import time
from dataclasses import dataclass
from pathlib import Path

import duckdb
import numpy as np
import polars as pl

from healthcare.features.aminoacids import ORDER, substitution_features
from healthcare.features.esm_features import esm_features

KEY = ["chrom", "pos", "ref", "alt"]
CATEGORICAL_FEATURES = ["f_aa_ref", "f_aa_alt"]
AMINO_ACID_CODES = {aa: i for i, aa in enumerate(ORDER + "U*X")}
FREQUENCY_FEATURES = [
    "f_gnomad_af",
    "f_gnomad_ac",
    "f_gnomad_nhomalt",
    "f_gnomad_af_grpmax",
    "f_gnomad_filter_pass",
    "f_gnomad_covered",
    "f_gnomad_exome_an",
]
UNIPROT_GROUPS = {
    "f_uniprot_domain": ["DOMAIN", "REPEAT", "ZN_FING", "DNA_BIND"],
    "f_uniprot_active_site": ["ACT_SITE"],
    "f_uniprot_binding_site": ["BINDING", "SITE"],
    "f_uniprot_motif": ["MOTIF"],
    "f_uniprot_region": ["REGION"],
    "f_uniprot_coiled_coil": ["COILED"],
    "f_uniprot_transmembrane": ["TRANSMEM", "INTRAMEM"],
    "f_uniprot_signal_or_propeptide": ["SIGNAL", "TRANSIT", "PROPEP", "PEPTIDE"],
    "f_uniprot_disulfide": ["DISULFID"],
    "f_uniprot_ptm_site": ["MOD_RES", "CARBOHYD", "LIPID", "CROSSLNK"],
    "f_uniprot_helix": ["HELIX"],
    "f_uniprot_strand": ["STRAND"],
    "f_uniprot_turn": ["TURN"],
}
CONSTRAINT_FEATURES = {
    "loeuf": "f_gene_loeuf",
    "lof_oe": "f_gene_lof_oe",
    "pli": "f_gene_pli",
    "mis_oe": "f_gene_mis_oe",
    "mis_oe_upper": "f_gene_mis_oe_upper",
    "mis_z": "f_gene_mis_z",
    "syn_z": "f_gene_syn_z",
}
ID_COLUMNS = [*KEY, "gene_name", "transcript_id", "protein_position", "aa_ref", "aa_alt"]


@dataclass(frozen=True)
class Inputs:
    transcripts: Path
    gnomad_exomes: Path
    coverage: Path
    conservation: Path
    constraint: Path
    uniprot: Path
    esm_store: Path


def resolve_transcripts(missense: pl.DataFrame) -> tuple[pl.DataFrame, int]:
    counts = missense.group_by(KEY).len()
    single = counts.filter(pl.col("len") == 1).select(KEY)
    return missense.join(single, on=KEY), counts.height - single.height


def amino_acid_features(table: pl.DataFrame) -> pl.DataFrame:
    rows = [substitution_features(a, b) for a, b in table.select("aa_ref", "aa_alt").iter_rows()]
    return table.hstack(pl.DataFrame(rows, schema={name: pl.Float64 for name in rows[0]}))


def region_filter(variants: pl.DataFrame) -> str:
    bounds = variants.group_by("chrom").agg(pl.col("pos").min().alias("lo"), pl.col("pos").max())
    terms = [
        f"(chrom = '{chrom}' AND pos BETWEEN {lo} AND {hi})"
        for chrom, lo, hi in bounds.sort("chrom").iter_rows()
    ]
    return " OR ".join(terms) or "FALSE"


def join_sources(variants: pl.DataFrame, inputs: Inputs) -> pl.DataFrame:
    con = duckdb.connect()
    con.register("variants", variants.to_arrow())
    regions = region_filter(variants)
    group_columns = ",\n".join(
        f"bool_or(u.feature_type IN ({', '.join(repr(t) for t in types)})) AS {name}"
        for name, types in UNIPROT_GROUPS.items()
    )
    query = f"""
    WITH g AS (
        SELECT * FROM read_parquet('{(inputs.gnomad_exomes / "*.parquet").as_posix()}')
        WHERE {regions}
    ),
    c AS (
        SELECT * FROM read_parquet('{inputs.coverage.as_posix()}') WHERE {regions}
    ),
    k AS (
        SELECT * FROM read_parquet('{(inputs.conservation / "*.parquet").as_posix()}')
        WHERE {regions}
    ),
    u AS (
        SELECT v.chrom, v.pos, v.ref, v.alt, {group_columns}
        FROM variants v
        JOIN read_parquet('{inputs.uniprot.as_posix()}') u
          ON u.transcript_id = v.transcript_id
         AND v.protein_position BETWEEN u.start AND u."end"
        GROUP BY v.chrom, v.pos, v.ref, v.alt
    ),
    mapped AS (
        SELECT DISTINCT transcript_id FROM read_parquet('{inputs.uniprot.as_posix()}')
    )
    SELECT
        v.*,
        g.af AS gnomad_af_raw, g.ac AS gnomad_ac_raw, g.nhomalt AS gnomad_nhomalt_raw,
        g.af_grpmax AS gnomad_af_grpmax_raw, g.filter_pass AS gnomad_filter_pass_raw,
        c.exome_an AS f_gnomad_exome_an,
        k.phylop AS f_phylop, k.phastcons AS f_phastcons,
        t.cds_length / 3 - 1 AS f_protein_length,
        m.transcript_id IS NOT NULL AS uniprot_mapped,
        {", ".join(f"u.{name}" for name in UNIPROT_GROUPS)},
        {", ".join(f"q.{src} AS {dst}" for src, dst in CONSTRAINT_FEATURES.items())}
    FROM variants v
    LEFT JOIN g USING (chrom, pos, ref, alt)
    LEFT JOIN c USING (chrom, pos)
    LEFT JOIN k USING (chrom, pos)
    LEFT JOIN read_parquet('{inputs.transcripts.as_posix()}') t USING (transcript_id)
    LEFT JOIN mapped m USING (transcript_id)
    LEFT JOIN u USING (chrom, pos, ref, alt)
    LEFT JOIN read_parquet('{inputs.constraint.as_posix()}') q
      ON q.transcript_id = split_part(v.transcript_id, '.', 1)
    """
    return con.execute(query).pl()


def shared_columns() -> list[pl.Expr]:
    covered = pl.col("f_gnomad_exome_an").fill_null(0) > 0
    in_gnomad = pl.col("gnomad_ac_raw").is_not_null()
    return [
        covered.cast(pl.Float64).alias("f_gnomad_covered"),
        *[
            pl.when(in_gnomad)
            .then(pl.col(f"gnomad_{name}_raw"))
            .when(covered)
            .then(0.0)
            .otherwise(None)
            .cast(pl.Float64)
            .alias(f"f_gnomad_{name}")
            for name in ("af", "ac", "nhomalt", "af_grpmax")
        ],
        pl.col("gnomad_filter_pass_raw").cast(pl.Float64).alias("f_gnomad_filter_pass"),
        (pl.col("protein_position") / pl.col("f_protein_length")).alias("f_relative_position"),
        *[
            pl.when(pl.col("uniprot_mapped"))
            .then(pl.col(name).fill_null(False))
            .otherwise(None)
            .cast(pl.Float64)
            .alias(name)
            for name in UNIPROT_GROUPS
        ],
        pl.col("uniprot_mapped").cast(pl.Float64).alias("f_uniprot_mapped"),
    ]


def finalize(joined: pl.DataFrame) -> pl.DataFrame:
    table = joined.with_columns(
        *shared_columns(),
        pl.col("aa_ref").alias("f_aa_ref"),
        pl.col("aa_alt").alias("f_aa_alt"),
    )
    table = amino_acid_features(table)
    features = sorted(c for c in table.columns if c.startswith("f_"))
    labels = [c for c in ("label",) if c in table.columns]
    return table.select(*ID_COLUMNS, *labels, *features).sort(KEY)


def to_matrix(table: pl.DataFrame, columns: list[str]) -> np.ndarray:
    encoded = table.select(
        [
            pl.col(c).replace_strict(AMINO_ACID_CODES, default=None, return_dtype=pl.Float64)
            if c in CATEGORICAL_FEATURES
            else pl.col(c).cast(pl.Float64)
            for c in columns
        ]
    )
    return encoded.to_numpy()


def missing_rates(table: pl.DataFrame) -> dict:
    features = [c for c in table.columns if c.startswith("f_")]
    return {c: round(table[c].null_count() / table.height, 4) for c in features}


def build(annotations_dir: Path, clinvar_dir: Path, inputs: Inputs, out_dir: Path) -> dict:
    started = time.monotonic()
    out_dir.mkdir(parents=True, exist_ok=True)
    summary: dict = {}
    for name in ("train_labels", "temporal_labels", "vus"):
        annotations = pl.read_parquet(annotations_dir / f"{name}.parquet").filter(
            pl.col("consequence") == "missense"
        )
        clinvar = pl.read_parquet(clinvar_dir / f"{name}.parquet")
        clinvar = clinvar.select(*KEY, *(["label"] if "label" in clinvar.columns else []))
        missense = annotations.join(clinvar, on=KEY)
        resolved, dropped = resolve_transcripts(missense)
        table = finalize(esm_features(join_sources(resolved, inputs), inputs.esm_store))
        table.write_parquet(out_dir / f"missense_{name}.parquet")
        entry = {"rows": table.height, "dropped_multi_gene": dropped}
        if "label" in table.columns:
            entry["pathogenic"] = int(table["label"].sum())
            entry["genes"] = table["gene_name"].n_unique()
        summary[name] = entry
        if name == "train_labels":
            summary["feature_columns"] = len([c for c in table.columns if c.startswith("f_")])
            summary["missing_rates"] = missing_rates(table)
    summary["seconds"] = round(time.monotonic() - started, 1)
    return summary
