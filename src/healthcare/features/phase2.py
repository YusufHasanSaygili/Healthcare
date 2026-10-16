import json
import time
from collections.abc import Callable
from pathlib import Path

import duckdb
import polars as pl
from pyfaidx import Fasta

from healthcare.annotate.indel import edit
from healthcare.annotate.snv import Transcript
from healthcare.annotate.truncation import (
    boundary_distance,
    cds_edit,
    splice_features,
    start_loss_features,
    truncation_features,
    utr3_sequence,
)
from healthcare.data.annotations import load_transcripts
from healthcare.data.protein_hgvs import residues
from healthcare.esm import store
from healthcare.features import missense
from healthcare.features.esm_features import EMPTY_SPAN, SPAN_FEATURES, span_features

KEY = missense.KEY
PHASE2_TYPES = [
    "synonymous",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "canonical_splice",
    "frameshift",
    "inframe_deletion",
    "inframe_insertion",
    "inframe_indel",
]
VARIANT_TYPES = ["missense", *PHASE2_TYPES]
TRUNCATING = {"stop_gained", "stop_lost", "frameshift"}
IN_FRAME = {"inframe_deletion", "inframe_insertion", "inframe_indel"}
TYPE_FEATURES = [
    "f_variant_type",
    "f_indel_length",
    "f_boundary_distance",
    "f_no_stop",
    "f_protein_retained",
    "f_new_length_ratio",
    "f_novel_residues",
    "f_ptc_junction_distance",
    "f_ptc_in_last_exon",
    "f_nmd_escape",
    "f_next_met_fraction",
    "f_splice_donor",
    "f_splice_intron_position",
    "f_skipped_exon_in_frame",
    "f_skipped_exon_cds_fraction",
    "f_skipped_exon_rank",
    "f_residues_removed",
    "f_residues_added",
]
ID_COLUMNS = [*KEY, "gene_name", "transcript_id", "consequence", "protein_position", "hgvs_p"]


def resolve(annotations: pl.DataFrame) -> tuple[pl.DataFrame, int]:
    scored = annotations.filter(pl.col("consequence").is_in(VARIANT_TYPES))
    genes = scored.group_by(KEY).agg(pl.col("gene_name").n_unique().alias("genes"))
    single = genes.filter(pl.col("genes") == 1).select(KEY)
    phase2 = scored.filter(pl.col("consequence").is_in(PHASE2_TYPES))
    kept = phase2.join(single, on=KEY)
    return kept, phase2.select(KEY).unique().height - kept.select(KEY).unique().height


def variant_features(
    row: dict, transcript: Transcript, utr3: Callable[[Transcript], str]
) -> dict[str, float | int | None]:
    start, deleted, inserted = edit(row["pos"], row["ref"], row["alt"])
    span_start = start if deleted else start - 1
    span_end = start + len(deleted) - 1 if deleted else start
    consequence = row["consequence"]
    values: dict[str, float | int | None] = {
        "span_start": span_start,
        "span_end": span_end,
        "f_variant_type": float(VARIANT_TYPES.index(consequence)),
        "f_indel_length": float(len(inserted) - len(deleted)),
        "f_boundary_distance": boundary_distance(transcript, span_start),
    }
    if consequence in TRUNCATING:
        change = cds_edit(transcript, row["pos"], row["ref"], row["alt"])
        if change is not None:
            values |= truncation_features(transcript, change, utr3(transcript))
    elif consequence == "start_lost":
        values |= start_loss_features(transcript)
    elif consequence == "canonical_splice":
        values |= splice_features(transcript, span_start, span_end)
    if consequence in IN_FRAME:
        removed = len(row["aa_ref"] or "")
        values["f_residues_removed"] = float(removed)
        values["f_residues_added"] = removed + (len(inserted) - len(deleted)) / 3
    return values


def type_features(
    table: pl.DataFrame, transcripts: dict[str, Transcript], utr3: Callable[[Transcript], str]
) -> pl.DataFrame:
    rows = [
        variant_features(row, transcripts[row["transcript_id"]], utr3)
        for row in table.select(*KEY, "transcript_id", "consequence", "aa_ref").iter_rows(
            named=True
        )
    ]
    schema = {"span_start": pl.Int64, "span_end": pl.Int64} | {
        name: pl.Float64 for name in TYPE_FEATURES
    }
    return table.hstack(pl.DataFrame(rows, schema=schema))


def inserted_residues(hgvs: str | None, protein: str, position: int, added: int) -> str:
    if not hgvs or added <= 0:
        return ""
    if hgvs.endswith("dup"):
        return protein[position - 1 - added : position - 1]
    for marker in ("delins", "ins"):
        if marker in hgvs:
            return residues(hgvs.rsplit(marker, 1)[1])
    return residues(hgvs[-3:])


def indel_esm_features(
    table: pl.DataFrame, transcripts: dict[str, Transcript], store_dir: Path
) -> pl.DataFrame:
    values: list[tuple] = [EMPTY_SPAN] * table.height
    in_frame = (
        table.with_row_index()
        .filter(pl.col("consequence").is_in(IN_FRAME))
        .select(
            "index", "transcript_id", "protein_position", "aa_ref", "hgvs_p", "f_residues_added"
        )
    )
    for (transcript_id,), group in in_frame.group_by("transcript_id"):
        if not store.is_done(store_dir, transcript_id):
            continue
        log_probs, _ = store.load(store_dir, transcript_id)
        protein = transcripts[transcript_id].protein
        for row, position, removed, hgvs, added in group.drop("transcript_id").iter_rows():
            inserted = inserted_residues(hgvs, protein, position, int(added))
            values[row] = span_features(log_probs, position, removed or "", inserted)
    return table.hstack(
        pl.DataFrame(values, schema={name: pl.Float64 for name in SPAN_FEATURES}, orient="row")
    )


def span_conservation(variants: pl.DataFrame, conservation: Path) -> pl.DataFrame:
    con = duckdb.connect()
    con.register("variants", variants.select(*KEY, "span_start", "span_end").to_arrow())
    regions = " OR ".join(
        f"(chrom = '{chrom}' AND pos BETWEEN {lo} AND {hi})"
        for chrom, lo, hi in variants.group_by("chrom")
        .agg(pl.col("span_start").min(), pl.col("span_end").max())
        .sort("chrom")
        .iter_rows()
    )
    query = f"""
    WITH k AS (
        SELECT * FROM read_parquet('{(conservation / "*.parquet").as_posix()}')
        WHERE {regions or "FALSE"}
    )
    SELECT v.chrom, v.pos, v.ref, v.alt,
           avg(k.phylop) AS f_phylop, max(k.phylop) AS f_phylop_max,
           avg(k.phastcons) AS f_phastcons
    FROM variants v
    JOIN k ON k.chrom = v.chrom AND k.pos BETWEEN v.span_start AND v.span_end
    GROUP BY v.chrom, v.pos, v.ref, v.alt
    """
    return con.execute(query).pl()


def finalize(joined: pl.DataFrame, conservation: pl.DataFrame) -> pl.DataFrame:
    position_known = pl.col("protein_position").is_not_null()
    table = (
        joined.drop("f_phylop", "f_phastcons")
        .join(conservation, on=KEY, how="left")
        .with_columns(*missense.shared_columns())
        .with_columns(
            pl.when(position_known).then(pl.col(name)).otherwise(None).alias(name)
            for name in missense.UNIPROT_GROUPS
        )
    )
    features = sorted(c for c in table.columns if c.startswith("f_"))
    labels = [c for c in ("label",) if c in table.columns]
    return table.select(*ID_COLUMNS, *labels, *features).sort(KEY)


def missing_by_type(table: pl.DataFrame) -> dict:
    features = [c for c in table.columns if c.startswith("f_")]
    return {
        consequence: {c: round(group[c].null_count() / group.height, 4) for c in features}
        for (consequence,), group in table.sort("consequence").group_by(
            "consequence", maintain_order=True
        )
    }


def utr3_lookup(genome: Fasta) -> Callable[[Transcript], str]:
    cache: dict[str, str] = {}

    def lookup(transcript: Transcript) -> str:
        if transcript.transcript_id not in cache:
            sequence = genome[transcript.chrom]
            cache[transcript.transcript_id] = utr3_sequence(
                transcript, lambda start, end: sequence[start - 1 : end].seq.upper()
            )
        return cache[transcript.transcript_id]

    return lookup


def build(
    annotations_dir: Path,
    clinvar_dir: Path,
    inputs: missense.Inputs,
    genome_path: Path,
    out_dir: Path,
    report_path: Path,
    limit: int | None,
) -> dict:
    started = time.monotonic()
    transcripts = {t.transcript_id: t for t in load_transcripts(inputs.transcripts)}
    utr3 = utr3_lookup(Fasta(str(genome_path)))
    out_dir.mkdir(parents=True, exist_ok=True)
    summary: dict = {"limit": limit}
    for name in ("train_labels", "temporal_labels", "vus"):
        annotations = pl.read_parquet(annotations_dir / f"{name}.parquet")
        clinvar = pl.read_parquet(clinvar_dir / f"{name}.parquet")
        clinvar = clinvar.select(*KEY, *(["label"] if "label" in clinvar.columns else []))
        resolved, dropped = resolve(annotations)
        resolved = resolved.join(clinvar, on=KEY)
        resolved = resolved.head(limit) if limit else resolved
        typed = indel_esm_features(
            type_features(resolved, transcripts, utr3), transcripts, inputs.esm_store
        )
        joined = missense.join_sources(typed, inputs)
        table = finalize(joined, span_conservation(typed, inputs.conservation))
        table.write_parquet(out_dir / f"phase2_{name}.parquet")
        counts = table.group_by("consequence").agg(
            pl.len().alias("rows"),
            *([pl.col("label").sum().alias("pathogenic")] if "label" in table.columns else []),
        )
        summary[name] = {
            "rows": table.height,
            "dropped_multi_gene": dropped,
            "by_type": {row.pop("consequence"): row for row in counts.sort("rows").to_dicts()},
        }
        if name == "train_labels":
            report = missing_by_type(table)
            report_path.parent.mkdir(parents=True, exist_ok=True)
            report_path.write_text(json.dumps(report, indent=1), encoding="utf-8")
            summary["feature_columns"] = len([c for c in table.columns if c.startswith("f_")])
            summary["missing_report"] = report_path
    summary["seconds"] = round(time.monotonic() - started, 1)
    return summary
