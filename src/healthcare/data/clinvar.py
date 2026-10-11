import gzip
import re
import time
from pathlib import Path

import polars as pl
from pyfaidx import Fasta

from healthcare.data.mane import PRIMARY_CHROMS

PATHOGENIC = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
BENIGN = {"Benign", "Likely_benign", "Benign/Likely_benign"}
UNCERTAIN = {"Uncertain_significance"}
REVIEW_STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
}
MIN_STARS = 1
KEY = ["chrom", "pos", "ref", "alt"]
VCF_COLUMNS = ["chrom", "pos", "variation_id", "ref", "alt", "qual", "filter", "info"]
HGVS_P = re.compile(r"\((p\.[^)]+)\)")
HGVS_TRANSCRIPT = re.compile(r"^(NM_[0-9.]+)")


def info_field(name: str) -> pl.Expr:
    return pl.col("info").str.extract(rf"(?:^|;){name}=([^;]*)", 1).alias(name)


def category_of(clnsig: str | None) -> str:
    if clnsig in PATHOGENIC:
        return "pathogenic"
    if clnsig in BENIGN:
        return "benign"
    if clnsig in UNCERTAIN:
        return "uncertain"
    return "other"


def read_vcf(path: Path) -> pl.DataFrame:
    raw = pl.read_csv(
        path,
        separator="\t",
        comment_prefix="#",
        has_header=False,
        new_columns=VCF_COLUMNS,
        schema_overrides={column: pl.String for column in VCF_COLUMNS},
        quote_char=None,
    )
    fields = ["CLNSIG", "CLNREVSTAT", "CLNVC", "GENEINFO", "MC", "ALLELEID"]
    return (
        raw.with_columns(info_field(name) for name in fields)
        .select(
            pl.col("chrom"),
            pl.col("pos").cast(pl.Int64),
            "ref",
            "alt",
            pl.col("variation_id").cast(pl.Int64).alias("clinvar_variation_id"),
            pl.col("ALLELEID").cast(pl.Int64).alias("clinvar_allele_id"),
            pl.col("GENEINFO").str.extract(r"^([^:|]+)", 1).alias("clinvar_gene"),
            pl.col("CLNVC").alias("clinvar_variant_type"),
            pl.col("CLNSIG").alias("clinvar_clnsig"),
            pl.col("CLNREVSTAT").alias("clinvar_review_status"),
            pl.col("MC").alias("validation_mc"),
        )
        .with_columns(
            pl.col("clinvar_review_status")
            .replace_strict(REVIEW_STARS, default=0, return_dtype=pl.Int8)
            .alias("clinvar_stars"),
            pl.col("clinvar_clnsig")
            .fill_null("")
            .map_elements(category_of, return_dtype=pl.String)
            .alias("clinvar_category"),
        )
    )


def keep_valid(table: pl.DataFrame) -> pl.DataFrame:
    return table.filter(
        pl.col("chrom").is_in(PRIMARY_CHROMS)
        & pl.col("ref").str.contains(r"^[ACGT]+$")
        & pl.col("alt").str.contains(r"^[ACGT]+$")
    )


def ref_matches(table: pl.DataFrame, genome: Fasta) -> pl.Series:
    matches = []
    for chrom, group in table.group_by("chrom", maintain_order=True):
        sequence = str(genome[chrom[0]][:]).upper()
        matches.append(
            group.select(
                pl.struct("pos", "ref")
                .map_elements(
                    lambda row, seq=sequence: (
                        seq[row["pos"] - 1 : row["pos"] - 1 + len(row["ref"])] == row["ref"]
                    ),
                    return_dtype=pl.Boolean,
                )
                .alias("ref_ok"),
                pl.col("row_id"),
            )
        )
    return pl.concat(matches).sort("row_id")["ref_ok"]


def deduplicate(table: pl.DataFrame) -> tuple[pl.DataFrame, int]:
    conflicting = (
        table.group_by(KEY)
        .agg(pl.col("clinvar_category").n_unique().alias("n"))
        .filter(pl.col("n") > 1)
    )
    cleaned = table.join(conflicting.select(KEY), on=KEY, how="anti").unique(
        KEY, keep="first", maintain_order=True
    )
    return cleaned, conflicting.height


def clean_release(path: Path, genome: Fasta) -> tuple[pl.DataFrame, dict]:
    raw = read_vcf(path)
    raw_keys = raw.select(KEY).unique()
    valid = keep_valid(raw).with_row_index("row_id")
    checked = valid.with_columns(ref_ok=ref_matches(valid, genome))
    ref_ok = checked.filter("ref_ok").drop("row_id", "ref_ok")
    deduplicated, conflicting_keys = deduplicate(ref_ok)
    stats = {
        "records": raw.height,
        "valid_primary_alleles": valid.height,
        "ref_mismatches": checked.height - ref_ok.height,
        "conflicting_duplicate_keys": conflicting_keys,
        "unique_variants": deduplicated.height,
    }
    return deduplicated, stats | {"raw_keys": raw_keys}


def labeled(table: pl.DataFrame) -> pl.DataFrame:
    return table.filter(
        pl.col("clinvar_category").is_in(["pathogenic", "benign"])
        & (pl.col("clinvar_stars") >= MIN_STARS)
    ).with_columns((pl.col("clinvar_category") == "pathogenic").cast(pl.Int8).alias("label"))


def uncertain(table: pl.DataFrame) -> pl.DataFrame:
    return table.filter(
        (pl.col("clinvar_category") == "uncertain") & (pl.col("clinvar_stars") >= MIN_STARS)
    )


def temporal_and_reclassified(
    train_all: pl.DataFrame, test_all: pl.DataFrame, train_raw_keys: pl.DataFrame | None = None
) -> tuple[pl.DataFrame, pl.DataFrame]:
    test_labeled = labeled(test_all)
    seen = train_all.select(KEY) if train_raw_keys is None else train_raw_keys
    temporal = test_labeled.join(seen, on=KEY, how="anti")
    previous = train_all.select(*KEY, pl.col("clinvar_category").alias("previous_category"))
    reclassified = (
        test_labeled.join(previous, on=KEY, how="inner")
        .filter(pl.col("previous_category") != pl.col("clinvar_category"))
        .filter(pl.col("previous_category").is_in(["pathogenic", "benign", "uncertain"]))
        .with_columns(
            pl.when(pl.col("previous_category") == "uncertain")
            .then(pl.lit("vus_resolved"))
            .otherwise(pl.lit("label_flipped"))
            .alias("change")
        )
    )
    return temporal, reclassified


def read_hgvs(path: Path) -> pl.DataFrame:
    rows = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = handle.readline().lstrip("#").rstrip("\n").split("\t")
        assembly, name, variation = (header.index(c) for c in ("Assembly", "Name", "VariationID"))
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[assembly] != "GRCh38":
                continue
            protein = HGVS_P.search(fields[name])
            transcript = HGVS_TRANSCRIPT.search(fields[name])
            rows[int(fields[variation])] = (
                transcript.group(1) if transcript else None,
                protein.group(1) if protein else None,
            )
    return pl.DataFrame(
        {
            "clinvar_variation_id": list(rows),
            "validation_transcript": [value[0] for value in rows.values()],
            "validation_hgvs_p": [value[1] for value in rows.values()],
        },
        schema={
            "clinvar_variation_id": pl.Int64,
            "validation_transcript": pl.String,
            "validation_hgvs_p": pl.String,
        },
    )


def distribution(table: pl.DataFrame) -> list[dict]:
    return (
        table.group_by("clinvar_category", "clinvar_stars")
        .len()
        .sort("clinvar_category", "clinvar_stars")
        .to_dicts()
    )


def top_types(table: pl.DataFrame, n: int = 6) -> dict:
    counts = table.group_by("clinvar_variant_type", "label").len().sort("len", descending=True)
    return {
        f"{row['clinvar_variant_type']}|{row['label']}": row["len"]
        for row in counts.head(n * 2).to_dicts()
    }


def build(
    train_vcf: Path, test_vcf: Path, summary_tsv: Path, reference_fa: Path, out_dir: Path
) -> dict:
    started = time.monotonic()
    genome = Fasta(str(reference_fa))
    train_all, train_stats = clean_release(train_vcf, genome)
    test_all, test_stats = clean_release(test_vcf, genome)
    train_raw_keys = train_stats.pop("raw_keys")
    test_stats.pop("raw_keys")
    hgvs = read_hgvs(summary_tsv)

    train = labeled(train_all).join(hgvs, on="clinvar_variation_id", how="left")
    vus = uncertain(train_all).join(hgvs, on="clinvar_variation_id", how="left")
    temporal, reclassified = temporal_and_reclassified(train_all, test_all, train_raw_keys)
    newly_labeled_from_unused = (
        labeled(test_all)
        .join(train_all.filter(pl.col("clinvar_stars") == 0).select(KEY), on=KEY, how="semi")
        .height
        + labeled(test_all)
        .join(
            train_all.filter(pl.col("clinvar_category") == "other").select(KEY), on=KEY, how="semi"
        )
        .height
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "train_labels": train,
        "vus": vus,
        "temporal_labels": temporal,
        "reclassified": reclassified,
    }
    for name, table in outputs.items():
        table.write_parquet(out_dir / f"{name}.parquet")

    return {
        "out_dir": out_dir,
        "train_release": train_stats,
        "temporal_release": test_stats,
        "train_labels": train.height,
        "train_pathogenic": int(train["label"].sum()),
        "train_benign": int(train.height - train["label"].sum()),
        "train_category_by_stars": distribution(train_all),
        "train_top_types": top_types(train),
        "vus": vus.height,
        "temporal_labels": temporal.height,
        "temporal_pathogenic": int(temporal["label"].sum()),
        "reclassified": dict(reclassified["change"].value_counts().iter_rows()),
        "dropped_previously_unusable_now_labeled": newly_labeled_from_unused,
        "hgvs_p_coverage_snv": round(
            train.filter(pl.col("clinvar_variant_type") == "single_nucleotide_variant")[
                "validation_hgvs_p"
            ]
            .is_not_null()
            .mean(),
            4,
        ),
        "seconds": round(time.monotonic() - started, 1),
    }
