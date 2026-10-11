import time
from pathlib import Path

import polars as pl

from healthcare.annotate.snv import Annotation, Transcript, TranscriptIndex

ANNOTATION_FIELDS = list(Annotation.__dataclass_fields__)
KEY = ["chrom", "pos", "ref", "alt"]
MC_TO_CONSEQUENCE = {
    "missense_variant": "missense",
    "synonymous_variant": "synonymous",
    "nonsense": "stop_gained",
    "stop_lost": "stop_lost",
    "initiator_codon_variant": "start_lost",
    "start_lost": "start_lost",
    "splice_donor_variant": "canonical_splice",
    "splice_acceptor_variant": "canonical_splice",
    "intron_variant": "intronic",
    "5_prime_UTR_variant": "utr",
    "3_prime_UTR_variant": "utr",
}
MAX_EXAMPLES = 15


def load_transcripts(table_path: Path) -> list[Transcript]:
    transcripts = []
    for row in pl.read_parquet(table_path).filter("translation_matches").iter_rows(named=True):
        transcripts.append(
            Transcript(
                transcript_id=row["transcript_id"],
                gene_name=row["gene_name"],
                refseq_transcript=row["refseq_transcript"],
                chrom=row["chrom"],
                strand=row["strand"],
                exons=tuple(zip(row["exon_starts"], row["exon_ends"], strict=True)),
                cds=tuple(zip(row["cds_starts"], row["cds_ends"], strict=True)),
                cds_sequence=row["cds_sequence"],
                protein=row["protein"],
            )
        )
    return transcripts


def snvs_only(variants: pl.DataFrame) -> pl.DataFrame:
    return variants.filter(
        (pl.col("ref").str.len_chars() == 1) & (pl.col("alt").str.len_chars() == 1)
    )


def annotate_table(index: TranscriptIndex, variants: pl.DataFrame) -> pl.DataFrame:
    rows = []
    for chrom, pos, ref, alt in snvs_only(variants).select(KEY).iter_rows():
        for annotation in index.annotate(chrom, pos, ref, alt):
            rows.append(
                (chrom, pos, ref, alt)
                + tuple(getattr(annotation, field) for field in ANNOTATION_FIELDS)
            )
    schema = {
        "chrom": pl.String,
        "pos": pl.Int64,
        "ref": pl.String,
        "alt": pl.String,
        "transcript_id": pl.String,
        "gene_name": pl.String,
        "refseq_transcript": pl.String,
        "consequence": pl.String,
        "cds_position": pl.Int64,
        "protein_position": pl.Int64,
        "codon_ref": pl.String,
        "codon_alt": pl.String,
        "aa_ref": pl.String,
        "aa_alt": pl.String,
        "hgvs_p": pl.String,
    }
    return pl.DataFrame(rows, schema=schema, orient="row")


def hgvs_agreement(labels: pl.DataFrame, annotations: pl.DataFrame) -> dict:
    compared = (
        labels.select(*KEY, "validation_transcript", "validation_hgvs_p")
        .join(annotations, on=KEY)
        .filter(
            (pl.col("validation_transcript") == pl.col("refseq_transcript"))
            & (pl.col("consequence") == "missense")
            & pl.col("validation_hgvs_p").is_not_null()
        )
    )
    agrees = compared.filter(pl.col("validation_hgvs_p") == pl.col("hgvs_p"))
    examples = (
        compared.join(agrees.select(KEY), on=KEY, how="anti")
        .select("gene_name", *KEY, "validation_hgvs_p", "hgvs_p")
        .head(MAX_EXAMPLES)
        .to_dicts()
    )
    return {
        "missense_compared": compared.height,
        "missense_agreement": round(agrees.height / compared.height, 5)
        if compared.height
        else None,
        "missense_disagreement_examples": examples,
    }


def consequence_agreement(labels: pl.DataFrame, annotations: pl.DataFrame) -> dict:
    mapped = labels.select(
        *KEY,
        pl.col("validation_mc")
        .str.extract_all(r"\|([A-Za-z0-9_']+)")
        .list.eval(
            pl.element().str.strip_prefix("|").replace_strict(MC_TO_CONSEQUENCE, default=None)
        )
        .list.drop_nulls()
        .alias("clinvar_consequences"),
    ).filter(pl.col("clinvar_consequences").list.len() > 0)
    single = (
        annotations.group_by(KEY)
        .agg(pl.col("consequence"))
        .filter(pl.col("consequence").list.len() == 1)
    )
    joined = mapped.join(single.explode("consequence"), on=KEY)
    hit = joined.filter(pl.col("clinvar_consequences").list.contains(pl.col("consequence")))
    return {
        "consequence_compared": joined.height,
        "consequence_agreement": round(hit.height / joined.height, 5) if joined.height else None,
    }


def build(transcripts_path: Path, clinvar_dir: Path, out_dir: Path, limit: int | None) -> dict:
    started = time.monotonic()
    index = TranscriptIndex(load_transcripts(transcripts_path))
    out_dir.mkdir(parents=True, exist_ok=True)
    summary: dict = {"out_dir": out_dir, "limit": limit}
    for name in ("train_labels", "temporal_labels", "vus"):
        variants = pl.read_parquet(clinvar_dir / f"{name}.parquet")
        variants = variants.head(limit) if limit else variants
        annotations = annotate_table(index, variants)
        annotations.write_parquet(out_dir / f"{name}.parquet")
        per_variant = annotations.group_by(KEY).len()
        summary[name] = {
            "snv_inputs": snvs_only(variants).height,
            "annotated_variants": per_variant.height,
            "multi_transcript_variants": per_variant.filter(pl.col("len") > 1).height,
            "consequences": dict(
                annotations["consequence"].value_counts().sort("count", descending=True).iter_rows()
            ),
        }
        if name == "train_labels":
            summary[name] |= hgvs_agreement(variants, annotations) | consequence_agreement(
                variants, annotations
            )
    summary["seconds"] = round(time.monotonic() - started, 1)
    return summary
