import time
from collections.abc import Callable
from pathlib import Path

import polars as pl
from pyfaidx import Fasta

from healthcare.annotate.indel import Fetch, annotate_variant, normalize
from healthcare.annotate.snv import Annotation, Transcript, TranscriptIndex
from healthcare.data.protein_hgvs import equivalent

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
    "frameshift_variant": "frameshift",
    "inframe_deletion": "inframe_deletion",
    "inframe_insertion": "inframe_insertion",
    "inframe_indel": "inframe_indel",
}
HGVS_TYPES = [
    "missense",
    "stop_gained",
    "stop_lost",
    "frameshift",
    "inframe_deletion",
    "inframe_insertion",
    "inframe_indel",
]
MAX_EXAMPLES = 15
TYPE_EXAMPLES = 6


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


def is_indel() -> pl.Expr:
    return (pl.col("ref").str.len_chars() != 1) | (pl.col("alt").str.len_chars() != 1)


def genome_fetch(genome: Fasta) -> Callable[[str], Fetch]:
    def fetch_for(chrom: str) -> Fetch:
        sequence = genome[chrom]
        return lambda start, end: sequence[start - 1 : end].seq.upper()

    return fetch_for


def annotate_table(
    index: TranscriptIndex, variants: pl.DataFrame, fetch_for: Callable[[str], Fetch]
) -> tuple[pl.DataFrame, int]:
    rows, shifted = [], 0
    for chrom, pos, ref, alt in variants.select(KEY).iter_rows():
        position, reference, alternate = pos, ref, alt
        if len(ref) != 1 or len(alt) != 1:
            position, reference, alternate = normalize(pos, ref, alt, fetch_for(chrom))
            shifted += (position, reference, alternate) != (pos, ref, alt)
        for annotation in annotate_variant(index, chrom, position, reference, alternate):
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
    return pl.DataFrame(rows, schema=schema, orient="row"), shifted


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


def hgvs_by_type(labels: pl.DataFrame, annotations: pl.DataFrame, proteins: dict[str, str]) -> dict:
    compared = (
        labels.select(*KEY, "validation_transcript", "validation_hgvs_p")
        .join(annotations, on=KEY)
        .filter(
            (pl.col("validation_transcript") == pl.col("refseq_transcript"))
            & pl.col("consequence").is_in(HGVS_TYPES)
            & pl.col("validation_hgvs_p").is_not_null()
        )
        .with_columns(
            pl.concat_str(
                "consequence",
                pl.when(is_indel()).then(pl.lit("/indel")).otherwise(pl.lit("/snv")),
            ).alias("type"),
            (pl.col("validation_hgvs_p") == pl.col("hgvs_p")).alias("agrees"),
        )
    )
    compared = compared.with_columns(
        pl.struct("hgvs_p", "validation_hgvs_p", "transcript_id")
        .map_elements(
            lambda row: equivalent(
                row["hgvs_p"], row["validation_hgvs_p"], proteins[row["transcript_id"]]
            ),
            return_dtype=pl.Boolean,
        )
        .alias("equivalent")
    )
    report = {}
    for (kind,), group in compared.sort("type").group_by("type", maintain_order=True):
        report[kind] = {
            "compared": group.height,
            "agreement": round(group["agrees"].mean(), 5),
            "equivalent": round(group["equivalent"].mean(), 5),
            "disagreement_examples": group.filter(~pl.col("equivalent"))
            .select("gene_name", *KEY, "validation_hgvs_p", "hgvs_p")
            .head(TYPE_EXAMPLES)
            .to_dicts(),
        }
    return report


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
    joined = mapped.join(single.explode("consequence"), on=KEY).with_columns(
        pl.col("clinvar_consequences").list.contains(pl.col("consequence")).alias("hit"),
        is_indel().alias("indel"),
    )
    snvs, indels = joined.filter(~pl.col("indel")), joined.filter("indel")
    by_type = (
        indels.group_by("consequence")
        .agg(pl.len().alias("compared"), pl.col("hit").mean().round(5).alias("agreement"))
        .sort("compared", descending=True)
    )
    return {
        "consequence_compared": snvs.height,
        "consequence_agreement": round(snvs["hit"].mean(), 5) if snvs.height else None,
        "indel_consequence_compared": indels.height,
        "indel_consequence_agreement": round(indels["hit"].mean(), 5) if indels.height else None,
        "indel_consequence_by_type": {
            row["consequence"]: {"compared": row["compared"], "agreement": row["agreement"]}
            for row in by_type.to_dicts()
        },
    }


def build(
    transcripts_path: Path, genome_path: Path, clinvar_dir: Path, out_dir: Path, limit: int | None
) -> dict:
    started = time.monotonic()
    transcripts = load_transcripts(transcripts_path)
    index = TranscriptIndex(transcripts)
    proteins = {t.transcript_id: t.protein for t in transcripts}
    fetch_for = genome_fetch(Fasta(str(genome_path)))
    out_dir.mkdir(parents=True, exist_ok=True)
    summary: dict = {"out_dir": out_dir, "limit": limit}
    for name in ("train_labels", "temporal_labels", "vus"):
        variants = pl.read_parquet(clinvar_dir / f"{name}.parquet")
        variants = variants.head(limit) if limit else variants
        annotations, shifted = annotate_table(index, variants, fetch_for)
        annotations.write_parquet(out_dir / f"{name}.parquet")
        per_variant = annotations.group_by(KEY).len()
        summary[name] = {
            "snv_inputs": snvs_only(variants).height,
            "indel_inputs": variants.filter(is_indel()).height,
            "indels_moved_by_normalization": shifted,
            "annotated_variants": per_variant.height,
            "multi_transcript_variants": per_variant.filter(pl.col("len") > 1).height,
            "consequences": dict(
                annotations["consequence"].value_counts().sort("count", descending=True).iter_rows()
            ),
        }
        if name == "train_labels":
            summary[name] |= (
                hgvs_agreement(variants, annotations)
                | consequence_agreement(variants, annotations)
                | {"hgvs_by_type": hgvs_by_type(variants, annotations, proteins)}
            )
    summary["seconds"] = round(time.monotonic() - started, 1)
    return summary
