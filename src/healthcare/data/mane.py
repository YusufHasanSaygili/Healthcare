import gzip
import shutil
import time
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path

import polars as pl
from pyfaidx import Fasta

from healthcare.annotate.sequence import reverse_complement, translate

CODING_FEATURES = {"CDS", "stop_codon"}
MAX_REPORTED_MISMATCHES = 20
PRIMARY_CHROMS = {*map(str, range(1, 23)), "X", "Y"}


def normalize_chrom(chrom: str) -> str:
    chrom = chrom.removeprefix("chr")
    return "MT" if chrom == "M" else chrom


def parse_attributes(field: str) -> dict[str, list[str]]:
    attributes = defaultdict(list)
    for item in field.strip().rstrip(";").split(";"):
        key, _, value = item.strip().partition(" ")
        attributes[key].append(value.strip('"'))
    return attributes


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def parse_mane_gtf(lines: Iterable[str]) -> list[dict]:
    transcripts: dict[str, dict] = {}
    for line in lines:
        if line.startswith("#"):
            continue
        chrom, _, feature, start, end, _, strand, _, attribute_field = line.rstrip("\n").split("\t")
        if feature not in {"transcript", "exon", *CODING_FEATURES}:
            continue
        attributes = parse_attributes(attribute_field)
        if "MANE_Select" not in attributes["tag"]:
            continue
        transcript_id = attributes["transcript_id"][0]
        record = transcripts.setdefault(
            transcript_id,
            {
                "gene_name": attributes["gene_name"][0],
                "gene_id": attributes["gene_id"][0],
                "transcript_id": transcript_id,
                "protein_id": None,
                "refseq_transcript": None,
                "chrom": normalize_chrom(chrom),
                "strand": strand,
                "exons": [],
                "coding": [],
            },
        )
        interval = (int(start), int(end))
        if feature == "transcript":
            refseq = [x for x in attributes["db_xref"] if x.startswith("RefSeq:NM_")]
            record["refseq_transcript"] = refseq[0].removeprefix("RefSeq:") if refseq else None
        elif feature == "exon":
            record["exons"].append(interval)
        else:
            record["coding"].append(interval)
        if attributes["protein_id"]:
            record["protein_id"] = attributes["protein_id"][0]
    return [record for record in transcripts.values() if record["coding"]]


def read_protein_fasta(path: Path) -> dict[str, str]:
    proteins: dict[str, list[str]] = {}
    current = None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                current = line[1:].split()[0]
                proteins[current] = []
            elif current:
                proteins[current].append(line.strip())
    return {key: "".join(parts) for key, parts in proteins.items()}


def coding_sequence(genome: Fasta, chrom: str, strand: str, coding: list[tuple[int, int]]) -> str:
    sequence = "".join(genome[chrom][start - 1 : end].seq.upper() for start, end in coding)
    return reverse_complement(sequence) if strand == "-" else sequence


ACCEPTED_NOTES = {"exact", "non_atg_start", "selenocysteine"}


def translation_note(translated: str, expected: str | None) -> str:
    if expected is None:
        return "missing official protein"
    if not translated.endswith("*"):
        return "no stop codon"
    body = translated[:-1]
    if len(body) != len(expected):
        return "length differs"
    notes = []
    if body[0] != expected[0] and expected[0] == "M":
        notes.append("non_atg_start")
        body = "M" + body[1:]
    if "U" in expected and body.replace("*", "U") == expected and body != expected:
        notes.append("selenocysteine")
        body = expected
    if body != expected:
        return "sequence differs"
    return "+".join(notes) or "exact"


def is_accepted(note: str) -> bool:
    return set(note.split("+")) <= ACCEPTED_NOTES


def build_transcript_table(
    genome: Fasta, transcripts: list[dict], proteins: dict[str, str]
) -> tuple[pl.DataFrame, list[dict]]:
    rows, mismatches = [], []
    for record in transcripts:
        coding = merge_intervals(record["coding"])
        exons = merge_intervals(record["exons"])
        cds = coding_sequence(genome, record["chrom"], record["strand"], coding)
        official = proteins.get(record["protein_id"])
        note = translation_note(translate(cds), official)
        if not is_accepted(note):
            mismatches.append(
                {
                    "transcript_id": record["transcript_id"],
                    "gene": record["gene_name"],
                    "reason": note,
                }
            )
        rows.append(
            {
                "gene_name": record["gene_name"],
                "gene_id": record["gene_id"],
                "transcript_id": record["transcript_id"],
                "protein_id": record["protein_id"],
                "refseq_transcript": record["refseq_transcript"],
                "chrom": record["chrom"],
                "strand": record["strand"],
                "exon_starts": [s for s, _ in exons],
                "exon_ends": [e for _, e in exons],
                "cds_starts": [s for s, _ in coding],
                "cds_ends": [e for _, e in coding],
                "cds_length": len(cds),
                "protein": official,
                "translation_note": note,
                "translation_matches": is_accepted(note),
            }
        )
    table = pl.DataFrame(rows).sort("chrom", pl.col("cds_starts").list.first())
    return table, mismatches


def decompress(source: Path, dest: Path) -> Path:
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        part = dest.with_name(dest.name + ".part")
        with gzip.open(source, "rb") as src, part.open("wb") as out:
            shutil.copyfileobj(src, out, 16 * 1024 * 1024)
        part.replace(dest)
    return dest


def build(
    fasta_gz: Path,
    gtf_gz: Path,
    protein_gz: Path,
    reference_fa: Path,
    output: Path,
    chrom: str | None = None,
    limit: int | None = None,
) -> dict:
    started = time.monotonic()
    genome = Fasta(str(decompress(fasta_gz, reference_fa)))
    with gzip.open(gtf_gz, "rt") as handle:
        parsed = parse_mane_gtf(handle)
    transcripts = [t for t in parsed if t["chrom"] in PRIMARY_CHROMS]
    primary_genes = {t["gene_name"] for t in transcripts}
    patch_only = sorted({t["gene_name"] for t in parsed} - primary_genes)
    if chrom:
        transcripts = [t for t in transcripts if t["chrom"] == chrom]
    transcripts = transcripts[:limit] if limit else transcripts
    table, mismatches = build_transcript_table(genome, transcripts, read_protein_fasta(protein_gz))
    output.parent.mkdir(parents=True, exist_ok=True)
    table.write_parquet(output)
    return {
        "output": output,
        "chrom_filter": chrom,
        "limit": limit,
        "transcripts": table.height,
        "excluded_patch_transcripts": len(parsed) - len(primary_genes),
        "patch_only_genes": len(patch_only),
        "patch_only_gene_examples": patch_only[:10],
        "genes": table["gene_name"].n_unique(),
        "translation_match_rate": round(table["translation_matches"].mean(), 5),
        "translation_notes": dict(table["translation_note"].value_counts().iter_rows()),
        "mismatch_examples": mismatches[:MAX_REPORTED_MISMATCHES],
        "seconds": round(time.monotonic() - started, 1),
    }
