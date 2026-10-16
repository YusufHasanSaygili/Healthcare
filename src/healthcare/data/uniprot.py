import gzip
import re
import time
from collections import Counter
from collections.abc import Iterable, Iterator
from pathlib import Path

import polars as pl

KEPT_TYPES = {
    "DOMAIN",
    "REPEAT",
    "ZN_FING",
    "DNA_BIND",
    "MOTIF",
    "REGION",
    "COILED",
    "COMPBIAS",
    "ACT_SITE",
    "BINDING",
    "SITE",
    "MOD_RES",
    "CARBOHYD",
    "LIPID",
    "CROSSLNK",
    "DISULFID",
    "TRANSMEM",
    "INTRAMEM",
    "TOPO_DOM",
    "SIGNAL",
    "TRANSIT",
    "PROPEP",
    "PEPTIDE",
    "HELIX",
    "STRAND",
    "TURN",
}
EXCLUDED_TYPES = {
    "VARIANT",
    "MUTAGEN",
    "CONFLICT",
    "VAR_SEQ",
    "NON_STD",
    "NON_TER",
    "INIT_MET",
    "CHAIN",
}
AUTOMATIC_EVIDENCE = {"ECO:0000256", "ECO:0000259", "ECO:0000313", "ECO:0000213"}
LOCATION = re.compile(r"^[<?]?(\d+)(?:\.\.[>?]?(\d+))?$")
EVIDENCE = re.compile(r"ECO:\d{7}")
NOTE = re.compile(r'/note="([^"]*)"')
CITATION = re.compile(r"\s*\((?:PubMed|Ref\.)[^)]*\)|\s*\{ECO:[^}]*\}")
SENTENCE_END = re.compile(r"(?<=\.)\s+(?=[A-Z])")
FUNCTION_SENTENCES = 2


def is_curated(evidence: set[str]) -> bool:
    return not evidence or not evidence <= AUTOMATIC_EVIDENCE


def feature_rows(feature_type: str, location: str, qualifiers: str) -> list[tuple]:
    match = LOCATION.match(location)
    evidence = set(EVIDENCE.findall(qualifiers))
    if feature_type not in KEPT_TYPES or not match or not is_curated(evidence):
        return []
    start = int(match.group(1))
    end = int(match.group(2) or start)
    note_match = NOTE.search(qualifiers)
    note = note_match.group(1) if note_match else None
    codes = ",".join(sorted(evidence))
    if feature_type == "DISULFID" and end != start:
        return [(feature_type, start, start, note, codes), (feature_type, end, end, note, codes)]
    return [(feature_type, start, end, note, codes)]


def clean_function(text: str) -> str | None:
    plain = " ".join(CITATION.sub("", text).split())
    if not plain:
        return None
    return " ".join(SENTENCE_END.split(plain)[:FUNCTION_SENTENCES])


def parse_entries(lines: Iterable[str]) -> Iterator[dict]:
    accession, gene, sequence, features = None, None, [], []
    current: list[str] | None = None
    function: list[str] = []
    in_function = False

    def flush_feature():
        if current:
            features.extend(feature_rows(current[0], current[1], " ".join(current[2:])))

    for line in lines:
        tag = line[:2]
        if tag == "AC" and accession is None:
            accession = line[5:].split(";")[0].strip()
        elif tag == "GN" and gene is None and "Name=" in line:
            gene = line.split("Name=")[1].split(";")[0].split(" {")[0].strip()
        elif tag == "CC":
            body = line[5:].rstrip()
            if body.startswith("-!-"):
                in_function = body.startswith("-!- FUNCTION:") and not function
                if in_function:
                    function.append(body.removeprefix("-!- FUNCTION:"))
            elif in_function and not body.startswith("---"):
                function.append(body)
            else:
                in_function = False
        elif tag == "FT":
            if line[5] != " ":
                flush_feature()
                parts = line[5:].split()
                current = [parts[0], parts[1] if len(parts) > 1 else ""]
            elif current is not None:
                current.append(line[21:].strip())
        elif tag == "  ":
            sequence.append(line.replace(" ", "").strip())
        elif tag == "//":
            flush_feature()
            yield {
                "accession": accession,
                "gene": gene,
                "sequence": "".join(sequence),
                "features": features,
                "function": clean_function(" ".join(function)),
            }
            accession, gene, sequence, features, current = None, None, [], [], None
            function, in_function = [], False


def build_functions(dat_gz: Path, transcripts: pl.DataFrame, output: Path) -> dict:
    proteins = dict(transcripts.select("protein", "transcript_id").iter_rows())
    genes = dict(transcripts.select("transcript_id", "gene_name").iter_rows())
    rows = []
    with gzip.open(dat_gz, "rt", encoding="utf-8") as handle:
        for entry in parse_entries(handle):
            transcript_id = proteins.get(entry["sequence"])
            if transcript_id is not None and entry["function"]:
                rows.append((transcript_id, genes[transcript_id], entry["function"]))
    table = pl.DataFrame(
        rows,
        schema={"transcript_id": pl.String, "gene_name": pl.String, "function": pl.String},
        orient="row",
    ).sort("transcript_id")
    output.parent.mkdir(parents=True, exist_ok=True)
    table.write_parquet(output)
    return {"proteins_with_function": table.height}


def build(dat_gz: Path, transcripts: pl.DataFrame, output: Path) -> dict:
    started = time.monotonic()
    proteins = dict(transcripts.select("protein", "transcript_id").iter_rows())
    genes = dict(transcripts.select("transcript_id", "gene_name").iter_rows())
    rows, matched, entries = [], set(), 0
    with gzip.open(dat_gz, "rt", encoding="utf-8") as handle:
        for entry in parse_entries(handle):
            entries += 1
            transcript_id = proteins.get(entry["sequence"])
            if transcript_id is None:
                continue
            matched.add(transcript_id)
            rows.extend(
                (transcript_id, genes[transcript_id], entry["accession"], *feature)
                for feature in entry["features"]
            )
    table = pl.DataFrame(
        rows,
        schema={
            "transcript_id": pl.String,
            "gene_name": pl.String,
            "uniprot_accession": pl.String,
            "feature_type": pl.String,
            "start": pl.Int64,
            "end": pl.Int64,
            "note": pl.String,
            "evidence": pl.String,
        },
        orient="row",
    ).sort("transcript_id", "start", "feature_type")
    output.parent.mkdir(parents=True, exist_ok=True)
    table.write_parquet(output)
    assert not set(table["feature_type"].unique()) & EXCLUDED_TYPES
    return {
        "version": "UniProtKB/Swiss-Prot 2026_03 (human)",
        "swissprot_entries": entries,
        "mane_proteins": transcripts.height,
        "mane_proteins_matched": len(matched),
        "match_rate": round(len(matched) / transcripts.height, 4),
        "features": table.height,
        "features_by_type": dict(Counter(table["feature_type"].to_list()).most_common()),
        "seconds": round(time.monotonic() - started, 1),
    }
