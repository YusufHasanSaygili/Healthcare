import gzip
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import polars as pl

from healthcare.data.tabix import read_index

KEY = ["chrom", "pos", "ref", "alt"]
TOOLS = {
    "SIFT_converted_rankscore": "SIFT",
    "SIFT4G_converted_rankscore": "SIFT4G",
    "Polyphen2_HDIV_rankscore": "PolyPhen-2 HDIV",
    "Polyphen2_HVAR_rankscore": "PolyPhen-2 HVAR",
    "MutationTaster_rankscore": "MutationTaster",
    "MutationAssessor_rankscore": "MutationAssessor",
    "PROVEAN_converted_rankscore": "PROVEAN",
    "VEST4_rankscore": "VEST4",
    "MetaSVM_rankscore": "MetaSVM",
    "MetaLR_rankscore": "MetaLR",
    "MetaRNN_rankscore": "MetaRNN",
    "M-CAP_rankscore": "M-CAP",
    "REVEL_rankscore": "REVEL",
    "MutPred2_rankscore": "MutPred2",
    "MVP_rankscore": "MVP",
    "gMVP_rankscore": "gMVP",
    "MisFit_D_rankscore": "MisFit",
    "MPC_rankscore": "MPC",
    "PrimateAI_rankscore": "PrimateAI",
    "DEOGEN2_rankscore": "DEOGEN2",
    "BayesDel_addAF_rankscore": "BayesDel (with AF)",
    "BayesDel_noAF_rankscore": "BayesDel (no AF)",
    "ClinPred_rankscore": "ClinPred",
    "LIST-S2_rankscore": "LIST-S2",
    "VARITY_R_rankscore": "VARITY_R",
    "VARITY_ER_rankscore": "VARITY_ER",
    "ESM1b_converted_rankscore": "ESM1b",
    "AlphaMissense_rankscore": "AlphaMissense",
    "PHACTboost_rankscore": "PHACTboost",
    "MutFormer_rankscore": "MutFormer",
    "MutScore_rankscore": "MutScore",
    "popEVE_converted_rankscore": "popEVE",
    "CADD_raw_rankscore": "CADD",
    "DANN_rankscore": "DANN",
    "fathmm-XF_coding_rankscore": "FATHMM-XF",
    "Eigen-raw_coding_rankscore": "Eigen",
    "GPN_MSA_converted_rankscore": "GPN-MSA",
}


def header(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        return handle.readline().rstrip("\n").lstrip("#").split("\t")


def chromosome_starts(index_path: Path) -> dict[str, int]:
    return {
        name: min(offset for offset in reference.linear if offset) if any(reference.linear) else 0
        for name, reference in read_index(index_path).items()
    }


def scan(path: Path, chrom: str, start: int, wanted: set, columns: list[int]) -> list[list]:
    found = []
    with path.open("rb") as raw:
        raw.seek(start >> 16)
        with gzip.GzipFile(fileobj=raw) as stream:
            stream.read(start & 0xFFFF)
            for line in stream:
                fields = line.split(b"\t", 4)
                if fields[0].decode() != chrom:
                    break
                key = (int(fields[1]), fields[2].decode(), fields[3].decode())
                if key in wanted:
                    values = line.rstrip(b"\r\n").split(b"\t")
                    found.append([chrom, *key, *(values[i].decode() for i in columns)])
    return found


def as_score(column: str) -> pl.Expr:
    return pl.col(column).replace(".", None).cast(pl.Float64, strict=False)


def extract(path: Path, index_path: Path, keys: pl.DataFrame, workers: int = 8) -> pl.DataFrame:
    names = header(path)
    columns = [names.index(column) for column in TOOLS]
    starts = chromosome_starts(index_path)
    jobs = []
    for (chrom,), group in keys.group_by("chrom"):
        if chrom in starts:
            wanted = set(group.select("pos", "ref", "alt").iter_rows())
            jobs.append((path, chrom, starts[chrom], wanted, columns))
    with ProcessPoolExecutor(workers) as pool:
        rows = [row for result in pool.map(scan, *zip(*jobs, strict=True)) for row in result]
    table = pl.DataFrame(rows, schema=[*KEY, *TOOLS], orient="row")
    return table.with_columns(pl.col("pos").cast(pl.Int64), *(as_score(c) for c in TOOLS))
