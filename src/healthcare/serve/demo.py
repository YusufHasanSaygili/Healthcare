import hashlib
import json
import shutil
import time
import zipfile
from dataclasses import dataclass
from pathlib import Path

import duckdb
import polars as pl
from pyfaidx import Fasta

from healthcare.esm import store
from healthcare.features import missense
from healthcare.serve.pipeline import BACKGROUND_SEED, Predictor

DEFAULT_GENES = (
    "APC",
    "BRCA1",
    "BRCA2",
    "CFTR",
    "FBN1",
    "KCNQ1",
    "LDLR",
    "MLH1",
    "MSH2",
    "MYH7",
    "PTEN",
    "RYR1",
    "SCN5A",
    "TP53",
)
PADDING = 2000
BACKGROUND_SIZE = 32
MODEL_NAME = "all_types_v2.onnx"


class OutsideDemoPack(ValueError):
    pass


@dataclass(frozen=True)
class Sequence:
    seq: str


class RegionSequence:
    def __init__(self, chrom: str, length: int, regions: list[tuple[int, str]], genes: str):
        self.chrom = chrom
        self.length = length
        self.regions = regions
        self.genes = genes

    def __len__(self) -> int:
        return self.length

    def __getitem__(self, key: slice) -> Sequence:
        for start, sequence in self.regions:
            if start <= key.start and key.stop <= start + len(sequence):
                return Sequence(sequence[key.start - start : key.stop - start])
        raise OutsideDemoPack(
            f"{self.chrom}:{key.start + 1} is outside the demo pack, which covers {self.genes}"
        )


class RegionGenome:
    def __init__(self, path: Path, lengths: dict[str, int], genes: list[str]):
        table = pl.read_parquet(path).sort("chrom", "start")
        listed = ", ".join(genes)
        self.chromosomes = {
            chrom: RegionSequence(
                chrom,
                length,
                [
                    (start, sequence)
                    for start, sequence in table.filter(pl.col("chrom") == chrom)
                    .select("start", "sequence")
                    .iter_rows()
                ],
                listed,
            )
            for chrom, length in lengths.items()
        }

    def __getitem__(self, chrom: str) -> RegionSequence:
        return self.chromosomes[chrom]


def spans(transcripts: pl.DataFrame) -> pl.DataFrame:
    return transcripts.select(
        "gene_name",
        "transcript_id",
        "chrom",
        (pl.col("exon_starts").list.min() - PADDING).clip(1).alias("start"),
        (pl.col("exon_ends").list.max() + PADDING).alias("end"),
    )


def merge_regions(regions: pl.DataFrame) -> list[tuple[str, int, int]]:
    merged: list[tuple[str, int, int]] = []
    for chrom, start, end in regions.sort("chrom", "start").select("chrom", "start", "end").rows():
        if merged and merged[-1][0] == chrom and start <= merged[-1][2] + 1:
            merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append((chrom, start, end))
    return merged


def pack_transcripts(transcripts: pl.DataFrame, genes: list[str]) -> pl.DataFrame:
    unknown = sorted(set(genes) - set(transcripts["gene_name"]))
    if unknown:
        raise ValueError(f"not MANE Select genes: {', '.join(unknown)}")
    all_spans = spans(transcripts)
    chosen = all_spans.filter(pl.col("gene_name").is_in(genes))
    while True:
        overlapping = all_spans.join(chosen.select("chrom", "start", "end"), on="chrom").filter(
            (pl.col("start") <= pl.col("end_right")) & (pl.col("end") >= pl.col("start_right"))
        )
        grown = all_spans.filter(
            pl.col("transcript_id").is_in(overlapping["transcript_id"].to_list())
        )
        if grown.height == chosen.height:
            return grown.sort("chrom", "start")
        chosen = grown


def region_filter(regions: list[tuple[str, int, int]]) -> str:
    return " OR ".join(
        f"(chrom = '{chrom}' AND pos BETWEEN {start} AND {end})" for chrom, start, end in regions
    )


def copy_rows(source: str, target: Path, regions: list[tuple[str, int, int]]) -> int:
    table = duckdb.sql(
        f"SELECT * FROM read_parquet('{source}') WHERE {region_filter(regions)} ORDER BY chrom, pos"
    ).pl()
    target.parent.mkdir(parents=True, exist_ok=True)
    table.write_parquet(target)
    return table.height


def background_rows(train: pl.DataFrame, feature_names: list[str]) -> pl.DataFrame:
    return pl.concat(
        group.sample(min(BACKGROUND_SIZE, group.height), seed=BACKGROUND_SEED)
        for _, group in train.select("consequence", *feature_names).group_by("consequence")
    )


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def layout(root: Path) -> dict[str, Path]:
    data = root / "data"
    return {
        "transcripts": data / "transcripts.parquet",
        "gnomad_exomes": data / "gnomad_exomes",
        "coverage": data / "coverage.parquet",
        "conservation": data / "conservation",
        "constraint": data / "constraint.parquet",
        "uniprot": data / "uniprot_features.parquet",
        "functions": data / "uniprot_functions.parquet",
        "genome": data / "genome_regions.parquet",
        "background": data / "background.parquet",
        "esm_store": root / "esm",
        "model": root / MODEL_NAME,
        "manifest": root / "manifest.json",
    }


def build(
    genes: list[str],
    inputs: missense.Inputs,
    functions_path: Path,
    genome_path: Path,
    model_path: Path,
    train: pl.DataFrame,
    feature_names: list[str],
    out_dir: Path,
    zip_path: Path,
) -> dict:
    started = time.monotonic()
    if out_dir.exists():
        shutil.rmtree(out_dir)
    paths = layout(out_dir)
    paths["esm_store"].mkdir(parents=True)
    paths["transcripts"].parent.mkdir(parents=True)

    transcripts = pl.read_parquet(inputs.transcripts)
    packed = pack_transcripts(transcripts, genes)
    regions = merge_regions(packed)
    transcript_ids = packed["transcript_id"].to_list()
    by_id = pl.col("transcript_id").is_in(transcript_ids)

    shutil.copyfile(inputs.transcripts, paths["transcripts"])
    gnomad_rows = copy_rows(
        (inputs.gnomad_exomes / "*.parquet").as_posix(),
        paths["gnomad_exomes"] / "variants.parquet",
        regions,
    )
    coverage_rows = copy_rows(inputs.coverage.as_posix(), paths["coverage"], regions)
    conservation_rows = copy_rows(
        (inputs.conservation / "*.parquet").as_posix(),
        paths["conservation"] / "positions.parquet",
        regions,
    )
    unversioned = [transcript_id.split(".")[0] for transcript_id in transcript_ids]
    pl.read_parquet(inputs.constraint).filter(
        pl.col("transcript_id").is_in(unversioned)
    ).write_parquet(paths["constraint"])
    pl.read_parquet(inputs.uniprot).filter(by_id).write_parquet(paths["uniprot"])
    pl.read_parquet(functions_path).filter(by_id).write_parquet(paths["functions"])

    reference = Fasta(str(genome_path))
    pl.DataFrame(
        [
            (chrom, start - 1, reference[chrom][start - 1 : end].seq.upper())
            for chrom, start, end in regions
        ],
        schema={"chrom": pl.String, "start": pl.Int64, "sequence": pl.String},
        orient="row",
    ).write_parquet(paths["genome"])
    lengths = {name: len(record) for name, record in reference.records.items()}

    esm_proteins = 0
    for transcript_id in transcript_ids:
        if not store.is_done(inputs.esm_store, transcript_id):
            continue
        for kind, source in store.paths(inputs.esm_store, transcript_id).items():
            target = store.paths(paths["esm_store"], transcript_id)[kind]
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
        esm_proteins += 1

    background_rows(train, feature_names).write_parquet(paths["background"])
    shutil.copyfile(model_path, paths["model"])

    manifest = {
        "genes": sorted(genes),
        "pack_genes": sorted(set(packed["gene_name"])),
        "transcripts": len(transcript_ids),
        "regions": [{"chrom": c, "start": s, "end": e} for c, s, e in regions],
        "chrom_lengths": lengths,
        "padding_bp": PADDING,
        "model_sha256": file_sha256(model_path),
    }
    paths["manifest"].write_text(json.dumps(manifest, indent=2))

    zip_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(out_dir.rglob("*")):
            if path.is_file():
                archive.write(path, Path("demo") / path.relative_to(out_dir))

    return {
        "genes": manifest["genes"],
        "pack_genes": manifest["pack_genes"],
        "transcripts": len(transcript_ids),
        "regions": len(regions),
        "region_bp": sum(e - s + 1 for _, s, e in regions),
        "gnomad_rows": gnomad_rows,
        "coverage_rows": coverage_rows,
        "conservation_rows": conservation_rows,
        "esm_proteins": esm_proteins,
        "background_rows": pl.read_parquet(paths["background"]).height,
        "pack_bytes": sum(p.stat().st_size for p in out_dir.rglob("*") if p.is_file()),
        "zip_path": str(zip_path),
        "zip_bytes": zip_path.stat().st_size,
        "zip_sha256": file_sha256(zip_path),
        "seconds": round(time.monotonic() - started, 1),
    }


def load(root: Path) -> Predictor:
    if not root.exists():
        raise FileNotFoundError(
            f"no demo pack at {root}: download healthcare-demo-v2.zip from the GitHub release "
            "and unzip it in the repository root"
        )
    paths = layout(root)
    missing = [str(p) for p in paths.values() if not p.exists()]
    if missing:
        raise FileNotFoundError(f"incomplete demo pack, missing: {', '.join(missing)}")
    manifest = json.loads(paths["manifest"].read_text())
    inputs = missense.Inputs(
        transcripts=paths["transcripts"],
        gnomad_exomes=paths["gnomad_exomes"],
        coverage=paths["coverage"],
        conservation=paths["conservation"],
        constraint=paths["constraint"],
        uniprot=paths["uniprot"],
        esm_store=paths["esm_store"],
    )
    functions = dict(
        pl.read_parquet(paths["functions"]).select("transcript_id", "function").iter_rows()
    )
    genome = RegionGenome(paths["genome"], manifest["chrom_lengths"], manifest["genes"])
    return Predictor(
        inputs, paths["model"], pl.read_parquet(paths["background"]), genome, functions
    )


def same(a, b) -> bool:
    return a == b or (isinstance(a, float) and isinstance(b, float) and a != a and b != b)


def parity(pack: Predictor, full: Predictor, variants: list[tuple[str, int, str, str]]) -> dict:
    largest, feature_mismatches, decision_mismatches, types = 0.0, 0, 0, set()
    for small, reference in zip(
        pack.predict_many(variants), full.predict_many(variants), strict=True
    ):
        if type(small) is not type(reference):
            decision_mismatches += 1
            continue
        if not hasattr(small, "probability"):
            continue
        types.add(small.consequence)
        largest = max(largest, abs(small.probability - reference.probability))
        decision_mismatches += small.decision != reference.decision
        feature_mismatches += any(
            not same(value, reference.features[name]) for name, value in small.features.items()
        )
    return {
        "variants": len(variants),
        "types": sorted(types),
        "max_probability_diff": largest,
        "feature_mismatches": feature_mismatches,
        "decision_mismatches": decision_mismatches,
    }
