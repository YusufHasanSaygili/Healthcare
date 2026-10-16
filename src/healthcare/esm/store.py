import time
from pathlib import Path

import numpy as np
import polars as pl

MODEL_ID = "facebook/esm2_t33_650M_UR50D"
REVISION = "08e4846e537177426273712802403f7ba8261b6c"
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
KINDS = ("log_probs", "embeddings")


def paths(store_dir: Path, transcript_id: str) -> dict[str, Path]:
    return {kind: store_dir / kind / f"{transcript_id}.npy" for kind in KINDS}


def is_done(store_dir: Path, transcript_id: str) -> bool:
    return all(path.exists() for path in paths(store_dir, transcript_id).values())


def save(
    store_dir: Path, transcript_id: str, log_probs: np.ndarray, embeddings: np.ndarray
) -> None:
    for kind, array in zip(KINDS, (log_probs, embeddings), strict=True):
        target = paths(store_dir, transcript_id)[kind]
        target.parent.mkdir(parents=True, exist_ok=True)
        part = target.with_name(target.name + ".part")
        with part.open("wb") as handle:
            np.save(handle, array)
        part.replace(target)


def load(store_dir: Path, transcript_id: str) -> tuple[np.ndarray, np.ndarray]:
    files = paths(store_dir, transcript_id)
    return tuple(np.load(files[kind], mmap_mode="r") for kind in KINDS)


def build(transcripts_path: Path, store_dir: Path, cache_dir: Path, limit: int | None) -> dict:
    from healthcare.esm.model import Esm

    started = time.monotonic()
    proteins = pl.read_parquet(transcripts_path).filter("translation_matches")
    proteins = proteins.select("transcript_id", "protein").sort(pl.col("protein").str.len_chars())
    if limit:
        proteins = proteins.head(limit)
    todo = [(t, p) for t, p in proteins.iter_rows() if not is_done(store_dir, t)]
    esm = Esm(cache_dir) if todo else None
    errors = []
    computed = 0
    for transcript_id, protein in todo:
        try:
            output = esm.run_protein(protein)
        except RuntimeError as error:
            errors.append(
                {"transcript_id": transcript_id, "length": len(protein), "error": str(error)[:200]}
            )
            continue
        if output.log_probs.shape[0] != len(protein):
            errors.append({"transcript_id": transcript_id, "error": "length mismatch"})
            continue
        save(store_dir, transcript_id, output.log_probs, output.embeddings)
        computed += 1
    done = sum(is_done(store_dir, t) for t in proteins["transcript_id"])
    sizes = {
        kind: sum(f.stat().st_size for f in (store_dir / kind).glob("*.npy")) for kind in KINDS
    }
    return {
        "model": MODEL_ID,
        "revision": REVISION,
        "proteins_requested": proteins.height,
        "computed_now": computed,
        "proteins_done": done,
        "errors": len(errors),
        "error_examples": errors[:10],
        "residues_total": int(proteins["protein"].str.len_chars().sum()),
        "gigabytes": {kind: round(size / 1e9, 2) for kind, size in sizes.items()},
        "seconds": round(time.monotonic() - started, 1),
    }
