import time
from pathlib import Path

import numpy as np
import polars as pl

from healthcare.esm import store

BLOCKS = ("position", "window3", "window7", "protein_mean")
KEY = ["chrom", "pos", "ref", "alt"]


def variant_blocks(embeddings: np.ndarray, position: int) -> np.ndarray:
    index = position - 1
    length = embeddings.shape[0]
    window3 = embeddings[max(0, index - 3) : min(length, index + 4)]
    window7 = embeddings[max(0, index - 7) : min(length, index + 8)]
    return np.stack(
        [
            embeddings[index],
            window3.mean(axis=0, dtype=np.float32),
            window7.mean(axis=0, dtype=np.float32),
            embeddings.mean(axis=0, dtype=np.float32),
        ]
    ).astype(np.float16)


def build(
    table: pl.DataFrame, store_dir: Path, output: Path, columns: list[int] | None = None
) -> dict:
    started = time.monotonic()
    kept = columns or list(range(len(BLOCKS)))
    output.parent.mkdir(parents=True, exist_ok=True)
    blocks = np.lib.format.open_memmap(
        output, mode="w+", dtype=np.float16, shape=(table.height, len(kept), 1280)
    )
    missing = 0
    indexed = table.select("transcript_id", "protein_position").with_row_index()
    for (transcript_id,), group in indexed.group_by("transcript_id"):
        if not store.is_done(store_dir, transcript_id):
            missing += group.height
            continue
        _, embeddings = store.load(store_dir, transcript_id)
        embeddings = np.asarray(embeddings)
        for row, position in group.select("index", "protein_position").iter_rows():
            if position is not None and 1 <= position <= embeddings.shape[0]:
                blocks[row] = variant_blocks(embeddings, position)[kept]
            else:
                missing += 1
    blocks.flush()
    del blocks
    table.select(KEY).write_parquet(output.with_suffix(".keys.parquet"))
    return {
        "rows": table.height,
        "missing": missing,
        "gigabytes": round(output.stat().st_size / 1e9, 2),
        "seconds": round(time.monotonic() - started, 1),
    }


def load(output: Path, table: pl.DataFrame, columns: list[int] | None = None) -> np.ndarray:
    keys = pl.read_parquet(output.with_suffix(".keys.parquet")).with_row_index("row")
    rows = table.select(KEY).join(keys, on=KEY, how="left", maintain_order="left")["row"]
    if rows.null_count():
        raise RuntimeError(f"{rows.null_count()} variants have no precomputed ESM blocks")
    blocks = np.load(output, mmap_mode="r")
    if columns is None:
        return blocks[rows.to_numpy()]
    return blocks[rows.to_numpy()[:, None], np.asarray(columns)]
