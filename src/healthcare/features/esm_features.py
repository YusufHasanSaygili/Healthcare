from pathlib import Path

import numpy as np
import polars as pl

from healthcare.esm import store
from healthcare.esm.model import AMINO_ACIDS

ESM_FEATURES = [
    "f_esm_llr",
    "f_esm_wt_log_prob",
    "f_esm_mut_log_prob",
    "f_esm_entropy",
    "f_esm_mut_rank",
]
INDEX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}
EMPTY = (None,) * len(ESM_FEATURES)


def position_features(log_probs: np.ndarray, position: int, ref: str, alt: str) -> tuple:
    if ref not in INDEX or alt not in INDEX or not 1 <= position <= log_probs.shape[0]:
        return EMPTY
    row = np.asarray(log_probs[position - 1], dtype=np.float32)
    wild_type, mutant = row[INDEX[ref]], row[INDEX[alt]]
    entropy = float(-(np.exp(row) * row).sum())
    rank = int((row > mutant).sum()) + 1
    return (float(mutant - wild_type), float(wild_type), float(mutant), entropy, rank)


def esm_features(table: pl.DataFrame, store_dir: Path) -> pl.DataFrame:
    values: list[tuple] = [EMPTY] * table.height
    indexed = table.select("transcript_id", "protein_position", "aa_ref", "aa_alt").with_row_index()
    for (transcript_id,), group in indexed.group_by("transcript_id"):
        if not store.is_done(store_dir, transcript_id):
            continue
        log_probs, _ = store.load(store_dir, transcript_id)
        for row, position, ref, alt in group.select(
            "index", "protein_position", "aa_ref", "aa_alt"
        ).iter_rows():
            values[row] = position_features(log_probs, position, ref, alt)
    features = pl.DataFrame(
        values, schema={name: pl.Float64 for name in ESM_FEATURES}, orient="row"
    )
    return table.hstack(features)
