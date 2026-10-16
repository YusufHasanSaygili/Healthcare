from pathlib import Path

import numpy as np
import polars as pl

from healthcare.esm import store
from healthcare.esm.store import AMINO_ACIDS

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


SPAN_FEATURES = [
    "f_esm_removed_log_prob",
    "f_esm_inserted_log_prob",
    "f_esm_span_llr",
    "f_esm_site_entropy",
]
EMPTY_SPAN = (None,) * len(SPAN_FEATURES)


def span_features(log_probs: np.ndarray, position: int, removed: str, inserted: str) -> tuple:
    length = log_probs.shape[0]
    if not 2 <= position <= length + 1:
        return EMPTY_SPAN
    rows = np.asarray(log_probs, dtype=np.float32)
    site = range(position - 1, min(position - 1 + max(len(removed), 1), length))
    flank = [position - 2, *site] if not removed else list(site)
    entropy = float(np.mean([-(np.exp(rows[i]) * rows[i]).sum() for i in flank]))
    removed_scores = [
        rows[position - 1 + i, INDEX[aa]]
        for i, aa in enumerate(removed)
        if aa in INDEX and position - 1 + i < length
    ]
    anchors = list(site) if removed else [position - 2]
    inserted_scores = [
        rows[anchors[min(i, len(anchors) - 1)], INDEX[aa]]
        for i, aa in enumerate(inserted)
        if aa in INDEX
    ]
    removed_mean = float(np.mean(removed_scores)) if removed_scores else None
    inserted_mean = float(np.mean(inserted_scores)) if inserted_scores else None
    llr = (
        inserted_mean - removed_mean
        if removed_mean is not None and inserted_mean is not None
        else None
    )
    return (removed_mean, inserted_mean, llr, entropy)
