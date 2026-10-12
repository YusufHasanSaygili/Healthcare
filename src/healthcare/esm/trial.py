import time
from pathlib import Path

import numpy as np
import polars as pl
import torch

from healthcare.annotate.snv import THREE_LETTER
from healthcare.esm.model import (
    AMINO_ACIDS,
    MODEL_ID,
    REVISION,
    Esm,
    entropy,
    log_likelihood_ratio,
)
from healthcare.split.genes import load_split

THREE_TO_ONE = {three: one for one, three in THREE_LETTER.items()}


def labeled_llr(output, protein: str, labeled: pl.DataFrame) -> dict:
    standard = list(AMINO_ACIDS)
    rows = labeled.filter(pl.col("aa_ref").is_in(standard) & pl.col("aa_alt").is_in(standard))
    values = {0: [], 1: []}
    for position, ref, alt, label in rows.select(
        "protein_position", "aa_ref", "aa_alt", "label"
    ).iter_rows():
        if protein[position - 1] == ref:
            values[label].append(log_likelihood_ratio(output.log_probs, position, ref, alt))
    return {
        "benign_n": len(values[0]),
        "benign_median_llr": round(float(np.median(values[0])), 3) if values[0] else None,
        "pathogenic_n": len(values[1]),
        "pathogenic_median_llr": round(float(np.median(values[1])), 3) if values[1] else None,
    }


def run(
    transcripts_path: Path,
    splits_dir: Path,
    cache_dir: Path,
    gene: str,
    variant: tuple[str, int, str] | None,
) -> dict:
    started = time.monotonic()
    transcripts = pl.read_parquet(transcripts_path)
    protein = transcripts.filter(pl.col("gene_name") == gene)["protein"][0]
    longest = transcripts.sort(pl.col("protein").str.len_chars(), descending=True).row(
        0, named=True
    )

    torch.cuda.reset_peak_memory_stats()
    load_started = time.monotonic()
    esm = Esm(cache_dir)
    load_seconds = time.monotonic() - load_started
    esm.run_protein("MKTAYIAKQR")
    output = esm.run_protein(protein)
    long_output = esm.run_protein(longest["protein"])

    open_rows = pl.concat([load_split(splits_dir, "train"), load_split(splits_dir, "validation")])
    gene_rows = open_rows.filter(pl.col("gene_name") == gene)
    summary = {
        "model": MODEL_ID,
        "revision": REVISION,
        "gpu": torch.cuda.get_device_name(0),
        "load_seconds": round(load_seconds, 1),
        "gene": gene,
        "protein_length": len(protein),
        "windows": output.windows,
        "seconds": round(output.seconds, 2),
        "log_probs_shape": list(output.log_probs.shape),
        "embeddings_shape": list(output.embeddings.shape),
        "mean_position_entropy": round(float(entropy(output.log_probs).mean()), 3),
        "longest_protein": {
            "gene": longest["gene_name"],
            "length": len(longest["protein"]),
            "windows": len(long_output.windows),
            "seconds": round(long_output.seconds, 2),
        },
        "peak_gpu_memory_gib": round(torch.cuda.max_memory_allocated() / 1024**3, 2),
        "labeled_variants_in_open_splits": labeled_llr(output, protein, gene_rows)
        if gene_rows.height
        else "gene is in the locked test split; skipped",
        "total_seconds": round(time.monotonic() - started, 1),
    }
    if variant:
        ref, position, alt = variant
        summary["variant"] = {
            "change": f"{ref}{position}{alt}",
            "llr": round(log_likelihood_ratio(output.log_probs, position, ref, alt), 3),
            "wild_type_probability": round(
                float(np.exp(output.log_probs[position - 1][AMINO_ACIDS.index(ref)])), 4
            ),
        }
    return summary


def parse_protein_change(text: str) -> tuple[str, int, str]:
    body = text.removeprefix("p.")
    return THREE_TO_ONE[body[:3]], int(body[3:-3]), THREE_TO_ONE[body[-3:]]
