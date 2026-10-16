import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch

from healthcare.esm.store import AMINO_ACIDS, MODEL_ID, REVISION

MAX_RESIDUES = 1022
WINDOW_STEP = 511


@dataclass
class ProteinOutput:
    log_probs: np.ndarray
    embeddings: np.ndarray
    windows: list[tuple[int, int]]
    seconds: float


def windows(
    length: int, size: int = MAX_RESIDUES, step: int = WINDOW_STEP
) -> list[tuple[int, int]]:
    if length <= size:
        return [(0, length)]
    starts = list(range(0, length - size, step)) + [length - size]
    return [(start, start + size) for start in sorted(set(starts))]


def window_assignment(length: int, spans: list[tuple[int, int]]) -> np.ndarray:
    positions = np.arange(length)
    centers = np.array([(start + end - 1) / 2 for start, end in spans])
    inside = np.array([(positions >= start) & (positions < end) for start, end in spans])
    distance = np.where(inside, np.abs(positions[None, :] - centers[:, None]), np.inf)
    return distance.argmin(axis=0)


def log_likelihood_ratio(
    log_probs: np.ndarray, position: int, wild_type: str, mutant: str
) -> float:
    row = log_probs[position - 1]
    return float(row[AMINO_ACIDS.index(mutant)] - row[AMINO_ACIDS.index(wild_type)])


def entropy(log_probs: np.ndarray) -> np.ndarray:
    probabilities = np.exp(log_probs.astype(np.float32))
    return -(probabilities * log_probs.astype(np.float32)).sum(axis=1)


class Esm:
    def __init__(self, cache_dir: Path, device: str = "cuda"):
        from transformers import AutoTokenizer, EsmForMaskedLM

        self.device = device
        self.tokenizer = AutoTokenizer.from_pretrained(
            MODEL_ID, revision=REVISION, cache_dir=cache_dir
        )
        self.model = (
            EsmForMaskedLM.from_pretrained(
                MODEL_ID,
                revision=REVISION,
                cache_dir=cache_dir,
                use_safetensors=True,
                dtype=torch.float16,
            )
            .to(device)
            .eval()
        )
        self.amino_acid_ids = torch.tensor(
            self.tokenizer.convert_tokens_to_ids(list(AMINO_ACIDS)), device=device
        )

    @torch.inference_mode()
    def run_window(self, sequence: str) -> tuple[np.ndarray, np.ndarray]:
        tokens = self.tokenizer(sequence, return_tensors="pt").to(self.device)
        output = self.model(**tokens, output_hidden_states=True)
        logits = output.logits[0, 1:-1].float().index_select(-1, self.amino_acid_ids)
        log_probs = torch.log_softmax(logits, dim=-1)
        embeddings = output.hidden_states[-1][0, 1:-1]
        return log_probs.half().cpu().numpy(), embeddings.half().cpu().numpy()

    def run_protein(self, sequence: str) -> ProteinOutput:
        started = time.monotonic()
        sequence = sequence.replace("U", "C")
        spans = windows(len(sequence))
        results = [self.run_window(sequence[start:end]) for start, end in spans]
        chosen = window_assignment(len(sequence), spans)
        log_probs = np.empty((len(sequence), len(AMINO_ACIDS)), dtype=np.float16)
        embeddings = np.empty((len(sequence), results[0][1].shape[1]), dtype=np.float16)
        for i, (start, _) in enumerate(spans):
            mask = chosen == i
            log_probs[mask] = results[i][0][np.flatnonzero(mask) - start]
            embeddings[mask] = results[i][1][np.flatnonzero(mask) - start]
        if self.device == "cuda":
            torch.cuda.synchronize()
        return ProteinOutput(log_probs, embeddings, spans, time.monotonic() - started)
