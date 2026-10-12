import numpy as np
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    matthews_corrcoef,
    roc_auc_score,
)

ECE_BINS = 10


def expected_calibration_error(labels: np.ndarray, probabilities: np.ndarray) -> float:
    bins = np.minimum((probabilities * ECE_BINS).astype(int), ECE_BINS - 1)
    error = 0.0
    for b in range(ECE_BINS):
        mask = bins == b
        if mask.any():
            error += mask.mean() * abs(labels[mask].mean() - probabilities[mask].mean())
    return float(error)


def classification_metrics(
    labels: np.ndarray, probabilities: np.ndarray, threshold: float = 0.5
) -> dict:
    predicted = probabilities >= threshold
    return {
        "roc_auc": float(roc_auc_score(labels, probabilities)),
        "pr_auc": float(average_precision_score(labels, probabilities)),
        "mcc": float(matthews_corrcoef(labels, predicted)),
        "balanced_accuracy": float(balanced_accuracy_score(labels, predicted)),
        "brier": float(brier_score_loss(labels, probabilities)),
        "ece": expected_calibration_error(labels, probabilities),
        "threshold": threshold,
        "n": int(labels.size),
        "positives": int(labels.sum()),
    }


def gene_bootstrap(
    labels: np.ndarray,
    probabilities: np.ndarray,
    genes: np.ndarray,
    rounds: int = 200,
    seed: int = 0,
) -> dict:
    rng = np.random.default_rng(seed)
    unique, inverse = np.unique(genes, return_inverse=True)
    members = [np.flatnonzero(inverse == i) for i in range(unique.size)]
    values = {"roc_auc": [], "pr_auc": []}
    for _ in range(rounds):
        index = np.concatenate([members[i] for i in rng.integers(0, unique.size, unique.size)])
        y, p = labels[index], probabilities[index]
        if 0 < y.sum() < y.size:
            values["roc_auc"].append(roc_auc_score(y, p))
            values["pr_auc"].append(average_precision_score(y, p))
    return {
        name: [round(float(np.percentile(v, 2.5)), 4), round(float(np.percentile(v, 97.5)), 4)]
        for name, v in values.items()
    }


def rounded(metrics: dict, digits: int = 4) -> dict:
    return {k: round(v, digits) if isinstance(v, float) else v for k, v in metrics.items()}
