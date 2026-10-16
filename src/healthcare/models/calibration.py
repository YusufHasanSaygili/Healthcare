import json
import random
import time
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
import torch
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression

from healthcare.esm import variant_inputs
from healthcare.evaluate.metrics import classification_metrics, rounded
from healthcare.features.missense import to_matrix
from healthcare.models import dual_tower
from healthcare.split.genes import load_split

SEED = 20261013
FOLDS = 5
THRESHOLDS = {"likely_benign_max": 0.10, "likely_pathogenic_min": 0.90}
EPSILON = 1e-6
DECISIONS = ("likely_benign", "uncertain", "likely_pathogenic")


def logit(probabilities: np.ndarray) -> np.ndarray:
    clipped = np.clip(probabilities, EPSILON, 1 - EPSILON)
    return np.log(clipped / (1 - clipped))


def fit_calibrator(method: str, probabilities: np.ndarray, labels: np.ndarray) -> dict:
    if method == "platt":
        model = LogisticRegression(C=1e6).fit(logit(probabilities)[:, None], labels)
        return {
            "method": method,
            "slope": float(model.coef_[0, 0]),
            "bias": float(model.intercept_[0]),
        }
    model = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
    model.fit(probabilities, labels)
    return {
        "method": method,
        "x": model.X_thresholds_.tolist(),
        "y": model.y_thresholds_.tolist(),
    }


def calibrate(calibrator: dict, probabilities: np.ndarray) -> np.ndarray:
    if calibrator["method"] == "platt":
        z = calibrator["slope"] * logit(probabilities) + calibrator["bias"]
        return 1 / (1 + np.exp(-z))
    return np.interp(probabilities, calibrator["x"], calibrator["y"])


def decide(probabilities: np.ndarray) -> np.ndarray:
    return np.where(
        probabilities >= THRESHOLDS["likely_pathogenic_min"],
        2,
        np.where(probabilities <= THRESHOLDS["likely_benign_max"], 0, 1),
    )


def decision_summary(probabilities: np.ndarray, labels: np.ndarray) -> dict:
    decisions = decide(probabilities)
    summary = {}
    for code, name in enumerate(DECISIONS):
        mask = decisions == code
        summary[name] = {
            "share": round(float(mask.mean()), 4),
            "n": int(mask.sum()),
            "pathogenic_fraction": round(float(labels[mask].mean()), 4) if mask.any() else None,
        }
    for name, label in (("pathogenic", 1), ("benign", 0)):
        of_class = labels == label
        summary[f"{name}_left_uncertain"] = round(float((decisions[of_class] == 1).mean()), 4)
    return summary


def gene_folds(genes: np.ndarray, folds: int = FOLDS, seed: int = SEED) -> np.ndarray:
    unique = sorted(set(genes.tolist()))
    random.Random(seed).shuffle(unique)
    fold_of = {gene: i % folds for i, gene in enumerate(unique)}
    return np.array([fold_of[g] for g in genes.tolist()])


def cross_fit(
    method: str, probabilities: np.ndarray, labels: np.ndarray, folds: np.ndarray
) -> np.ndarray:
    calibrated = np.empty_like(probabilities)
    for fold in np.unique(folds):
        held = folds == fold
        calibrator = fit_calibrator(method, probabilities[~held], labels[~held])
        calibrated[held] = calibrate(calibrator, probabilities[held])
    return calibrated


def reliability_bins(labels: np.ndarray, probabilities: np.ndarray, bins: int = 10) -> list[dict]:
    index = np.minimum((probabilities * bins).astype(int), bins - 1)
    return [
        {
            "bin": b,
            "n": int((index == b).sum()),
            "mean_predicted": round(float(probabilities[index == b].mean()), 4),
            "observed": round(float(labels[index == b].mean()), 4),
        }
        for b in range(bins)
        if (index == b).any()
    ]


def plot_reliability(curves: dict, path: Path, title: str = "validation, missense") -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.plot([0, 1], [0, 1], color="grey", linestyle="--", linewidth=1, label="perfect")
    for name, bins in curves.items():
        x = [b["mean_predicted"] for b in bins]
        y = [b["observed"] for b in bins]
        ax.plot(x, y, marker="o", label=name)
    ax.set_xlabel("Predicted pathogenic probability")
    ax.set_ylabel("Observed pathogenic fraction")
    ax.set_title(f"Reliability ({title})")
    ax.legend(loc="upper left")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=120)
    plt.close(fig)


def load_dual_tower(path: Path, device: str) -> tuple[dual_tower.DualTower, dict]:
    checkpoint = torch.load(path, map_location=device, weights_only=True)
    fields = checkpoint["encoder"]
    encoder = dual_tower.TabularEncoder(
        columns=fields["columns"],
        medians=np.array(fields["medians"]),
        means=np.array(fields["means"]),
        stds=np.array(fields["stds"]),
        flag_columns=fields["flag_columns"],
    )
    config = checkpoint["config"]
    model = dual_tower.DualTower(
        encoder.width, len(dual_tower.ESM_FORMS[config["esm_form"]]), config
    )
    model.load_state_dict(checkpoint["state_dict"])
    return model.to(device), {"encoder": encoder, "config": config}


def tower_probabilities(
    model: dual_tower.DualTower,
    encoder: dual_tower.TabularEncoder,
    esm_form: str,
    table: pl.DataFrame,
    esm_blocks: Path,
    device: str = "cuda",
) -> np.ndarray:
    numeric, amino_acids = encoder.encode(table)
    blocks = variant_inputs.load(esm_blocks, table)[:, dual_tower.ESM_FORMS[esm_form]]
    tensors = dual_tower.Tensors(
        torch.from_numpy(numeric).to(device),
        torch.from_numpy(amino_acids).to(device),
        torch.from_numpy(np.ascontiguousarray(blocks.reshape(len(table), -1))).to(device),
        torch.from_numpy(table["label"].to_numpy().astype(np.float32)).to(device),
    )
    return dual_tower.predict(model.to(device), tensors)


def dual_tower_probabilities(
    checkpoint: Path, esm_blocks: Path, table: pl.DataFrame, device: str = "cuda"
) -> np.ndarray:
    model, meta = load_dual_tower(checkpoint, device)
    return tower_probabilities(
        model, meta["encoder"], meta["config"]["esm_form"], table, esm_blocks, device
    )


def lightgbm_probabilities(model_path: Path, table: pl.DataFrame) -> np.ndarray:
    booster = lgb.Booster(model_file=str(model_path))
    return booster.predict(to_matrix(table, booster.feature_name()))


def run(splits_dir: Path, esm_blocks: Path, models_dir: Path, reports_dir: Path) -> dict:
    started = time.monotonic()
    validation = load_split(splits_dir, "validation")
    labels = validation["label"].to_numpy()
    genes = validation["gene_name"].to_numpy()
    members = {
        "lightgbm": lightgbm_probabilities(models_dir / "lightgbm_full.txt", validation),
        "dual_tower": dual_tower_probabilities(
            models_dir / "dual_tower.pt", esm_blocks, validation
        ),
    }
    ensemble = (members["lightgbm"] + members["dual_tower"]) / 2
    folds = gene_folds(genes)

    methods = {}
    out_of_fold = {}
    for method in ("platt", "isotonic"):
        calibrated = cross_fit(method, ensemble, labels, folds)
        out_of_fold[method] = calibrated
        methods[method] = {
            "metrics": rounded(classification_metrics(labels, calibrated)),
            "decisions": decision_summary(calibrated, labels),
        }
    chosen = min(methods, key=lambda m: methods[m]["metrics"]["brier"])
    calibrator = fit_calibrator(chosen, ensemble, labels)
    calibrated = calibrate(calibrator, ensemble)

    curves = {
        "ensemble (uncalibrated)": reliability_bins(labels, ensemble),
        f"ensemble + {chosen} (out-of-fold)": reliability_bins(labels, out_of_fold[chosen]),
    }
    plot_reliability(curves, reports_dir / "figures" / "calibration_reliability.png")

    model = {
        "members": ["lightgbm_full.txt", "dual_tower.pt"],
        "combination": "mean of member probabilities",
        "calibrator": calibrator,
        "thresholds": THRESHOLDS,
        "decisions": list(DECISIONS),
    }
    (models_dir / "calibration.json").write_text(json.dumps(model, indent=2))

    report = {
        "model": "Missense ensemble (LightGBM + dual tower), calibrated",
        "split": "validation (gene-disjoint)",
        "test_set_opened": False,
        "seed": SEED,
        "folds": FOLDS,
        "rows": int(labels.size),
        "pathogenic_fraction": round(float(labels.mean()), 4),
        "uncalibrated": {
            name: rounded(classification_metrics(labels, p))
            for name, p in (members | {"ensemble": ensemble}).items()
        },
        "methods_out_of_fold": methods,
        "chosen_method": chosen,
        "selection": "lowest out-of-fold Brier score (gene-grouped 5-fold cross-fitting)",
        "final_calibrator_fit_on": "all validation rows",
        "thresholds": THRESHOLDS,
        "threshold_rule": "per-variant calibrated probability (ACMG Bayesian posterior: "
        "likely pathogenic >= 0.90, likely benign <= 0.10)",
        "decisions_in_sample": decision_summary(calibrated, labels),
        "reliability_out_of_fold": curves[f"ensemble + {chosen} (out-of-fold)"],
        "note": "calibrated probabilities assume the validation pathogenic fraction; the final "
        "calibrator is fit on all validation rows, so out-of-fold numbers and the locked test set "
        "(task 6.1) are the honest estimates",
        "seconds": round(time.monotonic() - started, 1),
    }
    (reports_dir / "calibration.json").write_text(json.dumps(report, indent=2))
    return report
