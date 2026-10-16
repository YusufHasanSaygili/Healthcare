import json
import random
import time
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl

from healthcare.esm.model import MODEL_ID as ESM_MODEL_ID
from healthcare.evaluate.metrics import classification_metrics, gene_bootstrap, rounded
from healthcare.features.missense import CATEGORICAL_FEATURES, FREQUENCY_FEATURES
from healthcare.models.baseline import (
    BASE_PARAMS,
    ESM_PREFIX,
    SEED,
    Data,
    sample_params,
    train_model,
)
from healthcare.split.genes import KEY, load_split

MIN_CLASS = 20
CATEGORICAL = (*CATEGORICAL_FEATURES, "f_variant_type")
ABLATIONS = {
    "full": (),
    "no_frequency": tuple(FREQUENCY_FEATURES),
    "no_esm": ESM_PREFIX,
}


def load_joint(splits_dir: Path, name: str, allow_locked: bool = False) -> pl.DataFrame:
    multi = pl.read_parquet(splits_dir / "multi_gene_keys.parquet")
    missense = load_split(splits_dir, name, allow_locked, "missense").with_columns(
        pl.lit("missense").alias("consequence"), pl.lit(0.0).alias("f_variant_type")
    )
    phase2 = load_split(splits_dir, name, allow_locked, "phase2")
    joined = pl.concat([missense, phase2], how="diagonal_relaxed")
    return joined.join(multi, on=KEY, how="anti")


def joint_data(splits_dir: Path, dropped: tuple[str, ...] = ()) -> Data:
    return Data(splits_dir, dropped, loader=load_joint, categorical=CATEGORICAL)


def eligible(labels: np.ndarray) -> bool:
    positives = int(labels.sum())
    return positives >= MIN_CLASS and labels.size - positives >= MIN_CLASS


def per_type(table: pl.DataFrame, probabilities: np.ndarray, bootstrap: bool = False) -> dict:
    report = {}
    consequences = table["consequence"].to_numpy()
    labels = table["label"].to_numpy()
    genes = table["gene_name"].to_numpy()
    for consequence in sorted(set(consequences)):
        mask = consequences == consequence
        y, p = labels[mask], probabilities[mask]
        if not eligible(y):
            report[consequence] = {"n": int(y.size), "positives": int(y.sum()), "eligible": False}
            continue
        entry = rounded(classification_metrics(y, p)) | {"eligible": True}
        if bootstrap:
            entry["gene_bootstrap_95ci"] = gene_bootstrap(y, p, genes[mask])
        report[consequence] = entry
    return report


def selection_score(report: dict) -> float:
    return float(np.mean([m["roc_auc"] for m in report.values() if m["eligible"]]))


def predict(model: lgb.Booster, data: Data, name: str) -> np.ndarray:
    x, _ = data.matrix(name)
    return model.predict(x, num_iteration=model.best_iteration)


def sweep(data: Data, trials: int, results_path: Path) -> list[dict]:
    rng = random.Random(SEED)
    results = []
    results_path.parent.mkdir(parents=True, exist_ok=True)
    with results_path.open("w", encoding="utf-8") as handle:
        for trial in range(trials):
            params = sample_params(rng)
            started = time.monotonic()
            model = train_model(data, params)
            report = per_type(data.validation, predict(model, data, "validation"))
            row = {
                "trial": trial,
                "params": params,
                "best_iteration": model.best_iteration,
                "mean_type_roc_auc": round(selection_score(report), 4),
                "missense_roc_auc": report["missense"]["roc_auc"],
                "seconds": round(time.monotonic() - started, 1),
            }
            handle.write(json.dumps(row) + "\n")
            handle.flush()
            results.append(row)
    return sorted(results, key=lambda r: -r["mean_type_roc_auc"])


def run(splits_dir: Path, artifacts_dir: Path, reports_dir: Path, trials: int) -> dict:
    started = time.monotonic()
    data = joint_data(splits_dir)
    ranked = sweep(data, trials, artifacts_dir / "joint_lightgbm_sweep.jsonl")
    best = ranked[0]
    model = train_model(data, best["params"])
    model.save_model(str(artifacts_dir / "joint_lightgbm.txt"), num_iteration=model.best_iteration)
    ablations = {}
    for name, dropped in ABLATIONS.items():
        ablation_data = data if not dropped else joint_data(splits_dir, dropped)
        ablation_model = model if not dropped else train_model(ablation_data, best["params"])
        report = per_type(
            ablation_data.validation, predict(ablation_model, ablation_data, "validation")
        )
        ablations[name] = {
            "features": len(ablation_data.columns),
            "mean_type_roc_auc": round(selection_score(report), 4),
            "roc_auc_by_type": {t: m.get("roc_auc") for t, m in report.items()},
        }
    validation = per_type(data.validation, predict(model, data, "validation"), bootstrap=True)
    phase1 = json.loads((reports_dir / "baseline_lightgbm.json").read_text())
    report = {
        "model": "LightGBM baseline (all coding types, phase 2)",
        "lightgbm_version": lgb.__version__,
        "esm_model": ESM_MODEL_ID,
        "seed": SEED,
        "base_params": BASE_PARAMS,
        "rows": {
            "fit": data.fit.height,
            "early_stopping": data.stop.height,
            "validation": data.validation.height,
        },
        "test_set_opened": False,
        "selection": f"mean validation ROC-AUC over types with at least {MIN_CLASS} per class",
        "trials": trials,
        "best_params": best["params"],
        "best_iteration": model.best_iteration,
        "validation_by_type": validation,
        "missense_vs_phase1": {
            "phase1_lightgbm": phase1["ablations"]["full"]["validation"]["roc_auc"],
            "joint_lightgbm": validation["missense"]["roc_auc"],
        },
        "ablations": ablations,
        "sweep_top5": [
            {k: r[k] for k in ("trial", "mean_type_roc_auc", "missense_roc_auc")}
            for r in ranked[:5]
        ],
        "seconds": round(time.monotonic() - started, 1),
    }
    (reports_dir / "joint_lightgbm.json").write_text(json.dumps(report, indent=2))
    return report
