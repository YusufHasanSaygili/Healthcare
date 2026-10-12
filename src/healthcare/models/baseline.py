import json
import random
import time
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl

from healthcare.evaluate.metrics import classification_metrics, gene_bootstrap, rounded
from healthcare.features.aminoacids import ORDER
from healthcare.features.missense import CATEGORICAL_FEATURES, FREQUENCY_FEATURES
from healthcare.split.genes import load_split

SEED = 20261012
EARLY_STOPPING_GENE_FRACTION = 0.1
MAX_ROUNDS = 5000
AMINO_ACID_CODES = {aa: i for i, aa in enumerate(ORDER + "U*X")}
GENE_FEATURE_PREFIXES = ("f_gene_", "f_protein_length")
ABLATIONS = {
    "full": (),
    "no_frequency": tuple(FREQUENCY_FEATURES),
    "no_gene_level": GENE_FEATURE_PREFIXES,
    "no_frequency_no_gene_level": (*FREQUENCY_FEATURES, *GENE_FEATURE_PREFIXES),
}
BASE_PARAMS = {
    "objective": "binary",
    "metric": "binary_logloss",
    "verbosity": -1,
    "seed": SEED,
    "deterministic": True,
    "force_row_wise": True,
    "num_threads": 12,
}


def feature_columns(table: pl.DataFrame, dropped: tuple[str, ...] = ()) -> list[str]:
    return [
        c for c in table.columns if c.startswith("f_") and not any(c.startswith(d) for d in dropped)
    ]


def to_matrix(table: pl.DataFrame, columns: list[str]) -> np.ndarray:
    encoded = table.select(
        [
            pl.col(c).replace_strict(AMINO_ACID_CODES, default=None, return_dtype=pl.Float64)
            if c in CATEGORICAL_FEATURES
            else pl.col(c).cast(pl.Float64)
            for c in columns
        ]
    )
    return encoded.to_numpy()


def early_stopping_genes(train: pl.DataFrame, seed: int = SEED) -> set[str]:
    genes = sorted(train["gene_name"].unique())
    random.Random(seed).shuffle(genes)
    return set(genes[: int(len(genes) * EARLY_STOPPING_GENE_FRACTION)])


def sample_params(rng: random.Random) -> dict:
    return {
        "learning_rate": rng.choice([0.02, 0.03, 0.05, 0.08]),
        "num_leaves": rng.choice([15, 31, 63, 127, 255]),
        "min_child_samples": rng.choice([10, 20, 50, 100, 200]),
        "feature_fraction": rng.choice([0.6, 0.75, 0.9, 1.0]),
        "bagging_fraction": rng.choice([0.7, 0.85, 1.0]),
        "bagging_freq": 1,
        "lambda_l2": rng.choice([0.0, 1.0, 5.0, 20.0]),
        "max_depth": rng.choice([-1, 6, 8, 12]),
    }


class Data:
    def __init__(self, splits_dir: Path, dropped: tuple[str, ...] = ()):
        train = load_split(splits_dir, "train")
        stop_genes = early_stopping_genes(train)
        self.fit = train.filter(~pl.col("gene_name").is_in(stop_genes))
        self.stop = train.filter(pl.col("gene_name").is_in(stop_genes))
        self.validation = load_split(splits_dir, "validation")
        self.columns = feature_columns(train, dropped)
        self.categorical = [i for i, c in enumerate(self.columns) if c in CATEGORICAL_FEATURES]

    def matrix(self, name: str) -> tuple[np.ndarray, np.ndarray]:
        table = getattr(self, name)
        return to_matrix(table, self.columns), table["label"].to_numpy()


def train_model(data: Data, params: dict) -> lgb.Booster:
    x_fit, y_fit = data.matrix("fit")
    x_stop, y_stop = data.matrix("stop")
    fit = lgb.Dataset(x_fit, y_fit, feature_name=data.columns, categorical_feature=data.categorical)
    stop = lgb.Dataset(x_stop, y_stop, reference=fit)
    return lgb.train(
        BASE_PARAMS | params,
        fit,
        num_boost_round=MAX_ROUNDS,
        valid_sets=[stop],
        callbacks=[lgb.early_stopping(100, verbose=False)],
    )


def evaluate(model: lgb.Booster, data: Data, name: str) -> dict:
    x, y = data.matrix(name)
    return classification_metrics(y, model.predict(x, num_iteration=model.best_iteration))


def importance(model: lgb.Booster, data: Data) -> list[dict]:
    x, _ = data.matrix("validation")
    contributions = model.predict(x, num_iteration=model.best_iteration, pred_contrib=True)
    mean_abs = np.abs(contributions[:, :-1]).mean(axis=0)
    order = np.argsort(-mean_abs)
    return [
        {"feature": data.columns[i], "mean_abs_shap": round(float(mean_abs[i]), 4)} for i in order
    ]


def plot_importance(ranking: list[dict], path: Path, top: int = 20) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = ranking[:top][::-1]
    fig, ax = plt.subplots(figsize=(9, 7), dpi=150)
    ax.barh(
        [r["feature"].removeprefix("f_") for r in rows],
        [r["mean_abs_shap"] for r in rows],
        color="#1F7A7A",
    )
    ax.set_xlabel("mean |SHAP| on validation (log-odds)")
    ax.set_title("LightGBM baseline: feature importance")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def sweep(data: Data, trials: int, results_path: Path) -> list[dict]:
    rng = random.Random(SEED)
    results = []
    results_path.parent.mkdir(parents=True, exist_ok=True)
    with results_path.open("w", encoding="utf-8") as handle:
        for trial in range(trials):
            params = sample_params(rng)
            started = time.monotonic()
            model = train_model(data, params)
            row = {
                "trial": trial,
                "params": params,
                "best_iteration": model.best_iteration,
                "validation": rounded(evaluate(model, data, "validation")),
                "seconds": round(time.monotonic() - started, 1),
            }
            handle.write(json.dumps(row) + "\n")
            handle.flush()
            results.append(row)
    return sorted(results, key=lambda r: -r["validation"]["roc_auc"])


def run(splits_dir: Path, artifacts_dir: Path, reports_dir: Path, trials: int) -> dict:
    started = time.monotonic()
    data = Data(splits_dir)
    ranked = sweep(data, trials, artifacts_dir / "lightgbm_sweep.jsonl")
    best = ranked[0]
    model = train_model(data, best["params"])
    model.save_model(str(artifacts_dir / "lightgbm_full.txt"), num_iteration=model.best_iteration)
    ranking = importance(model, data)
    plot_importance(ranking, reports_dir / "figures" / "lightgbm_importance.png")
    x_val, y_val = data.matrix("validation")
    val_probabilities = model.predict(x_val, num_iteration=model.best_iteration)

    ablations = {}
    for name, dropped in ABLATIONS.items():
        ablation_data = data if not dropped else Data(splits_dir, dropped)
        ablation_model = model if not dropped else train_model(ablation_data, best["params"])
        ablations[name] = {
            "features": len(ablation_data.columns),
            "fit": rounded(evaluate(ablation_model, ablation_data, "fit")),
            "validation": rounded(evaluate(ablation_model, ablation_data, "validation")),
        }

    report = {
        "model": "LightGBM baseline (missense, phase 1)",
        "lightgbm_version": lgb.__version__,
        "seed": SEED,
        "rows": {
            "fit": data.fit.height,
            "early_stopping": data.stop.height,
            "validation": data.validation.height,
        },
        "test_set_opened": False,
        "trials": trials,
        "best_params": best["params"],
        "best_iteration": model.best_iteration,
        "validation": rounded(classification_metrics(y_val, val_probabilities)),
        "validation_gene_bootstrap_95ci": gene_bootstrap(
            y_val, val_probabilities, data.validation["gene_name"].to_numpy()
        ),
        "ablations": ablations,
        "top_features": ranking[:12],
        "sweep_top5": [
            {
                "trial": r["trial"],
                "roc_auc": r["validation"]["roc_auc"],
                "pr_auc": r["validation"]["pr_auc"],
            }
            for r in ranked[:5]
        ],
        "seconds": round(time.monotonic() - started, 1),
    }
    (reports_dir / "baseline_lightgbm.json").write_text(json.dumps(report, indent=2))
    return report
