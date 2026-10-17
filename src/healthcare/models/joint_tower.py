import json
import time
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
import torch

from healthcare.esm import variant_inputs
from healthcare.esm.model import MODEL_ID as ESM_MODEL_ID
from healthcare.features.esm_features import ESM_FEATURES, SPAN_FEATURES
from healthcare.features.missense import FREQUENCY_FEATURES, to_matrix
from healthcare.models import dual_tower
from healthcare.models.joint import load_joint, per_type, selection_score

SEEDS = (0, 1, 2)
ESM_COLUMNS = [0, 1]
NO_ESM = (*ESM_FEATURES, *SPAN_FEATURES)
ABLATIONS = {
    "no_frequency": (tuple(FREQUENCY_FEATURES), None),
    "no_esm": (NO_ESM, "none"),
}


def prepare(splits_dir: Path, store_dir: Path, output: Path) -> dict:
    table = pl.concat(
        [load_joint(splits_dir, "train"), load_joint(splits_dir, "validation")],
        how="diagonal_relaxed",
    )
    return variant_inputs.build(table, store_dir, output, ESM_COLUMNS)


def tower_validation(data: dual_tower.Data, result: dict, dropped: tuple, esm_form: str):
    _, tensors = data.tensors(dropped, esm_form)
    return dual_tower.predict(result["model"], tensors["validation"])


def lightgbm_validation(model_path: Path, validation: pl.DataFrame) -> np.ndarray:
    model = lgb.Booster(model_file=str(model_path))
    return model.predict(to_matrix(validation, model.feature_name()))


def summarize(reports: list[dict]) -> dict:
    types = [t for t, m in reports[0].items() if m["eligible"]]
    return {
        t: {
            "roc_auc_mean": round(float(np.mean([r[t]["roc_auc"] for r in reports])), 4),
            "roc_auc_std": round(float(np.std([r[t]["roc_auc"] for r in reports])), 4),
            "pr_auc_mean": round(float(np.mean([r[t]["pr_auc"] for r in reports])), 4),
        }
        for t in types
    }


def run(
    splits_dir: Path,
    blocks_path: Path,
    models_dir: Path,
    reports_dir: Path,
    device: str = "cuda",
) -> dict:
    started = time.monotonic()
    config = json.loads((reports_dir / "dual_tower.json").read_text())["winner"]["params"]
    data = dual_tower.Data(splits_dir, blocks_path, device, load_joint, ESM_COLUMNS)
    validation = data.tables["validation"]
    seed_reports, seed_probabilities, runs = [], [], []
    for seed in SEEDS:
        result = dual_tower.train(data, config, seed)
        probabilities = tower_validation(data, result, (), config["esm_form"])
        seed_reports.append(per_type(validation, probabilities))
        seed_probabilities.append(probabilities)
        runs.append(
            {"seed": seed, "best_epoch": result["best_epoch"], "seconds": result["seconds"]}
        )
        if seed == SEEDS[0]:
            torch.save(
                {
                    "state": result["model"].state_dict(),
                    "config": config,
                    "encoder": result["encoder"],
                },
                models_dir / "joint_dual_tower.pt",
            )
    tower = seed_probabilities[0]
    lightgbm = lightgbm_validation(models_dir / "joint_lightgbm.txt", validation)
    ensemble = (tower + lightgbm) / 2
    ablations = {}
    for name, (dropped, form) in ABLATIONS.items():
        esm_form = form or config["esm_form"]
        result = dual_tower.train(data, config | {"esm_form": esm_form}, SEEDS[0], dropped)
        report = per_type(validation, tower_validation(data, result, dropped, esm_form))
        ablations[name] = {
            "mean_type_roc_auc": round(selection_score(report), 4),
            "roc_auc_by_type": {t: m.get("roc_auc") for t, m in report.items()},
        }
    joint_lgbm = json.loads((reports_dir / "joint_lightgbm.json").read_text())
    report = {
        "model": "Dual tower (all coding types, phase 2)",
        "config": config,
        "esm_model": ESM_MODEL_ID,
        "torch_version": torch.__version__,
        "seeds": list(SEEDS),
        "test_set_opened": False,
        "rows": {name: table.height for name, table in data.tables.items()},
        "protein_tower_masked_rows": {
            name: int(table["protein_position"].is_null().sum())
            for name, table in data.tables.items()
        },
        "validation_by_type_seeds": summarize(seed_reports),
        "mean_type_roc_auc_seed0": round(selection_score(seed_reports[0]), 4),
        "joint_lightgbm_by_type": {
            t: m.get("roc_auc") for t, m in joint_lgbm["validation_by_type"].items()
        },
        "ensemble_by_type": per_type(validation, ensemble, bootstrap=True),
        "ensemble_mean_type_roc_auc": round(selection_score(per_type(validation, ensemble)), 4),
        "ablations": ablations,
        "runs": runs,
        "seconds": round(time.monotonic() - started, 1),
    }
    (reports_dir / "joint_dual_tower.json").write_text(json.dumps(report, indent=2))
    return report
