import hashlib
import json
import time
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
import torch
from sklearn.linear_model import LogisticRegression

from healthcare.esm import variant_inputs
from healthcare.features.missense import to_matrix
from healthcare.models import dual_tower
from healthcare.models.baseline import train_model
from healthcare.models.calibration import THRESHOLDS, decide
from healthcare.models.joint import (
    ABLATIONS,
    MIN_CLASS,
    joint_data,
    load_joint,
    per_type,
)
from healthcare.models.joint_tower import ESM_COLUMNS

POOLED = "pooled_small_types"
DECISIONS = ("likely_pathogenic", "uncertain", "likely_benign")
EPSILON = 1e-6


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def logit(probabilities: np.ndarray) -> np.ndarray:
    clipped = np.clip(probabilities, EPSILON, 1 - EPSILON)
    return np.log(clipped / (1 - clipped))


def lightgbm_predict(model_path: Path, table: pl.DataFrame) -> np.ndarray:
    model = lgb.Booster(model_file=str(model_path))
    return model.predict(to_matrix(table, model.feature_name()))


def tower_predict(checkpoint_path: Path, table: pl.DataFrame, blocks_path: Path) -> np.ndarray:
    checkpoint = torch.load(checkpoint_path, weights_only=False)
    encoder, config = checkpoint["encoder"], checkpoint["config"]
    model = dual_tower.DualTower(
        encoder.width, len(dual_tower.ESM_FORMS[config["esm_form"]]), config
    )
    model.load_state_dict(checkpoint["state"])
    numeric, amino_acids = encoder.encode(table)
    blocks = variant_inputs.load(blocks_path, table, ESM_COLUMNS)
    device = "cuda"
    tensors = dual_tower.Tensors(
        torch.from_numpy(numeric).to(device),
        torch.from_numpy(amino_acids).to(device),
        torch.from_numpy(np.ascontiguousarray(blocks.reshape(len(table), -1))).to(device),
        torch.from_numpy(table["label"].to_numpy().astype(np.float32)).to(device),
        torch.from_numpy(table["protein_position"].is_not_null().to_numpy().astype(np.float32)).to(
            device
        ),
    )
    return dual_tower.predict(model.to(device), tensors)


def calibration_group(consequence: str, groups: dict[str, str]) -> str:
    return groups.get(consequence, POOLED)


def blend(table: pl.DataFrame, lightgbm: np.ndarray, tower: np.ndarray, routes: dict) -> np.ndarray:
    averaged = table["consequence"].replace_strict(routes, default="lightgbm") == "average"
    return np.where(averaged.to_numpy(), (lightgbm + tower) / 2, lightgbm)


def calibrate(table: pl.DataFrame, raw: np.ndarray, platt: dict, groups: dict) -> np.ndarray:
    keys = [calibration_group(c, groups) for c in table["consequence"].to_list()]
    slopes = np.array([platt[k][0] for k in keys])
    intercepts = np.array([platt[k][1] for k in keys])
    return 1 / (1 + np.exp(-(slopes * logit(raw) + intercepts)))


def fit_platt(raw: np.ndarray, labels: np.ndarray) -> list[float]:
    model = LogisticRegression(C=1e6, max_iter=1000).fit(logit(raw)[:, None], labels)
    return [float(model.coef_[0, 0]), float(model.intercept_[0])]


def choose(splits_dir: Path, blocks_path: Path, models_dir: Path, reports_dir: Path) -> dict:
    started = time.monotonic()
    validation = load_joint(splits_dir, "validation")
    lightgbm = lightgbm_predict(models_dir / "joint_lightgbm.txt", validation)
    tower = tower_predict(models_dir / "joint_dual_tower.pt", validation, blocks_path)
    averaged = (lightgbm + tower) / 2
    by_lightgbm, by_average = per_type(validation, lightgbm), per_type(validation, averaged)
    routes = {
        consequence: "average"
        if metrics["eligible"] and by_average[consequence]["roc_auc"] > metrics["roc_auc"]
        else "lightgbm"
        for consequence, metrics in by_lightgbm.items()
    }
    raw = blend(validation, lightgbm, tower, routes)
    labels = validation["label"].to_numpy()
    consequences = validation["consequence"].to_numpy()
    groups = {c: c for c, m in by_lightgbm.items() if m["eligible"]}
    platt = {}
    for group in [*groups, POOLED]:
        mask = np.array([calibration_group(c, groups) == group for c in consequences])
        if mask.any():
            platt[group] = fit_platt(raw[mask], labels[mask])
    calibrated = calibrate(validation, raw, platt, groups)
    ablation_models = {}
    params = json.loads((reports_dir / "joint_lightgbm.json").read_text())["best_params"]
    for name, dropped in ABLATIONS.items():
        if not dropped:
            continue
        data = joint_data(splits_dir, dropped)
        path = models_dir / f"joint_lightgbm_{name}.txt"
        model = train_model(data, params)
        model.save_model(str(path), num_iteration=model.best_iteration)
        ablation_models[name] = path.name
    final = {
        "routes": routes,
        "calibration_groups": groups,
        "platt": platt,
        "thresholds": THRESHOLDS,
        "members": {
            "lightgbm": "joint_lightgbm.txt",
            "dual_tower": "joint_dual_tower.pt",
        },
        "ablation_models": ablation_models,
        "chosen_on": "validation only",
        "min_class_per_type": MIN_CLASS,
    }
    final_path = models_dir / "v2_final.json"
    final_path.write_text(json.dumps(final, indent=2))
    return {
        "final": final_path,
        "routes": routes,
        "validation_calibrated": per_type(validation, calibrated),
        "seconds": round(time.monotonic() - started, 1),
    }


def decision_shares(probabilities: np.ndarray) -> dict:
    codes = decide(probabilities)
    total = max(len(codes), 1)
    return {
        name: round(float((codes == 2 - i).sum()) / total, 4) for i, name in enumerate(DECISIONS)
    }


def type_report(table: pl.DataFrame, probabilities: np.ndarray) -> dict:
    report = per_type(table, probabilities, bootstrap=True)
    consequences = table["consequence"].to_numpy()
    for consequence, entry in report.items():
        entry["decisions"] = decision_shares(probabilities[consequences == consequence])
    return report


def predictions(table: pl.DataFrame, models_dir: Path, blocks_path: Path, final: dict) -> dict:
    lightgbm = lightgbm_predict(models_dir / final["members"]["lightgbm"], table)
    tower = tower_predict(models_dir / final["members"]["dual_tower"], table, blocks_path)
    raw = blend(table, lightgbm, tower, final["routes"])
    calibrated = calibrate(table, raw, final["platt"], final["calibration_groups"])
    ablations = {
        name: lightgbm_predict(models_dir / path, table)
        for name, path in final["ablation_models"].items()
    }
    return {"final": calibrated, "lightgbm": lightgbm, "dual_tower": tower} | ablations


def roc_by_type(table: pl.DataFrame, probabilities: np.ndarray) -> dict:
    return {t: m.get("roc_auc") for t, m in per_type(table, probabilities).items()}


def evaluate_set(table: pl.DataFrame, scores: dict) -> dict:
    return {
        "variants": table.height,
        "final": type_report(table, scores["final"]),
        "compare_roc_auc": {
            name: roc_by_type(table, values) for name, values in scores.items() if name != "final"
        },
    }


class AlreadyOpenedError(RuntimeError):
    pass


def evaluate(
    splits_dir: Path,
    store_dir: Path,
    locked_blocks: Path,
    models_dir: Path,
    reports_dir: Path,
    force: bool = False,
) -> dict:
    started = time.monotonic()
    report_path = reports_dir / "all_types_v2.json"
    if report_path.exists() and not force:
        raise AlreadyOpenedError(f"{report_path} exists: the test set was already opened")
    final = json.loads((models_dir / "v2_final.json").read_text())
    pinned = ["v2_final.json", *final["members"].values(), *final["ablation_models"].values()]
    hashes_before = {name: file_hash(models_dir / name) for name in pinned}
    test = load_joint(splits_dir, "test", allow_locked=True)
    temporal = load_joint(splits_dir, "temporal", allow_locked=True)
    locked = pl.concat([test, temporal], how="diagonal_relaxed").unique(
        ["chrom", "pos", "ref", "alt"], maintain_order=True
    )
    blocks = variant_inputs.build(locked, store_dir, locked_blocks, ESM_COLUMNS)
    sets = {
        "test": test,
        "temporal_all": temporal,
        "temporal_test_gene": temporal.filter(pl.col("temporal_subset") == "test_gene"),
        "temporal_validation_gene": temporal.filter(pl.col("temporal_subset") == "validation_gene"),
    }
    results = {
        name: evaluate_set(table, predictions(table, models_dir, locked_blocks, final))
        for name, table in sets.items()
    }
    phase1 = json.loads((reports_dir / "missense_v1.json").read_text())
    hashes_after = {name: file_hash(models_dir / name) for name in hashes_before}
    report = {
        "model": "v2: joint LightGBM, averaged with the joint dual tower per type, per-type Platt",
        "routes": final["routes"],
        "thresholds": final["thresholds"],
        "test_set_opened": True,
        "rerun": report_path.exists(),
        "model_hashes": hashes_before,
        "model_changed_after_test": hashes_before != hashes_after,
        "locked_esm_blocks": blocks,
        "temporal_subsets": dict(
            temporal.group_by("temporal_subset").len().sort("temporal_subset").iter_rows()
        ),
        "temporal_sharing_a_train_position": int(temporal["shares_position_with_train"].sum()),
        "results": results,
        "missense_v1_test_roc_auc": phase1["results"]["test"]["models"]["ensemble_calibrated"][
            "roc_auc"
        ],
        "seconds": round(time.monotonic() - started, 1),
    }
    report_path.write_text(json.dumps(report, indent=2, default=str))
    (reports_dir / "all_types_v2.md").write_text(markdown(report), encoding="utf-8")
    return report


def number(value) -> str:
    return "–" if value is None else f"{value:.4f}"


TITLES = {
    "test": "Gene-disjoint test",
    "temporal_all": "Temporal (added to ClinVar after October 2025), all genes",
    "temporal_test_gene": "Temporal, test genes only (never seen in training or model choice)",
    "temporal_validation_gene": "Temporal, validation genes (used for model choice)",
}


def markdown(report: dict) -> str:
    subsets = ", ".join(f"{k}: {v:,}" for k, v in report["temporal_subsets"].items())
    lines = [
        "# All coding types v2 — final evaluation",
        "",
        "For research and educational use only; not for clinical diagnosis or treatment decisions.",
        "",
        f"Test opened once. Model changed after opening: {report['model_changed_after_test']}.",
        "",
        "Routes chosen on validation: "
        + ", ".join(f"{t} → {r}" for t, r in sorted(report["routes"].items())),
        "",
        f"Temporal variants by gene set: {subsets}; "
        f"{report['temporal_sharing_a_train_position']:,} share a training position.",
        "",
        "Every table is per variant type; no number mixes types. Comparison columns are ROC-AUC of "
        "the joint LightGBM, the joint dual tower, and LightGBM trained without allele frequency "
        "or without ESM features.",
        "",
    ]
    for name, title in TITLES.items():
        result = report["results"][name]
        lines += [
            f"## {title} ({result['variants']:,} variants)",
            "",
            "| Type | n | Pathogenic | ROC-AUC [95% CI] | PR-AUC | MCC | Bal. acc. | Brier | ECE "
            "| Uncertain | LightGBM | Dual tower | LightGBM no-AF | LightGBM no-ESM |",
            "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
        ]
        compare = result["compare_roc_auc"]
        for consequence, entry in result["final"].items():
            if not entry["eligible"]:
                lines.append(
                    f"| {consequence} | {entry['n']} | {entry['positives']} | too few of one class "
                    "| | | | | | | | | | |"
                )
                continue
            ci = entry["gene_bootstrap_95ci"]["roc_auc"]
            lines.append(
                f"| {consequence} | {entry['n']} | {entry['positives']} | "
                f"{entry['roc_auc']:.4f} [{ci[0]:.4f}–{ci[1]:.4f}] | {entry['pr_auc']:.4f} | "
                f"{entry['mcc']:.4f} | {entry['balanced_accuracy']:.4f} | {entry['brier']:.4f} | "
                f"{entry['ece']:.4f} | {entry['decisions']['uncertain']:.1%} | "
                f"{number(compare['lightgbm'].get(consequence))} | "
                f"{number(compare['dual_tower'].get(consequence))} | "
                f"{number(compare['no_frequency'].get(consequence))} | "
                f"{number(compare['no_esm'].get(consequence))} |"
            )
        lines.append("")
    lines += [
        "Phase 1 missense model on the same test genes: ROC-AUC "
        f"{report['missense_v1_test_roc_auc']:.4f}.",
        "",
        "Allele frequency: for loss-of-function types most ClinVar benign labels follow population "
        "frequency, so the no-AF column is the honest measure of what the model learned beyond it.",
        "",
    ]
    return "\n".join(lines)
