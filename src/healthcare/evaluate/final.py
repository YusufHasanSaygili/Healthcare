import json
import time
from pathlib import Path

import lightgbm as lgb
import numpy as np
import onnx
import polars as pl

from healthcare.data.sources import load_sources
from healthcare.esm import variant_inputs
from healthcare.evaluate.metrics import classification_metrics, gene_bootstrap, rounded
from healthcare.features.missense import FREQUENCY_FEATURES, to_matrix
from healthcare.models import baseline, calibration, dual_tower, export
from healthcare.models.baseline import ESM_PREFIX
from healthcare.models.sweep import CONFIRM_SEEDS
from healthcare.split.genes import file_hash, load_split

REPORT_NAME = "missense_v1"
MEMBERS = ("ensemble_calibrated", "lightgbm", "dual_tower")
CURVE_SETS = ("test", "temporal_test_gene")


def ablation_specs(esm_form: str) -> dict:
    return {
        "no_frequency": {
            "lightgbm": tuple(FREQUENCY_FEATURES),
            "dual_tower": (dual_tower.NO_FREQUENCY, esm_form),
        },
        "no_esm": {
            "lightgbm": ESM_PREFIX,
            "dual_tower": (dual_tower.ESM_SCALARS, "none"),
        },
    }


def train_ablations(splits_dir: Path, esm_blocks: Path, reports_dir: Path) -> dict:
    lightgbm_params = json.loads((reports_dir / "baseline_lightgbm.json").read_text())[
        "best_params"
    ]
    tower_params = json.loads((reports_dir / "dual_tower.json").read_text())["winner"]["params"]
    tower_data = dual_tower.Data(splits_dir, esm_blocks)
    models = {}
    for name, spec in ablation_specs(tower_params["esm_form"]).items():
        booster = baseline.train_model(baseline.Data(splits_dir, spec["lightgbm"]), lightgbm_params)
        dropped, form = spec["dual_tower"]
        trained = dual_tower.train(
            tower_data, tower_params | {"esm_form": form}, CONFIRM_SEEDS[0], dropped
        )
        models[name] = {
            "lightgbm": booster,
            "dual_tower": (trained["model"], trained["encoder"], form),
        }
    del tower_data
    return models


def locked_blocks(table: pl.DataFrame, name: str, esm_store: Path, blocks_dir: Path) -> tuple:
    path = blocks_dir / f"locked_{name}.npy"
    if not path.exists():
        variant_inputs.build(table, esm_store, path)
    blocks = variant_inputs.load(path, table)
    missing = int((~blocks.reshape(table.height, -1).any(axis=1)).sum())
    return path, missing


def previous_onnx_hash(report_path: Path) -> str | None:
    if not report_path.exists():
        return None
    return json.loads(report_path.read_text()).get("onnx_sha256")


def onnx_outputs(model: onnx.ModelProto, table: pl.DataFrame, blocks_path: Path) -> dict:
    metadata = {p.key: json.loads(p.value) for p in model.metadata_props}
    indices = [variant_inputs.BLOCKS.index(b) for b in metadata["esm_blocks"]]
    protein = variant_inputs.load(blocks_path, table)[:, indices].reshape(table.height, -1)
    inputs = {
        "features": to_matrix(table, metadata["feature_names"]).astype(np.float32),
        "protein": np.ascontiguousarray(protein, dtype=np.float32),
    }
    outputs = export.predict(export.session(model), inputs)
    return {
        "ensemble_calibrated": outputs["probability"].astype(np.float64),
        "lightgbm": outputs["lightgbm"].astype(np.float64),
        "dual_tower": outputs["dual_tower"].astype(np.float64),
    }


def ablation_outputs(ablations: dict, table: pl.DataFrame, blocks_path: Path) -> dict:
    outputs = {}
    for name, models in ablations.items():
        booster = models["lightgbm"]
        boosted = booster.predict(
            to_matrix(table, booster.feature_name()), num_iteration=booster.best_iteration
        )
        model, encoder, form = models["dual_tower"]
        tower = calibration.tower_probabilities(model, encoder, form, table, blocks_path)
        outputs[f"{name}/lightgbm"] = boosted
        outputs[f"{name}/dual_tower"] = tower.astype(np.float64)
        outputs[f"{name}/ensemble_uncalibrated"] = (boosted + tower) / 2
    return outputs


def evaluation_sets(test: pl.DataFrame, temporal: pl.DataFrame) -> dict:
    sets = {
        "test": ("test", np.ones(test.height, dtype=bool)),
        "temporal_all": ("temporal", np.ones(temporal.height, dtype=bool)),
    }
    subsets = temporal["temporal_subset"]
    for subset in sorted(subsets.unique()):
        sets[f"temporal_{subset}"] = ("temporal", (subsets == subset).to_numpy())
    return sets


def score(labels: np.ndarray, probabilities: np.ndarray, genes: np.ndarray) -> dict:
    if labels.min() == labels.max():
        return {"n": int(labels.size), "positives": int(labels.sum()), "note": "single class"}
    return rounded(classification_metrics(labels, probabilities)) | {
        "gene_bootstrap_95ci": gene_bootstrap(labels, probabilities, genes)
    }


def plot_curves(labels: np.ndarray, outputs: dict, title: str, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import precision_recall_curve, roc_curve

    fig, (roc_ax, pr_ax) = plt.subplots(1, 2, figsize=(11, 5))
    for name in MEMBERS:
        fpr, tpr, _ = roc_curve(labels, outputs[name])
        precision, recall, _ = precision_recall_curve(labels, outputs[name])
        roc_ax.plot(fpr, tpr, label=name)
        pr_ax.plot(recall, precision, label=name)
    roc_ax.plot([0, 1], [0, 1], color="grey", linestyle="--", linewidth=1)
    roc_ax.set(xlabel="False positive rate", ylabel="True positive rate", title=f"ROC ({title})")
    pr_ax.axhline(labels.mean(), color="grey", linestyle="--", linewidth=1)
    pr_ax.set(xlabel="Recall", ylabel="Precision", title=f"Precision-recall ({title})")
    roc_ax.legend(loc="lower right")
    pr_ax.legend(loc="lower left")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=120)
    plt.close(fig)


def run(
    splits_dir: Path, esm_store: Path, esm_blocks: Path, models_dir: Path, reports_dir: Path
) -> dict:
    started = time.monotonic()
    ablations = train_ablations(splits_dir, esm_blocks, reports_dir)
    onnx_path = models_dir / export.ONNX_NAME
    model = onnx.load(str(onnx_path))
    onnx_hash = file_hash(onnx_path)
    previous_hash = previous_onnx_hash(reports_dir / f"{REPORT_NAME}.json")

    test = load_split(splits_dir, "test", allow_locked=True)
    temporal = load_split(splits_dir, "temporal", allow_locked=True)
    blocks = {}
    missing_blocks = {}
    for name, table in (("test", test), ("temporal", temporal)):
        blocks[name], missing_blocks[name] = locked_blocks(
            table, name, esm_store, esm_blocks.parent
        )

    predictions = {}
    for name, table in (("test", test), ("temporal", temporal)):
        predictions[name] = onnx_outputs(model, table, blocks[name]) | ablation_outputs(
            ablations, table, blocks[name]
        )

    results = {}
    figures = reports_dir / "figures"
    reliability = {}
    tables = {"test": test, "temporal": temporal}
    for set_name, (source, mask) in evaluation_sets(test, temporal).items():
        table = tables[source].filter(pl.Series(mask))
        labels = table["label"].to_numpy()
        genes = table["gene_name"].to_numpy()
        outputs = {name: values[mask] for name, values in predictions[source].items()}
        results[set_name] = {
            "variants": int(labels.size),
            "genes": int(np.unique(genes).size),
            "pathogenic": int(labels.sum()),
            "shares_position_with_train": int(
                table.get_column("shares_position_with_train", default=pl.Series([False])).sum()
            ),
            "models": {name: score(labels, values, genes) for name, values in outputs.items()},
        }
        if labels.min() != labels.max():
            results[set_name]["decisions"] = calibration.decision_summary(
                outputs["ensemble_calibrated"], labels
            )
            reliability[set_name] = calibration.reliability_bins(
                labels, outputs["ensemble_calibrated"]
            )
        if set_name in CURVE_SETS:
            plot_curves(labels, outputs, set_name, figures / f"{REPORT_NAME}_curves_{set_name}.png")
    calibration.plot_reliability(
        {name: reliability[name] for name in ("test", "temporal_all") if name in reliability},
        figures / f"{REPORT_NAME}_reliability.png",
        title="final model, missense",
    )

    report = {
        "model": "Missense ensemble (LightGBM + dual tower, Platt), ONNX",
        "variant_type": "missense",
        "test_set_opened": True,
        "model_changed_after_test": previous_hash is not None and previous_hash != onnx_hash,
        "onnx": export.ONNX_NAME,
        "onnx_sha256": onnx_hash,
        "previous_onnx_sha256": previous_hash,
        "seeds": {"split": 20261012, "lightgbm": baseline.SEED, "dual_tower": CONFIRM_SEEDS[0]},
        "data_versions": {name: source.version for name, source in load_sources().items()},
        "library_versions": {
            "lightgbm": lgb.__version__,
            "onnx": onnx.__version__,
            "onnxruntime": export.ort.__version__,
            "torch": export.torch.__version__,
        },
        "ablations": {
            "note": "retrained with one seed from the frozen hyperparameters on training data "
            "only, before the locked sets were read; the full model is the exported ONNX (float32 "
            "trees); ensembles of ablated members are uncalibrated, so compare their ROC/PR, not "
            "Brier/ECE",
            "dropped": {
                name: {
                    "lightgbm": list(spec["lightgbm"]),
                    "dual_tower_esm_form": spec["dual_tower"][1],
                }
                for name, spec in ablation_specs(ablations["no_esm"]["dual_tower"][2]).items()
            },
        },
        "missing_esm_blocks": missing_blocks,
        "temporal_note": "temporal_train_gene and temporal_validation_gene share genes with data "
        "the model saw (see shares_position_with_train per set); temporal_test_gene is the "
        "gene-disjoint temporal estimate",
        "results": results,
        "reliability": reliability,
        "seconds": round(time.monotonic() - started, 1),
    }
    (reports_dir / f"{REPORT_NAME}.json").write_text(json.dumps(report, indent=2))
    (reports_dir / f"{REPORT_NAME}.md").write_text(markdown(report), encoding="utf-8")
    return report


def metric_row(name: str, metrics: dict) -> str:
    if "roc_auc" not in metrics:
        return f"| {name} | – | – | – | – | – | – |"
    low, high = metrics["gene_bootstrap_95ci"]["roc_auc"]
    return (
        f"| {name} | {metrics['roc_auc']:.4f} [{low:.4f}–{high:.4f}] | {metrics['pr_auc']:.4f} "
        f"| {metrics['mcc']:.4f} | {metrics['balanced_accuracy']:.4f} | {metrics['brier']:.4f} "
        f"| {metrics['ece']:.4f} |"
    )


def markdown(report: dict) -> str:
    lines = [
        "# Missense model v1 — final evaluation",
        "",
        "Variant type: **missense SNVs only** (GRCh38, MANE Select). "
        "For research and educational use only; not for clinical diagnosis or treatment decisions.",
        "",
        f"- Model: {report['model']} (`{report['onnx']}`)",
        "- Test set opened in this report for the first time; model changed after testing: "
        f"**{'yes' if report['model_changed_after_test'] else 'no'}**",
        "- Thresholds and calibration were fit on the validation set only "
        "(likely pathogenic ≥ 0.90, likely benign ≤ 0.10).",
        f"- Ablations: {report['ablations']['note']}.",
        f"- Temporal sets: {report['temporal_note']}.",
        "",
    ]
    for set_name, result in report["results"].items():
        lines += [
            f"## {set_name}",
            "",
            f"{result['variants']} variants, {result['genes']} genes, "
            f"{result['pathogenic']} pathogenic.",
            "",
            "| Model | ROC-AUC [95% gene CI] | PR-AUC | MCC | Bal. acc. | Brier | ECE |",
            "|---|---|---|---|---|---|---|",
        ]
        lines += [metric_row(name, m) for name, m in result["models"].items()]
        if "decisions" in result:
            decisions = result["decisions"]
            lines += [
                "",
                "| Decision | Share | Pathogenic fraction |",
                "|---|---|---|",
            ]
            for name in calibration.DECISIONS:
                d = decisions[name]
                fraction = (
                    "–" if d["pathogenic_fraction"] is None else f"{d['pathogenic_fraction']:.3f}"
                )
                lines.append(f"| {name} | {d['share']:.3f} | {fraction} |")
        lines.append("")
    lines += [
        "## Figures",
        "",
        f"- `figures/{REPORT_NAME}_curves_test.png`, "
        f"`figures/{REPORT_NAME}_curves_temporal_test_gene.png`, "
        f"`figures/{REPORT_NAME}_reliability.png`",
        "",
    ]
    return "\n".join(lines)
