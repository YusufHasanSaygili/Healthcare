import io
import json
import time
import warnings
from pathlib import Path

import lightgbm as lgb
import numpy as np
import onnx
import onnxruntime as ort
import torch
from torch import nn

from healthcare.esm import variant_inputs
from healthcare.evaluate import all_types
from healthcare.features.missense import CATEGORICAL_FEATURES, to_matrix
from healthcare.features.phase2 import VARIANT_TYPES
from healthcare.models import calibration, dual_tower
from healthcare.models.export import (
    LIGHTGBM_INPUT,
    OPSET,
    latency,
    lightgbm_graph,
    merge,
    session,
)
from healthcare.models.joint import load_joint, per_type
from healthcare.models.joint_tower import ESM_COLUMNS

ONNX_NAME = "all_types_v2.onnx"
OUTPUTS = ["probability", "lightgbm", "dual_tower"]
INPUTS = ["features", "protein", "protein_mask"]


def type_parameters(final: dict) -> tuple[list[float], list[float], list[float]]:
    weights, slopes, biases = [], [], []
    for consequence in VARIANT_TYPES:
        weights.append(0.5 if final["routes"].get(consequence) == "average" else 0.0)
        group = all_types.calibration_group(consequence, final["calibration_groups"])
        slope, bias = final["platt"][group]
        slopes.append(slope)
        biases.append(bias)
    return weights, slopes, biases


class RoutedHead(nn.Module):
    def __init__(
        self,
        model: dual_tower.DualTower,
        encoder: dual_tower.TabularEncoder,
        feature_names: list[str],
        final: dict,
    ):
        super().__init__()
        self.model = model
        self.numeric_index = [feature_names.index(c) for c in encoder.columns]
        self.flag_index = [encoder.columns.index(c) for c in encoder.flag_columns]
        self.amino_acid_index = [feature_names.index(c) for c in CATEGORICAL_FEATURES]
        self.type_index = feature_names.index(dual_tower.VARIANT_TYPE)
        self.types = encoder.variant_types
        for name, columns in (
            ("log10_mask", dual_tower.LOG10_COLUMNS),
            ("log1p_mask", dual_tower.LOG1P_COLUMNS),
        ):
            self.register_buffer(name, torch.tensor([c in columns for c in encoder.columns]))
        for name in ("medians", "means", "stds"):
            values = np.asarray(getattr(encoder, name))
            self.register_buffer(name, torch.tensor(values, dtype=torch.float32))
        weights, slopes, biases = type_parameters(final)
        for name, values in (("weights", weights), ("slopes", slopes), ("biases", biases)):
            self.register_buffer(name, torch.tensor(values, dtype=torch.float32))

    def forward(self, features, protein, protein_mask, lightgbm_probabilities):
        raw = features[:, self.numeric_index]
        transformed = torch.where(
            self.log10_mask,
            torch.log10(raw + 1e-6),
            torch.where(self.log1p_mask, torch.log1p(raw), raw),
        )
        missing = torch.isnan(transformed)
        scaled = (torch.where(missing, self.medians, transformed) - self.means) / self.stds
        codes = features[:, self.type_index].long()
        one_hot = nn.functional.one_hot(codes, self.types).float()
        numeric = torch.cat([scaled, missing[:, self.flag_index].float(), one_hot], 1)
        amino = features[:, self.amino_acid_index]
        amino_acids = torch.where(
            torch.isnan(amino), torch.full_like(amino, dual_tower.AA_VOCAB - 1), amino
        ).long()
        tower = torch.sigmoid(self.model(numeric, amino_acids, protein, protein_mask))
        boosted = lightgbm_probabilities[:, 1]
        weight = self.weights[codes]
        blended = (weight * tower + (1 - weight) * boosted).clamp(
            all_types.EPSILON, 1 - all_types.EPSILON
        )
        logit = torch.log(blended / (1 - blended))
        probability = torch.sigmoid(self.slopes[codes] * logit + self.biases[codes])
        return probability, boosted, tower


def head_graph(head: RoutedHead, width: int, protein_width: int) -> onnx.ModelProto:
    head.eval()
    example = (
        torch.zeros(2, width),
        torch.zeros(2, protein_width),
        torch.ones(2),
        torch.full((2, 2), 0.5),
    )
    buffer = io.BytesIO()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        torch.onnx.export(
            head,
            example,
            buffer,
            dynamo=False,
            opset_version=OPSET,
            input_names=[*INPUTS, LIGHTGBM_INPUT],
            output_names=OUTPUTS,
            dynamic_axes={name: {0: "batch"} for name in [*INPUTS, LIGHTGBM_INPUT, *OUTPUTS]},
        )
    return onnx.load_from_string(buffer.getvalue())


def load_tower(path: Path) -> tuple[dual_tower.DualTower, dict]:
    checkpoint = torch.load(path, weights_only=False)
    encoder, config = checkpoint["encoder"], checkpoint["config"]
    model = dual_tower.DualTower(
        encoder.width, len(dual_tower.ESM_FORMS[config["esm_form"]]), config
    )
    model.load_state_dict(checkpoint["state"])
    return model.float().cpu().eval(), checkpoint


def build(booster: lgb.Booster, model, encoder, config: dict, final: dict) -> onnx.ModelProto:
    feature_names = booster.feature_name()
    missing = set(encoder.columns) - set(feature_names)
    if missing:
        raise ValueError(f"dual tower columns absent from LightGBM: {sorted(missing)}")
    head = RoutedHead(model, encoder, feature_names, final)
    protein_width = len(dual_tower.ESM_FORMS[config["esm_form"]]) * 1280
    merged = merge(head_graph(head, len(feature_names), protein_width), lightgbm_graph(booster))
    metadata = {
        "feature_names": json.dumps(feature_names),
        "esm_blocks": json.dumps([variant_inputs.BLOCKS[i] for i in ESM_COLUMNS]),
        "thresholds": json.dumps(final["thresholds"]),
        "decisions": json.dumps(calibration.DECISIONS),
        "variant_types": json.dumps(VARIANT_TYPES),
        "routes": json.dumps(final["routes"]),
    }
    for key, value in metadata.items():
        merged.metadata_props.add(key=key, value=value)
    return merged


def onnx_inputs(feature_names: list[str], table, blocks_path: Path) -> dict:
    protein = variant_inputs.load(blocks_path, table, ESM_COLUMNS)
    return {
        "features": to_matrix(table, feature_names).astype(np.float32),
        "protein": np.ascontiguousarray(protein.reshape(table.height, -1), dtype=np.float32),
        "protein_mask": table["protein_position"].is_not_null().to_numpy().astype(np.float32),
    }


@torch.inference_mode()
def tower_reference(model, encoder, inputs: dict, table) -> np.ndarray:
    numeric, amino_acids = encoder.encode(table)
    outputs = []
    for start in range(0, table.height, 8192):
        end = start + 8192
        logits = model(
            torch.from_numpy(numeric[start:end]),
            torch.from_numpy(amino_acids[start:end]),
            torch.from_numpy(inputs["protein"][start:end]),
            torch.from_numpy(inputs["protein_mask"][start:end]),
        )
        outputs.append(torch.sigmoid(logits).numpy().astype(np.float64))
    return np.concatenate(outputs)


def run(splits_dir: Path, blocks_path: Path, models_dir: Path, reports_dir: Path) -> dict:
    started = time.monotonic()
    final = json.loads((models_dir / "v2_final.json").read_text())
    booster = lgb.Booster(model_file=str(models_dir / final["members"]["lightgbm"]))
    model, checkpoint = load_tower(models_dir / final["members"]["dual_tower"])
    exported = build(booster, model, checkpoint["encoder"], checkpoint["config"], final)
    path = models_dir / ONNX_NAME
    onnx.save(exported, str(path))

    validation = load_joint(splits_dir, "validation")
    lightgbm = all_types.lightgbm_predict(models_dir / final["members"]["lightgbm"], validation)
    inputs = onnx_inputs(booster.feature_name(), validation, blocks_path)
    tower = tower_reference(model, checkpoint["encoder"], inputs, validation)
    raw = all_types.blend(validation, lightgbm, tower, final["routes"])
    expected = {
        "probability": all_types.calibrate(
            validation, raw, final["platt"], final["calibration_groups"]
        ),
        "lightgbm": lightgbm,
        "dual_tower": tower,
    }
    sess = session(exported)
    actual = dict(zip(OUTPUTS, sess.run(OUTPUTS, inputs), strict=True))
    decisions = int(
        (
            calibration.decide(expected["probability"]) != calibration.decide(actual["probability"])
        ).sum()
    )
    report = {
        "model": "v2 all coding types (LightGBM + dual tower, per-type route + Platt), ONNX",
        "test_set_opened_before_export": True,
        "model_unchanged_since_test": True,
        "file": ONNX_NAME,
        "megabytes": round(path.stat().st_size / 1e6, 2),
        "versions": {"onnx": onnx.__version__, "onnxruntime": ort.__version__},
        "inputs": {i.name: str(i.shape) for i in sess.get_inputs()},
        "parity_rows": validation.height,
        "parity_max_abs_diff": {
            name: round(float(np.abs(expected[name] - actual[name]).max()), 7) for name in OUTPUTS
        },
        "parity_decision_mismatches": decisions,
        "validation_onnx_by_type": {
            t: {k: m.get(k) for k in ("n", "roc_auc", "ece")}
            for t, m in per_type(validation, actual["probability"]).items()
        },
        "latency_cpu_single_variant": latency(sess, inputs),
        "seconds": round(time.monotonic() - started, 1),
    }
    (reports_dir / "onnx_export_v2.json").write_text(json.dumps(report, indent=2))
    return report
