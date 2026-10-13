import io
import json
import statistics
import time
import warnings
from pathlib import Path

import lightgbm as lgb
import numpy as np
import onnx
import onnxruntime as ort
import polars as pl
import torch
from onnxmltools.convert import convert_lightgbm
from onnxmltools.convert.common.data_types import FloatTensorType
from torch import nn

from healthcare.esm import variant_inputs
from healthcare.evaluate.metrics import classification_metrics, rounded
from healthcare.features.missense import CATEGORICAL_FEATURES
from healthcare.models import calibration, dual_tower
from healthcare.models.baseline import to_matrix
from healthcare.split.genes import load_split

OPSET = 15
ONNX_NAME = "missense_v1.onnx"
OUTPUTS = ["probability", "lightgbm", "dual_tower"]
LIGHTGBM_INPUT = "lightgbm_probabilities"
LATENCY_RUNS = 200


class EnsembleHead(nn.Module):
    def __init__(
        self,
        model: dual_tower.DualTower,
        encoder: dual_tower.TabularEncoder,
        feature_names: list[str],
        calibrator: dict,
    ):
        super().__init__()
        if calibrator["method"] != "platt":
            raise ValueError("only a Platt calibrator can be exported")
        self.model = model
        self.numeric_index = [feature_names.index(c) for c in encoder.columns]
        self.flag_index = [encoder.columns.index(c) for c in encoder.flag_columns]
        self.amino_acid_index = [feature_names.index(c) for c in CATEGORICAL_FEATURES]
        for name, columns in (
            ("log10_mask", dual_tower.LOG10_COLUMNS),
            ("log1p_mask", dual_tower.LOG1P_COLUMNS),
        ):
            self.register_buffer(name, torch.tensor([c in columns for c in encoder.columns]))
        for name in ("medians", "means", "stds"):
            values = np.asarray(getattr(encoder, name))
            self.register_buffer(name, torch.tensor(values, dtype=torch.float32))
        self.slope = float(calibrator["slope"])
        self.bias = float(calibrator["bias"])

    def forward(self, features, protein, lightgbm_probabilities):
        raw = features[:, self.numeric_index]
        transformed = torch.where(
            self.log10_mask,
            torch.log10(raw + 1e-6),
            torch.where(self.log1p_mask, torch.log1p(raw), raw),
        )
        missing = torch.isnan(transformed)
        scaled = (torch.where(missing, self.medians, transformed) - self.means) / self.stds
        numeric = torch.cat([scaled, missing[:, self.flag_index].float()], 1)
        codes = features[:, self.amino_acid_index]
        amino_acids = torch.where(
            torch.isnan(codes), torch.full_like(codes, dual_tower.AA_VOCAB - 1), codes
        ).long()
        ones = torch.ones_like(features[:, 0])
        tower = torch.sigmoid(self.model(numeric, amino_acids, protein, ones))
        boosted = lightgbm_probabilities[:, 1]
        mean = ((tower + boosted) / 2).clamp(calibration.EPSILON, 1 - calibration.EPSILON)
        logit = torch.log(mean / (1 - mean))
        return torch.sigmoid(self.slope * logit + self.bias), boosted, tower


def lightgbm_graph(booster: lgb.Booster) -> onnx.ModelProto:
    width = booster.num_feature()
    model = convert_lightgbm(
        booster,
        initial_types=[("features", FloatTensorType([None, width]))],
        zipmap=False,
        target_opset=OPSET,
    )
    probabilities = [o.name for o in model.graph.output if o.name != "label"]
    if len(probabilities) != 1:
        raise RuntimeError(f"unexpected LightGBM ONNX outputs: {probabilities}")
    for node in model.graph.node:
        node.output[:] = [LIGHTGBM_INPUT if o == probabilities[0] else o for o in node.output]
    for output in model.graph.output:
        if output.name == probabilities[0]:
            output.name = LIGHTGBM_INPUT
    return model


def head_graph(head: EnsembleHead, width: int, protein_width: int) -> onnx.ModelProto:
    head.eval()
    example = (
        torch.zeros(2, width),
        torch.zeros(2, protein_width),
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
            input_names=["features", "protein", LIGHTGBM_INPUT],
            output_names=OUTPUTS,
            dynamic_axes={
                name: {0: "batch"} for name in ["features", "protein", LIGHTGBM_INPUT, *OUTPUTS]
            },
        )
    return onnx.load_from_string(buffer.getvalue())


def merge(head: onnx.ModelProto, boosted: onnx.ModelProto) -> onnx.ModelProto:
    graph = head.graph
    taken = {i.name for i in graph.initializer} | {o for n in graph.node for o in n.output}
    clashes = {i.name for i in boosted.graph.initializer} | {
        o for n in boosted.graph.node for o in n.output if o != LIGHTGBM_INPUT
    }
    if taken & clashes:
        raise RuntimeError(f"name clash while merging graphs: {sorted(taken & clashes)[:5]}")
    inputs = [i for i in graph.input if i.name != LIGHTGBM_INPUT]
    del graph.input[:]
    graph.input.extend(inputs)
    nodes = list(boosted.graph.node) + list(graph.node)
    del graph.node[:]
    graph.node.extend(nodes)
    graph.initializer.extend(boosted.graph.initializer)
    domains = {o.domain for o in head.opset_import}
    head.opset_import.extend(o for o in boosted.opset_import if o.domain not in domains)
    head.ir_version = max(head.ir_version, boosted.ir_version)
    onnx.checker.check_model(head, full_check=True)
    return head


def build(
    booster: lgb.Booster,
    model: dual_tower.DualTower,
    encoder: dual_tower.TabularEncoder,
    config: dict,
    calibrator: dict,
) -> onnx.ModelProto:
    feature_names = booster.feature_name()
    missing = set(encoder.columns) - set(feature_names)
    if missing:
        raise ValueError(f"dual tower columns absent from LightGBM: {sorted(missing)}")
    head = EnsembleHead(model.float().cpu(), encoder, feature_names, calibrator)
    protein_width = len(dual_tower.ESM_FORMS[config["esm_form"]]) * 1280
    merged = merge(
        head_graph(head, len(feature_names), protein_width),
        lightgbm_graph(booster),
    )
    metadata = {
        "feature_names": json.dumps(feature_names),
        "esm_blocks": json.dumps(
            [variant_inputs.BLOCKS[i] for i in dual_tower.ESM_FORMS[config["esm_form"]]]
        ),
        "thresholds": json.dumps(calibration.THRESHOLDS),
        "decisions": json.dumps(calibration.DECISIONS),
    }
    for key, value in metadata.items():
        merged.metadata_props.add(key=key, value=value)
    return merged


@torch.inference_mode()
def reference(
    booster: lgb.Booster,
    model: dual_tower.DualTower,
    encoder: dual_tower.TabularEncoder,
    calibrator: dict,
    table: pl.DataFrame,
    protein: np.ndarray,
) -> dict:
    boosted = booster.predict(to_matrix(table, booster.feature_name()))
    numeric, amino_acids = encoder.encode(table)
    model = model.float().cpu().eval()
    logits = model(
        torch.from_numpy(numeric),
        torch.from_numpy(amino_acids),
        torch.from_numpy(protein),
        torch.ones(len(table)),
    )
    tower = torch.sigmoid(logits).numpy().astype(np.float64)
    return {
        "probability": calibration.calibrate(calibrator, (boosted + tower) / 2),
        "lightgbm": boosted,
        "dual_tower": tower,
    }


def session(model: onnx.ModelProto) -> ort.InferenceSession:
    options = ort.SessionOptions()
    options.intra_op_num_threads = 1
    return ort.InferenceSession(
        model.SerializeToString(), options, providers=["CPUExecutionProvider"]
    )


def onnx_inputs(booster: lgb.Booster, table: pl.DataFrame, protein: np.ndarray) -> dict:
    return {
        "features": to_matrix(table, booster.feature_name()).astype(np.float32),
        "protein": np.ascontiguousarray(protein, dtype=np.float32),
    }


def predict(sess: ort.InferenceSession, inputs: dict) -> dict:
    return dict(zip(OUTPUTS, sess.run(OUTPUTS, inputs), strict=True))


def parity(expected: dict, actual: dict, decisions: bool = True) -> dict:
    result = {
        name: round(float(np.abs(expected[name] - actual[name]).max()), 7) for name in OUTPUTS
    }
    if decisions:
        result["decision_mismatches"] = int(
            (
                calibration.decide(expected["probability"])
                != calibration.decide(actual["probability"])
            ).sum()
        )
    return result


def latency(sess: ort.InferenceSession, inputs: dict) -> dict:
    single = {name: value[:1] for name, value in inputs.items()}
    for _ in range(20):
        sess.run(OUTPUTS, single)
    timings = []
    for _ in range(LATENCY_RUNS):
        started = time.perf_counter()
        sess.run(OUTPUTS, single)
        timings.append((time.perf_counter() - started) * 1000)
    timings.sort()
    return {
        "runs": LATENCY_RUNS,
        "median_ms": round(statistics.median(timings), 3),
        "p95_ms": round(timings[int(0.95 * len(timings)) - 1], 3),
        "threads": 1,
    }


def run(splits_dir: Path, esm_blocks: Path, models_dir: Path, reports_dir: Path) -> dict:
    started = time.monotonic()
    booster = lgb.Booster(model_file=str(models_dir / "lightgbm_full.txt"))
    model, meta = calibration.load_dual_tower(models_dir / "dual_tower.pt", "cpu")
    calibrator = json.loads((models_dir / "calibration.json").read_text())["calibrator"]
    exported = build(booster, model, meta["encoder"], meta["config"], calibrator)
    path = models_dir / ONNX_NAME
    onnx.save(exported, str(path))

    validation = load_split(splits_dir, "validation")
    blocks = dual_tower.ESM_FORMS[meta["config"]["esm_form"]]
    protein = variant_inputs.load(esm_blocks, validation)[:, blocks].reshape(len(validation), -1)
    protein = protein.astype(np.float32)
    expected = reference(booster, model, meta["encoder"], calibrator, validation, protein)
    sess = session(exported)
    inputs = onnx_inputs(booster, validation, protein)
    actual = predict(sess, inputs)
    labels = validation["label"].to_numpy()

    report = {
        "model": "Missense ensemble (LightGBM + dual tower + Platt), ONNX",
        "test_set_opened": False,
        "file": ONNX_NAME,
        "megabytes": round(path.stat().st_size / 1e6, 2),
        "opset": {o.domain or "ai.onnx": o.version for o in exported.opset_import},
        "versions": {
            "onnx": onnx.__version__,
            "onnxruntime": ort.__version__,
            "torch": torch.__version__,
            "lightgbm": lgb.__version__,
        },
        "inputs": {i.name: str(i.shape) for i in sess.get_inputs()},
        "parity_rows": len(validation),
        "parity_max_abs_diff": parity(expected, actual),
        "parity_rows_over_1e-4": {
            name: int((np.abs(expected[name] - actual[name]) > 1e-4).sum()) for name in OUTPUTS
        },
        "parity_note": "LightGBM differences come from float32 tree comparisons in ONNX: "
        "f_delta_* values that differ only by float64 rounding noise collapse to one float32 value",
        "validation_onnx": {
            "metrics": rounded(classification_metrics(labels, actual["probability"])),
            "decisions": calibration.decision_summary(actual["probability"], labels),
        },
        "latency_cpu_single_variant": latency(sess, inputs),
        "seconds": round(time.monotonic() - started, 1),
    }
    (reports_dir / "onnx_export.json").write_text(json.dumps(report, indent=2))
    return report
