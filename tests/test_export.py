import json

import lightgbm as lgb
import numpy as np
import polars as pl
import pytest
import torch

from healthcare.models import calibration, dual_tower, export
from healthcare.models.baseline import AMINO_ACID_CODES, to_matrix

CONFIG = {
    "hidden": 16,
    "tabular_depth": 2,
    "fusion_depth": 1,
    "aa_dim": 4,
    "dropout": 0.0,
    "esm_form": "window3",
}
CALIBRATOR = {"method": "platt", "slope": 1.3, "bias": -0.2}


def table(n=400, seed=0):
    rng = np.random.default_rng(seed)
    letters = list("ACDEFGHIKLMNPQRSTVWY") + ["*", "U", "B"]
    labels = rng.integers(0, 2, n)
    af = np.where(rng.random(n) < 0.3, None, rng.random(n) * 0.01 * (1 - labels))
    phylop = rng.normal(2 + 3 * labels, 1.5)
    phylop = np.where(rng.random(n) < 0.1, None, phylop)
    return pl.DataFrame(
        {
            "f_gnomad_af": af.tolist(),
            "f_phylop": phylop.tolist(),
            "f_protein_length": rng.integers(50, 3000, n).astype(float),
            "f_aa_ref": rng.choice(letters, n).tolist(),
            "f_aa_alt": rng.choice(letters, n).tolist(),
            "label": labels,
        },
        schema_overrides={"f_gnomad_af": pl.Float64, "f_phylop": pl.Float64},
    )


@pytest.fixture(scope="module")
def models():
    data = table()
    columns = [c for c in data.columns if c.startswith("f_")]
    categorical = [columns.index(c) for c in ("f_aa_ref", "f_aa_alt")]
    booster = lgb.train(
        {"objective": "binary", "verbosity": -1, "num_leaves": 7, "min_child_samples": 5},
        lgb.Dataset(
            to_matrix(data, columns),
            data["label"].to_numpy(),
            feature_name=columns,
            categorical_feature=categorical,
        ),
        num_boost_round=20,
    )
    torch.manual_seed(0)
    encoder = dual_tower.TabularEncoder.fit(data)
    model = dual_tower.DualTower(encoder.width, 2, CONFIG).eval()
    return booster, model, encoder


def test_onnx_matches_native_lightgbm_and_pytorch(models):
    booster, model, encoder = models
    exported = export.build(booster, model, encoder, CONFIG, CALIBRATOR)
    data = table(300, seed=1)
    protein = np.random.default_rng(2).normal(size=(len(data), 2 * 1280)).astype(np.float32)
    expected = export.reference(booster, model, encoder, CALIBRATOR, data, protein)
    actual = export.predict(export.session(exported), export.onnx_inputs(booster, data, protein))
    diffs = export.parity(expected, actual, decisions=False)
    assert diffs["lightgbm"] < 1e-5
    assert diffs["dual_tower"] < 1e-5
    assert diffs["probability"] < 1e-5
    assert actual["probability"].shape == (len(data),)


def test_onnx_metadata_carries_feature_order_and_thresholds(models):
    booster, model, encoder = models
    exported = export.build(booster, model, encoder, CONFIG, CALIBRATOR)
    metadata = {p.key: json.loads(p.value) for p in exported.metadata_props}
    assert metadata["feature_names"] == booster.feature_name()
    assert metadata["thresholds"] == calibration.THRESHOLDS
    assert metadata["esm_blocks"] == ["position", "window3"]


def test_unknown_amino_acid_uses_the_same_code_in_both_members(models):
    booster, model, encoder = models
    data = table(5, seed=3).with_columns(pl.lit("B").alias("f_aa_alt"))
    assert "B" not in AMINO_ACID_CODES
    protein = np.zeros((len(data), 2 * 1280), dtype=np.float32)
    exported = export.build(booster, model, encoder, CONFIG, CALIBRATOR)
    expected = export.reference(booster, model, encoder, CALIBRATOR, data, protein)
    actual = export.predict(export.session(exported), export.onnx_inputs(booster, data, protein))
    assert export.parity(expected, actual, decisions=False)["dual_tower"] < 1e-5


def test_only_platt_calibrators_export(models):
    booster, model, encoder = models
    with pytest.raises(ValueError):
        export.build(booster, model, encoder, CONFIG, {"method": "isotonic", "x": [], "y": []})
