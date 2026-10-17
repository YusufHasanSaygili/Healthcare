import numpy as np
import polars as pl
import pytest

from healthcare.evaluate import all_types


def test_blend_averages_only_routed_types():
    table = pl.DataFrame({"consequence": ["missense", "frameshift", "start_lost"]})
    lightgbm, tower = np.array([0.2, 0.9, 0.5]), np.array([0.4, 0.1, 0.1])
    routes = {"missense": "average", "frameshift": "lightgbm"}
    assert np.allclose(all_types.blend(table, lightgbm, tower, routes), [0.3, 0.9, 0.5])


def test_calibrate_uses_each_type_group_and_pools_the_rest():
    table = pl.DataFrame({"consequence": ["missense", "start_lost"]})
    platt = {"missense": [1.0, 0.0], all_types.POOLED: [1.0, 2.0]}
    calibrated = all_types.calibrate(table, np.array([0.5, 0.5]), platt, {"missense": "missense"})
    assert calibrated[0] == 0.5
    assert calibrated[1] == 1 / (1 + np.exp(-2.0))


def test_fit_platt_recovers_a_monotone_mapping():
    rng = np.random.default_rng(0)
    raw = rng.uniform(0.01, 0.99, 2000)
    labels = (rng.uniform(size=2000) < raw).astype(int)
    slope, _ = all_types.fit_platt(raw, labels)
    assert 0.7 < slope < 1.3


def test_decision_shares_follow_the_thresholds():
    shares = all_types.decision_shares(np.array([0.95, 0.5, 0.05, 0.02]))
    assert shares == {"likely_pathogenic": 0.25, "uncertain": 0.25, "likely_benign": 0.5}


def test_evaluate_refuses_a_second_opening(tmp_path):
    (tmp_path / "all_types_v2.json").write_text("{}")
    with pytest.raises(all_types.AlreadyOpenedError):
        all_types.evaluate(tmp_path, tmp_path, tmp_path / "b.npy", tmp_path, tmp_path)
