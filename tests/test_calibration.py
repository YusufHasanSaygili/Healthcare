import numpy as np
import pytest

from healthcare.models import calibration


def synthetic(n=4000, seed=0):
    rng = np.random.default_rng(seed)
    labels = (rng.random(n) < 0.3).astype(int)
    scores = np.clip(rng.normal(0.25 + 0.5 * labels, 0.15), 0.001, 0.999)
    genes = np.array([f"G{i % 400}" for i in range(n)])
    return scores, labels, genes


@pytest.mark.parametrize("method", ["platt", "isotonic"])
def test_calibrator_is_monotonic_and_reduces_brier(method):
    scores, labels, _ = synthetic()
    calibrator = calibration.fit_calibrator(method, scores, labels)
    calibrated = calibration.calibrate(calibrator, scores)
    order = np.argsort(scores)
    assert np.all(np.diff(calibrated[order]) >= -1e-12)
    assert ((calibrated - labels) ** 2).mean() < ((scores - labels) ** 2).mean()
    edges = calibration.calibrate(calibrator, np.array([0.0, 1.0]))
    assert np.all((edges >= 0) & (edges <= 1))


def test_decide_maps_bands_inclusively():
    decisions = calibration.decide(np.array([0.05, 0.1, 0.5, 0.9, 0.95]))
    assert decisions.tolist() == [0, 0, 1, 2, 2]


def test_gene_folds_keep_each_gene_in_one_fold():
    _, _, genes = synthetic()
    folds = calibration.gene_folds(genes)
    for gene in np.unique(genes):
        assert np.unique(folds[genes == gene]).size == 1
    assert np.unique(folds).size == calibration.FOLDS


def test_cross_fit_predicts_every_row_from_other_genes():
    scores, labels, genes = synthetic()
    folds = calibration.gene_folds(genes)
    calibrated = calibration.cross_fit("platt", scores, labels, folds)
    assert calibrated.shape == scores.shape and np.isfinite(calibrated).all()
    held = folds == 0
    alone = calibration.fit_calibrator("platt", scores[~held], labels[~held])
    assert np.allclose(calibrated[held], calibration.calibrate(alone, scores[held]))


def test_calibrated_bands_hold_their_per_variant_precision():
    scores, labels, _ = synthetic(20000)
    calibrator = calibration.fit_calibrator("isotonic", scores, labels)
    summary = calibration.decision_summary(calibration.calibrate(calibrator, scores), labels)
    assert summary["likely_pathogenic"]["pathogenic_fraction"] >= 0.9
    assert summary["likely_benign"]["pathogenic_fraction"] <= 0.1
    shares = sum(summary[name]["share"] for name in calibration.DECISIONS)
    assert shares == pytest.approx(1.0, abs=1e-3)
