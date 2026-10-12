import numpy as np
import polars as pl
import pytest

from healthcare.evaluate import metrics
from healthcare.models import baseline


def test_metrics_on_perfect_and_random_scores():
    labels = np.array([0, 0, 1, 1] * 50)
    perfect = metrics.classification_metrics(labels, labels.astype(float) * 0.9 + 0.05)
    assert perfect["roc_auc"] == 1.0 and perfect["mcc"] == 1.0
    assert metrics.expected_calibration_error(labels, np.full(labels.size, 0.5)) == pytest.approx(
        0.0
    )
    noisy = np.random.default_rng(0).random(labels.size)
    assert 0.3 < metrics.classification_metrics(labels, noisy)["roc_auc"] < 0.7


def test_gene_bootstrap_interval_contains_point_estimate():
    rng = np.random.default_rng(1)
    labels = rng.integers(0, 2, 400)
    scores = (labels * 0.6 + rng.random(400) * 0.6) / 1.2
    genes = np.repeat(np.arange(40), 10)
    low, high = metrics.gene_bootstrap(labels, scores, genes, rounds=100)["roc_auc"]
    point = metrics.classification_metrics(labels, scores)["roc_auc"]
    assert low <= point <= high


def test_feature_columns_and_ablation_drops():
    table = pl.DataFrame(
        {
            "gene_name": ["A"],
            "label": [1],
            "f_gnomad_af": [0.1],
            "f_gene_pli": [0.5],
            "f_protein_length": [100],
            "f_phylop": [2.0],
            "f_aa_ref": ["R"],
        }
    )
    assert baseline.feature_columns(table) == [
        "f_gnomad_af",
        "f_gene_pli",
        "f_protein_length",
        "f_phylop",
        "f_aa_ref",
    ]
    kept = baseline.feature_columns(table, baseline.ABLATIONS["no_frequency_no_gene_level"])
    assert kept == ["f_phylop", "f_aa_ref"]
    assert not {"gene_name", "label"} & set(baseline.feature_columns(table))


def test_amino_acids_are_encoded_with_a_fixed_table():
    table = pl.DataFrame({"f_aa_ref": ["A", "W", "U"], "f_x": [1, None, 3]})
    matrix = baseline.to_matrix(table, ["f_aa_ref", "f_x"])
    assert matrix[:, 0].tolist() == [0.0, 17.0, 20.0]
    assert np.isnan(matrix[1, 1])


def test_early_stopping_genes_are_deterministic_and_gene_grouped():
    train = pl.DataFrame({"gene_name": [f"G{i}" for i in range(100)]})
    first = baseline.early_stopping_genes(train)
    assert first == baseline.early_stopping_genes(train)
    assert len(first) == 10
