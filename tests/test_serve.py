import subprocess
import sys

import numpy as np
import polars as pl
import pytest

from healthcare.annotate.snv import Annotation
from healthcare.serve import pipeline

THRESHOLDS = {"likely_benign_max": 0.1, "likely_pathogenic_min": 0.9}


def annotation(consequence: str, gene: str = "GENE") -> Annotation:
    return Annotation("T1", gene, None, consequence, aa_ref="R", aa_alt="W", protein_position=5)


def test_serve_pipeline_never_imports_torch():
    code = "import sys, healthcare.serve.pipeline; print('torch' in sys.modules)"
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "False"


@pytest.mark.parametrize(
    ("annotations", "consequence", "reason"),
    [
        ([], "intergenic", "no MANE Select transcript"),
        ([annotation("missense", "A"), annotation("missense", "B")], "multiple_genes", "more than"),
        ([annotation("synonymous")], "synonymous", "phase 2"),
        ([annotation("canonical_splice")], "canonical_splice", "phase 2"),
        ([annotation("intronic")], "intronic", "permanently out of scope"),
    ],
)
def test_classify_explains_out_of_scope_variants(annotations, consequence, reason):
    result = pipeline.classify("1:100 A>G", annotations)
    assert isinstance(result, pipeline.OutOfScope)
    assert result.consequence == consequence and reason in result.reason


def test_classify_accepts_missense_and_rejects_ref_mismatch():
    missense = annotation("missense")
    assert pipeline.classify("1:100 A>G", [missense]) is missense
    with pytest.raises(ValueError, match="REF does not match"):
        pipeline.classify("1:100 A>G", [annotation("ref_mismatch")])


def test_decide_uses_inclusive_thresholds():
    decisions = [pipeline.decide(p, THRESHOLDS) for p in (0.05, 0.1, 0.5, 0.9, 0.99)]
    assert decisions == [
        "likely_benign",
        "likely_benign",
        "uncertain",
        "likely_pathogenic",
        "likely_pathogenic",
    ]


def test_explanation_replaces_one_group_at_a_time_with_background_values():
    names = sorted(c for columns in pipeline.FEATURE_GROUPS.values() for c in columns)
    features = np.zeros(len(names), dtype=np.float32)
    background = np.ones((3, len(names)), dtype=np.float32)
    rows = pipeline.explanation_rows(features, background, names)
    assert rows.shape == (1 + 3 * len(pipeline.FEATURE_GROUPS), len(names))
    assert not rows[0].any()
    for g, columns in enumerate(pipeline.FEATURE_GROUPS.values()):
        block = rows[1 + 3 * g : 4 + 3 * g]
        replaced = [names.index(c) for c in columns]
        assert (block[:, replaced] == 1).all()
        assert block.sum() == 3 * len(columns)


def test_contributions_are_log_odds_changes_against_the_background_mean():
    groups = len(pipeline.FEATURE_GROUPS)
    probabilities = np.full(1 + 2 * groups, 0.5)
    probabilities[0] = 0.8
    probabilities[1:3] = [0.2, 0.5]
    probabilities[3:5] = [0.9, 0.9]
    by_group = dict(pipeline.group_contributions(probabilities, 2))
    names = list(pipeline.FEATURE_GROUPS)
    high, low, mid, top = pipeline.log_odds(np.array([0.8, 0.2, 0.5, 0.9]))
    assert by_group[names[0]] == pytest.approx(high - (low + mid) / 2, abs=1e-3)
    assert by_group[names[1]] == pytest.approx(high - top, abs=1e-3)
    assert by_group[names[2]] == pytest.approx(high - mid, abs=1e-3)


def test_contributions_stay_informative_when_the_probability_saturates():
    groups = len(pipeline.FEATURE_GROUPS)
    probabilities = np.full(1 + groups, 1e-7)
    probabilities[1] = 0.3
    contributions = pipeline.group_contributions(probabilities, 1)
    assert contributions[0][0] == list(pipeline.FEATURE_GROUPS)[0]
    assert contributions[0][1] < -5
    assert all(np.isfinite(delta) for _, delta in contributions)


def test_background_sample_is_deterministic_and_bounded():
    table = pl.DataFrame({"f_aa_ref": ["R", "W", "A"], "f_x": [1.0, None, 3.0]})
    first = pipeline.background_sample(table, ["f_aa_ref", "f_x"], size=2)
    second = pipeline.background_sample(table, ["f_aa_ref", "f_x"], size=2)
    assert first.shape == (2, 2) and np.array_equal(first, second, equal_nan=True)
    assert pipeline.background_sample(table, ["f_x"], size=10).shape == (3, 1)
