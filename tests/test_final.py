import numpy as np
import polars as pl

from healthcare.evaluate import final


def test_evaluation_sets_cover_test_temporal_and_each_subset():
    test = pl.DataFrame({"label": [0, 1, 1]})
    temporal = pl.DataFrame({"temporal_subset": ["train_gene", "test_gene", "test_gene", "x"]})
    sets = final.evaluation_sets(test, temporal)
    assert list(sets) == [
        "test",
        "temporal_all",
        "temporal_test_gene",
        "temporal_train_gene",
        "temporal_x",
    ]
    assert sets["test"][0] == "test" and sets["test"][1].all()
    assert sets["temporal_test_gene"][1].tolist() == [False, True, True, False]


def test_score_handles_single_class_subsets():
    result = final.score(np.array([1, 1]), np.array([0.2, 0.9]), np.array(["A", "B"]))
    assert result["note"] == "single class"


def test_markdown_states_variant_type_test_status_and_every_set():
    metrics = final.score(
        np.array([0, 1] * 20), np.linspace(0.05, 0.95, 40), np.array(list("ABCDEFGHIJ") * 4)
    )
    decisions = {
        name: {"share": 0.5, "n": 1, "pathogenic_fraction": None}
        for name in ("likely_benign", "uncertain", "likely_pathogenic")
    }
    report = {
        "model": "m",
        "onnx": "m.onnx",
        "model_changed_after_test": False,
        "ablations": {"note": "n"},
        "temporal_note": "t",
        "results": {
            "test": {
                "variants": 40,
                "genes": 10,
                "pathogenic": 20,
                "models": {"ensemble_calibrated": metrics, "single": {"n": 2}},
                "decisions": decisions,
            }
        },
    }
    text = final.markdown(report)
    assert "missense SNVs only" in text and "**no**" in text
    assert "## test" in text and "| single | – |" in text
    assert "not for clinical diagnosis" in text
