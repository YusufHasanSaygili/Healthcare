import numpy as np
import polars as pl
import pytest

from healthcare.esm import store
from healthcare.esm.model import AMINO_ACIDS
from healthcare.features.esm_features import ESM_FEATURES, esm_features, position_features


def log_probs(length=5):
    probabilities = np.full((length, 20), 0.01, dtype=np.float32)
    probabilities[:, AMINO_ACIDS.index("R")] = 0.81
    return np.log(probabilities).astype(np.float16)


def test_position_features_values():
    llr, wild_type, mutant, entropy, rank = position_features(log_probs(), 2, "R", "W")
    assert llr == pytest.approx(np.log(0.01 / 0.81), abs=0.01)
    assert wild_type == pytest.approx(np.log(0.81), abs=0.01)
    assert rank > 1 and entropy > 0
    assert position_features(log_probs(), 2, "R", "R")[4] == 1
    assert position_features(log_probs(), 9, "R", "W") == (None,) * 5
    assert position_features(log_probs(), 2, "U", "W") == (None,) * 5


def test_esm_features_keep_row_order_and_handle_missing_proteins(tmp_path):
    store.save(tmp_path, "T1", log_probs(), np.zeros((5, 4), dtype=np.float16))
    table = pl.DataFrame(
        {
            "transcript_id": ["T1", "MISSING", "T1"],
            "protein_position": [1, 3, 2],
            "aa_ref": ["R", "R", "R"],
            "aa_alt": ["W", "W", "R"],
        }
    )
    result = esm_features(table, tmp_path)
    assert result.columns[-5:] == ESM_FEATURES
    assert result["f_esm_llr"][0] < 0
    assert result["f_esm_llr"][1] is None
    assert result["f_esm_llr"][2] == 0
