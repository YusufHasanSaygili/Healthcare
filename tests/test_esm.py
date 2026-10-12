import numpy as np
import pytest

from healthcare.esm import model
from healthcare.esm.trial import parse_protein_change


@pytest.mark.parametrize("length", [10, 1022, 1023, 2500, 34350])
def test_windows_cover_every_residue_within_model_limit(length):
    spans = model.windows(length)
    assert spans[0][0] == 0 and spans[-1][1] == length
    assert all(end - start <= model.MAX_RESIDUES for start, end in spans)
    chosen = model.window_assignment(length, spans)
    for position, window in enumerate(chosen):
        start, end = spans[window]
        assert start <= position < end


def test_long_protein_uses_window_with_most_context():
    spans = model.windows(2000)
    chosen = model.window_assignment(2000, spans)
    start, end = spans[chosen[1000]]
    assert min(1000 - start, end - 1000) > 400


def test_log_likelihood_ratio_and_entropy():
    log_probs = np.log(np.full((3, 20), 0.01, dtype=np.float32))
    log_probs[1, model.AMINO_ACIDS.index("R")] = np.log(0.81)
    assert model.log_likelihood_ratio(log_probs, 2, "R", "W") == pytest.approx(np.log(0.01 / 0.81))
    assert model.log_likelihood_ratio(log_probs, 1, "A", "A") == 0
    assert model.entropy(log_probs).shape == (3,)


def test_parse_protein_change():
    assert parse_protein_change("p.Arg1699Trp") == ("R", 1699, "W")
