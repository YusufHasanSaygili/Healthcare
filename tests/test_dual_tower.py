import random

import numpy as np
import optuna
import polars as pl
import pytest
import torch

from healthcare.esm.variant_inputs import variant_blocks
from healthcare.models import dual_tower, sweep

CONFIG = {
    "hidden": 16,
    "tabular_depth": 2,
    "fusion_depth": 1,
    "aa_dim": 4,
    "dropout": 0.1,
    "esm_form": "window3",
}


def table():
    return pl.DataFrame(
        {
            "f_gnomad_af": [0.0, 0.01, None, 0.5],
            "f_phylop": [1.0, 2.0, 3.0, None],
            "f_aa_ref": ["R", "W", "U", "A"],
            "f_aa_alt": ["W", "R", "C", "*"],
            "label": [1, 0, 1, 0],
        }
    )


def test_encoder_fits_on_given_rows_and_flags_missing_values():
    encoder = dual_tower.TabularEncoder.fit(table())
    numeric, amino_acids = encoder.encode(table())
    assert encoder.columns == ["f_gnomad_af", "f_phylop"]
    assert encoder.flag_columns == ["f_gnomad_af", "f_phylop"]
    assert numeric.shape == (4, encoder.width) and not np.isnan(numeric).any()
    assert numeric[2, 2] == 1.0 and numeric[3, 3] == 1.0
    assert amino_acids.shape == (4, 2)
    dropped = dual_tower.TabularEncoder.fit(table(), dropped=("f_gnomad_af",))
    assert dropped.columns == ["f_phylop"]


def test_forward_shapes_and_protein_mask():
    torch.manual_seed(0)
    model = dual_tower.DualTower(5, 2, CONFIG)
    numeric = torch.randn(3, 5)
    amino_acids = torch.randint(0, dual_tower.AA_VOCAB, (3, 2))
    protein = torch.randn(3, 2 * 1280)
    model.eval()
    on = model(numeric, amino_acids, protein, torch.ones(3))
    off = model(numeric, amino_acids, protein, torch.zeros(3))
    off_other_protein = model(numeric, amino_acids, torch.randn(3, 2 * 1280), torch.zeros(3))
    assert on.shape == (3,)
    assert not torch.allclose(on, off)
    assert torch.allclose(off, off_other_protein)
    no_esm = dual_tower.DualTower(5, 0, CONFIG | {"esm_form": "none"})
    assert no_esm(numeric, amino_acids, torch.zeros(3, 0), torch.ones(3)).shape == (3,)


def test_variant_blocks_clip_windows_at_protein_edges():
    embeddings = np.arange(10 * 1280, dtype=np.float32).reshape(10, 1280).astype(np.float16)
    blocks = variant_blocks(embeddings, 1)
    assert blocks.shape == (4, 1280)
    assert np.allclose(blocks[0], embeddings[0])
    assert np.allclose(blocks[1], embeddings[:4].astype(np.float32).mean(0), rtol=1e-2)
    assert np.allclose(blocks[3], embeddings.astype(np.float32).mean(0), rtol=1e-2)


def test_random_configs_are_inside_the_search_space():
    rng = random.Random(0)
    for _ in range(50):
        config = sweep.random_config(rng)
        for name, choices in sweep.SPACE.items():
            assert config[name] in choices
        for name, (low, high) in sweep.RANGES.items():
            assert low <= config[name] <= high


def test_enqueued_random_config_is_what_the_trial_suggests(tmp_path):
    study = sweep.open_study(tmp_path / "journal.log")
    config = sweep.random_config(random.Random(1))
    study.enqueue_trial(config)
    trial = study.ask()
    assert sweep.suggest(trial) == pytest.approx(config)
    study.tell(trial, 0.5)
    assert study.trials[0].state == optuna.trial.TrialState.COMPLETE
    assert len(sweep.open_study(tmp_path / "journal.log").trials) == 1


def test_esm_ablation_drops_every_esm_scalar():
    frame = table().with_columns(pl.lit(0.5).alias("f_esm_llr"), pl.lit(1.0).alias("f_esm_entropy"))
    encoder = dual_tower.TabularEncoder.fit(frame, dropped=dual_tower.ESM_SCALARS)
    assert not any(c.startswith("f_esm_") for c in encoder.columns)


def test_checkpoint_loads_with_weights_only(tmp_path):
    encoder = dual_tower.TabularEncoder.fit(table())
    payload = {
        key: value.tolist() if isinstance(value, np.ndarray) else value
        for key, value in vars(encoder).items()
    }
    torch.save({"encoder": payload}, tmp_path / "model.pt")
    loaded = torch.load(tmp_path / "model.pt", weights_only=True)
    assert loaded["encoder"]["columns"] == encoder.columns


def test_enqueued_trials_are_flagged_random(tmp_path):
    study = sweep.open_study(tmp_path / "journal.log")
    study.enqueue_trial(sweep.random_config(random.Random(2)))
    assert sweep.is_enqueued(study, study.ask())
    assert not sweep.is_enqueued(study, study.ask())
