import numpy as np

from healthcare.esm import store


def test_save_and_memory_mapped_load(tmp_path):
    log_probs = np.random.default_rng(0).random((7, 20)).astype(np.float16)
    embeddings = np.ones((7, 1280), dtype=np.float16)
    assert not store.is_done(tmp_path, "ENST1")
    store.save(tmp_path, "ENST1", log_probs, embeddings)
    assert store.is_done(tmp_path, "ENST1")
    loaded_log_probs, loaded_embeddings = store.load(tmp_path, "ENST1")
    assert isinstance(loaded_embeddings, np.memmap)
    assert np.array_equal(loaded_log_probs, log_probs)
    assert loaded_embeddings.shape == (7, 1280)
    assert not list(tmp_path.rglob("*.part"))
