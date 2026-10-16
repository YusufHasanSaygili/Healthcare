import numpy as np
import polars as pl

from healthcare.models import joint


def test_load_joint_tags_missense_and_drops_multi_gene_keys(tmp_path):
    key = {"chrom": ["1", "1"], "ref": ["A", "A"], "alt": ["G", "T"]}
    missense = pl.DataFrame(
        key
        | {"pos": [10, 20], "gene_name": ["G1", "G1"], "label": [1, 0], "split": "train"}
        | {"f_esm_llr": [-9.0, -1.0]}
    )
    phase2 = pl.DataFrame(
        {"chrom": ["1"], "pos": [30], "ref": ["A"], "alt": ["AT"], "gene_name": ["G1"]}
        | {"label": [1], "split": "train", "consequence": "frameshift", "f_variant_type": [6.0]}
    )
    missense.write_parquet(tmp_path / "missense_open.parquet")
    phase2.write_parquet(tmp_path / "phase2_open.parquet")
    missense.head(1).select("chrom", "pos", "ref", "alt").write_parquet(
        tmp_path / "multi_gene_keys.parquet"
    )
    table = joint.load_joint(tmp_path, "train")
    assert sorted(table["pos"].to_list()) == [20, 30]
    rows = dict(table.select("consequence", "f_variant_type").iter_rows())
    assert rows == {"missense": 0.0, "frameshift": 6.0}
    assert table.filter(pl.col("pos") == 30)["f_esm_llr"][0] is None


def test_per_type_skips_types_without_enough_of_each_class():
    labels = [1] * 30 + [0] * 30 + [1] * 25 + [0] * 2
    table = pl.DataFrame(
        {
            "consequence": ["missense"] * 60 + ["frameshift"] * 27,
            "label": labels,
            "gene_name": [f"G{i % 9}" for i in range(87)],
        }
    )
    probabilities = np.array(labels, dtype=float) * 0.8 + 0.1
    report = joint.per_type(table, probabilities)
    assert report["missense"]["roc_auc"] == 1.0
    assert report["frameshift"] == {"n": 27, "positives": 25, "eligible": False}
    assert joint.selection_score(report) == 1.0
