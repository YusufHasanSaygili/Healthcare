import json

import polars as pl
import pytest

from healthcare.split import genes as gene_split


def variants(rows):
    return pl.DataFrame(
        rows, schema=["chrom", "pos", "ref", "alt", "gene_name", "label"], orient="row"
    )


def transcripts(genes):
    return pl.DataFrame(
        {
            "gene_name": [g for g, _, _ in genes],
            "chrom": ["1"] * len(genes),
            "cds_starts": [[s] for _, s, _ in genes],
            "cds_ends": [[e] for _, _, e in genes],
        }
    )


def synthetic(n_genes=300):
    genes = [(f"G{g}", g * 1000, g * 1000 + 100) for g in range(n_genes)]
    genes += [("UNLABELED", 999_000, 999_100)]
    rows = [
        ("1", g * 1000 + v, "A", "G", f"G{g}", int((g + v) % 3 == 0))
        for g in range(n_genes)
        for v in range(1 + g % 7)
    ]
    return transcripts(genes), variants(rows)


def test_gene_groups_merge_cds_overlap_and_shared_positions():
    tx = transcripts(
        [("UGT1A1", 100, 200), ("UGT1A4", 150, 300), ("FAR", 500, 600), ("FAR2", 700, 800)]
    )
    shared = variants([("1", 550, "A", "G", "FAR", 0), ("1", 550, "A", "C", "FAR2", 1)])
    groups = gene_split.gene_groups(tx, shared)
    assert groups["UGT1A1"] == groups["UGT1A4"]
    assert groups["FAR"] == groups["FAR2"]
    assert groups["UGT1A1"] != groups["FAR"]


def test_split_is_disjoint_balanced_deterministic_and_covers_all_genes():
    tx, labeled = synthetic()
    first = gene_split.build_split(tx, labeled)
    assert first.equals(gene_split.build_split(tx, labeled))
    assert set(first["gene_name"]) == set(tx["gene_name"])
    joined = labeled.join(first.select("gene_name", "split"), on="gene_name")
    assert set(gene_split.overlaps(joined).values()) == {0}
    summary = gene_split.split_summary(joined)
    for name, target in gene_split.TARGETS.items():
        assert summary[name]["variant_share"] == pytest.approx(target, abs=0.03)
    assert gene_split.cds_overlap_violations(tx, first) == 0


def write_inputs(tmp_path):
    tx, labeled = synthetic()
    features = tmp_path / "features"
    features.mkdir()
    labeled.write_parquet(features / "missense_train_labels.parquet")
    temporal = variants([("1", 5, "A", "T", "G0", 1), ("1", 999_050, "A", "T", "BRANDNEW", 0)])
    temporal.write_parquet(features / "missense_temporal_labels.parquet")
    tx_path = tmp_path / "transcripts.parquet"
    tx.write_parquet(tx_path)
    return features, tx_path, tmp_path / "splits", tmp_path / "report.json"


def test_build_pins_split_and_locks_test(tmp_path):
    features, tx_path, splits, report = write_inputs(tmp_path)
    summary = gene_split.build(features, tx_path, splits, report, force=False)
    assert summary["regenerated"]
    assert summary["temporal_variants_in_train_release"] == 0
    assert not gene_split.build(features, tx_path, splits, report, force=False)["regenerated"]
    assert set(gene_split.load_split(splits, "train")["split"]) == {"train"}
    with pytest.raises(gene_split.LockedSplitError):
        gene_split.load_split(splits, "test")
    with pytest.raises(gene_split.LockedSplitError):
        gene_split.load_split(splits, "temporal")
    assert set(gene_split.load_split(splits, "test", allow_locked=True)["split"]) == {"test"}
    assert "test" not in set(pl.read_parquet(splits / "missense_open.parquet")["split"])


def test_build_refuses_missing_or_tampered_pin(tmp_path):
    features, tx_path, splits, report = write_inputs(tmp_path)
    gene_split.build(features, tx_path, splits, report, force=False)
    genes = pl.read_parquet(splits / "genes.parquet")
    genes.head(5).write_parquet(splits / "genes.parquet")
    with pytest.raises(RuntimeError):
        gene_split.build(features, tx_path, splits, report, force=False)
    (splits / "genes.parquet").unlink()
    with pytest.raises(RuntimeError):
        gene_split.build(features, tx_path, splits, report, force=False)
    gene_split.build(features, tx_path, splits, report, force=True)
    pl.read_parquet(features / "missense_train_labels.parquet").head(50).write_parquet(
        features / "missense_train_labels.parquet"
    )
    with pytest.raises(RuntimeError):
        gene_split.build(features, tx_path, splits, report, force=False)
    assert json.loads(report.read_text())["regenerated"]


def test_temporal_subsets_and_shared_positions():
    genes = pl.DataFrame(
        {"gene_name": ["A", "B", "C"], "group": [0, 1, 2], "split": ["train", "validation", "test"]}
    )
    temporal = variants(
        [
            ("1", 10, "A", "G", "A", 1),
            ("1", 20, "A", "G", "B", 0),
            ("1", 30, "A", "G", "C", 1),
            ("1", 40, "A", "G", "NEW", 0),
        ]
    )
    train_positions = pl.DataFrame({"chrom": ["1"], "pos": [10]})
    labeled = gene_split.label_temporal(temporal, genes, train_positions).sort("pos")
    assert labeled["temporal_subset"].to_list() == [
        "train_gene",
        "validation_gene",
        "test_gene",
        "new_gene",
    ]
    assert labeled["shares_position_with_train"].to_list() == [True, False, False, False]


def test_pin_allows_new_feature_columns_but_not_label_changes(tmp_path):
    features, tx_path, splits, report = write_inputs(tmp_path)
    gene_split.build(features, tx_path, splits, report, force=False)
    train_path = features / "missense_train_labels.parquet"
    original = pl.read_parquet(train_path)
    original.with_columns(pl.lit(1.0).alias("f_new_feature")).write_parquet(train_path)
    assert not gene_split.build(features, tx_path, splits, report, force=False)["regenerated"]
    assert "f_new_feature" in gene_split.load_split(splits, "train").columns
    flipped = original.with_columns(
        pl.when(pl.int_range(pl.len()) == 0)
        .then(1 - pl.col("label"))
        .otherwise(pl.col("label"))
        .alias("label")
    )
    flipped.write_parquet(train_path)
    with pytest.raises(RuntimeError):
        gene_split.build(features, tx_path, splits, report, force=False)


def write_phase2(tmp_path, features):
    phase2 = pl.DataFrame(
        [
            ("1", g * 1000 + 50, "A", "AT", f"G{g}", int(g % 2 == 0), "frameshift")
            for g in range(0, 300, 3)
        ],
        schema=["chrom", "pos", "ref", "alt", "gene_name", "label", "consequence"],
        orient="row",
    )
    phase2.write_parquet(features / "phase2_train_labels.parquet")
    phase2.head(2).with_columns(pl.col("pos") + 1).write_parquet(
        features / "phase2_temporal_labels.parquet"
    )
    annotations = tmp_path / "annotations"
    annotations.mkdir()
    shared = pl.DataFrame(
        [("1", 1000, "A", "G", "G1", "missense"), ("1", 1000, "A", "G", "G2", "canonical_splice")],
        schema=["chrom", "pos", "ref", "alt", "gene_name", "consequence"],
        orient="row",
    )
    shared.write_parquet(annotations / "train_labels.parquet")
    shared.head(0).write_parquet(annotations / "temporal_labels.parquet")
    return annotations


def test_phase2_reuses_the_pinned_split_and_locks_test_rows(tmp_path):
    features, tx_path, splits, report = write_inputs(tmp_path)
    gene_split.build(features, tx_path, splits, report, force=False)
    annotations = write_phase2(tmp_path, features)
    phase2_report = tmp_path / "phase2.json"
    summary = gene_split.build_phase2(features, annotations, splits, report, phase2_report, False)
    assert set(summary["overlaps"].values()) == {0}
    assert summary["missense_cross_table"]["genes_in_multiple_splits"] == 0
    assert summary["multi_gene_keys"] == 1
    genes = dict(pl.read_parquet(splits / "genes.parquet").select("gene_name", "split").iter_rows())
    train = gene_split.load_split(splits, "train", table="phase2")
    assert all(genes[gene] == "train" for gene in train["gene_name"])
    with pytest.raises(gene_split.LockedSplitError):
        gene_split.load_split(splits, "test", table="phase2")
    assert "test" not in set(pl.read_parquet(splits / "phase2_open.parquet")["split"])


def test_phase2_refuses_changed_variants_or_a_tampered_gene_split(tmp_path):
    features, tx_path, splits, report = write_inputs(tmp_path)
    gene_split.build(features, tx_path, splits, report, force=False)
    annotations = write_phase2(tmp_path, features)
    phase2_report = tmp_path / "phase2.json"
    gene_split.build_phase2(features, annotations, splits, report, phase2_report, False)
    pl.read_parquet(features / "phase2_train_labels.parquet").head(10).write_parquet(
        features / "phase2_train_labels.parquet"
    )
    with pytest.raises(RuntimeError):
        gene_split.build_phase2(features, annotations, splits, report, phase2_report, False)
    pl.read_parquet(splits / "genes.parquet").head(5).write_parquet(splits / "genes.parquet")
    with pytest.raises(RuntimeError):
        gene_split.build_phase2(features, annotations, splits, report, phase2_report, True)


def test_cross_table_conflicts_are_counted():
    first = pl.DataFrame(
        {"chrom": ["1", "1"], "pos": [10, 20], "ref": "A", "alt": "G"}
        | {"gene_name": ["G1", "G2"], "split": ["train", "validation"]}
    )
    second = pl.DataFrame(
        {"chrom": ["1", "1"], "pos": [10, 30], "ref": ["A", "A"], "alt": ["G", "T"]}
        | {"gene_name": ["G1", "G2"], "split": ["test", "train"]}
    )
    assert gene_split.cross_table_conflicts(first, second) == {
        "shared_variants": 1,
        "positions_in_multiple_splits": 1,
        "genes_in_multiple_splits": 2,
    }
