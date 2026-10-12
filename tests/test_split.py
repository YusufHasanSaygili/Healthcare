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
