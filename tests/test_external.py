import gzip
import zipfile

import numpy as np
import polars as pl
import pytest

from healthcare.evaluate import external

KEYS = pl.DataFrame(
    {
        "chrom": ["17", "17", "X"],
        "pos": [100, 200, 300],
        "ref": ["G", "C", "A"],
        "alt": ["A", "T", "G"],
        "transcript_id": ["ENST1.4", "ENST2.1", "ENST3.2"],
    }
)


def test_alphamissense_prefers_the_mane_transcript_and_strips_chr(tmp_path):
    path = tmp_path / "am.tsv.gz"
    rows = [
        "# Copyright",
        "#",
        "# Licensed under CC BY-NC-SA 4.0 license",
        "#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class",
        "chr17\t100\tG\tA\thg38\tP1\tENST9.1\tR1W\t0.95\tlikely_pathogenic",
        "chr17\t100\tG\tA\thg38\tP1\tENST1.3\tR1W\t0.40\tambiguous",
        "chr17\t200\tC\tT\thg38\tP2\tENST7.1\tA2V\t0.10\tlikely_benign",
        "chr1\t300\tA\tG\thg38\tP3\tENST3.2\tK3E\t0.80\tlikely_pathogenic",
    ]
    with gzip.open(path, "wt") as handle:
        handle.write("\n".join(rows) + "\n")
    scores = dict(
        external.alphamissense_scores(path, KEYS)
        .select(pl.concat_str("chrom", "pos", separator=":"), "alphamissense")
        .iter_rows()
    )
    assert scores == {"17:100": 0.40, "17:200": 0.10}


def test_revel_is_extracted_once_and_matched_on_grch38_positions(tmp_path):
    archive = tmp_path / "revel.zip"
    rows = [
        "chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid",
        "17,90,100,G,A,R,W,0.900,ENST1;ENST5",
        "17,100,999,C,T,A,V,0.300,ENST2",
        "X,250,300,A,G,K,E,0.200,ENST3",
        "X,250,.,A,G,K,E,0.700,ENST8",
    ]
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr(external.REVEL_MEMBER, "\n".join(rows) + "\n")
    csv = external.revel_csv(archive, tmp_path / "out")
    assert external.revel_csv(archive, tmp_path / "out") == csv
    scores = dict(
        external.revel_scores(csv, KEYS)
        .select(pl.concat_str("chrom", "pos", separator=":"), "revel")
        .iter_rows()
    )
    assert scores == {"17:100": 0.9, "X:300": 0.2}


def test_alphamissense_classes_follow_the_published_cutoffs():
    scores = np.array([0.1, 0.34, 0.5, 0.564, 0.6])
    assert external.alphamissense_classes(scores).tolist() == [0, 1, 1, 1, 2]


def test_three_way_reports_precision_of_each_confident_call():
    labels = np.array([1, 1, 0, 0, 1, 0])
    codes = np.array([2, 2, 2, 0, 1, 1])
    result = external.three_way(labels, codes)
    assert result["decided_share"] == round(4 / 6, 4)
    assert result["likely_pathogenic_precision"] == round(2 / 3, 4)
    assert result["likely_benign_precision"] == 1.0
    assert result["pathogenic_found"] == round(2 / 3, 4)
    assert result["benign_found"] == round(1 / 3, 4)


def test_dbnsfp_scan_keeps_only_wanted_variants_and_stops_at_the_next_chromosome(tmp_path):
    from healthcare.evaluate import dbnsfp

    path = tmp_path / "db.gz"
    rows = [
        "#chr\tpos(1-based)\tref\talt\tREVEL_rankscore",
        "17\t100\tG\tA\t0.9",
        "17\t100\tG\tC\t0.1",
        "17\t200\tC\tT\t.",
        "18\t100\tG\tA\t0.5",
    ]
    with gzip.open(path, "wb") as handle:
        handle.write(("\n".join(rows) + "\n").encode())
    start = len(rows[0]) + 1
    found = dbnsfp.scan(path, "17", start, {(100, "G", "A"), (200, "C", "T")}, [4])
    assert found == [["17", 100, "G", "A", "0.9"], ["17", 200, "C", "T", "."]]


def test_broad_scores_ours_on_exactly_the_variants_each_tool_covers():
    from healthcare.evaluate import dbnsfp

    rng = np.random.default_rng(1)
    n = 400
    labels = rng.integers(0, 2, n)
    table = pl.DataFrame(
        {
            "label": labels,
            "gene_name": [f"G{i % 40}" for i in range(n)],
            "v2": labels * 0.8 + rng.random(n) * 0.2,
            "v2_no_frequency": rng.random(n),
            **{column: None for column in dbnsfp.TOOLS},
        }
    ).with_columns(
        pl.Series("REVEL_rankscore", np.where(np.arange(n) < 300, rng.random(n), np.nan)).fill_nan(
            None
        )
    )
    rows = external.broad(table, rounds=200)
    assert [row["tool"] for row in rows] == ["REVEL"]
    assert rows[0]["variants"] == 300
    low, high = rows[0]["difference_95ci"]
    assert low <= rows[0]["difference"] <= high
    assert rows[0]["difference"] > 0.3
    assert rows[0]["verdict"] == "ahead"


def test_holm_adjusts_in_rank_order_and_stays_monotone():
    assert external.holm([0.01, 0.04, 0.03, 0.005]) == pytest.approx([0.03, 0.06, 0.06, 0.02])


def test_verdict_needs_a_holm_adjusted_p_below_alpha():
    assert external.verdict(0.01, 0.01) == "ahead"
    assert external.verdict(-0.01, 0.01) == "behind"
    assert external.verdict(-0.01, 0.2) == "tie"


def test_bootstrap_p_value_is_two_sided_and_never_zero():
    assert external.bootstrap_p_value(np.full(999, 0.02)) == pytest.approx(2 / 1000)
    assert external.bootstrap_p_value(np.array([-1.0, 1.0])) == 1.0
