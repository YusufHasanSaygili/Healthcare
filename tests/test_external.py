import gzip
import zipfile

import numpy as np
import polars as pl

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
