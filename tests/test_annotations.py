import polars as pl

from healthcare.annotate.snv import TranscriptIndex
from healthcare.data import annotations
from tests.test_snv import MINUS, PLUS

LABELS = pl.DataFrame(
    {
        "chrom": ["1", "1", "1", "1", "1"],
        "pos": [114, 201, 33, 50, 114],
        "ref": ["G", "A", "A", "C", "GC"],
        "alt": ["A", "C", "C", "CA", "G"],
        "validation_transcript": ["NM_1.1", "NM_1.1", "NM_2.0", None, "NM_1.1"],
        "validation_hgvs_p": ["p.Ala2Thr", "p.Lys4Asn", "p.Phe2Val", None, "p.Trp3fs"],
        "validation_mc": [
            "SO:0001583|missense_variant",
            "SO:0001583|missense_variant,SO:0001627|intron_variant",
            "SO:0001819|synonymous_variant",
            None,
            "SO:0001589|frameshift_variant",
        ],
    }
)


def annotated():
    table, _ = annotations.annotate_table(
        TranscriptIndex([PLUS, MINUS]), LABELS, lambda chrom: lambda start, end: "N"
    )
    return table


def test_annotate_table_covers_snvs_and_indels_and_keeps_refseq():
    table = annotated()
    assert table.height == 4
    assert table.filter(pl.col("ref") == "GC")["consequence"][0] == "frameshift"
    assert table.filter(pl.col("pos") == 114)["refseq_transcript"][0] == "NM_1.1"


def test_annotate_table_counts_indels_moved_by_normalization():
    labels = pl.DataFrame({"chrom": ["1"], "pos": [205], "ref": ["T"], "alt": ["TT"]})
    _, moved = annotations.annotate_table(
        TranscriptIndex([PLUS]), labels, lambda chrom: lambda start, end: "A"
    )
    assert moved == 1


def test_hgvs_by_type_separates_snvs_from_indels():
    report = annotations.hgvs_by_type(
        LABELS, annotated(), {"TPLUS": PLUS.protein, "TMINUS": MINUS.protein}
    )
    assert report["missense/snv"]["compared"] == 2
    assert report["frameshift/indel"] == {
        "compared": 1,
        "agreement": 1.0,
        "equivalent": 1.0,
        "disagreement_examples": [],
    }


def test_hgvs_agreement_requires_same_transcript_version():
    result = annotations.hgvs_agreement(LABELS, annotated())
    assert result["missense_compared"] == 2
    assert result["missense_agreement"] == 0.5
    assert result["missense_disagreement_examples"][0]["validation_hgvs_p"] == "p.Lys4Asn"


def test_consequence_agreement_counts_each_variant_once():
    result = annotations.consequence_agreement(LABELS, annotated())
    assert result["consequence_compared"] == 3
    assert round(result["consequence_agreement"], 3) == 0.667
    assert result["indel_consequence_by_type"] == {"frameshift": {"compared": 1, "agreement": 1.0}}
