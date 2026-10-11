import polars as pl

from healthcare.annotate.snv import TranscriptIndex
from healthcare.data import annotations
from tests.test_snv import MINUS, PLUS

LABELS = pl.DataFrame(
    {
        "chrom": ["1", "1", "1", "1"],
        "pos": [114, 201, 33, 50],
        "ref": ["G", "A", "A", "C"],
        "alt": ["A", "C", "C", "CA"],
        "validation_transcript": ["NM_1.1", "NM_1.1", "NM_2.0", None],
        "validation_hgvs_p": ["p.Ala2Thr", "p.Lys4Asn", "p.Phe2Val", None],
        "validation_mc": [
            "SO:0001583|missense_variant",
            "SO:0001583|missense_variant,SO:0001627|intron_variant",
            "SO:0001819|synonymous_variant",
            None,
        ],
    }
)


def annotated():
    return annotations.annotate_table(TranscriptIndex([PLUS, MINUS]), LABELS)


def test_annotate_table_skips_non_snvs_and_keeps_refseq():
    table = annotated()
    assert table.height == 3
    assert table.filter(pl.col("pos") == 114)["refseq_transcript"][0] == "NM_1.1"


def test_hgvs_agreement_requires_same_transcript_version():
    result = annotations.hgvs_agreement(LABELS, annotated())
    assert result["missense_compared"] == 2
    assert result["missense_agreement"] == 0.5
    assert result["missense_disagreement_examples"][0]["validation_hgvs_p"] == "p.Lys4Asn"


def test_consequence_agreement_counts_each_variant_once():
    result = annotations.consequence_agreement(LABELS, annotated())
    assert result["consequence_compared"] == 3
    assert round(result["consequence_agreement"], 3) == 0.667
