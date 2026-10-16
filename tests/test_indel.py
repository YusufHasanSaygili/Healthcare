import pytest

from healthcare.annotate.indel import annotate_indel, annotate_variant, normalize, parse_variant
from healthcare.annotate.snv import Transcript, TranscriptIndex
from tests.test_snv import MINUS, PLUS


@pytest.mark.parametrize(
    ("pos", "ref", "alt", "consequence", "hgvs"),
    [
        (114, "GC", "G", "frameshift", "p.Trp3fs"),
        (117, "T", "TA", "frameshift", "p.Trp3Ter"),
        (205, "TGGA", "T", "inframe_deletion", "p.Gly6del"),
        (208, "A", "ATTC", "inframe_insertion", "p.Phe7dup"),
        (202, "A", "AGGG", "inframe_insertion", "p.Lys4_Arg5insGly"),
        (203, "CGT", "TGGTTA", "inframe_indel", "p.Arg5delinsTrpLeu"),
        (202, "A", "ATAA", "stop_gained", "p.Lys4_Arg5insTer"),
        (111, "ATGG", "A", "start_lost", "p.Met1?"),
    ],
)
def test_plus_strand_indels(pos, ref, alt, consequence, hgvs):
    result = annotate_indel(PLUS, pos, ref, alt)
    assert (result.consequence, result.hgvs_p) == (consequence, hgvs)


def test_in_frame_deletion_of_the_stop_codon_is_stop_lost():
    assert annotate_indel(PLUS, 214, "GTAA", "G").consequence == "stop_lost"


def test_in_frame_insertion_inside_a_kept_stop_codon_leaves_the_protein_unchanged():
    result = annotate_indel(PLUS, 216, "A", "AAAA")
    assert (result.consequence, result.hgvs_p) == ("synonymous", "p.Ter9=")


def test_frameshift_after_a_kept_stop_codon_leaves_the_protein_unchanged():
    result = annotate_indel(PLUS, 216, "A", "AGC")
    assert (result.consequence, result.hgvs_p) == ("synonymous", "p.Ter9=")


@pytest.mark.parametrize(
    ("pos", "ref", "alt", "consequence", "hgvs"),
    [
        (33, "A", "AT", "frameshift", "p.Phe2fs"),
        (16, "CGGG", "C", "inframe_deletion", "p.Pro3del"),
    ],
)
def test_minus_strand_indels(pos, ref, alt, consequence, hgvs):
    result = annotate_indel(MINUS, pos, ref, alt)
    assert (result.consequence, result.hgvs_p) == (consequence, hgvs)


@pytest.mark.parametrize(
    ("pos", "ref", "alt", "consequence"),
    [
        (120, "ANN", "A", "canonical_splice"),
        (121, "N", "NA", "canonical_splice"),
        (120, "A", "AT", "splice_region"),
        (150, "NN", "N", "intronic"),
        (104, "N", "NA", "utr"),
        (114, "GC", "AT", "multi_nucleotide"),
    ],
)
def test_indels_outside_the_coding_frame(pos, ref, alt, consequence):
    assert annotate_indel(PLUS, pos, ref, alt).consequence == consequence


def test_deletion_that_does_not_match_the_cds_is_flagged():
    assert annotate_indel(PLUS, 114, "GA", "G").consequence == "ref_mismatch"


def test_deletion_across_the_start_codon_boundary_is_start_lost():
    assert annotate_indel(PLUS, 109, "NNAT", "N").consequence == "start_lost"


GENOME = "GGCACACATTTTTGCAG"


def fetch(start: int, end: int) -> str:
    return GENOME[start - 1 : end]


@pytest.mark.parametrize(
    ("pos", "ref", "alt", "expected"),
    [
        (5, "CAC", "C", (2, "GCA", "G")),
        (12, "T", "TT", (8, "A", "AT")),
        (9, "TTTT", "TTT", (8, "AT", "A")),
        (3, "CAC", "CTC", (4, "A", "T")),
        (14, "GCA", "GA", (14, "GC", "G")),
    ],
)
def test_normalize_left_aligns_and_trims(pos, ref, alt, expected):
    assert normalize(pos, ref, alt, fetch) == expected


def test_parse_variant_accepts_indels_and_validates():
    assert parse_variant("chr17:43045712 ga>g") == ("17", 43045712, "GA", "G")
    with pytest.raises(ValueError):
        parse_variant("17:100 GA>GA")
    with pytest.raises(ValueError):
        parse_variant("M:100 A>AT")


def test_annotate_variant_dispatches_snvs_and_indels():
    index = TranscriptIndex([PLUS, MINUS])
    assert [a.hgvs_p for a in annotate_variant(index, "1", 114, "G", "A")] == ["p.Ala2Thr"]
    assert [a.hgvs_p for a in annotate_variant(index, "1", 205, "TGGA", "T")] == ["p.Gly6del"]
    assert [a.consequence for a in annotate_variant(index, "1", 120, "ANN", "A")] == [
        "canonical_splice"
    ]


SELENO = Transcript(
    transcript_id="TSEC",
    gene_name="SEC",
    refseq_transcript=None,
    chrom="1",
    strand="+",
    exons=((1, 15),),
    cds=((1, 15),),
    cds_sequence="ATGTGAGCCTGGTAA",
    protein="MUAW",
)


def test_selenocysteine_codons_are_not_read_as_stops():
    result = annotate_indel(SELENO, 6, "AGCC", "A")
    assert (result.consequence, result.hgvs_p) == ("inframe_deletion", "p.Ala3del")


def test_normalize_rejects_an_allele_that_cannot_be_anchored():
    with pytest.raises(ValueError):
        normalize(1, "A", "AA", fetch)


def test_frameshift_inside_the_stop_codon_is_a_stop_loss():
    result = annotate_indel(PLUS, 215, "TA", "T")
    assert (result.consequence, result.hgvs_p) == ("stop_lost", "p.Ter9?")
