import pytest

from healthcare.annotate.snv import (
    Transcript,
    TranscriptIndex,
    annotate_snv,
    cds_offset,
    genomic_position,
    parse_protein_change,
    parse_snv,
    protein_change_snvs,
)

PLUS = Transcript(
    transcript_id="TPLUS",
    gene_name="PLUS",
    refseq_transcript="NM_1.1",
    chrom="1",
    strand="+",
    exons=((101, 120), (201, 230)),
    cds=((111, 120), (201, 217)),
    cds_sequence="ATGGCCTGGAAACGTGGATTCCAGTAA",
    protein="MAWKRGFQ",
)
MINUS = Transcript(
    transcript_id="TMINUS",
    gene_name="MINUS",
    refseq_transcript="NM_2.1",
    chrom="1",
    strand="-",
    exons=((5, 19), (31, 40)),
    cds=((11, 19), (31, 36)),
    cds_sequence="ATGTTTCCCGGGTGA",
    protein="MFPG",
)


@pytest.mark.parametrize(
    ("pos", "ref", "alt", "consequence", "hgvs"),
    [
        (111, "A", "G", "start_lost", "p.Met1?"),
        (114, "G", "A", "missense", "p.Ala2Thr"),
        (116, "C", "T", "synonymous", "p.Ala2="),
        (119, "G", "A", "stop_gained", "p.Trp3Ter"),
        (201, "A", "C", "missense", "p.Lys4Thr"),
        (215, "T", "C", "stop_lost", "p.Ter9GlnextTer?"),
        (217, "A", "G", "synonymous", "p.Ter9="),
    ],
)
def test_plus_strand_coding(pos, ref, alt, consequence, hgvs):
    result = annotate_snv(PLUS, pos, ref, alt)
    assert (result.consequence, result.hgvs_p) == (consequence, hgvs)


def test_codon_split_across_exons_uses_both_exons():
    result = annotate_snv(PLUS, 201, "A", "C")
    assert (result.cds_position, result.codon_ref, result.codon_alt) == (11, "AAA", "ACA")


@pytest.mark.parametrize(
    ("pos", "ref", "alt", "hgvs", "cds_position"),
    [
        (33, "A", "C", "p.Phe2Val", 4),
        (19, "G", "A", "p.Pro3Ser", 7),
    ],
)
def test_minus_strand_coding(pos, ref, alt, hgvs, cds_position):
    result = annotate_snv(MINUS, pos, ref, alt)
    assert (result.consequence, result.hgvs_p, result.cds_position) == (
        "missense",
        hgvs,
        cds_position,
    )


@pytest.mark.parametrize(
    ("pos", "consequence"),
    [
        (121, "canonical_splice"),
        (122, "canonical_splice"),
        (123, "intronic"),
        (199, "canonical_splice"),
        (200, "canonical_splice"),
        (198, "intronic"),
        (105, "utr"),
        (225, "utr"),
    ],
)
def test_non_coding_positions(pos, consequence):
    assert annotate_snv(PLUS, pos, "A", "G").consequence == consequence


def test_ref_mismatch_is_flagged():
    assert annotate_snv(PLUS, 114, "C", "T").consequence == "ref_mismatch"


def test_rejects_non_snv():
    with pytest.raises(ValueError):
        annotate_snv(PLUS, 114, "GC", "A")


def test_index_finds_overlapping_transcripts_only():
    index = TranscriptIndex([PLUS, MINUS])
    assert [a.gene_name for a in index.annotate("1", 114, "G", "A")] == ["PLUS"]
    assert [a.gene_name for a in index.annotate("1", 33, "A", "C")] == ["MINUS"]
    assert index.annotate("1", 500, "A", "C") == []
    assert index.annotate("2", 114, "G", "A") == []


def test_parse_snv_normalizes_and_validates():
    assert parse_snv(" chr17:43045712 g>a ") == ("17", 43045712, "G", "A")
    assert parse_snv("X-100 C>T") == ("X", 100, "C", "T")
    for bad in ("23:5 A>G", "0:5 A>G", "1:0 A>G", "1:5 A>A", "1:5 AC>G", "MT:5 A>G", "hello"):
        with pytest.raises(ValueError):
            parse_snv(bad)


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("BRCA1 p.Arg1699Trp", ("BRCA1", "R", 1699, "W")),
        ("brca1 R1699W", ("BRCA1", "R", 1699, "W")),
        ("TP53 p.ARG175his", ("TP53", "R", 175, "H")),
    ],
)
def test_parse_protein_change_accepts_one_and_three_letter_codes(text, expected):
    assert parse_protein_change(text) == expected


@pytest.mark.parametrize("text", ["BRCA1", "BRCA1 p.Xyz12Trp", "p.Arg1699Trp", "BRCA1 Arg-5Trp"])
def test_parse_protein_change_rejects_malformed_input(text):
    with pytest.raises(ValueError):
        parse_protein_change(text)


@pytest.mark.parametrize("transcript", [PLUS, MINUS])
def test_genomic_position_inverts_cds_offset_on_both_strands(transcript):
    for offset in range(len(transcript.cds_sequence)):
        assert cds_offset(transcript, genomic_position(transcript, offset)) == offset
    with pytest.raises(ValueError):
        genomic_position(transcript, len(transcript.cds_sequence))


def test_protein_change_maps_back_to_genomic_snvs():
    assert protein_change_snvs(PLUS, "A", 2, "T") == [("1", 114, "G", "A")]
    assert protein_change_snvs(PLUS, "K", 4, "T") == [("1", 201, "A", "C")]
    assert protein_change_snvs(MINUS, "F", 2, "V") == [("1", 33, "A", "C")]
    snvs = protein_change_snvs(MINUS, "P", 3, "S")
    assert snvs
    for _, pos, ref, alt in snvs:
        assert annotate_snv(MINUS, pos, ref, alt).hgvs_p == "p.Pro3Ser"
    assert protein_change_snvs(PLUS, "A", 2, "W") == []


def test_protein_change_validates_reference_residue_and_length():
    with pytest.raises(ValueError, match="is Ala, not Arg"):
        protein_change_snvs(PLUS, "R", 2, "W")
    with pytest.raises(ValueError, match="has 8 amino acids"):
        protein_change_snvs(PLUS, "A", 9, "T")
