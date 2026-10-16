from dataclasses import replace

import pytest

from healthcare.annotate.snv import Transcript
from healthcare.annotate.truncation import (
    boundary_distance,
    cds_edit,
    cds_start,
    last_junction,
    splice_features,
    start_loss_features,
    truncation_features,
    utr3_sequence,
)
from tests.test_snv import MINUS, PLUS

LONG = Transcript(
    transcript_id="TLONG",
    gene_name="LONG",
    refseq_transcript=None,
    chrom="1",
    strand="+",
    exons=((1, 100), (201, 300), (401, 500)),
    cds=((51, 100), (201, 300), (401, 448)),
    cds_sequence="ATG" + "GCT" * 64 + "TAA",
    protein="M" + "A" * 64,
)


def test_mrna_coordinates_of_cds_start_and_last_junction():
    assert (cds_start(LONG), last_junction(LONG)) == (50, 200)
    assert (cds_start(PLUS), last_junction(PLUS)) == (10, 20)
    assert (cds_start(MINUS), last_junction(MINUS)) == (4, 10)


def test_utr3_follows_the_stop_codon_on_both_strands():
    def fetch(start, end):
        return "".join("ACGT"[pos % 4] for pos in range(start, end + 1))

    assert len(utr3_sequence(PLUS, fetch)) == 13
    assert utr3_sequence(MINUS, fetch) == "CGTACG"


def test_cds_edit_maps_snvs_and_indels_to_transcript_coordinates():
    assert cds_edit(PLUS, 114, "G", "A") == (3, 3, "A")
    assert cds_edit(PLUS, 114, "GC", "G") == (4, 4, "")
    assert cds_edit(MINUS, 33, "A", "AT") == (3, 2, "A")
    assert cds_edit(PLUS, 150, "NN", "N") is None


def test_early_stop_far_from_the_last_junction_triggers_nmd():
    features = truncation_features(LONG, (27, 29, "TAA"), "")
    assert features["f_nmd_escape"] == 0.0
    assert features["f_ptc_junction_distance"] == 123.0
    assert features["f_protein_retained"] == pytest.approx(9 / 65)
    assert features["f_novel_residues"] == 0.0


def test_stop_in_the_last_exon_escapes_nmd():
    features = truncation_features(LONG, (159, 161, "TAA"), "")
    assert (features["f_ptc_in_last_exon"], features["f_nmd_escape"]) == (1.0, 1.0)


def test_frameshift_reads_into_the_3_prime_utr():
    features = truncation_features(LONG, (190, 190, ""), "CTAGG")
    assert features["f_no_stop"] == 0.0
    assert features["f_novel_residues"] == 3.0
    assert truncation_features(LONG, (190, 190, ""), "CCCCCC")["f_no_stop"] == 1.0


def test_start_loss_finds_the_next_methionine():
    assert start_loss_features(PLUS) == {"f_next_met_fraction": None}
    assert start_loss_features(replace(PLUS, protein="MAWMRGFQ")) == {"f_next_met_fraction": 3 / 8}


@pytest.mark.parametrize(
    ("transcript", "pos", "expected"),
    [
        (PLUS, 121, (1.0, 1.0, 0.0, 0.0)),
        (PLUS, 199, (0.0, 2.0, 0.0, 1.0)),
        (MINUS, 30, (1.0, 1.0, 1.0, 0.0)),
    ],
)
def test_splice_features_name_the_site_and_the_skipped_exon(transcript, pos, expected):
    features = splice_features(transcript, pos, pos)
    assert (
        features["f_splice_donor"],
        features["f_splice_intron_position"],
        features["f_skipped_exon_in_frame"],
        features["f_skipped_exon_rank"],
    ) == expected


def test_splice_features_are_empty_away_from_splice_sites():
    assert splice_features(PLUS, 150, 150) == {}


def test_boundary_distance_is_to_the_nearest_exon_edge():
    assert boundary_distance(PLUS, 118) == 2.0
