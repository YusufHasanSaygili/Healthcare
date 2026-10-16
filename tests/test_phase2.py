import polars as pl

from healthcare.features import phase2
from tests.test_snv import PLUS


def annotation(pos, ref, alt, gene, consequence, aa_ref=None):
    return {
        "chrom": "1",
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "gene_name": gene,
        "transcript_id": f"T{gene}",
        "consequence": consequence,
        "aa_ref": aa_ref,
    }


def test_resolve_drops_variants_scored_in_two_genes_and_keeps_phase2_types():
    table = pl.DataFrame(
        [
            annotation(10, "A", "G", "A1", "synonymous"),
            annotation(20, "A", "G", "A1", "canonical_splice"),
            annotation(20, "A", "G", "B1", "missense"),
            annotation(30, "A", "G", "A1", "frameshift"),
            annotation(30, "A", "G", "B1", "intronic"),
            annotation(40, "A", "G", "A1", "missense"),
        ]
    )
    kept, dropped = phase2.resolve(table)
    assert sorted(kept["pos"].to_list()) == [10, 30]
    assert dropped == 1


def test_variant_features_follow_the_consequence():
    def utr3(_):
        return "TTAG"

    frameshift = phase2.variant_features(
        {"pos": 114, "ref": "GC", "alt": "G", "consequence": "frameshift", "aa_ref": None},
        PLUS,
        utr3,
    )
    assert frameshift["f_indel_length"] == -1.0
    assert frameshift["f_variant_type"] == float(phase2.VARIANT_TYPES.index("frameshift"))
    assert frameshift["f_nmd_escape"] == 1.0
    assert (frameshift["span_start"], frameshift["span_end"]) == (115, 115)

    deletion = phase2.variant_features(
        {"pos": 205, "ref": "TGGA", "alt": "T", "consequence": "inframe_deletion", "aa_ref": "G"},
        PLUS,
        utr3,
    )
    assert (deletion["f_residues_removed"], deletion["f_residues_added"]) == (1.0, 0.0)
    assert "f_nmd_escape" not in deletion

    splice = phase2.variant_features(
        {"pos": 121, "ref": "G", "alt": "A", "consequence": "canonical_splice", "aa_ref": None},
        PLUS,
        utr3,
    )
    assert splice["f_splice_donor"] == 1.0
    assert splice["f_boundary_distance"] == 1.0


def test_missing_by_type_reports_each_consequence_separately():
    table = pl.DataFrame(
        {
            "consequence": ["frameshift", "frameshift", "synonymous"],
            "f_nmd_escape": [1.0, None, None],
        }
    )
    assert phase2.missing_by_type(table) == {
        "frameshift": {"f_nmd_escape": 0.5},
        "synonymous": {"f_nmd_escape": 1.0},
    }


def test_inserted_residues_are_read_from_each_hgvs_form():
    protein = "MAWKRGFQ"
    assert phase2.inserted_residues("p.Phe7dup", protein, 8, 1) == "F"
    assert phase2.inserted_residues("p.Lys4_Arg5insGly", protein, 5, 1) == "G"
    assert phase2.inserted_residues("p.Arg5delinsTrpLeu", protein, 5, 2) == "WL"
    assert phase2.inserted_residues("p.Gly6del", protein, 6, 0) == ""
