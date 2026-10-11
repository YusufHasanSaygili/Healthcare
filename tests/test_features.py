import polars as pl
import pytest

from healthcare.features import missense
from healthcare.features.aminoacids import BLOSUM62, grantham, substitution_features


def test_amino_acid_tables():
    assert grantham("L", "I") == pytest.approx(5, abs=1)
    assert grantham("C", "W") == pytest.approx(215, abs=1)
    assert grantham("U", "A") is None
    assert BLOSUM62[("W", "W")] == 11
    assert all(BLOSUM62[(a, b)] == BLOSUM62[(b, a)] for a, b in BLOSUM62)
    features = substitution_features("D", "K")
    assert features["f_delta_charge"] == 2
    assert substitution_features("U", "C")["f_delta_charge"] is None


def test_resolve_transcripts_drops_multi_gene_variants():
    rows = pl.DataFrame(
        {
            "chrom": ["1", "1", "1"],
            "pos": [10, 20, 20],
            "ref": ["A"] * 3,
            "alt": ["G"] * 3,
            "transcript_id": ["T1", "T1", "T2"],
        }
    )
    resolved, dropped = missense.resolve_transcripts(rows)
    assert resolved["pos"].to_list() == [10]
    assert dropped == 1


def joined_row(**overrides):
    row = {
        "chrom": "1",
        "pos": 10,
        "ref": "A",
        "alt": "G",
        "gene_name": "G",
        "transcript_id": "T",
        "protein_position": 50,
        "aa_ref": "D",
        "aa_alt": "K",
        "label": 1,
        "gnomad_af_raw": None,
        "gnomad_ac_raw": None,
        "gnomad_nhomalt_raw": None,
        "gnomad_af_grpmax_raw": None,
        "gnomad_filter_pass_raw": None,
        "f_gnomad_exome_an": 1000,
        "f_phylop": 5.0,
        "f_phastcons": 1.0,
        "f_protein_length": 100,
        "uniprot_mapped": True,
        **{name: None for name in missense.UNIPROT_GROUPS},
        **{dst: 0.5 for dst in missense.CONSTRAINT_FEATURES.values()},
    }
    return row | overrides


def test_finalize_separates_not_observed_from_not_measured():
    table = missense.finalize(
        pl.DataFrame(
            [
                joined_row(),
                joined_row(pos=11, f_gnomad_exome_an=None),
                joined_row(
                    pos=12, gnomad_af_raw=0.01, gnomad_ac_raw=5, gnomad_filter_pass_raw=True
                ),
                joined_row(pos=13, uniprot_mapped=False, f_uniprot_domain=None),
            ]
        )
    )
    by_pos = {row["pos"]: row for row in table.iter_rows(named=True)}
    assert by_pos[10]["f_gnomad_af"] == 0.0 and by_pos[10]["f_gnomad_covered"] == 1.0
    assert by_pos[11]["f_gnomad_af"] is None and by_pos[11]["f_gnomad_covered"] == 0.0
    assert by_pos[12]["f_gnomad_af"] == 0.01
    assert by_pos[10]["f_uniprot_domain"] == 0.0
    assert by_pos[13]["f_uniprot_domain"] is None
    assert by_pos[10]["f_relative_position"] == 0.5


def test_finalize_exposes_only_whitelisted_columns():
    table = missense.finalize(pl.DataFrame([joined_row()]))
    allowed = {*missense.ID_COLUMNS, "label"}
    assert all(c in allowed or c.startswith("f_") for c in table.columns)
    assert not any(c.startswith(("clinvar_", "validation_")) for c in table.columns)
    assert set(missense.FREQUENCY_FEATURES) <= set(table.columns)
