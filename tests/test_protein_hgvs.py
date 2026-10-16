import pytest

from healthcare.data.protein_hgvs import equivalent, protein_effect

PROTEIN = "MAEERERELKKLQ"


@pytest.mark.parametrize(
    ("ours", "theirs"),
    [
        ("p.Glu3dup", "p.Glu3_Glu4insGlu"),
        ("p.Arg5_Glu6dup", "p.4_5ER[3]"),
        ("p.Arg5_Glu6del", "p.4_5ER[1]"),
        ("p.Lys10fs", "p.Lys10delinsThrTer"),
        ("p.Lys11Ter", "p.Lys10_Lys11insTer"),
        ("p.Glu3_Glu4delinsGln", "p.Glu3_Glu4delinsGln"),
        ("p.Ter14?", "p.Ter14TyrextTer?"),
    ],
)
def test_equivalent_descriptions_of_the_same_change(ours, theirs):
    assert equivalent(ours, theirs, PROTEIN)


@pytest.mark.parametrize(
    ("ours", "theirs"),
    [
        ("p.Glu3dup", "p.Glu3del"),
        ("p.Lys10fs", "p.Lys11fs"),
        ("p.Lys10fs", "p.Lys10del"),
        ("p.Ala2Thr", "p.Ala2Ser"),
    ],
)
def test_different_changes_are_not_equivalent(ours, theirs):
    assert not equivalent(ours, theirs, PROTEIN)


def test_descriptions_that_contradict_the_reference_are_rejected():
    assert protein_effect("p.Trp3del", PROTEIN) is None
    assert protein_effect("p.Ter9?", PROTEIN) is None
    assert protein_effect("p.Ter13?", PROTEIN) is None
    assert not equivalent(None, "p.Ter14TyrextTer?", PROTEIN)
