from healthcare.annotate.sequence import reverse_complement, translate


def test_reverse_complement():
    assert reverse_complement("ATGCN") == "NGCAT"


def test_translate_standard_code():
    assert translate("ATGGCCTGGTAA") == "MAW*"
    assert translate("ATGGC") == "M"
    assert translate("ATGNNN") == "MX"
