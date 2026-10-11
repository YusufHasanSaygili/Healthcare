import polars as pl

from healthcare.data import gnomad

INFO = "AC=3;AN=1000;AF=0.003;cadd_phred=25.1;nhomalt=1;AF_grpmax=0.01;grpmax=afr;revel_max=0.9"


def test_info_value_extracts_exact_keys():
    assert gnomad.info_value(INFO, "AC") == "3"
    assert gnomad.info_value(INFO, "AF") == "0.003"
    assert gnomad.info_value(INFO, "grpmax") == "afr"
    assert gnomad.info_value("AN=5;AC_XX=7", "AC") is None
    assert gnomad.info_value("AC=.;AN=5", "AC") is None


def test_parse_line_keeps_only_allowed_fields():
    line = f"chr22\t100\t.\tA\tG\t.\tPASS\t{INFO}\n"
    assert gnomad.parse_line(line) == ("22", 100, "A", "G", True, 3, 1000, 0.003, 1, 0.01, "afr")
    assert gnomad.parse_line("##INFO=<ID=AC>\n") is None


def test_written_parquet_has_no_forbidden_columns(tmp_path):
    lines = [f"chr22\t{100 + i}\t.\tA\tG\t.\tAC0\t{INFO}\n" for i in range(5)]
    output = tmp_path / "chr22.parquet"
    assert gnomad.write_batches(gnomad.parse_lines(lines), output) == 5
    table = pl.read_parquet(output)
    assert table.height == 5
    assert not table["filter_pass"].any()
    assert not any(gnomad.FORBIDDEN.search(column) for column in table.columns)


def test_constraint_drops_polyphen_columns_and_refseq_rows(tmp_path):
    path = tmp_path / "constraint.tsv"
    header = [
        "gene",
        "transcript",
        "mane_select",
        "mis_pphen.oe",
        *list(gnomad.CONSTRAINT_COLUMNS)[2:],
    ]
    rows = [
        ["BRCA1", "ENST00000357654", "true", "0.5", *["0.1"] * 7],
        ["BRCA1", "NM_007294.4", "true", "0.5", *["0.1"] * 7],
        ["OTHER", "ENST00000000001", "false", "0.5", *["NA"] * 7],
    ]
    path.write_text("\n".join("\t".join(r) for r in [header, *rows]) + "\n")
    table = gnomad.read_constraint(path)
    assert table["transcript_id"].to_list() == ["ENST00000357654"]
    assert "mis_pphen.oe" not in table.columns


def test_lookup_frequency_reads_parquet(tmp_path):
    lines = [
        f"chr22\t{100 + i}\t.\tA\tG\t.\tPASS\tAC={i};AN=1000;AF={i / 1000}\n" for i in range(3)
    ]
    gnomad.write_batches(gnomad.parse_lines(lines), tmp_path / "chr22.parquet")
    hit = gnomad.lookup_frequency(tmp_path, "22", 102, "A", "G")
    assert hit["ac"] == 2
    assert hit["af"] == 0.002
    assert gnomad.lookup_frequency(tmp_path, "22", 102, "A", "T") is None
    assert gnomad.lookup_frequency(tmp_path, "1", 102, "A", "G") is None
