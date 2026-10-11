import gzip

import polars as pl
from pyfaidx import Fasta

from healthcare.data import clinvar

GENOME = "ACGTACGTAC" * 10
SINGLE = "criteria_provided,_single_submitter"
EXPERT = "reviewed_by_expert_panel"
NONE = "no_assertion_criteria_provided"


def vcf_line(pos, ref, alt, vid, clnsig, revstat=SINGLE, chrom="1", vc="single_nucleotide_variant"):
    info = (
        f"ALLELEID={vid + 1000};CLNREVSTAT={revstat};CLNSIG={clnsig};CLNVC={vc};"
        f"GENEINFO=GENE{vid}:{vid};MC=SO:0001583|missense_variant"
    )
    return f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t{info}\n"


def write_vcf(path, lines):
    with gzip.open(path, "wt") as handle:
        handle.write("##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        handle.writelines(lines)
    return path


def genome(tmp_path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">1\n" + GENOME + "\n")
    return Fasta(str(fasta))


def test_category_rules():
    assert clinvar.category_of("Pathogenic") == "pathogenic"
    assert clinvar.category_of("Pathogenic/Likely_pathogenic") == "pathogenic"
    assert clinvar.category_of("Likely_benign") == "benign"
    assert clinvar.category_of("Uncertain_significance") == "uncertain"
    assert clinvar.category_of("Conflicting_classifications_of_pathogenicity") == "other"
    assert clinvar.category_of("Pathogenic|risk_factor") == "other"
    assert clinvar.category_of("Pathogenic,_low_penetrance") == "other"
    assert clinvar.category_of("risk_factor") == "other"
    assert clinvar.category_of(None) == "other"


def test_clean_release_applies_filters(tmp_path):
    lines = [
        vcf_line(1, "A", "G", 1, "Pathogenic"),
        vcf_line(2, "C", "T", 2, "Benign", revstat=EXPERT),
        vcf_line(3, "G", "A", 3, "Likely_pathogenic", revstat=NONE),
        vcf_line(4, "T", "C", 4, "Uncertain_significance"),
        vcf_line(5, "G", "T", 5, "Pathogenic"),
        vcf_line(6, "C", "N", 6, "Pathogenic"),
        vcf_line(7, "G", "A", 7, "Pathogenic", chrom="MT"),
        vcf_line(8, "T", "A", 8, "Pathogenic"),
        vcf_line(8, "T", "A", 9, "Benign"),
        vcf_line(10, "C", "CA", 10, "Pathogenic", vc="Insertion"),
    ]
    table, stats = clinvar.clean_release(write_vcf(tmp_path / "a.vcf.gz", lines), genome(tmp_path))
    assert stats["ref_mismatches"] == 1
    assert stats["conflicting_duplicate_keys"] == 1
    labels = clinvar.labeled(table)
    assert dict(zip(labels["pos"], labels["label"], strict=True)) == {1: 1, 2: 0, 10: 1}
    assert labels.filter(pl.col("pos") == 2)["clinvar_stars"][0] == 3
    assert clinvar.uncertain(table)["pos"].to_list() == [4]


def test_temporal_and_reclassified(tmp_path):
    train_lines = [
        vcf_line(1, "A", "G", 1, "Pathogenic"),
        vcf_line(2, "C", "T", 2, "Benign"),
        vcf_line(4, "T", "C", 4, "Uncertain_significance"),
    ]
    test_lines = [
        vcf_line(1, "A", "G", 1, "Pathogenic"),
        vcf_line(2, "C", "T", 2, "Pathogenic"),
        vcf_line(4, "T", "C", 4, "Likely_benign"),
        vcf_line(9, "A", "T", 11, "Likely_pathogenic"),
    ]
    g = genome(tmp_path)
    train, _ = clinvar.clean_release(write_vcf(tmp_path / "train.vcf.gz", train_lines), g)
    test, _ = clinvar.clean_release(write_vcf(tmp_path / "test.vcf.gz", test_lines), g)
    temporal, reclassified = clinvar.temporal_and_reclassified(train, test)
    assert temporal["pos"].to_list() == [9]
    changes = dict(zip(reclassified["pos"], reclassified["change"], strict=True))
    assert changes == {2: "label_flipped", 4: "vus_resolved"}
    assert temporal.join(train, on=clinvar.KEY, how="inner").is_empty()


def test_metadata_columns_are_prefixed(tmp_path):
    table, _ = clinvar.clean_release(
        write_vcf(tmp_path / "a.vcf.gz", [vcf_line(1, "A", "G", 1, "Pathogenic")]), genome(tmp_path)
    )
    unprefixed = set(clinvar.labeled(table).columns) - {*clinvar.KEY, "label"}
    assert all(c.startswith(("clinvar_", "validation_")) for c in unprefixed)


def test_read_hgvs_keeps_grch38_protein_change(tmp_path):
    path = tmp_path / "summary.txt.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("#AlleleID\tName\tAssembly\tVariationID\n")
        handle.write("1\tNM_007294.4(BRCA1):c.5096G>A (p.Arg1699Gln)\tGRCh38\t55\n")
        handle.write("1\tNM_007294.4(BRCA1):c.5096G>A (p.Arg1699Gln)\tGRCh37\t55\n")
        handle.write("2\tNM_000059.4(BRCA2):c.68-7T>A\tGRCh38\t56\n")
    hgvs = clinvar.read_hgvs(path).sort("clinvar_variation_id")
    assert "clinvar_variation_id" in hgvs.columns
    assert hgvs["validation_hgvs_p"].to_list() == ["p.Arg1699Gln", None]
    assert hgvs["validation_transcript"].to_list() == ["NM_007294.4", "NM_000059.4"]
