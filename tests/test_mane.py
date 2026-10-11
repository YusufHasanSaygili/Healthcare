import gzip

import polars as pl
from pyfaidx import Fasta

from healthcare.annotate.sequence import reverse_complement
from healthcare.data import mane

PLUS_CDS = "ATGGCC" + "TGGAAA" + "TAG"
MINUS_CDS = "ATGTTTCCC" + "GGGTGA"
GENOME = (
    "NNNNN"
    + PLUS_CDS[:6]
    + "GTAAGTNNNNAG"
    + PLUS_CDS[6:]
    + "NNNNN"
    + reverse_complement(MINUS_CDS[9:])
    + "CTTACNNNAC"
    + reverse_complement(MINUS_CDS[:9])
    + "NNNNN"
)


def gtf_line(feature, start, end, strand, transcript, gene, protein, tag="MANE_Select"):
    attributes = (
        f'gene_id "G{gene}"; transcript_id "{transcript}"; gene_name "{gene}"; '
        f'tag "{tag}"; protein_id "{protein}"; db_xref "RefSeq:NM_{gene}.1";'
    )
    return f"chr1\ttest\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n"


def write_fixtures(tmp_path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">1\n" + GENOME + "\n")
    p1 = GENOME.index(PLUS_CDS[:6]) + 1
    p2 = GENOME.index(PLUS_CDS[6:]) + 1
    m2 = GENOME.index(reverse_complement(MINUS_CDS[9:])) + 1
    m1 = GENOME.index(reverse_complement(MINUS_CDS[:9])) + 1
    lines = [
        gtf_line("transcript", p1, p2 + 8, "+", "T1", "PLUS", "P1"),
        gtf_line("exon", p1, p1 + 5, "+", "T1", "PLUS", "P1"),
        gtf_line("exon", p2, p2 + 8, "+", "T1", "PLUS", "P1"),
        gtf_line("CDS", p1, p1 + 5, "+", "T1", "PLUS", "P1"),
        gtf_line("CDS", p2, p2 + 5, "+", "T1", "PLUS", "P1"),
        gtf_line("stop_codon", p2 + 6, p2 + 8, "+", "T1", "PLUS", "P1"),
        gtf_line("exon", m2, m2 + 5, "-", "T2", "MINUS", "P2"),
        gtf_line("exon", m1, m1 + 8, "-", "T2", "MINUS", "P2"),
        gtf_line("CDS", m2 + 3, m2 + 5, "-", "T2", "MINUS", "P2"),
        gtf_line("stop_codon", m2, m2 + 2, "-", "T2", "MINUS", "P2"),
        gtf_line("CDS", m1, m1 + 8, "-", "T2", "MINUS", "P2"),
        gtf_line("CDS", p1, p1 + 5, "+", "T3", "PLUSCLIN", "P3", tag="MANE_Plus_Clinical"),
    ]
    gtf = tmp_path / "mane.gtf.gz"
    with gzip.open(gtf, "wt") as handle:
        handle.writelines(lines)
    proteins = tmp_path / "protein.faa.gz"
    with gzip.open(proteins, "wt") as handle:
        handle.write(">P1 desc\nMAWK\n>P2\nMFPG\n")
    return fasta, gtf, proteins


def test_parses_only_mane_select(tmp_path):
    _, gtf, _ = write_fixtures(tmp_path)
    with gzip.open(gtf, "rt") as handle:
        transcripts = mane.parse_mane_gtf(handle)
    assert [t["transcript_id"] for t in transcripts] == ["T1", "T2"]
    assert transcripts[0]["refseq_transcript"] == "NM_PLUS.1"


def test_translates_both_strands(tmp_path):
    fasta, gtf, proteins = write_fixtures(tmp_path)
    with gzip.open(gtf, "rt") as handle:
        transcripts = mane.parse_mane_gtf(handle)
    table, mismatches = mane.build_transcript_table(
        Fasta(str(fasta)), transcripts, mane.read_protein_fasta(proteins)
    )
    assert mismatches == []
    assert table["translation_matches"].all()
    minus = table.filter(pl.col("transcript_id") == "T2").row(0, named=True)
    assert minus["cds_length"] == len(MINUS_CDS)
    assert minus["cds_starts"][0] < minus["cds_starts"][1]


def test_translation_notes():
    assert mane.translation_note("MAWK*", "MAWK") == "exact"
    assert mane.translation_note("LAWK*", "MAWK") == "non_atg_start"
    assert mane.translation_note("MA*K*", "MAUK") == "selenocysteine"
    assert mane.translation_note("VA*K*", "MAUK") == "non_atg_start+selenocysteine"
    assert mane.translation_note("MAWK", "MAWK") == "no stop codon"
    assert mane.translation_note("MAW*", "MAWK") == "length differs"
    assert mane.translation_note("MAYK*", "MAWK") == "sequence differs"
    assert mane.translation_note("MAWK*", None) == "missing official protein"
    assert mane.is_accepted("non_atg_start+selenocysteine")
    assert not mane.is_accepted("sequence differs")


def test_merges_adjacent_intervals():
    assert mane.merge_intervals([(10, 12), (1, 5), (6, 8)]) == [(1, 8), (10, 12)]
