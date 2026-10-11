import gzip

import polars as pl

from healthcare.data import uniprot

ENTRY = """ID   TEST_HUMAN              Reviewed;          10 AA.
AC   P99999; Q00000;
GN   Name=TESTG {ECO:0000303};
FT   DOMAIN          2..6
FT                   /note="RING-type"
FT                   /evidence="ECO:0000255|PROSITE-ProRule:PRU00175"
FT   VARIANT         3
FT                   /note="R -> W (in disease; pathogenic)"
FT                   /evidence="ECO:0000269|PubMed:1"
FT   MUTAGEN         4
FT                   /note="C->A: Loss of activity."
FT   REGION          1..10
FT                   /note="Disordered"
FT                   /evidence="ECO:0000256|SAM:MobiDB-lite"
FT   DISULFID        2..9
FT                   /evidence="ECO:0000269|PubMed:2"
FT   BINDING         <1..>3
FT   ACT_SITE        Q11111:5
SQ   SEQUENCE   10 AA;  1000 MW;  0 CRC64;
     MACDEFGHCK
//
"""


def test_parse_entries_keeps_curated_types_only():
    (entry,) = list(uniprot.parse_entries(ENTRY.splitlines(keepends=True)))
    assert entry["accession"] == "P99999"
    assert entry["gene"] == "TESTG"
    assert entry["sequence"] == "MACDEFGHCK"
    types = [f[0] for f in entry["features"]]
    assert types == ["DOMAIN", "DISULFID", "DISULFID", "BINDING"]
    assert entry["features"][0][1:4] == (2, 6, "RING-type")
    assert [f[1] for f in entry["features"] if f[0] == "DISULFID"] == [2, 9]


def test_is_curated_drops_purely_automatic_evidence():
    assert uniprot.is_curated({"ECO:0000269"})
    assert uniprot.is_curated(set())
    assert uniprot.is_curated({"ECO:0000256", "ECO:0000269"})
    assert not uniprot.is_curated({"ECO:0000256"})


def test_build_maps_by_exact_sequence(tmp_path):
    dat = tmp_path / "sprot.dat.gz"
    with gzip.open(dat, "wt") as handle:
        handle.write(ENTRY + ENTRY.replace("MACDEFGHCK", "MACDEFGHCR").replace("P99999", "P88888"))
    transcripts = pl.DataFrame(
        {
            "transcript_id": ["T1", "T2"],
            "gene_name": ["TESTG", "OTHER"],
            "protein": ["MACDEFGHCK", "MXXX"],
        }
    )
    summary = uniprot.build(dat, transcripts, tmp_path / "features.parquet")
    table = pl.read_parquet(tmp_path / "features.parquet")
    assert summary["mane_proteins_matched"] == 1
    assert set(table["uniprot_accession"]) == {"P99999"}
    assert not set(table["feature_type"]) & uniprot.EXCLUDED_TYPES
