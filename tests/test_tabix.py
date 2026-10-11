import gzip
import struct

import polars as pl

from healthcare.data import gnomad, tabix


def test_region_bins_cover_all_levels():
    bins = tabix.region_bins(0, 1)
    assert bins == [0, 1, 9, 73, 585, 4681]
    assert 4681 + 1 in tabix.region_bins(16384, 16385)


def test_merge_chunks_joins_nearby_and_keeps_far_apart():
    near = [tabix.Chunk(100 << 16, 200 << 16), tabix.Chunk(150 << 16, 300 << 16)]
    far = tabix.Chunk(10_000_000 << 16, 10_000_100 << 16)
    merged = tabix.merge_chunks([far, *near])
    assert merged == [tabix.Chunk(100 << 16, 300 << 16), far]


def test_inflate_blocks_skips_truncated_tail():
    blocks = gzip.compress(b"a\tb\n") + gzip.compress(b"c\td\n")
    truncated = blocks + gzip.compress(b"lost\n")[:10]
    assert tabix.inflate_blocks(truncated) == b"a\tb\nc\td\n"


def test_chunk_lines_starts_at_virtual_offset_and_drops_partial_line():
    data = gzip.compress(b"#header\nline1\nline2\npartial")
    chunk = tabix.Chunk(len(b"#header\n"), len(b"#header\nline1\n"))
    assert tabix.chunk_lines(data, chunk) == ["line1", "line2"]


def test_read_index_parses_bins_and_linear_index(tmp_path):
    names = b"chr22\x00"
    header = b"TBI\x01" + struct.pack("<8i", 1, 2, 1, 2, 0, ord("#"), 0, len(names)) + names
    body = struct.pack("<i", 1) + struct.pack("<Ii", 4681, 1) + struct.pack("<2Q", 5 << 16, 9 << 16)
    body += struct.pack("<i", 2) + struct.pack("<2Q", 5 << 16, 7 << 16)
    path = tmp_path / "x.tbi"
    path.write_bytes(gzip.compress(header + body))
    index = tabix.read_index(path)["chr22"]
    assert index.bins[4681] == [tabix.Chunk(5 << 16, 9 << 16)]
    assert index.linear == [5 << 16, 7 << 16]
    assert tabix.region_chunks(index, 0, 100) == [tabix.Chunk(5 << 16, 9 << 16)]


def test_coding_regions_pad_and_merge():
    transcripts = pl.DataFrame(
        {"chrom": ["1", "1"], "cds_starts": [[100, 150], [5]], "cds_ends": [[120, 160], [10]]}
    )
    assert gnomad.coding_regions(transcripts)["1"] == [(1, 30), (80, 180)]
    assert gnomad.in_regions([1, 80], [30, 180], 179)
    assert not gnomad.in_regions([1, 80], [30, 180], 50)


def test_build_coverage_keeps_only_region_positions(tmp_path):
    tsv = tmp_path / "an.tsv.bgz"
    with gzip.open(tsv, "wt") as handle:
        handle.write("locus\tAN\toutside\n")
        handle.writelines(f"chr1:{pos}\t{an}\tfalse\n" for pos, an in ((70, 5), (85, 10), (200, 0)))
    transcripts = pl.DataFrame({"chrom": ["1"], "cds_starts": [[100]], "cds_ends": [[150]]})
    summary = gnomad.build_coverage(tsv, transcripts, tmp_path / "coverage.parquet")
    table = pl.read_parquet(tmp_path / "coverage.parquet")
    assert table["pos"].to_list() == [85]
    assert summary["positions"] == 1
