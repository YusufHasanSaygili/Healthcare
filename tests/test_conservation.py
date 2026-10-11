import gzip

import polars as pl
import pytest

from healthcare.data import conservation

WIG = [
    "fixedStep chrom=chr1 start=8 step=1\n",
    *[f"{v}\n" for v in (0.8, 0.9, 1.0, 1.1, 1.2)],
    "fixedStep chrom=chr1 start=20 step=1\n",
    *[f"{v}\n" for v in (2.0, 2.1, 2.2)],
]


def test_parse_wigfix_keeps_only_region_positions():
    positions, values = conservation.parse_wigfix(WIG, [(10, 11), (21, 30)])
    assert positions == [10, 11, 21, 22]
    assert values == pytest.approx([1.0, 1.1, 2.1, 2.2])


def test_parse_wigfix_rejects_non_unit_step():
    with pytest.raises(ValueError):
        conservation.parse_wigfix(["fixedStep chrom=chr1 start=1 step=2\n", "1.0\n"], [(1, 5)])


def test_convert_chromosome_joins_tracks_and_lookup(tmp_path):
    paths = {}
    for track, shift in (("phylop", 0.0), ("phastcons", 10.0)):
        path = tmp_path / f"{track}.wigFix.gz"
        lines = WIG[:6] if track == "phastcons" else WIG
        with gzip.open(path, "wt") as handle:
            handle.writelines(
                line if line.startswith("f") else f"{float(line) + shift}\n" for line in lines
            )
        paths[track] = path
    output = tmp_path / "out" / "chr1.parquet"
    result = conservation.convert_chromosome("1", paths, [(10, 11), (21, 22)], output)
    table = pl.read_parquet(output)
    assert result["rows"] == 4
    assert table.columns == ["chrom", "pos", "phylop", "phastcons"]
    assert table.filter(pl.col("pos") == 21)["phastcons"][0] is None
    hit = conservation.lookup_conservation(output.parent, "1", 10)
    assert hit["phylop"] == pytest.approx(1.0)
    assert hit["phastcons"] == pytest.approx(11.0)
    assert conservation.lookup_conservation(output.parent, "1", 15) is None
