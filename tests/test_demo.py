import subprocess
import sys

import polars as pl
import pytest

from healthcare.serve import demo


def transcripts(rows: list[tuple[str, str, str, int, int]]) -> pl.DataFrame:
    return pl.DataFrame(
        {
            "gene_name": [r[0] for r in rows],
            "transcript_id": [r[1] for r in rows],
            "chrom": [r[2] for r in rows],
            "exon_starts": [[r[3], r[3] + 100] for r in rows],
            "exon_ends": [[r[3] + 50, r[4]] for r in rows],
        }
    )


def test_demo_module_never_imports_torch():
    code = "import sys, healthcare.serve.demo; print('torch' in sys.modules)"
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "False"


def test_pack_pulls_in_every_gene_overlapping_the_chosen_regions():
    table = transcripts(
        [
            ("A", "T1", "1", 10_000, 20_000),
            ("B", "T2", "1", 21_000, 30_000),
            ("C", "T3", "1", 31_500, 40_000),
            ("D", "T4", "1", 100_000, 110_000),
            ("E", "T5", "2", 10_000, 20_000),
        ]
    )
    packed = demo.pack_transcripts(table, ["A"])
    assert packed["gene_name"].to_list() == ["A", "B", "C"]
    assert demo.merge_regions(packed) == [("1", 10_000 - demo.PADDING, 40_000 + demo.PADDING)]


def test_pack_rejects_unknown_genes():
    with pytest.raises(ValueError, match="not MANE Select genes: NOPE"):
        demo.pack_transcripts(transcripts([("A", "T1", "1", 10_000, 20_000)]), ["A", "NOPE"])


def test_merge_keeps_separate_chromosomes_and_distant_regions_apart():
    regions = pl.DataFrame({"chrom": ["2", "1", "1"], "start": [5, 300, 1], "end": [9, 400, 100]})
    assert demo.merge_regions(regions) == [("1", 1, 100), ("1", 300, 400), ("2", 5, 9)]


def test_region_genome_serves_slices_like_pyfaidx_and_refuses_outside(tmp_path):
    path = tmp_path / "genome.parquet"
    pl.DataFrame({"chrom": ["17"], "start": [100], "sequence": ["ACGTACGTAC"]}).write_parquet(path)
    genome = demo.RegionGenome(path, {"17": 5_000, "1": 9_000}, ["BRCA1"])
    sequence = genome["17"]
    assert len(sequence) == 5_000
    assert sequence[100:104].seq == "ACGT"
    assert sequence[106:110].seq == "GTAC"
    with pytest.raises(demo.OutsideDemoPack, match="17:110 is outside the demo pack"):
        sequence[109:111]
    with pytest.raises(demo.OutsideDemoPack, match="covers BRCA1"):
        genome["1"][0:1]


def test_outside_pack_errors_are_value_errors_so_serving_reports_them():
    assert issubclass(demo.OutsideDemoPack, ValueError)


def test_load_names_missing_pack_files(tmp_path):
    with pytest.raises(FileNotFoundError, match="manifest.json"):
        demo.load(tmp_path)
    with pytest.raises(FileNotFoundError, match="unzip it in the repository root"):
        demo.load(tmp_path / "demo")


def test_same_treats_missing_values_as_equal():
    assert demo.same(float("nan"), float("nan"))
    assert demo.same(None, None)
    assert not demo.same(None, 0.5)
    assert not demo.same(float("nan"), 0.5)
