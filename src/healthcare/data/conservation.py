import gzip
import time
from bisect import bisect_right
from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import duckdb
import polars as pl

from healthcare.data.gnomad import coding_regions

TRACKS = ("phylop", "phastcons")


def parse_wigfix(lines: Iterable[str], regions: list[tuple[int, int]]) -> tuple[list, list]:
    starts = [start for start, _ in regions]
    ends = [end for _, end in regions]
    positions: list[int] = []
    values: list[float] = []
    region = 0
    pos = 0
    for line in lines:
        if line[0] == "f":
            fields = dict(item.split("=") for item in line.split()[1:])
            if fields.get("step", "1") != "1":
                raise ValueError(f"unsupported wigFix step: {line.strip()}")
            pos = int(fields["start"]) - 1
            region = bisect_right(ends, pos)
            continue
        pos += 1
        while region < len(ends) and pos > ends[region]:
            region += 1
        if region < len(ends) and pos >= starts[region]:
            positions.append(pos)
            values.append(float(line))
    return positions, values


def convert_track(path: Path, regions: list[tuple[int, int]]) -> pl.DataFrame:
    with gzip.open(path, "rt", encoding="ascii") as handle:
        positions, values = parse_wigfix(handle, regions)
    return pl.DataFrame(
        {"pos": positions, "value": values}, schema={"pos": pl.Int64, "value": pl.Float32}
    )


def convert_chromosome(
    chrom: str, track_paths: dict[str, Path], regions: list[tuple[int, int]], output: Path
) -> dict:
    started = time.monotonic()
    tables = [
        convert_track(path, regions).rename({"value": track}) for track, path in track_paths.items()
    ]
    joined = tables[0]
    for table in tables[1:]:
        joined = joined.join(table, on="pos", how="full", coalesce=True)
    joined = (
        joined.sort("pos")
        .with_columns(pl.lit(chrom).alias("chrom"))
        .select("chrom", "pos", *track_paths)
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    joined.write_parquet(output)
    return {"rows": joined.height, "seconds": round(time.monotonic() - started, 1)}


def build(
    track_paths: dict[str, dict[str, Path]],
    transcripts: pl.DataFrame,
    out_dir: Path,
    workers: int,
    force: bool = False,
) -> dict:
    started = time.monotonic()
    regions = coding_regions(transcripts)
    todo = {
        chrom: paths
        for chrom, paths in track_paths.items()
        if all(p.exists() for p in paths.values())
        and (force or not (out_dir / f"chr{chrom}.parquet").exists())
    }
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {
            chrom: pool.submit(
                convert_chromosome,
                chrom,
                paths,
                regions.get(chrom, []),
                out_dir / f"chr{chrom}.parquet",
            )
            for chrom, paths in todo.items()
        }
        converted = {chrom: future.result()["rows"] for chrom, future in futures.items()}
    files = sorted(out_dir.glob("chr*.parquet"))
    stats = (
        duckdb.execute(
            "SELECT count(*), avg((phylop IS NULL)::int), avg((phastcons IS NULL)::int) "
            "FROM read_parquet(?)",
            [[str(f) for f in files]],
        ).fetchone()
        if files
        else (0, None, None)
    )
    region_bases = sum(end - start + 1 for rs in regions.values() for start, end in rs)
    return {
        "version": "UCSC hg38 phyloP100way + phastCons100way",
        "converted_now": converted,
        "chromosomes_done": len(files),
        "rows_total": stats[0],
        "region_bases": region_bases,
        "phylop_missing_fraction": round(stats[1], 5) if stats[1] is not None else None,
        "phastcons_missing_fraction": round(stats[2], 5) if stats[2] is not None else None,
        "bytes_total": sum(f.stat().st_size for f in files),
        "seconds": round(time.monotonic() - started, 1),
    }


def lookup_conservation(out_dir: Path, chrom: str, pos: int) -> dict | None:
    path = out_dir / f"chr{chrom}.parquet"
    if not path.exists():
        return None
    row = duckdb.execute(
        "SELECT phylop, phastcons FROM read_parquet(?) WHERE pos = ?", [str(path), pos]
    ).fetchone()
    return None if row is None else {"phylop": row[0], "phastcons": row[1]}
