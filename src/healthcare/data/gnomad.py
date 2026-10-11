import gzip
import re
import time
from collections.abc import Iterable, Iterator
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import duckdb
import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq

INFO_FIELDS = {"AC": int, "AN": int, "AF": float, "nhomalt": int, "AF_grpmax": float, "grpmax": str}
FORBIDDEN = re.compile(r"cadd|revel|spliceai|pangolin|polyphen|pphen|sift", re.IGNORECASE)
SCHEMA = pa.schema(
    [
        ("chrom", pa.string()),
        ("pos", pa.int64()),
        ("ref", pa.string()),
        ("alt", pa.string()),
        ("filter_pass", pa.bool_()),
        ("ac", pa.int64()),
        ("an", pa.int64()),
        ("af", pa.float64()),
        ("nhomalt", pa.int64()),
        ("af_grpmax", pa.float64()),
        ("grpmax", pa.string()),
    ]
)
BATCH_ROWS = 500_000
VERSION = "gnomAD v4.1 exomes"
assert not any(FORBIDDEN.search(name) for name in SCHEMA.names)
CONSTRAINT_COLUMNS = {
    "gene": "gene_name",
    "transcript": "transcript_id",
    "lof.oe": "lof_oe",
    "lof.oe_ci.upper": "loeuf",
    "lof.pLI": "pli",
    "mis.oe": "mis_oe",
    "mis.oe_ci.upper": "mis_oe_upper",
    "mis.z_score": "mis_z",
    "syn.z_score": "syn_z",
}


def info_value(info: str, key: str) -> str | None:
    marker = key + "="
    if info.startswith(marker):
        start = len(marker)
    else:
        index = info.find(";" + marker)
        if index == -1:
            return None
        start = index + 1 + len(marker)
    end = info.find(";", start)
    value = info[start:] if end == -1 else info[start:end]
    return value if value not in ("", ".") else None


def parse_line(line: str) -> tuple | None:
    if line.startswith("#"):
        return None
    chrom, pos, _, ref, alt, _, filter_, info = line.split("\t", 8)[:8]
    values = []
    for key, cast in INFO_FIELDS.items():
        raw = info_value(info, key)
        values.append(None if raw is None else cast(raw))
    return (chrom.removeprefix("chr"), int(pos), ref, alt, filter_ == "PASS", *values)


def parse_lines(lines: Iterable[str]) -> Iterator[tuple]:
    for line in lines:
        record = parse_line(line)
        if record is not None:
            yield record


def to_table(batch: list[tuple]) -> pa.Table:
    columns = list(zip(*batch, strict=True))
    return pa.Table.from_arrays(
        [pa.array(column, type=field.type) for column, field in zip(columns, SCHEMA, strict=True)],
        schema=SCHEMA,
    )


def write_batches(records: Iterator[tuple], output: Path, limit: int | None = None) -> int:
    output.parent.mkdir(parents=True, exist_ok=True)
    part = output.with_name(output.name + ".part")
    count = 0
    with pq.ParquetWriter(part, SCHEMA, compression="zstd") as writer:
        batch: list[tuple] = []
        for record in records:
            batch.append(record)
            count += 1
            if len(batch) == BATCH_ROWS:
                writer.write_table(to_table(batch))
                batch = []
            if limit and count >= limit:
                break
        if batch:
            writer.write_table(to_table(batch))
    part.replace(output)
    return count


def convert_vcf(vcf: Path, output: Path, limit: int | None = None) -> dict:
    started = time.monotonic()
    with gzip.open(vcf, "rt", encoding="utf-8") as handle:
        rows = write_batches(parse_lines(handle), output, limit)
    return {
        "output": str(output),
        "rows": rows,
        "bytes": output.stat().st_size,
        "seconds": round(time.monotonic() - started, 1),
    }


def read_constraint(path: Path) -> pl.DataFrame:
    table = pl.read_csv(path, separator="\t", infer_schema_length=0, null_values=["NA", ""])
    selected = (
        table.filter(
            (pl.col("mane_select") == "true") & pl.col("transcript").str.starts_with("ENST")
        )
        .select(list(CONSTRAINT_COLUMNS))
        .rename(CONSTRAINT_COLUMNS)
        .with_columns(pl.exclude("gene_name", "transcript_id").cast(pl.Float64))
        .unique("transcript_id", keep="first", maintain_order=True)
    )
    assert not any(FORBIDDEN.search(column) for column in selected.columns)
    return selected


def build(
    vcfs: dict[str, Path],
    constraint_tsv: Path,
    out_dir: Path,
    workers: int,
    limit: int | None = None,
    force: bool = False,
) -> dict:
    started = time.monotonic()
    exome_dir = out_dir / "exomes"
    todo = {
        chrom: vcf
        for chrom, vcf in vcfs.items()
        if vcf.exists() and (force or not (exome_dir / f"chr{chrom}.parquet").exists())
    }
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {
            chrom: pool.submit(convert_vcf, vcf, exome_dir / f"chr{chrom}.parquet", limit)
            for chrom, vcf in todo.items()
        }
        converted = {chrom: future.result() for chrom, future in futures.items()}
    constraint = read_constraint(constraint_tsv)
    constraint.write_parquet(out_dir / "constraint.parquet")
    done = sorted(path.stem.removeprefix("chr") for path in exome_dir.glob("chr*.parquet"))
    return {
        "version": VERSION,
        "converted_now": {chrom: result["rows"] for chrom, result in converted.items()},
        "chromosomes_done": done,
        "chromosomes_missing_raw": sorted(c for c, v in vcfs.items() if not v.exists()),
        "exome_rows_total": sum(
            pq.ParquetFile(path).metadata.num_rows for path in exome_dir.glob("chr*.parquet")
        ),
        "exome_bytes_total": sum(path.stat().st_size for path in exome_dir.glob("chr*.parquet")),
        "constraint_genes": constraint.height,
        "seconds": round(time.monotonic() - started, 1),
    }


def lookup_frequency(exome_dir: Path, chrom: str, pos: int, ref: str, alt: str) -> dict | None:
    path = exome_dir / f"chr{chrom}.parquet"
    if not path.exists():
        return None
    row = duckdb.execute(
        "SELECT filter_pass, ac, an, af, nhomalt, af_grpmax, grpmax FROM read_parquet(?) "
        "WHERE pos = ? AND ref = ? AND alt = ?",
        [str(path), pos, ref, alt],
    ).fetchone()
    if row is None:
        return None
    return dict(
        zip(["filter_pass", "ac", "an", "af", "nhomalt", "af_grpmax", "grpmax"], row, strict=True)
    )
