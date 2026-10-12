import hashlib
import json
import random
from pathlib import Path

import polars as pl

SPLITS = ("train", "validation", "test")
OPEN_SPLITS = ("train", "validation")
TARGETS = {"train": 0.70, "validation": 0.15, "test": 0.15}
SEED = 20261012
KEY = ["chrom", "pos", "ref", "alt"]


class LockedSplitError(RuntimeError):
    pass


class UnionFind:
    def __init__(self, items):
        self.parent = {item: item for item in items}

    def find(self, item):
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, a, b):
        self.parent[self.find(b)] = self.find(a)


def gene_groups(transcripts: pl.DataFrame, variants: pl.DataFrame | None = None) -> dict[str, int]:
    genes = sorted(transcripts["gene_name"].unique())
    union = UnionFind(genes)
    intervals = transcripts.select("chrom", "gene_name", "cds_starts", "cds_ends").explode(
        "cds_starts", "cds_ends"
    )
    for _, chrom_intervals in intervals.group_by("chrom"):
        reach_end, reach_gene = -1, None
        for gene, start, end in (
            chrom_intervals.sort("cds_starts")
            .select("gene_name", "cds_starts", "cds_ends")
            .iter_rows()
        ):
            if start <= reach_end:
                union.union(reach_gene, gene)
            if end > reach_end:
                reach_end, reach_gene = end, gene
    if variants is not None:
        shared = variants.group_by("chrom", "pos").agg(pl.col("gene_name").unique())
        for names in shared["gene_name"]:
            for name in names[1:]:
                union.union(names[0], name)
    roots = sorted({union.find(gene) for gene in genes})
    index = {root: i for i, root in enumerate(roots)}
    return {gene: index[union.find(gene)] for gene in genes}


def assign_groups(group_stats: pl.DataFrame, seed: int = SEED) -> dict[int, str]:
    rows = group_stats.sort("group").to_dicts()
    rng = random.Random(seed)
    rng.shuffle(rows)
    labeled = sorted((r for r in rows if r["variants"]), key=lambda r: -r["variants"])
    unlabeled = [r for r in rows if not r["variants"]]
    totals = {split: [0, 0] for split in SPLITS}
    all_variants = max(1, sum(r["variants"] for r in labeled))
    all_pathogenic = max(1, sum(r["pathogenic"] for r in labeled))
    assignment = {}
    for row in labeled:

        def cost(split: str, row: dict = row) -> float:
            score = 0.0
            for name in SPLITS:
                variants, pathogenic = totals[name]
                if name == split:
                    variants += row["variants"]
                    pathogenic += row["pathogenic"]
                score += (variants / all_variants - TARGETS[name]) ** 2
                score += (pathogenic / all_pathogenic - TARGETS[name]) ** 2
            return score

        best = min(SPLITS, key=cost)
        assignment[row["group"]] = best
        totals[best][0] += row["variants"]
        totals[best][1] += row["pathogenic"]
    weights = [TARGETS[name] for name in SPLITS]
    for row in unlabeled:
        assignment[row["group"]] = rng.choices(SPLITS, weights=weights)[0]
    return assignment


def build_split(transcripts: pl.DataFrame, labeled: pl.DataFrame, seed: int = SEED) -> pl.DataFrame:
    groups = gene_groups(transcripts, labeled.select("chrom", "pos", "gene_name"))
    gene_table = pl.DataFrame({"gene_name": list(groups), "group": list(groups.values())})
    counts = labeled.group_by("gene_name").agg(
        pl.len().alias("variants"), pl.col("label").sum().alias("pathogenic")
    )
    stats = (
        gene_table.join(counts, on="gene_name", how="left")
        .fill_null(0)
        .group_by("group")
        .agg(pl.col("variants").sum(), pl.col("pathogenic").sum())
    )
    assignment = assign_groups(stats, seed)
    return gene_table.with_columns(
        pl.col("group").replace_strict(assignment, return_dtype=pl.String).alias("split")
    ).sort("gene_name")


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def split_summary(table: pl.DataFrame) -> dict:
    stats = (
        table.group_by("split")
        .agg(
            pl.col("gene_name").n_unique().alias("genes"),
            pl.len().alias("variants"),
            pl.col("label").sum().alias("pathogenic"),
        )
        .sort("split")
    )
    total = table.height
    return {
        row["split"]: {
            "genes": row["genes"],
            "variants": row["variants"],
            "variant_share": round(row["variants"] / total, 4),
            "pathogenic": row["pathogenic"],
            "pathogenic_fraction": round(row["pathogenic"] / row["variants"], 4),
        }
        for row in stats.to_dicts()
    }


def overlaps(table: pl.DataFrame) -> dict:
    def count(columns: list[str]) -> int:
        grouped = table.group_by(columns).agg(pl.col("split").n_unique().alias("n"))
        return grouped.filter(pl.col("n") > 1).height

    return {
        "genes_in_multiple_splits": count(["gene_name"]),
        "positions_in_multiple_splits": count(["chrom", "pos"]),
        "variants_in_multiple_splits": count(KEY),
    }


def cds_overlap_violations(transcripts: pl.DataFrame, genes: pl.DataFrame) -> int:
    groups = gene_groups(transcripts)
    split = dict(genes.select("gene_name", "split").iter_rows())
    by_group: dict[int, set] = {}
    for gene, group in groups.items():
        by_group.setdefault(group, set()).add(split[gene])
    return sum(1 for splits in by_group.values() if len(splits) > 1)


def label_temporal(
    temporal: pl.DataFrame, genes: pl.DataFrame, train_positions: pl.DataFrame
) -> pl.DataFrame:
    return (
        temporal.join(genes.select("gene_name", "split"), on="gene_name", how="left")
        .with_columns(
            pl.when(pl.col("split").is_null())
            .then(pl.lit("new_gene"))
            .otherwise(pl.col("split") + "_gene")
            .alias("temporal_subset")
        )
        .drop("split")
        .join(
            train_positions.unique().with_columns(pl.lit(True).alias("shares_position_with_train")),
            on=["chrom", "pos"],
            how="left",
        )
        .with_columns(pl.col("shares_position_with_train").fill_null(False))
    )


def load_split(splits_dir: Path, name: str, allow_locked: bool = False) -> pl.DataFrame:
    if name in OPEN_SPLITS:
        return pl.read_parquet(splits_dir / "missense_open.parquet").filter(pl.col("split") == name)
    if name not in ("test", "temporal"):
        raise ValueError(f"unknown split: {name}")
    if not allow_locked:
        raise LockedSplitError(f"'{name}' is locked; it is opened only in final evaluation tasks")
    return pl.read_parquet(splits_dir / "locked" / f"missense_{name}.parquet")


def identity_hash(labeled: pl.DataFrame) -> str:
    identity = labeled.select(*KEY, "gene_name", "label").sort(KEY)
    return hashlib.sha256(identity.write_csv().encode()).hexdigest()


def check_pin(genes_path: Path, report_path: Path, train_hash: str, force: bool) -> bool:
    report = json.loads(report_path.read_text()) if report_path.exists() else None
    if force or (report is None and not genes_path.exists()):
        return True
    if not genes_path.exists():
        raise RuntimeError(f"{genes_path} is missing but {report_path} pins it; use --force")
    if report is None or report.get("genes_file_sha256") != file_hash(genes_path):
        raise RuntimeError(f"{genes_path} does not match the hash in {report_path}")
    if report.get("train_identity_sha256") != train_hash:
        raise RuntimeError(
            "training variants, genes or labels changed since the split was pinned; use --force"
        )
    return False


def build(
    features_dir: Path, transcripts_path: Path, splits_dir: Path, report_path: Path, force: bool
) -> dict:
    train_path = features_dir / "missense_train_labels.parquet"
    labeled = pl.read_parquet(train_path)
    train_hash = identity_hash(labeled)
    genes_path = splits_dir / "genes.parquet"
    regenerate = check_pin(genes_path, report_path, train_hash, force)
    transcripts = pl.read_parquet(transcripts_path)
    if regenerate:
        splits_dir.mkdir(parents=True, exist_ok=True)
        build_split(transcripts, labeled).write_parquet(genes_path)
    genes = pl.read_parquet(genes_path)

    labeled = labeled.join(genes.select("gene_name", "split"), on="gene_name", how="left")
    if labeled["split"].null_count():
        raise RuntimeError("some training genes have no split assignment")
    locked = splits_dir / "locked"
    locked.mkdir(parents=True, exist_ok=True)
    labeled.filter(pl.col("split").is_in(OPEN_SPLITS)).write_parquet(
        splits_dir / "missense_open.parquet"
    )
    labeled.filter(pl.col("split") == "test").write_parquet(locked / "missense_test.parquet")
    temporal = label_temporal(
        pl.read_parquet(features_dir / "missense_temporal_labels.parquet"),
        genes,
        labeled.filter(pl.col("split") == "train").select("chrom", "pos"),
    )
    temporal.write_parquet(locked / "missense_temporal.parquet")

    temporal_stats = (
        temporal.group_by("temporal_subset")
        .agg(
            pl.len().alias("variants"),
            pl.col("label").sum().alias("pathogenic"),
            pl.col("shares_position_with_train").sum().alias("shares_position_with_train"),
        )
        .sort("temporal_subset")
        .to_dicts()
    )
    summary = {
        "seed": SEED,
        "targets": TARGETS,
        "regenerated": regenerate,
        "genes_file_sha256": file_hash(genes_path),
        "train_identity_sha256": train_hash,
        "mane_genes_assigned": genes.height,
        "gene_groups": genes["group"].n_unique(),
        "merged_gene_groups": genes.group_by("group").len().filter(pl.col("len") > 1).height,
        "cds_overlap_groups_split_apart": cds_overlap_violations(transcripts, genes),
        "genes_per_split_all_mane": dict(genes["split"].value_counts().sort("split").iter_rows()),
        "splits": split_summary(labeled),
        "overlaps": overlaps(labeled),
        "temporal_variants_in_train_release": temporal.join(labeled, on=KEY, how="semi").height,
        "temporal": {row.pop("temporal_subset"): row for row in temporal_stats},
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(summary, indent=2))
    return summary
