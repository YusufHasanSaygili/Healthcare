import json
import shutil
import time
import zipfile
from pathlib import Path

import duckdb
import numpy as np
import polars as pl
from sklearn.metrics import roc_auc_score

from healthcare.data.download import ChecksumMismatch, digests_of
from healthcare.evaluate import dbnsfp
from healthcare.evaluate.all_types import file_hash, predictions
from healthcare.evaluate.metrics import classification_metrics, gene_bootstrap, rounded
from healthcare.models.calibration import decide
from healthcare.models.joint import load_joint

KEY = ["chrom", "pos", "ref", "alt"]
REVEL_MEMBER = "revel_with_transcript_ids"
ALPHAMISSENSE_BENIGN_BELOW = 0.34
ALPHAMISSENSE_PATHOGENIC_ABOVE = 0.564
BOOTSTRAP_ROUNDS = 2000
RANKING_METRICS = ("roc_auc", "pr_auc", "n", "positives")
SCORERS = {
    "v2": "Ours: v2 (all features)",
    "v2_no_frequency": "Ours: LightGBM without allele frequency",
    "alphamissense": "AlphaMissense (2023)",
    "revel": "REVEL (2016)",
}


def alphamissense_scores(path: Path, keys: pl.DataFrame) -> pl.DataFrame:
    con = duckdb.connect()
    con.register("keys", keys.select(*KEY, "transcript_id").to_arrow())
    return con.execute(
        f"""
        WITH am AS (
            SELECT replace(column0, 'chr', '') AS chrom, column1::BIGINT AS pos,
                   column2 AS ref, column3 AS alt,
                   split_part(column6, '.', 1) AS transcript, column8::DOUBLE AS score
            FROM read_csv('{path.as_posix()}', delim = '\t', header = false, skip = 4,
                          all_varchar = true)
        )
        SELECT k.chrom, k.pos, k.ref, k.alt,
               coalesce(
                   max(am.score) FILTER (WHERE am.transcript = split_part(k.transcript_id, '.', 1)),
                   max(am.score)
               ) AS alphamissense,
               count(*) FILTER (
                   WHERE am.transcript = split_part(k.transcript_id, '.', 1)
               ) > 0 AS alphamissense_on_mane
        FROM keys k JOIN am USING (chrom, pos, ref, alt)
        GROUP BY k.chrom, k.pos, k.ref, k.alt
        """
    ).pl()


def revel_csv(zip_path: Path, out_dir: Path) -> Path:
    target = out_dir / f"{REVEL_MEMBER}.csv"
    if not target.exists():
        out_dir.mkdir(parents=True, exist_ok=True)
        part = target.with_name(target.name + ".part")
        with (
            zipfile.ZipFile(zip_path) as archive,
            archive.open(REVEL_MEMBER) as source,
            part.open("wb") as sink,
        ):
            shutil.copyfileobj(source, sink, 1 << 20)
        part.replace(target)
    return target


def revel_scores(path: Path, keys: pl.DataFrame) -> pl.DataFrame:
    con = duckdb.connect()
    con.register("keys", keys.select(*KEY, "transcript_id").to_arrow())
    return con.execute(
        f"""
        WITH rv AS (
            SELECT chr AS chrom, TRY_CAST(grch38_pos AS BIGINT) AS pos, ref, alt,
                   TRY_CAST(REVEL AS DOUBLE) AS score, Ensembl_transcriptid AS transcripts
            FROM read_csv('{path.as_posix()}', header = true, all_varchar = true)
        )
        SELECT k.chrom, k.pos, k.ref, k.alt,
               coalesce(
                   max(rv.score) FILTER (
                       WHERE contains(rv.transcripts, split_part(k.transcript_id, '.', 1))
                   ),
                   max(rv.score)
               ) AS revel,
               count(*) FILTER (
                   WHERE contains(rv.transcripts, split_part(k.transcript_id, '.', 1))
               ) > 0 AS revel_on_mane
        FROM keys k JOIN rv USING (chrom, pos, ref, alt)
        GROUP BY k.chrom, k.pos, k.ref, k.alt
        """
    ).pl()


def alphamissense_classes(scores: np.ndarray) -> np.ndarray:
    return np.where(
        scores > ALPHAMISSENSE_PATHOGENIC_ABOVE,
        2,
        np.where(scores < ALPHAMISSENSE_BENIGN_BELOW, 0, 1),
    )


def three_way(labels: np.ndarray, codes: np.ndarray) -> dict:
    pathogenic_calls, benign_calls = codes == 2, codes == 0
    return {
        "decided_share": round(float((codes != 1).mean()), 4),
        "uncertain_share": round(float((codes == 1).mean()), 4),
        "likely_pathogenic_precision": round(float(labels[pathogenic_calls].mean()), 4)
        if pathogenic_calls.any()
        else None,
        "likely_benign_precision": round(float(1 - labels[benign_calls].mean()), 4)
        if benign_calls.any()
        else None,
        "pathogenic_found": round(float(pathogenic_calls[labels == 1].mean()), 4),
        "benign_found": round(float(benign_calls[labels == 0].mean()), 4),
    }


def compare(table: pl.DataFrame) -> dict:
    labels = table["label"].to_numpy().astype(int)
    genes = table["gene_name"].to_numpy()
    metrics = {}
    for name in SCORERS:
        scores = table[name].to_numpy().astype(float)
        ranking = classification_metrics(labels, scores)
        metrics[name] = rounded({k: ranking[k] for k in RANKING_METRICS}) | {
            "gene_bootstrap_95ci": gene_bootstrap(labels, scores, genes)
        }
    return {
        "metrics": metrics,
        "three_way": {
            "v2": three_way(labels, decide(table["v2"].to_numpy())),
            "alphamissense": three_way(
                labels, alphamissense_classes(table["alphamissense"].to_numpy())
            ),
        },
    }


def paired_differences(
    labels: np.ndarray, ours: np.ndarray, theirs: np.ndarray, genes: np.ndarray, rounds: int
) -> np.ndarray:
    rng = np.random.default_rng(0)
    unique, inverse = np.unique(genes, return_inverse=True)
    members = [np.flatnonzero(inverse == i) for i in range(unique.size)]
    differences = []
    for _ in range(rounds):
        index = np.concatenate([members[i] for i in rng.integers(0, unique.size, unique.size)])
        y = labels[index]
        if 0 < y.sum() < y.size:
            differences.append(roc_auc_score(y, ours[index]) - roc_auc_score(y, theirs[index]))
    return np.asarray(differences)


def bootstrap_p_value(differences: np.ndarray) -> float:
    below = (np.sum(differences <= 0) + 1) / (differences.size + 1)
    above = (np.sum(differences >= 0) + 1) / (differences.size + 1)
    return float(min(1.0, 2 * min(below, above)))


def holm(p_values: list[float]) -> list[float]:
    order = np.argsort(p_values)
    adjusted = np.empty(len(p_values))
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (len(p_values) - rank) * p_values[i]))
        adjusted[i] = running
    return adjusted.tolist()


def verdict(difference: float, p_holm: float, alpha: float = 0.05) -> str:
    if p_holm >= alpha:
        return "tie"
    return "ahead" if difference > 0 else "behind"


def broad(table: pl.DataFrame, rounds: int = BOOTSTRAP_ROUNDS) -> list[dict]:
    rows = []
    for column, title in dbnsfp.TOOLS.items():
        subset = table.filter(pl.col(column).is_not_null())
        labels = subset["label"].to_numpy().astype(int)
        if not 0 < labels.sum() < labels.size:
            continue
        theirs = subset[column].to_numpy()
        ours = subset["v2"].to_numpy()
        genes = subset["gene_name"].to_numpy()
        tool_roc, our_roc = roc_auc_score(labels, theirs), roc_auc_score(labels, ours)
        differences = paired_differences(labels, ours, theirs, genes, rounds)
        rows.append(
            {
                "tool": title,
                "column": column,
                "variants": subset.height,
                "tool_roc_auc": round(float(tool_roc), 4),
                "v2_roc_auc": round(float(our_roc), 4),
                "v2_no_frequency_roc_auc": round(
                    float(roc_auc_score(labels, subset["v2_no_frequency"].to_numpy())), 4
                ),
                "difference": round(float(our_roc - tool_roc), 4),
                "difference_95ci": [
                    round(float(np.percentile(differences, q)), 4) for q in (2.5, 97.5)
                ],
                "p_value": bootstrap_p_value(differences),
            }
        )
    for row, p_holm in zip(rows, holm([row["p_value"] for row in rows]), strict=True):
        row["p_holm"] = p_holm
        row["verdict"] = verdict(row["difference"], p_holm)
    return sorted(rows, key=lambda row: -row["tool_roc_auc"])


def dbnsfp_scores(path: Path, md5: str, keys: pl.DataFrame, cache: Path) -> pl.DataFrame:
    if not cache.exists():
        if digests_of(path)[1] != md5:
            raise ChecksumMismatch(f"{path.name}: md5 does not match {md5}")
        cache.parent.mkdir(parents=True, exist_ok=True)
        dbnsfp.extract(path, path.with_name(path.name + ".tbi"), keys).write_parquet(cache)
    return pl.read_parquet(cache).group_by(KEY).agg(pl.col(column).max() for column in dbnsfp.TOOLS)


def scored(table: pl.DataFrame, models_dir: Path, blocks_path: Path, final: dict) -> pl.DataFrame:
    scores = predictions(table, models_dir, blocks_path, final)
    return table.with_columns(
        pl.Series("v2", scores["final"]), pl.Series("v2_no_frequency", scores["no_frequency"])
    )


def run(
    splits_dir: Path,
    models_dir: Path,
    blocks_path: Path,
    alphamissense_path: Path,
    revel_zip: Path,
    interim_dir: Path,
    reports_dir: Path,
    dbnsfp_path: Path | None = None,
    dbnsfp_md5: str = "",
) -> dict:
    started = time.monotonic()
    final = json.loads((models_dir / "v2_final.json").read_text())
    pinned = ["v2_final.json", *final["members"].values(), *final["ablation_models"].values()]
    hashes_before = {name: file_hash(models_dir / name) for name in pinned}
    missense = pl.col("consequence") == "missense"
    test = load_joint(splits_dir, "test", allow_locked=True).filter(missense)
    temporal = load_joint(splits_dir, "temporal", allow_locked=True).filter(missense)
    sets = {
        "test": test,
        "temporal_all": temporal,
        "temporal_test_gene": temporal.filter(pl.col("temporal_subset") == "test_gene"),
    }
    keys = pl.concat([t.select(*KEY, "transcript_id") for t in sets.values()]).unique(KEY)
    external = alphamissense_scores(alphamissense_path, keys).join(
        revel_scores(revel_csv(revel_zip, interim_dir), keys), on=KEY, how="full", coalesce=True
    )
    broad_scores = (
        dbnsfp_scores(
            dbnsfp_path,
            dbnsfp_md5,
            keys,
            interim_dir.parent / "dbnsfp" / "scores.parquet",
        )
        if dbnsfp_path
        else None
    )
    results = {}
    for name, table in sets.items():
        joined = scored(table, models_dir, blocks_path, final).join(external, on=KEY, how="left")
        common = joined.filter(
            pl.col("alphamissense").is_not_null() & pl.col("revel").is_not_null()
        )
        results[name] = {
            "variants": table.height,
            "with_alphamissense": int(joined["alphamissense"].is_not_null().sum()),
            "with_revel": int(joined["revel"].is_not_null().sum()),
            "common": common.height,
            "common_pathogenic": int(common["label"].sum()),
            "common_genes": common["gene_name"].n_unique(),
            "alphamissense_from_other_transcript": int((~common["alphamissense_on_mane"]).sum()),
            "revel_from_other_transcript": int((~common["revel_on_mane"]).sum()),
            "missing_by_class": {
                tool: {
                    str(label): int(count)
                    for label, count in joined.filter(pl.col(tool).is_null())
                    .group_by("label")
                    .len()
                    .iter_rows()
                }
                for tool in ("alphamissense", "revel")
            },
        } | compare(common)
        if broad_scores is not None:
            results[name]["dbnsfp"] = broad(joined.join(broad_scores, on=KEY, how="left"))
    hashes_after = {name: file_hash(models_dir / name) for name in pinned}
    report = {
        "scope": "missense SNVs; every scorer evaluated on the same variants",
        "external_scores_used_for": "evaluation only; never an input, feature, label or filter",
        "sources": {
            "alphamissense": {
                "file": alphamissense_path.name,
                "sha256": file_hash(alphamissense_path),
                "license": "CC BY-NC-SA 4.0",
            },
            "revel": {
                "file": revel_zip.name,
                "sha256": file_hash(revel_zip),
                "license": "free for non-commercial use",
            },
        },
        "dbnsfp": {"file": dbnsfp_path.name, "license": "CC BY-NC-ND 4.0 (academic branch)"}
        if dbnsfp_path
        else None,
        "alphamissense_thresholds": [ALPHAMISSENSE_BENIGN_BELOW, ALPHAMISSENSE_PATHOGENIC_ABOVE],
        "model_hashes": hashes_before,
        "model_changed": hashes_before != hashes_after,
        "results": results,
        "seconds": round(time.monotonic() - started, 1),
    }
    (reports_dir / "external_comparison.json").write_text(json.dumps(report, indent=2))
    (reports_dir / "external_comparison.md").write_text(markdown(report), encoding="utf-8")
    return report


THREE_WAY_TITLES = {
    "v2": "Ours: v2 (0.10 / 0.90)",
    "alphamissense": "AlphaMissense (0.34 / 0.564)",
}
VERDICT_TITLES = {"ahead": "ours ahead", "tie": "no clear difference", "behind": "tool ahead"}
SET_TITLES = {
    "test": "Gene-disjoint test",
    "temporal_all": "Temporal (added to ClinVar after October 2025), all genes",
    "temporal_test_gene": "Temporal, test genes only",
}


def fmt(value) -> str:
    return "–" if value is None else f"{value:.1%}"


def markdown(report: dict) -> str:
    lines = [
        "# External comparison: AlphaMissense and REVEL",
        "",
        "For research and educational use only; not for clinical diagnosis or treatment decisions.",
        "",
        "AlphaMissense and REVEL scores were downloaded only to evaluate them next to our model "
        "on the same locked missense variants. They are never a model input, feature, label or "
        f"filter, and our models did not change (model changed: {report['model_changed']}).",
        "",
        "Read the numbers with these caveats:",
        "",
        "- Our v2 model uses gnomAD allele frequency as an input. AlphaMissense and REVEL do not, "
        "but AlphaMissense learned its benign signal from human and primate population frequency, "
        "so it is frequency-informed through its training labels. The frequency-free row removes "
        "our direct use of it.",
        "- Our test genes were never seen by our models. REVEL (trained on HGMD and rare ESP "
        "variants, from predictors trained on HGMD, SwissVar and ClinVar) and AlphaMissense "
        "(calibrated on ClinVar) may have seen these variants or genes; gene-disjointness does not "
        "apply to them.",
        "- Both tools predate the temporal labels, but since the 2022 ClinGen guidance submitters "
        "use REVEL and AlphaMissense as ACMG PP3/BP4 evidence, so some newer ClinVar labels partly "
        "rest on these very scores. Both effects favour the external tools.",
        "- When a tool has no score on our MANE transcript, its highest score on another "
        "transcript is used; the count is given under each set.",
        "- dbNSFP tables use each tool's rankscore (oriented so higher means more damaging, taken "
        "from its most damaging transcript). Several tools there, such as ClinPred, MetaRNN and "
        "BayesDel with AF, use allele frequency and were trained on ClinVar labels, so on the "
        "gene-disjoint test set they may have seen our test variants. The temporal sets guard "
        "against that only if a tool's training data predates October 2025, which we could not "
        "verify tool by tool.",
    ]
    for name, result in report["results"].items():
        lines += [
            "",
            f"## {SET_TITLES[name]}",
            "",
            f"{result['common']:,} of {result['variants']:,} variants scored by all four "
            f"({result['common_pathogenic']:,} pathogenic, {result['common_genes']:,} genes). "
            f"Scores from a non-MANE transcript: AlphaMissense "
            f"{result['alphamissense_from_other_transcript']:,}, REVEL "
            f"{result['revel_from_other_transcript']:,}.",
            "",
            "| Scorer | ROC-AUC [95% gene CI] | PR-AUC [95% gene CI] |",
            "|---|---|---|",
        ]
        for key, title in SCORERS.items():
            m = result["metrics"][key]
            roc, pr = m["gene_bootstrap_95ci"]["roc_auc"], m["gene_bootstrap_95ci"]["pr_auc"]
            lines.append(
                f"| {title} | {m['roc_auc']:.4f} [{roc[0]:.4f}–{roc[1]:.4f}] | "
                f"{m['pr_auc']:.4f} [{pr[0]:.4f}–{pr[1]:.4f}] |"
            )
        lines += [
            "",
            "| Three-way call | Decided | Likely pathogenic correct | Likely benign correct | "
            "Pathogenic found | Benign found |",
            "|---|---|---|---|---|---|",
        ]
        for key, title in THREE_WAY_TITLES.items():
            t = result["three_way"][key]
            cells = [
                t["decided_share"],
                t["likely_pathogenic_precision"],
                t["likely_benign_precision"],
                t["pathogenic_found"],
                t["benign_found"],
            ]
            lines.append(f"| {title} | " + " | ".join(fmt(c) for c in cells) + " |")
        if "dbnsfp" in result:
            lines += [
                "",
                f"### {SET_TITLES[name]}: every predictor in dbNSFP 5.4a",
                "",
                "Each tool is scored on the variants it covers, and our model on exactly those "
                "variants. The difference is ours minus the tool, with a paired 95% gene-bootstrap "
                "interval (unadjusted). The verdict uses the two-sided bootstrap p-value after "
                f"Holm correction across all {len(result['dbnsfp'])} tools in this set "
                f"(alpha 0.05). With {BOOTSTRAP_ROUNDS:,} resamples a p-value cannot fall below "
                f"2/{BOOTSTRAP_ROUNDS + 1:,}, so the smallest Holm p here is "
                f"{min(1.0, len(result['dbnsfp']) * 2 / (BOOTSTRAP_ROUNDS + 1)):.4f}: "
                "no resample crossed zero.",
                "",
                "Verdicts: "
                + ", ".join(
                    f"{label} {sum(row['verdict'] == key for row in result['dbnsfp'])}"
                    for key, label in VERDICT_TITLES.items()
                )
                + ".",
                "",
                "| Tool | Variants | Tool ROC-AUC | Ours v2 | Difference [95% CI] "
                "| Holm p | Verdict | Ours without AF |",
                "|---|---|---|---|---|---|---|---|",
            ]
            for row in result["dbnsfp"]:
                low, high = row["difference_95ci"]
                lines.append(
                    f"| {row['tool']} | {row['variants']:,} | {row['tool_roc_auc']:.4f} | "
                    f"{row['v2_roc_auc']:.4f} | {row['difference']:+.4f} "
                    f"[{low:+.4f}, {high:+.4f}] | {row['p_holm']:.4f} | "
                    f"{VERDICT_TITLES[row['verdict']]} | {row['v2_no_frequency_roc_auc']:.4f} |"
                )
    return "\n".join(lines) + "\n"
