import argparse
import sys
import time

import polars as pl
from pyfaidx import Fasta

from healthcare import __version__, doctor
from healthcare.annotate.indel import annotate_variant, normalize, parse_variant
from healthcare.annotate.snv import TranscriptIndex, parse_snv
from healthcare.data import annotations, clinvar, conservation, gnomad, mane, uniprot
from healthcare.data.download import ChecksumMismatch, download
from healthcare.data.sources import load_sources
from healthcare.features import missense, phase2
from healthcare.models import baseline
from healthcare.paths import ARTIFACTS_DIR, DATA_DIR, PROJECT_ROOT, REPORTS_DIR
from healthcare.runlog import start_log, write_summary
from healthcare.split import genes as gene_split


def run_download(names: list[str], force: bool) -> int:
    sources = load_sources()
    unknown = sorted(set(names) - sources.keys())
    if unknown:
        print(f"unknown sources: {', '.join(unknown)}; known: {', '.join(sorted(sources))}")
        return 2
    log_path = start_log("download")
    results = {}
    failed = False
    for name in names or sorted(sources):
        source = sources[name]
        try:
            results[name] = {"version": source.version} | download(
                source.url, source.path, source.sha256, force, source.md5
            )
        except (ChecksumMismatch, OSError) as error:
            results[name] = {"status": "failed", "error": str(error)}
            failed = True
    summary_path = write_summary("download", {"log": log_path, "sources": results})
    for name, result in results.items():
        print(f"{name}: {result['status']}")
    print(f"summary: {summary_path}")
    return 1 if failed else 0


def run_build_transcripts(chrom: str | None, limit: int | None, force: bool) -> int:
    sources = load_sources()
    parts = [f"chr{chrom}" if chrom else "", f"n{limit}" if limit else ""]
    suffix = "".join(f".{part}" for part in parts if part)
    output = DATA_DIR / "processed" / f"transcripts{suffix}.parquet"
    if output.exists() and not force:
        print(f"skipped: {output} exists (use --force)")
        return 0
    required = ("grch38_fasta", "mane_gtf", "mane_protein")
    missing = [name for name in required if not sources[name].path.exists()]
    if missing:
        print(f"missing raw files, run: healthcare download {' '.join(missing)}")
        return 2
    log_path = start_log("build-transcripts")
    summary = {"log": log_path} | mane.build(
        sources["grch38_fasta"].path,
        sources["mane_gtf"].path,
        sources["mane_protein"].path,
        DATA_DIR / "interim" / "reference" / "GRCh38.primary_assembly.fa",
        output,
        chrom,
        limit,
    )
    summary_path = write_summary("build-transcripts", summary)
    print(f"transcripts: {summary['transcripts']}, match rate: {summary['translation_match_rate']}")
    print(f"summary: {summary_path}")
    return 0


def run_build_clinvar(force: bool) -> int:
    sources = load_sources()
    out_dir = DATA_DIR / "processed" / "clinvar"
    reference = DATA_DIR / "interim" / "reference" / "GRCh38.primary_assembly.fa"
    if (out_dir / "train_labels.parquet").exists() and not force:
        print(f"skipped: {out_dir} exists (use --force)")
        return 0
    required = ("clinvar_2025", "clinvar_2026", "clinvar_summary_2025")
    missing = [name for name in required if not sources[name].path.exists()]
    if missing or not reference.exists():
        print("missing inputs: run `healthcare download` and `healthcare build-transcripts` first")
        return 2
    log_path = start_log("build-clinvar")
    summary = {"log": log_path} | clinvar.build(
        sources["clinvar_2025"].path,
        sources["clinvar_2026"].path,
        sources["clinvar_summary_2025"].path,
        reference,
        out_dir,
    )
    summary_path = write_summary("build-clinvar", summary)
    print(f"train labels: {summary['train_labels']}, temporal labels: {summary['temporal_labels']}")
    print(f"summary: {summary_path}")
    return 0


TRANSCRIPTS = DATA_DIR / "processed" / "transcripts.parquet"
REFERENCE = DATA_DIR / "interim" / "reference" / "GRCh38.primary_assembly.fa"


def run_build_annotations(limit: int | None, force: bool) -> int:
    suffix = f".n{limit}" if limit else ""
    out_dir = DATA_DIR / "processed" / f"annotations{suffix}"
    if (out_dir / "train_labels.parquet").exists() and not force:
        print(f"skipped: {out_dir} exists (use --force)")
        return 0
    clinvar_dir = DATA_DIR / "processed" / "clinvar"
    if not TRANSCRIPTS.exists() or not (clinvar_dir / "train_labels.parquet").exists():
        print("missing inputs: run `healthcare build-transcripts` and `healthcare build-clinvar`")
        return 2
    log_path = start_log("build-annotations")
    summary = {"log": log_path} | annotations.build(
        TRANSCRIPTS, REFERENCE, clinvar_dir, out_dir, limit
    )
    summary_path = write_summary("build-annotations", summary)
    print(f"missense agreement: {summary['train_labels'].get('missense_agreement')}")
    print(f"summary: {summary_path}")
    return 0


def run_annotate(variant: str) -> int:
    try:
        chrom, pos, ref, alt = parse_variant(variant)
    except ValueError as error:
        print(error)
        return 2
    sequence = Fasta(str(REFERENCE))[chrom]
    if pos + len(ref) - 1 > len(sequence):
        print(f"position {pos} is beyond the end of chromosome {chrom} ({len(sequence)})")
        return 2
    genome_ref = sequence[pos - 1 : pos - 1 + len(ref)].seq.upper()
    if genome_ref != ref:
        print(f"REF mismatch: the reference genome has {genome_ref} at {chrom}:{pos}, not {ref}")
        return 2
    normalized = normalize(pos, ref, alt, lambda a, b: sequence[a - 1 : b].seq.upper())
    if normalized != (pos, ref, alt):
        pos, ref, alt = normalized
        print(f"left-normalized to {chrom}:{pos} {ref}>{alt}")
    index = TranscriptIndex(annotations.load_transcripts(TRANSCRIPTS))
    results = annotate_variant(index, chrom, pos, ref, alt)
    if not results:
        print("no MANE Select transcript overlaps this position (out of scope)")
    for result in results:
        details = f" {result.hgvs_p} (c.{result.cds_position})" if result.hgvs_p else ""
        print(f"{result.gene_name} {result.transcript_id}: {result.consequence}{details}")
    return 0


def run_build_gnomad(chrom: str | None, limit: int | None, force: bool, workers: int) -> int:
    sources = load_sources()
    names = {
        name.removeprefix("gnomad_exomes_chr"): source.path
        for name, source in sources.items()
        if name.startswith("gnomad_exomes_chr")
    }
    if chrom:
        names = {chrom: names[chrom]}
    required = ("gnomad_constraint", "gnomad_exomes_an_all_sites")
    if not all(sources[name].path.exists() for name in required) or not TRANSCRIPTS.exists():
        print(
            f"missing inputs: run `healthcare download {' '.join(required)}` and build-transcripts"
        )
        return 2
    suffix = f".n{limit}" if limit else ""
    log_path = start_log("build-gnomad")
    summary = {"log": log_path} | gnomad.build(
        names,
        sources["gnomad_constraint"].path,
        sources["gnomad_exomes_an_all_sites"].path,
        pl.read_parquet(TRANSCRIPTS),
        DATA_DIR / "processed" / f"gnomad{suffix}",
        workers,
        limit,
        force,
    )
    summary_path = write_summary("build-gnomad", summary)
    done, rows = len(summary["chromosomes_done"]), summary["exome_rows_total"]
    print(f"chromosomes done: {done}/24, rows: {rows}")
    print(f"summary: {summary_path}")
    return 0


def run_frequency(variant: str) -> int:
    try:
        chrom, pos, ref, alt = parse_snv(variant)
    except ValueError as error:
        print(error)
        return 2
    started = time.perf_counter()
    result = gnomad.lookup_frequency(
        DATA_DIR / "processed" / "gnomad" / "exomes", chrom, pos, ref, alt
    )
    elapsed_ms = (time.perf_counter() - started) * 1000
    print(result if result else "not observed in gnomAD v4.1 exomes")
    print(f"query: {elapsed_ms:.1f} ms")
    return 0


CONSERVATION_DIR = DATA_DIR / "processed" / "conservation"
ESM_STORE = ARTIFACTS_DIR / "esm"


def run_build_conservation(chrom: str | None, force: bool, workers: int) -> int:
    sources = load_sources()
    chroms = [chrom] if chrom else [*map(str, range(1, 23)), "X", "Y"]
    paths = {c: {t: sources[f"{t}_chr{c}"].path for t in conservation.TRACKS} for c in chroms}
    if not TRANSCRIPTS.exists():
        print("missing inputs: run `healthcare build-transcripts` first")
        return 2
    log_path = start_log("build-conservation")
    summary = {"log": log_path} | conservation.build(
        paths, pl.read_parquet(TRANSCRIPTS), CONSERVATION_DIR, workers, force
    )
    summary_path = write_summary("build-conservation", summary)
    print(f"chromosomes done: {summary['chromosomes_done']}/24, rows: {summary['rows_total']}")
    print(f"summary: {summary_path}")
    return 0


def run_conservation(variant: str) -> int:
    try:
        chrom, pos, _, _ = parse_snv(variant)
    except ValueError as error:
        print(error)
        return 2
    started = time.perf_counter()
    result = conservation.lookup_conservation(CONSERVATION_DIR, chrom, pos)
    elapsed_ms = (time.perf_counter() - started) * 1000
    print(result if result else "position is outside MANE CDS +/-20 bp")
    print(f"query: {elapsed_ms:.1f} ms")
    return 0


UNIPROT_FEATURES = DATA_DIR / "processed" / "uniprot" / "features.parquet"
UNIPROT_FUNCTIONS = DATA_DIR / "processed" / "uniprot" / "functions.parquet"


def run_build_uniprot(force: bool) -> int:
    source = load_sources()["uniprot_human"]
    if not source.path.exists() or not TRANSCRIPTS.exists():
        print("missing inputs: run `healthcare download uniprot_human` and build-transcripts")
        return 2
    transcripts = pl.read_parquet(TRANSCRIPTS)
    if force or not UNIPROT_FUNCTIONS.exists():
        functions = uniprot.build_functions(source.path, transcripts, UNIPROT_FUNCTIONS)
        write_summary("build-uniprot-functions", functions)
        print(f"proteins with a function summary: {functions['proteins_with_function']}")
    if UNIPROT_FEATURES.exists() and not force:
        print(f"skipped: {UNIPROT_FEATURES} exists (use --force)")
        return 0
    log_path = start_log("build-uniprot")
    summary = {"log": log_path} | uniprot.build(source.path, transcripts, UNIPROT_FEATURES)
    summary_path = write_summary("build-uniprot", summary)
    print(f"matched MANE proteins: {summary['match_rate']:.1%}, features: {summary['features']}")
    print(f"summary: {summary_path}")
    return 0


def run_protein_features(gene: str) -> int:
    table = pl.read_parquet(UNIPROT_FEATURES).filter(pl.col("gene_name") == gene.upper())
    if table.is_empty():
        print(f"no curated UniProt features mapped for {gene}")
        return 0
    with pl.Config(tbl_rows=60, tbl_formatting="ASCII_MARKDOWN", fmt_str_lengths=40):
        print(
            table.filter(~pl.col("feature_type").is_in(["HELIX", "STRAND", "TURN"])).select(
                "feature_type", "start", "end", "note"
            )
        )
    return 0


def feature_inputs() -> missense.Inputs:
    processed = DATA_DIR / "processed"
    return missense.Inputs(
        transcripts=TRANSCRIPTS,
        gnomad_exomes=processed / "gnomad" / "exomes",
        coverage=processed / "gnomad" / "coverage.parquet",
        conservation=CONSERVATION_DIR,
        constraint=processed / "gnomad" / "constraint.parquet",
        uniprot=UNIPROT_FEATURES,
        esm_store=ESM_STORE,
    )


def missing_inputs(inputs: missense.Inputs) -> list[str]:
    return [str(p) for p in vars(inputs).values() if not p.exists()]


def run_build_features(force: bool) -> int:
    out_dir = DATA_DIR / "processed" / "features"
    if (out_dir / "missense_train_labels.parquet").exists() and not force:
        print(f"skipped: {out_dir} exists (use --force)")
        return 0
    processed = DATA_DIR / "processed"
    inputs = feature_inputs()
    missing = missing_inputs(inputs)
    if missing:
        print(f"missing inputs: {', '.join(missing)}")
        return 2
    log_path = start_log("build-features")
    summary = {"log": log_path} | missense.build(
        processed / "annotations", processed / "clinvar", inputs, out_dir
    )
    summary_path = write_summary("build-features", summary)
    train = summary["train_labels"]
    print(f"train rows: {train['rows']} ({train['pathogenic']} pathogenic)")
    print(f"features: {summary['feature_columns']}")
    print(f"summary: {summary_path}")
    return 0


def run_build_phase2_features(limit: int | None, force: bool) -> int:
    suffix = f".n{limit}" if limit else ""
    out_dir = DATA_DIR / "processed" / f"features{suffix}"
    if (out_dir / "phase2_train_labels.parquet").exists() and not force:
        print(f"skipped: {out_dir} exists (use --force)")
        return 0
    inputs = feature_inputs()
    missing = missing_inputs(inputs)
    if missing:
        print(f"missing inputs: {', '.join(missing)}")
        return 2
    log_path = start_log("build-phase2-features")
    processed = DATA_DIR / "processed"
    summary = {"log": log_path} | phase2.build(
        processed / "annotations",
        processed / "clinvar",
        inputs,
        REFERENCE,
        out_dir,
        REPORTS_DIR / f"phase2_missing_values{suffix}.json",
        limit,
    )
    summary_path = write_summary("build-phase2-features", summary)
    print(f"train rows: {summary['train_labels']['rows']}, features: {summary['feature_columns']}")
    print(f"summary: {summary_path}")
    return 0


def run_build_split(force: bool) -> int:
    features_dir = DATA_DIR / "processed" / "features"
    if not (features_dir / "missense_train_labels.parquet").exists():
        print("missing inputs: run `healthcare build-features` first")
        return 2
    log_path = start_log("build-split")
    try:
        summary = gene_split.build(
            features_dir,
            TRANSCRIPTS,
            DATA_DIR / "processed" / "splits",
            REPORTS_DIR / "split_summary.json",
            force,
        )
    except RuntimeError as error:
        print(error)
        return 1
    write_summary("build-split", {"log": log_path} | summary)
    for name, stats in summary["splits"].items():
        share = f"{stats['pathogenic_fraction']:.1%}"
        print(f"{name}: {stats['genes']} genes, {stats['variants']} variants, {share} P")
    print(f"overlaps: {summary['overlaps']}")
    print(f"report: {REPORTS_DIR / 'split_summary.json'}")
    return 0


def run_build_phase2_split(force: bool) -> int:
    processed = DATA_DIR / "processed"
    if not (processed / "features" / "phase2_train_labels.parquet").exists():
        print("missing inputs: run `healthcare build-phase2-features` first")
        return 2
    log_path = start_log("build-phase2-split")
    try:
        summary = gene_split.build_phase2(
            processed / "features",
            processed / "annotations",
            processed / "splits",
            REPORTS_DIR / "split_summary.json",
            REPORTS_DIR / "phase2_split_summary.json",
            force,
        )
    except RuntimeError as error:
        print(error)
        return 1
    write_summary("build-phase2-split", {"log": log_path} | summary)
    for name, stats in summary["splits"].items():
        print(f"{name}: {stats['genes']} genes, {stats['variants']} variants")
    print(
        f"overlaps: {summary['overlaps']}, missense cross-table: {summary['missense_cross_table']}"
    )
    print(f"report: {REPORTS_DIR / 'phase2_split_summary.json'}")
    return 0


def run_train_joint_baseline(trials: int) -> int:
    from healthcare.models import joint

    splits_dir = DATA_DIR / "processed" / "splits"
    if not (splits_dir / "phase2_open.parquet").exists():
        print("missing inputs: run `healthcare build-phase2-split` first")
        return 2
    log_path = start_log("train-joint-baseline")
    artifacts = ARTIFACTS_DIR / "models"
    artifacts.mkdir(parents=True, exist_ok=True)
    report = {"log": log_path} | joint.run(splits_dir, artifacts, REPORTS_DIR, trials)
    write_summary("train-joint-baseline", report)
    for consequence, metrics in report["validation_by_type"].items():
        print(f"{consequence}: ROC-AUC {metrics.get('roc_auc', 'n/a')} (n={metrics['n']})")
    print(f"report: {REPORTS_DIR / 'joint_lightgbm.json'}")
    return 0


def run_choose_final_v2() -> int:
    from healthcare.evaluate import all_types

    log_path = start_log("choose-final-v2")
    summary = {"log": log_path} | all_types.choose(
        DATA_DIR / "processed" / "splits", JOINT_ESM_BLOCKS, ARTIFACTS_DIR / "models", REPORTS_DIR
    )
    write_summary("choose-final-v2", summary)
    print(f"routes: {summary['routes']}")
    for consequence, metrics in summary["validation_calibrated"].items():
        print(f"{consequence}: validation ECE {metrics.get('ece', 'n/a')}")
    return 0


def run_evaluate_final_v2(force: bool) -> int:
    from healthcare.evaluate import all_types

    if not (ARTIFACTS_DIR / "models" / "v2_final.json").exists():
        print("missing inputs: run `healthcare choose-final-v2` first")
        return 2
    log_path = start_log("evaluate-final-v2")
    try:
        report = {"log": log_path} | all_types.evaluate(
            DATA_DIR / "processed" / "splits",
            ESM_STORE,
            ARTIFACTS_DIR / "esm_variants" / "joint_locked.npy",
            ARTIFACTS_DIR / "models",
            REPORTS_DIR,
            force,
        )
    except all_types.AlreadyOpenedError as error:
        print(error)
        return 1
    write_summary(
        "evaluate-final-v2", {k: report[k] for k in ("log", "model_changed_after_test", "seconds")}
    )
    print(f"model changed after test: {report['model_changed_after_test']}")
    print(f"report: {REPORTS_DIR / 'all_types_v2.md'}")
    return 0


def run_export_onnx_v2() -> int:
    from healthcare.models import export_v2

    log_path = start_log("export-onnx-v2")
    report = {"log": log_path} | export_v2.run(
        DATA_DIR / "processed" / "splits", JOINT_ESM_BLOCKS, ARTIFACTS_DIR / "models", REPORTS_DIR
    )
    write_summary("export-onnx-v2", report)
    print(f"{report['file']}: {report['megabytes']} MB, parity {report['parity_max_abs_diff']}")
    print(f"decision mismatches: {report['parity_decision_mismatches']}")
    return 0


def run_train_baseline(trials: int) -> int:
    splits_dir = DATA_DIR / "processed" / "splits"
    if not (splits_dir / "missense_open.parquet").exists():
        print("missing inputs: run `healthcare build-split` first")
        return 2
    log_path = start_log("train-baseline")
    report = baseline.run(splits_dir, ARTIFACTS_DIR / "models", REPORTS_DIR, trials)
    write_summary("train-baseline", {"log": log_path} | report)
    validation = report["validation"]
    print(f"validation ROC-AUC {validation['roc_auc']:.4f}, PR-AUC {validation['pr_auc']:.4f}")
    print(f"report: {REPORTS_DIR / 'baseline_lightgbm.json'}")
    return 0


def run_esm_trial(gene: str, change: str | None) -> int:
    from healthcare.esm import trial

    log_path = start_log("esm-trial")
    summary = {"log": log_path} | trial.run(
        TRANSCRIPTS,
        DATA_DIR / "processed" / "splits",
        ARTIFACTS_DIR / "hf",
        gene.upper(),
        trial.parse_protein_change(change) if change else None,
    )
    summary_path = write_summary("esm-trial", summary)
    print(f"{summary['gene']}: {summary['protein_length']} aa in {summary['seconds']} s")
    if "variant" in summary:
        print(f"{summary['variant']['change']}: LLR {summary['variant']['llr']}")
    print(f"peak GPU memory: {summary['peak_gpu_memory_gib']} GiB")
    print(f"summary: {summary_path}")
    return 0


def run_build_esm(limit: int | None) -> int:
    from healthcare.esm import store

    log_path = start_log("build-esm")
    summary = {"log": log_path} | store.build(TRANSCRIPTS, ESM_STORE, ARTIFACTS_DIR / "hf", limit)
    summary_path = write_summary("build-esm", summary)
    done, total = summary["proteins_done"], summary["proteins_requested"]
    print(f"proteins done: {done}/{total}, errors: {summary['errors']}")
    print(f"summary: {summary_path}")
    return 0


SWEEP_DIR = ARTIFACTS_DIR / "sweeps"
ESM_BLOCKS = ARTIFACTS_DIR / "esm_variants" / "open.npy"


def run_prepare_dual_tower() -> int:
    import polars as pl

    from healthcare.esm import variant_inputs

    splits_dir = DATA_DIR / "processed" / "splits"
    table = pl.read_parquet(splits_dir / "missense_open.parquet")
    log_path = start_log("prepare-dual-tower")
    summary = {"log": log_path} | variant_inputs.build(table, ESM_STORE, ESM_BLOCKS)
    write_summary("prepare-dual-tower", summary)
    rows, size, missing = summary["rows"], summary["gigabytes"], summary["missing"]
    print(f"ESM blocks: {rows} variants, {size} GB, missing {missing}")
    return 0


JOINT_ESM_BLOCKS = ARTIFACTS_DIR / "esm_variants" / "joint_open.npy"


def run_prepare_joint_dual_tower() -> int:
    from healthcare.models import joint_tower

    log_path = start_log("prepare-joint-dual-tower")
    summary = {"log": log_path} | joint_tower.prepare(
        DATA_DIR / "processed" / "splits", ESM_STORE, JOINT_ESM_BLOCKS
    )
    write_summary("prepare-joint-dual-tower", summary)
    rows, size, missing = summary["rows"], summary["gigabytes"], summary["missing"]
    print(f"ESM blocks: {rows} variants, {size} GB, without a protein position {missing}")
    return 0


def run_train_joint_dual_tower() -> int:
    from healthcare.models import joint_tower

    if not JOINT_ESM_BLOCKS.exists():
        print("missing inputs: run `healthcare prepare-joint-dual-tower` first")
        return 2
    log_path = start_log("train-joint-dual-tower")
    report = {"log": log_path} | joint_tower.run(
        DATA_DIR / "processed" / "splits", JOINT_ESM_BLOCKS, ARTIFACTS_DIR / "models", REPORTS_DIR
    )
    write_summary("train-joint-dual-tower", report)
    for consequence, metrics in report["validation_by_type_seeds"].items():
        print(f"{consequence}: ROC-AUC {metrics['roc_auc_mean']} ± {metrics['roc_auc_std']}")
    print(f"report: {REPORTS_DIR / 'joint_dual_tower.json'}")
    return 0


def dual_tower_data():
    from healthcare.models import dual_tower

    return dual_tower.Data(DATA_DIR / "processed" / "splits", ESM_BLOCKS)


def run_sweep_dual_tower(hours: float, trials: int | None) -> int:
    from healthcare.keepawake import prevent_idle_sleep
    from healthcare.models import sweep

    prevent_idle_sleep()

    log_path = start_log("sweep-dual-tower")
    summary = {"log": log_path} | sweep.run_sweep(
        dual_tower_data(),
        SWEEP_DIR / "dual_tower.journal",
        SWEEP_DIR / "dual_tower_trials.jsonl",
        hours,
        trials,
    )
    summary_path = write_summary("sweep-dual-tower", summary)
    total, now = summary["trials_completed_total"], summary["trials_this_run"]
    print(f"trials completed: {total} (this run {now})")
    if summary["best"]:
        print(f"best validation ROC-AUC so far: {summary['best'][0]['roc_auc']:.4f}")
    print(f"summary: {summary_path}")
    return 0


def run_confirm_dual_tower() -> int:
    from healthcare.keepawake import prevent_idle_sleep
    from healthcare.models import sweep

    prevent_idle_sleep()

    log_path = start_log("confirm-dual-tower")
    report = sweep.run_confirmation(
        dual_tower_data(), SWEEP_DIR / "dual_tower.journal", ARTIFACTS_DIR / "models", REPORTS_DIR
    )
    write_summary("confirm-dual-tower", {"log": log_path} | report)
    roc = report["winner"]["seeds"]["roc_auc"]
    print(f"winner validation ROC-AUC {roc['mean']:.4f} +/- {roc['std']:.4f}")
    print(f"report: {REPORTS_DIR / 'dual_tower.json'}")
    return 0


def run_calibrate() -> int:
    from healthcare.models import calibration

    log_path = start_log("calibrate")
    report = calibration.run(
        DATA_DIR / "processed" / "splits", ESM_BLOCKS, ARTIFACTS_DIR / "models", REPORTS_DIR
    )
    write_summary("calibrate", {"log": log_path} | report)
    chosen = report["chosen_method"]
    out_of_fold = report["methods_out_of_fold"][chosen]
    print(f"method: {chosen}, out-of-fold ECE {out_of_fold['metrics']['ece']:.4f}")
    print(f"uncertain share (out-of-fold): {out_of_fold['decisions']['uncertain']['share']:.1%}")
    print(f"report: {REPORTS_DIR / 'calibration.json'}")
    return 0


def run_export_onnx() -> int:
    from healthcare.models import export

    log_path = start_log("export-onnx")
    report = export.run(
        DATA_DIR / "processed" / "splits", ESM_BLOCKS, ARTIFACTS_DIR / "models", REPORTS_DIR
    )
    write_summary("export-onnx", {"log": log_path} | report)
    diff, timing = report["parity_max_abs_diff"], report["latency_cpu_single_variant"]
    print(f"onnx: {ARTIFACTS_DIR / 'models' / report['file']} ({report['megabytes']} MB)")
    mismatches = diff["decision_mismatches"]
    print(f"parity max |diff| {diff['probability']:.2e}, decision mismatches {mismatches}")
    print(f"single-variant CPU latency: median {timing['median_ms']} ms, p95 {timing['p95_ms']} ms")
    return 0


def run_evaluate_final() -> int:
    from healthcare.evaluate import final

    log_path = start_log("evaluate-final")
    report = final.run(
        DATA_DIR / "processed" / "splits",
        ESM_STORE,
        ESM_BLOCKS,
        ARTIFACTS_DIR / "models",
        REPORTS_DIR,
    )
    write_summary("evaluate-final", {"log": log_path, "seconds": report["seconds"]})
    for name in ("test", "temporal_all", "temporal_test_gene"):
        metrics = report["results"][name]["models"]["ensemble_calibrated"]
        print(f"{name}: ROC-AUC {metrics['roc_auc']:.4f}, PR-AUC {metrics['pr_auc']:.4f}")
    print(f"report: {REPORTS_DIR / (final.REPORT_NAME + '.md')}")
    return 0


def load_predictor():
    from healthcare.serve.pipeline import Predictor

    inputs = feature_inputs()
    model = ARTIFACTS_DIR / "models" / "all_types_v2.onnx"
    missing = missing_inputs(inputs) + ([] if model.exists() else [str(model)])
    if missing:
        print(f"missing inputs: {', '.join(missing)}")
        return None
    from healthcare.models.joint import load_joint

    train = load_joint(DATA_DIR / "processed" / "splits", "train")
    functions = (
        dict(pl.read_parquet(UNIPROT_FUNCTIONS).select("transcript_id", "function").iter_rows())
        if UNIPROT_FUNCTIONS.exists()
        else {}
    )
    return Predictor(inputs, model, train, REFERENCE, functions)


def print_result(result) -> None:
    from healthcare.serve.pipeline import Prediction

    if not isinstance(result, Prediction):
        named = result.gene is None or result.variant.upper().startswith(result.gene.upper())
        label = result.variant if named else f"{result.gene} {result.variant}"
        print(f"{label}: out of scope - {result.reason}")
        return
    features = result.features
    print(f"{result.gene} {result.hgvs_p}  {result.variant}  ({result.transcript_id})")
    print(
        f"  {result.decision.replace('_', ' ')}: probability {result.probability:.3f} "
        f"(LightGBM {result.lightgbm:.3f}, dual tower {result.dual_tower:.3f})"
    )
    shown = {
        "gnomAD AF": features["f_gnomad_af"],
        "phyloP": features["f_phylop"],
        "ESM LLR": features["f_esm_llr"],
        "Grantham": features["f_grantham"],
    }
    print(
        "  evidence: "
        + ", ".join(f"{k} {'n/a' if v is None else f'{v:.3g}'}" for k, v in shown.items())
    )
    print("  most influential evidence groups (log-odds vs. typical training variants):")
    for group, delta in result.contributions[:4]:
        direction = "toward pathogenic" if delta > 0 else "toward benign"
        print(f"    {group.replace('_', ' '):<20} {delta:+.2f}  {direction}")


def run_predict(words: list[str]) -> int:
    predictor = load_predictor()
    if predictor is None:
        return 2
    started = time.perf_counter()
    try:
        results = predictor.predict(" ".join(words))
    except ValueError as error:
        print(error)
        return 2
    for result in results:
        print_result(result)
    milliseconds = (time.perf_counter() - started) * 1000
    print(f"{milliseconds:.0f} ms; research and educational use only, not for clinical use")
    return 0


def run_serve(port: int) -> int:
    import json

    import uvicorn

    from healthcare.serve.api import create_app, model_card

    report_path = REPORTS_DIR / "all_types_v2.json"
    if not report_path.exists():
        print(f"missing inputs: {report_path}")
        return 2
    predictor = load_predictor()
    if predictor is None:
        return 2
    versions = json.loads((REPORTS_DIR / "missense_v1.json").read_text())["data_versions"]
    card = model_card(json.loads(report_path.read_text()), versions, predictor.thresholds)
    print(f"API docs: http://127.0.0.1:{port}/docs")
    uvicorn.run(create_app(predictor, card), host="127.0.0.1", port=port)
    return 0


def run_openapi() -> int:
    import json

    from healthcare.serve.api import openapi_schema

    path = PROJECT_ROOT / "web" / "openapi.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(openapi_schema(), indent=2) + "\n", encoding="utf-8")
    print(f"schema: {path}")
    return 0


def run_serve_check(count: int) -> int:
    from healthcare.serve.pipeline import EXAMPLES, Prediction, parity_check

    predictor = load_predictor()
    if predictor is None:
        return 2
    log_path = start_log("serve-check")
    from healthcare.models.joint import load_joint

    validation = load_joint(DATA_DIR / "processed" / "splits", "validation")
    per_type = max(1, count // validation["consequence"].n_unique())
    table = pl.concat(
        [
            group.sample(min(per_type, group.height), seed=20261016)
            for _, group in validation.group_by("consequence", maintain_order=True)
        ]
    )
    summary = {"log": log_path} | parity_check(predictor, table, JOINT_ESM_BLOCKS)
    summary["examples"] = {}
    for text in EXAMPLES:
        for result in predictor.predict(text):
            summary["examples"][text] = (
                {"variant": result.variant, "probability": round(result.probability, 4)}
                | {"decision": result.decision}
                if isinstance(result, Prediction)
                else {"variant": result.variant, "out_of_scope": result.reason}
            )
    summary_path = write_summary("serve-check", summary)
    print(
        f"feature mismatches: {summary['feature_mismatches']} over {summary['variants']} variants"
    )
    print(f"max probability diff vs batch ONNX: {summary['max_probability_diff']}")
    print(f"latency ms: {summary['latency_ms']}")
    print(f"summary: {summary_path}")
    return 0


def add_sample_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--chrom", help="Process a single chromosome, e.g. 22")
    parser.add_argument("--limit", type=int, help="Process at most N records")
    parser.add_argument("--force", action="store_true", help="Rebuild even if output exists")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="healthcare")
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("doctor", help="Check Python, PyTorch, CUDA, GPU and disk space")
    download_parser = commands.add_parser("download", help="Download pinned data sources")
    download_parser.add_argument("sources", nargs="*", help="Source names (default: all)")
    download_parser.add_argument("--force", action="store_true", help="Redownload even if present")
    transcripts_parser = commands.add_parser(
        "build-transcripts", help="Build the MANE Select transcript table and verify translations"
    )
    add_sample_args(transcripts_parser)
    clinvar_parser = commands.add_parser(
        "build-clinvar", help="Build ClinVar training, VUS, temporal and reclassification tables"
    )
    clinvar_parser.add_argument(
        "--force", action="store_true", help="Rebuild even if output exists"
    )
    annotations_parser = commands.add_parser(
        "build-annotations", help="Annotate ClinVar SNVs with our MANE-based annotator"
    )
    annotations_parser.add_argument("--limit", type=int, help="Annotate at most N rows per table")
    annotations_parser.add_argument("--force", action="store_true", help="Rebuild even if present")
    annotate_parser = commands.add_parser(
        "annotate", help="Annotate one SNV or indel, e.g. 17:43045712 G>A or 17:43124096 GA>G"
    )
    annotate_parser.add_argument("variant")
    gnomad_parser = commands.add_parser(
        "build-gnomad", help="Convert downloaded gnomAD exome VCFs and constraint to Parquet"
    )
    add_sample_args(gnomad_parser)
    gnomad_parser.add_argument("--workers", type=int, default=4, help="Parallel chromosomes")
    frequency_parser = commands.add_parser(
        "frequency", help="Look up a variant's gnomAD exome frequency, e.g. 17:43063931 G>A"
    )
    frequency_parser.add_argument("variant")
    conservation_parser = commands.add_parser(
        "build-conservation", help="Convert UCSC phyloP/phastCons wigFix to Parquet (CDS +/-20 bp)"
    )
    conservation_parser.add_argument("--chrom", help="Process a single chromosome, e.g. 22")
    conservation_parser.add_argument("--force", action="store_true", help="Rebuild even if present")
    conservation_parser.add_argument("--workers", type=int, default=6, help="Parallel chromosomes")
    lookup_parser = commands.add_parser(
        "conservation", help="Look up phyloP/phastCons at a position, e.g. 17:43063931 G>A"
    )
    lookup_parser.add_argument("variant")
    uniprot_parser = commands.add_parser(
        "build-uniprot", help="Map curated UniProt Swiss-Prot features onto MANE proteins"
    )
    uniprot_parser.add_argument("--force", action="store_true", help="Rebuild even if present")
    features_parser = commands.add_parser(
        "protein-features", help="Show curated UniProt regions for a gene, e.g. BRCA1"
    )
    features_parser.add_argument("gene")
    features_build_parser = commands.add_parser(
        "build-features", help="Build the missense feature tables (train, temporal, VUS)"
    )
    features_build_parser.add_argument("--force", action="store_true", help="Rebuild if present")
    phase2_parser = commands.add_parser(
        "build-phase2-features", help="Build feature tables for the Phase 2 variant types"
    )
    phase2_parser.add_argument("--limit", type=int, help="Use at most N variants per table")
    phase2_parser.add_argument("--force", action="store_true", help="Rebuild if present")
    split_parser = commands.add_parser(
        "build-split", help="Create (once) the gene-disjoint train/validation/test split"
    )
    split_parser.add_argument("--force", action="store_true", help="Regenerate the pinned split")
    phase2_split_parser = commands.add_parser(
        "build-phase2-split", help="Place Phase 2 variants into the pinned gene split"
    )
    phase2_split_parser.add_argument("--force", action="store_true", help="Re-pin if changed")
    joint_parser = commands.add_parser(
        "train-joint-baseline", help="Train LightGBM on missense and Phase 2 types together"
    )
    joint_parser.add_argument("--trials", type=int, default=20, help="Random-search trials")
    baseline_parser = commands.add_parser(
        "train-baseline", help="Train the LightGBM baseline with a small sweep (validation only)"
    )
    baseline_parser.add_argument("--trials", type=int, default=30, help="Random-search trials")
    esm_parser = commands.add_parser("esm-trial", help="Run ESM-2 on one MANE protein")
    esm_parser.add_argument("--gene", default="BRCA1")
    esm_parser.add_argument("--change", default="p.Arg1699Trp", help="Protein change to score")
    build_esm_parser = commands.add_parser(
        "build-esm", help="Precompute ESM-2 log-probabilities and embeddings for MANE proteins"
    )
    build_esm_parser.add_argument("--limit", type=int, help="Only the N shortest proteins")
    commands.add_parser("export-onnx-v2", help="Export the v2 all-types model to one ONNX file")
    commands.add_parser(
        "choose-final-v2", help="Choose per-type routes and calibration on validation only"
    )
    evaluate_v2_parser = commands.add_parser(
        "evaluate-final-v2", help="Open the locked test once and write reports/all_types_v2.md"
    )
    evaluate_v2_parser.add_argument("--force", action="store_true", help="Re-run a past opening")
    commands.add_parser(
        "prepare-joint-dual-tower", help="Precompute ESM blocks for the joint open splits"
    )
    commands.add_parser(
        "train-joint-dual-tower", help="Train the dual tower on missense and Phase 2 types"
    )
    commands.add_parser("prepare-dual-tower", help="Precompute ESM input blocks for open splits")
    sweep_parser = commands.add_parser("sweep-dual-tower", help="Resumable dual tower sweep")
    sweep_parser.add_argument("--hours", type=float, default=12.0, help="Time budget")
    sweep_parser.add_argument("--trials", type=int, help="Stop after N trials")
    commands.add_parser("confirm-dual-tower", help="Re-train top configs with 5 seeds")
    commands.add_parser("calibrate", help="Calibrate the ensemble on validation")
    commands.add_parser("export-onnx", help="Export the calibrated ensemble to one ONNX file")
    predict_parser = commands.add_parser(
        "predict", help="Predict one variant, e.g. 17:43057062 C>T or BRCA1 p.Arg1699Trp"
    )
    predict_parser.add_argument("variant", nargs="+")
    serve_parser = commands.add_parser("serve", help="Run the local API on 127.0.0.1")
    serve_parser.add_argument("--port", type=int, default=8000)
    commands.add_parser("openapi", help="Write the API schema to web/openapi.json")
    serve_check_parser = commands.add_parser(
        "serve-check", help="Check served features and latency against validation variants"
    )
    serve_check_parser.add_argument("--n", type=int, default=100, help="Validation variants")
    commands.add_parser(
        "evaluate-final", help="Open the locked test sets and write the final missense report"
    )
    args = parser.parse_args(argv)
    if args.command == "doctor":
        return doctor.run()
    if args.command == "download":
        return run_download(args.sources, args.force)
    if args.command == "build-transcripts":
        return run_build_transcripts(args.chrom, args.limit, args.force)
    if args.command == "build-clinvar":
        return run_build_clinvar(args.force)
    if args.command == "build-annotations":
        return run_build_annotations(args.limit, args.force)
    if args.command == "annotate":
        return run_annotate(args.variant)
    if args.command == "build-gnomad":
        return run_build_gnomad(args.chrom, args.limit, args.force, args.workers)
    if args.command == "frequency":
        return run_frequency(args.variant)
    if args.command == "build-conservation":
        return run_build_conservation(args.chrom, args.force, args.workers)
    if args.command == "conservation":
        return run_conservation(args.variant)
    if args.command == "build-uniprot":
        return run_build_uniprot(args.force)
    if args.command == "protein-features":
        return run_protein_features(args.gene)
    if args.command == "build-features":
        return run_build_features(args.force)
    if args.command == "build-phase2-features":
        return run_build_phase2_features(args.limit, args.force)
    if args.command == "export-onnx-v2":
        return run_export_onnx_v2()
    if args.command == "choose-final-v2":
        return run_choose_final_v2()
    if args.command == "evaluate-final-v2":
        return run_evaluate_final_v2(args.force)
    if args.command == "prepare-joint-dual-tower":
        return run_prepare_joint_dual_tower()
    if args.command == "train-joint-dual-tower":
        return run_train_joint_dual_tower()
    if args.command == "train-joint-baseline":
        return run_train_joint_baseline(args.trials)
    if args.command == "build-phase2-split":
        return run_build_phase2_split(args.force)
    if args.command == "build-split":
        return run_build_split(args.force)
    if args.command == "train-baseline":
        return run_train_baseline(args.trials)
    if args.command == "esm-trial":
        return run_esm_trial(args.gene, args.change)
    if args.command == "build-esm":
        return run_build_esm(args.limit)
    if args.command == "prepare-dual-tower":
        return run_prepare_dual_tower()
    if args.command == "sweep-dual-tower":
        return run_sweep_dual_tower(args.hours, args.trials)
    if args.command == "confirm-dual-tower":
        return run_confirm_dual_tower()
    if args.command == "calibrate":
        return run_calibrate()
    if args.command == "export-onnx":
        return run_export_onnx()
    if args.command == "predict":
        return run_predict(args.variant)
    if args.command == "serve":
        return run_serve(args.port)
    if args.command == "openapi":
        return run_openapi()
    if args.command == "serve-check":
        return run_serve_check(args.n)
    if args.command == "evaluate-final":
        return run_evaluate_final()
    return 2


if __name__ == "__main__":
    sys.exit(main())
