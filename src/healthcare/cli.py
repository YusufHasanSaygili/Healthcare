import argparse
import sys
import time

import polars as pl
from pyfaidx import Fasta

from healthcare import __version__, doctor
from healthcare.annotate.snv import TranscriptIndex, parse_snv
from healthcare.data import annotations, clinvar, conservation, gnomad, mane, uniprot
from healthcare.data.download import ChecksumMismatch, download
from healthcare.data.sources import load_sources
from healthcare.features import missense
from healthcare.paths import DATA_DIR
from healthcare.runlog import start_log, write_summary


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
    summary = {"log": log_path} | annotations.build(TRANSCRIPTS, clinvar_dir, out_dir, limit)
    summary_path = write_summary("build-annotations", summary)
    print(f"missense agreement: {summary['train_labels'].get('missense_agreement')}")
    print(f"summary: {summary_path}")
    return 0


def run_annotate(variant: str) -> int:
    try:
        chrom, pos, ref, alt = parse_snv(variant)
    except ValueError as error:
        print(error)
        return 2
    sequence = Fasta(str(REFERENCE))[chrom]
    if pos > len(sequence):
        print(f"position {pos} is beyond the end of chromosome {chrom} ({len(sequence)})")
        return 2
    genome_ref = sequence[pos - 1 : pos].seq.upper()
    if genome_ref != ref:
        print(f"REF mismatch: the reference genome has {genome_ref} at {chrom}:{pos}, not {ref}")
        return 2
    index = TranscriptIndex(annotations.load_transcripts(TRANSCRIPTS))
    results = index.annotate(chrom, pos, ref, alt)
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


def run_build_uniprot(force: bool) -> int:
    source = load_sources()["uniprot_human"]
    if UNIPROT_FEATURES.exists() and not force:
        print(f"skipped: {UNIPROT_FEATURES} exists (use --force)")
        return 0
    if not source.path.exists() or not TRANSCRIPTS.exists():
        print("missing inputs: run `healthcare download uniprot_human` and build-transcripts")
        return 2
    log_path = start_log("build-uniprot")
    summary = {"log": log_path} | uniprot.build(
        source.path, pl.read_parquet(TRANSCRIPTS), UNIPROT_FEATURES
    )
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


def run_build_features(force: bool) -> int:
    out_dir = DATA_DIR / "processed" / "features"
    if (out_dir / "missense_train_labels.parquet").exists() and not force:
        print(f"skipped: {out_dir} exists (use --force)")
        return 0
    processed = DATA_DIR / "processed"
    inputs = missense.Inputs(
        transcripts=TRANSCRIPTS,
        gnomad_exomes=processed / "gnomad" / "exomes",
        coverage=processed / "gnomad" / "coverage.parquet",
        conservation=CONSERVATION_DIR,
        constraint=processed / "gnomad" / "constraint.parquet",
        uniprot=UNIPROT_FEATURES,
    )
    missing = [str(p) for p in vars(inputs).values() if not p.exists()]
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
    annotate_parser = commands.add_parser("annotate", help="Annotate one SNV, e.g. 17:43045712 G>A")
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
    return 2


if __name__ == "__main__":
    sys.exit(main())
