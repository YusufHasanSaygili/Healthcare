import argparse
import sys

from healthcare import __version__, doctor
from healthcare.data import clinvar, mane
from healthcare.data.download import ChecksumMismatch, download
from healthcare.data.sources import load_sources
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
    args = parser.parse_args(argv)
    if args.command == "doctor":
        return doctor.run()
    if args.command == "download":
        return run_download(args.sources, args.force)
    if args.command == "build-transcripts":
        return run_build_transcripts(args.chrom, args.limit, args.force)
    if args.command == "build-clinvar":
        return run_build_clinvar(args.force)
    return 2


if __name__ == "__main__":
    sys.exit(main())
