import argparse
import sys

from healthcare import __version__, doctor
from healthcare.download import ChecksumMismatch, download
from healthcare.runlog import start_log, write_summary
from healthcare.sources import load_sources


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
                source.url, source.path, source.sha256, force
            )
        except (ChecksumMismatch, OSError) as error:
            results[name] = {"status": "failed", "error": str(error)}
            failed = True
    summary_path = write_summary("download", {"log": log_path, "sources": results})
    for name, result in results.items():
        print(f"{name}: {result['status']}")
    print(f"summary: {summary_path}")
    return 1 if failed else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="healthcare")
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("doctor", help="Check Python, PyTorch, CUDA, GPU and disk space")
    download_parser = commands.add_parser("download", help="Download pinned data sources")
    download_parser.add_argument("sources", nargs="*", help="Source names (default: all)")
    download_parser.add_argument("--force", action="store_true", help="Redownload even if present")
    args = parser.parse_args(argv)
    if args.command == "doctor":
        return doctor.run()
    if args.command == "download":
        return run_download(args.sources, args.force)
    return 2


if __name__ == "__main__":
    sys.exit(main())
