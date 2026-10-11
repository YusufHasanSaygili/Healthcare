import argparse

from healthcare import __version__, doctor


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="healthcare")
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("doctor", help="Check Python, PyTorch, CUDA, GPU and disk space")
    args = parser.parse_args(argv)
    if args.command == "doctor":
        return doctor.run()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
