import json
import logging
import time
from pathlib import Path

from healthcare.paths import LOGS_DIR

SUMMARY_DIR = LOGS_DIR / "summaries"


def start_log(name: str, logs_dir: Path = LOGS_DIR) -> Path:
    logs_dir.mkdir(parents=True, exist_ok=True)
    path = logs_dir / f"{name}-{time.strftime('%Y%m%d-%H%M%S')}.log"
    logging.basicConfig(
        filename=path,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        force=True,
    )
    return path


def write_summary(name: str, summary: dict, summary_dir: Path = SUMMARY_DIR) -> Path:
    summary_dir.mkdir(parents=True, exist_ok=True)
    path = summary_dir / f"{name}.json"
    path.write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    return path
