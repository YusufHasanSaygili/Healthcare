import tomllib
from dataclasses import dataclass
from pathlib import Path

from healthcare.paths import CONFIG_DIR, DATA_DIR

SOURCES_FILE = CONFIG_DIR / "sources.toml"


@dataclass(frozen=True)
class Source:
    name: str
    url: str
    version: str
    sha256: str | None = None
    md5: str | None = None

    @property
    def path(self) -> Path:
        return DATA_DIR / "raw" / self.name / self.url.rsplit("/", 1)[-1]


def load_sources(path: Path = SOURCES_FILE) -> dict[str, Source]:
    entries = tomllib.loads(path.read_text(encoding="utf-8")).get("sources", {})
    return {name: Source(name=name, **fields) for name, fields in entries.items()}
