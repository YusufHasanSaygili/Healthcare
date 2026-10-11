import hashlib
import logging
import shutil
import time
import urllib.request
from pathlib import Path

log = logging.getLogger(__name__)

CHUNK_SIZE = 1024 * 1024


class ChecksumMismatch(Exception):
    pass


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(CHUNK_SIZE):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, dest: Path, sha256: str | None = None, force: bool = False) -> dict:
    started = time.monotonic()
    if dest.exists() and not force:
        log.info("skip %s (exists)", dest)
        return {"status": "skipped", "path": str(dest), "bytes": dest.stat().st_size}

    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    if force:
        part.unlink(missing_ok=True)
    offset = part.stat().st_size if part.exists() else 0

    request = urllib.request.Request(url, headers={"Range": f"bytes={offset}-"} if offset else {})
    with urllib.request.urlopen(request, timeout=60) as response:
        resumed = offset > 0 and response.status == 206
        if offset and not resumed:
            log.info("server ignored range for %s, restarting", url)
        with part.open("ab" if resumed else "wb") as handle:
            shutil.copyfileobj(response, handle, CHUNK_SIZE)

    actual = sha256_of(part)
    if sha256 and actual != sha256:
        part.unlink()
        raise ChecksumMismatch(f"{dest.name}: expected {sha256}, got {actual}")

    part.replace(dest)
    log.info("downloaded %s (%d bytes, resumed=%s)", dest, dest.stat().st_size, resumed)
    return {
        "status": "resumed" if resumed else "downloaded",
        "path": str(dest),
        "bytes": dest.stat().st_size,
        "sha256": actual,
        "verified": sha256 is not None,
        "seconds": round(time.monotonic() - started, 2),
    }
