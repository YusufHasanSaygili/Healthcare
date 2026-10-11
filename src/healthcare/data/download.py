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


def digests_of(path: Path) -> tuple[str, str]:
    sha256, md5 = hashlib.sha256(), hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(CHUNK_SIZE):
            sha256.update(chunk)
            md5.update(chunk)
    return sha256.hexdigest(), md5.hexdigest()


def download(
    url: str, dest: Path, sha256: str | None = None, force: bool = False, md5: str | None = None
) -> dict:
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

    actual, actual_md5 = digests_of(part)
    for expected, got in ((sha256, actual), (md5, actual_md5)):
        if expected and got != expected:
            part.unlink()
            raise ChecksumMismatch(f"{dest.name}: expected {expected}, got {got}")

    part.replace(dest)
    log.info("downloaded %s (%d bytes, resumed=%s)", dest, dest.stat().st_size, resumed)
    return {
        "status": "resumed" if resumed else "downloaded",
        "path": str(dest),
        "bytes": dest.stat().st_size,
        "sha256": actual,
        "verified": bool(sha256 or md5),
        "seconds": round(time.monotonic() - started, 2),
    }
