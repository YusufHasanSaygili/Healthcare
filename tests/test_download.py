import hashlib
import json
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest

from healthcare.data.download import ChecksumMismatch, download
from healthcare.data.sources import load_sources
from healthcare.runlog import write_summary

PAYLOAD = bytes(range(256)) * 4096
PAYLOAD_SHA256 = hashlib.sha256(PAYLOAD).hexdigest()


class Handler(BaseHTTPRequestHandler):
    honor_range = True
    range_requests: list[str] = []

    def do_GET(self):
        requested = self.headers.get("Range")
        if requested:
            Handler.range_requests.append(requested)
        if requested and self.honor_range:
            start = int(requested.removeprefix("bytes=").rstrip("-"))
            body = PAYLOAD[start:]
            self.send_response(206)
        else:
            body = PAYLOAD
            self.send_response(200)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *args):
        pass


@pytest.fixture
def server():
    Handler.honor_range = True
    Handler.range_requests = []
    httpd = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{httpd.server_port}/file.bin"
    httpd.shutdown()


def test_downloads_and_verifies(server, tmp_path):
    dest = tmp_path / "file.bin"
    result = download(server, dest, sha256=PAYLOAD_SHA256)
    assert result["status"] == "downloaded"
    assert result["verified"]
    assert dest.read_bytes() == PAYLOAD


def test_skips_existing_file_unless_forced(server, tmp_path):
    dest = tmp_path / "file.bin"
    dest.write_bytes(b"old")
    assert download(server, dest)["status"] == "skipped"
    assert download(server, dest, force=True)["status"] == "downloaded"
    assert dest.read_bytes() == PAYLOAD


def test_resumes_partial_download(server, tmp_path):
    dest = tmp_path / "file.bin"
    (tmp_path / "file.bin.part").write_bytes(PAYLOAD[:1000])
    result = download(server, dest, sha256=PAYLOAD_SHA256)
    assert result["status"] == "resumed"
    assert Handler.range_requests == ["bytes=1000-"]
    assert dest.read_bytes() == PAYLOAD


def test_restarts_when_server_ignores_range(server, tmp_path):
    Handler.honor_range = False
    dest = tmp_path / "file.bin"
    (tmp_path / "file.bin.part").write_bytes(b"garbage")
    result = download(server, dest, sha256=PAYLOAD_SHA256)
    assert result["status"] == "downloaded"
    assert dest.read_bytes() == PAYLOAD


def test_checksum_mismatch_leaves_no_file(server, tmp_path):
    dest = tmp_path / "file.bin"
    with pytest.raises(ChecksumMismatch):
        download(server, dest, sha256="0" * 64)
    assert not dest.exists()
    assert not (tmp_path / "file.bin.part").exists()


def test_summary_is_json(tmp_path):
    path = write_summary("demo", {"count": 3, "path": tmp_path}, summary_dir=tmp_path)
    assert json.loads(path.read_text())["count"] == 3


def test_load_sources(tmp_path):
    config = tmp_path / "sources.toml"
    config.write_text(
        '[sources.demo]\nurl = "https://example.org/a/b.txt.gz"\nversion = "1"\n', encoding="utf-8"
    )
    source = load_sources(config)["demo"]
    assert source.sha256 is None
    assert source.path.parts[-3:] == ("raw", "demo", "b.txt.gz")
