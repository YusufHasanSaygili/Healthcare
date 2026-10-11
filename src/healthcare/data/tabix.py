import gzip
import struct
import urllib.request
import zlib
from dataclasses import dataclass
from pathlib import Path

LINEAR_SHIFT = 14
MAX_BLOCK = 65536
MERGE_GAP = 65536


@dataclass(frozen=True)
class Chunk:
    start: int
    end: int

    @property
    def first_block(self) -> int:
        return self.start >> 16

    @property
    def last_block(self) -> int:
        return self.end >> 16


@dataclass
class ReferenceIndex:
    bins: dict[int, list[Chunk]]
    linear: list[int]


def read_index(path: Path) -> dict[str, ReferenceIndex]:
    data = gzip.decompress(path.read_bytes())
    if data[:4] != b"TBI\x01":
        raise ValueError(f"not a tabix index: {path}")
    n_ref, *_, names_length = struct.unpack_from("<8i", data, 4)
    offset = 36
    names = data[offset : offset + names_length].split(b"\x00")[:n_ref]
    offset += names_length
    references = {}
    for name in names:
        (n_bin,) = struct.unpack_from("<i", data, offset)
        offset += 4
        bins = {}
        for _ in range(n_bin):
            bin_id, n_chunk = struct.unpack_from("<Ii", data, offset)
            offset += 8
            pairs = struct.unpack_from(f"<{2 * n_chunk}Q", data, offset)
            offset += 16 * n_chunk
            bins[bin_id] = [Chunk(pairs[i], pairs[i + 1]) for i in range(0, len(pairs), 2)]
        (n_intv,) = struct.unpack_from("<i", data, offset)
        offset += 4
        linear = list(struct.unpack_from(f"<{n_intv}Q", data, offset))
        offset += 8 * n_intv
        references[name.decode()] = ReferenceIndex(bins, linear)
    return references


def region_bins(begin: int, end: int) -> list[int]:
    end -= 1
    bins = [0]
    for shift, base in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        bins.extend(range(base + (begin >> shift), base + (end >> shift) + 1))
    return bins


def region_chunks(index: ReferenceIndex, begin: int, end: int) -> list[Chunk]:
    window = begin >> LINEAR_SHIFT
    minimum = index.linear[window] if window < len(index.linear) else 0
    return [
        chunk
        for bin_id in region_bins(begin, end)
        for chunk in index.bins.get(bin_id, [])
        if chunk.end > minimum
    ]


def merge_chunks(chunks: list[Chunk]) -> list[Chunk]:
    merged: list[Chunk] = []
    for chunk in sorted(chunks, key=lambda c: c.start):
        if merged and chunk.first_block <= merged[-1].last_block + MERGE_GAP:
            merged[-1] = Chunk(merged[-1].start, max(merged[-1].end, chunk.end))
        else:
            merged.append(chunk)
    return merged


def byte_range(chunk: Chunk) -> tuple[int, int]:
    return chunk.first_block, chunk.last_block + MAX_BLOCK - 1


def inflate_blocks(data: bytes) -> bytes:
    output = bytearray()
    while data:
        decompressor = zlib.decompressobj(31)
        try:
            block = decompressor.decompress(data)
        except zlib.error:
            break
        if not decompressor.eof:
            break
        output += block
        data = decompressor.unused_data
    return bytes(output)


def chunk_lines(data: bytes, chunk: Chunk) -> list[str]:
    text = inflate_blocks(data)[chunk.start & 0xFFFF :]
    lines = text.split(b"\n")
    return [line.decode() for line in lines[:-1] if line and not line.startswith(b"#")]


def fetch(url: str, chunk: Chunk, timeout: int = 60) -> bytes:
    first, last = byte_range(chunk)
    request = urllib.request.Request(url, headers={"Range": f"bytes={first}-{last}"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        if response.status != 206:
            raise OSError(f"range request not honored for {url}")
        return response.read()
