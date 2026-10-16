from healthcare.annotate.indel import Fetch, edit
from healthcare.annotate.sequence import reverse_complement
from healthcare.annotate.snv import Transcript, cds_offset

STOP_CODONS = {"TAA", "TAG", "TGA"}
NMD_DISTANCE = 55

Edit = tuple[int, int, str]


def in_transcript_order(transcript: Transcript) -> tuple[tuple[int, int], ...]:
    return transcript.exons if transcript.strand == "+" else transcript.exons[::-1]


def mrna_offset(transcript: Transcript, pos: int) -> int | None:
    offset = 0
    for start, end in in_transcript_order(transcript):
        if start <= pos <= end:
            return offset + (pos - start if transcript.strand == "+" else end - pos)
        offset += end - start + 1
    return None


def last_junction(transcript: Transcript) -> int | None:
    exons = in_transcript_order(transcript)
    if len(exons) < 2:
        return None
    return sum(end - start + 1 for start, end in exons[:-1])


def cds_start(transcript: Transcript) -> int:
    first_base = transcript.cds[0][0] if transcript.strand == "+" else transcript.cds[-1][1]
    offset = mrna_offset(transcript, first_base)
    if offset is None:
        raise ValueError(f"{transcript.transcript_id}: CDS start is outside its exons")
    return offset


def utr3_sequence(transcript: Transcript, fetch: Fetch) -> str:
    plus = transcript.strand == "+"
    stop_base = transcript.cds[-1][1] if plus else transcript.cds[0][0]
    pieces = []
    for start, end in in_transcript_order(transcript):
        low, high = (max(start, stop_base + 1), end) if plus else (start, min(end, stop_base - 1))
        if low <= high:
            bases = fetch(low, high)
            pieces.append(bases if plus else reverse_complement(bases))
    return "".join(pieces)


def cds_edit(transcript: Transcript, pos: int, ref: str, alt: str) -> Edit | None:
    start, deleted, inserted = edit(pos, ref, alt)
    oriented = inserted if transcript.strand == "+" else reverse_complement(inserted)
    if deleted:
        first, last = (
            cds_offset(transcript, start),
            cds_offset(transcript, start + len(deleted) - 1),
        )
        if first is None or last is None or abs(first - last) != len(deleted) - 1:
            return None
        return min(first, last), max(first, last), oriented
    left, right = cds_offset(transcript, start - 1), cds_offset(transcript, start)
    if left is None or right is None or abs(left - right) != 1:
        return None
    point = max(left, right)
    return point, point - 1, oriented


def next_stop(sequence: str, start: int) -> int | None:
    for index in range(start, len(sequence) - 2, 3):
        if sequence[index : index + 3] in STOP_CODONS:
            return index
    return None


def truncation_features(transcript: Transcript, change: Edit, utr3: str) -> dict[str, float | None]:
    first, last, inserted = change
    mutant = transcript.cds_sequence[:first] + inserted + transcript.cds_sequence[last + 1 :] + utr3
    length = len(transcript.protein)
    stop = next_stop(mutant, first - first % 3)
    retained = (first // 3) / length
    if stop is None:
        return {"f_no_stop": 1.0, "f_protein_retained": retained}
    shift = len(inserted) - (last - first + 1)
    reference_stop = stop if stop < first else stop - shift
    junction = last_junction(transcript)
    distance = None if junction is None else junction - (cds_start(transcript) + reference_stop)
    escapes = distance is None or distance <= NMD_DISTANCE
    return {
        "f_no_stop": 0.0,
        "f_protein_retained": retained,
        "f_new_length_ratio": (stop // 3) / length,
        "f_novel_residues": float(stop // 3 - first // 3),
        "f_ptc_junction_distance": None if distance is None else float(distance),
        "f_ptc_in_last_exon": float(distance is None or distance <= 0),
        "f_nmd_escape": float(escapes),
    }


def start_loss_features(transcript: Transcript) -> dict[str, float | None]:
    index = transcript.protein.find("M", 1)
    return {"f_next_met_fraction": None if index < 0 else index / len(transcript.protein)}


def cds_bases(transcript: Transcript, exon: tuple[int, int]) -> int:
    return sum(max(0, min(exon[1], end) - max(exon[0], start) + 1) for start, end in transcript.cds)


def splice_features(transcript: Transcript, start: int, end: int) -> dict[str, float | None]:
    exons = in_transcript_order(transcript)
    plus = transcript.strand == "+"
    for index, (upstream, downstream) in enumerate(zip(exons, exons[1:], strict=False)):
        boundary_out = upstream[1] if plus else upstream[0]
        boundary_in = downstream[0] if plus else downstream[1]
        step = 1 if plus else -1
        sites = [(True, index, offset, boundary_out + step * offset) for offset in (1, 2)] + [
            (False, index + 1, offset, boundary_in - step * offset) for offset in (1, 2)
        ]
        hits = [site for site in sites if start <= site[3] <= end]
        if hits:
            donor, skipped, offset, _ = min(hits, key=lambda site: site[2])
            coding = cds_bases(transcript, exons[skipped])
            return {
                "f_splice_donor": float(donor),
                "f_splice_intron_position": float(offset),
                "f_skipped_exon_in_frame": None if coding == 0 else float(coding % 3 == 0),
                "f_skipped_exon_cds_fraction": coding / len(transcript.cds_sequence),
                "f_skipped_exon_rank": skipped / (len(exons) - 1),
            }
    return {}


def boundary_distance(transcript: Transcript, pos: int) -> float:
    return float(min(min(abs(pos - start), abs(pos - end)) for start, end in transcript.exons))
