import re
from collections.abc import Callable

from healthcare.annotate.sequence import reverse_complement, translate
from healthcare.annotate.snv import (
    THREE_LETTER,
    VALID_CHROMS,
    Annotation,
    Transcript,
    TranscriptIndex,
    cds_offset,
    in_exon,
    splice_site,
)

VARIANT_PATTERN = re.compile(
    r"^(?:chr)?([0-9]{1,2}|X|Y)[:-](\d+)\s*([ACGT]+)\s*>\s*([ACGT]+)$", re.IGNORECASE
)

Fetch = Callable[[int, int], str]


def parse_variant(text: str) -> tuple[str, int, str, str]:
    match = VARIANT_PATTERN.match(text.strip())
    if not match:
        raise ValueError("expected a variant like 17:43045712 G>A or 17:43045712 GA>G")
    chrom, pos, ref, alt = match.groups()
    chrom, ref, alt = chrom.upper(), ref.upper(), alt.upper()
    if chrom not in VALID_CHROMS or int(pos) < 1:
        raise ValueError(f"chromosome must be 1-22, X or Y and position >= 1: {text}")
    if ref == alt:
        raise ValueError(f"REF and ALT are identical: {text}")
    return chrom, int(pos), ref, alt


def normalize(pos: int, ref: str, alt: str, fetch: Fetch) -> tuple[int, str, str]:
    while True:
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        elif (not ref or not alt) and pos > 1:
            pos -= 1
            base = fetch(pos, pos)
            ref, alt = base + ref, base + alt
        else:
            break
    if not ref or not alt:
        raise ValueError("the variant cannot be anchored at the start of the chromosome")
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def edit(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    shared = 0
    while shared < min(len(ref), len(alt)) and ref[shared] == alt[shared]:
        shared += 1
    ref, alt, pos = ref[shared:], alt[shared:], pos + shared
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    return pos, ref, alt


def three(residues: str) -> str:
    return "".join(THREE_LETTER.get(residue, "Xaa") for residue in residues)


def span(residues: str, start: int) -> str:
    first = f"{THREE_LETTER.get(residues[0], 'Xaa')}{start}"
    if len(residues) == 1:
        return first
    return f"{first}_{THREE_LETTER.get(residues[-1], 'Xaa')}{start + len(residues) - 1}"


def common_prefix(a: str, b: str) -> int:
    size = 0
    while size < min(len(a), len(b)) and a[size] == b[size]:
        size += 1
    return size


def common_suffix(a: str, b: str, floor: int) -> int:
    size = 0
    while size < min(len(a), len(b)) - floor and a[-1 - size] == b[-1 - size]:
        size += 1
    return size


def stop_retained(ref_tail: str, first: int, codon: int) -> tuple[str, str, int, str] | None:
    stop = ref_tail.find("*")
    if 0 <= stop < first:
        return "synonymous", "", codon + stop + 1, f"p.Ter{codon + stop + 1}="
    return None


def frameshift(ref_tail: str, mutant_tail: str, codon: int) -> tuple[str, str, int, str]:
    stop = mutant_tail.find("*")
    mutant_tail = mutant_tail[: stop + 1] if stop >= 0 else mutant_tail
    first = common_prefix(ref_tail, mutant_tail)
    if retained := stop_retained(ref_tail, first, codon):
        return retained
    position = codon + first + 1
    ref_residue = ref_tail[first] if first < len(ref_tail) else "*"
    if ref_residue == "*":
        new = f"{three(mutant_tail[first])}extTer?" if first < len(mutant_tail) else "?"
        return "stop_lost", "", position, f"p.Ter{position}{new}"
    if first < len(mutant_tail) and mutant_tail[first] == "*":
        return "frameshift", "", position, f"p.{three(ref_residue)}{position}Ter"
    return "frameshift", "", position, f"p.{three(ref_residue)}{position}fs"


def in_frame(ref_tail: str, mutant_tail: str, codon: int) -> tuple[str, str, int, str | None]:
    first = common_prefix(ref_tail, mutant_tail)
    if retained := stop_retained(ref_tail, first, codon):
        return retained
    last = common_suffix(ref_tail, mutant_tail, first)
    removed = ref_tail[first : len(ref_tail) - last]
    added = mutant_tail[first : len(mutant_tail) - last]
    position = codon + first + 1
    if "*" in removed:
        new = f"{THREE_LETTER.get(added[0], 'Xaa')}extTer?" if added else "?"
        return "stop_lost", removed, position, f"p.Ter{codon + ref_tail.index('*') + 1}{new}"
    if "*" in added:
        added = added[: added.index("*") + 1]
    consequence = "stop_gained" if "*" in added else None
    if not removed and not added:
        return "synonymous", "", position, None
    if not removed:
        before = ref_tail[first - len(added) : first] if first >= len(added) else ""
        if before == added and "*" not in added:
            return "inframe_insertion", "", position, f"p.{span(added, position - len(added))}dup"
        flank = f"{three(ref_tail[first - 1])}{position - 1}_{three(ref_tail[first])}{position}"
        return consequence or "inframe_insertion", "", position, f"p.{flank}ins{three(added)}"
    if not added:
        return "inframe_deletion", removed, position, f"p.{span(removed, position)}del"
    if len(removed) == len(added) == 1:
        hgvs = f"p.{three(removed)}{position}{three(added)}"
        return consequence or "inframe_indel", removed, position, hgvs
    return (
        consequence or "inframe_indel",
        removed,
        position,
        f"p.{span(removed, position)}delins{three(added)}",
    )


def selenocysteine(
    protein: str, codon: int, ref_tail: str, mutant_tail: str, before: int, offset: int | None
) -> tuple[str, str]:
    ref, mutant = list(ref_tail), list(mutant_tail)
    for index in range(len(ref)):
        if codon + index >= len(protein) or protein[codon + index] != "U":
            continue
        ref[index] = "U"
        target = index if index < before else None if offset is None else index + offset
        if target is not None and 0 <= target < len(mutant) and mutant[target] == "*":
            mutant[target] = "U"
    return "".join(ref), "".join(mutant)


def coding_effect(transcript: Transcript, first: int, last: int, inserted: str) -> Annotation:
    cds = transcript.cds_sequence
    mutant = cds[:first] + inserted + cds[last + 1 :]
    if first < 3:
        consequence, removed, position, hgvs = "start_lost", "M", 1, "p.Met1?"
    else:
        codon = max(first // 3 - 1 - len(inserted) // 3, 0)
        change = len(inserted) - (last - first + 1)
        shift = change % 3
        ref_tail, mutant_tail = selenocysteine(
            transcript.protein,
            codon,
            translate(cds[codon * 3 :]),
            translate(mutant[codon * 3 :]),
            first // 3 - codon,
            None if shift else change // 3,
        )
        consequence, removed, position, hgvs = (
            frameshift(ref_tail, mutant_tail, codon)
            if shift
            else in_frame(ref_tail, mutant_tail, codon)
        )
    return Annotation(
        transcript_id=transcript.transcript_id,
        gene_name=transcript.gene_name,
        refseq_transcript=transcript.refseq_transcript,
        consequence=consequence,
        cds_position=first + 1,
        protein_position=position,
        aa_ref=removed or None,
        hgvs_p=hgvs,
    )


def flagged(transcript: Transcript, consequence: str) -> Annotation:
    return Annotation(
        transcript.transcript_id, transcript.gene_name, transcript.refseq_transcript, consequence
    )


def splice_overlap(transcript: Transcript, start: int, end: int) -> bool:
    for (_, left_end), (right_start, _) in zip(
        transcript.exons, transcript.exons[1:], strict=False
    ):
        if start <= left_end + 2 and left_end + 1 <= end:
            return True
        if start <= right_start - 1 and right_start - 2 <= end:
            return True
    return False


def deletion_effect(transcript: Transcript, start: int, deleted: str, inserted: str) -> Annotation:
    end = start + len(deleted) - 1
    if splice_overlap(transcript, start, end):
        return flagged(transcript, "canonical_splice")
    ends = cds_offset(transcript, start), cds_offset(transcript, end)
    if None not in ends and abs(ends[0] - ends[1]) == end - start:
        coding = list(ends)
    else:
        offsets = (cds_offset(transcript, pos) for pos in range(start, end + 1))
        coding = [offset for offset in offsets if offset is not None]
    if not coding:
        exonic = any(s <= end and start <= e for s, e in transcript.exons)
        return flagged(transcript, "utr" if exonic else "intronic")
    if len(coding) < end - start + 1 and abs(coding[0] - coding[-1]) != end - start:
        if min(coding) < 3:
            return flagged(transcript, "start_lost")
        if max(coding) >= len(transcript.cds_sequence) - 3:
            return flagged(transcript, "stop_lost")
        return flagged(transcript, "splice_region")
    first, last = min(coding), max(coding)
    plus = transcript.strand == "+"
    reference = deleted if plus else reverse_complement(deleted)
    if transcript.cds_sequence[first : last + 1] != reference:
        return flagged(transcript, "ref_mismatch")
    return coding_effect(
        transcript, first, last, inserted if plus else reverse_complement(inserted)
    )


def insertion_effect(transcript: Transcript, start: int, inserted: str) -> Annotation:
    before, after = start - 1, start
    if splice_site(transcript, before) and splice_site(transcript, after):
        return flagged(transcript, "canonical_splice")
    left, right = cds_offset(transcript, before), cds_offset(transcript, after)
    if left is not None and right is not None and abs(left - right) == 1:
        point = max(left, right)
        plus = transcript.strand == "+"
        return coding_effect(
            transcript, point, point - 1, inserted if plus else reverse_complement(inserted)
        )
    exonic = [in_exon(transcript, before), in_exon(transcript, after)]
    if all(exonic):
        return flagged(transcript, "utr")
    return flagged(transcript, "splice_region" if any(exonic) else "intronic")


def annotate_indel(transcript: Transcript, pos: int, ref: str, alt: str) -> Annotation:
    start, deleted, inserted = edit(pos, ref, alt)
    if len(deleted) == len(inserted):
        return flagged(transcript, "multi_nucleotide")
    if not deleted:
        return insertion_effect(transcript, start, inserted)
    return deletion_effect(transcript, start, deleted, inserted)


def affected_span(pos: int, ref: str) -> tuple[int, int]:
    return pos, pos + max(len(ref), 2) - 1


def annotate_variant(
    index: TranscriptIndex, chrom: str, pos: int, ref: str, alt: str
) -> list[Annotation]:
    if len(ref) == 1 and len(alt) == 1:
        return index.annotate(chrom, pos, ref, alt)
    start, end = affected_span(pos, ref)
    return [annotate_indel(t, pos, ref, alt) for t in index.overlapping_span(chrom, start, end)]
