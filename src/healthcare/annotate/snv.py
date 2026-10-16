import re
from collections import defaultdict
from dataclasses import dataclass

from healthcare.annotate.sequence import CODON_TABLE, reverse_complement

THREE_LETTER = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "U": "Sec",
    "*": "Ter",
}
CODING_CONSEQUENCES = {"missense", "synonymous", "stop_gained", "stop_lost", "start_lost"}
BIN_SIZE = 100_000
SNV_PATTERN = re.compile(
    r"^(?:chr)?([0-9]{1,2}|X|Y)[:-](\d+)\s*([ACGT])\s*>\s*([ACGT])$", re.IGNORECASE
)
VALID_CHROMS = {*map(str, range(1, 23)), "X", "Y"}
ONE_LETTER = {three.upper(): one for one, three in THREE_LETTER.items()}
PROTEIN_PATTERN = re.compile(
    r"^([A-Za-z0-9-]+)\s+(?:p\.)?([A-Za-z]{3}|[A-Za-z*])(\d+)([A-Za-z]{3}|[A-Za-z*])$"
)


@dataclass(frozen=True)
class Transcript:
    transcript_id: str
    gene_name: str
    refseq_transcript: str | None
    chrom: str
    strand: str
    exons: tuple[tuple[int, int], ...]
    cds: tuple[tuple[int, int], ...]
    cds_sequence: str
    protein: str

    @property
    def start(self) -> int:
        return self.exons[0][0]

    @property
    def end(self) -> int:
        return self.exons[-1][1]


@dataclass(frozen=True)
class Annotation:
    transcript_id: str
    gene_name: str
    refseq_transcript: str | None
    consequence: str
    cds_position: int | None = None
    protein_position: int | None = None
    codon_ref: str | None = None
    codon_alt: str | None = None
    aa_ref: str | None = None
    aa_alt: str | None = None
    hgvs_p: str | None = None


def parse_snv(text: str) -> tuple[str, int, str, str]:
    match = SNV_PATTERN.match(text.strip())
    if not match:
        raise ValueError("expected an SNV like 17:43045712 G>A")
    chrom, pos, ref, alt = match.groups()
    chrom, ref, alt = chrom.upper(), ref.upper(), alt.upper()
    if chrom not in VALID_CHROMS or int(pos) < 1:
        raise ValueError(f"chromosome must be 1-22, X or Y and position >= 1: {text}")
    if ref == alt:
        raise ValueError(f"REF and ALT are identical: {text}")
    return chrom, int(pos), ref, alt


def residue(token: str) -> str:
    one = ONE_LETTER.get(token.upper()) if len(token) == 3 else token.upper()
    if one not in THREE_LETTER:
        raise ValueError(f"unknown amino acid: {token}")
    return one


def parse_protein_change(text: str) -> tuple[str, str, int, str]:
    match = PROTEIN_PATTERN.match(text.strip())
    if not match:
        raise ValueError("expected a protein change like BRCA1 p.Arg1699Trp")
    gene, ref, position, alt = match.groups()
    return gene.upper(), residue(ref), int(position), residue(alt)


def genomic_position(transcript: Transcript, offset: int) -> int:
    intervals = transcript.cds if transcript.strand == "+" else transcript.cds[::-1]
    for start, end in intervals:
        length = end - start + 1
        if offset < length:
            return start + offset if transcript.strand == "+" else end - offset
        offset -= length
    raise ValueError("CDS offset is beyond the coding sequence")


def protein_change_snvs(
    transcript: Transcript, aa_ref: str, position: int, aa_alt: str
) -> list[tuple[str, int, str, str]]:
    if not 1 <= position <= len(transcript.protein):
        raise ValueError(
            f"{transcript.gene_name} has {len(transcript.protein)} amino acids, not {position}"
        )
    actual = transcript.protein[position - 1]
    if actual != aa_ref:
        raise ValueError(
            f"{transcript.gene_name} residue {position} is {THREE_LETTER[actual]}, "
            f"not {THREE_LETTER[aa_ref]}"
        )
    variants = []
    for within in range(3):
        offset = 3 * (position - 1) + within
        pos = genomic_position(transcript, offset)
        base = transcript.cds_sequence[offset]
        ref = base if transcript.strand == "+" else reverse_complement(base)
        for alt in "ACGT":
            if alt != ref and annotate_snv(transcript, pos, ref, alt).aa_alt == aa_alt:
                variants.append((transcript.chrom, pos, ref, alt))
    return sorted(variants, key=lambda v: (v[1], v[3]))


def cds_offset(transcript: Transcript, pos: int) -> int | None:
    offset = 0
    intervals = transcript.cds if transcript.strand == "+" else transcript.cds[::-1]
    for start, end in intervals:
        if start <= pos <= end:
            return offset + (pos - start if transcript.strand == "+" else end - pos)
        offset += end - start + 1
    return None


def splice_site(transcript: Transcript, pos: int) -> bool:
    for (_, left_end), (right_start, _) in zip(
        transcript.exons, transcript.exons[1:], strict=False
    ):
        if pos in (left_end + 1, left_end + 2, right_start - 2, right_start - 1):
            return True
    return False


def in_exon(transcript: Transcript, pos: int) -> bool:
    return any(start <= pos <= end for start, end in transcript.exons)


def protein_hgvs(aa_ref: str, position: int, aa_alt: str, consequence: str) -> str:
    ref = THREE_LETTER[aa_ref]
    if consequence == "synonymous":
        return f"p.{ref}{position}="
    if consequence == "start_lost":
        return "p.Met1?"
    if consequence == "stop_lost":
        return f"p.{ref}{position}{THREE_LETTER[aa_alt]}extTer?"
    return f"p.{ref}{position}{THREE_LETTER[aa_alt]}"


def coding_annotation(transcript: Transcript, offset: int, alt: str) -> Annotation:
    codon_index = offset // 3
    codon_ref = transcript.cds_sequence[codon_index * 3 : codon_index * 3 + 3]
    alt_base = alt if transcript.strand == "+" else reverse_complement(alt)
    within = offset % 3
    codon_alt = codon_ref[:within] + alt_base + codon_ref[within + 1 :]
    last_codon = len(transcript.cds_sequence) // 3 - 1
    aa_ref = "*" if codon_index == last_codon else transcript.protein[codon_index]
    aa_alt = CODON_TABLE.get(codon_alt, "X")

    if codon_index == 0:
        consequence = "start_lost"
    elif codon_index == last_codon:
        consequence = "synonymous" if aa_alt == "*" else "stop_lost"
    elif aa_alt == "*" and aa_ref != "U":
        consequence = "stop_gained"
    elif aa_alt == aa_ref:
        consequence = "synonymous"
    else:
        consequence = "missense"

    position = codon_index + 1
    return Annotation(
        transcript_id=transcript.transcript_id,
        gene_name=transcript.gene_name,
        refseq_transcript=transcript.refseq_transcript,
        consequence=consequence,
        cds_position=offset + 1,
        protein_position=position,
        codon_ref=codon_ref,
        codon_alt=codon_alt,
        aa_ref=aa_ref,
        aa_alt=aa_alt,
        hgvs_p=protein_hgvs(aa_ref, position, aa_alt, consequence),
    )


def annotate_snv(transcript: Transcript, pos: int, ref: str, alt: str) -> Annotation:
    if len(ref) != 1 or len(alt) != 1 or ref == alt:
        raise ValueError(f"not an SNV: {ref}>{alt}")
    offset = cds_offset(transcript, pos)
    if offset is not None:
        ref_base = ref if transcript.strand == "+" else reverse_complement(ref)
        if transcript.cds_sequence[offset] != ref_base:
            return Annotation(
                transcript.transcript_id,
                transcript.gene_name,
                transcript.refseq_transcript,
                "ref_mismatch",
            )
        return coding_annotation(transcript, offset, alt)
    if splice_site(transcript, pos):
        consequence = "canonical_splice"
    elif in_exon(transcript, pos):
        consequence = "utr"
    else:
        consequence = "intronic"
    return Annotation(
        transcript.transcript_id, transcript.gene_name, transcript.refseq_transcript, consequence
    )


class TranscriptIndex:
    def __init__(self, transcripts: list[Transcript]):
        self.bins: dict[tuple[str, int], list[Transcript]] = defaultdict(list)
        for transcript in transcripts:
            for bin_id in range(transcript.start // BIN_SIZE, transcript.end // BIN_SIZE + 1):
                self.bins[(transcript.chrom, bin_id)].append(transcript)

    def overlapping(self, chrom: str, pos: int) -> list[Transcript]:
        return [t for t in self.bins.get((chrom, pos // BIN_SIZE), []) if t.start <= pos <= t.end]

    def annotate(self, chrom: str, pos: int, ref: str, alt: str) -> list[Annotation]:
        return [annotate_snv(t, pos, ref, alt) for t in self.overlapping(chrom, pos)]
