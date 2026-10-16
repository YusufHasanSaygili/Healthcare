import re

from healthcare.annotate.snv import ONE_LETTER

AA = r"[A-Z][a-z]{2}"
SUBSTITUTION = re.compile(rf"^p\.({AA})(\d+)({AA})(?:extTer\?)?$")
FRAMESHIFT = re.compile(rf"^p\.({AA})(\d+)fs$")
RANGE_EDIT = re.compile(rf"^p\.({AA})(\d+)(?:_({AA})(\d+))?(del|dup|delins((?:{AA})+))$")
INSERTION = re.compile(rf"^p\.({AA})(\d+)_({AA})(\d+)ins((?:{AA})+)$")
EXTENSION = re.compile(r"^p\.Ter(\d+)\?$")
REPEAT = re.compile(r"^p\.(\d+)(?:_(\d+))?([A-Z]+)\[(\d+)\]$")

Effect = tuple[str, int | str]


def residues(three_letter: str) -> str:
    return "".join(ONE_LETTER.get(code.upper(), "?") for code in re.findall(AA, three_letter))


def mutant(reference: str, start: int, end: int, replacement: str) -> Effect:
    changed = reference[: start - 1] + replacement + reference[end:]
    stop = changed.find("*")
    return ("protein", changed[: stop + 1] if stop >= 0 else changed)


def matches(reference: str, position: int, code: str) -> bool:
    return 1 <= position <= len(reference) and reference[position - 1] == residues(code)


def repeat_effect(reference: str, match: re.Match) -> Effect | None:
    start, unit, count = int(match.group(1)), match.group(3), int(match.group(4))
    copies = 0
    while reference[start - 1 + copies * len(unit) :].startswith(unit):
        copies += 1
    if copies == 0:
        return None
    return mutant(reference, start, start - 1 + copies * len(unit), unit * count)


def protein_effect(hgvs: str | None, protein: str) -> Effect | None:
    if hgvs is None:
        return None
    reference = protein + "*"
    if match := EXTENSION.match(hgvs):
        position = int(match.group(1))
        return ("stop_lost", position) if matches(reference, position, "Ter") else None
    if match := REPEAT.match(hgvs):
        return repeat_effect(reference, match)
    if match := FRAMESHIFT.match(hgvs):
        position = int(match.group(2))
        return ("frameshift", position) if matches(reference, position, match.group(1)) else None
    if match := SUBSTITUTION.match(hgvs):
        position = int(match.group(2))
        if not matches(reference, position, match.group(1)):
            return None
        if match.group(1) == "Ter":
            return ("stop_lost", position)
        return mutant(reference, position, position, residues(match.group(3)))
    if match := INSERTION.match(hgvs):
        position = int(match.group(2))
        if not matches(reference, position, match.group(1)):
            return None
        return mutant(reference, position + 1, position, residues(match.group(5)))
    if match := RANGE_EDIT.match(hgvs):
        start = int(match.group(2))
        end = int(match.group(4)) if match.group(4) else start
        if not matches(reference, start, match.group(1)):
            return None
        if match.group(5) == "del":
            return mutant(reference, start, end, "")
        if match.group(5) == "dup":
            return mutant(reference, end + 1, end, reference[start - 1 : end])
        return mutant(reference, start, end, residues(match.group(6)))
    return None


def first_change(effect: Effect, protein: str) -> int:
    if effect[0] != "protein":
        return int(effect[1])
    reference, changed = protein + "*", str(effect[1])
    size = 0
    while size < min(len(reference), len(changed)) and reference[size] == changed[size]:
        size += 1
    return size + 1


def equivalent(ours: str | None, theirs: str, protein: str) -> bool:
    a, b = protein_effect(ours, protein), protein_effect(theirs, protein)
    if a is None or b is None:
        return False
    if a == b:
        return True
    truncating = {a[0], b[0]} == {"frameshift", "protein"} and "Ter" in (
        theirs if b[0] == "protein" else ours
    )
    return truncating and first_change(a, protein) == first_change(b, protein)
