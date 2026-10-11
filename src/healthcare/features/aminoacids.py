import math

ORDER = "ARNDCQEGHILKMFPSTWYV"

BLOSUM62_ROWS = [
    [4],
    [-1, 5],
    [-2, 0, 6],
    [-2, -2, 1, 6],
    [0, -3, -3, -3, 9],
    [-1, 1, 0, 0, -3, 5],
    [-1, 0, 0, 2, -4, 2, 5],
    [0, -2, 0, -1, -3, -2, -2, 6],
    [-2, 0, 1, -1, -3, 0, 0, -2, 8],
    [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4],
    [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4],
    [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5],
    [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5],
    [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6],
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7],
    [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4],
    [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5],
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11],
    [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7],
    [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4],
]
BLOSUM62 = {
    (a, b): value
    for i, row in enumerate(BLOSUM62_ROWS)
    for j, value in enumerate(row)
    for a, b in ((ORDER[i], ORDER[j]), (ORDER[j], ORDER[i]))
}

GRANTHAM_PROPERTIES = {
    "S": (1.42, 9.2, 32.0),
    "R": (0.65, 10.5, 124.0),
    "L": (0.0, 4.9, 111.0),
    "P": (0.39, 8.0, 32.5),
    "T": (0.71, 8.6, 61.0),
    "A": (0.0, 8.1, 31.0),
    "V": (0.0, 5.9, 84.0),
    "G": (0.74, 9.0, 3.0),
    "I": (0.0, 5.2, 111.0),
    "F": (0.0, 5.2, 132.0),
    "Y": (0.20, 6.2, 136.0),
    "C": (2.75, 5.5, 55.0),
    "H": (0.58, 10.4, 96.0),
    "Q": (0.89, 10.5, 85.0),
    "N": (1.33, 11.6, 56.0),
    "K": (0.33, 11.3, 119.0),
    "D": (1.38, 13.0, 54.0),
    "E": (0.92, 12.3, 83.0),
    "M": (0.0, 5.7, 105.0),
    "W": (0.13, 5.4, 170.0),
}
KYTE_DOOLITTLE = dict(
    zip(
        ORDER,
        [
            1.8,
            -4.5,
            -3.5,
            -3.5,
            2.5,
            -3.5,
            -3.5,
            -0.4,
            -3.2,
            4.5,
            3.8,
            -3.9,
            1.9,
            2.8,
            -1.6,
            -0.8,
            -0.7,
            -0.9,
            -1.3,
            4.2,
        ],
        strict=True,
    )
)
CHARGE = {aa: 0 for aa in ORDER} | {"D": -1, "E": -1, "K": 1, "R": 1}


def grantham(a: str, b: str) -> float | None:
    if a not in GRANTHAM_PROPERTIES or b not in GRANTHAM_PROPERTIES:
        return None
    (c1, p1, v1), (c2, p2, v2) = GRANTHAM_PROPERTIES[a], GRANTHAM_PROPERTIES[b]
    distance = math.sqrt(
        1.833 * (c1 - c2) ** 2 + 0.1018 * (p1 - p2) ** 2 + 0.000399 * (v1 - v2) ** 2
    )
    return round(50.723 * distance, 1)


def substitution_features(a: str, b: str) -> dict:
    known = a in GRANTHAM_PROPERTIES and b in GRANTHAM_PROPERTIES
    return {
        "f_grantham": grantham(a, b),
        "f_blosum62": BLOSUM62.get((a, b)),
        "f_delta_hydropathy": KYTE_DOOLITTLE[b] - KYTE_DOOLITTLE[a] if known else None,
        "f_delta_volume": GRANTHAM_PROPERTIES[b][2] - GRANTHAM_PROPERTIES[a][2] if known else None,
        "f_delta_polarity": GRANTHAM_PROPERTIES[b][1] - GRANTHAM_PROPERTIES[a][1]
        if known
        else None,
        "f_delta_charge": CHARGE[b] - CHARGE[a] if known else None,
    }
