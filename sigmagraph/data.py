"""Small test data.

MULTIPERFECT lists the first 26 multiperfect numbers after 1 (the start of
OEIS A007691).  Each entry is checked by the test suite to satisfy
n | sigma(n).  For the full list of known multiperfect numbers, download a
table such as Achim Flammenkamp's "Multiply Perfect Numbers" page and load
it with :func:`load_multiperfect`.
"""

from __future__ import annotations

MULTIPERFECT = [
    6, 28, 120, 496, 672, 8128, 30240, 32760, 523776, 2178540,
    23569920, 33550336, 45532800, 142990848, 459818240, 1379454720,
    1476304896, 8589869056, 14182439040, 31998395520, 43861478400,
    51001180160, 66433720320, 137438691328, 153003540480, 403031236608,
]

# Descartes' spoof odd perfect number: 22021 = 19**2 * 61 is treated as prime.
DESCARTES_SPOOF = {3: 2, 7: 2, 11: 2, 13: 2, 22021: 1}


def load_multiperfect(path: str) -> list[int]:
    """Read one integer per line (blank lines and '#' comments ignored)."""
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if line:
                out.append(int(line.split()[0]))
    return out
