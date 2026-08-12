#!/usr/bin/env python3
"""Generate ``qcdf_wide.f`` from ``qcdf.f`` by raising the colour-array caps.

``qcdf.f`` is left untouched so it keeps reproducing the 1990-era output
byte-for-byte.  This script produces a variant whose colour-contraction arrays
are large enough for high Fock sectors -- see ``docs/fortran-color-overflow.md``.

The fix is purely dimensional: no statement is added, removed or reordered, and
no algorithm changes.  A colour matrix element between two ``L``-parton states
spans

    LNG = LRT + NOPS + LLT  =  L + 4 + L  =  2L + 4          (qcdf.f:5032)

colour indices, so every array indexed by that running length must hold at
least ``2*MXNP + 4 = 54`` entries.  ``qcdf.f`` sizes them 25 (and ``NCR`` only
13), which is why runs reaching ``L >= 11`` corrupt silently.

Two other arrays are dimensioned 25 for unrelated reasons and are deliberately
NOT touched:

    MSTATE(27608,25)   second index is the parton slot (MXNP), not LNG
    IBRPM(9523,25)     second index is the colour count N

The total extra memory is about 3 MB.
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

MXLNG_OLD, MXLNG_NEW = 25, 54          # 2*MXNP + 4

# (old text, new text, expected occurrences, why)
SUBSTITUTIONS = [
    ("INTEGER IDEL0(12552,25),IDELT(12552,25)",
     "INTEGER IDEL0(12552,54),IDELT(12552,54)", 1,
     "delta-function index tables, indexed 1..LNG"),

    ("INTEGER IDEL0(12552,25)",
     "INTEGER IDEL0(12552,54)", 1,
     "same table, declared again in NWTERM"),

    ("INTEGER IDELT(12552,25)",
     "INTEGER IDELT(12552,54)", 1,
     "same table, declared again in CLSUM"),

    ("COMMON/CLR/RESLT(12552),IDELT(12552,25)",
     "COMMON/CLR/RESLT(12552),IDELT(12552,54)", 2,
     "inline COMMON declarations in CNTRCT and BREDCE"),

    ("      COMMON/BAR/LNKB(25),ICYCL(25)",
     "      COMMON/BAR/LNKB(54),ICYCL(54)", 3,
     "baryon link and cycle arrays, indexed 1..LNG"),

    ("      INTEGER MX(4,25)",
     "      INTEGER MX(4,54)", 2,
     "concatenated (right | operators | left) state, length LNG"),

    ("      INTEGER MY(4,25)",
     "      INTEGER MY(4,54)", 1,
     "working copy of MX"),

    ("      INTEGER NCR(2,13)",
     "      INTEGER NCR(2,54)", 1,
     "quark/antiquark index lists; 13 was the tightest cap in the code"),

    ("      INTEGER NPERM(8,25)",
     "      INTEGER NPERM(54,54)", 1,
     "permutation groups; second index runs to 2*count+1"),

    ("      MXLNG=25",
     "      MXLNG=54", 2,
     "the loop bound itself"),
]

# Guard: these must NOT change.
UNTOUCHED = [
    ("COMMON/STATE/MSTATE(27608,25),MSTINF(6902,8)", "parton slots, not LNG"),
    ("COMMON/DLPM/IBRPM(9523,25),NPRMS", "colour count N, not LNG"),
]


def widen(text: str) -> tuple[str, list[str]]:
    log = []
    for old, new, expected, why in SUBSTITUTIONS:
        found = text.count(old)
        if found != expected:
            raise SystemExit(
                f"qcdf.f does not look as expected: found {found} occurrence(s) "
                f"of {old!r}, expected {expected}. Refusing to patch."
            )
        text = text.replace(old, new)
        log.append(f"  {found}x  {old.strip()}  ->  {new.strip()}    ({why})")

    for pattern, why in UNTOUCHED:
        if pattern not in text:
            raise SystemExit(f"expected to leave {pattern!r} untouched ({why})")

    return text, log


def main():
    src = HERE / "qcdf.f"
    dst = HERE / "qcdf_wide.f"

    text = src.read_text()
    widened, log = widen(text)

    header = (
        "c     **** GENERATED FILE -- DO NOT EDIT ****\n"
        "c     **** produced by fortran/widen.py from qcdf.f ****\n"
        "c     **** colour-array caps raised 25 -> 54 (= 2*MXNP+4), NCR 13 -> 54 ****\n"
        "c     **** see docs/fortran-color-overflow.md ****\n"
    )
    dst.write_text(header + widened)

    print(f"wrote {dst}")
    for line in log:
        print(line)
    print(f"\nline count unchanged: "
          f"{len(text.splitlines())} -> {len(widened.splitlines())}")


if __name__ == "__main__":
    sys.exit(main())
