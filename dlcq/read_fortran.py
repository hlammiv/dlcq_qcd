"""Parse the Fortran solver's ``qcdf.out`` / ``qcdf.ham`` into a DLCQResult.

``qcdf.f`` writes a human-readable log, not a data file, so this module carries
the knowledge needed to read it safely.  Four hazards, all of which silently
produce wrong physics if mishandled:

1. **Repeated eigenvalue blocks.**  ``PRTEIG`` is called once per weeding
   iteration with ``IVECT=0``, then once at the end with ``IVECT=IVCTMX``.  A
   K=21 run emits 27 ``COMPUTED EIGENVALUES:`` blocks and only the last is the
   physical spectrum -- the first 26 are *norm-matrix* eigenvalues.  We anchor
   on the final block, the one followed by ``COMPUTED EIGENVECTORS:``.

2. **The IVCTMX=75 eigenvector cap.**  ``PRTEIG`` prints ``N1=NUMSTA`` rows but
   only ``min(IVCTMX, NUMSTA)`` columns.  A run reports 189 eigenvalues and 75
   eigenvectors.  We record both counts so consumers fail loudly.

3. **Ragged ``PRTSPM`` rows.**  The final partial row uses a format with five
   ``(i,j)`` groups but supplies fewer, so Fortran emits a dangling ``" ("``
   before running out of data.  The index regex must tolerate it.

4. **``Z`` is reused.**  ``NUZ(W,Z)`` makes Z the orthonormalizing basis change
   (printed as ``MATRIX Z WHICH GIVES NEW BASIS``); ``DIAG(HNU,W,Z,IERR)`` then
   overwrites it with HNU's eigenvectors (printed as ``COMPUTED EIGENVECTORS``).
   The Fock-basis coefficients are the product: ``c_orig = Z_basis @ z_eig``.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from .dataset import DLCQResult

__all__ = ["read_out", "read_ham", "parse_eigenvalues"]


# Fortran reals: 0.123E+02, -.5D-03, 1.5e+00 ... accept D/d/E/e exponents.
_FLOAT_RE = re.compile(r"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][+-]?\d+)?")
# Index pairs from PRTSPM; tolerates the dangling "(" on the ragged last row.
_PAIR_RE = re.compile(r"\(\s*(\d+)\s*,\s*(\d+)\s*\)")


def _floats(line: str) -> list[float]:
    """Every Fortran real on a line, with D exponents normalized."""
    return [float(m.group(0).replace("D", "E").replace("d", "e"))
            for m in _FLOAT_RE.finditer(line)]


def _find(lines: list[str], needle: str, start: int = 0) -> int:
    for i in range(start, len(lines)):
        if needle in lines[i]:
            return i
    return -1


def _rfind(lines: list[str], needle: str) -> int:
    for i in range(len(lines) - 1, -1, -1):
        if needle in lines[i]:
            return i
    return -1


# ──────────────────────────────────────────────────────────────────────────
# Header
# ──────────────────────────────────────────────────────────────────────────

def _parse_header(lines: list[str]) -> dict:
    """Read the run parameters echoed at the top of qcdf.out."""
    meta: dict = {}

    i = _find(lines, "N,NF,B,K,LAMBDA:")
    if i < 0:
        raise ValueError("no 'N,NF,B,K,LAMBDA:' header found; not a qcdf.out?")
    vals = _floats(lines[i].split(":", 1)[1])
    if len(vals) < 5:
        raise ValueError(f"malformed header line: {lines[i]!r}")
    meta["N"], meta["NF"], meta["B"], meta["K_code"] = (int(v) for v in vals[:4])
    meta["rlamb"] = float(vals[4])

    i = _find(lines, "LPN:")
    meta["LPN"] = int(_floats(lines[i].split(":", 1)[1])[0]) if i >= 0 else 0

    i = _find(lines, "CUTOFF/MASS SQD:")
    meta["cutoff"] = float(_floats(lines[i].split(":", 1)[1])[0]) if i >= 0 else -1.0

    i = _find(lines, "ratio of quark masses")
    if i >= 0:
        meta["rmq"] = np.array(_floats(lines[i + 1])[: meta["NF"]])
    else:
        meta["rmq"] = np.ones(meta["NF"])

    i = _find(lines, "flavor numbers")
    if i >= 0:
        meta["iflv"] = np.array([int(v) for v in _floats(lines[i + 1])][: meta["NF"]])
    else:
        meta["iflv"] = np.zeros(meta["NF"], dtype=int)

    i = _find(lines, "NUMBER OF STATES BEFORE WEEDING")
    meta["numsta_pre"] = int(_floats(lines[i])[0]) if i >= 0 else 0

    i = _find(lines, "NUMBER OF STATES AFTER WEEDING")
    meta["numsta_post"] = int(_floats(lines[i])[0]) if i >= 0 else meta["numsta_pre"]

    return meta


# ──────────────────────────────────────────────────────────────────────────
# Eigenvalues / eigenvectors
# ──────────────────────────────────────────────────────────────────────────

def parse_eigenvalues(lines: list[str], expected: int | None = None) -> np.ndarray:
    """The *physical* spectrum: the last COMPUTED EIGENVALUES block.

    Equivalent in intent to ``fortran/extract_eigenvalues.sh`` but without its
    re-sort, so the solver's own ordering is preserved.
    """
    start = _rfind(lines, "COMPUTED EIGENVALUES")
    if start < 0:
        raise ValueError("no 'COMPUTED EIGENVALUES' block found")

    values: list[float] = []
    for line in lines[start + 1:]:
        if "COMPUTED EIGENVECTORS" in line or "IERR" in line:
            break
        if "HAVE DROPPED" in line or "MDROP" in line or "NUMSTA=" in line:
            break
        values.extend(_floats(line))
        if expected is not None and len(values) >= expected:
            break

    arr = np.array(values[:expected] if expected else values, dtype=float)
    return arr


def _parse_eigenvectors(lines: list[str], n_rows: int) -> np.ndarray | None:
    """The COMPUTED EIGENVECTORS block: column-blocks of <=4, ``' ------'``-separated.

    ``PRTEIG`` emits, for each block of up to 4 columns, ``n_rows`` lines of
    ``4e24.16`` followed by a separator.  Column count is capped at IVCTMX.
    """
    start = _rfind(lines, "COMPUTED EIGENVECTORS")
    if start < 0:
        return None

    blocks: list[np.ndarray] = []
    current: list[list[float]] = []

    for line in lines[start + 1:]:
        if line.strip().startswith("------"):
            if current:
                blocks.append(np.array(current, dtype=float))
                current = []
            continue
        vals = _floats(line)
        if vals:
            current.append(vals)
    if current:
        blocks.append(np.array(current, dtype=float))

    usable = [b for b in blocks if b.ndim == 2 and b.shape[0] == n_rows]
    if not usable:
        return None
    return np.hstack(usable)


# ──────────────────────────────────────────────────────────────────────────
# Sparse matrices (PRTSPM)
# ──────────────────────────────────────────────────────────────────────────

def _parse_sparse(lines: list[str], start: int, dim: int,
                  symmetrize: bool, n_cols: int | None = None) -> np.ndarray:
    """Read a PRTSPM block into a dense array.

    ``PRTSPM`` alternates an index line -- up to five ``(row, col)`` pairs --
    with a value line holding the matching values.  It stops at the blank
    padding emitted by ``FORMAT(//)``.  Indices are 1-based Fortran.
    """
    n_cols = dim if n_cols is None else n_cols
    mat = np.zeros((dim, n_cols), dtype=float)

    i = start
    blanks = 0
    while i < len(lines):
        line = lines[i]
        pairs = _PAIR_RE.findall(line)

        if not pairs:
            if not line.strip():
                blanks += 1
                # FORMAT(//) writes the terminating blank padding.
                if blanks >= 2:
                    break
                i += 1
                continue
            # A non-blank, non-pair line means the block is over.
            break

        blanks = 0
        vals = _floats(lines[i + 1]) if i + 1 < len(lines) else []
        for (r, c), v in zip(pairs, vals):
            r, c = int(r) - 1, int(c) - 1
            if 0 <= r < dim and 0 <= c < n_cols:
                mat[r, c] = v
                if symmetrize and 0 <= c < dim and 0 <= r < n_cols:
                    mat[c, r] = v
        i += 2

    return mat


# ──────────────────────────────────────────────────────────────────────────
# Fock basis (PRNTST)
# ──────────────────────────────────────────────────────────────────────────

def _parse_states(lines: list[str], numsta: int):
    """Read the ``TOTAL NUM STATES`` / ``STATE NUMBER:`` blocks.

    Each state prints ``LNGTH,NMES,NBRB,NBRD,KMES,KBRB,KBRD`` then four rows of
    ``25I6``: quark flags, antiquark flags, 2k^+ momenta, flavors.
    """
    start = _find(lines, "TOTAL NUM STATES")
    if start < 0:
        return None, None, None, None

    idx = [i for i in range(start, len(lines)) if "STATE NUMBER:" in lines[i]]
    if not idx:
        return None, None, None, None
    idx = idx[:numsta]

    lengths, types, moms, flavors = [], [], [], []
    for i in idx:
        info = _floats(lines[i + 2])          # the 7I5 line
        lng = int(info[0])
        # i+3 is the ' ------- ' rule; the four data rows follow.
        rows = [[int(v) for v in _floats(lines[i + 4 + r])[:lng]] for r in range(4)]
        lengths.append(lng)
        types.append(rows[0])
        moms.append(rows[2])
        flavors.append(rows[3])

    maxlen = max(lengths) if lengths else 0
    n = len(lengths)

    def pad(rows):
        out = np.zeros((n, maxlen), dtype=int)
        for r, row in enumerate(rows):
            out[r, : len(row)] = row
        return out

    return np.array(lengths, dtype=int), pad(types), pad(moms), pad(flavors)


# ──────────────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────────────

def read_out(path, with_matrices: bool = True) -> DLCQResult:
    """Parse ``qcdf.out`` into a :class:`DLCQResult`.

    Set ``with_matrices=False`` to skip the norm and Z blocks when only the
    spectrum is wanted -- much faster on the large runs.
    """
    path = Path(path)
    lines = path.read_text(errors="replace").splitlines()

    meta = _parse_header(lines)
    n = meta["numsta_post"]

    eigenvalues = parse_eigenvalues(lines, expected=n)
    z_eig = _parse_eigenvectors(lines, n_rows=n)

    state_len = state_types = state_moms = state_flavors = None
    norm = Z = c_orig = None

    if with_matrices:
        state_len, state_types, state_moms, state_flavors = _parse_states(lines, n)

        i = _find(lines, "NORM MATRIX")
        if i >= 0:
            norm = _parse_sparse(lines, i + 1, n, symmetrize=True)

        i = _find(lines, "MATRIX Z WHICH GIVES NEW BASIS")
        if i >= 0:
            j = _find(lines, "ORIG. STATE NUM.", i)
            Z = _parse_sparse(lines, (j if j >= 0 else i) + 1, n, symmetrize=False)

        if Z is not None and z_eig is not None:
            # c_orig[:, k] = coefficients of eigenstate k in the original Fock
            # basis.  Z maps orthonormal -> original; z_eig lives in the
            # orthonormal basis.
            c_orig = Z @ z_eig

    return DLCQResult(
        N=meta["N"], NF=meta["NF"], B=meta["B"], K_code=meta["K_code"],
        rlamb=meta["rlamb"], cutoff=meta["cutoff"], LPN=meta["LPN"],
        rmq=meta["rmq"], iflv=meta["iflv"],
        source="fortran",
        provenance={"file": str(path.resolve()), "reader": "dlcq.read_fortran"},
        numsta_pre=meta["numsta_pre"], numsta_post=n,
        state_len=state_len, state_types=state_types,
        state_moms=state_moms, state_flavors=state_flavors,
        norm=norm, Z=Z, eigenvalues=eigenvalues,
        eigenvectors=z_eig, c_orig=c_orig,
    )


def read_ham(path):
    """Parse ``qcdf.ham`` -> ``(meta, HNU0, HNU)``.

    The stored Hamiltonian: ``HNU0`` is the free part (diagonal, per flavor),
    ``HNU`` the interacting part as ``(i, j, value)`` triplets terminated by a
    ``-1 -1 0.0`` sentinel.  Because the coupling only enters as a lambda^2
    reweighting at the end, one of these files serves an entire m/g sweep.
    """
    path = Path(path)
    lines = path.read_text(errors="replace").splitlines()

    hdr = _floats(lines[1])
    N, NF, B, LPN, K_code, numsta = (int(v) for v in hdr[:6])
    cutoff = float(hdr[6]) if len(hdr) > 6 else -1.0
    meta = dict(N=N, NF=NF, B=B, LPN=LPN, K_code=K_code,
                numsta=numsta, cutoff=cutoff)

    i = _find(lines, "HNU0")
    hnu0 = np.zeros((NF, numsta))
    row = 0
    i += 1
    while row < numsta and i < len(lines):
        vals = _floats(lines[i])
        if vals:
            hnu0[:, row] = vals[:NF]
            row += 1
        i += 1

    i = _find(lines, "HNU:", i - 1)
    hnu = np.zeros((numsta, numsta))
    for line in lines[i + 1:]:
        toks = _floats(line)
        for k in range(0, len(toks) - 2, 3):
            a, b, v = int(toks[k]), int(toks[k + 1]), toks[k + 2]
            if a == -1:
                return meta, hnu0, hnu
            hnu[a - 1, b - 1] = v
            hnu[b - 1, a - 1] = v

    return meta, hnu0, hnu
