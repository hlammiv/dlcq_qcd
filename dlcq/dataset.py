"""The Fortran/Python data contract.

``DLCQResult`` is the single structure every producer emits and every consumer
reads.  The figure code depends on this and never on a solver, so a figure
built from Fortran output and one built from the Python port go through
identical arithmetic -- any difference in the plot is a difference in the
physics, not in the plotting.

Layout notes that are easy to get wrong:

``state_types``   1 = quark (b^dag), 0 = antiquark (d^dag)
``state_moms``    2k^+, i.e. odd integers under antiperiodic boundary conditions
``eigenvalues``   M^2_code = K_code * w / 2, ascending, in units of m^2 + g^2/pi
``eigenvectors``  columns are states, expressed in the *orthonormal* basis
``c_orig``        Z @ eigenvectors -- coefficients in the original Fock basis,
                  which is what the structure-function code actually needs
``norm``          the Gram matrix of the original (non-orthogonal) basis

The Fortran prints at most ``IVCTMX = 75`` eigenvectors regardless of how many
eigenvalues it computed, so ``n_eigenvectors`` can be smaller than
``n_eigenvalues``.  Consumers must check :meth:`DLCQResult.has_eigenvector`
rather than assuming.
"""

from __future__ import annotations

import json
import os
import tempfile
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

import numpy as np

__all__ = ["DLCQResult", "save", "load"]

SCHEMA_VERSION = 1


@dataclass
class BlockDiagonal:
    """A block-diagonal matrix kept as its blocks, never materialized.

    The norm is *exactly* block diagonal in parton configuration -- off-block
    entries are ``0.0``, not 1e-16 -- so the dense form is overwhelmingly zeros
    that still have to be allocated.  At LPN=7, 2K=71 the weeded norm is
    142352 x 142352, which is **151 GB** dense and a few MB as blocks.  That
    array was the last thing standing between this solver and the large-K runs.

    Supports ``@`` and nothing else, because ``@`` is all any consumer needs:
    :func:`dlcq.observables.structure_function` uses the norm solely as the
    matvec ``norm @ c``.  Indexing raises rather than silently densifying -- the
    whole point is that the dense form does not fit, so a stray ``norm[i, j]``
    must be a loud failure rather than an out-of-memory much later.

    Not ``scipy.sparse``: that would give ``@`` for free but discard the block
    structure, which is load-bearing three separate times -- weeding runs per
    block, ``orthonormalize_blockwise`` needs every column of Z inside one
    block, and the ``H0 = D N`` shortcut needs D constant on a block.
    """

    n_rows: int
    n_cols: int
    rows: list = field(default_factory=list)     # index arrays, one per block
    blocks: list = field(default_factory=list)   # dense square-ish blocks

    @property
    def shape(self):
        return (self.n_rows, self.n_cols)

    @property
    def nbytes(self):
        return int(sum(b.nbytes for b in self.blocks))

    def __matmul__(self, x):
        x = np.asarray(x)
        out = np.zeros((self.n_rows,) + x.shape[1:], dtype=float)
        for idx, blk in zip(self.rows, self.blocks):
            out[idx] = blk @ x[idx]
        return out

    def diagonal(self):
        d = np.zeros(min(self.n_rows, self.n_cols))
        for idx, blk in zip(self.rows, self.blocks):
            d[idx] = np.diag(blk)
        return d

    def toarray(self):
        """The dense form, explicitly.  Never called implicitly."""
        out = np.zeros((self.n_rows, self.n_cols))
        for idx, blk in zip(self.rows, self.blocks):
            out[np.ix_(idx, idx)] = blk
        return out

    def __getitem__(self, key):
        raise TypeError(
            "BlockDiagonal does not support indexing -- it exists because the "
            "dense form does not fit (151 GB at LPN=7, 2K=71). Use `@` for the "
            "matvec, .diagonal(), or .toarray() if you can afford it.")

    def __array__(self, dtype=None):
        raise TypeError(
            "refusing to densify a BlockDiagonal implicitly; call .toarray() "
            "if that is really what you want.")


@dataclass
class DLCQResult:
    """One DLCQ diagonalization, from either solver."""

    # ── run parameters ────────────────────────────────────────────────────
    N: int                      # number of colors
    NF: int                     # number of flavors
    B: int                      # baryon number
    K_code: int                 # = 2 * the paper's K
    rlamb: float                # dimensionless coupling lambda
    cutoff: float = -1.0        # Pauli-Villars invariant-mass cutoff (<=0: none)
    LPN: int = 0                # particle-number truncation (0: none)
    rmq: np.ndarray = field(default_factory=lambda: np.array([1.0]))
    iflv: np.ndarray = field(default_factory=lambda: np.array([0]))

    # ── provenance ────────────────────────────────────────────────────────
    source: str = "unknown"     # "fortran" | "python"
    provenance: dict = field(default_factory=dict)

    # ── basis ─────────────────────────────────────────────────────────────
    numsta_pre: int = 0         # states before weeding
    numsta_post: int = 0        # states after weeding (rows of the Fock basis)
    n_orth: int = 0             # columns of Z = dimension of the orthonormal
                                # space, and so the length of a *full* spectrum.
                                # 0 means "not recorded" (schema-1 files, which
                                # load() repairs from Z).
    state_len: Optional[np.ndarray] = None      # (numsta_post,) partons per state
    state_types: Optional[np.ndarray] = None    # (numsta_post, maxlen) 1=q 0=qbar
    state_moms: Optional[np.ndarray] = None     # (numsta_post, maxlen) odd ints
    state_flavors: Optional[np.ndarray] = None  # (numsta_post, maxlen)

    # ── matrices ──────────────────────────────────────────────────────────
    norm: Optional[np.ndarray] = None           # (numsta_post, numsta_post)
    Z: Optional[np.ndarray] = None              # (numsta_post, n_orth)
    eigenvalues: Optional[np.ndarray] = None    # (n_eig,) ascending M^2_code
    eigenvectors: Optional[np.ndarray] = None   # (n_orth, n_vec) orthonormal basis
    c_orig: Optional[np.ndarray] = None         # (numsta_post, n_vec) Fock basis

    # ── derived conveniences ──────────────────────────────────────────────

    @property
    def K_paper(self) -> float:
        return self.K_code / 2.0

    @property
    def n_eigenvalues(self) -> int:
        return 0 if self.eigenvalues is None else int(self.eigenvalues.size)

    @property
    def n_eigenvectors(self) -> int:
        if self.c_orig is not None:
            return int(self.c_orig.shape[1])
        if self.eigenvectors is None:
            return 0
        return int(self.eigenvectors.shape[1])

    @property
    def spectrum_is_truncated(self) -> bool:
        """True when ``eigenvalues`` holds the lowest few of a larger spectrum.

        ``eigenvalues[i]`` is still the i-th lowest level of the full spectrum;
        what is missing is everything above.  So indexing near the bottom stays
        correct and indexing deep silently runs off the end -- which is exactly
        the failure a consumer cannot detect from the array alone, since a run
        with 40 levels and a run that *computed* 40 of 15235 look identical.

        Derived from ``n_orth`` rather than asserted by a flag, so it cannot
        disagree with the data it describes.
        """
        return self.n_orth > 0 and self.n_eigenvalues < self.n_orth

    def has_eigenvector(self, idx: int) -> bool:
        """Whether eigenvector ``idx`` is actually present.

        Guards the Fortran's IVCTMX=75 print cap: a run can report 189
        eigenvalues but only 75 eigenvectors.
        """
        return 0 <= idx < self.n_eigenvectors

    def require_eigenvector(self, idx: int) -> np.ndarray:
        """Return ``c_orig[:, idx]``, failing loudly if it was never printed."""
        if not self.has_eigenvector(idx):
            raise IndexError(
                f"eigenvector {idx} unavailable: this run carries "
                f"{self.n_eigenvectors} eigenvectors for {self.n_eigenvalues} "
                f"eigenvalues (the Fortran caps printed eigenvectors at "
                f"IVCTMX=75). Raise IVCTMX and rerun, or pick a lower state."
            )
        if self.c_orig is None:
            raise ValueError("c_orig not populated")
        return self.c_orig[:, idx]

    def describe(self) -> str:
        return (
            f"DLCQResult(source={self.source}, N={self.N}, NF={self.NF}, "
            f"B={self.B}, K_code={self.K_code} (2K), lambda={self.rlamb:.6g}, "
            f"states {self.numsta_pre}->{self.numsta_post}, "
            f"{self.n_eigenvalues} eigenvalues, {self.n_eigenvectors} eigenvectors)"
        )

    # ── consistency ───────────────────────────────────────────────────────

    def validate(self) -> list[str]:
        """Return a list of structural problems; empty means self-consistent."""
        problems = []
        n = self.numsta_post

        if self.state_types is not None and self.state_types.shape[0] != n:
            problems.append(
                f"state_types has {self.state_types.shape[0]} rows, expected {n}"
            )
        if self.norm is not None and self.norm.shape != (n, n):
            problems.append(f"norm is {self.norm.shape}, expected {(n, n)}")
        if self.Z is not None and self.Z.shape[0] != n:
            problems.append(f"Z has {self.Z.shape[0]} rows, expected {n}")
        if self.c_orig is not None and self.c_orig.shape[0] != n:
            problems.append(f"c_orig has {self.c_orig.shape[0]} rows, expected {n}")

        if self.eigenvalues is not None and self.eigenvalues.size > 1:
            if np.any(np.diff(self.eigenvalues) < -1e-9):
                problems.append("eigenvalues are not ascending")

        # A truncated spectrum is self-consistent and must NOT be flagged here;
        # what is checked is that n_orth agrees with the arrays it describes.
        if self.n_orth:
            if self.Z is not None and self.Z.shape[1] != self.n_orth:
                problems.append(
                    f"Z has {self.Z.shape[1]} columns, expected "
                    f"n_orth={self.n_orth}")
            if (self.eigenvectors is not None
                    and self.eigenvectors.shape[0] != self.n_orth):
                problems.append(
                    f"eigenvectors has {self.eigenvectors.shape[0]} rows, "
                    f"expected n_orth={self.n_orth}")
            if self.n_eigenvalues > self.n_orth:
                problems.append(
                    f"{self.n_eigenvalues} eigenvalues exceeds "
                    f"n_orth={self.n_orth}")

        if (self.eigenvectors is not None and self.c_orig is not None
                and self.eigenvectors.shape[1] != self.c_orig.shape[1]):
            problems.append(
                "eigenvectors and c_orig disagree on the number of states")

        if self.state_moms is not None and self.state_len is not None:
            for s in range(min(n, self.state_moms.shape[0])):
                moms = self.state_moms[s, : self.state_len[s]]
                if np.any(moms % 2 == 0):
                    problems.append(
                        f"state {s} has even momenta {moms.tolist()}; "
                        "antiperiodic BCs require odd"
                    )
                    break
                if int(moms.sum()) != self.K_code:
                    problems.append(
                        f"state {s} momenta sum to {int(moms.sum())}, "
                        f"expected K_code={self.K_code}"
                    )
                    break

        return problems


# ──────────────────────────────────────────────────────────────────────────
# HDF5 serialization
# ──────────────────────────────────────────────────────────────────────────

_SCALARS = ("N", "NF", "B", "K_code", "rlamb", "cutoff", "LPN",
            "numsta_pre", "numsta_post", "source")

#: Written like ``_SCALARS`` but read with a default, because ``load`` takes
#: every ``_SCALARS`` entry unconditionally -- adding ``n_orth`` there would
#: raise ``KeyError`` on every file written before it existed, of which there
#: are several hundred in ``runs/python_cache`` representing hours of compute.
_SCALARS_OPTIONAL = (("n_orth", 0),)
_ARRAYS = ("rmq", "iflv", "state_len", "state_types", "state_moms",
           "state_flavors", "norm", "Z", "eigenvalues", "eigenvectors",
           "c_orig")


@dataclass(frozen=True)
class Spectrum:
    """Just the levels of a cached run, without inflating the rest of it.

    ``Z`` is stored dense.  It is 0.014% non-zero at 2K=71, so an 8.1 MB file
    on disk becomes 689 MB in memory the instant it is read -- and the Richardson
    sweeps, ``figures.table1`` and ``tools/plot_holdout.py`` alike, use nothing
    from a result but ``eigenvalues`` and ``rlamb``.  Reading one Table I panel
    was costing ~700 MB per K.

    ``spurious_zero_modes`` needs exactly these two fields, so this duck-types
    into :func:`~dlcq.observables.physical_indices` unchanged.

    The real fix is to store ``Z`` sparsely, which is a schema change; this
    removes the cost from every caller that never wanted ``Z`` in the first
    place.
    """
    eigenvalues: np.ndarray
    rlamb: float
    n_orth: int = 0

    @property
    def n_eigenvalues(self) -> int:
        return 0 if self.eigenvalues is None else int(self.eigenvalues.size)

    @property
    def spectrum_is_truncated(self) -> bool:
        return bool(self.n_orth) and self.n_eigenvalues < self.n_orth


def read_spectrum(path) -> Spectrum:
    """Read only the spectrum of a cached result.  See :class:`Spectrum`."""
    import h5py

    with h5py.File(Path(path), "r") as f:
        ev = f["eigenvalues"][:] if "eigenvalues" in f else np.zeros(0)
        return Spectrum(eigenvalues=np.asarray(ev),
                        rlamb=float(f.attrs["rlamb"]),
                        n_orth=int(f.attrs.get("n_orth", 0)))


def save(result: DLCQResult, path) -> Path:
    """Write a :class:`DLCQResult` to HDF5, atomically.

    Written to a temporary file in the same directory and then ``os.replace``d
    into place, which is atomic on POSIX.  Writing straight to the destination
    is what left a 1.34 GB unreadable file in ``runs/python_cache`` ("bad object
    header version number" -- a partial write, killed mid-flush by the OOM
    reaper during a large run).  Once such a file exists nothing removes it:
    ``PythonProvider.get`` sees ``path.exists()`` and hands it to ``load``,
    which raises, so every future request for that parameter set fails.

    Same directory rather than the system temp dir, so the rename cannot cross
    a filesystem boundary and silently degrade to a copy.
    """
    import h5py

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix=f".{path.name}.",
                               suffix=".tmp")
    os.close(fd)
    tmp = Path(tmp)
    try:
        with h5py.File(tmp, "w") as f:
            f.attrs["schema_version"] = SCHEMA_VERSION
            for name in _SCALARS:
                f.attrs[name] = getattr(result, name)
            for name, _default in _SCALARS_OPTIONAL:
                f.attrs[name] = getattr(result, name)
            f.attrs["provenance"] = json.dumps(result.provenance, default=str)
            for name in _ARRAYS:
                value = getattr(result, name)
                if value is None:
                    continue
                if isinstance(value, BlockDiagonal):
                    # Blocks concatenated, plus their row indices and sizes.
                    # Written under a distinct name so a reader that predates
                    # this simply does not find the field, rather than reading
                    # it as something it is not.
                    g = f.create_group(f"{name}_blocks")
                    g.attrs["n_rows"] = value.n_rows
                    g.attrs["n_cols"] = value.n_cols
                    g.create_dataset("sizes",
                                     data=np.array([len(r) for r in value.rows],
                                                   dtype=np.int64))
                    g.create_dataset("rowindex",
                                     data=np.concatenate(value.rows).astype(
                                         np.int64) if value.rows
                                     else np.zeros(0, dtype=np.int64),
                                     compression="gzip")
                    g.create_dataset("data",
                                     data=np.concatenate(
                                         [b.ravel() for b in value.blocks])
                                     if value.blocks else np.zeros(0),
                                     compression="gzip")
                    continue
                f.create_dataset(name, data=np.asarray(value),
                                 compression="gzip")
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    return path


def load(path) -> DLCQResult:
    """Read a :class:`DLCQResult` back from HDF5."""
    import h5py

    with h5py.File(Path(path), "r") as f:
        kwargs = {}
        for name in _SCALARS:
            value = f.attrs[name]
            if name == "source":
                value = value.decode() if isinstance(value, bytes) else str(value)
            elif name in ("rlamb", "cutoff"):
                value = float(value)
            else:
                value = int(value)
            kwargs[name] = value
        kwargs["provenance"] = json.loads(f.attrs.get("provenance", "{}"))
        for name, default in _SCALARS_OPTIONAL:
            kwargs[name] = int(f.attrs.get(name, default))
        for name in _ARRAYS:
            if name in f:
                kwargs[name] = f[name][...]
            elif f"{name}_blocks" in f:
                g = f[f"{name}_blocks"]
                sizes = g["sizes"][...]
                rowindex = g["rowindex"][...]
                data = g["data"][...]
                rows, blocks, ri, di = [], [], 0, 0
                for b in sizes:
                    b = int(b)
                    rows.append(rowindex[ri:ri + b])
                    blocks.append(data[di:di + b * b].reshape(b, b))
                    ri += b
                    di += b * b
                kwargs[name] = BlockDiagonal(int(g.attrs["n_rows"]),
                                             int(g.attrs["n_cols"]),
                                             rows, blocks)
        # Repair schema-1 files rather than leaving them ambiguous: Z's column
        # count *is* n_orth for every result written before the field existed,
        # since all of them kept the full spectrum.  Without this they would
        # report n_orth=0, and spectrum_is_truncated would answer False for the
        # wrong reason -- indistinguishable from a genuinely complete run.
        if not kwargs.get("n_orth") and "Z" in f:
            kwargs["n_orth"] = int(f["Z"].shape[1])
        return DLCQResult(**kwargs)
