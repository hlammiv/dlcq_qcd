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
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

import numpy as np

__all__ = ["DLCQResult", "save", "load"]

SCHEMA_VERSION = 1


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
_ARRAYS = ("rmq", "iflv", "state_len", "state_types", "state_moms",
           "state_flavors", "norm", "Z", "eigenvalues", "eigenvectors",
           "c_orig")


def save(result: DLCQResult, path) -> Path:
    """Write a :class:`DLCQResult` to HDF5."""
    import h5py

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(path, "w") as f:
        f.attrs["schema_version"] = SCHEMA_VERSION
        for name in _SCALARS:
            f.attrs[name] = getattr(result, name)
        f.attrs["provenance"] = json.dumps(result.provenance, default=str)
        for name in _ARRAYS:
            value = getattr(result, name)
            if value is None:
                continue
            f.create_dataset(name, data=np.asarray(value), compression="gzip")
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
        for name in _ARRAYS:
            if name in f:
                kwargs[name] = f[name][...]
        return DLCQResult(**kwargs)
