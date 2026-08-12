"""Run the Python solver and emit the same :class:`DLCQResult` contract.

This is the second producer feeding ``dlcq.figures``.  It deliberately mirrors
the Fortran main program's sequence so the two can be compared step for step:

    qcdsta -> build norm -> weed -> Z = z_norm / sqrt(w_norm)  (NUZ)
           -> build Hamiltonian -> HNU = Z^T H Z  (NUHAM)
           -> HNU = lambda^2 HNU + (1-lambda^2) rmq^2 HNU0 on the diagonal
           -> eigh -> M^2 = K * w / 2

Weeding policy
--------------
Three mutually inconsistent policies existed in this codebase:

===================================  ==========================================
``qcdf.weedr2``                      drops on ``abs(w) < eps``   (matches Fortran)
``qcdf_opt.weed``                    drops on ``w < eps``        (also eats large
                                     negative eigenvalues)
``reproduce_figures.run_dlcq``       keeps ``w > 0.5`` instead of weeding at all
===================================  ==========================================

That is the direct cause of Fortran retaining 189 states where Python retained
191 at N=3, B=1, 2K=21.  ``policy="fortran"`` reproduces ``WEEDR``/``WEEDR2``
faithfully and is the default; ``policy="spectral"`` keeps the better
conditioned projection as an explicit opt-in.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from .dataset import DLCQResult

__all__ = ["run_python", "weed_fortran", "weed_spectral"]

_PYDIR = Path(__file__).resolve().parent.parent / "python"


def _import_solver(prefer_opt: bool = True):
    """Import the solver, preferring the optimized build."""
    if str(_PYDIR) not in sys.path:
        sys.path.insert(0, str(_PYDIR))
    if prefer_opt:
        try:
            import qcdf_opt
            import qcdf
            return qcdf_opt, qcdf, True
        except ImportError:
            pass
    import qcdf
    return None, qcdf, False


# ──────────────────────────────────────────────────────────────────────────
# Weeding
# ──────────────────────────────────────────────────────────────────────────

def weed_fortran(hnorm, mstinf, numsta, eps=1e-4, max_iter=2000):
    """Faithful port of Fortran ``WEEDR`` + ``WEEDR2``.

    Stage 1 (``WEEDR``) removes states whose norm row is proportional to an
    earlier row -- literal linear redundancy.

    Stage 2 (``WEEDR2``) diagonalizes the norm matrix, counts eigenvalues with
    ``abs(w) < eps`` (note: *absolute* value, unlike ``qcdf_opt.weed``), and for
    each null direction discards one state from its support, marking the rest
    so the discarded states stay independent.  Repeats until no null directions
    remain.
    """
    hnorm = np.array(hnorm, dtype=float, copy=True)
    mstinf = np.array(mstinf, copy=True)

    # ── WEEDR: linear redundancy ──
    loc = 0
    while loc < numsta - 1:
        drops = []
        for j in range(loc + 1, numsta):
            if abs(hnorm[j, loc]) > eps:
                r = hnorm[j, loc] / hnorm[loc, loc]
                if np.all(np.abs(r * hnorm[loc, :numsta] - hnorm[j, :numsta]) <= eps):
                    drops.append(j)
        for d in sorted(drops, reverse=True):
            hnorm = np.delete(np.delete(hnorm, d, 0), d, 1)
            mstinf = np.delete(mstinf, d, 0)
            numsta -= 1
        loc += 1

    # ── WEEDR2: null directions of the norm matrix ──
    for _ in range(max_iter):
        w, z = eigh(hnorm[:numsta, :numsta])
        nzer = int(np.sum(np.abs(w) < eps))       # abs(), matching Fortran
        if nzer == 0:
            break
        drops, used = [], set()
        for i in range(nzer):
            dropped = False
            # Fortran scans I2 = NUMSTA..1 descending.
            for j in range(numsta - 1, -1, -1):
                if abs(z[j, i]) > eps:
                    if not dropped and j not in used:
                        drops.append(j)
                        dropped = True
                    used.add(j)
        if not drops:
            break
        for d in sorted(set(drops), reverse=True):
            hnorm = np.delete(np.delete(hnorm, d, 0), d, 1)
            mstinf = np.delete(mstinf, d, 0)
            numsta -= 1

    return hnorm, mstinf, numsta


def weed_spectral(hnorm, mstinf, numsta, threshold=0.5):
    """Keep only norm eigen-directions above ``threshold``; no state deletion.

    Better conditioned than iterative weeding and useful in its own right, but
    it produces a different basis from the Fortran, so it must not be used for
    tight Fortran-vs-Python comparison.
    """
    w, z = eigh(hnorm[:numsta, :numsta])
    good = w > threshold
    return hnorm, mstinf, numsta, z[:, good] / np.sqrt(w[good])[np.newaxis, :]


# ──────────────────────────────────────────────────────────────────────────
# Driver
# ──────────────────────────────────────────────────────────────────────────

def run_python(N, NF, B, K_code, rlamb, cutoff=-1.0, LPN=0,
               rmq=None, iflv=None, ncpus=1, policy="fortran",
               assembly="exact", prefer_opt=True) -> DLCQResult:
    """Run the Python DLCQ solver and return a :class:`DLCQResult`.

    Parameters
    ----------
    assembly : {"exact", "fortran"}
        How the free part of the Hamiltonian is combined with the interacting
        part in the orthonormal basis.

        ``"exact"`` (default) adds the full ``Z^T H0 Z`` matrix.  This is
        basis independent -- two different valid orthonormalizations agree to
        6e-13 -- and is therefore the physically meaningful answer.

        ``"fortran"`` adds only its diagonal, reproducing what ``qcdf.f`` does.
        That step is **basis dependent**: ``Z^T H0 Z`` has off-diagonal
        elements up to ~0.8, and the norm matrix is heavily degenerate (26
        distinct eigenvalues among 189 states at 2K=21), so the result depends
        on which eigenvectors the diagonalization returns.

        This is not a hypothetical.  Recompiling the *unmodified* ``qcdf.f``
        with ``-O2`` instead of ``-O0`` changes its own answer: 190 retained
        states instead of 189, and the ground state moves from 10.390380 to
        10.390084.  The Fortran spectrum is therefore reproducible only to
        ~1e-4, by anyone, including itself.  See ``docs/basis-dependence.md``.

        The shift is far below the precision the paper quotes (Table I gives
        M_bar/g = 10.71(2)), so the published physics is unaffected.

    ``rlamb`` is taken literally.  Do not derive it from m/g when comparing
    against a Fortran run: the input files use 0.3325 while
    ``mg_to_lambda(1.6)`` is 0.3325495, and that 1.5e-5 difference moves
    eigenvalues well above a 1e-8 comparison tolerance.
    """
    if assembly not in ("exact", "fortran"):
        raise ValueError(f"unknown assembly {assembly!r}")
    opt, base, use_opt = _import_solver(prefer_opt)
    os.environ["QCDF_NCPUS"] = str(ncpus)

    if iflv is None:
        iflv = [0] * max(NF, 1)
        iflv[0] = N * B                      # matches qcdf.f line 250

    p = base.Params()
    p.N, p.NF, p.B, p.K, p.rlamb = N, NF, B, K_code, rlamb
    p.cutoff, p.LPN = cutoff, LPN
    if rmq is not None:
        for i, v in enumerate(rmq[:NF]):
            p.rmq[i] = v
    for i, v in enumerate(iflv[:NF]):
        p.iflv[i] = v

    states = base.StateData()
    perm = base.PermTables()
    flav = base.FlavorTables()
    selfen = opt.compute_selfen(N) if use_opt else base.compute_selfen(p)

    base.qcdsta(p, states, perm, flav)
    numsta_pre = states.numsta
    if numsta_pre == 0:
        return DLCQResult(N=N, NF=NF, B=B, K_code=K_code, rlamb=rlamb,
                          cutoff=cutoff, LPN=LPN, source="python",
                          numsta_pre=0, numsta_post=0,
                          eigenvalues=np.array([]))

    # ── norm matrix ──
    if use_opt:
        _, _, hnorm = opt.build_matrices(0, states.mstate, states.mstinf,
                                         numsta_pre, N, NF, B, K_code,
                                         selfen, p.cbreak, ncpus)
    else:
        _, _, hnorm = base.clrdis(0, p, states, selfen, ncpus=ncpus)
    mstinf = states.mstinf[:numsta_pre].copy()

    # ── weed ──
    if policy == "fortran":
        hnorm_w, mstinf_w, numsta = weed_fortran(hnorm, mstinf, numsta_pre)
        w_n, z_n = eigh(hnorm_w[:numsta, :numsta])
        Z = z_n / np.sqrt(w_n)[np.newaxis, :]          # NUZ
    elif policy == "spectral":
        hnorm_w, mstinf_w, numsta, Z = weed_spectral(hnorm, mstinf, numsta_pre)
    else:
        raise ValueError(f"unknown weeding policy {policy!r}")

    if numsta == 0:
        return DLCQResult(N=N, NF=NF, B=B, K_code=K_code, rlamb=rlamb,
                          cutoff=cutoff, LPN=LPN, source="python",
                          numsta_pre=numsta_pre, numsta_post=0,
                          eigenvalues=np.array([]))

    # ── Hamiltonian on the surviving basis ──
    states.numsta = numsta
    states.mstinf[:numsta] = mstinf_w[:numsta]
    if use_opt:
        ham0, ham, _ = opt.build_matrices(1, states.mstate, mstinf_w[:numsta],
                                          numsta, N, NF, B, K_code,
                                          selfen, p.cbreak, ncpus)
    else:
        ham0, ham, _ = base.clrdis(1, p, states, selfen, ncpus=ncpus)

    # ── NUHAM, then combine free + interacting ──
    hnu = Z.T @ ham[:numsta, :numsta] @ Z
    free = [Z.T @ ham0[ifl, :numsta, :numsta] @ Z for ifl in range(NF)]

    rlmsq = rlamb ** 2
    hnu = hnu * rlmsq
    for ifl in range(NF):
        rmq_val = p.rmq[ifl]
        coeff = (1 - rlmsq) * rmq_val ** 2
        if assembly == "exact":
            hnu = hnu + coeff * free[ifl]
        else:                                    # "fortran": diagonal only
            hnu[np.diag_indices(hnu.shape[0])] += coeff * np.diag(free[ifl])

    w_eig, z_eig = eigh(hnu)
    eigenvalues = K_code * w_eig / 2.0
    c_orig = Z @ z_eig

    # ── Fock content of the surviving basis ──
    lengths = mstinf_w[:numsta, 1].astype(int)
    maxlen = int(lengths.max()) if numsta else 0
    types = np.zeros((numsta, maxlen), dtype=int)
    moms = np.zeros((numsta, maxlen), dtype=int)
    flavors = np.zeros((numsta, maxlen), dtype=int)
    for s in range(numsta):
        loc = int(mstinf_w[s, 0]) - 1
        L = int(lengths[s])
        types[s, :L] = states.mstate[loc, :L]
        moms[s, :L] = states.mstate[loc + 2, :L]
        flavors[s, :L] = states.mstate[loc + 3, :L]

    return DLCQResult(
        N=N, NF=NF, B=B, K_code=K_code, rlamb=rlamb, cutoff=cutoff, LPN=LPN,
        rmq=np.asarray(p.rmq[:NF]), iflv=np.asarray(p.iflv[:NF]),
        source="python",
        provenance={"solver": "qcdf_opt" if use_opt else "qcdf",
                    "weeding_policy": policy, "assembly": assembly,
                    "ncpus": ncpus},
        numsta_pre=numsta_pre, numsta_post=numsta,
        state_len=lengths, state_types=types, state_moms=moms,
        state_flavors=flavors,
        norm=hnorm_w[:numsta, :numsta], Z=Z,
        eigenvalues=eigenvalues, eigenvectors=z_eig, c_orig=c_orig,
    )
