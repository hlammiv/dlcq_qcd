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

__all__ = ["run_python", "weed_fortran", "weed_spectral",
           "orthonormalize_blockwise"]

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

def _weed_indices(A, eps=1e-4, max_iter=2000):
    """Which rows of ``A`` survive ``WEEDR`` + ``WEEDR2``.

    The algorithm itself, factored out so it can be run on the whole norm
    matrix or on one block of it.  Returns surviving indices into ``A``,
    ascending.
    """
    M = np.array(A, dtype=float, copy=True)
    idx = list(range(M.shape[0]))
    n = len(idx)

    # ── WEEDR: literal linear redundancy, a row proportional to an earlier one ──
    loc = 0
    while loc < n - 1:
        drops = []
        for j in range(loc + 1, n):
            if abs(M[j, loc]) > eps:
                r = M[j, loc] / M[loc, loc]
                if np.all(np.abs(r * M[loc, :n] - M[j, :n]) <= eps):
                    drops.append(j)
        for d in sorted(drops, reverse=True):
            M = np.delete(np.delete(M, d, 0), d, 1)
            idx.pop(d)
            n -= 1
        loc += 1

    # ── WEEDR2: null directions of the norm matrix ──
    for _ in range(max_iter):
        w, z = eigh(M[:n, :n])
        nzer = int(np.sum(np.abs(w) < eps))       # abs(), matching Fortran
        if nzer == 0:
            break
        drops, used = [], set()
        for i in range(nzer):
            dropped = False
            # Fortran scans I2 = NUMSTA..1 descending.
            for j in range(n - 1, -1, -1):
                if abs(z[j, i]) > eps:
                    if not dropped and j not in used:
                        drops.append(j)
                        dropped = True
                    used.add(j)
        if not drops:
            break
        for d in sorted(set(drops), reverse=True):
            M = np.delete(np.delete(M, d, 0), d, 1)
            idx.pop(d)
            n -= 1

    return idx


def weed_fortran(hnorm, mstinf, numsta, eps=1e-4, max_iter=2000, blocked=True):
    """Faithful port of Fortran ``WEEDR`` + ``WEEDR2``.

    Stage 1 (``WEEDR``) removes states whose norm row is proportional to an
    earlier row -- literal linear redundancy.

    Stage 2 (``WEEDR2``) diagonalizes the norm matrix, counts eigenvalues with
    ``abs(w) < eps`` (note: *absolute* value, unlike ``qcdf_opt.weed``), and for
    each null direction discards one state from its support, marking the rest
    so the discarded states stay independent.  Repeats until no null directions
    remain.

    ``blocked=True`` runs that identical algorithm **per connected block of the
    norm matrix** instead of on all states at once, which is a large speedup and
    was the dominant cost of a run: 76% of the wall time at 2K=23, against 16%
    for the Hamiltonian build and 6% for the norm.

    It is exact, not an approximation.  The norm is block-diagonal in momentum
    configuration -- verified by explicit colour enumeration, see
    ``tools/colour_norm.py`` -- so a state can only be linearly dependent on
    states in its own block, and both stages decompose:

    * ``WEEDR`` compares row j against row loc only when ``|N[j,loc]| > eps``,
      which already restricts it to a single block;
    * ``WEEDR2``'s null space is the direct sum of the blocks' null spaces.

    The one place the two can differ is if ``eigh`` returns null eigenvectors
    that mix several blocks, which it may when their zero eigenvalues are
    degenerate.  Then the two runs drop a different member of the same null
    space -- the surviving span, and so the spectrum, is unchanged.
    ``tests/test_fortran_python.py`` checks both the counts and the spectrum.

    At 2K=23 the 897 states form 180 blocks of median size 4 and maximum 19, so
    the O(n^2) work drops by about 120x.
    """
    A = np.asarray(hnorm, dtype=float)[:numsta, :numsta]
    mstinf = np.array(mstinf, copy=True)

    if blocked and numsta > 1:
        from scipy.sparse import csr_matrix
        from scipy.sparse.csgraph import connected_components

        adj = csr_matrix((np.abs(A) > 1e-9).astype(np.int8))
        nblocks, labels = connected_components(adj, directed=False)
        keep = []
        for b in range(nblocks):
            members = np.flatnonzero(labels == b)
            sub = _weed_indices(A[np.ix_(members, members)], eps, max_iter)
            keep.extend(int(members[i]) for i in sub)
        keep.sort()
    else:
        keep = _weed_indices(A, eps, max_iter)

    keep = np.asarray(keep, dtype=int)
    return A[np.ix_(keep, keep)].copy(), mstinf[keep], int(keep.size)


def orthonormalize_blockwise(hnorm, numsta, tol=1e-8):
    """Build Z block by block, which makes the Fortran's assembly well-posed.

    ``qcdf.f`` adds the free Hamiltonian to the diagonal only of ``Z^T H0 Z``,
    which is exact only if that matrix is diagonal.  With a *global*
    orthonormalization it is not -- the norm matrix is massively degenerate, so
    an eigensolver may mix states from different Fock momentum classes, and the
    discarded off-diagonal terms make the answer depend on which eigenvectors
    came back (see docs/basis-dependence.md).

    But ``H0 = D N`` with ``D`` **constant on each block** of mutually
    overlapping states: overlapping states share their parton content, hence
    their free energy.  Measured at N=3, B=1, 2K=21 the norm matrix splits into
    121 blocks (largest 6) and the free energy is constant across every one.

    Orthonormalizing inside each block therefore keeps every column of Z within
    a single block, so ``Z^T H0 Z`` comes out diagonal to machine precision
    (3.4e-14, against ~0.8 for a global Z) and the diagonal-only rule stops
    being an approximation.  The residual freedom -- rotations *within* a block
    -- does not matter, because D is a multiple of the identity there.

    Returns ``(Z, kept_columns)``.
    """
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components

    A = sp.csr_matrix((np.abs(hnorm[:numsta, :numsta]) > 1e-9).astype(np.int8))
    nblocks, labels = connected_components(A, directed=False)

    cols = []
    for b in range(nblocks):
        idx = np.flatnonzero(labels == b)
        w, v = eigh(hnorm[np.ix_(idx, idx)])
        keep = w > tol
        if not np.any(keep):
            continue
        Zb = v[:, keep] / np.sqrt(w[keep])
        for j in range(Zb.shape[1]):
            col = np.zeros(numsta)
            col[idx] = Zb[:, j]
            cols.append(col)

    Z = np.column_stack(cols) if cols else np.zeros((numsta, 0))
    return Z, Z.shape[1]


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
               assembly="exact", prefer_opt=True, backend=None) -> DLCQResult:
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

    backend : {"thread", "process", None}
        How the matrix builds are parallelised.  ``"thread"`` (the default)
        runs the ``nogil`` numba kernels of ``qcdf_kernels.py`` in a thread
        pool; ``"process"`` runs the interpreted reference routines under
        multiprocessing.  The two produce bit-identical matrices -- that is
        asserted in ``tests/test_kernels.py`` -- so this is purely a speed
        knob and results stay cache-compatible across it.  ``None`` takes
        ``QCDF_BACKEND`` from the environment, else the default.

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
                                         selfen, p.cbreak, ncpus,
                                         backend=backend)
    else:
        _, _, hnorm = base.clrdis(0, p, states, selfen, ncpus=ncpus)
    mstinf = states.mstinf[:numsta_pre].copy()

    # ── weed ──
    if policy == "fortran":
        hnorm_w, mstinf_w, numsta = weed_fortran(hnorm, mstinf, numsta_pre)
        w_n, z_n = eigh(hnorm_w[:numsta, :numsta])
        Z = z_n / np.sqrt(w_n)[np.newaxis, :]          # NUZ
    elif policy == "blockwise":
        hnorm_w, mstinf_w, numsta = weed_fortran(hnorm, mstinf, numsta_pre)
        Z, _ = orthonormalize_blockwise(hnorm_w, numsta)
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
                                          selfen, p.cbreak, ncpus,
                                          backend=backend)
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
                    "ncpus": ncpus,
                    "backend": backend or os.environ.get("QCDF_BACKEND")
                    or getattr(opt, "DEFAULT_BACKEND", "process")},
        numsta_pre=numsta_pre, numsta_post=numsta,
        state_len=lengths, state_types=types, state_moms=moms,
        state_flavors=flavors,
        norm=hnorm_w[:numsta, :numsta], Z=Z,
        eigenvalues=eigenvalues, eigenvectors=z_eig, c_orig=c_orig,
    )
