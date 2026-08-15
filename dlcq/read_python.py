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
from scipy import sparse
from scipy.linalg import eigh

from .dataset import BlockDiagonal, DLCQResult

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


def _labels_from_matrix(A):
    """Block labels via connected components of ``|A| > 1e-9``.

    The fallback for callers that cannot supply the state list.  Materializes
    dense ``n x n`` temporaries, so prefer :func:`config_block_labels`.
    """
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components

    adj = csr_matrix((np.abs(A) > 1e-9).astype(np.int8))
    return connected_components(adj, directed=False)[1]


def config_block_labels(mstate, mstinf, numsta):
    """Block labels for the norm matrix, straight from the state list.

    The norm couples two states only if their parton content is identical, so
    the multiset of ``(type, momentum, flavour)`` labels the blocks -- in
    ``O(n L)``, without touching the matrix.  This is the same grouping
    ``qcdf_opt._config_keys`` uses to schedule the norm build.

    It is *coarser* than the connected components of ``|N| > 1e-9``: a few
    configuration groups split further because their states happen to be
    colour-orthogonal (119 groups vs 120 components at 2K=21, 567 vs 577 at
    2K=29).  Coarser is safe -- weeding a superset of a block is still exact,
    just marginally more work -- and measured, both retain exactly the rank of
    the Gram matrix at 2K=21/25/29 and B=0/1.

    The alternative, ``csr_matrix(np.abs(A) > 1e-9)``, has to materialize two
    more dense ``n x n`` arrays to recover labels it could have had for free.
    At 2K=37 the norm is already 8.6 GB, so that is what runs a 15 GB machine
    out of memory.
    """
    keys = {}
    labels = np.empty(numsta, dtype=np.int64)
    for s in range(numsta):
        loc = int(mstinf[s, 0]) - 1
        L = int(mstinf[s, 1])
        key = tuple(sorted(
            (int(mstate[loc, j]), int(mstate[loc + 2, j]), int(mstate[loc + 3, j]))
            for j in range(L)))
        labels[s] = keys.setdefault(key, len(keys))
    return labels


def weed_fortran(hnorm, mstinf, numsta, eps=1e-4, max_iter=2000,
                 blocked=True, labels=None, return_kept=False):
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
        if labels is None:
            labels = _labels_from_matrix(A)
        keep = []
        for b in np.unique(labels[:numsta]):
            members = np.flatnonzero(labels[:numsta] == b)
            sub = _weed_indices(A[np.ix_(members, members)], eps, max_iter)
            keep.extend(int(members[i]) for i in sub)
        keep.sort()
    else:
        keep = _weed_indices(A, eps, max_iter)

    keep = np.asarray(keep, dtype=int)
    out = (A[np.ix_(keep, keep)].copy(), mstinf[keep], int(keep.size))
    # `return_kept` lets the caller carry block labels across the deletion
    # instead of rediscovering them from the weeded matrix.
    return out + (keep,) if return_kept else out


def free_energy_diagonal(mstate, mstinf, numsta, NF, cbreak):
    """``D`` such that ``H0 = D N``, straight from the parton content.

    The free part is diagonal in parton content, so ``H0`` and the norm differ
    only by a factor that depends on the state's momenta and flavours:

        D_s(f) = sum over partons p of flavour f of (2 + sgn_p * cbreak) / k_p

    with ``sgn = +1`` for a quark and ``-1`` for an antiquark -- read off the
    H8/H9 self-energy vertices in ``qcdf_kernels.hamqcd_nb``, which contribute
    ``cf * (2 + sgn*cbreak) / k`` to ``elem0``.  Verified against
    ``diag(H0)/diag(N)`` to 3.6e-15 at N=2/3/4, mesons and baryons.

    Why this matters: a configuration block is *defined* by the multiset of
    ``(type, momentum, flavour)``, and D depends on nothing else, so D is
    constant on a block **by construction** rather than by measurement.  With a
    blockwise Z every column then satisfies ``Z^T H0 Z = D I`` on its block, so
    the whole projection is ``diag(D)`` -- checked at 1e-14 off-diagonal.

    So ``ham0`` never has to be built or projected: two ``n_post**2`` arrays
    that were 1.86 GB each at 2K=41, for a length-``n`` vector.
    """
    D = np.zeros((numsta, max(NF, 1)))
    for s in range(numsta):
        loc = int(mstinf[s, 0]) - 1
        L = int(mstinf[s, 1])
        for j in range(L):
            typ = int(mstate[loc, j])
            k = int(mstate[loc + 2, j])
            f = int(mstate[loc + 3, j]) - 1
            if k == 0 or not (0 <= f < D.shape[1]):
                continue
            D[s, f] += (2.0 + (1.0 if typ == 1 else -1.0) * cbreak) / k
    return D


def weed_and_orthonormalize_blocks(blocks, numsta_pre, eps=1e-4, tol=1e-8):
    """Weed and orthonormalize straight from the norm blocks.

    The same arithmetic as ``weed_fortran(blocked=True)`` followed by
    :func:`orthonormalize_blockwise`, with the dense ``numsta_pre**2`` norm
    never allocated.  That array is the binding memory constraint of the solver
    -- 8.6 GB at 2K=37 and 51.5 GB at 2K=41, against tens of MB for the blocks
    -- and it is mostly zeros, since the norm is exactly block diagonal in
    parton configuration.

    Both stages already decompose over blocks, which is what makes this a pure
    change of container: ``WEEDR`` only ever compares rows inside one block, and
    ``WEEDR2``'s null space is the direct sum of the blocks' null spaces.  The
    per-block matrices are bit-identical to slicing the dense build, so the
    surviving indices and every column of ``Z`` are unchanged.

    Column order matches :func:`orthonormalize_blockwise` exactly -- blocks in
    ascending label, members ascending within a block -- so ``Z`` comes out
    column for column the same, not merely spanning the same space.

    ``Z`` comes back as CSR.  It is provably block diagonal here, so it is
    almost entirely zeros -- 32,135 non-zeros in a 6915 x 6915 array at 2K=37,
    0.39 MB against 0.38 GB dense -- and the triple product ``Z^T H Z`` is 11x
    faster sparse (0.94 s against 10.69 s).  That reorders the summation, so the
    spectrum moves by ~1e-14; the algorithm's own reproducibility floor is 1e-4
    (docs/basis-dependence.md), so this is four orders below the noise it lives
    in, but it does mean the result is *not* bit-identical to a dense triple
    product, only to everything upstream of it.

    Returns ``(Z, kept, weeded_blocks)``.  ``kept`` indexes the original states,
    ascending; ``weeded_blocks`` is ``[(members_in_kept_indexing, block), ...]``
    so a caller can assemble the weeded norm, or not.
    """
    from scipy.sparse import csr_matrix
    keep_parts, z_parts, weeded = [], [], []
    for members, M in blocks:
        sub = _weed_indices(M, eps)
        if not sub:
            continue
        members_kept = np.asarray([int(members[i]) for i in sub])
        Mw = M[np.ix_(sub, sub)]
        keep_parts.append(members_kept)
        weeded.append((members_kept, Mw))

        w, v = eigh(Mw)
        if w.min() < -tol:
            raise ValueError(
                f"norm block is not positive semidefinite: most negative "
                f"eigenvalue {w.min():.6g}. A Gram matrix cannot do this, so "
                f"the build is wrong rather than ill-conditioned -- see "
                f"docs/fortran-color-overflow.md.")
        sel = w > tol
        if np.any(sel):
            z_parts.append((members_kept, v[:, sel] / np.sqrt(w[sel])))

    if not keep_parts:
        return csr_matrix((0, 0)), np.zeros(0, dtype=int), []

    kept = np.sort(np.concatenate(keep_parts))
    pos = -np.ones(numsta_pre, dtype=np.int64)
    pos[kept] = np.arange(kept.size)

    # Assembled straight into COO triplets: the dense intermediate is the very
    # object being avoided (1.86 GB at 2K=41 for 90k non-zeros).
    ncol = sum(Zb.shape[1] for _, Zb in z_parts)
    rows_i, cols_i, vals = [], [], []
    c0 = 0
    for members_kept, Zb in z_parts:
        r = pos[members_kept]
        rr, cc = np.meshgrid(r, np.arange(c0, c0 + Zb.shape[1]), indexing="ij")
        rows_i.append(rr.ravel())
        cols_i.append(cc.ravel())
        vals.append(Zb.ravel())
        c0 += Zb.shape[1]

    Z = csr_matrix((np.concatenate(vals),
                    (np.concatenate(rows_i), np.concatenate(cols_i))),
                   shape=(kept.size, ncol))
    weeded = [(pos[m], M) for m, M in weeded]
    return Z, kept, weeded


#: Eigenpairs kept when a count is requested but not given.
#:
#: Not 30.  ``docs/next-steps.md`` says "the figures use at most 30 eigenvalues
#: and 11 eigenvectors", and both are true but they are two *different* thirties:
#: ``figures.figure2`` slices ``eigenvalues[:30]`` raw, while ``figures.figure5``
#: wants ``physical_indices(r)[10]`` -- the 11th *physical* level, which needs
#: 11 + however many spurious modes sit below it.  30 has no margin under either
#: reading, and the two readings disagree.  40 costs 7.4 MB of ``c_orig`` at
#: 2K=41 and covers both.
DEFAULT_NEV = 40

#: Extra eigenpairs requested from the iterative solver, then filtered.
#:
#: ``which="SA"`` returns the algebraically smallest, and non-positive M^2
#: artifacts live at the very bottom -- measured, N=4 B=1 2K=20 on the shipped
#: blockwise/exact path returns M^2 = -41.93 as its lowest eigenvalue, despite
#: observables.py having claimed the Python solver was immune.  Those would eat
#: the slots the caller asked for.
_NEV_SLACK = 8


def _solve_eigen(hnu, solver="dense", nev=None, tol=1e-8):
    """``(w, z)`` -- eigenvalues ascending, eigenvectors in columns.

    ``nev=None`` on the dense path keeps the whole spectrum, which is the
    historical behaviour.  ``nev=k`` keeps the lowest k however they were
    obtained, so a *dense* run with ``nev`` set is the shape-exact reference for
    a sparse one and the two can be compared without slicing.

    Numerics, all measured rather than assumed (see docs/next-steps.md):

    * ``which="SA"``, no ``sigma``, no folding.  Shift-invert would cut the
      iteration count 394 -> 160, but ``splu`` on this matrix fills in to ~83%
      of dense -- O(n^3) time and O(n^2) memory to replace a 3.6 s solve.
      Folding squares the condition number for a window already at the bottom
      of the spectrum, where ``SA`` is naturally efficient.
    * ``tol=1e-8`` delivers ~7e-15 relative -- machine precision, far past the
      1e-10 asked of it -- at 10-20% fewer matvecs than ``tol=0``.
    * ``ncv`` left to scipy's ``2k+1``: 40 is much worse (783 matvecs against
      426) and 150-250 cost more wall time for nothing.
    * A **fixed** ``v0``.  ``eigsh`` defaults to a random start vector, so two
      runs of the same parameters return different bits; combined with a cache
      keyed on those parameters, whichever landed first would win silently.
      Determinism here is correctness, not tidiness.

    The iteration count barely grows with n -- 220 at n=193, 432 at n=15235 --
    so the cost is O(nnz) per solve with a flat ~400-iteration constant.  That
    is why this pays: on a sparse ``hnu`` at 2K=41 it is 3.4 s against 123 s for
    a dense ``eigh(subset_by_index)``.
    """
    n = hnu.shape[0]
    if solver not in ("dense", "sparse"):
        raise ValueError(f"unknown solver {solver!r}; use 'dense' or 'sparse'")

    if solver == "dense":
        w, z = eigh(np.asarray(hnu.todense()) if sparse.issparse(hnu) else hnu)
        if nev is not None:
            w, z = w[:nev], z[:, :nev]
        return w, z

    from scipy.sparse.linalg import eigsh

    k = (DEFAULT_NEV if nev is None else nev) + _NEV_SLACK
    if k >= n - 1:
        # ARPACK needs k < n; below that threshold dense is faster anyway.
        w, z = eigh(np.asarray(hnu.todense()) if sparse.issparse(hnu) else hnu)
        keep = n if nev is None else min(nev, n)
        return w[:keep], z[:, :keep]

    v0 = np.full(n, 1.0 / np.sqrt(n))
    w, z = eigsh(hnu, k=k, which="SA", v0=v0, tol=tol)
    order = np.argsort(w)
    w, z = w[order], z[:, order]
    keep = k if nev is None else min(nev, k)
    return w[:keep], z[:, :keep]


def _require_positive_norm(w, N, B, K_code, tol=1e-8):
    """Fail clearly if the weeded norm is not positive definite.

    ``Z = z / sqrt(w)`` is where a negative eigenvalue first does damage, and
    it does it quietly: ``np.sqrt`` of a negative returns ``nan`` with only a
    ``RuntimeWarning``, and the failure then surfaces several steps later as
    ``ValueError: array must not contain infs or NaNs`` from inside ``eigh``,
    pointing at the wrong place entirely.

    A weeded norm cannot legitimately have a negative eigenvalue: it is a Gram
    matrix, and weeding has already removed everything with ``|w| < 1e-4``.  So
    this only fires on a real defect -- it is what silently truncated colour
    sums look like by the time they reach here (see
    docs/fortran-color-overflow.md, "The Python port's own overflow").
    """
    bad = np.flatnonzero(w < -tol)
    if bad.size:
        raise ValueError(
            f"norm matrix is not positive definite at N={N} B={B} "
            f"2K={K_code}: {bad.size} eigenvalue(s) below {-tol:g}, most "
            f"negative {w.min():.6g}. A Gram matrix cannot do this, so the "
            f"matrix build is wrong rather than merely ill-conditioned -- see "
            f"docs/fortran-color-overflow.md.")


def orthonormalize_blockwise(hnorm, numsta, tol=1e-8, labels=None):
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
    if labels is None:
        labels = _labels_from_matrix(np.asarray(hnorm)[:numsta, :numsta])

    cols = []
    worst = 0.0
    for b in np.unique(labels[:numsta]):
        idx = np.flatnonzero(labels[:numsta] == b)
        w, v = eigh(hnorm[np.ix_(idx, idx)])
        keep = w > tol
        worst = min(worst, float(w.min()))
        if not np.any(keep):
            continue
        Zb = v[:, keep] / np.sqrt(w[keep])
        for j in range(Zb.shape[1]):
            col = np.zeros(numsta)
            col[idx] = Zb[:, j]
            cols.append(col)

    # A column short of ``numsta`` means weeding and this step disagreed about
    # the rank.  They test different things -- weeding drops on ``|w| < 1e-4``,
    # matching Fortran WEEDR2, while this keeps ``w > tol`` -- and the one
    # eigenvalue class they disagree on is *negative with magnitude above
    # 1e-4*: weeding keeps such a state, this discards its direction.
    #
    # That is not a threshold to be tuned.  A weeded Gram matrix has no such
    # eigenvalue, so its appearance means the matrix is wrong.  Measured: it
    # showed up as n_post=10332 against n_orth=10328 at N=3, B=1, 2K=39, which
    # is exactly where the colour sums begin overrunning MXTRM.  Silently
    # returning a narrower Z there leaves every consumer that assumes
    # ``n_orth == n_post`` -- hnu's shape, c_orig's columns, validate() --
    # quietly wrong.
    if worst < -tol:
        raise ValueError(
            f"norm block matrix is not positive semidefinite: most negative "
            f"eigenvalue {worst:.6g} (tol {-tol:g}). A Gram matrix cannot do "
            f"this, so the build is wrong rather than ill-conditioned -- see "
            f"docs/fortran-color-overflow.md.")

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
               assembly="exact", prefer_opt=True, backend=None,
               solver="dense", nev=None, keep_norm=True,
               hamiltonian="standard") -> DLCQResult:
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

    hamiltonian : {"standard", "improved"}
        ``"standard"`` (default) is the reference path and is bit-identical to
        what this code has always produced.

        ``"improved"`` applies van de Sande's endpoint subtraction
        (`hep-ph/9605409 <https://arxiv.org/abs/hep-ph/9605409>`_) to the
        self-induced inertia, which removes the ``K^-2b`` error that makes the
        weak-coupling limit unreachable -- see ``docs/weak-coupling-limit.md``
        and :mod:`dlcq.endpoint`.  Measured against his exact ``M^2/g^2 =
        0.779141``: standard DLCQ reaches 0.413 at K = 100, improved reaches
        0.771 at K = **10**.

        **Two-parton sectors only**, and it raises otherwise.  There the
        partner momentum is fixed by conservation, so the subtraction stays a
        function of one momentum and is applied by swapping the ``selfen``
        table.  With three or more partons it depends on *which partner* a
        parton is paired with, which a table indexed by one momentum cannot
        express -- a limit of the representation, not of the physics.  Use
        ``LPN=2`` for a meson.

    ``rlamb`` is taken literally.  Do not derive it from m/g when comparing
    against a Fortran run: the input files use 0.3325 while
    ``mg_to_lambda(1.6)`` is 0.3325495, and that 1.5e-5 difference moves
    eigenvalues well above a 1e-8 comparison tolerance.
    """
    # Validate before importing the solver or generating states: an unknown
    # policy used to surface only after the norm build, which at large K is
    # minutes of work to learn about a typo.
    if assembly not in ("exact", "fortran"):
        raise ValueError(f"unknown assembly {assembly!r}")
    if policy not in ("fortran", "blockwise", "spectral"):
        raise ValueError(f"unknown weeding policy {policy!r}")
    if solver not in ("dense", "sparse"):
        raise ValueError(f"unknown solver {solver!r}; use 'dense' or 'sparse'")
    if nev is not None and nev < 1:
        raise ValueError(f"nev must be at least 1, got {nev}")
    if hamiltonian not in ("standard", "improved"):
        raise ValueError(
            f"unknown hamiltonian {hamiltonian!r}; use 'standard' or "
            f"'improved' (dlcq/endpoint.py, docs/weak-coupling-limit.md)")
    if solver == "sparse" and (assembly, policy) != ("exact", "blockwise"):
        raise ValueError(
            f"solver='sparse' supports assembly='exact' with "
            f"policy='blockwise' only; got assembly={assembly!r}, "
            f"policy={policy!r}. The Fortran-compatible diagonal-only assembly "
            f"is inherently dense and basis dependent "
            f"(docs/basis-dependence.md), and a global Z is dense by "
            f"construction, so neither can reach the sizes this path exists "
            f"for.")
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

    # qcdsta_fast produces element-identical mstate/mstinf (tests/test_states.py)
    # and grows its own tables, so it is not bound by LKPXMX/NMSTMX.
    if use_opt:
        opt.qcdsta_fast(p, states, perm, flav)
    else:
        base.qcdsta(p, states, perm, flav)
    numsta_pre = states.numsta
    if numsta_pre == 0:
        return DLCQResult(N=N, NF=NF, B=B, K_code=K_code, rlamb=rlamb,
                          cutoff=cutoff, LPN=LPN, source="python",
                          numsta_pre=0, numsta_post=0,
                          eigenvalues=np.array([]))

    mstinf = states.mstinf[:numsta_pre].copy()

    if hamiltonian == "improved":
        # Applied by swapping the ``selfen`` table rather than touching the
        # kernel, which only works while the partner momentum is fixed by
        # conservation -- i.e. two partons.  Checked here rather than assumed,
        # because LPN=0 means *no* truncation, not valence only, so an
        # innocent-looking call lands in a mixed-Fock basis.
        parton_counts = sorted({int(mstinf[s, 1]) for s in range(numsta_pre)})
        if parton_counts != [2]:
            raise ValueError(
                f"hamiltonian='improved' supports two-parton states only; this "
                f"basis has parton counts {parton_counts}. Use LPN=2 for a "
                f"meson. Beyond two partons the subtraction depends on which "
                f"partner a parton is paired with, which a selfen table indexed "
                f"by a single momentum cannot express; it needs a partner-"
                f"indexed table instead. See dlcq/endpoint.py and "
                f"docs/next-steps.md.")
        from .endpoint import improved_selfen
        from .units import lambda_to_mg, endpoint_exponent
        selfen = improved_selfen(N, K_code,
                                 endpoint_exponent(lambda_to_mg(rlamb), N),
                                 mxslfn=len(selfen))

    # Block labels from the state list, not by thresholding the matrix: the
    # latter costs two more dense n x n temporaries on top of a norm that is
    # already 8.6 GB at 2K=37.  See config_block_labels.
    labels = config_block_labels(states.mstate, mstinf, numsta_pre)

    # ── norm matrix ──
    #
    # ``policy="blockwise"`` never forms the pre-weeding norm.  It is the
    # binding memory constraint of the solver -- 51.5 GB at 2K=41 against
    # 27.9 MB of blocks -- and it is mostly zeros, the norm being exactly block
    # diagonal in parton configuration.  The blocks are bit-identical to
    # slicing the dense build, so this is a change of container only.
    blocks = None
    hnorm = None
    used_blocks = False
    # getattr, not a bare call: the interpreted reference solver has no block
    # builder, and a caller may be running against an older qcdf_opt.  Falling
    # back to the dense build is always correct, just heavier.
    _blocks_builder = getattr(opt, "build_norm_blocks", None) if use_opt else None
    if policy == "blockwise" and callable(_blocks_builder):
        blocks = _blocks_builder(states.mstate, mstinf, numsta_pre,
                                 N, NF, B, K_code, labels, nthreads=ncpus)
    elif use_opt:
        _, _, hnorm = opt.build_matrices(0, states.mstate, states.mstinf,
                                         numsta_pre, N, NF, B, K_code,
                                         selfen, p.cbreak, ncpus,
                                         backend=backend)
    else:
        _, _, hnorm = base.clrdis(0, p, states, selfen, ncpus=ncpus)

    # ── weed ──
    if policy == "fortran":
        hnorm_w, mstinf_w, numsta, kept = weed_fortran(
            hnorm, mstinf, numsta_pre, labels=labels, return_kept=True)
        w_n, z_n = eigh(hnorm_w[:numsta, :numsta])
        _require_positive_norm(w_n, N, B, K_code)
        Z = z_n / np.sqrt(w_n)[np.newaxis, :]          # NUZ
    elif policy == "blockwise" and blocks is not None:
        Z, kept, weeded_blocks = weed_and_orthonormalize_blocks(
            blocks, numsta_pre)
        numsta = int(kept.size)
        mstinf_w = mstinf[kept]
        used_blocks = True
        del blocks
        # The *weeded* norm is what DLCQResult carries, and it is far smaller
        # than the pre-weeding one (n_post ~ 0.3 n_pre at LPN=0).  Assembling it
        # densely keeps the data contract unchanged; making it block-sparse is
        # the next step, and is what removes the remaining n_post**2.
        if keep_norm:
            # Kept as blocks, not assembled.  The dense weeded norm was the last
            # n^2 object in the pipeline -- 151 GB at LPN=7, 2K=71 -- and it is
            # used only as the matvec ``norm @ c`` by structure_function, which
            # BlockDiagonal provides.  So the structure functions (Figs. 4-6,
            # and the unresolved Fig. 6a) stay available at every K the solver
            # can now reach, rather than being traded away for the memory.
            hnorm_w = BlockDiagonal(
                n_rows=numsta, n_cols=numsta,
                rows=[np.asarray(r) for r, _ in weeded_blocks],
                blocks=[M for _, M in weeded_blocks])
        else:
            hnorm_w = None
        del weeded_blocks
    elif policy == "blockwise":
        hnorm_w, mstinf_w, numsta, kept = weed_fortran(
            hnorm, mstinf, numsta_pre, labels=labels, return_kept=True)
        Z, _ = orthonormalize_blockwise(hnorm_w, numsta,
                                        labels=labels[kept])
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
    # H0 = D*N with D fixed by parton content, and a blockwise Z keeps every
    # column inside one configuration block, so Z^T H0 Z is exactly diag(D) --
    # verified to 1e-14 off-diagonal, 1e-13 on the diagonal.  Building ham0 and
    # projecting it is then pure waste: two n_post**2 arrays, 1.86 GB each at
    # 2K=41.  Only on the blockwise path, where the diagonality argument holds.
    use_h0_shortcut = used_blocks
    # H is the next n_post**2 array after the norm: 1.86 GB at 2K=41 for 6.85 M
    # non-zeros, 2.95% dense.  Density falls 26% -> 3% over 2K = 21 -> 41
    # because n grows 1.58x per +2 while nnz/row grows only 1.25x.  The CSR
    # build stores bit-identical values -- same kernel, same traversal, only the
    # zeros dropped -- and keeping hnu sparse is also what makes the iterative
    # solver worth using: on a dense operator every ARPACK matvec is an O(n^2)
    # GEMV, which is why solver="sparse" was no faster before this.
    _sparse_ham = getattr(opt, "build_ham_sparse", None) if use_opt else None
    use_sparse_ham = used_blocks and callable(_sparse_ham)
    ham0 = None
    if use_sparse_ham:
        ham = _sparse_ham(states.mstate, mstinf_w[:numsta], numsta, N, NF, B,
                          K_code, selfen, p.cbreak, nthreads=ncpus)
    elif use_opt:
        _bm_kw = {} if not use_h0_shortcut else {"want_h0": False}
        ham0, ham, _ = opt.build_matrices(1, states.mstate, mstinf_w[:numsta],
                                          numsta, N, NF, B, K_code,
                                          selfen, p.cbreak, ncpus,
                                          backend=backend, **_bm_kw)
        ham = ham[:numsta, :numsta]
    else:
        use_h0_shortcut = False
        ham0, ham, _ = base.clrdis(1, p, states, selfen, ncpus=ncpus)
        ham = ham[:numsta, :numsta]

    # ── NUHAM, then combine free + interacting ──
    #
    # ``Z`` is CSR on the blockwise path, where it is provably block diagonal
    # and so almost all zeros.  ``(Z^T (Z^T H)^T)`` keeps both products sparse x
    # dense and never forms a dense Z: 0.94 s against 10.69 s at 2K=37, and it
    # is what stops the 1.86 GB dense Z appearing at 2K=41.
    def _project(M):
        if sparse.issparse(Z) and sparse.issparse(M):
            # Both sparse: the result stays sparse, which is the point.  hnu is
            # only 4% denser than H (3.08% against 2.95% at 2K=41), so the
            # triple product does not fill in -- measured, not assumed.
            return (Z.T @ M @ Z).tocsr()
        if sparse.issparse(Z):
            return np.asarray(Z.T @ (Z.T @ M).T)
        return Z.T @ M @ Z

    hnu = _project(ham)
    del ham

    if use_h0_shortcut:
        # diag(D) per column, D read off the block each column occupies.
        Dstate = free_energy_diagonal(states.mstate, mstinf_w[:numsta], numsta,
                                      NF, p.cbreak)
        Zc = Z.tocsc()
        first = Zc.indices[Zc.indptr[:-1]]       # a representative row per column
        free_diag = [Dstate[first, ifl] for ifl in range(NF)]
        free = None
    else:
        free = [_project(ham0[ifl, :numsta, :numsta]) for ifl in range(NF)]
        free_diag = [np.diag(f) for f in free]

    rlmsq = rlamb ** 2
    hnu = hnu * rlmsq

    if sparse.issparse(hnu):
        # Add the free part as a diagonal rather than densifying: a dia_matrix
        # of the accumulated D keeps hnu at ~3% density all the way to the
        # eigensolver, which is the only reason the iterative path pays.
        acc = np.zeros(hnu.shape[0])
        for ifl in range(NF):
            acc += (1 - rlmsq) * p.rmq[ifl] ** 2 * free_diag[ifl]
        hnu = (hnu + sparse.diags(acc)).tocsr()
        free = None
        NF_loop = range(0)
    else:
        NF_loop = range(NF)

    for ifl in NF_loop:
        rmq_val = p.rmq[ifl]
        coeff = (1 - rlmsq) * rmq_val ** 2
        if assembly == "exact" and free is not None:
            hnu = hnu + coeff * free[ifl]
        else:
            # Diagonal only -- exact here rather than an approximation, because
            # Z^T H0 Z *is* diagonal in a blockwise basis (this is the whole
            # point of docs/basis-dependence.md), and it is what qcdf.f does
            # unconditionally under assembly="fortran".
            hnu[np.diag_indices(hnu.shape[0])] += coeff * free_diag[ifl]

    # Recorded before truncation: it is the length a *full* spectrum would have,
    # and it is what lets a consumer tell "this run has 40 levels" from "this
    # run computed 40 of 15235".
    n_orth = hnu.shape[0]
    w_eig, z_eig = _solve_eigen(hnu, solver=solver, nev=nev)
    eigenvalues = K_code * w_eig / 2.0
    c_orig = np.asarray(Z @ z_eig)
    # DLCQResult's contract is ndarray, and every consumer -- structure_function,
    # the figures, dataset.save -- assumes that.  Densifying Z here confines the
    # sparse form to the arithmetic above, where it earns its keep; making the
    # stored form block-sparse is a schema change and belongs with Stage 2.
    if sparse.issparse(Z):
        Z = Z.toarray()

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
        # "solver" here has always meant the engine module, so the new flag
        # gets its own key rather than colliding with it.
        provenance={"solver": "qcdf_opt" if use_opt else "qcdf",
                    "weeding_policy": policy, "assembly": assembly,
                    "ncpus": ncpus,
                    "eigensolver": solver, "nev": nev,
                    "hamiltonian": hamiltonian,
                    "backend": backend or os.environ.get("QCDF_BACKEND")
                    or getattr(opt, "DEFAULT_BACKEND", "process")},
        numsta_pre=numsta_pre, numsta_post=numsta, n_orth=int(n_orth),
        state_len=lengths, state_types=types, state_moms=moms,
        state_flavors=flavors,
        norm=hnorm_w, Z=Z,
        eigenvalues=eigenvalues, eigenvectors=z_eig, c_orig=c_orig,
    )
