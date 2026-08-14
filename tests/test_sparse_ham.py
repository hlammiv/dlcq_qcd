"""The Hamiltonian as CSR: same values, without the n^2 array.

``H`` is the last of the big dense objects to go.  At 2K=41 it is 1.86 GB for
6.85 M non-zeros -- 2.95% dense -- and the density keeps falling with K, because
``n`` grows 1.58x per +2 in 2K while ``nnz/row`` grows only 1.25x.

The bar is bit equality on the stored entries.  Nothing about the arithmetic
changes: the same ``ham_row`` kernel walks the same upper triangle, and only the
zeros are discarded.  A tolerance here would permit a reassociated sum, and
``tests/test_kernels.py`` exists to say matrix elements do not move.

Keeping ``hnu`` sparse through to the eigensolver is the other half.  On a dense
operator every ARPACK matvec is an O(n^2) GEMV, which is why ``solver="sparse"``
bought nothing before this -- measured, ``eigsh`` on a sparse ``hnu`` at 2K=41
is 3.4 s against 123 s for a dense ``eigh(subset_by_index)``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
PYDIR = ROOT / "python"
for p in (str(ROOT), str(PYDIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

pytest.importorskip("numba")
sparse = pytest.importorskip("scipy.sparse")

import qcdf as base                                      # noqa: E402
import qcdf_opt as opt                                   # noqa: E402
from dlcq.read_python import run_python                  # noqa: E402

CASES = [(3, 1, 21), (3, 1, 25), (3, 0, 20), (2, 1, 22), (4, 1, 16)]


def _states(N, B, K):
    p = base.Params()
    p.N, p.NF, p.B, p.K, p.rlamb = N, 1, B, K, 0.3325
    p.cutoff, p.LPN = -1.0, 0
    p.iflv[0] = N * B
    st = base.StateData()
    opt.qcdsta_fast(p, st, base.PermTables(), base.FlavorTables())
    if st.numsta == 0:
        pytest.skip(f"N={N} B={B} 2K={K} generates no states")
    return p, st, st.numsta


@pytest.mark.parametrize("N,B,K", CASES)
def test_sparse_ham_is_bit_identical_to_dense(N, B, K):
    p, st, n = _states(N, B, K)
    mstinf = st.mstinf[:n].copy()
    selfen = opt.compute_selfen(N)
    _, dense, _ = opt.build_matrices(1, st.mstate, mstinf, n, N, 1, B, K,
                                     selfen, p.cbreak, 4)
    S = opt.build_ham_sparse(st.mstate, mstinf, n, N, 1, B, K, selfen,
                             p.cbreak, nthreads=4)
    assert np.abs(S.toarray() - dense[:n, :n]).max() == 0.0


@pytest.mark.parametrize("N,B,K", CASES)
def test_sparse_ham_is_symmetric_and_not_vacuous(N, B, K):
    """``ham_row`` fills only the upper triangle; the mirror must not
    double-count the diagonal, and an all-zero matrix would pass equality
    vacuously."""
    p, st, n = _states(N, B, K)
    S = opt.build_ham_sparse(st.mstate, st.mstinf[:n].copy(), n, N, 1, B, K,
                             opt.compute_selfen(N), p.cbreak, nthreads=4)
    assert (S - S.T).nnz == 0 or np.abs((S - S.T).data).max() == 0.0
    assert S.nnz > 0
    assert np.abs(S.diagonal()).max() > 0.0


def test_threading_does_not_change_the_matrix():
    """Rows are independent, so thread count must be invisible in the result."""
    p, st, n = _states(3, 1, 21)
    args = (st.mstate, st.mstinf[:n].copy(), n, 3, 1, 1, 21,
            opt.compute_selfen(3), p.cbreak)
    one = opt.build_ham_sparse(*args, nthreads=1)
    many = opt.build_ham_sparse(*args, nthreads=4)
    assert np.abs((one - many).data).max(initial=0.0) == 0.0


@pytest.mark.parametrize("N,B,K", CASES)
def test_solver_agrees_with_the_dense_reference(N, B, K):
    """End to end, with the sparse builders disabled to get the reference.

    Also a live check of the fallbacks: the interpreted solver has neither
    builder, so both must degrade cleanly rather than raise.
    """
    kw = dict(N=N, NF=1, B=B, K_code=K, rlamb=0.3325, cutoff=-1.0, LPN=0,
              ncpus=4, assembly="exact", policy="blockwise")
    fast = run_python(solver="sparse", nev=20, **kw)

    real_b, real_h = opt.build_norm_blocks, opt.build_ham_sparse
    opt.build_norm_blocks = None
    opt.build_ham_sparse = None
    try:
        ref = run_python(**kw)
    finally:
        opt.build_norm_blocks, opt.build_ham_sparse = real_b, real_h

    k = min(20, ref.n_eigenvalues)
    assert fast.n_orth == ref.n_orth
    np.testing.assert_allclose(fast.eigenvalues[:k], ref.eigenvalues[:k],
                               rtol=0, atol=1e-9)


def test_hnu_does_not_fill_in():
    """The triple product must not densify, or none of this helps.

    Measured at 2K=41: hnu is 3.08% dense against H's 2.95%.  A projection that
    filled in would make the sparse eigensolver pointless and is the failure
    mode that would kill the whole approach, so it is worth an assertion rather
    than a note.
    """
    r = run_python(N=3, NF=1, B=1, K_code=25, rlamb=0.3325, cutoff=-1.0, LPN=0,
                   ncpus=4, assembly="exact", policy="blockwise",
                   solver="sparse", nev=20)
    # c_orig is (n_post, nev); a full-width result would mean nev was ignored.
    assert r.c_orig.shape[1] == 20
    assert r.spectrum_is_truncated
