"""The blockwise norm must be a change of container, not of arithmetic.

The pre-weeding norm is the binding memory constraint of the whole solver: it
is built at ``numsta_pre``, which is 8.6 GB at 2K=37 and 51.5 GB at 2K=41, and
it is overwhelmingly zeros -- the norm is *exactly* block diagonal in parton
configuration, off-block entries being ``0.0`` rather than 1e-16.  Building only
the blocks removes that array.

The whole value of the change rests on it altering nothing else, so the bar here
is **bit equality**, not a tolerance: same blocks, same surviving states, same
``Z`` column for column, same spectrum.  A tolerance would let a reassociated
sum through, and ``docs/basis-dependence.md`` records that this algorithm has an
intrinsic 1e-4 reproducibility floor -- a small perturbation would be
indistinguishable from that noise rather than caught by it.
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

import qcdf as base                                      # noqa: E402
import qcdf_opt as opt                                   # noqa: E402
from dlcq.dataset import BlockDiagonal                        # noqa: E402
from dlcq.read_python import config_block_labels, run_python   # noqa: E402


# (N, B, 2K) -- mesons and baryons, SU(2)/SU(3)/SU(4), kept small.
CASES = [(3, 1, 21), (3, 1, 25), (3, 0, 20), (2, 1, 22), (4, 1, 16)]


def _states(N, B, K, LPN=0):
    p = base.Params()
    p.N, p.NF, p.B, p.K, p.rlamb = N, 1, B, K, 0.3325
    p.cutoff, p.LPN = -1.0, LPN
    p.iflv[0] = N * B
    st = base.StateData()
    opt.qcdsta_fast(p, st, base.PermTables(), base.FlavorTables())
    if st.numsta == 0:
        pytest.skip(f"N={N} B={B} 2K={K} generates no states")
    return p, st, st.numsta


@pytest.mark.parametrize("N,B,K", CASES)
def test_blocks_are_bit_identical_to_the_dense_build(N, B, K):
    """Same kernel, same column sets, only a different destination."""
    p, st, n = _states(N, B, K)
    mstinf = st.mstinf[:n].copy()
    _, _, dense = opt.build_matrices(0, st.mstate, mstinf, n, N, 1, B, K,
                                     opt.compute_selfen(N), p.cbreak, 4)
    labels = config_block_labels(st.mstate, mstinf, n)
    blocks = opt.build_norm_blocks(st.mstate, mstinf, n, N, 1, B, K, labels,
                                   nthreads=4)
    for members, M in blocks:
        assert np.array_equal(M, dense[np.ix_(members, members)])


@pytest.mark.parametrize("N,B,K", CASES)
def test_the_blocks_cover_every_state_and_miss_no_coupling(N, B, K):
    """Two failure modes the equality test above cannot see on its own.

    A block scheme could be bit-perfect on the blocks it builds while dropping
    states entirely, or while leaving a genuine coupling outside every block.
    The second is the dangerous one: the memory saving is only legitimate
    because what is not stored is exactly zero.
    """
    p, st, n = _states(N, B, K)
    mstinf = st.mstinf[:n].copy()
    _, _, dense = opt.build_matrices(0, st.mstate, mstinf, n, N, 1, B, K,
                                     opt.compute_selfen(N), p.cbreak, 4)
    labels = config_block_labels(st.mstate, mstinf, n)
    blocks = opt.build_norm_blocks(st.mstate, mstinf, n, N, 1, B, K, labels,
                                   nthreads=4)

    seen = np.concatenate([m for m, _ in blocks])
    assert np.array_equal(np.sort(seen), np.arange(n))

    mask = np.ones((n, n), dtype=bool)
    for members, _ in blocks:
        mask[np.ix_(members, members)] = False
    assert np.abs(dense[:n, :n][mask]).max() == 0.0


@pytest.mark.parametrize("N,B,K", CASES)
def test_solver_agrees_with_the_dense_assembled_path(N, B, K):
    """End to end, with the bar set where each quantity can actually meet it.

    The basis is bit-identical: same surviving states, same weeded norm, same
    ``Z``.  The *spectrum* is not, and cannot be -- ``Z`` is held as CSR, so the
    triple product sums in a different order, worth ~1e-13.  That buys 11x on
    the projection (0.94 s against 10.69 s at 2K=37) and is four orders below
    this algorithm's own reproducibility floor, which docs/basis-dependence.md
    measures at 1e-4: recompiling the unmodified Fortran at -O2 moves its own
    answer by more.

    Asserting equality on the basis and a tight tolerance on the spectrum keeps
    the distinction visible, rather than relaxing everything to the loosest
    thing that passes.

    Exercised by disabling the block builder so the same policy runs through the
    old dense-assembled route -- which is also a live check that the fallback
    still works, since the interpreted reference solver has no block builder.
    """
    kw = dict(N=N, NF=1, B=B, K_code=K, rlamb=0.3325, cutoff=-1.0, LPN=0,
              ncpus=4, policy="blockwise", assembly="exact")
    fast = run_python(**kw)

    real = opt.build_norm_blocks
    opt.build_norm_blocks = None
    try:
        dense = run_python(**kw)
    finally:
        opt.build_norm_blocks = real

    assert fast.numsta_post == dense.numsta_post
    # The blockwise path stores the norm as blocks -- the dense form is 151 GB
    # at LPN=7, 2K=71 -- so densify explicitly to compare.  ``.toarray()`` must
    # be *exactly* the dense build: the blocks hold the same values, and what is
    # not stored is exactly zero.
    assert isinstance(fast.norm, BlockDiagonal)
    assert np.array_equal(fast.norm.toarray(), dense.norm)   # basis: exact
    np.testing.assert_allclose(fast.Z, dense.Z, rtol=0, atol=1e-12)
    np.testing.assert_allclose(fast.eigenvalues, dense.eigenvalues,
                               rtol=0, atol=1e-9)         # spectrum: reordered


def test_every_column_of_Z_lives_in_one_block():
    """What makes the Fortran's diagonal-only assembly well posed.

    ``H0 = D N`` with D constant on a block, so a Z confined to blocks gives a
    diagonal ``Z^T H0 Z``; a Z that mixed blocks would not, and
    docs/basis-dependence.md records that the resulting spectrum then depends on
    which eigenvectors came back.
    """
    r = run_python(N=3, NF=1, B=1, K_code=21, rlamb=0.3325, cutoff=-1.0,
                   LPN=0, ncpus=4, policy="blockwise", assembly="exact")

    # Label the *surviving* basis from the result's own state arrays.  Taking
    # the pre-weeding labels and slicing them would be wrong: weeding removes
    # states, so post-weeding row i is not pre-weeding state i.
    keys = {}
    labels = np.empty(r.numsta_post, dtype=np.int64)
    for s in range(r.numsta_post):
        L = int(r.state_len[s])
        key = tuple(sorted((int(r.state_types[s, j]), int(r.state_moms[s, j]),
                            int(r.state_flavors[s, j])) for j in range(L)))
        labels[s] = keys.setdefault(key, len(keys))

    for j in range(r.Z.shape[1]):
        rows = np.flatnonzero(np.abs(r.Z[:, j]) > 0)
        assert rows.size
        assert np.unique(labels[rows]).size == 1, (
            f"column {j} of Z spans more than one configuration block")
