"""The threaded numba kernels must reproduce the reference build exactly.

``qcdf_kernels.py`` is a hand transcription of ``clfact`` and ``hamqcd`` from
``qcdf_opt.py``, so nothing but a direct comparison establishes it is faithful.
These tests build both matrices both ways and demand **bit-identical** results,
not merely close ones: the kernels do the same arithmetic in the same order, so
any difference at all is a transcription error rather than rounding.

That matters more than usual here.  ``docs/basis-dependence.md`` records that
the diagonal-only assembly step is basis dependent at the 1e-4 level, so a
small perturbation of the matrices would be indistinguishable from the noise
the Fortran already has against itself.  Exact equality leaves no such gap.
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


# (N, NF, B, 2K, iflv) -- mesons and baryons, SU(2)/SU(3)/SU(4), one and two
# flavours.  Kept small so the interpreted reference side stays affordable.
CASES = [
    (3, 1, 1, 11, None),
    (3, 1, 1, 15, None),
    (3, 1, 0, 14, None),
    (2, 1, 1, 12, None),
    (2, 1, 0, 12, None),
    (4, 1, 1, 16, None),
    (3, 2, 0, 10, [0, 0]),
    (3, 2, 1, 11, [3, 0]),
]

# (4, 1, 1, 16) was silently wrong for a while: its Hamiltonian needs more than
# the 12552 colour terms ``MXTRM`` used to allow, and *both* backends carry the
# same limit, so they agreed bit-for-bit on a truncated sum and this file's
# equality assertions passed.  MXTRM is now 50000 and the values are exact
# against brute-force colour enumeration, so it is back to being an ordinary
# case -- but only tests/test_colour_overflow.py can tell the difference, since
# backend equality cannot see a defect both backends share.
#
# These still exceed the limit and must fail loudly rather than truncate.  Note
# it is the *Hamiltonian* that overruns here while the norm is fine, so parton
# count alone does not predict which runs are affected.
OVERFLOW_CASES = [
    (4, 1, 1, 20, None),
    (4, 1, 1, 24, None),
]


def _states(N, NF, B, K, iflv=None):
    p = base.Params()
    p.N, p.NF, p.B, p.K, p.rlamb = N, NF, B, K, 0.3325
    p.cutoff, p.LPN = -1.0, 0
    if iflv is None:
        iflv = [0] * NF
        iflv[0] = N * B
    for i, v in enumerate(iflv[:NF]):
        p.iflv[i] = v
    st = base.StateData()
    base.qcdsta(p, st, base.PermTables(), base.FlavorTables())
    return p, st


def _both(norh, N, NF, B, K, iflv):
    p, st = _states(N, NF, B, K, iflv)
    n = st.numsta
    if n == 0:
        pytest.skip(f"N={N} NF={NF} B={B} 2K={K} generates no states")
    selfen = opt.compute_selfen(N)
    mstinf = st.mstinf[:n].copy()
    args = (norh, st.mstate, mstinf, n, N, NF, B, K, selfen, p.cbreak)
    ref = opt.build_matrices(*args, 4, backend="process")
    got = opt.build_matrices(*args, 4, backend="thread")
    return ref, got, n


@pytest.mark.parametrize("N,NF,B,K,iflv", CASES)
def test_norm_matrix_identical(N, NF, B, K, iflv):
    """Threaded norm build equals the process build, element for element."""
    ref, got, n = _both(0, N, NF, B, K, iflv)
    assert np.array_equal(ref[2], got[2]), (
        f"norm differs, max |delta| = {np.max(np.abs(ref[2] - got[2])):.3e}")
    # A norm matrix with a zero diagonal would pass a comparison vacuously.
    assert np.count_nonzero(np.diag(got[2])) == n


@pytest.mark.parametrize("N,NF,B,K,iflv", CASES)
def test_hamiltonian_identical(N, NF, B, K, iflv):
    """Threaded Hamiltonian build equals the process build, including H0."""
    ref, got, _ = _both(1, N, NF, B, K, iflv)
    assert np.array_equal(ref[1], got[1]), (
        f"H differs, max |delta| = {np.max(np.abs(ref[1] - got[1])):.3e}")
    assert np.array_equal(ref[0], got[0]), (
        f"H0 differs, max |delta| = {np.max(np.abs(ref[0] - got[0])):.3e}")
    assert np.any(got[1] != 0.0)


@pytest.mark.parametrize("N,NF,B,K,iflv", OVERFLOW_CASES)
def test_colour_overflow_raises_in_both_backends(N, NF, B, K, iflv):
    """A configuration past MXTRM must fail, and fail the same way both ways.

    The failure has to be symmetric: ``qcdf_kernels.py`` and ``qcdf_opt.py``
    carry independent copies of MXTRM, and if they ever drift apart one backend
    would return a number where the other raises -- which is a far more
    confusing bug than the one this replaces.
    """
    p, st = _states(N, NF, B, K, iflv)
    n = st.numsta
    selfen = opt.compute_selfen(N)
    args = (1, st.mstate, st.mstinf[:n].copy(), n, N, NF, B, K, selfen, p.cbreak)

    outcomes = {}
    for backend in ("process", "thread"):
        try:
            opt.build_matrices(*args, 4, backend=backend)
            outcomes[backend] = "returned"
        except OverflowError:
            outcomes[backend] = "raised"

    assert outcomes["process"] == outcomes["thread"] == "raised", outcomes


def test_epsilon_table_is_the_signed_permutations():
    """``epsb_table`` reproduces what ``clfact`` used to rebuild per call.

    This is the SU(N) epsilon reduction BREDCE applies: the permutations of
    ``1..N-1`` with their signatures.  ``qcdf_opt.clfact`` derived it inside
    its innermost loop -- 15409 times in one profiled row batch -- for a table
    that depends only on N.  Hoisting it is the one algorithmic change in the
    kernels, so pin the result.
    """
    import qcdf_kernels as kern
    from qcdf import gprbig_array, psign, MXNP

    nprms3, ib3 = kern.epsb_table(3)
    assert nprms3 == 2
    assert ib3[:2, :3].tolist() == [[1, 2, 1], [2, 1, -1]]

    # N=2 has no extra terms; clfact's own branch set nprms_ep = 0.
    assert kern.epsb_table(2)[0] == 0

    for N in (3, 4, 5):
        nmin = np.ones(MXNP, dtype=np.int64)
        nmax = np.ones(MXNP, dtype=np.int64) * (N - 1)
        want_n, want = gprbig_array(nmin, nmax, N - 1, 1, 0, 9523)
        for i in range(want_n):
            want[i, N - 1] = psign(N - 1, want[i, :N - 1].tolist())
        got_n, got = kern.epsb_table(N)
        assert got_n == want_n
        assert np.array_equal(got[:got_n, :N], want[:want_n, :N])


@pytest.mark.parametrize("nthreads", [1, 2, 8, 20])
def test_thread_count_does_not_change_the_answer(nthreads):
    """Scratch buffers are per-thread, so results must not depend on the pool.

    ``qcdf_opt`` keeps one module-global ``ClfactBuffers``; that is safe under
    ``fork`` and would silently corrupt results under threads.  Oversubscribing
    (20 threads on 14 physical cores) makes any sharing show up.
    """
    p, st = _states(3, 1, 1, 17)
    n = st.numsta
    selfen = opt.compute_selfen(3)
    mstinf = st.mstinf[:n].copy()
    args = (1, st.mstate, mstinf, n, 3, 1, 1, 17, selfen, p.cbreak)
    one = opt.build_matrices(*args, 1, backend="thread")
    many = opt.build_matrices(*args, nthreads, backend="thread")
    assert np.array_equal(one[1], many[1])
    assert np.array_equal(one[0], many[0])


def test_unknown_backend_is_an_error():
    p, st = _states(3, 1, 1, 11)
    with pytest.raises(ValueError, match="backend"):
        opt.build_matrices(0, st.mstate, st.mstinf[:st.numsta], st.numsta,
                           3, 1, 1, 11, opt.compute_selfen(3), p.cbreak, 1,
                           backend="magic")


@pytest.mark.parametrize("B", [0, 1])
def test_run_python_agrees_across_backends(B):
    """End to end: the spectrum is identical whichever backend built it.

    This is what licenses ``PythonProvider`` to leave the backend out of its
    cache tag.
    """
    from dlcq.read_python import run_python

    kw = dict(N=3, NF=1, B=B, K_code=15, rlamb=0.3325, ncpus=4)
    a = run_python(backend="process", **kw)
    b = run_python(backend="thread", **kw)
    assert a.numsta_post == b.numsta_post
    assert np.array_equal(a.norm, b.norm)
    np.testing.assert_array_equal(a.eigenvalues, b.eigenvalues)
