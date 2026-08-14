"""The norm matrix is a Gram matrix, so it must be positive semidefinite.

This is the invariant that catches the colour-contraction overflow, and it is
the only one that can.  ``MXTRM = 12552`` bounds the colour term tables in
*both* ``qcdf_kernels.clfact_nb`` and the interpreted ``qcdf_opt.clfact``.
Because the two paths share the bound they agreed bit-for-bit on the *wrong*
answer, so every ``array_equal`` assertion in ``test_kernels.py`` passed while
the matrix was corrupt.  Backend equality cannot see a defect both backends
have.

What it looked like, measured at N=4, B=1, 2K=20 against
``tools/colour_norm.norm_bruteforce`` (which builds the true Gram matrix by
explicit colour summation and is therefore PSD by construction):

    25 wrong entries, every one between L=12 states
    the worst was a *diagonal*: 150048 against a true 331776
    the assembled matrix had an eigenvalue of -1.8e5

A state whose own norm is under half its true value is not a rounding problem.
The two guards now raise instead of returning ``0.0`` / breaking out of the sum,
because a wrong answer that is returned quietly is worse than no answer.

Scope, measured: this is an **N >= 4 baryon** defect, not a parton-count
threshold -- N=3 is clean at L=15 (2K=41) while N=4 breaks at L=12.  No
published result is affected: Table I and Figs. 7-8 run under ``sweep_lpn``
(L <= 6) and Figs. 4-6 are N=3.

The PSD check is done **per configuration block**.  The norm is exactly
block-diagonal in parton content -- off-block entries are 0.0 exactly, not
1e-16 -- so this is equivalent to checking the whole matrix, and it costs 46 s
at 2K=41 where the dense norm would be 51 GB.
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
from dlcq.read_python import config_block_labels         # noqa: E402


def _norm(N, B, K, LPN=0, backend="thread", ncpus=4):
    """The pre-weeding norm, plus the state list it was built from."""
    p = base.Params()
    p.N, p.NF, p.B, p.K, p.rlamb = N, 1, B, K, 0.3325
    p.cutoff, p.LPN = -1.0, LPN
    p.iflv[0] = N * B
    st = base.StateData()
    opt.qcdsta_fast(p, st, base.PermTables(), base.FlavorTables())
    n = st.numsta
    if n == 0:
        pytest.skip(f"N={N} B={B} 2K={K} generates no states")
    _, _, hnorm = opt.build_matrices(
        0, st.mstate, st.mstinf[:n].copy(), n, N, 1, B, K,
        opt.compute_selfen(N), p.cbreak, ncpus, backend=backend)
    return hnorm[:n, :n], st, n


def _worst_block_eigenvalue(hnorm, st, n):
    """Most negative eigenvalue over the config blocks, relative to its scale.

    Relative, because the norm entries span several orders of magnitude with K
    and an absolute threshold would either miss real corruption at small K or
    fire on rounding at large K.
    """
    labels = config_block_labels(st.mstate, st.mstinf[:n].copy(), n)
    worst = 0.0
    for b in np.unique(labels):
        idx = np.flatnonzero(labels == b)
        block = hnorm[np.ix_(idx, idx)]
        w = np.linalg.eigvalsh(block)
        scale = max(float(np.abs(np.diag(block)).max()), 1.0)
        worst = min(worst, float(w.min()) / scale)
    return worst


# (N, B, 2K) -- N=3 mesons and baryons, small enough for the default tier.
PSD_CASES = [(3, 1, 21), (3, 1, 25), (3, 0, 20), (3, 0, 24), (2, 1, 22)]


@pytest.mark.parametrize("N,B,K", PSD_CASES)
def test_norm_is_positive_semidefinite(N, B, K):
    """A Gram matrix has no negative eigenvalues; a truncated colour sum does."""
    hnorm, st, n = _norm(N, B, K)
    worst = _worst_block_eigenvalue(hnorm, st, n)
    assert worst > -1e-9, (
        f"norm is not positive semidefinite at N={N} B={B} 2K={K}: worst "
        f"relative block eigenvalue {worst:.3g}. The colour sum has been "
        f"truncated -- see this module's docstring.")


def test_the_configuration_that_was_silently_wrong_is_now_exact():
    """N=4, B=1, 2K=20 -- the case that exposed the defect.

    Its norm used to come back with a diagonal of 150048 against a true 331776
    and an eigenvalue of -1.8e5.  With MXTRM at 50000 it must now agree with
    explicit colour enumeration *exactly*: this is a sum of integer-weighted
    contractions, so anything other than bit equality means terms are still
    being dropped.
    """
    pytest.importorskip("scipy")
    sys.path.insert(0, str(ROOT / "tools"))
    from colour_norm import norm_bruteforce

    hnorm, st, n = _norm(4, 1, 20)
    G = norm_bruteforce(st, n, 4)
    assert np.abs(hnorm - G).max() == 0.0
    assert _worst_block_eigenvalue(hnorm, st, n) > -1e-9


def test_both_backends_overflow_alike():
    """The guards must stay in lockstep across the two paths.

    ``qcdf_kernels.py`` and ``qcdf_opt.py`` carry separate copies of MXTRM.  If
    they ever diverge, one backend raises where the other returns a number, and
    the bit-identity contract in ``test_kernels.py`` breaks in a confusing way.
    Uses a Hamiltonian build at N=4, 2K=20, which still exceeds the limit.
    """
    def outcome(backend):
        p = base.Params()
        p.N, p.NF, p.B, p.K, p.rlamb = 4, 1, 1, 20, 0.3325
        p.cutoff, p.LPN = -1.0, 0
        p.iflv[0] = 4
        st = base.StateData()
        opt.qcdsta_fast(p, st, base.PermTables(), base.FlavorTables())
        n = st.numsta
        try:
            opt.build_matrices(1, st.mstate, st.mstinf[:n].copy(), n, 4, 1, 1,
                               20, opt.compute_selfen(4), p.cbreak, 4,
                               backend=backend)
            return "ok"
        except OverflowError:
            return "overflow"

    assert outcome("thread") == outcome("process") == "overflow"


@pytest.mark.slow
@pytest.mark.parametrize("N,B,K", [(3, 1, 31), (3, 1, 35), (3, 0, 30)])
def test_norm_is_positive_semidefinite_at_large_K(N, B, K):
    """Where the parton count is highest and the term tables come closest."""
    hnorm, st, n = _norm(N, B, K)
    worst = _worst_block_eigenvalue(hnorm, st, n)
    assert worst > -1e-9, (
        f"norm is not positive semidefinite at N={N} B={B} 2K={K}: {worst:.3g}")


@pytest.mark.slow
def test_norm_matches_bruteforce_colour_sum():
    """The diagrammatic norm against explicit colour enumeration.

    ``tools/colour_norm`` expands every state over (type, k, flavour, colour)
    and takes inner products, so it is the true Gram matrix and is independent
    of the diagrammatic machinery entirely.
    """
    sys.path.insert(0, str(ROOT / "tools"))
    from colour_norm import norm_bruteforce

    hnorm, st, n = _norm(3, 1, 21)
    G = norm_bruteforce(st, n, 3)
    assert np.abs(hnorm - G).max() < 1e-9
