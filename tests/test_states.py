"""The compiled state generator must reproduce the interpreted one exactly.

``qcdf_states.py`` is a transcription of ``qcdf.brmsgn``/``prmx``/``prmm``, so
only a direct comparison establishes it is faithful.  These tests demand
**element-identical** ``mstate`` and ``mstinf``, not merely the same set of
states, because state *order* is load-bearing: ``WEEDR``/``WEEDR2`` drop
linearly dependent states by index, so a permuted basis would silently change
which states survive and, through them, the spectrum.

One deviation needs pinning in its own right.  The compiled path hoists the
momentum-parity test (``oddchk``) out of the candidate loop: parity depends
only on the ``kpx``/``kpm`` row selected by the quotient ``(ist[l,j]-1)//dim``,
never on the flavour remainder, and ``oddchk`` is a conjunction over the
baryon/antibaryon/meson groups, so filtering each group's rows once is
equivalent to testing every candidate.  That is what collapses a 46:1 candidate
funnel -- 209,152 built for 4,529 kept at 2K=29 -- and it is the one place a
subtle error would not show up as a crash.
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


# (N, NF, B, 2K, iflv).  Mesons and baryons, SU(2)/SU(3)/SU(4), one and two
# flavours, B=2, and both K parities.
CASES = [
    (3, 1, 1, 11, None),
    (3, 1, 1, 21, None),
    (3, 1, 1, 25, None),
    (3, 1, 0, 12, None),
    (3, 1, 0, 24, None),
    (2, 1, 1, 12, None),
    (2, 1, 0, 16, None),
    (4, 1, 1, 16, None),
    (3, 1, 2, 22, None),
    (3, 2, 0, 12, [0, 0]),
    (3, 2, 1, 11, [3, 0]),
]


def _fresh(N, NF, B, K, iflv=None):
    p = base.Params()
    p.N, p.NF, p.B, p.K, p.rlamb = N, NF, B, K, 0.3325
    p.cutoff, p.LPN = -1.0, 0
    flv = list(iflv) if iflv is not None else [0] * NF
    if iflv is None:
        flv[0] = N * B
    for i, v in enumerate(flv[:NF]):
        p.iflv[i] = v
    return p, base.StateData(), base.PermTables(), base.FlavorTables()


def _both(N, NF, B, K, iflv=None):
    ref = _fresh(N, NF, B, K, iflv)
    base.qcdsta(*ref)
    got = _fresh(N, NF, B, K, iflv)
    opt.qcdsta_fast(*got)
    return ref, got


@pytest.mark.parametrize("N,NF,B,K,iflv", CASES)
def test_states_are_element_identical(N, NF, B, K, iflv):
    """Same states, in the same order, to the last element."""
    (_, s_ref, _, _), (_, s_new, _, _) = _both(N, NF, B, K, iflv)
    n = s_ref.numsta
    assert s_new.numsta == n
    if n == 0:
        pytest.skip(f"N={N} NF={NF} B={B} 2K={K} generates no states")
    assert np.array_equal(s_ref.mstinf[:n], s_new.mstinf[:n]), "mstinf differs"
    assert np.array_equal(s_ref.mstate[:4 * n], s_new.mstate[:4 * n]), (
        "mstate differs -- state content or ORDER changed")


@pytest.mark.parametrize("N,NF,B,K,iflv", CASES)
def test_permutation_tables_are_identical(N, NF, B, K, iflv):
    """``prmx``/``prmm`` outputs match, including the locator arrays.

    ``kpxloc``/``kpmloc`` use 0 as a "no partitions" sentinel that ``brmsgn``
    branches on, so they must match everywhere, not just where populated.
    """
    (_, _, p_ref, _), (_, _, p_new, _) = _both(N, NF, B, K, iflv)
    assert p_new.numprm == p_ref.numprm
    assert p_new.nprmm == p_ref.nprmm
    assert np.array_equal(p_ref.kpx[:p_ref.numprm], p_new.kpx[:p_new.numprm])
    assert np.array_equal(p_ref.kpm[:p_ref.nprmm], p_new.kpm[:p_new.nprmm])
    assert np.array_equal(p_ref.kpxloc, p_new.kpxloc)
    assert np.array_equal(p_ref.kpmloc, p_new.kpmloc)


@pytest.mark.parametrize("N,NF,B,K,iflv", CASES)
def test_every_generated_state_has_odd_momenta(N, NF, B, K, iflv):
    """The property the hoisted ``oddchk`` is responsible for.

    Antiperiodic boundary conditions put every parton on an odd momentum.  If
    the row-level parity filter were wrong in the permissive direction, states
    with an even momentum would survive here even though the identity test
    above might still pass on a case where the reference agreed.
    """
    _, (_, s_new, _, _) = _both(N, NF, B, K, iflv)
    n = s_new.numsta
    if n == 0:
        pytest.skip("no states")
    for s in range(n):
        loc = int(s_new.mstinf[s, 0]) - 1
        L = int(s_new.mstinf[s, 1])
        moms = s_new.mstate[loc + 2, :L]
        assert np.all(moms % 2 == 1), f"state {s} has an even momentum: {moms}"


def test_stachk_matches_the_dict_version():
    """The occurrence-count rewrite of ``stachk`` is exactly equivalent.

    The original keys two dicts by ``val+1`` and ``-val+1``; both keyings are
    injective on their own domain and the dicts are separate, so the function
    reduces to "does any value occur more than ``indx`` times".
    """
    import qcdf_states as ks

    rng = np.random.default_rng(0)
    for _ in range(400):
        n = int(rng.integers(1, 12))
        vals = rng.integers(-4, 6, size=n).astype(np.int64)
        for indx in (1, 2, 3):
            assert (ks.stachk_nb(n, indx, vals)
                    is base.stachk(n, indx, vals.tolist()))


def test_growable_tables_are_not_capped_by_the_module_constants():
    """The compiled path sizes ``kpx``/``mstate`` to the run, not to a constant.

    ``qcdf.py``'s ``LKPXMX``/``NMSTMX`` bound the interpreted generator only.
    2K=37 needs 85,093 momentum partitions and 32,816 states, both past the
    defaults, and must still run.
    """
    p, s, e, f = _fresh(3, 1, 1, 37)
    opt.qcdsta_fast(p, s, e, f)
    assert s.numsta == 32816
    assert e.numprm == 85093
    assert e.numprm > base.LKPXMX or s.numsta > base.NMSTMX, (
        "raise the case until it exceeds a default, or this proves nothing")


@pytest.mark.slow
@pytest.mark.parametrize("B,K", [(1, 33), (1, 35), (0, 30), (0, 34)])
def test_identical_at_large_K(B, K):
    """Same check at sizes where generation used to take tens of seconds."""
    (_, s_ref, p_ref, _), (_, s_new, p_new, _) = _both(3, 1, B, K)
    n = s_ref.numsta
    assert s_new.numsta == n
    assert p_new.numprm == p_ref.numprm
    assert np.array_equal(s_ref.mstinf[:n], s_new.mstinf[:n])
    assert np.array_equal(s_ref.mstate[:4 * n], s_new.mstate[:4 * n])
