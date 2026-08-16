"""Exact checks on the equal-time Schwinger construction.

These are all *analytic* targets, not regression values, which is the point:
this project's dominant failure mode is convention errors, and a regression test
locks in whatever convention was in force when it was written.  Every assertion
here is something that can be worked out on paper.

The strong-coupling limit is the workhorse.  At ``x = 0`` the Hamiltonian is
diagonal, so the whole spectrum is exact by hand, and it independently pins the
electric weights, the electric constant and the mass term.
"""
from __future__ import annotations

import numpy as np
import pytest

from equaltime.schwinger_ed import (SQRT_PI_INV, charge_sector,
                                    electric_weights, gap_to_M_over_g,
                                    hamiltonian, mu_from, spectrum)


@pytest.mark.parametrize("n", [6, 8, 10, 12])
def test_strong_coupling_vacuum_is_exactly_zero(n):
    """``x=0``: the staggered vacuum has ``L_n = 0`` on every link, so ``W0 = 0``.

    This is what catches a dropped constant.  ``sum_n L_n^2`` expands as
    ``Q_n^2 + 2 Q_n c_n + c_n^2``; omitting the last term leaves every *gap*
    correct while shifting the ground state by exactly ``-N/8``.  Masses would
    never notice.  The condensate, being ``-d(epsilon)/dm`` on the ground-state
    energy itself, would be wrong everywhere.
    """
    assert spectrum(n, 0.0, 0.0, k=1)[0] == pytest.approx(0.0, abs=1e-10)


@pytest.mark.parametrize("n", [6, 8, 10, 12])
def test_strong_coupling_gap_is_exactly_one(n):
    """``x=0, m=0``: the cheapest excitation puts ``L = +-1`` on one link."""
    ev = spectrum(n, 0.0, 0.0, k=2)
    assert float(ev[1] - ev[0]) == pytest.approx(1.0, abs=1e-10)


@pytest.mark.parametrize("mu", [0.3, 1.0, 2.5])
def test_strong_coupling_gap_with_mass(mu):
    """``x=0``: flipping a pair costs 1 (electric) + ``2 mu`` (staggered mass)."""
    import scipy.sparse as sp
    from equaltime.schwinger_ed import _ops
    n = 8
    H = hamiltonian(n, 0.0, 0.0)
    SZ, _ = _ops(n)
    for k in range(n):
        H = H + mu * ((-1) ** k) * SZ[k]
    sel = charge_sector(n)
    ev = np.linalg.eigvalsh(H[sel][:, sel].toarray())[:2]
    assert float(ev[1] - ev[0]) == pytest.approx(1.0 + 2.0 * mu, abs=1e-10)


def test_electric_weights_are_the_link_count():
    """``w[k,l] = N-1-max(k,l)``: the links lying to the right of both sites."""
    w, C, E0 = electric_weights(8)
    assert w[0, 0] == 7 and w[3, 5] == 2 and w[7, 7] == 0
    assert w[2, 6] == w[6, 2]
    # E0 = sum_n c_n^2 with c_n = 1/2 on even links: N/2 links x 1/4
    assert E0 == pytest.approx(8 / 8.0)


def test_charge_sector_size_is_the_central_binomial():
    """Total ``S^z = 0`` on N sites is ``C(N, N/2)`` states."""
    from math import comb
    for n in (6, 8, 10):
        assert len(charge_sector(n)) == comb(n, n // 2)


def test_unit_conversions_round_trip():
    """``mu`` and the gap conversion are the two places conventions live."""
    for mg in (0.01, 0.1, 1.0):
        for x in (4.0, 25.0, 100.0):
            assert mu_from(mg, x) == pytest.approx(2.0 * mg * np.sqrt(x))
            dW = 2.0 * np.sqrt(x) * 0.5641896
            assert gap_to_M_over_g(dW, x) == pytest.approx(SQRT_PI_INV, rel=1e-6)


def test_hamiltonian_is_hermitian():
    for n in (6, 8):
        H = hamiltonian(n, 3.0, 0.2)
        assert abs((H - H.T.conj())).max() < 1e-12


@pytest.mark.parametrize("L,x,mg", [(8, 2.0, 0.0), (8, 4.0, 0.3),
                                    (10, 2.0, 0.0), (10, 4.0, 0.3)])
def test_mpo_reproduces_ed(L, x, mg):
    """The MPO must equal the ED build exactly.

    This is the hopping term's only check -- the analytic ``x = 0`` tests above
    pin the electric weights, the electric constant and the mass term, but say
    nothing about the hop.  Machine precision or it is wrong.
    """
    from equaltime.schwinger_mps import ground_state
    ed = float(spectrum(L, x, mg, k=1)[0])
    dmrg, _, _ = ground_state(L, x, mg, chi=256)
    assert dmrg == pytest.approx(ed, abs=1e-9)


@pytest.mark.parametrize("L,x", [(10, 4.0), (12, 4.0), (12, 9.0)])
def test_excited_state_dmrg_matches_ed(L, x):
    """The gap from orthogonality-constrained DMRG must equal ED's.

    Guards a silent failure: ``orthogonal_to`` is a keyword-only constructor
    argument, and passing it in the options dict is ignored -- the second solve
    re-converges to the ground state and the gap is exactly 0.  A blanket
    ``warnings.filterwarnings("ignore")`` also hides the unused-option warning
    that would otherwise flag it.
    """
    from equaltime.schwinger_mps import mass_gap
    from equaltime.schwinger_ed import gap_to_M_over_g
    ed = spectrum(L, x, 0.0, k=2)
    want = gap_to_M_over_g(float(ed[1] - ed[0]), x)
    assert mass_gap(L, x, 0.0, chi=200) == pytest.approx(want, abs=1e-8)
