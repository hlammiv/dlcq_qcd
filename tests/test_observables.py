"""Tier 1 -- sum rules and analytic anchors.

These are the strongest tests available: they hold at every K and every
coupling, follow from the paper's own equations, and need no digitized data.
The sum rules in particular gate everything downstream -- if they fail, every
wavefunction figure is wrong.
"""

import numpy as np
import pytest

from dlcq.observables import (structure_function, momentum_sum_rule,
                              number_sum_rule, thooft_valence_limit,
                              richardson_extrapolate, valence_parton_count)
from dlcq.units import endpoint_exponent


pytestmark = pytest.mark.fortran


# ── sum rules ─────────────────────────────────────────────────────────────

def test_momentum_sum_rule_ground_state(fortran_K21):
    """int x [q + qbar] dx = 1 -- exact, since every state carries K_code."""
    assert momentum_sum_rule(fortran_K21, 0) == pytest.approx(1.0, abs=1e-10)


def test_number_sum_rule_ground_state(fortran_K21):
    """int [q - qbar] dx = N*B = 3."""
    assert number_sum_rule(fortran_K21, 0) == pytest.approx(3.0, abs=1e-10)


def test_sum_rules_hold_for_every_available_eigenvector(fortran_K21):
    """Not just the ground state -- all 75 printed eigenvectors."""
    worst_p = max(abs(momentum_sum_rule(fortran_K21, i) - 1.0)
                  for i in range(fortran_K21.n_eigenvectors))
    worst_n = max(abs(number_sum_rule(fortran_K21, i) - 3.0)
                  for i in range(fortran_K21.n_eigenvectors))
    assert worst_p < 1e-10, f"momentum sum rule worst deviation {worst_p:.2e}"
    assert worst_n < 1e-10, f"number sum rule worst deviation {worst_n:.2e}"


def test_fock_sectors_partition_the_number_sum_rule(fortran_K21):
    """Summing the sectors must recover the total exactly."""
    dx = 2.0 / fortran_K21.K_code
    total = 0.0
    for npart in np.unique(fortran_K21.state_len):
        _, q, qbar = structure_function(fortran_K21, 0, nparton=int(npart))
        total += float(np.sum(q - qbar) * dx)
    assert total == pytest.approx(3.0, abs=1e-10)


def test_higher_fock_is_strongly_suppressed(fortran_K21):
    """The paper: higher Fock states are suppressed by 2 to 4 orders of magnitude."""
    dx = 2.0 / fortran_K21.K_code
    frac = {}
    for npart in (3, 5, 7):
        _, q, qbar = structure_function(fortran_K21, 0, nparton=npart)
        frac[npart] = float(np.sum(q - qbar) * dx) / 3.0
    assert frac[3] > 0.99                      # valence dominates
    assert 1e-5 < frac[5] < 1e-2               # next sector, 2-4 orders down
    assert frac[5] > frac[7]                   # monotone suppression


def test_structure_function_lives_on_the_odd_momentum_grid(fortran_K21):
    x, q, qbar = structure_function(fortran_K21, 0)
    K = fortran_K21.K_code
    np.testing.assert_allclose(x, np.arange(1, K, 2) / K)
    assert np.all(q >= -1e-12)                 # a number density


def test_valence_parton_count():
    assert valence_parton_count(3, 0) == 2     # meson: q qbar
    assert valence_parton_count(3, 1) == 3     # baryon: qqq
    assert valence_parton_count(3, 2) == 6
    assert valence_parton_count(4, 1) == 4


# ── analytic limits ───────────────────────────────────────────────────────

def test_thooft_valence_limit_formula():
    """Eq. (22): q(x) = N(N-1)(1-x)^(N-2)."""
    x = np.linspace(0, 1, 11)
    np.testing.assert_allclose(thooft_valence_limit(x, 3), 6 * (1 - x))
    np.testing.assert_allclose(thooft_valence_limit(x, 2), np.full_like(x, 2.0))
    # Normalization: integrates to N over [0,1].
    xs = np.linspace(0, 1, 200001)
    for N in (2, 3, 4):
        assert np.trapezoid(thooft_valence_limit(xs, N), xs) == pytest.approx(N, rel=1e-4)


# ── Richardson extrapolation ──────────────────────────────────────────────

def test_richardson_recovers_a_known_continuum_limit():
    """Synthetic M(K) built from Eq. (27) must extrapolate back to M(0)."""
    N, mg, M0_true = 3, 1.6, 4.618
    a = endpoint_exponent(mg, N)
    K_codes = np.arange(16, 25, 2)
    Kp = K_codes / 2.0
    masses = M0_true + 0.7 / Kp + 0.3 / Kp ** (1 + a) - 0.2 / Kp ** 2
    M0, last = richardson_extrapolate(K_codes, masses, mg, N, n_terms=4)
    assert M0 == pytest.approx(M0_true, abs=1e-6)
    assert last >= 0.0


def test_richardson_is_exact_for_a_constant():
    K = np.arange(16, 25, 2)
    M0, last = richardson_extrapolate(K, np.full(K.size, 3.21), 0.8, 3)
    assert M0 == pytest.approx(3.21, abs=1e-9)
    assert last < 1e-9


def test_richardson_handles_the_chiral_limit_degeneracy():
    """At m/g = 0, a = 0 makes 1/K^(1+a) degenerate with 1/K; must not blow up."""
    K = np.arange(16, 25, 2)
    Kp = K / 2.0
    M0, _ = richardson_extrapolate(K, 1.0 + 0.5 / Kp, 0.0, 3)
    assert M0 == pytest.approx(1.0, abs=1e-6)
