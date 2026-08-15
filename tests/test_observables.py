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
                              richardson_extrapolate, valence_parton_count,
                              physical_indices)
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


# ──────────────────────────────────────────────────────────────────────────
# Reach beyond the paper
# ──────────────────────────────────────────────────────────────────────────

@pytest.mark.slow
@pytest.mark.parametrize("B,K_code,min_pre", [(1, 31, 7569), (0, 30, 4141)])
def test_sum_rules_hold_past_the_papers_largest_K(B, K_code, min_pre):
    """2K = 30/31 exceeds every basis the paper used, and still closes exactly.

    The paper extrapolates "for 2K in the range of roughly 16-24", so 2K = 31
    is a basis it never reached.  Getting there needed no new physics, only two
    caps raised in ``qcdf.py``: ``LKPXMX`` (25001 momentum permutations, 32072
    needed) and ``NMSTMX`` (6902 states, 7569 needed).  This pins that they are
    in fact large enough, by checking the two sum rules that must hold exactly
    at any K -- if the state generator had silently truncated, they would not.
    """
    from dlcq.read_python import run_python
    from dlcq.units import mg_to_lambda

    r = run_python(N=3, NF=1, B=B, K_code=K_code, rlamb=mg_to_lambda(1.6),
                   ncpus=8)
    assert r.numsta_pre >= min_pre, (
        f"only {r.numsta_pre} states generated; a cap truncated the basis")
    i = int(physical_indices(r)[0])
    assert momentum_sum_rule(r, i) == pytest.approx(1.0, abs=1e-9)
    assert number_sum_rule(r, i) == pytest.approx(3 * B, abs=1e-9)


@pytest.mark.slow
def test_baryon_mass_still_converges_past_2K_24():
    """M^2(K) must keep rising monotonically toward its continuum limit.

    A cap overflow that quietly dropped high-Fock states would show up as the
    sequence turning over, since those states lower the ground state.
    """
    from dlcq.read_python import run_python
    from dlcq.units import mg_to_lambda

    lam = mg_to_lambda(1.6)
    masses = []
    for K in (27, 29, 31):
        r = run_python(N=3, NF=1, B=1, K_code=K, rlamb=lam, ncpus=8)
        masses.append(r.eigenvalues[int(physical_indices(r)[0])])
    assert masses[0] < masses[1] < masses[2], masses
    # Converging, not drifting: successive gaps must shrink.
    assert masses[2] - masses[1] < masses[1] - masses[0]


def test_richardson_stability_is_tiny_for_an_exact_series():
    """On data generated from the fit form itself, every sub-window agrees."""
    from dlcq.observables import richardson_stability
    from dlcq.units import endpoint_exponent

    mg, N = 1.6, 3
    a = endpoint_exponent(mg, N)
    Ks = np.array([25, 27, 29, 31, 33, 35])
    Kp = Ks / 2.0
    exact = 5.0 + 0.7 / Kp + 0.3 / Kp ** (1 + a)
    M0, spread, nwin = richardson_stability(Ks, exact, mg, N, n_terms=2)
    # Contiguous sub-windows [i, j) with at least min_points=4 points, the full
    # window excluded: sizes 4 and 5 give 3 + 2 = 5.  Pinned exactly rather than
    # bounded, because the count *is* the definition -- this used to enumerate
    # every subset (21 here, 2^n in general, which never returns at the 24
    # points of a 2K = 25-71 series).
    assert nwin == 5
    assert abs(M0 - 5.0) < 1e-8
    assert spread < 1e-8


def test_richardson_stability_needs_enough_points():
    from dlcq.observables import richardson_stability

    M0, spread, nwin = richardson_stability([21, 23, 25, 27],
                                            [1.0, 1.1, 1.15, 1.18], 1.6, 3)
    assert nwin < 2
    assert np.isnan(spread)


@pytest.mark.slow
def test_last_term_rule_overstates_the_error_at_weak_coupling():
    """The paper's error definition degrades as the Eq. (27) basis degenerates.

    Eq. (26)'s exponent ``a`` goes to zero in the chiral limit, so
    ``1/K^(1+a)`` collapses onto ``1/K``.  The fit is unaffected -- it still
    reproduces M(K) to five digits -- but the coefficients grow large and
    cancel, inflating the last term.  Refitting on sub-windows measures the
    actual stability and shows the gap widening as m/g falls.
    """
    from dlcq.figures import sweep_lpn
    from dlcq.observables import richardson_extrapolate, richardson_stability
    from dlcq.providers import PythonProvider
    from dlcq.units import code_to_M_over_g, mg_to_lambda

    prov = PythonProvider(ncpus=8)
    ratios = {}
    for mg in (1.6, 0.1):
        lam = float(mg_to_lambda(mg))
        Ks, ms = [], []
        for K in (25, 27, 29, 31, 33, 35):
            r = prov.get(3, 1, 1, K, lam, -1.0, sweep_lpn(3, 1))
            ph = physical_indices(r)
            if ph.size:
                Ks.append(K)
                ms.append(code_to_M_over_g(float(r.eigenvalues[ph[0]]), lam))
        _, last = richardson_extrapolate(Ks, ms, mg, 3, n_terms=2)
        _, spread, _ = richardson_stability(Ks, ms, mg, 3, n_terms=2)
        ratios[mg] = last / spread

    # Both overstate, and the weak-coupling one overstates far more.
    assert ratios[1.6] > 2
    assert ratios[0.1] > 50
    assert ratios[0.1] > 10 * ratios[1.6]
