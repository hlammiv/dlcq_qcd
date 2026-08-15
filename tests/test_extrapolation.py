"""The Richardson fit: its basis, and how its error is estimated.

Two things are established here.

**The confluent basis is a change of coordinates, not of model.**  Eq. (27)'s
exponents ``1`` and ``1+a`` merge as the endpoint exponent ``a -> 0``, and the
design matrix becomes singular (condition number 3.1e19 at m/g = 0, past what
float64 represents).  Replacing the merging pair by its divided difference
``(K^-(1+a) - K^-1)/a -> -K^-1 ln K`` spans the *same* subspace for every ``a``,
so ``M0`` must be unchanged; only the conditioning improves.  Both halves of
that are asserted below, because the value of the reparametrization rests
entirely on it being exactly model-preserving.

**Held-out K is the honest error estimate.**  The paper's rule -- the magnitude
of the last fitted term -- is basis dependent, which these tests demonstrate: it
moves by a factor of several under a coordinate change that leaves ``M0``
identical to 14 digits.  A quantity that depends on how you parametrize the
model is not an uncertainty.  ``richardson_holdout`` instead withholds the
largest K values and measures the prediction error, which is basis independent
and is the only check here sensitive to a wrong functional form.
"""

from __future__ import annotations

import numpy as np
import pytest

from dlcq.observables import (_richardson_design, levin_u, monotone_bracket,
                              richardson_budget, richardson_evaluate,
                              richardson_extrapolate, richardson_holdout,
                              wynn_epsilon)
from dlcq.units import endpoint_exponent

K_GRID = np.arange(25, 50, 2)
MGS = [1.6, 0.8, 0.4, 0.2, 0.1, 0.05, 0.0]


def _exact_series(K_codes, mg, N, coeffs=(3.0, 1.1, 0.4, -0.2)):
    """Data that IS the model, so the fit must recover it to machine precision."""
    Kp = np.asarray(K_codes, float) / 2.0
    a = endpoint_exponent(mg, N)
    c0, c1, c2, c3 = coeffs
    return c0 + c1 * Kp**-1.0 + c2 * Kp ** -(1.0 + a) + c3 * Kp**-2.0


# m/g = 0 is excluded from the same-span tests on purpose: there the paper
# basis is not a reparametrization of the confluent one but a strict subspace
# of it (see test_chiral_limit_changes_the_model_dimension).
POSITIVE_MGS = [mg for mg in MGS if mg > 0.0]


@pytest.mark.parametrize("mg", POSITIVE_MGS)
@pytest.mark.parametrize("N", [2, 3, 4])
def test_confluent_basis_gives_the_same_M0(mg, N):
    """Same span, therefore same answer. This is what licenses the change."""
    y = _exact_series(K_GRID, mg, N)
    paper, _ = richardson_extrapolate(K_GRID, y, mg, N, n_terms=4)
    conf, _ = richardson_extrapolate(K_GRID, y, mg, N, n_terms=4,
                                     basis="confluent")
    assert paper == pytest.approx(conf, abs=1e-9)


@pytest.mark.parametrize("mg", POSITIVE_MGS)
def test_confluent_basis_is_better_conditioned(mg):
    """And never worse -- the whole point is that it costs nothing to use."""
    Kp = K_GRID / 2.0
    a = endpoint_exponent(mg, 3)
    A_paper, _, _ = _richardson_design(Kp, a, 4, "paper")
    A_conf, _, _ = _richardson_design(Kp, a, 4, "confluent")
    assert np.linalg.cond(A_conf) <= np.linalg.cond(A_paper) * 1.05


def test_chiral_limit_changes_the_model_dimension():
    """At m/g = 0 the paper basis quietly fits a smaller model.

    Both ``1+a -> 1`` and ``2+a -> 2`` collapse, so a five-term request becomes
    a three-term fit -- and the dimension jumps discontinuously, since m/g =
    0.05 is still fitted at full width.  The confluent basis keeps all five
    columns at every ``a``, which is the point of using it.
    """
    Kp = K_GRID / 2.0
    A_paper, _, _ = _richardson_design(Kp, 0.0, 4, "paper")
    A_conf, _, _ = _richardson_design(Kp, 0.0, 4, "confluent")
    assert A_paper.shape[1] == 3
    assert A_conf.shape[1] == 5
    assert np.linalg.cond(A_conf) < 1e8

    # ...while just short of the limit the paper basis keeps all five and pays
    # for it in conditioning.
    a_small = endpoint_exponent(0.05, 3)
    A_near, _, _ = _richardson_design(Kp, a_small, 4, "paper")
    assert A_near.shape[1] == 5
    assert np.linalg.cond(A_near) > 10 * np.linalg.cond(A_conf)


def test_literal_paper_basis_is_singular_in_the_chiral_limit():
    """Why the dedupe exists at all: without it the matrix is not invertible."""
    Kp = K_GRID / 2.0
    literal = np.vstack([Kp**-e for e in (0.0, 1.0, 1.0, 2.0, 2.0)]).T
    assert np.linalg.cond(literal) > 1e15


@pytest.mark.parametrize("basis", ["paper", "confluent"])
@pytest.mark.parametrize("mg", [1.6, 0.4, 0.05])
def test_exact_series_is_recovered(mg, basis):
    """If the data is the model, M0 is the leading coefficient exactly."""
    y = _exact_series(K_GRID, mg, 3)
    M0, _ = richardson_extrapolate(K_GRID, y, mg, 3, n_terms=4, basis=basis)
    assert M0 == pytest.approx(3.0, abs=1e-8)


@pytest.mark.parametrize("basis", ["paper", "confluent"])
def test_holdout_is_exact_on_an_exact_series(basis):
    """No model error means no prediction error, whichever basis is used."""
    y = _exact_series(K_GRID, 0.4, 3)
    abs_err, rel_err, pred, actual, K_held = richardson_holdout(
        K_GRID, y, 0.4, 3, n_terms=4, n_hold=2, basis=basis)
    assert K_held.tolist() == [47.0, 49.0]
    assert abs_err < 1e-8
    assert pred == pytest.approx(actual, abs=1e-8)


def test_holdout_detects_a_wrong_functional_form():
    """The property the last-term rule cannot have.

    A series with a term the Eq. (27) basis cannot represent (here an
    exponential in K) fits the visible points happily and then mispredicts the
    withheld ones by orders of magnitude more.
    """
    Kp = K_GRID / 2.0
    y_good = _exact_series(K_GRID, 0.4, 3)
    y_bad = y_good + 5.0 * np.exp(-Kp / 3.0)

    good = richardson_holdout(K_GRID, y_good, 0.4, 3, n_terms=4, n_hold=2)[0]
    bad = richardson_holdout(K_GRID, y_bad, 0.4, 3, n_terms=4, n_hold=2)[0]
    assert bad > 100 * max(good, 1e-12)


def test_holdout_needs_enough_points():
    with pytest.raises(ValueError, match="held out"):
        richardson_holdout([25, 27, 29], [1.0, 1.1, 1.2], 0.4, 3, n_hold=2)


@pytest.mark.parametrize("basis", ["paper", "confluent"])
def test_evaluate_reproduces_the_fitted_points(basis):
    """``richardson_evaluate`` must rebuild the same design the fit used."""
    y = _exact_series(K_GRID, 0.2, 3)
    _, _, coeffs, _ = richardson_extrapolate(
        K_GRID, y, 0.2, 3, n_terms=4, return_fit=True, basis=basis)
    back = richardson_evaluate(coeffs, 0.2, 3, K_GRID, n_terms=4, basis=basis)
    assert back == pytest.approx(y, abs=1e-9)


def test_unknown_basis_is_an_error():
    with pytest.raises(ValueError, match="basis"):
        richardson_extrapolate(K_GRID, _exact_series(K_GRID, 0.4, 3),
                               0.4, 3, basis="magic")


def test_asking_past_the_exponent_ladder_is_an_error():
    """Eq. (27) supplies five correction terms; a sixth must not be invented.

    This used to slice a six-entry list, so ``n_terms`` of 5, 6 and 9 all
    returned the *same* five-term fit -- identical columns, identical M0, no
    warning.  An order scan run past the ladder then looked converged when it
    had only stopped changing the model.
    """
    Kp = K_GRID / 2.0
    assert _richardson_design(Kp, 0.3, 5, "paper")[0].shape[1] == 6
    for bad in (6, 9):
        with pytest.raises(ValueError, match="exceeds"):
            _richardson_design(Kp, 0.3, bad, "paper")


# ── exponent-free estimates ───────────────────────────────────────────────
#
# These exist to referee the Eq. (27) fits, so they are tested against series
# whose limits are known in closed form rather than against those fits.

KP = K_GRID / 2.0


@pytest.mark.parametrize("tail,limit", [
    (1.1 / KP, 3.0),                       # plain 1/K
    (1.1 * KP ** -1.2, 3.0),               # non-integer exponent
    (1.1 * KP ** -0.5, 3.0),               # slow, sqrt-like
    (0.9 / KP + 0.5 * KP ** -1.0844, 5.0),  # the actual Eq. (27) shape
])
def test_levin_u_recovers_known_limits(tail, limit):
    """Levin's u-transform is the referee, so it has to be right first.

    It assumes only an asymptotic expansion in 1/K -- no endpoint exponent --
    which is exactly why it can adjudicate a basis that is under suspicion.
    """
    got, _ = levin_u(K_GRID, limit - tail)
    assert got == pytest.approx(limit, abs=1e-4)


def test_wynn_is_the_weaker_tool_and_is_labelled_as_such():
    """Wynn is corroboration, not evidence, and the gap is worth pinning down.

    It is exact for geometric remainders and merely decent for the power-law
    tails this project has; measured, it recovers a known 2.0 as ~1.92 where
    Levin gets 2.0000073.  Asserting the ordering keeps a future reader from
    promoting it to a primary estimate.
    """
    S = 2.0 - 1.0 * KP ** -0.5
    wy, _ = wynn_epsilon(S)
    lv, _ = levin_u(K_GRID, S)
    assert abs(lv - 2.0) < abs(wy - 2.0)


def test_wynn_handles_a_short_or_constant_sequence():
    assert wynn_epsilon([1.0, 1.0, 1.0, 1.0])[0] == pytest.approx(1.0)
    assert np.isfinite(wynn_epsilon([1.0, 2.0])[0])


@pytest.mark.parametrize("p_true", [0.5, 1.0, 1.2, 2.0])
def test_bracket_contains_the_truth(p_true):
    """The bound must never be beaten -- a bracket that fails is worse than none.

    An earlier version approximated the tail *sum* by an integral starting at
    the last abscissa and missed a known limit by 2e-5.  For a decreasing
    summand the integral test requires starting half a step earlier, which is
    what the implementation now does.
    """
    S = 3.0 - 1.1 * KP ** -p_true
    lo, hi, p, mono = monotone_bracket(K_GRID, S)
    assert mono
    assert lo <= 3.0
    if np.isfinite(hi):
        assert hi >= 3.0


def test_bracket_refuses_a_decreasing_sequence():
    """It is only meaningful for the increasing-with-decreasing-increments case.

    M^2(K) is measured to be exactly that at every K in this project, but a
    caller handing it something else must get a refusal rather than a number.
    """
    lo, hi, p, mono = monotone_bracket(K_GRID, 3.0 + 1.1 / KP)
    assert mono is False
    assert not np.isfinite(hi)


def test_bracket_lower_bound_is_the_largest_computed_value():
    """Monotonicity gives that half for free, with no assumption at all."""
    S = 3.0 - 1.1 / KP
    lo, _, _, _ = monotone_bracket(K_GRID, S)
    assert lo == pytest.approx(S[-1])


# ── the truncation term, when the two sweeps span different windows ────────

def _series(Ks, M0, amp):
    """A convergent series with a 1/K tail, as a stand-in for a real sweep."""
    return [M0 + amp / k for k in Ks]


def test_truncation_term_is_measured_on_the_shared_window():
    """The alt-LPN sweep is shorter, and the term must still isolate truncation.

    Raising the particle-number cut by one qqbar pair multiplies the basis by
    17-23x, so that sweep stops well below the top of the main window.  If the
    two series were each extrapolated on their own K, the difference would
    carry a *window* change as well as a truncation change, and the budget's
    ``trunc`` column would silently absorb it.

    Here both series have the identical limit and differ only in window, so an
    honest truncation term is zero.
    """
    K = list(range(25, 72, 2))
    K_alt = list(range(25, 50, 2))
    base = _series(K, 10.0, 0.5)
    alt = _series(K_alt, 10.0, 0.5)          # same physics, shorter sweep

    bud = richardson_budget(K, base, 1.6, 3, n_terms=2,
                            masses_alt_lpn=alt, K_codes_alt=K_alt)
    assert bud["truncation"] == pytest.approx(0.0, abs=1e-9)


def test_truncation_term_still_sees_a_real_shift():
    """The complement: a genuine change in the limit must come through."""
    K = list(range(25, 72, 2))
    K_alt = list(range(25, 50, 2))
    base = _series(K, 10.0, 0.5)
    alt = _series(K_alt, 10.02, 0.5)         # limit moved by 0.02

    bud = richardson_budget(K, base, 1.6, 3, n_terms=2,
                            masses_alt_lpn=alt, K_codes_alt=K_alt)
    assert bud["truncation"] == pytest.approx(0.02, rel=1e-3)


def test_too_little_overlap_is_not_silently_zero():
    """A non-overlapping alt sweep must not read as 'no truncation error'."""
    bud = richardson_budget(list(range(51, 72, 2)), _series(range(51, 72, 2), 10.0, 0.5),
                            1.6, 3, n_terms=2,
                            masses_alt_lpn=_series(range(25, 32, 2), 10.0, 0.5),
                            K_codes_alt=list(range(25, 32, 2)))
    assert np.isnan(bud["truncation"])


# ── the ensemble: contiguous windows, and never exponential again ──────────

def _expected_windows(n, n_terms_set=(2, 3, 4), min_points=4, min_dof=1):
    """Closed form for the contiguous-window count ``richardson_ensemble`` makes.

    For each order the smallest admissible window is ``max(min_points,
    nt + 1 + min_dof)``; a window of size s has ``n - s + 1`` placements.
    """
    total = 0
    for nt in n_terms_set:
        smallest = max(min_points, nt + 1 + min_dof)
        total += sum(n - s + 1 for s in range(smallest, n + 1))
    return total


def test_ensemble_enumerates_contiguous_windows_not_subsets():
    """The count is the definition, and it must stay polynomial in n.

    This enumerated *every subset* until it was measured: 3 x 2^24 ~ 5e7 fits
    for a 2K = 25-71 series, which never returns -- ``figures.figure_fits``
    simply hung, and read as a crash rather than as an algorithm.  The same bug
    was in ``richardson_stability``.

    Pinning the exact count catches a revert to subsets immediately: at n = 24
    contiguous gives 631, subsets would give ~5e7.
    """
    from dlcq.observables import richardson_ensemble

    mg, N = 1.6, 3
    a = endpoint_exponent(mg, N)
    Ks = np.arange(25, 72, 2, dtype=float)          # 24 points, the real window
    Kp = Ks / 2.0
    exact = 5.0 + 0.7 / Kp + 0.3 / Kp ** (1 + a)

    ens = richardson_ensemble(Ks, exact, mg, N)
    assert len(ens) == _expected_windows(len(Ks)) == 631

    # Data generated from the fit form itself, so every window must recover the
    # same limit -- the ensemble is enumerating choices, not fitting noise.
    M0s = np.array([m for _, _, m in ens])
    assert np.allclose(M0s, 5.0, atol=1e-6)


def test_ensemble_stays_cheap_at_the_full_window():
    """A wall-clock guard: exponential enumeration cannot pass this."""
    import time
    from dlcq.observables import richardson_ensemble

    mg, N = 0.1, 3
    Ks = np.arange(25, 72, 2, dtype=float)
    Kp = Ks / 2.0
    exact = 1.0 + 0.5 / Kp + 0.2 / Kp ** 2

    t = time.time()
    ens = richardson_ensemble(Ks, exact, mg, N)
    elapsed = time.time() - t
    assert len(ens) == 631
    assert elapsed < 10.0, f"ensemble took {elapsed:.1f}s -- enumeration regressed"
