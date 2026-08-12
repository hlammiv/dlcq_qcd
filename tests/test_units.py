"""Tier 0 -- unit conversions and the K-recovery machinery.

Fast, no solver required.
"""

import numpy as np
import pytest

from dlcq.units import (lambda_to_mg, mg_to_lambda, code_to_M_over_g,
                        M_over_g_to_code, allowed_momenta, allowed_x,
                        infer_K_from_x_grid, infer_K_from_normalization,
                        endpoint_exponent, meson_baryon_ratio_bosonization,
                        K_code_from_paper, K_paper_from_code, thooft_rescale)


@pytest.mark.parametrize("mg", [0.0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.2, 1.6, 5.0])
def test_coupling_roundtrip(mg):
    assert lambda_to_mg(mg_to_lambda(mg)) == pytest.approx(mg, abs=1e-12)


def test_lambda_for_paper_coupling():
    """m/g = 1.6 is the paper's weak-coupling point; inputs use 0.3325."""
    assert mg_to_lambda(1.6) == pytest.approx(0.3325494921702995, rel=1e-12)
    # The rounded value the runs actually used maps back to m/g ~ 1.6002.
    assert lambda_to_mg(0.3325) == pytest.approx(1.6, abs=1e-3)


def test_mass_conversion_roundtrip():
    lam = mg_to_lambda(1.6)
    for Msq in [0.5, 10.39038, 137.0]:
        assert M_over_g_to_code(code_to_M_over_g(Msq, lam), lam) == pytest.approx(Msq)


def test_negative_eigenvalue_clamps():
    """The Fortran's own output has small negative eigenvalues at 2K=25, 29."""
    assert code_to_M_over_g(-0.11, 0.5) == 0.0


def test_K_conventions():
    assert K_paper_from_code(24) == 12.0
    assert K_code_from_paper(10.5) == 21
    with pytest.raises(ValueError):
        K_code_from_paper(10.3)


def test_momenta_are_odd_and_bounded():
    """Antiperiodic BCs: odd momenta only; spectators take the rest."""
    assert allowed_momenta(24, B=0, N=3).tolist() == list(range(1, 24, 2))
    # A B=1, N=3 baryon leaves 2 spectators, so k <= K-2.
    assert allowed_momenta(21, B=1, N=3).tolist() == list(range(1, 20, 2))
    for K, B in [(24, 0), (21, 1), (22, 2)]:
        assert np.all(allowed_momenta(K, B=B, N=3) % 2 == 1)


@pytest.mark.parametrize("K,B", [(24, 0), (21, 1), (16, 0), (13, 1), (10, 0)])
def test_infer_K_from_exact_grid(K, B):
    """The x-grid pins K: x_min = 1/K and dx = 2/K."""
    K_hat, diag = infer_K_from_x_grid(allowed_x(K, B=B, N=3))
    assert K_hat == K
    assert diag["consistent"]
    assert diag["max_grid_residual"] == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("K,B", [(24, 0), (21, 1)])
def test_infer_K_survives_digitization_noise(K, B):
    """1% jitter is well beyond what a 612 ppi trace should produce."""
    rng = np.random.default_rng(12345)
    x = allowed_x(K, B=B, N=3)
    for _ in range(20):
        K_hat, _ = infer_K_from_x_grid(x * (1 + rng.normal(0, 0.01, x.size)))
        assert K_hat == K


def test_infer_K_from_normalization():
    """The vertical scale pins K independently: sum q(x_i) / K_paper = quark number."""
    K, B, N = 21, 1, 3
    x = allowed_x(K, B=B, N=N)
    # A normalized q(x) whose integral with dx = 2/K is exactly N*B = 3.
    shape = (1 - x) ** 2
    q = shape * (N * B) / (np.sum(shape) * (2.0 / K))
    K_hat, diag = infer_K_from_normalization(q, quark_number=N * B)
    assert K_hat == K


def test_endpoint_exponent_limits():
    """Eq. (26): a runs from 0 in the chiral limit to 1 in the free limit."""
    assert endpoint_exponent(0.0, 3) == 0.0
    assert endpoint_exponent(1e-6, 3) < 1e-2
    assert endpoint_exponent(1000.0, 3) > 0.99
    # Monotone in m/g.
    vals = [endpoint_exponent(mg, 3) for mg in [0.05, 0.1, 0.4, 0.8, 1.6]]
    assert all(b > a for a, b in zip(vals, vals[1:]))


def test_bosonization_ratio():
    """2 sin[pi/(2(2N-1))]; the paper quotes agreement to about 10%."""
    assert meson_baryon_ratio_bosonization(2) == pytest.approx(1.0, abs=1e-12)
    assert meson_baryon_ratio_bosonization(3) == pytest.approx(0.618034, abs=1e-6)
    assert meson_baryon_ratio_bosonization(4) == pytest.approx(0.445042, abs=1e-6)


def test_thooft_rescale():
    """Fig. 8(b) axes; this paper's unit is g^2 N/2pi, not 't Hooft's g^2 N/pi."""
    assert thooft_rescale(3) == pytest.approx(np.sqrt(2 * np.pi / 3))
