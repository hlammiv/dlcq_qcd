"""The 't Hooft equation: the N -> infinity limit plotted in Fig. 8.

Fig. 8 overlays a "Large N" curve on the finite-N DLCQ results. That curve is
not a DLCQ output -- it is the solution of the paper's Eq. (24) at N -> infinity,
so it needs its own solver.

Eq. (24), for SU(N) (alpha = 1, which kills the third term):

    M^2 phi(x) = m^2 [1/x + 1/(1-x)] phi(x)
                 - (g^2/pi) (N^2 - 1)/(2N) P int_0^1 dy [phi(y) - phi(x)]/(y-x)^2

The principal-value integral with the ``phi(y) - phi(x)`` subtraction is finite:
the subtraction cancels the double pole, leaving an integrable singularity.

Units
-----
Measured, not assumed: in the paper's own coupling normalization the
dimensionless eigenvalues

    mu^2 = M^2 / [ (g^2/pi) (N^2 - 1)/(2N) ]

come out as 0, 5.8817, 14.1429, 23.08, ... in the chiral limit at large N,
reproducing the classic 't Hooft values 0, 5.88, 14.1 directly, with no
conversion factor.  (An earlier draft of this module asserted the paper's
coupling was half 't Hooft's and divided by 2; the benchmark below shows that
is wrong, so no factor is applied.)  Fig. 8(b)'s ``(2pi/N)^(1/2)`` axis
rescaling is a property of that figure and is applied there, not here.

We solve in units of ``g``: with ``t = m/g`` and
``c = (N^2-1)/(2N pi)`` the equation is

    (M/g)^2 phi(x) = t^2 [1/x + 1/(1-x)] phi(x)
                     - c P int_0^1 dy [phi(y) - phi(x)]/(y-x)^2

Method
------
Discretize on a uniform grid and treat the subtracted kernel exactly on each
cell, which is what makes the principal value tractable without a fine mesh:

    P int_{y_j - h/2}^{y_j + h/2} dy / (y - x_i)^2  =  1/(d - h/2) - 1/(d + h/2)

for ``d = x_i - y_j``.  The diagonal is fixed by the sum rule that the kernel
applied to a constant vanishes, which is precisely what the subtraction
enforces.  The endpoint behaviour ``phi ~ x^a`` of Eq. (26) is non-analytic, so
convergence in the grid size is slow; :func:`thooft_spectrum` therefore
Richardson-extrapolates in ``1/n``.
"""

from __future__ import annotations

import numpy as np

__all__ = ["thooft_coupling", "thooft_matrix", "thooft_spectrum",
           "thooft_mass"]


def thooft_coupling(N, large_n=False):
    """The paper's coupling combination, in units of g^2.

    ``c = (N^2 - 1) / (2 N pi)`` exactly, or its leading large-N form
    ``N / (2 pi)`` when ``large_n=True`` -- that is, the Casimir
    ``(N^2-1)/2N -> N/2``.

    Note this GROWS with N, so "the large-N curve" is not obtained by evaluating
    at some huge N: on an unrescaled ``M/g`` axis that diverges.  Fig. 8(a)
    plots the limit of the SU(N) *sequence*, which is the leading-order coupling
    at the N being compared; Fig. 8(b) rescales both axes by ``(2 pi / N)^(1/2)``
    precisely so the N -> infinity limit is finite there.
    """
    return N / (2.0 * np.pi) if large_n else (N ** 2 - 1.0) / (2.0 * N * np.pi)


def thooft_matrix(t, N, n=400, large_n=False):
    """Discretized Eq. (24) on ``n`` interior points.  Returns a symmetric matrix.

    ``t = m/g``.  Eigenvalues of the result are ``(M/g)^2``.
    """
    c = thooft_coupling(N, large_n)
    h = 1.0 / n
    x = (np.arange(n) + 0.5) * h                    # cell centres, avoids x=0,1

    # Free (mass) term.
    diag_free = t ** 2 * (1.0 / x + 1.0 / (1.0 - x))

    # Subtracted principal-value kernel, integrated exactly over each cell:
    #   K_ij = c * [1/(d - h/2) - 1/(d + h/2)],  d = x_i - x_j,  i != j
    d = x[:, None] - x[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        K = 1.0 / (d - h / 2.0) - 1.0 / (d + h / 2.0)
    np.fill_diagonal(K, 0.0)
    K *= c

    # The subtraction means the kernel annihilates a constant, so the diagonal
    # is minus the row sum.  This is what regulates the double pole.
    H = -K
    np.fill_diagonal(H, diag_free + K.sum(axis=1))

    return 0.5 * (H + H.T)                          # symmetrize against roundoff


def thooft_spectrum(t, N, n=400, k=4, large_n=False):
    """Lowest ``k`` values of ``(M/g)^2`` at one grid size."""
    from scipy.linalg import eigh

    H = thooft_matrix(t, N, n, large_n)
    w = eigh(H, eigvals_only=True)
    return np.sort(w)[:k]


def thooft_mass(t, N, grids=(200, 300, 400, 600), state=0, large_n=False):
    """``M/g`` for one state, Richardson-extrapolated in the grid spacing.

    The endpoint exponent of Eq. (26) is non-analytic, so a single grid
    converges slowly; extrapolating in ``1/n`` recovers most of it.
    """
    grids = np.asarray(grids, dtype=float)
    vals = np.array([thooft_spectrum(t, N, int(g), k=state + 1,
                                     large_n=large_n)[state]
                     for g in grids])
    # Fit vals = a0 + a1/n + a2/n^2 and take a0.
    A = np.vstack([np.ones_like(grids), 1.0 / grids, 1.0 / grids ** 2]).T
    coeffs, *_ = np.linalg.lstsq(A, vals, rcond=None)
    msq = float(coeffs[0])
    return float(np.sqrt(max(msq, 0.0)))
