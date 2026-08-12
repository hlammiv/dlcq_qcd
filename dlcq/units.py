"""Unit and convention conversions for DLCQ QCD(1+1).

Single source of truth for the conventions in Hornbostel, Brodsky & Pauli,
Phys. Rev. D 41, 3814 (1990).  Every other module imports from here rather
than re-deriving; getting any of these wrong silently rescales entire figures.

The K convention is the one that bites hardest:

    K_paper = (L/2pi) P^+     is half-integral (momenta are half-odd-integers)
    K_code  = 2 * K_paper     is the integer the Fortran and Python solvers take

So the paper's "2K = 24" means ``K_code = 24``, and the paper's own K is 12.
Confirmed independently by the cache filename produced by the existing Python
code (``...K24_lam0.332549...``) sitting alongside Fig. 5's "m/g = 1.6, 2K = 24".
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "lambda_to_mg",
    "mg_to_lambda",
    "code_to_M_over_g",
    "M_over_g_to_code",
    "K_paper_from_code",
    "K_code_from_paper",
    "allowed_momenta",
    "allowed_x",
    "infer_K_from_x_grid",
    "infer_K_from_normalization",
    "thooft_rescale",
    "endpoint_exponent",
    "meson_baryon_ratio_bosonization",
]


# ──────────────────────────────────────────────────────────────────────────
# Coupling
# ──────────────────────────────────────────────────────────────────────────

def lambda_to_mg(lam):
    """Dimensionless coupling lambda -> quark-mass ratio m/g.

    Inverse of :func:`mg_to_lambda`.  lambda = 1 is infinite coupling (m/g = 0,
    the chiral limit); lambda -> 0 is the free limit (m/g -> infinity).

    m/g = sqrt((1 - lambda^2) / (pi lambda^2))
    """
    lam = np.asarray(lam, dtype=float)
    return np.sqrt((1.0 - lam**2) / (np.pi * lam**2))


def mg_to_lambda(mg):
    """Quark-mass ratio m/g -> dimensionless coupling lambda.

    lambda = 1 / sqrt(1 + pi (m/g)^2)

    This is the x-axis of Fig. 2, labelled ``1/(1 + pi m^2/g^2)^(1/2)``.
    m/g = 1.6 gives lambda = 0.3325489..., the ``0.3325`` in every input file.
    """
    mg = np.asarray(mg, dtype=float)
    return 1.0 / np.sqrt(1.0 + np.pi * mg**2)


# ──────────────────────────────────────────────────────────────────────────
# Mass eigenvalues
# ──────────────────────────────────────────────────────────────────────────

def code_to_M_over_g(Msq_code, lam):
    """Solver eigenvalue -> M/g.

    The solvers report ``M^2_code = K_code * w / 2`` (qcdf.f:703-707), which is
    M^2 in units of (m^2 + g^2/pi) -- the y-axis of Fig. 2.  Converting to the
    M/g of Fig. 7, Fig. 8 and Table I:

        M/g = sqrt(M^2_code / (pi lambda^2))

    Negative eigenvalues (a genuine finite-K artifact at strong coupling, and
    present in the Fortran's own output at K=25 and K=29) clamp to 0 rather
    than producing NaN.
    """
    Msq_code = np.asarray(Msq_code, dtype=float)
    lam = np.asarray(lam, dtype=float)
    return np.sqrt(np.maximum(Msq_code, 0.0) / (np.pi * lam**2))


def M_over_g_to_code(M_over_g, lam):
    """Inverse of :func:`code_to_M_over_g`."""
    M_over_g = np.asarray(M_over_g, dtype=float)
    lam = np.asarray(lam, dtype=float)
    return (M_over_g**2) * np.pi * lam**2


# ──────────────────────────────────────────────────────────────────────────
# Harmonic resolution
# ──────────────────────────────────────────────────────────────────────────

def K_paper_from_code(K_code):
    """K_code -> the paper's K (may be half-integral)."""
    return K_code / 2.0


def K_code_from_paper(K_paper):
    """The paper's K -> K_code (must come out integral)."""
    K_code = 2.0 * K_paper
    if abs(K_code - round(K_code)) > 1e-9:
        raise ValueError(f"K_paper={K_paper} does not give an integer K_code")
    return int(round(K_code))


def allowed_momenta(K_code, B=0, N=3):
    """Odd integer momenta ``k`` a single parton may carry, in code units.

    Antiperiodic boundary conditions force half-odd-integer momenta, i.e. odd
    integers once doubled.  The largest a single parton can take is limited by
    the other partons, each of which must carry at least 1:

      * meson (B=0): at least one other parton   -> k <= K_code - 1
      * baryon (B=1, N colors): at least N-1     -> k <= K_code - (N-1)

    Returns the ascending array of allowed k.
    """
    n_spectators = 1 if B == 0 else N * abs(B) - 1
    k_max = K_code - n_spectators
    if k_max < 1:
        return np.array([], dtype=int)
    return np.arange(1, k_max + 1, 2, dtype=int)


def allowed_x(K_code, B=0, N=3):
    """The discrete ``x = k / K_code`` grid a structure function lives on."""
    return allowed_momenta(K_code, B=B, N=N) / float(K_code)


def infer_K_from_x_grid(x_values, atol=0.02, K_max=60, min_inliers=0.6):
    """Recover K_code from the x-positions of plotted structure-function points.

    Because momenta are odd integers, ``x = k/K_code`` for k = 1, 3, 5, ...
    The grid therefore pins K two independent ways:

        leftmost point   x_min = 1 / K_code
        point spacing    dx    = 2 / K_code

    That matters because the paper never states K for Figs. 3 and 4.  Digitize
    the marker positions and K falls out.  Figs. 5 and 6 *do* state K (2K = 24
    and 2K = 21), so run this on them first to validate the digitization before
    trusting it where K is unknown.

    Method: fit the lattice rather than trust a single spacing.  For each
    candidate K, score how well the observed x sit on ``{k/K : k odd}`` and take
    the **smallest** K that fits within ``atol``.  Smallest matters -- a finer
    lattice always fits at least as well (K=41 contains K=21's points), so
    minimising residual alone would run away to large K.

    Digitized traces are messy: markers can be missed where curves cross, and
    spurious blobs appear at intersections.  Scoring each point against its
    nearest lattice site tolerates both, where a median-spacing estimate does
    not.

    Returns ``(K_code, diagnostics)``.
    """
    x = np.sort(np.asarray(x_values, dtype=float))
    if x.size == 0:
        raise ValueError("no x values supplied")

    # A point at x = 0 is a plotted endpoint, not a momentum site (k >= 1).
    x = x[x > 1e-9]
    if x.size == 0:
        raise ValueError("no positive x values supplied")

    def score_for(K):
        """(inlier fraction, median residual) of the observed x against K's lattice."""
        grid = np.arange(1, K, 2) / float(K)
        if grid.size == 0:
            return 0.0, np.inf
        dist = np.array([np.min(np.abs(grid - xi)) for xi in x])
        return float(np.mean(dist <= atol)), float(np.median(dist))

    scores = {K: score_for(K) for K in range(4, K_max + 1)}

    # Smallest K that explains most of the points.  Max-residual scoring would
    # be defeated by a single spurious blob; requiring only a majority of
    # inliers tolerates the intersection artifacts a real trace contains.
    fitting = [K for K in sorted(scores) if scores[K][0] >= min_inliers]
    K_code = (fitting[0] if fitting
              else max(scores, key=lambda K: (scores[K][0], -scores[K][1])))

    inlier_frac, median_res = scores[K_code]
    grid = np.arange(1, K_code, 2) / float(K_code)
    max_res = float(np.max([np.min(np.abs(grid - xi)) for xi in x]))

    # Diagnostics that do not drive the choice but help judge a bad trace.
    K_from_min = 1.0 / x[0]
    spacings = np.diff(x)
    big = spacings[spacings > 0.5 * np.median(spacings)] if spacings.size else spacings
    K_from_spacing = 2.0 / float(np.median(big)) if big.size else np.nan

    diagnostics = {
        "K_code": K_code,
        "inlier_fraction": inlier_frac,
        "median_grid_residual": median_res,
        "max_grid_residual": max_res,
        "consistent": bool(inlier_frac >= min_inliers),
        "K_from_x_min": K_from_min,
        "K_from_spacing": K_from_spacing,
        "n_points": int(x.size),
    }
    return K_code, diagnostics


def infer_K_from_normalization(q_values, quark_number):
    """Recover K_code from the *vertical* scale of a digitized structure function.

    Independent of :func:`infer_K_from_x_grid`, which uses only the horizontal
    axis.  Eq. (12) defines ``q(x) = K_paper <phi| b_k^dag b_k |phi>`` and the
    grid spacing is ``dx = 1 / K_paper``, so

        integral q dx = sum_i q(x_i) / K_paper = quark number

    giving ``K_paper = sum_i q(x_i) / quark_number``.  Pass ``quark_number = 1``
    for a meson's quark distribution, ``N*B`` for a baryon's.

    Agreement between this and the x-grid estimate is a strong check that a
    digitized curve was traced and calibrated correctly.
    """
    total = float(np.sum(np.asarray(q_values, dtype=float)))
    if quark_number == 0:
        raise ValueError("quark_number must be non-zero")
    K_paper = total / float(quark_number)
    return K_code_from_paper(round(K_paper * 2) / 2.0), {
        "K_paper_raw": K_paper,
        "sum_q": total,
        "quark_number": quark_number,
    }


# ──────────────────────────────────────────────────────────────────────────
# Large-N and analytic results
# ──────────────────────────────────────────────────────────────────────────

def thooft_rescale(N):
    """The ``(2 pi / N)^(1/2)`` factor on both axes of Fig. 8(b).

    Careful: this paper's natural unit is ``g^2 N / 2pi``, not 't Hooft's
    ``g^2 N / pi``.  That factor-of-2 redefinition of g^2 is the single most
    likely source of a spurious sqrt(2) when comparing against the standard
    't Hooft literature.
    """
    return np.sqrt(2.0 * np.pi / N)


def endpoint_exponent(mg, N, bracket=None):
    """Solve Eq. (26) for the endpoint exponent ``a`` in ``phi(x) ~ x^a``.

        1 - a pi cot(a pi) = pi m^2 / (g^2 (N^2 - 1) / (2N))

    ``a`` runs from 0 in the chiral limit (m/g = 0) to 1 in the free limit.
    It is needed for the Richardson series of Eq. (27), whose non-analytic
    terms go as ``K^-(1+a)`` and ``K^-(2+a)``.
    """
    from scipy.optimize import brentq

    if bracket is None:
        bracket = (1e-12, 1.0 - 1e-12)

    casimir = (N**2 - 1.0) / (2.0 * N)
    rhs = np.pi * mg**2 / casimir

    def f(a):
        return 1.0 - a * np.pi / np.tan(a * np.pi) - rhs

    if mg == 0.0:
        return 0.0
    return float(brentq(f, *bracket, xtol=1e-14, rtol=1e-14))


def meson_baryon_ratio_bosonization(N):
    """``2 sin[pi / (2(2N-1))]`` -- the strong-coupling M_mes/M_bar ratio.

    Derived by bosonization; the paper reports agreement to roughly 10%.
    N = 2, 3, 4 give 1.0, 0.618034, 0.445042.
    """
    return 2.0 * np.sin(np.pi / (2.0 * (2 * N - 1)))
