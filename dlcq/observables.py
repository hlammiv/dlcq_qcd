"""Physical observables extracted from a :class:`DLCQResult`.

Structure functions, sum rules, and the Richardson extrapolation that turns
finite-K spectra into the continuum numbers of Table I.

``structure_function`` is a **reconstruction of a lost program.**  The x-space
conversion in the original work was done by a separate code, ``wf``/``wfbig``,
which does not survive: ``qcdf.f`` emits only eigenvectors, Fock content and the
basis-change matrix, and refers to ``wf`` in comments.  Every structure function
in the paper passed through it.

That also means this is the one piece of the pipeline a solver-vs-solver
comparison cannot check -- the Python port reproduces the Fortran's matrix
elements to 1e-14, so an error here would be shared by both.  It is instead
validated against the thesis's own definition, the two sum rules (machine
precision, both solvers), a brute-force colour enumeration of the norm it uses,
and the published curves.  See docs/baryon-higher-fock.md section 1.3.

The structure function follows Eq. (12) of the paper:

    q(x) = K_paper <phi| b_k^dag b_k |phi>,     x = k / K_code

Basis states are number eigenstates, so ``n_k |s> = n_k(s) |s>`` and the
expectation value in the *non-orthogonal* Fock basis collapses to

    <phi| n_k |phi> = sum_s n_k(s) c_s (N c)_s        with  c^T N c = 1

which is why the weight is ``c_s * (N c)_s`` rather than ``c_s^2``.  Using
``c_s^2`` would silently break every sum rule.

The measure is ``dx = 2 / K_code = 1 / K_paper`` because momenta step by 2 in
code units (odd integers only).
"""

from __future__ import annotations

import numpy as np

from .dataset import DLCQResult
from .units import endpoint_exponent

__all__ = [
    "structure_function",
    "momentum_sum_rule",
    "number_sum_rule",
    "valence_parton_count",
    "richardson_extrapolate",
    "richardson_stability",
    "thooft_valence_limit",
    "spurious_zero_modes",
    "physical_indices",
    "color_slots_required",
    "fortran_overflow_risk",
    "FORTRAN_COLOR_SLOTS",
]


def color_slots_required(max_partons):
    """Color-index slots a matrix element needs: ``2 * max_partons + 4``.

    Bra and ket each contribute their partons, plus up to four operators from a
    four-point vertex.  ``qcdf.f`` dimensions ``IDELT(12552, 25)`` and sets
    ``MXLNG = 25``, so any run reaching 11 partons overflows that array.
    ``qcdf_opt.py`` uses ``MXLNG = 2*MXP + 4 = 54`` and is unaffected.
    """
    return 2 * int(max_partons) + 4


FORTRAN_COLOR_SLOTS = 25


def fortran_overflow_risk(result):
    """Whether this configuration exceeds ``qcdf.f``'s 25-slot color array."""
    if result.state_len is None or len(result.state_len) == 0:
        return False
    return color_slots_required(int(np.max(result.state_len))) > FORTRAN_COLOR_SLOTS


def spurious_zero_modes(result, tol=1e-12):
    """Indices of unphysical eigenstates -- non-positive M^2 at non-zero mass.

    ``M^2 >= 0`` for any physical state, and the paper states that an exactly
    massless state occurs *only* in the chiral limit: "when m/g = 0 identically,
    the lightest state for any N or B is exactly zero, independently of K".  So
    at m/g > 0 any eigenvalue at or below zero is an artifact.

    These are not hypothetical, and they sit at the *bottom* of the spectrum --
    precisely where Figs. 7, 8 and Table I read off the lightest state.  Their
    cause is an array-bounds overflow in ``qcdf.f``: a matrix element between
    two L-parton states needs ``2L + 4`` color-index slots, but ``IDELT`` is
    dimensioned with only 25, so runs reaching L >= 11 corrupt the color
    contraction.  Observed:

    ===============  ====  ======  ==========================================
    run              L     slots   symptom
    ===============  ====  ======  ==========================================
    2K=21, B=1        9     22     clean
    2K=10, B=0        6     16     clean
    2K=24, B=0       12     28     one decoupled 12-parton state at M^2 = 0
    2K=25, B=1       11     26     four negative M^2, all 11-parton dominated
    ===============  ====  ======  ==========================================

    See ``docs/fortran-color-overflow.md``.  The Python solver is unaffected.

    Returns an array of indices into ``result.eigenvalues``.
    """
    if result.eigenvalues is None or result.eigenvalues.size == 0:
        return np.array([], dtype=int)
    # rlamb == 1 is the chiral limit, where a zero mode is genuine.
    if result.rlamb >= 1.0 - 1e-12:
        return np.array([], dtype=int)
    return np.flatnonzero(result.eigenvalues <= tol)


def physical_indices(result, tol=1e-12):
    """Eigenstate indices with the spurious zero modes removed.

    Use this wherever the paper says "the first three states"; indexing
    ``eigenvalues`` directly can silently pick up a decoupled mode.
    """
    bad = set(spurious_zero_modes(result, tol).tolist())
    n = result.n_eigenvalues
    return np.array([i for i in range(n) if i not in bad], dtype=int)


def valence_parton_count(N: int, B: int) -> int:
    """Partons in the minimal (valence) Fock sector: 2 for a meson, N*|B| otherwise."""
    return 2 if B == 0 else N * abs(B)


def structure_function(result: DLCQResult, state_idx: int = 0,
                       nparton: int | None = None):
    """Quark and antiquark structure functions for one eigenstate.

    Parameters
    ----------
    result
        A parsed run.  Needs ``c_orig``, ``norm``, and the basis arrays.
    state_idx
        Which eigenstate, 0 = ground state.
    nparton
        Restrict to the Fock sector with exactly this many partons (the
        decomposition plotted in Figs. 4-6).  ``None`` sums all sectors.

    Returns
    -------
    x, q_quark, q_antiquark
        ``x`` is the odd-momentum grid ``k / K_code``.  The measure that goes
        with it is ``dx = 2 / K_code``.
    """
    if result.c_orig is None or result.norm is None:
        raise ValueError("result lacks c_orig/norm; parse with with_matrices=True")

    c = result.require_eigenvector(state_idx)
    Nc = result.norm @ c
    weights = c * Nc                      # sums to 1 over all states

    K = result.K_code
    K_paper = result.K_paper

    q = np.zeros(K + 1)
    qbar = np.zeros(K + 1)

    for s in range(result.numsta_post):
        if nparton is not None and result.state_len[s] != nparton:
            continue
        w = weights[s]
        if w == 0.0:
            continue
        L = result.state_len[s]
        for p in range(L):
            k = result.state_moms[s, p]
            if result.state_types[s, p] == 1:
                q[k] += w
            else:
                qbar[k] += w

    ks = np.arange(1, K, 2)
    x = ks / float(K)
    return x, K_paper * q[ks], K_paper * qbar[ks]


def momentum_sum_rule(result: DLCQResult, state_idx: int = 0) -> float:
    """``integral x [q(x) + qbar(x)] dx``, which must equal 1 exactly.

    Every basis state carries total momentum ``K_code``, so this holds
    state-by-state and is independent of the basis conditioning -- the
    strongest available check that the structure function is built correctly.
    """
    x, q, qbar = structure_function(result, state_idx)
    dx = 2.0 / result.K_code
    return float(np.sum(x * (q + qbar)) * dx)


def number_sum_rule(result: DLCQResult, state_idx: int = 0) -> float:
    """``integral [q(x) - qbar(x)] dx``, which must equal ``N * B``."""
    x, q, qbar = structure_function(result, state_idx)
    dx = 2.0 / result.K_code
    return float(np.sum(q - qbar) * dx)


# ──────────────────────────────────────────────────────────────────────────
# Continuum extrapolation
# ──────────────────────────────────────────────────────────────────────────

def richardson_extrapolate(K_codes, masses, mg, N, n_terms=4, return_fit=False):
    """Extrapolate M(K) to the continuum using Eq. (27).

        M(1/K) = M(0) + c1/K + c2/K^(1+a) + c3/K^2 + c4/K^(2+a) + ...

    The non-analytic exponents come from the endpoint behaviour phi(x) ~ x^a of
    Eq. (26), with ``a`` solved per (N, m/g).

    This replaces the heuristic previously used in ``reproduce_figures.py``
    (take whichever K gives the highest positive lightest eigenvalue), which is
    not what the paper does and cannot reproduce Table I.

    Parameters
    ----------
    K_codes, masses
        Matched arrays. ``K_codes`` in code units (= 2 * paper K); the fit uses
        the paper's K internally, matching the paper's "2K in the range 16-24".
    mg, N
        Needed for the endpoint exponent.
    n_terms
        Number of correction terms. The paper estimates the error in M(0) as
        the magnitude of the last retained term -- that is what Table I's
        parentheses hold, so it is returned rather than a statistical error.

    Returns
    -------
    M0, last_term
        Continuum mass and the paper's error estimate.  With ``return_fit``,
        also the coefficient vector and the exponents used.
    """
    K_codes = np.asarray(K_codes, dtype=float)
    masses = np.asarray(masses, dtype=float)

    order = np.argsort(K_codes)
    K_codes, masses = K_codes[order], masses[order]

    a = endpoint_exponent(mg, N)
    # Exponents in 1/K_paper, ordered by how fast they vanish.
    exponents = [0.0, 1.0, 1.0 + a, 2.0, 2.0 + a, 3.0][: n_terms + 1]
    # a == 0 (chiral limit) makes 1/K^(1+a) degenerate with 1/K; drop dupes.
    seen, uniq = set(), []
    for e in exponents:
        key = round(e, 12)
        if key not in seen:
            seen.add(key)
            uniq.append(e)
    exponents = uniq

    if len(K_codes) < len(exponents):
        exponents = exponents[: len(K_codes)]

    Kp = K_codes / 2.0                       # the paper's K
    A = np.vstack([Kp ** (-e) for e in exponents]).T
    coeffs, *_ = np.linalg.lstsq(A, masses, rcond=None)

    M0 = float(coeffs[0])
    # "the numbers in parentheses give the magnitude of the last term in the
    # series fit", evaluated at the largest K used.
    last_term = float(abs(coeffs[-1]) * Kp[-1] ** (-exponents[-1]))

    if return_fit:
        return M0, last_term, coeffs, exponents
    return M0, last_term


def richardson_stability(K_codes, masses, mg, N, n_terms=2, min_points=4):
    """How much ``M(0)`` moves when the fit window is changed.

    The paper reports "the magnitude of the last term in the series fit" as its
    uncertainty, which :func:`richardson_extrapolate` returns.  That rule breaks
    down at weak coupling, and not because K is too small.  The Eq. (27) basis
    is ``1, 1/K, 1/K^(1+a), ...`` with the Eq. (26) exponent ``a``, and ``a -> 0``
    in the chiral limit, so ``1/K^(1+a)`` collapses onto ``1/K``.  The fit stays
    excellent -- maximum residual 1e-5 or better at every coupling, five digits
    on M(K) -- but the coefficients blow up with near-total cancellation: 0.892
    and 0.098 at m/g = 1.6, against 8.383 and 8.800 at m/g = 0.1.  The last term
    is then large because the basis is degenerate, not because the answer is
    uncertain.

    This refits on every sub-window of at least ``min_points`` points and
    returns the standard deviation of the resulting ``M(0)``.  Measured against
    the last-term rule for the N = 3 baryon over 2K = 25-35, the latter
    overstates by 6x at m/g = 1.6 rising to **341x at m/g = 0.1**, where the
    extrapolation is stable to 0.2% rather than the 72% the rule reports.

    **What this does not measure.**  It is the fit's sensitivity to the window,
    not a total uncertainty.  Systematics that shift every point in the same
    direction are invisible to it -- in particular the ``sweep_lpn``
    particle-number truncation, and any finite-K effect the Eq. (27) form does
    not capture.  At weak coupling the residual disagreement with the published
    values is larger than either estimate, so something systematic remains.

    Returns ``(M0, spread, n_windows)``; ``spread`` is NaN when there are too
    few points to form more than one sub-window.
    """
    from itertools import combinations

    K_codes = list(K_codes)
    masses = list(masses)
    M0, _ = richardson_extrapolate(K_codes, masses, mg, N, n_terms=n_terms)
    n = len(K_codes)
    fits = []
    for size in range(min_points, n + 1):
        for idx in combinations(range(n), size):
            if size == n:
                continue
            sub = richardson_extrapolate([K_codes[i] for i in idx],
                                         [masses[i] for i in idx],
                                         mg, N, n_terms=n_terms)[0]
            fits.append(sub)
    if len(fits) < 2:
        return M0, float("nan"), len(fits)
    return M0, float(np.std(fits)), len(fits)


def thooft_valence_limit(x, N, B=1):
    """Analytic chiral-limit baryon structure function, Eq. (22).

        q(x) = N (N-1) (1-x)^(N-2)

    N=3 gives 6(1-x); N=2 gives the constant 2, identical to the meson
    distribution.  A hard target for Fig. 3(b) as m/g -> 0.
    """
    x = np.asarray(x, dtype=float)
    return N * (N - 1) * (1.0 - x) ** (N - 2)
