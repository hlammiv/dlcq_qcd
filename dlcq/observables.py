"""Physical observables extracted from a :class:`DLCQResult`.

Structure functions, sum rules, and the Richardson extrapolation that turns
finite-K spectra into the continuum numbers of Table I.

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
    "thooft_valence_limit",
]


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


def thooft_valence_limit(x, N, B=1):
    """Analytic chiral-limit baryon structure function, Eq. (22).

        q(x) = N (N-1) (1-x)^(N-2)

    N=3 gives 6(1-x); N=2 gives the constant 2, identical to the meson
    distribution.  A hard target for Fig. 3(b) as m/g -> 0.
    """
    x = np.asarray(x, dtype=float)
    return N * (N - 1) * (1.0 - x) ** (N - 2)
