"""Per-state Fock-sector weights, and correlated binding energies.

The exotics questions — is any B=0 state *dominated* by the four-parton
sector, is the B=2 state bound against two baryons — need two primitives
that deliberately live outside :mod:`dlcq.observables`:

``fock_weights``
    The probability an eigenstate carries in each parton-number sector.
    The norm is exactly block diagonal in parton configuration (verified by
    brute-force colour enumeration, ``tools/colour_norm.py``), so the
    ``c * (N c)`` weight that :func:`dlcq.observables.structure_function`
    already uses splits exactly by sector: summing it over the states of one
    parton count *is* ``c_L^T N_L c_L``, with no cross terms to argue away.

``binding_series``
    ``Delta(K) = M_state(K) - n * M_threshold(K)`` from same-K correlated
    differences.  Subtracting two *separately extrapolated* masses double
    counts every K-systematic the two channels share; the correlated
    difference cancels them instead, and the difference series is what gets
    extrapolated.  The function refuses mismatched K grids and mismatched
    Hamiltonians outright — at N=3 the B=2 grid (even 2K) and the single-
    baryon grid (odd 2K) can never match, which is a physics fact to design
    around (run the hexaquark study at N=4, or bracket and say so), not to
    paper over with interpolation here.
"""
from __future__ import annotations

import numpy as np

from .dataset import DLCQResult


def fock_weights(result: DLCQResult, state_idx: int = 0) -> dict[int, float]:
    """Weight of each parton-number sector in one eigenstate.

    Returns ``{parton_count: weight}`` with the weights summing to 1 (the
    eigenvector is norm-normalized).  Needs ``c_orig`` and ``norm``, like
    :func:`dlcq.observables.structure_function`.
    """
    if result.c_orig is None or result.norm is None:
        raise ValueError("result lacks c_orig/norm; parse with with_matrices=True")
    c = result.require_eigenvector(state_idx)
    w = c * (result.norm @ c)
    lengths = np.asarray(result.state_len[:result.numsta_post], dtype=int)
    out: dict[int, float] = {}
    for L in np.unique(lengths):
        out[int(L)] = float(w[lengths == L].sum())
    return out


def dominant_sector(result: DLCQResult, state_idx: int = 0) -> tuple[int, float]:
    """The parton count carrying the largest weight, and that weight."""
    weights = fock_weights(result, state_idx)
    L = max(weights, key=weights.get)
    return L, weights[L]


def binding_series(K_state, M_state, K_thresh, M_thresh, n_thresh: int = 2,
                   hamiltonian_state: str | None = None,
                   hamiltonian_thresh: str | None = None):
    """Correlated binding-energy series ``Delta(K) = M_state - n * M_thresh``.

    Parameters
    ----------
    K_state, M_state
        The candidate's mass series (K in code units, matched arrays).
    K_thresh, M_thresh
        The threshold constituent's series **on the same K grid**.
    n_thresh
        How many constituents the threshold holds (2 for two mesons or two
        baryons, adjust for meson+baryon by passing the summed series with
        ``n_thresh=1``).
    hamiltonian_state, hamiltonian_thresh
        Pass both when known; a standard/improved mix is refused because the
        two solve different operators and their difference is meaningless.

    Returns
    -------
    K, delta
        The common grid and the difference series, ready for extrapolation
        (in the fit basis matching the Hamiltonian — plain 1/K for improved,
        the Eq. 27 ladder for standard; never a mix).
    """
    if hamiltonian_state is not None or hamiltonian_thresh is not None:
        if hamiltonian_state != hamiltonian_thresh:
            raise ValueError(
                f"Hamiltonian mismatch: state={hamiltonian_state!r} vs "
                f"threshold={hamiltonian_thresh!r}. A binding energy across "
                f"operators is meaningless.")
    K_state = np.asarray(K_state, dtype=int)
    K_thresh = np.asarray(K_thresh, dtype=int)
    if K_state.shape != K_thresh.shape or not np.array_equal(K_state, K_thresh):
        raise ValueError(
            f"K grids differ (state {K_state.tolist()} vs threshold "
            f"{K_thresh.tolist()}): a correlated difference needs identical "
            f"grids. At N=3 the B=2 and B=1 parities are incompatible — use "
            f"N=4 for exact differences, or bracket explicitly and label it.")
    M_state = np.asarray(M_state, dtype=float)
    M_thresh = np.asarray(M_thresh, dtype=float)
    return K_state.copy(), M_state - n_thresh * M_thresh
