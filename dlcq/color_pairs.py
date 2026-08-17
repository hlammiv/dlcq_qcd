"""Color-resolved pair correlations: the diquark discriminator.

The color-blind pair correlator (:func:`dlcq.fock_weights.pair_correlation`)
separates molecular from democratic momentum sharing but cannot see diquark
substructure, because a diquark is a COLOR correlation: a quark pair in the
antisymmetric irrep (the antitriplet at N=3).  This module measures it.

For two quarks, ``T_i . T_j`` has exactly two eigenvalues,

    t_A = -(N+1)/(2N)   (antisymmetric pair irrep — the "diquark" channel)
    t_S = +(N-1)/(2N)   (symmetric),

so the diquark probability of a momentum-resolved pair is a linear readout,

    P_A(k1,k2) = (<T.T>(k1,k2)/C(k1,k2) - t_S) / (t_A - t_S),

with ``C`` the color-blind pair count.  On explicit color labels the operator
is the Fierz form ``T_i.T_j = (swap colors of i,j)/2 - 1/(2N)``, so the whole
computation runs on the brute-force expansions of ``tools/colour_norm.py``:
expand the c-weighted eigenstate over explicit ``(type, k, flavour, colour)``
labels, and for every quark pair look up the color-swapped label.

Two structural checks come free.  Every quark pair inside one
epsilon-contracted baryon cluster is EXACTLY antisymmetric, so a valence
baryon must return ``P_A = 1`` on every bin.  And for an all-quark singlet
sector of n quarks, ``sum_pairs <T.T> = -n C_F / 2`` times the sector
weight, an operator identity.

Brute force scales as the color expansion (N!^clusters per state), which is
tiny for the states this question is about; this is an analysis tool, not a
solver path.
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from .dataset import DLCQResult

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))


def _expansion(result: DLCQResult, s: int, N: int):
    from colour_norm import state_expansion

    L = int(result.state_len[s])
    types = result.state_types[s, :L]
    moms = result.state_moms[s, :L]
    flavs = result.state_flavors[s, :L]
    # cluster counts: quarks in baryons, antiquarks in antibaryons, mesons
    nq = int((types == 1).sum())
    nqbar = L - nq
    net = nq - nqbar
    nbrb = max(net // N, 0)
    nbrd = max(-net // N, 0)
    nmes = (L - N * (nbrb + nbrd)) // 2
    return state_expansion(types, moms, flavs, nbrb, nbrd, nmes, N)


def diquark_map(result: DLCQResult, state_idx: int = 0,
                nparton: int | None = None):
    """``(ks, C, TT)``: pair count and pair ``<T.T>``, momentum-resolved.

    Restricted to quark-quark pairs, optionally within one Fock sector.
    ``P_A = (TT/C - t_S)/(t_A - t_S)`` wherever ``C > 0``.
    """
    if result.c_orig is None or result.norm is None:
        raise ValueError("result lacks c_orig/norm")
    N = result.N
    c = result.require_eigenvector(state_idx)

    # c-weighted expansion of the (sector-restricted) eigenstate over
    # explicit labels.  Labels from different momentum configurations never
    # meet, so one flat dictionary is safe.
    D: dict = defaultdict(float)
    for s in range(result.numsta_post):
        if abs(c[s]) < 1e-14:
            continue
        L = int(result.state_len[s])
        if nparton is not None and L != nparton:
            continue
        for label, amp in _expansion(result, s, N).items():
            D[label] += c[s] * amp

    K = result.K_code
    C = np.zeros((K + 1, K + 1))
    TT = np.zeros((K + 1, K + 1))
    for label, amp in D.items():
        if amp == 0.0:
            continue
        ops = list(label)
        for i in range(len(ops)):
            if ops[i][0] != 1:
                continue
            for j in range(i + 1, len(ops)):
                if ops[j][0] != 1:
                    continue
                k1, k2 = ops[i][1], ops[j][1]
                # color-blind count
                C[k1, k2] += amp * amp
                C[k2, k1] += amp * amp
                # Fierz: T_i.T_j = swap/2 - 1/(2N)
                swapped = ops.copy()
                swapped[i] = (ops[i][0], ops[i][1], ops[i][2], ops[j][3])
                swapped[j] = (ops[j][0], ops[j][1], ops[j][2], ops[i][3])
                if len(set(swapped)) != len(swapped):
                    swap_amp = 0.0
                else:
                    order = sorted(range(len(swapped)),
                                   key=lambda t: swapped[t])
                    inv = sum(1 for a in range(len(order))
                              for b in range(a + 1, len(order))
                              if order[a] > order[b])
                    key = tuple(swapped[t] for t in order)
                    swap_amp = (-1 if inv % 2 else 1) * D.get(key, 0.0)
                val = amp * (0.5 * swap_amp - amp / (2.0 * N))
                TT[k1, k2] += val
                TT[k2, k1] += val
    ks = np.arange(1, K, 2)
    return ks, C[np.ix_(ks, ks)], TT[np.ix_(ks, ks)]


def diquark_fraction(result: DLCQResult, state_idx: int = 0,
                     nparton: int | None = None):
    """``(P_A_map, P_A_mean)``: bin-wise and pair-averaged diquark odds."""
    N = result.N
    t_A = -(N + 1) / (2.0 * N)
    t_S = (N - 1) / (2.0 * N)
    ks, C, TT = diquark_map(result, state_idx, nparton)
    with np.errstate(invalid="ignore", divide="ignore"):
        P = (TT / C - t_S) / (t_A - t_S)
    P[C <= 1e-12] = np.nan
    mean = float((TT.sum() / C.sum() - t_S) / (t_A - t_S))
    return ks, P, mean
