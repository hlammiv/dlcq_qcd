#!/usr/bin/env python3
"""Compute the Fock-basis norm matrix by explicit colour summation.

`qcdf.f` builds the norm diagrammatically.  Its own header records when that
happened::

    **** 6/24/88 MODIFIED QCD2A2 SO THAT COLOR ****
    **** SUMS ARE PERFORMED DIAGRAMMATICALLY   ****
    **** RATHER THAN ITERATIVELY               ****

The thesis is 1988 and the article 1990, so the published higher-Fock curves may
predate that rewrite.  The pre-1988 path is gone from the file, but the quantity
it computed is not ambiguous: the norm is a sum over colour indices, and for
N = 3 with at most ten partons that sum is small enough to do by brute force.
This module does exactly that, which makes it an independent check of the
diagrammatic result rather than a reimplementation of it.

The basis states are products of colour-singlet clusters, laid out by
``qcdsta`` in a fixed order:

    partons [0, N*nbrb)                  nbrb baryons, N quarks each,
                                         antisymmetrized with epsilon
    partons [N*nbrb, N*(nbrb+nbrd))      nbrd antibaryons, likewise
    partons [N*(nbrb+nbrd), ... )        nmes mesons, a quark and an antiquark
                                         contracted with delta

so a state is

    |s> = prod_baryons  sum eps_{a1..aN} b+(k1,f1,a1) ... b+(kN,fN,aN)
        x prod_antibar. sum eps_{a1..aN} d+(...) ...
        x prod_mesons   sum_a b+(k,f,a) d+(k',f',a)   |0>

Expanding each cluster gives N! signed terms (epsilon) or N terms (delta).  Each
term is a product of creation operators in the order above; sorting that product
into a canonical order gives a Fock basis label and a fermionic sign.  The
overlap is then just the inner product of the two expansions.
"""

from __future__ import annotations

import itertools
from collections import defaultdict

import numpy as np


def _levi_civita_terms(n):
    """(colour tuple, sign) for every non-zero term of the rank-n epsilon."""
    out = []
    for perm in itertools.permutations(range(n)):
        sign = 1
        seen = list(perm)
        # parity by counting inversions
        inv = sum(1 for i in range(n) for j in range(i + 1, n) if seen[i] > seen[j])
        sign = -1 if inv % 2 else 1
        out.append((perm, sign))
    return out


def state_expansion(types, moms, flavs, nbrb, nbrd, nmes, N):
    """Expand one basis state over the (type, k, flavour, colour) Fock basis.

    Returns ``{canonical_label: amplitude}``.  The label is a sorted tuple of
    ``(type, momentum, flavour, colour)``; the amplitude carries the fermionic
    sign of sorting the operators into that order.
    """
    clusters = []
    p = 0
    for _ in range(nbrb + nbrd):
        clusters.append(("eps", list(range(p, p + N))))
        p += N
    for _ in range(nmes):
        clusters.append(("delta", [p, p + 1]))
        p += 2

    eps = _levi_civita_terms(N)
    choices = []
    for kind, slots in clusters:
        if kind == "eps":
            choices.append([(dict(zip(slots, colours)), sign) for colours, sign in eps])
        else:
            choices.append([({slots[0]: c, slots[1]: c}, 1) for c in range(N)])

    out = defaultdict(float)
    for combo in itertools.product(*choices):
        colour = {}
        amp = 1
        for assign, sign in combo:
            colour.update(assign)
            amp *= sign
        ops = [(int(types[i]), int(moms[i]), int(flavs[i]), colour[i])
               for i in range(p)]
        if len(set(ops)) != len(ops):        # Pauli: repeated operator
            continue
        order = sorted(range(len(ops)), key=lambda i: ops[i])
        inv = sum(1 for i in range(len(order)) for j in range(i + 1, len(order))
                  if order[i] > order[j])
        out[tuple(ops[i] for i in order)] += amp * (-1 if inv % 2 else 1)
    return out


def norm_bruteforce(states, numsta, N):
    """Gram matrix of the basis, by explicit colour summation."""
    exps = []
    for s in range(numsta):
        loc = states.mstinf[s, 0] - 1
        L = states.mstinf[s, 1]
        nmes, nbrb, nbrd = (states.mstinf[s, 2], states.mstinf[s, 3],
                            states.mstinf[s, 4])
        types = states.mstate[loc, :L]
        moms = states.mstate[loc + 2, :L]
        flavs = states.mstate[loc + 3, :L]
        exps.append(state_expansion(types, moms, flavs, nbrb, nbrd, nmes, N))

    G = np.zeros((numsta, numsta))
    for s in range(numsta):
        for t in range(s, numsta):
            a, b = exps[s], exps[t]
            if len(b) < len(a):
                a, b = b, a
            v = sum(amp * b.get(key, 0.0) for key, amp in a.items())
            G[s, t] = G[t, s] = v
    return G
