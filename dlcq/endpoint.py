"""van de Sande's endpoint subtraction for the DLCQ self-induced inertia.

Standard DLCQ converges as ``K^-2b`` at small quark mass, where ``b`` is the
Eq. (26) endpoint exponent, because the bound-state equation never uses the
known ``phi ~ x^b`` behaviour of the wavefunction near ``x = 0`` and ``x = 1``.
Since ``b -> 0`` in the chiral limit that is arbitrarily slow, and it is why
Table I's weak-coupling column measures a grid artifact rather than the
continuum limit -- see docs/weak-coupling-limit.md.

Brett van de Sande, `hep-ph/9605409 <https://arxiv.org/abs/hep-ph/9605409>`_,
fixes it by adding and subtracting a term so the kernel vanishes when *both*
momenta are near an endpoint.  With ``w(x) = (x(1-x))^b`` his Eq. (10) is

    M^2 Psi(x) = mu^2 Psi(x) [1/x + 1/(1-x)]
               + g^2 Psi(x) I(x)
               + g^2 int dy [Psi(x) w(y)/w(x) - Psi(y)] / (x-y)^2

    I(x) = PV int_0^1 dy [1 - w(y)/w(x)] / (x-y)^2

Adding the two ``g^2`` pieces reproduces ``int [Psi(x) - Psi(y)]/(x-y)^2``
exactly, so this is an identity in the continuum and differs only under
discretisation.  The **exchange term is untouched**; the whole modification is
to the self-induced inertia.

Why this module is so small
---------------------------
In the two-particle sector the partner momentum is fixed by conservation --
``k`` pairs with ``K_code - k`` and nothing else -- so the improved self-inertia
is still a function of a *single* momentum, exactly like ``compute_selfen``'s
table.  The improvement is therefore a **different ``selfen`` array**, not a
change to the kernel: ``qcdf_kernels.hamqcd_nb`` takes ``selfen`` as an
argument and never needs to know.

That is why :func:`improved_selfen` returns something with the same shape and
units as ``qcdf_opt.compute_selfen``, and why ``hamiltonian="improved"`` is
restricted to sectors where every state has two partons.  For three or more
partons the subtraction genuinely depends on which partner a parton is paired
with, and the reduction that makes this work no longer holds -- measured, the
identity ``diag(self-energy) = -sum_j H_exchange[i,j]`` is exact to 1e-16 for
two partons and fails by 8.3e-01 for three.

Validation
----------
Against van de Sande's own numbers at ``b = 0.1`` (his ``mu = 0.181981``), where
the exact ground state is ``M^2/g^2 = 0.779141``: standard DLCQ reaches only
0.413 at K = 100, while improved DLCQ gives 0.771 at K = **10**.  Ordinary
Richardson in ``1/K`` on the improved series returns 0.779315, i.e. the exact
answer to 2e-4 -- his Eq. (14) -- where the same fit on the standard series
returns 0.4536.  The chiral exponent, 2 in standard DLCQ and 1 in reality
(his Eq. 7, ``M^2 = 2 pi g mu / sqrt 3``), comes out 1.05.
"""
from __future__ import annotations

import numpy as np

__all__ = ["endpoint_integral", "improved_selfen"]


def endpoint_integral(x, b, limit=400):
    """``PV int_0^1 dy [1 - (y(1-y)/(x(1-x)))^b] / (x-y)^2`` -- van de Sande's I(x).

    Evaluated as a quadrature rather than approximated by the DLCQ sum, which
    is the entire point: this is the piece that carries the endpoint region,
    and computing it exactly is what removes the ``K^-2b`` error.

    The numerator vanishes linearly at ``y = x``, so the integrand has a simple
    pole and the principal value is finite.  Writing ``g(y) = 1 - w(y)/w(x)``,

        int g(y)/(x-y)^2 dy = int [g(y)/(y-x)] / (y-x) dy

    which is exactly SciPy's ``weight="cauchy"`` applied to ``g/(y-x)``.  That
    factor has a removable singularity at ``y = x``, filled analytically with
    ``g'(x) = -w'(x)/w(x) = -b(1-2x)/(x(1-x))``.
    """
    from scipy.integrate import quad

    if b == 0.0:
        return 0.0
    wx = (x * (1.0 - x)) ** b
    at_x = -b * (1.0 - 2.0 * x) / (x * (1.0 - x))

    def f(y):
        dy = y - x
        if abs(dy) < 1e-9:
            return at_x
        return (1.0 - ((y * (1.0 - y)) ** b) / wx) / dy

    val, _ = quad(f, 0.0, 1.0, weight="cauchy", wvar=x, limit=limit)
    return float(val)


def improved_selfen(N, K_code, b, mxslfn=100):
    """``selfen`` with the endpoint subtraction applied, for the two-parton sector.

    Drop-in replacement for ``qcdf_opt.compute_selfen(N)``: same shape, same
    indexing (``selfen[k-1]`` for a parton of odd momentum ``k``), same units.

    Construction.  The standard table satisfies, for the pair ``(k, K-k)``,

        selfen[k-1] + selfen[K-k-1] = cfact * sum_{j odd != k} 4/(k-j)^2

    because ``k`` and ``j`` are both odd, so ``|k-j|`` runs over the even
    integers ``{2..k-1} u {2..K-k-1}`` -- the two one-body self-inertias.  The
    improvement inserts the weight ``w_j/w_k`` into that sum and adds ``I(x)``:

        W(k)  = cfact * sum_{j odd != k} (w_j/w_k) * 4/(k-j)^2
        Ic(k) = cfact * I(k/K_code) / K_paper

    Both are symmetric under ``k -> K_code - k``, so splitting them evenly
    between the pair's two partons reproduces the code's per-parton
    accumulation:

        improved[k-1] = (W(k) + Ic(k)) / 2

    ``I`` carries no factor of K in Eq. (12) while the sum carries one, and the
    solver multiplies ``elem`` by ``K_paper`` on the way to ``M^2``
    (``eigenvalues = K_code * w / 2``), hence the ``/ K_paper`` on ``Ic``.

    At ``b = 0`` the weight is 1 and ``I`` vanishes, so this returns the
    standard table exactly -- the reduction that
    ``tests/test_endpoint.py`` pins.
    """
    b = float(b)
    K_code = int(K_code)
    K_paper = K_code / 2.0
    cfact = (N * N - 1.0) / N

    out = np.zeros(mxslfn)
    odd = np.arange(1, K_code, 2)
    if odd.size == 0:
        return out

    x = odd / float(K_code)
    w = (x * (1.0 - x)) ** b if b else np.ones_like(x)

    for idx, k in enumerate(odd):
        if k - 1 >= mxslfn:
            continue
        d = (k - odd).astype(float)
        keep = d != 0.0
        wsum = float(np.sum((w[keep] / w[idx]) * 4.0 / d[keep] ** 2))
        ic = (endpoint_integral(float(x[idx]), b) / K_paper) if b else 0.0
        out[k - 1] = cfact * (wsum + ic) / 2.0
    return out
