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
restricted to sectors where every state has two partons: with three or more,
``sigma`` depends on *which partner* a parton is paired with, and a table
indexed by a single momentum cannot express that.  The restriction is a limit
of this representation, not of the physics.

**It is not that anything breaks at three partons.**  An earlier version of
this module said the underlying identity failed there, on the strength of a
test that compared the self-energy against an *unweighted* row sum of the
exchange.  That test was wrong: it is only valid where the norm matrix is the
identity, which happens to be true in the meson valence sector and nowhere
else.  What actually holds, measured to 1e-16 relative at L = 2, 3, 4, 5 and in
mixed sectors up to L = 8 including the untruncated basis, is the
*norm-weighted* statement

    D  ==  diag(sigma_std) @ Norm,    sigma_std(s) = C_F sum_{partons p} S(k_p)

with ``C_F = (N^2-1)/(2N)`` and ``S(k) = sum_{n=1}^{(k-1)/2} 1/n^2``.  So the
self-energy stays a one-body scalar at every L, and ``Norm_ij`` is non-zero
only between states of identical parton content.  Extending to L >= 3 therefore
needs a partner-indexed table passed alongside ``selfen``, not a different
``selfen``; see ``docs/next-steps.md``.

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


def directed_pair_kernel(k, l, b):
    """``J^(k; l)`` -- one parton's share of the endpoint sum for a pair.

    ``k`` and ``l`` are odd code momenta; ``P = k + l`` is the pair's total.

        J^(k; l) = sum_{j odd, 0 < j < k} (j(P-j)/(k l))^b * 4/(k-j)^2
                 + I_b(k/P) / P_paper

    Directed rather than symmetric: ``J^(k;l) + J^(l;k)`` is van de Sande's
    symmetric pair value ``J(k,l)``, and splitting it this way is what lets a
    parton's share be assigned to the parton itself.

    At ``b = 0`` the weight is 1 and, since ``k`` and ``j`` are both odd,
    ``k - j`` runs over the even integers ``2..k-1``, so the sum collapses to
    ``S(k) = sum_{n=1}^{(k-1)/2} 1/n^2`` -- *independent of the partner*.  That
    is the reduction that makes this a safe drop-in, and it is exact rather
    than asymptotic.
    """
    k = int(k)
    l = int(l)
    P = k + l
    js = np.arange(1, k, 2, dtype=float)          # odd, 0 < j < k
    if b == 0.0:
        s = float(np.sum(4.0 / (k - js) ** 2)) if js.size else 0.0
        return s
    w = (js * (P - js) / (k * l)) ** b
    s = float(np.sum(w * 4.0 / (k - js) ** 2)) if js.size else 0.0
    # Half the counterterm per direction.  ``I`` carries no factor of K in
    # Eq. (12) while the sum does, and the solver multiplies the interaction by
    # K_paper on the way to M^2, so a *pair* wants ``I / P_paper``.  This is one
    # of the two directions and ``I_b`` is symmetric about 1/2, so
    # ``J^(k;l) + J^(l;k)`` restores it: 2 * I/P_code = I/P_paper.
    return s + endpoint_integral(k / float(P), b) / float(P)


def pair_kernel_table(K_code, b):
    """``J^(k; l)`` for every odd ``k, l`` with ``k + l <= K_code``.

    Indexed ``[k, l]`` so the kernel can look it up by raw momentum.  Built
    once per (K_code, b): ~1250 entries at 2K = 99, each a short sum plus one
    quadrature.  Never evaluate the principal value inside a matrix loop.
    """
    K_code = int(K_code)
    T = np.zeros((K_code + 1, K_code + 1))
    for k in range(1, K_code, 2):
        for l in range(1, K_code - k + 1, 2):
            T[k, l] = directed_pair_kernel(k, l, b)
    return T


def apply_sigma_correction(ham, norm, delta):
    """``ham + Norm @ diag(delta)``, preserving whatever container each is in.

    The self-energy enters the Hamiltonian as ``Norm @ diag(sigma)`` -- measured
    exactly, at every parton number -- so swapping ``sigma_std`` for
    ``sigma_imp`` is this one addition.  ``Norm_ij`` is non-zero only between
    states of identical parton content, on which ``sigma`` is constant, so the
    result stays symmetric and left/right multiplication coincide.

    ``norm`` may be a dense array or a :class:`~dlcq.dataset.BlockDiagonal`;
    ``ham`` may be dense or sparse.  Blocks are kept blocks: densifying the norm
    is what the block container exists to avoid (151 GB at LPN=7, 2K=71).
    """
    from scipy import sparse

    delta = np.asarray(delta, dtype=float)
    rows = getattr(norm, "rows", None)
    if rows is not None:                       # BlockDiagonal
        data, ri, ci = [], [], []
        for idx, blk in zip(norm.rows, norm.blocks):
            idx = np.asarray(idx)
            sub = blk * delta[idx][np.newaxis, :]
            rr, cc = np.nonzero(sub)
            if rr.size:
                data.append(sub[rr, cc])
                ri.append(idx[rr])
                ci.append(idx[cc])
        n = norm.n_rows
        corr = sparse.csr_matrix(
            (np.concatenate(data) if data else np.zeros(0),
             (np.concatenate(ri) if ri else np.zeros(0, int),
              np.concatenate(ci) if ci else np.zeros(0, int))), shape=(n, n))
    else:
        corr = np.asarray(norm) * delta[np.newaxis, :]

    if sparse.issparse(ham):
        return (ham + sparse.csr_matrix(corr)).tocsr()
    if sparse.issparse(corr):
        corr = corr.toarray()
    return ham + corr


def state_sigmas(mstate, mstinf, numsta, N, b, K_code, table=None):
    """``(sigma_std, sigma_imp)`` per basis state, for the matrix-level form.

    The standard self-inertia is a one-body scalar at every parton number --
    measured, ``H(selfen) - H(selfen=0) == diag(sigma_std) @ Norm`` to 1e-16 for
    L = 2..8 including untruncated bases::

        sigma_std(s) = C_F * sum_p S(k_p)

    The improvement replaces each parton's ``S(k_a)`` by its share of the
    endpoint sums with every partner::

        sigma_imp(s) = sum_a  C_F/(L-1) * sum_{c != a} J^(k_a; k_c)

    The ``C_F/(L-1)`` weight is not an arbitrary choice wherever the pairs are
    colour-equivalent, which covers every valence sector.  For a colour singlet
    ``sum_a T_a = 0`` gives ``sum_{c != a} (-T_a . T_c) = C_F`` per parton, and
    when all pairs sit in the same channel that forces ``C_F/(L-1)`` each.
    Checked against the measured exchange coefficients: ``c_qq = (N+1)/2N`` for
    N = 2..6 equals ``C_F/(L-1)`` at ``L = N``.

    Where pairs are *not* colour-equivalent -- Fock-extended sectors carrying
    more than one colour singlet, which is exactly Table I's ``LPN =
    valence+2`` -- the scalar stands in for the matrix ``-T_a . T_c`` and is an
    ansatz.  Measured effect there: ~3e-3 in M/g.  See docs/next-steps.md.

    ``sigma_imp`` reduces to ``sigma_std`` exactly at ``b = 0``, per parton, not
    merely per state, because ``J^(k;l)|_{b=0} = S(k)`` independently of ``l``.
    """
    if table is None:
        table = pair_kernel_table(K_code, b)
    CF = (N * N - 1.0) / (2.0 * N)
    std = np.zeros(numsta)
    imp = np.zeros(numsta)
    for s in range(numsta):
        loc = int(mstinf[s, 0]) - 1
        L = int(mstinf[s, 1])
        ks = [int(mstate[loc + 2, j]) for j in range(L)]
        std[s] = CF * sum(_S(k) for k in ks)
        if L < 2:
            imp[s] = std[s]
            continue
        tot = 0.0
        for a, ka in enumerate(ks):
            share = sum(table[ka, kc] for c, kc in enumerate(ks) if c != a)
            tot += share / (L - 1.0)
        imp[s] = CF * tot
    return std, imp


def _S(k):
    """``sum_{n=1}^{(k-1)/2} 1/n^2`` -- the standard one-body self-inertia sum."""
    return float(sum(1.0 / n ** 2 for n in range(1, (int(k) - 1) // 2 + 1)))


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
