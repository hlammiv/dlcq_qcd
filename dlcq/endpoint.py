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


#: Vertex patterns for the *number-conserving* instantaneous exchange, read off
#: ``qcdf.hamqcd``.  Each entry is (quark occupations, antiquark occupations,
#: clfact index pattern) for the four operator slots.
#:
#: H3-H6 are pair creation/annihilation, not scattering, so they carry no
#: self-inertia partner and are absent here.  H7 has two flavour structures; the
#: one that pairs with the self-energy is the *t-channel* A
#: (``ka = k4-k3, kb = k2-k1``, vanishing at zero momentum transfer), not the
#: annihilation structure B, which survives there but is a different graph.
_VERTEX = {
    "qq":   ([-1, -1, 1, 1], [0, 0, 0, 0], (0, 3, 1, 2)),      # H1
    "qbqb": ([0, 0, 0, 0], [-1, -1, 1, 1], (2, 1, 3, 0)),      # H2
    "qqb":  ([0, 0, -1, 1], [-1, 1, 0, 0], (1, 3, 2, 0)),      # H7 structure A
}


def pair_colour_weights(base, params, mstate, loc, L, selfen):
    """``<-T_a . T_c>`` for every parton pair of one state, normalised to C_F.

    Returns ``w[a][c]`` with ``sum_{c != a} w[a][c] == C_F`` exactly, which is
    what ``sum_a T_a = 0`` requires of a colour singlet and what makes the
    ``b = 0`` reduction to ``sigma_std`` automatic.

    No new colour machinery: ``clfact``'s ``nops=4`` path already *is*
    ``sum_a T^a_ij T^a_kl`` (``nct=1`` carries the 0.5, ``nct=2`` the -0.5/N).
    Calling it with the annihilated and created momenta equal gives the forward
    operator.

    Two bookkeeping facts, both found the hard way and both cancelled by the
    normalisation rather than corrected by hand:

    * ``clfact`` returns an *un-normalised* matrix element, so it carries the
      state norm ``Norm_ii``;
    * it sums over indistinguishable partons -- including at the queried
      momentum -- so it carries the multiplicity of the parton's own
      ``(type, momentum, flavour)``.

    Dividing by the row sum removes both, which is why this needs neither the
    norm nor a multiplicity count.  The identity that they *were* the only two
    factors is what ``tests/test_endpoint`` pins.
    """
    import numpy as _np

    keys = [(int(mstate[loc, j]), int(mstate[loc + 2, j]),
             int(mstate[loc + 3, j])) for j in range(L)]
    mat = mstate[loc:loc + 4, :L].astype(int)
    matc = base.conj_state(L, mat)

    def raw(a, c):
        ta, tc = keys[a][0], keys[c][0]
        ka, kc = keys[a][1], keys[c][1]
        fa, fc = keys[a][2], keys[c][2]
        if ta == 1 and tc == 1:
            occ_b, occ_d, idx = _VERTEX["qq"]
            moms, flavs = [ka, kc, ka, kc], [fa, fc, fa, fc]
        elif ta != 1 and tc != 1:
            occ_b, occ_d, idx = _VERTEX["qbqb"]
            moms, flavs = [ka, kc, ka, kc], [fa, fc, fa, fc]
        else:
            occ_b, occ_d, idx = _VERTEX["qqb"]
            kq, kqb = (ka, kc) if ta == 1 else (kc, ka)
            fq, fqb = (fa, fc) if ta == 1 else (fc, fa)
            moms, flavs = [kqb, kqb, kq, kq], [fqb, fqb, fq, fq]
        lng = 4 + 2 * L
        mx = _np.zeros((4, lng), dtype=int)
        mx[:, :L] = mat
        mx[:, L + 4:] = matc
        mx[0, L:L + 4] = occ_b
        mx[1, L:L + 4] = occ_d
        mx[2, L:L + 4] = moms
        mx[3, L:L + 4] = flavs
        return base.clfact_compute(*idx, mx, L, L, 4, params, selfen)

    CF = (params.N * params.N - 1.0) / (2.0 * params.N)
    w = _np.zeros((L, L))
    for a in range(L):
        seen, tot = {}, 0.0
        for c in range(L):
            if c == a:
                continue
            if keys[c] not in seen:
                seen[keys[c]] = raw(a, c)
                tot += seen[keys[c]]
        if tot == 0.0:
            continue
        cnt = {}
        for c in range(L):
            if c != a:
                cnt[keys[c]] = cnt.get(keys[c], 0) + 1
        for c in range(L):
            if c != a:
                w[a, c] = CF * seen[keys[c]] / (cnt[keys[c]] * tot)
    return w


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


def state_sigmas(mstate, mstinf, numsta, N, b, K_code, table=None,
                 base=None, params=None, selfen=None):
    """``(sigma_std, sigma_imp)`` per basis state, for the matrix-level form.

    The standard self-inertia is a one-body scalar at every parton number --
    measured, ``H(selfen) - H(selfen=0) == diag(sigma_std) @ Norm`` to 1e-16 for
    L = 2..8 including untruncated bases::

        sigma_std(s) = C_F * sum_p S(k_p)

    The improvement replaces each parton's ``S(k_a)`` by its share of the
    endpoint sums with every partner::

        sigma_imp(s) = sum_a  C_F/(L-1) * sum_{c != a} J^(k_a; k_c)

    ``w_ac`` is the exact ``<-T_a . T_c>`` when ``base``/``params``/``selfen``
    are supplied (:func:`pair_colour_weights`), and falls back to the scalar
    ``C_F/(L-1)`` otherwise.

    **Use the exact weights.** The scalar is right only where every pair is
    colour-equivalent, which covers valence sectors but not Table I's
    ``LPN = valence+2``.  Measured there, weight schemes that all satisfy the
    ``b = 0`` reduction spread the answer by up to 144% -- the Fock-extended
    numbers are undetermined, not approximate, without the true operator.

    ``sigma_imp`` reduces to ``sigma_std`` exactly at ``b = 0``, per parton, not
    merely per state, because ``J^(k;l)|_{b=0} = S(k)`` independently of ``l``.
    """
    if table is None:
        table = pair_kernel_table(K_code, b)
    CF = (N * N - 1.0) / (2.0 * N)
    std = np.zeros(numsta)
    imp = np.zeros(numsta)

    # sigma is a function of the state's (type, momentum, flavour) multiset,
    # which is exactly what config_block_labels keys the norm blocks on -- so
    # compute it once per block, not once per state.  That is what keeps the
    # exact colour weights cheap: L <= 6 across all of Table I
    # (sweep_lpn = min(valence+2, 10)), and there are far fewer blocks than
    # states (119 against 189 at 2K=21).
    cache = {}
    for s in range(numsta):
        loc = int(mstinf[s, 0]) - 1
        L = int(mstinf[s, 1])
        key = tuple(sorted(
            (int(mstate[loc, j]), int(mstate[loc + 2, j]),
             int(mstate[loc + 3, j])) for j in range(L)))
        hit = cache.get(key)
        if hit is not None:
            std[s], imp[s] = hit
            continue

        ks = [int(mstate[loc + 2, j]) for j in range(L)]
        sd = CF * sum(_S(k) for k in ks)
        if L < 2:
            si = sd
        else:
            if params is None:
                w = np.full((L, L), CF / (L - 1.0))
                np.fill_diagonal(w, 0.0)
            else:
                w = pair_colour_weights(base, params, mstate, loc, L, selfen)
            si = float(sum(w[a, c] * table[ks[a], ks[c]]
                           for a in range(L) for c in range(L) if c != a))
        cache[key] = (sd, si)
        std[s], imp[s] = sd, si
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
