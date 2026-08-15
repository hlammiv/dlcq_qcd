#!/usr/bin/env python3
"""The lattice Schwinger model by exact diagonalisation -- the transparent build.

This is deliberately the *slow, obvious* implementation.  It exists to be the
thing the MPO is checked against, because the dominant failure mode in this
project has been convention errors, not algorithmic ones, and those are
invisible in an MPO and glaring in a 2^N matrix.

Conventions, stated once and used everywhere (Banuls et al. arXiv:1305.3765,
Dempsey-Klebanov-Pufu-Zan arXiv:2206.05308):

Staggered fermions, Jordan-Wigner to spins, gauge field eliminated by Gauss'
law with no background field.  The dimensionless Hamiltonian is

    W = (2 / (a g^2)) H
      = x * sum_n [S+_n S-_{n+1} + h.c.]        hopping
      + mu * sum_n (-1)^n S^z_n                 staggered mass
      + sum_n L_n^2                             electric

with

    x   = 1 / (a g)^2
    mu  = 2 m / (a g^2) = 2 (m/g) sqrt(x)
    L_n = sum_{k<=n} [ S^z_k + (-1)^k / 2 ]     (Gauss law)

and physical energies recovered as

    E / g = (a g / 2) W = W / (2 sqrt(x))       so   M/g = dW / (2 sqrt(x)).

**The single number that validates all of this**: at ``m = 0`` the continuum
boson mass is ``M/g = 1/sqrt(pi) = 0.5641896``.  Anything wrong in the factors
above misses it.

Expanding the electric term (``S^z S^z = 1/4`` on spin-1/2) gives a long-range
Ising coupling whose weight is ``N-1-max(k,l)`` -- linear in the site index, not
exponential, which is why the MPO for this stays small.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

SQRT_PI_INV = 1.0 / np.sqrt(np.pi)      # M/g at m=0, the continuum target


def mu_from(mg: float, x: float) -> float:
    """``mu = 2 (m/g) sqrt(x)``.  The one place the mass convention lives."""
    return 2.0 * mg * np.sqrt(x)


def gap_to_M_over_g(dW: float, x: float) -> float:
    """``M/g = dW / (2 sqrt(x))``.  The one place the energy convention lives."""
    return dW / (2.0 * np.sqrt(x))


def _ops(n_sites: int):
    """Sparse S^z, S^+ , S^- on the full 2^N space, site by site."""
    sz = sp.csr_matrix(np.array([[0.5, 0.0], [0.0, -0.5]]))
    sp_ = sp.csr_matrix(np.array([[0.0, 1.0], [0.0, 0.0]]))   # |up><down|
    I2 = sp.identity(2, format="csr")

    def embed(op, k):
        out = sp.identity(1, format="csr")
        for j in range(n_sites):
            out = sp.kron(out, op if j == k else I2, format="csr")
        return out

    return ([embed(sz, k) for k in range(n_sites)],
            [embed(sp_, k) for k in range(n_sites)])


def electric_weights(n_sites: int):
    """``(w, C)`` for ``sum_n L_n^2`` with ``n`` running over the N-1 links.

    ``sum_{n} L_n^2 = sum_{k,l} w[k,l] S^z_k S^z_l + 2 sum_k C[k] S^z_k + E0``

    with ``w[k,l] = N-1-max(k,l)`` (the number of links to the right of both
    sites), ``C[k] = sum_{n>=k} c_n``, ``c_n = (1/2) sum_{j<=n} (-1)^j`` which is
    1/2 for even ``n`` and 0 for odd, and the constant ``E0 = sum_n c_n^2``.

    **``E0`` is not optional.**  Dropping it leaves gaps correct -- constants
    cancel in differences -- while shifting the ground-state energy by exactly
    ``-N/8``, which the x=0 check catches immediately (the staggered vacuum has
    ``L_n = 0`` on every link, so ``W0`` must be exactly 0, and it came out
    -0.75, -1.00, -1.25, -1.50 at N = 6, 8, 10, 12).  It matters because the
    condensate is ``-d(epsilon)/dm`` on the ground-state energy itself, so a
    dropped constant is harmless for masses and corrupting for Phase 3.
    """
    N = n_sites
    w = np.zeros((N, N))
    for k in range(N):
        for l in range(N):
            m = max(k, l)
            w[k, l] = max(N - 1 - m, 0)
    c = np.array([0.5 if n % 2 == 0 else 0.0 for n in range(N - 1)])
    C = np.array([c[k:].sum() if k < N - 1 else 0.0 for k in range(N)])
    E0 = float(np.sum(c ** 2))
    return w, C, E0


def hamiltonian(n_sites: int, x: float, mg: float):
    """Sparse ``W`` for ``n_sites`` staggered sites, open boundaries."""
    mu = mu_from(mg, x)
    SZ, SP = _ops(n_sites)
    w, C, E0 = electric_weights(n_sites)
    dim = 2 ** n_sites
    H = sp.csr_matrix((dim, dim))

    for n in range(n_sites - 1):                       # hopping
        H = H + x * (SP[n] @ SP[n + 1].T.conj() + SP[n + 1] @ SP[n].T.conj())
    for n in range(n_sites):                           # staggered mass
        H = H + mu * ((-1) ** n) * SZ[n]
    H = H + E0 * sp.identity(dim, format="csr")        # electric constant
    for k in range(n_sites):                           # electric
        H = H + 2.0 * C[k] * SZ[k]
        for l in range(n_sites):
            H = H + w[k, l] * (SZ[k] @ SZ[l])
    return H.tocsr()


def charge_sector(n_sites: int, total_sz: float = 0.0):
    """Indices of the basis states with a given total ``S^z``.

    **Not optional.**  Gauss' law with no incoming field fixes the physical
    sector to total charge zero, ``sum_n [S^z_n + (-1)^n/2] = 0``, which for even
    ``n_sites`` is ``sum_n S^z_n = 0``.  Diagonalising the full 2^N space instead
    returns a "gap" between states in *different* charge sectors, which is not a
    mass: measured that way the m=0 gap moves AWAY from 1/sqrt(pi) as N grows
    (0.363, 0.336, 0.323 at x=4 for N=8,10,12) rather than toward it.
    """
    idx = np.arange(2 ** n_sites)
    sz = np.zeros(2 ** n_sites)
    for k in range(n_sites):
        bit = (idx >> (n_sites - 1 - k)) & 1        # bit 0 = up = +1/2
        sz += np.where(bit == 0, 0.5, -0.5)
    return np.flatnonzero(np.isclose(sz, total_sz))


def spectrum(n_sites: int, x: float, mg: float, k: int = 4,
             total_sz: float = 0.0):
    """Lowest ``k`` eigenvalues of ``W`` in the physical charge sector."""
    H = hamiltonian(n_sites, x, mg)
    sel = charge_sector(n_sites, total_sz)
    H = H[sel][:, sel]
    if H.shape[0] <= 4096:
        ev = np.linalg.eigvalsh(H.toarray())
        return ev[:k]
    return np.sort(spla.eigsh(H, k=min(k, H.shape[0] - 2), which="SA",
                              return_eigenvectors=False))


def mass_gap(n_sites: int, x: float, mg: float):
    """``M/g`` from the first excitation of ``W``.

    Physical only in the double limit ``N -> inf`` then ``x -> inf``; at finite
    ``N`` this carries the whole finite-volume error, which is why the
    production runs use VUMPS instead.  Here it is a convention check.
    """
    ev = spectrum(n_sites, x, mg, k=2)
    return gap_to_M_over_g(float(ev[1] - ev[0]), x)


if __name__ == "__main__":
    print("  lattice Schwinger model, ED.  Target at m=0: M/g -> "
          f"1/sqrt(pi) = {SQRT_PI_INV:.7f}\n")
    print(f"  {'N':>3} {'x':>6} {'m/g':>6} {'dim':>7} {'W0':>12} {'gap W':>10} {'M/g':>9}")
    for n in (8, 10, 12, 14):
        for x in (4.0, 9.0, 16.0):
            ev = spectrum(n, x, 0.0, k=2)
            print(f"  {n:>3} {x:>6} {0.0:>6} {len(charge_sector(n)):>7} "
                  f"{ev[0]:12.6f} {ev[1]-ev[0]:10.6f} "
                  f"{gap_to_M_over_g(float(ev[1]-ev[0]), x):9.5f}")
