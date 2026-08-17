#!/usr/bin/env python3
"""Multi-flavour lattice Schwinger model.  The target is the two-flavour `m^{2/3}`.

Why this is the important one.  The single-flavour model validates conventions
against `M = e/sqrt(pi)`, but that is an *ordinary* mass scale.  The two-flavour
model has

    M ~ A m^{2/3} g^{1/3}

which is the same `M ~ m^{1/(2-Delta)}` logic with `Delta[psi-bar psi] = 1/2` that
gives the disputed SU(N) exponent `2N/(2N-1)`.  It is the framework's one clean
equal-time pass.  **If this pipeline cannot recover 2/3, nothing it later says
about SU(2) is interpretable** -- so this is a capability test, not a physics
result.

Construction.  Same conventions as ``schwinger_ed`` with the flavour sum added:

    Q_k = sum_f S^z_{k,f} + N_f (-1)^k / 2          (staggered charge)
    L_n = sum_{k<=n} Q_k                            (Gauss law)
    W   = x sum_{n,f} [S+_{n,f} S-_{n+1,f} + h.c.]
        + mu sum_{n,f} (-1)^n S^z_{n,f}
        + sum_n L_n^2

One chain site per *lattice* site carrying all `N_f` flavours (a TeNPy
``GroupedSite``, local dimension `2^N_f`), which keeps the hop
nearest-neighbour.  Interleaving flavours along the chain instead would make it
next-nearest and buy nothing.

The electric term keeps the running-sum structure of ``schwinger_fast``: the
weight ``N-1-max(k,l)`` depends only on the right index, so it costs one MPO bond
state whatever `N_f` is.
"""
from __future__ import annotations

import numpy as np
from tenpy.models.lattice import Chain
from tenpy.models.model import MPOModel
from tenpy.networks.mpo import MPO
from tenpy.networks.site import GroupedSite, SpinHalfSite

from .schwinger_ed import gap_to_M_over_g, mu_from


def electric_weights_nf(n_sites: int, n_flavour: int):
    """``(w, C, E0)`` for ``sum_n L_n^2`` with ``N_f`` flavours per site.

    Identical to the one-flavour case except that ``c_n`` carries ``N_f/2`` per
    site rather than ``1/2`` -- the staggered background charge is per flavour.
    """
    N = n_sites
    w = np.zeros((N, N))
    for k in range(N):
        for l in range(N):
            w[k, l] = max(N - 1 - max(k, l), 0)
    c = np.array([(n_flavour / 2.0) * (1.0 if n % 2 == 0 else 0.0)
                  for n in range(N - 1)])
    C = np.array([c[k:].sum() if k < N - 1 else 0.0 for k in range(N)])
    return w, C, float(np.sum(c ** 2))


class SchwingerNf(MPOModel):
    """``W`` for ``n_flavour`` degenerate flavours, hand-built MPO."""

    def __init__(self, L: int, x: float, mg: float, n_flavour: int = 2,
                 charges: str = "same"):
        """``charges="independent"`` conserves each flavour's number separately.

        The Hamiltonian conserves them anyway -- the hop and the mass are
        flavour diagonal, and the electric term sees only the total -- so this
        costs nothing and buys the thing that matters: the isotriplet pion then
        lives in a **different symmetry sector** from the vacuum and DMRG finds
        it directly.  With ``"same"`` only the total is conserved, the pion
        shares a sector with the vacuum, and an orthogonality-constrained sweep
        lands on the heavier flavour-singlet eta instead (measured: M/g -> 0.970
        as m -> 0, against the pion's required M/g -> 0).
        """
        nf = n_flavour
        site = GroupedSite([SpinHalfSite(conserve="Sz", sort_charge=True)
                            for _ in range(nf)], charges=charges)
        lat = Chain(L, site, bc="open", bc_MPS="finite")
        mu = mu_from(mg, x)
        w, C, E0 = electric_weights_nf(L, nf)

        # bond states: Id, {Sp_f}, {Sm_f}, Q, H
        n_hop = 2 * nf
        ID, H = 0, 1 + n_hop + 1
        SPf = [1 + f for f in range(nf)]
        SMf = [1 + nf + f for f in range(nf)]
        Q = 1 + n_hop
        dim = H + 1

        grids = []
        for n in range(L):
            g = [[None] * dim for _ in range(dim)]
            g[ID][ID] = "Id"
            g[H][H] = "Id"
            for f in range(nf):
                sp, sm, sz = f"Sp{f}", f"Sm{f}", f"Sz{f}"
                g[ID][SPf[f]] = sp
                g[ID][SMf[f]] = sm
                g[SPf[f]][H] = [(sm, x)]
                g[SMf[f]][H] = [(sp, x)]
            # running total charge: sum over flavours on this site
            g[ID][Q] = [(f"Sz{f}", 1.0) for f in range(nf)]
            g[Q][Q] = "Id"
            if n > 0 and w[0, n] != 0.0:
                g[Q][H] = [(f"Sz{f}", 2.0 * float(w[0, n])) for f in range(nf)]
            h = mu * (-1) ** n + 2.0 * float(C[n])
            if h != 0.0:
                g[ID][H] = [(f"Sz{f}", h) for f in range(nf)]
            grids.append(g)

        Hmpo = MPO.from_grids(lat.mps_sites(), grids, bc="finite",
                              IdL=ID, IdR=H, mps_unit_cell_width=L)
        self._charges = charges
        # Sz_f Sz_f = 1/4 each, and the same-site cross terms Sz_f Sz_f' (f != f')
        # both come from the diagonal of w; plus the c_n^2 constant
        self._const = E0 + 0.25 * nf * float(np.trace(w))
        self._nf = nf
        MPOModel.__init__(self, lat, Hmpo)

    @property
    def energy_offset(self) -> float:
        return self._const


def _opts(chi):
    return {"trunc_params": {"chi_max": chi, "svd_min": 1e-12},
            "max_sweeps": 60, "mixer": True, "combine": True}


def _staggered(model, L, nf):
    from tenpy.networks.mps import MPS
    # staggered vacuum: every flavour filled on odd sites, empty on even.
    # A GroupedSite labels its states by joining the constituents with a space
    # ("down_0 down_1"), so the one-flavour labels are not accepted directly.
    lo = " ".join(f"down_{f}" for f in range(nf))
    hi = " ".join(f"up_{f}" for f in range(nf))
    prod = [lo if i % 2 == 0 else hi for i in range(L)]
    return MPS.from_product_state(model.lat.mps_sites(), prod,
                                  bc=model.lat.bc_MPS, unit_cell_width=L)


def _pion_state(model, L, nf):
    """Staggered vacuum with one unit of flavour moved 0 -> 1.

    Total charge is unchanged, so the state is physical, but the individual
    flavour numbers differ from the vacuum by (+1, -1): the ``I_3 = +-1`` member
    of the isotriplet, degenerate with the neutral pion.  With
    ``charges="independent"`` this is a *different sector*, so plain DMRG in it
    returns the pion without any orthogonality constraint or penalty term.
    """
    from tenpy.networks.mps import MPS
    lo = [f"down_{f}" for f in range(nf)]
    hi = [f"up_{f}" for f in range(nf)]
    # Anchor the flips to site PARITY, not to L//2.  Even sites are empty
    # ("down") and odd sites filled ("up") in the staggered vacuum, so the flips
    # must land on an even and an odd site respectively.  Keying off L//2 works
    # only when L//2 happens to be even: at L=100 it is (site 50), at L=150 it is
    # not (site 75, already "up"), and BOTH flips silently become no-ops -- the
    # "pion" state is then literally the vacuum and the gap comes out exactly
    # 0.00000, which is what x=36 produced.
    c = 2 * (L // 4)                          # an even site near the middle
    prod = []
    for i in range(L):
        st = list(lo if i % 2 == 0 else hi)
        if i == c:                            # even site: add a flavour-0 particle
            st[0] = "up_0"
        if i == c + 1:                        # odd site: remove a flavour-1 one
            st[1] = "down_1"
        prod.append(" ".join(st))
    return MPS.from_product_state(model.lat.mps_sites(), prod,
                                  bc=model.lat.bc_MPS, unit_cell_width=L)


def pion_gap(L: int, x: float, mg: float, n_flavour: int = 2, chi: int = 120):
    """``M/g`` of the isotriplet pion, by targeting its flavour sector."""
    from tenpy.algorithms.dmrg import TwoSiteDMRGEngine
    m = SchwingerNf(L, x, mg, n_flavour, charges="independent")
    o = _opts(chi)
    E0, _ = TwoSiteDMRGEngine(_staggered(m, L, n_flavour), m, o).run()
    E1, _ = TwoSiteDMRGEngine(_pion_state(m, L, n_flavour), m, dict(o)).run()
    return gap_to_M_over_g(float(E1 - E0), x)


def mass_gap(L: int, x: float, mg: float, n_flavour: int = 2, chi: int = 120):
    """``M/g`` from the first excitation, orthogonality-constrained DMRG.

    **This finds the flavour SINGLET**, not the pion -- see ``pion_gap``.
    """
    from tenpy.algorithms.dmrg import TwoSiteDMRGEngine
    m = SchwingerNf(L, x, mg, n_flavour)
    o = _opts(chi)
    E0, psi0 = TwoSiteDMRGEngine(_staggered(m, L, n_flavour), m, o).run()
    E1, _ = TwoSiteDMRGEngine(_staggered(m, L, n_flavour), m, dict(o),
                              orthogonal_to=[psi0]).run()
    return gap_to_M_over_g(float(E1 - E0), x)


if __name__ == "__main__":
    import logging
    logging.getLogger("tenpy").setLevel(logging.ERROR)
    from .schwinger_ed import spectrum
    from .schwinger_fast import ground_state as gs1
    print("  Nf=1 through the Nf machinery must reproduce the one-flavour build\n")
    print(f"  {'L':>3} {'x':>5} {'m/g':>5} {'1-flavour':>14} {'via Nf=1':>14} {'diff':>10}")
    from tenpy.algorithms.dmrg import TwoSiteDMRGEngine
    for L, x, mg in ((8, 2.0, 0.0), (10, 4.0, 0.3)):
        a, _, _ = gs1(L, x, mg, chi=200)
        m = SchwingerNf(L, x, mg, 1)
        E, _ = TwoSiteDMRGEngine(_staggered(m, L, 1), m, _opts(200)).run()
        b = float(E) + m.energy_offset
        ed = float(spectrum(L, x, mg, k=1)[0])
        print(f"  {L:>3} {x:>5} {mg:>5} {a:14.9f} {b:14.9f} {abs(a-b):10.2e}"
              f"   (ED {ed:.9f})")
