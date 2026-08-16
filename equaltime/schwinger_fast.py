#!/usr/bin/env python3
"""The Schwinger MPO written by hand: bond dimension 5, independent of L.

``schwinger_mps.SchwingerChain`` hands TeNPy ``L(L-1)/2`` separate ``Sz Sz``
coupling terms, so its MPO bond dimension grows with ``L`` and DMRG slows to a
crawl -- L=64 at x=16 took 46 minutes.  That is an artefact of *how* the term was
handed over, not of the physics.

The structure that fixes it.  The electric weight is ``w[k,l] = N-1-max(k,l)``,
which for ``k < l`` is ``N-1-l`` -- it depends only on the RIGHT index.  So

    sum_{k<l} 2 (N-1-l) Sz_k Sz_l  =  2 sum_l (N-1-l) Sz_l * Q_{l-1},
    Q_{l-1} = sum_{k<l} Sz_k

i.e. each site multiplies its own ``Sz`` by the running charge to its left, with
a site-dependent coefficient.  A running sum needs exactly one MPO bond state, so
the whole electric term costs one, and the total is

    [ Id | Sp | Sm | Q | H ]      -- bond dimension 5, for any L.

The MPO here is validated against ``schwinger_mps`` (which is validated against
``schwinger_ed``, which is validated analytically at x=0), so the chain of
custody runs all the way down to numbers worked out on paper.

Still finite-only: the coefficient ``N-1-l`` is N-dependent, which is the same
reason gauge elimination cannot go to VUMPS.  See ``schwinger_mps`` docstring.
"""
from __future__ import annotations

import numpy as np
from tenpy.models.lattice import Chain
from tenpy.models.model import MPOModel
from tenpy.networks.mpo import MPO
from tenpy.networks.site import SpinHalfSite

from .schwinger_ed import electric_weights, gap_to_M_over_g, mu_from

# bond-state indices
_ID, _SP, _SM, _Q, _H = 0, 1, 2, 3, 4
_DIM = 5


class SchwingerFast(MPOModel):
    """``W`` as a hand-built MPO.  Same conventions as ``schwinger_ed``."""

    def __init__(self, L: int, x: float, mg: float):
        site = SpinHalfSite(conserve="Sz", sort_charge=True)
        lat = Chain(L, site, bc="open", bc_MPS="finite")
        mu = mu_from(mg, x)
        w, C, E0 = electric_weights(L)

        grids = []
        for n in range(L):
            g = [[None] * _DIM for _ in range(_DIM)]
            g[_ID][_ID] = "Id"
            g[_H][_H] = "Id"
            # start a hop, and close the one started on the site before
            g[_ID][_SP] = "Sp"
            g[_ID][_SM] = "Sm"
            g[_SP][_H] = [("Sm", x)]
            g[_SM][_H] = [("Sp", x)]
            # running charge Q = sum of Sz to the left
            g[_ID][_Q] = "Sz"
            g[_Q][_Q] = "Id"
            # electric quadratic piece: Sz_n * Q_{n-1} with weight 2(N-1-n)
            if n > 0 and w[0, n] != 0.0:
                g[_Q][_H] = [("Sz", 2.0 * float(w[0, n]))]
            # onsite: staggered mass + electric linear piece
            h = mu * (-1) ** n + 2.0 * float(C[n])
            if h != 0.0:
                g[_ID][_H] = [("Sz", h)]
            grids.append(g)

        H = MPO.from_grids(lat.mps_sites(), grids, bc="finite",
                           IdL=_ID, IdR=_H, mps_unit_cell_width=L)
        self._const = E0 + 0.25 * float(np.trace(w))
        MPOModel.__init__(self, lat, H)

    @property
    def energy_offset(self) -> float:
        """``Sz^2 = 1/4`` terms plus ``sum_n c_n^2``: cancels in gaps, not in
        the ground-state energy the condensate differentiates."""
        return self._const


def _opts(chi: int):
    return {"trunc_params": {"chi_max": chi, "svd_min": 1e-12},
            "max_sweeps": 60, "mixer": True, "combine": True}


def _staggered(model, L):
    from tenpy.networks.mps import MPS
    prod = ["down" if i % 2 == 0 else "up" for i in range(L)]
    return MPS.from_product_state(model.lat.mps_sites(), prod,
                                  bc=model.lat.bc_MPS, unit_cell_width=L)


def ground_state(L: int, x: float, mg: float, chi: int = 128):
    """``(W0, psi, model)``, ``W0`` including :attr:`energy_offset`."""
    from tenpy.algorithms.dmrg import TwoSiteDMRGEngine
    m = SchwingerFast(L, x, mg)
    E, psi = TwoSiteDMRGEngine(_staggered(m, L), m, _opts(chi)).run()
    return float(E) + m.energy_offset, psi, m


def mass_gap(L: int, x: float, mg: float, chi: int = 128):
    """``M/g`` from the first excitation.

    ``orthogonal_to`` is a keyword-only constructor argument; in the options
    dict it is silently ignored and the gap comes out exactly 0.
    """
    from tenpy.algorithms.dmrg import TwoSiteDMRGEngine
    m = SchwingerFast(L, x, mg)
    o = _opts(chi)
    E0, psi0 = TwoSiteDMRGEngine(_staggered(m, L), m, o).run()
    E1, _ = TwoSiteDMRGEngine(_staggered(m, L), m, dict(o),
                              orthogonal_to=[psi0]).run()
    return gap_to_M_over_g(float(E1 - E0), x)


if __name__ == "__main__":
    import logging
    import time
    logging.getLogger("tenpy").setLevel(logging.ERROR)
    from .schwinger_ed import spectrum
    print("  hand-built MPO vs ED, and the speedup over the all-pairs build\n")
    print(f"  {'L':>4} {'x':>5} {'m/g':>5} {'ED W0':>14} {'fast W0':>14} "
          f"{'diff':>9} {'bond dim':>9} {'sec':>6}")
    for L, x, mg in ((8, 2.0, 0.0), (10, 4.0, 0.3), (12, 4.0, 0.0)):
        ed = float(spectrum(L, x, mg, k=1)[0])
        t = time.time()
        w0, _, m = ground_state(L, x, mg, chi=200)
        chi_mpo = max(m.H_MPO.chi)
        print(f"  {L:>4} {x:>5} {mg:>5} {ed:14.9f} {w0:14.9f} {abs(w0-ed):9.2e} "
              f"{chi_mpo:>9} {time.time()-t:6.1f}")
