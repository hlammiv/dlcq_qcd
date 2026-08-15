#!/usr/bin/env python3
"""The lattice Schwinger model as a TeNPy MPO, checked against ``schwinger_ed``.

Same conventions as ``equaltime.schwinger_ed`` -- that module is the definition,
this one must reproduce it exactly.  Nothing here is trusted until
``tests/test_equaltime.py`` shows MPO and ED agreeing to machine precision on
N = 8..12, because the hopping term is the one piece the analytic ``x = 0``
checks cannot reach.

**A structural point that decides Phase 2.**  Eliminating the gauge field by
Gauss' law is what makes this model tractable on a finite open chain, but it
produces a coupling ``w[k,l] = N-1-max(k,l)`` that is *long-range and
N-dependent*.  That is fine for finite DMRG and impossible for VUMPS, which has
no N.  So the elimination trick and the infinite-volume algorithm are mutually
exclusive, and the way out is to keep the gauge field local -- which is exactly
what the loop-string-hadron basis does, and is the real reason Phase 2 uses it:
LSH is nearest-neighbour, hence VUMPS-compatible, where the eliminated basis is
not.

So this module is finite-DMRG only, by construction.  It exists to validate the
conventions and the machinery, not to produce the physics numbers.

``conserve="Sz"`` is not a performance option: it *is* the Gauss-law charge
sector.  See ``schwinger_ed.charge_sector`` for what goes wrong without it.
"""
from __future__ import annotations

import numpy as np
from tenpy.models.lattice import Chain
from tenpy.models.model import CouplingModel, MPOModel
from tenpy.networks.site import SpinHalfSite

from .schwinger_ed import electric_weights, gap_to_M_over_g, mu_from


class SchwingerChain(CouplingModel, MPOModel):
    """``W = x*hop + mu*staggered + sum_n L_n^2`` on ``L`` staggered sites.

    Parameters mirror ``schwinger_ed.hamiltonian``: ``L`` sites, ``x =
    1/(ag)^2``, ``mg = m/g``.  Energies are ``W``; convert with
    ``schwinger_ed.gap_to_M_over_g``.
    """

    def __init__(self, L: int, x: float, mg: float, bc_MPS: str = "finite"):
        site = SpinHalfSite(conserve="Sz", sort_charge=True)
        lat = Chain(L, site, bc="open", bc_MPS=bc_MPS)
        CouplingModel.__init__(self, lat)

        mu = mu_from(mg, x)
        w, C, E0 = electric_weights(L)

        # hopping -- the only term the x=0 analytic checks cannot reach, and
        # therefore the whole reason the ED comparison exists
        for i in range(L - 1):
            self.add_coupling_term(x, i, i + 1, "Sp", "Sm", plus_hc=True)

        # staggered mass, and the electric term's linear piece
        for i in range(L):
            self.add_onsite_term(mu * (-1) ** i + 2.0 * C[i], i, "Sz")

        # electric quadratic piece.  Long-range and N-dependent by construction;
        # see the module docstring for why that rules out VUMPS here.
        for k in range(L):
            for l in range(k, L):
                if w[k, l] == 0.0:
                    continue
                if k == l:
                    # Sz^2 = 1/4 on spin-1/2: a constant, folded in below
                    continue
                self.add_coupling_term(2.0 * w[k, l], k, l, "Sz", "Sz")
        self._const = E0 + 0.25 * float(np.trace(w))

        MPOModel.__init__(self, lat, self.calc_H_MPO())

    @property
    def energy_offset(self) -> float:
        """Constant part of ``W`` that the MPO does not carry.

        ``Sz_k Sz_k = 1/4`` and the ``sum_n c_n^2`` term are both pure numbers.
        They cancel in every gap, and they do **not** cancel in the ground-state
        energy -- which is what the condensate differentiates.  Add this to any
        energy before comparing with ``schwinger_ed``.
        """
        return self._const


def ground_state(L: int, x: float, mg: float, chi: int = 128,
                 bc_MPS: str = "finite", verbose: bool = False):
    """DMRG ground state.  Returns ``(W0, psi, model)`` with ``W0`` including
    :attr:`SchwingerChain.energy_offset`."""
    from tenpy.algorithms.dmrg import TwoSiteDMRGEngine
    from tenpy.networks.mps import MPS

    model = SchwingerChain(L, x, mg, bc_MPS=bc_MPS)
    # start from the staggered vacuum: the exact x=0 ground state, and the
    # physical charge sector (total Sz = 0) in one choice
    prod = ["down" if i % 2 == 0 else "up" for i in range(L)]
    psi = MPS.from_product_state(model.lat.mps_sites(), prod, bc=model.lat.bc_MPS)
    eng = TwoSiteDMRGEngine(psi, model, {
        "trunc_params": {"chi_max": chi, "svd_min": 1e-12},
        "max_sweeps": 60, "mixer": True,
        "combine": True,
    })
    E, psi = eng.run()
    return float(E) + model.energy_offset, psi, model


if __name__ == "__main__":
    from .schwinger_ed import spectrum
    print("  MPO vs ED -- the hopping term's only check\n")
    print(f"  {'L':>3} {'x':>6} {'m/g':>6} {'ED W0':>14} {'DMRG W0':>14} {'diff':>10}")
    for L in (8, 10, 12):
        for x, mg in ((2.0, 0.0), (4.0, 0.3)):
            ed = float(spectrum(L, x, mg, k=1)[0])
            dm, _, _ = ground_state(L, x, mg, chi=256)
            print(f"  {L:>3} {x:>6} {mg:>6} {ed:14.9f} {dm:14.9f} {abs(dm-ed):10.2e}")
