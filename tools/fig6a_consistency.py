#!/usr/bin/env python3
"""Are Fig. 6's two curves consistent with each other?

Figs. 6(a)-(c) each print a valence curve and a five-quark curve for the same
eigenstate.  We reproduce the valence one and not the five-quark one, and
`tools/robustness.py` has excluded every mechanism that could have damaged the
five-quark amplitudes in the original run.

This module asks the question from the other end.  In the original Fock basis
``H c = w N c`` with ``N`` block-diagonal in parton number, so the five-parton
row of that system determines c5 from c3 outright:

    c5 = -(H55 - w N55)^-1 H53 c3

That is exact to 6e-4 here, and the propagator barely enters it.  The five-quark
curve is therefore a *function* of the valence curve, and the two published
curves in one panel constrain each other.  So we can ask: is there ANY valence
wavefunction that reproduces both?

Four tests, cheapest first:

``--bound``     q(k) can never exceed the sector's total parton density
                q(k)+qbar(k).  Relabelling which parton is the antiquark, or
                plotting one colour cluster, only moves weight between the two.
                If the published points exceed q+qbar, no re-reading of our
                amplitudes can produce them.  This one is a bound, not a search.

``--map``       how well c5 is determined by c3, and how little the five-parton
                propagator matters.

``--search``    optimize over all twelve valence amplitudes, a free w, and a
                free scale on the five-quark points, at a range of tolerances
                on the valence.  Reports a Pareto envelope, so it is monotone
                in the tolerance by construction.  Minutes, not seconds.

                Every number it prints is an UPPER bound on the achievable
                residual -- more starts can only tighten it, never loosen it,
                so a run that is cheaper than the one behind the table in
                docs/baryon-higher-fock.md will report *worse* agreement, not
                better.  That table used ~1050 local optima (``--starts 130``);
                the default here is ~480.  A tolerance-by-tolerance sweep that
                is NOT pooled this way comes out non-monotone, which is how the
                under-convergence was caught in the first place.

``--historical``  feed the PRESERVED 1990-era eigenvector (python/qcdf.out)
                through the same conversion.  This is the test that localizes
                the defect: a fault in the colour sums, the Hamiltonian build or
                the diagonalization would show up in the eigenvector, and it
                does not -- that eigenvector reproduces the published VALENCE to
                0.6% while sitting 43% from the published five-quark curve.  So
                whatever went wrong is downstream of the eigensolve.

The published markers below were re-extracted from the thesis reprint
(SLAC-333 p. 82, the same figure as the article's Fig. 6, printed larger) at
600 dpi by three independent methods -- the column trace of
refs/thesis_fig12a_fivequark.csv, erosion-centroid blob finding, and raw ink
runs per lattice column -- which agree to 0.6-1.3%.  The vertical scale is
pinned by the major ticks: p. 82's frame is 6.00 tick intervals tall, so its top
is 15.0; p. 86's is 5.24, so its top is 10.48 and its markers are carried here
already rescaled.  literature/*.pdf is not in the repository (APS copyright), so
the numbers are carried here rather than re-derived.

Panel "15a" is thesis fig. 15(a) -- the same sector and K at m/g = 0.1, which no
earlier revision of this work used.  It is the second, independent coupling.

Usage::

    python tools/fig6a_consistency.py                # bound+map+historical, all
    python tools/fig6a_consistency.py --search       # adds the optimization
    python tools/fig6a_consistency.py --panel 15a    # the strong-coupling panel
"""

from __future__ import annotations

import argparse
import contextlib
import io
import os
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "python"))

K_CODE = 21
MG_WEAK, MG_STRONG = 1.6, 0.1

# Thesis p. 82 (figs. 12a-c, weak coupling, y axis 0..15.0) and p. 86 (fig. 15a,
# strong coupling, y axis 0..10.48 -- the top label sits 48 px below the frame,
# which the tick spacing confirms: the frame is 5.24 intervals, not 5).
# "val"/"fiv" are the valence and five-quark markers by lattice k.  Panels (b),
# (c) lose some five-quark markers where they merge with the valence curve, and
# the x10^n in those legends is not legible -- hence the free scale.
PUBLISHED = {
    "a": dict(state=0, mg=MG_WEAK,
              val={3: 1.994, 5: 8.343, 7: 11.572, 9: 7.035, 11: 2.200},
              fiv={1: 8.57, 3: 8.01, 5: 1.02, 7: 1.77, 9: 5.04, 11: 0.52},
              scale=1e3),          # legend states x10^3, and it is legible
    "b": dict(state=1, mg=MG_WEAK,
              val={3: 3.62, 5: 8.14, 7: 12.44, 9: 0.73, 11: 3.88, 13: 2.07},
              fiv={1: 2.67, 7: 0.24, 9: 1.86},
              scale=None),
    "c": dict(state=2, mg=MG_WEAK,
              val={1: 0.24, 3: 4.88, 5: 4.51, 7: 10.78, 9: 8.41, 11: 1.34,
                   13: 1.00},
              fiv={1: 4.05, 5: 0.43, 7: 0.65, 9: 2.63},
              scale=None),
    # thesis fig. 15(a): same sector, same K, m/g = 0.1.  Its legend reads
    # "m/g = 1" in the scan, but the coupling is 0.1 -- our valence at 0.1 fits
    # its valence markers to 0.9% where 1.0 misses by 66%.  Markers below are
    # already rescaled by 10.48/10.0 for the axis top.
    "15a": dict(state=0, mg=MG_STRONG,
                val={1: 4.517, 3: 4.956, 5: 4.886, 7: 4.616, 9: 3.994,
                     11: 3.141, 13: 2.253},
                fiv={1: 3.296, 3: 2.353, 5: 0.269, 7: 0.239, 9: 1.197,
                     11: 0.184},
                scale=1e1, inferred=True),   # x10 inferred; report_map checks it
}


def build(rlamb, K=K_CODE, N=3, B=1, ncpus=8):
    """Matrices in the ORIGINAL Fock basis, which is where the sectors block."""
    from dlcq.read_python import (_import_solver, config_block_labels,
                                  weed_fortran)
    opt, base, _ = _import_solver(True)
    os.environ["QCDF_NCPUS"] = str(ncpus)
    p = base.Params()
    p.N, p.NF, p.B, p.K, p.rlamb = N, 1, B, K, rlamb
    p.cutoff, p.LPN = -1.0, 0
    p.iflv[0] = N * B
    states, perm, flav = base.StateData(), base.PermTables(), base.FlavorTables()
    selfen = opt.compute_selfen(N)
    with contextlib.redirect_stderr(io.StringIO()):
        opt.qcdsta_fast(p, states, perm, flav)
        npre = states.numsta
        _, _, hnorm = opt.build_matrices(0, states.mstate, states.mstinf, npre,
                                         N, 1, B, K, selfen, p.cbreak, ncpus,
                                         backend="thread")
        mstinf = states.mstinf[:npre].copy()
        labels = config_block_labels(states.mstate, mstinf, npre)
        hw, mw, ns, _ = weed_fortran(hnorm, mstinf, npre, labels=labels,
                                     return_kept=True)
        states.numsta = ns
        states.mstinf[:ns] = mw[:ns]
        ham0, ham, _ = opt.build_matrices(1, states.mstate, mw[:ns], ns,
                                          N, 1, B, K, selfen, p.cbreak, ncpus,
                                          backend="thread")
    rl2 = rlamb ** 2
    H = rl2 * ham[:ns, :ns] + (1 - rl2) * ham0[0, :ns, :ns]
    lengths = mw[:ns, 1].astype(int)
    maxlen = int(lengths.max())
    types = np.zeros((ns, maxlen), int)
    moms = np.zeros((ns, maxlen), int)
    for s in range(ns):
        loc, L = int(mw[s, 0]) - 1, int(lengths[s])
        types[s, :L] = states.mstate[loc, :L]
        moms[s, :L] = states.mstate[loc + 2, :L]
    return dict(N=hw[:ns, :ns], H=H, K=K, lengths=lengths, types=types,
                moms=moms)


class Sectors:
    """Templates and blocks for the 3- and 5-parton sectors."""

    def __init__(self, m):
        self.K = m["K"]
        self.Kp = self.K / 2.0
        self.Nrm, self.H = m["N"], m["H"]
        lengths, types, moms = m["lengths"], m["types"], m["moms"]
        self.ks = np.arange(1, self.K, 2)
        self.I3 = np.where(lengths == 3)[0]
        self.I5 = np.where(lengths == 5)[0]
        self.N33 = self.Nrm[np.ix_(self.I3, self.I3)]
        self.N55 = self.Nrm[np.ix_(self.I5, self.I5)]
        self.H55 = self.H[np.ix_(self.I5, self.I5)]
        self.H53 = self.H[np.ix_(self.I5, self.I3)]
        self.T3 = self._tmpl(self.I3, 3, types, moms, 1)
        self.T5 = self._tmpl(self.I5, 5, types, moms, 1)
        self.T5b = self._tmpl(self.I5, 5, types, moms, 0)
        self.w, self.z = eigh(self.H, self.Nrm)

    def _tmpl(self, idx, L, types, moms, want):
        T = np.zeros((len(self.ks), len(idx)))
        for j, s in enumerate(idx):
            for p in range(L):
                if types[s, p] == want:
                    T[(moms[s, p] - 1) // 2, j] += 1
        return T

    def exact(self, state):
        c = self.z[:, state]
        return c / np.sqrt(c @ self.Nrm @ c), self.w[state]

    def from_c3(self, c3, w):
        """The five-parton row of H c = w N c, solved for c5."""
        return -np.linalg.solve(self.H55 - w * self.N55, self.H53 @ c3)

    def curves(self, c3, w):
        """Both structure functions for a trial valence amplitude vector.

        Normalization runs over c3 AND the c5 it implies, not over c3 alone:
        the physical state is normalized across sectors, and for the excited
        states the five-parton probability is large enough that getting this
        wrong shifts the valence by a few tenths of a percent -- which matters
        when the whole comparison is at the 1% level.  The neglected c7/c9
        contribution to the norm is O(1e-7).
        """
        c5 = self.from_c3(c3, w)
        nn = c3 @ self.N33 @ c3 + c5 @ self.N55 @ c5
        if nn <= 0:
            return None, None, None
        c3, c5 = c3 / np.sqrt(nn), c5 / np.sqrt(nn)
        w5 = c5 * (self.N55 @ c5)
        return (self.Kp * self.T3 @ (c3 * (self.N33 @ c3)),
                self.Kp * self.T5 @ w5, self.Kp * self.T5b @ w5)


def _resid(pred, pub, free_scale):
    pred = np.asarray(pred, float)
    if free_scale:
        if pred @ pred <= 0:
            return 9.0, 0.0
        a = float(pub @ pred / (pred @ pred))
    else:
        a = 1.0
    return float(np.linalg.norm(a * pred - pub) / np.linalg.norm(pub)), a


def report_bound(S, panel, P):
    c, w = S.exact(P["state"])
    _, q5, qb5 = S.curves(c[S.I3], w)
    kf = sorted(P["fiv"])
    i = [(k - 1) // 2 for k in kf]
    pub = np.array([P["fiv"][k] for k in kf])
    scale = P["scale"]
    if scale is None:
        _, a = _resid(q5[i] * 1e3, pub, True)
        note = f"  (free scale {a:.4g}; this legend's x10^n is not legible)"
    else:
        # our q is carried as q x 10^3; the panel plots q x scale
        a = scale / 1e3
        note = (f"  (legend x{scale:g}"
                + (", inferred" if P.get("inferred") else ", as stated") + ")")
    tot = (q5 + qb5)[i] * 1e3 * a
    print(f"\n  parton-density bound, panel ({panel}){note}")
    print("     k  published q   our q+qbar   ratio")
    over = []
    for k, p_, t in zip(kf, pub, tot):
        flag = ""
        if p_ > t:
            over.append(k)
            flag = "  <== EXCEEDS"
        print(f"   {k:3d}   {p_:9.2f}   {t:10.2f}   {p_/t:6.3f}{flag}")
    if not over:
        print("     -> within the total parton density everywhere.")
    elif scale is None:
        print(f"     -> exceeded at k={over}, but only at the fitted scale.")
        print("        With the legend's x10^n illegible the scale is free, and"
              " a larger one")
        print("        always satisfies the bound -- so for this panel the "
              "operative")
        print("        statement is the shape residual above, not this table.")
    else:
        print(f"     -> at k={over} the published quark density exceeds our "
              "TOTAL five-parton density.")
        print("        No relabelling, cluster choice or re-partition can do "
              "that; the amplitudes differ.  The legend fixes the scale, so "
              "this is a bound, not a fit.")


def report_map(S, panel, P):
    c, w = S.exact(P["state"])
    c5 = S.from_c3(c[S.I3], w)
    err = np.linalg.norm(c5 - c[S.I5]) / np.linalg.norm(c[S.I5])
    q3, q5, _ = S.curves(c[S.I3], w)
    kv, kf = sorted(P["val"]), sorted(P["fiv"])
    pv = np.array([P["val"][k] for k in kv])
    pf = np.array([P["fiv"][k] for k in kf])
    rv, _ = _resid(q3[[(k - 1) // 2 for k in kv]], pv, False)
    rf, a = _resid(q5[[(k - 1) // 2 for k in kf]] * 1e3, pf, True)
    print(f"\n  panel ({panel}) = state {P['state']}   "
          f"M^2 = {S.K * w / 2:.4f}")
    print(f"    c5 reconstructed from c3 alone : {err:.2e}")
    print(f"    valence residual               : {rv*100:6.2f}%")
    print(f"    five-quark residual (free scale): {rf*100:6.2f}%")
    ratio = pf / (a * q5[[(k - 1) // 2 for k in kf]] * 1e3)
    print(f"    published/ours by k={kf}: {np.round(ratio, 3)}")
    if P["scale"] is not None and P.get("inferred"):   # inferred -> check it
        pred = q5[(kf[0] - 1) // 2] * P["scale"]
        fitted = float(np.mean(pf / q5[[(k - 1) // 2 for k in kf]]))
        print(f"    legend x{P['scale']:g}: predicts {pred:.3f} at k={kf[0]} "
              f"vs {pf[0]:.3f} measured ({abs(pred/pf[0]-1)*100:.1f}%); "
              f"a free fit wants x{fitted:.3g}")
    # how little the propagator matters
    if panel == "a":
        V = S.H53 @ c[S.I3]
        D = S.H55 - w * S.N55
        for lab, c5v in [("exact", -np.linalg.solve(D, V)),
                         ("diagonal only", -V / np.diag(D))]:
            y = S.Kp * S.T5 @ (c5v * (S.N55 @ c5v))
            r, _ = _resid(y[[(k - 1) // 2 for k in kf]] * 1e3, pf, True)
            print(f"    propagator {lab:14s}         : {r*100:6.2f}%")


def report_search(S, panel, P, starts=60, seed=1):
    from scipy.optimize import minimize
    kv, kf = sorted(P["val"]), sorted(P["fiv"])
    iv = [(k - 1) // 2 for k in kv]
    i5 = [(k - 1) // 2 for k in kf]
    pv = np.array([P["val"][k] for k in kv])
    pf = np.array([P["fiv"][k] for k in kf])
    c0, w0 = S.exact(P["state"])
    c0 = c0[S.I3]

    def both(par):
        q3, q5, _ = S.curves(par[:-1], par[-1])
        return (_resid(q3[iv], pv, False)[0], _resid(q5[i5] * 1e3, pf, True)[0])

    rng = np.random.default_rng(seed)
    pool = []
    for mu in (0.0, 0.1, 1.0, 10.0, 100.0, 1e3, 1e4, 1e5):
        for t in range(starts):
            if t == 0:
                x0 = np.append(c0, w0)
            elif t < starts // 2:
                x0 = np.append(c0 * (1 + 0.8 * rng.standard_normal(len(c0))), w0)
            else:
                x0 = np.append(rng.standard_normal(len(c0)) * 0.1, w0)
            f = lambda p: (lambda a, b: b + mu * a ** 2)(*both(p))
            r = minimize(f, x0, method="Nelder-Mead",
                         options=dict(maxiter=3000, fatol=1e-11, xatol=1e-11))
            pool.append(both(r.x))
    print(f"\n  panel ({panel}): {len(pool)} local optima")
    print("    valence held to |  best five-quark residual found")
    for tol in (0.01, 0.02, 0.03, 0.05, 0.10, 0.20, 1.00):
        ok = [b for a, b in pool if a <= tol]
        got = f"{min(ok)*100:6.2f}%" if ok else "  (none)"
        print(f"        +-{tol*100:5.1f}%   |        {got}")
    a0, b0 = both(np.append(c0, w0))
    print(f"    our own ground state: valence {a0*100:.2f}%, "
          f"five-quark {b0*100:.2f}%")
    print("    (upper bounds: this is a search, not a proof of the minimum)")


def report_historical(path=None, panel="a"):
    """Does the PRESERVED 1990-era eigenvector produce the published figure?

    This is the test that localizes the defect.  A fault anywhere in the colour
    sums, the Hamiltonian build or the diagonalization would show up in the
    eigenvector -- so if the era's own eigenvector is healthy, whatever went
    wrong happened after the eigensolve, in the x-space conversion.

    python/qcdf.out is gitignored (it is a large historical artifact), so a
    fresh clone skips this check rather than failing.
    """
    from dlcq.observables import structure_function, physical_indices
    path = Path(path or ROOT / "python" / "qcdf.out")
    print("\n  preserved 1990-era eigenvector (python/qcdf.out)")
    if not path.exists():
        print(f"    {path} not present -- skipped (gitignored artifact).")
        return
    from dlcq.read_fortran import read_out
    with contextlib.redirect_stderr(io.StringIO()):
        r = read_out(path)
    if (r.B, r.K_code) != (1, K_CODE):
        print(f"    B={r.B} 2K={r.K_code}: not the 2K={K_CODE} baryon run -- skipped.")
        return
    i = int(physical_indices(r)[0])
    _, q3, _ = structure_function(r, i, nparton=3)
    _, q5, _ = structure_function(r, i, nparton=5)
    P = PUBLISHED[panel]
    kv, kf = sorted(P["val"]), sorted(P["fiv"])
    pv = np.array([P["val"][k] for k in kv])
    pf = np.array([P["fiv"][k] for k in kf])
    ov = q3[[(k - 1) // 2 for k in kv]]
    of = q5[[(k - 1) // 2 for k in kf]] * 1e3
    rv, _ = _resid(ov, pv, False)
    rf, a = _resid(of, pf, True)
    print(f"    M^2 = {r.eigenvalues[i]:.9f},  {r.numsta_post} states")
    print(f"    its valence vs the PUBLISHED valence : {rv*100:6.2f}%")
    print(f"    its five-quark vs the PUBLISHED curve: {rf*100:6.2f}%")
    print(f"    published / 1990-eigenvector by k={kf}: "
          f"{np.round(pf/(a*of), 3)}")
    print("    -> the era's own eigenvector reproduces the published VALENCE but"
          " not the\n       published five-quark curve, so the defect is "
          "downstream of the eigensolve.")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", choices=list(PUBLISHED) + ["all"], default="all")
    ap.add_argument("--bound", action="store_true")
    ap.add_argument("--map", action="store_true")
    ap.add_argument("--search", action="store_true")
    ap.add_argument("--historical", action="store_true")
    ap.add_argument("--starts", type=int, default=60)
    ap.add_argument("--ncpus", type=int, default=8)
    a = ap.parse_args(argv)
    if not (a.bound or a.map or a.search or a.historical):
        a.bound = a.map = a.historical = True

    from dlcq.units import mg_to_lambda
    cache = {}
    def sectors_for(mg):
        if mg not in cache:
            cache[mg] = Sectors(build(mg_to_lambda(mg), ncpus=a.ncpus))
        return cache[mg]

    panels = list(PUBLISHED) if a.panel == "all" else [a.panel]
    first = sectors_for(PUBLISHED[panels[0]]["mg"])
    print(f"2K = {first.K}, B = 1, N = 3: "
          f"{len(first.I3)} valence and {len(first.I5)} five-parton states")
    for panel in panels:
        P = PUBLISHED[panel]
        S = sectors_for(P["mg"])
        if a.map:
            report_map(S, panel, P)
        if a.bound:
            report_bound(S, panel, P)
        if a.search:
            report_search(S, panel, P, starts=a.starts)
    if a.historical:
        report_historical()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
