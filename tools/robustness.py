#!/usr/bin/env python3
"""Can degraded arithmetic reproduce the published higher-Fock curves?

Figs. 6(a)-(c) and thesis Fig. 18(b) disagree with our higher-Fock structure
functions.  A natural suspicion is that the 1988/1990 runs suffered numerical
damage the current code does not -- lower precision, a different weeding
threshold, or one of the array overflows qcdf.f is known to be vulnerable to
(it dimensions IDELT with 25 colour slots while an element between states of L
partons needs 2L+4, so runs reaching 12 partons are genuinely corrupted).

That suspicion is testable, and the constraints are tight: whatever moved the
five-quark distribution by ~40% at 2K=21 must have left M^2 at 1e-6, the
valence curves under 1%, and the sector total within 1.5%.

This module perturbs each stage of the solve and reports all four quantities,
so a candidate mechanism either shows up or is excluded.  As of writing every
mechanism tried is excluded -- see docs/baryon-higher-fock.md for the table.

Usage::

    python tools/robustness.py                 # 2K=21 baryon, all probes
    python tools/robustness.py --K 15 --B 1
"""

from __future__ import annotations

import argparse
import contextlib
import io
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "python"))


def build(B, K, rlamb, eps=1e-4, ncpus=8):
    """States, weeded norm, and both Hamiltonian pieces."""
    import qcdf as base

    from dlcq.read_python import weed_fortran

    p = base.Params()
    p.N, p.NF, p.B, p.K, p.rlamb = 3, 1, B, K, rlamb
    p.cutoff, p.LPN = -1.0, 0
    p.iflv[0] = 3 * B
    st, pm, fl = base.StateData(), base.PermTables(), base.FlavorTables()
    se = base.compute_selfen(p)
    with contextlib.redirect_stdout(io.StringIO()):
        base.qcdsta(p, st, pm, fl)
        _, _, hnorm = base.clrdis(0, p, st, se, ncpus=ncpus)
        n0 = st.numsta
        hn, mst, n = weed_fortran(hnorm, st.mstinf[:n0].copy(), n0, eps=eps)
        st.numsta = n
        st.mstinf[:n] = mst[:n]
        ham0, ham, _ = base.clrdis(1, p, st, se, ncpus=ncpus)
    return p, st, np.asarray(hn[:n, :n], float), mst, n, \
        np.array(ham0[0, :n, :n]), np.array(ham[:n, :n])


def solve(p, A, n, H0, H, dtype=np.float64, assembly="exact"):
    from scipy.linalg import eigh

    from dlcq.read_python import orthonormalize_blockwise

    Z, _ = orthonormalize_blockwise(A, n)
    cast = lambda M: np.asarray(M, dtype=dtype).astype(np.float64)
    hnu = cast(Z.T @ cast(H) @ Z)
    free = cast(Z.T @ cast(H0) @ Z)
    r2 = p.rlamb ** 2
    hnu = hnu * r2
    coeff = (1 - r2) * p.rmq[0] ** 2
    if assembly == "exact":
        hnu = hnu + coeff * free
    else:
        hnu[np.diag_indices(n)] += coeff * np.diag(free)
    w, z = eigh(cast(hnu))
    ev = p.K * w / 2.0
    i = int(np.argmin(np.where(ev > 1e-6, ev, np.inf)))
    return ev[i], (Z @ z)[:, i]


def sectors(st, mst, n, K, c, A, wanted):
    w = c * (A @ c)
    out = {}
    for npart in wanted:
        q = np.zeros(K + 1)
        for s in range(n):
            if int(mst[s, 1]) != npart:
                continue
            loc = int(mst[s, 0]) - 1
            for j in range(npart):
                if st.mstate[loc, j] == 1:
                    q[int(st.mstate[loc + 2, j])] += w[s]
        out[npart] = np.array([q[k] * (K / 2.0) for k in range(1, K, 2)])
    return out


def main(argv=None):
    from dlcq.figures import paper_lambda

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--B", type=int, default=1)
    ap.add_argument("--K", type=int, default=21)
    ap.add_argument("--ncpus", type=int, default=8)
    args = ap.parse_args(argv)

    lam = paper_lambda(1.6)
    p, st, A, mst, n, H0, H = build(args.B, args.K, lam, ncpus=args.ncpus)
    L = np.array([int(mst[s, 1]) for s in range(n)])
    val, hf = (3, 5) if args.B == 1 else (2, 4)

    print(f"2K={args.K} B={args.B}: {n} states, sectors "
          f"{ {int(v): int((L == v).sum()) for v in sorted(set(L))} }")
    print(f"\n{'probe':<34}{'M^2':>13}{'valence pk':>12}"
          f"{f'{hf}q total':>11}   {hf}q x1e3")

    def row(label, ev, c):
        s = sectors(st, mst, n, args.K, c, A, (val, hf))
        a = s[hf] * 1e3
        m = min(6, len(a))
        print(f"{label:<34}{ev:>13.6f}{s[val].max():>12.4f}"
              f"{a[:m].sum():>11.3f}   {np.round(a[:m], 2)}")

    ev, c = solve(p, A, n, H0, H)
    row("reference (float64, exact)", ev, c)
    for dt, nm in ((np.float32, "float32"), (np.float16, "float16")):
        ev, c = solve(p, A, n, H0, H, dtype=dt)
        row(f"arithmetic in {nm}", ev, c)
    for eps in (1e-2, 1e-6, 1e-10):
        p2, st2, A2, mst2, n2, H02, H2 = build(args.B, args.K, lam, eps=eps,
                                               ncpus=args.ncpus)
        ev, c = solve(p2, A2, n2, H02, H2)
        s = sectors(st2, mst2, n2, args.K, c, A2, (val, hf))
        a = s[hf] * 1e3
        m = min(6, len(a))
        print(f"{'weeding eps ' + format(eps, '.0e') + f' ({n2} states)':<34}"
              f"{ev:>13.6f}{s[val].max():>12.4f}{a[:m].sum():>11.3f}"
              f"   {np.round(a[:m], 2)}")
    for thresh in sorted(set(L))[1:]:
        for f, lab in ((0.0, "zeroed"), (2.0, "doubled"), (-1.0, "sign-flipped")):
            m = (L[:, None] >= thresh) | (L[None, :] >= thresh)
            Hc, H0c = H.copy(), H0.copy()
            Hc[m] *= f
            H0c[m] *= f
            ev, c = solve(p, A, n, H0c, Hc)
            row(f"L>={thresh} elements {lab}", ev, c)
    ev, c = solve(p, A, n, H0, H, assembly="fortran")
    row("ill-posed diagonal-only assembly", ev, c)
    return 0


if __name__ == "__main__":
    sys.exit(main())
