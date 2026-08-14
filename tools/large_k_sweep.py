#!/usr/bin/env python3
"""Push the Richardson window past 2K=49, at two truncations at once.

Stages 1-3 make 2K = 51-99 affordable at the Table I truncation, which is the
range that can narrow the weak-coupling bracket: at m/g = 0.05 the increments
decay as K^-1.02, the bracket spans a factor of ~10, and the only cure is data
further out.

**Both LPN=5 and LPN=7, deliberately.**  The truncation error is not
K-independent.  Measured on the ground state at N=3, B=1:

    m/g = 1.6:  3.4e-06 -> 4.4e-06 over 2K = 25 -> 49, growth ~ K^0.38
    m/g = 0.1:  5.2e-05 -> 1.7e-04 over 2K = 25 -> 41, growth ~ K^2.4

So the strong-coupling reassurance in ``figures.sweep_lpn`` -- and every
truncation check in docs/table1-units.md, all run at 2K <= 35 -- does not
transfer to weak coupling at large K.  The error is one-signed (more Fock space
lowers the mass, as it must) and grows with the same kind of K-dependence the
Richardson fit is trying to model, so it distorts the fitted convergence rather
than merely offsetting the answer.  Running both truncations and comparing the
two extrapolations is the only honest way to use the extra reach.

Writes a CSV per (quantity, N, m/g, LPN) so the fits can be redone without
recomputing.
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import logging
logging.disable(logging.INFO)

from dlcq.figures import _K_grid, sweep_lpn
from dlcq.observables import monotone_bracket, physical_indices
from dlcq.read_python import run_python
from dlcq.units import mg_to_lambda


FIELDS = ["N", "B", "mg", "LPN", "K_code", "n_pre", "n_orth", "msq", "seconds"]


def _flush(rows, out):
    """Rewrite the CSV after every series.

    These sweeps are long and get interrupted; writing only at the end means an
    interruption costs everything computed so far.
    """
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)


def series(N, B, mg, K_lo, K_hi, lpn, nev, ncpus, out_rows, verbose=True):
    lam = float(mg_to_lambda(mg))
    Ks, ms = [], []
    for K in _K_grid(B, N, K_lo, K_hi):
        t = time.time()
        try:
            r = run_python(N=N, NF=1, B=B, K_code=K, rlamb=lam, cutoff=-1.0,
                           LPN=lpn, ncpus=ncpus, assembly="exact",
                           policy="blockwise", solver="sparse", nev=nev,
                           keep_norm=False)
        except Exception as exc:                    # OOM, overflow, whatever
            if verbose:
                print(f"    2K={K:<4} FAILED {type(exc).__name__}: "
                      f"{str(exc)[:60]}", flush=True)
            break
        phys = physical_indices(r)
        if phys.size == 0:
            continue
        m = float(r.eigenvalues[phys[0]])
        Ks.append(K)
        ms.append(m)
        out_rows.append(dict(N=N, B=B, mg=mg, LPN=lpn, K_code=K,
                             n_pre=r.numsta_pre, n_orth=r.n_orth, msq=m,
                             seconds=round(time.time() - t, 2)))
        if verbose:
            print(f"    2K={K:<4} n={r.numsta_pre:>7} -> {r.n_orth:<7} "
                  f"M2={m:.9f}   [{time.time()-t:.1f}s]", flush=True)
    return np.array(Ks, float), np.array(ms, float)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--K-lo", type=int, default=25)
    ap.add_argument("--K-hi", type=int, default=71)
    ap.add_argument("--mg", type=float, nargs="+", default=[1.6, 0.4, 0.1])
    ap.add_argument("--N", type=int, nargs="+", default=[3])
    ap.add_argument("--B", type=int, nargs="+", default=[1])
    ap.add_argument("--lpn", type=int, nargs="+", default=None,
                    help="default: sweep_lpn and sweep_lpn+2")
    ap.add_argument("--nev", type=int, default=12)
    ap.add_argument("--ncpus", type=int, default=8)
    ap.add_argument("--out", default=str(ROOT / "runs" / "large_k_sweep.csv"))
    args = ap.parse_args(argv)

    rows = []
    for N in args.N:
        for B in args.B:
            lpns = args.lpn or [sweep_lpn(N, B), sweep_lpn(N, B) + 2]
            for mg in args.mg:
                print(f"\nN={N} B={B} m/g={mg}", flush=True)
                got = {}
                for lpn in lpns:
                    print(f"  LPN={lpn}", flush=True)
                    got[lpn] = series(N, B, mg, args.K_lo, args.K_hi, lpn,
                                      args.nev, args.ncpus, rows)
                    _flush(rows, args.out)
                # The comparison this exists for.
                for lpn, (Ks, ms) in got.items():
                    if ms.size >= 4:
                        lo, hi, p, mono = monotone_bracket(Ks, ms)
                        br = (f"[{lo:.4f}, {hi:.4f}]" if mono and
                              np.isfinite(hi) else "n/a")
                        print(f"  LPN={lpn}: {ms.size} points, "
                              f"M2(K_max)={ms[-1]:.6f}, bracket {br}, p={p:.2f}",
                              flush=True)
                a, b = lpns[0], lpns[-1]
                if got[a][1].size and got[b][1].size:
                    n = min(got[a][1].size, got[b][1].size)
                    d = got[b][1][:n] - got[a][1][:n]
                    print(f"  truncation LPN{a}->LPN{b}: at 2K={int(got[a][0][n-1])} "
                          f"the difference is {d[-1]:.3e} "
                          f"({abs(d[-1])/abs(got[a][1][n-1]):.2e} relative)",
                          flush=True)

    _flush(rows, args.out)
    print(f"\nsaved {args.out}  ({len(rows)} points)")


if __name__ == "__main__":
    main()
