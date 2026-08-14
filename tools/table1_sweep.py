#!/usr/bin/env python3
"""Populate the solver cache for a Table I window, in parallel.

``figures.table1`` and ``tools/plot_holdout.py`` both read through
``PythonProvider``, so the points they need have to be *in the cache*.  Running
``figures.py --table1`` alone would compute them, but serially: ~800 solves for
2K = 25-71 across five entries and seven couplings.  This does the same solves
across processes and then leaves ``figures.py`` nothing to do but read.

``tools/large_k_sweep.py`` looks similar and is not a substitute: it calls
``run_python`` directly and writes a CSV, so it never touches the cache the
figure code reads.

Each worker builds its own provider.  Concurrent writes are safe because
``DLCQResult.save`` is atomic (tempfile + ``os.replace``) -- a race on the same
key costs a duplicate computation, not a corrupt file.

Memory is the local constraint, not cores: the largest entry (N=4 baryon at
2K=70) peaks at 1.6 GB, so ``--jobs`` x 1.6 GB has to fit.  Threads per job is
``--ncpus``; the product is the real core count in use.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dlcq.figures import MG_TABLE, _K_grid, sweep_lpn
from dlcq.units import mg_to_lambda

# (label, N, B) -- the columns Table I actually has.  No N=2 baryon.
ENTRIES = [("mes_N2", 2, 0), ("mes_N3", 3, 0), ("mes_N4", 4, 0),
           ("bar_N3", 3, 1), ("bar_N4", 4, 1)]


def _one(item):
    """Solve and cache a single (entry, K, coupling).  Returns a status tuple."""
    label, N, B, K, mg, lpn, ncpus, assembly, policy, solver = item
    os.environ.setdefault("OMP_NUM_THREADS", str(ncpus))
    from dlcq.providers import PythonProvider
    prov = PythonProvider(ncpus=ncpus, assembly=assembly, policy=policy,
                          solver=solver)
    t = time.time()
    try:
        r = prov.get(N, 1, B, K, float(mg_to_lambda(mg)), -1.0, lpn)
        return (label, K, mg, "ok", r.numsta_pre, r.n_eigenvalues,
                time.time() - t, "")
    except Exception as exc:                                  # noqa: BLE001
        # Overflow is the one that matters: the MXTRM guards raise rather than
        # truncate, so a colour-table overflow surfaces here instead of
        # silently corrupting a high-Fock sector.
        return (label, K, mg, type(exc).__name__, 0, 0, time.time() - t,
                str(exc)[:90])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--K-lo", type=int, default=25)
    ap.add_argument("--K-hi", type=int, default=71)
    ap.add_argument("--lpn-delta", type=int, default=0,
                    help="0 for the Table I sweep; 2 for the alt-LPN sweep the "
                         "error budget's truncation term needs")
    ap.add_argument("--entries", nargs="+", default=[e[0] for e in ENTRIES],
                    help=" ".join(e[0] for e in ENTRIES))
    ap.add_argument("--mg", type=float, nargs="+", default=list(MG_TABLE))
    ap.add_argument("--jobs", type=int, default=4, help="concurrent processes")
    ap.add_argument("--ncpus", type=int, default=4, help="threads per process")
    ap.add_argument("--assembly", default="exact")
    ap.add_argument("--policy", default="blockwise")
    ap.add_argument("--solver", default="sparse")
    args = ap.parse_args(argv)

    chosen = [e for e in ENTRIES if e[0] in args.entries]
    if not chosen:
        ap.error(f"no entries matched {args.entries}")

    work = []
    for label, N, B in chosen:
        lpn = sweep_lpn(N, B) + args.lpn_delta
        for K in _K_grid(B, N, args.K_lo, args.K_hi):
            for mg in args.mg:
                work.append((label, N, B, K, mg, lpn, args.ncpus,
                             args.assembly, args.policy, args.solver))

    # Biggest first: the long pole starts immediately instead of last, and a
    # memory blowup shows up in the first minute rather than the last.
    work.sort(key=lambda w: -w[3])

    print(f"{len(work)} solves | 2K = {args.K_lo}-{args.K_hi} | "
          f"LPN+{args.lpn_delta} | {args.jobs} jobs x {args.ncpus} threads | "
          f"{args.assembly}/{args.policy}/{args.solver}", flush=True)

    t0 = time.time()
    done = failed = 0
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(_one, w): w for w in work}
        for fut in as_completed(futs):
            label, K, mg, status, npre, nev, secs, msg = fut.result()
            done += 1
            if status != "ok":
                failed += 1
                print(f"  [{done}/{len(work)}] {label} 2K={K} m/g={mg} "
                      f"** {status}: {msg}", flush=True)
            elif done % 25 == 0 or secs > 30:
                print(f"  [{done}/{len(work)}] {label} 2K={K} m/g={mg} "
                      f"n={npre} nev={nev} {secs:.1f}s "
                      f"({time.time()-t0:.0f}s elapsed)", flush=True)

    print(f"done: {done - failed} ok, {failed} failed, "
          f"{time.time() - t0:.0f}s", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
