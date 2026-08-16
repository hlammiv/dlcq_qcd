#!/usr/bin/env python3
"""Warm the solver cache for the Phase 4a figure panels.

Runs the pinned list of phase4a_config.py through PythonProvider, in
parallel, writing into runs/python_cache/ — the same store
paper/make_figures.py reads.  Safe to re-run: cached keys return instantly.

Usage (on Lenore, inside tmux):
    ~/venvs/dlcq/bin/python paper/jobs/phase4a_lenore.py --jobs 6
"""
from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "paper" / "jobs"))

from phase4d_config import RUNS  # noqa: E402


def one(job):
    label, cfg, (kind, val) = job
    from dlcq.figures import paper_lambda
    from dlcq.providers import PythonProvider

    lam = float(val) if kind == "lambda" else float(paper_lambda(val))
    prov = PythonProvider(ncpus=4, assembly="exact", policy="blockwise",
                          solver="sparse", nev=cfg["nev"])
    t = time.time()
    r = prov.get(cfg["N"], cfg["NF"], cfg["B"], cfg["K_code"], lam,
                 LPN=cfg["LPN"])
    return label, cfg["K_code"], lam, r.numsta_post, time.time() - t


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--jobs", type=int, default=6)
    args = ap.parse_args(argv)

    print(f"{len(RUNS)} runs, {args.jobs} jobs", flush=True)
    fails = 0
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(one, job): job for job in RUNS}
        for i, f in enumerate(as_completed(futs), 1):
            label, cfg, _ = futs[f]
            try:
                lab, K, lam, n, dt = f.result()
                print(f"[{i}/{len(RUNS)}] {lab} 2K={K} lam={lam:.4f} "
                      f"n={n} {dt:.1f}s", flush=True)
            except Exception as exc:
                fails += 1
                print(f"[{i}/{len(RUNS)}] {label} FAILED "
                      f"{type(exc).__name__}: {str(exc)[:100]}", flush=True)
    print(f"done, {fails} failures", flush=True)
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
