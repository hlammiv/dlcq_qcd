#!/usr/bin/env python3
"""Reduce the Phase 4d exotics runs to one compact CSV.

For every cached run in phase4d_config.RUNS and every physical level in it,
write: label, N, B, mg, K_code, level, msq, dominant parton count, its
weight, and the 2/4/6/8-parton sector weights.  The heavy lifting (loading
446k-state style results) happens wherever this runs — Lenore, in practice
— and the paper's number extractor reads only the small CSV this emits.

Usage:
    python paper/jobs/exotics_analysis.py --out data/paper/exotics_levels.csv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "paper" / "jobs"))

from phase4d_config import RUNS  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(ROOT / "data" / "paper" /
                                         "exotics_levels.csv"))
    args = ap.parse_args(argv)

    from dlcq.figures import paper_lambda
    from dlcq.fock_weights import fock_weights
    from dlcq.observables import physical_indices
    from dlcq.providers import PythonProvider

    rows = []
    done = set()
    for label, cfg, (kind, val) in RUNS:
        key = (label, cfg["N"], cfg["B"], cfg["K_code"], val)
        if key in done:
            continue
        done.add(key)
        lam = float(paper_lambda(val)) if kind == "mg" else float(val)
        prov = PythonProvider(ncpus=4, assembly="exact", policy="blockwise",
                              solver="sparse", nev=cfg["nev"])
        r = prov.get(cfg["N"], cfg["NF"], cfg["B"], cfg["K_code"], lam,
                     LPN=cfg["LPN"])
        phys = physical_indices(r)
        for lev, idx in enumerate(phys[:cfg["nev"]]):
            try:
                w = fock_weights(r, int(idx))
            except (IndexError, ValueError):
                break
            dom = max(w, key=w.get)
            rows.append(dict(
                label=label, N=cfg["N"], B=cfg["B"], mg=val,
                K_code=cfg["K_code"], level=lev,
                msq=float(r.eigenvalues[idx]),
                dominant=dom, w_dom=round(w[dom], 6),
                **{f"w{L}": round(w.get(L, 0.0), 6)
                   for L in (2, 3, 4, 5, 6, 8, 10)}))
        del r
        print(f"{label} N={cfg['N']} B={cfg['B']} 2K={cfg['K_code']} "
              f"mg={val}: {len(phys)} levels", flush=True)

    out = Path(args.out)
    out.parent.mkdir(exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} level rows -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
