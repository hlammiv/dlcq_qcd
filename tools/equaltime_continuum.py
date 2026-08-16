#!/usr/bin/env python3
"""Continuum limit of the lattice Schwinger boson mass, at FIXED physical volume.

The first attempt varied ``L`` and ``x`` independently, which is wrong: physical
volume is ``Lg = L/sqrt(x)``, so at fixed ``L`` the finest lattice has the
*smallest* volume -- 60 at x=16 down to 24 at x=100.  The finite-volume error
was therefore largest exactly where the continuum extrapolation needed it
smallest, the ``A + B/L`` fits carried residuals growing with x (0.004 -> 0.019),
and two estimators of the same L -> inf limit disagreed by 9%.  The apparent
0.36% agreement with ``1/sqrt(pi)`` at x=100 was a coincidence of two
uncontrolled extrapolations.

Here ``L = round(Lg * sqrt(x))`` instead, so every lattice spacing is compared at
the same physical size.  Then the limits separate: ``Lg -> inf`` at fixed ``x``
is infinite volume at fixed spacing, and ``x -> inf`` afterwards is the continuum.

Target: ``M/g -> 1/sqrt(pi) = 0.5641896`` at ``m = 0``.

Usage:
    python tools/equaltime_continuum.py --out data/equaltime/continuum.csv
    python tools/equaltime_continuum.py --lg 30 45 --x 16 36 --chi 80
"""
from __future__ import annotations

import argparse
import csv
import logging
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
logging.getLogger("tenpy").setLevel(logging.ERROR)

from equaltime.schwinger_ed import SQRT_PI_INV          # noqa: E402
from equaltime.schwinger_fast import mass_gap           # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lg", type=float, nargs="+", default=[30.0, 45.0, 60.0],
                    help="physical volumes L*g = L/sqrt(x)")
    ap.add_argument("--x", type=float, nargs="+", default=[16.0, 36.0, 64.0, 100.0])
    ap.add_argument("--mg", type=float, default=0.0)
    ap.add_argument("--chi", type=int, default=100)
    ap.add_argument("--out", default=str(ROOT / "data" / "equaltime" / "continuum.csv"))
    args = ap.parse_args(argv)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    print(f"  target 1/sqrt(pi) = {SQRT_PI_INV:.7f}   m/g = {args.mg}")
    print(f"  {'Lg':>6} {'x':>7} {'L':>5} {'M/g':>9} {'sec':>7}", flush=True)
    for lg in args.lg:
        for x in args.x:
            L = int(round(lg * np.sqrt(x)))
            L += L % 2                      # even: mesons need even total momentum
            t = time.time()
            try:
                v = mass_gap(L, x, args.mg, chi=args.chi)
            except Exception as e:
                print(f"  {lg:>6} {x:>7} {L:>5}  {type(e).__name__}: {str(e)[:40]}",
                      flush=True)
                continue
            dt = time.time() - t
            rows.append({"Lg": lg, "x": x, "L": L, "mg": args.mg,
                         "chi": args.chi, "M_over_g": v, "seconds": round(dt, 1)})
            print(f"  {lg:>6} {x:>7} {L:>5} {v:9.5f} {dt:7.0f}", flush=True)
            with open(out, "w", newline="") as fh:
                wtr = csv.DictWriter(fh, fieldnames=list(rows[0]))
                fh.write("# Lattice Schwinger boson mass at fixed physical volume\n")
                fh.write("# Lg = L/sqrt(x) held fixed so the infinite-volume and\n")
                fh.write("# continuum limits separate.  m=0 target: M/g = 1/sqrt(pi).\n")
                fh.write("# regenerate: python tools/equaltime_continuum.py\n")
                wtr.writeheader()
                wtr.writerows(rows)
    print(f"\n  saved {out}  ({len(rows)} points)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
