#!/usr/bin/env python3
"""M_meson / M_baryon against bosonization, and against just counting quarks.

Non-abelian bosonization (Date-Frishman-Sonnenschein NP B283 365) predicts a
strong-coupling ratio ``2 sin[pi nu / 2]`` with ``nu = 1/(2N-1)``
(``units.meson_baryon_ratio_bosonization``).  The 1990 paper reports agreement
to roughly 10% and takes that as support.

**It is not support, because a second law fits at least as well.**  A meson is
two quarks and a baryon is N, so naive constituent counting gives ``2/N`` -- and
that must hold trivially at large m/g, where the masses are dominated by the
quark masses rather than by any dynamics.  Where the paper has data the two are
nearly degenerate, and counting actually wins:

    N=2   measured 1.0 (forced by pseudo-reality)   boson 1.0000   2/N 1.0000
    N=3   measured 0.6657                           boson 0.6180   2/N 0.6667
    N=4   measured 0.4779                           boson 0.4450   2/N 0.5000

They separate as N grows -- 7.3% at N=3, but 13.2 / 14.6 / 15.6 / 16.4% at
N = 5,6,7,8 -- which is why SU(5)-SU(8) baryons (reachable since the NCOL_MAX,
MXTRM and MXTRMS caps were raised) are what actually decides this.

Read the output this way.  At m/g >= 0.4 the ratio *must* approach 2/N; that is
the trivial baseline and confirms the extrapolations rather than testing physics.
A genuine bosonization signal has to show up as a **departure from 2/N toward
2 sin(pi nu / 2) as m/g falls**.  A ratio that simply tracks 2/N at every
coupling means the published agreement was counting all along.

Both extrapolations use the basis appropriate to the Hamiltonian: plain 1/K for
improved (van de Sande Eq. 14), the Eq. (27) ladder for standard.  Bars are the
same fit-order x sub-window ensemble used everywhere else in this repo, carried
through the ratio.

Usage:
    python tools/meson_baryon_ratio.py
    python tools/meson_baryon_ratio.py --hamiltonian standard
    python tools/meson_baryon_ratio.py --mes-dir data/gmor_scan --bar-dir data/bar_scan
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dlcq.observables import richardson_extrapolate  # noqa: E402
from dlcq.units import meson_baryon_ratio_bosonization  # noqa: E402


def discover(directory: Path, pattern: str):
    """{N: path} for files matching e.g. 'improved_N*.csv' or 'bar_improved_N*.csv'."""
    out = {}
    if not directory.is_dir():
        return out
    for p in sorted(directory.glob(pattern)):
        m = re.search(r"_N(\d+)\.csv$", p.name)
        if m:
            out[int(m.group(1))] = p
    return dict(sorted(out.items()))


def read_series(path: Path):
    """-> {(N, mg): {K_code: msq}}, skipping header-only files."""
    out = defaultdict(dict)
    with open(path) as fh:
        for r in csv.DictReader(l for l in fh if not l.lstrip().startswith("#")):
            out[(int(r["N"]), float(r["mg"]))][int(r["K_code"])] = float(r["msq"])
    return out


def _fit_1k(ks, y, order):
    Kp = np.array(ks, float) / 2.0
    y = np.asarray(y, float)
    A = np.vstack([np.ones_like(Kp)] + [Kp ** -(i + 1.0) for i in range(order)]).T
    sc = np.linalg.norm(A, axis=0)
    c, *_ = np.linalg.lstsq(A / sc, y, rcond=None)
    return c / sc


def ensemble(series, mg, N, hamiltonian, min_pts=5):
    """(median, 68% half-width) of M^2(K->inf) over orders x contiguous windows.

    Improved output is extrapolated in a plain 1/K series and standard in the
    Eq. (27) ladder -- using either basis on the other's output is wrong in both
    directions (docs/weak-coupling-limit.md).
    """
    ks = sorted(series)
    y = [series[k] for k in ks]
    n = len(ks)
    if n < min_pts:
        return None
    vals = []
    if hamiltonian == "improved":
        for order in (2, 3, 4):
            for i in range(n):
                for j in range(i + max(min_pts, order + 2), n + 1):
                    c = _fit_1k(ks[i:j], y[i:j], order)
                    x = np.array(ks[i:j], float) / 2.0
                    A = np.vstack([np.ones_like(x)] +
                                  [x ** -(t + 1.0) for t in range(order)]).T
                    r = np.max(np.abs(A @ c - np.array(y[i:j])))
                    if r / max(np.max(np.abs(y[i:j])), 1e-300) < 1e-3 and c[0] > 0:
                        vals.append(c[0])
    else:
        for nt in (2, 3, 4):
            for i in range(n):
                for j in range(i + max(min_pts, nt + 2), n + 1):
                    try:
                        m0 = richardson_extrapolate(ks[i:j], y[i:j], mg, N,
                                                    n_terms=nt)[0]
                    except Exception:
                        continue
                    if m0 > 0:
                        vals.append(m0)
    if not vals:
        return None
    a = np.array(vals)
    return (float(np.median(a)),
            float(0.5 * (np.percentile(a, 84) - np.percentile(a, 16))))


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mes-dir", default=str(ROOT / "data" / "gmor_scan"))
    ap.add_argument("--bar-dir", default=str(ROOT / "data" / "bar_scan"))
    ap.add_argument("--hamiltonian", default="improved",
                    choices=["improved", "standard"])
    args = ap.parse_args(argv)

    ham = args.hamiltonian
    mes_files = discover(Path(args.mes_dir), f"{ham}_N*.csv")
    mes_files.update(discover(Path(args.mes_dir), f"sc_mes_N*.csv"))
    bar_files = discover(Path(args.bar_dir), f"bar_{ham}_N*.csv")
    bar_files.update(discover(Path(args.bar_dir), f"sc_bar_N*.csv"))

    mes, bar = {}, {}
    for p in mes_files.values():
        mes.update(read_series(p))
    for p in bar_files.values():
        bar.update(read_series(p))
    if not mes or not bar:
        print(f"need both channels: {len(mes)} meson and {len(bar)} baryon series"
              f"\n  meson dir {args.mes_dir}\n  baryon dir {args.bar_dir}")
        return 1

    both = sorted(set(mes) & set(bar))
    NS = sorted({N for N, _ in both})
    MGS = sorted({mg for _, mg in both}, reverse=True)
    print(f"  {ham}:  N = {NS}   m/g = {MGS}   ({len(both)} matched (N, m/g) pairs)\n")

    print("  M_meson / M_baryon    [bars = ensemble, propagated through the ratio]")
    hdr = "".join(f"{f'm/g={mg}':>19}" for mg in MGS)
    print(f"  {'N':>2}{hdr}{'2sin(pi.nu/2)':>15}{'2/N':>8}")
    table = {}
    for N in NS:
        cells = []
        for mg in MGS:
            if (N, mg) not in mes or (N, mg) not in bar:
                cells.append(f"{'--':>19}")
                continue
            em = ensemble(mes[(N, mg)], mg, N, ham)
            eb = ensemble(bar[(N, mg)], mg, N, ham)
            if not em or not eb or em[0] <= 0 or eb[0] <= 0:
                cells.append(f"{'--':>19}")
                continue
            # ratio of MASSES, so sqrt of the M^2 values; the (m^2+g^2/pi)
            # conversion factor is common to both channels and cancels exactly.
            r = np.sqrt(em[0] / eb[0])
            rel = 0.5 * np.hypot(em[1] / em[0], eb[1] / eb[0])
            table[(N, mg)] = (r, r * rel)
            cells.append(f"{r:12.4f}+-{r*rel:.4f}")
        print(f"  {N:>2}" + "".join(cells) +
              f"{meson_baryon_ratio_bosonization(N):15.4f}{2.0/N:8.4f}")

    print("\n  which law is closer, in units of the measurement's own bar:")
    print(f"  {'N':>2} {'m/g':>7} {'vs bosonization':>17} {'vs 2/N':>10}  verdict")
    for N in NS:
        for mg in MGS:
            if (N, mg) not in table:
                continue
            r, e = table[(N, mg)]
            b = meson_baryon_ratio_bosonization(N)
            c = 2.0 / N
            if abs(b - c) / max(c, 1e-12) < 0.02:
                verd = "DEGENERATE (laws agree here)"
            elif abs(r - b) < abs(r - c):
                verd = "bosonization"
            else:
                verd = "counting"
            e = max(e, 1e-12)
            print(f"  {N:>2} {mg:>7} {(r-b)/e:+17.1f}s {(r-c)/e:+10.1f}s  {verd}")
    print("\n  Expected if bosonization is right: tracks 2/N at m/g >= 0.4 and")
    print("  departs toward 2sin(pi.nu/2) as m/g falls, most visibly at N >= 5.")
    print("  Expected if it is not: tracks 2/N everywhere.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
