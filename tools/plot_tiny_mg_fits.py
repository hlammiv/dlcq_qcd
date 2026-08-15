#!/usr/bin/env python3
"""The K-extrapolations at very small m/g, standard against improved.

Why this figure exists.  At `m/g = 1.95e-4` the improved meson mass is 1462x the
standard one, and the natural question is whether standard is merely slow -- a
larger K would get there -- or converging to the wrong place.

The endpoint exponent answers it.  Increments fall as `K^-(1+a)`, so the
*remainder* falls as `K^-a`, and `a -> 0` in the chiral limit:

    m/g        a          captured at 2K=70    2K to halve the remainder
    0.8        0.465          98.1%                      4.4
    0.05       0.0315         23.5%                      3.5e9
    0.0125     0.0079          6.5%                      1.5e38
    1.95e-4    1.2e-4          0.105%                    1e2443

So standard is not wrong in principle -- its K -> infinity limit is the right
one -- but at these couplings it has captured a tenth of a percent of the answer
and the rest is unreachable by any conceivable K.  `1/captured = 956` accounts
for the measured 1462x to within 1.5x.

What to look for in each panel:

* **improved (blue)** is nearly flat in `1/K`: the extrapolation adds 0.3-0.6%
  to the last computed meson point, ~6% for the baryon.
* **standard (black)** adds ~20-22%, with visible curvature -- a bigger
  extension, but by no means a broken-looking one.
* **and that is the point worth taking away.** The standard fits look
  *healthy*: smooth, well determined, stable in shape across a factor of 8 in
  `m/g` (the last/M(0) percentages are identical column to column). Nothing in
  the K-convergence reveals that the answer is ~1500x too small. The failure is
  invisible from inside the extrapolation, which is precisely why it survived
  1990-2026 and why the error budget built on these fits cannot detect it.

Usage:
    python tools/plot_tiny_mg_fits.py
    python tools/plot_tiny_mg_fits.py --N 6 --out figures/tiny_fits_N6
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from dlcq.units import endpoint_exponent  # noqa: E402


def read(path: Path):
    out = defaultdict(dict)
    if not path.exists():
        return out
    with open(path) as fh:
        for r in csv.DictReader(l for l in fh if not l.lstrip().startswith("#")):
            out[float(r["mg"])][int(r["K_code"])] = float(r["msq"])
    return out


def fit_1k(ks, y, order=3):
    Kp = np.array(ks, float) / 2.0
    y = np.asarray(y, float)
    A = np.vstack([np.ones_like(Kp)] + [Kp ** -(i + 1.0) for i in range(order)]).T
    sc = np.linalg.norm(A, axis=0)
    c, *_ = np.linalg.lstsq(A / sc, y, rcond=None)
    return c / sc


def curve(c, grid):
    Kp = np.where(grid > 0, 1.0 / np.maximum(grid, 1e-300), np.inf)
    return c[0] + sum(c[i + 1] * Kp ** -(i + 1.0) for i in range(len(c) - 1))


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--N", type=int, default=5)
    ap.add_argument("--data-dir", default=str(ROOT / "data" / "tiny_scan"))
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)

    N = args.N
    d = Path(args.data_dir)
    rows = [("meson, improved", read(d / f"mes_improved_N{N}.csv"), "C0", "s"),
            ("baryon, improved", read(d / f"bar_improved_N{N}.csv"), "C0", "s"),
            ("meson, standard", read(d / f"mes_standard_N{N}.csv"), "k", "o"),
            ("baryon, standard", read(d / f"bar_standard_N{N}.csv"), "k", "o")]
    mgs = sorted(set().union(*[set(r[1]) for r in rows]), reverse=True)
    if not mgs:
        print(f"no data under {d}")
        return 1
    if len(mgs) > 4:                      # span the range, not just the top
        idx = [0, len(mgs)//3, 2*len(mgs)//3, len(mgs)-1]
        mgs = [mgs[i] for i in sorted(set(idx))]

    fig, axes = plt.subplots(len(rows), len(mgs),
                             figsize=(3.1 * len(mgs), 2.5 * len(rows)),
                             squeeze=False)
    for ir, (label, data, col, mk) in enumerate(rows):
        for ic, mg in enumerate(mgs):
            ax = axes[ir][ic]
            g = data.get(mg)
            if not g or len(g) < 5:
                ax.text(0.5, 0.5, "no data", ha="center", va="center",
                        transform=ax.transAxes, fontsize=8, color="0.5")
                ax.set_xticks([]); ax.set_yticks([])
                continue
            ks = sorted(g)
            y = [g[k] for k in ks]
            inv = 2.0 / np.array(ks, float)          # 1/K_paper
            c = fit_1k(ks, y)
            grid = np.linspace(0, inv.max() * 1.03, 300)
            ax.plot(grid, curve(c, grid), "-", color=col, lw=1.0, alpha=0.8, zorder=2)
            ax.plot(inv, y, mk, color=col, ms=3.5, zorder=4)
            ax.plot([0], [c[0]], mk, color=col, ms=8, markeredgecolor="w",
                    markeredgewidth=1.0, zorder=6)
            # how much of the extrapolated answer the computed points actually reach
            frac = max(y) / c[0] if c[0] > 0 else np.nan
            ax.set_title(f"{label}\nm/g={mg:.3g}   last/M(0)={100*frac:.1f}%",
                         fontsize=7.5)
            ax.tick_params(labelsize=6.5)
            ax.set_xlim(left=-0.02 * inv.max())
            if ir == len(rows) - 1:
                ax.set_xlabel(r"$1/K$", fontsize=8)
            if ic == 0:
                ax.set_ylabel(r"$M^2/(m^2+g^2/\pi)$", fontsize=8)

    a_small = endpoint_exponent(min(mgs), N)
    fig.suptitle(
        f"K-extrapolations at small m/g, SU({N}).  Both look like healthy convergent "
        f"series — improved (blue) adds 0.3–6% past the last point, standard (black) "
        f"~20–22%,\nand the shape is stable across a factor of 8 in m/g.  Nothing here "
        f"reveals that the standard answer is ~1500x too small.\n"
        f"At m/g={min(mgs):.3g}: a={a_small:.2e}, the grid captures "
        f"{100*(1-70.0**(-2*a_small)):.3f}% of the endpoint weight, and halving the "
        f"remainder needs 2K ~ 1e{np.log10(2.0)/a_small:.0f}.", fontsize=9)
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)
    out = args.out or str(ROOT / "figures" / f"tiny_mg_fits_N{N}")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=150)
        print(f"  saved {out}.{ext}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
