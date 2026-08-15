#!/usr/bin/env python3
"""Table I's extrapolations under both Hamiltonians, on one set of axes.

The companion to ``figures.figure_fits``, which plots the standard path alone.
Every Table I number is an extrapolation, and the whole weak-coupling story is
about how much of each one is fitted rather than computed -- so the two paths
have to be seen on the same axes to be judged.

What to look for, per panel:

* the **black** points (standard) climb steeply and are still climbing at the
  left edge, so their star sits a long way past the last datum -- 15-35% of the
  answer at weak coupling;
* the **blue** points (improved) are nearly flat, so their star is a short
  extension -- <=3.5%;
* the **green** band is the published value with its own quoted last term;
* the **orange** line, on meson panels at m/g <= 0.1, is van de Sande's GMOR law
  ``M^2 = 2 pi g mu / sqrt(3)`` mapped into these units.  It is a parameter-free
  prediction of the continuum theory and neither fit sees it.

Standard is extrapolated in the Eq. (27) basis it needs; improved in a plain
1/K series, per van de Sande Eq. (14).  Using the Eq. (27) basis on improved
output, or a 1/K basis on standard output, is wrong in both directions -- see
docs/weak-coupling-limit.md.
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

from dlcq.observables import richardson_budget, richardson_extrapolate  # noqa: E402

MGS = [1.6, 0.8, 0.4, 0.2, 0.1, 0.05]
COLS = [("mes", 2, 0), ("mes", 3, 0), ("mes", 4, 0),
        ("bar", 3, 1), ("bar", 4, 1)]


def _read(path, key="msq"):
    out = defaultdict(dict)
    with open(path) as fh:
        rows = csv.DictReader(l for l in fh if not l.lstrip().startswith("#"))
        for r in rows:
            out[(int(r["N"]), int(r["B"]), float(r["mg"]))][
                int(r["K_code"])] = float(r[key])
    return out


def _paper():
    out = {}
    p = ROOT / "refs" / "table1.csv"
    with open(p) as fh:
        for r in csv.DictReader(l for l in fh if not l.startswith("#")):
            out[(r["quantity"], int(r["N"]), float(r["mg"]))] = (
                float(r["value"]), float(r["last_term"]))
    return out


def _fit_1k(ks, y, order=3):
    """Plain 1/K series -- the right basis for improved output (his Eq. 14)."""
    Kp = np.array(ks, float) / 2.0
    y = np.asarray(y)
    A = np.vstack([np.ones_like(Kp)] + [Kp ** -(i + 1.0)
                                        for i in range(order)]).T
    sc = np.linalg.norm(A, axis=0)
    c, *_ = np.linalg.lstsq(A / sc, y, rcond=None)
    return c / sc


def _spread_1k(ks, y):
    """Bar from an ensemble over orders x contiguous sub-windows."""
    Kp = np.array(ks, float) / 2.0
    y = np.asarray(y)
    n = len(ks)
    vals = []
    for order in (2, 3, 4):
        for i in range(n):
            for j in range(i + max(6, order + 2), n + 1):
                c = _fit_1k(ks[i:j], y[i:j], order)
                x, yy = Kp[i:j], y[i:j]
                A = np.vstack([np.ones_like(x)] +
                              [x ** -(t + 1.0) for t in range(order)]).T
                if np.max(np.abs(A @ c - yy)) / np.max(np.abs(yy)) < 1e-3 \
                        and c[0] > 0:
                    vals.append(c[0])
    a = np.array(vals)
    return float(np.median(a)), float(
        0.5 * (np.percentile(a, 84) - np.percentile(a, 16)))


def gmor(N, mg):
    """van de Sande Eq. (7) mapped into the repo's units. Mesons, large N."""
    c = (N * N - 1.0) / (2.0 * N * np.pi)
    return 2 * np.pi * np.sqrt(c) * mg / np.sqrt(3.0) / (mg ** 2 + 1.0 / np.pi)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--std", default=str(ROOT / "data" / "chiral_grid_msq.csv"))
    ap.add_argument("--imp", default=str(ROOT / "data" / "improved_table1_msq.csv"))
    ap.add_argument("--out", default=str(ROOT / "figures" / "table1_fits_improved_vs_standard"))
    args = ap.parse_args(argv)

    std, imp, pap = _read(args.std), _read(args.imp), _paper()

    fig, axes = plt.subplots(len(MGS), len(COLS),
                             figsize=(3.0 * len(COLS), 2.3 * len(MGS)),
                             squeeze=False)
    for ir, mg in enumerate(MGS):
        for ic, (q, N, B) in enumerate(COLS):
            ax = axes[ir][ic]
            gs, gi = std.get((N, B, mg)), imp.get((N, B, mg))
            if not gs and not gi:
                ax.text(0.5, 0.5, "no data", ha="center", va="center",
                        transform=ax.transAxes, fontsize=8, color="0.5")
                ax.set_xticks([]); ax.set_yticks([])
                continue

            hi = 0.0
            # ── standard: Eq. (27) basis, the one it needs ──
            if gs:
                ks = sorted(gs); ys = [gs[k] for k in ks]
                inv = 2.0 / np.array(ks, float)          # 1/K_paper
                bud = richardson_budget(ks, ys, mg, N, n_terms=4)
                ax.plot(inv, ys, "o", color="k", ms=3.5, zorder=4,
                        label="standard")
                ax.errorbar([0], [bud["M0"]], yerr=[bud["total"]], fmt="*",
                            color="k", ms=13, capsize=3, zorder=5)
                hi = max(hi, max(ys), bud["M0"])

            # ── improved: plain 1/K, per his Eq. (14) ──
            if gi:
                ki = sorted(gi); yi = [gi[k] for k in ki]
                invi = 2.0 / np.array(ki, float)
                m0, sd = _spread_1k(ki, yi)
                c = _fit_1k(ki, yi)
                grid = np.linspace(0, invi.max() * 1.02, 200)
                Kp = np.where(grid > 0, 2.0 / np.maximum(grid, 1e-12), np.inf)
                curve = c[0] + sum(c[i + 1] * Kp ** -(i + 1.0)
                                   for i in range(len(c) - 1))
                ax.plot(grid, curve, "-", color="C0", lw=1.0, alpha=0.8,
                        zorder=3)
                ax.plot(invi, yi, "s", color="C0", ms=3.5, zorder=4,
                        label="improved")
                ax.errorbar([0], [m0], yerr=[sd], fmt="*", color="C0", ms=13,
                            capsize=3, zorder=6)
                hi = max(hi, max(yi), m0)

            # ── published value ──
            pv, pe = pap.get((q, N, mg), (None, None))
            if pv is not None:
                ax.axhspan(pv - pe, pv + pe, color="C2", alpha=0.18, zorder=0)
                ax.axhline(pv, color="C2", lw=1.0, ls="--", zorder=2,
                           label="paper")
                hi = max(hi, pv + pe)

            # ── GMOR: parameter-free, mesons only, small m/g ──
            if B == 0 and mg <= 0.1:
                g = gmor(N, mg)
                ax.axhline(g, color="C1", lw=1.4, ls=":", zorder=2,
                           label="GMOR (Eq. 7)")
                hi = max(hi, g)

            ax.set_ylim(0, hi * 1.12)
            ax.set_xlim(left=-0.015 * max(2.0 / 25.0, 0.01))
            ax.set_title(f"{q} SU({N}),  m/g={mg}", fontsize=9)
            ax.tick_params(labelsize=7)
            if ir == len(MGS) - 1:
                ax.set_xlabel(r"$1/K$", fontsize=8)
            if ic == 0:
                ax.set_ylabel(r"$M^2/(m^2+g^2/\pi)$", fontsize=8)
            # Legend goes on a panel that actually carries every element --
            # the GMOR line only exists on meson panels at m/g <= 0.1, so the
            # top-left panel would label three of the four.
            if ir == len(MGS) - 1 and ic == 0:
                ax.legend(fontsize=6.5, loc="center right", framealpha=0.9)

    fig.suptitle("TABLE I extrapolations: improved vs standard DLCQ\n"
                 "stars are $M^2(K\\to\\infty)$; standard fitted in the Eq. (27) "
                 "basis, improved in a plain $1/K$ series (van de Sande Eq. 14)",
                 fontsize=11)
    fig.tight_layout()
    fig.subplots_adjust(top=0.93)
    for ext in ("png", "pdf"):
        p = f"{args.out}.{ext}"
        fig.savefig(p, dpi=150)
        print(f"  saved {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
