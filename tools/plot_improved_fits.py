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

from dlcq.observables import (richardson_budget, richardson_curve,  # noqa: E402
                              richardson_extrapolate)

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


def _spread_1k(ks, y, grid=None):
    """Median, 68% half-width, and the ensemble's curves on ``grid``.

    The bar and the band come from the same ensemble -- orders x contiguous
    sub-windows -- so the band is the visual form of the quoted uncertainty
    rather than a separate construction.
    """
    Kp = np.array(ks, float) / 2.0
    y = np.asarray(y)
    n = len(ks)
    vals = []
    curves = []
    Kg = None if grid is None else np.where(grid > 0,
                                            1.0 / np.maximum(grid, 1e-12), np.inf)
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
                    # Band from FULL-window fits only.  A fit on a short
                    # sub-window is meaningless outside its own range, and
                    # drawing it across the panel makes the band flare at the
                    # right edge into something that is an artefact of
                    # extrapolating a fit backwards, not an uncertainty.
                    if Kg is not None and i == 0 and j == n:
                        curves.append(c[0] + sum(
                            c[t + 1] * Kg ** -(t + 1.0) for t in range(order)))
    a = np.array(vals)
    return (float(np.median(a)),
            float(0.5 * (np.percentile(a, 84) - np.percentile(a, 16))),
            np.array(curves) if curves else None)


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

            # Track the actual span of everything drawn.  Anchoring the axis at
            # zero buries the whole panel in a sliver at the top wherever the
            # values are large -- at m/g = 1.6 the two paths differ in the
            # fourth decimal and are indistinguishable on a 0-to-5 axis.
            seen = []
            # ── standard: Eq. (27) basis, the one it needs ──
            if gs:
                ks = sorted(gs); ys = [gs[k] for k in ks]
                inv = 2.0 / np.array(ks, float)          # 1/K_paper
                bud = richardson_budget(ks, ys, mg, N, n_terms=4)
                gridS = np.linspace(0, inv.max() * 1.02, 200)
                cs = []
                for nt in (2, 3, 4, 5):
                    try:
                        _, _, co, ex = richardson_extrapolate(
                            ks, ys, mg, N, n_terms=nt, return_fit=True)
                    except Exception:
                        continue
                    cs.append(richardson_curve(co, ex, gridS))
                if cs:
                    cs = np.array(cs)
                    ax.fill_between(gridS, cs.min(axis=0), cs.max(axis=0),
                                    color="k", alpha=0.13, lw=0, zorder=1)
                # The central fit -- the one whose M(0) is the black marker.
                # Without it the standard path has a band and no curve, while
                # improved has both, which makes them hard to compare.
                try:
                    _, _, co4, ex4 = richardson_extrapolate(
                        ks, ys, mg, N, n_terms=4, return_fit=True)
                    ax.plot(gridS, richardson_curve(co4, ex4, gridS), "-",
                            color="k", lw=1.0, alpha=0.75, zorder=3)
                except Exception:
                    pass
                ax.plot(inv, ys, "o", color="k", ms=3.5, zorder=4,
                        label="standard")
                ax.errorbar([0], [bud["M0"]], yerr=[bud["total"]], fmt="o",
                            color="k", ms=6, capsize=4, elinewidth=1.6,
                            markeredgecolor="w", markeredgewidth=0.8, zorder=7)
                seen += [min(ys), max(ys), bud["M0"] - bud["total"],
                         bud["M0"] + bud["total"]]

            # ── improved: plain 1/K, per his Eq. (14) ──
            if gi:
                ki = sorted(gi); yi = [gi[k] for k in ki]
                invi = 2.0 / np.array(ki, float)
                grid = np.linspace(0, invi.max() * 1.02, 200)
                m0, sd, band = _spread_1k(ki, yi, grid)
                c = _fit_1k(ki, yi)
                # ``invi`` is 2/K_code, i.e. already 1/K_paper, and the fit is
                # in K_paper -- so this inverts to K_paper, not 2*K_paper.
                Kp = np.where(grid > 0, 1.0 / np.maximum(grid, 1e-12), np.inf)
                curve = c[0] + sum(c[i + 1] * Kp ** -(i + 1.0)
                                   for i in range(len(c) - 1))
                if band is not None:
                    ax.fill_between(grid, band.min(axis=0), band.max(axis=0),
                                    color="C0", alpha=0.20, lw=0, zorder=1)
                ax.plot(grid, curve, "-", color="C0", lw=1.0, alpha=0.9,
                        zorder=3)
                ax.plot(invi, yi, "s", color="C0", ms=3.5, zorder=4,
                        label="improved")
                ax.errorbar([0], [m0], yerr=[sd], fmt="s", color="C0", ms=6,
                            capsize=4, elinewidth=1.6, markeredgecolor="w",
                            markeredgewidth=0.8, zorder=7)
                seen += [min(yi), max(yi), m0 - sd, m0 + sd]

            # ── published value ──
            pv, pe = pap.get((q, N, mg), (None, None))
            if pv is not None:
                ax.axhspan(pv - pe, pv + pe, color="C2", alpha=0.18, zorder=0)
                ax.axhline(pv, color="C2", lw=1.0, ls="--", zorder=2,
                           label="paper")
                seen += [pv - pe, pv + pe]

            # ── GMOR: parameter-free, mesons only, small m/g ──
            if B == 0 and mg <= 0.1:
                g = gmor(N, mg)
                ax.axhline(g, color="C1", lw=1.4, ls=":", zorder=2,
                           label="GMOR (Eq. 7)")
                seen.append(g)

            if seen:
                lo, hi = min(seen), max(seen)
                pad = 0.08 * (hi - lo) if hi > lo else max(abs(hi), 1.0) * 0.05
                ax.set_ylim(lo - pad, hi + pad)
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
                ax.legend(fontsize=6.5, loc="lower right", framealpha=0.95)

    fig.suptitle("TABLE I extrapolations: improved vs standard DLCQ\n"
                 "filled markers at $1/K=0$ are $M^2(K\\to\\infty)$ with their "
                 "ensemble bar; bands are the same ensemble.  Standard fitted "
                 "in the Eq. (27) "
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
