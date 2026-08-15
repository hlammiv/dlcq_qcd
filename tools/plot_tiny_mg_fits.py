#!/usr/bin/env python3
"""The K-extrapolations at very small m/g, both Hamiltonians on one set of axes.

Companion to ``tools/plot_improved_fits.py`` and follows its conventions: black
circles for standard, blue squares for improved, the filled marker at ``1/K = 0``
is ``M²(K→∞)`` with its bar, the shaded band is the same ensemble that produced
that bar, and the green line is the published value where one exists.

Why this figure exists.  At ``m/g = 1.95e-4`` the improved meson mass is 1462x
the standard one, which raises the obvious question: is standard merely slow, so
a larger K would get there, or is it converging to the wrong place?

The endpoint integral answers it.  The ``m²`` term contributes
``m² ∫ x^{2a-1} dx`` near the endpoint, whose true value is ``1/(2a)``.  A grid
whose smallest momentum fraction is ``x_min = 1/K_code`` evaluates
``(1 - K^{-2a})/(2a)`` instead — and as ``a → 0`` that **saturates at ``ln K``**
while the true value diverges:

    m/g       a         true 1/(2a)   grid sees   captured
    0.8       0.465        1.076        1.055      98.07%
    0.05      0.0315      15.86         3.727      23.50%
    0.0125    0.0079      63.41         4.109       6.48%
    1.95e-4   1.23e-4   4058            4.246       0.105%

so ``M²_true ~ m²/a ~ m·g`` (α = 1) against ``M²_grid ~ m² ln K`` (α = 2).  A
larger K does not rescue it: the remainder falls as ``K^{-a}``, so halving it at
``m/g = 1.95e-4`` needs ``2K ~ 10^2443``.

**What the picture shows, and it is not what you would guess.** Both fits look
healthy.  Standard's extrapolation adds only ~20% past its last computed point,
its curves are smooth and well determined, and their shape is identical across a
factor of 32 in ``m/g``.  Nothing inside the K-convergence reveals that the
answer is three orders of magnitude low — which is why fit-quality diagnostics
and the error budget built on them cannot detect this, and why it took an
external anchor to find.

Error bars.  Each computed point carries the intrinsic reproducibility floor of
the original algorithm, ~1e-4 relative (``docs/basis-dependence.md``: ``qcdf.f``
adds ``H₀`` to the diagonal only, so the answer is basis dependent at that
level).  The bar on ``M(0)`` propagates **both** sources in quadrature: the
fit-form ensemble (orders x contiguous sub-windows) and a Monte-Carlo over the
per-point uncertainty.  Neither alone is the whole error.

Usage:
    python tools/plot_tiny_mg_fits.py
    python tools/plot_tiny_mg_fits.py --N 6 --channel bar
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

from dlcq.observables import (richardson_curve,  # noqa: E402
                              richardson_extrapolate)
from dlcq.units import endpoint_exponent  # noqa: E402

POINT_REL = 1.0e-4          # docs/basis-dependence.md reproducibility floor
RNG = np.random.default_rng(20260815)


def read(path: Path):
    out = defaultdict(dict)
    if not path.exists():
        return out
    with open(path) as fh:
        for r in csv.DictReader(l for l in fh if not l.lstrip().startswith("#")):
            out[float(r["mg"])][int(r["K_code"])] = float(r["msq"])
    return out


def paper_values():
    """Published Table I entries, keyed (quantity, N, m/g)."""
    out = {}
    p = ROOT / "refs" / "table1.csv"
    if not p.exists():
        return out
    with open(p) as fh:
        for r in csv.DictReader(l for l in fh if not l.startswith("#")):
            out[(r["quantity"], int(r["N"]), float(r["mg"]))] = (
                float(r["value"]), float(r["last_term"]))
    return out


def fit_1k(ks, y, order=3):
    Kp = np.array(ks, float) / 2.0
    A = np.vstack([np.ones_like(Kp)] + [Kp ** -(i + 1.0) for i in range(order)]).T
    sc = np.linalg.norm(A, axis=0)
    c, *_ = np.linalg.lstsq(A / sc, np.asarray(y, float), rcond=None)
    return c / sc


def curve(c, grid):
    Kp = np.where(grid > 0, 1.0 / np.maximum(grid, 1e-300), np.inf)
    return c[0] + sum(c[i + 1] * Kp ** -(i + 1.0) for i in range(len(c) - 1))


def extrapolate(series, mg, N, ham, grid=None, draws=200):
    """``M(0)``, its total bar, and the ensemble's full-window curves.

    **Each Hamiltonian is extrapolated in the basis it needs.**  Improved
    converges in a plain ``1/K`` series (van de Sande Eq. 14); standard needs the
    Eq. (27) ladder carrying ``K^-(1+a)``.  Using ``1/K`` on standard output is
    documented as wrong in both directions, and it is not a harmless choice: at
    ``m/g = 1.95e-4`` it gives an 8.4% order spread against the ladder's 17.1%,
    so it silently halves the quoted uncertainty and biases ``M(0)`` low.

    Two error sources, combined in quadrature:
      * **form** -- spread over fit orders x contiguous sub-windows;
      * **point** -- Monte-Carlo over the per-datum reproducibility floor.
    """
    ks = sorted(series)
    y = np.array([series[k] for k in ks])
    n = len(ks)
    if n < 5:
        return None
    improved = (ham == "improved")
    orders = (2, 3, 4) if improved else (2, 3, 4, 5)

    def one(kk, yy, order):
        if improved:
            return fit_1k(kk, yy, order)[0]
        return richardson_extrapolate(kk, yy, mg, N, n_terms=order)[0]

    vals, curves = [], []
    for order in orders:
        for i in range(n):
            for j in range(i + max(5, order + 2), n + 1):
                try:
                    v = one(ks[i:j], y[i:j], order)
                except Exception:
                    continue
                if v <= 0:
                    continue
                vals.append(v)
                if grid is not None and i == 0 and j == n and order == 3:
                    if improved:
                        curves.append(curve(fit_1k(ks, y, 3), grid))
                    else:
                        try:
                            _, _, co, ex = richardson_extrapolate(
                                ks, y, mg, N, n_terms=3, return_fit=True)
                            curves.append(richardson_curve(co, ex, grid))
                        except Exception:
                            pass
    if not vals:
        return None
    a = np.array(vals)
    m0 = float(np.median(a))
    form = float(0.5 * (np.percentile(a, 84) - np.percentile(a, 16)))
    pert = []
    for _ in range(draws):
        try:
            pert.append(one(ks, y * (1.0 + RNG.normal(0, POINT_REL, n)), 3))
        except Exception:
            pass
    point = float(np.std(pert)) if pert else 0.0
    return m0, float(np.hypot(form, point)), form, point, np.array(curves)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--N", type=int, default=5)
    ap.add_argument("--channel", default="mes", choices=["mes", "bar"])
    ap.add_argument("--data-dir", default=str(ROOT / "data" / "tiny_scan"))
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)

    N, ch = args.N, args.channel
    d = Path(args.data_dir)
    imp = read(d / f"{ch}_improved_N{N}.csv")
    std = read(d / f"{ch}_standard_N{N}.csv")
    pap = paper_values()
    mgs = sorted(set(imp) | set(std), reverse=True)
    if not mgs:
        print(f"no data under {d}")
        return 1

    nsel = min(len(mgs), 6)
    ncol = 3
    nrow = int(np.ceil(nsel / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.3 * ncol, 3.9 * nrow),
                             squeeze=False)
    for ic, mg in enumerate(mgs[:nsel]):
        ax = axes[ic // ncol][ic % ncol]
        seen = []
        for data, col, mk, lbl in ((std, "k", "o", "standard"),
                                   (imp, "C0", "s", "improved")):
            g = data.get(mg)
            if not g or len(g) < 5:
                continue
            ks = sorted(g)
            y = np.array([g[k] for k in ks])
            inv = 2.0 / np.array(ks, float)
            grid = np.linspace(0, inv.max() * 1.04, 300)
            res = extrapolate(g, mg, N, lbl, grid)
            if res is None:
                continue
            m0, tot, form, point, curves = res
            if len(curves):
                lo = np.min(curves, axis=0); hi = np.max(curves, axis=0)
                ax.fill_between(grid, lo, hi, color=col, alpha=0.15, lw=0, zorder=1)
            if len(curves):
                ax.plot(grid, curves[0], "-", color=col, lw=1.0, alpha=0.85,
                        zorder=3)
            ax.errorbar(inv, y, yerr=y * POINT_REL, fmt=mk, color=col, ms=3.4,
                        elinewidth=0.8, capsize=1.5, zorder=4, label=lbl)
            ax.errorbar([0], [m0], yerr=[tot], fmt=mk, color=col, ms=8,
                        capsize=4, elinewidth=1.6, markeredgecolor="w",
                        markeredgewidth=0.9, zorder=7)
            seen += [y.min(), m0 + tot]
        q = "mes" if ch == "mes" else "bar"
        pv = pap.get((q, N, mg))
        if pv:
            ax.axhspan(pv[0] - pv[1], pv[0] + pv[1], color="C2", alpha=0.18, zorder=0)
            ax.axhline(pv[0], color="C2", lw=1.1, ls="--", zorder=2, label="paper")
            seen += [pv[0]]
        ax.set_yscale("log")
        if seen:
            ax.set_ylim(min(seen) * 0.6, max(seen) * 1.8)
        ax.set_xlim(left=-0.003)
        a = endpoint_exponent(mg, N)
        ax.set_title(f"m/g={mg:.3g}\n$a$={a:.2e}, grid sees "
                     f"{100*(1-70.0**(-2*a)):.2f}%", fontsize=8)
        ax.tick_params(labelsize=7)
        if ic // ncol == nrow - 1:
            ax.set_xlabel(r"$1/K$", fontsize=9)
        if ic % ncol == 0:
            ax.set_ylabel(r"$M^2/(m^2+g^2/\pi)$", fontsize=9)
        if ic == 0:
            ax.legend(fontsize=8, loc="center right", framealpha=0.95)
    for k in range(nsel, nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")

    chan = "meson" if ch == "mes" else "baryon"
    note = ("Table I has no SU(5) entries (its colours are N=2,3,4) and its "
            "smallest coupling is m/g=0.05, so no published value exists at any "
            "of these points."
            if not any(pap.get((ch, N, m)) for m in mgs[:ncol])
            else "green = published Table I value with its quoted last term")
    fig.suptitle(
        f"SU({N}) {chan}: K-extrapolation at small m/g — standard (black) vs "
        f"improved (blue), log scale, they differ by up to 1500x.\n"
        f"Point bars = 1e-4 reproducibility floor; the $M(0)$ bar adds the "
        f"fit-form ensemble and that floor in quadrature.\n{note}", fontsize=10)
    fig.tight_layout()
    fig.subplots_adjust(top=0.88 if nrow > 1 else 0.80)
    out = args.out or str(ROOT / "figures" / f"tiny_mg_fits_{ch}_N{N}")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=150)
        print(f"  saved {out}.{ext}")

    # numeric companion, so the bars are readable as numbers too
    print(f"\n  {'m/g':>10} {'ham':>9} {'M(0)':>13} {'total':>11} "
          f"{'form':>11} {'point':>11}")
    for mg in mgs[:ncol]:
        for data, lbl in ((std, "standard"), (imp, "improved")):
            g = data.get(mg)
            if not g or len(g) < 5:
                continue
            r = extrapolate(g, mg, N, lbl)
            if r:
                print(f"  {mg:10.3e} {lbl:>9} {r[0]:13.6e} {r[1]:11.2e} "
                      f"{r[2]:11.2e} {r[3]:11.2e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
