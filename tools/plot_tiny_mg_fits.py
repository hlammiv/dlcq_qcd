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
    rows = [("improved", read(d / f"{ch}_improved_N{N}.csv"), "C0", "s"),
            ("standard", read(d / f"{ch}_standard_N{N}.csv"), "k", "o")]
    # improved's M(0) per coupling, drawn on the standard panels as the target
    # the one-sided endpoint bar is supposed to reach
    imp_ref = {}
    for mg_, g_ in rows[0][1].items():
        if len(g_) >= 5:
            r_ = extrapolate(g_, mg_, N, "improved")
            if r_:
                imp_ref[mg_] = r_[0]
    pap = paper_values()
    mgs = sorted(set().union(*[set(r[1]) for r in rows]), reverse=True)
    if not mgs:
        print(f"no data under {d}")
        return 1
    mgs = mgs[:6]

    # One row per Hamiltonian, linear y scaled to that panel.  They differ by up
    # to 1500x, so a shared axis compresses each fit into a flat line and hides
    # the only thing the figure is for -- the shape of the extrapolation.
    fig, axes = plt.subplots(len(rows), len(mgs),
                             figsize=(3.05 * len(mgs), 3.7 * len(rows)),
                             squeeze=False)
    for ir, (ham, data, col, mk) in enumerate(rows):
        for ic, mg in enumerate(mgs):
            ax = axes[ir][ic]
            g = data.get(mg)
            if not g or len(g) < 5:
                ax.text(0.5, 0.5, "no data", ha="center", va="center",
                        transform=ax.transAxes, fontsize=8, color="0.5")
                ax.set_xticks([]); ax.set_yticks([])
                continue
            ks = sorted(g)
            y = np.array([g[k] for k in ks])
            inv = 2.0 / np.array(ks, float)
            grid = np.linspace(0, inv.max() * 1.04, 300)
            res = extrapolate(g, mg, N, ham, grid)
            if res is None:
                continue
            m0, tot, form, point, curves = res
            # endpoint systematic: one-sided, upward, no free parameters --
            # M(0)/(1 - K^-2a) is where the captured-fraction model says the
            # true answer sits.  Only drawn for standard; improved already
            # subtracts the endpoint, so it does not carry this term.
            a_end = endpoint_exponent(mg, N)
            cap = 1.0 - float(max(ks)) ** (-2.0 * a_end)
            endp = m0 * (1.0 / cap - 1.0) if (ham == "standard" and cap > 0) else 0.0
            if len(curves):
                ax.plot(grid, curves[0], "-", color=col, lw=1.0, alpha=0.85, zorder=3)
            ax.errorbar(inv, y, yerr=y * POINT_REL, fmt=mk, color=col, ms=3.6,
                        elinewidth=0.8, capsize=1.5, zorder=4)
            if endp > 0:
                ax.errorbar([0], [m0], yerr=[[0.0], [endp]], fmt="none",
                            ecolor="C3", elinewidth=2.4, capsize=6, alpha=0.85,
                            zorder=5)
            ax.errorbar([0], [m0], yerr=[tot], fmt=mk, color=col, ms=8, capsize=4,
                        elinewidth=1.7, markeredgecolor="w", markeredgewidth=0.9,
                        zorder=7)
            seen = [y.min(), y.max(), m0 - tot, m0 + tot + endp]
            if endp > 0:
                gi = imp_ref.get(mg)
                if gi is not None:
                    ax.axhline(gi, color="C0", lw=1.2, ls=":", zorder=6)
                    seen.append(gi)
            pv = pap.get((ch, N, mg))
            if pv:
                ax.axhspan(pv[0] - pv[1], pv[0] + pv[1], color="C2", alpha=0.18,
                           zorder=0)
                ax.axhline(pv[0], color="C2", lw=1.1, ls="--", zorder=2)
                seen += [pv[0] - pv[1], pv[0] + pv[1]]
            lo, hi = min(seen), max(seen)
            pad = 0.10 * (hi - lo) if hi > lo else abs(hi) * 0.05 + 1e-300
            ax.set_ylim(lo - pad, hi + pad)
            ax.set_xlim(left=-0.035 * inv.max())
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 3))
            a = endpoint_exponent(mg, N)
            ttl = (f"{ham}   m/g={mg:.3g}\n$M(0)$={m0:.4g}$\\pm${100*tot/m0:.1f}%"
                   f"   grid sees {100*(1-70.0**(-2*a)):.2f}%")
            if endp > 0:
                ttl += f"\nendpoint syst. $+${endp/m0:.0f}x"
            ax.set_title(ttl, fontsize=7.5)
            ax.tick_params(labelsize=6.5)
            if endp > 0:
                ax.set_yscale("log")
                ax.set_ylim(min(seen) * 0.7, max(seen) * 1.6)
            if ir == len(rows) - 1:
                ax.set_xlabel(r"$1/K$", fontsize=8.5)
            if ic == 0:
                ax.set_ylabel(r"$M^2/(m^2+g^2/\pi)$", fontsize=8.5)

    chan = "meson" if ch == "mes" else "baryon"
    note = ("Table I has no SU(%d) entries (its colours are N=2,3,4, smallest "
            "coupling m/g=0.05), so no published value exists here." % N
            if not any(pap.get((ch, N, m)) for m in mgs)
            else "green dashed = published Table I value with its quoted last term")
    fig.suptitle(
        f"SU({N}) {chan}: K-extrapolation at small m/g.  Rows split because the two "
        f"differ by up to 1500x.\n"
        f"Improved: plain 1/K basis (van de Sande Eq. 14).  Standard: Eq. (27) "
        f"ladder, log axis, RED = one-sided endpoint systematic "
        f"$M(0)/(1-K^{{-2a}})$ — no free parameters — and blue dotted = improved "
        f"$M(0)$.\n{note}", fontsize=9.5)
    fig.tight_layout()
    fig.subplots_adjust(top=0.79, hspace=0.42)
    out = args.out or str(ROOT / "figures" / f"tiny_mg_fits_{ch}_N{N}")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{ext}".join([out + ".", ""]), dpi=150)
        print(f"  saved {out}.{ext}")

    print(f"\n  {'m/g':>10} {'ham':>9} {'M(0)':>13} {'total':>11} {'rel':>7} "
          f"{'form':>11} {'point':>11}")
    for mg in mgs:
        for ham, data, _, _ in rows:
            g = data.get(mg)
            if not g or len(g) < 5:
                continue
            r = extrapolate(g, mg, N, ham)
            if r:
                print(f"  {mg:10.3e} {ham:>9} {r[0]:13.6e} {r[1]:11.2e} "
                      f"{100*r[1]/r[0]:6.2f}% {r[2]:11.2e} {r[3]:11.2e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
