#!/usr/bin/env python3
"""Plot each Richardson fit together with a held-out-K test of it.

``figures.figure_fits`` already shows the fit and where it lands.  What it
cannot show is whether the fit *deserves* to be believed, because every point on
screen was used to make it.  This adds the missing check: refit on the low-K
points only, then predict the two highest K -- which the fit never saw -- and
draw the prediction against the truth.

Two numbers are annotated per panel.

``extrap``  how much of the answer comes from the fit rather than the data:
            ``(M0 - M2(K_max)) / M0``.  At strong coupling this is ~1%; at weak
            coupling it is ~25%, which is the real content of "the extrapolation
            is doing the work here".

``holdout`` the relative error of the withheld prediction.  It is basis
            independent, unlike the paper's last-term rule, which moves by a
            factor of several under a change of coordinates that leaves M0
            identical to 14 digits (see observables._richardson_design).

Read the two together.  A small holdout error with a large extrap fraction --
the weak-coupling panels -- means the fit interpolates its own neighbourhood
well while still extrapolating a quarter of the answer on trust.  Holdout is a
necessary check, not a sufficient one.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from dlcq.figures import _K_grid, _mass_series, _load_paper_table1
from dlcq.observables import (richardson_evaluate, richardson_extrapolate,
                              richardson_holdout)
from dlcq.providers import PythonProvider

# Categorical slots 1 and 2 of the reference palette.  Documented as passing
# the all-pairs gates in both modes (CVD dE 9.2 light, normal-vision 24.0), and
# neither sits in the sub-3:1 contrast set.  Identity is carried by marker shape
# as well as hue, so the panels survive greyscale and CVD.
C_FIT = "#2a78d6"     # the full-window fit: what Table I actually reports
C_HOLD = "#eb6834"    # the reduced fit and the points withheld from it
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#d8d7d2"

N_HOLD = 2
N_TERMS = 2           # figures.py:334 -- what Table I uses


def panel(ax, prov, q, N, B, mg, paper, lo=25, hi=49):
    Ks, ms = _mass_series(prov, N, B, mg, _K_grid(B, N, lo, hi), msq_units=True)
    if len(Ks) < N_HOLD + 3:
        ax.set_visible(False)
        return None
    Ks = np.asarray(Ks, float)
    ms = np.asarray(ms, float)
    x = 2.0 / Ks                                   # 1 / K_paper

    M0, _ = richardson_extrapolate(Ks, ms, mg, N, n_terms=N_TERMS)
    _, rel, pred, actual, K_held = richardson_holdout(
        Ks, ms, mg, N, n_terms=N_TERMS, n_hold=N_HOLD)

    # Full-window fit, continued to the axis.
    _, _, coeffs, _ = richardson_extrapolate(
        Ks, ms, mg, N, n_terms=N_TERMS, return_fit=True)
    grid_K = np.linspace(Ks.min(), 4000.0, 600)
    ax.plot(2.0 / grid_K, richardson_evaluate(coeffs, mg, N, grid_K,
                                              n_terms=N_TERMS),
            "-", color=C_FIT, lw=2, zorder=2)

    # Fit on the low-K points only, extended over the gap it must predict.
    K_fit = Ks[:-N_HOLD]
    _, _, c_hold, _ = richardson_extrapolate(
        K_fit, ms[:-N_HOLD], mg, N, n_terms=N_TERMS, return_fit=True)
    gap = np.linspace(K_fit.min(), Ks.max(), 200)
    ax.plot(2.0 / gap, richardson_evaluate(c_hold, mg, N, gap,
                                           n_terms=N_TERMS),
            "--", color=C_HOLD, lw=2, zorder=3)

    if paper is not None:
        ax.axhline(paper, color=INK_2, lw=1, ls=":", zorder=1)

    ax.plot(x[:-N_HOLD], ms[:-N_HOLD], "o", color=INK, ms=5,
            mfc="white", mew=1.5, zorder=4)
    ax.plot(x[-N_HOLD:], actual, "s", color=C_HOLD, ms=8,
            mfc="white", mew=2, zorder=6)
    ax.plot(2.0 / K_held, pred, "x", color=C_HOLD, ms=9, mew=2, zorder=5)
    ax.plot([0], [M0], "*", color=C_FIT, ms=15, zorder=7)

    extrap = (M0 - ms[-1]) / abs(M0) if M0 else float("nan")
    ax.set_title(f"{q} SU({N})   m/g = {mg:g}", fontsize=9, color=INK, pad=4)
    # Upper right: the data runs lower-right to upper-left, so this is the one
    # corner no panel occupies.  Boxed, because the published-value line passes
    # near the top in the strong-coupling panels.
    ax.text(0.97, 0.94,
            f"extrap {100*extrap:.1f}%\nheld-out {rel:.1e}",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
            color=INK_2, linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="none",
                      alpha=0.85))

    ax.set_xlim(-0.004, x.max() * 1.06)
    ax.grid(True, color=GRID, lw=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(labelsize=7.5, colors=INK_2, length=3)
    return M0


def main():
    prov = PythonProvider(ncpus=8, assembly="exact", policy="fortran")
    paper = _load_paper_table1()
    entries = [("mes", 2, 0, 1.6), ("mes", 3, 0, 1.6), ("bar", 3, 1, 1.6),
               ("mes", 2, 0, 0.1), ("mes", 3, 0, 0.1), ("bar", 3, 1, 0.1)]

    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.2))
    for ax, (q, N, B, mg) in zip(axes.ravel(), entries):
        ref = paper.get((q, N, mg))
        panel(ax, prov, q, N, B, mg, ref[0] if ref else None)

    for ax in axes[-1]:
        ax.set_xlabel(r"$1/K_{\rm paper}$", fontsize=8.5, color=INK_2)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$M^2/(m^2+g^2/\pi)$", fontsize=8.5, color=INK_2)

    handles = [
        plt.Line2D([], [], marker="o", ls="none", mfc="white", mec=INK,
                   mew=1.5, ms=5, label="used in fit"),
        plt.Line2D([], [], marker="s", ls="none", mfc="white", mec=C_HOLD,
                   mew=2, ms=8, label="held out (truth)"),
        plt.Line2D([], [], marker="x", ls="none", color=C_HOLD, mew=2, ms=9,
                   label="held out (predicted)"),
        plt.Line2D([], [], color=C_FIT, lw=2, label="fit on all K"),
        plt.Line2D([], [], color=C_HOLD, lw=2, ls="--", label="fit on low K only"),
        plt.Line2D([], [], marker="*", ls="none", color=C_FIT, ms=13,
                   label=r"$M^2(K\to\infty)$"),
        plt.Line2D([], [], color=INK_2, lw=1, ls=":", label="published value"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=7, frameon=False,
               fontsize=8, bbox_to_anchor=(0.5, -0.005),
               labelcolor=INK_2, handletextpad=0.5, columnspacing=1.4)

    fig.suptitle("Richardson fits, tested against K the fit never saw  "
                 "(2K = 25-49)", fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))

    outdir = ROOT / "figures"
    outdir.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        p = outdir / f"extrapolation_holdout_python.{ext}"
        fig.savefig(p, dpi=170 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
        print(f"  saved {p}")


if __name__ == "__main__":
    main()
