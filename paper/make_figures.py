#!/usr/bin/env python3
"""Single figure driver for the paper. Reads, never solves.

Contract (see the plan's "Repo mechanics"):
  * Inputs are exclusively ``data/**``, ``runs/python_cache/**`` (through
    ``PythonProvider``), ``refs/digitized/**`` and ``refs/*.csv``.  Before
    any provider call the driver verifies the cache file exists; on a miss
    it exits listing the missing keys as ``paper/jobs/`` commands.  A paper
    build is therefore deterministic and fast, and missing data surfaces as
    a to-do list rather than a surprise solve.
  * Standard-Hamiltonian series come from the cache; improved series come
    from ``data/*.csv`` — the never-cross-bases rule at the I/O layer.
  * Every output lands in ``paper/figs/`` and embeds provenance (git SHA)
    in the PDF metadata.

Reproduction protocol (uniform): each 1990 figure renders as two rows —
row (a) at the 1990 parameters with the digitized published markers
overlaid, row (b) at the 2026 resolutions of paper/jobs/phase4a_config.py.

Palette (checked: all OKLab pairs >= 15.9, ColorBrewer RdBu extremes, and
every series also differs by linestyle or marker fill):
  this-work valence/primary  #2166ac solid
  this-work secondary        #b2182b dashed
  digitized 1990             #7f7f7f hollow markers

Usage:
    python paper/make_figures.py [--only fig2 fig3 ...] [--list]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIGS = Path(__file__).resolve().parent / "figs"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "paper" / "jobs"))

import matplotlib                                    # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402
import numpy as np                                   # noqa: E402

C_PRIMARY, C_SECONDARY, C_1990 = "#2166ac", "#b2182b", "#7f7f7f"

REGISTRY: dict = {}
_MISSING: list = []


def register(name):
    def deco(fn):
        REGISTRY[name] = fn
        return fn
    return deco


def _git_sha() -> str:
    out = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                         cwd=ROOT, capture_output=True, text=True)
    return out.stdout.strip() or "unknown"


# ── cache access that never solves ─────────────────────────────────────────

def provider_1990():
    """The settings the reproduction cache was built with (exact:fortran)."""
    from dlcq.providers import PythonProvider
    return PythonProvider(ncpus=1)


def provider_2026(nev):
    """The Phase 4a settings (exact:blockwise:sparse:k<nev>)."""
    from dlcq.providers import PythonProvider
    return PythonProvider(ncpus=1, assembly="exact", policy="blockwise",
                          solver="sparse", nev=nev)


def require_cached(prov, N, NF, B, K_code, rlamb, LPN=0):
    """Return the cached DLCQResult, or record the miss and return None."""
    from dlcq.providers import _tag
    tag = _tag(N, NF, B, K_code, rlamb, -1.0, LPN, extra=prov._extra())
    path = prov.cache_dir / f"{tag}.h5"
    if not path.exists():
        _MISSING.append(f"N={N} B={B} 2K={K_code} lam={rlamb:.6g} LPN={LPN} "
                        f"[{prov._extra()}]")
        return None
    return prov.get(N, NF, B, K_code, rlamb, LPN=LPN)


def digitized(panel):
    """Digitized 1990 markers: (x, y, filled) arrays from refs/digitized."""
    import csv
    xs, ys, fill = [], [], []
    with open(ROOT / "refs" / "digitized" / f"{panel}.csv") as fh:
        for r in csv.DictReader(l for l in fh if not l.startswith("#")):
            xs.append(float(r["x"])); ys.append(float(r["y"]))
            fill.append(r.get("filled") == "1")
    return np.array(xs), np.array(ys), np.array(fill)


def finish(fig, name):
    FIGS.mkdir(exist_ok=True)
    fig.savefig(FIGS / f"{name}.pdf",
                metadata={"Creator": f"make_figures.py @ {_git_sha()}"})
    plt.close(fig)
    print(f"  {name}.pdf")


def style(ax):
    ax.grid(alpha=0.15, linewidth=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# ── figures ────────────────────────────────────────────────────────────────

@register("fig2")
def fig2():
    """Spectra vs coupling: (a) 1990 resolutions + digitized, (b) 2026."""
    from dlcq.observables import physical_indices
    from phase4a_config import LAMBDAS

    # The 1990-resolution sweep in the cache uses the legacy grid that
    # matches the preserved Fortran runs (runs/fig2), not the modern list.
    LEGACY = [0.066667, 0.216667, 0.366667, 0.516667, 0.666667, 0.816667,
              0.966667]

    rows = [
        ("a", provider_1990(), [(0, 10, 0, 50, "fig2a"), (1, 13, 0, 60, "fig2b"),
                                (2, 22, 0, 100, "fig2c")], None, LEGACY),
        ("b", None, [(0, 50, 6, 50, None), (1, 51, 5, 60, None),
                     (2, 50, 8, 100, None)], 40, LAMBDAS),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    NLEV = 14           # levels drawn in the 2026 row: connected, not a cloud
    for (row, prov, panels, nev, grid), axrow in zip(rows, axes):
        for (B, K, LPN, ymax, dig), ax in zip(panels, axrow):
            p = prov or provider_2026(nev)
            if row == "a":
                # Full spectra, as the 1990 figure drew them.
                for lam in grid:
                    r = require_cached(p, 3, 1, B, K, lam, LPN=LPN)
                    if r is None:
                        continue
                    ev = r.eigenvalues[physical_indices(r)]
                    ax.plot([lam] * len(ev[ev <= ymax]), ev[ev <= ymax], ".",
                            color=C_PRIMARY, markersize=2.5, alpha=0.8)
                ax.set_ylim(0, ymax)
            else:
                # nev bounds what a sparse solve computed, so a full-spectrum
                # cloud would really be the bottom-nev envelope (an artifact:
                # we drew it, saw the arch, and replaced it).  The 2026 panel
                # shows the lowest NLEV levels as curves instead.
                levels = [[] for _ in range(NLEV)]
                lams = []
                for lam in grid:
                    r = require_cached(p, 3, 1, B, K, lam, LPN=LPN)
                    if r is None:
                        continue
                    ev = r.eigenvalues[physical_indices(r)][:NLEV]
                    lams.append(lam)
                    for i in range(NLEV):
                        levels[i].append(ev[i] if i < len(ev) else np.nan)
                for i in range(NLEV):
                    ax.plot(lams, levels[i], "-", color=C_PRIMARY,
                            linewidth=0.9, alpha=0.85)
            if dig:
                x, y, _ = digitized(dig)
                ax.plot(x, y, "o", mfc="none", mec=C_1990, markersize=4.5,
                        markeredgewidth=0.9)
            ax.set_xlim(0, 1)
            ax.set_title(f"({row}) $B={B}$, $2K={K}$"
                         + ("" if LPN == 0 else f", LPN$={LPN}$")
                         + ("" if row == "a" else f", lowest {NLEV}"),
                         fontsize=10)
            style(ax)
            if B == 0:
                ax.set_ylabel(r"$M^2/(m^2+g^2/\pi)$")
            ax.set_xlabel(r"$\lambda$")
    fig.suptitle("FIG. 2: spectra vs coupling — 1990 reproduction (top, "
                 "markers digitized from HBP) and 2026 (bottom)", fontsize=11)
    fig.tight_layout()
    finish(fig, "fig2_spectra")


@register("fig3")
def fig3():
    """Valence structure functions: (a) 2K=14/15 + digitized, (b) 2K=70/71."""
    from dlcq.figures import paper_lambda
    from dlcq.observables import structure_function
    from phase4a_config import MG_STRONG, MG_WEAK

    rows = [("a", provider_1990(), {0: (14, 0), 1: (15, 0)},
             {0: "fig3a", 1: "fig3b"}),
            ("b", provider_2026(4), {0: (70, 4), 1: (71, 5)}, {})]
    fig, axes = plt.subplots(2, 2, figsize=(9, 7.5))
    for (row, prov, cfg, digs), axrow in zip(rows, axes):
        for B, ax in zip((0, 1), axrow):
            K, LPN = cfg[B]
            for mg, mfc in ((MG_WEAK, C_PRIMARY), (MG_STRONG, "none")):
                r = require_cached(prov, 3, 1, B, K, float(paper_lambda(mg)),
                                   LPN=LPN)
                if r is None:
                    continue
                x, q, qbar = structure_function(r, 0)
                ax.plot(x, q, "o-" if mfc != "none" else "s--",
                        color=C_PRIMARY, mfc=mfc, markersize=4,
                        linewidth=1.2, label=f"$m/g={mg}$")
            if B in digs:
                x, y, fill = digitized(digs[B])
                for f, mk in ((True, "o"), (False, "s")):
                    m = fill == f
                    # On top and larger, so the 1990 ring shows around our
                    # marker instead of vanishing underneath it.
                    ax.plot(x[m], y[m], mk, mfc="none", mec=C_1990,
                            markersize=7.5, markeredgewidth=1.1, zorder=5,
                            label="HBP (dig.)" if f else None)
            ax.set_title(f"({row}) {'meson' if B == 0 else 'baryon'}, "
                         f"$2K={K}$", fontsize=10)
            ax.set_xlabel("$x$")
            if B == 0:
                ax.set_ylabel("$q(x)$")
            ax.legend(fontsize=8, frameon=False)
            style(ax)
    fig.suptitle("FIG. 3: valence structure functions, $N=3$ — 1990 (top) "
                 "and 2026 (bottom)", fontsize=11)
    fig.tight_layout()
    finish(fig, "fig3_valence")


# ── main ───────────────────────────────────────────────────────────────────

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args(argv)

    if args.list:
        print(f"{len(REGISTRY)} registered figure(s):")
        for name in sorted(REGISTRY):
            print(f"  {name}")
        return 0

    names = args.only or sorted(REGISTRY)
    unknown = [n for n in names if n not in REGISTRY]
    if unknown:
        ap.error(f"unknown figure(s): {unknown}; use --list")

    print(f"[{_git_sha()}] building: {' '.join(names)}")
    for name in names:
        REGISTRY[name]()

    if _MISSING:
        print(f"\n{len(_MISSING)} cache miss(es) — run these on Lenore via "
              f"paper/jobs/phase4a_lenore.py (or extend phase4a_config):")
        for m in _MISSING:
            print(f"  {m}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
