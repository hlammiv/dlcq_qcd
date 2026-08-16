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

    # 2026 row: resolution traded for Fock depth so the panels span the SAME
    # y-range as the untruncated 1990 originals (a state of L partons sits at
    # ~L^2 at small lambda, so the panel top needs LPN 8/9/10).  nev=250
    # covers the range; the coverage check below verifies that per panel.
    rows = [
        ("a", provider_1990(),
         [(0, 10, 0, 50, "fig2a", None), (1, 13, 0, 60, "fig2b", None),
          (2, 22, 0, 100, "fig2c", None)], LEGACY),
        ("b", None,
         [(0, 34, 8, 50, None, 250), (1, 35, 9, 60, None, 400),
          (2, 34, 10, 100, None, 350)], LAMBDAS),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    for (row, prov, panels, grid), axrow in zip(rows, axes):
        for (B, K, LPN, ymax, dig, nev), ax in zip(panels, axrow):
            p = prov or provider_2026(nev)
            # Near the chiral point the number of states below any fixed
            # M^2 grows toward the whole basis, so "compute everything in
            # the window" stops being a meaningful target for the 2026 row.
            # Instead the panel shows its own coverage: the region above
            # the highest computed level is shaded, so truncation can
            # never read as emptiness (it did, twice, before this).
            ceil_x, ceil_y = [], []
            for lam in grid:
                r = require_cached(p, 3, 1, B, K, lam, LPN=LPN)
                if r is None:
                    continue
                ev = r.eigenvalues[physical_indices(r)]
                if row == "b":
                    ceil_x.append(lam)
                    ceil_y.append(min(float(ev.max()), ymax))
                ax.plot([lam] * len(ev[ev <= ymax]), ev[ev <= ymax], ".",
                        color=C_PRIMARY, markersize=2.0 if row == "b" else 2.5,
                        alpha=0.7 if row == "b" else 0.8)
            if row == "b" and ceil_x and min(ceil_y) < ymax:
                ax.fill_between(ceil_x, ceil_y, ymax, color="0.85",
                                alpha=0.55, linewidth=0, zorder=0)
                ax.plot(ceil_x, ceil_y, "-", color="0.6", linewidth=0.7)
            if dig:
                x, y, _ = digitized(dig)
                ax.plot(x, y, "o", mfc="none", mec=C_1990, markersize=4.5,
                        markeredgewidth=0.9)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, ymax)
            ax.set_title(f"({row}) $B={B}$, $2K={K}$"
                         + ("" if LPN == 0 else f", LPN$={LPN}$"),
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


@register("fig5")
def fig5():
    """Meson excitations: valence + 4q sectors, 1990 (2K=24) vs 2026 (2K=60).

    Rescale factors are the 1990 panel multipliers, applied identically in
    both rows so the digitized markers overlay directly.
    """
    from dlcq.figures import paper_lambda
    from dlcq.observables import require_physical_index, structure_function
    from phase4a_config import MG_WEAK

    from dlcq.fock_weights import fock_weights

    lam = float(paper_lambda(MG_WEAK))
    # (state selector, multiplier, digitized panel).  In panels 1-3 the 4q
    # sector is suppressed and carries the multiplier; in the last panel the
    # state IS 4q-dominant (the two-meson continuum), so the roles swap and
    # the multiplier sits on the valence.  At 2026 resolutions "the 11th
    # state" is not the same physical state, so the last panel selects the
    # first 4q-dominant state by its Fock weight instead of by index.
    # fig5d's digitized markers are calibrated to the panel's RIGHT-hand
    # 0-12 axis (the dual-axis trap of docs/figure-validation.md §2), so no
    # overlay until that frame mapping is re-validated against our units.
    states = [(0, 1e3, "fig5a"), (1, 1e2, "fig5b"), (2, 1e2, "fig5c"),
              ("first-4q", 1e4, None)]
    rows = [("a", provider_1990(), 24, 0), ("b", provider_2026(40), 100, 6)]
    fig, axes = plt.subplots(2, 4, figsize=(14, 7))
    for (row, prov, K, LPN), axrow in zip(rows, axes):
        r = require_cached(prov, 3, 1, 0, K, lam, LPN=LPN)
        for (w, mult, dig), ax in zip(states, axrow):
            idx, label_w = None, w
            if r is not None:
                if w == "first-4q":
                    for lev in range(r.n_eigenvectors):
                        try:
                            i = require_physical_index(r, lev)
                        except IndexError:
                            break
                        fw = fock_weights(r, i)
                        if fw.get(4, 0.0) > 0.5:
                            idx, label_w = i, lev
                            break
                else:
                    try:
                        idx = require_physical_index(r, w)
                    except IndexError:
                        idx = None
            if idx is not None:
                x, q2, _ = structure_function(r, idx, nparton=2)
                _, q4, q4b = structure_function(r, idx, nparton=4)
                if w == "first-4q":
                    ax.plot(x, q4 + q4b, "o-", color=C_PRIMARY, markersize=3,
                            linewidth=1.1, label="4q")
                    ax.plot(x, mult * q2, "s--", mfc="none",
                            color=C_SECONDARY, markersize=3, linewidth=1.1,
                            label=rf"valence $\times{mult:g}$")
                else:
                    ax.plot(x, q2, "o-", color=C_PRIMARY, markersize=3,
                            linewidth=1.1, label="valence")
                    ax.plot(x, mult * (q4 + q4b) / 2, "s--", mfc="none",
                            color=C_SECONDARY, markersize=3, linewidth=1.1,
                            label=rf"4q $\times{mult:g}$")
            if dig and row == "a":
                x, y, fill = digitized(dig)
                for f in (True, False):
                    m = fill == f
                    ax.plot(x[m], y[m], "o" if f else "s", mfc="none",
                            mec=C_1990, markersize=6.5, markeredgewidth=1.0,
                            zorder=5, label="HBP (dig.)" if f else None)
            slabel = ("4q-dominant" if w == "first-4q"
                      else f"state {w + 1}")
            if w == "first-4q" and isinstance(label_w, int):
                slabel += f" (level {label_w + 1})"
            ax.set_title(f"({row}) {slabel}, $2K={K}$"
                         + ("" if LPN == 0 else f", LPN$={LPN}$"),
                         fontsize=9)
            ax.set_xlabel("$x$")
            style(ax)
        axrow[0].set_ylabel("$q(x)$")
        axrow[0].legend(fontsize=7, frameon=False)
    fig.suptitle("FIG. 5: $N=3$ meson states at $m/g=1.6$ — valence and "
                 "four-quark content, 1990 (top) and 2026 (bottom)",
                 fontsize=11)
    fig.tight_layout()
    finish(fig, "fig5_meson_wf")


@register("fig6")
def fig6():
    """Baryon states + the two-baryon state.  The (a)-(c) five-quark panels
    show the known 25-43% disagreement with the 1990 curves deliberately:
    the defect is downstream of the 1990 eigensolve (lost conversion code),
    and panel (d)'s 0.4% agreement is the control."""
    from dlcq.figures import paper_lambda
    from dlcq.observables import require_physical_index, structure_function
    from phase4a_config import MG_WEAK

    lam = float(paper_lambda(MG_WEAK))
    fig, axes = plt.subplots(2, 4, figsize=(14, 7))
    specs = [
        # (state, B, mult, digitized panel), valence nparton / hf nparton
        (0, 1, 1e3, "fig6a", 3, 5), (1, 1, 1e2, "fig6b", 3, 5),
        (2, 1, 1e2, "fig6c", 3, 5), (0, 2, 5e2, "fig6d", 6, 8),
    ]
    rows = [("a", provider_1990(), {1: (21, 0), 2: (24, 0)}, True),
            ("b", None, {1: (99, 5), 2: (50, 8)}, False)]
    for (row, prov, kmap, overlay), axrow in zip(rows, axes):
        for (w, B, mult, dig, npv, nph), ax in zip(specs, axrow):
            K, LPN = kmap[B]
            p = prov or provider_2026(12 if B == 1 else 8)
            r = require_cached(p, 3, 1, B, K, lam, LPN=LPN)
            if r is not None:
                try:
                    idx = require_physical_index(r, w)
                except IndexError:
                    idx = None
                if idx is not None:
                    x, qv, _ = structure_function(r, idx, nparton=npv)
                    _, qh, qhb = structure_function(r, idx, nparton=nph)
                    ax.plot(x, qv, "o-", color=C_PRIMARY, markersize=3,
                            linewidth=1.1, label="valence")
                    ax.plot(x, mult * (qh + qhb) / 2, "s--", mfc="none",
                            color=C_SECONDARY, markersize=3, linewidth=1.1,
                            label=rf"{nph}q $\times{mult:g}$")
            if overlay and dig:
                x, y, fill = digitized(dig)
                for f in (True, False):
                    m = fill == f
                    ax.plot(x[m], y[m], "o" if f else "s", mfc="none",
                            mec=C_1990, markersize=6.5, markeredgewidth=1.0,
                            zorder=5, label="HBP (dig.)" if f else None)
            ax.set_title(f"({row}) $B={B}$ state {w + 1}, $2K={K}$"
                         + ("" if LPN == 0 else f", LPN$={LPN}$"),
                         fontsize=9)
            ax.set_xlabel("$x$")
            style(ax)
        axrow[0].set_ylabel("$q(x)$")
        axrow[0].legend(fontsize=7, frameon=False)
    fig.suptitle("FIG. 6: $N=3$ baryons and the $B=2$ state at $m/g=1.6$ — "
                 "1990 (top; the five-quark disagreement is the documented "
                 "1990 conversion defect) and 2026 (bottom)", fontsize=11)
    fig.tight_layout()
    finish(fig, "fig6_baryon_wf")


@register("fig22")
def fig22():
    """New figure: N=2 meson four-quark probability vs resolution.

    The only published K-scan of higher-Fock content is thesis fig. 22
    (2K = 14 and 20); this extends it to 2K = 100 at both couplings using
    the Fock-sector weights of dlcq/fock_weights.py.
    """
    from dlcq.figures import paper_lambda
    from dlcq.fock_weights import fock_weights
    from phase4a_config import MG_STRONG, MG_WEAK

    prov = provider_2026(4)
    fig, ax = plt.subplots(figsize=(5.2, 4))
    for mg, color, mk in ((MG_WEAK, C_PRIMARY, "o"),
                          (MG_STRONG, C_SECONDARY, "s")):
        lam = float(paper_lambda(mg))
        Ks, P4 = [], []
        for K in (14, 20, 28, 40, 56, 80, 100):
            r = require_cached(prov, 2, 1, 0, K, lam, LPN=4)
            if r is None:
                continue
            w = fock_weights(r, 0)
            Ks.append(K)
            P4.append(w.get(4, 0.0))
        ax.semilogy(Ks, P4, mk + "-", color=color, markersize=5,
                    linewidth=1.2, label=f"$m/g={mg}$")
    ax.set_xlabel("$2K$")
    ax.set_ylabel("four-quark probability $P_4$")
    ax.legend(fontsize=9, frameon=False)
    style(ax)
    fig.suptitle("N=2 meson: four-quark content vs resolution", fontsize=10)
    fig.tight_layout()
    finish(fig, "fig22_p4_vs_K")


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
