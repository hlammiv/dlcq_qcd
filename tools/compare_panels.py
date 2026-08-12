#!/usr/bin/env python3
"""Paper | Fortran | Python, side by side, for every reproducible panel.

Three columns that are supposed to be the same picture:

    left    the panel cropped from the scanned page, exactly as published
    middle  the same quantity computed by fortran/qcdf.f
    right   the same quantity computed by the Python port

The computed panels are drawn on the axis limits read off the paper's own
frame, so a shape or scale difference is visible rather than hidden by
autoscaling.  Digitized markers are overlaid on both computed panels as hollow
grey circles, which is what the numerical comparison in
docs/figure-validation.md actually measures.

Usage::

    python tools/compare_panels.py                 # all panels
    python tools/compare_panels.py --panel fig6a
    python tools/compare_panels.py --outdir figures/compare
"""

from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))


# panel -> what physics it shows.
#   B, K       the run
#   state      index into the PHYSICAL spectrum (spurious modes removed)
#   series     list of (parton count, multiplier, marker kind, m/g, parton,
#              label)
# The multipliers are the paper's own x10^n legend entries.  Note Fig. 5(d)
# puts its x10^4 on the FILLED (valence) curve, not the open one.
#
# ``marker kind`` names the shape the paper draws the series with, which is
# what ties a computed curve to the right rows of the digitized CSV -- several
# panels use more than two shapes.  And Fig. 3 and Fig. 4 do not show two Fock
# sectors at one mass but one sector at two masses, which is what ``m/g``
# carries; ``parton`` picks q(x) or qbar(x).
PANEL_PHYSICS = {
    "fig3a": dict(B=0, K=14, state=0,
                  series=[(2, 1.0, "filled", 1.6, "q", r"$q\bar q$, $m/g=1.6$"),
                          (2, 1.0, "open", 0.1, "q", r"$q\bar q$, $m/g=0.1$")]),
    "fig3b": dict(B=1, K=15, state=0,
                  series=[(3, 1.0, "filled", 1.6, "q", "qqq, $m/g=1.6$"),
                          (3, 1.0, "open", 0.1, "q", "qqq, $m/g=0.1$")]),
    "fig4a": dict(B=0, K=14, state=0,
                  series=[(4, 1e4, "filled", 1.6, "q",
                           r"$q\bar qq\bar q\ (\times10^4)$, $m/g=1.6$"),
                          (4, 1e4, "open", 0.1, "q",
                           r"$q\bar qq\bar q\ (\times10^4)$, $m/g=0.1$")]),
    "fig4b": dict(B=1, K=15, state=0,
                  series=[(5, 1e3, "filled", 1.6, "q",
                           r"$q\ (\times10^3)$, $m/g=1.6$"),
                          (5, 1e2, "open", 0.1, "q",
                           r"$q\ (\times10^2)$, $m/g=0.1$"),
                          (5, 1e3, "triangle_down", 1.6, "qbar",
                           r"$\bar q\ (\times10^3)$, $m/g=1.6$"),
                          (5, 1e2, "triangle", 0.1, "qbar",
                           r"$\bar q\ (\times10^2)$, $m/g=0.1$")]),
    "fig4c": dict(B=1, K=15, state=0,
                  series=[(7, 1e4, "filled", 1.6, "q",
                           r"$7$-parton $(\times10^4)$, $m/g=1.6$"),
                          (7, 1e7, "open", 0.1, "q",
                           r"$7$-parton $(\times10^7)$, $m/g=0.1$")]),
    "fig5a": dict(B=0, K=24, state=0,
                  series=[(2, 1.0, "filled", 1.6, "q", r"$q\bar q$"),
                          (4, 1e3, "open", 1.6, "q",
                           r"$q\bar qq\bar q\ (\times10^3)$")]),
    "fig5b": dict(B=0, K=24, state=1,
                  series=[(2, 1.0, "filled", 1.6, "q", r"$q\bar q$"),
                          (4, 1e2, "open", 1.6, "q",
                           r"$q\bar qq\bar q\ (\times10^2)$")]),
    "fig5c": dict(B=0, K=24, state=2,
                  series=[(2, 1.0, "filled", 1.6, "q", r"$q\bar q$"),
                          (4, 1e2, "open", 1.6, "q",
                           r"$q\bar qq\bar q\ (\times10^2)$")]),
    "fig5d": dict(B=0, K=24, state=10,
                  series=[(2, 1e4, "filled", 1.6, "q", r"$q\bar q\ (\times10^4)$"),
                          (4, 1.0, "open", 1.6, "q", r"$q\bar qq\bar q$")]),
    "fig6a": dict(B=1, K=21, state=0,
                  series=[(3, 1.0, "filled", 1.6, "q", "qqq"),
                          (5, 1e3, "open", 1.6, "q",
                           r"$qqqq\bar q\ (\times10^3)$")]),
    "fig6b": dict(B=1, K=21, state=1,
                  series=[(3, 1.0, "filled", 1.6, "q", "qqq"),
                          (5, 1e2, "open", 1.6, "q",
                           r"$qqqq\bar q\ (\times10^2)$")]),
    "fig6c": dict(B=1, K=21, state=2,
                  series=[(3, 1.0, "filled", 1.6, "q", "qqq"),
                          (5, 1e2, "open", 1.6, "q",
                           r"$qqqq\bar q\ (\times10^2)$")]),
    # 2K = 24, not the 22 of Fig. 2(c): the caption never states it, and the
    # marker positions, the quark-number sum rule and the mean momentum all
    # say 24.  See the panel comment in tools/digitize.py.
    "fig6d": dict(B=2, K=24, state=0,
                  series=[(6, 1.0, "filled", 1.6, "q", "6q"),
                          (8, 5e2, "open", 1.6, "q",
                           r"$6qq\bar q\ (\times5\!\times\!10^2)$"),
                          (8, 1e3, "triangle", 1.6, "qbar",
                           r"$6qq\bar q$, $\bar q(x)\ (\times10^3)$")]),
}


# How each marker kind is drawn on the computed panels, so that a curve and
# the digitized points it should pass through look alike: (style, face, colour).
MARKER_STYLE = {
    "filled": ("-o", "k", "k"),
    "open": ("--s", "none", "crimson"),
    "triangle": ("--^", "none", "darkorchid"),
    "triangle_down": ("--v", "none", "teal"),
}


def paper_crop(panel, dpi=600, pages_dir=None):
    """The published panel as an image array."""
    from PIL import Image
    from render_pages import DEFAULT_OUT, PDF, render

    pages_dir = Path(pages_dir or DEFAULT_OUT)
    hits = sorted(pages_dir.glob(f"p{panel.page:02d}_{dpi}dpi*.png"))
    if not hits:
        if not PDF.exists():
            return None
        hits = [render(PDF, panel.page, dpi, pages_dir)]
    img = Image.open(hits[0]).convert("L")
    w, h = img.size
    l, t, r, b = panel.bbox
    return np.asarray(img.crop((int(l * w), int(t * h), int(r * w), int(b * h))))


def digitized(name):
    """Digitized markers for a panel, keyed by marker kind.

    Splitting on the ``kind`` column rather than on filled/open matters wherever
    a panel plots more than two series: Fig. 6(d) has circles *and* triangles,
    all of which are "not filled".
    """
    path = ROOT / "refs" / "digitized" / f"{name}.csv"
    out = {}
    if not path.exists():
        return out
    header = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        if line.startswith("x,"):
            header = line.split(",")
            continue
        parts = line.split(",")
        if len(parts) < 3:
            continue
        row = dict(zip(header, parts))
        kind = row.get("kind") or ("filled" if parts[2] == "1" else "open")
        out.setdefault(kind, []).append(
            (float(parts[0]), float(parts[1]), row.get("conf", "probe")))
    return {k: sorted(v) for k, v in out.items()}


def computed(provider, name, phys):
    """(x, [(y, kind, label)]) for one panel from one solver.

    One run per distinct m/g in the panel: Figs. 3 and 4 overlay two masses.
    """
    from dlcq.observables import structure_function, physical_indices
    from dlcq.figures import paper_lambda

    runs, out, xs = {}, [], None
    for npart, scale, kind, mg, parton, label in phys["series"]:
        if mg not in runs:
            runs[mg] = provider.get(3, 1, phys["B"], phys["K"],
                                    paper_lambda(mg))
        r = runs[mg]
        idx = physical_indices(r)
        if phys["state"] >= idx.size:
            return None, []
        x, q, qbar = structure_function(r, int(idx[phys["state"]]),
                                        nparton=npart)
        xs = x
        out.append(((qbar if parton == "qbar" else q) * scale, kind, label))
    return xs, out


def series_stats(name, provider, phys, include_coincident=False):
    """Per-series agreement between the trace and one solver.

    Reported against the series' own peak, not point by point: these curves run
    down to zero, where a relative error is meaningless.  ``scale`` is the
    least-squares multiplier that would best map ours onto the paper's, so a
    mis-stated power of ten shows up as a number far from 1 while the shape
    stays visible in ``corr``.
    """
    x, series = computed(provider, name, phys)
    if x is None:
        return []
    dig = digitized(name)
    rows = []
    for (q, kind, label) in series:
        pts = [p for p in dig.get(kind, [])
               if include_coincident or p[2] != "coincident"]
        if len(pts) < 3:
            rows.append(dict(series=label, kind=kind, n=len(pts)))
            continue
        px = np.array([p[0] for p in pts])
        py = np.array([p[1] for p in pts])
        j = [int(np.argmin(np.abs(x - v))) for v in px]
        ours = q[j]
        peak = max(abs(py).max(), 1e-30)
        scale = float((py @ ours) / (ours @ ours)) if ours @ ours else float("nan")
        corr = float(np.corrcoef(py, ours)[0, 1]) if py.std() and ours.std() \
            else float("nan")
        rows.append(dict(series=label, kind=kind, n=len(pts),
                         corr=corr, scale=scale,
                         max_dev=float(np.abs(ours - py).max() / peak),
                         rms_dev=float(np.sqrt(((ours - py) ** 2).mean()) / peak)))
    return rows


def print_table(names, include_coincident=False):
    from dlcq.providers import PythonProvider

    pp = PythonProvider(ncpus=8)
    print(f"{'panel':7s} {'series':44s} {'n':>3s} {'corr':>7s} "
          f"{'scale':>8s} {'rms/peak':>9s} {'max/peak':>9s}")
    for name in names:
        phys = PANEL_PHYSICS.get(name)
        if phys is None:
            continue
        try:
            rows = series_stats(name, pp, phys, include_coincident)
        except Exception as exc:                       # noqa: BLE001
            print(f"{name:7s} unavailable: {exc}")
            continue
        for r in rows:
            if "corr" not in r:
                print(f"{name:7s} {r['series'][:44]:44s} {r['n']:3d}  "
                      "too few digitized points")
                continue
            print(f"{name:7s} {r['series'][:44]:44s} {r['n']:3d} "
                  f"{r['corr']:7.4f} {r['scale']:8.4f} "
                  f"{100*r['rms_dev']:8.2f}% {100*r['max_dev']:8.2f}%")


def make_figure(name, outdir, dpi=600, pages_dir=None):
    from digitize import PANELS

    panel = PANELS.get(name)
    phys = PANEL_PHYSICS.get(name)
    if panel is None or phys is None:
        return None

    from dlcq.providers import FortranProvider, PythonProvider

    fp = FortranProvider(allow_run=True, extra_search=[ROOT / "python"])
    pp = PythonProvider(ncpus=8)

    img = paper_crop(panel, dpi=dpi, pages_dir=pages_dir)
    dig = digitized(name)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    # ── the published panel ──
    if img is not None:
        axes[0].imshow(img, cmap="gray", aspect="auto")
    else:
        axes[0].text(0.5, 0.5, "paper PDF not available\n(see CITATION.md)",
                     ha="center", va="center", transform=axes[0].transAxes)
    axes[0].set_xticks([]); axes[0].set_yticks([])
    axes[0].set_title("paper (as published)", fontsize=11)

    # ── the two solvers, on the paper's own axis limits ──
    for ax, prov, label in ((axes[1], fp, "fortran"), (axes[2], pp, "python")):
        try:
            x, series = computed(prov, name, phys)
        except Exception as exc:
            # Wrapped, because an unwrapped subprocess error is one enormous
            # line and matplotlib's tight_layout squeezes the sibling panels
            # down to slivers to make room for it.
            msg = textwrap.fill(f"{type(exc).__name__}: {exc}", 46,
                                max_lines=6, placeholder=" ...")
            ax.text(0.5, 0.5, f"unavailable\n\n{msg}", ha="center",
                    va="center", transform=ax.transAxes, fontsize=7)
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(label, fontsize=11)
            continue
        if x is None:
            ax.text(0.5, 0.5, "state unavailable", ha="center", va="center",
                    transform=ax.transAxes)
        else:
            for q, kind, lab in series:
                ax.plot(x, q, MARKER_STYLE[kind][0],
                        markerfacecolor=MARKER_STYLE[kind][1],
                        color=MARKER_STYLE[kind][2],
                        markersize=4, lw=1.1, label=lab)
        first = True
        for kind, pts in dig.items():
            ax.plot([p[0] for p in pts], [p[1] for p in pts],
                    MARKER_STYLE.get(kind, ("o", "none", "0.55"))[0][-1],
                    mfc="none", mec="0.55", ms=9, mew=1.2, ls="none",
                    label="paper, digitized" if first else None, zorder=0)
            first = False
        ax.set_xlim(*panel.xlim)
        ax.set_ylim(*panel.ylim)
        ax.set_xlabel("x = k/K"); ax.set_ylabel("q(x)")
        ax.set_title(label, fontsize=11)
        ax.legend(fontsize=7)

    fig.suptitle(f"{panel.description}   —   these three should be the same plot",
                 fontsize=12)
    fig.tight_layout()
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"compare_{name}.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", nargs="+", default=sorted(PANEL_PHYSICS))
    ap.add_argument("--outdir", type=Path, default=ROOT / "figures" / "compare")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--pages-dir", type=Path, default=None)
    ap.add_argument("--table", action="store_true",
                    help="print the per-series agreement instead of drawing")
    ap.add_argument("--with-coincident", action="store_true",
                    help="include points recorded as buried under another "
                         "series (conf=coincident)")
    args = ap.parse_args(argv)

    if args.table:
        print_table(args.panel, include_coincident=args.with_coincident)
        return 0

    for name in args.panel:
        path = make_figure(name, args.outdir, dpi=args.dpi,
                           pages_dir=args.pages_dir)
        print(f"  {name}: {path}" if path else f"  {name}: skipped")
    return 0


if __name__ == "__main__":
    sys.exit(main())
