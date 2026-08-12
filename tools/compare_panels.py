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
#   series     list of (parton count, multiplier, marker-is-filled, label)
# The multipliers are the paper's own x10^n legend entries.  Note Fig. 5(d)
# puts its x10^4 on the FILLED (valence) curve, not the open one.
PANEL_PHYSICS = {
    "fig3a": dict(B=0, K=14, state=0,
                  series=[(2, 1.0, True, r"$q\bar q$ valence")]),
    "fig3b": dict(B=1, K=15, state=0,
                  series=[(3, 1.0, True, "qqq valence")]),
    "fig5a": dict(B=0, K=24, state=0,
                  series=[(2, 1.0, True, r"$q\bar q$"),
                          (4, 1e3, False, r"$q\bar qq\bar q\ (\times10^3)$")]),
    "fig5b": dict(B=0, K=24, state=1,
                  series=[(2, 1.0, True, r"$q\bar q$"),
                          (4, 1e2, False, r"$q\bar qq\bar q\ (\times10^2)$")]),
    "fig5c": dict(B=0, K=24, state=2,
                  series=[(2, 1.0, True, r"$q\bar q$"),
                          (4, 1e2, False, r"$q\bar qq\bar q\ (\times10^2)$")]),
    "fig5d": dict(B=0, K=24, state=10,
                  series=[(2, 1e4, True, r"$q\bar q\ (\times10^4)$"),
                          (4, 1.0, False, r"$q\bar qq\bar q$")]),
    "fig6a": dict(B=1, K=21, state=0,
                  series=[(3, 1.0, True, "qqq"),
                          (5, 1e3, False, r"$qqqq\bar q\ (\times10^3)$")]),
    "fig6b": dict(B=1, K=21, state=1,
                  series=[(3, 1.0, True, "qqq"),
                          (5, 1e2, False, r"$qqqq\bar q\ (\times10^2)$")]),
    "fig6c": dict(B=1, K=21, state=2,
                  series=[(3, 1.0, True, "qqq"),
                          (5, 1e2, False, r"$qqqq\bar q\ (\times10^2)$")]),
    "fig6d": dict(B=2, K=22, state=0,
                  series=[(6, 1.0, True, "6q"),
                          (8, 5e2, False, r"$6qq\bar q\ (\times5\!\times\!10^2)$")]),
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
    """Digitized markers for a panel, split into filled and open."""
    path = ROOT / "refs" / "digitized" / f"{name}.csv"
    if not path.exists():
        return [], []
    filled, open_ = [], []
    for line in path.read_text().splitlines():
        if line.startswith("#") or line.startswith("x,"):
            continue
        parts = line.split(",")
        if len(parts) < 3:
            continue
        x, y, f = float(parts[0]), float(parts[1]), parts[2] == "1"
        (filled if f else open_).append((x, y))
    return sorted(filled), sorted(open_)


def computed(provider, name, phys):
    """(x, [(y, filled, label)]) for one panel from one solver."""
    from dlcq.observables import structure_function, physical_indices
    from dlcq.figures import paper_lambda

    r = provider.get(3, 1, phys["B"], phys["K"], paper_lambda(1.6))
    idx = physical_indices(r)
    if phys["state"] >= idx.size:
        return None, []
    out, xs = [], None
    for npart, scale, filled, label in phys["series"]:
        x, q, _ = structure_function(r, int(idx[phys["state"]]), nparton=npart)
        xs = x
        out.append((q * scale, filled, label))
    return xs, out


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
    dig_f, dig_o = digitized(name)

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
            ax.text(0.5, 0.5, f"unavailable:\n{exc}", ha="center", va="center",
                    transform=ax.transAxes, fontsize=8)
            ax.set_title(label, fontsize=11)
            continue
        if x is None:
            ax.text(0.5, 0.5, "state unavailable", ha="center", va="center",
                    transform=ax.transAxes)
        else:
            for q, filled, lab in series:
                ax.plot(x, q, "-o" if filled else "--s",
                        markerfacecolor="k" if filled else "none",
                        color="k" if filled else "crimson",
                        markersize=4, lw=1.1, label=lab)
        if dig_f:
            ax.plot(*zip(*dig_f), "o", mfc="none", mec="0.55", ms=9, mew=1.2,
                    label="paper, digitized", zorder=0)
        if dig_o:
            ax.plot(*zip(*dig_o), "s", mfc="none", mec="0.75", ms=9, mew=1.2,
                    zorder=0)
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


# Figures that cannot be compared marker-by-marker -- dense level bundles
# (Fig. 2) or vertically offset families (Figs. 7, 8(b)).  For these we place
# the scanned original beside our two generated figures so the shapes can be
# judged directly.  Page/crop is the article; the thesis twin is noted.
WHOLE_FIGURE = {
    "fig2": dict(page=3, bbox=(0.58, 0.08, 0.86, 0.49), stem="fig2_spectra",
                 title="FIG. 2  spectra vs coupling, B = 0, 1, 2 "
                       "(thesis Fig. 4 confirms 2K = 10, 13, 22)"),
    "fig7": dict(page=7, bbox=(0.06, 0.06, 0.50, 0.34), stem="fig7_masses",
                 title="FIG. 7  extrapolated masses, N = 2, 3, 4"),
    "fig8": dict(page=7, bbox=(0.55, 0.05, 0.97, 0.46), stem="fig8_large_N",
                 title="FIG. 8  comparison with large N "
                       "(thesis Fig. 7 is the same plot)"),
}


def make_whole_figure(key, outdir, dpi=600, pages_dir=None):
    """Scanned figure beside our fortran and python versions of it."""
    from PIL import Image
    from render_pages import DEFAULT_OUT, PDF, render

    spec = WHOLE_FIGURE[key]
    pages_dir = Path(pages_dir or DEFAULT_OUT)
    hits = sorted(pages_dir.glob(f"p{spec['page']:02d}_{dpi}dpi*.png"))
    if not hits:
        if not PDF.exists():
            return None
        hits = [render(PDF, spec["page"], dpi, pages_dir)]
    img = Image.open(hits[0]).convert("L")
    w, h = img.size
    l, t, r, b = spec["bbox"]
    crop = np.asarray(img.crop((int(l * w), int(t * h), int(r * w), int(b * h))))

    figdir = ROOT / "figures"
    panels = [("paper (as published)", crop)]
    for src in ("fortran", "python"):
        png = figdir / f"{spec['stem']}_{src}.png"
        panels.append((src, np.asarray(Image.open(png)) if png.exists() else None))

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    for ax, (label, arr) in zip(axes, panels):
        if arr is None:
            ax.text(0.5, 0.5, "not generated\n(run dlcq.figures)", ha="center",
                    va="center", transform=ax.transAxes)
        else:
            ax.imshow(arr, cmap="gray" if arr.ndim == 2 else None, aspect="auto")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(label, fontsize=11)
    fig.suptitle(spec["title"] + "   —   these three should show the same physics",
                 fontsize=12)
    fig.tight_layout()
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"compare_{key}.png"
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", nargs="+", default=sorted(PANEL_PHYSICS))
    ap.add_argument("--outdir", type=Path, default=ROOT / "figures" / "compare")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--pages-dir", type=Path, default=None)
    ap.add_argument("--whole", action="store_true",
                    help="also emit whole-figure comparisons for Figs 2, 7, 8")
    args = ap.parse_args(argv)

    for name in args.panel:
        path = make_figure(name, args.outdir, dpi=args.dpi,
                           pages_dir=args.pages_dir)
        print(f"  {name}: {path}" if path else f"  {name}: skipped")

    if args.whole:
        for key in WHOLE_FIGURE:
            path = make_whole_figure(key, args.outdir, dpi=args.dpi,
                                     pages_dir=args.pages_dir)
            print(f"  {key}: {path}" if path else f"  {key}: skipped")
    return 0


if __name__ == "__main__":
    sys.exit(main())
