#!/usr/bin/env python3
"""Trace curve data out of the scanned figures.

The paper carries no vector content (see ``tools/render_pages.py``), so the only
way to per-point numbers is pixel tracing.  The scan is 1-bit at ~612 ppi, which
helps: axis frames are solid runs, ticks are short perpendicular marks, and the
paper's filled/open circle markers separate cleanly by fill ratio.

Calibration is semi-automatic and auditable.  You supply the panel's bounding
box and the *values* of its axis ticks; the code locates the tick *pixels* and
fits the mapping.  Everything it used is written to a provenance JSON beside the
CSV, so any number can be traced back to the pixels it came from.

Self-validation
---------------
Figs. 5 and 6 state their K (2K = 24 and 21).  Because momenta are odd
integers, a structure-function panel's marker positions obey ``x_min = 1/K`` and
``dx = 2/K``, so ``dlcq.units.infer_K_from_x_grid`` recovers K from the trace
alone.  Running that on a panel whose K is known is the check that the
calibration is right -- and once it passes there, the same procedure determines
the K for Figs. 3 and 4, which the paper never states.

Usage::

    python tools/digitize.py --list
    python tools/digitize.py --panel fig6a --check
    python tools/digitize.py --panel fig6a --out refs/digitized/
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))


# ──────────────────────────────────────────────────────────────────────────
# Panel description
# ──────────────────────────────────────────────────────────────────────────

@dataclass
class Panel:
    """One plot panel and how to read its axes.

    ``bbox`` is a fractional (left, top, right, bottom) crop of the rendered
    page -- fractions rather than pixels so it survives a change of dpi.
    ``xticks``/``yticks`` are the tick *values* printed on the axes, ascending.
    """

    name: str
    page: int                      # PDF page (1-based)
    bbox: tuple                    # fractional (l, t, r, b)
    xlim: tuple                    # data values at the left and right frame edges
    ylim: tuple                    # data values at the bottom and top frame edges
    xticks: tuple = ()             # major tick values, for cross-checking only
    yticks: tuple = ()
    description: str = ""
    # Expected K, when the paper states it -- used by --check.
    expected_K: int | None = None
    B: int = 0
    N: int = 3
    # If set, the vertical scale is fixed by the exact sum rule
    # int q dx = quark_number rather than by reading tick labels.  Non-circular:
    # it uses only an identity, never the solver.  1 for a meson's quark
    # distribution, N*B for a baryon's.
    quark_number: float | None = None
    sector: str = ""


PANELS = {
    # ── FIG. 3 (page 4): valence structure functions, K NOT stated ──────────
    "fig3a": Panel(name="fig3a", page=4, bbox=(0.1149, 0.7900, 0.2765, 0.8813),
                   xlim=(0.0, 1.0), ylim=(0.0, 3.5),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   description="FIG. 3(a) SU(3) meson valence; K not stated",
                   B=0, N=3, quark_number=1.0, sector="valence"),
    "fig3b": Panel(name="fig3b", page=4, bbox=(0.3061, 0.7906, 0.4677, 0.8818),
                   xlim=(0.0, 1.0), ylim=(0.0, 11.25),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   description="FIG. 3(b) SU(3) baryon valence; K not stated",
                   B=1, N=3, quark_number=3.0, sector="valence"),

    # ── FIG. 4 (page 4): higher-Fock sectors, scaled by 10^n ────────────────
    "fig4a": Panel(name="fig4a", page=4, bbox=(0.5376, 0.0949, 0.6999, 0.1902),
                   xlim=(0.0, 1.0), ylim=(0.0, 21.0),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   description="FIG. 4(a) meson, one extra qqbar pair (x10^4/10^2)",
                   B=0, N=3, sector="higher-fock"),
    "fig4b": Panel(name="fig4b", page=4, bbox=(0.7304, 0.0961, 0.8929, 0.1911),
                   xlim=(0.0, 1.0), ylim=(0.0, 22.5),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   description="FIG. 4(b) baryon, one extra pair, incl. antiquarks",
                   B=1, N=3, sector="higher-fock"),
    "fig4c": Panel(name="fig4c", page=4, bbox=(0.5372, 0.1811, 0.6994, 0.2765),
                   xlim=(0.0, 1.0), ylim=(0.0, 11.25),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   description="FIG. 4(c) baryon, two extra pairs",
                   B=1, N=3, sector="higher-fock"),

    # ── FIG. 5 (page 4): meson spectrum, 2K = 24 ───────────────────────────
    "fig5a": Panel(name="fig5a", page=4, bbox=(0.5478, 0.7005, 0.7149, 0.7973),
                   xlim=(0.0, 1.0), ylim=(0.0, 3.6),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   description="FIG. 5(a) 1st meson state, m/g=1.6, 2K=24",
                   expected_K=24, B=0, N=3, quark_number=1.0, sector="valence"),
    "fig5b": Panel(name="fig5b", page=4, bbox=(0.7037, 0.7005, 0.8708, 0.7973),
                   xlim=(0.0, 1.0), ylim=(0.0, 12.0),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   description="FIG. 5(b) 2nd meson state, m/g=1.6, 2K=24",
                   expected_K=24, B=0, N=3, quark_number=1.0, sector="valence"),
    "fig5c": Panel(name="fig5c", page=4, bbox=(0.5478, 0.7873, 0.7106, 0.8836),
                   xlim=(0.0, 1.0), ylim=(0.0, 3.6),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   description="FIG. 5(c) 3rd meson state, m/g=1.6, 2K=24",
                   expected_K=24, B=0, N=3, quark_number=1.0, sector="valence"),
    "fig5d": Panel(name="fig5d", page=4, bbox=(0.6994, 0.7873, 0.8665, 0.8836),
                   xlim=(0.0, 1.0), ylim=(0.0, 12.0),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   description="FIG. 5(d) 11th meson state, m/g=1.6, 2K=24",
                   expected_K=24, B=0, N=3, quark_number=1.0, sector="valence"),

    # ── FIG. 6 (page 5): baryon spectrum, 2K = 21 ──────────────────────────
    # The y frame top is NOT a labelled value: ticks sit at height-fractions
    # 0.184/0.365/0.540/0.720/0.895 and the labelled 10.5 is the one at 0.720,
    # so the step is 10.5/4 = 2.625 and the top is 10.5/0.720 = 14.77.
    # Reading the top as 21 inflates every y by 1.42x.
    "fig6a": Panel(name="fig6a", page=5, bbox=(0.1257, 0.1108, 0.2853, 0.2030),
                   xlim=(0.0, 1.0), ylim=(0.0, 14.77),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   yticks=(2.625, 5.25, 7.875, 10.5, 13.125),
                   description="FIG. 6(a) 1st baryon state, m/g=1.6, 2K=21",
                   expected_K=21, B=1, N=3, quark_number=3.0, sector="valence"),
    "fig6b": Panel(name="fig6b", page=5, bbox=(0.2737, 0.1108, 0.4337, 0.2030),
                   xlim=(0.0, 1.0), ylim=(0.0, 14.77),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   description="FIG. 6(b) 2nd baryon state, m/g=1.6, 2K=21",
                   expected_K=21, B=1, N=3, quark_number=3.0, sector="valence"),
    "fig6c": Panel(name="fig6c", page=5, bbox=(0.1257, 0.2141, 0.2853, 0.3059),
                   xlim=(0.0, 1.0), ylim=(0.0, 14.77),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   description="FIG. 6(c) 3rd baryon state, m/g=1.6, 2K=21",
                   expected_K=21, B=1, N=3, quark_number=3.0, sector="valence"),
    "fig6d": Panel(name="fig6d", page=5, bbox=(0.2737, 0.2141, 0.4337, 0.3059),
                   xlim=(0.0, 0.6), ylim=(0.0, 48.0),
                   description="FIG. 6(d) 1st B=2 state, m/g=1.6, 2K=22",
                   expected_K=22, B=2, N=3, quark_number=6.0, sector="valence"),

    # ── FIG. 8 (page 7): meson mass vs m/g, with Hamer's SU(2) lattice points ──
    "fig8a": Panel(name="fig8a", page=7, bbox=(0.5943, 0.1381, 0.8324, 0.2760),
                   # 8 x ticks every 0.25 (0..1.75) and 9 y ticks every 0.5
                   # (0..4.0); the frame edges sit just beyond the last of each.
                   xlim=(0.0, 1.756), ylim=(0.0, 4.02),
                   xticks=(0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75),
                   yticks=(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0),
                   description="FIG. 8(a) meson mass M/g vs m/g; filled circles are "
                               "Hamer's SU(2) lattice results (Nucl. Phys. B195, 503)",
                   B=0, N=2),

    # ── FIG. 2 (page 3): eigenvalue spectra vs coupling ─────────────────────
    # Dense level bundles rather than markers on a curve; traced as a scatter
    # and flagged low confidence.
    "fig2a": Panel(name="fig2a", page=3, bbox=(0.5985, 0.0921, 0.8402, 0.2250),
                   xlim=(0.0, 1.0), ylim=(0.0, 50.0),
                   description="FIG. 2(a) SU(3) B=0, 2K=10", B=0, N=3),
    "fig2b": Panel(name="fig2b", page=3, bbox=(0.5985, 0.2189, 0.8402, 0.3511),
                   xlim=(0.0, 1.0), ylim=(0.0, 60.0),
                   description="FIG. 2(b) SU(3) B=1, 2K=13", B=1, N=3),
    "fig2c": Panel(name="fig2c", page=3, bbox=(0.5985, 0.3451, 0.8402, 0.4783),
                   xlim=(0.0, 1.0), ylim=(0.0, 100.0),
                   description="FIG. 2(c) SU(3) B=2, 2K=22", B=2, N=3),
}


# ──────────────────────────────────────────────────────────────────────────
# Image helpers
# ──────────────────────────────────────────────────────────────────────────

def load_panel_image(panel: Panel, dpi: int = 600, pages_dir: Path | None = None):
    """Return the cropped panel as a boolean 'ink' array (True = dark)."""
    from PIL import Image

    from render_pages import DEFAULT_OUT, PDF, render

    pages_dir = pages_dir or DEFAULT_OUT
    candidates = sorted(Path(pages_dir).glob(f"p{panel.page:02d}_{dpi}dpi*.png"))
    if not candidates:
        if not PDF.exists():
            raise SystemExit(
                f"{PDF} not found; the article is not redistributed with this "
                "repo. See CITATION.md."
            )
        candidates = [render(PDF, panel.page, dpi, Path(pages_dir))]

    img = Image.open(candidates[0]).convert("L")
    w, h = img.size
    l, t, r, b = panel.bbox
    box = (int(l * w), int(t * h), int(r * w), int(b * h))
    crop = np.asarray(img.crop(box))
    return crop < 128, box, (w, h)


def find_axis_frame(ink: np.ndarray):
    """Locate the plot frame: the strongest long horizontal and vertical runs.

    Returns ``(left, right, top, bottom)`` in crop pixel coordinates.
    """
    rows = ink.sum(axis=1)
    cols = ink.sum(axis=0)
    h, w = ink.shape

    # A frame line spans most of the panel.
    row_cand = np.flatnonzero(rows > 0.5 * w)
    col_cand = np.flatnonzero(cols > 0.5 * h)
    if row_cand.size < 1 or col_cand.size < 1:
        raise ValueError("could not find a plot frame; check the bbox")

    top, bottom = int(row_cand.min()), int(row_cand.max())
    left, right = int(col_cand.min()), int(col_cand.max())
    if bottom - top < 10 or right - left < 10:
        raise ValueError("frame too small; check the bbox")
    return left, right, top, bottom


def find_ticks(ink: np.ndarray, frame, axis: str, depth=26, frac=0.6):
    """Major-tick pixel positions, read off the frame edge that has no data on it.

    Ticks point inward.  Minor ticks are short and major ticks long, so a band
    ``depth`` pixels deep selects only the majors.  For x we read the **top**
    frame: curves live in the lower part of a panel, so the top edge is clean,
    whereas a band above the bottom axis is contaminated by the data itself.
    """
    left, right, top, bottom = frame
    if axis == "x":
        band = ink[top + 4:top + 4 + depth, left:right + 1]
        profile = band.sum(axis=0)
        origin = left
    else:
        band = ink[top:bottom + 1, right - 3 - depth:right - 3]
        profile = band.sum(axis=1)
        origin = top

    hot = profile >= depth * frac
    clusters, run = [], []
    for i, on in enumerate(hot):
        if on:
            run.append(i)
        elif run:
            clusters.append(run)
            run = []
    if run:
        clusters.append(run)
    return [origin + float(np.mean(c)) for c in clusters]


def fit_axis(pixels, values):
    """Least-squares linear pixel -> value map."""
    pixels = np.asarray(pixels, dtype=float)
    values = np.asarray(values, dtype=float)
    A = np.vstack([pixels, np.ones_like(pixels)]).T
    (slope, intercept), *_ = np.linalg.lstsq(A, values, rcond=None)
    predicted = slope * pixels + intercept
    residual = float(np.max(np.abs(predicted - values))) if pixels.size else np.inf
    return float(slope), float(intercept), residual


def detect_markers(ink: np.ndarray, frame, stroke=5, min_area=150):
    """Find plot markers, in two channels.

    Markers sit *on* continuous curves, so plain connected components merge
    everything into one blob.  Instead:

    * **filled** discs survive a morphological opening with a disc of radius
      ``stroke``, while the curves (about ``stroke`` px wide) do not;
    * **open** circles are recovered as enclosed holes -- ``fill_holes(ink) ^
      ink`` -- which is also why we must *not* fill holes before the opening:
      a curve closed off by the axis encloses a large region that would fill
      into a spurious blob.
    """
    from scipy import ndimage

    left, right, top, bottom = frame
    inner = ink[top + 4:bottom - 3, left + 4:right - 3]

    def blobs(mask, lo, hi, filled):
        lab, n = ndimage.label(mask)
        out = []
        for i in range(1, n + 1):
            ys, xs = np.where(lab == i)
            area = ys.size
            if not (lo <= area <= hi):
                continue
            hh = ys.max() - ys.min() + 1
            ww = xs.max() - xs.min() + 1
            if hh == 0 or ww == 0 or not (0.5 < ww / hh < 2.0):
                continue
            if area / (hh * ww) < 0.5:            # roughly round and solid
                continue
            out.append(dict(px=float(xs.mean()) + left + 4,
                            py=float(ys.mean()) + top + 4,
                            area=int(area), filled=filled))
        return out

    r = stroke
    yy, xx = np.ogrid[-r:r + 1, -r:r + 1]
    disc = (xx * xx + yy * yy) <= r * r

    solid = blobs(ndimage.binary_opening(inner, structure=disc), min_area, 1500, True)
    holes = blobs(ndimage.binary_fill_holes(inner) ^ inner, 60, 1200, False)

    markers = solid + holes
    markers.sort(key=lambda m: m["px"])
    return markers


# ──────────────────────────────────────────────────────────────────────────
# Driver
# ──────────────────────────────────────────────────────────────────────────

def digitize(panel: Panel, dpi=600, pages_dir=None):
    """Trace one panel.  Returns ``(records, provenance)``."""
    ink, box, page_size = load_panel_image(panel, dpi=dpi, pages_dir=pages_dir)
    frame = find_axis_frame(ink)

    left, right, top, bottom = frame

    # Frame-edge calibration: deterministic, and independent of whether the
    # tick finder happened to resolve every mark.
    xs, xi, xres = fit_axis([left, right], list(panel.xlim))
    ys, yi, yres = fit_axis([bottom, top], list(panel.ylim))

    notes = []

    # Cross-check: the detected major ticks should land on their stated values.
    xt = find_ticks(ink, frame, "x")
    if panel.xticks and len(xt) == len(panel.xticks):
        pred = [xs * t + xi for t in xt]
        err = max(abs(a - b) for a, b in zip(pred, panel.xticks))
        notes.append(f"x ticks agree to {err:.4f}")
        if err > 0.02:
            notes.append("WARNING: x tick mismatch, check xlim")
    else:
        notes.append(f"x ticks: found {len(xt)}, expected {len(panel.xticks)} "
                     "(calibration still from frame edges)")

    markers = detect_markers(ink, frame)
    records = []
    for m in markers:
        records.append(dict(x=xs * m["px"] + xi, y=ys * m["py"] + yi,
                            filled=m["filled"], area=m["area"],
                            px=m["px"], py=m["py"]))

    provenance = dict(
        panel=panel.name, description=panel.description,
        source="Phys. Rev. D 41, 3814 (1990)",
        pdf_page=panel.page, dpi=dpi, page_size_px=list(page_size),
        crop_box_px=list(box), bbox_fraction=list(panel.bbox),
        frame_px=list(frame),
        xlim=list(panel.xlim), ylim=list(panel.ylim),
        xtick_values=list(panel.xticks), xtick_pixels=list(xt),
        x_fit=dict(slope=xs, intercept=xi, max_residual=xres),
        y_fit=dict(slope=ys, intercept=yi, max_residual=yres),
        n_markers=len(records), notes=notes,
        method="tools/digitize.py connected-component marker tracing",
    )
    return records, provenance


def snap_to_lattice(records, K, atol=0.02):
    """Keep only markers that sit on the ``x = k/K`` lattice, and snap x to it.

    A trace of a plot whose curves cross picks up blobs at the intersections.
    They are recognisable because they do not land on an allowed momentum site
    -- momenta are odd integers, so every genuine marker has ``x = k/K``.
    Dropping the rest removes exactly the artifacts without any hand-tuning.
    """
    grid = np.arange(1, K, 2) / float(K)
    kept, dropped = [], []
    for r in records:
        j = int(np.argmin(np.abs(grid - r["x"])))
        if abs(grid[j] - r["x"]) <= atol:
            out = dict(r)
            out["x_raw"] = r["x"]
            out["x"] = float(grid[j])
            out["k"] = int(2 * j + 1)
            kept.append(out)
        else:
            dropped.append(r)
    return kept, dropped


def check(panel: Panel, records, provenance):
    """Validate a trace by recovering K from the marker x-positions."""
    from dlcq.units import infer_K_from_x_grid

    xs = np.array([r["x"] for r in records])
    if xs.size < 2:
        return {"ok": False, "reason": f"only {xs.size} markers found"}

    K, diag = infer_K_from_x_grid(xs)
    result = {"K_inferred": K, **diag}
    if panel.expected_K is not None:
        result["K_expected"] = panel.expected_K
        result["ok"] = (K == panel.expected_K)
    else:
        result["ok"] = bool(diag["consistent"])
        result["note"] = "paper does not state K for this panel; inferred"
    return result


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", help="panel key (see --list)")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--check", action="store_true",
                    help="recover K from the trace and compare with the paper")
    ap.add_argument("--out", type=Path, help="directory for the CSV + provenance")
    ap.add_argument("--pages-dir", type=Path, default=None)
    args = ap.parse_args(argv)

    if args.list or not args.panel:
        for key, p in PANELS.items():
            k = p.expected_K if p.expected_K else "not stated"
            print(f"  {key:8s} page {p.page}  2K={k:>10}  {p.description}")
        return 0

    panel = PANELS[args.panel]
    records, provenance = digitize(panel, dpi=args.dpi, pages_dir=args.pages_dir)

    print(f"{panel.name}: {len(records)} markers")
    print(f"  frame px      : {provenance['frame_px']}")
    print(f"  x fit residual: {provenance['x_fit']['max_residual']:.4g}")
    print(f"  y fit residual: {provenance['y_fit']['max_residual']:.4g}")
    for note in provenance["notes"]:
        print(f"  note: {note}")

    verdict = check(panel, records, provenance)
    provenance["check"] = verdict
    if args.check:
        print(f"  K inferred    : {verdict['K_inferred']} "
              f"(expected {verdict.get('K_expected', 'not stated')}) "
              f"inliers {verdict['inlier_fraction']:.2f}")

    K = verdict["K_inferred"]
    records, dropped = snap_to_lattice(records, K)

    # Sum-rule calibration of the vertical scale.  int q dx is exactly the
    # quark number (1 for a meson, N*B for a baryon), so the valence series
    # fixes its own normalization -- no tick labels needed, and non-circular
    # because it uses only that identity.  This catches the trap that the frame
    # top is often NOT the last labelled tick: reading Fig 6's top as 21
    # instead of 14.77 inflates every y by 1.42x, and Fig 3(b) needs 1.20x.
    if panel.quark_number and records:
        valence = [r for r in records if r["filled"]] or records
        total = sum(r["y"] for r in valence) * (2.0 / K)
        if total > 1e-9:
            scale = panel.quark_number / total
            for r in records:
                r["y_raw"] = r["y"]
                r["y"] *= scale
            provenance["sum_rule_scale"] = scale
            provenance["sum_rule_quark_number"] = panel.quark_number
            print(f"  sum rule      : rescaled y by {scale:.4f} so "
                  f"int q dx = {panel.quark_number:g}")
    provenance["n_on_lattice"] = len(records)
    provenance["n_dropped_off_lattice"] = len(dropped)
    print(f"  on lattice    : {len(records)} kept, {len(dropped)} dropped "
          f"(curve-crossing artifacts)")

    if args.out:
        args.out.mkdir(parents=True, exist_ok=True)
        csv_path = args.out / f"{panel.name}.csv"
        with open(csv_path, "w") as fh:
            fh.write(f"# {panel.description}\n")
            fh.write(f"# digitized from {provenance['source']}, "
                     f"PDF page {panel.page} at {args.dpi} dpi\n")
            fh.write("# see the companion .json for full calibration provenance\n")
            fh.write("x,y,filled,k,x_raw\n")
            for r in records:
                fh.write(f"{r['x']:.6f},{r['y']:.6f},{int(r['filled'])},"
                         f"{r.get('k','')},{r.get('x_raw', r['x']):.6f}\n")
        json_path = args.out / f"{panel.name}.json"
        json_path.write_text(json.dumps(provenance, indent=2))
        print(f"  wrote {csv_path}")
        print(f"  wrote {json_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
