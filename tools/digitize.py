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
    # Which document the panel is traced from.  The thesis (SLAC-333) reprints
    # the article's figures at better print quality, so it is the preferred
    # target where a panel appears in both.
    source: str = "article"


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

    # ── THESIS Figs. 11 and 12 -- the same panels as the article's Figs. 5
    # and 6, but far more legibly printed.  These are the preferred targets.
    "t12a": Panel(name="t12a", page=82, source="thesis",
                  bbox=(0.330, 0.089, 0.690, 0.300),
                  xlim=(0.0, 1.0), ylim=(0.0, 15.0),
                  xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                  description="THESIS Fig 12(a) = article Fig 6(a): 1st baryon, 2K=21",
                  expected_K=21, B=1, N=3, quark_number=3.0, sector="valence"),
    "t12b": Panel(name="t12b", page=82, source="thesis",
                  bbox=(0.330, 0.313, 0.690, 0.523),
                  xlim=(0.0, 1.0), ylim=(0.0, 15.0),
                  xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                  description="THESIS Fig 12(b) = article Fig 6(b): 2nd baryon, 2K=21",
                  expected_K=21, B=1, N=3, quark_number=3.0, sector="valence"),
    "t12c": Panel(name="t12c", page=82, source="thesis",
                  bbox=(0.330, 0.540, 0.690, 0.750),
                  xlim=(0.0, 1.0), ylim=(0.0, 15.0),
                  xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                  description="THESIS Fig 12(c) = article Fig 6(c): 3rd baryon, 2K=21",
                  expected_K=21, B=1, N=3, quark_number=3.0, sector="valence"),

    # THESIS Fig 4: the article's Fig. 2, and it states the K on each panel
    # ("K = 10/2", "13/2", "22/2"), confirming 2K = 10, 13, 22.
    "t4a": Panel(name="t4a", page=64, source="thesis",
                 bbox=(0.325, 0.083, 0.690, 0.240),
                 xlim=(0.0, 1.0), ylim=(0.0, 55.0),
                 description="THESIS Fig 4(a) = article Fig 2(a): SU(3) B=0, 2K=10",
                 B=0, N=3, sector="spectrum"),
    "t4b": Panel(name="t4b", page=64, source="thesis",
                 bbox=(0.325, 0.243, 0.690, 0.395),
                 xlim=(0.0, 1.0), ylim=(0.0, 70.0),
                 description="THESIS Fig 4(b) = article Fig 2(b): SU(3) B=1, 2K=13",
                 B=1, N=3, sector="spectrum"),
    "t4c": Panel(name="t4c", page=64, source="thesis",
                 bbox=(0.325, 0.398, 0.690, 0.556),
                 xlim=(0.0, 1.0), ylim=(0.0, 125.0),
                 description="THESIS Fig 4(c) = article Fig 2(c): SU(3) B=2, 2K=22",
                 B=2, N=3, sector="spectrum"),

    # THESIS Fig 18(a): the five-quark contribution to the lightest N=3 baryon
    # structure function, shown ALONE with q and qbar separated, at 2K=15.
    # This is the article's Fig. 4(b), and it is a far better target than
    # Fig. 12(a)/Fig. 6(a): those overlay the valence curve on the higher-Fock
    # one and the two cross three times, so point-to-series assignment is
    # unreliable -- and at scan resolution the dashes of the dashed curve read
    # as markers.  Here there are no crossings.
    "t18a": Panel(name="t18a", page=89, source="thesis",
                  bbox=(0.220, 0.078, 0.780, 0.405),
                  xlim=(0.0, 1.0), ylim=(0.0, 25.0),
                  xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                  description="THESIS Fig 18(a) = article Fig 4(b): five-quark "
                              "contribution to the lightest N=3 baryon, 2K=15",
                  expected_K=15, B=1, N=3, sector="higher-fock"),

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

    from render_pages import DEFAULT_OUT, SOURCES, render

    src = SOURCES.get(panel.source)
    prefix = "p" if panel.source == "article" else "t"
    pages_dir = pages_dir or DEFAULT_OUT
    candidates = sorted(Path(pages_dir).glob(
        f"{prefix}{panel.page:02d}_{dpi}dpi*.png"))
    if not candidates:
        if src is None or not src.exists():
            raise SystemExit(
                f"{src} not found; source documents are not redistributed with "
                "this repo. See CITATION.md."
            )
        candidates = [render(src, panel.page, dpi, Path(pages_dir),
                             prefix=prefix)]

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

def digitize(panel: Panel, dpi=600, pages_dir=None, use_lattice_probe=True):
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

    # Two stages.  Blob detection is general enough to establish K on a panel
    # whose K the paper never states; the column probe then re-reads the plot
    # at that K, which is far more robust where curves cross or markers touch.
    markers = detect_markers(ink, frame)
    markers, legend_hits = suppress_legend(markers, frame)
    if legend_hits:
        notes.append(f"suppressed {len(legend_hits)} in-plot legend/text blobs")
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
        n_markers=len(records), n_legend_suppressed=len(legend_hits),
        notes=notes,
        method="tools/digitize.py connected-component marker tracing",
    )
    # Second pass.  Blob detection is general enough to establish K even when
    # the paper never states it; the column probe then re-reads the plot at that
    # K, which is far more robust where curves cross or markers touch.
    if use_lattice_probe:
        K_probe = panel.expected_K
        if not K_probe and records:
            try:
                from dlcq.units import infer_K_from_x_grid
                K_probe, _ = infer_K_from_x_grid([r["x"] for r in records])
            except Exception:
                K_probe = None
        if K_probe:
            probed = trace_at_lattice(ink, frame, K_probe)
            probed, probe_legend = suppress_legend(probed, frame)
            # The probe samples only lattice columns, so it recovers too few
            # legend glyphs per row for the row rule to fire.  Reuse the box
            # the blob pass established.
            probed, boxed = drop_in_box(probed, frame,
                                        legend_bbox(legend_hits, frame))
            probe_legend = list(probe_legend) + list(boxed)
            # Prefer the probe whenever it found a plausible number of
            # markers.  Requiring it to match the blob count was backwards:
            # the blob count is inflated by legend text, so a clean probe was
            # being rejected for finding *fewer* things.
            if len(probed) >= 3:
                provenance["probe_records"] = [
                    dict(x=xs * m["px"] + xi, y=ys * m["py"] + yi,
                         filled=m["filled"], area=m["area"],
                         px=m["px"], py=m["py"], k=m.get("k"))
                    for m in probed]
                provenance["detector"] = "lattice column probe"
                provenance["K_probe"] = K_probe
                provenance["n_legend_suppressed"] += len(probe_legend)
                notes.append(
                    f"column probe at 2K={K_probe}: {len(probed)} markers "
                    f"({sum(m['filled'] for m in probed)} filled, "
                    f"{sum(not m['filled'] for m in probed)} open)")
            else:
                provenance["detector"] = "blob (probe found too few)"
        else:
            provenance["detector"] = "blob (no K)"
    else:
        provenance["detector"] = "blob"

    return records, provenance


def trace_at_lattice(ink, frame, K, stroke=5, marker_min=14, marker_max=46,
                     pad=3):
    """Probe the columns where markers *must* be, instead of hunting blobs.

    Momenta are odd integers, so a marker can only sit at ``x = k/K``.  That
    turns 2D detection into a handful of 1D problems and fixes the failure mode
    of blob detection on these plots: where curves cross or markers touch,
    connected components merge and morphology loses them.  Here a merged blob
    is irrelevant -- we only ask what lies in a known column.

    The subtlety is open circles.  They are drawn with a white interior that
    masks the connecting curve, so a vertical slice through one hits the upper
    arc, a gap, then the lower arc: **two** runs, not one.  A filled disc gives
    a single solid run.  So we group runs whose combined span is marker-sized
    and read the grouping itself:

        one run, solid            -> filled marker
        two runs around a gap     -> open marker

    which also gives the fill classification for free, rather than by
    thresholding an interior density that overlaps between the two classes.

    Returns marker dicts with the same keys as :func:`detect_markers`.
    """
    left, right, top, bottom = frame
    fw = float(right - left)
    out = []

    for k in range(1, K, 2):
        xc = left + fw * (k / float(K))
        c0, c1 = int(round(xc - pad)), int(round(xc + pad)) + 1
        c0, c1 = max(c0, left + 1), min(c1, right)
        if c1 <= c0:
            continue

        strip = ink[top + 2:bottom - 1, c0:c1]
        if strip.size == 0:
            continue
        # A row counts as ink if most of the strip's width is dark, which keeps
        # a steeply climbing curve from masquerading as a tall marker.
        col = strip.mean(axis=1) > 0.5

        runs, start = [], None
        for i, v in enumerate(col):
            if v and start is None:
                start = i
            elif not v and start is not None:
                runs.append([start, i - 1])
                start = None
        if start is not None:
            runs.append([start, len(col) - 1])

        # Group runs that together span a marker: a ring is two arcs about a
        # hollow centre, so allow one interior gap.
        used = [False] * len(runs)
        for i, (a, b) in enumerate(runs):
            if used[i]:
                continue
            span_a, span_b, members = a, b, [i]
            for j in range(i + 1, len(runs)):
                if used[j]:
                    continue
                gap = runs[j][0] - span_b - 1
                new_span = runs[j][1] - span_a + 1
                if gap <= 0 or new_span > marker_max:
                    break
                # The hollow of a ring is comparable to its arcs; a gap far
                # larger than that belongs to a different series.
                if gap > 0.8 * new_span:
                    break
                span_b = runs[j][1]
                members.append(j)
                break                      # at most one interior gap
            height = span_b - span_a + 1
            if not (marker_min <= height <= marker_max):
                continue

            # Height alone is not enough: a steeply climbing curve rises as far
            # across the strip as a marker is tall.  A marker is also WIDE --
            # round, roughly as broad as it is high -- while a curve stroke is
            # only ``stroke`` px across.  Measure the horizontal ink span at the
            # run's mid-height and require it to be marker-like.
            row = top + 2 + (span_a + span_b) // 2
            lo = max(int(xc - marker_max), left + 1)
            hi = min(int(xc + marker_max) + 1, right)
            band = ink[row, lo:hi]
            idx = np.flatnonzero(band)
            width = (idx[-1] - idx[0] + 1) if idx.size else 0
            if width < 0.55 * height:
                continue

            for m in members:
                used[m] = True
            is_open = len(members) > 1
            out.append(dict(px=xc, py=top + 2 + 0.5 * (span_a + span_b),
                            area=int(strip[span_a:span_b + 1].sum()),
                            height=int(height), filled=not is_open, k=k))
    return out


def legend_bbox(dropped, frame, margin=0.02):
    """Bounding box (in frame fractions) of the blobs :func:`suppress_legend` cut.

    The blob pass sees the whole legend -- its markers *and* the letter counters
    of its text -- so it localizes the region reliably.  The column probe only
    samples lattice columns and therefore finds too few legend glyphs per row
    for the >=3-in-a-row rule to fire, so it is given this box instead.
    """
    if not dropped:
        return None
    left, right, top, bottom = frame
    fw, fh = float(right - left), float(bottom - top)
    xs = [(m["px"] - left) / fw for m in dropped]
    ys = [(bottom - m["py"]) / fh for m in dropped]
    return (min(xs) - margin, min(ys) - margin,
            max(xs) + margin, max(ys) + margin)


def drop_in_box(markers, frame, box):
    """Remove markers inside a frame-fraction bounding box."""
    if box is None:
        return markers, []
    left, right, top, bottom = frame
    fw, fh = float(right - left), float(bottom - top)
    x0, y0, x1, y1 = box
    keep, cut = [], []
    for m in markers:
        fx = (m["px"] - left) / fw
        fy = (bottom - m["py"]) / fh
        (cut if (x0 <= fx <= x1 and y0 <= fy <= y1) else keep).append(m)
    return keep, cut


def trace_spectrum(ink, frame, x_fracs, min_gap=3):
    """Read a dense level-trajectory plot column by column.

    Figs. 2 / thesis Fig. 4 show the whole eigenvalue spectrum as ~25-40
    trajectories against the coupling.  Individual curves cannot be followed
    through the crossings, and there are no markers to find -- but the thing we
    actually want is the *set* of eigenvalues at a given coupling, which is
    exactly what one column contains.

    So we scan the requested columns and return the centre of every ink run.
    Runs closer than ``min_gap`` are merged, since two trajectories that touch
    at this resolution cannot be separated.  That merging is the honest
    limitation: a traced level is a lower bound on the number of levels there.

    Returns ``{x_frac: [y_frac_from_bottom, ...]}`` with y ascending.
    """
    left, right, top, bottom = frame
    fw, fh = float(right - left), float(bottom - top)
    out = {}

    for xf in x_fracs:
        col = int(round(left + xf * fw))
        if not (left < col < right):
            continue
        # Stay clear of the frame lines themselves, and take ANY ink in the
        # column: a trajectory is a thin stroke, so requiring most of a
        # multi-pixel strip to be dark loses it.
        pad = max(4, int(0.012 * fh))
        strip = ink[top + pad:bottom - pad, max(col - 1, left + 1):min(col + 2, right)]
        dark = strip.any(axis=1)

        runs, start = [], None
        for i, v in enumerate(dark):
            if v and start is None:
                start = i
            elif not v and start is not None:
                runs.append((start, i - 1))
                start = None
        if start is not None:
            runs.append((start, len(dark) - 1))

        merged = []
        for a, b in runs:
            if merged and a - merged[-1][1] <= min_gap:
                merged[-1] = (merged[-1][0], b)
            else:
                merged.append((a, b))

        ys = [1.0 - (0.5 * (a + b) + pad) / fh for a, b in merged]
        out[xf] = sorted(y for y in ys if 0.0 <= y <= 1.0)
    return out


def suppress_legend(markers, frame, y_tol=0.015, min_row=3, x_span=0.15,
                    y_min=0.55):
    """Drop markers belonging to an in-plot legend.

    These panels put their legends *inside* the axes, and two things there get
    mistaken for data:

    * the legend's own sample markers;
    * the enclosed counters of letters like q, b, o, d, which the hole channel
      sees as open circles.

    Neither is separable by size -- measured on Fig. 6(a), legend blobs span
    124-416 px and data blobs 164-458 px, overlapping.  Nor by position alone,
    since the valence peak reaches 79% of the frame height.

    What does separate them is *shape of the arrangement*: a legend is a
    horizontal row.  Three or more markers sharing a y to within ``y_tol`` of
    the frame height and spanning ``x_span`` of its width cannot be a peaked
    structure function, which is single-valued and falls away from its maximum.
    Restricting to the upper part of the frame (``y_min``) keeps a genuine flat
    tail near zero safe.

    Returns ``(kept, dropped)``.
    """
    left, right, top, bottom = frame
    fw, fh = float(right - left), float(bottom - top)
    if fw <= 0 or fh <= 0 or not markers:
        return markers, []

    enriched = [(m, (m["px"] - left) / fw, (bottom - m["py"]) / fh)
                for m in markers]
    dropped_ids = set()

    for i, (_, _, fy) in enumerate(enriched):
        if fy < y_min:
            continue
        row = [j for j, (_, _, fy2) in enumerate(enriched)
               if abs(fy2 - fy) <= y_tol and enriched[j][2] >= y_min]
        if len(row) < min_row:
            continue
        xs = [enriched[j][1] for j in row]
        if max(xs) - min(xs) >= x_span:
            dropped_ids.update(row)

    kept = [m for i, (m, _, _) in enumerate(enriched) if i not in dropped_ids]
    dropped = [m for i, (m, _, _) in enumerate(enriched) if i in dropped_ids]
    return kept, dropped


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

    # Snap to the K the paper STATES when it states one.  The inference exists
    # to validate the trace and to settle Figs. 3 and 4, where K is unstated --
    # it must not override a published fact.  With few surviving points a
    # coarser lattice can fit as well as the true one (fig6b/c returned 20
    # instead of 21 after legend suppression thinned them), and silently
    # adopting that would corrupt every x.
    K = panel.expected_K if panel.expected_K else verdict["K_inferred"]
    provenance["K_used"] = K
    provenance["K_source"] = "stated in paper" if panel.expected_K else "inferred"

    # digitize() already re-read the plot with the column probe if it could;
    # this just records which detector produced the records we are snapping.
    if provenance.get("probe_records") is not None:
        records = provenance.pop("probe_records")

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
