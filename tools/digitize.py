#!/usr/bin/env python3
"""Trace curve data out of the scanned figures.

The paper carries no vector content (see ``tools/render_pages.py``), so the only
way to per-point numbers is pixel tracing.  The scan is 1-bit at ~612 ppi, which
helps: axis frames are solid runs and ticks are short perpendicular marks.

Detection
---------
Momenta are odd integers, so a marker on a structure-function panel can only sit
at ``x = k/K``.  Probing those columns turns 2D detection into a few 1D problems.
In each column three concentric templates are scored against the ink: an annulus
on the marker's stroke, which is inked for every series and says a marker is
*there*; and a disc and a core inside it, which say *which* series, since the
paper draws both with the same circle and distinguishes them only by whether the
middle is filled.  Open triangles, where a panel has them, are found by the
shape of their hole.

Two things make this work where reading vertical runs did not.  Frame lines are
excluded from the template averages rather than counted as ink, so a marker
sitting *on* the x axis -- which is every marker in the large-x tail -- still
reads as a full disc or a full ring.  And the two classes are searched
independently, so where the two series overlap, as they do throughout that
tail, a ring resting on a disc is two detections instead of one merged blob.

Calibration is semi-automatic and auditable.  You supply the panel's bounding
box, the values of its x ticks, and the values of its y ticks by index; the code
locates the tick *pixels*, snaps them to a uniform ladder and fits the mapping.
Everything it used is written to a provenance JSON beside the CSV, so any number
can be traced back to the pixels it came from.

Self-validation
---------------
Three checks, none of which involves a solver:

* Figs. 5 and 6(a-c) state their K.  ``dlcq.units.infer_K_from_x_grid`` recovers
  it from the marker positions alone, which is what licenses the same procedure
  on Figs. 3, 4 and 6(d), where K is never stated.
* ``int q dx`` is the quark number exactly, so a valence series that integrates
  to something else has a vertical scale that is wrong.
* ``tools/trace_overlay.py`` draws the detections back onto the scan, which is
  what makes a *missing* marker visible.

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
    # The y axis is calibrated from its own tick marks when these are given:
    # pairs of (index of the tick counting up from the bottom, printed value).
    # Reading ``ylim`` off the topmost printed number is a trap -- in Figs.
    # 3(b), 4(b), 4(c) and 6 that number sits well below the frame edge, and
    # taking it for the edge shrinks every value by 13-19%.
    ytick_anchor: tuple = ()
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
    # Marker shapes the panel plots, one per series.  Used to report which
    # momentum sites came back short, so a series that was never found is
    # visible instead of silently absent.
    shapes: tuple = ("filled", "open")


PANELS = {
    # ── FIG. 3 (page 4): valence structure functions, K NOT stated ──────────
    "fig3a": Panel(name="fig3a", page=4, bbox=(0.1149, 0.7900, 0.2765, 0.8813),
                   xlim=(0.0, 1.0), ylim=(0.0, 3.5),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   ytick_anchor=((0, 0.0), (2, 0.7), (4, 1.4), (6, 2.1),
                                 (8, 2.8), (10, 3.5)),
                   description="FIG. 3(a) SU(3) meson valence at m/g=1.6 "
                               "(filled) and 0.1 (open); K not stated",
                   B=0, N=3, quark_number=1.0, sector="valence"),
    "fig3b": Panel(name="fig3b", page=4, bbox=(0.3061, 0.7906, 0.4677, 0.8818),
                   xlim=(0.0, 1.0), ylim=(0.0, 13.5),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   ytick_anchor=((0, 0.0), (2, 3.75), (4, 7.5), (6, 11.25)),
                   description="FIG. 3(b) SU(3) baryon valence at m/g=1.6 "
                               "(filled) and 0.1 (open); K not stated",
                   B=1, N=3, quark_number=3.0, sector="valence"),

    # ── FIG. 4 (page 4): higher-Fock sectors, scaled by 10^n ────────────────
    "fig4a": Panel(name="fig4a", page=4, bbox=(0.5376, 0.0949, 0.6999, 0.1902),
                   xlim=(0.0, 1.0), ylim=(0.0, 21.0),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   ytick_anchor=((0, 0.0), (2, 7.0), (4, 14.0), (6, 21.0)),
                   description="FIG. 4(a) meson, one extra qqbar pair, at "
                               "m/g=1.6 (filled) and 0.1 (open), both x10^4",
                   B=0, N=3, sector="higher-fock"),
    "fig4b": Panel(name="fig4b", page=4, bbox=(0.7304, 0.0961, 0.8929, 0.1911),
                   xlim=(0.0, 1.0), ylim=(0.0, 26.6),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   ytick_anchor=((0, 0.0), (2, 7.5), (4, 15.0), (6, 22.5)),
                   description="FIG. 4(b) baryon, one extra pair: q(x) at "
                               "m/g=1.6 (filled, x10^3) and 0.1 (open, x10^2), "
                               "qbar(x) at 1.6 (down triangle) and 0.1 (up)",
                   B=1, N=3, sector="higher-fock",
                   shapes=("filled", "open", "triangle", "triangle_down")),
    "fig4c": Panel(name="fig4c", page=4, bbox=(0.5372, 0.1811, 0.6994, 0.2765),
                   xlim=(0.0, 1.0), ylim=(0.0, 12.7),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   ytick_anchor=((0, 0.0), (2, 3.75), (4, 7.5), (6, 11.25)),
                   description="FIG. 4(c) baryon, two extra pairs, at m/g=1.6 "
                               "(filled, x10^4) and 0.1 (open, x10^7)",
                   B=1, N=3, sector="higher-fock"),

    # ── FIG. 5 (page 4): meson spectrum, 2K = 24 ───────────────────────────
    "fig5a": Panel(name="fig5a", page=4, bbox=(0.5478, 0.7005, 0.7149, 0.7973),
                   xlim=(0.0, 1.0), ylim=(0.0, 3.62),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   ytick_anchor=((0, 0.0), (2, 1.2), (4, 2.4), (6, 3.6)),
                   description="FIG. 5(a) 1st meson state, m/g=1.6, 2K=24",
                   expected_K=24, B=0, N=3, quark_number=1.0, sector="valence"),
    "fig5b": Panel(name="fig5b", page=4, bbox=(0.7037, 0.7005, 0.8708, 0.7973),
                   # (b) carries no y labels of its own.  Its tick ladder is
                   # (a)'s to a pixel, so it is on (a)'s 0..3.6 scale, not the
                   # 0..12 of the right-hand axis -- which belongs to (d).  The
                   # quark-number sum rule then lands on 1.00; reading it as
                   # 0..12 misses by the ratio 12/3.6 = 3.33.
                   xlim=(0.0, 1.0), ylim=(0.0, 3.62),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   ytick_anchor=((0, 0.0), (2, 1.2), (4, 2.4), (6, 3.6)),
                   description="FIG. 5(b) 2nd meson state, m/g=1.6, 2K=24",
                   expected_K=24, B=0, N=3, quark_number=1.0, sector="valence"),
    "fig5c": Panel(name="fig5c", page=4, bbox=(0.5478, 0.7873, 0.7106, 0.8836),
                   xlim=(0.0, 1.0), ylim=(0.0, 3.62),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   ytick_anchor=((0, 0.0), (2, 1.2), (4, 2.4)),
                   description="FIG. 5(c) 3rd meson state, m/g=1.6, 2K=24",
                   expected_K=24, B=0, N=3, quark_number=1.0, sector="valence"),
    "fig5d": Panel(name="fig5d", page=4, bbox=(0.6994, 0.7873, 0.8665, 0.8836),
                   # (d) owns the right-hand 0..12 axis.  No sum rule here:
                   # this is the one panel that puts its x10^n on the *filled*
                   # series, so the valence curve as printed integrates to
                   # 10^4 x P(2 partons), not to the quark number.
                   xlim=(0.0, 1.0), ylim=(0.0, 12.06),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   ytick_anchor=((0, 0.0), (2, 2.4), (4, 4.8), (6, 7.2),
                                 (8, 9.6), (10, 12.0)),
                   description="FIG. 5(d) 11th meson state, m/g=1.6, 2K=24",
                   expected_K=24, B=0, N=3, sector="valence"),

    # ── FIG. 6 (page 5): baryon spectrum, 2K = 21 ──────────────────────────
    # Only 0 and 10.5 are printed, and 10.5 is the 5th tick from the bottom,
    # not the frame top: the ladder runs 0.003/0.185/0.365/0.540/0.720/0.895/
    # 0.997 of the frame height, so the top is 14.6.  (b) and (c) carry no
    # labels but their ladders match (a)'s to a pixel.  Reading the top as
    # 10.5 would shrink every value by 30%, and as 21 inflate it by 44%.
    "fig6a": Panel(name="fig6a", page=5, bbox=(0.1257, 0.1108, 0.2853, 0.2030),
                   xlim=(0.0, 1.0), ylim=(0.0, 14.6),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   ytick_anchor=((0, 0.0), (4, 10.5)),
                   description="FIG. 6(a) 1st baryon state, m/g=1.6, 2K=21",
                   expected_K=21, B=1, N=3, quark_number=3.0, sector="valence"),
    "fig6b": Panel(name="fig6b", page=5, bbox=(0.2737, 0.1108, 0.4337, 0.2030),
                   xlim=(0.0, 1.0), ylim=(0.0, 14.6),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   ytick_anchor=((0, 0.0), (4, 10.5)),
                   description="FIG. 6(b) 2nd baryon state, m/g=1.6, 2K=21",
                   expected_K=21, B=1, N=3, quark_number=3.0, sector="valence"),
    "fig6c": Panel(name="fig6c", page=5, bbox=(0.1257, 0.2141, 0.2853, 0.3059),
                   xlim=(0.0, 1.0), ylim=(0.0, 14.6),
                   xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                   ytick_anchor=((0, 0.0), (4, 10.5)),
                   description="FIG. 6(c) 3rd baryon state, m/g=1.6, 2K=21",
                   expected_K=21, B=1, N=3, quark_number=3.0, sector="valence"),
    # The caption gives 2K = 21 for (a)-(c) and says nothing at all about (d),
    # so its K is inferred -- see docs/inferred-K.md.  It is 24, not the 22 of
    # Fig. 2(c), on three independent counts: the markers sit on x = k/24 to
    # 2 px where k/22 misses by 23; the valence series then integrates to
    # 6.03 against the exact quark number 6 (6.48 at K = 22); and its mean
    # momentum is 0.1677 against the exact 1/6 = 0.1667 (0.1810 at K = 22).
    # Only the first of those uses the x calibration at all.
    "fig6d": Panel(name="fig6d", page=5, bbox=(0.2737, 0.2141, 0.4337, 0.3059),
                   xlim=(0.0, 0.6), ylim=(0.0, 47.9),
                   ytick_anchor=((0, 0.0), (2, 12.0), (4, 24.0), (6, 36.0),
                                 (8, 48.0)),
                   description="FIG. 6(d) 1st B=2 state, m/g=1.6, 2K=24 (inferred)",
                   expected_K=24, B=2, N=3, quark_number=6.0, sector="valence",
                   shapes=("filled", "open", "triangle")),

    # ── FIG. 8 (page 7): meson mass vs m/g, with Hamer's SU(2) lattice points ──
    "fig8a": Panel(name="fig8a", page=7, bbox=(0.5943, 0.1381, 0.8324, 0.2760),
                   # 8 x ticks every 0.25 (0..1.75) and 9 y ticks every 0.5
                   # (0..4.0); the frame edges sit just beyond the last of each.
                   xlim=(0.0, 1.756), ylim=(0.0, 4.02),
                   xticks=(0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75),
                   ytick_anchor=((0, 0.0), (2, 1.0), (4, 2.0), (6, 3.0),
                                 (8, 4.0)),
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

def load_panel_image(panel: Panel, dpi: int = 600, pages_dir: Path | None = None,
                     bottom_margin: float = 0.005):
    """Return the cropped panel as a boolean 'ink' array (True = dark).

    The crop is extended ``bottom_margin`` (a fraction of the page height)
    below the stated bbox.  A marker at ``q = 0`` is centred *on* the x axis,
    so half of it lies below the frame; without that margin the detector sees
    a cap instead of a disc and the whole large-x tail is lost.  The margin is
    small enough to stay inside the gap between stacked panels -- ~70 px at
    600 dpi for Fig. 6 -- and ``find_axis_frame`` still sees exactly one frame.
    """
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
    stated_bottom = int(b * h)
    box = (int(l * w), int(t * h), int(r * w),
           min(h, stated_bottom + int(bottom_margin * h)))
    crop = np.asarray(img.crop(box))
    margin_rows = box[3] - stated_bottom
    return crop < 128, box, (w, h), margin_rows


def find_axis_frame(ink: np.ndarray, margin_rows: int = 0):
    """Locate the plot frame: the strongest long horizontal and vertical runs.

    ``margin_rows`` are the extra rows :func:`load_panel_image` appended below
    the stated bbox; they are excluded from the search so that the frame of a
    panel stacked underneath can never be mistaken for this one's x axis.

    Returns ``(left, right, top, bottom)`` in crop pixel coordinates.
    """
    if margin_rows:
        ink = ink[:-margin_rows]
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


def find_tick_ladder(ink, frame, R, depth_frac=2.3, frac=0.6):
    """Rows of the y tick marks on the left frame, bottom first.

    The ticks are evenly spaced by construction, so the detections are snapped
    to the best-fitting uniform ladder and anything off it dropped.  That is
    what makes the reading robust where a curve runs into the left frame and
    adds a spurious "tick", as in Fig. 4(b).
    """
    left, right, top, bottom = frame
    depth = int(round(depth_frac * R))
    best = []
    # Either frame edge carries the same ladder.  Try both and keep the better
    # reading: in Fig. 4(b) four curves pile into the left edge and drown its
    # ticks, while the right edge is clean.
    for band in (ink[top:bottom + 1, left + 3:left + 3 + depth],
                 ink[top:bottom + 1, right - 2 - depth:right - 2]):
        rows = _band_runs(band, top, depth * frac)
        rungs = _uniform_ladder(rows)
        if len(rungs) > len(best):
            best = rungs
    return best


def _band_runs(band, origin, threshold):
    """Row centres of the ink runs in a band along a frame edge."""
    hot = band.sum(axis=1) >= threshold
    runs, cur = [], []
    for i, on in enumerate(hot):
        if on:
            cur.append(i)
        elif cur:
            runs.append(cur)
            cur = []
    if cur:
        runs.append(cur)
    return sorted(origin + float(np.mean(c)) for c in runs)


def _uniform_ladder(rows):
    """Snap tick candidates to the uniform ladder that explains most of them.

    Ticks are evenly spaced by construction, so anything off the grid is not a
    tick: a curve running into the frame, or the frame line itself broadened
    into two or three detections.  Only the candidate closest to each grid
    position is kept, which is what stops those extra detections near the axis
    from shifting every index by one and mis-labelling the whole axis.

    Returns the rungs bottom-first, so index 0 is the tick at the bottom.
    """
    if len(rows) < 3:
        return []
    step = float(np.median(np.diff(rows)))
    if step <= 0:
        return []
    best = {}
    for anchor in rows:
        idx = np.round((np.asarray(rows) - anchor) / step).astype(int)
        keep = {}
        for r, i in zip(rows, idx):
            err = abs(r - (anchor + i * step))
            if err > 0.25 * step:
                continue
            if i not in keep or err < keep[i][1]:
                keep[i] = (r, err)
        if len(keep) > len(best):
            best = keep
    return sorted((r for r, _ in best.values()), reverse=True)


def calibrate_y(ink, frame, panel, R):
    """Pixel -> value map for the y axis, from the ticks when it can.

    Falls back to stretching ``panel.ylim`` between the frame edges, which is
    all that is possible on a panel whose ticks are unlabelled.  Returns
    ``(slope, intercept, residual, source, implied_top)``.
    """
    left, right, top, bottom = frame
    if panel.ytick_anchor:
        ladder = find_tick_ladder(ink, frame, R)
        idx = [i for i, _ in panel.ytick_anchor]
        if ladder and max(idx) < len(ladder):
            pix = [ladder[i] for i in idx]
            vals = [v for _, v in panel.ytick_anchor]
            # Zero is the *centre* of the axis line, not the rung the ladder
            # reports there.  A tick at zero is drawn inward from the line and
            # merges with it, so that run's centre sits above the line's --
            # 3 px on Fig. 8(a), whose axis is 7 px thick.  Interior ticks are
            # isolated and unbiased, so they fix the slope; the axis fixes zero.
            axis = _axis_centre(ink, frame)
            pix = [axis if v == 0.0 else p for p, v in zip(pix, vals)]
            slope, inter, res = fit_axis(pix, vals)
            return slope, inter, res, "y ticks", slope * top + inter
    slope, inter, res = fit_axis([bottom, top], list(panel.ylim))
    return slope, inter, res, "frame edges", panel.ylim[1]


def _axis_centre(ink, frame):
    """Row of the middle of the bottom frame line, which is where y = 0 is."""
    left, right, top, bottom = frame
    rows = np.flatnonzero(ink.sum(axis=1) > 0.5 * ink.shape[1])
    line = rows[rows > (top + bottom) // 2]
    return float(line.mean()) if line.size else float(bottom)


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
    ink, box, page_size, margin = load_panel_image(panel, dpi=dpi,
                                                   pages_dir=pages_dir)
    frame = find_axis_frame(ink, margin_rows=margin)

    left, right, top, bottom = frame

    # x: frame-edge calibration, deterministic and independent of whether the
    # tick finder resolved every mark.  y: from the ticks where their values
    # are known, because on these panels the top printed number is usually not
    # the frame edge.
    R = MARKER_R_600 * dpi / 600.0
    xs, xi, xres = fit_axis([left, right], list(panel.xlim))
    ys, yi, yres, ysrc, ytop = calibrate_y(ink, frame, panel, R)

    notes = []
    notes.append(f"y calibrated from {ysrc}; frame top = {ytop:.3f}")
    if ysrc == "y ticks" and abs(ytop - panel.ylim[1]) > 0.02 * abs(ytop):
        notes.append(f"WARNING: y ticks put the frame top at {ytop:.3f}, "
                     f"but ylim says {panel.ylim[1]:.3f}")

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
                            px=m["px"], py=m["py"],
                            kind="filled" if m["filled"] else "open",
                            conf="blob"))

    provenance = dict(
        panel=panel.name, description=panel.description,
        source="Phys. Rev. D 41, 3814 (1990)",
        pdf_page=panel.page, dpi=dpi, page_size_px=list(page_size),
        crop_box_px=list(box), bbox_fraction=list(panel.bbox),
        frame_px=list(frame),
        xlim=list(panel.xlim), ylim=list(panel.ylim),
        xtick_values=list(panel.xticks), xtick_pixels=list(xt),
        x_fit=dict(slope=xs, intercept=xi, max_residual=xres),
        y_fit=dict(slope=ys, intercept=yi, max_residual=yres,
                   source=ysrc, frame_top=ytop,
                   anchors=[list(a) for a in panel.ytick_anchor]),
        n_markers=len(records), n_legend_suppressed=len(legend_hits),
        notes=notes,
        method="tools/digitize.py connected-component marker tracing",
    )
    # Second pass.  Blob detection is general enough to establish K even when
    # the paper never states it; the column probe then re-reads the plot at that
    # K, which is far more robust where curves cross or markers touch.
    # Only on panels that plot a structure function: the probe's whole premise
    # is that markers sit on x = k/K.  Fig. 2's spectra and Fig. 8's mass curve
    # have no momentum lattice, and probing them invents one.
    if use_lattice_probe and panel.sector:
        K_probe = panel.expected_K
        if not K_probe and records:
            try:
                from dlcq.units import infer_K_from_x_grid
                K_probe, _ = infer_K_from_x_grid([r["x"] for r in records])
            except Exception:
                K_probe = None
        if K_probe:
            labels, big = curve_components(ink, MARKER_R_600 * dpi / 600.0)
            xpix = lattice_columns(K_probe, xs, xi, panel.xlim)
            probed = trace_at_lattice(ink, frame, K_probe, xpix, dpi=dpi,
                                      labels=labels, big=big,
                                      shapes=panel.shapes)
            # No legend box here.  The blob pass needs one because its notion
            # of a marker is loose enough to match a letter's counter, but the
            # probe rejects the legend by connectivity, and the box is a blunt
            # instrument: in Fig. 6(c) it reaches down over the valence peak,
            # which is the tallest marker in the panel.
            if len(probed) >= max(3, len(records) // 2):
                provenance["probe_records"] = [
                    dict(x=xs * m["px"] + xi, y=ys * m["py"] + yi,
                         filled=m["filled"], area=m["area"],
                         px=m["px"], py=m["py"], k=m.get("k"),
                         kind=m.get("kind"), conf=m.get("conf"),
                         ring=m.get("ring"), core=m.get("core"),
                         disc=m.get("disc"))
                    for m in probed]
                provenance["detector"] = "lattice column probe"
                provenance["K_probe"] = K_probe
                found = {(m["k"], m["kind"]) for m in probed}
                missing = [k for k in sorted(xpix)
                           if not {(k, s) for s in panel.shapes} <= found]
                coincident = sorted({m["k"] for m in probed
                                     if m.get("conf") == "coincident"})
                provenance["sites"] = len(xpix)
                provenance["sites_incomplete"] = missing
                provenance["sites_coincident"] = coincident
                by_shape = {s: sum(m.get("kind") == s for m in probed)
                            for s in panel.shapes}
                notes.append(
                    f"column probe at 2K={K_probe}: {len(probed)} markers over "
                    f"{len(xpix)} momentum sites ("
                    + ", ".join(f"{v} {s}" for s, v in by_shape.items()) + ")")
                if coincident:
                    notes.append(
                        "series buried under another at k = "
                        + ", ".join(str(k) for k in coincident)
                        + " (recorded with conf=coincident)")
                if missing:
                    notes.append("sites still missing a series: k = "
                                 + ", ".join(str(k) for k in missing))
            else:
                provenance["detector"] = "blob (probe found too few)"
        else:
            provenance["detector"] = "blob (no K)"
    else:
        provenance["detector"] = "blob"

    return records, provenance


# ──────────────────────────────────────────────────────────────────────────
# Marker templates
# ──────────────────────────────────────────────────────────────────────────
#
# Measured on the 600 dpi render, over every panel of Figs. 3-6: a filled disc
# is 24 x 23 px and a ring's white hole 15 x 13 px.  So the outer radius is
# 11.5 px, the hole radius 7 px, and the stroke between them 4.5 px.  Both
# series use the same circle, which is why the *only* thing that distinguishes
# them is whether the middle is inked.
MARKER_R_600 = 11.5

# Acceptance thresholds on the template scores.  They are far apart from the
# values the other class produces -- a ring scores disc ~ 0.05-0.2 where a disc
# scores 1.0 -- so nothing here is finely tuned.
FILLED_DISC_MIN = 0.93
OPEN_RING_MIN = 0.80
OPEN_CORE_MAX = 0.45     # a tick through a hole reaches ~0.3; a disc gives 1.0
# ``core`` is what separates the classes; this only rules out a disc whose
# middle happens to drop out.  Set tight it rejects real rings drawn with a
# heavy stroke, whose inner disc is more than half ink even though the hole
# itself is empty -- Fig. 5(a)'s last marker scores disc 0.59, core 0.00.
OPEN_DISC_MAX = 0.80
# How far a marker may sit off its nominal column, as a fraction of the lattice
# spacing.  Most land within 4 px, but the print misplaces the odd one by more
# -- Fig. 6(c)'s open marker at x = 3/7 is 17 px right of its site.  A quarter
# of a spacing catches those while staying well clear of the next site.
SHIFT_FRAC = 0.25


def _template_offsets(R: float):
    """Pixel offsets of the three test regions, relative to a marker centre.

    ``ring`` straddles the stroke: inked for *both* series, so it answers "is a
    marker here".  ``core`` is well inside the hole and ``disc`` well inside the
    outer edge; together they answer "which series", with a wide margin either
    side of the stroke so that a pixel or two of registration error is harmless.
    """
    n = int(np.ceil(1.2 * R)) + 1
    yy, xx = np.mgrid[-n:n + 1, -n:n + 1]
    d = np.hypot(yy, xx)
    regions = {
        "core": d <= 0.44 * R,                       # r <= 5.0 px: deep in the hole
        "disc": d <= 0.74 * R,                       # r <= 8.5 px: inside the edge
        "ring": (d >= 0.68 * R) & (d <= 1.00 * R),   # 7.8..11.5 px: the stroke
    }
    return {k: (yy[m], xx[m]) for k, m in regions.items()}


def _score_column(ink, known, cx, rows, offsets):
    """Template scores for every candidate centre row in one column.

    Returns ``{region: array over rows}`` holding the inked fraction of that
    region, averaged over *known* pixels only.  Frame lines and tick marks are
    excluded via ``known`` rather than counted as ink: an axis line through a
    ring's hole would otherwise read as a filled disc, and a marker sitting on
    the axis -- which is where every large-x point sits -- would be misread.
    """
    H, W = ink.shape
    out = {}
    for name, (dy, dx) in offsets.items():
        yy = rows[:, None] + dy[None, :]
        xx = int(round(cx)) + dx[None, :]
        good = (yy >= 0) & (yy < H) & (xx >= 0) & (xx < W)
        yc = np.clip(yy, 0, H - 1)
        xc = np.clip(np.broadcast_to(xx, yy.shape), 0, W - 1)
        usable = good & known[yc, xc]
        n = usable.sum(axis=1)
        hit = (ink[yc, xc] & usable).sum(axis=1)
        out[name] = np.where(n > 0, hit / np.maximum(n, 1), 0.0)
        out[name + "_n"] = n
    return out


def lattice_columns(K, x_slope, x_intercept, xlim):
    """Pixel column of every allowed momentum site, ``{k: pixel}``.

    Momenta are odd integers, so ``x = k/K``.  The pixel comes from the panel's
    own axis fit, not from a fraction of the frame width -- Fig. 6(d) plots
    only 0 <= x <= 0.6, and treating its frame as a full unit of x put every
    probe column in the wrong place.
    """
    out = {}
    for k in range(1, K, 2):
        x = k / float(K)
        if not (min(xlim) - 1e-9 <= x <= max(xlim) + 1e-9):
            continue
        out[k] = (x - x_intercept) / x_slope
    return out


def trace_at_lattice(ink, frame, K, xpix, dpi=600, labels=None, big=None,
                     shapes=("filled", "open")):
    """Probe the columns where markers *must* be, instead of hunting blobs.

    Momenta are odd integers, so a marker can only sit at ``x = k/K``; ``xpix``
    supplies the pixel column for each such x, taken from the panel's own axis
    calibration (so a panel whose x range stops at 0.6 is handled like any
    other).  That turns 2D detection into a handful of 1D problems.

    Within a column the two series are found by **template score**, not by
    reading vertical runs.  Runs fail exactly where these panels are hardest:

    * every large-x marker sits *on* the x axis, so its run merges with the
      axis line and its lower half falls outside the old crop;
    * where the valence value has fallen to ~0 the higher-Fock marker sits
      almost on top of it, and two overlapping circles form one tall run that
      passes neither the marker-height nor the ring-gap test.

    A template does not care.  ``ring`` is inked for both series and locates a
    marker; ``disc`` and ``core`` then say which series it is.  The two classes
    are searched independently, so a ring resting on a disc is two detections
    at almost the same y rather than one confused blob.

    ``labels``/``min_component`` reject the legend: a data marker is drawn on
    its curve and so belongs to a huge connected component, while a legend
    sample is an island of its own size.

    Returns marker dicts with the same keys as :func:`detect_markers`.
    """
    left, right, top, bottom = frame
    R = MARKER_R_600 * dpi / 600.0
    offsets = _template_offsets(R)
    # A centre may sit anywhere from the top frame down to the axis itself;
    # q >= 0, so nothing is centred more than a pixel or two below y = 0.
    rows = np.arange(top + 2, min(bottom + 3, ink.shape[0] - int(R) - 1))
    known = _known_mask(ink, frame, R)

    spacing = (max(xpix.values()) - min(xpix.values())) / max(1, len(xpix) - 1)
    shift_px = max(3, int(round(SHIFT_FRAC * spacing)))

    out = []
    hollow = {}      # site -> best evidence of a hollow marker anywhere in it
    for k, xc in sorted(xpix.items()):
        if not (left < xc < right + R):
            continue
        best = {}
        hollow[k] = 0.0
        # Print and scan leave a few pixels of registration error; try the
        # neighbouring columns and keep whichever explains the ink best.  The
        # search stays well inside half a lattice spacing, so it can drift onto
        # a smudge but never onto the next momentum site.
        for shift in range(-shift_px, shift_px + 1):
            s = _score_column(ink, known, xc + shift, rows, offsets)
            for i, py in enumerate(rows):
                if s["ring_n"][i] < 0.4 * offsets["ring"][0].size:
                    continue
                filled = s["disc"][i] >= FILLED_DISC_MIN
                open_ = (s["ring"][i] >= OPEN_RING_MIN
                         and s["core"][i] <= OPEN_CORE_MAX
                         and s["disc"][i] <= OPEN_DISC_MAX)
                # Track the strongest hollow-looking spot in the column even
                # when it misses the cut.  Its *absence* is what later licenses
                # calling an undetected series buried rather than missed.
                if s["core"][i] <= OPEN_CORE_MAX:
                    hollow[k] = max(hollow[k], float(s["ring"][i]))
                if not (filled or open_):
                    continue
                kind = "filled" if filled else "open"
                key = (kind, int(py))
                # ``ring`` localizes both classes: it is a matched filter for
                # the stroke and falls off as the template slides off centre.
                # ``disc``/``core`` only classify -- both saturate over a range
                # of offsets, so ranking on them puts the centre wherever the
                # search happened to start.
                score = s["ring"][i]
                if key not in best or score > best[key]["score"]:
                    best[key] = dict(px=float(xc + shift), py=float(py),
                                     filled=filled, k=k, score=float(score),
                                     ring=float(s["ring"][i]),
                                     core=float(s["core"][i]),
                                     disc=float(s["disc"][i]))
        # A structure function is single-valued, so a momentum site carries at
        # most one marker per series: take the best-scoring candidate rather
        # than every local maximum.  That is what stops a steep curve from
        # yielding the same marker twice at neighbouring centres.
        for kind in ("filled", "open"):
            cands = [c for key, c in best.items() if key[0] == kind
                     and (labels is None or _on_a_curve(labels, big, c, R))]
            if not cands:
                continue
            c = max(cands, key=lambda c: c["score"])
            c["area"] = int(round(np.pi * R * R * (1.0 if c["filled"] else 0.45)))
            c["conf"] = "probe"
            c["kind"] = kind
            out.append(c)

    want = tuple(s for s in shapes if s.startswith("triangle"))
    if want:
        out += find_triangles(ink, frame, xpix, R, labels, big, shift_px, want)

    out += _recover_coincident(out, xpix, shapes, R, frame, hollow)
    out.sort(key=lambda m: (m["px"], m["py"]))
    return out


def _recover_coincident(found, xpix, shapes, R, frame, hollow):
    """Account for a series whose marker is buried under another series'.

    Two series can land on the same point, and then the paper's own plot holds
    only one visible marker.  It happens throughout the large-x tail, where
    every series has fallen to ~0 and their markers pile up on the axis -- at
    k = 17 and 19 of Fig. 6(a) the open circle is *inside* the filled disc and
    the ink there is solid black -- and again wherever two curves cross, as at
    k = 3 of Fig. 6(c).  Nothing can be detected there, by this method or any
    other.  Dropping those sites is not neutral either: it drops exactly the
    places where the two series agree and keeps the places where they differ.

    The site is therefore recorded at the covering marker's value with
    ``conf = "coincident"``, which says the two centres are within a marker of
    each other -- 0.3 in q at Fig. 6's scale.  It takes one of two licences,
    either of which pins the hidden value; without one, the site is reported
    missing rather than papered over.

    * The covering marker sits **on the axis**.  Then q = 0 there for it, and
      the buried series, being non-negative and not visible above it, is inside
      one marker of zero as well.  This is the whole large-x tail.
    * The covering marker is a **filled disc** and the column holds no hollow
      spot anywhere.  A disc is opaque, so it can hide a ring exactly; that no
      ring shows up elsewhere in the column is what says the ring is under it
      rather than somewhere the detector failed to look.  (The converse is real
      too -- an open marker drawn later paints its white middle over a disc --
      but that case only arises on the axis, where the first licence covers it.)
    """
    left, right, top, bottom = frame
    by_site = {}
    for m in found:
        by_site.setdefault(m["k"], []).append(m)

    out = []
    for k, marks in by_site.items():
        have = {m["kind"] for m in marks}
        missing = [s for s in shapes if s not in have]
        if not missing:
            continue
        on_axis = [m for m in marks if abs(m["py"] - bottom) <= 1.2 * R]
        opaque = [m for m in marks if m["kind"] == "filled"
                  and hollow.get(k, 1.0) < OPEN_RING_MIN]
        if on_axis:
            host = min(on_axis, key=lambda m: abs(m["py"] - bottom))
        elif opaque:
            host = opaque[0]
        else:
            continue
        for shape in missing:
            out.append(dict(px=host["px"], py=host["py"],
                            filled=(shape == "filled"), kind=shape, k=k,
                            score=0.0, area=host["area"], ring=None,
                            core=None, disc=None, conf="coincident"))
    return out


def _known_mask(ink, frame, R):
    """Pixels whose colour is evidence about a marker.

    The frame lines are drawn *through* the plotting area and belong to no
    marker, so they are excluded from every template average.  Excluding them,
    rather than counting them as ink, is what lets a marker centred on the x
    axis still read as a full disc or a full ring: the region simply has fewer
    usable pixels, not wrong ones.  Counting them instead is the bug that made
    the whole large-x tail unreadable -- an axis line through a ring's hole
    reads as a filled disc.

    Tick marks are deliberately *not* masked.  They are the same width as a
    ring's arc, so any rule narrow enough to catch a tick also erases the sides
    of every marker resting on the axis, which is precisely the population this
    is meant to rescue.  Left alone, a tick inside a ring's hole raises ``core``
    to about 0.3 against a filled disc's 1.0, which the threshold clears
    comfortably.
    """
    H, W = ink.shape
    known = np.ones_like(ink, dtype=bool)
    for r in np.flatnonzero(ink.sum(axis=1) > 0.5 * W):
        known[max(0, r - 1):r + 2, :] = False
    for c in np.flatnonzero(ink.sum(axis=0) > 0.5 * H):
        known[:, max(0, c - 1):c + 2] = False
    return known


def find_triangles(ink, frame, xpix, R, labels, big, shift_px,
                   want=("triangle",)):
    """Locate open triangle markers -- Fig. 6(d)'s third series, Fig. 4(b)'s
    two antiquark series.

    A triangle is found by its hole, not its outline: the white middle survives
    the axis line and the curves crossing it, and its *shape* is unmistakable.
    An apex-up triangle's hole widens steadily from top to bottom and an
    apex-down one narrows, where a circle's is widest across the middle -- so
    comparing the top third of the hole with the bottom third sorts all three
    kinds of open marker with nothing to tune.  A triangle's hole is also about
    half a ring's area, which is the cross-check rather than the test.
    """
    from scipy import ndimage

    left, right, top, bottom = frame
    holes = ndimage.binary_fill_holes(ink) ^ ink
    lab, n = ndimage.label(holes)
    out = []
    for i in range(1, n + 1):
        ys, xs_ = np.where(lab == i)
        if not (0.1 * np.pi * R * R <= ys.size <= 0.45 * np.pi * R * R):
            continue
        y0, y1, x0, x1 = ys.min(), ys.max(), xs_.min(), xs_.max()
        h, w = y1 - y0 + 1, x1 - x0 + 1
        if not (0.5 * R <= h <= 1.4 * R and 0.5 * R <= w <= 1.6 * R):
            continue
        third = max(1, h // 3)
        wide_top = np.ptp(xs_[ys < y0 + third]) + 1 if (ys < y0 + third).any() else 0
        wide_bot = np.ptp(xs_[ys > y1 - third]) + 1 if (ys > y1 - third).any() else 0
        if wide_top < 0.62 * wide_bot:
            kind = "triangle"                       # apex up
        elif wide_bot < 0.62 * wide_top:
            kind = "triangle_down"                  # apex down
        else:
            continue                                # symmetric: a ring
        if kind not in want:
            continue
        cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        k = min(xpix, key=lambda kk: abs(xpix[kk] - cx))
        if abs(xpix[k] - cx) > shift_px:
            continue
        m = dict(px=float(cx), py=float(cy), filled=False, kind=kind,
                 k=k, score=1.0, area=int(ys.size), ring=None, core=None,
                 disc=None, conf="probe")
        if not _on_a_curve(labels, big, m, R):
            continue
        out.append(m)
    # One per site per kind: a structure function is single-valued.
    best = {}
    for m in out:
        key = (m["k"], m["kind"])
        if key not in best or m["area"] > best[key]["area"]:
            best[key] = m
    return sorted(best.values(), key=lambda m: m["px"])


def curve_components(ink, R):
    """Label connected ink and flag the components big enough to be a curve.

    Returns ``(labels, big)`` where ``big[label]`` is True when that component
    is larger than a few markers' worth of ink -- in practice the frame with
    every curve welded to it, which is one component of tens of thousands of
    pixels, against a few hundred for anything free-standing.
    """
    from scipy import ndimage

    labels, n = ndimage.label(ink, structure=np.ones((3, 3), bool))
    sizes = np.bincount(labels.ravel(), minlength=n + 1)
    big = sizes > 6.0 * np.pi * R * R
    big[0] = False
    return labels, big


def _on_a_curve(labels, big, marker, R):
    """True if this marker is joined to a curve, false if it is a legend sample.

    Every plotted point is drawn on top of the line through its series, and
    every line runs into the frame, so a data marker's connected component is
    the whole curve-and-frame: thousands of pixels.  A legend sample sits in
    white space and forms a component the size of one marker.  That separates
    them far more sharply than position or size, which overlap -- measured on
    Fig. 6(a), legend blobs span 124-416 px and data blobs 164-458 px.
    """
    yy, xx = int(round(marker["py"])), int(round(marker["px"]))
    H, W = labels.shape
    # Sample the stroke, not the centre: a ring has no ink in the middle.
    seen = []
    for ang in np.linspace(0, 2 * np.pi, 32, endpoint=False):
        y = int(round(yy + 0.95 * R * np.sin(ang)))
        x = int(round(xx + 0.95 * R * np.cos(ang)))
        if 0 <= y < H and 0 <= x < W and labels[y, x]:
            seen.append(labels[y, x])
    if 0 <= yy < H and 0 <= xx < W and labels[yy, xx]:
        seen.append(labels[yy, xx])
    if not seen:
        return False
    return bool(big[np.bincount(np.asarray(seen)).argmax()])


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
        if r.get("k"):
            # The probe already knows which site it was looking at, and the
            # print misplaces the occasional marker by more than ``atol`` --
            # Fig. 6(c)'s open marker at k = 9 sits 17 px right of its column.
            # Re-deriving the site from x would throw exactly those away.
            out = dict(r)
            out["x_raw"] = r["x"]
            out["x"] = r["k"] / float(K)
            kept.append(out)
            continue
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

    # Only a structure function lives on a momentum lattice.  Fig. 2 plots
    # eigenvalues against a coupling and Fig. 8 a mass against m/g; snapping
    # those to k/K would not clean them up, it would overwrite the very
    # coordinate being read (Fig. 8's m/g = 0.020 became 0.040).
    if panel.sector:
        records, dropped = snap_to_lattice(records, K)
    else:
        dropped = []
        provenance["K_source"] += " (not applied: no momentum lattice here)"

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
            fh.write("# conf=probe: traced from the marker's own ink.  "
                     "conf=coincident: this series' marker is buried under\n"
                     "# another series' at the same site, so the value is that "
                     "marker's, good to about a marker radius.\n")
            fh.write("x,y,filled,kind,conf,k,x_raw\n")
            for r in records:
                fh.write(f"{r['x']:.6f},{r['y']:.6f},{int(r['filled'])},"
                         f"{r.get('kind', '')},{r.get('conf', 'blob')},"
                         f"{r.get('k','')},{r.get('x_raw', r['x']):.6f}\n")
        json_path = args.out / f"{panel.name}.json"
        json_path.write_text(json.dumps(provenance, indent=2))
        print(f"  wrote {csv_path}")
        print(f"  wrote {json_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
