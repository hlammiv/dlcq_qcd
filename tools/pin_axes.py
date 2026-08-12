#!/usr/bin/env python3
"""Pin each panel's vertical scale from the positions of its y-axis labels.

The recurring error in this reproduction has been assuming the frame top equals
the last labelled tick.  It usually does not.  Five panels in the article have a
frame top above their highest label, and in every case reading it naively
produced something that looked exactly like a physics discrepancy -- Fig. 5(b)
by 233%, Fig. 6 by 44%, Fig. 4(b) by 16%.

Tick *spacing* is not a reliable substitute either: on these scans the tick scan
picks up curve and marker ink, and on Fig. 6 it gave 14.63 where the labels give
14.77.

What is reliable is where the labels are printed.  Each label is centred on the
value it marks, so given two or more label rows and the values they carry, the
map from height fraction to data value is a straight-line fit, and the frame top
follows.  This module locates the label rows; the values come from
``Panel.ylabels``, read off the render once and recorded there.

Usage::

    python tools/pin_axes.py                # every panel that declares ylabels
    python tools/pin_axes.py --panel fig6a
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))


def label_rows(ink, left, top, bottom, width=260, pad=60, min_area=25,
               row_gap=22):
    """Vertical centres of the text rows in the margin left of the frame.

    Returns ``[(y_pixel, n_glyphs), ...]`` top to bottom.  Multi-glyph rows are
    the numeric labels; single stray glyphs are usually the axis title or a
    piece of the x-axis labels below the frame.
    """
    from scipy import ndimage

    strip = ink[max(0, top - pad):bottom + pad, max(0, left - width):left - 8]
    lab, _ = ndimage.label(strip)
    ys = []
    for i, sl in enumerate(ndimage.find_objects(lab)):
        if (lab[sl] == i + 1).sum() < min_area:
            continue
        ys.append(0.5 * (sl[0].start + sl[0].stop) + max(0, top - pad))
    ys.sort()
    groups = []
    for y in ys:
        if groups and y - groups[-1][-1] <= row_gap:
            groups[-1].append(y)
        else:
            groups.append([y])
    return [(float(np.mean(g)), len(g)) for g in groups]


def fit_top(rows, values, top, bottom):
    """Least-squares frame top, given label rows and the values they carry.

    ``rows`` and ``values`` are matched in order, highest value first.  Returns
    ``(frame_top_value, value_at_frame_bottom, rms_residual)``.
    """
    fh = float(bottom - top)
    frac = np.array([(bottom - y) / fh for y, _ in rows[:len(values)]])
    vals = np.array(values, dtype=float)
    a, b = np.polyfit(frac, vals, 1)
    resid = float(np.sqrt(np.mean((a * frac + b - vals) ** 2)))
    return float(a * 1.0 + b), float(b), resid


def page_frame(panel, dpi=600):
    """Panel image and its frame, in FULL-PAGE coordinates.

    The panel crops used elsewhere start at the frame, so the label margin is
    outside them; the labels can only be found on the whole page.
    """
    import numpy as np
    from PIL import Image

    from render_pages import DEFAULT_OUT, SOURCES, render

    src = getattr(panel, "source", "article")
    pdf = SOURCES.get(src, SOURCES["article"])
    prefix = "p" if src == "article" else "t"
    pages = Path(DEFAULT_OUT)
    hits = sorted(pages.glob(f"{prefix}{panel.page:02d}_{dpi}dpi*.png"))
    if not hits:
        hits = [render(pdf, panel.page, dpi, pages, prefix=prefix)]
    ink = np.asarray(Image.open(hits[0]).convert("L")) < 128
    H, W = ink.shape
    l, t, r, b = panel.bbox
    sub = ink[int(t * H):int(b * H), int(l * W):int(r * W)]
    # 0.5, not 0.6: several top borders carry outward tick marks and the panel
    # letter, which drop their density just below a stricter cut and collapse
    # the detected frame onto the bottom border alone.
    cols = np.flatnonzero(sub.sum(axis=0) > 0.5 * sub.shape[0])
    rows = np.flatnonzero(sub.sum(axis=1) > 0.5 * sub.shape[1])
    if not cols.size or not rows.size:
        raise RuntimeError("frame not found")
    return ink, (int(cols[0] + l * W), int(cols[-1] + l * W),
                 int(rows[0] + t * H), int(rows[-1] + t * H))


def best_fit(rows, values, top, bottom):
    """Fit using the subset of label rows that matches the values best.

    A margin also holds the axis title and the top of the x-axis labels, so
    there are usually more rows than values; take the contiguous run that fits
    a straight line best.
    """
    import itertools

    # Every axis in this figure set starts at zero, so a fit whose value at the
    # frame bottom is not ~0 has picked up the axis title or a neighbouring
    # panel's labels.  Requiring it rejects those outright -- on Fig. 4(b) the
    # unconstrained best fit puts the bottom at -18.4 and the top at 39.
    scale = max(abs(v) for v in values) or 1.0
    best = feasible = None
    for pick in itertools.combinations(range(len(rows)), len(values)):
        sel = [rows[i] for i in pick]
        t_, b_, rms = fit_top(sel, values, top, bottom)
        cand = (t_, b_, rms, pick)
        if best is None or rms < best[2]:
            best = cand
        if abs(b_) < 0.03 * scale and (feasible is None or rms < feasible[2]):
            feasible = cand
    return feasible or best


def main(argv=None):
    from digitize import PANELS

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", nargs="+")
    ap.add_argument("--dpi", type=int, default=600)
    args = ap.parse_args(argv)

    names = args.panel or sorted(PANELS)
    print(f"{'panel':<8}{'declared':>10}{'fitted top':>12}{'bottom':>9}"
          f"{'rms':>8}  labels used")
    for name in names:
        panel = PANELS.get(name)
        if panel is None or not getattr(panel, "ylabels", None):
            continue
        try:
            ink, frame = page_frame(panel, dpi=args.dpi)
        except Exception as exc:
            print(f"{name:<8} unavailable ({exc})")
            continue
        left, _right, top, bottom = frame
        fh = float(bottom - top)
        rows = [r for r in label_rows(ink, left, top, bottom)
                if -0.03 <= (bottom - r[0]) / fh <= 1.05]
        vals = list(panel.ylabels)
        if len(rows) < len(vals):
            print(f"{name:<8} found only {len(rows)} label rows for "
                  f"{len(vals)} values")
            continue
        fitted, bot, rms, _pick = best_fit(rows, vals, top, bottom)
        flag = "" if abs(fitted - panel.ylim[1]) / max(panel.ylim[1], 1e-9) < 0.01 \
            else "   <-- differs from declared"
        print(f"{name:<8}{panel.ylim[1]:>10.4g}{fitted:>12.4g}{bot:>9.3f}"
              f"{rms:>8.3f}  {vals}{flag}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
