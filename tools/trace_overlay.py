#!/usr/bin/env python3
"""Draw a trace back onto the panel it came from, to check it by eye.

The numbers in ``refs/digitized/`` are only as good as the marker detection,
and the failure that matters is silent: a marker that was never found leaves
no trace in the CSV.  This renders what the detector decided on top of the
scan -- one ring per detection, one cross per momentum site where a series is
missing -- so an omission is as visible as a mistake.

``tests/test_digitized.py`` measures that the committed traces are missing
markers on all thirteen structure-function panels.  This is the tool for
seeing *which* ones, and why the detector passed them over.

Colours::

    green circle    filled marker (valence series)
    blue circle     open marker (higher-Fock series)
    red cross       a momentum site where fewer series were found than the
                    panel carries
    grey column     the lattice column that was probed

Usage::

    python tools/trace_overlay.py --panel fig6a --outdir /tmp/overlays
    python tools/trace_overlay.py --all --outdir /tmp/overlays

Requires the paper PDF that ``tools/digitize.py`` reads, which is not committed
to this repository -- point ``--pages-dir`` at the rendered pages.  Everything
below is deliberately self-contained: a diagnostic must not require editing the
thing it diagnoses, so no change to ``digitize.py`` is needed to run it.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

import digitize as D

# Marker radius in pixels at 600 dpi.  ``digitize`` expresses the same scale as
# the area window ``marker_min``/``marker_max``; this is the radius that window
# brackets, and it is used here only to size the drawn rings.
MARKER_R_600 = 11.5


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


def incomplete_sites(marks, columns):
    """Momentum sites carrying fewer series than the panel does.

    Derived from the trace rather than from a per-panel table of expected
    series counts: the kinds a panel draws are whichever kinds the detector
    found *somewhere* on it, and a site is short if it is missing one of them.
    A site with no markers at all counts as short of every kind.

    That self-calibrates, which is the point -- a hard-coded expectation would
    have to be right about panels the detector is currently getting wrong.
    """
    kinds = {_kind(m) for m in marks}
    if not kinds:
        return sorted(columns)
    found = {}
    for m in marks:
        if m.get("k") is not None:
            found.setdefault(int(m["k"]), set()).add(_kind(m))
    return sorted(k for k in columns if found.get(k, set()) != kinds)


def _kind(m):
    """``filled``/``open``, from a ``kind`` field if the record has one.

    The committed CSVs carry a binary ``filled`` flag, so a panel drawing a
    third glyph -- Fig. 4(b)'s triangles, Fig. 6(d)'s -- has it folded in with
    one of these two rather than drawn in its own colour.
    """
    k = m.get("kind")
    if k:
        return k
    return "filled" if m["filled"] else "open"


def draw_circle(rgb, cy, cx, colour, r, width=2):
    h, w = rgb.shape[:2]
    for t in np.linspace(0, 2 * np.pi, 400):
        for dr in range(width):
            y = int(round(cy + (r + dr) * np.sin(t)))
            x = int(round(cx + (r + dr) * np.cos(t)))
            if 0 <= y < h and 0 <= x < w:
                rgb[y, x] = colour


def draw_cross(rgb, cy, cx, colour, r):
    h, w = rgb.shape[:2]
    for d in range(-r, r + 1):
        for off in (-1, 0, 1):
            for y, x in ((cy + d, cx + d + off), (cy + d, cx - d + off)):
                if 0 <= y < h and 0 <= x < w:
                    rgb[y, x] = colour


def overlay(panel_name, dpi=600, pages_dir=None):
    panel = D.PANELS[panel_name]
    records, prov = D.digitize(panel, dpi=dpi, pages_dir=pages_dir)
    ink, box, page_size, margin = D.load_panel_image(panel, dpi=dpi,
                                                     pages_dir=pages_dir)
    rgb = np.where(ink[:, :, None], np.uint8([0, 0, 0]), np.uint8([255, 255, 255]))
    rgb = np.ascontiguousarray(rgb)

    frame = prov["frame_px"]
    R = MARKER_R_600 * dpi / 600.0
    probed = prov.get("probe_records")
    marks = probed if probed is not None else records

    columns = {}
    K = prov.get("K_probe")
    if K:
        xs = prov["x_fit"]["slope"]
        xi = prov["x_fit"]["intercept"]
        columns = lattice_columns(K, xs, xi, panel.xlim)
        for k, xc in columns.items():
            c = int(round(xc))
            if 0 <= c < rgb.shape[1]:
                col = rgb[frame[2]:frame[3] + 1, c]
                rgb[frame[2]:frame[3] + 1, c] = np.where(
                    (col == 255).all(axis=1)[:, None],
                    np.uint8([215, 215, 215]), col)

    colours = {"filled": np.uint8([0, 160, 0]), "open": np.uint8([0, 80, 255]),
               "triangle": np.uint8([200, 0, 200])}
    for m in marks:
        colour = colours.get(_kind(m), np.uint8([0, 160, 0]))
        # A marker recovered as buried under another gets a wider, orange ring:
        # it is a recorded site, not something the detector saw on its own.
        if m.get("conf") == "coincident":
            draw_circle(rgb, m["py"], m["px"], np.uint8([255, 150, 0]),
                        r=R + 8, width=3)
        else:
            draw_circle(rgb, m["py"], m["px"], colour, r=R + 3)

    short = incomplete_sites(marks, columns)
    for k in short:
        draw_cross(rgb, frame[3] + int(1.6 * R), int(round(columns[k])),
                   np.uint8([230, 0, 0]), r=int(0.5 * R))
    return rgb, prov, marks, short


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--pages-dir", type=Path, default=None)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args(argv)

    if not args.all and not args.panel:
        ap.error("give --panel NAME or --all")

    from PIL import Image
    names = list(D.PANELS) if args.all else [args.panel]
    args.outdir.mkdir(parents=True, exist_ok=True)
    for name in names:
        try:
            rgb, prov, marks, short = overlay(name, dpi=args.dpi,
                                              pages_dir=args.pages_dir)
        except Exception as exc:                       # noqa: BLE001
            print(f"  {name}: FAILED {exc}")
            continue
        path = args.outdir / f"{name}_overlay.png"
        Image.fromarray(rgb).save(path)
        print(f"  {name}: {len(marks)} markers, "
              f"{len(short)} incomplete sites {short or ''} -> {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
