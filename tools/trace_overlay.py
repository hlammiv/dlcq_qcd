#!/usr/bin/env python3
"""Draw a trace back onto the panel it came from, to check it by eye.

The numbers in ``refs/digitized/`` are only as good as the marker detection,
and the failure that matters is silent: a marker that was never found leaves
no trace in the CSV.  This renders what the detector decided on top of the
scan -- one ring per detection, one cross per momentum site where a series is
missing -- so an omission is as visible as a mistake.

Colours::

    green circle    filled marker (valence series)
    blue circle     open marker (higher-Fock series)
    magenta circle  open triangle (Fig. 6(d) third series)
    red cross       a momentum site where fewer markers were found than the
                    panel has series
    grey column     the lattice column that was probed

Usage::

    python tools/trace_overlay.py --panel fig6a
    python tools/trace_overlay.py --all --outdir /tmp/overlays
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
    R = D.MARKER_R_600 * dpi / 600.0
    probed = prov.get("probe_records")
    marks = probed if probed is not None else records

    K = prov.get("K_probe")
    if K:
        xs = prov["x_fit"]["slope"]
        xi = prov["x_fit"]["intercept"]
        for k, xc in D.lattice_columns(K, xs, xi, panel.xlim).items():
            c = int(round(xc))
            if 0 <= c < rgb.shape[1]:
                col = rgb[frame[2]:frame[3] + 1, c]
                rgb[frame[2]:frame[3] + 1, c] = np.where(
                    (col == 255).all(axis=1)[:, None],
                    np.uint8([215, 215, 215]), col)

    colours = {"filled": np.uint8([0, 160, 0]), "open": np.uint8([0, 80, 255]),
               "triangle": np.uint8([200, 0, 200])}
    for m in marks:
        kind = m.get("kind") or ("filled" if m["filled"] else "open")
        colour = colours.get(kind, np.uint8([0, 160, 0]))
        # A marker recovered as buried under another gets a wider, orange ring:
        # it is a recorded site, not something the detector saw on its own.
        if m.get("conf") == "coincident":
            draw_circle(rgb, m["py"], m["px"], np.uint8([255, 150, 0]),
                        r=R + 8, width=3)
        else:
            draw_circle(rgb, m["py"], m["px"], colour, r=R + 3)

    for k in prov.get("sites_incomplete", []):
        xs = prov["x_fit"]["slope"]
        xi = prov["x_fit"]["intercept"]
        xc = (k / float(K) - xi) / xs
        draw_cross(rgb, frame[3] + int(1.6 * R), int(round(xc)),
                   np.uint8([230, 0, 0]), r=int(0.5 * R))
    return rgb, prov


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--pages-dir", type=Path, default=None)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args(argv)

    from PIL import Image
    names = list(D.PANELS) if args.all else [args.panel]
    args.outdir.mkdir(parents=True, exist_ok=True)
    for name in names:
        try:
            rgb, prov = overlay(name, dpi=args.dpi, pages_dir=args.pages_dir)
        except Exception as exc:                       # noqa: BLE001
            print(f"  {name}: FAILED {exc}")
            continue
        path = args.outdir / f"{name}_overlay.png"
        Image.fromarray(rgb).save(path)
        n = len(prov.get("probe_records") or [])
        print(f"  {name}: {n} probed markers, "
              f"{len(prov.get('sites_incomplete', []))} incomplete sites "
              f"-> {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
