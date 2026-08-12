#!/usr/bin/env python3
"""Render pages of the scanned paper to PNG for digitization.

The article is an APS scan: one JBIG2 bitonal image per page at ~612 ppi, with
an invisible OCR text layer and **zero vector content**.  ``mutool trace``
reports only ``fill_image`` and ``ignore_text`` on every page, so there are no
curve coordinates to extract -- the only route to numbers is tracing pixels.

That the source is 1-bit is actually favourable: curves are clean strokes and
the filled/open circle markers segment cleanly by fill ratio.

Renders go to the scratch directory, never into the repository: they are
derived from copyrighted material.  Only the extracted CSVs are committed.

Usage::

    python tools/render_pages.py                 # all figure pages at 600 dpi
    python tools/render_pages.py --pages 4 --dpi 900
    python tools/render_pages.py --list          # what is on which page
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PDF = ROOT / "literature" / "PhysRevD.41.3814.pdf"
# Hornbostel's thesis (SLAC-333) reprints the same figures at better print
# quality; see CITATION.md.  Its Figs. 11 and 12 are the article's Figs. 5 and 6.
THESIS_PDF = ROOT / "literature" / "SLAC-333_Hornbostel_thesis.pdf"

SOURCES = {"article": PDF, "thesis": THESIS_PDF}

THESIS_PAGE_CONTENT = {
    81: ("thesis 68", "FIG. 11 meson wavefunctions = article FIG. 5"),
    82: ("thesis 70", "FIG. 12 baryon wavefunctions = article FIG. 6"),
}

# PDF page -> journal page -> content.  Journal pages run 3814-3821.
PAGE_CONTENT = {
    3: (3816, "FIG. 1 (vertices), FIG. 2 (spectra vs coupling)"),
    4: (3817, "FIG. 3 (valence), FIG. 4 (higher Fock), FIG. 5 (meson wf)"),
    5: (3818, "FIG. 6 (baryon wf)"),
    7: (3820, "FIG. 7 (extrapolated masses), FIG. 8 (large N), TABLE I"),
}

DEFAULT_OUT = Path("/tmp/claude-1000/-home-hlamm-Desktop-Ideas-hornbostel/"
                   "14a87109-c0c0-46a4-a343-37f4b9c69905/scratchpad/pages")


def render(pdf: Path, page: int, dpi: int, outdir: Path, prefix="p") -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    stem = outdir / f"{prefix}{page:02d}_{dpi}dpi"
    subprocess.run(
        ["pdftoppm", "-r", str(dpi), "-f", str(page), "-l", str(page),
         "-png", str(pdf), str(stem)],
        check=True, capture_output=True,
    )
    produced = sorted(outdir.glob(f"{stem.name}*.png"))
    if not produced:
        raise RuntimeError(f"pdftoppm produced nothing for page {page}")
    return produced[0]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pages", nargs="+", type=int, default=sorted(PAGE_CONTENT))
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--list", action="store_true", help="show page contents and exit")
    args = ap.parse_args(argv)

    if args.list:
        for page, (journal, what) in sorted(PAGE_CONTENT.items()):
            print(f"  PDF page {page}  (journal {journal}):  {what}")
        return 0

    if not PDF.exists():
        print(f"error: {PDF} not found.\n"
              "The article is copyrighted and not redistributed with this repo; "
              "place your own copy there. See CITATION.md.", file=sys.stderr)
        return 1

    for page in args.pages:
        path = render(PDF, page, args.dpi, args.outdir)
        journal, what = PAGE_CONTENT.get(page, ("?", "?"))
        print(f"  page {page} (journal {journal}) -> {path}")
        print(f"      {what}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
