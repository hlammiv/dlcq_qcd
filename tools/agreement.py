#!/usr/bin/env python3
"""How far is each digitized marker from our computed curve?

`tools/compare_panels.py` renders the three-way images; this puts a number on
them.  For every marker in `refs/digitized/`, find the computed value at the
same lattice site in whichever series of that panel comes closest, and report
the relative deviation.  Matching to the nearest series rather than a declared
one is deliberate: the trace does not reliably know which curve a marker came
from, so forcing an assignment would measure the classifier, not the physics.

Valence and higher-Fock are summarized separately.  They fail for different
reasons and at different scales -- the higher-Fock curves are two to four orders
of magnitude smaller and are the ones the paper draws crossing each other.

Usage::

    python tools/agreement.py
    python tools/agreement.py --panel fig6d --verbose
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))


def series_for(provider, phys):
    """[(x, y, label, is_valence)] for one panel, from one solver."""
    from dlcq.figures import paper_lambda
    from dlcq.observables import physical_indices, structure_function

    out = []
    for n, entry in enumerate(phys["series"]):
        npart, scale, _filled, label = entry[:4]
        mg = entry[4] if len(entry) > 4 else 1.6
        r = provider.get(3, 1, phys["B"], phys["K"], paper_lambda(mg))
        idx = physical_indices(r)
        if phys["state"] >= idx.size:
            raise IndexError(f"state {phys['state']} absent")
        x, q, qbar = structure_function(r, int(idx[phys["state"]]),
                                        nparton=abs(npart))
        # The first series of a panel is normally its valence curve.  Two
        # exceptions: Fig. 3 draws two valence curves at different couplings,
        # and Figs. 4 and 18 draw NO valence curve at all -- they are dedicated
        # higher-Fock plots, so every series there is higher-Fock.
        valence = (not phys.get("no_valence")) and (n == 0 or len(entry) > 4)
        out.append((x, (qbar if npart < 0 else q) * scale, label, valence))
    return out


def markers(name):
    path = ROOT / "refs" / "digitized" / f"{name}.csv"
    if not path.exists():
        return []
    body = [ln for ln in path.read_text().splitlines() if not ln.startswith("#")]
    rows = csv.DictReader(io.StringIO("\n".join(body)))
    return [(float(r["x"]), float(r["y"])) for r in rows]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", nargs="+")
    ap.add_argument("--source", choices=("python", "fortran"), default="python")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args(argv)

    from compare_panels import PANEL_PHYSICS
    from dlcq.providers import FortranProvider, PythonProvider

    provider = (PythonProvider(ncpus=8) if args.source == "python"
                else FortranProvider(allow_run=True,
                                     extra_search=[ROOT / "python"]))
    names = args.panel or sorted(PANEL_PHYSICS)
    rows = []
    for name in names:
        phys = PANEL_PHYSICS.get(name)
        if phys is None:
            continue
        try:
            ser = series_for(provider, phys)
        except Exception as exc:
            print(f"{name}: unavailable ({exc})")
            continue
        val, high = [], []
        for xd, yd in markers(name):
            best = None
            for x, y, label, is_val in ser:
                i = int(np.argmin(np.abs(x - xd)))
                if abs(x[i] - xd) > 1e-3:
                    continue
                dev = abs(y[i] - yd) / max(abs(yd), 1e-12)
                if best is None or dev < best[0]:
                    best = (dev, y[i], label, is_val)
            if best is None:
                continue
            (val if best[3] else high).append(best[0])
            if args.verbose:
                print(f"  {name} x={xd:.4f}  paper {yd:9.3f}  ours "
                      f"{best[1]:9.3f}  {100 * best[0]:6.1f}%  {best[2]}")
        rows.append((name, val, high, bool(phys.get("partial"))))

    def med(v):
        return f"{100 * float(np.median(v)):.1f}% ({len(v)})" if v else "--"

    print(f"\nmedian |deviation| vs the published panel  [{args.source}]\n")
    print(f"{'panel':<8}{'valence':>16}{'higher-Fock':>18}")
    for name, val, high, partial in rows:
        mark = "  (partial: a published series is not computed)" if partial else ""
        print(f"{name:<8}{med(val):>16}{med(high):>18}{mark}")
    if any(r[3] for r in rows):
        print("\nPartial panels draw a second coupling on the same lattice "
              "sites, so their\nmarkers cannot be attributed by a column probe. "
              "Those panels are measured\nby hand instead -- see "
              "refs/article_fig4.csv and refs/thesis_fig18*.csv.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
