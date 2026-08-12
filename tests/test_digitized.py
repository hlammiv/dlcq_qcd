"""Tier 2 -- the committed traces of the published figures.

These assert on ``refs/digitized/*.csv`` alone, so they run without the paper
PDF.  What they check is the property that actually failed in practice: not
whether a value is a little off, but whether a *marker is missing*.  A silently
incomplete trace reads as clean data, and it biases exactly the region where
two series overlap -- the large-x tail, where the valence curve has fallen to
zero and the rescaled higher-Fock curve sits on top of it.

The quark-number sum rule is the independent check on the vertical scale: it
is an identity, so it never involves the solver.
"""

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
DIGITIZED = ROOT / "refs" / "digitized"

# panel -> (2K, momentum sites, marker kinds, quark number or None)
STRUCTURE_PANELS = {
    "fig3a": (14, 7, ("filled", "open"), 1.0),
    "fig3b": (15, 7, ("filled", "open"), 3.0),
    "fig4a": (14, 7, ("filled", "open"), None),
    "fig4b": (15, 7, ("filled", "open", "triangle", "triangle_down"), None),
    "fig4c": (15, 7, ("filled", "open"), None),
    "fig5a": (24, 12, ("filled", "open"), 1.0),
    "fig5b": (24, 12, ("filled", "open"), 1.0),
    "fig5c": (24, 12, ("filled", "open"), 1.0),
    "fig5d": (24, 12, ("filled", "open"), None),
    "fig6a": (21, 10, ("filled", "open"), 3.0),
    "fig6b": (21, 10, ("filled", "open"), 3.0),
    "fig6c": (21, 10, ("filled", "open"), 3.0),
    "fig6d": (24, 7, ("filled", "open", "triangle"), 6.0),
}

# Panels where a marker is genuinely unreadable: the paper draws two series on
# top of one another and one is not recoverable even in principle.  Each entry
# is the number of (site, series) slots allowed to be absent.
ALLOWED_GAPS = {"fig4a": 1, "fig4b": 1, "fig5d": 1}


def read_trace(name):
    path = DIGITIZED / f"{name}.csv"
    if not path.exists():
        pytest.skip(f"{path} not present")
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            if line.startswith("x,"):
                header = line.strip().split(",")
                continue
            rows.append(dict(zip(header, line.strip().split(","))))
    return rows


@pytest.mark.parametrize("name", sorted(STRUCTURE_PANELS))
def test_every_series_covers_every_momentum_site(name):
    """One marker per series per site, and no site skipped.

    A structure function is single-valued, so a duplicate is as much a defect
    as an omission -- it means a steep curve was read twice at one x.
    """
    K, nsites, kinds, _ = STRUCTURE_PANELS[name]
    rows = read_trace(name)
    seen = {}
    for r in rows:
        key = (int(r["k"]), r["kind"])
        assert key not in seen, f"{name}: two {r['kind']} markers at k={r['k']}"
        seen[key] = r
    sites = sorted({int(r["k"]) for r in rows})
    assert len(sites) == nsites, f"{name}: {len(sites)} sites, expected {nsites}"
    for k in sites:
        assert k % 2 == 1, f"{name}: k={k} is not an odd momentum"

    missing = [(k, kind) for k in sites for kind in kinds
               if (k, kind) not in seen]
    assert len(missing) <= ALLOWED_GAPS.get(name, 0), \
        f"{name}: series missing at {missing}"


@pytest.mark.parametrize("name", sorted(STRUCTURE_PANELS))
def test_markers_sit_on_the_momentum_lattice(name):
    """Momenta are odd integers, so every x must be k/2K exactly."""
    K = STRUCTURE_PANELS[name][0]
    for r in read_trace(name):
        assert abs(float(r["x"]) - int(r["k"]) / K) < 1e-6, \
            f"{name}: x={r['x']} is not k/{K} for k={r['k']}"


@pytest.mark.parametrize("name", sorted(STRUCTURE_PANELS))
def test_valence_series_carries_the_quark_number(name):
    """``int q dx = N*B`` (1 for a meson) -- an identity, not a fit.

    This is what catches a wrong vertical calibration.  Reading Fig. 6's
    y axis off its topmost printed number instead of its tick ladder, for
    instance, lands here as a 30% error.
    """
    K, _, _, quark_number = STRUCTURE_PANELS[name]
    if quark_number is None:
        pytest.skip("panel rescales its valence series by a power of ten")
    rows = [r for r in read_trace(name) if r["kind"] == "filled"]
    total = sum(float(r["y"]) for r in rows) * (2.0 / K)
    assert total == pytest.approx(quark_number, rel=0.05), \
        f"{name}: int q dx = {total:.3f}, expected {quark_number}"


@pytest.mark.parametrize("name", sorted(STRUCTURE_PANELS))
def test_traced_values_are_inside_the_panel(name):
    """No marker may sit outside its own axes.

    Legend text used to leak into the trace as data; those blobs land near the
    top of the frame, well above any curve, which is how they were caught.
    """
    rows = read_trace(name)
    ys = [float(r["y"]) for r in rows]
    # A marker drawn at q = 0 straddles the axis line, so its measured centre
    # can land a pixel or two under it.  The allowance is set by the height of
    # the panel, not by how tall its data happen to be: a marker radius is
    # about 2% of a frame.  Anything further below the axis is not a rounding
    # error.
    top = json.loads((DIGITIZED / f"{name}.json").read_text())["y_fit"]["frame_top"]
    assert min(ys) > -0.03 * top, \
        f"{name}: y = {min(ys)} is below the axis of a 0..{top:.3g} panel"
    # The tallest marker sets the scale; nothing should tower over it.
    assert max(ys) < 3.0 * sorted(ys)[-2], \
        f"{name}: one marker at y = {max(ys)} dwarfs the rest"
