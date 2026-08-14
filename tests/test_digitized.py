"""Tier 2 -- the committed traces of the published figures.

These assert on ``refs/digitized/*.csv`` alone, so they run without the paper
PDF.  What they check is the property that actually failed in practice: not
whether a value is a little off, but whether a *marker is missing*.  A silently
incomplete trace reads as clean data, and it biases exactly the region where
two series overlap -- the large-x tail, where the valence curve has fallen to
zero and the rescaled higher-Fock curve sits on top of it.

The quark-number sum rule is the independent check on the vertical scale: it
is an identity, so it never involves the solver.

Origin, and why most of this file is expected to fail
----------------------------------------------------
Written against the ``worktree-fig6-digitize`` trace, which finds every marker
on all thirteen panels.  The traces committed here do not: twelve of thirteen
are short of markers, and the sum rule fails on five as a direct consequence.
That is a real defect in the committed data, not a defect in these tests, and
the repository already records it -- ``refs/digitized/fig6b.json`` carries

    "sum_rule_rejected_scale": 1.7042917132580284

i.e. the digitizer computed that its trace would need a 70% rescale to satisfy
an identity, and correctly declined to apply it.  A 70% correction is a missing
series, not a miscalibration.

So the incomplete panels are marked ``xfail(strict=True)`` rather than skipped
or deleted.  Strict matters: each one flips to a failure the moment its trace is
completed, which is the signal that the panel is fixed.  Deleting the assertions
would leave the defect exactly as invisible as it was before.

The two checks that pass on every panel -- the momentum lattice and the panel
bounds -- are live, unmarked, and guard the calibration that *is* right here.
"""

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
DIGITIZED = ROOT / "refs" / "digitized"

# panel -> (2K, momentum sites, quark number or None)
#
# Marker kinds are not in this table.  The trace this file came from carried a
# ``kind`` column naming the glyph (filled/open/triangle/triangle_down); the
# committed CSVs carry a binary ``filled`` flag instead, which cannot express a
# third or fourth series.  ``read_trace`` synthesises ``kind`` from that flag,
# so on panels that really do draw triangles -- fig4b, fig6d -- the triangle
# series is folded in with one of the other two rather than checked separately.
STRUCTURE_PANELS = {
    "fig3a": (14, 7, 1.0),
    "fig3b": (15, 7, 3.0),
    "fig4a": (14, 7, None),
    "fig4b": (15, 7, None),
    "fig4c": (15, 7, None),
    "fig5a": (24, 12, 1.0),
    "fig5b": (24, 12, 1.0),
    "fig5c": (24, 12, 1.0),
    "fig5d": (24, 12, None),
    "fig6a": (21, 10, 3.0),
    "fig6b": (21, 10, 3.0),
    "fig6c": (21, 10, 3.0),
    "fig6d": (24, 7, 6.0),
}

# Every panel: none of the committed traces is complete.  Measured, not
# assumed -- sites found against the counts above, then per-series gaps within
# the sites that were found:
#
#   fig3a 7/7  but 2 gaps    fig3b 6/7    fig4a 6/7    fig4b 5/7
#   fig4c 6/7                fig5a 10/12  fig5b 10/12  fig5c 10/12
#   fig5d 7/12 and 7 gaps    fig6a 8/10   fig6b 6/10   fig6c 6/10
#   fig6d 4/7
#
# fig3a is the instructive one: it reaches all seven momentum sites and passes
# the sum rule, yet is missing the *filled* marker at k=1 and k=13.  A site
# count alone would have called it clean, which is why the per-series check
# below is separate from the site count.
INCOMPLETE = set(STRUCTURE_PANELS)

# Panels where the missing markers move the sum rule outside 5%.  A subset of
# INCOMPLETE: a panel can lose a marker from the tail, where q has already
# fallen to nearly zero, and still integrate correctly -- fig5a, fig5c and
# fig6d each drop sites while staying inside tolerance.  Listing these
# separately keeps the sum rule a live check on the panels it can still police.
SUM_RULE_BROKEN = {"fig3b", "fig5b", "fig6a", "fig6b", "fig6c"}


def _params(marked):
    """Parametrise over the panels, xfailing the ones in ``marked``."""
    return [
        pytest.param(
            name,
            marks=pytest.mark.xfail(
                strict=True,
                reason="committed trace is missing markers; see this module's "
                       "docstring and sum_rule_rejected_scale in the json",
            ),
        )
        if name in marked else pytest.param(name)
        for name in sorted(STRUCTURE_PANELS)
    ]


def read_trace(name):
    """Rows of a committed trace, with ``kind`` synthesised from ``filled``."""
    path = DIGITIZED / f"{name}.csv"
    if not path.exists():
        pytest.skip(f"{path} not present")
    rows, header = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            if line.startswith("x,"):
                header = line.strip().split(",")
                continue
            row = dict(zip(header, line.strip().split(",")))
            row["kind"] = "filled" if row["filled"] == "1" else "open"
            rows.append(row)
    return rows


def frame_top(name):
    """The top of the panel in data units.

    The trace this file came from stored it as ``y_fit.frame_top``.  Here the
    vertical scale is pinned by fitting the y-label positions, so the frame top
    is the top of ``ylim`` and ``y_fit`` holds the fit itself (slope,
    intercept, max_residual) rather than the resulting bound.
    """
    return json.loads((DIGITIZED / f"{name}.json").read_text())["ylim"][1]


@pytest.mark.parametrize("name", _params(INCOMPLETE))
def test_every_series_covers_every_momentum_site(name):
    """One marker per series per site, and no site skipped.

    A structure function is single-valued, so a duplicate is as much a defect
    as an omission -- it means a steep curve was read twice at one x.
    """
    _, nsites, _ = STRUCTURE_PANELS[name]
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

    missing = [(k, kind) for k in sites for kind in ("filled", "open")
               if (k, kind) not in seen]
    assert not missing, f"{name}: series missing at {missing}"


@pytest.mark.parametrize("name", _params(set()))
def test_markers_sit_on_the_momentum_lattice(name):
    """Momenta are odd integers, so every x must be k/2K exactly.

    Live on every panel: what the committed traces get wrong is which markers
    they found, not where they placed the ones they did.
    """
    K = STRUCTURE_PANELS[name][0]
    for r in read_trace(name):
        assert abs(float(r["x"]) - int(r["k"]) / K) < 1e-6, \
            f"{name}: x={r['x']} is not k/{K} for k={r['k']}"


@pytest.mark.parametrize("name", _params(SUM_RULE_BROKEN))
def test_valence_series_carries_the_quark_number(name):
    """``int q dx = N*B`` (1 for a meson) -- an identity, not a fit.

    This is what catches a wrong vertical calibration.  Reading Fig. 6's
    y axis off its topmost printed number instead of its tick ladder, for
    instance, lands here as a 30% error.
    """
    K, _, quark_number = STRUCTURE_PANELS[name]
    if quark_number is None:
        pytest.skip("panel rescales its valence series by a power of ten")
    rows = [r for r in read_trace(name) if r["kind"] == "filled"]
    total = sum(float(r["y"]) for r in rows) * (2.0 / K)
    assert total == pytest.approx(quark_number, rel=0.05), \
        f"{name}: int q dx = {total:.3f}, expected {quark_number}"


@pytest.mark.parametrize("name", _params(set()))
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
    top = frame_top(name)
    assert min(ys) > -0.03 * top, \
        f"{name}: y = {min(ys)} is below the axis of a 0..{top:.3g} panel"
    # The tallest marker sets the scale; nothing should tower over it.
    assert max(ys) < 3.0 * sorted(ys)[-2], \
        f"{name}: one marker at y = {max(ys)} dwarfs the rest"
