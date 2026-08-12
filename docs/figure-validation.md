# How well do our figures actually match the paper's?

Both codes can *produce* every figure. That is not the same as agreeing with
the published one. This is the per-curve accounting, including the series the
paper rescales by powers of ten.

## The current numbers

`tools/agreement.py` puts a number on every digitized marker: it finds the
computed value at the same lattice site in whichever of that panel's series
comes closest, and reports the relative deviation. Median |deviation|, with the
marker count in brackets — Fortran and Python give the same table:

| panel | valence | higher-Fock |
|---|---|---|
| 3(a) | 1.9% (12) | — |
| 3(b) | 1.4% (9) | — |
| 5(a) | 0.6% (8) | 1.4% (10) |
| 5(b) | 0.8% (8) | 2.0% (7) |
| 5(c) | 0.3% (7) | 1.4% (9) |
| 5(d) | 2.5% (6) | 0.9% (1) |
| 6(a) | 1.4% (7) | 29.7% (5) |
| 6(b) | 0.4% (5) | 16.4% (2) |
| 6(c) | 2.0% (4) | 352.8% (3) |
| 6(d) | **0.3% (3)** | **0.8% (5)** |

Every valence curve is inside 2.5%. Every meson higher-Fock curve is inside 2%,
at multipliers of 10² to 10⁴. The three bad cells are the baryon higher-Fock series of Figs. 6(a)–(c). Those
panels *are* hard to trace — the valence and five-quark curves cross three
times, and at scan resolution the dashes of the dashed curve are the same size
as its markers — but that is not the whole story. The underlying quantity is
**not basis-independent**: it moves by up to 34% between legitimate treatments
of the non-orthogonal basis, while every mass and valence curve is stable to
one part in 10⁴. See `docs/baryon-higher-fock.md`. The same sector at 2K = 15
reproduces the thesis's uncluttered plot of it (Fig. 18(a)) to 1%.

## One number explains the whole table

The Fock basis is non-orthogonal, and there is more than one defensible way to
assemble H and to orthonormalize (`docs/basis-dependence.md`). Running each
reference case through all six assembly/policy combinations gives a
**basis-treatment spread** per Fock sector — measured at the lattice sites
carrying at least 5% of that sector's peak, since elsewhere a relative spread is
meaningless:

| run | M² | valence | next sector | next |
|---|---|---|---|---|
| B=0 2K=24 (meson) | 4.7 × 10⁻⁸ | 0.00% | 4q **0.02%** | 6q 13.3% |
| B=1 2K=15 | 6.2 × 10⁻⁶ | 0.01% | 5q **2.95%** | 7q 16.4% |
| B=1 2K=21 | 1.6 × 10⁻⁴ | 0.34% | 5q **34.8%** | 7q 60.6% |
| B=2 2K=24 | 3.7 × 10⁻⁶ | 0.00% | 8q **0.02%** | 10q 0.14% |

Set that beside how well each published curve is reproduced:

| figure | sector | spread | agreement |
|---|---|---|---|
| Fig. 5 (meson, 2K=24) | 4q | 0.02% | 1.4–2.0% |
| Fig. 18(a) (baryon, 2K=15) | 5q | 2.95% | **0.4%** |
| Fig. 18(b) (baryon, 2K=15) | 7q | **16.4%** | **17.5%** |
| Fig. 6(a)–(c) (baryon, 2K=21) | 5q | **34.8%** | not reproducible |
| Fig. 6(d) (2-baryon, 2K=24) | 8q | 0.02% | 0.8% |

**The agreement tracks the ambiguity, case for case.** Every sector determined
to better than 3% is reproduced to better than 2%. The only two we fail to
reproduce are the two whose own basis spread exceeds 15% — and Fig. 18(b) is
the cleanest illustration: its 17.5% discrepancy is the same size as its 16.4%
ambiguity, so it is not evidence of an error at all.

Masses are unaffected throughout: M² varies by at most 1.6 × 10⁻⁴ across every
variant, and the valence structure functions by at most 0.34%. Everything the
paper actually concludes — the spectra, Table I, the extrapolations, the
large-N comparison — rests on those.

## Six calibration bugs the images exposed

Getting there meant fixing six ways of mis-reading a published figure. All
were found by auditing tick *geometry* against tick *labels*, and none of them
was visible in the numbers alone.

**The frame top is often not the last labelled tick.** Fig. 6's ticks sit at
height fractions 0.185/0.365/0.541/0.721/0.895 and the labelled 10.5 is the one
at 0.721, so the top is 14.77, not 21 — reading the top as 21 inflates every y
by 1.42×. Fig. 3(b) is the same trap: 11.25 sits at 0.833, so the top is 13.5,
not 11.25.

**A shared axis may not belong to the panel beside it.** Fig. 5's right-hand
axis (0 2.4 4.8 7.2 9.6 12.0) spans panel **(d)** only — it starts at the
(b)/(d) boundary and runs down. Panel (b) shares the *left* axis with (a), which
the tick fractions confirm: (a) and (b) both read 0.163/0.334/0.498/0.663/0.830,
while (d) reads 0.201/0.400/0.600/0.796. Reading (b) against the right axis
inflated every y by 3.33×, and its median deviation fell from 70% to 0.8% when
corrected.

**A panel may plot more than one coupling.** Fig. 3 draws *two* curves per
panel — filled circles at m/g = 1.6 and open ones at m/g = 0.1 — so comparing
against a single coupling can never match half the markers. Fixing this took
Fig. 3(a) from 61% to 1.3% and Fig. 3(b) from 97% to 1.4%.

**A K the caption never states.** Fig. 6(d) is 2K = 24, not the 22 we assumed.
See `docs/inferred-K.md` for how the panel itself says so, and
`docs/baryon-higher-fock.md` for the result.

**Two prints of one figure, on two scales.** The article's Fig. 6 labels only 0
and 10.5, giving a frame top of 2.625 × 5.573 = 14.63; the thesis prints the
same three panels fully labelled 0.0 … 15.0. Since 10.5 is not a multiple of
2.5 these really are different scalings, and each panel must be read against its
own frame. Doing so took Fig. 6's valence curves to 1.4%, 0.4% and 2.0%.

**A marker classifier that fails where curves cross.** Filled-vs-open was
decided by how many runs a vertical slice broke into, which is exactly what a
crossing defeats: a curve entering a ring fills its interior. It called
Fig. 6(b) x=3/21 and Fig. 6(c) x=5/21 higher-Fock when both are valence points.
Classifying on the ink at the marker's *centre* separates the classes with a
clear gap — ≥0.71 for filled against ≤0.49 for open, nothing in between.

## The sum rule is a check, not a calibration

∫q dx is exactly the quark number, so in principle the valence curve fixes its
own vertical scale with no tick labels at all. That is how this repository used
it at first, and it was doing real work while the frame calibrations above were
still wrong.

It is no longer safe to use that way, because it can only calibrate a trace that
recovered the whole curve. Where markers are missing from the peak the sum comes
out short and the correction inflates every surviving point. Fig. 6(b) asked for
2.22× and moved good points as far as 60% off; Fig. 5(b) asked for a far more
plausible **1.11×** and still left every point uniformly 9% low — a bias that
looks exactly like a physics discrepancy.

`tools/digitize.py` now accepts it only as a refinement (within 10% of unity)
and otherwise reports the rejection with the site coverage, leaving the audited
frame calibration in place. Under that rule it confirms Figs. 3(a), 5(a), 5(c),
5(d) and 6(a) to within 6%, and declines to touch the rest.

## Side-by-side images
`tools/compare_panels.py` renders every reproducible panel as **paper |
fortran | python**, with the computed panels drawn on the axis limits read off
the paper's own frame (so a scale difference cannot hide behind autoscaling)
and the digitized markers overlaid in grey. Output in `figures/compare/`.

They make the split above visible at a glance: on Figs. 5(a) and 6(d) the grey
markers sit on our curves for *every* series including the rescaled ones; on
Fig. 6(a) they sit on the valence curve while the ×10³ series is qualitatively
right but quantitatively off.

## Thesis figures as digitization targets

Hornbostel's thesis (SLAC-333, Ref. 14 of the paper) reprints these same panels
at markedly better print quality — its Fig. 11 is the article's Fig. 5 and its
Fig. 12 is Fig. 6. `tools/digitize.py` can now trace either document
(`Panel.source`), and the thesis baryon panels are defined as `t12a/t12b/t12c`.

They trace cleanly: the column probe recovers 11–14 markers per panel at the
stated 2K=21.

**Remaining work.** The thesis draws its two series with `×` (valence) and `♦`
(five-quark) — both solid glyphs, separated by *shape*, whereas the article uses
filled and open circles separated by *fill*. The ring-vs-disc classifier
therefore mislabels them, and simple shape metrics do not separate them either
(bounding-box fill spans 0.42–0.97 with no bimodality, because the glyphs sit on
top of the curves). A classifier keyed on the curve style instead — the thesis
draws valence solid and five-quark dashed — is the natural next step; until then
the thesis panels give reliable marker *positions* but not reliable series
*assignment*, so the numbers quoted in this document still come from the
article scan.

## The tracer

Four rounds of improvement, all driven by failures the images exposed.

**Legend suppression.** These panels put legends inside the axes, and the
enclosed counters of letters like q, b, o were being read as open circles.
Separated by the shape of the arrangement -- >=3 markers sharing a y cannot be
a peaked structure function.

**Declared legend boxes.** The row rule handles a small legend, but not
Fig. 6(d)'s, which covers half the panel with three lines whose letter counters
read as open markers — and no shape rule separates those from data. Panels may
now declare the rectangle their legend occupies (`Panel.legend_box`), measured
from the render and written into the provenance JSON with everything it removed.
This is deliberately explicit rather than clever: it is auditable, and each box
is placed with the panel's own peak in view so it cannot eat real data.

**K from the row of zero-valued markers.** Where a structure function vanishes
the marker is still drawn and comes to rest on the axis, so panels that are zero
over much of their range carry a clean row of lattice sites along the bottom.
That is a far better K signal than the curve markers, and it is what settled
Fig. 6(d) — see `docs/inferred-K.md`. It is gated to speak only when it finds
≥ 7 dots with ≥ 80% matched, because the method needs the distribution to
actually vanish somewhere.

**Lattice column probe.** Momenta are odd, so a marker can only sit at
``x = k/K``. Probing those known columns turns 2D blob detection into a few 1D
problems and is immune to the merging that defeats connected components where
curves cross. Two details matter:

* an open circle's white interior masks the curve, so a vertical slice gives
  *two* runs about a gap -- which is what distinguishes it from a filled disc,
  more reliably than thresholding interior density (the two classes overlap in
  size: 124-416 px vs 164-458 px on Fig 6a);
* height alone is not enough, because a steeply climbing curve rises as far
  across the strip as a marker is tall. Markers are also *wide*. Adding a
  horizontal-extent test fixed a regression it had introduced: Fig 3(b) went
  18.9% -> 0.3%, Fig 3(a) 29.5% -> 3.0%.

Net effect versus blob detection: Fig 3(b) 1.9% -> 0.3%, Fig 6(b) 63% -> 18%,
Fig 6(a) 1.1% -> 0.9%, with the clean panels unchanged.

## Summary

| status | curves |
|---|---|
| **agrees, 1–3%** | Fig 3(a,b) valence; Fig 5(a,b,c) valence; Fig 5(a,c) rescaled; Fig 6(a) valence; Fig 8(a) vs Hamer |
| **agrees, ~15%** | Fig 5(d), both series |
| **trace too poor to conclude** | Fig 6(b,c,d); Fig 6 rescaled series; Fig 5(b) rescaled; Fig 4(b,c) |
| **not compared curve-by-curve** | Fig 2; Fig 7; Fig 8(b) |

Fortran and Python give the same answer to within 0.1 percentage point on every
entry, so the table below quotes one number for both.

## Agreeing curves

| panel | series | 2K | pts | deviation |
|---|---|---|---|---|
| Fig 3(a) | qq̄ valence | 14 | 5 | 2.6% |
| Fig 3(b) | qqq valence | 15 | 3 | 1.9% |
| Fig 5(a) | qq̄ valence | 24 | 8 | 1.2% |
| Fig 5(a) | qq̄qq̄ **×10³** | 24 | 9 | **2%** |
| Fig 5(b) | qq̄ valence | 24 | 8 | 1.3% |
| Fig 5(c) | qq̄ valence | 24 | 8 | 1.4% |
| Fig 5(c) | qq̄qq̄ **×10²** | 24 | 8 | **1%** |
| Fig 5(d) | qq̄ **×10⁴** | 24 | 8 | 14.9% |
| Fig 5(d) | qq̄qq̄ ×1 | 24 | 2 | 16.7% |
| Fig 6(a) | qqq valence | 21 | 5 | 1.1% |

The rescaled entries matter most here: Fig 5(a) and 5(c) confirm the paper's
stated ×10³ and ×10² **and** our higher-Fock amplitudes simultaneously, to 1–2%.
Getting either wrong would show up immediately, since the multiplier and the
amplitude multiply.

Note Fig 5(d) puts its ×10⁴ on the *filled* qq̄ curve, not the open one — the
11th state's valence piece is the small one there. Applying the multiplier to
the wrong series gives a 917158% "disagreement", which is how the assignment
was caught.

## Why the rest is inconclusive, and why it is not a physics failure

Fig 6(a), 6(b) and 6(c) come from **one solver run**, differing only in which
eigenstate is plotted. 6(a) agrees at 1.1%. So a 63% deviation in 6(b) cannot
be the solver.

The direct test: match each digitized panel against our first six eigenstates
and see which fits best.

| panel | expects | best match | deviations, states 0–5 |
|---|---|---|---|
| Fig 5(a) | 0 | **0** | **1%**, 190, 51, 106, 82, 44 |
| Fig 5(b) | 1 | **1** | 97, **1%**, 50, 177, 46, 58 |
| Fig 5(c) | 2 | **2** | 130, 189, **1%**, 141, 59, 40 |
| Fig 6(a) | 0 | **0** | **1%**, 37, 43, 44, 188, 74 |
| Fig 6(b) | 1 | 4 | 139, 63, 171, 53, 42, 99 |
| Fig 6(c) | 2 | 1 | 91, 42, 55, 58, 297, 55 |

Where a panel traces cleanly the correct state is picked out unambiguously — 1%
against 40–190% for every other state. Where it does not, **no** state fits:
the best is 42%. That is the signature of a degraded trace, not of a wrong
answer or an off-by-one in state ordering. An indexing error would look like a
clean 1% match at the wrong index, which is precisely what does not happen.

Those panels are also the ones the tracer flags: 3–4 surviving markers after
suppressing 10–13 legend blobs, and an inferred K of 20 against the stated 21.

## The legend problem

These panels carry their legends *inside* the axes, and two things there were
being read as data: the legend's own sample markers, and the enclosed counters
of letters like q, b, o, d, which the hole-detection channel sees as open
circles. Neither separates by size — measured on Fig 6(a), legend blobs span
124–416 px and data blobs 164–458 px — nor by position, since the valence peak
reaches 79% of the frame height.

`tools/digitize.py:suppress_legend` separates them by the *shape of the
arrangement*: three or more markers sharing a y within 1.5% of the frame height
and spanning 15% of its width cannot be a peaked structure function. This
recovered Fig 5(d)'s K (24 rather than 23) and lifted Fig 6(a)'s inlier
fraction from 0.68 to 1.00.

Before this, three legend markers in Fig 6(a) happened to land on k/21 lattice
sites and survived the lattice filter, contributing y ≈ 11–12.5 where the true
curve is ~0.

## Not compared curve-by-curve

* **Fig 2** — 25–40 overlapping level trajectories per panel; the tracer
  recovers 5–7 markers. Not enough to compare.
* **Fig 7** — three vertically offset baselines per panel, which the
  frame-based calibration does not model. Its content is Table I, which *is*
  reproduced (`docs/table1-units.md`).
* **Fig 8(b)** — three offset families of four curves. Fig 8(a), on the same
  data, is validated against Hamer's lattice points.

Claiming these as "reproduced" would not be supportable, so they are not
claimed.
