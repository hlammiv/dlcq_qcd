# How well do our figures actually match the paper's?

Both codes can *produce* every figure. That is not the same as agreeing with
the published one. This is the per-curve accounting, including the series the
paper rescales by powers of ten.

## Side-by-side images

`tools/compare_panels.py` renders every reproducible panel as **paper |
fortran | python**, with the computed panels drawn on the axis limits read off
the paper's own frame (so a scale difference cannot hide behind autoscaling)
and the digitized markers overlaid in grey. Output in `figures/compare/`.

They make the split below visible at a glance: on Fig 5(a) the grey markers sit
on our curves for *both* series including the x10^3; on Fig 6(a) they sit on the
valence curve while the x10^3 series is qualitatively right but quantitatively
off.

## The tracer

Two rounds of improvement, both driven by failures the images exposed.

**Legend suppression.** These panels put legends inside the axes, and the
enclosed counters of letters like q, b, o were being read as open circles.
Separated by the shape of the arrangement -- >=3 markers sharing a y cannot be
a peaked structure function.

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
