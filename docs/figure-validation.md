# How well do our figures actually match the paper's?

Both codes can *produce* every figure. That is not the same as agreeing with
the published one. This is the per-curve accounting, including the series the
paper rescales by powers of ten.

Reproduce the table with:

```
python tools/compare_panels.py --table
```

## Side-by-side images

`tools/compare_panels.py` renders every reproducible panel as **paper |
fortran | python**, with the computed panels drawn on the axis limits read off
the paper's own frame (so a scale difference cannot hide behind autoscaling)
and the digitized markers overlaid in grey. Output in `figures/compare/`.

`tools/trace_overlay.py` draws the detections back onto the scan itself — one
ring per marker found, one cross per momentum site still short of a series.
That image is what turns a *missing* marker into something you can see, and it
is how every tracer bug below was caught.

## The tracer

The panels are 1-bit scans at ~612 ppi. Markers are circles of outer radius
11.5 px with a 7 px hole and a 4.5 px stroke, so both series are drawn with the
*same* shape and the only thing distinguishing them is whether the middle is
inked.

Detection probes the columns where markers must be — momenta are odd, so a
marker can only sit at `x = k/K` — and scores three concentric templates there:

| region | radius | filled disc | open ring |
|---|---|---|---|
| `ring` | 0.68–1.00 R | inked | inked |
| `disc` | ≤ 0.74 R | inked | mostly clear |
| `core` | ≤ 0.44 R | inked | clear |

`ring` says a marker is present and localizes it; `core` says which series it
belongs to. The two classes are searched independently, so a ring resting on a
disc is two detections rather than one confused blob.

### What was wrong before

Four bugs, all of which made markers *vanish* rather than move — the failure
mode that leaves no trace in the output:

**Frame lines counted as ink.** Every large-x marker sits *on* the x axis,
which is where a structure function has fallen to zero and where the rescaled
higher-Fock series sits on top of the valence one. The axis line running
through a ring's hole made it read as a filled disc, and the crop stopped at
the frame so the lower half of any such marker was outside the image. Frame
lines are now excluded from the template averages — treated as unknown, not as
background — and the crop extends below the axis. Tick marks are deliberately
*not* excluded: they are a stroke wide, so any rule sharp enough to catch a
tick also erases the sides of the markers this is meant to rescue.

**Two markers merged.** The old detector read vertical runs and classified by
whether a run had a gap in it. Where the two series overlap — which is most of
the large-x tail — two circles form one tall run that passes neither the height
test nor the gap test. Template scores do not care.

**The legend box ate data.** Legend suppression drew a bounding box around
everything it took for legend text and dropped detections inside it. In
Fig. 6(c) that box reaches down over the valence peak, the tallest marker in
the panel. Legend samples are now rejected by connectivity instead: a data
marker is drawn on its curve and so belongs to a component of tens of thousands
of pixels, while a legend sample is an island of its own size. Size and
position both overlap (legend blobs 124–416 px, data blobs 164–458 px);
connectivity does not.

**Panels read on the wrong axes.** Fig. 6(d) plots `0 ≤ x ≤ 0.6`, but the
probe placed its columns as fractions of the frame width, so every column was
in the wrong place. And the y axis was taken from the topmost printed number in
five panels where that number is *not* the frame edge — see below.

### Calibration

x comes from the frame edges, cross-checked against the printed tick values.

y is fitted to the panel's own **tick ladder**: the ticks are evenly spaced, so
the detections are snapped to the uniform ladder that explains most of them and
anything off it is discarded (a curve running into the left frame, or the frame
line broadened into two detections). The printed values are then attached to
their rungs by index. Reading `ylim` off the topmost printed number instead is
a trap that had five panels wrong:

| panel | printed top | its height fraction | true frame top |
|---|---|---|---|
| Fig 3(b) | 11.25 | 0.832 | **13.51** |
| Fig 4(b) | 22.5 | 0.863 | **26.62** |
| Fig 4(c) | 11.25 | 0.886 | **12.67** |
| Fig 6(a,b,c) | 10.5 | 0.720 | **14.60** |
| Fig 5(b) | — | — | **3.62**, not the 12.0 of the right-hand axis |

Fig. 5(b) is the odd one out: it carries no y labels at all, and the 0–12 axis
printed down the right-hand side of the figure belongs to Fig. 5(d), not to it.
Its tick ladder matches Fig. 5(a)'s to a pixel, which is the giveaway.

### The check that catches all of it

`∫ q dx` is the quark number — 1 for a meson, N·B for a baryon — exactly, as an
identity that never involves the solver. With the y axis calibrated from ticks
it is an *independent* test rather than a correction, and every panel it applies
to now passes. The column below is what the tick calibration alone gives,
*before* the sum rule is applied; the committed CSVs are then rescaled by that
residual factor, so they integrate to the exact value:

| panel | ∫ q dx | exact |
|---|---|---|
| Fig 3(a) | 0.993 | 1 |
| Fig 3(b) | 3.005 | 3 |
| Fig 5(a) | 0.999 | 1 |
| Fig 5(b) | 0.999 | 1 |
| Fig 5(c) | 0.999 | 1 |
| Fig 6(a) | 2.977 | 3 |
| Fig 6(b) | 2.933 | 3 |
| Fig 6(c) | 2.947 | 3 |
| Fig 6(d) | 5.923 | 6 |

Before the fixes the same check read 1.61× off on Fig. 6(c) and 1.18× on
Fig. 3(b).

### Coverage

Every marker the paper draws, found, against what the old tracer returned:

| panel | sites × series | old | now |
|---|---|---|---|
| Fig 6(a) | 10 × 2 | 14 | **20** |
| Fig 6(b) | 10 × 2 | 9 | **20** |
| Fig 6(c) | 10 × 2 | 9 | **20** |
| Fig 6(d) | 7 × 3 | 8 | **21** |
| Fig 5(a–d) | 12 × 2 | 21–23 | **23–24** |
| Fig 4(b) | 7 × 4 | 0 (fell back to blobs) | **27** |

Fig. 6(d) has three series, one of them open triangles, and Fig. 4(b) has four,
including both triangle orientations. Neither was traced at all before; a
triangle is now told from a circle by the shape of its hole, which widens
steadily from top to bottom where a circle's is widest across the middle.

`tests/test_digitized.py` asserts the coverage, the lattice positions and the
sum rule against the committed CSVs, so a regression here fails the suite
rather than quietly degrading a number.

### Markers that cannot be recovered

Where two series land on the same point the paper's own plot holds one visible
marker, and no method recovers the other. Those sites are recorded at the
covering marker's value with `conf = coincident` in the CSV — good to about a
marker radius, which is 0.3 in q at Fig. 6's scale — and never silently
dropped: dropping them is not neutral, because it drops exactly the region
where the two series agree and keeps the region where they differ. Two
conditions license it, either of which pins the hidden value: the covering
marker sits on the axis (so both series are within a marker of zero), or it is
a filled disc and the column holds no hollow spot anywhere.

## Summary

| status | curves |
|---|---|
| **agrees, ≤3%** | Fig 3(a,b) both masses; Fig 5(a,b,c) both series; Fig 6(a,b,c) valence; **Fig 6(d) all three series**; Fig 4(a) at m/g=1.6; Fig 4(b) q̄ at m/g=1.6; Fig 8(a) vs Hamer |
| **agrees, ~5%** | Fig 4(b) q at m/g=1.6 |
| **shape agrees, scale does not** | Fig 4(a,b,c) at m/g=0.1; Fig 4(c) at m/g=1.6 |
| **disagrees** | Fig 6(a,b,c) higher-Fock; Fig 5(d) |
| **not compared curve-by-curve** | Fig 2; Fig 7; Fig 8(b) |

Fortran and Python give the same answer to within 0.1 percentage point on every
entry, so the table below quotes one number for both.

## The numbers

`corr` is the Pearson correlation over the plotted points, `scale` the
least-squares multiplier that would best map ours onto the paper's (so 1.000
means the paper's stated ×10ⁿ is confirmed), and the deviations are RMS and
worst-case as a fraction of that series' own peak.

| panel | series | corr | scale | rms | max |
|---|---|---|---|---|---|
| Fig 3(a) | qq̄ valence, m/g=1.6 | 0.9999 | 0.997 | 0.4% | 0.8% |
| Fig 3(a) | qq̄ valence, m/g=0.1 | 0.9764 | 0.998 | 2.2% | 4.1% |
| Fig 3(b) | qqq valence, m/g=1.6 | 0.9998 | 0.997 | 0.7% | 1.4% |
| Fig 3(b) | qqq valence, m/g=0.1 | 0.9996 | 0.984 | 1.9% | 2.9% |
| Fig 4(a) | qq̄qq̄ ×10⁴, m/g=1.6 | 1.0000 | 1.003 | 0.4% | 0.6% |
| Fig 4(b) | q ×10³, m/g=1.6 | 0.9995 | 1.030 | 4.5% | 6.6% |
| Fig 4(b) | q̄ ×10³, m/g=1.6 | 0.9994 | 0.992 | 2.3% | 3.6% |
| Fig 5(a) | qq̄ valence | 0.9999 | 1.002 | 0.5% | 1.4% |
| Fig 5(a) | qq̄qq̄ **×10³** | 0.9996 | 0.990 | 1.1% | 2.2% |
| Fig 5(b) | qq̄ valence | 0.9999 | 1.000 | 0.6% | 1.2% |
| Fig 5(b) | qq̄qq̄ **×10²** | 0.9990 | 0.965 | 2.9% | 6.0% |
| Fig 5(c) | qq̄ valence | 1.0000 | 1.006 | 0.5% | 0.9% |
| Fig 5(c) | qq̄qq̄ **×10²** | 0.9999 | 0.977 | 1.3% | 2.3% |
| Fig 6(a) | qqq valence | 0.9999 | 1.004 | 0.5% | 0.9% |
| Fig 6(b) | qqq valence | 1.0000 | 1.018 | 1.0% | 2.1% |
| Fig 6(c) | qqq valence | 0.9998 | 1.012 | 0.9% | 2.1% |
| Fig 6(d) | 6q valence | 0.9999 | 1.009 | 1.0% | 1.2% |
| Fig 6(d) | 6q qq̄ q(x) **×5×10²** | 1.0000 | 1.011 | 0.9% | 1.4% |
| Fig 6(d) | 6q qq̄ q̄(x) **×10³** | 0.9993 | 1.009 | 1.9% | 3.2% |
| Fig 6(a) | qqqqq̄ ×10³ | 0.811 | 1.055 | 22.8% | 35.6% |
| Fig 6(b) | qqqqq̄ ×10² | 0.387 | 0.712 | 36.4% | 68.9% |
| Fig 6(c) | qqqqq̄ ×10² | 0.628 | 0.645 | 30.5% | 68.0% |
| Fig 5(d) | qq̄ ×10⁴ | −0.011 | 1.002 | 30.1% | 99.9% |
| Fig 5(d) | qq̄qq̄ ×1 | 0.450 | 0.988 | 40.8% | 100.0% |

The rescaled entries matter most: Fig. 5(a), 5(c) and all of Fig. 6(d) confirm
the paper's stated multipliers **and** our higher-Fock amplitudes
simultaneously. Getting either wrong would show up immediately, since the
multiplier and the amplitude multiply.

## What is left

**Fig 6(a,b,c) higher-Fock.** The one substantive physics disagreement.
Localized and analysed in `docs/baryon-higher-fock.md`; note that Fig. 6(d)'s
higher-Fock series now agree to 1%, which kills the standing explanation.

**Fig 5(d).** The 11th meson state. Correlation −0.01 on the valence series
means we are not plotting the same state the paper is; the state index or the
spurious-mode filter needs revisiting at that depth in the spectrum. This
panel's y scale had been forced by a sum rule that does not apply to it (its
×10⁴ sits on the *valence* curve, so the printed series does not integrate to
the quark number), which masked the problem behind a plausible-looking scale.

**Fig 4 at m/g = 0.1.** Newly compared, since the tracer previously could not
separate the two masses these panels overlay. Shapes agree (corr 0.93–0.99),
scales do not (factors of 4–40). Either the stated multipliers are being misread
off the scan or the weak-coupling runs need a different K; not yet resolved.

**Fig 2, Fig 7, Fig 8(b).** Not compared curve-by-curve. Fig. 2 is 25–40
overlapping level trajectories per panel with no momentum lattice to probe;
Fig. 7's content is Table I, which *is* reproduced (`docs/table1-units.md`);
Fig. 8(b) is three offset families of four curves, and Fig. 8(a) on the same
data is validated against Hamer's lattice points. Claiming these as
"reproduced" would not be supportable, so they are not claimed.
