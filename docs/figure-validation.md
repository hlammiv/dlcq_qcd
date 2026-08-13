# How well do our figures match the paper's?

**Current results first; method second; the history of how the numbers moved is
the appendix.**

Both solvers produce every figure. That is not the same as agreeing with the
published one, and this document is the per-curve accounting — including the
series the paper rescales by powers of ten, where getting either the multiplier
or the amplitude wrong shows up immediately because the two multiply.

Fortran and Python agree to within 0.1 percentage point on every entry below, so
one number is quoted for both.

---

## 1. Results

### 1.1 Per-panel, from the automated trace

`tools/agreement.py` finds, for every digitized marker, the computed value at the
same lattice site in whichever of that panel's series comes closest, and reports
the relative deviation. Matching to the nearest series rather than a declared one
is deliberate: the trace does not reliably know which curve a marker came from,
so forcing an assignment would measure the classifier, not the physics.

| panel | valence | higher-Fock |
|---|---|---|
| 3(a) | 1.9% (12) | — |
| 3(b) | 1.4% (9) | — |
| 5(a) | 0.6% (8) | 1.4% (10) |
| 5(b) | 0.8% (8) | 2.0% (7) |
| 5(c) | 0.3% (7) | 1.5% (9) |
| 5(d) | 2.5% (6) | 0.9% (1) |
| 6(a) | 1.4% (7) | 29.7% (5) |
| 6(b) | 0.4% (5) | 16.4% (2) |
| 6(c) | 2.0% (4) | 352.8% (3) |
| 6(d) | **0.3% (3)** | **0.8% (5)** |

Median |deviation|, marker count in brackets. Every valence curve is inside
2.5%; every meson higher-Fock curve inside 2%. The three bad cells are the
baryon higher-Fock series of Figs. 6(a)–(c) — see
`docs/baryon-higher-fock.md`, which characterizes that disagreement in full.

Figs. 4 and 18 are omitted from this table: they each draw a **second coupling**
(m/g = 0.1) whose markers sit on the same lattice sites as the m/g = 1.6 ones, so
a column probe cannot attribute them and half the traced markers can never match.
`agreement.py` flags those panels rather than reporting a number that really
measures the unmatched series. They are measured by hand instead — §1.2.

### 1.2 Every published higher-Fock curve

Fig. 4 is the paper's dedicated higher-Fock figure — the sectors plotted alone,
with no valence curve crossing them. It answers whether the baryon higher-Fock is
wrong in general or only in particular places. **Only in particular places.**

| figure | sector | 2K | agreement |
|---|---|---|---|
| Fig. 4(a) | meson qq̄qq̄ | 14 | **0.6, 0.9, 0.7, 2.5%** |
| Fig. 5(a)–(d) | meson qq̄qq̄ | 24 | **1.4–2.0%** |
| Fig. 4(b) | baryon qqqqq̄ | 15 | **2.0, 1.5, 4.0%** (11.6% in the tail) |
| thesis Fig. 18(a) | baryon qqqqq̄ | 15 | **0.4%** |
| Fig. 6(d) | two-baryon 6q+qq̄ | 24 | **0.8%** |
| Fig. 4(c) | baryon +2 pairs | 15 | 3.8, 16.9, 8.1% |
| thesis Fig. 18(b) | baryon +2 pairs | 15 | 17.5% |
| Figs. 6(a)–(c) | baryon qqqqq̄ | **21** | ~40% |

The baryon five-quark sector — the one Figs. 6(a)–(c) plot — is reproduced to
1.5–4% at 2K = 15 by two independent published panels.

### 1.3 Fig. 7 / Table I: the extrapolated masses

Fig. 7 plots exactly what Table I tabulates — the Richardson extrapolation of the
lightest state to the continuum, Eq. (27), over the paper's own 2K = 16–24 window
— so the table is the better target: numbers, rather than three curves in a
two-inch panel with a broken y axis.

All 30 non-trivial entries are compared, in M² units (`docs/table1-units.md`).
The paper states its own reliability cut, so the tolerance follows it:

> For m/g ≳ 0.2, these are reasonable estimates of the actual error. Beyond this,
> the largest K employed is likely not large enough for these to be more than a
> rough guide.

| band | entries | criterion | worst |
|---|---|---|---|
| m/g ≥ 0.8 | 10 | max(last term, 5%) | 4.2% |
| 0.2 ≤ m/g < 0.8 | 10 | max(last term, 12%) | 19.0% |
| m/g < 0.2 | 10 | max(last term, 30%) | 52.1% |

`last_term` is the magnitude of the last term in the paper's series fit, not an
error bar, and is as tight as ±0.002 — hence the combined criterion.

**Where the paper calls the extrapolation reliable, this is the sharpest
comparison in the whole reproduction.** At m/g = 1.6 the mesons are

| N | Table I | ours | |
|---|---|---|---|
| 2 | 4.314 | 4.313 | 0.02% |
| 3 | 4.618 | 4.615 | 0.07% |
| 4 | 4.845 | 4.844 | 0.02% |

and the baryons 10.710 → 10.753 (0.4%) and 21.200 → 20.312 (4.2%). Median over
the m/g ≥ 0.8 band: **0.37%**.

Below m/g ≈ 0.2 the two extrapolations drift apart, exactly where the paper says
its own are no more than a rough guide — the largest available K is not big
enough for the series to have converged, and the residual measures that, not a
disagreement about physics.

Fig. 7(b) was also read directly as a cross-check: tracking its three curves
across 191 columns and sampling at the tabulated couplings gives 4.5–12.5%,
consistent with the table and limited by panel size. Fig. 7(a) cannot be read
that way — its y axis is **broken into three offset baselines**, one per N.

### 1.4 Fig. 2: the spectra

Fig. 2 plots ~20 eigenvalue trajectories against the **coupling**, which is
continuous. It has no momentum lattice. `Panel.sector == "spectrum"` routes it to
a column scan: at each of 59 coupling columns, take the centre of every ink run
and report the *set* of levels there. No lattice, no marker classification, no
sum rule.

| panel | B | 2K | traced levels | computed | median gap | 90th pct |
|---|---|---|---|---|---|---|
| 2(a) | 0 | 10 | 8–19 | 19–21 | 1.5% | 5.4% |
| 2(b) | 1 | 13 | 11–18 | 18–21 | 1.8% | 4.4% |
| 2(c) | 2 | 22 | 10–14 | 17–21 | 1.0% | 6.2% |

Gaps as a fraction of the panel's y range. This is a weaker check than the
structure-function panels get, and honestly so: with twenty trajectories in one
frame, a traced level landing near *some* eigenvalue is not as sharp as matching
a marker to a curve.

What gives it teeth is the direction of the count. Trajectories that touch at
scan resolution merge, so the trace can only ever find **fewer** levels than
exist — never more. A trace reporting more levels than we compute would mean we
are missing states. It never does, and that is an assertion in the suite.

---

## 2. Reading a published axis

Misreading a vertical scale caused every false alarm in this project, so the
scales are not read by hand or inferred from tick spacing.
`tools/pin_axes.py` locates the y-axis **labels**, fits the height-fraction →
value map through them, and reads the frame top off the fit. Each panel records
the label values it carries in `Panel.ylabels`, so the assumption is visible and
the fit reproducible.

Two things make it trustworthy. The fit is over-determined wherever a panel
carries three or more labels, giving residuals of 0.01–0.22 in data units. And
every axis here starts at zero, so the fitted value at the frame *bottom* is an
independent check never used in the fit: it lands between −0.14 and +0.01 on all
fifteen panels. That constraint also rejects bad fits automatically — on
Fig. 4(b) the unconstrained best fit puts the bottom at −18.4 and the top at 39.

| panel | highest label | fitted top | value at bottom | rms |
|---|---|---|---|---|
| 3(a) | 3.5 | 3.588 | −0.047 | 0.049 |
| 3(b) | 11.25 | **13.67** | −0.056 | 0.106 |
| 4(a) | 21 | 20.81 | −0.139 | 0.032 |
| 4(b) | 22.5 | **26.17** | −0.572 | 0.222 |
| 4(c) | 11.25 | **12.62** | −0.051 | 0.092 |
| 5(a) | 3.6 | 3.656 | −0.039 | 0.027 |
| 5(c) | 3.6 | 3.601 | −0.023 | 0.036 |
| 6(a) | 10.5 | **14.77** | −0.068 | 0.000 |
| 6(c) | 10.5 | **14.55** | +0.013 | 0.000 |
| t12(a)–(c) | 15.0 | 15.12 / 15.07 / 15.10 | −0.03 … +0.01 | ≤0.045 |
| t13(a) | 30 | 30.17 | −0.138 | 0.035 |
| t18(a) | 25 | 25.07 | −0.043 | 0.026 |
| t18(b) | 12 | 12.08 | +0.001 | 0.011 |

### Six ways this went wrong

| panel | labelled | actual | effect if misread |
|---|---|---|---|
| 3(b) | 11.25 | 13.67 | 20% low |
| 4(b) | 22.5 | 26.17 | 16% low |
| 4(c) | 11.25 | 12.62 | 11% low |
| 5(b) | — | uses (a)'s axis, 3.6 not 12.0 | 233% high |
| 6(a)–(c) | 10.5 | 14.77 | 44% low |
| 6(d) | K not stated | 2K = 24, not 22 | shape wrong, not just scale |

Every one of those looks exactly like a physics result. Three further traps:

**A shared axis may not belong to the panel beside it.** Fig. 5's right-hand axis
(0 2.4 4.8 7.2 9.6 12.0) spans panel **(d)** only — it starts at the (b)/(d)
boundary and runs down. Panel (b) shares the *left* axis with (a), which the tick
fractions confirm: (a) and (b) both read 0.163/0.334/0.498/0.663/0.830, while (d)
reads 0.201/0.400/0.600/0.796.

**A panel may plot more than one coupling.** Fig. 3 draws *two* curves per panel —
filled circles at m/g = 1.6 and open ones at m/g = 0.1 — so comparing against a
single coupling can never match half the markers. Fixing this took Fig. 3(a) from
61% to 1.3% and Fig. 3(b) from 97% to 1.4%.

**A multiplier may sit on the series you least expect.** Fig. 5(d) puts its ×10⁴
on the *filled* qq̄ curve, not the open one — the 11th state's valence piece is
the small one there. Applying it to the wrong series gives a 917158%
"disagreement", which is how the assignment was caught.

### Two prints, two scales

The article's Fig. 6 labels only 0 and 10.5, giving a frame top of 14.77; the
thesis prints the same three panels fully labelled 0.0 … 15.0. Since 10.5 is not
a multiple of 2.5 these really are different scalings, and each panel must be
read against its own frame.

---

## 3. The sum rule is a check, not a calibration

∫q dx is exactly the quark number, so in principle the valence curve fixes its
own vertical scale with no tick labels at all. That is how this repository used
it at first, and it was doing real work while the frame calibrations above were
still wrong.

It is no longer applied. Now that every frame has been audited against its own
label positions, applying it makes things **worse** almost everywhere, because it
can only calibrate a trace that recovered the whole curve and these traces
recover about half:

| panel | sum rule applied | frame only |
|---|---|---|
| 5(a) | 3.9% | **0.6%** |
| 5(c) | 3.4% | **0.3%** |
| 6(a) | 7.9% | **2.3%** |

The failures are not obvious ones either. Fig. 6(b) asked for 2.22×, which is
visibly absurd, but Fig. 6(a) asked for 1.063× and Fig. 5(b) for 1.110× — both
entirely plausible, and both leaving every point uniformly biased in a way that
reads as a physics discrepancy.

So it stays as a diagnostic: a scale near 1 is evidence the frame was read
correctly, and a scale far from 1 says the trace is incomplete.

---

## 4. The tracer

The paper carries no vector content — `mutool trace` reports only `fill_image`
and `ignore_text` on every page — so the only route to numbers is pixel tracing.
The scan is 1-bit at ~612 ppi, which helps: curves are clean strokes and the
markers segment reasonably.

**Lattice column probe.** Momenta are odd, so a marker can only sit at `x = k/K`.
Probing those known columns turns 2D blob detection into a few 1D problems and is
immune to the merging that defeats connected components where curves cross.
Height alone is not enough — a steeply climbing curve rises as far across the
strip as a marker is tall — so a horizontal-extent test is applied too.

**Filled vs open, from the marker's centre.** The original test counted runs in a
vertical slice: an open circle's white interior masks the curve, giving two runs
about a gap. That fails exactly where curves cross, because a curve entering a
ring fills its interior. Measuring the ink in a small disc at the marker's centre
separates the classes with a clear gap — ≥0.71 for filled against ≤0.49 for open,
nothing in between — and fixed two valence points that had been called
higher-Fock (Fig. 6(b) at x = 3/21, Fig. 6(c) at x = 5/21).

**Legend suppression, two layers.** These panels put legends *inside* the axes,
and two things there get read as data: the legend's own sample markers, and the
enclosed counters of letters like q, b, o, d, which the hole channel sees as open
circles. Neither separates by size — on Fig. 6(a), legend blobs span 124–416 px
and data blobs 164–458 px — nor by position, since the valence peak reaches 79%
of the frame height. What separates them is the *shape of the arrangement*: three
or more markers sharing a y cannot be a peaked structure function. For legends
too large for that rule — Fig. 6(d)'s covers half the panel — a panel may declare
the rectangle its legend occupies (`Panel.legend_box`), measured from the render
and written into the provenance JSON with everything it removed. Deliberately
explicit rather than clever: it is auditable, and each box is placed with the
panel's own peak in view so it cannot eat real data.

**K from the row of zero-valued markers.** Where a structure function vanishes
the marker is still drawn and comes to rest on the axis, so panels that are zero
over much of their range carry a clean row of lattice sites along the bottom.
That is a far better K signal than the curve markers, and it settled Fig. 6(d) —
see `docs/inferred-K.md`. Gated to speak only with ≥ 7 dots and ≥ 80% matched,
because the method needs the distribution to actually vanish somewhere.

**Spectrum scan** for panels with no lattice — §1.4.

### Side-by-side images

`tools/compare_panels.py` renders every reproducible panel as **paper | fortran |
python**, with the computed panels drawn on the axis limits read off the paper's
own frame (so a scale difference cannot hide behind autoscaling) and the
digitized markers overlaid in grey. Output in `figures/compare/` — 22 panels.
Panels are cropped from whichever document they name, so thesis and article
versions can be compared side by side.

### Thesis panels

Hornbostel's thesis (SLAC-333, Ref. 14) reprints these panels at better print
quality — its Fig. 11 is the article's Fig. 5, Fig. 12 is Fig. 6, Fig. 13 is
Fig. 6(d). It also prints two figures the article does not: Figs. 18(a) and (b),
the five- and seven-quark baryon sectors **plotted alone**, which are the
cleanest higher-Fock targets in either document.

One limitation remains. The thesis draws its two series with `×` (valence) and
`♦`/`●` (higher-Fock) — both solid glyphs, separated by *shape*, where the
article uses filled and open circles separated by *fill*. Neither the ring-vs-disc
test nor simple shape metrics separate them (bounding-box fill spans 0.42–0.97
with no bimodality, because the glyphs sit on the curves). So the thesis panels
give reliable marker *positions* but not reliable series *assignment*, and the
numbers quoted here for them come from hand measurement recorded in
`refs/thesis_fig18a.csv`, `refs/thesis_fig18b.csv` and
`refs/thesis_fig12a_fivequark.csv`.

---

# Appendix: how the numbers moved

## A.1 Withdrawn: "a basis spread explains the table"

An earlier revision reported a basis-treatment spread per Fock sector and
concluded that how well we reproduce a published curve is set by how well that
sector is determined — 0.02% spread ↔ 1.4% agreement, 16.4% spread ↔ 17.5%
agreement, and so on across every case.

**That was withdrawn.** The spread was measured over six assembly/policy
combinations, two of which are ill posed rather than alternative conventions.
Among well-posed treatments the sectors agree to 10⁻¹² and there is no spread to
explain anything with. See `docs/baryon-higher-fock.md` §A.1.

## A.2 Superseded: the state-identification diagnostic

Before the tracer and axes were fixed, several panels could not be compared and
the question was whether that indicated a wrong answer or a degraded trace. The
test: match each digitized panel against our first six eigenstates and see which
fits best.

| panel | expects | best match | deviations, states 0–5 |
|---|---|---|---|
| Fig. 5(a) | 0 | **0** | **1%**, 190, 51, 106, 82, 44 |
| Fig. 5(b) | 1 | **1** | 97, **1%**, 50, 177, 46, 58 |
| Fig. 5(c) | 2 | **2** | 130, 189, **1%**, 141, 59, 40 |
| Fig. 6(a) | 0 | **0** | **1%**, 37, 43, 44, 188, 74 |
| Fig. 6(b) | 1 | 4 | 139, 63, 171, 53, 42, 99 |
| Fig. 6(c) | 2 | 1 | 91, 42, 55, 58, 297, 55 |

Where a panel traced cleanly the correct state was picked out unambiguously — 1%
against 40–190% for every other state. Where it did not, **no** state fitted: the
best was 42%. That is the signature of a degraded trace, not of a wrong answer or
an off-by-one in state ordering; an indexing error would look like a clean 1%
match at the wrong index, which is precisely what did not happen.

The diagnostic did its job and the panels it flagged now agree at 0.4–2.0%
(§1.1). It is kept because the reasoning generalizes.

## A.3 Superseded status table

For reference, the position before Figs. 2, 4 and 7 were compared and before the
axes were pinned:

| status | curves |
|---|---|
| agrees, 1–3% | Fig. 3(a,b) valence; Fig. 5(a,b,c) valence; Fig. 5(a,c) rescaled; Fig. 6(a) valence; Fig. 8(a) vs Hamer |
| agrees, ~15% | Fig. 5(d), both series |
| trace too poor to conclude | Fig. 6(b,c,d); Fig. 6 rescaled; Fig. 5(b) rescaled; Fig. 4(b,c) |
| not compared curve-by-curve | Fig. 2; Fig. 7; Fig. 8(b) |

Every row of that table has since been superseded except the last entry of the
third row, which became `docs/baryon-higher-fock.md`.
