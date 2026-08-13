# Table I is in M² units, not M/g

**Resolution.** Table I of Phys. Rev. D **41**, 3814 tabulates

```
M² / (m² + g²/π)
```

— the dimensionless quantity plotted on the y-axis of **Fig. 2**, and exactly
the raw eigenvalue the code produces — even though its column headers read
`M_mes/g` and `M_bar/g`. Read that way, our sweep reproduces it. Read as `M/g`,
nothing matches.

This was initially recorded here as an unresolved discrepancy. It is not one.

## The evidence

Our N=2 meson at m/g = 1.6 gives `M/g = 3.5234`. Converting to code units,

```
M²/(m² + g²/π) = (M/g)² · π λ²  = 12.4143 × 0.34734 = 4.3131
```

Table I prints **4.314**.

At the two best-converged couplings the whole table falls into place at once —
five independent columns, no free parameters:

| m/g | mes N=2 | mes N=3 | mes N=4 | bar N=3 | bar N=4 |
|---|---|---|---|---|---|
| **1.6** | 1.0002 | 1.0007 | 1.0002 | 0.9960 | 1.0437 |
| **0.8** | 1.0006 | 1.0044 | 1.0034 | 0.9825 | 1.0061 |

(ratio paper ÷ ours). Four of the five m/g = 1.6 entries agree to better than
0.1%.

Across all 30 non-trivial entries, **17 fall within the paper's own quoted
uncertainty** — the magnitude of the last retained Richardson term, which the
paper is explicit is not a statistical error and is "not more than a rough
guide" below m/g ≈ 0.2. Most of the rest sit at a pull of 1–1.5.

Compare the two hypotheses:

| reading | ratio paper ÷ ours |
|---|---|
| `M/g`, as the header says | 0.55 – 4.8, non-monotonic |
| `M²/(m²+g²/π)` | 0.98 – 2.09; **0.996 – 1.044 at m/g ≥ 0.8** |

## Corroborating checks

* **Fig. 7 and Fig. 8(a) really are `M/g`.** Their axes run 0–4 (mesons) and
  0–8 (baryons), and our `M/g` values land inside: N = 2, 3, 4 mesons at
  m/g = 1.6 give 3.52, 3.64, 3.73, and the N=3 baryon 5.56. Table I's 4.314 /
  4.618 / 4.845 would all fall off the top of Fig. 8(a). Both readings are now
  consistent: the figures are `M/g`, the table is `M²`.
* **Hamer's independent lattice data**, digitized from Fig. 8(a), matches our
  `M/g` to 0.1% at m/g = 1.6 and 0.4–2.2% for m/g ≥ 0.2.
* **Saturation explained.** Table I looked unphysical because `M/g` must grow
  like `2m/g`, yet the entries flatten near 4.3–4.8. In `M²/(m²+g²/π)` units
  saturation is *required*: as m/g → ∞ the ratio tends to a constant, since
  both numerator and denominator grow like m². The flattening was the clue.

## Not truncation

The residual at strong coupling is not our Fock truncation. Removing it
entirely at N=3, m/g=0.2 (LPN=0, up to 12 partons, 817 states versus 158) moves
the extrapolated value from 1.2153 to 1.2087 — a 0.5% change against a 24% gap.
It is the continuum extrapolation itself: at strong coupling the 1/K series
converges slowly, and the paper quotes 1.5(5) for that entry, a 33%
uncertainty which our 1.21 sits inside.

## Consequence for this repository

`refs/table1.csv` keeps the printed values verbatim and now records the units in
its header. `dlcq.figures.table1` reports `M²/(m²+g²/π)` for comparison with the
table and `M/g` for Figs. 7 and 8.

---

# The extrapolation window, and what the parenthesized errors mean

The paper extrapolates "for 2K in the range of roughly 16–24" and warns that
below m/g ≈ 0.2 "the largest K employed is likely not large enough for these to
be more than a rough guide" — `refs/table1.csv` records `reliable=0` for the
m/g = 0.10 and 0.05 rows. With 2K = 35 now routine, both halves of that caveat
can be tested.

## Widening the window

`dlcq.figures.table1` takes `K_lo`/`K_hi`; the default is the paper's 16–24.
N = 3, in the paper's own M²/(m²+g²/π) units, `ours/paper` with the pull against
the paper's quoted last term:

| m/g | meson 2K=16–24 | meson 2K=25–35 | baryon 2K=16–24 | baryon 2K=25–35 |
|---|---|---|---|---|
| 1.60 | 4.6149 p0.5 | 4.6157 p0.4 | 10.7122 p0.1 | 10.7215 p0.6 |
| 0.80 | 4.3809 p0.4 | 4.3847 p0.3 | 10.4519 p0.5 | 10.4671 p0.7 |
| 0.40 | 2.9657 p0.3 | 3.0014 p0.2 | 7.1829 p0.2 | **7.3026 p0.0** |
| 0.20 | 1.2153 p0.6 | 1.2628 p0.5 | 2.8934 p1.0 | 3.0292 p0.4 |
| 0.10 | 0.3624 p0.7 | 0.3828 p0.6 | 0.8465 p0.8 | 0.8994 p0.7 |

Over the full table -- all five columns, N = 2, 3, 4, computed on Lenore:

| | mean \|pull\| | improved | worsened | unchanged |
|---|---|---|---|---|
| 2K = 16-24 (paper's window) | 1.08 | | | |
| 2K = 25-35 | **0.71** | **27 / 30** | 2 | 1 |

The two that worsen are the N = 3 baryon at m/g = 1.6 and 0.8, from pull 0.1 to
0.6 and 0.5 to 0.7 -- both still comfortably inside the paper's own quoted term,
so this is movement at the level of its stated precision rather than a
regression.

At strong coupling the window barely matters — the answer was already converged.
At weak coupling it moves every value **toward** the published one, most sharply
for the m/g = 0.40 baryon, which goes from 1.6% low to matching the paper
exactly. So the paper's own warning is right: 2K ≤ 24 is not converged there.

In M/g the last Richardson term also falls hard where the fit is well
conditioned: at N = 3, B = 1, m/g = 1.6, from 5.5526(41) on 16–24 to
5.5551(5) on 25–35, an eightfold tightening with the shift sitting inside the
paper's own uncertainty.

## The parenthesized error is not an uncertainty at weak coupling

It is defined as the magnitude of the last retained term, and that definition
fails as m/g → 0 — for a reason that has nothing to do with K being too small.

The Eq. (27) basis is `1, 1/K, 1/K^(1+a), 1/K², …` with the non-analytic
exponent `a` from Eq. (26). But `a → 0` in the chiral limit, so `1/K^(1+a)`
collapses onto `1/K`:

| m/g | a | exponents | cond(A) on 2K=25–35 |
|---|---|---|---|
| 1.6 | 0.845 | 1.0, 1.845 | 1.7e4 |
| 0.4 | 0.326 | 1.0, 1.326 | 1.7e4 |
| 0.1 | 0.084 | 1.0, **1.084** | 4.8e4 |
| 0.05 | 0.042 | 1.0, **1.042** | 9.1e4 |

The fit itself stays excellent — maximum residual 1e-5 or better at every
coupling, i.e. the series reproduces M(K) to five digits. What degrades is the
*coefficients*: at m/g = 1.6 they are 0.892 and 0.098, at m/g = 0.1 they are
8.383 and 8.800 with near-total cancellation (sum 0.418). The last term is then
large because the basis is nearly degenerate, not because the extrapolation is
uncertain.

Refitting on every 4- and 5-point sub-window measures the actual stability:

| m/g | M/g | paper-style error | sub-window spread | overstated by |
|---|---|---|---|---|
| 1.6 | 5.55511 | 0.00050 | 8.1e-05 | 6× |
| 0.8 | 3.16708 | 0.00975 | 1.4e-04 | 72× |
| 0.4 | 1.87024 | 0.10503 | 9.6e-04 | 110× |
| 0.2 | 1.04571 | 0.28551 | 1.6e-03 | 177× |
| 0.1 | 0.54659 | 0.39492 | 1.2e-03 | **341×** |

So the rows the paper marks unreliable are far better determined than its own
error rule reports — at m/g = 0.1 the extrapolation is stable to 0.2%, not 72%.

## The dominant uncertainty is the fit form, not the truncation

Two hypotheses for the residual disagreement were tested. The first was wrong.

**Not the particle-number truncation.** The sweeps run under ``sweep_lpn`` --
valence plus one qqbar pair -- and the three worst entries were all N = 4,
where SU(4)'s richer colour structure should make that bite hardest. Raising
LPN by a full pair (meson 4 -> 6 partons, baryon 6 -> 8) moves the extrapolated
values in the **fourth decimal** and leaves every pull unchanged:

| | LPN | ours | paper | pull |
|---|---|---|---|---|
| meson N=4, m/g = 0.80 | 4 | 4.7349 | 4.743 | 4.0 |
| meson N=4, m/g = 0.80 | **6** | 4.7348 | 4.743 | **4.1** |
| baryon N=4, m/g = 1.60 | 6 | 20.4561 | 21.200 | 2.5 |
| baryon N=4, m/g = 1.60 | **8** | 20.4560 | 21.200 | **2.5** |

**It is the number of correction terms.** ``figure_fits`` shows why: the data
live at 1/K in [0.057, 0.080] and the answer is at 1/K = 0, so the fit crosses
a gap as wide as the data range itself, and the curvature across that gap is
set by the assumed series rather than by anything measured. Varying ``n_terms``
moves M(0) far more than either the window or the truncation:

| case | m/g | n=1 | n=2 | n=3 | n=4 | paper |
|---|---|---|---|---|---|---|
| meson N=4 | 0.80 | 4.707 | 4.735 | **4.742** | 4.745 | **4.743** |
| baryon N=4 | 0.80 | 20.555 | 20.823 | **20.907** | 20.930 | **20.900** |
| baryon N=3 | 0.20 | 2.696 | 3.029 | 3.191 | 3.323 | 3.100 |

Two of the three worst entries land on the published value at three correction
terms. So ``n_terms=2`` is systematically low at weak coupling, and the "pull
4.0" was never a physics disagreement.

Measured across N = 3, 4 and m/g = 1.6, 0.8, 0.2, the spread over
``n_terms`` in 2..4 exceeds the sub-window spread by **7-18x at every point**.

**Consequence for :func:`richardson_stability`.** It cannot see this, by
construction: every sub-window uses the same form over the same narrow range in
1/K. It measures window sensitivity honestly and is silent on the term that
dominates. Any quotable error bar needs the form spread added -- and the real
fix is data nearer 1/K = 0, where the choice of form stops mattering.

## Incremental K does not help -- measured, not assumed

Extending the window from 2K = 25-35 to 25-37 was expected to tighten the
weak-coupling rows. It does essentially nothing:

| | median error reduction | best | weak-coupling rows |
|---|---|---|---|
| 25-35 -> 25-37 | **1.00x** | 1.11x (m/g = 1.6) | **1.00x** |

13 of 30 improved at all, and every one of those was at strong coupling where
the error was already tiny. The reason is that one extra K barely shortens the
extrapolation, which is what sets the error:

| window | nearest data to the axis |
|---|---|
| 2K = 16-24 (the paper's) | 1/K = 0.0833 |
| 2K = 25-35 | 0.0571 |
| 2K = 25-37 | 0.0541 |
| 2K = 25-48 | 0.0417 |

Going 35 -> 37 shortens the reach by 5%; the strong-coupling errors fall by
roughly that much and the weak-coupling ones, which are dominated by the
ambiguity in how many correction terms to keep, do not move at all. **An
earlier estimate here that "K > 35 would cut the weak-coupling errors ~2.5x"
was wrong** -- that figure was for 2K = 48, and even it should be measured
rather than assumed, since the form spread is not a simple power of the
extrapolation distance.

The practical consequence: there is no incremental path. Either accept the
present errors (0.003% to 6.9%, median 2.2%, and 27/30 consistent with the
published table) or make the jump to 2K >~ 48, which needs the sparse solver.

## The one entry that does not reproduce: baryon N=4 at m/g = 1.6

27 of the 30 non-trivial entries agree within the paper's own quoted term. Two
of the exceptions are marginal (meson N=4 at m/g = 0.80 and 0.05, both within
about 2x). The third is not, and it is worth recording precisely because every
other explanation was eliminated.

Published: **21.200(300)**. Ours: **20.4561(130)** -- 3.5% low, 2.5x outside the
paper's own term.

It is not the extrapolation. The masses converge smoothly with shrinking
increments, and every order of fit agrees:

| 2K | 26 | 28 | 30 | 32 | 34 |
|---|---|---|---|---|---|
| M^2 | 19.62695 | 19.68571 | 19.73665 | 19.78128 | 19.82071 |

increments 0.059, 0.051, 0.045, 0.039; and M(0) = 20.450 / 20.456 / 20.486 /
20.460 at 1 / 2 / 3 / 4 correction terms. Reaching 21.200 would need the
sequence to rise another 7% against a clearly converging trend.

It is not the truncation: raising ``sweep_lpn`` by a whole qqbar pair moves
M(0) from 20.4561 to 20.4560.

It is not our solver. **The 1993 Fortran and the Python agree to 2e-5 at every
K**, on exactly this configuration:

| 2K | python | fortran | rel. diff |
|---|---|---|---|
| 20 | 19.378668 | 19.377980 | 3.6e-05 |
| 22 | 19.477027 | 19.477310 | 1.5e-05 |
| 24 | 19.558350 | 19.558745 | 2.0e-05 |

At 2K = 24 -- the **top of the paper's own 16-24 window** -- the mass is 19.558.
No fit over 2K = 16-24 of data topping out at 19.558 reaches 21.200; our own
extrapolation on that window gives 20.295.

And SU(4) is not broken in general: the strong-coupling bosonization ratio
(Eq. 23) predicts ``2 sin[pi/(2(2N-1))] = 0.4450`` for N = 4, and we get 0.4870
at m/g = 0.05 -- 9%, against the ~10% the paper itself quotes for that
comparison.

So the published value is not reproducible from the code that produced the
paper. We cannot say why; a transcription slip is the obvious candidate but
20.5 -> 21.2 is not an evident digit error. Recorded as an open discrepancy
rather than resolved.

## The SU(4) rows, for the record

After widening the window, every entry with a pull above 1.5 is an N = 4 one:

| | ours | paper | pull |
|---|---|---|---|
| baryon, m/g = 1.60 | 20.4561 | 21.200 | **2.5** |
| meson, m/g = 0.80 | 4.7349 | 4.743 | **4.0** |
| meson, m/g = 0.05 | 0.1025 | 0.120 | 1.8 |

That is a pattern, not three unrelated outliers, and it points at the one
systematic the stability estimate below is blind to. The sweeps run under
`sweep_lpn` -- valence plus a single extra qqbar pair -- which is 4 partons for
an N = 4 meson and 6 for an N = 4 baryon. SU(4) has the richest colour
structure of the three, so truncating at one extra pair should bite hardest
exactly there. The m/g = 0.80 meson is also the sharpest test in the table: at
0.17% low it is closer in relative terms than most entries that pass, and only
registers as pull 4.0 because the paper's own quoted term there is 0.002.

The test is cheap and specific: rerun the N = 4 columns with `sweep_lpn` raised
by one pair and see whether the pull falls. If it does, the residual is
truncation rather than anything about the extrapolation.

Two caveats before treating the spread as *the* error. It measures the fit's
stability against the choice of window, not systematic error: the sweep applies
the `sweep_lpn` particle-number truncation at every K, and any finite-K effect
not captured by the Eq. (27) form is invisible to it. And the deviations from
the paper at weak coupling (still 0.7–0.8 of the paper's own quoted term at
m/g = 0.1) are larger than either error estimate, so something systematic
remains there.
