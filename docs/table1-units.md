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

Two caveats before treating the spread as *the* error. It measures the fit's
stability against the choice of window, not systematic error: the sweep applies
the `sweep_lpn` particle-number truncation at every K, and any finite-K effect
not captured by the Eq. (27) form is invisible to it. And the deviations from
the paper at weak coupling (still 0.7–0.8 of the paper's own quoted term at
m/g = 0.1) are larger than either error estimate, so something systematic
remains there.
