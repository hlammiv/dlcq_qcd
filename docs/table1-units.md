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
