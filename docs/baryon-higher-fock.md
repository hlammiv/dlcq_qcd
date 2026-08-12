# The baryon higher-Fock discrepancy

Both codes reproduce every valence curve and the *meson* higher-Fock curves to
~1%. Of the baryon higher-Fock series, **Fig. 6(d) now also agrees to ~1%**;
what does not agree is Fig. 6(a)–(c). This is what the investigation
established, and what it ruled out.

Everything below is against the corrected trace. The earlier version of this
document was written on a trace that was missing half the markers in the panels
it was drawing conclusions from — see `docs/figure-validation.md`.

## Where exactly it breaks

Correlation between the digitized paper series and our computed curve, with the
scale fitted by least squares (so shape is tested independently of any ×10ⁿ):

| series | corr | fitted scale | paper's stated scale |
|---|---|---|---|
| meson Fig 5(a) valence | 0.9999 | 1.002 | ×1 |
| meson Fig 5(a) higher-Fock | 0.9996 | 0.990 | ×10³ |
| meson Fig 5(c) higher-Fock | 0.9999 | 0.977 | ×10² |
| baryon Fig 6(a) valence | 0.9999 | 1.004 | ×1 |
| baryon Fig 6(b) valence | 1.0000 | 1.018 | ×1 |
| baryon Fig 6(c) valence | 0.9998 | 1.012 | ×1 |
| **B=2 Fig 6(d) valence** | **0.9999** | **1.009** | ×1 |
| **B=2 Fig 6(d) higher-Fock q** | **1.0000** | **1.011** | **×5×10²** |
| **B=2 Fig 6(d) higher-Fock q̄** | **0.9993** | **1.009** | **×10³** |
| **baryon Fig 6(a) higher-Fock** | **0.811** | 1.055 | ×10³ |
| **baryon Fig 6(b) higher-Fock** | **0.387** | 0.712 | ×10² |
| **baryon Fig 6(c) higher-Fock** | **0.628** | 0.645 | ×10² |

Fig. 6(d) is the entry that changed everything. Its K is not stated in the
paper and had been assumed to be the 2K = 22 of Fig. 2(c); it is 24
(`docs/inferred-K.md`), and at 24 all three of its series — including the
antiquark distribution the old comparison did not even plot — fall on ours:

```
x        : 0.042   0.125   0.208     (8-parton q, paper's x5x10^2)
paper    : 18.441  19.595  26.254
ours     : 18.296  19.231  26.063

x        : 0.042   0.125             (8-parton qbar, paper's x10^3)
paper    : 16.266   2.194
ours     : 16.133   2.077
```

So the disagreement is not "the baryon higher-Fock sector". It is Fig. 6(a),
(b) and (c) specifically: B = 1, 2K = 21, the 5-parton sector.

## The decisive constraint: amplitude right, shape wrong

Integrating each 5-parton series over the plotted points:

| panel | paper | ours | ratio |
|---|---|---|---|
| Fig 6(a) | 0.0024201 | 0.0024209 | 1.0003 |
| Fig 6(b) | 0.0075197 | 0.0081806 | 1.088 |
| Fig 6(c) | 0.0116710 | 0.0123544 | 1.059 |

The **total weight in the sector is right to 0.03–9%**. Only its distribution
over x differs, and it differs the same way in all three panels: the published
curve is high at k = 1, 3, dives to near zero at k = 5, 7, and comes back in a
bump at k = 9. Ours is smooth.

```
Fig 6(a) x : 0.048  0.143  0.238  0.333  0.429  0.524  0.619  0.714
   paper   : 8.592  8.078  1.042  1.718  4.858  0.582  0.176  0.365
   ours    : 6.296  5.673  3.860  4.773  3.483  1.144  0.178  0.012
```

## The new diagnostic: the three panels are too alike

Fig. 6(a), (b) and (c) are three different eigenstates of one run. Correlating
the paper's own traced series against each other, and ours against each other:

| pair | valence, paper | valence, ours | higher-Fock, paper | higher-Fock, ours |
|---|---|---|---|---|
| (a)–(b) | 0.842 | 0.844 | 0.979 | 0.889 |
| (a)–(c) | 0.924 | 0.926 | 0.990 | 0.930 |
| (b)–(c) | 0.732 | 0.739 | 0.990 | 0.987 |

The valence column is the control, and it is about as good as a control gets:
the paper's state-to-state relationships come out at 0.842 / 0.924 / 0.732 and
ours at 0.844 / 0.926 / 0.739. Three decimal places, on quantities nobody
tuned. The trace is faithful and the states are correctly identified.

Against that, the published higher-Fock curves for three visibly different
states are mutually correlated at 0.979–0.990 — appreciably more alike than the
states themselves are, and more alike than our calculation says they should be.

## Ruled out

**Fock truncation.** LPN = 0 (untruncated, 189 states, up to 9 partons) and
LPN = 5 (93 states) give an identical 5-parton shape. The truncation level does
not affect it.

**A different function.** A search over eigenstates 0–11 × sectors
{3,5,6,7,8,9,10} × {q, q̄, q+q̄} finds nothing sensible above corr 0.93 for
Figs. 6(a)–(c); the leaders are physically absurd (state 6's 7-parton sector at
scale 1.8×10⁵). Run against Fig. 6(d) the *same* search picks the physical
candidate outright — state 0, 8-parton, q(x), corr 1.0000, fitted scale 505.3
against the paper's stated 5×10². So the search works; it simply finds nothing
for the B = 1 panels.

**The state weighting.** `structure_function` weights basis state *s* by
`c_s (N c)_s`, which is what makes the sum rules exact in a non-orthogonal
basis; the lost `wf`/`wfbig` program might have used something else. Trying
`c²`, `|c(Nc)|`, `(Nc)²` and `|c|` in its place changes nothing: on Fig. 6(a)
they span corr 0.80–0.86 against our 0.84, and on Fig. 6(d) every one of them
gives 0.99+, because that sector is dominated by a few basis states. No
weighting convention reproduces the published shape.

**The baryon/meson cluster partition.** A 5-parton state is a 3-quark baryon
cluster plus a qq̄ meson cluster, and the code tracks that split. But plotting
only one cluster would change the sector integral, and for Fig. 6(a) the
integral agrees to 0.03%.

**The ε-tensor-plus-extra-pair hypothesis.** This was the standing explanation:
that our reimplementation apportions weight wrongly in the one configuration
where the baryon ε-tensor and an extra qq̄ pair both appear. Fig. 6(d)'s
8-parton sector *is* that configuration — six quarks in two colour-singlet
clusters plus an extra pair — and it reproduces to 1% in both q and q̄. The
hypothesis is dead.

**A tracer artifact.** Three separate tracer bugs were found and fixed while
chasing this (`docs/figure-validation.md`). What remains was checked against
overlays of the detections on the scan, and the valence control above shows the
trace reproduces state-to-state structure to three decimals.

## Status

Reported, not resolved, but much more sharply localized than before. The
affected curves are the higher-Fock series of Fig. 6(a), (b) and (c) — B = 1,
2K = 21, 5 partons. Every mass, every eigenvalue, every valence structure
function, all meson higher-Fock results and all three of Fig. 6(d)'s series are
unaffected.

The remaining possibilities are (i) something specific to the B = 1 5-parton
sector that survives all the tests above, or (ii) an error in the original
figure — for which the mutual-correlation table is the evidence, since three
different eigenstates should not share a higher-Fock shape to 0.99 while their
valence shapes differ at 0.73. Note the sector carries ~0.06% of the wave
function, which is the regime where a plotting slip is easiest to make and
hardest to notice.
