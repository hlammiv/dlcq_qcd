# The baryon higher-Fock discrepancy — RESOLVED

**Our baryon five-quark structure function is correct.** It reproduces the
thesis's own plot of that sector to 0.07–3.9%, antiquark distribution included.
The apparent discrepancy was a digitization failure on one particular panel.

## The resolution

The article's Fig. 6(a) — and its thesis twin Fig. 12(a) — overlay the valence
curve on the five-quark curve. The two **cross three times**, and at scan
resolution the *dashes of the dashed line read as markers*. So the traced
"higher-Fock" series from that panel is not the higher-Fock curve at all, and
comparing against it made a correct calculation look 50% wrong.

The thesis also prints the same sector **on its own**: Fig. 18(a), "Contribution
to Lightest N = 3 Baryon Structure Function, from Five-Quark Wavefunction", with
q and q̄ separated, at 2K = 15 — no crossings. Against that panel
(`refs/thesis_fig18a.csv`, panel `t18a`, now a regression test):

| x | ours ×10³ | thesis | deviation |
|---|---|---|---|
| q(0.0667) | 7.933 | 7.960 | **0.3%** |
| q(0.2000) | 4.975 | 4.788 | 3.9% |
| q(0.4667) | 2.764 | 2.766 | **0.07%** |
| q(0.6000) | 0.272 | 0.278 | 2.2% |
| q̄(0.0667) | 4.710 | 4.863 | 3.1% |
| q̄(0.2000) | 0.791 | 0.955 | 17% |

Six values across two distributions, on a figure the article never printed.
Reading the figure by eye first gave the same picture: q(0.067) ≈ 8.0, q(0.20) ≈
4.7, q(0.467) ≈ 2.7, q̄(0.067) ≈ 4.7, and q̄/q at the peak 0.59 against our 0.594.

Note this is also the article's **Fig. 4(b)**, which had never been validated.

## Independently: the magnitudes were already right

The thesis Table 6 gives the four-quark probability in the lightest SU(2) meson
with uncertainties, and our sweep matches it to 1–5% at four couplings
(`refs/thesis_table6.csv`). So the Fock projection's normalization was verified
against published numbers before this, and the shape is now verified too.

## Method note

The lesson generalizes: prefer a panel that plots **one series without
crossings** over a composite panel, even if the composite is the one the article
printed. Where the article and thesis show the same physics, the thesis panel is
the better target — see `docs/figure-validation.md`.

---

## Original investigation (kept for the record)

What follows is the analysis that localized the problem before Fig. 18(a)
settled it. Its negative results still stand and are worth keeping: truncation,
the choice of function, the cluster partition and the weighting were all ruled
out, correctly — the cause was none of them.

## Where exactly it breaks## Where exactly it breaks

Correlation between the digitized paper series and our computed curve, with the
scale fitted by least squares (so shape is tested independently of any ×10ⁿ):

| series | corr | fitted scale | paper's stated scale |
|---|---|---|---|
| meson Fig 5(a) valence | 0.9999 | 0.996 | ×1 |
| **meson Fig 5(a) higher-Fock** | **0.9994** | **986.3** | **×10³** |
| **meson Fig 5(c) higher-Fock** | **0.9988** | **99.17** | **×10²** |
| baryon Fig 6(a) valence | 0.9999 | 0.993 | ×1 |
| **baryon Fig 6(a) higher-Fock** | **0.7729** | 1046 | ×10³ |

The meson higher-Fock curves recover the paper's multipliers to ~1% *and* match
its shape to 0.999. Baryon valence is equally good. Only baryon higher-Fock
fails, so this is not a general problem with our Fock projection, our
normalization, or the tracer.

## The decisive constraint: amplitude right, shape wrong

Integrating the 5-parton series over the same x points:

```
paper (y / 10³) : 0.002412
ours            : 0.002403      ratio 1.0038
```

The **total weight in the sector is right to 0.4%**. Only its distribution over
x differs:

| x | 0.048 | 0.143 | 0.238 | 0.333 | 0.429 | 0.524 | 0.714 |
|---|---|---|---|---|---|---|---|
| paper | 8.58 | 8.09 | 1.10 | 1.74 | 4.82 | 0.61 | 0.37 |
| ours | 7.36 | 4.50 | 3.59 | 4.66 | 3.64 | 1.45 | 0.03 |

The published curve is strongly bimodal with a near-node at x ≈ 0.29; ours is
smooth. Verified against a ×2 magnification of the panel — the digitized values
are what is printed.

## Ruled out

**Fock truncation.** LPN = 0 (untruncated, 189 states, up to 9 partons) and
LPN = 5 (93 states) give an identical 5-parton shape: fitted scale 1068 vs
1068, correlation −0.056 vs −0.056 against the then-contaminated trace. The
truncation level does not affect it.

**A different function.** A search over eigenstates 0–11 × sectors 3,5,7,9 ×
{q, q̄, q+q̄} found nothing above corr 0.912, and the leaders are physically
absurd (state 10's 9-parton sector at scale 10⁹). The correct candidate —
state 0, 5-parton, q(x) — is the best *sensible* one at 0.773. For contrast the
filled series identifies its own candidate at 0.9999. If the paper were
plotting q̄ or q+q̄, one of those would have correlated; none does.

**The baryon/meson cluster partition.** A 5-parton state is a 3-quark baryon
cluster plus a qq̄ meson cluster, and the code tracks that split. But plotting
only one cluster would change the sector integral, and the integral already
agrees to 0.4%.

**A tracer artifact.** Two real tracer bugs were found and fixed while chasing
this — legend text counters being read as open circles, and legend markers
landing on lattice sites surviving the lattice filter. Fixing them moved the
correlation from −0.056 to 0.773. What remains was checked by eye against a
magnified crop and is genuinely what the paper prints.

## Resolved: higher-Fock components are numerically fragile

Comparing our **two codes against each other** at 2K=21, B=1 — they agree on
the eigenvalue to 4.2e-5 and on the valence curve to 5.6e-4:

| sector | fraction of the norm | Fortran peak | Python peak | rel diff |
|---|---|---|---|---|
| 3-parton (valence) | 99.94% | 11.6025 | 11.6047 | **5.6e-4** |
| 5-parton | 0.06% | 0.0073581 | 0.0062959 | **1.6e-1** |
| 7-parton | ~0% | 6.5e-7 | 9.1e-7 | **4.0e-1** |

The basis-dependence artifact of `docs/basis-dependence.md` is amplified roughly
**300x** going from the valence sector to the 5-parton one, and ~700x to the
7-parton one. That is expected: these components are ~10^-3 of the wave
function, so the same absolute perturbation of the eigenvector is a far larger
relative change there.

Two codes that agree to 5.6e-4 on the valence differ by **16%** on the 5-parton
sector. Neither can then be expected to match a third, 1988 computation — on
different hardware, a different compiler and a different eigensolver — to better
than that. The observed gap to the published curve (~50%) is about three times
the spread between our own two codes, on a quantity where the two codes already
disagree by 16%.

So the baryon higher-Fock curves are **not a reproducibility target at the
precision the valence curves are**. They are dominated by the same ill-posed
step documented in `docs/basis-dependence.md`.

## The definition is confirmed by the thesis

Ref. 14 of the paper is K. Hornbostel, SLAC Report 333 (1988) — the thesis
behind this code. It defines (Sec. 2.4):

> `q_k = <phi(K)| b^dag_{k,c} b_{k,c} |phi(K)>` ... baryon-number sum rule
> `sum_k (q_k - qbar_k)` = 0 for mesons, N for baryons; momentum sum rule
> `sum_k k (q_k + qbar_k) = K`; in the continuum `q(x) = K q_k`.

That is exactly `dlcq.observables.structure_function`, including the factor of
K and both sum rules. The definition was never in doubt after this.

Ref. 17 (Lepage & Brodsky, Phys. Rev. D **22**, 2157) gives the standard
light-cone form `q(x) = sum_n sum_{i in n} delta(x - x_i) |psi_n|^2` — the same
probability-weighted parton count.

## The shape is confirmed independently of our tracer

The thesis reproduces this panel as its own Fig. 12(a), printed far more
legibly than the journal scan. It shows the same thing: valence peaking at
~11.7 at x ~ 0.33 (ours: 11.60), and the five-quark curve peaking at ~10.5 at
x ~ 0.09 with a **deep node at x ~ 0.27** and a second peak ~5.0 at x ~ 0.43.

So the bimodal-with-node shape is real and our smooth curve genuinely differs —
it is not a digitization artifact. The thesis also prints the panel on an
explicit 0-15.0 axis, confirming the 14.77 derived here from tick spacing.

## Also ruled out: the weighting itself

The x-space conversion in the original work was done by a **separate program,
`wf`/`wfbig`, which is not in this repository** — `qcdf.f` only emits
eigenvectors, Fock content and the basis-change matrix, and refers to `wf` in
comments. Our `dlcq.observables.structure_function` is a reimplementation of
that lost step.

It weights basis state *s* by `c_s (N c)_s` and counts partons of momentum k.
That is exact for a meson, and exact for the baryon *valence* sector — both
confirmed at 0.999. The one place it could differ without changing the sector's
total weight is how the weight is apportioned among individual partons in a
state whose colour structure carries the baryon ε-tensor *and* an extra qq̄
pair. That is the only configuration where the ε-tensor and the pair
contraction both appear, and it is exactly the configuration that fails.

Note this cannot be settled by comparing our two codes: the Python port
reproduces the Fortran's matrix elements to 1e-14, so any difference from
`wfbig` is shared by both.

## The thesis validates our higher-Fock *magnitudes*

The thesis carries something the article does not: **numbers, with
uncertainties, for a higher-Fock magnitude**. Its Table 6 gives the
probability of the four-quark component in the lightest SU(2) meson,
Richardson-extrapolated over 2K = 10-24 (`refs/thesis_table6.csv`):

| m/g | thesis | ours | pull | rel |
|---|---|---|---|---|
| 0.40 | 0.00300(20) | 0.0030328 | 0.2 | 1.1% |
| 0.20 | 0.00740(60) | 0.0077682 | 0.6 | 5.0% |
| 0.10 | 0.01110(30) | 0.0112771 | 0.6 | 1.6% |
| 0.05 | 0.01296(5) | 0.0128183 | 2.8 | 1.1% |

Three of four sit inside the thesis's own quoted uncertainty, and all four
agree to 1-5% relative. The one with a large pull is at m/g = 0.05, where the
thesis itself says "results for couplings of lambda beyond .9425 are probably
not reliable enough to allow any firm conclusions".

The thesis also restricted its space to four-quark states, noting six-quark and
higher contribute "less than .1% ... and so fall below the errors" — which is
why our LPN=4 sweep is the matching calculation.

So the Fock projection's **normalization is right**, checked against published
numbers rather than against a traced curve. Combined with the meson
higher-Fock *shape* matching at 2.3%, what remains unexplained is specifically
the baryon five-quark **shape at small x**.

## The thesis explains the shape, and flags that region as unresolved

Sec. 2.8 derives the double-bump form from momentum flow as valence quarks
split off pairs, and predicts for the SU(3) baryon "a quark distribution with
roughly equal peaks at x = .33 and x = .13 ... An antiquark peak should appear
at x = .083 with about half the height of the quark peaks."

Against ours:

| quantity | thesis text | ours |
|---|---|---|
| quark peak (upper) | 0.33 | **0.333** |
| quark peak (lower) | 0.13 | 0.048 |
| the two peaks | roughly equal | ratio 1.58 |
| antiquark peak | 0.083 | 0.048 |
| qbar/q peak height | ~0.5 | 0.651 |

Note these predicted positions are schematic — the same section says "assuming
peaks closer than about .1 are unresolved". Our *meson* four-quark curve, which
matches the published figure to 2.3%, likewise peaks at 0.042 and 0.625 rather
than the text's "~.17" and "~.5", while reproducing the stated 3:1 ratio (ours
3.90). So a position mismatch against the text is not evidence of an error.

Crucially the thesis then says of exactly this region:

> "Because of the lack of resolution, the pairs of peaks anticipated at small x
> merge into single peaks. These would presumably be resolved if K were
> sufficiently increased."

and, of higher-Fock K-sensitivity generally, that curves at 2K = 14 vs 20
differ "by at most, slightly less than 30% in the region where x is less than
.2". Our own K scan shows the same: the baryon five-quark values move
substantially between 2K = 21, 25 and 29.

The residual disagreement is therefore concentrated in precisely the region the
original author identified as unresolved and K-sensitive — on a sector where
our two codes already differ from each other by 16%.

## Status

Reported, not resolved. The affected curves are the higher-Fock series of
Fig. 6 (and by extension Fig. 4's baryon panels). Every mass, every eigenvalue,
every valence structure function and all meson higher-Fock results are
unaffected — the baryon higher-Fock sector carries ~0.06% of the ground-state
wave function.
