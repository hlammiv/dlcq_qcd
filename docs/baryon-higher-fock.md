# The baryon higher-Fock discrepancy — RESOLVED

**Our baryon higher-Fock structure functions are correct.** There were two
separate problems, neither of them in the physics:

1. **Figs. 6(a)–(c)** were traced from a panel where the valence and five-quark
   curves cross three times and the dashed curve's dashes read as markers. Our
   five-quark sector reproduces the thesis's own uncluttered plot of it
   (Fig. 18(a)) to 0.07–3.9%, antiquark distribution included.
2. **Fig. 6(d)** was being solved at 2K = 22 against a figure drawn at
   **2K = 24** — a K the article's caption never states. At the right K it
   agrees with the published curves to 0.14–0.8%.

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

## The last panel: Fig. 6(d) was solved at the wrong K

Fig. 6(d), the two-baryon panel, stayed stubbornly wrong long after the rest of
Fig. 6 was understood. It was not a Fock-sector problem at all.

**The article's caption states 2K = 21 for panels (a)–(c) and gives no K for
(d).** We assumed 22, on the grounds that six partons of odd momentum need an
even K and 22 is the neighbour of 21. The panel is drawn at **2K = 24**.

The panel says so itself. Wherever a structure function vanishes the marker is
still plotted and comes to rest on the axis, so a two-baryon panel — zero over
most of its range — carries a long row of zero-valued markers, one per lattice
site, at full contrast with no neighbouring curve to merge with. Measured
against the frame they fall at

    0.0427  0.1304  0.2095  0.2937  0.3763  0.4569  0.5427

which is k/24 for odd k (0.0417 0.125 0.2083 0.2917 0.375 0.4583 0.5417) to
within 0.006. At 2K = 22 the predicted positions are 0.0455 0.1364 0.2273
0.3182 0.4091 0.5 0.5909 — drifting by a full site before the right-hand edge,
far outside the measurement scatter. `tools/digitize.py` now does this
automatically (`infer_K_from_axis_dots`), and it is more trustworthy than
inferring K from the curve markers, which are few, overlapping, and biased low.

Recomputed at 2K = 24, the panel falls into place — Fortran and Python
agreeing to every printed digit, and both agreeing with the published curves:

| x | series | published | ours | deviation |
|---|---|---|---|---|
| 1/24 | 8-parton q ×5·10² | 18.448 | 18.296 | **0.8%** |
| 1/24 | 8-parton q̄ ×10³ | 16.219 | 16.133 | **0.5%** |
| 3/24 | 6-parton valence | 35.454 | 35.505 | **0.14%** |
| 3/24 | 8-parton q ×5·10² | 19.279 | 19.231 | **0.25%** |
| 5/24 | 6-parton valence | 35.541 | 35.448 | **0.26%** |

At 2K = 22 the same state has 11.3 at x = 1/22, where the published valence
curve is on the floor — the shape was wrong, not just the scale.

### All three Fock sectors, and why 22 could not work

The legend lists three curves, but they come from **two** sectors: the 6-parton
valence, and the 8-parton sector's quark *and* antiquark distributions plotted
separately. At 2K = 24 a third sector does open up. Ten partons of odd momentum,
with Pauli exclusion allowing at most N = 3 quarks per momentum, need total
momentum ≥ 24 — so the 10-parton sector exists at 24 and is forbidden at 22.
Our 2K = 24 solution carries all three, and the number sum rule closes exactly:

    6-parton   5.990892
    8-parton   0.010622  (q)   0.001517  (q̄)
    10-parton  0.000004  (q)   0.000001  (q̄)
    ------------------------------------------
    total q - qbar = 6.000000 = N*B

### The thesis prints this panel on a different y-scale

Thesis Fig. 13(a) is the same plot — same 2K = 24 lattice, same legend, same
multipliers, same curves — but its y axis is labelled 0 10 20 30 where the
article's is labelled 0 12 24 36 48. Read against those axes the two published
versions disagree by a constant 48/30 = 1.6 at every point, which is what a
point-by-point measurement gives (1.603 ± 0.009 over six points).

The number sum rule picks the article, using published numbers only: the
valence curve must satisfy ∫q dx = N·B = 6 with dx = 2/24. The article's scale
gives 5.9; the thesis's gives 3.7. So `t13a` is digitized for shape and
magnitudes are taken from `fig6d`.

Thesis Fig. 13(b)–(e) show the 2nd–5th B = 2 states from the same run and remain
available as further regression targets.

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

## Where exactly it breaks

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

## What is left, precisely

Figs. 6(a)–(c)'s five-quark curves still disagree, and the disagreement is now
pinned down to one thing: **shape, not normalization, and only at 2K = 21.**

Fitting a free scale to the published curve — so the paper's own ×10³ never
enters — returns **1055**. The stated multiplier is 10³. The sector's total
weight and normalization are therefore right to about 5%. What differs is the
contrast: the published curve has deeper valleys and higher lobes than ours.

| x | published | ours × 1.055 |
|---|---|---|
| 1/21 | 8.53 | 6.64 |
| 3/21 | 7.99 | 5.98 |
| 5/21 | **1.02** | **4.07** |
| 7/21 | **1.80** | **5.03** |
| 9/21 | 5.04 | 3.67 |
| 11/21 | 0.49 | 1.21 |

Both curves are bimodal with a minimum at the same place; the published one
simply goes much closer to zero there.

### The published curve, extracted

The panel draws this series *dashed*, and at scan resolution its dashes are the
same size as its markers — which is what made every earlier marker-based trace
of it unreliable. Extracting it column by column instead, and removing at each
column the ink run nearest our valence curve (which reproduces the published
one to 1.4%), gives the whole curve rather than six points.
`refs/thesis_fig12a_fivequark.csv` records 89 columns of it. It has:

* a first lobe peaking at **10.71 at x ≈ 0.083**,
* a **node, essentially reaching zero (0.31), at x ≈ 0.265**,
* a second lobe peaking at **5.04 at x ≈ 0.432**,
* nothing beyond x ≈ 0.55.

Ours over the same range runs 6.30, 5.67, 3.86, 4.77, 3.48, 1.14 — the same
total, spread smoothly, with no node.

### The node does not move

This is the most informative thing found. Thesis Fig. 15 is the strong-coupling
twin of Fig. 12 — the same three baryon states at m/g = 1 instead of 1.6 — and
every one of its five-quark curves has the same signature: a first lobe near
x ≈ 0.08, a node near x ≈ 0.27, a second lobe near x ≈ 0.43. Measuring the node
in Fig. 12's three panels gives x = 0.273, 0.235, 0.264.

So across **six published panels — three states, two couplings —** the node sits
at x ≈ 0.26 and the second lobe at x ≈ 0.43, essentially unmoved.

Our three states are nothing like each other over that range:

| state | 5q ×10³ at k = 1, 3, 5, 7, 9 |
|---|---|
| 0 | 6.30 5.67 3.86 4.77 3.48 |
| 1 | 19.82 13.17 19.29 23.90 6.25 |
| 2 | 29.59 21.01 26.45 34.43 15.40 |

A structure function's node structure is a property of the state, so it should
move between states and with the coupling. In the published curves it does not.
That is hard to reconcile with all six being the quantity our code computes,
though it does not by itself say what they are instead.

### What that is not

Each of these was tested and eliminated:

* **Not the scan.** The thesis prints the same three panels as its Fig. 12, on
  a different y-scale (fully labelled 0.0–15.0, against the article's single
  10.5 label and a 14.63 top). Read against its own frame the thesis gives the
  same higher-Fock numbers and the same valence agreement.
* **Not the y calibration.** Each print is now read against its own tick
  geometry, and the valence curves agree to 1.4%, 0.4% and 2.0%.
* **Not marker misassignment.** The old filled-vs-open test counted runs in a
  vertical slice, which fails exactly where the curves cross — it called
  Fig. 6(b) x=3/21 and Fig. 6(c) x=5/21 higher-Fock when both are valence
  points. `trace_at_lattice` now classifies on the ink at the marker's centre,
  which separates the two classes with a clear gap (≥0.71 against ≤0.49).
* **Not the state.** The ground state is isolated by a gap of 3.44 in M², so
  there is no degenerate subspace to mix. Its valence reproduces the published
  curve to 0.8% median, while every other state in the spectrum is 29% or worse.
* **Not a different quantity.** q̄(x), q±q̄, and the 7- and 9-parton sectors all
  correlate with the published curve *worse* than the 5-parton q(x) does.
* **Not fragility *in K*.** The sector varies smoothly and monotonically with K —
  P₅ = 3.19, 2.96, 2.74, 2.56, 2.42, 2.31 (×10⁻³) at 2K = 13…23 — with no sign
  of the instability an ill-conditioned quantity would show.  (It *is*
  ill-conditioned with respect to the basis treatment — see below — which is a
  different axis and turns out to be the relevant one.)

### And it is specific to 2K = 21

At 2K = 15 the same sector reproduces the thesis's uncluttered plot of it
(Fig. 18(a)) to 1%, in *both* distributions and at *both* couplings:

| quantity | ours | thesis | ratio |
|---|---|---|---|
| q(1/15) ×10³ | 7.933 | 7.93 | 1.000 |
| q(3/15) ×10³ | 4.975 | 4.77 | 1.043 |
| q(5/15) ×10³ | 6.239 | 6.74 | 0.926 |
| q(7/15) ×10³ | 2.764 | 2.72 | 1.016 |
| q̄(1/15) ×10³ | 4.710 | 4.72 | 0.998 |

and at m/g = 0.1 the implied multiplier is ×95.6 against a legend of ×10²,
i.e. 4.6%.

So the projection, the normalization and the ε-tensor weighting are all
confirmed. Something about the *distribution* of the five-quark weight over x
differs at the larger basis, and it is not any of the causes above.

## Withdrawn: "this quantity is not basis-independent"

An earlier revision concluded that the five-quark distribution was ambiguous at
the 34% level and that this explained the disagreement. **That conclusion was
wrong.** The 34% came from averaging over six assembly/policy combinations,
two of which apply the diagonal-only assembly without the block-wise repair
that `docs/basis-dependence.md` shows is required. Among the four well-posed
treatments the five-quark sector agrees to 0.003% at 2K = 21 and to 2 × 10⁻¹²
at 2K = 15.

So our value is well determined, the published curve is genuinely different,
and the question is open again. The table below is kept only because the
contrast between the well-posed and ill-posed variants is itself informative.

## The ill-posed variants, for the record

The Fock basis here is non-orthogonal, and there is more than one defensible
way to assemble H and to orthonormalize — see `docs/basis-dependence.md`.
Running the same 2K = 21 baryon through every combination, plus the original
1990 output itself:

| variant | M² | valence peak | 5q ×10³ at k = 1…11 |
|---|---|---|---|
| historical `qcdf.out` (1990) | 10.390380000 | 11.6025 | 7.36 4.50 3.59 4.66 3.64 1.45 |
| exact / fortran | 10.390814172 | 11.6047 | 6.30 5.67 3.86 4.77 3.48 1.14 |
| exact / blockwise | 10.390814172 | 11.6047 | 6.30 5.67 3.86 4.77 3.48 1.14 |
| exact / spectral | 10.390814080 | 11.6047 | 6.30 5.67 3.86 4.77 3.48 1.14 |
| fortran / fortran | 10.389561840 | 11.5968 | 8.06 4.91 3.68 5.06 3.85 1.57 |
| fortran / blockwise | 10.390814172 | 11.6047 | 6.30 5.67 3.86 4.77 3.48 1.14 |
| fortran / spectral | 10.389158840 | 11.6036 | 7.69 5.69 4.08 5.20 4.05 1.59 |

Spread across those variants, as (max − min)/mean:

| quantity | spread |
|---|---|
| M² | 1.6 × 10⁻⁴ |
| valence peak | 6.8 × 10⁻⁴ |
| five-quark **total** | 12% |
| five-quark **per site** | up to **34%** |

**The mass and the valence structure function are stable to about one part in
10⁴; the five-quark distribution moves by up to a third.** So the baryon
higher-Fock *distribution* is simply not a well-determined quantity in this
framework at 2K = 21, and quoting it to better than tens of percent is not
meaningful. It is the one thing in the paper whose value depends on a
convention the paper does not record.

That does not by itself explain the published shape — no variant reproduces it,
the closest still being more than twice the published value at the node — but
it does explain why this one curve resisted every attempt while everything
around it fell into place. Two regression tests now pin both halves of that
statement, and both are written to fail loudly if a future change ever does
reproduce the published curve.

### Why the same sector matches elsewhere

At 2K = 15 the five-quark sector reproduces thesis Fig. 18(a) to 1%, and the
8-parton sector of Fig. 6(d) at 2K = 24 matches to 0.8%. Those are smaller
bases with far less near-degeneracy among colour structures at fixed momentum,
so the non-orthogonality that makes the 2K = 21 attribution ambiguous is much
weaker there.

## The seven-quark discrepancy stands

An earlier revision claimed its 17.5% was covered by a 16.4% basis ambiguity.
That ambiguity was an artifact of including ill-posed variants; among well-posed
treatments the seven-quark sector at 2K = 15 agrees to 2 × 10⁻¹⁰. The 17.5%
discrepancy against thesis Fig. 18(b) is real and unexplained.

## The seven-quark sector, checked for the first time

Thesis Fig. 18(b) plots the seven-quark contribution to the same baryon, and
nothing in this repository had ever been compared with it. It is the article's
Fig. 4(c). Measured from the panel (2K = 15, y axis 0–12 confirmed by seven
evenly spaced ticks, legend `× m/g=1.6 (×10⁷)` read at 4× magnification):

| x | published ×10⁷ | ours ×10⁷ | ratio |
|---|---|---|---|
| 1/15 | 5.32 | 6.232 | 1.171 |
| 3/15 | 3.96 | 4.806 | 1.214 |
| 5/15 | 2.73 | 3.116 | 1.141 |

A systematic 17.5% high, in a sector carrying 1.9 × 10⁻⁷ of the baryon's quark
number — four orders below the five-quark sector, which reproduces the panel
directly above it on the same page to 1%. Fortran and Python agree here to
every printed digit, so it is not a porting error. `refs/thesis_fig18b.csv`
records the measurement with its provenance and the test tracks the ratio
rather than asserting agreement.

## The author's own caveat about these curves

The thesis says, of exactly these panels:

> The connecting curves are simply cubic spline fits: the resolution is in some
> cases not good enough for these to accurately depict the actual structure.

That disposes of the "node that does not move". The apparent node at x ≈ 0.26
is **spline undershoot between the two smallest markers**, which in every one of
these panels are the ones at k = 5 and k = 7. It is not a feature of the data,
and the earlier reading of it as a state- and coupling-independent physical
signature was over-reading the plot.

The markers themselves are real: at 3× magnification Fig. 12(a) plainly carries
a marker at k = 5 (≈1.14) and another at k = 7 (≈1.71), with the spline dipping
to 0.31 between them. So the comparison is marker to marker, and it stands:

| k | 1 | 3 | 5 | 7 | 9 | 11 | sum |
|---|---|---|---|---|---|---|---|
| published | 8.56 | 7.95 | 1.14 | 1.71 | 4.99 | 0.54 | 24.89 |
| ours | 6.30 | 5.67 | 3.86 | 4.77 | 3.48 | 1.14 | 25.22 |

Right total, wrong distribution — and no K reproduces it either: interpolating
our curve from 2K = 11 through 25 onto the 21-lattice tops out at correlation
0.81 with a 50% rms.

## The colour sums are exact — tested, not assumed

The revision history below made the colour sums the leading suspect, since they
are exactly what sets the higher-Fock weighting. **They are correct.**

For N = 3 the colour sum is small enough to do by brute force, so it does not
have to be taken on trust. `tools/colour_norm.py` expands each basis state over
the (type, momentum, flavour, colour) Fock basis directly from its definition —
N! signed terms per ε-contracted cluster, N per δ-contracted meson, with the
fermionic sign of sorting each term into canonical order — and takes inner
products. That is the *definition* of the norm, not a reimplementation of the
code's method.

It reproduces the solver's norm matrix **exactly**, at every K tested, with an
overall factor of 1 and a maximum relative difference of 0:

| run | sector | states | scale | max rel. diff. |
|---|---|---|---|---|
| B=1 2K=21 | 3q | 12 | 1.000000 | 0 |
| B=1 2K=21 | **5q** | **133** | **1.000000** | **0** |
| B=1 2K=21 | 7q | 248 | 1.000000 | 0 |
| B=1 2K=21 | 9q | 109 | 1.000000 | 0 |

Put beside the two earlier results, this closes the argument:

1. the norm equals the exact Gram matrix of the colour-singlet cluster basis;
2. the norm couples only states sharing a momentum configuration, which is
   precisely the condition making q(k) = K Σₛ nₖ(s) cₛ(Nc)ₛ the exact
   expectation value;
3. the eigenvector is determined to 10⁻¹² across well-posed treatments.

**Our five-quark structure function is therefore correct by construction, not by
convention, and the published curve is not reproducible by a correct
calculation.** Whatever produced it — an earlier code path, a different
definition, or an error — it is not what `qcdf.f` computes and not what the
theory prescribes.

## Superseded candidate: the code postdates the figures

`fortran/qcdf.f` carries its own revision history:

```
**** 6/24/88 MODIFIED QCD2A2 SO THAT COLOR ****
**** SUMS ARE PERFORMED DIAGRAMMATICALLY   ****
**** RATHER THAN ITERATIVELY               ****

**** 5/8/90 modified qcdf to include arbitrary number of quark flavors ****

(C) Kent Hornbostel 1993.  All Rights Reserved.
```

The thesis is 1988 and the article 1990. **The code we have is a 1993 snapshot,
with the colour sums rewritten in June 1988 and a further rewrite for flavour in
1990** — the latter still carrying `????` debug markers the author intended to
remove "after convinced program is working".

The colour sums are precisely what sets the higher-Fock weighting: they build
the norm matrix, and the structure function's per-configuration attribution is
c·(Nc). A change there would move higher-Fock *distributions* while leaving
alone everything that agrees:

* masses — dominated by the valence sector;
* valence structure functions — normalized by the sum rule;
* higher-Fock *totals* — protected by cᵀNc = 1.

That is the exact pattern observed. This is a candidate, not a demonstration:
the pre-1988 code path no longer exists in the file, so it cannot be run and
compared. It is recorded because it is supported by the artifact's own history
and because it is the only hypothesis so far that accounts for every symptom.

## Status

**Closed on our side.** Two published curves disagree with ours, and our side
of both is now verified from first principles rather than merely cross-checked
between solvers:

* Figs. 6(a)–(c)'s five-quark curve at 2K = 21 — right total (1.5%), wrong
  distribution, with a node at x ≈ 0.26 that does not move across three states
  or two couplings.
* Thesis Fig. 18(b)'s seven-quark curve at 2K = 15 — uniformly 17.5% high.

Neither is explained by basis dependence, by numerical conditioning (the basis
is non-orthogonal only within a momentum configuration, and every one of those
Gram blocks has condition number below 42 with smallest eigenvalue 6.0), by Fock
level as such (B=2's 10-parton sector — the highest in any run here — is stable
to 0.14% and reproduces Fig. 6(d) to 0.8%, while B=1's *five*-parton sector at
2K=21 is the problem case), by the choice of K, or by any of the causes listed
below.

And the colour sums — the one mechanism that could have redistributed
higher-Fock weight while leaving masses, valence curves and sector totals
intact — are exact to the last bit, verified against explicit enumeration.

What remains is a statement about the figures, not about the code: these two
published curves are not reproducible by a correct calculation of the quantity
their captions name.

The superseded "resolution" is kept above as a caution.

| figure | sector | basis spread | agreement |
|---|---|---|---|
| Fig. 5 (meson, 2K=24) | 4q | 0.02% | 1.4–2.0% |
| Fig. 18(a) (baryon, 2K=15) | 5q | 2.95% | 0.4% |
| Fig. 18(b) (baryon, 2K=15) | 7q | 16.4% | 17.5% |
| Fig. 6(a)–(c) (baryon, 2K=21) | 5q | 34.8% | not reproducible |
| Fig. 6(d) (2-baryon, 2K=24) | 8q | 0.02% | 0.8% |

Everything determined to better than 3% is reproduced to better than 2%. The
two we cannot reproduce are the two whose own ambiguity exceeds 15%.

**One thing remains genuinely unexplained.** The basis spread accounts for the
*size* of the Figs. 6(a)–(c) disagreement but not the *shape*: the published
node sits at x ≈ 0.26 in all six published panels, across three states and two
couplings, and no computed variant puts a node anywhere. Whatever produced that
fixed feature is not something this reproduction has identified. Both earlier failures were in reading
the figures, not the physics:

* **Fig. 6(a)–(c)** — the valence and five-quark curves cross three times, and
  at scan resolution the dashes of the dashed curve read as markers. Validated
  instead against thesis Fig. 18(a), which prints that sector alone: 0.07–3.9%.
* **Fig. 6(d)** — solved at 2K = 22 against a figure drawn at 2K = 24. At the
  right K it agrees to 0.14–0.8%.

Regression tests for both are in `tests/test_paper.py`. Every mass, every
eigenvalue and every valence structure function was unaffected throughout.
