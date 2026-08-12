# The baryon higher-Fock discrepancy

Both codes reproduce every valence curve and the *meson* higher-Fock curves to
~1%. The **baryon** higher-Fock series in Fig. 6 does not match. This is what
the investigation established, and what it ruled out.

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

## The remaining hypothesis

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

## Status

Reported, not resolved. The affected curves are the higher-Fock series of
Fig. 6 (and by extension Fig. 4's baryon panels). Every mass, every eigenvalue,
every valence structure function and all meson higher-Fock results are
unaffected — the baryon higher-Fock sector carries ~0.06% of the ground-state
wave function.
