# Baryon higher-Fock structure functions

**Conclusions first. The investigation history, including two conclusions that
were published here and later withdrawn, is the appendix.**

Two published higher-Fock curves are not reproduced by either solver. Everything
else in the same family is, and our side of both exceptions is verified from
first principles rather than merely cross-checked between the two codes.

| figure | sector | 2K | agreement |
|---|---|---|---|
| Fig. 4(a) | meson qq̄qq̄ | 14 | **0.0–1.6%** |
| Fig. 5(a)–(d) | meson qq̄qq̄ | 24 | **0.9–2.0%** |
| Fig. 4(b) | baryon qqqqq̄ | 15 | **1.4–3.9%** |
| thesis Fig. 18(a) | baryon qqqqq̄ | 15 | **0.4%** |
| Fig. 6(d) | two-baryon 6q+qq̄ | 24 | **0.8%** |
| Fig. 4(c) / thesis Fig. 18(b) | baryon +2 pairs | 15 | ✗ 17.5% high |
| Figs. 6(a)–(c) | baryon qqqqq̄ | 21 | ✗ **incompatible with the same panels' own valence curves** |

The five-quark sector — the one Figs. 6(a)–(c) plot — is reproduced to 0.4–3.9%
at 2K = 15 by two independent published panels. The disagreement is specific to
2K = 21, not general to baryon higher-Fock.

For Figs. 6(a)–(c) the question is now settled in one direction and closed in
the other. **It is not our observable and not our reading of the figure**: the
published five-quark points exceed a bound our amplitudes cannot violate
(§2.1.2), and they are incompatible with the valence points printed in the same
panel (§2.1.3), in all three panels (§2.1.4). **What the 1988–90 code did
instead is not recoverable** from the surviving artifacts (§2.5).

---

## 1. What is established

### 1.1 Fig. 6(d) is 2K = 24, and it agrees

The caption gives 2K = 21 for panels (a)–(c) and **states no K at all for (d)**.
We assumed 22. It is 24, and the panel says so itself: wherever a structure
function vanishes the marker is still plotted and comes to rest on the axis, so a
two-baryon panel carries a long row of zero-valued markers, one per lattice site.
They fall at

    0.0427  0.1304  0.2095  0.2937  0.3763  0.4569  0.5427

which is k/24 for odd k to within 0.006; at 2K = 22 the predicted positions drift
by a full site before the right-hand edge. See `docs/inferred-K.md`.

Recomputed at 2K = 24 — Fortran and Python agreeing to every printed digit:

| x | series | published | ours | dev |
|---|---|---|---|---|
| 1/24 | 8-parton q ×5·10² | 18.448 | 18.296 | 0.8% |
| 1/24 | 8-parton q̄ ×10³ | 16.219 | 16.133 | 0.5% |
| 3/24 | 6-parton valence | 35.454 | 35.505 | **0.14%** |
| 3/24 | 8-parton q ×5·10² | 19.279 | 19.231 | 0.25% |
| 5/24 | 6-parton valence | 35.541 | 35.448 | 0.26% |

The panel's three curves come from **two** Fock sectors: the 6-parton valence,
plus the 8-parton sector's quark *and* antiquark distributions. A third sector
does open at 2K = 24 — ten partons of odd momentum, with at most N = 3 quarks per
momentum, need total momentum ≥ 24, so the 10-parton sector exists at 24 and is
forbidden at 22. All three are present and the number sum rule closes exactly:

    6-parton   5.990892
    8-parton   0.010622 (q)   0.001517 (q̄)
    10-parton  0.000004 (q)   0.000001 (q̄)
    ------------------------------------------
    total q − q̄ = 6.000000 = N·B

The thesis prints this panel as its Fig. 13(a) on a different y scale — labelled
0 10 20 30 where the article uses 0 12 24 36 48 — so the two published versions
differ by a constant 48/30 = 1.6 everywhere (measured: 1.603 ± 0.009 over six
points). The number sum rule picks the article using published numbers only:
∫q dx must be N·B = 6, which the article's scale gives (5.9) and the thesis's
misses (3.7). Thesis Fig. 13(b)–(e) show the 2nd–5th B = 2 states and remain
available as further regression targets.

### 1.2 The five-quark sector is correct at 2K = 15

Thesis Fig. 18(a) plots it alone, with q and q̄ separated and no curve crossings
— the cleanest published target for this sector:

| quantity | ours | thesis | ratio |
|---|---|---|---|
| q(1/15) ×10³ | 7.933 | 7.93 | **1.000** |
| q(3/15) ×10³ | 4.975 | 4.77 | 1.043 |
| q(5/15) ×10³ | 6.239 | 6.74 | 0.926 |
| q(7/15) ×10³ | 2.764 | 2.72 | 1.016 |
| q̄(1/15) ×10³ | 4.710 | 4.72 | **0.998** |

and at m/g = 0.1 the implied multiplier is ×95.6 against a legend of ×10², i.e.
4.6%. The article's own Fig. 4(b) is the same data and gives 1.4–3.9%
independently.

### 1.3 `structure_function` is a reconstruction of the lost `wfbig`

**The x-space conversion in the original work was done by a separate program,
`wf`/`wfbig`, which is not in this repository and appears not to survive.**
`qcdf.f` emits only eigenvectors, Fock content and the basis-change matrix, and
refers to `wf` in comments. Every structure function in the paper passed through
that program.

`dlcq.observables.structure_function` is a reconstruction of it. That is worth
stating plainly rather than as a caveat: it is the one piece of the original
pipeline that had to be rebuilt from the physics rather than ported, and it is
now better validated than the code it replaces.

It cannot be checked the way the rest of the port was. The Python solver
reproduces the Fortran's matrix elements to 10⁻¹⁴, so any error in the
conversion would be shared by both and invisible to a solver-vs-solver
comparison. Five independent routes were used instead:

1. **The definition is the thesis's.** Sec. 2.4 gives
   `q_k = <phi(K)| b^dag_{k,c} b_{k,c} |phi(K)>`, the baryon-number sum rule
   `sum_k (q_k − qbar_k)` = 0 for mesons and N for baryons, the momentum sum
   rule `sum_k k (q_k + qbar_k) = K`, and the continuum limit `q(x) = K q_k` —
   which is what `structure_function` implements, factor of K included. Ref. 17
   (Lepage & Brodsky, Phys. Rev. D **22**, 2157) gives the same
   probability-weighted parton count in standard light-cone form.
2. **The weighting is provably the exact expectation value.** The norm couples
   only basis states sharing a momentum configuration (§1.4), so the number
   operator acts as a scalar within every coupled block and
   `<phi|n_k|phi>` collapses to `sum_s n_k(s) c_s (Nc)_s` with `c^T N c = 1`.
   The `c_s (Nc)_s` weight is not a convention; `c_s^2` would silently break
   every sum rule.
3. **The norm it uses is exact** against brute-force colour enumeration (§1.4).
4. **Both sum rules hold to machine precision**, across four runs and both
   solvers:

   | run | \|momentum − 1\| | \|number − N·B\| |
   |---|---|---|
   | N=3 B=0 2K=24 meson | 0.00e+00 | 3.4e−16 |
   | N=3 B=1 2K=21 baryon | 2.2e−16 | 0.00e+00 |
   | N=3 B=1 2K=15 baryon | 2.2e−16 | 4.4e−16 |
   | N=3 B=2 2K=24 two-baryon | 2.2e−16 | 0.00e+00 |
   | the same, read from Fortran output | 3.6e−15 | 1.1e−14 |

5. **It reproduces the published curves** wherever they can be read cleanly:
   0.0–2.0% for every meson sector, 0.4–3.9% for the baryon five-quark sector at
   2K = 15, 0.8% for the two-baryon eight-parton sector, and 1–5% against the
   thesis's Table 6 higher-Fock probabilities — the only place the thesis quotes
   a higher-Fock magnitude with uncertainties.

The two curves in §2 are the exceptions, and the same five checks are what make
it hard to attribute them to the reconstruction.

### 1.4 The observable is exact by construction

Three facts, each tested rather than assumed:

* **The norm equals the true Gram matrix.** For N = 3 the colour sum is small
  enough to do by brute force. `tools/colour_norm.py` expands each basis state
  over the (type, momentum, flavour, colour) Fock basis straight from the
  definition — N! signed terms per ε-contracted cluster, N per δ-contracted
  meson, with the fermionic sign of sorting into canonical order — and takes
  inner products. It reproduces the solver's norm **exactly** at every K tested,
  overall factor 1, maximum relative difference 0. At 2K = 21: 3q (12 states),
  5q (133), 7q (248), 9q (109), all exact.
* **The norm couples only states sharing a momentum configuration** — 184 such
  pairs at 2K = 21, and zero others. That is precisely the condition making
  q(k) = K Σₛ nₖ(s) cₛ(Nc)ₛ the exact expectation value.
* **The eigenvector is determined to 10⁻¹²** across well-posed treatments of the
  non-orthogonal basis.

So the five-quark structure function is correct by construction, not by
convention.

### 1.5 The magnitudes were independently confirmed

The thesis carries something the article does not: numbers, with uncertainties,
for a higher-Fock *magnitude*. Its Table 6 gives the probability of the
four-quark component in the lightest SU(2) meson, Richardson-extrapolated over
2K = 10–24. Our sweep matches it to 1–5% at four couplings
(`refs/thesis_table6.csv`).

---

## 2. What remains open

### 2.1 Figs. 6(a)–(c), five-quark at 2K = 21

Marker for marker, at the published lattice. The published row is the mean of
**three independent extractions** of thesis p. 82 at 600 dpi — the column trace
in `refs/thesis_fig12a_fivequark.csv`, erosion-centroid blob finding, and raw
ink runs per lattice column with no morphology at all — with their spread as
the uncertainty. `tools/fig6a_consistency.py` carries them and every test below:

| k | 1 | 3 | 5 | 7 | 9 | 11 | sum |
|---|---|---|---|---|---|---|---|
| published | 8.57(4) | 8.01(5) | 1.02(5) | 1.77(6) | 5.04(0) | 0.52(5) | 24.9 |
| ours | 6.30 | 5.68 | 3.86 | 4.77 | 3.49 | 1.14 | 25.2 |

**Reading the panel is a 1% operation and the disagreement is 42%.** The three
methods agree with each other to 0.6–1.3% in the same best-scale metric that
puts each of them 42.2–42.4% from ours.

The vertical scale is not in question either. It is pinned by the label
positions, per §3, at 15.0; the frames of panels (b) and (c) then measure 6.00
major-tick intervals, confirming it. (Panel (a)'s own tick scan does *not*
confirm it — it returns 173.5, 175.5, 175.0, 32.5, 35.0, 106.5 px, the curve-ink
contamination §3 warns about, which is why label positions and not ticks are
the primary calibration.) The decisive check is that the same calibration lands
the *valence* markers on ours to 0.6%, and §2.1.1 removes the axis from the
argument altogether.

**The sector total agrees to 1.5%; only its distribution over x differs.** That
is the constraint any explanation has to satisfy, and it is what rules out every
mechanism in §2.3. Fitting a free scale — so the paper's own ×10³ never enters —
returns 1055 against a stated 10³, confirming the normalization independently.

#### 2.1.1 The disagreement survives with no axis at all

Both curves are printed on the same axis, so the five-quark/valence ratio
*within* a panel is independent of how that axis is pinned — every error in the
calibration cancels exactly:

| k | 3 | 5 | 7 | 9 | 11 |
|---|---|---|---|---|---|
| published 5q/valence ×10⁻³ | 4.042 | 0.128 | 0.148 | 0.717 | 0.258 |
| ours | 2.877 | 0.466 | 0.412 | 0.494 | 0.538 |
| **ratio** | **1.41** | **0.28** | **0.36** | **1.45** | **0.48** |

while the valence itself agrees to 1.6% rms. No axis-pinning error of any size
produces that table.

#### 2.1.2 It is not the observable: a bound our amplitudes cannot violate

`q(k)` can never exceed the sector's **total** parton density `q(k) + q̄(k)`.
Relabelling which parton is the antiquark, plotting one colour cluster,
partitioning the sector any way at all — every such choice only moves weight
between `q` and `q̄`, so none can push `q` above the sum. Against our
amplitudes:

| k | 1 | 3 | 5 | 7 | 9 | 11 |
|---|---|---|---|---|---|---|
| published q | 8.57 | 8.01 | 1.02 | 1.77 | 5.04 | 0.52 |
| our q + q̄ | 10.72 | **7.27** | 4.16 | 4.81 | **3.49** | 1.14 |
| ratio | 0.80 | **1.10** | 0.25 | 0.37 | **1.44** | 0.45 |

At k = 3 and k = 9 the published quark density exceeds our entire five-parton
parton density. **No re-reading of our amplitudes can produce the published
points**; the amplitudes themselves would have to differ. This is a bound, not
a search — it needs the ×10³ the legend states, and that is legible in panel (a)
and corroborated by the sector total.

Measured directly rather than argued: neither colour cluster alone matches
either. Splitting the five-parton states into their (qqq) and (qq̄) clusters —
verified state by state against `mstinf`'s own `kmes`/`kbrb` — gives free-scale
residuals of 40.5% (qqq cluster), 90.1% (sea quark) and 54.6% (antiquark)
against 42.4% for the full distribution. Nothing beats plotting all four quarks,
which is what `structure_function` does. (The sea-quark/cluster split is
partly an enumeration convention — no basis state carries its meson-cluster
quark at k = 1 — so read those rows only as the elimination they are.)

#### 2.1.3 The five-quark points are a function of the valence points

In the original Fock basis `H c = w N c` with `N` block-diagonal in parton
number, so the five-parton row determines `c₅` from `c₃` outright:

    c₅ = −(H₅₅ − w N₅₅)⁻¹ H₅₃ c₃

Measured against the full eigenvector this is exact to **5.8 × 10⁻⁴** (2K = 21
ground state), 2.1 × 10⁻³ and 2.7 × 10⁻³ for the two excited states, and
5.5 × 10⁻⁴ at 2K = 15. And the propagator hardly enters: replacing
`(H₅₅ − wN₅₅)⁻¹` by its diagonal, or by free masses, moves the residual from
42.41% to 44.24% and 42.57%. **The five-quark curve is fixed by the valence
curve**, with no freedom worth the name.

So the two curves in the panel constrain each other. Optimizing over all twelve
valence amplitudes *and* a free `w` *and* a free overall scale on the five-quark
points — every benefit of the doubt to the published side — over 1050 local
optima:

| valence held to | best five-quark residual found |
|---|---|
| ±1% | 23.9% |
| ±2% | 18.5% |
| ±3% | 15.7% |
| ±20% | 7.9% |
| ±100% | 0.0% |

Our actual ground state sits at valence 0.61%, five-quark 42.4%. Read against a
1% reading error: **with the valence pinned where it demonstrably agrees, no
wavefunction was found that fits the five-quark points to better than 16–24%**,
and driving that residual to zero costs a valence that is ~100% wrong. The
published five-quark points *are* a legitimate five-parton distribution of some
state — they are not kinematically impossible — but not of the state whose
valence the same panel prints.

Two honest limits on this one. It is a search, so those are upper bounds on the
true minima, not proofs. And the ±3% optimum partly exploits the k = 1 valence
point, left unconstrained because that marker sits on the axis and is the one
valence point that cannot be read.

#### 2.1.4 All three panels, not one

The same test on panels (b) and (c) — the 2nd and 3rd baryons, same K, same
basis, same extraction. Markers re-read from thesis p. 82; the five-quark
residual carries a free scale, because the ×10ⁿ in those two legends is **not**
legible in the scan, which also means the §2.1.2 bound does not apply to them:

| panel | state | valence residual | five-quark residual |
|---|---|---|---|
| (a) | 1st | **0.61%** | **42.4%** |
| (b) | 2nd | **1.69%** | **73.3%** |
| (c) | 3rd | **2.93%** | **68.9%** |

and the deviation has the same sign structure in all three — suppressed at
k = 5, 7 and enhanced at k = 1, 9:

| panel | k=1 | 3 | 5 | 7 | 9 | 11 |
|---|---|---|---|---|---|---|
| (a) | 1.29 | 1.34 | 0.25 | 0.35 | 1.37 | 0.43 |
| (b) | 1.92 | — | — | 0.14 | 4.24 | — |
| (c) | 2.11 | — | 0.25 | 0.29 | 2.64 | — |

(dashes are markers merged with the valence curve and not separable.) The
valence peaks at k = 7 in all three states, so **the suppression sits exactly
under the valence peak and the enhancement away from it.** That is a property of
the figure, not of one panel, and it is the sharpest surviving clue to whatever
produced it.

The apparent "node" at x ≈ 0.26 is **not** a data feature. The thesis says of
exactly these panels:

> The connecting curves are simply cubic spline fits: the resolution is in some
> cases not good enough for these to accurately depict the actual structure.

It is spline undershoot between the two smallest markers, which in every one of
these panels are those at k = 5 and k = 7. At 3× magnification both markers are
plainly there, at 1.14 and 1.71, with the spline dipping to 0.31 between them.

#### Fock truncation cannot be the mechanism

Re-verified against the current code, and the result is sharper than "LPN makes
no difference". The five-quark distribution at 2K = 21 is **identical to four
decimals** across truncations that change the basis by a factor of two:

| LPN | states retained | q(x) x10^3 at k = 1..11 | corr. with published |
|---|---|---|---|
| 0 (none) | 193 | 6.30 5.67 3.86 4.77 3.48 1.14 | +0.7661 |
| 5 | **95** | 6.30 5.67 3.86 4.77 3.48 1.14 | +0.7662 |
| 7 | 175 | 6.30 5.67 3.86 4.77 3.48 1.14 | +0.7661 |

Halving the basis moves nothing. The five-quark sector is effectively decoupled
from everything above it -- unsurprising given the sector probabilities
(2.4x10^-3 at five partons, 2.1x10^-7 at seven), but worth pinning, because a
truncation in the original run is otherwise the most natural explanation for a
panel that disagrees at one K and not another.

#### The published 2K = 21 panel is inconsistent with the paper's own 2K = 15 panel

The five-quark quark distribution is published twice, at two different K, and we
reproduce one of them. Tracking how the shape evolves with K settles which of
the two is the outlier.

Our q(x) for this sector, ×10³, ground state, m/g = 1.6:

| 2K | k=1 | 3 | 5 | 7 | 9 | 11 | 13 |
|---|---|---|---|---|---|---|---|
| 15 | 7.93 | 4.97 | 6.24 | 2.76 | 0.27 | — | — |
| 17 | 7.28 | 5.26 | 5.11 | 4.46 | 1.08 | 0.08 | — |
| 19 | 6.75 | 5.51 | 4.26 | 5.01 | 2.36 | 0.43 | 0.03 |
| 21 | 6.30 | 5.67 | 3.86 | 4.77 | 3.48 | 1.14 | 0.18 |
| 23 | 5.91 | 5.77 | 3.76 | 4.28 | 4.10 | 2.06 | 0.55 |
| 25 | 5.57 | 5.82 | 3.83 | 3.82 | 4.25 | 2.90 | 1.15 |

The distribution is two-lobed at every K, with a shallow interior minimum that
drifts right as the lattice refines but stays near x ≈ 0.22–0.28. Its **depth is
stable**: the minimum sits at 55–63% of the peak throughout.

Against the two published panels:

| 2K | source | peak | interior minimum | min/peak |
|---|---|---|---|---|
| 15 | ours | 7.93 | 4.97 | 0.63 |
| 15 | **published** (thesis Fig. 18(a)) | 7.93 | 4.77 | **0.60** |
| 21 | ours | 6.30 | 3.48 | 0.55 |
| 21 | **published** (Fig. 6(a)) | 8.56 | 1.14 | **0.13** |

At 2K = 15 the two agree on the depth to 5%. At 2K = 21 they differ by a factor
of four, and the published minimum is deeper than anything our calculation
produces at any K in 15–25.

So the disagreement is not that our five-quark distribution has the wrong shape
— at 2K = 15 it has the published shape, marker for marker. It is that the
published shape **changes qualitatively between the paper's own two panels**,
from a shallow 60% minimum to a deep 13% one, while ours evolves smoothly. Note
also that the apparent node in Fig. 6(a) is *not* the spline undershoot
discussed above: the markers themselves, at 1.14 and 1.71, are what make it
deep.

This does not identify the cause, and it is not proof that the 2K = 21 panel is
wrong — only that one of the two published panels is inconsistent with a smooth
K dependence, and it is not the one we fail to reproduce. Combined with §1.4
(colour sums exact against brute-force enumeration, norm exact, eigenvector
determined to 10⁻¹², both sum rules at machine precision) and with
`docs/basis-dependence.md`'s finding that the seven-parton sector is
basis-independent to 2×10⁻¹², an explanation on our side would have to be
K-specific and leave every one of those checks intact.

### 2.2 Fig. 4(c) / thesis Fig. 18(b), seven-quark at 2K = 15

| x | published ×10⁷ | ours ×10⁷ | ratio |
|---|---|---|---|
| 1/15 | 5.32 | 6.232 | 1.171 |
| 3/15 | 3.96 | 4.806 | 1.214 |
| 5/15 | 2.73 | 3.116 | 1.141 |

Systematically 17.5% high, in a sector carrying 1.9 × 10⁻⁷ of the baryon's quark
number — four orders below the five-quark sector, which reproduces the panel
directly above it on the same page to 0.4%.

### 2.3 Everything eliminated

| candidate | verdict |
|---|---|
| the K of the panel | 2K = 21 confirmed by the caption *and* by the valence matching to 0.4–2.0% |
| the y-axis calibration | fitted from label positions; Fig. 6's top is 14.77, and a 1% shift leaves a 40% disagreement |
| two curves on different y-scales | no right-hand axis exists; the panel row was checked out to its margins |
| marker misassignment | fixed by classifying on ink at the marker's *centre*; two valence points had been called higher-Fock |
| the spline "node" | the author's own caveat, §2.1 |
| the state / degeneracy | ground state isolated by a gap of 3.44 in M²; its valence matches to 0.8% median, every other state ≥ 29% off |
| a different quantity | q̄, q ± q̄ and every Fock sector of the first twelve states — ~100 candidates — top out at correlation 0.88 |
| a different K | 2K = 11…25 interpolated onto the 21-lattice tops out at 0.81, rms 50% |
| the cluster partition | plotting one cluster of a baryon+meson state would change the sector integral, which already agrees to 1.5% |
| basis treatment | well-posed variants agree to 0.003% at 2K = 21 and 2 × 10⁻¹² at 2K = 15 |
| numerical conditioning | every momentum-configuration Gram block has condition number ≤ 42 and smallest eigenvalue 6.0 |
| Fock level as such | B = 2's 10-parton sector — the highest in any run here — is stable to 0.14% and reproduces Fig. 6(d) to 0.8% |
| basis truncation (LPN) | cutting at 5, 7, 9 partons leaves the amplitudes identical to three decimals |
| the non-orthogonal weighting | c²·diag(N), c², \|c(Nc)\|, (Nc)² all correlate *worse* than c·(Nc) |
| the colour sums | exact against brute-force enumeration, §1.4 |
| K-dependence / fragility | P₅ runs 3.19, 2.96, 2.74, 2.56, 2.42, 2.31 (×10⁻³) across 2K = 13…23 — smooth and monotonic |
| **numerical damage** | §2.4 |
| the y-axis calibration, again | the in-panel five-quark/valence ratio needs no axis at all and still disagrees by 0.28–1.45× (§2.1.1) |
| the marker reading | three independent extractions agree to 0.6–1.3%, against a 42% disagreement (§2.1) |
| the colour cluster plotted | neither cluster alone, nor the antiquark, beats the full distribution: 40.5 / 90.1 / 54.6% vs 42.4% (§2.1.2) |
| the sign structure of c₃ | **all 2048** sign patterns of the twelve valence amplitudes scanned; best 34.8%, unflipped 42.3%. Worth doing precisely because the valence curve is blind to these signs and the five-quark curve is not — it was the last free thing in c₃ |
| the five-quark propagator | exact / diagonal-only / free-mass denominators give 42.3 / 44.1 / 42.4%. The sector's own dynamics barely enter (§2.1.3) |
| a truncated state list | a fixed-array overflow in the 1988 code was the best mechanical candidate left — K-specific, and it would leave the valence untouched. Every prefix of the generation order solved; best 31.8% at 61 of 193 states. It fails |
| the coupling | every cached m/g from 0.05 to 1.6: 36.6–42.3%, none reproducing the shape |
| **any wavefunction at all** | with the valence pinned to its own published points, nothing found below 16–24% (§2.1.3) |

### 2.4 Numerical damage cannot produce it

`qcdf.f` is genuinely vulnerable to array overflow — it dimensions `IDELT` with
25 colour slots while an element between L-parton states needs 2L + 4, so runs
reaching 12 partons really are corrupted (`docs/fortran-color-overflow.md`) — and
part of it was single precision before the `**** modified for double precision
****` change. So the suspicion is fair. `tools/robustness.py` tests it:

| probe | M² | valence pk | 5q total | 5q distribution |
|---|---|---|---|---|
| arithmetic in float32 | exact | exact | exact | **0.00%** |
| arithmetic in float16 | 3×10⁻⁴ | 0.02% | 0.05% | **<0.9%** |
| weeding eps 10⁻²…10⁻¹⁰ (189–193 states survive) | exact | exact | exact | **0.00%** |
| L ≥ 7 elements zeroed / doubled / sign-flipped | <10⁻⁶ | <0.01% | <0.06% | **<0.2%** |
| L ≥ 9 elements zeroed / doubled / sign-flipped | exact | exact | exact | **0.00%** |
| L ≥ 5 elements zeroed | 0.2% | 0.4% | **25.2 → 0.0** | — |
| L ≥ 5 elements doubled | 0.2% | 0.3% | **25.2 → 19.0** | — |
| L ≥ 5 elements sign-flipped | 0.3% | 0.7% | **25.2 → 9.6** | — |

**Every perturbation that leaves the sector total intact leaves the distribution
intact; the only ones that move the distribution destroy the total.** The
published curve has the right total. So no damage of this kind produced it.

At 2K = 21 nothing overflows in any case: the longest state is 9 partons, needing
22 of the 25 available slots.

### 2.5 The most plausible remaining account

`fortran/qcdf.f` postdates the figures. Its own header records:

```
**** 6/24/88 MODIFIED QCD2A2 SO THAT COLOR ****
**** SUMS ARE PERFORMED DIAGRAMMATICALLY   ****
**** RATHER THAN ITERATIVELY               ****

**** 5/8/90 modified qcdf to include arbitrary number of quark flavors ****

(C) Kent Hornbostel 1993.  All Rights Reserved.
```

The thesis is 1988, the article 1990, the file a 1993 snapshot. The colour sums
set the higher-Fock weighting, and they were rewritten mid-thesis.

This is *not* a demonstration, and §1.4 shows it does not implicate the code we
have: the current colour sums are exact, so if the pre-1988 path differed it was
the one in error. That path no longer exists in the file and cannot be run.
Short of finding an intermediate version of `qcdf.f`, this is as far as the
artifacts allow. The lost `wf`/`wfbig` conversion program (§1.3) is the other
candidate and is equally unavailable.

What §2.1 adds is a specification any such account now has to meet. It must
act on the five-parton **amplitudes** — the observable and every re-partition of
it are excluded by the bound in §2.1.2 — while leaving the valence amplitudes
alone to better than 3%, since the same panels' valence curves reproduce. It
must be inactive at 2K = 15, where the same sector reproduces and the two
published curves *are* mutually consistent. And it must suppress k = 5, 7 and
enhance k = 1, 9 across three different eigenstates (§2.1.4). Nothing in the
surviving code does that, which is why this section stays an account rather
than a finding.

---

## 3. Method notes worth keeping

**Prefer a panel that plots one series without crossings**, even when the
composite is the one the article printed. Figs. 6(a)–(c) overlay the valence and
five-quark curves, which cross three times, and at scan resolution the *dashes*
of the dashed curve are the same size as its markers — which is what made this
sector look irreproducible for a long time. Thesis Figs. 18(a)/(b) plot the same
sectors alone.

**Fit the vertical scale from the label positions**, never from tick spacing
(contaminated by curve ink) and never by assuming the frame top is the last
label. See `docs/figure-validation.md`.

**Classify filled vs open markers on the ink at the marker's centre**, not on how
many runs a vertical slice breaks into: a curve entering a ring fills its
interior, so the run count fails exactly where curves cross.

**Extract each panel three ways and quote the spread.** Column trace,
erosion-centroid blob finding, and raw ink runs per lattice column fail
differently — the trace picks up the wrong curve where two cross, erosion drops
markers that merge with a steep curve, ink runs merge a marker with the curve
passing through it — so their agreement is real evidence and their spread is an
honest error bar. On Fig. 6(a) they agree to 0.6–1.3% against a 42%
disagreement with the calculation, which is what makes that 42% a finding
rather than a digitization artifact (§2.1).

**Read the ink-run widths, not just the positions.** A marker core is 20–45 px
at 600 dpi and a curve crossing a column is 4–8 px, so the width alone
separates them, and a run far outside that range is a merge to be flagged
rather than a datum. Listing every run in a lattice column — rather than only
the ones a detector accepts — is also what lets you say where a marker is *not*:
in panel (b) the five-quark markers at k = 3 and 5 could not be separated (at
k = 3 the dashed lobe crosses the valence marker; at k = 5 the curve has come
down onto the axis), but the column shows no ink at all where our calculation
puts them, 0.93 and 1.36 in that panel's units. Those two were left out of
§2.1.4 as unmeasured rather than counted as agreeing or disagreeing.

---

# Appendix: investigation history

Kept because two conclusions were published in this file and later withdrawn. A
reader who saw them deserves the correction rather than a silent edit.

## A.1 Withdrawn: "the five-quark distribution is basis-dependent"

An earlier revision reported a 34% "basis-treatment spread", concluded that it
explained the disagreement, and presented a table in which the size of the spread
tracked the published agreement across every sector. **That was wrong.**

The spread was measured over all six assembly/policy combinations, two of which
apply the diagonal-only assembly *without* the block-wise repair that
`docs/basis-dependence.md` shows is required — H₀ = D·N with D block-constant, so
reading off only the diagonal does not determine it. Those two are ill posed, not
alternative conventions, so the spread was measuring the historical code's error
rather than any ambiguity in the physics.

Restricted to the four well-posed treatments the sectors agree to 8 × 10⁻¹³ at
2K = 15 and 0.003% at 2K = 21.

What the ill-posed variants do show is worth keeping, and is why the error was
easy to make: they shift the 2K = 21 five-quark curve by ~30% while moving M²
only in its sixth digit — exactly what a genuine physical ambiguity would look
like.

| variant | M² | valence pk | 5q ×10³ at k = 1…11 |
|---|---|---|---|
| historical `qcdf.out` (1990) | 10.390380000 | 11.6025 | 7.36 4.50 3.59 4.66 3.64 1.45 |
| exact / fortran | 10.390814172 | 11.6047 | 6.30 5.67 3.86 4.77 3.48 1.14 |
| exact / blockwise | 10.390814172 | 11.6047 | 6.30 5.67 3.86 4.77 3.48 1.14 |
| exact / spectral | 10.390814080 | 11.6047 | 6.30 5.67 3.86 4.77 3.48 1.14 |
| *fortran / fortran* | *10.389561840* | *11.5968* | *8.06 4.91 3.68 5.06 3.85 1.57* |
| fortran / blockwise | 10.390814172 | 11.6047 | 6.30 5.67 3.86 4.77 3.48 1.14 |
| *fortran / spectral* | *10.389158840* | *11.6036* | *7.69 5.69 4.08 5.20 4.05 1.59* |

Italicised rows are the ill-posed ones.

An earlier version of this section also argued that the higher-Fock sectors are
"numerically fragile", from the ~16% spread between our own two codes at 2K = 21.
That spread is the same artifact: it is the gap between the historical
diagonal-only run and the repaired one, not an intrinsic sensitivity. §2.4 shows
the quantity is in fact robust to half-precision arithmetic.

## A.2 Withdrawn: "the seven-quark discrepancy is within its own ambiguity"

Followed from A.1 and falls with it. The 16.4% "spread" was the same artifact;
among well-posed treatments the seven-quark sector at 2K = 15 agrees to
2 × 10⁻¹⁰. The 17.5% discrepancy is real (§2.2).

## A.3 Over-read: "the node does not move"

An earlier revision made much of the published node sitting at x ≈ 0.26 across
three states and two couplings, and treated its refusal to move as a physical
signature no computed variant reproduced. The thesis's own caveat (§2.1) shows it
is spline undershoot between the two smallest markers, which sit at k = 5 and
k = 7 in all of those panels. Not a data feature.

## A.4 Where the trail started

The original symptom: Fig. 6(a)'s higher-Fock series correlated with our curve at
only 0.77, where the meson equivalents managed 0.999 *and* recovered the paper's
×10ⁿ multipliers to ~1%.

| series | corr | fitted scale | paper's stated scale |
|---|---|---|---|
| meson Fig. 5(a) valence | 0.9999 | 0.996 | ×1 |
| meson Fig. 5(a) higher-Fock | 0.9994 | 986.3 | ×10³ |
| meson Fig. 5(c) higher-Fock | 0.9988 | 99.17 | ×10² |
| baryon Fig. 6(a) valence | 0.9999 | 0.993 | ×1 |
| **baryon Fig. 6(a) higher-Fock** | **0.7729** | 1046 | ×10³ |

That the fitted scale was right while the shape was not is the same signature
that survives today: total right, distribution wrong.

Two things confused matters for a long time and were resolved along the way.
Fig. 6(d) was being solved at the wrong K (§1.1), which made a second, unrelated
panel look broken. And the five-quark sector was being judged against
Figs. 6(a)–(c) alone, before thesis Figs. 18(a)/(b) were found — those plot the
same sectors with no crossings and settled that the sector itself is right at
2K = 15.

Two real tracer bugs were also found and fixed while chasing this: legend text
counters being read as open circles, and legend markers landing on lattice sites
surviving the lattice filter. Fixing them moved the Fig. 6(a) correlation from
−0.056 to 0.773.
