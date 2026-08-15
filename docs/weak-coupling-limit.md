# The weak-coupling limit is a grid artifact, not the continuum

**Summary.** At `m/g ≲ 0.2` the DLCQ momentum grid cannot resolve the region
where the physics lives. The wavefunction goes as `φ ~ x^a` at the endpoints
with `a → 0` in the chiral limit, so the endpoint region is `x ≲ e^{-1/a}` —
`5 × 10⁻¹¹` at `m/g = 0.05`, `9 × 10⁻⁴²` at `m/g = 0.0125` — against a grid
whose smallest momentum fraction is `x = 1/K ≈ 1/35`. Everything anomalous
about the weak-coupling column follows from that one fact: the slow K
convergence, the fit-order ambiguity that dominates the error budget, and a
chiral scaling `M² ∝ m²` where the physical answer is `M² ∝ m`.

This is **not** a defect of `qcdf.f` or of the port. It is a property of
uniform-grid methods applied to this problem, and it is reproduced here in a
*continuum* solver that shares no code with the DLCQ path. It does not affect
the paper's figures or its strong-coupling results. It does mean the
weak-coupling entries of Table I converge to something other than the physical
continuum limit, and that the uncertainty quoted on them cannot be reduced by
better fitting.

## What is measured

The increments of `M²(K)` are a clean power law, and the power is not the one
Eq. (27) assumes. Fitting `dM²/dK ~ K^-p` locally and extrapolating the local
exponent as `p(K) = p_∞ + c/K`:

| m/g | `a` | `1+a` | `p_∞` | closer to |
|---|---|---|---|---|
| 1.60 | 0.845 | 1.845 | 1.974 | **2** |
| 0.80 | 0.586 | 1.586 | 1.893 | **2** |
| 0.40 | 0.326 | 1.326 | 1.554 | 1+a |
| 0.20 | 0.168 | 1.168 | 1.215 | **1+a** |
| 0.10 | 0.084 | 1.084 | 1.067 | **1+a** |
| 0.05 | 0.042 | 1.042 | 1.017 | **1+a** |

(N = 3 baryon; the other four channels agree.) At strong coupling the
increments decay as `K^-2`, which is what Eq. (27)'s `1, K⁻¹, K^-(1+a), K⁻², …`
basis is built for. At weak coupling they decay as `K^-(1+a)`, matching
Eq. (26)'s `a` to 0.003–0.047 across all thirty Table I entries — so the
remainder goes as `K^-a`, and with `a = 0.036` that factor moves **3.6% across
the entire window**.

The gap being extrapolated across is correspondingly large. At `m/g ≤ 0.2`,
**15–35% of `M(0)` lies beyond 2K = 70** — the extrapolation adds 17–54% to the
last computed point. For mes N = 4 at `m/g = 0.05`, `M²(2K=70) = 0.08554`
against `M(0) = 0.1226`.

## The cause, and the control that proves it

Taking the `K^-a` remainder literally gives limits 3–4× the published values
(see "what does not work" below), so it cannot be the whole story. The
resolution is that `p_∞ = 1+a` is a *local* exponent that has not reached its
asymptote — the true asymptotic regime sits at `K ~ e^{1/a}`.

The control is `dlcq/thooft.py`, which solves the continuum 't Hooft equation
on a uniform grid in `x`. It shares no code with the DLCQ path, has no `K`, no
Fock-space truncation, and no boundary conditions. **It fails the same way.**
At fixed `m/g = 0.05`, refining its grid:

```
n =  200   1/n = 5.0e-3   (M/g)² = 0.03585
n =  400         2.5e-3            0.03912
n =  800         1.2e-3            0.04234
n = 1600         6.3e-4            0.04550
n = 3200         3.1e-4            0.04861
```

~7% per doubling, no plateau — the DLCQ `M²(K)` pathology, reproduced in the
continuum. Both methods are blind to `x ≲ e^{-1/a}`, and neither is converging
to the physical answer at accessible resolution.

## Consequence 1 — the chiral exponent is wrong, and it is not zero modes

`M²` in Table I's units scales as `(m/g)^α` with **α → 2**, measured at fixed K
and essentially K-independent (1.9537, 1.9497, 1.9450 at 2K = 35, 49, 69 for
the N = 3 baryon; `tools/joint_fit.py` reports 1.9968 at `m/g = 0.0125`).

The physical answer is `α = 1`. The 't Hooft model obeys a GMOR relation,
`M_π² ∝ m`, and the mechanism is visible in the equation: the `m²` term
contributes `m² ∫dx φ²(1/x + 1/(1−x))`, log-divergent at the endpoints and cut
off by the endpoint behaviour itself, `∫₀¹ x^{2a−1}dx = 1/(2a)`. So
`M² ~ m²/a ~ m·g`, since `a ∝ m/g`.

A grid that cannot resolve the endpoint never sees the `1/a` enhancement and
returns `m²` instead. The continuum control above returns α = 2 as well
(1.90, 1.97, 1.993, 1.998 down to `m/g = 0.00625`), which rules out the
obvious suspect: **this is not the missing zero-mode sector.** DLCQ here uses
antiperiodic boundary conditions and so has no fermion zero modes
(`python/qcdf_states._filter_parity_group`), and that is a genuine limitation
for vacuum structure — but it is not the cause of this, because a solver with
no such truncation reproduces the same exponent.

### Eq. (26) is not the source either

The endpoint relation contributes no non-analyticity. Expanding
`1 − πa·cot(πa) = (πa)²/3 + (πa)⁴/45 + … = μ` with `μ = π(m/g)²/C` and
`C = (N²−1)/2N`, and inverting:

```
a = (1/π)·√(3π/C)·(m/g)·[1 − π(m/g)²/(10C) + …]
```

an odd analytic series in `m/g`, verified against `units.endpoint_exponent` to
`1e-10` for `m/g ≤ 0.05`. The leading constants are 1.128379 (N=2), 0.846284
(N=3), 0.713650 (N=4).

## Consequence 2 — the error budget's `form` term is irreducible

`richardson_budget`'s `form` component (spread over the number of correction
terms) dominates the weak-coupling budget, and it cannot be argued away. `M(0)`
does not converge in the order: mes N = 4 at `m/g = 0.1` gives 0.409, 0.435,
0.458, 0.474 for `n_terms` = 2, 3, 4, 5, each term still moving it ~3%.

Nor does any criterion select the order. The two sound ones disagree:

* **Held-out K** picks the *largest order offered*, every time — 4 of {2,3,4} on
  2K = 25–49, 5 of {2,3,4,5} on 25–71, 30/30 both. It holds out 2K = 67–71 while
  fitting 25–65, which is near-interpolation, so flexibility always wins.
* **The chiral limit** runs the other way and prefers 2. Mean `|α − 2|` of the
  extrapolated masses over the 0.05–0.1 pair: 0.087, 0.098, 0.110, 0.119 for
  `n_terms` = 2, 3, 4, 5.

Two defensible criteria pointing at opposite ends of the range is what an
undetermined model order looks like. `form` therefore runs over the whole
Eq. (27) ladder, 2 through 5.

### What does not work — measured, do not re-try

| approach | result |
|---|---|
| More K | remainder falls as `K^-a`; halving the bracket at `m/g = 0.1` needs K × 84,000 |
| Confluent basis | spans the same space; `M(0)` identical to `1.1e-14` |
| Committing to `n_terms = 4` | the next admissible term still moves `M(0)` by ~3% |
| Validity floor on the window set | floor 2K ≥ 25 → ≥ 49 shrinks the spread 3.5× but walks the central value monotonically by ~7× the shrunken bar (bar N=3 at `m/g=0.1`: 1.0668 → 1.1045 against 0.0203 → 0.0057) |
| Free-exponent fit `M0 − C·K^-p` | unidentifiable: `p ≈ 0.01` makes `K^-p` nearly constant, `M0` and `C` degenerate; lands 10× off with a `6e-7` residual |
| Increment tail-sum (Hurwitz ζ) | validates at strong coupling, diverges at weak (`q → 1`); 5/30 inside the paper's bars against 26/30 for Eq. (27) |
| Explicit `K^-a` column, `a` fixed | limits 3–4× the published values; `\|α−2\|` = 0.87 |
| `K^-a` column in the joint fit | 300–1400% variant spreads, negative masses |

One thing the budget does **not** double-count: `form` and `window` are
strongly correlated (`window/form` = 0.42–0.48 across all fifteen weak entries),
but measured against the full (order × sub-window) ensemble, quadrature
reproduces the joint spread — `quadrature/joint` = 0.997–1.05 at `m/g ≤ 0.2`.
It *does* over-count at strong coupling (1.29 at `m/g = 1.6`).

## Consequence 3 — the independent fit is biased high

A joint `(m/g, K)` fit anchored at the chiral point disagrees with the
per-coupling Eq. (27) fit, systematically and in the same direction in every
channel. `tools/joint_fit.py` was previously scored over `m/g ≤ 0.45`, where the
`(m/g)²` form is already breaking down (93% variant spread at the top coupling);
restricted to the chiral region it behaves quite differently — `|α − 2|` = 0.113,
0.052, 0.016 at `m/g_max` = 0.45, 0.35, 0.25.

Fitted at `m/g ≤ 0.15` on a grid down to `m/g = 0.025`, the joint value sits
**below** the independent one at all three of the weakest couplings, in all five
channels — 21/25 below, 15/25 outside the independent total bar — and the gap
grows at nearly the same rate everywhere: ~1.1 bars at `m/g = 0.075`, ~1.5 at
0.05, ~1.9–2.1 at 0.025. At `m/g = 0.15` the two agree, which is the check
against the anchor simply pulling everything down.

The chiral referee backs the joint value in every channel:

| channel | joint `\|α−2\|` | independent `\|α−2\|` |
|---|---|---|
| mes N=2 | 0.037 | 0.146 |
| mes N=3 | 0.007 | 0.126 |
| mes N=4 | 0.037 | 0.097 |
| bar N=3 | 0.021 | 0.101 |
| bar N=4 | 0.055 | 0.081 |
| **mean** | **0.031** | **0.110** |

Read this carefully. It does **not** establish that the joint values are right:
by §"the cause" above, both are extrapolating within the pre-asymptotic regime,
and the `M² ∝ (m/g)²` anchor is an anchor to the artifact. What it does
establish is that the `form + window` budget is **missing a systematic** that a
second, independent method exposes at 1–2× the quoted bar.

## A limit form that reaches all `m/g`

`joint_fit`'s limit series is a polynomial in `m/g`. It has the correct
small-`m/g` behaviour built in (the `mg²` prefactor) and *no* large-`m/g`
behaviour at all, which is why it fails above `m/g ≈ 0.25`: `R = M²/(m²+g²/π)`
saturates, and a polynomial cannot.

At large `m/g` the quarks are free, `M → n_q·m` with `n_q = 2` for a meson and
`N` for a baryon, so `R → n_q²`. Verified in the raw data, not assumed — `R/n_q²`
at 2K = 68 runs 1.0603, 1.0320, 1.0145, **1.0060** (mes N=2) and 1.1147,
1.0553, 1.0239, **1.0096** (bar N=3) at `m/g` = 3, 6, 12, 24.

In `t = (m/g)²` a rational function carries both limits exactly:

```
R(t) = t·P(t)/Q(t),   deg Q = deg P + 1,   lead(P)/lead(Q) = n_q²
```

with the ratio **imposed, not fitted**. Max relative residual, best polynomial
against best Padé at comparable parameter count:

| channel | range | polynomial | Padé |
|---|---|---|---|
| mes N=2 | 0.025–1.6 | 9.3e-2 | **6.1e-3** |
| mes N=3 | 0.025–1.6 | 8.4e-2 | **7.4e-3** |
| mes N=4 | 0.025–1.6 | 1.5e-1 | **6.7e-3** |
| bar N=3 | 0.025–1.6 | 1.3e-1 | **5.3e-4** |
| mes N=2 | 0.025–**24** | 6.96 *(696%)* | **8.1e-3** |
| bar N=3 | 0.025–**24** | 6.75 *(675%)* | **1.2e-2** |

Three decades in `m/g` at ~1%, where no polynomial works at all. This is not
yet wired into `joint_fit`, which is a linear SVD solve — the Padé denominator
makes it nonlinear, so it wants a two-stage fit or a rewrite.

## The chiral expansion does not converge either

Fitting `M0/mg² = g₀ + g₁·mg + …` to the continuum values, decoupled from the
K-fit, pins **g₀ to ~4%** per channel (mes N=2 45.9–47.9, N=3 50.6–52.7, N=4
50.4–52.6, bar N=3 117.8–123.4). But the coefficients swing wildly with order
while the residual falls — mes N=2 gives `g₁` = −88.9, −62.9, −22.5 and `g₂` =
−115, −553 and `g₃` = +1315. That is a polynomial chasing a non-polynomial
function, the same trap as `n_terms` in the K-fit, and going to higher order
buys residual rather than knowledge.

A `g₀ + h·ln(mg) + g₁·mg` form fits better at equal parameter count in all four
channels tested (3.8e-3 vs 8.2e-3; 1.3e-3 vs 6.1e-3; 2.6e-3 vs 4.6e-3; 3.0e-3 vs
5.1e-3). **Treat this as fit selection, not a derivation** — with `h > 0` it
sends `M0/mg² → −∞` as `mg → 0`, and 0.025–0.2 is under one decade, too short to
separate a log from a small fractional power. A log *could* arise: the endpoint
integrals carry `1/a` poles (§"Consequence 1"), and logs appear where those
cancel, the same confluent structure that produces `K⁻¹ln K`. Whether the
coefficient is nonzero needs the perturbation theory, which has not been done.

## Defect: the Eq. (27) basis degenerates at *both* ends of `a`

The chiral end is documented in `observables._richardson_design`: `a → 0` sends
`1+a → 1` and `2+a → 2`. The heavy end is not, and it is worse in practice. As
`m/g` grows, `a → 1`, so `1+a → 2` and `2+a → 3` go near-parallel to columns
already present. Condition number runs `1.7e5` at `m/g = 3` to `1.0e7` at 24.
Against the known free limit for bar N = 3 (truth 9) at `m/g = 24`:

```
raw 2K=68   9.086
n_terms=2   9.157      <- the only reliable one
n_terms=3   8.900
n_terms=4  10.752
n_terms=5  -5.438      <- negative M²
```

and mes N = 2 at `m/g = 24` (truth 4): `nt=2` 4.023, `nt=4` 3.660, `nt=5` 6.619.

**Use `n_terms = 2` whenever `a > 0.95`.** `figures._terms_for` keys only off the
point count and would hand back 4. Table I's own range is safe — `a = 0.909` at
`m/g = 1.6`, where `n_terms` 2 through 5 agree to five digits — but any sweep
past `m/g ≈ 3` is not. The confluent basis does not help here: its
`partner_of = {2:1, 4:3}` pairs `1+a` with `1`, which is the `a → 0` merger.

## This was known in 1996, and there is a published fix

Everything above was found independently, but it is not new. Brett van de Sande,
*Convergence of Discretized Light Cone Quantization in the small mass limit*
([hep-ph/9605409](https://arxiv.org/abs/hep-ph/9605409)), diagnoses the same
thing and proposes a cure. Four points from it bear directly on this repo.

**The mechanism, confirmed.** "The source of this slow convergence comes from
the end-point behavior of the wavefunction which is handled poorly by DLCQ" —
and "the DLCQ bound state equation does not make use of our knowledge of the
endpoint behavior."

**The convergence law.** DLCQ eigenvalues converge as powers of `1/K^{2β}` from
the endpoint regions, *plus* powers of `1/K` from the `1/(x−y)²` singularity.
His fit at β = 0.1 is

```
M² = 0.779141 − 1.25943/K^{2β} + 1.15609/K^{4β} − 1.01505/K^{6β} + …
```

against an exact 0.77914. **Eq. (27)'s basis has the Coulomb `K⁻¹` piece but no
`K^{-2β}` term at all** — which is exactly the missing column identified above,
now with its correct exponent (`2β`, not the `β` that the local increment slope
suggests).

**The true chiral law.** His Eq. (7): `M² = 2πgμ/√3 + O(μ²)` — *linear* in the
quark mass, confirming GMOR and confirming that DLCQ's `M² ∝ m²` is an artifact.
He says so outright: the DLCQ results "have an incorrect functional form for M²
versus μ in the small μ limit." Also worth noting, his Eq. (6) gives
`β = (μ√3/gπ)(1 − μ/(10g²) + …)`, the same odd-series structure derived above
from Eq. (26).

**The fix — "improved DLCQ."** Add and subtract a term so the integrand vanishes
when both `x` and `y` are near an endpoint. The discretized form, his Eq. (12),
replaces the off-diagonal sum with

```
g² Ψᵢ I(i/K)  +  g² K Σ_{j≠i} [ Ψᵢ ((j(K−j))/(i(K−i)))^β − Ψⱼ ] / (i−j)²
```

with `I(x) = ∫dy [1 − (y(1−y)/(x(1−x)))^β] / (x−y)²` evaluated exactly. This is
a change to the *Hamiltonian*, not to the fit.

### Do not use the convergence law as a fitting basis on this K range

The obvious shortcut — refit our existing data in `{1, K^{-2β}, K^{-4β}, K^{-6β}}`
— appears to work spectacularly, moving the chiral exponent from 1.85–1.95 to
0.95–1.29 across all five channels, against a truth of 1. **It is an artifact.**
Fed a synthetic series with chiral exponent exactly 2 and a pure `K⁻¹`
K-dependence carrying no endpoint physics at all, the same fit returns **0.33**
(N=2) and **−0.54** (N=3), where Eq. (27) correctly recovers 2.0000.

The reason is degeneracy. At `2β = 0.042`, `K^{-2β}` runs only from 0.898 to
0.861 across the whole window, so `{1, K^{-2β}, K^{-4β}, …}` is a set of
near-identical constants; condition numbers reach `1e7`, and the fitted
intercept drifts with `β` and hence with `m/g`. Where `β` is O(1) the basis is
well conditioned and behaves: at `m/g` = 0.8 and 1.6 it agrees with Eq. (27) to
~0.2% and both match the published values.

The leverage boundary is quantitative. For `K^{-2β}` to fall to 0.5 — the point
at which a fit can see it at all — needs `K`:

| m/g | 0.4 | 0.2 | 0.1 | 0.05 | 0.025 |
|---|---|---|---|---|---|
| K needed | 2.9 | 7.9 | **61** | 3.6e3 | 1.3e7 |

against our `K_paper ≤ 34.5`. Van de Sande demonstrates at β = 0.1 — comparable
to our `m/g ≈ 0.2` — with K up to ~100. So his law is right, and it is precisely
*why* this K range cannot determine `M(0)` at weak coupling; it is not a way to
extract it.

## What would actually work, and what it would cost

Not more K, and not a better fit.

The cost is quantifiable, because the failure has a one-parameter model. The
endpoint integral is `∫₀¹ x^{2a−1}dx = 1/(2a)`, and a mesh reaching down to
`x_min = δ` captures `1 − δ^{2a}` of it. For a uniform mesh `δ = 1/n`, so the
captured fraction is `1 − n^{-2a}` — which approaches 1 only for `n ≫ e^{1/2a}`.

That model accounts for essentially all of the observed non-convergence.
Dividing the `thooft.py` values by their captured fraction at `m/g = 0.05`:

```
n =  200   M² = 0.03585   captured 0.3612   M²/captured = 0.09927
n =  400        0.03912            0.3975                 0.09840
n =  800        0.04234            0.4319                 0.09804
n = 1600        0.04550            0.4642                 0.09802
n = 3200        0.04861            0.4947                 0.09826
```

`M²` moves 36% across that range; `M²/captured` moves **1.3%**. That collapse is
a real test of the *n*-dependence, and it passes.

**Do not use it to "correct" the chiral exponent.** An earlier revision of this
document reported that applying the same rescaling moved the exponent from 1.979
to 1.151, and read that as evidence for the mechanism. It is not. For small `a`,
`1 − n^{-2a} ≈ 2a·ln n`, so the correction factor goes as `1/a ∝ g/m` and
mechanically subtracts about 1 from *any* exponent. Feeding it a synthetic
series that is exactly `C·(m/g)²` with no `n`-dependence at all returns a
"corrected exponent" of 1.17 — indistinguishable from what the real data give.
The rescaling tests convergence in the grid parameter and nothing about the
quark-mass dependence.

The scaling then says what each fix buys. Nodes needed to capture 99% of the
endpoint integral, N = 3 baryon:

| m/g | `a` | `x_min` needed | uniform | geometric (10%/cell) | `x^a` basis |
|---|---|---|---|---|---|
| 0.40 | 0.326 | 8.6e-4 | 1.2e3 | 74 | ~30 |
| 0.20 | 0.168 | 1.1e-6 | 9.2e5 | 144 | ~30 |
| 0.10 | 0.084 | 1.4e-12 | 7.0e11 | 286 | ~30 |
| 0.05 | 0.042 | 2.3e-24 | 4.4e23 | 571 | ~30 |
| 0.025 | 0.021 | 5.3e-48 | 1.9e47 | 1142 | ~30 |
| 0.0125 | 0.011 | 2.9e-95 | 3.4e94 | 2284 | ~30 |

A **graded mesh** turns a cost exponential in `1/a` into one linear in it — a few
thousand nodes spans ninety-five decades — and makes the problem merely hard.
A **basis with `x^a` in it** removes the requirement entirely: after factoring
the singular behaviour the remainder is analytic, so convergence is spectral and
the endpoint is exact by construction rather than resolved.

This is why the established fix is a basis change and not a finer grid. Anand,
Fitzpatrick, Katz and Xin
([arXiv:2111.00021](https://arxiv.org/abs/2111.00021)) restore convergence at
quark masses far below the strong-coupling scale in lightcone conformal
truncation by modifying the basis according to the quark mass, motivated by
't Hooft's endpoint analysis.

**Note what does not transfer to DLCQ.** The graded-mesh option exists only for
`thooft.py`. In DLCQ `x = k/K` is not a numerical choice — it follows from
quantizing momentum in a periodic box, which is exactly what makes the Fock
basis finite at fixed `K⁺`. There is no non-uniform spacing of `k` to be had, so
on that path the basis change is the *only* option.

Nor does the captured-fraction *correction* transfer. DLCQ has the same
`x_min = 1/K_code` structure, so `M²/(1 − K_code^{-2a})` is the obvious thing to
try, and it does flatten the K-dependence (15–19% raw spread down to 0.3–5%).
But it fails its own falsification test: fitting `M² = A(1 − K_code^{-b})` with
`b` free returns `b/2a` anywhere from **0.32 to 4.2**, with no sign of settling
on 1. The flattening is a two-parameter fit absorbing a smooth monotone curve,
not the model being right. Do not use it to correct DLCQ numbers.

## What the DLCQ code can still do

Not reach smaller `m/g` — but two things worth having.

**A validity boundary, from two independent directions.** The endpoint integral
is captured to O(1) only when `2a·ln K_code ≳ 1`. At 2K = 70 that puts the
threshold at `m/g` = 0.105 (mes N=2), 0.140 (N=3, bar N=3), 0.166 (N=4, bar
N=4). Independently, the chirally-anchored joint fit and the per-coupling fit
agree at `m/g = 0.15` (0.02–0.64 bars) and diverge below it (1.0–2.1 bars). Two
unrelated criteria landing on `m/g ≈ 0.15` is the boundary worth quoting:
**Table I's `m/g` = 0.1 and 0.05 entries are grid-limited in every channel**,
and the 0.2 entry is marginal.

**More K, but only above it.** "More K is dead" is a weak-coupling statement and
does not generalize. The remainder goes as `K^-(p−1)` with `p` the measured
increment exponent, so the K factor needed to halve it is:

| m/g | `p_∞` | remainder | K factor to halve |
|---|---|---|---|
| 1.60 | 1.974 | `K^-0.974` | **2.0** |
| 0.80 | 1.893 | `K^-0.893` | **2.2** |
| 0.40 | 1.554 | `K^-0.554` | **3.5** |
| 0.20 | 1.215 | `K^-0.215` | 25 |
| 0.10 | 1.067 | `K^-0.067` | 3.1e4 |
| 0.05 | 1.017 | `K^-0.017` | 5.1e17 |

So at `m/g ≥ 0.4` the sparse solver's reach is a real lever — a factor of 2–3.5
in K halves the bar, and 2K = 71 at the Table I truncation costs 21 s
(`next-steps.md` §1). At `m/g = 0.2` it is marginal. Below that it is hopeless,
and in exactly the region where the grid does not resolve the physics anyway.

## Zero modes: a separate question, and more tractable than it first looks

Zero modes are **not** the cause of anything in §"Consequence 1" — the continuum
control settles that. But they are a real omission with a measured effect, and
the difficulty should not be overstated.

The naive picture — allow `k = 0` in the Fock basis — is wrong and is what makes
it look impossible: at fixed `K⁺` one could then add arbitrarily many
zero-momentum quanta, so the space is infinite-dimensional. That is not how it
is done. Zero modes come in two kinds:

* **constrained** — not independent degrees of freedom; each obeys a constraint
  equation making it a complicated operator-valued function of all other modes.
  This is the genuinely hard part, and it is routinely dropped.
* **dynamical** — a small number of *collective* variables that cannot be gauged
  away. In QCD(1+1) the light-cone gauge `A₋ = 0` is unavailable and the gauge
  zero mode is dynamical and gauge invariant
  ([hep-th/9504026](https://arxiv.org/abs/hep-th/9504026)).

Dropping the constrained modes leaves the dynamical ones governed by an infinite
square-well potential — an ordinary quantum-mechanical problem. Its wavefunctions
tensor onto the existing Fock basis, so the matrix dimension multiplies by the
number of zero-mode levels retained, a modest factor rather than an unbounded
tower.

This has been carried through numerically. Müller, Kalloniatis and Pauli
([hep-th/9803204](https://arxiv.org/abs/hep-th/9803204)) compute the full
bound-state spectrum of two-dimensional SU(2) Yang–Mills coupled to massive
adjoint scalar matter in DLCQ, and find **21% shifts in the masses of the lowest
lying states** once the zero-mode wavefunctions are included. The formalism is
Pauli–Kalloniatis–Pinsky ([hep-th/9403038](https://arxiv.org/abs/hep-th/9403038),
[hep-th/9509020](https://arxiv.org/abs/hep-th/9509020)).

So the honest assessment is: **doable, and worth doing for its own sake**, but
not free here. What survives untouched is the eigensolver, the colour algebra and
the state generator. What is new is (a) a zero-mode sector with its own quantum
mechanics and (b) the Hamiltonian matrix elements coupling it to the non-zero
modes — and those have to be *derived* for this case, because the published
implementations are SU(2) pure glue and adjoint scalars, not SU(N) with
fundamental quarks and baryon number. Gauge fixing also has to be revisited,
since `A₋ = 0` is not available.

Two smaller practical notes. Antiperiodicity for the quarks is baked into
`qcdf_states._filter_parity_group` as a hoisted "all momenta odd" test that
rejects 96.3% of candidates, so it is simultaneously the boundary condition and
the main pruning step — but it can stay, since antiperiodic fermions have no zero
mode and the gluon zero mode is the one at issue. And any such run would no
longer be comparable to `fortran/qcdf.f`, which is antiperiodic and zero-mode-free
throughout, so the bit-identical anchor would not cover it.

## The data

Every number above comes from **`data/chiral_grid_msq.csv`** — all five
channels, `m/g` = 0.025 to 1.6 (and to 24 for mes N=2 and bar N=3), 2K = 25–69
at the Table I truncation, 1446 rows. It supersets `data/joint_fit_msq.csv`,
which is kept because `tools/joint_fit.py` takes it by name.

To regenerate or extend:

```
python tools/large_k_sweep.py --N <N> --B <B> --K-lo 25 --K-hi 69 \
    --mg 0.025 0.05 0.075 0.1 0.15 0.2 --lpn <sweep_lpn(N,B)> --nev 12
```

`--nev 12` is verified to give the identical ground state to a full solve
(0.3012878408 at mes N=2, 2K = 68, `m/g = 0.1`), and it is cheap: bar N=4 at
2K = 68 is ~10 s, the mesons ~1 s.
