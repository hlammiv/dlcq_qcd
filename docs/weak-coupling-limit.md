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

## Implemented: improved DLCQ

`hamiltonian="improved"` on `run_python` / `PythonProvider`, `--hamiltonian
improved` on the CLI. Off by default; the reference path is untouched and
bit-identical.

**It needed no kernel change.** The self-energy enters the Hamiltonian as
`Norm · diag(σ)` — measured exactly, at every parton number — so replacing
`σ_std` by `σ_imp` is a single matrix addition after the existing build.
Nothing in the hot path, the colour algebra or the bit-identity wall moves.

The improved per-state scalar is built from a *directed* pair kernel, one
parton's share of the endpoint sum with a given partner (`P = k + l`):

```
J^(k; l) = Σ_{j odd, 0<j<k} (j(P−j)/(k·l))^b · 4/(k−j)²  +  I_b(k/P)/P
σ_imp(s) = Σ_a  C_F/(L−1) · Σ_{c≠a} J^(k_a; k_c)
```

At `b = 0` the weight is 1 and, since `k` and `j` are both odd, `k−j` runs over
the even integers `2..k−1`, so `J^(k;l) → S(k)` **independently of the
partner** — and `σ_imp → σ_std` exactly, per parton rather than merely per
state. That is what makes it a safe drop-in, and it is the gate the tests run
first.

The `C_F/(L−1)` weight is not arbitrary in valence sectors: for a colour singlet
`Σ_a T_a = 0` gives `Σ_{c≠a}(−T_a·T_c) = C_F` per parton, and when all pairs sit
in the same channel that forces `C_F/(L−1)` each — matching the measured
`c_qq = (N+1)/2N` for N = 2…6. Where pairs are *not* colour-equivalent
(Fock-extended sectors carrying more than one colour singlet, i.e. Table I's
`LPN = valence+2`) it is a scalar stand-in for the matrix `−T_a·T_c`; measured
effect ~3e-3 in M/g. See `docs/next-steps.md`.

**Measured, N = 3 meson at m/g = 0.1**, over 2K = 10–70, with bars from an
ensemble of fit form × contiguous sub-window (each family using the exponents
appropriate to it — `K^{-2β}` powers for standard per his Eq. (13), plain `1/K`
for improved per Eq. (14)):

| | M(0) | bar | at 2K=70 | extrapolation supplies |
|---|---|---|---|---|
| standard | 0.7792 | ±0.0679 (8.7%) | 0.329 | **58%** |
| improved | 0.8296 | **±0.0003 (0.04%)** | 0.827 | **0.3%** |

They **agree** (0.7σ), and the improved bar is ~200× tighter. That is the whole
point: the gain is not a moved central value, it is that almost nothing is left
to extrapolate. Standard drifts 48% across that K range; improved drifts 2.4%.

**And it fixes the chiral exponent.** `M² ∝ (m/g)^α` with α = 1 physically
(GMOR; his Eq. 7). Standard returns 1.954, 1.951, 1.946 at 2K = 30, 40, 60 —
the artifact this document is about. Improved returns **1.0696, 1.0688,
1.0678**: near 1 and *stable in K*, so it is a property of the computed masses
rather than of a fit. This is the one non-circular demonstration that the `m²`
artifact becomes the physical `m¹` — the operator changed, not the fit basis
(contrast the trap recorded above).

**It runs at every parton number**, valence baryons included. The N=3 valence
baryon (`LPN=3`) at m/g = 0.1 over 2K = 15–45: standard drifts **25.8%**,
improved **4.5%**; a plain 1/K fit has residual 1.1e-5 against standard's
3.5e-4; and the chiral exponent goes **1.958 → 1.045** at 2K=21 and
**1.952 → 1.037** at 2K=33 — closer to the physical 1 than the meson manages.

An earlier revision of this section claimed the underlying identity "fails by
8.3e-01 at three partons", and restricted the flag on that basis. **That was
wrong**, and the error was mine: the test behind it compared the self-energy
against an *unweighted* row sum of the exchange, which is only valid where the
norm matrix is the identity — true in the meson valence sector and nowhere
else. The correct, norm-weighted statement holds at every L:

```
D ≡ H(selfen) − H(selfen=0)  =  diag(σ_std) · Norm
σ_std(s) = C_F · Σ_{partons p} S(k_p),  C_F = (N²−1)/2N,  S(k) = Σ_{n≤(k−1)/2} 1/n²
```

Verified here to 1e-16…5e-18 relative at L = 2, 3, 4, 5 and in mixed sectors up
to L = 8, including the untruncated basis. Two consequences: the self-energy
stays a **one-body scalar** at every L, and `Norm_ij ≠ 0` only between states of
identical parton content — so any replacement enters as `Norm · diag(σ_new)` and
needs only a per-state scalar. The `8.3e-01` was the `(5,5,5)` state's `1 − 1/6`
with `6 = Norm_ii/L!`: a normalisation number, not physics.

### Table I under the improved Hamiltonian

Run over Table I's own grid — five channels, six couplings, 2K = 25–71, at the
published truncation `LPN = sweep_lpn = valence+2`. Data in
`data/improved_table1_msq.csv`; the standard-Hamiltonian companion is
`data/chiral_grid_msq.csv`. Improved values are a plain 1/K extrapolation with
the bar from an ensemble over orders × contiguous sub-windows; standard values
are `richardson_budget` at four Eq. (27) terms.

| channel | m/g | improved | standard | paper | imp vs std |
|---|---|---|---|---|---|
| mes N=2 | 1.60 | 4.3139 ± 0.0001 | 4.3139 ± 0.0001 | 4.314 | 0.0σ |
| mes N=2 | 0.10 | 0.5725 ± 0.0037 | 0.3938 ± 0.0198 | 0.380 | 8.9σ |
| mes N=2 | 0.05 | 0.2595 ± 0.0023 | 0.1089 ± 0.0070 | 0.100 | 20.4σ |
| mes N=4 | 0.05 | 0.4502 ± 0.0013 | 0.1224 ± 0.0091 | 0.120 | 35.8σ |
| bar N=3 | 0.05 | 0.9694 ± 0.0045 | 0.2858 ± 0.0242 | 0.310 | 27.7σ |
| bar N=4 | 0.05 | 2.2482 ± 0.0063 | 0.5322 ± 0.0457 | 0.420 | 37.2σ |

**They agree exactly where standard is trustworthy and diverge exactly where it
is not.** At m/g ≥ 0.8 every channel agrees to 0.0–0.6σ. Below that the gap
opens as it must if the chiral exponents differ (1.93 standard, 1.03 improved):
1.45× at m/g = 0.1 rising to 2.38× at 0.05.

Improved also barely moves across the window — drift 0.0–2.6% against standard's
14–21% — so ≤3.5% of its answer is extrapolated, against 15–35% for standard.

#### GMOR in the chiral limit, N = 2..5: improved converges, standard does not

The GMOR check was first run at a single coupling, where a ratio away from 1
mixes two effects — the `O(μ²)` correction to GMOR and finite `N`. They separate
by pushing `m/g` down, because GMOR is exact only as `m/g → 0`. Scanning
`m/g = 0.1, 0.05, 0.025, 0.0125` at `N = 2,3,4,5` (`data/gmor_scan/`,
2K = 26–70, LPN = 4), `M²(K→∞) / GMOR`:

| m/g | N=2 | N=3 | N=4 | N=5 |
|---|---|---|---|---|
| 0.1 | 1.0595 | 1.0893 | 1.0875 | 1.0834 |
| 0.05 | 0.9386 | 1.0219 | 1.0311 | 1.0333 |
| 0.025 | 0.8840 | 0.9992 | 1.0095 | 1.0127 |
| **0.0125** | 0.8610 | **0.9946** | **1.0023** | **1.0048** |

**Read with the extrapolation bar, this is a two-parameter law, not a ratio
sitting on 1** — and it is worth being careful here, because "converges to 1
within 0.5%" is what this table looks like and is not what it says.

At fixed `m/g` the ratio does **not** sit at 1. With the ensemble bar (fit
orders × contiguous sub-windows, 68% half-width) at `m/g = 0.0125`:
N=3 → 0.9942 ± 0.0027, N=4 → 1.0020 ± 0.0017, N=5 → **1.0046 ± 0.0012**, i.e.
N=5 is **3.8σ above 1**, and the N-trend is still rising. Extrapolating in
`1/N²` over N = 3,4,5 gives `r_∞ = 1.0105 ± 0.0024` — 4.4σ high. Taken alone
that reads as a failure to converge.

It is not, because the residual is proportional to `m/g`:

| m/g | large-N ratio − 1 (N = 3…8) | dev/(m/g) |
|---|---|---|
| 0.1 | 0.0712 ± 0.0014 | 0.712 |
| 0.05 | 0.0347 ± 0.0013 | 0.694 |
| 0.025 | 0.0172 ± 0.0010 | 0.688 |
| 0.0125 | 0.0086 ± 0.0007 | 0.688 |

`dev/(m/g)` is flat to 3% across a factor of 8 in coupling — the linearity is
visible before any fit is run.

Fitting `dev = B + A·(m/g)` over the full N = 3…8 range gives
**B = −0.00048 ± 0.00075, 0.6σ from zero**, with χ²/dof = 0.15 and slope 0.714.
The quadratic alternative `B + A·(m/g)²` is excluded — intercept **18.4σ** from
zero, χ²/dof = 33.5. (With only N = 3,4,5 in the large-N fit these read
0.0005 ± 0.0025 at 0.2σ and 6.9σ / 4.35; extending to SU(8) tightened the
intercept ~3× and sharpened the linear-vs-quadratic call from 6.9σ to 18.4σ.) So

```
M² / GMOR  =  1 + 0.71 (m/g)      (large N, N = 3..8)
```

and the double limit `N → ∞` then `m/g → 0` reaches 1 to within 0.25%. A
*relative* `O(μ)` correction is exactly what should be there, since Eq. (7) is
only the leading term of the mass expansion. The chiral exponent agrees
independently: 1.005, 1.008, 1.009 at N = 3,4,5 on the last coupling pair.

Standard, on the same data:

| m/g | N=2 | N=3 | N=4 | N=5 |
|---|---|---|---|---|
| 0.1 | 0.7303 | 0.6175 | 0.5370 | 0.4819 |
| **0.0125** | 0.1021 | 0.0847 | 0.0715 | **0.0632** |

The ratio **halves every time `m/g` halves** — the signature of `M² ∝ m²`
against GMOR's `M² ∝ m` — with exponent 1.986–1.993. It is not converging to
GMOR; it is diverging from it at a fixed rate.

**N = 2 is the one exception, and it separates cleanly into exponent and
coefficient.** Its ratio moves *away* from 1 (1.0605 → 0.8619) and does not
converge. But its chiral *exponent* still goes to 1 (1.036 on the last pair) —
so what fails at N = 2 is Eq. (7)'s coefficient `2π/√3`, which is a large-N
result, not the GMOR scaling itself. N = 2 is also special on its own terms, its
"meson" and "baryon" being the same state by pseudo-reality — which this code
reproduces to four digits and which Hamer observed independently. The −14% is a
finite-N deviation of the formula, not evidence against the Hamiltonian.

This is the third independent meson anchor, and the sharpest, because what is
confirmed is a *law with an intercept* — `1 + 0.71(m/g)`, intercept consistent
with zero at 0.6σ, linear form preferred over quadratic at 18.4σ — rather than a
single ratio eyeballed against 1.

#### External anchors: what exists, and what each one can settle

Four independent checks now bear on the improved Hamiltonian. Two are decisive,
two are consistent but not discriminating, and one more exists that cannot be
applied here at all. Listing them together because the weak ones were twice
mistaken for strong ones during this work.

| anchor | what it is | verdict |
|---|---|---|
| van de Sande exact, two-body | polynomial-basis solution of the 't Hooft equation, `M²/g² = 0.779141` at β = 0.1 | improved to **5.5e-4** — decisive |
| GMOR (his Eq. 7) | `M² = 2πgμ/√3`, parameter-free, large-N **meson** | improved **1.03**, standard **0.20** — decisive |
| Hamer 1982, meson | SU(2) finite-lattice, **with quoted errors** | both ≤0.7σ; direction favours improved |
| Hamer 1982, baryon | same table, SU(2) baryon | both ≤0.7σ; direction favours improved |
| Date–Frishman–Sonnenschein 1987 | bosonized mass formula, `α = 2N_c/(2N_c−1)` | **failed at finite N_c** (1.01 measured vs 1.20) — see below |

**Hamer's actual error bars matter, and are much wider than the figure suggests.**
Table 1 of Nucl. Phys. B195, 503 gives both columns on our own coupling grid:

| m/g | M_B/g | M_meson/g |
|---|---|---|
| 0.05 | 0.25 ± 0.2 | 0.30 ± 0.2 |
| 0.10 | 0.40 ± 0.2 | 0.45 ± 0.2 |
| 0.20 | 0.65 ± 0.15 | 0.65 ± 0.15 |
| 0.40 | 1.05 ± 0.10 | 1.08 ± 0.07 |
| 0.80 | 1.90 ± 0.05 | 1.92 ± 0.03 |
| 1.60 | 3.50 ± 0.05 | 3.51 ± 0.02 |

±0.2 on 0.30 is ±67%. An earlier revision of this work read points off Fig. 8(a)
instead and assigned them a ~1.8% *digitisation* bar, which manufactured an
apparent 5σ tension at m/g = 0.2. With Hamer's own errors everything sits inside
0.7σ and the tension does not exist. **Use Table 1, not the figure.**

Worth noting as an independent structural check: our SU(2) baryon and SU(2)
meson come out identical to four digits at every coupling, reproducing Hamer's
own observation that "the ratio of the baryon mass to the meson mass is
consistent with 1 everywhere". For N = 2 both are two-quark states.

**The DFS baryon formula: three transcription errors, now corrected.** Date,
Frishman & Sonnenschein (Nucl. Phys. B283, 365) derive, by non-abelian
bosonization at strong coupling, their Eq. (4.1). An earlier revision of this
document transcribed it as

```
WRONG:  E = 4m√(2N_c/π) + √(½)·m·√( (π/N_c)³ [ C₂ − N_c²·N/(2(N+1)) ] )
```

and concluded from it that DFS predicts baryon exponent **2** against GMOR's
**1** for mesons. That conclusion is **withdrawn**. Three things were wrong:

1. **The Casimir bracket sits OUTSIDE the radical**, not inside:
   `E = 4m√(2N_c/π) + √(½)·m·√((π/N_c)³)·[C₂ − N_c²N/(2(N+1))]`.
   This is settled by DFS themselves, not a judgement call: they quote a
   decuplet/octet ratio of **1.41** on the same page, and bracket-outside gives
   **1.4112** while bracket-inside gives 1.2374.
2. **`N` in DFS is `N_f − 1`**, not the flavour count. Their flavour index runs
   `i = 1..N+1`; and they describe `N = 2, N_c = 3` as "decuplet to octet",
   which is SU(3) *flavour*, i.e. three flavours.
3. **`m` is not the quark mass.** It is a dynamically generated scale defined
   under their Eq. (2.11), `m ∝ m_q^{1/(1+P)}` with
   `P(N_c,N_f) = (N_c²−1)/(N_c(N_c+N_f))`. So "linear in `m`" is *not*
   "linear in `m_q`". A one-line check settles it: were `m` the quark mass, the
   strong-coupling baryon mass `4m√(2N_c/π)` would carry no coupling at all.

**Consequences of the correction.** NF = 1 is perfectly well defined — the
bracket is exactly `0`, giving `E = 4m√(2N_c/π)`, not an imaginary root. And
since after bosonization the effective action carries exactly *one* dimensionful
parameter `m`, every mass is a pure number times `m`: **DFS predicts no
meson/baryon split at all.** The prediction is a single channel-universal
exponent

```
α_DFS = 2/(1 + P(N_c,N_f))   =   2N_c/(2N_c − 1)   at N_f = 1
```

— 1.333 at N_c = 2, **1.200** at N_c = 3, 1.111 at N_c = 5, and **→ 1 as
N_c → ∞**, i.e. it *reduces to GMOR* in the large-N limit where GMOR is derived.
Independent support: Frishman & Sonnenschein, Phys. Rept. 223 (1993) 309 quote
`M ∝ m_q^{N_c/(2N_c−1)}`, and Steinhardt's `ν = 2/3` at N_c = 2 is the same
`α = 4/3`.

**The corrected anchor is then failed in the meson channel, where improved was
supposed to be validated.** Measured on the improved data, in physical units
(`M²/g² = M²_code·((m/g)² + 1/π)`), against the corrected law:

| N_c | α (0.0125→0.025) | vs GMOR's 1 | α_DFS = 2N_c/(2N_c−1) | vs DFS |
|---|---|---|---|---|
| 2 | 1.0381 ± 0.0169 | +2.3σ | 1.3333 | −17.5σ |
| 3 | 1.0064 ± 0.0065 | +1.0σ | 1.2000 | **−30.0σ** |
| 4 | 1.0102 ± 0.0042 | +2.4σ | 1.1429 | −31.6σ |
| 5 | 1.0112 ± 0.0031 | +3.6σ | 1.1111 | −32.3σ |
| 6 | 1.0114 ± 0.0024 | +4.7σ | 1.0909 | −32.8σ |
| 7 | 1.0114 ± 0.0020 | +5.8σ | 1.0769 | −33.3σ |
| 8 | 1.0112 ± 0.0017 | +6.8σ | 1.0667 | **−33.6σ** |

Bars are the ensemble 68% half-width propagated through the log-ratio, treating
the two couplings as independent; they are positively correlated (same K window,
same fit family), so this *overstates* them.

Two readings, both worth stating. **DFS's finite-N_c law is excluded at ~30σ.**
And α also sits above 1 by up to 6.8σ — but that residual is already accounted
for: the separately measured transient `M²/GMOR = 1 + 0.71(m/g)` implies a
*local* exponent of **1.013** on this coupling pair, against 1.011 measured. So
the excess is the known O(m/g) contamination of a finite-coupling difference,
not a third exponent, and α → 1 underneath it.

**DFS requires `α` to fall by 0.133 from N_c = 3 to N_c = 8. The measured change
is +0.005** — flat. The gap closes only because DFS itself tends to 1; improved
agrees with it exactly where DFS coincides with GMOR, and nowhere else.

**So the state of the baryon question is different from what was recorded, and
worse.** There is no meson/baryon split to fail, so "improved fails DFS for
baryons" is withdrawn along with the prediction that generated it. What replaces
it is a sharper and more uncomfortable question: **two analytic results disagree
at finite N_c** — GMOR says `α = 1`, bosonization says `α = 2N_c/(2N_c−1)` — they
agree only as `N_c → ∞`, and improved DLCQ lands on 1 in every channel. Nothing
in this work adjudicates which is right at N_c = 3.

One structural argument, offered and not refuted, cuts against DLCQ rather than
against the subtraction: the anomalous exponent `1/(2N_c−1)` arises from the
condensate/WZW sector, and this code uses antiperiodic boundary conditions with
no fermion zero modes and no vacuum structure. If that is where the correction
lives, **no** DLCQ variant here — standard or improved — can produce it, and the
flat 1.01 is a property of the framework, not evidence about the Hamiltonian.

**And the zero-mode dismissal recorded earlier in this document does not cover
this.** "Consequence 1" argues *this is not the missing zero-mode sector*, on the
grounds that `dlcq/thooft.py` reproduces α = 2 with no truncation. That argument
is sound for what it addressed — the **standard** Hamiltonian's α = 2 grid
artefact — but `thooft.py` solves the **large-N** 't Hooft equation, which cannot
exhibit a finite-N_c anomalous exponent by construction. It is therefore no
evidence at all on the present question, and must not be reused for it.

#### The blindness is real, but it is not about zero modes

The hypothesis floated above — that antiperiodic DLCQ cannot see `1/(2N_c−1)`
because it has no zero modes and no condensate — is **wrong on its stated
grounds**, and the true reason is narrower, provable, and about this code rather
than about light-front quantisation.

Wrong on two counts. Anand, Fitzpatrick, Katz & Xin (arXiv:2111.00021) quantise
on a null plane at N_c = 3 with a trivial Fock vacuum and no zero modes and
*do* resolve the sine-Gordon IR; and the condensate attribution was backwards —
the anomalous exponent corresponds to a condensate that **vanishes** as
`m_q^{1/(2N_c−1)}`, whereas `α = 1` is the statement that needs a *constant*
condensate.

The real reason, verified here to `1e-13`: **in the valence meson sector `N_c`
enters the Hamiltonian only through `C_F`.** Writing `H = m²A + g²C_F B` with
`A`, `B` independent of `N_c` and `N_f` on an `N_c`-free basis,

```
M² / (g² C_F)  =  F(u),        u = (m/g)² / C_F
```

is a single universal curve. Measured at 2K = 20, LPN = 2, for
N_c = 2, 3, 4, 5, 8 and both Hamiltonians, at u = 0.05, 0.01, 0.002, the spread
across all five colours is **2e-14 to 4e-13**. Eq. (26)'s `a` does depend on
`N_c`, but only through the same `u` (`units.endpoint_exponent` has
`rhs = π(m/g)²/C_F`), so it is a reparametrisation, not independent information.

A single universal curve admits exactly **one** chiral-limit exponent. So the
valence meson sector cannot return an `N_c`-dependent α whatever the physics is,
and **our flat α = 1.01 is forced by the construction rather than measured.**

**Baryons do not collapse, and that is the important corollary.** The same test
at valence baryon (B=1, LPN=N), improved, N = 3,4,5:

| u | mesons | baryons |
|---|---|---|
| 0.05 | 0.63619302 ×3, spread **1.9e-14** | 1.561 / 2.832 / 4.416, spread **1.8** |
| 0.01 | 0.23817308 ×3, spread **5.2e-14** | 0.591 / 1.072 / 1.668, spread **1.8** |

Nor is the baryon spread trivial quark-counting: `M_B ∝ N` would give successive
`M²/(g²C_F)` ratios of 1.26 and 1.22, against 1.81 and 1.56 measured. So the
**baryon sector is the only part of this calculation carrying genuine finite-`N_c`
information**, and the meson α is silent by construction.

This cuts both ways. The improved baryons are the most valuable numbers here —
the only ones that can see the physics in dispute — and the least validated:
no external anchor exists for N ≥ 3 (Hamer's SU(2) "baryon" is the same state as
his meson by pseudo-reality), and the corrections are largest there (bar N=4 at
`m/g = 0.05` moves 0.53 → 2.25). The first genuine constraint on them is the
strong-coupling limit of the meson/baryon ratio, which **must** approach `2/N` at
`m/g = 1.6` where masses are quark-mass dominated. It
is not evidence against DFS. What breaks the collapse is only the higher Fock
content (LPN = 4 shifts M² by 1–2.5%), and that moves α by **+0.011** across
N_c = 3→8 against DFS's required **−0.133** — wrong sign, and an order of
magnitude too small to discriminate.

This also retires an argument made earlier in this session: that "the code sees
finite-N structure in the coefficient, so any blindness would have to be specific
to the anomalous dimension". The coefficient's `N_c` dependence is real but comes
from the sub-leading Fock sectors, and it is far too weak to carry the exponent.

#### The ratio and the exponent are the same parameter

`units.meson_baryon_ratio_bosonization` returns `2 sin[π/(2(2N−1))]`. With
`ν = 1/(2N_c−1)` that is `2 sin(πν/2)`, and the bosonization exponent is
`α = 1 + ν`. Verified identical at N_c = 2,3,4,5. **These are two readings of one
number, not two independent tests of bosonization** — agreement of the mass ratio
to ~10% (which the paper reports) and the chiral exponent stand or fall together.

#### What the finite-N_c literature actually says

Two independent determinations exist at N_c = 3, and both agree with α → 1:

* **AFKX (arXiv:2111.00021)**, lightcone conformal truncation: local exponents
  1.235 → 1.050 → **1.002** as `m_q` falls 0.5 → 0.1 → 0.02 → 0.001, descending
  monotonically toward 1. At `m_q = 0.02`, α = 1.2 predicts M² = 0.11487 against
  0.06333 measured — a factor 1.81.
* **Kochergin (arXiv:2405.04031)**, 1/N-corrected 't Hooft perturbation theory,
  states in print that bosonization does not give the correct scaling — and that
  he does not know why.

So the numerical consensus is against DFS's finite-N_c exponent, the mechanism of
its failure is unexplained by anyone including the authors who found it, and
**every** determination in existence — AFKX, Kochergin, this repo — is quantised
on a null plane with a trivial vacuum. That single shared assumption is the one
the anomalous exponent is about, which is why "unresolved" remains the honest
verdict rather than "DFS is wrong".

#### The meson sector is blind to `nu` in *every* observable, not just the exponent

The `C_F` collapse is a statement about the Hamiltonian (`H = m²A + g²C_F B`), so
it forces the **whole spectrum** to collapse, not just the ground state. Then no
meson observable — scaling or ratio — can carry `N_c`. Tested on the
excited-state ratio at matched `u = (m/g)²/C_F`:

| u | N=3 | N=5 | N=8 | sine-Gordon `2cos(πν/2)` |
|---|---|---|---|---|
| 0.02 | 2.624551 | 2.624551 | 2.624551 | 1.902 / 1.970 / 1.989 |
| 0.005 | 3.619510 | 3.619510 | 3.619510 | 1.902 / 1.970 / 1.989 |

Identical to seven digits across three colours, where sine-Gordon needs a 4.6%
spread. And the DLCQ ratio **diverges** as `u → 0` (2.62 → 3.62) where the
sine-Gordon breather ratio is finite: the valence meson is the 't Hooft tower,
whose pion goes massless while the excited states stay finite, whereas
sine-Gordon breathers all scale together. That is exactly right, because the
valence meson sector *is* the 't Hooft equation with `g² → g²C_F`.

**Why this matters for reading AFKX.** Their sine-Gordon evidence at N_c = 3 is
`m₂/m₁ = 1.90421` against `2cos(πν/2) = 1.90211`. That cannot come from a
valence-truncated spectrum — LCT carries the full Fock space. It matches this
repo's own finding that the ratio needs the higher sectors: the first qq̄ pair
moves the meson/baryon ratio 6.2%, and only then does it land on
`2 sin(πν/2)`.

So the unified picture is: **valence = 't Hooft = no `nu`; the bosonization
structure lives in the higher Fock sectors**, and the baryon channel is the only
place this code can see it, because the baryon is the only sector that does not
collapse.

#### The meson/baryon ratio recovers the bosonization parameter

`M_meson/M_baryon` needs the baryon, which (unlike the meson) does **not**
collapse in `C_F` — so unlike the chiral exponent it can carry `N_c`. Improved
Hamiltonian, `tools/meson_baryon_ratio.py`:

| N | m/g=1.6 | 0.4 | 0.05 | 0.0125 | `2sin(πν/2)` | `2/N` |
|---|---|---|---|---|---|---|
| 2 | 1.0000 | 1.0006 | — | — | 1.0000 | 1.0000 |
| 3 | 0.6562 | 0.6364 | — | — | 0.6180 | 0.6667 |
| 4 | 0.4869 | 0.4640 | — | — | 0.4450 | 0.5000 |
| 5 | 0.3865 | 0.3643 | 0.3476 | **0.3477** | 0.3473 | 0.4000 |
| 6 | 0.3192 | 0.2993 | 0.2848 | **0.2844** | 0.2846 | 0.3333 |
| 7 | 0.2764 | 0.2538 | 0.2412 | **0.2406** | 0.2411 | 0.2857 |
| 8 | — | 0.2204 | 0.2091 | **0.2089** | 0.2091 | 0.2500 |

A monotone descent from constituent counting at strong coupling — where the
quark masses dominate and nothing else is possible — onto the bosonization value
in the chiral limit. `N = 2` returns `1.0000 ± 0.0006`, the pseudo-reality
identity, as a free pipeline check. (The `m/g = 1.6` entries at N = 7, 8 carry
bars of ±0.08 and ±0.098 and are not usable; the strong-coupling baryon
extrapolation degrades at large N.)

**The standard Hamiltonian sits on `2/N` at every coupling** (0.7–1.0σ) and never
departs. So this is a direct consequence of the van de Sande endpoint
subtraction, in the sector where its corrections are largest.

**Fock convergence**, the check that could have killed it — deviation from
`2sin(πν/2)` at `m/g = 0.025`:

| N | valence | valence+2 | valence+4 |
|---|---|---|---|
| 5 | −6.31% | +0.02% | **−0.41%** |
| 6 | −6.26% | −0.02% | **−0.16%** |
| 7 | −6.19% | −0.03% | **+0.04%** |

The first qq̄ pair moves the ratio ~6.2% at every N; the second moves it 0.43%,
0.14%, 0.07% — converging, and converging *faster* as N grows. The sub-0.1%
agreements at valence+2 (which the scan uses) are fortuitously precise; the
converged statement is **within 0.4%, improving with N**, against 15–19% for
counting.

**This is the first external-facing check the improved baryons have passed.**
They have no anchor at N ≥ 3 — Hamer's SU(2) "baryon" is his meson by
pseudo-reality — and their corrections are the largest anywhere (bar N=4 at
`m/g = 0.05` moves 0.53 → 2.25). The `m/g = 1.6 → 2/N` limit is forced by
physics alone, and they reproduce it.

#### Why our two readings of the same parameter disagreed by 10x

Bosonization is one parameter, `nu = 1/(2N-1)`, and this repo measures it two
ways: the chiral exponent `alpha = 1 + nu`, and the meson/baryon mass ratio
`r = 2 sin(pi nu / 2)`. They disagreed violently — and the resolution is not
that either is wrong.

| N | nu from ratio | nu from exponent | `1/(2N-1)` |
|---|---|---|---|
| 5 | 0.1112 | 0.0112 | 0.1111 |
| 6 | 0.0908 | 0.0114 | 0.0909 |
| 7 | 0.0768 | 0.0114 | 0.0769 |
| 8 | 0.0666 | 0.0112 | 0.0667 |

The ratio recovers `1/(2N-1)` to four decimals. The exponent column does not
move with N at all, which is the giveaway: any `1/N`-corrected law must vary.

**The exponent's "nu" converges to zero, not to 0.1111.** At N=5, as `m/g` falls
it reads 0.0683, 0.0290, 0.0113 — roughly `∝ m/g`. That is the separately
measured O(m/g) transient (`M²/GMOR = 1 + 0.71(m/g)`), not `nu`. Meanwhile
ratio-nu reads 0.1118, 0.1112, 0.1111 and is converged by `m/g = 0.05`.

**So the two are not both measurements of `nu`.** `r = 2 sin(pi nu / 2)` is a
statement about the *shape* of the spectrum, fixed by `β²` throughout the scaling
regime. `alpha = 1 + nu` is a statement about the *overall scale*, and only holds
below the crossover at `m_q ~ 1.7e-4` — 70x below our smallest coupling. Above
it, `alpha -> 1` is the correct 't Hooft answer, which is what we measure, more
emphatically the smaller `m/g` gets.

The `C_F` collapse is a *second*, independent reason the meson exponent cannot
show `N_c` dependence — but it is not the main one. Even an unblinded exponent
would read ~1 in this window.

**AFKX see the same split, at the same data point.** Their `m₂/m₁ = 1.90421`
matches the sine-Gordon `sin(π/5)/sin(π/10) = 1.90211` to 0.1%, while their mass
gap sits on the 't Hooft line to 0.8%. Spectrum ratios sine-Gordon, absolute
scaling 't Hooft — two independent methods, one pattern.

**Consequence for how this document reads its own results.** The "~30σ exclusion
of DFS" recorded earlier came from the exponent, i.e. from a quantity that in
this window measures the transient rather than `nu`. It is withdrawn as evidence
about bosonization. What survives is that in the accessible window ratios have
converged and scaling has not, and no measurement by any method has yet gone
below `m_q ~ 10⁻³`.

#### The equal-time measurement exists, and it favours bosonization

Bañuls, Cichy, Cirac, Jansen & Kühn, PRX **7**, 041046 (2017),
[arXiv:1707.06434](https://arxiv.org/abs/1707.06434) — (1+1)d SU(2), N_f = 1
staggered, gauge field eliminated exactly by Gauss' law, MPS/DMRG. Their Table 1
fits the vector mass gap against the **quark mass** over `m/g ∈ [0.1, 0.4]`:

| j_max | ν in `M ~ (m/g)^ν` |
|---|---|
| 1 | 0.781(93)(65) |
| 3/2 | **0.700(29)(11)** |
| 2 | 0.700(29)(12) |

So `α = 2ν = 1.40`: **+1.07σ from bosonization's 2/3, +6.45σ from α = 1**. They
cite Steinhardt for the 2/3. (Verified against the paper's full text — an
earlier reading of this work as *lattice-spacing* scaling rather than chiral
scaling was wrong.)

#### Three corrections to arguments made earlier in this document

1. **There is no U(1)_A anomaly for SU(N) in 2d.** The abelian axial anomaly is
   `∝ ε^{μν} tr F_{μν}`, which vanishes for SU(N). So the N_f = 1 meson *is*
   exactly massless at `m_q = 0` and the chiral-exponent framing is sound.
   Measured: `c_IR = 1.04` against the free-boson `c = 1` (Fujikura & Hidaka,
   arXiv:2605.17183, equal-time VUMPS), the flow to WZW confirmed on the lattice
   with overlap fermions for `(N_c,N_f) = (2,1), (3,1)`
   (Karthik–Narayanan–Narayanan, arXiv:2312.13929), and a `c = 1` plateau in
   AFKX's own c-function.
2. **The Schwinger model is not the N = 1 member of this family.** `2N/(2N−1)`
   at N=1 gives 2, and the Schwinger boson is *massive* at `m_q = 0`
   (`M = e/√π`) — but that is the *gauged abelian* anomaly generating
   `−(e²/2π)φ²`, a term absent when `tr F = 0`. So this is not a falsification.
   **The correct structural check is the two-flavour Schwinger model**, where
   the identical logic `M ~ m^{1/(2−Δ)}` gives `m^{2/3}` — and *that* is
   confirmed by equal-time DMRG (Dempsey–Klebanov–Pufu–Zan arXiv:2206.05308 and
   several tensor-network follow-ups). The framework passes its one clean
   equal-time test.
3. **Coleman's theorem makes `α = 1` untenable asymptotically.** GMOR
   (`M²f_π² ∝ m_q⟨ψ̄ψ⟩`) with `α = 1` *requires* `⟨ψ̄ψ⟩ → const ≠ 0`, which
   Coleman forbids at finite N (large N evades it; Zhitnitsky's
   `⟨ψ̄ψ⟩ = −N_c/√12` is a large-N statement). So the exponent **must** leave 1
   at small enough mass.

#### The reconciliation: a crossover below everything computed

The light-front numerics are not wrong either. AFKX's `M² = 0.0031251` at
`m_q = 2^-10` sits on the 't Hooft line `M² = (2π/√3)m_q·(1−1/N²) = 0.003149`
to **0.8%**. Both facts survive if the anomalous exponent takes over below the
computed window. With Kochergin's normalisation at N=3, `M²_bos = 18.2 m^1.2`
crosses `M²_tHooft = 3.22 m` at

```
m_q ~ 1.7e-4        -- a factor ~6 below the smallest mass anyone has computed
```

Above that crossing bosonization *overshoots* (33% in M² at `m_q = 10⁻³`), which
is what the measured local slope 1.002 excludes. **Nobody has looked below
`m_q ~ 10⁻³` at finite N by any method.** That is the experiment that decides it.

The cleanest discriminator remains unmeasured at finite N: `⟨ψ̄ψ⟩(m_q)` must be
flat if `α = 1` and go as `m_q^{1/(2N−1)}` — at N=2, `m^{1/3}` — if the RG is
right. The DMRG machinery exists (Fujikura–Hidaka; Hayata–Hidaka–Nishimura,
arXiv:2311.11643, SU(2) *and* SU(3), N_f=1) and simply never went below
`m/g = 0.5`.

**A caution on the measured baryon exponent itself.** The improved eigenvalue in
the chiral region is, at the K this repo reaches, close to a linear map of the
Eq. (26) exponent `b` that is *fed in as input* and satisfies `b ∝ m/g`. To that
extent `α ≈ 1` is that input propagated, not an independent dynamical
measurement — which is a further reason not to read the baryon 1.04 as
confirmation of anything.

#### The external check, which is what makes this assertable

Everything above is internal consistency or the two-body anchor. Van de Sande's
Eq. (7) is neither: `M² = 2πgμ/√3` is the GMOR law of the continuum theory, with
no free parameters. Mapping his units onto the repo (`M²/g² = (M/g)²/c`,
`μ²/g² = (m/g)²/c`, `c = (N²−1)/2Nπ`) predicts the repo eigenvalue outright:

| channel | m/g | GMOR | improved | ratio | standard | ratio |
|---|---|---|---|---|---|---|
| mes N=2 | 0.05 | 0.2762 | 0.2593 | **0.94** | 0.0792 | 0.29 |
| mes N=3 | 0.05 | 0.3683 | 0.3764 | **1.02** | 0.0848 | 0.23 |
| mes N=4 | 0.05 | 0.4368 | 0.4504 | **1.03** | 0.0851 | 0.20 |
| mes N=4 | 0.10 | 0.8536 | 0.9283 | **1.09** | 0.3277 | 0.38 |

Improved lands within 3–9% of the continuum law with nothing fitted; standard is
off by 3–5× and worsening as m/g falls. This is the first external anchor at
Table I's own truncation.

**So the published weak-coupling entries are low by roughly 2.5–5×**, because
standard DLCQ converges to the grid artifact this document is about.

Three caveats worth keeping attached. GMOR is a large-N result, so N=4 is the
meaningful comparison and N=2's 0.94 is expected to drift. It is a *meson* law
and says nothing directly about the baryon entries, which carry the largest
corrections and still have no independent anchor. And the improved chiral
exponent is 1.03–1.14 rather than exactly 1, consistent with the few-percent
GMOR gap and worth understanding rather than dismissing.

### The colour weight, and why it is now exact

**Resolved.** `state_sigmas` uses the true `w_ac = ⟨−T_a·T_c⟩`, extracted from
`clfact` — whose `nops=4` path already *is* `Σ_a T^a_{ij}T^a_{kl}` — by calling
it at zero momentum transfer per parton pair. The gate is
`Σ_{c≠a}⟨−T_a·T_c⟩ = C_F`, which `Σ_a T_a = 0` forces on a colour singlet, and
it holds to **0.00e+00** across N = 2…6, B = 0 and 1, L = 2…6, valence and
Fock-extended, on 1,300+ partons. Because the identity holds exactly, the b=0
reduction to `σ_std` is automatic rather than arranged.

Three things had to be right together, and each was found by the gate rather
than by reasoning:

* **the vertex patterns** — H1 for qq, H2 for q̄q̄, and H7's *t-channel*
  structure A for qq̄ (`clfact(1,3,2,0)`, slots `[d̄_ann, d̄_cre, b_ann, b_cre]`).
  H3–H6 are pair creation/annihilation and carry no self-inertia partner; an
  invented qq̄ pattern returned the `−1/N` Fierz piece of a graph that does not
  exist;
* **normalisation** — `clfact` returns an un-normalised element carrying
  `Norm_ii`;
* **multiplicity keyed on `(type, momentum, flavour)`** — `clfact` sums over
  indistinguishable partons, including at the queried momentum. Keying on
  momentum alone passes every baryon and fails exactly the meson states where a
  quark and an antiquark share a momentum but stay distinguishable.

Normalising each parton's weights to sum to `C_F` makes both `Norm_ii` and the
multiplicity cancel, so the implementation needs neither.

**The ansatz turned out to be nearly right, which was not knowable in advance.**
Against the exact operator, `C_F/(L−1)` is off by 0.06% (bar N=3 at m/g = 0.10
and 0.40) and 0.21% (mes N=2 at m/g = 0.10). An earlier revision of this section
called the Fock-extended numbers "undetermined" on the strength of a 144% spread
across weight schemes. That spread is real, but it measures the width of the
*admissible* space, not the uncertainty: the schemes probing it
(momentum-weighted, inverse-momentum) are arbitrary, and the physical answer sits
essentially on uniform. The honest statement is that the space was wide and the
answer was near its middle.

### Pair creation needs no subtraction, and do not "test" that with a K^-2β fit

Only the number-conserving exchange carries the endpoint singularity. `brack`
returns `4/L²` **only when `L + M = 0`**, so momentum conservation pins each
kernel's denominator:

| vertex | kind | written as | forced to equal |
|---|---|---|---|
| H1, H2, H7A | exchange | `k4−k2`, `k2−k4`, `k4−k3` | another *difference* — can reach one grid spacing |
| H3, H5, H7B | pair creation | `k4+k1`, `k1+k3`, `k4+k2` | a *sum* of momenta |
| H4, H6 | pair creation | `k1−k3`, `k2−k3` | `−(k4+k2)`, `−(k1+k4)` — **also sums** |

H4 and H6 look like differences but the constraint pins them to sums, so every
pair-creation kernel is bounded below by a sum of positive momenta. Only the
exchange reaches the grid spacing, and only it has a self-energy partner to
cancel against. That is the whole argument, and it needs no numerics.

**The tempting numerical check is a trap.** Refitting the improved series with a
`K^{-2β}` component allowed, to see whether any endpoint error survives, looks
damning: the residual falls 7–25× and `M(0)` shifts 1.5% in valence sectors but
7.7% (bar N=3, LPN=5) and 22% (mes N=2, LPN=4 at m/g = 0.05) in Fock-extended
ones — exactly where pair creation is active.

It is an artifact. Run the same test in the two-body sector, where the exact
answer is known:

```
1/K only        M0 = 0.995553   resid 4.4e-06   err vs exact  +6.3e-04
+ K^-2b allowed M0 = 0.991663   resid 1.8e-07   err vs exact  -3.3e-03
```

The residual improves **25×** while the answer gets **5× worse**. At small β the
`K^{-2β}` column is nearly constant and so nearly degenerate with the intercept
— the same degeneracy that makes standard DLCQ unextrapolable in the first
place. It absorbs residual and corrupts `M(0)`. Extrapolate improved output in a
**plain 1/K series**, per van de Sande's Eq. (14), and nothing else.

### Historical: why the scalar weight was only ever justified in valence sectors

**The `C_F/(L−1)` weight is derived only where every pair is colour-equivalent.**
For a colour singlet `Σ_a T_a = 0` forces `Σ_{c≠a}(−T_a·T_c) = C_F` per parton,
and when all pairs sit in the same channel that fixes each at `C_F/(L−1)`. That
covers the valence sectors — `LPN=2` meson, `LPN=3` baryon — where there is no
freedom and the results above stand.

It does **not** cover `LPN = valence+2`, which is Table I's truncation for all
thirty entries. There the true weight is the matrix `−T_a·T_c`, and a scalar
cannot be right. Measured, by running weight schemes that all satisfy the b=0
reduction exactly (pair weights summing to 1 per parton, so nothing in the
derivation distinguishes them):

| channel | LPN | m/g | uniform | momentum-weighted | inverse-momentum | spread |
|---|---|---|---|---|---|---|
| bar N=3 | 5 | 0.10 | 2.047 | 0.492 | 3.446 | **144%** |
| bar N=3 | 5 | 0.40 | 7.725 | 5.206 | 9.663 | 58% |
| mes N=2 | 4 | 0.10 | 0.591 | 0.465 | 0.612 | 25% |
| mes N=2 | 4 | 0.40 | 2.628 | 2.583 | 2.646 | 2.4% |

So the Fock-extended numbers are not merely uncertain, they are **undetermined**
— a factor of seven apart between two equally admissible choices at m/g = 0.1.
Note also what this implies about the improved-vs-standard gap on those sectors
(2.7× at bar N=3): that gap is a property of the chosen weight, not a physics
prediction.

Why the sensitivity is so large when `σ_imp/σ_std` never exceeds 1.28: the
interaction annihilates the chiral zero mode, so the lowest eigenvalue is a
small residual of large cancelling terms and a modest change in the diagonal
moves it a long way. That is the same near-cancellation that makes the
weak-coupling limit hard in the first place.

**The prerequisite is therefore the exact `−T_a·T_c`**, which needs no new
colour machinery: `clfact`'s `nops=4` path already *is* `Σ_a T^a_{ij}T^a_{kl}`
(`nct=1` carries the `0.5`, `nct=2` the `−0.5/N`), so what is required is
invoking it at zero momentum transfer per parton pair. The gate is ready and
sharp — `Σ_{c≠a}(−T_a·T_c) = C_F` per parton must come out exactly, and it
already rejected a first attempt whose operator patterns were guessed rather
than read (mesons returned `−1/N` instead of `+C_F`, baryons carried an
uncorrected state norm).

There is also no independently known L ≥ 3 answer, so correctness there will
rest on standard and improved meeting in the continuum rather than on an exact
reference. The two-body sector, which does have one, is what anchors the whole
construction.

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
