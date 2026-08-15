# Open threads

Six, in the order they are worth taking. Each is written to be started cold.
§1 has since been **built**, and is kept because what it measured overturned
several of its own premises — including one this file asserted twice. §4-§6
came out of `weak-coupling-limit.md`. §4 is the published fix and the one to
start with; §5 is the heavier alternative; §6 is a separate question that is
often confused with both.

Everything here was established by measurement; the supporting numbers are in
`performance.md`, `table1-units.md` and `baryon-higher-fock.md`. Note
`performance.md` predates §1 and its profile is now inverted: the eigensolve is
no longer the dominant line, and its claim that 2K = 37 does not fit in memory
is superseded by the table below.

---

## 1. Sparse solver — **built**; what is left is a smaller list

Landed behind flags, dense still the default and still the reference:
`run_python(solver="dense"|"sparse", nev=None)`, `PythonProvider(solver=, nev=)`,
`--solver` / `--nev` / `--policy` on the CLI. `solver="sparse"` is legal only
with `assembly="exact"` **and** `policy="blockwise"` — one cell of six, rejected
early. The cache tag extends only on non-default values, because appending
unconditionally orphans ~7 GB of existing cache.

Reach, measured end to end:

| case | before | after |
|---|---|---|
| 2K = 37, LPN = 0 | OOM on a 15 GB box | 29.9 s, 1.12 GB |
| 2K = 41, LPN = 0 | OOM | 105 s, 4.93 GB |
| 2K = 71, LPN = 5 (Table I truncation) | unreachable | 21 s |

### The premise this section used to carry was wrong

It said **"there is no incremental path — 2K ≳ 48 is the only lever"**, on the
strength of 25–35 → 25–37 giving a median error reduction of 1.00×. That
measurement was made at `n_terms=4` against a Table I that then ran at
`n_terms=2`. At the model in use *at the time*, widening 25–35 → 25–49 improves
**26 of 30** entries, median **1.208×**. (Both now run at `N_TERMS = 4`; the
`n_terms=2` default in `table1_budget` was a defect, fixed in f3c5bf4.)

What replaces it is a sharper statement. The binding constraint at weak coupling
is not K, it is the Eq. (27) basis missing the term that actually controls the
convergence as `a → 0`:

- the endpoint term is `K^{-2a}`, and Eq. (27)'s basis
  `{1, K⁻¹, K^-(1+a), …}` does not contain it. That is van de Sande's result
  ([hep-ph/9605409](https://arxiv.org/abs/hep-ph/9605409)) and it is what
  `weak-coupling-limit.md` is about. The local increment slope suggests `1+a`
  rather than `1+2a`, but the window is too short to separate them.
- so at m/g = 0.1 the tail bound falls only as `K^-0.06`, and **halving the
  bracket would need K × 84,000**
- a different fitting *basis* does not help: the confluent parametrization spans
  the same space and leaves `M₀` identical to 1.1e-14, and refitting in
  `K^{-2a}` powers is circular — measured, see the same doc

**More K is still worth having, but only above `m/g ≈ 0.4`**, where the
remainder goes as `K^-0.55` or steeper and a factor of 2–3.5 in K halves it.
Below `m/g ≈ 0.15` the grid does not resolve the endpoint at all, and no amount
of K closes it. `weak-coupling-limit.md` has the boundary and the per-coupling
factors.

### What the staging turned out to be

The old A → B → C → D order was wrong in three places.

**B was the unlock, not A.** The norm is built pre-weeding and is *exactly*
block diagonal in parton configuration — off-block entries are `0.0`, not 1e-16.
Storing only the blocks:

| 2K | n_pre | dense norm | blocks |
|---|---|---|---|
| 37 | 32,816 | 8.62 GB | 7.56 MB |
| 41 | 83,167 | **55.33 GB** | **27.86 MB** |

Plus the piece of D that was worth keeping: `H0 = D·N` with `D` constant per
block by construction, so `Zᵀ H0 Z` is diagonal and need never be built.

**A is real but modest, and pays only after C.** `eigsh` beats
`eigh(subset_by_index)` 40× at 2K = 41 (3.56 s against 143.80 s) — but on the
*sparse* operator. Iterative methods pay off because of sparsity, not instead of
it. Settled numerics: `which='SA'`, no `sigma`, no folding, `tol=1e-8`, fixed
`v0` (default `eigsh` is non-deterministic, which a cache turns into whichever
run landed first). Matvec count does not grow with n — 220 → 432 across
n = 193 → 15,235.

**D's headline is dead.** `hnu = ZᵀHZ` is only **4.3%** denser than `H`
(3.079% against 2.952% at 2K = 41). The triple product does not fill in, so
matrix-free saves nothing and adds two sparse products per matvec.

**`n_post` is far smaller than assumed** — the ratio is ≈ 8/(2K), not 0.28. At
2K = 41 it is 15,235, so the dense objects to size for are 15,235², not 23,000².
Meson 2K = X has the same `n_post` as baryon 2K = X+3, so "meson 48" is
"baryon 51" — five steps past what now runs, not one.

### Not built, deliberately

**C2, candidate enumeration.** The **|ΔL| ≤ 2** rule this section used to
recommend admits **81–86%** of matrix entries at 2K = 21–31, against true
densities of 26% → 8.5% — a 1.2× filter that does not improve with K. The rule
that does work: every surviving parton must contract with one of identical
`(type, momentum, flavour)`, so the two states share a residual key — a hash
join on "state minus j partons", `j ∈ 0..3`. Verified a superset against every
non-zero of the dense build at 2K = 21/23/25/29/31, overshoot 1.04–1.06×.

It is simply not needed: the existing bare scan costs 20–26 ns/pair, 1.36 s at
2K = 39 out of ~90 s. Revisit only if a profile puts it above 20% of a run.

### Constraints that must survive

- The sparse path supports **`assembly="exact"` + blockwise `Z` only**. The
  Fortran-compatible diagonal-only assembly is inherently dense and basis
  dependent; it stays on the dense path for small-K cross-validation.
- **Matrix element values must not move.** `tests/test_kernels.py` carries seven
  `array_equal` assertions plus one `assert_array_equal`, and
  `tests/test_fortran_python.py` asserts the norm is bit-identical (`== 0.0`) to
  the 1990 Fortran. Storage, blocking and traversal order across independent
  rows are unconstrained; **reassociating the arithmetic is not** — which is why
  the blockwise gate is bit-equality on stored entries but a *tolerance* on
  eigenvalues: sparse `Z` sums the triple product in a different order, worth
  ~1e-13 against this algorithm's own 1e-4 reproducibility floor.
- `structure_function` needs `norm @ c` as a matvec, and the weight is
  `c * (Nc)`, **not** `c²` — `c²` silently breaks every sum rule. The blockwise
  norm is therefore a `BlockDiagonal` supporting `@`, not a `csr_matrix`.

### A correctness bug found while measuring

`MXTRM = 12552` capped the colour-contraction tables in **both** the numba
kernels and the interpreted path, and neither raised. Because both shared the
cap they agreed bit-for-bit, so backend-equality tests **structurally could not
detect it**. At N=4, B=1, 2K=20 the assembled norm had an eigenvalue of −1.8e5,
where a Gram matrix cannot go below 0.

Both guards now raise and `MXTRM` is 200000, which makes the sums exact —
N=4/2K=20 now matches `tools/colour_norm.norm_bruteforce` with difference **0**.
Nothing published was affected (Table I and Figs. 7/8 go through `sweep_lpn`,
Figs. 4–6 are N=3 at 2K ≤ 29). The regression test is that the norm is positive
semidefinite per config block — which catches this class where backend equality
provably cannot.

---

## 2. QCD + QED with Standard-Model charges

### The structural fact that makes it cheap

`clfact` already computes **both** Fierz contractions, from
`Σ_a T^a_{ij}T^a_{kl} = ½(δ_il δ_jk − (1/N) δ_ij δ_kl)`:

```python
nctf = 2 if nops == 4 else 1
if nct == 1:  reslt *= 0.5        # δ_il δ_jk
else:         reslt *= -0.5/N     # δ_ij δ_kl
```

For U(1) the generator is the charge, so `Σ T T → Q_a Q_b δ_ij δ_kl` — **the
same index structure as the −1/N term.** QED is therefore a reweighting of a
contraction already being computed, not a new one:

```
coeff(nct=2) = −0.5·α/N + c·(e²/g²)·Q_a·Q_b
```

`α = 1` gives SU(N) and `α = 0` gives U(N) — the paper's own parameter under
Eq. (24) ("For SU(N), α=1, and α=0 for U(N)"), so the abelian case is already
anticipated in the formalism.

`H = H_QCD + (e²/g²)·H_QED` is **linear** in the coupling ratio, so both
matrices are built once and every ratio is one `eigh`. QCD enters at O(e⁰) and
QED at O(e²), so they separate **algebraically** rather than statistically.

### Changes

1. `Params`: `chg[NF]` alongside `rmq`/`iflv`, plus one `e2g2` ratio. In 1+1
   `g²` has mass dimension, so `e²/g²` is dimensionless and is a genuine free
   parameter — it is *not* fixed by α_em.
2. `clfact_nb`: the coefficient above.
3. `hamqcd_nb`: pass the two currents' charges. **The fiddly part** — the
   flavour pairing differs per vertex (`_V4_FL` is `[0,1,0,1]` for four and
   `[0,1,1,0]` for two), so the charge product must follow the actual
   contraction, not the slot order. Derive once, test per vertex.
4. Fix `c` by matching an exact result, not by convention-chasing.

### Validation ladder

1. `e2g2 = 0` reproduces current results **bit-identically** — the existing
   `array_equal` tests give this for free.
2. **Schwinger model**: `α = 0`, `N = 1`, massless → `M² = g²/π` exactly, i.e.
   **exactly 1.0** in the paper's `M²/(m²+g²/π)`. Validates the abelian path
   *and* fixes `c`. N=1 already generates states (16 at 2K=12 — the neutral
   qq̄ sector, which is the right one).
3. **Isospin restoration**: `m_u = m_d`, `Q_u = Q_d` → recover the
   flavour-symmetric degeneracy already tested.
4. Sum rules unchanged (charge-blind).
5. The e²-linear part of any mass must vanish for a neutral, isospin-symmetric
   state.

### Two landmines, both already documented

- **`qcdf.f:253` overwrites `iflv(1)` with `N·B`** (`flavour.md`). Charge
  sectors hit this immediately.
- **Zero modes.** Antiperiodic BCs exclude k=0, where the θ-vacuum and chiral
  condensate live. You get the spectrum, not the vacuum structure.

### Charged states

In 1+1 the Coulomb potential is linear, so an isolated charge has energy ∝ L —
charged states genuinely diverge in infinite volume. That is not a
compactification artifact and no boundary condition removes it. **But** the
standard lattice prescription (QED_L: delete the photon zero mode) makes them
finite-volume-computable, and this code *already* drops the zero-momentum
exchange — `brack(L, M)` returns 0 for `L == 0`. So a charged sector is
probably reachable.

The cost is that QED_L finite-volume artifacts are a series in 1/L — the same
variable as the DLCQ extrapolation. That is **not** fatal: both vanish as
K → ∞ and extra fitted terms absorb them. The real risk is exponents outside
the Eq. (27) basis (log K / K in particular), which is a paper-and-pencil check
on the light-cone constraint equation `∂_-² A⁻ = g J⁺`, worth doing before
writing code.

With SM charges, Gauss's law admits the **neutron** (udd, Q=0) and the **π⁰**,
but not the proton or π⁺, unless the zero-mode subtraction above is used.

---

## 3. Fig. 6(a)–(c), the five-quark series at 2K = 21

**Closed as far as the surviving artifacts allow — do not reopen it as a
reproduction problem.** See `baryon-higher-fock.md` §2.1 for the accounting.

What is settled: the published five-quark points are not reproducible, and the
failure is not on our side in any way a better reading or a different convention
could fix. They exceed a bound our amplitudes cannot violate (§2.1.2); they are
incompatible with the valence points printed in the same panel, since
`c₅ = −(H₅₅ − wN₅₅)⁻¹H₅₃c₃` holds to 5.8 × 10⁻⁴ (§2.1.3); all three panels fail
the same way while their valence curves reproduce to 0.6–2.9% (§2.1.4); and it
recurs at a second coupling in thesis fig. 15(a) (§2.1.5).

**The defect is downstream of the eigensolve** (§2.1.6). The preserved 1990 run's
own eigenvector reproduces the published *valence* to 1% and sits 43.2% from the
published five-quark curve — so the colour sums, the Hamiltonian build and the
diagonalization are all exonerated, and with them any story about numerical
instability in the solve. What is left is the step from eigenvector to plotted
curve: the lost `wf`/`wfbig`.

What is not settled: **what that program did.** It cannot be run, and no
intermediate version survives.

**If you pick this up**, the one high-value move is archival, not computational:
find `wf`/`wfbig`, or an intermediate `qcdf.f`. Everything computational has
been tried; §2.3 is the list.

### Published data still unused

Found while chasing this, and worth more than it cost. The thesis prints far
more than the article, and most of it has never been checked:

| source | what it is | why it is worth doing |
|---|---|---|
| **thesis fig. 22 + Table 5** | **SU(2)** meson four-quark contribution at 2K = 14 and 20 | the only published K-scan of a higher-Fock structure function anywhere in this work. It measures the K-dependence of the conversion instead of inferring it. Even-channel, so it cannot settle the odd-channel question above — but it is the only direct probe there is. **N = 2 only** |
| **thesis figs. 19, 20** | four-quark contribution as a function of **N** — fig. 20 is **U(N)**, not SU(N) | an untouched published *N*-scan and U(N)-scan, not the K-scan an earlier revision of this table claimed. Corrected after reading the thesis pages: the captions read "as a Function of N", and the body text puts the K-scan in fig. 22 alone ("at 2K equal to 14 and 20 for both the weakly and strongly coupled SU(2) meson"). The thesis OCR renders both N and K as `A'`, which is how the two got conflated |
| thesis figs. 14, 16 | meson 4q and two-baryon 8q at **strong** coupling | free regression targets; the weak-coupling versions already reproduce |
| thesis fig. 15(b), (c) | 2nd and 3rd baryon, five-quark, m/g = 0.1, 2K = 21 | two more states at the second coupling; 15(a) is done (§2.1.5) |
| thesis figs. 13(b)–(e) | 2nd–5th B = 2 states | already flagged in §1.1 |
| thesis figs. 21, 23, 24 | momentum splitting, incl. **SU(3) baryons** (fig. 24) | untouched; fig. 24 is the only other baryon higher-Fock content in the thesis |

`literature/` now also holds **SLAC-PUB-5811** (= Brodsky, McCartor, Pauli &
Pinsky, *Particle World* **3** (1993) 109), whose figs. 2 and 3 redraw this
sector cleanly at both couplings — see §1.3 item 6. Note CERN's copy is behind
an anti-bot challenge and the SLAC copy's TLS chain does not validate against
older CA bundles; fetch it in a browser.

### Qualifications

Two, both of which cost earlier revisions of this file some credibility:

- The f(k) reweighting is coupling-independent in **sign** but not in magnitude
  (1–4% at k = 3, 5; 15–37% elsewhere). Do not quote it as a fixed function.
- Every "higher K reproduces" data point is **even-channel** (B = 0 or 2), and
  `K_code ≡ L (mod 2)` locks parton parity to K parity, so the only within-B=1
  comparison is 2K = 15 vs 21. A K-growing effect confined to the odd channel is
  not excluded by published data — no baryon-versus-K panel was ever printed.

Do **not** re-derive from the drawn curves. They are cubic splines through very
few points and the thesis says so itself; panel (c)'s valence curve has two
spline "zeros" that are not in its own data. Work from the markers.

---

## 4. Improved DLCQ — the published fix, and it applies to the DLCQ path

**Start here, not with the basis rewrite below.** Van de Sande
([hep-ph/9605409](https://arxiv.org/abs/hep-ph/9605409)) diagnosed this exact
problem in 1996 and published a cure that modifies the *Hamiltonian* rather than
the basis, so it keeps everything DLCQ is good for — orthogonal basis, sparse
matrix, matrix elements with no integrals to evaluate.

The move is to add and subtract a term so the kernel vanishes when both momenta
are near an endpoint. His Eq. (12) replaces the off-diagonal sum with

```
g² Ψᵢ I(i/K) + g² K Σ_{j≠i} [ Ψᵢ ((j(K−j))/(i(K−i)))^β − Ψⱼ ] / (i−j)²
I(x) = ∫dy [1 − (y(1−y)/(x(1−x)))^β] / (x−y)²
```

with `β` the endpoint exponent — our `a`, already available exactly from
`units.endpoint_exponent`.

**Why this is promising for precisely our case.** Van de Sande notes that for
small μ improved DLCQ is near-exact for the *lowest* state already at K = 10,
and calls it "a fortuitous accident: for small μ, the lowest eigenfunction is
Ψ(x) ≈ x^β(1−x)^β," so the correction term nearly annihilates. That accident is
exactly Table I — the lightest state at weak coupling. The entries this repo
cannot currently determine are the ones the fix should nail. (Excited states
improve too, but by less; see his Fig. 4.)

**What is genuinely new work.** His derivation is for the two-body 't Hooft
equation at large N. Ours is a multi-particle Fock space at finite N with
baryons, so the endpoint weight has to be generalized — each parton carries its
own `x`, and what plays the role of `((j(K−j))/(i(K−i)))^β` for an L-parton
state needs deriving. The change lands in the instantaneous/Coulomb matrix
elements, i.e. `hamqcd` and `clfact` in `python/qcdf.py`, which is the hot path.
He claims the method "can be easily applied to other theories and large
numerical calculations" but demonstrates only on 't Hooft.

**Do not** try to shortcut this by refitting existing data in a `K^{-2β}` basis.
That is measured, circular, and documented under "Do not use the convergence law
as a fitting basis" in `weak-coupling-limit.md`.

## 5. An `x^a`-adapted basis — the heavier alternative

Read `weak-coupling-limit.md` first; this entry is its consequence. The
weak-coupling column is limited by the momentum grid failing to resolve
`φ ~ x^a` near the endpoints, where the relevant scale is `x ≲ e^{-1/a}`
(`5e-11` at `m/g = 0.05`). Nothing about K, the fit form, or the error budget
addresses that — six approaches were tried and measured dead. The fix is to
*represent* the endpoint rather than resolve it: build `x^a` into the basis.

This is established practice elsewhere. Anand, Fitzpatrick, Katz and Xin
([arXiv:2111.00021](https://arxiv.org/abs/2111.00021)) restore convergence at
quark masses far below the strong-coupling scale in lightcone conformal
truncation by exactly this move, modifying the basis according to the quark
mass following 't Hooft's endpoint analysis.

The hard part here is that DLCQ's basis is a *momentum grid*, not a function
basis, so this is not a drop-in: it changes what a "state" is. And note that the
graded-mesh alternative is **not** available on the DLCQ path — `x = k/K` follows
from quantizing momentum in a periodic box, not from a discretization choice.

So start in `dlcq/thooft.py`, where a mesh *is* a free choice. Two steps, in
order:

1. **Graded mesh.** Replace the uniform grid with a geometric one reaching
   `x_min ~ 10⁻²⁵`. The cell-wise exact treatment of the principal value
   generalizes directly — `P∫ dy/(y−x)²` over a cell is
   `1/(x−y⁺) − 1/(x−y⁻)` for any cell width — and the diagonal is still fixed by
   the sum rule that the kernel annihilates a constant. `weak-coupling-limit.md`
   costs this: ~571 nodes at `m/g = 0.05` against `4e23` for a uniform mesh, and
   2284 nodes at `m/g = 0.0125` against `3e94`. Expected outcome: `M²` plateaus
   under refinement where it currently climbs ~7% per doubling, and the chiral
   exponent falls from 2 toward 1.
2. **`x^a`-weighted basis**, if step 1 confirms the mechanism. ~30 coefficients,
   spectral convergence, endpoint exact by construction.

That solver is 120 lines, has an exact benchmark (0, 5.8817, 14.1429 at
`m = 0`), and its failure to converge is documented and reproducible — so step 1
is cheap and falsifiable. **Do it before touching the DLCQ path.**

## 6. Zero modes — real, measurable, and not what §4 or §5 is about

Keep this separate from §4 and §5. Zero modes are *not* the cause of the
weak-coupling problem: the continuum solver in `thooft.py` has no zero-mode truncation and
reproduces the same artifact. But they are a genuine omission — DLCQ here uses
antiperiodic quarks and drops the gluon zero mode — and the effect is not small.
Müller, Kalloniatis and Pauli
([hep-th/9803204](https://arxiv.org/abs/hep-th/9803204)) find **21% shifts** in
the lowest bound-state masses of 2D SU(2) Yang–Mills with adjoint scalars once
zero-mode wavefunctions are included.

It is more tractable than the naive picture suggests. Do not add `k = 0` to the
Fock basis — at fixed `K⁺` that is infinite-dimensional. The dynamical zero
modes are a few *collective* variables; drop the constrained ones and they are
governed by an infinite square-well potential, whose wavefunctions tensor onto
the existing basis, multiplying the matrix dimension by the levels retained.
Formalism: Pauli–Kalloniatis–Pinsky
([hep-th/9403038](https://arxiv.org/abs/hep-th/9403038),
[hep-th/9509020](https://arxiv.org/abs/hep-th/9509020)).

What is genuinely new work: the published implementations are SU(2) pure glue
and adjoint scalar matter, **not** SU(N) with fundamental quarks and baryon
number, so the zero-mode/non-zero-mode Hamiltonian matrix elements have to be
derived for this case. Gauge fixing needs revisiting too — `A₋ = 0` is not
available. The eigensolver, colour algebra and state generator all survive
unchanged, and the quark antiperiodicity can stay.

---

## A method note worth keeping

An apparent disagreement between two extrapolations can be **entirely** an
artifact of fitting different K windows. Refitting on the other side's window is
cheap and should be the first check whenever a published extrapolation fails to
reproduce. It is what dissolved most of the baryon N=4 discrepancy at
m/g = 1.6 (see `table1-units.md`).
