# An equal-time calculation: what it would take

## Why

Every finite-`N_c` determination of the 2d chiral exponent in existence — this
repo's two channels, Anand–Fitzpatrick–Katz–Xin (arXiv:2111.00021), Kochergin
(arXiv:2405.04031) — is quantised on a null plane with a trivial vacuum. They all
give `alpha -> 1`. The one equal-time measurement, Bañuls, Cichy, Cirac, Jansen &
Kühn (PRX **7**, 041046, arXiv:1707.06434), gives `nu = 0.700(31)` i.e.
`alpha = 1.40`, which is 1.1σ from bosonization's 4/3 and 6.4σ from 1.

That single shared assumption is exactly what the anomalous exponent is about, so
"unresolved" is the honest verdict and an independent equal-time calculation is
the thing that would settle it.

**Two things make it worth doing rather than citing.**

1. **Bañuls' window is a factor of 4** (`m/g ∈ [0.1, 0.4]`), not a decade.
   Separating 0.700 from 0.667 from 0.500 on that lever arm is what limits them.
2. ~~No additive mass renormalisation appears to be applied... plausibly
   explains 0.700 against 2/3.~~ **Retracted — see "Corrections" below.** The
   SU(N) additive shift is *zero*; the sign of such a bias runs the wrong way
   in any case (a positive shift pushes the measured exponent **down**, and
   faking 0.700 from a true 2/3 would need `δ < 0`, i.e. `ag < 0`); and Bañuls
   extrapolate `ag → 0` at fixed bare `m/g` anyway, which removes an `O(g²a)`
   shift by construction. **Their 0.700 needs a different explanation, and this
   was not it.**

And it measures **both** observables with one machinery, so it tests this repo's
meson/baryon ratio result (`2 sin(pi nu / 2)` to 0.4%) as well as the exponent —
not just the half we disagree with.

## What already exists, and the gap it leaves

The model has been simulated three times. **None of them computed the thing in
question**, and none ships code.

| work | method | what it covers | quark masses |
|---|---|---|---|
| Bañuls et al., arXiv:1707.06434 | MPS, SU(2), N_f=1 | the only published chiral-exponent fit, `nu = 0.700(31)` | `m/g ∈ [0.1, 0.4]` — a factor of **4** |
| Fujikura & Hidaka, arXiv:2605.17183 | **VUMPS**, SU(2), N_f=1 | `c_IR = 1.04` at `m = 0`; baryon masses | `m/g = 0, 0.5, 1.0` |
| Hayata, Hidaka & Nishimura, arXiv:2311.11643 | MPS, SU(2) **and SU(3)**, N_f=1 | finite-density EoS, condensate vs density | `m/g = 0.5, 1.0` |

So the chiral region is essentially untouched: the only non-zero masses ever run
below `m/g = 0.1` are Bañuls' four-fold window, and the only `m = 0` point is a
CFT identification (`c_IR`), not a scaling measurement. **Nobody has run a decade
in `m/g`**, which is exactly where the crossover question lives.

### Public code: there is none for the SU(N) basis

Checked directly rather than inferred from citations:

| repo | state |
|---|---|
| `falquez/GCB_SUN_LGT` — "Gauss constrained basis for SU(N) LGT with Tensor Networks" | **not an implementation.** Six files: README, two plotting scripts, two data files, a PDF. The 1.3 MB is the paper. |
| `epv180502/openSU2LGT` | Julia, 0 stars, first pushed 2026-08-11. Too new to rely on; worth reading. |
| `gcataldi96/ed-lgt` | **404** — gone or private. Only the Zenodo archive (record 11145318) survives. ED only in any case. |
| `tenpy/tenpy` | healthy and active (493 stars). General library, no gauge-theory model. |

**A better basis than the one above: loop-string-hadron (LSH).** Raychowdhury &
Stryker, **arXiv:1912.06133**, PRD **101**, 114502 — a *local*, manifestly
gauge-invariant description in terms of loop, string and hadron degrees of
freedom. It applies to **SU(2) and SU(3)**, to **open and periodic** boundaries,
and in 1+1d and higher, and it has been built into an MPS ansatz already:
Mathew, Gupta, Kadam, Bapat, Stryker, Davoudi & Raychowdhury, arXiv:2501.18301
("Tensor-network toolbox for probing dynamics of non-Abelian gauge theories"),
plus arXiv:2603.24698 for (1+1)d SU(2) string breaking and **arXiv:2212.04490** (PRD 107, 094513)
giving (1+1)d SU(3) LSH with dynamical quarks.

This matters because LSH is **local**: it removes the non-local recoupling that
makes the Gauss-constrained irrep basis a week of hand algebra, and it gives a
route to SU(3) — i.e. to N_c = 3, which is where Table I and AFKX live — instead
of stopping at SU(2). It carries a bosonic occupation cutoff to extrapolate in,
which the exact-elimination basis does not, but that is a far more tractable
systematic than hand-written 6j symbols.

**And TeNPy already implements VUMPS** (`tenpy.algorithms.vumps`,
`SingleSiteVUMPSEngine` / `TwoSiteVUMPSEngine`), so the calculation can be done
directly in infinite volume. That **eliminates rather than mitigates** the two
open-boundary concerns: the left-boundary electric-field choice (which is a
physical background-field/θ choice, not a convention) and the finite-chain
sensitivity of the baryon, which is a whole-chain colour singlet and therefore
the observable most exposed to boundary effects. Fujikura & Hidaka use VUMPS on
exactly this model for exactly this reason. Alternatives if needed:
`MPSKit.jl`, `ITensorInfiniteMPS.jl` (both Julia).

### Verified by reading the papers, not the search results

* **No code is released by any of them.** arXiv:2501.18301 has no code- or
  data-availability statement in body, footnotes or acknowledgements; its only
  tooling remark is *"We rely upon existing methods defined in the ITensors.jl
  package [49] to construct and optimize the MPS and MPOs needed for this
  study."* Same for Bañuls and for Fujikura & Hidaka.
* **LSH is nearest-neighbour**, quoting arXiv:2501.18301 directly: *"The LSH
  Hamiltonian comprises at most nearest-neighbor interactions, allowing a compact
  matrix-product-operator (MPO) representation."* That is the decisive fact for
  cost — it makes the model a stock TeNPy `CouplingMPOModel`, days of work rather
  than a week of hand-written 6j recoupling.
* **Its cutoff is one parameter**: *"the infinite-dimensional Hilbert space
  associated with the gauge bosons must be truncated... restrict the loop
  quantum number to `n_l(r) ≤ n_l,max`, leading to a local physical
  Hilbert-space dimension of `4(n_l,max+1)`."* A single extrapolation, against
  the irrep basis's non-local recoupling.

### The three choices, which are not alternatives

`LSH` is a *basis*; `VUMPS` is an *algorithm*. The plan uses both:

| choice | decision | why |
|---|---|---|
| basis | **LSH** | nearest-neighbour ⇒ compact MPO; extends to **SU(3)**; published MPS precedent. Costs one cutoff `n_l,max` |
| algorithm | **VUMPS** (`tenpy.algorithms.vumps`) | infinite volume by construction, so the boundary-field choice and the baryon's whole-chain sensitivity never arise |
| build vs borrow | **build**, after asking | nothing public exists; but LSH being local makes building cheap |

### Decided

* **Build it here.** No approach to the authors — nothing public exists, and LSH
  being nearest-neighbour makes building it cheap enough that borrowing is not
  worth the dependency on someone else's timeline.
* **`physics-tenpy` is an accepted dependency**, the first beyond
  numpy/scipy/numba/matplotlib/h5py/pillow/pytest. It carries both the DMRG and
  the VUMPS engines, so one dependency covers finite and infinite volume.


## Corrections from the derivation pass

Six errors were found by an adversarially-checked derivation round, four of them
in material a coder would have typed. The three load-bearing ones were
re-verified independently here.

### 1. An OBC chain needs boundary Gauss constraints, not just the internal one

A finite open chain has no link off either end, so beyond the internal Abelian
Gauss law it needs

```
N_R(r=1) = 0        and        N_L(r=N) = 0
```

Without them the sector contains states with a unit of dangling colour flux at
each open end — unphysical static sources charged **zero** electric energy,
because `H_E` sums only internal links. An AGL penalty does not catch them:
their internal defect is exactly zero.

Measured against a truncation-free gauge-fixed reference at N=4, x=μ=1: the
AGL-only sector has dimension 206 against the correct 20, and its first excited
state is −3.767 against the true −2.856 — **a gap 69% of the right one**. At
N=6, `x=4, μ=0.2` (i.e. `m/g = 0.05`, this project's target corner) the gap is
2.745 against 3.490.

**The ground state is unaffected**, since it satisfies the boundary condition
anyway — which is exactly why this is dangerous: every `E_0` validation passes.

**Mitigating:** the production path is VUMPS in infinite volume, where there is
no boundary and none of this arises. It bites only the ED/finite-DMRG validation
ladder — which is where the plan already puts the gate, so it would have been
caught there.

### 2. The electric term is easy to write twice too large

`H_E = (g²a/2) Σ N̂_L(N̂_L/2 + 1)` glossed as "= C₂(j) = j(j+1)" is false:
with `j = N_L/2`, `N_L(N_L/2+1) = 2j(j+1)` — verified exactly for
`N_L = 1..5`, ratio 2.000 throughout. The correct form (arXiv:1912.06133 Eq. 41;
arXiv:2501.18301 Eq. 1) carries **¼** per site-end, giving `(g²a/2) Σ j(j+1)` on
the AGL surface.

**Why it matters:** doubling `H_E` alone is identical to running at
`g'² = 2g²`, so every `(ag, m/g)` label is wrong by `√2` — a run intended at
`(0.10, 0.010)` is physically at `(0.1414, 0.00707)`. The *exponent* survives a
pure rescaling of the x-axis, so it would not fake `α`; the amplitude, the
reproduction of Bañuls' window, the `(ag)²` extrapolation, and the U(N)-vs-SU(N)
null test all break silently.

### 3. The additive mass shift is zero for SU(N), and the sign argument was backwards

Three independent reasons the original framing was wrong:

1. **There is no shift to subtract.**
2. **The sign runs the wrong way.** With `M = A(m+δ)^ν` the local slope is
   `ν·m/(m+δ)`, so a *positive* shift pushes the measured exponent **down**.
   Faking 0.700 from a true 2/3 needs `δ < 0`, i.e. `ag < 0`. Verified:
   at `m = 0.02`, `δ = +0.005` gives slope 0.533, not 0.700.
3. **Bañuls remove it anyway**, extrapolating `ag → 0` at fixed bare `m/g` (their
   App. E, "a polynomial in `1/√x` up to second order").

So **Bañuls' 0.700 needs a different explanation.** Live candidates, none
examined: the short lever arm (a factor of 4 in `m/g`), residual `j_max`
truncation, contamination of the vector by the near-degenerate diquark — or it
is simply right and 2/3 is not the answer.

### 4. The `−2/π` free-fermion check, right number, wrong reasoning

Setting `U → 1` gives two *decoupled* tight-binding chains with hopping
coefficient 1, so each half-filled chain gives `−2/π` and a colour **pair**
gives `−4/π`. Bañuls' published `ω₀ = E₀/(2Nx) → −2/π` is nonetheless right,
because the factor 2 is already in their denominator. Stated as "two colours
give `−2/π`", a correct implementation fails this check by exactly 2 — on the
test billed as the factor-2 audit.

### 5. Provenance, and it was my error

`arXiv:1812.07554` is *"Solving Gauss's Law on Digital Quantum Computers with
LSH Digitization"* — circuits, no Hamiltonian. The LSH Hamiltonian and operator
algebra is **arXiv:1912.06133** (PRD 101, 114502). Both appeared as separate
hits in the search that produced this plan, and I merged them. Likewise the
(1+1)d SU(3) LSH with dynamical quarks is **arXiv:2212.04490** (PRD 107,
094513); `arXiv:2407.19181` is the trivalent-vertex `d>1` paper and carries an
unsolved outer-multiplicity problem that does *not* affect 1+1d.

## The physics being discretised

Kogut–Susskind, one spatial dimension, staggered fermions, open boundaries:

```
H = -(i/2a) Σ_n [ψ†_n U_n ψ_{n+1} − h.c.]  +  m Σ_n (−1)^n ψ†_n ψ_n
    + (g²a/2) Σ_n E_n^a E_n^a
```

In 1+1d with open boundaries **Gauss' law determines every link from the matter
to its left**, so the gauge field can be eliminated exactly — no electric-flux
truncation, no `j_max` extrapolation. That is what Bañuls' "efficient basis" and
the gaugeless formulation of arXiv:2511.00154 both exploit.

**The one hard part, and it is specific to non-abelian.** For U(1) the
elimination is trivial: `E_n = Σ_{k≤n} q_k`, and the electric term becomes a
long-range Coulomb interaction in the charges — a standard MPO. For SU(N) the
link carries an *irrep*, not a number: `E_n^a E_n^a = C₂(j_n)` where `j_n` is the
total colour representation of sites `1..n`, which depends on the **recoupling
path**, not merely on the total charge. So the irrep label has to live in the
basis, and the state is `(fermion occupations, j_0 … j_L)` with Clebsch–Gordan
weights on the hops.

**Consequence for tooling:** TeNPy's `np_conserved` implements *abelian* charges
(U(1), Z_N). Non-abelian symmetric tensors are not available, so the SU(2) irrep
structure has to be carried explicitly in the local Hilbert space rather than by
the symmetry machinery. That is what makes Phase 2 a real piece of work and
Phase 1 nearly free.

## Phase 1 result: the pipeline reproduces `1/sqrt(pi)` to 0.06%

Twelve points at fixed physical volume `Lg = L/sqrt(x)`, `Lg` in {30,45,60},
`x` in {16,36,64,100}, `L` up to 600, plus `x` = 144 and 196 at `Lg = 45`.

**Continuum `M/g = 0.56451` against `1/sqrt(pi) = 0.56419`, +0.06%.**

Both limits are controlled, and both forms were *measured* rather than assumed:

* **volume: `A + B/Lg^2`.** Successive differences at fixed `x` have ratio
  0.354, 0.353, 0.353, 0.353 across all four `x` — against 0.350 predicted by
  `1/Lg^2` and 0.500 by `1/Lg`. This is just the boson's momentum quantisation
  in a box, `sqrt(M^2 + (pi/L)^2) ~ M + pi^2/(2 M L^2)`.
* **spacing: `A + B/sqrt(x)`**, the staggered form. It fits with max residual
  3.8e-5 where `1/x` fits 70x worse, and the intercept is stable as the range
  extends: 0.56878, 0.56875, 0.56873 at `Lg = 45` for `x` up to 100, 144, 196.

Systematics checked and excluded: **MPS truncation** — `M/g` = 0.609922,
0.609919, 0.609919, 0.609919 at `chi` = 60, 100, 160, 240, converged to six
digits; and **the `x` range** — extending to `x = 196` moves the extrapolation by
5e-5.

**Two wrong readings on the way there, both recorded because the failure modes
recur.** An earlier scan varied `L` and `x` independently, which put the largest
finite-volume error on the *finest* lattice and made two estimators of the same
limit disagree by 9%; its apparent 0.36% agreement was a coincidence of two
uncontrolled extrapolations. Then, having fixed the volume, fitting `A + B/Lg`
instead of `A + B/Lg^2` over-corrected to −0.86%, and the power-law drift was
misread as an open-boundary pathology requiring VUMPS. It was ordinary
particle-in-a-box quantisation. **The form of an extrapolation is a measurement,
not a modelling choice** — in both cases the data said which power, and reading
it off took one line.

So the conventions are confirmed end to end: `x = 1/(ag)^2`,
`mu = 2(m/g)sqrt(x)`, `L_n = sum_{k<=n}[S^z_k + (-1)^k/2]`, and
`M/g = dW/(2 sqrt(x))` — the last of which no fixed-`x` analytic check can see,
and which is what this scan was for.

## The two-flavour test found a state-targeting problem, not the exponent

The two-flavour capability test — does the pipeline resolve an *anomalous*
exponent, `M ~ A m^{2/3} g^{1/3}` — came back measuring the wrong state.

Four points at `x=16, Lg=25, L=100`, `m/g` = 0.4, 0.2, 0.1, 0.05:

| m/g | M/g | local exponent |
|---|---|---|
| 0.4 | 2.47845 | |
| 0.2 | 1.72350 | 0.524 |
| 0.1 | 1.34856 | 0.354 |
| 0.05 | 1.15819 | 0.220 |

The local exponent *decreases toward zero*, so `M/g` tends to a constant rather
than vanishing: `M/g = M0 + c(m/g)` fits with `M0 = 0.970`, and successive
differences halve exactly (ratios 1.97, 2.01). That is a state which stays
massive in the chiral limit — the **flavour-singlet eta**, continuum mass
`sqrt(N_f/pi) = 0.798` (at `ag = 0.25` a ~20% lattice shift is expected). The
massless pion would need `M/g -> 0`; at `m/g = 0.05` the 2/3 law predicts 0.273
against 1.158 measured.

**This is not a Hamiltonian error.** The `N_f = 1` reduction reproduces the
one-flavour build to `0.00e+00` and ED to the same, so the flavour sum, the
staggered background charge `N_f(-1)^k/2` and the constant are all right.

It is a **state-targeting failure**, and a general one:
orthogonality-constrained DMRG converged to the eta even at `m/g = 0.4`, where
the pion should be lighter (1.09 against 2.48). Orthogonality to the ground
state is not enough to reach a state in a different flavour sector when the
Hamiltonian, the initial product state and the conserved charge are all flavour
symmetric. This is exactly the contamination Bañuls needed an explicit penalty
term for, and it is why their vector-mass extraction required one.

**Diagnosed, then re-diagnosed — the first diagnosis was wrong.** Four controlled tests at
`m/g = 0.1` (target `M/g = 0.4326`, baseline 1.25452):

| variable changed | M/g |
|---|---|
| chi 90 → 200 | 1.25450 |
| volume Lg 25 → 50 (L 100 → 200) | **1.25452** |
| lattice x 16 → 36 | **0.00000** |

Doubling the volume changes the answer by *nothing* — identical to seven digits,
which no propagating particle's mass can be — and the finer lattice collapses the
gap to exactly zero. Both point the same way, and a direct check confirms it:

```
charges="same"        : vacuum Q=[0]    pion Q=[0]
charges="independent" : vacuum Q=[0 0]  pion Q=[0 0]
per-flavour Sz        : vacuum (0, 0)   pion (+1, -1)
```

That check was made with the wrong call. `get_total_charge()` returns the
*gauged* charge, and for finite boundary conditions the total can be gauged away
into the trivial end legs. With `only_physical_legs=True`:

```
charges="same"        : vacuum Q=[0]    pion Q=[0]      same sector
charges="independent" : vacuum Q=[0 0]  pion Q=[2 -2]   DIFFERENT sectors
```

So the sector targeting **did** work, and the bookkeeping was never the problem.

**The actual bug was site parity.** `_pion_state` anchored its two flips to
`L//2`. Even sites are empty in the staggered vacuum and odd sites filled, so the
flips must land on an even and an odd site respectively — which holds only when
`L//2` is even. At `L = 100` it is (site 50); at `L = 150` it is not (site 75,
already filled), so **both flips silently became no-ops and the "pion" state was
literally the vacuum**, giving a gap of exactly 0.00000. Anchoring to parity
instead (`c = 2*(L//4)`) gives `Q = [2, -2]` at every `L`, and `x = 36` returns
1.09327 rather than 0.

With that fixed the state responds to the lattice: `M/g` = 1.25452 at `x = 16`
and 1.09327 at `x = 36`, decreasing as the spacing falls, which is the right
direction for a lattice artefact. Two points extrapolate in `1/sqrt(x)` to ~0.77
against the target 0.4326 — not there, but no longer inert.

**The volume-independence was also misread.** Seven-digit agreement between
`L = 100` and `L = 200` is exactly what a *converged massive state* should give:
corrections go as `exp(-M Lg)`, and `M Lg ~ 31` here makes them `~1e-14`. It was
evidence the calculation was working, and it was read as evidence of a defect.

**Resolved: the sector is right, the state is a taste partner.** Two probes,
both conclusive.

*Continuum limit at fixed `m/g = 0.1`* (target 0.4326 if this were the pion):

| x | 16 | 36 | 64 | 100 | extrapolated |
|---|---|---|---|---|---|
| M/g | 1.25452 | 1.09327 | 0.99256 | 0.92936 | **0.722** (`1/sqrt(x)`) / 0.896 (`1/x`) |

*Mass dependence at `x = 36`*: local exponents **0.425, 0.308, 0.205** as `m/g`
falls 0.4 → 0.05 — decreasing, so `M/g` tends to a constant near 0.85 rather
than to zero.

A state with non-zero mass in the chiral limit is not the pion, and no
extrapolation in `a` or `a²` brings it to 0.4326. **The cause is taste
breaking.** Staggered fermions in 1+1d preserve exactly *one* chiral (shift)
symmetry; only its partner pion is protected and goes massless at `m = 0`. The
other taste partners acquire a lattice mass that survives the chiral limit and
vanishes only as `a -> 0`. A generic flavour rotation of the vacuum — which is
what `(Sz0, Sz1) = (+1,-1)` is — lands on a taste partner, and that is precisely
what these numbers describe.

So the measurement machinery works and the sector bookkeeping is correct; the
**operator** is wrong. Recovering `m^{2/3}` requires the state associated with
the exact staggered shift symmetry, not an arbitrary flavour sector. That is a
different construction, not a parameter change.

### Status: the two-flavour exponent is unresolved, and needs a fine lattice

Three routes tried, each failing for a different and now-understood reason:

1. **Orthogonality-constrained DMRG** → the flavour singlet (eta), `M/g -> 0.970`.
   Orthogonality does not select a flavour sector.
2. **Flavour-sector targeting** (`charges="independent"`, `Q = [2,-2]`) → a
   staggered **taste partner**, which keeps a lattice mass in the chiral limit.
   The sector is right; the operator is not.
3. **Isovector correlator** → reads *higher* still (3.77 at `m/g = 0.4`, `x = 16`,
   against 1.94 from route 2), because the fit window cannot reach the
   asymptotic regime at a coarse lattice.

Route 3's failure is quantitative and worth keeping, because it says where the
calculation has to live:

| x | M·a | decay length | signal at r=4 | at r=24 |
|---|---|---|---|---|
| 16 | 0.485 | 2.1 sites | 1.4e-1 | 9e-6 |
| 36 | 0.182 | 5.5 sites | 4.8e-1 | 1e-2 |
| 100 | 0.060 | 16.7 sites | 7.9e-1 | 2e-1 |

A correlator needs `M·a << 1` to have an exponential to fit at all — a **fine**
lattice, which is the expensive corner, and the opposite of where these runs are
cheap. At `x = 16` the window runs from marginal signal straight into numerical
noise.

**Route 3 tested on fine lattices, and it fails too.** `A = (M/g)/(m/g)^{2/3}`
should be flat at ~2.008:

| x | m/g=0.2 | 0.1 | 0.05 |
|---|---|---|---|
| 64 | 6.58 | 8.41 | 11.51 |
| 144 | 5.92 | 7.49 | 10.16 |

It climbs at both spacings, so the state stays massive in the chiral limit. And
the correlator reads *heavier* than the variational search (1.61 against 1.09 at
comparable spacing) — which is the diagnostic that matters: **a correlator cannot
be heavier than a state a variational method already found unless its operator
barely couples to that state.** So resolution was not the whole story, and the
`(-1)^n` isovector density is not the protected pion's interpolating operator
here.

**So the remaining work is not a parameter tweak.** It needs the taste-exact
interpolating operator (the one tied to the exact staggered shift symmetry, not a
generic flavour rotation) *and* a lattice fine enough for its correlator to
decay over many sites *and* the volume to hold `M·Lg >> 1` at the same time.
Those pull in opposite directions on cost, which is why this is a piece of work
rather than a fix.

**What is established regardless:** the Hamiltonian is verified four ways, the
continuum mass check passes at +0.06%, and the machinery reaches the states it
targets. The open item is narrow and well-posed: *which operator creates the
protected pion*, answered before more compute is spent.

**Consequences for Phase 2, which are larger than for Phase 1.** The SU(2)
measurement wants the vector meson, sitting near a diquark of similar mass, in a
theory with more sectors than this one. Nothing about "take the first excited
state" is safe there. The fix is to target the sector deliberately — a quantum
number the MPS conserves, or a penalty term — and the two-flavour Schwinger
model is now the right place to develop and validate that, because the answer is
known (`2/3`, amplitude `~2.008`) and the wrong answer is unmistakable.

## Phases

### Phase 1 — the Schwinger model, U(1). Validates everything, and is itself a result.

Elimination is exact and trivial, the model is a long-range spin chain, and TeNPy
handles it with a stock `CouplingMPOModel`. Everything downstream is exercised
here before any non-abelian work starts: the MPO, the DMRG driver, the condensate
measurement, the mass-shift determination, the exponent fit, the ED cross-check.

Validation targets, all published:
* `m = 0`: boson mass `M = e/√π` — the sharpest single check of the whole chain.
* the additive mass shift `m_r/g = m_lat/g + 1/(8√x)` (arXiv:2206.05308),
  reproduced by locating the point where the discrete chiral/shift symmetry is
  restored — the same technique that must then be *derived* for SU(2).
* **two-flavour Schwinger: `M ~ m^{2/3}`.** This is not just a test. It is the
  framework's one clean equal-time pass — the identical logic
  `M ~ m^{1/(2−Δ)}` with `Δ[ψ̄ψ] = 1/2` — and reproducing it establishes that
  our pipeline can *see* an anomalous exponent at all before we ask whether SU(2)
  has one. If Phase 1 cannot recover 2/3, nothing in Phase 2 is interpretable.

### Phase 2 — SU(2), `N_f = 1`

* Local Hilbert space: staggered site holds 0, 1 (two colours), or 2 quarks →
  4 states/site, times the incoming irrep label.
* Two routes, and the first is preferred:
  * **loop-string-hadron** (**arXiv:1912.06133**, MPS version arXiv:2501.18301) —
    *preferred*: local, manifestly gauge invariant, extends to SU(3), and has a
    published MPS implementation to follow.  Costs a bosonic cutoff to
    extrapolate in;
  * **exact elimination in the irrep basis** (Bañuls / arXiv:2511.00154) — no
    truncation error at all, at the cost of writing the recoupling by hand;
  * **Kogut–Susskind with `j ≤ j_max`** — simpler to write, but Bañuls found
    `j_max ≥ 3/2` necessary (their `j_max = 1` gave `nu = 0.781(93)(65)`, visibly
    corrupted), so it buys an extrapolation we do not need.
* Cross-check every MPO against a ~20-site exact diagonalisation built with
  `scipy.sparse` / `qutip 5.2.3` — both already installed — before trusting DMRG.
* Reproduce Bañuls' `nu = 0.700` **in their own window** before quoting anything
  from a wider one.

### Phase 3 — the measurements

* **`⟨ψ̄ψ⟩(m)` first, not the gap.** A ground-state one-point function: free once
  DMRG converges, no excited state, no correlator fit, no plateau identification,
  and no penalty term to separate the near-degenerate diquark (which Bañuls
  needed). The discrimination is also starker — `m^{1/(2N−1)} = m^{1/3}` at N=2
  against flat, a factor 2.15 over `m/g ∈ [0.01, 0.1]`.
* **Then the gap**, and check the two against Feynman–Hellmann
  (`⟨ψ̄ψ⟩ = −∂ε_vac/∂m`), which locks the two exponents together and is a free
  internal consistency test.
* **Then the meson/baryon ratio**, which is what tests this repo's own result.
* Over a genuine decade, `m/g ∈ [0.01, 0.2]`, with the mass shift subtracted.
* **U(N_c) vs SU(N_c) at `N_f = 1` is a near-free null test**: the anomaly gives
  U(N_c) a nonzero `⟨ψ̄ψ⟩(0)` while SU(N_c) must vanish as `m^{1/(2N_c−1)}`.
  Berruto et al. (hep-lat/0201010) report the U side. If our machinery cannot
  separate those two, it cannot address the question.

## Codebase

New, self-contained, no changes to the DLCQ path:

```
equaltime/
  __init__.py
  schwinger.py     Phase 1 model: U(1) staggered, gauge eliminated, long-range
                   Coulomb MPO.  CouplingMPOModel subclass.
  su2.py           Phase 2 model: irrep-basis SU(2), recoupling coefficients.
  massshift.py     locate the additive shift by restoration of the discrete
                   chiral/shift symmetry; the SU(2) coefficient is new work
  observables.py   condensate, gap, meson/baryon ratio, Feynman-Hellmann check
  ed.py            small-lattice exact diagonalisation (scipy.sparse/qutip) used
                   only to validate the MPOs
tools/
  equaltime_scan.py    driver; writes CSVs in the house format (see
                       tools/large_k_sweep.py) into data/equaltime/
tests/
  test_equaltime.py    M = e/sqrt(pi) at m=0; MPO vs ED on 8-12 sites;
                       two-flavour m^(2/3); Banuls nu=0.700 in their window
```

Reuse from this repo: nothing computational (the DLCQ path shares no physics with
a Hamiltonian lattice), but **the conventions matter** —
`docs/table1-units.md` and `dlcq/units.py` define how anything gets compared back,
and exponents are convention-free while prefactors are not. The lattice `g` and
this repo's `lambda ↔ m/g` normalisation differ; fix that before comparing any
coefficient, only after comparing exponents.

## Cost and the risk that could sink it

* **Dependency:** `pip install physics-tenpy` (★493, actively maintained). Adds a
  dependency the repo does not currently have; `requirements.txt` pins only
  numpy/scipy/numba/matplotlib/h5py/pillow/pytest, so this would be the first.
  Alternative if that is unwanted: `quimb`, or hand-rolled DMRG on `scipy.sparse`
  (feasible in 1d but a week of avoidable work).
* **Compute:** trivial. The system is gapped for `m > 0` with `S ~ (c/6)ln ξ`,
  `c = 1`, so `χ` grows polylogarithmically. N = 200–400 sites, χ = 100–400 is
  seconds to minutes per ground state; a full grid of ~20 masses × 4 spacings ×
  4 volumes × 4 bond dimensions is a few CPU-days, embarrassingly parallel.
  Lenore (32 cores, 125 GB) absorbs it overnight.
* **Effort:** Phase 1 ~3–5 days including validation. Phase 2 ~1 week, dominated
  by the recoupling algebra. Phase 3 ~2 days plus run time.
* **The mass shift is no longer the killer risk** — it is zero for SU(N), and
  the bias would run the wrong way regardless. Keep the numerical symmetry-
  restoration check as a cheap null test, not as a mitigation.
* **What replaces it as the leading risk is the boundary Gauss law** in the
  finite-volume validation ladder (see Corrections §1). It is two lines to fix
  and it is *invisible in the ground state*, so every `E_0` check passes while
  the gap is wrong by ~30%.
* **And the truncation `n_l,max`**, which is now the only extrapolation the
  method carries.
* **Secondary risk:** contamination of the meson by the near-degenerate
  diquark/baryon, which Bañuls handled with an explicit penalty term. Budget for
  needing one.
