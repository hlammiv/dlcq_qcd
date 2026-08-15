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
2. **No additive mass renormalisation appears to be applied.** For the Schwinger
   model, Hamiltonian staggered has `m_r/g = m_lat/g + 1/(8√x)`
   (arXiv:2206.05308, confirmed with MPS in arXiv:2303.11016). **No published
   SU(N) analogue was found.** The shift is `O(g²a)` — the same order as the `m`
   whose exponent is being fitted — so leaving it out biases `nu` systematically,
   and plausibly explains 0.700 against 2/3. Deriving or measuring it is the
   single most valuable new piece of work available here.

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

**Therefore, before writing anything: ask.** Requesting code from Bañuls/Kühn or
Fujikura/Hidaka is ordinary practice and would collapse Phase 2 — the week of
recoupling algebra — into running an existing implementation at smaller `m/g`.
Their papers answer a different question with the same machinery, so there is no
competitive reason for them to decline. That is the cheapest path by a wide
margin and should be exhausted first.

If that fails, the phases below are the from-scratch route.

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
* **The killer risk is the mass shift.** It is `O(g²a)`, the same order as the
  quantity whose exponent is being fitted, no SU(N) value is published, and
  getting it wrong biases `nu` in exactly the direction that would fake or hide
  the effect. Mitigations: work at several `ag` and extrapolate; locate the shift
  numerically via symmetry restoration; and lead with the condensate, whose
  leading `m^{1/3}` is less sensitive to an `O(g²a)` shift in `m` than a log-log
  slope fit is.
* **Secondary risk:** contamination of the meson by the near-degenerate
  diquark/baryon, which Bañuls handled with an explicit penalty term. Budget for
  needing one.
