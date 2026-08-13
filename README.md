# DLCQ QCD in 1+1 dimensions

Validating two solvers — Kent Hornbostel's 1993 Fortran 77 code and a Python
port of it — against **Hornbostel, Brodsky & Pauli, Phys. Rev. D 41, 3814
(1990)**, and freezing that validation into an automated regression suite.

The article PDF is not redistributed here; see [CITATION.md](CITATION.md).

## Layout

```
fortran/   qcdf.f (6397 lines, F77) -- the historical code, deliberately unpatched
python/    qcdf.py / qcdf_opt.py    -- the Python port
dlcq/      the shared pipeline both solvers feed
refs/      Table I transcribed from the paper
docs/      findings that change how results must be read (docs/README.md indexes them)
tests/     the regression suite
```

## The idea

Both solvers produce the same `DLCQResult`, and the figure code reads only that.
A figure built from Fortran output and one built from Python go through
identical arithmetic, so any difference in the plot is a difference in physics.

```
fortran/qcdf.f ──> qcdf.out + qcdf.ham ──┐
                                         ├──> DLCQResult ──> observables ──> figures
python/qcdf_opt.py ──> arrays ───────────┘                                    └──> tests
```

| module | role |
|---|---|
| `dlcq/units.py` | λ ↔ m/g, `K_code = 2·K_paper`, `M²_code` → `M/g`, Eq. 26 exponent |
| `dlcq/dataset.py` | `DLCQResult` + HDF5 — the data contract |
| `dlcq/read_fortran.py` | parses `qcdf.out` / `qcdf.ham` |
| `dlcq/read_python.py` | runs the Python solver |
| `dlcq/providers.py` | `PythonProvider` / `FortranProvider`, both cached |
| `dlcq/observables.py` | structure functions, sum rules, Richardson (Eq. 27) |
| `dlcq/figures.py` | Figs. 1–8 and Table I, from `DLCQResult` only |

## Quick start

```bash
pip install -r requirements.txt
make -C fortran                                   # needs gfortran

# run a case (qcdf.f hard-codes its output filenames, so runs are isolated)
bash fortran/run_case.sh runs/K21 3 1 1 0.3325 -1.0 0 21

# figures, from either solver
python -m dlcq.figures --source fortran --fig 5 6
python -m dlcq.figures --source python  --fig 5 6 --ncpus 8

pytest                    # Tier 0-1, under a second
pytest -m slow            # the expensive sweeps
```

## Conventions worth reading once

| quantity | definition |
|---|---|
| `K_code` | **2×** the paper's `K = (L/2π)P⁺`. The paper's "2K = 24" is `K_code = 24` |
| λ | `1/(1 + πm²/g²)^½`; `m/g = [(1−λ²)/(πλ²)]^½` |
| eigenvalue | `M²_code = K_code·w/2`, in units of `m² + g²/π` |
| `M/g` | `sqrt(M²_code/(πλ²))` |
| structure fn | `q(x) = K_paper⟨φ\|b†_k b_k\|φ⟩`, `x = k/K_code`, `dx = 2/K_code` |
| momenta | odd integers (antiperiodic boundary conditions) |

Two traps. **λ must be passed literally**: the historical runs use `0.3325`
while `mg_to_lambda(1.6)` is `0.33254949`, and that 1.5×10⁻⁵ swamps any tight
tolerance. And this paper's natural unit is `g²N/2π`, not 't Hooft's `g²N/π` —
hence the `(2π/N)^½` rescaling on Fig. 8(b).

## Does it reproduce the paper?

Every figure and table, from **both** solvers. `docs/figure-validation.md` has
the per-curve accounting; this is the summary.

| | quantity | agreement |
|---|---|---|
| **Table I** | 30 extrapolated masses, M² units | 0.02–4% at m/g ≥ 0.8 |
| **Fig. 2** | eigenvalue spectra vs coupling | 1.0–1.8% of the axis |
| **Fig. 3** | valence structure functions | 1.4–1.9% |
| **Fig. 4** | higher-Fock, plotted alone | 0.0–3.9% (meson, 5q) |
| **Fig. 5** | meson wavefunctions, 2K=24 | 0.3–2.5% valence, 0.9–2.0% higher-Fock |
| **Fig. 6** | baryon wavefunctions, 2K=21 | 0.4–2.0% valence |
| **Fig. 6(d)** | two-baryon, 2K=24 | 0.3% valence, 0.8% higher-Fock |
| **Fig. 7** | extrapolated masses | see Table I — it plots the same numbers |
| **Fig. 8** | large-N and lattice comparison | 't Hooft eigenvalues to 4 digits |

The sharpest single comparison is Table I at m/g = 1.6, where the paper quotes
four significant figures: mesons 4.314 / 4.618 / 4.845 for N = 2, 3, 4 against
our 4.313 / 4.615 / 4.844.

**Two curves are not reproduced**, both baryon higher-Fock: Figs. 6(a)–(c)'s
five-quark series at 2K = 21, and the two-extra-pair sector (Fig. 4(c) /
thesis Fig. 18(b), 17.5% high). `docs/baryon-higher-fock.md` records what was
eliminated: the colour sums are exact against brute-force enumeration, the norm
equals the true Gram matrix, the eigenvector is determined to 10⁻¹², and the
distribution is invariant under precision loss down to float16, under the
weeding threshold over eight orders of magnitude, and under corrupting every
high-parton matrix element. The same sector at 2K = 15 is reproduced to 0.4%.

### Six ways a published axis can be misread

Worth its own line, because it caused every false alarm in this project. Five
panels have a frame top **above** their highest printed label, and Fig. 5(b)
uses its neighbour's axis entirely. Read naively these give apparent
discrepancies of 16%, 20%, 44% and 233% — each of which looks exactly like a
physics result. `tools/pin_axes.py` fits the scale from the label positions
instead, and checks itself against the value at the frame bottom, which must
come out zero. And Fig. 6(d)'s K is never stated in its caption; it is 24, not
the 22 we first assumed, which alone accounted for that panel's disagreement.

## Speed

The Python port is **7–11× faster than `qcdf.f` as shipped** and 3–4× faster
than an `-O2` rebuild of it — numba-compiled and parallel, against a serial
unoptimized F77. (The shipped binary carries no optimization flag deliberately:
rebuilding at `-O2` changes its own answer, see
[docs/basis-dependence.md](docs/basis-dependence.md).) Details and the remaining
opportunities are in [docs/performance.md](docs/performance.md).

## Reach

The paper extrapolates over 2K = 16–24. That window is set by cost, and the cost
was dominated by an avoidable step: **weeding** — discarding linearly dependent
basis states — was 76% of a run's wall time, against 16% for the Hamiltonian
build and 6% for the norm.

The norm is block-diagonal in momentum configuration (verified by explicit
colour enumeration, `tools/colour_norm.py`), so a state can only be linearly
dependent on states in its own block and both weeding stages decompose. At
2K = 23 the 897 states form 180 blocks of median size 4, so weeding per block is
~120× less work. Measured, on the weeding step alone: **528× at 2K = 21, 890× at
2K = 23**.

| 2K | states | before | after |
|---|---|---|---|
| 21 | 502 → 193 | 18 s | **3.6 s** |
| 23 | 897 → 319 | 53 s | **8.6 s** |
| 25 | 1559 → 510 | 283 s | **16.5 s** |
| 27 | 2692 → 818 | did not finish in 5 min | **37 s** |
| 29 | 4529 → 1274 | — | ~**180 s** |

So 2K = 29 is now routine where 2K = 25 had been the practical limit, and the
extrapolation window can be widened past the paper's. The two passes are not
bit-identical — see `weed_fortran`'s docstring — but the ground state agrees to
9×10⁻⁸ and the low-lying spectrum to 10⁻⁶, a thousandfold inside this
algorithm's own 10⁻⁴ reproducibility floor.

## One piece had to be rebuilt, not ported

`qcdf.f` emits eigenvectors, Fock content and the basis-change matrix — but not
structure functions. Those were produced in the original work by a **separate
program, `wf`/`wfbig`, which does not survive**; the source refers to it only in
comments. Every structure function in the paper passed through it.

`dlcq/observables.py:structure_function` is a reconstruction of it, and it is the
one part of the pipeline a solver-vs-solver comparison cannot validate — the
Python port matches the Fortran's matrix elements to `1e-14`, so an error in the
conversion would be shared by both. It is checked five other ways instead: the
thesis's own definition (Sec. 2.4), a proof that the `c_s (Nc)_s` weighting is
the exact expectation value given the norm's block structure, a brute-force
colour enumeration confirming that norm exactly, both sum rules to machine
precision across four runs and both solvers, and agreement with the published
curves at 0.0–2.0% (mesons), 0.4–3.9% (baryon five-quark) and 0.8% (two-baryon).

See [docs/baryon-higher-fock.md](docs/baryon-higher-fock.md) §1.3.

## Port validation

Verified to machine precision, against the Fortran's own basis and its own `Z`:

| quantity | agreement |
|---|---|
| norm matrix (189×189) | **exact**, max diff `0.0` |
| interacting `HNU = ZᵀHZ` | `1.1e-14` |
| free `HNU0` | `8.9e-15` |
| `ZᵀNZ = I` | `2.7e-14` |
| momentum sum rule `∫x[q+q̄]dx = 1` | `1.4e-14` |
| number sum rule `∫[q−q̄]dx = N·B` | `4.0e-14` |

So the state generator, the diagrammatic colour-contraction engine, all seven
four-point vertices, both self-energies and the structure-function machinery are
confirmed correct. The rebuilt Fortran also reproduces the preserved 1990-era
`qcdf.out`/`qcdf.ham` **byte-for-byte** at 2K = 21 and 2K = 25.

Reproduced from the paper's own text: at 2K=24 the 11th meson state has
`M/g` ratio **1.989** against the stated "twice that of the first", and its
qq̄qq̄ component peaks at `x = 5/24 = 0.208`, the nearest odd-momentum grid
point to the stated `x = 1/4`.

## Multiple flavours

`NF > 1` was added to `qcdf.f` in May 1990 and **no figure in the paper uses
it** — it shipped with no coverage. It works, and is the natural direction for
new physics: unequal quark masses give kaon and strange-baryon analogues the
paper never computes. `tests/test_flavour.py` now covers it via both sum rules,
flavour symmetry at degenerate masses, two variational inequalities, and
cross-solver agreement to `3e-7`.

Two things to know first, both in [docs/flavour.md](docs/flavour.md): `qcdf.f`
**cannot reach flavour non-singlet channels at all** (it overwrites `iflv(1)`
with `N*B`, so a kaon-like channel silently generates zero states), and its
colour-array overflow arrives at `2K = 12` rather than `2K = 24`.

## Two defects found in `qcdf.f`

Both are documented in full under [docs/](docs/README.md). Neither changes the paper's physics
conclusions, but both change how its output must be read.

**1. The spectrum is basis-dependent** —
[docs/basis-dependence.md](docs/basis-dependence.md).
`qcdf.f` adds the free Hamiltonian to the *diagonal only* of the
orthonormal-basis matrix. `Zᵀ H₀ Z` is not diagonal, and the norm matrix has
just 26 distinct eigenvalues among 189 states, so the result depends on which
eigenvectors the diagonalization returns. Recompiling the **unmodified** source
at `-O2` instead of `-O0` changes its own answer: 190 retained states instead of
189, ground state 10.390084 instead of 10.390380. ~10⁻⁴ is the intrinsic
reproducibility floor, for anyone. The end-to-end test derives its tolerance
from this rather than hard-coding one.

**2. A silent array-bounds overflow** —
[docs/fortran-color-overflow.md](docs/fortran-color-overflow.md).
A matrix element between two `L`-parton states needs `2L + 4` colour-index
slots; `IDELT` is dimensioned with 25. Runs reaching `L ≥ 11` are corrupted, and
the damage lands as non-positive `M²` at the **bottom** of the spectrum — where
Figs. 7, 8 and Table I read the lightest state. At 2K=24, B=0 a naive
`eigenvalues[0]` returns `M/g = 0` instead of 3.617.
`dlcq.observables.physical_indices` removes these; the reader warns on load.

Neither is patched in `qcdf.f`, which stays byte-faithful to the code that
produced the paper. The corrected physics lives in the Python path, where
`assembly="exact"` is the basis-independent reference.

## Reproduction status by figure

| Fig | status | notes |
|---|---|---|
| 1 | schematic | interaction vertices; nothing to validate |
| 2 | **both** | λ sweep at 2K = 10/13/22, 21 couplings each |
| 3 | **agrees, 2.6% / 1.9%** | K recovered as 2K = 14 / 15 (never stated) |
| 4 | **(b) agrees, 0.07–3.9%** | via the thesis twin Fig 18(a); q and q̄ both |
| 5 | **agrees, 1.2–1.4%** | incl. the ×10³ and ×10² rescaled series, to 1–2% |
| 6 | **(a) agrees, 1.1%** | higher-Fock verified via thesis Fig 18(a) |
| 7 | **both** | Richardson over 2K = 16–24, N = 2,3,4, both sectors |
| 8 | **both** | incl. the 't Hooft large-N curve; matches Hamer's lattice data to 0.1–2.2% |
| Table I | **both** | reproduced once read in M² units — see `docs/table1-units.md` |

Every figure is generated from **both** solvers; `figures/` carries a
`_fortran` and a `_python` copy of each. Producing a figure is not the same as
agreeing with the published one — see
[docs/figure-validation.md](docs/figure-validation.md) for the per-curve
accounting, including the series the paper rescales by powers of ten. Fortran
and Python agree with each other to within 0.1 percentage point on every curve
compared.

Cross-checked over the whole 174-configuration sweep, the two codes agree on
the lightest M² to a **median of 3.9×10⁻⁶** (max 4.0×10⁻³, which is the
basis-dependence floor). After Richardson extrapolation the meson Table I
columns agree exactly; the baryon columns differ by up to 2.5%, since
extrapolation amplifies the per-K spread — and there the Python value is the
one closer to the paper.

`dlcq/thooft.py` solves Eq. (24) for Fig. 8's large-N curve. It reproduces the
classic chiral-limit eigenvalues 0, 5.88, 14.1 in the paper's own coupling
normalization, and agrees with the finite-N DLCQ sweep to 0.07% at m/g = 1.6 —
an independent equation, discretization and code path.

The paper never states `K` for Figs. 3 and 4; it was **recovered from the plots
themselves** as 2K = 14 (meson) and 15 (baryon). Momenta are odd, so markers lie
on `x = k/K`, and fitting that lattice reproduces the stated 2K = 24 and 21 of
Figs. 5 and 6 before being trusted where K is unknown. Every meson panel
returned an even K and every baryon panel an odd one — the parity each sector
requires, which the fit knows nothing about. See `docs/inferred-K.md`.

## Table I is in M² units, not M/g

Table I tabulates `M²/(m² + g²/π)` — the y-axis of **Fig. 2**, and the code's
raw eigenvalue — despite its column headers reading `M_mes/g` and `M_bar/g`.
Read that way our sweep reproduces it; read as `M/g` nothing matches.

At the best-converged coupling, m/g = 1.6, four of the five columns agree to
better than **0.1%**, and 17 of the 30 non-trivial entries fall within the
paper's own quoted uncertainty. The `M/g` reading gives ratios scattered over
0.55–4.8.

Three things corroborate it. Figs. 7 and 8(a) *are* `M/g` — their axes run 0–4
and 0–8, and our values sit inside while Table I's would not. **Hamer's
independent lattice data** matches our `M/g` to 0.1% at m/g = 1.6. And Table I's
apparent "unphysical saturation" is required in M² units, since numerator and
denominator both grow like m². See `docs/table1-units.md`.
