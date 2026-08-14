# Open threads

Three, in the order they are worth taking. Each is written to be started cold.

Everything here was established by measurement; the supporting numbers are in
`performance.md`, `table1-units.md` and `baryon-higher-fock.md`.

---

## 1. Sparse / Lanczos solver, behind a flag

**Do not replace the dense path.** It stays the default and stays the
reference: it is what every bit-exactness test compares against.

Flag: `run_python(solver="dense"|"sparse")`, `PythonProvider(solver=)`,
`--solver` on the CLI. **The cache tag must include it** — unlike `backend=`,
which is excluded because the two backends produce bit-identical matrices. A
Lanczos path will not, so caching across it would silently mix results.

### Why this and not more K

Measured, not assumed: extending the Richardson window from 2K = 25–35 to
25–37 gave a **median error reduction of 1.00×**. Only 13 of 30 entries moved
at all and every one was at strong coupling, where the error was already
negligible. One extra K shortens the extrapolation by 5%; the weak-coupling
errors are set by how many correction terms to keep, which one point does not
resolve. **There is no incremental path — 2K ≳ 48 is the only lever**, and even
that should be measured rather than assumed.

### Stages

**A — swap the eigensolver only.** Keep the dense build; replace `eigh` with
`scipy.sparse.linalg.eigsh` for the lowest ~30. This alone removes the O(n³)
eigensolve, **48% of a 2K = 35 run**, and is small enough to validate
completely. Do it first: it establishes the flag and the validation pattern
before anything structural moves. Check against dense at 2K = 21–29 to 1e-10,
and confirm `test_reader`'s three published eigenvalues still pass at rtol
1e-12.

Nothing downstream needs more than this: the figures use at most 30 eigenvalues
and 11 eigenvectors, and `DLCQResult` already tolerates fewer eigenvectors than
eigenvalues (that is what `has_eigenvector`/`require_eigenvector` exist for).

**B — never materialize the norm.** Build it per block from
`dlcq.read_python.config_block_labels`, which is already exact and already used
by weeding. Blocks are tiny (max 6 → 11 over 2K = 21 → 31). Store `Z` as blocks:
`orthonormalize_blockwise` already *proves* it is block-diagonal and then
assembles it dense.

**C — sparse H.** Do **not** scan all pairs: that is 1.15e10 pair tests at
n = 150k. Enumerate candidates from configuration keys plus the exact
**|ΔL| ≤ 2** rule — measured over every non-zero off-diagonal at 2K = 21 and 23,
not one couples states differing by more than one qq̄ pair, and it is still
unused. Density falls 26% → 8.5% across 2K = 21 → 31 while `n` grows 1.58× per
+2 and `nnz/row` only 1.25×, projecting to ~0.7% at 2K = 48: **~2 GB sparse
against 185 GB dense.**

**D — matrix-free triple product.** `y = Zᵀ(H(Zx))`. `Zᵀ H0 Z` need never be
formed at all: `H0 = D·N` with `D` constant on each norm block (off-diagonal of
`N⁻¹H0` measured at 7e-15), so it is diagonal by construction in a blockwise
basis.

### Constraints that must survive

- The sparse path supports **`assembly="exact"` + blockwise `Z` only**. The
  Fortran-compatible diagonal-only assembly is inherently dense and basis
  dependent; it stays on the dense path for small-K cross-validation.
- **Matrix element values must not move.** `tests/test_kernels.py` has six
  `array_equal` assertions and `tests/test_fortran_python.py` asserts the norm
  is bit-identical (`== 0.0`) to the 1990 Fortran. Storage, blocking and
  traversal order across independent rows are unconstrained; **reassociating
  the arithmetic is not.**
- `structure_function` needs `norm @ c` as a matvec, and the weight is
  `c * (Nc)`, **not** `c²` — `c²` silently breaks every sum rule.

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
| **thesis figs. 19, 20, 22** | meson four-quark contribution **as a function of K** | the only published K-scan of a higher-Fock structure function anywhere in this work. It measures the K-dependence of the conversion instead of inferring it. Even-channel, so it cannot settle the odd-channel question above — but it is the only direct probe there is |
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

## A method note worth keeping

An apparent disagreement between two extrapolations can be **entirely** an
artifact of fitting different K windows. Refitting on the other side's window is
cheap and should be the first check whenever a published extrapolation fails to
reproduce. It is what dissolved most of the baryon N=4 discrepancy at
m/g = 1.6 (see `table1-units.md`).
