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

The one genuinely unresolved item. See `baryon-higher-fock.md` for the full
accounting. Eliminated so far: Fock truncation (identical to four decimals at
LPN 0/5/7, i.e. across a 2× change in basis size), numerical precision (the
7-parton sector is basis-independent to 2e-12), state identity, sector
identification, and the tracer.

The strongest remaining lead is the ratio pattern: published/ours is
**1.36, 1.40, 1.43 at k = 1, 3, 9** and **0.30, 0.36, 0.47 at k = 5, 7, 11** —
near-constant within each group, alternating between them. That looks mechanical
rather than physical.

Note also that the same sector is published at 2K = 15 and **we reproduce it
marker for marker**, with a dip depth of 0.63 against 0.60. At 2K = 21 the
published depth is 0.13 where ours is 0.55, and ours varies smoothly (0.55–0.63)
across 2K = 15–25.

---

## A method note worth keeping

An apparent disagreement between two extrapolations can be **entirely** an
artifact of fitting different K windows. Refitting on the other side's window is
cheap and should be the first check whenever a published extrapolation fails to
reproduce. It is what dissolved most of the baryon N=4 discrepancy at
m/g = 1.6 (see `table1-units.md`).
