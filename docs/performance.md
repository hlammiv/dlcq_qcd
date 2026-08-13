# Performance

Two questions, both measured on a 14-physical-core machine, N = 3, B = 1
baryon, m/g = 1.6.

## Is Python faster than the Fortran?

**Yes — 7–11× as the Fortran is shipped, 3–4× against an `-O2` rebuild.**

| 2K | Fortran `-O0` (as shipped) | Fortran `-O2` | Python |
|---|---|---|---|
| 15 | 4.1 s | 1.7 s | **0.51 s** |
| 19 | 21.1 s | 7.6 s | **1.47 s** |
| 21 | 43.8 s | 16.2 s | **3.6 s** |
| 23 | 112.8 s | 37.8 s | **8.6 s** |

`fortran/makefile` passes only `-mcmodel=medium`, no optimization flag, so the
shipped binary is effectively `-O0`. That is deliberate and should stay:
`docs/basis-dependence.md` records that rebuilding at `-O2` **changes qcdf.f's
own answer** — 190 retained states instead of 189, ground state 10.390084
instead of 10.390380 — because the diagonal-only assembly step is basis
dependent. The `-O2` column above is for a fair speed comparison only, not a
recommendation.

Those Python numbers predate the threading work below; the current figures are
in "Threading" further down. The port wins for two reasons that have nothing to
do with the language: the hot path is numba-compiled, and it is parallel where
`qcdf.f` is strictly serial.

## Where the time goes, and what was done about it

Profiling a 2K = 23 run, before any of this work:

| phase | time | share |
|---|---|---|
| weeding | 40.8 s | **76%** |
| Hamiltonian build | 8.3 s | 16% |
| norm build | 3.0 s | 6% |
| triple product + `eigh` | 0.9 s | 2% |
| state generation | 0.5 s | 1% |

Two changes followed, both exploiting the same structural fact — **the norm is
block-diagonal in momentum configuration**, verified by explicit colour
enumeration in `tools/colour_norm.py`.

**Weeding per block** (`weed_fortran(blocked=True)`). A state can only be
linearly dependent on states in its own block, so both `WEEDR` and `WEEDR2`
decompose. At 2K = 23 the 897 states form 180 blocks of median size 4, so the
O(n²) work drops ~120×. Measured on the step: **528× at 2K = 21, 890× at 2K = 23.**

**Norm build by configuration group.** The worker already declined to compute a
colour factor for mismatched states, so the waste was not the contraction — it
was the O(n²) × O(L²) scan that discovers the mismatch. Hashing each state's
parton multiset once, in O(nL), and handing each row only its own group's
columns replaces that scan entirely. At 2K = 23 the norm build went 2.75 s →
0.55 s.

**Pool scheduling.** Row *r* computes *n − r* elements, so `pool.map`'s default
contiguous chunking is badly unbalanced. Tasks are now ordered longest-first and
dispatched with `chunksize=1`.

End to end:

| 2K | states | before this work | now |
|---|---|---|---|
| 21 | 502 → 193 | 18 s | **3.6 s** |
| 23 | 897 → 319 | 53 s | **8.6 s** |
| 25 | 1559 → 510 | 283 s | **16.5 s** |
| 27 | 2692 → 818 | did not finish in 5 min | **37 s** |
| 29 | 4529 → 1274 | — | ~180 s |

So 2K = 29 is routine where 2K = 25 had been the practical limit. The paper
extrapolates over 2K = 16–24; that window can now be widened.

## Threading

The process pool plateaued at ~5×, and the reason turned out **not** to be the
one recorded here previously. IPC was measured, not assumed: one Hamiltonian
build at 2K = 23 moves **2.4 MB** across the process boundary in 319 tasks, and
accumulating all of it into the destination matrix costs **2 ms** out of a 4.4 s
build. Pool creation is 10–19 ms. Effectively all the time was inside
`pool.map`, doing arithmetic.

The actual reason is simpler: `clfact` and `hamqcd` were never compiled. Only
the micro-kernels (`brack`, `conj_state`, `selfen`) carried `@njit`. Profiling
one row batch:

| routine | tottime | share |
|---|---|---|
| `clfact` | 4.62 s | **79%** |
| `hamqcd` | 0.50 s | 7% |
| `numpy.zeros` (called 140831×) | 0.65 s | 9% |
| `gprbig_array` (called 15409×) | 0.11 s | 2% |

So the work was interpreter overhead, which more processes can only divide, not
remove — and dividing it across a 6 P-core + 8 E-core machine is what capped
the useful thread count at ~8.

`python/qcdf_kernels.py` compiles both routines with **`nogil=True`**, which
fixes both halves at once: the serial path stops being interpreted, and because
the kernels hold no lock, a plain `ThreadPoolExecutor` gives real parallelism
with no fork, no pickling, and each row written straight into the destination
matrix. Three things had to change to make it compile, all behaviour-preserving:

- **The epsilon table is hoisted.** `clfact` rebuilt the SU(N) permutation table
  via `gprbig_array` on *every call* — 15409 times in the batch above — for a
  table that depends only on `N`. It is now built once per run by
  `epsb_table(N)`.
- **Scratch buffers are per-thread.** The module-global `ClfactBuffers` was safe
  only because `fork` gave each process a copy.
- **The vertex table is arrays, not lambdas.** `hamqcd`'s `V4` closures became
  integer tables: the `BRACK` arguments as coefficient vectors over `(k₀..k₃)`,
  the flavour assignment as a 0/1 selector over `(f₁, f₂)`.

**Both backends produce bit-identical matrices.** That is asserted, not assumed:
`tests/test_kernels.py` builds the norm and the Hamiltonian both ways for eight
configurations — mesons and baryons, SU(2)/SU(3)/SU(4), one and two flavours —
and requires `array_equal`, not a tolerance. Exact equality is the right bar
here, because `docs/basis-dependence.md` records that the assembly step is basis
dependent at the 1e-4 level, so anything looser would be indistinguishable from
noise the Fortran already has against itself. `PythonProvider` therefore leaves
the backend out of its cache tag.

Hamiltonian build alone, 14 threads:

| 2K | states | process pool | threads | |
|---|---|---|---|---|
| 21 | 193 | 2.01 s | **0.09 s** | 22× |
| 23 | 319 | 4.95 s | **0.12 s** | 41× |
| 25 | 510 | 11.44 s | **0.14 s** | 82× |
| 27 | 818 | 24.17 s | **0.20 s** | 121× |
| 29 | 1274 | 50.39 s | **0.31 s** | 163× |

End to end (`run_python`, N = 3, B = 1, m/g = 1.6, medians of 3):

| 2K | before | after | |
|---|---|---|---|
| 21 | 2.43 s | **0.66 s** | 3.7× |
| 23 | 6.09 s | **1.07 s** | 5.7× |
| 25 | 13.27 s | **1.51 s** | 8.8× |
| 27 | 27.81 s | **1.99 s** | 14× |
| 29 | 59.92 s | **4.87 s** | 12× |

Thread scaling on the build itself, 2K = 29: 1.34 s → 0.42 (4) → 0.26 (8) →
0.21 (14) → 0.18 (20), i.e. **7.3×**. Past 14 the gain is small, as expected
when 6 of the 20 logical CPUs are hyperthreads and 8 are E-cores. `ncpus` above
14 is not harmful, just not useful.

Select with `--backend`, `backend=` on `run_python`/`PythonProvider`, or
`QCDF_BACKEND`. The process backend is kept because it runs the *interpreted
reference* routines, which is what makes the equality test meaningful.

## State generation

The same treatment, applied to the next bottleneck. `qcdf_states.py` compiles
`brmsgn`, `prmx`, `prmm` and their helpers with `nogil=True`; `qcdf_opt.
qcdsta_fast` drives them. `qcdf.py` keeps its versions as the interpreted
reference, and `tests/test_states.py` asserts the two produce
**element-identical** `mstate`/`mstinf` — not merely the same set of states,
because state *order* is load-bearing: `WEEDR`/`WEEDR2` drop linearly dependent
states by index, so a permuted basis would silently change which states survive.

The profile said the work was not where it looked. At 2K = 29, `brmsgn`'s own
bytecode was 73% of generation — not a callee, but the `istate` fill loops doing
~44 boxed `int64` item-assignments per candidate, ~9.2M per run. Underneath that
sat three compounding wastes:

- **A 46:1 candidate funnel.** 209,152 candidates were fully constructed for
  4,529 kept; 96.3% died on `oddchk`. But `oddchk` is a *momentum-parity* test,
  and momenta come from the `kpx`/`kpm` row picked by the quotient
  `(ist[l,j]-1)//dim` while flavour comes from the remainder — so parity depends
  only on the row, and `oddchk` is a conjunction over the baryon/antibaryon/meson
  groups. Filtering each group's rows once is equivalent and collapses the funnel.
- **1.9 MB allocated per call, 1,129 times.** `gprbig_array` allocated a fresh
  `(ISTMAX_BG, MXNP)` table each call and wrote ~67 rows on average. Worse,
  absent species got the same 1.9 MB and never read it — ~2.7 GB of zeroed pages
  per run. Buffers are now caller-owned.
- **`stachk` hashed Python dicts** — 364,389 lookups per run. Its two dicts are
  keyed injectively on disjoint domains, so the function reduces to "does any
  value occur more than `indx` times"; with `nptcl <= 25` a pairwise count is
  equivalent and allocation-free.

Measured, N = 3, B = 1 (interpreted vs compiled, both element-identical):

| 2K | states | interpreted | compiled | |
|---|---|---|---|---|
| 29 | 4,529 | 1.4 s | **0.004 s** | 313× |
| 33 | 12,471 | 9.4 s | 0.08 s | 125× |
| 35 | 20,353 | 17.4 s | **0.02 s** | 736× |

**No cap needs hand-tuning again.** `kpx` and the state arrays now double on
demand, so `LKPXMX`/`NMSTMX` bound only the interpreted reference. Reach of the
generator alone, which is no longer a constraint anywhere in range:

| 2K | 37 | 39 | 41 | 43 | | 2K (meson) | 42 | 46 | **48** |
|---|---|---|---|---|---|---|---|---|---|
| states | 32,816 | 52,519 | 83,167 | 130,794 | | | 56,375 | 127,500 | **190,171** |
| time | 0.05 s | 0.13 s | 0.23 s | 0.47 s | | | 0.83 s | 3.4 s | **5.8 s** |

Threading it over momentum sectors was planned and then **measured to be
unnecessary**: generation is now 0.2% of a run. The sectors are independent and
the split is easy to add if the sparse work below ever makes it matter again.

## What is left

**The profile has inverted twice.** The Hamiltonian build was 88% of a run and
is now 7%; state generation was 41% and is now 0.2%. What dominates at 2K = 35:

| phase | time | share |
|---|---|---|
| dense `eigh` (all eigenpairs) | 13.98 s | **48%** |
| NUZ (`eigh` of the norm + scaling) | 8.68 s | **30%** |
| weeding | 2.99 s | 10% |
| Hamiltonian build | 2.05 s | 7% |
| norm build | 1.36 s | 5% |
| state generation | 0.05 s | 0.2% |

**It is now dense linear algebra, and the wall is memory, not time.** 2K = 37 is
OOM-killed on a 15 GB machine. The norm is built at the *pre*-weeding size, so:

| case | n_pre | dense n² |
|---|---|---|
| baryon 2K = 35 | 20,353 | 3.3 GB |
| baryon 2K = 37 | 32,816 | 8.6 GB |
| baryon 2K = 39 | 52,519 | 22.1 GB |
| baryon 2K = 41 | 83,167 | 55.3 GB |

`_symmetrize` used to multiply that by four (`np.triu(a) + np.triu(a,1).T`
touches three extra n² arrays); it now mirrors in place. The immediate remaining
offender is `weed_fortran`, which densifies `|N| > 1e-9` into a CSR matrix
(`dlcq/read_python.py:158`) — `np.abs(A)` alone is another full n² — to recover
block labels that `qcdf_opt._config_keys` already computes in O(n·L) without
touching the matrix. That is what kills 2K = 37 today.

**The sparse rewrite, now the only thing between here and 2K = 48.** Density
falls steadily, which is what makes the target reachable:

| 2K | n_post | nnz/row | density |
|---|---|---|---|
| 21 | 193 | 50.3 | 26.05% |
| 25 | 510 | 85.0 | 16.67% |
| 29 | 1,274 | 136.0 | 10.68% |
| 31 | 1,983 | 169.0 | 8.52% |

`n` grows 1.58× per +2 in K; `nnz/row` only 1.25×. Extrapolated to 2K = 48:
~150k states at ~0.7% density, so **sparse H is ~2 GB where dense is 185 GB**.
The norm is better still — block-diagonal, maximum block size 6 → 11.

The pieces: block-diagonal norm built from `_config_keys` labels and never
materialized; `Z` stored as blocks (it is provably block-diagonal under
`orthonormalize_blockwise`, yet assembled dense today); sparse H with candidate
pairs enumerated from configuration keys plus the exact **|ΔL| ≤ 2** rule
(measured over every non-zero off-diagonal at 2K = 21 and 23 — not one couples
states differing by more than one qq̄ pair — and still unused); and Lanczos for
the lowest ~30 eigenvalues instead of a full `eigh`. Nothing downstream needs
more: the figures use at most 30 eigenvalues and 11 eigenvectors, and
`DLCQResult` already tolerates fewer eigenvectors than eigenvalues.

`Z^T H0 Z` need never be formed at all: `H0 = D·N` with `D` constant on each
norm block (off-diagonal of `N⁻¹H0` measured at 7e-15), so in a blockwise basis
it is diagonal by construction.
