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

## What is left

**The profile has inverted.** The Hamiltonian build was 88% of a 2K = 29 run;
it is now 6%. What dominates now:

| phase | 2K = 29 | share |
|---|---|---|
| state generation (`qcdsta`) | 2.02 s | **41%** |
| dense `eigh` | 1.79 s | **37%** |
| Hamiltonian build | 0.31 s | 6% |
| norm build | 0.26 s | 5% |
| weeding | 0.21 s | 4% |
| NUZ | 0.19 s | 4% |

**State generation is now the single largest cost** and is still pure
interpreted Python in `qcdf.py` (`qcdsta`, `prmx`, `gprbig_array`). Compiling it
is the next step, and it is independent of everything above.

It *was* also what blocked going higher — 2K = 31 died in `prmx` with `index
25001 is out of bounds` — but that was a fixed cap, not a time limit, and the
caps are now raised. Requirements measured on the N = 3, B = 1 baryon:

| 2K | partons | KPX perms | states pre → post | run |
|---|---|---|---|---|
| 29 | 11 | 18801 | 4529 → 1274 | 2.9 s |
| 31 | 13 | 32072 | 7569 → 1983 | 7.8 s |
| 33 | 13 | 45013 | 12471 → 3032 | 16 s |
| 35 | 13 | 62290 | 20353 → 4610 | 46 s |

The old caps were `LKPXMX = 25001` and `NMSTMX = 6902`, both inherited from
`qcdf.f`; they now sit at 80000 and 25000, costing ~37 MB of zeroed array per
run instead of ~11 MB. `MXNP = 25` is **not** binding — `lpnsub` caps this
channel at 13 partons — and `MXKMX`/`MAXK_CLF` have ample headroom.
`ISTMAX_SM` is dead code.

Raising `NMSTMX` also required fixing a latent allocation in `qcdf.py`'s own
`main()`: `z_full` was sized `NMSTMX × NMSTMX` regardless of the problem, which
was 381 MB wasted at the old cap and would have been 5 GB at the new one. `nuz`
and `nuham` only ever touch `[:numsta]`.

Past 2K = 35, raise the caps again; each overflow path now names the constant.

**The dramatic option, not attempted.** Every eigenvalue is computed, by dense
`eigh`, when the figures and Table I need only the lowest few. A matrix-free
formulation — applying H to a vector directly and using Lanczos — would replace
O(n²) stored elements and O(n³) diagonalization with O(nnz) per matvec. That is
a real rewrite and would need validating state by state, but with the build cost
gone it is now the other half of what stands between here and 2K ≫ 30.

**A selection rule still unused.** |ΔL| ≤ 2 holds exactly: measured over every
non-zero off-diagonal at 2K = 21 and 23, not one couples states differing by
more than one qq̄ pair. The four-mismatch filter already captures it; skipping
those pairs before the matching loop would save ~16% of a now-small cost.
