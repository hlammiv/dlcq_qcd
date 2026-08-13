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

The Python port wins for two reasons that have nothing to do with the language:
`qcdf_opt.py` is numba-compiled, and it is parallel where `qcdf.f` is strictly
serial.

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

## What is left

**Parallel efficiency is the weakest point.** At 2K = 21 the scaling was 3.4× on
8 processes (43%) and *negative* past that — 16 processes were slower than 8.
Longest-first dispatch helps, but each task still returns arrays of length
`numsta`, so O(n²) data crosses the process boundary. Threads with `nogil`
numba kernels, or accumulating into shared memory, would avoid that. On a
14-core machine there is roughly a 3× gain still on the table.

**The Hamiltonian build is now the dominant cost** (~80% of a run). It is
genuinely denser than the norm — 21–26% non-zero — so there is no equivalent
structural shortcut. The one selection rule that does hold exactly is
**|ΔL| ≤ 2**: measured over every non-zero off-diagonal at 2K = 21 and 23, not
one couples states differing by more than one qq̄ pair. The worker's existing
four-mismatch filter already captures this, and skipping those pairs before the
matching loop would save about 16%.

**The dramatic option, not attempted.** Every eigenvalue is currently computed,
by dense `eigh`, when the figures and Table I need only the lowest few. A
matrix-free formulation — applying H to a vector directly and using Lanczos —
would replace O(n²) stored elements and O(n³) diagonalization with O(nnz) per
matvec. That is a real rewrite and would need validating against the present
code state by state, but it is what would make 2K ≫ 30 reachable.
