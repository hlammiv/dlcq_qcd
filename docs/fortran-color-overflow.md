# An array-bounds overflow corrupts high Fock sectors in `qcdf.f`

**Summary.** A colour matrix element between two `L`-parton states needs
`2L + 4` colour-index slots. `qcdf.f` dimensions that array with 25. Any run
reaching **11 or more partons in a single Fock state** therefore overflows it,
silently, and the corrupted states appear as **non-positive `M²` eigenvalues at
the bottom of the spectrum** — exactly where Figs. 7, 8 and Table I read off the
lightest mass.

The Python port is not affected **by this particular overflow**: it uses
`MXLNG = 2*MXP + 4 = 54`, and its author had already documented the correct
size. It has a *different* colour-array defect of its own, at the term count
rather than the slot count — see [The Python port's own overflow](#the-python-ports-own-overflow)
below. Earlier revisions of this document said flatly that the Python solver was
unaffected; that was wrong.

## The arithmetic

```fortran
      INTEGER IDEL0(12552,25),IDELT(12552,25)      ! qcdf.f:4994
      MXLNG=25                                     ! qcdf.f:5018, 5672
```

```python
MXTRM = 12552; MXLNG = 54; MXP = 25
# MXLNG = 2*MXP + 4 (max right state + max operators + max left state)
#                                                  qcdf_opt.py:24-25
```

A four-point vertex contracts the bra's `L` partons, the ket's `L` partons and
up to 4 operators, so the index row must hold `2L + 4` entries. With 25 slots
the limit is `L = 10`; at `L = 11` the write runs past the row.

Note that `qcdf.f` *does* guard the other dimension — `NTRMS.GT.MXTRMS` prints
`NTRMS EXCEEDED MXTRMS IN FN CLFACT` and stops. That guard never fires here.
There is no guard on the second dimension, so the overflow is silent.

## Measured symptoms

| run | max partons | slots needed | result |
|---|---|---|---|
| 2K=10, B=0 | 6 | 16 | clean |
| 2K=21, B=1 | 9 | 22 | clean, lowest `M²` = 10.3904 |
| 2K=24, B=0 | 12 | 28 | **one decoupled 12-parton state with an identically zero Hamiltonian row, giving a spurious `M² = 0` as the ground state** |
| 2K=25, B=1 | 11 | 26 | **four negative `M²`, every one dominated by an 11-parton state** |

The correlation is exact: runs that fit in 25 slots have no non-positive
eigenvalues; runs that do not have them, and the offending eigenvectors live in
precisely the sector that overflows.

### The 2K=24 case in detail

The Fortran retains a single 12-parton state — three quarks and three
antiquarks at `k = 1`, and the same at `k = 3`, so its norm diagonal is
`(3!)⁴ = 1296`. Its entire Hamiltonian row comes back zero:

```
HNU0[816]        = 0.000000000000e+00
max |HNU[816,:]| = 0.000000000000e+00
```

A state of massive quarks must have free energy `Σ 2m²/kᵢ > 0`, so that is
wrong. Running the same configuration through the Python solver gives the same
state a free energy of `HAM0[i,i] = 6912`, i.e. `6912/1296 = 5.333` in the
orthonormal basis, which would place it near `M² ≈ 57` — a heavy state, not a
massless one.

## Why the paper is probably unaffected

`input_parameters.pdf` documents the particle-number truncation `LPN`:

> "the number of states can grow rapidly; this allows you to limit by hand the
> number of particles in higher-Fock states. For mesons, you might start with 2
> or 4, then allow for more later until it gets out of hand."

Any `LPN ≤ 10` keeps the run inside the 25-slot array. The repo's own
`input_small.txt` uses `LPN = 4`. The authors very likely ran their production
cases truncated, which both bounds the basis size and incidentally avoids this
bug. Our reproductions use `LPN = 0` (no truncation) and therefore hit it.

The physics impact on the *low-lying* spectrum is small either way — the
11- and 12-parton sectors carry ~10⁻⁸ of the ground-state wave function — but
the artifacts land at the bottom of the eigenvalue list, so anything that reads
`eigenvalues[0]` is corrupted outright.

## How this repository handles it

- `dlcq.observables.spurious_zero_modes` flags eigenvalues `≤ 0` at `m/g > 0`.
- `dlcq.observables.physical_indices` removes them; every figure that speaks of
  "the first three states" or "the lightest state" indexes through it.
- `dlcq.observables.fortran_overflow_risk` reports whether a given run exceeded
  the 25-slot limit, and `dlcq.read_fortran` warns on load.
- The Python solver (`assembly="exact"`) is the scientific reference and is not
  subject to *this* bug, but see the section below for one of its own.

`qcdf.f` is deliberately left unpatched so it continues to reproduce the
historical output byte-for-byte.

## The Python port's own overflow

Found 2026-08-13. Different array, same shape of failure, and it went unnoticed
for a specific reason worth recording.

`qcdf.f` runs out of **index slots** (`IDELT(12552,25)`, needing `2L+4`). The
Python port sized that correctly at 54, but inherited the *other* dimension
verbatim: `MXTRM = 12552`, the number of colour contraction **terms**. Two
guards enforced it, and neither said anything:

| site | what it did |
|---|---|
| `qcdf_kernels.py` `clfact_nb` | `return 0.0` — the whole colour factor becomes zero |
| `qcdf_kernels.py` `hamqcd_nb` | `break` — leaves a **partial sum** |

with identical copies in the interpreted `qcdf_opt.clfact`/`hamqcd`. The partial
sum is the dangerous one: it returns a plausible non-zero number rather than an
obvious zero.

**Why the tests could not catch it.** `tests/test_kernels.py` establishes
correctness by demanding the threaded and interpreted builds agree bit-for-bit.
Both carry the same `MXTRM`, so they agreed bit-for-bit on the wrong answer.
Backend equality is blind to any defect both backends share.

**What it looked like**, at N=4, B=1, 2K=20, against
`tools/colour_norm.norm_bruteforce` — which sums colours explicitly and is
therefore the true Gram matrix:

```
25 wrong entries, every one between L=12 states
worst was a DIAGONAL: 150048 against a true 331776
assembled norm had an eigenvalue of -1.8e5 (a Gram matrix cannot)
```

**Scope.** The term count grows with both `N` and `L`, so higher `N` overflows
at lower parton number. Measured at `MXTRM = 12552`:

| configuration | maxL | outcome |
|---|---|---|
| N=3, B=1, LPN=0, 2K ≤ 37 | 13 | ok |
| N=3, B=1, LPN=0, 2K ≥ 39 | 15 | **overflows** |
| N=3, B=0, LPN=0, 2K=40 | 14 | ok |
| N=4, B=1, LPN=0, 2K=20 (norm) | 12 | **overflows** |
| N=4, B=1, LPN=0, 2K=16 (**Hamiltonian**) | 10 | **overflows** |

Two things that are easy to get wrong. It is not an N ≥ 4 problem — N=3
overflows too, just later. And parton count alone does not predict it: at N=4,
2K=16 the *norm* is fine while the *Hamiltonian* overruns, because the
Hamiltonian's four-point vertices generate far more terms per element.

**No published result is affected.** Table I and Figs. 7–8 run under
`sweep_lpn`, whose cap exists precisely to keep `2L+4` inside the Fortran's
array and which therefore also holds the term count far below `MXTRM` (L ≤ 5 at
every K); Figs. 4–6 are N=3 at 2K ≤ 29, well inside the limit.

**Fix, in two parts.**

1. Both guards now raise `OverflowError` naming `MXTRM`. A wrong answer returned
   quietly is worse than no answer.
2. `MXTRM` was raised from the inherited Fortran value of 12552 to **200000**.
   That is not a workaround — it makes the sums *exact*. Rechecked against
   `tools/colour_norm.norm_bruteforce` afterwards:

   | configuration | max \|diagrammatic − brute force\| |
   |---|---|
   | N=4, B=1, 2K=20 | **0** (was 181728) |
   | N=4, B=1, 2K=16 | **0** |
   | N=3, B=1, 2K=21 | 0 |
   | N=3, B=1, 2K=25 | 0 |

   and the N=4, 2K=20 norm is positive semidefinite again (−3.4e-11, against
   −1.8e5 before).

   **Size it against the Hamiltonian, not the norm.** The norm at N=3, LPN=0
   needs only 25000 to reach 2K=41, and sizing on that alone is a trap: the
   Hamiltonian's four-point vertices generate far more terms per element, and
   at 2K=39 it still overran at 50000. Measured requirements:

   | build | 2K=39 | 2K=41 |
   |---|---|---|
   | norm | fits at 25000 | fits at 25000 |
   | **Hamiltonian** | needs > 50000 | needs > 50000 |

   200000 covers both (2K=39 in 24.3 s at 2.0 GB peak, 2K=41 in 20.3 s at
   4.1 GB, the peak being the dense Hamiltonian rather than the scratch).

   Cost is much smaller than the reservation suggests: `np.zeros` gives lazy
   zero pages and only rows `0..ntrms` are ever touched, so residency tracks the
   real term count. Measured after the increase, 2K=29 runs in **1.52 s**
   against the **4.87 s** in `performance.md` from before it. And since the
   guard now fails loudly, under-sizing is safe — it can only refuse to run,
   never corrupt.

**The invariant that does catch it**, now in `tests/test_colour_overflow.py`:
the norm is a Gram matrix, so it is positive semidefinite. Checked per
configuration block — the norm is exactly block-diagonal in parton content, so
that is equivalent to checking the whole matrix and costs 46 s at 2K=41, where
the dense norm would be 51 GB.
