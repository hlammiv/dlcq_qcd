# An array-bounds overflow corrupts high Fock sectors in `qcdf.f`

**Summary.** A colour matrix element between two `L`-parton states needs
`2L + 4` colour-index slots. `qcdf.f` dimensions that array with 25. Any run
reaching **11 or more partons in a single Fock state** therefore overflows it,
silently, and the corrupted states appear as **non-positive `M²` eigenvalues at
the bottom of the spectrum** — exactly where Figs. 7, 8 and Table I read off the
lightest mass.

The Python port is not affected: it uses `MXLNG = 2*MXP + 4 = 54`, and its
author had already documented the correct size.

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
  subject to this bug.

`qcdf.f` is deliberately left unpatched so it continues to reproduce the
historical output byte-for-byte.
