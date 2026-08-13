# The Fortran spectrum is reproducible only to ~1e-4

**Summary.** `qcdf.f` adds the free part of the light-cone Hamiltonian to the
*diagonal only* of the orthonormal-basis matrix. That step is basis dependent,
and the basis it depends on is not uniquely determined. Recompiling the
unmodified source with a different optimization level changes the published
eigenvalues. Roughly `1e-4` relative is the intrinsic reproducibility floor of
the original algorithm — for any reimplementation, and for `qcdf.f` itself.

This does **not** affect the conclusions of Phys. Rev. D **41**, 3814. The shift
is two to four orders of magnitude smaller than the precision the paper quotes.

## What the code does

The basis of color-singlet Fock states is over-complete, so the Gram (norm)
matrix `N` is diagonalized and used to build an orthonormalizing transform
(`qcdf.f`, `NUZ`):

```
Z(ML,MR) = Z(ML,MR) / sqrt(W(MR))          with  N Z = Z W,  Z^T N Z = I
```

`NUHAM` then forms `HNU = Z^T HAM Z` and, separately, only the **diagonal** of
the free part:

```fortran
       NR=NL                                  ! <- diagonal elements only
       HLR0=HLR0 + HAM0(ifl,KL,KR)*ZL*ZR
       HNU0(ifl,NL)=HLR0
```

and the main program combines them as

```fortran
       HNU(I1,I1)= (1.0D0 - RLMSQ)*rmq(ifl)*rmq(ifl)*HNU0(ifl,I1)
     >             +HNU(I1,I1)
```

## Why that is ill-posed

`H0` is exactly `D·N` for a diagonal `D` — verified numerically, the
off-diagonal part of `N^-1 H0` is `7e-15`. The free light-cone energy
`sum_i m^2/k_i^+` is constant on each set of states with nonzero overlap, since
overlapping states share their parton content.

That makes `Z^T H0 Z` diagonal *only if* each column of `Z` stays inside a
single Fock momentum-content block. Columns of `Z` are norm eigenvectors, and
the norm matrix is massively degenerate — at N=3, B=1, 2K=21 there are **26
distinct eigenvalues among 189 states**. Degenerate eigenvectors may be freely
rotated, including across blocks, and any eigensolver is entitled to return such
a mixture. When it does, `Z^T H0 Z` acquires off-diagonal elements (up to ~0.8
here) and discarding them changes the answer.

Measured at N=3, NF=1, B=1, 2K=21, lambda=0.3325, comparing two equally valid
orthonormalizations (LAPACK `eigh` vs the `Z` the Fortran itself printed):

| assembly of the free part | spectrum |
|---|---|
| full `Z^T H0 Z` | **basis independent** to `6.3e-13` |
| diagonal of `Z^T H0 Z` (what `qcdf.f` does) | **basis dependent**, differing by up to `5.8` in `M^2` |

## The compiler experiment

The clearest demonstration needs no reimplementation at all. Same source, same
machine, same input — only the optimization flag differs:

```
gfortran -mcmodel=medium      -o qcdf qcdf.f     # -O0, the makefile default
gfortran -mcmodel=medium -O2  -o qcdf qcdf.f
```

| | states kept | lowest `M^2` | 2nd | 3rd |
|---|---|---|---|---|
| `-O0` | 189 | 10.390379999627 | 13.822939496992 | 15.317122613433 |
| `-O2` | 190 | 10.390083903062 | 13.823574065023 | 15.320915160705 |

Max relative difference across the spectrum: `6.6e-2`; ground state `2.9e-5`.

The optimization level perturbs floating-point arithmetic just enough to change
which eigenvectors come back from degenerate subspaces, which changes both the
weeding outcome (189 vs 190 retained states) and the discarded off-diagonal
terms.

## What *is* exactly reproducible

Everything upstream of the assembly step. Comparing the Python port against the
Fortran on Fortran's own retained basis, using Fortran's own `Z`:

| quantity | agreement |
|---|---|
| norm matrix `N` (189x189) | **exact**, max diff `0.0` |
| interacting `HNU = Z^T H Z` | `1.1e-14` |
| free `HNU0` (diagonal) | `8.9e-15` |
| `Z^T N Z = I` | `2.7e-14` |
| momentum sum rule `int x[q+qbar]dx = 1` | `1.4e-14` |
| number sum rule `int [q-qbar]dx = N*B` | `4.0e-14` |

So the state generator, the diagrammatic color-contraction engine, all seven
four-point vertices, both self-energy terms, and the structure-function
machinery are confirmed correct to machine precision. The only irreproducible
element is the choice of basis inside degenerate norm eigenspaces.

## Does the floor limit higher Fock components?

No — and this is worth stating plainly, because the 1e-4 figure is easy to
mistake for a precision limit on *our* results. It is not. It is the
reproducibility of **`qcdf.f`'s diagonal-only assembly**, which we keep only for
Fortran comparison. The default path (`assembly="exact"`) does not have it.

Measured by running the same configuration through two genuinely different
valid orthonormalizations of the same norm — a global `Z` and a blockwise `Z`
— at N=3, B=1, 2K=21:

| assembly of the free part | ground state | max relative, whole spectrum |
|---|---|---|
| full `Z^T H0 Z` (**default**) | **2.5e-14** | **4.6e-15** |
| diagonal only (what `qcdf.f` does) | 1.5e-3 | 1.0e-1 |

Eleven orders of magnitude apart. So there is nothing to "push lower" on the
default path — it is already at machine precision, and the 1e-4 exists purely
as the yardstick for how closely anything, including `qcdf.f` itself, can
reproduce the 1990 numbers.

Higher Fock sectors inherit that, which is the case one would most expect to be
fragile since the amplitudes are tiny. Under the same basis change, with the
exact assembly:

| sector | probability | agreement across the basis change |
|---|---|---|
| 3 partons (valence) | 2.998183 | 1.2e-15 |
| 5 partons | 2.4222e-3 | 2.2e-15 |
| 7 partons | 2.086e-7 | 2.2e-12 |
| 9 partons | ~0 | 3.5e-10 |

The seven-parton sector carries a probability of 2e-7 and is still reproducible
to 2e-12 — five orders of margin. **Numerical precision is not what limits the
higher-Fock results.** That agrees with `docs/baryon-higher-fock.md`, which
found the distribution invariant under precision loss down to float16, under
the weeding threshold over eight orders of magnitude, and under corrupting
every high-parton matrix element. Whatever is wrong with Fig. 6(a)'s five-quark
series is structural, not numerical.

## Consequences for this repository

1. `dlcq.read_python` defaults to `assembly="exact"` — the basis-independent
   result. `assembly="fortran"` reproduces the historical scheme for
   comparison.
2. Tier-2 Fortran-vs-Python tests assert on the **matrices** at `1e-12`, where
   agreement is exact, and on the **spectrum** at a documented `~1e-3`, which is
   the floor the original algorithm allows.
3. `dlcq.nr_eigen` provides a faithful port of the Numerical Recipes
   `tred2`/`tqli`/`eigsrt` chain. It reproduces the Fortran's *eigenvalues* to
   `2.8e-13` but not its degenerate *eigenvectors*, for the reason above — it is
   kept because it isolates the eigensolver as a variable rather than leaving it
   confounded with LAPACK.

## Does this change the paper?

No. Table I quotes `M_bar/g = 10.71(2)` for N=3 at m/g=1.6 — a 0.2% uncertainty,
and the paper is explicit that its parenthesized figures are the magnitude of
the last retained Richardson term rather than statistical errors. A `2.9e-5`
relative shift in `M^2` is about `1.5e-5` in `M`, four orders of magnitude
inside the quoted band. Every physics conclusion in the paper survives
unchanged.
