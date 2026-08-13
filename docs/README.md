# Findings

Six documents. Each records something that **changes how the paper's output must
be read** — not implementation notes, and not a log of what was tried.

Start with `figure-validation.md` if you want to know whether the reproduction
worked; with `basis-dependence.md` or `fortran-color-overflow.md` if you intend
to run `qcdf.f` yourself.

## What was reproduced

**[figure-validation.md](figure-validation.md)** — the per-curve accounting for
every figure and Table I, from both solvers, with the measured deviations. Also
the method: how a published axis is pinned from its label positions, and the six
ways doing that wrong produced apparent physics results of 16–233%.

**[baryon-higher-fock.md](baryon-higher-fock.md)** — despite the name, also the
home of §1.3, which documents the one part of the original pipeline that had to
be **rebuilt rather than ported**: the x-space conversion, done in 1990 by a
separate program `wf`/`wfbig` that does not survive. And the two published curves
that are *not* reproduced, and the full list of what was eliminated for them.
Our side is closed from every direction available: the colour sums are exact
against brute-force enumeration, the norm equals the true Gram matrix, the
eigenvector is determined to 10⁻¹², and the distribution is invariant under
precision loss down to float16.

## Facts about the paper

**[table1-units.md](table1-units.md)** — Table I tabulates `M²/(m² + g²/π)`
despite column headers reading `M_mes/g`. Figs. 7 and 8 plot `M/g`. Two
conventions, one paper; comparing across them without converting is a factor of
several.

**[inferred-K.md](inferred-K.md)** — the paper never states K for Figs. 3 and 4,
and never for Fig. 6(d) either. All three can be read off the plots, because
antiperiodic boundary conditions put markers only at odd `k/K`. Fig. 6(d) is
2K = 24, not the 22 that seemed natural, and that alone accounted for its
apparent disagreement.

## Ready for new work

**[performance.md](performance.md)** — is the Python faster than the Fortran?
Yes, decisively: 2K = 29 runs in 4.9 s. Records the measurement that redirected
the parallel work — the process pool's bottleneck was *not* IPC (2.4 MB per
build, 2 ms to accumulate) but the fact that `clfact` and `hamqcd` were never
compiled — and what compiling them `nogil` and threading them bought. Then the
same for state generation, which the profile showed was 96% wasted work: 209,152
candidates built for 4,529 kept, because a momentum-parity test ran per candidate
instead of per momentum row. Also where the time goes now that the build is 7% of
a run and generation 0.2%: dense linear algebra, and a memory wall at 2K = 37.

**[flavour.md](flavour.md)** — NF > 1 was added to `qcdf.f` in 1990 and no figure
in the paper uses it, so it shipped untested. It is the obvious direction for new
work (kaon and strange-baryon analogues), and it works — validated by both sum
rules, flavour symmetry at degenerate masses, two variational inequalities, and
cross-solver agreement to 3 × 10⁻⁷. Two things to know before using it: `qcdf.f`
cannot reach flavour non-singlet channels at all, and the colour-array overflow
arrives at half the K it does for one flavour.

## Defects in `qcdf.f`

Neither changes the paper's physics conclusions; both change how its output must
be read. A third, specific to multiple flavours, is in
[flavour.md](flavour.md). `fortran/qcdf.f` is kept unpatched — corrected behaviour lives in the
Python path.

**[basis-dependence.md](basis-dependence.md)** — `qcdf.f` adds the free
Hamiltonian to the *diagonal only* of the orthonormal-basis matrix. `Zᵀ H₀ Z` is
not diagonal, so the result depends on which eigenvectors the diagonalization
happens to return. Recompiling the unmodified source at `-O2` instead of `-O0`
changes its own answer. About `1e-4` relative is the intrinsic reproducibility
floor of the original algorithm, for anyone.

**[fortran-color-overflow.md](fortran-color-overflow.md)** — a colour matrix
element between two `L`-parton states needs `2L + 4` index slots; `qcdf.f`
dimensions that array with 25. Runs reaching 11+ partons in one Fock state
overflow it silently, and the damage lands as non-positive `M²` at the **bottom**
of the spectrum — exactly where Figs. 7, 8 and Table I read the lightest state.
At 2K = 24, B = 0 a naive `eigenvalues[0]` returns `M/g = 0` instead of 3.617.

## Conventions

`docs/` assumes the conventions table in the top-level [README](../README.md):
`K_code = 2K_paper`, odd momenta, `M²_code = K_code·w/2` in units of `m² + g²/π`.
`dlcq/units.py` is the single source of truth for all of them.
