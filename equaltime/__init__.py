"""Equal-time (Hamiltonian-lattice) calculations, independent of the DLCQ path.

Why this package exists: every finite-`N_c` determination of the 2d chiral
exponent -- this repo's meson and baryon channels, Anand-Fitzpatrick-Katz-Xin,
Kochergin -- is quantised on a null plane with a trivial vacuum, and all give
``alpha -> 1``.  The one equal-time measurement (Banuls et al., PRX 7, 041046)
gives ``alpha = 1.40``.  That shared assumption is exactly what the anomalous
exponent is about.  See ``docs/equal-time-plan.md``.

Nothing here imports from ``dlcq`` or ``python``: it is a different
quantisation of a different discretisation, and the only thing the two share is
the physics they are supposed to agree on.  Units are the one place they must be
reconciled -- see ``docs/table1-units.md`` -- and only after exponents, which
are convention-free, rather than before.
"""
