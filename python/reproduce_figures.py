#!/usr/bin/env python3
"""Deprecated shim -- the figure pipeline now lives in ``dlcq.figures``.

The previous version of this file computed everything from the Python solver
and could not read Fortran output at all, so it could not be used to compare
the two codes.  It also carried several approximations that quietly changed the
physics:

* ``get_lightest_mass`` took whichever K gave the *highest positive* lightest
  eigenvalue.  The paper uses Richardson extrapolation (Eq. 27) over 2K ~ 16-24;
  the heuristic cannot reproduce Table I.  Replaced by
  ``dlcq.observables.richardson_extrapolate``, which fits Eq. (27) with the
  non-analytic exponent from Eq. (26) and reports the last retained term as the
  paper's own error estimate.

* ``if lam >= 0.998: return 0.0`` hard-coded the chiral-limit masslessness.
  That is one of the paper's *exact results* and belongs in the test suite, not
  in the solver; see ``tests/test_paper.py``.

* It read ``eigenvalues[0]`` directly.  At N=3, B=0, 2K=24 that is a spurious
  M^2 = 0 caused by an array-bounds overflow in ``qcdf.f``, so the "lightest
  mass" came back as 0 instead of M/g = 3.617.  Use
  ``dlcq.observables.physical_indices``; see ``docs/fortran-color-overflow.md``.

* Bare ``except:`` blocks turned solver failures into silently missing curves.

Use instead::

    python -m dlcq.figures --source fortran --fig 5 6
    python -m dlcq.figures --source python  --fig 5 6 --ncpus 8
    python -m dlcq.figures --source python  --table1

The original implementation is preserved in git history (commit e426ea4).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

if __name__ == "__main__":
    print(__doc__)
    print("Forwarding to dlcq.figures ...\n")
    from dlcq.figures import main

    main()
