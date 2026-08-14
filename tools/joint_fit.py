#!/usr/bin/env python3
"""Extrapolate in K and m/g at once, anchored at the chiral point.

Why a joint fit at all
----------------------
Fitting each coupling separately leaves weak coupling undetermined, and no
amount of computing fixes it.  The convergence rate is the endpoint exponent --
measured, ``p ~ 1 + a``: 1.985 against 1+a=1.845 at m/g=1.6, and 1.061 against
1.084 at m/g=0.1.  So as ``a -> 0`` the series converges as slowly as it can,
the tail bound falls only as ``K^-0.06``, and halving the bracket at m/g=0.1
would need K x 84,000.  A different fitting *basis* cannot help either: the
confluent basis spans the same space and leaves M0 identical to 1.1e-14.  The
only remaining lever is a change of model.

The anchor, and why it is allowed
---------------------------------
Eq. (16): at m/g = 0 the lightest state is exactly zero **for any N, B and K**.
That is stronger than a statement about the continuum limit -- every finite-K
curve passes through the origin too, so the K-corrections must vanish there as
well.

Non-analyticity was the obvious objection, and it was measured rather than
waved away.  At fixed K, ``d ln M^2 / d ln(m/g)`` runs 1.400, 1.813, 1.950,
1.987, **1.9968** across m/g = 0.4 down to 0.0125, and agrees to four digits at
2K=49 and 2K=71.  Integer power, no logarithm, no fractional exponent, no K
dependence in the exponent.  ``a/(m/g) -> 0.8463`` and is flat below m/g=0.05.
Both analytic at the origin, so the anchor is legitimate.

The model
---------
``M^2 / (m/g)^2`` is fitted -- factoring out the measured power is what makes
the target O(1) instead of spanning 0.011 to 10.7, and an unweighted fit over
that range is dominated by its largest values::

    M^2(mg, K) / mg^2 = sum_j g_j mg^j + sum_k sum_i b_ki mg^i Kp^-e_k(mg)

with ``e_k`` the Eq. (27) exponents ``1, 1+a, 2, 2+a`` at each row's own
``a(mg)``.  ``j`` starts at 0 (so the limit starts at ``mg^2`` overall) and the
correction terms start at ``i=0`` likewise.

Conditioning
------------
Columns are equilibrated to unit norm and the solve is an SVD with explicit
rank truncation.  Without that the 16-parameter variant reached condition
number 1.6e18 -- numerically singular, interpolating rather than fitting, so
its 1.7e-05 residual meant nothing.  The rank actually used is reported.

Validation
----------
Two hold-outs, because they fail differently:

* **highest K** -- the ordinary extrapolation check, present for every fit.
* **an entire coupling** -- only a joint fit can even attempt this, and it is
  the one that tests whether tying the couplings together is doing real work or
  just adding parameters.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dlcq.units import endpoint_exponent


def design(mg, Kp, N, n_mass, n_corr, n_corr_mass):
    """Columns of the joint model for ``M^2/mg^2``, and a label for each."""
    mg = np.asarray(mg, float)
    Kp = np.asarray(Kp, float)
    a = np.array([endpoint_exponent(float(m), N) for m in mg])

    cols, labels = [], []
    for j in range(n_mass):
        cols.append(mg ** j)
        labels.append(f"limit mg^{j}")
    exps = [(np.ones_like(a), "K^-1"), (a + 1.0, "K^-(1+a)"),
            (np.full_like(a, 2.0), "K^-2"), (a + 2.0, "K^-(2+a)")]
    for ex, name in exps[:n_corr]:
        for i in range(n_corr_mass):
            cols.append(mg ** i * Kp ** (-ex))
            labels.append(f"{name} mg^{i}")
    return np.vstack(cols).T, labels, n_mass


def solve(A, y, rcond=1e-10):
    """Equilibrated least squares with explicit rank truncation.

    Scaling each column to unit norm removes the part of the ill-conditioning
    that is only a choice of units -- ``mg^3`` and ``K^-2`` differ by orders of
    magnitude in size before they differ in direction -- so what the rank cut
    then removes is genuine degeneracy rather than scale.
    """
    scale = np.linalg.norm(A, axis=0)
    scale[scale == 0] = 1.0
    As = A / scale
    U, s, Vt = np.linalg.svd(As, full_matrices=False)
    keep = s > rcond * s[0]
    rank = int(keep.sum())
    coef_s = Vt[keep].T @ ((U[:, keep].T @ y) / s[keep])
    return coef_s / scale, rank, (s[0] / s[keep][-1] if rank else np.inf)


def limit_at(mg_q, coef, n_mass):
    """Continuum M^2 at a coupling: the limit series, times the mg^2 taken out."""
    return sum(coef[j] * mg_q ** j for j in range(n_mass)) * mg_q ** 2


def fit_all(mg, Kp, y, N, n_mass, n_corr, n_corr_mass, rcond=1e-10):
    A, labels, nm = design(mg, Kp, N, n_mass, n_corr, n_corr_mass)
    coef, rank, cond = solve(A, y / mg ** 2, rcond)
    pred = (A @ coef) * mg ** 2
    return coef, rank, cond, labels, pred


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--N", type=int, default=3)
    ap.add_argument("--mg-max", type=float, default=0.45,
                    help="the (m/g)^2 form is asymptotic; large m/g breaks it")
    ap.add_argument("--rcond", type=float, default=1e-10)
    args = ap.parse_args(argv)

    rows = list(csv.DictReader(open(args.csv)))
    mg = np.array([float(r["mg"]) for r in rows])
    Kp = np.array([float(r["K_code"]) for r in rows]) / 2.0
    y = np.array([float(r["msq"]) for r in rows])
    sel = mg <= args.mg_max
    mg, Kp, y = mg[sel], Kp[sel], y[sel]
    couplings = sorted(set(mg), reverse=True)
    print(f"{len(y)} points, m/g <= {args.mg_max}: "
          f"{[round(c, 4) for c in couplings]}\n")

    # Several model sizes, so the spread between them is reported as the
    # uncertainty rather than one of them being presented as the answer.
    variants = [(nm, nc, ncm) for nm in (2, 3, 4)
                for nc in (2, 4) for ncm in (2, 3)]

    print("%-14s %5s %5s %10s %11s | %s" % (
        "model", "par", "rank", "cond", "max rel res",
        "  ".join("%9.4g" % c for c in couplings)))
    keep = []
    for nm, nc, ncm in variants:
        coef, rank, cond, labels, pred = fit_all(mg, Kp, y, args.N, nm, nc, ncm,
                                                 args.rcond)
        rel = np.abs(pred - y) / np.abs(y)
        lims = [limit_at(c, coef, nm) for c in couplings]
        flag = "" if rank == len(coef) else "  (rank deficient)"
        print("%-14s %5d %5d %10.2e %11.2e | %s%s" % (
            f"{nm}/{nc}/{ncm}", len(coef), rank, cond, rel.max(),
            "  ".join("%9.5f" % v for v in lims), flag))
        if rel.max() < 5e-3:
            keep.append(lims)

    if keep:
        arr = np.array(keep)
        print("\nspread over %d acceptable variants (max rel residual < 5e-3):"
              % len(arr))
        print("%10s %11s %11s %11s %9s" % ("m/g", "median", "min", "max",
                                           "spread"))
        for i, c in enumerate(couplings):
            col = arr[:, i]
            print("%10.4g %11.5f %11.5f %11.5f %8.1f%%" % (
                c, np.median(col), col.min(), col.max(),
                100 * (col.max() - col.min()) / abs(np.median(col))))

    # ── hold out the highest K ────────────────────────────────────────────
    Kmax = Kp.max()
    tr = Kp < Kmax
    print(f"\nheld-out K = {Kmax:g} (fit without it, predict it):")
    print("%10s %13s %13s %10s" % ("m/g", "predicted", "actual", "rel err"))
    coef, rank, cond, _, _ = fit_all(mg[tr], Kp[tr], y[tr], args.N, 3, 4, 2,
                                     args.rcond)
    A_h, _, _ = design(mg[~tr], Kp[~tr], args.N, 3, 4, 2)
    pred = (A_h @ coef) * mg[~tr] ** 2
    for m, p, t in zip(mg[~tr], pred, y[~tr]):
        print("%10.4g %13.6f %13.6f %10.2e" % (m, p, t, abs(p - t) / abs(t)))

    # ── hold out a whole coupling ─────────────────────────────────────────
    print("\nheld-out coupling (fit without it entirely, then predict its "
          "continuum value):")
    print("%10s %15s %15s %10s" % ("m/g", "from held-out", "from all", "diff"))
    coef_all, *_ = fit_all(mg, Kp, y, args.N, 3, 4, 2, args.rcond)
    for c in couplings:
        tr = mg != c
        if len(set(mg[tr])) < 3:
            continue
        coef_h, *_ = fit_all(mg[tr], Kp[tr], y[tr], args.N, 3, 4, 2, args.rcond)
        a_, b_ = limit_at(c, coef_h, 3), limit_at(c, coef_all, 3)
        print("%10.4g %15.6f %15.6f %10.2e" % (c, a_, b_, abs(a_ - b_)))


if __name__ == "__main__":
    main()
