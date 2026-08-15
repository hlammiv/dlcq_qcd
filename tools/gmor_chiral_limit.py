#!/usr/bin/env python3
"""The chiral limit of the improved Hamiltonian, with uncertainties, at every N present.

Reproduces every GMOR/chiral-exponent number quoted in
``docs/weak-coupling-limit.md``.  The colour range is **discovered from the data
directory**, not hardcoded -- the scan grew from N<=5 to N<=8 mid-analysis and a
hardcoded list silently kept reporting the old answer.

Four blocks, each answering a question the previous one raises:

1. ``M^2(K->inf) / GMOR`` per (N, m/g), with the ensemble bar.  Quoting this
   ratio bare is what made an earlier revision claim "converges to 1 within
   0.5%" when at fixed m/g it does not converge to 1 at all.
2. That same ratio in units of its own bar, so "close to 1" becomes a number.
3. The large-N limit, ``ratio = r_inf + c/N^2``.  At fixed m/g the ratio rises
   with N and is significantly above 1 -- the question is whether it rises
   toward a limit, which needs more than three N to answer.
4. Whether the surviving large-N residual vanishes as m/g -> 0.  It does,
   **linearly**: ``M^2/GMOR = 1 + A(m/g)`` with an intercept consistent with
   zero, the quadratic alternative excluded.  This is the real content of the
   GMOR anchor -- a law with an intercept, not a ratio eyeballed against 1.

Then the chiral exponent itself, in physical units, against the two competing
analytic laws: GMOR's ``alpha = 1`` (exact only as N -> infinity) and the
bosonization result ``alpha = 2N/(2N-1)`` of Date-Frishman-Sonnenschein, which
contains GMOR as its large-N limit.  See docs for why DFS is a hypothesis under
test rather than a validated anchor.

Usage:
    python tools/gmor_chiral_limit.py
    python tools/gmor_chiral_limit.py --hamiltonian standard
    python tools/gmor_chiral_limit.py --data-dir data/gmor_scan --min-n 3
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))


def discover(data_dir: Path, ham: str):
    """Find every ``{ham}_N*.csv`` and return {N: path}, sorted by N."""
    out = {}
    for p in sorted(data_dir.glob(f"{ham}_N*.csv")):
        m = re.search(rf"{re.escape(ham)}_N(\d+)\.csv$", p.name)
        if m:
            out[int(m.group(1))] = p
    return dict(sorted(out.items()))


def read_series(path: Path):
    """-> {(N, mg): {K_code: msq}}.  Skips files with only a header."""
    out = defaultdict(dict)
    with open(path) as fh:
        for r in csv.DictReader(l for l in fh if not l.lstrip().startswith("#")):
            out[(int(r["N"]), float(r["mg"]))][int(r["K_code"])] = float(r["msq"])
    return out


def fit_1k(ks, y, order):
    """Plain 1/K series -- the right basis for improved output (van de Sande Eq. 14)."""
    Kp = np.array(ks, float) / 2.0            # K_code -> K_paper
    y = np.asarray(y, float)
    A = np.vstack([np.ones_like(Kp)] + [Kp ** -(i + 1.0) for i in range(order)]).T
    sc = np.linalg.norm(A, axis=0)
    c, *_ = np.linalg.lstsq(A / sc, y, rcond=None)
    return c / sc


def ensemble(series, min_pts=6):
    """Median and 68% half-width of M^2(K->inf) over orders x contiguous sub-windows.

    The bar and any band drawn from it come from one ensemble, so the quoted
    uncertainty and the picture of it cannot disagree.
    """
    ks = sorted(series)
    y = [series[k] for k in ks]
    n = len(ks)
    vals = []
    for order in (2, 3, 4):
        for i in range(n):
            for j in range(i + max(min_pts, order + 2), n + 1):
                c = fit_1k(ks[i:j], y[i:j], order)
                x = np.array(ks[i:j], float) / 2.0
                A = np.vstack([np.ones_like(x)] +
                              [x ** -(t + 1.0) for t in range(order)]).T
                resid = np.max(np.abs(A @ c - np.array(y[i:j])))
                if resid / max(np.max(np.abs(y[i:j])), 1e-300) < 1e-3 and c[0] > 0:
                    vals.append(c[0])
    if not vals:
        return None
    a = np.array(vals)
    return (float(np.median(a)),
            float(0.5 * (np.percentile(a, 84) - np.percentile(a, 16))),
            len(a))


def gmor(N, mg):
    """van de Sande Eq. (7), M^2 = 2 pi g mu / sqrt(3), in the repo's units."""
    c = (N * N - 1.0) / (2.0 * N * np.pi)
    return 2 * np.pi * np.sqrt(c) * mg / np.sqrt(3.0) / (mg ** 2 + 1.0 / np.pi)


def dfs_alpha(N, NF=1):
    """Bosonization exponent: alpha = 2/(1+P), P = (N^2-1)/(N(N+NF)).

    At NF=1 this is 2N/(2N-1) -- 1.200 at N=3 -- and tends to GMOR's 1 as
    N -> infinity.  Derived from DFS Eq. (4.1) with the Casimir bracket OUTSIDE
    the radical (their own decuplet/octet = 1.41 fixes that) and with m the
    generated scale m ~ m_q^{1/(1+P)}, not the quark mass.
    """
    P = (N * N - 1.0) / (N * (N + NF))
    return 2.0 / (1.0 + P)


def to_physical(msq_code, mg):
    """Table I units -> M^2/g^2.  Matters only at the largest couplings."""
    return msq_code * (mg ** 2 + 1.0 / np.pi)


def wls_intercept(x, y, e):
    """Weighted least squares y = b + a x.  Returns (b, sigma_b, a, chi2/dof)."""
    A = np.vstack([np.ones_like(x), x]).T
    W = np.diag(1.0 / e)
    c, *_ = np.linalg.lstsq(W @ A, W @ y, rcond=None)
    cov = np.linalg.inv(A.T @ np.diag(1.0 / e ** 2) @ A)
    dof = max(len(x) - 2, 1)
    chi2 = float(np.sum(((A @ c - y) / e) ** 2))
    return float(c[0]), float(np.sqrt(cov[0, 0])), float(c[1]), chi2 / dof


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default=str(ROOT / "data" / "gmor_scan"))
    ap.add_argument("--hamiltonian", default="improved",
                    choices=["improved", "standard"])
    ap.add_argument("--min-n", type=int, default=3,
                    help="smallest N included in the large-N fits (N=2 is far "
                         "from large-N and is reported but excluded)")
    args = ap.parse_args(argv)

    data_dir = Path(args.data_dir)
    files = discover(data_dir, args.hamiltonian)
    if not files:
        print(f"no {args.hamiltonian}_N*.csv under {data_dir}")
        return 1

    d = {}
    for N, p in files.items():
        d.update(read_series(p))
    NS = sorted({N for (N, _) in d})
    MGS = sorted({mg for (_, mg) in d}, reverse=True)
    if not NS or not MGS:
        print("files present but contain no rows")
        return 1

    print(f"  {args.hamiltonian}: N = {NS}   m/g = {MGS}")
    print(f"  from {len(files)} files in {data_dir}\n")

    # ── 1. the ratio, with its bar ──────────────────────────────────────────
    R = {}
    print("  M^2(K->inf) / GMOR   +- ensemble 68%")
    print("   m/g   " + "".join(f"{'N=' + str(N):>18}" for N in NS))
    for mg in MGS:
        cells = []
        for N in NS:
            g = d.get((N, mg))
            res = ensemble(g) if g else None
            if not res:
                cells.append(f"{'--':>18}")
                continue
            m0, sd, _ = res
            gm = gmor(N, mg)
            R[(N, mg)] = (m0 / gm, sd / gm, m0, sd)
            cells.append(f"{R[(N, mg)][0]:11.4f}+-{R[(N, mg)][1]:.4f}")
        print(f"  {mg:<7}" + "".join(cells))

    # ── 2. in units of its own bar ──────────────────────────────────────────
    print("\n  deviation from 1, in units of its own bar:")
    for mg in MGS:
        cells = [f"N={N}:{(R[(N, mg)][0] - 1) / R[(N, mg)][1]:+7.1f}s"
                 if (N, mg) in R and R[(N, mg)][1] > 0 else f"N={N}:     --"
                 for N in NS]
        print(f"  {mg:<8}" + " ".join(cells))

    # ── 3. the large-N limit ────────────────────────────────────────────────
    big = {}
    fitns = [N for N in NS if N >= args.min_n]
    print(f"\n  large-N:  ratio = r_inf + c/N^2   over N = {fitns}")
    for mg in MGS:
        ns = [N for N in fitns if (N, mg) in R]
        if len(ns) < 3:
            continue
        x = np.array([1.0 / N ** 2 for N in ns])
        y = np.array([R[(N, mg)][0] for N in ns])
        e = np.array([max(R[(N, mg)][1], 1e-9) for N in ns])
        b, sb, a, chi2 = wls_intercept(x, y, e)
        big[mg] = (b, sb)
        print(f"   m/g={mg:<7} r_inf={b:.4f}+-{sb:.4f}  c={a:+.4f}  "
              f"chi2/dof={chi2:.2f}  [{len(ns)} pts]")

    # ── 4. does the large-N residual vanish as m/g -> 0? ────────────────────
    if len(big) >= 3:
        mgs = np.array(sorted(big))
        dv = np.array([big[m][0] - 1.0 for m in mgs])
        ev = np.array([max(big[m][1], 1e-9) for m in mgs])
        print("\n  does the large-N residual vanish as m/g -> 0?")
        for lbl, xx in (("B + A*(m/g)  ", mgs), ("B + A*(m/g)^2", mgs ** 2)):
            b, sb, a, chi2 = wls_intercept(xx, dv, ev)
            print(f"   {lbl}: B={b:+.5f}+-{sb:.5f} ({abs(b) / sb:.1f}s from 0)"
                  f"  A={a:+.4f}  chi2/dof={chi2:.2f}")
        print("   (linear preferred => the residual is a relative O(mu) correction"
              " to GMOR,\n    which is only the leading term of the mass expansion)")

    # ── 5. per-N form of the same law ───────────────────────────────────────
    print("\n  the same law fitted per N   (ratio - 1 = B + A*(m/g)):")
    for N in NS:
        pts = [mg for mg in MGS if (N, mg) in R]
        if len(pts) < 3:
            continue
        x = np.array(sorted(pts))
        y = np.array([R[(N, m)][0] - 1.0 for m in x])
        e = np.array([max(R[(N, m)][1], 1e-9) for m in x])
        b, sb, a, _ = wls_intercept(x, y, e)
        print(f"   N={N}:  B={b:+.5f}+-{sb:.5f}   A={a:+.4f}")

    # ── 6. the chiral exponent, against both analytic laws ──────────────────
    print("\n  chiral exponent  d ln(M^2/g^2) / d ln(m/g), physical units")
    print("  bars propagate the ensemble bar through the log-ratio treating the")
    print("  two couplings as independent; they are positively correlated (same K")
    print("  window, same fit family), so this OVERstates them.")
    asc = sorted(MGS)
    pairs = list(zip(asc[:-1], asc[1:]))
    head = "".join(f"{f'{a}->{b}':>20}" for a, b in pairs)
    print(f"\n   {'N':>2}{head}{'GMOR':>8}{'DFS 2N/(2N-1)':>15}")
    for N in NS:
        cells, chiral = [], None
        for a, b in pairs:
            if (N, a) not in R or (N, b) not in R:
                cells.append(f"{'--':>20}")
                continue
            Ma, Sa = to_physical(R[(N, a)][2], a), to_physical(R[(N, a)][3], a)
            Mb, Sb = to_physical(R[(N, b)][2], b), to_physical(R[(N, b)][3], b)
            al = np.log(Mb / Ma) / np.log(b / a)
            sa = np.hypot(Sa / Ma, Sb / Mb) / abs(np.log(b / a))
            cells.append(f"{al:13.4f}+-{sa:.4f}")
            if chiral is None:
                chiral = (al, sa)
        line = f"   {N:>2}" + "".join(cells) + f"{1.0:8.3f}{dfs_alpha(N):15.4f}"
        if chiral:
            line += (f"   [{(chiral[0] - 1) / chiral[1]:+.1f}s vs GMOR, "
                     f"{(chiral[0] - dfs_alpha(N)) / chiral[1]:+.1f}s vs DFS]")
        print(line)
    print("\n  GMOR's alpha=1 is exact only as N -> infinity; DFS's 2N/(2N-1)"
          " contains it as\n  the large-N limit.  See docs/weak-coupling-limit.md"
          " for why DFS is a hypothesis\n  under test here rather than a validated"
          " anchor.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
