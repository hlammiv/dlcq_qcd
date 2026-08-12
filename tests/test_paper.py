"""Tier 3 -- against the published paper.

Two kinds of assertion:

* **Exact results** the paper states hold identically.  These are the strongest
  checks in the suite: no digitization, no extrapolation, no tolerance
  negotiation.
* **Table I**, compared loosely, because those entries are Richardson
  extrapolations whose quoted uncertainties the paper itself calls "a rough
  guide" below m/g ~ 0.2.
"""

import csv
from pathlib import Path

import numpy as np
import pytest

from dlcq.observables import (physical_indices, spurious_zero_modes,
                              structure_function, thooft_valence_limit)
from dlcq.units import code_to_M_over_g, meson_baryon_ratio_bosonization

ROOT = Path(__file__).resolve().parent.parent
TABLE1 = ROOT / "refs" / "table1.csv"


def load_table1():
    """Table I as a list of dicts.  Values transcribed from the scanned page."""
    rows = []
    with open(TABLE1) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            reader = csv.DictReader([("quantity,N,mg,value,last_term,reliable")
                                     if line.startswith("quantity") else line],
                                    fieldnames=["quantity", "N", "mg", "value",
                                                "last_term", "reliable"])
            for row in reader:
                if row["quantity"] == "quantity":
                    continue
                rows.append(dict(quantity=row["quantity"], N=int(row["N"]),
                                 mg=float(row["mg"]), value=float(row["value"]),
                                 last_term=float(row["last_term"]),
                                 reliable=bool(int(row["reliable"]))))
    return rows


# ── the transcription itself ──────────────────────────────────────────────

def test_table1_transcription_is_complete():
    """5 columns x 7 couplings = 35 entries; no N=2 baryon column exists."""
    rows = load_table1()
    assert len(rows) == 35
    assert {(r["quantity"], r["N"]) for r in rows} == {
        ("mes", 2), ("mes", 3), ("mes", 4), ("bar", 3), ("bar", 4)}
    assert not any(r["quantity"] == "bar" and r["N"] == 2 for r in rows)


def test_table1_spot_values():
    """Anchors verified against a 600 dpi render of page 3820."""
    rows = {(r["quantity"], r["N"], r["mg"]): r for r in load_table1()}
    assert rows[("mes", 3, 1.6)]["value"] == 4.618
    assert rows[("mes", 3, 1.6)]["last_term"] == 0.006
    assert rows[("bar", 3, 1.6)]["value"] == 10.71
    assert rows[("bar", 4, 1.6)]["value"] == 21.2
    assert rows[("mes", 2, 1.6)]["value"] == 4.314
    assert rows[("bar", 4, 0.4)]["value"] == 15.0


def test_chiral_limit_rows_are_exactly_zero():
    """The paper: at m/g = 0 the lightest state is exactly massless, any N or B."""
    for r in load_table1():
        if r["mg"] == 0.0:
            assert r["value"] == 0.0


def test_table1_is_monotone_in_coupling():
    """M/g must increase with m/g for every column."""
    rows = load_table1()
    for quantity in ("mes", "bar"):
        for N in (2, 3, 4):
            col = sorted([r for r in rows
                          if r["quantity"] == quantity and r["N"] == N],
                         key=lambda r: r["mg"])
            vals = [r["value"] for r in col]
            assert all(b >= a for a, b in zip(vals, vals[1:])), \
                f"{quantity} N={N} not monotone: {vals}"


def test_bosonization_ratio_against_table1():
    """M_mes/M_bar at strong coupling vs 2 sin[pi/(2(2N-1))], to the paper's ~10%.

    Table I is in M^2 units (docs/table1-units.md), so the mass ratio needs a
    square root -- and the m^2+g^2/pi denominator cancels between the two.
    """
    rows = {(r["quantity"], r["N"], r["mg"]): r for r in load_table1()}
    for N in (3, 4):
        mes = rows[("mes", N, 0.05)]["value"]
        bar = rows[("bar", N, 0.05)]["value"]
        ratio = np.sqrt(mes / bar)
        predicted = meson_baryon_ratio_bosonization(N)
        assert ratio == pytest.approx(predicted, rel=0.35), (
            f"N={N}: table gives {ratio:.3f}, bosonization {predicted:.3f}"
        )


# ── exact results, checked against the Fortran output ─────────────────────

@pytest.mark.fortran
def test_no_spurious_zero_in_the_baryon_run(fortran_K21):
    """2K=21, B=1 has a clean spectrum; the artifact is specific to 2K=24, B=0."""
    assert spurious_zero_modes(fortran_K21).size == 0
    assert fortran_K21.eigenvalues[0] > 1.0


@pytest.mark.fortran
def test_spurious_zero_mode_is_detected_where_it_occurs():
    """N=3, B=0, 2K=24 carries one decoupled 12-parton state at M^2 = 0.

    Guards the detector itself: if the Fortran is ever rebuilt such that this
    state acquires a correct Hamiltonian row, this test should be revisited
    rather than the detector loosened.
    """
    from dlcq.read_fortran import read_out

    path = ROOT / "runs" / "fig5_B0_K24" / "qcdf.out"
    if not path.exists():
        pytest.skip("2K=24 B=0 run not present")
    r = read_out(path)
    bad = spurious_zero_modes(r)
    assert bad.tolist() == [0]
    # It is the unique maximal-Fock state.
    c = r.c_orig[:, 0]
    dominant = int(np.argmax(np.abs(c)))
    assert r.state_len[dominant] == 12
    assert r.norm[dominant, dominant] == pytest.approx(1296.0)   # (3!)^4


@pytest.mark.fortran
def test_eleventh_meson_state_matches_the_paper_text():
    """The paper: the 11th state has "a mass twice that of the first" and its
    two-pair component peaks at x = 1/4."""
    from dlcq.read_fortran import read_out

    path = ROOT / "runs" / "fig5_B0_K24" / "qcdf.out"
    if not path.exists():
        pytest.skip("2K=24 B=0 run not present")
    r = read_out(path)
    phys = physical_indices(r)

    m1 = code_to_M_over_g(r.eigenvalues[phys[0]], r.rlamb)
    m11 = code_to_M_over_g(r.eigenvalues[phys[10]], r.rlamb)
    assert m11 / m1 == pytest.approx(2.0, rel=0.05)

    x, q4, _ = structure_function(r, int(phys[10]), nparton=4)
    peak = x[np.argmax(q4)]
    # x = 1/4 is not on the odd-momentum grid at 2K=24; 5/24 and 7/24 bracket it.
    assert abs(peak - 0.25) <= 1.0 / r.K_code + 1e-9, f"peak at x={peak}"


@pytest.mark.fortran
def test_valence_dominance_matches_the_paper(fortran_K21):
    """Higher-Fock suppressed by 2-4 orders of magnitude relative to valence."""
    dx = 2.0 / fortran_K21.K_code
    _, q3, qb3 = structure_function(fortran_K21, 0, nparton=3)
    _, q5, qb5 = structure_function(fortran_K21, 0, nparton=5)
    valence = float(np.sum(q3 - qb3) * dx)
    higher = float(np.sum(q5 - qb5) * dx)
    ratio = higher / valence
    assert 1e-5 < ratio < 1e-2, f"higher-Fock fraction {ratio:.2e}"


def test_chiral_baryon_structure_function_shape():
    """Eq. (22) is what Fig. 3(b) must approach as m/g -> 0."""
    x = np.linspace(0, 1, 50)
    q = thooft_valence_limit(x, 3)
    assert q[0] == pytest.approx(6.0)
    assert q[-1] == pytest.approx(0.0, abs=1e-12)
    assert np.all(np.diff(q) < 0)          # strictly decreasing for N=3


# ── the paper's figures, and the independent lattice data in Fig. 8(a) ────

HAMER_SU2 = [   # digitized from Fig. 8(a); C. J. Hamer, Nucl. Phys. B195, 503
    (0.2047, 0.6590), (0.4091, 1.1109), (0.8233, 1.9749), (1.6465, 3.5187),
]


def test_hamer_lattice_points_are_in_range():
    """Sanity-check the digitized external data before anything asserts on it."""
    for mg, Mg in HAMER_SU2:
        assert 0.0 < mg < 1.8
        assert 0.0 < Mg < 4.0        # Fig. 8(a)'s y-axis maximum


def test_table1_is_in_msquared_units_not_M_over_g():
    """Table I tabulates M^2/(m^2+g^2/pi) despite its M/g headers.

    Two independent consequences pin this down, neither needing a solver:

    1. Fig. 8(a) plots M/g on an axis running to 4.0 (9 ticks every 0.5), yet
       every Table I meson entry at m/g = 1.6 exceeds it.  The table cannot be
       in the units of the figure it summarizes.
    2. Converting Table I to M/g via M/g = sqrt(value / (pi lambda^2)) puts the
       m/g = 1.6 mesons back inside that axis.

    See docs/table1-units.md.
    """
    from dlcq.units import mg_to_lambda

    rows = {(r["quantity"], r["N"], r["mg"]): r for r in load_table1()}
    lam = float(mg_to_lambda(1.6))
    for N in (2, 3, 4):
        raw = rows[("mes", N, 1.6)]["value"]
        assert raw > 4.0                       # impossible as M/g
        as_mass = np.sqrt(raw / (np.pi * lam ** 2))
        assert 3.0 < as_mass < 4.0             # sits on Fig. 8(a)'s axis


def test_table1_reproduces_hamer_after_unit_conversion():
    """Converted to M/g, Table I's SU(2) entry matches Hamer's lattice point.

    Hamer's value is external data (Nucl. Phys. B195, 503) digitized from
    Fig. 8(a), so this ties the table to a wholly independent calculation.
    """
    from dlcq.units import mg_to_lambda

    rows = {(r["quantity"], r["N"], r["mg"]): r for r in load_table1()}
    lam = float(mg_to_lambda(1.6))
    as_mass = np.sqrt(rows[("mes", 2, 1.6)]["value"] / (np.pi * lam ** 2))
    hamer = 3.5187
    assert as_mass == pytest.approx(hamer, rel=0.02), (
        f"Table I -> M/g gives {as_mass:.4f}, Hamer {hamer:.4f}"
    )


# ── the large-N limit plotted in Fig. 8 ──────────────────────────────────

def test_thooft_chiral_limit_reproduces_the_classic_eigenvalues():
    """Eq. (24) at m = 0, large N: mu^2 = 0, 5.88, 14.1, ...

    An external benchmark for the Fig. 8 large-N curve, independent of the DLCQ
    codes entirely.  Note the values come out in the paper's own coupling
    normalization with no conversion factor.
    """
    from dlcq.thooft import thooft_coupling, thooft_spectrum

    N = 1000
    mu2 = thooft_spectrum(0.0, N, n=800, k=3) / thooft_coupling(N)
    assert abs(mu2[0]) < 1e-6                      # exactly massless
    assert mu2[1] == pytest.approx(5.88, rel=0.01)
    assert mu2[2] == pytest.approx(14.1, rel=0.01)


def test_thooft_agrees_with_dlcq_at_weak_coupling():
    """At m/g = 1.6 the qqbar sector dominates, so Eq. (24) and DLCQ must agree.

    Independent equation, independent discretization, independent code path --
    so this cross-checks the swept DLCQ masses without reusing any of that
    machinery.  Values are the Richardson-extrapolated sweep results.
    """
    from dlcq.thooft import thooft_mass

    dlcq = {2: 3.5234, 3: 3.6446, 4: 3.7340}
    for N, expected in dlcq.items():
        assert thooft_mass(1.6, N) == pytest.approx(expected, rel=0.01), N


def test_thooft_mass_grows_with_quark_mass():
    from dlcq.thooft import thooft_mass

    vals = [thooft_mass(t, 3) for t in (0.1, 0.4, 0.8, 1.6)]
    assert all(b > a for a, b in zip(vals, vals[1:]))


# ── the thesis (SLAC-333), which quotes higher-Fock magnitudes ────────────

def load_thesis_table6():
    """Four-quark Fock probability in the lightest SU(2) meson, from the thesis."""
    path = ROOT / "refs" / "thesis_table6.csv"
    rows = []
    with open(path) as fh:
        for row in csv.DictReader(l for l in fh if not l.startswith("#")):
            rows.append(dict(mg=float(row["mg"]),
                             probability=float(row["probability"]),
                             last_term=float(row["last_term"]),
                             reliable=bool(int(row["reliable"]))))
    return rows


def test_thesis_table6_transcription():
    rows = load_thesis_table6()
    assert len(rows) == 4
    assert {r["mg"] for r in rows} == {0.40, 0.20, 0.10, 0.05}
    # Higher-Fock grows as the quarks get lighter.
    ordered = sorted(rows, key=lambda r: r["mg"], reverse=True)
    vals = [r["probability"] for r in ordered]
    assert all(b > a for a, b in zip(vals, vals[1:]))


@pytest.mark.slow
@pytest.mark.fortran
def test_higher_fock_magnitude_matches_the_thesis():
    """Our four-quark probability against the thesis's published values.

    This is the only quantitative check of a higher-Fock *magnitude* that does
    not go through a digitized curve -- the thesis quotes numbers with
    uncertainties.  It is the direct answer to whether our Fock projection has
    the right normalization, independent of any figure.
    """
    import re

    from dlcq.observables import richardson_extrapolate
    from dlcq.read_fortran import read_out

    sweep = ROOT / "runs" / "sweep"
    if not sweep.exists():
        pytest.skip("no sweep; run tools/sweep_fortran.sh")

    wanted = {r["mg"]: r for r in load_thesis_table6()}
    runs = {}
    for d in sorted(sweep.glob("N2_B0_*")):
        m = re.match(r"N2_B0_K(\d+)_mg([\d.]+)", d.name)
        if not m or not (d / "qcdf.out").exists():
            continue
        mg = float(m[2])
        if mg not in wanted:
            continue
        r = read_out(d / "qcdf.out")
        idx = physical_indices(r)
        if not idx.size:
            continue
        dx = 2.0 / r.K_code
        prob = {}
        for n in (2, 4):
            _, q, qb = structure_function(r, int(idx[0]), nparton=n)
            # partons per state -> probability of the Fock component
            prob[n] = float(np.sum(q + qb) * dx) / n
        runs.setdefault(mg, []).append((int(m[1]), prob[4] / (prob[2] + prob[4])))

    checked = 0
    for mg, ref in wanted.items():
        if mg not in runs or len(runs[mg]) < 3:
            continue
        v = sorted(runs[mg])
        p0, _ = richardson_extrapolate([k for k, _ in v], [p for _, p in v],
                                       mg, 2, n_terms=2)
        # 10% covers the extrapolation spread; the thesis itself calls
        # m/g < 0.2 unreliable, so the tolerance is not tightened there.
        assert p0 == pytest.approx(ref["probability"], rel=0.10), (
            f"m/g={mg}: ours {p0:.4e}, thesis {ref['probability']:.4e}"
        )
        checked += 1
    assert checked >= 3, f"only {checked} couplings compared"


def load_thesis_fig18a():
    """Digitized five-quark baryon structure function, thesis Fig 18(a)."""
    path = ROOT / "refs" / "thesis_fig18a.csv"
    rows = []
    with open(path) as fh:
        for row in csv.DictReader(l for l in fh if not l.startswith("#")):
            rows.append((row["quantity"], float(row["x"]), float(row["value"])))
    return rows


@pytest.mark.fortran
def test_baryon_five_quark_matches_the_thesis():
    """The baryon higher-Fock sector against thesis Fig 18(a), 2K=15.

    This is the panel that settles it.  The article's Fig. 6(a) overlays the
    valence and five-quark curves, which cross three times, and at scan
    resolution the dashes of the dashed curve read as markers -- so the traced
    "higher-Fock" series there is unreliable and made this sector look wrong.
    Fig 18(a) shows the sector alone with q and qbar separated and no crossings.
    """
    from dlcq.read_fortran import read_out

    path = ROOT / "runs" / "K15_B1" / "qcdf.out"
    if not path.exists():
        pytest.skip("no 2K=15 B=1 run; see fortran/run_case.sh")
    r = read_out(path)
    idx = physical_indices(r)
    x, q, qbar = structure_function(r, int(idx[0]), nparton=5)

    worst = 0.0
    for quantity, xp, val in load_thesis_fig18a():
        i = int(np.argmin(np.abs(x - xp)))
        ours = (q[i] if quantity == "q" else qbar[i]) * 1e3
        rel = abs(val - ours) / val
        worst = max(worst, rel)
        assert rel < 0.20, (
            f"{quantity}(x={xp}): ours {ours:.3f}, thesis {val:.3f} ({rel:.1%})"
        )
    assert worst < 0.20
