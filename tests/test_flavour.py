"""Multiple flavours (NF > 1).

Flavour was added to ``qcdf.f`` on 5/8/90 and **no figure in the paper uses
it**, so this path shipped with no coverage at all.  It is also the most likely
direction for new work -- unequal quark masses give the analogue of a kaon, and
of a strange baryon, neither of which the paper computes.

The checks here are the ones that do not need a published number to compare
against: exact identities (the two sum rules, flavour symmetry at degenerate
masses), variational inequalities (adding a flavour can only lower the ground
state; a heavier quark can only raise the mass), and agreement between the two
solvers where both can run.

One thing they establish is a defect: ``qcdf.f`` cannot reach flavour
non-singlet channels at all.  See ``test_fortran_cannot_reach_non_singlet``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent

from dlcq.observables import (momentum_sum_rule, number_sum_rule,  # noqa: E402
                              physical_indices)
from dlcq.units import mg_to_lambda  # noqa: E402

LAM = float(mg_to_lambda(1.6))


def solve(NF, B, K, rmq, iflv, ncpus=4):
    from dlcq.read_python import run_python

    return run_python(N=3, NF=NF, B=B, K_code=K, rlamb=LAM,
                      rmq=rmq, iflv=iflv, ncpus=ncpus)


def ground(r):
    idx = physical_indices(r)
    assert idx.size, "no physical states"
    return float(r.eigenvalues[idx[0]]), int(idx[0])


# ── the exact identities have to survive the larger basis ─────────────────

@pytest.mark.slow
@pytest.mark.parametrize("label,NF,B,K,rmq,iflv", [
    ("NF=1 meson",              1, 0, 12, None,       None),
    ("NF=2 neutral",            2, 0, 12, [1.0, 1.0], [0, 0]),
    ("NF=2 neutral, m2=2m1",    2, 0, 12, [1.0, 2.0], [0, 0]),
    ("NF=2 non-singlet",        2, 0, 12, [1.0, 1.0], [1, -1]),
    ("NF=2 non-singlet, m2=2m1", 2, 0, 12, [1.0, 2.0], [1, -1]),
    ("NF=2 baryon uud",         2, 1, 13, [1.0, 1.0], [2, 1]),
    ("NF=2 baryon uus",         2, 1, 13, [1.0, 2.0], [2, 1]),
])
def test_sum_rules_hold_with_flavour(label, NF, B, K, rmq, iflv):
    """Momentum and number sum rules, which flavour must not disturb.

    Every basis state carries total momentum K regardless of how its partons
    are flavoured, so the momentum rule holds state by state and is the
    strongest available check that the larger basis is being built and weighted
    correctly.
    """
    r = solve(NF, B, K, rmq, iflv)
    assert r.numsta_post > 0, f"{label}: no states"
    _, i = ground(r)
    assert momentum_sum_rule(r, i) == pytest.approx(1.0, abs=1e-12), label
    assert number_sum_rule(r, i) == pytest.approx(3 * B, abs=1e-12), label


# ── symmetry and variational inequalities ─────────────────────────────────

@pytest.mark.slow
def test_flavour_symmetry_at_degenerate_masses():
    """With equal masses, flavour is a symmetry and the channels are degenerate.

    The flavour-neutral channel (iflv = [0, 0], which admits u-ubar and d-dbar
    and their mixing) and the non-singlet channel (iflv = [1, -1], a u-dbar
    "kaon") must give the same ground state when m1 = m2.  They are built from
    different numbers of basis states -- 443 against 316 -- so this is a real
    check on the flavour bookkeeping rather than a tautology.
    """
    neutral, _ = ground(solve(2, 0, 12, [1.0, 1.0], [0, 0]))
    nonsinglet, _ = ground(solve(2, 0, 12, [1.0, 1.0], [1, -1]))
    assert neutral == pytest.approx(nonsinglet, rel=1e-9), (
        f"degenerate masses should not split the channels: "
        f"{neutral:.9f} vs {nonsinglet:.9f}")


@pytest.mark.slow
def test_a_second_flavour_can_only_lower_the_ground_state():
    """Variational: a larger basis cannot raise the lowest eigenvalue.

    NF=2 contains every NF=1 configuration plus states with the second flavour,
    so its ground state must sit at or below the NF=1 one -- even when the
    second quark is heavier, since the extra states can simply not be used.
    """
    one, _ = ground(solve(1, 0, 12, None, None))
    both, _ = ground(solve(2, 0, 12, [1.0, 1.0], [0, 0]))
    heavier, _ = ground(solve(2, 0, 12, [1.0, 2.0], [0, 0]))
    assert both <= one * (1 + 1e-9), f"{both:.6f} should be <= {one:.6f}"
    assert heavier <= one * (1 + 1e-9), f"{heavier:.6f} should be <= {one:.6f}"
    assert both <= heavier * (1 + 1e-9), "degenerate should be lowest of the two"


@pytest.mark.slow
@pytest.mark.parametrize("B,K,iflv", [(0, 12, [1, -1]), (1, 13, [2, 1])])
def test_a_heavier_quark_raises_the_mass(B, K, iflv):
    """Monotonicity in the quark mass, for a channel that contains that quark.

    The non-singlet meson and the "uus" baryon each carry exactly one quark of
    flavour 2, so doubling its mass must raise the state.  This is the physics
    a flavour calculation is for, and the paper never does it.
    """
    light, _ = ground(solve(2, B, K, [1.0, 1.0], iflv))
    heavy, _ = ground(solve(2, B, K, [1.0, 2.0], iflv))
    assert heavy > light, f"m2=2m1 gave {heavy:.6f}, not above {light:.6f}"


# ── both solvers, where both can run ──────────────────────────────────────

@pytest.mark.slow
@pytest.mark.fortran
def test_nf2_agrees_between_the_solvers():
    """The NF=2 path in both codes, compared the way the NF=1 path is.

    Agreement is held to the same standard as everywhere else: the ~1e-4
    reproducibility floor that docs/basis-dependence.md establishes for this
    algorithm.  Observed at 2K=12: identical state counts (510 -> 443) and a
    ground state matching to 3e-7.
    """
    from dlcq.providers import FortranProvider

    try:
        f = FortranProvider(allow_run=True).get(3, 2, 0, 12, LAM,
                                                rmq=[1.0, 2.0], iflv=[0, 0])
    except Exception as exc:                       # no gfortran, no binary
        pytest.skip(f"Fortran unavailable: {exc}")

    p = solve(2, 0, 12, [1.0, 2.0], [0, 0])
    assert f.numsta_pre == p.numsta_pre, "state generation differs"
    assert f.numsta_post == p.numsta_post, "weeding differs"

    ef, _ = ground(f)
    ep, _ = ground(p)
    assert abs(ef - ep) / ef < 1e-4, (
        f"ground state: fortran {ef:.9f}, python {ep:.9f}")


# ── a defect this coverage turned up ──────────────────────────────────────

@pytest.mark.slow
@pytest.mark.fortran
def test_fortran_cannot_reach_non_singlet_channels():
    """`qcdf.f` silently overwrites the flavour number of flavour 1.

    Line 253, immediately before ``CALL QCDSTA`` and flagged with the author's
    own ``??????``::

        iflv(1) = N*B

    is unconditional, so whatever the user entered at the ``input flavor qm
    number for flavor 1`` prompt is discarded.  Only the channel with
    ``iflv(1) = N*B`` is reachable; every other one asks for a total quark
    number the code cannot satisfy and generates **zero states rather than an
    error**.

    So the Fortran cannot compute a flavour non-singlet -- no kaon analogue, no
    flavoured baryon other than the one whose first flavour carries all of the
    baryon number.  The Python solver honours the value it is given, and the
    same channel there yields 316 states.

    This test pins the defect: if a future `qcdf.f` ever produces states here,
    the workaround in the Python path can be reconsidered.
    """
    import subprocess

    binary = ROOT / "fortran" / "qcdf"
    if not binary.exists():
        pytest.skip("fortran/qcdf not built")

    outdir = ROOT / "runs" / "nf2_nonsinglet_probe"
    script = ROOT / "fortran" / "run_case.sh"
    subprocess.run(
        ["bash", str(script), str(outdir), "3", "2", "0", f"{LAM:.10g}",
         "-1.0", "0", "12", "1.0", "1", "-1"],
        check=True, capture_output=True,
    )
    text = (outdir / "qcdf.out").read_text(errors="replace")
    assert "NUMBER OF STATES BEFORE WEEDING:      0" in text, (
        "qcdf.f now generates states for a non-singlet channel -- "
        "the iflv(1) override at qcdf.f:253 may have been fixed"
    )

    # and the same channel is perfectly well defined for the Python solver
    p = solve(2, 0, 12, [1.0, 1.0], [1, -1])
    assert p.numsta_post > 300
