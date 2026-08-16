"""Ceiling and parity guards on run_python.

Commit 9487675 measured the failure modes these guards close: past 2K=100 the
colour path returns spurious negative eigenvalues below an otherwise intact
spectrum with no error raised, and a wrong-parity request returns zero states
in ~0.1 s, which reads exactly like a fast successful run.  The slow test here
doubles as the ceiling-validation record: 2K=100 must keep reproducing the
values recorded in that commit.
"""
import numpy as np
import pytest

from dlcq.read_python import run_python
from dlcq.units import mg_to_lambda


def test_beyond_ceiling_raises():
    with pytest.raises(ValueError, match="ceiling"):
        run_python(N=5, NF=1, B=0, K_code=102, rlamb=float(mg_to_lambda(0.1)),
                   LPN=4, assembly="exact", policy="blockwise",
                   solver="sparse", nev=6)


def test_wrong_parity_meson_raises():
    # Mesons need even K_code; 21 is a valid baryon grid point, not a meson one.
    with pytest.raises(ValueError, match="parity"):
        run_python(N=3, NF=1, B=0, K_code=21, rlamb=0.3325)


def test_wrong_parity_baryon_raises():
    # N=3 baryons need odd K_code.
    with pytest.raises(ValueError, match="parity"):
        run_python(N=3, NF=1, B=1, K_code=22, rlamb=0.3325)


def test_right_parity_not_rejected_early():
    # A correct-parity call must get past the guards; 2K=6 meson is a
    # subsecond solve, so run it outright and check it produced states.
    r = run_python(N=2, NF=1, B=0, K_code=6, rlamb=0.3325)
    assert r.numsta_post > 0


@pytest.mark.slow
def test_ceiling_2K100_reproduces_recorded_values():
    """2K=100 is the validated boundary: it must reproduce commit 9487675."""
    from dlcq.observables import physical_indices

    lam = float(mg_to_lambda(0.1))
    expected = {"standard": 0.351212, "improved": 1.051555}
    for ham, want in expected.items():
        r = run_python(N=5, NF=1, B=0, K_code=100, rlamb=lam, LPN=4,
                       assembly="exact", policy="blockwise", solver="sparse",
                       nev=6, keep_norm=False, hamiltonian=ham)
        phys = physical_indices(r)
        got = float(r.eigenvalues[phys[0]])
        assert got == pytest.approx(want, rel=1e-4), (ham, got)


@pytest.mark.slow
def test_beyond_ceiling_probe_warns_and_shows_corruption():
    """With the override, 2K=102 runs — and the corruption is visible."""
    lam = float(mg_to_lambda(0.1))
    with pytest.warns(UserWarning, match="negative eigenvalue"):
        r = run_python(N=5, NF=1, B=0, K_code=102, rlamb=lam, LPN=4,
                       assembly="exact", policy="blockwise", solver="sparse",
                       nev=6, keep_norm=False, allow_beyond_ceiling=True)
    assert float(np.min(r.eigenvalues)) < 0.0
