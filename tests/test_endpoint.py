"""van de Sande's endpoint subtraction: the reduction, the contract, the physics.

The load-bearing test here is :func:`test_reduces_to_standard_at_zero_exponent`.
The improved self-inertia is a *rewriting* of the standard one with a weight
inserted, so at zero exponent the weight is 1, the compensating integral
vanishes, and the two Hamiltonians must agree to machine precision.  Anything
wrong with the pair decomposition, the colour normalisation or the units shows
up there first, which is why it gates the rest.

See ``dlcq/endpoint.py`` and ``docs/weak-coupling-limit.md``.
"""
from __future__ import annotations

import inspect
import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / "python"))

from dlcq.endpoint import endpoint_integral, improved_selfen
from dlcq.providers import PythonProvider, _tag
from dlcq.read_python import run_python
from dlcq.units import endpoint_exponent, mg_to_lambda

import qcdf as base  # noqa: E402
import qcdf_opt as opt  # noqa: E402


def _ham(N, B, K, LPN, selfen):
    """The interacting Hamiltonian for one sector, with a given selfen table."""
    p = base.Params()
    p.N, p.NF, p.B, p.K, p.rlamb = N, 1, B, K, 0.3325
    p.cutoff, p.LPN = -1.0, LPN
    p.iflv[0] = N * B
    st = base.StateData()
    base.qcdsta(p, st, base.PermTables(), base.FlavorTables())
    n = st.numsta
    if n == 0:
        pytest.skip(f"N={N} B={B} 2K={K} LPN={LPN} generates no states")
    lengths = sorted({int(st.mstinf[s, 1]) for s in range(n)})
    ham = opt.build_matrices(1, st.mstate, st.mstinf[:n].copy(), n, N, 1, B, K,
                             selfen, p.cbreak, 4, backend="thread")[1]
    return ham, lengths


# ── the gate ──────────────────────────────────────────────────────────────

@pytest.mark.parametrize("N,K", [(3, 10), (3, 16), (3, 24), (2, 12), (2, 30),
                                 (4, 14)])
def test_reduces_to_standard_at_zero_exponent(N, K):
    """At b = 0 the improved Hamiltonian *is* the standard one.

    Not merely close: the weight is identically 1 and ``I(x)`` is identically
    zero, so the only difference is how the same total is distributed between
    a pair's two partons -- and the solver only ever sees the sum.
    """
    std, lens = _ham(N, 0, K, 2, opt.compute_selfen(N))
    imp, _ = _ham(N, 0, K, 2, improved_selfen(N, K, 0.0))
    assert lens == [2]
    scale = max(np.max(np.abs(std)), 1e-300)
    assert np.max(np.abs(std - imp)) / scale < 1e-13


def test_the_weight_actually_changes_something():
    """Guard against the reduction test passing vacuously."""
    N, K = 3, 24
    b = endpoint_exponent(0.1, N)
    assert b > 0
    std, _ = _ham(N, 0, K, 2, opt.compute_selfen(N))
    imp, _ = _ham(N, 0, K, 2, improved_selfen(N, K, b))
    scale = np.max(np.abs(std))
    assert np.max(np.abs(std - imp)) / scale > 1e-3


# ── the integral ──────────────────────────────────────────────────────────

def test_endpoint_integral_vanishes_at_zero_exponent():
    for x in (0.1, 0.3, 0.5, 0.9):
        assert endpoint_integral(x, 0.0) == 0.0


def test_endpoint_integral_is_symmetric():
    """I(x) is built from w(x) = (x(1-x))^b, which is symmetric about 1/2."""
    for x in (0.1, 0.25, 0.4):
        a = endpoint_integral(x, 0.3)
        c = endpoint_integral(1.0 - x, 0.3)
        assert a == pytest.approx(-c, rel=1e-6, abs=1e-9) or \
               a == pytest.approx(c, rel=1e-6, abs=1e-9)


# ── the flag contract ─────────────────────────────────────────────────────

def test_standard_is_the_default_everywhere():
    assert (inspect.signature(run_python)
            .parameters["hamiltonian"].default == "standard")
    assert (inspect.signature(PythonProvider.__init__)
            .parameters["hamiltonian"].default == "standard")
    assert PythonProvider().hamiltonian == "standard"


def test_the_default_cache_tag_is_unchanged():
    """The 31 GB on disk must not re-key because a new flag exists.

    md5("exact:fortran")[:8] = 31275504 and md5("exact:blockwise:sparse")[:8]
    = 2900464f are the suffixes on 811 and 1246 cached files respectively.
    """
    before = _tag(3, 1, 1, 21, 0.3325, -1.0, 0, extra="exact:fortran")
    assert _tag(3, 1, 1, 21, 0.3325, -1.0, 0,
                extra=PythonProvider()._extra()) == before
    sparse = PythonProvider(assembly="exact", policy="blockwise",
                            solver="sparse")
    assert sparse._extra() == "exact:blockwise:sparse"


def test_improved_gets_its_own_cache_namespace():
    a = PythonProvider()._extra()
    b = PythonProvider(hamiltonian="improved")._extra()
    assert a != b and b.endswith(":improved")


def test_get_and_spectrum_agree_on_the_tag():
    """They built the discriminator separately once; a flag added to one only
    would have sent ``spectrum`` to a different file than ``get`` wrote."""
    src = inspect.getsource(PythonProvider)
    assert src.count("self._extra()") >= 2


def test_unknown_hamiltonian_is_rejected_before_any_work():
    with pytest.raises(ValueError, match="hamiltonian"):
        run_python(N=3, NF=1, B=0, K_code=10, rlamb=0.3325, LPN=2,
                   hamiltonian="magic")


def test_improved_refuses_more_than_two_partons():
    """A valence baryon has no partner momentum fixed by conservation."""
    with pytest.raises(ValueError, match="two-parton"):
        run_python(N=3, NF=1, B=1, K_code=15, rlamb=0.3325, LPN=3, ncpus=2,
                   assembly="exact", policy="blockwise",
                   hamiltonian="improved")


def test_improved_refuses_an_untruncated_basis():
    """LPN=0 means *no* truncation, not valence only -- an easy mistake."""
    with pytest.raises(ValueError, match="two-parton"):
        run_python(N=3, NF=1, B=0, K_code=12, rlamb=0.3325, LPN=0, ncpus=2,
                   assembly="exact", policy="blockwise",
                   hamiltonian="improved")


# ── the physics ───────────────────────────────────────────────────────────

@pytest.mark.slow
def test_improved_fixes_the_chiral_exponent():
    """The point of the whole exercise.

    ``M^2`` scales as ``(m/g)^alpha``.  The physical answer is alpha = 1 -- GMOR,
    and van de Sande's Eq. (7), ``M^2 = 2 pi g mu / sqrt 3``.  Standard DLCQ
    returns ~2 at every K, which is the artifact ``docs/weak-coupling-limit.md``
    traces.  Improved DLCQ must land near 1, and must do so *stably in K*, since
    a value that drifted would just be a fit artefact.
    """
    def M2(K, mg, ham):
        r = run_python(N=3, NF=1, B=0, K_code=K, rlamb=float(mg_to_lambda(mg)),
                       LPN=2, ncpus=4, assembly="exact", policy="blockwise",
                       hamiltonian=ham)
        return float(r.eigenvalues[0])

    for K in (30, 40):
        std = np.log(M2(K, 0.1, "standard") / M2(K, 0.05, "standard")) / np.log(2)
        imp = np.log(M2(K, 0.1, "improved") / M2(K, 0.05, "improved")) / np.log(2)
        assert 1.85 < std < 2.05, f"standard alpha moved: {std}"
        assert 0.95 < imp < 1.20, f"improved alpha off: {imp}"


@pytest.mark.slow
def test_improved_converges_far_faster_in_K():
    """Standard drifts ~50% over 2K = 10..60; improved should barely move."""
    def M2(K):
        out = {}
        for ham in ("standard", "improved"):
            r = run_python(N=3, NF=1, B=0, K_code=K,
                           rlamb=float(mg_to_lambda(0.1)), LPN=2, ncpus=4,
                           assembly="exact", policy="blockwise",
                           hamiltonian=ham)
            out[ham] = float(r.eigenvalues[0])
        return out

    lo, hi = M2(10), M2(60)
    drift = {h: abs(hi[h] - lo[h]) / hi[h] for h in lo}
    assert drift["standard"] > 0.25
    assert drift["improved"] < 0.05
    assert drift["improved"] < drift["standard"] / 5
