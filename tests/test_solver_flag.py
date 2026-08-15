"""The ``solver``/``nev`` flag, and the contract a truncated spectrum needs.

Two separate things are being protected.

**The default must not move.**  Every bit-exactness guarantee in the suite --
seven ``array_equal`` assertions in test_kernels, the norm matching the 1990
Fortran to exactly 0.0, the element-identical state generation -- rests on
``solver="dense"`` staying the default in all three places it can be set.  That
is one assertion, and it is stronger and cheaper than auditing the tests that
depend on it.

**A truncated spectrum must not pass for a complete one.**  ``eigenvalues[i]``
is still the i-th lowest level either way, so nothing looks wrong near the
bottom; the difference only shows when a consumer indexes deep, and by then it
has already drawn the wrong state.  ``n_orth`` records what a full spectrum
would have been so the two cases stay distinguishable, and
``require_physical_index`` turns the deep case into an error rather than a
plausible-looking panel.
"""

from __future__ import annotations

import inspect
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
PYDIR = ROOT / "python"
for p in (str(ROOT), str(PYDIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

pytest.importorskip("numba")

from dlcq.dataset import DLCQResult, load, save                # noqa: E402
from dlcq.observables import require_physical_index            # noqa: E402
from dlcq.providers import PythonProvider, _tag                # noqa: E402
from dlcq.read_python import _solve_eigen, run_python          # noqa: E402

# Every test here drives the solver end to end -- 10 run_python/provider calls
# across this file -- so the whole module is the slow tier.  Without this,
# ``pytest -m "not slow"`` still runs them: it took 30 minutes and was killed by
# a timeout, against seconds for the genuinely fast tests.
pytestmark = pytest.mark.slow


BASE = dict(N=3, NF=1, B=1, K_code=15, rlamb=0.3325, cutoff=-1.0, LPN=0,
            ncpus=4, assembly="exact", policy="blockwise")


# ── the default ───────────────────────────────────────────────────────────

def test_dense_is_the_default_everywhere():
    """The single fact the whole bit-exactness surface rests on."""
    assert inspect.signature(run_python).parameters["solver"].default == "dense"
    assert (inspect.signature(PythonProvider.__init__)
            .parameters["solver"].default == "dense")
    assert PythonProvider().solver == "dense"
    assert PythonProvider().nev is None


def test_the_default_cache_tag_is_unchanged():
    """Adding the flag must not orphan the existing cache.

    ``md5("exact:fortran")[:8]`` is the suffix on the great majority of the
    files in ``runs/python_cache`` -- several GB, hours of compute.  Appending
    ":dense" unconditionally would have changed every one of them.
    """
    before = _tag(3, 1, 1, 21, 0.3325, -1.0, 0, extra="exact:fortran")
    p = PythonProvider()
    extra = f"{p.assembly}:{p.policy}"
    assert p.solver == "dense" and p.nev is None
    assert _tag(3, 1, 1, 21, 0.3325, -1.0, 0, extra=extra) == before


# ── validation ────────────────────────────────────────────────────────────

def test_unknown_solver_is_an_error():
    with pytest.raises(ValueError, match="solver"):
        run_python(solver="magic", **BASE)


def test_unknown_policy_is_rejected_before_any_work():
    """It used to surface only after the norm build -- minutes, at large K."""
    with pytest.raises(ValueError, match="policy"):
        run_python(**{**BASE, "policy": "magic"})


@pytest.mark.parametrize("assembly,policy", [
    ("fortran", "blockwise"),
    ("exact", "fortran"),
    ("exact", "spectral"),
    ("fortran", "fortran"),
    ("fortran", "spectral"),
])
def test_sparse_rejects_every_unsupported_combination(assembly, policy):
    """One legal cell of six, and the other five must say so rather than run.

    A silent fallback to dense would be worse than an error: the run would
    succeed, take the memory the sparse path exists to avoid, and the caller
    would conclude the flag does not work.
    """
    kw = {**BASE, "assembly": assembly, "policy": policy, "solver": "sparse"}
    with pytest.raises(ValueError, match="sparse"):
        run_python(**kw)


def test_nev_must_be_positive():
    with pytest.raises(ValueError, match="nev"):
        run_python(**{**BASE, "nev": 0})


# ── the eigensolver itself ────────────────────────────────────────────────

def _clustered(n=200, seed=0):
    """Symmetric, bounded below, with a deliberately tight low end."""
    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.standard_normal((n, n)))
    w = np.concatenate([np.linspace(1.0, 1.02, 12), np.linspace(2.0, 40.0,
                                                                n - 12)])
    return (q * w) @ q.T


def test_sparse_matches_dense_on_a_clustered_spectrum():
    """The failure that matters is a *skipped* level, not an imprecise one.

    A missed level shifts every subsequent index by one, so the ground state a
    figure plots becomes the second state -- and the Richardson fit through it
    still looks perfectly smooth.  Comparing the whole sorted vector, rather
    than just the lowest value, is what catches that.
    """
    A = _clustered()
    w_dense, _ = _solve_eigen(A, solver="dense", nev=20)
    w_sparse, _ = _solve_eigen(A, solver="sparse", nev=20)
    np.testing.assert_allclose(w_sparse, w_dense, rtol=0, atol=1e-10)


def test_sparse_is_deterministic():
    """``eigsh`` defaults to a random v0; with a cache that is a correctness bug.

    Two runs of identical parameters would return different bits, and whichever
    reached the cache first would win silently thereafter.
    """
    A = _clustered(seed=1)
    a, _ = _solve_eigen(A, solver="sparse", nev=15)
    b, _ = _solve_eigen(A, solver="sparse", nev=15)
    np.testing.assert_array_equal(a, b)


def test_dense_with_nev_is_the_full_spectrum_sliced():
    """Which is what makes it the shape-exact reference for the sparse path."""
    A = _clustered(seed=2)
    w_all, z_all = _solve_eigen(A, solver="dense")
    w_cut, z_cut = _solve_eigen(A, solver="dense", nev=25)
    np.testing.assert_array_equal(w_cut, w_all[:25])
    np.testing.assert_array_equal(z_cut, z_all[:, :25])


def test_sparse_falls_back_when_k_approaches_n():
    """ARPACK needs k < n, and below that size dense is faster anyway."""
    A = _clustered(n=30)
    w, z = _solve_eigen(A, solver="sparse", nev=25)
    assert w.size == 25
    np.testing.assert_allclose(w, _solve_eigen(A, solver="dense", nev=25)[0],
                               rtol=0, atol=1e-10)


# ── the truncation contract ───────────────────────────────────────────────

def test_n_orth_records_the_full_dimension():
    full = run_python(**BASE)
    cut = run_python(**{**BASE, "solver": "sparse", "nev": 12})
    assert full.n_orth == cut.n_orth == full.n_eigenvalues
    assert not full.spectrum_is_truncated
    assert cut.spectrum_is_truncated
    assert cut.n_eigenvalues == 12
    np.testing.assert_allclose(cut.eigenvalues, full.eigenvalues[:12],
                               rtol=0, atol=1e-9)


def test_a_truncated_result_is_still_self_consistent():
    """``validate`` reports structural problems; truncation is not one."""
    cut = run_python(**{**BASE, "solver": "sparse", "nev": 12})
    assert cut.validate() == []


def test_validate_catches_more_eigenvalues_than_n_orth():
    r = DLCQResult(N=3, NF=1, B=1, K_code=15, rlamb=0.3325, n_orth=5,
                   eigenvalues=np.arange(9.0))
    assert any("exceeds" in p for p in r.validate())


def test_legacy_files_get_n_orth_repaired_on_load(tmp_path):
    """Several hundred cached results predate the field.

    Left at 0, ``spectrum_is_truncated`` would answer False for the wrong
    reason -- indistinguishable from a genuinely complete run.  Every one of
    them kept the full spectrum, so Z's column count is the answer.
    """
    import h5py
    from dlcq.dataset import _ARRAYS, _SCALARS

    r = run_python(**BASE)
    path = tmp_path / "legacy.h5"
    save(r, path)
    with h5py.File(path, "a") as f:          # strip the field, as a v1 file
        del f.attrs["n_orth"]
    back = load(path)
    assert back.n_orth == r.Z.shape[1]
    assert not back.spectrum_is_truncated


def test_require_physical_index_fails_loudly_past_the_truncation():
    cut = run_python(**{**BASE, "solver": "sparse", "nev": 6})
    assert require_physical_index(cut, 0) == 0
    with pytest.raises(IndexError, match="Raise nev"):
        require_physical_index(cut, 50)


def test_require_physical_index_distinguishes_short_from_truncated():
    """A complete-but-small run gets a different message from a truncated one."""
    full = run_python(**BASE)
    with pytest.raises(IndexError) as exc:
        require_physical_index(full, full.n_eigenvalues + 5)
    assert "Raise nev" not in str(exc.value)
