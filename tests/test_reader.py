"""Tier 0 -- the Fortran reader, including the four qcdf.out parsing hazards."""

import numpy as np
import pytest

from dlcq.dataset import DLCQResult, save, load


pytestmark = pytest.mark.fortran


def test_header_and_structure(fortran_K21):
    r = fortran_K21
    assert (r.N, r.NF, r.B, r.K_code) == (3, 1, 1, 21)
    assert r.rlamb == pytest.approx(0.3325)
    assert r.cutoff == pytest.approx(-1.0)
    assert r.LPN == 0
    assert r.source == "fortran"
    assert (r.numsta_pre, r.numsta_post) == (502, 189)


def test_result_self_consistent(fortran_K21):
    """Momenta odd, summing to K_code; eigenvalues ascending; shapes agree."""
    assert fortran_K21.validate() == []


def test_last_eigenvalue_block_is_the_physical_one(fortran_K21):
    """Hazard 1: 27 COMPUTED EIGENVALUES blocks; the first 26 are norm-matrix.

    The norm-matrix blocks are all ~1e-12.  Picking one of them would give a
    spectrum of essentially zeros.
    """
    ev = fortran_K21.eigenvalues
    assert ev.size == 189
    assert ev[0] == pytest.approx(10.390379999627, rel=1e-12)
    assert ev[1] == pytest.approx(13.822939496992, rel=1e-12)
    assert ev[2] == pytest.approx(15.317122613433, rel=1e-12)
    assert np.all(np.diff(ev) >= -1e-9)


def test_eigenvector_cap_is_recorded(fortran_K21):
    """Hazard 2: IVCTMX=75 caps printed eigenvectors below the eigenvalue count."""
    r = fortran_K21
    assert r.n_eigenvalues == 189
    assert r.n_eigenvectors == 75
    assert r.has_eigenvector(74)
    assert not r.has_eigenvector(75)
    with pytest.raises(IndexError, match="IVCTMX"):
        r.require_eigenvector(75)


def test_orthonormalization_identity(fortran_K21):
    """Hazard 3+4 jointly: ragged PRTSPM rows and the reused Z.

    Z^T N Z = I only holds if the norm matrix and the basis-change matrix were
    both parsed correctly, including the dangling "(" on the ragged final row.
    """
    r = fortran_K21
    gram = r.Z.T @ r.norm @ r.Z
    assert np.abs(gram - np.eye(gram.shape[0])).max() < 1e-10


def test_eigenvectors_orthonormal_under_the_norm_metric(fortran_K21):
    r = fortran_K21
    c = r.c_orig
    assert np.abs(c.T @ r.norm @ c - np.eye(c.shape[1])).max() < 1e-10


def test_norm_matrix_is_symmetric_and_positive_definite(fortran_K21):
    r = fortran_K21
    assert np.allclose(r.norm, r.norm.T, atol=1e-12)
    # Weeding should leave a well-conditioned basis; the minimum here is 6.
    assert np.linalg.eigvalsh(r.norm).min() > 0.5


def test_fock_sectors_are_physical(fortran_K21):
    """A B=1, N=3 baryon has 3 quarks plus qqbar pairs: 3, 5, 7, 9 partons."""
    r = fortran_K21
    assert set(np.unique(r.state_len).tolist()) <= {3, 5, 7, 9}
    for s in range(r.numsta_post):
        L = r.state_len[s]
        nq = int((r.state_types[s, :L] == 1).sum())
        nqbar = int((r.state_types[s, :L] == 0).sum())
        assert nq - nqbar == r.N * r.B


def test_ham_file_matches_out_file():
    """qcdf.ham and qcdf.out must describe the same run."""
    from dlcq.read_fortran import read_ham
    from conftest import find_run

    out = find_run("K21_B1")
    if out is None:
        pytest.skip("no 2K=21 run")
    ham = out.with_name("qcdf.ham")
    if not ham.exists():
        pytest.skip("no qcdf.ham beside the qcdf.out")
    meta, hnu0, hnu = read_ham(ham)
    assert (meta["N"], meta["NF"], meta["B"], meta["K_code"]) == (3, 1, 1, 21)
    assert meta["numsta"] == 189
    assert hnu.shape == (189, 189)
    assert np.allclose(hnu, hnu.T, atol=1e-12)
    assert hnu0.shape == (1, 189)


def test_hdf5_roundtrip(fortran_K21, tmp_path):
    """The data contract must survive serialization intact."""
    path = tmp_path / "r.h5"
    save(fortran_K21, path)
    back = load(path)
    assert back.validate() == []
    for field in ("N", "NF", "B", "K_code", "LPN", "numsta_pre", "numsta_post"):
        assert getattr(back, field) == getattr(fortran_K21, field)
    assert back.rlamb == pytest.approx(fortran_K21.rlamb)
    assert back.source == fortran_K21.source
    np.testing.assert_allclose(back.eigenvalues, fortran_K21.eigenvalues)
    np.testing.assert_allclose(back.c_orig, fortran_K21.c_orig)
    np.testing.assert_allclose(back.norm, fortran_K21.norm)


def test_spectrum_only_mode_agrees(fortran_K21, fortran_K21_nomat):
    np.testing.assert_allclose(fortran_K21_nomat.eigenvalues,
                               fortran_K21.eigenvalues, rtol=0, atol=0)
    assert fortran_K21_nomat.norm is None
