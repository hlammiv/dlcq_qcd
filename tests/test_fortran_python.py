"""Tier 2 -- Python solver against the Fortran.

Two tolerance levels, and the reason for the split is the whole point:

* The **matrices** agree to machine precision.  Restricted to the Fortran's own
  retained basis and contracted with the Fortran's own Z, the norm matrix is
  bit-identical and the Hamiltonian agrees to ~1e-14.  This is where the real
  verification happens: it exercises the state generator, the diagrammatic
  color-contraction engine, all seven four-point vertices and both self-energy
  terms.

* The **spectrum** can only be compared at ~1e-3, because qcdf.f adds the free
  Hamiltonian to the diagonal only of the orthonormal-basis matrix, which is a
  basis-dependent step, and the basis is not uniquely determined.  Recompiling
  the unmodified qcdf.f with -O2 instead of -O0 shifts its own ground state by
  2.9e-5 relative and changes the retained-state count from 189 to 190.

  See docs/basis-dependence.md.  Tightening this tolerance would not be
  testing the Python port; it would be testing which way the compiler rounded.
"""

import numpy as np
import pytest

from dlcq.units import mg_to_lambda


pytestmark = pytest.mark.fortran

# The floor imposed by the Fortran's basis-dependent assembly, measured by
# recompiling qcdf.f at -O0 vs -O2 (2.9e-5 on the ground state, 6.6e-2 worst
# case across the full spectrum).  Applied to the low-lying states that the
# figures actually use.
# (the end-to-end spectrum tolerance is derived at runtime by _basis_spread)
MATRIX_ATOL = 1e-10


@pytest.fixture(scope="module")
def matched(fortran_K21):
    """Python matrices restricted to the Fortran's exact retained basis.

    States are matched on (parton types, momenta), which is unique, so this
    compares like with like rather than relying on both codes happening to
    order or weed identically.
    """
    import os
    import sys
    from pathlib import Path

    root = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(root / "python"))
    os.environ["QCDF_NCPUS"] = "4"
    import qcdf
    import qcdf_opt

    f = fortran_K21
    p = qcdf.Params()
    p.N, p.NF, p.B, p.K, p.rlamb = f.N, f.NF, f.B, f.K_code, f.rlamb
    p.cutoff, p.LPN = f.cutoff, f.LPN
    p.iflv[0] = f.N * f.B

    st = qcdf.StateData()
    qcdf.qcdsta(p, st, qcdf.PermTables(), qcdf.FlavorTables())
    selfen = qcdf_opt.compute_selfen(f.N)

    def key(s):
        loc = int(st.mstinf[s, 0]) - 1
        L = int(st.mstinf[s, 1])
        return (tuple(st.mstate[loc, :L]), tuple(st.mstate[loc + 2, :L]))

    table = {}
    for s in range(st.numsta):
        table.setdefault(key(s), []).append(s)

    idx = []
    for s in range(f.numsta_post):
        L = f.state_len[s]
        k = (tuple(f.state_types[s, :L]), tuple(f.state_moms[s, :L]))
        assert k in table and table[k], f"Fortran state {s} absent from Python basis"
        idx.append(table[k].pop(0))
    idx = np.array(idx)

    _, _, hnorm = qcdf_opt.build_matrices(0, st.mstate, st.mstinf, st.numsta,
                                          f.N, f.NF, f.B, f.K_code,
                                          selfen, p.cbreak, 4)
    ham0, ham, _ = qcdf_opt.build_matrices(1, st.mstate, st.mstinf, st.numsta,
                                           f.N, f.NF, f.B, f.K_code,
                                           selfen, p.cbreak, 4)
    sel = np.ix_(idx, idx)
    return dict(numsta_pre=st.numsta, norm=hnorm[sel],
                ham=ham[sel], ham0=ham0[0][sel], fortran=f)


def test_state_generation_agrees(matched, fortran_K21):
    """Both codes generate the same pre-weeding basis, and every Fortran
    post-weeding state is found in it."""
    assert matched["numsta_pre"] == fortran_K21.numsta_pre == 502


def test_norm_matrix_is_bit_identical(matched):
    """The color-singlet overlap engine, to the last bit."""
    diff = np.abs(matched["norm"] - matched["fortran"].norm).max()
    assert diff == 0.0, f"norm matrix differs by {diff:.3e}"


def test_interacting_hamiltonian_matches(matched):
    """HNU = Z^T H Z in the Fortran's own basis: the seven four-point vertices."""
    f = matched["fortran"]
    Z = f.Z
    py = Z.T @ matched["ham"] @ Z
    from dlcq.read_fortran import read_ham
    from conftest import find_run

    out = find_run("K21_B1")
    ham_file = out.with_name("qcdf.ham") if out else None
    if ham_file is None or not ham_file.exists():
        pytest.skip("no qcdf.ham to compare against")
    _, _, fort = read_ham(ham_file)
    assert np.abs(py - fort).max() < MATRIX_ATOL


def test_free_hamiltonian_matches(matched):
    """HNU0: the mass term and both self-energy contributions."""
    from dlcq.read_fortran import read_ham
    from conftest import find_run

    out = find_run("K21_B1")
    ham_file = out.with_name("qcdf.ham") if out else None
    if ham_file is None or not ham_file.exists():
        pytest.skip("no qcdf.ham to compare against")
    _, fort0, _ = read_ham(ham_file)
    Z = matched["fortran"].Z
    py0 = np.diag(Z.T @ matched["ham0"] @ Z)
    assert np.abs(py0 - fort0[0]).max() < MATRIX_ATOL


def test_free_hamiltonian_is_proportional_to_the_norm(matched):
    """H0 = D*N for diagonal D -- the free energy is a per-state scalar.

    This is what makes the Fortran's diagonal-only assembly *look* harmless,
    and understanding it is what explains why it is not.
    """
    R = np.linalg.solve(matched["norm"], matched["ham0"])
    off = R - np.diag(np.diag(R))
    assert np.abs(off).max() < 1e-8


def test_exact_assembly_is_basis_independent(matched):
    """The full Z^T H0 Z assembly gives the same spectrum for any valid Z.

    Contrast with test_diagonal_assembly_is_basis_dependent below.  Together
    they are the evidence for choosing assembly="exact" as the reference.
    """
    from scipy.linalg import eigh

    f = matched["fortran"]
    lam2 = f.rlamb ** 2
    w_n, z_n = eigh(matched["norm"])
    Z_lapack = z_n / np.sqrt(w_n)[None, :]

    def spectrum(Z):
        M = lam2 * (Z.T @ matched["ham"] @ Z) + (1 - lam2) * (Z.T @ matched["ham0"] @ Z)
        return f.K_code * eigh(M, eigvals_only=True) / 2.0

    assert np.abs(spectrum(f.Z) - spectrum(Z_lapack)).max() < 1e-9


def test_diagonal_assembly_is_basis_dependent(matched):
    """Regression guard on the documented artifact.

    If a future change made this pass, the basis-dependence would have been
    fixed -- update docs/basis-dependence.md and the tolerances rather than
    deleting the test.
    """
    from scipy.linalg import eigh

    f = matched["fortran"]
    lam2 = f.rlamb ** 2
    w_n, z_n = eigh(matched["norm"])
    Z_lapack = z_n / np.sqrt(w_n)[None, :]

    def spectrum(Z):
        M = lam2 * (Z.T @ matched["ham"] @ Z)
        M[np.diag_indices(M.shape[0])] += (1 - lam2) * np.diag(Z.T @ matched["ham0"] @ Z)
        return f.K_code * eigh(M, eigvals_only=True) / 2.0

    spread = np.abs(spectrum(f.Z) - spectrum(Z_lapack)).max()
    assert spread > 1e-3, (
        "the diagonal-only assembly now looks basis independent; "
        "see docs/basis-dependence.md"
    )


def test_fortran_scheme_reproduces_fortran_spectrum(matched):
    """With the Fortran's own Z and its diagonal-only rule, agreement is exact."""
    from scipy.linalg import eigh

    f = matched["fortran"]
    lam2 = f.rlamb ** 2
    Z = f.Z
    M = lam2 * (Z.T @ matched["ham"] @ Z)
    M[np.diag_indices(M.shape[0])] += (1 - lam2) * np.diag(Z.T @ matched["ham0"] @ Z)
    ev = f.K_code * eigh(M, eigvals_only=True) / 2.0
    assert np.abs(ev - f.eigenvalues).max() < 1e-8


def _basis_spread(matched, n):
    """How far the Fortran's own scheme moves under an equally valid basis.

    Measured in-repo rather than hard-coded: run qcdf.f's diagonal-only
    assembly twice, once with the Z it printed and once with LAPACK's, and take
    the relative spread over the lowest ``n`` states.  This is the precision
    the original algorithm actually defines, so it is the only defensible
    yardstick for an independent reimplementation.
    """
    from scipy.linalg import eigh

    f = matched["fortran"]
    lam2 = f.rlamb ** 2
    w_n, z_n = eigh(matched["norm"])
    Z_lapack = z_n / np.sqrt(w_n)[None, :]

    def spectrum(Z):
        M = lam2 * (Z.T @ matched["ham"] @ Z)
        M[np.diag_indices(M.shape[0])] += (1 - lam2) * np.diag(Z.T @ matched["ham0"] @ Z)
        return f.K_code * eigh(M, eigvals_only=True) / 2.0

    a, b = spectrum(f.Z)[:n], spectrum(Z_lapack)[:n]
    return float(np.max(np.abs(a - b) / np.maximum(np.abs(a), 1e-30)))


@pytest.mark.slow
def test_python_spectrum_within_algorithmic_floor(python_K21, fortran_K21, matched):
    """End-to-end: an independent Python run vs the Fortran, low-lying states.

    The tolerance is *derived*, not assumed: Python must agree with the Fortran
    at least as well as the Fortran's own algorithm agrees with itself when
    handed a different but equally valid orthonormalization.  A hard-coded
    constant here would be testing floating-point luck.
    """
    n = min(10, python_K21.n_eigenvalues, fortran_K21.n_eigenvalues)
    floor = _basis_spread(matched, n)
    tol = max(3.0 * floor, 1e-6)          # 3x headroom over the measured floor

    a = fortran_K21.eigenvalues[:n]
    b = python_K21.eigenvalues[:n]
    rel = np.abs(a - b) / np.maximum(np.abs(a), 1e-30)

    assert rel.max() < tol, (
        f"Python differs from Fortran by {rel.max():.2e} on the lowest {n} "
        f"states, exceeding 3x the algorithm's own basis-dependence floor "
        f"({floor:.2e}). That is a real regression, not the known artifact."
    )


@pytest.mark.slow
def test_measured_floor_is_in_the_documented_range(matched):
    """Pin the artifact's size so a change in its magnitude is noticed.

    docs/basis-dependence.md quotes ~1e-4 on low-lying states, corroborated by
    the -O0 vs -O2 recompilation experiment.
    """
    floor = _basis_spread(matched, 10)
    assert 1e-6 < floor < 1e-2, f"basis-dependence floor moved to {floor:.2e}"


@pytest.mark.slow
def test_python_sum_rules(python_K21):
    """The Python path must satisfy the same exact sum rules."""
    from dlcq.observables import momentum_sum_rule, number_sum_rule

    assert momentum_sum_rule(python_K21, 0) == pytest.approx(1.0, abs=1e-10)
    assert number_sum_rule(python_K21, 0) == pytest.approx(3.0, abs=1e-10)


# ── removing the basis dependence ─────────────────────────────────────────

def test_blockwise_Z_makes_the_free_part_diagonal(matched):
    """Block-wise orthonormalization makes qcdf.f's assembly well-posed.

    ``H0 = D N`` with D constant on each block of mutually overlapping states,
    so a Z built inside blocks keeps every column in one block and ``Z^T H0 Z``
    comes out diagonal -- which is exactly the condition the Fortran's
    diagonal-only rule silently assumes.
    """
    from dlcq.read_python import orthonormalize_blockwise

    Z, ncol = orthonormalize_blockwise(matched["norm"], matched["norm"].shape[0])
    assert ncol > 0

    gram = Z.T @ matched["norm"] @ Z
    assert np.abs(gram - np.eye(ncol)).max() < 1e-9

    free = Z.T @ matched["ham0"] @ Z
    off = free - np.diag(np.diag(free))
    assert np.abs(off).max() < 1e-9, (
        f"Z^T H0 Z off-diagonal {np.abs(off).max():.2e}; a global Z gives ~0.8"
    )


def test_blockwise_Z_removes_the_assembly_ambiguity(matched):
    """With a block-wise Z the two assemblies agree, so the answer is unique."""
    from scipy.linalg import eigh

    from dlcq.read_python import orthonormalize_blockwise

    f = matched["fortran"]
    lam2 = f.rlamb ** 2
    Z, _ = orthonormalize_blockwise(matched["norm"], matched["norm"].shape[0])

    full = lam2 * (Z.T @ matched["ham"] @ Z) + (1 - lam2) * (Z.T @ matched["ham0"] @ Z)
    diag = lam2 * (Z.T @ matched["ham"] @ Z)
    diag[np.diag_indices(diag.shape[0])] += (1 - lam2) * np.diag(Z.T @ matched["ham0"] @ Z)

    a = f.K_code * eigh(full, eigvals_only=True) / 2.0
    b = f.K_code * eigh(diag, eigvals_only=True) / 2.0
    assert np.abs(a - b).max() < 1e-9, (
        "block-wise Z should make the diagonal-only and full assemblies identical"
    )


# ── blocked weeding: the same algorithm, per norm block ───────────────────

@pytest.mark.slow
def test_blocked_weeding_preserves_the_spectrum():
    """Weeding per norm block gives the same physics as weeding globally.

    The norm is block-diagonal in momentum configuration (verified by explicit
    colour enumeration, `tools/colour_norm.py`), so a state can only be linearly
    dependent on states in its own block and both weeding stages decompose.
    Running them per block instead of on all states at once was worth ~500-900x
    on the weeding step, which had been 76% of a run's wall time.

    The two are not bit-identical.  WEEDR2 diagonalizes the norm, and when
    several blocks have degenerate zero eigenvalues `eigh` can return null
    vectors that mix them; the global pass then drops a different member of the
    same null space, and can drop one or two states the blocked pass keeps
    (191 vs 193 at 2K=21).  The surviving span -- and so the spectrum -- is what
    must agree, and it does, far inside the ~1e-4 floor that
    docs/basis-dependence.md establishes for this algorithm.
    """
    import numpy as np

    import dlcq.read_python as rp
    from dlcq.observables import physical_indices
    from dlcq.units import mg_to_lambda

    lam = float(mg_to_lambda(1.6))
    original = rp.weed_fortran
    spectra = {}
    try:
        for blocked in (False, True):
            rp.weed_fortran = (
                lambda *a, _b=blocked, **k: original(*a, **{**k, "blocked": _b}))
            r = rp.run_python(N=3, NF=1, B=1, K_code=21, rlamb=lam, ncpus=4,
                              policy="fortran")
            spectra[blocked] = (r.numsta_post,
                                np.sort(r.eigenvalues[physical_indices(r)]))
    finally:
        rp.weed_fortran = original

    (n_global, e_global), (n_blocked, e_blocked) = spectra[False], spectra[True]
    assert n_blocked >= n_global, "blocked weeding should not drop more states"

    # The ground state is the number every figure and Table I entry reads.
    assert abs(e_blocked[0] - e_global[0]) / e_global[0] < 1e-6, (
        f"ground state moved: {e_global[0]:.9f} -> {e_blocked[0]:.9f}")

    # And the low-lying spectrum, well inside the algorithm's own 1e-4 floor.
    for j in range(6):
        assert abs(e_blocked[j] - e_global[j]) / e_global[j] < 1e-5, (
            f"level {j}: {e_global[j]:.9f} -> {e_blocked[j]:.9f}")


def test_the_norm_really_is_block_diagonal():
    """The premise blocked weeding rests on, checked directly.

    If this ever fails, blocked weeding is unsound and `blocked=False` is the
    correct setting -- the speedup is not worth a wrong basis.
    """
    import numpy as np
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components

    from dlcq.figures import paper_lambda
    from dlcq.providers import PythonProvider

    r = PythonProvider(ncpus=4).get(3, 1, 1, 21, paper_lambda(1.6))
    N = r.norm
    n = r.numsta_post
    lengths = np.asarray(r.state_len[:n])

    # no norm element couples different Fock sectors ...
    for a in sorted(set(lengths.tolist())):
        for b in sorted(set(lengths.tolist())):
            if a == b:
                continue
            block = N[np.ix_(lengths == a, lengths == b)]
            if block.size:
                assert np.abs(block).max() == 0.0, f"{a}q-{b}q coupling"

    # ... and the blocks are many and small, which is what makes it pay
    ncomp, _ = connected_components(
        csr_matrix((np.abs(N[:n, :n]) > 1e-9).astype(np.int8)), directed=False)
    assert ncomp > n / 4, f"only {ncomp} blocks among {n} states"


def test_config_labels_match_the_matrix_derived_blocks():
    """Block labels from the state list agree with those from the norm matrix.

    ``weed_fortran`` used to recover its blocks by thresholding the assembled
    norm -- ``csr_matrix(np.abs(A) > 1e-9)`` -- which costs two more dense n x n
    temporaries on top of a matrix that is already 8.6 GB at 2K=37.  That is
    what ran a 15 GB machine out of memory.  ``config_block_labels`` gets the
    same blocks in O(n L) from the parton content, which is what makes them
    block-diagonal in the first place.

    The two are not identical: configuration groups are *coarser*, because a
    few groups split further when their states happen to be colour-orthogonal.
    Coarser is safe -- weeding a superset of a block is still exact -- and what
    must hold is that no component straddles two configuration groups, and that
    both retain the same basis.
    """
    import numpy as np

    import qcdf as base
    import qcdf_opt as opt
    from dlcq.read_python import (_labels_from_matrix, config_block_labels,
                                  weed_fortran)

    for B, K in ((1, 21), (0, 24)):
        p = base.Params()
        p.N, p.NF, p.B, p.K, p.rlamb = 3, 1, B, K, 0.3325
        p.cutoff, p.LPN = -1.0, 0
        p.iflv[0] = 3 * B
        st = base.StateData()
        opt.qcdsta_fast(p, st, base.PermTables(), base.FlavorTables())
        n = st.numsta
        _, _, hn = opt.build_matrices(0, st.mstate, st.mstinf, n,
                                      3, 1, B, K, opt.compute_selfen(3),
                                      p.cbreak, 4)
        A = hn[:n, :n]
        mi = st.mstinf[:n].copy()

        from_matrix = _labels_from_matrix(A)
        from_states = config_block_labels(st.mstate, mi, n)

        # Every matrix-derived component must sit inside one config group.
        for c in np.unique(from_matrix):
            assert len(set(from_states[from_matrix == c].tolist())) == 1, (
                f"component {c} straddles two configuration groups")

        # And the weeding outcome must be the same basis size -- which is the
        # rank of the Gram matrix, so neither labelling loses a direction.
        rank = int(np.sum(np.abs(np.linalg.eigvalsh(A)) > 1e-8))
        _, _, n_matrix = weed_fortran(A, mi, n, labels=from_matrix)
        _, _, n_states = weed_fortran(A, mi, n, labels=from_states)
        assert n_matrix == n_states == rank, (
            f"B={B} 2K={K}: matrix {n_matrix}, states {n_states}, rank {rank}")


def test_weeding_removes_exactly_the_null_space():
    """Weeding is redundancy removal, not truncation.

    The Fock basis is non-orthogonal -- different colour contractions of the
    same partons are linear combinations of each other -- so the Gram matrix is
    singular and ``Z = N^-1/2`` does not exist until the dependent states go.
    What must hold is that the retained count equals the *rank*: no independent
    direction is discarded, so the span, and with it the physics, is unchanged.

    The threshold has enormous margin: the smallest retained norm eigenvalue is
    6.0 and the largest discarded is ~1e-12.
    """
    import numpy as np

    import qcdf as base
    import qcdf_opt as opt
    from dlcq.read_python import config_block_labels, weed_fortran

    p = base.Params()
    p.N, p.NF, p.B, p.K, p.rlamb = 3, 1, 1, 21, 0.3325
    p.cutoff, p.LPN = -1.0, 0
    p.iflv[0] = 3
    st = base.StateData()
    opt.qcdsta_fast(p, st, base.PermTables(), base.FlavorTables())
    n = st.numsta
    _, _, hn = opt.build_matrices(0, st.mstate, st.mstinf, n, 3, 1, 1, 21,
                                  opt.compute_selfen(3), p.cbreak, 4)
    A = hn[:n, :n]
    w = np.sort(np.abs(np.linalg.eigvalsh(A)))
    rank = int(np.sum(w > 1e-8))

    labels = config_block_labels(st.mstate, st.mstinf[:n], n)
    _, _, kept = weed_fortran(A, st.mstinf[:n].copy(), n, labels=labels)
    assert kept == rank, f"kept {kept} but rank is {rank}"

    # The gap the eps=1e-4 threshold sits in, spanning ~12 orders of magnitude.
    assert w[n - rank - 1] < 1e-9 < 1.0 < w[n - rank]
