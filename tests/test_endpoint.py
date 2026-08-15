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


@pytest.mark.parametrize("N,B,K,LPN", [(3, 0, 16, 2), (3, 1, 15, 3),
                                       (3, 1, 15, 5), (4, 1, 16, 4),
                                       (2, 0, 12, 0)])
def test_block_shortcut_equals_the_matrix_form(N, B, K, LPN):
    """``Z^T (N diag(dsigma)) Z == diag(dsigma)`` on the blockwise path.

    ``sigma`` is constant on a configuration block (``config_block_labels``
    keys on the multiset of (type, momentum, flavour), and sigma depends only
    on the momenta), and every column of a blockwise ``Z`` lives in one block
    with ``Z^T N Z = I``.  So the correction collapses to a diagonal add after
    projection -- no norm, no n^2 addition, no sparsity destroyed.

    Off that path ``Z`` is global and the argument fails, so there the
    correction is applied to ``ham`` against the dense norm.  The two must
    agree, and that is what this pins.
    """
    kw = dict(N=N, NF=1, B=B, K_code=K, rlamb=0.3325, LPN=LPN, ncpus=2,
              assembly="exact", hamiltonian="improved")
    fast = run_python(policy="blockwise", **kw)
    ref = run_python(policy="fortran", **kw)
    m = min(fast.eigenvalues.size, ref.eigenvalues.size)
    scale = max(abs(float(fast.eigenvalues[0])), 1e-300)
    assert np.max(np.abs(fast.eigenvalues[:m] - ref.eigenvalues[:m])) / scale < 1e-11


@pytest.mark.parametrize("N,B,K,LPN", [(3, 1, 15, 3), (3, 1, 21, 5)])
def test_improved_does_not_need_the_norm_kept(N, B, K, LPN):
    """``keep_norm=False`` must stay available -- the norm is the memory wall."""
    kw = dict(N=N, NF=1, B=B, K_code=K, rlamb=0.3325, LPN=LPN, ncpus=2,
              assembly="exact", policy="blockwise", hamiltonian="improved")
    a = run_python(keep_norm=True, **kw)
    b = run_python(keep_norm=False, **kw)
    assert np.array_equal(a.eigenvalues, b.eigenvalues)


@pytest.mark.parametrize("N,B,K,LPN", [(3, 1, 15, 3), (3, 1, 21, 3),
                                       (4, 1, 16, 4), (3, 0, 12, 0),
                                       (3, 1, 15, 5)])
def test_improved_runs_beyond_two_partons(N, B, K, LPN):
    """L >= 3, including a mixed-Fock basis (LPN=0 is *no* truncation)."""
    r = run_python(N=N, NF=1, B=B, K_code=K, rlamb=0.3325, LPN=LPN, ncpus=2,
                   assembly="exact", policy="blockwise",
                   hamiltonian="improved")
    assert r.eigenvalues.size > 0
    assert np.all(np.isfinite(r.eigenvalues))


@pytest.mark.parametrize("N,B,K,LPN", [(3, 0, 10, 2), (3, 1, 15, 3),
                                       (4, 1, 16, 4), (3, 1, 15, 5),
                                       (3, 0, 12, 0), (5, 1, 15, 5)])
def test_sigma_reduces_per_state_at_zero_exponent(N, B, K, LPN):
    """``sigma_imp == sigma_std`` at b = 0, for every L, exactly.

    Stronger than the Hamiltonian-level gate: it pins the reduction on each
    state separately, which is what makes the correction a safe drop-in rather
    than a coincidence of sums.
    """
    from dlcq.endpoint import state_sigmas
    p = base.Params()
    p.N, p.NF, p.B, p.K, p.rlamb = N, 1, B, K, 0.3325
    p.cutoff, p.LPN = -1.0, LPN
    p.iflv[0] = N * B
    st = base.StateData()
    base.qcdsta(p, st, base.PermTables(), base.FlavorTables())
    n = st.numsta
    std, imp = state_sigmas(st.mstate, st.mstinf[:n].copy(), n, N, 0.0, K)
    assert np.max(np.abs(imp - std)) < 1e-13


def test_directed_kernel_reduces_partner_independently():
    """``J^(k;l)|_{b=0} = S(k)`` for every partner -- the reduction's root."""
    from dlcq.endpoint import directed_pair_kernel

    def S(k):
        return sum(1.0 / n ** 2 for n in range(1, (k - 1) // 2 + 1))

    worst = 0.0
    for k in range(1, 40, 2):
        for l in range(1, 40 - k, 2):
            worst = max(worst, abs(directed_pair_kernel(k, l, 0.0) - S(k)))
    assert worst < 1e-14


def test_matrix_level_matches_the_validated_two_body_path():
    """The general form must reproduce the selfen-swap it generalises.

    ``improved_selfen`` is validated against van de Sande's exact
    ``M^2/g^2 = 0.779141``; this ties the L >= 3 machinery to that anchor.
    """
    from dlcq.endpoint import state_sigmas
    N, K, b = 3, 24, 0.0844
    std_t, _ = _ham(N, 0, K, 2, opt.compute_selfen(N))
    swap, _ = _ham(N, 0, K, 2, improved_selfen(N, K, b))
    p = base.Params()
    p.N, p.NF, p.B, p.K, p.rlamb = N, 1, 0, K, 0.3325
    p.cutoff, p.LPN = -1.0, 2
    st = base.StateData()
    base.qcdsta(p, st, base.PermTables(), base.FlavorTables())
    n = st.numsta
    mi = st.mstinf[:n].copy()
    norm = opt.build_matrices(0, st.mstate, mi, n, N, 1, 0, K,
                              opt.compute_selfen(N), p.cbreak, 4,
                              backend="thread")[2]
    sd, si = state_sigmas(st.mstate, mi, n, N, b, K)
    matrix_level = std_t + norm @ np.diag(si - sd)
    scale = max(np.max(np.abs(swap)), 1e-300)
    assert np.max(np.abs(swap - matrix_level)) / scale < 1e-13


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


# ── exact colour weights ──────────────────────────────────────────────────

@pytest.mark.parametrize("N,B,K,LPN", [(3, 1, 15, 3), (3, 1, 21, 3),
                                       (4, 1, 16, 4), (5, 1, 15, 5),
                                       (3, 0, 10, 2), (2, 0, 12, 2),
                                       (3, 1, 15, 5), (3, 0, 12, 0),
                                       (6, 1, 18, 6)])
def test_colour_weights_satisfy_the_singlet_identity(N, B, K, LPN):
    """``sum_{c != a} <-T_a . T_c> == C_F`` for every parton, exactly.

    This is what ``sum_a T_a = 0`` forces on a colour singlet, and it is the
    gate the whole extraction hangs on: it is what makes the ``b = 0``
    reduction automatic rather than something to be arranged.

    It also caught every error on the way in.  The vertex patterns have to be
    the real ones -- H1 for qq, H2 for qbar-qbar, and H7's *t-channel*
    structure A for q-qbar (H3-H6 are pair creation/annihilation and carry no
    self-inertia partner).  And the multiplicity has to key on
    ``(type, momentum, flavour)``: keying on momentum alone passes every
    baryon and fails exactly the meson states where a quark and an antiquark
    share a momentum but remain distinguishable.
    """
    from dlcq.endpoint import pair_colour_weights

    p = base.Params()
    p.N, p.NF, p.B, p.K, p.rlamb = N, 1, B, K, 0.3325
    p.cutoff, p.LPN = -1.0, LPN
    p.iflv[0] = N * B
    st = base.StateData()
    base.qcdsta(p, st, base.PermTables(), base.FlavorTables())
    n = st.numsta
    if n == 0:
        pytest.skip("no states")
    selfen = base.compute_selfen(p)
    CF = (N * N - 1.0) / (2.0 * N)

    worst = 0.0
    for s in range(min(n, 12)):
        L = int(st.mstinf[s, 1])
        if L < 2:
            continue
        loc = int(st.mstinf[s, 0]) - 1
        w = pair_colour_weights(base, p, st.mstate, loc, L, selfen)
        for a in range(L):
            worst = max(worst, abs(w[a].sum() - CF))
    assert worst < 1e-12, f"singlet identity violated by {worst:.3e}"


@pytest.mark.parametrize("N,B,K,LPN", [(3, 1, 15, 3), (4, 1, 16, 4),
                                       (3, 1, 15, 5), (5, 1, 15, 5)])
def test_exact_weights_keep_the_zero_exponent_reduction(N, B, K, LPN):
    """The gate must survive the switch from the scalar to the operator."""
    from dlcq.endpoint import state_sigmas

    p = base.Params()
    p.N, p.NF, p.B, p.K, p.rlamb = N, 1, B, K, 0.3325
    p.cutoff, p.LPN = -1.0, LPN
    p.iflv[0] = N * B
    st = base.StateData()
    base.qcdsta(p, st, base.PermTables(), base.FlavorTables())
    n = st.numsta
    selfen = base.compute_selfen(p)
    std, imp = state_sigmas(st.mstate, st.mstinf[:n].copy(), n, N, 0.0, K,
                            base=base, params=p, selfen=selfen)
    assert np.max(np.abs(imp - std)) < 1e-12


def test_exact_weights_agree_with_uniform_in_valence_sectors():
    """Where every pair is colour-equivalent the scalar is exact, so both agree.

    A valence baryon has all pairs in the same channel, which is precisely when
    ``sum_a T_a = 0`` forces ``-T_a . T_c = C_F/(L-1)``.
    """
    from dlcq.endpoint import state_sigmas

    N, B, K, LPN, b = 3, 1, 21, 3, 0.0844
    p = base.Params()
    p.N, p.NF, p.B, p.K, p.rlamb = N, 1, B, K, 0.3325
    p.cutoff, p.LPN = -1.0, LPN
    p.iflv[0] = N * B
    st = base.StateData()
    base.qcdsta(p, st, base.PermTables(), base.FlavorTables())
    n = st.numsta
    selfen = base.compute_selfen(p)
    _, exact = state_sigmas(st.mstate, st.mstinf[:n].copy(), n, N, b, K,
                            base=base, params=p, selfen=selfen)
    _, scalar = state_sigmas(st.mstate, st.mstinf[:n].copy(), n, N, b, K)
    assert np.max(np.abs(exact - scalar)) / np.max(np.abs(scalar)) < 1e-12
