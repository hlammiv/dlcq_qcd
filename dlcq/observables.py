"""Physical observables extracted from a :class:`DLCQResult`.

Structure functions, sum rules, and the Richardson extrapolation that turns
finite-K spectra into the continuum numbers of Table I.

``structure_function`` is a **reconstruction of a lost program.**  The x-space
conversion in the original work was done by a separate code, ``wf``/``wfbig``,
which does not survive: ``qcdf.f`` emits only eigenvectors, Fock content and the
basis-change matrix, and refers to ``wf`` in comments.  Every structure function
in the paper passed through it.

That also means this is the one piece of the pipeline a solver-vs-solver
comparison cannot check -- the Python port reproduces the Fortran's matrix
elements to 1e-14, so an error here would be shared by both.  It is instead
validated against the thesis's own definition, the two sum rules (machine
precision, both solvers), a brute-force colour enumeration of the norm it uses,
and the published curves.  See docs/baryon-higher-fock.md section 1.3.

The structure function follows Eq. (12) of the paper:

    q(x) = K_paper <phi| b_k^dag b_k |phi>,     x = k / K_code

Basis states are number eigenstates, so ``n_k |s> = n_k(s) |s>`` and the
expectation value in the *non-orthogonal* Fock basis collapses to

    <phi| n_k |phi> = sum_s n_k(s) c_s (N c)_s        with  c^T N c = 1

which is why the weight is ``c_s * (N c)_s`` rather than ``c_s^2``.  Using
``c_s^2`` would silently break every sum rule.

The measure is ``dx = 2 / K_code = 1 / K_paper`` because momenta step by 2 in
code units (odd integers only).
"""

from __future__ import annotations

import numpy as np

from .dataset import DLCQResult
from .units import endpoint_exponent

__all__ = [
    "structure_function",
    "momentum_sum_rule",
    "number_sum_rule",
    "valence_parton_count",
    "richardson_extrapolate",
    "richardson_evaluate",
    "richardson_holdout",
    "richardson_stability",
    "wynn_epsilon",
    "levin_u",
    "monotone_bracket",
    "richardson_ensemble",
    "richardson_budget",
    "thooft_valence_limit",
    "spurious_zero_modes",
    "physical_indices",
    "require_physical_index",
    "color_slots_required",
    "fortran_overflow_risk",
    "FORTRAN_COLOR_SLOTS",
]


def color_slots_required(max_partons):
    """Color-index slots a matrix element needs: ``2 * max_partons + 4``.

    Bra and ket each contribute their partons, plus up to four operators from a
    four-point vertex.  ``qcdf.f`` dimensions ``IDELT(12552, 25)`` and sets
    ``MXLNG = 25``, so any run reaching 11 partons overflows that array.
    ``qcdf_opt.py`` uses ``MXLNG = 2*MXP + 4 = 54`` and is unaffected.
    """
    return 2 * int(max_partons) + 4


FORTRAN_COLOR_SLOTS = 25


def fortran_overflow_risk(result):
    """Whether this configuration exceeds ``qcdf.f``'s 25-slot color array."""
    if result.state_len is None or len(result.state_len) == 0:
        return False
    return color_slots_required(int(np.max(result.state_len))) > FORTRAN_COLOR_SLOTS


def spurious_zero_modes(result, tol=1e-12):
    """Indices of unphysical eigenstates -- non-positive M^2 at non-zero mass.

    ``M^2 >= 0`` for any physical state, and the paper states that an exactly
    massless state occurs *only* in the chiral limit: "when m/g = 0 identically,
    the lightest state for any N or B is exactly zero, independently of K".  So
    at m/g > 0 any eigenvalue at or below zero is an artifact.

    These are not hypothetical, and they sit at the *bottom* of the spectrum --
    precisely where Figs. 7, 8 and Table I read off the lightest state.  Their
    cause is an array-bounds overflow in ``qcdf.f``: a matrix element between
    two L-parton states needs ``2L + 4`` color-index slots, but ``IDELT`` is
    dimensioned with only 25, so runs reaching L >= 11 corrupt the color
    contraction.  Observed:

    ===============  ====  ======  ==========================================
    run              L     slots   symptom
    ===============  ====  ======  ==========================================
    2K=21, B=1        9     22     clean
    2K=10, B=0        6     16     clean
    2K=24, B=0       12     28     one decoupled 12-parton state at M^2 = 0
    2K=25, B=1       11     26     four negative M^2, all 11-parton dominated
    ===============  ====  ======  ==========================================

    See ``docs/fortran-color-overflow.md``.

    The Python solver has its **own** version of this defect, and the claim
    that it is unaffected -- previously made here -- was wrong.  ``MXTRM =
    12552`` bounds the colour term tables in both ``qcdf_kernels.clfact_nb``
    and the interpreted ``qcdf_opt.clfact``; overrunning it used to return
    ``0.0`` from one guard and truncate to a partial sum at the other.  Since
    both paths share the bound they agreed bit-for-bit on the wrong answer, so
    the ``array_equal`` backend tests could not see it.  Measured at N=4, B=1,
    2K=20 against ``tools/colour_norm.norm_bruteforce``: a state's own norm came
    back 150048 against a true 331776, and the assembled matrix had an
    eigenvalue of -1.8e5.

    Unlike the Fortran's, it is an **N >= 4 baryon** defect rather than a
    parton-count one -- N=3 is clean at L=15 (2K=41) while N=4 breaks at L=12.
    Both guards now raise ``OverflowError``, so it can no longer reach a
    spectrum silently.  ``tests/test_colour_overflow.py`` holds the invariant
    that catches it: the norm is a Gram matrix, so it must be positive
    semidefinite.

    No published result is affected either way: Table I and Figs. 7-8 run under
    ``sweep_lpn``, whose cap exists precisely to keep ``2L + 4`` inside the
    array, and Figs. 4-6 are N=3.

    Returns an array of indices into ``result.eigenvalues``.
    """
    if result.eigenvalues is None or result.eigenvalues.size == 0:
        return np.array([], dtype=int)
    # rlamb == 1 is the chiral limit, where a zero mode is genuine.
    if result.rlamb >= 1.0 - 1e-12:
        return np.array([], dtype=int)
    return np.flatnonzero(result.eigenvalues <= tol)


def physical_indices(result, tol=1e-12):
    """Eigenstate indices with the spurious zero modes removed.

    Use this wherever the paper says "the first three states"; indexing
    ``eigenvalues`` directly can silently pick up a decoupled mode.
    """
    bad = set(spurious_zero_modes(result, tol).tolist())
    n = result.n_eigenvalues
    return np.array([i for i in range(n) if i not in bad], dtype=int)


def require_physical_index(result, level: int, tol: float = 1e-12) -> int:
    """Index of the ``level``-th physical eigenstate, failing loudly if absent.

    :func:`physical_indices` cannot tell "this run has few levels" from "this
    run *computed* few levels", so indexing it directly is safe only near the
    bottom of the spectrum.  Fig. 5 wants the 11th physical state, which a
    40-level solve supplies only if fewer than 29 spurious modes sit below it.

    Callers that index deep should come through here.  The failure mode this
    replaces is worse than a crash: ``figures._wavefunction_panel`` fell back to
    indexing the *raw* array when ``physical_indices`` came up short, which
    draws a plausible-looking panel of the wrong state.
    """
    phys = physical_indices(result, tol)
    if level < phys.size:
        return int(phys[level])
    if getattr(result, "spectrum_is_truncated", False):
        raise IndexError(
            f"physical level {level} unavailable: this run computed the lowest "
            f"{result.n_eigenvalues} of {result.n_orth} eigenvalues, of which "
            f"{phys.size} are physical. Raise nev and rerun.")
    raise IndexError(
        f"physical level {level} unavailable: this run has only {phys.size} "
        f"physical states among {result.n_eigenvalues} eigenvalues.")


def valence_parton_count(N: int, B: int) -> int:
    """Partons in the minimal (valence) Fock sector: 2 for a meson, N*|B| otherwise."""
    return 2 if B == 0 else N * abs(B)


def structure_function(result: DLCQResult, state_idx: int = 0,
                       nparton: int | None = None):
    """Quark and antiquark structure functions for one eigenstate.

    Parameters
    ----------
    result
        A parsed run.  Needs ``c_orig``, ``norm``, and the basis arrays.
    state_idx
        Which eigenstate, 0 = ground state.
    nparton
        Restrict to the Fock sector with exactly this many partons (the
        decomposition plotted in Figs. 4-6).  ``None`` sums all sectors.

    Returns
    -------
    x, q_quark, q_antiquark
        ``x`` is the odd-momentum grid ``k / K_code``.  The measure that goes
        with it is ``dx = 2 / K_code``.
    """
    if result.c_orig is None or result.norm is None:
        raise ValueError("result lacks c_orig/norm; parse with with_matrices=True")

    c = result.require_eigenvector(state_idx)
    Nc = result.norm @ c
    weights = c * Nc                      # sums to 1 over all states

    K = result.K_code
    K_paper = result.K_paper

    q = np.zeros(K + 1)
    qbar = np.zeros(K + 1)

    for s in range(result.numsta_post):
        if nparton is not None and result.state_len[s] != nparton:
            continue
        w = weights[s]
        if w == 0.0:
            continue
        L = result.state_len[s]
        for p in range(L):
            k = result.state_moms[s, p]
            if result.state_types[s, p] == 1:
                q[k] += w
            else:
                qbar[k] += w

    ks = np.arange(1, K, 2)
    x = ks / float(K)
    return x, K_paper * q[ks], K_paper * qbar[ks]


def momentum_sum_rule(result: DLCQResult, state_idx: int = 0) -> float:
    """``integral x [q(x) + qbar(x)] dx``, which must equal 1 exactly.

    Every basis state carries total momentum ``K_code``, so this holds
    state-by-state and is independent of the basis conditioning -- the
    strongest available check that the structure function is built correctly.
    """
    x, q, qbar = structure_function(result, state_idx)
    dx = 2.0 / result.K_code
    return float(np.sum(x * (q + qbar)) * dx)


def number_sum_rule(result: DLCQResult, state_idx: int = 0) -> float:
    """``integral [q(x) - qbar(x)] dx``, which must equal ``N * B``."""
    x, q, qbar = structure_function(result, state_idx)
    dx = 2.0 / result.K_code
    return float(np.sum(q - qbar) * dx)


# ──────────────────────────────────────────────────────────────────────────
# Continuum extrapolation
# ──────────────────────────────────────────────────────────────────────────

def _richardson_design(Kp, a, n_terms, basis="paper"):
    """Columns of the Eq. (27) fit at the paper's ``K``, plus a label for each.

    ``basis="paper"`` uses the exponents literally: ``1, K^-1, K^-(1+a), ...``.
    That is what the paper wrote and it is the default, but it degenerates in
    the chiral limit, where ``a -> 0`` sends ``K^-(1+a)`` onto ``K^-1``.  Taken
    literally the design matrix is then singular -- its condition number over
    2K = 25-49 at N = 3 goes 3.2e7 at m/g = 1.6, 5.6e7 at 0.1, 1.1e8 at 0.05,
    and **3.1e19 at 0**, past what float64 can represent.

    This function never evaluates that singular matrix, because it drops
    duplicate exponents first.  But the cure has its own problem: at ``a = 0``
    *both* ``1+a -> 1`` and ``2+a -> 2`` collapse, so a five-term request
    silently becomes a **three-term fit**.  The model dimension therefore jumps
    discontinuously as ``m/g -> 0``: five terms at m/g = 0.05, three at m/g = 0,
    with nothing in between.  The near-singular region just short of the limit
    (cond 1.1e8 at m/g = 0.05) is fitted at full width and is where the
    coefficient blow-up actually bites.

    ``basis="confluent"`` fixes it properly.  Two exponents merging is a
    *confluent* limit -- the same situation as a repeated root of a
    characteristic polynomial -- and the correct limiting basis is the function
    together with its derivative in the exponent::

        d/da K^-(1+a) = -ln(K) K^-(1+a)

    so the merging pair is replaced by the divided difference

        (K^-(1+a) - K^-1) / a   ->   -K^-1 ln K   as a -> 0

    This **spans exactly the same space for every a > 0**, so ``M0`` is
    unchanged (verified identical to 1.1e-14); only the coordinates change.
    Measured conditioning gain against the literal basis: 1.2x at m/g = 1.6,
    14x at 0.1, 29x at 0.05.

    At ``a = 0`` exactly it is *not* the same space -- it is strictly larger,
    because it keeps the two log columns the paper path deletes.  That is the
    intended behaviour: the confluent basis is continuous in ``a`` where the
    paper one is not.

    That matters for the *error estimate* rather than the answer.  The paper
    reports "the magnitude of the last term in the series fit"; in the paper
    basis at weak coupling that term is large because two near-parallel columns
    carry huge cancelling coefficients (8.383 and 8.800 at m/g = 0.1, against
    0.892 and 0.098 at m/g = 1.6 -- see :func:`richardson_stability`).  In the
    confluent basis the columns are not near-parallel, so the last term measures
    a genuine independent contribution.

    Returns ``(A, exponents, labels)``.  ``A`` has shape
    ``(len(Kp), len(exponents))``.  ``exponents`` are the numeric exponents the
    columns derive from -- for ``basis="paper"`` the columns *are* ``K**-e``, so
    they can be fed to :func:`richardson_curve`; for ``"confluent"`` they are
    only the provenance of each column, and the curve must be evaluated with
    :func:`richardson_evaluate` instead.
    """
    Kp = np.asarray(Kp, dtype=float)
    base = [0.0, 1.0, 1.0 + a, 2.0, 2.0 + a, 3.0][: n_terms + 1]

    if basis == "paper":
        seen, uniq = set(), []
        for e in base:
            key = round(e, 12)
            if key not in seen:
                seen.add(key)
                uniq.append(e)
        return (np.vstack([Kp ** (-e) for e in uniq]).T, uniq,
                [f"K^-{e:g}" for e in uniq])

    if basis != "confluent":
        raise ValueError(f"unknown basis {basis!r}; use 'paper' or 'confluent'")

    # Position, not value: at a == 0 the exponents 1+a and 1 are numerically
    # equal, so keying on the value turns *both* columns into log columns and
    # loses the plain K^-1 the pair is supposed to be built around.
    #   base index:  0    1      2      3      4      5
    #   exponent:    0    1    1+a      2    2+a      3
    partner_of = {2: 1, 4: 3}
    cols, labels = [], []
    for i, e in enumerate(base):
        j = partner_of.get(i)
        if j is None:
            cols.append(Kp ** (-e))
            labels.append(f"K^-{e:g}")
            continue
        p = base[j]
        if a < 1e-8:
            cols.append(-np.log(Kp) * Kp ** (-p))
            labels.append(f"K^-{p:g} lnK")
        else:
            cols.append((Kp ** (-e) - Kp ** (-p)) / a)
            labels.append(f"(K^-{e:g} - K^-{p:g})/a")
    return np.vstack(cols).T, base, labels


def richardson_extrapolate(K_codes, masses, mg, N, n_terms=4, return_fit=False,
                           basis="paper"):
    """Extrapolate M(K) to the continuum using Eq. (27).

        M(1/K) = M(0) + c1/K + c2/K^(1+a) + c3/K^2 + c4/K^(2+a) + ...

    The non-analytic exponents come from the endpoint behaviour phi(x) ~ x^a of
    Eq. (26), with ``a`` solved per (N, m/g).

    This replaces the heuristic previously used in ``reproduce_figures.py``
    (take whichever K gives the highest positive lightest eigenvalue), which is
    not what the paper does and cannot reproduce Table I.

    Parameters
    ----------
    K_codes, masses
        Matched arrays. ``K_codes`` in code units (= 2 * paper K); the fit uses
        the paper's K internally, matching the paper's "2K in the range 16-24".
    mg, N
        Needed for the endpoint exponent.
    n_terms
        Number of correction terms. The paper estimates the error in M(0) as
        the magnitude of the last retained term -- that is what Table I's
        parentheses hold, so it is returned rather than a statistical error.

    Returns
    -------
    M0, last_term
        Continuum mass and the paper's error estimate.  With ``return_fit``,
        also the coefficient vector and the exponents used.
    """
    K_codes = np.asarray(K_codes, dtype=float)
    masses = np.asarray(masses, dtype=float)

    order = np.argsort(K_codes)
    K_codes, masses = K_codes[order], masses[order]

    a = endpoint_exponent(mg, N)
    Kp = K_codes / 2.0                       # the paper's K

    # Keep the fit over-determined by at least one point.  ``n_terms`` costs
    # ``n_terms + 1`` columns, so the original ``len(Kp) - 1`` clamp permitted
    # exact interpolation: zero residual, and an extrapolation driven entirely
    # by rounding.
    #
    # This is a floor against nonsense, not a model-selection rule.  How many
    # terms the *data* supports is a question for the caller, which knows what
    # window it is fitting -- see ``figures.N_TERMS``.  Deciding it here would
    # silently rewrite the model underneath :func:`richardson_stability`, whose
    # whole method is refitting small sub-windows at a fixed order.
    n_terms = min(n_terms, max(len(Kp) - 2, 1))
    A, exponents, _labels = _richardson_design(Kp, a, n_terms, basis)
    coeffs, *_ = np.linalg.lstsq(A, masses, rcond=None)

    M0 = float(coeffs[0])
    # "the numbers in parentheses give the magnitude of the last term in the
    # series fit", evaluated at the largest K used.  Taken from the design
    # matrix rather than as Kp**-e so it stays correct for a basis whose last
    # column is not a bare power (K_codes was sorted, so row -1 is the largest
    # K).  For basis="paper" this is identical to the original expression.
    last_term = float(abs(coeffs[-1] * A[-1, -1]))

    if return_fit:
        return M0, last_term, coeffs, exponents
    return M0, last_term


def richardson_evaluate(coeffs, mg, N, K_codes, n_terms=None, basis="paper"):
    """The fitted series evaluated at arbitrary ``K_codes``.

    Rebuilds the same design the fit used, so it works for either basis.  Used
    by :func:`richardson_holdout` to predict points the fit never saw.
    """
    coeffs = np.asarray(coeffs, dtype=float)
    Kp = np.asarray(K_codes, dtype=float) / 2.0
    if n_terms is None:
        n_terms = len(coeffs) - 1
    A, _, _ = _richardson_design(Kp, endpoint_exponent(mg, N), n_terms, basis)
    return A[:, : coeffs.size] @ coeffs


def richardson_holdout(K_codes, masses, mg, N, n_terms=2, n_hold=2,
                       basis="paper"):
    """Fit on the low-K points, predict the highest ``n_hold``, report the error.

    The honest error estimate, and the one thing neither the paper's last-term
    rule nor :func:`richardson_stability` provides: **a measurement of the
    fit's actual extrapolating ability.**

    The last-term rule reports the size of the final fitted coefficient, which
    shrinks like ``K^-e`` as the window widens whether or not the answer is
    better determined -- measured over the full Table I, widening 2K = 25-35 to
    25-49 improves the last term by 1.65x while ``richardson_stability`` gets
    2x *worse* and agreement with the published values is unchanged.
    ``richardson_stability`` refits sub-windows, which detects instability but
    still never confronts the fit with data it has not seen.

    This does.  Holding out the largest K values is the right split rather than
    a random one: extrapolation to ``K -> infinity`` is precisely a prediction
    beyond the data, so the held-out points must lie on that side.  It is also
    the only criterion here that can catch a wrong *functional form*, which no
    amount of refitting the same basis will reveal.

    Returns ``(max_abs_err, max_rel_err, predicted, actual, K_held)``.
    """
    K_codes = np.asarray(K_codes, dtype=float)
    masses = np.asarray(masses, dtype=float)
    order = np.argsort(K_codes)
    K_codes, masses = K_codes[order], masses[order]

    if n_hold < 1 or len(K_codes) - n_hold < 2:
        raise ValueError(
            f"need at least 2 fitting points and 1 held out; got "
            f"{len(K_codes)} points with n_hold={n_hold}")

    K_fit, K_held = K_codes[:-n_hold], K_codes[-n_hold:]
    m_fit, m_held = masses[:-n_hold], masses[-n_hold:]

    n_terms = min(n_terms, max(len(K_fit) - 1, 0))
    _, _, coeffs, _ = richardson_extrapolate(
        K_fit, m_fit, mg, N, n_terms=n_terms, return_fit=True, basis=basis)

    pred = richardson_evaluate(coeffs, mg, N, K_held, n_terms=n_terms,
                               basis=basis)
    err = np.abs(pred - m_held)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(m_held != 0, err / np.abs(m_held), np.inf)
    return (float(err.max()), float(np.max(rel)), pred, m_held, K_held)


def wynn_epsilon(values):
    """Limit of a sequence by Wynn's epsilon algorithm.

    Assumes **nothing** about the exponents.  That is the point: every estimate
    in this module other than this one and :func:`levin_u` is built on the
    Eq. (27) basis, so when that basis is itself under suspicion -- and near the
    chiral limit it is, because ``1/K^(1+a)`` collapses onto ``1/K`` -- they
    cannot referee each other.  These two can.

    The recurrence is ``e[k+1][n] = e[k-1][n+1] + 1/(e[k][n+1] - e[k][n])``,
    with the even columns holding the estimates.  It is exact for a sequence
    whose remainder is a sum of geometrics and works well far beyond that, at
    the cost of being numerically delicate: the reciprocal of a small difference
    amplifies rounding, so the table is walked only while the differences stay
    meaningful, and the last stable even column is returned.

    Returns ``(limit, order)``; ``order`` is the column it came from, so a
    caller can tell a well-accelerated estimate from one that stalled at k=0.
    """
    s = np.asarray(values, dtype=float)
    n = s.size
    if n < 3:
        return (float(s[-1]) if n else float("nan")), 0

    prev = np.zeros(n)          # e[k-1]
    cur = s.copy()              # e[k]
    best, best_order = float(s[-1]), 0
    scale = max(float(np.abs(s).max()), 1e-300)

    for k in range(1, n):
        m = n - k
        nxt = np.empty(m)
        ok = True
        for i in range(m):
            d = cur[i + 1] - cur[i]
            # A difference at the level of rounding makes 1/d meaningless; the
            # table has converged as far as this data can take it.
            if abs(d) < 1e-14 * scale:
                ok = False
                break
            nxt[i] = prev[i + 1] + 1.0 / d
        if not ok:
            break
        prev, cur = cur[:m], nxt
        if k % 2 == 0 and np.all(np.isfinite(cur)):
            best, best_order = float(cur[-1]), k
        if m <= 1:
            break
    return best, best_order


def levin_u(K_codes, masses):
    """Limit of ``M^2(K)`` by the Levin u-transform.

    Like :func:`wynn_epsilon` this assumes no exponents, but unlike it the
    u-transform is built for exactly our situation: a sequence whose remainder
    has an asymptotic expansion in ``1/x``.  It needs a remainder estimate, and
    the "u" choice is ``omega_n = x_n * (S_n - S_{n-1})`` -- the increment scaled
    by the abscissa -- which is the right one when the tail behaves like a power
    series in ``1/x``.

    ``x`` is the paper's ``K``, not the sequence index, because the expansion is
    in ``1/K`` and the code's ``K_codes`` step by 2.

    Returns ``(limit, order)``.
    """
    from math import comb

    K = np.asarray(K_codes, dtype=float)
    S = np.asarray(masses, dtype=float)
    order_idx = np.argsort(K)
    K, S = K[order_idx] / 2.0, S[order_idx]
    n = S.size
    if n < 3:
        return (float(S[-1]) if n else float("nan")), 0

    dS = np.diff(S)
    omega = K[1:] * dS                      # u-transform remainder estimate
    S_ = S[1:]
    m = S_.size
    if np.any(omega == 0.0):
        return float(S[-1]), 0

    best, best_order = float(S[-1]), 0
    for k in range(1, m):
        num = den = 0.0
        for j in range(k + 1):
            w = ((-1) ** j) * comb(k, j)
            # Levin's (n+j+1)^(k-1) weight, written with the actual abscissae.
            w *= (K[j] / K[k]) ** (k - 1)
            num += w * S_[j] / omega[j]
            den += w / omega[j]
        if den == 0.0 or not np.isfinite(num / den):
            break
        best, best_order = float(num / den), k
    return best, best_order


def monotone_bracket(K_codes, masses, n_tail=4):
    """A genuine interval for ``M^2(K -> infinity)``, not a heuristic error bar.

    Measured at every K in this project, ``M^2(K)`` increases with K and its
    increments decrease.  Monotonicity alone gives the lower bound for free:
    the limit cannot be below the largest ``M^2`` computed.

    The upper bound needs one stated assumption -- that the increments continue
    to decay as a power law ``dM ~ C K^-p``, which is what the Eq. (26)/(27)
    structure predicts and what the data shows.  Then the remaining tail is
    bounded by the integral ``C K_max^(1-p) / (p-1)``.

    ``p`` is fitted on the last ``n_tail`` increments in log-log.  To keep the
    upper bound conservative the *shallowest* local slope over that window is
    used rather than the least-squares one: a smaller ``p`` decays more slowly
    and so over-estimates the tail, which is the safe direction for an upper
    bound.

    The tail is a *sum*, not an integral, and the two differ by about half a
    step -- enough to matter.  Approximating it by the integral from the last
    abscissa gave an upper bound that failed by 2e-5 on a series whose limit is
    known exactly, and a bound that can be beaten is worse than no bound.  For a
    decreasing summand the integral test gives the safe direction: the sum from
    the next increment on is at most the integral started half a step *before*
    the last point.

    Returns ``(lo, hi, p, monotone)``.  ``monotone`` is False if the sequence
    does not actually increase with decreasing increments, in which case the
    bracket is not meaningful and the caller must not use it.
    """
    K = np.asarray(K_codes, dtype=float)
    S = np.asarray(masses, dtype=float)
    o = np.argsort(K)
    K, S = K[o] / 2.0, S[o]
    if S.size < 4:
        return float("nan"), float("nan"), float("nan"), False

    d = np.diff(S)
    monotone = bool(np.all(d > 0) and np.all(np.diff(d) < 0))
    lo = float(S[-1])
    if not monotone:
        return lo, float("nan"), float("nan"), False

    mid = 0.5 * (K[1:] + K[:-1])
    t = min(n_tail, d.size - 1)
    lx, ly = np.log(mid[-t - 1:]), np.log(d[-t - 1:])
    slopes = np.diff(ly) / np.diff(lx)
    p = float(-np.max(slopes))              # shallowest decay = safest bound
    if p <= 1.0:
        # Tail not summable under this assumption; no finite upper bound.
        return lo, float("inf"), p, True

    C = d[-1] * mid[-1] ** p
    step = float(K[-1] - K[-2])
    # Integral test for a decreasing summand: start half a step early.
    tail = C * (K[-1] - 0.5 * step) ** (1.0 - p) / (p - 1.0)
    return lo, lo + float(tail), p, True


def richardson_stability(K_codes, masses, mg, N, n_terms=2, min_points=4):
    """How much ``M(0)`` moves when the fit window is changed.

    The paper reports "the magnitude of the last term in the series fit" as its
    uncertainty, which :func:`richardson_extrapolate` returns.  That rule breaks
    down at weak coupling, and not because K is too small.  The Eq. (27) basis
    is ``1, 1/K, 1/K^(1+a), ...`` with the Eq. (26) exponent ``a``, and ``a -> 0``
    in the chiral limit, so ``1/K^(1+a)`` collapses onto ``1/K``.  The fit stays
    excellent -- maximum residual 1e-5 or better at every coupling, five digits
    on M(K) -- but the coefficients blow up with near-total cancellation: 0.892
    and 0.098 at m/g = 1.6, against 8.383 and 8.800 at m/g = 0.1.  The last term
    is then large because the basis is degenerate, not because the answer is
    uncertain.

    This refits on every sub-window of at least ``min_points`` points and
    returns the standard deviation of the resulting ``M(0)``.  Measured against
    the last-term rule for the N = 3 baryon over 2K = 25-35, the latter
    overstates by 6x at m/g = 1.6 rising to **341x at m/g = 0.1**, where the
    extrapolation is stable to 0.2% rather than the 72% the rule reports.

    **What this does not measure.**  It is the fit's sensitivity to the window,
    not a total uncertainty.  Systematics that shift every point in the same
    direction are invisible to it -- in particular the ``sweep_lpn``
    particle-number truncation, and any finite-K effect the Eq. (27) form does
    not capture.  At weak coupling the residual disagreement with the published
    values is larger than either estimate, so something systematic remains.

    Returns ``(M0, spread, n_windows)``; ``spread`` is NaN when there are too
    few points to form more than one sub-window.
    """
    from itertools import combinations

    K_codes = list(K_codes)
    masses = list(masses)
    M0, _ = richardson_extrapolate(K_codes, masses, mg, N, n_terms=n_terms)
    n = len(K_codes)
    fits = []
    for size in range(min_points, n + 1):
        for idx in combinations(range(n), size):
            if size == n:
                continue
            sub = richardson_extrapolate([K_codes[i] for i in idx],
                                         [masses[i] for i in idx],
                                         mg, N, n_terms=n_terms)[0]
            fits.append(sub)
    if len(fits) < 2:
        return M0, float("nan"), len(fits)
    return M0, float(np.std(fits)), len(fits)


def richardson_ensemble(K_codes, masses, mg, N, n_terms_set=(2, 3, 4),
                        min_points=4, min_dof=1):
    """Every defensible Eq. (27) fit of this data: each order x each sub-window.

    An extrapolated mass is not one number but a choice -- how many correction
    terms to keep, and which K to fit.  Enumerating those choices is what turns
    "the fit residual is 1e-5" into an honest uncertainty, because the residual
    measures how well a *given* form describes the data and is blind to the
    form being wrong.

    ``min_dof`` requires at least that many degrees of freedom, so an
    exactly-determined fit (which has zero residual by construction and no
    predictive content) is never admitted.

    Returns a list of ``(coeffs, exponents, M0)``.
    """
    from itertools import combinations

    K_codes = list(K_codes)
    masses = list(masses)
    n = len(K_codes)
    out = []
    for nt in n_terms_set:
        for size in range(max(min_points, nt + 1 + min_dof), n + 1):
            for idx in combinations(range(n), size):
                kk = [K_codes[i] for i in idx]
                mm = [masses[i] for i in idx]
                try:
                    M0, _, coeffs, exps = richardson_extrapolate(
                        kk, mm, mg, N, n_terms=nt, return_fit=True)
                except Exception:
                    continue
                if np.isfinite(M0):
                    out.append((coeffs, exps, float(M0)))
    return out


def richardson_curve(coeffs, exponents, inv_K):
    """Evaluate a fitted Eq. (27) series at ``1/K_paper`` (0 = the continuum)."""
    inv_K = np.asarray(inv_K, dtype=float)
    y = np.zeros_like(inv_K)
    for c, e in zip(coeffs, exponents):
        y = y + (c if e == 0 else c * np.where(inv_K > 0, inv_K, 0.0) ** e)
    return y


def richardson_budget(K_codes, masses, mg, N, n_terms=2, n_terms_set=(2, 3, 4),
                      min_points=4, masses_alt_lpn=None):
    """Full error budget for one extrapolated mass.

    Four components, combined in quadrature, each measured rather than assumed:

    ``form``
        Spread of M(0) over the number of correction terms kept, on the full
        window.  **This dominates** -- 7-18x the window term at every point
        tested -- because the data sit at 1/K ~ 0.06 while the answer is at
        1/K = 0, so the curvature across the gap is set by the assumed series.
    ``window``
        Spread over sub-windows at fixed order (:func:`richardson_stability`).
        Honest but narrow: every sub-window shares the same form over the same
        range, so it cannot see the term above.
    ``truncation``
        Sensitivity to the particle-number cut, from ``masses_alt_lpn`` if
        given.  Measured at ~1e-4 -- raising ``sweep_lpn`` by a whole qqbar
        pair moves M(0) in the fourth decimal.
    ``numerical``
        The solver floor.  ``assembly="exact"`` is basis independent to 6e-13
        (docs/basis-dependence.md), so this is negligible and carried only so
        the budget is complete rather than selectively quoted.

    Returns a dict with ``M0`` and the components.
    """
    M0, last = richardson_extrapolate(K_codes, masses, mg, N, n_terms=n_terms)

    fits = []
    for nt in n_terms_set:
        try:
            fits.append(richardson_extrapolate(K_codes, masses, mg, N,
                                               n_terms=nt)[0])
        except Exception:
            pass
    fits = [f for f in fits if np.isfinite(f)]
    err_form = float(np.std(fits)) if len(fits) > 1 else float("nan")

    _, err_window, _ = richardson_stability(K_codes, masses, mg, N,
                                            n_terms=n_terms,
                                            min_points=min_points)

    err_trunc = 0.0
    if masses_alt_lpn is not None:
        try:
            M0b, _ = richardson_extrapolate(K_codes, masses_alt_lpn, mg, N,
                                            n_terms=n_terms)
            err_trunc = abs(M0b - M0)
        except Exception:
            err_trunc = float("nan")

    err_num = abs(M0) * 1e-12

    parts = [p for p in (err_form, err_window, err_trunc, err_num)
             if p is not None and np.isfinite(p)]
    total = float(np.sqrt(np.sum(np.square(parts)))) if parts else float("nan")
    return {"M0": float(M0), "last_term": float(last), "form": err_form,
            "window": float(err_window), "truncation": float(err_trunc),
            "numerical": err_num, "total": total}


def thooft_valence_limit(x, N, B=1):
    """Analytic chiral-limit baryon structure function, Eq. (22).

        q(x) = N (N-1) (1-x)^(N-2)

    N=3 gives 6(1-x); N=2 gives the constant 2, identical to the meson
    distribution.  A hard target for Fig. 3(b) as m/g -> 0.
    """
    x = np.asarray(x, dtype=float)
    return N * (N - 1) * (1.0 - x) ** (N - 2)
