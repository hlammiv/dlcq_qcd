"""qcdf_kernels.py -- nogil numba kernels for the DLCQ matrix build.

``qcdf_opt.py`` already JIT-compiles its micro-kernels, but ``clfact`` and
``hamqcd`` -- 79% and 7% of a run -- stayed interpreted, so the only way to use
more than one core was ``multiprocessing``.  That worked but plateaued at ~5x:
not because of IPC (measured at 2.4 MB per build, with accumulation costing
2 ms), but because a process pool cannot share the interpreter's work and this
class of code is dominated by interpreter overhead in the first place.

Compiling the two hot routines with ``nogil=True`` fixes both halves at once.
The serial path gets much faster, and because the kernels hold no lock a plain
``ThreadPoolExecutor`` gives true parallelism with no fork, no pickling, and
rows written straight into the destination matrix.

Everything here is a transcription of the routines in ``qcdf_opt.py``, which
remain in place as the reference implementation.  ``tests/test_kernels.py``
asserts the two agree exactly, element by element.  Three changes were needed
to make the code compilable, and all three are behaviour-preserving:

1. **The epsilon permutation table is hoisted.**  ``clfact`` rebuilt it via
   ``gprbig_array`` on *every* call -- 15409 times in one profiled row batch --
   for a table that depends only on ``N``.  It is now built once by
   :func:`epsb_table` (still using the original ``gprbig_array``/``psign``, so
   it is the same table) and passed in.
2. **Scratch buffers are per-thread.**  ``qcdf_opt`` keeps one module-global
   ``ClfactBuffers``, which is safe only because each forked process got its
   own copy.  Threads share memory, so each worker now owns a :class:`Scratch`.
3. **The vertex table is arrays, not lambdas.**  ``hamqcd``'s ``V4`` list holds
   closures; the same information is encoded below as integer tables, with the
   momentum arguments as coefficient vectors over ``(k0..k3)`` and the flavour
   assignment as a 0/1 selector over ``(f1, f2)``.
"""

from __future__ import annotations

import numpy as np
from numba import njit

from qcdf import gprbig_array, psign, MXNP, MXTRM

MXLNG = 54          # 2 * MXP + 4
MXP = 25

__all__ = ["Scratch", "epsb_table", "clfact_nb", "hamqcd_nb",
           "norm_row", "ham_row"]


# ──────────────────────────────────────────────────────────────────────────
# Vertex tables (were lambdas in qcdf_opt.hamqcd)
# ──────────────────────────────────────────────────────────────────────────
# Per four-point vertex: the operator signature checked by OPMTCH, the b/d
# occupation rows written into the four operator slots, the colour indices
# handed to clfact, the two BRACK arguments as coefficients over (k0..k3), the
# flavour pattern as a selector over (f1, f2), and the overall coefficient.

_V4_MATCH = np.array([[2, 2, 0, 0],
                      [0, 0, 2, 2],
                      [1, 2, 0, 1],
                      [0, 1, 1, 2],
                      [2, 1, 1, 0],
                      [1, 0, 2, 1]], dtype=np.int32)

_V4_ROW0 = np.array([[-1, -1, 1, 1],
                     [0, 0, 0, 0],
                     [0, -1, 1, 1],
                     [0, 0, 0, 1],
                     [0, -1, -1, 1],
                     [0, 0, 0, -1]], dtype=np.int32)

_V4_ROW1 = np.array([[0, 0, 0, 0],
                     [-1, -1, 1, 1],
                     [1, 0, 0, 0],
                     [-1, 1, 1, 0],
                     [-1, 0, 0, 0],
                     [-1, -1, 1, 0]], dtype=np.int32)

_V4_IC = np.array([[0, 3, 1, 2],
                   [2, 1, 3, 0],
                   [1, 3, 0, 2],
                   [1, 0, 2, 3],
                   [1, 0, 2, 3],
                   [3, 1, 2, 0]], dtype=np.int32)

# ka = _V4_KA . (k0,k1,k2,k3),  kb = _V4_KB . (k0,k1,k2,k3)
_V4_KA = np.array([[0, -1, 0, 1],       # k3 - k1
                   [0, 1, 0, -1],       # k1 - k3
                   [1, 0, 0, 1],        # k3 + k0
                   [1, 0, -1, 0],       # k0 - k2
                   [1, 0, 1, 0],        # k0 + k2
                   [0, 1, -1, 0]],      # k1 - k2
                  dtype=np.int32)

_V4_KB = np.array([[-1, 0, 1, 0],       # k2 - k0
                   [1, 0, -1, 0],       # k0 - k2
                   [0, -1, 1, 0],       # k2 - k1
                   [0, -1, 0, -1],      # -k3 - k1
                   [0, 1, 0, -1],       # k1 - k3
                   [1, 0, 0, 1]],       # k0 + k3
                  dtype=np.int32)

# 0 selects f1, 1 selects f2
_V4_FL = np.array([[0, 1, 0, 1],
                   [0, 1, 0, 1],
                   [0, 1, 1, 0],
                   [0, 1, 0, 1],
                   [0, 1, 0, 1],
                   [0, 1, 1, 0]], dtype=np.int32)

_V4_COEF = np.array([-0.5, -0.5, 1.0, 1.0, 1.0, 1.0])

_H7_ROW0 = np.array([0, 0, -1, 1], dtype=np.int32)
_H7_ROW1 = np.array([-1, 1, 0, 0], dtype=np.int32)


# ──────────────────────────────────────────────────────────────────────────
# Per-thread scratch
# ──────────────────────────────────────────────────────────────────────────

class Scratch:
    """Work buffers for one worker thread.

    ``qcdf_opt`` allocates one set at module scope; that is safe under ``fork``
    and unsafe under threads, so each thread gets its own.  The two ``MXTRM x
    MXLNG`` term tables dominate: ~90 MB per thread at ``MXTRM = 200000``,
    against ~5.5 MB at the inherited Fortran value of 12552, which was found to
    be silently truncating colour sums (docs/fortran-color-overflow.md).

    The reservation costs far less than it looks.  ``np.zeros`` hands back lazy
    zero pages, and the kernels only ever touch rows ``0..ntrms``, so what is
    resident tracks the actual term count rather than the cap.  Measured: after
    the increase 2K=29 runs in 1.52 s, against the 4.87 s recorded in
    docs/performance.md before it.
    """

    __slots__ = ("idel0", "resl0", "idelt", "reslt", "my", "lnkb", "icycl",
                 "ncr", "nperm", "nb1", "nb2", "mxw", "ltc", "mqn", "mqn_n",
                 "elem0", "fwd", "rev", "used")

    def __init__(self, NF: int, K: int):
        self.idel0 = np.zeros((MXTRM, MXLNG), dtype=np.int32)
        self.resl0 = np.zeros(MXTRM)
        self.idelt = np.zeros((MXTRM, MXLNG), dtype=np.int32)
        self.reslt = np.zeros(MXTRM)
        self.my = np.zeros((4, MXLNG), dtype=np.int32)
        self.lnkb = np.zeros(MXLNG, dtype=np.int32)
        self.icycl = np.ones(MXLNG, dtype=np.int32)
        self.ncr = np.zeros((2, MXLNG), dtype=np.int32)
        self.nperm = np.zeros((MXLNG, MXLNG), dtype=np.int32)
        self.nb1 = np.zeros(MXLNG, dtype=np.int32)
        self.nb2 = np.zeros(MXLNG, dtype=np.int32)
        self.mxw = np.zeros((4, MXLNG), dtype=np.int32)
        self.ltc = np.zeros((4, MXP), dtype=np.int32)
        self.mqn = np.zeros((4, MXLNG), dtype=np.int32)
        self.mqn_n = np.zeros(4, dtype=np.int32)
        self.elem0 = np.zeros(max(NF, 1))
        # MOCHK counters: momenta run to K, flavours 1..NF.
        self.fwd = np.zeros((K + 2, NF + 1), dtype=np.int32)
        self.rev = np.zeros((K + 2, NF + 1), dtype=np.int32)
        self.used = np.zeros(MXP, dtype=np.bool_)

    def clfact_args(self):
        """Buffers used by :func:`clfact_nb`, in order."""
        return (self.idel0, self.resl0, self.idelt, self.reslt, self.my,
                self.lnkb, self.icycl, self.ncr, self.nperm, self.nb1,
                self.nb2, self.fwd, self.rev)

    def norm_args(self):
        return self.clfact_args() + (self.mxw, self.ltc, self.used)

    def ham_args(self):
        return self.clfact_args() + (self.mxw, self.ltc, self.mqn,
                                     self.mqn_n, self.elem0, self.used)


def epsb_table(N: int):
    """The ``(nprms, ibrpm)`` epsilon-reduction table for SU(N).

    ``qcdf_opt.clfact`` computed this inline, per call.  It depends only on
    ``N``: the signed permutations of ``1..N-1`` used by BREDCE to expand a
    contracted pair of epsilon tensors into antisymmetrized deltas.  For N=3
    that is ``[[1,2,+1], [2,1,-1]]``.

    ``N=2`` has no extra terms, matching the ``nprms_ep = 0`` branch.
    """
    if N <= 2:
        return 0, np.zeros((1, MXNP), dtype=np.int64)
    nmin = np.ones(MXNP, dtype=np.int64)
    nmax = np.ones(MXNP, dtype=np.int64) * (N - 1)
    nprms, ibrpm = gprbig_array(nmin, nmax, N - 1, 1, 0, 9523)
    ibrpm = np.ascontiguousarray(ibrpm, dtype=np.int64)
    for i in range(nprms):
        ibrpm[i, N - 1] = psign(N - 1, ibrpm[i, :N - 1].tolist())
    return nprms, ibrpm


# ──────────────────────────────────────────────────────────────────────────
# Kernels
# ──────────────────────────────────────────────────────────────────────────

@njit(nogil=True, cache=True)
def _conj_state(lst, ist, out):
    for i in range(lst):
        j = lst - 1 - i
        out[0, j] = -ist[0, i]
        out[1, j] = -ist[1, i]
        out[2, j] = ist[2, i]
        out[3, j] = ist[3, i]


@njit(nogil=True, cache=True)
def _brack(L, M):
    if L == 0 or M == 0 or L + M != 0:
        return 0.0
    return 4.0 / (float(L) * float(L))


@njit(nogil=True, cache=True)
def _mi(u, v):
    """Slot index for an operator with b-occupation ``u`` and d-occupation ``v``.

    Maps (-1,0)->0 (annihilate quark), (0,-1)->1 (annihilate antiquark),
    (1,0)->2 (create quark), (0,1)->3 (create antiquark), selecting which of
    the four momentum lists the operator's momentum must come from.
    """
    return 2 * u * u + 3 * v * v + u + v - 1


@njit(nogil=True, cache=True)
def clfact_nb(ic1, ic2, ic3, ic4, mx, lng, llt, lrt, nops, N, NF, B, K,
              nprms_ep, ibrpm,
              idel0, resl0, idelt, reslt, my, lnkb, icycl, ncr, nperm,
              nb1, nb2, fwd, rev):
    """Colour factor for one operator string.  See ``qcdf_opt.clfact``."""
    for i in range(lng):
        idel0[0, i] = 0
        my[0, i] = mx[0, i]
        my[1, i] = mx[1, i]
        my[2, i] = mx[2, i]
        my[3, i] = mx[3, i]
    resl0[0] = 1.0

    # ── MOCHK: operator balance at each momentum+flavour ──
    for k3 in range(2):
        for a in range(fwd.shape[0]):
            for b in range(fwd.shape[1]):
                fwd[a, b] = 0
                rev[a, b] = 0
        for i3 in range(lng):
            ri = lng - 1 - i3
            mom_f = mx[2, i3]
            flv_f = mx[3, i3]
            mom_r = mx[2, ri]
            flv_r = mx[3, ri]
            fwd[mom_f, flv_f] += mx[k3, i3]
            rev[mom_r, flv_r] += mx[k3, ri]
            if (fwd[mom_f, flv_f] < 0 or fwd[mom_f, flv_f] > N or
                    rev[mom_r, flv_r] > 0 or rev[mom_r, flv_r] < -N):
                return 0.0

    ipow = 0
    for ja in range(1, lng):
        if my[0, ja] + my[1, ja] != -1:
            continue
        for ka in range(ja - 1, -1, -1):
            if (my[0, ka] == -my[0, ja] and my[1, ka] == -my[1, ja]
                    and my[2, ka] == my[2, ja] and my[3, ka] == my[3, ja]):
                idel0[0, ja] = ka
                idel0[0, ka] = ja
                my[0, ka] = 0
                my[1, ka] = 0
                my[0, ja] = 0
                my[1, ja] = 0
                break
            elif my[0, ka] != 0 or my[1, ka] != 0:
                ipow += 1
    resl0[0] = 1.0 if ipow % 2 == 0 else -1.0

    for i in range(lng):
        lnkb[i] = 0
        icycl[i] = 1
    if B:
        for nbs in range(B):
            for i in range(N * nbs, N * (nbs + 1) - 1):
                lnkb[i] = i + 1
            lnkb[N * (nbs + 1) - 1] = N * nbs
            for j in range(lng - N * nbs - 1, lng - N * (nbs + 1), -1):
                lnkb[j] = j - 1
            lnkb[lng - N * (nbs + 1)] = lng - N * nbs - 1
        if N % 2 == 0:
            for nbs in range(B):
                s = 1
                for i in range(N * nbs, N * (nbs + 1)):
                    icycl[i] = s
                    s *= -1
                s = 1
                for j in range(lng - N * nbs - 1, lng - N * (nbs + 1) - 1, -1):
                    icycl[j] = s
                    s *= -1

    for i in range(lng):
        ncr[0, i] = 0
        ncr[1, i] = 0
        for j in range(lng):
            nperm[i, j] = 0

    lbcr = 0
    ldcr = 0
    lprm = 0
    lpo = 1
    ntrms = 1

    for i in range(lng):
        for bd in range(2):
            is_b = (bd == 0 and mx[1, i] == 0)
            is_d = (bd == 1 and mx[0, i] == 0 and mx[1, i] != 0)
            if not (is_b or is_d):
                continue
            cr = bd
            lcr = lbcr if bd == 0 else ldcr
            if mx[bd, i] == -1:
                lprm += 1
                lpchk = False
                jfst = -1
                for j in range(lcr - 1, -1, -1):
                    if (not lpchk and mx[2, ncr[cr, j]] == mx[2, i]
                            and mx[3, ncr[cr, j]] == mx[3, i]):
                        jfst = j
                        for ll in range(j - 1, -1, -1):
                            if (mx[2, ncr[cr, ll]] == mx[2, ncr[cr, j]]
                                    and mx[3, ncr[cr, ll]] == mx[3, ncr[cr, j]]):
                                nperm[lprm - 1, 0] += 1
                                mp = 2 * nperm[lprm - 1, 0]
                                nperm[lprm - 1, mp - 1] = ncr[cr, j]
                                nperm[lprm - 1, mp] = ncr[cr, ll]
                                lpchk = True
                        break
                if not lpchk:
                    lprm -= 1
                if jfst >= 0:
                    for m in range(jfst, lcr - 1):
                        ncr[cr, m] = ncr[cr, m + 1]
                    ncr[cr, lcr - 1] = 0
                    if bd == 0:
                        lbcr -= 1
                    else:
                        ldcr -= 1
            elif mx[bd, i] == 1:
                ncr[cr, lcr] = i
                if bd == 0:
                    lbcr += 1
                else:
                    ldcr += 1

    if lprm > 0:
        for pid in range(lprm - 1, -1, -1):
            nsbar = 1
            for jd in range(nperm[pid, 0]):
                kd = 2 * (jd + 1)
                ip1 = nperm[pid, kd - 1]
                ip2 = nperm[pid, kd]
                ismb = False
                if B > 0 and lnkb[ip1] != 0:
                    nxt = lnkb[ip1]
                    for _ in range(N - 1):
                        if ip2 == nxt:
                            ismb = True
                            break
                        nxt = lnkb[nxt]
                if ismb:
                    nsbar += 1
                else:
                    for j in range(ntrms):
                        if lpo >= idel0.shape[0]:
                            # Guard on the buffer's ACTUAL length, not on the
                            # MXTRM global.  numba freezes globals at compile
                            # time and ``cache=True`` does not notice MXTRM
                            # changing, so a kernel compiled at a high cap can be
                            # reused against buffers Python allocated at a low
                            # one.  njit does not bounds check, so that combination
                            # writes past the end of the array -- a far worse bug
                            # than the overflow this guard exists to catch.
                            # Reading .shape[0] is free and always correct.
                            # Silently returning 0.0 here produced a norm matrix
                            # that was not positive semidefinite (N=4, B=1,
                            # 2K=20: a state's own norm came out 150048 against
                            # a true 331776, and the matrix had an eigenvalue of
                            # -1.8e5).  Both backends share MXTRM, so they agreed
                            # bit-for-bit on the wrong answer and the
                            # backend-equality tests could not see it.
                            raise OverflowError(
                                "colour contraction exceeded MXTRM terms; the "
                                "result would be silently wrong. Raise it with "
                                "the DLCQ_MXTRM environment variable (see MXTRM "
                                "in qcdf.py), or reduce the cutoff LPN.")
                        for t in range(lng):
                            idel0[lpo, t] = idel0[j, t]
                        resl0[lpo] = -resl0[j]
                        idel0[lpo, ip1] = idel0[j, ip2]
                        idel0[lpo, ip2] = idel0[j, ip1]
                        v1 = idel0[j, ip2]
                        v2 = idel0[j, ip1]
                        if 0 <= v1 < lng:
                            idel0[lpo, v1] = idel0[j, v2]
                        if 0 <= v2 < lng:
                            idel0[lpo, v2] = idel0[j, v1]
                        lpo += 1
            if nsbar > 1:
                for im in range(ntrms):
                    resl0[im] *= nsbar
            ntrms *= (nperm[pid, 0] + 2 - nsbar)

    ntrms0 = ntrms
    nctf = 2 if nops == 4 else 1
    total = 0.0

    for nct in range(1, nctf + 1):
        for t in range(ntrms):
            for j in range(lng):
                idelt[t, j] = 0
            reslt[t] = resl0[t]
        for j in range(lng):
            if mx[0, j] - mx[1, j] == 1:
                for t in range(ntrms):
                    idelt[t, j] = idel0[t, j]
        if nops == 4:
            if nct == 1:
                for t in range(ntrms):
                    idelt[t, lrt + ic1] = lrt + ic2
                    idelt[t, lrt + ic3] = lrt + ic4
                    reslt[t] *= 0.5
            else:
                for t in range(ntrms):
                    idelt[t, lrt + ic1] = lrt + ic4
                    idelt[t, lrt + ic3] = lrt + ic2
                    reslt[t] *= -0.5 / N
        elif nops == 2:
            for t in range(ntrms):
                idelt[t, lrt + ic1] = lrt + ic2

        for il in range(ntrms):
            for im in range(N * B + 1, lrt, 2):
                idelt[il, im] = im - 1
            for in_ in range(lng - N * B - 1, lng - llt, -2):
                idelt[il, in_] = in_ - 1

        for nt in range(ntrms):
            for lt in range(lng):
                ic = lt
                for _ in range(2 * lng):
                    nx = idelt[nt, ic]
                    if nx == 0:
                        break
                    nx2 = idelt[nt, nx]
                    if nx2 == 0:
                        break
                    idelt[nt, ic] = nx2
                    idelt[nt, nx] = 0
                    if nx2 == ic:
                        idelt[nt, ic] = 0
                        reslt[nt] *= N
                        break

        # ── Baryon epsilon reduction (BREDCE) ──
        if B:
            for _iter in range(2 * B + 1):
                for nt in range(ntrms):
                    for lt in range(lng):
                        ic = lt
                        for __ in range(2 * lng):
                            nx = idelt[nt, ic]
                            if nx == 0:
                                break
                            nx2 = idelt[nt, nx]
                            if nx2 == 0:
                                break
                            idelt[nt, ic] = nx2
                            idelt[nt, nx] = 0
                            if nx2 == ic:
                                idelt[nt, ic] = 0
                                reslt[nt] *= N
                                break

                ntpr = ntrms
                has_terms = False
                for nt in range(ntrms):
                    for lt in range(lng):
                        if idelt[nt, lt] != 0:
                            has_terms = True
                            inb1 = lt
                            inb2 = idelt[nt, lt]
                            reslt[nt] *= icycl[inb1] * icycl[inb2]
                            c1 = inb1
                            c2 = inb2
                            for i1 in range(N - 1):
                                nb1[i1] = lnkb[c1]
                                nb2[i1] = lnkb[c2]
                                c1 = lnkb[c1]
                                c2 = lnkb[c2]
                            idelt[nt, lt] = 0
                            for j1 in range(N - 1):
                                idelt[nt, nb2[j1]] = nb1[j1]
                            for j2 in range(1, nprms_ep):
                                if ntpr >= idelt.shape[0]:
                                    # Buffer length, not the MXTRM global -- see
                                    # the note on the other guard above.
                                    # ``break`` here left a *partial sum*, which
                                    # is how wrong non-zero values arose rather
                                    # than obvious zeros.  See the note on the
                                    # other MXTRM guard above.
                                    raise OverflowError(
                                        "colour contraction exceeded MXTRM "
                                        "terms; the result would be silently "
                                        "wrong. Raise it with the DLCQ_MXTRM "
                                        "environment variable (see MXTRM in "
                                        "qcdf.py), or reduce LPN.")
                                reslt[ntpr] = reslt[nt] * float(ibrpm[j2, N - 1])
                                for t in range(lng):
                                    idelt[ntpr, t] = idelt[nt, t]
                                for j4 in range(N - 1):
                                    idelt[ntpr, nb2[ibrpm[j2, j4] - 1]] = nb1[j4]
                                ntpr += 1
                            break
                ntrms = ntpr
                if not has_terms:
                    break

        ct = 0.0
        for i in range(ntrms):
            ct += reslt[i]
        total += ct
        ntrms = ntrms0
    return total


@njit(nogil=True, cache=True)
def hamqcd_nb(lrt, llt, matrt, matlt, N, NF, B, K, selfen, cbreak,
              nprms_ep, ibrpm,
              idel0, resl0, idelt, reslt, my, lnkb, icycl, ncr, nperm,
              nb1, nb2, fwd, rev, mxw, ltc, mqn, mqn_n, elem0):
    """One Hamiltonian matrix element.  See ``qcdf_opt.hamqcd``.

    ``elem0`` is filled in place (one entry per flavour); the interacting part
    is returned.
    """
    _conj_state(llt, matlt, ltc)
    for i in range(NF):
        elem0[i] = 0.0
    elem = 0.0

    for i in range(4):
        mqn_n[i] = 0
    for i in range(lrt):
        idx = 0 if matrt[0, i] == 1 else 1
        m = matrt[2, i]
        seen = False
        for t in range(mqn_n[idx]):
            if mqn[idx, t] == m:
                seen = True
                break
        if not seen:
            mqn[idx, mqn_n[idx]] = m
            mqn_n[idx] += 1
    for j in range(llt):
        idx = 2 if ltc[0, j] == -1 else 3
        m = ltc[2, j]
        seen = False
        for t in range(mqn_n[idx]):
            if mqn[idx, t] == m:
                seen = True
                break
        if not seen:
            mqn[idx, mqn_n[idx]] = m
            mqn_n[idx] += 1

    nbr = 0
    ndr = 0
    for i in range(lrt):
        if matrt[0, i] == 1:
            nbr += 1
        if matrt[1, i] == 1:
            ndr += 1
    nbl = 0
    ndl = 0
    for i in range(llt):
        if matlt[0, i] == 1:
            nbl += 1
        if matlt[1, i] == 1:
            ndl += 1

    # ── four-point vertices ──
    for v in range(6):
        a = _V4_MATCH[v, 0]
        b = _V4_MATCH[v, 1]
        c = _V4_MATCH[v, 2]
        d = _V4_MATCH[v, 3]
        if not (a <= nbr and b <= nbl and c <= ndr and d <= ndl
                and a + nbl == b + nbr and c + ndl == d + ndr):
            continue
        for i in range(lrt):
            for r in range(4):
                mxw[r, i] = matrt[r, i]
        for j in range(llt):
            for r in range(4):
                mxw[r, lrt + 4 + j] = ltc[r, j]
        for j in range(4):
            mxw[0, lrt + j] = _V4_ROW0[v, j]
            mxw[1, lrt + j] = _V4_ROW1[v, j]
        lng = 4 + lrt + llt

        m0 = _mi(mxw[0, lrt], mxw[1, lrt])
        m1 = _mi(mxw[0, lrt + 1], mxw[1, lrt + 1])
        m2 = _mi(mxw[0, lrt + 2], mxw[1, lrt + 2])
        m3 = _mi(mxw[0, lrt + 3], mxw[1, lrt + 3])

        ic1 = _V4_IC[v, 0]
        ic2 = _V4_IC[v, 1]
        ic3 = _V4_IC[v, 2]
        ic4 = _V4_IC[v, 3]
        coeff = _V4_COEF[v]
        h = 0.0
        for l1 in range(mqn_n[m0]):
            for l2 in range(mqn_n[m1]):
                for l3 in range(mqn_n[m2]):
                    for l4 in range(mqn_n[m3]):
                        k0 = mqn[m0, l1]
                        k1 = mqn[m1, l2]
                        k2 = mqn[m2, l3]
                        k3 = mqn[m3, l4]
                        mxw[2, lrt] = k0
                        mxw[2, lrt + 1] = k1
                        mxw[2, lrt + 2] = k2
                        mxw[2, lrt + 3] = k3
                        tot = 0
                        for j in range(4):
                            tot += mxw[2, lrt + j] * (mxw[0, lrt + j]
                                                      + mxw[1, lrt + j])
                        if tot != 0:
                            continue
                        ka = (_V4_KA[v, 0] * k0 + _V4_KA[v, 1] * k1
                              + _V4_KA[v, 2] * k2 + _V4_KA[v, 3] * k3)
                        kb = (_V4_KB[v, 0] * k0 + _V4_KB[v, 1] * k1
                              + _V4_KB[v, 2] * k2 + _V4_KB[v, 3] * k3)
                        bk = _brack(ka, kb)
                        if bk == 0.0:
                            continue
                        for f1 in range(1, NF + 1):
                            for f2 in range(1, NF + 1):
                                for j in range(4):
                                    mxw[3, lrt + j] = f1 if _V4_FL[v, j] == 0 else f2
                                h += coeff * clfact_nb(
                                    ic1, ic2, ic3, ic4, mxw, lng, llt, lrt, 4,
                                    N, NF, B, K, nprms_ep, ibrpm,
                                    idel0, resl0, idelt, reslt, my, lnkb,
                                    icycl, ncr, nperm, nb1, nb2, fwd, rev) * bk
        elem += h

    # ── H7 (exchange, two flavour structures) ──
    if (1 <= nbr and 1 <= nbl and 1 <= ndr and 1 <= ndl
            and 1 + nbl == 1 + nbr and 1 + ndl == 1 + ndr):
        for i in range(lrt):
            for r in range(4):
                mxw[r, i] = matrt[r, i]
        for j in range(llt):
            for r in range(4):
                mxw[r, lrt + 4 + j] = ltc[r, j]
        for j in range(4):
            mxw[0, lrt + j] = _H7_ROW0[j]
            mxw[1, lrt + j] = _H7_ROW1[j]
        lng = 4 + lrt + llt

        m0 = _mi(mxw[0, lrt], mxw[1, lrt])
        m1 = _mi(mxw[0, lrt + 1], mxw[1, lrt + 1])
        m2 = _mi(mxw[0, lrt + 2], mxw[1, lrt + 2])
        m3 = _mi(mxw[0, lrt + 3], mxw[1, lrt + 3])

        h = 0.0
        for l1 in range(mqn_n[m0]):
            for l2 in range(mqn_n[m1]):
                for l3 in range(mqn_n[m2]):
                    for l4 in range(mqn_n[m3]):
                        k0 = mqn[m0, l1]
                        k1 = mqn[m1, l2]
                        k2 = mqn[m2, l3]
                        k3 = mqn[m3, l4]
                        mxw[2, lrt] = k0
                        mxw[2, lrt + 1] = k1
                        mxw[2, lrt + 2] = k2
                        mxw[2, lrt + 3] = k3
                        tot = 0
                        for j in range(4):
                            tot += mxw[2, lrt + j] * (mxw[0, lrt + j]
                                                      + mxw[1, lrt + j])
                        if tot != 0:
                            continue
                        bka = _brack(k3 - k2, k1 - k0)
                        bkc = _brack(k3 + k1, -k0 - k2)
                        if bka != 0.0:
                            for f1 in range(1, NF + 1):
                                for f2 in range(1, NF + 1):
                                    mxw[3, lrt] = f1
                                    mxw[3, lrt + 1] = f1
                                    mxw[3, lrt + 2] = f2
                                    mxw[3, lrt + 3] = f2
                                    h -= clfact_nb(
                                        1, 3, 2, 0, mxw, lng, llt, lrt, 4,
                                        N, NF, B, K, nprms_ep, ibrpm,
                                        idel0, resl0, idelt, reslt, my, lnkb,
                                        icycl, ncr, nperm, nb1, nb2, fwd, rev) * bka
                        if bkc != 0.0:
                            for f1 in range(1, NF + 1):
                                for f2 in range(1, NF + 1):
                                    mxw[3, lrt] = f1
                                    mxw[3, lrt + 1] = f2
                                    mxw[3, lrt + 2] = f1
                                    mxw[3, lrt + 3] = f2
                                    h += clfact_nb(
                                        2, 3, 1, 0, mxw, lng, llt, lrt, 4,
                                        N, NF, B, K, nprms_ep, ibrpm,
                                        idel0, resl0, idelt, reslt, my, lnkb,
                                        icycl, ncr, nperm, nb1, nb2, fwd, rev) * bkc
        elem += h

    # ── H8, H9 (diagonal self-energy) ──
    for bd in range(2):
        if bd == 0:
            o00 = -1; o01 = 1; o10 = 0; o11 = 0
            c1 = 0; c2 = 1; sgn = 1.0
            a = 1; b = 1; c = 0; d = 0
        else:
            o00 = 0; o01 = 0; o10 = -1; o11 = 1
            c1 = 1; c2 = 0; sgn = -1.0
            a = 0; b = 0; c = 1; d = 1
        if not (a <= nbr and b <= nbl and c <= ndr and d <= ndl
                and a + nbl == b + nbr and c + ndl == d + ndr):
            continue
        for i in range(lrt):
            for r in range(4):
                mxw[r, i] = matrt[r, i]
        for j in range(llt):
            for r in range(4):
                mxw[r, lrt + 2 + j] = ltc[r, j]
        mxw[0, lrt] = o00
        mxw[0, lrt + 1] = o01
        mxw[1, lrt] = o10
        mxw[1, lrt + 1] = o11
        lng = 2 + lrt + llt
        s0 = _mi(o00, o10)
        s1 = _mi(o01, o11)

        for l1 in range(mqn_n[s0]):
            for l2 in range(mqn_n[s1]):
                mxw[2, lrt] = mqn[s0, l1]
                mxw[2, lrt + 1] = mqn[s1, l2]
                if (mxw[2, lrt] * (mxw[0, lrt] + mxw[1, lrt])
                        + mxw[2, lrt + 1] * (mxw[0, lrt + 1]
                                             + mxw[1, lrt + 1])) != 0:
                    continue
                k1 = mxw[2, lrt]
                k2 = mxw[2, lrt + 1]
                if k1 == 0 or k1 != k2:
                    continue
                for ifl in range(1, NF + 1):
                    mxw[3, lrt] = ifl
                    mxw[3, lrt + 1] = ifl
                    cf = clfact_nb(c1, c2, 0, 0, mxw, lng, llt, lrt, 2,
                                   N, NF, B, K, nprms_ep, ibrpm,
                                   idel0, resl0, idelt, reslt, my, lnkb,
                                   icycl, ncr, nperm, nb1, nb2, fwd, rev)
                    elem0[ifl - 1] += cf * (2.0 + sgn * cbreak) / float(k1)
                    if 0 < k1 <= 99:
                        elem += cf * 0.5 * selfen[k1 - 1]
    return elem


# ──────────────────────────────────────────────────────────────────────────
# Row workers -- these write straight into the destination matrices
# ──────────────────────────────────────────────────────────────────────────

@njit(nogil=True, cache=True)
def _overlap(matrt, lstr, matlt, lstl, used):
    """Partons the two states share, split into quarks (sb) and antiquarks (sd)."""
    for i in range(lstl):
        used[i] = False
    sb = 0
    sd = 0
    for ir in range(lstr):
        for il in range(lstl):
            if used[il]:
                continue
            if (matrt[0, ir] == matlt[0, il] and matrt[2, ir] == matlt[2, il]
                    and matrt[3, ir] == matlt[3, il]):
                if matlt[0, il] == 1:
                    sb += 1
                else:
                    sd += 1
                used[il] = True
                break
    return sb, sd


@njit(nogil=True, cache=True)
def norm_row(nstrt, cols, mstate, mstinf, N, NF, B, K, out_n,
             nprms_ep, ibrpm,
             idel0, resl0, idelt, reslt, my, lnkb, icycl, ncr, nperm,
             nb1, nb2, fwd, rev, mxw, ltc, used):
    """Norm-matrix row ``nstrt``, restricted to ``cols``, written to ``out_n``."""
    lr = mstinf[nstrt, 0] - 1
    lstr = mstinf[nstrt, 1]
    matrt = mstate[lr:lr + 4, :MXP]
    for ci in range(len(cols)):
        nstlt = cols[ci]
        ll = mstinf[nstlt, 0] - 1
        lstl = mstinf[nstlt, 1]
        matlt = mstate[ll:ll + 4, :MXP]
        sb, sd = _overlap(matrt, lstr, matlt, lstl, used)
        nbr = 0
        ndr = 0
        for i in range(lstr):
            if matrt[0, i] == 1:
                nbr += 1
            if matrt[1, i] == 1:
                ndr += 1
        nbl = 0
        ndl = 0
        for i in range(lstl):
            if matlt[0, i] == 1:
                nbl += 1
            if matlt[1, i] == 1:
                ndl += 1
        if abs(nbr - sb) + abs(nbl - sb) + abs(ndr - sd) + abs(ndl - sd) != 0:
            continue
        _conj_state(lstl, matlt, ltc)
        for i in range(lstr):
            for r in range(4):
                mxw[r, i] = matrt[r, i]
        for j in range(lstl):
            for r in range(4):
                mxw[r, lstr + j] = ltc[r, j]
        out_n[nstlt] = clfact_nb(0, 0, 0, 0, mxw, lstr + lstl, lstl, lstr, 0,
                                 N, NF, B, K, nprms_ep, ibrpm,
                                 idel0, resl0, idelt, reslt, my, lnkb, icycl,
                                 ncr, nperm, nb1, nb2, fwd, rev)


@njit(nogil=True, cache=True)
def ham_row(nstrt, numsta, mstate, mstinf, N, NF, B, K, selfen, cbreak,
            out_h0, out_h,
            nprms_ep, ibrpm,
            idel0, resl0, idelt, reslt, my, lnkb, icycl, ncr, nperm,
            nb1, nb2, fwd, rev, mxw, ltc, mqn, mqn_n, elem0, used):
    """Hamiltonian row ``nstrt`` (upper triangle) written into ``out_h0``/``out_h``."""
    lr = mstinf[nstrt, 0] - 1
    lstr = mstinf[nstrt, 1]
    matrt = mstate[lr:lr + 4, :MXP]
    for nstlt in range(nstrt, numsta):
        ll = mstinf[nstlt, 0] - 1
        lstl = mstinf[nstlt, 1]
        matlt = mstate[ll:ll + 4, :MXP]
        sb, sd = _overlap(matrt, lstr, matlt, lstl, used)
        nbr = 0
        ndr = 0
        for i in range(lstr):
            if matrt[0, i] == 1:
                nbr += 1
            if matrt[1, i] == 1:
                ndr += 1
        nbl = 0
        ndl = 0
        for i in range(lstl):
            if matlt[0, i] == 1:
                nbl += 1
            if matlt[1, i] == 1:
                ndl += 1
        if abs(nbr - sb) + abs(nbl - sb) + abs(ndr - sd) + abs(ndl - sd) > 4:
            continue
        e = hamqcd_nb(lstr, lstl, matrt, matlt, N, NF, B, K, selfen, cbreak,
                      nprms_ep, ibrpm,
                      idel0, resl0, idelt, reslt, my, lnkb, icycl, ncr, nperm,
                      nb1, nb2, fwd, rev, mxw, ltc, mqn, mqn_n, elem0)
        for ifl in range(NF):
            out_h0[ifl, nstlt] = elem0[ifl]
        out_h[nstlt] = e
