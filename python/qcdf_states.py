"""qcdf_states.py -- nogil numba kernels for Fock-state generation.

After the matrix build was compiled (``qcdf_kernels.py``), state generation
became the largest single cost in a run -- 41% at 2K=29, and 228 s at 2K=41.
This is the same treatment applied to ``qcdf.py``'s ``brmsgn`` and its helpers.
``qcdf.py`` keeps its versions as the interpreted reference; ``tests/
test_states.py`` asserts the two produce **element-identical** ``mstate`` and
``mstinf``, not merely the same set of states.

That distinction matters: state *order* is load-bearing.  ``WEEDR``/``WEEDR2``
drop linearly dependent states in index order, so a permutation of the basis
would silently change which states survive and, through that, the spectrum.

Three things changed beyond compilation, and only the third alters any
intermediate quantity (never a stored one):

1. **``stachk`` counts by comparison, not by dict.**  The original builds two
   Python dicts keyed ``val+1`` (non-negative) and ``-val+1`` (negative).  Those
   keyings are injective on their own domains and live in separate dicts, so the
   function is exactly "does any value occur more than ``indx`` times".  With
   ``nptcl <= MXNP = 25`` a pairwise count is equivalent, allocation-free, and
   removes 364,389 dict lookups per run.

2. **Work buffers are caller-owned.**  ``gprbig_array`` allocated a fresh
   ``(ISTMAX_BG, MXNP)`` array -- 1.9 MB -- on each of its 1129 calls per run,
   writing ~67 rows on average.  ``qcdf.py:658,664,670`` allocated the same
   1.9 MB for *absent* species and never read it.  A :class:`Scratch` grows on
   demand instead, which is what makes raising ``ISTMAX_BG`` for large K
   affordable at all: the cap needs ~150,000 at 2K=48, i.e. 30 MB per call.

3. **The momentum-parity test is hoisted.**  ``oddchk`` rejected 96.3% of
   candidates *after* each was fully built -- 209,152 constructed for 4,529 kept
   at 2K=29.  But a candidate's momenta come from the ``kpx``/``kpm`` row picked
   by the quotient ``(ist[l,j]-1)//dim``, and the flavour from the remainder, so
   parity is a property of ``l`` alone and is independent across the baryon,
   antibaryon and meson groups.  Testing each group's rows once and iterating
   only over survivors is exactly equivalent -- ``oddchk`` is a conjunction over
   those same groups -- and collapses the 46:1 funnel.
"""

from __future__ import annotations

import numpy as np
from numba import njit

MXNP = 25

__all__ = ["Scratch", "stachk_nb", "gprbig_nb", "brmsgn_nb",
           "OVERFLOW_PERM", "OVERFLOW_STATES"]

#: ``brmsgn_nb`` return codes.  Both are recoverable: the driver grows the
#: offending buffer and re-runs the sector from the same start index.
OVERFLOW_PERM = -1        # istb/istd/istm too small
OVERFLOW_STATES = -2      # output mstate/mstinf too small


class Scratch:
    """Per-thread work buffers for one ``brmsgn_nb`` call.

    ``istmax`` sizes the three permutation tables and grows on demand via
    :meth:`grown`, so a small run never pays for a large one.
    """

    __slots__ = ("istb", "istd", "istm", "nuterm", "nminbf", "nmaxbf",
                 "nmindf", "nmaxdf", "nminmf", "nmaxmf", "istate",
                 "ichkb", "ichkd", "fct", "okb", "okd", "okm", "istmax")

    def __init__(self, NF: int, istmax: int = 16384):
        self.istmax = int(istmax)
        z = lambda *s: np.zeros(s, dtype=np.int64)          # noqa: E731
        self.istb = z(self.istmax, MXNP)
        self.istd = z(self.istmax, MXNP)
        self.istm = z(self.istmax, MXNP)
        self.okb = z(self.istmax)
        self.okd = z(self.istmax)
        self.okm = z(self.istmax)
        self.nuterm = z(MXNP)
        self.nminbf, self.nmaxbf = z(MXNP), z(MXNP)
        self.nmindf, self.nmaxdf = z(MXNP), z(MXNP)
        self.nminmf, self.nmaxmf = z(MXNP), z(MXNP)
        self.istate = z(4, MXNP)
        self.ichkb, self.ichkd = z(MXNP), z(MXNP)
        self.fct = z(max(NF, 1))

    def grown(self, NF: int) -> "Scratch":
        """A fresh Scratch with twice the permutation capacity."""
        return Scratch(NF, self.istmax * 2)

    def args(self):
        """The buffers, in the order ``brmsgn_nb`` expects."""
        return (self.istb, self.istd, self.istm, self.okb, self.okd, self.okm,
                self.nuterm, self.nminbf, self.nmaxbf, self.nmindf,
                self.nmaxdf, self.nminmf, self.nmaxmf, self.istate,
                self.ichkb, self.ichkd, self.fct)


# ──────────────────────────────────────────────────────────────────────────
# Leaf kernels
# ──────────────────────────────────────────────────────────────────────────

@njit(nogil=True, cache=True)
def stachk_nb(nptcl, indx, nuterm):
    """Fermi statistics: no value may occur more than ``indx`` times.

    Equivalent to ``qcdf.stachk``; see this module's docstring for why the
    two-dict formulation reduces to a plain occurrence count.
    """
    for i in range(nptcl):
        c = 0
        v = nuterm[i]
        for j in range(nptcl):
            if nuterm[j] == v:
                c += 1
                if c > indx:
                    return False
    return True


@njit(nogil=True, cache=True)
def _pavill_nb(istate, lngsta, cutoff, K, rmq):
    """Pauli-Villars cutoff.  A no-op when ``cutoff <= 0``, which is the default."""
    if cutoff <= 0.0:
        return True
    eps = 1.0e-6
    elc = 0.0
    for ic in range(lngsta):
        x = float(istate[2, ic]) / float(K)
        if x < eps:
            return False
        rm = rmq[istate[3, ic] - 1]
        elc += rm * rm / x
    return elc < cutoff


@njit(nogil=True, cache=True)
def _flvchk_nb(istate, lngsta, NF, iflv, fct):
    """Flavour quantum numbers must match ``iflv`` exactly."""
    for i in range(NF):
        fct[i] = 0
    for ic in range(lngsta):
        mf = istate[3, ic] - 1
        if istate[0, ic] == 1:
            fct[mf] += 1
        else:
            fct[mf] -= 1
    for i in range(NF):
        if fct[i] != iflv[i]:
            return False
    return True


@njit(nogil=True, cache=True)
def gprbig_nb(nmin, nmax, lngth, indx, idbl, ist, nuterm):
    """Generalized permutation generator, writing into caller-owned ``ist``.

    A transcription of ``qcdf.gprbig_array`` with its odometer control flow
    intact.  Returns the number of rows written, or ``OVERFLOW_PERM`` if ``ist``
    is too small -- the caller grows it and retries rather than raising, since
    the required size is only knowable at runtime and grows ~1.6x per +2 in K.
    """
    istmax = ist.shape[0]
    nprmg = 1
    idirec = 1
    jj = 0

    for ll in range(lngth):
        nuterm[ll] = nmin[ll]

    while True:
        km = nmax[jj] if jj < lngth else 0
        if idirec == 1 or nuterm[jj] < km:
            if idirec != 1:
                nuterm[jj] += 1
            if (idbl == 1 and jj > 0 and nmin[jj] == nmin[jj - 1]
                    and (nuterm[jj] - 1) < nuterm[jj - 1]):
                nuterm[jj] = nuterm[jj - 1]
            jj += 1
            idirec = 1
            if jj >= lngth - 1:
                if (idbl == 1 and lngth > 1 and nmin[lngth - 1] == nmin[lngth - 2]
                        and (nuterm[lngth - 1] - 1) < nuterm[lngth - 2]):
                    nuterm[lngth - 1] = nuterm[lngth - 2]
                for ii in range(nuterm[lngth - 1], nmax[lngth - 1] + 1):
                    nuterm[lngth - 1] = ii
                    if indx != 0:
                        ixok = stachk_nb(lngth, indx, nuterm)
                    else:
                        ixok = True
                    if ixok:
                        if nprmg > istmax:
                            return OVERFLOW_PERM
                        for j in range(lngth):
                            ist[nprmg - 1, j] = nuterm[j]
                        nprmg += 1
                nuterm[lngth - 1] = nmin[lngth - 1]
                jj = lngth - 2
                if jj < 0:
                    break
                idirec = -1
        else:
            nuterm[jj] = nmin[jj]
            jj -= 1
            if jj < 0:
                break
            idirec = -1

    nprmg -= 1
    if nprmg > istmax:
        return OVERFLOW_PERM
    return nprmg


@njit(nogil=True, cache=True)
def prmx_nb(kmx, nptmx, N_eff, IX0, IXN, kpx, kpxloc, nuterm):
    """Momentum partitions with Fermi statistics -> ``kpx``, ``kpxloc``.

    Transcription of ``qcdf.prmx``.  The two closures the original defined
    *inside* its innermost loop -- recreated ~640 times per run, each slicing a
    Python list to recompute ``sum(nuterm[:mm-1])`` -- are inlined.  The sum is
    kept as an explicit loop rather than a running total: ``nptcl <= 25`` makes
    it free, and a running total is subtly wrong at ``nptcl == 1``, where the
    original's ``mm > 1`` guard returns 0 rather than the accumulated value.

    Returns the number of partitions, or :data:`OVERFLOW_PERM` if ``kpx`` is
    too small; the caller grows it and retries.
    """
    lkpx = kpx.shape[0]
    for i in range(lkpx):
        for j in range(kpx.shape[1]):
            kpx[i, j] = 0
    for a in range(kpxloc.shape[0]):
        for b in range(kpxloc.shape[1]):
            for c in range(kpxloc.shape[2]):
                kpxloc[a, b, c, 0] = 0
                kpxloc[a, b, c, 1] = 0

    numprm = 1
    for ix in range(IX0, IXN + 1):
        indx = ix - 1
        if ix == IXN:
            indx = N_eff
        for nptcl in range(1, nptmx + 1):
            for kptcl in range(nptcl, kmx + 1):
                for i in range(MXNP):
                    nuterm[i] = 1
                kptpr = kptcl + 1
                inmprm = numprm

                idirec = 1
                jj = 0
                while True:
                    # _kmax_prmx(jj+1, ...): kleft = sum(nuterm[:jj])
                    kleft = 0
                    for t in range(jj):
                        kleft += nuterm[t]
                    nrt = nptcl - (jj + 1)
                    if nrt + 1 == 0:
                        km = kptcl - kleft
                    else:
                        km = (kptcl - kleft) // (nrt + 1)
                    if idirec == 1 or nuterm[jj] < km:
                        if idirec != 1:
                            nuterm[jj] += 1
                        jj += 1
                        idirec = 1
                        if jj >= nptcl - 1:
                            # _kleft_prmx(nptcl, ...): sum(nuterm[:nptcl-1])
                            klft = 0
                            for t in range(nptcl - 1):
                                klft += nuterm[t]
                            nuterm[nptcl - 1] = kptcl - klft
                            if indx != 0:
                                ixok = stachk_nb(nptcl, indx, nuterm)
                            else:
                                ixok = True
                            if ixok:
                                if numprm > lkpx:
                                    return OVERFLOW_PERM
                                for j in range(nptcl):
                                    kpx[numprm - 1, j] = nuterm[j]
                                kpxloc[ix - 1, nptcl - 1, kptpr - 1, 0] = inmprm
                                kpxloc[ix - 1, nptcl - 1, kptpr - 1, 1] = numprm
                                numprm += 1
                            jj -= 1
                            idirec = -1
                        else:
                            nuterm[jj] = nuterm[jj - 1]
                    else:
                        nuterm[jj] = 1
                        jj -= 1
                        if jj < 0:
                            break
                        idirec = -1
    return numprm - 1


@njit(nogil=True, cache=True)
def prmm_nb(kmx, kpm, kpmloc):
    """Two-body momentum splits -> ``kpm``, ``kpmloc``.  See ``qcdf.prmm``."""
    for i in range(kpm.shape[0]):
        kpm[i, 0] = 0
        kpm[i, 1] = 0
    for i in range(kpmloc.shape[0]):
        kpmloc[i, 0] = 0
        kpmloc[i, 1] = 0
    nprmm = 0
    for k1 in range(2, kmx + 1):
        kpmloc[k1, 0] = nprmm + 1
        for k2 in range(1, k1):
            nprmm += 1
            if nprmm > kpm.shape[0]:
                return OVERFLOW_PERM
            kpm[nprmm - 1, 0] = k2
            kpm[nprmm - 1, 1] = k1 - k2
        kpmloc[k1, 1] = nprmm
    return nprmm


# ──────────────────────────────────────────────────────────────────────────
# Parity pre-filters (the hoisted oddchk)
# ──────────────────────────────────────────────────────────────────────────

@njit(nogil=True, cache=True)
def _filter_parity_group(ist, ntrms, ngrp, dim, table, width, ok):
    """Rows of ``ist`` whose momenta are all odd, in ascending order.

    ``table`` is ``kpx`` (baryons, ``width = N``) or ``kpm`` (mesons,
    ``width = 2``); the row index is the quotient ``(ist[l,j]-1)//dim`` and the
    flavour is the remainder, so parity does not depend on the flavour half.
    Returns how many rows survived.
    """
    nok = 0
    for l in range(ntrms):
        good = True
        for j in range(ngrp):
            row = (ist[l, j] - 1) // dim
            for k in range(width):
                if table[row, k] % 2 == 0:
                    good = False
                    break
            if not good:
                break
        if good:
            ok[nok] = l
            nok += 1
    return nok


# ──────────────────────────────────────────────────────────────────────────
# BRMSGN
# ──────────────────────────────────────────────────────────────────────────

@njit(nogil=True, cache=True)
def brmsgn_nb(N, NF, nbrb, nbrd, nmes, kbrb, kbrd, kmes,
              LPN, cutoff, K, rmq, iflv,
              kpx, kpxloc, kpm, kpmloc, IX0, IXN,
              mbarfl, mmesfl, ibfdim, imfdim,
              out_mstate, out_mstinf, ns0,
              istb, istd, istm, okb, okd, okm, nuterm,
              nminbf, nmaxbf, nmindf, nmaxdf, nminmf, nmaxmf,
              istate, ichkb, ichkd, fct):
    """States for one momentum sector, appended starting at index ``ns0``.

    Returns the number written, or :data:`OVERFLOW_PERM` / :data:`OVERFLOW_STATES`.
    A transcription of ``qcdf.brmsgn``; see the module docstring for the three
    deviations.
    """
    nnf = N * NF
    lngsta = N * (nbrb + nbrd) + 2 * nmes

    if LPN != 0 and lngsta > LPN:
        return 0

    ixbar = nnf if N % 2 == 0 else 1

    if nbrb > 0:
        ibmn = kpxloc[IX0 - 1, nbrb - 1, kbrb, 0]
        ibmx = kpxloc[IX0 - 1, nbrb - 1, kbrb, 1]
    else:
        ibmn, ibmx = -1, -1
    if nbrd > 0:
        idmn = kpxloc[IX0 - 1, nbrd - 1, kbrd, 0]
        idmx = kpxloc[IX0 - 1, nbrd - 1, kbrd, 1]
    else:
        idmn, idmx = -1, -1
    if nmes > 0:
        immn = kpxloc[IX0 - 1, nmes - 1, kmes, 0]
        immx = kpxloc[IX0 - 1, nmes - 1, kmes, 1]
    else:
        immn, immx = -1, -1

    if ibmn == 0 or idmn == 0 or immn == 0:
        return 0

    ib_lo, ib_hi = (ibmn - 1, ibmx) if ibmn > 0 else (0, 1)
    id_lo, id_hi = (idmn - 1, idmx) if idmn > 0 else (0, 1)
    im_lo, im_hi = (immn - 1, immx) if immn > 0 else (0, 1)

    ns = ns0
    nmax_states = out_mstinf.shape[0]

    for ib in range(ib_lo, ib_hi):
        for id_ in range(id_lo, id_hi):
            for im in range(im_lo, im_hi):
                skip = False
                if nbrb > 0:
                    for jb in range(nbrb):
                        kv = kpx[ib, jb]
                        nm = kpxloc[IXN - 1, N - 1, kv, 0]
                        nx = kpxloc[IXN - 1, N - 1, kv, 1]
                        if nm == 0:
                            skip = True
                            break
                        nminbf[jb] = (nm - 1) * ibfdim + 1
                        nmaxbf[jb] = nx * ibfdim
                if skip:
                    continue

                if nbrd > 0:
                    for jd in range(nbrd):
                        kv = kpx[id_, jd]
                        nm = kpxloc[IXN - 1, N - 1, kv, 0]
                        nx = kpxloc[IXN - 1, N - 1, kv, 1]
                        if nm == 0:
                            skip = True
                            break
                        nmindf[jd] = (nm - 1) * ibfdim + 1
                        nmaxdf[jd] = nx * ibfdim
                if skip:
                    continue

                if nmes > 0:
                    for jm in range(nmes):
                        kv = kpx[im, jm]
                        nm = kpmloc[kv, 0]
                        nx = kpmloc[kv, 1]
                        if nm == 0:
                            skip = True
                            break
                        nminmf[jm] = (nm - 1) * imfdim + 1
                        nmaxmf[jm] = nx * imfdim
                if skip:
                    continue

                if nbrb > 0:
                    ntrmsb = gprbig_nb(nminbf, nmaxbf, nbrb, ixbar, 1,
                                       istb, nuterm)
                    if ntrmsb == OVERFLOW_PERM:
                        return OVERFLOW_PERM
                else:
                    ntrmsb = 1
                if nbrd > 0:
                    ntrmsd = gprbig_nb(nmindf, nmaxdf, nbrd, ixbar, 1,
                                       istd, nuterm)
                    if ntrmsd == OVERFLOW_PERM:
                        return OVERFLOW_PERM
                else:
                    ntrmsd = 1
                if nmes > 0:
                    ntrmsm = gprbig_nb(nminmf, nmaxmf, nmes, nnf, 1,
                                       istm, nuterm)
                    if ntrmsm == OVERFLOW_PERM:
                        return OVERFLOW_PERM
                else:
                    ntrmsm = 1

                # ── hoisted parity: reject whole rows before building states ──
                if nbrb > 0:
                    nb_ok = _filter_parity_group(istb, ntrmsb, nbrb, ibfdim,
                                                 kpx, N, okb)
                else:
                    nb_ok = 1
                    okb[0] = 0
                if nb_ok == 0:
                    continue
                if nbrd > 0:
                    nd_ok = _filter_parity_group(istd, ntrmsd, nbrd, ibfdim,
                                                 kpx, N, okd)
                else:
                    nd_ok = 1
                    okd[0] = 0
                if nd_ok == 0:
                    continue
                if nmes > 0:
                    nm_ok = _filter_parity_group(istm, ntrmsm, nmes, imfdim,
                                                 kpm, 2, okm)
                else:
                    nm_ok = 1
                    okm[0] = 0
                if nm_ok == 0:
                    continue

                for a1 in range(nb_ok):
                    l1 = okb[a1]
                    for a2 in range(nd_ok):
                        l2 = okd[a2]
                        for a3 in range(nm_ok):
                            l3 = okm[a3]

                            if nbrb > 0:
                                for j1 in range(nbrb):
                                    ibk = (istb[l1, j1] - 1) // ibfdim + 1
                                    ibf = istb[l1, j1] - (ibk - 1) * ibfdim
                                    for k1 in range(N):
                                        i1 = N * j1 + k1
                                        istate[0, i1] = 1
                                        istate[1, i1] = 0
                                        istate[2, i1] = kpx[ibk - 1, k1]
                                        istate[3, i1] = mbarfl[ibf - 1, k1]

                            if nbrd > 0:
                                for j2 in range(nbrd):
                                    idk = (istd[l2, j2] - 1) // ibfdim + 1
                                    idf = istd[l2, j2] - (idk - 1) * ibfdim
                                    for k2 in range(N):
                                        i2 = N * nbrb + N * j2 + k2
                                        istate[0, i2] = 0
                                        istate[1, i2] = 1
                                        istate[2, i2] = kpx[idk - 1, k2]
                                        istate[3, i2] = mbarfl[idf - 1, k2]

                            if nmes > 0:
                                for j3 in range(nmes):
                                    imk = (istm[l3, j3] - 1) // imfdim + 1
                                    imf = istm[l3, j3] - (imk - 1) * imfdim
                                    for k3 in range(2):
                                        i3 = N * (nbrb + nbrd) + 2 * j3 + k3
                                        istate[0, i3] = -k3 + 1
                                        istate[1, i3] = k3
                                        istate[2, i3] = kpm[imk - 1, k3]
                                        istate[3, i3] = mmesfl[imf - 1, k3]

                            # Parity already guaranteed by the row filter.

                            nb = 0
                            nd = 0
                            for mch in range(lngsta):
                                mom = (istate[2, mch] + 1) // 2
                                if istate[0, mch] == 1:
                                    ichkb[nb] = istate[3, mch] + (mom - 1) * NF
                                    nb += 1
                                else:
                                    ichkd[nd] = istate[3, mch] + (mom - 1) * NF
                                    nd += 1

                            jxok = True
                            if nb > 0:
                                jxok = stachk_nb(nb, N, ichkb)
                            if jxok and nd > 0:
                                jxok = stachk_nb(nd, N, ichkd)
                            if not jxok:
                                continue

                            if not _pavill_nb(istate, lngsta, cutoff, K, rmq):
                                continue
                            if not _flvchk_nb(istate, lngsta, NF, iflv, fct):
                                continue

                            if ns >= nmax_states:
                                return OVERFLOW_STATES
                            out_mstinf[ns, 0] = 4 * ns + 1
                            out_mstinf[ns, 1] = lngsta
                            out_mstinf[ns, 2] = nmes
                            out_mstinf[ns, 3] = nbrb
                            out_mstinf[ns, 4] = nbrd
                            out_mstinf[ns, 5] = kmes
                            out_mstinf[ns, 6] = kbrb
                            out_mstinf[ns, 7] = kbrd
                            for iw in range(4):
                                for iz in range(lngsta):
                                    out_mstate[4 * ns + iw, iz] = istate[iw, iz]
                            ns += 1

    return ns - ns0
