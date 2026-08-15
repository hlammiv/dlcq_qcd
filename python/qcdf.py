#!/usr/bin/env python3
"""
qcdf.py — Discretized Light-Cone Quantization for 1+1D QCD with multiple flavors.

Translated from Fortran 77 code by Kent Hornbostel (1993).
Python translation with multiprocessing parallelization.

This program computes hadron mass spectra in 1+1 dimensional QCD using DLCQ.
It generates Fock-space states constrained by baryon number, momentum conservation,
Pauli exclusion, flavor quantum numbers, and optional Pauli-Villars cutoffs.
It constructs norm and Hamiltonian matrices with diagrammatic color-singlet
projections, orthonormalizes the basis, and diagonalizes to obtain mass² eigenvalues.

Usage:
    python qcdf.py              # interactive mode (prompts for parameters)
    python qcdf.py input.json   # batch mode from JSON parameter file

Parallelization:
    The O(N²) Hamiltonian matrix construction loop (CLRDIS) is parallelized
    across CPU cores using multiprocessing. Set QCDF_NCPUS env var to control.
"""

import os
import sys
import json
import math
import time
import numpy as np
from scipy.linalg import eigh
from multiprocessing import Pool, cpu_count, shared_memory
from dataclasses import dataclass, field
from typing import Optional
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants / dimension limits (mirroring Fortran PARAMETER-like values)
# ---------------------------------------------------------------------------
# The three K-sensitive caps.  qcdf.f fixes them at NMSTMX = 6902 and
# LKPXMX = 25001, which is what the paper's 2K <= 24 needed; they are raised
# here because nothing else stands between this code and larger K.  Measured
# requirements for the N=3, B=1 baryon (the channel that grows fastest):
#
#     2K     partons    KPX perms    states
#     29        11         18801       4529     <- fits the original caps
#     31        13         32072       7569     <- needs both raised
#     33        13         45013      12471
#     35        13         62290      20353
#
# The cost is ~37 MB of zeroed array per run instead of ~11 MB.  Beyond 2K=35
# raise these again; both overflow paths raise RuntimeError naming the cap, and
# MXNP is not the binding constraint (lpnsub caps the baryon at 13 partons).
NMSTMX   = 25000   # max number of states
NMXFR    = 100000  # 4 * NMSTMX  (state storage)
NMFLMX   = 3       # max number of flavors
MXNP     = 25      # max length of states (particles per state)
LKPXMX   = 80000   # max permutations in KPX
MXKMX    = 100     # max total momentum in KPX
ISTMAX_SM = 6902   # array size for gprmx intermediate arrays
ISTMAX_BG = 9523   # array size for gprbig intermediate arrays
# Colour-contraction capacity for the pure-Python clfact/bredce path.  The
# improved (van de Sande) Hamiltonian reaches this path through
# dlcq.endpoint.pair_colour_weights and needs far more terms than the inherited
# Fortran value at large N: an N=7 baryon overflows it at 2K=7.  Override with
# DLCQ_MXTRMS.  Distinct from MXTRM below, which serves the njit/opt kernels.
MXTRMS   = int(os.environ.get("DLCQ_MXTRMS", "12552"))

# Colour-contraction capacity for the njit/opt kernels.  Distinct from MXTRMS
# above: that one belongs to the pure-Python reference path, this one to the
# hot kernels, and they are NOT interchangeable.  Defined here so both
# qcdf_kernels and qcdf_opt read one value -- they used to declare it
# separately, and the overflow message had to tell you to raise it in *both*.
#
# Cost is small even when large: np.zeros hands back lazy zero pages and the
# kernels only touch rows 0..ntrms, so resident memory tracks the actual term
# count rather than the cap.  Raise via DLCQ_MXTRM for large-N baryons, where
# the epsilon tensor needs N! permutations (8! = 40320 against 7! = 5040).
MXTRM = int(os.environ.get("DLCQ_MXTRM", "200000"))
MAXEP    = 9523    # max epsilon permutations
IVCTMX   = 75      # max eigenvectors to print
NCOL_MAX = 12      # max colors for baryon flavor perms (a baryon has N quarks)
NPRM_BAR = 1001    # max baryon flavor perms
NPRM_MES = 101     # max meson flavor perms
MAXK_CLF = 130     # max momentum in CLFACT

# ---------------------------------------------------------------------------
# Global parameter container
# ---------------------------------------------------------------------------
@dataclass
class Params:
    """Physics and numerical parameters."""
    rlamb: float = 0.0       # coupling λ
    N: int = 3               # number of colors
    NF: int = 1              # number of flavors
    B: int = 0               # baryon number
    K: int = 4               # 2× total LC momentum
    MASS: int = 1            # massive quarks flag
    rmq: np.ndarray = field(default_factory=lambda: np.ones(NMFLMX))
    iflv: np.ndarray = field(default_factory=lambda: np.zeros(NMFLMX, dtype=int))
    cutoff: float = -1.0     # Pauli-Villars cutoff (≤0 = none)
    LPN: int = 0             # particle number limit (0 = none)
    IB0: int = 0             # allow zero-momentum B's
    ID0: int = 0             # allow zero-momentum D's
    IBAB: int = 0            # allow extra bar-antibar pairs
    INTPRT: int = 0          # print intermediate Hamiltonian values
    cbreak: float = 1.0e-8   # charge conjugation breaking parameter

# ---------------------------------------------------------------------------
# State storage
# ---------------------------------------------------------------------------
@dataclass
class StateData:
    """Fock-space state storage."""
    numsta: int = 0
    # mstate[i, j]: state data; i indexes into 4*state_number + row, j is particle
    mstate: np.ndarray = field(default_factory=lambda: np.zeros((NMXFR, MXNP), dtype=int))
    # mstinf[i, j]: state info (location, length, nmes, nbrb, nbrd, kmes, kbrb, kbrd)
    mstinf: np.ndarray = field(default_factory=lambda: np.zeros((NMSTMX, 8), dtype=int))

# ---------------------------------------------------------------------------
# Momentum permutation tables
# ---------------------------------------------------------------------------
@dataclass
class PermTables:
    """KPX, KPXLOC, KPM, KPMLOC tables for momentum distributions."""
    kpx: np.ndarray = field(default_factory=lambda: np.zeros((LKPXMX, MXNP), dtype=int))
    kpxloc: np.ndarray = field(default_factory=lambda: np.zeros((3, MXNP, MXKMX, 2), dtype=int))
    kpm: np.ndarray = field(default_factory=lambda: np.zeros((LKPXMX, 2), dtype=int))
    kpmloc: np.ndarray = field(default_factory=lambda: np.zeros((MXKMX, 2), dtype=int))
    numprm: int = 1
    nprmm: int = 0
    IX0: int = 1
    IXN: int = 2

# ---------------------------------------------------------------------------
# Flavor permutation tables
# ---------------------------------------------------------------------------
@dataclass
class FlavorTables:
    """Baryon and meson flavor permutation tables."""
    mbarfl: np.ndarray = field(default_factory=lambda: np.zeros((NPRM_BAR, NCOL_MAX), dtype=int))
    ibfdim: int = 0
    mmesfl: np.ndarray = field(default_factory=lambda: np.zeros((NPRM_MES, 2), dtype=int))
    imfdim: int = 0

# ---------------------------------------------------------------------------
# Self-energy table
# ---------------------------------------------------------------------------
def compute_selfen(params):
    """Compute self-energies (SUB SLFN)."""
    mxslfn = 100
    selfen = np.zeros(mxslfn)
    selfen[0] = 0.0
    for i in range(0, mxslfn - 3, 2):
        selfen[i + 2] = selfen[i] + 4.0 / ((i + 2) ** 2)
    # Apply color factor to odd-indexed entries (1-based odd = 0-based even)
    cfact = (params.N * params.N - 1.0) / params.N
    for i in range(0, mxslfn - 1, 2):
        selfen[i] *= cfact
    return selfen

# ---------------------------------------------------------------------------
# Small utility functions
# ---------------------------------------------------------------------------
def iodev(i):
    """Return 1 if i is even, -1 if odd (for anti-periodic BCs)."""
    return 1 if (i % 2 == 0) else -1

def nfact(m):
    """Compute m! (factorial)."""
    result = 1
    for j in range(1, m + 1):
        result *= j
    return result

def minmom(nquant, N, NF):
    """
    Given nquant particles, compute minimum momentum they carry.
    Assumes anti-periodic BCs and accounts for flavor.
    """
    nnf = N * NF
    if nnf == 0:
        return 0
    nefbr = nquant // nnf
    nqrem = nquant - nnf * nefbr
    return nnf * nefbr * nefbr + nqrem * (2 * nefbr + 1)

def brack(L, M):
    """
    Compute the square bracket [L, M] used in the Hamiltonian.
    Returns 4/L² if L≠0, M≠0, and L+M=0; else 0.
    """
    if L == 0 or M == 0 or (L + M) != 0:
        return 0.0
    return 4.0 / (float(L) ** 2)

def jdelta(i, j):
    """Kronecker delta."""
    return 1 if i == j else 0

# ---------------------------------------------------------------------------
# State check: Fermi statistics (SUB STACHK)
# ---------------------------------------------------------------------------
def stachk(nptcl, indx, nuterm):
    """
    Check if state satisfies Fermi statistics with index 'indx'.
    Returns True if OK, False if violated.
    nuterm is a list/array of length nptcl.
    """
    mcount = {}
    negmct = {}
    for nn in range(nptcl):
        val = nuterm[nn]
        if val >= 0:
            kp1 = val + 1
            mcount[kp1] = mcount.get(kp1, 0) + 1
            if mcount[kp1] > indx:
                return False
        else:
            negkp1 = -val + 1
            negmct[negkp1] = negmct.get(negkp1, 0) + 1
            if negmct[negkp1] > indx:
                return False
    return True

# ---------------------------------------------------------------------------
# Odd-momentum check (SUB ODDCHK)
# ---------------------------------------------------------------------------
def oddchk(istate, lngsta):
    """Check if all momenta in state are odd. Returns True if all odd."""
    for nd in range(lngsta):
        mom = istate[2, nd]
        if mom % 2 == 0:
            return False
    return True

# ---------------------------------------------------------------------------
# Pauli-Villars cutoff check (SUB PAVILL)
# ---------------------------------------------------------------------------
def pavill(istate, lngsta, params):
    """Check if state satisfies Pauli-Villars cutoff. Returns True if passes."""
    if params.cutoff <= 0:
        return True
    eps = 1.0e-6
    elc = 0.0
    for ic in range(lngsta):
        x = float(istate[2, ic]) / float(params.K)
        if x < eps:
            return False
        rm = params.rmq[istate[3, ic] - 1]  # flavor index is 1-based in Fortran
        elc += rm * rm / x
    return elc < params.cutoff

# ---------------------------------------------------------------------------
# Flavor quantum number check (SUB FLVCHK)
# ---------------------------------------------------------------------------
def flvchk(istate, lngsta, params):
    """Check if state has correct flavor quantum numbers. Returns True if OK."""
    fct = np.zeros(params.NF, dtype=int)
    for ic in range(lngsta):
        mf = istate[3, ic] - 1  # 0-based flavor index
        if istate[0, ic] == 1:
            fct[mf] += 1
        else:
            fct[mf] -= 1
    for id_ in range(params.NF):
        if fct[id_] != params.iflv[id_]:
            return False
    return True

# ---------------------------------------------------------------------------
# Conjugate state (SUB CONJ)
# ---------------------------------------------------------------------------
def conj_state(lst, ist):
    """Conjugate state: reverse order, negate rows 0 and 1."""
    istc = np.zeros_like(ist)
    for i1 in range(lst):
        j1 = lst - 1 - i1
        for i2 in range(4):
            if i2 <= 1:
                istc[i2, j1] = -ist[i2, i1]
            else:
                istc[i2, j1] = ist[i2, i1]
    return istc

# ---------------------------------------------------------------------------
# Permutation sign (SUB PSIGN)
# ---------------------------------------------------------------------------
def psign(lngth, iperm):
    """Compute sign of antisymmetric permutation."""
    if lngth <= 1:
        return 1
    # Initialize so that sgn(1,2,3,...,n) = +1
    i0 = 0
    for l1 in range(1, lngth):
        for l2 in range(l1):
            i0 += 1
    iexp = 1 if ((-1) ** i0) < 0 else 0
    for i1 in range(1, lngth):
        for i2 in range(i1):
            if iperm[i1] > iperm[i2]:
                iexp += 1
    return (-1) ** iexp

# ---------------------------------------------------------------------------
# LPNSUB: compute max quarks/antiquarks for given N, NF, B, K
# ---------------------------------------------------------------------------
def lpnsub(N, NF, B, K):
    """
    Given N, NF, B, K compute max number of b-adj (LNB) and d-adj (LND)
    that can appear in a state. For anti-periodic BCs.
    Returns (LNB, LND).
    """
    nbmin = N * B
    nbmax = 25
    for nsubb in range(nbmin, nbmax + 2):
        nsubd = nsubb - N * B
        if nsubb == nbmax + 1:
            raise RuntimeError("NSUBB exceeds NBMAX in lpnsub")
        kbmin = minmom(nsubb, N, NF)
        kdmin = minmom(nsubd, N, NF)
        kmin = kdmin + kbmin
        if kmin > K:
            lnb = nsubb - 1
            lnd = lnb - N * B
            if lnb < N * B:
                raise RuntimeError("K insufficient in lpnsub")
            return lnb, lnd
    raise RuntimeError("Failed to find bounds in lpnsub")

# ---------------------------------------------------------------------------
# PRMX: generate momentum distributions with fermionic statistics
# ---------------------------------------------------------------------------
def prmx(kmx, nptmx, N_eff, perm):
    """
    Generate momentum distributions with fermionic statistics.
    N_eff is the effective fermionic index (N*NF or N depending on context).
    Populates perm.kpx, perm.kpxloc, perm.numprm.
    """
    perm.kpx[:] = 0
    perm.kpxloc[:] = 0
    numprm = 1

    for ix in range(perm.IX0, perm.IXN + 1):
        indx = ix - 1
        if ix == perm.IXN:
            indx = N_eff
        for nptcl in range(1, nptmx + 1):
            # No zero momenta allowed (IB0=0, ID0=0)
            kptcli = nptcl
            kptclf = kmx
            for kptcl in range(kptcli, kptclf + 1):
                nuterm = [1] * MXNP
                kptpr = kptcl + 1
                inmprm = numprm

                def _kmax_prmx(mm, nptcl_l, kptcl_l, nuterm_l):
                    nrt = nptcl_l - mm
                    kleft = sum(nuterm_l[:mm - 1]) if mm > 1 else 0
                    if nrt + 1 == 0:
                        return kptcl_l - kleft
                    return (kptcl_l - kleft) // (nrt + 1)

                def _kleft_prmx(mm, nuterm_l):
                    return sum(nuterm_l[:mm - 1]) if mm > 1 else 0

                # Iterative permutation generator (mirrors Fortran logic)
                idirec = 1
                jj = 0  # 0-based
                while True:
                    km = _kmax_prmx(jj + 1, nptcl, kptcl, nuterm)
                    if idirec == 1 or nuterm[jj] < km:
                        if idirec != 1:
                            nuterm[jj] += 1
                        jj += 1
                        idirec = 1
                        if jj >= nptcl - 1:
                            klft = _kleft_prmx(nptcl, nuterm)
                            nuterm[nptcl - 1] = kptcl - klft
                            if indx != 0:
                                ixok = stachk(nptcl, indx, nuterm[:nptcl])
                            else:
                                ixok = True
                            if ixok:
                                for j in range(nptcl):
                                    perm.kpx[numprm - 1, j] = nuterm[j]
                                perm.kpxloc[ix - 1, nptcl - 1, kptpr - 1, 0] = inmprm
                                perm.kpxloc[ix - 1, nptcl - 1, kptpr - 1, 1] = numprm
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

    perm.numprm = numprm - 1
    if perm.numprm >= LKPXMX:
        raise RuntimeError(
            f"Overflow in KPX: {perm.numprm} momentum permutations needed, "
            f"LKPXMX = {LKPXMX}. Raise LKPXMX in qcdf.py.")

# ---------------------------------------------------------------------------
# PRMM: generate meson momentum distributions
# ---------------------------------------------------------------------------
def prmm(kmx, perm):
    """Generate momentum distributions for individual mesons."""
    perm.kpm[:] = 0
    perm.kpmloc[:] = 0
    nprmm = 0
    # No zero momenta allowed
    for k1 in range(2, kmx + 1):
        perm.kpmloc[k1, 0] = nprmm + 1
        for k2 in range(1, k1):
            nprmm += 1
            perm.kpm[nprmm - 1, 0] = k2
            perm.kpm[nprmm - 1, 1] = k1 - k2
        perm.kpmloc[k1, 1] = nprmm
    perm.nprmm = nprmm
    if nprmm >= LKPXMX:
        raise RuntimeError(
            f"Overflow in KPM: {nprmm} meson permutations needed, "
            f"LKPXMX = {LKPXMX}. Raise LKPXMX in qcdf.py.")

# ---------------------------------------------------------------------------
# IBARFL: generate baryon flavor permutations
# ---------------------------------------------------------------------------
def ibarfl(params, flavtab):
    """Generate flavor permutations for baryons."""
    N = params.N
    NF = params.NF
    # A baryon carries N quarks, so mbarfl's second axis must reach N.  This
    # runs for meson jobs too, so the cap binds at every B.
    if N > NCOL_MAX:
        raise RuntimeError(
            f"Overflow in MBARFL: N = {N} colors needed, NCOL_MAX = {NCOL_MAX}. "
            f"Raise NCOL_MAX in qcdf.py.")
    flavtab.mbarfl[:] = 0

    # Generate all NF^N permutations (no Fermi restriction on flavor alone)
    nuterm = [1] * N
    nprmg = 0

    def _recurse(depth):
        nonlocal nprmg
        if depth == N:
            for j in range(N):
                flavtab.mbarfl[nprmg, j] = nuterm[j]
            nprmg += 1
            if nprmg >= NPRM_BAR:
                raise RuntimeError("nprmg exceeds nprmmx in ibarfl")
            return
        for val in range(1, NF + 1):
            nuterm[depth] = val
            _recurse(depth + 1)

    _recurse(0)
    flavtab.ibfdim = nprmg

# ---------------------------------------------------------------------------
# IMESFL: generate meson flavor permutations
# ---------------------------------------------------------------------------
def imesfl(params, flavtab):
    """Generate flavor permutations for mesons."""
    NF = params.NF
    nprm = 0
    for i in range(1, NF + 1):
        for j in range(1, NF + 1):
            flavtab.mmesfl[nprm, 0] = i
            flavtab.mmesfl[nprm, 1] = j
            nprm += 1
            if nprm >= NPRM_MES:
                raise RuntimeError("nprm exceeds nprmmx in imesfl")
    flavtab.imfdim = nprm

# ---------------------------------------------------------------------------
# GPRBIG: permutation generator (general, large arrays)
# ---------------------------------------------------------------------------
def gprbig(nmin, nmax, lngth, indx, idbl):
    """
    Generate permutations of indices nmin[0..lngth-1] to nmax[0..lngth-1]
    respecting Fermi statistics with index 'indx' and double-counting
    prevention if idbl=1.
    Returns (nprmg, ist) where ist is list of tuples.
    """
    results = []
    nuterm = list(nmin[:lngth])

    def _gen(pos):
        if pos == lngth:
            # Check Fermi statistics
            if indx != 0:
                if not stachk(lngth, indx, nuterm[:lngth]):
                    return
            results.append(tuple(nuterm[:lngth]))
            return
        start = nmin[pos]
        if idbl == 1 and pos > 0 and nmin[pos] == nmin[pos - 1]:
            start = max(start, nuterm[pos - 1])
        for val in range(start, nmax[pos] + 1):
            nuterm[pos] = val
            _gen(pos + 1)

    _gen(0)
    return len(results), results

# ---------------------------------------------------------------------------
# Generalized permutation for state generation (matches Fortran GPRBIG logic)
# ---------------------------------------------------------------------------
def gprbig_array(nmin, nmax, lngth, indx, idbl, istmax):
    """
    Same as gprbig but returns numpy array matching Fortran IST(istmax, 25).
    Uses iterative approach matching Fortran control flow exactly.
    """
    ist = np.zeros((istmax, MXNP), dtype=int)
    nuterm = np.zeros(MXNP, dtype=int)
    nprmg = 1
    idirec = 1
    jj = 0  # 0-based

    for ll in range(lngth):
        nuterm[ll] = nmin[ll]

    while True:
        km = nmax[jj] if jj < lngth else 0
        if idirec == 1 or nuterm[jj] < km:
            if idirec != 1:
                nuterm[jj] += 1
            # Prevent double counting
            if (idbl == 1 and jj > 0 and nmin[jj] == nmin[jj - 1]
                    and (nuterm[jj] - 1) < nuterm[jj - 1]):
                nuterm[jj] = nuterm[jj - 1]
            jj += 1
            idirec = 1
            if jj >= lngth - 1:
                # Prevent double counting for last element
                if (idbl == 1 and lngth > 1 and nmin[lngth - 1] == nmin[lngth - 2]
                        and (nuterm[lngth - 1] - 1) < nuterm[lngth - 2]):
                    nuterm[lngth - 1] = nuterm[lngth - 2]
                for ii in range(nuterm[lngth - 1], nmax[lngth - 1] + 1):
                    nuterm[lngth - 1] = ii
                    if indx != 0:
                        ixok = stachk(lngth, indx, nuterm[:lngth].tolist())
                    else:
                        ixok = True
                    if ixok:
                        for j in range(lngth):
                            ist[nprmg - 1, j] = nuterm[j]
                        nprmg += 1
                        if nprmg > istmax:
                            raise RuntimeError(f"nprmg exceeds istmax ({istmax}) in gprbig_array")
                nuterm[lngth - 1] = nmin[lngth - 1]
                jj = lngth - 2
                if jj < 0:
                    break
                idirec = -1
            # else: nuterm[jj] = nuterm[jj-1] already handled by loop
        else:
            nuterm[jj] = nmin[jj]
            jj -= 1
            if jj < 0:
                break
            idirec = -1

    nprmg -= 1
    return nprmg, ist


# ===================================================================
# BRMSGN: given baryon/meson counts and momenta, generate states
# ===================================================================
def brmsgn(params, states, perm, flavtab, N, NF, nbrb, nbrd, nmes, kbrb, kbrd, kmes):
    """
    Given number and momentum of baryons, antibaryons, and mesons,
    compute possible Fock states and store them in states.
    """
    nnf = N * NF
    lngsta = N * (nbrb + nbrd) + 2 * nmes
    ibfdim = flavtab.ibfdim
    imfdim = flavtab.imfdim

    # Particle number limit check
    if params.LPN != 0 and lngsta > params.LPN:
        return

    # Determine baryon exchange symmetry
    ixbar = nnf if ((-1) ** N > 0) else 1

    # Get permutation ranges from KPX tables
    if nbrb > 0:
        ibmn = perm.kpxloc[perm.IX0 - 1, nbrb - 1, kbrb, 0]
        ibmx = perm.kpxloc[perm.IX0 - 1, nbrb - 1, kbrb, 1]
    else:
        ibmn, ibmx = -1, -1

    if nbrd > 0:
        idmn = perm.kpxloc[perm.IX0 - 1, nbrd - 1, kbrd, 0]
        idmx = perm.kpxloc[perm.IX0 - 1, nbrd - 1, kbrd, 1]
    else:
        idmn, idmx = -1, -1

    if nmes > 0:
        immn = perm.kpxloc[perm.IX0 - 1, nmes - 1, kmes, 0]
        immx = perm.kpxloc[perm.IX0 - 1, nmes - 1, kmes, 1]
    else:
        immn, immx = -1, -1

    if ibmn == 0 or idmn == 0 or immn == 0:
        return

    # Ranges for iteration (convert to 0-based, handle -1 sentinel)
    ib_range = range(ibmn - 1, ibmx) if ibmn > 0 else range(1)
    id_range = range(idmn - 1, idmx) if idmn > 0 else range(1)
    im_range = range(immn - 1, immx) if immn > 0 else range(1)

    for ib in ib_range:
        for id_ in id_range:
            for im in im_range:
                # Get momentum sub-distributions and their flavor ranges
                nminbf = np.zeros(MXNP, dtype=int)
                nmaxbf = np.zeros(MXNP, dtype=int)
                nmindf = np.zeros(MXNP, dtype=int)
                nmaxdf = np.zeros(MXNP, dtype=int)
                nminmf = np.zeros(MXNP, dtype=int)
                nmaxmf = np.zeros(MXNP, dtype=int)

                skip = False
                if nbrb > 0:
                    for jb in range(nbrb):
                        kpx_val = perm.kpx[ib, jb]
                        nm = perm.kpxloc[perm.IXN - 1, N - 1, kpx_val, 0]
                        nx = perm.kpxloc[perm.IXN - 1, N - 1, kpx_val, 1]
                        if nm == 0:
                            skip = True
                            break
                        nminbf[jb] = (nm - 1) * ibfdim + 1
                        nmaxbf[jb] = nx * ibfdim

                if skip:
                    continue

                if nbrd > 0:
                    for jd in range(nbrd):
                        kpx_val = perm.kpx[id_, jd]
                        nm = perm.kpxloc[perm.IXN - 1, N - 1, kpx_val, 0]
                        nx = perm.kpxloc[perm.IXN - 1, N - 1, kpx_val, 1]
                        if nm == 0:
                            skip = True
                            break
                        nmindf[jd] = (nm - 1) * ibfdim + 1
                        nmaxdf[jd] = nx * ibfdim

                if skip:
                    continue

                if nmes > 0:
                    for jm in range(nmes):
                        kpx_val = perm.kpx[im, jm]
                        nm = perm.kpmloc[kpx_val, 0]
                        nx = perm.kpmloc[kpx_val, 1]
                        if nm == 0:
                            skip = True
                            break
                        nminmf[jm] = (nm - 1) * imfdim + 1
                        nmaxmf[jm] = nx * imfdim

                if skip:
                    continue

                # Generate baryon momentum+flavor combinations
                if nbrb > 0:
                    ntrmsb, istb = gprbig_array(nminbf, nmaxbf, nbrb, ixbar, 1, ISTMAX_BG)
                else:
                    ntrmsb = 1
                    istb = np.zeros((ISTMAX_BG, MXNP), dtype=int)

                if nbrd > 0:
                    ntrmsd, istd = gprbig_array(nmindf, nmaxdf, nbrd, ixbar, 1, ISTMAX_BG)
                else:
                    ntrmsd = 1
                    istd = np.zeros((ISTMAX_BG, MXNP), dtype=int)

                if nmes > 0:
                    ntrmsm, istm = gprbig_array(nminmf, nmaxmf, nmes, nnf, 1, ISTMAX_BG)
                else:
                    ntrmsm = 1
                    istm = np.zeros((ISTMAX_BG, MXNP), dtype=int)

                # Iterate over all combinations
                for l1 in range(ntrmsb):
                    for l2 in range(ntrmsd):
                        for l3 in range(ntrmsm):
                            istate = np.zeros((4, MXNP), dtype=int)

                            # Fill baryons
                            if nbrb > 0:
                                for j1 in range(nbrb):
                                    ibk = (istb[l1, j1] - 1) // ibfdim + 1
                                    ibf = istb[l1, j1] - (ibk - 1) * ibfdim
                                    for k1 in range(N):
                                        i1 = N * j1 + k1
                                        istate[0, i1] = 1
                                        istate[1, i1] = 0
                                        istate[2, i1] = perm.kpx[ibk - 1, k1]
                                        istate[3, i1] = flavtab.mbarfl[ibf - 1, k1]

                            # Fill antibaryons
                            if nbrd > 0:
                                for j2 in range(nbrd):
                                    idk = (istd[l2, j2] - 1) // ibfdim + 1
                                    idf = istd[l2, j2] - (idk - 1) * ibfdim
                                    for k2 in range(N):
                                        i2 = N * nbrb + N * j2 + k2
                                        istate[0, i2] = 0
                                        istate[1, i2] = 1
                                        istate[2, i2] = perm.kpx[idk - 1, k2]
                                        istate[3, i2] = flavtab.mbarfl[idf - 1, k2]

                            # Fill mesons
                            if nmes > 0:
                                for j3 in range(nmes):
                                    imk = (istm[l3, j3] - 1) // imfdim + 1
                                    imf = istm[l3, j3] - (imk - 1) * imfdim
                                    for k3 in range(2):
                                        i3 = N * (nbrb + nbrd) + 2 * j3 + k3
                                        istate[0, i3] = -k3 + 1  # 1 for quark, 0 for antiquark
                                        istate[1, i3] = k3        # 0 for quark, 1 for antiquark
                                        istate[2, i3] = perm.kpm[imk - 1, k3]
                                        istate[3, i3] = flavtab.mmesfl[imf - 1, k3]

                            # --- Validation checks ---
                            # Odd momenta check
                            if not oddchk(istate, lngsta):
                                continue

                            # Fermi statistics check (N quarks/antiquarks max per momentum+flavor)
                            numbs = N * nbrb + nmes
                            numds = N * nbrd + nmes
                            ichkb = []
                            ichkd = []
                            for mch in range(lngsta):
                                if istate[0, mch] == 1:
                                    momb = (istate[2, mch] + 1) // 2
                                    ichkb.append(istate[3, mch] + (momb - 1) * params.NF)
                                else:
                                    momd = (istate[2, mch] + 1) // 2
                                    ichkd.append(istate[3, mch] + (momd - 1) * params.NF)

                            jxok = True
                            if len(ichkb) > 0:
                                jxok = stachk(len(ichkb), N, ichkb)
                            if jxok and len(ichkd) > 0:
                                jxok = stachk(len(ichkd), N, ichkd)

                            if not jxok:
                                continue

                            # Pauli-Villars cutoff
                            if not pavill(istate, lngsta, params):
                                continue

                            # Flavor quantum numbers
                            if not flvchk(istate, lngsta, params):
                                continue

                            # Store the state
                            ns = states.numsta
                            if ns >= NMSTMX:
                                raise RuntimeError(
                                    f"numsta exceeds NMSTMX = {NMSTMX}. "
                                    f"Raise NMSTMX (and NMXFR = 4*NMSTMX) "
                                    f"in qcdf.py.")
                            states.mstinf[ns, 0] = 4 * ns + 1  # 1-based location
                            states.mstinf[ns, 1] = lngsta
                            states.mstinf[ns, 2] = nmes
                            states.mstinf[ns, 3] = nbrb
                            states.mstinf[ns, 4] = nbrd
                            states.mstinf[ns, 5] = kmes
                            states.mstinf[ns, 6] = kbrb
                            states.mstinf[ns, 7] = kbrd

                            for iw in range(4):
                                for iz in range(lngsta):
                                    states.mstate[4 * ns + iw, iz] = istate[iw, iz]

                            states.numsta += 1

# ===================================================================
# QCDSTA: QCD state generator (main)
# ===================================================================
def qcdsta(params, states, perm, flavtab):
    """Generate all valid Fock states."""
    N, NF, B, K = params.N, params.NF, params.B, params.K
    states.numsta = 0

    lnb, lnd = lpnsub(N, NF, B, K)
    log.info(f"  LPNSUB: LNB={lnb}, LND={lnd}")

    # Compute max particle number
    if params.LPN != 0 and params.LPN < lnb:
        nptmx = params.LPN
    else:
        nptmx = lnb

    if params.LPN > (lnb + lnd):
        params.LPN = lnb + lnd

    # Generate momentum permutation tables
    nnf = N * NF
    prmx(K, nptmx, nnf, perm)
    log.info("  Called PRMX")

    prmm(K, perm)
    log.info("  Called PRMM")

    # Generate flavor permutation tables
    ibarfl(params, flavtab)
    log.info(f"  Called IBARFL (ibfdim={flavtab.ibfdim})")

    imesfl(params, flavtab)
    log.info(f"  Called IMESFL (imfdim={flavtab.imfdim})")

    # Iterate over Fock-space sectors
    lmbf = lnb - N * B
    for lmb in range(lmbf + 1):
        nxbmx = lmb // N if params.IBAB != 0 else 0
        for nxb in range(nxbmx + 1):
            nbar = B + 2 * nxb
            nmes = lmb - N * nxb

            nbb = N * (B + nxb)
            ndb = N * nxb
            kbmn = minmom(nbb, N, NF)
            kdmn = minmom(ndb, N, NF)
            kmmn = minmom(nmes, N, NF)
            kbarmn = kbmn + kdmn
            kmesmn = 2 * kmmn

            kbarmx = K - kmesmn
            kmesmx = K - kbarmn

            if nbar == 0:
                kmesmn, kmesmx = K, K
                kbarmn, kbarmx = 0, 0
            if nmes == 0:
                kmesmn, kmesmx = 0, 0
                kbarmn, kbarmx = K, K

            if ((nmes != 0 or nbar != 0) and
                    ((N * B + lmb) <= lnb and lmb <= lnd)):

                nbrb = B + nxb
                nbrd = nxb
                noe = iodev(N)
                nbrboe = iodev(nbrb)
                nbrdoe = iodev(nbrd)

                for kbar in range(kbarmn, kbarmx + 1):
                    kmes_val = K - kbar
                    kbrbmn = minmom(N * nbrb, N, NF)
                    kbrdmn = minmom(N * nbrd, N, NF)
                    kbrbmx = kbar - kbrdmn
                    kbrdmx = kbar - kbrbmn

                    if nbrd == 0:
                        kbrbmn, kbrbmx = kbar, kbar
                        kbrdmn, kbrdmx = 0, 0

                    for kbrb in range(kbrbmn, kbrbmx + 1):
                        kbrd = kbar - kbrb

                        # APBC check
                        kmesoe = iodev(kmes_val)
                        kbrboe = iodev(kbrb)
                        kbrdoe = iodev(kbrd)

                        if kmesoe == 1:
                            if ((noe == 1 and kbrboe == 1 and kbrdoe == 1) or
                                    (noe == -1 and kbrboe == nbrboe and kbrdoe == nbrdoe)):
                                brmsgn(params, states, perm, flavtab,
                                        N, NF, nbrb, nbrd, nmes, kbrb, kbrd, kmes_val)

    log.info(f"  QCDSTA complete: {states.numsta} states generated")


# ===================================================================
# COLOR FACTOR ENGINE
# ===================================================================
def clfact_compute(ic1, ic2, ic3, ic4, mx, llt, lrt, nops, params, selfen):
    """
    Evaluate color factor for Hamiltonian matrix element.
    This is a simplified but faithful translation of the Fortran CLFACT
    + CLSUM + CNTRCT + BREDCE + NWTERM + EPSB subroutines.

    Operates on the combined operator/state matrix MX.
    """
    N, NF, B, K = params.N, params.NF, params.B, params.K
    lng = nops + llt + lrt
    maxfl = NMFLMX

    # Check momentum balance
    my = mx[:, :lng].copy()
    mochk = np.zeros((2, 2, K + 2, NF + 1), dtype=int)

    for i3 in range(lng):
        for k3 in range(2):
            mom_idx = mx[2, i3]
            flv_idx = mx[3, i3]
            mochk[k3, 0, mom_idx, flv_idx] += mx[k3, i3]
            ri = lng - 1 - i3
            mochk[k3, 1, mx[2, ri], mx[3, ri]] += mx[k3, ri]
            # Quick kill check
            if (mochk[k3, 0, mom_idx, flv_idx] < 0 or
                    mochk[k3, 0, mom_idx, flv_idx] > N or
                    mochk[k3, 1, mx[2, ri], mx[3, ri]] > 0 or
                    mochk[k3, 1, mx[2, ri], mx[3, ri]] < -N):
                return 0.0

    # --- Build initial delta functions by contracting operators ---
    idel0 = np.zeros((MXTRMS, lng), dtype=int)
    resl0 = np.zeros(MXTRMS)
    resl0[0] = 1.0

    ipow = 0
    for ja in range(1, lng):
        if (my[0, ja] + my[1, ja]) != -1:
            continue
        for ka in range(ja - 1, -1, -1):
            if (my[0, ka] == -my[0, ja] and my[1, ka] == -my[1, ja] and
                    my[2, ka] == my[2, ja] and my[3, ka] == my[3, ja]):
                idel0[0, ja] = ka
                idel0[0, ka] = ja
                my[0, ka] = 0
                my[1, ka] = 0
                my[0, ja] = 0
                my[1, ja] = 0
                break
            else:
                if my[0, ka] != 0 or my[1, ka] != 0:
                    ipow += 1

    resl0[0] = float((-1) ** ipow)

    # --- Build baryon index linking ---
    lnkb = np.zeros(lng, dtype=int)
    icycl = np.ones(lng, dtype=int)
    if B != 0:
        for nbs in range(B):
            for inb in range(N * nbs, N * (nbs + 1) - 1):
                lnkb[inb] = inb + 1
            lnkb[N * (nbs + 1) - 1] = N * nbs

            for jnb in range(lng - N * nbs - 1, lng - N * (nbs + 1), -1):
                lnkb[jnb] = jnb - 1
            lnkb[lng - N * (nbs + 1)] = lng - N * nbs - 1

        if (-1) ** N > 0:
            for nbs in range(B):
                ii = 1
                for inb in range(N * nbs, N * (nbs + 1)):
                    icycl[inb] = ii
                    ii *= -1
                jj = 1
                for jnb in range(lng - N * nbs - 1, lng - N * (nbs + 1) - 1, -1):
                    icycl[jnb] = jj
                    jj *= -1

    # --- Build permutations for degenerate momenta ---
    ncr = np.zeros((2, lng // 2 + 1), dtype=int)
    nperm_arr = np.zeros((lng // 2 + 1, lng), dtype=int)
    lbcr = 0
    ldcr = 0
    ntrms = 1
    lprm = 0
    lpo = 1

    for i in range(lng):
        if mx[1, i] == 0:
            # B operator
            if mx[0, i] == -1:
                lprm += 1
                lpchk = False
                jfstb = -1
                for j in range(lbcr - 1, -1, -1):
                    if not lpchk:
                        if (mx[2, ncr[0, j]] == mx[2, i] and
                                mx[3, ncr[0, j]] == mx[3, i]):
                            jfstb = j
                            for l in range(j - 1, -1, -1):
                                if (mx[2, ncr[0, l]] == mx[2, ncr[0, j]] and
                                        mx[3, ncr[0, l]] == mx[3, ncr[0, j]]):
                                    nperm_arr[lprm - 1, 0] += 1
                                    mp = 2 * nperm_arr[lprm - 1, 0]
                                    nperm_arr[lprm - 1, mp - 1] = ncr[0, j]
                                    nperm_arr[lprm - 1, mp] = ncr[0, l]
                                    lpchk = True
                            break
                if not lpchk:
                    lprm -= 1
                if jfstb >= 0:
                    for m in range(jfstb, lbcr - 1):
                        ncr[0, m] = ncr[0, m + 1]
                    ncr[0, lbcr - 1] = 0
                    lbcr -= 1
            else:
                ncr[0, lbcr] = i
                lbcr += 1
        else:
            # D operator
            if mx[1, i] == -1:
                lprm += 1
                lpchk = False
                jfstd = -1
                for j in range(ldcr - 1, -1, -1):
                    if not lpchk:
                        if (mx[2, ncr[1, j]] == mx[2, i] and
                                mx[3, ncr[1, j]] == mx[3, i]):
                            jfstd = j
                            for l in range(j - 1, -1, -1):
                                if (mx[2, ncr[1, l]] == mx[2, ncr[1, j]] and
                                        mx[3, ncr[1, l]] == mx[3, ncr[1, j]]):
                                    nperm_arr[lprm - 1, 0] += 1
                                    mp = 2 * nperm_arr[lprm - 1, 0]
                                    nperm_arr[lprm - 1, mp - 1] = ncr[1, j]
                                    nperm_arr[lprm - 1, mp] = ncr[1, l]
                                    lpchk = True
                            break
                if not lpchk:
                    lprm -= 1
                if jfstd >= 0:
                    for m in range(jfstd, ldcr - 1):
                        ncr[1, m] = ncr[1, m + 1]
                    ncr[1, ldcr - 1] = 0
                    ldcr -= 1
            else:
                ncr[1, ldcr] = i
                ldcr += 1

    # --- Generate new terms from permutations ---
    if lprm > 0:
        for id_ in range(lprm - 1, -1, -1):
            nsbar = 1
            for jd in range(nperm_arr[id_, 0]):
                kd = 2 * (jd + 1)
                ip1 = nperm_arr[id_, kd - 1]
                ip2 = nperm_arr[id_, kd]

                # Check if permuted indices are in same baryon
                ismb = False
                if B > 0:
                    nxtb = lnkb[ip1]
                    if nxtb != 0 or (ip1 < N * B and ip1 >= 0):
                        for _ in range(N - 1):
                            if ip2 == nxtb:
                                ismb = True
                                break
                            nxtb = lnkb[nxtb]

                if ismb:
                    nsbar += 1
                else:
                    # Create new terms by permuting
                    for j in range(ntrms):
                        lpo += 1
                        if lpo > idel0.shape[0]:
                            raise RuntimeError("NTRMS exceeded MXTRMS in clfact")
                        idel0[lpo - 1, :lng] = idel0[j, :lng]
                        resl0[lpo - 1] = -resl0[j]
                        idel0[lpo - 1, ip1] = idel0[j, ip2]
                        idel0[lpo - 1, ip2] = idel0[j, ip1]
                        if idel0[j, ip2] < lng:
                            idel0[lpo - 1, idel0[j, ip2]] = idel0[j, idel0[j, ip1]]
                        if idel0[j, ip1] < lng:
                            idel0[lpo - 1, idel0[j, ip1]] = idel0[j, idel0[j, ip2]]

            if nsbar > 1:
                for imlt in range(ntrms):
                    resl0[imlt] *= float(nsbar)
            ntrms *= (nperm_arr[id_, 0] + 2 - nsbar)
            if ntrms > idel0.shape[0]:
                raise RuntimeError("NTRMS exceeded MXTRMS in clfact")

    ntrms0 = ntrms

    # --- Sum over color factor terms ---
    nmsr = (lrt - N * B) // 2
    nmsl = (llt - N * B) // 2

    nctf = 2 if nops == 4 else 1
    clfact_val = 0.0

    for nct in range(1, nctf + 1):
        idelt = np.zeros((MXTRMS, lng), dtype=int)
        reslt = np.zeros(MXTRMS)

        # Copy delta functions for upper indices
        for j in range(lng):
            if (mx[0, j] - mx[1, j]) == 1:
                for icl in range(ntrms):
                    idelt[icl, j] = idel0[icl, j]

        for m1 in range(ntrms):
            reslt[m1] = resl0[m1]

        # Fill in color vertex deltas
        if nops == 4:
            if nct == 1:
                for jc in range(ntrms):
                    idelt[jc, lrt + ic1] = lrt + ic2
                    idelt[jc, lrt + ic3] = lrt + ic4
                    reslt[jc] *= 0.5
            else:
                for jc in range(ntrms):
                    idelt[jc, lrt + ic1] = lrt + ic4
                    idelt[jc, lrt + ic3] = lrt + ic2
                    reslt[jc] *= -0.5 / float(N)
        elif nops == 2:
            for jc in range(ntrms):
                idelt[jc, lrt + ic1] = lrt + ic2

        # Fill in meson deltas
        if nmsr > 0:
            for il in range(ntrms):
                for im in range(N * B + 1, lrt, 2):
                    idelt[il, im] = im - 1

        if nmsl > 0:
            for il in range(ntrms):
                for in_ in range(lng - N * B - 1, lng - llt, -2):
                    idelt[il, in_] = in_ - 1

        # --- Contract deltas (SUB CNTRCT) ---
        for nt in range(ntrms):
            for lt in range(lng):
                icurr = lt
                while True:
                    inxt = idelt[nt, icurr]
                    if inxt == 0:
                        break
                    inxt2 = idelt[nt, inxt]
                    if inxt2 == 0:
                        break
                    idelt[nt, icurr] = inxt2
                    idelt[nt, inxt] = 0
                    if inxt2 == icurr:
                        idelt[nt, icurr] = 0
                        reslt[nt] *= float(N)
                        break

        # --- Baryon epsilon reduction (SUB BREDCE) ---
        # Fortran: EPSB called only for B≠0 AND N>2 (generates permutation table)
        #          BREDCE called for ALL B≠0 (reduces epsilon pairs to deltas)
        if B != 0:
            # Generate baryon epsilon perms (only needed for N>2)
            if N > 2:
                nmin_ep = np.ones(MXNP, dtype=int)
                nmax_ep = np.ones(MXNP, dtype=int) * (N - 1)
                nprms_ep, ibrpm = gprbig_array(nmin_ep, nmax_ep, N - 1, 1, 0, MAXEP)
                # Assign signs
                for is1 in range(nprms_ep):
                    iprm = ibrpm[is1, :N - 1].tolist()
                    ibrpm[is1, N - 1] = psign(N - 1, iprm)
            else:
                # For N≤2: no extra permutation terms (Fortran: NPRMS=0)
                nprms_ep = 0
                ibrpm = None

            mxiter = 2 * B + 1
            for _ in range(mxiter):
                # Contract (SUB CNTRCT)
                for nt in range(ntrms):
                    for lt in range(lng):
                        icurr = lt
                        while True:
                            inxt = idelt[nt, icurr]
                            if inxt == 0:
                                break
                            inxt2 = idelt[nt, inxt]
                            if inxt2 == 0:
                                break
                            idelt[nt, icurr] = inxt2
                            idelt[nt, inxt] = 0
                            if inxt2 == icurr:
                                idelt[nt, icurr] = 0
                                reslt[nt] *= float(N)
                                break

                # Reduce epsilon pairs (SUB BREDCE)
                ntpr = ntrms
                has_terms = False
                for nt in range(ntrms):
                    for lt in range(lng):
                        if idelt[nt, lt] != 0:
                            has_terms = True
                            inb1 = lt
                            inb2 = idelt[nt, lt]
                            reslt[nt] *= icycl[inb1] * icycl[inb2]
                            nb1 = []
                            nb2 = []
                            curr1, curr2 = inb1, inb2
                            for _ in range(N - 1):
                                nb1.append(lnkb[curr1])
                                nb2.append(lnkb[curr2])
                                curr1 = lnkb[curr1]
                                curr2 = lnkb[curr2]
                            idelt[nt, lt] = 0
                            for j1 in range(N - 1):
                                idelt[nt, nb2[j1]] = nb1[j1]
                            # Extra permutation terms (only for N>2)
                            for j2 in range(1, nprms_ep):
                                ntpr += 1
                                if ntpr > idelt.shape[0]:
                                    raise RuntimeError("NTPR exceeded MXTRMS in bredce")
                                reslt[ntpr - 1] = reslt[nt] * float(ibrpm[j2, N - 1])
                                idelt[ntpr - 1, :lng] = idelt[nt, :lng]
                                for j4 in range(N - 1):
                                    perm_idx = ibrpm[j2, j4] - 1
                                    idelt[ntpr - 1, nb2[perm_idx]] = nb1[j4]
                            break
                ntrms = ntpr
                if not has_terms:
                    break

        # Sum color contributions
        clrtot = sum(reslt[:ntrms])
        clfact_val += clrtot

        # Reset for second color term
        ntrms = ntrms0

    return clfact_val


# ===================================================================
# HAMQCD: Evaluate Hamiltonian matrix element between two states
# ===================================================================
def hamqcd(lrt, llt, matrt, matlt, params, selfen):
    """
    Compute Hamiltonian (and free Hamiltonian) matrix elements
    between right state matrt and left state matlt.
    Returns (elem0[NF], elem) where elem0 are free (mass) terms per flavor,
    elem is the interaction term.
    """
    N, NF, B, K = params.N, params.NF, params.B, params.K
    nf = NF

    # Conjugate left state
    matltc = conj_state(llt, matlt)

    elem0 = np.zeros(nf)
    elem = 0.0

    # Build MX: [right state | ham ops | conjugated left state]
    # We'll build it fresh for each Hamiltonian term

    # Collect possible quantum numbers for operators
    nbdrl = [0, 0, 0, 0]
    mqnpos = [[], [], [], []]

    for i9 in range(lrt):
        if matrt[0, i9] == 1:
            mom = matrt[2, i9]
            if mom not in mqnpos[0]:
                mqnpos[0].append(mom)
                nbdrl[0] += 1
        else:
            mom = matrt[2, i9]
            if mom not in mqnpos[1]:
                mqnpos[1].append(mom)
                nbdrl[1] += 1

    for j9 in range(llt):
        if matltc[0, j9] == -1:
            mom = matltc[2, j9]
            if mom not in mqnpos[2]:
                mqnpos[2].append(mom)
                nbdrl[2] += 1
        else:
            mom = matltc[2, j9]
            if mom not in mqnpos[3]:
                mqnpos[3].append(mom)
                nbdrl[3] += 1

    # Helper: count excess operators
    # (mirrors OPMTCH logic)
    nbr = sum(1 for i in range(lrt) if matrt[0, i] == 1)
    ndr = sum(1 for i in range(lrt) if matrt[1, i] == 1)
    nbl = sum(1 for i in range(llt) if matlt[0, i] == 1)
    ndl = sum(1 for i in range(llt) if matlt[1, i] == 1)

    def opmtch(nhb, nhbaj, nhd, nhdaj):
        """Check if Hamiltonian operators can match state operators."""
        if nhb > nbr or nhbaj > nbl or nhd > ndr or nhdaj > ndl:
            return False
        if (nhb + nbl) != (nhbaj + nbr) or (nhd + ndl) != (nhdaj + ndr):
            return False
        # Count unmatched when comparing states
        itempl = matlt[:4, :llt].copy()
        nbsame = 0
        ndsame = 0
        for ir in range(lrt):
            for il in range(llt):
                if (matrt[0, ir] == itempl[0, il] and
                        matrt[2, ir] == itempl[2, il] and
                        matrt[3, ir] == itempl[3, il]):
                    if itempl[0, il] == 1:
                        itempl[0, il] = 2
                        nbsame += 1
                        break
                    else:
                        itempl[1, il] = 2
                        ndsame += 1
                        break
        nbajx = abs(nbr - nbsame)
        nbx = abs(nbl - nbsame)
        ndajx = abs(ndr - ndsame)
        ndx = abs(ndl - ndsame)
        if nhb < nbajx or nhbaj < nbx or nhd < ndajx or nhdaj < ndx:
            return False
        return True

    def build_mx_base(nops_val):
        """Build base MX array with right state and conjugated left state."""
        lng = nops_val + lrt + llt
        mx = np.zeros((4, lng), dtype=int)
        for i in range(lrt):
            for j in range(4):
                mx[j, i] = matrt[j, i]
        for i in range(llt):
            for j in range(4):
                mx[j, lrt + nops_val + i] = matltc[j, i]
        return mx

    def m_index(mu, md):
        """Map operator type to NBDRL index (1-based)."""
        return 2 * mu * mu + 3 * md * md + mu + md - 1  # 0-based

    # === NORM (norh=0): H0 ===
    # (This is handled separately in the calling code)

    # === HAMILTONIAN TERMS ===

    # --- H1: b†b†bb (4-point, quark-quark) ---
    if opmtch(2, 2, 0, 0):
        mx = build_mx_base(4)
        mx[0, lrt:lrt + 4] = [-1, -1, 1, 1]
        mx[1, lrt:lrt + 4] = [0, 0, 0, 0]
        m_ids = [m_index(mx[0, lrt + im], mx[1, lrt + im]) for im in range(4)]

        h1 = 0.0
        for l1 in range(nbdrl[m_ids[0]]):
            for l2 in range(nbdrl[m_ids[1]]):
                for l3 in range(nbdrl[m_ids[2]]):
                    for l4 in range(nbdrl[m_ids[3]]):
                        mx[2, lrt] = mqnpos[m_ids[0]][l1]
                        mx[2, lrt + 1] = mqnpos[m_ids[1]][l2]
                        mx[2, lrt + 2] = mqnpos[m_ids[2]][l3]
                        mx[2, lrt + 3] = mqnpos[m_ids[3]][l4]
                        inth = sum(mx[2, lrt + ll] * (mx[0, lrt + ll] + mx[1, lrt + ll]) for ll in range(4))
                        if inth == 0:
                            k1, k2, k3, k4 = mx[2, lrt], mx[2, lrt + 1], mx[2, lrt + 2], mx[2, lrt + 3]
                            ka, kb = k4 - k2, k3 - k1
                            for if1 in range(1, nf + 1):
                                for jf1 in range(1, nf + 1):
                                    mx[3, lrt:lrt + 4] = [if1, jf1, if1, jf1]
                                    cf = clfact_compute(0, 3, 1, 2, mx, llt, lrt, 4, params, selfen)
                                    h1 += 0.5 * (-cf * brack(ka, kb))
        elem += h1

    # --- H2: d†d†dd (4-point, antiquark-antiquark) ---
    if opmtch(0, 0, 2, 2):
        mx = build_mx_base(4)
        mx[0, lrt:lrt + 4] = [0, 0, 0, 0]
        mx[1, lrt:lrt + 4] = [-1, -1, 1, 1]
        m_ids = [m_index(mx[0, lrt + im], mx[1, lrt + im]) for im in range(4)]

        h2 = 0.0
        for l1 in range(nbdrl[m_ids[0]]):
            for l2 in range(nbdrl[m_ids[1]]):
                for l3 in range(nbdrl[m_ids[2]]):
                    for l4 in range(nbdrl[m_ids[3]]):
                        mx[2, lrt] = mqnpos[m_ids[0]][l1]
                        mx[2, lrt + 1] = mqnpos[m_ids[1]][l2]
                        mx[2, lrt + 2] = mqnpos[m_ids[2]][l3]
                        mx[2, lrt + 3] = mqnpos[m_ids[3]][l4]
                        inth = sum(mx[2, lrt + ll] * (mx[0, lrt + ll] + mx[1, lrt + ll]) for ll in range(4))
                        if inth == 0:
                            k1, k2, k3, k4 = mx[2, lrt], mx[2, lrt + 1], mx[2, lrt + 2], mx[2, lrt + 3]
                            ka, kb = k2 - k4, k1 - k3
                            for if2 in range(1, nf + 1):
                                for jf2 in range(1, nf + 1):
                                    mx[3, lrt:lrt + 4] = [if2, jf2, if2, jf2]
                                    cf = clfact_compute(2, 1, 3, 0, mx, llt, lrt, 4, params, selfen)
                                    h2 += 0.5 * (-cf * brack(ka, kb))
        elem += h2

    # --- H3-H7: quark-antiquark vertices (pair creation/annihilation) ---
    # H3: d†b†b†b
    if opmtch(1, 2, 0, 1):
        mx = build_mx_base(4)
        mx[0, lrt:lrt + 4] = [0, -1, 1, 1]
        mx[1, lrt:lrt + 4] = [1, 0, 0, 0]
        m_ids = [m_index(mx[0, lrt + im], mx[1, lrt + im]) for im in range(4)]
        h3 = 0.0
        for l1 in range(nbdrl[m_ids[0]]):
            for l2 in range(nbdrl[m_ids[1]]):
                for l3 in range(nbdrl[m_ids[2]]):
                    for l4 in range(nbdrl[m_ids[3]]):
                        mx[2, lrt] = mqnpos[m_ids[0]][l1]
                        mx[2, lrt + 1] = mqnpos[m_ids[1]][l2]
                        mx[2, lrt + 2] = mqnpos[m_ids[2]][l3]
                        mx[2, lrt + 3] = mqnpos[m_ids[3]][l4]
                        inth = sum(mx[2, lrt + ll] * (mx[0, lrt + ll] + mx[1, lrt + ll]) for ll in range(4))
                        if inth == 0:
                            k1, k2, k3, k4 = mx[2, lrt], mx[2, lrt + 1], mx[2, lrt + 2], mx[2, lrt + 3]
                            kc, kd = k4 + k1, k3 - k2
                            for if3 in range(1, nf + 1):
                                for jf3 in range(1, nf + 1):
                                    mx[3, lrt:lrt + 4] = [if3, jf3, jf3, if3]
                                    cf = clfact_compute(1, 3, 0, 2, mx, llt, lrt, 4, params, selfen)
                                    h3 += cf * brack(kc, kd)
        elem += h3

    # H4: d†dd†b (pair annihilation + creation)
    if opmtch(0, 1, 1, 2):
        mx = build_mx_base(4)
        mx[0, lrt:lrt + 4] = [0, 0, 0, 1]
        mx[1, lrt:lrt + 4] = [-1, 1, 1, 0]
        m_ids = [m_index(mx[0, lrt + im], mx[1, lrt + im]) for im in range(4)]
        h4 = 0.0
        for l1 in range(nbdrl[m_ids[0]]):
            for l2 in range(nbdrl[m_ids[1]]):
                for l3 in range(nbdrl[m_ids[2]]):
                    for l4 in range(nbdrl[m_ids[3]]):
                        mx[2, lrt] = mqnpos[m_ids[0]][l1]
                        mx[2, lrt + 1] = mqnpos[m_ids[1]][l2]
                        mx[2, lrt + 2] = mqnpos[m_ids[2]][l3]
                        mx[2, lrt + 3] = mqnpos[m_ids[3]][l4]
                        inth = sum(mx[2, lrt + ll] * (mx[0, lrt + ll] + mx[1, lrt + ll]) for ll in range(4))
                        if inth == 0:
                            k1, k2, k3, k4 = mx[2, lrt], mx[2, lrt + 1], mx[2, lrt + 2], mx[2, lrt + 3]
                            ka, kb = k1 - k3, -k4 - k2
                            for if4 in range(1, nf + 1):
                                for jf4 in range(1, nf + 1):
                                    mx[3, lrt:lrt + 4] = [if4, jf4, if4, jf4]
                                    cf = clfact_compute(1, 0, 2, 3, mx, llt, lrt, 4, params, selfen)
                                    h4 += cf * brack(ka, kb)
        elem += h4

    # H5: d†b†bb† (inverted pair)
    if opmtch(2, 1, 1, 0):
        mx = build_mx_base(4)
        mx[0, lrt:lrt + 4] = [0, -1, -1, 1]
        mx[1, lrt:lrt + 4] = [-1, 0, 0, 0]
        m_ids = [m_index(mx[0, lrt + im], mx[1, lrt + im]) for im in range(4)]
        h5 = 0.0
        for l1 in range(nbdrl[m_ids[0]]):
            for l2 in range(nbdrl[m_ids[1]]):
                for l3 in range(nbdrl[m_ids[2]]):
                    for l4 in range(nbdrl[m_ids[3]]):
                        mx[2, lrt] = mqnpos[m_ids[0]][l1]
                        mx[2, lrt + 1] = mqnpos[m_ids[1]][l2]
                        mx[2, lrt + 2] = mqnpos[m_ids[2]][l3]
                        mx[2, lrt + 3] = mqnpos[m_ids[3]][l4]
                        inth = sum(mx[2, lrt + ll] * (mx[0, lrt + ll] + mx[1, lrt + ll]) for ll in range(4))
                        if inth == 0:
                            k1, k2, k3, k4 = mx[2, lrt], mx[2, lrt + 1], mx[2, lrt + 2], mx[2, lrt + 3]
                            kc, kd = k1 + k3, k2 - k4
                            for if5 in range(1, nf + 1):
                                for jf5 in range(1, nf + 1):
                                    mx[3, lrt:lrt + 4] = [if5, jf5, if5, jf5]
                                    cf = clfact_compute(1, 0, 2, 3, mx, llt, lrt, 4, params, selfen)
                                    h5 += cf * brack(kc, kd)
        elem += h5

    # H6: d†d†d b† (quark-antiquark annihilation)
    if opmtch(1, 0, 2, 1):
        mx = build_mx_base(4)
        mx[0, lrt:lrt + 4] = [0, 0, 0, -1]
        mx[1, lrt:lrt + 4] = [-1, -1, 1, 0]
        m_ids = [m_index(mx[0, lrt + im], mx[1, lrt + im]) for im in range(4)]
        h6 = 0.0
        for l1 in range(nbdrl[m_ids[0]]):
            for l2 in range(nbdrl[m_ids[1]]):
                for l3 in range(nbdrl[m_ids[2]]):
                    for l4 in range(nbdrl[m_ids[3]]):
                        mx[2, lrt] = mqnpos[m_ids[0]][l1]
                        mx[2, lrt + 1] = mqnpos[m_ids[1]][l2]
                        mx[2, lrt + 2] = mqnpos[m_ids[2]][l3]
                        mx[2, lrt + 3] = mqnpos[m_ids[3]][l4]
                        inth = sum(mx[2, lrt + ll] * (mx[0, lrt + ll] + mx[1, lrt + ll]) for ll in range(4))
                        if inth == 0:
                            k1, k2, k3, k4 = mx[2, lrt], mx[2, lrt + 1], mx[2, lrt + 2], mx[2, lrt + 3]
                            ka, kb = k2 - k3, k1 + k4
                            for if6 in range(1, nf + 1):
                                for jf6 in range(1, nf + 1):
                                    mx[3, lrt:lrt + 4] = [if6, jf6, jf6, if6]
                                    cf = clfact_compute(3, 1, 2, 0, mx, llt, lrt, 4, params, selfen)
                                    h6 += cf * brack(ka, kb)
        elem += h6

    # H7: d†d b†b (quark-antiquark exchange)
    if opmtch(1, 1, 1, 1):
        mx = build_mx_base(4)
        mx[0, lrt:lrt + 4] = [0, 0, -1, 1]
        mx[1, lrt:lrt + 4] = [-1, 1, 0, 0]
        m_ids = [m_index(mx[0, lrt + im], mx[1, lrt + im]) for im in range(4)]
        h7 = 0.0
        for l1 in range(nbdrl[m_ids[0]]):
            for l2 in range(nbdrl[m_ids[1]]):
                for l3 in range(nbdrl[m_ids[2]]):
                    for l4 in range(nbdrl[m_ids[3]]):
                        mx[2, lrt] = mqnpos[m_ids[0]][l1]
                        mx[2, lrt + 1] = mqnpos[m_ids[1]][l2]
                        mx[2, lrt + 2] = mqnpos[m_ids[2]][l3]
                        mx[2, lrt + 3] = mqnpos[m_ids[3]][l4]
                        inth = sum(mx[2, lrt + ll] * (mx[0, lrt + ll] + mx[1, lrt + ll]) for ll in range(4))
                        if inth == 0:
                            k1, k2, k3, k4 = mx[2, lrt], mx[2, lrt + 1], mx[2, lrt + 2], mx[2, lrt + 3]
                            ka, kb = k4 - k3, k2 - k1
                            kc, kd = k4 + k2, -k1 - k3
                            # Two flavor structures for H7
                            for if7a in range(1, nf + 1):
                                for jf7a in range(1, nf + 1):
                                    mx[3, lrt:lrt + 4] = [if7a, if7a, jf7a, jf7a]
                                    cf = clfact_compute(1, 3, 2, 0, mx, llt, lrt, 4, params, selfen)
                                    h7 += 0.5 * 2.0 * (-cf * brack(ka, kb))
                            for if7b in range(1, nf + 1):
                                for jf7b in range(1, nf + 1):
                                    mx[3, lrt:lrt + 4] = [if7b, jf7b, if7b, jf7b]
                                    cf = clfact_compute(2, 3, 1, 0, mx, llt, lrt, 4, params, selfen)
                                    h7 += 0.5 * 2.0 * cf * brack(kc, kd)
        elem += h7

    # --- Diagonal (2-point) vertices ---
    cbreak = params.cbreak

    # H8: b†b (quark self-energy)
    if opmtch(1, 1, 0, 0):
        mx = build_mx_base(2)
        mx[0, lrt:lrt + 2] = [-1, 1]
        mx[1, lrt:lrt + 2] = [0, 0]
        m_ids = [m_index(-1, 0), m_index(1, 0)]

        h80 = np.zeros(nf)
        h8 = 0.0
        for l1 in range(nbdrl[m_ids[0]]):
            for l2 in range(nbdrl[m_ids[1]]):
                mx[2, lrt] = mqnpos[m_ids[0]][l1]
                mx[2, lrt + 1] = mqnpos[m_ids[1]][l2]
                inth = sum(mx[2, lrt + ll] * (mx[0, lrt + ll] + mx[1, lrt + ll]) for ll in range(2))
                if inth == 0:
                    k1 = mx[2, lrt]
                    k2 = mx[2, lrt + 1]
                    for ifl8 in range(1, nf + 1):
                        mx[3, lrt:lrt + 2] = [ifl8, ifl8]
                        if k1 != 0:
                            cf = clfact_compute(0, 1, 0, 0, mx, llt, lrt, 2, params, selfen)
                            delta = float(jdelta(k1, k2))
                            h80[ifl8 - 1] += cf * delta * (2.0 + cbreak) / k1
                            h8 += cf * delta * 0.5 * selfen[k1 - 1] if k1 <= len(selfen) else 0.0
        elem0 += h80
        elem += h8

    # H9: d†d (antiquark self-energy)
    if opmtch(0, 0, 1, 1):
        mx = build_mx_base(2)
        mx[0, lrt:lrt + 2] = [0, 0]
        mx[1, lrt:lrt + 2] = [-1, 1]
        m_ids = [m_index(0, -1), m_index(0, 1)]

        h90 = np.zeros(nf)
        h9 = 0.0
        for l1 in range(nbdrl[m_ids[0]]):
            for l2 in range(nbdrl[m_ids[1]]):
                mx[2, lrt] = mqnpos[m_ids[0]][l1]
                mx[2, lrt + 1] = mqnpos[m_ids[1]][l2]
                inth = sum(mx[2, lrt + ll] * (mx[0, lrt + ll] + mx[1, lrt + ll]) for ll in range(2))
                if inth == 0:
                    k1 = mx[2, lrt]
                    k2 = mx[2, lrt + 1]
                    for ifl9 in range(1, nf + 1):
                        mx[3, lrt:lrt + 2] = [ifl9, ifl9]
                        if k1 != 0:
                            cf = clfact_compute(1, 0, 0, 0, mx, llt, lrt, 2, params, selfen)
                            delta = float(jdelta(k1, k2))
                            h90[ifl9 - 1] += cf * delta * (2.0 - cbreak) / k1
                            h9 += cf * delta * 0.5 * selfen[k1 - 1] if k1 <= len(selfen) else 0.0
        elem0 += h90
        elem += h9

    return elem0, elem


# ===================================================================
# CLRDIS: Build norm or Hamiltonian matrix (PARALLELIZED)
# ===================================================================
def _compute_element(args):
    """Worker function for parallel Hamiltonian construction."""
    (nstrt, nstlt, mstate_flat, mstinf_flat, numsta,
     params_dict, selfen, norh, N, NF, B, K) = args

    # Reconstruct arrays from flat data
    mstate = mstate_flat.reshape(NMXFR, MXNP)
    mstinf = mstinf_flat.reshape(NMSTMX, 8)

    # Reconstruct params
    params = Params()
    params.N = params_dict['N']
    params.NF = params_dict['NF']
    params.B = params_dict['B']
    params.K = params_dict['K']
    params.MASS = params_dict['MASS']
    params.rmq = np.array(params_dict['rmq'])
    params.iflv = np.array(params_dict['iflv'])
    params.rlamb = params_dict['rlamb']
    params.cutoff = params_dict['cutoff']
    params.LPN = params_dict['LPN']
    params.cbreak = params_dict['cbreak']

    # Extract states
    locstr = mstinf[nstrt, 0] - 1  # 0-based
    lstrt = mstinf[nstrt, 1]
    locstl = mstinf[nstlt, 0] - 1
    lstlt = mstinf[nstlt, 1]

    matrt = np.zeros((4, MXNP), dtype=int)
    matlt = np.zeros((4, MXNP), dtype=int)
    for j1 in range(4):
        for j2 in range(lstrt):
            matrt[j1, j2] = mstate[locstr + j1, j2]
    for k1 in range(4):
        for k2 in range(lstlt):
            matlt[k1, k2] = mstate[locstl + k1, k2]

    # Quick check: can any Hamiltonian term connect these states?
    itempl = matlt[:4, :lstlt].copy()
    nbsame = ndsame = 0
    for ir in range(lstrt):
        for il in range(lstlt):
            if (matrt[0, ir] == itempl[0, il] and
                    matrt[2, ir] == itempl[2, il] and
                    matrt[3, ir] == itempl[3, il]):
                if itempl[0, il] == 1:
                    itempl[0, il] = 2
                    nbsame += 1
                    break
                else:
                    itempl[1, il] = 2
                    ndsame += 1
                    break

    nbr_loc = sum(1 for i in range(lstrt) if matrt[0, i] == 1)
    ndr_loc = sum(1 for i in range(lstrt) if matrt[1, i] == 1)
    nbl_loc = sum(1 for i in range(lstlt) if matlt[0, i] == 1)
    ndl_loc = sum(1 for i in range(lstlt) if matlt[1, i] == 1)

    nbajx = abs(nbr_loc - nbsame)
    nbx = abs(nbl_loc - nbsame)
    ndajx = abs(ndr_loc - ndsame)
    ndx = abs(ndl_loc - ndsame)
    nbdx = nbajx + nbx + ndajx + ndx

    if nbdx > 4:  # max ops in any Hamiltonian term
        return nstrt, nstlt, np.zeros(NF), 0.0, 0.0

    if norh == 0:
        # Norm matrix: just overlap (H0 with no operators)
        if nbdx == 0:
            mx = np.zeros((4, lstrt + lstlt), dtype=int)
            matltc = conj_state(lstlt, matlt)
            for i in range(lstrt):
                for j in range(4):
                    mx[j, i] = matrt[j, i]
            for i in range(lstlt):
                for j in range(4):
                    mx[j, lstrt + i] = matltc[j, i]
            h0 = clfact_compute(0, 0, 0, 0, mx, lstlt, lstrt, 0, params, selfen)
            return nstrt, nstlt, np.zeros(NF), 0.0, h0
        else:
            return nstrt, nstlt, np.zeros(NF), 0.0, 0.0
    else:
        # Hamiltonian
        elem0, elem = hamqcd(lstrt, lstlt, matrt, matlt, params, selfen)
        return nstrt, nstlt, elem0, elem, 0.0


def clrdis(norh, params, states, selfen, ncpus=None):
    """
    Build norm matrix (norh=0) or Hamiltonian (norh=1).
    PARALLELIZED using multiprocessing.
    Returns (ham0, ham, hnorm) as appropriate.
    """
    numsta = states.numsta
    N, NF = params.N, params.NF

    ham0 = np.zeros((NF, numsta, numsta))
    ham = np.zeros((numsta, numsta))
    hnorm = np.zeros((numsta, numsta))

    # Flatten data for multiprocessing
    params_dict = {
        'N': params.N, 'NF': params.NF, 'B': params.B, 'K': params.K,
        'MASS': params.MASS, 'rmq': params.rmq.tolist(),
        'iflv': params.iflv.tolist(), 'rlamb': params.rlamb,
        'cutoff': params.cutoff, 'LPN': params.LPN, 'cbreak': params.cbreak
    }
    mstate_flat = states.mstate.flatten()
    mstinf_flat = states.mstinf.flatten()

    # Build task list (upper triangle only)
    tasks = []
    for nstrt in range(numsta):
        for nstlt in range(nstrt, numsta):
            tasks.append((nstrt, nstlt, mstate_flat, mstinf_flat, numsta,
                          params_dict, selfen, norh, params.N, params.NF,
                          params.B, params.K))

    log.info(f"  CLRDIS({'norm' if norh == 0 else 'ham'}): {len(tasks)} matrix elements, "
             f"{numsta} states")

    if ncpus is None:
        ncpus = int(os.environ.get('QCDF_NCPUS', min(cpu_count(), 8)))

    t0 = time.time()

    if ncpus > 1 and len(tasks) > 100:
        log.info(f"  Using {ncpus} parallel workers")
        with Pool(ncpus) as pool:
            results = pool.map(_compute_element, tasks, chunksize=max(1, len(tasks) // (ncpus * 4)))
    else:
        log.info("  Running serially")
        results = [_compute_element(t) for t in tasks]

    # Assemble matrices
    for (nstrt, nstlt, elem0, elem, norm_val) in results:
        if norh == 0:
            hnorm[nstlt, nstrt] = hnorm[nstlt, nstrt] + norm_val
            hnorm[nstrt, nstlt] = hnorm[nstlt, nstrt]
        else:
            for ifl in range(NF):
                ham0[ifl, nstlt, nstrt] += elem0[ifl]
                ham0[ifl, nstrt, nstlt] = ham0[ifl, nstlt, nstrt]
            ham[nstlt, nstrt] += elem
            ham[nstrt, nstlt] = ham[nstlt, nstrt]

    dt = time.time() - t0
    log.info(f"  CLRDIS complete in {dt:.1f}s")

    return ham0, ham, hnorm


# ===================================================================
# WEEDR: Weed redundant states from norm matrix
# ===================================================================
def weedr(states, hnorm):
    """
    Search norm matrix to weed out redundant states.
    Modifies states.numsta, states.mstinf, and hnorm in place.
    """
    eps = 1.0e-4
    nstpr = states.numsta

    loc = 0
    while loc < nstpr - 1:
        ndump = 0
        mdrop = []
        j = loc + 1
        while j < nstpr:
            if abs(hnorm[j, loc]) > eps:
                idump = True
                r = hnorm[j, loc] / hnorm[loc, loc]
                for ir in range(nstpr):
                    rh = r * hnorm[loc, ir]
                    diff = rh - hnorm[j, ir]
                    if abs(diff) > eps:
                        idump = False
                        break
                if idump:
                    ndump += 1
                    mdrop.append(j)
            j += 1

        if ndump > 0:
            _dropr(states, hnorm, nstpr, mdrop)
            nstpr -= ndump
        loc += 1

    states.numsta = nstpr


def weedr2(states, hnorm, w, z):
    """Weed states with zero norm eigenvalues."""
    eps = 1.0e-4
    numsta = states.numsta
    nzer = sum(1 for i in range(numsta) if abs(w[i]) < eps)

    if nzer > 0:
        idr = 0
        mdrop = []
        indchk = np.zeros(numsta, dtype=int)
        for i1 in range(nzer):
            idis = False
            for i2 in range(numsta - 1, -1, -1):
                if abs(z[i2, i1]) > eps and indchk[i2] == 0:
                    if not idis:
                        idr += 1
                        mdrop.append(i2)
                        idis = True
                    indchk[i2] = 1

        if idr > 0:
            log.info(f"  Dropping {idr} states after norm diagonalization")
            _dropr(states, hnorm, states.numsta, mdrop)
            return True
    return False


def _dropr(states, hnorm, nstpr, mdrop):
    """Drop states from hnorm and mstinf."""
    for idx in sorted(mdrop, reverse=True):
        if idx == nstpr - 1:
            hnorm[nstpr - 1, :nstpr] = 0
            hnorm[:nstpr, nstpr - 1] = 0
            states.mstinf[nstpr - 1, :] = 0
        else:
            # Shift rows/cols
            for jdr in range(idx, nstpr - 1):
                hnorm[jdr, :nstpr] = hnorm[jdr + 1, :nstpr]
                hnorm[jdr + 1, :nstpr] = 0
            for ldr in range(idx, nstpr - 1):
                hnorm[:nstpr, ldr] = hnorm[:nstpr, ldr + 1]
                hnorm[:nstpr, ldr + 1] = 0
            for ns in range(idx, nstpr - 1):
                states.mstinf[ns, :] = states.mstinf[ns + 1, :]
                states.mstinf[ns + 1, :] = 0
        nstpr -= 1
    states.numsta = nstpr


# ===================================================================
# NUZ: Normalize Z by eigenvalue square roots
# ===================================================================
def nuz(w, z, numsta):
    """Divide eigenvectors by sqrt of eigenvalues to normalize states."""
    for mr in range(numsta):
        z[:numsta, mr] /= np.sqrt(w[mr])
    return z


# ===================================================================
# NUHAM: Compute Hamiltonian in orthonormal basis
# ===================================================================
def nuham(ham0, ham, z, params, numsta):
    """
    Compute new Hamiltonian in orthonormal basis.
    Returns (hnu0[NF, numsta], hnu[numsta, numsta]).
    """
    NF = params.NF
    eps = 1.0e-10

    # HNU0: diagonal free Hamiltonian per flavor
    hnu0 = np.zeros((NF, numsta))
    for ifl in range(NF):
        for nl in range(numsta):
            hlr0 = 0.0
            for kl in range(numsta):
                zl = z[kl, nl]
                if abs(zl) > eps:
                    for kr in range(numsta):
                        zr = z[kr, nl]
                        if abs(zr) > eps:
                            hlr0 += ham0[ifl, kl, kr] * zl * zr
            hnu0[ifl, nl] = hlr0

    # HNU: interaction Hamiltonian
    hnu = np.zeros((numsta, numsta))
    for nl in range(numsta):
        for nr in range(nl + 1):
            hlr = 0.0
            for kl in range(numsta):
                zl = z[kl, nl]
                if abs(zl) > eps:
                    for kr in range(numsta):
                        zr = z[kr, nr]
                        if abs(zr) > eps:
                            hlr += ham[kl, kr] * zl * zr
            hnu[nl, nr] = hlr
            hnu[nr, nl] = hlr

    return hnu0, hnu


# ===================================================================
# MAIN PROGRAM
# ===================================================================
def main(input_file=None):
    """Main DLCQ QCD program."""
    log.info("=" * 60)
    log.info(" DLCQ QCD with Multiple Flavors (Python/Parallel)")
    log.info(" Original Fortran by K. Hornbostel (1993)")
    log.info("=" * 60)

    params = Params()
    states = StateData()
    perm = PermTables()
    flavtab = FlavorTables()

    # --- Read parameters ---
    if input_file and os.path.exists(input_file):
        with open(input_file) as f:
            cfg = json.load(f)
        ihfile = cfg.get('ihfile', 0)
        params.N = cfg['N']
        params.NF = cfg['NF']
        params.B = cfg['B']
        params.rlamb = cfg['rlamb']
        params.K = cfg['K']
        params.cutoff = cfg.get('cutoff', -1.0)
        params.LPN = cfg.get('LPN', 0)
        if 'rmq' in cfg:
            for i, v in enumerate(cfg['rmq']):
                params.rmq[i] = v
        if 'iflv' in cfg:
            for i, v in enumerate(cfg['iflv']):
                params.iflv[i] = v
    else:
        ihfile = int(input("DOES A HAMILTONIAN FILE EXIST? YES: 1; NO: 0\n> "))
        if ihfile != 1:
            params.N = int(input("INPUT N: "))
            params.NF = int(input("INPUT NF: "))
            params.B = int(input("INPUT B: "))
            params.rlamb = float(input("INPUT LAMBDA (0.0<LAMBDA<1.0): "))

            params.rmq[0] = 1.0
            if params.NF > 1:
                for nqk in range(1, params.NF):
                    params.rmq[nqk] = float(input(
                        f"Input ratio of quark mass {nqk + 1} to quark mass 1: "))
                for nfl in range(params.NF):
                    params.iflv[nfl] = int(input(
                        f"Input flavor qm number for flavor {nfl + 1}: "))

            params.cutoff = float(input("CUTOFF/MASS SQD (LE.0 = NO CUTOFF): "))
            params.LPN = int(input("LIMIT PARTICLE NUMBER TO (0 = NO LIMIT): "))

            nqmin = params.N * params.B
            kmin = minmom(nqmin, params.N, params.NF) if params.B != 0 else 2
            params.K = int(input(f"INPUT K (GE {kmin}): "))
        else:
            params.rlamb = float(input("INPUT LAMBDA (0.0<LAMBDA<1.0): "))

    log.info(f"\n  N={params.N}, NF={params.NF}, B={params.B}, K={params.K}, "
             f"lambda={params.rlamb}")
    log.info(f"  LPN={params.LPN}, cutoff={params.cutoff}")
    log.info(f"  rmq={params.rmq[:params.NF]}")
    log.info(f"  iflv={params.iflv[:params.NF]}")

    # --- Output file ---
    outf = open('qcdf_py.out', 'w')
    outf.write(f" N,NF,B,K,LAMBDA:  {params.N:5d}{params.NF:5d}"
               f"{params.B:5d}{params.K:5d} {params.rlamb:24.16e}\n")
    outf.write(f" LPN:  {params.LPN:5d}\n")
    outf.write(f" CUTOFF/MASS SQD:  {params.cutoff:24.16e}\n")

    if ihfile != 1:
        # --- COMPUTE STATES AND HAMILTONIAN ---
        selfen = compute_selfen(params)
        log.info("\n  Computing self-energies...")

        log.info("  Generating states...")
        qcdsta(params, states, perm, flavtab)
        outf.write(f" NUMBER OF STATES BEFORE WEEDING: {states.numsta:6d}\n")
        log.info(f"  States before weeding: {states.numsta}")

        if states.numsta == 0:
            log.info("  No states generated. Stopping.")
            outf.close()
            return

        # --- Compute norm matrix ---
        log.info("\n  Computing norm matrix...")
        ham0, ham, hnorm = clrdis(0, params, states, selfen)

        # --- Weed redundant states ---
        log.info("  Weeding redundant states...")
        weedr(states, hnorm)
        log.info(f"  States after linear weeding: {states.numsta}")

        # --- Diagonalize norm and weed zero eigenvalues ---
        log.info("  Diagonalizing norm matrix...")
        numsta = states.numsta
        for iw in range(2000):
            hnorm_sub = hnorm[:numsta, :numsta].copy()
            w, z_mat = eigh(hnorm_sub)
            # Sort by eigenvalue (ascending)
            idx = np.argsort(w)
            w = w[idx]
            z_mat = z_mat[:, idx]

            changed = weedr2(states, hnorm, w, z_mat)
            numsta = states.numsta
            if not changed:
                break

        outf.write(f" NUMBER OF STATES AFTER WEEDING: {states.numsta:6d}\n")
        log.info(f"  States after norm weeding: {states.numsta}")

        if states.numsta == 0:
            log.info("  No states remain. Stopping.")
            outf.close()
            return

        # Final norm diag for Z
        numsta = states.numsta
        hnorm_sub = hnorm[:numsta, :numsta].copy()
        w_norm, z_norm = eigh(hnorm_sub)
        idx = np.argsort(w_norm)
        w_norm = w_norm[idx]
        z_norm = z_norm[:, idx]

        # Normalize Z.  Sized to numsta, not NMSTMX: nuz and nuham only ever
        # touch [:numsta], so allocating the cap wasted 381 MB at the old
        # NMSTMX and would be 5 GB at the raised one.
        z_full = np.zeros((numsta, numsta))
        z_full[:numsta, :numsta] = z_norm
        z_full = nuz(w_norm, z_full, numsta)

        # --- Compute Hamiltonian ---
        log.info("\n  Computing Hamiltonian matrix...")
        ham0, ham, _ = clrdis(1, params, states, selfen)

        # Transform to orthonormal basis
        log.info("  Transforming to orthonormal basis...")
        hnu0, hnu = nuham(ham0, ham, z_full, params, numsta)

        # --- Save Hamiltonian ---
        log.info("  Saving Hamiltonian to qcdf_py.ham...")
        with open('qcdf_py.ham', 'w') as hf:
            hf.write(f"     N     NF    B    LPN    K   NUMSTA   CUTOFF\n")
            hf.write(f"{params.N:6d}{params.NF:6d}{params.B:6d}"
                     f"{params.LPN:6d}{params.K:6d}{numsta:6d}      "
                     f"{params.cutoff:24.16e}\n")
            hf.write(f"      HNU0(nf):   diagonal terms for each flavor\n")
            for ir1 in range(numsta):
                vals = " ".join(f"{hnu0[ifl, ir1]:24.16e}" for ifl in range(params.NF))
                hf.write(f"{vals}\n")
            hf.write(f"      HNU:   non-zero elements of interaction hamilt.\n")
            eps_h = 1.0e-15
            nbuff = 0
            buf = []
            for ir1 in range(numsta):
                for ir2 in range(ir1 + 1):
                    if abs(hnu[ir1, ir2]) > eps_h:
                        buf.append((ir1 + 1, ir2 + 1, hnu[ir1, ir2]))
                        if len(buf) == 3:
                            line = "".join(f"{b[0]:6d}{b[1]:6d} {b[2]:24.16e} " for b in buf)
                            hf.write(line + "\n")
                            buf = []
            # Write remaining
            while len(buf) < 3:
                buf.append((-1, -1, 0.0))
            line = "".join(f"{b[0]:6d}{b[1]:6d} {b[2]:24.16e} " for b in buf)
            hf.write(line + "\n")
    else:
        # --- READ HAMILTONIAN FROM FILE ---
        log.info("  Reading Hamiltonian from qcdf_py.ham...")
        with open('qcdf_py.ham', 'r') as hf:
            hf.readline()  # header
            parts = hf.readline().split()
            params.N = int(parts[0])
            params.NF = int(parts[1])
            params.B = int(parts[2])
            params.LPN = int(parts[3])
            params.K = int(parts[4])
            numsta = int(parts[5])
            params.cutoff = float(parts[6])
            states.numsta = numsta

            hf.readline()  # HNU0 header
            hnu0 = np.zeros((params.NF, numsta))
            for ir1 in range(numsta):
                vals = list(map(float, hf.readline().split()))
                for ifl in range(params.NF):
                    hnu0[ifl, ir1] = vals[ifl]

            hf.readline()  # HNU header
            hnu = np.zeros((numsta, numsta))
            while True:
                line = hf.readline()
                if not line:
                    break
                parts = line.split()
                for k in range(0, len(parts), 3):
                    if k + 2 < len(parts):
                        i = int(parts[k])
                        j = int(parts[k + 1])
                        v = float(parts[k + 2])
                        if i == -1:
                            break
                        hnu[i - 1, j - 1] = v
                        hnu[j - 1, i - 1] = v

        # Input quark masses
        params.rmq[0] = 1.0
        if params.NF > 1:
            for nqk in range(1, params.NF):
                params.rmq[nqk] = float(input(
                    f"Input ratio of quark mass {nqk + 1} to quark mass 1: "))

    # --- Combine free and interacting Hamiltonians ---
    numsta = states.numsta
    log.info(f"\n  Combining free + interacting Hamiltonians (numsta={numsta})...")

    rlmsq = params.rlamb ** 2

    # Scale interaction by coupling
    for i1 in range(numsta):
        for i2 in range(i1 + 1):
            hnu[i1, i2] *= rlmsq

    # Add diagonal (free) terms
    for i1 in range(numsta):
        for ifl in range(params.NF):
            hnu[i1, i1] += (1.0 - rlmsq) * params.rmq[ifl] ** 2 * hnu0[ifl, i1]

    # Symmetrize
    for i1 in range(numsta):
        for i2 in range(i1):
            hnu[i2, i1] = hnu[i1, i2]

    # --- Diagonalize ---
    log.info("  Diagonalizing Hamiltonian...")
    hnu_sub = hnu[:numsta, :numsta].copy()
    w_eig, z_eig = eigh(hnu_sub)

    # Convert to mass² eigenvalues: M² = K * eigenvalue / 2
    w_mass = params.K * w_eig / 2.0

    # Sort by mass²
    idx = np.argsort(w_mass)
    w_mass = w_mass[idx]
    z_eig = z_eig[:, idx]

    # --- Output ---
    outf.write(f"\n IERR =     0\n\n")
    outf.write(f" COMPUTED EIGENVALUES (MASS SQUARED):\n\n")
    for i in range(0, numsta, 3):
        vals = w_mass[i:min(i + 3, numsta)]
        line = "".join(f"{v:24.16e}" for v in vals)
        outf.write(line + "\n")

    outf.write(f"\n\n COMPUTED EIGENVECTORS:\n\n")
    n_print = min(numsta, IVCTMX)
    for l in range(0, n_print, 4):
        for i in range(numsta):
            nn = min(l + 4, n_print)
            vals = z_eig[i, l:nn]
            line = "".join(f"{v:24.16e}" for v in vals)
            outf.write(line + "\n")
        outf.write(" ------\n")

    outf.close()

    # --- Print summary ---
    log.info("\n" + "=" * 60)
    log.info(" MASS² EIGENVALUES (lowest 20):")
    log.info("=" * 60)
    for i in range(min(20, numsta)):
        log.info(f"  {i + 1:4d}:  {w_mass[i]:20.12e}")
    log.info("=" * 60)
    log.info(f" Output written to qcdf_py.out")
    log.info(f" Hamiltonian saved to qcdf_py.ham")


# ===================================================================
if __name__ == '__main__':
    input_file = sys.argv[1] if len(sys.argv) > 1 else None
    main(input_file)
