#!/usr/bin/env python3
"""
qcdf_opt.py — DLCQ QCD 1+1D, optimized Python.

Key optimizations vs qcdf.py:
  1. Pre-allocated work buffers (eliminates ~8000 numpy.zeros/call)
  2. Vectorized basis transform: Z.T @ H @ Z  (replaces O(N⁴) loops)
  3. Row-batched parallelism (1 IPC call/row instead of 1/element)
  4. Dynamic matrix sizing (allocate to actual numsta, not 6902)
  5. Numba JIT for hot micro-kernels (brack, stachk, conj_state, selfen)
  6. Table-driven vertex loop in hamqcd (less code, same logic)
"""
import sys, os, json, time, logging, threading
import numpy as np
from scipy.linalg import eigh
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Pool, cpu_count
from numba import njit

# Import everything except the hot-path functions from the original
from qcdf import (
    Params, StateData, PermTables, FlavorTables,
    lpnsub, qcdsta, NMSTMX, NMXFR,
    prmx, prmm, ibarfl, imesfl, minmom, iodev,
)
MXTRM = 200000; MXLNG = 54; MXP = 25
# MXLNG = 2*MXP + 4 (max right state + max operators + max left state)

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger(__name__)

# ─── Numba micro-kernels ────────────────────────────────────────
@njit(cache=True)
def brack(L, M):
    if L == 0 or M == 0 or L + M != 0: return 0.0
    return 4.0 / (float(L) * float(L))

@njit(cache=True)
def conj_state_nb(lst, ist, out):
    for i in range(lst):
        j = lst - 1 - i
        out[0,j]=-ist[0,i]; out[1,j]=-ist[1,i]; out[2,j]=ist[2,i]; out[3,j]=ist[3,i]

@njit(cache=True)
def compute_selfen(N):
    s = np.zeros(100)
    for i in range(0, 97, 2): s[i+2] = s[i] + 4.0/float((i+2)**2)
    cf = (N*N - 1.0)/N
    for i in range(0, 99, 2): s[i] *= cf
    return s

# ─── Color factor with PRE-ALLOCATED buffers ────────────────────
# These are allocated once, reused on every clfact call (NOT thread-safe,
# but each process in multiprocessing gets its own copy via fork)
class ClfactBuffers:
    def __init__(self):
        self.idel0 = np.zeros((MXTRM, MXLNG), dtype=np.int32)
        self.resl0 = np.zeros(MXTRM)
        self.idelt = np.zeros((MXTRM, MXLNG), dtype=np.int32)
        self.reslt = np.zeros(MXTRM)
        self.my    = np.zeros((4, MXLNG), dtype=np.int32)
        self.lnkb  = np.zeros(MXLNG, dtype=np.int32)
        self.icycl = np.ones(MXLNG, dtype=np.int32)
        self.ncr   = np.zeros((2, MXLNG), dtype=np.int32)
        self.nperm = np.zeros((MXLNG, MXLNG), dtype=np.int32)

_buf = ClfactBuffers()

def clfact(ic1, ic2, ic3, ic4, mx, lng, llt, lrt, nops, N, NF, B, K):
    """Color factor using pre-allocated work buffers."""
    b = _buf
    idel0=b.idel0; resl0=b.resl0; my=b.my
    lnkb=b.lnkb; icycl=b.icycl; ncr=b.ncr; nperm=b.nperm

    # Reset only what we need
    idel0[0, :lng] = 0
    my[:, :lng] = mx[:, :lng]
    resl0[0] = 1.0

    # ── MOCHK: verify operator balance at each momentum+flavor ──
    # Scan L→R and R→L simultaneously; if B or D count at any
    # momentum exceeds N or goes negative, color factor is 0.
    maxmom = K + 2
    for k3 in range(2):
        fwd = np.zeros((maxmom, NF + 1), dtype=np.int32)
        rev = np.zeros((maxmom, NF + 1), dtype=np.int32)
        for i3 in range(lng):
            ri = lng - 1 - i3
            mom_f = mx[2, i3]; flv_f = mx[3, i3]
            mom_r = mx[2, ri]; flv_r = mx[3, ri]
            fwd[mom_f, flv_f] += mx[k3, i3]
            rev[mom_r, flv_r] += mx[k3, ri]
            if (fwd[mom_f, flv_f] < 0 or fwd[mom_f, flv_f] > N or
                    rev[mom_r, flv_r] > 0 or rev[mom_r, flv_r] < -N):
                return 0.0

    ipow = 0
    for ja in range(1, lng):
        if my[0,ja]+my[1,ja] != -1: continue
        for ka in range(ja-1, -1, -1):
            if (my[0,ka]==-my[0,ja] and my[1,ka]==-my[1,ja]
                    and my[2,ka]==my[2,ja] and my[3,ka]==my[3,ja]):
                idel0[0,ja]=ka; idel0[0,ka]=ja
                my[0,ka]=my[1,ka]=my[0,ja]=my[1,ja]=0
                break
            elif my[0,ka]!=0 or my[1,ka]!=0:
                ipow += 1
    resl0[0] = (-1.0)**ipow

    lnkb[:lng]=0; icycl[:lng]=1
    if B:
        for nbs in range(B):
            for i in range(N*nbs, N*(nbs+1)-1): lnkb[i]=i+1
            lnkb[N*(nbs+1)-1]=N*nbs
            for j in range(lng-N*nbs-1, lng-N*(nbs+1), -1): lnkb[j]=j-1
            lnkb[lng-N*(nbs+1)]=lng-N*nbs-1
        if (-1)**N > 0:
            for nbs in range(B):
                s=1
                for i in range(N*nbs, N*(nbs+1)): icycl[i]=s; s*=-1
                s=1
                for j in range(lng-N*nbs-1, lng-N*(nbs+1)-1, -1): icycl[j]=s; s*=-1

    ncr[:,:lng]=0; nperm[:lng,:lng]=0
    lbcr=ldcr=lprm=0; lpo=1; ntrms=1

    for i in range(lng):
        for bd in range(2):
            is_b=(bd==0 and mx[1,i]==0)
            is_d=(bd==1 and mx[0,i]==0 and mx[1,i]!=0)
            if not(is_b or is_d): continue
            cr=bd; lcr=lbcr if bd==0 else ldcr
            if mx[bd,i]==-1:
                lprm+=1; lpchk=False; jfst=-1
                for j in range(lcr-1,-1,-1):
                    if not lpchk and mx[2,ncr[cr,j]]==mx[2,i] and mx[3,ncr[cr,j]]==mx[3,i]:
                        jfst=j
                        for ll in range(j-1,-1,-1):
                            if mx[2,ncr[cr,ll]]==mx[2,ncr[cr,j]] and mx[3,ncr[cr,ll]]==mx[3,ncr[cr,j]]:
                                nperm[lprm-1,0]+=1; mp=2*nperm[lprm-1,0]
                                nperm[lprm-1,mp-1]=ncr[cr,j]; nperm[lprm-1,mp]=ncr[cr,ll]
                                lpchk=True
                        break
                if not lpchk: lprm-=1
                if jfst>=0:
                    for m in range(jfst,lcr-1): ncr[cr,m]=ncr[cr,m+1]
                    ncr[cr,lcr-1]=0
                    if bd==0: lbcr-=1
                    else: ldcr-=1
            elif mx[bd,i]==1:
                ncr[cr,lcr]=i
                if bd==0: lbcr+=1
                else: ldcr+=1

    if lprm>0:
        for pid in range(lprm-1,-1,-1):
            nsbar=1
            for jd in range(nperm[pid,0]):
                kd=2*(jd+1); ip1=nperm[pid,kd-1]; ip2=nperm[pid,kd]
                ismb=False
                if B>0 and lnkb[ip1]!=0:
                    nxt=lnkb[ip1]
                    for _ in range(N-1):
                        if ip2==nxt: ismb=True; break
                        nxt=lnkb[nxt]
                if ismb: nsbar+=1
                else:
                    for j in range(ntrms):
                        # Must stay in lockstep with the same guard in
                        # qcdf_kernels.clfact_nb: the two paths sharing MXTRM is
                        # exactly why they agreed bit-for-bit on a wrong answer.
                        if lpo>=MXTRM:
                            raise OverflowError(
                                "colour contraction exceeded MXTRM terms; the "
                                "result would be silently wrong. Raise MXTRM in "
                                "BOTH qcdf_opt.py and qcdf_kernels.py (they must "
                                "match), or reduce the particle-number cutoff LPN.")
                        idel0[lpo,:lng]=idel0[j,:lng]; resl0[lpo]=-resl0[j]
                        idel0[lpo,ip1]=idel0[j,ip2]; idel0[lpo,ip2]=idel0[j,ip1]
                        v1,v2=idel0[j,ip2],idel0[j,ip1]
                        if 0<=v1<lng: idel0[lpo,v1]=idel0[j,v2]
                        if 0<=v2<lng: idel0[lpo,v2]=idel0[j,v1]
                        lpo+=1
            if nsbar>1:
                for im in range(ntrms): resl0[im]*=nsbar
            ntrms*=(nperm[pid,0]+2-nsbar)

    ntrms0=ntrms; nmsr=(lrt-N*B)//2; nmsl=(llt-N*B)//2
    nctf=2 if nops==4 else 1; total=0.0
    idelt=b.idelt; reslt=b.reslt

    for nct in range(1,nctf+1):
        idelt[:ntrms,:lng]=0; reslt[:ntrms]=resl0[:ntrms]
        for j in range(lng):
            if mx[0,j]-mx[1,j]==1: idelt[:ntrms,j]=idel0[:ntrms,j]
        if nops==4:
            if nct==1:
                idelt[:ntrms,lrt+ic1]=lrt+ic2; idelt[:ntrms,lrt+ic3]=lrt+ic4; reslt[:ntrms]*=0.5
            else:
                idelt[:ntrms,lrt+ic1]=lrt+ic4; idelt[:ntrms,lrt+ic3]=lrt+ic2; reslt[:ntrms]*=-0.5/N
        elif nops==2:
            idelt[:ntrms,lrt+ic1]=lrt+ic2
        for il in range(ntrms):
            for im in range(N*B+1,lrt,2): idelt[il,im]=im-1
            for in_ in range(lng-N*B-1,lng-llt,-2): idelt[il,in_]=in_-1
        for nt in range(ntrms):
            for lt in range(lng):
                ic=lt
                for _ in range(2*lng):
                    nx=idelt[nt,ic]
                    if nx==0: break
                    nx2=idelt[nt,nx]
                    if nx2==0: break
                    idelt[nt,ic]=nx2; idelt[nt,nx]=0
                    if nx2==ic: idelt[nt,ic]=0; reslt[nt]*=N; break
        # ── Baryon epsilon reduction (BREDCE) ──
        # For B≠0: reduce epsilon pairs to antisymmetrized delta products
        if B:
            # Generate permutation table for N>2 (EPSB)
            if N > 2:
                from qcdf import gprbig_array, psign as _psign, MXNP
                nmin_ep = np.ones(MXNP, dtype=int)
                nmax_ep = np.ones(MXNP, dtype=int) * (N - 1)
                nprms_ep, ibrpm = gprbig_array(nmin_ep, nmax_ep, N-1, 1, 0, 9523)
                for is1 in range(nprms_ep):
                    ibrpm[is1, N-1] = _psign(N-1, ibrpm[is1, :N-1].tolist())
            else:
                nprms_ep = 0  # N=2: no extra permutation terms

            for _iter in range(2*B+1):
                # CNTRCT
                for nt in range(ntrms):
                    for lt in range(lng):
                        ic=lt
                        for __ in range(2*lng):
                            nx=idelt[nt,ic]
                            if nx==0: break
                            nx2=idelt[nt,nx]
                            if nx2==0: break
                            idelt[nt,ic]=nx2; idelt[nt,nx]=0
                            if nx2==ic: idelt[nt,ic]=0; reslt[nt]*=N; break

                # BREDCE: reduce epsilon pairs to deltas
                ntpr = ntrms; has_terms = False
                for nt in range(ntrms):
                    for lt in range(lng):
                        if idelt[nt, lt] != 0:
                            has_terms = True
                            inb1 = lt; inb2 = idelt[nt, lt]
                            reslt[nt] *= icycl[inb1] * icycl[inb2]
                            # Follow LNKB to get the other N-1 indices of each epsilon
                            nb1 = [0]*(N-1); nb2 = [0]*(N-1)
                            c1, c2 = inb1, inb2
                            for i1 in range(N-1):
                                nb1[i1] = lnkb[c1]; nb2[i1] = lnkb[c2]
                                c1 = lnkb[c1]; c2 = lnkb[c2]
                            # Identity permutation: replace contracted pair with delta
                            idelt[nt, lt] = 0
                            for j1 in range(N-1):
                                idelt[nt, nb2[j1]] = nb1[j1]
                            # Additional permutations (only for N>2)
                            for j2 in range(1, nprms_ep):
                                if ntpr >= MXTRM:
                                    raise OverflowError(
                                        "colour contraction exceeded MXTRM "
                                        "terms; the result would be silently "
                                        "wrong. Raise MXTRM in BOTH qcdf_opt.py "
                                        "and qcdf_kernels.py (they must match), "
                                        "or reduce LPN.")
                                reslt[ntpr] = reslt[nt] * float(ibrpm[j2, N-1])
                                idelt[ntpr, :lng] = idelt[nt, :lng]
                                for j4 in range(N-1):
                                    idelt[ntpr, nb2[ibrpm[j2, j4]-1]] = nb1[j4]
                                ntpr += 1
                            break  # one epsilon pair per term per iteration
                ntrms = ntpr
                if not has_terms: break

        # Sum result (avoid numpy overhead on small array)
        ct = 0.0
        for i in range(ntrms): ct += reslt[i]
        total += ct; ntrms=ntrms0
    return total

# ─── Hamiltonian element (table-driven vertices) ─────────────────
_mx_work = np.zeros((4, MXLNG), dtype=np.int32)
_ltc_work = np.zeros((4, MXP), dtype=np.int32)

def hamqcd(lrt, llt, matrt, matlt, N, NF, B, K, selfen, cbreak):
    matltc = _ltc_work; conj_state_nb(llt, matlt, matltc)
    elem0=np.zeros(NF); elem=0.0
    mqn=[[], [], [], []]
    for i in range(lrt):
        idx=0 if matrt[0,i]==1 else 1; m=matrt[2,i]
        if m not in mqn[idx]: mqn[idx].append(m)
    for j in range(llt):
        idx=2 if matltc[0,j]==-1 else 3; m=matltc[2,j]
        if m not in mqn[idx]: mqn[idx].append(m)
    nbr=sum(matrt[0,i]==1 for i in range(lrt)); ndr=sum(matrt[1,i]==1 for i in range(lrt))
    nbl=sum(matlt[0,i]==1 for i in range(llt)); ndl=sum(matlt[1,i]==1 for i in range(llt))
    def opmtch(a,b,c,d): return a<=nbr and b<=nbl and c<=ndr and d<=ndl and a+nbl==b+nbr and c+ndl==d+ndr
    def mi(u,v): return 2*u*u+3*v*v+u+v-1

    def do_mx(nops):
        mx=_mx_work; mx[:,:lrt]=matrt[:,:lrt]; mx[:,lrt+nops:lrt+nops+llt]=matltc[:,:llt]
        return mx

    # 4-point vertices: (nhb,nhba,nhd,nhda, row0,row1, ic, ka_fn, flav_fn, coeff)
    V4 = [
        (2,2,0,0, [-1,-1,1,1],[0,0,0,0], (0,3,1,2), lambda k:(k[3]-k[1],k[2]-k[0]), lambda a,b:(a,b,a,b), -0.5),
        (0,0,2,2, [0,0,0,0],[-1,-1,1,1], (2,1,3,0), lambda k:(k[1]-k[3],k[0]-k[2]), lambda a,b:(a,b,a,b), -0.5),
        (1,2,0,1, [0,-1,1,1],[1,0,0,0], (1,3,0,2), lambda k:(k[3]+k[0],k[2]-k[1]), lambda a,b:(a,b,b,a), 1.0),
        (0,1,1,2, [0,0,0,1],[-1,1,1,0], (1,0,2,3), lambda k:(k[0]-k[2],-k[3]-k[1]), lambda a,b:(a,b,a,b), 1.0),
        (2,1,1,0, [0,-1,-1,1],[-1,0,0,0], (1,0,2,3), lambda k:(k[0]+k[2],k[1]-k[3]), lambda a,b:(a,b,a,b), 1.0),
        (1,0,2,1, [0,0,0,-1],[-1,-1,1,0], (3,1,2,0), lambda k:(k[1]-k[2],k[0]+k[3]), lambda a,b:(a,b,b,a), 1.0),
    ]
    for nhb,nhba,nhd,nhda, r0,r1, (c1,c2,c3,c4), ka_fn, fl_fn, coeff in V4:
        if not opmtch(nhb,nhba,nhd,nhda): continue
        mx=do_mx(4); mx[0,lrt:lrt+4]=r0; mx[1,lrt:lrt+4]=r1; lng=4+lrt+llt
        mids=[mi(mx[0,lrt+j],mx[1,lrt+j]) for j in range(4)]
        h=0.0
        for l1 in range(len(mqn[mids[0]])):
            for l2 in range(len(mqn[mids[1]])):
                for l3 in range(len(mqn[mids[2]])):
                    for l4 in range(len(mqn[mids[3]])):
                        ks=[mqn[mids[0]][l1], mqn[mids[1]][l2], mqn[mids[2]][l3], mqn[mids[3]][l4]]
                        for j in range(4): mx[2,lrt+j]=ks[j]
                        if sum(mx[2,lrt+j]*(mx[0,lrt+j]+mx[1,lrt+j]) for j in range(4))!=0: continue
                        ka,kb=ka_fn(ks); bk=brack(ka,kb)
                        if bk==0.0: continue
                        for f1 in range(1,NF+1):
                            for f2 in range(1,NF+1):
                                fl=fl_fn(f1,f2)
                                for j in range(4): mx[3,lrt+j]=fl[j]
                                h+=coeff*clfact(c1,c2,c3,c4,mx,lng,llt,lrt,4,N,NF,B,K)*bk
        elem+=h

    # H7 (exchange, two flavor structures)
    if opmtch(1,1,1,1):
        mx=do_mx(4); mx[0,lrt:lrt+4]=[0,0,-1,1]; mx[1,lrt:lrt+4]=[-1,1,0,0]; lng=4+lrt+llt
        mids=[mi(mx[0,lrt+j],mx[1,lrt+j]) for j in range(4)]
        h=0.0
        for l1 in range(len(mqn[mids[0]])):
            for l2 in range(len(mqn[mids[1]])):
                for l3 in range(len(mqn[mids[2]])):
                    for l4 in range(len(mqn[mids[3]])):
                        ks=[mqn[mids[0]][l1], mqn[mids[1]][l2], mqn[mids[2]][l3], mqn[mids[3]][l4]]
                        for j in range(4): mx[2,lrt+j]=ks[j]
                        if sum(mx[2,lrt+j]*(mx[0,lrt+j]+mx[1,lrt+j]) for j in range(4))!=0: continue
                        bka=brack(ks[3]-ks[2],ks[1]-ks[0]); bkc=brack(ks[3]+ks[1],-ks[0]-ks[2])
                        if bka!=0:
                            for f1 in range(1,NF+1):
                                for f2 in range(1,NF+1):
                                    mx[3,lrt:lrt+4]=[f1,f1,f2,f2]
                                    h-=clfact(1,3,2,0,mx,lng,llt,lrt,4,N,NF,B,K)*bka
                        if bkc!=0:
                            for f1 in range(1,NF+1):
                                for f2 in range(1,NF+1):
                                    mx[3,lrt:lrt+4]=[f1,f2,f1,f2]
                                    h+=clfact(2,3,1,0,mx,lng,llt,lrt,4,N,NF,B,K)*bkc
        elem+=h

    # H8, H9 (diagonal self-energy)
    for bd,ops0,ops1,c12,sgn in [(0,[-1,1],[0,0],(0,1),1),(1,[0,0],[-1,1],(1,0),-1)]:
        nhb,nhba,nhd,nhda=(1,1,0,0) if bd==0 else (0,0,1,1)
        if not opmtch(nhb,nhba,nhd,nhda): continue
        mx=do_mx(2); mx[0,lrt:lrt+2]=ops0; mx[1,lrt:lrt+2]=ops1; lng=2+lrt+llt
        m0,m1=mi(ops0[0],ops1[0]),mi(ops0[1],ops1[1])
        for l1 in range(len(mqn[m0])):
            for l2 in range(len(mqn[m1])):
                mx[2,lrt]=mqn[m0][l1]; mx[2,lrt+1]=mqn[m1][l2]
                if mx[2,lrt]*(mx[0,lrt]+mx[1,lrt])+mx[2,lrt+1]*(mx[0,lrt+1]+mx[1,lrt+1])!=0: continue
                k1,k2=mx[2,lrt],mx[2,lrt+1]
                if k1==0 or k1!=k2: continue
                for ifl in range(1,NF+1):
                    mx[3,lrt]=ifl; mx[3,lrt+1]=ifl
                    cf=clfact(c12[0],c12[1],0,0,mx,lng,llt,lrt,2,N,NF,B,K)
                    elem0[ifl-1]+=cf*(2.0+sgn*cbreak)/k1
                    if 0<k1<=99: elem+=cf*0.5*selfen[k1-1]
    return elem0, elem

# ─── Row-batched matrix construction ─────────────────────────────
# Global state for forked workers (set by parent before Pool creation)
_g_mstate = None
_g_mstinf = None
_g_numsta = None
_g_N = _g_NF = _g_B = _g_K = 0
_g_selfen = None
_g_cbreak = 1e-8
_g_norh = 0

def _init_worker(mstate, mstinf, numsta, N, NF, B, K, selfen, cbreak, norh):
    """Initialize global state in forked worker processes."""
    global _g_mstate, _g_mstinf, _g_numsta, _g_N, _g_NF, _g_B, _g_K, _g_selfen, _g_cbreak, _g_norh
    _g_mstate = mstate; _g_mstinf = mstinf; _g_numsta = numsta
    _g_N = N; _g_NF = NF; _g_B = B; _g_K = K
    _g_selfen = selfen; _g_cbreak = cbreak; _g_norh = norh

def _config_keys(mstate, mstinf, numsta):
    """One hashable key per state: its multiset of (type, momentum, flavour).

    Two states can have a non-zero norm element only if these agree exactly --
    the worker's own filter already demands it, and tools/colour_norm.py
    confirms the norm is block-diagonal in exactly this grouping.  Computing
    the key once per state, in O(n L), replaces an O(n^2 L^2) scan over pairs
    that almost all fail.
    """
    keys = []
    for s in range(numsta):
        loc = mstinf[s, 0] - 1
        L = mstinf[s, 1]
        keys.append(tuple(sorted(
            (int(mstate[loc, j]), int(mstate[loc + 2, j]), int(mstate[loc + 3, j]))
            for j in range(L))))
    return keys


def _worker_row_global(task):
    """Compute one row using global (forked) state — no pickling of arrays.

    ``task`` is either a row index (visit every column at or after it) or a
    ``(row, columns)`` pair naming exactly the columns that can be non-zero.
    """
    if isinstance(task, tuple):
        nstrt, cols = task
    else:
        nstrt, cols = task, None
    mstate=_g_mstate; mstinf=_g_mstinf; numsta=_g_numsta
    N=_g_N; NF=_g_NF; B=_g_B; K=_g_K; selfen=_g_selfen; cbreak=_g_cbreak; norh=_g_norh

    lr=mstinf[nstrt,0]-1; lstr=mstinf[nstrt,1]; matrt=mstate[lr:lr+4,:MXP].copy()
    rh0=np.zeros((NF,numsta)); rh=np.zeros(numsta); rn=np.zeros(numsta)
    for nstlt in (range(nstrt, numsta) if cols is None else cols):
        ll=mstinf[nstlt,0]-1; lstl=mstinf[nstlt,1]; matlt=mstate[ll:ll+4,:MXP].copy()
        used=np.zeros(lstl,dtype=bool); sb=sd=0
        for ir in range(lstr):
            for il in range(lstl):
                if used[il]: continue
                if matrt[0,ir]==matlt[0,il] and matrt[2,ir]==matlt[2,il] and matrt[3,ir]==matlt[3,il]:
                    if matlt[0,il]==1: sb+=1
                    else: sd+=1
                    used[il]=True; break
        nbr=np.sum(matrt[0,:lstr]==1); ndr=np.sum(matrt[1,:lstr]==1)
        nbl=np.sum(matlt[0,:lstl]==1); ndl=np.sum(matlt[1,:lstl]==1)
        if abs(nbr-sb)+abs(nbl-sb)+abs(ndr-sd)+abs(ndl-sd)>4: continue
        if norh==0:
            if abs(nbr-sb)+abs(nbl-sb)+abs(ndr-sd)+abs(ndl-sd)==0:
                mx=np.zeros((4,MXLNG),dtype=np.int32); ltc=np.zeros((4,MXP),dtype=np.int32)
                conj_state_nb(lstl,matlt,ltc); mx[:,:lstr]=matrt[:,:lstr]; mx[:,lstr:lstr+lstl]=ltc[:,:lstl]
                rn[nstlt]=clfact(0,0,0,0,mx,lstr+lstl,lstl,lstr,0,N,NF,B,K)
        else:
            e0,e=hamqcd(lstr,lstl,matrt,matlt,N,NF,B,K,selfen,cbreak)
            rh0[:,nstlt]=e0; rh[nstlt]=e
    return nstrt, rh0, rh, rn

def _row_tasks(norh, mstate, mstinf, numsta):
    """Work items for one matrix build, ordered longest-first for balance.

    The norm couples only states with identical parton content, so group by
    that and hand each row just its own group's columns.  For the Hamiltonian a
    four-point vertex allows up to four mismatches, so every pair has to be
    visited and the full upper triangle stands -- and since row *r* then costs
    ``numsta - r``, contiguous chunks would be badly unbalanced.
    """
    if norh == 0:
        groups = {}
        for s, key in enumerate(_config_keys(mstate, mstinf, numsta)):
            groups.setdefault(key, []).append(s)
        tasks = [(r, np.array([c for c in groups[k] if c >= r], dtype=np.int64))
                 for k in groups for r in groups[k]]
        tasks.sort(key=lambda t: -len(t[1]))
        return tasks
    return sorted(range(numsta), key=lambda r: -(numsta - r))


@njit(nogil=True, cache=True)
def _mirror_upper(a):
    """Copy the strict upper triangle of ``a`` into its lower triangle."""
    n = a.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            a[j, i] = a[i, j]


def _symmetrize(a):
    """Mirror an upper-triangular build into a full symmetric matrix, in place.

    This was ``np.triu(a) + np.triu(a, 1).T``, which is correct but touches
    three extra n x n arrays.  The norm is built at the *pre*-weeding size, so
    at 2K=39 that turned a 22 GB matrix into 88 GB of live allocation and was
    the binding memory wall -- ahead of anything about speed.  The workers fill
    only ``j >= i``, so mirroring in place is the same result with no
    temporaries, and being a pure copy it cannot perturb a value (which the
    bit-exactness tests in ``tests/test_kernels.py`` would catch).
    """
    _mirror_upper(a)
    return a


def build_norm_blocks(mstate, mstinf, numsta, N, NF, B, K, labels,
                      nthreads=1):
    """The norm as a list of dense blocks, never allocating ``numsta**2``.

    The norm is *exactly* block diagonal in parton configuration -- off-block
    entries are ``0.0``, not 1e-16, verified by explicit colour enumeration --
    so the dense array was mostly zeros that had to be allocated anyway.  It is
    also built at the **pre**-weeding size, which makes it the binding memory
    constraint of the whole solver: 8.6 GB at 2K=37, 51.5 GB at 2K=41, against
    7.6 MB and 27.9 MB for the blocks alone.

    Nothing about the arithmetic changes.  ``kern.norm_row`` is called with the
    same arguments over the same column sets that ``_row_tasks(norh=0)`` already
    groups by configuration key; only the destination differs, so the stored
    values are bit-identical to slicing the dense build.  Each thread keeps one
    full-length scratch row (0.7 MB at 2K=41) and copies out just its block's
    columns.

    ``labels`` comes from :func:`dlcq.read_python.config_block_labels`, which
    derives the grouping from the state list in O(n L) without touching a
    matrix.

    Returns ``[(members, block), ...]`` with ``members`` the ascending state
    indices of that block and ``block`` its dense symmetric submatrix.
    """
    import qcdf_kernels as kern

    nprms, ibrpm = kern.epsb_table(N)
    labels = np.asarray(labels)[:numsta]
    order = np.argsort(labels, kind="stable")
    groups = np.split(order, np.flatnonzero(np.diff(labels[order])) + 1)
    groups = [np.ascontiguousarray(np.sort(g)) for g in groups if g.size]

    blocks = [None] * len(groups)
    tls = threading.local()

    def run(gi):
        members = groups[gi]
        b = members.size
        sc = getattr(tls, "buf", None)
        if sc is None:
            sc = tls.buf = (kern.Scratch(NF, K), np.zeros(numsta))
        scratch, row = sc
        M = np.zeros((b, b))
        for i in range(b):
            cols = members[i:]
            row[cols] = 0.0
            kern.norm_row(members[i], cols, mstate, mstinf, N, NF, B, K, row,
                          nprms, ibrpm, *scratch.norm_args())
            M[i, i:] = row[cols]
        _mirror_upper(M)
        blocks[gi] = M

    if groups:
        # Compile on this thread first; JIT holds the GIL, so racing workers
        # into it would serialize startup (same reason as _build_threaded).
        run(0)
        rest = range(1, len(groups))
        if nthreads > 1 and len(groups) > 2:
            with ThreadPoolExecutor(nthreads) as ex:
                for _ in ex.map(run, rest):
                    pass
        else:
            for gi in rest:
                run(gi)
    return list(zip(groups, blocks))


def build_ham_sparse(mstate, mstinf, numsta, N, NF, B, K, selfen, cbreak,
                     nthreads=1, drop=0.0):
    """The Hamiltonian as CSR, never allocating ``numsta**2``.

    Same traversal and same kernel as the dense build -- ``ham_row`` writes the
    upper triangle of one row -- so the stored values are bit-identical; only
    the zeros are not kept.  Each thread reuses one length-``numsta`` row buffer
    and compresses it before moving on, which is what removes the n^2 array:
    1.86 GB at 2K=41 for 6.85 M non-zeros, i.e. 82 MB as CSR.

    Density falls steadily with K -- 26.05% at 2K=21, 10.68% at 29, 8.52% at 31,
    2.95% at 41 -- because ``n`` grows 1.58x per +2 while ``nnz/row`` grows only
    1.25x.  That is what makes the sparse form worth its overhead, and it is
    also why the sparse *eigensolver* only pays here: on a dense operator every
    ARPACK matvec is an O(n^2) GEMV, and measured, ``eigsh`` on this matrix
    sparse is 3.4 s against 123 s for a dense ``eigh(subset_by_index)``.

    ``drop`` discards |value| <= drop.  It defaults to 0.0 -- exact zeros only.
    A threshold here would silently change matrix elements, and
    ``tests/test_kernels.py`` demands they do not move.

    The buffer must be zeroed per row: ``ham_row`` leaves entries that fail its
    ``|delta| <= 4`` pre-filter untouched rather than writing zero.
    """
    import qcdf_kernels as kern
    import scipy.sparse as sp

    nprms, ibrpm = kern.epsb_table(N)
    rows_out = [None] * numsta
    tls = threading.local()

    def run(r):
        buf = getattr(tls, "buf", None)
        if buf is None:
            buf = tls.buf = (kern.Scratch(NF, K), np.zeros(numsta),
                             np.zeros((NF, numsta)))
        scratch, row, h0_row = buf
        row[r:] = 0.0
        kern.ham_row(r, numsta, mstate, mstinf, N, NF, B, K, selfen, cbreak,
                     h0_row, row, nprms, ibrpm, *scratch.ham_args())
        seg = row[r:]
        nz = np.flatnonzero(np.abs(seg) > drop) if drop else np.flatnonzero(seg)
        rows_out[r] = (nz + r, seg[nz].copy())

    if numsta:
        run(0)                       # compile on this thread; JIT holds the GIL
        rest = range(1, numsta)
        if nthreads > 1 and numsta > 2:
            with ThreadPoolExecutor(nthreads) as ex:
                for _ in ex.map(run, rest):
                    pass
        else:
            for r in rest:
                run(r)

    indptr = np.zeros(numsta + 1, dtype=np.int64)
    for r, (cols, _vals) in enumerate(rows_out):
        indptr[r + 1] = indptr[r] + cols.size
    indices = np.empty(indptr[-1], dtype=np.int64)
    data = np.empty(indptr[-1])
    for r, (cols, vals) in enumerate(rows_out):
        indices[indptr[r]:indptr[r + 1]] = cols
        data[indptr[r]:indptr[r + 1]] = vals
    upper = sp.csr_matrix((data, indices, indptr), shape=(numsta, numsta))

    # Mirror, without double-counting the diagonal the upper triangle already
    # holds.  The dense path does the same thing via _symmetrize.
    diag = sp.dia_matrix((upper.diagonal()[None, :], [0]),
                         shape=(numsta, numsta))
    return (upper + upper.T - diag).tocsr()


def _build_threaded(norh, mstate, mstinf, numsta, N, NF, B, K, selfen, cbreak,
                    nthreads, want_h0=True):
    """Threaded build over the nogil kernels in ``qcdf_kernels``.

    Each row is written straight into the destination matrix, so nothing is
    pickled and nothing is copied back; the workers touch disjoint rows, and
    only the upper triangle is filled, so there is no write conflict and no
    lock.  Per-thread :class:`~qcdf_kernels.Scratch` replaces the module-global
    buffers, which were safe only because ``fork`` gave each process a copy.
    """
    import qcdf_kernels as kern

    nprms, ibrpm = kern.epsb_table(N)
    # Allocate only the matrix this path fills.  Both were allocated
    # unconditionally before, so a norm build reserved (NF+1) n^2 it never
    # touched -- lazily faulted, so mostly address space rather than RSS, but
    # there is no reason to reserve it.
    if norh == 0:
        ham0 = np.zeros((NF, 0, 0))
        ham = np.zeros((0, 0))
        hnorm = np.zeros((numsta, numsta))
    else:
        # ``want_h0=False`` gives each thread a throwaway row instead of a slice
        # of an (NF, n, n) array.  H0 = D*N with D fixed by parton content, so a
        # blockwise Z makes Z^T H0 Z exactly diag(D) and the matrix is redundant
        # -- see dlcq.read_python.free_energy_diagonal.  At 2K=41 that is
        # 1.86 GB per flavour not allocated, and another 1.86 GB of projection
        # not performed.  The kernel still needs somewhere to write.
        ham0 = (np.zeros((NF, numsta, numsta)) if want_h0
                else np.zeros((NF, 0, 0)))
        ham = np.zeros((numsta, numsta))
        hnorm = np.zeros((0, 0))
    tasks = _row_tasks(norh, mstate, mstinf, numsta)

    tls = threading.local()

    def run(task):
        sc = getattr(tls, "scratch", None)
        if sc is None:
            sc = tls.scratch = kern.Scratch(NF, K)
        if norh == 0:
            r, cols = task
            kern.norm_row(r, cols, mstate, mstinf, N, NF, B, K, hnorm[r],
                          nprms, ibrpm, *sc.norm_args())
        else:
            r = task
            if want_h0:
                h0_row = ham0[:, r, :]
            else:
                h0_row = getattr(tls, "h0_row", None)
                if h0_row is None or h0_row.shape[1] != numsta:
                    h0_row = tls.h0_row = np.zeros((NF, numsta))
            kern.ham_row(r, numsta, mstate, mstinf, N, NF, B, K, selfen,
                         cbreak, h0_row, ham[r], nprms, ibrpm,
                         *sc.ham_args())

    if tasks:
        # Compile on this thread first: JIT compilation holds the GIL, so
        # letting several workers race into it would serialize the startup.
        run(tasks[0])
        rest = tasks[1:]
        if nthreads > 1 and len(rest) > 1:
            with ThreadPoolExecutor(nthreads) as ex:
                for _ in ex.map(run, rest):
                    pass
        else:
            for t in rest:
                run(t)

    if norh == 0:
        return ham0, ham, _symmetrize(hnorm)
    if want_h0:
        for ifl in range(NF):
            ham0[ifl] = _symmetrize(ham0[ifl])
    return ham0, _symmetrize(ham), hnorm


def _build_process(norh, mstate, mstinf, numsta, N, NF, B, K, selfen, cbreak,
                   ncpus):
    """The original process-pool build over the interpreted reference routines.

    Kept as the cross-check for ``_build_threaded`` -- ``tests/test_kernels.py``
    asserts the two produce identical matrices -- and as a fallback if numba is
    unavailable.
    """
    if norh == 0:
        ham0 = np.zeros((NF, 0, 0))
        ham = np.zeros((0, 0))
        hnorm = np.zeros((numsta, numsta))
    else:
        ham0 = np.zeros((NF, numsta, numsta))
        ham = np.zeros((numsta, numsta))
        hnorm = np.zeros((0, 0))
    tasks = _row_tasks(norh, mstate, mstinf, numsta)
    if norh == 0:                       # this path wants plain lists
        tasks = [(r, cols.tolist()) for r, cols in tasks]

    if ncpus > 1 and numsta > 10:
        # Use initializer to share state via fork — no pickling of arrays!
        with Pool(ncpus, initializer=_init_worker,
                  initargs=(mstate, mstinf, numsta, N, NF, B, K, selfen, cbreak, norh)) as pool:
            # chunksize=1: dynamic scheduling, since per-row cost varies a lot.
            results = pool.map(_worker_row_global, tasks, chunksize=1)
    else:
        _init_worker(mstate, mstinf, numsta, N, NF, B, K, selfen, cbreak, norh)
        results = [_worker_row_global(t) for t in tasks]

    for r,rh0,rh,rn in results:
        if norh==0: hnorm[r,r:numsta]+=rn[r:]; hnorm[r:numsta,r]=hnorm[r,r:numsta]
        else:
            ham0[:,r,r:numsta]+=rh0[:,r:]; ham0[:,r:numsta,r]=ham0[:,r,r:numsta]
            ham[r,r:numsta]+=rh[r:]; ham[r:numsta,r]=ham[r,r:numsta]
    return ham0, ham, hnorm


#: Default parallel backend.  ``"thread"`` runs the nogil numba kernels in a
#: thread pool; ``"process"`` runs the interpreted reference routines under
#: multiprocessing.  Override per call, or globally with ``QCDF_BACKEND``.
DEFAULT_BACKEND = "thread"


def build_matrices(norh, mstate, mstinf, numsta, N, NF, B, K, selfen, cbreak,
                   ncpus, backend=None, want_h0=True):
    """Build the norm (``norh=0``) or Hamiltonian (``norh=1``) matrix.

    ``backend`` selects how the work is spread; both backends return identical
    matrices, bit for bit.  See ``docs/performance.md``.

    ``want_h0=False`` skips storing the free part, which a caller using the
    ``H0 = D N`` shortcut does not need; the returned ``ham0`` is then empty.
    Only the threaded backend honours it -- the process backend is the
    bit-exactness reference and is left alone.
    """
    backend = backend or os.environ.get("QCDF_BACKEND") or DEFAULT_BACKEND
    if backend not in ("thread", "process"):
        raise ValueError(f"unknown backend {backend!r}")
    t0 = time.time()
    label = "norm" if norh == 0 else "ham"
    log.info(f"  {label} ({numsta}², {ncpus} {backend}s)...")
    if backend == "thread":
        ham0, ham, hnorm = _build_threaded(norh, mstate, mstinf, numsta, N, NF,
                                           B, K, selfen, cbreak, ncpus,
                                           want_h0=want_h0)
    else:
        ham0, ham, hnorm = _build_process(norh, mstate, mstinf, numsta, N, NF,
                                          B, K, selfen, cbreak, ncpus)
    log.info(f"  {label}: {time.time()-t0:.2f}s")
    return ham0, ham, hnorm

# ─── State generation ────────────────────────────────────────────

def _sector_list(params, lnb, lnd):
    """The ``(nbrb, nbrd, nmes, kbrb, kbrd, kmes)`` sectors ``qcdsta`` visits.

    Factored out of ``qcdf.qcdsta``'s four nested loops verbatim, so the
    compiled driver walks them in exactly the same order.  That order is the
    state order, and the state order is load-bearing: ``WEEDR``/``WEEDR2`` drop
    linearly dependent states by index.
    """
    N, NF, B, K = params.N, params.NF, params.B, params.K
    sectors = []
    lmbf = lnb - N * B
    for lmb in range(lmbf + 1):
        nxbmx = lmb // N if params.IBAB != 0 else 0
        for nxb in range(nxbmx + 1):
            nbar = B + 2 * nxb
            nmes = lmb - N * nxb

            kbmn = minmom(N * (B + nxb), N, NF)
            kdmn = minmom(N * nxb, N, NF)
            kmmn = minmom(nmes, N, NF)
            kbarmn = kbmn + kdmn
            kmesmn = 2 * kmmn
            kbarmx = K - kmesmn

            if nbar == 0:
                kbarmn, kbarmx = 0, 0
            if nmes == 0:
                kbarmn, kbarmx = K, K

            if not ((nmes != 0 or nbar != 0)
                    and (N * B + lmb) <= lnb and lmb <= lnd):
                continue

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
                    if iodev(kmes_val) != 1:
                        continue
                    kbrboe = iodev(kbrb)
                    kbrdoe = iodev(kbrd)
                    if ((noe == 1 and kbrboe == 1 and kbrdoe == 1)
                            or (noe == -1 and kbrboe == nbrboe
                                and kbrdoe == nbrdoe)):
                        sectors.append((nbrb, nbrd, nmes, kbrb, kbrd, kmes_val))
    return sectors


def qcdsta_fast(params, states, perm, flavtab, ncpus=1):
    """Compiled equivalent of ``qcdf.qcdsta``.

    Produces **element-identical** ``mstate``/``mstinf`` to the interpreted
    version -- asserted in ``tests/test_states.py`` -- because state order
    feeds weeding, which feeds the basis, which feeds the spectrum.

    The permutation tables (``prmx``/``prmm``) and flavour tables
    (``ibarfl``/``imesfl``) are built once per run and stay in Python; only
    ``brmsgn``, which is 98% of the cost, is compiled.
    """
    import qcdf_states as ks

    N, NF, B, K = params.N, params.NF, params.B, params.K
    states.numsta = 0

    lnb, lnd = lpnsub(N, NF, B, K)
    log.info(f"  LPNSUB: LNB={lnb}, LND={lnd}")
    nptmx = params.LPN if (params.LPN != 0 and params.LPN < lnb) else lnb
    if params.LPN > (lnb + lnd):
        params.LPN = lnb + lnd

    # The momentum-partition table is the one array whose size is not knowable
    # in advance and which grows ~1.6x per +2 in K (62k rows at 2K=35, 140k at
    # 2K=39, ~1.5M projected at 2K=48).  Rather than hard-code an ever-larger
    # LKPXMX, size it to the run: start from whatever PermTables allocated and
    # double until prmx fits.  This removes LKPXMX as a ceiling entirely.
    # ``kpxloc``/``kpmloc`` are indexed by total momentum, up to K, but
    # ``PermTables`` sizes their momentum axis with the fixed ``MXKMX = 100``.
    # Numba does not bounds-check, so beyond 2K = 99 ``prmx_nb`` and ``prmm_nb``
    # wrote past the end of the allocation and smashed the heap: generation
    # *completed*, returned a plausible state count, and the process then died
    # with "corrupted size vs. prev_size" whenever the allocator next walked
    # its metadata -- at interpreter shutdown, typically, far from the cause.
    # Measured: clean at 2K=99, corrupt from 2K=119, exactly where K first
    # exceeds MXKMX.
    #
    # Sized to the run, like ``kpx`` below, rather than by raising the constant,
    # so there is no new ceiling to rediscover later.
    need = K + 2
    if perm.kpxloc.shape[2] < need:
        s = perm.kpxloc.shape
        perm.kpxloc = np.zeros((s[0], s[1], need, s[3]), dtype=np.int64)
    if perm.kpmloc.shape[0] < need:
        perm.kpmloc = np.zeros((need, perm.kpmloc.shape[1]), dtype=np.int64)

    nuterm = np.zeros(ks.MXNP, dtype=np.int64)
    kpx, kpxloc = perm.kpx, perm.kpxloc
    while True:
        numprm = ks.prmx_nb(K, nptmx, N * NF, perm.IX0, perm.IXN,
                            kpx, kpxloc, nuterm)
        if numprm != ks.OVERFLOW_PERM:
            break
        kpx = np.zeros((kpx.shape[0] * 2, ks.MXNP), dtype=np.int64)
    perm.kpx, perm.numprm = kpx, numprm

    # kpm needs only K^2/2 rows; the default LKPXMX sizing was vastly generous.
    kpm = np.zeros((K * K // 2 + 2, 2), dtype=np.int64)
    nprmm = ks.prmm_nb(K, kpm, perm.kpmloc)
    if nprmm == ks.OVERFLOW_PERM:            # cannot happen; K^2/2 is exact
        raise RuntimeError("kpm sizing bound violated")
    perm.kpm, perm.nprmm = kpm, nprmm

    ibarfl(params, flavtab)
    imesfl(params, flavtab)

    sectors = _sector_list(params, lnb, lnd)
    log.info(f"  {len(sectors)} momentum sectors")

    rmq = np.ascontiguousarray(params.rmq, dtype=np.float64)
    iflv = np.ascontiguousarray(params.iflv, dtype=np.int64)
    tables = (perm.kpx, perm.kpxloc, perm.kpm, perm.kpmloc, perm.IX0, perm.IXN,
              flavtab.mbarfl, flavtab.mmesfl, flavtab.ibfdim, flavtab.imfdim)

    scratch = ks.Scratch(NF)
    ns = 0
    for (nbrb, nbrd, nmes, kbrb, kbrd, kmes) in sectors:
        while True:
            r = ks.brmsgn_nb(N, NF, nbrb, nbrd, nmes, kbrb, kbrd, kmes,
                             params.LPN, params.cutoff, K, rmq, iflv,
                             *tables,
                             states.mstate, states.mstinf, ns,
                             *scratch.args())
            if r == ks.OVERFLOW_PERM:
                # Required size is only knowable at runtime and grows ~1.6x
                # per +2 in K; grow and redo the sector from the same index.
                scratch = scratch.grown(NF)
                continue
            if r == ks.OVERFLOW_STATES:
                # Same treatment as kpx: the state count grows ~1.6x per +2 in
                # K, so grow rather than making the caller guess NMSTMX.  The
                # sector is redone from ns, overwriting its partial writes.
                old = states.mstinf.shape[0]
                new = max(old * 2, 1024)
                mstinf = np.zeros((new, 8), dtype=np.int64)
                mstinf[:old] = states.mstinf
                mstate = np.zeros((4 * new, ks.MXNP), dtype=np.int64)
                mstate[:4 * old] = states.mstate
                states.mstinf, states.mstate = mstinf, mstate
                continue
            ns += r
            break

    states.numsta = ns
    log.info(f"  QCDSTA complete: {ns} states generated")


# ─── Weeding ─────────────────────────────────────────────────────
def weed(hnorm, mstinf, numsta):
    eps=1e-4; loc=0
    # Step 1: linear weeding (proportional rows)
    while loc<numsta-1:
        drops=[]
        for j in range(loc+1,numsta):
            if abs(hnorm[j,loc])>eps:
                r=hnorm[j,loc]/hnorm[loc,loc]
                if np.all(np.abs(r*hnorm[loc,:numsta]-hnorm[j,:numsta])<eps): drops.append(j)
        for d in sorted(drops,reverse=True):
            hnorm=np.delete(np.delete(hnorm,d,0),d,1); mstinf=np.delete(mstinf,d,0); numsta-=1
        loc+=1
    log.info(f"    linear weeding → {numsta}")
    # Step 2: eigenvalue weeding (remove non-positive-definite directions)
    for iteration in range(2000):
        w,z=eigh(hnorm[:numsta,:numsta])
        # Count eigenvalues that are zero OR negative (not just near-zero)
        nzer=int(np.sum(w < eps))
        if nzer==0: break
        drops=[]; used=set()
        for i in range(nzer):
            for j in range(numsta-1,-1,-1):
                if abs(z[j,i])>eps and j not in used:
                    drops.append(j)
                    for k in range(numsta):
                        if abs(z[k,i])>eps: used.add(k)
                    break
        if not drops: break
        for d in sorted(set(drops),reverse=True):
            hnorm=np.delete(np.delete(hnorm,d,0),d,1); mstinf=np.delete(mstinf,d,0); numsta-=1
        log.info(f"    iter {iteration}: dropped {len(drops)} → {numsta}")
    return hnorm, mstinf, numsta

# ─── MAIN ────────────────────────────────────────────────────────
def main(input_file=None):
    log.info("="*60+"\n DLCQ QCD 1+1D (optimized)\n"+"="*60)
    ncpus=int(os.environ.get('QCDF_NCPUS',min(cpu_count(),8)))
    cfg=json.load(open(input_file)) if input_file else {}
    p=Params(); p.N=cfg.get('N',3); p.NF=cfg.get('NF',1); p.B=cfg.get('B',0)
    p.rlamb=cfg.get('rlamb',0.3325); p.K=cfg.get('K',6); p.cutoff=cfg.get('cutoff',-1.0)
    p.LPN=cfg.get('LPN',0)
    for i,v in enumerate(cfg.get('rmq',[1.0])[:p.NF]): p.rmq[i]=v
    for i,v in enumerate(cfg.get('iflv',[0])[:p.NF]): p.iflv[i]=v
    log.info(f"  N={p.N} NF={p.NF} B={p.B} K={p.K} λ={p.rlamb}")

    selfen=compute_selfen(p.N)
    states=StateData(); perm=PermTables(); flav=FlavorTables()
    qcdsta(p, states, perm, flav)
    numsta=states.numsta; log.info(f"  {numsta} states before weeding")

    _,_,hnorm=build_matrices(0,states.mstate,states.mstinf,numsta,p.N,p.NF,p.B,p.K,selfen,p.cbreak,ncpus)
    hnorm,mstinf,numsta=weed(hnorm,states.mstinf[:numsta],numsta)
    log.info(f"  {numsta} states after weeding")

    w_n,z_n=eigh(hnorm[:numsta,:numsta]); Z=z_n/np.sqrt(w_n)[np.newaxis,:]
    ham0,ham,_=build_matrices(1,states.mstate,mstinf,numsta,p.N,p.NF,p.B,p.K,selfen,p.cbreak,ncpus)

    # Vectorized basis transform  (replaces O(N⁴) Python loop)
    hnu=Z.T@ham[:numsta,:numsta]@Z
    hnu0=np.array([np.diag(Z.T@ham0[ifl,:numsta,:numsta]@Z) for ifl in range(p.NF)])

    rlmsq=p.rlamb**2; hnu*=rlmsq
    for ifl in range(p.NF): hnu[np.diag_indices(numsta)]+=(1-rlmsq)*p.rmq[ifl]**2*hnu0[ifl]

    w_eig,_=eigh(hnu); w_mass=p.K*w_eig/2.0

    with open('qcdf_py.out','w') as f:
        f.write(f" N,NF,B,K,LAMBDA: {p.N:5d}{p.NF:5d}{p.B:5d}{p.K:5d} {p.rlamb:24.16e}\n")
        f.write(f" NUMSTA: {numsta}\n\n COMPUTED EIGENVALUES:\n\n")
        for i in range(0,numsta,3):
            f.write("".join(f"{v:24.16e}" for v in w_mass[i:min(i+3,numsta)])+"\n")

    log.info(f"\n{'='*60}\n MASS² (lowest 20):")
    for i in range(min(20,numsta)): log.info(f"  {i+1:4d}: {w_mass[i]:20.12e}")
    log.info(f"{'='*60}\n  → qcdf_py.out")

if __name__=='__main__':
    main(sys.argv[1] if len(sys.argv)>1 else None)
