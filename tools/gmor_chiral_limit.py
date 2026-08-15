"""GMOR ratio WITH the extrapolation uncertainty, and the large-N behaviour.

The bar is the same ensemble the plots use: fit orders x contiguous sub-windows,
median and 68% half-width.  Without it a 0.5% residual cannot be called either
significant or not.
"""
import csv, sys
from collections import defaultdict
import numpy as np
sys.path.insert(0, ".")

def read(p):
    out = defaultdict(dict)
    for r in csv.DictReader(l for l in open(p) if not l.lstrip().startswith("#")):
        out[(int(r["N"]), float(r["mg"]))][int(r["K_code"])] = float(r["msq"])
    return out

def fit1k(ks, y, n):
    Kp = np.array(ks, float)/2; y = np.asarray(y)
    A = np.vstack([np.ones_like(Kp)]+[Kp**-(i+1.) for i in range(n)]).T
    sc = np.linalg.norm(A, axis=0)
    c, *_ = np.linalg.lstsq(A/sc, y, rcond=None)
    return c/sc

def spread(g):
    ks = sorted(g); y = [g[k] for k in ks]; n = len(ks); vals = []
    for order in (2,3,4):
        for i in range(n):
            for j in range(i+max(6,order+2), n+1):
                c = fit1k(ks[i:j], y[i:j], order)
                x = np.array(ks[i:j],float)/2
                A = np.vstack([np.ones_like(x)]+[x**-(t+1.) for t in range(order)]).T
                if np.max(np.abs(A@c-np.array(y[i:j])))/np.max(np.abs(y[i:j])) < 1e-3 and c[0] > 0:
                    vals.append(c[0])
    a = np.array(vals)
    return float(np.median(a)), float(0.5*(np.percentile(a,84)-np.percentile(a,16))), len(a)

def gmor(N, mg):
    c = (N*N-1.)/(2.*N*np.pi)
    return 2*np.pi*np.sqrt(c)*mg/np.sqrt(3.)/(mg**2+1./np.pi)

MGS=[0.1,0.05,0.025,0.0125]; NS=[2,3,4,5]
d={}
for N in NS: d.update(read(f"data/gmor_scan/improved_N{N}.csv"))
R={}
print("improved:  M^2(K->inf)/GMOR  +- ensemble (68%)   [n fits]")
for mg in MGS:
    cells=[]
    for N in NS:
        g=d.get((N,mg))
        if not g: cells.append("        --"); continue
        m0,sd,nf = spread(g); r=m0/gmor(N,mg); e=sd/gmor(N,mg)
        R[(N,mg)]=(r,e)
        cells.append(f"{r:7.4f}+-{e:6.4f}")
    print(f"  m/g={mg:<7}" + "  ".join(cells))

print("\n  deviation from 1, in units of its own bar (n_sigma):")
for mg in MGS:
    print(f"  m/g={mg:<7}" + "  ".join(
        f"N={N}:{(R[(N,mg)][0]-1)/R[(N,mg)][1]:+7.1f}" if (N,mg) in R and R[(N,mg)][1]>0
        else f"N={N}:    --" for N in NS))

print("\n  large-N at the smallest coupling: fit ratio = r_inf + c/N^2 on N=3,4,5")
for mg in MGS:
    ns=[N for N in (3,4,5) if (N,mg) in R]
    if len(ns)<3: continue
    x=np.array([1.0/N**2 for N in ns]); y=np.array([R[(N,mg)][0] for N in ns])
    e=np.array([R[(N,mg)][1] for N in ns])
    A=np.vstack([np.ones_like(x),x]).T
    W=np.diag(1/np.maximum(e,1e-9))
    c,*_=np.linalg.lstsq(W@A, W@y, rcond=None)
    # propagate: bar on r_inf from the fit covariance
    cov=np.linalg.inv(A.T@np.diag(1/np.maximum(e,1e-9)**2)@A)
    print(f"  m/g={mg:<7} r_inf = {c[0]:.4f} +- {np.sqrt(cov[0,0]):.4f}   (c={c[1]:+.4f})")
