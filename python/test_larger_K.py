#!/usr/bin/env python3
"""Compare Python vs Fortran structure functions at any K.

Usage: python test_larger_K.py K qcdf_KNN.ham qcdf_KNN.out [ncpus]

Strategy for Fortran: parse the post-weeding state list from qcdf.out,
match each to a Python state, extract the sub-matrices for exactly those
states, orthonormalize, diagonalize, and compute structure functions.
This avoids the Z matrix entirely (which uses pre-weeding indices).
"""
import os, sys, re
import numpy as np
from scipy.linalg import eigh
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

K = int(sys.argv[1]) if len(sys.argv) > 1 else 25
ham_file = sys.argv[2] if len(sys.argv) > 2 else 'qcdf.ham'
out_file = sys.argv[3] if len(sys.argv) > 3 else 'qcdf.out'
ncpus = int(sys.argv[4]) if len(sys.argv) > 4 else 8

os.environ['QCDF_NCPUS'] = str(ncpus)

from reproduce_figures import run_dlcq, extract_structure_function, mg_to_lambda, code_to_M_over_g
from qcdf_opt import build_matrices, compute_selfen
from qcdf import Params, StateData, PermTables, FlavorTables, qcdsta, compute_selfen as cs1

lam = mg_to_lambda(1.6)
print(f"K={K}, lam={lam:.4f}, ham={ham_file}, out={out_file}, ncpus={ncpus}")

# ═══════════════════════════════════════════════════════════════
# PYTHON: normal run
# ═══════════════════════════════════════════════════════════════
w_p, c_p, ns_p, mi_p, st_p, nso_p, Nm_p = run_dlcq(3, 1, 1, K, lam, ncpus=ncpus)
print(f"Python: {nso_p} states after weeding, {ns_p} orthonormal")

# ═══════════════════════════════════════════════════════════════
# FORTRAN: parse post-weeding states, match to Python, rebuild
# ═══════════════════════════════════════════════════════════════
with open(out_file, 'r') as f:
    olines = f.readlines()

# Parse post-weeding state list
f_states = []
i = 0
while i < len(olines):
    if 'STATE NUMBER:' in olines[i]:
        info = olines[i+2].split()
        lng = int(info[0])
        types = [int(x) for x in olines[i+4].split()[:lng]]
        moms = [int(x) for x in olines[i+6].split()[:lng]]
        f_states.append({'lng': lng, 'types': types, 'moms': moms})
    i += 1
NS_f = len(f_states)
print(f"Fortran: {NS_f} states after weeding")

# Generate all Python states (same as Fortran's pre-weeding set)
p = Params(); p.N=3; p.NF=1; p.B=1; p.K=K; p.rlamb=lam; p.iflv[0]=3
st_all = StateData(); perm = PermTables(); flav = FlavorTables()
se = cs1(p)
qcdsta(p, st_all, perm, flav)
n_all = st_all.numsta
print(f"Python generated: {n_all} states (should match Fortran pre-weeding)")

# Match each Fortran post-weeding state to a Python state index
matched_py_indices = []
unmatched = 0
for fi, fs in enumerate(f_states):
    found = False
    for s in range(n_all):
        loc = st_all.mstinf[s, 0] - 1
        lng = st_all.mstinf[s, 1]
        if lng != fs['lng']:
            continue
        py_types = [int(st_all.mstate[loc, pp]) for pp in range(lng)]
        py_moms = [int(st_all.mstate[loc+2, pp]) for pp in range(lng)]
        if py_types == fs['types'] and py_moms == fs['moms']:
            matched_py_indices.append(s)
            found = True
            break
    if not found:
        unmatched += 1
        if unmatched <= 5:
            print(f"  WARNING: Fortran state {fi+1} not found: types={fs['types']}, moms={fs['moms']}")

print(f"Matched: {len(matched_py_indices)}/{NS_f}, unmatched: {unmatched}")

if unmatched > 0:
    print("Cannot proceed with Fortran comparison - state matching failed")
    print("Plotting Python-only results")
    matched_py_indices = None

# Build norm and Hamiltonian for the Fortran's exact set of states
if matched_py_indices is not None and len(matched_py_indices) == NS_f:
    # Create a sub-mstinf for just the matched states (in Fortran order)
    sub_mstinf = np.zeros((NS_f, 8), dtype=int)
    for i, si in enumerate(matched_py_indices):
        sub_mstinf[i, :] = st_all.mstinf[si, :]

    # Build norm matrix for these states
    se2 = compute_selfen(3)
    _, _, Nm_sub = build_matrices(0, st_all.mstate, sub_mstinf, NS_f,
                                   3, 1, 1, K, se2, 1e-8, ncpus)
    Nm_sub = Nm_sub[:NS_f, :NS_f]

    # Orthonormalize
    w_n, z_n = eigh(Nm_sub)
    good = w_n > 0.5
    n_good = int(np.sum(good))
    print(f"Fortran-matched: {NS_f} states, {n_good} orthonormal (good norm eigenvalues)")
    Z_sub = z_n[:, good] / np.sqrt(w_n[good])[np.newaxis, :]

    # Build Hamiltonian
    h0_sub, hi_sub, _ = build_matrices(1, st_all.mstate, sub_mstinf, NS_f,
                                        3, 1, 1, K, se2, 1e-8, ncpus)

    # Transform to orthonormal basis
    hnu = Z_sub.T @ hi_sub[:NS_f, :NS_f] @ Z_sub
    hnu0 = np.diag(Z_sub.T @ h0_sub[0, :NS_f, :NS_f] @ Z_sub)

    rl2 = lam**2
    hnu_full = rl2 * hnu
    hnu_full[np.diag_indices(n_good)] += (1 - rl2) * hnu0

    w_f, z_f = eigh(hnu_full)
    w_mass_f = K * w_f / 2.0

    # Eigenvectors in original (post-weeding) basis
    c_f = Z_sub @ z_f

    # Also compare with Fortran's own eigenvalues from qcdf.ham
    with open(ham_file, 'r') as f:
        hlines = f.readlines()
    header = hlines[1].split()
    fNS_ham = int(header[5])
    f_hnu0_ham = np.zeros(fNS_ham)
    for ii in range(fNS_ham):
        f_hnu0_ham[ii] = float(hlines[3+ii].strip())
    f_hnu_ham = np.zeros((fNS_ham, fNS_ham))
    idx = 4 + fNS_ham
    while idx < len(hlines):
        parts = hlines[idx].split()
        if len(parts) < 3: idx += 1; continue
        for chunk in range(0, len(parts), 3):
            if chunk+2 >= len(parts): break
            try: i_f=int(parts[chunk]); j_f=int(parts[chunk+1]); v_f=float(parts[chunk+2])
            except: break
            if i_f==-1: break
            f_hnu_ham[i_f-1,j_f-1]=v_f; f_hnu_ham[j_f-1,i_f-1]=v_f
        idx += 1
    f_full_ham = rl2 * f_hnu_ham + np.diag((1-rl2)*f_hnu0_ham)
    f_eig_ham = np.sort(eigh(f_full_ham, eigvals_only=True))
    f_mass_ham = K / 2.0 * f_eig_ham

    print(f"\nEigenvalue comparison (M²_code):")
    print(f"{'#':>3s}  {'F(ham)':>12s}  {'F(rebuilt)':>12s}  {'Python':>12s}  {'ham-reb':>10s}  {'ham-py':>10s}")
    for i in range(min(5, n_good, len(w_p))):
        fm = f_mass_ham[i]; fr = w_mass_f[i]; pp = w_p[i]
        print(f"{i:3d}  {fm:12.4f}  {fr:12.4f}  {pp:12.4f}  {abs(fm-fr):10.4f}  {abs(fm-pp):10.4f}")

    # ═══════════════════════════════════════════════════════════
    # Structure functions from Fortran-matched basis
    # ═══════════════════════════════════════════════════════════
    K_paper = K / 2.0

    def fortran_sf(state_idx, nparton_filter=None):
        evec = c_f[:, state_idx]
        Nc = Nm_sub @ evec
        q_q = np.zeros(K+1)
        q_a = np.zeros(K+1)
        for s in range(NS_f):
            loc = sub_mstinf[s, 0] - 1
            lng = sub_mstinf[s, 1]
            if nparton_filter is not None and lng != nparton_filter:
                continue
            wt = evec[s] * Nc[s]
            for pp in range(lng):
                mom = st_all.mstate[loc+2, pp]
                if st_all.mstate[loc, pp] == 1:
                    q_q[mom] += wt
                else:
                    q_a[mom] += wt
        x = np.arange(1, K, 2) / float(K)
        return x, np.array([K_paper*q_q[n] for n in range(1,K,2)]), \
                  np.array([K_paper*q_a[n] for n in range(1,K,2)])

    has_fortran = True
else:
    has_fortran = False

# ═══════════════════════════════════════════════════════════════
# PLOT
# ═══════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 3, figsize=(15, 9))

for si in range(3):
    Mg_py = code_to_M_over_g(max(w_p[si], 0), lam)

    # Python
    x_p, qq3_p, _ = extract_structure_function(si, c_p, nso_p, mi_p, st_p.mstate, Nm_p, K, 3, 1, nparton_filter=3)
    x_p, qq5_p, qa5_p = extract_structure_function(si, c_p, nso_p, mi_p, st_p.mstate, Nm_p, K, 3, 1, nparton_filter=5)

    # Top: valence
    ax = axes[0, si]
    ax.plot(x_p, qq3_p, 'k-o', ms=3, lw=1, label='Python qqq')

    # Bottom: 5-parton
    ax2 = axes[1, si]
    sq = [1e3, 1e2, 1e2][si]
    ax2.plot(x_p, qq5_p*sq, 'r-o', ms=3, lw=1, label=f'Py q (x{sq:.0e})')
    ax2.plot(x_p, qa5_p*sq, 'b-o', ms=3, lw=1, label=f'Py qbar (x{sq:.0e})')

    if has_fortran:
        Mg_f = code_to_M_over_g(max(w_mass_f[si], 0), lam)
        x_f, qq3_f, _ = fortran_sf(si, nparton_filter=3)
        x_f, qq5_f, qa5_f = fortran_sf(si, nparton_filter=5)

        ax.plot(x_f, qq3_f, 'g--^', ms=3, lw=1, label='Fortran-matched qqq')
        ax2.plot(x_f, qq5_f*sq, 'g--^', ms=3, lw=1, label=f'F q (x{sq:.0e})')
        ax2.plot(x_f, qa5_f*sq, 'c--v', ms=3, lw=1, label=f'F qbar (x{sq:.0e})')
        ax.set_title(f"State {si+1}: Py={Mg_py:.2f}, F={Mg_f:.2f}")

        f_sr = np.sum(qq3_f + qq5_f) * 2/K
    else:
        ax.set_title(f"State {si+1}: Py={Mg_py:.2f}")
        f_sr = 0

    p_sr = np.sum(qq3_p + qq5_p) * 2/K
    print(f"State {si}: Py sr={p_sr:.4f}" + (f", F sr={f_sr:.4f}" if has_fortran else ""))

    ax.set_xlabel("x"); ax.set_ylabel("q(x)"); ax.set_xlim(0, 1); ax.legend(fontsize=7)
    ax2.set_xlabel("x"); ax2.set_ylabel("q(x)"); ax2.set_xlim(0, 1); ax2.legend(fontsize=6)
    ax2.set_title(f"5-parton state {si+1}")

title = f"K={K}: Python ({nso_p} st)"
if has_fortran:
    title += f" vs Fortran-matched ({NS_f} st, {n_good} orth)"
fig.suptitle(title, fontsize=13)
fig.tight_layout()
outname = f'fig6_K{K}_comparison.pdf'
fig.savefig(outname, dpi=150)
print(f"Saved {outname}")
