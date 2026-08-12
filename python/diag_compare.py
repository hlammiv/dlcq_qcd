#!/usr/bin/env python3
"""
Diagnostic: compare Python vs Fortran for N=3, B=1, K=7 baryon.

Usage:
  python diag_compare.py

For Fortran comparison, you need to patch qcdf.f:
  1. Find "CALL QCDSTA" (around line 255)
  2. Add this line BEFORE it:
       iflv(1) = N*B
  3. Recompile
  4. Run with:  echo "0 3 1 1 0.3325 -1.0 0 7" | tr ' ' '\n' | ./a.out
  5. The Fortran will write qcdf.ham and qcdf.out
  6. Copy qcdf.ham and qcdf.out to this directory
  7. Re-run this script with --compare to diff the matrices
"""

import os, sys
import numpy as np
from scipy.linalg import eigh

os.environ['QCDF_NCPUS'] = '1'

# ── Parameters ──
N, NF, B, K = 3, 1, 1, 21
lam = 0.3325

print(f"N={N}, NF={NF}, B={B}, K={K}, λ={lam}")
print("=" * 60)

# ── Build using qcdf_opt.py ──
from qcdf_opt import (build_matrices, compute_selfen, weed,
                       Params, StateData, PermTables, FlavorTables, qcdsta)

p = Params(); p.N=N; p.NF=NF; p.B=B; p.K=K; p.rlamb=lam; p.iflv[0]=N*B
st = StateData(); perm = PermTables(); flav = FlavorTables()
se = compute_selfen(N)
qcdsta(p, st, perm, flav)
n = st.numsta
print(f"\nStates generated: {n}")

# Show states
for i in range(n):
    loc = st.mstinf[i,0]-1; lng = st.mstinf[i,1]
    moms = list(st.mstate[loc+2,:lng])
    types = list(st.mstate[loc,:lng])
    print(f"  state {i}: types={types} moms={moms}")

# Build norm
_, _, Nm = build_matrices(0, st.mstate, st.mstinf, n, N, NF, B, K, se, 1e-8, 1)
Nm = Nm[:n,:n]
print(f"\nNORM MATRIX ({n}×{n}):")
np.set_printoptions(precision=10, linewidth=120)
for i in range(n):
    print(f"  [{', '.join(f'{Nm[i,j]:14.8f}' for j in range(n))}]")

# Norm eigenvalues
wn = np.linalg.eigvalsh(Nm)
print(f"\nNorm eigenvalues: {wn}")

# Weed
Nm_w, mi_w, n_w = weed(Nm.copy(), st.mstinf[:n].copy(), n)
print(f"\nAfter weeding: {n_w} states (from {n})")

# Orthonormalize
w_n, z_n = eigh(Nm_w[:n_w,:n_w])
good = w_n > 0.5; ng = int(np.sum(good))
Z = z_n[:, good] / np.sqrt(w_n[good])[np.newaxis, :]
print(f"Orthonormal states: {ng}")

# Build Hamiltonian
h0, hi, _ = build_matrices(1, st.mstate, mi_w, n_w, N, NF, B, K, se, 1e-8, 1)
h0 = h0[:,:n_w,:n_w]; hi = hi[:n_w,:n_w]

print(f"\nHAM0 (free) diagonal in original basis:")
for i in range(n_w):
    print(f"  state {i}: H0 = {h0[0,i,i]:.10f}")

print(f"\nHAM_I (interaction) in original basis ({n_w}×{n_w}):")
for i in range(n_w):
    print(f"  [{', '.join(f'{hi[i,j]:14.8f}' for j in range(n_w))}]")

# Transform to orthonormal basis
hnu = Z.T @ hi @ Z
hnu0_diag = np.diag(Z.T @ h0[0] @ Z)
print(f"\nHNU0 (free) diagonal in orthonormal basis:")
for i in range(ng):
    print(f"  state {i}: HNU0 = {hnu0_diag[i]:.10f}")

print(f"\nHNU (interaction) in orthonormal basis ({ng}×{ng}):")
for i in range(ng):
    print(f"  [{', '.join(f'{hnu[i,j]:14.8f}' for j in range(ng))}]")

# Combine and diagonalize
rl2 = lam**2
hnu_full = rl2 * hnu
hnu_full[np.diag_indices(ng)] += (1 - rl2) * hnu0_diag
w_eig, z_eig = eigh(hnu_full)
w_mass = K * w_eig / 2.0

print(f"\nFINAL EIGENVALUES (M²_code = K/2 * eigenvalue):")
for i in range(ng):
    print(f"  state {i}: eigenvalue = {w_eig[i]:.12f},  M²_code = {w_mass[i]:.12f}")

# ── Compare with Fortran if qcdf.ham exists ──
if '--compare' in sys.argv and os.path.exists('qcdf.ham'):
    print("\n" + "=" * 60)
    print("COMPARING WITH FORTRAN qcdf.ham")
    print("=" * 60)
    
    with open('qcdf.ham', 'r') as f:
        lines = f.readlines()
    
    # Parse header
    header = lines[1].split()
    f_N, f_NF, f_B, f_LPN, f_K, f_NUMSTA = [int(x) for x in header[:6]]
    f_CUTOFF = float(header[6])
    print(f"Fortran: N={f_N}, NF={f_NF}, B={f_B}, K={f_K}, NUMSTA={f_NUMSTA}")
    
    # Parse HNU0 diagonal
    f_hnu0 = []
    for i in range(f_NUMSTA):
        f_hnu0.append(float(lines[3 + i].strip()))
    
    print(f"\nHNU0 comparison (Fortran vs Python):")
    for i in range(min(ng, f_NUMSTA)):
        diff = abs(f_hnu0[i] - hnu0_diag[i])
        match = "✓" if diff < 1e-8 else f"✗ diff={diff:.2e}"
        print(f"  state {i}: Fortran={f_hnu0[i]:.10f}  Python={hnu0_diag[i]:.10f}  {match}")
    
    # Parse HNU non-zero elements
    f_hnu = np.zeros((f_NUMSTA, f_NUMSTA))
    idx = 4 + f_NUMSTA  # skip header + label + HNU0 entries
    while idx < len(lines):
        parts = lines[idx].split()
        if len(parts) < 3:
            idx += 1
            continue
        for chunk in range(0, len(parts), 3):
            if chunk + 2 >= len(parts):
                break
            i_f = int(parts[chunk])
            j_f = int(parts[chunk + 1])
            v_f = float(parts[chunk + 2])
            if i_f == -1:
                break
            f_hnu[i_f - 1, j_f - 1] = v_f
            f_hnu[j_f - 1, i_f - 1] = v_f
        idx += 1
    
    print(f"\nHNU comparison (Fortran vs Python):")
    max_diff = 0
    for i in range(min(ng, f_NUMSTA)):
        for j in range(i + 1):
            diff = abs(f_hnu[i, j] - hnu[i, j])
            if diff > max_diff:
                max_diff = diff
            if diff > 1e-6 or abs(f_hnu[i, j]) > 1e-10:
                match = "✓" if diff < 1e-8 else f"✗ diff={diff:.2e}"
                print(f"  HNU[{i},{j}]: Fortran={f_hnu[i,j]:14.8f}  Python={hnu[i,j]:14.8f}  {match}")
    print(f"  Max difference: {max_diff:.2e}")
    
    # Compare eigenvalues from qcdf.out
    if os.path.exists('qcdf.out'):
        print(f"\nEigenvalue comparison from qcdf.out:")
        with open('qcdf.out', 'r') as f:
            out_lines = f.readlines()
        # Find eigenvalue section
        for line in out_lines:
            if 'EIGENVAL' in line.upper():
                print(f"  {line.strip()}")
else:
    if '--compare' in sys.argv:
        print("\nNo qcdf.ham found. Run the Fortran first:")
        print('  echo "0 3 1 1 0.3325 -1.0 0 7" | tr \' \' \'\\n\' | ./a.out')
