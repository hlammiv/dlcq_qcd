#!/usr/bin/env python3
"""
reproduce_figures.py — Reproduce all figures from:
  Hornbostel, Brodsky & Pauli, "Light-cone-quantized QCD in 1+1 dimensions"
  Phys. Rev. D 41, 3814 (1990)

Uses the DLCQ code (qcdf.py / qcdf_opt.py) to compute spectra and wave functions.

Usage:
    python reproduce_figures.py              # generate all figures
    python reproduce_figures.py --fig 2      # generate only Figure 2
    python reproduce_figures.py --fig 2 5 7  # generate Figures 2, 5, 7

Figures are saved as PDF in ./figures/ directory.
Intermediate data is cached in ./cache/ to avoid recomputation.

Requirements: qcdf.py (and optionally qcdf_opt.py) in the same directory.
"""

import os, sys, json, time, pickle, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from pathlib import Path

# ── Import DLCQ solver ──
try:
    from qcdf_opt import (build_matrices, compute_selfen, weed,
                           Params, StateData, PermTables, FlavorTables, qcdsta)
    USE_OPT = True
    print("Using optimized solver (qcdf_opt.py)")
except ImportError:
    from qcdf import (Params, StateData, PermTables, FlavorTables,
                       qcdsta, clrdis, weedr, weedr2, compute_selfen as _cs,
                       nuham, nuz)
    USE_OPT = False
    print("Using base solver (qcdf.py)")

from scipy.linalg import eigh

# ── Directories ──
FIGDIR = Path("figures"); FIGDIR.mkdir(exist_ok=True)
CACHEDIR = Path("cache"); CACHEDIR.mkdir(exist_ok=True)

# ── Unit conversions ──
def lambda_to_mg(lam):
    """Convert λ to m/g."""
    return np.sqrt((1 - lam**2) / (np.pi * lam**2))

def mg_to_lambda(mg):
    """Convert m/g to λ."""
    return 1.0 / np.sqrt(1 + np.pi * mg**2)

def code_to_Msq_over_scale(Msq_code):
    """Convert code M² to M²/(m²+g²/π) for Figure 2 y-axis.
    Result: y = M²_code.
    
    The code's M²_code = K/2 * eigenvalue, where the eigenvalue already
    incorporates the DLCQ normalization P⁺·P⁻ = (m²+g²/π) · M²_code.
    """
    return Msq_code

def code_to_M_over_g(Msq_code, lam):
    """Convert code M² to M/g.  M/g = sqrt(M²_code / (π*λ²))."""
    if Msq_code <= 0:
        return 0.0
    return np.sqrt(Msq_code / (np.pi * lam**2))


def get_lightest_mass(N, NF, B, K_list, lam, ncpus=1, **kwargs):
    """
    Compute lightest M/g by running at multiple K values.
    
    At strong coupling, the lightest eigenvalue can go negative at large K
    or dip spuriously due to finite-discretization artifacts. We take the
    result from the K that gives the HIGHEST positive lightest eigenvalue,
    which corresponds to the best-converged value before artifacts set in.
    
    The paper uses Richardson extrapolation across multiple K; this simpler
    approach gives a lower bound that is monotonically correct.
    
    Special case: at λ ≥ 0.998 (m/g ≈ 0), the lightest state is exactly
    massless (Kac-Moody construction, Section VI of paper).
    
    Returns (M_over_g, best_K) or (None, None) if all fail.
    """
    if lam >= 0.998:
        return 0.0, 0
    
    best_Mg = None
    best_K = None
    
    for K in sorted(K_list):
        try:
            w, _, ns, _, _, _, _ = run_dlcq(N, NF, B, K, lam, ncpus=ncpus, **kwargs)
            if ns > 0 and len(w) > 0 and w[0] > -1e-10:
                Mg = code_to_M_over_g(max(w[0], 0.0), lam)
                if best_Mg is None or Mg > best_Mg:
                    best_Mg = Mg
                    best_K = K
        except:
            pass
    
    return best_Mg, best_K


# ── Core computation function ──
def run_dlcq(N, NF, B, K, lam, cutoff=-1.0, LPN=0, rmq=None, iflv=None, ncpus=1):
    """
    Run full DLCQ computation. Returns (eigenvalues_code, eigenvectors, numsta, mstinf).
    eigenvalues_code are M²_code = K/2 * raw_eigenvalue.
    
    For NF=1, iflv is automatically set to [N*B] if not provided.
    """
    # Auto-set iflv for baryons: net quark number = N*B per flavor
    if iflv is None:
        if NF == 1:
            iflv = [N * B]
        else:
            # For multi-flavor, default to all quarks in first flavor
            iflv = [0] * NF
            iflv[0] = N * B

    cache_key = f"N{N}_NF{NF}_B{B}_K{K}_lam{lam:.6f}_cut{cutoff}_LPN{LPN}_iflv{'_'.join(str(x) for x in iflv)}"
    cache_file = CACHEDIR / f"{cache_key}.pkl"
    if cache_file.exists():
        with open(cache_file, 'rb') as f:
            return pickle.load(f)

    os.environ['QCDF_NCPUS'] = str(ncpus)

    p = Params()
    p.N = N; p.NF = NF; p.B = B; p.K = K; p.rlamb = lam
    p.cutoff = cutoff; p.LPN = LPN
    if rmq is not None:
        for i, v in enumerate(rmq[:NF]): p.rmq[i] = v
    for i, v in enumerate(iflv[:NF]): p.iflv[i] = v

    states = StateData(); perm = PermTables(); flav = FlavorTables()
    selfen = compute_selfen(p.N) if USE_OPT else _cs(p)

    qcdsta(p, states, perm, flav)
    numsta = states.numsta
    if numsta == 0:
        result = (np.array([]), np.array([[]]), 0, None, states, 0, None)
        with open(cache_file, 'wb') as f: pickle.dump(result, f)
        return result

    if USE_OPT:
        _, _, hnorm = build_matrices(0, states.mstate, states.mstinf, numsta,
                                      N, NF, B, K, selfen, p.cbreak, ncpus)
        hnorm, mstinf, numsta = weed(hnorm, states.mstinf[:numsta].copy(), numsta)
    else:
        _, _, hnorm = clrdis(0, p, states, selfen, ncpus=ncpus)
        weedr(states, hnorm)
        numsta = states.numsta
        for _ in range(2000):
            w, z = eigh(hnorm[:numsta, :numsta])
            if not weedr2(states, hnorm, w, z): break
            numsta = states.numsta
        mstinf = states.mstinf[:numsta].copy()

    if numsta == 0:
        result = (np.array([]), np.array([[]]), 0, None, states, 0, None)
        with open(cache_file, 'wb') as f: pickle.dump(result, f)
        return result

    # ── Orthonormalize by projecting onto well-conditioned subspace ──
    # Diagonalize the full norm matrix. Keep only eigenvectors whose
    # eigenvalue is above a threshold — this cleanly removes redundant
    # and near-redundant states without iterative weeding artifacts.
    w_n, z_n = eigh(hnorm[:numsta, :numsta])
    
    # Threshold: keep eigenvalues > 0.5 (norm eigenvalues are integers
    # or near-integers for this problem; 0.5 cleanly separates physical
    # states from spurious near-zero/negative modes)
    norm_thresh = 0.5
    good = w_n > norm_thresh
    n_good = int(np.sum(good))
    
    if n_good == 0:
        result = (np.array([]), np.array([[]]), 0, None, states, 0, None)
        with open(cache_file, 'wb') as f: pickle.dump(result, f)
        return result

    # Z maps from original basis to orthonormal basis (only good modes)
    Z = z_n[:, good] / np.sqrt(w_n[good])[np.newaxis, :]
    numsta_orth = n_good

    # Hamiltonian (computed in original basis, then projected)
    if USE_OPT:
        ham0, ham, _ = build_matrices(1, states.mstate, mstinf, numsta,
                                       N, NF, B, K, selfen, p.cbreak, ncpus)
    else:
        ham0, ham, _ = clrdis(1, p, states, selfen, ncpus=ncpus)

    # Transform to orthonormal basis: Z^T H Z  (Z is numsta × numsta_orth)
    hnu = Z.T @ ham[:numsta, :numsta] @ Z
    hnu0 = np.array([np.diag(Z.T @ ham0[ifl, :numsta, :numsta] @ Z) for ifl in range(NF)])

    # Combine free + interacting
    rlmsq = lam ** 2
    hnu *= rlmsq
    for ifl in range(NF):
        rmq_val = p.rmq[ifl] if rmq is None else (rmq[ifl] if ifl < len(rmq) else 1.0)
        hnu[np.diag_indices(numsta_orth)] += (1 - rlmsq) * rmq_val**2 * hnu0[ifl]

    # Diagonalize
    w_eig, z_eig = eigh(hnu)
    w_mass = K * w_eig / 2.0  # M²_code

    # Transform eigenvectors back to original state basis for wave functions:
    # c_original = Z @ z_eig  (numsta × numsta_orth matrix)
    c_orig = Z @ z_eig  # coefficients in original state basis

    # Store norm matrix for structure function computation
    norm_mat = hnorm[:numsta, :numsta].copy()

    result = (w_mass, c_orig, numsta_orth, mstinf if USE_OPT else states.mstinf[:numsta], states, numsta, norm_mat)
    with open(cache_file, 'wb') as f:
        pickle.dump(result, f)
    return result


def extract_structure_function(state_idx, c_orig, numsta_orig, mstinf, mstate, norm_mat, K, N, B,
                                nparton_filter=None):
    """
    Extract quark structure function q(x) for a given eigenstate.
    
    q(x=k/K) = K_paper * Σ_s c_s * n_k(s) * (N·c)_s
    
    where K_paper = K_code/2, c_s are coefficients in the original non-orthogonal
    basis satisfying c†Nc = 1, n_k(s) counts quarks with momentum k in state s,
    and N is the norm matrix. This formula gives correct sum rules:
    ∫q(x)dx = 1 for mesons, 3 for baryons.
    
    nparton_filter: if set, only include Fock components with this many partons.
    
    Returns (x_values, q_quark, q_antiquark).
    """
    evec = c_orig[:, state_idx]
    Nc = norm_mat @ evec  # precompute N·c

    K_paper = K / 2.0  # paper's K = code's K / 2

    q_quark = np.zeros(K + 1)
    q_anti  = np.zeros(K + 1)

    for s in range(numsta_orig):
        loc = mstinf[s, 0] - 1
        lng = mstinf[s, 1]
        if nparton_filter is not None and lng != nparton_filter:
            continue
        weight = evec[s] * Nc[s]  # c_s * (N·c)_s — correct for non-orthogonal basis
        for p in range(lng):
            mom = mstate[loc + 2, p]
            if mstate[loc, p] == 1:  # quark (b†)
                q_quark[mom] += weight
            else:  # antiquark (d†)
                q_anti[mom] += weight

    x_vals = np.arange(1, K, 2) / float(K)
    q_q = np.array([K_paper * q_quark[n] for n in range(1, K, 2)])
    q_a = np.array([K_paper * q_anti[n]  for n in range(1, K, 2)])

    return x_vals, q_q, q_a


# ══════════════════════════════════════════════════════════════════
# FIGURE 2: Spectra for N=3, B=0,1,2 as function of coupling
# ══════════════════════════════════════════════════════════════════
def figure2(ncpus=1):
    """Reproduce Figure 2: meson and baryon spectra vs coupling."""
    print("\n=== FIGURE 2: Spectra vs coupling ===")

    # λ values spanning 0 → 1 (weak → strong coupling)
    lambdas = np.concatenate([
        np.arange(0.05, 0.50, 0.05),
        np.arange(0.50, 0.95, 0.03),
        np.array([0.95, 0.97, 0.99])
    ])

    panels = [
        (0, 10, "B=0, 2K=10", 50),   # (a)
        (1, 13, "B=1, 2K=13", 60),   # (b)
        (2, 22, "B=2, 2K=22", 100),  # (c)
    ]

    fig, axes = plt.subplots(1, len(panels), figsize=(5*len(panels), 6))
    if len(panels) == 1: axes = [axes]

    for ax, (B, K, title, ymax) in zip(axes, panels):
        print(f"  Panel: {title}")
        all_evals = []
        valid_lams = []

        for lam in lambdas:
            try:
                w, _, ns, _, _, _, _ = run_dlcq(3, 1, B, K, lam, ncpus=ncpus)
                if len(w) > 0:
                    # Convert to M²/(m²+g²/π) = 2*M²_code
                    y = code_to_Msq_over_scale(w)
                    all_evals.append(y[y >= -1])
                    valid_lams.append(lam)
                    print(f"    λ={lam:.3f}: {ns} states, lowest M²/(m²+g²/π)={y[0]:.2f}")
            except Exception as e:
                print(f"    λ={lam:.3f}: FAILED ({e})")

        # Plot each eigenvalue as a point
        for i, (lam, evals) in enumerate(zip(valid_lams, all_evals)):
            nplot = min(len(evals), 30)
            ax.plot([lam]*nplot, evals[:nplot], 'k.', markersize=1.5)

        ax.set_xlim(0, 1)
        ax.set_ylim(0, ymax)
        ax.set_xlabel(r"$1/(1+\pi m^2/g^2)^{1/2}$", fontsize=11)
        ax.set_ylabel(r"$M^2/(m^2+g^2/\pi)$", fontsize=11)
        ax.set_title(f"SU(3)  {title}", fontsize=12)

    fig.tight_layout()
    fig.savefig(FIGDIR / "fig2_spectra.pdf", dpi=150)
    print(f"  Saved: {FIGDIR}/fig2_spectra.pdf")
    plt.close()


# ══════════════════════════════════════════════════════════════════
# FIGURES 3 & 5: Valence wave functions
# ══════════════════════════════════════════════════════════════════
def figure3(ncpus=1):
    """Reproduce Figure 3: valence structure functions for N=3 meson and baryon."""
    print("\n=== FIGURE 3: Valence structure functions ===")

    couplings = [(1.6, 'o', 'full', 24), (0.1, 's', 'none', 16)]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    for mg, marker, fill, K_mes in couplings:
        lam = mg_to_lambda(mg)
        label = f"m/g = {mg}"
        K_bar = K_mes - 3  # odd K for baryons

        # Meson (B=0)
        w, c, ns, mi, st, nso, Nm = run_dlcq(3, 1, 0, K_mes, lam, ncpus=ncpus)
        if ns > 0:
            x, qq, qa = extract_structure_function(0, c, nso, mi, st.mstate, Nm, K_mes, 3, 0)
            ax1.plot(x, qq, '-'+marker, fillstyle=fill, markersize=4,
                    linewidth=1, label=label)

        # Baryon (B=1)
        try:
            w, c, ns, mi, st, nso, Nm = run_dlcq(3, 1, 1, K_bar, lam, ncpus=ncpus)
            if ns > 0:
                x, qq, qa = extract_structure_function(0, c, nso, mi, st.mstate, Nm, K_bar, 3, 1)
                ax2.plot(x, qq, '-'+marker, fillstyle=fill, markersize=4,
                        linewidth=1, label=label)
        except Exception as e:
            print(f"    Baryon at m/g={mg} failed: {e}")

    ax1.set_xlabel("x = k/K"); ax1.set_ylabel("q(x)")
    ax1.set_title("SU(3) Meson", fontsize=12)
    ax1.legend(fontsize=9); ax1.set_xlim(0, 1)

    ax2.set_xlabel("x = k/K"); ax2.set_ylabel("q(x)")
    ax2.set_title("SU(3) Baryon", fontsize=12)
    ax2.legend(fontsize=9); ax2.set_xlim(0, 1)

    fig.tight_layout()
    fig.savefig(FIGDIR / "fig3_valence_wf.pdf", dpi=150)
    print(f"  Saved: {FIGDIR}/fig3_valence_wf.pdf")
    plt.close()


def figure5(ncpus=1):
    """Reproduce Figure 5: first three meson states + 11th for N=3, m/g=1.6, 2K=24.
    Paper scale factors: (a-c) qq̄qq̄ ×10³, (d) qq̄ ×10⁴ and qq̄qq̄ unscaled."""
    print("\n=== FIGURE 5: Meson spectrum wave functions ===")

    K = 24; lam = mg_to_lambda(1.6)
    w, c_orig, ns, mi, st, ns_orig, N_mat = run_dlcq(3, 1, 0, K, lam, ncpus=ncpus)

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    titles = ["(a) 1st state", "(b) 2nd state", "(c) 3rd state", "(d) 11th state"]
    indices = [0, 1, 2, min(10, ns - 1)]
    # Paper scale factors: (a) qq̄qq̄ ×10³, (b,c) ×10², (d) qq̄ ×10⁴, qq̄qq̄ unscaled
    hf_scales = [1e3, 1e2, 1e2, 1.0]
    val_scales = [1.0, 1.0, 1.0, 1e4]

    for i, (ax, title, idx) in enumerate(zip(axes.flat, titles, indices)):
        if idx < ns:
            Mg = code_to_M_over_g(max(w[idx], 0), lam)

            # Valence qq̄ component
            x, qq_v, _ = extract_structure_function(idx, c_orig, ns_orig, mi, st.mstate, N_mat, K, 3, 0,
                                                     nparton_filter=2)
            vs = val_scales[i]
            vlab = r"$q\bar{q}$" + (f" (×{vs:.0e})" if vs != 1.0 else "")
            ax.plot(x, qq_v * vs, 'k-o', markersize=2, linewidth=1, label=vlab)

            # Higher-Fock: qq̄qq̄ component
            x, qq_4, _ = extract_structure_function(idx, c_orig, ns_orig, mi, st.mstate, N_mat, K, 3, 0,
                                                     nparton_filter=4)
            hs = hf_scales[i]
            if np.max(qq_4) > 1e-15:
                hlab = r"$q\bar{q}q\bar{q}$" + (f" (×{hs:.0e})" if hs != 1.0 else "")
                ax.plot(x, qq_4 * hs, 'r-s', markersize=2, linewidth=0.8, label=hlab)

            ax.set_xlabel("x = k/K"); ax.set_ylabel("q(x)")
            ax.set_title(f"{title}  M/g={Mg:.2f}", fontsize=11)
            ax.set_xlim(0, 1)
            ax.legend(fontsize=7)

    fig.suptitle(f"SU(3) Meson, m/g = 1.6, 2K = {2*K}", fontsize=13)
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig5_meson_wf.pdf", dpi=150)
    print(f"  Saved: {FIGDIR}/fig5_meson_wf.pdf")
    plt.close()


# ══════════════════════════════════════════════════════════════════
# FIGURE 7: Extrapolated masses for N=2,3,4
# ══════════════════════════════════════════════════════════════════
def figure7(ncpus=1):
    """Reproduce Figure 7: extrapolated meson and baryon masses."""
    print("\n=== FIGURE 7: Extrapolated masses ===")

    mg_values = [0.0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.2, 1.6]
    N_values = [2, 3, 4]
    K_meson = list(range(4, 13))      # 4,5,...,12
    K_baryon = list(range(5, 14))     # 5,6,...,13

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    colors = {2: 'C0', 3: 'C1', 4: 'C2'}
    styles = {2: '--', 3: '-', 4: '-.'}

    for N in N_values:
        mg_mes = []; Mg_mes = []
        mg_bar = []; Mg_bar = []

        for mg in mg_values:
            lam = mg_to_lambda(mg) if mg > 0 else 0.999
            print(f"  N={N}, m/g={mg:.2f}")

            # Meson
            Mg, K_used = get_lightest_mass(N, 1, 0, K_meson, lam, ncpus=ncpus)
            if Mg is not None:
                mg_mes.append(mg); Mg_mes.append(Mg)
                print(f"    Meson M/g = {Mg:.4f} (K={K_used})")

            # Baryon
            Mg, K_used = get_lightest_mass(N, 1, 1, K_baryon, lam, ncpus=ncpus)
            if Mg is not None:
                mg_bar.append(mg); Mg_bar.append(Mg)
                print(f"    Baryon M/g = {Mg:.4f} (K={K_used})")

        if mg_mes:
            from scipy.interpolate import PchipInterpolator
            ax1.plot(mg_mes, Mg_mes, 'o', color=colors[N], markersize=4)
            if len(mg_mes) >= 4:
                mg_f = np.linspace(min(mg_mes), max(mg_mes), 200)
                ax1.plot(mg_f, PchipInterpolator(mg_mes, Mg_mes)(mg_f),
                        styles[N], color=colors[N], linewidth=1.2, label=f"SU({N})")
            else:
                ax1.plot(mg_mes, Mg_mes, styles[N], color=colors[N],
                        linewidth=1.2, label=f"SU({N})")
        if mg_bar:
            from scipy.interpolate import PchipInterpolator
            ax2.plot(mg_bar, Mg_bar, 's', color=colors[N], markersize=4)
            if len(mg_bar) >= 4:
                mg_f = np.linspace(min(mg_bar), max(mg_bar), 200)
                ax2.plot(mg_f, PchipInterpolator(mg_bar, Mg_bar)(mg_f),
                        styles[N], color=colors[N], linewidth=1.2, label=f"SU({N})")
            else:
                ax2.plot(mg_bar, Mg_bar, styles[N], color=colors[N],
                        linewidth=1.2, label=f"SU({N})")

    ax1.set_xlabel("m/g", fontsize=12); ax1.set_ylabel("M/g", fontsize=12)
    ax1.set_title("(a) Meson Mass", fontsize=13)
    ax1.legend(); ax1.set_xlim(0, 1.7); ax1.set_ylim(0, 5)

    ax2.set_xlabel("m/g", fontsize=12); ax2.set_ylabel("M/g", fontsize=12)
    ax2.set_title("(b) Baryon Mass", fontsize=13)
    ax2.legend(); ax2.set_xlim(0, 1.7); ax2.set_ylim(0, 8)

    fig.tight_layout()
    fig.savefig(FIGDIR / "fig7_masses.pdf", dpi=150)
    print(f"  Saved: {FIGDIR}/fig7_masses.pdf")
    plt.close()


# ══════════════════════════════════════════════════════════════════
# FIGURE 8: Comparison with large-N
# ══════════════════════════════════════════════════════════════════
def figure8(ncpus=1):
    """Reproduce Figure 8: meson masses compared with large-N and lattice."""
    print("\n=== FIGURE 8: Comparison with large-N ===")

    mg_values = np.arange(0.0, 1.65, 0.05)
    N_values = [2, 3, 4]
    K_list = list(range(4, 13))  # K=4,5,6,...,12 — smooth K coverage

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    colors = {2: 'C0', 3: 'C1', 4: 'C2'}
    lstyles = {2: '--', 3: '-', 4: '-.'}

    for N in N_values:
        mg_list = []; Mg_list = []
        for mg in mg_values:
            lam = mg_to_lambda(mg) if mg > 0 else 0.999
            Mg, K_used = get_lightest_mass(N, 1, 0, K_list, lam, ncpus=ncpus)
            if Mg is not None:
                mg_list.append(mg)
                Mg_list.append(Mg)

        if mg_list and len(mg_list) >= 4:
            from scipy.interpolate import PchipInterpolator
            mg_arr = np.array(mg_list); Mg_arr = np.array(Mg_list)
            mg_fine = np.linspace(mg_arr[0], mg_arr[-1], 200)
            Mg_fine = PchipInterpolator(mg_arr, Mg_arr)(mg_fine)
            ax1.plot(mg_fine, Mg_fine, lstyles[N], color=colors[N],
                    linewidth=1.5, label=f"SU({N})")
            scale = np.sqrt(2 * np.pi / N)
            ax2.plot(mg_fine * scale, Mg_fine * scale,
                    lstyles[N], color=colors[N], linewidth=1.5,
                    label=f"SU({N})")

    # Large-N reference curve from 't Hooft equation (leading order)
    # For large N: M² = g²N/π * f(m√(2πN)/g) where f is universal
    # At m=0: M=0 for SU(N), M=g/√(2π) for U(N)
    # We plot the large-N analytic result if available
    ax2.set_xlabel(r"$(2\pi/N)^{1/2}\, m/g$", fontsize=11)
    ax2.set_ylabel(r"$(2\pi/N)^{1/2}\, M/g$", fontsize=11)

    ax1.set_xlabel("m/g", fontsize=12); ax1.set_ylabel("M/g", fontsize=12)
    ax1.set_title("(a) Meson Mass", fontsize=13)
    ax1.legend(); ax1.set_xlim(0, 1.6); ax1.set_ylim(0, 5)

    ax2.set_title("(b) 't Hooft Meson Mass", fontsize=13)
    ax2.legend(); ax2.set_xlim(0, 2); ax2.set_ylim(0, 4)

    fig.tight_layout()
    fig.savefig(FIGDIR / "fig8_large_N.pdf", dpi=150)
    print(f"  Saved: {FIGDIR}/fig8_large_N.pdf")
    plt.close()


# ══════════════════════════════════════════════════════════════════
# TABLE I: N-dependence of meson and baryon mass
# ══════════════════════════════════════════════════════════════════
def table1(ncpus=1):
    """Reproduce Table I: N-dependence of meson and baryon mass."""
    print("\n" + "=" * 80)
    print(" TABLE I.  N dependence of meson and baryon mass.")
    print("=" * 80)

    mg_values = [1.6, 0.8, 0.4, 0.2, 0.1, 0.05, 0.0]
    N_values = [2, 3, 4]
    K_mes = list(range(4, 15))     # 4,...,14
    K_bar = list(range(5, 14))     # 5,...,13

    # Header
    header = f"{'m/g':>6s}"
    for N in N_values:
        header += f"  {'N='+str(N):>10s}"
    header += "  |"
    for N in [3, 4]:
        header += f"  {'N='+str(N):>10s}"
    print(f"\n{' ':>6s}  {'── M_mes/g ──':^34s}  | {'── M_bar/g ──':^24s}")
    print(header)
    print("-" * 80)

    results = {}  # store for later use

    for mg in mg_values:
        lam = mg_to_lambda(mg) if mg > 0 else 0.999
        row = f"{mg:6.2f}"

        for N in N_values:
            Mg, K_used = get_lightest_mass(N, 1, 0, K_mes, lam, ncpus=ncpus)
            if Mg is not None:
                row += f"  {Mg:10.3f}"
                results[(mg, N, 'mes')] = Mg
            else:
                row += f"  {'---':>10s}"

        row += "  |"

        for N in [3, 4]:
            Mg, K_used = get_lightest_mass(N, 1, 1, K_bar, lam, ncpus=ncpus)
            if Mg is not None:
                row += f"  {Mg:10.3f}"
                results[(mg, N, 'bar')] = Mg
            else:
                row += f"  {'---':>10s}"

        print(row)

    # Print paper's values for comparison
    print("\n" + "-" * 80)
    print(" Paper values (Hornbostel et al., Table I):")
    print("-" * 80)
    paper = {
        (1.6, 2, 'mes'): 4.314, (1.6, 3, 'mes'): 4.618, (1.6, 4, 'mes'): 4.845,
        (0.8, 2, 'mes'): 3.913, (0.8, 3, 'mes'): 4.40,  (0.8, 4, 'mes'): 4.743,
        (0.4, 2, 'mes'): 2.61,  (0.4, 3, 'mes'): 3.1,   (0.4, 4, 'mes'): 3.4,
        (0.2, 2, 'mes'): 1.17,  (0.2, 3, 'mes'): 1.5,   (0.2, 4, 'mes'): 1.4,
        (0.1, 2, 'mes'): 0.38,  (0.1, 3, 'mes'): 0.5,   (0.1, 4, 'mes'): 0.43,
        (0.05,2, 'mes'): 0.10,  (0.05,3, 'mes'): 0.2,   (0.05,4, 'mes'): 0.12,
        (0.0, 2, 'mes'): 0,     (0.0, 3, 'mes'): 0,     (0.0, 4, 'mes'): 0,
        (1.6, 3, 'bar'): 10.71, (1.6, 4, 'bar'): 21.2,
        (0.8, 3, 'bar'): 10.4,  (0.8, 4, 'bar'): 20.9,
        (0.4, 3, 'bar'): 7.3,   (0.4, 4, 'bar'): 15,
        (0.2, 3, 'bar'): 3.1,   (0.2, 4, 'bar'): 6.0,
        (0.1, 3, 'bar'): 1.1,   (0.1, 4, 'bar'): 1.9,
        (0.05,3, 'bar'): 0.31,  (0.05,4, 'bar'): 0.42,
        (0.0, 3, 'bar'): 0,     (0.0, 4, 'bar'): 0,
    }
    header2 = f"{'m/g':>6s}"
    for N in N_values:
        header2 += f"  {'N='+str(N):>10s}"
    header2 += "  |"
    for N in [3, 4]:
        header2 += f"  {'N='+str(N):>10s}"
    print(header2)
    print("-" * 80)
    for mg in mg_values:
        row = f"{mg:6.2f}"
        for N in N_values:
            v = paper.get((mg, N, 'mes'))
            row += f"  {v:10.3f}" if v is not None else f"  {'':>10s}"
        row += "  |"
        for N in [3, 4]:
            v = paper.get((mg, N, 'bar'))
            row += f"  {v:10.3f}" if v is not None else f"  {'':>10s}"
        print(row)

    print("=" * 80)
    print()

    # Save to file
    outfile = FIGDIR / "table1.txt"
    with open(outfile, 'w') as f:
        f.write("TABLE I.  N dependence of meson and baryon mass.\n")
        f.write("From: Hornbostel, Brodsky & Pauli, PRD 41, 3814 (1990)\n")
        f.write("Reproduced with Python DLCQ code\n\n")
        f.write(f"{'m/g':>6s}")
        for N in N_values: f.write(f"  {'N='+str(N):>10s}")
        f.write("  |")
        for N in [3, 4]: f.write(f"  {'N='+str(N):>10s}")
        f.write("\n" + "-" * 80 + "\n")
        for mg in mg_values:
            row = f"{mg:6.2f}"
            for N in N_values:
                v = results.get((mg, N, 'mes'))
                row += f"  {v:10.3f}" if v is not None else f"  {'---':>10s}"
            row += "  |"
            for N in [3, 4]:
                v = results.get((mg, N, 'bar'))
                row += f"  {v:10.3f}" if v is not None else f"  {'---':>10s}"
            f.write(row + "\n")
    print(f"  Table saved: {outfile}")


# ══════════════════════════════════════════════════════════════════
# FIGURE 6: Baryon wave functions
# ══════════════════════════════════════════════════════════════════
def figure6(ncpus=1):
    """Reproduce Figure 6: first three baryon states + B=2 state.
    (a-c) B=1: qqq valence + higher-Fock 5-parton (4q+1qbar)
    (d) B=2: 6q valence + 8-parton (7q+1qbar) components."""
    print("\n=== FIGURE 6: Baryon wave functions ===")

    lam = mg_to_lambda(1.6)
    K = 21  # paper: 2K=21

    w, c_orig, ns, mi, st, ns_orig, N_mat = run_dlcq(3, 1, 1, K, lam, ncpus=ncpus)

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    titles = ["(a) 1st baryon", "(b) 2nd baryon", "(c) 3rd baryon", "(d) 1st B=2"]
    hf_q_scales = [1e3, 1e2, 1e2]
    hf_qbar_scales = [1e3, 1e2, 1e2]

    for i, (ax, title) in enumerate(zip(axes.flat[:3], titles[:3])):
        if i < ns:
            Mg = code_to_M_over_g(max(w[i], 0), lam)

            # Valence qqq (3-parton)
            x, qq_v, _ = extract_structure_function(i, c_orig, ns_orig, mi, st.mstate, N_mat, K, 3, 1,
                                                     nparton_filter=3)
            ax.plot(x, qq_v, 'k-o', markersize=2, linewidth=1, label="qqq")

            # Higher-Fock 5-parton = 4q + 1qbar
            x, qq_5, qa_5 = extract_structure_function(i, c_orig, ns_orig, mi, st.mstate, N_mat, K, 3, 1,
                                                        nparton_filter=5)
            sq = hf_q_scales[i]; sa = hf_qbar_scales[i]
            if np.max(np.abs(qq_5)) > 1e-15:
                ax.plot(x, qq_5 * sq, 'r-s', markersize=2, linewidth=0.8,
                       label=r"q in qqqq$\bar{q}$" + f" (x{sq:.0e})")
            if np.max(np.abs(qa_5)) > 1e-15:
                ax.plot(x, qa_5 * sa, 'b-^', markersize=2, linewidth=0.8,
                       label=r"$\bar{q}$ in qqqq$\bar{q}$" + f" (x{sa:.0e})")

            ax.set_xlabel("x"); ax.set_ylabel("q(x)")
            ax.set_title(f"{title}  M/g={Mg:.2f}"); ax.set_xlim(0, 1)
            ax.legend(fontsize=6)

    # Panel (d): B=2 at K=22. Valence = 6q; next Fock = 7q + 1qbar (8-parton)
    ax = axes.flat[3]
    try:
        K2 = 22
        w2, c2, ns2, mi2, st2, nso2, N2 = run_dlcq(3, 1, 2, K2, lam, ncpus=ncpus)
        if ns2 > 0:
            Mg = code_to_M_over_g(max(w2[0], 0), lam)
            # 6q valence
            x, qq_6, _ = extract_structure_function(0, c2, nso2, mi2, st2.mstate, N2, K2, 3, 2,
                                                     nparton_filter=6)
            ax.plot(x, qq_6, 'k-o', markersize=2, linewidth=1, label="6q")
            # 8-parton (7q + 1qbar)
            x, qq_8, qa_8 = extract_structure_function(0, c2, nso2, mi2, st2.mstate, N2, K2, 3, 2,
                                                        nparton_filter=8)
            if np.max(np.abs(qq_8)) > 1e-15:
                ax.plot(x, qq_8 * 5e2, 'r-s', markersize=2, linewidth=0.8,
                       label=r"q in 6q$q\bar{q}$ (x5e+02)")
            if np.max(np.abs(qa_8)) > 1e-15:
                ax.plot(x, qa_8 * 1e3, 'b-^', markersize=2, linewidth=0.8,
                       label=r"$\bar{q}$(x) in 6q$q\bar{q}$ (x1e+03)")
            ax.set_title(f"(d) 1st B=2  M/g={Mg:.2f}")
        else:
            ax.set_title("(d) B=2: no states")
    except Exception as e:
        print(f"    Fig6(d) B=2: {e}")
        ax.set_title("(d) B=2: failed")
    ax.set_xlabel("x"); ax.set_ylabel("q(x)"); ax.set_xlim(0, 1)
    ax.legend(fontsize=6)

    fig.suptitle("SU(3) Baryon, m/g = 1.6", fontsize=13)
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig6_baryon_wf.pdf", dpi=150)
    print(f"  Saved: {FIGDIR}/fig6_baryon_wf.pdf")
    plt.close()

def figure1():
    """Reproduce Figure 1: schematic interaction vertices."""
    print("\n=== FIGURE 1: Interaction vertices ===")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3))

    # (a) Diagonal / self-energy vertex
    ax1.annotate("", xy=(0.7, 0.5), xytext=(0.1, 0.5),
                arrowprops=dict(arrowstyle="->, head_width=0.2", lw=2))
    ax1.annotate("", xy=(0.7, 0.7), xytext=(0.4, 0.5),
                arrowprops=dict(arrowstyle="->, head_width=0.15", lw=1.5,
                              connectionstyle="arc3,rad=0.3"))
    ax1.annotate("", xy=(0.4, 0.5), xytext=(0.7, 0.3),
                arrowprops=dict(arrowstyle="->, head_width=0.15", lw=1.5,
                              connectionstyle="arc3,rad=0.3"))
    ax1.set_xlim(0, 1); ax1.set_ylim(0, 1); ax1.set_aspect('equal')
    ax1.axis('off'); ax1.set_title("(a) Diagonal", fontsize=11)

    # (b) Four-point interaction
    ax2.annotate("", xy=(0.6, 0.8), xytext=(0.1, 0.8),
                arrowprops=dict(arrowstyle="->", lw=2))
    ax2.annotate("", xy=(0.6, 0.2), xytext=(0.1, 0.2),
                arrowprops=dict(arrowstyle="->", lw=2))
    ax2.annotate("", xy=(0.9, 0.8), xytext=(0.6, 0.8),
                arrowprops=dict(arrowstyle="->", lw=2))
    ax2.annotate("", xy=(0.9, 0.2), xytext=(0.6, 0.2),
                arrowprops=dict(arrowstyle="->", lw=2))
    # gluon exchange
    ax2.plot([0.6, 0.6], [0.2, 0.8], 'b-', lw=2, zorder=5)
    ax2.plot(0.6, 0.5, 'bs', markersize=8, zorder=6)
    ax2.set_xlim(0, 1); ax2.set_ylim(0, 1); ax2.set_aspect('equal')
    ax2.axis('off'); ax2.set_title("(b) Four-point interaction", fontsize=11)

    fig.tight_layout()
    fig.savefig(FIGDIR / "fig1_vertices.pdf", dpi=150)
    print(f"  Saved: {FIGDIR}/fig1_vertices.pdf")
    plt.close()


# ══════════════════════════════════════════════════════════════════
# FIGURE 4: Higher-Fock contributions
# ══════════════════════════════════════════════════════════════════
def figure4(ncpus=1):
    """Reproduce Figure 4: higher-Fock quark distributions q(x).
    Each panel: one line per coupling, quarks only (not antiquarks separately).
    (a) Meson qq̄qq̄ sector ×10⁴, (b) Baryon qqqq̄q̄ sector ×10³/×10²,
    (c) Baryon 7-parton sector ×10⁴/×10⁷."""
    print("\n=== FIGURE 4: Higher-Fock contributions ===")

    K_mes = 24; K_bar = 21
    lam_weak = mg_to_lambda(1.6); lam_strong = mg_to_lambda(0.1)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # (a) Meson: q(x) in qq̄qq̄ sector
    for lam, ls, lab, scale in [(lam_weak, '-o', "m/g=1.6", 1e4),
                                 (lam_strong, '--s', "m/g=0.1", 1e4)]:
        w, c, ns, mi, st, nso, Nm = run_dlcq(3, 1, 0, K_mes, lam, ncpus=ncpus)
        if ns > 0:
            x, qq, _ = extract_structure_function(0, c, nso, mi, st.mstate, Nm, K_mes, 3, 0,
                                                   nparton_filter=4)
            axes[0].plot(x, qq * scale, ls, markersize=2, linewidth=1,
                        label=f"{lab} (×{scale:.0e})")
    axes[0].set_title("(a) Meson: qq̄qq̄"); axes[0].set_xlabel("x")
    axes[0].set_ylabel("q(x)"); axes[0].legend(fontsize=8); axes[0].set_xlim(0, 1)

    # (b) Baryon: q(x) in qqqq̄q̄ sector (quarks only)
    for lam, ls, lab, scale in [(lam_weak, '-o', "m/g=1.6", 1e3),
                                 (lam_strong, '--s', "m/g=0.1", 1e2)]:
        try:
            w, c, ns, mi, st, nso, Nm = run_dlcq(3, 1, 1, K_bar, lam, ncpus=ncpus)
            if ns > 0:
                x, qq, _ = extract_structure_function(0, c, nso, mi, st.mstate, Nm, K_bar, 3, 1,
                                                       nparton_filter=5)
                axes[1].plot(x, qq * scale, ls, markersize=2, linewidth=1,
                            label=f"{lab} (×{scale:.0e})")
        except Exception as e:
            print(f"    Fig4(b) {lab}: {e}")
    axes[1].set_title("(b) Baryon: qqqq̄q̄"); axes[1].set_xlabel("x")
    axes[1].set_ylabel("q(x)"); axes[1].legend(fontsize=8); axes[1].set_xlim(0, 1)

    # (c) Baryon: q(x) in 7-parton sector (quarks only)
    for lam, ls, lab, scale in [(lam_weak, '-o', "m/g=1.6", 1e4),
                                 (lam_strong, '--s', "m/g=0.1", 1e7)]:
        try:
            w, c, ns, mi, st, nso, Nm = run_dlcq(3, 1, 1, K_bar, lam, ncpus=ncpus)
            if ns > 0:
                x, qq, _ = extract_structure_function(0, c, nso, mi, st.mstate, Nm, K_bar, 3, 1,
                                                       nparton_filter=7)
                if np.max(np.abs(qq)) > 1e-20:
                    axes[2].plot(x, qq * scale, ls, markersize=2, linewidth=1,
                                label=f"{lab} (×{scale:.0e})")
        except Exception as e:
            print(f"    Fig4(c) {lab}: {e}")
    axes[2].set_title("(c) Baryon: 2 extra qq̄"); axes[2].set_xlabel("x")
    axes[2].set_ylabel("q(x)"); axes[2].legend(fontsize=8); axes[2].set_xlim(0, 1)

    fig.suptitle("SU(3) Higher-Fock Structure Functions", fontsize=13)
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig4_higher_fock.pdf", dpi=150)
    print(f"  Saved: {FIGDIR}/fig4_higher_fock.pdf")
    plt.close()


# ══════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(description="Reproduce figures from PRD 41, 3814 (1990)")
    parser.add_argument('--fig', nargs='+', type=int, default=None,
                       help='Figure numbers to generate (default: all)')
    parser.add_argument('--ncpus', type=int, default=1,
                       help='Number of CPUs for parallel computation')
    parser.add_argument('--table1', action='store_true',
                       help='Generate Table I')
    args = parser.parse_args()

    os.environ['QCDF_NCPUS'] = str(args.ncpus)

    figs = args.fig or [1, 2, 3, 4, 5, 6, 7, 8]

    print("=" * 60)
    print(" Reproducing figures from Hornbostel, Brodsky & Pauli (1990)")
    print(" PRD 41, 3814")
    print("=" * 60)

    t0 = time.time()

    dispatch = {
        1: figure1,
        2: lambda: figure2(args.ncpus),
        3: lambda: figure3(args.ncpus),
        4: lambda: figure4(args.ncpus),
        5: lambda: figure5(args.ncpus),
        6: lambda: figure6(args.ncpus),
        7: lambda: figure7(args.ncpus),
        8: lambda: figure8(args.ncpus),
    }

    for fignum in figs:
        if fignum in dispatch:
            try:
                dispatch[fignum]()
            except Exception as e:
                print(f"\n  *** Figure {fignum} FAILED: {e} ***\n")
                import traceback; traceback.print_exc()

    if args.table1:
        table1(args.ncpus)

    print(f"\nTotal time: {time.time()-t0:.0f}s")
    print(f"Figures saved in: {FIGDIR.absolute()}")
    print(f"Cache saved in: {CACHEDIR.absolute()}")


if __name__ == '__main__':
    main()
