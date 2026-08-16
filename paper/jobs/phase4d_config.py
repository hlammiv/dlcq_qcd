"""Phase 4d: the exotics run list — tetraquark and hexaquark searches.

Design constraints from the plan:
  * Fock weights need eigenvectors and the norm, so every run keeps both
    (the provider does by default) with nev deep enough to cover the
    threshold region.
  * Binding energies come from same-K correlated differences ONLY
    (dlcq.fock_weights.binding_series refuses anything else), so every
    candidate sector runs on the same K grid as its threshold hadrons.
  * The hexaquark study runs at N=4: B=2, B=1 and B=0 all live on even
    K_code there, so state and threshold share grids exactly.  At N=3 the
    B=1 grid is odd and a same-K difference cannot exist.
  * DLCQ_MXTRM must be raised in the environment for the deep-LPN N=4
    sectors (the runner's shell exports it).
"""

MGS = [1.6, 0.8, 0.4, 0.1]
KS_EVEN = [20, 24, 28, 32, 36, 40, 44, 48]

RUNS = []

# Tetraquark search, B=0: spectra deep enough in Fock space to hold a
# 4q-dominant state, plus the same-K meson ground state as the threshold.
for N in (2, 3):
    for K in KS_EVEN:
        for mg in MGS:
            RUNS.append((f"tetra_N{N}", dict(N=N, NF=1, B=0, K_code=K,
                                             LPN=6, nev=32), ("mg", mg)))

# Hexaquark study at N=4: the B=2 state, the single baryon, and the meson,
# all on one even-K grid.  LPN = valence + 2 throughout.
for K in KS_EVEN:
    for mg in MGS:
        RUNS.append(("hexa_B2", dict(N=4, NF=1, B=2, K_code=K, LPN=10,
                                     nev=12), ("mg", mg)))
        RUNS.append(("hexa_B1", dict(N=4, NF=1, B=1, K_code=K, LPN=6,
                                     nev=8), ("mg", mg)))
        RUNS.append(("hexa_B0", dict(N=4, NF=1, B=0, K_code=K, LPN=6,
                                     nev=8), ("mg", mg)))
