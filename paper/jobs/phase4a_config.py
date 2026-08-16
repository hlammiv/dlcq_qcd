"""Phase 4a: the pinned run list for the paper's large-K figure panels.

One source of truth, imported by both the cache-warming wrapper
(phase4a_lenore.py) and, later, paper/make_figures.py — so the runs that
warm the cache and the keys the figures read can never drift apart.

Conventions: every run is sparse+blockwise+exact (the scaling path), with
nev eigenvectors kept for structure functions.  Parities follow
figures._K_grid; resolutions are the "2026 panels" of the reproduction
protocol.  LPN choices: mesons valence+4 where higher-Fock content is
plotted (figs 2, 5), valence+2 where only the valence is (fig 3); baryons
valence+2; the B=2 state valence+2 (LPN=8).
"""

# The 1990 figure-2 coupling grid, restated from dlcq.figures.figure2's
# default (explicit, not arange, for exact cache matching).
LAMBDAS = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
           0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95,
           0.97, 0.99]

MG_STRONG, MG_WEAK = 0.1, 1.6          # strong/weak COUPLING (small/large m/g)

# (label, dict(N, NF, B, K_code, LPN, nev), coupling spec)
# coupling spec: ("lambda", value) uses rlamb directly; ("mg", value) goes
# through figures.paper_lambda.
RUNS = []

# Fig. 2, 2026 panels: spectra vs coupling.  Resolution is traded for Fock
# depth here on purpose: at small lambda a state of L partons sits at
# M^2/(m^2+g^2/pi) ~ L^2, and the 1990 panels are untruncated, so matching
# their y-range (50/60/100) needs LPN >= 8/9/10 — a starved LPN at larger K
# cannot contain the upper half of the original panel at all (first drawn
# that way at 2K=50, LPN=6, and misread immediately).  nev=250 covers the
# full range: the 250th level clears the panel top (probed: 52.5 at
# 2K=34, LPN=8, lambda=0.5).
for lam in LAMBDAS:
    RUNS.append(("fig2_B0", dict(N=3, NF=1, B=0, K_code=34, LPN=8, nev=250),
                 ("lambda", lam)))
    RUNS.append(("fig2_B1", dict(N=3, NF=1, B=1, K_code=35, LPN=9, nev=400),
                 ("lambda", lam)))
    RUNS.append(("fig2_B2", dict(N=3, NF=1, B=2, K_code=34, LPN=10, nev=350),
                 ("lambda", lam)))

# Fig. 3, 2026 panels: valence structure functions at 2K = 70/71.
for mg in (MG_WEAK, MG_STRONG):
    RUNS.append(("fig3_mes", dict(N=3, NF=1, B=0, K_code=70, LPN=4, nev=4),
                 ("mg", mg)))
    RUNS.append(("fig3_bar", dict(N=3, NF=1, B=1, K_code=71, LPN=5, nev=4),
                 ("mg", mg)))

# Fig. 5, 2026 panel: meson excitations with 4q content, 2K = 60.
RUNS.append(("fig5", dict(N=3, NF=1, B=0, K_code=100, LPN=6, nev=40),
             ("mg", MG_WEAK)))

# Fig. 6, 2026 panels: baryon states (also serves thesis fig 24 splitting),
# and the two-baryon state with its 8-parton sector.
RUNS.append(("fig6_bar", dict(N=3, NF=1, B=1, K_code=99, LPN=5, nev=12),
             ("mg", MG_WEAK)))
RUNS.append(("fig6_B2", dict(N=3, NF=1, B=2, K_code=50, LPN=8, nev=8),
             ("mg", MG_WEAK)))

# Thesis fig. 22: the K-scan of the N=2 meson four-quark content.
for K in (14, 20, 28, 40, 56, 80, 100):
    for mg in (MG_WEAK, MG_STRONG):
        RUNS.append(("fig22", dict(N=2, NF=1, B=0, K_code=K, LPN=4, nev=4),
                     ("mg", mg)))
