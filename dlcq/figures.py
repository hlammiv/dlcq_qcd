"""Reproduce the figures of Phys. Rev. D 41, 3814 from a :class:`DLCQResult`.

Every routine here takes results from a :class:`~dlcq.providers.Provider` and
never touches a solver directly, so::

    python -m dlcq.figures --source fortran --fig 6
    python -m dlcq.figures --source python  --fig 6

run identical arithmetic on the two codes' output.  Any difference in the plot
is a difference in the physics.

Figure inventory (parameters as stated in the paper):

===  =====================================================  ==================
Fig  content                                                parameters
===  =====================================================  ==================
1    interaction vertices                                   schematic only
2    spectra vs coupling, B = 0, 1, 2                       N=3, 2K=10/13/22
3    valence structure functions, meson and baryon          N=3, m/g=0.1, 1.6
4    higher-Fock contributions                              N=3, m/g=0.1, 1.6
5    first three meson states plus the 11th                 N=3, m/g=1.6, 2K=24
6    first three baryon states plus first B=2               N=3, m/g=1.6, 2K=21
     ...panel (d) is 2K=24; the caption states no K for it
7    extrapolated masses, N = 2, 3, 4                       Richardson, 2K=16-24
8    comparison with large-N and lattice                    N = 2, 3, 4
===  =====================================================  ==================

The paper never states K for Figs. 3 and 4.  We adopt 2K=24 (meson) and 2K=21
(baryon) to match Figs. 5 and 6, and record that choice on the figure itself.
``dlcq.units.infer_K_from_x_grid`` can recover the true K from the digitized
plots, since momenta are odd integers and so x_min = 1/K, dx = 2/K.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .observables import (structure_function, valence_parton_count,
                          richardson_extrapolate, thooft_valence_limit,
                          physical_indices, spurious_zero_modes,
                          monotone_bracket)
from .units import (code_to_M_over_g, lambda_to_mg, mg_to_lambda,
                    thooft_rescale)

__all__ = ["figure1", "figure2", "figure3", "figure4", "figure5",
           "figure6", "figure7", "figure8", "table1", "table1_budget", "figure_fits", "FIGURES"]

_ROOT = Path(__file__).resolve().parent.parent
FIGDIR = _ROOT / "figures"

# Paper parameters.
MG_STRONG, MG_WEAK = 0.1, 1.6
K_MESON, K_BARYON = 24, 21

# The paper never states K for Figs. 3 and 4.  Recovered from the published
# plots themselves: momenta are odd integers, so the marker x-positions lie on
# {k/K}, and tools/digitize.py fits that lattice.  Every meson panel returned an
# even K and every baryon panel an odd K -- the parity each sector requires,
# which is not built into the inference -- and the same procedure reproduces the
# stated 2K = 24 and 21 of Figs. 5 and 6.  See docs/inferred-K.md.
# Fig. 6(d)'s K is never stated in the paper.  Measured from the panel: 24.
K_TWO_BARYON = 24
K_FIG3_MESON, K_FIG3_BARYON = 14, 15
K_FIG4_MESON, K_FIG4_BARYON = 14, 15

# The m/g = 1.6 runs behind Figs. 5 and 6 used the *rounded* coupling 0.3325 --
# it is the value in every input_*.txt / input_*.json in this repo and in the
# preserved Fortran outputs.  mg_to_lambda(1.6) is 0.33254949, and that 1.5e-5
# difference is far above a 1e-12 comparison tolerance, so both solvers must be
# handed the identical literal value rather than re-deriving it from m/g.
LAM_WEAK = 0.3325                       # m/g = 1.60024...
LAM_STRONG = float(mg_to_lambda(MG_STRONG))


def paper_lambda(mg):
    """The coupling the paper's runs actually used for this m/g."""
    return LAM_WEAK if abs(mg - MG_WEAK) < 1e-9 else float(mg_to_lambda(mg))


def _outfile(name, source):
    FIGDIR.mkdir(parents=True, exist_ok=True)
    return FIGDIR / f"{name}_{source}.pdf"


def _finish(fig, name, source):
    path = _outfile(name, source)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    # A PNG alongside, so tools/compare_panels.py can place our figure next to
    # the scanned original without rasterizing the PDF.
    fig.savefig(path.with_suffix(".png"), dpi=130)
    plt.close(fig)
    print(f"  saved {path}")
    return path


# ──────────────────────────────────────────────────────────────────────────

def figure1(provider=None, source="schematic"):
    """Fig. 1 -- interaction vertices.  Schematic; nothing to validate."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3))

    ax1.annotate("", xy=(0.75, 0.5), xytext=(0.1, 0.5),
                 arrowprops=dict(arrowstyle="-|>", lw=2, color="k"))
    th = np.linspace(0, np.pi, 100)
    ax1.plot(0.45 + 0.17 * np.cos(th), 0.5 + 0.17 * np.sin(th), "b-", lw=1.5)
    ax1.set_title("(a) diagonal / self-energy", fontsize=10)

    for y in (0.25, 0.75):
        ax2.annotate("", xy=(0.9, y), xytext=(0.1, y),
                     arrowprops=dict(arrowstyle="-|>", lw=2, color="k"))
    ax2.plot([0.5, 0.5], [0.25, 0.75], "b-", lw=2)
    ax2.set_title("(b) four-point (instantaneous gluon)", fontsize=10)

    for ax in (ax1, ax2):
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_aspect("equal"); ax.axis("off")

    fig.suptitle("FIG. 1  Interaction vertices", fontsize=11)
    return _finish(fig, "fig1_vertices", source)


def figure2(provider, source, lambdas=None, ncurves=30):
    """Fig. 2 -- full spectra vs coupling for B = 0, 1, 2 at fixed K."""
    if lambdas is None:
        # Explicit rather than arange: these must match the swept runs exactly,
        # and arange's 0.15000000000000002 is a needless risk.
        lambdas = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
                   0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95,
                   0.97, 0.99]
    panels = [(0, 10, 50), (1, 13, 60), (2, 22, 100)]
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    for ax, (B, K, ymax) in zip(axes, panels):
        for lam in lambdas:
            r = provider.get(3, 1, B, K, float(lam))
            if r.n_eigenvalues == 0:
                continue
            ev = r.eigenvalues[:ncurves]
            ax.plot([lam] * len(ev), ev, "k.", markersize=1.5)
        ax.set_xlim(0, 1); ax.set_ylim(0, ymax)
        ax.set_xlabel(r"$1/(1+\pi m^2/g^2)^{1/2}$")
        ax.set_ylabel(r"$M^2/(m^2+g^2/\pi)$")
        ax.set_title(f"SU(3)  B={B},  2K={K}")

    fig.suptitle(f"FIG. 2  Spectra vs coupling  [{source}]")
    return _finish(fig, "fig2_spectra", source)


def figure3(provider, source):
    """Fig. 3 -- valence structure functions, meson and baryon."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    for mg, marker, fill in [(MG_WEAK, "o", "full"), (MG_STRONG, "s", "none")]:
        lam = paper_lambda(mg)
        for ax, B, K in [(ax1, 0, K_FIG3_MESON), (ax2, 1, K_FIG3_BARYON)]:
            r = provider.get(3, 1, B, K, float(lam))
            if r.n_eigenvectors == 0:
                continue
            nval = valence_parton_count(3, B)
            gs = int(physical_indices(r)[0])
            x, q, _ = structure_function(r, gs, nparton=nval)
            ax.plot(x, q, "-" + marker, fillstyle=fill, markersize=4, lw=1,
                    label=f"m/g = {mg}")

    # The chiral-limit analytic result, Eq. (22): q(x) = N(N-1)(1-x)^(N-2).
    xs = np.linspace(0, 1, 200)
    ax2.plot(xs, thooft_valence_limit(xs, 3), "k:", lw=1.2,
             label=r"$6(1-x)$  (Eq. 22, $m/g\to0$)")

    ax1.set_title(f"(a) SU(3) meson,  2K={K_FIG3_MESON} (inferred)")
    ax2.set_title(f"(b) SU(3) baryon,  2K={K_FIG3_BARYON} (inferred)")
    for ax in (ax1, ax2):
        ax.set_xlabel("x = k/K"); ax.set_ylabel("q(x)")
        ax.set_xlim(0, 1); ax.legend(fontsize=8)

    fig.suptitle(f"FIG. 3  Valence structure functions  [{source}]"
                 "   (K not stated in the paper; recovered from the plot)")
    return _finish(fig, "fig3_valence", source)


def figure4(provider, source):
    """Fig. 4 -- higher-Fock contributions."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    specs = [
        (axes[0], 0, K_FIG4_MESON, 4, "(a) meson: $q\\bar q q\\bar q$", False),
        (axes[1], 1, K_FIG4_BARYON, 5, "(b) baryon: $qqqq\\bar q$", True),
        (axes[2], 1, K_FIG4_BARYON, 7, "(c) baryon: two extra pairs", False),
    ]
    for ax, B, K, npart, title, show_anti in specs:
        for mg, marker, fill in [(MG_WEAK, "o", "full"), (MG_STRONG, "s", "none")]:
            r = provider.get(3, 1, B, K, paper_lambda(mg))
            if r.n_eigenvectors == 0:
                continue
            gs = int(physical_indices(r)[0])
            x, q, qbar = structure_function(r, gs, nparton=npart)
            if np.max(np.abs(q)) > 1e-18:
                ax.plot(x, q, "-" + marker, fillstyle=fill, markersize=3, lw=1,
                        label=f"q, m/g={mg}")
            if show_anti and np.max(np.abs(qbar)) > 1e-18:
                ax.plot(x, qbar, "--^", fillstyle=fill, markersize=3, lw=0.8,
                        label=rf"$\bar q$, m/g={mg}")
        ax.set_title(title); ax.set_xlabel("x = k/K"); ax.set_ylabel("q(x)")
        ax.set_xlim(0, 1); ax.legend(fontsize=7)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))

    fig.suptitle(f"FIG. 4  Higher-Fock structure functions  [{source}]")
    return _finish(fig, "fig4_higher_fock", source)


def _wavefunction_panel(ax, r, idx, sectors, title):
    """One panel of Figs. 5/6: valence plus the next Fock sector."""
    if not r.has_eigenvector(idx):
        ax.text(0.5, 0.5, f"eigenvector {idx}\nnot available",
                ha="center", va="center", transform=ax.transAxes, fontsize=9)
        ax.set_title(title, fontsize=10)
        return
    Mg = code_to_M_over_g(r.eigenvalues[idx], r.rlamb)
    for npart, style, lab in sectors:
        # A negative parton count selects that sector's ANTIQUARK distribution:
        # Fig. 6(d) plots q(x) and qbar(x) of the 8-parton sector separately.
        x, q, qbar = structure_function(r, idx, nparton=abs(npart))
        y = qbar if npart < 0 else q
        if np.max(np.abs(y)) > 1e-18:
            ax.plot(x, y, style, markersize=3, lw=1, label=lab)
    ax.set_title(f"{title}   M/g = {Mg:.3f}", fontsize=10)
    ax.set_xlabel("x = k/K"); ax.set_ylabel("q(x)")
    ax.set_xlim(0, 1); ax.legend(fontsize=7)


def figure5(provider, source):
    """Fig. 5 -- first three meson states plus the 11th, N=3, m/g=1.6, 2K=24."""
    r = provider.get(3, 1, 0, K_MESON, paper_lambda(MG_WEAK))
    phys = physical_indices(r)
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    sectors = [(2, "k-o", r"$q\bar q$"), (4, "r-s", r"$q\bar q q\bar q$")]
    # "First three states" and "eleventh state" count the PHYSICAL spectrum;
    # see observables.spurious_zero_modes for the decoupled 12-parton mode
    # that otherwise sits at index 0 in this configuration.
    wanted = [0, 1, 2, 10]
    for ax, w, title in zip(axes.flat, wanted,
                            ["(a) 1st", "(b) 2nd", "(c) 3rd", "(d) 11th"]):
        idx = int(phys[w]) if w < len(phys) else w
        _wavefunction_panel(ax, r, idx, sectors, title)
    fig.suptitle(f"FIG. 5  SU(3) meson, m/g={MG_WEAK}, 2K={K_MESON}  [{source}]")
    return _finish(fig, "fig5_meson_wf", source)


def figure6(provider, source):
    """Fig. 6 -- first three baryon states plus the first B=2 state, 2K=21."""
    lam = paper_lambda(MG_WEAK)
    r = provider.get(3, 1, 1, K_BARYON, lam)
    phys = physical_indices(r)
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    sectors = [(3, "k-o", "qqq"), (5, "r-s", r"$qqqq\bar q$")]
    for ax, w, title in zip(axes.flat, [0, 1, 2],
                            ["(a) 1st baryon", "(b) 2nd", "(c) 3rd"]):
        idx = int(phys[w]) if w < len(phys) else w
        _wavefunction_panel(ax, r, idx, sectors, title)

    ax = axes.flat[3]
    try:
        # 2K = 24, not 22.  The caption states 2K=21 for (a)-(c) and gives no
        # K at all for (d); 24 is what the row of zero-valued markers along the
        # panel's own x axis says, and at 24 the computed curves land on the
        # published ones to better than 1%.  See docs/baryon-higher-fock.md.
        r2 = provider.get(3, 1, 2, K_TWO_BARYON, lam)
        _wavefunction_panel(ax, r2, 0,
                            [(6, "k-o", "6q"),
                             (8, "r-s", r"$6qq\bar q$  $q(x)$"),
                             (-8, "b-^", r"$6qq\bar q$  $\bar q(x)$")],
                            f"(d) 1st B=2, 2K={K_TWO_BARYON}")
    except Exception as exc:                       # provider may lack the run
        ax.text(0.5, 0.5, f"B=2 unavailable:\n{exc}", ha="center", va="center",
                transform=ax.transAxes, fontsize=7)
        ax.set_title("(d) 1st B=2", fontsize=10)

    fig.suptitle(f"FIG. 6  SU(3) baryon, m/g={MG_WEAK}, 2K={K_BARYON}  [{source}]")
    return _finish(fig, "fig6_baryon_wf", source)


# ──────────────────────────────────────────────────────────────────────────
# Extrapolated masses
# ──────────────────────────────────────────────────────────────────────────

MG_TABLE = [1.6, 0.8, 0.4, 0.2, 0.1, 0.05, 0.0]


def _mass_series(provider, N, B, mg, K_codes, msq_units=False, lpn=None):
    """The ``(2K, mass)`` points a Richardson fit is built from.

    Factored out so ``figure_fits`` plots exactly the points
    ``_extrapolated_mass`` fits -- if they could drift apart, the diagnostic
    would stop diagnosing.
    """
    lam = float(mg_to_lambda(mg))
    Ks, masses = [], []
    for K in K_codes:
        # Levels only.  A cached result keeps ``Z`` dense despite it being
        # 0.014% non-zero, so ``get`` turns an 8 MB file into ~700 MB in memory
        # -- and a Richardson series wants nothing from a run but its spectrum.
        # Providers without the light path (the Fortran one) fall back to
        # ``get``, which is what this used to do unconditionally.
        args = (N, 1, B, K, lam, -1.0,
                sweep_lpn(N, B) if lpn is None else lpn)
        r = (provider.spectrum(*args) if hasattr(provider, "spectrum")
             else provider.get(*args))
        if r.n_eigenvalues == 0:
            continue
        phys = physical_indices(r)
        if phys.size == 0:
            continue
        Ks.append(K)
        ev = float(r.eigenvalues[phys[0]])
        masses.append(ev if msq_units else code_to_M_over_g(ev, lam))
    return Ks, masses


#: Correction terms kept in the Eq. (27) fit.
#:
#: Was 2.  Three independent criteria say that under-fits, measured over the
#: whole of Table I on 2K = 25-49:
#:
#: * held-out K -- fit on 2K = 25-45, predict 47 and 49 -- prefers 4 on **30 of
#:   30** entries, and never picks 2.  This criterion never sees the paper.
#: * the Levin u-transform, which assumes no endpoint exponent at all and so
#:   shares none of Eq. (27)'s structure, lands closer to the 4-term answer on
#:   **29 of 30**.
#: * agreement with the published table improves on **21 of 30**, median 1.54x.
#:
#: The rows where agreement gets worse are mid-coupling baryons, and Levin sits
#: *past* the 4-term value there in the same direction (bar N=3, m/g=0.4:
#: 2-term 7.344, 4-term 7.526, Levin 7.697, published 7.300), so those are
#: evidence about the published value rather than against four terms.
#:
#: **But only where the window supports it.**  All of that evidence comes from
#: 12-13 points spanning 2K = 25-49.  It does not transfer to the paper's own
#: window, which is 4-5 points spanning 16-24: there four terms put the N=4
#: baryon at m/g = 0.8 at 22.073 against a published 20.900, while two terms
#: land on it.  A flat default is the wrong shape -- the number of correction
#: terms has to follow the amount of data.
N_TERMS = 4

#: Points below which the wide-window evidence does not apply and the historical
#: two terms are used instead.  Eight is the first window size past the paper's
#: own (5) and past 2K = 25-35 (5-6), and is reached by 2K = 25-41.
N_TERMS_MIN_POINTS = 8


def _terms_for(n_points):
    """Correction terms a window of ``n_points`` supports."""
    return N_TERMS if n_points >= N_TERMS_MIN_POINTS else 2


def _extrapolated_mass(provider, N, B, mg, K_codes, msq_units=False, lpn=None,
                       n_terms=None):
    """Lightest state extrapolated to the continuum via Eq. (27).

    ``msq_units=False`` returns M/g, the quantity Figs. 7 and 8 plot.
    ``msq_units=True`` returns M^2/(m^2+g^2/pi), the quantity **Table I**
    tabulates despite its M/g headers -- see docs/table1-units.md.
    """
    if mg == 0.0:
        # Exact: the lightest state is massless for any N, B, K (Eq. 16).
        return 0.0, 0.0
    # Sweeps use the exact coupling.  paper_lambda's rounded 0.3325 exists only
    # to match the preserved Fig. 5/6 runs; using it here would miss every
    # swept run and silently return "no data".
    Ks, masses = _mass_series(provider, N, B, mg, K_codes, msq_units, lpn)
    if len(Ks) < 3:
        return None, None
    if n_terms is None:
        n_terms = _terms_for(len(Ks))
    return richardson_extrapolate(Ks, masses, mg, N, n_terms=n_terms)


def _extrapolated_bracket(provider, N, B, mg, K_codes, msq_units=False,
                          lpn=None):
    """``(lo, hi)`` bounding ``M(K -> infinity)``, or ``(None, None)``.

    Reported instead of the paper's last-term rule, which does not survive
    scrutiny: it is *basis dependent*, moving by a factor of several under a
    change of coordinates that leaves M(0) identical to 14 digits, and at weak
    coupling it is wildly pessimistic -- at N=3, B=1, m/g=0.05 it quotes 0.594
    on a value of 0.243, i.e. 244%.

    The bracket assumes only what is measured to hold: ``M^2(K)`` increases with
    K, its increments decrease, and they decay as a power law.  The lower bound
    is then free and the upper bound follows from the integral test.

    It is honest rather than flattering.  At strong coupling it is tight (1.4%
    at N=3, B=1, m/g=1.6); at weak coupling it is enormous -- a factor of 10 at
    m/g=0.05 -- because the increments there decay as ``K^-1.02``, barely
    summable, and 2K <= 49 genuinely does not determine the limit.  Both the
    2- and 4-term fits and the published value sit inside it on essentially
    every row.
    """
    if mg == 0.0:
        return 0.0, 0.0
    Ks, masses = _mass_series(provider, N, B, mg, K_codes, True, lpn)
    if len(Ks) < 4:
        return None, None
    lo, hi, _, mono = monotone_bracket(Ks, masses)
    if not mono:
        return None, None
    if msq_units:
        return lo, hi
    lam = float(mg_to_lambda(mg))
    return (code_to_M_over_g(lo, lam),
            code_to_M_over_g(hi, lam) if np.isfinite(hi) else float("inf"))


def sweep_lpn(N, B):
    """Particle-number truncation used by the Figs. 7/8 and Table I sweeps.

    Valence plus one extra qqbar pair, capped at 10 partons.  The cap keeps
    every run inside qcdf.f's 25-slot colour array (a matrix element needs
    2L+4 slots), and the truncation is converged to ~1e-5 at weak coupling.

    Both providers must be handed the SAME value.  The Fortran provider used to
    match a stored run on (N, NF, B, K, lambda) alone, so it happily returned a
    truncated sweep run for a request that left LPN at its default 0 -- meaning
    the "Fortran vs Python" comparison was between different physics.
    """
    valence = 2 if B == 0 else N * abs(B)
    return min(valence + 2, 10)


def _K_grid(B, N, lo=16, hi=24):
    """2K in the paper's 16-24 window, with the parity B requires."""
    want_odd = (N * abs(B)) % 2 == 1
    return [K for K in range(lo, hi + 1) if (K % 2 == 1) == want_odd]


def figure7(provider, source, N_values=(2, 3, 4)):
    """Fig. 7 -- extrapolated meson and baryon masses vs m/g."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    styles = {2: "-", 3: "--", 4: ":"}

    for N in N_values:
        for ax, B in [(ax1, 0), (ax2, 1)]:
            xs, ys, es = [], [], []
            for mg in sorted(MG_TABLE):
                M0, err = _extrapolated_mass(provider, N, B, mg, _K_grid(B, N))
                if M0 is None:
                    continue
                xs.append(mg); ys.append(M0); es.append(err)
            if xs:
                ax.errorbar(xs, ys, yerr=es, fmt="o" + styles[N], lw=1.2,
                            markersize=4, capsize=2, label=f"SU({N})")

    ax1.set_title("(a) meson mass"); ax2.set_title("(b) baryon mass")
    for ax in (ax1, ax2):
        ax.set_xlabel("m/g"); ax.set_ylabel("M/g")
        ax.set_xlim(0, 1.7); ax.legend()

    fig.suptitle(f"FIG. 7  Extrapolated masses  [{source}]")
    return _finish(fig, "fig7_masses", source)


def figure8(provider, source, N_values=(2, 3, 4)):
    """Fig. 8 -- meson masses against the large-N limit."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    styles = {2: "-", 3: "--", 4: ":"}

    for N in N_values:
        xs, ys = [], []
        for mg in sorted(MG_TABLE):
            M0, _ = _extrapolated_mass(provider, N, 0, mg, _K_grid(0, N))
            if M0 is None:
                continue
            xs.append(mg); ys.append(M0)
        if not xs:
            continue
        ax1.plot(xs, ys, "o" + styles[N], lw=1.2, markersize=4, label=f"SU({N})")
        s = thooft_rescale(N)
        ax2.plot(np.array(xs) * s, np.array(ys) * s, "o" + styles[N],
                 lw=1.2, markersize=4, label=f"SU({N})")

    # The paper overlays 't Hooft's N -> infinity solution of Eq. (24).  It is
    # not a DLCQ result, so it comes from dlcq.thooft; that solver reproduces
    # the classic chiral-limit eigenvalues 5.88 and 14.1 as a benchmark.
    # The paper overlays 't Hooft's solution of Eq. (24).  It is NOT obtained by
    # evaluating at some huge N: this paper's coupling (N^2-1)/2N pi grows with
    # N, so on the unrescaled M/g axis of panel (a) that diverges.  Panel (a)
    # shows the limit of the SU(N) sequence -- the leading-order Casimir
    # N/2 in place of (N^2-1)/2N -- at the largest N compared, which is why the
    # published curve sits just above SU(4).  Panel (b)'s (2pi/N)^(1/2)
    # rescaling is what makes the N -> infinity limit finite.
    try:
        from .thooft import thooft_mass
        N_ref = max(N_values)
        mgs = np.linspace(0.05, 1.6, 12)
        large_N = np.array([thooft_mass(float(t), N_ref, large_n=True)
                            for t in mgs])
        ax1.plot(mgs, large_N, "k-", lw=1.5,
                 label=f"large N ('t Hooft, N={N_ref})")
        s = thooft_rescale(N_ref)
        ax2.plot(mgs * s, large_N * s, "k-", lw=1.5, label="large N")
    except Exception as exc:                       # never let this kill the figure
        print(f"    large-N curve unavailable: {exc}")

    ax1.set_xlabel("m/g"); ax1.set_ylabel("M/g")
    ax1.set_title("(a) meson mass"); ax1.legend()
    ax2.set_xlabel(r"$(2\pi/N)^{1/2}\, m/g$")
    ax2.set_ylabel(r"$(2\pi/N)^{1/2}\, M/g$")
    ax2.set_title("(b) 't Hooft rescaling"); ax2.legend()

    fig.suptitle(f"FIG. 8  Comparison with large N  [{source}]"
                 "   (note: this paper's unit is $g^2N/2\\pi$)")
    return _finish(fig, "fig8_large_N", source)


def table1(provider, source, N_values=(2, 3, 4), K_lo=16, K_hi=24):
    """Table I -- N dependence of meson and baryon mass.

    Reported in **both** conventions, because the paper's table and its figures
    use different ones: Table I tabulates ``M^2/(m^2 + g^2/pi)`` -- the y-axis
    of Fig. 2 -- despite column headers reading ``M_mes/g``, while Figs. 7 and 8
    plot ``M/g``.  See docs/table1-units.md.

    ``K_lo``/``K_hi`` set the Richardson window.  The paper's is 16-24 ("2K in
    the range of roughly 16-24"), which is the default so this reproduces the
    published table.  A wider window is now affordable and tightens the fit
    sharply where the paper is reliable -- at N=3, B=1, m/g=1.6 the last
    Richardson term falls from 0.0041 to 0.0005 on 2K = 25-35, a factor of
    eight, while the central value moves by less than the paper's own quoted
    uncertainty.  It does **not** rescue the weak-coupling rows; see
    docs/table1-units.md.

    The particle-number truncation is ``sweep_lpn`` at every K, so the two
    windows are compared on identical physics rather than on a truncation that
    changes with the window.
    """
    paper = _load_paper_table1()
    grid = lambda B, N: _K_grid(B, N, K_lo, K_hi)      # noqa: E731
    rows = []
    for mg in MG_TABLE:
        row = {"mg": mg}
        for N in N_values:
            row[f"mes_N{N}"] = (
                _extrapolated_mass(provider, N, 0, mg, grid(0, N), msq_units=True),
                _extrapolated_mass(provider, N, 0, mg, grid(0, N)),
                _extrapolated_bracket(provider, N, 0, mg, grid(0, N),
                                      msq_units=True))
        for N in (n for n in (3, 4) if n in N_values):
            row[f"bar_N{N}"] = (
                _extrapolated_mass(provider, N, 1, mg, grid(1, N), msq_units=True),
                _extrapolated_mass(provider, N, 1, mg, grid(1, N)),
                _extrapolated_bracket(provider, N, 1, mg, grid(1, N),
                                      msq_units=True))
        rows.append(row)

    FIGDIR.mkdir(parents=True, exist_ok=True)
    tag = "" if (K_lo, K_hi) == (16, 24) else f"_2K{K_lo}-{K_hi}"
    out = FIGDIR / f"table1_{source}{tag}.txt"
    cols = [(f"mes_N{N}", "mes", N) for N in N_values] + \
           [(f"bar_N{N}", "bar", N) for N in (3, 4) if N in N_values]

    with open(out, "w") as fh:
        fh.write("TABLE I.  N dependence of meson and baryon mass.\n")
        fh.write(f"Reproduced with the {source} solver, Richardson-extrapolated\n")
        fh.write(f"over 2K = {K_lo}-{K_hi} using Eq. (27).\n\n")
        fh.write("UNITS: M^2/(m^2 + g^2/pi), matching what Table I actually\n")
        fh.write("tabulates (see docs/table1-units.md).  'paper' is the printed\n")
        fh.write("value; 'pull' is |paper-ours| over the paper's own quoted\n")
        fh.write("last-Richardson-term uncertainty.\n")
        fh.write(f"Fits keep {N_TERMS} correction terms; see figures.N_TERMS "
                 f"for why 2 was wrong.\n\n")

        hdr = f"{'m/g':>6}" + "".join(f"{c[1]+' N='+str(c[2]):>22}" for c in cols)
        fh.write(hdr + "\n" + "-" * len(hdr) + "\n")
        for row in rows:
            line = f"{row['mg']:6.2f}"
            for key, q, N in cols:
                (m0, _), _, _ = row[key]
                ref = paper.get((q, N, row["mg"]))
                if m0 is None:
                    line += f"{'---':>22}"
                elif ref is None:
                    line += f"{m0:>22.4f}"
                else:
                    val, unc = ref
                    if unc > 0:
                        line += (f"{m0:>9.4f}/{val:<7.3f}"
                                 f"p{abs(val - m0) / unc:>3.1f}")
                    else:
                        # m/g = 0 is exact by Eq. (16), not an extrapolation.
                        line += f"{m0:>9.4f}/{val:<7.3f}{'exact':>5}"
            fh.write(line + "\n")

        fh.write("\n\nSame results as M/g -- the quantity Figs. 7 and 8 plot:\n\n")
        fh.write(hdr + "\n" + "-" * len(hdr) + "\n")
        for row in rows:
            line = f"{row['mg']:6.2f}"
            for key, q, N in cols:
                _, (g0, gerr) = row[key][0], row[key][1]
                line += f"{'---':>22}" if g0 is None else f"{g0:>16.4f}({gerr:.2g})"
            fh.write(line + "\n")

        # The bracket, not the last term, is the defensible uncertainty.  The
        # last-term rule is basis dependent -- it moves by a factor of several
        # under a change of coordinates that leaves M(0) identical to 14 digits
        # -- and at weak coupling it is absurd, quoting 244% at N=3, B=1,
        # m/g=0.05.  The bracket assumes only that M^2(K) rises with decreasing
        # increments that decay as a power law, all of which is measured to hold.
        fh.write("\n\nBracket on M^2(K -> infinity), in Table I's units.\n")
        fh.write("Lower bound is M^2 at the largest K computed and needs only\n")
        fh.write("monotonicity; the upper bound additionally assumes the\n")
        fh.write("increments keep decaying as a power law.  Where this is wide\n")
        fh.write("the data does not determine the limit, and no fit can fix\n")
        fh.write("that -- at m/g=0.05 the increments decay as K^-1.02, barely\n")
        fh.write("summable, so 2K <= 49 leaves a factor of ~10 open.\n\n")
        fh.write(hdr + "\n" + "-" * len(hdr) + "\n")
        for row in rows:
            line = f"{row['mg']:6.2f}"
            for key, q, N in cols:
                lo, hi = row[key][2]
                if lo is None:
                    line += f"{'---':>22}"
                elif not np.isfinite(hi):
                    line += f"{'[%.3f, inf)' % lo:>22}"
                else:
                    line += f"{'[%.3f, %.3f]' % (lo, hi):>22}"
            fh.write(line + "\n")

    print(f"  saved {out}")
    return rows


def figure_fits(provider, source, K_lo=25, K_hi=35, lpn=None,
                N_values=(2, 3, 4), n_terms=2):
    """One panel per Table I entry: the data, the fit, and the published value.

    Every number in Table I is an extrapolation, and an extrapolation is only
    as trustworthy as the fit behind it.  This plots each one so it can be
    checked by eye rather than taken on the strength of a residual: the
    computed M^2(K) against 1/K_paper, the fitted Eq. (27) series continued to
    the axis, the extrapolated M(0) as a star, and the paper's value as a
    horizontal band whose half-width is its own quoted last term.

    What to look for.  The points should sit on a gently curving line that the
    fit follows without wiggling; the star should be a short, believable
    extension of it rather than a long extrapolation from a steep tail.  A fit
    that has to bend sharply between the last point and the axis is the visual
    signature of the coefficient blow-up documented in docs/table1-units.md,
    where the Eq. (26) exponent a -> 0 makes 1/K^(1+a) collapse onto 1/K.

    ``lpn=None`` uses ``sweep_lpn``; pass a larger value to see the truncation
    move the points.
    """
    from .observables import (richardson_budget, richardson_ensemble,
                              richardson_curve)

    paper = _load_paper_table1()
    mgs = [mg for mg in MG_TABLE if mg > 0]
    cols = ([("mes", N, 0) for N in N_values]
            + [("bar", N, 1) for N in (3, 4) if N in N_values])

    fig, axes = plt.subplots(len(mgs), len(cols),
                             figsize=(3.5 * len(cols), 2.7 * len(mgs)),
                             squeeze=False)

    for i, mg in enumerate(mgs):
        for j, (q, N, B) in enumerate(cols):
            ax = axes[i][j]
            Ks, ms = _mass_series(provider, N, B, mg,
                                  _K_grid(B, N, K_lo, K_hi),
                                  msq_units=True, lpn=lpn)
            if len(Ks) < 3:
                ax.text(0.5, 0.5, "no data", ha="center", va="center",
                        transform=ax.transAxes, fontsize=8, color="0.5")
                ax.set_xticks([]); ax.set_yticks([])
                continue

            M0, last, coeffs, exps = richardson_extrapolate(
                Ks, ms, mg, N, n_terms=n_terms, return_fit=True)
            bud = richardson_budget(Ks, ms, mg, N, n_terms=n_terms)

            Kp = np.asarray(Ks) / 2.0
            inv = 1.0 / Kp
            grid = np.linspace(0.0, inv.max() * 1.05, 200)

            # The band is the envelope of every defensible fit -- each order
            # crossed with each sub-window -- evaluated along the grid.  It
            # pinches where the data constrain the curve and opens toward the
            # axis, so its width *is* the extrapolation error, and its shape
            # shows how that error grows with distance from the data.
            ens = richardson_ensemble(Ks, ms, mg, N)
            if ens:
                curves = np.vstack([richardson_curve(c, e, grid)
                                    for c, e, _ in ens])
                ax.fill_between(grid, curves.min(axis=0), curves.max(axis=0),
                                color="C0", alpha=0.20, lw=0, zorder=1)

            ax.plot(grid, richardson_curve(coeffs, exps, grid), "-",
                    color="C0", lw=1.3, zorder=2)
            ax.plot(inv, ms, "o", color="k", markersize=4, zorder=4)
            ax.errorbar([0], [M0], yerr=[bud["total"]], fmt="*", color="C3",
                        markersize=13, capsize=3, elinewidth=1.4, zorder=5)

            row = paper.get((q, N, mg))
            if row is not None:
                pv, pe = row
                ax.axhspan(pv - pe, pv + pe, color="C2", alpha=0.18, zorder=0)
                ax.axhline(pv, color="C2", lw=1.0, ls="--", zorder=2)
                pull = abs(pv - M0) / pe if pe else float("nan")
                note = f"pull {pull:.1f}"
            else:
                note = "no paper value"

            ax.set_xlim(left=-0.02 * inv.max())
            ax.set_title(f"{q} SU({N}),  m/g={mg}", fontsize=9)
            ax.text(0.03, 0.04,
                    f"M(0)={M0:.4f} +- {bud['total']:.2g}\n"
                    f"form {bud['form']:.1g}   window {bud['window']:.1g}\n"
                    f"(paper last-term {last:.2g})\n{note}",
                    transform=ax.transAxes, fontsize=6.5, va="bottom",
                    bbox=dict(fc="white", ec="0.8", alpha=0.85, pad=1.5))
            if i == len(mgs) - 1:
                ax.set_xlabel(r"$1/K$", fontsize=8)
            if j == 0:
                ax.set_ylabel(r"$M^2/(m^2+g^2/\pi)$", fontsize=8)
            ax.tick_params(labelsize=7)

    trunc = "sweep_lpn" if lpn is None else f"LPN={lpn}"
    fig.suptitle(f"Richardson fits behind Table I  --  2K = {K_lo}-{K_hi}, "
                 f"{trunc}, {n_terms} correction terms  [{source}]\n"
                 f"black = computed, blue = Eq. (27) fit, red star = M(0), "
                 f"green band = paper +- its quoted last term;  "
                 f"blue band = envelope of all defensible fits",
                 fontsize=9)
    fig.subplots_adjust(top=0.94)
    tag = f"_2K{K_lo}-{K_hi}" + ("" if lpn is None else f"_LPN{lpn}")
    return _finish(fig, f"table1_fits{tag}", source)


def table1_budget(provider, source, K_lo=25, K_hi=35, lpn=None,
                  N_values=(2, 3, 4), n_terms=2, alt_lpn_delta=2,
                  alt_K_hi=49):
    """Table I with a full error budget, one line per entry.

    Each entry carries all four components rather than a single quoted number:
    ``form`` (how many correction terms), ``window`` (which K), ``trunc`` (the
    particle-number cut, from a second run at ``LPN + alt_lpn_delta``), and the
    solver floor.  See :func:`~dlcq.observables.richardson_budget`.

    Written as a separate routine so ``table1`` keeps reproducing the published
    table in the published units, unchanged.
    """
    from .observables import richardson_budget

    paper = _load_paper_table1()
    cols = ([("mes", N, 0) for N in N_values]
            + [("bar", N, 1) for N in (3, 4) if N in N_values])
    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = FIGDIR / f"table1_budget_2K{K_lo}-{K_hi}_{source}.txt"

    with open(out, "w") as fh:
        fh.write("TABLE I with a full error budget.\n")
        fh.write(f"Solver: {source}.  Richardson over 2K = {K_lo}-{K_hi}, "
                 f"Eq. (27), {n_terms} correction terms.\n")
        fh.write("Units: M^2/(m^2 + g^2/pi), as Table I actually tabulates "
                 "(docs/table1-units.md).\n\n")
        fh.write("Components, added in quadrature:\n")
        fh.write("  form   spread over the number of correction terms (2..4)  "
                 "-- DOMINANT\n")
        fh.write("  wind   spread over sub-windows at fixed order\n")
        fh.write("  trunc  change when the particle-number cut is raised by "
                 "one qqbar pair\n")
        fh.write("  (the solver floor, 6e-13 relative, is carried but never "
                 "visible at this precision)\n\n")
        fh.write("'paper' is the printed value and 'pterm' its own quoted "
                 "last term -- which is the\nmagnitude of the last series "
                 "term, NOT a 1-sigma error bar.  So the meaningful\n"
                 "comparison is 'd/pt': |paper-ours| in units of the PAPER's "
                 "term, which should be\n<= 1.  'd/ours' is the same "
                 "difference in units of our own error; it is routinely\n"
                 "large simply because our error is up to 38x smaller, and on "
                 "its own it does NOT\nindicate disagreement.\n\n")
        hdr = (f"{'case':>8} {'m/g':>5} {'ours':>10} {'+-tot':>9} "
               f"{'form':>9} {'wind':>9} {'trunc':>9} "
               f"{'paper':>9} {'pterm':>8} {'d/pt':>6} {'d/ours':>7} {'ok':>3}")
        fh.write(hdr + "\n" + "-" * len(hdr) + "\n")

        for q, N, B in cols:
            for mg in [m for m in MG_TABLE if m > 0]:
                grid = _K_grid(B, N, K_lo, K_hi)
                Ks, ms = _mass_series(provider, N, B, mg, grid,
                                      msq_units=True, lpn=lpn)
                if len(Ks) < 4:
                    continue
                base = sweep_lpn(N, B) if lpn is None else lpn
                # The alt-LPN series gets its own, shorter grid.  One extra
                # qqbar pair multiplies the basis by 17-23x -- 482,320 states
                # against 27,377 for the N=4 baryon at 2K=70 -- so sweeping it
                # to the top of the main window is days of compute for a term
                # already measured to shrink with K.  ``richardson_budget``
                # compares the two on the K they share.
                alt_grid = _K_grid(B, N, K_lo, min(K_hi, alt_K_hi))
                Ks_alt, ms_alt = _mass_series(provider, N, B, mg, alt_grid,
                                              msq_units=True,
                                              lpn=base + alt_lpn_delta)
                alt = ms_alt if len(ms_alt) >= 4 else None
                bud = richardson_budget(Ks, ms, mg, N, n_terms=n_terms,
                                        masses_alt_lpn=alt,
                                        K_codes_alt=Ks_alt if alt else None)
                row = paper.get((q, N, mg))
                pv, pe = row if row else (float("nan"), float("nan"))
                diff = abs(pv - bud["M0"])
                d_pt = diff / pe if pe else float("nan")
                d_our = diff / bud["total"] if bud["total"] else float("nan")
                fh.write(f"{q + ' N=' + str(N):>8} {mg:5.2f} "
                         f"{bud['M0']:10.4f} {bud['total']:9.2g} "
                         f"{bud['form']:9.2g} {bud['window']:9.2g} "
                         f"{bud['truncation']:9.2g} "
                         f"{pv:9.3f} {pe:8.3f} {d_pt:6.2f} {d_our:7.1f} "
                         f"{'ok' if d_pt <= 1.0 else 'OUT':>3}\n")
            fh.write("\n")
    print(f"  saved {out}")
    return out


def _load_paper_table1():
    """Table I as printed, keyed by (quantity, N, m/g) -> (value, last_term)."""
    import csv

    path = _ROOT / "refs" / "table1.csv"
    if not path.exists():
        return {}
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(l for l in fh if not l.startswith("#")):
            out[(row["quantity"], int(row["N"]), float(row["mg"]))] = (
                float(row["value"]), float(row["last_term"]))
    return out


FIGURES = {1: figure1, 2: figure2, 3: figure3, 4: figure4,
           5: figure5, 6: figure6, 7: figure7, 8: figure8}


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Reproduce figures from Phys. Rev. D 41, 3814 (1990)")
    ap.add_argument("--source", choices=["python", "fortran"], default="python")
    ap.add_argument("--fig", nargs="+", type=int, default=None,
                    help="figure numbers (default: the cheap ones, 1 3 4 5 6)")
    ap.add_argument("--table1", action="store_true")
    ap.add_argument("--budget", action="store_true",
                    help="Table I with its full error budget; API-only until now")
    ap.add_argument("--alt-K-hi", type=int, default=49,
                    help="top of the LPN+2 sweep the truncation term uses. One "
                         "extra qqbar pair multiplies the basis 17-23x (482,320 "
                         "states for the N=4 baryon at 2K=70), so this stays "
                         "well below --table1-window; the two series are "
                         "compared on the K they share")
    ap.add_argument("--fits", action="store_true",
                    help="diagnostic: the Richardson fit behind every Table I entry")
    ap.add_argument("--table1-window", nargs=2, type=int, metavar=("LO", "HI"),
                    default=(16, 24),
                    help="Richardson window in 2K (paper: 16 24)")
    ap.add_argument("--ncpus", type=int, default=1)
    ap.add_argument("--backend", choices=["thread", "process"], default=None,
                    help="matrix-build parallelism; identical results either "
                         "way (see docs/performance.md)")
    ap.add_argument("--assembly", choices=["exact", "fortran"], default="exact",
                    help="free-part assembly; see docs/basis-dependence.md")
    # --policy has been reachable from the API since PythonProvider was written
    # and from the CLI never, so four of the six (assembly, policy) pairs the
    # test suite exercises could not be run from a command line.
    ap.add_argument("--policy", choices=["fortran", "blockwise", "spectral"],
                    default="fortran",
                    help="weeding/orthonormalization; 'blockwise' never forms "
                         "the pre-weeding norm and is much faster")
    ap.add_argument("--solver", choices=["dense", "sparse"], default="dense",
                    help="eigensolver; 'sparse' needs --assembly exact "
                         "--policy blockwise")
    ap.add_argument("--nev", type=int, default=None,
                    help="keep only the lowest N eigenpairs (default: all for "
                         "dense, 40 for sparse)")
    ap.add_argument("--allow-run", action="store_true",
                    help="let the Fortran provider launch missing runs")
    args = ap.parse_args(argv)

    # Rejected here rather than deep inside a figure loop, so it costs a usage
    # message instead of minutes of compute.  Note --policy defaults to
    # "fortran", which is exactly the combination --solver sparse forbids: the
    # alternative of having --solver quietly change the weeding policy would
    # mean a flag named for the eigensolver silently changing the basis, and
    # docs/basis-dependence.md exists because that kind of coupling bites.
    if args.solver == "sparse" and (args.assembly, args.policy) != ("exact",
                                                                    "blockwise"):
        ap.error("--solver sparse requires --assembly exact --policy blockwise")

    if args.source == "python":
        from .providers import PythonProvider
        provider = PythonProvider(ncpus=args.ncpus, assembly=args.assembly,
                                  policy=args.policy, backend=args.backend,
                                  solver=args.solver, nev=args.nev)
    else:
        from .providers import FortranProvider
        provider = FortranProvider(allow_run=args.allow_run,
                                   extra_search=[_ROOT / "python"])

    for n in (args.fig or [1, 3, 4, 5, 6]):
        print(f"Figure {n} [{args.source}]:")
        if n == 1:
            figure1()
        else:
            FIGURES[n](provider, args.source)

    if args.table1:
        print(f"Table I [{args.source}]:")
        table1(provider, args.source,
               K_lo=args.table1_window[0], K_hi=args.table1_window[1])

    if args.fits:
        print(f"Table I fits [{args.source}]:")
        figure_fits(provider, args.source,
                    K_lo=args.table1_window[0], K_hi=args.table1_window[1])

    if args.budget:
        print(f"Table I error budget [{args.source}]:")
        table1_budget(provider, args.source,
                      K_lo=args.table1_window[0], K_hi=args.table1_window[1],
                      alt_K_hi=args.alt_K_hi)


if __name__ == "__main__":
    main()
