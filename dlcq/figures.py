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
                          physical_indices, spurious_zero_modes)
from .units import (code_to_M_over_g, lambda_to_mg, mg_to_lambda,
                    thooft_rescale)

__all__ = ["figure1", "figure2", "figure3", "figure4", "figure5",
           "figure6", "figure7", "figure8", "table1", "FIGURES"]

_ROOT = Path(__file__).resolve().parent.parent
FIGDIR = _ROOT / "figures"

# Paper parameters.
MG_STRONG, MG_WEAK = 0.1, 1.6
K_MESON, K_BARYON = 24, 21

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
        lambdas = np.concatenate([np.arange(0.05, 0.50, 0.05),
                                  np.arange(0.50, 0.95, 0.05),
                                  [0.95, 0.97, 0.99]])
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
        for ax, B, K in [(ax1, 0, K_MESON), (ax2, 1, K_BARYON)]:
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

    ax1.set_title(f"(a) SU(3) meson,  2K={K_MESON}")
    ax2.set_title(f"(b) SU(3) baryon,  2K={K_BARYON}")
    for ax in (ax1, ax2):
        ax.set_xlabel("x = k/K"); ax.set_ylabel("q(x)")
        ax.set_xlim(0, 1); ax.legend(fontsize=8)

    fig.suptitle(f"FIG. 3  Valence structure functions  [{source}]"
                 "   (K not stated in paper; adopted from Figs. 5-6)")
    return _finish(fig, "fig3_valence", source)


def figure4(provider, source):
    """Fig. 4 -- higher-Fock contributions."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    specs = [
        (axes[0], 0, K_MESON, 4, "(a) meson: $q\\bar q q\\bar q$", False),
        (axes[1], 1, K_BARYON, 5, "(b) baryon: $qqqq\\bar q$", True),
        (axes[2], 1, K_BARYON, 7, "(c) baryon: two extra pairs", False),
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
        x, q, qbar = structure_function(r, idx, nparton=npart)
        if np.max(np.abs(q)) > 1e-18:
            ax.plot(x, q, style, markersize=3, lw=1, label=lab)
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
        r2 = provider.get(3, 1, 2, 22, lam)
        _wavefunction_panel(ax, r2, 0, [(6, "k-o", "6q"),
                                        (8, "r-s", r"$6qq\bar q$")],
                            "(d) 1st B=2, 2K=22")
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


def _extrapolated_mass(provider, N, B, mg, K_codes):
    """Lightest M/g extrapolated to the continuum via Eq. (27)."""
    if mg == 0.0:
        # Exact: the lightest state is massless for any N, B, K (Eq. 16).
        return 0.0, 0.0
    lam = paper_lambda(mg)
    Ks, masses = [], []
    for K in K_codes:
        r = provider.get(N, 1, B, K, lam)
        if r.n_eigenvalues == 0:
            continue
        phys = physical_indices(r)
        if phys.size == 0:
            continue
        Ks.append(K)
        masses.append(code_to_M_over_g(r.eigenvalues[phys[0]], lam))
    if len(Ks) < 3:
        return None, None
    return richardson_extrapolate(Ks, masses, mg, N)


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

    ax1.set_xlabel("m/g"); ax1.set_ylabel("M/g")
    ax1.set_title("(a) meson mass"); ax1.legend()
    ax2.set_xlabel(r"$(2\pi/N)^{1/2}\, m/g$")
    ax2.set_ylabel(r"$(2\pi/N)^{1/2}\, M/g$")
    ax2.set_title("(b) 't Hooft rescaling"); ax2.legend()

    fig.suptitle(f"FIG. 8  Comparison with large N  [{source}]"
                 "   (note: this paper's unit is $g^2N/2\\pi$)")
    return _finish(fig, "fig8_large_N", source)


def table1(provider, source, N_values=(2, 3, 4)):
    """Table I -- N dependence of meson and baryon mass, with Richardson errors."""
    rows = []
    for mg in MG_TABLE:
        row = {"mg": mg}
        for N in N_values:
            M0, err = _extrapolated_mass(provider, N, 0, mg, _K_grid(0, N))
            row[f"mes_N{N}"] = (M0, err)
        for N in (3, 4):
            M0, err = _extrapolated_mass(provider, N, 1, mg, _K_grid(1, N))
            row[f"bar_N{N}"] = (M0, err)
        rows.append(row)

    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = FIGDIR / f"table1_{source}.txt"
    with open(out, "w") as fh:
        fh.write("TABLE I.  N dependence of meson and baryon mass.\n")
        fh.write(f"Reproduced with the {source} solver.\n")
        fh.write("Parenthesized value = magnitude of the last Richardson term,\n")
        fh.write("matching the paper's convention (not a statistical error).\n\n")
        hdr = f"{'m/g':>6}" + "".join(f"{'mes N='+str(N):>16}" for N in N_values)
        hdr += "".join(f"{'bar N='+str(N):>16}" for N in (3, 4))
        fh.write(hdr + "\n" + "-" * len(hdr) + "\n")
        for row in rows:
            line = f"{row['mg']:6.2f}"
            for key in [f"mes_N{N}" for N in N_values] + [f"bar_N{N}" for N in (3, 4)]:
                M0, err = row.get(key, (None, None))
                line += f"{'---':>16}" if M0 is None else f"{M0:>11.3f}({err:.2g})".rjust(16)
            fh.write(line + "\n")
    print(f"  saved {out}")
    return rows


FIGURES = {1: figure1, 2: figure2, 3: figure3, 4: figure4,
           5: figure5, 6: figure6, 7: figure7, 8: figure8}


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Reproduce figures from Phys. Rev. D 41, 3814 (1990)")
    ap.add_argument("--source", choices=["python", "fortran"], default="python")
    ap.add_argument("--fig", nargs="+", type=int, default=None,
                    help="figure numbers (default: the cheap ones, 1 3 4 5 6)")
    ap.add_argument("--table1", action="store_true")
    ap.add_argument("--ncpus", type=int, default=1)
    ap.add_argument("--assembly", choices=["exact", "fortran"], default="exact",
                    help="free-part assembly; see docs/basis-dependence.md")
    ap.add_argument("--allow-run", action="store_true",
                    help="let the Fortran provider launch missing runs")
    args = ap.parse_args(argv)

    if args.source == "python":
        from .providers import PythonProvider
        provider = PythonProvider(ncpus=args.ncpus, assembly=args.assembly)
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
        table1(provider, args.source)


if __name__ == "__main__":
    main()
