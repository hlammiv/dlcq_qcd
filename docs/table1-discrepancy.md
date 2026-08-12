# Table I disagrees with the paper's own Fig. 8(a)

Our sweep reproduces the paper's **figures** but not its **Table I**. The
evidence says the figures are right.

## The independent arbiter

Fig. 8(a) overlays SU(2) Hamiltonian-lattice results from C. J. Hamer,
Nucl. Phys. **B195**, 503 (1982) — data computed by a different method, by
different authors, published separately. Digitizing those points from the
paper's own figure and comparing with our DLCQ sweep (Richardson-extrapolated
over 2K = 16–24):

| m/g | Hamer lattice | our DLCQ | diff |
|---|---|---|---|
| 0.205 | 0.6590 | 0.6561 | 0.4% |
| 0.409 | 1.1109 | 1.1348 | 2.2% |
| 0.823 | 1.9749 | 1.9826 | 0.4% |
| 1.647 | 3.5187 | 3.5236 | **0.1%** |

(The two smallest-m/g points differ more; that is the strong-coupling regime
where our LPN = 4 Fock truncation bites hardest and where the paper itself says
its own error estimates are "not more than a rough guide".)

## The disagreement

At m/g = 1.6, SU(2) meson `M/g`:

| source | value |
|---|---|
| Hamer lattice, read off the paper's Fig. 8(a) | 3.519 |
| our DLCQ | 3.524 |
| the paper's **Table I** | **4.314** |

**Fig. 8(a)'s y-axis maximum is 4.0**, verified from its tick structure (9 ticks
every 0.5). Table I's 4.314 cannot be plotted on the figure it supposedly
summarizes. The same holds across the table: every N = 2, 3, 4 meson entry at
m/g = 1.6 (4.314, 4.618, 4.845) exceeds that axis.

Two further symptoms:

* **Table I saturates.** The N=3 meson goes 3.1 → 4.40 → 4.618 for
  m/g = 0.4 → 0.8 → 1.6, i.e. nearly flat once m/g > 0.8. `M/g` must grow like
  `2m/g` in the free limit, so saturation is unphysical.
* **The ratio to our values is not constant** — it runs 0.55 to 4.8 across the
  table — so this is not a single units or normalization factor.

## What we did about it

Nothing is "corrected". `refs/table1.csv` transcribes Table I exactly as
printed (verified against a 600 dpi render, all 35 entries). The comparison is
reported, not reconciled, because the resolution is a physics-judgement call
for the user rather than something to infer.

What the repository asserts instead:

* the paper's **exact** results, which hold identically (N=2 meson/baryon
  degeneracy is reproduced to 2×10⁻⁹ relative at every m/g and K);
* the paper's **figures**, reproduced by both solvers at 1–3%;
* **Hamer's independent lattice data**, matched to 0.1–2.2%.

Table I comparisons are marked `xfail` pending a resolution.
