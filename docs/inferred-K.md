# Recovering the unstated K of Figs. 3 and 4

The paper gives K for Figs. 5 (2K = 24) and 6 (2K = 21) but never for Figs. 3
and 4. It can be read off the published plots.

## Why the plot knows

Antiperiodic boundary conditions force half-odd-integer momenta, so in code
units every parton carries an **odd integer** k and a structure function is
sampled only at

```
x = k / K ,  k = 1, 3, 5, ...     =>   x_min = 1/K ,  Δx = 2/K
```

The marker positions therefore *are* K. `dlcq.units.infer_K_from_x_grid` fits
that lattice: for each candidate K it scores how many traced points fall within
tolerance of a site, and takes the smallest K explaining ≥60% of them. Smallest
matters — a finer lattice always contains a coarser one.

A second, independent handle is the vertical scale. Eq. (12) carries a factor
`K_paper = K/2`, and the grid spacing is `1/K_paper`, so
`Σᵢ q(xᵢ)/K_paper` is the quark number. `infer_K_from_normalization` uses it.

## Validation before use

Run the procedure on the panels whose K the paper *does* state:

| panel | paper | recovered |
|---|---|---|
| fig5a, fig5b, fig5c | 2K = 24 | **24** |
| fig6a, fig6b, fig6c | 2K = 21 | **21** |
| fig6d (B=2) | not stated | 22 |

## Result

| panel | sector | recovered 2K | required parity | agrees? |
|---|---|---|---|---|
| fig3a | meson valence | **14** | even | yes |
| fig3b | baryon valence | **15** | odd | yes |
| fig4a | meson, one extra pair | **14** | even | yes |
| fig4b | baryon, one extra pair | **15** | odd | yes |
| fig4c | baryon, two extra pairs | **15** | odd | yes |

The parity check is the strongest part of this. A state's momenta are odd
integers summing to K, so an even parton count forces K even and an odd count
forces K odd: mesons need even K, N=3 baryons odd. **Nothing in the lattice fit
knows this** — yet every meson panel returned an even K and every baryon panel
an odd one. Five independent coin flips landing correctly.

Confirmed against the physics: at the recovered K, both solvers reproduce the
published curves to 2.6% (Fig. 3a) and 1.9% (Fig. 3b).

`dlcq.figures` uses 2K = 14 (meson) and 15 (baryon) for Figs. 3 and 4.
