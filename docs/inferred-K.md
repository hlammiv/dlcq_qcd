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

## Fig. 6(d): K from the row of zero-valued markers

The caption of Fig. 6 states 2K = 21 for panels (a)–(c) and **no K at all for
(d)**. Inferring it from the curve markers, as everywhere else on this page,
gives the wrong answer: those markers are few, they sit where three curves
overlap, and `infer_K_from_x_grid` returns the *smallest* K explaining a quorum,
so their scatter biases it low. It answers 21 for the article's panel and 22 for
the thesis's reprint of the same plot. Both are wrong; it is **24**.

There is a much better signal on exactly this kind of panel. Wherever a
structure function vanishes, the marker is still drawn and comes to rest on the
axis — and a two-baryon distribution is zero over most of its range, so the
panel carries a long row of zero-valued markers along the bottom, one per
lattice site, at full contrast with nothing to merge with. Measured against the
frame:

| site | measured | k/24 | k/22 |
|---|---|---|---|
| 1 | 0.0427 | 0.0417 | 0.0455 |
| 3 | 0.1304 | 0.1250 | 0.1364 |
| 5 | 0.2095 | 0.2083 | 0.2273 |
| 7 | 0.2937 | 0.2917 | 0.3182 |
| 9 | 0.3763 | 0.3750 | 0.4091 |
| 11 | 0.4569 | 0.4583 | 0.5000 |
| 13 | 0.5427 | 0.5417 | 0.5909 |

Worst deviation: 0.006 for K = 24, 0.048 for K = 22 — a full site's drift by the
right-hand edge, far outside the scatter.

Two facts make the fit well posed. Momenta are odd, so the lattices for K and
2K are **disjoint** rather than nested and a factor-of-two error cannot hide.
And K's parity is fixed outright: a state of n partons, each of odd momentum,
has K_code ≡ n (mod 2) — even for mesons and for the B = 2 baryon pair, odd for
a single N = 3 baryon. That alone rules out the K = 23 the thesis panel's
noisier dot row otherwise prefers.

`tools/digitize.py:infer_K_from_axis_dots` implements this. It is deliberately
gated to speak only when it finds ≥ 7 dots with ≥ 80% matched, because the
method needs the distribution to actually vanish somewhere: on a valence meson
panel, nonzero nearly everywhere, there is no dot row and the band picks up tick
marks instead — Fig. 3(a) yields a confident-looking 5-of-6 for K = 18, which is
wrong. Under that gate it reports on Figs. 6(b) and 6(d) only, and is correct on
both.

