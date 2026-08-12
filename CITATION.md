# Citation

## The paper this repository reproduces

> K. Hornbostel, S. J. Brodsky, and H. C. Pauli,
> *"Light-cone-quantized QCD in 1+1 dimensions"*,
> Phys. Rev. D **41**, 3814–3821 (1990).
> DOI: [10.1103/PhysRevD.41.3814](https://doi.org/10.1103/PhysRevD.41.3814)

```bibtex
@article{Hornbostel:1988fb,
    author    = "Hornbostel, K. and Brodsky, S. J. and Pauli, H. C.",
    title     = "{Light-cone quantized QCD in (1+1)-dimensions}",
    journal   = "Phys. Rev. D",
    volume    = "41",
    pages     = "3814--3821",
    year      = "1990",
    doi       = "10.1103/PhysRevD.41.3814"
}
```

**The article PDF is not redistributed in this repository.** It is copyrighted by the
American Physical Society. `literature/*.pdf` is listed in `.gitignore`. To run the
digitization tooling in `tools/`, place your own copy at
`literature/PhysRevD.41.3814.pdf`.

## The original Fortran code

`fortran/qcdf.f` carries the notice:

> `(C) Kent Hornbostel 1993.  All Rights Reserved.`
> `The routines ESRTR8 TQR8 TRR8 are property of the authors of Numerical Recipes.`

It is included here as the historical artifact being validated. The Numerical Recipes
routines `TRR8` / `TQR8` / `ESRTR8` (`tred2` / `tqli` / eigenvalue sort) are replaced by
`scipy.linalg.eigh` in the Python port.

## The thesis behind the code

Ref. 14 of the paper is the thesis this Fortran was written for:

> K. Hornbostel, *"The application of light-cone quantization to quantum
> chromodynamics in one-plus-one dimensions"*, Ph.D. thesis, Stanford
> University; SLAC Report No. 333 (1988).

It is a DOE technical report, freely available from OSTI:
<https://www.osti.gov/biblio/6783753> (full text:
<https://www.osti.gov/servlets/purl/6783753>).

It matters here for two reasons. Its Sec. 2.4 gives the structure-function
definition the code implements, and it reprints the same figures as the journal
article at noticeably better print quality — its Figs. 11 and 12 are the
article's Figs. 5 and 6. `tools/digitize.py` therefore traces the **thesis**
panels in preference to the journal scan.

Place a copy at `literature/SLAC-333_Hornbostel_thesis.pdf` to run that
tooling. Like the article, it is gitignored rather than redistributed here.

## External data appearing in the paper's figures

Figure 8(a) overlays SU(2) Hamiltonian-lattice results from:

> C. J. Hamer, *Nucl. Phys.* **B195**, 503 (1982).

The large-N reference curves in Fig. 8 come from 't Hooft's solution:

> G. 't Hooft, *Nucl. Phys.* **B75**, 461 (1974).
