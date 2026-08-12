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

## External data appearing in the paper's figures

Figure 8(a) overlays SU(2) Hamiltonian-lattice results from:

> C. J. Hamer, *Nucl. Phys.* **B195**, 503 (1982).

The large-N reference curves in Fig. 8 come from 't Hooft's solution:

> G. 't Hooft, *Nucl. Phys.* **B75**, 461 (1974).
