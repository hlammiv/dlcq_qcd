# DLCQ QCD — Python/Parallel Translation

## Overview

This is a Python translation of Kent Hornbostel's 1993 Fortran 77 program `qcdf.f`, which computes hadron mass spectra in 1+1 dimensional QCD using **Discretized Light-Cone Quantization (DLCQ)** with support for multiple quark flavors.

The translation uses:
- **NumPy** for array operations and matrix storage
- **SciPy** (`scipy.linalg.eigh`) for symmetric eigenvalue decomposition (replacing the Numerical Recipes routines TRR8/TQR8/ESRTR8)
- **multiprocessing** to parallelize the O(N²) Hamiltonian matrix construction across CPU cores

## Files

| File | Description |
|------|-------------|
| `qcdf.py` | Main Python program (~1800 lines) |
| `input_small.json` | Small test case (K=6, LPN=4, ~5 states) |
| `input_medium.json` | Medium test case (K=10, ~21 states) |
| `input_original.json` | Original Fortran test case (K=24, ~818 states after weeding) |

## Usage

### Interactive mode
```bash
python3 qcdf.py
```
Prompts for parameters just like the original Fortran program.

### Batch mode
```bash
python3 qcdf.py input_small.json
```

### Control parallelism
```bash
export QCDF_NCPUS=8
python3 qcdf.py input_original.json
```

## JSON Input Format

```json
{
    "ihfile": 0,         // 0 = compute from scratch, 1 = read saved Hamiltonian
    "N": 3,              // number of colors
    "NF": 1,             // number of flavors
    "B": 0,              // baryon number
    "rlamb": 0.3325,     // coupling λ = sqrt(1/(1 + π m²/g²))
    "K": 24,             // 2× total light-cone momentum (discretization parameter)
    "cutoff": -1.0,      // Pauli-Villars cutoff (≤0 = none)
    "LPN": 0,            // particle number limit (0 = none)
    "rmq": [1.0],        // quark mass ratios (relative to first flavor)
    "iflv": [0]          // flavor quantum numbers for states
}
```

## Output

- `qcdf_py.out` — eigenvalues (mass²) and eigenvectors
- `qcdf_py.ham` — saved Hamiltonian for re-use with different λ values

## Physics

The program solves the light-cone Hamiltonian eigenvalue problem:

    P⁻ P⁺ |ψ⟩ = M² |ψ⟩

where P⁻ is the light-cone Hamiltonian for 1+1D QCD with SU(N) gauge group, discretized on a periodic spatial box with anti-periodic boundary conditions for fermions. States are Fock-space configurations of quarks and antiquarks, projected onto color singlets.

Key features:
- Anti-periodic boundary conditions (momenta are odd integers)
- Diagrammatic color-singlet projection (no explicit color indices stored)
- Arbitrary number of colors N and flavors NF
- Optional Pauli-Villars cutoff on Fock states
- Particle number truncation

## Parallelization Details

The dominant computational cost is the double loop over state pairs to compute Hamiltonian matrix elements (subroutine `CLRDIS` in the original Fortran). In Python this is parallelized with `multiprocessing.Pool.map()`:

- Each (i, j) matrix element is computed independently
- State data is passed as flat arrays to worker processes
- Results are assembled into the full symmetric matrix

For the original test case (818 states → ~335k matrix elements), this scales well to 4-8 cores. The eigenvalue decomposition itself uses SciPy's highly optimized LAPACK routines, which are internally threaded.

## Correspondence to Fortran Subroutines

| Fortran | Python |
|---------|--------|
| MAIN PROG | `main()` |
| QCDSTA | `qcdsta()` |
| BRMSGN | `brmsgn()` |
| PRMX | `prmx()` |
| PRMM | `prmm()` |
| GPRBIG | `gprbig_array()` |
| IBARFL | `ibarfl()` |
| IMESFL | `imesfl()` |
| CLRDIS | `clrdis()` (parallelized) |
| HAMQCD | `hamqcd()` |
| CLFACT + CLSUM + CNTRCT + BREDCE + NWTERM | `clfact_compute()` |
| WEEDR / WEEDR2 / DROPR | `weedr()` / `weedr2()` / `_dropr()` |
| NUHAM | `nuham()` |
| NUZ | `nuz()` |
| DIAG (TRR8/TQR8/ESRTR8) | `scipy.linalg.eigh()` |
| STACHK | `stachk()` |
| PAVILL | `pavill()` |
| FLVCHK | `flvchk()` |
| ODDCHK | `oddchk()` |
| LPNSUB | `lpnsub()` |
| SLFN | `compute_selfen()` |
| BRACK | `brack()` |
| CONJ | `conj_state()` |
| PSIGN | `psign()` |

## Requirements

```
numpy
scipy
```

## Notes

- The Fortran code uses extensive `COMMON` blocks for global state; the Python version uses dataclass containers passed explicitly
- Array indexing has been translated from 1-based (Fortran) to 0-based (Python) throughout
- The color contraction engine (`clfact_compute`) faithfully reproduces the diagrammatic delta-function contraction and epsilon-tensor reduction algorithms from the original
- For large K (K≥24), the computation may take significant time; the parallelization helps but the color factor computation remains the bottleneck
