# DLCQ QCD — Python/Parallel Translation

## Overview

This is a Python translation of Kent Hornbostel's 1993 Fortran 77 program `qcdf.f`, which computes hadron mass spectra in 1+1 dimensional QCD using **Discretized Light-Cone Quantization (DLCQ)** with support for multiple quark flavors.

The translation uses:
- **NumPy** for array operations and matrix storage
- **SciPy** (`scipy.linalg.eigh`) for symmetric eigenvalue decomposition (replacing the Numerical Recipes routines TRR8/TQR8/ESRTR8)
- **numba** `nogil` kernels (`qcdf_kernels.py`) for the colour factor and Hamiltonian element, driven from a **thread pool** so the O(N²) matrix construction parallelizes without forking or pickling. A **multiprocessing** path over the interpreted reference routines is kept for cross-checking; see [../docs/performance.md](../docs/performance.md)

## Files

| File | Description |
|------|-------------|
| `qcdf.py` | Main Python program (~1800 lines), and the reference for state generation |
| `qcdf_opt.py` | Optimized driver: matrix builds, weeding, backend selection |
| `qcdf_kernels.py` | `nogil` numba kernels for `clfact`/`hamqcd` and the row workers |
| `qcdf_states.py` | `nogil` numba kernels for Fock-state generation (`brmsgn`, `prmx`, `prmm`) |
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
export QCDF_NCPUS=14          # threads; past the physical core count adds little
python3 qcdf.py input_original.json

export QCDF_BACKEND=process   # interpreted reference path, for cross-checking
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

The double loop over state pairs (subroutine `CLRDIS` in the original Fortran) used to be the dominant cost. It is parallelized one row at a time:

- `qcdf_kernels.py` compiles the colour factor and the Hamiltonian element with `@njit(nogil=True)`, so worker **threads** run genuinely concurrently
- each thread owns a `Scratch` — the module-global buffers in `qcdf_opt.py` were safe only because `fork` gave each process a copy
- each row is written straight into the destination matrix, so nothing is pickled and nothing is copied back; workers touch disjoint rows and fill only the upper triangle, so there is no lock
- rows are dispatched longest-first, since row *r* costs `numsta - r`

Set `QCDF_BACKEND=process` to run the interpreted reference routines under `multiprocessing` instead. The two produce **bit-identical** matrices — `tests/test_kernels.py` asserts it — so the choice is purely about speed.

Scaling is ~7× on 14 physical cores; past that the machine has only hyperthreads and E-cores left, so more threads add little. The eigenvalue decomposition uses SciPy's LAPACK routines, which are internally threaded, and is now a *larger* share of a run than the matrix build.

## Correspondence to Fortran Subroutines

| Fortran | Python |
|---------|--------|
| MAIN PROG | `main()` |
| QCDSTA | `qcdsta()`; `qcdf_opt.qcdsta_fast()` |
| BRMSGN | `brmsgn()` |
| PRMX | `prmx()` |
| PRMM | `prmm()` |
| GPRBIG | `gprbig_array()` |
| IBARFL | `ibarfl()` |
| IMESFL | `imesfl()` |
| CLRDIS | `clrdis()`; `qcdf_opt.build_matrices()` / `qcdf_kernels.norm_row()`, `ham_row()` |
| HAMQCD | `hamqcd()`; `qcdf_kernels.hamqcd_nb()` |
| CLFACT + CLSUM + CNTRCT + BREDCE + NWTERM | `clfact_compute()`; `qcdf_kernels.clfact_nb()` |
| EPSB (epsilon permutation table) | `qcdf_kernels.epsb_table()` |
| BRMSGN | `brmsgn()`; `qcdf_states.brmsgn_nb()` |
| PRMX / PRMM | `prmx()` / `prmm()`; `qcdf_states.prmx_nb()` / `prmm_nb()` |
| STACHK | `stachk()`; `qcdf_states.stachk_nb()` |
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
numba
```

## Notes

- The Fortran code uses extensive `COMMON` blocks for global state; the Python version uses dataclass containers passed explicitly
- Array indexing has been translated from 1-based (Fortran) to 0-based (Python) throughout
- The color contraction engine (`clfact_compute`) faithfully reproduces the diagrammatic delta-function contraction and epsilon-tensor reduction algorithms from the original
- The colour factor was the bottleneck for a long time, then state generation was; with both compiled, neither is. A full N=3, B=1 run at 2K=29 takes ~1.8 s, of which the Hamiltonian build is 0.28 s and state generation 0.01 s. The dense `eigh` and the norm orthonormalization now dominate
- The compiled generator sizes `kpx` and the state arrays to the run, so `LKPXMX`/`NMSTMX` bound only `qcdf.py`'s interpreted reference. It reaches 2K=48 (190,171 states) in 5.8 s
- What limits reach now is memory: the norm is built dense at the pre-weeding size, 8.6 GB at 2K=37
