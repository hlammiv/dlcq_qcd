"""Shared fixtures and markers for the DLCQ regression suite."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

FIXTURES = ROOT / "tests" / "fixtures"

# Fortran outputs are gitignored (they are large and regenerable), so a fresh
# checkout -- CI included -- has none until fortran/run_case.sh is run. Each
# entry lists the locations to try: the preserved historical files first, then
# where run_case.sh puts a fresh run.
HISTORICAL = {
    "K21_B1": [ROOT / "python" / "qcdf.out",
               ROOT / "runs" / "K21" / "qcdf.out",
               ROOT / "runs" / "K21_B1" / "qcdf.out"],
    "K25_B1": [ROOT / "python" / "qcdf_K25.out",
               ROOT / "runs" / "K25" / "qcdf.out"],
    "K29_B1": [ROOT / "python" / "qcdf_K29.out",
               ROOT / "runs" / "K29" / "qcdf.out"],
}


def _available(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def find_run(key: str) -> Path | None:
    """First existing output for ``key``, or None."""
    for path in HISTORICAL[key]:
        if _available(path):
            return path
    return None


@pytest.fixture(scope="session")
def fortran_K21():
    """The 2K=21, N=3, B=1 reference run -- the workhorse for these tests."""
    from dlcq.read_fortran import read_out

    path = find_run("K21_B1")
    if path is None:
        pytest.skip("no 2K=21 B=1 Fortran run; see fortran/run_case.sh")
    return read_out(path)


@pytest.fixture(scope="session")
def fortran_K21_nomat():
    """Same run, spectrum only -- much cheaper when matrices are not needed."""
    from dlcq.read_fortran import read_out

    path = find_run("K21_B1")
    if path is None:
        pytest.skip("no 2K=21 B=1 Fortran run; see fortran/run_case.sh")
    return read_out(path, with_matrices=False)


@pytest.fixture(scope="session")
def python_K21():
    """Python solver on the 2K=21 configuration, Fortran-compatible weeding.

    Uses the literal lambda from the Fortran run rather than re-deriving it
    from m/g -- 0.3325 vs mg_to_lambda(1.6) = 0.33254949 differ by enough to
    swamp any tight tolerance.
    """
    from dlcq.read_python import run_python

    return run_python(N=3, NF=1, B=1, K_code=21, rlamb=0.3325,
                      cutoff=-1.0, LPN=0, ncpus=4,
                      policy="fortran", assembly="exact")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: expensive; excluded by default")
    config.addinivalue_line("markers", "fortran: needs Fortran output or the binary")
