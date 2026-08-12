"""Result providers: where a figure gets its :class:`DLCQResult` from.

A figure asks for "N=3, B=1, 2K=21, lambda=0.3325" and a provider returns it,
either by running the Python solver or by locating/launching a Fortran run.
The figure code never learns which.  That is what lets the same plotting
routine be driven from either solver and the outputs compared meaningfully.

Both providers cache to HDF5 so an expensive sweep is paid for once.
"""

from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path

from .dataset import DLCQResult, load, save

__all__ = ["Provider", "PythonProvider", "FortranProvider"]

_ROOT = Path(__file__).resolve().parent.parent


def _tag(N, NF, B, K_code, rlamb, cutoff, LPN, extra="") -> str:
    base = f"N{N}_NF{NF}_B{B}_K{K_code}_lam{rlamb:.10g}_cut{cutoff:g}_LPN{LPN}"
    if extra:
        base += "_" + hashlib.md5(extra.encode()).hexdigest()[:8]
    return base


class Provider:
    """Interface: return a :class:`DLCQResult` for one parameter set."""

    name = "abstract"

    def get(self, N, NF, B, K_code, rlamb, cutoff=-1.0, LPN=0) -> DLCQResult:
        raise NotImplementedError


class PythonProvider(Provider):
    """Run the Python solver, caching results to HDF5.

    ``assembly="exact"`` (the default) is the basis-independent scheme and the
    scientific reference; ``"fortran"`` reproduces the historical diagonal-only
    step.  See ``docs/basis-dependence.md``.
    """

    name = "python"

    def __init__(self, ncpus=1, assembly="exact", policy="fortran",
                 cache_dir=None):
        self.ncpus = ncpus
        self.assembly = assembly
        self.policy = policy
        self.cache_dir = Path(cache_dir or _ROOT / "runs" / "python_cache")

    def get(self, N, NF, B, K_code, rlamb, cutoff=-1.0, LPN=0) -> DLCQResult:
        from .read_python import run_python

        tag = _tag(N, NF, B, K_code, rlamb, cutoff, LPN,
                   extra=f"{self.assembly}:{self.policy}")
        path = self.cache_dir / f"{tag}.h5"
        if path.exists():
            return load(path)

        result = run_python(N=N, NF=NF, B=B, K_code=K_code, rlamb=rlamb,
                            cutoff=cutoff, LPN=LPN, ncpus=self.ncpus,
                            policy=self.policy, assembly=self.assembly)
        save(result, path)
        return result


class FortranProvider(Provider):
    """Locate an existing Fortran run, or launch one via ``fortran/run_case.sh``.

    Set ``allow_run=False`` to make missing runs an error instead of a several-
    minute compute (2K=25 takes ~5 minutes, 2K=29 substantially longer).
    """

    name = "fortran"

    # Tolerance on matching lambda to a stored run.  It has to thread a needle:
    #   * loose enough for the 8-digit lambda a sweep driver writes into the
    #     Fortran input (0.99609604 vs mg_to_lambda(0.05) = 0.9960959907,
    #     a 5e-8 difference that is physically irrelevant);
    #   * tight enough to keep 0.3325 and mg_to_lambda(1.6) = 0.33254949
    #     distinct -- they differ by 1.5e-5 and DO give different eigenvalues,
    #     which is the whole reason paper_lambda exists.
    # 1e-6 sits between them by two orders of magnitude on each side.
    LAMBDA_ATOL = 1e-6

    def __init__(self, run_root=None, allow_run=True, extra_search=(),
                 lam_atol=None):
        self.run_root = Path(run_root or _ROOT / "runs")
        self.allow_run = allow_run
        self.lam_atol = self.LAMBDA_ATOL if lam_atol is None else lam_atol
        # Historical outputs preserved in python/ carry K-tagged names.
        self.extra_search = [Path(p) for p in extra_search]

    # The Python solver writes a header in the same format as the Fortran, so
    # matching on parameters alone would silently accept qcdf_py.out as a
    # "Fortran" run.  These sections are emitted only by qcdf.f (PRNTST and
    # PRTEIG with IVECT>0) and never by the Python solver.
    _FORTRAN_MARKERS = ("TOTAL NUM STATES", "NUMBER OF STATES BEFORE WEEDING")

    @staticmethod
    def _is_fortran_output(text: str) -> bool:
        return all(m in text for m in FortranProvider._FORTRAN_MARKERS)

    def _find_existing(self, N, NF, B, K_code, rlamb, LPN=None) -> Path | None:
        from .read_fortran import _parse_header

        # rglob, not glob: sweep runs live at runs/sweep/<tag>/qcdf.out
        candidates = list(self.run_root.rglob("qcdf.out"))
        for extra in self.extra_search:
            # Exclude *_py.out explicitly as well as by content.
            candidates.extend(sorted(p for p in extra.glob("qcdf*.out")
                                     if not p.name.endswith("_py.out")))

        for path in candidates:
            try:
                text = path.read_text(errors="replace")
            except OSError:
                continue
            if not self._is_fortran_output(text):
                continue
            try:
                meta = _parse_header(text.splitlines())
            except Exception:
                continue
            if (meta["N"] == N and meta["NF"] == NF and meta["B"] == B
                    and meta["K_code"] == K_code
                    and abs(meta["rlamb"] - rlamb) < self.lam_atol
                    and (LPN is None or meta["LPN"] == LPN)):
                return path
        return None

    def get(self, N, NF, B, K_code, rlamb, cutoff=-1.0, LPN=0) -> DLCQResult:
        from .read_fortran import read_out

        path = self._find_existing(N, NF, B, K_code, rlamb, LPN)
        if path is None:
            if not self.allow_run:
                raise FileNotFoundError(
                    f"no Fortran run for N={N} NF={NF} B={B} 2K={K_code} "
                    f"lambda={rlamb}; run fortran/run_case.sh or pass allow_run=True"
                )
            tag = _tag(N, NF, B, K_code, rlamb, cutoff, LPN)
            outdir = self.run_root / tag
            script = _ROOT / "fortran" / "run_case.sh"
            subprocess.run(
                ["bash", str(script), str(outdir), str(N), str(NF), str(B),
                 f"{rlamb:.10g}", f"{cutoff:g}", str(LPN), str(K_code)],
                check=True, capture_output=True,
            )
            path = outdir / "qcdf.out"

        return read_out(path)
