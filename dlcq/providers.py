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

from .dataset import DLCQResult, Spectrum, load, read_spectrum, save

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

    ``backend`` picks how the matrix builds are parallelised.  It is absent
    from the cache tag on purpose: the two backends produce bit-identical
    matrices, so a cached run stays valid whichever built it.
    """

    name = "python"

    def __init__(self, ncpus=1, assembly="exact", policy="fortran",
                 cache_dir=None, backend=None, solver="dense", nev=None):
        self.ncpus = ncpus
        self.assembly = assembly
        self.policy = policy
        self.backend = backend
        self.solver = solver
        self.nev = nev
        self.cache_dir = Path(cache_dir or _ROOT / "runs" / "python_cache")

    def get(self, N, NF, B, K_code, rlamb, cutoff=-1.0, LPN=0) -> DLCQResult:
        from .read_python import run_python

        # Extend the tag only for non-default values.  ``solver`` must be in it
        # -- unlike ``backend``, an iterative solve is not bit-identical, so
        # caching across it would silently mix results -- but appending
        # ":dense" unconditionally would change the hash for every existing
        # file: md5("exact:fortran")[:8] is "31275504", the suffix on the great
        # majority of runs/python_cache, several GB representing hours of
        # compute.  ``nev`` likewise, since a 40-level result is a different
        # object from a full one.  And no "a cached k=100 satisfies a k=40
        # request" rule: the lowest 40 of a 100-vector Lanczos are not the same
        # bits as a 40-vector one, which is the same argument that puts
        # ``solver`` in the tag at all.
        extra = f"{self.assembly}:{self.policy}"
        if self.solver != "dense":
            extra += f":{self.solver}"
        if self.nev is not None:
            extra += f":k{self.nev}"
        tag = _tag(N, NF, B, K_code, rlamb, cutoff, LPN, extra=extra)
        path = self.cache_dir / f"{tag}.h5"
        if path.exists():
            return load(path)

        result = run_python(N=N, NF=NF, B=B, K_code=K_code, rlamb=rlamb,
                            cutoff=cutoff, LPN=LPN, ncpus=self.ncpus,
                            policy=self.policy, assembly=self.assembly,
                            backend=self.backend, solver=self.solver,
                            nev=self.nev)
        save(result, path)
        return result

    def spectrum(self, N, NF, B, K_code, rlamb, cutoff=-1.0, LPN=0):
        """The levels only, reading ~700 MB less per cached point than ``get``.

        A cached result stores ``Z`` densely even though it is 0.014% non-zero,
        so ``get`` inflates an 8 MB file to ~700 MB to hand back a 48-element
        spectrum.  Every Richardson consumer -- the sweeps, ``table1``,
        ``plot_holdout`` -- wants only the levels.

        Falls back to a full solve on a cache miss, so this is a read
        optimisation and never changes what gets computed.
        """
        extra = f"{self.assembly}:{self.policy}"
        if self.solver != "dense":
            extra += f":{self.solver}"
        if self.nev is not None:
            extra += f":k{self.nev}"
        path = self.cache_dir / (
            _tag(N, NF, B, K_code, rlamb, cutoff, LPN, extra=extra) + ".h5")
        if path.exists():
            return read_spectrum(path)
        r = self.get(N, NF, B, K_code, rlamb, cutoff, LPN)
        return Spectrum(eigenvalues=r.eigenvalues, rlamb=r.rlamb,
                        n_orth=getattr(r, "n_orth", 0))


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
        # Parsed headers, keyed by (path, mtime_ns, size).  Without this every
        # get() re-read every qcdf.out in full -- some are 12 MB -- just to
        # recover a header.  With ~160 run directories that turned a sweep into
        # tens of GB of I/O: the first full test run on a 16-core box took 50
        # minutes, the second (warm page cache) two.  The value is the parsed
        # metadata, or None for "not a Fortran output".
        self._header_cache: dict = {}

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
                st = path.stat()
            except OSError:
                continue
            key = (str(path), st.st_mtime_ns, st.st_size)
            if key in self._header_cache:
                meta = self._header_cache[key]
            else:
                meta = None
                try:
                    text = path.read_text(errors="replace")
                except OSError:
                    text = None
                # The whole file has to be read: the discriminating marker,
                # TOTAL NUM STATES, is emitted by PRNTST after the state list,
                # so its offset scales with the basis size.  The other marker
                # is NOT sufficient on its own -- qcdf.py:1956 writes
                # "NUMBER OF STATES BEFORE WEEDING" too.  Reading once per
                # file and caching is what makes that affordable.
                if text is not None and self._is_fortran_output(text):
                    try:
                        meta = _parse_header(text.splitlines())
                    except Exception:
                        meta = None
                self._header_cache[key] = meta
            if meta is None:
                continue
            if (meta["N"] == N and meta["NF"] == NF and meta["B"] == B
                    and meta["K_code"] == K_code
                    and abs(meta["rlamb"] - rlamb) < self.lam_atol
                    and (LPN is None or meta["LPN"] == LPN)):
                return path
        return None

    def get(self, N, NF, B, K_code, rlamb, cutoff=-1.0, LPN=0,
            rmq=None, iflv=None) -> DLCQResult:
        """``rmq``/``iflv`` are only meaningful for NF > 1.

        ``qcdf.f`` fixes ``rmq(1) = 1`` and reads ``rmq(2..NF)`` as ratios to
        it, then one flavour quantum number per flavour.  Note that it then
        **overwrites** ``iflv(1)`` with ``N*B`` regardless of what was entered
        (qcdf.f:253, flagged there with the author's own ``??????``), so any
        channel needing a different flavour-1 number is unreachable from the
        Fortran and generates zero states rather than an error.  The Python
        solver honours the value it is given.  See docs/flavour.md.
        """
        from .read_fortran import read_out

        path = self._find_existing(N, NF, B, K_code, rlamb, LPN)
        if path is None:
            if not self.allow_run:
                raise FileNotFoundError(
                    f"no Fortran run for N={N} NF={NF} B={B} 2K={K_code} "
                    f"lambda={rlamb}; run fortran/run_case.sh or pass allow_run=True"
                )
            extra = []
            if NF > 1:
                ratios = list(rmq[1:NF]) if rmq is not None else [1.0] * (NF - 1)
                numbers = list(iflv[:NF]) if iflv is not None else \
                    [N * B] + [0] * (NF - 1)
                extra = [f"{v:g}" for v in ratios] + [str(int(v)) for v in numbers]
            tag = _tag(N, NF, B, K_code, rlamb, cutoff, LPN,
                       extra="_".join(extra))
            outdir = self.run_root / tag
            script = _ROOT / "fortran" / "run_case.sh"
            subprocess.run(
                ["bash", str(script), str(outdir), str(N), str(NF), str(B),
                 f"{rlamb:.10g}", f"{cutoff:g}", str(LPN), str(K_code)] + extra,
                check=True, capture_output=True,
            )
            path = outdir / "qcdf.out"

        return read_out(path)
