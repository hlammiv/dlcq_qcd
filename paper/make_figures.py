#!/usr/bin/env python3
"""Single figure driver for the paper. Reads, never solves.

Contract (see the plan's "Repo mechanics"):
  * Inputs are exclusively ``data/**``, ``runs/python_cache/**`` (through
    ``PythonProvider(cache_dir=...)``), ``refs/digitized/**`` and
    ``refs/*.csv``.  This driver never calls the solver: before any provider
    call it verifies the cache file exists for every required key, and on a
    miss it exits listing the missing keys as ``paper/runs/`` commands to
    execute.  A paper build is therefore deterministic and fast, and missing
    data surfaces as a to-do list rather than a surprise multi-hour solve.
  * Standard-Hamiltonian series come from the cache; improved series come
    from ``data/*.csv`` — which enforces the never-cross-bases rule at the
    I/O layer (the cache holds no improved runs, by construction).
  * Peer-owned tools (plot_tiny_mg_fits.py, meson_baryon_ratio.py,
    gmor_chiral_limit.py) are invoked as subprocesses with ``--*-dir``
    arguments, never imported-and-monkeypatched, never edited.
  * Every output lands in ``paper/figs/`` (NOT ``figures/`` — that name is
    a gitignore glob at any depth) and embeds provenance (git SHA + inputs)
    in the PDF metadata.

Usage:
    python paper/make_figures.py [--only fig2 fig5 ratio ...] [--check] [--list]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIGS = Path(__file__).resolve().parent / "figs"
sys.path.insert(0, str(ROOT))

# ── registry ───────────────────────────────────────────────────────────────
# name -> (builder, [required cache keys or data files])
# Builders are added phase by phase; --list shows what exists so far.
REGISTRY: dict = {}


def register(name):
    def deco(fn):
        REGISTRY[name] = fn
        return fn
    return deco


def _git_sha() -> str:
    out = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                         cwd=ROOT, capture_output=True, text=True)
    return out.stdout.strip() or "unknown"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", nargs="*", default=None,
                    help="build only these registered figures")
    ap.add_argument("--check", action="store_true",
                    help="run caption/staleness checks instead of building")
    ap.add_argument("--list", action="store_true",
                    help="list registered figures and exit")
    args = ap.parse_args(argv)

    if args.list or not REGISTRY:
        print(f"{len(REGISTRY)} registered figure(s):")
        for name in sorted(REGISTRY):
            print(f"  {name}")
        if not REGISTRY:
            print("  (registry fills in as drafting phases land)")
        return 0

    names = args.only or sorted(REGISTRY)
    unknown = [n for n in names if n not in REGISTRY]
    if unknown:
        ap.error(f"unknown figure(s): {unknown}; use --list")

    FIGS.mkdir(exist_ok=True)
    for name in names:
        print(f"[{_git_sha()}] building {name}")
        REGISTRY[name]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
