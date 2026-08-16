#!/usr/bin/env python3
"""Frozen-number extractor: every number the paper quotes, from data, once.

Emits into ``paper/frozen/``:
  * ``vN.json``       — the frozen values with their inputs, one file per freeze
  * ``numbers.tex``   — LaTeX macros; prose uses ONLY these, never literals
  * ``provenance.json`` — git HEAD, input-file SHA256s, and the pinned blob
    SHA of docs/weak-coupling-limit.md that all quotes cite

Basis-pairing gate (G3): every entry carries its (hamiltonian, fit_basis)
pair, and this script refuses to emit a mismatch — improved pairs only with
plain 1/K (van de Sande Eq. 14), standard only with the Eq. (27) ladder.

Freeze discipline: v1 (pre-production, writing basis) and v2 (post-production,
publication basis).  Between freezes the macro file is read-only; a v1→v2
diff with a sign flip or a >1-bar migration is investigated, not absorbed.

Usage:
    python paper/make_numbers.py --freeze v1
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FROZEN = Path(__file__).resolve().parent / "frozen"

VALID_PAIRS = {("improved", "1/K"), ("standard", "eq27")}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def _git(*args) -> str:
    out = subprocess.run(["git", *args], cwd=ROOT, capture_output=True,
                         text=True)
    return out.stdout.strip()


def gather() -> dict:
    """Compute every frozen number from data/ + runs/python_cache/.

    Entries: name -> {value, err, hamiltonian, fit_basis, inputs, note}.
    Filled in per drafting phase; the gate below runs on whatever exists.
    """
    numbers: dict = {}
    # TODO(phase 2): Table I entries with full budget, improved Table I,
    # ratio table (converged 0.4% statement), GMOR intercept, alpha
    # exponents, captured fractions, timing/reach table.
    return numbers


def check_pairs(numbers: dict) -> None:
    for name, entry in numbers.items():
        ham = entry.get("hamiltonian")
        basis = entry.get("fit_basis")
        if ham is None and basis is None:
            continue  # not an extrapolated quantity
        if (ham, basis) not in VALID_PAIRS:
            raise SystemExit(
                f"G3 violation in {name!r}: hamiltonian={ham!r} with "
                f"fit_basis={basis!r}. improved pairs with '1/K', standard "
                f"with 'eq27'; nothing else is emittable.")


def emit(numbers: dict, tag: str) -> None:
    FROZEN.mkdir(exist_ok=True)
    (FROZEN / f"{tag}.json").write_text(json.dumps(numbers, indent=2))
    lines = [f"% frozen {tag} — {date.today().isoformat()} — generated file, "
             f"do not edit"]
    for name, entry in sorted(numbers.items()):
        macro = "".join(w.capitalize() for w in name.split("_"))
        lines.append(f"\\newcommand{{\\num{macro}}}{{{entry['value']}}}")
    (FROZEN / "numbers.tex").write_text("\n".join(lines) + "\n")
    prov = {
        "freeze": tag,
        "git_head": _git("rev-parse", "HEAD"),
        "weak_coupling_doc_blob": _git("rev-parse",
                                       "HEAD:docs/weak-coupling-limit.md"),
        "inputs": {str(p.relative_to(ROOT)): _sha256(p)
                   for entry in numbers.values()
                   for p in (ROOT / i for i in entry.get("inputs", []))
                   if p.exists()},
    }
    (FROZEN / "provenance.json").write_text(json.dumps(prov, indent=2))
    print(f"froze {len(numbers)} numbers as {tag} at {prov['git_head'][:8]}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--freeze", required=True, help="freeze tag, e.g. v1")
    args = ap.parse_args(argv)
    numbers = gather()
    check_pairs(numbers)
    emit(numbers, args.freeze)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
