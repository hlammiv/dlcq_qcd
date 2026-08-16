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


BUDGET_TXT = ROOT / "figures" / "table1_budget_2K25-71_python.txt"
IMPROVED_CSV = ROOT / "data" / "improved_table1_msq.csv"


def improved_fits() -> dict:
    """Improved-Hamiltonian Table I extrapolates, plain-1/K ensemble.

    Reuses the exact fit machinery of tools/plot_improved_fits.py (orders
    2--4 x contiguous sub-windows; median and 68% half-width), so the frozen
    numbers and the figure can never disagree.  Basis is plain 1/K per van
    de Sande Eq. (14) — enforced by the G3 pair on every entry.
    """
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "tools"))
    from plot_improved_fits import _read, _spread_1k

    data = _read(IMPROVED_CSV)
    out = {}
    for (N, B, mg), series in sorted(data.items()):
        ks = sorted(series)
        ys = [series[k] for k in ks]
        med, err, _ = _spread_1k(ks, ys)
        out[("mes" if B == 0 else "bar", N, mg)] = (med, err)
    return out


def parse_budget(path: Path = BUDGET_TXT) -> list[dict]:
    """Parse the table1_budget artifact into row dicts.

    Columns: case, m/g, ours, tot, form, wind, trunc, paper, pterm,
    d/pt, d/ours, flag.  Units are M^2/(m^2+g^2/pi) (the table's native
    units, see docs/table1-units.md); conversion to M/g happens where a
    figure or sentence needs it, never here.
    """
    rows = []
    for line in path.read_text().splitlines():
        parts = line.split()
        # data rows look like: mes N=2  1.60  4.3139  0.0019 ... 0.02 0.0 ok|OUT
        if len(parts) == 13 and parts[0] in ("mes", "bar") and "=" in parts[1]:
            rows.append({
                "channel": parts[0],
                "N": int(parts[1].split("=")[1]),
                "mg": float(parts[2]),
                "ours": float(parts[3]),
                "tot": float(parts[4]),
                "form": float(parts[5]),
                "wind": float(parts[6]),
                "trunc": float(parts[7]),
                "paper": float(parts[8]),
                "pterm": float(parts[9]),
                "d_pt": float(parts[10]),
                "d_ours": float(parts[11]),
                "flag": parts[12],
            })
    if len(rows) != 30:
        raise SystemExit(
            f"parse_budget: expected 30 Table I entries, got {len(rows)} "
            f"from {path} — regenerate the artifact before freezing.")
    return rows


def gather() -> dict:
    """Compute every frozen number from data/ + regenerated artifacts.

    Entries: name -> {value, err, hamiltonian, fit_basis, inputs, note}.
    Filled in per drafting phase; the gate below runs on whatever exists.
    """
    numbers: dict = {}
    budget_rel = str(BUDGET_TXT.relative_to(ROOT))
    rows = parse_budget()

    # Headline scalars quoted in prose (standard Hamiltonian, Eq. 27 ladder).
    weak = [r for r in rows if r["mg"] <= 0.1]
    worst = max(weak, key=lambda r: r["tot"] / r["ours"])
    numbers["tableone_n_entries"] = {
        "value": len(rows), "inputs": [budget_rel],
        "note": "entries in the reproduced Table I"}
    numbers["tableone_worst_weak_bar_pct"] = {
        "value": round(100 * worst["tot"] / worst["ours"]),
        "hamiltonian": "standard", "fit_basis": "eq27",
        "inputs": [budget_rel],
        "note": f"largest relative bar at m/g<=0.1: {worst['channel']} "
                f"N={worst['N']} m/g={worst['mg']}"}
    numbers["tableone_n_outside_paper_term"] = {
        "value": sum(1 for r in rows if r["flag"] != "ok"),
        "hamiltonian": "standard", "fit_basis": "eq27",
        "inputs": [budget_rel],
        "note": "entries with |paper-ours| > paper's own last-term magnitude"}

    # Improved Table I, plain-1/K ensemble (G3: improved <-> 1/K).
    imp = improved_fits()
    imp_rel = str(IMPROVED_CSV.relative_to(ROOT))
    for r in rows:
        key = (r["channel"], r["N"], r["mg"])
        if key in imp:
            r["imp"], r["imp_err"] = imp[key]
    n_disagree = sum(
        1 for r in rows if "imp" in r
        and abs(r["imp"] - r["ours"]) > 3 * (r["tot"] ** 2 + r["imp_err"] ** 2) ** 0.5)
    numbers["tableone_improved_n_beyond_three_sigma"] = {
        "value": n_disagree, "hamiltonian": "improved", "fit_basis": "1/K",
        "inputs": [budget_rel, imp_rel],
        "note": "entries where improved and standard disagree past 3 sigma "
                "of the combined TOTAL bar (endpoint included): the coverage "
                "check on the endpoint systematic"}
    # The same divergence in fit-only sigmas is the headline the endpoint
    # component exists to fix; freeze it from the same rows.
    fitonly = [
        abs(r["imp"] - r["ours"])
        / (r["form"] ** 2 + r["wind"] ** 2 + r["trunc"] ** 2
           + r["imp_err"] ** 2) ** 0.5
        for r in rows if "imp" in r]
    numbers["tableone_improved_worst_fitonly_sigma"] = {
        "value": round(max(fitonly)), "hamiltonian": "improved",
        "fit_basis": "1/K", "inputs": [budget_rel, imp_rel],
        "note": "largest improved-vs-standard divergence in fit-only sigmas "
                "(form+wind+trunc quadrature, no endpoint)"}

    # TODO(next passes): ratio table (converged 0.4% statement), GMOR
    # intercept, alpha exponents, captured fractions, timing table.
    numbers["_table1_rows"] = {
        "value": rows, "hamiltonian": "standard", "fit_basis": "eq27",
        "inputs": [budget_rel, imp_rel],
        "note": "full rows incl. improved columns; emitted as "
                "frozen/table1_body.tex, not macros"}
    return numbers


def emit_table1_body(rows: list[dict]) -> None:
    """Write the Table I body in HBP's own layout.

    Columns are his five channels (meson N=2,3,4; baryon N=3,4); rows are
    m/g blocks, largest first, each with three lines: the published value
    with its parenthesis, our standard value with the full budget, and our
    improved value with its ensemble spread.  The exact m/g = 0 row of the
    original is restated last.
    """
    import math

    cols = [("mes", 2), ("mes", 3), ("mes", 4), ("bar", 3), ("bar", 4)]
    by = {(r["channel"], r["N"], r["mg"]): r for r in rows}
    mgs = sorted({r["mg"] for r in rows}, reverse=True)

    def fmt(v: float, e: float) -> str:
        """PDG-style compact notation: 1.23446(8), error to 1-2 sig figs."""
        if e <= 0:
            return f"${v:.4f}$"
        exp = int(math.floor(math.log10(e)))
        l3 = int(round(e / 10 ** exp * 100))        # leading 3 digits, 100-999
        if l3 >= 950:
            sig, exp = 1, exp + 1
        elif l3 >= 355:
            sig = 1
        else:
            sig = 2
        dec = max(0, -(exp - (sig - 1)))
        par = int(round(e * 10 ** dec)) if dec else int(round(e))
        return f"${v:.{dec}f}({par})$" if dec else f"${round(v)}({par})$"

    lines = ["% generated by make_numbers.py — do not edit"]
    for mg in mgs:
        rs = [by.get((c, N, mg)) for c, N in cols]
        hbp = " & ".join(
            fmt(r["paper"], r["pterm"]) if r else "---" for r in rs)
        std = " & ".join(
            fmt(r["ours"], r["tot"]) if r else "---" for r in rs)
        imp = " & ".join(
            fmt(r["imp"], r["imp_err"]) if r and "imp" in r else "---"
            for r in rs)
        lines.append(f"\\multirow{{3}}{{*}}{{{mg:g}}} & HBP & {hbp} \\\\")
        lines.append(f" & standard & {std} \\\\")
        lines.append(f" & improved & {imp} \\\\")
        lines.append("\\hline")
    # No trailing \\ on the final row: it would render as an empty row of
    # whitespace between "exact" and the bottom rule.
    lines.append("0 & exact & $0$ & $0$ & $0$ & $0$ & $0$")
    (FROZEN / "table1_body.tex").write_text("\n".join(lines) + "\n")


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
        if name.startswith("_"):
            continue                       # emitted as a body file instead
        macro = "".join(w.capitalize() for w in name.split("_"))
        lines.append(f"\\newcommand{{\\num{macro}}}{{{entry['value']}}}")
    (FROZEN / "numbers.tex").write_text("\n".join(lines) + "\n")
    if "_table1_rows" in numbers:
        emit_table1_body(numbers["_table1_rows"]["value"])
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
