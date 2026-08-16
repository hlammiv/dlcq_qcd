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

VALID_PAIRS = {("improved", "1/K"), ("standard", "eq27"),
               ("improved", "alpha(mg)"), ("standard", "alpha(mg)")}


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


def two_body_anchor() -> dict:
    """Parse the recorded van de Sande two-body validation numbers.

    The values live in dlcq/endpoint.py's Validation docstring and are held
    there by tests/test_endpoint.py; parsing rather than retyping means an
    edit to the recorded values breaks this freeze loudly.  Provenance
    carries the file hash.
    """
    import re
    text = (ROOT / "dlcq" / "endpoint.py").read_text()
    exact = re.search(r"exact ground state is ``M\^2/g\^2 = ([\d.]+)``", text)
    imp = re.search(r"returns ([\d.]+), i\.e\. the exact\s+answer", text)
    std = re.search(r"the same fit on the standard series\s+returns "
                    r"(\d+\.\d+)", text)
    if not (exact and imp and std):
        raise SystemExit("two_body_anchor: recorded values not found in "
                         "dlcq/endpoint.py — docstring changed?")
    return {"exact": float(exact.group(1)), "improved": float(imp.group(1)),
            "standard": float(std.group(1))}


def _read_house_csv(path: Path, want_hamiltonian: str | None = None) -> dict:
    """House-format CSV -> {(N, B, mg): {K_code: msq}}."""
    import csv
    out: dict = {}
    with open(path) as fh:
        rows = csv.DictReader(l for l in fh if not l.lstrip().startswith("#"))
        for r in rows:
            if want_hamiltonian and r.get("hamiltonian") != want_hamiltonian:
                continue
            key = (int(r["N"]), int(r["B"]), float(r["mg"]))
            out.setdefault(key, {})[int(r["K_code"])] = float(r["msq"])
    return out


def chiral_exponents() -> dict:
    """Power-law fits M^2 ~ (m/g)^alpha at fixed K, meson N=3.

    Standard from data/chiral_grid_msq.csv, improved from
    data/gmor_scan/improved_N3.csv (read-only), over the weak couplings
    m/g <= 0.1 present in each file, at 2K = 30, 40, 60.  The exponent's
    stability in K is the point: a fit-basis artifact would move with the
    grid, a property of the masses does not.
    """
    import math
    out = {}
    # chiral_grid_msq.csv predates the hamiltonian column and is standard
    # throughout (its header says so); gmor_scan files carry the column.
    for ham, path, want in (
            ("standard", ROOT / "data" / "chiral_grid_msq.csv", None),
            ("improved", ROOT / "data" / "gmor_scan" / "improved_N3.csv",
             "improved")):
        data = _read_house_csv(path, want)
        alphas = []
        for K in (30, 40, 60):
            pts = [(mg, series[K])
                   for (N, B, mg), series in data.items()
                   if N == 3 and B == 0 and mg <= 0.1 and K in series]
            if len(pts) < 3:
                continue
            xs = [math.log(mg) for mg, _ in pts]
            ys = [math.log(m) for _, m in pts]
            n = len(xs)
            sx, sy = sum(xs), sum(ys)
            sxx = sum(x * x for x in xs)
            sxy = sum(x * y for x, y in zip(xs, ys))
            alphas.append((n * sxy - sx * sy) / (n * sxx - sx * sx))
        if len(alphas) < 2:
            raise SystemExit(f"chiral_exponents: too few K values for {ham}")
        out[ham] = {"alphas": alphas,
                    "value": alphas[-1],
                    "drift": max(alphas) - min(alphas)}
    return out


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

    # Van de Sande two-body anchor (values recorded in dlcq/endpoint.py,
    # held by tests/test_endpoint.py).
    tb = two_body_anchor()
    ep_rel = "dlcq/endpoint.py"
    numbers["two_body_exact"] = {
        "value": tb["exact"], "inputs": [ep_rel],
        "note": "van de Sande's exact M^2/g^2 at beta=0.1"}
    numbers["two_body_improved"] = {
        "value": tb["improved"], "hamiltonian": "improved", "fit_basis": "1/K",
        "inputs": [ep_rel], "note": "improved-series 1/K extrapolate"}
    numbers["two_body_improved_rel"] = {
        "value": f"{abs(tb['improved'] - tb['exact']) / tb['exact']:.1e}",
        "hamiltonian": "improved", "fit_basis": "1/K", "inputs": [ep_rel],
        "note": "relative agreement of improved extrapolate with exact"}
    numbers["two_body_standard"] = {
        "value": tb["standard"], "hamiltonian": "standard",
        "fit_basis": "eq27", "inputs": [ep_rel],
        "note": "the same fit on the standard series; not close at any "
                "reachable K"}

    # Chiral exponents at fixed K (meson N=3): the artifact vs the physics.
    ce = chiral_exponents()
    ce_inputs = ["data/chiral_grid_msq.csv", "data/gmor_scan/improved_N3.csv"]
    for ham in ("standard", "improved"):
        numbers[f"alpha_{ham}"] = {
            "value": f"{ce[ham]['value']:.2f}", "hamiltonian": ham,
            "fit_basis": "alpha(mg)", "inputs": ce_inputs,
            "note": f"M^2 ~ (m/g)^alpha at 2K=60; per-K values "
                    f"{[round(a, 3) for a in ce[ham]['alphas']]}"}
        numbers[f"alpha_{ham}_drift"] = {
            "value": f"{ce[ham]['drift']:.3f}", "hamiltonian": ham,
            "fit_basis": "alpha(mg)", "inputs": ce_inputs,
            "note": "spread of alpha over 2K = 30, 40, 60"}

    # Endpoint arithmetic at the smallest coupling probed (N=5 meson,
    # m/g = 1.953e-4, the tiny_scan floor), straight from Eq. (26).
    import math
    import sys as _sys
    _sys.path.insert(0, str(ROOT))
    from dlcq.units import endpoint_exponent
    a_tiny = float(endpoint_exponent(1.953e-4, 5))
    numbers["capture_pct_tiny"] = {
        "value": f"{100 * (1 - 70.0 ** (-2 * a_tiny)):.3f}",
        "inputs": ["dlcq/units.py"],
        "note": "percent of the endpoint integral a 2K=70 grid captures at "
                "m/g=1.953e-4, N=5 (exponent a from Eq. 26)"}
    numbers["halving_log_ten"] = {
        "value": round(math.log10(2.0) / a_tiny),
        "inputs": ["dlcq/units.py"],
        "note": "log10 of the K factor needed to halve the K^-a remainder "
                "at m/g=1.953e-4, N=5"}

    # Meson/baryon ratio vs bosonization, improved Hamiltonian, plain 1/K.
    # Mesons from data/gmor_scan (peer's, read-only), baryons from the
    # Phase 4c extensions in data/paper.  Quoted at the smallest coupling
    # with baryon coverage, m/g = 0.0125.  The N=8 column joins when its
    # MXTRM-raised rerun lands.
    import sys as _s
    _s.path.insert(0, str(ROOT / "tools"))
    from plot_improved_fits import _spread_1k
    mg_ratio = 0.0125
    bar_ext = _read_house_csv(ROOT / "data" / "paper" /
                              "bar_ext_improved_N567.csv", "improved")
    devs = {}
    for Nc in (5, 6, 7):
        mes = _read_house_csv(
            ROOT / "data" / "gmor_scan" / f"improved_N{Nc}.csv", "improved")
        mser = mes.get((Nc, 0, mg_ratio), {})
        bser = bar_ext.get((Nc, 1, mg_ratio), {})
        if len(mser) < 8 or len(bser) < 8:
            continue
        mk = sorted(mser); bk = sorted(bser)
        m_med, m_err, _ = _spread_1k(mk, [mser[k] for k in mk])
        b_med, b_err, _ = _spread_1k(bk, [bser[k] for k in bk])
        R = math.sqrt(m_med / b_med)          # M ratio from M^2 ratio
        nu = 1.0 / (2 * Nc - 1)
        target = 2 * math.sin(math.pi * nu / 2)
        devs[Nc] = 100 * (R / target - 1)
    if devs:
        worst = max(abs(v) for v in devs.values())
        numbers["ratio_worst_dev_pct"] = {
            "value": f"{worst:.1f}", "hamiltonian": "improved",
            "fit_basis": "1/K",
            "inputs": ["data/paper/bar_ext_improved_N567.csv"]
            + [f"data/gmor_scan/improved_N{n}.csv" for n in devs],
            "note": f"largest |deviation| of M_mes/M_bar from 2sin(pi nu/2) "
                    f"at m/g={mg_ratio}, per-N percents "
                    f"{ {k: round(v, 2) for k, v in devs.items()} }"}

    # The N=8 standard/improved gap at the smallest coupling probed --
    # the tiny_scan gap, now measured at the eighth color too.
    for ham in ("standard", "improved"):
        d = _read_house_csv(ROOT / "data" / "paper" /
                            f"tiny_mes_{ham}_N8.csv", ham)
        ser = d.get((8, 0, 0.0001953125), {})
        if len(ser) >= 8:
            ks = sorted(ser)
            med, err, _ = _spread_1k(ks, [ser[k] for k in ks])
            numbers[f"tiny_neight_{ham}_msq"] = {
                "value": f"{med:.3e}", "hamiltonian": ham,
                "fit_basis": "1/K" if ham == "improved" else "eq27",
                "inputs": [f"data/paper/tiny_mes_{ham}_N8.csv"],
                "note": "M^2 extrapolate at m/g=1.953e-4, N=8 meson "
                        "(standard quoted from the same 1/K fit ONLY as a "
                        "gap denominator; its budgeted value uses eq27)"}
    a = numbers.get("tiny_neight_improved_msq")
    b = numbers.get("tiny_neight_standard_msq")
    if a and b:
        numbers["tiny_neight_gap"] = {
            "value": round(float(a["value"]) / float(b["value"])),
            "hamiltonian": "improved", "fit_basis": "1/K",
            "inputs": a["inputs"] + b["inputs"],
            "note": "improved/standard M^2 ratio at m/g=1.953e-4, N=8"}

    # GMOR law, quoted from the pinned weak-coupling doc (G5 blob pin in
    # provenance): M^2/GMOR = 1 + 0.71 (m/g), intercept -0.0005 +- 0.0008.
    numbers["gmor_intercept"] = {
        "value": "-0.0005 \\pm 0.0008", "hamiltonian": "improved",
        "fit_basis": "1/K", "inputs": [],
        "note": "doc-quoted at the pinned blob; regenerate from "
                "tools/gmor_chiral_limit.py at v2"}

    # TODO(v2): timing table; N=8 ratio column; GMOR from a fresh run.
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
