#!/usr/bin/env python3
"""Scriptable verification gates (the plan's G1, G2, G3, G8, G9).

Exit nonzero on any failure; `make check` runs this.

G1  No bare numeric literal in section prose: every number in running text
    comes from a frozen macro.  Whitelisted: years, \\cite/\\ref/\\label
    arguments, equation/section numbers in refs like "Eq.~(27)", counts
    written as words, values inside comments or \\num macros, figure
    multipliers of the form \\times10, resolutions "2K = nn" (parameters,
    not results — captions state them by the caption gate), m/g values,
    LPN/N/B parameters, and percent literals spelled in prose words.
G2  Every \\caption states N, B (or channel names), 2K, m/g, LPN (or
    "untruncated"), and the Hamiltonian.
G3  make_numbers.py runs clean (its internal basis-pairing gate).
G8  Withdrawn-claim signatures never appear as assertions.
G9  The equal-time slot (10e file) is either absent or referenced only by
    an \\input in main.tex — deleting it must leave no dangling refs.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SECTIONS = sorted((HERE / "sections").glob("*.tex"))
fails: list[str] = []


def strip_protected(t: str) -> str:
    t = re.sub(r"(?<!\\)%.*", "", t)                      # comments
    t = re.sub(r"\\(cite|ref|eqref|label|input|includegraphics)"
               r"(\[[^]]*\])?\{[^}]*\}", " ", t)
    t = re.sub(r"\\num[A-Za-z]+", " ", t)                 # frozen macros
    t = re.sub(r"\\caption\{", " ", t)                    # caption braces
    return t


# ── G1: bare digits in prose ───────────────────────────────────────────────
ALLOW = [
    r"19\d\d|20\d\d",                       # years
    r"Eq\.~?\(?\d+\)?", r"Sec\.~?[IVX\d]+", r"Fig[s]?\.~?\d+\w?",
    r"Table~?[I1]+",
    r"2K\s*[=~]?\s*\$?\d+", r"\$2K[^$]*\$",
    r"m/g\s*=\s*[\d.e\\times^{}\-]+", r"\$m/g[^$]*\$",
    r"LPN[^$ ]*\d+", r"N\s*=\s*\d", r"\$N[^$]*\$", r"B\s*=\s*\d",
    r"SU\(\d\)", r"U\(\d\)", r"\d\+1", r"1\+1",
    r"\\times10", r"10\^", r"\^\{?-?\d",
    r"one|two|three|four|five|six|seven|eight|nine|ten",
    r"logical qubits?",                     # analytic 2NcK ledger estimates
    r"IBM \d+",                            # proper noun
    r"lowest \d+(--\d+)? states?",        # nev parameters
    r"level \d+",                          # level-index identifications
    r"runs? from \d+ to \d+",             # K-window parameters
    r"factor of [\d.]+",                   # parameter arithmetic
]
allow_re = re.compile("|".join(ALLOW))
num_re = re.compile(r"\d+\.\d+|\d{2,}")

for f in SECTIONS:
    body = strip_protected(f.read_text())
    # drop inline math entirely: equations and parameter math are not prose
    body = re.sub(r"\$[^$]*\$", " ", body)
    body = re.sub(r"\\begin\{equation\}.*?\\end\{equation\}", " ",
                  body, flags=re.S)
    flat = re.sub(r"\s+", " ", body)     # phrases must survive line wraps
    for m in num_re.finditer(flat):
        ctx = flat[max(0, m.start() - 60):m.end() + 25]
        if allow_re.search(ctx):
            continue
        fails.append(f"G1 {f.name}: bare number "
                     f"'{m.group()}' in: ...{ctx.strip()}...")

# ── G2: caption fields ─────────────────────────────────────────────────────
REQUIRED = [
    (re.compile(r"N\s*=|N\{?=\}?|SU\("), "N"),
    (re.compile(r"2K"), "2K"),
    (re.compile(r"m/g|\\lambda"), "m/g"),
    (re.compile(r"LPN|untruncated|no fits|Fock truncation"), "LPN"),
    (re.compile(r"[Hh]amiltonian"), "Hamiltonian"),
]
for f in SECTIONS:
    t = f.read_text()
    for m in re.finditer(r"\\caption\{(.*?)\}\s*\n\\label", t, flags=re.S):
        cap = m.group(1)
        for rx, name in REQUIRED:
            if not rx.search(cap):
                fails.append(f"G2 {f.name}: caption missing '{name}' "
                             f"field: {cap[:60]}...")

# ── G3: make_numbers runs clean (basis pairing inside it) ─────────────────
r = subprocess.run([sys.executable, str(HERE / "make_numbers.py"),
                    "--freeze", "gatecheck"], capture_output=True, text=True)
if r.returncode != 0:
    fails.append(f"G3 make_numbers failed: {r.stderr.strip()[:200]}")
for tmp in (HERE / "frozen").glob("gatecheck*"):
    tmp.unlink()

# ── G8: withdrawn-claim signatures ────────────────────────────────────────
FORBIDDEN = ["30\\sigma", "30 sigma", "0.02--0.2\\%", "5--7\\%",
             "confirms bosonization to"]
for f in SECTIONS:
    t = strip_protected(f.read_text())
    for sig in FORBIDDEN:
        if sig in t:
            fails.append(f"G8 {f.name}: withdrawn-claim signature '{sig}'")

# ── G9: additivity of the equal-time slot ─────────────────────────────────
slot = HERE / "sections" / "10e-equal-time-crosscheck.tex"
main = (HERE / "main.tex").read_text()
if slot.exists():
    refs = re.findall(r"10e-equal-time-crosscheck", main)
    if len(refs) != 1:
        fails.append("G9: 10e slot exists but is referenced other than by "
                     "its single \\input")
else:
    if "10e-equal-time-crosscheck" in main:
        fails.append("G9: main.tex references the absent 10e slot")

print(f"{len(SECTIONS)} section files checked")
if fails:
    for x in fails:
        print("FAIL", x)
    sys.exit(1)
print("all gates pass")
