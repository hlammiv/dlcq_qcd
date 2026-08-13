#!/usr/bin/env bash
# Run one DLCQ case and collect its outputs under a named directory.
#
# qcdf.f hard-codes its output filenames (unit 15 -> ./qcdf.out, unit 14 ->
# ./qcdf.ham) and writes them to the current working directory.  Runs must
# therefore be isolated, which is what this wrapper does: it executes the
# solver inside the output directory itself.
#
# Usage:
#   run_case.sh OUTDIR N NF B LAMBDA CUTOFF LPN K [RMQ2..RMQNF IFLV1..IFLVNF]
#   run_case.sh OUTDIR --reuse HAMFILE LAMBDA        # IHFILE=1 fast lambda sweep
#
# Examples:
#   run_case.sh runs/K21 3 1 1 0.3325 -1.0 0 21
#   run_case.sh runs/K21_l07 --reuse runs/K21/qcdf.ham 0.7
#
#   # NF=2: one mass ratio (rmq2/rmq1) then one flavour number per flavour.
#   # A flavour-neutral meson with a doubly heavy second flavour:
#   run_case.sh runs/nf2 3 2 0 0.3325 -1.0 0 12   2.0   0 0
#   # ...and the non-singlet ("kaon") channel, one quark of each flavour:
#   run_case.sh runs/nf2ns 3 2 0 0.3325 -1.0 0 12   2.0   1 -1
#
# qcdf.f fixes rmq(1) = 1 and reads rmq(2..NF) as RATIOS to it; iflv(f) is the
# net quark number of flavour f, so a meson is flavour-neutral overall and a
# baryon's numbers must sum to N*B.  For NF=1 the code sets iflv(1) = N*B
# itself and prompts for neither.
#
# Produces OUTDIR/{qcdf.out,qcdf.ham,stdout.log,case.json}.
#
# The IHFILE=1 path exists because the Hamiltonian is lambda-independent: the
# coupling only enters as a lambda^2 reweighting at the very end (qcdf.f:600).
# One Hamiltonian per (N,B,K) therefore serves an entire m/g sweep, which is
# what makes Figs. 2, 7, 8 and Table I affordable.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# QCDF_BIN lets callers select the widened variant (see widen.py).
QCDF="${QCDF_BIN:-$SCRIPT_DIR/qcdf}"

if [ ! -x "$QCDF" ]; then
    echo "error: $QCDF not built. Run 'make -C $SCRIPT_DIR' first." >&2
    exit 1
fi

if [ $# -lt 2 ]; then
    sed -n '2,25p' "${BASH_SOURCE[0]}" >&2
    exit 1
fi

OUTDIR="$1"; shift
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

if [ "${1:-}" = "--reuse" ]; then
    # IHFILE=1: read a stored Hamiltonian, only lambda is prompted for.
    HAMFILE="$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"; shift 2
    LAMBDA="$1"
    cp "$HAMFILE" "$OUTDIR/qcdf.ham"
    printf '1\n%s\n' "$LAMBDA" > "$OUTDIR/input.txt"
    cat > "$OUTDIR/case.json" <<EOF
{"mode": "reuse", "hamfile": "$HAMFILE", "rlamb": $LAMBDA}
EOF
else
    if [ $# -lt 7 ]; then
        echo "error: expected at least 7 parameters (N NF B LAMBDA CUTOFF LPN K), got $#" >&2
        exit 1
    fi
    N="$1"; NF="$2"; B="$3"; LAMBDA="$4"; CUTOFF="$5"; LPN="$6"; K="$7"
    shift 7
    # For NF>1 qcdf.f prompts for NF-1 mass ratios then NF flavour numbers,
    # between RLAMB and CUTOFF.  For NF=1 it prompts for neither.
    NEXTRA=$(( NF > 1 ? 2 * NF - 1 : 0 ))
    if [ $# -ne "$NEXTRA" ]; then
        echo "error: NF=$NF needs $NEXTRA extra parameters" \
             "($(( NF - 1 )) mass ratios + $NF flavour numbers), got $#" >&2
        exit 1
    fi
    FLAVOUR=""
    for v in "$@"; do FLAVOUR="${FLAVOUR}${v}\n"; done
    # Read order in qcdf.f: IHFILE, N, NF, B, RLAMB, [rmq(2..NF), iflv(1..NF)],
    # CUTOFF, LPN, K.
    printf "0\n%s\n%s\n%s\n%s\n${FLAVOUR}%s\n%s\n%s\n" \
        "$N" "$NF" "$B" "$LAMBDA" "$CUTOFF" "$LPN" "$K" > "$OUTDIR/input.txt"
    EXTRA_JSON=""
    if [ "$NEXTRA" -gt 0 ]; then
        EXTRA_JSON=", \"flavour_params\": [$(printf '%s, ' "$@" | sed 's/, $//')]"
    fi
    cat > "$OUTDIR/case.json" <<EOF
{"mode": "compute", "N": $N, "NF": $NF, "B": $B, "rlamb": $LAMBDA,
 "cutoff": $CUTOFF, "LPN": $LPN, "K": $K$EXTRA_JSON}
EOF
fi

cd "$OUTDIR"
# stdout is extremely chatty (QCDSTA prints from an inner momentum loop), so
# it goes to a log rather than the terminal.
"$QCDF" < input.txt > stdout.log 2>&1

echo "$OUTDIR/qcdf.out"
