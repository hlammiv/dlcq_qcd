#!/usr/bin/env bash
# Fortran sweep for Figs. 2, 7, 8 and Table I.
#
# Grid
# ----
# K parity is fixed by the sector: a state's momenta are odd integers, so the
# number of partons determines whether K is even or odd.
#   B=0  (q qbar + pairs, even parton count)  -> K even
#   B=1, N odd  (N quarks + pairs)            -> K odd
#   B=1, N even                               -> K even
# The paper extrapolates over "2K in the range of roughly 16-24", so we use the
# five (or four) admissible K in that window.
#
# LPN truncation caps the parton number. It is what keeps every run inside
# qcdf.f's 25-slot colour array (needs 2L+4 <= 25, so L <= 10) and it is almost
# certainly what the authors did -- input_parameters.pdf recommends starting at
# 2 or 4. It is also well converged: at 2K=21, B=1 the lightest mass is 5.46952
# truncated at LPN=5 against 5.469519 untruncated.
#
# Runs are independent, so they go in parallel up to the core count.

set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

JOBS=${JOBS:-16}
OUT=runs/sweep
mkdir -p "$OUT"

# m/g -> lambda = 1/sqrt(1 + pi (m/g)^2).  m/g = 0 is exactly massless by
# Eq. (16) and needs no run.
declare -A LAM=( [1.6]=0.33254949 [0.8]=0.57633109 [0.4]=0.81577472
                 [0.2]=0.94253081 [0.1]=0.98465268 [0.05]=0.99609604 )

launch () {   # N B LPN K mg
    local N=$1 B=$2 LPN=$3 K=$4 mg=$5
    local lam=${LAM[$mg]}
    local tag="N${N}_B${B}_K${K}_mg${mg}"
    [ -f "$OUT/$tag/qcdf.out" ] && return 0
    bash fortran/run_case.sh "$OUT/$tag" "$N" 1 "$B" "$lam" -1.0 "$LPN" "$K" >/dev/null 2>&1 \
        && echo "  ok   $tag" || echo "  FAIL $tag"
}

running=0
for N in 2 3 4; do
  for B in 0 1; do
    # parton count of the valence sector fixes K parity
    if [ "$B" -eq 0 ]; then valence=2; else valence=$(( N * B )); fi
    if [ $(( valence % 2 )) -eq 0 ]; then parity=0; else parity=1; fi
    # one extra qqbar pair beyond valence, capped at 10 partons
    LPN=$(( valence + 2 )); [ "$LPN" -gt 10 ] && LPN=10
    for K in 16 17 18 19 20 21 22 23 24; do
      [ $(( K % 2 )) -ne $parity ] && continue
      for mg in 1.6 0.8 0.4 0.2 0.1 0.05; do
        launch "$N" "$B" "$LPN" "$K" "$mg" &
        running=$(( running + 1 ))
        if [ "$running" -ge "$JOBS" ]; then wait -n; running=$(( running - 1 )); fi
      done
    done
  done
done
wait
echo "Fortran sweep complete: $(ls -d $OUT/*/ 2>/dev/null | wc -l) runs"
