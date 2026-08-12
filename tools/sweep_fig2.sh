#!/usr/bin/env bash
# Fig. 2 needs a dense lambda scan at three FIXED K, which the Table I sweep
# (2K = 16-24) does not cover: B=0 at 2K=10, B=1 at 2K=13, B=2 at 2K=22.
# These bases are small, so no Fock truncation is needed.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
JOBS=${JOBS:-18}
OUT=runs/fig2; mkdir -p "$OUT"
running=0
for spec in "0 10" "1 13" "2 22"; do
  set -- $spec; B=$1; K=$2
  for lam in 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 \
             0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 0.97 0.99; do
    tag="B${B}_K${K}_lam${lam}"
    [ -f "$OUT/$tag/qcdf.out" ] && continue
    ( bash fortran/run_case.sh "$OUT/$tag" 3 1 "$B" "$lam" -1.0 0 "$K" >/dev/null 2>&1 \
        && echo "  ok $tag" || echo "  FAIL $tag" ) &
    running=$((running+1))
    [ $running -ge $JOBS ] && { wait -n; running=$((running-1)); }
  done
done
wait
echo "Fig 2 sweep complete: $(ls -d $OUT/*/ 2>/dev/null | wc -l) runs"
