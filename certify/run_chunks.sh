#!/usr/bin/env bash
# Run a chunked certification step over primes in [2, Zimmert bound] in parallel.
# Usage: certify/run_chunks.sh <step>   where <step> is phase1 or generators
#   phase1     -> certify/phase1.gp      (PARI bnftestprimes, range-restricted)
#   generators -> certify/generators.gp  (explicit verified generators)
# Resumable: chunks whose output already ends in "STATUS OK" are skipped.
set -euo pipefail
cd "$(dirname "$0")/.."
GP=.tools/pari/bin/gp
STEP=${1:?usage: run_chunks.sh phase1|generators}
case $STEP in phase1|generators) ;; *) echo "unknown step $STEP" >&2; exit 2;; esac
JOBS=${JOBS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}
CHUNKS=${CHUNKS:-64}
grep -q "STATUS OK" results/setup.out || { echo "run certify/setup.gp first" >&2; exit 1; }
BOUND=$(awk '/^zimmert_bound/{print $2}' results/setup.out)
mkdir -p results/$STEP

run_chunk() {
  local a=$1 b=$2 out="results/$STEP/$1-$2.out"
  if [[ -f $out ]] && tail -1 "$out" | grep -q "STATUS OK"; then return 0; fi
  local err="results/$STEP/$1-$2.err"
  P1_A=$a P1_B=$b P1_BOUND=$BOUND P1_OUT=$out \
    "$GP" -q -D parisizemax=1000000000 "certify/$STEP.gp" < /dev/null > "$err" 2>&1
  if grep -v "Warning" "$err" | grep -q '\*\*\*'; then
    echo "ERROR in chunk $a-$b, see $err" >&2; return 1
  fi
  tail -1 "$out" | grep -q "STATUS OK" || { echo "chunk $a-$b incomplete" >&2; return 1; }
  rm -f "$err"
  echo "done $a-$b"
}
export -f run_chunk
export GP BOUND STEP

python3 - "$BOUND" "$CHUNKS" <<'PY' | xargs -P "$JOBS" -n 2 bash -c 'run_chunk "$0" "$1"'
import sys
bound, n = int(sys.argv[1]), int(sys.argv[2])
edges = [2] + [bound * i // n for i in range(1, n)] + [bound + 1]
for lo, hi in zip(edges, edges[1:]):
    print(lo, hi - 1)
PY

python3 certify/check_results.py "$STEP"
