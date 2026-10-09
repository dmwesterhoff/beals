#!/usr/bin/env bash
# Run Stoll's criterion for each prime given on the command line (in parallel).
# Usage: JOBS=2 stoll/run_range.sh 7 11 13 ...
set -euo pipefail
cd "$(dirname "$0")/.."
JOBS=${JOBS:-2}
mkdir -p results/stoll

run_one() {
  local l=$1 err="results/stoll/l$1.err"
  STOLL_L=$l .tools/pari/bin/gp -q -D parisizemax=2000000000 stoll/run.gp < /dev/null > "$err" 2>&1
  if grep -v Warning "$err" | grep -q '\*\*\*'; then echo "l=$l ERROR (see $err)"; return 1; fi
  tail -1 "results/stoll/l$l.out" | grep -q "STATUS OK" || { echo "l=$l incomplete"; return 1; }
  rm -f "$err"
  echo "l=$l $(grep -E '^(PASS|bnf_source|ms)' "results/stoll/l$l.out" | tr '\n' ' ')"
}
export -f run_one

printf '%s\n' "$@" | xargs -P "$JOBS" -n 1 bash -c 'run_one "$0"'
