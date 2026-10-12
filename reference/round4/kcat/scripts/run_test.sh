#!/bin/bash
# ABOUTME: Runs analyze.jl on the small test sets under the shared Julia lock, one at a time.
K=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat
L=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/julia.lock
for s in "$@"; do
  (time timeout -s KILL 1800 flock -w 14400 $L julia --project=$K/env --heap-size-hint=1500M $K/scripts/analyze.jl $s) > $K/logs/analyze_$s.log 2>&1
  echo "exit $?" >> $K/logs/analyze_$s.log
done
