#!/bin/bash
# ABOUTME: Runs analyze.jl on each label in turn under the shared Julia lock (KILL timeout per
# ABOUTME: run); logs in logs/analyze_<label>.log, rows in data/rows_<label>.jls.
K=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat
L=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/julia.lock
for s in "$@"; do
  (time timeout -s KILL 3600 flock -w 14400 $L julia --project=$K/env --heap-size-hint=1500M $K/scripts/analyze.jl $s) > $K/logs/analyze_$s.log 2>&1
  echo "exit $?" >> $K/logs/analyze_$s.log
done
echo finished > $K/logs/run_all_$1.done
