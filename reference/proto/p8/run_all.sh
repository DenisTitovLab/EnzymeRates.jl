#!/bin/bash
# ABOUTME: Runs the P8 key driver on each argument (a population, or several joined by '+') in turn
# ABOUTME: under the shared Julia lock; logs go to logs/run_<arg>.log.
cd "$(dirname "$0")"
LOCK=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/julia.lock
for arg in "$@"; do
  flock -w 14400 "$LOCK" julia --startup-file=no --heap-size-hint=1500M run_pop.jl ${arg//+/ } > "logs/run_${arg}.log" 2>&1
  echo "done $arg $(date +%T)" >> logs/run_all.progress
done
