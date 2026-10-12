#!/bin/bash
# ABOUTME: Runs the P8 reparametrization estimate on each population and variant in turn
# ABOUTME: under the shared Julia lock; logs go to logs/reparam_<variant>_<pops>.log.
cd "$(dirname "$0")"
LOCK=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/julia.lock
for arg in "$@"; do
  for v in redesign today; do
    NPAIRS=300 flock -w 14400 "$LOCK" julia --startup-file=no --heap-size-hint=1500M reparam.jl $v ${arg//+/ } > "logs/reparam_${v}_${arg}.log" 2>&1
    echo "done $v $arg $(date +%T)" >> logs/reparam.progress
  done
done
