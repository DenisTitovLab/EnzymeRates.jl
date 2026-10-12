#!/bin/bash
# ABOUTME: Queues a retry of the first large-law timing script, to see whether its GC crash
# ABOUTME: reproduces.
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
K=$S/slice/skeptic
cd $K
until [ -f logs/queue8.done ]; do sleep 10; done
g() { log=$1; shift; LIMIT_MB=${LIMIT_MB:-3500} MAX_S=${MAX_S:-3600} $K/rerun/run_guarded.sh $log "$@"; tail -2 $log | head -1; }
NTOP=10 g logs/bigbench_orig_retry.log scripts/bigbench_orig.jl
echo finished > logs/queue9.done
