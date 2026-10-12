#!/bin/bash
# ABOUTME: Queues the large-law timings without the law derivation that precedes the GC crash.
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
K=$S/slice/skeptic
cd $K
until [ -f logs/queue12.done ]; do sleep 10; done
g() { log=$1; shift; LIMIT_MB=${LIMIT_MB:-3500} MAX_S=${MAX_S:-3600} $K/rerun/run_guarded.sh $log "$@"; tail -2 $log | head -1; }
NTOP=20 g logs/bigbench_bbG.log scripts/inv/bbG.jl
echo finished > logs/queue13.done
