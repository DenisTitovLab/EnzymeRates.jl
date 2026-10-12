#!/bin/bash
# ABOUTME: Queues the P2 tie investigation and the P2 spot check with the chemistry-aware
# ABOUTME: alignment.
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
K=$S/slice/skeptic
cd $K
until [ -f logs/queue10.done ]; do sleep 10; done
g() { log=$1; shift; LIMIT_MB=${LIMIT_MB:-3500} MAX_S=${MAX_S:-3600} $K/rerun/run_guarded.sh $log "$@"; tail -2 $log | head -1; }
g logs/p2ties.log scripts/inv/p2ties.jl
g logs/spot_p2adv_rerun2.log scripts/spot.jl p2adv_rerun2 43 p2_adv_pkg p2_adv_pkg2 p2_fuzz_dump p2_fuzz_dump_flip p2_hand
echo finished > logs/queue11.done
