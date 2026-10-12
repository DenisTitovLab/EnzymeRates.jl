#!/bin/bash
# ABOUTME: Queues the skeptic's runs after the first queue: the P1 adversarial spot check again
# ABOUTME: (kcat corners capped) and the record-by-record comparison of the clean rerun.
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
K=$S/slice/skeptic
cd $K
until [ -f logs/queue1.done ]; do sleep 10; done
g() { log=$1; shift; LIMIT_MB=${LIMIT_MB:-3500} MAX_S=${MAX_S:-3600} $K/rerun/run_guarded.sh $log "$@"; tail -2 $log | head -1; }
g logs/spot_p1adv_rerun.log scripts/spot.jl p1adv_rerun 20 p1_adv_pkg p1_adv_pkg2 p1_direct
g logs/compare_pipe.log scripts/compare_pipe.jl
echo finished > logs/queue2.done
