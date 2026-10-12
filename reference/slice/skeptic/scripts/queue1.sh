#!/bin/bash
# ABOUTME: Queues the skeptic's runs after the spot checks: the second rerun batch, the
# ABOUTME: adversarial set, the key attacks and the large-law timings, one at a time.
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
K=$S/slice/skeptic
cd $K
until [ -f logs/spot_all.done ]; do sleep 10; done
g() { log=$1; shift; LIMIT_MB=${LIMIT_MB:-3500} MAX_S=${MAX_S:-3600} $K/rerun/run_guarded.sh $log "$@"; tail -2 $log | head -1; }
g logs/kcat_specs40.log scripts/inv/kcat40.jl specs/40 specs
g logs/kcat_AD3b.log scripts/inv/kcat40.jl AD3b_copy_square_tied p5a_adv_cases
g logs/spot_p5a_rerun.log scripts/spot.jl p5a_rerun 20 p5a_adv_cases p5a_adv_cases2 p5a_adv_cases3 p5a_adv_cases4 p5a_cases
g logs/adv.log scripts/adv.jl
g logs/keytest.log scripts/keytest.jl
g logs/bigbench.log scripts/bigbench.jl
rerun/rerun_B.sh > logs/rerun_B.out 2>&1
echo finished > logs/queue1.done
