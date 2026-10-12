#!/bin/bash
# ABOUTME: Queues the adversarial set (cube roots rebuilt) and the spot checks whose kcat oracle
# ABOUTME: changed (products at 1e-80, substrates at 1e24).
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
K=$S/slice/skeptic
cd $K
g() { log=$1; shift; LIMIT_MB=${LIMIT_MB:-3500} MAX_S=${MAX_S:-3600} $K/rerun/run_guarded.sh $log "$@"; tail -2 $log | head -1; }
g logs/adv2.log scripts/adv.jl
g logs/spot_specs2.log scripts/spot.jl specs2 42 specs
g logs/spot_p2adv3.log scripts/spot.jl p2adv3 43 p2_adv_pkg p2_adv_pkg2 p2_fuzz_dump p2_fuzz_dump_flip p2_hand
g logs/spot_p5a3.log scripts/spot.jl p5a3 40 p5a_adv_cases p5a_adv_cases2 p5a_adv_cases3 p5a_adv_cases4 p5a_cases
echo finished > logs/queue12.done
