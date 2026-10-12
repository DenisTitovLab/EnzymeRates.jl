#!/bin/bash
# ABOUTME: Reruns every allosteric-stage check under the shared Julia lock, one at a time, each
# ABOUTME: log in logs/ (about 25 minutes without waiting for the lock).
cd /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/allo
LOCK=../../proto/julia.lock
run() { ( time timeout -s KILL 2400 flock -w 14400 $LOCK julia --startup-file=no --heap-size-hint=1500M "$1" ) > logs/${1%.jl}.log 2>&1 ; }
for f in test_hand.jl check_oracle.jl check_verdicts.jl check_g1.jl check_limit.jl check_print.jl check_kcat.jl check_extra.jl check_lists_p5b.jl enum_admit.jl explore_print.jl; do
    run $f
done
echo finished > logs/run_all.done
