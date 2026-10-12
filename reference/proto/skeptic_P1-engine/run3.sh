#!/bin/bash
# ABOUTME: Third batch: adversarial checks with the corrected rank helper, then the population rank census.
# ABOUTME: Each run takes the shared Julia lock; logs go to skeptic_P1-engine/logs/.
P=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
S=$P/proto/skeptic_P1-engine
run() {
    local t0=$(date +%s)
    flock -w 14400 $P/proto/julia.lock julia --startup-file=no --heap-size-hint=1500M $1 > $S/logs/$2.log 2>&1
    echo "$(date +%T) $2 exit $? ($(( $(date +%s) - t0 )) s incl. lock wait)" >> $S/logs/run3.status
}
run $S/adv_check.jl adv_check
run $S/rank_pop.jl rank_pop
echo "$(date +%T) all done" >> $S/logs/run3.status
