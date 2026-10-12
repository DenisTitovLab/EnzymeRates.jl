#!/bin/bash
# ABOUTME: Fourth batch: the free-form package case, and the enumerated dead-end-inhibitor populations.
# ABOUTME: Each run takes the shared Julia lock; logs go to skeptic_P1-engine/logs/.
P=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
S=$P/proto/skeptic_P1-engine
run() {  # log project? script args...
    local log=$1 proj=$2; shift 2
    local t0=$(date +%s)
    if [ "$proj" = pkg ]; then
        flock -w 14400 $P/proto/julia.lock julia --startup-file=no --project=$P/er_oracle --heap-size-hint=2500M "$@" > $S/logs/$log.log 2>&1
    else
        flock -w 14400 $P/proto/julia.lock julia --startup-file=no --heap-size-hint=1500M "$@" > $S/logs/$log.log 2>&1
    fi
    echo "$(date +%T) $log exit $? ($(( $(date +%s) - t0 )) s incl. lock wait)" >> $S/logs/run4.status
}
run adv_export2 pkg $S/adv_export.jl extra
run adv_check2 std $S/adv_check.jl adv_pkg2.jls
run enum_export pkg $S/enum_export.jl
run enum_check std $S/enum_check.jl
echo "$(date +%T) all done" >> $S/logs/run4.status
