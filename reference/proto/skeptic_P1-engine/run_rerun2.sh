#!/bin/bash
# ABOUTME: Second batch: probe scripts that read results.jls (patched to load Engine first) and my own checks.
# ABOUTME: Each run takes the shared Julia lock; logs go to skeptic_P1-engine/logs/.
P=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
S=$P/proto/skeptic_P1-engine
run() {  # script log [project]
    local t0=$(date +%s)
    if [ -n "$3" ]; then
        flock -w 14400 $P/proto/julia.lock julia --startup-file=no --project=$P/er_oracle --heap-size-hint=2500M $1 > $S/logs/$2.log 2>&1
    else
        flock -w 14400 $P/proto/julia.lock julia --startup-file=no --heap-size-hint=1500M $1 > $S/logs/$2.log 2>&1
    fi
    echo "$(date +%T) $2 exit $? ($(( $(date +%s) - t0 )) s incl. lock wait)" >> $S/logs/run_rerun2.status
}
run $S/rerun/phantom_groups.jl rerun_phantom_groups
run $S/rerun/blocks_peek.jl rerun_blocks_peek
run $S/rerun/coverage.jl rerun_coverage
run $S/adv_direct.jl adv_direct
run $S/adv_export.jl adv_export pkg
run $S/adv_check.jl adv_check
echo "$(date +%T) all done" >> $S/logs/run_rerun2.status
