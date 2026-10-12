#!/bin/bash
# ABOUTME: Re-runs the P1-engine probe's own scripts (copied into rerun/) one Julia process at a time.
# ABOUTME: Each run takes the shared Julia lock; logs go to skeptic_P1-engine/logs/.
P=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
S=$P/proto/skeptic_P1-engine
run() {  # script log
    local t0=$(date +%s)
    flock -w 14400 $P/proto/julia.lock julia --startup-file=no --heap-size-hint=1500M $1 > $S/logs/$2.log 2>&1
    echo "$(date +%T) $2 exit $? ($(( $(date +%s) - t0 )) s incl. lock wait)" >> $S/logs/run_rerun.status
}
run $S/rerun/validate_engine.jl rerun_validate_engine
run $S/rerun/synthetic.jl rerun_synthetic
run $S/rerun/extras.jl rerun_extras
run $S/rerun/compare_lean.jl rerun_compare_lean
run $S/rerun/bfs_and_shift.jl rerun_bfs_and_shift
run $S/rerun/phantom_groups.jl rerun_phantom_groups
run $S/rerun/sizes_grouped.jl rerun_sizes_grouped
run $S/rerun/peek_deps.jl rerun_peek_deps
echo "$(date +%T) all done" >> $S/logs/run_rerun.status
