#!/bin/bash
# ABOUTME: Exports every population after specs, one Julia run per population or chunk.
# ABOUTME: Each run takes the shared Julia lock; logs go to common/logs/.
P=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
C=$P/proto/common
run_pkg() {  # population chunk nchunks
    local log=$C/logs/$1${3:+.part$2}.log
    flock -w 14400 $P/proto/julia.lock julia --startup-file=no --project=$P/er_oracle \
        --heap-size-hint=2500M $C/export/export_pop.jl "$@" > $log 2>&1
    echo "$(date +%T) done $* (exit $?)"
}
run_merge() {
    flock -w 14400 $P/proto/julia.lock julia --startup-file=no --heap-size-hint=1500M \
        $C/export/merge_parts.jl "$@" >> $C/logs/$1.merge.log 2>&1
    echo "$(date +%T) merged $*"
}
for pop in specs init_uniuni init_unibi init_bibi init_bibi_pp handwritten allo_uniuni; do
    run_pkg $pop
done
run_pkg expand_bibi 1 2; run_pkg expand_bibi 2 2; run_merge expand_bibi 2
run_pkg allo_bibi 1 2; run_pkg allo_bibi 2 2; run_merge allo_bibi 2
echo "$(date +%T) all done"
