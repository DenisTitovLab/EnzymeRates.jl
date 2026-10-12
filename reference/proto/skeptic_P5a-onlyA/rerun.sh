#!/bin/bash
# ABOUTME: Skeptic P5a: reruns the probe's own validation and study scripts (copied here) so the
# ABOUTME: outputs land in this directory; every Julia process goes through the shared lock.
cd "$(dirname "$0")"
J="flock -w 14400 ../julia.lock julia --startup-file=no --heap-size-hint=1500M"
for step in "validate.jl:validate" "run_cases.jl cases:run_cases" "analyze_cases.jl cases:analyze_cases" \
            "run_cases.jl allo_uniuni:run_allo_uniuni" "run_cases.jl allo_bibi_KV:run_allo_bibi_KV" \
            "run_pop.jl expand_bibi:pop_expand_bibi" "run_pop.jl allo_bibi:pop_allo_bibi" \
            "summary.jl:summary" "kcat.jl:kcat"; do
    cmd=${step%%:*}; log=${step##*:}
    echo "$(date +%T) start $cmd" >> logs/rerun.progress
    /usr/bin/time -v $J $cmd > logs/$log.log 2> logs/$log.err
    echo "$(date +%T) end $cmd rc=$?" >> logs/rerun.progress
done
echo done > logs/rerun.done
