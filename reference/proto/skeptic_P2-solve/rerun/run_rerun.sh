#!/bin/bash
# ABOUTME: Reruns the P2 probe's stdlib validation scripts (copied, include paths made absolute).
# ABOUTME: Holds the shared Julia lock for the whole batch; one log per script in skeptic logs/.
cd /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/rerun
for s in smoke census determinism correctness extra specs_compare ident ident_zf ident_hw17 pins nontext haldane_count zf splitcheck; do
    t0=$(date +%s)
    julia --startup-file=no --heap-size-hint=1500M $s.jl > /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/logs/rerun_$s.log 2>&1
    echo "$s exit $? $(($(date +%s)-t0)) s"
done
