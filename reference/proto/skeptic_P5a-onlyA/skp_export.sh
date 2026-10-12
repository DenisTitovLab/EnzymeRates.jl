#!/bin/bash
# ABOUTME: Skeptic P5a: runs skp_export.jl for each extra population through the shared Julia lock.
# ABOUTME: Logs go to logs/export_<pop>.log.
cd "$(dirname "$0")"
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
for pop in allo_uniuni_inh allo_bibi_pp allo_unibi allo_terbi; do
    echo "$(date +%T) start $pop" >> logs/export.progress
    flock -w 14400 ../julia.lock julia --project=$S/er_oracle --heap-size-hint=2500M skp_export.jl $pop 400 > logs/export_$pop.log 2>&1
    echo "$(date +%T) end $pop rc=$?" >> logs/export.progress
done
echo done > logs/export.done
