#!/bin/bash
# ABOUTME: Follow-up B1: runs full_eval.jl over every population in chunks through the shared Julia lock,
# ABOUTME: so other agents can take the lock between chunks. Logs to logs_full/.
cd "$(dirname "$0")"
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
J="flock -w 14400 $S/proto/julia.lock julia --startup-file=no --heap-size-hint=1500M"
mkdir -p logs_full
B1_MAXDIRS=32 $J full_eval.jl specs > logs_full/specs.log 2>&1
for p in allo_uniuni_inh allo_unibi allo_terbi; do $J full_eval.jl $p > logs_full/$p.log 2>&1; done
for r in "1 250" "251 460"; do set -- $r; $J full_eval.jl expand_bibi $1 $2 > logs_full/expand_bibi_$1.log 2>&1; done
for r in "1 200" "201 400"; do set -- $r; $J full_eval.jl allo_bibi_pp $1 $2 > logs_full/allo_bibi_pp_$1.log 2>&1; done
for r in "1 250" "251 500" "501 750" "751 1000" "1001 1250" "1251 1477"; do set -- $r; $J full_eval.jl allo_bibi $1 $2 > logs_full/allo_bibi_$1.log 2>&1; done
echo done > logs_full/run_full.done
