#!/bin/bash
# ABOUTME: Third phase after the second: today's package on the slowest ter-ter laws, the worst
# ABOUTME: oracle errors at a faster rapid-equilibrium scale, the result tables, line counts.
cd "$(dirname "$0")"
until [ -f logs/run_phase2.done ]; do sleep 20; done
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
g() { log=$1; shift; LIMIT_MB=${LIMIT_MB:-3500} MAX_S=${MAX_S:-3600} ./run_guarded.sh $log "$@"; tail -2 $log | head -1; }
( exec 9>$S/proto/julia.lock; flock -w 14400 9 && timeout -s KILL 3000 julia --project=$S/er_oracle --heap-size-hint=2500M pkg_bench_terter.jl > logs/pkg_bench_terter.log 2>&1 9>&- )
g logs/check_oracle_worst.log check_oracle_worst.jl
g logs/final_summary.log final_summary.jl
python3 linecounts.py > logs/linecounts.log 2>&1
echo finished > logs/run_phase3.done
