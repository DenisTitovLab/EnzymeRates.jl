#!/bin/bash
# ABOUTME: Reruns every solve-stage check script under the shared Julia lock, one at a time,
# ABOUTME: writing each log to logs/; the order matches the report's checks (a)-(f).
cd /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/solve
for s in check_plain check_determinism check_rational check_against check_allo check_cost; do
    ( time timeout 3000 flock -w 14400 /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/julia.lock julia --startup-file=no --heap-size-hint=1500M $s.jl ) > logs/$s.log 2>&1
done
( time LAMBDA_EXP=60 NRANDOM=60 timeout 3000 flock -w 14400 /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/julia.lock julia --startup-file=no --heap-size-hint=1500M check_rank.jl ) > logs/check_rank.log 2>&1
echo finished > logs/run_all.done
