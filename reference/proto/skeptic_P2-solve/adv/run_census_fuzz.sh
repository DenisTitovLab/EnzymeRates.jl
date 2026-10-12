#!/bin/bash
# ABOUTME: Runs the skeptic's rank census and fuzzer under one hold of the shared Julia lock.
# ABOUTME: Logs go to skeptic_P2-solve/logs/.
t0=$(date +%s)
julia --startup-file=no --heap-size-hint=1500M /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/adv/fuzz.jl 40 > /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/logs/fuzz_smoke.log 2>&1
echo "fuzz smoke exit $? $(($(date +%s)-t0)) s"
t0=$(date +%s)
julia --startup-file=no --heap-size-hint=1500M /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/adv/rank_census.jl > /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/logs/rank_census.log 2>&1
echo "rank_census exit $? $(($(date +%s)-t0)) s"
t0=$(date +%s)
julia --startup-file=no --heap-size-hint=1500M /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/adv/fuzz.jl 600 > /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/skeptic_P2-solve/logs/fuzz.log 2>&1
echo "fuzz exit $? $(($(date +%s)-t0)) s"
