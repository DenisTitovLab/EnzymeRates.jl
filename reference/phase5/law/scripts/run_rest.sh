#!/bin/bash
# ABOUTME: Runs the law study's timing chunks, the engine checks (Q9, Q10, Q14, bounds) and the Q11
# ABOUTME: package check, one Julia at a time through the research lock; logs under law/logs.
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
L=$S/phase45/research5/law
J=$S/phase45/jl_research.sh
for sets in spec,twosites,uniuni bibi pingpong terbi,prefix; do
    SETS=$sets TAG=${sets//,/_} $J $L/logs/timing_${sets//,/_}.log 600 --startup-file=no --heap-size-hint=2000M $L/scripts/timing.jl
done
$J $L/logs/checks.log 600 --startup-file=no --heap-size-hint=2000M $L/scripts/checks.jl
$J $L/logs/q11.log 600 --startup-file=no --project=$L/env --heap-size-hint=2000M $L/scripts/q11.jl
echo "rest done" > $L/logs/run_rest.done
