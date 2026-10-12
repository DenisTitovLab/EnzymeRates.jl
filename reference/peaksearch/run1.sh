#!/bin/bash
# ABOUTME: Run 1 of the peak search: package cross-check of the solver, then a short pilot search.
# ABOUTME: Each Julia process takes the machine-wide lock; the timeout covers the lock wait.
D=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
L=$D/proto/julia.lock
P=$D/peaksearch
echo "crosscheck start $(date +%T)"
timeout -s KILL 4800 flock -w 14400 $L julia --project=/home/denis.linux/.julia/dev/EnzymeRates --heap-size-hint=2500M $P/crosscheck.jl > $P/crosscheck.log 2>&1
echo "crosscheck exit $? $(date +%T)"
timeout -s KILL 4800 flock -w 14400 $L julia --project=$P --heap-size-hint=1500M $P/search.jl pilot 20000 T,g,d1,f zero,Q 0 > $P/pilot.log 2>&1
echo "pilot exit $? $(date +%T)"
