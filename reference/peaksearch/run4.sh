#!/bin/bash
# ABOUTME: Run 4 of the peak search: package rate_equation evaluated at the optima and examples.
# ABOUTME: One Julia process under the machine-wide lock; the timeout covers the lock wait.
D=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
P=$D/peaksearch
echo "best start $(date +%T)"
timeout -s KILL 3600 flock -w 14400 $D/proto/julia.lock julia --project=/home/denis.linux/.julia/dev/EnzymeRates --heap-size-hint=2500M $P/crosscheck_best.jl > $P/crosscheck_best.log 2>&1
echo "best exit $? $(date +%T)"
