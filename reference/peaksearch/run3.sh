#!/bin/bash
# ABOUTME: Run 3 of the peak search: multi-seed searches with the least-cancellation turnover.
# ABOUTME: One Julia process under the machine-wide lock; the timeout covers the lock wait.
D=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
P=$D/peaksearch
echo "multiseed start $(date +%T)"
timeout -s KILL 5400 flock -w 14400 $D/proto/julia.lock julia --project=$P --heap-size-hint=1500M $P/search2.jl > $P/multiseed.log 2>&1
echo "multiseed exit $? $(date +%T)"
