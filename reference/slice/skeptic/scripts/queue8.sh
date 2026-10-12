#!/bin/bash
# ABOUTME: Queues the package-side reproduction of the init_terter/77789 expansion failure.
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
K=$S/slice/skeptic
cd $K
until [ -f logs/queue7.done ]; do sleep 10; done
( exec 9>$S/proto/julia.lock; flock -w 14400 9 && timeout -s KILL 1800 julia --project=$S/er_oracle --heap-size-hint=2500M scripts/pkg_expand77789.jl > logs/pkg_expand77789.log 2>&1 9>&- )
echo finished > logs/queue8.done
