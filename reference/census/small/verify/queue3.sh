#!/bin/bash
# ABOUTME: Third verification queue: waits for the second to finish, then runs each script under the
# ABOUTME: shared Julia lock with a hard timeout, logging to logs/<script>.log.
V=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/small/verify
L=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/julia.lock
P=/home/denis.linux/.julia/dev/EnzymeRates
cd $V
until grep -q finished logs/queue2.status 2>/dev/null; do sleep 5; done
for s in "$@"; do
  ( time timeout -s KILL 3600 flock -w 14400 $L julia --project=$P --heap-size-hint=2500M $s.jl ) > logs/$s.log 2>&1
  echo "$s exit $?" >> logs/queue3.status
done
echo finished >> logs/queue3.status
