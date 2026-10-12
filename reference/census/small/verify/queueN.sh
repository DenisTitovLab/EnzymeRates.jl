#!/bin/bash
# ABOUTME: Generic follow-on queue: waits until the status file named first contains "finished", then
# ABOUTME: runs each named script under the shared Julia lock, logging to logs/<script>.log.
V=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/small/verify
L=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/julia.lock
P=/home/denis.linux/.julia/dev/EnzymeRates
cd $V
waitfor=$1; mine=$2; shift 2
until grep -q finished logs/$waitfor 2>/dev/null; do sleep 5; done
for s in "$@"; do
  ( time timeout -s KILL 3600 flock -w 14400 $L julia --project=$P --heap-size-hint=2500M $s.jl ) > logs/$s.log 2>&1
  echo "$s exit $?" >> logs/$mine
done
echo finished >> logs/$mine
