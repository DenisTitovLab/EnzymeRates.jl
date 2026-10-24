#!/bin/bash
# ABOUTME: Runs analyze.jl over one reaction's records in resumable 540-second chunks until every
# ABOUTME: record has a row; usage: run_analyze.sh <reaction> <rows> (one Julia at a time via the lock).
P=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45
cd $P/research5/phantom
for i in 1 2 3 4 5 6; do
  $P/jl_research.sh logs/analyze_$1_$i.log 600 --project=$P/../phase3/census/env --heap-size-hint=2500M \
      scripts/analyze.jl data $1 2 $2 data/facts_$1.tsv 1
  grep -q "stopping for time" logs/analyze_$1_$i.log || break
done
echo "all chunks done"
