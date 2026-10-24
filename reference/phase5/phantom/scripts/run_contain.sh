#!/bin/bash
# ABOUTME: Runs contain.jl over the selection in resumable chunks of at most 600 seconds each,
# ABOUTME: through the research lock, until every selected pair has a row.
P=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45
cd $P/research5/phantom
for i in $(seq 1 12); do
  $P/jl_research.sh logs/contain_$i.log 600 --project=$P/../phase3/census/env --heap-size-hint=2500M \
      scripts/contain.jl data data/sel.tsv data/contain.tsv
  grep -q "stopping for time" logs/contain_$i.log || break
done
echo "all chunks done"
