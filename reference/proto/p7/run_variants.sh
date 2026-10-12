#!/bin/bash
# ABOUTME: Runs the P7 benchmark for every variant (one fresh Julia process each) under the shared lock.
# ABOUTME: Usage: run_variants.sh "variant:mode ..."
cd "$(dirname "$0")"
for vm in $1; do
  v=${vm%%:*}; m=${vm##*:}
  flock -w 14400 ../julia.lock julia --startup-file=no --heap-size-hint=1500M bench.jl $v $m > logs/bench_${v}_${m}.log 2>&1
  echo "finished $v $m: $(tail -1 logs/bench_${v}_${m}.log)"
done
