#!/bin/bash
# ABOUTME: Fourth batch of P7 runs: segment-only factoring (no blocks) check, benchmark and report.
cd "$(dirname "$0")"
L="flock -w 14400 ../julia.lock julia --startup-file=no --heap-size-hint=1500M"
BLOCKS=0 $L check.jl > logs/check_all_noblocks.log 2>&1; echo "check noblocks done"
$L bench.jl segh unique > logs/bench_segh_unique.log 2>&1; echo "segh: $(tail -1 logs/bench_segh_unique.log)"
$L report.jl > logs/report.log 2>&1; echo "report done"
