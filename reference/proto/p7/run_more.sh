#!/bin/bash
# ABOUTME: Second batch of P7 runs (Haldane-only factoring with Horner, sharing, accuracy, type size).
cd "$(dirname "$0")"
L="flock -w 14400 ../julia.lock julia --startup-file=no --heap-size-hint=1500M"
$L bench.jl halh unique > logs/bench_halh_unique.log 2>&1; echo "halh: $(tail -1 logs/bench_halh_unique.log)"
$L bench.jl hal unique > logs/bench_hal_unique.log 2>&1; echo "hal: $(tail -1 logs/bench_hal_unique.log)"
$L share.jl facth > logs/share_facth.log 2>&1; echo "share facth: $(tail -2 logs/share_facth.log)"
$L share.jl flat > logs/share_flat.log 2>&1; echo "share flat: $(tail -2 logs/share_flat.log)"
$L accuracy.jl > logs/accuracy.log 2>&1; echo "accuracy: $(tail -1 logs/accuracy.log)"
$L typesize.jl > logs/typesize.log 2>&1; echo "typesize: $(tail -1 logs/typesize.log)"
