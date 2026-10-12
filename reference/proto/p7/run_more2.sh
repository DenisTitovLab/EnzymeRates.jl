#!/bin/bash
# ABOUTME: Third batch of P7 runs: accuracy, type size, Cleland-form sizes, report.
cd "$(dirname "$0")"
L="flock -w 14400 ../julia.lock julia --startup-file=no --heap-size-hint=1500M"
$L typesize.jl > logs/typesize.log 2>&1; echo "typesize done"
$L cleland.jl > logs/cleland.log 2>&1; echo "cleland done"
$L report.jl > logs/report.log 2>&1; echo "report done"
$L accuracy.jl > logs/accuracy.log 2>&1; echo "accuracy done"
