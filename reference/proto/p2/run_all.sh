#!/bin/bash
# ABOUTME: Runs every stdlib-only P2 check in sequence under the shared Julia lock, one log each.
# ABOUTME: Usage: bash run_all.sh (from proto/p2).
cd "$(dirname "$0")"
for s in smoke census determinism correctness extra specs_compare ident; do
    flock -w 14400 ../julia.lock julia --startup-file=no --heap-size-hint=1500M $s.jl > logs/$s.log 2>&1
    echo "$s exit $?"
done
