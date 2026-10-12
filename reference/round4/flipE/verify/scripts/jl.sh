#!/bin/bash
# ABOUTME: Runs one Julia script under the machine-wide Julia lock with a hard kill timeout.
# ABOUTME: Usage: jl.sh <seconds> <heap> <project> <log> <script> [args...]
LOCK=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/julia.lock
T=$1; H=$2; P=$3; LOG=$4; shift 4
echo "start $(date +%T) project=$P" > "$LOG"
timeout -s KILL "$T" flock -w 14400 "$LOCK" julia --project="$P" --heap-size-hint="$H" "$@" >> "$LOG" 2>&1
echo "exit $? $(date +%T)" >> "$LOG"
