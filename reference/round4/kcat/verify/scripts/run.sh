#!/bin/bash
# ABOUTME: Runs one verifier script under the shared machine-wide Julia lock with a hard
# ABOUTME: KILL timeout; usage: run.sh <seconds> <log name> <script> [args...]
V=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat/verify
L=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/julia.lock
T=$1; LOG=$V/logs/$2.log; shift 2
echo "start $(date +%T)" > "$LOG"
timeout -s KILL "$T" flock -w 14400 "$L" julia --project=$V/env --heap-size-hint=1500M "$V/scripts/$@" >> "$LOG" 2>&1
echo "exit $? $(date +%T)" >> "$LOG"
