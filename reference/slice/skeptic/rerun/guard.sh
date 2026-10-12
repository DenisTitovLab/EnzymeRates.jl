#!/bin/bash
# ABOUTME: Runs a Julia script under the shared lock with a memory guard: kills it when its
# ABOUTME: resident set exceeds LIMIT_MB (default 3000); usage: guard.sh log script args...
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
LOG=$1; shift
LIMIT=${LIMIT_MB:-3000}
exec 9>$S/proto/julia.lock
flock -w 14400 9 || exit 1
julia --startup-file=no --heap-size-hint=1500M "$@" > "$LOG" 2>&1 9>&- &
PID=$!
while kill -0 $PID 2>/dev/null; do
  RSS=$(ps -o rss= -p $PID 2>/dev/null | tr -d ' ')
  if [ -n "$RSS" ] && [ "$RSS" -gt $((LIMIT * 1024)) ]; then
    echo "GUARD: killed at ${RSS} KB" >> "$LOG"; kill -9 $PID; break
  fi
  sleep 1
done
wait $PID 2>/dev/null
echo "exit $?" >> "$LOG"
