#!/bin/bash
# ABOUTME: Runs one Julia job under the shared lock with a hard kill timeout; logs to $1.
# ABOUTME: Usage: run.sh <log> <timeout_s> <julia args...>  (JL env picks the julia binary/channel)
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
LOG=$1; T=$2; shift 2
START=$(date +%s)
timeout -s KILL $T flock -w 14400 $S/proto/julia.lock /usr/bin/time -f "maxrss_kb %M" julia ${JLCH:+$JLCH} "$@" > "$LOG" 2>&1
RC=$?
echo "exit $RC elapsed $(( $(date +%s) - START ))s" >> "$LOG"
echo "$LOG: exit $RC, GCerr=$(grep -c "GC error" $LOG), $(grep maxrss_kb $LOG), $(( $(date +%s) - START ))s"
