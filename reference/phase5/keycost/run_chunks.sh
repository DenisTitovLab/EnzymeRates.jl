#!/bin/bash
# ABOUTME: Runs keytime.jl over an input TSV in resumable chunks (480 s of work, 600 s cap) through
# ABOUTME: the research lock. Args: <in.tsv> <out.tsv> <key|code> <log prefix> [julia version arg].
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
K=$S/phase45/research5/keycost
IN=$1; OUT=$2; MODE=$3; PRE=$4; VER=$5
i=1; c=1
while true; do
    log=$K/logs/${PRE}_$c.log
    $S/phase45/jl_research.sh $log 600 $VER --startup-file=no --heap-size-hint=2500M \
        $K/keytime.jl $IN $i 480 $OUT $MODE
    grep -q '^DONE' $log && { echo "done $PRE"; exit 0; }
    next=$(sed -n 's/^NEXT //p' $log)
    [ -n "$next" ] && [ "$next" -gt "$i" ] || { echo "chunk $PRE $c failed at $i"; exit 1; }
    i=$next; c=$((c + 1))
done
