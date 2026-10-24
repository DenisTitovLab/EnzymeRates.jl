#!/bin/bash
# ABOUTME: Runs the unit-rule census through the research lock in time-budgeted chunks (480 s of
# ABOUTME: work, 600 s cap), resuming each population where its last chunk stopped. Args: pops.
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
U=$S/phase45/research5/unit
for pop in "$@"; do
    n=$(wc -l < $U/data/$pop.tsv)
    i=$(ls $U/data/ | sed -n "s/^${pop}_[0-9]*_\([0-9]*\)\.jls$/\1/p" | sort -n | tail -1)
    i=$(( ${i:-0} + 1 ))
    while [ $i -le $n ]; do
        log=$U/logs/census_${pop}_$i.log
        $S/phase45/jl_research.sh $log 600 --startup-file=no --heap-size-hint=2500M \
            $U/census.jl $pop $i $n 480
        next=$(sed -n 's/^NEXT //p' $log)
        [ -n "$next" ] && [ "$next" -gt "$i" ] || { echo "chunk $pop $i failed"; exit 1; }
        i=$next
    done
done
echo "census done"
