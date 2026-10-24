#!/bin/bash
# ABOUTME: Summarizes the remark-scan logs: per log, hoistscan.py's counts, the rounds that ran,
# ABOUTME: whether a GC error appeared, the exit status and peak memory.
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
L=$S/phase45/research5/law
for f in "$@"; do
    h=$(python3 -I $S/phase3/gcprobe/hoistscan.py $f | sed -n 2p | sed 's/^ *//')
    d=$(grep -c "^derivation full gc" $f); i=$(grep -c "^identify full gc" $f)
    printf "%s | %s | derivation rounds %s, identify rounds %s | GC error lines %s | %s | %s\n" \
        $(basename $f) "$h" $d $i $(grep -c "GC error" $f) "$(grep '^exit' $f)" "$(grep maxrss $f)"
done
