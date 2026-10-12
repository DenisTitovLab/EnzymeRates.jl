#!/bin/bash
# ABOUTME: Runs the census analyses of the saved runs, one locked Julia process at a time,
# ABOUTME: each with a hard timeout; logs go to ../logs/<step>.log.
D=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/large
L=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/julia.lock
P=/home/denis.linux/.julia/dev/EnzymeRates
R="ter_bi_rxn ter_ter_rxn pyruvate_carboxylase_rxn"
run() {  # run <log> <timeout> <script> [args...]
    local log=$D/logs/$1.log t=$2; shift 2
    echo "start $(date +%T)" > $log
    (time timeout -s KILL $t flock -w 14400 $L julia --project=$P --heap-size-hint=2500M $D/scripts/"$@") >> $log 2>&1
    echo "exit $? $(date +%T)" >> $log
}
run report_all 900 report_saved.jl $R
run moves 600 moves.jl $R
run examples 600 examples.jl 3 $R
for r in $R; do run samelaw_$r 900 samelaw.jl $r 12; done
run costs 1500 costs.jl 120 $R
run degenerate 900 degenerate.jl $R
echo "all done $(date +%T)" >> $D/logs/analyze.done
