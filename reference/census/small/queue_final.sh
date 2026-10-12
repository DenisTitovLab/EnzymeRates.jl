#!/bin/bash
# ABOUTME: Runs the final census summaries one after another, each under the shared Julia lock
# ABOUTME: with a hard timeout: tables, :OnlyA examples, and the report of the supplementary runs.
CEN=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/small
L=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/julia.lock
P=/home/denis.linux/.julia/dev/EnzymeRates
cd $CEN
ALL=uni_uni_rxn,uni_uni_with_reg,uni_uni_allo,uni_uni_allo_reg,uni_uni_allo_2reg,uni_uni_reg_and_inhibitor,uni_bi_rxn,bi_bi_rxn,bi_bi_pp_rxn,bi_bi_rxn.r2_100000,bi_bi_pp_rxn.r2_1200
( time timeout -s KILL 10800 flock -w 14400 $L julia --project=$P --heap-size-hint=2500M tables.jl $ALL ) > logs/tables.log 2>&1
( time timeout -s KILL 10800 flock -w 14400 $L julia --project=$P --heap-size-hint=2500M wexamples.jl bi_bi_pp_rxn,bi_bi_pp_rxn.r2_1200,bi_bi_rxn.r2_100000 ) > logs/wexamples.log 2>&1
( time timeout -s KILL 10800 flock -w 14400 $L julia --project=$P --heap-size-hint=2500M report.jl bi_bi_rxn.r2_100000,bi_bi_pp_rxn.r2_1200 ) > logs/report_supp.log 2>&1
echo finished > logs/queue_final.done
