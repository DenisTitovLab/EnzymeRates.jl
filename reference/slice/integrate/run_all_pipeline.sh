#!/bin/bash
# ABOUTME: Runs the pipeline driver on every population, one guarded Julia process per batch
# ABOUTME: of sets; logs in logs/pipe_<label>.log, rows in data/pipe_<label>.jls.
cd "$(dirname "$0")"
run() { label=$1; shift; LIMIT_MB=3500 timeout -s KILL 3600 ./guard.sh logs/pipe_$label.log run_pipeline.jl $label "$@"; tail -2 logs/pipe_$label.log | head -1; }
run small specs handwritten init_uniuni init_unibi p1_adv_pkg p1_adv_pkg2 p1_direct p2_adv_pkg p2_adv_pkg2 p2_fuzz_dump p2_fuzz_dump_flip p2_hand
run init init_bibi init_bibi_pp
run expand expand_bibi
run enuminh p1_enum_inh
run allobibi allo_bibi
run allo allo_uniuni allo_uniuni_inh allo_unibi allo_bibi_pp allo_terbi p5a_adv_cases p5a_adv_cases2 p5a_adv_cases3 p5a_adv_cases4 p5a_cases
run expand_terter expand_terter
for k in $(seq 1 11); do run terter$k init_terter.part$k; done
echo finished > logs/run_all_pipeline.done
