#!/bin/bash
# ABOUTME: Compiles and checks the generated code of every derived law of the core sets, one
# ABOUTME: guarded Julia process per batch; logs in logs/bench_<label>.log.
cd "$(dirname "$0")"
run() { label=$1; shift; LIMIT_MB=3500 timeout -s KILL 3600 ./guard.sh logs/bench_$label.log bench.jl $label all "$@"; tail -2 logs/bench_$label.log | head -1; }
run small specs handwritten init_uniuni init_unibi p1_adv_pkg p1_adv_pkg2 p1_direct p2_adv_pkg p2_adv_pkg2 p2_fuzz_dump p2_fuzz_dump_flip p2_hand
run init init_bibi init_bibi_pp
run expand expand_bibi
run enuminh p1_enum_inh
run allobibi allo_bibi
run allo allo_uniuni allo_uniuni_inh allo_unibi allo_bibi_pp allo_terbi p5a_adv_cases p5a_adv_cases2 p5a_adv_cases3 p5a_adv_cases4 p5a_cases
echo finished > logs/run_bench_core.done
