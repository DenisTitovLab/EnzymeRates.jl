#!/bin/bash
# ABOUTME: Second phase after the pipeline runs: ter-ter compile checks, key invariances, K1
# ABOUTME: classes, the ter-ter moves, the stage checks rerun on the integrated sources.
cd "$(dirname "$0")"
until [ -f logs/run_all_pipeline.done ]; do sleep 20; done
g() { log=$1; shift; LIMIT_MB=${LIMIT_MB:-3500} MAX_S=${MAX_S:-3600} ./run_guarded.sh $log "$@"; tail -2 $log | head -1; }
# parts 1-4 ran before the engine's cap counted parallel steps: rerun them
for k in 1 2 3 4; do g logs/pipe_terter$k.log run_pipeline.jl terter$k init_terter.part$k; done
g logs/select_terter.log select_terter.jl
for b in $(seq 1 9); do g logs/bench_terter$b.log bench.jl terter$b list:terter_bench_batch$b.jls $(for k in $(seq 1 11); do echo init_terter.part$k; done); done
g logs/bench_expand_terter.log bench.jl expand_terter all expand_terter
g logs/check_enum_terter.log check_enum_terter.jl
g logs/keys_core.log check_keys.jl core specs handwritten init_uniuni init_unibi init_bibi init_bibi_pp expand_bibi allo_uniuni allo_bibi allo_uniuni_inh allo_unibi allo_bibi_pp allo_terbi p5a_adv_cases p5a_adv_cases2 p5a_adv_cases3 p5a_adv_cases4 p5a_cases p1_adv_pkg p1_adv_pkg2 p1_enum_inh p1_direct p2_adv_pkg p2_adv_pkg2 p2_fuzz_dump p2_fuzz_dump_flip p2_hand
g logs/keys_expand_terter.log check_keys.jl expand_terter expand_terter
for k in $(seq 1 11); do g logs/keys_terter$k.log check_keys.jl terter$k init_terter.part$k; done
g logs/classes_core.log check_classes.jl core specs handwritten init_uniuni init_unibi init_bibi init_bibi_pp expand_bibi init_bibi_pp+expand_bibi allo_uniuni allo_bibi allo_uniuni_inh allo_unibi allo_bibi_pp allo_terbi p5a_adv_cases p5a_adv_cases2 p5a_adv_cases3 p5a_adv_cases4 p5a_cases p1_adv_pkg p1_adv_pkg2 p1_enum_inh p1_direct p2_adv_pkg p2_adv_pkg2 p2_fuzz_dump p2_fuzz_dump_flip p2_hand
LIMIT_MB=4500 g logs/classes_terter.log check_classes.jl terter init_terter.part1+init_terter.part2+init_terter.part3+init_terter.part4+init_terter.part5+init_terter.part6+init_terter.part7+init_terter.part8+init_terter.part9+init_terter.part10+init_terter.part11+expand_terter
g logs/stage_engine_check_kcat.log stage_engine_check_kcat.jl
g logs/stage_allo_check_kcat.log ../allo/check_kcat.jl
g logs/stage_allo_check_verdicts.log ../allo/check_verdicts.jl
g logs/bench_small.log bench.jl small all specs handwritten init_uniuni init_unibi p1_adv_pkg p1_adv_pkg2 p1_direct p2_adv_pkg p2_adv_pkg2 p2_fuzz_dump p2_fuzz_dump_flip p2_hand
g logs/peek_zero.log peek_zero.jl
g logs/check_bugs.log check_bugs.jl
echo finished > logs/run_phase2.done
