#!/bin/bash
# ABOUTME: Skeptic's clean rerun, second batch: compile checks of the core sets and expand_terter,
# ABOUTME: key invariances and K1 classes of the core sets, the bug checks and two stage checks.
cd "$(dirname "$0")"
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
g() { log=$1; shift; LIMIT_MB=${LIMIT_MB:-3500} MAX_S=${MAX_S:-3600} ./run_guarded.sh $log "$@"; tail -2 $log | head -1; }
g logs/bench_small.log bench.jl small all specs handwritten init_uniuni init_unibi p1_adv_pkg p1_adv_pkg2 p1_direct p2_adv_pkg p2_adv_pkg2 p2_fuzz_dump p2_fuzz_dump_flip p2_hand
g logs/bench_init.log bench.jl init all init_bibi init_bibi_pp
g logs/bench_expand.log bench.jl expand all expand_bibi
g logs/bench_enuminh.log bench.jl enuminh all p1_enum_inh
g logs/bench_allobibi.log bench.jl allobibi all allo_bibi
g logs/bench_allo.log bench.jl allo all allo_uniuni allo_uniuni_inh allo_unibi allo_bibi_pp allo_terbi p5a_adv_cases p5a_adv_cases2 p5a_adv_cases3 p5a_adv_cases4 p5a_cases
g logs/bench_expand_terter.log bench.jl expand_terter all expand_terter
g logs/keys_core.log check_keys.jl core specs handwritten init_uniuni init_unibi init_bibi init_bibi_pp expand_bibi allo_uniuni allo_bibi allo_uniuni_inh allo_unibi allo_bibi_pp allo_terbi p5a_adv_cases p5a_adv_cases2 p5a_adv_cases3 p5a_adv_cases4 p5a_cases p1_adv_pkg p1_adv_pkg2 p1_enum_inh p1_direct p2_adv_pkg p2_adv_pkg2 p2_fuzz_dump p2_fuzz_dump_flip p2_hand
g logs/classes_core.log check_classes.jl core specs handwritten init_uniuni init_unibi init_bibi init_bibi_pp expand_bibi init_bibi_pp+expand_bibi allo_uniuni allo_bibi allo_uniuni_inh allo_unibi allo_bibi_pp allo_terbi p5a_adv_cases p5a_adv_cases2 p5a_adv_cases3 p5a_adv_cases4 p5a_cases p1_adv_pkg p1_adv_pkg2 p1_enum_inh p1_direct p2_adv_pkg p2_adv_pkg2 p2_fuzz_dump p2_fuzz_dump_flip p2_hand
g logs/check_bugs.log check_bugs.jl
g logs/stage_engine_check_kcat.log stage_engine_check_kcat.jl
g logs/check_enum_terter.log check_enum_terter.jl
g logs/compare_regress_before.log compare_regress.jl regress_before.jls regress_skeptic.jls
echo finished > logs/rerun_B.done
