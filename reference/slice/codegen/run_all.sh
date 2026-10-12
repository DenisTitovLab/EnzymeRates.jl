#!/bin/bash
# ABOUTME: Reruns every codegen-stage check under the shared Julia lock, one process per batch
# ABOUTME: of sets (values, MWC values, benchmarks, the package's benchmarks, compile budget).
cd /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/codegen
LOCK=../../proto/julia.lock
run() { ( time timeout -s KILL 3000 flock -w 14400 $LOCK julia --startup-file=no --heap-size-hint=1500M "$@" ) ; }
pkg() { ( time timeout -s KILL 3000 flock -w 14400 $LOCK julia --project=../../er_oracle --heap-size-hint=2500M "$@" ) ; }
SMALL="specs handwritten init_uniuni init_unibi p1_adv_pkg p1_adv_pkg2 p1_direct p2_adv_pkg p2_adv_pkg2 p2_fuzz_dump p2_fuzz_dump_flip p2_hand"
ALLO="allo_uniuni allo_uniuni_inh allo_unibi allo_bibi_pp allo_terbi p5a_adv_cases p5a_adv_cases2 p5a_adv_cases3 p5a_adv_cases4 p5a_cases"
run check_values.jl $SMALL > logs/check_values_small.log 2>&1
run check_values.jl init_bibi init_bibi_pp > logs/check_values_init.log 2>&1
run check_values.jl expand_bibi > logs/check_values_expand.log 2>&1
run check_values.jl p1_enum_inh > logs/check_values_enum_inh.log 2>&1
run check_allo.jl specs handwritten expand_bibi > logs/check_allo_shared.log 2>&1
run check_allo.jl allo_bibi > logs/check_allo_bibi.log 2>&1
run check_allo.jl $ALLO > logs/check_allo_other.log 2>&1
run bench.jl small specs handwritten init_uniuni init_unibi > logs/bench_small.log 2>&1
run bench.jl init init_bibi init_bibi_pp > logs/bench_init.log 2>&1
run bench.jl expand expand_bibi > logs/bench_expand.log 2>&1
run bench.jl adv p1_adv_pkg p1_adv_pkg2 p1_direct p2_adv_pkg p2_adv_pkg2 p2_fuzz_dump p2_fuzz_dump_flip p2_hand > logs/bench_adv.log 2>&1
run bench.jl enuminh p1_enum_inh > logs/bench_enuminh.log 2>&1
run bench.jl allobibi allo_bibi > logs/bench_allobibi.log 2>&1
run bench.jl allo $ALLO > logs/bench_allo.log 2>&1
pkg bench_pkg.jl small specs handwritten init_uniuni init_unibi > logs/bench_pkg_small.log 2>&1
pkg bench_pkg.jl init init_bibi init_bibi_pp > logs/bench_pkg_init.log 2>&1
pkg bench_pkg.jl expand expand_bibi > logs/bench_pkg_expand.log 2>&1
pkg bench_pkg.jl allo allo_uniuni allo_bibi > logs/bench_pkg_allo.log 2>&1
pkg bench_pkg_inh.jl > logs/bench_pkg_inh.log 2>&1
run budget.jl > logs/budget.log 2>&1
rm -f data/body_hashes.jls
run check_determinism.jl > logs/check_determinism1.log 2>&1
run check_determinism.jl > logs/check_determinism2.log 2>&1
run summary.jl > logs/summary.log 2>&1
echo finished > logs/run_all.done
