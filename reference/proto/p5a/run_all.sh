#!/bin/bash
# ABOUTME: Reruns every P5a study in sequence through the shared Julia lock.
# ABOUTME: Logs go to logs/; data records to data/.
cd "$(dirname "$0")"
J="flock -w 14400 ../julia.lock julia --startup-file=no --heap-size-hint=1500M"
$J run_cases.jl cases > logs/run_cases.log 2>&1
$J analyze_cases.jl cases > logs/analyze_cases.log 2>&1
$J run_cases.jl allo_uniuni > logs/run_allo_uniuni.log 2>&1
$J run_cases.jl allo_bibi_KV > logs/run_allo_bibi_KV.log 2>&1
$J details.jl > logs/details.log 2>&1
$J run_pop.jl expand_bibi > logs/pop_expand_bibi.log 2>&1
$J run_pop.jl allo_bibi > logs/pop_allo_bibi.log 2>&1
$J pop_analyze.jl expand_bibi i 2 > logs/pop_analyze_expand_i.log 2>&1
$J pop_analyze.jl expand_bibi iu 2 > logs/pop_analyze_expand_iu.log 2>&1
$J pop_analyze.jl allo_bibi i 2 > logs/pop_analyze_allo_i.log 2>&1
$J summary.jl > logs/summary.log 2>&1
echo done > logs/run_all.done
