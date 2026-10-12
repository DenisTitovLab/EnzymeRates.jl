#!/bin/bash
# ABOUTME: Skeptic's clean rerun of the integration scripts (copies in this directory), one
# ABOUTME: guarded Julia process each under the shared lock; outputs stay in this directory.
cd "$(dirname "$0")"
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
g() { log=$1; shift; LIMIT_MB=${LIMIT_MB:-3500} MAX_S=${MAX_S:-3600} ./run_guarded.sh $log "$@"; tail -2 $log | head -1; }
g logs/pipe_small.log run_pipeline.jl small specs handwritten init_uniuni init_unibi p1_adv_pkg p1_adv_pkg2 p1_direct p2_adv_pkg p2_adv_pkg2 p2_fuzz_dump p2_fuzz_dump_flip p2_hand
g logs/pipe_init.log run_pipeline.jl init init_bibi init_bibi_pp
g logs/pipe_expand.log run_pipeline.jl expand expand_bibi
g logs/pipe_enuminh.log run_pipeline.jl enuminh p1_enum_inh
g logs/pipe_allobibi.log run_pipeline.jl allobibi allo_bibi
g logs/pipe_allo.log run_pipeline.jl allo allo_uniuni allo_uniuni_inh allo_unibi allo_bibi_pp allo_terbi p5a_adv_cases p5a_adv_cases2 p5a_adv_cases3 p5a_adv_cases4 p5a_cases
g logs/pipe_expand_terter.log run_pipeline.jl expand_terter expand_terter
g logs/pipe_terter11.log run_pipeline.jl terter11 init_terter.part11
g logs/pipe_terter1.log run_pipeline.jl terter1 init_terter.part1
( exec 9>$S/proto/julia.lock; flock -w 14400 9 && timeout -s KILL 1800 julia --startup-file=no --heap-size-hint=1500M regress.jl $S/slice/src data/regress_skeptic.jls > logs/regress_skeptic.log 2>&1 9>&- ); tail -1 logs/regress_skeptic.log
g logs/compare_regress.log compare_regress.jl regress_after2.jls regress_skeptic.jls
echo finished > logs/rerun_A.done
