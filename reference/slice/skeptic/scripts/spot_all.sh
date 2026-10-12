#!/bin/bash
# ABOUTME: Runs the skeptic's independent spot checks on every population, 20 derived laws each,
# ABOUTME: one guarded Julia process per population under the shared lock; logs in skeptic/logs.
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
K=$S/slice/skeptic
cd $K
g() { log=$1; shift; LIMIT_MB=${LIMIT_MB:-3500} MAX_S=${MAX_S:-2400} $K/rerun/run_guarded.sh $log "$@"; tail -2 $log | head -1; }
N=20
g logs/spot_specs.log scripts/spot.jl specs $N specs
g logs/spot_handwritten.log scripts/spot.jl handwritten $N handwritten
g logs/spot_init_uniuni.log scripts/spot.jl init_uniuni $N init_uniuni
g logs/spot_init_unibi.log scripts/spot.jl init_unibi $N init_unibi
g logs/spot_init_bibi.log scripts/spot.jl init_bibi $N init_bibi
g logs/spot_init_bibi_pp.log scripts/spot.jl init_bibi_pp $N init_bibi_pp
g logs/spot_expand_bibi.log scripts/spot.jl expand_bibi $N expand_bibi
g logs/spot_allo_uniuni.log scripts/spot.jl allo_uniuni $N allo_uniuni
g logs/spot_allo_bibi.log scripts/spot.jl allo_bibi $N allo_bibi
g logs/spot_allo_uniuni_inh.log scripts/spot.jl allo_uniuni_inh $N allo_uniuni_inh
g logs/spot_allo_unibi.log scripts/spot.jl allo_unibi $N allo_unibi
g logs/spot_allo_bibi_pp.log scripts/spot.jl allo_bibi_pp $N allo_bibi_pp
g logs/spot_allo_terbi.log scripts/spot.jl allo_terbi $N allo_terbi
g logs/spot_p5a.log scripts/spot.jl p5a $N p5a_adv_cases p5a_adv_cases2 p5a_adv_cases3 p5a_adv_cases4 p5a_cases
g logs/spot_p1adv.log scripts/spot.jl p1adv $N p1_adv_pkg p1_adv_pkg2 p1_direct
g logs/spot_p1_enum_inh.log scripts/spot.jl p1_enum_inh $N p1_enum_inh
g logs/spot_p2adv.log scripts/spot.jl p2adv $N p2_adv_pkg p2_adv_pkg2 p2_fuzz_dump p2_fuzz_dump_flip p2_hand
g logs/spot_expand_terter.log scripts/spot.jl expand_terter $N expand_terter
g logs/spot_init_terter.log scripts/spot.jl init_terter $N init_terter
echo finished > logs/spot_all.done
