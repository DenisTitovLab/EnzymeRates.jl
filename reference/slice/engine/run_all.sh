#!/bin/bash
# ABOUTME: Reruns every engine-stage check script under the shared Julia lock, one at a time,
# ABOUTME: writing each log to logs/ (the package export needs the er_oracle project).
cd /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/engine
LOCK=../../proto/julia.lock
run() { ( time timeout -s KILL 3000 flock -w 14400 $LOCK julia --startup-file=no --heap-size-hint=1500M "$@" ) ; }
[ -f data/pkg_kcat_times.jls ] || ( timeout -s KILL 3000 flock -w 14400 $LOCK julia --project=../../er_oracle --heap-size-hint=2500M export_pkg.jl > logs/export_pkg.log 2>&1 )
run check_engine.jl > logs/check_engine.log 2>&1
rm -f data/derived.jls
run check_determinism.jl > logs/check_determinism.log 2>&1
run check_determinism.jl > logs/check_determinism2.log 2>&1
run check_kcat.jl > logs/check_kcat.log 2>&1
run check_cost.jl > logs/check_cost.log 2>&1
run check_extra.jl > logs/check_extra.log 2>&1
run peek_inactive.jl > logs/peek_inactive.log 2>&1
run show_laws.jl specs/1 specs/5 init_bibi/55 specs/7 adv2/1 specs/11 handwritten/21 > logs/show_laws.log 2>&1
echo finished > logs/run_all.done
