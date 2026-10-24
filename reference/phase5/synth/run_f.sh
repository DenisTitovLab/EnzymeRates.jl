#!/bin/bash
# ABOUTME: Runs variant F's agreement and timing check, then the typed GC workload for C, E and F on Julia
# ABOUTME: 1.12 and 1.13 with JuliaLICM remarks and --pkgimages=no, one Julia at a time; hoistscan at the end.
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
D=$S/phase45/research5/synth
J=$S/phase45/jl_research.sh
# agreement ran separately (logs/agree_f.log)
export JULIA_LLVM_ARGS="-pass-remarks=julia-licm -pass-remarks-analysis=julia-licm"
for ver in 112 113; do for var in F E C; do
  v=""; [ "$ver" = "113" ] && v="+1.13"
  $J $D/logs/scan_${var}_${ver}.log 600 $v --startup-file=no --pkgimages=no --heap-size-hint=2000M $D/scan_f.jl $var
done; done
unset JULIA_LLVM_ARGS
cd $D/logs && python3 -I $S/phase3/gcprobe/hoistscan.py scan_*.log > hoistscan.log 2>&1
echo "exit 0" >> $D/logs/run_f.done
