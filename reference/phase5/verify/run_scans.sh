#!/bin/bash
# ABOUTME: Runs the verifier's GC workload for law types U (control) and C on Julia 1.12 and 1.13
# ABOUTME: with JuliaLICM remarks and --pkgimages=no, one at a time through the research lock.
export JULIA_LLVM_ARGS="-pass-remarks=julia-licm -pass-remarks-analysis=julia-licm"
for spec in "C 112" "U 112" "C 113" "U 113"; do
  set -- $spec
  ver=""; [ "$2" = "113" ] && ver="+1.13"
  /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/jl_research.sh /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research5/verify/logs/scan_$1_$2.log 600 $ver --startup-file=no --pkgimages=no --heap-size-hint=2000M /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research5/verify/scan_v.jl $1
done
