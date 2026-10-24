#!/bin/bash
# ABOUTME: Runs crash_c.jl for law types C and D on Julia 1.12 and 1.13 with JuliaLICM remarks,
# ABOUTME: one Julia at a time; the logs keep the exit code (134 = GC abort).
export JULIA_LLVM_ARGS="-pass-remarks=julia-licm -pass-remarks-analysis=julia-licm"
for ver in 112 113; do for var in C D; do
  v=""; [ "$ver" = "113" ] && v="+1.13"
  /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/jl_research.sh /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research5/verify/logs/crash_${var}_${ver}.log 300 $v --startup-file=no --heap-size-hint=2000M /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research5/verify/crash_c.jl $var
done; done
