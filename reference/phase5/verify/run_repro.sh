#!/bin/bash
# ABOUTME: Runs repro_c.jl in each mode on Julia 1.12 and 1.13 with JuliaLICM remarks, one Julia at
# ABOUTME: a time through the research lock; the logs keep the exit code (134 = GC abort).
export JULIA_LLVM_ARGS="-pass-remarks=julia-licm -pass-remarks-analysis=julia-licm"
for ver in 112 113; do for mode in vec nothing alias; do
  v=""; [ "$ver" = "113" ] && v="+1.13"
  /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/jl_research.sh /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research5/verify/logs/repro_${mode}_${ver}.log 120 $v --startup-file=no /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/research5/verify/repro_c.jl $mode
done; done
