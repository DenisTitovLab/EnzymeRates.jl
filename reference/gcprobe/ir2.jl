# ABOUTME: Dumps the optimized LLVM IR of lawT(Int, ...) and wplain, the callee and caller
# ABOUTME: of the crashing struct return, for comparison of how the sret result is rooted.
ENV["IR"] = "1"
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/gcprobe/gc3.jl")
using InteractiveUtils
code_llvm(stdout, lawT, (Type{Int}, S.Mechanism, S.Solution, BitVector, Bool); debuginfo = :none)
