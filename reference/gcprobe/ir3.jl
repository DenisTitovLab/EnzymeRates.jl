# ABOUTME: Dumps the optimized LLVM IR of Slice.law(Int, ...), the function that fills the
# ABOUTME: struct result through sret, to see when its stores to the sret buffer happen.
ENV["IR"] = "1"
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/gcprobe/gc3.jl")
using InteractiveUtils
code_llvm(stdout, S.law, (Type{Int}, S.Mechanism, S.Solution, BitVector, Bool); debuginfo = :none)
