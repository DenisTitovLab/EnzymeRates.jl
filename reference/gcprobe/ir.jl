# ABOUTME: Dumps the typed code and optimized LLVM IR of the crashing wrapper `wplain` and of
# ABOUTME: `lawT(Int, ...)`'s signature, to see how the struct result is boxed for the Union.
ENV["IR"] = "1"
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/gcprobe/gc3.jl")
using InteractiveUtils
println(Base.return_types(wplain, (S.Mechanism, S.Solution, BitVector, Bool)))
code_llvm(stdout, wplain, (S.Mechanism, S.Solution, BitVector, Bool); debuginfo = :none)
