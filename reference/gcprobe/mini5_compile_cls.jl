# ABOUTME: Compiles the reproducer's main without running it, under the JuliaLICM detector, then
# ABOUTME: classifies the hoisted allocations, to check the classifier names P{Int} and P{BigInt}.
ENV["IR"] = "1"
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/gcprobe/mini5.jl")
precompile(main, (Int, Bool)) || error("precompile failed")
using InteractiveUtils; code_llvm(devnull, main, (Int, Bool))
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/gcprobe/classify.jl")
