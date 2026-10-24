# ABOUTME: Dumps the optimized LLVM IR of gc6's main_b, the smallest crashing caller of the
# ABOUTME: Union{Law{Int}, Law{BigInt}}-returning wdirect.
ENV["IR"] = "1"
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/gcprobe/gc6.jl")
using InteractiveUtils
code_typed(main_b, ()) |> only |> first |> x -> show(stdout, MIME"text/plain"(), x)
println()
code_llvm(stdout, main_b, (); debuginfo = :none)
