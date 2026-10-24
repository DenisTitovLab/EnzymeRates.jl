# ABOUTME: Dumps the unoptimized LLVM IR of gc6's main_b, to see where codegen itself places
# ABOUTME: the box for the Union{Poly{Int}, Poly{BigInt}} phi before any LLVM pass runs.
ENV["IR"] = "1"
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/gcprobe/gc6.jl")
using InteractiveUtils
code_llvm(stdout, main_b, (); debuginfo = :none, optimize = false)
