# ABOUTME: Compiles gc6's main_b (optimized) so that LLVM's -print-changed output, set through
# ABOUTME: JULIA_LLVM_ARGS, shows which pass moves the Poly box allocation out of the loop.
ENV["IR"] = "1"
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/gcprobe/gc6.jl")
using InteractiveUtils
code_llvm(devnull, main_b, (); debuginfo = :none)
