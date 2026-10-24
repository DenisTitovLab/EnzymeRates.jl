# ABOUTME: Prints the optimized LLVM IR of a reproducer's main (script path in MINI) to check
# ABOUTME: whether the union phi's box allocation sits outside the loop.
ENV["IR"] = "1"
include(ENV["MINI"])
using InteractiveUtils
code_llvm(stdout, main, (); debuginfo = :none)
