# ABOUTME: Prints the typed code and unoptimized LLVM IR of a reproducer's main (script in
# ABOUTME: MINI), to compare its union handling with the crashing gc6 main_b.
ENV["IR"] = "1"
include(ENV["MINI"])
using InteractiveUtils
sig = applicable(main) ? () : (Int,)
show(stdout, MIME"text/plain"(), first(only(code_typed(main, sig)))); println()
code_llvm(stdout, main, sig; debuginfo = :none, optimize = false)
