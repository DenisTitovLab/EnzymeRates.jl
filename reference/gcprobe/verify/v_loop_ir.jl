# ABOUTME: Prints the optimized LLVM IR of v_loop.jl's main to locate the box allocations
# ABOUTME: relative to the loop (preheader vs body).
include(joinpath(@__DIR__, "v_common.jl"))
function main(rounds, flag::Bool)
    for r in 1:rounds
        println(r, " ", nterms(make(flag, 1000).D)); flush(stdout)
        GC.gc(true)
        println("full gc ", r, " ok"); flush(stdout)
    end
    println("finished")
end
using InteractiveUtils
code_llvm(stdout, main, (Int, Bool); optimize = true, debuginfo = :none)
