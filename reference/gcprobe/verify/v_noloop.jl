# ABOUTME: Same body as v_loop.jl but unrolled into straight-line code (no loop for LICM to
# ABOUTME: hoist from); everything else, including the Union box and full GCs, is identical.
include(joinpath(@__DIR__, "v_common.jl"))
function main(flag::Bool)
    println(1, " ", nterms(make(flag, 1000).D)); flush(stdout); GC.gc(true); println("full gc 1 ok"); flush(stdout)
    println(2, " ", nterms(make(flag, 1000).D)); flush(stdout); GC.gc(true); println("full gc 2 ok"); flush(stdout)
    println(3, " ", nterms(make(flag, 1000).D)); flush(stdout); GC.gc(true); println("full gc 3 ok"); flush(stdout)
    println("finished")
end
main(false)
