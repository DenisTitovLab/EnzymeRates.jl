# ABOUTME: repro_upstream.jl with flushed output, so the log shows which round aborts.
# ABOUTME: Same loop: `make(flag, n).D` is a Union{P{Int}, P{Float64}} read inside a loop.
include(joinpath(@__DIR__, "v_common.jl"))
function main(rounds, flag::Bool)
    for r in 1:rounds
        println(r, " ", nterms(make(flag, 1000).D)); flush(stdout)
        GC.gc(true)
        println("full gc ", r, " ok"); flush(stdout)
    end
    println("finished")
end
main(3, false)
