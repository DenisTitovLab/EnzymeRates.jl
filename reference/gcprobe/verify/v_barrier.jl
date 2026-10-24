# ABOUTME: v_loop.jl with the field read moved behind a @noinline function barrier, so the
# ABOUTME: caller's loop never boxes a Union{P{Int}, P{Float64}} (the report's mitigation).
include(joinpath(@__DIR__, "v_common.jl"))
@noinline dterms(x) = nterms(x.D)
function main(rounds, flag::Bool)
    for r in 1:rounds
        println(r, " ", dterms(make(flag, 1000))); flush(stdout)
        GC.gc(true)
        println("full gc ", r, " ok"); flush(stdout)
    end
    println("finished")
end
main(3, false)
