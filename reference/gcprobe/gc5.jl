# ABOUTME: Deterministic F6 trigger: one rate law through a non-inlined function returning
# ABOUTME: Union{Law{Int}, Law{BigInt}}, then a full GC; MODE picks how the caller uses it.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/src/Slice.jl")
const S = Slice
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/gcprobe/mechs_src.jl")
const MODE = get(ENV, "MODE", "split")

@noinline function wdirect(m, sol, mask, flag::Bool)
    flag ? S.law(BigInt, m, sol, mask, false) : S.law(Int, m, sol, mask, false)
end

function run_split(m)
    sol = S.solve(m)
    L = wdirect(m, sol, trues(length(m.edges)), false)
    println("nD ", S.nterms(L.D)); flush(stdout)
    nothing
end
function run_barrier(m)
    sol = S.solve(m)
    L = Base.inferencebarrier(wdirect(m, sol, trues(length(m.edges)), false))
    println("nD ", S.nterms(L.D)); flush(stdout)
    nothing
end
function run_drop(m)
    sol = S.solve(m)
    wdirect(m, sol, trues(length(m.edges)), false)
    println("dropped"); flush(stdout)
    nothing
end

function main()
    m = S.Mechanism(RAWS[1][1]...; RAWS[1][2]...)
    f = MODE == "split" ? run_split : MODE == "barrier" ? run_barrier : run_drop
    for i in 1:3
        f(m)
        println("full gc $i"); flush(stdout); GC.gc(true); println("ok $i"); flush(stdout)
    end
    println("finished")
end
get(ENV, "IR", "0") == "1" || main()
