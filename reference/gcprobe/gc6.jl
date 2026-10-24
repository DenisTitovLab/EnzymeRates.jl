# ABOUTME: gc3's wdirect path alone (six mechanisms built, a full GC after each law), with
# ABOUTME: switches that remove one difference from the non-crashing gc5 at a time.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/src/Slice.jl")
const S = Slice
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase3/gcprobe/mechs_src.jl")
const NM = parse(Int, get(ENV, "NM", "6"))
const V = get(ENV, "V", "a")

@noinline function wdirect(m, sol, mask, flag::Bool)
    flag ? S.law(BigInt, m, sol, mask, false) : S.law(Int, m, sol, mask, false)
end

function main_a()   # as gc3: loop over ms, union-split result
    ms = [S.Mechanism(a...; kw...) for (a, kw) in RAWS[1:NM]]
    println("built ", length(ms)); flush(stdout)
    for r in 1:2, (k, m) in enumerate(ms)
        sol = S.solve(m)
        println(r, " ", k, " ", S.nterms(wdirect(m, sol, trues(length(m.edges)), false).D))
        println("full gc"); flush(stdout); GC.gc(true); println("full gc ok"); flush(stdout)
    end
    println("finished")
end
function main_b()   # one mechanism built in main, same loop body
    m = S.Mechanism(RAWS[1][1]...; RAWS[1][2]...)
    for r in 1:2
        sol = S.solve(m)
        println(r, " ", S.nterms(wdirect(m, sol, trues(length(m.edges)), false).D))
        println("full gc"); flush(stdout); GC.gc(true); println("full gc ok"); flush(stdout)
    end
    println("finished")
end
get(ENV, "IR", "0") == "1" || (V == "a" ? main_a() : main_b())
