# ABOUTME: Bisects the GC crash of the timing script: each ENV flag switches one ingredient off
# ABOUTME: (ROWS, POP, LAW, TIME) while the rest repeats bigbench.jl's expand_bibi steps.
include(joinpath(@__DIR__, "..", "indep.jl"))
flag(s) = get(ENV, s, "1") == "1"
function protocol(em, concs, params)
    S.rate_equation(em, concs, params)
    allocs = @allocated S.rate_equation(em, concs, params)
    best, acc = Inf, 0.0
    for _ in 1:5
        acc = 0.0
        t = @elapsed for _ in 1:10_000
            acc += S.rate_equation(em, concs, params)
        end
        best = min(best, t / 10_000)
    end
    allocs, best
end
function main()
    if flag("ROWS")
        rows = [x for k in 1:11 for x in deserialize(joinpath(IDATA_IN, "pipe_terter$k.jls")) if x.ok]
        println("rows ", length(rows))
        rows = nothing
        GC.gc()
    end
    recs = flag("POP") ? Dict(r.c.index => r for r in skeptic_records("expand_bibi")) :
        Dict(c.index => (; c) for c in deserialize(joinpath(CORE, "expand_bibi.jls")))
    rng = Xoshiro(3)
    out = NamedTuple[]
    for index in (626, 642, 658)
        m = recs[index].c.mech
        flag("LAW") && S.rate_law(m)
        em = S.compile(m)
        concs = (; (x => 10.0^(4 * rand(rng) - 2) for x in S.metabolites(m))...)
        params = (; (x => (x === :E_total ? 1.0 : 10.0^(4 * rand(rng) - 2))
                     for x in S.parameters(m))...)
        tc = flag("ELAPSED") ? (@elapsed S.rate_equation(em, concs, params)) :
             (S.rate_equation(em, concs, params); 0.0)
        res = flag("TIME") ? [protocol(em, concs, params) for _ in 1:3] : [(0, 0.0)]
        ns = [1e9 * r[2] for r in res]
        flag("OUT") && push!(out, (; index, allocs = maximum(first, res), ns, tc))
        println(index, " done")
    end
    recs = nothing
    GC.gc()
    println("full GC ok")
end
main()
