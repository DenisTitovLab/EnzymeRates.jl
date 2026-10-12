# ABOUTME: Re-times the largest ter-ter laws and expand_bibi/626, 642, 658 in a clean process
# ABOUTME: with the package's protocol (direct calls inside a function), three repeats each.
include(joinpath(@__DIR__, "..", "indep.jl"))

const NTOP = parse(Int, get(ENV, "NTOP", "40"))

"""The package's protocol: warm call, `@allocated`, minimum over 5 batches of 10,000 direct
calls (seconds per call)."""
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
    isfinite(acc) || error("non-finite")
    allocs, best
end

function main()
    rows = [x for k in 1:11
            for x in deserialize(joinpath(IDATA_IN, "pipe_terter$k.jls")) if x.ok]
    top = sort(rows; by = x -> -x.nD)[1:NTOP]
    rng = Xoshiro(3)
    pick = top
    out = NamedTuple[]
    todo = [[("expand_bibi", "expand_bibi", i) for i in (626, 642, 658)];
            [(x.id, x.set, x.index) for x in pick]]
    rows = nothing
    GC.gc()
    for set in unique(t[2] for t in todo)        # one population file in memory at a time
        ps = set == "expand_bibi" ? nothing : deserialize(joinpath(IDATA_IN, "pop_$set.jls"))
        recs = set == "expand_bibi" ? Dict(r.c.index => r for r in skeptic_records(set)) : nothing
        for (label, _, index) in filter(t -> t[2] == set, todo)
            c = ps === nothing ? recs[index].c : S.Converted(ps[index], index)
            m = c.mech
            nD = 0
            em = S.compile(m)
            concs = (; (x => 10.0^(4 * rand(rng) - 2) for x in S.metabolites(m))...)
            params = (; (x => (x === :E_total ? 1.0 : 10.0^(4 * rand(rng) - 2))
                         for x in S.parameters(m))...)
            S.rate_equation(em, concs, params)
            println(label, " ", nD)
            flush(stdout)
        end
        ps = recs = nothing
        GC.gc()
    end
    println("finished")
end

main()
