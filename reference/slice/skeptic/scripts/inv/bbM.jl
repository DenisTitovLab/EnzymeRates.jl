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
    big = filter(x -> 1001 <= x.nD <= 5000, rows)
    pick = [top; big[randperm(rng, length(big))[1:NTOP]]]
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
            L = S.rate_law(m)
            nD = L isa S.Law ? S.nterms(L.D) : S.nterms(L.active.D)
            em = S.compile(m)
            concs = (; (x => 10.0^(4 * rand(rng) - 2) for x in S.metabolites(m))...)
            params = (; (x => (x === :E_total ? 1.0 : 10.0^(4 * rand(rng) - 2))
                         for x in S.parameters(m))...)
            S.rate_equation(em, concs, params); tc = 0.0
            res = [protocol(em, concs, params) for _ in 1:3]
            ns = [1e9 * r[2] for r in res]
            allocs = maximum(first, res)
            push!(out, (; label, nD, allocs, ns, tc))
            @printf("%-26s D terms %5d  allocs %d  ns %6.1f %6.1f %6.1f  first call %.2f s\n",
                    label, nD, allocs, ns..., tc)
            flush(stdout)
        end
        ps = recs = nothing
        GC.gc()
    end
    serialize(joinpath(SCR, "slice", "skeptic", "data", "bigbench.jls"), out)
    over(f) = count(x -> f(x.ns) >= 120, out)
    @printf("laws %d; allocations > 0: %d; over 120 ns by the best of 3: %d, by the worst of 3: %d\n",
            length(out), count(x -> x.allocs > 0, out), over(minimum), over(maximum))
end

main()
