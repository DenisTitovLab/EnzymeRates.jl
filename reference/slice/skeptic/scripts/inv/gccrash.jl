# ABOUTME: Minimal attempt to reproduce the GC corruption seen in the timing runs: compile one
# ABOUTME: large law, time it with direct calls and @allocated, then force full collections.
const SCR = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(SCR, "proto", "common", "Population.jl"))
include(joinpath(SCR, "slice", "src", "Slice.jl"))
using Serialization, Random, Printf
const S = Slice
const MODE = get(ENV, "MODE", "direct")

function protocol(em, concs, params)
    S.rate_equation(em, concs, params)
    allocs = MODE == "noalloc" ? 0 : @allocated S.rate_equation(em, concs, params)
    best, acc = Inf, 0.0
    for _ in 1:5
        acc = 0.0
        t = @elapsed for _ in 1:10_000
            acc += S.rate_equation(em, concs, params)
        end
        best = min(best, t / 10_000)
    end
    allocs, best, acc
end

function main()
    cs = deserialize(joinpath(SCR, "slice", "core", "data", "expand_bibi.jls"))
    rng = Xoshiro(3)
    for idx in (626, 642, 658)
        m = only(c for c in cs if c.index == idx).mech
        MODE == "lawfirst" && S.rate_law(m)
        em = S.compile(m)
        concs = (; (x => 10.0^(4 * rand(rng) - 2) for x in S.metabolites(m))...)
        params = (; (x => (x === :E_total ? 1.0 : 10.0^(4 * rand(rng) - 2))
                     for x in S.parameters(m))...)
        for _ in 1:3
            a, b, acc = protocol(em, concs, params)
            @printf("%d: allocs %d, %.1f ns\n", idx, a, 1e9 * b)
        end
        L = S.rate_law(m)
        GC.gc()
        println("  full GC ok after ", idx)
    end
    for _ in 1:3
        GC.gc()
    end
    println("done")
end

main()
