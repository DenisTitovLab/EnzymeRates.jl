# ABOUTME: Local identifiable rank of the package's Reduced law vs its fitted-parameter count on every plain
# ABOUTME: population mechanism, to test whether "phantoms" (speed-free groups) explain every rank deficiency.
include(joinpath(@__DIR__, "SkLib.jl"))
using Random, Printf, Serialization

const POPS = ["specs", "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi",
              "handwritten"]
rows = []
t0 = time()
for pop in POPS
    ms = [m for m in Pop.load_population(pop) if !m.allosteric && m.ok]
    ndef = 0; defs = Int[]
    for (n, m) in enumerate(ms)
        nθ, rk = try
            Sk.identifiable_rank(m, Random.Xoshiro(31 * n + length(pop)))
        catch e
            println("  ", m.id, " rank ERROR ", first(sprint(showerror, e), 200)); continue
        end
        # second draw guards against an unlucky point
        nθ2, rk2 = Sk.identifiable_rank(m, Random.Xoshiro(97 * n + 5))
        rk = max(rk, rk2)
        push!(rows, (id = m.id, source = m.source, nθ = nθ, rank = rk, fitted = m.fitted_params))
        nθ > rk && (ndef += 1; push!(defs, nθ - rk))
    end
    @printf("%-13s n=%3d  rank < fitted in %3d (deficiency total %d, max %d)  [%.0f s]\n", pop,
            length(ms), ndef, sum(defs; init = 0), maximum(defs; init = 0), time() - t0)
    flush(stdout)
end
serialize(joinpath(@__DIR__, "data", "rank_pop.jls"), rows)
println("\nexamples of deficient mechanisms:")
for r in filter(r -> r.nθ > r.rank, rows)[1:min(end, 25)]
    println("  ", rpad(r.id, 18), " fitted ", r.nθ, " rank ", r.rank, "  ", first(r.source, 70))
end
