# ABOUTME: Re-times the slowest core laws (from both bench runs) five times each with the
# ABOUTME: package's protocol, to see how close the core sets sit to the 120 ns gate.
include(joinpath(@__DIR__, "indep.jl"))

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
    allocs, 1e9 * best
end

function main()
    mine = joinpath(SCR, "slice", "skeptic", "rerun", "data")
    rows = [x for dir in (IDATA_IN, mine) for f in readdir(dir)
            if startswith(f, "bench_") && !startswith(f, "bench_terter") &&
               f != "bench_expand_terter.jls"
            for x in deserialize(joinpath(dir, f))]
    slow = sort(rows; by = x -> -x.t1n_ns)
    seen, pick = Set{Tuple{String, Int}}(), []
    for x in slow
        (x.set, x.index) in seen && continue
        push!(seen, (x.set, x.index)); push!(pick, x)
        length(pick) == 6 && break
    end
    rng = Xoshiro(9)
    for x in pick
        r = only(r for r in skeptic_records(x.set) if r.c.index == x.index)
        m = r.c.mech
        em = S.compile(m)
        concs = (; (y => 10.0^(4 * rand(rng) - 2) for y in S.metabolites(m))...)
        params = (; (y => (y === :E_total ? 1.0 : 10.0^(4 * rand(rng) - 2))
                     for y in S.parameters(m))...)
        res = [protocol(em, concs, params) for _ in 1:5]
        @printf("%-14s %-20s D terms %4d ops %5d | bench run %.1f ns | now %s ns, allocs %d\n",
                x.set, x.id, x.nD, x.ops, x.t1n_ns, join((@sprintf("%.1f", r[2]) for r in res), " "),
                maximum(first, res))
    end
end

main()
