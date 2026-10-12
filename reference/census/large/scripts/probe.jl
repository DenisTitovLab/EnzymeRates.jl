# ABOUTME: Smoke test of the census on bi-bi (small rounds, package-rank check) and timings of
# ABOUTME: init_mechanisms and of the per-candidate stages on ter-bi and the pyruvate reactions.

include(joinpath(@__DIR__, "driver.jl"))
include(joinpath(@__DIR__, "summarize.jl"))

info, rounds = run_census(:bi_bi_rxn; N0 = 2000, N1 = 12, N2 = 6, pkgrank_n = 40)
print_tables(info, rounds; nexamples = 2)
print_costs(rounds)
flush(stdout)

rng = Xoshiro(5)
for name in (:ter_bi_rxn, :pyruvate_carboxylase_rxn, :pyruvate_dehydrogenase_rxn)
    rxn = getfield(Main, name)
    t = @elapsed seeds = unique!(collect(ER.init_mechanisms(rxn)))
    @printf("%s: %d seeds in %.1f s; maxrss %.0f MB\n", name, length(seeds), t,
            Sys.maxrss() / 2^20)
    pick = seeds[sort!(randperm(rng, length(seeds))[1:min(40, length(seeds))])]
    tp, ts = Float64[], Float64[]
    st = Dict{Symbol, Int}()
    for m in pick
        t0 = time_ns()
        r = pkg_stage(m)
        push!(tp, secs(t0))
        st[r.status] = get(st, r.status, 0) + 1
        r.status === :ok || continue
        t0 = time_ns()
        slice_stage(m, rng)
        push!(ts, secs(t0))
    end
    @printf("  40 seeds: statuses %s; package stage median %.1f ms max %.0f ms; slice \
            stage median %.1f ms max %.0f ms\n", st, 1e3 * median(tp), 1e3 * maximum(tp),
            isempty(ts) ? NaN : 1e3 * median(ts), isempty(ts) ? NaN : 1e3 * maximum(ts))
    flush(stdout)
end
