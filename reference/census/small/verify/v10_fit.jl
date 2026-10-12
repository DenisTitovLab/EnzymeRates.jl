# ABOUTME: Fits the package law of :OnlyA-rejected and admitted candidates to mass-action steady states:
# ABOUTME: the deletion model (the package's), and for islands the :OnlyA family with the island populated.
include(joinpath(@__DIR__, "vlib.jl"))
include(joinpath(@__DIR__, "moracle.jl"))
const SAMPLES = deserialize(joinpath(VER, "samples.jls"))

rng = Xoshiro(101010)
tally = Dict{String, Vector{Int}}()
bump!(k, hit) = (v = get!(tally, k, [0, 0]); v[1] += 1; v[2] += hit)
census_kind(msg) = msg === nothing ? "none" : occursin("still join", msg) ? "G1_inside" :
                   occursin("close a cycle", msg) ? "G1_cycle" : occursin("reach the free", msg) ? "island" : "other"

function check(lab, name, d, mechs, g)
    r = d.recs[g]
    ms = mstruct(mechs[g])
    pl = plaw(mechs[g])
    A = aconsts(ms, rng)
    Ip = iconsts(ms, A, rng)
    xp = Ip.draw(rng)
    lnkp = Ip.lnk(xp)
    L = big(lu10(rng, -1, 1))
    σ = [big(lu10(rng, -1, 1)) for _ in 1:(length(ms.steps) + 1)]
    keq = Float64(exp(BigFloat(sum(A.mu[x] for x in ms.subs) - sum(A.mu[x] for x in ms.prods))))
    npts = 2 * length(pl.names) + 14
    pts = [Dict{Symbol, Float64}(x => lu10(rng, -1, 1) for x in ms.mets) for _ in 1:npts]
    bpts = [Dict{Symbol, BigFloat}(x => big(v) for (x, v) in c) for c in pts]
    vdel = [Float64(oracle_rate(ms, A, Ip, lnkp, L, c, big(10)^30, σ)) for c in bpts]
    rdel = fit_values(pl, keq, pts, vdel, rng)
    line = @sprintf("  [%s] %s #%d census %s np %d | package law fitted to the deletion-model steady state: %.1e",
                    lab, name, g, census_kind(r[:admit]), length(pl.names), rdel)
    bump!("$lab: package law = deletion-model ODE (fit < 1e-8)", rdel < 1e-8)
    if lab == "island"
        If = iconsts(ms, A, rng; family = true, pin = Ip.pins(xp))
        if If !== nothing
            rs = Float64[]
            for _ in 1:2
                l = If.lnk(If.draw(rng))
                vf = [Float64(oracle_rate(ms, A, If, l, L, c, big(10)^30, σ)) for c in bpts]
                spread = maximum(abs.(vf .- vdel) ./ abs.(vdel))
                push!(rs, fit_values(pl, keq, pts, vf, rng))
                line *= @sprintf(" | island-populated draw: differs from deletion by %.1e, package-law fit %.1e", spread, rs[end])
            end
            bump!("island: populated variant inside the package family (fit < 1e-8)", all(<(1e-8), rs))
        end
    end
    println(line)
    flush(stdout)
end

for name in ["uni_bi_rxn", "bi_bi_rxn", "bi_bi_pp_rxn", "bi_bi_rxn.r2_100000", "bi_bi_pp_rxn.r2_1200"]
    d, mechs = loadrun(name)
    s = SAMPLES[name]
    println("\n== $name")
    for g in s[:W_island][1:min(end, 10)]
        check("island", name, d, mechs, g)
    end
    for g in s[:W_G1][1:min(end, 10)]
        check("G1", name, d, mechs, g)
    end
    for g in s[:W_unflagged_allo][1:min(end, 8)]
        check("admitted", name, d, mechs, g)
    end
    d = mechs = nothing
    GC.gc()
end
println("\nTally (n, hits):")
for (k, v) in sort!(collect(tally))
    println("  $k: $(v[2]) of $(v[1])")
end
