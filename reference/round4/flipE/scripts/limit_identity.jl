# ABOUTME: For R3 and R4 the kinetic-block ε-limit equals the law of the same mechanism with the
# ABOUTME: :OnlyA isomerization tagged :NonequalAI and its inactive ratio set to the I-copy's.
include(joinpath(@__DIR__, "lib.jl"))
include(joinpath(@__DIR__, "fixtures_def.jl"))

println("package: ", pathof(EnzymeRates))
for (label, subs, prods, inh, steps) in fixtures
    (startswith(label, "R3") || startswith(label, "R4")) || continue
    m = build(dsl(subs, prods, steps; inhibitors = inh))
    twin_steps = [(s, t == "OnlyA" ? "NonequalAI" : t) for (s, t) in steps]
    twin = build(dsl(subs, prods, twin_steps; inhibitors = inh))
    am, amt = ER.AllostericMechanism(m), ER.AllostericMechanism(twin)
    println("\n=== $label; twin fitted_params: ", ER.fitted_params(twin))
    rng = MersenneTwister(3)
    for draw in 1:3, base in (big(1), big(7))
        vals = all_values(am, draw_fitted(rng, ER.fitted_params(m)),
                          exp(BigFloat(log(0.5) + log(40) * rand(rng))))
        # The twin's inactive constants of the isomerization: the active ones, the
        # equilibrium ratio times `base` (the I-copy the ε network uses).
        tv = copy(vals)
        for p in ER.fitted_params(twin)
            haskey(tv, p) && continue
            s = String(p)
            a = Symbol(replace(s, "_I_" => "_A_"))
            tv[p] = vals[a] * (startswith(s, "K_") || occursin("E_to_F", s) ? base : big(1))
        end
        tvals = all_values(amt, Dict(p => tv[p] for p in ER.fitted_params(twin)), vals[:Keq])
        concs = Dict{Symbol, BigFloat}(Symbol(x) => draw_conc(rng) for x in [subs; prods; inh])
        sAe, sIe, _ = oracle_steps(am, vals; mode = :literal, ε = big(10)^-12, base)
        ve, _, _ = oracle_v(sAe, sIe, concs, vals[:L], Symbol(subs[1]); strict = true)
        vt = pkg_rate(twin, tvals, concs)
        vd = pkg_rate(m, vals, concs)
        @printf("  draw %d base %d: ε=1e-12 network vs :NonequalAI twin %.1e | vs deletion law %.1e\n",
                draw, Int(base), relerr(ve, vt), relerr(ve, vd))
    end
end
