# ABOUTME: Whether deleting the stranded inactive forms is the ε → 0 limit when every :OnlyA
# ABOUTME: step between components loses only its rate away from free E (mode :away).
include(joinpath(@__DIR__, "lib.jl"))
include(joinpath(@__DIR__, "fixtures_def.jl"))

const EPS = (big(10)^-4, big(10)^-6, big(10)^-8)
fmt(x) = @sprintf("%.1e", x)
println("package: ", pathof(EnzymeRates))
for (label, subs, prods, inh, steps) in fixtures
    m = build(dsl(subs, prods, steps; n = 1, inhibitors = inh))
    am = ER.AllostericMechanism(m)
    allforms = Set(ER.name(f) for g in ER.steps(am) for s in g
                   for f in (ER.from_species(s), ER.to_species(s)))
    isempty(setdiff(allforms, inactive_forms(am))) && continue
    ia = island_analysis(am)
    println("\n=== $label (n = 1); distances from E's component: ",
            join(["$(join(sort(ia.comps[c]), ",")) $d" for (c, d) in ia.dist], " | "))
    rng = MersenneTwister(hash(label) % UInt32)
    worst = Dict{String, BigFloat}()
    for draw in 1:4
        vals = all_values(am, draw_fitted(rng, ER.fitted_params(m)),
                          exp(BigFloat(log(0.5) + log(40) * rand(rng))))
        gen() = Dict{Symbol, BigFloat}(Symbol(x) => draw_conc(rng) for x in [subs; prods; inh])
        pts = ["generic" => gen(), "products 0" => merge(gen(), Dict(Symbol(x) => big(0) for x in prods)),
               "P = 1e-6" => merge(gen(), Dict(Symbol(prods[end]) => big(10)^-6))]
        for (pname, concs) in pts
            v = pkg_rate(m, vals, concs)
            errs = String[]
            for ε in EPS
                sAe, sIe, miss = oracle_steps(am, vals; mode = :away, ε)
                isempty(miss) || error("missing inactive constants $miss")
                ve, _, _ = oracle_v(sAe, sIe, concs, vals[:L], Symbol(subs[1]); strict = true)
                e = relerr(ve, v)
                push!(errs, fmt(e))
                worst[pname * " ε=" * fmt(ε)] = max(get(worst, pname * " ε=" * fmt(ε), big(0)), e)
            end
            println("  draw $draw, $(rpad(pname, 10)): v = $(@sprintf("%+.4e", v)) | away ε=1e-4,-6,-8: ",
                    join(errs, " "))
        end
        if draw == 1
            ce = eq_concs(rng, am, vals[:Keq], Dict{Symbol, BigFloat}(Symbol(x) => draw_conc(rng) for x in inh))
            sAe, sIe, _ = oracle_steps(am, vals; mode = :away, ε = EPS[1])
            ve, _, _ = oracle_v(sAe, sIe, ce, vals[:L], Symbol(subs[1]); strict = true)
            println("  away ε=1e-4 network at equilibrium: |v| = ", fmt(abs(ve)))
        end
    end
    println("  worst: ", join(["$k $(fmt(v))" for (k, v) in sort(collect(worst))], "; "))
end
