# ABOUTME: Splits the degenerate seeds of each census sample by the seed screen's conditions
# ABOUTME: (V: no Vmax forward or reverse; C: chemistry in equilibrium with both sides).
include(joinpath(@__DIR__, "census_lib.jl"))
for name in Symbol.(ARGS)
    rxn = getfield(Main, name)
    seeds = unique!(collect(ER.init_mechanisms(rxn)))
    n = length(seeds)
    idx = n > 2000 ? sort!(randperm(Xoshiro(20261011), n)[1:2000]) : collect(1:n)
    tally = Dict{Tuple{Bool, Bool, Bool}, Int}()
    ex = Dict{Tuple{Bool, Bool, Bool}, String}()
    for m in seeds[idx]
        ER._degenerate(m) || continue
        g = ER.steps(m)
        k = (!ER._has_vmax(g, rxn, ER.Substrate), !ER._has_vmax(g, rxn, ER.Product),
             ER._chemistry_equilibrates_both_sides(g, rxn))
        tally[k] = get(tally, k, 0) + 1
        haskey(ex, k) || (ex[k] = mech_text(to_slice(m)))
    end
    println("$name: degenerate seeds in the 2000-seed sample by (no fwd Vmax, no rev Vmax, condition C): ", tally)
    for (k, t) in ex
        println("  e.g. $k: ", first(t, 700))
    end
    flush(stdout)
end
