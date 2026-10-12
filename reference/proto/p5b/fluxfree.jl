# ABOUTME: Checks whether rank deficits not caused by L come from flux-free steady-state groups
# ABOUTME: (P2's Thermo.zero_flux_ss) in either conformation's graph; probe P5b.
using Serialization
include(joinpath(@__DIR__, "Allo.jl"))
for name in ("allo_bibi", "specs", "handwritten", "expand_bibi", "allo_uniuni")
    rs = deserialize(joinpath(@__DIR__, "data", "results_$(name).jls"))
    P = Dict(m.id => m for m in Pop.load_population(name) if m.allosteric)
    tab = Dict{Any, Int}()
    for r in rs
        m = P[r[:id]]; md = Allo.derive(m)
        ff = 0
        for keep in (trues(length(m.steps)), md.keepI)
            any(keep) || continue
            pm, js = Allo.submech(m, keep)
            z = Thermo.zero_flux_ss(pm)
            g = Pop.group_of_step(m)
            # whole kinetic groups (within this state) that are flux-free SS
            for grp in unique(g[js])
                mem = [k for k in eachindex(js) if g[js[k]] == grp]
                all(k -> !pm.steps[k].re && z[k], mem) && (ff += 1)
            end
        end
        k = (deficit = r[:count] - r[:rank_mine][1], fluxfree_groups = ff)
        tab[k] = get(tab, k, 0) + 1
    end
    println(name, ": ", sort!(collect(tab); by = string ∘ first))
end
