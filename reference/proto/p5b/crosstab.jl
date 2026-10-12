# ABOUTME: Cross-tabulates P5b count differences against the package's phantom names and fold hazards.
# ABOUTME: Stdlib only.
using Serialization
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
for name in ("allo_bibi", "expand_bibi")
    rs = deserialize(joinpath(@__DIR__, "data", "results_$(name).jls"))
    tab = Dict{Any, Int}()
    for r in rs
        k = (diff = r[:pkg_count] - r[:count], undeclared = length(r[:undeclared_pkg]),
             opp = r[:pkg_opp], one = r[:pkg_one], myfold = r[:folds] > 0)
        tab[k] = get(tab, k, 0) + 1
    end
    println("== $name (pkg-mine, #undeclared pkg names, opposite fold, one-state fold, my fold) => n")
    for (k, n) in sort!(collect(tab); by = x -> Tuple(x.first))
        println("   ", k, " => ", n)
    end
end
# one example of each hazard
P = Dict(m.id => m for m in Pop.load_population("allo_bibi"))
rs = deserialize(joinpath(@__DIR__, "data", "results_allo_bibi.jls"))
for (lab, f) in (("opposite", r -> r[:pkg_opp]), ("one-state", r -> r[:pkg_one] && !r[:pkg_opp]))
    r = rs[findfirst(f, rs)]; m = P[r[:id]]
    println("== $lab example $(r[:id]): p2A=", m.pass2_rename_A, " p2I=", m.pass2_rename_I)
    println("   pkg fitted ($(r[:pkg_count])): ", m.fitted_params)
    println("   my fitted  ($(r[:count])): ", r[:fitted])
    g = Pop.group_of_step(m)
    for (j, s) in enumerate(m.steps)
        println("   step $j g$(g[j]) [$(m.group_tags[g[j]])] ", Pop.step_text(m, s), " A=", m.step_names_A[j], " Ipkg=", m.step_names_I_pkg[j])
    end
    println(r[:text][1:findfirst("v = ", r[:text])[1]-1])
end
