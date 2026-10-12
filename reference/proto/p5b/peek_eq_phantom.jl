# ABOUTME: Shows allo_bibi mechanisms whose package fit has undeclared names at an equal count (P5b).
# ABOUTME: Stdlib only.
using Serialization
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
P = Dict(m.id => m for m in Pop.load_population("allo_bibi"))
rs = deserialize(joinpath(@__DIR__, "data", "results_allo_bibi.jls"))
for r in filter(r -> r[:count] == r[:pkg_count] && !isempty(r[:undeclared_pkg]), rs)[1:2]
    m = P[r[:id]]
    println("== $(r[:id]) undeclared ", r[:undeclared_pkg])
    println("   pkg fitted: ", m.fitted_params, "\n   pkg deps: ", m.dependents)
    println("   my fitted:  ", r[:fitted])
    g = Pop.group_of_step(m)
    for (j, s) in enumerate(m.steps)
        println("   step $j g$(g[j]) [$(m.group_tags[g[j]])] ", Pop.step_text(m, s), " A=", m.step_names_A[j], " Ipkg=", m.step_names_I_pkg[j])
    end
    println(r[:text][1:findfirst("v = ", r[:text])[1]-1])
end
