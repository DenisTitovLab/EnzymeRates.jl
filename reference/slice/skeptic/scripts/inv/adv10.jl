# ABOUTME: P2's adv/10: the package's stored steps and names of the mixed-tier group, against the
# ABOUTME: slice's oriented members, to see which tie convention the package uses.
include(joinpath(@__DIR__, "..", "indep.jl"))
r = only(x for x in skeptic_records("p2_adv_pkg") if x.p.id == "adv/10")
p = r.p
for grp in p.groups
    length(grp) > 1 || continue
    for j in grp
        s = p.steps[j]
        println("  stored: ", p.species[s.from].name, " -> ", p.species[s.to].name, "  names ", p.step_names[j])
    end
end
println("  package fitted: ", p.fitted_params)
println("  slice fitted:   ", S.fitted_params(r.c.mech))
