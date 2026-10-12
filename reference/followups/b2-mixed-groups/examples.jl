# ABOUTME: Prints example mechanisms with mixed groups (RE and SS, dead-end and catalytic plain
# ABOUTME: member, :OnlyA) with their names, fitted parameters and Haldane/Wegscheider dependents.
include(joinpath(@__DIR__, "..", "..", "proto", "common", "Population.jl"))
include(joinpath(@__DIR__, "mixed_lib.jl"))

function show_full(m)
    print_mech(m)
    names = m.allosteric ? m.step_names_A : m.step_names
    for g in mixed_groups(m)
        println("     mixed group constants: ", names[first(m.groups[g])],
                m.allosteric ? "  (inactive: $(m.step_names_I[first(m.groups[g])]))" : "")
    end
    println("     fitted: ", m.fitted_params)
    for (k, v) in m.dependents
        println("     dependent ", k, " = ", v)
    end
    isempty(m.pass2_rename) || println("     pass-2 rename: ", m.pass2_rename)
end

pops = Dict(n => Pop.load_population(n) for n in ("init_bibi", "init_bibi_pp", "expand_bibi",
                                                   "allo_bibi"))
byid(id) = only(filter(x -> x.id == id, pops[split(id, "/")[1]]))
foreach(id -> show_full(byid(id)), ("init_bibi/56", "init_bibi/57", "allo_bibi/826"))

# First mechanism of each further shape: RE random-order mixed group (catalytic plain member),
# an :OnlyA RE mixed group in expand_bibi, and a mechanism with an RE isomerization.
first_where(pop, f) = (i = findfirst(f, pops[pop]); i === nothing ? nothing : pops[pop][i])
re_random = first_where("init_bibi", m -> any(g -> m.steps[first(m.groups[g])].re &&
    any(j -> step_class(m, m.steps[j]) === :plain && !deadend(m.species[m.steps[j].to]),
        m.groups[g]), mixed_groups(m)))
onlya_re = first_where("expand_bibi", m -> m.allosteric && any(g -> m.group_tags[g] === :OnlyA &&
    m.steps[first(m.groups[g])].re, mixed_groups(m)))
re_iso = first_where("init_bibi_pp", m -> any(s -> s.re && step_class(m, s) === :iso, m.steps))
for m in (re_random, onlya_re, re_iso)
    m === nothing ? println("-- none") : show_full(m)
end
