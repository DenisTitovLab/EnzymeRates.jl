# ABOUTME: Prints a few fields of chosen mechanisms of a population for a quick look.
# ABOUTME: Usage: julia peek.jl <population> <index...>
include(joinpath(@__DIR__, "..", "Population.jl"))
pop = Pop.load_population(ARGS[1])
for i in parse.(Int, ARGS[2:end])
    m = pop[i]
    println("== ", m.id, "  ", m.source, "  ok=", m.ok, "  allo=", m.allosteric, "  family=", m.family)
    println("subs=", m.substrates, " prods=", m.products, " inh=", m.inhibitors, " regs=", m.regulators, " mets=", m.metabolites)
    for (j, s) in enumerate(m.steps)
        println("  step $j: ", Pop.step_text(m, s), "   names=",
                m.allosteric ? (m.step_names_A[j], m.step_names_I[j], isempty(m.step_names_I_pkg) ? "-" : m.step_names_I_pkg[j], isempty(m.inactive_kept) ? "-" : m.inactive_kept[j]) : m.step_names[j])
    end
    println("groups=", m.groups, " tags=", m.group_tags, " sites=", m.sites)
    println("fitted=", m.fitted_params, "\nfull=", m.full_params)
    println("deps=", m.dependents)
    println("pass2=", m.pass2_rename, " changes=", m.pass2_changes_polys, " ipc=", m.independent_param_count)
    println("raw_num terms=", length(m.raw_num), " raw_den terms=", length(m.raw_den), " dfree=", m.raw_dfree)
    println("errors=", m.stage_errors, " ctor_ns=", m.ctor_ns)
    for s in m.samples
        println("  v=", s.v, " v_full=", s.v_full, " kcat=", s.kcat, " ndep=", length(s.dep_values))
    end
    println(m.rate_equation_string)
end
