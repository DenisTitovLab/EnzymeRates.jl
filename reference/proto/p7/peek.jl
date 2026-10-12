# ABOUTME: Prints records of allosteric specs and dependent exponent statistics for the P7 probe.
# ABOUTME: Throwaway inspection script.
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
specs = Pop.load_population("specs")
for id in ("specs/27", "specs/28", "specs/35", "specs/40")
    m = only(filter(x -> x.id == id, specs))
    println("=== ", m.id, " ", m.source)
    println("fitted: ", m.fitted_params)
    println("metabolites: ", m.metabolites, " regs: ", m.regulators, " sites: ", m.sites, " n=", m.cat_multiplicity)
    println("tags: ", m.group_tags)
    println("deps: ", m.dependent_monomials)
    println("pass2 A/I: ", m.pass2_rename_A, " / ", m.pass2_rename_I)
    for (st, P) in (("A", m.state_polys_A), ("I", m.state_polys_I))
        println(" state $st: num ", length(P[1]), " den ", length(P[2]), " dfree ", P[3])
        println("   den syms: ", sort(unique([s for (mono, c) in P[2] for (s, e) in mono])))
        println("   num: ", P[1])
    end
    println("sample params: ", m.samples[1].params)
    println("rate string:\n", m.rate_equation_string)
end
# dependent exponents across all populations
for pop in Pop.populations()
    ms = Pop.load_population(pop)
    nonint = 0; tot = 0; maxabs = 0
    for m in ms, (k, mono) in m.dependent_monomials, (s, e) in mono
        tot += 1; isinteger(e) || (nonint += 1); maxabs = max(maxabs, abs(numerator(e)))
    end
    coefs = Set{Rational{Int}}()
    for m in ms
        for P in (m.pkg_num, m.pkg_den), (mono, c) in P; push!(coefs, c); end
        for SP in (m.state_polys_A, m.state_polys_I), P in SP, (mono, c) in P; push!(coefs, c); end
    end
    println(pop, ": dep factors $tot nonint $nonint max|e| $maxabs; distinct coefs ", sort(collect(coefs))[1:min(end, 12)], " .. ", length(coefs))
end
