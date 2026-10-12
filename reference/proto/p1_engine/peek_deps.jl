# ABOUTME: Checks the plain populations for non-integer exponents in the package's dependent monomials.
# ABOUTME: Also counts Keq exponents other than 0 and -1 (cycles turning over more than once).
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
for p in ["specs", "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi", "handwritten"]
    ms = [m for m in Pop.load_population(p) if !m.allosteric]
    nonint = count(m -> any(!isinteger(e) for (_, mono) in m.dependent_monomials for (_, e) in mono), ms)
    keqs = unique(e for m in ms for (_, mono) in m.dependent_monomials for (s, e) in mono if s == :Keq)
    println(rpad(p, 14), " plain ", length(ms), "  with non-integer dependent exponents ", nonint,
            "  Keq exponents ", keqs)
end
