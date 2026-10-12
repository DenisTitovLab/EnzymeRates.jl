# ABOUTME: How the candidate normalizations behave on the redesign's fitted-name laws: whether the
# ABOUTME: denominator's constant term and its Π[S] coefficient are single monomials (plain laws).
include(joinpath(@__DIR__, "Laws.jl"))
include(joinpath(@__DIR__, "Factor.jl"))
for pop in ("init_bibi", "init_bibi_pp", "expand_bibi")
    ms = filter(m -> !m.allosteric, Pop.load_population(pop))
    n, c0, cs, zs, nterms0, ntermsS = 0, 0, 0, 0, Int[], Int[]
    for m in ms
        b = Laws.plain_redesign(m)
        law, _, _ = Factor.cancel(b.law)
        den = Canon.normalize(law)[2]
        concs = Set(m.metabolites)
        conc_part(mo) = sort([x for x in mo if x.first in concs]; by = first)
        S = sort([x => 1 for x in m.substrates]; by = first)
        t0 = count(mo -> isempty(conc_part(mo)), keys(den))
        tS = count(mo -> conc_part(mo) == S, keys(den))
        n += 1; c0 += t0 == 1; cs += tS == 1; zs += t0 == 0
        push!(nterms0, t0); push!(ntermsS, tS)
    end
    println(pop, ": ", n, " plain laws; denominator constant term a single monomial in ", c0,
            " (absent in ", zs, ", terms max ", maximum(nterms0), "); Π[S] coefficient a single monomial in ",
            cs, " (terms max ", maximum(ntermsS), ")")
end
