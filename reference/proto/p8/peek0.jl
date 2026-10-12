# ABOUTME: First look at the populations for the P8 key probe: eq_hash duplicates and record fields.
# ABOUTME: Stdlib only; prints counts and one record per kind.
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
using SHA
println(VERSION)
for name in ("init_bibi", "init_bibi_pp", "expand_bibi", "allo_bibi")
    ms = Pop.load_population(name)
    plain = filter(m -> !m.allosteric, ms); allo = filter(m -> m.allosteric, ms)
    for (lab, sub) in (("plain", plain), ("allo", allo))
        isempty(sub) && continue
        nk = length(unique(m.eq_hash for m in sub))
        nstr = length(unique(m.rate_equation_string for m in sub))
        nfp = [length(m.fitted_params) for m in sub]
        println(rpad(name, 14), rpad(lab, 6), " n=", length(sub), " eq_hash classes=", nk,
                " string classes=", nstr, " fitted range=", extrema(nfp),
                " ok=", count(m -> m.ok, sub))
    end
end
ms = Pop.load_population("allo_bibi")
m = ms[1]
for f in (:id, :source, :family, :substrates, :products, :metabolites, :fitted_params, :dependents,
          :dependent_monomials, :pass2_rename_A, :pass2_rename_I, :sites, :cat_multiplicity,
          :group_tags)
    println(f, " = ", getfield(m, f))
end
println("state_polys_A sizes ", length.(m.state_polys_A), " I ", length.(m.state_polys_I))
println("A dfree: ", m.state_polys_A[3]); println("I dfree: ", m.state_polys_I[3])
println("sample params: ", m.samples[1].params); println("v=", m.samples[1].v)
println(m.rate_equation_string[1:min(end, 1500)])
