# ABOUTME: Package side of the split-difference examples: the fitted parameters and dependents
# ABOUTME: the package derives for each exported parent and child.
using EnzymeRates, Serialization
const S = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(S, "proto/common/Population.jl"))
include(joinpath(@__DIR__, "pkg_lib.jl"))
for p in deserialize(joinpath(S, "slice/enum/data/diag_examples.jls"))
    m = to_er(p, bi_bi_pp_rxn)
    dep, indep = ER._dependent_param_exprs(m)
    println("== ", p.id, " package indep (", length(indep), "): ", join(indep, " "))
    for (k, v) in sort!(collect(dep); by = string ∘ first)
        println("     ", k, " = ", v)
    end
end
