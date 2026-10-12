# ABOUTME: Split diagnosis, step 2 (package side): the package's independent-parameter count of
# ABOUTME: every exported parent and unit-set child (minus L and the regulator constants).
using EnzymeRates, Serialization, Printf
const S = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(S, "proto/common/Population.jl"))
include(joinpath(@__DIR__, "pkg_lib.jl"))
ps = deserialize(joinpath(S, "slice/enum/data/diag_split_export.jls"))
counts = map(ps) do p
    rxn = startswith(p.id, "expand_bibi") ? bi_bi_pp_rxn : bi_bi_allo_rxn
    try
        m = to_er(p, rxn)
        nreg = m isa ER.AllostericMechanism ? 1 + length(ER._kreg_params(m, :A)) +
               length(ER._kreg_params(m, :I)) : 0
        ER._independent_param_count(m) - nreg
    catch e
        println(p.id, ": ", first(sprint(showerror, e), 200))
        -1
    end
end
serialize(joinpath(S, "slice/enum/data/diag_split_pkgcounts.jls"), counts)
@printf("%d counts, %d errors\n", length(counts), count(==(-1), counts))
