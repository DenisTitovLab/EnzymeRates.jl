# ABOUTME: Package side: flags ROADMAP known issues 1 (an untagged inactive-state name the active
# ABOUTME: state lacks) and 2 (a name folded in one state but used by the other) per mechanism.
using EnzymeRates, Serialization, Printf
const S = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(S, "proto/common/Population.jl"))
include(joinpath(@__DIR__, "pkg_lib.jl"))
ps = deserialize(joinpath(S, "slice/enum/data/diag_split_export.jls"))
flags = map(ps) do p
    am = to_er(p, startswith(p.id, "expand_bibi") ? bi_bi_pp_rxn : bi_bi_allo_rxn)
    names, renames = Dict{Symbol, Vector{Symbol}}(), Dict{Symbol, Vector{Pair{Symbol, Symbol}}}()
    for st in (:A, :I)
        cm, sp = ER._state_parts(am, st)
        names[st] = unique(Symbol[ER.name(x, cm) for ps in sp for x in ps])
        renames[st] = collect(ER._build_wegscheider_rename_map(cm; step_params = sp))
    end
    tagged(n) = startswith(string(n), "K_I_") || startswith(string(n), "k_I_")
    ki1 = any(n -> !tagged(n) && !(n in names[:A]), names[:I])
    ki2 = any(any(r -> r.first in names[st == :A ? :I : :A], renames[st]) for st in (:A, :I))
    (; ki1, ki2)
end
serialize(joinpath(S, "slice/enum/data/diag_split_issues.jls"), flags)
@printf("%d mechanisms: known issue 1 in %d, known issue 2 in %d\n", length(flags),
        count(f -> f.ki1, flags), count(f -> f.ki2, flags))
