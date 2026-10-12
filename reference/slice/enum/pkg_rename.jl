# ABOUTME: Package side: the per-state Wegscheider rename maps and inactive-state names of the
# ABOUTME: split-difference examples, to tie each difference to a ROADMAP known issue.
using EnzymeRates, Serialization
const S = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(S, "proto/common/Population.jl"))
include(joinpath(@__DIR__, "pkg_lib.jl"))
for p in deserialize(joinpath(S, "slice/enum/data/diag_examples.jls"))
    am = to_er(p, bi_bi_pp_rxn)
    println("== ", p.id)
    for st in (:A, :I)
        cm, sp = ER._state_parts(am, st)
        rn = ER._build_wegscheider_rename_map(cm; step_params = sp)
        names = unique(Symbol[ER.name(x, cm) for ps in sp for x in ps])
        println("   $st rename ", collect(rn), "\n   $st names ", names)
    end
end
