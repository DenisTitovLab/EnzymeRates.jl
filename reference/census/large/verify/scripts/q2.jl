# ABOUTME: Queries on the saved records: v = 0 laws among fitted, the steady-state-local P cases
# ABOUTME: (their NonequalAI steady-state groups and whether the inactive state turns over).
include(joinpath(@__DIR__, "vlib.jl"))
for name in (:ter_bi_rxn, :ter_ter_rxn, :pyruvate_carboxylase_rxn)
    data = load_run(name)
    fit = [r for recs in data.rounds for r in recs if r[:fitted]]
    println(name, ": fitted with v0 = true: ", count(r -> get(r, :v0, false), fit),
            "; with a dead face: ", count(r -> !isempty(get(r, :dead, String[])), fit),
            "; with N zero on a face: ", count(r -> !isempty(get(r, :nzero, String[])), fit))
    ss = [r for r in fit if haskey(get(r, :explain, Dict()), :ss_local)]
    for r in ss
        m = r[:smech]
        neq = [g for g in 1:S.ngroups(m) if m.tags[g] === :NonequalAI]
        neq_ss = [g for g in neq if any(k -> !m.edges[k].re, S.members(m, g))]
        imask = S.inactive_mask(m)
        # chemistry edges: not a plain binding in their physical-forward direction
        fwd = S.orient(m)
        chem_in = [k for k in findall(imask) if !S.is_plain_binding(m, S.oriented(m, k, fwd))]
        println("   ", r[:id], " NonequalAI groups ", neq, " (steady-state ", neq_ss, "), chemistry edges in the inactive conformation ", length(chem_in),
                ", rank ", r[:rank], "/", r[:n_slice])
    end
end
