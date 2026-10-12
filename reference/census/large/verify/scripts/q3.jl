# ABOUTME: Query: step lists of the records the verification report cites.
include(joinpath(@__DIR__, "vlib.jl"))
for (name, ids) in ((:ter_bi_rxn, [(1, 639), (0, 526)]), (:pyruvate_carboxylase_rxn, [(1, 1955), (0, 2251)]))
    data = load_run(name)
    byid = Dict(r[:id] => r for recs in data.rounds for r in recs)
    for id in ids
        r = byid[id]
        println(name, " ", id, " parent ", r[:parent], " labels ", r[:labels], " slice ", r[:rank], "/", r[:n_slice], " pkg ", r[:pkg_rank_used], "/", r[:n_params])
        println("   ", r[:text])
    end
end
