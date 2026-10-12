# ABOUTME: Query: labels and admit verdicts of the records in the both-ways same-law sample.
include(joinpath(@__DIR__, "vlib.jl"))
for (name, ids) in ((:ter_bi_rxn, [(1,48),(1,305),(1,450),(1,1255),(2,138),(2,563),(2,208),(2,396),(2,398)]),
                    (:ter_ter_rxn, [(1,281),(1,321),(1,419),(1,725),(1,889),(1,1212),(1,2099),(1,2149),(1,2159)]),
                    (:pyruvate_carboxylase_rxn, [(1,363),(1,509),(1,763),(1,932),(1,1519),(1,1986),(1,1955)]))
    data = load_run(name)
    byid = Dict(r[:id] => r for recs in data.rounds for r in recs)
    for id in ids
        r = byid[id]
        println(name, " ", id, " labels ", r[:labels], " explain ", get(r, :explain, ""), " admit ", r[:admit] === nothing ? "ok" : first(r[:admit], 80))
    end
end
