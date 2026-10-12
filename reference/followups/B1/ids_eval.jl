# ABOUTME: Follow-up B1: runs full_eval on a list of mechanism ids (any populations) and prints each verdict
# ABOUTME: under (ii), (iu), (iiu), (Dk) with the direction counts. Usage: julia ids_eval.jl id...
include(joinpath(@__DIR__, "full_eval.jl"))
bypop = Dict{String, Vector{String}}()
for id in ARGS; push!(get!(bypop, split(id, "/")[1], String[]), id); end
for (name, want) in sort!(collect(bypop); by = first)
    pop = Dict(m.id => m for m in load_allo(name))
    for id in want, conv in (:ii, :iu, :iiu, :Dk)
        r = full_eval(pop[id], conv)
        @printf("%-18s %-4s %-8s dim=%d n=%-4d eval=%-4s distinct=%-3s spread=%.1e spread16=%.1e vs_del=%.1e seq=%s zero=%s keep=%d %.2fs\n",
                id, conv, verdict(r), r[:cone_dim], r[:cone_n], get(r, :n_eval, 0), get(r, :n_distinct, "-"),
                get(r, :cone_spread, NaN), get(r, :spread16, NaN), get(r, :cone_vs_del, NaN), get(r, :cone_seq, "-"),
                get(r, :any_zero, "-"), r[:nextra], r[:time])
    end
end
