# ABOUTME: Exports the large census's per-candidate classes (census/large/scripts/summarize.jl
# ABOUTME: rules), sizes and keys to one TSV, so the critic can weight them by fit cost.
const LS = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/large/scripts"
include(joinpath(LS, "census_lib.jl"))
include(joinpath(LS, "summarize.jl"))
const OUT = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/critic/large.tsv"
kh(r) = haskey(r, :key) ? string(hash(r[:key]); base = 16) : ""
open(OUT, "w") do io
    println(io, join(["run", "round", "idx", "parent", "fitted", "allo", "np", "np_slice",
                      "dterms", "rank", "pkg_rank", "labels", "explain", "key", "eqh"], '\t'))
    for name in ("ter_bi_rxn", "ter_ter_rxn", "pyruvate_carboxylase_rxn")
        data = deserialize(joinpath(LS, "..", "data", "$(name).jls"))
        label!(data.rounds)
        for (ri, recs) in enumerate(data.rounds), (i, r) in enumerate(recs)
            r[:status] === :ok || continue
            ex = get(r, :explain, nothing)
            exs = ex === nothing ? "" : join(("$k=$v" for (k, v) in ex), ",")
            println(io, join(Any[name, ri - 1, i, get(r, :parent, 0), r[:fitted],
                                 get(r, :allo, false), r[:n_params], get(r, :n_slice, -1),
                                 get(r, :nD, -1), get(r, :rank, -1),
                                 get(r, :pkg_rank_used, -1),
                                 join(string.(get(r, :labels, Symbol[])), ","), exs, kh(r),
                                 string(get(r, :eq_hash, UInt64(0)); base = 16)], '\t'))
        end
        println(name, " done"); flush(stdout)
        data = nothing; GC.gc()
    end
end
println("maxrss MB ", round(Int, Sys.maxrss() / 2^20))
