# ABOUTME: Compares the end-to-end recomputation with the census: the fitted sets per round by
# ABOUTME: structural hash, and the D row recomputed on the recomputed mechanisms in run order.
include(joinpath(@__DIR__, "vlib.jl"))

for (name, file) in ((:ter_bi_rxn, "v5_ter_bi_rxn.jls"), (:ter_ter_rxn, "v5_ter_ter_rxn.jls"),
                     (:pyruvate_carboxylase_rxn, "v5_pyruvate_carboxylase_rxn.jls"))
    path = joinpath(VER, "data", file)
    isfile(path) || (println("$name: no recomputation saved"); continue)
    mine = deserialize(path)
    data = load_run(name)
    println("="^100, "\n$name")
    seen_keys = Set{String}()
    for (ri, R) in enumerate(mine)
        theirs = Set(r[:hash] for r in data.rounds[ri] if r[:fitted])
        ours = Set(R.hashes)
        @printf("round %d: recomputed fitted %d, census fitted %d, in both %d\n", ri - 1,
                length(ours), length(theirs), length(intersect(ours, theirs)))
        # the D row on the recomputed mechanisms, in their run order
        nd, nerr = 0, 0
        for m in R.fitted
            k = try
                S.canonical(S.rate_law(to_slice(m))).key
            catch
                nerr += 1
                continue
            end
            k in seen_keys && (nd += 1)
            push!(seen_keys, k)
        end
        @printf("   D recomputed: %d of %d (%.1f%%); conversion or law errors %d; W_admit (own island rule): %d\n",
                nd, length(R.fitted), 100nd / length(R.fitted), nerr, count(>(0), R.lost))
    end
end
