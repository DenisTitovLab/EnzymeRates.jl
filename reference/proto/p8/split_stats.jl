# ABOUTME: Splits the saved P8 class counts of a mixed population into plain and MWC mechanisms.
# ABOUTME: Reads data/<pop>.jls and prints eq_hash, exact-name and rename-proof class counts per kind.
using Serialization
for pop in ARGS
    d = deserialize(joinpath(@__DIR__, "data", "$pop.jls"))
    for allo in (false, true)
        rs = [r for r in d.recs if r[:allosteric] == allo]
        isempty(rs) && continue
        println(pop, allo ? " MWC" : " plain", ": n = ", length(rs), ", eq_hash ",
                length(unique(r[:eq_hash] for r in rs)))
        for variant in (:today, :redesign), kk in (:k1, :k2)
            ks = [r[(variant, kk, :key)] for r in rs if haskey(r, (variant, kk, :key))]
            ex = [r[(variant, kk, :exact)] for r in rs if haskey(r, (variant, kk, :exact))]
            println("   ", variant, " ", kk, ": exact-name ", length(unique(ex)), ", rename-proof ",
                    length(unique(ks)))
        end
    end
end
