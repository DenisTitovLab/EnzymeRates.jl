# ABOUTME: Prints allo_terbi/131 (the one island rejection whose limit equals deletion) and two
# ABOUTME: other rejected ping-pong mechanisms with B1's step table and the admission message.
include(joinpath(@__DIR__, "lib.jl"))
module B1
import Main.Pop
include(joinpath(Main.SCRATCH, "followups", "B1", "b1_lib.jl"))
end
recs = load_records(("allo_terbi", "expand_bibi", "allo_bibi_pp"))
for id in ("allo_terbi/131", "allo_terbi/130", "expand_bibi/308", "allo_bibi_pp/314")
    r = only(filter(x -> x.p.id == id, recs))
    B1.print_mech(stdout, r.p, :iu)
    println("   admit: ", S.admit(r.c.mech))
    full = deserialize(joinpath(SCRATCH, "followups", "B1", "data",
                                id == "allo_terbi/131" || id == "allo_terbi/130" ? "full_allo_terbi.jls" :
                                id == "expand_bibi/308" ? "full_expand_bibi_251.jls" : "full_allo_bibi_pp_201.jls"))
    for x in filter(x -> x[:id] == id && x[:conv] === :iu, full)
        println("   B1 full (iu): ", B1.verdict(x), " dirs ", get(x, :n_eval, 0), " spread ",
                get(x, :cone_spread, NaN), " vs deletion ", get(x, :cone_vs_del, NaN))
    end
end
