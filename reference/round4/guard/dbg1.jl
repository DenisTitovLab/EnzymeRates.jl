# ABOUTME: Debug: prints a record's steps and its face oracle values under several RE speeds.
include(joinpath(@__DIR__, "lib.jl"))
out = deserialize(joinpath(R, "data", "slice_smoke.jls"))
conv(set) = deserialize(joinpath(SP, "slice/core/data/$set.jls"))
for key in ARGS
    set, idx = split(key, "/")
    c = only(filter(c -> c.index == parse(Int, idx), conv(set)))
    m = c.mech
    println("== $key  subs=$(m.substrates) prods=$(m.products) allo=$(S.is_allosteric(m)) tags=$(m.tags)")
    for (k, e) in enumerate(m.edges)
        println("  [$k g$(m.group[k])] ", S.step_text(m, e))
    end
    r = only(filter(r -> r.key == key, out))
    for f in (:pos, :fwd, :rev); println("  $f: ", r.num[f]); end
end
