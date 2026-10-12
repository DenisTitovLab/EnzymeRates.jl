# ABOUTME: Prints every record of saved census runs (one line of facts plus the step list),
# ABOUTME: for checking the census by eye; usage: peek.jl <name,...> [max records]
using Serialization, Printf
const CEN = @__DIR__
maxn = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 60
for name in split(ARGS[1], ",")
    d = deserialize(joinpath(CEN, "data", "$name.jls"))
    println("== ", d.label)
    for r in d.recs[1:min(end, maxn)]
        g(k) = get(r, k, "-")
        println(@sprintf("#%d r%d par=%s %s %s np=%s/%s rank=%s/%s fit=%s key=%s v0=%s dvan=%s face=%s admit=%s ret=%s vmax=%s/%s unb=%s/%s C=%s kcat=%s inv=%s parts=%s",
            r[:gid], r[:round], r[:parent], r[:move], r[:status], g(:np), g(:np_slice),
            g(:rank_pkg), g(:rank_slice), g(:fitted), first(string(g(:key)), 8), g(:v0),
            g(:dvan), g(:face_kind), g(:admit) === nothing ? "ok" : first(string(g(:admit)), 40),
            g(:re_turnover) === nothing ? "-" : "YES", g(:vmax_fwd), g(:vmax_rev), g(:unb_fwd),
            g(:unb_rev), g(:condC), g(:kcat_err) === nothing ? "ok" : first(string(g(:kcat_err)), 30),
            g(:involved), g(:deficit_parts)))
        haskey(r, :slice_err) && println("    slice_err: ", r[:slice_err])
        haskey(r, :pkg_err) && println("    pkg_err: ", r[:pkg_err])
        haskey(r, :rank_err) && println("    rank_err: ", r[:rank_err])
        println("    ", g(:text))
    end
end
