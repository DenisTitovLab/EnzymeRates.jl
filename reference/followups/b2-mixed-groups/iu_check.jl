# ABOUTME: Which P5a convention-(iu) failures in expand_bibi carry an :OnlyA mixed group, and
# ABOUTME: whether every (i) failure with a mixed group becomes unique = deletion under (iu).
include(joinpath(@__DIR__, "..", "..", "proto", "common", "Population.jl"))
include(joinpath(@__DIR__, "mixed_lib.jl"))
using Serialization
const P5A = joinpath(@__DIR__, "..", "..", "proto", "p5a", "data")
for name in ("expand_bibi", "allo_bibi")
    byid = Dict(m.id => m for m in Pop.load_population(name))
    recs = deserialize(joinpath(P5A, "pop_$name.jls"))
    ok(r) = !haskey(r, :error) && haskey(r, :cone_v) && r[:cone_spread] < 1e-20 && r[:cone_vs_del] < 1e-20
    fails(conv) = Set(r[:id] for r in recs if r[:conv] === conv && !ok(r))
    has_mixed(id) = any(g -> byid[id].group_tags[g] === :OnlyA, mixed_groups(byid[id]))
    fi, fiu = fails(:i), fails(:iu)
    mi = filter(has_mixed, fi)
    println(name, ": (i) failures ", length(fi), ", with mixed ", length(mi),
            "; (iu) failures ", length(fiu), ", with mixed ", count(has_mixed, fiu),
            "; mixed (i) failures still failing under (iu): ", sort!(collect(intersect(mi, fiu))),
            "; ping-pong among them: ", count(id -> :pingpong in byid[id].family, intersect(mi, fiu)))
end
