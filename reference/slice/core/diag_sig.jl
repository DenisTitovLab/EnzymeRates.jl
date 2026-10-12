# ABOUTME: Lists the method instances that encoding a Sig of a new length compiles, to find
# ABOUTME: what in `sig` specializes on the tuple length.
include(joinpath(@__DIR__, "..", "src", "Slice.jl"))
using Serialization
function snapshot()
    out = Set{Core.MethodInstance}()
    Base.visit(m -> foreach(mi -> push!(out, mi), Base.specializations(m)), Core.methodtable)
    out
end
ms = [c.mech for c in deserialize(joinpath(@__DIR__, "data", "expand_bibi.jls"))]
Slice.sig(ms[1])
seen = Set([length(Slice.sig(ms[1]))])
n = 0
for m in ms
    L = length(Slice.sig(m))
    L in seen && continue
    before = snapshot()
    Slice.sig(m)
    new = setdiff(snapshot(), before)
    push!(seen, L)
    println("new length $L: $(length(new)) instances")
    foreach(mi -> println("    ", first(string(mi.specTypes), 160)), new)
    global n += 1
    n >= 2 && break
end
