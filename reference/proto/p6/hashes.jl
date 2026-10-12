# ABOUTME: Writes hash(m) of every population mechanism as a P6 GMech, plus control hashes, so two
# ABOUTME: separate processes can be compared. Usage: julia hashes.jl <out-file> <fwd|rev>
include(joinpath(@__DIR__, "harness.jl"))

struct _testhelper_Ctrl
    x::Int
end
module _testhelper_Tmp
struct C
    x::Int
end
end

out, order = ARGS[1], ARGS[2]
if order == "rev"
    # Different process history: intern many unrelated symbols first, build in reverse order.
    rng = Xoshiro(round(Int, time()))
    syms = [Symbol(randstring(rng, 8)) for _ in 1:20_000]
end
ms = _testhelper_load_all()
order == "rev" && reverse!(ms)
lines = String[]
for m in ms
    g = _testhelper_gm(_testhelper_toraw(m))
    push!(lines, string(m.id, " ", hash(g)))
end
sort!(lines)
push!(lines, string("control hash(:A) ", hash(:A)))
push!(lines, string("control hash(GM.Met(:A, false)) ", hash(GM.Met(:A, false))))
push!(lines, string("control default hash of a Main struct ", hash(_testhelper_Ctrl(1))))
push!(lines, string("control default hash of a runtime-module struct ", hash(_testhelper_Tmp.C(1))))
push!(lines, string("control objectid(GM.Form) ", objectid(GM.Form)))
write(out, join(lines, "\n") * "\n")
println("wrote ", length(ms), " hashes to ", out)
