# ABOUTME: Prints examples of the :OnlyA admission rejections (W) by kind, with the admission
# ABOUTME: message, and the share of W candidates that are also D or P; usage: wexamples.jl <name>
using Serialization, Printf
const CEN = @__DIR__
kind(msg) = occursin("has none", msg) ? "no free form" :
            occursin("flips at one free form", msg) ? "second free form" :
            occursin("still join the two sides", msg) ? "G1: arc inside a component" :
            occursin("close a cycle", msg) ? "G1: arcs close a cycle" :
            occursin("reach the free enzyme only", msg) ? "island" : "other"
for name in split(ARGS[1], ",")
    d = deserialize(joinpath(CEN, "data", "$name.jls"))
    fitted = [r for r in d.recs if r[:status] === :ok && r[:fitted]]
    W = [r for r in fitted if get(r, :admit, nothing) !== nothing]
    println("== $name: $(length(W)) :OnlyA rejections among $(length(fitted)) fitted")
    byk = Dict{String, Vector{Any}}()
    for r in W
        push!(get!(byk, kind(r[:admit]), Any[]), r)
    end
    for (k, rs) in sort!(collect(byk); by = first)
        println("  [$k] $(length(rs)); rounds ", join(("r$x $(count(r -> r[:round] == x, rs))" for x in 0:2), ", "),
                "; tags: OnlyA groups per candidate median ",
                sort([count("::OnlyA", r[:text]) for r in rs])[cld(length(rs), 2)])
        for r in rs[1:min(end, 2)]
            println("    #$(r[:gid]) r$(r[:round]) $(r[:move]) np=$(r[:np]) rank=$(r[:rank_pkg])\n      $(r[:text])\n      admit: $(first(r[:admit], 300))")
        end
    end
    # do the rejected candidates use chemistry tags only (V-type) or bindings
    pp = count(r -> occursin("_res_", r[:text]), W)
    println("  ping-pong (residual forms) among them: $pp")
end
