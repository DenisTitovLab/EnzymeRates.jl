# ABOUTME: Coverage of the moves' filters (each switched off in turn by a method override in this
# ABOUTME: process only) and the split's time with thermo.jl's cycles fixed to Int indices.
include(joinpath(@__DIR__, "load.jl"))

const SETS = ["init_bibi_pp", "dead1", "expand_bibi", "allo_bibi", "p1_enum_inh"]
const DECLARED = Dict("p1_enum_inh" => [:I])
rows = Dict(set => deserialize(joinpath(S, "slice/enum/data/cmp_$set.jls")) for set in SETS)
declared(set) = get(DECLARED, set, [:I, :A, :Q])

# ── Shapes of the children ────────────────────────────────────────────────────────────────
for set in SETS
    flipped, split, rev, dual = Dict{Int, Int}(), Dict{Int, Int}(), 0, 0
    for r in rows[set]
        m = r.mech
        for c in r.flip.mine
            n = length(unique(m.group[k] for k in eachindex(m.edges)
                              if m.edges[k].re != c.edges[k].re))
            flipped[n] = get(flipped, n, 0) + 1
        end
        for c in r.split.mine
            n = Slice.ngroups(c) - Slice.ngroups(m)
            split[n] = get(split, n, 0) + 1
            rev += any(k -> m.edges[k].re != c.edges[k].re, eachindex(m.edges))
        end
        dual += count(c -> any(x -> Slice.is_dual(c, x), c.mets) >
                           any(x -> Slice.is_dual(m, x), m.mets), r.dead.mine)
    end
    @printf("%-12s flip children by groups flipped %s; split children by groups split %s, \
            with a part sent to RE %d; dead-end children binding a reactant copy %d\n",
            set, sort!(collect(flipped)), sort!(collect(split)), rev, dual)
end

# ── The split with cycles fixed (isassigned with an Int index) ───────────────────────────
src = read(joinpath(Slice.SRC, "thermo.jl"), String)
a = findfirst("function cycles(", src)[1]
b = findnext("\nend\n", src, a)[end]
fixed = replace(src[a:b], "isassigned(path, e.a)" => "isassigned(path, Int(e.a))",
                "isassigned(path, e.b)" => "isassigned(path, Int(e.b))")
occursin("Int(e.a)", fixed) || error("thermo.jl's cycles no longer matches the patch")
Core.eval(Slice, Meta.parse(fixed))
timed(f) = (t = @elapsed r = f(); t < 0.05 && for _ in 1:2; t = min(t, @elapsed f()); end; (r, t))
for set in SETS
    ts, same = Float64[], 0
    for r in rows[set]
        cs, t = timed(() -> Slice.split_children(r.mech))
        push!(ts, t)
        same += cs == r.split.mine
    end
    tp = [r.split.t_pkg for r in rows[set]]
    @printf("%-12s split with fixed cycles: same children %d/%d; median %.1f µs (p90 %.1f, \
            total %.3f s); package median %.1f µs; median ratio %.1f\n", set, same,
            length(ts), 1e6 * median(ts), 1e6 * quantile(ts, 0.9), sum(ts), 1e6 * median(tp),
            median(tp ./ ts))
end

# ── Each filter switched off in turn ──────────────────────────────────────────────────────
function extra(set, mv, f)
    n = 0
    for r in rows[set]
        n += length(f(r.mech)) - length(getfield(r, mv).mine)
    end
    n
end
@eval Slice hyperbolic(m::Mechanism) = true
@printf("children the hyperbolic filter removes (allo_bibi, n = 2): %d\n",
        extra("allo_bibi", :flip, Slice.flip_children))
@eval Slice chain_flanks(m::Mechanism) = Int[]
for set in SETS
    @printf("%-12s flip children the chain-flank rule removes: %d\n", set,
            extra(set, :flip, Slice.flip_children))
end
@eval Slice redundant_copies(m::Mechanism) = Int[]
for set in SETS
    @printf("%-12s children the copy rule removes: split %d, dead end %d\n", set,
            extra(set, :split, Slice.split_children),
            extra(set, :dead, m -> Slice.dead_end_children(m, declared(set))))
end
