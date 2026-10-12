# ABOUTME: Skeptic P5a: (1) whether "unique" cone limits depend on the hat values of the vanishing
# ABOUTME: constants; (2) the population verdicts re-drawn at an independent second point (seed 12).
include(joinpath(@__DIR__, "limits_lib.jl"))
using Printf

"""Cone limit at a point drawn with `seed`, along every picked direction; optionally at a second
hat point that moves only the :OnlyA inactive constants (`shift` along the last direction)."""
function probe(m, conv; seed = 11, ncone = 16, shift = false)
    untie = conv === :iu
    J = MWC.joint_solve(m; order = :pkg, untie)
    rng = Random.Xoshiro(seed)
    hat = MWC.random_hat(J, rng); keq, L = logu(rng), logu(rng)
    rA, rI = draw_reg(m, rng)
    concs = Dict(x => logu(rng) for x in m.metabolites)
    dirs, d = MWC.cone_directions(J, :i)
    isempty(dirs) && return (empty = true,)
    pick = unique([1; round.(Int, range(1, length(dirs); length = min(ncone, length(dirs))))])
    S1 = sym_setup(J, hat, keq, L, rA, rI, concs)
    vs = [sym_eval(S1, dirs[k]) for k in pick]
    vdel = try v_deletion(J, hat, keq, L, rA, rI, concs) catch; big(NaN) end
    out = (empty = false, v = vs[1].v, spread = maximum(reldiff(s.v, vs[1].v, 0) for s in vs),
           seq = any(s -> s.ordI[2] !== nothing && s.ordI[3] !== nothing && s.ordI[2] < s.ordI[3], vs),
           del = maximum(reldiff(s.v, vdel, 0) for s in vs), hatdep = NaN, deldep = NaN)
    shift || return out
    # move only the :OnlyA inactive constants: values along the last direction at ε0
    x = MWC.path_values(J, hat, dirs[pick[end]], MWC.big1024(0.05), keq)
    hat2 = Dict(i => x[i] for i in J.indep)
    S2 = sym_setup(J, hat2, keq, L, rA, rI, concs)
    v2 = sym_eval(S2, dirs[1]).v
    vdel2 = try v_deletion(J, hat2, keq, L, rA, rI, concs) catch; big(NaN) end
    merge(out, (hatdep = reldiff(v2, vs[1].v, 0), deldep = reldiff(vdel2, vdel, 0)))
end

recs = Dict(name => deserialize(joinpath(@__DIR__, "data", "pop_$name.jls")) for name in ("expand_bibi", "allo_bibi"))
pops = Dict(name => Dict(m.id => m for m in Pop.load_population(name)) for name in ("expand_bibi", "allo_bibi"))
println("== (1) hat dependence of the (i) 'unique ≠ deletion' limits (and a control of = deletion ones)")
for name in ("expand_bibi", "allo_bibi")
    rs = filter(r -> r[:conv] === :i && haskey(r, :cone_v) && r[:cone_spread] < 1e-20, recs[name])
    bad = filter(r -> r[:cone_vs_del] >= 1e-20, rs)
    good = filter(r -> r[:cone_vs_del] < 1e-20, rs)[1:20]
    for (label, set) in (("unique≠del", bad), ("unique=del (control, first 20)", good))
        hd = Float64[]; dd = Float64[]; ex = String[]
        for r in set
            p = probe(pops[name][r[:id]], :i; shift = true)
            push!(hd, p.hatdep); push!(dd, p.deldep)
            p.hatdep > 1e-20 && length(ex) < 6 && push!(ex, @sprintf("%s(%.1e)", r[:id], p.hatdep))
        end
        @printf("  %-12s %-32s n=%3d: limit moves with the vanishing constants' hat values in %d (max %.1e); deletion moves in %d %s\n",
                name, label, length(set), count(>(1e-20), hd), maximum(hd; init = 0.0), count(>(1e-20), dd), join(ex, " "))
    end
end
adv = deserialize(joinpath(@__DIR__, "data", "adv_cases.jls"))
for id in ("AD11_mixed_two_plain",)
    m = only(filter(x -> x.id == id, adv))
    p = probe(m, :i; shift = true)
    @printf("  %s (i): limit vs deletion %.1e; limit moves with vanishing hats %.1e; deletion moves %.1e\n", id, p.del, p.hatdep, p.deldep)
end

println("\n== (2) verdicts re-drawn at an independent point (seed 12) vs the probe's (seed 11)")
for name in ("expand_bibi", "allo_bibi"), conv in (:i, :iu)
    rs = filter(r -> r[:conv] === conv && haskey(r, :cone_v), recs[name])
    cls(spread, del, seq) = spread >= 1e-20 ? "nonunique" : del < 1e-20 ? "=del" : "≠del"
    flips = String[]; n = 0
    t = @elapsed for r in rs
        p = probe(pops[name][r[:id]], conv; seed = 12)
        n += 1
        c1 = cls(r[:cone_spread], r[:cone_vs_del], r[:cone_seq]); c2 = p.empty ? "empty" : cls(p.spread, p.del, p.seq)
        c1 == c2 || push!(flips, "$(r[:id]):$c1→$c2")
    end
    @printf("  %-12s conv %-2s: %d records re-drawn; verdict changed in %d %s (%.0f s)\n", name, conv, n,
            length(flips), join(flips[1:min(end, 8)], " "), t)
end
