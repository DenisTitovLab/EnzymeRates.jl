# ABOUTME: Skeptic P5a: two-free-form cases. Compares the package's rate with deletion oracles that flip
# ABOUTME: at different free forms, and the probe's cone limit with both readings of deletion.
include(joinpath(@__DIR__, "limits_lib.jl"))
using Printf
const Λ0 = MWC.big1024(10)^40
rel(a, b) = Float64(abs(a - b) / max(abs(a), abs(b), MWC.big1024(1e-300)))

"""Forms of the inactive graph as the package builds it: forms that appear in an active step
whose two forms are reachable from a free form, in first-seen order of `m.steps`
(reversed steps, `m.inactive_reversed`, list `to` first)."""
function inactive_forms(m, active, keep)
    out = Int[]
    for (j, s) in enumerate(m.steps)
        (active[j] && s.from in keep && s.to in keep) || continue
        rev = !isempty(m.inactive_reversed) && m.inactive_reversed[j]
        for f in (rev ? (s.to, s.from) : (s.from, s.to)); f in out || push!(out, f); end
    end
    out
end

"""Single-conformation steady state over `forms` with the flip at `free`: (π(free), flux);
a flip form outside `forms` is an inert conformation (π = 1, no flux)."""
function cstate(m, k, re, active, concs, Λ, forms, free)
    free in forms || return MWC.big1024(1), MWC.big1024(0)
    loc = Dict(f => i for (i, f) in enumerate(forms))
    edges = MWC.cat_edges(m, k, re, active, concs, Λ, f -> get(loc, f, nothing))
    length(forms) == 1 && return MWC.big1024(1), MWC.big1024(0)
    x, J = MWC.solve_net(length(forms), edges)
    x[loc[free]], J
end

"""Deletion rate with the inactive flip at `freeI` (:pkg = first free form of the inactive
graph in first-seen order, :E = the mechanism's first free form)."""
function v_del2(m, kA, kI, L, regA, regI, concs, n, freeI; Λ = Λ0)
    setprecision(BigFloat, MWC.PREC) do
        tags = MWC.step_tags(m)
        reA = [s.re for s in m.steps]
        nf = length(m.species)
        fA = MWC.free_form(m)
        πA, JA = cstate(m, kA, reA, trues(length(m.steps)), concs, Λ, collect(1:nf), fA)
        active = [t !== :OnlyA for t in tags]
        keep = MWC.reachable(m, active)
        forms = inactive_forms(m, active, keep)
        fI = freeI === :E ? fA : something(findfirst(f -> Pop.is_free(m.species[f]), forms) |>
                                           (i -> i === nothing ? nothing : forms[i]), fA)
        πI, JI = cstate(m, kI, reA, active, concs, Λ, forms, fI)
        wA = MWC.W_factor(m, regA, concs) / πA^n
        wI = BigFloat(L) * MWC.W_factor(m, regI, concs) / πI^n
        (wA * JA + wI * JI) / (wA + wI), fI
    end
end

for file in ("adv_cases.jls", "adv_cases2.jls")
    ms = deserialize(joinpath(@__DIR__, "data", file))
    for m in ms
        (m.ok && length([sp for sp in m.species if Pop.is_free(sp)]) > 1) || continue
        println("── ", m.id, "  tags ", m.group_tags, "  free forms ",
                [sp.name for sp in m.species if Pop.is_free(sp)])
        # (1) package value vs deletion read two ways, at the package's own samples
        for (si, s) in enumerate(m.samples)
            vals = pkg_values(m, s)
            kA = pkg_step_constants(m, vals, :A)
            kI = [x === nothing ? [big(1), big(1)] : x for x in pkg_step_constants(m, vals, :I)]
            regA, regI = pkg_regK(m, vals)
            c = Dict(k => MWC.big1024(v) for (k, v) in s.concs)
            vp, fp = v_del2(m, kA, kI, vals[:L], regA, regI, c, m.cat_multiplicity, :pkg)
            ve, _ = v_del2(m, kA, kI, vals[:L], regA, regI, c, m.cat_multiplicity, :E)
            si <= 2 && @printf("   sample %d: package v %+.6e | deletion, I flips at %s (inactive graph's first free) %+.6e rel %.1e | deletion, I flips at %s %+.6e rel %.1e\n",
                    si, s.v, m.species[fp].name, Float64(vp), isnan(s.v) ? NaN : rel(vp, big(s.v)),
                    m.species[MWC.free_form(m)].name, Float64(ve), isnan(s.v) ? NaN : rel(ve, big(s.v)))
        end
        # (2) the cone limit at a joint-solve point vs both deletion readings at that point
        for conv in (:i, :iu)
            untie = conv === :iu
            rng = Random.Xoshiro(11)
            J = MWC.joint_solve(m; order = :pkg, untie)
            hat = MWC.random_hat(J, rng); keq, L = logu(rng), logu(rng)
            rA, rI = draw_reg(m, rng)
            concs = Dict(x => logu(rng) for x in m.metabolites)
            dirs, d = MWC.cone_directions(J, :i)
            if isempty(dirs)
                println("   conv $conv: cone EMPTY"); continue
            end
            S2 = sym_setup(J, hat, keq, L, rA, rI, concs)
            vs = [sym_eval(S2, a) for a in dirs[unique(round.(Int, range(1, length(dirs); length = min(16, length(dirs)))))]]
            x = MWC.all_values(J, hat, keq)
            kA = MWC.stored(J, x, :A); kI = MWC.stored(J, x, :I)
            cc = Dict(k => MWC.big1024(v) for (k, v) in concs)
            vp, fp = v_del2(m, kA, kI, L, rA, rI, cc, m.cat_multiplicity, :pkg)
            ve, _ = v_del2(m, kA, kI, L, rA, rI, cc, m.cat_multiplicity, :E)
            spread = maximum(reldiff(s.v, vs[1].v, 0) for s in vs)
            seq = any(s -> s.ordI[2] !== nothing && s.ordI[3] !== nothing && s.ordI[2] < s.ordI[3], vs)
            e = setprecision(() -> BigFloat(10)^(-max(6, 30 ÷ max(1, Int(maximum(abs, dirs[1]))))), BigFloat, MWC.PREC)
            vo = v_path(J, hat, dirs[1], e, keq, L, rA, rI, concs)
            @printf("   conv %-2s: cone d=%d, %d directions, spread %.1e, seq=%s; limit %+.6e (oracle on the path %+.6e, rel %.1e) | deletion flip-%s %+.6e (Δ %.1e) | deletion flip-%s %+.6e (Δ %.1e)\n",
                    conv, d, length(dirs), spread, seq, Float64(vs[1].v), Float64(vo), rel(vo, vs[1].v),
                    m.species[fp].name, Float64(vp), rel(vp, vs[1].v),
                    m.species[MWC.free_form(m)].name, Float64(ve), rel(ve, vs[1].v))
        end
    end
end
