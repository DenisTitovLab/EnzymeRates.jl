# ABOUTME: Skeptic P5a: the probe's V4 check and light limit study (conv i, ii, iu) on the extra
# ABOUTME: populations exported by skp_export.jl; prints the probe's summary counts per population.
include(joinpath(@__DIR__, "limits_lib.jl"))
using Printf
const Λ = MWC.big1024(10)^40
relerr_(a, b) = abs(a - b) / max(abs(b), MWC.big1024(1e-300))
function inputs(m, s)
    vals = pkg_values(m, s)
    kA = pkg_step_constants(m, vals, :A)
    kI = pkg_step_constants(m, vals, :I)
    tags = MWC.step_tags(m)
    activeI = [tags[j] !== :OnlyA && kI[j] !== nothing for j in eachindex(m.steps)]
    kIf = [kI[j] === nothing ? [big(1), big(1)] : kI[j] for j in eachindex(m.steps)]
    keepI = MWC.reachable(m, [tags[j] !== :OnlyA for j in eachindex(m.steps)])
    regA, regI = pkg_regK(m, vals)
    concs = Dict(k => MWC.big1024(v) for (k, v) in s.concs)
    vals, kA, kIf, activeI, keepI, regA, regI, concs
end
function light(m, conv; seed = 11, ncone = 16)
    rng = Random.Xoshiro(seed)
    untie = conv === :iu
    c = untie ? :i : conv
    J = MWC.joint_solve(m; order = :pkg, untie)
    hat = MWC.random_hat(J, rng); keq, L = logu(rng), logu(rng)
    rA, rI = draw_reg(m, rng)
    concs = Dict(x => logu(rng) for x in m.metabolites)
    r = Dict{Symbol, Any}(:id => m.id, :conv => conv)
    v1 = v_path(J, hat, zeros(MWC.Q, length(J.cols)), big(1), keq, L, rA, rI, concs)
    r[:v_del] = try v_deletion(J, hat, keq, L, rA, rI, concs) catch; big(NaN) end
    dirs, d = MWC.cone_directions(J, c)
    r[:cone_dim], r[:cone_n] = d, length(dirs)
    isempty(dirs) && return r
    S = sym_setup(J, hat, keq, L, rA, rI, concs)
    pick = unique([1; round.(Int, range(1, length(dirs); length = min(ncone, length(dirs))))])
    res = [sym_eval(S, dirs[k]) for k in pick]
    r[:cone_v] = res[1].v
    r[:cone_spread] = maximum(reldiff(s.v, res[1].v, 0) for s in res)
    r[:cone_seq] = any(s -> s.ordI[2] !== nothing && s.ordI[3] !== nothing && s.ordI[2] < s.ordI[3], res)
    r[:cone_cancel] = any(s -> s.cancel, res)
    r[:cone_vs_del] = maximum(reldiff(s.v, r[:v_del], 0) for s in res)
    e = setprecision(() -> BigFloat(10)^(-max(6, 30 ÷ max(1, Int(maximum(abs, dirs[1]))))), BigFloat, MWC.PREC)
    r[:oracle_err] = try reldiff(v_path(J, hat, dirs[1], e, keq, L, rA, rI, concs), res[1].v, v1) catch; NaN end
    r
end
for name in ARGS
    f = joinpath(@__DIR__, "data", "$name.jls")
    isfile(f) || (println("missing $name"); continue)
    pop = deserialize(f)
    println("== $name: $(length(pop)) mechanisms; families: ",
            join(("$k=$(count(m -> k in m.family, pop))" for k in (:pingpong, :inhibitor, :reg_sites, :fused_binding, :theorell_chance, :multi_metabolite_side, :dead_end_leaf, :tag_NonequalAI)), ", "))
    # V4
    worst = 0.0; nbad = 0; nerr = 0; bad = String[]
    for m in pop
        m.ok || continue
        e = 0.0
        for s in m.samples
            try
                vals, kA, kI, act, keep, rA, rI, c = inputs(m, s)
                v = MWC.v_factorized(m, kA, kI, [st.re for st in m.steps], act, keep, vals[:L], rA, rI, c;
                                     n = m.cat_multiplicity, Λ = Λ)
                e = max(e, Float64(relerr_(v, MWC.big1024(s.v))))
            catch err
                e = Inf
            end
        end
        isinf(e) && (nerr += 1)
        worst = max(worst, isinf(e) ? 0.0 : e)
        e > 1e-10 && (nbad += 1; length(bad) < 6 && push!(bad, @sprintf("%s(%.1e)", m.id, e)))
    end
    @printf("  V4 deletion oracle vs package v: worst finite rel %.1e; above 1e-10: %d (oracle errors %d) %s\n",
            worst, nbad, nerr, join(bad, " "))
    recs = Dict{Symbol, Any}[]
    t = @elapsed for m in pop, conv in (:i, :ii, :iu)
        r = try light(m, conv) catch e
            Dict{Symbol, Any}(:id => m.id, :conv => conv, :error => first(sprint(showerror, e), 160))
        end
        r[:admitted] = m.ok; push!(recs, r)
    end
    serialize(joinpath(@__DIR__, "data", "res_$name.jls"), recs)
    for conv in (:i, :ii, :iu)
        rs = filter(r -> r[:conv] === conv, recs)
        ok = filter(r -> !haskey(r, :error), rs)
        feas = filter(r -> haskey(r, :cone_v), ok)
        uniq = filter(r -> r[:cone_spread] < 1e-20, feas)
        udel = filter(r -> r[:cone_vs_del] < 1e-20, uniq)
        @printf("  conv %-2s: %d mechanisms, errors %d; cone empty %d; unique %d (= deletion %d, ≠ deletion %d); non-unique %d; sequestration possible %d; face cancellation %d; oracle err max %.1e\n",
                conv, length(rs), length(rs) - length(ok), length(ok) - length(feas), length(uniq), length(udel),
                length(uniq) - length(udel), length(feas) - length(uniq), count(r -> r[:cone_seq], feas),
                count(r -> r[:cone_cancel], feas), maximum((r[:oracle_err] for r in feas if !isnan(r[:oracle_err])); init = 0.0))
        for r in filter(r -> haskey(r, :error), rs)[1:min(end, 3)]; println("     error ", r[:id], ": ", r[:error]); end
        em = [r[:id] for r in ok if !haskey(r, :cone_v)]
        isempty(em) || println("     cone empty: ", join(em[1:min(end, 10)], " "), length(em) > 10 ? " …" : "")
        g = [r[:id] for r in uniq if r[:cone_vs_del] >= 1e-20 && !r[:cone_seq]]
        isempty(g) || println("     unique, no sequestration, ≠ deletion (G1–G3 pass): ", join(g[1:min(end, 10)], " "))
    end
    @printf("  (%.0f s)\n", t)
end
