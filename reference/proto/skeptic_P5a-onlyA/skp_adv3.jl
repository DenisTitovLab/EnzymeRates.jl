# ABOUTME: Skeptic P5a: the skp_adv.jl checks (V4 and the limit study, conv i, ii, iu) on the third batch
# ABOUTME: on the adversarial mechanisms built by skp_build.jl, and prints a verdict table.
include(joinpath(@__DIR__, "limits_lib.jl"))
using Printf

const Λ = MWC.big1024(10)^40
relerr_(a, b) = abs(a - b) / max(abs(b), MWC.big1024(1e-300))

"""validate.jl's oracle inputs of a package sample."""
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

"""run_pop.jl's light study (cone directions, no ε sweep)."""
function light(m, conv; seed = 11, ncone = 16)
    rng = Random.Xoshiro(seed)
    untie = conv === :iu
    conv = untie ? :i : conv
    J = MWC.joint_solve(m; order = :pkg, untie)
    hat = MWC.random_hat(J, rng)
    keq, L = logu(rng), logu(rng)
    rA, rI = draw_reg(m, rng)
    concs = Dict(x => logu(rng) for x in m.metabolites)
    r = Dict{Symbol, Any}(:id => m.id, :conv => untie ? :iu : conv)
    S = sym_setup(J, hat, keq, L, rA, rI, concs)
    v1 = v_path(J, hat, zeros(MWC.Q, length(J.cols)), big(1), keq, L, rA, rI, concs)
    r[:v_del] = try v_deletion(J, hat, keq, L, rA, rI, concs) catch e; "ERR: " * first(sprint(showerror, e), 80) end
    dirs, d = MWC.cone_directions(J, conv)
    r[:cone_dim], r[:cone_n] = d, length(dirs)
    if !isempty(dirs)
        pick = unique([1; round.(Int, range(1, length(dirs); length = min(ncone, length(dirs))))])
        res = [sym_eval(S, dirs[k]) for k in pick]
        r[:cone_v] = res[1].v
        r[:cone_spread] = maximum(reldiff(s.v, res[1].v, 0) for s in res)
        r[:cone_seq] = any(s -> s.ordI[2] !== nothing && s.ordI[3] !== nothing &&
                                s.ordI[2] < s.ordI[3], res)
        r[:cone_vs_del] = r[:v_del] isa String ? NaN :
                          maximum(reldiff(s.v, r[:v_del], 0) for s in res)
        e = setprecision(() -> BigFloat(10)^(-max(6, 30 ÷ max(1, Int(maximum(abs, dirs[1]))))),
                         BigFloat, MWC.PREC)
        r[:oracle_err] = try
            reldiff(v_path(J, hat, dirs[1], e, keq, L, rA, rI, concs), res[1].v, v1)
        catch err
            NaN
        end
        r[:dirs] = dirs[pick]; r[:J] = J; r[:hat] = hat; r[:S] = S
    end
    r
end

verdict(r) = haskey(r, :error) ? "ERROR " * r[:error] :
    !haskey(r, :cone_v) ? "cone EMPTY" :
    r[:cone_spread] >= 1e-20 ? @sprintf("NON-UNIQUE(spread %.1e)%s", r[:cone_spread], r[:cone_seq] ? "+seq" : "") :
    (r[:cone_vs_del] isa Number && r[:cone_vs_del] < 1e-20) ? "unique = deletion" :
    @sprintf("unique ≠ deletion (Δ %.1e)%s", r[:cone_vs_del], r[:cone_seq] ? "+seq" : "")

ms = deserialize(joinpath(@__DIR__, "data", "adv_cases3.jls"))
println("== V4/V1 on the adversarial cases: deletion oracle (package semantics) vs package v")
for m in ms
    m.ok || continue
    e = 0.0; note = ""
    for s in m.samples
        try
            vals, kA, kI, act, keep, rA, rI, c = inputs(m, s)
            reI = [st.re for st in m.steps]
            v = MWC.v_factorized(m, kA, kI, reI, act, keep, vals[:L], rA, rI, c;
                                 n = m.cat_multiplicity, Λ = Λ)
            e = max(e, isnan(s.v) ? NaN : Float64(relerr_(v, MWC.big1024(s.v))))
            isnan(s.v) && (note = @sprintf(" package v NaN; oracle v = %.4g", Float64(v)))
        catch err
            note = " oracle error: " * first(sprint(showerror, err), 120); e = NaN
        end
    end
    @printf("  %-28s free forms %-8s  first free (Pop order) %-3s  rel %.1e%s\n", m.id,
            join([sp.name for sp in m.species if Pop.is_free(sp)], ","),
            m.species[MWC.free_form(m)].name, e, note)
end

println("\n== limit study (conv i, ii with ε sweep via study(); iu via light())")
recs = Dict{Symbol, Any}[]
for m in ms, conv in (:i, :ii, :iu)
    r = try
        conv === :iu ? light(m, :iu) : study(m, conv; seed = 7)
    catch e
        Dict{Symbol, Any}(:id => m.id, :conv => conv, :error => first(sprint(showerror, e), 200))
    end
    r[:conv] = conv; r[:admitted] = m.ok
    push!(recs, r)
    extra = conv === :iu || haskey(r, :error) ? "" :
        @sprintf(" | RREF: partner=%d wrong=%d cons=%.0e sym-vs-oracle=%.0e orderSpread=%.1e cone-oracle=%.0e",
                 length(r[:partner]), length(r[:wrong]), r[:consistency], r[:sym_vs_oracle], r[:order_spread],
                 get(r, :cone_oracle_err, NaN))
    @printf("  %-28s %-3s %-4s %s%s\n", m.id, conv, m.ok ? "adm" : "REJ", verdict(r), extra)
end
serialize(joinpath(@__DIR__, "data", "adv_results3.jls"),
          [Dict(k => v for (k, v) in r if !(k in (:J, :S))) for r in recs])
