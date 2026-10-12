# ABOUTME: Follow-up B1: convention Dka = untied R_D with the ratio kept (kf_I, kr_I same exponent) on every
# ABOUTME: kinetic door that thermodynamics does not force; run where kinetic cycle doors exist. Usage: julia dka_eval.jl pop...
include(joinpath(@__DIR__, "full_eval.jl"))

"""Full-direction verdict under Dka."""
function full_eval_dka(m; maxdirs = 3000)
    rng = Random.Xoshiro(11)
    J = MWC.joint_solve(m; order = :pkg, untie = true)
    hat = MWC.random_hat(J, rng); keq, L = logu(rng), logu(rng)
    rA, rI = draw_reg(m, rng)
    concs = Dict(x => logu(rng) for x in m.metabolites)
    a = analyze(m, :D)
    keep = [j for j in a.doors if !a.aff[j] && !(j in a.forced)]
    extra = [(J.idx[(:I, j, 1)], J.idx[(:I, j, 2)]) for j in keep]
    r = Dict{Symbol, Any}(:id => m.id, :conv => :Dka, :nextra => length(extra))
    dirs, d = cone_dirs(J, :ii; extra)
    r[:cone_dim], r[:cone_n] = d, length(dirs)
    r[:v_del] = v_deletion(J, hat, keq, L, rA, rI, concs)
    if !isempty(dirs)
        Sx = sym_setup(J, hat, keq, L, rA, rI, concs)
        res = [sym_eval(Sx, x) for x in dirs[1:min(end, maxdirs)]]
        vs = [s.v for s in res]
        r[:cone_v] = vs[1]
        r[:cone_spread] = maximum(reldiff(v, vs[1], 0) for v in vs)
        r[:cone_seq] = any(s -> s.ordI[2] !== nothing && s.ordI[3] !== nothing && s.ordI[2] < s.ordI[3], res)
        r[:cone_vs_del] = maximum(reldiff(v, r[:v_del], 0) for v in vs)
        r[:n_eval] = length(vs)
    end
    r
end

recs = Dict{Symbol, Any}[]
for name in ARGS
    pop = filter(m -> :OnlyA in m.group_tags, load_allo(name))
    sel = filter(m -> (a = analyze(m, :D); a.cyc_main || a.vanish_cycle || a.cyc_hang), pop)
    tab = Dict{Tuple{String, String}, Int}()
    t = @elapsed for m in sel
        r = try full_eval_dka(m) catch e
            Dict{Symbol, Any}(:id => m.id, :conv => :Dka, :error => first(sprint(showerror, e), 200))
        end
        push!(recs, r)
        k = (analyze(m, :D).cls, verdict(r)); tab[k] = get(tab, k, 0) + 1
    end
    @printf("%s: %d mechanisms with kinetic cycle doors (%.0f s)\n", name, length(sel), t)
    for ((c, v), n) in sort!(collect(tab); by = first); @printf("   %-12s %-8s %4d\n", c, v, n); end
end
serialize(joinpath(@__DIR__, "data", "dka_" * join(ARGS, "+") * ".jls"), recs)
