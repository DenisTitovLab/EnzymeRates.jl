# ABOUTME: Details of the key P5a cases: per-order RREF-path limits with their partner columns,
# ABOUTME: cone limits, deletion, inactive orders, and kcat of the limit law by the same engine.
include(joinpath(@__DIR__, "limits_lib.jl"))

recs = deserialize(joinpath(@__DIR__, "data", "results_cases.jls"))
byid = Dict((r[:id], r[:conv]) => r for r in recs)
f(x) = @sprintf("%.4g", Float64(x))
for id in ("uni/RE/RE/EqualAI,NonequalAI,OnlyA", "uni/RE/RE/OnlyA,EqualAI,EqualAI",
           "uni/RE/RE/OnlyA,OnlyA,EqualAI", "familyA_cut", "multiA", "pp_err1", "pp_dead",
           "pp_err1_ss", "pp_re_chem", "handwritten/34", "specs/41"), conv in (:i, :ii)
    r = byid[(id, conv)]
    println("── $id  conv ($conv)  ", r[:admitted] ? "admitted" : "REJECTED by the guard")
    println("   RREF path (pkg order): v_lim = ", f(r[:v_lim_sym]), "  partner moves ", r[:partner],
            "  wrong-sign ", r[:wrong], "  ord(N_I,Q_I,D_I) = ", [o === nothing ? "-" : string(o[1]) for o in r[:ordI]])
    for (o, p) in sort!(collect(r[:per_order]); by = x -> string(x[1]))
        println("      order ", rpad(string(o), 10), " v_lim = ", rpad(f(p.v), 12), " partner ", p.partner,
                " wrong ", p.wrong)
    end
    println("   deletion v = ", f(r[:v_del]))
    if haskey(r, :cone)
        println("   cone d=", r[:cone_dim], ", ", r[:cone_n], " feasible directions; picked:")
        for c in r[:cone]
            println("      a = ", c.a, "  v_lim = ", f(c.v), "  ord(N_I,Q_I,D_I) = ",
                    [o === nothing ? "-" : string(o[1]) for o in c.ordI], c.seq ? "  SEQUESTERS" : "")
        end
    else
        println("   cone EMPTY (no path moves only the :OnlyA constants)")
    end
end

# kcat of the limit law: the same leading-order routine with three stages
# (ε → 0, then the corner variables, then the substrates), against the oracle.
println("\n── kcat of the conv (i) limit law by the leading-order engine vs the oracle")
cases = Dict(m.id => retag(m) for m in deserialize(joinpath(@__DIR__, "data", "cases.jls")))
for id in ("uni/RE/RE/OnlyA,OnlyA,EqualAI", "uni/SS/SS/EqualAI,OnlyA,NonequalAI", "familyA_cut",
           "multiA", "pp_err1", "pp_err1_ss", "pp_re_chem")
    m = cases[id]
    rng = Random.Xoshiro(3)
    J = MWC.joint_solve(m)
    dirs, _ = MWC.cone_directions(J, :i)
    a = dirs[1]
    hat = MWC.random_hat(J, rng); keq, L = logu(rng), logu(rng)
    rA, rI = draw_reg(m, rng)
    tab = MWC.Table(Poly.Vars())
    lawA = MWC.state_law(J, :A, tab); lawI = MWC.state_law(J, :I, tab)
    T = tab.T
    Lv = Poly.var!(T, :L)
    one_ = one(Poly.P{Int})
    num, den = MWC.combine(lawA, lawI, one_, one_, Poly.term(Poly.P{Int}, Poly.mono(Lv => 1)),
                           m.cat_multiplicity)
    prods = Set(T.index[p] for p in m.products if haskey(T.index, p))
    drop(p) = Poly.P{Int}(Dict(mn => c for (mn, c) in p.t if !any(v -> v.first in prods, mn)))
    mapx(e) = e[1] === :leaf ? (:leaf, drop(e[2])) : (e[1], Any[mapx(x) for x in e[2]])
    dε = Dict{Int, MWC.Q}(); vals = fill(big(1.0), length(T.names))
    for i in J.indep
        k = get(T.index, Symbol(J.cols[i].name), 0); k == 0 && continue
        dε[k] = a[i]; vals[k] = hat[i]
    end
    vals[T.index[:Keq]] = keq; vals[Lv] = L
    dS = Dict(T.index[s] => MWC.Q(-1) for s in m.substrates if haskey(T.index, s))
    r = MWC.lead_value(mapx(num), mapx(den), [dε, dS], vals)
    # oracle: tiny ε, huge substrates, zero products
    c = Dict{Symbol, BigFloat}(x => (x in m.substrates ? big(10)^12 : big(0)) for x in m.metabolites)
    e = setprecision(() -> BigFloat(10)^(-max(8, 40 ÷ max(1, Int(maximum(abs, a))))), BigFloat, MWC.PREC)
    vo = v_path(J, hat, a, e, keq, L, rA, rI, c)
    @printf("   %-36s kcat(engine) = %-12s oracle(ε small, S = 1e12, P = 0) = %-12s rel %.1e\n",
            id, f(r.v), f(vo), reldiff(r.v, vo, vo))
end
