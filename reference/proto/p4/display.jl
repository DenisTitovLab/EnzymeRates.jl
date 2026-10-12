# ABOUTME: P4 parts 2-3 per family: normalized display forms after the joint solve, the checks
# ABOUTME: (division-free evaluation, display = evaluated = oracle, units) and the Hill quotient.
include(joinpath(@__DIR__, "P4.jl"))
include(joinpath(@__DIR__, "Build.jl"))
include(joinpath(@__DIR__, "..", "common", "Oracle.jl"))
using Random, Printf
using .Poly: Mono, mono, mmul, mpow, mexp, P, nterms, evaluate, term, lowest_shift

const PI = P{Int}
setprecision(BigFloat, 512)
rng = Xoshiro(5)
logu(r) = exp(log(1e-2) + rand(r) * (log(1e2) - log(1e-2)))

byid(id) = (pi = split(id, "/"); Pop.load_population(pi[1])[parse(Int, pi[2])])
re_random = Build.mech("re_random_bibi", [:A, :B], [:P, :Q],
    ["E + A ⇌ E(A)", "E + B ⇌ E(B)", "E(A) + B ⇌ E(A,B)", "E(B) + A ⇌ E(A,B)",
     "E(A,B) <--> E(P,Q)", "E(P,Q) ⇌ E(Q) + P", "E(P,Q) ⇌ E(P) + Q", "E(P) ⇌ E + P",
     "E(Q) ⇌ E + Q"])
families = [
    "SS uni-uni (Uni-Uni)" => byid("specs/1"),
    "SS uni-uni (Segel IX-8)" => byid("specs/2"),
    "RE uni-uni" => byid("specs/18"),
    "iso uni-uni (Segel IX-45)" => byid("specs/3"),
    "ordered bi-bi SS (Segel IX-87)" => byid("specs/5"),
    "ordered bi-bi RE bindings" => byid("specs/19"),
    "random bi-bi SS" => byid("specs/11"),
    "random bi-bi RE (textbook, built)" => re_random,
    "random bi-bi RE on free E (spec)" => byid("specs/20"),
    "Theorell-Chance (Segel IX-122)" => byid("specs/6"),
    "ping-pong, conformation F (Segel IX-140)" => byid("specs/7"),
    "ping-pong, residual, SS chemistry (handwritten/11)" => byid("handwritten/11"),
    "ping-pong, residual, RE (handwritten/23)" => byid("handwritten/23"),
    "competitive inhibitor" => byid("specs/21"),
    "non-competitive inhibitor" => byid("specs/22"),
    "ordered ter-ter (Segel IX-261)" => byid("specs/9"),
    "hexa uni ping-pong (Segel IX-308)" => byid("specs/17"),
]

"""The reference monomial of rule A: the constant term when present, else ΠS, else the lowest
monomial in graded order (degree, then text)."""
function reference(L, D)
    haskey(D, Mono()) && return Mono(), "1"
    S = P4.cmono(L, L.m.substrates)
    haskey(D, S) && return S, "ΠS"
    first(sort(collect(keys(D)); by = mo -> (P4.degree(mo), P4.fmt_mono(mo, L.vars)))), "graded"
end

"""Display of `p/(c x^m*)` grouped by concentration monomial; each coefficient is `c_m/c` with
common monomial content and common linear factors cancelled."""
function display(L, p, c::PI, ms::Mono)
    W = L.vars
    parts = String[]
    for (mo, k) in sort(collect(P4.cmap(L, p)); by = x -> (P4.degree(x[1]), P4.fmt_mono(x[1], W)))
        x = mmul(mo, mpow(ms, -1))
        xs = isempty(x) ? "" : P4.fmt_mono(x, W)
        a, b = P4.reduce_ratio(k, c)
        sa = P4.fmt_factored(a, W); sb = P4.fmt_factored(b, W)
        coef = b == one(PI) ? sa : sa * "/(" * sb * ")"
        push!(parts, isempty(xs) ? coef : coef == "1" ? xs : coef * "*" * xs)
    end
    replace(join(parts, " + "), "+ -" => "- ")
end

for (label, m) in families
    L = P4.law(m)
    W = L.vars
    D, N = P4.cmap(L, L.den), P4.cmap(L, L.num)
    ms, rule = reference(L, D)
    c = D[ms]
    cr = P4.criteria(m)
    println("\n=== ", label, "  [", m.id, "]  forms ", length(m.species), ", steps ", length(m.steps),
            " (RE ", count(s -> s.re, m.steps), "), terminal SCCs (c3) ", cr.c3,
            ", reference ", rule, " = ", P4.fmt_mono(ms, W), ", coefficient terms ", nterms(c))
    for d in L.sol.deps
        d.label === :haldane || d.label === :wegscheider || continue
        println("  ", d.label, ": ", L.sol.cols[d.col].name, " (class ", P4.depclass(L.sol, d), ")")
    end
    dd = display(L, L.den, c, ms); nd = display(L, L.num, c, ms)
    println("  v/E_t = [", length(nd) > 400 ? first(nd, 400) * " …" : nd, "]")
    println("        / [", length(dd) > 700 ? first(dd, 700) * " …" : dd, "]")
    # (a) evaluated form = (N/c, D/c): concentrations with non-negative exponents only
    nonneg = all(mexp(t, v) >= 0 for p in (L.num, L.den) for t in keys(p.t) for v in values(L.conc))
    finite_all_zero = haskey(D, Mono())
    finite_one_zero = all(x -> any(mo -> mexp(mo, L.conc[x]) == 0, keys(D)), m.metabolites)
    # (b) display (Laurent, divided by c x^m*) = evaluated (divided by c) = oracle
    worst = big(0)
    for _ in 1:3
        yv = Dict(v => logu(rng) for v in values(L.colvar)); keq = logu(rng)
        concs = Dict(x => logu(rng) for x in m.metabolites)
        vals = zeros(BigFloat, length(W.names))
        for (v, y) in yv; vals[v] = y; end
        vals[L.keq] = keq
        for (x, cx) in concs; vals[L.conc[x]] = cx; end
        cv = evaluate(c, vals); xs = evaluate(term(PI, ms), vals)
        v_eval = (evaluate(L.num, vals) / cv) / (evaluate(L.den, vals) / cv)
        v_disp = sum(evaluate(k, vals) / cv * evaluate(term(PI, mmul(mo, mpow(ms, -1))), vals) for (mo, k) in N) /
                 sum(evaluate(k, vals) / cv * evaluate(term(PI, mmul(mo, mpow(ms, -1))), vals) for (mo, k) in D)
        cols = Dict(i => BigFloat(yv[L.colvar[i]])^L.colsgn[i] for i in L.sol.indep)
        k = Thermo.stored_constants(L.sol, Thermo.all_values(L.sol, cols, keq))
        v_or = Oracle.steady_state_v(m, k, concs)
        worst = max(worst, abs(v_disp - v_eval) / abs(v_eval), abs(v_eval - v_or) / abs(v_or))
    end
    # (c) units of each display coefficient: (conc exponent, time exponent)
    ud = Set((P4.unit(L, t) .+ (P4.degree(mo), 0)) .- (P4.unit(L, first(keys(c.t))) .+ (P4.degree(ms), 0))
             for (mo, k) in D for t in keys(k.t))
    un = Set((P4.unit(L, t) .+ (P4.degree(mo), 0)) .- (P4.unit(L, first(keys(c.t))) .+ (P4.degree(ms), 0))
             for (mo, k) in N for t in keys(k.t))
    @printf("  checks: (a) division-free %s, finite at all-zero %s, finite with any one zero %s; (b) max rel diff %.1e; (c) display-term units D %s N %s; homogeneous %s\n",
            nonneg, finite_all_zero, finite_one_zero, Float64(worst), ud, un, P4.homogeneous(L))
    h = P4.hill(L)
    qs = h.Q === nothing ? "none" : P4.fmt_poly(h.Q, W)
    println("  Hill: N = (ΠS − ΠP/Keq)·Q, Q ", h.Q === nothing ? "absent" : "terms $(nterms(h.Q))",
            ", polynomial ", h.poly, ", positive ", h.positive, ", Keq-free ", h.keq_free,
            "; Q = ", length(qs) > 300 ? first(qs, 300) * " …" : qs)
end
