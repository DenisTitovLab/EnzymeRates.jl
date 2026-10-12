# ABOUTME: Compares EngineLean with Engine on every plain mechanism, times it, and runs the reduced pipeline
# ABOUTME: (group names, Pass-2 renames, package dependents, (1 − Γ/Keq) division) against the package's v.
const C = joinpath(@__DIR__, "..", "common")
include(joinpath(C, "Oracle.jl")); include(joinpath(C, "Engine.jl")); include(joinpath(@__DIR__, "EngineLean.jl"))
using Printf, Statistics, Serialization
setprecision(BigFloat, Oracle.PREC)
const POPS = ["specs", "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi",
              "handwritten"]
mechs = [(p, m) for p in POPS for m in Pop.load_population(p) if !m.allosteric]
for (_, m) in mechs[1:40]; EngineLean.bare_rate_law(m); end

same = 0; tl = Float64[]; nfac = npos = 0; verr = Float64[]; halve = 0
for (pop, m) in mechs
    L = Engine.bare_rate_law(m)
    V, sv, num, den, dfree = EngineLean.bare_rate_law(m)
    global same += (Poly.symbolic(num, V), Poly.symbolic(den, V), Poly.symbolic(dfree, V)) ==
                   (Poly.symbolic(L.num, L.vars), Poly.symbolic(L.den, L.vars), Poly.symbolic(L.dfree, L.vars))
    push!(tl, minimum(@elapsed(EngineLean.bare_rate_law(m)) for _ in 1:3))

    # Reduced pipeline: step variables → group names (Pass-2 folded) → package dependents.
    rename = Dict(m.pass2_rename)
    nm(s) = get(rename, s, s)
    sub = Dict(v => Poly.mono(Poly.var!(V, nm(m.step_names[j][i])) => 1)
               for (j, vs) in enumerate(sv) for (i, v) in enumerate(vs))
    dep = Dict(Poly.var!(V, k) => Poly.mono((Poly.var!(V, s) => Int(e) for (s, e) in mono)...)
               for (k, mono) in m.dependent_monomials)
    numr, denr = (Poly.subst(Poly.subst(p, v -> get(sub, v, nothing)), v -> get(dep, v, nothing))
                  for p in (num, den))
    shift = Poly.lowest_shift((numr, denr), [V.index[x] for x in m.metabolites])
    numr, denr = numr * shift, denr * shift
    kv = Poly.var!(V, :Keq)
    x = Poly.mono((Poly.var!(V, p) => 1 for p in m.products)...,
                  (Poly.var!(V, s) => -1 for s in m.substrates)..., kv => -1)
    q = EngineLean.divide_one_minus(numr, x, kv)
    q === nothing && (println("not divisible: ", m.id); continue)
    global nfac += 1
    global npos += all(>(0), values(q.t))
    global halve += 2 * Poly.nterms(q) == Poly.nterms(numr)
    for s in m.samples
        vals = ones(BigFloat, length(V.names))
        for (k, v) in [s.params; s.concs]
            haskey(V.index, k) && (vals[V.index[k]] = big(v))
        end
        xv = Poly.evaluate(Poly.term(Poly.P{Int}, x), vals)
        v = (1 - xv) * Poly.evaluate(q, vals) / Poly.evaluate(denr, vals)
        push!(verr, Float64(abs(v - s.v) / abs(s.v)))
    end
end
n = length(mechs)
println("lean == full engine (num, den, dfree identical): $same of $n")
@printf("lean time per derivation: median %.1f µs, max %.1f µs, total %.3f s\n",
        1e6 * median(tl), 1e6 * maximum(tl), sum(tl))
println("reduced numerator divisible by (1 − Γ/Keq): $nfac of $n; positive quotient $npos; ",
        "quotient has exactly half the numerator's terms in $halve")
@printf("reduced v = (1 − Γ/Keq)·Q/den vs the package's Reduced rate_equation: max rel %.1e over %d samples\n",
        maximum(verr), length(verr))
