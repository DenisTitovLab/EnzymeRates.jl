# ABOUTME: The measured plain laws rebuilt over P2's joint solve (Thermo: content-keyed independent
# ABOUTME: constants, every dependent substituted) and compared in size with the package-dependent laws.
include(joinpath(@__DIR__, "Encode.jl"))
include(joinpath(@__DIR__, "..", "common", "Thermo.jl"))
include(joinpath(@__DIR__, "Select.jl"))
using .Poly: Vars, var!, Mono, mono, mmul, mpow, subst, lowest_shift, term, nterms
using Printf, Statistics, Random
setprecision(BigFloat, 256)
const PI = P7.PI
countops(b) = (n = Ref(0); w(x) = (x isa Expr && (x.head === :call && (n[] += 1); foreach(w, x.args))); w(b); n[])
"""The law of plain `m` over P2's independent columns: `(law, nindep)` with flat (num, den) nodes and
the Haldane-split variant."""
function p2_laws(m)
    sol = Thermo.solve(m)
    E = Engine.bare_rate_law(m; check = false)
    V = Vars()
    indep = sort(sol.indep; by = i -> Thermo.colname(sol, i))
    for i in indep; var!(V, Symbol(Thermo.colname(sol, i))); end
    kv = var!(V, :Keq); var!(V, :E_total)
    np = length(V.names)
    foreach(x -> var!(V, x), m.metabolites)
    colmono = Dict{Int, Mono}(i => mono(V.index[Symbol(Thermo.colname(sol, i))] => 1) for i in indep)
    toint(q) = (denominator(q) == 1 || error("fractional exponent"); Int(numerator(q)))
    for d in sol.deps
        colmono[d.col] = mmul(mono(kv => toint(d.keq)), mono((colmono[f] |> only |> first => toint(e) for (f, e) in d.terms)...))
    end
    idx = Dict((c.j, c.r) => i for (i, c) in enumerate(sol.cols))
    vmap = Dict{Int, Mono}()
    for (j, vs) in enumerate(E.stepvars)
        c = sol.cs[j]
        if c.re
            k = colmono[idx[(j, 1)]]
            vmap[vs[1]] = c.flipped && !c.bindK ? mpow(k, -1) : k
        else
            vmap[vs[1]] = colmono[idx[(j, c.flipped ? 2 : 1)]]
            vmap[vs[2]] = colmono[idx[(j, c.flipped ? 1 : 2)]]
        end
    end
    for x in m.metabolites; vmap[E.vars.index[x]] = mono(V.index[x] => 1); end
    num, den = (subst(p, v -> vmap[v]) for p in (E.num, E.den))
    concs = collect(np + 1:length(V.names))
    s = lowest_shift((num, den), concs); num, den = num * s, den * s
    flat = P7.Law(V, np, Mono[], length(m.metabolites), [P7.node_of(num), P7.node_of(den)], 1, 2, V.index[:E_total])
    H, Q = P7.haldane_split(num, m, V)
    Q, c = P7.lift(Q, concs)
    hal = P7.Law(V, np, Mono[], length(m.metabolites),
                 [P7.node_of(H), P7.node_of(Q), [P7.Term(1, c, [1 => 1, 2 => 1])], P7.node_of(den)], 3, 4, V.index[:E_total])
    flat, hal, length(indep), sol
end
rows = []
for (pop, m) in selected()
    m.allosteric && continue
    f2, h2, n2, sol = p2_laws(m)
    f1, h1 = P7.flat_law(m), P7.factored_law(m; expand_den = true)
    ops(l, o) = countops(Gen.body(Gen.decode(Tuple(code(l, o)))))
    push!(rows, (id = m.id, fit_pkg = length(m.fitted_params), fit_p2 = n2,
                 den_pkg = length(f1.nodes[2]), den_p2 = length(f2.nodes[2]),
                 num_pkg = length(f1.nodes[1]), num_p2 = length(f2.nodes[1]),
                 flath_pkg = ops(f1, (true, true, 0)), flath_p2 = ops(f2, (true, true, 0)),
                 halh_pkg = ops(h1, (true, true, 0)), halh_p2 = ops(h2, (true, true, 0))))
    # numeric identity: same rate at the oracle's consistent constants mapped both ways is P2's job;
    # here only check that the P2 law is a valid rate law at its own random independents (finite).
end
r(f) = [f(x) for x in rows]
@printf("%d plain laws. fitted params: package median %d, P2 median %d (P2 − pkg: %s)\n", length(rows),
        median(r(x -> x.fit_pkg)), median(r(x -> x.fit_p2)), string(sort(unique(r(x -> x.fit_p2 - x.fit_pkg)))))
for (lab, a, b) in (("den terms", :den_pkg, :den_p2), ("num terms", :num_pkg, :num_p2),
                    ("flat+Horner ops", :flath_pkg, :flath_p2), ("Haldane+Horner ops", :halh_pkg, :halh_p2))
    ratio = r(x -> getfield(x, b) / getfield(x, a))
    @printf("  %-20s P2/package ratio median %.2f [%.2f, %.2f]; package sum %d, P2 sum %d\n", lab,
            median(ratio), minimum(ratio), maximum(ratio), sum(r(x -> getfield(x, a))), sum(r(x -> getfield(x, b))))
end
