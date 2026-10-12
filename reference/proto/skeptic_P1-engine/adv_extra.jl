# ABOUTME: Extra checks on the adversarial package mechanisms: the package's Reduced rate against the oracle,
# ABOUTME: EngineLean against the package where the Hill engine fails, and RE-tree (BFS) order dependence.
include(joinpath(@__DIR__, "SkLib.jl"))
using Random, Printf, Serialization
setprecision(BigFloat, Oracle.PREC)
pop = deserialize(joinpath(@__DIR__, "data", "adv_pkg.jls"))

"""EngineLean law with steps renamed to the package's group names, concentrations to lowest terms."""
function lean_named(m; order = eachindex(m.steps), rename = Dict{Symbol, Symbol}())
    V, sv, num, den, dfree = EngineLean.bare_rate_law(m; order)
    sub = Dict(v => Poly.mono(Poly.var!(V, get(rename, m.step_names[j][i], m.step_names[j][i])) => 1)
               for (j, vs) in enumerate(sv) for (i, v) in enumerate(vs))
    n, d, f = (Poly.subst(p, v -> get(sub, v, nothing)) for p in (num, den, dfree))
    shift = Poly.lowest_shift((n, d), [V.index[x] for x in m.metabolites])
    V, n * shift, d * shift, f * shift
end

for m in pop
    @printf("%-52s", m.source)
    # (1) package Reduced v against the oracle at the package's own dependents
    if m.ok
        errs = Float64[]
        for s in m.samples
            k = Oracle.constants_from_names(m, Pop.all_values(s); rename = Dict(m.pass2_rename))
            push!(errs, Sk.relerr(big(s.v), Oracle.steady_state_v(m, k, Pop.concs(s))))
        end
        @printf(" pkg Reduced v vs oracle %.1e", maximum(errs))
    else
        print(" pkg not ok                     ")
    end
    # (2) lean law vs package raw polynomials
    V, n, d, f = lean_named(m)
    same = (Sk.sym(n, V), Sk.sym(d, V), Sk.sym(f, V)) ==
           (Sk.topoly(m.raw_num), Sk.topoly(m.raw_den), Sk.topoly(m.raw_dfree))
    print(" | lean == pkg raw ", isempty(m.raw_den) ? "n/a" : same)
    # (3) BFS order: reversed step scan, Full mode (group names) via cross-multiplication
    V2, n2, d2, _ = lean_named(m; order = reverse(eachindex(m.steps)))
    conv(p) = Poly.subst(p, v -> Poly.mono(Poly.var!(V, V2.names[v]) => 1))
    full_dep = n * conv(d2) != conv(n2) * d
    print(" | BFS order changes Full law ", full_dep)
    println()
end
