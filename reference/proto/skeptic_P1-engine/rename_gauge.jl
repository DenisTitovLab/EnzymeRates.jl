# ABOUTME: Shows that the engine's polynomials depend on metabolite names through the segment-root tie-break:
# ABOUTME: one mechanism with its product named P or T gives different (num, den) for the same rate function.
include(joinpath(@__DIR__, "SkLib.jl"))
using Printf

F(name, mets...) = Pop.PSpecies([Pop.PMet(x, r) for (x, r) in mets], :E, Symbol[], Symbol[], name)
M(x, r = :S) = Pop.PMet(x, r)
"""Uni-uni with an RE isomerization between the two bound forms; `p` names the product."""
function mech(p)
    bp = Symbol("E", p)
    forms = [F(:E), F(:ES, (:S, :S)), F(bp, (p, :P))]
    steps = [Pop.PStep(1, 2, [M(:S)], [], false), Pop.PStep(2, 3, [], [], true),
             Pop.PStep(3, 1, [], [M(p, :P)], false)]
    Pop.PMech(; id = "rename_$p", source = "rename_$p", substrates = [:S], products = [p],
              metabolites = [:S, p], species = forms, steps = steps, groups = [[1], [2], [3]],
              step_names = [[:kf1, :kr1], [:K2], [:kf3, :kr3]])
end
laws = Dict{Symbol, Any}()
for p in (:P, :T)
    m = mech(p)
    L = Engine.bare_rate_law(m)
    seg, _, _ = Sk.segment_graph(m)
    roots = [argmin(f -> (length(m.species[f].bound), string(m.species[f].name)),
                    [x for x in eachindex(m.species) if seg[x] == seg[2]])]
    # Rename the product variable to :P so both laws live in the same symbols.
    V = L.vars
    sub = Dict(V.index[p] => Poly.mono(Poly.var!(V, :P) => 1))
    num, den = (Poly.subst(q, v -> get(sub, v, nothing)) for q in (L.num, L.den))
    laws[p] = (Sk.sym(num, V), Sk.sym(den, V))
    println("product named $p: root of the RE segment {ES, E$p} = ", m.species[roots[1]].name)
    println("   den = ", sort([join(["$(s)^$(e)" for (s, e) in k], "·") for k in keys(laws[p][2])]))
end
same = laws[:P] == laws[:T]
println("identical (num, den) after renaming T → P: ", same)
