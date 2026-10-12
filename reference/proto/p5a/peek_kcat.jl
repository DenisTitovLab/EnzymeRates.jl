# ABOUTME: Prints the structure and leading faces of mechanisms whose kcat differs (P5a).
# ABOUTME: Usage: julia peek_kcat.jl <population> <index...>
include(joinpath(@__DIR__, "pkgvals.jl"))
function topoly_(p, T)
    d = Dict{Poly.Mono, Int}()
    for (mono, c) in p
        Poly.addto!(d, Poly.mono((Poly.var!(T, s) => e for (s, e) in mono)...), Int(c))
    end
    Poly.P{Int}(d)
end
pop = Pop.load_population(ARGS[1])
for i in parse.(Int, ARGS[2:end])
    m = pop[i]
    println("== ", m.id, " ", m.source, " family=", m.family)
    for s in m.steps; println("   ", Pop.step_text(m, s)); end
    T = Poly.Vars()
    num = topoly_(m.pkg_num, T); den = topoly_(m.pkg_den, T)
    prods = Set(T.index[p] for p in m.products if haskey(T.index, p))
    keep(p) = Poly.P{Int}(Dict(mn => c for (mn, c) in p.t if !any(v -> v.first in prods, mn)))
    n2, d2 = keep(num), keep(den)
    ds = [Dict(T.index[s] => MWC.Q(-1) for s in m.substrates if haskey(T.index, s))]
    on, pn = MWC.initial(n2, ds); od, pd = MWC.initial(d2, ds)
    pr(p) = join(sort([string(c, "·", join([string(T.names[v], "^", e) for (v, e) in mn], "·")) for (mn, c) in p.t]), "  +  ")
    println("   num face order ", on, ": ", first(pr(pn), 600))
    println("   den face order ", od, ": ", first(pr(pd), 600))
end
