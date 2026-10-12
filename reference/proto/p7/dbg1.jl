include(joinpath(@__DIR__, "P7.jl"))
ms = Pop.load_population("handwritten")
for id in ("handwritten/21", "handwritten/23")
    m = only(filter(x -> x.id == id, ms))
    println("=== $id $(m.source) family=$(m.family)")
    for s in m.steps; println("   ", Pop.step_text(m, s)); end
    println("groups ", m.groups, " names ", m.step_names)
    println("pass2 ", m.pass2_rename, " deps ", m.dependent_monomials)
    fx = P7.factored_law(m)
    V = fx.vars
    show_node(k) = println("  node $k: ", join([string(t.coef, "*", join(["$(V.names[v])^$e" for (v, e) in t.mono], "*"), " ", join(["n$(f)^$p" for (f, p) in t.facs], "*")) for t in fx.nodes[k]], " + "))
    for k in eachindex(fx.nodes); show_node(k); end
    println(" pkg_den: ", m.pkg_den)
    println(" pkg_num: ", m.pkg_num)
    E = Engine.bare_rate_law(m; check = false)
    println(" sigma: ", length(E.stats[:sigma]), " blocks: ", E.stats[:block_sizes], " attach: ", E.stats[:attach])
end
