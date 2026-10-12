# ABOUTME: Follows the traps the first pass flagged as unclean (an exit other than releasing the
# ABOUTME: substrate or taking up a product) to larger substrate concentrations, dominant forms only.
include(joinpath(@__DIR__, "vlib.jl"))

cases = [("expand_terter", 1695, :A), ("expand_terter", 1697, :D), ("expand_terter", 1700, :A),
         ("expand_terter", 1820, :B)]
rows = deserialize(joinpath(VD, "data", "zero_rows.jls"))
odd = [(x.set, x.index, t.X) for x in rows if x.ok && x.cls === :zero && !x.allo
       for t in x.traps if :other in t.kinds || !t.all_hold_P || !t.all_hold_X || t.share < 0.5]
println("unclean traps from the first pass: ", odd)
recs = Dict((r.set, r.index) => r.m for r in pop_records("expand_terter"))
picks = deserialize(joinpath(KDATA, "terter_sample.jls"))
for set in unique(s for (s, _, _) in odd if s != "expand_terter")
    foreach(r -> recs[(r.set, r.index)] = r.m, pop_records(set, [i for (s, i, _) in odd if s == set]))
end
rng = Xoshiro(5)
for (set, index, X) in unique(odd)
    m = recs[(set, index)]
    ok, why, sol, L = derive_law(m)
    names, subs, prods, mods = roles(m)
    xi = findfirst(==(X), names)
    pd = draw_params(L.vars, rng)
    kv = edge_constants(m, sol, pd)
    base = zeros(length(names))
    foreach(i -> base[i] = 10.0^(-1 + 2 * rand(rng)), subs)
    nl = NL(L, pd)
    println("\n", set, "/", index, " X = ", X)
    for g in (1e3, 1e6, 1e9, 1e12, 1e15)
        c = big.(base); c[xi] *= g
        v, x = oracle(m, kv, c)
        dom = [(S.form_text(m, i), round(Float64(x[i]); sigdigits = 3)) for i in sortperm(x; rev = true)
               if x[i] > 1e-3]
        @printf("  X x %.0e: v oracle %.3e law %.3e; forms above 1e-3: %s\n", g, Float64(v),
                rate(nl, c), string(dom))
    end
    # exits of the dominant form at the largest concentration
    c = big.(base); c[xi] *= 1e15
    _, x = oracle(m, kv, c)
    T = [i for i in eachindex(x) if x[i] > 1e-3]
    for (k, e) in enumerate(m.edges)
        (e.a in T) == (e.b in T) && continue
        println("     exit: ", S.step_text(m, e), " (from ", S.form_text(m, e.a in T ? e.a : e.b), ")")
    end
end
