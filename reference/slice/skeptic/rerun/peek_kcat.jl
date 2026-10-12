# ABOUTME: Shows one ter-ter law whose kcat depends on how the substrates grow: its steps, the
# ABOUTME: denominator's leading face and the ratio on each face monomial; usage: peek_kcat.jl
include(joinpath(@__DIR__, "lib.jl"))
rows = deserialize(joinpath(IDATA, "pipe_terter2.jls"))
bad = filter(x -> x.ok && x.kcat isa String && occursin("grow relative", x.kcat), rows)
x = first(sort(bad; by = x -> x.nD))
r = only(filter(r -> r.c.index == x.index, records_of("init_terter.part2")))
m = r.c.mech
println(x.id, " (", length(bad), " such laws in part 2; this one has ", x.nD, " D terms)")
foreach(k -> println("  ", S.step_text(m, m.edges[k]), "   group ", m.group[k]), eachindex(m.edges))
core = S.catalytic_core(m)
println("  dead-end edges removed: ", count(!, core))
L = S.law(m, S.solve(m), core)
names = L.vars.names
subs = [L.vars.index[s] for s in m.substrates]
N = L.Q * S.conc_mono(L.vars, m.substrates)
keep(p) = S.split_conc(S.Poly(filter(t -> !any(v -> v.first <= length(S.metabolites(m)) &&
                                         !(v.first in subs), t.first), p.terms)), L.vars)
n, d = keep(N), keep(L.D)
deg(mo) = sum(e for (v, e) in mo if v in subs; init = 0 // 1)
top = maximum(deg, keys(d))
for mo in sort(collect(keys(d)); by = string)
    deg(mo) == top || continue
    println("  face ", join(("$(names[v])^$e" for (v, e) in mo), " "), ": numerator ",
            haskey(n, mo) ? S.sum_text(n[mo], names) : "0", "  over  ", S.sum_text(d[mo], names))
end
