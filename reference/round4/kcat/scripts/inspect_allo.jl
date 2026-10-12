# ABOUTME: Cross-tabulates sup > 1.01 kcat in the allosteric core laws against how substrates bind
# ABOUTME: the two conformations (tags of substrate-binding and chemistry steps) and the multiplicity.
include(joinpath(@__DIR__, "klib.jl"))

rows = filter(x -> x.ok && x.allo && x.cls === :finite, deserialize(joinpath(KDATA, "rows_core.jls")))
recs = Dict((r.set, r.c.index) => r for set in CORE_SETS for r in records_of(set))

"""Flags of an allosteric mechanism: a substrate-binding step in a :NonequalAI or :OnlyA group
(`subtag`); a steady-state substrate-binding step together with chemistry in such a group
(`sschem`); regulatory sites (`sites`)."""
function allo_flags(m)
    fwd = S.orient(m)
    names(p) = S.mets_of(m, p)
    takes_sub(o) = any(x -> x.name in m.substrates && !x.modifier, names(o.cons))
    chem(o) = !any(x -> x.name in m.products, names(o.cons)) &&
              (any(x -> x.name in m.products, names(o.rel)) ||
               (holds_product(m, o.b) && !holds_product(m, o.a)))
    tag(k) = m.tags[m.group[k]]
    os = [S.oriented(m, k, fwd) for k in eachindex(m.edges)]
    subtag = any(k -> takes_sub(os[k]) && tag(k) in (:NonequalAI, :OnlyA), eachindex(os))
    sschem = any(k -> takes_sub(os[k]) && !os[k].re, eachindex(os)) &&
             any(k -> chem(os[k]) && tag(k) in (:NonequalAI, :OnlyA), eachindex(os))
    (n = m.multiplicity, subtag, sschem, sites = !isempty(m.sites))
end

tab = Dict{Any, Vector{Int}}()
for x in rows
    f = allo_flags(recs[(x.set, x.index)].c.mech)
    k = (f.n, f.subtag, f.sschem)
    v = get!(tab, k, [0, 0])
    v[1] += 1
    v[2] += x.ratio > 1.01
end
println("allosteric finite laws by (multiplicity, substrate binding conformation-dependent, SS substrate binding with conformation-dependent chemistry): [laws, sup > 1.01 kcat]")
foreach(((k, v),) -> println("  ", k, " => ", v), sort!(collect(tab); by = string))
