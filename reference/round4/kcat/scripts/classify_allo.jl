# ABOUTME: Cross-tabulates sup > 1.01 kcat in the allosteric laws against SS steps at abortive forms
# ABOUTME: on non-turnover cycles, the multiplicity and conformation-dependent substrate binding.
include(joinpath(@__DIR__, "klib.jl"))

rows = [x for l in ("core", "expand_terter", "terter_sample")
        for x in deserialize(joinpath(KDATA, "rows_$l.jls")) if x.ok && x.allo && x.cls === :finite]
need = Set((x.set, x.index) for x in rows)
recs = Dict{Tuple{String, Int}, Any}()
for set in unique(first.(need))
    if set in CORE_SETS || set == "expand_terter"
        foreach(r -> (r.set, r.c.index) in need && (recs[(r.set, r.c.index)] = r.c.mech), records_of(set))
    else
        ps = deserialize(joinpath(IDATA, "pop_$set.jls"))
        for (s, i) in need
            s == set && (recs[(s, i)] = S.Converted(ps[i], i).mech)
        end
        ps = nothing
        GC.gc()
    end
end

"""Whether form `i` is an abortive complex: it holds a substrate (or its modifier copy) and a
product, or a substrate that its residual already added."""
function abortive(m, i)
    f = m.forms[i]
    b = [x.name for x in S.mets_of(m, f.bound)]
    res = [x.name for x in S.mets_of(m, f.res)]
    any(in(m.substrates), b) && (any(in(m.products), b) || any(x -> x in res, filter(in(m.substrates), b)))
end

"""Support of the non-turnover cycle space (as `features`), restricted to `mask`."""
function nt_support(m, mask)
    zs, cs = S.cycles(m, S.orient(m))
    ref = findfirst(!=(0), cs)
    sup = falses(length(m.edges))
    for (z, c) in zip(zs, cs)
        v = c == 0 ? z : ref === nothing || z === zs[ref] ? nothing : cs[ref] .* z .- c .* zs[ref]
        v === nothing || (sup .|= v .!= 0)
    end
    sup .& mask
end

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
    m = recs[(x.set, x.index)]
    core = S.catalytic_core(m)
    sup = nt_support(m, core)
    ss_ab = any(k -> sup[k] && !m.edges[k].re && (abortive(m, m.edges[k].a) || abortive(m, m.edges[k].b)),
                eachindex(m.edges))
    f = allo_flags(m)
    k = (n = f.n, ss_abortive_cycle = ss_ab, conf_dependent = f.subtag || f.sschem)
    v = get!(tab, k, [0, 0])
    v[1] += 1
    v[2] += x.ratio > 1.01
    if x.ratio > 1.01 && !ss_ab && !(f.subtag || f.sschem)
        println("  neither flag: ", x.id, " (", x.set, "/", x.index, ") ratio ", round(x.ratio; digits = 3),
                " (substrates only ", round(x.ratios; digits = 3), ") modifiers ", x.mods, " sites ", m.sites)
        foreach(l -> println("      ", l), step_lines(m))
    end
end
println("allosteric finite laws (core sets, expand_terter) by (multiplicity, SS step at an abortive form on a non-turnover cycle, substrate binding conformation-dependent or SS with conformation-dependent chemistry): [laws, sup > 1.01 kcat]")
foreach(((k, v),) -> println("  ", k, " => ", v), sort!(collect(tab); by = string))
