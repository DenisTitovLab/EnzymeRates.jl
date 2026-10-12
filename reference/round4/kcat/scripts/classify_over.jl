# ABOUTME: Cross-tabulates sup > 1.01 kcat in the plain laws against abortive complexes (a form
# ABOUTME: holding a substrate with a product, or a substrate its residual already holds) on SS cycles.
include(joinpath(@__DIR__, "klib.jl"))

rows = [x for l in ("core", "expand_terter", "terter_sample")
        for x in deserialize(joinpath(KDATA, "rows_$l.jls")) if x.ok && !x.allo && x.cls === :finite]
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

tab = Dict{Tuple{Bool, Bool, Bool}, Vector{Int}}()
examples = Dict{Tuple{Bool, Bool, Bool}, Any}()
for x in rows
    m = recs[(x.set, x.index)]
    core = S.catalytic_core(m)
    sup = nt_support(m, core)
    ss_on = any(k -> sup[k] && !m.edges[k].re, eachindex(m.edges))
    ab_on = any(k -> sup[k] && (abortive(m, m.edges[k].a) || abortive(m, m.edges[k].b)), eachindex(m.edges))
    ss_ab = any(k -> sup[k] && !m.edges[k].re && (abortive(m, m.edges[k].a) || abortive(m, m.edges[k].b)),
                eachindex(m.edges))
    k = (ss_on, ab_on, ss_ab)
    v = get!(tab, k, [0, 0])
    v[1] += 1
    v[2] += x.ratio > 1.01
    x.ratio > 1.01 && !haskey(examples, k) && (examples[k] = (x, m))
end
println("plain finite laws by (SS step on a non-turnover cycle, abortive form on such a cycle, SS step at an abortive form on such a cycle): [laws, sup > 1.01 kcat]")
foreach(((k, v),) -> println("  ", k, " => ", v), sort!(collect(tab); by = string))
for (k, (x, m)) in sort!(collect(examples); by = string)
    println("  example of ", k, ": ", x.id, " (", x.set, "/", x.index, ") ratio ", round(x.ratio; digits = 3))
    foreach(l -> println("      ", l), step_lines(m))
end
