# ABOUTME: How often distinct mechanisms produce identical positional laws (P7 probe): per population
# ABOUTME: and across the plain bi-bi populations, against the package's eq_hash classes.
include(joinpath(@__DIR__, "Encode.jl"))
using Printf
classes(keys) = (d = Dict{Any, Vector{Int}}(); for (i, k) in enumerate(keys); push!(get!(d, k, Int[]), i); end; d)
function report(label, ms)
    flat = [code(P7.flat_law(m), (false, false, 0)) for m in ms]
    fact = [code(P7.factored_law(m), (true, false, 0)) for m in ms]
    eqh = [m.eq_hash for m in ms]
    cf, cx, ce = classes(flat), classes(fact), classes(eqh)
    merged = count(v -> length(unique(eqh[v])) > 1, values(cf))       # positional class spanning eq_hash classes
    split = count(v -> length(unique(flat[v])) > 1, values(ce))        # eq_hash class spanning positional classes
    shared = sum(length(v) for v in values(cf) if length(v) > 1; init = 0)
    @printf("%-28s n=%4d eq_hash=%4d positional(flat)=%4d positional(factored)=%4d  mechanisms in a shared class=%4d (%.0f%%)  pos-classes merging eq_hash classes=%d  eq_hash classes split=%d  largest class=%d\n",
            label, length(ms), length(ce), length(cf), length(cx), shared, 100 * shared / length(ms), merged, split,
            maximum(length, values(cf)))
    for v in sort(collect(values(cf)); by = length, rev = true)[1:min(3, end)]
        length(v) > 1 && println("    e.g. ", join((ms[i].id for i in v[1:min(6, end)]), ", "), length(v) > 6 ? ", …" : "")
    end
    cf
end
allms = Dict(p => [m for m in Pop.load_population(p) if m.ok] for p in Pop.populations())
for p in sort(collect(keys(allms)))
    report(p, allms[p])
end
plain_bibi = [m for p in ("init_bibi", "init_bibi_pp", "expand_bibi") for m in allms[p] if !m.allosteric]
report("plain bi-bi (3 populations)", plain_bibi)
allo_bibi = [m for p in ("expand_bibi", "allo_bibi") for m in allms[p] if m.allosteric]
report("allosteric bi-bi (2 pops)", allo_bibi)
