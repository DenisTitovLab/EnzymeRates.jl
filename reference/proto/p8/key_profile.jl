# ABOUTME: Splits the P8 key cost into normalization, refinement and encoding+hash on a few laws
# ABOUTME: (the largest plain laws and a median one), minimum over repeated runs.
include(joinpath(@__DIR__, "Laws.jl"))
using Printf
ms = filter(m -> !m.allosteric, Pop.load_population("expand_bibi"))
bs = [(m.id, Laws.plain_redesign(m).law) for m in ms]
sort!(bs; by = x -> -sum(length, x[2].comps))
picks = [bs[1], bs[2], bs[div(end, 2)]]
for (id, law) in picks
    Canon.canonical(law)
    tmin(f) = minimum(@elapsed(f()) for _ in 1:5) * 1e6
    tn = tmin(() -> Canon.normalize(law))
    tc = tmin(() -> Canon.canonical(law))
    comps = Canon.normalize(law)
    present = sort!(unique(s for c in comps for m in keys(c) for (s, _) in m if !(s in law.fixed)))
    pidx = Dict(p => i for (i, p) in enumerate(present))
    labels, inc = String[], Vector{Tuple{Int, Int}}[]
    for (ci, c) in enumerate(comps), (m, coef) in c
        push!(labels, string(ci, "|", coef, "|", join(("$s^$e" for (s, e) in m if s in law.fixed), "*")))
        push!(inc, [(pidx[s], e) for (s, e) in m if !(s in law.fixed)])
    end
    occ = [Tuple{Int, Int}[] for _ in present]
    for (t, ps) in enumerate(inc), (p, e) in ps; push!(occ[p], (t, e)); end
    tr = tmin(() -> Canon.refine(ones(Int, length(present)), Canon.ranks(labels), inc, occ))
    pcol, _ = Canon.refine(ones(Int, length(present)), Canon.ranks(labels), inc, occ)
    te = tmin(() -> Canon.bytes2hex(Canon.sha256(Canon.encode(pcol, inc, labels))))
    @printf("%-16s terms %4d params %2d | total %7.0f µs = normalize %6.0f + labels/incidence %6.0f + refine %6.0f + encode+sha %6.0f\n",
            id, length(inc), length(present), tc, tn, tc - tn - tr - te, tr, te)
end
