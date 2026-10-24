# ABOUTME: Law-derivation time and live bytes of law types C (immutable Poly) and E (mutable Poly)
# ABOUTME: and D (Dict-alias Poly) on the scan population, interleaved per round, best of 5.
const SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
const VD = joinpath(SP, "phase45/research5/verify")
include(joinpath(SP, "phase45/research5/law/variants/C/SliceC.jl"))
include(joinpath(VD, "variantE/SliceE.jl"))
include(joinpath(VD, "variantD/SliceD.jl"))
include(joinpath(VD, "sliceside.jl"))
raws = [eval(Meta.parse(l)) for l in eachline(joinpath(VD, "data", "scan_pop.tsv"))]
function build(S)
    out = S.Mechanism[]
    for r in raws
        try
            push!(out, mkmech(S, r))
        catch
        end
    end
    out
end
pops = Dict(:C => (SliceC, build(SliceC)), :E => (SliceE, build(SliceE)), :D => (SliceD, build(SliceD)))
function derive_all(S, ms, sols)
    n = 0
    for (m, sol) in zip(ms, sols)
        L = S.rate_law(m, sol)
        n += L isa S.Law ? length(L.D isa AbstractDict ? L.D : L.D.terms) : 1
    end
    n
end
sols = Dict(k => [S.solve(m) for m in ms] for (k, (S, ms)) in pops)
best = Dict(k => Inf for k in keys(pops))
for rd in 1:5, k in (:C, :E, :D)
    S, ms = pops[k]
    GC.gc()
    t = @elapsed derive_all(S, ms, sols[k])
    best[k] = min(best[k], t)
end
for k in (:C, :E, :D)
    println(k, ": ", length(pops[k][2]), " laws, best ", round(1e6 * best[k] / length(pops[k][2]); digits = 1),
            " µs per law (", round(best[k] / best[:C]; digits = 3), " of C)")
end
