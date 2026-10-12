# ABOUTME: Checks the bare-constant convention and the edge <-> Pop step constant map against the
# ABOUTME: oracle's thermodynamically consistent constants on every converted record.
const PROTO = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto"
const DATA = joinpath(@__DIR__, "data")
include(joinpath(PROTO, "common", "Oracle.jl"))
include(joinpath(@__DIR__, "..", "src", "Slice.jl"))
using Serialization, Printf, Random
const S = Slice
setprecision(BigFloat, 512)

rng = Xoshiro(11)
nrec, nedge, roundtrip_bad, worst = 0, 0, 0, big(0)
for r in [split(l, '\t') for l in readlines(joinpath(DATA, "MANIFEST.tsv"))[2:end]]
    pops = deserialize(r[4])
    for c in deserialize(joinpath(DATA, r[1] * ".jls"))
        p = pops[c.index]
        m = c.mech
        k, G, mu = Oracle.consistent_constants(p, rng; Keq = 3.0)
        x = S.edge_constants(c, k)
        S.pop_constants(c, x) == k || (global roundtrip_bad += 1)
        # each form's free energy, read through the steps that touch it
        Gm = zeros(BigFloat, length(m.forms))
        for (j, s) in enumerate(p.steps)
            e = m.edges[c.edge_of[j]]
            pa, pb = c.flipped[j] ? (s.to, s.from) : (s.from, s.to)
            Gm[e.a], Gm[e.b] = G[pa], G[pb]
        end
        fwd = S.orient(m)
        μ(pk) = sum((mu[y.name] for y in S.mets_of(m, pk)); init = big(0))
        for kk in eachindex(m.edges)
            o = S.oriented(m, kk, fwd)
            Ka = exp(Gm[o.a] + μ(o.cons) - Gm[o.b] - μ(o.rel))   # association, forward
            got = length(x[kk]) == 2 ? x[kk][1] / x[kk][2] :
                  S.is_binding(o) ? 1 / x[kk][1] : x[kk][1]
            global worst = max(worst, abs(got - Ka) / Ka)
            global nedge += 1
        end
        global nrec += 1
    end
end
@printf("constant map: %d records, %d edges; Pop round trip fails %d; worst relative error \
         of the forward equilibrium ratio (kf/kr, 1/K for a binding, K otherwise) %.1e\n",
        nrec, nedge, roundtrip_bad, Float64(worst))
