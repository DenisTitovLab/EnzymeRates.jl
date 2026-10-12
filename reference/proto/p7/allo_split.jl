# ABOUTME: How often the Haldane factor comes out of the whole MWC numerator (P7 probe).
include(joinpath(@__DIR__, "P7.jl"))
n = 0; ok = 0; deadI = 0
for pop in Pop.populations(), m in Pop.load_population(pop)
    (m.allosteric && m.ok) || continue
    global n += 1
    law = P7.factored_law(m)
    t = only(law.nodes[law.num])
    H = length(t.facs) == 2 && length(law.nodes[t.facs[1].first]) == 2 &&
        any(x -> x.coef < 0, law.nodes[t.facs[1].first])
    global ok += H
    global deadI += isempty(m.state_polys_I[1])
end
println("allosteric laws: $n; Haldane factor out of the MWC numerator: $ok; inactive numerator 0 (dead I-state): $deadI")
