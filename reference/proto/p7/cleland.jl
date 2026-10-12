# ABOUTME: Size of the concentration-monomial (Cleland coefficient) form of the measured laws: how many
# ABOUTME: distinct concentration monomials remain per row if parameter-only coefficients were hoisted.
include(joinpath(@__DIR__, "Encode.jl"))
include(joinpath(@__DIR__, "Select.jl"))
using Printf, Statistics
rows = []
for (pop, m) in selected()
    m.allosteric && continue
    law = P7.factored_law(m; expand_den = true)
    nc0 = law.np + length(law.deps)
    cm(node) = length(unique([[(v, e) for (v, e) in t.mono if v > nc0] for t in node]))
    den = law.nodes[law.den]; numQ = law.nodes[law.num - 1]
    push!(rows, (m.id, length(den), cm(den), length(numQ), cm(numQ)))
end
@printf("plain laws %d: den terms median %d max %d -> distinct conc monomials median %d max %d; Q terms median %d max %d -> conc monomials median %d max %d\n",
        length(rows), median(r[2] for r in rows), maximum(r[2] for r in rows), median(r[3] for r in rows), maximum(r[3] for r in rows),
        median(r[4] for r in rows), maximum(r[4] for r in rows), median(r[5] for r in rows), maximum(r[5] for r in rows))
for r in rows
    r[2] >= 150 && println("  ", r)
end
