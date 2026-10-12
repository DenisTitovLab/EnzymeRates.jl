# ABOUTME: Counts Haldane-labelled dependents per plain mechanism (P2 probe); stdlib only.
# ABOUTME: Throwaway script.
include(joinpath("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common", "Thermo.jl"))
h = Dict{Int, Int}(); ex = String[]
for p in Pop.populations(), m in Pop.load_population(p)
    m.allosteric && continue
    sol = Thermo.solve(m)
    n = count(d -> d.label === :haldane, sol.deps)
    h[n] = get(h, n, 0) + 1
    n > 1 && length(ex) < 6 && push!(ex, m.id * " " * m.source)
end
println("Haldane-labelled dependents per mechanism: ", sort(collect(h)))
foreach(println, ex)
