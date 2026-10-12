# ABOUTME: Lists plain mechanisms whose flux-free SS steps form whole kinetic groups (true phantoms).
# ABOUTME: Throwaway P2 probe script; stdlib only.
include(joinpath("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common", "Thermo.jl"))
const T = Thermo
n_any = 0; n_whole = 0; whole = String[]
for pop in Pop.populations(), m in Pop.load_population(pop)
    m.allosteric && continue
    z = T.zero_flux_ss(m)
    any(z) || continue
    global n_any += 1
    wg = [g for g in m.groups if all(j -> !m.steps[j].re && z[j], g)]
    if !isempty(wg)
        global n_whole += 1
        push!(whole, string(m.id, " [", m.source, "] whole flux-free SS groups: ",
              join([join([Pop.step_text(m, m.steps[j]) for j in g], " | ") for g in wg], " ; ")))
    end
end
println("mechanisms with a flux-free SS step: $n_any; with a whole flux-free SS group: $n_whole")
foreach(println, whole)
