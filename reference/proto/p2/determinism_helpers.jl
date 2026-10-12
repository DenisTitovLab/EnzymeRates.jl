# ABOUTME: Helpers shared by the P2 determinism scripts: permuted and regrouped mechanisms, signatures.
# ABOUTME: Included after Thermo.jl (expects const T = Thermo); stdlib only.
using Random
"""`m` with species, step order, group order, member order and step directions shuffled."""
function permuted(m::Pop.PMech, rng)
    ns = length(m.species); π = randperm(rng, ns); inv = invperm(π)   # new i holds old π[i]
    species = m.species[π]
    flip = rand(rng, Bool, length(m.steps))
    σ = randperm(rng, length(m.steps)); σinv = invperm(σ)               # new j holds old σ[j]
    steps = map(σ) do j
        s = m.steps[j]
        flip[j] ? Pop.PStep(inv[s.to], inv[s.from], s.released, s.consumed, s.re) :
                  Pop.PStep(inv[s.from], inv[s.to], s.consumed, s.released, s.re)
    end
    groups = [shuffle(rng, [σinv[j] for j in g]) for g in shuffle(rng, m.groups)]
    Pop.PMech(; id = m.id, source = m.source, substrates = m.substrates, products = m.products,
              inhibitors = m.inhibitors, metabolites = m.metabolites, species = species,
              steps = steps, groups = groups)
end

regrouped(m, groups) = Pop.PMech(; id = m.id, source = m.source, substrates = m.substrates,
    products = m.products, inhibitors = m.inhibitors, metabolites = m.metabolites,
    species = m.species, steps = m.steps, groups = groups)

signature(sol) = (T.independent_names(sol), T.dependent_table(sol),
                  sort!([T.colname(sol, d.col) => d.fold for d in sol.deps]; by = first))
labels(sol) = sort!([T.colname(sol, d.col) => d.label for d in sol.deps]; by = first)

