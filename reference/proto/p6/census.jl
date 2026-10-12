# ABOUTME: Size census of every population for the P6 graph-type encoding limits.
# ABOUTME: Counts distinct metabolite roles, multiplicities, forms, steps, groups and sites.
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
using Printf

for pop in Pop.populations()
    ms = Pop.load_population(pop)
    nmet = Int[]; mult = Int[]; nsp = Int[]; nst = Int[]; ngr = Int[]; nconf = Int[]
    nres = Int[]; side = Int[]; nsite = Int[]; nlig = Int[]
    for m in ms
        mets = Set{Tuple{Symbol, Bool}}()
        mx = 0
        for sp in m.species
            for x in sp.bound; push!(mets, (x.name, x.role === :Inh)); end
            for x in sp.res_added; push!(mets, (x, false)); end
            for x in sp.res_subtracted; push!(mets, (x, false)); end
            for x in unique(sp.bound); mx = max(mx, count(==(x), sp.bound)); end
            push!(nres, length(sp.res_added) + length(sp.res_subtracted))
        end
        for s in m.steps, v in (s.consumed, s.released)
            for x in v; push!(mets, (x.name, x.role === :Inh)); end
            for x in unique(v); mx = max(mx, count(==(x), v)); end
            push!(side, length(v))
        end
        push!(nmet, length(mets)); push!(mult, mx); push!(nsp, length(m.species))
        push!(nst, length(m.steps)); push!(ngr, length(m.groups))
        push!(nconf, length(unique(sp.conformation for sp in m.species)))
        push!(nsite, length(m.sites)); push!(nlig, isempty(m.sites) ? 0 : maximum(length(s.ligands) for s in m.sites))
    end
    @printf("%-12s n=%4d  mets max %2d  mult max %d  forms max %2d  steps max %2d  groups max %2d  confs max %d  residual max %d  side max %d  sites max %d  ligands/site max %d\n",
            pop, length(ms), maximum(nmet), maximum(mult), maximum(nsp), maximum(nst), maximum(ngr),
            maximum(nconf), maximum(nres), maximum(side), maximum(nsite), maximum(nlig))
end
