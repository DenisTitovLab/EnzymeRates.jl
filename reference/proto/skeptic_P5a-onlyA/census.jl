# ABOUTME: Census of the allosteric populations by tag pattern, multiplicity, sites and step kinds (P5a).
# ABOUTME: Prints counts used to choose the probe's cases.
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
using Printf
for name in ("specs", "handwritten", "allo_uniuni", "allo_bibi", "expand_bibi")
    pop = filter(m -> m.allosteric, Pop.load_population(name))
    nonly = count(m -> :OnlyA in m.group_tags, pop)
    ns = Dict{Int,Int}(); for m in pop; ns[m.cat_multiplicity] = get(ns, m.cat_multiplicity, 0) + 1; end
    nsites = count(m -> !isempty(m.sites), pop)
    # :OnlyA on RE chemistry, on SS binding, on RE binding, on SS chemistry
    cnt = Dict{Symbol,Int}()
    kinds = Dict{Symbol,Int}()
    for m in pop
        gos = Pop.group_of_step(m)
        seen = Set{Symbol}()
        for (j, s) in enumerate(m.steps)
            m.group_tags[gos[j]] === :OnlyA || continue
            fs, ts = m.species[s.from], m.species[s.to]
            plain = Pop.is_binding(s) && sort([string(x.name, x.role) for x in [fs.bound; s.consumed]]) == sort([string(x.name, x.role) for x in ts.bound]) && fs.res_added == ts.res_added && fs.res_subtracted == ts.res_subtracted
            k = Symbol(s.re ? "RE_" : "SS_", plain ? "binding" : "chem")
            push!(seen, k)
        end
        for k in seen; cnt[k] = get(cnt, k, 0) + 1; end
        kind = isempty(seen) ? :none : (any(k -> endswith(string(k), "binding"), seen) ? (any(k -> endswith(string(k), "chem"), seen) ? :K_and_V : :K_only) : :V_only)
        kinds[kind] = get(kinds, kind, 0) + 1
    end
    @printf("%-12s allo=%4d OnlyA=%4d n=%s sites=%d  OnlyA step kinds(mechs)=%s  class=%s\n", name, length(pop), nonly, ns, nsites, cnt, kinds)
    noA = [m.id for m in pop if !(:OnlyA in m.group_tags)]
    println("   no-OnlyA ids: ", noA)
end
