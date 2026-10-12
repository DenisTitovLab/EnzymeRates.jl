# ABOUTME: Skeptic P5a: counts :OnlyA mixed groups (binding + chemistry members) in the populations and
# ABOUTME: whether each plain member lies on a cycle of the full graph (a bridge cannot constrain it).
include(joinpath(@__DIR__, "MWC.jl"))
"""Whether step j is a bridge of the undirected step graph (removing it disconnects its forms)."""
function is_bridge(m, j)
    s = m.steps[j]
    adj = Dict{Int, Vector{Int}}()
    for (k, t) in enumerate(m.steps)
        k == j && continue
        push!(get!(adj, t.from, Int[]), t.to); push!(get!(adj, t.to, Int[]), t.from)
    end
    seen = Set([s.from]); q = [s.from]
    while !isempty(q)
        u = popfirst!(q)
        for v in get(adj, u, Int[]); v in seen || (push!(seen, v); push!(q, v)); end
    end
    !(s.to in seen)
end
for name in ("specs", "allo_uniuni", "expand_bibi", "allo_bibi")
    ms = filter(m -> m.allosteric && :OnlyA in m.group_tags, Pop.load_population(name))
    nmix = 0; ncyc = 0; ex = String[]
    for m in ms
        kind = [c.kind for c in MWC.Thermo.canonical_steps(m)]
        hit = false; cyc = false
        for (g, grp) in enumerate(m.groups)
            m.group_tags[g] === :OnlyA || continue
            any(kind[j] != 1 for j in grp) && any(kind[j] == 1 for j in grp) || continue
            hit = true
            for j in grp
                kind[j] == 1 && !is_bridge(m, j) && (cyc = true)
            end
        end
        nmix += hit; ncyc += cyc
        cyc && length(ex) < 5 && push!(ex, m.id)
    end
    println(rpad(name, 12), " :OnlyA mechanisms ", length(ms), "; with a mixed :OnlyA group ", nmix,
            "; with a plain member of a mixed :OnlyA group on a cycle ", ncyc, " ", ex)
end
