# ABOUTME: For mechanisms whose Hill quotient carries Keq, tries every column as the first-eliminated
# ABOUTME: one and reports whether some content-based pivot gives a Keq-free quotient; P4 probe.
using Serialization
include(joinpath(@__DIR__, "P4.jl"))
rows = deserialize(joinpath(@__DIR__, "census.jls"))
byid = Dict{String, Any}()
for pop in Pop.populations(), m in Pop.load_population(pop)
    m.allosteric || (byid[m.id] = m)
end
bad = [r[:id] for r in rows if !r[:hill_keqfree]]
fixable, how = 0, Dict{Any, Int}()
unfixable = String[]
for id in bad
    m = byid[id]
    L = Engine.bare_rate_law(m; check = false)
    sol0 = Thermo.solve(m)
    found = false
    for c in sol0.cols
        NL = try P4.law(m; front = [c.name], L) catch; continue end
        h = P4.hill(NL)
        if h.keq_free
            global fixable += 1; found = true
            hd = [P4.depclass(NL.sol, d) for d in NL.sol.deps if d.keq != 0]
            k = (P4.depclass(NL.sol, only(d for d in NL.sol.deps if NL.sol.cols[d.col].name == c.name)))
            how[k] = get(how, k, 0) + 1
            break
        end
    end
    found || push!(unfixable, id)
end
println("Keq inside Q: ", length(bad), "; a single front column makes Q Keq-free in ", fixable)
println("class of that column (re, kind, against): ", how)
println("no single column works: ", length(unfixable), " ", first(unfixable, 20))
