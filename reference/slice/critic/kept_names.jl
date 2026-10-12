# ABOUTME: Counts how many fitted names change if each independent is named after its own edge
# ABOUTME: (the member the solve keeps) instead of its kinetic group's naming representative.
const SCR = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(SCR, "proto", "common", "Population.jl"))
include(joinpath(SCR, "slice", "src", "Slice.jl"))
using Serialization, Printf
const SL = Slice

"""The names edge `k` renders for itself in a state, as `constant_names` renders a group."""
function own_names(m, k, fwd, inactive)
    e = SL.oriented(m, k, fwd)
    a, b = SL.side_name(m, e.a, e.cons), SL.side_name(m, e.b, e.rel)
    t = SL.state_tag(m, m.group[k], inactive)
    e.re ? [Symbol("K_", t, SL.is_binding(e) ? "$(b)_to_$a" : "$(a)_to_$b")] :
           [Symbol("k_", t, a, "_to_", b), Symbol("k_", t, b, "_to_", a)]
end

data = joinpath(SCR, "slice", "core", "data")
sets = [String(split(l, '\t')[1]) for l in readlines(joinpath(data, "MANIFEST.tsv"))[2:end]]
tot = Dict{String, NTuple{4, Int}}()
for set in sets
    recs, changed, names_changed, nnames = 0, 0, 0, 0
    for c in deserialize(joinpath(data, set * ".jls"))
        m = c.mech
        sol = try SL.solve(m) catch; continue end
        fwd = SL.orient(m)
        today = SL.column_names(m, sol)[sol.indep]
        kept = [own_names(m, sol.cols[i].edge, fwd, sol.cols[i].inactive)[sol.cols[i].slot]
                for i in sol.indep]
        recs += 1
        d = count(i -> today[i] != kept[i], eachindex(today))
        changed += d > 0
        names_changed += d
        nnames += length(today)
    end
    tot[set] = (recs, changed, names_changed, nnames)
    @printf("%-18s records %5d, fitted sets that change %5d, names that change %5d of %6d\n",
            set, recs, changed, names_changed, nnames)
    flush(stdout)
end
r = sum(first, values(tot)); c = sum(x -> x[2], values(tot))
n = sum(x -> x[3], values(tot)); N = sum(x -> x[4], values(tot))
@printf("all: records %d, fitted sets that change %d (%.1f%%), names %d of %d (%.1f%%)\n",
        r, c, 100c / r, n, N, 100n / N)
