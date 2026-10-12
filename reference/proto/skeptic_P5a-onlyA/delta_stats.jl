# ABOUTME: Size of the departure from deletion (P5a): for the mechanisms whose limit is not
# ABOUTME: unique or not deletion, the spread of the maximum relative difference to deletion.
using Serialization, Printf, Statistics
for name in ("expand_bibi", "allo_bibi")
    recs = deserialize(joinpath(@__DIR__, "data", "pop_$name.jls"))
    for conv in (:i, :ii, :iu)
        rs = filter(r -> r[:conv] === conv && haskey(r, :cone_v), recs)
        bad = filter(r -> !(r[:cone_spread] < 1e-20 && r[:cone_vs_del] < 1e-20), rs)
        isempty(bad) && (@printf("%-12s %-3s none\n", name, conv); continue)
        d = sort([r[:cone_vs_del] for r in bad])
        q(p) = d[clamp(round(Int, p * length(d)), 1, length(d))]
        @printf("%-12s %-3s n=%3d  max rel diff to deletion over directions: min %.1e  median %.1e  90%% %.1e  max %.1e;  ≥1e-2: %d\n",
                name, conv, length(d), d[1], q(0.5), q(0.9), d[end], count(>=(1e-2), d))
    end
end
