# ABOUTME: Tabulates the classes of fitted expansion children by the move that made them (flip,
# ABOUTME: split, to-allosteric, allosteric-state change), from saved runs.

include(joinpath(@__DIR__, "census_lib.jl"))
include(joinpath(@__DIR__, "summarize.jl"))

"""The expansion move that turned `p` into `c`, read off their structure."""
function move_of(p, c)
    pa, ca = p isa ER.AllostericMechanism, c isa ER.AllostericMechanism
    !pa && ca && return "to-allosteric"
    ng(m) = length(ER.steps(m))
    nre(m) = count(ER.is_equilibrium, (s for g in ER.steps(m) for s in g))
    ng(c) > ng(p) && return "split"
    nre(c) < nre(p) && return "flip RE->SS"
    pa && ca && ER.cat_allo_states(c) != ER.cat_allo_states(p) && return "allo state change"
    "other"
end

for name in Symbol.(ARGS)
    data = deserialize(joinpath(@__DIR__, "..", "data", "$(name).jls"))
    rounds = data.rounds
    label!(rounds)
    byid = Dict(r[:id] => r for recs in rounds for r in recs)
    println("="^100, "\n", name)
    for (ri, recs) in enumerate(rounds)
        ri == 1 && continue
        tab = Dict{String, Dict{String, Int}}()
        for r in recs
            r[:status] === :ok || continue
            mv = move_of(byid[r[:parent]][:mech], r[:mech])
            t = get!(tab, mv, Dict{String, Int}())
            t["compiled"] = get(t, "compiled", 0) + 1
            r[:fitted] || continue
            t["fitted"] = get(t, "fitted", 0) + 1
            ls = r[:labels]
            for (k, hit) in (("W", any(isW, ls)), ("D", :D in ls), ("P", :P in ls),
                             ("Ptoday", :Ptoday in ls), ("N", :N in ls), ("Z", :Z in ls),
                             ("Xlaw", :Xlaw in ls),
                             ("OK", isempty(setdiff(ls, [:Xlaw]))),
                             ("useless(W|D|P|Z)", any(l -> isW(l) || l in (:D, :P, :Z), ls)))
                hit && (t[k] = get(t, k, 0) + 1)
            end
        end
        println("round $(ri - 1): fitted children by move (shares of the move's fitted)")
        for (mv, t) in sort!(collect(tab); by = x -> -get(x[2], "fitted", 0))
            f = get(t, "fitted", 0)
            parts = join(("$k $(get(t, k, 0)) ($(pct(get(t, k, 0), f)))" for k in
                          ("W", "D", "P", "Ptoday", "Z", "N", "Xlaw", "useless(W|D|P|Z)", "OK")), ", ")
            @printf("  %-18s compiled %5d fitted %5d: %s\n", mv, get(t, "compiled", 0), f, parts)
        end
    end
    flush(stdout)
end
