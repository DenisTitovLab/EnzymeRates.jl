# ABOUTME: Prints split-difference examples: the parent, the unit, and both counts' independents
# ABOUTME: on the slice side; exports the same mechanisms for the package's fitted lists.
include(joinpath(@__DIR__, "load.jl"))
include(joinpath(@__DIR__, "diag_lib.jl"))
rows = deserialize(joinpath(S, "slice/enum/data/cmp_expand_bibi.jls"))
out = Pop.PMech[]
for (index, set) in ((24, [1]), (163, [2]), (32, [1]))
    r = only(filter(r -> r.index == index, rows))
    m = r.mech
    units, child = _testhelper_units(m)
    c = child(set)
    fwd = Slice.orient(m)
    println("== expand_bibi/$index, unit set $set: group $(units[set[1]][1]) split into ",
            units[set[1]][2], " | ", units[set[1]][3], " reverted ", units[set[1]][4])
    for (k, e) in enumerate(m.edges)
        println("   e$k g$(m.group[k]) $(m.tags[m.group[k]])  ", Slice.step_text(m, Slice.oriented(m, k, fwd)))
    end
    println("   inactive mask: ", findall(Slice.inactive_mask(m)))
    for (lbl, x) in (("parent", m), ("child", c))
        sol = Slice.solve(x)
        names = [(cl.inactive ? Slice.constant_names(x; inactive = true) : Slice.constant_names(x))[cl.edge][cl.slot] for cl in sol.cols]
        println("   $lbl slice indep ($(length(sol.indep))): ", join(names[sol.indep], " "))
    end
    push!(out, _testhelper_to_pmech(m, "expand_bibi/$index"), _testhelper_to_pmech(c, "expand_bibi/$index/" * join(set, ",")))
end
serialize(joinpath(S, "slice/enum/data/diag_examples.jls"), out)
