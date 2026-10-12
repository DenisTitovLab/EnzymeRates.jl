# ABOUTME: Prints an allosteric mechanism's steps, tags, moved columns and cone constraints (P5a).
# ABOUTME: Usage: julia peek_cone.jl <population> <index...>
include(joinpath(@__DIR__, "MWC.jl"))
pop = Pop.load_population(ARGS[1])
for i in parse.(Int, ARGS[2:end])
    m = pop[i]
    println("== ", m.id, " ", m.source)
    gos = Pop.group_of_step(m); ch = MWC.is_chem(m)
    for (j, s) in enumerate(m.steps)
        println("   g", gos[j], " [", m.group_tags[gos[j]], "] ", ch[j] ? "chem " : "bind ", Pop.step_text(m, s))
    end
    J = MWC.joint_solve(m)
    for conv in (:i, :ii)
        sig = MWC.cone_signs(J, conv)
        O = sort!(collect(keys(sig)))
        println("  conv $conv moved columns: ", [(J.cols[i].name, sig[i]) for i in O])
        rows = [J.C[r, O] for r in axes(J.C, 1) if any(!iszero, J.C[r, O])]
        for r in unique(rows); println("     row ", r); end
        dirs, d = MWC.cone_directions(J, conv)
        println("     d=$d feasible=$(length(dirs))")
    end
end
