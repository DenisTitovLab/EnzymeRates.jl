# ABOUTME: Prints the steps and untied-variant cone constraints of chosen mechanisms (P5a).
# ABOUTME: Usage: julia peek_untie.jl <population> <index...>
include(joinpath(@__DIR__, "MWC.jl"))
pop = Pop.load_population(ARGS[1])
for i in parse.(Int, ARGS[2:end])
    m = pop[i]
    println("== ", m.id, " family=", m.family)
    gos = Pop.group_of_step(m); ch = MWC.is_chem(m; per_group = false)
    for (j, s) in enumerate(m.steps)
        println("   g", gos[j], " [", m.group_tags[gos[j]], "] ", ch[j] ? "chem " : "bind ", s.re ? "RE " : "SS ", Pop.step_text(m, s))
    end
    J = MWC.joint_solve(m; untie = true)
    sig = MWC.cone_signs(J, :i)
    O = sort!(collect(keys(sig)))
    println("  moved: ", [(J.cols[i].name, sig[i]) for i in O])
    for r in unique([J.C[r, O] for r in axes(J.C, 1) if any(!iszero, J.C[r, O])]); println("     row ", r); end
end
