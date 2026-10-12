# ABOUTME: Prints the plain records the solve rejects, with their steps, to judge whether the
# ABOUTME: rejection is right (a cycle whose net change is no multiple of the net reaction).
include(joinpath(@__DIR__, "lib.jl"))
for r in load_records()
    m = r.c.mech
    r.p.ok || println("not ok: ", r.set, " ", r.p.id, " ", r.p.source, ": ", first(r.p.error, 200))
    try
        S.solve(m)
    catch err
        println("== ", r.set, " ", r.p.id, " ", r.p.source, ": ", sprint(showerror, err))
        println("   substrates ", m.substrates, " products ", m.products)
        for e in m.edges; println("   ", S.step_text(m, e)); end
    end
end
