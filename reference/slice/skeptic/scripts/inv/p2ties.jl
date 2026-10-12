# ABOUTME: The slice's tied group members of P2's adv/9 and adv/10 (isomerization groups), and
# ABOUTME: the independent check with chemistry-aware alignment.
include(joinpath(@__DIR__, "..", "spotlib.jl"))
recs = skeptic_records("p2_adv_pkg")
for id in ("adv/9", "adv/10")
    r = only(x for x in recs if x.p.id == id)
    m = r.c.mech
    fwd = S.orient(m)
    println("== ", id, " ", r.p.source)
    for g in 1:S.ngroups(m)
        ks = S.members(m, g)
        length(ks) > 1 || continue
        println("  group $g: ", join((S.step_text(m, S.oriented(m, k, fwd)) for k in ks), "  ||  "))
    end
    x = check(r, S.derive(m), Xoshiro(1))
    @printf("  check: count %d/%d, law vs oracle %.1e, dependents %.1e, kcat %s\n", x.mine,
            x.slice, x.errB, x.worst_dep, x.kc)
end
