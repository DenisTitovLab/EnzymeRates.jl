# ABOUTME: Allosteric check (c): the graph form of G1 against an exact rational LP (phase-1 simplex,
# ABOUTME: Bland's rule) on the joint system's rows, and both against today's guard on hand cases.
include(joinpath(@__DIR__, "lib.jl"))

include(joinpath(@__DIR__, "g1_lp.jl"))

records = [r for r in load_records() if S.is_allosteric(r.c.mech)]
agree, n, lp_t, mismatch = 0, 0, Float64[], String[]
hand = Dict{Any, Int}()
for r in records
    m, p = r.c.mech, r.p
    t = @elapsed lp = g1_lp(m)
    push!(lp_t, t)
    gr = g1_graph(m)
    global n += 1
    lp == gr ? (global agree += 1) : push!(mismatch, "$(r.set) $(p.id): LP $lp graph $gr")
    if !(r.set in ("specs", "handwritten", "allo_uniuni", "allo_bibi", "expand_bibi",
                   "allo_uniuni_inh", "allo_unibi", "allo_bibi_pp", "allo_terbi"))
        k = (startswith(r.set, "p5a_cases") ? (startswith(p.source, "uni:") ? "grid" :
             "P5a special") : "AD", gr ? "G1 feasible" : "G1 empty",
             p.ok ? "guard admits" : "guard rejects")
        hand[k] = get(hand, k, 0) + 1
        startswith(r.set, "p5a_adv") && (gr != p.ok) &&
            println("  G1 and today's guard differ: $(p.id) (G1 $(gr ? "feasible" : "empty"), guard $(p.ok ? "admits" : "rejects"))")
    end
end
@printf("G1 graph test = exact LP in %d of %d allosteric records (LP median %.2f ms, max %.0f ms)\n",
        agree, n, 1e3 * median(lp_t), 1e3 * maximum(lp_t))
foreach(println, mismatch)
println("hand-built cases (kind, G1, today's guard) => count:")
foreach(((k, v),) -> println("  ", k, " => ", v), sort!(collect(hand); by = string))
