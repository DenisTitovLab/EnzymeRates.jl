# ABOUTME: Re-tests faces with the long-time state from the free enzyme: every G1-rejected candidate,
# ABOUTME: the face-flagged ones and unflagged samples; examines the two odd :OnlyA-oracle cases.
include(joinpath(@__DIR__, "vlib.jl"))
include(joinpath(@__DIR__, "moracle.jl"))
const SAMPLES = deserialize(joinpath(VER, "samples.jls"))
census_kind(msg) = msg === nothing ? "none" : occursin("still join", msg) ? "G1_inside" :
                   occursin("close a cycle", msg) ? "G1_cycle" : occursin("reach the free", msg) ? "island" : "other"
census_face(r) = [string(r[:face_kind][i], "@", r[:faces][i]) for i in eachindex(r[:faces]) if r[:face_kind][i] !== :none]

rng = Xoshiro(111111)
tab = Dict{Tuple{String, Bool, Bool}, Int}()
function ftest(lab, d, mechs, g; draws = 2, show = true)
    r = d.recs[g]
    ms = mstruct(mechs[g])
    res = [face_test_free(ms, rng) for _ in 1:draws]
    bad = unique([f.face for ft in res for f in ft if f.fails])
    worst = maximum(max(f.gap, f.pattern) for ft in res for f in ft)
    nr = maximum(f.nreach for ft in res for f in ft)
    nt = maximum(f.ntotal for ft in res for f in ft)
    cf = census_face(r)
    k = (lab, !isempty(cf), !isempty(bad))
    tab[k] = get(tab, k, 0) + 1
    show && @printf("  [%s] #%d %s admit %s | census %s | from-free oracle fails %s (worst %.1e; terminal classes reachable ≤%d, in all ≤%d)\n      %s\n",
                    lab, g, r[:move], census_kind(r[:admit]), isempty(cf) ? "-" : join(cf, ","),
                    isempty(bad) ? "none" : join(bad, ", "), worst, nr, nt, r[:text])
    flush(stdout)
end

for name in ["bi_bi_rxn.r2_100000", "bi_bi_pp_rxn.r2_1200"]
    d, mechs = loadrun(name)
    println("\n== $name")
    g1 = [r[:gid] for r in d.recs if r[:status] === :ok && r[:fitted] &&
          census_kind(r[:admit]) in ("G1_inside", "G1_cycle")]
    for g in g1
        ftest(census_kind(d.recs[g][:admit]), d, mechs, g)
    end
    for g in SAMPLES[name][:W_face]
        ftest("face-flagged", d, mechs, g; draws = 1, show = false)
    end
    for g in [SAMPLES[name][:W_unflagged_plain]; SAMPLES[name][:W_unflagged_allo]]
        ftest("unflagged", d, mechs, g; draws = 1, show = false)
    end
    if name == "bi_bi_rxn.r2_100000"
        # the admitted candidate for which no :OnlyA family was found
        g = 19519
        ms = mstruct(mechs[g])
        println("\n  #19519 ", d.recs[g][:text], "\n   my_admit ", my_admit(ms), "; census admit ", d.recs[g][:admit])
        A = aconsts(ms, rng)
        Ip = iconsts(ms, A, rng)
        pin = Ip.pins(Ip.draw(rng))
        nb = count(k -> !Ip.kept[k] && plain_binding(ms, ms.steps[k]), eachindex(ms.steps))
        found = 0
        for t in 1:200
            expo = t == 1 ? nothing : rand(rng, 1:8, nb)
            If = iconsts(ms, A, rng; family = true, expo, pin)
            If === nothing || (found = t; break)
        end
        println("   :OnlyA plain bindings $nb; family found at try $found (0 = none in 200, exponents 1:8)")
        If0 = iconsts(ms, A, rng; family = true, expo = rand(rng, 1:8, nb))
        If1 = iconsts(ms, A, rng; family = true)
        println("   without the pins to the package model: family exists for exponents all 1: $(If1 !== nothing)")
        for (k, s) in enumerate(ms.steps)
            println("     step $k group $(s.g) tag $(ms.tags[s.g]) re $(s.re) plain-binding $(plain_binding(ms, s)) kept $(Ip.kept[k]) inI $(Ip.present[k]): ",
                    ER.name(ms.forms[s.a]), " + ", join(s.cons, "+"), " → ", ER.name(ms.forms[s.b]), " + ", join(s.rel, "+"))
        end
    end
    d = mechs = nothing
    GC.gc()
end

# the island whose two family draws agreed
d, mechs = loadrun("bi_bi_pp_rxn")
g = 3273
ms = mstruct(mechs[g])
println("\n  #3273 ", d.recs[g][:text])
for k in 1:3
    ot = onlya_test(ms, rng)
    println("   draw $k: nfree $(ot.nfree); positive pkg $(ot.points[1].pkg) family $(ot.points[1].fam)")
end
for (k, s) in enumerate(ms.steps)
    println("     step $k group $(s.g) tag $(ms.tags[s.g]) re $(s.re) plain-binding $(plain_binding(ms, s)): ",
            ER.name(ms.forms[s.a]), " + ", join(s.cons, "+"), " → ", ER.name(ms.forms[s.b]), " + ", join(s.rel, "+"))
end

println("\n(label, census face flag, from-free oracle fails) counts:")
for (k, v) in sort!(collect(tab))
    println("  $k: $v")
end
