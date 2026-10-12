# ABOUTME: Verifies class W: admission re-implemented on every fitted allosteric candidate against the
# ABOUTME: census verdicts, the :OnlyA-limit oracle and the face oracle on flagged and unflagged samples.
include(joinpath(@__DIR__, "vlib.jl"))
include(joinpath(@__DIR__, "moracle.jl"))

const SAMPLES = deserialize(joinpath(VER, "samples.jls"))
census_kind(msg) = msg === nothing ? "none" :
    occursin("has none", msg) ? "nofree" : occursin("flips at one free form", msg) ? "secondfree" :
    occursin("still join the two sides", msg) ? "G1_inside" : occursin("close a cycle", msg) ? "G1_cycle" :
    occursin("reach the free enzyme only", msg) ? "island" : "other"
const ORDER = ["nofree", "secondfree", "G1_inside", "G1_cycle", "island"]
first_kind(v) = (i = findfirst(in(v), ORDER); i === nothing ? "none" : ORDER[i])

wflag(r) = get(r, :v0, false) || r[:re_turnover] !== nothing || any(!=(:none), r[:face_kind]) ||
           any(r[:dvan]) || r[:admit] !== nothing

function admission_all(name)
    d, mechs = loadrun(name)
    tab = Dict{Tuple{String, String}, Int}()
    multi = Dict{String, Int}()
    ex = Dict{Tuple{String, String}, Int}()
    nallo = 0
    t = @elapsed for r in d.recs
        r[:status] === :ok && r[:fitted] && r[:allo] || continue
        nallo += 1
        v = my_admit(mstruct(mechs[r[:gid]]))
        k = (census_kind(r[:admit]), first_kind(v))
        tab[k] = get(tab, k, 0) + 1
        k[1] != k[2] && (ex[k] = r[:gid])
        length(v) > 1 && (s = join(sort(v), "+"); multi[s] = get(multi, s, 0) + 1)
    end
    @printf("\n[%s] admission on %d fitted allosteric candidates (%.1f s): (census, mine) %s\n", name,
            nallo, t, join(("$(k[1])/$(k[2]) $v" for (k, v) in sort!(collect(tab))), ", "))
    isempty(multi) || println("   candidates with several violations: ",
                              join(("$k $v" for (k, v) in sort!(collect(multi))), ", "))
    for (k, g) in ex
        println("   disagreement $(k) #$g: ", d.recs[g][:text])
    end
    d, mechs
end

fmtp(xs) = join((@sprintf("%.5g", x) for x in xs), "/")

function onlya_line(lab, name, d, mechs, g, rng; ndraws = 2)
    r = d.recs[g]
    ms = mstruct(mechs[g])
    res = [onlya_test(ms, rng) for _ in 1:ndraws]
    fam = all(x -> x.family, res)
    spread = fam ? maximum(p.spread for x in res for p in x.points) : NaN
    vsp = fam ? maximum(p.vs_pkg for x in res for p in x.points) : NaN
    pspread = fam ? maximum(x.points[1].spread for x in res) : NaN
    @printf("  [%s] %s #%d census %s mine %s | family %s nfree %s | spread pos %.1e all %.1e | vs package model %.1e\n",
            lab, name, g, census_kind(r[:admit]), join(my_admit(ms), "+"), fam,
            join((x.nfree for x in res), "/"), pspread, spread, vsp)
    (fam = fam, spread = spread, vsp = vsp)
end

function face_line(lab, name, d, mechs, g, rng; ndraws = 2, lawcheck = true)
    r = d.recs[g]
    m = mechs[g]
    ms = mstruct(m)
    fts = [face_test(ms, rng) for _ in 1:ndraws]
    bad = unique([f.face for ft in fts for f in ft.faces if f.fails])
    worst_pat = maximum(f.pattern for ft in fts for f in ft.faces)
    worst_gap = maximum(f.gap for ft in fts for f in ft.faces)
    census = [string(r[:face_kind][i], "@", r[:faces][i]) for i in eachindex(r[:faces]) if r[:face_kind][i] !== :none]
    any(r[:dvan]) && push!(census, "dvan")
    lv = ""
    if lawcheck && !ms.allo
        pl = plaw(m)
        lo = law_vs_oracle(m, ms, pl, rng)
        pos = lo[1]
        fb = [p.point for p in lo[2:end] if p.law_vs_exact > 1e-8 || p.exact_spread > 1e-8]
        lv = @sprintf(" | law vs oracle: positive %.1e; faces law≠exact or spread: %s (law vs interior max %.1e)",
                      max(pos.law_vs_interior, pos.law_vs_exact), isempty(fb) ? "none" : join(fb, ", "),
                      maximum(p.law_vs_interior for p in lo))
    end
    @printf("  [%s] %s #%d r%d %s %s admit %s census %s | oracle fails %s (max pattern %.1e, max gap %.1e)%s\n",
            lab, name, g, r[:round], r[:move], ms.allo ? "allo" : "plain", census_kind(r[:admit]),
            isempty(census) ? "-" : join(census, ","),
            isempty(bad) ? "none" : join(bad, ", "), worst_pat, worst_gap, lv)
    (bad = !isempty(bad), census = !isempty(census))
end

rng = Xoshiro(20261011)
tally = Dict{String, Vector{Int}}()     # label => [n, flagged-by-oracle]
bump!(lab, hit) = (v = get!(tally, lab, [0, 0]); v[1] += 1; v[2] += hit)

for name in ["uni_uni_rxn", "uni_uni_with_reg", "uni_bi_rxn", "bi_bi_rxn", "bi_bi_pp_rxn",
             "bi_bi_rxn.r2_100000", "bi_bi_pp_rxn.r2_1200"]
    d, mechs = admission_all(name)
    s = SAMPLES[name]
    println("  -- :OnlyA-limit oracle")
    for g in s[:W_island]
        x = onlya_line("island", name, d, mechs, g, rng)
        bump!("onlya island: family ok and spread>1e-6", x.fam && x.spread > 1e-6)
        y = face_line("island-face", name, d, mechs, g, rng)
        bump!("island: package model fails a face (oracle)", y.bad)
    end
    for g in s[:W_G1]
        x = onlya_line("G1", name, d, mechs, g, rng)
        bump!("onlya G1: no family", !x.fam)
        y = face_line("G1-face", name, d, mechs, g, rng)
        bump!("G1: package model fails a face (oracle)", y.bad)
        bump!("G1: census face flag", y.census)
        bump!("G1: oracle face verdict = census face verdict", y.bad == y.census)
    end
    for g in s[:W_unflagged_allo]
        x = onlya_line("admitted", name, d, mechs, g, rng)
        bump!("onlya admitted: family ok, spread<1e-12, = package model", x.fam && x.spread < 1e-12 && x.vsp < 1e-12)
    end
    println("  -- face oracle")
    for g in s[:W_face]
        x = face_line("face-flagged", name, d, mechs, g, rng)
        bump!("face flagged: oracle fails", x.bad)
    end
    for g in s[:W_unflagged_plain]
        x = face_line("unflagged-plain", name, d, mechs, g, rng)
        bump!("face unflagged plain: oracle fails", x.bad)
    end
    for g in s[:W_unflagged_allo][1:min(end, 6)]
        x = face_line("unflagged-allo", name, d, mechs, g, rng)
        bump!("face unflagged allo: oracle fails", x.bad)
    end
    d = mechs = nothing
    GC.gc()
end
println("\nTally (n, hits):")
for (k, v) in sort!(collect(tally))
    println("  $k: $(v[2]) of $(v[1])")
end
