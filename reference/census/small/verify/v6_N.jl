# ABOUTME: Verifies class N (never saturates) by growth scans of the compiled package law along joint
# ABOUTME: and single-reactant rays, on samples, the numeric-probe hits and whole small runs.
include(joinpath(@__DIR__, "vlib.jl"))
const SAMPLES = deserialize(joinpath(VER, "samples.jls"))

rng = Xoshiro(616161)
global nscan = 0
global nbad = 0
worst = Tuple{Float64, String, Int}[]

function scan(lab, name, d, mechs, g; show = false)
    r = d.recs[g]
    ms = mstruct(mechs[g])
    pl = plaw(mechs[g])
    gr = growth(pl, ms, rng; ndraw = 3)
    s = max(gr[:fwd], gr[:rev])
    global nscan += 1
    push!(worst, (s, name, g))
    bad = s > 0.05
    bad && (global nbad += 1)
    (show || bad) && @printf("  [%s] %s #%d vmax %s/%s probe %s/%s → max log-slope fwd %.2e rev %.2e %s\n",
                             lab, name, g, r[:vmax_fwd], r[:vmax_rev], r[:unb_fwd], r[:unb_rev],
                             gr[:fwd], gr[:rev], bad ? "GROWS" : "")
end

for name in ["uni_uni_rxn", "uni_uni_with_reg", "uni_uni_allo", "uni_uni_allo_reg",
             "uni_uni_allo_2reg", "uni_uni_reg_and_inhibitor", "uni_bi_rxn", "bi_bi_rxn",
             "bi_bi_pp_rxn", "bi_bi_rxn.r2_100000", "bi_bi_pp_rxn.r2_1200"]
    d, mechs = loadrun(name)
    s = SAMPLES[name]
    fitted = [r[:gid] for r in d.recs if r[:status] === :ok && r[:fitted]]
    println("\n== $name")
    for g in s[:N_probe]
        scan("probe hit", name, d, mechs, g; show = true)
    end
    whole = length(fitted) <= 400 ? fitted :
            fitted[sort!(randperm(rng, length(fitted))[1:250])]
    t = @elapsed for g in whole
        scan(length(fitted) <= 400 ? "all fitted" : "250 random fitted", name, d, mechs, g)
    end
    @printf("  scanned %d (%s) in %.0f s\n", length(whole), length(fitted) <= 400 ? "every fitted" : "random sample", t)
    d = mechs = nothing
    GC.gc()
end
sort!(worst; rev = true)
println("\nScanned $nscan candidates; growing (max log-slope > 0.05): $nbad")
println("Largest log-slopes: ", join((@sprintf("%.2e (%s #%d)", w...) for w in worst[1:min(end, 8)]), ", "))
