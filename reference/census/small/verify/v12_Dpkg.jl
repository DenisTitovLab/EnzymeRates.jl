# ABOUTME: Examines D pairs whose package laws did not match: equations side by side, a harder family
# ABOUTME: fit, and the thermodynamically correct law (oracle) of Z children fitted by the parent's law.
include(joinpath(@__DIR__, "vlib.jl"))
include(joinpath(@__DIR__, "moracle.jl"))
const SAMPLES = deserialize(joinpath(VER, "samples.jls"))
rng = Xoshiro(121212)

const ODD = [("bi_bi_rxn", 3642, 3637), ("bi_bi_pp_rxn", 3671, 3663),
             ("bi_bi_rxn.r2_100000", 10285, 10277), ("bi_bi_rxn.r2_100000", 10040, 10004),
             ("bi_bi_rxn.r2_100000", 4129, 4121), ("bi_bi_pp_rxn.r2_1200", 13267, 13251)]

function eqlines(m)
    s = ER.rate_equation_string(ER.compile_mechanism(m))
    filter(l -> !startswith(strip(l), "(; S") && !startswith(strip(l), "(; A"), split(s, '\n'))
end

cache = Dict{String, Any}()
getrun(name) = get!(() -> loadrun(name), cache, name)
for (name, x, y) in ODD
    d, mechs = getrun(name)
    println("\n== $name #$x vs #$y")
    println("  #$x ", d.recs[x][:text]); println("  #$y ", d.recs[y][:text])
    a, b = plaw(mechs[x]), plaw(mechs[y])
    println("  names #$x: ", join(a.names, ", ")); println("  names #$y: ", join(b.names, ", "))
    for l in eqlines(mechs[x]); println("   x| ", first(l, 600)); end
    for l in eqlines(mechs[y]); println("   y| ", first(l, 600)); end
    keq = lu10(rng, -1, 1)
    best_ab = minimum(fit_into(a, [lu10(rng, -1, 1) for _ in a.names], b, keq, rng; starts = 40, iters = 800) for _ in 1:3)
    best_ba = minimum(fit_into(b, [lu10(rng, -1, 1) for _ in b.names], a, keq, rng; starts = 40, iters = 800) for _ in 1:3)
    worst_ab = maximum(fit_into(a, [lu10(rng, -1, 1) for _ in a.names], b, keq, rng; starts = 40, iters = 800) for _ in 1:3)
    worst_ba = maximum(fit_into(b, [lu10(rng, -1, 1) for _ in b.names], a, keq, rng; starts = 40, iters = 800) for _ in 1:3)
    @printf("  harder fit (40 starts, 3 draws each way): x-in-y best %.1e worst %.1e; y-in-x best %.1e worst %.1e\n",
            best_ab, worst_ab, best_ba, worst_ba)
    flush(stdout)
end

# the thermodynamically correct law of D/Z pairs whose package counts differ: the oracle of the
# later candidate (its own kinetic groups, consistent constants) fitted by the earlier one's law
println("\n== correct-law test: oracle of the later candidate fitted by the earlier candidate's package law")
global nyes = 0
global ntot = 0
for name in ["uni_bi_rxn", "bi_bi_rxn", "bi_bi_pp_rxn", "bi_bi_rxn.r2_100000", "bi_bi_pp_rxn.r2_1200"]
    d, mechs = getrun(name)
    prs = SAMPLES[name][:D_pairs_npdiff]
    pick = prs[sort!(randperm(rng, length(prs))[1:min(8, length(prs))])]
    for (x, y) in pick
        msx = mstruct(mechs[x])
        msx.allo && msx.n == 1 || continue
        px, py = plaw(mechs[x]), plaw(mechs[y])
        A = aconsts(msx, rng)
        Ip = iconsts(msx, A, rng)
        lnkp = Ip.lnk(Ip.draw(rng))
        L = big(lu10(rng, -1, 1))
        σ = [big(lu10(rng, -1, 1)) for _ in 1:(length(msx.steps) + 1)]
        keq = Float64(exp(BigFloat(sum(A.mu[z] for z in msx.subs) - sum(A.mu[z] for z in msx.prods))))
        npts = 2 * max(length(px.names), length(py.names)) + 14
        pts = [Dict{Symbol, Float64}(z => lu10(rng, -1, 1) for z in msx.mets) for _ in 1:npts]
        vx = [Float64(oracle_rate(msx, A, Ip, lnkp, L, Dict{Symbol, BigFloat}(k => big(v) for (k, v) in c),
                                  big(10)^30, σ)) for c in pts]
        ry = fit_values(py, keq, pts, vx, rng; starts = 20)
        rx = fit_values(px, keq, pts, vx, rng; starts = 20)
        global ntot += 1
        ry < 1e-8 && (global nyes += 1)
        @printf("  %s #%d (np %d) earlier #%d (np %d): oracle of #%d fitted by its own law %.1e, by #%d's law %.1e\n",
                name, x, length(px.names), y, length(py.names), x, rx, y, ry)
        flush(stdout)
    end
end
println("  later candidate's correct law inside the earlier one's package family: $nyes of $ntot")
