# ABOUTME: Fit-cost calibration, step 6: how the fit cost scales with the number of data points (48 vs
# ABOUTME: 480) for 12 candidates, and a direct saturation check of the candidate the numeric probe flagged.
using EnzymeRates, Serialization, Random, Printf, Statistics, LinearAlgebra
const ER = EnzymeRates
const FCD = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/fitcost"
include(joinpath(FCD, "fitlib.jl"))
const KEQ = 10.0

sels = [deserialize(joinpath(FCD, "data", f)) for f in ("selection.jls", "selection_big.jls")]
mech = Dict{Tuple{Symbol, Int}, Any}()
rows = Dict{Tuple{Symbol, Int}, Any}()
for sk in sels, s in sk, (r, m) in s.sel
    mech[(s.name, r.gid)] = m
    rows[(s.name, r.gid)] = r
end

# ── 1. saturation of bi_bi #1858 along rays (substrates up, products at zero) ──────────────
let em = ER.compile_mechanism(mech[(:bi_bi_rxn, 1858)]), pn = ER.fitted_params(em)
    println("bi_bi #1858: v as substrates grow (products 0), BigFloat (Float64 in brackets):")
    rng = Xoshiro(11)
    for trial in 1:4
        θ = [10.0^(2rand(rng) - 1) for _ in pn]
        for (label, conc) in (("A = B = s", s -> (A = s, B = s, P = zero(s), Q = zero(s))),
                              ("A = s, B = 1", s -> (A = s, B = one(s), P = zero(s), Q = zero(s))),
                              ("B = s, A = 1", s -> (A = one(s), B = s, P = zero(s), Q = zero(s))))
            vb = setprecision(BigFloat, 256) do
                pb = merge(NamedTuple{pn}(Tuple(big.(θ))), (Keq = big(KEQ), E_total = big(1.0)))
                [Float64(ER.rate_equation(em, conc(big(10.0)^e), pb)) for e in 0:3:12]
            end
            pf = merge(NamedTuple{pn}(Tuple(θ)), (Keq = KEQ, E_total = 1.0))
            vf = [ER.rate_equation(em, conc(10.0^e), pf) for e in 0:3:12]
            @printf("  θ draw %d, %-13s s = 1e0,1e3,1e6,1e9,1e12: %s  [%s]\n", trial, label,
                    join((@sprintf("%.6g", v) for v in vb), ", "),
                    join((@sprintf("%.6g", v) for v in vf), ", "))
        end
    end
end
flush(stdout)

# ── 2. data size: the same design repeated 10 times (480 points, fresh noise) ─────────────
"""`make_data` with every design row repeated `nrep` times (independent noise per copy)."""
function make_data_rep(name, em, θ, keq, rng, nrep; noise = 0.05, gscale = 0.3)
    mets = ER.metabolites(em)
    rows = design(name)
    scale = Dict(g => 10.0^(gscale * randn(rng)) for g in unique(first.(rows)))
    params = merge(θ, (Keq = keq, E_total = 1.0))
    cols = Dict(x => Float64[] for x in mets)
    G, R = String[], Float64[]
    for _ in 1:nrep, (g, c) in rows
        concs = NamedTuple{mets}(Tuple(get(c, x, 0.0) for x in mets))
        v = ER.rate_equation(em, concs, params)
        push!(G, g)
        push!(R, v * scale[g] * exp(noise * randn(rng)))
        for x in mets
            push!(cols[x], get(c, x, 0.0))
        end
    end
    merge((group = G, Rate = R), NamedTuple{mets}(Tuple(cols[x] for x in mets)))
end

const PICKS = [(:bi_bi_rxn, 1365), (:bi_bi_rxn, 2068), (:bi_bi_rxn, 2573), (:bi_bi_rxn, 1864),
               (:bi_bi_rxn, 2823), (:bi_bi_rxn, 2790), (:ter_bi_rxn, 5208), (:ter_bi_rxn, 1777),
               (:ter_bi_rxn, 5416), (:ter_bi_rxn, 5932), (:ter_bi_rxn, 5447), (:ter_bi_rxn, 5621)]
opt = CMAEvolutionStrategyOpt()
datas = Dict{Tuple{Symbol, Int}, Any}()
for s in sels[1]
    (rrow, rm) = s.ref
    em = ER.compile_mechanism(rm)
    pn = ER.fitted_params(em)
    rng = Xoshiro(s.name === :bi_bi_rxn ? 101 : 202)
    θ = NamedTuple{pn}(Tuple(10.0 .^ (rand(rng, length(pn)) .- 0.5)))
    datas[(s.name, 48)] = make_data(s.name, em, θ, KEQ, rng)[1]     # the main runs' data
    datas[(s.name, 480)] = make_data_rep(s.name, em, θ, KEQ, Xoshiro(303), 10)
end
println("\ndata size: rxn gid np ops | points: t_loss us, generations/20, t_fit s (warm), us/gen, best loss")
open(joinpath(FCD, "data", "datasize.tsv"), "w") do io
    println(io, "rxn\tgid\tnp\tops\tnpts\tt_loss\tgens\tt_fit\tbest_loss")
    first_call = true
    for (name, gid) in PICKS
        m = mech[(name, gid)]
        em = ER.compile_mechanism(m)
        np = length(ER.fitted_params(em))
        ops = rows[(name, gid)].pkg_ops
        seed = gid + (name === :bi_bi_rxn ? 0 : 1_000_000)
        for npts in (48, 480)
            fp = ER.FittingProblem(em, datas[(name, npts)]; Keq = KEQ, scale_k_to_kcat = 1.0)
            ER.loss!(zeros(np), fp)                       # compile before timing
            if first_call                                  # compile the restart loop once
                fit_instrumented(fp, opt; n_restarts = 1)
                first_call = false
            end
            tl, _, _ = time_loss(fp, np, Xoshiro(seed))
            Random.seed!(seed)
            t = time_ns()
            f = fit_instrumented(fp, opt; n_restarts = 20, maxtime = 600.0)
            tf = (time_ns() - t) / 1e9
            g = sum(r.gens for r in f.restarts)
            @printf("  %-10s #%-5d np %2d ops %5d | %3d: %7.2f us, %6.0f, %6.2f s, %6.1f us/gen, %.4g\n",
                    name, gid, np, ops, npts, tl * 1e6, g / 20, tf, tf / g * 1e6, f.loss)
            println(io, join((name, gid, np, ops, npts, tl, g, tf, f.loss), '\t'))
            flush(stdout)
        end
    end
end
