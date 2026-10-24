# ABOUTME: Fit cost of ter-ter candidates as identify fits them (CMA-ES restarts, relative loss, 48
# ABOUTME: synthetic points): per-restart time, evaluations and per-loss cost. Args: <first> <last> <cap s>.
using EnzymeRates, OptimizationCMAEvolutionStrategy, Random, Printf, Statistics, Serialization
const ER = EnzymeRates
const Optimization = ER.Optimization
println("EnzymeRates ", pkgversion(EnzymeRates), "; Julia ", VERSION)

const SERIES = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
const KEQ = 1000.0

"""48 points in 4 groups over the ter-ter metabolites A, B, D (substrates) and P, Q, R
(products), the design of census/fitcost/fitlib.jl:9-45 extended to three substrates."""
function design()
    rows = Tuple{String, Dict{Symbol, Float64}}[]
    add!(g; kw...) = push!(rows, (g, Dict{Symbol, Float64}(kw)))
    for x in (0.3, 3.0), a in SERIES
        add!("G1"; A = a, B = x, D = x)
    end
    for b in SERIES
        add!("G2"; A = 1.0, B = b, D = 1.0)
    end
    for d in SERIES
        add!("G2"; A = 1.0, B = 1.0, D = d)
    end
    for a in SERIES
        add!("G3"; A = a, B = 1.0, D = 1.0, P = 0.1)
    end
    for d in SERIES
        add!("G4"; A = 1.0, B = 1.0, D = d, R = 0.1)
    end
    rows
end

"""Relative-mode data from `em` at constants `θ` (fitcost's make_data: group factors
10^(0.3·N(0,1)), noise exp(0.05·N(0,1)))."""
function make_data(em, θ, rng)
    mets = ER.metabolites(em)
    rows = design()
    scale = Dict(g => 10.0^(0.3 * randn(rng)) for g in unique(first.(rows)))
    params = merge(θ, (Keq = KEQ, E_total = 1.0))
    cols = Dict(x => Float64[] for x in mets)
    G, R = String[], Float64[]
    for (g, c) in rows
        concs = NamedTuple{mets}(Tuple(get(c, x, 0.0) for x in mets))
        v = ER.rate_equation(em, concs, params)
        (isfinite(v) && v > 0) || error("reference rate $v")
        push!(G, g)
        push!(R, v * scale[g] * exp(0.05 * randn(rng)))
        foreach(x -> push!(cols[x], get(c, x, 0.0)), mets)
    end
    merge((group = G, Rate = R), NamedTuple{mets}(Tuple(cols[x] for x in mets)))
end

"""Median seconds per `loss!` call at 5 random points."""
function time_loss(fp, np, rng)
    xs = [clamp.(randn(rng, np) .* 2.0, -15.0, 15.0) for _ in 1:5]
    acc, ts = 0.0, Float64[]
    ER.loss!(xs[1], fp)
    for x in xs
        n = 200
        t = @elapsed for _ in 1:n
            acc += ER.loss!(x, fp)
        end
        push!(ts, t / n)
    end
    median(ts)
end

const FIRST, LAST, CAP = parse(Int, ARGS[1]), parse(Int, ARGS[2]), parse(Float64, ARGS[3])
picks = deserialize(joinpath(@__DIR__, "data", "fitpick.jls"))
opt = CMAEvolutionStrategyOpt()
for (k, (set, label, id, m)) in enumerate(picks)
    FIRST <= k <= LAST || continue
    rng = Xoshiro(k)
    t0 = time()
    em = ER.compile_mechanism(m)
    names = ER.fitted_params(em)
    np = length(names)
    θ = NamedTuple{names}(Tuple(10.0^(rand(rng) - 0.5) for _ in 1:np))
    data = try
        make_data(em, θ, rng)
    catch e
        @printf("fit %d %s %s: no data (%s)\n", k, set, id, sprint(showerror, e)[1:min(end, 80)])
        continue
    end
    fp = ER.FittingProblem(em, data; Keq = KEQ)
    tfixed = time() - t0
    tl = time_loss(fp, np, rng)
    lb, ub = fill(-15.0, np), fill(15.0, np)
    rs = NamedTuple[]
    total = 0.0
    Random.seed!(k)
    for r in 1:20
        x0 = clamp.(randn(np) .* 2.0, lb, ub)
        prob = ER._optimization_problem(fp, x0, lb, ub)
        t = @elapsed sol = Optimization.solve(prob, opt; maxtime = 600.0, maxiters = 10_000_000)
        gens = sol.stats.iterations
        lam = sol.original.p.λ
        push!(rs, (t = t, gens = gens, fevals = gens * lam, loss = sol.objective))
        total += t
        total > CAP && break
    end
    ts, fe = [r.t for r in rs], [r.fevals for r in rs]
    @printf("fit %d %s %s: np %d, fixed %.2f s (compile, data, problem), loss %.2f µs, %d restarts in %.1f s: median %.2f s, %.0f evaluations, %.0f generations; overhead per evaluation %.2f µs; best loss %.4g\n",
            k, set, id, np, tfixed, 1e6 * tl, length(rs), total, median(ts), median(fe),
            median([r.gens for r in rs]), 1e6 * median((ts .- fe .* tl) ./ fe),
            minimum(r.loss for r in rs))
    flush(stdout)
end
