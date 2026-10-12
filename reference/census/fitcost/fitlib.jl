# ABOUTME: Fit-cost calibration helpers: synthetic relative-mode data from a reference mechanism, and
# ABOUTME: identify's multi-start CMA-ES fit replicated with per-restart wall time, generations and stop reason.
using Optimization, OptimizationCMAEvolutionStrategy, Random, Printf, Statistics, LinearAlgebra

const SERIES = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]

"""The data design of a reaction: `(group, concentrations)` rows, 48 per reaction in 4
groups, each series spanning 0.05 to 10 (the reference constants lie in [0.32, 3.2])."""
function design(name::Symbol)
    rows = Tuple{String, Dict{Symbol, Float64}}[]
    add!(g; kw...) = push!(rows, (g, Dict{Symbol, Float64}(kw)))
    if name === :bi_bi_rxn
        for b in (0.3, 3.0), a in SERIES
            add!("G1"; A = a, B = b)
        end
        for a in (0.3, 3.0), b in SERIES
            add!("G2"; A = a, B = b)
        end
        for a in SERIES
            add!("G3"; A = a, B = 1.0, P = 1.0)
        end
        for b in SERIES
            add!("G4"; A = 1.0, B = b, Q = 1.0)
        end
    elseif name === :ter_bi_rxn
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
            add!("G3"; A = a, B = 1.0, D = 1.0, P = 1.0)
        end
        for d in SERIES
            add!("G4"; A = 1.0, B = 1.0, D = d, Q = 1.0)
        end
    else
        error("no data design for $name")
    end
    rows
end

"""Relative-mode data from the reference mechanism `em` at constants `θ`: each group scaled
by its own factor 10^(0.3·N(0,1)) (an unknown E_total) and each rate by exp(0.05·N(0,1))."""
function make_data(name, em, θ::NamedTuple, keq, rng; noise = 0.05, gscale = 0.3)
    mets = ER.metabolites(em)
    rows = design(name)
    scale = Dict(g => 10.0^(gscale * randn(rng)) for g in unique(first.(rows)))
    params = merge(θ, (Keq = keq, E_total = 1.0))
    cols = Dict(x => Float64[] for x in mets)
    G, R, V = String[], Float64[], Float64[]
    for (g, c) in rows
        concs = NamedTuple{mets}(Tuple(get(c, x, 0.0) for x in mets))
        v = ER.rate_equation(em, concs, params)
        (isfinite(v) && v > 0) || error("reference rate $v at $c")
        push!(G, g)
        push!(V, v)
        push!(R, v * scale[g] * exp(noise * randn(rng)))
        for x in mets
            push!(cols[x], get(c, x, 0.0))
        end
    end
    data = merge((group = G, Rate = R), NamedTuple{mets}(Tuple(cols[x] for x in mets)))
    data, V, scale
end

"""
`fit_rate_equation`'s restart loop with each restart timed: the same draws from the global
RNG (`x0 = clamp.(randn(np) .* 2, lb, ub)`, then CMA-ES's own seed), the same
`_optimization_problem` and the same `solve` options. Stops early (`truncated`) once the
restarts so far took more than `fit_cap` seconds. `fp` is `@nospecialize`, as in
`fit_rate_equation`, so this loop compiles once rather than once per mechanism.
"""
function fit_instrumented(@nospecialize(fp), opt; n_restarts = 20, maxtime = 600.0, maxiters = 10_000_000,
                          fit_cap = Inf)
    np = length(ER.fitted_params(fp.mechanism))
    lb, ub = fill(-15.0, np), fill(15.0, np)
    best_loss, best_x, best_ret = Inf, zeros(np), :NoFiniteLoss
    rs = NamedTuple[]
    total = 0.0
    truncated = false
    for r in 1:n_restarts
        x0 = clamp.(randn(np) .* 2.0, lb, ub)
        prob = ER._optimization_problem(fp, x0, lb, ub)
        t = time_ns()
        sol = Optimization.solve(prob, opt; maxtime = maxtime, maxiters = maxiters)
        dt = (time_ns() - t) / 1e9
        total += dt
        o = sol.original
        gens = sol.stats.iterations
        push!(rs, (restart = r, t = dt, gens = gens, lambda = o.p.λ, fevals = gens * o.p.λ,
                   retcode = Symbol(sol.retcode), reason = o.stop.reason,
                   loss = sol.objective))
        if sol.objective < best_loss
            best_loss = sol.objective
            best_x .= sol.u
            best_ret = Symbol(sol.retcode)
        end
        if total > fit_cap && r < n_restarts
            truncated = true
            break
        end
    end
    (loss = best_loss, x = best_x, retcode = best_ret, restarts = rs, truncated = truncated)
end

"""Median and minimum seconds per `loss!` call at `nx` random points (each timed over
enough calls to take about `target` seconds)."""
function time_loss(fp, np, rng; nx = 5, target = 0.02)
    xs = [clamp.(randn(rng, np) .* 2.0, -15.0, 15.0) for _ in 1:nx]
    acc = 0.0
    t1 = @elapsed for _ in 1:50
        acc += ER.loss!(xs[1], fp)
    end
    n = clamp(round(Int, target / max(t1 / 50, 1e-9)), 50, 500_000)
    ts = Float64[]
    for x in xs
        t = time_ns()
        for _ in 1:n
            acc += ER.loss!(x, fp)
        end
        push!(ts, (time_ns() - t) / 1e9 / n)
    end
    median(ts), minimum(ts), acc
end

"""CMA-ES cost per generation apart from the objective: `gens` generations of a sphere in
`n` dimensions (far from its minimum, so no stop rule fires), through Optimization.jl."""
function cma_overhead(opt; ns = (2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 18, 20),
                      gens = 400, reps = 3)
    sphere(x, p) = (s = 0.0; for i in eachindex(x); s += (x[i] - p[i])^2; end; s)
    f = Optimization.OptimizationFunction(sphere)
    out = NamedTuple[]
    for n in ns, rep in 1:reps
        prob = Optimization.OptimizationProblem(f, fill(5.0, n), fill(-5.0, n);
                                                lb = fill(-15.0, n), ub = fill(15.0, n))
        t = time_ns()
        sol = Optimization.solve(prob, opt; maxiters = gens, maxtime = 60.0)
        dt = (time_ns() - t) / 1e9
        push!(out, (n = n, rep = rep, t = dt, gens = sol.stats.iterations,
                    lambda = sol.original.p.λ, reason = sol.original.stop.reason))
    end
    out
end

ncalls(x) = x isa Expr ? (x.head === :call) + sum(ncalls, x.args; init = 0) : 0
