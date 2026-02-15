# ─── Fitting rate equations to experimental data ─────────────────────

"""
    FittingProblem{M, D}

Holds pre-processed experimental data and mechanism info for fitting rate constants.

# Fields
- `mechanism`: the `EnzymeMechanism` instance
- `data`: `NamedTuple` of column vectors (via `Tables.columntable`)
- `fig_point_indexes`: row indices grouped by unique `(Article, Fig)` pairs
- `Keq`: fixed equilibrium constant
- `log_abs_rates`: pre-computed `log.(abs.(Rate))`
- `log_ratios_buffer`: pre-allocated working buffer for loss computation
"""
struct FittingProblem{M<:EnzymeMechanism, D<:NamedTuple}
    mechanism::M
    data::D
    fig_point_indexes::Vector{Vector{Int}}
    Keq::Float64
    log_abs_rates::Vector{Float64}
    log_ratios_buffer::Vector{Float64}
end

"""
    FittingProblem(mechanism::EnzymeMechanism, table; Keq::Real)

Construct a `FittingProblem` from an enzyme mechanism and tabular data.

The table must have columns: `Article`, `Fig`, `Rate`, and one column per metabolite
in `metabolites(mechanism)`. Uses `Tables.columntable` for conversion.

Rate values must be nonzero (zero rates produce `-Inf` in log space).
"""
function FittingProblem(mechanism::EnzymeMechanism, table; Keq::Real)
    data = Tables.columntable(table)

    mnames = metabolite_names(mechanism)

    # Validate required columns
    col_names = keys(data)
    for req in (:Article, :Fig, :Rate)
        req in col_names || error("Missing required column: $req")
    end
    for m in mnames
        m in col_names || error("Missing metabolite column: $m")
    end

    # Validate no zero rates
    rates = data.Rate
    n = length(rates)
    for i in 1:n
        rates[i] == 0 && error("Zero rate at row $i: log(0) is undefined")
    end

    # Pre-compute log(abs(rates))
    log_abs_rates = log.(abs.(rates))

    # Build fig_point_indexes by grouping on (Article, Fig) pairs
    articles = data.Article
    figs = data.Fig
    group_map = Dict{Tuple{eltype(articles),eltype(figs)}, Vector{Int}}()
    for i in 1:n
        key = (articles[i], figs[i])
        if haskey(group_map, key)
            push!(group_map[key], i)
        else
            group_map[key] = [i]
        end
    end
    fig_point_indexes = collect(values(group_map))

    # Allocate working buffer
    log_ratios_buffer = Vector{Float64}(undef, n)

    FittingProblem{typeof(mechanism), typeof(data)}(
        mechanism, data, fig_point_indexes, Float64(Keq),
        log_abs_rates, log_ratios_buffer
    )
end

"""
    loss!(x::AbstractVector, fp::FittingProblem)

Compute the per-figure-centered log-ratio loss. Zero-allocation on the hot path.

Parameters in `x` are in log-space (i.e., actual rate constants = `exp.(x)`).
Each figure's log-ratios are centered (mean-subtracted) before squaring,
making the loss invariant to per-figure E_total scaling.

Sign mismatches (predicted vs measured rate sign) incur a penalty of 10.0.
"""
function loss!(x::AbstractVector, fp::FittingProblem{M,D}) where {M,D}
    buf = fp.log_ratios_buffer
    ParamNames = fitted_params(fp.mechanism)
    MetNames = metabolite_names(fp.mechanism)
    N = length(ParamNames)
    K = length(MetNames)
    n_data = length(fp.log_abs_rates)

    # Build params NamedTuple from log-space x
    fitted = NamedTuple{ParamNames}(ntuple(i -> exp(x[i]), Val(N)))
    params = merge(fitted, (Keq = fp.Keq, E_total = 1.0))

    # Pass 1: fill log_ratios_buffer
    @inbounds for i in 1:n_data
        concs = NamedTuple{MetNames}(ntuple(
            j -> getproperty(fp.data, MetNames[j])[i], Val(K),
        ))
        pred = rate_equation(fp.mechanism, params, concs)
        meas_sign = sign(fp.data.Rate[i])
        if sign(pred) != meas_sign || pred == 0.0
            buf[i] = 10.0
        else
            buf[i] = log(abs(pred)) - fp.log_abs_rates[i]
        end
    end

    # Pass 2: per-figure centered loss
    total_loss = 0.0
    @inbounds for fig_idx in fp.fig_point_indexes
        n_fig = length(fig_idx)
        mean_lr = 0.0
        for j in fig_idx
            mean_lr += buf[j]
        end
        mean_lr /= n_fig
        for j in fig_idx
            d = buf[j] - mean_lr
            total_loss += d * d
        end
    end

    return total_loss / n_data
end

"""
    fit_rate_equation(fp::FittingProblem, optimizer;
        n_restarts=10, maxtime=60.0,
        lb=fill(-15.0, length(fitted_params(fp.mechanism))),
        ub=fill(15.0, length(fitted_params(fp.mechanism))),
        kwargs...)

Fit rate constants by minimizing `loss!` using Optimization.jl.

Runs `n_restarts` independent optimizations from random initial points and returns
the best result.

Returns a NamedTuple `(params, loss, x)` where:
- `params`: fitted rate constants as a NamedTuple (in real space, i.e., `exp.(x)`)
- `loss`: the best loss value achieved
- `x`: the best parameter vector in log-space
"""
function fit_rate_equation(fp::FittingProblem, optimizer;
    n_restarts::Int=10,
    maxtime::Real=60.0,
    lb=fill(-15.0, length(fitted_params(fp.mechanism))),
    ub=fill(15.0, length(fitted_params(fp.mechanism))),
    kwargs...
)
    obj = Optimization.OptimizationFunction((x, p) -> loss!(x, p))
    np = length(fitted_params(fp.mechanism))

    best_x = zeros(np)
    best_loss = Inf

    for _ in 1:n_restarts
        x0 = clamp.(randn(np) .* 2.0, lb, ub)
        prob = Optimization.OptimizationProblem(obj, x0, fp; lb=lb, ub=ub)
        sol = Optimization.solve(prob, optimizer; maxtime=maxtime, kwargs...)
        if sol.objective < best_loss
            best_loss = sol.objective
            best_x .= sol.u
        end
    end

    pnames = fitted_params(fp.mechanism)
    result_params = NamedTuple{pnames}(ntuple(i -> exp(best_x[i]), Val(length(pnames))))
    return (params = result_params, loss = best_loss, x = copy(best_x))
end
