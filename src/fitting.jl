# ABOUTME: Fitting rate equations to experimental data.
# ABOUTME: Constructs FittingProblem, computes loss, and runs multi-start optimization.

"""
    FittingProblem{M, D}

Holds pre-processed experimental data and mechanism info for fitting
rate constants.

# Fields
- `mechanism`: an `AbstractEnzymeMechanism` instance
- `data`: `NamedTuple` of column vectors (via `Tables.columntable`)
- `group_point_indexes`: row indices grouped by unique `group` values
- `Keq`: fixed equilibrium constant
- `scale_k_to_kcat`: a positive `Real` selects relative mode (per-group-centered
  loss); `nothing` selects absolute per-enzyme-turnover mode (uncentered loss)
- `log_abs_rates`: pre-computed `log.(abs.(Rate))`
- `log_ratios_buffer`: pre-allocated working buffer for loss computation
"""
struct FittingProblem{M<:AbstractEnzymeMechanism, D<:NamedTuple}
    mechanism::M
    data::D
    group_point_indexes::Vector{Vector{Int}}
    Keq::Float64
    scale_k_to_kcat::Union{Float64,Nothing}
    log_abs_rates::Vector{Float64}
    log_ratios_buffer::Vector{Float64}
end

"""
    FittingProblem(mechanism::AbstractEnzymeMechanism, table;
        Keq::Real, scale_k_to_kcat::Union{Real,Nothing}=1.0)

Construct a `FittingProblem` from an enzyme mechanism and tabular data.

The table must have columns: `group`, `Rate`, and one column per
metabolite matching `metabolites(mechanism)`. Uses
`Tables.columntable` for conversion.

`scale_k_to_kcat` selects the loss mode: a positive `Real` (default `1.0`)
treats the data as relative (per-group-centered loss); `nothing` treats it as
absolute per-enzyme turnover (uncentered loss).

Rate values must be nonzero (zero rates produce `-Inf` in log space).
"""
function FittingProblem(mechanism::AbstractEnzymeMechanism, table;
        Keq::Real, scale_k_to_kcat::Union{Real,Nothing}=1.0)
    scale_k_to_kcat !== nothing && scale_k_to_kcat <= 0 && error(
        "scale_k_to_kcat must be positive (or nothing); got $scale_k_to_kcat")
    data = Tables.columntable(table)

    mnames = metabolites(mechanism)

    # Validate required columns
    col_names = keys(data)
    for req in (:group, :Rate)
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

    # Build group_point_indexes by grouping on group column
    groups = data.group
    group_map = Dict{eltype(groups), Vector{Int}}()
    for i in 1:n
        key = groups[i]
        if haskey(group_map, key)
            push!(group_map[key], i)
        else
            group_map[key] = [i]
        end
    end
    group_point_indexes = collect(values(group_map))

    # Allocate working buffer
    log_ratios_buffer = Vector{Float64}(undef, n)

    sk = scale_k_to_kcat === nothing ? nothing : Float64(scale_k_to_kcat)
    FittingProblem{typeof(mechanism), typeof(data)}(
        mechanism, data, group_point_indexes, Float64(Keq), sk,
        log_abs_rates, log_ratios_buffer
    )
end

# Accept the concrete working-representation mechanism: compile to the
# singleton once at construction so `loss!`'s hot path operates on the
# @generated `EnzymeMechanism` / `AllostericEnzymeMechanism` (0-alloc).
FittingProblem(mechanism::Union{Mechanism, AllostericMechanism}, table;
        Keq::Real, scale_k_to_kcat::Union{Real,Nothing}=1.0) =
    FittingProblem(compile_mechanism(mechanism), table;
        Keq=Keq, scale_k_to_kcat=scale_k_to_kcat)

"""
    loss!(x::AbstractVector, fp::FittingProblem)

Compute the log-ratio loss. Zero-allocation on the hot path.

Parameters in `x` are in log-space (i.e., actual rate constants =
`exp.(x)`). When `fp.scale_k_to_kcat` is a `Real`, each group's
log-ratios are centered (mean-subtracted) before squaring, making the
loss invariant to per-group E_total scaling (relative data). When it is
`nothing`, the log-ratios are squared without centering, so the absolute
magnitude is scored (absolute per-enzyme turnover data).

Sign mismatches (predicted vs measured rate sign) incur a flat penalty
of 100.0 per point, accumulated after the per-point loop. In the
centered mode this prevents all-mismatch groups from contributing zero
loss (a uniform sentinel would cancel under mean-subtraction).
"""
function loss!(x::AbstractVector, fp::FittingProblem{M,D}) where {M,D}
    buf = fp.log_ratios_buffer
    ParamNames = fitted_params(fp.mechanism)
    MetNames = metabolites(fp.mechanism)
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
        pred = rate_equation(fp.mechanism, concs, params)
        meas_sign = sign(fp.data.Rate[i])
        if sign(pred) != meas_sign || pred == 0.0
            buf[i] = 10.0
        else
            buf[i] = log(abs(pred)) - fp.log_abs_rates[i]
        end
    end

    # Pass 2: loss. Relative (scale_k_to_kcat isa Real) → per-group mean-
    # centered, removing each group's arbitrary scale. Absolute
    # (scale_k_to_kcat === nothing) → uncentered: the y-axis is absolute
    # per-enzyme turnover, so the absolute magnitude is meaningful.
    total_loss = 0.0
    if fp.scale_k_to_kcat === nothing
        @inbounds for i in 1:n_data
            total_loss += buf[i] * buf[i]
        end
    else
        @inbounds for grp_idx in fp.group_point_indexes
            n_grp = length(grp_idx)
            mean_lr = 0.0
            for j in grp_idx
                mean_lr += buf[j]
            end
            mean_lr /= n_grp
            for j in grp_idx
                d = buf[j] - mean_lr
                total_loss += d * d
            end
        end
    end

    # Flat penalty for sign mismatches. In centered mode an all-mismatch group
    # would otherwise contribute zero loss (the uniform 10.0 sentinel cancels
    # under mean-subtraction); the post-hoc penalty keeps it positive. In
    # uncentered mode the sentinel already contributes 100.0 per point, so this
    # adds a second 100.0 (stronger steering away from sign flips, intentional).
    n_mismatch = 0
    @inbounds for i in 1:n_data
        buf[i] == 10.0 && (n_mismatch += 1)
    end

    return (total_loss + 100.0 * n_mismatch) / n_data
end

"""
    fit_rate_equation(fp::FittingProblem, optimizer;
        n_restarts=20, maxtime=60.0, maxiters=10_000_000,
        abstol=nothing, reltol=nothing, callback=nothing,
        lb=fill(-15.0, length(fitted_params(fp.mechanism))),
        ub=fill(15.0, length(fitted_params(fp.mechanism))),
        solver_kwargs=(;))

Fit rate constants by minimizing `loss!` using Optimization.jl.

Runs `n_restarts` independent optimizations from random initial points and returns
the best result.

`maxtime` and `maxiters` are Optimization.jl common solver options forwarded to
every `Optimization.solve` call. `abstol`, `reltol`, and `callback` are common
options forwarded only when set (otherwise each solver uses its own default).
`solver_kwargs` is a `NamedTuple` of solver-specific options forwarded verbatim
to `solve` (e.g. `(; popsize=200)` for a CMA-ES solver that supports it). A key
present in both a named common option and `solver_kwargs` takes the
`solver_kwargs` value. `lb`/`ub` are the log-space bounds passed to the
`OptimizationProblem`.

Rescaling is driven by `fp.scale_k_to_kcat`: when it is a `Real`, the returned
parameters are rescaled so that `_kcat_forward(mechanism, params) ≈
fp.scale_k_to_kcat`; when it is `nothing`, the raw (unrescaled) parameters are
returned (the data fixes the absolute scale).

Returns a NamedTuple `(params, loss, retcode)` where:
- `params`: fitted rate constants as a NamedTuple
- `loss`: the best loss value achieved
- `retcode`: the `Symbol` form of the best restart's `sol.retcode`. Only
  `:Success` indicates the optimizer converged on its own criteria; any other
  value (e.g. `:MaxTime` — hit the time budget; `:Failure`; `:Default` — ran but
  did not flag success) means the fit should be treated as un-converged (check
  `retcode !== :Success`).
"""
function fit_rate_equation(fp::FittingProblem, optimizer;
    n_restarts::Int=20,
    maxtime::Real=60.0,
    maxiters::Integer=10_000_000,
    abstol::Union{Real,Nothing}=nothing,
    reltol::Union{Real,Nothing}=nothing,
    callback=nothing,
    lb=fill(-15.0, length(fitted_params(fp.mechanism))),
    ub=fill(15.0, length(fitted_params(fp.mechanism))),
    solver_kwargs=(;),
)
    obj = Optimization.OptimizationFunction((x, p) -> loss!(x, p))
    np = length(fitted_params(fp.mechanism))

    # Common solver options: maxtime/maxiters always forwarded; the optional
    # ones only when set, so each solver keeps its own default otherwise.
    # solver_kwargs is forwarded verbatim and wins on key conflicts.
    common = (; maxtime, maxiters)
    abstol   === nothing || (common = (; common..., abstol))
    reltol   === nothing || (common = (; common..., reltol))
    callback === nothing || (common = (; common..., callback))
    solve_kwargs = merge(common, solver_kwargs)

    best_x = zeros(np)
    best_loss = Inf
    # Sentinel for "no restart produced a finite objective" (loss stays Inf, so
    # the fit is dropped by the beam's non-finite filter). Deliberately NOT a
    # real SciMLBase ReturnCode name — `Symbol(ReturnCode.Default) === :Default`,
    # so `:NoFiniteLoss` stays distinct from a genuine `:Default` solver return
    # (which only fires when a restart achieves a finite objective).
    best_retcode = :NoFiniteLoss

    for _ in 1:n_restarts
        x0 = clamp.(randn(np) .* 2.0, lb, ub)
        prob = Optimization.OptimizationProblem(obj, x0, fp; lb=lb, ub=ub)
        sol = Optimization.solve(prob, optimizer; solve_kwargs...)
        if sol.objective < best_loss
            best_loss = sol.objective
            best_x .= sol.u
            best_retcode = Symbol(sol.retcode)
        end
    end

    pnames = fitted_params(fp.mechanism)
    result_params = NamedTuple{pnames}(ntuple(i -> exp(best_x[i]), Val(length(pnames))))
    if fp.scale_k_to_kcat !== nothing
        full = merge(result_params, (Keq = fp.Keq, E_total = 1.0))
        rp = rescale_parameter_values(
            fp.mechanism, full; scale_k_to_kcat=fp.scale_k_to_kcat,
        )
        result_params = NamedTuple{pnames}(
            ntuple(i -> rp[pnames[i]], Val(length(pnames))),
        )
    end
    return (params = result_params, loss = best_loss, retcode = best_retcode)
end
