# ABOUTME: Fitting rate equations to experimental data and scaling them to a reference
# ABOUTME: scale: FittingProblem, the loss, multi-start fits and test_rate_equation.

"""
    FittingProblem{M, D}

Holds pre-processed experimental data and mechanism info for fitting
rate constants.

# Fields
- `mechanism`: an `AbstractEnzymeMechanism` instance
- `rate_data`: the validated table, `Keq` (fixed equilibrium constant), `absolute` and
  `reference`, as a `RateData`. `absolute = false` selects relative mode (per-group-centered
  loss, fit scaled to the reference turnovers or, without any, so the largest predicted
  rate over the data rows at `E_total = 1` is 1); `true` selects absolute per-enzyme
  turnover mode (uncentered loss, no scale).
  `reference` holds the anchors `concentrations => turnover` that set a relative fit's
  scale; empty in absolute mode and when a relative table gives none
- `group_point_indexes`: row indices grouped by unique `group` values
- `log_abs_rates`: pre-computed `log.(abs.(Rate))`
- `log_ratios_buffer`: pre-allocated working buffer for loss computation
"""
struct FittingProblem{M<:Mechanisms.AbstractEnzymeMechanism, D<:NamedTuple}
    mechanism::M
    rate_data::Data.RateData{D}
    group_point_indexes::Vector{Vector{Int}}
    log_abs_rates::Vector{Float64}
    log_ratios_buffer::Vector{Float64}
end

"""
    FittingProblem(mechanism::AbstractEnzymeMechanism, table;
        Keq::Real, absolute::Bool = false, reference = nothing)

Construct a `FittingProblem` from an enzyme mechanism and tabular data.

The table must have columns: `group`, `Rate`, and one column per
metabolite matching `metabolites(mechanism)`, so no metabolite may be
named `group` or `Rate`. Uses `Tables.columntable` for conversion.

With `absolute = false` (the default) the data are relative: the loss centers each
group's log-ratios, which removes each group's scale, so the data fix the rate constants
only up to one common factor. `reference`, optional in this mode, sets that factor: a
vector of anchors `concentrations => turnover`, each a `NamedTuple` of concentrations in
the data's units, where a metabolite left out is zero, and a turnover per unit `E_total`
measured there (a kcat, or a specific activity in turnover units).
`fit_rate_equation` scales the fitted rate constants so the mechanism reproduces the
reference turnovers. Without a `reference` it scales them so the largest absolute rate
predicted over the data rows, at `E_total = 1`, is 1: the rate constants are then in
units of the fastest measured condition, not kcats. With `absolute = true` the data
are absolute per-enzyme turnover: the loss is uncentered, the fit is returned unscaled,
and passing `reference` is an error.

Rate values must be finite and nonzero (the loss takes their log; a zero rate gives
`-Inf`), every concentration must be a finite number ≥ 0 (zero is valid), and `Keq`
must be positive and finite. A reference anchor must name metabolites of the mechanism
at finite concentrations ≥ 0, with a finite, nonzero turnover.

```julia
fp = FittingProblem(m, data; Keq = 2.0, reference = [(A = 10.0, B = 5.0) => 42.0])
```
"""
function FittingProblem(@nospecialize(mechanism::Mechanisms.AbstractEnzymeMechanism),
        table; Keq::Real, absolute::Bool=false, reference=nothing)
    mets = Mechanisms.metabolites(mechanism)
    rate_data = Data.RateData(table, mets; Keq, absolute, reference, owner = "mechanism")
    data = rate_data.data
    group_map = Dict{eltype(data.group), Vector{Int}}()
    for (i, g) in enumerate(data.group)
        push!(get!(() -> Int[], group_map, g), i)
    end
    FittingProblem{typeof(mechanism), typeof(data)}(
        mechanism, rate_data, collect(values(group_map)),
        log.(abs.(data.Rate)), Vector{Float64}(undef, length(data.Rate)))
end

"""
Accept the concrete working-representation mechanism: compile to the
singleton once at construction so `loss!`'s hot path operates on the
@generated `EnzymeMechanism` / `AllostericEnzymeMechanism` (0-alloc).
"""
FittingProblem(mechanism::Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism},
               table; kwargs...) =
    FittingProblem(Mechanisms.compile_mechanism(mechanism), table; kwargs...)

"""
    loss!(x::AbstractVector, fp::FittingProblem)

Compute the log-ratio loss. Zero-allocation on the hot path.

Parameters in `x` are in log-space (i.e., actual rate constants =
`exp.(x)`). In relative mode (`fp.rate_data.absolute == false`) each group's
log-ratios are centered (mean-subtracted) before squaring, making the
loss invariant to per-group E_total scaling (relative data). In absolute
mode (`fp.rate_data.absolute`) the log-ratios are squared without centering, so the
absolute magnitude is scored (absolute per-enzyme turnover data).

Sign mismatches (predicted vs measured rate sign) incur a flat penalty
of 100.0 per point, counted in the per-point loop and added after. In the
centered mode this prevents all-mismatch groups from contributing zero
loss (a uniform sentinel would cancel under mean-subtraction).
"""
function loss!(x::AbstractVector, fp::FittingProblem)
    buf = fp.log_ratios_buffer
    ParamNames = Derivation.parameters(fp.mechanism, Derivation.FittedParameters())
    MetNames = Mechanisms.metabolites(fp.mechanism)
    N = length(ParamNames)
    K = length(MetNames)
    n_data = length(fp.log_abs_rates)

    # Build params NamedTuple from log-space x
    fitted = NamedTuple{ParamNames}(ntuple(i -> exp(x[i]), Val(N)))
    params = merge(fitted, (Keq = fp.rate_data.Keq, E_total = 1.0))

    # Pass 1: fill log_ratios_buffer and count sign mismatches. A zero or NaN
    # prediction counts as one: its sign never equals a nonzero rate's ±1.
    n_mismatch = 0
    @inbounds for i in 1:n_data
        concs = NamedTuple{MetNames}(ntuple(
            j -> getproperty(fp.rate_data.data, MetNames[j])[i], Val(K),
        ))
        pred = Derivation.rate_equation(fp.mechanism, concs, params)
        meas_sign = sign(fp.rate_data.data.Rate[i])
        if sign(pred) != meas_sign
            buf[i] = 10.0
            n_mismatch += 1
        else
            buf[i] = log(abs(pred)) - fp.log_abs_rates[i]
        end
    end

    # Pass 2: loss. Relative (fp.rate_data.absolute false) → per-group mean-centered,
    # removing each group's arbitrary scale. Absolute → uncentered: the y-axis is absolute
    # per-enzyme turnover, so the absolute magnitude is meaningful.
    total_loss = 0.0
    if fp.rate_data.absolute
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
    return (total_loss + 100.0 * n_mismatch) / n_data
end

"""
The objective `fit_rate_equation` hands the optimizer. `loss` holds `loss!(x, fp)` as a
function of `x` behind a `FunctionWrapper`, which hides the mechanism and data types:
every fit gives `Optimization.solve` one problem type, so the solver compiles once per
optimizer rather than once per mechanism. `x_buf` takes a copy of each point before
`loss` scores it, because CMA-ES passes a view and the wrapper takes a `Vector{Float64}`.
"""
struct FitObjective
    loss::FunctionWrapper{Float64, Tuple{Vector{Float64}}}
    x_buf::Vector{Float64}
end

"""The optimizer's objective: `obj.loss` at a copy of `x` (see `FitObjective`)."""
_fit_loss(x, obj::FitObjective) = obj.loss(copyto!(obj.x_buf, x))

"""
The `OptimizationProblem` minimizing `loss!` of `fp` from `x0` within `lb`/`ub`. Its
type is the same for every mechanism and data table (see `FitObjective`).
"""
function _optimization_problem(@nospecialize(fp::FittingProblem), x0, lb, ub)
    # `Base.Fix2` takes the type of `fp` from its value. A closure over the
    # `@nospecialize` argument may store it untyped (Julia 1.13 does), and then each
    # call dispatches `loss!` at run time and boxes its result.
    obj = FitObjective(Base.Fix2(loss!, fp), Vector{Float64}(undef, length(x0)))
    Optimization.OptimizationProblem(Optimization.OptimizationFunction(_fit_loss), x0, obj;
                                     lb=lb, ub=ub)
end

# Distance in log units from a box edge within which a fitted parameter counts as at bound.
const AT_BOUND_TOLERANCE = 0.01

"""
    fit_rate_equation(fp::FittingProblem, optimizer;
        n_restarts=20, maxtime=60.0, maxiters=10_000_000,
        abstol=nothing, reltol=nothing, callback=nothing,
        lb=centre .- 15.0, ub=centre .+ 15.0,
        solver_kwargs=(;))

Fit rate constants by minimizing `loss!` using Optimization.jl.

Runs `n_restarts` independent optimizations from random initial points and returns
the best result. Each start draws every log-space parameter from a normal distribution
around its `centre`, with a standard deviation of 2, and clamps it to the box. The centre
comes from the geometric means of the table's nonzero concentrations and of its `|Rate|`:
a binding's K centres on the concentrations it binds, a rate constant on the rate over
the concentrations it takes up, a regulator's K on its ligand's concentrations and `L`
on 1.

`maxtime` and `maxiters` are Optimization.jl common solver options forwarded to
every `Optimization.solve` call. `abstol`, `reltol`, and `callback` are common
options forwarded only when set (otherwise each solver uses its own default).
`solver_kwargs` is a `NamedTuple` of solver-specific options forwarded verbatim
to `solve` (e.g. `(; popsize=200)` for a CMA-ES solver that supports it). A key
present in both a named common option and `solver_kwargs` takes the
`solver_kwargs` value. `lb`/`ub` are the log-space bounds passed to the
`OptimizationProblem`. They hold absolute log values, and by default lie 15 log units
either side of each centre, so the box follows the data.

In relative mode (`fp.rate_data.absolute == false`) the returned steady-state rate
constants are multiplied by one factor so the mechanism reproduces the reference turnovers
`fp.rate_data.reference` (see `FittingProblem`): the geometric mean of measured over
predicted turnover. Without a reference the factor is 1 over the largest absolute rate
the fit predicts at the data's concentrations, at `E_total = 1`. A fit whose predicted
turnover at a reference anchor is zero, non-finite or of the other sign, or whose largest
predicted rate is zero or non-finite, has no such factor: it returns its parameters
unscaled, with `loss = Inf` and `retcode = :ReferenceFailed`, and does not throw.
In absolute mode, or when no restart reaches a finite loss (`loss = Inf`,
`retcode = :NoFiniteLoss`), the rate constants are returned as fitted.

Returns a NamedTuple `(params, loss, retcode, at_bound)` where:
- `params`: fitted rate constants as a NamedTuple
- `loss`: the best loss value achieved
- `retcode`: the `Symbol` form of the best restart's `sol.retcode`. Only
  `:Success` indicates the optimizer converged on its own criteria; any other
  value (e.g. `:MaxTime` — hit the time budget; `:Failure`; `:Default` — ran but
  did not flag success) means the fit should be treated as un-converged (check
  `retcode !== :Success`); `:ReferenceFailed` marks a failed reference scale.
- `at_bound`: the names of the fitted parameters whose best point lies within 0.01 log
  units of its `lb` or `ub`, read before the relative-mode scale, which can move a
  reported rate constant outside the box. `()` when no restart reaches a finite loss. "At
  bound" does not mean wrong: a parameter no data constrain wanders across the box.
"""
function fit_rate_equation(@nospecialize(fp::FittingProblem), optimizer;
    n_restarts::Int=20,
    maxtime::Real=60.0,
    maxiters::Integer=10_000_000,
    abstol::Union{Real,Nothing}=nothing,
    reltol::Union{Real,Nothing}=nothing,
    callback=nothing,
    lb=_log_centres(fp) .- 15.0,
    ub=_log_centres(fp) .+ 15.0,
    solver_kwargs=(;),
)
    # `fp` is @nospecialize, so `parameters(fp.mechanism, FittedParameters())` infers as
    # `Any`; asserting a tuple of Symbols lets `np` infer as an `Int` and the NamedTuple
    # built from the names as a `NamedTuple`.
    pnames = Derivation.parameters(fp.mechanism,
                                   Derivation.FittedParameters())::Tuple{Vararg{Symbol}}
    np = length(pnames)

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
    # `_ingest!` drops it from the beam search). Deliberately NOT a
    # real SciMLBase ReturnCode name — `Symbol(ReturnCode.Default) === :Default`,
    # so `:NoFiniteLoss` stays distinct from a genuine `:Default` solver return
    # (which only fires when a restart achieves a finite objective).
    best_retcode = :NoFiniteLoss

    centre = _log_centres(fp)
    for _ in 1:n_restarts
        x0 = clamp.(centre .+ randn(np) .* 2.0, lb, ub)
        prob = _optimization_problem(fp, x0, lb, ub)
        sol = Optimization.solve(prob, optimizer; solve_kwargs...)
        if sol.objective < best_loss
            best_loss = sol.objective
            best_x .= sol.u
            best_retcode = Symbol(sol.retcode)
        end
    end

    result_params = NamedTuple{pnames}(ntuple(i -> exp(best_x[i]), Val(np)))
    at_bound = isfinite(best_loss) ?
        Tuple(pnames[i] for i in 1:np
              if min(best_x[i] - lb[i], ub[i] - best_x[i]) <= AT_BOUND_TOLERANCE) : ()
    (fp.rate_data.absolute || !isfinite(best_loss)) &&
        return (params = result_params, loss = best_loss, retcode = best_retcode, at_bound)
    scale = _reference_scale(fp, merge(result_params, (Keq = fp.rate_data.Keq,)))
    scale === nothing &&
        return (params = result_params, loss = Inf, retcode = :ReferenceFailed, at_bound)
    return (params = _scale_rate_constants(fp.mechanism, result_params, scale),
            loss = best_loss, retcode = best_retcode, at_bound)
end

"""
    _concentration_scales(fp) -> Dict{Symbol, Float64}

Each metabolite of `fp.mechanism` to the geometric mean of its nonzero concentrations in
the table. A metabolite that is always zero takes the geometric mean of all positive
concentrations of those metabolites in the table, or 1 when there are none.
"""
function _concentration_scales(@nospecialize(fp::FittingProblem))
    data = fp.rate_data.data
    mets = Mechanisms.metabolites(fp.mechanism)
    # `getfield` reads a column by its runtime name without compiling for the table's type.
    positive = Vector{Float64}[Float64[c for c in getfield(data, x) if c > 0] for x in mets]
    pooled = reduce(vcat, positive)
    fallback = isempty(pooled) ? 1.0 : _geometric_mean(pooled)
    Dict{Symbol, Float64}(x => isempty(v) ? fallback : _geometric_mean(v)
                          for (x, v) in zip(mets, positive))
end

_geometric_mean(v) = exp(mean(log, v))

"""The geometric mean of the table's `|Rate|`, which `RateData` guarantees nonzero."""
_rate_scale(@nospecialize(fp::FittingProblem)) = exp(mean(fp.log_abs_rates))

"""
    _log_centre(p, m, s, r) -> Float64

The log value of parameter `p` of the concrete mechanism `m` at which it matches the
concentration scales `s` (`_concentration_scales`) and the rate scale `r` (`_rate_scale`).
A rate constant turns the scales of what it takes up into `r`. A binding's K, a dissociation
constant, is the scale of what binds, and any other rapid-equilibrium K is the scales given
off over those taken up. A regulator's K is the scale of its ligand, and `L` centres on 1.
"""
function _log_centre(p::Derivation.StepParameter, m, s::Dict{Symbol, Float64}, r::Real)
    cm = Mechanisms.catalytic(m)
    e = cm.steps[p.step]
    e.rapid_equilibrium ||
        return log(r) - _sum_log(s, p.slot == 1 ? e.taken_up : e.given_off)
    cm.step_kinds[p.step] in (:binding, :fused_binding) ? _sum_log(s, e.taken_up) :
        _sum_log(s, e.given_off) - _sum_log(s, e.taken_up)
end
_log_centre(p::Derivation.RegulatorParameter, m, s::Dict{Symbol, Float64}, r::Real) =
    log(s[p.ligand])
_log_centre(::Derivation.AllostericParameter, m, s::Dict{Symbol, Float64}, r::Real) = 0.0

_sum_log(s, mets) = sum(log(s[x.name]) for x in mets; init = 0.0)

"""
    _log_centres(fp) -> Vector{Float64}

The centre of `fit_rate_equation`'s search: `_log_centre` of each parameter of
`Derivation.parameter_info`, in the order of `parameters(fp.mechanism, FittedParameters())`.
"""
function _log_centres(@nospecialize(fp::FittingProblem))
    m = _concrete(fp.mechanism)::Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}
    s = _concentration_scales(fp)
    r = _rate_scale(fp)
    Float64[_log_centre(p, m, s, r) for p in Derivation.parameter_info(m)]
end

"""
    _scale_rate_constants(mechanism, params, factor) -> NamedTuple

`params` with every fitted rate constant on a steady-state step of `mechanism` (a concrete
or a compiled one) multiplied by `factor`, every other entry unchanged. A common factor on
those constants multiplies every rate by the same factor and leaves the binding constants,
`Keq`, `E_total`, `L` and the regulatory constants alone.
"""
function _scale_rate_constants(@nospecialize(mechanism), params::NamedTuple, factor::Real)
    m = _concrete(mechanism)
    steady_state = Set(Derivation.name(p, m) for p in Derivation.parameter_info(m)
                       if _on_steady_state_step(p, m))
    NamedTuple{keys(params)}(map((k, v) -> k in steady_state ? v * factor : v,
                                 keys(params), values(params)))
end

_on_steady_state_step(p::Derivation.StepParameter, m) =
    !Mechanisms.catalytic(m).steps[p.step].rapid_equilibrium
_on_steady_state_step(p, m) = false

_concrete(m::Mechanisms.EnzymeMechanism) = Mechanisms.Mechanism(m)
_concrete(m::Mechanisms.AllostericEnzymeMechanism) = Mechanisms.AllostericMechanism(m)
_concrete(m) = m

"""
    _reference_scale(fp, params)

The factor on the rate constants `params` (which hold `Keq`) of `fp.mechanism` that sets
a relative fit's scale, or `nothing` if no factor does. With reference anchors it is the
geometric mean of measured over predicted turnover per unit `E_total`, where a metabolite
an anchor leaves out is zero, and no factor reproduces an anchor whose predicted turnover
is zero, non-finite or of the other sign than the measured one. Without anchors it is
1 over the largest absolute rate predicted at the data's concentrations, which must be
positive and finite.
"""
function _reference_scale(@nospecialize(fp::FittingProblem),
                          @nospecialize(params::NamedTuple))
    mets = Mechanisms.metabolites(fp.mechanism)
    turnover_params = merge(params, (E_total = 1.0,))
    reference = fp.rate_data.reference
    if isempty(reference)
        peak = maximum(
            abs(Derivation.rate_equation(fp.mechanism,
                NamedTuple{mets}(Tuple(Float64(fp.rate_data.data[k][i]) for k in mets)),
                turnover_params))
            for i in eachindex(fp.rate_data.data.Rate))
        return isfinite(peak) && peak > 0 ? 1 / peak : nothing
    end
    log_ratios = Float64[]
    for (concs, measured) in reference
        predicted = Derivation.rate_equation(fp.mechanism,
            NamedTuple{mets}(Tuple(Float64(get(concs, k, 0.0)) for k in mets)),
            turnover_params)
        isfinite(predicted) && predicted != 0 && sign(predicted) == sign(measured) ||
            return nothing
        push!(log_ratios, log(measured / predicted))
    end
    exp(sum(log_ratios) / length(log_ratios))
end

"""
    test_rate_equation(m, params::NamedTuple, table; Keq, absolute = false) -> Float64

The fitting loss of mechanism `m` at the rate constants `params` on `table`: the loss
`fit_rate_equation` minimizes, for scoring a fit on data it did not see. `params` must
hold every name of `parameters(m, FittedParameters())`, as a fit's `params` does; other
names are ignored. `table` has the columns `FittingProblem` needs. `Keq` is the fixed
equilibrium constant. With `absolute = false` (the default) the loss centers each group's
log-ratios, so it does not depend on the scale of `params`; with `absolute = true` it
scores the absolute turnover.

```julia
fit = fit_rate_equation(FittingProblem(m, train; Keq = 2.0), optimizer)
EnzymeRates.test_rate_equation(m, fit.params, test; Keq = 2.0)
```
"""
function test_rate_equation(@nospecialize(m), @nospecialize(params::NamedTuple), table;
                            Keq::Real, absolute::Bool = false)
    fp = FittingProblem(m, table; Keq, absolute)
    names = Derivation.parameters(fp.mechanism, Derivation.FittedParameters())
    for name in names
        haskey(params, name) || error("params is missing the fitted parameter $name")
    end
    loss!([log(params[name]) for name in names], fp)
end
