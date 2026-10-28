# ABOUTME: Fitting rate equations to experimental data and scaling them to reference
# ABOUTME: turnovers: FittingProblem, the loss, multi-start fits and rescaled parameters.

"""
    FittingProblem{M, D}

Holds pre-processed experimental data and mechanism info for fitting
rate constants.

# Fields
- `mechanism`: an `AbstractEnzymeMechanism` instance
- `data`: `NamedTuple` of column vectors (via `Tables.columntable`)
- `group_point_indexes`: row indices grouped by unique `group` values
- `Keq`: fixed equilibrium constant
- `absolute`: `false` selects relative mode (per-group-centered loss, fit scaled to the
  reference turnovers); `true` selects absolute per-enzyme turnover mode (uncentered
  loss, no scale)
- `reference`: the anchors `concentrations => turnover` that set a relative fit's scale;
  empty in absolute mode
- `log_abs_rates`: pre-computed `log.(abs.(Rate))`
- `log_ratios_buffer`: pre-allocated working buffer for loss computation
"""
struct FittingProblem{M<:Mechanisms.AbstractEnzymeMechanism, D<:NamedTuple}
    mechanism::M
    data::D
    group_point_indexes::Vector{Vector{Int}}
    Keq::Float64
    absolute::Bool
    reference::Vector{Pair{NamedTuple, Float64}}
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
only up to one common factor. `reference`, required in this mode, sets that factor: a
vector of anchors `concentrations => turnover`, each a `NamedTuple` of concentrations in
the data's units, where a metabolite left out is zero, and a turnover per unit `E_total`
measured there (a kcat, or a specific activity in turnover units).
`fit_rate_equation` scales the fitted rate constants so the mechanism reproduces the
reference turnovers, as `rescale_parameter_values` does. With `absolute = true` the data
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
function FittingProblem(mechanism::Mechanisms.AbstractEnzymeMechanism, table;
        Keq::Real, absolute::Bool=false, reference=nothing)
    mets = Mechanisms.metabolites(mechanism)
    (; data, reference) = rate_data(table, mets; Keq, absolute, reference,
                                    owner = "mechanism")
    group_map = Dict{eltype(data.group), Vector{Int}}()
    for (i, g) in enumerate(data.group)
        push!(get!(() -> Int[], group_map, g), i)
    end
    # The typed constructor converts Keq to the field type.
    FittingProblem{typeof(mechanism), typeof(data)}(
        mechanism, data, collect(values(group_map)), Keq, absolute, reference,
        log.(abs.(data.Rate)), Vector{Float64}(undef, length(data.Rate)))
end

"""
    rate_data(table, metabolite_names; Keq, absolute, reference, owner) → (; data, reference)

Validate a rate table and its reference anchors. `data` is `Tables.columntable(table)`:
`Keq` must be positive and finite, the table needs a `group` column, a `Rate` column and
one column per name in `metabolite_names`, no metabolite may be named `group` or `Rate`,
every concentration must be a finite number ≥ 0 (zero is valid), and every rate must be a
finite, nonzero number (the loss takes its log). `reference` is what a problem stores: in
relative mode (`absolute = false`) the anchors, which it requires, checked against
`metabolite_names` of the `owner` (`"mechanism"` or `"reaction"`), with Float64 turnovers;
in absolute mode, where `reference` must be `nothing`, an empty vector.
"""
function rate_data(table, metabolite_names; Keq, absolute::Bool, reference, owner::String)
    (; data = _rate_table(table, metabolite_names, Keq),
       reference = _reference(reference, absolute, metabolite_names, owner))
end

# The validated `Tables.columntable(table)`; see `rate_data`.
function _rate_table(table, mnames, Keq)
    0 < Keq < Inf || error("Keq must be positive and finite; got $Keq")
    data = Tables.columntable(table)
    for req in (:group, :Rate)
        req in mnames && error("Metabolite $req has the name of the required $req " *
                               "column; rename the metabolite")
        req in keys(data) || error("Missing required column: $req")
    end
    for m in mnames
        m in keys(data) || error("Missing metabolite column: $m")
        i = findfirst(c -> !(c isa Real && isfinite(c) && c >= 0), data[m])
        i === nothing || error("Concentration $m at row $i must be a finite number ≥ 0; " *
                               "got $(repr(data[m][i]))")
    end
    i = findfirst(r -> !(r isa Real && isfinite(r)), data.Rate)
    i === nothing ||
        error("Rate at row $i must be a finite number; got $(repr(data.Rate[i]))")
    i = findfirst(iszero, data.Rate)
    i === nothing || error("Zero rate at row $i: log(0) is undefined")
    data
end

# The reference anchors a problem stores; see `rate_data`.
function _reference(reference, absolute::Bool, mets, owner::String)
    if absolute
        reference === nothing || error(
            "reference applies only to relative mode; absolute turnover data set the " *
            "rate scale themselves")
        return Pair{NamedTuple, Float64}[]
    end
    reference === nothing && error(
        "relative mode needs reference = [(concentrations...) => turnover, ...] to set " *
        "the rate scale; pass absolute = true for absolute per-enzyme turnover data")
    _check_anchors(reference, mets, "reference", owner)
    Pair{NamedTuple, Float64}[concs => Float64(t) for (concs, t) in reference]
end

"""
    _check_anchors(anchors, mets, name, owner)

Error unless `anchors`, the argument `name` (`"reference"` or `"kcats"`), is a nonempty
vector of `concentrations => turnover` pairs whose concentrations are a `NamedTuple` of
finite numbers ≥ 0 keyed by names in `mets`, the metabolites of the `owner`, and whose
turnover is a finite, nonzero number.
"""
function _check_anchors(anchors, mets, name::String, owner::String)
    anchors isa AbstractVector{<:Pair} || error("$name must be a vector of " *
        "(concentrations => turnover) pairs; got $(repr(anchors))")
    isempty(anchors) &&
        error("$name is empty; pass at least one (concentrations => turnover) anchor")
    for (concs, measured) in anchors
        concs isa NamedTuple ||
            error("$name: concentrations must be a NamedTuple; got $(repr(concs))")
        for (k, c) in pairs(concs)
            k in mets || error("$name: $k is not a metabolite of the $owner " *
                               "(metabolites: $(join(mets, ", ")))")
            c isa Real && isfinite(c) && c >= 0 || error(
                "$name: concentration of $k must be a finite number ≥ 0; got $(repr(c))")
        end
        measured isa Real && isfinite(measured) && measured != 0 || error(
            "$name: turnover must be a finite, nonzero number; got $(repr(measured))")
    end
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
    _fold_problem(mechanism, data, mask; keywords...) → FittingProblem

The `FittingProblem` of `mechanism` on the rows of `data` that `mask` selects, with the
`FittingProblem` `keywords`. The fold copies its rows rather than viewing them, so its
problem has the full fit's type and `loss!` does not compile a second time for it; the
copy is small next to a fit.
"""
function _fold_problem(@nospecialize(mechanism::Mechanisms.AbstractEnzymeMechanism),
                       data::NamedTuple, mask; keywords...)
    FittingProblem(mechanism, map(col -> col[mask], data); keywords...)
end

"""
    cross_validation_loss(mechanism, data, held_out; Keq, absolute, reference = nothing,
        optimizer, kwargs...) → Float64

The loss on the rows of `data` whose `group` is `held_out` of `mechanism` fitted to the
other rows: `fit_rate_equation(problem, optimizer; kwargs...)` on a `FittingProblem` with
`Keq`, `absolute` and `reference`, which absolute mode leaves `nothing`. A training fit
that fails its reference (`:ReferenceFailed`) keeps its parameters unscaled, and the
relative loss, blind to the scale, still scores them. A non-finite loss raises (naming the
held-out group): a corrupted fold must abort model selection rather than propagate a bad
score.

This function and `_fold_problem` take `mechanism` `@nospecialize`. Cross-validation runs
a candidate's folds after its full fit, which compiled `loss!` for the folds' problem type,
so a fold infers and compiles nothing for the mechanism.
"""
function cross_validation_loss(
    @nospecialize(mechanism::Mechanisms.AbstractEnzymeMechanism),
    data::NamedTuple, held_out;
    Keq, absolute, reference = nothing, optimizer, kwargs...)
    held = data.group .== held_out
    problem_keywords = (; Keq, absolute, reference)
    fit = fit_rate_equation(_fold_problem(mechanism, data, .!held; problem_keywords...),
                            optimizer; kwargs...)
    # `fit.params` is keyed by `fitted_params(mechanism)`, the order `loss!` reads.
    test_loss = loss!([log(v) for v in fit.params],
                      _fold_problem(mechanism, data, held; problem_keywords...))
    # A non-finite fold loss means the fit is unusable; aborting model
    # selection is correct (re-run CV from the saved CSVs after fixing
    # the fit).
    isfinite(test_loss) || error(
        "LOOCV produced a non-finite test loss for held-out group " *
        "$held_out — the fit is unusable; aborting model selection.")
    test_loss
end

"""
    loss!(x::AbstractVector, fp::FittingProblem)

Compute the log-ratio loss. Zero-allocation on the hot path.

Parameters in `x` are in log-space (i.e., actual rate constants =
`exp.(x)`). In relative mode (`fp.absolute == false`) each group's
log-ratios are centered (mean-subtracted) before squaring, making the
loss invariant to per-group E_total scaling (relative data). In absolute
mode (`fp.absolute`) the log-ratios are squared without centering, so the absolute
magnitude is scored (absolute per-enzyme turnover data).

Sign mismatches (predicted vs measured rate sign) incur a flat penalty
of 100.0 per point, counted in the per-point loop and added after. In the
centered mode this prevents all-mismatch groups from contributing zero
loss (a uniform sentinel would cancel under mean-subtraction).
"""
function loss!(x::AbstractVector, fp::FittingProblem)
    buf = fp.log_ratios_buffer
    ParamNames = Derivation.fitted_params(fp.mechanism)
    MetNames = Mechanisms.metabolites(fp.mechanism)
    N = length(ParamNames)
    K = length(MetNames)
    n_data = length(fp.log_abs_rates)

    # Build params NamedTuple from log-space x
    fitted = NamedTuple{ParamNames}(ntuple(i -> exp(x[i]), Val(N)))
    params = merge(fitted, (Keq = fp.Keq, E_total = 1.0))

    # Pass 1: fill log_ratios_buffer and count sign mismatches. A zero or NaN
    # prediction counts as one: its sign never equals a nonzero rate's ±1.
    n_mismatch = 0
    @inbounds for i in 1:n_data
        concs = NamedTuple{MetNames}(ntuple(
            j -> getproperty(fp.data, MetNames[j])[i], Val(K),
        ))
        pred = Derivation.rate_equation(fp.mechanism, concs, params)
        meas_sign = sign(fp.data.Rate[i])
        if sign(pred) != meas_sign
            buf[i] = 10.0
            n_mismatch += 1
        else
            buf[i] = log(abs(pred)) - fp.log_abs_rates[i]
        end
    end

    # Pass 2: loss. Relative (fp.absolute false) → per-group mean-centered, removing
    # each group's arbitrary scale. Absolute → uncentered: the y-axis is absolute
    # per-enzyme turnover, so the absolute magnitude is meaningful.
    total_loss = 0.0
    if fp.absolute
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

In relative mode (`fp.absolute == false`) the returned steady-state rate constants are
multiplied by one factor so the mechanism reproduces the reference turnovers
`fp.reference` (see `FittingProblem`): the geometric mean of measured over predicted
turnover, as `rescale_parameter_values` computes it. A fit whose predicted turnover at a
reference anchor is zero, non-finite or of the other sign has no such factor: it returns
its parameters unscaled, with `loss = Inf` and `retcode = :ReferenceFailed`, and does not
throw. In absolute mode the rate constants are returned as fitted.

Returns a NamedTuple `(params, loss, retcode)` where:
- `params`: fitted rate constants as a NamedTuple
- `loss`: the best loss value achieved
- `retcode`: the `Symbol` form of the best restart's `sol.retcode`. Only
  `:Success` indicates the optimizer converged on its own criteria; any other
  value (e.g. `:MaxTime` — hit the time budget; `:Failure`; `:Default` — ran but
  did not flag success) means the fit should be treated as un-converged (check
  `retcode !== :Success`); `:ReferenceFailed` marks a failed reference scale.
"""
function fit_rate_equation(@nospecialize(fp::FittingProblem), optimizer;
    n_restarts::Int=20,
    maxtime::Real=60.0,
    maxiters::Integer=10_000_000,
    abstol::Union{Real,Nothing}=nothing,
    reltol::Union{Real,Nothing}=nothing,
    callback=nothing,
    lb=fill(-15.0, length(Derivation.fitted_params(fp.mechanism))),
    ub=fill(15.0, length(Derivation.fitted_params(fp.mechanism))),
    solver_kwargs=(;),
)
    # `fp` is @nospecialize, so `fitted_params(fp.mechanism)` infers as `Any`; asserting a
    # tuple of Symbols lets `np` infer as an `Int` and the NamedTuple built from the names
    # as a `NamedTuple`.
    pnames = Derivation.fitted_params(fp.mechanism)::Tuple{Vararg{Symbol}}
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

    for _ in 1:n_restarts
        x0 = clamp.(randn(np) .* 2.0, lb, ub)
        prob = _optimization_problem(fp, x0, lb, ub)
        sol = Optimization.solve(prob, optimizer; solve_kwargs...)
        if sol.objective < best_loss
            best_loss = sol.objective
            best_x .= sol.u
            best_retcode = Symbol(sol.retcode)
        end
    end

    result_params = NamedTuple{pnames}(ntuple(i -> exp(best_x[i]), Val(np)))
    fp.absolute && return (params = result_params, loss = best_loss, retcode = best_retcode)
    ratios = _anchor_log_ratios(fp.mechanism, merge(result_params, (Keq = fp.Keq,)),
                                fp.reference)
    ratios === nothing &&
        return (params = result_params, loss = Inf, retcode = :ReferenceFailed)
    scale = exp(sum(ratios) / length(ratios))
    return (params = Derivation.scale_rate_constants(fp.mechanism, result_params, scale),
            loss = best_loss, retcode = best_retcode)
end

"""
    rescale_parameter_values(m, params, kcats) -> (params, spread)

Re-express the parameters `params` of mechanism `m` on another rate scale by multiplying
every steady-state rate constant (the `k_…` constants) by one factor. The binding
constants `K_…`, `Keq`, `E_total`, `L` and the regulatory constants stay as they are, and
every rate changes by that factor. `params` must include `Keq`, which the dependent rate
constants need; merge it into a fit's `params`.

The factor matches measured turnovers `kcats` — a kcat, or a specific activity in turnover
units — at stated concentrations: a vector of `concentrations => turnover` pairs, the
concentrations a `NamedTuple` per anchor; metabolites an anchor leaves out are zero.
Anchor concentrations are in the data's concentration units. Turnover is per unit
`E_total`, the active-site concentration, so divide a kcat quoted per oligomer by its
number of active sites first. Several anchors combine by the geometric mean of measured
over predicted turnover, the least-squares scale in log space. `spread` is the standard
deviation of the anchors' `log(measured/predicted)`, 0 for one anchor; a large spread says
the anchors disagree with the model or with each other. An anchor where the model predicts
zero or non-finite turnover, or turnover of the other sign, is an error. A relative fit
comes back already scaled to its problem's `reference` this way (see `FittingProblem`);
rescaling it re-expresses it on other anchors.

```julia
fit = fit_rate_equation(
    FittingProblem(m, data; Keq = 2.0, reference = [(A = 10.0, B = 5.0) => 1.0]),
    optimizer)
kcats = [(A = 10.0, B = 5.0) => 42.0, (A = 20.0, B = 5.0) => 44.0]
params, spread = rescale_parameter_values(m, merge(fit.params, (; Keq = 2.0)), kcats)
```
"""
function rescale_parameter_values(@nospecialize(m::Mechanisms.AbstractEnzymeMechanism),
                                  params::NamedTuple, kcats::AbstractVector{<:Pair})
    _check_has_keq(params)
    _check_anchors(kcats, Mechanisms.metabolites(m), "kcats", "mechanism")
    log_ratios = _anchor_log_ratios(m, params, kcats)
    log_ratios === nothing && error(
        "kcats: the mechanism predicts zero, non-finite or opposite-sign turnover at an " *
        "anchor; an anchor needs a finite, nonzero predicted turnover of the measured " *
        "turnover's sign")
    center = sum(log_ratios) / length(log_ratios)
    spread = length(log_ratios) == 1 ? 0.0 :
        sqrt(sum(r -> (r - center)^2, log_ratios) / (length(log_ratios) - 1))
    (Derivation.scale_rate_constants(m, params, exp(center)), spread)
end

"""
    _anchor_log_ratios(m, params, anchors) -> Union{Vector{Float64}, Nothing}

`log(measured / predicted)` for each anchor `concentrations => measured` of `anchors`,
where `predicted` is the turnover of `m` at `params` (which hold `Keq`) per unit
`E_total`, with every metabolite the anchor leaves out at zero; `nothing` when a
predicted turnover is zero, non-finite or of the other sign than the measured one, which
no factor on the rate constants reproduces.
"""
function _anchor_log_ratios(@nospecialize(m::Mechanisms.AbstractEnzymeMechanism),
                            params::NamedTuple, anchors)
    mets = Mechanisms.metabolites(m)
    turnover_params = merge(params, (E_total = 1.0,))
    log_ratios = Float64[]
    for (concs, measured) in anchors
        predicted = Derivation.rate_equation(m,
            NamedTuple{mets}(Tuple(Float64(get(concs, k, 0.0)) for k in mets)),
            turnover_params)
        isfinite(predicted) && predicted != 0 && sign(predicted) == sign(measured) ||
            return nothing
        push!(log_ratios, log(measured / predicted))
    end
    log_ratios
end

"""Errors unless `params` holds `Keq`, which the dependent rate constants need."""
_check_has_keq(params::NamedTuple) = haskey(params, :Keq) || error(
    "params must include Keq, which the dependent rate constants need; merge it into a " *
    "fit's params")
