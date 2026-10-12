# ABOUTME: Fitting rate equations to experimental data and scaling them to peak turnover.
# ABOUTME: Constructs FittingProblem, computes loss, runs multi-start fits, rescales fits.

"""
    FittingProblem{M, D}

Holds pre-processed experimental data and mechanism info for fitting
rate constants.

# Fields
- `mechanism`: an `AbstractEnzymeMechanism` instance
- `data`: `NamedTuple` of column vectors (via `Tables.columntable`)
- `group_point_indexes`: row indices grouped by unique `group` values
- `Keq`: fixed equilibrium constant
- `absolute`: `false` selects relative mode (per-group-centered loss, fit scaled to peak
  turnover 1); `true` selects absolute per-enzyme turnover mode (uncentered loss, no
  scale)
- `peak_range`: the concentration range `(lo, hi)` of the peak-turnover scale
- `log_abs_rates`: pre-computed `log.(abs.(Rate))`
- `log_ratios_buffer`: pre-allocated working buffer for loss computation
"""
struct FittingProblem{M<:AbstractEnzymeMechanism, D<:NamedTuple}
    mechanism::M
    data::D
    group_point_indexes::Vector{Vector{Int}}
    Keq::Float64
    absolute::Bool
    peak_range::Tuple{Float64,Float64}
    log_abs_rates::Vector{Float64}
    log_ratios_buffer::Vector{Float64}
end

"""
    FittingProblem(mechanism::AbstractEnzymeMechanism, table;
        Keq::Real, absolute::Bool = false, peak_range = (1e-12, 1e3))

Construct a `FittingProblem` from an enzyme mechanism and tabular data.

The table must have columns: `group`, `Rate`, and one column per
metabolite matching `metabolites(mechanism)`, so no metabolite may be
named `group` or `Rate`. Uses `Tables.columntable` for conversion.

With `absolute = false` (the default) the data are relative: the loss centers each
group's log-ratios, and `fit_rate_equation` scales the fitted rate constants so the
mechanism's peak turnover — its largest net forward rate per enzyme over every
metabolite's concentration within `peak_range`, products and regulators included — is 1.
With `absolute = true` the data are absolute per-enzyme turnover: the loss is uncentered
and the fit is returned unscaled. `peak_range = (lo, hi)` is in the data's concentration
units; the default covers physiological concentrations in M or mM, and data in µM need a
higher upper bound.

Rate values must be finite and nonzero (the loss takes their log; a zero rate gives
`-Inf`), every concentration must be a finite number ≥ 0 (zero is valid), `Keq`
must be positive and finite, and `peak_range` must be a tuple `(lo, hi)` with
`0 < lo < hi < Inf`.
"""
function FittingProblem(mechanism::AbstractEnzymeMechanism, table;
        Keq::Real, absolute::Bool=false, peak_range=DEFAULT_PEAK_RANGE)
    data = _rate_table(table, metabolites(mechanism), Keq)
    group_map = Dict{eltype(data.group), Vector{Int}}()
    for (i, g) in enumerate(data.group)
        push!(get!(() -> Int[], group_map, g), i)
    end
    # The typed constructor converts Keq to the field type.
    FittingProblem{typeof(mechanism), typeof(data)}(
        mechanism, data, collect(values(group_map)), Keq, absolute,
        _peak_range(peak_range), log.(abs.(data.Rate)),
        Vector{Float64}(undef, length(data.Rate)))
end

"""
    _rate_table(table, mnames, Keq) → NamedTuple

Validate a rate table and return it as `Tables.columntable(table)`: `Keq` must be
positive and finite, the table needs a `group` column, a `Rate` column and one column per
name in `mnames`, no name in `mnames` may be `group` or `Rate`, every concentration must
be a finite number ≥ 0 (zero is valid), and every rate must be a finite, nonzero number
(the loss takes its log).
"""
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

"""
Accept the concrete working-representation mechanism: compile to the
singleton once at construction so `loss!`'s hot path operates on the
@generated `EnzymeMechanism` / `AllostericEnzymeMechanism` (0-alloc).
"""
FittingProblem(mechanism::Union{Mechanism, AllostericMechanism}, table; kwargs...) =
    FittingProblem(compile_mechanism(mechanism), table; kwargs...)

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
    ParamNames = fitted_params(fp.mechanism)
    MetNames = metabolites(fp.mechanism)
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
        pred = rate_equation(fp.mechanism, concs, params)
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

In relative mode (`fp.absolute == false`) the returned rate constants are scaled so the
mechanism's peak turnover within `fp.peak_range` is 1 (see `FittingProblem`); a law with
no positive turnover in the range is returned unscaled. In absolute mode they are returned
as fitted.

Returns a NamedTuple `(params, loss, retcode, peak_at_range_edge)` where:
- `params`: fitted rate constants as a NamedTuple
- `loss`: the best loss value achieved
- `retcode`: the `Symbol` form of the best restart's `sol.retcode`. Only
  `:Success` indicates the optimizer converged on its own criteria; any other
  value (e.g. `:MaxTime` — hit the time budget; `:Failure`; `:Default` — ran but
  did not flag success) means the fit should be treated as un-converged (check
  `retcode !== :Success`).
- `peak_at_range_edge`: `true` when the scaled law still rises at the upper edge of
  `fp.peak_range` — by more than 10% over the range's last decade — so a wider range would
  change the scale; always `false` in absolute mode
"""
function fit_rate_equation(@nospecialize(fp::FittingProblem), optimizer;
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
    # `fp` is @nospecialize, so `fitted_params(fp.mechanism)` infers as `Any`; asserting a
    # tuple of Symbols lets `np` infer as an `Int` and the NamedTuple built from the names
    # as a `NamedTuple`.
    pnames = fitted_params(fp.mechanism)::Tuple{Vararg{Symbol}}
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
    peak_at_range_edge = false
    if !fp.absolute
        peak = _peak_turnover(fp.mechanism, merge(result_params, (Keq = fp.Keq,)),
                              fp.peak_range)
        peak_at_range_edge = peak.at_edge
        0 < peak.value < Inf && (result_params =
            _scale_rate_constants(fp.mechanism, result_params, 1 / peak.value))
    end
    return (params = result_params, loss = best_loss, retcode = best_retcode,
            peak_at_range_edge = peak_at_range_edge)
end

"""
`params` with every steady-state rate constant of `m` multiplied by `s`, every other
entry unchanged. A common factor on those constants multiplies every rate by the same
factor and leaves the binding constants, `Keq`, `E_total`, `L` and the regulatory
constants alone.
"""
function _scale_rate_constants(@nospecialize(m::AbstractEnzymeMechanism),
                               params::NamedTuple, s::Real)
    rate_constants = _ss_rate_constant_names(m)
    NamedTuple{keys(params)}(map((k, v) -> k in rate_constants ? v * s : v,
                                 keys(params), values(params)))
end

"""
The default `peak_range`: the concentrations, in the data's units, over which the peak
turnover that scales a relative fit is sought. It covers physiological concentrations in M
or in mM.
"""
const DEFAULT_PEAK_RANGE = (1e-12, 1e3)

"""
`peak_range` as a `Tuple{Float64,Float64}`, after checking that it is a tuple `(lo, hi)`
with `0 < lo < hi < Inf`.
"""
function _peak_range(peak_range)
    peak_range isa Tuple{Real,Real} && 0 < peak_range[1] < peak_range[2] < Inf ||
        error("peak_range must be a tuple (lo, hi) with 0 < lo < hi < Inf; got " *
              repr(peak_range))
    (Float64(peak_range[1]), Float64(peak_range[2]))
end

"""
The turnover `v/E_total` of `mechanism` at `params` (whose `E_total` is 1) as a function
of the log10 concentrations of `metabolites(mechanism)`, in that order; `N` counts them.
"""
struct LogConcentrationTurnover{M<:AbstractEnzymeMechanism, P<:NamedTuple, N}
    mechanism::M
    params::P
end

LogConcentrationTurnover(mechanism::AbstractEnzymeMechanism, params::NamedTuple) =
    LogConcentrationTurnover{typeof(mechanism), typeof(params),
                             length(metabolites(mechanism))}(mechanism, params)

(t::LogConcentrationTurnover{M, P, N})(u::Vector{Float64}) where {M, P, N} =
    rate_equation(t.mechanism,
        NamedTuple{metabolites(t.mechanism)}(ntuple(i -> exp10(u[i]), Val(N))), t.params)

"""Decades between the points of one coordinate's scan in the peak search."""
const PEAK_SCAN_STEP = 0.25

"""`v`, or `-Inf` when `v` is NaN or infinite, so the peak search never picks it."""
_finite_or_minus_inf(v) = isfinite(v) ? v : -Inf

"""
The largest value of `f` that a heuristic search finds over the box `[a, b]^n` of log10
concentrations, the point that reaches it, and whether `f` still rises at the box's upper
edge. Points where `f` is NaN or infinite never win.

The search is deterministic. It evaluates every corner of the box (the two extreme
corners only, above 12 coordinates) and its center, then climbs from the four best of
them one coordinate at a time (`_coordinate_ascent!`); it can stop at a local maximum.
`f` still rises at the upper edge when some coordinate of the best point sits at `b` and
`f` there exceeds 1.1 times its value with that coordinate lowered by one decade (to `a`
at most); a best value that is not positive never counts as rising.
"""
function _maximize_on_box(f::FunctionWrapper{Float64, Tuple{Vector{Float64}}}, n::Int,
                          a::Float64, b::Float64)
    corners = n <= 12 ? (0:(2^n - 1)) : (0, 2^n - 1)
    starts = [Float64[(k >> (j - 1)) & 1 == 1 ? b : a for j in 1:n] for k in corners]
    push!(starts, fill((a + b) / 2, n))
    values = [_finite_or_minus_inf(f(u)) for u in starts]
    order = sort(eachindex(starts); by = i -> (-values[i], i))
    best_u, best_v = starts[order[1]], values[order[1]]
    for i in order[1:min(4, length(order))]
        u = copy(starts[i])
        v = _coordinate_ascent!(f, u, values[i], a, b)
        v > best_v && ((best_u, best_v) = (u, v))
    end
    at_edge = false
    if best_v > 0
        for j in 1:n
            best_u[j] >= b - 1e-6 || continue
            lowered = copy(best_u)
            lowered[j] = max(a, b - 1)
            if best_v > 1.1 * f(lowered)
                at_edge = true
                break
            end
        end
    end
    best_u, best_v, at_edge
end

"""
Climb `f` from `u`, where its value is `v`, one coordinate at a time (`_line_maximum!`),
in place, until a full sweep gains less than one part in 10¹²; returns the final value.
"""
function _coordinate_ascent!(f::FunctionWrapper{Float64, Tuple{Vector{Float64}}},
                             u::Vector{Float64}, v::Float64, a::Float64, b::Float64)
    for _ in 1:100
        v_start = v
        for j in eachindex(u)
            v = _line_maximum!(f, u, j, v, a, b)
        end
        v > v_start && v - v_start > 1e-12 * abs(v) || break
    end
    v
end

"""
Move `u[j]` to the best point of `f` along coordinate `j` within `[a, b]`, in place, and
return the value there; `v` is the value at the current `u`. Scans the coordinate every
`PEAK_SCAN_STEP` decades, then refines the best scan point by golden-section search
between its neighbours. The current point stays unless a candidate beats it.
"""
function _line_maximum!(f::FunctionWrapper{Float64, Tuple{Vector{Float64}}},
                        u::Vector{Float64}, j::Int, v::Float64, a::Float64, b::Float64)
    n_steps = max(1, ceil(Int, (b - a) / PEAK_SCAN_STEP))
    h = (b - a) / n_steps
    best_t, best_v = u[j], v
    for k in 0:n_steps
        u[j] = k == n_steps ? b : a + k * h
        w = _finite_or_minus_inf(f(u))
        w > best_v && ((best_t, best_v) = (u[j], w))
    end
    lo, hi = max(a, best_t - h), min(b, best_t + h)
    φ = (sqrt(5) - 1) / 2
    x1, x2 = hi - φ * (hi - lo), lo + φ * (hi - lo)
    u[j] = x1; f1 = _finite_or_minus_inf(f(u))
    u[j] = x2; f2 = _finite_or_minus_inf(f(u))
    while hi - lo > 1e-9
        if f1 >= f2
            hi, x2, f2 = x2, x1, f1
            x1 = hi - φ * (hi - lo)
            u[j] = x1; f1 = _finite_or_minus_inf(f(u))
        else
            lo, x1, f1 = x1, x2, f2
            x2 = lo + φ * (hi - lo)
            u[j] = x2; f2 = _finite_or_minus_inf(f(u))
        end
    end
    f1 > best_v && ((best_t, best_v) = (x1, f1))
    f2 > best_v && ((best_t, best_v) = (x2, f2))
    u[j] = best_t
    best_v
end

"""
    _peak_turnover(m, params, peak_range) -> (; value, concs, at_edge)

The peak turnover of mechanism `m` at `params`, which must hold `Keq`: the largest net
forward rate per enzyme, `v/E_total`, over every metabolite's concentration within
`peak_range` — substrates, products, inhibitors and regulators alike
(`_maximize_on_box`). `concs` is the point that reaches it, and `at_edge` says whether
the turnover still rises at the range's upper edge. The `E_total` of `params` is
ignored. `value` is not positive when the mechanism has no forward turnover in the range.
The turnover sits behind a `FunctionWrapper`, so the search compiles once, shared by all
mechanisms.
"""
function _peak_turnover(@nospecialize(m::AbstractEnzymeMechanism),
                        @nospecialize(params::NamedTuple),
                        peak_range::Tuple{Float64, Float64})
    mets = metabolites(m)
    f = FunctionWrapper{Float64, Tuple{Vector{Float64}}}(
        LogConcentrationTurnover(m, merge(params, (E_total = 1.0,))))
    u, value, at_edge = _maximize_on_box(f, length(mets), log10(peak_range[1]),
                                         log10(peak_range[2]))
    (value = value, concs = NamedTuple{mets}(Tuple(exp10.(u))), at_edge = at_edge)
end

"""
    rescale_parameter_values(m, params; peak_range = (1e-12, 1e3)) -> (params, nothing)
    rescale_parameter_values(m, params, kcats) -> (params, spread)

Re-express the parameters `params` of mechanism `m` on another rate scale by multiplying
every steady-state rate constant (the `k_…` constants) by one factor. The binding
constants `K_…`, `Keq`, `E_total`, `L` and the regulatory constants stay as they are, and
every rate changes by that factor. `params` must include `Keq`, which the dependent rate
constants need; merge it into a fit's `params`.

With `kcats`, a vector of `concentrations => turnover` pairs, the factor matches measured
turnovers — a kcat, or a specific activity in turnover units — at stated concentrations,
given as a `NamedTuple` per anchor; metabolites an anchor leaves out are zero. Anchor
concentrations are in the data's concentration units. Turnover is per unit `E_total`, the
active-site concentration, so divide a kcat quoted per oligomer by its number of active
sites first. Several anchors combine by the geometric mean of measured over predicted
turnover, the least-squares scale in log space. `spread` is the standard deviation of the
anchors' `log(measured/predicted)`, 0 for one anchor; a large spread says the anchors
disagree with the model or with each other. An anchor where the model predicts zero or
non-finite turnover, or turnover of the other sign, is an error.

Without `kcats`, the factor makes the mechanism's peak turnover — its largest net forward
rate per enzyme over every metabolite's concentration within `peak_range`, in the data's
concentration units — equal 1, the scale `fit_rate_equation` gives relative fits, and
`spread` is `nothing`. The default range covers physiological concentrations in M or mM;
data in µM need a higher upper bound.

```julia
fit = fit_rate_equation(FittingProblem(m, data; Keq = 2.0), optimizer)
kcats = [(A = 10.0, B = 5.0) => 42.0, (A = 20.0, B = 5.0) => 44.0]
params, spread = rescale_parameter_values(m, merge(fit.params, (; Keq = 2.0)), kcats)
```
"""
function rescale_parameter_values(@nospecialize(m::AbstractEnzymeMechanism),
                                  params::NamedTuple; peak_range=DEFAULT_PEAK_RANGE)
    _check_has_keq(params)
    range = _peak_range(peak_range)
    peak = _peak_turnover(m, params, range).value
    0 < peak < Inf || error("the mechanism has no positive forward turnover within " *
                            "peak_range $range; the peak scale is undefined")
    (_scale_rate_constants(m, params, 1 / peak), nothing)
end

function rescale_parameter_values(@nospecialize(m::AbstractEnzymeMechanism),
                                  params::NamedTuple, kcats::AbstractVector{<:Pair})
    _check_has_keq(params)
    isempty(kcats) &&
        error("kcats is empty; pass at least one (concentrations => turnover) anchor")
    mets = metabolites(m)
    turnover_params = merge(params, (E_total = 1.0,))
    log_ratios = map(kcats) do (concs, measured)
        concs isa NamedTuple ||
            error("kcats: concentrations must be a NamedTuple; got $(repr(concs))")
        for (k, c) in pairs(concs)
            k in mets || error("kcats: $k is not a metabolite of the mechanism " *
                               "(metabolites: $(join(mets, ", ")))")
            c isa Real && isfinite(c) && c >= 0 || error(
                "kcats: concentration of $k must be a finite number ≥ 0; got $(repr(c))")
        end
        measured isa Real && isfinite(measured) && measured != 0 || error(
            "kcats: turnover must be a finite, nonzero number; got $(repr(measured))")
        predicted = rate_equation(m,
            NamedTuple{mets}(Tuple(Float64(get(concs, k, 0.0)) for k in mets)),
            turnover_params)
        isfinite(predicted) && predicted != 0 && sign(predicted) == sign(measured) || error(
            "kcats: at $concs the mechanism predicts turnover $predicted; an anchor " *
            "needs a finite, nonzero predicted turnover of the measured turnover's sign")
        log(measured / predicted)
    end
    center = sum(log_ratios) / length(log_ratios)
    spread = length(log_ratios) == 1 ? 0.0 :
        sqrt(sum(r -> (r - center)^2, log_ratios) / (length(log_ratios) - 1))
    (_scale_rate_constants(m, params, exp(center)), spread)
end

"""Errors unless `params` holds `Keq`, which the dependent rate constants need."""
_check_has_keq(params::NamedTuple) = haskey(params, :Keq) || error(
    "params must include Keq, which the dependent rate constants need; merge it into a " *
    "fit's params")
