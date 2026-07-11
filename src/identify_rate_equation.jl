# ABOUTME: Beam-search pipeline to identify the best rate equation.
# ABOUTME: Enumerates mechanisms, fits each, selects via CV; a comment-
# ABOUTME: stripped rate-equation string key dedups equivalent equations.

using DataFrames
using CSV
using Statistics

"""
    IdentifyRateEquationProblem{R, D}

Holds the reaction, experimental data, and equilibrium
constant for rate equation identification.

# Fields
- `reaction`: the `EnzymeReaction` instance
- `data`: `NamedTuple` of column vectors with `:group`,
  `:Rate`, and metabolite columns
- `Keq`: fixed equilibrium constant
- `scale_k_to_kcat`: target kcat for SS-rate rescaling
  before fitting (`nothing` = no rescaling, positive Float64 = target)
"""
struct IdentifyRateEquationProblem{
    R<:EnzymeReaction, D<:NamedTuple
}
    reaction::R
    data::D
    Keq::Float64
    scale_k_to_kcat::Union{Float64,Nothing}
end

function IdentifyRateEquationProblem(
    reaction::EnzymeReaction, table; Keq::Real,
    scale_k_to_kcat::Union{Real,Nothing}=1.0
)
    scale_k_to_kcat !== nothing && scale_k_to_kcat <= 0 && error(
        "scale_k_to_kcat must be positive (or nothing); got $scale_k_to_kcat")
    data = Tables.columntable(table)
    col_names = keys(data)

    for req in (:group, :Rate)
        req in col_names ||
            error("Missing required column: $req")
    end
    # Extract metabolite names from reaction (struct accessors return
    # Substrate/Product/RegulatorMults; pull the underlying Symbol via
    # `name()`).
    mnames = tuple(
        (name(s) for s in substrates(reaction))...,
        (name(p) for p in products(reaction))...,
        (name(regulator(r)) for r in regulators(reaction))...,
    )
    for m in mnames
        m in col_names ||
            error(
                "Missing metabolite column: $m")
    end

    # Validate non-zero rates
    for i in eachindex(data.Rate)
        data.Rate[i] == 0 &&
            error(
                "Zero rate at row $i: " *
                "log(0) is undefined")
    end

    # Validate at least 2 groups for CV
    n_groups = length(unique(data.group))
    n_groups >= 2 || error(
        "Need at least 2 unique groups for " *
        "cross-validation, got $n_groups")

    sk = scale_k_to_kcat === nothing ? nothing : Float64(scale_k_to_kcat)
    IdentifyRateEquationProblem{
        typeof(reaction),typeof(data)
    }(reaction, data, Float64(Keq), sk)
end

"""
    IdentifyRateEquationResults

Results from `identify_rate_equation`.

# Fields
- `best`: the best `AbstractEnzymeMechanism`
  (lowest loss at optimal param count)
- `cv_results`: `DataFrame` with LOOCV results for
  top candidates per param count
"""
struct IdentifyRateEquationResults
    best::AbstractEnzymeMechanism
    cv_results::DataFrame
end

"""
    identify_rate_equation(prob; optimizer,
        min_beam_width=50, loss_rel_threshold=2.0, loss_abs_threshold=0.01,
        loss_parsimony_threshold=1.01,
        max_param_count=20, n_restarts=20, maxtime=60.0, maxiters=10_000_000,
        abstol=nothing, reltol=nothing, callback=nothing, solver_kwargs=(;),
        n_cv_candidates=5, se_threshold=1.0, perm_p_threshold=0.16,
        save_dir=_default_save_dir(), show_progress=true)

Find the best rate equation for the given reaction
and data using beam search.

# Keyword Arguments
- `min_beam_width::Int = 50`: minimum mechanisms
  to keep per param-count tier
- `loss_rel_threshold::Float64 = 2.0`: relative tolerance
  for beam selection (see "Beam selection" below)
- `loss_abs_threshold::Float64 = 0.01`: absolute tolerance
  for beam selection
- `loss_parsimony_threshold::Float64 = 1.01`: a mechanism
  keeps expanding only if its loss is within this factor of
  the best model of any smaller parameter count — an added
  parameter must earn its keep. Combined with the other loss
  thresholds via `min`; `min_beam_width` is a cumulative
  per-count budget (see "Beam selection" below). `Inf` disables it.
- `max_param_count::Int = 20`: stop expanding beyond
- `eq_complexity_filter::Int = 337`: skip (before fitting, before any
  derivation) any mechanism whose rate equation is more complex than this —
  roughly the number of terms in the rate equation's denominator: V×τ, the
  number of RE-segments times the spanning-tree count of the catalytic segment
  graph, i.e. the products the compiled equation evaluates on every call. The
  default 337 is one above a fully steady-state random-order bi-bi (V×τ = 336),
  so such a mechanism passes and anything more complex is skipped — equations
  denser than that are impractical to fit and can blow up derivation/codegen.
  Computed from the mechanism graph alone.
- `optimizer`: Optimization.jl optimizer (required).
  Recommended: `CMAEvolutionStrategyOpt()` from OptimizationCMAEvolutionStrategy.
- `n_restarts::Int = 20`: multi-start restarts per fit
- `maxtime::Real = 600.0`: max time per fit (seconds; common solver
  option, forwarded to `Optimization.solve`)
- `maxiters::Integer = 10_000_000`: max iterations per optimizer run
  (common solver option, forwarded to `Optimization.solve`)
- `abstol`/`reltol`/`callback = nothing`: Optimization.jl common solver
  options, forwarded to `Optimization.solve` only when set
- `solver_kwargs::NamedTuple = (;)`: solver-specific options forwarded
  verbatim to `Optimization.solve` (e.g. `(; popsize=200)` for a CMA-ES
  solver that supports it); the caller matches its contents to `optimizer`
- `n_cv_candidates::Int = 5`: LOOCV top N
  **unique-rate-equation** candidates per param count
- `se_threshold::Float64 = 1.0`: paired 1-SE multiplier for
  model selection. Simpler-model bucket accepted iff its mean
  paired log-loss difference vs the best bucket is `≤
  se_threshold * std(diffs)/sqrt(n_folds)`. Default 1.0 is the
  textbook "1-SE rule".
- `perm_p_threshold::Float64 = 0.16`: minimum one-sided
  permutation p-value for model selection. Simpler-model
  bucket accepted iff `p > perm_p_threshold` under the
  sign-flip null. Default 0.16 matches paired 1-SE empirically.
- `save_dir::String = _default_save_dir()`: output directory for the
  search CSVs (`initial_mechanisms.csv` + `equation_search_iteration_N.csv`),
  plus `loocv_results.csv` (the LOOCV table for every cross-validated
  candidate) and `best_equation.csv` (the selected equation and its
  fitted parameters)

# Beam selection

A mechanism at parameter count `n` qualifies for the next-level
beam if either:
- its loss ≤ `min(loss_rel_threshold * best(n) + loss_abs_threshold,
  loss_parsimony_threshold * best(<n))`, where `best(n)` is the
  lowest loss at parameter count `n` and `best(<n)` is the lowest
  loss over all smaller counts; the second term is dropped at the
  base count (no smaller level exists yet),
- OR it falls under the `min_beam_width` budget: a cumulative
  per-count allowance that expands at least `min_beam_width`
  mechanisms at each count over the whole search, then stops —
  spent once, not re-granted each time the count is revisited.

The additive term protects against `best_loss` approaching zero
(simulated / very-low-loss data) where a purely multiplicative
threshold would collapse the beam to the single best mechanism.

# Model selection (LOOCV)

For each `n_params` bucket below `n_min` (lowest mean
log-fold-loss) the rule computes paired log-loss differences
between the bucket's representative and `n_min`'s, then
accepts the simpler bucket iff BOTH:

1. **Paired 1-SE rule**: `mean(diffs) ≤ se_threshold *
   std(diffs)/sqrt(n_folds)`.
2. **One-sided permutation test**:
   `permutation_p > perm_p_threshold`, where p is computed by
   exact enumeration when `n_folds ≤ 20` and Monte Carlo
   (10⁶ samples) otherwise.

Returns the smallest passing `n_params`; falls through to
`n_min` if none pass. Within the chosen bucket the mechanism
with lowest training loss wins. Per-bucket representative =
the row with the lowest `cv_score` in that bucket. Diagnostic
columns `mean_log_loss_diff`, `se_paired`, `permutation_p`
are surfaced in `cv_results`; the `n_min` bucket has 0.0 in
all three.
"""
function identify_rate_equation(
    prob::IdentifyRateEquationProblem;
    # Beam search
    min_beam_width::Int = 50,
    loss_rel_threshold::Float64 = 2.0,
    loss_abs_threshold::Float64 = 0.01,
    loss_parsimony_threshold::Float64 = 1.01,
    max_param_count::Int = 20,
    eq_complexity_filter::Int = 337,
    # Fitting
    optimizer,
    n_restarts::Int = 20,
    maxtime::Real = 600.0,
    maxiters::Integer = 10_000_000,
    abstol::Union{Real,Nothing} = nothing,
    reltol::Union{Real,Nothing} = nothing,
    callback = nothing,
    solver_kwargs = (;),
    # Model selection
    n_cv_candidates::Int = 5,
    se_threshold::Float64 = 1.0,
    perm_p_threshold::Float64 = 0.16,
    # Output & parallelism
    save_dir::String = _default_save_dir(),
    show_progress::Bool = true,
)
    fitting_kwargs = (;
        n_restarts, maxtime, maxiters,
        abstol, reltol, callback, solver_kwargs)

    if isdir(save_dir)
        existing = filter(
            f -> endswith(f, ".csv"),
            readdir(save_dir))
        isempty(existing) || error(
            "save_dir already contains CSV " *
            "files. Use an empty directory " *
            "to avoid mixing results.")
    end

    mechanisms, df = _beam_search(prob;
        min_beam_width, loss_rel_threshold,
        loss_abs_threshold, loss_parsimony_threshold,
        max_param_count, eq_complexity_filter, save_dir, show_progress,
        optimizer, n_cv_candidates,
        fitting_kwargs...)

    result = _cv_model_selection(
        mechanisms, df, prob;
        n_cv_candidates, se_threshold, perm_p_threshold,
        optimizer, save_dir, show_progress,
        fitting_kwargs...)
    _progress(save_dir, show_progress, "Done. Results saved to $save_dir")
    return result
end

# Write result rows to `<save_dir>/<filename>`, creating `save_dir` if absent.
function _write_rows_csv(save_dir::String, filename::String, rows)
    isdir(save_dir) || mkpath(save_dir)
    CSV.write(joinpath(save_dir, filename), _rows_to_dataframe(rows))
end

"""Save the base-tier fit (all init mechanisms) to `initial_mechanisms.csv`."""
_save_initial_csv(save_dir::String, rows) =
    _write_rows_csv(save_dir, "initial_mechanisms.csv", rows)

"""
Save one expansion iteration to `equation_search_iteration_<iteration>.csv`.
`iteration` is a 1-based sequential counter, NOT a parameter count — the
real fitted count is the `n_params` column of each row.
"""
_save_iteration_csv(save_dir::String, rows, iteration::Int) =
    _write_rows_csv(
        save_dir, "equation_search_iteration_$(iteration).csv", rows)

"""
Convert result row NamedTuples to a DataFrame.
Row order is preserved (no sorting) to maintain
alignment with the mechanism vector.
"""
function _rows_to_dataframe(rows)
    isempty(rows) && return DataFrame()
    all_pnames = Set{Symbol}()
    for row in rows
        for p in row.fitted_param_names
            push!(all_pnames, p)
        end
    end
    sorted_pnames = sort(collect(all_pnames))

    df = DataFrame(
        n_params = [r.n_params for r in rows],
        parent_n_params =
            Union{Missing,Int}[get(r, :parent_n_params, missing) for r in rows],
        loss = [r.loss for r in rows],
        mechanism_type = [r.mechanism_type for r in rows],
        parent_mechanism_type = Union{Missing,String}[
            get(r, :parent_mechanism_type, missing) for r in rows],
        rate_equation = [r.rate_equation for r in rows],
        retcode = Union{Missing,String}[r.retcode for r in rows],
        error = Union{Missing,String}[r.error for r in rows],
        eq_hash = [r.eq_hash for r in rows],
        fit_inherited = Union{Missing,Bool}[r.fit_inherited for r in rows],
    )
    for pn in sorted_pnames
        df[!, pn] = [
            pn in r.fitted_param_names ?
                r.fitted_param_values[
                    findfirst(==(pn), r.fitted_param_names)] :
                missing
            for r in rows
        ]
    end
    df
end

"""
Equation-identity key: the rendered rate-equation string with provenance
removed — `# …` header lines and Wegscheider `(substituted into v)` lines
(the choice of which dependent K was eliminated is cosmetic; it is already
substituted into v). Two mechanisms with the same key compute the identical
rate function. Used as a CSV tag and the LOOCV distinct-equation key.

Exact-render identity only: it hashes the (provenance-stripped) equation string
verbatim, so it detects textually identical equations, not algebraic equivalence.
"""
function _rate_eq_dedup_key(eq_text::AbstractString)
    kept = Iterators.filter(split(eq_text, '\n')) do ln
        l = strip(ln)
        !startswith(l, "#") && !occursin(ANNOTATION_SUBSTITUTED, l)
    end
    hash(join(kept, '\n'))
end

"""
Return indices into `losses` for mechanisms that qualify for the
beam at this level. A mechanism qualifies if either:
  • its loss ≤ cutoff, where
    cutoff = min(loss_rel_threshold * best_loss + loss_abs_threshold,
                 parsimony_cutoff) and the parsimony term is dropped
    when `parsimony_cutoff === nothing`,
  • OR its rank (1-indexed by ascending loss) ≤ min_beam_width.

`parsimony_cutoff` (the loss-parsimony threshold times the best loss
over all smaller parameter counts) only tightens the loss cutoff.
`min_beam_width` here is the number kept by the width floor for this
call; `_select_count!` passes the remaining cumulative per-count
budget, not a fixed per-sweep floor.

Mechanisms with non-finite losses (`Inf`, `NaN`) are excluded
unconditionally — they represent failed or non-converging fits
that should not propagate to the next level.
"""
function _select_beam(
    losses::AbstractVector{<:Real};
    loss_rel_threshold::Float64,
    loss_abs_threshold::Float64,
    min_beam_width::Int,
    best_override::Union{Nothing,Float64}=nothing,
    parsimony_cutoff::Union{Nothing,Float64}=nothing,
)
    finite_idx = [i for i in eachindex(losses) if isfinite(losses[i])]
    isempty(finite_idx) && return Int[]

    perm = sort(finite_idx; by=i -> losses[i])
    best = best_override === nothing ? losses[perm[1]] : best_override
    cutoff = loss_rel_threshold * best + loss_abs_threshold
    parsimony_cutoff !== nothing && (cutoff = min(cutoff, parsimony_cutoff))
    selected = Int[]
    for (rank, idx) in enumerate(perm)
        if losses[idx] <= cutoff || rank <= min_beam_width
            push!(selected, idx)
        end
    end
    # Return indices in original (input) order so callers don't
    # rely on the by-loss sort order, which is a side-effect of
    # the rank computation rather than part of the contract.
    sort!(selected)
end

"""
Select this sweep's parents at one parameter count under a *cumulative* floor
budget, and advance the budget. `expanded[c]` tracks how many count-`c`
mechanisms the whole search has expanded so far; the width floor may add at most
`min_beam_width - expanded[c]` more. Once the budget is spent, only the loss
cutoff admits at that count. Returns the selected indices (input order).
"""
function _select_count!(
    expanded::Dict{Int,Int}, c::Int, losses::AbstractVector{<:Real};
    loss_rel_threshold::Float64, loss_abs_threshold::Float64,
    min_beam_width::Int,
    best_override::Union{Nothing,Float64}=nothing,
    parsimony_cutoff::Union{Nothing,Float64}=nothing,
)
    budget = max(0, min_beam_width - get(expanded, c, 0))
    sel = _select_beam(losses;
        loss_rel_threshold, loss_abs_threshold,
        min_beam_width=budget, best_override, parsimony_cutoff)
    expanded[c] = get(expanded, c, 0) + length(sel)
    sel
end

"""One fitted mechanism: its own params + retcode + eq_hash + the CSV row."""
struct BatchEntry
    mech::Union{Mechanism, AllostericMechanism}
    n_params::Int
    loss::Float64
    retcode::Symbol
    eq_hash::UInt64
    row::NamedTuple
end

"""A mechanism that failed to compile or fit, with the exception message
captured as a string. Kept for the CSV + summary."""
struct FitFailure
    mech::Union{Mechanism, AllostericMechanism}
    error::String
end

# Compact, CSV-safe rendering of a thrown exception: type + truncated message.
_exc_string(e) = first(sprint(showerror, e), 200)

# CSV row for a mechanism that threw. Same NamedTuple schema as a fitted row,
# with `missing` wherever the value is unavailable (compile/fit never produced it).
# `mechanism_type` is the round-trippable parametric `EnzymeMechanism{Sig}` string when
# the mechanism compiles; falls back to the bare concrete type name
# (`"EnzymeRates.Mechanism"` / `"EnzymeRates.AllostericMechanism"`) when compilation
# itself fails, so the row still identifies the mechanism family.
function _failure_row(f::FitFailure)
    (n_params = missing,
     parent_n_params = missing,
     loss = missing,
     mechanism_type = try
         string(typeof(compile_mechanism(f.mech)))
     catch
         string(typeof(f.mech))
     end,
     parent_mechanism_type = missing,
     rate_equation = missing,
     retcode = missing,
     error = f.error,
     fitted_param_names = (),
     fitted_param_values = (),
     eq_hash = missing,
     fit_inherited = missing)
end

"""
Emit one progress line to flushed stdout AND append it to
`<save_dir>/progress.log`. Gated by `show_progress`. The explicit flush makes
lines appear in a redirected cluster job log (otherwise stdout is block-buffered
when redirected and withholds output until the process exits).
"""
function _progress(save_dir::AbstractString, show_progress::Bool, msg::AbstractString)
    show_progress || return nothing
    println(msg)
    flush(stdout)
    isdir(save_dir) || mkpath(save_dir)
    open(joinpath(save_dir, "progress.log"), "a") do io
        println(io, msg)
    end
    nothing
end

"""
One-line batch summary reconciling every child mechanism into six buckets that
sum to the child count: `new fits` (a genuine optimizer run, `fit_inherited ==
false`), `inherited` (a memoized fit reused, `fit_inherited == true`), three
`skipped` buckets dropped before fitting — one for a structure already fit in
an earlier batch, one for exceeding `max_param_count`, one for exceeding
`eq_complexity_filter` — and `errored` (compile/render/fit threw). Trailing
`Success` / `non-Success retcode` percentages are over the fitted set (`new +
inherited`).
"""
function _batch_summary(
    entries::Vector{BatchEntry}, failures::Vector{FitFailure};
    n_param_skipped::Int, n_complexity_skipped::Int, n_fitted_skipped::Int,
    max_param_count::Int, eq_complexity_filter::Int)
    n_fit   = length(entries)
    n_err   = length(failures)
    n_new   = count(e -> !e.row.fit_inherited, entries)
    n_inh   = n_fit - n_new
    n_succ  = count(e -> e.retcode === :Success, entries)
    n_other = n_fit - n_succ
    pct(x)  = n_fit == 0 ? 0.0 : round(100 * x / n_fit; digits=1)
    string(n_new, " new fits + ", n_inh, " inherited + ",
           n_fitted_skipped, " skipped (already fit) + ",
           n_param_skipped, " skipped (>", max_param_count, " params) + ",
           n_complexity_skipped, " skipped (>", eq_complexity_filter,
           " complexity) + ", n_err, " errored | Success ", pct(n_succ),
           "% | non-Success retcode ", pct(n_other), "%")
end

"""
Render the running best loss per parameter count, ascending, marking the counts
that improved this iteration with `*`. Counts read from `best_loss_by_count`;
`improved` is the set of counts whose best strictly dropped (or first appeared)
this iteration.
"""
function _best_loss_line(best_loss_by_count::Dict{Int,Float64}, improved::Set{Int})
    parts = [string(c, ":", round(best_loss_by_count[c]; sigdigits=4),
                    c in improved ? "*" : "")
             for c in sort(collect(keys(best_loss_by_count)))]
    annotation = isempty(improved) ? "   (no improvement)" : "   (* improved)"
    string("best loss by n_params: ", join(parts, " "), annotation)
end

"""
Compile + cap-check + fit every mechanism in `mechs`, deduplicating by `eq_hash`
so each distinct rate equation is fit once. Returns `(entries, failures)`:
`entries::Vector{BatchEntry}` are the mechanisms with a usable fit (each keeping
its own row, `retcode`, and `eq_hash`); `failures::Vector{FitFailure}` are
mechanisms that threw at compile/render/fit — captured WITH the exception text,
never silently swallowed. A mechanism whose fitted-param count exceeds
`max_param_count` is dropped before fitting (a cap skip, not a failure). A
mechanism whose structural `hash` is already in `fitted` (produced in an earlier
batch) is dropped before PASS-1 (an already-fit skip, not a failure); `fitted`
gets every structure on first production regardless of outcome (so a capped or
errored repeat also counts as already-fit — the common case is a genuine refit).

Two passes. PASS-1 (`pmap`) compiles, caps, and renders every mechanism to get its
`eq_hash`. PASS-2 (`pmap`) fits ONE representative per `eq_hash` not already in
`memo`, in parallel across all workers; the fit is stored in `memo` and copied to
every mechanism sharing that `eq_hash` — same equation ⟹ same fitted params and
same rescaling, so no refit is needed. `memo` persists across iterations (threaded
from `_beam_search`), so an equation fit in an earlier iteration is never refit. A
representative whose fit throws fails ALL of its duplicates (all-or-nothing per
equation). `fit_inherited` is `false` for the representative actually fit this
batch, `true` for every reused row. `mechs` is already structurally deduped by the
caller (`unique!`).
"""
function _process_batch(
    mechs, prob::IdentifyRateEquationProblem;
    optimizer, max_param_count, eq_complexity_filter::Int = typemax(Int),
    memo::Dict{UInt64,NamedTuple}=Dict{UInt64,NamedTuple}(),
    fitted::Set{UInt64}=Set{UInt64}(),
    parent_of::AbstractDict = Dict(), kwargs...
)
    # Skip structures already produced in an earlier batch — expand each once.
    # Added to `fitted` on first sight regardless of outcome (fit / cap / error).
    fresh = empty(mechs)
    n_fitted_skip = 0
    for m in mechs
        h = hash(m)
        if h in fitted
            n_fitted_skip += 1
        else
            push!(fitted, h)
            push!(fresh, m)
        end
    end

    # PASS 1 (workers): complexity-cap + compile + param-cap + render.
    # `:complexity_skip` = over eq_complexity_filter (checked first, before any
    # derivation); `nothing` = over max_param_count; `FitFailure` = threw; else a
    # record with everything the row needs + `eq_hash`.
    compiled = pmap(fresh) do m
        try
            _eq_complexity(m) > eq_complexity_filter && return :complexity_skip
            em = compile_mechanism(m)
            fkeys = fitted_params(em)
            length(fkeys) > max_param_count && return nothing
            eq_text = rate_equation_string(em)
            (mech = m, orig = m, n_params = length(fkeys),
             mechanism_type = string(typeof(em)),
             eq_text = eq_text, eq_hash = _rate_eq_dedup_key(eq_text),
             fitted_param_names = fkeys)
        catch e
            FitFailure(m, _exc_string(e))
        end
    end

    # Pick one representative per `eq_hash` not already fit, then fit them in
    # PARALLEL (`pmap`) — fitting dominates cost, so it must run across all workers.
    rep_idx = Dict{UInt64,Int}()
    for (i, c) in enumerate(compiled)
        c isa NamedTuple || continue
        (haskey(memo, c.eq_hash) || haskey(rep_idx, c.eq_hash)) && continue
        rep_idx[c.eq_hash] = i
    end
    reps = [(mech = compiled[i].mech, eq_hash = compiled[i].eq_hash)
            for i in values(rep_idx)]
    rep_fits = pmap(reps) do r
        try
            fp = FittingProblem(compile_mechanism(r.mech), prob.data;
                Keq=prob.Keq, scale_k_to_kcat=prob.scale_k_to_kcat)
            (eq_hash = r.eq_hash, fit = fit_rate_equation(fp, optimizer; kwargs...),
             error = nothing)
        catch e
            (eq_hash = r.eq_hash, fit = nothing, error = _exc_string(e))
        end
    end
    fit_error = Dict{UInt64,String}()   # eq_hash → representative fit threw
    for r in rep_fits
        r.error === nothing ? (memo[r.eq_hash] = r.fit) : (fit_error[r.eq_hash] = r.error)
    end

    # Build a row per compiled mechanism by copying its equation's fit.
    entries  = BatchEntry[]
    failures = FitFailure[]
    emitted_eq_hashes = Set{UInt64}()
    for c in compiled
        (c === nothing || c === :complexity_skip) && continue    # cap skip
        c isa FitFailure && (push!(failures, c); continue)
        if haskey(fit_error, c.eq_hash)              # representative fit threw
            push!(failures, FitFailure(c.orig, fit_error[c.eq_hash]))
            continue
        end
        fit = memo[c.eq_hash]
        inherited = !haskey(rep_idx, c.eq_hash) || (c.eq_hash in emitted_eq_hashes)
        push!(emitted_eq_hashes, c.eq_hash)
        parent = get(parent_of, c.orig, nothing)
        row = (
            n_params = c.n_params,
            parent_n_params = parent === nothing ? missing : parent.n_params,
            loss = fit.loss,
            mechanism_type = c.mechanism_type,
            parent_mechanism_type =
                parent === nothing ? missing : parent.mechanism_type,
            rate_equation = c.eq_text,
            retcode = string(fit.retcode),
            error = missing,
            fitted_param_names = c.fitted_param_names,
            fitted_param_values = Tuple(fit.params[k] for k in c.fitted_param_names),
            eq_hash = string(c.eq_hash, base=16, pad=16),
            fit_inherited = inherited,
        )
        push!(entries, BatchEntry(c.mech, c.n_params, fit.loss, fit.retcode,
                                  c.eq_hash, row))
    end
    return entries, failures,
           count(x -> x === nothing, compiled),          # param-count skips
           count(x -> x === :complexity_skip, compiled), # complexity skips
           n_fitted_skip
end

"""
Fold a batch of `BatchEntry`s into the search state: every entry joins the
`frontier` (the unexpanded work queue — ALL structurally-distinct
mechanisms, no eq-dedup); `best_loss_by_count` tracks the per-count running
min (the beam-cutoff reference); `cv_pool` keeps the top `n_cv_candidates`
DISTINCT equations (by `eq_hash`, lowest loss each) per param count.
"""
function _ingest!(frontier, cv_pool, best_loss_by_count, entries;
                  n_cv_candidates)
    for e in entries
        push!(get!(frontier, e.n_params, BatchEntry[]), e)
        if !haskey(best_loss_by_count, e.n_params) ||
                e.loss < best_loss_by_count[e.n_params]
            best_loss_by_count[e.n_params] = e.loss
        end
        _offer_cv!(get!(cv_pool, e.n_params, BatchEntry[]),
                   e, n_cv_candidates)
    end
    nothing
end

"""
Keep `pool` at the top `n` distinct-`eq_hash` entries by loss. A repeat
`eq_hash` only ever updates its own slot (to the lower loss); it never
consumes a second slot.
"""
function _offer_cv!(pool::Vector{BatchEntry}, e::BatchEntry, n::Int)
    n == 0 && return pool
    idx = findfirst(p -> p.eq_hash == e.eq_hash, pool)
    if idx !== nothing
        e.loss < pool[idx].loss && (pool[idx] = e)
        return pool
    end
    if length(pool) < n
        push!(pool, e)
    else
        worst = argmax([p.loss for p in pool])
        e.loss < pool[worst].loss && (pool[worst] = e)
    end
    pool
end

"""
Expand one parent into its children, catching a per-mechanism expansion error
(e.g. a canonicalization that fails to reach a fixed point) so it is recorded as
a failure rather than aborting the whole search. Returns `(children, failure)`
with `failure === nothing` on success, else a `FitFailure` carrying the parent.
"""
function _expand_parent(m::Union{Mechanism, AllostericMechanism},
                        rxn::EnzymeReaction)
    try
        (expand_mechanisms(Union{Mechanism, AllostericMechanism}[m], rxn), nothing)
    catch e
        (Union{Mechanism, AllostericMechanism}[], FitFailure(m, _exc_string(e)))
    end
end

function _beam_search(
    prob::IdentifyRateEquationProblem;
    min_beam_width, loss_rel_threshold, loss_abs_threshold,
    loss_parsimony_threshold,
    max_param_count, eq_complexity_filter, save_dir, show_progress,
    optimizer, n_cv_candidates, kwargs...
)
    frontier = Dict{Int, Vector{BatchEntry}}()
    cv_pool  = Dict{Int, Vector{BatchEntry}}()
    best_loss_by_count = Dict{Int, Float64}()
    # Mechanisms expanded so far per parameter count — the cumulative floor
    # budget. Spent once over the whole search, never re-granted per sweep.
    expanded_by_count = Dict{Int, Int}()
    # Raw pre-rescale fits keyed by `eq_hash`, shared across iterations so each
    # distinct equation is fit exactly once over the whole search.
    memo = Dict{UInt64,NamedTuple}()
    # Structures already produced — each is expanded at most once (termination).
    # This preserves the selected model: expansion is deterministic, so a
    # re-produced structure yields identical children (no new reachable
    # mechanism), and its FIRST production is its best beam-selection chance —
    # the width-floor budget (`expanded_by_count`) only shrinks and the loss
    # cutoff (`best_loss_by_count` running min) only tightens, so anything the
    # old code selected on a later re-production it would already have selected
    # on the first. If `_select_count!` ever gains a non-monotone budget, this
    # invariant breaks.
    fitted = Set{UInt64}()

    # ── Base tier: fit ALL init mechanisms (no bucketing — siblings) ──
    _progress(save_dir, show_progress, "Enumerating initial mechanisms…")
    base = unique!(collect(init_mechanisms(prob.reaction)))
    _progress(save_dir, show_progress,
        "Fitting $(length(base)) initial mechanisms…")
    base_entries, base_failures, n_base_param_skip, n_base_cx_skip, n_base_fitted_skip =
        _process_batch(base, prob;
            optimizer, max_param_count, eq_complexity_filter, memo, fitted, kwargs...)
    if isempty(base_entries)
        isempty(base_failures) && return (
            Union{Mechanism, AllostericMechanism}[],
            _rows_to_dataframe(NamedTuple[]))
        _save_initial_csv(save_dir, [_failure_row(f) for f in base_failures])
        error("Every base-tier fit failed ($(length(base_failures)) " *
              "mechanisms; failure rows written to " *
              "$(joinpath(save_dir, "initial_mechanisms.csv"))). This usually " *
              "indicates an optimizer/solver configuration problem (e.g. an " *
              "unsupported kwarg). First failure: $(base_failures[1].error)")
    end
    _save_initial_csv(save_dir,
        vcat([e.row for e in base_entries],
             [_failure_row(f) for f in base_failures]))
    pre_best = copy(best_loss_by_count)
    _ingest!(frontier, cv_pool, best_loss_by_count,
             base_entries; n_cv_candidates)
    improved = Set(c for c in keys(best_loss_by_count)
                   if !haskey(pre_best, c) || best_loss_by_count[c] < pre_best[c])
    _progress(save_dir, show_progress, string(
        "Base tier: ",
        _batch_summary(base_entries, base_failures;
                       n_param_skipped=n_base_param_skip,
                       n_complexity_skipped=n_base_cx_skip,
                       n_fitted_skipped=n_base_fitted_skip,
                       max_param_count, eq_complexity_filter),
        "\n  ", _best_loss_line(best_loss_by_count, improved)))

    # ── Advancing-target sweep over actual param counts ──
    iteration = 0
    target = minimum(keys(frontier))
    while !isempty(frontier)
        # Sweep this tier plus any same-or-lower-count stragglers.
        swept = BatchEntry[]
        for c in collect(keys(frontier))
            c <= target && append!(swept, pop!(frontier, c))
        end

        to_expand = BatchEntry[]
        for c in unique(e.n_params for e in swept)
            entries_at_count = [e for e in swept if e.n_params == c]
            sel = _select_count!(expanded_by_count, c,
                [e.loss for e in entries_at_count];
                loss_rel_threshold, loss_abs_threshold,
                min_beam_width, best_override = best_loss_by_count[c],
                parsimony_cutoff = _parsimony_cutoff(
                    best_loss_by_count, c, loss_parsimony_threshold))
            append!(to_expand, entries_at_count[sel])
        end

        if !isempty(to_expand)
            # Expand each parent and record which parent produced each child
            # (first parent wins on structural dedup, matching `unique!`), so the
            # saved CSV can carry the parent's round-trippable mechanism type and
            # parameter count for diagnosing per-move parameter changes. The
            # parent's `mechanism_type` is already on its `BatchEntry.row`, so no
            # recompile. Typed for dispatch: expand_mechanisms needs a concrete
            # Vector{<:Union{Mechanism, AllostericMechanism}} eltype.
            parent_of = Dict{Union{Mechanism, AllostericMechanism},
                             @NamedTuple{mechanism_type::String, n_params::Int}}()
            children = Union{Mechanism, AllostericMechanism}[]
            expand_failures = FitFailure[]
            for pe in to_expand
                kids, failure = _expand_parent(pe.mech, prob.reaction)
                failure === nothing || push!(expand_failures, failure)
                for child in kids
                    haskey(parent_of, child) && continue
                    parent_of[child] = (mechanism_type = pe.row.mechanism_type,
                                        n_params = pe.n_params)
                    push!(children, child)
                end
            end
            child_entries, child_failures, n_child_param_skip, n_child_cx_skip,
                n_child_fitted_skip =
                _process_batch(children, prob;
                    optimizer, max_param_count, eq_complexity_filter, memo,
                    fitted, parent_of, kwargs...)
            # A per-parent expansion error is recorded like a fit failure (CSV row
            # + the errored bucket), so a canonicalization bug flags itself in the
            # search output instead of aborting the whole run.
            append!(child_failures, expand_failures)
            if !isempty(child_entries) || !isempty(child_failures)
                # Count only iterations that produced rows, so the
                # equation_search_iteration_N CSVs are gap-free.
                iteration += 1
                _save_iteration_csv(save_dir,
                    vcat([e.row for e in child_entries],
                         [_failure_row(f) for f in child_failures]),
                    iteration)
                pre_best = copy(best_loss_by_count)
                !isempty(child_entries) && _ingest!(
                    frontier, cv_pool, best_loss_by_count,
                    child_entries; n_cv_candidates)
                improved = Set(c for c in keys(best_loss_by_count)
                    if !haskey(pre_best, c) || best_loss_by_count[c] < pre_best[c])
                np_range = isempty(child_entries) ? "n/a" :
                    let lo = minimum(e -> e.n_params, child_entries),
                        hi = maximum(e -> e.n_params, child_entries)
                        lo == hi ? string(lo) : "$lo-$hi"
                    end
                _progress(save_dir, show_progress, string(
                    "Iteration $iteration (child n_params $np_range): ",
                    length(to_expand), " parents → ", length(children),
                    " children\n  ",
                    _batch_summary(child_entries, child_failures;
                                   n_param_skipped=n_child_param_skip,
                                   n_complexity_skipped=n_child_cx_skip,
                                   n_fitted_skipped=n_child_fitted_skip,
                                   max_param_count, eq_complexity_filter),
                    "\n  ", _best_loss_line(best_loss_by_count, improved)))
            elseif !isempty(children)
                # Whole batch cap-skipped (param count and/or complexity) or
                # already fit, so it produced no rows and no CSV. Report it
                # anyway; don't bump the iteration counter, since there are no
                # rows to save.
                _progress(save_dir, show_progress, string(
                    "Expanded ", length(to_expand), " parents → ",
                    length(children), " children | all skipped (",
                    n_child_fitted_skip, " already fit, ",
                    n_child_param_skip, " >", max_param_count, " params, ",
                    n_child_cx_skip, " >", eq_complexity_filter, " complexity)"))
            end
        end

        isempty(frontier) && break
        target = max(target + 1, minimum(keys(frontier)))
    end

    pool_entries = BatchEntry[e for v in values(cv_pool) for e in v]
    mechs = Union{Mechanism, AllostericMechanism}[
        e.mech for e in pool_entries]
    df = _rows_to_dataframe([e.row for e in pool_entries])
    return mechs, df
end

# Parsimony reference = threshold × best loss over ALL counts strictly below c
# (not just c-1): an added parameter must beat the best simpler model of any size.
# Returns nothing when no simpler tier has been fit yet.
function _parsimony_cutoff(best_loss_by_count::Dict{Int,Float64}, c::Int,
                           loss_parsimony_threshold::Float64)
    prev = [best_loss_by_count[k] for k in keys(best_loss_by_count) if k < c]
    isempty(prev) && return nothing
    loss_parsimony_threshold * minimum(prev)
end

"""
Subset a columnar NamedTuple by a boolean mask, returning views
to avoid per-fold copying. Called G times per LOOCV; copying
every column ×2 (train/test) ×G grew O(G·N·ncols).
"""
function _subset_data(data::NamedTuple, mask)
    idx = findall(mask)
    return map(col -> view(col, idx), data)
end

"""
Evaluate loss of a mechanism on data with given
params and scaling mode.
"""
function _evaluate_loss(
    mechanism, data, params, Keq, scale_k_to_kcat
)
    pnames = fitted_params(mechanism)
    x = [log(params[p]) for p in pnames]
    fp = FittingProblem(mechanism, data; Keq=Keq, scale_k_to_kcat=scale_k_to_kcat)
    return loss!(x, fp)
end

"""
    _onesided_permutation_p(diffs; exact_threshold=20,
                             mc_samples=10^6,
                             rng=Random.default_rng()) → Float64

One-sided p-value `Pr(perm_mean ≥ observed)` for paired-difference vector
`diffs` under the sign-flip null. Exact enumeration of `2^n` sign patterns
when `length(diffs) ≤ exact_threshold` (default 20 → up to ~10^6 perms);
Monte Carlo with `mc_samples` random sign-flips otherwise.

Both default branches do ~10^6 inner iterations. The `exact_threshold` and
`mc_samples` kwargs are exposed primarily for tests (forcing the MC branch
on small fixtures and using seeded RNGs).

Errors on empty `diffs` (caller's invariant: a bucket comparison always has
at least one fold).
"""
function _onesided_permutation_p(
    diffs::Vector{Float64};
    exact_threshold::Int = 20,
    mc_samples::Int = 10^6,
    rng = Random.default_rng(),
)
    n = length(diffs)
    n == 0 && error("_onesided_permutation_p: empty diffs vector")
    # Compute `observed` with the same sequential reduction as the inner
    # permutation loop, so the identity permutation reproduces it
    # bit-identically and `s/n >= observed` always counts it. Using
    # `mean(diffs)` here would invoke pairwise summation for n ≥ 16,
    # which differs from the loop at 1 ULP and can drop the identity
    # from the count.
    observed_sum = 0.0
    for i in 1:n
        observed_sum += diffs[i]
    end
    observed = observed_sum / n

    if n <= exact_threshold
        total = 1 << n   # 2^n
        count_ge = 0
        for mask in 0:(total - 1)
            s = 0.0
            @inbounds for i in 1:n
                bit = (mask >> (i - 1)) & 1
                s += bit == 1 ? -diffs[i] : diffs[i]
            end
            count_ge += (s / n >= observed)
        end
        return count_ge / total
    else
        count_ge = 0
        @inbounds for _ in 1:mc_samples
            s = 0.0
            for i in 1:n
                s += rand(rng, Bool) ? -diffs[i] : diffs[i]
            end
            count_ge += (s / n >= observed)
        end
        return count_ge / mc_samples
    end
end

"""
    _select_best_n_params(cv_df; se_threshold=1.0,
                          perm_p_threshold=0.16) → NamedTuple

Paired 1-SE rule AND-combined with a one-sided sign-flip permutation test
on log-transformed per-fold LOOCV scores. Returns:

  best_n::Int                — selected `n_params`
  n_min::Int                 — bucket with lowest mean log-fold-loss
  diagnostics::Dict{Int, NamedTuple{(:mean_log_loss_diff, :se_paired,
                                     :permutation_p)}}
                              — `n_min` bucket has all three = 0.0

Per-bucket representative = the row with the lowest `cv_score` in that
`n_params` bucket. For each non-`n_min` bucket, computes paired diffs vs
the rep of `n_min`'s log-folds. The simpler bucket is accepted iff BOTH:

  mean(diffs) ≤ se_threshold * std(diffs)/sqrt(n_folds)   (paired 1-SE)
  permutation_p > perm_p_threshold                         (perm test)

Iterates smaller buckets in ascending `n_params`; returns the first that
passes. Falls through to `n_min` if none pass. Larger-than-`n_min` buckets
have diagnostics computed but are never selected.

Tiebreak: when two buckets tie on `mean(log_scores)`, `n_min` resolves to
the smallest `n_params` (parsimony).

Errors if any bucket's fold-score length differs from `n_min`'s. (Every
`cv_fold_scores` row is non-empty — the LOOCV grid scores every
`(candidate, fold)` pair or raises — so there is no empty-input case to handle.)

When `n_folds_min == 1` the SE is undefined; the selection loop is
skipped and `n_min` is returned. Diagnostics are still populated with
`se_paired = 0.0`. When the input has only one `n_params` value, returns
it as both `n_min` and `best_n` with empty smaller/larger comparisons.
"""
function _select_best_n_params(
    cv_df::DataFrame;
    se_threshold::Float64 = 1.0,
    perm_p_threshold::Float64 = 0.16,
)
    sorted = sort(cv_df, [:n_params, :cv_score])
    reps = combine(groupby(sorted, :n_params), first)

    log_scores = Dict(
        row.n_params => log.(row.cv_fold_scores)
        for row in eachrow(reps))
    log_means = Dict(n => mean(ls) for (n, ls) in log_scores)
    # Tie-break on log-mean by smallest n_params (parsimony). Without
    # this, argmin over Dict keys is iteration-order-dependent.
    n_min = argmin(n -> (log_means[n], n), keys(log_means))
    n_folds_min = length(log_scores[n_min])

    diagnostics = Dict{Int, @NamedTuple{
        mean_log_loss_diff::Float64,
        se_paired::Float64,
        permutation_p::Float64}}()
    diagnostics[n_min] = (
        mean_log_loss_diff = 0.0,
        se_paired = 0.0,
        permutation_p = 0.0,
    )

    smaller_ns = sort([n for n in keys(log_means) if n < n_min])

    for n in keys(log_means)
        n == n_min && continue
        ls = log_scores[n]
        length(ls) == n_folds_min || error(
            "fold-count mismatch for n_params=$n: " *
            "got $(length(ls)), expected $n_folds_min " *
            "(n_min=$n_min)")
        diffs = ls .- log_scores[n_min]
        md  = mean(diffs)
        sep = n_folds_min == 1 ? 0.0 :
              std(diffs) / sqrt(n_folds_min)
        p   = _onesided_permutation_p(diffs)
        diagnostics[n] = (
            mean_log_loss_diff = md,
            se_paired = sep,
            permutation_p = p,
        )
    end

    best_n = n_min
    if n_folds_min > 1
        for n in smaller_ns
            d = diagnostics[n]
            if d.mean_log_loss_diff <= se_threshold * d.se_paired &&
               d.permutation_p > perm_p_threshold
                best_n = n
                break
            end
        end
    end

    return (
        best_n = best_n,
        n_min = n_min,
        diagnostics = diagnostics,
    )
end

function _cv_model_selection(
    mechs::Vector, df::DataFrame,
    prob::IdentifyRateEquationProblem;
    n_cv_candidates, optimizer,
    se_threshold::Float64,
    perm_p_threshold::Float64,
    save_dir, show_progress,
    kwargs...
)
    isempty(mechs) && error(
        "No mechanisms were successfully " *
        "fitted during beam search")

    # LOOCV candidates: the top `n_cv_candidates` DISTINCT equations (by
    # `eq_hash`, lowest `loss` each) per `n_params` bucket. Deduping by `eq_hash`
    # keeps textually identical equations out of LOOCV — at most one row per
    # `eq_hash` per param count — so folds are never wasted on, or biased by, them.
    candidate_indices = Int[]
    df_idx = DataFrame(
        row_idx = 1:nrow(df), n_params = df.n_params,
        loss = df.loss, eq_hash = df.eq_hash,
    )
    for gdf in groupby(df_idx, :n_params)
        seen_hashes = Set{String}()
        for row in eachrow(sort(gdf, :loss))
            row.eq_hash in seen_hashes && continue
            push!(seen_hashes, row.eq_hash)
            push!(candidate_indices, row.row_idx)
            length(seen_hashes) >= n_cv_candidates && break
        end
    end
    candidate_mechs = mechs[candidate_indices]
    candidate_rows = df[candidate_indices, :]

    _progress(save_dir, show_progress,
        "Cross-validating $(length(candidate_mechs)) candidate equations (LOOCV)…")
    # Flatten LOOCV to a (candidate, fold) grid so all folds of all candidates
    # run across every worker, not one candidate per worker with serial folds.
    groups = unique(prob.data.group)
    tasks = [(ci, g) for ci in eachindex(candidate_mechs) for g in groups]
    flat = pmap(tasks) do task
        ci, g = task
        m = compile_mechanism(candidate_mechs[ci])
        (ci, g, _cv_fold_loss(m, prob, g; optimizer, kwargs...))
    end
    fold_scores_per_candidate = _scatter_fold_scores(
        flat, length(candidate_mechs), groups)

    cv_df = copy(candidate_rows)
    cv_df.cv_fold_scores = collect(fold_scores_per_candidate)
    cv_df.cv_score = [mean(log.(v)) for v in cv_df.cv_fold_scores]

    sel = _select_best_n_params(
        cv_df;
        se_threshold = se_threshold,
        perm_p_threshold = perm_p_threshold,
    )

    # Best mechanism = lowest training `loss` within sel.best_n.
    at_best_pc_idx = findall(==(sel.best_n), cv_df.n_params)
    isempty(at_best_pc_idx) && error(
        "internal: best_n=$(sel.best_n) has no rows in cv_df")
    sort_perm = sortperm(cv_df.loss[at_best_pc_idx])
    best_row_idx = at_best_pc_idx[sort_perm[1]]
    best_mech = candidate_mechs[best_row_idx]
    best_mechanism = compile_mechanism(best_mech)

    # Populate diagnostic columns. A bucket may be absent from
    # diagnostics only if every row in it had empty fold scores.
    for fld in (:mean_log_loss_diff, :se_paired, :permutation_p)
        cv_df[!, fld] = [
            haskey(sel.diagnostics, n) ?
                sel.diagnostics[n][fld] : missing
            for n in cv_df.n_params
        ]
    end

    # Flatten per-fold scores into one column per held-out group.
    # Group order matches the `groups = unique(prob.data.group)` iteration above.
    for (i, g) in enumerate(groups)
        col = Symbol("cv_fold_$g")
        cv_df[!, col] = [v[i] for v in cv_df.cv_fold_scores]
    end

    _progress(save_dir, show_progress,
        "Selected: $(nameof(typeof(best_mechanism))) " *
        "(eq_hash=$(cv_df.eq_hash[best_row_idx])), n_params=$(sel.best_n)")
    select!(cv_df, Not(:cv_fold_scores))

    # Save the LOOCV table and the selected best equation alongside the
    # per-iteration fit CSVs, so cluster runs persist the model-selection
    # outcome (recoverable from the iteration CSVs, but wasteful to omit).
    isdir(save_dir) || mkpath(save_dir)
    CSV.write(joinpath(save_dir, "loocv_results.csv"), cv_df)
    CSV.write(joinpath(save_dir, "best_equation.csv"), cv_df[[best_row_idx], :])

    return IdentifyRateEquationResults(best_mechanism, cv_df)
end

"""
One LOOCV fold: fit `mechanism` on every group except `held_out`, score it on
`held_out`, and return the test loss floored at `eps(Float64)`. A non-finite
test loss raises (naming the held-out group) — a corrupted fold must abort model
selection rather than propagate a bad score.
"""
function _cv_fold_loss(
    mechanism::AbstractEnzymeMechanism,
    prob::IdentifyRateEquationProblem, held_out;
    optimizer, kwargs...)
    train_mask = prob.data.group .!= held_out
    test_mask  = prob.data.group .== held_out
    train_data = _subset_data(prob.data, train_mask)
    test_data  = _subset_data(prob.data, test_mask)

    fp_train = FittingProblem(mechanism, train_data;
        Keq=prob.Keq, scale_k_to_kcat=prob.scale_k_to_kcat)
    fit = fit_rate_equation(fp_train, optimizer; kwargs...)

    test_loss = _evaluate_loss(mechanism, test_data,
        fit.params, prob.Keq, prob.scale_k_to_kcat)
    # A non-finite fold loss means the fit is unusable; aborting model
    # selection is correct (re-run CV from the saved CSVs after fixing
    # the fit). max(NaN, eps) === NaN, so the floor below would not catch it.
    isfinite(test_loss) || error(
        "LOOCV produced a non-finite test loss for held-out group " *
        "$held_out — the fit is unusable; aborting model selection.")
    # Floor at eps so log(score) is finite. The centered-residuals loss
    # can be exactly 0 (e.g. a single-row held-out group).
    max(test_loss, eps(Float64))
end

"""
Scatter flat `(candidate_index, group, score)` triples into one fold-score
vector per candidate, each ordered by `groups`. Every `(ci, g)` in the grid
appears exactly once, so every slot is written.
"""
function _scatter_fold_scores(flat, n_candidates::Int, groups)
    gi = Dict(g => i for (i, g) in enumerate(groups))
    out = [Vector{Float64}(undef, length(groups)) for _ in 1:n_candidates]
    for (ci, g, s) in flat
        out[ci][gi[g]] = s
    end
    out
end

"""
Default results directory: the first non-existent `<date>_results[_N]`
directory in the cwd (e.g. `YYYY_MM_DD_results`, then `…_results_2`, `_3`).
"""
function _default_save_dir()
    base = string(Dates.format(Dates.today(), "yyyy_mm_dd"), "_results")
    isdir(base) || return base
    n = 2
    while isdir(string(base, "_", n)); n += 1; end
    string(base, "_", n)
end
