# ABOUTME: Beam-search pipeline to identify the best rate equation.
# ABOUTME: Enumerates mechanisms, fits each, selects via CV; a comment-
# ABOUTME: stripped rate-equation string key dedups equivalent equations.

using DataFrames
using CSV
using Statistics

"""
    IdentifyRateEquationProblem{D}

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
struct IdentifyRateEquationProblem{D<:NamedTuple}
    reaction::EnzymeReaction
    data::D
    Keq::Float64
    scale_k_to_kcat::Union{Float64,Nothing}
end

function IdentifyRateEquationProblem(
    reaction::EnzymeReaction, table; Keq::Real,
    scale_k_to_kcat::Union{Real,Nothing}=1.0
)
    # Every metabolite the reaction declares needs a concentration column.
    data = _rate_table(table, _metabolite_names(reaction), scale_k_to_kcat, Keq)

    # Validate at least 2 groups for CV
    n_groups = length(unique(data.group))
    n_groups >= 2 || error(
        "Need at least 2 unique groups for " *
        "cross-validation, got $n_groups")

    IdentifyRateEquationProblem{typeof(data)}(reaction, data, Keq, scale_k_to_kcat)
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
        min_beam_width=50, loss_rel_threshold=1.3, loss_abs_threshold=0.001,
        loss_parsimony_threshold=0.99,
        max_param_count=20, n_restarts=20, maxtime=600.0, maxiters=10_000_000,
        abstol=nothing, reltol=nothing, callback=nothing, solver_kwargs=(;),
        n_cv_candidates=5, se_threshold=1.0,
        save_dir=_default_save_dir(), show_progress=true)

Find the best rate equation for the given reaction
and data using beam search.

# Keyword Arguments
- `min_beam_width::Int = 50`: minimum mechanisms
  to keep per param-count tier
- `loss_rel_threshold::Float64 = 1.3`: relative tolerance
  for beam selection (see "Beam selection" below)
- `loss_abs_threshold::Float64 = 0.001`: absolute tolerance
  for beam selection
- `loss_parsimony_threshold::Float64 = 0.99`: a mechanism
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
- `optional_allosteric_regulators::Vector{Symbol} = Symbol[]`: declared
  allosteric regulators the seed may omit. A reaction that declares regulators
  seeds the beam fully regulated by default; naming a regulator here lets the
  beam also start from seeds that do not bind it (added back later by
  expansion). Listing every declared allosteric regulator restores the
  unregulated `init_mechanisms` seed.
- `optional_competitive_inhibitors::Vector{Symbol} = Symbol[]`: declared
  competitive inhibitors the seed may omit, as above.
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
- `se_threshold::Float64 = 1.0`: 1-SE multiplier for model selection.
  A simpler equation is accepted iff its `cv_score` is `≤` the best
  equation's `cv_score + se_threshold * cv_score_se`. Default 1.0 is
  the textbook "1-SE rule".
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

Every LOOCV candidate gets a `cv_score` (mean of its per-fold losses)
and a `cv_score_se` (`std(fold losses)/sqrt(n_folds)`). The best
equation is the one with the lowest `cv_score`; its cutoff is
`cv_score + se_threshold * cv_score_se`. The selected equation is the
lowest-`cv_score` candidate at the smallest `n_params` that has a
candidate at or below the cutoff, so a simpler equation wins whenever
it predicts held-out groups as well as the best one to within the
best one's fold-to-fold standard error.
"""
function identify_rate_equation(
    prob::IdentifyRateEquationProblem;
    # Beam search
    min_beam_width::Int = 50,
    loss_rel_threshold::Float64 = 1.3,
    loss_abs_threshold::Float64 = 0.001,
    loss_parsimony_threshold::Float64 = 0.99,
    max_param_count::Int = 20,
    eq_complexity_filter::Int = 337,
    optional_allosteric_regulators::Vector{Symbol} = Symbol[],
    optional_competitive_inhibitors::Vector{Symbol} = Symbol[],
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
    # Output & parallelism
    save_dir::String = _default_save_dir(),
    show_progress::Bool = true,
)
    fitting_kwargs = (;
        n_restarts, maxtime, maxiters,
        abstol, reltol, callback, solver_kwargs)

    isdir(save_dir) &&
        any(f -> endswith(f, ".csv") || f == "progress.log", readdir(save_dir)) &&
        error("save_dir already contains results (CSV files or progress.log). " *
              "Use an empty directory to avoid mixing results.")

    mechanisms, df = _beam_search(prob;
        min_beam_width, loss_rel_threshold,
        loss_abs_threshold, loss_parsimony_threshold,
        max_param_count, eq_complexity_filter, save_dir, show_progress,
        optimizer, n_cv_candidates,
        optional_allosteric_regulators, optional_competitive_inhibitors,
        fitting_kwargs...)

    result = _cv_model_selection(
        mechanisms, df, prob;
        se_threshold, optimizer, save_dir, show_progress,
        fitting_kwargs...)
    _progress(save_dir, show_progress, "Done. Results saved to $save_dir")
    return result
end

"""Write result rows to `<save_dir>/<filename>`, creating `save_dir` if absent."""
function _write_rows_csv(save_dir::String, filename::String, rows)
    mkpath(save_dir)
    CSV.write(joinpath(save_dir, filename), _rows_to_dataframe(rows))
end

"""
Convert result row NamedTuples to a DataFrame, one column per fitted parameter
name across `rows` (`missing` where a row lacks it).
Row order is preserved (no sorting) to maintain
alignment with the mechanism vector.
"""
function _rows_to_dataframe(rows)
    isempty(rows) && return DataFrame()
    sorted_pnames = sort!(unique(k for r in rows for k in keys(r.params)))

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
        df[!, pn] = [get(r.params, pn, missing) for r in rows]
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
Select this sweep's parents among the count-`c` mechanisms with `losses`, and advance
the cumulative floor budget. Returns the selected indices in input order. A mechanism
qualifies if either:
  • its loss ≤ cutoff, where
    cutoff = min(loss_rel_threshold * best_loss_by_count[c] + loss_abs_threshold,
                 loss_parsimony_threshold * best(<c))
    and best(<c) is the best loss over ALL counts strictly below c (not just c-1): an
    added parameter must beat the best simpler model of any size. The parsimony term
    only tightens the cutoff, and is dropped while no smaller count has been fit,
  • OR its rank (1-indexed by ascending loss) is within the *cumulative* floor budget.
    `expanded[c]` tracks how many count-`c` mechanisms the whole search has expanded
    so far; the width floor may add at most `min_beam_width - expanded[c]` more. Once
    the budget is spent, only the loss cutoff admits at that count.

Mechanisms with non-finite losses (`Inf`, `NaN`) are excluded
unconditionally — they represent failed or non-converging fits
that should not propagate to the next level.
"""
function _select_count!(
    expanded::Dict{Int,Int}, best_loss_by_count::Dict{Int,Float64}, c::Int,
    losses::AbstractVector{<:Real};
    loss_rel_threshold::Float64, loss_abs_threshold::Float64,
    loss_parsimony_threshold::Float64, min_beam_width::Int,
)
    cutoff = loss_rel_threshold * best_loss_by_count[c] + loss_abs_threshold
    smaller = [l for (k, l) in best_loss_by_count if k < c]
    isempty(smaller) ||
        (cutoff = min(cutoff, loss_parsimony_threshold * minimum(smaller)))
    budget = max(0, min_beam_width - get(expanded, c, 0))
    perm = sort([i for i in eachindex(losses) if isfinite(losses[i])]; by=i -> losses[i])
    # Return indices in original (input) order so callers don't
    # rely on the by-loss sort order, which is a side-effect of
    # the rank computation rather than part of the contract.
    sel = sort!([i for (rank, i) in enumerate(perm)
                 if losses[i] <= cutoff || rank <= budget])
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

"""Compact, CSV-safe rendering of a thrown exception: type + truncated message."""
_exc_string(e) = first(sprint(showerror, e), 200)

"""
CSV row for a mechanism that threw. Same NamedTuple schema as a fitted row,
with `missing` wherever the value is unavailable (compile/fit never produced it).
`mechanism_type` is the round-trippable parametric `EnzymeMechanism{Sig}` string when
the mechanism compiles; falls back to the bare concrete type name
(`"EnzymeRates.Mechanism"` / `"EnzymeRates.AllostericMechanism"`) when compilation
itself fails, so the row still identifies the mechanism family.
"""
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
     params = (;),
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
    mkpath(save_dir)
    open(joinpath(save_dir, "progress.log"), "a") do io
        println(io, msg)
    end
    nothing
end

"""
Pre-fit half of a batch summary: the five buckets known before fitting —
`new fits` (distinct new equations being fit this batch), `inherited` (rows
reusing a memoized or in-batch fit), and the three `skipped` buckets
(already-fit, over `max_param_count`, over `eq_complexity_filter`).
"""
function _prefit_summary(n_new::Int, n_inherited::Int,
        n_param_skip::Int, n_cx_skip::Int, n_fitted_skip::Int;
        max_param_count::Int, eq_complexity_filter::Int)
    string(n_new, " new fits + ", n_inherited, " inherited + ",
           n_fitted_skip, " skipped (already fit) + ",
           n_param_skip, " skipped (>", max_param_count, " params) + ",
           n_cx_skip, " skipped (>", eq_complexity_filter, " complexity)")
end

"""
Post-fit half of a batch summary: `errored` (`length(failures)` — compile, fit,
and any expansion failures) and the `Success` / `non-Success retcode`
percentages over the fitted set (`entries`).
"""
function _postfit_summary(entries::Vector{BatchEntry}, failures::Vector{FitFailure})
    n_fit  = length(entries)
    n_succ = count(e -> e.retcode === :Success, entries)
    pct(x) = n_fit == 0 ? 0.0 : round(100 * x / n_fit; digits=1)
    string(length(failures), " errored | Success ", pct(n_succ),
           "% | non-Success retcode ", pct(n_fit - n_succ), "%")
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
Fit one batch of mechanisms. Returns `(entries, failures, n_param_skip, n_cx_skip,
n_seen_skip)`.

PASS-1: skip mechanisms already produced, in an earlier batch or earlier in this one (a
mechanism whose structural `hash` is already in `seen` is dropped before compiling;
`seen` gets every structure on first sight regardless of outcome). Compile, cap-check
(`max_param_count`), and render every fresh mechanism in parallel (`pmap`) to get its
`eq_hash`; a mechanism over `eq_complexity_filter` or `max_param_count` is recorded as a
skip, not fit. Then pick one fit representative per `eq_hash` not already in `memo`, and
pass the pre-fit summary (`_prefit_summary`) to `log`. `mechs` is already structurally
deduped by the caller (`unique!`).

PASS-2: fit ONE representative per `eq_hash` not already in `memo`, in parallel across
all workers (`pmap`); the fit is stored in `memo` and copied to every mechanism sharing
that `eq_hash` — same equation ⟹ same fitted params and same rescaling, so no refit is
needed. `memo` persists across iterations (threaded from `_beam_search`), so an equation
fit in an earlier iteration is never refit. A representative whose fit throws fails ALL
of its duplicates (all-or-nothing per equation). `fit_inherited` is `false` for the
representative actually fit this batch, `true` for every reused row. `parent_of` maps a
child to its parent's `BatchEntry`, which fills the row's parent columns.
`entries::Vector{BatchEntry}` are the mechanisms with a usable fit (each keeping its own
row, `retcode`, and `eq_hash`); `failures::Vector{FitFailure}` are mechanisms that threw
at compile/render/fit — captured WITH the exception text, never silently swallowed.
"""
function _process_batch(
    mechs, prob::IdentifyRateEquationProblem;
    optimizer, max_param_count, eq_complexity_filter::Int = typemax(Int),
    memo::Dict{UInt64,NamedTuple}=Dict{UInt64,NamedTuple}(),
    seen::Set{UInt64}=Set{UInt64}(),
    parent_of::AbstractDict = Dict(), log = msg -> nothing, kwargs...
)
    # Skip structures already produced — expand each once. Added to `seen` on first
    # sight regardless of outcome (fit / cap / error).
    fresh = empty(mechs)
    for m in mechs
        h = hash(m)
        h in seen || (push!(seen, h); push!(fresh, m))
    end

    # PASS 1 (workers): complexity-cap + compile + param-cap + render.
    # `:complexity_skip` = over eq_complexity_filter (checked first, before any
    # derivation); `:param_skip` = over max_param_count; `FitFailure` = threw; else a
    # record with everything the row needs + `eq_hash`.
    compiled = pmap(fresh) do m
        try
            _eq_complexity(m) > eq_complexity_filter && return :complexity_skip
            em = compile_mechanism(m)
            fkeys = fitted_params(em)
            length(fkeys) > max_param_count && return :param_skip
            eq_text = rate_equation_string(em)
            (mech = m, n_params = length(fkeys),
             mechanism_type = string(typeof(em)),
             eq_text = eq_text, eq_hash = _rate_eq_dedup_key(eq_text))
        catch e
            FitFailure(m, _exc_string(e))
        end
    end

    # Pick one representative per `eq_hash` not already fit.
    rep_idx = Dict{UInt64,Int}()
    for (i, c) in enumerate(compiled)
        c isa NamedTuple && !haskey(memo, c.eq_hash) && get!(rep_idx, c.eq_hash, i)
    end
    reps = [(mech = compiled[i].mech, eq_hash = compiled[i].eq_hash)
            for i in values(rep_idx)]
    n_param_skip = count(==(:param_skip), compiled)
    n_cx_skip = count(==(:complexity_skip), compiled)
    n_seen_skip = length(mechs) - length(fresh)
    log(_prefit_summary(length(reps),
        count(c -> c isa NamedTuple, compiled) - length(reps),
        n_param_skip, n_cx_skip, n_seen_skip; max_param_count, eq_complexity_filter))

    # PASS 2: fit them in PARALLEL (`pmap`) — fitting dominates cost, so it must run
    # across all workers.
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
        c isa Symbol && continue                     # cap skip
        c isa FitFailure && (push!(failures, c); continue)
        if haskey(fit_error, c.eq_hash)              # representative fit threw
            push!(failures, FitFailure(c.mech, fit_error[c.eq_hash]))
            continue
        end
        fit = memo[c.eq_hash]
        inherited = !haskey(rep_idx, c.eq_hash) || (c.eq_hash in emitted_eq_hashes)
        push!(emitted_eq_hashes, c.eq_hash)
        parent = get(parent_of, c.mech, nothing)
        row = (
            n_params = c.n_params,
            parent_n_params = parent === nothing ? missing : parent.n_params,
            loss = fit.loss,
            mechanism_type = c.mechanism_type,
            parent_mechanism_type =
                parent === nothing ? missing : parent.row.mechanism_type,
            rate_equation = c.eq_text,
            retcode = string(fit.retcode),
            error = missing,
            params = fit.params,
            eq_hash = string(c.eq_hash, base=16, pad=16),
            fit_inherited = inherited,
        )
        push!(entries, BatchEntry(c.mech, c.n_params, fit.loss, fit.retcode,
                                  c.eq_hash, row))
    end
    (entries, failures, n_param_skip, n_cx_skip, n_seen_skip)
end

"""
Fold a batch of `BatchEntry`s into the search state: every entry joins the
`frontier` (the unexpanded work queue — ALL structurally-distinct
mechanisms, no eq-dedup); `best_loss_by_count` tracks the per-count running
min (the beam-cutoff reference); `cv_pool` keeps the top `n_cv_candidates`
DISTINCT equations (by `eq_hash`, lowest loss each) per param count. Returns the
set of counts whose best loss strictly dropped (or first appeared) in this batch.
"""
function _ingest!(frontier, cv_pool, best_loss_by_count, entries;
                  n_cv_candidates)
    improved = Set{Int}()
    for e in entries
        push!(get!(frontier, e.n_params, BatchEntry[]), e)
        if !haskey(best_loss_by_count, e.n_params) ||
                e.loss < best_loss_by_count[e.n_params]
            best_loss_by_count[e.n_params] = e.loss
            push!(improved, e.n_params)
        end
        _offer_cv!(get!(cv_pool, e.n_params, BatchEntry[]),
                   e, n_cv_candidates)
    end
    improved
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
Expand every selected parent into its children across the workers, then merge
serially. `pmap` preserves input order, so iterating `zip(to_expand, results)`
reproduces the serial loop's first-parent-wins dedup, child order, and failure
order exactly. Returns `(children, parent_of, expand_failures)`: `parent_of` maps
each child to the `BatchEntry` of the first parent that produced it, and
`expand_failures` holds a `FitFailure` carrying the parent for each parent whose
expansion threw.
"""
function _expand_parents(to_expand::Vector{BatchEntry},
                         reaction::EnzymeReaction)
    results = pmap([pe.mech for pe in to_expand]) do m
        # Record a per-parent expansion error (e.g. `expand_mechanisms`' atom-conservation
        # assertion, `_assert_atom_conserving`) as a failure; never abort the search.
        try
            (expand_mechanisms(Union{Mechanism, AllostericMechanism}[m], reaction),
             nothing)
        catch e
            (Union{Mechanism, AllostericMechanism}[], FitFailure(m, _exc_string(e)))
        end
    end
    parent_of = Dict{Union{Mechanism, AllostericMechanism}, BatchEntry}()
    children = Union{Mechanism, AllostericMechanism}[]
    expand_failures = FitFailure[]
    for (pe, (kids, failure)) in zip(to_expand, results)
        failure === nothing || push!(expand_failures, failure)
        for child in kids
            haskey(parent_of, child) && continue
            parent_of[child] = pe
            push!(children, child)
        end
    end
    (children, parent_of, expand_failures)
end

"""
The mechanisms the beam fits first: every seed of `mechs` that is not degenerate
(`_degenerate`), and in place of each one that is, its flip children (`_expand_re_to_ss`)
that are not. A degenerate seed's law ignores a substrate or never saturates, so it is not
worth a fit, but it is the only parent of mechanisms that are not degenerate (such as the
ping-pong seeds' children). Only a flip can cure it: the other moves keep every catalytic
step's rapid-equilibrium flag, and the steps they add bind regulators, which neither
degeneracy condition reads. An expansion error is returned as a failure, one per seed,
carrying the seed.
"""
function _base_tier(mechs::Vector, rxn::EnzymeReaction)
    base = Union{Mechanism, AllostericMechanism}[m for m in mechs if !_degenerate(m)]
    failures = FitFailure[]
    for m in mechs
        _degenerate(m) || continue
        try
            for child in _expand_re_to_ss(m)
                _assert_atom_conserving(child)
                _degenerate(child) || push!(base, child)
            end
        catch e
            push!(failures, FitFailure(m, _exc_string(e)))
        end
    end
    unique!(base), failures
end

"""
    _required_regulators(rxn, optional_allosteric_regulators,
                         optional_competitive_inhibitors)
        -> (required_allo::Set{Symbol}, required_comp::Set{Symbol})

The regulators the beam seed must bind: every `AllostericRegulator` and every
`CompetitiveInhibitor` declared in `rxn`, minus the names the caller marked
optional. Both sets empty means the beam takes its seeds from `init_mechanisms`; a
non-empty set means it takes them from `seed_mechanisms`. Either way the base tier
fits the seeds that are not degenerate and, in place of those that are, their flip
children that are not (`_base_tier`).
"""
function _required_regulators(rxn::EnzymeReaction,
                              optional_allosteric_regulators::Vector{Symbol},
                              optional_competitive_inhibitors::Vector{Symbol})
    required(T, optional) = setdiff(
        Set{Symbol}(name(regulator(rm)) for rm in regulators(rxn)
                    if regulator(rm) isa T),
        optional)
    (required(AllostericRegulator, optional_allosteric_regulators),
     required(CompetitiveInhibitor, optional_competitive_inhibitors))
end

function _beam_search(
    prob::IdentifyRateEquationProblem;
    min_beam_width, loss_rel_threshold, loss_abs_threshold,
    loss_parsimony_threshold,
    max_param_count, eq_complexity_filter, save_dir, show_progress,
    optimizer, n_cv_candidates,
    optional_allosteric_regulators, optional_competitive_inhibitors,
    kwargs...
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
    seen = Set{UInt64}()

    progress(msg) = _progress(save_dir, show_progress, msg)
    # Fit one batch and record it: log `header` with the pre-fit summary, write the
    # rows to `csv`, fold the entries into the search state, and log the post-fit
    # summary. Returns `(entries, failures)`, or `nothing`, writing no CSV, when the
    # batch has no rows.
    function run_batch!(mechs, extra_failures, header, csv; parent_of = Dict())
        entries, failures = _process_batch(mechs, prob; optimizer, max_param_count,
            eq_complexity_filter, memo, seen, parent_of,
            log = msg -> progress(header * "\n  " * msg), kwargs...)
        append!(failures, extra_failures)
        isempty(entries) && isempty(failures) && return nothing
        _write_rows_csv(save_dir, csv,
            vcat([e.row for e in entries], [_failure_row(f) for f in failures]))
        improved = _ingest!(frontier, cv_pool, best_loss_by_count, entries;
                            n_cv_candidates)
        progress(string("  ", _postfit_summary(entries, failures),
                        "\n  ", _best_loss_line(best_loss_by_count, improved)))
        (entries, failures)
    end

    # ── Base tier: fit every seed, with each degenerate seed replaced by its
    # non-degenerate flip children (`_base_tier`; no bucketing — siblings) ──
    progress("Enumerating initial mechanisms…")
    required_allo, required_comp = _required_regulators(
        prob.reaction, optional_allosteric_regulators,
        optional_competitive_inhibitors)
    seeds = (isempty(required_allo) && isempty(required_comp)) ?
        unique!(collect(init_mechanisms(prob.reaction))) :
        unique!(collect(seed_mechanisms(
            prob.reaction, required_allo, required_comp)))
    base, base_expand_failures = _base_tier(seeds, prob.reaction)
    # A degenerate seed's expansion error is recorded like a fit failure (a row of
    # initial_mechanisms.csv and the errored bucket).
    base_batch = run_batch!(base, base_expand_failures,
        "Fitting $(length(base)) initial mechanisms…", "initial_mechanisms.csv")
    base_batch === nothing && return (
        Union{Mechanism, AllostericMechanism}[], _rows_to_dataframe(NamedTuple[]))
    base_entries, base_failures = base_batch
    isempty(base_entries) && error(
        "Every base-tier fit failed ($(length(base_failures)) " *
        "mechanisms; failure rows written to " *
        "$(joinpath(save_dir, "initial_mechanisms.csv"))). This usually " *
        "indicates an optimizer/solver configuration problem (e.g. an " *
        "unsupported kwarg). First failure: $(base_failures[1].error)")

    # ── Advancing-target sweep over actual param counts ──
    iteration = 0
    target = 0
    while !isempty(frontier)
        # Sweep this tier plus any same-or-lower-count stragglers, bucket by bucket
        # in the frontier's key order (each bucket holds one count).
        target = max(target + 1, minimum(keys(frontier)))
        to_expand = BatchEntry[]
        for c in [k for k in keys(frontier) if k <= target]
            entries_at_count = pop!(frontier, c)
            sel = _select_count!(expanded_by_count, best_loss_by_count, c,
                [e.loss for e in entries_at_count];
                loss_rel_threshold, loss_abs_threshold, loss_parsimony_threshold,
                min_beam_width)
            append!(to_expand, entries_at_count[sel])
        end
        isempty(to_expand) && continue

        # Expand each parent and record which parent produced each child
        # (first parent wins on structural dedup, matching `unique!`), so the
        # saved CSV can carry the parent's round-trippable mechanism type and
        # parameter count for diagnosing per-move parameter changes. The
        # parent's `mechanism_type` is already on its `BatchEntry.row`, so no
        # recompile. Typed for dispatch: expand_mechanisms needs a concrete
        # Vector{<:Union{Mechanism, AllostericMechanism}} eltype.
        children, parent_of, expand_failures =
            _expand_parents(to_expand, prob.reaction)
        # A per-parent expansion error is recorded like a fit failure (CSV row
        # + the errored bucket), so a bug in an expansion move flags itself in
        # the search output instead of aborting the whole run. Count only
        # iterations that produced rows, so the equation_search_iteration_N CSVs
        # are gap-free: a batch with no rows (every child skipped) still logs its
        # header and pre-fit summary, but writes no CSV and keeps the counter.
        # `iteration` is a 1-based sequential counter, NOT a parameter count —
        # the real fitted count is the `n_params` column of each row.
        run_batch!(children, expand_failures,
            "Iteration $(iteration + 1): $(length(to_expand)) parents → " *
            "$(length(children)) children",
            "equation_search_iteration_$(iteration + 1).csv"; parent_of) === nothing ||
            (iteration += 1)
    end

    pool_entries = BatchEntry[e for v in values(cv_pool) for e in v]
    mechs = Union{Mechanism, AllostericMechanism}[
        e.mech for e in pool_entries]
    df = _rows_to_dataframe([e.row for e in pool_entries])
    return mechs, df
end

"""
    _select_best_row(cv_df; se_threshold=1.0) → Int

1-SE rule on the best equation's fold scores. The best row is the one with
the lowest `cv_score` (ties resolve to the smallest `n_params`); its cutoff
is `cv_score + se_threshold * cv_score_se`. Returns the index of the
lowest-`cv_score` row at the smallest `n_params` that has a row at or below
the cutoff. The best row always passes its own cutoff, so a row with more
parameters than the best row is never returned.
"""
function _select_best_row(cv_df::DataFrame; se_threshold::Float64 = 1.0)
    by_score = sortperm(collect(zip(cv_df.cv_score, cv_df.n_params)))
    i_min = by_score[1]
    cutoff = cv_df.cv_score[i_min] + se_threshold * cv_df.cv_score_se[i_min]
    passing = [i for i in by_score if cv_df.cv_score[i] <= cutoff]
    best_n = minimum(cv_df.n_params[passing])
    return first(i for i in passing if cv_df.n_params[i] == best_n)
end

function _cv_model_selection(
    mechs::Vector, df::DataFrame,
    prob::IdentifyRateEquationProblem;
    optimizer, se_threshold::Float64,
    save_dir, show_progress,
    kwargs...
)
    isempty(mechs) && error(
        "No mechanisms were successfully " *
        "fitted during beam search")

    # LOOCV candidates: every row of `df`, ordered by `n_params`, then `loss`. The
    # rows are the cv pool, where `_offer_cv!` keeps the top `n_cv_candidates`
    # DISTINCT equations (by `eq_hash`, lowest `loss` each) per `n_params` bucket.
    # That keeps textually identical equations out of LOOCV — at most one row per
    # `eq_hash` per param count — so folds are never wasted on, or biased by, them.
    candidate_indices = sortperm(collect(zip(df.n_params, df.loss)))
    candidate_mechs = mechs[candidate_indices]
    candidate_rows = df[candidate_indices, :]

    _progress(save_dir, show_progress,
        "Cross-validating $(length(candidate_mechs)) candidate equations (LOOCV)…")
    # Flatten LOOCV to a (candidate, fold) grid so all folds of all candidates
    # run across every worker, not one candidate per worker with serial folds.
    # `pmap` keeps the grid order (the fold varies fastest), so column ci of
    # `scores` holds candidate ci's fold losses in `groups` order.
    groups = unique(prob.data.group)
    fold_losses = pmap([(m, g) for m in candidate_mechs for g in groups]) do (m, g)
        _cv_fold_loss(compile_mechanism(m), prob, g; optimizer, kwargs...)
    end
    scores = reshape(Float64.(fold_losses), length(groups), :)

    cv_df = copy(candidate_rows)
    cv_df.cv_score = [mean(c) for c in eachcol(scores)]
    cv_df.cv_score_se = [std(c) / sqrt(length(groups)) for c in eachcol(scores)]

    best_row_idx = _select_best_row(cv_df; se_threshold)
    best_mechanism = compile_mechanism(candidate_mechs[best_row_idx])

    # Flatten per-fold scores into one column per held-out group.
    # Group order matches the `groups = unique(prob.data.group)` iteration above.
    for (i, g) in enumerate(groups)
        cv_df[!, Symbol("cv_fold_$g")] = scores[i, :]
    end

    _progress(save_dir, show_progress,
        "Selected: $(nameof(typeof(best_mechanism))) " *
        "(eq_hash=$(cv_df.eq_hash[best_row_idx])), " *
        "n_params=$(cv_df.n_params[best_row_idx])")

    # Save the LOOCV table and the selected best equation alongside the
    # per-iteration fit CSVs, so cluster runs persist the model-selection
    # outcome (recoverable from the iteration CSVs, but wasteful to omit).
    mkpath(save_dir)
    CSV.write(joinpath(save_dir, "loocv_results.csv"), cv_df)
    CSV.write(joinpath(save_dir, "best_equation.csv"), cv_df[[best_row_idx], :])

    return IdentifyRateEquationResults(best_mechanism, cv_df)
end

"""
One LOOCV fold: fit `mechanism` on every group except `held_out`, score it on
`held_out`, and return the finite test loss. A non-finite
test loss raises (naming the held-out group) — a corrupted fold must abort model
selection rather than propagate a bad score.
"""
function _cv_fold_loss(
    mechanism::AbstractEnzymeMechanism,
    prob::IdentifyRateEquationProblem, held_out;
    optimizer, kwargs...)
    held = prob.data.group .== held_out
    # Each fold's columns are views, so no fold copies the data: copying every
    # column for train and test in each of the G folds would cost O(G·N·ncols).
    fold(mask) = (idx = findall(mask);
        FittingProblem(mechanism, map(col -> view(col, idx), prob.data);
            Keq=prob.Keq, scale_k_to_kcat=prob.scale_k_to_kcat))
    fit = fit_rate_equation(fold(.!held), optimizer; kwargs...)
    # `fit.params` is keyed by `fitted_params(mechanism)`, the order `loss!` reads.
    test_loss = loss!([log(v) for v in fit.params], fold(held))
    # A non-finite fold loss means the fit is unusable; aborting model
    # selection is correct (re-run CV from the saved CSVs after fixing
    # the fit).
    isfinite(test_loss) || error(
        "LOOCV produced a non-finite test loss for held-out group " *
        "$held_out — the fit is unusable; aborting model selection.")
    test_loss
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
