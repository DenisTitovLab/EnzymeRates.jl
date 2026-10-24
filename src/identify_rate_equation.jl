# ABOUTME: Beam-search pipeline to identify the best rate equation.
# ABOUTME: Enumerates mechanisms, fits one representative per rate-law key
# ABOUTME: (`_canonical`), and selects among the fits by cross-validation.

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
- `absolute`: `false` (the default) treats the data as relative and scales every fit to
  peak turnover 1; `true` treats it as absolute per-enzyme turnover (see `FittingProblem`)
- `peak_range`: the concentration range `(lo, hi)` of the peak-turnover scale (see
  `FittingProblem`), in the data's concentration units; the default covers physiological
  concentrations in M or mM, and data in µM need a higher upper bound

Construct it with `IdentifyRateEquationProblem(reaction, table; Keq, absolute = false,
peak_range = (1e-12, 1e3))`.
"""
struct IdentifyRateEquationProblem{D<:NamedTuple}
    reaction::EnzymeReaction
    data::D
    Keq::Float64
    absolute::Bool
    peak_range::Tuple{Float64,Float64}
end

function IdentifyRateEquationProblem(
    reaction::EnzymeReaction, table; Keq::Real,
    absolute::Bool=false, peak_range=DEFAULT_PEAK_RANGE
)
    mnames = _metabolite_names(reaction)
    # The search fits allosteric mechanisms, and their rate equations take the
    # conformational constant L, so a metabolite named L would fail every one of them.
    # The reaction's constructor rejects the name L (`_RESERVED_NAMES`), so this check
    # cannot fire for a reaction the constructor builds.
    :L in mnames && error(
        "Metabolite L has the name of the conformational constant L of the allosteric " *
        "mechanisms the search fits; rename the metabolite")
    # Every metabolite the reaction declares needs a concentration column.
    data = _rate_table(table, mnames, Keq)

    # Validate at least 2 groups for CV
    n_groups = length(unique(data.group))
    n_groups >= 2 || error(
        "Need at least 2 unique groups for " *
        "cross-validation, got $n_groups")

    IdentifyRateEquationProblem{typeof(data)}(reaction, data, Keq, absolute,
                                              _peak_range(peak_range))
end

"""
    IdentifyRateEquationResults

Results from `identify_rate_equation`.

# Fields
- `best`: the selected `AbstractEnzymeMechanism` (1-SE rule on the LOOCV
  scores; see Model selection in `identify_rate_equation`)
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
        max_param_count=20, eq_complexity_filter=337,
        optional_allosteric_regulators=Symbol[],
        optional_competitive_inhibitors=Symbol[],
        n_restarts=20, maxtime=600.0, maxiters=10_000_000,
        abstol=nothing, reltol=nothing, callback=nothing, solver_kwargs=(;),
        n_cv_candidates=5, se_threshold=1.0,
        save_dir=_default_save_dir(), show_progress=true)

Find the best rate equation for the given reaction
and data using beam search.

# Keyword Arguments
- `min_beam_width::Int = 50`: minimum rate-law keys to expand per
  param-count tier; a sweep expands every mechanism it meets whose key is
  selected
- `loss_rel_threshold::Float64 = 1.3`: relative tolerance
  for beam selection (see "Beam selection" below)
- `loss_abs_threshold::Float64 = 0.001`: absolute tolerance
  for beam selection
- `loss_parsimony_threshold::Float64 = 0.99`: the loss cutoff
  admits a key only if its loss is within this factor of
  the best model of any smaller parameter count — an added
  parameter must earn its keep. Combined with the other loss
  thresholds via `min`; `min_beam_width` is a cumulative
  per-count budget (see "Beam selection" below). `Inf` disables it.
- `max_param_count::Int = 20`: skip, before fitting, any mechanism with
  more fitted constants than this (`fitted_params`), so the search stops
  at this parameter count
- `eq_complexity_filter::Int = 337`: skip (before fitting, before any
  derivation) any mechanism whose rate equation is more complex than this —
  roughly the number of terms in the rate equation's denominator: V×τ, the
  number of RE-segments times the spanning-tree count of the segment graph,
  i.e. the products the compiled equation evaluates on every call. An
  allosteric mechanism counts both conformations, V×τ of the active one plus
  V×τ of the inactive one. The default 337 is one above a fully steady-state
  random-order bi-bi (V×τ = 336), so such a mechanism passes and anything more
  complex is skipped — equations denser than that are impractical to fit and
  can blow up derivation/codegen. A one-subunit variant of that bi-bi whose
  inactive conformation keeps its substrate bindings counts 336 + 16 = 352 and
  is skipped. Computed from the mechanism graph alone.
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
- `maxtime::Real = 600.0`: max time per optimizer run in seconds (common solver
  option, forwarded to `Optimization.solve`); each of a fit's `n_restarts` runs
  gets it, so a fit can take up to `n_restarts × maxtime`
- `maxiters::Integer = 10_000_000`: max iterations per optimizer run
  (common solver option, forwarded to `Optimization.solve`)
- `abstol`/`reltol`/`callback = nothing`: Optimization.jl common solver
  options, forwarded to `Optimization.solve` only when set
- `solver_kwargs::NamedTuple = (;)`: solver-specific options forwarded
  verbatim to `Optimization.solve` (e.g. `(; popsize=200)` for a CMA-ES
  solver that supports it); the caller matches its contents to `optimizer`
- `n_cv_candidates::Int = 5`: LOOCV top N candidates with **distinct rate-law
  keys** per param count (at least 1)
- `se_threshold::Float64 = 1.0`: 1-SE multiplier for model selection.
  A simpler equation is accepted iff its `cv_score` is `≤` the best
  equation's `cv_score + se_threshold * cv_score_se`. Default 1.0 is
  the textbook "1-SE rule".
- `save_dir::String = _default_save_dir()`: output directory for the
  search CSVs (`initial_mechanisms.csv` + `equation_search_iteration_N.csv`),
  plus `loocv_results.csv` (the LOOCV table for every cross-validated
  candidate) and `best_equation.csv` (the selected equation and its
  fitted parameters). A `save_dir` that already holds a `.csv` file or a
  `progress.log` is refused, so the results of two runs never mix.

# Beam selection

The beam selects rate-law keys: the mechanisms of one key share one fit
(an L phantom's key aside, see below), so a key counts once in the ranks
and the budget, and a sweep expands every mechanism it meets of a selected
key. A key at parameter count `n` qualifies for the next-level beam if
either:
- its loss ≤ `min(loss_rel_threshold * best(n) + loss_abs_threshold,
  loss_parsimony_threshold * best(<n))`, where `best(n)` is the
  lowest loss at parameter count `n` and `best(<n)` is the lowest
  loss over all smaller counts; the second term is dropped at the
  base count (no smaller level exists yet),
- OR it falls under the `min_beam_width` budget: a cumulative
  per-count allowance that expands at least `min_beam_width`
  keys at each count over the whole search, then stops —
  spent once, not re-granted each time the count is revisited.

A selected key keeps its verdict: a mechanism found later with a
key already selected expands. One found later with a key already
rejected is judged again and, at the same loss, rejected again,
since the cutoff only tightens and a rejection by rank spends the
budget. Only an L phantom's key can be refused and later selected: a
member met before the fit of its law at `L = 0` is fitted itself and
one met after reads that fit, so the key can carry two losses, and a
member refused at one has left the frontier when the other selects
the key.

The additive term protects against `best_loss` approaching zero
(simulated / very-low-loss data) where a purely multiplicative
threshold would collapse the beam to the single best key.

# Model selection (LOOCV)

Every LOOCV candidate gets a `cv_score` (mean of its per-fold losses)
and a `cv_score_se` (`std(fold losses)/sqrt(n_folds)`). The best
equation is the one with the lowest `cv_score`; its cutoff is
`cv_score + se_threshold * cv_score_se`. The selected equation is the
lowest-`cv_score` candidate at the smallest `n_params` that has a
candidate at or below the cutoff, so a simpler equation wins whenever
it predicts held-out groups as well as the best one to within the
best one's fold-to-fold standard error.

A candidate's fold losses come from fitting, on each fold, the
mechanism whose fit it reads. An L phantom (`L_phantom` in the CSVs:
a one-subunit variant whose `L` the data cannot determine) takes its
source's fold losses, so beside its source it loses the tie by its
extra parameter, and selected without it, it reports `L = 0`.
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
    n_cv_candidates >= 1 || error("n_cv_candidates must be ≥ 1; got $n_cv_candidates")
    fitting_kwargs = (;
        n_restarts, maxtime, maxiters,
        abstol, reltol, callback, solver_kwargs)

    isdir(save_dir) &&
        any(f -> endswith(f, ".csv") || f == "progress.log", readdir(save_dir)) &&
        error("save_dir already contains results (CSV files or progress.log). " *
              "Use an empty directory to avoid mixing results.")

    pool, memo = _beam_search(prob;
        min_beam_width, loss_rel_threshold,
        loss_abs_threshold, loss_parsimony_threshold,
        max_param_count, eq_complexity_filter, save_dir, show_progress,
        optimizer, n_cv_candidates,
        optional_allosteric_regulators, optional_competitive_inhibitors,
        fitting_kwargs...)

    result = _cv_model_selection(
        pool, memo, prob;
        se_threshold, optimizer, save_dir, show_progress,
        fitting_kwargs...)
    _progress(save_dir, show_progress, "Done. Results saved to $save_dir")
    return result
end

"""
Write result rows to `<save_dir>/<filename>`, creating `save_dir` if absent. `texts` maps
each `mechanism_id` the run has written to its mechanism's text, and gains the rows' ids.
A row whose id names another text is an error, raised before the file is written: the
run's CSVs could not tell the two mechanisms apart.
"""
function _write_rows_csv!(save_dir::String, filename::String, rows,
                          texts::Dict{String,String})
    for r in rows
        other = get!(texts, r.mechanism_id, r.mechanism)
        other == r.mechanism ||
            error("two mechanisms share the mechanism_id $(r.mechanism_id), the first " *
                  "16 hex digits of the SHA-256 of their texts, so the run's CSVs cannot " *
                  "tell them apart:\n$other\n$(r.mechanism)")
    end
    mkpath(save_dir)
    CSV.write(joinpath(save_dir, filename), _rows_to_dataframe(rows))
end

"""
Convert result row NamedTuples to a DataFrame, one column per fitted parameter
name across `rows` (`missing` where a row lacks it).
Row order is preserved (no sorting), so row i stays aligned with the i-th entry the rows
came from (the LOOCV candidates in `_cv_model_selection`).
"""
function _rows_to_dataframe(rows)
    isempty(rows) && return DataFrame()
    sorted_pnames = sort!(unique(k for r in rows for k in keys(r.params)))

    df = DataFrame(
        n_params = [r.n_params for r in rows],
        parent_n_params =
            Union{Missing,Int}[get(r, :parent_n_params, missing) for r in rows],
        loss = [r.loss for r in rows],
        mechanism = [r.mechanism for r in rows],
        mechanism_id = [r.mechanism_id for r in rows],
        parent_id = Union{Missing,String}[get(r, :parent_id, missing) for r in rows],
        rate_equation = [r.rate_equation for r in rows],
        retcode = Union{Missing,String}[r.retcode for r in rows],
        error = Union{Missing,String}[r.error for r in rows],
        key = Union{Missing,String}[r.key for r in rows],
        fit_inherited = Union{Missing,Bool}[r.fit_inherited for r in rows],
        fit_source_id = Union{Missing,String}[r.fit_source_id for r in rows],
        L_phantom = Union{Missing,Bool}[r.L_phantom for r in rows],
        peak_at_range_edge =
            Union{Missing,Bool}[get(r, :peak_at_range_edge, missing) for r in rows],
    )
    for pn in sorted_pnames
        df[!, pn] = [get(r.params, pn, missing) for r in rows]
    end
    df
end

"""
A mechanism's id in identify's CSVs: the first 16 hex digits of the SHA-256 of `text`,
the mechanism's DSL text (`show`), so every process and Julia version gives one mechanism
one id.
"""
_mechanism_id(text::AbstractString) = bytes2hex(sha256(text))[1:16]

"""
Select this sweep's parents among the count-`c` mechanisms with `losses` and rate-law
keys `law_keys`, and advance the cumulative floor budget. The beam ranks keys, not
mechanisms: the mechanisms of one key share one fit (`_process_batch`; an L phantom's key
aside, whose members met before and after its source's fit read two), so a key enters at its
lowest loss, takes one rank and spends one unit of budget, and duplicates cannot crowd out a
distinct key. A key not yet in `selected` is selected if either:
  • its loss ≤ cutoff, where
    cutoff = min(loss_rel_threshold * best_loss_by_count[c] + loss_abs_threshold,
                 loss_parsimony_threshold * best(<c))
    and best(<c) is the best loss over ALL counts strictly below c (not just c-1): an
    added parameter must beat the best simpler model of any size. The parsimony term
    only tightens the cutoff, and is dropped while no smaller count has been fit,
  • OR its rank (1-indexed by ascending loss among this sweep's keys not yet selected) is
    within the *cumulative* floor budget. `expanded[c]` counts the keys at count `c` the
    whole search has selected so far; the width floor may add at most
    `min_beam_width - expanded[c]` more. Once the budget is spent, only the loss cutoff
    admits at that count.
A selected key joins `selected` and keeps that verdict for the rest of the search: a
later member of it is returned without spending budget. A rejected key is judged again by
the same rule when a later sweep meets it, and at the same loss it is rejected again,
since the cutoff only tightens and a rejection by rank spends the budget. Returns the
indices, in input order, of every mechanism whose key is selected.

Every loss is finite: a failed or non-converging fit (`Inf`, `NaN`) should not
propagate to the next level, so `_ingest!` keeps it out of the frontier.
"""
function _select_count!(
    expanded::Dict{Int,Int}, best_loss_by_count::Dict{Int,Float64}, c::Int,
    losses::AbstractVector{<:Real}, law_keys::AbstractVector{String};
    loss_rel_threshold::Float64, loss_abs_threshold::Float64,
    loss_parsimony_threshold::Float64, min_beam_width::Int,
    selected::Set{String},
)
    cutoff = loss_rel_threshold * best_loss_by_count[c] + loss_abs_threshold
    smaller = [l for (k, l) in best_loss_by_count if k < c]
    isempty(smaller) ||
        (cutoff = min(cutoff, loss_parsimony_threshold * minimum(smaller)))
    budget = max(0, min_beam_width - get(expanded, c, 0))
    # The keys not yet selected, in input order, each at its lowest loss.
    unselected, lowest = String[], Dict{String,Float64}()
    for (l, k) in zip(losses, law_keys)
        k in selected && continue
        haskey(lowest, k) || push!(unselected, k)
        lowest[k] = min(get(lowest, k, Inf), l)
    end
    for (rank, j) in enumerate(sortperm([lowest[k] for k in unselected]))
        k = unselected[j]
        (lowest[k] <= cutoff || rank <= budget) && push!(selected, k)
    end
    expanded[c] = get(expanded, c, 0) + count(in(selected), unselected)
    # Indices in original (input) order, so callers do not rely on the by-loss sort
    # order, a side-effect of the rank computation rather than part of the contract.
    findall(in(selected), law_keys)
end

"""
A fit as identify stores it, once per rate-law key (`_canonical`): `values` holds the
fitted constants at their canonical positions (`_Canonical.order`), then `L` for an MWC
law, so every mechanism with that key reads them through its own canonical order;
`loss`, `retcode` and `peak_at_range_edge` are `fit_rate_equation`'s; `source` is the
mechanism that was fitted and `source_id` its `mechanism_id`.
"""
struct _Fit
    values::Vector{Float64}
    loss::Float64
    retcode::Symbol
    peak_at_range_edge::Bool
    source::Union{Mechanism, AllostericMechanism}
    source_id::String
end

"""
A mechanism PASS1 of `_process_batch` derived: its fitted constants in `fitted_params`
order, its text and the text's id (`_mechanism_id`), its printed law, its law's canonical
form (`canon`), whether it is an L phantom that reads another law's fit (`_is_phantom`),
and `canon0`, the canonical form of the law whose fit it reads: for a phantom its active
law, its law at `L = 0`; `canon` otherwise.
"""
struct _Candidate
    mech::Union{Mechanism, AllostericMechanism}
    fitted::Vector{Symbol}
    mechanism::String
    mechanism_id::String
    eq_text::String
    canon::_Canonical
    phantom::Bool
    canon0::_Canonical
end

"""The order in which candidate `c` reads its fit's values (`_Fit`): the canonical order of
the law whose fit it reads (`canon0`), then `L` for an MWC law."""
_fit_order(c::_Candidate) =
    c.mech isa AllostericMechanism ? [c.canon0.order; :L] : c.canon0.order

"""
The parameters of candidate `c` as its row holds them, in its `fitted_params` order: a
fit's `values`, stored at canonical positions, read through `order`, the candidate's own
order for those positions.
"""
function _row_params(c::_Candidate, order::Vector{Symbol}, values::Vector{Float64})
    at = Dict(zip(order, values))
    NamedTuple{Tuple(c.fitted)}(Tuple(at[x] for x in c.fitted))
end

"""One fitted mechanism: its parameter count, its fit's loss and retcode, its rate-law key
(`_canonical`), the key of the fit it reads (`fit_key`: its own, or for an L phantom the
key of its law at `L = 0`), and its CSV row."""
struct BatchEntry
    mech::Union{Mechanism, AllostericMechanism}
    n_params::Int
    loss::Float64
    retcode::Symbol
    key::String
    fit_key::String
    row::NamedTuple
end

"""A mechanism that failed at derivation, compilation, fit or expansion, with the
exception message captured as a string. Kept for the CSV + summary."""
struct FitFailure
    mech::Union{Mechanism, AllostericMechanism}
    error::String
end

"""Compact, CSV-safe rendering of a thrown exception: type + truncated message."""
_exc_string(e) = first(sprint(showerror, e), 200)

"""
CSV row for a mechanism that threw. Same NamedTuple schema as a fitted row; it records
only the mechanism's text and id and the error (`missing` elsewhere, empty `params`).
`mechanism` is the mechanism's text, rendered from the concrete mechanism, so a row
names its mechanism even when it does not compile; `parent_id` is missing.
"""
function _failure_row(f::FitFailure)
    text = sprint(show, f.mech)
    (n_params = missing,
     parent_n_params = missing,
     loss = missing,
     mechanism = text,
     mechanism_id = _mechanism_id(text),
     parent_id = missing,
     rate_equation = missing,
     retcode = missing,
     error = f.error,
     params = (;),
     key = missing,
     fit_inherited = missing,
     fit_source_id = missing,
     L_phantom = missing,
     peak_at_range_edge = missing)
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
`new fits` (distinct new laws, one per key, being fit this batch), `inherited` (rows
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
Post-fit half of a batch summary: `errored` (`length(failures)` — derivation,
compile, fit, and any expansion failures) and the `Success` / `non-Success retcode`
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
mechanism whose structural `hash` is already in `seen` is dropped before it is derived;
`seen` gets every structure on first sight regardless of outcome). For each fresh
mechanism, in parallel (`pmap`) and from the concrete mechanism alone, without compiling
it: skip it over `eq_complexity_filter`; check it with the guard (`_check_well_posed`);
skip it over `max_param_count` by its solve's fitted constants (`_solve`); else derive its
law (`_law`), check the law on the guard's faces (`_assert_support`), print it
(`_law_text`), take its text and the text's id (`_mechanism_id`) and label the law
canonically (`_canonical`). A law whose labelling reaches the leaf cap gets a key of its
own, and the pre-fit summary names its mechanism. Then pick one fit representative per key
not already in `memo`, and pass the pre-fit summary (`_prefit_summary`) to `log`. `mechs`
is already structurally deduped by the caller (`unique!`).

PASS-2: fit ONE representative per key not already in `memo`, in parallel across all
workers (`pmap`), and store the fit in `memo` at canonical positions (`_Fit`). Every
mechanism with that key reads it through its own canonical order, so its row carries its
own parameter names: same law ⟹ same fitted values and same rescaling, so no refit is
needed. `memo` persists across iterations (threaded from `_beam_search`), so a law fit in
an earlier iteration is never refit. A representative whose fit throws fails ALL the
mechanisms with its key (all-or-nothing per law). An L phantom (`_is_phantom`: a rescaling
certificate, which needs one catalytic subunit and an inactive conformation that turns
nothing over) is never a representative. It reads the fit of its law at `L = 0`, its active
law's key (`fit_key`), from `memo` or from this batch's representatives, through its own
order and with `L = 0`, and its row sets `L_phantom`. A phantom whose law at `L = 0` has no
fit there, or lacks one of its other fitted constants, is fitted as an ordinary candidate;
one whose source's fit throws fails with the source's error. A row's
`fit_inherited` is `false` when its mechanism is its fit's source and `true` otherwise, and
`fit_source_id` names the source. `parent_of` maps a child to its parent's `BatchEntry`,
which fills the row's parent columns. `entries::Vector{BatchEntry}` are the mechanisms with
a usable fit (each keeping its own row, `retcode` and key); `failures::Vector{FitFailure}`
are mechanisms that threw at derivation or fit — captured WITH the exception text,
never silently swallowed.
"""
function _process_batch(
    mechs, prob::IdentifyRateEquationProblem;
    optimizer, max_param_count, eq_complexity_filter::Int = typemax(Int),
    memo::Dict{String,_Fit}=Dict{String,_Fit}(),
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

    # PASS 1 (workers), from the concrete mechanism and compiling nothing: the complexity
    # filter (checked first, before any derivation), the guard, the solve and the
    # parameter cap, the law, its check on the guard's faces, the printed law and the
    # canonical form. `:complexity_skip` and `:param_skip` are skips; `FitFailure` = threw;
    # else a `_Candidate`.
    records = pmap(fresh) do m
        try
            _eq_complexity(m) > eq_complexity_filter && return :complexity_skip
            _check_well_posed(m)
            s = _solve(m)
            length(s.fitted) > max_param_count && return :param_skip
            _with_law(m, s) do m, s, L
                eq_text = _law_text(m, s, L; coefficients = true)
                text = sprint(show, m)
                id = _mechanism_id(text)
                canon = _canonical(L, id)
                # An L phantom reads the fit of its law at L = 0, its active law, when that
                # law holds every fitted constant but L. A regulatory site whose ligands are
                # all :EqualAI has one factor in both conformations, so the certificate
                # admits it, yet the active law lacks its constants: the check fits such a
                # phantom itself, so a row never reads names its fit lacks.
                phantom = m isa AllostericMechanism && _is_phantom(L)
                canon0 = phantom ? _canonical(L.active, id) : canon
                if phantom && Set(canon0.order) != Set(s.fitted[1:(end - 1)])
                    phantom, canon0 = false, canon
                end
                _Candidate(m, s.fitted, text, id, eq_text, canon, phantom, canon0)
            end
        catch e
            FitFailure(m, _exc_string(e))
        end
    end

    # Pick one representative per key not already fit; an L phantom is none.
    rep_idx = Dict{String,Int}()
    for (i, c) in enumerate(records)
        c isa _Candidate && !c.phantom && !haskey(memo, c.canon.key) &&
            get!(rep_idx, c.canon.key, i)
    end
    # A phantom inherits the fit of whoever holds its L → 0 key, in `memo` or among this
    # batch's representatives; any other phantom becomes an ordinary candidate.
    for (i, c) in enumerate(records)
        c isa _Candidate && c.phantom || continue
        (haskey(memo, c.canon0.key) || haskey(rep_idx, c.canon0.key)) && continue
        records[i] = _Candidate(c.mech, c.fitted, c.mechanism, c.mechanism_id,
                                c.eq_text, c.canon, false, c.canon)
        haskey(memo, c.canon.key) || get!(rep_idx, c.canon.key, i)
    end
    reps = _Candidate[records[i] for i in values(rep_idx)]
    n_param_skip = count(==(:param_skip), records)
    n_cx_skip = count(==(:complexity_skip), records)
    n_seen_skip = length(mechs) - length(fresh)
    capped = ["leaf cap reached for $(c.mechanism_id); fitted under its own key"
              for c in records if c isa _Candidate && c.canon.capped]
    log(join([_prefit_summary(length(reps),
                  count(c -> c isa _Candidate, records) - length(reps),
                  n_param_skip, n_cx_skip, n_seen_skip;
                  max_param_count, eq_complexity_filter); capped], "\n  "))

    # PASS 2: fit them in PARALLEL (`pmap`) — fitting dominates cost, so it must run
    # across all workers. Each representative gives its fit at canonical positions, or
    # the text of the error its fit threw.
    rep_fits = pmap(reps) do r
        try
            fp = FittingProblem(compile_mechanism(r.mech), prob.data;
                Keq=prob.Keq, absolute=prob.absolute, peak_range=prob.peak_range)
            fit = fit_rate_equation(fp, optimizer; kwargs...)
            _Fit(Float64[fit.params[x] for x in _fit_order(r)], fit.loss, fit.retcode,
                 fit.peak_at_range_edge, r.mech, r.mechanism_id)
        catch e
            _exc_string(e)
        end
    end
    fit_error = Dict{String,String}()   # key → its representative's fit threw
    for (r, f) in zip(reps, rep_fits)
        f isa _Fit ? (memo[r.canon.key] = f) : (fit_error[r.canon.key] = f::String)
    end

    # Build a row per derived mechanism from its key's fit, read through its own order.
    entries  = BatchEntry[]
    failures = FitFailure[]
    for c in records
        c isa FitFailure && (push!(failures, c); continue)
        c isa _Candidate || continue                       # a skip
        fit_key = c.canon0.key                    # its own key, or a phantom's L → 0 key
        if haskey(fit_error, fit_key)             # the fit's source threw
            push!(failures, FitFailure(c.mech, fit_error[fit_key]))
            continue
        end
        fit = memo[fit_key]
        params = _row_params(c, _fit_order(c), c.phantom ? [fit.values; 0.0] : fit.values)
        parent = get(parent_of, c.mech, nothing)
        row = (
            n_params = length(c.fitted),
            parent_n_params = parent === nothing ? missing : parent.n_params,
            loss = fit.loss,
            mechanism = c.mechanism,
            mechanism_id = c.mechanism_id,
            parent_id = parent === nothing ? missing : parent.row.mechanism_id,
            rate_equation = c.eq_text,
            retcode = string(fit.retcode),
            error = missing,
            params,
            key = c.canon.key,
            fit_inherited = fit.source_id != c.mechanism_id,
            fit_source_id = fit.source_id,
            L_phantom = c.phantom,
            peak_at_range_edge = fit.peak_at_range_edge,
        )
        push!(entries, BatchEntry(c.mech, length(c.fitted), fit.loss, fit.retcode,
                                  c.canon.key, fit_key, row))
    end
    (entries, failures, n_param_skip, n_cx_skip, n_seen_skip)
end

"""
Fold a batch of `BatchEntry`s into the search state: every entry with a finite loss
joins the `frontier` (the unexpanded work queue — ALL structurally-distinct
mechanisms, no key dedup); `best_loss_by_count` tracks the per-count running
min (the beam-cutoff reference); `cv_pool` keeps the top `n_cv_candidates`
entries with DISTINCT rate-law keys (lowest loss each) per param count. An entry whose
loss is not finite (`fit_rate_equation` returns `Inf` when no restart is finite) is
skipped: it is never expanded, never sets a best loss and never enters LOOCV; its row
is already in the batch's CSV. Returns the set of counts whose best loss strictly
dropped (or first appeared) in this batch.
"""
function _ingest!(frontier, cv_pool, best_loss_by_count, entries;
                  n_cv_candidates)
    improved = Set{Int}()
    for e in entries
        isfinite(e.loss) || continue
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
Keep `pool` at the top `n` entries with distinct rate-law keys, by loss. A repeat key
only ever updates its own slot (to the lower loss); it never consumes a second slot.
"""
function _offer_cv!(pool::Vector{BatchEntry}, e::BatchEntry, n::Int)
    idx = findfirst(p -> p.key == e.key, pool)
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
One sweep of the beam: pop every `frontier` bucket whose count is at most `target`, in the
order the Dict yields them (each bucket holds one count), and return the entries whose
rate-law key `_select_count!` selects. `expanded`, `best_loss_by_count` and `selected` are
the search's state, which `_select_count!` reads and advances.
"""
function _sweep_parents!(
    frontier::Dict{Int, Vector{BatchEntry}}, target::Int, expanded::Dict{Int, Int},
    best_loss_by_count::Dict{Int, Float64}, selected::Set{String};
    loss_rel_threshold, loss_abs_threshold, loss_parsimony_threshold, min_beam_width,
)
    parents = BatchEntry[]
    for c in filter(<=(target), collect(keys(frontier)))
        entries_at_count = pop!(frontier, c)
        sel = _select_count!(expanded, best_loss_by_count, c,
            [e.loss for e in entries_at_count], [e.key for e in entries_at_count];
            loss_rel_threshold, loss_abs_threshold, loss_parsimony_threshold,
            min_beam_width, selected)
        append!(parents, entries_at_count[sel])
    end
    parents
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
carrying the seed; a flip child whose rate law the guard rejects (`_is_well_posed`) is
dropped, as `expand_mechanisms` drops one.
"""
function _base_tier(mechs::Vector, rxn::EnzymeReaction)
    base = Union{Mechanism, AllostericMechanism}[m for m in mechs if !_degenerate(m)]
    failures = FitFailure[]
    for m in mechs
        _degenerate(m) || continue
        try
            children = filter!(_is_well_posed, _expand_re_to_ss(m))
            for child in children
                _assert_atom_conserving(child, rxn)
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
    # Keys selected so far per parameter count — the cumulative floor budget. Spent once
    # over the whole search, never re-granted per sweep.
    expanded_by_count = Dict{Int, Int}()
    # The keys selected so far (`_select_count!`), kept for the whole search: a late
    # member of a selected key expands without budget, and a late member of any other
    # key is judged again.
    selected_keys = Set{String}()
    # Fits by rate-law key, at canonical positions (`_Fit`; scaled to peak turnover 1 in
    # relative mode), shared across iterations so each distinct law is fit exactly once
    # over the whole search. LOOCV refits each key's source on its folds.
    memo = Dict{String,_Fit}()
    # Structures already produced — each is expanded at most once (termination).
    # This preserves the selected model: expansion is deterministic, so a
    # re-produced structure yields identical children (no new reachable
    # mechanism), and its FIRST production is its best beam-selection chance —
    # the width-floor budget (`expanded_by_count`) only shrinks, the loss cutoff
    # (`best_loss_by_count` running min) only tightens, and a selected key stays
    # selected (`selected_keys`), so a later re-production, which reads the loss its
    # first production read, cannot be selected where the first was not. An L phantom
    # first fitted under its own key is the exception: re-produced once its source's law
    # is fitted, it would read that fit's loss, which a later cutoff may admit. If
    # `_select_count!` ever gains a non-monotone budget, this invariant breaks.
    seen = Set{UInt64}()
    # Every mechanism_id the run's CSVs hold, with its mechanism's text
    # (`_write_rows_csv!`).
    texts = Dict{String,String}()

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
        _write_rows_csv!(save_dir, csv,
            vcat([e.row for e in entries], [_failure_row(f) for f in failures]), texts)
        improved = _ingest!(frontier, cv_pool, best_loss_by_count, entries;
                            n_cv_candidates)
        progress(string("  ", _postfit_summary(entries, failures),
                        "\n  ", _best_loss_line(best_loss_by_count, improved)))
        (entries, failures)
    end

    # The first log line names the package version that wrote the run's CSVs.
    progress("EnzymeRates v$(pkgversion(EnzymeRates))")

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
    base_batch === nothing && return (BatchEntry[], memo)
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
        to_expand = _sweep_parents!(frontier, target, expanded_by_count,
            best_loss_by_count, selected_keys; loss_rel_threshold, loss_abs_threshold,
            loss_parsimony_threshold, min_beam_width)
        isempty(to_expand) && continue

        # Expand each parent and record which parent produced each child
        # (first parent wins on structural dedup, matching `unique!`), so the
        # saved CSV can carry the parent's mechanism_id and parameter count for
        # diagnosing per-move parameter changes. The parent's `mechanism_id` is
        # already on its `BatchEntry.row`. Typed for dispatch: expand_mechanisms
        # needs a concrete Vector{<:Union{Mechanism, AllostericMechanism}} eltype.
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

    return BatchEntry[e for v in values(cv_pool) for e in v], memo
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

"""
    _cv_model_selection(pool, memo, prob; optimizer, se_threshold, save_dir,
                        show_progress, kwargs...) -> IdentifyRateEquationResults

Cross-validate the candidates of `pool` (the cv pool `_beam_search` returns) by leaving
one group of `prob.data` out at a time, select one by the 1-SE rule (`_select_best_row`)
and write `loocv_results.csv` and `best_equation.csv` to `save_dir`. A candidate's fold
losses come from fitting, on each fold, the source of the fit it reads
(`memo[e.fit_key].source`), so candidates that read one fit share their fold losses and
each fit's folds run once.
"""
function _cv_model_selection(
    pool::Vector{BatchEntry}, memo::Dict{String,_Fit},
    prob::IdentifyRateEquationProblem;
    optimizer, se_threshold::Float64,
    save_dir, show_progress,
    kwargs...
)
    isempty(pool) && error(
        "No mechanisms were successfully " *
        "fitted during beam search")

    # LOOCV candidates: the cv pool, ordered by `n_params`, then `loss`. `_offer_cv!`
    # keeps the top `n_cv_candidates` entries with DISTINCT rate-law keys (lowest `loss`
    # each) per `n_params` bucket, so one law never takes two LOOCV slots at a count and
    # folds are never wasted on, or biased by, its duplicates.
    candidates = pool[sortperm([(e.n_params, e.loss) for e in pool])]
    cv_df = _rows_to_dataframe([e.row for e in candidates])

    _progress(save_dir, show_progress,
        "Cross-validating $(length(candidates)) candidate equations (LOOCV)…")
    # The fold losses of each fit key, from fitting its source once per fold, flattened to
    # a (fit key, fold) grid so all folds run across every worker. `pmap` keeps the grid
    # order (the fold varies fastest), so column j of `by_fit` holds fit key j's fold
    # losses in `groups` order, and column ci of `scores` candidate ci's.
    groups = unique(prob.data.group)
    fit_keys = unique(e.fit_key for e in candidates)
    fold_losses = pmap([(memo[k].source, g) for k in fit_keys for g in groups]) do (m, g)
        _cv_fold_loss(compile_mechanism(m), prob, g; optimizer, kwargs...)
    end
    by_fit = reshape(Float64.(fold_losses), length(groups), :)
    scores = by_fit[:, [findfirst(==(e.fit_key), fit_keys) for e in candidates]]

    cv_df.cv_score = [mean(c) for c in eachcol(scores)]
    cv_df.cv_score_se = [std(c) / sqrt(length(groups)) for c in eachcol(scores)]

    best_row_idx = _select_best_row(cv_df; se_threshold)
    best_mechanism = compile_mechanism(candidates[best_row_idx].mech)

    # Flatten per-fold scores into one column per held-out group.
    # Group order matches the `groups = unique(prob.data.group)` iteration above.
    for (i, g) in enumerate(groups)
        cv_df[!, Symbol("cv_fold_$g")] = scores[i, :]
    end

    _progress(save_dir, show_progress,
        "Selected: $(nameof(typeof(best_mechanism))) " *
        "(key=$(cv_df.key[best_row_idx])), " *
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

This function and `_fold_problem` take `mechanism` `@nospecialize`. LOOCV runs a
candidate's folds after its full fit, which compiled `loss!` for the folds' problem type,
so a fold infers and compiles nothing for the mechanism.
"""
function _cv_fold_loss(
    @nospecialize(mechanism::AbstractEnzymeMechanism),
    prob::IdentifyRateEquationProblem, held_out;
    optimizer, kwargs...)
    held = prob.data.group .== held_out
    fit = fit_rate_equation(_fold_problem(mechanism, prob, .!held), optimizer; kwargs...)
    # `fit.params` is keyed by `fitted_params(mechanism)`, the order `loss!` reads.
    test_loss = loss!([log(v) for v in fit.params], _fold_problem(mechanism, prob, held))
    # A non-finite fold loss means the fit is unusable; aborting model
    # selection is correct (re-run CV from the saved CSVs after fixing
    # the fit).
    isfinite(test_loss) || error(
        "LOOCV produced a non-finite test loss for held-out group " *
        "$held_out — the fit is unusable; aborting model selection.")
    test_loss
end

"""
The `FittingProblem` of `mechanism` on the rows of `prob.data` that `mask` selects. The
fold copies its rows rather than viewing them, so its problem has the full fit's type and
`loss!` does not compile a second time for it; the copy is small next to a fit.
"""
function _fold_problem(@nospecialize(mechanism::AbstractEnzymeMechanism),
                       prob::IdentifyRateEquationProblem, mask)
    FittingProblem(mechanism, map(col -> col[mask], prob.data);
        Keq=prob.Keq, absolute=prob.absolute, peak_range=prob.peak_range)
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
