# ABOUTME: identify_rate_equation, its problem and its results: the beam search over the
# ABOUTME: candidates, one fit per rate-law key, selection by cross-validation, the CSVs.

using DataFrames
using CSV
using Statistics
using SHA: sha256

"""
    IdentifyRateEquationProblem{D}

The reaction, rate table and equilibrium constant that `identify_rate_equation` finds a
mechanism for.

# Fields
- `reaction`: the `EnzymeReaction`
- `rate_data`: the checked rate table, as `FittingProblem` takes it, with `Keq`, `absolute`
  and `reference`, as a `RateData`

Construct it with `IdentifyRateEquationProblem(reaction, table, Keq; absolute = false,
reference = nothing)`. The table, the keywords and their errors are those of
`FittingProblem`, with one column per metabolite of `reaction`, anchors that name the
reaction's metabolites, and at least two `group` values, for cross-validation. Each fit
reads an anchor's concentrations of the metabolites its mechanism binds, so a mechanism
that leaves out an optional regulator ignores that regulator's concentration.
[Normalized vs absolute rate](@ref) explains the two modes.
"""
struct IdentifyRateEquationProblem{D<:NamedTuple}
    reaction::Mechanisms.EnzymeReaction
    rate_data::Data.RateData{D}
end

"""The distinct metabolite names of `rxn`: substrates, then products, then regulators,
each in stored order, keeping a name's first occurrence."""
_metabolite_names(rxn::Mechanisms.EnzymeReaction) =
    unique!(Symbol[[r.name for r in rxn.substrates]; [r.name for r in rxn.products];
                   [reg.name for reg in rxn.regulators]])

function IdentifyRateEquationProblem(
    reaction::Mechanisms.EnzymeReaction, table, Keq::Real;
    absolute::Bool=false, reference=nothing
)
    mnames = _metabolite_names(reaction)
    # Every metabolite the reaction declares needs a concentration column.
    rate_data = Data.RateData(table, mnames, Keq; absolute, reference, owner = "reaction")

    # Validate at least 2 groups for CV
    n_groups = length(unique(rate_data.data.group))
    n_groups >= 2 || error(
        "Need at least 2 unique groups for " *
        "cross-validation, got $n_groups")

    IdentifyRateEquationProblem(reaction, rate_data)
end

"""
The `FittingProblem` keywords of `prob` for `mechanism`: `absolute` and, in relative mode
with a reference, `prob.rate_data.reference` with each anchor's concentrations limited to
the metabolites `mechanism` binds (it does not respond to the others). The reference keeps
the problem's element type, so every fit and fold passes keywords of one type.
"""
function _fit_keywords(mechanism::Mechanisms.Mechanism, prob::IdentifyRateEquationProblem)
    (; absolute, reference) = prob.rate_data
    (absolute || isempty(reference)) && return (; absolute)
    mets = Mechanisms.metabolites(mechanism)
    (; absolute,
       reference = Pair{NamedTuple, Float64}[
           NamedTuple{Tuple(filter(in(mets), keys(concs)))}(concs) => turnover
           for (concs, turnover) in reference])
end

"""
    IdentifyRateEquationResults

Results from `identify_rate_equation`.

# Fields
- `best`: the selected mechanism, a `Mechanism`, chosen by the one-standard-error rule on
  the cross-validation scores ([Best mechanism selection](@ref))
- `cv_results`: `DataFrame` with LOOCV results for
  top candidates per param count; the same table as `loocv_results.csv`, so its
  `rate_equation` column holds the law without its coefficient values
"""
struct IdentifyRateEquationResults
    best::Mechanisms.Mechanism
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

Find the best rate equation for the reaction and data of `prob` and return an
`IdentifyRateEquationResults`. A beam search fits candidate mechanisms in order of
parameter count and expands the mechanisms whose rate-law key enters its beam. A key enters
when its loss is at most the loss cutoff or when it ranks within the `min_beam_width` floor.
The cutoff is the first of the two bounds below, or the smaller of the two once a smaller
parameter count has a fitted loss. Leave-one-group-out cross-validation then selects one
mechanism by the one-standard-error rule. [Best mechanism selection](@ref) states both
rules and the files the search writes to `save_dir`.

# Keyword Arguments
- `optimizer`: Optimization.jl optimizer (required).
  Recommended: `CMAEvolutionStrategyOpt()` from OptimizationCMAEvolutionStrategy.
- `min_beam_width::Int = 50`: the floor, the fewest rate-law keys the beam takes at each
  parameter count over the whole search, lowest loss first
- `loss_rel_threshold::Float64 = 1.3`: the first bound of the loss cutoff at a parameter
  count is this factor times the best loss there, plus `loss_abs_threshold`
- `loss_abs_threshold::Float64 = 0.001`: the additive term of the first bound
- `loss_parsimony_threshold::Float64 = 0.99`: the second bound is this factor times the
  best loss at any smaller parameter count; `Inf` disables it
- `max_param_count::Int = 20`: the search proposes no mechanism with more fitted constants
  than this (`parameters`)
- `eq_complexity_filter::Int = 337`: the most terms a proposed rate law may have
  (`EnzymeRates.Derivation.rate_equation_term_count`). The count is roughly the number of
  terms in the law's denominator, both conformations counted for an allosteric mechanism.
  The default admits a fully steady-state random-order bi-bi (336 terms); denser laws are
  slow to derive and fit. Counting costs no derivation.
- `optional_allosteric_regulators::Vector{Symbol} = Symbol[]`: declared
  allosteric regulators the seed may omit. A reaction that declares regulators
  seeds the beam fully regulated by default; naming a regulator here lets the
  beam also start from seeds that do not bind it (added back later by
  expansion). Listing every declared regulator restores the unregulated
  seeds.
- `optional_competitive_inhibitors::Vector{Symbol} = Symbol[]`: declared
  competitive inhibitors the seed may omit, as above.
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
- `n_cv_candidates::Int = 5`: how many distinct rate-law keys at each parameter count, those
  of lowest loss, enter cross-validation (at least 1)
- `se_threshold::Float64 = 1.0`: the multiplier of the one-standard-error rule; larger
  values favor simpler equations
- `save_dir::String`: the directory the search writes its CSVs to, by default a new
  `yyyy_mm_dd_results` directory in the working directory
- `show_progress::Bool = true`: print progress lines and append them to
  `save_dir/progress.log`

Errors when
- `n_cv_candidates` is less than 1;
- `save_dir` already holds a `.csv` file or a `progress.log`, so that the results of two
  runs never mix;
- the starting mechanisms cannot be built:
  - the reaction's metabolite names spell one string two ways (such as `Ac`, `CoA` and
    `AcCoA`);
  - the reaction has no catalytic cycle: every route needs an isomerization converting
    more than three substrates, or more than three products counting a covalent residue
    as one;
  - no mechanism binds every required regulator;
- no mechanism is fitted to a finite loss (when every starting mechanism throws, the error
  quotes the first failure);
- a cross-validation fold's loss is not finite;
- an internal check fails, such as two mechanisms whose texts share a `mechanism_id`.
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

function _write_csv(save_dir::String, filename::String, df::DataFrame)
    mkpath(save_dir)
    CSV.write(joinpath(save_dir, filename), df)
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
    _write_csv(save_dir, filename, _rows_to_dataframe(rows))
end

"""
The result rows `rows` as a `DataFrame`, in their order, with one column per fitted
parameter name across them (`missing` where a row lacks it).
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
        at_bound = Union{Missing,String}[get(r, :at_bound, missing) for r in rows],
    )
    for pn in sorted_pnames
        df[!, pn] = [get(r.params, pn, missing) for r in rows]
    end
    df
end

"""
The indices, in input order, of the count-`c` mechanisms with `losses` and rate-law keys
`law_keys` whose key is selected after this sweep. A key not yet in `selected` is judged at
its lowest loss by the beam rule ("Running the search" on the Best mechanism selection
page). It enters when that loss is within the loss cutoff, or when its rank among this
sweep's unselected keys is within the floor's budget left at `c`,
`min_beam_width - expanded[c]`. Adds the keys it selects to `selected` and their number to
`expanded[c]`. Every loss is finite: `_ingest!` keeps the others out of the frontier.
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

"""A mechanism's id in identify's CSVs: the first 16 hex digits of the SHA-256 of its text
(`show`), so every process and Julia version gives one mechanism one id."""
_mechanism_id(text::AbstractString) = bytes2hex(sha256(text))[1:16]

"""
A fit as identify stores it, once per rate-law key (`rate_equation_dedup_key`): `values`
holds the fitted constants at their canonical positions (its `parameters`), then `L` for an
MWC law, so every mechanism with that key reads them through its own canonical order;
`loss` and `retcode` are `fit_rate_equation`'s; `source` is the
mechanism that was fitted and `source_id` its `mechanism_id`; `at_bound` flags, at the
same positions, the constants `fit_rate_equation` reports at a bound.
"""
struct _Fit
    values::Vector{Float64}
    loss::Float64
    retcode::Symbol
    source::Mechanisms.Mechanism
    source_id::String
    at_bound::Vector{Bool}
end

"""A rate law's key and its fitted constants in canonical order
(`rate_equation_dedup_key`)."""
const _LawKey = @NamedTuple{key::String, order::Vector{Symbol}}

"""
What `_pass1_record` derived from a candidate (`candidate`): its mechanism's fitted
constants in `parameters` order, its text and the text's id (`mechanism_id`), its printed
law, its law's key (`canon`), whether it is an L phantom that reads another law's fit
(`rate_equation_dedup_key` gives it a fit key other than its key), `canon0`, the key of the
law whose fit it reads (for a phantom its active law, its law at `L = 0`; `canon`
otherwise), and whether its law's labelling reached the leaf cap (`capped`).
"""
struct _Pass1Record
    candidate::Enumeration.CandidateMechanism
    fitted::Tuple{Vararg{Symbol}}
    text::String
    mechanism_id::String
    eq_text::String
    canon::_LawKey
    phantom::Bool
    canon0::_LawKey
    capped::Bool
end

"""Record `r` as an ordinary one: not an L phantom, it reads the fit of its own law
(`canon0 = r.canon`)."""
_unflagged(r::_Pass1Record) = _Pass1Record(r.candidate, r.fitted, r.text, r.mechanism_id,
                                           r.eq_text, r.canon, false, r.canon, r.capped)

"""The order in which record `r` reads its fit's values (`_Fit`): the canonical order of
the law whose fit it reads (`canon0`), then `L` for an MWC law."""
_fit_order(r::_Pass1Record) =
    isempty(r.candidate.mechanism.allosteric.allosteric_states) ? r.canon0.order :
    [r.canon0.order; :L]

"""
The parameters of record `r` as its row holds them, in its `parameters` order: a fit's
`values`, stored at canonical positions, read through `order`, the record's own order for
those positions.
"""
function _row_params(r::_Pass1Record, order::Vector{Symbol}, values::Vector{Float64})
    at = Dict(zip(order, values))
    NamedTuple{r.fitted}(Tuple(at[x] for x in r.fitted))
end

"""One fitted candidate: its parameter count, its fit's loss and retcode, its rate-law key
(`rate_equation_dedup_key`), the key of the fit it reads (`fit_key`: its own, or for an L
phantom the key of its law at `L = 0`), and its CSV row."""
struct BatchEntry
    candidate::Enumeration.CandidateMechanism
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
    mech::Mechanisms.Mechanism
    error::String
end

_exc_string(@nospecialize(e)) = first(sprint(showerror, e), 200)

"""`pmap(f, itr)` as a `Vector{T}`. The element type of `pmap`'s result depends on the
values `f` returned, so inference cannot read it."""
_pmap_to_vector(f, ::Type{T}, itr) where {T} = convert(Vector{T}, pmap(f, itr))

"""
CSV row for a mechanism that threw. Same NamedTuple schema as a fitted row; it records
only the mechanism's text and id and the error (`missing` elsewhere, empty `params`).
`mechanism` is the mechanism's printed text (`show`), so a row names its mechanism even
when it does not compile; `parent_id` is missing.
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
     at_bound = missing)
end

"""
Print `msg` to stdout and append it to `<save_dir>/progress.log` when `show_progress`
holds. The flush makes the line appear at once in a redirected cluster job log: redirected
stdout is block-buffered and would withhold it until the process exits.
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
Pre-fit half of a batch summary: the three buckets known before fitting — `new fits`
(distinct new laws, one per key, being fit this batch), `inherited` (rows reusing a
memoized or in-batch fit), and `skipped (already fit)`.
"""
_prefit_summary(n_new::Int, n_inherited::Int, n_fitted_skip::Int) =
    string(n_new, " new fits + ", n_inherited, " inherited + ",
           n_fitted_skip, " skipped (already fit)")

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
What pass 1 of `_process_batch` derives from candidate `c` on a worker, from its mechanism
alone and compiling nothing: its fitted constants (`parameters`, which runs the guard
first), its law printed in named coefficients without their values
(`rate_equation_string`), its text and the text's id (`mechanism_id`), and its law's key
and the key of the fit it reads (`rate_equation_dedup_key`; a law whose labelling reaches
the leaf cap gets a key of its own, and an L phantom reads the fit of its law at
`L = 0`). Returns a `FitFailure` when a step throws, else a `_Pass1Record`.
"""
function _pass1_record(c::Enumeration.CandidateMechanism)
    m = c.mechanism
    try
        fitted = Derivation.parameters(m)
        eq_text = Derivation.rate_equation_string(m; coefficients = true,
                                                    coefficient_values = false)
        k = Derivation.rate_equation_dedup_key(m)
        canon = (key = k.key, order = k.parameters)
        phantom = k.fit_key != k.key
        canon0 = phantom ? (key = k.fit_key, order = k.fit_parameters) : canon
        text = sprint(show, m)
        _Pass1Record(c, fitted, text, _mechanism_id(text), eq_text, canon, phantom,
                     canon0, k.capped)
    catch e
        FitFailure(m, _exc_string(e))
    end
end

"""
Fit one batch of candidates and return `(entries, failures, n_seen_skip)`.

Pass 1 skips each candidate whose structure (`hash`) `seen` holds, counting it in
`n_seen_skip`, adds the other structures to `seen` and derives their candidates on the
workers (`_pass1_record`). Pass 2 fits one mechanism per rate-law key not in `memo`, on the
workers, and stores the fit in `memo` (`_Fit`). Every mechanism with that key, in this
batch or a later one, reads that fit under its own names. A fit that throws fails every
mechanism with its key. An L phantom is a mechanism to which `rate_equation_dedup_key`
gives a fit key other than its key. It reads, from `memo` or from this batch, the fit of
the law it has at `L = 0`, and its row adds `L = 0` to that fit's values. When neither
holds that fit, or that fit throws, the phantom is fitted under its own key.

`log` receives the pre-fit summary (`_prefit_summary`) with a line per mechanism whose law
reached the leaf cap, then a line per phantom whose source's fit threw. `parent_of` maps a
child to its parent's `BatchEntry`, which fills the row's parent columns. `entries` holds a
`BatchEntry` per candidate with a fit, each with its own row, `retcode` and key;
`failures` holds a `FitFailure`, with the exception text, per candidate whose derivation or
fit threw.
"""
function _process_batch(
    candidates::Vector{Enumeration.CandidateMechanism}, prob::IdentifyRateEquationProblem;
    optimizer, memo::Dict{String,_Fit}=Dict{String,_Fit}(),
    seen::Set{UInt64}=Set{UInt64}(),
    parent_of::Dict{Enumeration.CandidateMechanism, BatchEntry} =
        Dict{Enumeration.CandidateMechanism, BatchEntry}(),
    log = msg -> nothing, kwargs...
)
    fresh = Enumeration.CandidateMechanism[]
    for c in candidates
        h = hash(c)
        h in seen || (push!(seen, h); push!(fresh, c))
    end

    records = _pmap_to_vector(_pass1_record, Union{FitFailure, _Pass1Record}, fresh)

    # Pick one representative per key not already fit; an L phantom is none.
    rep_idx = Dict{String,Int}()
    for (i, r) in enumerate(records)
        r isa _Pass1Record && !r.phantom && !haskey(memo, r.canon.key) &&
            get!(rep_idx, r.canon.key, i)
    end
    # A phantom inherits the fit of whoever holds its L → 0 key, in `memo` or among this
    # batch's representatives; any other phantom becomes an ordinary record.
    for (i, r) in enumerate(records)
        r isa _Pass1Record || continue
        r.phantom || continue
        (haskey(memo, r.canon0.key) || haskey(rep_idx, r.canon0.key)) && continue
        records[i] = _unflagged(r)
        haskey(memo, r.canon.key) || get!(rep_idx, r.canon.key, i)
    end
    reps = _Pass1Record[records[i] for i in values(rep_idx)]
    n_seen_skip = length(candidates) - length(fresh)
    capped = ["leaf cap reached for $(r.mechanism_id); fitted under its own key"
              for r in records if r isa _Pass1Record && r.capped]
    log(join([_prefit_summary(length(reps),
                  count(r -> r isa _Pass1Record, records) - length(reps), n_seen_skip);
              capped], "\n  "))

    # Pass 2 runs across the workers, since fitting dominates the cost. Each
    # representative gives its fit at canonical positions, or the text of the error its
    # fit threw.
    fit_error = Dict{String,String}()   # key → its representative's fit threw
    function fit!(reps::Vector{_Pass1Record})
        rep_fits = _pmap_to_vector(Union{_Fit, String}, reps) do r
            m = r.candidate.mechanism
            try
                fp = FittingProblem(m, prob.rate_data.data, prob.rate_data.Keq;
                                    _fit_keywords(m, prob)...)
                fit = fit_rate_equation(fp, optimizer; kwargs...)
                order = _fit_order(r)
                _Fit(Float64[fit.params[x] for x in order], fit.loss, fit.retcode,
                     m, r.mechanism_id, Bool[x in fit.at_bound for x in order])
            catch e
                _exc_string(e)
            end
        end
        for (r, f) in zip(reps, rep_fits)
            f isa _Fit ? (memo[r.canon.key] = f) : (fit_error[r.canon.key] = f)
        end
    end
    fit!(reps)
    # A phantom whose source's fit threw becomes an ordinary record, fitted under its own
    # key: one representative per such key not already fit.
    empty!(rep_idx)
    threw = String[]
    for (i, r) in enumerate(records)
        r isa _Pass1Record || continue
        r.phantom && haskey(fit_error, r.canon0.key) || continue
        records[i] = _unflagged(r)
        haskey(memo, r.canon.key) || get!(rep_idx, r.canon.key, i)
        push!(threw, "source fit threw for L phantom $(r.mechanism_id); " *
                     "it reads its own key's fit")
    end
    isempty(threw) || log(join(threw, "\n  "))
    fit!(_Pass1Record[records[i] for i in values(rep_idx)])

    # Build a row per derived candidate from its key's fit, read through its own order.
    entries  = BatchEntry[]
    failures = FitFailure[]
    for r in records
        r isa FitFailure && (push!(failures, r); continue)
        fit_key = r.canon0.key                    # its own key, or a phantom's L → 0 key
        if haskey(fit_error, fit_key)             # the fit's source threw
            push!(failures, FitFailure(r.candidate.mechanism, fit_error[fit_key]))
            continue
        end
        fit = memo[fit_key]
        order = _fit_order(r)
        params = _row_params(r, order, r.phantom ? [fit.values; 0.0] : fit.values)
        flags = r.phantom ? [fit.at_bound; false] : fit.at_bound
        parent = get(parent_of, r.candidate, nothing)
        row = (
            n_params = length(r.fitted),
            parent_n_params = parent === nothing ? missing : parent.n_params,
            loss = fit.loss,
            mechanism = r.text,
            mechanism_id = r.mechanism_id,
            parent_id = parent === nothing ? missing : parent.row.mechanism_id::String,
            rate_equation = r.eq_text,
            retcode = string(fit.retcode),
            error = missing,
            params,
            key = r.canon.key,
            fit_inherited = fit.source_id != r.mechanism_id,
            fit_source_id = fit.source_id,
            L_phantom = r.phantom,
            at_bound = join(order[flags], ";"),
        )
        push!(entries, BatchEntry(r.candidate, length(r.fitted), fit.loss, fit.retcode,
                                  r.canon.key, fit_key, row))
    end
    (entries, failures, n_seen_skip)
end

"""
Fold the batch `entries` into the search state and return the set of counts whose best loss
dropped or first appeared. Each entry with a finite loss joins the `frontier` at its count,
even when another entry has its key; it lowers `best_loss_by_count` there and is offered to
`cv_pool` (`_offer_cv!`). An entry whose loss is not finite (`fit_rate_equation` returns
`Inf` when no restart is finite or the fit cannot match the reference turnovers) is skipped:
it is never expanded, never sets a best loss and never enters cross-validation; its row is
already in the batch's CSV.
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
Expand every selected parent into its children (`Enumeration.candidate_mechanisms`, with
the size limits `max_param_count` and `eq_complexity_filter`) across the workers. Returns
`(children, parent_of, expand_failures)`, in the order of `to_expand` whatever the number
of workers:
- `children`: each child once;
- `parent_of`: each child's first parent, as its `BatchEntry`;
- `expand_failures`: a `FitFailure` carrying the parent's mechanism for each parent whose
  expansion threw.
"""
function _expand_parents(to_expand::Vector{BatchEntry},
                         reaction::Mechanisms.EnzymeReaction;
                         max_param_count::Int, eq_complexity_filter::Int)
    results = _pmap_to_vector(
            Tuple{Vector{Enumeration.CandidateMechanism}, Union{Nothing, FitFailure}},
            [pe.candidate for pe in to_expand]) do parent
        # Record a per-parent expansion error (e.g. the atom-conservation assertion) as a
        # failure; never abort the search.
        try
            (Enumeration.candidate_mechanisms(reaction, parent; max_param_count,
                                              eq_complexity_filter), nothing)
        catch e
            (Enumeration.CandidateMechanism[], FitFailure(parent.mechanism, _exc_string(e)))
        end
    end
    parent_of = Dict{Enumeration.CandidateMechanism, BatchEntry}()
    children = Enumeration.CandidateMechanism[]
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
    # Fits by rate-law key, at canonical positions (`_Fit`; scaled to the reference
    # turnovers in relative mode), shared across iterations so each distinct law is fit
    # exactly once over the whole search. LOOCV refits each key's source on its folds.
    memo = Dict{String,_Fit}()
    # Structures already produced; each is expanded at most once, which ends the search.
    # This loses no model. Expansion is deterministic, so a structure produced again has the
    # same children and reads the same loss. Its first production had its best chance at
    # the beam: the width floor's budget (`expanded_by_count`) only shrinks, the loss cutoff
    # (`best_loss_by_count`) only tightens, and a selected key stays selected
    # (`selected_keys`). An L phantom first fitted under its own key is the exception:
    # produced again once its source's law is fitted, it would read that fit's loss, which a
    # later cutoff may admit. A budget that could grow back in `_select_count!` would break
    # this argument.
    seen = Set{UInt64}()
    # Every mechanism_id the run's CSVs hold, with its mechanism's text
    # (`_write_rows_csv!`).
    texts = Dict{String,String}()

    progress(msg) = _progress(save_dir, show_progress, msg)
    # Fit one batch and record it: log `header` with the pre-fit summary, write the
    # rows to `csv`, fold the entries into the search state, and log the post-fit
    # summary. Returns `(entries, failures)`, or `nothing`, writing no CSV, when the
    # batch has no rows.
    function run_batch!(candidates, extra_failures, header, csv;
                        parent_of = Dict{Enumeration.CandidateMechanism, BatchEntry}())
        entries, failures = _process_batch(candidates, prob; optimizer, memo, seen,
            parent_of, log = msg -> progress(header * "\n  " * msg), kwargs...)
        append!(failures, extra_failures)
        isempty(entries) && isempty(failures) && return nothing
        _write_rows_csv!(save_dir, csv,
            vcat(NamedTuple[e.row for e in entries],
                 NamedTuple[_failure_row(f) for f in failures]), texts)
        improved = _ingest!(frontier, cv_pool, best_loss_by_count, entries;
                            n_cv_candidates)
        progress(string("  ", _postfit_summary(entries, failures),
                        "\n  ", _best_loss_line(best_loss_by_count, improved)))
        (entries, failures)
    end

    # The first log line names the package version that wrote the run's CSVs.
    progress("EnzymeRates v$(pkgversion(EnzymeRates))")

    # ── Base tier: fit every starting candidate (no bucketing — siblings) ──
    progress("Enumerating initial mechanisms…")
    base = Enumeration.candidate_mechanisms(prob.reaction; max_param_count,
        eq_complexity_filter, optional_allosteric_regulators,
        optional_competitive_inhibitors)
    base_batch = run_batch!(base, FitFailure[],
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

        children, parent_of, expand_failures = _expand_parents(to_expand, prob.reaction;
            max_param_count, eq_complexity_filter)
        # `iteration` numbers the batches that wrote rows, so the iteration files have no
        # gaps; it is not a parameter count.
        run_batch!(children, expand_failures,
            "Iteration $(iteration + 1): $(length(to_expand)) parents → " *
            "$(length(children)) children",
            "equation_search_iteration_$(iteration + 1).csv"; parent_of) === nothing ||
            (iteration += 1)
    end

    return BatchEntry[e for v in values(cv_pool) for e in v], memo
end

"""The index of the row of `cv_df` that the one-standard-error rule selects ("The
cross-validation selection rule" on the Best mechanism selection page)."""
function _select_best_row(cv_df::DataFrame; se_threshold::Float64 = 1.0)::Int
    by_score = sortperm(collect(zip(cv_df.cv_score, cv_df.n_params)))
    i_min = by_score[1]
    cutoff = cv_df.cv_score[i_min] + se_threshold * cv_df.cv_score_se[i_min]
    passing = [i for i in by_score if cv_df.cv_score[i] <= cutoff]
    best_n = minimum(cv_df.n_params[passing])
    return first(i for i in passing if cv_df.n_params[i] == best_n)
end

"""
    _fold_problem(mechanism, data, mask, Keq; keywords...) → FittingProblem

The `FittingProblem` of `mechanism` on the rows of `data` that `mask` selects, with `Keq`
and the `FittingProblem` `keywords`. The fold copies its rows rather than viewing them, so
its problem has the full fit's type and `loss!` does not compile a second time for it; the
copy is small next to a fit.
"""
function _fold_problem(mechanism::Mechanisms.Mechanism, data::NamedTuple, mask, Keq::Real;
                       keywords...)
    FittingProblem(mechanism, map(col -> col[mask], data), Keq; keywords...)
end

_single_row_dataframe(df::DataFrame, i::Int) = df[[i], :]

"""
    _cross_validation_loss(mechanism, data, held_out, Keq; absolute, reference = nothing,
        optimizer, kwargs...) → Float64

The loss on the rows of `data` whose `group` is `held_out` of `mechanism` fitted to the
other rows: `fit_rate_equation(problem, optimizer; kwargs...)` on a `FittingProblem` with
`Keq`, `absolute` and `reference`, which absolute mode leaves `nothing`, scored by
`test_rate_equation` with `Keq`. A training fit that fails its reference
(`:ReferenceFailed`) keeps its parameters unscaled, and the relative loss, blind to the
scale, still scores them. A non-finite loss raises (naming the held-out group): a
corrupted fold must abort model selection rather than propagate a bad score.
"""
function _cross_validation_loss(
    mechanism::Mechanisms.Mechanism, data::NamedTuple, held_out, Keq::Real;
    absolute, reference = nothing, optimizer, kwargs...)
    held = isequal.(data.group, held_out)
    fit = fit_rate_equation(
        _fold_problem(mechanism, data, .!held, Keq; absolute, reference),
        optimizer; kwargs...)
    test_loss = test_rate_equation(mechanism, fit.params, map(col -> col[held], data),
                                   Keq; absolute)
    isfinite(test_loss) || error(
        "LOOCV produced a non-finite test loss for held-out group " *
        "$held_out — the fit is unusable; aborting model selection.")
    test_loss
end

"""
    _cv_model_selection(pool, memo, prob; optimizer, se_threshold, save_dir,
                        show_progress, kwargs...) -> IdentifyRateEquationResults

Cross-validate the candidates of `pool` (the cv pool `_beam_search` returns) by leaving
one group of `prob.rate_data.data` out at a time, select one by the 1-SE rule
(`_select_best_row`) and write `loocv_results.csv` and `best_equation.csv` to
`save_dir`. A candidate's fold losses come from fitting, on each fold, the source of the
fit it reads (`memo[e.fit_key].source`), so candidates that read one fit share their fold
losses and each fit's folds run once.
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

    # The cv pool by `n_params`, then `loss`. It holds one entry per key (`_offer_cv!`),
    # so no law takes two cross-validation slots at a count or weighs twice in the folds.
    candidates = pool[sortperm([(e.n_params, e.loss) for e in pool])]
    cv_df = _rows_to_dataframe(NamedTuple[e.row for e in candidates])

    _progress(save_dir, show_progress,
        "Cross-validating $(length(candidates)) candidate equations (LOOCV)…")
    # The fold losses of each fit key, from fitting its source once per fold, flattened to
    # a (fit key, fold) grid so all folds run across every worker. `pmap` keeps the grid
    # order (the fold varies fastest), so column j of `by_fit` holds fit key j's fold
    # losses in `groups` order, and column ci of `scores` candidate ci's.
    groups = unique(prob.rate_data.data.group)
    fit_keys = unique(e.fit_key for e in candidates)
    folds = Tuple{Mechanisms.Mechanism, eltype(groups)}[(memo[k].source, g)
                                            for k in fit_keys for g in groups]
    fold_losses = _pmap_to_vector(Float64, folds) do (m, g)
        _cross_validation_loss(m, prob.rate_data.data, g, prob.rate_data.Keq;
            _fit_keywords(m, prob)..., optimizer, kwargs...)
    end
    by_fit = reshape(fold_losses, length(groups), :)
    scores = by_fit[:, [findfirst(==(e.fit_key), fit_keys)::Int for e in candidates]]

    cv_df.cv_score = [mean(c) for c in eachcol(scores)]
    cv_df.cv_score_se = [std(c) / sqrt(length(groups)) for c in eachcol(scores)]

    best_row_idx = _select_best_row(cv_df; se_threshold)
    best = candidates[best_row_idx]

    for (i, g) in enumerate(groups)
        cv_df[!, Symbol("cv_fold_$g")] = scores[i, :]
    end

    _progress(save_dir, show_progress,
        "Selected: key=$(best.key), n_params=$(best.n_params)")

    # Save the LOOCV table and the selected best equation alongside the
    # per-iteration fit CSVs, so cluster runs persist the model-selection
    # outcome (recoverable from the iteration CSVs, but wasteful to omit).
    _write_csv(save_dir, "loocv_results.csv", cv_df)
    best_df = _single_row_dataframe(cv_df, best_row_idx)
    best_df.rate_equation = [Derivation.rate_equation_string(
        best.candidate.mechanism; coefficients = true)]
    _write_csv(save_dir, "best_equation.csv", best_df)

    return IdentifyRateEquationResults(best.candidate.mechanism, cv_df)
end

"""
Default results directory: the first non-existent `<date>_results[_N]`
directory in the cwd (e.g. `YYYY_MM_DD_results`, then `…_results_2`, `_3`).
"""
function _default_save_dir()
    base = string(Dates.format(Dates.today(), dateformat"yyyy_mm_dd"), "_results")
    isdir(base) || return base
    n = 2
    while isdir(string(base, "_", n)); n += 1; end
    string(base, "_", n)
end
