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
"""
struct IdentifyRateEquationProblem{
    R<:EnzymeReaction, D<:NamedTuple
}
    reaction::R
    data::D
    Keq::Float64
end

function IdentifyRateEquationProblem(
    reaction::EnzymeReaction, table; Keq::Real
)
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

    IdentifyRateEquationProblem{
        typeof(reaction),typeof(data)
    }(reaction, data, Float64(Keq))
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
    identify_rate_equation(prob; kwargs...)

Find the best rate equation for the given reaction
and data using beam search.

# Keyword Arguments
- `min_beam_width::Int = 50`: minimum mechanisms
  to keep per param-count tier
- `loss_rel_threshold::Float64 = 2.0`: relative tolerance
  for beam selection (see "Beam selection" below)
- `loss_abs_threshold::Float64 = 0.01`: absolute tolerance
  for beam selection
- `max_param_count::Int = 20`: stop expanding beyond
- `optimizer`: Optimization.jl optimizer (required).
  Recommended: `PyCMAOpt()` from OptimizationPyCMAES.
- `n_restarts::Int = 20`: multi-start restarts per fit
- `maxtime::Real = 60.0`: max time per fit (seconds)
- `maxiters::Int = 10_000_000`: max iterations per
  optimizer run (forwarded to `Optimization.solve`)
- `popsize::Int = 200`: population size for optimizer
  (forwarded to `Optimization.solve`)
- `verbose::Int = -9`: optimizer verbosity
  (forwarded to `Optimization.solve`)
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
  search CSVs (`initial_mechanisms.csv` + `equation_search_iteration_N.csv`)
- `pmap_function::Function = pmap`: parallelism
  function (Distributed.pmap by default)
- Extra kwargs are forwarded to `fit_rate_equation`
  and then to `Optimization.solve`.

# Beam selection

A mechanism qualifies for the next-level beam if either:
- its loss ≤ `loss_rel_threshold * best_loss + loss_abs_threshold`,
- OR its rank by loss (ascending) ≤ `min_beam_width`.

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
    max_param_count::Int = 20,
    # Fitting
    optimizer,
    n_restarts::Int = 20,
    maxtime::Real = 60.0,
    maxiters::Int = 10_000_000,
    popsize::Int = 200,
    verbose::Int = -9,
    # Model selection
    n_cv_candidates::Int = 5,
    se_threshold::Float64 = 1.0,
    perm_p_threshold::Float64 = 0.16,
    # Output & parallelism
    save_dir::String = _default_save_dir(),
    pmap_function::Function = pmap,
    # Extra fitting/optimizer kwargs
    optim_kwargs...
)
    fitting_kwargs = (;
        n_restarts, maxtime,
        maxiters, popsize, verbose,
        optim_kwargs...)

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
        loss_abs_threshold,
        max_param_count, save_dir,
        pmap_function, optimizer, n_cv_candidates,
        fitting_kwargs...)

    return _cv_model_selection(
        mechanisms, df, prob;
        n_cv_candidates, se_threshold, perm_p_threshold,
        pmap_function, optimizer, fitting_kwargs...)
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
        loss = [r.loss for r in rows],
        mechanism_type = [r.mechanism_type for r in rows],
        rate_equation = [r.rate_equation for r in rows],
        eq_hash = [r.eq_hash for r in rows],
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
  • its loss ≤ loss_rel_threshold * best_loss + loss_abs_threshold,
  • OR its rank (1-indexed by ascending loss) ≤ min_beam_width.

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
)
    finite_idx = [i for i in eachindex(losses) if isfinite(losses[i])]
    isempty(finite_idx) && return Int[]

    perm = sort(finite_idx; by=i -> losses[i])
    best = best_override === nothing ? losses[perm[1]] : best_override
    cutoff = loss_rel_threshold * best + loss_abs_threshold
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

"""One fitted mechanism: its own params + eq_hash + the CSV row."""
struct BatchEntry
    mech::Union{Mechanism, AllostericMechanism}
    n_params::Int
    loss::Float64
    eq_hash::UInt64
    row::NamedTuple
end

"""
Compile + cap-check + fit every mechanism in `mechs`, one `pmap` pass with
compile and fit fused on the same worker. Returns one `BatchEntry` per
fitted mechanism (each keeping its OWN fitted params and `eq_hash`). A
mechanism whose actual fitted-param count exceeds `max_param_count` is
dropped BEFORE fitting; compile/fit failures are dropped. No dedup here —
`mechs` is already structurally deduped by the caller (`_dedup_flat`).
"""
function _process_batch(
    mechs, prob::IdentifyRateEquationProblem;
    pmap_function, optimizer, max_param_count, kwargs...
)
    results = pmap_function(mechs) do m
        try
            em = compile_mechanism(m)
            fkeys = fitted_params(em)
            n = length(fkeys)
            n > max_param_count && return nothing
            eq_text = rate_equation_string(em)
            key = _rate_eq_dedup_key(eq_text)
            fp = FittingProblem(em, prob.data; Keq=prob.Keq)
            fit = fit_rate_equation(fp, optimizer; kwargs...)
            row = (
                n_params = n,
                loss = fit.loss,
                mechanism_type = string(typeof(em)),
                rate_equation = eq_text,
                fitted_param_names = fkeys,
                fitted_param_values =
                    Tuple(fit.params[k] for k in fkeys),
                eq_hash = string(key, base=16, pad=16),
            )
            BatchEntry(m, n, fit.loss, key, row)
        catch e
            @debug("_process_batch: compile or fit failed",
                   exception=(e, catch_backtrace()))
            nothing
        end
    end
    BatchEntry[r for r in results if r !== nothing]
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

function _beam_search(
    prob::IdentifyRateEquationProblem;
    min_beam_width, loss_rel_threshold, loss_abs_threshold,
    max_param_count, save_dir, pmap_function,
    optimizer, n_cv_candidates, kwargs...
)
    frontier = Dict{Int, Vector{BatchEntry}}()
    cv_pool  = Dict{Int, Vector{BatchEntry}}()
    best_loss_by_count = Dict{Int, Float64}()

    # ── Base tier: fit ALL init mechanisms (no bucketing — siblings) ──
    base = _dedup_flat(collect(init_mechanisms(prob.reaction)))
    base_entries = _process_batch(base, prob;
        pmap_function, optimizer, max_param_count, kwargs...)
    isempty(base_entries) && return (
        Union{Mechanism, AllostericMechanism}[],
        _rows_to_dataframe(NamedTuple[]))
    _save_initial_csv(save_dir, [e.row for e in base_entries])
    _ingest!(frontier, cv_pool, best_loss_by_count,
             base_entries; n_cv_candidates)

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
            sel = _select_beam([e.loss for e in entries_at_count];
                loss_rel_threshold, loss_abs_threshold,
                min_beam_width, best_override = best_loss_by_count[c])
            append!(to_expand, entries_at_count[sel])
        end

        if !isempty(to_expand)
            # Typed for dispatch: expand_mechanisms needs a concrete
            # Vector{<:Union{Mechanism, AllostericMechanism}} eltype.
            parents = Union{Mechanism, AllostericMechanism}[
                e.mech for e in to_expand]
            children = _dedup_flat(
                expand_mechanisms(parents, prob.reaction))
            child_entries = _process_batch(children, prob;
                pmap_function, optimizer, max_param_count, kwargs...)
            if !isempty(child_entries)
                # Count only iterations that fit children, so the
                # equation_search_iteration_N CSVs are gap-free.
                iteration += 1
                _save_iteration_csv(save_dir,
                    [e.row for e in child_entries], iteration)
                _ingest!(frontier, cv_pool, best_loss_by_count,
                         child_entries; n_cv_candidates)
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
params.
"""
function _evaluate_loss(
    mechanism, data, params, Keq
)
    pnames = fitted_params(mechanism)
    x = [log(params[p]) for p in pnames]
    fp = FittingProblem(mechanism, data; Keq=Keq)
    return loss!(x, fp)
end

"""
    _loocv(mechanism, prob; optimizer, kwargs...)

Leave-one-group-out cross-validation. Returns `Vector{Float64}`
of per-fold test losses (one per held-out group). Each score is
floored at `eps(Float64)` so `log(score)` is finite — the
selection rules in `_select_best_n_params` operate in log space.
On any internal failure (compile error, fit failure, etc.),
returns an empty `Float64[]`.
"""
function _loocv(
    mechanism::AbstractEnzymeMechanism,
    prob::IdentifyRateEquationProblem;
    optimizer, kwargs...
)
    groups = unique(prob.data.group)
    scores = Float64[]

    try
        for held_out in groups
            train_mask =
                prob.data.group .!= held_out
            test_mask =
                prob.data.group .== held_out

            train_data = _subset_data(
                prob.data, train_mask)
            test_data = _subset_data(
                prob.data, test_mask)

            fp_train = FittingProblem(
                mechanism, train_data;
                Keq=prob.Keq)
            fit = fit_rate_equation(
                fp_train, optimizer;
                kwargs...)

            test_loss = _evaluate_loss(
                mechanism, test_data,
                fit.params, prob.Keq)
            # Treat NaN or Inf as a fold-failure and abort the
            # entire LOOCV. `max(NaN, eps) === NaN` (NaN is not
            # ordered), so a non-finite test_loss would otherwise
            # bypass the floor and silently corrupt downstream
            # `cv_score = mean(v)` and `argmin(...)` selection.
            isfinite(test_loss) || return Float64[]
            # Floor at eps so log(score) is finite. The centered-
            # residuals loss can be exactly 0 (e.g. single-row
            # held-out group, or all residuals equal post-centering).
            push!(scores,
                  max(test_loss, eps(Float64)))
        end
    catch e
        @debug("LOOCV failed",
            exception=(e, catch_backtrace()))
        return Float64[]
    end

    scores
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

Errors if:
  * no bucket has any non-empty `cv_fold_scores` row;
  * any bucket's fold-score length differs from `n_min`'s.

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
    valid = filter(row -> !isempty(row.cv_fold_scores), cv_df)
    isempty(valid) && error(
        "no finite LOOCV scores in cv_df")
    sorted = sort(valid, [:n_params, :cv_score])
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
    n_cv_candidates, pmap_function, optimizer,
    se_threshold::Float64,
    perm_p_threshold::Float64,
    kwargs...
)
    isempty(mechs) && error(
        "No mechanisms were successfully " *
        "fitted during beam search")

    # Pick top n_cv_candidates per (n_params, eq_hash) bucket.
    candidate_indices = Int[]
    df_idx = DataFrame(
        row_idx = 1:nrow(df),
        n_params = df.n_params,
        loss = df.loss,
        eq_hash = df.eq_hash,
    )
    for gdf in groupby(df_idx, :n_params)
        seen_hashes = Set{String}()
        sorted = sort(gdf, :loss)
        for row in eachrow(sorted)
            row.eq_hash in seen_hashes && continue
            push!(seen_hashes, row.eq_hash)
            push!(candidate_indices, row.row_idx)
            length(seen_hashes) >= n_cv_candidates && break
        end
    end

    candidate_mechs = mechs[candidate_indices]
    candidate_rows = df[candidate_indices, :]

    fold_scores_per_candidate = pmap_function(
        candidate_mechs
    ) do mech
        m = compile_mechanism(mech)
        _loocv(m, prob; optimizer, kwargs...)
    end

    cv_df = copy(candidate_rows)
    cv_df.cv_fold_scores = collect(fold_scores_per_candidate)
    cv_df.cv_score = [isempty(v) ? Inf : mean(log.(v))
                      for v in cv_df.cv_fold_scores]

    all(!isfinite, cv_df.cv_score) && error(
        "All LOOCV scores are non-finite — every fold's fit " *
        "failed. The pipeline cannot select a best mechanism. " *
        "Inspect optimizer settings (n_restarts, maxtime), " *
        "data quality, or compile failures (run with " *
        "ENV[\"JULIA_DEBUG\"] = \"EnzymeRates\").")

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
    # Group order matches `_loocv`'s iteration over
    # `unique(prob.data.group)`.
    groups = unique(prob.data.group)
    for (i, g) in enumerate(groups)
        col = Symbol("cv_fold_$g")
        cv_df[!, col] = [
            isempty(v) ? missing : v[i]
            for v in cv_df.cv_fold_scores
        ]
    end

    select!(cv_df, Not(:cv_fold_scores))
    return IdentifyRateEquationResults(best_mechanism, cv_df)
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
