# ABOUTME: Beam-search pipeline to identify the best rate equation.
# ABOUTME: Enumerates mechanisms, fits each, selects via CV.

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
    # Extract metabolite names from reaction
    # (substrates/products are (name, atoms) pairs)
    mnames = tuple(
        [s[1] for s in substrates(reaction)]...,
        [p[1] for p in products(reaction)]...,
        regulators(reaction)...,
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

"""Cached fit result keyed by canonical rate-equation hash.
- `first_seen_estimate`: the beam-search level (the `pc` loop
  iteration value, equal to `n_fit_params_estimate`) at which
  this hash's fit was first performed.
- `first_seen_n_actual`: `length(fitted_params(m))` at first fit.
- `first_seen_eq_hash`: 16-char hex display string of the hash.
- `canon_to_rep`: pre-inverted `canonical_token => rep_orig_key`
  map, computed once at cache-insert. Spec members of the same
  hash group reuse this; avoids O(N) re-inversion per spec.
"""
struct _CachedFitResult
    loss::Float64
    params::NamedTuple
    canon_to_rep::Dict{String,String}
    first_seen_estimate::Int
    first_seen_n_actual::Int
    first_seen_eq_hash::String
end

"""Stage 1 result: uniform per-spec record so the `pmap` return
is concretely-typed. The `spec` field is `AbstractMechanismSpec`
(level vectors mix MechanismSpec and AllostericMechanismSpec, so
we can't tighten the type without splitting the pipeline). On
failure, every non-spec field has a sentinel value and `ok=false`."""
struct _Stage1Result
    spec::AbstractMechanismSpec
    eq_text::String
    h_full::UInt64
    h_short::String
    n_actual::Int
    mech_type_str::String
    name_map::Dict{String,String}
    fitted_keys::Tuple{Vararg{Symbol}}
    ok::Bool
end

"""Empty-failure sentinel."""
_Stage1Failure(spec::AbstractMechanismSpec) =
    _Stage1Result(
        spec, "", zero(UInt64), "", 0, "",
        Dict{String,String}(), (), false)

"""
Build the canonical text + name_map. Internal helper exposed
to `_canonical_rate_eq_hash_data` so callers can also retrieve
the name_map (needed by Stage 1 to project cached params across
specs in the same hash group).

Strategy: walk `parameters(m, Full)` to discover every parameter
symbol the mechanism could mention (including dependents that
appear in the v-line and on constraint LHSes), scan the
rate-equation body to find each parameter's first-appearance
position, then rename them as `p_1, p_2, …` in first-appearance
order. `:E_total` is in `Full` but is excluded from renaming;
`:Keq` and metabolite names are not in `Full` and aren't renamed.
The constraint lines are KEPT in the body — they encode
parameterization, so two mechanisms with the same v-line but
different choice of which parameter is dependent must hash
differently.

`parameters(m, Full)` is defined for both `EnzymeMechanism` and
`AllostericEnzymeMechanism` (Phase G.0). Allosteric coverage
includes T-state names, regulator-site names, and the allosteric
coupling `L` automatically.
"""
function _canonicalize_rate_eq_with_map(m::AbstractEnzymeMechanism)
    raw_body = rate_equation_string(m)
    raw_body === nothing && error(
        "rate_equation_string returned nothing for $(typeof(m))")
    body = String(raw_body)

    # Strip ONLY the destructure header lines.
    body = join(
        filter(
            ln -> !occursin(
                r"^\s*\(; .* = (params|concs)$", ln),
            split(body, '\n')),
        '\n')

    skip = (:E_total,)
    pnames = String[String(p) for p in parameters(m, Full)
                    if p ∉ skip]

    first_pos = Dict{String,Int}()
    for name in pnames
        rx = Regex("\\b" * name * "\\b")
        m_pos = match(rx, body)
        m_pos === nothing && continue
        first_pos[name] = m_pos.offset
    end
    appearing = collect(keys(first_pos))

    ordered = sort(appearing; by=name -> (first_pos[name], name))
    name_map = Dict(name => "p_$i"
                    for (i, name) in enumerate(ordered))

    # Substitute longest first to prevent prefix collisions
    # (e.g., rename `K1_T` before `K1`).
    for name in sort(appearing; by=length, rev=true)
        body = replace(body,
            Regex("\\b" * name * "\\b") => name_map[name])
    end

    canonical = strip(replace(body, r"\s+" => " "))
    (canonical, name_map)
end

"""
Return `(UInt64 hash, 16-char hex display string, name_map)`.
The single source for canonical hashing — both `_hash` and
`_hash_pair` delegate here so the canonicalizer runs once. Used
by Stage 1 of `_beam_search` to keep the rename mapping for later
per-spec param projection.

Hash collision probability over 10⁴ mechanisms is ~10⁻¹² with
Julia's built-in `hash(::String)::UInt64`.
"""
function _canonical_rate_eq_hash_data(m::AbstractEnzymeMechanism)
    canonical, name_map = _canonicalize_rate_eq_with_map(m)
    h = hash(canonical)
    (h, string(h, base=16, pad=16), name_map)
end

"""
Hash a mechanism's canonicalized rate equation. Returns the
`UInt64` hash.
"""
function _canonical_rate_eq_hash(m::AbstractEnzymeMechanism)
    first(_canonical_rate_eq_hash_data(m))
end

"""
Project cached params (keyed by rep spec's `fitted_params`
symbols) onto a target spec's own `fitted_params` keys, preserving
canonical-position values. Two specs in the same hash group have
isomorphic rate equations modulo parameter renaming; this function
applies the canonical position bijection
(rep_fitted_key → canonical_token → spec_fitted_key) to relabel
values without changing them.

`canon_to_rep` is the pre-inverted `canonical_token => rep_orig_key`
map (computed once at cache-insert from the rep's name_map).
`spec_name_map` is the spec's `orig_string => canonical_token` Dict
produced by the canonicalizer over `parameters(m, Full)`. They
include BOTH independent and dependent parameter names. We
restrict the projection to FITTED (independent) keys only —
`cached_params` is keyed by `fitted_params(rep_m)`, which doesn't
contain dep names. Iterating `keys(spec_name_map)` directly would
cause `KeyError` for any dep name (e.g., `:k10r`, `:K1_T` for
`:EqualRT` mirrors).

The return is a NamedTuple keyed by `fitted_params(spec_m)`.
"""
function _project_cached_params(
    cached_params::NamedTuple,
    canon_to_rep::Dict{String,String},
    spec_name_map::Dict{String,String},
    spec_fitted_keys::Tuple{Vararg{Symbol}},
)
    # Defensive lookup: a fitted key may not appear in the body
    # (e.g., a structurally-unidentifiable ghost param on a
    # zeroed `:NonequalRT` path), in which case `spec_name_map`
    # has no entry. Fall back to the spec key itself in cached_params
    # if both maps lack the canonical token; if even that misses,
    # use NaN as a sentinel that downstream loss/CV will surface.
    function _proj(k::Symbol)
        s = String(k)
        canon = get(spec_name_map, s, nothing)
        if canon !== nothing && haskey(canon_to_rep, canon)
            rep_key = Symbol(canon_to_rep[canon])
            haskey(cached_params, rep_key) &&
                return cached_params[rep_key]
        end
        haskey(cached_params, k) && return cached_params[k]
        return NaN
    end

    NamedTuple{spec_fitted_keys}(
        Tuple(_proj(k) for k in spec_fitted_keys))
end

"""
    identify_rate_equation(prob; kwargs...)

Find the best rate equation for the given reaction
and data using beam search.

# Keyword Arguments
- `min_beam_width::Int = 50`: minimum mechanisms
  to keep per level
- `loss_rel_threshold::Float64 = 2.0`: relative tolerance
  for beam selection (see "Beam selection" below)
- `loss_abs_threshold::Float64 = 0.01`: absolute tolerance
  for beam selection
- `max_param_count::Int = 20`: stop expanding beyond
- `optimizer`: Optimization.jl optimizer (required).
  Recommended: `PyCMAOpt()` from OptimizationPyCMAES.
- `n_restarts::Int = 10`: multi-start restarts per fit
- `maxtime::Real = 60.0`: max time per fit (seconds)
- `maxiters::Int = 10_000_000`: max iterations per
  optimizer run (forwarded to `Optimization.solve`)
- `popsize::Int = 200`: population size for optimizer
  (forwarded to `Optimization.solve`)
- `verbose::Int = -9`: optimizer verbosity
  (forwarded to `Optimization.solve`)
- `n_cv_candidates::Int = 5`: LOOCV top N
  **unique-rate-equation** candidates per param count
- `p_value_threshold::Float64 = 0.4`: parsimony threshold for
  the Wilcoxon signed-rank test in model selection. Smaller
  values demand stronger evidence to accept simpler models.
  Default 0.4 matches DataDrivenEnzymeRateEqs.jl convention
  (parsimony-permissive).
- `save_dir`: directory for per-level CSV files
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

The best `n_params` is chosen by the more parsimonious of two
methods, both operating on log of per-fold LOOCV scores:

1. **1-SE rule**: smallest `n_params` whose representative
   bucket-mean log-loss is within one standard error of the
   lowest representative bucket-mean log-loss.
2. **Wilcoxon signed-rank test**: smallest `n_params` whose
   representative per-fold log-losses are NOT statistically
   significantly worse than the best bucket
   (`pvalue > p_value_threshold`).

Per-bucket representative = the row with the lowest mean
fold-loss in that `n_params` bucket. Final pick is
`min(n_1se, n_wilcoxon)`. Within the chosen `n_params`, the
mechanism with lowest training loss wins.
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
    n_restarts::Int = 10,
    maxtime::Real = 60.0,
    maxiters::Int = 10_000_000,
    popsize::Int = 200,
    verbose::Int = -9,
    # Model selection
    n_cv_candidates::Int = 5,
    p_value_threshold::Float64 = 0.4,
    # Output & parallelism
    save_dir::Union{Nothing,String} = nothing,
    pmap_function::Function = pmap,
    # Extra fitting/optimizer kwargs
    optim_kwargs...
)
    fitting_kwargs = (;
        n_restarts, maxtime,
        maxiters, popsize, verbose,
        optim_kwargs...)

    if save_dir !== nothing && isdir(save_dir)
        existing = filter(
            f -> endswith(f, ".csv"),
            readdir(save_dir))
        isempty(existing) || error(
            "save_dir already contains CSV " *
            "files. Use an empty directory " *
            "to avoid mixing results.")
    end

    specs, df = _beam_search(prob;
        min_beam_width, loss_rel_threshold,
        loss_abs_threshold,
        max_param_count, save_dir,
        pmap_function, optimizer,
        fitting_kwargs...)

    return _cv_model_selection(
        specs, df, prob;
        n_cv_candidates, p_value_threshold,
        pmap_function, optimizer, fitting_kwargs...)
end

"""
Convert result row NamedTuples to a DataFrame.
Row order is preserved (no sorting) to maintain
alignment with the specs vector.
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
        fit_inherited_from_estimate = [
            r.fit_inherited_from_estimate for r in rows],
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
Save results for one beam level to a CSV file. The filename
encodes the level's `n_fit_params_estimate`; the actual `n_params`
of each row may be smaller (Haldane reduction collapses some
declared kinetic groups). Users wanting one file per actual
`n_params` value can post-process by reading and re-grouping.
"""
function _save_level_csv(
    save_dir::String, rows, n_fit_params_estimate::Int
)
    isdir(save_dir) || mkpath(save_dir)
    path = joinpath(
        save_dir,
        "params_estimate_$(n_fit_params_estimate).csv")
    df = _rows_to_dataframe(rows)
    CSV.write(path, df)
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
)
    finite_idx = [i for i in eachindex(losses) if isfinite(losses[i])]
    isempty(finite_idx) && return Int[]

    perm = sort(finite_idx; by=i -> losses[i])
    best = losses[perm[1]]
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

function _beam_search(
    prob::IdentifyRateEquationProblem;
    min_beam_width, loss_rel_threshold, loss_abs_threshold,
    max_param_count, save_dir, pmap_function,
    optimizer, kwargs...
)
    # Persistent cross-level cache keyed by canonical hash.
    fit_cache = Dict{UInt64, _CachedFitResult}()

    cache = Dict{Int,Vector{AbstractMechanismSpec}}()
    for spec in init_mechanisms(prob.reaction)
        push!(get!(cache, spec.n_fit_params_estimate,
                   AbstractMechanismSpec[]),
              spec)
    end
    dedup!(cache)

    all_specs = AbstractMechanismSpec[]
    all_rows  = NamedTuple[]

    isempty(cache) && return (
        all_specs, _rows_to_dataframe(all_rows))

    min_pc = minimum(keys(cache))
    for pc in min_pc:max_param_count
        level = pop!(cache, pc, AbstractMechanismSpec[])
        isempty(level) && (isempty(cache) ? break : continue)

        # ── Stage 1 (parallel): compile + hash ──
        compiled = pmap_function(level) do spec
            try
                m = compile_mechanism(spec)
                eq_text = rate_equation_string(m)
                h_full, h_short, name_map =
                    _canonical_rate_eq_hash_data(m)
                fkeys = fitted_params(m)
                n_actual = length(fkeys)
                mech_type_str = string(typeof(m))
                _Stage1Result(spec, eq_text, h_full, h_short,
                              n_actual, mech_type_str, name_map,
                              fkeys, true)
            catch e
                @debug("Mechanism compilation failed",
                       exception=(e, catch_backtrace()))
                _Stage1Failure(spec)
            end
        end
        filter!(c -> c.ok, compiled)
        isempty(compiled) && continue

        new_hashes = Set{UInt64}()
        for c in compiled
            haskey(fit_cache, c.h_full) && continue
            push!(new_hashes, c.h_full)
        end

        reps_by_hash = Dict{UInt64, _Stage1Result}()
        for c in compiled
            c.h_full in new_hashes || continue
            haskey(reps_by_hash, c.h_full) && continue
            reps_by_hash[c.h_full] = c
        end

        # ── Stage 2 (parallel): worker-side recompile + fit ──
        rep_results = pmap_function(
            collect(values(reps_by_hash))
        ) do rep
            try
                m = compile_mechanism(rep.spec)
                fp = FittingProblem(m, prob.data; Keq=prob.Keq)
                fit = fit_rate_equation(
                    fp, optimizer; kwargs...)
                (h_full=rep.h_full, h_short=rep.h_short,
                 n_actual=rep.n_actual,
                 name_map=rep.name_map,
                 loss=fit.loss, params=fit.params, ok=true)
            catch e
                @debug("Rep fit failed",
                       exception=(e, catch_backtrace()))
                (h_full=rep.h_full, ok=false)
            end
        end

        for r in rep_results
            r.ok || continue
            canon_to_rep = Dict(v => k for (k, v) in r.name_map)
            fit_cache[r.h_full] = _CachedFitResult(
                r.loss, r.params, canon_to_rep,
                pc, r.n_actual, r.h_short)
        end

        # ── Stage 3 (master): build ONE row per spec member ──
        # Use Stage 1's captured `fitted_keys` (computed once on
        # the worker that compiled) instead of recompiling on
        # master — saves a serial compile per spec.
        level_rows = NamedTuple[]
        level_specs = AbstractMechanismSpec[]
        for c in compiled
            haskey(fit_cache, c.h_full) || continue
            cached = fit_cache[c.h_full]
            is_inherited = !(c.h_full in new_hashes)
            spec_params = _project_cached_params(
                cached.params, cached.canon_to_rep,
                c.name_map, c.fitted_keys)
            row = (
                n_params = c.n_actual,
                loss = cached.loss,
                mechanism_type = c.mech_type_str,
                rate_equation = c.eq_text,
                fitted_param_names = c.fitted_keys,
                fitted_param_values =
                    Tuple(values(spec_params)),
                eq_hash = cached.first_seen_eq_hash,
                fit_inherited_from_estimate =
                    is_inherited ? cached.first_seen_estimate :
                                   missing,
            )
            push!(level_rows, row)
            push!(level_specs, c.spec)
        end

        append!(all_specs, level_specs)
        append!(all_rows,  level_rows)

        if save_dir !== nothing && !isempty(level_rows)
            _save_level_csv(save_dir, level_rows, pc)
        end

        sel = _select_beam(
            [r.loss for r in level_rows];
            loss_rel_threshold=loss_rel_threshold,
            loss_abs_threshold=loss_abs_threshold,
            min_beam_width=min_beam_width)
        beam_specs = level_specs[sel]

        new_cache = expand_mechanisms(beam_specs, prob.reaction)
        for (target_pc, specs) in new_cache
            target_pc > max_param_count && continue
            append!(get!(cache, target_pc,
                         AbstractMechanismSpec[]),
                    specs)
        end
        dedup!(cache)
    end

    df = _rows_to_dataframe(all_rows)
    return all_specs, df
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
Per-bucket setup shared by `_find_best_n_params_1se` and
`_find_best_n_params_wilcoxon`. Drops LOOCV-failure rows,
selects one representative row per `n_params` bucket (lowest
`cv_score` = lowest log-mean fold-loss), and builds two
parallel Dicts keyed by `n_params`:

- `log_means[n]` = `mean(log.(rep.cv_fold_scores))`. Recomputed
  from raw fold scores rather than read from `cv_score` so the
  helpers stay self-consistent if cv_score ever drifts from
  log-mean (e.g. test fixtures).
- `log_scores[n]` = `log.(rep.cv_fold_scores)`, used by Wilcoxon's
  paired signed-rank test and by 1-SE's standard error.

Also returns `n_min`, the bucket with the lowest log-mean.
Errors if no valid bucket remains. Per-bucket representatives
preserve the fold-pairing the Wilcoxon test relies on.
"""
function _per_bucket_log_stats(cv_df::DataFrame)
    valid = filter(
        row -> !isempty(row.cv_fold_scores), cv_df)
    isempty(valid) && error(
        "no finite LOOCV scores in cv_df")
    sorted = sort(valid, [:n_params, :cv_score])
    reps = combine(groupby(sorted, :n_params), first)
    log_scores = Dict(row.n_params => log.(row.cv_fold_scores)
                      for row in eachrow(reps))
    log_means = Dict(n => mean(ls)
                     for (n, ls) in log_scores)
    n_min = argmin(n -> log_means[n], keys(log_means))
    (log_means, log_scores, n_min)
end

"""
    _find_best_n_params_1se(cv_df) → Int

1-SE rule on log-transformed per-fold LOOCV scores. Returns the
smallest `n_params` whose representative log-mean is within one
standard error of the bucket with the lowest representative
log-mean. Standard error is `std(log_losses_at_min) / sqrt(n_folds)`.

If `n_folds == 1` for the best bucket (SE undefined), returns
`n_min` (no widening possible).
"""
function _find_best_n_params_1se(cv_df::DataFrame)
    log_means, log_scores, n_min =
        _per_bucket_log_stats(cv_df)
    losses_at_min = log_scores[n_min]
    n_folds = length(losses_at_min)
    n_folds == 1 && return n_min
    se = std(losses_at_min) / sqrt(n_folds)
    threshold = log_means[n_min] + se
    candidates = [n for n in keys(log_means)
                  if n <= n_min && log_means[n] <= threshold]
    minimum(candidates)
end

"""
    _find_best_n_params_wilcoxon(cv_df, p_threshold) → Int

Wilcoxon signed-rank rule on log-transformed per-fold LOOCV
scores. For each `n_params` bucket strictly below `n_min`, runs
a paired signed-rank test comparing that bucket's representative
per-fold log-losses to the best bucket's. Returns the smallest
`n_params` whose `pvalue > p_threshold` (NOT significantly
worse). Returns `n_min` if no smaller bucket qualifies.

`p_threshold = 0.4` is the parsimony-permissive default matching
`DataDrivenEnzymeRateEqs.jl`. Lower thresholds (e.g. 0.05)
require stronger evidence to accept simpler models.

Pairing semantics: the i-th element of each bucket's per-fold
score vector corresponds to the same held-out group, so
`pair(losses_smaller[i], losses_at_min[i])` is meaningful.
Per-bucket representatives ensure both vectors come from a
single mechanism. Skips comparisons where fold-counts differ.
"""
function _find_best_n_params_wilcoxon(
    cv_df::DataFrame, p_threshold::Float64,
)
    log_means, log_scores, n_min =
        _per_bucket_log_stats(cv_df)
    losses_at_min = log_scores[n_min]
    n_folds = length(losses_at_min)
    smaller_ns = sort([n for n in keys(log_means)
                       if n < n_min])
    for n in smaller_ns
        losses = log_scores[n]
        length(losses) == n_folds || continue
        try
            p = pvalue(ExactSignedRankTest(
                losses, losses_at_min))
            p > p_threshold && return n
        catch e
            @debug("Wilcoxon test failed for n_params=$n",
                   exception=(e, catch_backtrace()))
            continue
        end
    end
    n_min
end

"""
    _select_best_n_params(cv_df, p_threshold) → Int

Parsimony-aware model selection. Returns the more parsimonious
of the 1-SE rule and the Wilcoxon signed-rank rule. Both
operate on log-transformed per-fold LOOCV scores via per-bucket
representative selection (lowest cv_score row in each bucket).

Replaces strict-`argmin` over CV means: tiny CV improvements
from added parameters no longer justify the complexity bump.
"""
function _select_best_n_params(
    cv_df::DataFrame, p_threshold::Float64,
)
    n_se = _find_best_n_params_1se(cv_df)
    n_wx = _find_best_n_params_wilcoxon(
        cv_df, p_threshold)
    min(n_se, n_wx)
end

function _cv_model_selection(
    specs::Vector, df::DataFrame,
    prob::IdentifyRateEquationProblem;
    n_cv_candidates, pmap_function,
    optimizer, p_value_threshold, kwargs...
)
    # specs[i] corresponds to df[i, :] (same append
    # order in _beam_search, df is NOT sorted)
    isempty(specs) && error(
        "No mechanisms were successfully " *
        "fitted during beam search")
    df_indexed = copy(df)
    df_indexed.spec_idx = 1:nrow(df_indexed)

    candidate_indices = Int[]
    for gdf in groupby(df_indexed, :n_params)
        seen_hashes = Set{String}()
        sorted = sort(gdf, :loss)
        for row in eachrow(sorted)
            row.eq_hash in seen_hashes && continue
            push!(seen_hashes, row.eq_hash)
            push!(candidate_indices, row.spec_idx)
            length(seen_hashes) >= n_cv_candidates && break
        end
    end

    candidate_specs = specs[candidate_indices]
    candidate_rows = df[candidate_indices, :]

    # LOOCV each candidate in parallel — each result is the
    # candidate's Vector{Float64} of per-fold scores (or empty
    # on failure).
    fold_scores_per_candidate = pmap_function(
        candidate_specs
    ) do spec
        m = compile_mechanism(spec)
        _loocv(m, prob; optimizer, kwargs...)
    end

    # Build CV results DataFrame. `cv_score` holds the LOG-mean
    # of fold losses (the metric the 1-SE and Wilcoxon rules
    # operate on); using log-space everywhere avoids the
    # arithmetic-vs-log mismatch in rep selection.
    # `cv_fold_scores` holds the raw per-fold vector — INTERNAL
    # only; dropped before returning to the user so
    # `IdentifyRateEquationResults.cv_results` stays
    # CSV-serialisable.
    cv_df = copy(candidate_rows)
    cv_df.cv_fold_scores =
        collect(fold_scores_per_candidate)
    cv_df.cv_score = [isempty(v) ? Inf : mean(log.(v))
                      for v in cv_df.cv_fold_scores]
    cv_df.spec_idx = candidate_indices

    # Surface a hard error when every LOOCV fold failed —
    # otherwise `argmin([Inf, Inf, ...])` silently returns 1
    # and the user gets a meaningless "best" mechanism.
    all(!isfinite, cv_df.cv_score) && error(
        "All LOOCV scores are non-finite — every fold's fit " *
        "failed. The pipeline cannot select a best mechanism. " *
        "Inspect optimizer settings (n_restarts, maxtime), " *
        "data quality, or compile failures (run with " *
        "ENV[\"JULIA_DEBUG\"] = \"EnzymeRates\").")

    # Parsimony-aware selection: 1-SE rule + Wilcoxon
    # signed-rank test, both in log-loss space; take the more
    # parsimonious answer.
    best_param_count = _select_best_n_params(
        cv_df, p_value_threshold)

    # Best mechanism = lowest loss at that
    # param count
    at_best_pc = filter(
        row -> row.n_params == best_param_count,
        cv_df)
    sort!(at_best_pc, :loss)
    best_idx = at_best_pc[1, :spec_idx]
    best_mechanism = compile_mechanism(
        specs[best_idx])

    # Drop internal columns so user-facing cv_results stays
    # CSV-serialisable. cv_fold_scores is a Vector column —
    # CSV.write would stringify each cell.
    select!(cv_df, Not([:spec_idx, :cv_fold_scores]))
    return IdentifyRateEquationResults(
        best_mechanism, cv_df)
end
