# Log the fit plan before fitting — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Print each search batch's parent/child/new/inherited/skip counts before
fitting, and the errored/success counts after.

**Architecture:** Split `_process_batch` at its existing PASS-1/PASS-2 seam into
`_compile_batch` and `_fit_batch`; `_process_batch` stays as a thin wrapper so
its tests are unchanged. Replace `_batch_summary` with `_prefit_summary` (the
five pre-fit buckets) and `_postfit_summary` (errored + percentages). Both
`_beam_search` call sites print the pre-fit line, fit, then print the post-fit
line.

**Tech Stack:** Julia; existing test suite (`Pkg.test()`).

## Global Constraints

- Changes only what prints and when — never what is enumerated or fit.
- `_process_batch` keeps its signature and 5-tuple return `(entries, failures,
  n_param_skip, n_cx_skip, n_fitted_skip)`; its tests
  (`test/test_identify_rate_equation.jl:1123-1191`) must pass unchanged.
- `new fits` printed pre-fit = the distinct new equations being fit this batch
  (the fit representatives). `errored` printed post-fit = `length(failures)`.
  `new fits = successful + errored`.
- `rate_equation` untouched — its 0-alloc / sub-120 ns gate stays green.
- 92-char line limit, 4-space indent, match surrounding style.
- Do NOT run the full `Pkg.test()` inside a subagent — it exceeds the Bash cap
  and orphans. Verify with a scratch driver (below); the controller runs the
  full suite as the gate.

Scratch-driver pattern (fast, foreground, well under the cap):
```
julia --project /tmp/.../scratchpad/verify.jl
# verify.jl: using EnzymeRates, Test, DataFrames; @testset ... end
```

---

### Task 1: Split `_process_batch`; add the two summaries

**Files:**
- Modify: `src/identify_rate_equation.jl` — split `_process_batch:541-644`; replace
  `_batch_summary:474-501`.
- Test: `test/test_identify_rate_equation.jl:1089-1108` (rewrite the summary test).

**Interfaces:**
- Produces: `_compile_batch(mechs, prob; max_param_count, eq_complexity_filter,
  memo, fitted) -> (compiled, reps, rep_idx, n_fitted_skip)`;
  `_fit_batch(compiled, reps, rep_idx, prob, memo; optimizer, parent_of, kwargs...)
  -> (entries, failures)`;
  `_prefit_summary(n_new, n_inherited, n_param_skip, n_cx_skip, n_fitted_skip;
  max_param_count, eq_complexity_filter) -> String`;
  `_postfit_summary(entries, failures) -> String`.
- Unchanged: `_process_batch(...)` still returns the 5-tuple.

- [ ] **Step 1: Rewrite the summary test (RED)**

Replace the `_batch_summary` block at `test/test_identify_rate_equation.jl:1089-1108`
with:

```julia
    # _prefit_summary: the five pre-fit buckets, no errored, no percentages.
    pre = EnzymeRates._prefit_summary(2, 0, 4, 2, 3;
        max_param_count=8, eq_complexity_filter=337)
    @test occursin("2 new fits + 0 inherited + 3 skipped (already fit) + " *
                   "4 skipped (>8 params) + 2 skipped (>337 complexity)", pre)
    @test !occursin("errored", pre)
    @test !occursin("Success", pre)

    # _postfit_summary: errored + success/non-Success over the fitted set.
    mech = first(EnzymeRates.init_mechanisms(@enzyme_reaction begin
        substrates: S[C]; products: P[C] end))
    row = (n_params=3, loss=0.5, mechanism_type="M", rate_equation="v",
           retcode="Success", error=missing, fitted_param_names=(:K,),
           fitted_param_values=(1.0,), eq_hash="abc", fit_inherited=false)
    e_succ = EnzymeRates.BatchEntry(mech, 3, 0.5, :Success, hash(:a), row)
    e_mt   = EnzymeRates.BatchEntry(mech, 3, 0.9, :MaxTime, hash(:b), row)
    f      = EnzymeRates.FitFailure(mech, "StackOverflowError: ")
    post = EnzymeRates._postfit_summary([e_succ, e_mt], [f])
    @test occursin("1 errored", post)
    @test occursin("Success 50.0%", post)              # 1 of 2 fitted
    @test occursin("non-Success retcode 50.0%", post)  # e_mt is :MaxTime
    @test !occursin("best loss", post)
```

Copy this same block into a scratch driver (`using EnzymeRates, Test, DataFrames`
+ the block wrapped in a `@testset`). Run it → RED: `UndefVarError: _prefit_summary`.

- [ ] **Step 2: Add the two summary functions**

Replace `_batch_summary` (`src/identify_rate_equation.jl:474-501`, the docstring
and the whole function) with:

```julia
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
```

- [ ] **Step 3: Split `_process_batch` into `_compile_batch` + `_fit_batch` + wrapper**

Replace the `_process_batch` function body (`src/identify_rate_equation.jl:541-644`,
keeping its docstring above) with the three functions below. `_compile_batch` is
lines 550-591 of the original verbatim; `_fit_batch` is lines 592-639 verbatim;
the wrapper reproduces the original 5-tuple return.

```julia
function _compile_batch(
    mechs, prob::IdentifyRateEquationProblem;
    max_param_count, eq_complexity_filter::Int = typemax(Int),
    memo::Dict{UInt64,NamedTuple}=Dict{UInt64,NamedTuple}(),
    fitted::Set{UInt64}=Set{UInt64}(),
)
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
    rep_idx = Dict{UInt64,Int}()
    for (i, c) in enumerate(compiled)
        c isa NamedTuple || continue
        (haskey(memo, c.eq_hash) || haskey(rep_idx, c.eq_hash)) && continue
        rep_idx[c.eq_hash] = i
    end
    reps = [(mech = compiled[i].mech, eq_hash = compiled[i].eq_hash)
            for i in values(rep_idx)]
    (compiled, reps, rep_idx, n_fitted_skip)
end

function _fit_batch(compiled, reps, rep_idx::Dict{UInt64,Int},
    prob::IdentifyRateEquationProblem, memo::Dict{UInt64,NamedTuple};
    optimizer, parent_of::AbstractDict = Dict(), kwargs...)
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
    fit_error = Dict{UInt64,String}()
    for r in rep_fits
        r.error === nothing ? (memo[r.eq_hash] = r.fit) : (fit_error[r.eq_hash] = r.error)
    end
    entries  = BatchEntry[]
    failures = FitFailure[]
    emitted_eq_hashes = Set{UInt64}()
    for c in compiled
        (c === nothing || c === :complexity_skip) && continue
        c isa FitFailure && (push!(failures, c); continue)
        if haskey(fit_error, c.eq_hash)
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
    (entries, failures)
end

function _process_batch(
    mechs, prob::IdentifyRateEquationProblem;
    optimizer, max_param_count, eq_complexity_filter::Int = typemax(Int),
    memo::Dict{UInt64,NamedTuple}=Dict{UInt64,NamedTuple}(),
    fitted::Set{UInt64}=Set{UInt64}(),
    parent_of::AbstractDict = Dict(), kwargs...
)
    compiled, reps, rep_idx, n_fitted_skip = _compile_batch(mechs, prob;
        max_param_count, eq_complexity_filter, memo, fitted)
    entries, failures = _fit_batch(compiled, reps, rep_idx, prob, memo;
        optimizer, parent_of, kwargs...)
    (entries, failures,
     count(x -> x === nothing, compiled),
     count(x -> x === :complexity_skip, compiled),
     n_fitted_skip)
end
```

- [ ] **Step 4: Run the scratch driver (GREEN) + the `_process_batch` scratch check**

Run the summary scratch driver → GREEN. Then add the existing `_process_batch`
testset body (`test/test_identify_rate_equation.jl:1123-1164`) to the scratch
driver and run it → it must stay green (the wrapper preserves the contract).
Record both.

- [ ] **Step 5: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Split _process_batch into compile/fit passes; add pre/post summaries"
```

---

### Task 2: Print pre-fit before fitting, post-fit after

**Files:**
- Modify: `src/identify_rate_equation.jl` — the base-tier block (`:782-796` and
  its `_ingest!`/summary print `:784-796`) and the iteration block (`:842-891`).

**Interfaces:**
- Consumes: `_compile_batch`, `_fit_batch`, `_prefit_summary`, `_postfit_summary`
  from Task 1.

- [ ] **Step 1: Rewire the base tier**

Replace the base-tier `_process_batch` call and its post-fit summary print.
After `base = … unique!(…)` (`identify_rate_equation.jl:787-790`), the sequence
becomes:

```julia
    compiled, reps, rep_idx, n_base_fitted_skip = _compile_batch(base, prob;
        max_param_count, eq_complexity_filter, memo, fitted)
    n_base_param_skip = count(x -> x === nothing, compiled)
    n_base_cx_skip    = count(x -> x === :complexity_skip, compiled)
    n_base_nt = count(c -> c isa NamedTuple, compiled)
    _progress(save_dir, show_progress, string(
        "Fitting $(length(base)) initial mechanisms…\n  ",
        _prefit_summary(length(reps), n_base_nt - length(reps),
            n_base_param_skip, n_base_cx_skip, n_base_fitted_skip;
            max_param_count, eq_complexity_filter)))
    base_entries, base_failures = _fit_batch(compiled, reps, rep_idx, prob, memo;
        optimizer, kwargs...)
```

Leave the `isempty(base_entries)` error handling and `_save_initial_csv` calls
(`:796-810`) unchanged. Then the post-fit print (`:791-796` in the original,
now after ingest) becomes:

```julia
    _progress(save_dir, show_progress, string(
        "Base tier: ", _postfit_summary(base_entries, base_failures),
        "\n  ", _best_loss_line(best_loss_by_count, improved)))
```

- [ ] **Step 2: Rewire the iteration block**

Replace the iteration block (`identify_rate_equation.jl:842-891`) — from the
`child_entries, … = _process_batch(children, …)` call through the closing of the
`elseif !isempty(children)` branch — with:

```julia
            compiled, reps, rep_idx, n_child_fitted_skip = _compile_batch(
                children, prob; max_param_count, eq_complexity_filter, memo, fitted)
            n_child_param_skip = count(x -> x === nothing, compiled)
            n_child_cx_skip    = count(x -> x === :complexity_skip, compiled)
            n_child_nt   = count(c -> c isa NamedTuple, compiled)
            n_child_fail = count(c -> c isa FitFailure, compiled)
            if n_child_nt > 0 || n_child_fail > 0 || !isempty(expand_failures)
                iteration += 1
                np_range = n_child_nt == 0 ? "n/a" :
                    let ns = [c.n_params for c in compiled if c isa NamedTuple],
                        lo = minimum(ns), hi = maximum(ns)
                        lo == hi ? string(lo) : "$lo-$hi"
                    end
                _progress(save_dir, show_progress, string(
                    "Iteration $iteration (child n_params $np_range): ",
                    length(to_expand), " parents → ", length(children),
                    " children\n  ",
                    _prefit_summary(length(reps),
                        n_child_nt - length(reps), n_child_param_skip,
                        n_child_cx_skip, n_child_fitted_skip;
                        max_param_count, eq_complexity_filter)))
                child_entries, child_failures = _fit_batch(compiled, reps, rep_idx,
                    prob, memo; optimizer, parent_of, kwargs...)
                append!(child_failures, expand_failures)
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
                _progress(save_dir, show_progress, string("  ",
                    _postfit_summary(child_entries, child_failures),
                    "\n  ", _best_loss_line(best_loss_by_count, improved)))
            elseif !isempty(children)
                _progress(save_dir, show_progress, string(
                    "Expanded ", length(to_expand), " parents → ",
                    length(children), " children | all skipped (",
                    n_child_fitted_skip, " already fit, ",
                    n_child_param_skip, " >", max_param_count, " params, ",
                    n_child_cx_skip, " >", eq_complexity_filter, " complexity)"))
            end
```

Note: `n_params` for a compiled `NamedTuple` equals the fitted mechanism's
`n_params` (fitting does not change it), so the pre-fit `np_range` matches the
old post-fit range.

- [ ] **Step 3: Verify ordering with a small search + record**

Write a scratch driver that runs `identify_rate_equation` on a tiny reaction
(uni-uni, low `max_param_count`, `min_beam_width` small) with `save_dir` a temp
dir and `show_progress=true`, and prints the resulting `progress.log`. Confirm
the order per batch: `Iteration…`/`Fitting…` line, then the `new fits + …` line,
BOTH before any fit-time output, then the `errored | Success…` line, then
`best loss…`. Record the captured log in the report.

- [ ] **Step 4: Commit**

```bash
git add src/identify_rate_equation.jl
git commit -m "Print the fit plan before fitting, results after, per batch"
```

---

## Notes for the implementer

- The base tier passes `optimizer, kwargs...` to `_fit_batch` (no `parent_of`);
  the iteration passes `optimizer, parent_of, kwargs...`. Match the original
  `_process_batch` calls' argument sets.
- The base-tier `pre_best`/`_ingest!`/`improved` lines between `_save_initial_csv`
  and the post-fit print are unchanged — the post-fit print consumes `improved`.
