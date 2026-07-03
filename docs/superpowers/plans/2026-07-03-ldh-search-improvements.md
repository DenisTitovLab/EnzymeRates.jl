# LDH Search Improvements Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the beam-search floor so the search terminates on convergence, make the progress log honest and per-`n_params`, parallelize LOOCV across all workers, correct the docs, and diagnose the RAM growth.

**Architecture:** All code changes live in `src/identify_rate_equation.jl`. The beam floor becomes a cumulative per-count budget (`_select_count!`), the progress line gains a reconciling four-bucket count and a per-`n_params` best-loss map (`_batch_summary` rewrite + `_best_loss_line`), and LOOCV flattens to a `(candidate, fold)` grid (`_cv_fold_loss` extracted from `_loocv`). Docs and a memory-profiling script follow.

**Tech Stack:** Julia, Distributed (`pmap`), Optimization.jl / OptimizationCMAEvolutionStrategy, DataFrames, CSV, Test.

## Global Constraints

- Line length ≤ 92 chars; 4-space indent; match surrounding style.
- `rate_equation` stays allocation-free / sub-100 ns — no task here touches it; the full suite (including `test_rate_equation_performance`) must stay green.
- `min_beam_width` default stays `50`. This plan changes what the floor *means*, not the default.
- TDD: write the failing test, watch it fail, implement, watch it pass, commit.
- Parameter-naming chokepoint guard: introduce no bare `Symbol("K…")`/`k…`/`V…`/`L…` literals.
- **Full suite before every commit:**
  ```bash
  julia --project=. -e 'using Pkg; Pkg.test()'
  ```
- **Focused helper check** (fast, no fitting) — run a single assertion against a pure helper:
  ```bash
  julia --project=. -e 'using EnzymeRates, Test, DataFrames; <assertion>'
  ```
- **Focused file run** (the whole identify test file, standalone):
  ```bash
  julia --project=. -e 'using EnzymeRates, Test, DataFrames, CSV, Random, Statistics, OptimizationCMAEvolutionStrategy, Optimization; using Optimization.SciMLBase: build_solution, ReturnCode, DefaultOptimizationCache; include("test/test_identify_rate_equation.jl")'
  ```

---

## Task 1: Cumulative per-count beam floor

Replace the per-sweep rank floor with a cumulative per-count budget. `_select_beam` is unchanged; a new `_select_count!` wraps it, tracking how many mechanisms each count has already expanded and shrinking the floor accordingly.

**Files:**
- Modify: `src/identify_rate_equation.jl` — add `_select_count!` after `_select_beam` (near line 360); use it in `_beam_search` (per-count loop, ~lines 631-640; add `expanded_by_count` near line 593).
- Test: `test/test_identify_rate_equation.jl` — new `@testset` after the `_select_beam parsimony_cutoff` testset (~line 972).

**Interfaces:**
- Produces: `_select_count!(expanded::Dict{Int,Int}, c::Int, losses::AbstractVector{<:Real}; loss_rel_threshold, loss_abs_threshold, min_beam_width, best_override=nothing, parsimony_cutoff=nothing) -> Vector{Int}` — returns the selected indices (input order, same as `_select_beam`) and advances `expanded[c]` by the number selected.
- Consumes: `_select_beam` (unchanged), `_parsimony_cutoff` (unchanged).

- [ ] **Step 1: Write the failing test**

Add to `test/test_identify_rate_equation.jl` after line 972:

```julia
@testset "_select_count! cumulative per-count floor" begin
    expanded = Dict{Int,Int}()
    # Sweep 1 at count 5: rel cutoff admits only the best (loss 1.0); the
    # floor budget (3) tops it up to the top 3 by loss. expanded[5] -> 3.
    sel1 = EnzymeRates._select_count!(expanded, 5, [1.0, 2.0, 3.0, 4.0, 5.0];
        loss_rel_threshold=1.0, loss_abs_threshold=0.0,
        min_beam_width=3, best_override=1.0)
    @test sort(sel1) == [1, 2, 3]
    @test expanded[5] == 3

    # Sweep 2 at count 5: budget spent (3 of 3). New mechanisms all above the
    # cutoff -> the floor admits NONE (unlike the old per-sweep floor, which
    # would grant a fresh 3). expanded[5] stays 3.
    sel2 = EnzymeRates._select_count!(expanded, 5, [10.0, 11.0, 12.0];
        loss_rel_threshold=1.0, loss_abs_threshold=0.0,
        min_beam_width=3, best_override=1.0)
    @test isempty(sel2)
    @test expanded[5] == 3

    # A cutoff-passer is still admitted after the floor is spent.
    sel3 = EnzymeRates._select_count!(expanded, 5, [1.0, 20.0];
        loss_rel_threshold=1.0, loss_abs_threshold=0.0,
        min_beam_width=3, best_override=1.0)
    @test sel3 == [1]
    @test expanded[5] == 4
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```bash
julia --project=. -e 'using EnzymeRates, Test; EnzymeRates._select_count!(Dict{Int,Int}(), 5, [1.0]; loss_rel_threshold=1.0, loss_abs_threshold=0.0, min_beam_width=1)'
```
Expected: FAIL — `UndefVarError: _select_count! not defined`.

- [ ] **Step 3: Implement `_select_count!`**

In `src/identify_rate_equation.jl`, immediately after the `_select_beam` function (after line 360), add:

```julia
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
```

- [ ] **Step 4: Wire it into `_beam_search`**

In `_beam_search`, add the budget dict beside `best_loss_by_count` (after line 593's `best_loss_by_count = Dict{Int, Float64}()`):

```julia
    # Mechanisms expanded so far per parameter count — the cumulative floor
    # budget. Spent once over the whole search, never re-granted per sweep.
    expanded_by_count = Dict{Int, Int}()
```

Replace the per-count selection loop body (lines 632-639) — change the `_select_beam` call to `_select_count!`:

```julia
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
```

- [ ] **Step 5: Run the focused test to verify it passes**

Run the Focused file run command (Global Constraints). Expected: PASS — the `_select_count! cumulative per-count floor` testset passes and every existing testset still passes.

- [ ] **Step 6: Run the full suite**

Run:
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
Expected: PASS — all tests green (the existing full-pipeline test still runs; the cumulative floor only tightens selection).

- [ ] **Step 7: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Cumulative per-count beam floor: min_beam_width spent once, not per sweep"
```

---

## Task 2: Progress-log reporting

Rewrite `_batch_summary` to report the four reconciling buckets, add `_best_loss_line` for the per-`n_params` best-loss map with improvement marks, and rewire the base-tier and iteration progress lines (ingesting before reporting so the map reflects the current iteration).

**Files:**
- Modify: `src/identify_rate_equation.jl` — rewrite `_batch_summary` (lines 422-439); add `_best_loss_line` after it; rewire base-tier progress (lines 613-619) and iteration progress (lines 642-672) in `_beam_search`.
- Test: `test/test_identify_rate_equation.jl` — update the `_batch_summary` assertions in the `_progress` testset (lines 1014-1028); add a `_best_loss_line` testset; add a progress-log-content check to the `CSV output (new schema)` testset (~line 342).

**Interfaces:**
- Produces:
  - `_batch_summary(entries::Vector{BatchEntry}, failures::Vector{FitFailure}; n_skipped::Int, max_param_count::Int) -> String` — `"<new> new fits + <inh> inherited + <skip> skipped (><cap> params) + <err> errored | Success <p>% | non-Success retcode <q>%"`.
  - `_best_loss_line(best_loss_by_count::Dict{Int,Float64}, improved::Set{Int}) -> String` — `"best loss by n_params: <c>:<loss>[*] … (<annotation>)"`.
- Consumes: `BatchEntry.row.fit_inherited`, `best_loss_by_count` (from Task 1's `_beam_search`), `max_param_count`.

- [ ] **Step 1: Write the failing `_batch_summary` test (update existing)**

Replace the `_batch_summary` assertions at the end of the `_progress` testset (lines 1023-1027) with:

```julia
    s = EnzymeRates._batch_summary([e_succ, e_mt], [f]; n_skipped=4, max_param_count=8)
    @test occursin("2 new fits + 0 inherited + 4 skipped (>8 params) + 1 errored", s)
    @test occursin("Success 50.0%", s)                 # 1 of 2 fitted
    @test occursin("non-Success retcode 50.0%", s)     # e_mt is :MaxTime
    @test !occursin("best loss", s)                    # best loss moved to its own line
```

- [ ] **Step 2: Run to verify it fails**

Run the Focused file run command. Expected: FAIL — `_batch_summary` has no `n_skipped`/`max_param_count` keywords (MethodError) and the old string format does not match.

- [ ] **Step 3: Rewrite `_batch_summary`**

Replace the whole `_batch_summary` function (lines 422-439) with:

```julia
"""
One-line batch summary reconciling every child mechanism into four buckets that
sum to the child count: `new fits` (a genuine optimizer run, `fit_inherited ==
false`), `inherited` (a memoized fit reused, `fit_inherited == true`), `skipped`
(dropped before fitting for exceeding `max_param_count`), and `errored`
(compile/render/fit threw). Trailing `Success` / `non-Success retcode`
percentages are over the fitted set (`new + inherited`).
"""
function _batch_summary(
    entries::Vector{BatchEntry}, failures::Vector{FitFailure};
    n_skipped::Int, max_param_count::Int)
    n_fit   = length(entries)
    n_err   = length(failures)
    n_new   = count(e -> !e.row.fit_inherited, entries)
    n_inh   = n_fit - n_new
    n_succ  = count(e -> e.retcode === :Success, entries)
    n_other = n_fit - n_succ
    pct(x)  = n_fit == 0 ? 0.0 : round(100 * x / n_fit; digits=1)
    string(n_new, " new fits + ", n_inh, " inherited + ",
           n_skipped, " skipped (>", max_param_count, " params) + ",
           n_err, " errored | Success ", pct(n_succ),
           "% | non-Success retcode ", pct(n_other), "%")
end
```

- [ ] **Step 4: Write the failing `_best_loss_line` test**

Add after the `_progress` testset (after line 1028):

```julia
@testset "_best_loss_line" begin
    line = EnzymeRates._best_loss_line(
        Dict(5 => 0.01751, 6 => 0.009316), Set([6]))
    @test occursin("best loss by n_params:", line)
    @test occursin("5:0.01751 ", line)          # unimproved: no star
    @test occursin("6:0.009316*", line)         # improved: starred
    @test occursin("(* improved)", line)

    quiet = EnzymeRates._best_loss_line(Dict(5 => 0.01751), Set{Int}())
    @test occursin("(no improvement)", quiet)
    @test !occursin("*", quiet)
end
```

- [ ] **Step 5: Run to verify it fails**

Run the Focused file run command. Expected: FAIL — `UndefVarError: _best_loss_line`.

- [ ] **Step 6: Implement `_best_loss_line`**

Add immediately after `_batch_summary`:

```julia
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
```

- [ ] **Step 7: Rewire the base-tier progress line**

Replace lines 613-619 (the base-tier save/ingest/progress) with:

```julia
    _save_initial_csv(save_dir,
        vcat([e.row for e in base_entries],
             [_failure_row(f) for f in base_failures]))
    pre_best = copy(best_loss_by_count)
    _ingest!(frontier, cv_pool, best_loss_by_count,
             base_entries; n_cv_candidates)
    improved = Set(c for c in keys(best_loss_by_count)
                   if !haskey(pre_best, c) || best_loss_by_count[c] < pre_best[c])
    n_skipped = length(base) - length(base_entries) - length(base_failures)
    _progress(save_dir, show_progress, string(
        "Base tier: ",
        _batch_summary(base_entries, base_failures; n_skipped, max_param_count),
        "\n  ", _best_loss_line(best_loss_by_count, improved)))
```

- [ ] **Step 8: Rewire the iteration progress line (ingest before report)**

Replace the iteration block (lines 651-672, from `if !isempty(child_entries) || !isempty(child_failures)` through the standalone `!isempty(child_entries) && _ingest!(...)`) with:

```julia
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
                n_skipped = length(children) - length(child_entries) -
                    length(child_failures)
                _progress(save_dir, show_progress, string(
                    "Iteration $iteration (child n_params $np_range): ",
                    length(parents), " parents → ", length(children),
                    " children\n  ",
                    _batch_summary(child_entries, child_failures;
                                   n_skipped, max_param_count),
                    "\n  ", _best_loss_line(best_loss_by_count, improved)))
            end
```

(The old standalone `!isempty(child_entries) && _ingest!(...)` line that followed the `if` block is now inside it — delete the trailing duplicate.)

- [ ] **Step 9: Add a progress-log content check**

In the `CSV output (new schema)` testset, after line 346 (`@test filesize(joinpath(save_dir, "progress.log")) > 0`), add:

```julia
        log_text = read(joinpath(save_dir, "progress.log"), String)
        @test occursin("new fits", log_text)
        @test occursin("skipped (>", log_text)
        @test occursin("best loss by n_params:", log_text)
```

- [ ] **Step 10: Run the focused file, then the full suite**

Run the Focused file run command. Expected: PASS. Then:
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
Expected: PASS.

- [ ] **Step 11: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Progress log: reconcile new/inherited/skipped/errored + per-n_params best-loss map"
```

---

## Task 3: Cross-validation parallelism

Extract the per-fold body of `_loocv` into `_cv_fold_loss`, keep `_loocv` as a serial loop over it, and flatten `_cv_model_selection`'s `pmap` to the full `(candidate, fold)` grid.

**Files:**
- Modify: `src/identify_rate_equation.jl` — add `_cv_fold_loss` before `_loocv` (~line 718); rewrite `_loocv` body (lines 730-763) as a loop over it; flatten the `pmap` in `_cv_model_selection` (lines 966-971).
- Test: `test/test_identify_rate_equation.jl` — add a `_cv_fold_loss` testset and a flatten-equivalence testset after the `_loocv is loud on fit failure` testset (~line 454, still inside the outer `identify_rate_equation` testset) or at file top-level after line 456.

**Interfaces:**
- Produces: `_cv_fold_loss(mechanism::AbstractEnzymeMechanism, prob::IdentifyRateEquationProblem, held_out; optimizer, kwargs...) -> Float64` — one fold: fit on the complement of `held_out`, score on `held_out`, error on non-finite (naming the group), floor at `eps(Float64)`.
- Consumes: `_subset_data`, `_evaluate_loss` (both unchanged), `compile_mechanism`.

- [ ] **Step 1: Write the failing `_cv_fold_loss` test**

Add after the `_loocv is loud on fit failure` testset (after line 454, before the closing `end` at line 456):

```julia
    @testset "_cv_fold_loss: one fold, floored + finite; _loocv loops it" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
        end
        m = EnzymeRates.EnzymeMechanism(
            first(EnzymeRates.init_mechanisms(rxn)))
        data = DataFrame(
            S = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            P = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
            Rate = [0.5, 0.8, 1.0, 1.1, 1.2, 1.3],
            group = [1, 1, 2, 2, 3, 3])
        prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)
        kw = (; optimizer=CMAEvolutionStrategyOpt(),
                n_restarts=2, maxtime=2.0, maxiters=500)

        one = EnzymeRates._cv_fold_loss(m, prob, 2; kw...)
        @test one isa Float64
        @test one >= eps(Float64) && isfinite(one)
    end
```

(Do NOT assert bitwise equality between two independent real-optimizer fits —
CMA-ES draws restarts from the unseeded global RNG under a wall-clock `maxtime`,
so `_loocv(m,prob) == [_cv_fold_loss(m,prob,g) for g]` is flaky. `_loocv`'s
faithfulness as a serial loop over `_cv_fold_loss` is proven deterministically
by the Step 5 flatten-equivalence test with the `_CountingStubOpt` stub, plus
the existing `_loocv returns per-fold scores` testset.)

- [ ] **Step 2: Run to verify it fails**

Run the Focused file run command. Expected: FAIL — `UndefVarError: _cv_fold_loss`.

- [ ] **Step 3: Extract `_cv_fold_loss` and refactor `_loocv`**

In `src/identify_rate_equation.jl`, add before `_loocv` (before line 730's docstring):

```julia
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
    isfinite(test_loss) || error(
        "LOOCV produced a non-finite test loss for held-out group " *
        "$held_out — the fit is unusable; aborting model selection.")
    max(test_loss, eps(Float64))
end
```

Replace the `_loocv` body (the `for held_out in groups … end; scores` block, lines 735-762) with a loop over `_cv_fold_loss`:

```julia
    groups = unique(prob.data.group)
    [_cv_fold_loss(mechanism, prob, g; optimizer, kwargs...) for g in groups]
```

(Keep the `_loocv` docstring at lines 718-729.)

- [ ] **Step 4: Run to verify the `_cv_fold_loss` test passes**

Run the Focused file run command. Expected: PASS — the new testset passes and the existing `_loocv returns per-fold scores` / `_loocv is loud on fit failure` testsets still pass.

- [ ] **Step 5: Write the failing flatten-equivalence test**

Add at file top-level **after the `_DEDUP_SIG1` / `_DEDUP_SIG2` / `_CountingStubOpt` definitions** (near the file end, after the `LOOCV eq_hash-uniqueness guard (§4)` testset) — this testset consumes those fixtures, and Julia evaluates the file top-to-bottom, so placing it earlier UndefVars:

```julia
@testset "_cv_model_selection flatten reproduces serial LOOCV" begin
    # Deterministic stub optimizer → identical fits whether folds run serially
    # or across the flattened (candidate, fold) grid.
    recon(sig) = EnzymeRates.Mechanism(Core.eval(EnzymeRates, Meta.parse(sig))())
    m1 = recon(_DEDUP_SIG1)
    em1 = EnzymeRates.compile_mechanism(m1)
    fkeys = EnzymeRates.fitted_params(em1)
    h = string(EnzymeRates._rate_eq_dedup_key(rate_equation_string(em1)),
               base=16, pad=16)
    data = (group = ["G1", "G1", "G2", "G2", "G3", "G3"],
            Rate = [0.5, 0.8, 1.0, 1.1, 0.9, 1.2],
            A = [1.0, 2.0, 1.0, 2.0, 1.5, 2.5], B = [0.5, 0.5, 1.0, 1.0, 0.7, 0.7],
            P = [0.1, 0.2, 0.1, 0.2, 0.15, 0.25], Q = [0.3, 0.3, 0.4, 0.4, 0.35, 0.35])
    prob = IdentifyRateEquationProblem(EnzymeRates.reaction(m1), data; Keq=2.0)
    mechs = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m1]
    mkrow(loss) = (n_params=length(fkeys), loss=loss, mechanism_type="M",
        rate_equation="v", retcode="Success", error=missing,
        fitted_param_names=Tuple(fkeys),
        fitted_param_values=Tuple(fill(1.0, length(fkeys))),
        eq_hash=h, fit_inherited=false)
    df = EnzymeRates._rows_to_dataframe([mkrow(0.5)])

    res = EnzymeRates._cv_model_selection(mechs, df, prob;
        n_cv_candidates=5, optimizer=_CountingStubOpt(; uval=log(5.0)),
        se_threshold=1.0, perm_p_threshold=1.0, save_dir=mktempdir(),
        show_progress=false, n_restarts=1, maxtime=1.0)

    groups = unique(prob.data.group)
    flat = [res.cv_results[1, Symbol("cv_fold_$g")] for g in groups]
    serial = EnzymeRates._loocv(EnzymeRates.compile_mechanism(m1), prob;
        optimizer=_CountingStubOpt(; uval=log(5.0)), n_restarts=1, maxtime=1.0)
    @test flat == serial
end
```

- [ ] **Step 6: Run to verify it fails**

Run the Focused file run command. Expected: FAIL — with the un-flattened `pmap`, the numbers still match, so this test may PASS before the change; if it passes, note it exercises the same primitive. Confirm it fails only if the reduction is wrong. (If it passes here, keep it as a regression lock and proceed — the flatten in Step 7 must keep it green.)

- [ ] **Step 7: Flatten the `pmap` in `_cv_model_selection`**

Replace the `pmap` block (lines 966-971):

```julia
    fold_scores_per_candidate = pmap(
        candidate_mechs
    ) do mech
        m = compile_mechanism(mech)
        _loocv(m, prob; optimizer, kwargs...)
    end
```

with the flattened grid:

```julia
    # Flatten LOOCV to a (candidate, fold) grid so all folds of all candidates
    # run across every worker, not one candidate per worker with serial folds.
    groups = unique(prob.data.group)
    tasks = [(ci, g) for ci in eachindex(candidate_mechs) for g in groups]
    flat = pmap(tasks) do task
        ci, g = task
        m = compile_mechanism(candidate_mechs[ci])
        (ci, g, _cv_fold_loss(m, prob, g; optimizer, kwargs...))
    end
    gi = Dict(g => i for (i, g) in enumerate(groups))
    fold_scores_per_candidate = [Vector{Float64}(undef, length(groups))
                                 for _ in candidate_mechs]
    for (ci, g, s) in flat
        fold_scores_per_candidate[ci][gi[g]] = s
    end
```

- [ ] **Step 8: Run the focused file, then the full suite**

Run the Focused file run command. Expected: PASS — the flatten-equivalence test and every existing CV test pass. Then:
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
Expected: PASS.

- [ ] **Step 9: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Flatten LOOCV to a (candidate, fold) grid so all workers stay busy"
```

---

## Task 4: Documentation

Correct the filter description everywhere it appears: `best(n-1)` → `best(<n)` (best over *all* smaller counts) and the width floor from per-rank/hard to a cumulative per-count budget.

**Files:**
- Modify: `docs/src/identify/model_selection.md` (lines 23-46).
- Modify: `src/identify_rate_equation.jl` — the `identify_rate_equation` docstring (`loss_parsimony_threshold` bullet ~lines 114-119; "Beam selection" section ~lines 150-164) and the `_select_beam` docstring (~lines 317-334).
- No test; the CI doctest/checkdocs run in the full suite.

- [ ] **Step 1: Fix `model_selection.md` — the cutoff formula**

Replace lines 29-34:

```
  ```
  min(loss_rel_threshold * best(n) + loss_abs_threshold,
      loss_parsimony_threshold * best(n-1))
  ```

  where `best(k)` is the lowest loss seen at parameter count `k`. The
  `best(n-1)` term is dropped at the base count (there is no `n-1` level).
```

with:

```
  ```
  min(loss_rel_threshold * best(n) + loss_abs_threshold,
      loss_parsimony_threshold * best(<n))
  ```

  where `best(n)` is the lowest loss at parameter count `n` and `best(<n)` is
  the lowest loss over *all* smaller counts — an added parameter must beat the
  best simpler model of any size. The `best(<n)` term is dropped at the base
  count (no smaller level exists yet).
```

- [ ] **Step 2: Fix `model_selection.md` — the width floor**

Replace lines 35-37:

```
- **Width floor** — its loss-rank (ascending) is `≤ min_beam_width`, which
  always keeps the top `min_beam_width` mechanisms even when the loss cutoff
  would admit fewer.
```

with:

```
- **Width floor** — a cumulative per-count budget. The floor keeps expanding
  mechanisms at a parameter count until `min_beam_width` of them have been
  expanded *over the whole search*, then stops; after that only the loss cutoff
  admits at that count. The budget is spent once, not re-granted each time the
  count is revisited.
```

- [ ] **Step 3: Fix `model_selection.md` — the knob table**

Replace the `loss_parsimony_threshold` and `min_beam_width` rows (lines 45-46):

```
| `loss_parsimony_threshold` | `1.01` | An added parameter must earn its keep: keep expanding only if the loss is within this factor of `best(n-1)`, the best model with one fewer parameter. Set to `Inf` to disable it. |
| `min_beam_width` | `50` | Hard floor on the number kept per parameter count. The loss thresholds can only tighten the beam below this floor, never widen it past it. |
```

with:

```
| `loss_parsimony_threshold` | `1.01` | An added parameter must earn its keep: keep expanding only if the loss is within this factor of `best(<n)`, the best model of any smaller parameter count. Set to `Inf` to disable it. |
| `min_beam_width` | `50` | Cumulative floor: expand at least this many mechanisms per parameter count over the whole search. Spent once, then only the loss cutoff admits at that count. |
```

- [ ] **Step 4: Fix the `identify_rate_equation` docstring**

In `src/identify_rate_equation.jl`, update the `loss_parsimony_threshold` bullet (lines 114-119) to say `best(<n)` / "the best model of any smaller parameter count", and the "Beam selection" prose (lines 150-164) so the second cutoff term reads `loss_parsimony_threshold * best(<n)` with `best(<n)` defined as the lowest loss over all counts below `n`, and the floor bullet describes the cumulative budget (expand at least `min_beam_width` per count over the whole search, spent once). Update the `_select_beam` docstring (lines 317-334) so `min_beam_width` is described as the *remaining* per-count budget passed by the caller, not a per-sweep rank floor.

Concretely, change every `best(n-1)` / "one fewer parameter" to `best(<n)` / "any smaller parameter count", and every "hard floor" / "top `min_beam_width` always qualify per sweep" to the cumulative-budget wording.

- [ ] **Step 5: Run the full suite (doctest/checkdocs included)**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
Expected: PASS — no doctest or `checkdocs` failure (the edited docstrings contain no executable examples).

- [ ] **Step 6: Commit**

```bash
git add docs/src/identify/model_selection.md src/identify_rate_equation.jl
git commit -m "Docs: filter is best(<n) parsimony + cumulative min_beam_width floor"
```

---

## Task 5: Memory diagnostic (profile, don't fix)

Run a reduced LDH search locally with a non-invasive memory sampler to isolate the dominant growth term (worker code-cache vs. master data retention), and record the finding. **No fix in this task** — the fix is a follow-up plan once the term is known.

**Files:**
- Create: `docs/ldh_hpc_results/profile_memory.jl` — the sampler + reduced run.
- Create: `docs/superpowers/specs/2026-07-03-ldh-memory-findings.md` — the recorded conclusion.

- [ ] **Step 1: Write the profiling script**

Create `docs/ldh_hpc_results/profile_memory.jl`:

```julia
# ABOUTME: Local memory profiler for identify_rate_equation on reduced LDH data.
# ABOUTME: Samples master + per-worker gc_live_bytes / maxrss on a timer to
# ABOUTME: isolate code-cache growth (GC-immune) from data retention.
using Pkg
Pkg.activate(joinpath(@__DIR__))
using Distributed, Dates, CSV, DataFrames

addprocs(4; exeflags = ["--project"])
@everywhere using EnzymeRates, OptimizationCMAEvolutionStrategy

raw = CSV.read(joinpath(@__DIR__, "Enzyme data", "LDH_data.csv"), DataFrame)
filter!(row -> row.Rate != 0.0, raw)
data = (group = String.(raw.Article .* "_" .* raw.Fig),
        Rate = Float64.(raw.Rate), NADH = Float64.(raw.NADH),
        Pyruvate = Float64.(raw.Pyruvate), Lactate = Float64.(raw.Lactate),
        NAD = Float64.(raw.NAD))
rxn = @enzyme_reaction begin
    substrates:NADH[C21H29N7O14P2], Pyruvate[C3H4O3]
    products:Lactate[C3H6O3], NAD[C21H27N7O14P2]
    oligomeric_state:4
end
prob = IdentifyRateEquationProblem(rxn, data; Keq=20000.0, scale_k_to_kcat=1.0)

# Sampler: master + each worker, every 5 s, to a CSV keyed by wall-seconds.
@everywhere _mem() = (gc_live = Base.gc_live_bytes(), maxrss = Sys.maxrss())
samples = DataFrame(t=Float64[], who=String[], gc_live=Int[], maxrss=Int[])
t0 = time()
sampling = Ref(true)
sampler = @async while sampling[]
    t = time() - t0
    m = _mem(); push!(samples, (t, "master", m.gc_live, m.maxrss))
    for w in workers()
        wm = remotecall_fetch(_mem, w)
        push!(samples, (t, "worker$w", wm.gc_live, wm.maxrss))
    end
    sleep(5)
end

# Reduced run: small enough to finish in a few minutes, large enough to grow.
results = identify_rate_equation(prob;
    optimizer=CMAEvolutionStrategyOpt(),
    max_param_count=9, min_beam_width=10,
    loss_rel_threshold=1.3, loss_abs_threshold=0.001,
    loss_parsimony_threshold=1.01,
    save_dir=joinpath(@__DIR__, "memprofile_results"))

sampling[] = false; wait(sampler)

# Does a forced GC reclaim? If maxrss barely drops, growth is code cache.
@everywhere GC.gc(true); GC.gc(true)
after = _mem(); push!(samples, (time() - t0, "master_postGC", after.gc_live, after.maxrss))
for w in workers()
    wm = remotecall_fetch(() -> (GC.gc(true); _mem()), w)
    push!(samples, (time() - t0, "worker$(w)_postGC", wm.gc_live, wm.maxrss))
end

CSV.write(joinpath(@__DIR__, "memprofile_samples.csv"), samples)
println("peak master maxrss: ", maximum(samples[samples.who .== "master", :maxrss]))
for w in workers()
    ws = samples[samples.who .== "worker$w", :]
    isempty(ws) && continue
    println("worker$w maxrss: first ", first(ws.maxrss), " peak ", maximum(ws.maxrss))
end
for p in workers(); rmprocs(p); end
```

- [ ] **Step 2: Run the profiler**

```bash
cd docs/ldh_hpc_results && julia --project profile_memory.jl 2>&1 | tee memprofile_stdout.txt
```
Expected: completes in a few minutes; writes `memprofile_samples.csv` and prints per-worker peak `maxrss`.

- [ ] **Step 3: Analyze and record the finding**

Inspect `memprofile_samples.csv`. Decide the dominant term:
- If per-worker `maxrss` grows monotonically while `gc_live` stays roughly flat, and the `*_postGC` rows do **not** drop `maxrss` → **code-cache growth** (`@generated rate_equation` per `Sig`); `GC.gc()` will not help.
- If `gc_live` grows and `*_postGC` drops `maxrss` → **data retention** (candidate: the wide per-iteration DataFrame or `memo`).

Write `docs/superpowers/specs/2026-07-03-ldh-memory-findings.md` recording: the peak master/worker `maxrss`, whether GC reclaimed, the identified dominant term, and a one-paragraph recommended fix direction (e.g. worker recycling for code cache; streaming/narrow CSV writing for the DataFrame). Do **not** implement the fix here.

- [ ] **Step 4: Commit the diagnostic**

```bash
git add docs/ldh_hpc_results/profile_memory.jl docs/superpowers/specs/2026-07-03-ldh-memory-findings.md
git commit -m "Local memory profiler for identify_rate_equation + recorded findings"
```

(Leave `memprofile_results/`, `memprofile_samples.csv`, and `memprofile_stdout.txt` untracked — do not `git add` them.)

---

## Self-Review

**Spec coverage:**
- Change 1 (cumulative floor) → Task 1. ✓
- Change 2 (reporting: new/inherited/skipped + Option C map) → Task 2. ✓
- Change 3 (CV flatten) → Task 3. ✓
- Change 5 (docs) → Task 4. ✓
- Change 4 (memory: diagnose first) → Task 5 (diagnostic only; fix is a follow-up, per spec). ✓

**Placeholder scan:** every code step shows the full code; the only deferred item is Change 4's fix, which the spec itself defers. No `TODO`/"handle edge cases"/"similar to" placeholders.

**Type consistency:** `_select_count!` returns `Vector{Int}` and mutates `expanded::Dict{Int,Int}`, matching `expanded_by_count` in `_beam_search`. `_batch_summary` keyword `n_skipped::Int`/`max_param_count::Int` matches both call sites. `_best_loss_line(::Dict{Int,Float64}, ::Set{Int})` matches the `improved` set and `best_loss_by_count` built in `_beam_search`. `_cv_fold_loss` returns `Float64`; `_loocv` returns `Vector{Float64}`; the flattened reducer fills `Vector{Float64}` per candidate — matches `cv_df.cv_fold_scores` expectations downstream.

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-07-03-ldh-search-improvements.md`.
