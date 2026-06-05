# Actual-Param-Count Beam Search Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the parameter-*estimate*-driven beam search with one keyed on the **actual** fitted-param count, bound its memory (mandatory CSVs + top-N-per-count CV pool), and replace the ~200-line canonical-rate-eq-hash machinery with a ~5-line comment-stripped string key.

**Architecture:** `init_mechanisms` are all fit up front (`initial_mechanisms.csv`); an ascending **advancing-target sweep** then expands mechanisms by actual param count — each iteration pops all unexpanded entries at count ≤ a monotonic `target`, beam-selects per count, expands as one batch, fits children, and writes `equation_search_iteration_N.csv`. Each worker does compile→cap-check→fit locally (one `pmap`, no cross-worker recompile). Dedup is structural-only (`_dedup_flat`); `eq_hash` is a CSV tag + the LOOCV distinct-equation key. Termination rests on move irreversibility + the param cap.

**Tech Stack:** Julia 1.9+, package `EnzymeRates` at `~/.julia/dev/EnzymeRates/`. Full suite: `julia --project -e 'using Pkg; Pkg.test()'` from the package root. Single-test-file runs are NOT supported standalone (shared fixtures via `runtests.jl`); for fast inner-loop checks use `julia --project -e '<snippet>'`.

**Design doc:** `docs/superpowers/specs/2026-06-04-actual-param-count-beam-search-design.md`

**Invariants you MUST NOT break (`.claude/CLAUDE.md`):**
- `rate_equation` stays allocation-free and `< 100e-9` s/call (`test/test_rate_eq_derivation.jl`). This plan does **not** touch the derivation. If a change would force it to allocate/slow, STOP and ask Denis.
- Compile-budget gate (`test/test_compile_budget.jl`, ~750 bi-bi init) and the 18-export count hold — no new `Sig` types, no export changes.
- Style: 92-char lines, 4-space indent. Each source file starts with two `# ABOUTME:` lines (don't remove). Remove unused features entirely (no dead functions).

---

## File Structure

- **`src/identify_rate_equation.jl`** — new: `_rate_eq_dedup_key`, `_default_save_dir`, `_save_initial_csv`, `_save_iteration_csv`, `BatchEntry`, `_process_batch`, `_ingest!`, `_offer_cv!`; rewritten: `_beam_search`, `_select_beam` (adds `best_override`), `_rows_to_dataframe` (drops a column), `identify_rate_equation` (save_dir default); deleted: canonical-hash machinery, `_project_cached_params`, `_CachedFitResult`, `_CompiledMechanismResult`, `_save_level_csv`.
- **`src/mechanism_enumeration.jl`** — new: `_dedup_flat`; rewritten: `expand_mechanisms`, `_add_expansions_mech!` (flat vector); deleted: `_push_mech!`, `_n_fit_params_estimate` (both overloads).
- **`src/EnzymeRates.jl`** — add `using Dates`.
- **`src/rate_eq_derivation.jl`** — fix one stale doc comment (`:1171`) referencing a deleted helper.
- **`test/test_identify_rate_equation.jl`** / **`test/test_mechanism_enumeration.jl`** — migrate (new CSV names, flat `expand_mechanisms`, drop estimate/canonical tests, add per-function unit tests).
- **`.claude/CLAUDE.md`** — doc sync.

---

## Task 0: Baseline

**Files:** none (measurement only)

- [ ] **Step 1: Record the baseline**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -40`
Record pass/fail/broken counts and wall time. Confirm `test_rate_eq_derivation.jl` (perf gate) and `test_compile_budget.jl` are green on the untouched tree.

- [ ] **Step 2: Confirm the working branch**

Run: `git -C ~/.julia/dev/EnzymeRates branch --show-current`
Expected: `actual-param-count-beam-search` (the spec lives here). If on `main`, create/switch: `git checkout actual-param-count-beam-search`.

---

## Task 1: `_rate_eq_dedup_key`

**Files:**
- Modify: `src/identify_rate_equation.jl` (add the function near the other `_rate_eq_*` helpers, e.g. just above `_select_beam`)
- Test: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Write the failing test**

Add to `test/test_identify_rate_equation.jl` inside a new `@testset`:

```julia
@testset "_rate_eq_dedup_key" begin
    base = "(; K_a, k_b) = params\n(; A) = concs\n" *
           "# Haldane constraints:\nk_r = (1/Keq)*K_a\nv = k_b*A/K_a"
    # differs only in a comment header + a substituted-into-v provenance line:
    a = "# Wegscheider constraints:\nK_x = K_a  (substituted into v)\n" * base
    b = "# Wegscheider constraints:\nK_y = K_a  (substituted into v)\n" * base
    @test EnzymeRates._rate_eq_dedup_key(a) ==
          EnzymeRates._rate_eq_dedup_key(b)
    # differs in a Haldane definition -> different key:
    c = replace(base, "k_r = (1/Keq)*K_a" => "k_r = (2/Keq)*K_a")
    @test EnzymeRates._rate_eq_dedup_key(base) !=
          EnzymeRates._rate_eq_dedup_key(c)
    # differs in the v= line -> different key:
    d = replace(base, "v = k_b*A/K_a" => "v = k_b*A/(K_a + A)")
    @test EnzymeRates._rate_eq_dedup_key(base) !=
          EnzymeRates._rate_eq_dedup_key(d)
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `julia --project -e 'using EnzymeRates; EnzymeRates._rate_eq_dedup_key("x")'`
Expected: FAIL — `UndefVarError: _rate_eq_dedup_key`.

- [ ] **Step 3: Implement**

Add to `src/identify_rate_equation.jl`:

```julia
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
        !startswith(l, "#") && !occursin("(substituted into v)", l)
    end
    hash(join(kept, '\n'))
end
```

- [ ] **Step 4: Verify it passes**

Run: `julia --project -e 'using EnzymeRates; const E=EnzymeRates;
b="(; K_a) = params\nv = K_a";
println(E._rate_eq_dedup_key("# h\n"*b) == E._rate_eq_dedup_key(b))'`
Expected: `true`.

- [ ] **Step 5: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Add _rate_eq_dedup_key (comment-stripped equation key)"
```

---

## Task 2: `_dedup_flat`

**Files:**
- Modify: `src/mechanism_enumeration.jl:1789-1798` (refactor `dedup!`, add `_dedup_flat`)
- Test: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write the failing test**

Add a `@testset` to `test/test_mechanism_enumeration.jl`:

```julia
@testset "_dedup_flat" begin
    rxn = @enzyme_reaction begin
        substrates:S[C]
        products:P[C]
    end
    ms = collect(EnzymeRates.init_mechanisms(rxn))
    dup = vcat(ms, deepcopy(ms))          # every mechanism twice
    out = EnzymeRates._dedup_flat(dup)
    @test length(out) == length(EnzymeRates._dedup_flat(collect(ms)))
    @test length(out) <= length(dup)
    @test EnzymeRates._dedup_flat(Union{EnzymeRates.Mechanism,
        EnzymeRates.AllostericMechanism}[]) == []
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `julia --project -e 'using EnzymeRates; EnzymeRates._dedup_flat([])'`
Expected: FAIL — `UndefVarError: _dedup_flat`.

- [ ] **Step 3: Implement**

In `src/mechanism_enumeration.jl`, replace the body of `dedup!` (currently `src/mechanism_enumeration.jl:1789-1798`) and add `_dedup_flat` just above it:

```julia
"""
    _dedup_flat(mechs::Vector)

Canonicalize each mechanism in place via `_canonicalize_mechanism!`, then
`unique!` so structurally-equivalent mechanisms collapse. The cheap,
pre-compile dedup used on both the base set and every expansion batch.
"""
function _dedup_flat(mechs::Vector)
    for m in mechs
        _canonicalize_mechanism!(m)
    end
    unique!(mechs)
    mechs
end

"""
    dedup!(cache::Dict{Int, <:Vector})

Apply `_dedup_flat` to each bucket; drop emptied buckets.
"""
function dedup!(cache::Dict{Int, <:Vector})
    for (pc, mechs) in cache
        _dedup_flat(mechs)
        isempty(mechs) && delete!(cache, pc)
    end
    cache
end
```

- [ ] **Step 4: Verify it passes**

Run: `julia --project -e 'using EnzymeRates; const E=EnzymeRates;
rxn=E.@enzyme_reaction begin; substrates:S[C]; products:P[C]; end;
ms=collect(E.init_mechanisms(rxn));
println(length(E._dedup_flat(vcat(ms,deepcopy(ms)))) == length(E._dedup_flat(copy(ms))))'`
Expected: `true`.

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Factor _dedup_flat out of dedup!"
```

---

## Task 3: `_default_save_dir` (+ `using Dates`)

**Files:**
- Modify: `src/EnzymeRates.jl` (add `using Dates`)
- Modify: `src/identify_rate_equation.jl` (add `_default_save_dir`)
- Test: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "_default_save_dir" begin
    mktempdir() do tmp
        cd(tmp) do
            d1 = EnzymeRates._default_save_dir()
            @test endswith(d1, "_results")
            mkpath(d1)
            d2 = EnzymeRates._default_save_dir()
            @test d2 == d1 * "_2"
            mkpath(d2)
            @test EnzymeRates._default_save_dir() == d1 * "_3"
        end
    end
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `julia --project -e 'using EnzymeRates; EnzymeRates._default_save_dir()'`
Expected: FAIL — `UndefVarError: _default_save_dir`.

- [ ] **Step 3: Implement**

Add `using Dates` to `src/EnzymeRates.jl` alongside the other `using` lines (near `using Random`).

Add to `src/identify_rate_equation.jl`:

```julia
"""
Default results directory: the first non-existent `<date>_results[_N]`
directory in the cwd (e.g. `2026_06_04_results`, then `…_results_2`, `_3`).
"""
function _default_save_dir()
    base = string(Dates.format(Dates.today(), "yyyy_mm_dd"), "_results")
    isdir(base) || return base
    n = 2
    while isdir(string(base, "_", n)); n += 1; end
    string(base, "_", n)
end
```

- [ ] **Step 4: Verify it passes**

Run: `julia --project -e 'using EnzymeRates; cd(mktempdir()) do;
println(endswith(EnzymeRates._default_save_dir(), "_results")); end'`
Expected: `true`.

- [ ] **Step 5: Commit**

```bash
git add src/EnzymeRates.jl src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Add _default_save_dir with date + numeric suffix"
```

---

## Task 4: CSV writers `_save_initial_csv` / `_save_iteration_csv`

**Files:**
- Modify: `src/identify_rate_equation.jl` (add both; leave `_save_level_csv` for now — deleted in Task 11)
- Test: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "csv writers" begin
    rows = [(
        n_params = 5, loss = 1.0, mechanism_type = "M",
        rate_equation = "v = 1", fitted_param_names = (:K_a,),
        fitted_param_values = (2.0,), eq_hash = "abc",
    )]
    mktempdir() do tmp
        EnzymeRates._save_initial_csv(tmp, rows)
        @test isfile(joinpath(tmp, "initial_mechanisms.csv"))
        EnzymeRates._save_iteration_csv(tmp, rows, 3)
        @test isfile(joinpath(tmp, "equation_search_iteration_3.csv"))
        df = CSV.read(joinpath(tmp, "equation_search_iteration_3.csv"), DataFrame)
        @test df.n_params == [5]
        @test "eq_hash" in names(df)
    end
end
```
(Ensure `using CSV, DataFrames` is in scope in this test file — it already is for existing CSV tests.)

- [ ] **Step 2: Run it to verify it fails**

Run: `julia --project -e 'using EnzymeRates; EnzymeRates._save_initial_csv(".", [])'`
Expected: FAIL — `UndefVarError: _save_initial_csv`.

- [ ] **Step 3: Implement**

Add to `src/identify_rate_equation.jl` (near `_rows_to_dataframe`):

```julia
"""Save the base-tier fit (all init mechanisms) to `initial_mechanisms.csv`."""
function _save_initial_csv(save_dir::String, rows)
    isdir(save_dir) || mkpath(save_dir)
    CSV.write(joinpath(save_dir, "initial_mechanisms.csv"),
              _rows_to_dataframe(rows))
end

"""
Save one expansion iteration to `equation_search_iteration_<iteration>.csv`.
`iteration` is a 1-based sequential counter, NOT a parameter count — the
real fitted count is the `n_params` column of each row.
"""
function _save_iteration_csv(save_dir::String, rows, iteration::Int)
    isdir(save_dir) || mkpath(save_dir)
    CSV.write(joinpath(save_dir,
              "equation_search_iteration_$(iteration).csv"),
              _rows_to_dataframe(rows))
end
```

Note: these call `_rows_to_dataframe`, which today still emits the `fit_inherited_from_estimate` column. That column is dropped in Task 9; this test does not assert its absence, so it passes either way.

- [ ] **Step 4: Verify it passes**

Run: `julia --project -e 'using EnzymeRates, CSV, DataFrames;
rows=[(n_params=5,loss=1.0,mechanism_type="M",rate_equation="v=1",
fitted_param_names=(:K_a,),fitted_param_values=(2.0,),eq_hash="abc")];
t=mktempdir(); EnzymeRates._save_iteration_csv(t,rows,3);
println(isfile(joinpath(t,"equation_search_iteration_3.csv")))'`
Expected: `true`.

- [ ] **Step 5: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Add initial/iteration CSV writers"
```

---

## Task 5: `_select_beam` gains `best_override`

**Files:**
- Modify: `src/identify_rate_equation.jl:277-299`
- Test: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "_select_beam best_override" begin
    losses = [1.0, 1.5, 3.0]
    kw = (loss_rel_threshold=1.2, loss_abs_threshold=0.0, min_beam_width=1)
    # without override: best = min = 1.0, cutoff = 1.2 -> only index 1
    @test EnzymeRates._select_beam(losses; kw...) == [1]
    # override best = 2.0 -> cutoff 2.4 -> indices 1 and 2
    @test EnzymeRates._select_beam(losses; kw..., best_override=2.0) == [1, 2]
    # min_beam_width still honored
    @test EnzymeRates._select_beam(losses;
        loss_rel_threshold=1.0, loss_abs_threshold=0.0,
        min_beam_width=2, best_override=0.0) == [1, 2]
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `julia --project -e 'using EnzymeRates;
EnzymeRates._select_beam([1.0]; loss_rel_threshold=1.0,
loss_abs_threshold=0.0, min_beam_width=1, best_override=2.0)'`
Expected: FAIL — `MethodError`/unknown kwarg `best_override`.

- [ ] **Step 3: Implement**

Edit `_select_beam` (`src/identify_rate_equation.jl:277-299`) — add the kwarg and use it for `best`:

```julia
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
    sort!(selected)
end
```

- [ ] **Step 4: Verify it passes**

Run: `julia --project -e 'using EnzymeRates;
println(EnzymeRates._select_beam([1.0,1.5,3.0]; loss_rel_threshold=1.2,
loss_abs_threshold=0.0, min_beam_width=1, best_override=2.0) == [1,2])'`
Expected: `true`.

- [ ] **Step 5: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Add best_override kwarg to _select_beam"
```

---

## Task 6: `BatchEntry` + `_process_batch`

**Files:**
- Modify: `src/identify_rate_equation.jl` (add struct + function; place near `_beam_search`)
- Test: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Write the failing test**

Self-contained (build a tiny uni-uni `prob` + `PyCMAOpt()`, mirroring the `_loocv` testset's pattern at `test/test_identify_rate_equation.jl:356-377`). `pmap_function=map` runs single-process. Add a new top-level `@testset`:

```julia
@testset "_process_batch" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    data = DataFrame(
        S = [1.0, 2.0, 3.0, 4.0],
        P = [0.1, 0.2, 0.3, 0.4],
        Rate = [0.5, 0.8, 1.0, 1.1],
        group = [1, 1, 2, 2],
    )
    prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)
    ms = EnzymeRates._dedup_flat(collect(EnzymeRates.init_mechanisms(rxn)))
    entries = EnzymeRates._process_batch(ms, prob;
        pmap_function=map, optimizer=PyCMAOpt(),
        max_param_count=20, n_restarts=1, maxtime=1.0)
    @test !isempty(entries)
    @test all(e -> e isa EnzymeRates.BatchEntry, entries)
    @test all(e -> e.n_params == length(e.row.fitted_param_names), entries)
    @test all(e -> e.row.eq_hash isa String, entries)
    # cap filter: nothing over the cap is fit
    capped = EnzymeRates._process_batch(ms, prob;
        pmap_function=map, optimizer=PyCMAOpt(),
        max_param_count=0, n_restarts=1, maxtime=1.0)
    @test isempty(capped)
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `julia --project -e 'using EnzymeRates; EnzymeRates.BatchEntry'`
Expected: FAIL — `UndefVarError: BatchEntry`.

- [ ] **Step 3: Implement**

Add to `src/identify_rate_equation.jl` (above `_beam_search`):

```julia
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
            @debug("process_batch failed",
                   exception=(e, catch_backtrace()))
            nothing
        end
    end
    BatchEntry[r for r in results if r !== nothing]
end
```

- [ ] **Step 4: Run the targeted test**

Run the full suite (single-file not supported): `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -20`
Expected: the `_process_batch` testset passes (other not-yet-migrated tests may still be green at this point — nothing has been removed yet).

- [ ] **Step 5: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Add BatchEntry and _process_batch (fused compile+fit)"
```

---

## Task 7: `_ingest!` + `_offer_cv!`

**Files:**
- Modify: `src/identify_rate_equation.jl` (add both; near `_process_batch`)
- Test: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "_ingest! and cv pool" begin
    mk(n, loss, h) = EnzymeRates.BatchEntry(
        first(EnzymeRates.init_mechanisms(@enzyme_reaction begin
            substrates:S[C]; products:P[C] end)),
        n, loss, hash(h),
        (n_params=n, loss=loss, mechanism_type="M",
         rate_equation="v", fitted_param_names=(:K,),
         fitted_param_values=(1.0,), eq_hash=string(hash(h),base=16)))
    frontier = Dict{Int,Vector{EnzymeRates.BatchEntry}}()
    cv_pool  = Dict{Int,Vector{EnzymeRates.BatchEntry}}()
    best     = Dict{Int,Float64}()
    # two distinct equations + one duplicate-eq with worse loss, n_cv=2
    EnzymeRates._ingest!(frontier, cv_pool, best,
        [mk(5,2.0,:a), mk(5,1.0,:b), mk(5,3.0,:a)]; n_cv_candidates=2)
    @test length(frontier[5]) == 3            # frontier keeps ALL
    @test best[5] == 1.0                       # running min
    @test length(cv_pool[5]) == 2              # bounded, distinct eq_hash
    # the kept :a entry is the lower-loss one (2.0, not 3.0):
    a = only(filter(e -> e.eq_hash == string(hash(:a),base=16), cv_pool[5]))
    @test a.loss == 2.0
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `julia --project -e 'using EnzymeRates;
EnzymeRates._ingest!(Dict(),Dict(),Dict(),[]; n_cv_candidates=2)'`
Expected: FAIL — `UndefVarError: _ingest!`.

- [ ] **Step 3: Implement**

Add to `src/identify_rate_equation.jl`:

```julia
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
```

- [ ] **Step 4: Verify it passes**

Run the full suite: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -20`
Expected: the `_ingest!` testset passes.

- [ ] **Step 5: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Add _ingest! / _offer_cv! (bounded CV pool)"
```

---

## Task 8: `expand_mechanisms` returns a flat vector

**Files:**
- Modify: `src/mechanism_enumeration.jl:1673-1726` (`expand_mechanisms`, `_add_expansions_mech!`; delete `_push_mech!`)
- Modify: `src/identify_rate_equation.jl:775-781` (adapt the OLD `_beam_search` call site to the flat return — temporary, replaced in Task 9)
- Modify: `test/test_mechanism_enumeration.jl` (the ~19 `expand_mechanisms` call sites that treat the result as a `Dict`)

This task keeps `_n_fit_params_estimate` alive (the old `_beam_search` and `test_mechanism_enumeration.jl` still bucket by it); it is deleted in Task 10.

- [ ] **Step 1: Update the enumeration tests to the flat contract first (they encode it)**

In `test/test_mechanism_enumeration.jl`, every place that consumes `expand_mechanisms(...)` as a `Dict` must switch to a flat vector. The integration helper around `:156-180` (it iterates buckets / calls `_n_fit_params_estimate` to re-bucket) and the per-move testsets (search for `expand_mechanisms(`) currently do e.g.:

```julia
result = EnzymeRates.expand_mechanisms([am], rxn)   # was Dict
all_children = reduce(vcat, values(result); init=Union{...}[])
```

Replace each with the flat form:

```julia
all_children = EnzymeRates.expand_mechanisms([am], rxn)   # now a Vector
```

For assertions that grouped by bucket key (`result[k]`), regroup locally by `_n_fit_params_estimate` only where the test genuinely checks bucketing — otherwise assert directly on the flat vector. (These tests still use `_n_fit_params_estimate`; that's fine, it's deleted in Task 10 along with those assertions.) Run a grep to find them all:

Run: `grep -n "expand_mechanisms" test/test_mechanism_enumeration.jl`
Edit each of the ~19 sites to the flat contract.

- [ ] **Step 2: Run to confirm the tests now fail against the current Dict implementation**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: FAIL in `test_mechanism_enumeration.jl` — the flat-contract assertions don't match the current `Dict` return.

- [ ] **Step 3: Implement the flat `expand_mechanisms`**

In `src/mechanism_enumeration.jl`, replace `expand_mechanisms`/`_add_expansions_mech!` (`:1686-1718`) and DELETE `_push_mech!` (`:1720-1726`):

```julia
"""
    expand_mechanisms(mechs, reaction) -> Vector{Union{Mechanism, AllostericMechanism}}

Apply all expansion moves (RE→SS, split kinetic group, add dead-end
regulator, to-allosteric, add allosteric regulator, change allo state) to
each input mechanism and return the children as a flat vector. Bucketing by
parameter count is the caller's job, not enumeration's.
"""
function expand_mechanisms(
    mechs::Vector{<:Union{Mechanism, AllostericMechanism}},
    rxn::EnzymeReaction)
    result = Union{Mechanism, AllostericMechanism}[]
    for m in mechs
        _add_expansions_mech!(result, m, rxn)
    end
    result
end

function _add_expansions_mech!(
    result::Vector{Union{Mechanism, AllostericMechanism}},
    m::Union{Mechanism, AllostericMechanism},
    rxn::EnzymeReaction)
    for s in _expand_re_to_ss(m);                push!(result, s); end
    for s in _expand_split_kinetic_group(m);     push!(result, s); end
    for s in _expand_add_dead_end_regulator(m, rxn); push!(result, s); end
    for s in _expand_to_allosteric(m, rxn);      push!(result, s); end
    for s in _expand_add_allosteric_regulator(m, rxn); push!(result, s); end
    for s in _expand_change_allo_state(m);       push!(result, s); end
end
```

- [ ] **Step 4: Adapt the OLD `_beam_search` call site (temporary)**

In `src/identify_rate_equation.jl`, the old `_beam_search` (`:775-781`) iterates the Dict. Replace those lines so it consumes the flat vector but keeps the old estimate-bucketing behavior (this whole function is rewritten in Task 9):

```julia
        for child in expand_mechanisms(beam_mechs, prob.reaction)
            pc = _n_fit_params_estimate(child)
            pc > max_param_count && continue
            push!(get!(cache, pc,
                       Union{Mechanism, AllostericMechanism}[]), child)
        end
        dedup!(cache)
```

- [ ] **Step 5: Run the full suite to green**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: PASS — enumeration tests match the flat contract; the old `_beam_search` still works via inline estimate-bucketing.

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl src/identify_rate_equation.jl test/test_mechanism_enumeration.jl
git commit -m "expand_mechanisms returns a flat vector"
```

---

## Task 9: Rewrite `_beam_search`; new CSV schema; mandatory `save_dir`

**Files:**
- Modify: `src/identify_rate_equation.jl` — `_beam_search` (`:639-787`), `_rows_to_dataframe` (`:228-236`, drop a column), `identify_rate_equation` (`:161-211`, save_dir default + pass `n_cv_candidates`)
- Modify: `test/test_identify_rate_equation.jl` — replace the integration/save_dir/`_save_level_csv` testsets with new-schema versions; write the integration test first

- [ ] **Step 1: Migrate the existing tests to the new schema first**

These edits all live in `test/test_identify_rate_equation.jl`. The pipeline is already run once as `results` (`:171-180`, `save_dir`, `pmap_function=map`, `optimizer=pycma_opt`, `max_param_count=8`). Reuse that run.

(a) **`_rows_to_dataframe` testset (`:137-155`):** remove the `fit_inherited_from_estimate = missing,` field from the row literal (`:146`) and replace the column assertion (`:154`) `@test "fit_inherited_from_estimate" in names(df)` with:
```julia
        @test !("fit_inherited_from_estimate" in names(df))
```

(b) **"CSV output" testset (`:276-335`):** replace its whole body (the `params_estimate_*` checks and the `fit_inherited_from_estimate` chain) with:
```julia
    @testset "CSV output (new schema)" begin
        files = sort(filter(f -> endswith(f, ".csv"), readdir(save_dir)))
        @test "initial_mechanisms.csv" in files
        @test !any(startswith(f, "params_estimate_") for f in files)
        iters = filter(f -> startswith(f, "equation_search_iteration_"), files)
        @test !isempty(iters)
        nums = sort(parse.(Int, replace.(iters,
            "equation_search_iteration_" => "", ".csv" => "")))
        @test nums == collect(1:length(nums))      # sequential, no gaps
        init_df = CSV.read(joinpath(save_dir, "initial_mechanisms.csv"), DataFrame)
        @test nrow(init_df) == length(EnzymeRates._dedup_flat(
            collect(EnzymeRates.init_mechanisms(prob.reaction))))
        @test "eq_hash" in names(init_df)
        @test !("fit_inherited_from_estimate" in names(init_df))
        for f in files
            df_file = CSV.read(joinpath(save_dir, f), DataFrame)
            @test "eq_hash" in names(df_file)
            @test all(.!ismissing.(df_file.eq_hash))
            @test all(length.(string.(df_file.eq_hash)) .== 16)
            @test all(<=(8), df_file.n_params)     # max_param_count=8
        end
    end
```
The "save_dir non-empty check" testset (`:337-354`) is unchanged — it re-runs `identify_rate_equation` with the now-populated `save_dir` and still expects the `ErrorException` guard. The "results structure" testset (`:204-269`) is unchanged and already covers the LOOCV distinct-`eq_hash` invariant (`allunique(gdf.eq_hash)`), i.e. "N different equations" — keep it. The `_save_level_csv` testset (`:391-409`) is deleted in Task 11.

- [ ] **Step 2: Run it to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: FAIL — old `_beam_search` writes `params_estimate_*`, not the new names; `fit_inherited_from_estimate` column still present.

- [ ] **Step 3: Drop the `fit_inherited_from_estimate` column**

In `_rows_to_dataframe` (`src/identify_rate_equation.jl:228-236`) remove the `fit_inherited_from_estimate` line so the `DataFrame(...)` call is:

```julia
    df = DataFrame(
        n_params = [r.n_params for r in rows],
        loss = [r.loss for r in rows],
        mechanism_type = [r.mechanism_type for r in rows],
        rate_equation = [r.rate_equation for r in rows],
        eq_hash = [r.eq_hash for r in rows],
    )
```

- [ ] **Step 4: Rewrite `_beam_search`**

Replace the entire `_beam_search` function (`src/identify_rate_equation.jl:639-787`) with:

```julia
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
        group = BatchEntry[]
        for c in collect(keys(frontier))
            c <= target && append!(group, pop!(frontier, c))
        end

        to_expand = BatchEntry[]
        for c in unique(e.n_params for e in group)
            ec  = [e for e in group if e.n_params == c]
            sel = _select_beam([e.loss for e in ec];
                loss_rel_threshold, loss_abs_threshold,
                min_beam_width, best_override = best_loss_by_count[c])
            append!(to_expand, ec[sel])
        end

        if !isempty(to_expand)
            children = _dedup_flat(expand_mechanisms(
                [e.mech for e in to_expand], prob.reaction))
            child_entries = _process_batch(children, prob;
                pmap_function, optimizer, max_param_count, kwargs...)
            if !isempty(child_entries)
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
```

- [ ] **Step 5: Make `save_dir` mandatory-with-default and pass `n_cv_candidates`**

In `identify_rate_equation` (`src/identify_rate_equation.jl:161-211`):

Change the signature default:
```julia
    save_dir::String = _default_save_dir(),
```
Change the guard (no longer `!== nothing`):
```julia
    if isdir(save_dir)
        existing = filter(
            f -> endswith(f, ".csv"),
            readdir(save_dir))
        isempty(existing) || error(
            "save_dir already contains CSV " *
            "files. Use an empty directory " *
            "to avoid mixing results.")
    end
```
Pass `n_cv_candidates` into `_beam_search`:
```julia
    mechanisms, df = _beam_search(prob;
        min_beam_width, loss_rel_threshold,
        loss_abs_threshold,
        max_param_count, save_dir,
        pmap_function, optimizer, n_cv_candidates,
        fitting_kwargs...)
```

- [ ] **Step 6: Run the full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -40`
Expected: PASS for `test_identify_rate_equation.jl` (new artifacts) — note `_n_fit_params_estimate`, the canonical-hash helpers, `_CachedFitResult`, `_CompiledMechanismResult`, `_save_level_csv`, and `_project_cached_params` are now UNUSED by `src/` but still defined (deleted in Tasks 10–11). They do not break anything yet.

- [ ] **Step 7: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Rewrite _beam_search: advancing-target sweep, bounded CV pool, new CSVs"
```

---

## Task 10: Delete `_n_fit_params_estimate` and its tests

**Files:**
- Modify: `src/mechanism_enumeration.jl` (delete both overloads `:1070-1119`; update the `expand_mechanisms` docstring at `:1680` if it still mentions the estimate)
- Modify: `test/test_mechanism_enumeration.jl` (delete every `_n_fit_params_estimate` testset / assertion and the stale comments)

- [ ] **Step 1: Delete the test references**

Run: `grep -n "_n_fit_params_estimate" test/test_mechanism_enumeration.jl`
Delete each: the per-move delta assertions (`est(r) == est(am)+1`, the `deltas = sort([est(r)-est(am) …])` blocks), the floor-applying sites (`max(_n_fit_params_estimate(…), floor_pc)`), the lower-bound/`== n_actual` property tests, and the explanatory comments that reference it. If a whole `@testset` exists only to check the estimate, delete the testset. If a testset checks a real expansion-move property AND incidentally used the estimate for the count, keep the testset but replace the estimate call with the exact count `length(EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(m)))` ONLY where the test still asserts something meaningful about actual params; otherwise delete the assertion.

Note on termination coverage: do NOT add an exact-count "+1 per move" monotonicity test — it would FAIL (Δ=0 is common, by design). Termination is covered by the Task 9 integration test (a full LDH-like run completes).

- [ ] **Step 2: Run to confirm the tests fail (estimate still defined, refs removed cleanly)**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -20`
Expected: PASS (the deletions are self-consistent; estimate is still defined but now unused). If any remaining test references `_n_fit_params_estimate`, the grep in Step 1 missed it — fix.

- [ ] **Step 3: Delete the function**

In `src/mechanism_enumeration.jl`, delete both `_n_fit_params_estimate` overloads (the `Mechanism` one and the `AllostericMechanism` one, `:1070-1119` including their shared docstring). If the `expand_mechanisms` docstring still says "bucket results by their `_n_fit_params_estimate` value," it was already rewritten in Task 8 — confirm no lingering mention.

- [ ] **Step 4: Verify nothing references it**

Run: `grep -rn "_n_fit_params_estimate" src/ test/`
Expected: no output.

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -20`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Delete _n_fit_params_estimate and its tests"
```

---

## Task 11: Delete the canonical-hash machinery, projection, and dead structs

**Files:**
- Modify: `src/identify_rate_equation.jl` — delete `_expr_canonical_via_name_map` (`:348`), `_canonicalize_for_hash` (`:419`), `_build_name_map` (`:446`), `_dep_exprs_canonical` (`:481`), `_synth_dep_a_names` (`:503`), `_canonical_rate_eq_hash_data` (`:532`), `_canonical_rate_eq_hash` (`:543`), `_project_cached_params` (`:611`), `_CachedFitResult` (`:557`), `_CompiledMechanismResult` + `_CompiledMechanismFailure` (`:571-587`), `_save_level_csv` (`:256-265`)
- Modify: `src/rate_eq_derivation.jl:1168-1173` (the doc comment referencing `_expr_canonical_via_name_map`)
- Modify: `test/test_identify_rate_equation.jl` — delete the `_canonical_rate_eq_hash` testset (`:840-854`)

- [ ] **Step 1: Delete the obsolete tests**

In `test/test_identify_rate_equation.jl`, delete:
- the testset that calls `EnzymeRates._canonical_rate_eq_hash` (around `:840-854`);
- the `@testset "save_level_csv uses estimate-level filename"` (around `:391-409`) — it tests the deleted `_save_level_csv` and uses the dropped `fit_inherited_from_estimate` field.

- [ ] **Step 2: Delete the source functions/structs**

Delete, in `src/identify_rate_equation.jl`, the functions and structs listed above. They form one contiguous-ish canonicalization block plus the two records and the old CSV writer. After deleting, the file should have no `name_map`, `canon_to_rep`, `_canonical*`, `_build_name_map`, `_dep_exprs_canonical`, `_synth_dep_a_names`, `_expr_canonical_via_name_map`, `_project_cached_params`, `_CachedFitResult`, `_CompiledMechanism*`, or `_save_level_csv`.

- [ ] **Step 3: Fix the stale doc comment in `rate_eq_derivation.jl`**

The comment at `src/rate_eq_derivation.jl:1168-1173` explains why over-emitting name_map entries is harmless "because downstream `_expr_canonical_via_name_map` only substitutes Symbols that actually appear…". That helper is gone. Rewrite the comment to describe the current behavior without referencing the deleted function (the enumeration still intentionally over-emits I-state names; just drop the `_expr_canonical_via_name_map` sentence, or replace it with: "Unused names are inert — nothing downstream consumes a name that does not appear in the rate-equation Exprs.").

- [ ] **Step 4: Verify no references remain**

Run: `grep -rn "_canonical_rate_eq_hash\|_project_cached_params\|_CachedFitResult\|_CompiledMechanism\|_save_level_csv\|_build_name_map\|_dep_exprs_canonical\|_synth_dep_a_names\|_canonicalize_for_hash\|_expr_canonical_via_name_map" src/ test/`
Expected: no output.

- [ ] **Step 5: Run the full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: PASS. Aqua's stale-deps / undefined-export checks stay green.

- [ ] **Step 6: Commit**

```bash
git add src/identify_rate_equation.jl src/rate_eq_derivation.jl test/test_identify_rate_equation.jl
git commit -m "Delete canonical-hash machinery, projection, and dead records"
```

---

## Task 12: Doc sync (`.claude/CLAUDE.md`)

**Files:**
- Modify: `.claude/CLAUDE.md`

- [ ] **Step 1: Update the enumeration + identify descriptions**

In `.claude/CLAUDE.md`:
- The `expand_mechanisms` description (the "Mechanism enumeration building blocks" section and the Source Layout entry for `mechanism_enumeration.jl`) says it returns `Dict{Int, Vector{...}}` keyed by `_n_fit_params_estimate`. Change to: returns `Vector{Union{Mechanism, AllostericMechanism}}`; bucketing by **actual** fitted-param count now lives in `identify_rate_equation.jl`'s beam search.
- Remove `_n_fit_params_estimate` from the building-blocks list and from `_push_mech!`/floor mentions.
- In the `identify_rate_equation.jl` Source Layout entry, replace the canonical-rate-equation-hashing description with: structural dedup (`_dedup_flat`) + a comment-stripped `_rate_eq_dedup_key` used as a CSV tag and the LOOCV distinct-equation key; the beam search is an advancing-target sweep over actual param counts with mandatory CSV output (`initial_mechanisms.csv` + `equation_search_iteration_N.csv`) and `max_param_count` capping **actual fitted** params.
- Remove references to the deleted `_canonical_rate_eq_hash` / `_project_cached_params` as the dedup/fit-reuse mechanism.

- [ ] **Step 2: Verify**

Run: `grep -n "_n_fit_params_estimate\|_canonical_rate_eq_hash\|_project_cached_params\|Dict{Int, Vector" .claude/CLAUDE.md`
Expected: no stale references (or only historical ones you intentionally keep — there should be none).

- [ ] **Step 3: Commit**

```bash
git add .claude/CLAUDE.md
git commit -m "Doc sync: flat expand_mechanisms, actual-count beam search, string dedup"
```

---

## Task 13: Full-suite green + end-to-end acceptance

**Files:** none (verification)

- [ ] **Step 1: Full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -40`
Expected: pass count ≥ Task-0 baseline minus the intentionally-deleted tests; 0 fail. Perf gate (`test_rate_eq_derivation.jl`) and compile-budget gate green.

- [ ] **Step 2: End-to-end artifact check on LDH**

Run:
```bash
julia --project -e '
using EnzymeRates, DataFrames, CSV; const ER=EnzymeRates
rxn = ER.@enzyme_reaction begin
    substrates:NADH[N], Pyruvate[C]
    products:Lactate[C], NAD[N]
    oligomeric_state:4
end
# tiny stub data: reuse a fixture or a few rows with a `group` column.
# Build prob + run with a small cap and pmap_function=map, then:
#   readdir(save_dir)  ==>  initial_mechanisms.csv + equation_search_iteration_1.csv, _2.csv, …
#   every CSV n_params <= cap ; initial row count == length(_dedup_flat(init_mechanisms(rxn)))
'
```
Expected: `initial_mechanisms.csv` present with all deduped init mechanisms at their real 5–6 params; `equation_search_iteration_N.csv` sequential with no gaps; no `params_estimate_*.csv`; `n_params` never exceeds the cap.

- [ ] **Step 3: Grep clean-up confirmation**

Run: `grep -rn "_n_fit_params_estimate\|_canonical_rate_eq_hash\|_project_cached_params\|params_estimate_\|_save_level_csv\|fit_inherited" src/ test/`
Expected: no output.

- [ ] **Step 4: Final commit (if any stray fixups)**

```bash
git add -A
git commit -m "Final cleanup for actual-param-count beam search" || echo "nothing to commit"
```

---

## Self-Review checklist (run before handing off)

- **Spec coverage:** every spec section maps to a task — `_rate_eq_dedup_key` (T1), `_dedup_flat` (T2), `_default_save_dir` (T3), CSV writers (T4), `_select_beam` override (T5), `_process_batch` (T6), `_ingest!`/cv_pool (T7), flat `expand_mechanisms` (T8), `_beam_search`/`_rows_to_dataframe`/`identify_rate_equation` (T9), delete estimate (T10), delete canonical machinery (T11), docs (T12), verify (T13).
- **No new `Sig` types** introduced — derivation untouched. Perf + compile-budget gates verified in T0 and T13.
- **Type consistency:** `BatchEntry` fields (`mech`, `n_params`, `loss`, `eq_hash::UInt64`, `row`) are used identically in `_process_batch`, `_ingest!`, `_offer_cv!`, and `_beam_search`; the row NamedTuple schema (`n_params`, `loss`, `mechanism_type`, `rate_equation`, `fitted_param_names`, `fitted_param_values`, `eq_hash`) matches `_rows_to_dataframe` after the column drop.
