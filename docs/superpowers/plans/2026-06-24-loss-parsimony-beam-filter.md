# Loss-parsimony beam filter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `loss_parsimony_threshold` beam knob that stops expanding a mechanism whose loss fails to come within a factor of the best model with one fewer parameter.

**Architecture:** The threshold tightens the existing per-count beam cutoff inside `_select_beam` via `min(...)`, leaving `min_beam_width` as a guaranteed OR floor. The keyword threads `identify_rate_equation → _beam_search → _select_beam` exactly like the two existing loss thresholds. The call site supplies the pre-multiplied cutoff `loss_parsimony_threshold * best_loss_by_count[c-1]`, which is already in scope. No new data structures, no sensitivity Jacobian.

**Tech Stack:** Julia, `Optimization.jl` / `OptimizationCMAEvolutionStrategy`, `Test`. All changes live in `src/identify_rate_equation.jl`, `test/test_identify_rate_equation.jl`, and two docs tutorial pages.

## Global Constraints

- Spec: `docs/superpowers/specs/2026-06-24-loss-parsimony-beam-filter-design.md`.
- Public keyword name is exactly `loss_parsimony_threshold`; default `1.01`.
- Internal `_select_beam` parameter name is exactly `parsimony_cutoff` (the pre-multiplied cutoff value, `Union{Nothing,Float64}`, default `nothing`).
- `min_beam_width` MUST stay a guaranteed floor: the top `min_beam_width` by loss always survive, even when zero mechanisms clear the loss cutoff.
- The filter is an expansion gate only — it must not touch `cv_pool` or model selection.
- 92-character line length, 4-space indentation; match surrounding style.
- Run the full suite green before declaring a code task done: `julia --project -e 'using Pkg; Pkg.test()'`.

---

### Task 1: `_select_beam` gains `parsimony_cutoff`

**Files:**
- Modify: `src/identify_rate_equation.jl` — `_select_beam` docstring (`296-305`), signature (`306-312`), body (`318`).
- Test: `test/test_identify_rate_equation.jl` — new testset after the `_select_beam best_override` testset (`888-899`).

**Interfaces:**
- Consumes: nothing new.
- Produces: `_select_beam(losses; loss_rel_threshold::Float64, loss_abs_threshold::Float64, min_beam_width::Int, best_override::Union{Nothing,Float64}=nothing, parsimony_cutoff::Union{Nothing,Float64}=nothing) -> Vector{Int}`. When `parsimony_cutoff !== nothing`, the qualifying cutoff becomes `min(loss_rel_threshold*best + loss_abs_threshold, parsimony_cutoff)`; the `min_beam_width` rank floor is unchanged.

- [ ] **Step 1: Write the failing tests**

In `test/test_identify_rate_equation.jl`, immediately after the closing `end` of the `@testset "_select_beam best_override"` block (currently line 899), add:

```julia
@testset "_select_beam parsimony_cutoff" begin
    # Floor guarantee: a parsimony_cutoff below every loss admits nothing via
    # the loss filter, yet min_beam_width still keeps the top-k by loss.
    losses = [1.0, 1.5, 2.5, 5.0, 10.0]
    @test EnzymeRates._select_beam(losses;
        loss_rel_threshold=2.0, loss_abs_threshold=0.0,
        min_beam_width=2, parsimony_cutoff=0.5) == [1, 2]

    # Tightening: a parsimony_cutoff stricter than the rel/abs cutoff drops the
    # mechanisms between the two cutoffs (min_beam_width=1 so the floor
    # re-admits only the single best). Without it, rel=10 would admit all four.
    losses = [1.0, 1.5, 2.5, 5.0]
    @test EnzymeRates._select_beam(losses;
        loss_rel_threshold=10.0, loss_abs_threshold=0.0,
        min_beam_width=1, parsimony_cutoff=2.0) == [1, 2]

    # No-op: parsimony_cutoff=nothing reproduces the parsimony-free selection.
    kw = (loss_rel_threshold=2.0, loss_abs_threshold=0.0, min_beam_width=1)
    @test EnzymeRates._select_beam(losses; kw..., parsimony_cutoff=nothing) ==
          EnzymeRates._select_beam(losses; kw...)

    # Interaction: min() picks the smaller cutoff. With best_override=2.0 the
    # rel cutoff is 2.4 (admits 1,2); a tighter parsimony_cutoff=1.0 overrides
    # it down to just the single best.
    losses = [1.0, 1.5, 3.0]
    ov = (loss_rel_threshold=1.2, loss_abs_threshold=0.0,
          min_beam_width=1, best_override=2.0)
    @test EnzymeRates._select_beam(losses; ov...) == [1, 2]
    @test EnzymeRates._select_beam(losses; ov..., parsimony_cutoff=1.0) == [1]
end
```

- [ ] **Step 2: Run the tests to verify they fail**

Run:
```bash
julia --project=. -e '
using EnzymeRates, Test
@testset "parsimony red" begin
    @test EnzymeRates._select_beam([1.0,1.5,2.5,5.0,10.0];
        loss_rel_threshold=2.0, loss_abs_threshold=0.0,
        min_beam_width=2, parsimony_cutoff=0.5) == [1, 2]
end'
```
Expected: FAIL — a `MethodError`/keyword error, because `_select_beam` does not yet accept `parsimony_cutoff`.

- [ ] **Step 3: Add the `parsimony_cutoff` keyword to the signature**

In `src/identify_rate_equation.jl`, change the `_select_beam` signature:

```julia
function _select_beam(
    losses::AbstractVector{<:Real};
    loss_rel_threshold::Float64,
    loss_abs_threshold::Float64,
    min_beam_width::Int,
    best_override::Union{Nothing,Float64}=nothing,
)
```

to:

```julia
function _select_beam(
    losses::AbstractVector{<:Real};
    loss_rel_threshold::Float64,
    loss_abs_threshold::Float64,
    min_beam_width::Int,
    best_override::Union{Nothing,Float64}=nothing,
    parsimony_cutoff::Union{Nothing,Float64}=nothing,
)
```

- [ ] **Step 4: Tighten the cutoff in the body**

In the same function, change:

```julia
    cutoff = loss_rel_threshold * best + loss_abs_threshold
    selected = Int[]
```

to:

```julia
    cutoff = loss_rel_threshold * best + loss_abs_threshold
    parsimony_cutoff !== nothing && (cutoff = min(cutoff, parsimony_cutoff))
    selected = Int[]
```

- [ ] **Step 5: Update the `_select_beam` docstring**

Replace the docstring above `function _select_beam` (currently lines `296-305`):

```julia
"""
Return indices into `losses` for mechanisms that qualify for the
beam at this level. A mechanism qualifies if either:
  • its loss ≤ loss_rel_threshold * best_loss + loss_abs_threshold,
  • OR its rank (1-indexed by ascending loss) ≤ min_beam_width.

Mechanisms with non-finite losses (`Inf`, `NaN`) are excluded
unconditionally — they represent failed or non-converging fits
that should not propagate to the next level.
"""
```

with:

```julia
"""
Return indices into `losses` for mechanisms that qualify for the
beam at this level. A mechanism qualifies if either:
  • its loss ≤ cutoff, where
    cutoff = min(loss_rel_threshold * best_loss + loss_abs_threshold,
                 parsimony_cutoff) and the parsimony term is dropped
    when `parsimony_cutoff === nothing`,
  • OR its rank (1-indexed by ascending loss) ≤ min_beam_width.

`parsimony_cutoff` (the loss-parsimony threshold times the best loss
at one fewer parameter) only tightens the loss cutoff. `min_beam_width`
stays a hard floor: the top `min_beam_width` always qualify, even when
the loss cutoff admits fewer.

Mechanisms with non-finite losses (`Inf`, `NaN`) are excluded
unconditionally — they represent failed or non-converging fits
that should not propagate to the next level.
"""
```

- [ ] **Step 6: Run the tests to verify they pass**

Run:
```bash
julia --project=. -e '
using EnzymeRates, Test
@testset "parsimony green" begin
    @test EnzymeRates._select_beam([1.0,1.5,2.5,5.0,10.0];
        loss_rel_threshold=2.0, loss_abs_threshold=0.0,
        min_beam_width=2, parsimony_cutoff=0.5) == [1, 2]
    @test EnzymeRates._select_beam([1.0,1.5,2.5,5.0];
        loss_rel_threshold=10.0, loss_abs_threshold=0.0,
        min_beam_width=1, parsimony_cutoff=2.0) == [1, 2]
    kw = (loss_rel_threshold=2.0, loss_abs_threshold=0.0, min_beam_width=1)
    @test EnzymeRates._select_beam([1.0,1.5,2.5,5.0]; kw..., parsimony_cutoff=nothing) ==
          EnzymeRates._select_beam([1.0,1.5,2.5,5.0]; kw...)
    ov = (loss_rel_threshold=1.2, loss_abs_threshold=0.0,
          min_beam_width=1, best_override=2.0)
    @test EnzymeRates._select_beam([1.0,1.5,3.0]; ov...) == [1, 2]
    @test EnzymeRates._select_beam([1.0,1.5,3.0]; ov..., parsimony_cutoff=1.0) == [1]
end'
```
Expected: PASS — `Test Summary: parsimony green | Pass 5`.

- [ ] **Step 7: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "feat: _select_beam honors a parsimony_cutoff that tightens the beam"
```

---

### Task 2: Thread `loss_parsimony_threshold` through the search

**Files:**
- Modify: `src/identify_rate_equation.jl` — docstring signature line (`97`), kwarg doc list (`109-113`), "Beam selection" section (`140-148`), `identify_rate_equation` signature (`176-178`), `_beam_search` call (`210-215`), `_beam_search` signature (`494-499`), `_select_beam` call site (`543-545`).
- Test: `test/test_identify_rate_equation.jl` — new testset after the `"identify runs on a solver that rejects popsize"` testset (`1019-1040`).

**Interfaces:**
- Consumes: `_select_beam(...; parsimony_cutoff=...)` from Task 1.
- Produces: `identify_rate_equation(prob; ..., loss_parsimony_threshold::Float64 = 1.01, ...)` and `_beam_search(prob; ..., loss_parsimony_threshold, ...)`. The call site passes `parsimony_cutoff = haskey(best_loss_by_count, c-1) ? loss_parsimony_threshold * best_loss_by_count[c-1] : nothing`.

- [ ] **Step 1: Write the failing threading test**

In `test/test_identify_rate_equation.jl`, immediately after the closing `end` of `@testset "identify runs on a solver that rejects popsize"` (currently line 1040), add:

```julia
@testset "loss_parsimony_threshold threads through identify_rate_equation" begin
    # An unknown keyword throws at the call boundary (see the removed-kwargs
    # test), so a clean end-to-end run with an explicit non-default value
    # proves the keyword is accepted and forwarded to the beam.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    data = (group = ["G1", "G1", "G2", "G2"],
            Rate = [0.5, 0.8, 1.0, 1.1],
            S = [1.0, 2.0, 3.0, 4.0],
            P = [0.1, 0.2, 0.3, 0.4])
    prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)
    results = identify_rate_equation(prob;
        optimizer=CMAEvolutionStrategyOpt(),
        min_beam_width=1, loss_rel_threshold=1.0, loss_abs_threshold=0.0,
        loss_parsimony_threshold=2.0,
        max_param_count=6, n_cv_candidates=1, n_restarts=1, maxtime=1.0,
        save_dir=mktempdir(), show_progress=false)
    @test results isa IdentifyRateEquationResults
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```bash
julia --project=. -e '
using EnzymeRates, Test, DataFrames, CSV, Random, Statistics, OptimizationCMAEvolutionStrategy
rxn = EnzymeRates.@enzyme_reaction begin
    substrates: S[C]
    products: P[C]
end
data = (group=["G1","G1","G2","G2"], Rate=[0.5,0.8,1.0,1.1],
        S=[1.0,2.0,3.0,4.0], P=[0.1,0.2,0.3,0.4])
prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)
identify_rate_equation(prob; optimizer=CMAEvolutionStrategyOpt(),
    min_beam_width=1, loss_rel_threshold=1.0, loss_abs_threshold=0.0,
    loss_parsimony_threshold=2.0, max_param_count=6, n_cv_candidates=1,
    n_restarts=1, maxtime=1.0, save_dir=mktempdir(), show_progress=false)'
```
Expected: FAIL — an error that `loss_parsimony_threshold` is not a recognized keyword of `identify_rate_equation`.

- [ ] **Step 3: Add the keyword to `identify_rate_equation`'s signature**

In `src/identify_rate_equation.jl`, change:

```julia
    loss_rel_threshold::Float64 = 2.0,
    loss_abs_threshold::Float64 = 0.01,
    max_param_count::Int = 20,
```

to:

```julia
    loss_rel_threshold::Float64 = 2.0,
    loss_abs_threshold::Float64 = 0.01,
    loss_parsimony_threshold::Float64 = 1.01,
    max_param_count::Int = 20,
```

- [ ] **Step 4: Forward it to `_beam_search`**

Change the `_beam_search` call:

```julia
    mechanisms, df = _beam_search(prob;
        min_beam_width, loss_rel_threshold,
        loss_abs_threshold,
        max_param_count, save_dir, show_progress,
        optimizer, n_cv_candidates,
        fitting_kwargs...)
```

to:

```julia
    mechanisms, df = _beam_search(prob;
        min_beam_width, loss_rel_threshold,
        loss_abs_threshold, loss_parsimony_threshold,
        max_param_count, save_dir, show_progress,
        optimizer, n_cv_candidates,
        fitting_kwargs...)
```

- [ ] **Step 5: Accept it in `_beam_search`'s signature**

Change:

```julia
function _beam_search(
    prob::IdentifyRateEquationProblem;
    min_beam_width, loss_rel_threshold, loss_abs_threshold,
    max_param_count, save_dir, show_progress,
    optimizer, n_cv_candidates, kwargs...
)
```

to:

```julia
function _beam_search(
    prob::IdentifyRateEquationProblem;
    min_beam_width, loss_rel_threshold, loss_abs_threshold,
    loss_parsimony_threshold,
    max_param_count, save_dir, show_progress,
    optimizer, n_cv_candidates, kwargs...
)
```

- [ ] **Step 6: Pass the pre-multiplied cutoff at the call site**

Change the per-count `_select_beam` call:

```julia
            sel = _select_beam([e.loss for e in entries_at_count];
                loss_rel_threshold, loss_abs_threshold,
                min_beam_width, best_override = best_loss_by_count[c])
```

to:

```julia
            sel = _select_beam([e.loss for e in entries_at_count];
                loss_rel_threshold, loss_abs_threshold,
                min_beam_width, best_override = best_loss_by_count[c],
                parsimony_cutoff = haskey(best_loss_by_count, c - 1) ?
                    loss_parsimony_threshold * best_loss_by_count[c - 1] :
                    nothing)
```

- [ ] **Step 7: Update the `identify_rate_equation` docstring**

(a) In the docstring signature block, change:

```julia
        min_beam_width=50, loss_rel_threshold=2.0, loss_abs_threshold=0.01,
        max_param_count=20, n_restarts=20, maxtime=60.0, maxiters=10_000_000,
```

to:

```julia
        min_beam_width=50, loss_rel_threshold=2.0, loss_abs_threshold=0.01,
        loss_parsimony_threshold=1.01,
        max_param_count=20, n_restarts=20, maxtime=60.0, maxiters=10_000_000,
```

(b) In the keyword list, change:

```julia
- `loss_abs_threshold::Float64 = 0.01`: absolute tolerance
  for beam selection
- `max_param_count::Int = 20`: stop expanding beyond
```

to:

```julia
- `loss_abs_threshold::Float64 = 0.01`: absolute tolerance
  for beam selection
- `loss_parsimony_threshold::Float64 = 1.01`: a mechanism
  keeps expanding only if its loss is within this factor of
  the best model with one fewer parameter — an added
  parameter must earn its keep. Combined with the other loss
  thresholds via `min`; `min_beam_width` stays a hard floor.
  `Inf` disables it.
- `max_param_count::Int = 20`: stop expanding beyond
```

(c) Replace the "Beam selection" section:

```markdown
# Beam selection

A mechanism qualifies for the next-level beam if either:
- its loss ≤ `loss_rel_threshold * best_loss + loss_abs_threshold`,
- OR its rank by loss (ascending) ≤ `min_beam_width`.

The additive term protects against `best_loss` approaching zero
(simulated / very-low-loss data) where a purely multiplicative
threshold would collapse the beam to the single best mechanism.
```

with:

```markdown
# Beam selection

A mechanism at parameter count `n` qualifies for the next-level
beam if either:
- its loss ≤ `min(loss_rel_threshold * best(n) + loss_abs_threshold,
  loss_parsimony_threshold * best(n-1))`, where `best(k)` is the
  lowest loss seen at parameter count `k`; the second term is
  dropped at the base count (no `n-1` level),
- OR its rank by loss (ascending) ≤ `min_beam_width`. This floor
  always keeps the top `min_beam_width` mechanisms, even when the
  loss cutoff admits fewer.

The additive term protects against `best_loss` approaching zero
(simulated / very-low-loss data) where a purely multiplicative
threshold would collapse the beam to the single best mechanism.
```

- [ ] **Step 8: Run the threading test to verify it passes**

Run the same command as Step 2. Expected: PASS — the run completes and prints no error; the wrapping `@testset` reports the `results isa IdentifyRateEquationResults` test passing when run through the suite.

- [ ] **Step 9: Run the full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS — all testsets green, including the pre-existing `"beam selection: loss thresholds + min_beam_width floor"`, `"_select_beam best_override"`, `"_select_beam parsimony_cutoff"`, and the end-to-end recovery tests. The new active default (`1.01`) must not break any existing test; with the greedy `min_beam_width=1` settings those tests use, the parsimony term never tightens below the floor, so their results are unchanged.

- [ ] **Step 10: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "feat: add loss_parsimony_threshold beam knob (default 1.01)"
```

---

### Task 3: Sync the tutorial default lists

**Files:**
- Modify: `docs/src/getting_started.md` (`153-154`).
- Modify: `docs/src/identify/tutorial.md` (`11-12`).

**Interfaces:**
- Consumes: the default `loss_parsimony_threshold=1.01` introduced in Task 2.
- Produces: nothing code-facing; the two prose lists of production defaults now name the new knob.

- [ ] **Step 1: Update `getting_started.md`**

Change:

```markdown
(`min_beam_width=50`, `loss_rel_threshold=2.0`, `loss_abs_threshold=0.01`,
`max_param_count=20`) and would often run for many hours and require a High
```

to:

```markdown
(`min_beam_width=50`, `loss_rel_threshold=2.0`, `loss_abs_threshold=0.01`,
`loss_parsimony_threshold=1.01`, `max_param_count=20`) and would often run for
many hours and require a High
```

- [ ] **Step 2: Update `identify/tutorial.md`**

Change:

```markdown
widens the beam to the defaults (`min_beam_width=50`, `loss_rel_threshold=2.0`,
`loss_abs_threshold=0.01`, `max_param_count=20`) and would often run for many
```

to:

```markdown
widens the beam to the defaults (`min_beam_width=50`, `loss_rel_threshold=2.0`,
`loss_abs_threshold=0.01`, `loss_parsimony_threshold=1.01`, `max_param_count=20`)
and would often run for many
```

- [ ] **Step 3: Verify the edits**

Run:
```bash
grep -rn "loss_parsimony_threshold" docs/src/
```
Expected: both `docs/src/getting_started.md` and `docs/src/identify/tutorial.md` now list `loss_parsimony_threshold=1.01` in their production-default sentences.

The two edits are prose only (not inside `@example` blocks), so they cannot change any rendered output. The runnable `@example` blocks already omit the keyword and pick up the new `1.01` default automatically; because the filter never discards a level's best (the `min_beam_width` floor guarantees it), the recovered mechanisms those blocks render are unchanged. The CI docs job (`docs/make.jl`, `doctest = true`) validates the full build; a local check is `julia --project=docs docs/make.jl` if desired (heavy — needs the docs environment).

- [ ] **Step 4: Commit**

```bash
git add docs/src/getting_started.md docs/src/identify/tutorial.md
git commit -m "docs: list loss_parsimony_threshold among the production beam defaults"
```

---

## Self-Review

**Spec coverage:**
- The rule (`min()` tightening) → Task 1 (Steps 3-4) and Task 2 (Step 6 call site).
- The floor stays → Task 1 keeps the `min_beam_width` OR untouched; Task 1 Step 1 floor-guarantee test asserts it.
- Expansion-only / `cv_pool` untouched → no edit touches `cv_pool`; the call site only feeds `to_expand`.
- Default `1.01` → Task 2 Step 3.
- Exemptions (base / missing `n-1`) → Task 2 Step 6 `haskey(best_loss_by_count, c-1)` guard.
- Code changes (four edits) → Task 1 (`_select_beam`), Task 2 (threading + docstrings).
- Safety → covered by the running-min cutoff plus the floor; no code beyond the above.
- Tests (floor / tightening / no-op / interaction + threading) → Task 1 Step 1, Task 2 Step 1.
- Doc defaults sync (beyond spec, keeps docs accurate) → Task 3.

**Placeholder scan:** none — every step shows the exact code and command.

**Type consistency:** public keyword `loss_parsimony_threshold::Float64` and internal `parsimony_cutoff::Union{Nothing,Float64}` are used identically across Task 1 and Task 2; the call site multiplies the threshold by `best_loss_by_count[c-1]` to produce the `parsimony_cutoff` value, matching `_select_beam`'s signature.
