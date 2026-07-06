# Canonical Kinetic-Group Partition Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the `Mechanism` (and `AllostericMechanism`) constructor merge kinetic groups whose binding-K representatives are single-symbol Wegscheider-tied, so same-graph rate-equivalent mechanisms collapse to one struct before derivation or fitting, and remove the now-redundant derivation-time Pass-2 rename.

**Architecture:** The tie machinery (`_build_wegscheider_rename_map` → `_dependent_param_exprs_kernel`) already runs on a plain `Mechanism`. The public constructor builds a provisional (direction/order-canonical, un-merged) mechanism through a private `Val(:_raw)` inner constructor, computes the binding-K tie relation on it, union-find-merges the tied kinetic groups, and re-canonicalizes order. `Base.hash`/`==` on `(reaction, steps)` then make split and merged variants identical, so the existing `unique!` and the `eq_hash` fit memo collapse them.

**Tech Stack:** Julia 1.12, EnzymeRates.jl. Tests via `Pkg.test()` / `julia --project`.

## Global Constraints

- **92-character line length, 4-space indentation.** (`.claude/CLAUDE.md`)
- **Match surrounding code style.** Names describe purpose, not history; no "new"/"old"/"improved". Every source file starts with two `# ABOUTME:` lines (existing files already have them — do not add duplicates).
- **TDD, always.** Write the failing test, watch it fail, implement minimally, watch it pass, commit. Commit frequently.
- **The `rate_equation` 0-alloc / sub-100 ns runtime contract is untouchable.** All changes here are compile-time construction/derivation. If any change would make `rate_equation` allocate or slow down, STOP and raise it with Denis.
- **Canonical Step Form is load-bearing.** This plan adds partition canonicalization as a third construction-time invariant alongside the existing step-direction and step/group-order canonicalization. Do not relax the existing two.
- **Run tests before every commit:** `julia --project -e 'using Pkg; Pkg.test()'` (cold; pays precompile). For a single file during iteration: `julia --project -e 'using EnzymeRates; include("test/<file>.jl")'` inside a `@testset`, or the project's usual per-file runner.

---

## File structure

- **Modify `src/types.jl`** — `Mechanism` and `AllostericMechanism` constructors: extract the existing direction+order canonicalization into a helper, add a private `Val(:_raw)` provisional inner constructor, and call the partition merge from the public constructor.
- **Modify `src/rate_eq_derivation.jl`** — add a `Mechanism`-accepting method of `_build_wegscheider_rename_map`; add `_merge_tied_kinetic_groups`; later delete Pass-2. Add the allosteric merge alongside `_state_wegscheider_rename_map`.
- **Modify `test/test_types.jl`** — new tests for the provisional constructor, the merge helper, and the partition-merge invariant.
- **Modify `test/test_identify_rate_equation.jl`** — re-verify / update the 55-class assertion.
- **Modify `test/test_rate_eq_derivation.jl`, `test/test_allosteric_golden.jl`** — re-baseline any golden that a merged mechanism changes.

---

### Task 1: Extract direction/order canonicalization; add the provisional constructor

**Files:**
- Modify: `src/types.jl:553-563` (the `Mechanism` inner constructor)
- Test: `test/test_types.jl`

**Interfaces:**
- Produces: `Mechanism(reaction, steps, ::Val{:_raw})` — a `Mechanism` that is direction- and order-canonical but does **not** merge kinetic groups. Used by later tasks and tests to build split forms.
- Produces: `_canon_step_groups(reaction, steps)::Vector{Vector{Step}}` — the existing direction + order canonicalization, factored out.

- [ ] **Step 1: Write the failing test** (in `test/test_types.jl`, inside a new `@testset "provisional Mechanism ctor"`):

```julia
@testset "provisional Mechanism ctor" begin
    e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
    e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
    e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)
    rxn = @enzyme_reaction(begin
        substrates: S[C]
        products:   P[C]
    end)
    groups = [
        [EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)],
        [EnzymeRates.Step(e_s, e_p, nothing, false)],
        [EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)],
    ]
    raw  = EnzymeRates.Mechanism(rxn, deepcopy(groups), Val(:_raw))
    full = EnzymeRates.Mechanism(rxn, deepcopy(groups))
    # With no tied binding-K split, provisional and public agree.
    @test EnzymeRates.steps(raw) == EnzymeRates.steps(full)
    # Provisional is a real, usable Mechanism.
    @test length(EnzymeRates.steps(raw)) == 3
end
```

- [ ] **Step 2: Run it, verify it fails** with `MethodError` (no `Val{:_raw}` method).

- [ ] **Step 3: Implement.** Replace the `Mechanism` inner constructor block (`src/types.jl:556-562`) with:

```julia
    # Direction- and order-canonical, no partition merge. Used to compute
    # ties before merging, and by tests that build split forms.
    function Mechanism(reaction::EnzymeReaction,
                       steps::Vector{Vector{Step}}, ::Val{:_raw})
        new(reaction, _canon_step_groups(reaction, steps))
    end
    function Mechanism(reaction::EnzymeReaction,
                       steps::Vector{Vector{Step}})
        prov = Mechanism(reaction, steps, Val(:_raw))
        merged = _merge_tied_kinetic_groups(prov)
        merged === prov.steps ?
            prov : Mechanism(reaction, merged, Val(:_raw))
    end
```

Add the helper just above the `struct Mechanism` block (near `_canonical_group_order!`, `src/types.jl:517`):

```julia
# Direction + order canonicalization shared by both Mechanism inner
# constructors: canonicalize iso-step direction, sort steps within groups,
# order groups. Returns fresh vectors.
function _canon_step_groups(reaction::EnzymeReaction,
                            steps::Vector{Vector{Step}})
    groups = _canonicalize_iso_groups(reaction, steps)
    permute!(groups, _canonical_group_order!(groups))
    _assert_no_re_ss_duplicate(groups)
    groups
end
```

`_merge_tied_kinetic_groups` does not exist yet (Task 3). Julia resolves it at call time, and no `Mechanism` is constructed during module load, so this compiles. Confirm no module-load-time `Mechanism(...)` exists: `grep -rn "Mechanism(" src/ | grep -v function | grep -v "::"` and check none run at top level.

- [ ] **Step 4: Add a temporary stub so the suite loads** until Task 3. At the bottom of `src/rate_eq_derivation.jl`, add:

```julia
# Replaced with the real merge in the partition-canonicalization work.
_merge_tied_kinetic_groups(mech::Mechanism) = mech.steps
```

- [ ] **Step 5: Run the test, verify it passes.**

Run: `julia --project -e 'using EnzymeRates; include("test/test_types.jl")'`
Expected: the new testset passes; no other testset regresses.

- [ ] **Step 6: Commit.**

```bash
git add src/types.jl src/rate_eq_derivation.jl test/test_types.jl
git commit -m "Add provisional Mechanism constructor and _canon_step_groups helper"
```

---

### Task 2: Give `_build_wegscheider_rename_map` a `Mechanism` method

**Files:**
- Modify: `src/rate_eq_derivation.jl:117-143`
- Test: `test/test_types.jl` (or `test/test_rate_eq_derivation.jl`)

**Interfaces:**
- Produces: `_build_wegscheider_rename_map(mech::Mechanism)::Dict{Symbol,Symbol}` — the single-symbol binding-K tie relation, computed directly on a struct.
- Consumes: `_dependent_param_exprs_kernel(mech::Mechanism, rename)`, `_step_parameters(mech)`, `_flat_steps(mech)`, `name(p, mech)`, `is_binding`, `is_equilibrium` — all already `Mechanism`-based.

- [ ] **Step 1: Write the failing test:**

```julia
@testset "wegscheider rename map on Mechanism" begin
    m = first(EnzymeRates.init_mechanisms(@enzyme_reaction(begin
        substrates: A[C], B[N]
        products:   P[C], Q[N]
    end)))
    from_type = EnzymeRates._build_wegscheider_rename_map(
        EnzymeRates.compile_mechanism(m))
    from_mech = EnzymeRates._build_wegscheider_rename_map(m)
    @test from_type == from_mech
end
```

- [ ] **Step 2: Run it, verify it fails** with `MethodError` (no `Mechanism` method).

- [ ] **Step 3: Implement.** Change the head of `_build_wegscheider_rename_map` (`src/rate_eq_derivation.jl:117`) from the type to the struct, drop the `mech = Mechanism(M())` line (118), and change the kernel call (130) from `M` to `mech`. Add a type-accepting delegate. Final shape:

```julia
function _build_wegscheider_rename_map(mech::Mechanism)
    rename = Dict{Symbol, Symbol}()
    step_params = _step_parameters(mech)
    binding_set = Set{Symbol}()
    for (idx, (s, _)) in enumerate(_flat_steps(mech))
        is_equilibrium(s) && is_binding(s) || continue
        push!(binding_set, name(step_params[idx][1], mech))
    end
    dep_raw, _ = _dependent_param_exprs_kernel(mech, rename)
    for (lhs, rhs) in dep_raw
        rhs isa Symbol || continue
        lhs in binding_set && rhs in binding_set || continue
        target = get(rename, rhs, rhs)
        rename[lhs] = target
        for k in collect(keys(rename))
            rename[k] == lhs && (rename[k] = target)
        end
    end
    rename
end

_build_wegscheider_rename_map(M::Type{<:EnzymeMechanism}) =
    _build_wegscheider_rename_map(Mechanism(M()))
_build_wegscheider_rename_map(m::EnzymeMechanism) =
    _build_wegscheider_rename_map(typeof(m))
```

Keep the existing docstring on the function.

- [ ] **Step 4: Run the test, verify it passes.**

Run: `julia --project -e 'using EnzymeRates; include("test/test_types.jl")'`
Expected: PASS.

- [ ] **Step 5: Run the full suite** to confirm the type path still behaves identically: `julia --project -e 'using Pkg; Pkg.test()'`. Expected: green (the type delegate reconstructs the same mechanism it did before).

- [ ] **Step 6: Commit.**

```bash
git add src/rate_eq_derivation.jl test/test_types.jl
git commit -m "Add Mechanism method for _build_wegscheider_rename_map"
```

---

### Task 3: Implement `_merge_tied_kinetic_groups`

**Files:**
- Modify: `src/rate_eq_derivation.jl` (replace the Task-1 stub)
- Test: `test/test_types.jl`

**Interfaces:**
- Produces: `_merge_tied_kinetic_groups(mech::Mechanism)::Vector{Vector{Step}}` — the canonical merged partition. Returns `mech.steps` unchanged (identity, `===`) when nothing merges.
- Consumes: `_build_wegscheider_rename_map(mech)` (Task 2), `_step_parameters(mech)`, `_flat_steps(mech)`, `name(p, mech)`, `is_binding`, `is_equilibrium`.

- [ ] **Step 1: Write the failing test.** This splits a real init-mechanism group and checks the helper re-merges the tied ones:

```julia
@testset "_merge_tied_kinetic_groups re-merges tied splits" begin
    rxn = @enzyme_reaction(begin
        substrates: A[C], B[N]
        products:   P[C], Q[N]
    end)
    merged_something = false
    for m in EnzymeRates.init_mechanisms(rxn)
        gs = EnzymeRates.steps(m)
        for gi in eachindex(gs)
            length(gs[gi]) >= 2 || continue
            # Split group gi into singletons via the provisional (non-merging) ctor.
            split_groups = Vector{Vector{EnzymeRates.Step}}()
            for (j, g) in enumerate(gs)
                j == gi ? append!(split_groups, [[s] for s in g]) :
                          push!(split_groups, copy(g))
            end
            raw = EnzymeRates.Mechanism(EnzymeRates.reaction(m),
                                        split_groups, Val(:_raw))
            # A tied split has the same rate function as the merged original.
            same_rate = EnzymeRates._rate_eq_dedup_key(
                            rate_equation_string(EnzymeRates.compile_mechanism(raw))) ==
                        EnzymeRates._rate_eq_dedup_key(
                            rate_equation_string(EnzymeRates.compile_mechanism(m)))
            same_rate || continue
            merged_something = true
            re = EnzymeRates._merge_tied_kinetic_groups(raw)
            @test length(re) < length(split_groups)   # tied groups collapsed
        end
    end
    @test merged_something   # the test actually exercised a merge
end
```

- [ ] **Step 2: Run it, verify it fails** — with the Task-1 stub, `_merge_tied_kinetic_groups` returns `raw.steps` unchanged, so `length(re) < length(split_groups)` fails.

- [ ] **Step 3: Implement.** Replace the stub at the bottom of `src/rate_eq_derivation.jl` with:

```julia
"""
Canonical kinetic-group partition: merge kinetic groups whose binding-K
representatives are single-symbol Wegscheider-tied (the relation
`_build_wegscheider_rename_map` finds), so split and merged encodings of the
same rate-equivalent graph collapse to one partition. Returns `mech.steps`
unchanged when nothing is tied.
"""
function _merge_tied_kinetic_groups(mech::Mechanism)
    rename = _build_wegscheider_rename_map(mech)
    isempty(rename) && return mech.steps
    groups = mech.steps
    step_params = _step_parameters(mech)
    grp_of_flat = Int[]
    for (gi, g) in enumerate(groups), _ in g
        push!(grp_of_flat, gi)
    end
    # Canonical representative binding-K name per group (nothing if no binding step).
    rep = Vector{Union{Symbol, Nothing}}(nothing, length(groups))
    for (idx, (s, _)) in enumerate(_flat_steps(mech))
        is_equilibrium(s) && is_binding(s) || continue
        k = name(step_params[idx][1], mech)
        rep[grp_of_flat[idx]] = get(rename, k, k)
    end
    byrep = Dict{Symbol, Vector{Int}}()
    for (gi, r) in enumerate(rep)
        r === nothing && continue
        push!(get!(byrep, r, Int[]), gi)
    end
    any(length(v) > 1 for v in values(byrep)) || return mech.steps
    merged = Vector{Vector{Step}}()
    done = falses(length(groups))
    for gi in eachindex(groups)
        done[gi] && continue
        r = rep[gi]
        if r !== nothing && length(byrep[r]) > 1
            gis = byrep[r]
            push!(merged, Step[s for j in gis for s in groups[j]])
            for j in gis
                done[j] = true
            end
        else
            push!(merged, copy(groups[gi]))
            done[gi] = true
        end
    end
    merged
end
```

- [ ] **Step 4: Run the test, verify it passes** (and `merged_something` is true).

Run: `julia --project -e 'using EnzymeRates; include("test/test_types.jl")'`
Expected: PASS. If `merged_something` is false, the bi_bi init set has no mergeable-and-tied group; widen the reaction (add a third substrate) or use the recovered LDH pair (see Task 4 note) — do NOT weaken the assertion.

- [ ] **Step 5: Commit.**

```bash
git add src/rate_eq_derivation.jl test/test_types.jl
git commit -m "Implement _merge_tied_kinetic_groups for canonical partition"
```

---

### Task 4: Wire the merge into the public `Mechanism` constructor (acceptance)

**Files:**
- Verify: `src/types.jl` (the public constructor from Task 1 already calls `_merge_tied_kinetic_groups`)
- Test: `test/test_types.jl`

**Interfaces:**
- Consumes: everything from Tasks 1-3.
- Produces: the invariant that a tied split-form and its merged form build to `==` `Mechanism`s with the same `eq_hash`.

- [ ] **Step 1: Write the failing test** — the acceptance invariant, over real init mechanisms:

```julia
@testset "canonical partition: tied splits build == merged" begin
    rxn = @enzyme_reaction(begin
        substrates: A[C], B[N]
        products:   P[C], Q[N]
    end)
    checked = 0
    for m in EnzymeRates.init_mechanisms(rxn)
        gs = EnzymeRates.steps(m)
        for gi in eachindex(gs)
            length(gs[gi]) >= 2 || continue
            split_groups = Vector{Vector{EnzymeRates.Step}}()
            for (j, g) in enumerate(gs)
                j == gi ? append!(split_groups, [[s] for s in g]) :
                          push!(split_groups, copy(g))
            end
            raw = EnzymeRates.Mechanism(EnzymeRates.reaction(m),
                                        split_groups, Val(:_raw))
            key(x) = EnzymeRates._rate_eq_dedup_key(
                         rate_equation_string(EnzymeRates.compile_mechanism(x)))
            key(raw) == key(m) || continue          # tied split
            checked += 1
            # Public constructor must re-merge the split back to the canonical form.
            rebuilt = EnzymeRates.Mechanism(EnzymeRates.reaction(m), split_groups)
            @test rebuilt == m
            @test key(rebuilt) == key(m)
        end
    end
    @test checked > 0
end
```

- [ ] **Step 2: Run it, verify it fails** — pre-wiring, the public constructor merge collapses correctly only if `_merge_tied_kinetic_groups` is wired; if you completed Task 1's constructor edit, this should already pass. If Task 1 left the public constructor unwired, wire it now (Step 3), then this fails→passes.

- [ ] **Step 3: Confirm the public constructor (from Task 1) calls the merge.** It should read:

```julia
    function Mechanism(reaction::EnzymeReaction,
                       steps::Vector{Vector{Step}})
        prov = Mechanism(reaction, steps, Val(:_raw))
        merged = _merge_tied_kinetic_groups(prov)
        merged === prov.steps ?
            prov : Mechanism(reaction, merged, Val(:_raw))
    end
```

- [ ] **Step 4: Run the test, verify it passes.**

Run: `julia --project -e 'using EnzymeRates; include("test/test_types.jl")'`
Expected: PASS with `checked > 0`.

Note: if you need the exact confirmed LDH pair (`3a2788df` / `f5f7e53b`, a bi-bi with `substrates: NADH, Pyruvate; products: Lactate, NAD`), its two mechanisms are recoverable from the `mechanism_type` column in `docs/ldh_hpc_results/2026_07_03_results/loocv_results.csv` (identical 13-step set, 8 vs 9 kinetic groups). The general bi_bi test above subsumes it; only reconstruct the specific pair if you want a named regression case.

- [ ] **Step 5: Run the full suite.** Some tests will now shift (Task 7 handles them). Record which fail; expect them to be golden/count baselines, not logic failures.

- [ ] **Step 6: Commit.**

```bash
git add src/types.jl test/test_types.jl
git commit -m "Merge Wegscheider-tied kinetic groups in the Mechanism constructor"
```

---

### Task 5: Extend to `AllostericMechanism` (parallel `cat_allo_states`)

**Files:**
- Modify: `src/types.jl:597-643` (the `AllostericMechanism` inner constructor)
- Modify: `src/rate_eq_derivation.jl` (add `_merge_tied_kinetic_groups(am::AllostericMechanism, ...)` using `_state_wegscheider_rename_map`)
- Test: `test/test_types.jl`

**Interfaces:**
- Produces: `AllostericMechanism(reaction, cat_steps, cat_allo_states, mult, sites, ::Val{:_raw})` and a merged public path.
- Consumes: `_state_wegscheider_rename_map(am, state)` (`rate_eq_derivation.jl:1201`), `_state_step_params(am, state)` (`1123`).

**Design note — tag reconciliation:** `cat_steps` and `cat_allo_states` are parallel; merging two catalytic groups must merge their state tags. A binding-K tie only ties groups that play the same catalytic role, so tied groups are expected to share a tag. The merge must assert this and error loudly if two tied groups carry different tags (a case to surface, not silently pick). The per-state binding-K ties come from `_state_wegscheider_rename_map(am, :A)` and `(:I)`; two groups merge only if tied in the state(s) that govern them.

- [ ] **Step 1: Write the failing test** — split a catalytic group of an allosteric init mechanism and assert re-merge:

```julia
@testset "allosteric canonical partition: tied splits build == merged" begin
    # Use an allosteric init mechanism with a multi-binding-step catalytic group.
    rxn = @enzyme_reaction(begin
        substrates: A[C], B[N]
        products:   P[C], Q[N]
    end)
    checked = 0
    for am in EnzymeRates.init_mechanisms(rxn)   # filter to AllostericMechanism if init mixes
        am isa EnzymeRates.AllostericMechanism || continue
        cs = EnzymeRates.steps(am)
        for gi in eachindex(cs)
            length(cs[gi]) >= 2 || continue
            split = Vector{Vector{EnzymeRates.Step}}()
            states = Symbol[]
            for (j, g) in enumerate(cs)
                if j == gi
                    for s in g
                        push!(split, [s]); push!(states, EnzymeRates.cat_allo_state(am, j))
                    end
                else
                    push!(split, copy(g)); push!(states, EnzymeRates.cat_allo_state(am, j))
                end
            end
            raw = EnzymeRates.AllostericMechanism(
                EnzymeRates.reaction(am), split, states,
                EnzymeRates.catalytic_multiplicity(am),
                copy(EnzymeRates.regulatory_sites(am)), Val(:_raw))
            key(x) = EnzymeRates._rate_eq_dedup_key(
                         rate_equation_string(EnzymeRates.compile_mechanism(x)))
            key(raw) == key(am) || continue
            checked += 1
            rebuilt = EnzymeRates.AllostericMechanism(
                EnzymeRates.reaction(am), split, states,
                EnzymeRates.catalytic_multiplicity(am),
                copy(EnzymeRates.regulatory_sites(am)))
            @test rebuilt == am
        end
    end
    # `checked` may be 0 if the chosen reaction yields no allosteric tied split;
    # if so, pick a reaction/spec from MECHANISM_TEST_SPECS that does and adapt.
end
```

- [ ] **Step 2: Run it, verify it fails** (`MethodError` on the `Val(:_raw)` allosteric ctor, or `rebuilt != am`).

- [ ] **Step 3: Implement** the provisional `Val(:_raw)` allosteric inner constructor (mirror Task 1: run `_canon_step_groups` on `cat_steps`, permute `cat_allo_states` with the same perm, validate, `new(...)` — the existing validation block stays), and a public path that calls the allosteric merge. Write `_merge_tied_kinetic_groups(am::AllostericMechanism)` returning `(cat_steps, cat_allo_states)`: compute per-state renames, union groups tied in their governing state, merge steps, and reconcile tags — erroring if two tied groups disagree on tag. Keep the merge structurally identical to the non-allosteric one; the only additions are the two-state rename union and the tag carry-through.

- [ ] **Step 4: Run the test, verify it passes.**

- [ ] **Step 5: Run the full suite;** record shifted goldens for Task 7.

- [ ] **Step 6: Commit.**

```bash
git add src/types.jl src/rate_eq_derivation.jl test/test_types.jl
git commit -m "Canonical kinetic-group partition for AllostericMechanism"
```

---

### Task 6: Assert Pass-2 is now a no-op, then remove it

**Files:**
- Modify: `src/rate_eq_derivation.jl` (delete `_build_wegscheider_rename_map` Pass-2 usage from `_raw_symbolic_rate_polys:391`; likewise `_state_wegscheider_rename_map` where the constructor merge subsumes it)
- Test: `test/test_rate_eq_derivation.jl`

**Interfaces:**
- Consumes: the constructor merge (Tasks 4-5).

- [ ] **Step 1: Write the guard test** — over every curated spec, Pass-2 finds nothing to absorb once mechanisms arrive pre-merged:

```julia
@testset "Pass-2 rename is empty after constructor merge" begin
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism   # adapt to the actual field name in MECHANISM_TEST_SPECS
        @test isempty(EnzymeRates._build_wegscheider_rename_map(
                          EnzymeRates.compile_mechanism(m)))
    end
end
```

Also add, if the flagged families are not already in `MECHANISM_TEST_SPECS`, a case each for the non-competitive-inhibitor and non-essential-activator families (`project_dedup_pass2_dead_code`), asserting the same emptiness.

- [ ] **Step 2: Run it.** If any spec yields a non-empty rename, the constructor merge did **not** cover it — STOP, investigate that mechanism, and do not proceed to deletion. (Expected: all empty.)

- [ ] **Step 3: Remove Pass-2.** With the guard green, replace the `rename_map = _build_wegscheider_rename_map(M)` at `src/rate_eq_derivation.jl:391` with an empty `Dict{Symbol,Symbol}()` (the polynomial rename is now a structural identity), and delete the now-unused `_build_wegscheider_rename_map` Pass-2 body. Do the same for `_state_wegscheider_rename_map` at its call sites (`1166`, `1242`) if the allosteric guard is likewise green. Keep the guard test — it now protects the constructor merge, asserting the derivation never needs to re-absorb.

- [ ] **Step 4: Run the full suite, verify green.**

Run: `julia --project -e 'using Pkg; Pkg.test()'`

- [ ] **Step 5: Commit.**

```bash
git add src/rate_eq_derivation.jl test/test_rate_eq_derivation.jl
git commit -m "Remove derivation-time Pass-2 Wegscheider absorption, subsumed by constructor merge"
```

---

### Task 7: Re-baseline moved tests

**Files:**
- Modify: `test/test_identify_rate_equation.jl:928` (55-class assertion), `test/test_rate_eq_derivation.jl`, `test/test_allosteric_golden.jl`

- [ ] **Step 1: Run the full suite; list every failing assertion.** Expected failures are golden strings / counts for mechanisms whose canonical partition changed, not logic errors.

- [ ] **Step 2: For the 55-class assertion** (`test/test_identify_rate_equation.jl:928`): re-run and confirm the count. Init has no same-rate-function dupes, so it likely stays 55; if a curated init mechanism was itself a tied-split form, the count or its comment may change. Update the number and the explanatory comment to record the constructor merge as the cause (do not just bump the number silently).

- [ ] **Step 3: For each golden** (`rate_equation_string` / Expr-shape / flat-string / `fitted_params`) that changed: verify the change is a merge (fewer kinetic groups / one constant where a tied pair stood), not a regression — the numeric rate law must be unchanged. Regenerate the golden and note in the test comment that the merged partition is the cause.

- [ ] **Step 4: Confirm the `rate_equation` performance test still passes** (`test_rate_equation_performance` in `test/test_rate_eq_derivation.jl`) — `allocs == 0`, `t < 100e-9`. If it fails, STOP: the construction-time change must not touch the runtime body.

- [ ] **Step 5: Run the full suite, verify fully green.**

- [ ] **Step 6: Commit.**

```bash
git add test/
git commit -m "Re-baseline goldens and partition counts for canonical kinetic-group merge"
```

---

### Task 8 (optional): End-to-end LDH cross-check

**Files:** none (verification only)

- [ ] **Step 1:** Re-run the LDH identify pipeline (`docs/ldh_hpc_results/identify_ldh.jl`, reduced `max_param_count` for speed) and re-run the independent fingerprint canonicalizer (`scratchpad/canon_ev.py` / `finger.py`) on the output.

- [ ] **Step 2:** Confirm the distinct-`eq_hash` count falls from ~68k toward the ~26k cross-graph floor, and that the renaming-dup rate drops from ~66% toward the ~7% residual.

- [ ] **Step 3:** Record the before/after numbers in the PR description. No commit.

---

## Self-Review

**Spec coverage:**
- Constructor partition merge (binding-K single-symbol ties) → Tasks 1-4. ✓
- Allosteric constructor → Task 5. ✓
- Provisional non-merging path to avoid recursion → Task 1 (`Val(:_raw)`). ✓
- `_build_wegscheider_rename_map` refactor to take a struct → Task 2. ✓
- Remove Pass-2 after asserting no-op → Task 6. ✓
- Blast radius (55-class, goldens, perf contract) → Task 7. ✓
- Validation (tied split builds == merged, same `eq_hash`) → Task 4 (+ Task 5 allosteric). ✓
- Defer 15.2% cross-graph residual → out of scope; no task, correct. ✓

**Correction to the spec's blast-radius claim:** the spec says the 55-class assertion "should fall." Init enumeration has 0% same-rate-function dupes (measured), so the count likely holds at 55; the fix's effect is on expanded/split forms, which the init-only test does not exercise. Task 7 Step 2 verifies rather than assumes.

**Placeholder scan:** no TBD/TODO; every code step shows code. The one genuine unknown — the exact `MECHANISM_TEST_SPECS` field name — is flagged inline in Task 6 Step 1 ("adapt to the actual field name").

**Type consistency:** `_merge_tied_kinetic_groups(mech)::Vector{Vector{Step}}` returns `mech.steps` (identity) on no-op; the public constructor checks `merged === prov.steps` — consistent. `Val(:_raw)` used identically in Tasks 1, 3, 4, 5.
