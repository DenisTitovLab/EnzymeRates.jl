# Concrete-types refactor — continuation implementation plan (Stages 6β–7e)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete the concrete-types refactor on branch `refactor-to-concrete-types-instead-of-symbols`: retire spec types, switch enumeration dispatch to the concrete `EnzymeReaction`, switch DSL to emit decomposed Species, delete the legacy parallel data path, and land in a single PR at ≤8,200 src LOC with all existing test coverage preserved.

**Architecture:** One concrete struct family (`Mechanism`/`Step`/`Species`/`Parameter`) shared between enumeration and derivation. The chokepoint `name(p::Parameter, m::Mechanism)` is the sole `Parameter → Symbol` mapping. The DSL emits `EnzymeMechanism(Mechanism(...))` directly with decomposed-Species notation (`E(S)`). Enumeration dispatches on the concrete `EnzymeReaction`; `EnzymeReactionLegacy` and the opaque-form helpers are deleted.

**Tech Stack:** Julia 1.10+, Test.jl, Aqua.jl, JET.jl. The codebase uses `@generated` functions for type-stable rate-equation derivation; perf-critical paths must stay 0-alloc/<100ns.

**Source spec:** `docs/superpowers/specs/2026-05-22-concrete-types-refactor-continuation-design.md` (commit `7f84916`).

---

## Conventions for every task

These apply to every commit in this plan:

- **TDD (per `.claude/CLAUDE.md`):** write the failing test, run-fail, implement, run-pass, refactor if needed.
- **Test integrity (spec §2 + continuation §4):** never delete, comment, weaken, or `@test_skip`/`@test_broken` a test. When migrating a testset, the replacement must land in the same commit as the deletion. If you can't fit migration + deletion in one commit, split: migration first, deletion second.
- **Run the test-integrity script after every commit:**
  ```bash
  bash scripts/check_test_integrity.sh main
  ```
  Expected: PASS. If FAIL: stop, restore deleted/weakened tests, recommit.
- **Run the full suite after every commit:**
  ```bash
  julia --project=test -e 'include("test/runtests.jl")'
  ```
  (This dev box OOMs on `Pkg.test()`; use the include path.) Expected: all tests pass. Test count fluctuates ±1 between runs due to CMA-ES non-determinism in `test_identify_rate_equation.jl` — don't over-interpret single-test count shifts.
- **Commit message format:**
  ```
  <stage>: <one-line summary>

  <optional body>

  src delta: -X / +Y net Z, cumulative: ±W
  ```
  Compute Δ via `wc -l src/*.jl` before vs after. Cumulative is vs main baseline (7,136).
- **No `--amend`.** Always create new commits.
- **No temporal-context comments** in code (no "Stage N", "previously", "legacy", "will be"). Documentation in plans/specs is fine; code comments must be evergreen.
- **Perf gates green at every commit:** `test_rate_equation_performance` in `test/test_rate_eq_derivation.jl` must show 0 allocs and < 100ns per call. The 3 compile-budget gates in `test/test_compile_budget.jl` must pass.
- **Branch:** stay on `refactor-to-concrete-types-instead-of-symbols`. Single PR at the end.

---

## File structure overview

**Files modified across the plan:**

| File | Stages touching it | Role |
|------|-------------------|------|
| `test/test_mechanism_enumeration.jl` | 6β | 50+ spec testsets migrated to Mechanism construction |
| `src/mechanism_enumeration.jl` | 6β, 7b, 7c, 7d | spec types + Legacy dispatch + opaque helpers deleted; `_expand_*` rewritten to use structured Species |
| `src/identify_rate_equation.jl` | 6, 7c | canonicalizer rewrite; dispatch on `EnzymeReaction` |
| `src/types.jl` | 7a, 7d | chokepoint usage; `EnzymeReactionLegacy` + dual-Sig branches + 2-arg `EnzymeMechanism` deleted |
| `src/rate_eq_derivation.jl` | 7a | 4 direct `Symbol("K$idx")` sites routed through chokepoint |
| `src/dsl.jl` | 7b | grammar extended to decomposed-Species; macros switch to emit `Mechanism(...)` |
| `test/mechanism_definitions_for_test_enzyme_derivation.jl` | 7b | 37 DSL invocations migrated to new grammar |
| `test/test_dsl.jl`, `test/test_types.jl`, `test/test_rate_eq_derivation.jl`, others | 7b | 70+ DSL invocations migrated |
| `test/test_compile_budget.jl` | 7c | direct `EnzymeReactionLegacy` construction switched to `EnzymeReaction`; budgets re-baselined if compile metrics shift |
| `README.md`, `.claude/CLAUDE.md` | 7e | docs updated to final architecture |

**Files created:**

| File | Stage | Purpose |
|------|-------|---------|
| `test/test_chokepoint.jl` | 7a | regression test asserting only `name(p, m)` renders parameter symbols |
| `docs/superpowers/refactor-deleted-tests.md` | as needed | log entries for any §2.1-exception test deletions |

---

# Stage 6β — Migrate spec testsets to Mechanism (opaque-form)

**Stage goal:** retire `MechanismSpec`/`StepSpec`/`AllostericMechanismSpec` from `test/test_mechanism_enumeration.jl` and from `src/mechanism_enumeration.jl` without losing assertion coverage. Use opaque-form Species (`Species([], :ES)`) so the existing `_expand_*` moves work unchanged.

**Stage LOC delta target:** −500 to −800 in `src/mechanism_enumeration.jl`.

**Migration pattern (refer back from each task):**

Translate this spec construction:
```julia
spec = MechanismSpec(rxn,
    [StepSpec([:E, :S], [:E_S], true, 1),
     StepSpec([:E, :P], [:E_P], true, 2),
     StepSpec([:E_S], [:E_P], false, 3)],
    3)
variants = EnzymeRates._expand_re_to_ss(spec)
```

To this Mechanism construction (opaque-form Species, structural grouping by Vector{Vector{Step}}):
```julia
m = Mechanism(rxn,
    [[Step(Species([], :E), Species([], :E_S), Substrate(:S), true)],
     [Step(Species([], :E), Species([], :E_P), Product(:P), true)],
     [Step(Species([], :E_S), Species([], :E_P), nothing, false)]])
variants = EnzymeRates._expand_re_to_ss(m)
```

Key translation rules:
- `MechanismSpec(rxn, [step1, step2, ...], n)` → `Mechanism(rxn, group_steps([step1, step2, ...]))` where `group_steps` puts each StepSpec into its own inner Vector unless multiple StepSpecs share `kinetic_group::Int` (then they're grouped).
- `StepSpec([:E, :S], [:E_S], true, kg)` → `Step(Species([], :E), Species([], :E_S), Substrate(:S), true)`. The bound metabolite is detected from `reactants ∪ products` difference; for binding (`E + S → ES`), it's `:S`. For dissociation (`ES → E + P`), the order is reversed and bound is `:P`.
- `StepSpec([:E_S], [:E_P], false, kg)` (iso, no metabolite) → `Step(Species([], :E_S), Species([], :E_P), nothing, false)`.
- For substrate vs product vs regulator classification, look at `rxn`'s declared substrates/products/regulators and dispatch the right metabolite constructor: `Substrate(:S)`, `Product(:P)`, `CompetitiveInhibitor(:I)`, `AllostericRegulator(:R)`.
- `spec.kinetic_group[i]` indexing → reach via structural lookup: `findfirst(g -> step in g, m.steps)`.
- `spec.n_fit_params_estimate` → `EnzymeRates._n_fit_params_estimate(m)`.
- `EnzymeRates._init_mechanism_specs(rxn)` → `EnzymeRates.init_mechanisms(rxn)`.
- `_assert_spec_invariants(spec)` → write a parallel `_assert_mechanism_invariants(m)` (Task 6β.1 creates this).
- AllostericMechanismSpec → AllostericMechanism (similar pattern; tags become Dict fields).

If a translation is ambiguous (you can't tell which Species role to use, or which group a step belongs to), STOP and ask — don't guess.

## Task 6β.1: Migrate test helpers + add `_assert_mechanism_invariants`

**Files:**
- Modify: `test/test_mechanism_enumeration.jl:1-500` (helpers section)

The current file has spec-form helpers near the top (`mechanism_spec_from_mechanism_and_rxn`, `allosteric_spec_from_mechanism_and_rxn`, `_assert_spec_invariants`, `enumerate_all`). This task adds parallel Mechanism-form helpers used by every subsequent 6β task. Spec helpers stay until Task 6β.10 (deletion).

- [ ] **Step 1: Read the current helper block**

Read `test/test_mechanism_enumeration.jl` lines 1–500 to understand existing helpers.

- [ ] **Step 2: Write a failing self-test for the new helper**

Add this testset at the end of the helpers block (around line 490, before `@testset "Mechanism Enumeration"`):

```julia
@testset "_assert_mechanism_invariants on uni-uni init" begin
    rxn = @enzyme_reaction begin
        substrates: S[C6H12O6]
        products:   P[C6H13O9P]
    end
    m = first(EnzymeRates.init_mechanisms(rxn))
    @test EnzymeRates._assert_mechanism_invariants(m) === nothing
end
```

- [ ] **Step 3: Run the failing test**

```bash
julia --project=test -e 'include("test/test_mechanism_enumeration.jl")' 2>&1 | tail -20
```

Expected: FAIL with "UndefVarError: `_assert_mechanism_invariants` not defined" (or similar).

- [ ] **Step 4: Implement `_assert_mechanism_invariants` in `src/mechanism_enumeration.jl`**

Add near the bottom of the file (after the existing `_assert_spec_invariants` helpers):

```julia
"""
    _assert_mechanism_invariants(m::Mechanism) -> Nothing

Structural invariants every valid Mechanism should satisfy:
- Every group is non-empty
- source_idx values are unique and dense (1 through n_steps)
- Each binding step's bound_metabolite is non-nothing AND iso steps have nothing
- from_species != to_species for every step
"""
function _assert_mechanism_invariants(m::Mechanism)
    flat = collect(Iterators.flatten(m.steps))
    isempty(flat) && error("empty steps in Mechanism")
    for g in m.steps
        isempty(g) && error("empty kinetic group in Mechanism")
    end
    src_indices = [s.source_idx for s in flat]
    sort(src_indices) == collect(1:length(flat)) ||
        error("source_idx values not dense 1..n: got $(sort(src_indices))")
    for s in flat
        if is_binding(s)
            s.bound_metabolite === nothing &&
                error("binding step has nothing bound_metabolite")
        else
            s.bound_metabolite === nothing ||
                error("iso step has non-nothing bound_metabolite")
        end
        s.from_species == s.to_species &&
            error("from_species == to_species in step $s")
    end
    nothing
end

function _assert_mechanism_invariants(m::AllostericMechanism)
    _assert_mechanism_invariants(m.base)
    # Dense tag storage: every kinetic group has a group_tag, every ligand has a reg_ligand_tag
    n_groups = length(m.base.steps)
    Set(keys(m.group_tags)) == Set(1:n_groups) ||
        error("group_tags not dense: keys=$(sort(collect(keys(m.group_tags))))")
    declared_ligands = Set(name(metabolite(r)) for s in m.allosteric_reg_sites for r in s)
    Set(keys(m.reg_ligand_tags)) == declared_ligands ||
        error("reg_ligand_tags not dense")
    nothing
end
```

(If `_assert_mechanism_invariants` already exists from prior commits, skim its body and skip to Step 5.)

- [ ] **Step 5: Run the test to verify it passes**

```bash
julia --project=test -e 'include("test/test_mechanism_enumeration.jl")' 2>&1 | tail -20
```

Expected: PASS.

- [ ] **Step 6: Add `enumerate_all_mechanism` helper to the test file**

Add this above the existing `enumerate_all` helper in `test/test_mechanism_enumeration.jl`:

```julia
"""
    enumerate_all_mechanism(rxn::EnzymeReaction) -> Dict{Int, Vector{<:AbstractMechanism}}

Mechanism-form parallel of `enumerate_all` (which works on specs).
Used by 6β-migrated testsets that exercise full enumeration via the
Mechanism API.
"""
function enumerate_all_mechanism(rxn)
    cache = Dict{Int, Vector{EnzymeRates.AbstractMechanism}}()
    for m in EnzymeRates.init_mechanisms(rxn)
        pc = EnzymeRates._n_fit_params_estimate(m)
        push!(get!(cache, pc, EnzymeRates.AbstractMechanism[]), m)
    end
    results = Dict{Int, Vector{EnzymeRates.AbstractMechanism}}()
    while !isempty(cache)
        pc = minimum(keys(cache))
        level = pop!(cache, pc, EnzymeRates.AbstractMechanism[])
        EnzymeRates.dedup!(level)
        results[pc] = level
        for m in level
            for child in EnzymeRates.expand_mechanisms_one(m)
                cpc = EnzymeRates._n_fit_params_estimate(child)
                cpc > pc || continue
                push!(get!(cache, cpc, EnzymeRates.AbstractMechanism[]), child)
            end
        end
    end
    results
end
```

(If `expand_mechanisms_one` doesn't exist yet as a single-mechanism entry point, locate the existing per-move dispatch in `src/mechanism_enumeration.jl` and either add a thin wrapper or call the moves directly inside this helper.)

- [ ] **Step 7: Verify the helper compiles**

```bash
julia --project=test -e 'include("test/test_mechanism_enumeration.jl"); println("compiled")' 2>&1 | tail -5
```

Expected: "compiled" with no errors.

- [ ] **Step 8: Run the full test-integrity check + suite**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -10
```

Expected: integrity PASS; test suite passes (26,904 ± a few).

- [ ] **Step 9: Commit**

```bash
git add test/test_mechanism_enumeration.jl src/mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Stage 6β.1: add Mechanism-form test helpers parallel to spec helpers

_assert_mechanism_invariants validates structural shape of Mechanism /
AllostericMechanism. enumerate_all_mechanism parallels enumerate_all
for full-enumeration integration tests. Both are used by subsequent
6β tasks that migrate spec testsets to Mechanism form.

src delta: -X / +Y net Z, cumulative: ±W
EOF
)"
```

Fill in actual LOC deltas after `wc -l src/*.jl`.

## Task 6β.2: Migrate `_expand_re_to_ss` testsets

**Files:**
- Modify: `test/test_mechanism_enumeration.jl:1601-2270` (the `@testset "_expand_re_to_ss"` block)

The block contains both `MechanismSpec — ...` testsets (lines ~1603–2180) and a few earlier-migrated `Mechanism — ...` testsets (lines ~2182–2270). The spec ones need parallel Mechanism testsets added; spec testsets stay until Task 6β.10.

- [ ] **Step 1: Inventory the spec testsets in the block**

```bash
sed -n '1601,2275p' test/test_mechanism_enumeration.jl | grep '@testset' 
```

Expected output: ~10 testset headers. Note each `MechanismSpec — ...` line.

- [ ] **Step 2: Pick the first spec testset to migrate — write the Mechanism-form parallel**

For `@testset "MechanismSpec — uni-uni: 2 RE binding groups → 2 variants"` (starts ~line 1603), write a parallel Mechanism testset *immediately after* it. Translation follows the pattern at the top of Stage 6β.

Worked example:

```julia
@testset "Mechanism — uni-uni: 2 RE binding groups → 2 variants" begin
    rxn = @enzyme_reaction begin
        substrates: S[C6H12O6]
        products:   P[C6H13O9P]
    end
    m = Mechanism(rxn,
        [[Step(Species([], :E), Species([], :E_S), Substrate(:S), true)],
         [Step(Species([], :E), Species([], :E_P), Product(:P), true)],
         [Step(Species([], :E_S), Species([], :E_P), nothing, false)]])
    EnzymeRates._assert_mechanism_invariants(m)

    variants = EnzymeRates._expand_re_to_ss(m)
    @test length(variants) == 2

    # Each variant converts one of the two RE binding groups to SS.
    # In the spec version this was asserted by counting RE vs SS step counts;
    # mirror that assertion structurally:
    for v in variants
        EnzymeRates._assert_mechanism_invariants(v)
        n_re = count(s -> s.is_equilibrium, Iterators.flatten(v.steps))
        n_ss = count(s -> !s.is_equilibrium, Iterators.flatten(v.steps))
        @test n_re + n_ss == 3
        @test n_ss == 2  # 1 original SS iso + 1 RE→SS-converted
    end
end
```

- [ ] **Step 3: Run the new testset to verify it passes alongside the existing spec one**

```bash
julia --project=test -e 'include("test/test_mechanism_enumeration.jl")' 2>&1 | grep -E "uni-uni: 2 RE|Test Failed|PASS|FAIL" | head -20
```

Expected: both spec and Mechanism testsets pass.

- [ ] **Step 4: Migrate the remaining spec testsets in the block**

Apply the same pattern to each remaining `MechanismSpec — ...` testset in the block (lines ~1677–2180). For each:
1. Read the spec testset's setup + assertions.
2. Write the parallel `Mechanism — ...` testset immediately after it.
3. Verify both pass.

The full list (from Step 1):
- `MechanismSpec — bi-bi sequential: 4 RE binding groups → 4 variants`
- `MechanismSpec — bi-bi multi-step kinetic group: atomic conversion`
- `MechanismSpec — bi-bi ping-pong: 5 RE groups → 5 variants`
- `MechanismSpec — ter-ter sequential`
- `MechanismSpec — all-SS catalytic seed: empty (negative)`
- `AllostericMechanismSpec — :EqualRT group: Δ=+1`
- `AllostericMechanismSpec — :OnlyR group: Δ=+1`
- `AllostericMechanismSpec — :NonequalRT group: Δ=+2`
- `Substrate-as-dead-end-inhibitor overlap (S used as both)`
- `AllostericMechanismSpec — substrate-as-dead-end-I overlap`

If you discover an assertion in a spec testset that doesn't translate cleanly to the Mechanism API (e.g., it inspects `spec.kinetic_group` field directly in a way no public Mechanism accessor exposes), STOP and ask — don't quietly weaken the assertion.

- [ ] **Step 5: Run the full suite + integrity check**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -10
```

Expected: integrity PASS; test count grew by ~10 (one per migrated testset).

- [ ] **Step 6: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Stage 6β.2: migrate _expand_re_to_ss spec testsets to Mechanism form

Adds Mechanism-form parallels of each MechanismSpec/AllostericMechanismSpec
testset in the _expand_re_to_ss block. Spec testsets retained until
6β.10 deletes the spec types entirely; the parallel ensures we have a
clean migration target with full assertion coverage.

src delta: 0 / 0 (test-only)
EOF
)"
```

## Task 6β.3: Migrate `_expand_split_kinetic_group` testsets

**Files:**
- Modify: `test/test_mechanism_enumeration.jl:2273-2670` (the `@testset "_expand_split_kinetic_group"` block)

Same pattern as Task 6β.2.

- [ ] **Step 1: Inventory spec testsets**

```bash
sed -n '2273,2670p' test/test_mechanism_enumeration.jl | grep '@testset'
```

- [ ] **Step 2: For each `MechanismSpec — ...` / `AllostericMechanismSpec — ...` testset, write the Mechanism-form parallel immediately after it**

Inventory (from prior survey):
- `MechanismSpec — bi-bi: 4 multi-step groups → 8 splits`
- `MechanismSpec — all singleton groups: empty (negative)`
- `MechanismSpec — mixed RE/SS group sizes: deltas differ`
- `AllostericMechanismSpec — split inherits parent's tag`
- `AllostericMechanismSpec — SS multi-step :NonequalRT split: Δ=+4`

Translation pattern: same as 6β.2, plus the split-specific assertion that the variant count = sum of (group_size - 1) over each multi-step group. The Mechanism version reads group size as `length(m.steps[i])` instead of counting `spec.kinetic_group` entries.

Worked example (one testset):

```julia
@testset "Mechanism — bi-bi: 4 multi-step groups → 8 splits" begin
    # Setup mirrors the MechanismSpec testset above, but with grouped steps:
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products:   P[C], Q[N]
    end
    m = Mechanism(rxn,
        [[Step(Species([], :E), Species([], :E_A), Substrate(:A), true),
          Step(Species([], :E_B), Species([], :E_AB), Substrate(:A), true)],
         [Step(Species([], :E), Species([], :E_B), Substrate(:B), true),
          Step(Species([], :E_A), Species([], :E_AB), Substrate(:B), true)],
         [Step(Species([], :E_AB), Species([], :E_PQ), nothing, false)],
         [Step(Species([], :E_PQ), Species([], :E_Q), Product(:P), true),
          Step(Species([], :E_P), Species([], :E), Product(:P), true)],
         [Step(Species([], :E_PQ), Species([], :E_P), Product(:Q), true),
          Step(Species([], :E_Q), Species([], :E), Product(:Q), true)]])
    EnzymeRates._assert_mechanism_invariants(m)

    variants = EnzymeRates._expand_split_kinetic_group(m)
    @test length(variants) == 4  # 4 multi-step groups, each can be split once
                                  # — assertion mirrors the original spec testset
    for v in variants
        EnzymeRates._assert_mechanism_invariants(v)
        @test length(v.steps) == length(m.steps) + 1  # one group became two
    end
end
```

(Verify the assertion counts against the existing spec testset — translate, don't invent.)

- [ ] **Step 3: Run the suite + integrity**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -10
```

- [ ] **Step 4: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Stage 6β.3: migrate _expand_split_kinetic_group testsets to Mechanism form

src delta: 0 / 0 (test-only)
EOF
)"
```

## Task 6β.4: Migrate `_expand_add_dead_end_regulator` testsets

**Files:**
- Modify: `test/test_mechanism_enumeration.jl:2675-3385` (the `@testset "_expand_add_dead_end_regulator"` block)

Pattern as before. Inventory ~12 spec testsets.

- [ ] **Step 1: Inventory + migrate per the 6β.2 pattern**

```bash
sed -n '2675,3385p' test/test_mechanism_enumeration.jl | grep '@testset'
```

Spec testsets to migrate (from prior survey):
- `Uni-uni + I: 1 variant (equivalence-style)` — not a `MechanismSpec —` prefix; check if this is pre-migrated or still spec-construction. If `MechanismSpec` appears in the setup, migrate; if not, skip.
- `Uni-uni no regulators → empty (negative)`
- `exclude_regs kwarg suppresses regulator addition (negative)`
- `Sequential bi-bi + I: 4 distinct form sets`
- `Bi-bi random + I: 9 variants (property-style)`
- `Bi-bi PP + I: 3 variants`
- `Two regulators chain: J__reg dummy naming has no numeric suffix`
- `Two regulators competition: 17 variants (property-style)`
- `Substrate-as-dead-end-inhibitor overlap (S used as both)`
- `AllostericMechanismSpec input: dead-end binding tagged :EqualRT`
- `Allosteric-only regulator → empty (negative)`
- `Plain MechanismSpec: I__reg bindings share one kinetic group`

For each, follow the 6β.2 pattern: read setup, write Mechanism parallel, verify both pass.

- [ ] **Step 2: Run suite + integrity + commit**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -10
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Stage 6β.4: migrate _expand_add_dead_end_regulator testsets to Mechanism form

src delta: 0 / 0 (test-only)
EOF
)"
```

## Task 6β.5: Migrate `_expand_to_allosteric` testsets

**Files:**
- Modify: `test/test_mechanism_enumeration.jl:3386-3640` (the `@testset "_expand_to_allosteric"` block)

- [ ] **Step 1: Inventory + migrate per the 6β.2 pattern**

```bash
sed -n '3386,3640p' test/test_mechanism_enumeration.jl | grep '@testset'
```

Spec testsets to migrate:
- `MechanismSpec — uni-uni: 4 variants (equivalence-style)`
- `AllostericMechanismSpec → empty (negative)`
- `oligomeric_state from reaction`
- `Bi-bi sequential: 5 groups → 6 variants`
- `Bi-bi ping-pong: 6 groups → 7 variants`

For each, follow the 6β.2 pattern.

- [ ] **Step 2: Run suite + integrity + commit**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -10
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Stage 6β.5: migrate _expand_to_allosteric testsets to Mechanism form

src delta: 0 / 0 (test-only)
EOF
)"
```

## Task 6β.6: Migrate `_expand_add_allosteric_regulator` testsets

**Files:**
- Modify: `test/test_mechanism_enumeration.jl:3640-4115` (the `@testset "_expand_add_allosteric_regulator"` block)

- [ ] **Step 1: Inventory + migrate**

```bash
sed -n '3640,4118p' test/test_mechanism_enumeration.jl | grep '@testset'
```

Spec testsets to migrate (inventory from prior survey):
- `Allosteric uni-uni + first allo regulator R: 3 variants`
- `existing_de exclusion prevents adding bound dead-end as allo regulator`
- `Non-allosteric MechanismSpec → empty (negative)`
- `Two regulators with site options: count = 7`
- `EqualRT ligand reachable at existing reg site`
- `Adding :EqualRT R2 at site with :OnlyT R1`
- `Substrate-as-allosteric-regulator overlap`
- `Two regulators at different sites`
- `Product-as-allosteric-regulator overlap`
- `All declared regs already present → empty (negative)`

(Some may already be in pre-migrated `AllostericMechanism — ...` form. Skip those.)

- [ ] **Step 2: Run + commit**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -10
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Stage 6β.6: migrate _expand_add_allosteric_regulator testsets to Mechanism form

src delta: 0 / 0 (test-only)
EOF
)"
```

## Task 6β.7: Migrate `_expand_change_allo_state` testsets

**Files:**
- Modify: `test/test_mechanism_enumeration.jl:4118-4440` (the `@testset "_expand_change_allo_state"` block)

- [ ] **Step 1: Inventory + migrate**

```bash
sed -n '4118,4440p' test/test_mechanism_enumeration.jl | grep '@testset'
```

Spec testsets to migrate:
- `Allosteric uni-uni all-:EqualRT: 3 group-tag relaxations`
- `Fully relaxed → empty (negative)`
- `MechanismSpec → empty`
- `Allosteric regulator tag removal delta`
- `:OnlyT regulator-ligand relaxation`
- `Multiple regulator ligands at independent tags`
- `Non-default site multiplicity (catalytic_n=4, reg site multiplicity=2)`

- [ ] **Step 2: Run + commit**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -10
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Stage 6β.7: migrate _expand_change_allo_state testsets to Mechanism form

src delta: 0 / 0 (test-only)
EOF
)"
```

## Task 6β.8: Migrate `_canonicalize!`/`_dedup_key`/Dedup testsets

**Files:**
- Modify: `test/test_mechanism_enumeration.jl:4439-4742` (the `@testset "_canonicalize!"`, `@testset "_dedup_key"`, `@testset "Dedup"` blocks)

The current `Dedup` block has BOTH spec and Mechanism testsets side by side; this task only adds Mechanism parallels to any spec-only testsets in the `_canonicalize!` + `_dedup_key` blocks.

- [ ] **Step 1: Inventory**

```bash
sed -n '4439,4742p' test/test_mechanism_enumeration.jl | grep '@testset'
```

- [ ] **Step 2: For each spec testset in the `_canonicalize!`/`_dedup_key` blocks, write a Mechanism-form parallel using the public dedup API**

The Mechanism replacement for `_canonicalize!(spec)` is `EnzymeRates._canonicalize_mechanism!(m)` (verify the function exists in `src/mechanism_enumeration.jl`; if not, the canonicalization is internal to `dedup!`). The Mechanism replacement for `_dedup_key(spec)` is whatever internal hash `dedup!` uses; for tests, you can usually substitute `(s.from_species, s.to_species, s.bound_metabolite, s.is_equilibrium)`-based tuple keys.

- [ ] **Step 3: Run + commit**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -10
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Stage 6β.8: migrate _canonicalize!/_dedup_key testsets to Mechanism form

src delta: 0 / 0 (test-only)
EOF
)"
```

## Task 6β.9: Migrate `expand_mechanisms` + `Integration` testsets

**Files:**
- Modify: `test/test_mechanism_enumeration.jl:4743-5100` (the `@testset "expand_mechanisms"` and `@testset "Integration"` blocks)

- [ ] **Step 1: Inventory**

```bash
sed -n '4743,5100p' test/test_mechanism_enumeration.jl | grep '@testset'
```

- [ ] **Step 2: For each spec testset, write a Mechanism-form parallel**

Use `enumerate_all_mechanism(rxn)` (added in 6β.1) for full enumeration tests. Replace `enumerate_all(rxn)` calls with the Mechanism variant in the new parallel testsets.

Worked example for `Uni-uni full enumeration`:

```julia
@testset "Mechanism — Uni-uni full enumeration" begin
    rxn = @enzyme_reaction begin
        substrates: S[C6H12O6]
        products:   P[C6H13O9P]
    end
    results = enumerate_all_mechanism(rxn)
    @test !isempty(results)
    min_pc = minimum(keys(results))
    @test length(results[min_pc]) == 1  # mirrors the spec testset's assertion
    for (pc, level) in results
        for m in level
            EnzymeRates._assert_mechanism_invariants(m)
        end
    end
end
```

- [ ] **Step 3: Run + commit**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -10
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Stage 6β.9: migrate expand_mechanisms + Integration testsets to Mechanism form

src delta: 0 / 0 (test-only)
EOF
)"
```

## Task 6β.10: Delete spec types + their overloads

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (delete spec testsets, helpers; possibly delete entire `@testset "AllostericMechanismSpec constructor density validation"` and `@testset "spec-from-mechanism helpers reject inconsistent inputs"` and `@testset "allosteric_spec_from_mechanism_and_rxn round-trip"`)
- Modify: `src/mechanism_enumeration.jl` (delete spec types + overloads)
- Modify: `test/runtests.jl` (drop `StepSpec, MechanismSpec, AllostericMechanismSpec, AbstractMechanismSpec` from the `using EnzymeRates: ...` imports if present)

This is the deletion commit. Every spec testset deleted here MUST have a corresponding Mechanism-form testset in the file from Tasks 6β.2–6β.9. The continuation rule (§4 of spec) requires migration + deletion to either be in the same commit OR migration-first/deletion-second commits.

Since Tasks 6β.2–6β.9 already migrated the testsets across 8 commits, the rule is satisfied: this commit is the "deletion-second" half.

- [ ] **Step 1: Verify migration completeness BEFORE deleting**

Grep for `MechanismSpec` / `StepSpec` / `AllostericMechanismSpec` usage outside the helpers block:

```bash
grep -n "MechanismSpec\|StepSpec\|AllostericMechanismSpec" test/test_mechanism_enumeration.jl | grep -v "^[0-9]*:#" | head -50
```

For each spec testset, search the file for a corresponding `Mechanism — ` or `AllostericMechanism — ` parallel testset. If any spec testset has no parallel, STOP — go back and complete the migration before deleting.

- [ ] **Step 2: Delete the spec testsets in `test/test_mechanism_enumeration.jl`**

Strategy: bottom-up. Start from the `@testset "Integration"` block, scroll up, and delete each `MechanismSpec — ...` / `AllostericMechanismSpec — ...` testset whose Mechanism parallel exists.

Also delete:
- `@testset "AllostericMechanismSpec constructor density validation"` (lines ~147–187) — its Mechanism parallel is the existing `AllostericMechanism` constructor validation in `test/test_types.jl`. Verify before deleting; if no parallel exists, write one in `test/test_types.jl` BEFORE deleting here.
- `@testset "spec-from-mechanism helpers reject inconsistent inputs"` (lines ~190–256) — these test `mechanism_spec_from_mechanism_and_rxn`, which is the spec-construction helper being deleted. After 6β.10, the function is gone, so the testset has no surviving referent. This is a §2.1 narrow-exception deletion — add a log entry (Step 4 below).
- `@testset "allosteric_spec_from_mechanism_and_rxn round-trip"` (lines ~257–...) — same as above. Log entry needed.

Delete the spec-form helpers:
- `mechanism_spec_from_mechanism_and_rxn` function
- `allosteric_spec_from_mechanism_and_rxn` function
- `_assert_spec_invariants` methods (both for `MechanismSpec` and `AllostericMechanismSpec`)
- The original `enumerate_all` function (the Mechanism version `enumerate_all_mechanism` from 6β.1 supersedes it).

- [ ] **Step 3: Delete the spec types + overloads in `src/mechanism_enumeration.jl`**

Bottom-up:
1. `StepSpec`, `MechanismSpec`, `AllostericMechanismSpec`, `AbstractMechanismSpec` struct definitions.
2. `_expand_re_to_ss(::AbstractMechanismSpec)`, `_expand_split_kinetic_group(::AbstractMechanismSpec)`, `_expand_add_dead_end_regulator(::AbstractMechanismSpec)`, `_expand_to_allosteric(::AbstractMechanismSpec)`, `_expand_add_allosteric_regulator(::AbstractMechanismSpec)`, `_expand_change_allo_state(::AbstractMechanismSpec)` — six legacy overloads.
3. `dedup!(::Dict{Int, Vector{<:AbstractMechanismSpec}})` legacy overload.
4. `_canonicalize!(::AbstractMechanismSpec)`.
5. `_dedup_key(::AbstractMechanismSpec)`.
6. `_init_mechanism_specs(::EnzymeReaction)`, `_init_mechanism_specs(::EnzymeReactionLegacy)`. (The 2-arg method on `EnzymeReactionLegacy` is referenced from `src/identify_rate_equation.jl`; if so, update the caller to use `init_mechanisms` directly.)
7. `_spec_from_mechanism`, `_mechanism_from_spec`.
8. `EnzymeMechanism(spec::MechanismSpec)` adapter.
9. `expand_mechanisms(::Vector{<:AbstractMechanismSpec}, ::EnzymeReactionLegacy)` legacy overload (the Mechanism version stays).
10. `_n_fit_params_estimate(::AbstractMechanismSpec)` if defined as a separate method.

For each deletion, also delete the docstring above it.

- [ ] **Step 4: Add §2.1 log entries for the spec-helper testset deletions**

```bash
mkdir -p docs/superpowers
cat > docs/superpowers/refactor-deleted-tests.md <<'EOF'
# Refactor — Deleted test log

Per spec §2.1 narrow exception: tests of deleted helpers whose
behavior is preserved by the replacement code path.

## Stage 6β.10 (commit TBD-after-commit)

### test_mechanism_enumeration.jl `@testset "spec-from-mechanism helpers reject inconsistent inputs"`
- Lines: ~190–256 (pre-deletion)
- Helper deleted: `mechanism_spec_from_mechanism_and_rxn`
- Replacement: direct `Mechanism(rxn, [[Step(...)]])` construction in
  6β-migrated testsets; the spec round-trip the helper exercised is
  no longer needed because callers don't go through specs.
- Integration coverage:
  - `test_mechanism_enumeration.jl @testset "_expand_re_to_ss"` (and
    each other `_expand_*` block) — the Mechanism inputs exercise the
    same construction paths the helper used to validate.

### test_mechanism_enumeration.jl `@testset "allosteric_spec_from_mechanism_and_rxn round-trip"`
- Lines: ~257–... (pre-deletion)
- Helper deleted: `allosteric_spec_from_mechanism_and_rxn`
- Replacement: direct `AllostericMechanism(...)` construction in 6β-migrated
  testsets.
- Integration coverage: same as above, plus
  `test_types.jl @testset "AllostericMechanism construction"`.
EOF
```

After the commit, edit `docs/superpowers/refactor-deleted-tests.md` to replace `TBD-after-commit` with the actual commit SHA.

- [ ] **Step 5: Run full suite + integrity check**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -30
wc -l src/*.jl
```

Expected: integrity PASS (with the §2.1 log entry covering the helper-testset deletions). Test suite passes (count drops by the number of deleted spec testsets, ~50). LOC line for the commit shows the deletion delta.

If `bash scripts/check_test_integrity.sh main` reports unexplained deletions, fix them: either the migration was incomplete (recover the missed parallel test) or the §2.1 log needs more entries.

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl docs/superpowers/refactor-deleted-tests.md test/runtests.jl
git commit -m "$(cat <<'EOF'
Stage 6β.10: delete spec types + overloads after Mechanism migration

Removes MechanismSpec, StepSpec, AllostericMechanismSpec,
AbstractMechanismSpec; six _expand_* legacy overloads on
AbstractMechanismSpec; dedup!/canonicalize/dedup_key spec branches;
_init_mechanism_specs / _spec_from_mechanism / _mechanism_from_spec
helpers; spec-construction testset helpers in
test/test_mechanism_enumeration.jl; and two spec-helper testsets per
spec §2.1 narrow exception (log entry in
docs/superpowers/refactor-deleted-tests.md).

Mechanism-form parallel testsets land in 6β.1–6β.9; this commit is
the deletion half of the migrate-then-delete pattern.

src delta: -X / +Y net Z, cumulative: ±W
EOF
)"
```

- [ ] **Step 7: Edit the log to fill in the commit SHA**

```bash
git rev-parse HEAD
# Edit docs/superpowers/refactor-deleted-tests.md, replace "TBD-after-commit"
# with the SHA. Then commit:
git add docs/superpowers/refactor-deleted-tests.md
git commit -m "Stage 6β.10: backfill deleted-tests log commit SHA

src delta: 0 / 0 (docs only)"
```

- [ ] **Step 8: Tag the stage**

```bash
git tag stage-6beta-complete
```

---

# Stage 6 — `identify_rate_equation` canonicalizer rewrite

**Stage goal:** replace the regex-pipeline canonical-hash in `src/identify_rate_equation.jl` with a structural hash over the `Mechanism` form.

**Stage LOC delta target:** −300 to −500 in `src/identify_rate_equation.jl`.

**Reference:** The original implementation plan §6 — `docs/superpowers/plans/2026-05-20-concrete-types-refactor.md` lines 3027–3450 — describes this stage in detail. Execute as written there.

- [ ] **Step 1: Re-read the original Stage 6 plan**

```bash
sed -n '3027,3450p' docs/superpowers/plans/2026-05-20-concrete-types-refactor.md
```

- [ ] **Step 2–N: Execute the original Stage 6 tasks (6.1, 6.2, ...) verbatim**

Each task includes its own TDD steps and commit. Tag after the final task:

```bash
git tag stage-6-complete
```

- [ ] **Step N+1: Run final verification**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -10
wc -l src/*.jl
```

Expected: all green; cumulative LOC delta ≈ +1,500 (start of stage was ≈ +1,900; this stage drops ≈ −400).

---

# Stage 7a — Chokepoint hygiene

**Stage goal:** every parameter-name rendering in `src/` flows through `name(p::Parameter, m::Mechanism)`. After this stage the chokepoint is the sole `Parameter → Symbol` mapping; a regression test enforces this for future commits.

**Stage LOC delta target:** ≈ 0 net.

## Task 7a.1: Route `Symbol("K$idx")` sites in `src/rate_eq_derivation.jl` through the chokepoint

**Files:**
- Modify: `src/rate_eq_derivation.jl` lines 132, 148, 627, 1632

Inventory:
- Line 132: `rename[Symbol("K$idx")] = Symbol("K$rep")` — builds a rename map from one positional K name to another.
- Line 148: `sym = Symbol("K$idx")` — looks up a K param.
- Line 627: `pass1_rename[Symbol("K$idx")] = Symbol("K$rep")` — another rename map entry.
- Line 1632: same pattern as 627.

These all build kinetic-group rename maps for the polynomial substitution passes. The Parameter struct that produces `Symbol("K$idx")` is `K(group_idx)`-style — verify by reading `src/types.jl` around the K Parameter definition and the chokepoint `name(::K, m)`.

- [ ] **Step 1: Read each call site to understand context**

```bash
sed -n '125,155p' src/rate_eq_derivation.jl
sed -n '620,635p' src/rate_eq_derivation.jl
sed -n '1625,1640p' src/rate_eq_derivation.jl
```

- [ ] **Step 2: Find the K Parameter constructor and chokepoint signature**

```bash
grep -n "^struct K\|name(::K\|name(p::K" src/types.jl
```

- [ ] **Step 3: Write a regression test that fails before the change**

In a new file `test/test_chokepoint.jl`:

```julia
using Test
using EnzymeRates

@testset "chokepoint: no direct Symbol(\"K\$idx\") construction outside name()" begin
    src_dir = joinpath(dirname(@__DIR__), "src")
    for f in readdir(src_dir; join=true)
        endswith(f, ".jl") || continue
        content = read(f, String)
        # Allow Symbol("...") calls inside the body of name() methods.
        # Strategy: a Symbol("K…") line outside a function whose signature
        # contains `::Parameter` or `name(p::` is a violation.
        # Quick first pass: just count the regex matches per file.
        matches = collect(eachmatch(r"Symbol\(\"K\d|Symbol\(\"k\d|Symbol\(\"V|Symbol\(\"L\"", content))
        # After 7a, the only acceptable matches should be inside name() bodies.
        # Use a stricter line-by-line check that skips name-method bodies.
        offending = String[]
        in_name_method = false
        depth = 0
        for line in eachline(IOBuffer(content))
            if occursin(r"^function name\(.*Parameter", line) || occursin(r"^name\(.*Parameter.*\) =", line)
                in_name_method = true
                depth = 1
            elseif in_name_method
                depth += count(==('('), line) - count(==(')'), line) +
                        count(==('{'), line) - count(==('}'), line)
                if occursin(r"^end\s*$", line) || depth <= 0
                    in_name_method = false
                end
            end
            if !in_name_method && occursin(r"Symbol\(\"K\d|Symbol\(\"k\d[fr]|Symbol\(\"K_|Symbol\(\"V|Symbol\(\"L\"", line)
                push!(offending, "$(basename(f)): $line")
            end
        end
        @test isempty(offending) || (println("offending:\n", join(offending, "\n")); false)
    end
end
```

(This is a heuristic test — refine if false positives appear. The key invariant: no top-level direct Symbol("K…") construction outside `name(p, m)` bodies.)

Add to `test/runtests.jl`:

```julia
include("test_chokepoint.jl")
```

- [ ] **Step 4: Run the test to confirm it fails before the fix**

```bash
julia --project=test -e 'include("test/test_chokepoint.jl")' 2>&1 | tail -20
```

Expected: FAIL with offending lines from `rate_eq_derivation.jl` and `types.jl`.

- [ ] **Step 5: Fix line 132 — rename map entry**

Read the surrounding context (lines 125–155). The variables `idx` and `rep` are integers identifying kinetic groups. The Mechanism `m` is in scope. The fix:

```julia
# Before:
rename[Symbol("K$idx")] = Symbol("K$rep")

# After:
rename[name(K(idx), m)] = name(K(rep), m)
```

(Verify `K(group_idx::Int)` is the correct constructor — read the K struct definition in `src/types.jl`. If the field is named differently, adjust.)

- [ ] **Step 6: Fix line 148 — symbol lookup**

```julia
# Before:
sym = Symbol("K$idx")

# After:
sym = name(K(idx), m)
```

- [ ] **Step 7: Fix lines 627, 1632 — rename map entries**

Same pattern as line 132.

- [ ] **Step 8: Re-run the regression test to confirm it passes**

```bash
julia --project=test -e 'include("test/test_chokepoint.jl")' 2>&1 | tail -10
```

Expected: PASS (offending lines now empty for `rate_eq_derivation.jl`; `types.jl` still has its 2 lines, fixed in 7a.2).

- [ ] **Step 9: Run the full suite + integrity check**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -10
```

Expected: full suite passes; the 4 changed lines should produce identical Symbols to before (verify by running the rate-equation tests — they compare exact Symbols).

- [ ] **Step 10: Commit**

```bash
git add src/rate_eq_derivation.jl test/test_chokepoint.jl test/runtests.jl
git commit -m "$(cat <<'EOF'
Stage 7a.1: route Symbol(\"K\$idx\") sites in rate_eq_derivation.jl through name(K, m)

Four call sites (lines 132, 148, 627, 1632) that built positional K
param names directly now call name(K(group_idx), m). Identical output;
removes string concatenation in favor of struct-driven rendering.

Adds test/test_chokepoint.jl regression test: no direct
Symbol(\"K...\") / Symbol(\"k...\") / Symbol(\"V...\") / Symbol(\"L\")
construction outside of name(p, m) bodies.

src delta: -X / +Y net Z, cumulative: ±W
EOF
)"
```

## Task 7a.2: Route regulator-K sites in `src/types.jl` through the chokepoint

**Files:**
- Modify: `src/types.jl` lines 1840, 1841

These two lines build `K_<lig>_T_reg<n>` / `K_<lig>_reg<n>` Symbols inline. The Parameter struct that produces these is `Kreg(site_idx, ligand_name, state)`.

- [ ] **Step 1: Read the call site**

```bash
sed -n '1830,1850p' src/types.jl
```

Identify the function this is inside — likely an accessor or rendering helper for an `AllostericEnzymeMechanism`. Note the variables in scope: `site_idx`, `lig_name`, `p.state` (`:T` or `:R`).

- [ ] **Step 2: Locate the Kreg struct + name method**

```bash
grep -n "^struct Kreg\|name(::Kreg\|name(p::Kreg" src/types.jl
```

- [ ] **Step 3: Replace the two lines**

```julia
# Before:
p.state === :T ? Symbol("K_$(lig_name)_T_reg$site_idx") :
                 Symbol("K_$(lig_name)_reg$site_idx")

# After:
name(Kreg(site_idx, lig_name, p.state), m)
```

(Verify the `Kreg` constructor signature — read the struct definition. The arguments may need reordering.)

- [ ] **Step 4: Run the regression test**

```bash
julia --project=test -e 'include("test/test_chokepoint.jl")' 2>&1 | tail -10
```

Expected: PASS (no offending lines anywhere).

- [ ] **Step 5: Run the full suite**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -10
```

Expected: full suite passes; allosteric rate-equation tests in particular verify the regulator K names match exactly.

- [ ] **Step 6: Commit**

```bash
git add src/types.jl
git commit -m "$(cat <<'EOF'
Stage 7a.2: route regulator-K Symbol construction through name(Kreg, m)

Two lines in src/types.jl that built K_<lig>_(T_)reg<n> names directly
now call name(Kreg(site, lig, state), m). Chokepoint is now the sole
Parameter → Symbol mapping in src/.

src delta: -X / +Y net Z, cumulative: ±W
EOF
)"
```

- [ ] **Step 7: Tag the stage**

```bash
git tag stage-7a-complete
```

---

# Stage 7b — DSL native emission with decomposed Species

**Stage goal:** macros `@enzyme_mechanism` and `@allosteric_mechanism` emit `EnzymeMechanism(Mechanism(...))` with decomposed-Species notation. Opaque-form helpers (`_form_name`, `_parse_bound`, etc.) deleted at the end.

**Stage LOC delta target:** −200 to −250 (DSL rewrite ≈ neutral; opaque-helper deletion ≈ −200; `_mechanism_from_legacy_sig` deletion ≈ −60).

## Task 7b.1: Extend DSL parser to accept decomposed-Species notation

**Files:**
- Modify: `src/dsl.jl` lines 391–620 (`@enzyme_mechanism` parser) and lines 669–880 (`@allosteric_mechanism` parser)

Today the step-grammar parses each step's LHS and RHS as plain Symbols representing opaque enzyme-form names (`:E`, `:ES`, `:EP`). The extension: also accept `Call` expressions like `E(S)` parsed as `Species([Substrate(:S)], :E)`.

Parser strategy:
- If the step entry is a `Symbol` (e.g., `:ES`): parse as opaque form (legacy path).
- If the step entry is an `Expr(:call, :E, :S)` (i.e., `E(S)`): parse as decomposed-Species `Species([Substrate(:S)], :E)`. The metabolite role (`Substrate`/`Product`/`CompetitiveInhibitor`/`AllostericRegulator`) is determined by looking up the metabolite name in the surrounding `substrates:` / `products:` / `competitive_inhibitors:` / `allosteric_regulators:` blocks.

- [ ] **Step 1: Read the current step parser**

```bash
sed -n '440,560p' src/dsl.jl
```

Identify the function that walks the `steps: begin ... end` block and locate the per-entry parsing.

- [ ] **Step 2: Write a failing parser test**

Add to `test/test_dsl.jl` (or wherever DSL tests live):

```julia
@testset "DSL: decomposed-Species notation parses to Mechanism" begin
    m = @enzyme_mechanism begin
        substrates: S[C6H12O6]
        products:   P[C6H13O9P]
        steps: begin
            E + S => E(S)
            E(S) => E(P)
            E(P) => E + P
        end
    end
    mech = Mechanism(m)
    @test length(mech.steps) == 3
    # Each step should have decomposed Species — bound list non-empty for E(S)/E(P)
    es_step = first(mech.steps)[1]  # E + S => E(S)
    @test name(es_step.to_species) == :E
    @test bound(es_step.to_species) == [Substrate(:S)]
end
```

- [ ] **Step 3: Run the test to verify it fails**

```bash
julia --project=test -e 'include("test/test_dsl.jl")' 2>&1 | tail -20
```

Expected: FAIL — current parser doesn't recognize `E(S)`.

- [ ] **Step 4: Implement the parser extension**

In the step-entry parser, add a dispatch on `Expr` head:

```julia
function _parse_step_entry(entry, declared_mets::Dict{Symbol, Metabolite})
    if entry isa Symbol
        # Legacy opaque-form path: name as-is.
        return Species(Metabolite[], entry)
    elseif entry isa Expr && entry.head === :call
        # Decomposed: E(S, ATP) → Species([Substrate(:S), Substrate(:ATP)], :E)
        conformation = entry.args[1]::Symbol
        ligand_names = Symbol[a for a in entry.args[2:end]]
        ligands = Metabolite[declared_mets[n] for n in ligand_names]
        return Species(ligands, conformation)
    else
        error("@enzyme_mechanism: unrecognized step entry: $entry")
    end
end
```

The `declared_mets::Dict{Symbol, Metabolite}` is built from the top-level `substrates:` / `products:` / etc. blocks of the macro body. Look up each name to get the right Metabolite subtype.

Update the step-parsing function to dispatch through `_parse_step_entry` and emit appropriate `Step(...)` calls with decomposed Species.

- [ ] **Step 5: Make the macros emit `EnzymeMechanism(Mechanism(...))` when all step entries are decomposed**

In the macro return statement, branch on whether any step entry was opaque-form. If all decomposed, emit `EnzymeMechanism(Mechanism($rxn, $grouped_steps_expr))`. Otherwise fall back to the legacy `EnzymeMechanism($mets_expr, $rxns_expr)` form.

- [ ] **Step 6: Run the test to verify it passes**

```bash
julia --project=test -e 'include("test/test_dsl.jl")' 2>&1 | tail -10
```

Expected: PASS.

- [ ] **Step 7: Run the full suite to confirm legacy fixtures still work**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -10
```

Expected: all existing tests pass — dual-grammar support means legacy fixtures continue to work alongside the new test.

- [ ] **Step 8: Commit**

```bash
git add src/dsl.jl test/test_dsl.jl
git commit -m "$(cat <<'EOF'
Stage 7b.1: DSL parser accepts decomposed-Species notation E(S)

@enzyme_mechanism and @allosteric_mechanism parsers now recognize
Expr(:call, conformation, ligands...) syntax as a Species with the
ligands as bound metabolites. Opaque-form Symbol entries continue to
parse as Species([], :name) — dual-grammar transition.

When ALL step entries in a macro are decomposed, the macro emits
EnzymeMechanism(Mechanism(...)) directly. Otherwise falls back to the
legacy (metabolites, reactions) emission shape.

src delta: -X / +Y net Z, cumulative: ±W
EOF
)"
```

## Task 7b.2: Migrate `test/mechanism_definitions_for_test_enzyme_derivation.jl` to decomposed-Species grammar

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`

The file contains 37 DSL invocations defining test fixtures (uni-uni, bi-bi, ter-ter, allosteric, etc.) with hand-derived analytical formulas attached. The DSL invocations migrate to decomposed grammar; the analytical formulas stay unchanged (they reference `:K1`-style names which belong to the deferred parameter-naming refactor).

- [ ] **Step 1: Inventory the DSL invocations**

```bash
grep -n "@enzyme_mechanism\|@allosteric_mechanism" test/mechanism_definitions_for_test_enzyme_derivation.jl
```

Expected: 37 occurrences.

- [ ] **Step 2: Migrate one fixture as a worked example**

Pick the uni-uni fixture (typically the first one). Translate:

```julia
# Before:
@enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S => ES
        ES => EP
        EP => E + P
    end
end

# After:
@enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S => E(S)
        E(S) => E(P)
        E(P) => E + P
    end
end
```

Translation rule: every opaque enzyme-form Symbol containing bound metabolites (e.g., `:ES`, `:EAB`, `:E_PQ`) becomes a decomposed-Species call `E(S)`, `E(A, B)`, `E(P, Q)` etc. Pure conformation Symbols (`:E`, `:Estar`) stay as bare names.

For dead-end inhibitors and allosteric regulators, decompose similarly: `:E_I` → `E(I)`, `:E_R` → `E(R)`.

For ping-pong residual forms (e.g., `:Estar`), keep the bare Symbol — that's not a bound-metabolite form, it's a conformation change.

- [ ] **Step 3: Run the test for that one fixture to confirm the migration is equivalent**

```bash
julia --project=test -e '
include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
include("test/test_rate_eq_derivation.jl")
' 2>&1 | tail -20
```

Expected: the analytical-formula comparison test for that fixture still passes (rate equation matches the hand-derived formula).

- [ ] **Step 4: Migrate the remaining 36 fixtures**

Apply the same translation rule to each. Some fixtures (allosteric, dead-end, multi-substrate) have more complex step lists — be careful to identify which Symbols are bound-metabolite forms vs pure conformations.

For each fixture, after migration:
1. Verify the rate-equation-string test for that fixture passes.
2. If the test fails, the migration is wrong — STOP and investigate before continuing.

- [ ] **Step 5: Run the full suite + integrity check**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -10
```

Expected: all pass. Test count unchanged.

- [ ] **Step 6: Commit**

```bash
git add test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -m "$(cat <<'EOF'
Stage 7b.2: migrate mechanism_definitions fixtures to decomposed-Species DSL grammar

37 @enzyme_mechanism / @allosteric_mechanism invocations switched from
opaque-form Symbol entries (:ES, :EAB) to decomposed-Species notation
(E(S), E(A, B)). Analytical formulas (which reference :K1-style param
names) unchanged — they belong to the deferred parameter-naming
refactor.

src delta: 0 / 0 (test-fixture only)
EOF
)"
```

## Task 7b.3: Migrate remaining DSL-using test files

**Files:**
- Modify: `test/test_dsl.jl` (36 DSL invocations)
- Modify: `test/test_types.jl` (15)
- Modify: `test/test_rate_eq_derivation.jl` (14)
- Modify: `test/test_accessors.jl` (3)
- Modify: `test/test_identify_rate_equation.jl` (3)
- Modify: `test/test_fitting.jl` (2)
- Modify: `test/test_compile_budget.jl` (1)

- [ ] **Step 1: Migrate one file at a time**

For each file, apply the same translation rule as Task 7b.2. After each file, run that file's tests to confirm green:

```bash
julia --project=test -e 'include("test/<filename>")' 2>&1 | tail -10
```

- [ ] **Step 2: Run the full suite after each file**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -10
```

- [ ] **Step 3: Commit per file (or per logical group of files)**

```bash
git add test/test_dsl.jl
git commit -m "Stage 7b.3a: migrate test_dsl.jl fixtures to decomposed grammar

src delta: 0 / 0 (test-only)"

git add test/test_types.jl test/test_accessors.jl
git commit -m "Stage 7b.3b: migrate test_types.jl + test_accessors.jl fixtures

src delta: 0 / 0 (test-only)"

# ... and so on
```

Group small files (`test_accessors.jl`, `test_identify_rate_equation.jl`, `test_fitting.jl`, `test_compile_budget.jl`) into one or two commits if they all change together cleanly. Keep each commit focused.

## Task 7b.4: Rewrite `_expand_*` moves to read decomposed Species fields directly

**Files:**
- Modify: `src/mechanism_enumeration.jl` — the 6 native `_expand_*` move implementations

After Tasks 7b.2 and 7b.3, every fixture in the test suite produces decomposed Species. The `_expand_*` moves can stop calling `_form_name` / `_parse_bound` / `_bound_mets_from_form_name` and instead read `from_species.bound` / `to_species.bound` directly.

- [ ] **Step 1: Audit each `_expand_*` move for opaque-helper usage**

```bash
grep -n "_form_name\|_parse_bound\|_bound_mets_from_form_name\|_dead_end_form_name" src/mechanism_enumeration.jl
```

For each match, identify the calling `_expand_*` function.

- [ ] **Step 2: For each `_expand_*` function, rewrite the opaque-helper call to read structured Species fields**

Worked example: a move that asks "which metabolites are bound in this Species?"

```julia
# Before:
bound_mets = _bound_mets_from_form_name(name(species), reaction)

# After:
bound_mets = bound(species)  # Species.bound is the Vector{Metabolite}
```

A move that builds a new Species name from a bound list:

```julia
# Before:
new_name = _form_name(bound_mets, products, reaction)

# After:
new_species = Species(bound_mets, :E)  # or appropriate conformation symbol
```

(For empty-residual ping-pong / Estar forms, the conformation symbol comes from the existing Species; pass it through.)

- [ ] **Step 3: Use existing enumeration tests as regression coverage**

This is the critical TDD step. You don't need a new test — `test/test_mechanism_enumeration.jl` already has `@testset "Integration"` blocks that pin full-enumeration counts for uni-uni, bi-bi, allosteric, and dead-end cases. Run those before AND after each `_expand_*` rewrite:

```bash
julia --project=test -e 'include("test/test_mechanism_enumeration.jl")' 2>&1 | grep -E "Integration|Test (Pass|Fail)" | head -20
```

Every count assertion must remain unchanged. If a count changes, the rewrite is wrong — STOP and find the structural-field reading bug. Don't update the count.

- [ ] **Step 4: Run the full suite after each move rewrite**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -10
```

- [ ] **Step 5: Commit per move (6 commits total)**

```bash
git add src/mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Stage 7b.4a: _expand_re_to_ss reads structured Species fields directly

Replaces _form_name / _parse_bound / _bound_mets_from_form_name calls
with direct access to species.bound. Identical enumeration output —
verified by full test suite.

src delta: -X / +Y net Z, cumulative: ±W
EOF
)"

# Repeat for: _expand_split_kinetic_group, _expand_add_dead_end_regulator,
# _expand_to_allosteric, _expand_add_allosteric_regulator, _expand_change_allo_state
```

## Task 7b.5: Delete the dual-grammar legacy branch + emission path

**Files:**
- Modify: `src/dsl.jl` — remove the opaque-form fallback emission

After 7b.2 and 7b.3, no fixture in the test suite produces opaque-form entries. The `@enzyme_mechanism` / `@allosteric_mechanism` macros' legacy emission branch can be deleted.

- [ ] **Step 1: Confirm no test fixture still uses opaque form**

```bash
grep -rn "@enzyme_mechanism\|@allosteric_mechanism" test/ | xargs grep -l "=> [A-Z][A-Z_]*$" 2>/dev/null | head
```

Expected: empty (no fixture has a bare-name multi-char RHS in a `=>` step).

If non-empty, those fixtures need migration first (back to Task 7b.2/7b.3).

- [ ] **Step 2: Delete the opaque-form parser branch in `_parse_step_entry`**

Remove the `if entry isa Symbol` arm (or simplify it to only accept pure-conformation single-character symbols like `:E`, `:Estar`). The new entry parser expects every bound-metabolite form to be a `Call` expression.

- [ ] **Step 3: Delete the legacy macro emission path**

Both macros now unconditionally emit `EnzymeMechanism(Mechanism(...))`.

- [ ] **Step 4: Delete the 2-arg `EnzymeMechanism(metabolites, reactions)` constructor in `src/types.jl`**

Lines 749–870 (approximate; verify). The new emission path doesn't call it.

- [ ] **Step 5: Delete `_mechanism_from_legacy_sig` and the `_is_new_sig` branch in `Mechanism(em::EnzymeMechanism{Sig})`**

`src/types.jl` lines 785–860 (approximate). After deletion, `Mechanism(em::EnzymeMechanism{Sig})` is a single-body function: `_mechanism_from_sig(Sig)`.

- [ ] **Step 6: Delete the opaque-form helpers**

In `src/mechanism_enumeration.jl`:
- `_form_name`
- `_parse_bound`
- `_bound_mets_from_form_name`
- `_dead_end_form_name`
- `_atoms_dict`
- `_is_estar_form`
- `_can_pingpong`
- `_subtract_atoms`

For each, verify zero remaining callers via:

```bash
grep -rn "<function_name>" src/ test/
```

If any caller remains, fix the caller first (likely a leftover `_expand_*` move that 7b.4 missed).

- [ ] **Step 7: Run the full suite + integrity check + perf gates**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -30
wc -l src/*.jl
```

Expected: all green. Stage 7b cumulative LOC delta ≈ −200 to −250.

- [ ] **Step 8: Commit**

```bash
git add src/dsl.jl src/types.jl src/mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Stage 7b.5: delete opaque-form helpers + legacy DSL emission path

Removes _form_name, _parse_bound, _bound_mets_from_form_name,
_dead_end_form_name, _atoms_dict, _is_estar_form, _can_pingpong,
_subtract_atoms from src/mechanism_enumeration.jl. Removes the
2-arg EnzymeMechanism(metabolites, reactions) constructor and
_mechanism_from_legacy_sig + the _is_new_sig branch from src/types.jl.
Removes the dual-grammar fallback in src/dsl.jl.

After this commit, the only path from DSL to Mechanism is
EnzymeMechanism(Mechanism(...)) with decomposed Species.

src delta: -X / +Y net Z, cumulative: ±W
EOF
)"
```

- [ ] **Step 9: Tag the stage**

```bash
git tag stage-7b-complete
```

---

# Stage 7c — Enumeration dispatches on `EnzymeReaction`

**Stage goal:** the ~20 sites in `src/mechanism_enumeration.jl` that dispatch on `EnzymeReactionLegacy` move to dispatch on `EnzymeReaction`. After this stage, no internal call path needs `_to_legacy_reaction`.

**Stage LOC delta target:** small net (dispatch swap; deletions land in 7d).

## Task 7c.1: Inventory `EnzymeReactionLegacy` dispatch sites

**Files:**
- Read: `src/mechanism_enumeration.jl`, `src/identify_rate_equation.jl`

- [ ] **Step 1: List all dispatch sites**

```bash
grep -nE "::EnzymeReactionLegacy|<:EnzymeReactionLegacy" src/*.jl
```

Expected output: ~20 method signatures in `src/mechanism_enumeration.jl`, plus 2 in `src/identify_rate_equation.jl` (struct type parameter + outer `identify_rate_equation` signature).

- [ ] **Step 2: Verify `EnzymeReaction` is structurally compatible**

For each dispatched function, check what fields it reads from the reaction argument. Read the same fields off `EnzymeReaction` (the concrete struct) and confirm equivalence. Common access patterns: `substrates(r)`, `products(r)`, `regulators(r)`, `oligomeric_state(r)`.

If any access pattern is `EnzymeReactionLegacy`-specific (e.g., reads a type parameter directly via `where {S, P, R, N}`), STOP — that's a deeper rewrite, ask before proceeding.

- [ ] **Step 3: Note the inventory in a working file (mental or scratch — no commit)**

## Task 7c.2: Swap dispatch from `EnzymeReactionLegacy` to `EnzymeReaction`

**Files:**
- Modify: `src/mechanism_enumeration.jl` (~20 method signatures)
- Modify: `src/identify_rate_equation.jl` (struct parametric type + 1 outer function signature)

- [ ] **Step 1: Write a failing test**

The current test suite calls `init_mechanisms(rxn::EnzymeReaction)` (the DSL-produced concrete reaction) which goes through `_to_legacy_reaction` → `init_mechanisms(::EnzymeReactionLegacy)`. After this task, the concrete `EnzymeReaction` overload is the direct dispatch target.

Write a test that asserts no `_to_legacy_reaction` call happens for the common path:

```julia
@testset "init_mechanisms(::EnzymeReaction) dispatches directly (no Legacy conversion)" begin
    rxn = @enzyme_reaction begin
        substrates: S[C6H12O6]
        products:   P[C6H13O9P]
    end
    # The method that takes EnzymeReaction should be the one called —
    # verify by checking that the dispatched method's signature is
    # EnzymeReaction, not EnzymeReactionLegacy:
    method = which(EnzymeRates.init_mechanisms, (typeof(rxn),))
    @test occursin("EnzymeReaction", string(method.sig))
    @test !occursin("EnzymeReactionLegacy", string(method.sig))
end
```

- [ ] **Step 2: Run the test to verify it fails**

Expected: FAIL — current dispatch goes through `EnzymeReactionLegacy`.

- [ ] **Step 3: Change each dispatch site from `::EnzymeReactionLegacy` to `::EnzymeReaction`**

For each method signature:

```julia
# Before:
function init_mechanisms(@nospecialize(reaction::EnzymeReactionLegacy); kwargs...)

# After:
function init_mechanisms(reaction::EnzymeReaction; kwargs...)
```

Drop `@nospecialize` — the concrete struct doesn't need it (no per-arity type explosion). Or keep it if the original had a reason (read the docstring comment if any).

Inside each function body, replace any `EnzymeReactionLegacy`-specific accessor calls with their `EnzymeReaction` equivalents (likely no changes needed because the accessors are polymorphic).

- [ ] **Step 4: Update `IdentifyRateEquationProblem` type parameter**

```julia
# Before (src/identify_rate_equation.jl line ~21):
struct IdentifyRateEquationProblem{R<:EnzymeReactionLegacy, D<:NamedTuple}

# After:
struct IdentifyRateEquationProblem{R<:EnzymeReaction, D<:NamedTuple}
```

And the outer constructor signature:

```julia
# Before:
function identify_rate_equation(reaction::EnzymeReactionLegacy, table; Keq::Real)

# After:
function identify_rate_equation(reaction::EnzymeReaction, table; Keq::Real)
```

- [ ] **Step 5: Update `_init_mechanism_specs` signatures**

```bash
grep -n "_init_mechanism_specs" src/mechanism_enumeration.jl
```

The `EnzymeReactionLegacy` overload (line ~3414) can be deleted; the `EnzymeReaction` overload (line ~3411) becomes the sole method.

- [ ] **Step 6: Delete remaining `_to_legacy_reaction` call sites**

```bash
grep -rn "_to_legacy_reaction" src/
```

Each remaining call should be removable — the dispatch swap above means callers can pass the concrete `EnzymeReaction` directly. The function itself stays for now (deleted in 7d).

- [ ] **Step 7: Re-baseline compile-budget tests if needed**

Read `test/test_compile_budget.jl` lines 77, 147 etc. — these construct `EnzymeReactionLegacy` directly to exercise per-arity specialization paths. Change them to construct `EnzymeReaction` instead:

```julia
# Before:
r = EnzymeRates.EnzymeReactionLegacy(((:S, ((:C, 6), (:H, 12), (:O, 6))),),
                                       ((:P, ((:C, 6), (:H, 13), (:O, 9), (:P, 1))),))

# After:
r = @enzyme_reaction begin
    substrates: S[C6H12O6]
    products:   P[C6H13O9P]
end
```

Run the compile-budget tests:

```bash
julia --project=test -e 'include("test/test_compile_budget.jl")' 2>&1 | tail -20
```

If the trace-compile counts are now significantly lower than the budget (e.g., baseline drops from 29 to 12), reduce the budgets to keep them tight at 2× the new baseline:

```julia
const INIT_TRACE_BUDGET = 24    # baseline 2026-05-22: 12; budget = 2× (rounded up)
```

If any budget would NEED to be raised (test_compile_budget gets WORSE after the dispatch swap), STOP and investigate — this would indicate that the concrete `EnzymeReaction` is somehow producing more specialized methods than the singleton, which would be unexpected. Don't raise the budget; redesign the change.

- [ ] **Step 8: Run the regression test + full suite**

```bash
julia --project=test -e 'include("test/test_chokepoint.jl")' 2>&1 | tail -5
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -30
```

Expected: all green. The dispatch-confirmation test passes; compile-budget gates pass at the new (possibly tightened) limits.

- [ ] **Step 9: Commit**

```bash
git add src/mechanism_enumeration.jl src/identify_rate_equation.jl test/test_compile_budget.jl test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Stage 7c: enumeration dispatches on EnzymeReaction (no EnzymeReactionLegacy)

~20 method signatures in src/mechanism_enumeration.jl moved from
::EnzymeReactionLegacy to ::EnzymeReaction. IdentifyRateEquationProblem
type parameter switched. test_compile_budget.jl now constructs
EnzymeReaction directly via @enzyme_reaction; budgets re-baselined to
the new (lower) trace-compile counts.

_to_legacy_reaction has no remaining callers in src/; the function
itself is deleted in 7d.

src delta: -X / +Y net Z, cumulative: ±W
EOF
)"
```

- [ ] **Step 10: Tag the stage**

```bash
git tag stage-7c-complete
```

---

# Stage 7d — Delete legacy infrastructure

**Stage goal:** delete `EnzymeReactionLegacy`, its accessors, `_to_legacy_reaction`, the dual-Sig branches in `src/types.jl` accessors, and any remaining adapter code.

**Stage LOC delta target:** −450 to −550.

## Task 7d.1: Delete `EnzymeReactionLegacy` struct + accessors + `_to_legacy_reaction`

**Files:**
- Modify: `src/types.jl`
- Modify: `test/test_types.jl` — delete `EnzymeReactionLegacy`-specific testsets (§2.1 narrow exception)
- Modify: `test/runtests.jl` — drop `EnzymeReactionLegacy` from `using EnzymeRates: ...`

- [ ] **Step 1: Confirm zero callers in `src/`**

```bash
grep -rn "EnzymeReactionLegacy\|_to_legacy_reaction" src/
```

Expected: zero matches (after Stage 7c).

If any matches remain, STOP — fix the caller (likely a stray leftover from 7c) before deleting.

- [ ] **Step 2: Identify the `EnzymeReactionLegacy`-specific testsets in `test/test_types.jl`**

```bash
grep -n "EnzymeReactionLegacy" test/test_types.jl
```

Expected: ~15 testsets exercising the singleton-typed reaction's behavior (atom balance, canonical ordering, regulator overlap, oligomeric_state, etc.).

These need either (a) parallel `EnzymeReaction` testsets if not already in `test/test_types.jl`, or (b) §2.1 log entries if their behavior is already covered.

- [ ] **Step 3: Audit coverage parity**

For each `EnzymeReactionLegacy` testset, search `test/test_types.jl` for an `EnzymeReaction` parallel. If missing, add one BEFORE the deletion commit. The atom-balance check, duplicate-name check, canonical ordering, etc. should all have `EnzymeReaction` parallels because `EnzymeReaction` enforces the same invariants.

For testsets that genuinely only make sense in the singleton-typed form (e.g., "type parameter encoding"), add a §2.1 log entry.

- [ ] **Step 4: Add §2.1 log entries for genuine narrow-exception deletions**

Edit `docs/superpowers/refactor-deleted-tests.md`:

```markdown
## Stage 7d.1 (commit TBD-after-commit)

### test_types.jl `@testset "EnzymeReactionLegacy ..."` (15 testsets)
- Helper deleted: `EnzymeReactionLegacy` struct + outer constructor + accessors
- Replacement: `EnzymeReaction` (concrete struct, same invariants)
- Coverage parity:
  - "EnzymeReactionLegacy" atom mandatory → "EnzymeReaction atom mandatory" (test_types.jl @testset)
  - "EnzymeReactionLegacy" atom balance → "EnzymeReaction atom balance" (test_types.jl @testset)
  - "EnzymeReactionLegacy duplicate substrate names" → covered by EnzymeReaction constructor's identical validation
  - ... (one row per deleted testset)
```

- [ ] **Step 5: Delete the `EnzymeReactionLegacy` struct + accessors + helpers**

In `src/types.jl`:
- Delete `struct EnzymeReactionLegacy{S,P,R,N}` (line 628) and its docstring.
- Delete the outer constructor `EnzymeReactionLegacy(subs, prods, regs; oligomeric_state)` (line 654–693).
- Delete `_to_legacy_reaction(r::EnzymeReaction)` (line 697–730).
- Delete the `Base.show(io, ::EnzymeReactionLegacy)` method (line 1221).
- Delete the accessor methods on `EnzymeReactionLegacy`:
  - `substrates(::EnzymeReactionLegacy)` (line 1451)
  - `products(::EnzymeReactionLegacy)` (line 1462)
  - `regulators(::EnzymeReactionLegacy)` (line 1472)
  - `regulator_roles(::EnzymeReactionLegacy)` (line 1476)
  - `oligomeric_state(::EnzymeReactionLegacy)` (line 1479)

Also delete `_sum_atoms` IF it's only called by `EnzymeReactionLegacy`'s constructor — check with `grep -n "_sum_atoms" src/`. If the new `EnzymeReaction` constructor uses it too, leave it.

In `EnzymeRates.jl`: remove `EnzymeReactionLegacy` from the `export` list if it's exported.

- [ ] **Step 6: Delete the `EnzymeReactionLegacy`-specific testsets in `test/test_types.jl`**

Bottom-up. For each `@testset "EnzymeReactionLegacy ..."`, delete the block.

- [ ] **Step 7: Update `test/runtests.jl`**

```julia
# Before:
using EnzymeRates: EnzymeReactionLegacy

# After: (remove the line entirely if EnzymeReactionLegacy was its only import)
```

- [ ] **Step 8: Run full suite + integrity check**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -30
wc -l src/*.jl
```

Expected: all green. The §2.1 log covers the EnzymeReactionLegacy testset deletions.

- [ ] **Step 9: Commit**

```bash
git add src/types.jl src/EnzymeRates.jl test/test_types.jl test/runtests.jl docs/superpowers/refactor-deleted-tests.md
git commit -m "$(cat <<'EOF'
Stage 7d.1: delete EnzymeReactionLegacy + accessors + _to_legacy_reaction

Removes the singleton-typed parametric reaction struct and its
accessor methods from src/types.jl. _to_legacy_reaction adapter is
gone — no remaining callers after Stage 7c. Removes EnzymeReactionLegacy
testsets from test/test_types.jl per §2.1 narrow exception (log entry
records coverage parity with EnzymeReaction testsets).

src delta: -X / +Y net Z, cumulative: ±W
EOF
)"
```

Backfill the commit SHA in `docs/superpowers/refactor-deleted-tests.md`:

```bash
git rev-parse HEAD
# Edit doc, then:
git add docs/superpowers/refactor-deleted-tests.md
git commit -m "Stage 7d.1: backfill deleted-tests log commit SHA

src delta: 0 / 0 (docs only)"
```

## Task 7d.2: Collapse dual-Sig branches in `EnzymeMechanism` accessors

**Files:**
- Modify: `src/types.jl` — 12 accessors with `_is_new_sig` branches

After 7b.5 deleted `_mechanism_from_legacy_sig`, the `_is_new_sig` distinction is moot — every `EnzymeMechanism{Sig}` now has a new-shape Sig. The dual-branch accessors collapse to single bodies.

- [ ] **Step 1: Inventory dual-Sig accessors**

```bash
grep -n "_is_new_sig" src/types.jl
```

Expected: ~12 functions with `_is_new_sig(Sig) ? ... : ...` ternaries.

- [ ] **Step 2: Write a regression test**

Pick a non-trivial accessor (e.g., `metabolites(::EnzymeMechanism)`) and assert it returns the same value before and after the collapse for a representative fixture. (Probably this is already covered by existing tests; if so, skip this step.)

- [ ] **Step 3: For each accessor, delete the legacy branch + keep the new-Sig body**

```julia
# Before:
function metabolites(::EnzymeMechanism{Sig}) where {Sig}
    _is_new_sig(Sig) ? <new_body> : <legacy_body>
end

# After:
function metabolites(::EnzymeMechanism{Sig}) where {Sig}
    <new_body>
end
```

Also delete the `_is_new_sig` helper itself (it has no remaining callers).

- [ ] **Step 4: Run full suite + integrity check + perf gates**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -10
```

Expected: all green.

- [ ] **Step 5: Commit**

```bash
git add src/types.jl
git commit -m "$(cat <<'EOF'
Stage 7d.2: collapse dual-Sig branches in EnzymeMechanism accessors

After 7b.5 retired _mechanism_from_legacy_sig, the _is_new_sig
distinction is moot — every EnzymeMechanism{Sig} carries new-shape
Sig. 12 accessor functions in src/types.jl collapse from
ternary-bodied to single-bodied. _is_new_sig helper deleted.

src delta: -X / +Y net Z, cumulative: ±W
EOF
)"
```

## Task 7d.3: Final cleanup of `_n_fit_params_estimate` + stragglers

**Files:**
- Modify: `src/mechanism_enumeration.jl` (or wherever `_n_fit_params_estimate` lives)

- [ ] **Step 1: Inventory `_n_fit_params_estimate` methods**

```bash
grep -n "function _n_fit_params_estimate\|^_n_fit_params_estimate" src/
```

Expected: at least one method on `Mechanism` and possibly leftover overloads on now-deleted spec types.

- [ ] **Step 2: Delete any leftover overloads on deleted spec types**

If `_n_fit_params_estimate(::AbstractMechanismSpec)` or similar remains, delete it. The struct itself is gone, so the overload can't dispatch anyway.

- [ ] **Step 3: Final src-side dead-code grep**

```bash
for f in src/*.jl; do
    echo "=== $f ==="
    grep -oE '^function [a-zA-Z_][a-zA-Z_0-9!]*' "$f" | sed 's/function //' | \
    while read fn; do
        count=$(grep -r "\b$fn\b" src/ | wc -l)
        [ "$count" -le 1 ] && echo "POSSIBLY UNUSED: $fn"
    done
done
```

For each "POSSIBLY UNUSED" entry: verify it's not called from tests (`grep -rn "$fn" test/`) and not part of the public API (`grep "$fn" src/EnzymeRates.jl`). If genuinely unused, delete.

- [ ] **Step 4: Run full suite + integrity**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -10
wc -l src/*.jl
```

Expected: all green. Cumulative LOC ≈ +800.

- [ ] **Step 5: Commit**

```bash
git add src/
git commit -m "$(cat <<'EOF'
Stage 7d.3: delete leftover overloads + unused helpers

src delta: -X / +Y net Z, cumulative: ±W
EOF
)"
```

- [ ] **Step 6: Tag the stage**

```bash
git tag stage-7d-complete
```

---

# Stage 7e — Cleanup + docs + PR

**Stage goal:** dead-code sweep, README + CLAUDE.md updates, final PR.

**Stage LOC delta target:** ≈ −100 (mostly cleanup; possibly 0 if everything's already tight).

## Task 7e.1: Final dead-code sweep

**Files:**
- Modify: any `src/*.jl` with detected dead code

- [ ] **Step 1: Run the dead-code detector**

```bash
for f in src/*.jl; do
    echo "=== $f ==="
    grep -oE '^function [a-zA-Z_][a-zA-Z_0-9!]*' "$f" | sed 's/function //' | sort -u | \
    while read fn; do
        callers=$(grep -rl "\b$fn\b" src/ test/ 2>/dev/null | grep -v "/$(basename $f)$" | wc -l)
        callers_in_file=$(grep -c "\b$fn\b" "$f")
        if [ "$callers" -eq 0 ] && [ "$callers_in_file" -le 1 ]; then
            echo "POSSIBLY UNUSED: $fn"
        fi
    done
done
```

- [ ] **Step 2: For each candidate, verify before deleting**

```bash
grep -rn "\b<fn>\b" src/ test/
```

Confirm zero callers outside the defining file. If it's only called from within its own file (self-recursive helper that's no longer needed), delete; otherwise keep.

Some genuine candidates after this refactor:
- Private helpers (`_*`) that were called only by now-deleted spec or legacy code.
- Inline single-use helpers — consider inlining instead of deleting.

- [ ] **Step 3: Inline single-use private helpers where it improves readability**

For each `_<helper>` called in exactly one place, consider inlining. Skip if the helper has a clear name that documents intent — readability trumps minor LOC reduction.

- [ ] **Step 4: Run full suite after each round of changes**

```bash
bash scripts/check_test_integrity.sh main
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -10
```

- [ ] **Step 5: Commit**

```bash
git add src/
git commit -m "$(cat <<'EOF'
Stage 7e.1: final dead-code sweep — delete unused / inline single-use helpers

src delta: -X / +Y net Z, cumulative: ±W
EOF
)"
```

## Task 7e.2: Update README to new DSL grammar

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Read current README**

```bash
head -200 README.md
```

- [ ] **Step 2: Update DSL examples**

Replace every `@enzyme_mechanism` and `@allosteric_mechanism` example using opaque-form (`E + S => ES`) with decomposed-Species notation (`E + S => E(S)`).

Update the architecture section to describe the final state: one `Mechanism`/`Step`/`Species` family, single `EnzymeReaction`, chokepoint `name(p, m)` for parameter rendering.

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "Stage 7e.2: update README to new DSL grammar + final architecture

src delta: 0 / 0 (docs only)"
```

## Task 7e.3: Update `.claude/CLAUDE.md`

**Files:**
- Modify: `.claude/CLAUDE.md` — "Source Layout" + "Key Architecture Decisions" sections

- [ ] **Step 1: Read current CLAUDE.md**

Already in your context. Focus on the "Source Layout", "Key Architecture Decisions", and "Mechanism enumeration building blocks" sections.

- [ ] **Step 2: Update to reflect final architecture**

Key updates:
- "Source Layout": remove references to `EnzymeReactionLegacy`, `_form_name`, `_parse_bound`, spec types. Add reference to decomposed-Species DSL grammar.
- "Key Architecture Decisions" → "EnzymeReaction{S,P,R,N}" subsection: rewrite to describe concrete `EnzymeReaction` struct with `reactants::Vector{Reactant}`, `regulators::Vector{RegulatorMult}`, etc.
- "Mechanism enumeration building blocks": remove references to `MechanismSpec`, `AllostericMechanismSpec`, `StepSpec`. Mechanism / AllostericMechanism / Step are the only forms.
- Add a new short subsection: "Parameter naming chokepoint" — explains that all Symbol rendering goes through `name(p::Parameter, m::Mechanism)` and that a future refactor will switch to semantic names.

- [ ] **Step 3: Commit**

```bash
git add .claude/CLAUDE.md
git commit -m "Stage 7e.3: update CLAUDE.md to final architecture

src delta: 0 / 0 (docs only)"
```

## Task 7e.4: Final sanity sweep + PR preparation

**Files:**
- Verify: everything

- [ ] **Step 1: Full test suite — clean run**

```bash
julia --project=test -e 'include("test/runtests.jl")' 2>&1 | tail -30
```

Expected: ALL tests pass; `test_rate_equation_performance` 0-alloc/<100ns gate green; ALL compile-budget gates green.

- [ ] **Step 2: Final test-integrity verification (NON-NEGOTIABLE per spec §2 + continuation §4)**

```bash
bash scripts/check_test_integrity.sh main
```

If FAIL: **STOP — DO NOT OPEN PR**. Audit each commit on the branch; restore the deleted/weakened tests; re-verify.

Confirm @testset count grew vs main:

```bash
n_main=0
for f in $(git ls-tree -r main --name-only test/ 2>/dev/null | grep "\.jl$"); do
    c=$(git show "main:$f" 2>/dev/null | grep -c "@testset" || true)
    n_main=$((n_main + c))
done
n_head=$(grep -rh "@testset" test/ | grep -c "@testset")
echo "@testset count: main=$n_main, HEAD=$n_head"
test "$n_head" -gt "$n_main" || { echo "FAIL: no new tests added across the entire refactor — suspicious"; exit 1; }
```

- [ ] **Step 3: Final LOC verification**

```bash
wc -l src/*.jl
```

Expected: total ≤ 8,200 (per continuation spec §3.2). If not, return to 7e.1 for more aggressive cleanup, or revise the expected target with Denis.

- [ ] **Step 4: Update memory entries**

In `/home/denis.linux/.claude/projects/-home-denis-linux--julia-dev-EnzymeRates/memory/MEMORY.md`, add an index entry for the refactor's completion. Replace stale entries (e.g., `project_structs_throughout_refactor.md` which described an earlier attempt, `project_mechanism_types_refactor_complete.md` which described the earlier ship of partial structs).

Use the auto-memory format from CLAUDE.md.

- [ ] **Step 5: Draft PR description**

Based on `git log main..HEAD --oneline`, draft a PR description covering:

- Motivation (spec link)
- Per-stage bullet-point summary
- Behavior changes (DSL grammar `:ES → E(S)`, structural kinetic groups, removal of `EnzymeReactionLegacy`)
- Migration path for users (the `:K1` naming preserved — point to the deferred refactor)
- Perf gates passed
- LOC delta (honest numbers — celebrate the architectural simplification)

- [ ] **Step 6: Push branch + open PR**

```bash
git push -u origin refactor-to-concrete-types-instead-of-symbols
gh pr create --title "Refactor EnzymeRates to concrete types" \
  --body "$(cat <<'EOF'
## Summary
- Unified concrete struct family (`Mechanism`/`Step`/`Species`/`Parameter`) shared between enumeration and derivation.
- Single chokepoint `name(p::Parameter, m::Mechanism)` for parameter-name rendering — no Symbol-string dispatch in src.
- DSL emits `EnzymeMechanism(Mechanism(...))` with decomposed-Species notation (`E(S)`).
- `EnzymeReactionLegacy` and opaque-form helpers removed.

## Per-stage changes
[fill in from git log]

## Behavior changes
- `@enzyme_mechanism` / `@allosteric_mechanism` grammar: `:ES` → `E(S)` for bound-metabolite forms; bare conformations (`:E`, `:Estar`) unchanged.
- `EnzymeReactionLegacy` removed; all dispatch now on the concrete `EnzymeReaction`.
- Parameter naming preserved (positional `:K1`, `:k10f`) — change deferred to a follow-up refactor.

## Perf gates passed
- `rate_equation` 0-alloc/<100ns per call: GREEN
- `loss!` runtime: no regression vs main baseline
- Compile-budget gates: GREEN (re-baselined to new EnzymeReaction trace counts)

## LOC
- main: 7,136 src LOC
- HEAD: <FILL_IN> src LOC (+<X>%)
- Per continuation spec §3.2 the LOC target is ≤8,200; reaching ≤3,600 (original goal) requires removing parameter-name indexing, deferred to a follow-up refactor.

## Test integrity (spec §2 NON-NEGOTIABLE)
- [x] `bash scripts/check_test_integrity.sh main` PASSES
- [x] @testset count grew vs main
- [x] No `@test_skip` / `@test_broken` introduced
- [x] No `@test` lines commented out
- [x] All test deletions logged in docs/superpowers/refactor-deleted-tests.md per §2.1 narrow exception

## Test plan
- [x] Full test suite green
- [x] All perf gates passing

Spec: docs/superpowers/specs/2026-05-22-concrete-types-refactor-continuation-design.md
Original spec: docs/superpowers/specs/2026-05-20-concrete-types-refactor-design.md
EOF
)"
```

- [ ] **Step 7: Tag final**

```bash
git tag refactor-complete
git push --tags
```

---

# Plan complete

Cumulative target met:
- src LOC: 7,136 → ≤8,200 (≈+15%)
- All tests preserved + adapted mechanically + new tests added for chokepoint regression
- All perf gates green
- One PR opened on `refactor-to-concrete-types-instead-of-symbols`
- Foundation laid for the deferred parameter-naming refactor: `name(p, m)` chokepoint exclusivity enforced, `Step.source_idx` isolated as presentation metadata.
