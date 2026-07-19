# Skip Redundant Regulatory-Site Merges — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stop the `_expand_merge_regulatory_sites` move from emitting the one merged child that duplicates an equation the search already holds — the all-keep merge of an all-`:OnlyA` site with an all-`:OnlyI` site.

**Architecture:** Two sites act on disjoint conformations (an all-`:OnlyA` site binds only the active state, an all-`:OnlyI` site only the inactive state). Merging them onto one site derives to the same rate equation as leaving them separate, so the all-keep merged child is redundant. Detect disjoint coverage per merged pair and drop that one assignment; keep the `:EqualAI` antagonist retags and every other pair.

**Tech Stack:** Julia, `Test` stdlib. Package sources in `src/`, tests in `test/` wired through `test/runtests.jl`.

## Global Constraints

- 92-character line limit, 4-space indentation. Match surrounding style.
- `rate_equation` stays allocation-free and sub-120ns; this change does not touch it, but the full suite's performance gate must stay green.
- Canonical Step Form is load-bearing; do not reorder steps/groups.
- All new `src` files start with two `# ABOUTME:` lines. (No new files here.)
- Names describe purpose, not history. No "new"/"old"/"improved".

## Test-run commands

- **Fast (enumeration file only, ~1–2 min):**
  ```bash
  julia --project -e 'using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_mechanism_enumeration.jl")'
  ```
- **Full suite (final gate, several min):**
  ```bash
  julia --project -e 'using Pkg; Pkg.test()'
  ```

---

### Task 1: `_site_active_states` helper

**Files:**
- Modify: `src/mechanism_enumeration.jl` (add helper just above `_expand_merge_regulatory_sites`, currently line 2081)
- Test: `test/test_mechanism_enumeration.jl` (new testset just above the `_expand_merge_regulatory_sites` testset, currently line 4065)

**Interfaces:**
- Produces: `_site_active_states(site::RegulatorySite) -> Set{Symbol}` returning a subset of `{:active, :inactive}`.

- [ ] **Step 1: Write the failing test.** Insert above line 4065 (`# ─── _expand_merge_regulatory_sites ───`):

```julia
@testset "_site_active_states" begin
    mk(states) = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(Symbol("R", i)) for i in eachindex(states)],
        4, collect(Symbol, states))
    @test EnzymeRates._site_active_states(mk([:OnlyA])) == Set([:active])
    @test EnzymeRates._site_active_states(mk([:OnlyI])) == Set([:inactive])
    @test EnzymeRates._site_active_states(mk([:EqualAI])) == Set([:active, :inactive])
    @test EnzymeRates._site_active_states(mk([:NonequalAI])) == Set([:active, :inactive])
    @test EnzymeRates._site_active_states(mk([:OnlyA, :OnlyA])) == Set([:active])
    @test EnzymeRates._site_active_states(mk([:OnlyA, :OnlyI])) ==
          Set([:active, :inactive])
end
```

- [ ] **Step 2: Run to verify it fails.** Run the fast command. Expected: `UndefVarError: _site_active_states` (function not defined).

- [ ] **Step 3: Implement.** Insert above line 2081 (`function _expand_merge_regulatory_sites`):

```julia
"""
    _site_active_states(site::RegulatorySite) -> Set{Symbol}

The conformations a regulatory site's ligands act on: `:active` for any ligand
binding the active state (`:OnlyA`/`:EqualAI`/`:NonequalAI`), `:inactive` for
any binding the inactive state (`:OnlyI`/`:EqualAI`/`:NonequalAI`). Two sites
with disjoint active states — an all-`:OnlyA` site and an all-`:OnlyI` site —
merge to a rate equation identical to keeping them separate, so that all-keep
merge is redundant.
"""
function _site_active_states(site::RegulatorySite)
    active = Set{Symbol}()
    for st in allo_states(site)
        st in (:OnlyA, :EqualAI, :NonequalAI) && push!(active, :active)
        st in (:OnlyI, :EqualAI, :NonequalAI) && push!(active, :inactive)
    end
    active
end
```

- [ ] **Step 4: Run to verify it passes.** Run the fast command. Expected: the `_site_active_states` testset passes.

- [ ] **Step 5: Commit.**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add _site_active_states: conformations a regulatory site acts on"
```

---

### Task 2: `drop_all_keep` on `_merged_site_state_assignments`

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_merged_site_state_assignments`, currently ~line 2117)
- Test: `test/test_mechanism_enumeration.jl` (new testset directly after the `_site_active_states` testset from Task 1)

**Interfaces:**
- Consumes: nothing from Task 1.
- Produces: `_merged_site_state_assignments(base_states::Vector{Symbol}; drop_all_keep::Bool=false) -> Vector{Vector{Symbol}}`. Default behavior unchanged; `drop_all_keep=true` omits the `copy(base_states)` all-keep entry, retaining the `:EqualAI` retags.

- [ ] **Step 1: Write the failing test.** Insert after the `_site_active_states` testset:

```julia
@testset "_merged_site_state_assignments drop_all_keep" begin
    base = [:OnlyA, :OnlyI]
    keep = EnzymeRates._merged_site_state_assignments(base)
    @test [:OnlyA, :OnlyI] in keep          # all-keep present by default
    @test [:EqualAI, :OnlyI] in keep
    @test [:OnlyA, :EqualAI] in keep
    @test length(keep) == 3
    dropped = EnzymeRates._merged_site_state_assignments(base; drop_all_keep=true)
    @test !([:OnlyA, :OnlyI] in dropped)    # all-keep omitted
    @test [:EqualAI, :OnlyI] in dropped     # antagonist retags retained
    @test [:OnlyA, :EqualAI] in dropped
    @test length(dropped) == 2
end
```

- [ ] **Step 2: Run to verify it fails.** Run the fast command. Expected: a `MethodError` / unknown-keyword error on the `drop_all_keep=true` call.

- [ ] **Step 3: Implement.** Replace the body of `_merged_site_state_assignments` (currently line ~2117). Also update its docstring's first line to mention the keyword. New version:

```julia
"""
    _merged_site_state_assignments(base_states::Vector{Symbol};
                                   drop_all_keep=false) -> Vector{Vector{Symbol}}

Δ0-valid allo-state assignments for a merged site's ligands: the all-keep
assignment (omitted when `drop_all_keep`), plus each assignment retagging
exactly one non-`:EqualAI` ligand to `:EqualAI`. The all-`:EqualAI` result is
dropped. `drop_all_keep` omits the all-keep entry for a redundant merge (see
`_site_active_states`) while keeping the antagonist retags.
"""
function _merged_site_state_assignments(base_states::Vector{Symbol};
                                        drop_all_keep::Bool=false)
    assignments = Vector{Symbol}[]
    drop_all_keep || push!(assignments, copy(base_states))
    for i in eachindex(base_states)
        base_states[i] == :EqualAI && continue
        retagged = copy(base_states)
        retagged[i] = :EqualAI
        all(==(:EqualAI), retagged) && continue
        push!(assignments, retagged)
    end
    assignments
end
```

- [ ] **Step 4: Run to verify it passes.** Run the fast command. Expected: the `drop_all_keep` testset passes.

- [ ] **Step 5: Commit.**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add drop_all_keep to _merged_site_state_assignments"
```

---

### Task 3: Wire the redundancy skip into `_expand_merge_regulatory_sites`

This is the red-green driver: first update the two existing tests to assert the new behavior (they fail against the unchanged move), then implement the skip.

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_expand_merge_regulatory_sites`, line 2081–2099)
- Test: `test/test_mechanism_enumeration.jl` (update the testsets at lines ~4096, ~4122, and ~4783)

**Interfaces:**
- Consumes: `_site_active_states` (Task 1), `_merged_site_state_assignments(...; drop_all_keep)` (Task 2).

- [ ] **Step 1: Update the existing move tests to the new behavior.**

  (a) Replace the testset currently at lines 4096–4104 (`@testset "co-binding + both antagonist forms appear"`) with:

```julia
    @testset "disjoint OnlyA/OnlyI co-binding skipped; antagonists kept" begin
        # A (:OnlyA) acts only on the active state, I (:OnlyI) only on the
        # inactive one, so co-binding them on one site derives to the same
        # equation as separate sites — that all-keep merge is redundant and
        # skipped. The two antagonist retags stay (each :EqualAI ligand now
        # acts on both states, a genuinely distinct mechanism).
        @test all(c -> single_site_names(c) == Set([:A, :I]), children)
        states = Set((merged_state(c, :A), merged_state(c, :I)) for c in children)
        @test !((:OnlyA, :OnlyI) in states)  # redundant co-binding skipped
        @test (:EqualAI, :OnlyI) in states   # activator → antagonist
        @test (:OnlyA, :EqualAI) in states   # inhibitor → antagonist
        @test !((:EqualAI, :EqualAI) in states)  # all-EqualAI dropped
        @test length(children) == 2
    end
```

  (b) In the testset `@testset "reg_type filter keeps all when activator/inhibitor differ in type"` (currently lines 4122–4129), change the final assertion from `@test length(kept) == 3` to:

```julia
        @test length(kept) == 2
```

  (c) In the testset `@testset "Merge move wired for two-site allosteric; no-op for Mechanism"` (currently line 4758), update the assertion at line ~4785 and its comment. Change:

```julia
        # The three Δ0 merge children (co-binding + two antagonist forms)
        # survive the reg_type filter and appear in the output.
        @test count(is_merged, children) == 3
```

  to:

```julia
        # The two Δ0 antagonist merge children survive the reg_type filter and
        # appear in the output; the redundant OnlyA/OnlyI co-binding is skipped.
        @test count(is_merged, children) == 2
```

- [ ] **Step 2: Run to verify the updated tests fail.** Run the fast command. Expected: the three updated assertions fail (the move still emits 3 children, so `length(children) == 2` and `count(is_merged, children) == 2` fail, and `!((:OnlyA, :OnlyI) in states)` fails).

- [ ] **Step 3: Implement the skip.** In `_expand_merge_regulatory_sites`, replace the single line

```julia
        for states in _merged_site_state_assignments(base_states)
```

  with:

```julia
        redundant = isempty(intersect(_site_active_states(sites[i]),
                                       _site_active_states(sites[j])))
        for states in _merged_site_state_assignments(base_states;
                                                     drop_all_keep=redundant)
```

- [ ] **Step 4: Run to verify it passes.** Run the fast command. Expected: the full `_expand_merge_regulatory_sites` testset and the `Merge move wired…` testset pass. In particular the same-direction (`:OnlyA` + `:OnlyA`) and 3-way merge testsets still pass (those pairs are not disjoint, so their all-keep is retained).

- [ ] **Step 5: Commit.**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Skip the redundant OnlyA/OnlyI all-keep regulatory-site merge"
```

---

### Task 4: Derivation-equivalence premise guard

Lock the claim the skip rests on: the skipped merged child derives to the same rate equation as the separate-site form. This test passes before and after the fix; it guards against a future change to the derivation that would make the skip unsound.

**Files:**
- Test: `test/test_mechanism_enumeration.jl` (new testset directly after the `_expand_merge_regulatory_sites` testset, i.e. after its closing `end` near line 4230)

- [ ] **Step 1: Write the guard test.**

```julia
@testset "OnlyA+OnlyI merge is derivation-redundant (premise guard)" begin
    # The reason the all-keep OnlyA/OnlyI merge is skipped: merged onto one
    # site it evaluates to the same rate as on separate sites. Confirm
    # numerically over random parameters and concentrations.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A::Activator, I::Inhibitor
        oligomeric_state: 4
    end
    base = first(EnzymeRates.init_mechanisms(rxn))
    cat = Symbol[:OnlyA for _ in 1:length(EnzymeRates.steps(base))]
    mkm(sites) = EnzymeRates.AllostericMechanism(
        rxn, copy(EnzymeRates.steps(base)), cat, 4, sites)
    A = EnzymeRates.AllostericRegulator(:A)
    I = EnzymeRates.AllostericRegulator(:I)
    separate = mkm([EnzymeRates.RegulatorySite([A], 4, [:OnlyA]),
                    EnzymeRates.RegulatorySite([I], 4, [:OnlyI])])
    merged = mkm([EnzymeRates.RegulatorySite([A, I], 4, [:OnlyA, :OnlyI])])
    fps = EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(separate))
    fpm = EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(merged))
    @test Set(fps) == Set(fpm)
    Random.seed!(42)
    for _ in 1:25
        vals = Dict(p => exp(randn()) for p in union(fps, fpm))
        ps = (; (p => vals[p] for p in fps)..., Keq=100.0, E_total=1.0)
        pm = (; (p => vals[p] for p in fpm)..., Keq=100.0, E_total=1.0)
        c = (; S=exp(randn()), P=exp(randn()), A=exp(randn()), I=exp(randn()))
        vs = real(EnzymeRates.rate_equation(separate, c, ps))
        vm = real(EnzymeRates.rate_equation(merged, c, pm))
        @test isapprox(vs, vm; rtol=1e-9)
    end
end
```

- [ ] **Step 2: Run to verify it passes.** Run the fast command. Expected: the guard testset passes (25 numeric checks plus the param-set check).

  If `rate_equation` on this uni-uni construction errors (e.g. a missing conc field), read the error: add every regulator/reactant name it names to the `c` NamedTuple. Do not weaken the `rtol`; the equations are equal to machine precision.

- [ ] **Step 3: Commit.**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Guard: OnlyA/OnlyI merge derives identically to separate sites"
```

---

### Task 5: Full-suite regression gate

**Files:** none (verification only).

- [ ] **Step 1: Run the full suite.**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass, including `test_mechanism_enumeration.jl`, the allosteric golden/ground-truth files, `test_identify_rate_equation.jl`, and the `rate_equation` performance gate in `test_rate_eq_derivation.jl`.

- [ ] **Step 2: Adjudicate any failure.** A failure here should be an integration test that asserts a specific allosteric enumeration count for a reaction carrying both an `:OnlyA` and an `:OnlyI` regulator. If one fails:
  1. Confirm the delta is exactly the skipped redundant merges (the count drops, and the dropped mechanisms are all-keep OnlyA/OnlyI co-bindings whose equation is already present via the separate-site form).
  2. Update the expected count and add a one-line comment explaining the redundant-merge skip.
  3. If a failure is *not* explained by the redundant-merge skip, stop — it is a real regression; do not adjust the test.

- [ ] **Step 3: Commit any test-count updates.**

```bash
git add test/
git commit -m "Update allosteric enumeration counts for skipped redundant merges"
```

---

## Self-review

- **Spec coverage:** `_site_active_states` (Task 1) ↔ spec change 1; `drop_all_keep` (Task 2) ↔ spec change 2; wiring + move-behavior test (Task 3) ↔ spec change 3 and test 1; derivation guard (Task 4) ↔ test 2; `_site_active_states` unit ↔ test 3; full suite (Task 5) ↔ test 4. All covered.
- **Placeholders:** none — every step has exact code and commands.
- **Type consistency:** `_site_active_states` returns `Set{Symbol}` of `{:active,:inactive}`, consumed by `isempty(intersect(...))` in Task 3. `_merged_site_state_assignments(...; drop_all_keep::Bool)` defined in Task 2, called in Task 3. Names consistent throughout.
