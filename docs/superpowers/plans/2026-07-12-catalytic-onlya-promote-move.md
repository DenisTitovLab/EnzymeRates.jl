# Catalytic-`:OnlyA` Promote Move Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make every multi-`:OnlyA` catalytic-state assignment reachable by the mechanism search, via one new Δ0 expansion move.

**Architecture:** Add `_expand_promote_catalytic_to_onlya(am::AllostericMechanism)` — for each `:EqualAI` catalytic group, emit a variant with it set to `:OnlyA` — and wire it into the beam's per-parent expansion (`_add_expansions_mech!`), not the seed build. Repeated application from single-`:OnlyA` seeds reaches every `:OnlyA` subset; the beam prunes as with every other structural axis.

**Tech Stack:** Julia, EnzymeRates.jl. Tests run under `julia --project`.

## Global Constraints

- **`rate_equation` performance is non-negotiable:** every mechanism in `MECHANISM_TEST_SPECS` must keep `rate_equation` at `allocs == 0` and `t < 120e-9` (`test/test_rate_eq_derivation.jl`). If a newly-reachable multi-`:OnlyA` mechanism violates this, **STOP and flag Denis** — do not work around it.
- **TDD:** failing test first, confirm it fails, minimal implementation, confirm it passes, commit.
- **Style:** 92-char line limit, 4-space indent. Match surrounding code. All source files start with two `# ABOUTME:` lines (already present in the files touched).
- **Canonical Step Form is load-bearing:** the `AllostericMechanism` constructor canonicalizes group/tag order. Tests must assert by *structure* (tag multisets, per-group tags), never by raw index position.
- **Commit trailers** on every commit:
  ```
  Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01VFCCyUsAM2fdaXLq5Pdtw1
  ```
- Branch: `catalytic-onlya-promote-move` (already checked out; spec committed at `ad7cc20`).

---

### Task 1: The `_expand_promote_catalytic_to_onlya` move + move-level tests

**Files:**
- Modify: `src/mechanism_enumeration.jl` — add the move next to `_expand_change_allo_state` (≈ line 1959).
- Test: `test/test_mechanism_enumeration.jl` — new `@testset` after the `_expand_change_allo_state` testset (≈ line 3772+).

**Interfaces:**
- Produces: `_expand_promote_catalytic_to_onlya(am::AllostericMechanism) -> Vector{AllostericMechanism}` and `_expand_promote_catalytic_to_onlya(::Mechanism) -> AllostericMechanism[]`.
- Consumes: existing `cat_allo_states`, `_with_cat_allo_states` (`src/types.jl:701`), `rep_step`, `is_iso`, `compile_mechanism`, `fitted_params`, `rate_equation_string`, the `@allosteric_mechanism` DSL, `AllostericMechanism`.

- [ ] **Step 1: Write the failing tests.**

Add to `test/test_mechanism_enumeration.jl`:

```julia
# ─── _expand_promote_catalytic_to_onlya ────────────────────────────────
@testset "_expand_promote_catalytic_to_onlya" begin

    # SEED: bi-uni ordered, one substrate group :OnlyA, the other catalytic
    # groups :EqualAI (one binding B, one iso catalytic, one binding P).
    biuni_seed() = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + A ⇌ E(A)              :: OnlyA
            E(A) + B ⇌ E(A, B)       :: EqualAI
            E(A, B) <--> E(P)        :: EqualAI
            E + P ⇌ E(P)             :: EqualAI
        end
    end)

    @testset "one variant per :EqualAI catalytic group, each Δ0-shaped" begin
        am = biuni_seed()
        EnzymeRates._assert_mechanism_invariants(am)
        n_eq = count(==(:EqualAI), am.cat_allo_states)
        result = EnzymeRates._expand_promote_catalytic_to_onlya(am)

        @test length(result) == n_eq          # 3 here
        for child in result
            @test child isa EnzymeRates.AllostericMechanism
            EnzymeRates._assert_mechanism_invariants(child)
            # exactly one :EqualAI became :OnlyA; nothing else changed shape
            @test count(==(:OnlyA), child.cat_allo_states) ==
                  count(==(:OnlyA), am.cat_allo_states) + 1
            @test count(==(:EqualAI), child.cat_allo_states) == n_eq - 1
            @test child.regulatory_sites == am.regulatory_sites
            @test EnzymeRates.catalytic_multiplicity(child) ==
                  EnzymeRates.catalytic_multiplicity(am)
        end
    end

    @testset "covers the iso/catalytic group, not just binding groups" begin
        am = biuni_seed()
        iso_g = only(g for g in EnzymeRates.kinetic_groups(am)
                     if EnzymeRates.is_iso(EnzymeRates.rep_step(am, g)) &&
                        am.cat_allo_states[g] == :EqualAI)
        result = EnzymeRates._expand_promote_catalytic_to_onlya(am)
        # some child has the iso group promoted to :OnlyA
        @test any(c -> c.cat_allo_states[iso_g] == :OnlyA, result)
    end

    @testset "no-op when no catalytic group is :EqualAI" begin
        am = biuni_seed()
        all_onlya = EnzymeRates._with_cat_allo_states(
            am, fill(:OnlyA, length(am.cat_allo_states)))
        @test isempty(EnzymeRates._expand_promote_catalytic_to_onlya(all_onlya))
    end

    @testset "no-op on a non-allosteric Mechanism" begin
        m = first(EnzymeRates.init_mechanisms(@enzyme_reaction begin
            substrates: S; products: P; oligomeric_state: 2
        end))
        @test isempty(EnzymeRates._expand_promote_catalytic_to_onlya(m))
    end

    @testset "order-independent: promote i then j == j then i" begin
        am = biuni_seed()
        # promote any child again; two distinct first-promotions that then
        # reach the same 3-:OnlyA structure must be ==.
        kids = EnzymeRates._expand_promote_catalytic_to_onlya(am)
        grandkids = reduce(vcat,
            (EnzymeRates._expand_promote_catalytic_to_onlya(k) for k in kids))
        # dedup by structural hash; two orders to the same OnlyA-set collapse
        @test length(unique(grandkids)) < length(grandkids)
    end

    @testset "Δ0: promoted children keep the parent's fitted-param count" begin
        am = biuni_seed()
        base = length(EnzymeRates.fitted_params(
            EnzymeRates.compile_mechanism(am)))
        for child in EnzymeRates._expand_promote_catalytic_to_onlya(am)
            @test length(EnzymeRates.fitted_params(
                EnzymeRates.compile_mechanism(child))) == base
        end
    end

    @testset "distinct: every promotion changes the rate equation" begin
        am = biuni_seed()
        parent_eq = EnzymeRates.rate_equation_string(am)
        for child in EnzymeRates._expand_promote_catalytic_to_onlya(am)
            @test EnzymeRates.rate_equation_string(child) != parent_eq
        end
    end
end
```

- [ ] **Step 2: Run the tests to verify they fail.**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep -A3 promote_catalytic` (or run just the file — see note below).
Expected: FAIL with `UndefVarError: _expand_promote_catalytic_to_onlya` (function not yet defined).

Note: running the whole suite is slow. To iterate on just this file:
`julia --project -e 'include("test/runtests.jl")'` runs everything; for a faster loop, temporarily `include` the single test file inside a scratch script that first does `using EnzymeRates, Test` plus the file's local helpers. Do NOT commit any scratch runner.

- [ ] **Step 3: Implement the move.**

In `src/mechanism_enumeration.jl`, immediately after `_expand_change_allo_state(::Mechanism)` (≈ line 1967), add:

```julia
"""
    _expand_promote_catalytic_to_onlya(am::AllostericMechanism)
        → Vector{AllostericMechanism}

Δ0 catalytic-state move. For each catalytic kinetic group tagged `:EqualAI`,
emit one variant with that group set to `:OnlyA` — binding (K-type) and
iso/catalytic (V-type) groups alike. On an already-allosteric mechanism `L`
is already observable, so promoting any group is distinguishable and never
degenerate. The catalytic steps, multiplicity, regulatory sites, and every
other tag pass through unchanged.
"""
function _expand_promote_catalytic_to_onlya(am::AllostericMechanism)
    results = AllostericMechanism[]
    for g in 1:length(cat_allo_states(am))
        cat_allo_states(am)[g] == :EqualAI || continue
        new_states = copy(cat_allo_states(am))
        new_states[g] = :OnlyA
        push!(results, _with_cat_allo_states(am, new_states))
    end
    results
end

"""
    _expand_promote_catalytic_to_onlya(m::Mechanism)
        → Vector{AllostericMechanism}

Non-allosteric input: no-op; this move only elaborates allosteric states.
"""
_expand_promote_catalytic_to_onlya(::Mechanism) =
    AllostericMechanism[]
```

- [ ] **Step 4: Run the tests to verify they pass.**

Run the new testset. Expected: PASS (all subtests green).

- [ ] **Step 5: Commit.**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit  # message: "Add Δ0 _expand_promote_catalytic_to_onlya move + tests" + trailers
```

---

### Task 2: Wire the move into the beam + reachability integration test

**Files:**
- Modify: `src/mechanism_enumeration.jl` — `_add_expansions_mech!` (≈ line 2131-2142) and the `expand_mechanisms` docstring (≈ line 2109-2116).
- Test: `test/test_mechanism_enumeration.jl` — reachability testset (place near the other `expand_mechanisms` integration tests).

**Interfaces:**
- Consumes: `_expand_promote_catalytic_to_onlya` (Task 1), `expand_mechanisms`, `init_mechanisms`, `cat_allo_states`, `rep_step`, `bound_metabolite`, `is_iso`.

- [ ] **Step 1: Write the failing reachability test.**

Add to `test/test_mechanism_enumeration.jl`:

```julia
@testset "expand_mechanisms reaches ≥2 distinct-metabolite catalytic :OnlyA" begin
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products:   P[C1N1]
        allosteric_regulators: R
        oligomeric_state: 2
    end

    # #distinct metabolites bound by an :OnlyA catalytic group (iso → skip)
    onlya_mets(am) = Set(EnzymeRates.name(EnzymeRates.bound_metabolite(
                            EnzymeRates.rep_step(am, g)))
        for g in EnzymeRates.kinetic_groups(am)
        if EnzymeRates.cat_allo_states(am)[g] === :OnlyA &&
           !EnzymeRates.is_iso(EnzymeRates.rep_step(am, g)))

    seen = Set{UInt64}()
    frontier = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[]
    for m in EnzymeRates.init_mechanisms(rxn)
        h = hash(m); h in seen || (push!(seen, h); push!(frontier, m))
    end
    maxdistinct = 0
    gen = 0
    while !isempty(frontier) && gen < 14
        nextf = eltype(frontier)[]
        for c in EnzymeRates.expand_mechanisms(frontier, rxn)
            h = hash(c); h in seen || (push!(seen, h); push!(nextf, c))
            c isa EnzymeRates.AllostericMechanism &&
                (maxdistinct = max(maxdistinct, length(onlya_mets(c))))
        end
        frontier = nextf; gen += 1
    end

    @test maxdistinct >= 2
end
```

- [ ] **Step 2: Run to verify it fails.**

Run the testset. Expected: FAIL — `maxdistinct` is 1 (move not wired into the beam yet).

- [ ] **Step 3: Wire the move in.**

In `src/mechanism_enumeration.jl`, in `_add_expansions_mech!`, add one line after the `_expand_change_allo_state` append (≈ line 2140):

```julia
    append!(result, _expand_change_allo_state(m))
    append!(result, _expand_promote_catalytic_to_onlya(m))
    append!(result, _expand_merge_regulatory_sites(m))
```

Update the `expand_mechanisms` docstring move list (≈ line 2112) to include the new move, e.g. change `…, change allo state, merge regulatory sites)` to `…, change allo state, promote catalytic to OnlyA, merge regulatory sites)`.

- [ ] **Step 4: Run to verify it passes.**

Run the testset. Expected: PASS (`maxdistinct >= 2`).

- [ ] **Step 5: Commit.**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit  # "Wire promote-catalytic-to-OnlyA into expand_mechanisms" + trailers
```

---

### Task 3: Performance + derivation golden coverage for multi-`:OnlyA`

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` — add a `MechanismTestSpec` inside `build_mechanism_test_specs()` (before `return specs`, ≈ line 2505).

**Interfaces:**
- Consumes: `@allosteric_mechanism`, `MechanismTestSpec`, `MECHANISM_TEST_SPECS`, the derivation + performance harness in `test/test_rate_eq_derivation.jl`.

- [ ] **Step 1: Add a multi-`:OnlyA` spec with placeholder golden counts.**

Insert into `build_mechanism_test_specs()`:

```julia
    # Two catalytic :OnlyA binding groups — the multi-:OnlyA family the
    # promote move makes reachable. Guards derivation + the rate_equation
    # 0-alloc / sub-120ns contract for that family. RE bindings, no ODE.
    push!(specs, MechanismTestSpec(
        name="multi-OnlyA bi-uni (A,B both OnlyA)",
        mechanism=(@allosteric_mechanism begin
            substrates: A, B
            products: P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + A ⇌ E(A)          :: OnlyA
                E(A) + B ⇌ E(A, B)    :: OnlyA
                E(A, B) <--> E(P)     :: EqualAI
                E + P ⇌ E(P)          :: EqualAI
            end
        end),
        metabolite_names=[:A, :B, :P],
        expected_n_states=0, expected_n_steps=0, expected_n_metabolites=3,
        expected_n_haldane_constraints=1, expected_n_mirror_constraints=0,
        expected_n_wegscheider_constraints=0, expected_n_independent_params=0,
        run_ode_test=false))
```

- [ ] **Step 2: Derive it once to read the real structural counts.**

Run a scratch script (do NOT commit it):

```julia
using EnzymeRates
const ER = EnzymeRates
em = @allosteric_mechanism begin
    substrates: A, B
    products: P
    catalytic_multiplicity: 2
    catalytic_steps: begin
        E + A ⇌ E(A)          :: OnlyA
        E(A) + B ⇌ E(A, B)    :: OnlyA
        E(A, B) <--> E(P)     :: EqualAI
        E + P ⇌ E(P)          :: EqualAI
    end
end
am = ER.AllostericMechanism(em)
cem = ER.compile_mechanism(am)
println("fitted_params: ", ER.fitted_params(cem))
println("rate_eq: ", ER.rate_equation_string(am))
```

Note the `@allosteric_mechanism` here lives in the test harness (`test/mechanism_definitions_for_test_enzyme_derivation.jl`), which `include`s the DSL; run the scratch from within `julia --project` after `include`-ing that file, or reuse the exact mechanism already added to the specs. Record `n_independent_params` from `length(fitted_params)`.

For `expected_n_states` / `expected_n_steps`: run the derivation test for this spec name and read the reported actual values from the failure message, then paste them back into the spec. `expected_n_haldane_constraints` MUST be `1` (one reversible reaction); if the derivation reports otherwise, STOP — the mechanism is malformed.

- [ ] **Step 3: Fill in the golden counts and re-run the spec's derivation test.**

Replace the `0` placeholders with the observed values. Run the derivation + performance tests for this spec.
Expected: PASS, including `allocs == 0` and `t < 120e-9`.

**If the performance assertion fails for this mechanism, STOP and flag Denis** — per the non-negotiable perf contract, do not adjust the bound or work around it.

- [ ] **Step 4: Commit.**

```bash
git add test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit  # "Add multi-OnlyA golden + perf-contract coverage" + trailers
```

---

### Final verification

- [ ] **Run the full suite:** `julia --project -e 'using Pkg; Pkg.test()'`
- [ ] Confirm **all green** — especially `test_rate_equation_performance` (0-alloc / sub-120ns) and the allosteric derivation/collapse/golden testsets. Enumeration count tests elsewhere may shift because the beam now reaches more mechanisms; if any *count* assertion changes, verify the new count reflects the reachable multi-`:OnlyA` mechanisms (expected) rather than a regression, and update the asserted count with a comment noting the promote move as the cause.
- [ ] If everything is green, the branch is ready for Denis's review.
