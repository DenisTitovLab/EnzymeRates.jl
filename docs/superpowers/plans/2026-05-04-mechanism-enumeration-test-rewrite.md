# Mechanism Enumeration Test Rewrite Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rewrite `test/test_mechanism_enumeration.jl` in pipeline-execution order, replacing brittle `init_mechanisms |> first` patterns with literal `@enzyme_mechanism` / `@allosteric_mechanism` seeds, and applying a 7-item per-move checklist with independent-derivation comments for every count, delta, and equivalence-set entry.

**Architecture:** Test-only refactor. The single new piece of code is `allosteric_spec_from_mechanism(m, rxn)`, a round-trip helper symmetric to the existing `mechanism_spec_from_mechanism`, lives at the top of the test file. The file is reorganized into 7 sections (0 infrastructure → 6 integration) reflecting pipeline depth. No `src/` changes are planned; if a test surfaces a bug, follow the bug-handling protocol in §6 of the design (`docs/superpowers/specs/2026-05-04-mechanism-enumeration-test-rewrite-design.md`).

**Tech Stack:** Julia 1.x, Test stdlib, `EnzymeRates.jl` package's enumeration pipeline (`src/mechanism_enumeration.jl`).

---

## Reference

- **Design doc:** `docs/superpowers/specs/2026-05-04-mechanism-enumeration-test-rewrite-design.md` — read sections 3 (checklist) and 6 (bug-handling protocol) before starting any task.
- **Pipeline source:** `src/mechanism_enumeration.jl` — function definitions, expansion-move semantics.
- **Existing tests being rewritten:** `test/test_mechanism_enumeration.jl`.
- **DSL macros:** `src/dsl.jl` — `@enzyme_reaction`, `@enzyme_mechanism`, `@allosteric_mechanism`.

## File Structure

The plan modifies a single test file plus produces ~13 commits. The end-state layout:

```
test/test_mechanism_enumeration.jl          (rewritten in place, ~3000 lines)
├── 0. Test infrastructure                   (helpers + shared @enzyme_reaction defs)
├── 1. Support functions                     (no spec input; mostly preserved)
├── 2. Initialization                        (compile_mechanism + init_mechanisms)
├── 3. Base-spec expansion moves             (re_to_ss, split, add_dead_end_reg)
├── 4. Allosteric expansion moves            (to_allosteric, add_allo_reg, change_allo_state)
├── 5. Composition                           (dedup!, expand_mechanisms)
└── 6. Integration                           (enumerate_all)
```

Three testsets that don't belong to enumeration may be moved out (Task 13, optional and confirmation-gated):

```
test/test_types.jl                           (+ AllostericEnzymeMechanism TR equivalence)
test/test_dsl.jl                             (+ test reaction atom balance)
test/test_rate_eq_derivation.jl              (+ Tagged groups exclude T-state params)
```

## Common verification

After every test edit, the engineer runs:

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: green. If red, follow the bug-handling protocol in design §6 — surface the failure with the structured report, decide test-vs-code, fix the right side, never weaken or skip the test.

The test command is slow (cold precompilation each invocation per CLAUDE.md). Run it once at end of each task before commit; do NOT run after every individual seed addition unless debugging.

## Important DSL invariants the engineer must know

- `@enzyme_mechanism` step arrows: `⇌` = RE (equilibrium), `<-->` = SS (steady state).
- `@enzyme_mechanism` kinetic_group assignment: each top-level step gets a fresh `kinetic_group` integer in source order (1, 2, 3, …). A parenthesized step group `(stepA; stepB)` shares one group.
- `@allosteric_mechanism` requires per-step or per-step-group `:: <Tag>` annotations within `site(:catalytic, N)`. Allowed tags: `:OnlyR`, `:OnlyT`, `:EqualRT`, `:NonequalRT`. Allosteric regulators require `name::Tag` per entry.
- `mechanism_spec_from_mechanism` (existing) and `allosteric_spec_from_mechanism` (Task 1 adds it) round-trip a compiled mechanism back to a spec. Both round-trips must be validated at every call site via `=== m_seed`.
- Compiled `EnzymeMechanism` / `AllostericEnzymeMechanism` are singleton types; `===` and `==` are equivalent. Step ORDER and `kinetic_group` numbering are part of type identity. So when writing expected-mechanism literals for equivalence-style assertions, the step order and grouping in the literal must match what the move produces.
- `_expand_re_to_ss`, `_expand_split_kinetic_group`, and other moves preserve step order from the input spec; they do NOT canonicalize. `dedup!` is the only function that canonicalizes.
- Helper-seeded specs have `n_fit_params_estimate == length(fitted_params(m))` (the *exact* fitted count), while init-seeded specs have `n_fit_params_estimate >= length(fitted_params(m))` (an *upper-bound estimate* that can be strictly larger when mirror cycles exist). Move-delta assertions (`r.n_fit_params_estimate == spec.n_fit_params_estimate + delta`) are correct under either baseline because deltas are baseline-independent — each move computes delta from tag / RE-vs-SS / multi-vs-singleton properties of the affected group. The upper-bound invariant is exercised separately in Task 4 (`init_mechanisms` testset, init-seeded); helper-seeded specs satisfy it trivially via equality.

---

## Task 1: Add `allosteric_spec_from_mechanism` helper

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (add helper near top, just below the existing `mechanism_spec_from_mechanism` at line 8)

**Pre-existing tests this absorbs:** none. Pure addition.

- [ ] **Step 1.1: Add the helper function**

Add after the existing `mechanism_spec_from_mechanism` definition (line 20):

```julia
# Helper: round-trip a compiled AllostericEnzymeMechanism back to an
# AllostericMechanismSpec. Symmetric to mechanism_spec_from_mechanism.
# AllostericMechanismSpec uses dense Dict storage — every kinetic group and
# every regulator ligand has an explicit entry, so this is a pure pass-through.
function allosteric_spec_from_mechanism(
    m::AllostericEnzymeMechanism,
    @nospecialize(rxn::EnzymeReaction))
    cm = EnzymeRates.catalytic_mechanism(m)
    base_spec = mechanism_spec_from_mechanism(cm, rxn)

    cat_n = EnzymeRates.catalytic_multiplicity(m)

    n_groups = length(unique(s.kinetic_group for s in base_spec.steps))
    group_tags = Dict{Int, Symbol}()
    for g in 1:n_groups
        group_tags[g] = EnzymeRates.cat_allo_state(m, g)
    end

    n_reg_sites = length(EnzymeRates.regulatory_sites(m))
    reg_sites = Vector{Symbol}[]
    multiplicities = Int[]
    reg_ligand_tags = Dict{Symbol, Symbol}()
    for i in 1:n_reg_sites
        ligs = collect(EnzymeRates.regulatory_site_ligands(m, i))
        push!(reg_sites, ligs)
        push!(multiplicities,
            EnzymeRates.regulatory_site_multiplicity(m, i))
        for lig in ligs
            reg_ligand_tags[lig] = EnzymeRates.reg_allo_state(m, i, lig)
        end
    end

    AllostericMechanismSpec(
        base_spec, cat_n, reg_sites, multiplicities,
        group_tags, reg_ligand_tags,
        length(EnzymeRates.fitted_params(m)))
end
```

- [ ] **Step 1.2: Verify `regulatory_sites` accessor exists, or use the right one**

Run:

```bash
grep -n "^regulatory_sites\b\|function regulatory_sites" /home/denis.linux/.julia/dev/EnzymeRates/src/types.jl
```

Expected: a definition shows up. If the accessor is named differently (e.g., `n_regulatory_sites`), update the `n_reg_sites = …` line in the helper to use the correct API. Names to check as fallbacks: `n_reg_sites(m)`, `regulatory_sites(m)` returning a tuple to take `length` of.

If neither exists, look at `AllostericEnzymeMechanism{CM,CS,RS}` — `RS` is the type-parameter tuple of reg sites. Use `length(RS)` via parametric dispatch in the helper, or extract `RS` via `typeof(m).parameters[3]`.

- [ ] **Step 1.3: Add round-trip validation testset for the helper**

Add immediately after the helper:

```julia
@testset "allosteric_spec_from_mechanism round-trip" begin
    # K-type uni-uni: catalytic 2-mer, all bindings :EqualRT, iso :EqualRT,
    # no regulators. Round-trip must be lossless: spec → AllostericEnzymeMechanism
    # rebuilds to the same singleton type as the macro produced.
    m1 = @allosteric_mechanism begin
        substrates: S
        products: P
        site(:catalytic, 2): begin
            steps: begin
                E + S ⇌ E_S       :: EqualRT
                E + P ⇌ E_P       :: EqualRT
                E_S <--> E_P      :: EqualRT
            end
        end
    end
    rxn1 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        oligomeric_state: 2
    end
    spec1 = allosteric_spec_from_mechanism(m1, rxn1)
    @test AllostericEnzymeMechanism(spec1) === m1

    # Mixed group tags: one :OnlyR, one :EqualRT, one :NonequalRT.
    # Dense storage: every group has an explicit entry in group_tags.
    m2 = @allosteric_mechanism begin
        substrates: S
        products: P
        site(:catalytic, 2): begin
            steps: begin
                E + S ⇌ E_S       :: OnlyR
                E + P ⇌ E_P       :: EqualRT
                E_S <--> E_P      :: NonequalRT
            end
        end
    end
    spec2 = allosteric_spec_from_mechanism(m2, rxn1)
    @test AllostericEnzymeMechanism(spec2) === m2
    @test spec2.group_tags == Dict(1 => :OnlyR, 2 => :EqualRT, 3 => :NonequalRT)

    # With one allosteric regulator at its own site, tag :OnlyT.
    rxn3 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R
        oligomeric_state: 2
    end
    m3 = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R::OnlyT
        site(:catalytic, 2): begin
            steps: begin
                E + S ⇌ E_S       :: EqualRT
                E + P ⇌ E_P       :: EqualRT
                E_S <--> E_P      :: EqualRT
            end
        end
    end
    spec3 = allosteric_spec_from_mechanism(m3, rxn3)
    @test AllostericEnzymeMechanism(spec3) === m3
    @test spec3.reg_ligand_tags == Dict(:R => :OnlyT)

    # Two regulators at the same site, one :OnlyR, one :NonequalRT.
    # Dense storage: both ligands appear in reg_ligand_tags.
    rxn4 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R1, R2
        oligomeric_state: 2
    end
    m4 = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R1::OnlyR, R2::NonequalRT
        site(:catalytic, 2): begin
            steps: begin
                E + S ⇌ E_S       :: EqualRT
                E + P ⇌ E_P       :: EqualRT
                E_S <--> E_P      :: EqualRT
            end
        end
        site(:regulatory, 2): begin
            ligands: R1, R2
        end
    end
    spec4 = allosteric_spec_from_mechanism(m4, rxn4)
    @test AllostericEnzymeMechanism(spec4) === m4
    @test spec4.reg_ligand_tags == Dict(:R1 => :OnlyR, :R2 => :NonequalRT)
end
```

- [ ] **Step 1.4: Run tests and fix any helper issues**

Run: `julia --project -e 'using Pkg; Pkg.test()'`

Expected: PASS. If a `===` round-trip fails, the helper has a bug — investigate the affected field (group_tags ordering, reg sites order, etc). If an accessor isn't found, fix per Step 1.2.

- [ ] **Step 1.5: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test: add allosteric_spec_from_mechanism round-trip helper

Symmetric to the existing mechanism_spec_from_mechanism. Round-trips a
compiled AllostericEnzymeMechanism back to an AllostericMechanismSpec
with dense Dict storage (every kinetic group and regulator ligand has an
explicit entry). Validates the round-trip across four shapes via ===
assertions.

Required for upcoming per-move tests that need to seed AllostericMechanismSpec
inputs from @allosteric_mechanism literals.
EOF
)"
```

---

## Task 2: Section 1 — reorganize support-function tests in pipeline order

**Note on stoichiometry-2:** the `EnzymeReaction` constructor at `src/types.jl:46-49` errors on duplicate substrate or product names, and the macros pass through to that constructor. Reactions like `2 ADP ↔ ATP + AMP` cannot be expressed in the current API. Adding stoichiometric-coefficient support is out of scope for this rewrite. No stoich-2 seeds appear in any task below.

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

**Pre-existing tests this absorbs:** the existing testsets that test support functions (no spec input). These move under a `─── 1. Support functions ───` header in pipeline order; assertions are preserved verbatim.

Source line ranges of testsets to move (verify before each move):

| Testset | Existing lines |
|---|---|
| `Catalytic topologies` (incl. all sub-testsets) | 245–507 |
| `C6 iso size limit blocks 4x4` | 2408–2426 |
| `Competition patterns` | 685–739 |
| `Inhibitor competition patterns` | 1101–1128 |
| `Forms with binding step` | 1130–1174 |
| `Dead-end filtering by competition` (sub-testsets that test `_substrate_product_dead_end_opportunities`) | 741–921 (need to split: `_substrate_product_dead_end_opportunities` portions) |
| `Dead-end substrate/product expansion` (testset for `_expand_substrate_product_dead_ends`) | 576–683 |

- [ ] **Step 2.1: Insert the section header at the top of the test scope**

Locate `@testset "Mechanism Enumeration" begin` (around line 196). Add immediately after:

```julia
# ═══════════════════════════════════════════════════════════════════════
# 1. Support functions (no spec input)
# ═══════════════════════════════════════════════════════════════════════
```

- [ ] **Step 2.2: Move testsets in pipeline order**

Cut and paste the following testsets into the new section, in this order. Within each cut/paste, verify the original line range matches before moving. Preserve all assertions verbatim.

1. `_catalytic_topologies` (lines 245–507) — already a `Catalytic topologies` testset.
2. `C6 iso size limit blocks 4x4` (lines 2408–2426) — fold into the `_catalytic_topologies` testset as a sub-testset. New name: `quad-quad: C6 forces ping-pong`.
3. `_competition_patterns` (lines 685–739).
4. `_inhibitor_competition_patterns` (lines 1101–1128).
5. `_forms_with_binding_step` (lines 1130–1174).
6. `_substrate_product_dead_end_opportunities` (extract the relevant sub-testset from lines 741–921 — specifically the `Ter-ter diagonal: 12 of 27 allowed` portion that exercises this function directly).
7. `_expand_substrate_product_dead_ends` (lines 576–683).

After moves, add a divider before each:

```julia
# ─── _catalytic_topologies ─────────────────────────────────────────────
@testset "_catalytic_topologies" begin
    …
end

# ─── _competition_patterns ─────────────────────────────────────────────
@testset "_competition_patterns" begin
    …
end
```

- [ ] **Step 2.3: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`

Expected: PASS (no assertion changes; reorganization only).

If FAIL, the rearrangement broke a `const` reference or testset boundary. The most common cause: the cut-paste lost a `const <reaction>_rxn = @enzyme_reaction …` definition. Move it back to the infrastructure section (`─── 0. Test infrastructure ───`).

- [ ] **Step 2.4: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test(reorg): place support-function tests under section 1 in pipeline order

Reorganization only — no assertion changes. Tests for _catalytic_topologies
(absorbing the C6 quad-quad subcase), _competition_patterns,
_inhibitor_competition_patterns, _forms_with_binding_step,
_substrate_product_dead_end_opportunities, and
_expand_substrate_product_dead_ends now live under section 1 of the file
in pipeline order.
EOF
)"
```

---

## Task 4: Section 2 — rewrite init_mechanisms and add compile_mechanism testset

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

**Pre-existing tests this absorbs:**
- `n_fit_params_estimate semantics` (lines 119–133) — move into init testset.
- `n_fit_params_estimate upper-bound for dead-end mirrors` (lines 135–162) — keep, place in init testset.
- `init_mechanisms` testset block (lines 509–937) — absorb the parts that ARE init tests (Param count invariant, All have exactly 1 SS step, Mirror steps share kinetic_group, Uni-uni: no dead-end forms, Round-trip: competition-filtered specs compile).
- `init_mechanisms drops unbound regulators from spec→type` (lines 2428–2449) — fold into init testset.

The compile_mechanism round-trip pattern is currently scattered as one-liners (`@test EnzymeMechanism(spec) === m`) at every call site. Promote to a dedicated testset that asserts the round-trip explicitly across 4–5 representative seeds.

- [ ] **Step 4.1: Insert section 2 header**

After the close of section 1, add:

```julia
# ═══════════════════════════════════════════════════════════════════════
# 2. Initialization (compile_mechanism + init_mechanisms)
# ═══════════════════════════════════════════════════════════════════════
```

- [ ] **Step 4.2: Add `compile_mechanism` round-trip testset**

```julia
# ─── compile_mechanism / EnzymeMechanism round-trip ────────────────────
@testset "compile_mechanism round-trip" begin
    # Round-trip lossless invariant: for any mechanism built via the DSL,
    # mechanism_spec_from_mechanism ∘ EnzymeMechanism (== compile_mechanism)
    # returns the same singleton type. Validates the helper AND the
    # constructor's bidirectional consistency. Same idea for the allosteric
    # round-trip (covered by the dedicated testset added in Task 1).

    # uni-uni
    m_uu = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ E_P
            E + S ⇌ E_S
            E_S <--> E_P
        end
    end
    spec_uu = mechanism_spec_from_mechanism(m_uu, uni_uni_rxn)
    @test EnzymeMechanism(spec_uu) === m_uu

    # bi-bi sequential
    m_seq = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E_A
            E_A + B ⇌ E_A_B
            E + Q ⇌ E_Q
            E_Q + P ⇌ E_P_Q
            E_A_B <--> E_P_Q
        end
    end
    spec_seq = mechanism_spec_from_mechanism(m_seq, bi_bi_rxn)
    @test EnzymeMechanism(spec_seq) === m_seq

    # bi-bi ping-pong
    m_pp = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E_A
            Estar + B ⇌ Estar_B
            E + Q ⇌ E_Q
            Estar + P ⇌ Estar_A_P
            E_A <--> Estar_A_P
            Estar_B ⇌ E_Q
        end
    end
    spec_pp = mechanism_spec_from_mechanism(m_pp, bi_bi_pp_rxn)
    @test EnzymeMechanism(spec_pp) === m_pp

    # uni-uni with dead-end inhibitor (regulator strip in the round-trip)
    m_uu_i = @enzyme_mechanism begin
        substrates: S
        products: P
        regulators: I
        steps: begin
            E + P ⇌ E_P
            E + S ⇌ E_S
            E_S <--> E_P
            E + I ⇌ E_I
        end
    end
    spec_uu_i = mechanism_spec_from_mechanism(m_uu_i, uni_uni_with_reg)
    @test EnzymeMechanism(spec_uu_i) === m_uu_i
end
```

- [ ] **Step 4.3: Replace existing init_mechanisms testset with the rewritten one**

Delete lines 509–937 (the `init_mechanisms` testset block) and replace with:

```julia
# ─── init_mechanisms ───────────────────────────────────────────────────
@testset "init_mechanisms" begin

    @testset "min param count floor: subs + prods + 1" begin
        # Floor invariant: n_fit_params_estimate ≥ n_subs + n_prods + 1.
        # Derivation: every init mechanism has 1 SS step (the iso) plus
        # n_subs RE binding groups for substrates and n_prods for products.
        # The kinetic-group count = n_subs + n_prods (RE) + 1 (SS), and
        # n_thermo subtracts off based on cycle counts. The minimum after
        # subtractions equals exactly n_subs + n_prods + 1 for the simplest
        # topology with no dead-ends.
        for (rxn, n_s, n_p) in [
            (uni_uni_rxn, 1, 1),
            (uni_bi_rxn, 1, 2),
            (bi_bi_rxn, 2, 2),
            (bi_bi_pp_rxn, 2, 2),
        ]
            specs = EnzymeRates.init_mechanisms(rxn)
            min_pc = n_s + n_p + 1
            for s in specs
                @test s.n_fit_params_estimate >= min_pc
            end
        end
    end

    @testset "n_fit_params_estimate matches fitted_params for uni-uni init" begin
        # Uni-uni: 3 forms (E, E_S, E_P), 3 steps. n_thermo = 3 - 3 + 1 = 1
        # (one independent thermodynamic constraint = Keq).
        # Formula: n_re_groups + 2*n_ss_groups - n_thermo = 2 + 2 - 1 = 3.
        # length(fitted_params(m)) for uni-uni init = 3 (K1, K2, k3f).
        # Estimate must equal actual on the simplest case.
        init_specs = EnzymeRates.init_mechanisms(uni_uni_rxn)
        @test !isempty(init_specs)
        spec = first(init_specs)
        m = EnzymeRates.compile_mechanism(spec)
        n_actual = length(EnzymeRates.fitted_params(m))
        @test spec.n_fit_params_estimate == n_actual
    end

    @testset "n_fit_params_estimate upper-bound for dead-end mirrors" begin
        # When dead-end mirror cycles exist, the formula can underestimate
        # the true thermodynamic-constraint count, so the floor in
        # _apply_equivalence_grouping ensures pc ≥ n_subs + n_prods + 1.
        # This guards the upper-bound invariant: estimate ≥ actual.
        # Cap compiled specs to keep @generated cost bounded.
        cap = 30
        init_specs = EnzymeRates.init_mechanisms(uni_uni_with_reg)
        for spec in init_specs[1:min(cap, end)]
            m = EnzymeRates.compile_mechanism(spec)
            @test spec.n_fit_params_estimate >=
                length(EnzymeRates.fitted_params(m))
        end
        expanded = EnzymeRates.expand_mechanisms(
            init_specs, uni_uni_with_reg)
        expanded_specs = EnzymeRates.AbstractMechanismSpec[]
        for (_, specs) in expanded
            append!(expanded_specs, specs)
        end
        for spec in expanded_specs[1:min(cap, end)]
            m = EnzymeRates.compile_mechanism(spec)
            @test spec.n_fit_params_estimate >=
                length(EnzymeRates.fitted_params(m))
        end
    end

    @testset "exactly 1 SS step per init spec" begin
        # init_mechanisms produces minimum-parameter mechanisms — exactly
        # one isomerization step, which is SS by construction. Subsequent
        # RE→SS expansions add more SS steps; init never does.
        for rxn in [uni_uni_rxn, uni_bi_rxn,
                    bi_bi_rxn, bi_bi_pp_rxn]
            specs = EnzymeRates.init_mechanisms(rxn)
            for s in specs
                @test count(!st.is_equilibrium for st in s.steps) == 1
            end
        end
    end

    @testset "Same-metabolite RE bindings share kinetic_group" begin
        # _apply_equivalence_grouping collapses all RE binding steps for
        # the same metabolite into one kinetic group (one shared K).
        # For uni-uni + dead-end inhibitor, the inhibitor's mirror cycles
        # mean :I binds at multiple forms — these mirror bindings must
        # share a single kinetic_group (one K_I).
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            dead_end_inhibitors: I
        end
        specs = EnzymeRates.init_mechanisms(rxn)
        @test !isempty(specs)
        for spec in specs
            by_metabolite = Dict{Symbol, Vector{EnzymeRates.StepSpec}}()
            for step in spec.steps
                step.is_equilibrium || continue
                length(step.reactants) == 2 || continue
                met = step.reactants[2]
                push!(get!(by_metabolite, met,
                           EnzymeRates.StepSpec[]), step)
            end
            for (_met, steps) in by_metabolite
                length(steps) >= 2 || continue
                groups = Set(s.kinetic_group for s in steps)
                @test length(groups) == 1
            end
        end
    end

    @testset "Uni-uni: exactly 1 init mechanism" begin
        # Uni-uni topology: 1 catalytic topology × 1 dead-end variant
        # (none possible — see test_expand_substrate_product_dead_ends
        # uni-uni case). Hence init produces exactly 1 spec.
        specs = EnzymeRates.init_mechanisms(uni_uni_rxn)
        @test length(specs) == 1
    end

    @testset "Init compiles for all small reactions" begin
        # Every init spec must compile to a valid EnzymeMechanism, and
        # the actual fitted-param count must respect the upper-bound
        # invariant. Tests first 5 specs per reaction to cap @generated cost.
        for rxn in [uni_uni_rxn, bi_bi_rxn, bi_bi_pp_rxn]
            specs = EnzymeRates.init_mechanisms(rxn)
            for spec in first(specs, 5)
                m = EnzymeMechanism(spec)
                @test m isa EnzymeMechanism
                @test length(EnzymeRates.fitted_params(m)) <=
                    spec.n_fit_params_estimate
            end
        end
    end

    @testset "Drops unbound regulators from spec→type" begin
        # init_mechanisms produces specs without dead-end regulators bound.
        # When compiled to EnzymeMechanism, the regulator must NOT appear
        # in the regulators tuple — only the catalytic mechanism is built.
        # After expand_mechanisms adds the dead-end regulator, it should
        # appear.
        init_specs = EnzymeRates.init_mechanisms(uni_uni_with_reg)
        @test !isempty(init_specs)
        for spec in init_specs
            m = EnzymeRates.EnzymeMechanism(spec)
            @test :I ∉ EnzymeRates.regulators(m)
        end

        expanded = EnzymeRates.expand_mechanisms(init_specs, uni_uni_with_reg)
        found_with_reg = false
        for (_, specs) in expanded
            for spec in specs
                m = EnzymeRates.EnzymeMechanism(spec)
                if :I in EnzymeRates.regulators(m)
                    found_with_reg = true
                    break
                end
            end
            found_with_reg && break
        end
        @test found_with_reg
    end
end
```

- [ ] **Step 4.4: Run tests and commit**

Run: `julia --project -e 'using Pkg; Pkg.test()'`

Expected: PASS (assertion content unchanged; reorganization only).

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test(reorg): rewrite init_mechanisms tests; add compile_mechanism testset

Promotes the scattered EnzymeMechanism(spec) === m round-trip assertions
into a dedicated compile_mechanism testset across 4 representative
shapes. Reorganizes the init_mechanisms testset under section 2 in
pipeline order, absorbing n_fit_params_estimate-semantics tests and the
"drops unbound regulators" testset (formerly at file end). Each count
and delta carries an independent-derivation comment.
EOF
)"
```

---

## Task 5: Section 3a — `_expand_re_to_ss` (template-validating commit)

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

**Pre-existing tests this absorbs:**
- Existing `RE→SS conversion` testset (lines 939–1034).

This task is the **template validator** — once these patterns work, subsequent move tasks (6–10) follow the same shape. If a pattern problem surfaces here, fix it BEFORE proceeding to Task 6.

- [ ] **Step 5.1: Insert section 3 header**

After the close of section 2, add:

```julia
# ═══════════════════════════════════════════════════════════════════════
# 3. Base-spec expansion moves (polymorphic over Mechanism/AllostericMechanismSpec)
# ═══════════════════════════════════════════════════════════════════════
```

- [ ] **Step 5.2: Delete the existing RE→SS testset**

Locate and delete lines 939–1034 (the `RE→SS conversion` testset block).

- [ ] **Step 5.3: Add the rewritten testset for `_expand_re_to_ss`**

Drop in below the section 3 header:

```julia
# ─── _expand_re_to_ss ──────────────────────────────────────────────────
@testset "_expand_re_to_ss" begin

    @testset "MechanismSpec — uni-uni: 2 RE binding groups → 2 variants" begin
        # SEED: uni-uni with 3 singleton kinetic groups.
        # Group 1 = E+P binding (RE), group 2 = E+S binding (RE),
        # group 3 = iso E_S↔E_P (SS).
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism(m_seed, uni_uni_rxn)
        @test EnzymeMechanism(spec) === m_seed

        result = EnzymeRates._expand_re_to_ss(spec)

        # 1. count: RE→SS fires per all-RE kinetic group atomically.
        # The seed has 2 all-RE groups (P-binding, S-binding). The iso
        # group is already SS so it's excluded. → 2 variants.
        @test length(result) == 2

        # 2. Δ params: each conversion replaces 1 RE param (K) with 2 SS
        # params (kf, kr). For a plain MechanismSpec, _re_to_ss_delta = +1
        # (ratchet of 1 K → kf + kr is a +1 net under the kinetic-group
        # accounting in _n_fit_params_estimate_from_steps).
        for r in result
            @test r.n_fit_params_estimate ==
                spec.n_fit_params_estimate + 1
        end

        # 3. compilability
        for r in result
            @test compile_mechanism(r) isa EnzymeMechanism
        end

        # 4. structural change — equivalence-style (N=2 ≤ 6).
        # Variant A: P-binding flipped to SS (group 1 RE→SS).
        # Variant B: S-binding flipped to SS (group 2 RE→SS).
        # No third variant exists because the iso group was already SS.
        v_p_flipped = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P <--> E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        v_s_flipped = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S <--> E_S
                E_S <--> E_P
            end
        end
        expected = Set([v_p_flipped, v_s_flipped])
        @test Set(compile_mechanism(r) for r in result) == expected

        # 5. preservation: reaction unchanged; non-flipped steps remain
        # in their original RE/SS state with the same kinetic_group.
        for r in result
            @test r.reaction === spec.reaction
            # Exactly one step is now SS-with-metabolite that was RE in seed.
            n_newly_ss = count(zip(spec.steps, r.steps)) do (s_old, s_new)
                s_old.is_equilibrium && !s_new.is_equilibrium &&
                    s_old.kinetic_group == s_new.kinetic_group
            end
            @test n_newly_ss == 1
        end
    end

    @testset "MechanismSpec — bi-bi sequential: 2 RE binding groups → 2 variants" begin
        # SEED: bi-bi sequential. 2 binding groups (one for A, one for B
        # via parens to share kinetic group; same for P, Q). But here
        # we use the simplest sequential bi-bi where each metabolite has
        # its own singleton group — that gives 4 RE groups + 1 SS iso.
        # Sequential ordered: E + A → E_A + B → E_A_B ↔ E_P_Q → E + P/Q
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E_A
                E_A + B ⇌ E_A_B
                E + Q ⇌ E_Q
                E_Q + P ⇌ E_P_Q
                E_A_B <--> E_P_Q
            end
        end
        spec = mechanism_spec_from_mechanism(m_seed, bi_bi_rxn)
        @test EnzymeMechanism(spec) === m_seed

        result = EnzymeRates._expand_re_to_ss(spec)

        # 1. count: 4 all-RE singleton groups (A, B, Q, P bindings).
        # Iso group is already SS. → 4 variants, one per RE group.
        @test length(result) == 4

        # 2. Δ params: +1 per variant (plain MechanismSpec, no allosteric).
        for r in result
            @test r.n_fit_params_estimate ==
                spec.n_fit_params_estimate + 1
        end

        # 3. compilability
        for r in result
            @test compile_mechanism(r) isa EnzymeMechanism
        end

        # 4. structural change — N=4 ≤ 6 → equivalence-style.
        # Each variant flips exactly one of the four RE binding groups.
        v_a = @enzyme_mechanism begin
            substrates: A, B; products: P, Q
            steps: begin
                E + A <--> E_A
                E_A + B ⇌ E_A_B
                E + Q ⇌ E_Q
                E_Q + P ⇌ E_P_Q
                E_A_B <--> E_P_Q
            end
        end
        v_b = @enzyme_mechanism begin
            substrates: A, B; products: P, Q
            steps: begin
                E + A ⇌ E_A
                E_A + B <--> E_A_B
                E + Q ⇌ E_Q
                E_Q + P ⇌ E_P_Q
                E_A_B <--> E_P_Q
            end
        end
        v_q = @enzyme_mechanism begin
            substrates: A, B; products: P, Q
            steps: begin
                E + A ⇌ E_A
                E_A + B ⇌ E_A_B
                E + Q <--> E_Q
                E_Q + P ⇌ E_P_Q
                E_A_B <--> E_P_Q
            end
        end
        v_p = @enzyme_mechanism begin
            substrates: A, B; products: P, Q
            steps: begin
                E + A ⇌ E_A
                E_A + B ⇌ E_A_B
                E + Q ⇌ E_Q
                E_Q + P <--> E_P_Q
                E_A_B <--> E_P_Q
            end
        end
        expected = Set([v_a, v_b, v_q, v_p])
        @test Set(compile_mechanism(r) for r in result) == expected

        # 5. preservation: each result has exactly one step's is_equilibrium
        # flipped from true to false, with kinetic_group preserved.
        for r in result
            @test r.reaction === spec.reaction
            n_newly_ss = count(zip(spec.steps, r.steps)) do (s_old, s_new)
                s_old.is_equilibrium && !s_new.is_equilibrium &&
                    s_old.kinetic_group == s_new.kinetic_group
            end
            @test n_newly_ss == 1
        end
    end

    @testset "MechanismSpec — bi-bi multi-step kinetic group: atomic conversion" begin
        # SEED: bi-bi random where A binds at two forms (E and E_B) and
        # those two RE binding steps share kinetic_group 1 (parenthesized).
        # B-binding shares group 2 (E and E_A). P shares group 3.
        # Q shares group 4. Iso = group 5 (SS).
        # When _expand_re_to_ss fires on group 1, BOTH A-binding steps
        # flip atomically (same group → same kinetic params).
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E_A
                 E_B + A ⇌ E_A_B)
                (E + B ⇌ E_B
                 E_A + B ⇌ E_A_B)
                (E + P ⇌ E_P
                 E_P + Q ⇌ E_P_Q)
                (E + Q ⇌ E_Q
                 E_Q + P ⇌ E_P_Q)
                E_A_B <--> E_P_Q
            end
        end
        spec = mechanism_spec_from_mechanism(m_seed, bi_bi_rxn)
        @test EnzymeMechanism(spec) === m_seed

        result = EnzymeRates._expand_re_to_ss(spec)

        # 1. count: 4 multi-step RE groups (A, B, P, Q each with 2 steps).
        # Iso group SS. → 4 variants, each flipping 2 steps atomically.
        @test length(result) == 4

        # 2. Δ params: +1 per variant.
        for r in result
            @test r.n_fit_params_estimate ==
                spec.n_fit_params_estimate + 1
        end

        # 3. compilability
        for r in result
            @test compile_mechanism(r) isa EnzymeMechanism
        end

        # 4. property-style check: in each variant, exactly one kinetic
        # group has ALL its steps now SS (atomic conversion). All other
        # groups retain their original state.
        for r in result
            groups = Dict{Int, Vector{Bool}}()
            for st in r.steps
                push!(get!(groups, st.kinetic_group, Bool[]), st.is_equilibrium)
            end
            # Exactly one group: all-false (newly SS, was multi-step RE).
            n_all_ss_multi = count(((_, vs),) ->
                length(vs) >= 2 && all(==(false), vs), groups)
            @test n_all_ss_multi == 1
        end

        # 5. preservation: reaction unchanged.
        for r in result
            @test r.reaction === spec.reaction
        end
    end

    @testset "MechanismSpec — all-SS catalytic seed: empty (negative)" begin
        # When every catalytic step is already SS, _expand_re_to_ss has no
        # all-RE group to fire on → empty result.
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P <--> E_P
                E + S <--> E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism(m_seed, uni_uni_rxn)
        @test EnzymeMechanism(spec) === m_seed
        @test isempty(EnzymeRates._expand_re_to_ss(spec))
    end

    @testset "AllostericMechanismSpec — :EqualRT group: Δ=+1" begin
        # SEED: uni-uni with all groups :EqualRT. Each catalytic group's
        # R/T tag is :EqualRT (one shared K_R = K_T). When RE→SS converts
        # an :EqualRT group, the new (kf, kr) pair is also state-shared,
        # so Δ = +1 (the EqualRT/OnlyR/OnlyT cheap-tag delta).
        m_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P       :: EqualRT
                    E + S ⇌ E_S       :: EqualRT
                    E_S <--> E_P      :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism(m_seed, uni_uni_allo)
        @test AllostericEnzymeMechanism(spec) === m_seed

        result = EnzymeRates._expand_re_to_ss(spec)

        # 1. count: 2 all-RE groups (P-binding, S-binding); iso is SS.
        # _expand_re_to_ss fires per group; same as plain. → 2 variants.
        @test length(result) == 2

        # 2. Δ params: :EqualRT is a cheap tag → +1 per variant.
        for r in result
            @test r.n_fit_params_estimate ==
                spec.n_fit_params_estimate + 1
        end

        # 3. compilability — must produce AllostericEnzymeMechanism.
        for r in result
            @test compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. property-style: exactly one group now SS; all group_tags
        # preserved including the converted group's :EqualRT tag (move
        # MUST NOT change R/T-state semantics).
        for r in result
            n_newly_ss = count(zip(spec.base.steps, r.base.steps)) do (s_old, s_new)
                s_old.is_equilibrium && !s_new.is_equilibrium
            end
            @test n_newly_ss == 1
            @test r.group_tags == spec.group_tags  # tags untouched
        end

        # 5. preservation: catalytic_n, reg sites, reg ligand tags untouched.
        for r in result
            @test r.catalytic_n == spec.catalytic_n
            @test r.allosteric_reg_sites == spec.allosteric_reg_sites
            @test r.allosteric_multiplicities == spec.allosteric_multiplicities
            @test r.reg_ligand_tags == spec.reg_ligand_tags
            @test r.base.reaction === spec.base.reaction
        end
    end

    @testset "AllostericMechanismSpec — :NonequalRT group: Δ=+2" begin
        # SEED: uni-uni with one :NonequalRT group, others :EqualRT.
        # When RE→SS converts the :NonequalRT group, BOTH the R-state K
        # and the T-state K_T must split into (kf, kr) and (kf_T, kr_T).
        # Δ for :NonequalRT = 2 × base = +2.
        m_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P       :: EqualRT
                    E + S ⇌ E_S       :: NonequalRT
                    E_S <--> E_P      :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism(m_seed, uni_uni_allo)
        @test AllostericEnzymeMechanism(spec) === m_seed

        result = EnzymeRates._expand_re_to_ss(spec)

        # 1. count: 2 RE groups (group 1 = P-binding :EqualRT,
        # group 2 = S-binding :NonequalRT). → 2 variants.
        @test length(result) == 2

        # 2. Δ params: depends on tag of the converted group.
        # P-binding :EqualRT → +1; S-binding :NonequalRT → +2.
        deltas = sort([r.n_fit_params_estimate -
                       spec.n_fit_params_estimate for r in result])
        @test deltas == [1, 2]

        # 3. compilability
        for r in result
            @test compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. property-style: exactly one RE group flipped per variant; tag
        # of the flipped group preserved.
        for r in result
            @test r.group_tags == spec.group_tags
        end

        # 5. preservation
        for r in result
            @test r.catalytic_n == spec.catalytic_n
            @test r.allosteric_reg_sites == spec.allosteric_reg_sites
        end
    end

    @testset "Substrate-as-dead-end-inhibitor overlap (S used as both)" begin
        # SEED: uni-uni where S is BOTH a substrate AND a dead-end inhibitor.
        # The reaction declares dead_end_inhibitors: S, and the seed has
        # been pre-expanded to bind :S as inhibitor (giving rise to S__reg
        # binding steps). The base RE→SS move shouldn't be confused by
        # the metabolite-overlap — it operates on kinetic groups, not names.
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: S
        end
        # Build via init + dead-end expansion to get the S/__reg overlap form.
        init_specs = EnzymeRates.init_mechanisms(rxn)
        @test length(init_specs) == 1   # uni-uni: 1 catalytic topology
        seed_spec = first(init_specs)
        de_specs = EnzymeRates._expand_add_dead_end_regulator(seed_spec, rxn)
        @test !isempty(de_specs)
        spec = first(de_specs)

        # Move
        result = EnzymeRates._expand_re_to_ss(spec)

        # 1. count: each all-RE kinetic group can flip. The seed after
        # add-dead-end has groups: substrate-binding (RE), product-binding
        # (RE), iso (SS), and dead-end-S__reg-binding (RE). → 3 RE groups → 3 variants.
        @test length(result) == 3

        # 2. Δ params: +1 each (plain MechanismSpec).
        for r in result
            @test r.n_fit_params_estimate ==
                spec.n_fit_params_estimate + 1
        end

        # 3. compilability
        for r in result
            @test compile_mechanism(r) isa EnzymeMechanism
        end

        # 5. preservation: reaction === spec.reaction.
        for r in result
            @test r.reaction === spec.reaction
        end
    end

end
```

- [ ] **Step 5.4: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`

Expected: PASS. If FAIL on a count, delta, or equivalence-set assertion, **stop and follow the bug-handling protocol** in design §6:

1. Diagnose: read the assertion, the seed, the actual output. Form ONE hypothesis (test wrong vs code wrong).
2. Surface to Denis with the structured report (file:line, seed, expected vs actual, hypothesis).
3. Do NOT modify the test to match buggy output. Do NOT use `@test_broken`.
4. If code is wrong: fix `src/` in a separate "fix:" commit FIRST, then this commit lands green.
5. If test expectation is wrong: fix the derivation comment to explain why I miscounted, then update the test.

- [ ] **Step 5.5: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test(re_to_ss): rewrite with literal seeds and 7-item checklist

Replaces the existing RE→SS conversion testset (lines 939-1034) with the
new pattern: literal @enzyme_mechanism / @allosteric_mechanism seeds,
independent-derivation comments for every count and delta, equivalence-style
structural assertions for ≤6-output cases. Adds:

- bi-bi sequential (4 variants, equivalence-style)
- bi-bi multi-step group (atomic conversion, property-style)
- AllostericMechanismSpec :EqualRT (Δ=+1) and :NonequalRT (Δ=+2)
- substrate-as-dead-end-inhibitor overlap

Acts as the template-validating commit for the section 3/4 rewrite.
EOF
)"
```

If a `src/` bug surfaced, the bug fix lands as a separate "fix:" commit BEFORE this one, per the design §6 protocol.

---

## Task 6: Section 3b — `_expand_split_kinetic_group`

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

**Pre-existing tests this absorbs:**
- `Split kinetic group` testset (lines 1036–1099).
- `Remove constraint on allosteric` sub-testset (lines 2309–2375) — fold into the cross-type allosteric sub-testsets.

The template is now established by Task 5 — apply the same pattern. Each seed below specifies the literal mechanism, the count derivation, the delta derivation, and whether equivalence-style or property-style is appropriate.

- [ ] **Step 6.1: Delete the existing `Split kinetic group` testset (lines 1036–1099) and the `Remove constraint on allosteric` sub-testset (lines 2309–2375)**

- [ ] **Step 6.2: Add the rewritten testset for `_expand_split_kinetic_group`**

Insert below `_expand_re_to_ss`:

```julia
# ─── _expand_split_kinetic_group ───────────────────────────────────────
@testset "_expand_split_kinetic_group" begin

    @testset "MechanismSpec — bi-bi: 4 multi-step groups → 8 splits" begin
        # SEED: bi-bi random with 4 multi-step kinetic groups (A, B, P, Q
        # each with 2 binding steps shared via parens). Iso is singleton.
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E_A
                 E_B + A ⇌ E_A_B)
                (E + B ⇌ E_B
                 E_A + B ⇌ E_A_B)
                (E + P ⇌ E_P
                 E_P + Q ⇌ E_P_Q)
                (E + Q ⇌ E_Q
                 E_Q + P ⇌ E_P_Q)
                E_A_B <--> E_P_Q
            end
        end
        spec = mechanism_spec_from_mechanism(m_seed, bi_bi_rxn)
        @test EnzymeMechanism(spec) === m_seed

        result = EnzymeRates._expand_split_kinetic_group(spec)

        # 1. count: split picks one step out of a multi-step group into
        # a fresh group. 4 groups × 2 steps each → 4 × 2 = 8 split variants.
        @test length(result) == 8

        # 2. Δ params: each split adds +1 (RE plain, no allosteric).
        # The new group inherits is_equilibrium=true; one extra K param.
        for r in result
            @test r.n_fit_params_estimate ==
                spec.n_fit_params_estimate + 1
        end

        # 3. compilability
        for r in result
            @test compile_mechanism(r) isa EnzymeMechanism
        end

        # 4. property-style (N=8 > 6): the result has exactly one new
        # kinetic_group integer (max+1) and exactly one step now belongs
        # to it. The remaining group has size n_old - 1.
        for r in result
            old_groups = Set(s.kinetic_group for s in spec.steps)
            new_groups = Set(s.kinetic_group for s in r.steps)
            extra = setdiff(new_groups, old_groups)
            @test length(extra) == 1
            new_g = first(extra)
            n_in_new = count(s.kinetic_group == new_g for s in r.steps)
            @test n_in_new == 1
        end

        # 5. preservation: reaction unchanged; total step count unchanged.
        for r in result
            @test r.reaction === spec.reaction
            @test length(r.steps) == length(spec.steps)
        end
    end

    @testset "MechanismSpec — all singleton groups: empty (negative)" begin
        # SEED: uni-uni init, every kinetic group is a singleton (size 1).
        # Splitting requires a multi-step group → no eligible group → empty.
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism(m_seed, uni_uni_rxn)
        @test EnzymeMechanism(spec) === m_seed
        @test isempty(EnzymeRates._expand_split_kinetic_group(spec))
    end

    @testset "MechanismSpec — mixed RE/SS group sizes: deltas differ" begin
        # SEED: bi-bi where one multi-step group is now SS (after a prior
        # RE→SS conversion). Splitting an RE group adds +1; splitting an
        # SS group adds +2 (kf and kr both split into a fresh group).
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A <--> E_A          # SS group
                 E_B + A <--> E_A_B)
                (E + B ⇌ E_B
                 E_A + B ⇌ E_A_B)
                E + P ⇌ E_P
                E + Q ⇌ E_Q
                E_A_B <--> E_P_Q
                E_P + Q ⇌ E_P_Q
                E_Q + P ⇌ E_P_Q
            end
        end
        spec = mechanism_spec_from_mechanism(m_seed, bi_bi_rxn)
        @test EnzymeMechanism(spec) === m_seed

        result = EnzymeRates._expand_split_kinetic_group(spec)

        # 1. count: 2 multi-step groups (A-binding SS×2, B-binding RE×2).
        # Each can split per member → 2 × 2 = 4 variants.
        @test length(result) == 4

        # 2. Δ params: +1 if RE-group split, +2 if SS-group split.
        deltas = sort([r.n_fit_params_estimate -
                       spec.n_fit_params_estimate for r in result])
        # 2 SS splits (+2 each) + 2 RE splits (+1 each) = [1, 1, 2, 2]
        @test deltas == [1, 1, 2, 2]

        # 3. compilability
        for r in result
            @test compile_mechanism(r) isa EnzymeMechanism
        end

        # 5. preservation
        for r in result
            @test r.reaction === spec.reaction
            @test length(r.steps) == length(spec.steps)
        end
    end

    @testset "AllostericMechanismSpec — split inherits parent's tag" begin
        # SEED: bi-bi allosteric where one multi-step group is :NonequalRT
        # and another is :EqualRT. Split must produce results where the
        # NEW group inherits the parent's tag — splitting is a parameter-
        # relaxation move and MUST NOT change R/T-state semantics.
        m_seed = @allosteric_mechanism begin
            substrates: A, B
            products: P, Q
            site(:catalytic, 2): begin
                steps: begin
                    (E + A ⇌ E_A
                     E_B + A ⇌ E_A_B)        :: NonequalRT
                    (E + B ⇌ E_B
                     E_A + B ⇌ E_A_B)        :: EqualRT
                    E + P ⇌ E_P             :: EqualRT
                    E_P + Q ⇌ E_P_Q         :: EqualRT
                    E + Q ⇌ E_Q             :: EqualRT
                    E_Q + P ⇌ E_P_Q         :: EqualRT
                    E_A_B <--> E_P_Q        :: EqualRT
                end
            end
        end
        bi_bi_allo_rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            oligomeric_state: 2
        end
        spec = allosteric_spec_from_mechanism(m_seed, bi_bi_allo_rxn)
        @test AllostericEnzymeMechanism(spec) === m_seed

        result = EnzymeRates._expand_split_kinetic_group(spec)

        # 1. count: 2 multi-step groups × 2 members each = 4 variants.
        @test length(result) == 4

        # 2. Δ params: depends on tag and is_equilibrium of parent group.
        # :NonequalRT RE split: +2 (base 1 × 2 for NonequalRT factor).
        # :EqualRT RE split: +1.
        deltas = sort([r.n_fit_params_estimate -
                       spec.n_fit_params_estimate for r in result])
        @test deltas == [1, 1, 2, 2]

        # 3. compilability
        for r in result
            @test compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. property-style: the new group's tag equals its parent group's tag.
        pre_groups = Set(s.kinetic_group for s in spec.base.steps)
        for r in result
            post_groups = Set(s.kinetic_group for s in r.base.steps)
            new_groups = setdiff(post_groups, pre_groups)
            @test length(new_groups) == 1
            new_g = first(new_groups)
            # Identify parent: the only pre-group whose count dropped.
            pre_counts = Dict(g => count(s -> s.kinetic_group == g,
                                         spec.base.steps)
                              for g in pre_groups)
            post_counts = Dict(g => count(s -> s.kinetic_group == g,
                                          r.base.steps)
                               for g in pre_groups)
            old_g = only(g for g in pre_groups
                         if post_counts[g] < pre_counts[g])
            @test r.group_tags[new_g] == spec.group_tags[old_g]
            for g in pre_groups
                @test r.group_tags[g] == spec.group_tags[g]
            end
        end

        # 5. preservation
        for r in result
            @test r.catalytic_n == spec.catalytic_n
            @test r.allosteric_reg_sites == spec.allosteric_reg_sites
            @test r.reg_ligand_tags == spec.reg_ligand_tags
        end
    end
end
```

- [ ] **Step 6.3: Run tests; surface bugs per protocol; commit**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test(split_kinetic_group): rewrite with literal seeds and 7-item checklist

Replaces the existing Split kinetic group testset (lines 1036-1099) and
absorbs the Remove-constraint-on-allosteric sub-testset. Independent
derivations for count, delta, and tag-inheritance behavior across 4 seeds:
bi-bi multi-step plain, all-singleton (negative), mixed RE/SS group sizes,
and bi-bi allosteric (split must inherit parent tag).
EOF
)"
```

---

## Task 7: Section 3c — `_expand_add_dead_end_regulator`

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

**Pre-existing tests this absorbs:**
- `Add dead-end regulator` testset (lines 1176–1622).
- `Regulator dummy naming stability` (lines 1624–1665).
- `Metabolite overlap: substrate as dead-end inhibitor` (lines 2203–2231).
- `Add dead-end reg on allosteric` sub-testset (lines 2377–2406).

This is the largest move-test in the file. Use the template from Task 5; apply the seed battery from design §3.

- [ ] **Step 7.1: Delete the absorbed testsets**

Locate each testset by line and delete:
- `Add dead-end regulator` (~1176–1622)
- `Regulator dummy naming stability` (~1624–1665)
- `Metabolite overlap: substrate as dead-end inhibitor` (~2203–2231)
- `Add dead-end reg on allosteric` (~2377–2406)

- [ ] **Step 7.2: Add the rewritten testset**

Insert below `_expand_split_kinetic_group`. Each sub-testset follows Task 5's template. Seeds (each its own `@testset`):

1. **Uni-uni + new regulator I**: count = 1 variant. Derivation: E is the only eligible form (E_S has all subs, E_P has all prods, both ineligible). Equivalence-style.
2. **No regulators in reaction → empty (negative)**.
3. **Already-bound regulator → excluded**: seed has I already bound; calling the move with `exclude_regs=Set([:I])` yields empty.
4. **Sequential bi-bi + I**: count = 4 (derive from the inhibitor-competition pattern × eligible forms × dedup of identical form sets — see existing comments in current test at line 1487 lifted with derivation).
5. **Bi-bi random + I**: count = 9 (bi-bi has 7 dead-end forms × inhibitor patterns; full pattern at this seed produces 9 distinct dead-end-form sets).
6. **Bi-bi PP + I**: count = 3 (per existing test line 1458).
7. **Two regulators chain**: add I, then add J on a spec that already has I → asserts the J__reg dummy naming has no numeric suffix (current `Regulator dummy naming stability`).
8. **Two regulators competition**: count = 17 (per existing test). Property-style.
9. **Substrate-as-dead-end-inhibitor overlap**: seed declares `dead_end_inhibitors: S` where S is also a substrate. Move produces `S__reg`-named binding steps that compile correctly. Folds in the existing overlap testset.
10. **AllostericMechanismSpec input**: dead-end on allosteric with mixed regs (some allosteric, some dead-end). The move must exclude allosteric ligands from `eligible_regs` and tag the new dead-end group `:EqualRT`.
11. **Allosteric-only regulator → empty**: when reaction has only `allosteric_regulators:` and the spec is allosteric, no eligible dead-end regs → empty.

For each seed, write the same six-item checklist. For brevity in the plan, here is the structure for ONE seed; replicate for the others using the seed-battery descriptions above.

```julia
@testset "_expand_add_dead_end_regulator" begin

    @testset "Uni-uni + I: 1 variant (equivalence-style)" begin
        # SEED: uni-uni init, no regulators bound yet.
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
        end
        spec = MechanismSpec(rxn,
            mechanism_spec_from_mechanism(m_seed, uni_uni_rxn).steps,
            mechanism_spec_from_mechanism(m_seed, uni_uni_rxn).n_fit_params_estimate)
        # Round-trip on the catalytic side (excluding the unbound :I)
        @test EnzymeMechanism(spec) === m_seed

        result = EnzymeRates._expand_add_dead_end_regulator(spec, rxn)

        # 1. count: eligible forms for I-binding = forms NOT bound by all subs
        # AND NOT bound by all prods. E is the only eligible form (E_S has
        # all subs → ineligible; E_P has all prods → ineligible). Inhibitor
        # competition patterns for uni-uni: 1 (S × P × no-existing-inh = 1).
        # → 1 form set → 1 variant.
        @test length(result) == 1

        # 2. Δ params: +1 (one new K_I parameter for the dead-end binding group).
        for r in result
            @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
        end

        # 3. compilability
        for r in result
            @test compile_mechanism(r) isa EnzymeMechanism
        end

        # 4. equivalence-style (N=1).
        # Expected: same uni-uni catalytic + a new RE binding step E + I__reg ⇌ E_I.
        expected = @enzyme_mechanism begin
            substrates: S
            products: P
            regulators: I
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
                E + I ⇌ E_I
            end
        end
        @test compile_mechanism(first(result)) === expected

        # 5. preservation
        @test first(result).reaction === spec.reaction
    end

    @testset "Uni-uni no regulators → empty (negative)" begin
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism(m_seed, uni_uni_rxn)
        @test EnzymeMechanism(spec) === m_seed
        @test isempty(EnzymeRates._expand_add_dead_end_regulator(spec, uni_uni_rxn))
    end

    @testset "Sequential bi-bi + I: 4 distinct form sets" begin
        # SEED: bi-bi sequential. Eligible forms: E, E_A, E_Q.
        # Inhibitor competition patterns enumerate which subs/prods I
        # competes with; combined with which forms have a binding step
        # for those mets, produces exactly 4 distinct form sets:
        #   {E, E_Q}: ({A},{P}), ({A,B},{P})
        #   {E, E_A}: ({B},{Q}), ({B},{P,Q})
        #   {E_A, E_Q}: ({B},{P})
        #   {E}:      ({A},{Q}), ({A},{P,Q}), ({A,B},{Q}), ({A,B},{P,Q})
        # → 4 unique form sets after dedup.
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E_A
                E_A + B ⇌ E_A_B
                E + Q ⇌ E_Q
                E_Q + P ⇌ E_P_Q
                E_A_B <--> E_P_Q
            end
        end
        rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end
        base = mechanism_spec_from_mechanism(m_seed, bi_bi_rxn)
        spec = MechanismSpec(rxn, base.steps, base.n_fit_params_estimate)
        result = EnzymeRates._expand_add_dead_end_regulator(spec, rxn)

        @test length(result) == 4
        for r in result
            @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
            @test compile_mechanism(r) isa EnzymeMechanism
            @test r.reaction === spec.reaction
        end
    end

    # … (replicate the structure for seeds 5–11 from the seed battery above)
end
```

For seeds 5–11, the engineer follows the SAME `@testset` structure shown in seeds 1–3 above: literal seed → round-trip validation → move call → checklist items 1–5. Below is the per-seed data needed (seed mechanism literal, expected count, count derivation, expected delta, assertion style):

**Seed 4 — bi-bi random + I (9 variants, property-style since N>6):**

```julia
m_seed = @enzyme_mechanism begin
    substrates: A, B; products: P, Q
    steps: begin
        (E + A ⇌ E_A
         E_B + A ⇌ E_A_B)
        (E + B ⇌ E_B
         E_A + B ⇌ E_A_B)
        (E + P ⇌ E_P
         E_P + Q ⇌ E_P_Q)
        (E + Q ⇌ E_Q
         E_Q + P ⇌ E_P_Q)
        E_A_B <--> E_P_Q
    end
end
# Reaction with dead_end_inhibitors: I, atoms as bi_bi_rxn.
# Count derivation: bi-bi has 4 sub-binding forms (E, E_A, E_B, E_P, E_Q
# minus those bound by all subs/all prods) — eligible forms = E, E_A, E_B,
# E_P, E_Q (5). Inhibitor competition patterns: 9 (3 sub subsets × 3 prod
# subsets, comp_inh=∅). Each pattern produces a distinct active-form set
# after dedup → exactly 9 variants.
@test length(result) == 9
# Property-style: each result has at least one I__reg-suffixed binding step,
# all I__reg bindings share one kinetic group (one K_I).
```

**Seed 5 — bi-bi PP + I (3 variants, equivalence-style):**

```julia
m_seed = @enzyme_mechanism begin
    substrates: A, B; products: P, Q
    steps: begin
        E + A ⇌ E_A
        Estar + B ⇌ Estar_B
        E + Q ⇌ E_Q
        Estar + P ⇌ Estar_A_P
        E_A <--> Estar_A_P
        Estar_B ⇌ E_Q
    end
end
# Reaction: bi_bi_pp_rxn + dead_end_inhibitors: I.
# Count derivation: ping-pong has forms E, E_A, Estar, Estar_A_P, Estar_B,
# E_Q. Estar_A_P has all prods → ineligible. Estar_B has all subs →
# ineligible. Eligible: E, E_A, Estar, E_Q (4 forms). Inhibitor competition
# patterns × dedup → 3 unique form sets.
@test length(result) == 3
# Each variant must compile to a valid EnzymeMechanism. Equivalence-style
# requires writing 3 expected mechanisms; this is reasonable at N=3.
```

**Seed 6 — Two regulators chain (regulator dummy naming stability):**

```julia
m_seed = @enzyme_mechanism begin
    substrates: S; products: P
    steps: begin
        E + P ⇌ E_P
        E + S ⇌ E_S
        E_S <--> E_P
    end
end
# Reaction: dead_end_inhibitors: I, J.
# Step A: add I via _expand_add_dead_end_regulator.
# Step B: pick the variant with I bound, then call the move again — adds J.
# Property-style: in J's binding steps, the dummy must be exactly :J__reg
# (no numeric suffix like :J__reg2). Catches the regulator-naming-stability
# regression (existing test at line 1624).
for s in j_specs
    j_syms = [sym for st in s.steps
              for sym in Iterators.flatten((st.reactants, st.products))
              if contains(string(sym), "J__reg")]
    for sym in j_syms
        @test !occursin(r"__reg\d", string(sym))
    end
end
```

**Seed 7 — Two regulators competition (17 variants, property-style):**

Same base mechanism as seed 4 (bi-bi random) but with `dead_end_inhibitors: I1, I2`. After I1 binds at multiple forms, calling `_expand_add_dead_end_regulator` again adds I2 with patterns of competing-vs-not against I1 — count derivation: 9 base patterns × variations from I1's existing presence yields 17 unique form sets after dedup. (Existing test at line 1597 has the count, but the derivation needs to be re-derived from the inhibitor-pattern formula, not lifted.)

**Seed 8 — Substrate-as-dead-end-inhibitor overlap:**

Same shape as the overlap testset in Task 5 (substrate :S also declared as `dead_end_inhibitors: S`). The move produces `S__reg`-named binding steps; assert each variant compiles correctly and exactly one new kinetic group exists.

**Seed 9 — AllostericMechanismSpec input (mixed allosteric + dead-end):**

```julia
m_seed = @allosteric_mechanism begin
    substrates: S; products: P
    site(:catalytic, 2): begin
        steps: begin
            E + P ⇌ E_P    :: EqualRT
            E + S ⇌ E_S    :: EqualRT
            E_S <--> E_P   :: EqualRT
        end
    end
end
# Reaction: dead_end_inhibitors: I, allosteric_regulators: R, oligomeric_state: 2.
# Spec is allosteric. Move excludes :R (allosteric ligand) from
# eligible_regs but allows :I.
# Count derivation: same as uni-uni + I plain (1 eligible form E),
# but result is AllostericMechanismSpec.
@test length(result) == 1
for r in result
    @test r isa AllostericMechanismSpec
    # New dead-end binding kinetic group must be tagged :EqualRT
    # (cheapest tag, per design §3 cross-type rule).
    new_groups = setdiff(
        Set(s.kinetic_group for s in r.base.steps),
        Set(s.kinetic_group for s in spec.base.steps))
    @test length(new_groups) == 1
    new_g = first(new_groups)
    @test r.group_tags[new_g] == :EqualRT
end
```

**Seed 10 — Allosteric-only regulator → empty (negative):**

Reaction has `allosteric_regulators: R` only (no dead-end); spec is allosteric with R bound. `_expand_add_dead_end_regulator` should return empty — no eligible dead-end regulators.

**Seed 11 — Same regulator binding steps share one kinetic group:**

Bi-bi random + I as in seed 4. Property-style: across the result, every variant's I__reg binding steps (potentially multiple, e.g., E + I__reg ⇌ E_I and E_A + I__reg ⇌ E_A_I) all share the same kinetic_group integer (one K_I parameter, not two).

**Independent derivation reminder:** counts above (9, 3, 17, 1) must each be re-derived from the seed and the move's specification at write-time, NOT lifted from the existing tests. If the engineer writes the test, runs it, and gets a different number, the existing test's count may be wrong (bug-handling protocol §6 applies) OR the new derivation is wrong — surface either way, do not silently match.

- [ ] **Step 7.3: Run tests; surface bugs per protocol; commit**

```bash
julia --project -e 'using Pkg; Pkg.test()'
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test(add_dead_end_regulator): rewrite with literal seeds; absorb 4 testsets

Replaces the existing Add dead-end regulator (lines 1176-1622),
Regulator dummy naming stability (1624-1665), Metabolite overlap:
substrate as dead-end inhibitor (2203-2231), and Add dead-end reg on
allosteric (2377-2406) testsets. 11 seeds with independent-derivation
comments per count and delta. Includes substrate-as-I overlap and the
allosteric cross-type case.
EOF
)"
```

---

## Task 8: Section 4a — `_expand_to_allosteric`

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

**Pre-existing tests this absorbs:**
- `Allosteric conversion` testset (lines 1667–1723).

- [ ] **Step 8.1: Insert section 4 header**

After the close of section 3, add:

```julia
# ═══════════════════════════════════════════════════════════════════════
# 4. Allosteric expansion moves (AllostericMechanismSpec only)
# ═══════════════════════════════════════════════════════════════════════
```

- [ ] **Step 8.2: Delete `Allosteric conversion` (lines 1667–1723) and add the rewritten testset**

```julia
# ─── _expand_to_allosteric ─────────────────────────────────────────────
@testset "_expand_to_allosteric" begin

    @testset "MechanismSpec — uni-uni: 4 variants (equivalence-style)" begin
        # SEED: uni-uni init, 3 singleton kinetic groups.
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism(m_seed, uni_uni_allo)
        @test EnzymeMechanism(spec) === m_seed

        result = EnzymeRates._expand_to_allosteric(spec, uni_uni_allo)

        # 1. count: _expand_to_allosteric emits the all-:EqualRT baseline
        # once plus one :OnlyR variant per kinetic group. 3 groups → 1 + 3 = 4.
        @test length(result) == 4

        # 2. Δ params: +1 per variant (just L, the conformation equilibrium).
        # All other tag deltas are zero relative to the all-:EqualRT baseline.
        for r in result
            @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
        end

        # 3. compilability — every variant must build to AllostericEnzymeMechanism.
        for r in result
            @test compile_mechanism(r) isa AllostericEnzymeMechanism
            @test r.catalytic_n == 2
        end

        # 4. equivalence-style (N=4 ≤ 6). 4 expected mechanisms:
        # all-:EqualRT baseline + one :OnlyR variant per group.
        v_baseline = @allosteric_mechanism begin
            substrates: S; products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        v_g1_OnlyR = @allosteric_mechanism begin
            substrates: S; products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: OnlyR
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        v_g2_OnlyR = @allosteric_mechanism begin
            substrates: S; products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: OnlyR
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        v_g3_OnlyR = @allosteric_mechanism begin
            substrates: S; products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: OnlyR
                end
            end
        end
        @test Set(compile_mechanism(r) for r in result) ==
            Set([v_baseline, v_g1_OnlyR, v_g2_OnlyR, v_g3_OnlyR])

        # 5. preservation
        for r in result
            @test r.base === spec || r.base == spec
        end
    end

    @testset "AllostericMechanismSpec → empty (negative)" begin
        # Already-allosteric specs cannot be re-converted. The
        # _expand_to_allosteric specialization on AllostericMechanismSpec
        # returns an empty vector.
        m_seed = @allosteric_mechanism begin
            substrates: S; products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism(m_seed, uni_uni_allo)
        @test isempty(EnzymeRates._expand_to_allosteric(spec, uni_uni_allo))
    end

    @testset "oligomeric_state from reaction" begin
        # The catalytic_n of the result is taken from the reaction's
        # oligomeric_state, not hardcoded to 2.
        rxn4 = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            oligomeric_state: 4
        end
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism(m_seed, rxn4)
        result = EnzymeRates._expand_to_allosteric(spec, rxn4)
        @test !isempty(result)
        for r in result
            @test r.catalytic_n == 4
        end
    end

    @testset "Bi-bi sequential: 5 groups → 6 variants" begin
        # _expand_to_allosteric emits 1 baseline + 1 :OnlyR per group.
        # 5 groups → 1 + 5 = 6 variants.
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E_A
                E_A + B ⇌ E_A_B
                E + Q ⇌ E_Q
                E_Q + P ⇌ E_P_Q
                E_A_B <--> E_P_Q
            end
        end
        bi_bi_allo_rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            oligomeric_state: 2
        end
        spec = mechanism_spec_from_mechanism(m_seed, bi_bi_allo_rxn)
        result = EnzymeRates._expand_to_allosteric(spec, bi_bi_allo_rxn)
        @test length(result) == 6
        for r in result
            @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
            @test compile_mechanism(r) isa AllostericEnzymeMechanism
        end
    end
end
```

- [ ] **Step 8.3: Run tests; surface bugs per protocol; commit**

```bash
julia --project -e 'using Pkg; Pkg.test()'
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test(to_allosteric): rewrite with literal seeds; equivalence-style for uni-uni

Replaces the existing Allosteric conversion testset (1667-1723). 4 seeds:
uni-uni equivalence (4 variants), already-allosteric (negative),
oligomeric_state propagation, bi-bi sequential (6 variants).
EOF
)"
```

---

## Task 9: Section 4b — `_expand_add_allosteric_regulator`

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

**Pre-existing tests this absorbs:**
- `Add allosteric regulator` testset (lines 1725–1794).
- `Metabolite overlap: substrate as allosteric regulator` (lines 2233–2277).

- [ ] **Step 9.1: Delete the absorbed testsets**

- [ ] **Step 9.2: Add the rewritten testset**

Each sub-testset follows the Task 5 template: literal seed, round-trip validation, move call, checklist items 1–5. Below is the first seed in full as the in-task template; subsequent seeds list the seed mechanism, count derivation, and assertion-style.

```julia
# ─── _expand_add_allosteric_regulator ──────────────────────────────────
@testset "_expand_add_allosteric_regulator" begin

    @testset "Allosteric uni-uni + first allo regulator R: 3 variants" begin
        # SEED: uni-uni allosteric, all groups :EqualRT, no allosteric
        # regulator added yet.
        m_seed = @allosteric_mechanism begin
            substrates: S; products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism(m_seed, uni_uni_allo_reg)
        @test AllostericEnzymeMechanism(spec) === m_seed

        result = EnzymeRates._expand_add_allosteric_regulator(
            spec, uni_uni_allo_reg)

        # 1. count: R is the only un-added allosteric regulator. There
        # are 0 existing reg sites. The move enumerates non-:EqualRT
        # tags {:OnlyR, :OnlyT, :NonequalRT} × site options (new site
        # only, since 0 existing): 3 × 1 = 3 variants. The :EqualRT
        # branch is gated to "existing site with at least one non-:EqualRT
        # ligand" → not applicable here (no existing sites). → 3.
        @test length(result) == 3

        # 2. Δ params: cost of new R-binding K (+1) plus per-tag delta vs
        # :EqualRT base. :OnlyR/:OnlyT cheap → +1 total. :NonequalRT → +2
        # (K_R + K_T). So deltas across the 3 variants: [1, 1, 2].
        deltas = sort([r.n_fit_params_estimate -
                       spec.n_fit_params_estimate for r in result])
        @test deltas == [1, 1, 2]

        # 3. compilability
        for r in result
            @test compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. equivalence-style (N=3 ≤ 6).
        v_onlyR = @allosteric_mechanism begin
            substrates: S; products: P
            allosteric_regulators: R::OnlyR
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        v_onlyT = @allosteric_mechanism begin
            substrates: S; products: P
            allosteric_regulators: R::OnlyT
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        v_neq = @allosteric_mechanism begin
            substrates: S; products: P
            allosteric_regulators: R::NonequalRT
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        @test Set(compile_mechanism(r) for r in result) ==
            Set([v_onlyR, v_onlyT, v_neq])

        # 5. preservation
        for r in result
            @test r.catalytic_n == spec.catalytic_n
            @test r.group_tags == spec.group_tags
            @test r.base.reaction === spec.base.reaction
        end
    end

    # … (additional sub-testsets per seed listing below)
end
```

**Remaining seeds — apply the same template:**

**Seed 2 — Non-allosteric MechanismSpec → empty (negative):**

Plain `MechanismSpec` (uni-uni init, mechanism_spec_from_mechanism with `uni_uni_allo_reg`) → `_expand_add_allosteric_regulator` returns empty. The function specializes on `AllostericMechanismSpec`; non-allosteric input dispatches to the empty fallback.

**Seed 3 — Two regulators with site options (count = 7):**

Two-step process: first add R1 (3 variants from seed 1), pick the :OnlyR variant. Then call `_expand_add_allosteric_regulator` again to add R2.

```julia
# Count derivation: 3 non-:EqualRT tag flavors × 2 site options (new site
# OR R1's existing site) = 6. Plus 1 variant for :EqualRT at R1's
# (non-:EqualRT) existing site = 1. Total: 7.
@test length(r2_added) == 7
```

The :EqualRT-at-existing-site branch only fires because R1 is :OnlyR (non-:EqualRT). If R1 were already :EqualRT, that branch would not fire and the count would be 6.

**Seed 4 — EqualRT at existing reg site:**

Verify (via `findfirst`) that for a seed where R1 is :OnlyR at site 1, calling the move with R2 produces at least one variant where R2 is :EqualRT at site 1 (same site as R1). Property-style assertion on the result list.

**Seed 5 — Substrate-as-allosteric-regulator overlap:**

Reaction declares both substrate :S and allosteric_regulator :S. Seed is uni-uni allosteric. The move adds :S as an allosteric ligand at a new site; spec's allosteric_reg_sites must include `[:S]`. Compile verifies the dual-role metabolite shows up in both substrates and allosteric ligands. Property-style: at least one result has `:S in any(site for site in r.allosteric_reg_sites)`. Independent count derivation: 3 tag variants × 1 site option = 3 variants.

- [ ] **Step 9.3: Run tests; surface bugs per protocol; commit**

```bash
julia --project -e 'using Pkg; Pkg.test()'
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test(add_allosteric_regulator): rewrite with literal seeds; absorb overlap

Replaces Add allosteric regulator (1725-1794) and Metabolite overlap:
substrate as allosteric regulator (2233-2277). 5 seeds with independent
derivations.
EOF
)"
```

---

## Task 10: Section 4c — `_expand_change_allo_state`

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

**Pre-existing tests this absorbs:**
- `Remove TR equivalence` testset (lines 1796–1879).

- [ ] **Step 10.1: Delete the absorbed testset**

- [ ] **Step 10.2: Add the rewritten testset**

The first seed below is in full as the in-task template; subsequent seeds list the seed and derivation.

```julia
# ─── _expand_change_allo_state ─────────────────────────────────────────
@testset "_expand_change_allo_state" begin

    @testset "Allosteric uni-uni all-:EqualRT: 3 group-tag relaxations" begin
        # SEED: uni-uni allosteric with all 3 groups tagged :EqualRT.
        # Each non-:NonequalRT entry contributes ONE relaxation variant
        # (flip its value to :NonequalRT in the dense Dict).
        m_seed = @allosteric_mechanism begin
            substrates: S; products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism(m_seed, uni_uni_allo)
        @test AllostericEnzymeMechanism(spec) === m_seed

        result = EnzymeRates._expand_change_allo_state(spec, uni_uni_allo)

        # 1. count: 3 group_tags entries non-:NonequalRT (one per group,
        # all :EqualRT) + 0 reg_ligand_tags entries (no regulators) →
        # 3 relaxation variants.
        @test length(result) == 3

        # 2. Δ params: each removal converts :EqualRT → :NonequalRT.
        # _allo_state_delta(:EqualRT, :NonequalRT, is_re):
        #   For RE (is_re=true): factor 1 × (cost(NonequalRT) - cost(EqualRT))
        #   = 1 × (2 - 1) = +1.
        # For SS iso group (is_re=false): factor 2 × (2 - 1) = +2.
        # Two RE binding groups → +1 each (2 variants). One SS iso group →
        # +2 (1 variant). Deltas: [1, 1, 2].
        deltas = sort([r.n_fit_params_estimate -
                       spec.n_fit_params_estimate for r in result])
        @test deltas == [1, 1, 2]

        # 3. compilability
        for r in result
            @test compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. property-style: in each result, exactly one group's tag flipped
        # from non-:NonequalRT to :NonequalRT (relaxation move).
        for r in result
            relaxed = [g for g in keys(spec.group_tags)
                       if spec.group_tags[g] != :NonequalRT &&
                          r.group_tags[g] == :NonequalRT]
            @test length(relaxed) == 1
        end

        # 5. preservation: catalytic_n, reg sites, base.reaction unchanged.
        for r in result
            @test r.catalytic_n == spec.catalytic_n
            @test r.allosteric_reg_sites == spec.allosteric_reg_sites
            @test r.base.reaction === spec.base.reaction
        end
    end

    # … (additional sub-testsets per seed listing below)
end
```

**Remaining seeds:**

**Seed 2 — Fully relaxed → empty (negative):**

Under dense storage, "fully relaxed" means every group_tag is `:NonequalRT`
and every reg_ligand_tag is `:NonequalRT`. Construct the seed directly via
`@allosteric_mechanism` with every step `:: NonequalRT` and no allosteric
regulators — no eligible entries → empty result.

```julia
@testset "Fully relaxed → empty (negative)" begin
    # SEED: every group_tag and reg_ligand_tag is :NonequalRT.
    m_seed = @allosteric_mechanism begin
        substrates: S
        products: P
        site(:catalytic, 2): begin
            steps: begin
                E + P ⇌ E_P    :: NonequalRT
                E + S ⇌ E_S    :: NonequalRT
                E_S <--> E_P   :: NonequalRT
            end
        end
    end
    spec = allosteric_spec_from_mechanism(m_seed, uni_uni_allo)
    # All group_tags == :NonequalRT, no reg_ligand_tags → no eligible
    # entries → empty result.
    @test isempty(EnzymeRates._expand_change_allo_state(spec, uni_uni_allo))
end
```

**Seed 3 — MechanismSpec → empty:**

Plain `MechanismSpec` (uni-uni init, no allosteric conversion) → `_expand_change_allo_state` returns empty. The function dispatches to the MechanismSpec specialization at line 2046 of `src/mechanism_enumeration.jl` which returns empty.

**Seed 4 — Allosteric regulator tag removal delta:**

Seed: uni-uni allosteric with one regulator R tagged :OnlyR.

```julia
m_seed = @allosteric_mechanism begin
    substrates: S; products: P
    allosteric_regulators: R::OnlyR
    site(:catalytic, 2): begin
        steps: begin
            E + P ⇌ E_P    :: EqualRT
            E + S ⇌ E_S    :: EqualRT
            E_S <--> E_P   :: EqualRT
        end
    end
end
```

Move yields `(3 group-tag removals) + (1 reg-ligand-tag removal)` = 4 variants.

For the R-removal variant: `_allo_lig_state_delta(:OnlyR, :NonequalRT)` = cost(:NonequalRT) - cost(:OnlyR) = 2 - 1 = +1. Filter for the reg-ligand-removal variant; assert delta = +1.

- [ ] **Step 10.3: Run tests; surface bugs per protocol; commit**

```bash
julia --project -e 'using Pkg; Pkg.test()'
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test(change_allo_state): rewrite with literal seeds and 7-item checklist

Replaces Remove TR equivalence (1796-1879). 4 seeds with independent
derivations including the MechanismSpec specialization (yields empty).
EOF
)"
```

---

## Task 11: Section 5 — composition (dedup! and expand_mechanisms)

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

**Pre-existing tests this absorbs:**
- `Dedup` testset (lines 1881–1937).
- `expand_mechanisms` testset (lines 1939–2018).

`dedup!` is already tested via direct `MechanismSpec` literals (no spec-from-init brittleness) — preserve as-is, just reorganize under the section header. `expand_mechanisms` testset uses `init_mechanisms |> first` patterns; replace with literal seeds.

- [ ] **Step 11.1: Insert section 5 header**

```julia
# ═══════════════════════════════════════════════════════════════════════
# 5. Composition (dedup!, expand_mechanisms)
# ═══════════════════════════════════════════════════════════════════════
```

- [ ] **Step 11.2: Move existing `Dedup` testset under the section header**

No changes to assertions.

- [ ] **Step 11.3: Rewrite `expand_mechanisms` testset to use literal seeds**

For each existing sub-testset (`Returns dict keyed by param count`, `Allosteric expansion included`, `No self-expansion to same param count`, `Allosteric rewrap preserves structure`, `Dead-end excludes allosteric regs`), replace `EnzymeRates.init_mechanisms(rxn) |> first` with a literal `@enzyme_mechanism` seed + `mechanism_spec_from_mechanism` round-trip. Keep the assertions intact.

- [ ] **Step 11.4: Run tests; commit**

```bash
julia --project -e 'using Pkg; Pkg.test()'
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test(composition): reorganize dedup! and expand_mechanisms under section 5

dedup! tests preserved verbatim. expand_mechanisms tests replace
init_mechanisms |> first patterns with literal @enzyme_mechanism seeds.
EOF
)"
```

---

## Task 12: Section 6 — integration (enumerate_all)

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

**Pre-existing tests this absorbs:**
- `Integration` testset (lines 2020–2118).

Integration tests are coordination-layer; preserve substantively, just move under section 6 header. The `enumerate_all` helper at line 165 is fine to keep.

- [ ] **Step 12.1: Insert section 6 header**

```julia
# ═══════════════════════════════════════════════════════════════════════
# 6. Integration (enumerate_all)
# ═══════════════════════════════════════════════════════════════════════
```

- [ ] **Step 12.2: Move the `enumerate_all` helper and `Integration` testset under section 6**

Preserve assertions verbatim.

- [ ] **Step 12.3: Run tests; commit**

```bash
julia --project -e 'using Pkg; Pkg.test()'
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test(integration): place enumerate_all tests under section 6

Reorganization only — no assertion changes.
EOF
)"
```

---

## Task 13 (optional, confirmation-gated): Move out-of-scope testsets

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (delete three testsets)
- Modify: `test/test_types.jl`, `test/test_dsl.jl`, `test/test_rate_eq_derivation.jl` (add the moved testsets)

This task is OPTIONAL and gated on Denis's per-row confirmation. Do NOT execute it without explicit approval per testset.

Three testsets are out of scope for the enumeration pipeline file:

1. `AllostericEnzymeMechanism TR equivalence` (lines 217–243): tests `cat_allo_state` accessor → belongs in `test_types.jl`.
2. `test reaction atom balance` (lines 198–215): tests `@enzyme_reaction` macro → belongs in `test_dsl.jl`.
3. `Tagged groups exclude T-state params` (lines 2120–2201): tests `parameters(m)` and the canonicalizer → belongs in `test_rate_eq_derivation.jl`.

- [ ] **Step 13.1: Ask Denis per-row whether to move each one**

Surface a question:

```
Per design §4 "What's absorbed/moved out", three testsets in
test_mechanism_enumeration.jl test something other than the enumeration
pipeline. Move each? Per row, please:

1. AllostericEnzymeMechanism TR equivalence → test_types.jl  [yes/no/defer]
2. test reaction atom balance → test_dsl.jl                   [yes/no/defer]
3. Tagged groups exclude T-state params → test_rate_eq_derivation.jl [yes/no/defer]
```

- [ ] **Step 13.2: For each "yes", cut from `test_mechanism_enumeration.jl` and paste into the target file**

If the target file does not have a `using EnzymeRates` import that covers the moved testset's symbol needs, add the necessary `using` lines. Run tests.

- [ ] **Step 13.3: Commit per-row**

One commit per testset moved. Commit message format:

```
test: move <testset name> from test_mechanism_enumeration.jl to <target>

<testset> tests <thing-tested> — out of scope for the enumeration
pipeline file. No assertion changes.
```

If Denis declines all three, skip this task entirely — no commit.

---

## Self-Review Checklist

After completing all tasks, run a final self-review:

- [ ] **Spec coverage:** every section of the design doc maps to at least one task. Specifically:
  - §2 helper added (Task 1). ✓
  - §3 checklist applied to every spec-consuming move (Tasks 5–10). ✓
  - §4 file structure realized (Tasks 3–12 each operate on the right section). ✓
  - §5 execution plan: 11 commits committed plus optional Task 13. ✓
  - §6 bug-handling protocol invoked at each "Run tests" step. ✓

- [ ] **Independent-derivation rule:** every count, delta, and equivalence-set entry across Tasks 5–10 has a derivation comment explaining WHY the prediction is what it is (not just "the code returned this number"). Spot-check 3–5 random testsets.

- [ ] **No `init_mechanisms |> first` in spec-consuming testsets:** grep the rewritten file:

  ```bash
  grep -n "init_mechanisms" test/test_mechanism_enumeration.jl | \
      grep -v "init_mechanisms\""    # filter out the testset header
  ```

  Acceptable remaining uses: inside `@testset "init_mechanisms"` itself, inside `@testset "expand_mechanisms"` for input batteries, and inside the integration tests. NOT acceptable: any base-spec or allosteric move test using `init_mechanisms |> first` for the SEED.

- [ ] **Compilability check:** every per-move test asserts `compile_mechanism(r) isa <Type>` for at least one r in result.

- [ ] **Final test run:**

  ```bash
  julia --project -e 'using Pkg; Pkg.test()'
  ```

  Green. If any test fails, follow design §6 protocol — do not weaken or skip.

- [ ] **Commit count check:** between Task 1's commit and Task 12's commit, expect 11–12 commits (plus 0–3 from Task 13). Each commit is reviewable on its own.
