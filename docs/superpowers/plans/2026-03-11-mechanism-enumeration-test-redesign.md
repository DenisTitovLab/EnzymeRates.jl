# Mechanism Enumeration Test Redesign — Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development
> (if subagents available) or superpowers:executing-plans to implement this plan.
> Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rewrite mechanism enumeration tests with hand-calculated expected values,
per-stage isolation, and compile_mechanism round-trip verification — exposing known
pipeline bugs via `@test_broken`.

**Architecture:** Two test files: `mechanism_enumeration_test_specs.jl` (reaction
definitions, helper functions, `mechanism_spec_from_mechanism`) and
`test_mechanism_enumeration.jl` (all stage testsets). Each stage testset defines
mechanisms via `@enzyme_mechanism`, converts to `MechanismSpec`, verifies round-trip,
then runs the stage function and checks counts/properties. Known bugs produce
`@test_broken` with comments explaining expected vs actual behavior.

**Tech Stack:** Julia, Test stdlib, EnzymeRates.jl internals (`_catalytic_topologies`,
`_expand_ress_variants`, `_expand_dead_end_inhibitors`, `_expand_equivalence_constraints`,
`_deduplicate`, `_expand_allosteric`, `_expand_tr_equivalence`,
`_deduplicate_allosteric`, `compile_mechanism`, `enumerate_mechanisms`)

**Spec:** `docs/superpowers/specs/2026-03-10-mechanism-enumeration-test-redesign.md`

---

## Important context for implementer

### Key files
- `src/mechanism_enumeration.jl` — All pipeline stage functions and types
  (`MechanismSpec`, `AllostericMechanismSpec`, `compile_mechanism`)
- `src/types.jl` — `EnzymeReaction`, `EnzymeMechanism`, `AllostericEnzymeMechanism`
- `src/dsl.jl` — `@enzyme_reaction`, `@enzyme_mechanism` macros
- `test/mechanism_enumeration_test_specs.jl` — Will be rewritten (reactions, helpers)
- `test/test_mechanism_enumeration.jl` — Will be rewritten (all test logic)
- `test/runtests.jl` — Includes test files; do NOT modify

### How to run tests
```bash
# Full suite (slow, ~7 min):
julia --project -e 'using Pkg; Pkg.test()'

# Just mechanism enumeration (faster, interactive):
julia --project -e '
    using Test, EnzymeRates, Random
    include("test/mechanism_enumeration_test_specs.jl")
    include("test/test_mechanism_enumeration.jl")
'
```

### How `mechanism_spec_from_mechanism` works
Converts a compiled `EnzymeMechanism` back to a `MechanismSpec` by:
1. Calling `enumerate_enzyme_forms(rxn)` to get the form index mapping
2. Walking the mechanism's reactions to map enzyme form names → indices
3. Extracting edges, equilibrium_steps, param_constraints from the mechanism
Use: `spec = mechanism_spec_from_mechanism(m, rxn)` where `m` was defined via
`@enzyme_mechanism` and `rxn` via `@enzyme_reaction`.

### Known bugs (tests should expose these)
1. **RE/SS expansion** (`_expand_ress_variants`): Only toggles binding edges, keeps
   isomerization always SS. Should toggle ALL steps. Current: `2^n_binding - 1`.
   Correct: `2^n_total - 1`. Mark correct expectations as `@test_broken`.
2. **Missing stage 2.5**: No substrate/product dead-end expansion exists.
   All stage 2.5 tests are `@test_broken`.
3. **Regulator dead-end** (`_expand_dead_end_inhibitors`): Allows inhibitors to bind
   ALL catalytic forms. Should only allow binding to forms not fully occupied (not all
   substrates bound AND not all products bound). Mark correct expectations `@test_broken`.

### Style rules (from CLAUDE.md)
- 92-char line limit, 4-space indentation
- All files start with 2-line `ABOUTME:` comment
- Match existing code style
- No unnecessary comments about changes or improvements

---

## Chunk 1: Foundation — reactions and helpers

### Task 1: Rewrite reaction definitions in mechanism_enumeration_test_specs.jl

**Files:**
- Modify: `test/mechanism_enumeration_test_specs.jl` (full rewrite)

This task replaces the entire contents of `mechanism_enumeration_test_specs.jl` with
clean reaction definitions and the `mechanism_spec_from_mechanism` helper. All the
old `StageExpansionTestSpec`, `EnumerationTestSpec`, and builder functions are removed.

- [ ] **Step 1: Read the current file to understand what's being replaced**

Read `test/mechanism_enumeration_test_specs.jl` in full.

- [ ] **Step 2: Write the new file with all reaction definitions**

Replace entire contents with:

```julia
# ABOUTME: Reaction definitions and helpers for mechanism enumeration tests
# ABOUTME: Defines test reactions, mechanism_spec_from_mechanism, and run helpers

using Random

# ── Helper: convert EnzymeMechanism → MechanismSpec ──────────

# [Keep the existing mechanism_spec_from_mechanism function exactly as-is,
#  lines 8-90 of the current file. Do NOT modify it.]

# ── Reaction definitions ─────────────────────────────────────
# Organized by complexity. Regulated reactions use R for allosteric,
# I/J for dead-end inhibitors.

# --- No regulators ---
const uni_uni = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
end

const uni_bi = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
end

const bi_bi = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end

const bi_bi_ping_pong = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
end

const ter_ter = @enzyme_reaction begin
    substrates: A[C], B[N], D[X]
    products: P[C], Q[N], R[X]
end

# --- Dead-end inhibitor variants ---
const uni_uni_dead_end_I = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    dead_end_inhibitors: I
end

const uni_bi_dead_end_I = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    dead_end_inhibitors: I
end

const bi_bi_dead_end_I = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    dead_end_inhibitors: I
end

const bi_bi_ping_pong_dead_end_I = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    dead_end_inhibitors: I
end

const uni_uni_dead_end_I_J = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    dead_end_inhibitors: I, J
end

# --- Allosteric regulator variants ---
const uni_uni_allosteric_R = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: R
end

const uni_bi_allosteric_R = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    allosteric_regulators: R
end

const bi_bi_ping_pong_allosteric_R = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    allosteric_regulators: R
end

const uni_bi_allosteric_R_cn2 = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    allosteric_regulators: R
end

# --- Mixed ---
const bi_bi_dead_end_I_allosteric_R = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    dead_end_inhibitors: I
    allosteric_regulators: R
end

const bi_bi_allosteric_R1_R2 = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    allosteric_regulators: R1, R2
end

# --- Unknown role (for end-to-end partitioning) ---
const uni_uni_reg_unknown = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    regulators: I
end

const uni_bi_reg_unknown = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    regulators: I
end

const bi_bi_ping_pong_reg_unknown = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    regulators: I
end

# ── Helper: run pipeline stage by stage (all partitions) ─────

"""
Run the full enumeration pipeline stage by stage, summing
counts across all regulator partitions. Returns a NamedTuple
with per-stage totals.
"""
function _run_full_pipeline_stages(rxn; catalytic_n::Int=0,
                                   max_re_groups::Int=7)
    catalytic = EnzymeRates._catalytic_topologies(rxn)

    roles = EnzymeRates.regulator_roles(rxn)
    fixed_de = Symbol[r[1] for r in roles
                       if r[2] == :dead_end]
    fixed_allo = Symbol[r[1] for r in roles
                         if r[2] == :allosteric]
    unknown = Symbol[r[1] for r in roles
                      if r[2] == :unknown]
    n_unknown = length(unknown)

    n_ress = 0; n_de = 0; n_eq = 0; n_dd = 0
    n_allo = 0; n_tr = 0; n_allo_dd = 0

    for reg_mask in 0:(1 << n_unknown) - 1
        de = Symbol[fixed_de;
            [unknown[i] for i in 1:n_unknown
             if (reg_mask >> (i - 1)) & 1 == 0]]
        al = Symbol[fixed_allo;
            [unknown[i] for i in 1:n_unknown
             if (reg_mask >> (i - 1)) & 1 == 1]]

        ress = EnzymeRates._expand_ress_variants(
            catalytic, rxn; max_re_groups)
        n_ress += length(ress)

        with_de = EnzymeRates._expand_dead_end_inhibitors(
            ress, rxn; dead_end_regs=de)
        n_de += length(with_de)

        with_eq =
            EnzymeRates._expand_equivalence_constraints(
                with_de, rxn)
        n_eq += length(with_eq)

        dd = EnzymeRates._deduplicate(with_eq, rxn)
        n_dd += length(dd)

        if !isempty(al)
            cn = catalytic_n > 0 ? catalytic_n : 1
            allo = EnzymeRates._expand_allosteric(
                dd, rxn; catalytic_n=cn,
                allosteric_regs=al)
            n_allo += length(allo)
            allo = EnzymeRates._expand_tr_equivalence(
                allo, rxn)
            n_tr += length(allo)
            allo = EnzymeRates._deduplicate_allosteric(
                allo, rxn)
            n_allo_dd += length(allo)
        end
    end

    (catalytic=length(catalytic), ress=n_ress,
     dead_end=n_de, equivalence=n_eq, dedup=n_dd,
     allosteric=n_allo, tr_equiv=n_tr,
     allosteric_dedup=n_allo_dd)
end
```

- [ ] **Step 3: Verify the file loads without errors**

Run:
```bash
julia --project -e '
    using Test, EnzymeRates, Random
    include("test/mechanism_enumeration_test_specs.jl")
    println("All reactions defined successfully")
    println("Reactions: uni_uni, uni_bi, bi_bi, bi_bi_ping_pong, ter_ter")
'
```
Expected: prints success message, no errors.

- [ ] **Step 4: Commit**

```bash
git add test/mechanism_enumeration_test_specs.jl
git commit -m "Rewrite mechanism_enumeration_test_specs.jl with clean reaction definitions"
```

---

## Chunk 2: Stage 1 and Stage 2 tests

### Task 2: Stage 1 — Catalytic topologies

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (full rewrite, starting fresh)

- [ ] **Step 1: Write the Stage 1 testset**

Replace entire contents of `test/test_mechanism_enumeration.jl` with Stage 1 tests.
Define ALL expected mechanisms for uni_uni (1) and uni_bi (3) via `@enzyme_mechanism`.
Define 3-4 representative bi_bi mechanisms. Define the ping-pong mechanism.
Verify round-trip for each hand-defined mechanism.

The file structure:
```julia
# ABOUTME: Tests for the staged mechanism enumeration pipeline
# ABOUTME: Organized by stage with hand-calculated expected values

@testset "Mechanism Enumeration Pipeline" begin

@testset "Stage 1: Catalytic topologies" begin
    # --- Uni-Uni: 1 topology ---
    @testset "Uni-Uni" begin
        topos = EnzymeRates._catalytic_topologies(uni_uni)
        @test length(topos) == 1

        # The single Uni-Uni mechanism: E⇌ES⇌EP⇌E
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, ES[C], EP[C]
            end
            steps: begin
                [E, S] ⇌ [ES]
                [E, P] ⇌ [EP]
                [ES] <--> [EP]
            end
        end
        spec_uu = mechanism_spec_from_mechanism(
            m_uu, uni_uni)
        # Round-trip
        @test compile_mechanism(spec_uu) === m_uu

        # Properties
        for t in topos
            @test t.n_catalytic_edges == length(t.edges)
            @test count(.!t.equilibrium_steps) >= 1
        end
    end

    # --- Uni-Bi: 3 topologies ---
    @testset "Uni-Bi" begin
        topos = EnzymeRates._catalytic_topologies(uni_bi)
        @test length(topos) == 3

        # Define all 3 expected mechanisms by hand:
        # 1. Sequential: release P first, then Q
        # 2. Sequential: release Q first, then P
        # 3. Random: both P and Q can release from EPQ

        # [Implementer: define each mechanism using
        #  @enzyme_mechanism with the correct enzyme forms
        #  and steps. Use enumerate_enzyme_forms(uni_bi)
        #  to identify the form names. Verify round-trip
        #  for each.]

        # Properties
        for t in topos
            @test t.n_catalytic_edges == length(t.edges)
            @test count(.!t.equilibrium_steps) >= 1
        end
    end

    # --- Bi-Bi: 9 topologies ---
    @testset "Bi-Bi" begin
        topos = EnzymeRates._catalytic_topologies(bi_bi)
        @test length(topos) == 9

        # Define representative mechanisms (one per category):
        # 1. Ordered sequential: A binds first, Q released last
        # 2. Random binding, sequential release
        # 3. Sequential binding, random release
        # 4. Random binding, random release

        # [Implementer: define 4 representative mechanisms.
        #  Verify round-trip for each. Check that each
        #  appears in the topology list.]

        for t in topos
            @test t.n_catalytic_edges == length(t.edges)
            @test count(.!t.equilibrium_steps) >= 1
        end
    end

    # --- Bi-Bi Ping-Pong: 10 topologies ---
    @testset "Bi-Bi Ping-Pong" begin
        topos = EnzymeRates._catalytic_topologies(
            bi_bi_ping_pong)
        @test length(topos) == 10

        # Define the ping-pong mechanism specifically
        # (the one topology not present in Bi-Bi)
        # [Implementer: define the ping-pong mechanism
        #  with E_X intermediate form]

        for t in topos
            @test t.n_catalytic_edges == length(t.edges)
            @test count(.!t.equilibrium_steps) >= 1
        end
    end

    # --- Ter-Ter: 49 topologies ---
    @testset "Ter-Ter" begin
        topos = EnzymeRates._catalytic_topologies(ter_ter)
        @test length(topos) == 49

        for t in topos
            @test t.n_catalytic_edges == length(t.edges)
            @test count(.!t.equilibrium_steps) >= 1
        end
    end
end

end # outer testset
```

**Critical**: The implementer MUST:
1. Run `EnzymeRates.enumerate_enzyme_forms(rxn)` for each reaction to learn the
   exact form names (e.g., `:E_0_0_0_0`, `:E_A_0_0_0`) before writing mechanisms
2. Check that the `@enzyme_mechanism` DSL uses the correct atom annotations
   matching the reaction (e.g., `ES[C]` for S[C])
3. Verify each round-trip: `compile_mechanism(spec) === m`

- [ ] **Step 2: Run tests to verify Stage 1**

```bash
julia --project -e '
    using Test, EnzymeRates, Random
    include("test/mechanism_enumeration_test_specs.jl")
    include("test/test_mechanism_enumeration.jl")
'
```
Expected: Stage 1 tests pass (catalytic topology counts should match current code).
If any counts don't match, investigate — these are independently calculated.

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Add Stage 1 catalytic topology tests with hand-defined mechanisms"
```

### Task 3: Stage 2 — RE/SS expansion

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (add Stage 2 testset)

- [ ] **Step 1: Add Stage 2 testset after Stage 1**

Add inside the outer `@testset "Mechanism Enumeration Pipeline"`:

```julia
@testset "Stage 2: RE/SS expansion" begin
    # Use the Uni-Uni mechanism from Stage 1 tests
    @testset "Uni-Uni: 2^3 - 1 = 7 variants" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, ES[C], EP[C]
            end
            steps: begin
                [E, S] ⇌ [ES]
                [E, P] ⇌ [EP]
                [ES] <--> [EP]
            end
        end
        spec = mechanism_spec_from_mechanism(m, uni_uni)
        result = EnzymeRates._expand_ress_variants(
            [spec], uni_uni)
        # Correct: all 3 steps toggleable, 2^3 - 1 = 7
        # BUG: current code only toggles binding edges
        #   (2 binding edges), gives 2^2 - 1 = 3
        @test_broken length(result) == 7
        # Document current (buggy) behavior:
        @test length(result) == 3

        # Property: every output has at least 1 SS step
        for s in result
            @test any(.!s.equilibrium_steps)
        end
        # Property: edges unchanged
        for s in result
            @test s.edges == spec.edges
            @test s.n_catalytic_edges ==
                spec.n_catalytic_edges
        end
    end

    @testset "Uni-Bi: 2^n - 1 variants" begin
        topo = EnzymeRates._catalytic_topologies(uni_bi)[1]
        n_steps = length(topo.edges)
        result = EnzymeRates._expand_ress_variants(
            [topo], uni_bi)
        # [Implementer: count n_steps, compute 2^n-1,
        #  write @test_broken for correct and @test for
        #  current behavior. Current code gives
        #  2^n_binding - 1.]

        for s in result
            @test any(.!s.equilibrium_steps)
        end
    end

    @testset "Bi-Bi: 2^n - 1 variants" begin
        topo = EnzymeRates._catalytic_topologies(bi_bi)[1]
        n_steps = length(topo.edges)
        result = EnzymeRates._expand_ress_variants(
            [topo], bi_bi)
        @test_broken length(result) == 2^n_steps - 1
        # [Implementer: verify current buggy count]

        for s in result
            @test any(.!s.equilibrium_steps)
        end
    end

    @testset "max_re_groups filtering" begin
        topo = EnzymeRates._catalytic_topologies(
            ter_ter)[1]
        n_steps = length(topo.edges)
        # With default max_re_groups=7, some assignments
        # filtered (all-SS creates n_forms groups > 7)
        result = EnzymeRates._expand_ress_variants(
            [topo], ter_ter)
        @test length(result) < 2^n_steps - 1

        # The all-SS assignment should be excluded
        all_ss = all(s -> all(.!s.equilibrium_steps),
            result)
        @test !all_ss  # no output should be all-SS
    end
end
```

- [ ] **Step 2: Run tests**

```bash
julia --project -e '
    using Test, EnzymeRates, Random
    include("test/mechanism_enumeration_test_specs.jl")
    include("test/test_mechanism_enumeration.jl")
'
```
Expected: `@test_broken` tests should show as broken (expected). `@test` should pass.

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Add Stage 2 RE/SS expansion tests exposing isomerization toggle bug"
```

---

## Chunk 3: Stages 2.5, 3, 4

### Task 4: Stage 2.5 — Substrate/product dead-end expansion (future)

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (add testset)

- [ ] **Step 1: Compute off-cycle form counts**

Run this to determine exact off-cycle form counts for each reaction:
```bash
julia --project -e '
    using EnzymeRates
    for (name, rxn) in [
        ("uni_uni", @enzyme_reaction(begin
            substrates: S[C]; products: P[C] end)),
        ("bi_bi", @enzyme_reaction(begin
            substrates: A[C], B[N]; products: P[C], Q[N] end)),
        ("bi_bi_pp", @enzyme_reaction(begin
            substrates: A[CX], B[N]; products: P[C], Q[NX] end)),
    ]
        _, forms = EnzymeRates.enumerate_enzyme_forms(rxn)
        topos = EnzymeRates._catalytic_topologies(rxn)
        # Collect all on-cycle form indices across topologies
        on_cycle = Set{Int}()
        for t in topos
            for (i, j) in t.edges
                push!(on_cycle, i, j)
            end
        end
        n_off = length(forms) - length(on_cycle)
        println("$name: $(length(forms)) forms, " *
            "$(length(on_cycle)) on-cycle, $n_off off-cycle")
    end
'
```

Use these counts to fill in exact expected values.

- [ ] **Step 2: Add Stage 2.5 testset — all @test_broken**

```julia
@testset "Stage 2.5: Substrate/product dead-end expansion" begin
    # This stage does not exist yet. All tests are
    # @test_broken to document the expected behavior
    # for when the stage is implemented.

    @testset "Uni-Uni passthrough (no off-cycle forms)" begin
        # Uni-Uni: all 3 forms (E, ES, EP) are on the
        # catalytic cycle → 0 off-cycle → passthrough
        # [Implementer: @test_broken that the future
        #  stage function returns input unchanged]
    end

    @testset "Bi-Bi: 4 off-cycle forms" begin
        # Forms EAP, EBP, EAQ, EBQ are off-cycle
        # Expected: 2^4 = 16 variants per input spec
        # [Implementer: document specific off-cycle forms
        #  and expected expansion count]
    end

    @testset "Bi-Bi Ping-Pong" begin
        # [Implementer: compute off-cycle count from
        #  step 1 above, document expected expansion]
    end
end
```

- [ ] **Step 3: Run and commit**

```bash
julia --project -e '
    using Test, EnzymeRates, Random
    include("test/mechanism_enumeration_test_specs.jl")
    include("test/test_mechanism_enumeration.jl")
'
git add test/test_mechanism_enumeration.jl
git commit -m "Add Stage 2.5 substrate/product dead-end tests (all test_broken)"
```

### Task 5: Stage 3 — Regulator dead-end expansion

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (add testset)

- [ ] **Step 1: Add Stage 3 testset**

```julia
@testset "Stage 3: Regulator dead-end expansion" begin
    @testset "No regulators: passthrough" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni)[1]
        result = EnzymeRates._expand_dead_end_inhibitors(
            [topo], uni_uni; dead_end_regs=Symbol[])
        @test length(result) == 1
    end

    @testset "Uni-Uni + I: only E is eligible" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_dead_end_I)[1]
        result = EnzymeRates._expand_dead_end_inhibitors(
            [topo], uni_uni_dead_end_I;
            dead_end_regs=[:I])
        # Correct: only E eligible (ES, EP fully occupied)
        # → (2^1)^1 = 2 variants (E, EI)
        # BUG: current code allows I to bind all 3 forms
        # → (2^1)^3 = 8
        @test_broken length(result) == 2
        @test length(result) == 8
    end

    @testset "Uni-Bi + I: E, E_P, E_Q eligible" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_bi_dead_end_I)[1]
        result = EnzymeRates._expand_dead_end_inhibitors(
            [topo], uni_bi_dead_end_I;
            dead_end_regs=[:I])
        # Correct: 3 eligible forms → (2^1)^3 = 8
        # BUG: current code uses all catalytic forms
        # [Implementer: count current buggy output and
        #  add @test for it]
        @test_broken length(result) == 8
    end

    @testset "Bi-Bi + I: 5 eligible forms" begin
        topo = EnzymeRates._catalytic_topologies(
            bi_bi_dead_end_I)[1]
        result = EnzymeRates._expand_dead_end_inhibitors(
            [topo], bi_bi_dead_end_I;
            dead_end_regs=[:I])
        # Correct: E, EA, EB, EP, EQ eligible (5 forms)
        # → (2^1)^5 = 32
        # [Implementer: determine current buggy count]
        @test_broken length(result) == 32
    end

    @testset "2 inhibitors: Uni-Uni + I, J" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_dead_end_I_J)[1]
        result = EnzymeRates._expand_dead_end_inhibitors(
            [topo], uni_uni_dead_end_I_J;
            dead_end_regs=[:I, :J])
        # Correct: only E eligible → (2^2)^1 = 4
        # BUG: current code uses all forms
        @test_broken length(result) == 4
    end

    @testset "Allosteric-only: passthrough" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_allosteric_R)[1]
        result = EnzymeRates._expand_dead_end_inhibitors(
            [topo], uni_uni_allosteric_R;
            dead_end_regs=Symbol[])
        @test length(result) == 1
    end

    # Properties for all non-passthrough results
    @testset "Properties" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_dead_end_I)[1]
        result = EnzymeRates._expand_dead_end_inhibitors(
            [topo], uni_uni_dead_end_I;
            dead_end_regs=[:I])
        for s in result
            @test length(s.edges) >= s.n_catalytic_edges
            @test s.n_catalytic_edges ==
                topo.n_catalytic_edges
        end
    end
end
```

- [ ] **Step 2: Run and commit**

Run tests, verify `@test_broken` show as broken and `@test` pass.
```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Add Stage 3 regulator dead-end tests exposing eligible-form bug"
```

### Task 6: Stage 4 — Equivalence constraints

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (add testset)

- [ ] **Step 1: Add Stage 4 testset**

```julia
@testset "Stage 4: Equivalence constraints" begin
    @testset "No equiv groups: passthrough" begin
        # Uni-Uni: S and P bind different atoms
        topo = EnzymeRates._catalytic_topologies(
            uni_uni)[1]
        result =
            EnzymeRates._expand_equivalence_constraints(
                [topo], uni_uni)
        @test length(result) == 1
    end

    @testset "1 equiv group from dead-end inhibitor" begin
        # Uni-Uni + dead-end I: if I binds to E (→EI),
        # and the same I binding appears in another edge,
        # they form an equivalence group.
        # [Implementer: build a mechanism with 2 I-binding
        #  edges via mechanism_spec_from_mechanism. Check
        #  that equivalence constraints expand to 2
        #  variants (unconstrained + constrained)]
    end

    @testset "Same metabolite as substrate and regulator" begin
        # Edge case: metabolite S is both a substrate and
        # a dead-end regulator. Substrate-binding edges
        # and regulator-binding edges of S must NOT be
        # in the same equivalence group.
        # [Implementer: create this scenario and verify
        #  that the groups are separated]
    end

    @testset "Properties" begin
        # Output count ≥ input count
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_dead_end_I)[1]
        de = EnzymeRates._expand_dead_end_inhibitors(
            [topo], uni_uni_dead_end_I;
            dead_end_regs=[:I])
        eq = EnzymeRates._expand_equivalence_constraints(
            de, uni_uni_dead_end_I)
        @test length(eq) >= length(de)

        # Constrained variants have lower param_count
        unconstrained = filter(
            s -> isempty(s.param_constraints), eq)
        constrained = filter(
            s -> !isempty(s.param_constraints), eq)
        if !isempty(constrained) && !isempty(unconstrained)
            @test minimum(s.param_count
                for s in constrained) <
                maximum(s.param_count
                    for s in unconstrained)
        end
    end
end
```

- [ ] **Step 2: Run and commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Add Stage 4 equivalence constraint tests"
```

---

## Chunk 4: Stages 5, 6, 7, 8

### Task 7: Stage 5 — Deduplication

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (add testset)

- [ ] **Step 1: Add Stage 5 testset**

```julia
@testset "Stage 5: Deduplication" begin
    @testset "Uni-Uni: 7 variants dedup to 1" begin
        # All RE/SS variants of the single Uni-Uni
        # topology produce the same concentration
        # fingerprint {1, [S], [P]}.
        # BUG: current code only produces 3 variants
        # (can't test 7→1 until RE/SS bug is fixed)
        topo = EnzymeRates._catalytic_topologies(
            uni_uni)[1]
        ress = EnzymeRates._expand_ress_variants(
            [topo], uni_uni)
        deduped = EnzymeRates._deduplicate(
            ress, uni_uni)
        # Current behavior: 3 → 1 (still deduplicates
        # correctly within the buggy set)
        @test length(deduped) == 1
        @test deduped[1].param_count <=
            minimum(s.param_count for s in ress)
    end

    @testset "Bi-Bi single-SS: all 5 distinct" begin
        # Each single-SS variant of the ordered Bi-Bi
        # has a different concentration fingerprint.
        # After dedup, all 5 should survive.
        # [Implementer: get the first Bi-Bi topology,
        #  create 5 single-SS specs by hand using
        #  @enzyme_mechanism for each, run dedup,
        #  verify count == 5]
    end

    @testset "Duplicate removal" begin
        # Feed 2 identical specs → verify dedup
        # removes one
        topo = EnzymeRates._catalytic_topologies(
            uni_uni)[1]
        result = EnzymeRates._deduplicate(
            [topo, deepcopy(topo)], uni_uni)
        @test length(result) == 1
    end

    @testset "Keeps lower param_count" begin
        # Feed 2 specs with same fingerprint but
        # different param_count → verify the lower
        # param_count one is kept
        topo = EnzymeRates._catalytic_topologies(
            uni_uni)[1]
        topo2 = deepcopy(topo)
        # Artificially increase param_count
        topo2_modified = EnzymeRates.MechanismSpec(
            topo2.reaction, topo2.edges,
            topo2.n_catalytic_edges,
            topo2.equilibrium_steps,
            topo2.param_constraints,
            topo2.param_count + 5)
        result = EnzymeRates._deduplicate(
            [topo, topo2_modified], uni_uni)
        @test length(result) == 1
        @test result[1].param_count == topo.param_count
    end
end
```

- [ ] **Step 2: Run and commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Add Stage 5 deduplication tests with hand-verified cases"
```

### Task 8: Stage 6 — Allosteric expansion

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (add testset)

- [ ] **Step 1: Add Stage 6 testset**

```julia
@testset "Stage 6: Allosteric expansion" begin
    @testset "No allosteric regs: passthrough" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni)[1]
        dd = EnzymeRates._deduplicate([topo], uni_uni)
        result = EnzymeRates._expand_allosteric(
            dd, uni_uni; catalytic_n=1,
            allosteric_regs=Symbol[])
        @test length(result) == 1
    end

    @testset "1 reg, catalytic_n=1" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_allosteric_R)[1]
        dd = EnzymeRates._deduplicate(
            [topo], uni_uni_allosteric_R)
        result = EnzymeRates._expand_allosteric(
            dd, uni_uni_allosteric_R; catalytic_n=1,
            allosteric_regs=[:R])
        @test length(result) == 1  # m=1 only
    end

    @testset "1 reg, catalytic_n=2" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_bi_allosteric_R_cn2)[1]
        dd = EnzymeRates._deduplicate(
            [topo], uni_bi_allosteric_R_cn2)
        result = EnzymeRates._expand_allosteric(
            dd, uni_bi_allosteric_R_cn2; catalytic_n=2,
            allosteric_regs=[:R])
        @test length(result) == 2  # m=1, m=2
    end

    @testset "1 reg, catalytic_n=3" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_allosteric_R)[1]
        dd = EnzymeRates._deduplicate(
            [topo], uni_uni_allosteric_R)
        result = EnzymeRates._expand_allosteric(
            dd, uni_uni_allosteric_R; catalytic_n=3,
            allosteric_regs=[:R])
        @test length(result) == 3  # m=1, m=2, m=3
    end

    @testset "2 regs, catalytic_n=1" begin
        topo = EnzymeRates._catalytic_topologies(
            bi_bi_allosteric_R1_R2)[1]
        dd = EnzymeRates._deduplicate(
            [topo], bi_bi_allosteric_R1_R2)
        result = EnzymeRates._expand_allosteric(
            dd, bi_bi_allosteric_R1_R2; catalytic_n=1,
            allosteric_regs=[:R1, :R2])
        # 2 set partitions: {R1},{R2} and {R1,R2}
        # Each with m=1 only → 2 variants
        @test length(result) == 2
    end

    @testset "2 regs, catalytic_n=2" begin
        topo = EnzymeRates._catalytic_topologies(
            bi_bi_allosteric_R1_R2)[1]
        dd = EnzymeRates._deduplicate(
            [topo], bi_bi_allosteric_R1_R2)
        result = EnzymeRates._expand_allosteric(
            dd, bi_bi_allosteric_R1_R2; catalytic_n=2,
            allosteric_regs=[:R1, :R2])
        # Separate sites: 2×2 = 4; same site: 2 → total 6
        @test length(result) == 6
    end

    @testset "Dead-end edges: passthrough" begin
        # Spec with dead-end I edges but no allosteric regs
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_dead_end_I)[1]
        de = EnzymeRates._expand_dead_end_inhibitors(
            [topo], uni_uni_dead_end_I;
            dead_end_regs=[:I])
        dd = EnzymeRates._deduplicate(
            de, uni_uni_dead_end_I)
        result = EnzymeRates._expand_allosteric(
            dd, uni_uni_dead_end_I; catalytic_n=1,
            allosteric_regs=Symbol[])
        # Dead-end edges pass through unchanged
        @test length(result) == length(dd)
    end

    @testset "Dead-end I + allosteric R" begin
        topo = EnzymeRates._catalytic_topologies(
            bi_bi_dead_end_I_allosteric_R)[1]
        dd = EnzymeRates._deduplicate(
            [topo], bi_bi_dead_end_I_allosteric_R)
        result = EnzymeRates._expand_allosteric(
            dd, bi_bi_dead_end_I_allosteric_R;
            catalytic_n=1, allosteric_regs=[:R])
        # Only R expands (m=1), I passes through → 1
        @test length(result) == 1
    end

    # Properties
    @testset "Properties" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_bi_allosteric_R_cn2)[1]
        dd = EnzymeRates._deduplicate(
            [topo], uni_bi_allosteric_R_cn2)
        result = EnzymeRates._expand_allosteric(
            dd, uni_bi_allosteric_R_cn2;
            catalytic_n=2, allosteric_regs=[:R])
        for s in result
            @test s.catalytic_n == 2
            @test !isempty(s.allosteric_reg_sites)
        end
    end
end
```

- [ ] **Step 2: Run and commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Add Stage 6 allosteric expansion tests"
```

### Task 9: Stage 7 — TR equivalence

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (add testset)

- [ ] **Step 1: Add Stage 7 testset**

```julia
@testset "Stage 7: TR equivalence" begin
    @testset "Uni-Uni + R: 2^3 = 8 variants" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_allosteric_R)[1]
        dd = EnzymeRates._deduplicate(
            [topo], uni_uni_allosteric_R)
        allo = EnzymeRates._expand_allosteric(
            dd, uni_uni_allosteric_R; catalytic_n=1,
            allosteric_regs=[:R])
        tr = EnzymeRates._expand_tr_equivalence(
            allo, uni_uni_allosteric_R)
        # 3 metabolites with T-state params: S, P, R
        @test length(tr) == 8
    end

    @testset "Uni-Bi + R: 2^4 = 16 variants" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_bi_allosteric_R)[1]
        dd = EnzymeRates._deduplicate(
            [topo], uni_bi_allosteric_R)
        allo = EnzymeRates._expand_allosteric(
            dd, uni_bi_allosteric_R; catalytic_n=1,
            allosteric_regs=[:R])
        tr = EnzymeRates._expand_tr_equivalence(
            allo, uni_bi_allosteric_R)
        # 4 metabolites: S, P, Q, R
        @test length(tr) == 16
    end

    @testset "Bi-Bi + R1, R2: 2^6 = 64 variants" begin
        topo = EnzymeRates._catalytic_topologies(
            bi_bi_allosteric_R1_R2)[1]
        dd = EnzymeRates._deduplicate(
            [topo], bi_bi_allosteric_R1_R2)
        allo = EnzymeRates._expand_allosteric(
            dd, bi_bi_allosteric_R1_R2;
            catalytic_n=1, allosteric_regs=[:R1, :R2])
        # Use first allosteric spec
        tr = EnzymeRates._expand_tr_equivalence(
            [allo[1]], bi_bi_allosteric_R1_R2)
        # 6 metabolites: A, B, P, Q, R1, R2
        @test length(tr) == 64
    end

    @testset "Properties" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_allosteric_R)[1]
        dd = EnzymeRates._deduplicate(
            [topo], uni_uni_allosteric_R)
        allo = EnzymeRates._expand_allosteric(
            dd, uni_uni_allosteric_R;
            catalytic_n=1, allosteric_regs=[:R])
        tr = EnzymeRates._expand_tr_equivalence(
            allo, uni_uni_allosteric_R)

        no_tr = filter(
            s -> isempty(s.tr_equiv_metabolites), tr)
        with_tr = filter(
            s -> !isempty(s.tr_equiv_metabolites), tr)
        @test !isempty(no_tr)
        @test !isempty(with_tr)

        # TR equiv reduces parameter count
        m_no = compile_mechanism(no_tr[1])
        m_with = compile_mechanism(with_tr[end])
        @test length(parameters(m_with)) <
            length(parameters(m_no))
    end
end
```

- [ ] **Step 2: Run and commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Add Stage 7 TR equivalence tests"
```

### Task 10: Stage 8 — Allosteric deduplication

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (add testset)

- [ ] **Step 1: Add Stage 8 testset**

```julia
@testset "Stage 8: Allosteric deduplication" begin
    @testset "T/R mirrors dedup" begin
        # Two specs with complementary tr_equiv_metabolites
        # are T↔R mirrors and should dedup to 1
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_allosteric_R)[1]
        dd = EnzymeRates._deduplicate(
            [topo], uni_uni_allosteric_R)
        allo = EnzymeRates._expand_allosteric(
            dd, uni_uni_allosteric_R;
            catalytic_n=1, allosteric_regs=[:R])
        tr = EnzymeRates._expand_tr_equivalence(
            allo, uni_uni_allosteric_R)
        deduped = EnzymeRates._deduplicate_allosteric(
            tr, uni_uni_allosteric_R)
        # Should be < length(tr) if mirrors exist
        @test length(deduped) <= length(tr)
    end

    @testset "Different base mechanisms survive" begin
        # Two allosteric specs from different base
        # mechanisms should both survive dedup
        topos = EnzymeRates._catalytic_topologies(
            uni_bi_allosteric_R)
        @test length(topos) >= 2
        dd = EnzymeRates._deduplicate(
            topos[1:2], uni_bi_allosteric_R)
        if length(dd) >= 2
            allo = EnzymeRates._expand_allosteric(
                dd, uni_bi_allosteric_R;
                catalytic_n=1, allosteric_regs=[:R])
            tr = EnzymeRates._expand_tr_equivalence(
                allo, uni_bi_allosteric_R)
            deduped = EnzymeRates._deduplicate_allosteric(
                tr, uni_bi_allosteric_R)
            @test length(deduped) >= 2
        end
    end
end
```

- [ ] **Step 2: Run and commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Add Stage 8 allosteric deduplication tests"
```

---

## Chunk 5: End-to-end, param_count, and final validation

### Task 11: End-to-end pipeline and param_count tests

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (add testsets)

- [ ] **Step 1: Add end-to-end testset**

```julia
@testset "End-to-end pipeline" begin
    # These test the full enumerate_mechanisms function.
    # Expected counts will need updating as pipeline
    # bugs are fixed. For now, document current behavior
    # and mark incorrect counts.

    @testset "Uni-Uni, no regs" begin
        result = collect(
            EnzymeRates.enumerate_mechanisms(uni_uni))
        # Current code produces 1 mechanism
        @test length(result) == 1
    end

    @testset "Uni-Bi, no regs" begin
        result = collect(
            EnzymeRates.enumerate_mechanisms(uni_bi))
        # [Implementer: run enumerate_mechanisms and
        #  record the current count. This count may
        #  change when RE/SS and dead-end bugs are fixed.]
        @test length(result) > 0
    end

    @testset "Bi-Bi, no regs" begin
        result = collect(
            EnzymeRates.enumerate_mechanisms(bi_bi))
        @test length(result) > 0
    end
end
```

- [ ] **Step 2: Add param_count accuracy testset**

```julia
@testset "param_count accuracy" begin
    @testset "All Uni-Bi specs" begin
        all_specs = collect(
            EnzymeRates.enumerate_mechanisms(uni_bi))
        @test length(all_specs) > 0
        for s in all_specs
            m = compile_mechanism(s)
            @test s.param_count == length(parameters(m))
        end
    end

    @testset "Sampled Bi-Bi specs" begin
        all_specs = collect(
            EnzymeRates.enumerate_mechanisms(bi_bi))
        base = filter(
            s -> s isa EnzymeRates.MechanismSpec,
            all_specs)
        rng = Random.MersenneTwister(42)
        n = min(20, length(base))
        sample = base[randperm(rng, length(base))[1:n]]
        for s in sample
            m = compile_mechanism(s)
            @test s.param_count == length(parameters(m))
        end
    end
end
```

- [ ] **Step 3: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

This runs ALL tests (not just enumeration). Verify no regressions in other test files.
Expected: some `@test_broken` in enumeration tests, everything else passes.

- [ ] **Step 4: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Add end-to-end pipeline and param_count accuracy tests"
```

### Task 12: Final validation and cleanup

- [ ] **Step 1: Run full test suite one more time**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Verify:
- No test regressions in other files (test_accessors, test_types, test_dsl, etc.)
- All `@test_broken` are genuinely broken (known bugs)
- All `@test` pass
- No stale references to old test spec types (`StageExpansionTestSpec`,
  `EnumerationTestSpec`, `build_stage_expansion_specs`, etc.)

- [ ] **Step 2: Review test file for completeness**

Read through `test/test_mechanism_enumeration.jl` and verify:
- Every stage has a testset
- Round-trip verification exists for all hand-defined mechanisms
- Properties are checked for each stage
- Comments explain why `@test_broken` tests are broken

- [ ] **Step 3: Final commit if any cleanup needed**

```bash
git add test/test_mechanism_enumeration.jl test/mechanism_enumeration_test_specs.jl
git commit -m "Final cleanup of mechanism enumeration tests"
```
