# Mechanism Enumeration Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the beam enumeration and old staged pipeline with a simpler mechanism enumeration module that grows mechanisms incrementally by parameter count.

**Architecture:** Three composable building blocks (`init_mechanisms`, `expand_mechanisms`, `dedup!`) — no monolithic pipeline. The caller owns the loop and cache. Parameter counts are estimated upper bounds during enumeration; true counts come from `parameters()` at compile time.

**Tech Stack:** Julia, EnzymeRates.jl type system (`EnzymeMechanism`, `AllostericEnzymeMechanism`, `EnzymeReaction` singleton types, `@generated` functions)

**Spec:** `docs/superpowers/specs/2026-03-27-mechanism-enumeration-refactor-design.md`

---

## File Map

| File | Action | Responsibility |
|------|--------|---------------|
| `src/mechanism_enumeration.jl` | Rename → `src/old_mechanism_enumeration.jl` | Old staged pipeline (preserved) |
| `src/beam_enumeration.jl` | Rename → `src/old_beam_enumeration.jl` | Old beam pipeline (preserved) |
| `test/test_beam_enumeration.jl` | Rename → `test/old_test_beam_enumeration.jl` | Old beam tests (preserved) |
| `src/mechanism_enumeration.jl` | Create new | All enumeration: types, init, expand, dedup, constructors |
| `test/test_mechanism_enumeration.jl` | Create new | Unit + integration tests for new enumeration |
| `src/types.jl` | Modify lines 7-56, 342-350, 429-441 | Add `oligomeric_state` type param to `EnzymeReaction` |
| `src/dsl.jl` | Modify lines 82-105 | `@enzyme_reaction` accepts `oligomeric_state:` label |
| `src/EnzymeRates.jl` | Modify lines 24, 35-36 | Remove `compile_mechanism` export, update includes |
| `test/runtests.jl` | Modify lines 14-15 | Point to new test files |
| `test/test_types.jl` | Modify lines 2-78 | Add `oligomeric_state` tests |
| `test/test_dsl.jl` | Modify lines 2-41 | Add `oligomeric_state` DSL tests |
| `.claude/CLAUDE.md` | Modify | Update architecture notes |

---

### Task 1: Rename old files and update includes

**Files:**
- Rename: `src/mechanism_enumeration.jl` → `src/old_mechanism_enumeration.jl`
- Rename: `src/beam_enumeration.jl` → `src/old_beam_enumeration.jl`
- Rename: `test/test_beam_enumeration.jl` → `test/old_test_beam_enumeration.jl`
- Modify: `src/EnzymeRates.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Rename source files**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
git mv src/mechanism_enumeration.jl src/old_mechanism_enumeration.jl
git mv src/beam_enumeration.jl src/old_beam_enumeration.jl
git mv test/test_beam_enumeration.jl test/old_test_beam_enumeration.jl
```

- [ ] **Step 2: Update includes in `src/EnzymeRates.jl`**

Change lines 35-36 from:
```julia
include("mechanism_enumeration.jl")
include("beam_enumeration.jl")
```
to:
```julia
include("old_mechanism_enumeration.jl")
include("old_beam_enumeration.jl")
include("mechanism_enumeration.jl")
```

Also remove `compile_mechanism` from the export on line 24.

- [ ] **Step 3: Update test includes in `test/runtests.jl`**

Change lines 14-15 from:
```julia
include("old_test_mechanism_enumeration.jl")
include("test_beam_enumeration.jl")
```
to:
```julia
include("old_test_mechanism_enumeration.jl")
include("old_test_beam_enumeration.jl")
include("test_mechanism_enumeration.jl")
```

- [ ] **Step 4: Create stub files so includes don't fail**

Create empty `src/mechanism_enumeration.jl`:
```julia
# ABOUTME: Mechanism enumeration by incremental parameter count growth
# ABOUTME: Provides init_mechanisms, expand_mechanisms, dedup! building blocks
```

Create empty `test/test_mechanism_enumeration.jl`:
```julia
# ABOUTME: Tests for mechanism enumeration pipeline
# ABOUTME: Unit tests per move + integration tests per reaction
```

- [ ] **Step 5: Run tests to verify nothing broke**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: All existing tests pass (old files just renamed, includes updated).

- [ ] **Step 6: Commit**

```bash
git add src/EnzymeRates.jl src/old_mechanism_enumeration.jl src/old_beam_enumeration.jl src/mechanism_enumeration.jl test/old_test_beam_enumeration.jl test/test_mechanism_enumeration.jl test/runtests.jl
git commit -m "Rename old enumeration files, prepare for new mechanism_enumeration.jl"
```

---

### Task 2: Add `oligomeric_state` type parameter to `EnzymeReaction`

**Files:**
- Modify: `src/types.jl:7-56, 342-350, 429-441`
- Modify: `test/test_types.jl`
- Modify: `test/test_dsl.jl`

- [ ] **Step 1: Write failing tests for `oligomeric_state` in `test/test_types.jl`**

Add to the existing `EnzymeReaction` testset in `test/test_types.jl`:

```julia
@testset "EnzymeReaction oligomeric_state" begin
    # Default oligomeric_state is 1
    rxn = EnzymeReaction(
        ((:S, ((:C, 1),)),),
        ((:P, ((:C, 1),)),),
    )
    @test EnzymeRates.oligomeric_state(rxn) == 1

    # Explicit oligomeric_state
    rxn2 = EnzymeReaction(
        ((:S, ((:C, 1),)),),
        ((:P, ((:C, 1),)),),
        ();
        oligomeric_state=4
    )
    @test EnzymeRates.oligomeric_state(rxn2) == 4

    # Different oligomeric_state = different type
    @test typeof(rxn) !== typeof(rxn2)
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `oligomeric_state` not defined.

- [ ] **Step 3: Implement `oligomeric_state` type parameter**

In `src/types.jl`, change line 14 from:
```julia
struct EnzymeReaction{Substrates, Products, Regulators} end
```
to:
```julia
struct EnzymeReaction{Substrates, Products, Regulators, OligomericState} end
```

Update the constructor (lines 28-55). Change signature from:
```julia
function EnzymeReaction(subs::Tuple, prods::Tuple, regs::Tuple=())
```
to:
```julia
function EnzymeReaction(subs::Tuple, prods::Tuple, regs::Tuple=(); oligomeric_state::Int=1)
```

And change the final construction line from:
```julia
EnzymeReaction{subs, prods, sorted_regs}()
```
to:
```julia
EnzymeReaction{subs, prods, sorted_regs, oligomeric_state}()
```

Update all accessor functions (lines 429-441) to include the 4th type parameter. Change pattern from `::EnzymeReaction{S,P,R}` to `::EnzymeReaction{S,P,R,N}`. Add accessor:
```julia
oligomeric_state(::EnzymeReaction{S,P,R,N}) where {S,P,R,N} = N
```

Update the `show` method (lines 342-350) to include oligomeric_state when it's not 1.

Update `regulator_roles` to use 4-param pattern: `::EnzymeReaction{S,P,R,N}`.

- [ ] **Step 4: Fix all existing pattern matches**

Search for `EnzymeReaction{S,P,R}` in `src/` files and update to `EnzymeReaction{S,P,R,N}` (or `EnzymeReaction{S,P,R}` with a free N if unused). Key locations:

- `src/types.jl`: accessors (lines 429-441)
- `src/rate_eq_derivation.jl`: `@generated` functions using `EnzymeReaction` type params
- `src/fitting.jl`: any EnzymeReaction dispatch
- `src/old_mechanism_enumeration.jl`: functions with `@nospecialize(reaction::EnzymeReaction)` — these use `@nospecialize` so they should work without changes
- `src/old_beam_enumeration.jl`: same — `@nospecialize` should work

Note: `@nospecialize(reaction::EnzymeReaction)` calls don't need updating — they dispatch on the abstract type. Only `where {S,P,R}` destructuring needs the 4th param.

- [ ] **Step 5: Run tests to verify everything passes**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: All tests pass including the new oligomeric_state test.

- [ ] **Step 6: Commit**

```bash
git add src/types.jl src/rate_eq_derivation.jl src/fitting.jl test/test_types.jl
git commit -m "Add oligomeric_state type parameter to EnzymeReaction"
```

---

### Task 3: Add `oligomeric_state:` to `@enzyme_reaction` DSL

**Files:**
- Modify: `src/dsl.jl:82-105`
- Modify: `test/test_dsl.jl`

- [ ] **Step 1: Write failing test in `test/test_dsl.jl`**

Add to the `@enzyme_reaction` testset:

```julia
@testset "@enzyme_reaction with oligomeric_state" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        oligomeric_state: 4
    end
    @test EnzymeRates.oligomeric_state(rxn) == 4

    # Without oligomeric_state defaults to 1
    rxn2 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    @test EnzymeRates.oligomeric_state(rxn2) == 1
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `oligomeric_state:` label not recognized by macro.

- [ ] **Step 3: Update `@enzyme_reaction` macro in `src/dsl.jl`**

In the `macro enzyme_reaction(block)` function (lines 82-105), add parsing for the `oligomeric_state:` label. The macro currently parses `substrates:`, `products:`, `regulators:`, `dead_end_inhibitors:`, `allosteric_regulators:`. Add `oligomeric_state:` as an integer label.

In the label parsing loop, add a branch:
```julia
elseif label == :oligomeric_state
    oligomeric_state = stmt.args[2]
```

Initialize `oligomeric_state = 1` before the loop.

Change the final return from:
```julia
return esc(:(EnzymeReaction($subs, $prods, $regs)))
```
to:
```julia
return esc(:(EnzymeReaction($subs, $prods, $regs; oligomeric_state=$oligomeric_state)))
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: All tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/dsl.jl test/test_dsl.jl
git commit -m "Add oligomeric_state label to @enzyme_reaction DSL"
```

---

### Task 4: New `mechanism_enumeration.jl` — types and helpers

**Files:**
- Modify: `src/mechanism_enumeration.jl` (the new stub)
- Modify: `test/test_mechanism_enumeration.jl` (the new stub)

- [ ] **Step 1: Write failing test for types and `mechanism_spec_from_mechanism` round-trip**

In `test/test_mechanism_enumeration.jl`:

```julia
# ABOUTME: Tests for mechanism enumeration pipeline
# ABOUTME: Unit tests per move + integration tests per reaction

# ── Helper: convert EnzymeMechanism → MechanismSpec ──────────
"""
Convert a compiled EnzymeMechanism back to a MechanismSpec.
Reconstructs StepSpec entries from the mechanism's reactions.
"""
function mechanism_spec_from_mechanism(
    m::EnzymeMechanism,
    @nospecialize(rxn::EnzymeReaction))
    rxns = EnzymeRates.reactions(m)
    eq_steps_tuple = EnzymeRates.equilibrium_steps(m)
    pc = EnzymeRates.param_constraints(m)

    steps = EnzymeRates.StepSpec[]
    for (i, (lhs, rhs)) in enumerate(rxns)
        push!(steps, EnzymeRates.StepSpec(
            collect(lhs), collect(rhs),
            eq_steps_tuple[i]))
    end

    constraints = [
        (t, c, [(s, sc) for (s, sc) in f])
        for (t, c, f) in pc
    ]

    EnzymeRates.MechanismSpec(
        rxn, steps, constraints,
        length(parameters(m)))
end

# ── Reaction definitions ─────────────────────────────────────

const uni_uni_rxn = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
end

const bi_bi_rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end

const bi_bi_pp_rxn = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
end

const uni_bi_rxn = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
end

@testset "Mechanism Enumeration" begin

@testset "Types and round-trip" begin
    # Define a mechanism via @enzyme_mechanism
    m_uu = @enzyme_mechanism begin
        species: begin
            substrates: S[C]
            products: P[C]
            enzymes: E, E_P[C], E_S[C]
        end
        steps: begin
            [E, P] ⇌ [E_P]
            [E, S] ⇌ [E_S]
            [E_S] <--> [E_P]
        end
    end

    # Convert to MechanismSpec
    spec = mechanism_spec_from_mechanism(m_uu, uni_uni_rxn)

    # Verify fields
    @test length(spec.steps) == 3
    @test spec.max_estimated_param_count == length(parameters(m_uu))

    # Round-trip: compile back and verify same mechanism
    m_compiled = EnzymeMechanism(spec)
    @test m_compiled === m_uu
end

end # top-level testset
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `MechanismSpec` not defined in new file, `EnzymeMechanism(spec::MechanismSpec)` not defined.

- [ ] **Step 3: Implement types and constructor in `src/mechanism_enumeration.jl`**

Copy the following from `src/old_mechanism_enumeration.jl` into the new `src/mechanism_enumeration.jl`:

```julia
# ABOUTME: Mechanism enumeration by incremental parameter count growth
# ABOUTME: Provides init_mechanisms, expand_mechanisms, dedup! building blocks

"""Elementary reaction step in canonical binding direction."""
struct StepSpec
    reactants::Vector{Symbol}
    products::Vector{Symbol}
    is_equilibrium::Bool
end

const ParamConstraint = Tuple{
    Symbol, Int, Vector{Tuple{Symbol, Int}}}

abstract type AbstractMechanismSpec end

"""Monomeric enzyme mechanism specification."""
struct MechanismSpec <: AbstractMechanismSpec
    reaction::Any
    steps::Vector{StepSpec}
    param_constraints::Vector{ParamConstraint}
    max_estimated_param_count::Int
end
```

Copy `AllostericMechanismSpec` from `src/old_mechanism_enumeration.jl` (lines 60+). The old struct does NOT have a `param_count` field — it derives param count from its base. Add `max_estimated_param_count::Int` as a new field (the 10th field) so expansion moves can track estimated param count directly on allosteric specs without recomputing.

Copy the `_compile_enzyme_mechanism` function from `src/old_mechanism_enumeration.jl` (the internal function that `compile_mechanism` calls). Add constructor:

```julia
"""Construct EnzymeMechanism from MechanismSpec."""
function EnzymeMechanism(spec::MechanismSpec)
    _compile_enzyme_mechanism(spec)
end
```

Copy the `compile_mechanism(spec::AllostericMechanismSpec)` logic and wrap as:

```julia
"""Construct AllostericEnzymeMechanism from AllostericMechanismSpec."""
function AllostericEnzymeMechanism(spec::AllostericMechanismSpec)
    # ... (copy logic from old compile_mechanism for allosteric)
end
```

Also copy required helpers: `step_metabolite`, `step_forms`, `all_form_names`, `_form_name`, `_dead_end_form_name`, and any other helpers that the constructor depends on.

- [ ] **Step 4: Run tests to verify they pass**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: All tests pass including the round-trip test.

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add types and EnzymeMechanism constructor to new mechanism_enumeration.jl"
```

---

### Task 5: `init_mechanisms` — catalytic topology enumeration

**Files:**
- Modify: `src/mechanism_enumeration.jl`
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing tests for catalytic topology counts and round-trips**

Add to the testset in `test/test_mechanism_enumeration.jl` (adapted from `test/old_test_mechanism_enumeration.jl` lines 224-430):

```julia
@testset "Catalytic topologies" begin

    @testset "Uni-Uni" begin
        topos = EnzymeRates._catalytic_topologies(uni_uni_rxn)
        @test length(topos) == 1

        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end

        @test EnzymeMechanism(topos[1]) === m_uu

        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps
            ) == 1
        end
    end

    @testset "Uni-Bi" begin
        topos = EnzymeRates._catalytic_topologies(uni_bi_rxn)
        @test length(topos) == 3

        # Topo 1: ordered release P-first
        m_ub1 = @enzyme_mechanism begin
            species: begin
                substrates: S[AB]
                products: P[A], Q[B]
                enzymes: E, E_P_Q[AB],
                    E_Q[B], E_S[AB]
            end
            steps: begin
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P_Q]
            end
        end

        # Topo 2: random release
        m_ub2 = @enzyme_mechanism begin
            species: begin
                substrates: S[AB]
                products: P[A], Q[B]
                enzymes: E, E_P[A], E_P_Q[AB],
                    E_Q[B], E_S[AB]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E_P, Q] ⇌ [E_P_Q]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P_Q]
            end
        end

        # Topo 3: ordered release Q-first
        m_ub3 = @enzyme_mechanism begin
            species: begin
                substrates: S[AB]
                products: P[A], Q[B]
                enzymes: E, E_P[A],
                    E_P_Q[AB], E_S[AB]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E_P, Q] ⇌ [E_P_Q]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P_Q]
            end
        end

        defined = [m_ub1, m_ub2, m_ub3]
        for (i, m) in enumerate(defined)
            @test EnzymeMechanism(topos[i]) === m
        end

        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps
            ) == 1
        end
    end

    @testset "Bi-Bi" begin
        topos = EnzymeRates._catalytic_topologies(bi_bi_rxn)
        @test length(topos) == 9
        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps
            ) == 1
        end
    end

    @testset "Bi-Bi Ping-Pong" begin
        topos = EnzymeRates._catalytic_topologies(bi_bi_pp_rxn)
        @test length(topos) == 10
        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps
            ) == 1
        end
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `_catalytic_topologies` not defined in new file.

- [ ] **Step 3: Implement `_catalytic_topologies` in `src/mechanism_enumeration.jl`**

Copy the `_catalytic_topologies` function and all its dependencies from `src/old_mechanism_enumeration.jl`. This includes:
- `_form_name` (line 356)
- `_atoms_dict` (line 370)
- `_add_atoms`, `_subtract_atoms` (atom arithmetic helpers)
- `_can_pingpong` (line ~600)
- `_catalytic_topologies` itself (lines ~430-780)
- Path deduplication by step-set identity

Make sure every topology has exactly 1 SS step (first isomerization) and all other steps RE. Add equivalence constraints for metabolites appearing in multiple RE steps.

- [ ] **Step 4: Run tests to verify they pass**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: All topology count and round-trip tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add _catalytic_topologies to new mechanism_enumeration.jl"
```

---

### Task 6: `init_mechanisms` — dead-end saturation and full init

**Files:**
- Modify: `src/mechanism_enumeration.jl`
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing tests for dead-end saturation and init**

Add to testset in `test/test_mechanism_enumeration.jl`:

```julia
@testset "init_mechanisms" begin

    @testset "Param count invariant" begin
        for (rxn, n_s, n_p) in [
            (uni_uni_rxn, 1, 1),
            (uni_bi_rxn, 1, 2),
            (bi_bi_rxn, 2, 2),
            (bi_bi_pp_rxn, 2, 2),
        ]
            specs = EnzymeRates.init_mechanisms(rxn)
            expected_pc = n_s + n_p + 3
            for s in specs
                @test s.max_estimated_param_count == expected_pc
            end
        end
    end

    @testset "All have exactly 1 SS step" begin
        for rxn in [uni_uni_rxn, uni_bi_rxn,
                     bi_bi_rxn, bi_bi_pp_rxn]
            specs = EnzymeRates.init_mechanisms(rxn)
            for s in specs
                @test count(
                    !st.is_equilibrium for st in s.steps
                ) == 1
            end
        end
    end

    @testset "Uni-Uni: no dead-end forms" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_rxn)
        # 1 topology, no dead-end opportunities
        @test length(specs) == 1
    end

    @testset "Bi-Bi random: dead-end saturation" begin
        # The fully-random bi-bi topology has 4 dead-end
        # opportunities → 2^4 = 16 variants per topology
        # Need to find the random topology among the 9
        m_bb_random = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN],
                    E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_B, A] ⇌ [E_A_B]
                [E, B] ⇌ [E_B]
                [E_A, B] ⇌ [E_A_B]
                [E, P] ⇌ [E_P]
                [E_P, Q] ⇌ [E_P_Q]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        spec_random = mechanism_spec_from_mechanism(
            m_bb_random, bi_bi_rxn)
        # Find init specs matching this topology
        all_specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        # The random topology should produce 16 variants
        # (including the bare topology with 0 dead-ends)
        # The random topology has 9 steps (most of any
        # bi-bi topology). Count specs by step count
        # to identify variants from the random topology.
        random_specs = filter(all_specs) do s
            # Count steps not involving __reg forms
            n_cat = count(s.steps) do st
                !any(contains(string(f), "__reg")
                    for f in vcat(st.reactants, st.products))
            end
            n_cat == 9
        end
        @test length(random_specs) == 16
    end

    @testset "Uni-Bi ordered: no dead-end forms" begin
        # Ordered uni-bi topologies have no mixed
        # substrate+product dead-end opportunities
        specs = EnzymeRates.init_mechanisms(uni_bi_rxn)
        # 2 ordered topologies produce 1 variant each,
        # 1 random topology may produce more
        ordered_count = count(specs) do s
            length(s.steps) == 4  # ordered has 4 steps
        end
        @test ordered_count == 2  # P-first and Q-first
    end

    @testset "Bi-Bi Ping-Pong: dead-end forms" begin
        # Ping-pong with 3 dead-end opportunities → 2^3 = 8
        specs = EnzymeRates.init_mechanisms(bi_bi_pp_rxn)
        @test length(specs) >= 8  # at least one topology has dead-ends
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `init_mechanisms` not defined.

- [ ] **Step 3: Implement `init_mechanisms`**

In `src/mechanism_enumeration.jl`, implement:

```julia
"""
    init_mechanisms(reaction) → Vector{MechanismSpec}

Enumerate all mechanisms at minimum parameter count.
For each catalytic topology: 1 SS step, all K's constrained
equal per metabolite, all substrate/product dead-end subsets
(2^n) with constrained K's.
"""
function init_mechanisms(
    @nospecialize(reaction::EnzymeReaction))
    topos = _catalytic_topologies(reaction)
    result = MechanismSpec[]
    for topo in topos
        _saturate_dead_ends!(result, topo, reaction)
    end
    result
end
```

Copy `_expand_substrate_product_dead_ends` logic from `src/old_mechanism_enumeration.jl` (lines 986-1100) and adapt it into `_saturate_dead_ends!` that:
1. Finds all dead-end opportunities for the topology
2. Enumerates all 2^n subsets
3. For each subset, adds dead-end binding steps with K constrained to catalytic counterpart
4. Adds mirror steps where both endpoints have dead-end extension with same metabolite
5. Computes `max_estimated_param_count = n_substrates + n_products + 3`
6. Pushes each variant to the result vector

Also copy helpers: `_bound_entities`, `_substrate_product_dead_end_opportunities`, `_dead_end_form_name`.

- [ ] **Step 4: Run tests to verify they pass**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: All init_mechanisms tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add init_mechanisms with dead-end saturation"
```

---

### Task 7: Expansion Move 1 — RE→SS conversion

**Files:**
- Modify: `src/mechanism_enumeration.jl`
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing tests for RE→SS move**

```julia
@testset "Move 1: RE→SS conversion" begin
    @testset "Multiple RE steps" begin
        # Uni-uni has 2 RE steps → 2 new specs
        m = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        spec = mechanism_spec_from_mechanism(m, uni_uni_rxn)
        result = EnzymeRates._expand_re_to_ss(spec)
        @test length(result) == 2
        for r in result
            @test r.max_estimated_param_count ==
                spec.max_estimated_param_count + 1
        end
    end

    @testset "Constrained RE steps skipped" begin
        # Bi-bi random with K_A constrained on 2 forms
        # Constrained RE steps should be skipped
        specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        constrained_spec = first(filter(
            s -> !isempty(s.param_constraints), specs))
        # Collect indices of constrained steps
        constrained_idxs = Set{Int}()
        for (_, idx, _) in constrained_spec.param_constraints
            push!(constrained_idxs, idx)
        end
        for (_, _, followers) in constrained_spec.param_constraints
            for (_, idx) in followers
                push!(constrained_idxs, idx)
            end
        end
        n_eligible = count(
            s.is_equilibrium && !(i in constrained_idxs)
            for (i, s) in enumerate(constrained_spec.steps))
        result = EnzymeRates._expand_re_to_ss(
            constrained_spec)
        # Should have exactly one result per eligible RE step
        @test length(result) == n_eligible
        # Verify no result has a newly-SS step that was constrained
        for r in result
            new_ss_idxs = [i for (i, s) in enumerate(r.steps)
                if !s.is_equilibrium &&
                   constrained_spec.steps[i].is_equilibrium]
            for idx in new_ss_idxs
                @test !(idx in constrained_idxs)
            end
        end
    end

    @testset "All SS → yields nothing" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] <--> [E_P]
                [E, S] <--> [E_S]
                [E_S] <--> [E_P]
            end
        end
        spec = mechanism_spec_from_mechanism(m, uni_uni_rxn)
        result = EnzymeRates._expand_re_to_ss(spec)
        @test isempty(result)
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `_expand_re_to_ss` not defined.

- [ ] **Step 3: Implement `_expand_re_to_ss`**

```julia
"""
    _expand_re_to_ss(spec::MechanismSpec) → Vector{MechanismSpec}

Convert one RE step to SS. Skip constrained RE steps.
Each eligible RE step produces one new mechanism at +1.
Dead-end mirror steps that correspond to the converted
catalytic step also inherit SS status.
"""
function _expand_re_to_ss(spec::MechanismSpec)
    result = MechanismSpec[]
    constrained_steps = Set{Int}()
    for (_, idx, _) in spec.param_constraints
        push!(constrained_steps, idx)
    end
    for (_, _, followers) in spec.param_constraints
        for (_, idx) in followers
            push!(constrained_steps, idx)
        end
    end

    for (i, s) in enumerate(spec.steps)
        s.is_equilibrium || continue
        i in constrained_steps && continue

        new_steps = copy(spec.steps)
        new_steps[i] = StepSpec(
            s.reactants, s.products, false)
        # Propagate SS to dead-end mirror steps:
        # find steps connecting dead-end forms that
        # correspond to the endpoints of this catalytic step
        from_form, to_form = step_forms(s)
        for (j, ms) in enumerate(new_steps)
            j == i && continue
            ms.is_equilibrium || continue
            mf, mt = step_forms(ms)
            # Mirror step if both endpoints contain __reg
            # and they correspond to from_form/to_form
            if _is_mirror_of(mf, mt, from_form, to_form)
                new_steps[j] = StepSpec(
                    ms.reactants, ms.products, false)
            end
        end
        push!(result, MechanismSpec(
            spec.reaction, new_steps,
            copy(spec.param_constraints),
            spec.max_estimated_param_count + 1))
    end
    result
end
```

Where `_is_mirror_of(mf, mt, from_form, to_form)` checks whether forms `mf` and `mt` are dead-end extensions of `from_form` and `to_form`. This applies to BOTH regulator dead-ends (forms with `__reg` suffix) AND substrate/product dead-ends (forms like `E_A_P` where an extra metabolite is bound). The check should verify that `mf` and `mt` are non-catalytic forms whose "base" catalytic forms (the forms they extend via dead-end binding) are `from_form` and `to_form` respectively. Use the bound-entities map to determine which catalytic form each dead-end form extends.

- [ ] **Step 4: Run tests to verify they pass**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add Move 1: RE→SS conversion"
```

---

### Task 8: Expansion Move 2 — Remove equivalence constraint

**Files:**
- Modify: `src/mechanism_enumeration.jl`
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing tests**

```julia
@testset "Move 2: Remove equivalence constraint" begin
    @testset "Mechanism with constraints" begin
        # Bi-bi random has K_A and K_B each constrained
        # on 2 forms → 2 removable constraints
        specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        constrained = filter(
            s -> !isempty(s.param_constraints), specs)
        @test !isempty(constrained)
        spec = first(constrained)
        n_constraints = length(spec.param_constraints)
        result = EnzymeRates._expand_remove_constraint(spec)
        @test length(result) == n_constraints
        for r in result
            @test r.max_estimated_param_count ==
                spec.max_estimated_param_count + 1
            @test length(r.param_constraints) ==
                n_constraints - 1
        end
    end

    @testset "No constraints → yields nothing" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_rxn)
        spec = first(specs)
        @test isempty(spec.param_constraints)
        result = EnzymeRates._expand_remove_constraint(spec)
        @test isempty(result)
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 3: Implement `_expand_remove_constraint`**

```julia
"""
    _expand_remove_constraint(spec::MechanismSpec) → Vector{MechanismSpec}

Remove one equivalence constraint (+1 estimated param).
"""
function _expand_remove_constraint(spec::MechanismSpec)
    result = MechanismSpec[]
    for i in eachindex(spec.param_constraints)
        new_constraints = [
            spec.param_constraints[j]
            for j in eachindex(spec.param_constraints)
            if j != i]
        push!(result, MechanismSpec(
            spec.reaction, copy(spec.steps),
            new_constraints,
            spec.max_estimated_param_count + 1))
    end
    result
end
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add Move 2: remove equivalence constraint"
```

---

### Task 9: Expansion Move 3 — Add dead-end regulator

**Files:**
- Modify: `src/mechanism_enumeration.jl`
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing tests**

```julia
const uni_uni_with_reg = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    dead_end_inhibitors: I
end

const bi_bi_with_sub_reg = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    dead_end_inhibitors: A
end

@testset "Move 3: Add dead-end regulator" begin
    @testset "Uni-uni + new regulator" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_with_reg)
        spec = first(specs)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, uni_uni_with_reg)
        # 3 forms (E, E_S, E_P), all eligible
        # 2^3 - 1 = 7 non-empty subsets
        @test length(result) == 7
        for r in result
            @test r.max_estimated_param_count ==
                spec.max_estimated_param_count + 1
        end
    end

    @testset "Mirror step inherits RE/SS" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_with_reg)
        spec = first(specs)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, uni_uni_with_reg)
        # Find variant where I binds to both E_S and E_P
        # (endpoints of the SS step E_S <--> E_P)
        # The mirror step should inherit SS status
        for r in result
            mirror_steps = filter(r.steps) do s
                any(contains(string(x), "__reg")
                    for x in s.reactants) &&
                any(contains(string(x), "__reg")
                    for x in s.products)
            end
            for ms in mirror_steps
                # Mirror of SS step should be SS
                cat_step = first(
                    s for s in spec.steps
                    if !s.is_equilibrium)
                @test ms.is_equilibrium ==
                    cat_step.is_equilibrium
            end
        end
    end

    @testset "No regulators → yields nothing" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_rxn)
        spec = first(specs)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, uni_uni_rxn)
        @test isempty(result)
    end

    @testset "Same metabolite as substrate and regulator" begin
        # A is both substrate and dead-end inhibitor
        specs = EnzymeRates.init_mechanisms(bi_bi_with_sub_reg)
        spec = first(specs)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, bi_bi_with_sub_reg)
        # A-as-regulator should bind to forms where
        # A-as-substrate is already bound (different site)
        @test !isempty(result)
        for r in result
            # Verify __reg suffix in form names
            reg_forms = filter(r.steps) do s
                any(contains(string(x), "__reg")
                    for x in vcat(s.reactants, s.products))
            end
            @test !isempty(reg_forms)
        end
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 3: Implement `_expand_add_dead_end_regulator`**

Copy and adapt from `_expand_dead_end` / `_regulator_dead_end_opportunities` in `src/old_mechanism_enumeration.jl`. Key differences:
- Takes a single spec and reaction (not a vector)
- Returns all non-empty subsets of eligible forms for each new regulator
- All K_R constrained equal within each variant
- Adds mirror steps with inherited RE/SS
- Uses `__reg` suffix for regulator form names
- Only adds regulators NOT already in the mechanism

```julia
"""
    _expand_add_dead_end_regulator(spec, reaction; exclude_regs=Set{Symbol}()) → Vector{MechanismSpec}

Add a new dead-end regulator to all non-empty subsets of
eligible forms. Each variant is +1 (one new K, constrained
equal across forms). `exclude_regs` specifies regulators
to skip (e.g., regulators already present as allosteric
in an AllostericMechanismSpec wrapper).
"""
function _expand_add_dead_end_regulator(
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReaction);
    exclude_regs::Set{Symbol}=Set{Symbol}())
    # ... implementation
    # When checking which regulators are "not yet in the
    # mechanism", also exclude those in exclude_regs
end
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add Move 3: add dead-end regulator"
```

---

### Task 10: Expansion Move 6 — Allosteric conversion (+2)

**Files:**
- Modify: `src/mechanism_enumeration.jl`
- Modify: `test/test_mechanism_enumeration.jl`

Note: Implementing Move 6 before Moves 4-5 because Moves 4-5 operate on `AllostericMechanismSpec` which Move 6 creates.

- [ ] **Step 1: Write failing tests**

```julia
const uni_uni_allo = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    oligomeric_state: 2
end

@testset "Move 6: Allosteric conversion" begin
    @testset "Non-allosteric → allosteric variants" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        # S × {r_only, t_only} + P × {r_only, t_only}
        # = 2 metabolites × 2 flavors = 4 variants
        @test length(result) == 4
        for r in result
            @test r isa EnzymeRates.AllostericMechanismSpec
            @test r.max_estimated_param_count ==
                spec.max_estimated_param_count + 2
        end
    end

    @testset "Already allosteric → yields nothing" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        for a in allo_specs
            result = EnzymeRates._expand_to_allosteric(
                a, uni_uni_allo)
            @test isempty(result)
        end
    end

    @testset "oligomeric_state from reaction" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        for a in allo_specs
            @test a.catalytic_n == 2
        end
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 3: Implement `_expand_to_allosteric`**

```julia
"""
    _expand_to_allosteric(spec::MechanismSpec, reaction) → Vector{AllostericMechanismSpec}
    _expand_to_allosteric(spec::AllostericMechanismSpec, reaction) → Vector{AllostericMechanismSpec}

Convert non-allosteric mechanism to allosteric (+2):
L + one metabolite r_only or t_only. All other metabolites
TR-equivalent. Returns empty for already-allosteric specs.
"""
function _expand_to_allosteric(
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReaction))
    cn = EnzymeRates.oligomeric_state(reaction)
    # Collect eligible metabolites:
    #   - substrates and products (by bare symbol name)
    #   - allosteric/unknown regulators not in mechanism
    #     as dead-end (by __reg suffixed symbol, e.g. :R__reg1)
    # For each eligible metabolite × {r_only, t_only}:
    #   Create AllostericMechanismSpec with:
    #   - all OTHER metabolites in tr_equiv_metabolites
    #     (substrates/products as bare symbols,
    #      regulators as __reg suffixed symbols)
    #   - the one differentiated metabolite in
    #     r_only_metabolites or t_only_metabolites
    #   - all non-binding SS catalytic steps in
    #     tr_equiv_cat_steps
    #   - catalytic_n = cn
    # ...
end

function _expand_to_allosteric(
    ::AllostericMechanismSpec,
    @nospecialize(::EnzymeReaction))
    AllostericMechanismSpec[]
end
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add Move 6: allosteric conversion (+2)"
```

---

### Task 11: Expansion Moves 4-5 — Allosteric regulator and TR equivalence removal

**Files:**
- Modify: `src/mechanism_enumeration.jl`
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing tests for Move 4 (add allosteric regulator)**

```julia
const uni_uni_allo_reg = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: R
    oligomeric_state: 2
end

const uni_uni_allo_2reg = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: R1, R2
    oligomeric_state: 2
end

@testset "Move 4: Add allosteric regulator" begin
    @testset "Add second regulator to allosteric spec" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo_2reg)
        spec = first(specs)
        # First create an allosteric spec with R1
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo_2reg)
        # Find one that used R1 as differentiated metabolite
        # is not right — R1 isn't in init. Let's use a
        # substrate-differentiated one and add R1
        allo = first(allo_specs)
        result = EnzymeRates._expand_add_allosteric_regulator(
            allo, uni_uni_allo_2reg)
        # R1 and R2 available (not yet added):
        # each × {r_only, t_only, tr_equiv} × {new_site}
        # = 2 regs × 3 flavors = 6 variants (min)
        @test length(result) >= 6
        for r in result
            @test r.max_estimated_param_count ==
                allo.max_estimated_param_count + 1
        end
    end

    @testset "Non-allosteric → yields nothing" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo_reg)
        spec = first(specs)
        result = EnzymeRates._expand_add_allosteric_regulator(
            spec, uni_uni_allo_reg)
        @test isempty(result)
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 3: Implement `_expand_add_allosteric_regulator`**

```julia
"""
    _expand_add_allosteric_regulator(spec::AllostericMechanismSpec, reaction) → Vector{AllostericMechanismSpec}

Add a new allosteric regulator: r_only, t_only, or tr_equiv,
on same site as existing or new site. Each variant is +1.
"""
function _expand_add_allosteric_regulator(
    spec::AllostericMechanismSpec,
    @nospecialize(reaction::EnzymeReaction))
    # Find regulators with :unknown or :allosteric role
    # not yet in spec's allosteric_reg_sites
    # For each: 3 flavors × site options
    # ...
end

function _expand_add_allosteric_regulator(
    ::MechanismSpec,
    @nospecialize(::EnzymeReaction))
    AllostericMechanismSpec[]
end
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 5: Write failing tests for Move 5 (remove TR equivalence)**

```julia
@testset "Move 5: Remove TR equivalence" begin
    @testset "Remove metabolite TR equiv" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        allo = first(allo_specs)
        n_tr = length(allo.tr_equiv_metabolites) +
               length(allo.tr_equiv_cat_steps)
        result = EnzymeRates._expand_remove_tr_equiv(
            allo, uni_uni_allo)
        @test length(result) == n_tr
        for r in result
            @test r.max_estimated_param_count ==
                allo.max_estimated_param_count + 1
        end
    end

    @testset "Remove catalytic k TR equiv" begin
        # Ping-pong has multiple SS steps → multiple
        # TR-equivalent k's removable in stages
        specs = EnzymeRates.init_mechanisms(bi_bi_pp_rxn)
        spec = first(specs)
        pp_allo_rxn = @enzyme_reaction begin
            substrates: A[CX], B[N]
            products: P[C], Q[NX]
            oligomeric_state: 2
        end
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, pp_allo_rxn)
        @test !isempty(allo_specs)
        allo = first(allo_specs)
        result = EnzymeRates._expand_remove_tr_equiv(
            allo, pp_allo_rxn)
        # Should include both metabolite and
        # catalytic step TR equiv removals
        @test !isempty(result)
    end

    @testset "No TR equivs left → yields nothing" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        allo = first(allo_specs)
        # Remove all TR equivs
        fully_relaxed = allo
        while true
            result = EnzymeRates._expand_remove_tr_equiv(
                fully_relaxed, uni_uni_allo)
            isempty(result) && break
            fully_relaxed = first(result)
        end
        @test isempty(
            EnzymeRates._expand_remove_tr_equiv(
                fully_relaxed, uni_uni_allo))
    end

    @testset "Same metabolite as substrate and regulator" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: S
            oligomeric_state: 2
        end
        specs = EnzymeRates.init_mechanisms(rxn)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, rxn)
        @test !isempty(allo_specs)
        allo = first(allo_specs)
        result = EnzymeRates._expand_remove_tr_equiv(
            allo, rxn)
        # S-as-substrate appears as :S in tr_equiv_metabolites
        # S-as-regulator appears as :S__reg1 in tr_equiv_metabolites
        # (using __reg suffix convention). Removing each
        # produces a separate result.
        # Plus P-as-substrate TR equiv removal.
        @test length(result) >= 2
        # Verify :S and :S__reg1 are distinct entries
        @test :S in allo.tr_equiv_metabolites
        @test Symbol("S__reg1") in allo.tr_equiv_metabolites
    end
end
```

- [ ] **Step 6: Run test to verify it fails**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 7: Implement `_expand_remove_tr_equiv`**

```julia
"""
    _expand_remove_tr_equiv(spec::AllostericMechanismSpec, reaction) → Vector{AllostericMechanismSpec}

Remove one TR equivalence: metabolite K or catalytic step k.
Each removable equivalence produces one variant at +1.
"""
function _expand_remove_tr_equiv(
    spec::AllostericMechanismSpec,
    @nospecialize(reaction::EnzymeReaction))
    result = AllostericMechanismSpec[]
    # For each tr_equiv metabolite: remove from
    # tr_equiv_metabolites list → +1
    for (i, met) in enumerate(spec.tr_equiv_metabolites)
        new_equiv = [spec.tr_equiv_metabolites[j]
            for j in eachindex(spec.tr_equiv_metabolites)
            if j != i]
        push!(result, AllostericMechanismSpec(
            spec.base, spec.catalytic_n,
            spec.allosteric_reg_sites,
            spec.allosteric_multiplicities,
            new_equiv,
            copy(spec.tr_equiv_cat_steps),
            copy(spec.r_only_metabolites),
            copy(spec.t_only_metabolites),
            copy(spec.r_only_cat_steps),
            spec.max_estimated_param_count + 1))
    end
    # For each tr_equiv catalytic step: remove → +1
    for (i, step_idx) in enumerate(
            spec.tr_equiv_cat_steps)
        new_steps = [spec.tr_equiv_cat_steps[j]
            for j in eachindex(spec.tr_equiv_cat_steps)
            if j != i]
        push!(result, AllostericMechanismSpec(
            spec.base, spec.catalytic_n,
            spec.allosteric_reg_sites,
            spec.allosteric_multiplicities,
            copy(spec.tr_equiv_metabolites),
            new_steps,
            copy(spec.r_only_metabolites),
            copy(spec.t_only_metabolites),
            copy(spec.r_only_cat_steps),
            spec.max_estimated_param_count + 1))
    end
    result
end

function _expand_remove_tr_equiv(
    ::MechanismSpec,
    @nospecialize(::EnzymeReaction))
    AllostericMechanismSpec[]
end
```

- [ ] **Step 8: Run tests to verify they pass**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 9: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add Moves 4-5: allosteric regulator and TR equivalence removal"
```

---

### Task 12: `expand_mechanisms` orchestrator and `dedup!`

**Files:**
- Modify: `src/mechanism_enumeration.jl`
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing tests for dedup**

```julia
@testset "Dedup" begin
    @testset "Same mechanism, different step order" begin
        spec1 = EnzymeRates.MechanismSpec(
            uni_uni_rxn,
            [EnzymeRates.StepSpec([:E, :S], [:E_S], true),
             EnzymeRates.StepSpec([:E, :P], [:E_P], true),
             EnzymeRates.StepSpec([:E_S], [:E_P], false)],
            EnzymeRates.ParamConstraint[],
            5)
        spec2 = EnzymeRates.MechanismSpec(
            uni_uni_rxn,
            [EnzymeRates.StepSpec([:E, :P], [:E_P], true),
             EnzymeRates.StepSpec([:E_S], [:E_P], false),
             EnzymeRates.StepSpec([:E, :S], [:E_S], true)],
            EnzymeRates.ParamConstraint[],
            5)
        cache = Dict(5 => EnzymeRates.AbstractMechanismSpec[
            spec1, spec2])
        EnzymeRates.dedup!(cache)
        @test length(cache[5]) == 1
    end

    @testset "Different mechanisms preserved" begin
        specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        n = length(specs)
        cache = Dict(specs[1].max_estimated_param_count =>
            copy(specs))
        EnzymeRates.dedup!(cache)
        pc = specs[1].max_estimated_param_count
        # Some may dedup, but not all to 1
        @test length(cache[pc]) >= 1
        @test length(cache[pc]) <= n
    end

    @testset "Idempotent" begin
        specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        pc = specs[1].max_estimated_param_count
        cache = Dict(pc => copy(specs))
        EnzymeRates.dedup!(cache)
        n1 = length(cache[pc])
        EnzymeRates.dedup!(cache)
        @test length(cache[pc]) == n1
    end
end
```

- [ ] **Step 2: Write failing tests for `expand_mechanisms`**

```julia
@testset "expand_mechanisms" begin
    @testset "Returns dict keyed by param count" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_rxn)
        spec = first(specs)
        result = EnzymeRates.expand_mechanisms(
            specs, uni_uni_rxn)
        @test result isa Dict{Int, Vector{
            EnzymeRates.AbstractMechanismSpec}}
        base_pc = spec.max_estimated_param_count
        # Should have entries at base_pc + 1 (from moves 1,2)
        @test haskey(result, base_pc + 1)
    end
end
```

- [ ] **Step 3: Run tests to verify they fail**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 4: Implement `dedup!` and `_canonicalize!`**

```julia
"""Canonical sort key for a StepSpec."""
function _step_sort_key(s::StepSpec)
    (sort(s.reactants), sort(s.products), s.is_equilibrium)
end

"""Canonicalize a MechanismSpec for dedup comparison."""
function _canonicalize!(spec::MechanismSpec)
    sort!(spec.steps, by=_step_sort_key)
    sort!(spec.param_constraints,
        by=c -> c[1])
    spec
end

"""Canonicalize an AllostericMechanismSpec."""
function _canonicalize!(spec::AllostericMechanismSpec)
    _canonicalize!(spec.base)
    sort!(spec.tr_equiv_metabolites)
    sort!(spec.tr_equiv_cat_steps)
    sort!(spec.r_only_metabolites)
    sort!(spec.t_only_metabolites)
    sort!(spec.r_only_cat_steps)
    for site in spec.allosteric_reg_sites
        sort!(site)
    end
    spec
end

"""Remove structural duplicates from cache buckets."""
function dedup!(
    cache::Dict{Int, Vector{AbstractMechanismSpec}})
    for (pc, specs) in cache
        for s in specs
            _canonicalize!(s)
        end
        unique!(s -> _dedup_key(s), specs)
        if isempty(specs)
            delete!(cache, pc)
        end
    end
    cache
end
```

`_dedup_key` extracts a hashable representation from the canonicalized spec:

```julia
function _dedup_key(spec::MechanismSpec)
    steps = Tuple(
        (Tuple(s.reactants), Tuple(s.products),
         s.is_equilibrium) for s in spec.steps)
    constraints = Tuple(
        (c[1], c[2], Tuple(c[3]))
        for c in spec.param_constraints)
    (steps, constraints)
end

function _dedup_key(spec::AllostericMechanismSpec)
    base_key = _dedup_key(spec.base)
    (base_key, spec.catalytic_n,
     Tuple(Tuple.(spec.allosteric_reg_sites)),
     Tuple(spec.allosteric_multiplicities),
     Tuple(spec.tr_equiv_metabolites),
     Tuple(spec.tr_equiv_cat_steps),
     Tuple(spec.r_only_metabolites),
     Tuple(spec.t_only_metabolites),
     Tuple(spec.r_only_cat_steps))
end
```

- [ ] **Step 5: Implement `expand_mechanisms`**

```julia
"""
    expand_mechanisms(specs, reaction) → Dict{Int, Vector{AbstractMechanismSpec}}

Apply all +1 and +2 expansion moves. Returns results
grouped by target max_estimated_param_count.
"""
function expand_mechanisms(
    specs::Vector{<:AbstractMechanismSpec},
    @nospecialize(reaction::EnzymeReaction))
    result = Dict{Int, Vector{AbstractMechanismSpec}}()

    for spec in specs
        _add_expansions!(result, spec, reaction)
    end
    result
end

function _add_expansions!(result, spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReaction))
    for new_spec in _expand_re_to_ss(spec)
        _push_to_dict!(result, new_spec)
    end
    for new_spec in _expand_remove_constraint(spec)
        _push_to_dict!(result, new_spec)
    end
    for new_spec in _expand_add_dead_end_regulator(
            spec, reaction)
        _push_to_dict!(result, new_spec)
    end
    for new_spec in _expand_to_allosteric(spec, reaction)
        _push_to_dict!(result, new_spec)
    end
end

function _add_expansions!(result,
    spec::AllostericMechanismSpec,
    @nospecialize(reaction::EnzymeReaction))
    # Move 4: add allosteric regulator
    for new_spec in _expand_add_allosteric_regulator(
            spec, reaction)
        _push_to_dict!(result, new_spec)
    end
    # Move 5: remove TR equivalence
    for new_spec in _expand_remove_tr_equiv(
            spec, reaction)
        _push_to_dict!(result, new_spec)
    end
    # Also apply base mechanism moves to the base,
    # rewrapping each result as AllostericMechanismSpec
    for new_base in _expand_re_to_ss(spec.base)
        _push_to_dict!(result, _rewrap_allosteric(
            spec, new_base))
    end
    for new_base in _expand_remove_constraint(spec.base)
        _push_to_dict!(result, _rewrap_allosteric(
            spec, new_base))
    end
    # For dead-end regulators, exclude regulators already
    # present as allosteric (in spec.allosteric_reg_sites)
    allo_regs = Set{Symbol}()
    for site in spec.allosteric_reg_sites
        for lig in site
            push!(allo_regs, lig)
        end
    end
    for new_base in _expand_add_dead_end_regulator(
            spec.base, reaction;
            exclude_regs=allo_regs)
        _push_to_dict!(result, _rewrap_allosteric(
            spec, new_base))
    end
end

"""
    _rewrap_allosteric(original, new_base) → AllostericMechanismSpec

Create a new AllostericMechanismSpec with updated base
mechanism, preserving all allosteric fields (reg_sites,
multiplicities, TR equiv lists, r/t-only lists). The
max_estimated_param_count is recomputed as:
  original.max_estimated_param_count +
  (new_base.max_estimated_param_count -
   original.base.max_estimated_param_count)
"""
function _rewrap_allosteric(
    original::AllostericMechanismSpec,
    new_base::MechanismSpec)
    delta = new_base.max_estimated_param_count -
        original.base.max_estimated_param_count
    AllostericMechanismSpec(
        new_base, original.catalytic_n,
        deepcopy(original.allosteric_reg_sites),
        copy(original.allosteric_multiplicities),
        copy(original.tr_equiv_metabolites),
        copy(original.tr_equiv_cat_steps),
        copy(original.r_only_metabolites),
        copy(original.t_only_metabolites),
        copy(original.r_only_cat_steps),
        original.max_estimated_param_count + delta)
end

function _push_to_dict!(
    result::Dict{Int, Vector{AbstractMechanismSpec}},
    spec::AbstractMechanismSpec)
    pc = spec.max_estimated_param_count
    push!(get!(result, pc, AbstractMechanismSpec[]), spec)
end
```

- [ ] **Step 6: Run tests to verify they pass**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 7: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add expand_mechanisms orchestrator and dedup!"
```

---

### Task 13: Integration tests

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write integration test helper and tests**

```julia
"""Collect all mechanisms by running the full enumeration loop."""
function enumerate_all(
    @nospecialize(reaction::EnzymeReaction);
    max_params::Int=20)
    cache = Dict{Int, Vector{
        EnzymeRates.AbstractMechanismSpec}}()

    init_specs = EnzymeRates.init_mechanisms(reaction)
    min_pc = first(init_specs).max_estimated_param_count
    cache[min_pc] = init_specs
    EnzymeRates.dedup!(cache)

    results = Dict{Int, Vector{
        EnzymeRates.AbstractMechanismSpec}}()

    for pc in min_pc:max_params
        level = pop!(cache, pc,
            EnzymeRates.AbstractMechanismSpec[])
        isempty(level) && (isempty(cache) ? break :
            continue)

        results[pc] = level

        new_specs = EnzymeRates.expand_mechanisms(
            level, reaction)
        for (target_pc, specs) in new_specs
            target_pc > max_params && continue
            append!(get!(cache, target_pc,
                EnzymeRates.AbstractMechanismSpec[]),
                specs)
        end
        EnzymeRates.dedup!(cache)
    end
    results
end

@testset "Integration" begin
    @testset "Uni-uni full enumeration" begin
        results = enumerate_all(uni_uni_rxn; max_params=10)
        @test !isempty(results)
        pcs = sort(collect(keys(results)))
        # Sorted (gaps possible if no mechanisms at some level)
        @test issorted(pcs)
        # Every mechanism compiles
        for (pc, specs) in results
            for spec in specs
                if spec isa EnzymeRates.MechanismSpec
                    m = EnzymeMechanism(spec)
                    @test length(parameters(m)) <= pc
                elseif spec isa EnzymeRates.AllostericMechanismSpec
                    m = AllostericEnzymeMechanism(spec)
                    @test length(parameters(m)) <= pc
                end
            end
        end
    end

    @testset "Bi-bi full enumeration" begin
        results = enumerate_all(bi_bi_rxn; max_params=12)
        @test !isempty(results)
        for (pc, specs) in results
            for spec in specs
                if spec isa EnzymeRates.MechanismSpec
                    m = EnzymeMechanism(spec)
                    @test length(parameters(m)) <= pc
                elseif spec isa EnzymeRates.AllostericMechanismSpec
                    m = AllostericEnzymeMechanism(spec)
                    @test length(parameters(m)) <= pc
                end
            end
        end
    end

    @testset "With allosteric regulators" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: R
            oligomeric_state: 2
        end
        results = enumerate_all(rxn; max_params=10)
        # Should have both MechanismSpec and
        # AllostericMechanismSpec
        has_allo = any(
            any(s isa EnzymeRates.AllostericMechanismSpec
                for s in specs)
            for (_, specs) in results)
        @test has_allo
    end
end
```

- [ ] **Step 2: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: All integration tests pass.

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Add integration tests for full enumeration loop"
```

---

### Task 14: Remove `compile_mechanism` export and update call sites

**Files:**
- Modify: `src/EnzymeRates.jl`
- Modify: `test/old_test_mechanism_enumeration.jl`
- Modify: `test/old_test_beam_enumeration.jl`
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`

- [ ] **Step 1: Remove `compile_mechanism` from exports in `src/EnzymeRates.jl`**

This was already done in Task 1 Step 2. Verify it's removed.

- [ ] **Step 2: Update old test files to use qualified `EnzymeRates.compile_mechanism`**

In `test/old_test_mechanism_enumeration.jl` and `test/old_test_beam_enumeration.jl`, replace all bare `compile_mechanism(` calls with `EnzymeRates.compile_mechanism(` since it's no longer exported. Use find-and-replace.

Similarly update `test/mechanism_definitions_for_test_enzyme_derivation.jl` if it uses `compile_mechanism`.

- [ ] **Step 3: Run tests to verify nothing broke**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 4: Commit**

```bash
git add src/EnzymeRates.jl test/old_test_mechanism_enumeration.jl test/old_test_beam_enumeration.jl test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -m "Remove compile_mechanism export, qualify old call sites"
```

---

### Task 15: Update CLAUDE.md architecture notes

**Files:**
- Modify: `.claude/CLAUDE.md`

- [ ] **Step 1: Update the mechanism enumeration pipeline section**

Update the "Mechanism enumeration staged pipeline" section in CLAUDE.md to reflect the new architecture:
- Replace 8-stage pipeline description with init/expand/dedup building blocks
- Update `MechanismSpec` field descriptions (`param_count` → `max_estimated_param_count`)
- Add `oligomeric_state` to `EnzymeReaction` description
- Note that `compile_mechanism` is replaced by `EnzymeMechanism(spec)` / `AllostericEnzymeMechanism(spec)` constructors
- Update the "Parameter count: estimated vs. actual" distinction

- [ ] **Step 2: Run tests to verify nothing broke**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 3: Commit**

```bash
git add .claude/CLAUDE.md
git commit -m "Update CLAUDE.md for new mechanism enumeration architecture"
```

---

## Self-Review Checklist

**Spec coverage:**
- ✅ File renames (Task 1)
- ✅ `oligomeric_state` type parameter (Task 2)
- ✅ `@enzyme_reaction` DSL update (Task 3)
- ✅ Types and constructors (Task 4)
- ✅ `init_mechanisms` with catalytic topologies (Task 5)
- ✅ `init_mechanisms` with dead-end saturation (Task 6)
- ✅ Move 1: RE→SS (Task 7)
- ✅ Move 2: Remove constraint (Task 8)
- ✅ Move 3: Add dead-end regulator (Task 9)
- ✅ Move 4: Add allosteric regulator (Task 11)
- ✅ Move 5: Remove TR equivalence (Task 11)
- ✅ Move 6: Allosteric conversion (Task 10)
- ✅ `expand_mechanisms` orchestrator (Task 12)
- ✅ `dedup!` (Task 12)
- ✅ Integration tests (Task 13)
- ✅ `compile_mechanism` removal (Task 14)
- ✅ CLAUDE.md update (Task 15)

**Type consistency:**
- `MechanismSpec` uses `max_estimated_param_count` everywhere ✅
- `EnzymeReaction{S,P,R,N}` pattern used consistently ✅
- `EnzymeMechanism(spec)` constructor used (not `compile_mechanism`) ✅

**TDD followed:**
- Every task writes failing test first, then implementation ✅
