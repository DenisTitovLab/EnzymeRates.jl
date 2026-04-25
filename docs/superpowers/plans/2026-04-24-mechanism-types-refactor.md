# Mechanism Types Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refactor `EnzymeMechanism` and `AllostericEnzymeMechanism` to simpler, non-redundant type parameters. Drop atom tracking at the mechanism level. Replace `constraints:` blocks with kinetic-group encoding via parenthesized step-groups. Add per-kinetic-group TR-mode tagging (`OnlyR`/`OnlyT`/`EqualRT`/`NonequalRT`) via DSL. Split `@enzyme_mechanism` into `@enzyme_mechanism` (plain) and `@allosteric_mechanism` (MWC). Eliminate magic-index type-parameter access. Delete dead code and duplicated R/T-state derivation paths. Substantial net deletion across `src/`.

**Architecture:** `EnzymeMechanism{Metabolites, Reactions}` — `Metabolites` is `((subs_names,), (prods_names,), (regs_names,))` (Symbols only, no atoms); each step is `(lhs, rhs, is_eq, kinetic_group::Int)`; steps with identical `kinetic_group` share kinetic parameters. `AllostericEnzymeMechanism{CatalyticMech, CatSites, RegSites}` — `CatSites = (multiplicity, group_tags)`; `RegSites` entries `(ligands, multiplicity, ligand_tags)`. Accessor-only read interface. No canonicalization in constructor (preserves user step order for predictable parameter naming). Stoichiometric feasibility checked via `rank(S) == rank([S | r])` on the full stoichiometry matrix.

**Tech Stack:** Julia 1.x, `@generated` functions for rate-equation derivation, `LinearAlgebra.rank` (no `Graphs.jl` for validation), `Test`/`Aqua`/`JET` for testing.

**Spec:** `docs/superpowers/specs/2026-04-23-mechanism-types-refactor-design.md`

**Branch:** `allosteric-refactor-spec`.

---

## Execution notes

- The refactor has significant breaking surface. Work in `allosteric-refactor-spec` branch. Main stays untouched until merge.
- Each task either keeps tests green or **explicitly marks RED-OK transitional state**. A "commit" step runs only after the relevant tests pass unless the step header says RED-OK.
- Run tests via `julia --project -e 'using Pkg; Pkg.test()'` (cold). For incremental work use a persistent REPL.
- Code reduction is a first-class success criterion. After major tasks, run `git diff --shortstat origin/main` and check against the targets in spec §9.
- Strict accessor-only rule: grep for `Species[`, `CS[`, `RS[`, `.parameters[` in `src/` after refactor. Zero hits expected.

---

## Phase 0: Preliminaries

### Task 0.1: Capture baseline

**Files:** none (validation only)

- [ ] **Step 1: Run full test suite on baseline**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -40
```

Expected: all tests pass.

- [ ] **Step 2: Record baseline line counts**

```bash
wc -l src/types.jl src/dsl.jl src/rate_eq_derivation.jl src/mechanism_enumeration.jl src/sym_poly_for_rate_eq_derivation.jl src/thermodynamic_constr_for_rate_eq_derivation.jl
```

Save the numbers in a scratch file `REFACTOR_BASELINE.txt` for the line-count check at the end.

---

## Phase 1: Dead-code Removal and Export Cleanup

These are independently committable, keep main green, and reduce surface area for the main refactor.

### Task 1.1: Delete `graph()` accessor

The `graph()` accessor is defined in `src/types.jl:508-521`, used only by `test/test_accessors.jl` (lines 21, 79) for an allocation-smoke test. Not called from any `src/` computation.

**Files:**
- Modify: `src/types.jl`
- Modify: `test/test_accessors.jl`

- [ ] **Step 1: Remove the test references**

In `test/test_accessors.jl`, delete the two lines that call `graph`:
- Line ~21: remove `EnzymeRates.graph(m);`. Other accessor calls on that line stay.
- Line ~79: delete the `@test (@allocated EnzymeRates.graph(m)) == 0` line.

- [ ] **Step 2: Delete the accessor definition**

In `src/types.jl`, delete the `graph` docstring and `@generated function graph(...)` body (currently lines ~504-521).

- [ ] **Step 3: Decide on `Graphs.jl` dependency**

```bash
grep -rn "using Graphs\|import Graphs\|SimpleDiGraph\|add_edge!" src/
```

If `graph()` was the only user (likely), remove `using Graphs` from `src/types.jl:1` and the `Graphs` entry from `Project.toml`. Otherwise leave the import alone.

- [ ] **Step 4: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: pass.

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_accessors.jl Project.toml
git commit -m "Delete unused graph() accessor

Only consumer was an allocation-smoke test; not called from src/
computation paths. Remove the accessor, its test, and the Graphs
dependency if unused elsewhere."
```

### Task 1.2: Delete `RegulatorRole` type hierarchy

Defined in `src/types.jl:18-27` but never dispatched on. All regulator-role logic uses `Symbol`s (`:unknown`, `:dead_end`, `:allosteric`).

**Files:**
- Modify: `src/types.jl`

- [ ] **Step 1: Confirm no dispatch exists**

```bash
grep -rn "::RegulatorRole\|::Allosteric\b\|::DeadEnd\b\|::UnconstrainedRegulator" src/ test/ --include="*.jl"
```

Expected: zero hits beyond the type definitions themselves. If any real dispatch shows up, pause and consult Denis.

- [ ] **Step 2: Delete the declarations**

In `src/types.jl`, delete lines ~17-27 (`abstract type RegulatorRole`, `struct Allosteric`, `struct DeadEnd`, `struct UnconstrainedRegulator`, plus their docstrings).

- [ ] **Step 3: Run tests**

Expected: pass.

- [ ] **Step 4: Commit**

```bash
git add src/types.jl
git commit -m "Delete unused RegulatorRole type hierarchy

Defined but never dispatched on; symbols :unknown / :dead_end /
:allosteric are used directly throughout."
```

### Task 1.3: Remove `compile_mechanism` from exports

`compile_mechanism` is the dispatcher used internally by `identify_rate_equation` and the enumeration pipeline. CLAUDE.md already states it's not user-facing. SPEC.md exports it (stale). Resolution: keep the function as internal, remove it from the export list.

**Files:**
- Modify: `src/EnzymeRates.jl`

- [ ] **Step 1: Verify call sites are all internal**

```bash
grep -rn "compile_mechanism" src/ test/ --include="*.jl"
```

Expected hits: definitions in `src/mechanism_enumeration.jl`, internal calls in `src/identify_rate_equation.jl`, the test file. No user-tutorial consumers.

- [ ] **Step 2: Remove from export list**

In `src/EnzymeRates.jl`, find the line exporting `compile_mechanism` and delete it. Adjacent exports stay.

- [ ] **Step 3: Run tests**

Expected: pass — internal call sites use `EnzymeRates.compile_mechanism` or unqualified within the package, neither requires the export.

- [ ] **Step 4: Commit**

```bash
git add src/EnzymeRates.jl
git commit -m "Remove compile_mechanism from exports

compile_mechanism is an internal dispatcher used by identify_rate_equation
and the enumeration pipeline; it has no user-facing role per CLAUDE.md.
SPEC.md exports table will be updated in the docs cleanup phase."
```

### Task 1.4: Clean stale `old_*.jl` references in CLAUDE.md

These files were cleaned up earlier; CLAUDE.md still references them.

**Files:**
- Modify: `.claude/CLAUDE.md`

- [ ] **Step 1: Confirm files don't exist**

```bash
ls src/old_*.jl test/old_*.jl 2>&1
# Expected: "No such file or directory"
```

- [ ] **Step 2: Delete stale references**

In `.claude/CLAUDE.md`, delete:
- Line ~282: "Old pipeline files preserved as `old_mechanism_enumeration.jl`, `old_beam_enumeration.jl`"
- Line ~293: "`src/old_mechanism_enumeration.jl` — Legacy 8-stage pipeline (preserved, still included for shared helpers)"
- Line ~294: "`src/old_beam_enumeration.jl` — Legacy beam-search pipeline (preserved, still included for shared helpers)"
- Line ~325: "Old enumeration tests preserved in `test/old_test_mechanism_enumeration.jl` and `test/old_test_beam_enumeration.jl`"

- [ ] **Step 3: Commit**

```bash
git add .claude/CLAUDE.md
git commit -m "Remove stale old_*.jl file references from CLAUDE.md

These files were removed in an earlier cleanup; the CLAUDE.md
references were never updated."
```

---

## Phase 2: New `EnzymeMechanism` Type and DSL

This phase replaces the `EnzymeMechanism` type signature and the `@enzyme_mechanism` macro grammar. Breaking change — every plain-mechanism call site migrates in this phase. `AllostericEnzymeMechanism` is touched minimally to keep compilation working; its full refactor is Phase 3.

### Task 2.1: Define new `EnzymeMechanism` struct and accessors

**Files:**
- Modify: `src/types.jl`
- Modify: `test/test_types.jl`

- [ ] **Step 1: Write failing test for new struct + accessors**

Add to `test/test_types.jl`:

```julia
@testset "EnzymeMechanism struct + accessors (new design)" begin
    mets = ((:S,), (:P,), ())
    rxns = (
        ((:E, :S), (:ES,), true,  1),
        ((:ES,),   (:EP,), false, 2),
        ((:EP,),   (:E, :P), true, 3),
    )
    m = EnzymeRates.EnzymeMechanism{mets, rxns}()

    @test EnzymeRates.substrates(m) == (:S,)
    @test EnzymeRates.products(m) == (:P,)
    @test EnzymeRates.regulators(m) == ()
    @test EnzymeRates.metabolites(m) == (:S, :P)
    @test EnzymeRates.reactions(m) == rxns
    @test EnzymeRates.equilibrium_steps(m) == (true, false, true)
    @test EnzymeRates.n_steps(m) == 3
    @test EnzymeRates.kinetic_group(m, 1) == 1
    @test EnzymeRates.kinetic_group(m, 2) == 2
    @test EnzymeRates.kinetic_group(m, 3) == 3
    @test EnzymeRates.kinetic_groups(m) == (1, 2, 3)
    @test EnzymeRates.steps_in_group(m, 1) == (1,)
end
```

- [ ] **Step 2: Run — expect failure**

- [ ] **Step 3: Replace `EnzymeMechanism` struct + accessors**

In `src/types.jl`, replace the existing `EnzymeMechanism` definition with:

```julia
"""
    EnzymeMechanism{Metabolites, Reactions}

Singleton type encoding an enzyme mechanism.

- `Metabolites`: 3-tuple `(substrates::Tuple{Symbol,...}, products::Tuple{Symbol,...},
  regulators::Tuple{Symbol,...})`. Plain symbol names — no atom content stored.
- `Reactions`: tuple of 4-tuples `(lhs_syms, rhs_syms, is_eq::Bool,
  kinetic_group::Int)`. Steps with identical `kinetic_group` share kinetic
  parameters (one `K` for RE groups, one `k_f` and one `k_r` for SS groups).
"""
struct EnzymeMechanism{Metabolites, Reactions} <: AbstractEnzymeMechanism end
```

Replace the existing accessors. Delete the now-invalid `Species[k]` indexing throughout `src/types.jl`. New accessor implementations:

```julia
substrates(::EnzymeMechanism{M}) where {M} = M[1]
products(::EnzymeMechanism{M}) where {M} = M[2]
regulators(::EnzymeMechanism{M}) where {M} = M[3]

@generated function metabolites(::EnzymeMechanism{M}) where {M}
    seen = Set{Symbol}()
    names = Symbol[]
    for group in M
        for name in group
            if name ∉ seen
                push!(seen, name)
                push!(names, name)
            end
        end
    end
    Tuple(names)
end

reactions(::EnzymeMechanism{M, R}) where {M, R} = R

@generated function equilibrium_steps(::EnzymeMechanism{M, R}) where {M, R}
    Tuple(step[3] for step in R)
end

n_steps(::EnzymeMechanism{M, R}) where {M, R} = length(R)

kinetic_group(::EnzymeMechanism{M, R}, idx::Int) where {M, R} = R[idx][4]

@generated function kinetic_groups(::EnzymeMechanism{M, R}) where {M, R}
    Tuple(sort(unique(step[4] for step in R)))
end

@generated function steps_in_group(
    ::EnzymeMechanism{M, R}, ::Val{G},
) where {M, R, G}
    Tuple(i for (i, step) in enumerate(R) if step[4] == G)
end
steps_in_group(m::EnzymeMechanism, g::Int) = steps_in_group(m, Val(g))

@generated function enzyme_forms(::EnzymeMechanism{M, R}) where {M, R}
    met_names = Set{Symbol}()
    for group in M; for name in group; push!(met_names, name); end; end
    seen = Set{Symbol}()
    forms = Symbol[]
    for (lhs, rhs, _, _) in R
        for s in lhs; s ∉ met_names && s ∉ seen && (push!(seen, s); push!(forms, s)); end
        for s in rhs; s ∉ met_names && s ∉ seen && (push!(seen, s); push!(forms, s)); end
    end
    Tuple(forms)
end

n_states(m::EnzymeMechanism) = length(enzyme_forms(m))
```

- [ ] **Step 4: Run tests**

The new accessor test passes; many other tests break (they use the old DSL or old type). **RED-OK from this point through end of Phase 2.**

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "Introduce new EnzymeMechanism{Metabolites, Reactions} type

Metabolites is ((subs,), (prods,), (regs,)) of plain Symbols — no
atoms. Each step tuple is (lhs, rhs, is_eq, kinetic_group). Steps
sharing kinetic_group share K (RE) or k_f/k_r (SS).

KNOWN RED: downstream code uses old 4-param signature."
```

### Task 2.2: Refactor `stoich_matrix` to return full matrix

The existing `stoich_matrix` returns metabolites-only rows. We extend it to return the full matrix (enzymes followed by metabolites) for the rank-based stoichiometric feasibility check, and migrate the one existing caller.

**Files:**
- Modify: `src/types.jl`
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl`

- [ ] **Step 1: Replace `stoich_matrix` definition**

In `src/types.jl`, replace the existing `@generated function stoich_matrix(...)` with:

```julia
"""
    stoich_matrix(m::EnzymeMechanism) → Matrix{Int}

Full stoichiometry matrix. Rows are species in the order
`(enzyme_forms..., metabolites...)` (use `enzyme_row_range(m)` and
`metabolite_row_range(m)` to slice). Columns are step indices.
Positive = produced; negative = consumed in the forward direction.

Enzyme-row columns sum to zero by construction (each step has one
enzyme on each side).
"""
@generated function stoich_matrix(::EnzymeMechanism{M, R}) where {M, R}
    met_names_set = Set{Symbol}()
    for group in M; for name in group; push!(met_names_set, name); end; end

    seen = Set{Symbol}()
    enz = Symbol[]
    for (lhs, rhs, _, _) in R
        for s in lhs; s ∉ met_names_set && s ∉ seen && (push!(seen, s); push!(enz, s)); end
        for s in rhs; s ∉ met_names_set && s ∉ seen && (push!(seen, s); push!(enz, s)); end
    end

    met_seen = Set{Symbol}()
    mets = Symbol[]
    for group in M
        for name in group
            name ∉ met_seen && (push!(met_seen, name); push!(mets, name))
        end
    end

    species = [enz; mets]
    sp_idx = Dict(s => i for (i, s) in enumerate(species))
    S = zeros(Int, length(species), length(R))
    for (j, (lhs, rhs, _, _)) in enumerate(R)
        for s in lhs; S[sp_idx[s], j] -= 1; end
        for s in rhs; S[sp_idx[s], j] += 1; end
    end
    S
end

enzyme_row_range(m::EnzymeMechanism) = 1:n_states(m)
metabolite_row_range(m::EnzymeMechanism) = (n_states(m) + 1):(n_states(m) + length(metabolites(m)))
```

- [ ] **Step 2: Migrate the one existing caller**

In `src/thermodynamic_constr_for_rate_eq_derivation.jl:116`, change:

```julia
# Before:
stoich_matrix(m), collect(metabolites(m)),
# After:
stoich_matrix(m)[metabolite_row_range(m), :], collect(metabolites(m)),
```

- [ ] **Step 3: Verify the slice returns the same shape and values**

Add a quick smoke test in `test/test_types.jl`:

```julia
@testset "stoich_matrix has expected enzyme/metabolite rows" begin
    mets = ((:S,), (:P,), ())
    rxns = (
        ((:E, :S), (:ES,), true,  1),
        ((:ES,),   (:EP,), false, 2),
        ((:EP,),   (:E, :P), true, 3),
    )
    m = EnzymeRates.EnzymeMechanism{mets, rxns}()
    S = EnzymeRates.stoich_matrix(m)

    enz_idx = EnzymeRates.enzyme_row_range(m)
    met_idx = EnzymeRates.metabolite_row_range(m)
    @test all(sum(S[enz_idx, j]) == 0 for j in 1:size(S, 2))   # enzyme conservation
    @test S[met_idx, :] |> size == (2, 3)                       # 2 metabolites × 3 steps
end
```

- [ ] **Step 4: RED-OK; commit**

```bash
git add src/types.jl src/thermodynamic_constr_for_rate_eq_derivation.jl test/test_types.jl
git commit -m "Refactor stoich_matrix to return full enzyme + metabolite matrix

Adds enzyme_row_range / metabolite_row_range slicing accessors.
Migrates the one existing caller (thermodynamic_constr) to slice
metabolite rows. Used by the new rank-based stoichiometric
feasibility check in the EnzymeMechanism constructor."
```

### Task 2.3: Implement `EnzymeMechanism(metabolites, reactions)` constructor

**Files:**
- Modify: `src/types.jl`
- Modify: `test/test_types.jl`

- [ ] **Step 1: Write failing constructor test**

Add to `test/test_types.jl`:

```julia
@testset "EnzymeMechanism constructor" begin
    mets = ((:S,), (:P,), ())
    rxns = (
        ((:E, :S), (:ES,), true,  1),
        ((:ES,),   (:EP,), false, 2),
        ((:EP,),   (:E, :P), true, 3),
    )
    m = EnzymeRates.EnzymeMechanism(mets, rxns)
    @test EnzymeRates.n_steps(m) == 3
    @test EnzymeRates.substrates(m) == (:S,)

    # Same-kinetics group test
    rxns_grouped = (
        ((:E, :S),   (:ES,),   true,  1),
        ((:ES, :S),  (:ESS,),  true,  1),  # same group as step 1
        ((:ESS,),    (:EP,),   false, 2),
        ((:EP,),     (:E, :P), true,  3),
    )
    m_g = EnzymeRates.EnzymeMechanism(mets, rxns_grouped)
    @test EnzymeRates.kinetic_group(m_g, 1) == EnzymeRates.kinetic_group(m_g, 2)

    # Stoichiometry violation: substrate not actually consumed
    bad_rxns = (
        ((:E, :S), (:ES,), true,  1),
        ((:ES,),   (:E,),  false, 2),    # S "vanishes" — no product
    )
    @test_throws ErrorException EnzymeRates.EnzymeMechanism(((:S,), (:P,), ()), bad_rxns)

    # Iso group with size > 1 should error
    bad_iso = (
        ((:E, :S), (:ES,), true,  1),
        ((:ES,),   (:EP,), false, 99),
        ((:EP,),   (:EQ,), false, 99),    # second iso step in same group → error
        ((:EQ,),   (:E, :P), true, 2),
    )
    @test_throws ErrorException EnzymeRates.EnzymeMechanism(((:S,), (:P,), ()), bad_iso)
end
```

- [ ] **Step 2: Run — expect failure**

- [ ] **Step 3: Implement the constructor**

In `src/types.jl`, replace the existing `function EnzymeMechanism(species::Tuple, reactions::Tuple, eq_steps, constraints=())` with the new 2-arg constructor. Validation steps in order:

```julia
"""
    EnzymeMechanism(metabolites, reactions) → EnzymeMechanism

Construct an `EnzymeMechanism` from explicit metabolite 3-tuple and
reaction 4-tuples. Step order is preserved (no canonicalization).
Validates structure, enzyme-form connectivity, kinetic-group
composition rules, and stoichiometric feasibility.
"""
function EnzymeMechanism(
    mets::Tuple{Tuple{Vararg{Symbol}}, Tuple{Vararg{Symbol}}, Tuple{Vararg{Symbol}}},
    rxns::Tuple,
)
    subs, prods, regs = mets

    # ---- Sort each species list alphabetically (canonical form for type
    #      uniqueness). Step ORDER is preserved per spec §4.1.2.
    subs = Tuple(sort(collect(subs)))
    prods = Tuple(sort(collect(prods)))
    regs = Tuple(sort(collect(regs)))
    mets = (subs, prods, regs)

    isempty(rxns) && error("Reactions tuple must not be empty")
    all(step[3] for step in rxns) &&
        error("At least one SS step required (not all steps can be RE)")

    # ---- Build metabolite set
    met_set = Set{Symbol}()
    for group in (subs, prods, regs)
        for name in group; push!(met_set, name); end
    end

    # ---- Validate each step shape
    for (i, step) in enumerate(rxns)
        length(step) == 4 ||
            error("Step $i must be (lhs, rhs, is_eq, kinetic_group); got $step")
        lhs, rhs, is_eq, gnum = step
        is_eq isa Bool || error("Step $i is_eq must be Bool")
        gnum isa Int || error("Step $i kinetic_group must be Int")
        # Each side has exactly one enzyme form (= one non-metabolite name)
        n_enz_lhs = count(s -> s ∉ met_set, lhs)
        n_enz_rhs = count(s -> s ∉ met_set, rhs)
        n_enz_lhs == 1 ||
            error("Step $i LHS must contain exactly one enzyme form; got $lhs")
        n_enz_rhs == 1 ||
            error("Step $i RHS must contain exactly one enzyme form; got $rhs")
        n_met_lhs = count(s -> s ∈ met_set, lhs)
        n_met_rhs = count(s -> s ∈ met_set, rhs)
        n_met_lhs <= 1 || error("Step $i LHS has more than one metabolite")
        n_met_rhs <= 1 || error("Step $i RHS has more than one metabolite")
    end

    # ---- Canonicalize RE step direction (metabolite on LHS for binding steps)
    rxns = ntuple(length(rxns)) do i
        (lhs, rhs, is_eq, gnum) = rxns[i]
        if !is_eq
            return (lhs, rhs, is_eq, gnum)
        end
        rhs_has_met = any(s in met_set for s in rhs)
        lhs_has_met = any(s in met_set for s in lhs)
        if rhs_has_met && !lhs_has_met
            (rhs, lhs, is_eq, gnum)
        else
            (lhs, rhs, is_eq, gnum)
        end
    end

    # ---- Each substrate / product / regulator must appear in some step
    appears = Set{Symbol}()
    for (lhs, rhs, _, _) in rxns
        for s in lhs; appears |= Set([s]); end
        for s in rhs; appears |= Set([s]); end
    end
    for name in vcat(collect(subs), collect(prods), collect(regs))
        name in appears ||
            error("Listed metabolite $name does not appear in any reaction step")
    end

    # ---- No unlisted metabolite-looking name (anything in steps not in met_set
    #      is treated as enzyme form; no extra check needed).

    # ---- Kinetic-group composition rules
    _validate_kinetic_groups(rxns, met_set)

    # ---- Build the singleton type and run remaining checks via accessors
    m = EnzymeMechanism{mets, rxns}()

    # ---- Enzyme-form graph weakly connected
    _validate_enzyme_connectivity(m)

    # ---- Stoichiometric feasibility (rank check)
    _validate_stoichiometry(m)

    m
end
```

Add helpers:

```julia
"""
Validate kinetic-group composition: 2+ groups must be all RE binding
or all SS binding (same metabolite); iso steps must be singletons.
"""
function _validate_kinetic_groups(rxns, met_set)
    groups = Dict{Int, Vector{Int}}()
    for (i, step) in enumerate(rxns)
        push!(get!(groups, step[4], Int[]), i)
    end
    for (g, idxs) in groups
        length(idxs) == 1 && continue
        kinds = map(idxs) do i
            lhs, rhs, is_eq, _ = rxns[i]
            mets_in = [s for s in lhs if s in met_set]
            mets_out = [s for s in rhs if s in met_set]
            isempty(mets_in) && isempty(mets_out) &&
                error("Iso step (no metabolite) at index $i must be a " *
                      "singleton kinetic group; found in group $g of size $(length(idxs))")
            length(mets_in) == 1 ||
                error("Step $i has $(length(mets_in)) metabolites on LHS; expected 1")
            (is_eq, mets_in[1])
        end
        first_kind = kinds[1]
        for (i, k) in zip(idxs[2:end], kinds[2:end])
            k[1] == first_kind[1] ||
                error("Kinetic group $g contains both RE and SS binding steps")
            k[2] == first_kind[2] ||
                error("Kinetic group $g binds different metabolites: " *
                      "$(first_kind[2]) and $(k[2])")
        end
    end
end

"""Verify the enzyme-form graph is weakly connected."""
function _validate_enzyme_connectivity(m::EnzymeMechanism)
    enz = enzyme_forms(m)
    isempty(enz) && error("Mechanism has no enzyme forms")
    name_set = Set(enz)
    adj = Dict(n => Set{Symbol}() for n in enz)
    for (lhs, rhs, _, _) in reactions(m)
        e_l = first(s for s in lhs if s in name_set)
        e_r = first(s for s in rhs if s in name_set)
        push!(adj[e_l], e_r)
        push!(adj[e_r], e_l)
    end
    visited = Set{Symbol}()
    queue = [first(enz)]
    while !isempty(queue)
        cur = popfirst!(queue)
        cur in visited && continue
        push!(visited, cur)
        for n in adj[cur]; n in visited || push!(queue, n); end
    end
    visited == name_set ||
        error("Enzyme-form graph not connected; orphan forms: " *
              "$(setdiff(name_set, visited))")
end

"""
Stoichiometric feasibility via `r ∈ col(S)` rank test on the full
stoichiometry matrix. The target vector r has 0 on enzyme rows and
on regulator rows, and ±count(M in subs/prods) on substrate/product rows.
"""
function _validate_stoichiometry(m::EnzymeMechanism)
    S = stoich_matrix(m)
    species = (enzyme_forms(m)..., metabolites(m)...)
    sp_idx = Dict(s => i for (i, s) in enumerate(species))
    r = zeros(Int, length(species))
    for s in substrates(m); r[sp_idx[s]] -= 1; end
    for p in products(m);   r[sp_idx[p]] += 1; end

    rs = Rational.(S)
    rr = Rational.(r)
    rank(rs) == rank(hcat(rs, rr)) ||
        error("Mechanism stoichiometry does not match the declared net reaction. " *
              "Check substrate / product multiplicities and that regulators " *
              "have net zero change. Declared: " *
              "$(_pretty_reaction(substrates(m), products(m)))")
end

_pretty_reaction(subs, prods) =
    "$(join(string.(subs), " + ")) → $(join(string.(prods), " + "))"
```

- [ ] **Step 4: Run tests**

The new constructor tests pass. Other plain-mechanism tests still failing (still RED).

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "Implement EnzymeMechanism 2-arg constructor

Validation: at least one SS step, each species appears in steps,
exactly one enzyme form per side, kinetic-group composition rules,
weakly-connected enzyme-form graph, rank-based stoichiometric
feasibility (r ∈ col(S)). No canonicalization of step order; user
DSL order preserved.

KNOWN RED: downstream code uses old DSL / old type."
```

### Task 2.4: Rewrite `@enzyme_mechanism` macro

**Files:**
- Modify: `src/dsl.jl`
- Modify: `test/test_dsl.jl`

- [ ] **Step 1: Write failing DSL test**

```julia
@testset "@enzyme_mechanism (new grammar)" begin
    m = @enzyme_mechanism begin
        substrates: S
        products:   P
        regulators: I

        steps: begin
            ([E, S] ⇌ [ES], [EP, S] ⇌ [EPS])
            [ES, I] ⇌ [ESI]
            [ES]   <--> [EP]
            [EP]   ⇌   [E, P]
        end
    end
    @test EnzymeRates.substrates(m) == (:S,)
    @test EnzymeRates.products(m) == (:P,)
    @test EnzymeRates.regulators(m) == (:I,)
    @test EnzymeRates.kinetic_group(m, 1) == EnzymeRates.kinetic_group(m, 2)
    @test EnzymeRates.kinetic_group(m, 3) != EnzymeRates.kinetic_group(m, 4)

    # Reject atom bracket syntax
    @test_throws Exception eval(:(@enzyme_mechanism begin
        substrates: S[C]
        products:   P
        steps: begin
            [E, S] ⇌ [ES]
            [ES] <--> [EP]
            [EP] ⇌ [E, P]
        end
    end))

    # Reject allosteric-only syntax
    @test_throws Exception eval(:(@enzyme_mechanism begin
        substrates: S
        products:   P
        site(:catalytic, 2): begin
            steps: begin
                [E, S] ⇌ [ES]
                [ES] <--> [EP]
                [EP] ⇌ [E, P]
            end
        end
    end))
end
```

- [ ] **Step 2: Run — expect failure**

- [ ] **Step 3: Replace `@enzyme_mechanism` body**

Replace the existing `@enzyme_mechanism` body in `src/dsl.jl` with new grammar. The full replacement:

```julia
"""
    @enzyme_mechanism begin
        substrates: S
        products:   P
        regulators: I

        steps: begin
            ([E, S] ⇌ [ES], [EP, S] ⇌ [EPS])    # parenthesized → shared kinetics
            [ES, I] ⇌ [ESI]
            [ES]   <--> [EP]
            [EP]   ⇌   [E, P]
        end
    end

Build a plain (non-allosteric) `EnzymeMechanism`.
- `substrates:`, `products:`, `regulators:` accept comma-separated bare symbols.
  Atom brackets (e.g. `S[C]`) are rejected — atom declarations belong in `@enzyme_reaction`.
- `enzymes:` block deleted (forms inferred from steps).
- `constraints:` block deleted (same-kinetics groups are expressed by parenthesizing
  the steps that share parameters).
- Allosteric-only constructs (`site(...)` blocks, `::Tag` annotations,
  `allosteric_regulators:`, `catalytic_inhibitors:`) are rejected with clear errors.
"""
macro enzyme_mechanism(block)
    _reject_allosteric_syntax!(block)
    mets_expr, rxns_expr = _parse_plain_mechanism_body(block)
    return esc(:(EnzymeMechanism($mets_expr, $rxns_expr)))
end

function _reject_allosteric_syntax!(block)
    for arg in block.args
        arg isa LineNumberNode && continue
        if arg isa Expr && arg.head == :call && arg.args[1] == :(:)
            label = arg.args[2]
            if label isa Expr && label.head == :call && label.args[1] == :site
                error("@enzyme_mechanism: `site(...)` belongs in @allosteric_mechanism")
            end
            label in (:allosteric_regulators, :catalytic_inhibitors) &&
                error("@enzyme_mechanism: `$label:` is allosteric-only; " *
                      "use @allosteric_mechanism instead")
        end
    end
end

function _parse_plain_mechanism_body(block)
    subs_list, prods_list, regs_list = Symbol[], Symbol[], Symbol[]
    steps_block = nothing
    for arg in block.args
        arg isa LineNumberNode && continue
        arg isa Expr && arg.head == :call && arg.args[1] == :(:) ||
            error("Unexpected expression: $arg")
        label, value = arg.args[2], arg.args[3]
        if label == :substrates
            append!(subs_list, _parse_bare_symbol_list(value, label))
        elseif label == :products
            append!(prods_list, _parse_bare_symbol_list(value, label))
        elseif label == :regulators
            append!(regs_list, _parse_bare_symbol_list(value, label))
        elseif label == :steps
            steps_block = value
        else
            error("Unknown @enzyme_mechanism label: $label")
        end
    end
    isempty(subs_list) && error("substrates: not specified")
    isempty(prods_list) && error("products: not specified")
    steps_block === nothing && error("steps: not specified")

    rxns_expr = _parse_steps_block_with_groups(steps_block)
    mets_expr = Expr(:tuple,
        Expr(:tuple, QuoteNode.(subs_list)...),
        Expr(:tuple, QuoteNode.(prods_list)...),
        Expr(:tuple, QuoteNode.(regs_list)...),
    )
    mets_expr, rxns_expr
end

"""
Parse `substrates: S, A` — comma-separated bare Symbols. Reject
atom brackets and tag annotations.
"""
function _parse_bare_symbol_list(value, label)
    syms = Symbol[]
    function push_one(arg)
        if arg isa Symbol
            push!(syms, arg)
        elseif arg isa Expr && arg.head == :ref
            error("@enzyme_mechanism: atom bracket syntax `$arg` is not allowed " *
                  "at the mechanism level; declare atoms in @enzyme_reaction.")
        elseif arg isa Expr && arg.head == :(::)
            error("@enzyme_mechanism: tag annotation `$arg` is not allowed; " *
                  "tags are only valid in @allosteric_mechanism.")
        else
            error("@enzyme_mechanism `$label:` expects bare Symbol names; got $arg")
        end
    end
    if value isa Expr && value.head == :tuple
        for a in value.args; push_one(a); end
    else
        push_one(value)
    end
    syms
end

"""
Parse the steps block. Each line is either a single step or a
parenthesized tuple of steps sharing kinetics. Returns a tuple-Expr of
4-tuples `(lhs, rhs, is_eq, kinetic_group)`.
"""
function _parse_steps_block_with_groups(steps_block)
    next_group = Ref(0)
    rxns = Expr(:tuple)
    for arg in steps_block.args
        arg isa LineNumberNode && continue
        if arg isa Expr && arg.head == :tuple
            # Parenthesized group: all share one kinetic group
            next_group[] += 1
            gnum = next_group[]
            for e in arg.args
                push!(rxns.args, _parse_single_step(e, gnum, allow_tag=false))
            end
        else
            next_group[] += 1
            gnum = next_group[]
            push!(rxns.args, _parse_single_step(arg, gnum, allow_tag=false))
        end
    end
    rxns
end

"""
Parse a single step `[lhs] ⇌ [rhs]` or `[lhs] <--> [rhs]`. With
`allow_tag=false` (plain mechanism), reject any `::Tag` postfix. With
`allow_tag=true` (allosteric), the caller handles the tag.
"""
function _parse_single_step(expr, gnum::Int; allow_tag::Bool=false)
    if allow_tag && expr isa Expr && expr.head == :(::)
        # Caller is responsible for unwrapping the tag in allosteric context.
        # In plain context this branch isn't reached (allow_tag=false).
    end
    if !allow_tag && expr isa Expr && expr.head == :(::)
        error("@enzyme_mechanism: tag annotation `$expr` is not allowed; " *
              "tags are only valid in @allosteric_mechanism.")
    end
    expr isa Expr && expr.head == :call ||
        error("Expected [lhs] ⇌ [rhs] or [lhs] <--> [rhs]; got $expr")
    op = expr.args[1]
    is_eq = op == :⇌
    is_eq || op == :(<-->) ||
        error("Expected ⇌ or <--> step operator; got $op")
    lhs = _parse_step_side_symbols(expr.args[2])
    rhs = _parse_step_side_symbols(expr.args[3])
    Expr(:tuple, lhs, rhs, is_eq, gnum)
end
```

Reuse the existing `_parse_step_side_symbols` (which expects `[a, b, c]` vector syntax and produces a tuple of QuoteNodes).

Delete the old `_parse_enzyme_mechanism`, `_parse_allosteric_mechanism`, `_is_allosteric_label`, `_parse_constraint_rhs`, `_walk_rhs!`, `_push_constraint!`, `_parse_constraints_block`, `_parse_species_block`, `_parse_steps_block`, `_parse_chemical_formula`, `_parse_species_tuple_expr`, `_parse_label_species_tuple`, `_parse_labeled_block`, `_regulator_tuple_to_symbols`, `_parse_reg_ligands_block`, `_parse_catalytic_block`, `_met_sym` functions to the extent they're only used by the deleted DSL.

Note: `@enzyme_reaction` macro at the top of `dsl.jl` uses `_parse_chemical_formula`, `_parse_species_tuple_expr`, `_parse_label_species_tuple`, `_parse_labeled_block`, `_regulator_tuple_to_symbols`. These must be preserved. Audit before deletion.

- [ ] **Step 4: Run tests**

The new DSL test passes; existing tests using old grammar still fail (Phase 2.6 migrates them).

- [ ] **Step 5: Commit**

```bash
git add src/dsl.jl test/test_dsl.jl
git commit -m "Rewrite @enzyme_mechanism with new flat grammar

substrates: / products: / regulators: at top level (bare Symbols only;
atom brackets rejected). steps: block accepts standalone step lines and
parenthesized step-groups (shared kinetics). Allosteric-only syntax
rejected with clear errors. Constraint-DSL parser (_walk_rhs!,
_parse_constraint_rhs, _push_constraint!) deleted.

KNOWN RED: existing test mechanisms use old grammar."
```

### Task 2.5: Add `@allosteric_mechanism` macro skeleton

This task adds the macro and its DSL parser, sufficient to construct `AllostericEnzymeMechanism` types using the OLD allosteric type signature. The new allosteric signature comes in Phase 3, after which the macro emits the new shape. Splitting this way lets us migrate all plain-mechanism tests (Task 2.6) before touching the allosteric type itself.

**Files:**
- Modify: `src/dsl.jl`
- Modify: `src/EnzymeRates.jl`
- Modify: `test/test_dsl.jl`

- [ ] **Step 1: Write a smoke test for the new macro**

```julia
@testset "@allosteric_mechanism (smoke)" begin
    m = @allosteric_mechanism begin
        substrates: F6P
        products:   F16BP
        allosteric_regulators: I::OnlyT

        site(:catalytic, 2): begin
            steps: begin
                [E, F6P] ⇌ [E_F6P]    :: EqualRT
                [E_F6P] <--> [E_F16BP] :: EqualRT
                [E_F16BP] ⇌ [E, F16BP] :: EqualRT
            end
        end
    end
    @test m isa EnzymeRates.AllostericEnzymeMechanism
    @test EnzymeRates.allosteric_regulators(m) ⊇ ((:I, :OnlyT),)
end
```

- [ ] **Step 2: Run — expect failure**

- [ ] **Step 3: Add the macro**

In `src/dsl.jl`, add `macro allosteric_mechanism(block)` and its parser, generating the OLD allosteric type signature (`AllostericEnzymeMechanism{Mets, CM, CS, RS}`) for now. The new signature is introduced in Phase 3.

Parser responsibilities:
- Parse `substrates:`, `products:`, `allosteric_regulators:` (with required `::Tag`), `catalytic_inhibitors:` (no tag).
- Parse `site(:catalytic, N): begin steps: ... end` (required exactly once).
- Parse `site(:regulatory, N): begin ligands: A, I end` (optional, multiple allowed).
- Each step or step-group has a required `:: Tag`.
- Reject single-ligand `::EqualRT` reg sites (validated when site is built).
- Build the inner `EnzymeMechanism` (catalytic mechanism) using the same step-grouping logic as `@enzyme_mechanism`.

Export `@allosteric_mechanism` in `src/EnzymeRates.jl`.

- [ ] **Step 4: Run tests**

Smoke test passes; the old allosteric type is still used internally so the rest of the suite is unaffected.

- [ ] **Step 5: Commit**

```bash
git add src/dsl.jl src/EnzymeRates.jl test/test_dsl.jl
git commit -m "Add @allosteric_mechanism macro skeleton

Generates old AllostericEnzymeMechanism type signature for now;
Phase 3 swaps to the new signature. Parses substrates/products/
allosteric_regulators (with required ::Tag)/catalytic_inhibitors,
site(:catalytic, N): with required tagged steps, optional
site(:regulatory, N): for competing ligands."
```

### Task 2.6: Migrate plain-mechanism call sites

**Files (per migration):**
- Modify: `test/test_dsl.jl`
- Modify: `test/test_types.jl`
- Modify: `test/test_enzyme_derivation.jl`
- Modify: `test/test_mechanism_enumeration.jl`
- Modify: `test/test_fitting.jl`
- Modify: `test/test_identify_rate_equation.jl`
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`
- Modify: `test/test_accessors.jl`
- Modify: `test/test_sym_poly.jl`

For each `@enzyme_mechanism` call site that uses old grammar, rewrite to new grammar.

**Translation rules:**
1. `species: begin substrates: S[C]; products: P[C]; regulators: I; enzymes: E, ES[C] end`
   → `substrates: S` / `products: P` / `regulators: I` (drop atoms; drop `enzymes:` block).
2. `steps: begin [E,S]⇌[ES]; ...; end` stays as-is (each step keeps its expression).
3. `constraints: begin K2 = K1 end` → wrap the corresponding two steps in a parenthesized group: `([E, S] ⇌ [ES], [EP, S] ⇌ [EPS])`. Multiple constraints `K2=K1, K3=K1, K4=K1` become one group containing all four steps.
4. **Delete the three flat-written-out homodimer tests entirely** (decision B): "MWC Dimer" (`test/mechanism_definitions_for_test_enzyme_derivation.jl:~1453`), "Homodimer + Non-competitive Inhibitor" (~1636), "MWC Dimer + Independent Inhibitor" (~1903). Their `[AllostericEnzymeMechanism]` siblings remain (and migrate to the new `@allosteric_mechanism` macro).

- [ ] **Step 1: Find all call sites and triage**

```bash
grep -rln "@enzyme_mechanism\|species: begin\|constraints: begin\|enzymes:" test/ --include="*.jl"
```

- [ ] **Step 2: Per-file migration**

Work through each test file. After each file, run:

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -20
```

Expected: per-file failures decrease as plain mechanisms migrate.

- [ ] **Step 3: Migrate analytical-formula references**

Some tests have hand-derived `analytical_rate_fn(p, c)` formulas referencing `K1`, `k2f`, etc. by step index. Step order is preserved (no canonicalization), so indices stay aligned with the user-DSL order — formulas should still work as long as the step ORDER in the new DSL matches the old. Verify by running the test; if `@test rate_equation(...) ≈ analytical_rate_fn(...)` fails, the migration likely reordered steps. Restore original step order.

- [ ] **Step 4: Delete the flat homodimer tests**

In `test/mechanism_definitions_for_test_enzyme_derivation.jl`, delete the three "flat" homodimer mechanism definitions and their `MechanismTestSpec` entries. Their `[AllostericEnzymeMechanism]` siblings stay (now using `@allosteric_mechanism`).

- [ ] **Step 5: Commit progressively**

One commit per file or logical batch:

```bash
git add test/<file>.jl
git commit -m "Migrate <file> to new @enzyme_mechanism DSL"
```

End state: plain-mechanism tests pass; allosteric tests still failing (handled in Phase 3).

### Task 2.7: Migrate rate-equation derivation for new types

**Files:**
- Modify: `src/rate_eq_derivation.jl`
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl`
- Modify: `src/sym_poly_for_rate_eq_derivation.jl`

The rate-equation machinery currently indexes `Species[k]` and reads `param_constraints(m)`. Update to use new accessors and kinetic-group fan-out via symbol renaming (decision G).

- [ ] **Step 1: Audit magic-index reads**

```bash
grep -n "Species\[\|.parameters\[" src/rate_eq_derivation.jl src/thermodynamic_constr_for_rate_eq_derivation.jl src/sym_poly_for_rate_eq_derivation.jl
```

Replace each with the appropriate accessor: `Species[1]` → `substrates(m)`, etc.

- [ ] **Step 2: Replace constraint fan-out with symbol renaming**

In `src/sym_poly_for_rate_eq_derivation.jl`, delete the `_apply_param_constraints` methods on `POLY` / `FactoredSigma` / `FactoredPoly` / `DenomTerm`. Replace with a single `_rename_symbols(poly, rename_map)` that performs a recursive symbol substitution:

```julia
"""
Rename symbols in a polynomial structure. `rename_map` is a `Dict{Symbol, Symbol}`;
absent keys are left unchanged.
"""
_rename_symbols(p::POLY, m) = ... # apply to monomial keys
_rename_symbols(fp::FactoredPoly, m) = ...
_rename_symbols(fs::FactoredSigma, m) = ...
_rename_symbols(dt::DenomTerm, m) = ...
```

In `src/rate_eq_derivation.jl`, where the old code computed:

```julia
for (target, coeff, factors) in param_constraints(m)
    # build substitution Expr
end
csigma = _apply_param_constraints(csigma, pc; ...)
```

Replace with:

```julia
function _build_kinetic_rename_map(m)
    rename = Dict{Symbol, Symbol}()
    for g in kinetic_groups(m)
        idxs = steps_in_group(m, g)
        length(idxs) == 1 && continue
        rep = first(idxs)
        for idx in idxs
            idx == rep && continue
            if equilibrium_steps(m)[idx]
                rename[Symbol("K$idx")] = Symbol("K$rep")
            else
                rename[Symbol("k$(idx)f")] = Symbol("k$(rep)f")
                rename[Symbol("k$(idx)r")] = Symbol("k$(rep)r")
            end
        end
    end
    rename
end

# Apply to polynomials before Haldane runs:
poly = _rename_symbols(poly, _build_kinetic_rename_map(m))
```

Update `parameters(m)`, `_raw_param_symbols`, and any other helper to skip non-representative steps' parameters (they're aliased to their group representative).

- [ ] **Step 3: Update Haldane derivation in `thermodynamic_constr_for_rate_eq_derivation.jl`**

The Haldane Gaussian-elimination machinery operates on the parameter symbol set. With kinetic-group renaming applied to polynomials, the parameter set passed to Haldane shrinks (one symbol per group). Verify Haldane's `_dependent_param_exprs` still produces correct output.

- [ ] **Step 4: Run tests**

Plain-mechanism rate-equation tests pass. Allosteric tests still failing (handled in Phase 3).

- [ ] **Step 5: Commit**

```bash
git add src/rate_eq_derivation.jl src/thermodynamic_constr_for_rate_eq_derivation.jl src/sym_poly_for_rate_eq_derivation.jl
git commit -m "Rewrite rate-equation derivation for new EnzymeMechanism

Replace Species[k] magic-index reads with accessors. Replace
_apply_param_constraints (general constraint handling) with
_rename_symbols (simple symbol substitution) driven by kinetic-group
representatives. parameters(m) returns only group representatives
(one K or k_f/k_r symbol per kinetic group). Haldane runs on the
reduced parameter set.

Allosteric paths still use old type / old derivation; Phase 3 fixes."
```

---

## Phase 3: New `AllostericEnzymeMechanism` Type and Rate-Equation Derivation

### Task 3.1: Redefine `AllostericEnzymeMechanism` struct + accessors

**Files:**
- Modify: `src/types.jl`
- Modify: `test/test_types.jl`

- [ ] **Step 1: Write failing test for new struct**

```julia
@testset "AllostericEnzymeMechanism (new design)" begin
    cm = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            [E, S] ⇌ [ES]
            [ES] <--> [EP]
            [EP] ⇌ [E, P]
        end
    end
    cat_sites = (2, ((2, :OnlyR),))
    reg_sites = ((((:I,), 2, ((:I, :OnlyT),)),),)
    m = EnzymeRates.AllostericEnzymeMechanism{typeof(cm), cat_sites, reg_sites[1]}()

    @test EnzymeRates.catalytic_mechanism(m) === cm
    @test EnzymeRates.catalytic_multiplicity(m) == 2
    @test EnzymeRates.group_tag(m, 1) == :NonequalRT   # default
    @test EnzymeRates.group_tag(m, 2) == :OnlyR
    @test EnzymeRates.regulatory_sites(m) == reg_sites[1]
end
```

- [ ] **Step 2: Run — expect failure**

- [ ] **Step 3: Replace struct + accessors**

In `src/types.jl`, replace `AllostericEnzymeMechanism` with:

```julia
"""
    AllostericEnzymeMechanism{CatalyticMech, CatSites, RegSites}

Singleton type encoding a multi-subunit MWC allosteric enzyme.

- `CatalyticMech`: an `EnzymeMechanism` type (single-subunit catalytic mech).
- `CatSites`: `(multiplicity::Int, group_tags::Tuple{Pair{Int,Symbol}...})`.
  Non-default-only storage; absent groups have tag `:NonequalRT`.
- `RegSites`: tuple of entries `((ligands, multiplicity, ligand_tags),)`.
  One entry per reg site.
"""
struct AllostericEnzymeMechanism{
    CatalyticMech, CatSites, RegSites,
} <: AbstractEnzymeMechanism end
```

Add accessors:

```julia
catalytic_mechanism(::AllostericEnzymeMechanism{CM}) where {CM} = CM()
catalytic_multiplicity(::AllostericEnzymeMechanism{CM, CS}) where {CM, CS} = CS[1]

function group_tag(m::AllostericEnzymeMechanism, g::Int)
    CS = typeof(m).parameters[2]
    for (k, t) in CS[2]
        k == g && return t
    end
    :NonequalRT
end

step_tag(m::AllostericEnzymeMechanism, idx::Int) =
    group_tag(m, kinetic_group(catalytic_mechanism(m), idx))

substrates(m::AllostericEnzymeMechanism)         = substrates(catalytic_mechanism(m))
products(m::AllostericEnzymeMechanism)           = products(catalytic_mechanism(m))
reactions(m::AllostericEnzymeMechanism)          = reactions(catalytic_mechanism(m))
equilibrium_steps(m::AllostericEnzymeMechanism)  = equilibrium_steps(catalytic_mechanism(m))
n_steps(m::AllostericEnzymeMechanism)            = n_steps(catalytic_mechanism(m))
enzyme_forms(m::AllostericEnzymeMechanism)       = enzyme_forms(catalytic_mechanism(m))
n_states(m::AllostericEnzymeMechanism)           = n_states(catalytic_mechanism(m))
kinetic_group(m::AllostericEnzymeMechanism, i::Int) = kinetic_group(catalytic_mechanism(m), i)
kinetic_groups(m::AllostericEnzymeMechanism)     = kinetic_groups(catalytic_mechanism(m))
steps_in_group(m::AllostericEnzymeMechanism, g)  = steps_in_group(catalytic_mechanism(m), g)
stoich_matrix(m::AllostericEnzymeMechanism)      = stoich_matrix(catalytic_mechanism(m))
enzyme_row_range(m::AllostericEnzymeMechanism)   = enzyme_row_range(catalytic_mechanism(m))
metabolite_row_range(m::AllostericEnzymeMechanism) = metabolite_row_range(catalytic_mechanism(m))

function regulators(m::AllostericEnzymeMechanism)
    cat_regs = regulators(catalytic_mechanism(m))
    RS = typeof(m).parameters[3]
    extra = Symbol[]
    seen = Set{Symbol}(cat_regs)
    for entry in RS
        for lig in entry[1]
            lig in seen || (push!(seen, lig); push!(extra, lig))
        end
    end
    (cat_regs..., extra...)
end

function metabolites(m::AllostericEnzymeMechanism)
    cat_mets = metabolites(catalytic_mechanism(m))
    RS = typeof(m).parameters[3]
    extra = Symbol[]
    seen = Set{Symbol}(cat_mets)
    for entry in RS
        for lig in entry[1]
            lig in seen || (push!(seen, lig); push!(extra, lig))
        end
    end
    (cat_mets..., extra...)
end

function allosteric_regulators(m::AllostericEnzymeMechanism)
    RS = typeof(m).parameters[3]
    result = Tuple{Symbol, Symbol}[]
    for (ligands, _, lig_tags) in RS
        tag_map = Dict(lig_tags)
        for lig in ligands
            push!(result, (lig, get(tag_map, lig, :NonequalRT)))
        end
    end
    Tuple(result)
end

function catalytic_inhibitors(m::AllostericEnzymeMechanism)
    RS = typeof(m).parameters[3]
    rs_names = Set{Symbol}()
    for (ligs, _, _) in RS
        for l in ligs; push!(rs_names, l); end
    end
    cat_regs = regulators(catalytic_mechanism(m))
    Tuple(r for r in cat_regs if r ∉ rs_names)
end

regulatory_sites(::AllostericEnzymeMechanism{CM, CS, RS}) where {CM, CS, RS} = RS
regulatory_site_ligands(m::AllostericEnzymeMechanism, i::Int)     = regulatory_sites(m)[i][1]
regulatory_site_multiplicity(m::AllostericEnzymeMechanism, i::Int) = regulatory_sites(m)[i][2]

function regulatory_ligand_tag(m::AllostericEnzymeMechanism, i::Int, lig::Symbol)
    for (k, t) in regulatory_sites(m)[i][3]
        k == lig && return t
    end
    :NonequalRT
end
```

- [ ] **Step 4: Run tests**

New struct test passes. `@allosteric_mechanism` macro from Task 2.5 still emits old type — the test from Task 2.5 may now fail. Mark RED-OK; Task 3.2 fixes.

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "Redefine AllostericEnzymeMechanism with new shape

CatSites = (multiplicity, group_tags). RegSites entries are
(ligands, multiplicity, ligand_tags). Non-default-only tag storage.
All access via named accessors.

KNOWN RED: @allosteric_mechanism still emits old type."
```

### Task 3.2: Implement `AllostericEnzymeMechanism` constructor and update `@allosteric_mechanism`

**Files:**
- Modify: `src/types.jl`
- Modify: `src/dsl.jl`

- [ ] **Step 1: Write failing test**

```julia
@testset "AllostericEnzymeMechanism constructor + DSL" begin
    cm = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            [E, S] ⇌ [ES]
            [ES] <--> [EP]
            [EP] ⇌ [E, P]
        end
    end

    # Single-ligand :EqualRT reg site → error
    @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
        cm, (2, ()), ((((:I,), 2, ((:I, :EqualRT),)),),)[1],
    )

    # Iso group :OnlyT → error
    @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
        cm, (2, ((2, :OnlyT),)), (),
    )

    # Build via DSL
    m = @allosteric_mechanism begin
        substrates: S
        products:   P
        allosteric_regulators: I::OnlyT

        site(:catalytic, 2): begin
            steps: begin
                [E, S] ⇌ [ES]    :: EqualRT
                [ES] <--> [EP]    :: OnlyR
                [EP] ⇌ [E, P]    :: EqualRT
            end
        end
    end
    @test EnzymeRates.catalytic_multiplicity(m) == 2
    @test EnzymeRates.group_tag(m, 1) == :EqualRT
    @test EnzymeRates.group_tag(m, 2) == :OnlyR
    @test EnzymeRates.allosteric_regulators(m) == ((:I, :OnlyT),)
end
```

- [ ] **Step 2: Run — expect failure**

- [ ] **Step 3: Implement constructor in `src/types.jl`**

```julia
function AllostericEnzymeMechanism(
    cm::EnzymeMechanism,
    cat_sites::Tuple{Int, <:Tuple},
    reg_sites::Tuple,
)
    multiplicity, group_tags = cat_sites
    valid_groups = Set(kinetic_groups(cm))
    eq_steps = equilibrium_steps(cm)
    rxns = reactions(cm)
    cat_mets = Set(metabolites(cm))

    # Validate group tags
    for (g, tag) in group_tags
        g in valid_groups ||
            error("group_tag references non-existent kinetic_group $g")
        tag in (:OnlyR, :OnlyT, :EqualRT, :NonequalRT) ||
            error("Invalid group tag: $tag")
        # Iso-only groups can't be :OnlyT
        any_iso = false
        for idx in steps_in_group(cm, g)
            lhs, rhs, is_eq, _ = rxns[idx]
            mets_in = any(s in cat_mets for s in (lhs..., rhs...))
            if !is_eq && !mets_in
                any_iso = true
            end
        end
        any_iso && tag == :OnlyT &&
            error("Iso group $g tagged :OnlyT is forbidden (R-inactive is a relabel)")
    end

    # Validate reg sites
    for (i, entry) in enumerate(reg_sites)
        ligands, mult, lig_tags = entry
        isempty(ligands) && error("Reg site $i has no ligands")
        tag_map = Dict(lig_tags)
        for (lig, tag) in lig_tags
            tag in (:OnlyR, :OnlyT, :EqualRT, :NonequalRT) ||
                error("Invalid reg-site tag $tag for ligand $lig")
        end
        all_eq = all(get(tag_map, l, :NonequalRT) == :EqualRT for l in ligands)
        all_eq &&
            error("Reg site $i with all `:EqualRT` ligands cancels identically " *
                  "(or single-ligand :EqualRT reg site); at least one ligand " *
                  "must have a non-:EqualRT tag. Ligands: $ligands")
    end

    # Sort group_tags by group number for canonical form
    sorted_tags = Tuple(sort(collect(group_tags); by=first))
    cat_sites_canon = (multiplicity, sorted_tags)

    AllostericEnzymeMechanism{typeof(cm), cat_sites_canon, reg_sites}()
end
```

- [ ] **Step 4: Update `@allosteric_mechanism` macro to emit the new type**

In `src/dsl.jl`, change the macro's expansion to call the new 3-arg `AllostericEnzymeMechanism(cm, cat_sites, reg_sites)` constructor.

The parser:
- Builds `cm` from the `site(:catalytic, N):` block's `steps:` body via the same step-grouping logic as `@enzyme_mechanism`.
- Captures `::Tag` annotations on each step or step-group; aggregates them into `group_tags::Tuple{Pair{Int,Symbol}...}`.
- Builds `reg_sites` from the species declarations (all listed allosteric_regulators) and explicit `site(:regulatory, N):` blocks. Unlisted allosteric regulators each get their own reg site with multiplicity = N (catalytic multiplicity).
- Tag distribution: for each ligand, look up its declared tag from `allosteric_regulators:`. Build `ligand_tags::Tuple{Pair{Symbol,Symbol}...}`.

- [ ] **Step 5: Run tests**

The new tests pass. The Task 2.5 smoke test now produces the new type. Other allosteric tests still fail (rate-eq derivation, kcat — Task 3.3 / 3.4 fix).

- [ ] **Step 6: Commit**

```bash
git add src/types.jl src/dsl.jl test/test_types.jl
git commit -m "Implement AllostericEnzymeMechanism constructor + new DSL emission

Validates group tags, iso-group :OnlyT prohibition, single-/all-EqualRT
reg-site rejection. @allosteric_mechanism now emits the new 3-param
type signature.

KNOWN RED: rate-equation derivation and _kcat_forward use old impl."
```

### Task 3.3: Rewrite allosteric rate-equation derivation

**Files:**
- Modify: `src/rate_eq_derivation.jl`
- Modify: `src/sym_poly_for_rate_eq_derivation.jl`

The current parallel R/T-state derivation is replaced with: derive R-state via the plain `EnzymeMechanism` machinery, then build T-state by symbol substitution at POLY level (zero `:OnlyR` symbols, rename `:NonequalRT` symbols to T-counterparts; for R-state, zero `:OnlyT` symbols).

- [ ] **Step 1: Migrate the existing reference allosteric mechanism tests to `@allosteric_mechanism`**

In `test/mechanism_definitions_for_test_enzyme_derivation.jl`, find the existing `[AllostericEnzymeMechanism]`-suffixed test specs (e.g., `MWC Dimer [AllostericEnzymeMechanism]` at ~1502, the homodimer-noncomp-inh sibling at ~1691, the MWC-dimer-inh sibling at ~1959). Rewrite each using `@allosteric_mechanism` with explicit `:: NonequalRT` tags on every step group (matching their existing default-NonequalRT behavior).

- [ ] **Step 2: Run tests — expect failure for these allosteric mechanisms**

The new type doesn't have rate-equation derivation yet.

- [ ] **Step 3: Delete old allosteric rate-equation code**

In `src/rate_eq_derivation.jl`, delete:
- `_is_tr_equiv_catalytic_K`, `_is_tr_equiv_catalytic_param`, `_is_r_only_catalytic_param`.
- The old `_binding_K_symbols(::Type{<:AllostericEnzymeMechanism{...}})` method.
- The old `_dependent_param_exprs(::Type{<:AllostericEnzymeMechanism{...}})` method (~100 lines).
- `_allosteric_dep_assignments`.
- `_allosteric_num_den_exprs`.
- `_build_allosteric_rate_body`.
- The old `@generated rate_equation(...,::AllostericEnzymeMechanism, ::ReducedMode)`.
- The old `rate_equation_string(::AllostericEnzymeMechanism, ::ReducedMode)`.
- The old `structural_identifiability_deficit(::AllostericEnzymeMechanism)`.

In `src/sym_poly_for_rate_eq_derivation.jl`, delete:
- `_rs_tr_equiv` / `_rs_r_only` / `_rs_t_only`.
- `_count_allosteric_rate_monomials`.

- [ ] **Step 4: Implement new derivation**

Add to `src/rate_eq_derivation.jl`:

```julia
# ═══════════════════════════════════════════════════════════════════
# AllostericEnzymeMechanism rate equation — unified path
# ═══════════════════════════════════════════════════════════════════

"""
    _onlyT_syms(m::AllostericEnzymeMechanism) → Set{Symbol}

Symbols (K, k_f, k_r) that are absent in R-state polynomial because
their kinetic group is tagged :OnlyT, plus reg-site ligand-name
symbols for ligands tagged :OnlyT (their concentration is "absent" from
R-state Q).
"""
function _onlyT_syms(m)
    cm = catalytic_mechanism(m)
    syms = Set{Symbol}()
    for g in kinetic_groups(cm)
        group_tag(m, g) == :OnlyT || continue
        rep = first(steps_in_group(cm, g))
        if equilibrium_steps(cm)[rep]
            push!(syms, Symbol("K$rep"))
        else
            push!(syms, Symbol("k$(rep)f"))
            push!(syms, Symbol("k$(rep)r"))
        end
    end
    # Reg-site OnlyT ligands are absent from R-state Q only on a per-site basis;
    # they are not zeroed at the catalytic-poly level — they're handled in reg_Q construction.
    syms
end

"""
    _onlyR_syms(m) → Set{Symbol}

Symmetric to _onlyT_syms but for :OnlyR-tagged kinetic groups (T-state zeroing).
"""
function _onlyR_syms(m)
    cm = catalytic_mechanism(m)
    syms = Set{Symbol}()
    for g in kinetic_groups(cm)
        group_tag(m, g) == :OnlyR || continue
        rep = first(steps_in_group(cm, g))
        if equilibrium_steps(cm)[rep]
            push!(syms, Symbol("K$rep"))
        else
            push!(syms, Symbol("k$(rep)f"))
            push!(syms, Symbol("k$(rep)r"))
        end
    end
    syms
end

"""
    _nonequalRT_rename(m) → Dict{Symbol, Symbol}

Symbol rename map for T-state derivation: each `:NonequalRT`-tagged group's
representative symbols are renamed to T-suffixed counterparts. `:EqualRT` and
`:OnlyR` groups keep R-state names (EqualRT shares with R; OnlyR has been zeroed
already).
"""
function _nonequalRT_rename(m)
    cm = catalytic_mechanism(m)
    rename = Dict{Symbol, Symbol}()
    for g in kinetic_groups(cm)
        tag = group_tag(m, g)
        tag == :NonequalRT || continue
        rep = first(steps_in_group(cm, g))
        if equilibrium_steps(cm)[rep]
            rename[Symbol("K$rep")] = Symbol("K$(rep)_T")
        else
            rename[Symbol("k$(rep)f")] = Symbol("k$(rep)f_T")
            rename[Symbol("k$(rep)r")] = Symbol("k$(rep)r_T")
        end
    end
    rename
end

"""
    _build_reg_Q(entry, conformation::Symbol) → Expr

Build the reg-site partition function Expr for the given site entry and
conformation (:R or :T).
- :R: skip ligands tagged :OnlyT.
- :T: skip ligands tagged :OnlyR; use T-suffixed K names for :NonequalRT ligands.
"""
function _build_reg_Q(site_idx::Int, entry, conf::Symbol)
    ligs, mult, lig_tags_raw = entry
    lig_tags = Dict(lig_tags_raw)
    terms = Any[1]
    for lig in ligs
        tag = get(lig_tags, lig, :NonequalRT)
        if conf == :R && tag == :OnlyT; continue; end
        if conf == :T && tag == :OnlyR; continue; end
        K_sym = if conf == :T && tag == :NonequalRT
            Symbol("K_$(lig)_T_reg$(site_idx)")
        else
            Symbol("K_$(lig)_reg$(site_idx)")
        end
        push!(terms, :($lig / $K_sym))
    end
    isempty(terms) ? 1 : foldl((a, b) -> :($a + $b), terms)
end

@generated function rate_equation(
    m::AllostericEnzymeMechanism, concs::NamedTuple, params::NamedTuple, ::ReducedMode,
)
    M = m
    cm_type = M.parameters[1]
    cat_n = M.parameters[2][1]
    reg_sites = M.parameters[3]

    # R-state and T-state polynomials via tag-driven substitution
    raw_num, raw_den = _raw_rate_polys(cm_type)
    m_inst = M()

    onlyT = _onlyT_syms(m_inst)
    onlyR = _onlyR_syms(m_inst)
    rename_T = _nonequalRT_rename(m_inst)

    # R-state: zero :OnlyT syms
    num_R = _zero_symbols_in_poly(raw_num, onlyT)
    den_R = _zero_symbols_in_poly(raw_den, onlyT)

    # T-state: zero :OnlyR syms, then rename :NonequalRT syms
    num_T_zeroed = _zero_symbols_in_poly(raw_num, onlyR)
    den_T_zeroed = _zero_symbols_in_poly(raw_den, onlyR)
    num_T = _rename_symbols(num_T_zeroed, rename_T)
    den_T = _rename_symbols(den_T_zeroed, rename_T)

    # Convert to Exprs
    binding_Ks_r = ... # accessor-based
    N_R = _poly_to_expr(num_R, ..., binding_Ks_r)
    Q_R = _poly_to_expr(den_R, ..., binding_Ks_r)
    N_T = _poly_to_expr(num_T, ..., binding_Ks_r)
    Q_T = _poly_to_expr(den_T, ..., binding_Ks_r)

    reg_Q_R = [_build_reg_Q(i, e, :R) for (i, e) in enumerate(reg_sites)]
    reg_Q_T = [_build_reg_Q(i, e, :T) for (i, e) in enumerate(reg_sites)]

    # Assemble num and den (MWC formula)
    # num = cat_n * (N_R * Q_R^(cat_n-1) * Π reg_Q_R^mult + L * N_T * Q_T^(cat_n-1) * Π reg_Q_T^mult)
    # den = Q_R^cat_n * Π reg_Q_R^mult + L * Q_T^cat_n * Π reg_Q_T^mult

    num_expr = ... # build as before, with the new R/T poly Exprs
    den_expr = ...
    rate_expr = :(E_total * ($num_expr) / ($den_expr))

    # Parameter destructuring
    params_list = ... # via parameters(m_inst, ReducedMode())
    Expr(:block,
        _destructuring_expr(params_list, :params),
        _destructuring_expr(metabolites(m_inst), :concs),
        :(return $rate_expr))
end
```

- [ ] **Step 5: Update `rate_equation_string` and `structural_identifiability_deficit` to use the same path**

Both share the substitution machinery. Extract common substring/expression-building into a helper.

- [ ] **Step 6: Run tests**

Allosteric rate-equation tests pass.

- [ ] **Step 7: Commit**

```bash
git add src/rate_eq_derivation.jl src/sym_poly_for_rate_eq_derivation.jl test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -m "Rewrite allosteric rate-equation derivation as POLY-level substitution

R-state and T-state polynomials are produced by symbol substitutions
over a shared raw polynomial (derived once via the plain-EnzymeMechanism
path on the embedded CatalyticMech). For R: zero :OnlyT symbols. For T:
zero :OnlyR symbols, then rename :NonequalRT symbols to T-suffixed
counterparts. Reg-site Q polynomials are built per-site, per-conformation,
filtering ligands by tag.

Deletes _is_tr_equiv_*, _is_r_only_*, _rs_tr_equiv/r_only/t_only,
_count_allosteric_rate_monomials, _allosteric_dep_assignments, and the
parallel R/T control flow throughout."
```

### Task 3.4: Migrate `_kcat_forward(::AllostericEnzymeMechanism)`

**Files:**
- Modify: `src/rate_eq_derivation.jl`
- Modify: `test/test_enzyme_derivation.jl`

The existing `_kcat_forward(::AllostericEnzymeMechanism)` (~150 lines, magic-index access throughout) is replaced with a function that uses the same POLY-level substitution machinery as Task 3.3, plus a shared `_kcat_from_poly(num_poly, den_poly, params)` helper that also serves the plain-mechanism path.

- [ ] **Step 1: Extract `_kcat_from_poly` helper from existing plain-mechanism code**

In `src/rate_eq_derivation.jl`, find the existing `_kcat_forward(m::EnzymeMechanism, params)` body. Factor out the polynomial-traversal kcat extraction into a standalone function:

```julia
"""
    _kcat_from_poly(num_poly, den_poly, params) → Float64

Compute kcat (forward) analytically from the polynomial structure:
group monomials by metabolite pattern, identify the saturation-limit
contribution, return the bounded kcat ratio.
"""
function _kcat_from_poly(num_poly::POLY, den_poly::POLY, params::NamedTuple)
    ... # extracted body
end
```

Update the plain-mechanism `_kcat_forward` to call this helper.

- [ ] **Step 2: Write failing test for allosteric kcat**

```julia
@testset "_kcat_forward(::AllostericEnzymeMechanism)" begin
    m = mwc_dimer_new_dsl   # the migrated reference
    params = (...)            # specific test point
    expected = ...            # hand-derived kcat at this point
    @test EnzymeRates._kcat_forward(m, params) ≈ expected

    # Verify rescale invariant holds
    rescaled = rescale_parameter_values(m, params; kcat=1.0)
    sat_concs = (S=1e6, P=0.0)
    @test rate_equation(m, sat_concs, rescaled) ≈ params.E_total
end
```

- [ ] **Step 3: Run — expect failure**

- [ ] **Step 4: Implement new `_kcat_forward(::AllostericEnzymeMechanism)`**

```julia
function _kcat_forward(m::AllostericEnzymeMechanism, params)
    cm = catalytic_mechanism(m)
    raw_num, raw_den = _raw_rate_polys(typeof(cm))

    onlyT = _onlyT_syms(m)
    onlyR = _onlyR_syms(m)
    rename_T = _nonequalRT_rename(m)

    num_R = _zero_symbols_in_poly(raw_num, onlyT)
    den_R = _zero_symbols_in_poly(raw_den, onlyT)
    num_T = _rename_symbols(_zero_symbols_in_poly(raw_num, onlyR), rename_T)
    den_T = _rename_symbols(_zero_symbols_in_poly(raw_den, onlyR), rename_T)

    kcat_R = _kcat_from_poly(num_R, den_R, params)
    kcat_T = _kcat_from_poly(num_T, den_T, params)

    # Iterate regulator corners: for each subset of regulators "saturating",
    # compute effective L factor and weighted kcat. Return max.
    n_lig = sum(length(regulatory_site_ligands(m, i)) for i in 1:length(regulatory_sites(m)))
    best = max(kcat_R, kcat_T)
    for mask in 0:(2^n_lig - 1)
        # ... iterate corners, compute corner-specific R/T weighting from L and tag-modified Q's
        # update best if a higher kcat is reachable at this corner
    end
    best
end
```

The corner-iteration logic preserves the semantics of the old code: at saturating regulator concentrations, certain regulators are "fully bound" (R-state or T-state, per their tag) and shift the L_effective. The function returns the max kcat across corners.

- [ ] **Step 5: Run tests**

Allosteric kcat test passes; `rescale_parameter_values` invariant holds.

- [ ] **Step 6: Commit**

```bash
git add src/rate_eq_derivation.jl test/test_enzyme_derivation.jl
git commit -m "Migrate _kcat_forward(::AllostericEnzymeMechanism) to new tag system

Uses _kcat_from_poly shared helper (factored from the plain-mechanism
path) and the same POLY-level substitution machinery as the allosteric
rate-equation derivation. Magic-index access (CS[2]/CS[5]/CS[7]) is gone;
all access via accessors. Reduces ~150 lines to ~60."
```

---

## Phase 4: Mechanism Enumeration Simplification

### Task 4.1: Update `MechanismSpec` / `AllostericMechanismSpec`

**Files:**
- Modify: `src/mechanism_enumeration.jl`
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Update `StepSpec`**

```julia
struct StepSpec
    reactants::Vector{Symbol}
    products::Vector{Symbol}
    is_equilibrium::Bool
    kinetic_group::Int
end

Base.:(==)(a::StepSpec, b::StepSpec) =
    a.reactants == b.reactants &&
    a.products == b.products &&
    a.is_equilibrium == b.is_equilibrium &&
    a.kinetic_group == b.kinetic_group

Base.hash(s::StepSpec, h::UInt) =
    hash(s.kinetic_group, hash(s.is_equilibrium,
        hash(s.products, hash(s.reactants, h))))
```

- [ ] **Step 2: Update `MechanismSpec`**

```julia
struct MechanismSpec <: AbstractMechanismSpec
    reaction::Any
    steps::Vector{StepSpec}
    param_count::Int
end
```

Delete the `param_constraints::Vector{ParamConstraint}` field. Delete `ParamConstraint` type alias.

- [ ] **Step 3: Update `AllostericMechanismSpec`**

```julia
struct AllostericMechanismSpec <: AbstractMechanismSpec
    base::MechanismSpec
    catalytic_n::Int
    allosteric_reg_sites::Vector{Vector{Symbol}}
    allosteric_multiplicities::Vector{Int}
    group_tags::Dict{Int, Symbol}           # kinetic_group → tag
    reg_ligand_tags::Dict{Symbol, Symbol}   # ligand → tag
    param_count::Int
end
```

Delete the old fields: `tr_equiv_metabolites`, `tr_equiv_cat_steps`, `r_only_metabolites`, `t_only_metabolites`, `r_only_cat_steps`.

- [ ] **Step 4: Update spec → mechanism constructors**

```julia
function EnzymeMechanism(spec::MechanismSpec)
    rxn = spec.reaction
    subs = Tuple(s for s in substrates(rxn))   # already Symbol tuple at reaction level
    prods = Tuple(s for s in products(rxn))
    regs = Tuple(r for r in regulators(rxn))
    mets = (subs, prods, regs)
    rxns = Tuple(
        (Tuple(s.reactants), Tuple(s.products), s.is_equilibrium, s.kinetic_group)
        for s in spec.steps
    )
    EnzymeMechanism(mets, rxns)
end

function AllostericEnzymeMechanism(spec::AllostericMechanismSpec)
    cm = EnzymeMechanism(spec.base)
    group_tags = Tuple(sort(collect(spec.group_tags); by=first))
    cat_sites = (spec.catalytic_n, group_tags)

    reg_sites = Tuple(
        (Tuple(group),
         spec.allosteric_multiplicities[i],
         Tuple((l, spec.reg_ligand_tags[l])
               for l in group if haskey(spec.reg_ligand_tags, l)))
        for (i, group) in enumerate(spec.allosteric_reg_sites)
    )

    AllostericEnzymeMechanism(cm, cat_sites, reg_sites)
end
```

- [ ] **Step 5: Migrate `init_mechanisms` and existing call sites**

Every `StepSpec(reactants, products, is_eq)` becomes `StepSpec(reactants, products, is_eq, gnum)`. `init_mechanisms` should put mirror dead-end steps in the same kinetic group as their catalytic counterpart (replaces the old `K_mirror = K_catalytic` constraint encoding).

- [ ] **Step 6: Run tests**

Some enumeration tests may still fail due to expansion-move changes (Task 4.2). Mark RED.

- [ ] **Step 7: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Update MechanismSpec / AllostericMechanismSpec for new types

StepSpec gains kinetic_group::Int. MechanismSpec drops ParamConstraint
vector. AllostericMechanismSpec drops 5 TR/r-only fields, gains
group_tags and reg_ligand_tags Dicts. Constructors EnzymeMechanism(spec)
and AllostericEnzymeMechanism(spec) updated. init_mechanisms puts
catalytic + dead-end mirror steps in the same kinetic group.

KNOWN RED: expansion moves not yet updated."
```

### Task 4.2: Unify expansion moves

**Files:**
- Modify: `src/mechanism_enumeration.jl`

- [ ] **Step 1: Unify `_expand_re_to_ss`**

Replace the existing `MechanismSpec` and `AllostericMechanismSpec` methods with a single method parametric over spec type. The expansion converts ONE WHOLE kinetic group from RE to SS (atomic; mirror propagation is implicit).

```julia
function _expand_re_to_ss(spec::AbstractMechanismSpec)
    results = typeof(spec)[]
    steps = _steps(spec)
    groups = Dict{Int, Vector{Int}}()
    for (i, s) in enumerate(steps)
        push!(get!(groups, s.kinetic_group, Int[]), i)
    end
    for (g, idxs) in groups
        all(steps[i].is_equilibrium for i in idxs) || continue
        new_steps = copy(steps)
        for i in idxs
            old = steps[i]
            new_steps[i] = StepSpec(old.reactants, old.products, false, g)
        end
        # +1 net: K (1 param) replaced by k_f, k_r (2 params), shared across group
        push!(results, _with_steps(spec, new_steps, _param_count(spec) + 1))
    end
    results
end

_steps(s::MechanismSpec) = s.steps
_steps(s::AllostericMechanismSpec) = s.base.steps
_param_count(s::AbstractMechanismSpec) = s.param_count

_with_steps(spec::MechanismSpec, new_steps, new_pc) =
    MechanismSpec(spec.reaction, new_steps, new_pc)
_with_steps(spec::AllostericMechanismSpec, new_steps, new_pc) = ...
```

- [ ] **Step 2: Unify other expansion moves**

Apply the same pattern to:
- `_expand_split_kinetic_group` (replaces `_expand_remove_constraint`).
- `_expand_add_dead_end_regulator`.
- `_expand_add_allosteric_regulator`.
- `_expand_change_group_tag` (replaces `_expand_remove_tr_equiv` and `_expand_to_allosteric`).

- [ ] **Step 3: Delete the K-type/V-type hardcoded enumeration**

Delete `_valid_allosteric_differentiations`. Replace with uniform per-group tag enumeration in `_expand_change_group_tag` and `_expand_to_allosteric`.

- [ ] **Step 4: Delete dead helpers**

Delete:
- `_is_mirror_of` (mirror propagation now via kinetic-group atomicity).
- `_constrained_step_indices`.
- `_rewrap_allosteric` (single methods via `_with_steps` dispatch).
- `_tr_equiv_met_delta`.
- The `::AllostericMechanismSpec` variants of every expansion move.

- [ ] **Step 5: Update `_canonicalize!` for new structure**

```julia
function _canonicalize!(spec::MechanismSpec)
    sort!(spec.steps, by=_step_sort_key)
    # Renumber kinetic_groups by first-occurrence order in sorted list
    old_to_new = Dict{Int, Int}()
    next_g = 1
    for (i, s) in enumerate(spec.steps)
        if !haskey(old_to_new, s.kinetic_group)
            old_to_new[s.kinetic_group] = next_g
            next_g += 1
        end
        spec.steps[i] = StepSpec(
            s.reactants, s.products, s.is_equilibrium,
            old_to_new[s.kinetic_group],
        )
    end
    old_to_new
end

function _canonicalize!(spec::AllostericMechanismSpec)
    remap = _canonicalize!(spec.base)
    new_gtags = Dict{Int, Symbol}()
    for (old_g, t) in spec.group_tags
        haskey(remap, old_g) && (new_gtags[remap[old_g]] = t)
    end
    empty!(spec.group_tags); merge!(spec.group_tags, new_gtags)
    sort!.(spec.allosteric_reg_sites)
    if length(spec.allosteric_reg_sites) >= 2
        perm = sortperm(spec.allosteric_reg_sites)
        spec.allosteric_reg_sites .= spec.allosteric_reg_sites[perm]
        spec.allosteric_multiplicities .= spec.allosteric_multiplicities[perm]
    end
    spec
end
```

- [ ] **Step 6: Update `_dedup_key` for new structure**

Make sure `Dict` fields are converted to sorted tuples for deterministic hashing.

- [ ] **Step 7: Run enumeration tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Counts may shift due to (a) per-group tag enumeration replacing K-type/V-type hardcoded subsets, and (b) atomic group RE→SS conversion replacing independent-mirror conversion. Record new counts; verify with Denis.

- [ ] **Step 8: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Unify expansion moves; delete K-type/V-type hardcoding

Each expansion move is now a single method parametric over spec type via
_steps/_with_steps accessors. Mirror propagation is implicit in
kinetic-group atomicity (a group's RE→SS conversion converts every
member). _is_mirror_of, _constrained_step_indices, _rewrap_allosteric,
_tr_equiv_met_delta, _valid_allosteric_differentiations all deleted.
Canonicalization renumbers kinetic_groups by first-occurrence in the
sorted step list."
```

### Task 4.3: Verify mirror+catalytic kinetic-group invariant

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Add unit test asserting init_mechanisms puts catalytic + mirror in same kinetic group**

```julia
@testset "init_mechanisms: mirror steps share kinetic_group with catalytic" begin
    rxn = @enzyme_reaction begin
        substrates: S
        products:   P
        regulators: I
        dead_end_inhibitors: I
    end
    specs = init_mechanisms(rxn)
    for spec in specs
        # Find catalytic and mirror dead-end steps
        ... # assert they share kinetic_group
    end
end
```

- [ ] **Step 2: Run tests**

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Verify init_mechanisms preserves mirror-catalytic kinetic-group sharing"
```

---

## Phase 5: PFK and HK Hand-Verified Tests

### Task 5.1: PFK-1

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`
- Modify: `test/test_enzyme_derivation.jl`

- [ ] **Step 1: Write the PFK mechanism via `@allosteric_mechanism`**

```julia
pfk_mechanism = @allosteric_mechanism begin
    substrates: F6P, ATP
    products:   F16BP, ADP
    allosteric_regulators:
        Pi::EqualRT, ATP::OnlyT, ADP::OnlyR, Citrate::OnlyT, F26BP::NonequalRT

    site(:catalytic, 4): begin
        steps: begin
            ([E, F6P] ⇌ [E_F6P], [E_ATP, F6P] ⇌ [E_F6P_ATP])      :: OnlyR
            ([E, ATP] ⇌ [E_ATP], [E_F6P, ATP] ⇌ [E_F6P_ATP])      :: EqualRT
            [E_F6P_ATP] <--> [E_F16BP_ADP]                         :: EqualRT
            ([E_F16BP_ADP] ⇌ [E_ADP, F16BP], [E_F16BP] ⇌ [E, F16BP]) :: EqualRT
            ([E_F16BP_ADP] ⇌ [E_F16BP, ADP], [E_ADP] ⇌ [E, ADP])     :: EqualRT
        end
    end

    site(:regulatory, 4): begin
        ligands: Pi, ATP
    end
end
```

- [ ] **Step 2: Write hand-derived analytical rate function**

`pfk_rate_analytical(params, concs)` — hand-derive following spec §10.2.1. Key:
- F6P :: OnlyR → no F6P term in T-state Q.
- Pi at reg site 1 (EqualRT) → 1 + Pi/K_Pi in both reg_Q_R and reg_Q_T.
- ATP at reg site 1 (OnlyT) → ATP/K_ATP_T_reg1 only in reg_Q_T.
- ADP, Citrate, F26BP each at own reg sites.

- [ ] **Step 3: Write test**

```julia
@testset "PFK rate equation matches analytical form" begin
    concs = (F6P=0.1, ATP=1.0, F16BP=0.001, ADP=0.1, Pi=1.0, Citrate=0.01, F26BP=0.01)
    params = (...; Keq=1000.0, E_total=1.0, L=1.0)
    @test rate_equation(pfk_mechanism, concs, params) ≈ pfk_rate_analytical(params, concs)

    # F6P → 0 implies rate → 0
    @test rate_equation(pfk_mechanism, merge(concs, (F6P=1e-10,)), params) < 1e-6

    # ADP shifts via R-only stabilization; Citrate / ATP via T-only
    rate_low_ADP  = rate_equation(pfk_mechanism, merge(concs, (ADP=0.0,)), params)
    rate_high_ADP = rate_equation(pfk_mechanism, merge(concs, (ADP=10.0,)), params)
    @test rate_high_ADP > rate_low_ADP

    # F26BP NonequalRT: differential effect when K_F26BP_R ≠ K_F26BP_T
    ...
end
```

- [ ] **Step 4: Run tests**

- [ ] **Step 5: Commit**

```bash
git add test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_enzyme_derivation.jl
git commit -m "Add PFK hand-verified rate-equation test

Tests F6P :: OnlyR (K-type), ATP both substrate and allosteric
regulator with different tags per context, Pi :: EqualRT at reg site
with OnlyT co-ligand (non-cancelling), ADP / Citrate / F26BP at own
reg sites with their respective tags. Validates the full allosteric
rate-equation machinery against a hand-derived analytical form."
```

### Task 5.2: HK

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`
- Modify: `test/test_enzyme_derivation.jl`

- [ ] **Step 1: Write the HK mechanism**

```julia
hk_mechanism = @allosteric_mechanism begin
    substrates: Glucose, ATP
    products:   G6P, ADP
    allosteric_regulators: G6P::OnlyT, Pi::EqualRT
    catalytic_inhibitors:  G6P

    site(:catalytic, 2): begin
        steps: begin
            ([E, Glucose] ⇌ [E_Glc], [E_ATP, Glucose] ⇌ [E_Glc_ATP])  :: EqualRT
            ([E, ATP] ⇌ [E_ATP], [E_Glc, ATP] ⇌ [E_Glc_ATP])           :: EqualRT
            [E_Glc_ATP] <--> [E_G6P_ADP]                                :: EqualRT
            ([E_G6P_ADP] ⇌ [E_ADP, G6P], [E_G6P] ⇌ [E, G6P])           :: EqualRT
            ([E_G6P_ADP] ⇌ [E_G6P, ADP], [E_ADP] ⇌ [E, ADP])           :: EqualRT
            [E_ATP, G6P] ⇌ [E_ATP_G6P]                                  :: EqualRT
            [E_ADP, G6P] ⇌ [E_ADP_G6P]                                  :: EqualRT
        end
    end

    site(:regulatory, 2): begin
        ligands: G6P, Pi
    end
end
```

- [ ] **Step 2: Write `hk_rate_analytical` function**

- [ ] **Step 3: Test that G6P appears in three independent contributions**

- [ ] **Step 4: Run + commit**

### Task 5.3: Single-feature edge-case tests

**Files:**
- Modify: `test/test_enzyme_derivation.jl`

- [ ] **Step 1: Add edge-case tests**

```julia
@testset "Allosteric edge cases" begin
    # OnlyT substrate
    onlyT_sub = @allosteric_mechanism begin
        substrates: S
        products:   P
        site(:catalytic, 2): begin
            steps: begin
                [E, S] ⇌ [ES]  :: OnlyT
                [ES] <--> [EP] :: EqualRT
                [EP] ⇌ [E, P]  :: EqualRT
            end
        end
    end
    # rate → 0 as K_T → ∞
    ...

    # V-type only
    vtype = @allosteric_mechanism begin ... end
    # T-state numerator zero
    ...

    # Single-ligand :EqualRT reg site → error
    @test_throws ErrorException @allosteric_mechanism begin
        substrates: S
        products:   P
        allosteric_regulators: I::EqualRT
        site(:catalytic, 2): begin
            steps: begin
                [E, S] ⇌ [ES]  :: EqualRT
                [ES] <--> [EP] :: EqualRT
                [EP] ⇌ [E, P]  :: EqualRT
            end
        end
    end

    # Atom bracket in mechanism DSL → error
    @test_throws ErrorException @enzyme_mechanism begin
        substrates: S[C]
        ...
    end

    # Stoichiometric infeasibility → error
    @test_throws ErrorException EnzymeRates.EnzymeMechanism(
        ((:S,), (:P, :Q), ()),    # P AND Q listed as products
        (((:E, :S), (:ES,), true, 1),
         ((:ES,), (:EP,), false, 2),
         ((:EP,), (:E, :P), true, 3)),    # only P produced; Q never produced → error
    )

    # same_kinetics group across different metabolites → error
    @test_throws ErrorException EnzymeRates.EnzymeMechanism(
        ((:S, :A), (:P,), ()),
        (((:E, :S), (:ES,), true, 1),
         ((:E, :A), (:EA,), true, 1),    # same group as S-binding but different metabolite → error
         ((:EA,), (:E,), false, 2),     # placeholder
         ((:EA,), (:E, :P), true, 3)),
    )
end
```

- [ ] **Step 2: Run + commit**

---

## Phase 6: Final Cleanup

### Task 6.1: Magic-index audit

**Files:** various src/

- [ ] **Step 1: Grep for magic-index access**

```bash
grep -rn "Species\[\|CS\[\|RS\[\|\.parameters\[" src/ --include="*.jl" | grep -v "^.*#"
```

Expected: zero hits. Fix any remaining hits with appropriate accessor calls.

- [ ] **Step 2: Grep for stale function names**

```bash
grep -rn "_is_tr_equiv_catalytic\|_is_r_only_catalytic\|_rs_tr_equiv\|_rs_r_only\|_rs_t_only\|_rewrap_allosteric\|_tr_equiv_met_delta\|_valid_allosteric_differentiations\|_is_mirror_of\|_constrained_step_indices\|_apply_param_constraints\|_walk_rhs!\|_parse_constraint_rhs" src/ --include="*.jl"
```

Expected: zero hits.

- [ ] **Step 3: Run full suite**

Expected: all green.

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "Final magic-index and stale-symbol audit pass"
```

### Task 6.2: Update SPEC.md

**Files:**
- Modify: `SPEC.md`

- [ ] **Step 1: Remove `compile_mechanism` from exports table** (line ~57)

- [ ] **Step 2: Remove `compile_mechanism` from "Complete Exported Symbol List"** (line ~383)

- [ ] **Step 3: Add `@allosteric_mechanism` to the macros table**

- [ ] **Step 4: Update the example workflows that reference `S[C]` atoms** at the mechanism level — strip atoms (or rewrite to use `@enzyme_reaction`).

- [ ] **Step 5: Commit**

```bash
git add SPEC.md
git commit -m "Update SPEC.md for new macros and refactored type system

- Remove compile_mechanism from exports (it's an internal dispatcher).
- Add @allosteric_mechanism to the macros table.
- Strip atom syntax from mechanism-level DSL examples."
```

### Task 6.3: Update CLAUDE.md

**Files:**
- Modify: `.claude/CLAUDE.md`

- [ ] **Step 1: Update "Key Architecture Decisions"**

Replace old `EnzymeMechanism` / `AllostericEnzymeMechanism` descriptions with new ones (matching spec §4).

- [ ] **Step 2: Update "Regulator representation"**

Rewrite for new `allosteric_regulators` / `catalytic_inhibitors` / DSL split.

- [ ] **Step 3: Remove obsolete sections**

- "Dead-end SS/RE propagation" — replaced by kinetic-group atomicity.
- "Dead-end parameter equivalence constraints" — replaced by kinetic-group encoding.
- References to `RegulatorRole` types.
- (Stale `old_*.jl` references already removed in Task 1.4.)

- [ ] **Step 4: Update "Mechanism enumeration building blocks"**

Match the unified expansion-move structure.

- [ ] **Step 5: Verify "Vmax Normalization" still accurate**

Should reference the new `_kcat_from_poly` shared helper.

- [ ] **Step 6: Commit**

```bash
git add .claude/CLAUDE.md
git commit -m "Update CLAUDE.md for refactored type system and DSL"
```

### Task 6.4: Line-count reduction check

**Files:** none (verification)

- [ ] **Step 1: Compare against baseline**

```bash
wc -l src/types.jl src/dsl.jl src/rate_eq_derivation.jl src/mechanism_enumeration.jl src/sym_poly_for_rate_eq_derivation.jl src/thermodynamic_constr_for_rate_eq_derivation.jl
cat REFACTOR_BASELINE.txt
```

Verify substantial reduction across `src/` per spec §1 / §9. If any file is flat or grew significantly, investigate.

- [ ] **Step 2: Clean up baseline file**

```bash
rm REFACTOR_BASELINE.txt
```

- [ ] **Step 3: Final full-suite check**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all green.

---

## Completion

Final git diff against `origin/main`:
- Substantial code deletion across `src/types.jl`, `src/dsl.jl`, `src/rate_eq_derivation.jl`, `src/mechanism_enumeration.jl`, `src/sym_poly_for_rate_eq_derivation.jl`.
- New PFK + HK test fixtures.
- Public API (`rate_equation`, `rate_equation_string`, `parameters`, `identify_rate_equation`, `fit_rate_equation`) functionally unchanged.
- `compile_mechanism` removed from exports (internal only).
- `@allosteric_mechanism` added to exports.

```bash
git push -u origin allosteric-refactor-spec
gh pr create --title "Mechanism types refactor" --body "$(cat <<'EOF'
## Summary
- Refactor EnzymeMechanism and AllostericEnzymeMechanism with simpler type parameters
- Drop atom tracking at mechanism level; rank-based stoichiometric feasibility check
- Per-kinetic-group TR-mode tagging via @allosteric_mechanism
- Unify step grouping into parenthesized DSL syntax (no constraints/enzymes blocks)
- Eliminate magic-index access; strict accessor-only internal interface
- Substantial code deletion (target ≥30% in touched src/ files)

## Design reference
docs/superpowers/specs/2026-04-23-mechanism-types-refactor-design.md

## Test plan
- [ ] Full test suite passes (Aqua + JET + all unit/integration)
- [ ] PFK hand-verified rate equation matches analytical form
- [ ] HK hand-verified rate equation matches analytical form (G6P in three roles)
- [ ] Stoichiometric feasibility check rejects malformed mechanisms
- [ ] Magic-index grep returns zero hits in src/
- [ ] Line-count reduction targets met

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

---

## Spec-coverage self-check

- **§1 Primary goal (code reduction):** Task 6.4.
- **§2 Motivation:** Phases 2-5 implement.
- **§4.1 EnzymeMechanism type:** Task 2.1, 2.2, 2.3.
- **§4.2 AllostericEnzymeMechanism type:** Task 3.1, 3.2.
- **§5 Accessors:** Task 2.1, 3.1.
- **§6 DSL:** Task 2.4, 2.5, 3.2.
- **§7 Tag semantics:** Task 3.3 embodies; Task 5 tests.
- **§8 Error cases:** Tasks 2.3, 2.4, 3.2 (validation), Task 5.3 (tests).
- **§9 Dead code:** Phase 1 + 2.4 + 2.7 + 3.3 + 3.4 + 4.2.
- **§10.1 DSL tests:** Tasks 2.4, 2.5, 3.2; Task 2.6 migrates existing.
- **§10.2 PFK / HK / edge cases:** Tasks 5.1, 5.2, 5.3.
- **§10.3 Enumeration invariants:** Task 4.2 step 7 + Task 4.3.
- **§10.4 kcat invariants:** Task 3.4.
- **§10.5 Aqua/JET:** ongoing across tasks.
- **§10.6 Graphs.jl removal:** Task 1.1 step 3.
- **§11 Migration notes:** Phase 1 + 2.6 + Phase 3 spread.
- **§12 Out-of-scope:** Confirmed not touching listed APIs.
