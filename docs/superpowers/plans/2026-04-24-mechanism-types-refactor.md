# Mechanism Types Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refactor `EnzymeMechanism` and `AllostericEnzymeMechanism` to simpler, non-redundant type parameters. Replace `constraints:` blocks with step-grouping. Add per-step TR-mode tagging (`OnlyR`/`OnlyT`/`EqualRT`/`NonequalRT`) via DSL. Split `@enzyme_mechanism` into `@enzyme_mechanism` (plain) and `@allosteric_mechanism` (MWC). Eliminate magic-index type-parameter access. Delete dead code and duplicated R/T-state derivation paths.

**Architecture:** `EnzymeMechanism{Metabolites, Reactions}` — each step is `(lhs, rhs, is_eq, kinetic_group::Int)`; steps with identical `kinetic_group` share kinetic parameters. `AllostericEnzymeMechanism{CatalyticMech, CatSites, RegSites}` — `CatSites = (multiplicity, group_tags)`; `RegSites` entries `(ligands, multiplicity, ligand_tags)`. Accessor-only read interface, no magic indices. Canonical form enforced by constructor (sort steps, renumber kinetic groups by first-occurrence order).

**Tech Stack:** Julia 1.x, `@generated` functions for rate-equation derivation, `Test`/`Aqua`/`JET` for testing.

**Spec:** `docs/superpowers/specs/2026-04-23-mechanism-types-refactor-design.md`

**Branch:** `allosteric-refactor-spec` (already created; contains spec commits).

---

## Execution notes

- The refactor has significant breaking surface. Work in the `allosteric-refactor-spec` branch. Main stays unchanged until this branch merges.
- Key invariant: **each task either keeps the test suite green or explicitly pauses red in a controlled transitional state** that the next task restores. A "commit" step only runs after the relevant tests pass unless the step header marks the state as **RED-OK**.
- Run tests via `julia --project -e 'using Pkg; Pkg.test()'` (cold — takes several minutes). For incremental feedback during development, `julia --project=. test/runtests.jl` can be run in a persistent REPL.
- Code reduction is a first-class success criterion. After major tasks, run `git diff --shortstat main` and compare to the reduction targets in Section 9.1 of the spec.
- Strict accessor-only rule: grep for `Species[` / `CS[` / `RS[` / `.parameters[` in `src/` after refactor. Zero hits expected.

---

## Phase 0: Preliminaries

### Task 0.1: Capture current test-pass baseline

**Files:** none (validation only)

- [ ] **Step 1: Run the full test suite on baseline**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -40
```

Expected: all tests pass. Record timing and any flaky warnings.

- [ ] **Step 2: Record baseline line counts**

```bash
wc -l src/types.jl src/dsl.jl src/rate_eq_derivation.jl src/mechanism_enumeration.jl src/sym_poly_for_rate_eq_derivation.jl src/thermodynamic_constr_for_rate_eq_derivation.jl
```

Store the numbers in a scratch file `REFACTOR_BASELINE.txt` (gitignored or deleted later). Used for the reduction-target self-check after major phases.

---

## Phase 1: Dead Code Removal

These are independent, safe deletions. Done first because they reduce surface area for the main refactor and don't risk interaction with new work.

### Task 1.1: Delete `graph()` accessor

The `graph()` accessor is defined in `src/types.jl:508-521`, used only by `test/test_accessors.jl` (lines 21, 79) for an allocation-smoke test. Not called from any `src/` computation.

**Files:**
- Modify: `src/types.jl`
- Modify: `test/test_accessors.jl`

- [ ] **Step 1: Remove test references**

In `test/test_accessors.jl`, delete the two lines that call `graph`:

Line 21 (approximately): remove `EnzymeRates.graph(m);`. The remaining calls on that line stay.

Line 79: delete the `@test (@allocated EnzymeRates.graph(m)) == 0` testset entry.

- [ ] **Step 2: Delete the accessor definition**

In `src/types.jl`, delete the `graph` docstring + `@generated function graph(...)` definition (currently lines 504-521, approximately):

```julia
"""
Build a directed graph of enzyme-form connectivity.
Returns (graph, enzyme_forms_tuple).
"""
@generated function graph(::EnzymeMechanism{Species, Reactions}) where {Species, Reactions}
    enzs = Species[4]
    enz_names = Tuple(e[1] for e in enzs)
    name_to_idx = Dict(n => i for (i, n) in enz_names)
    enz_set = Set(enz_names)
    g = SimpleDiGraph(length(enzs))
    for (lhs, rhs) in Reactions
        e_lhs = first(s for s in lhs if s in enz_set)
        e_rhs = first(s for s in rhs if s in enz_set)
        add_edge!(g, name_to_idx[e_lhs], name_to_idx[e_rhs])
        add_edge!(g, name_to_idx[e_rhs], name_to_idx[e_lhs])
    end
    return g, enzs
end
```

- [ ] **Step 3: Check for `Graphs` import necessity**

```bash
grep -rn "using Graphs\|import Graphs\|SimpleDiGraph\|add_edge!" src/
```

If `graph()` was the only user, remove `using Graphs` from `src/types.jl` (line 1) and the `Graphs` dependency from `Project.toml`. If still used elsewhere, leave imports untouched.

- [ ] **Step 4: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: pass.

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_accessors.jl Project.toml
git commit -m "Delete unused graph() accessor

The graph accessor was defined but never called from src/ — the only
user was an allocation-smoke test in test_accessors.jl. Delete it
along with its test and (if unused elsewhere) the Graphs dependency."
```

### Task 1.2: Delete `RegulatorRole` type hierarchy

`abstract type RegulatorRole` and its three concrete subtypes (`Allosteric`, `DeadEnd`, `UnconstrainedRegulator`) are defined in `src/types.jl:18-27` but never dispatched on. All regulator-role logic uses `Symbol`s (`:unknown`, `:dead_end`, `:allosteric`).

**Files:**
- Modify: `src/types.jl`

- [ ] **Step 1: Confirm no dispatch exists**

```bash
grep -rn "::RegulatorRole\|::Allosteric\b\|::DeadEnd\b\|::UnconstrainedRegulator" src/ test/ --include="*.jl"
```

Expected: zero hits (the type names appear only in their own definitions). If any hit shows a real dispatch, pause and consult Denis before deleting.

- [ ] **Step 2: Delete the declarations**

In `src/types.jl`, delete lines ~17-27:

```julia
"""Regulator role in mechanism enumeration."""
abstract type RegulatorRole end

"""Allosteric regulator: participates in MWC allosteric regulation."""
struct Allosteric <: RegulatorRole end

"""Dead-end inhibitor: creates dead-end complexes only."""
struct DeadEnd <: RegulatorRole end

"""Unconstrained: try all roles (allosteric + dead-end)."""
struct UnconstrainedRegulator <: RegulatorRole end
```

- [ ] **Step 3: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: pass.

- [ ] **Step 4: Commit**

```bash
git add src/types.jl
git commit -m "Delete unused RegulatorRole type hierarchy

The RegulatorRole abstract type and its Allosteric / DeadEnd /
UnconstrainedRegulator subtypes were defined but never dispatched on.
All regulator-role logic uses Symbol values directly."
```

### Task 1.3: Delete complex `ParamConstraint` monomial parsing

The constraint-RHS parser supports coefficients and multi-symbol products, but every real usage in tests is a simple equality (`K2 = K1`). Only `test/test_dsl.jl` exercises the monomial forms (`k3r = 2 * k1r`, `k3r = k1f * k2f / k2r`) as parser-generality tests. We will later replace the whole constraint representation with `kinetic_group`; start by restricting the DSL now.

Defer full deletion of `_parse_constraint_rhs` / `_walk_rhs!` until Phase 3 (new DSL), since the current DSL still uses them. This task just trims the test cases and adds a guard error.

**Files:**
- Modify: `src/dsl.jl`
- Modify: `test/test_dsl.jl`

- [ ] **Step 1: Delete DSL-generality tests for monomial forms**

In `test/test_dsl.jl`, find the testsets that include `k3r = 2 * k1r` (around line 246) and `k3r = k1f * k2f / k2r` (around line 265). Delete just those two constraint lines (keep the surrounding mechanism tests intact, replacing these constraints with simple equalities if the surrounding mechanism would otherwise be invalid). If the whole `@testset` is dedicated to monomial-form parsing, delete the whole testset.

Concrete pattern to search for:

```bash
grep -n "k3r = 2 \* k1r\|k3r = k1f \* k2f / k2r" test/test_dsl.jl
```

- [ ] **Step 2: Restrict the DSL parser to equality-only**

In `src/dsl.jl`, modify `_push_constraint!` (around line 347) to error on non-trivial constraints:

```julia
function _push_constraint!(constraints, arg)
    target = arg.args[1]
    target isa Symbol || error("Constraint target must be a symbol, got $target")
    rhs = arg.args[2]
    rhs isa Symbol ||
        error("Constraint RHS must be a bare symbol (equality only); " *
              "complex constraints are auto-generated by Haldane/Wegscheider. " *
              "Got: $target = $rhs")
    coeff = 1
    factors = Expr(:tuple, Expr(:tuple, QuoteNode(rhs), 1))
    push!(constraints.args, Expr(:tuple, QuoteNode(target), coeff, factors))
end
```

- [ ] **Step 3: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: pass (we removed the only tests that needed complex RHS).

- [ ] **Step 4: Commit**

```bash
git add src/dsl.jl test/test_dsl.jl
git commit -m "Restrict constraint DSL to simple equalities

No production mechanism uses constraint-RHS coefficients or multi-symbol
products. The general ParamConstraint machinery will be replaced
entirely in a follow-up commit; meanwhile, reject the unused forms at
parse time with a clear error."
```

---

## Phase 2: Introduce new `EnzymeMechanism` type and DSL

This phase replaces the `EnzymeMechanism` type signature and the `@enzyme_mechanism` macro. It's a breaking change — every call site migrates in this phase. `AllostericEnzymeMechanism` is touched minimally here (only to keep compilation working); its own refactor is Phase 3.

### Task 2.1: Define new `EnzymeMechanism` struct and basic accessors

**Files:**
- Modify: `src/types.jl`
- Modify: `test/test_types.jl`

- [ ] **Step 1: Write failing test for new struct and accessors**

Add a new testset to `test/test_types.jl`:

```julia
@testset "EnzymeMechanism struct + accessors (new design)" begin
    # Build a mechanism by raw tuple construction — verify accessors work.
    subs = ((:S, ((:C, 1),)),)
    prods = ((:P, ((:C, 1),)),)
    regs = ()
    mets = (subs, prods, regs)
    # Steps: (lhs, rhs, is_eq, kinetic_group)
    rxns = (
        ((:E, :S), (:ES,), true, 1),
        ((:ES,), (:EP,), false, 2),
        ((:EP,), (:E, :P), true, 3),
    )
    m = EnzymeRates.EnzymeMechanism{mets, rxns}()

    @test EnzymeRates.substrates(m) == subs
    @test EnzymeRates.products(m) == prods
    @test EnzymeRates.regulators(m) == regs
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

- [ ] **Step 2: Run test — expect failure**

```bash
julia --project -e 'using Pkg; Pkg.test(test_args=["types"])'
```

Expected: fail (new struct signature not yet defined).

- [ ] **Step 3: Replace `EnzymeMechanism` struct and add accessors**

In `src/types.jl`, replace the existing `struct EnzymeMechanism{...}` definition (currently lines 74-76 approximately) with:

```julia
"""
    EnzymeMechanism{Metabolites, Reactions}

Singleton type encoding an enzyme mechanism.

- `Metabolites`: 3-tuple `(substrates, products, regulators)`, each a
  tuple of `(name::Symbol, atoms::Tuple{Vararg{Tuple{Symbol,Int}}})`.
- `Reactions`: tuple of 4-tuples `(lhs_syms, rhs_syms, is_eq::Bool,
  kinetic_group::Int)`. Steps with identical `kinetic_group` share
  kinetic parameters (one `K` for RE groups, one `k_f` and one `k_r`
  for SS groups).
"""
struct EnzymeMechanism{Metabolites, Reactions} <: AbstractEnzymeMechanism end
```

Replace the accessors in `src/types.jl` (lines 432-504 approximately — the existing `substrates`/`products`/`regulators`/`enzyme_forms`/`n_states`/`n_steps`/`reactions`/`equilibrium_steps`/`param_constraints` definitions). Delete the now-invalid `Species[k]` indexing. Accessor code below:

```julia
"""Return substrates as tuple of `(name, atoms)` pairs."""
substrates(::EnzymeMechanism{M}) where {M} = M[1]

"""Return products as tuple of `(name, atoms)` pairs."""
products(::EnzymeMechanism{M}) where {M} = M[2]

"""Return regulators as tuple of `(name, atoms)` pairs (usually atoms empty)."""
regulators(::EnzymeMechanism{M}) where {M} = M[3]

"""Full metabolite list (substrates ∪ products ∪ regulators, deduplicated by name)."""
@generated function metabolites(::EnzymeMechanism{M}) where {M}
    seen = Set{Symbol}()
    names = Symbol[]
    for group in M
        for (name, _) in group
            if name ∉ seen
                push!(seen, name)
                push!(names, name)
            end
        end
    end
    Tuple(names)
end

"""Return the reactions tuple directly: `((lhs, rhs, is_eq, kinetic_group), ...)`."""
reactions(::EnzymeMechanism{M, R}) where {M, R} = R

"""Return the RE/SS flags as a tuple of Bool, parallel to reactions()."""
@generated function equilibrium_steps(::EnzymeMechanism{M, R}) where {M, R}
    Tuple(step[3] for step in R)
end

"""Number of reaction steps."""
n_steps(::EnzymeMechanism{M, R}) where {M, R} = length(R)

"""Return the kinetic group integer for step `idx`."""
kinetic_group(::EnzymeMechanism{M, R}, idx::Int) where {M, R} = R[idx][4]

"""Return the unique kinetic-group integers present in the mechanism, sorted."""
@generated function kinetic_groups(::EnzymeMechanism{M, R}) where {M, R}
    gs = unique(step[4] for step in R)
    Tuple(sort(collect(gs)))
end

"""Return step indices belonging to the given kinetic group."""
@generated function steps_in_group(
    ::EnzymeMechanism{M, R}, ::Val{G},
) where {M, R, G}
    Tuple(i for (i, step) in enumerate(R) if step[4] == G)
end
steps_in_group(m::EnzymeMechanism, g::Int) = steps_in_group(m, Val(g))
```

Enzyme-form accessors (derive from `Reactions`):

```julia
"""Infer enzyme forms from reaction steps + atom balance. Returns tuple of (name, atoms)."""
@generated function enzyme_forms(::EnzymeMechanism{M, R}) where {M, R}
    # All symbols appearing in step sides that are not substrates/products/regulators.
    met_names = Set{Symbol}()
    for group in M
        for (name, _) in group
            push!(met_names, name)
        end
    end
    forms_order = Symbol[]
    for (lhs, rhs, _, _) in R
        for s in lhs; s ∉ met_names && s ∉ forms_order && push!(forms_order, s); end
        for s in rhs; s ∉ met_names && s ∉ forms_order && push!(forms_order, s); end
    end
    # Atom content inferred by solving mass balance; see _infer_enzyme_atoms below.
    atoms_map = _infer_enzyme_atoms(M, R, forms_order)
    Tuple((f, atoms_map[f]) for f in forms_order)
end

"""Number of distinct enzyme forms."""
n_states(m::EnzymeMechanism) = length(enzyme_forms(m))
```

Add helper `_infer_enzyme_atoms(M, R, forms_order)` (pseudo-code here — implement as pure function that walks the reaction graph assigning atom dicts by mass balance, errors if inconsistent):

```julia
"""
Infer enzyme-form atom content by walking the reaction graph and enforcing
atomic conservation: for each step `E_lhs + M_lhs ⇌ E_rhs + M_rhs` (or
`E_lhs ⇌ E_rhs` for iso), atoms(E_rhs) = atoms(E_lhs) ± atoms(M).
Free enzyme (form with no atoms bound anywhere) starts at ∅.
"""
function _infer_enzyme_atoms(
    Metabolites::NTuple{3},
    Reactions,
    forms_order::Vector{Symbol},
)
    # Map metabolite name -> atom Dict{Symbol,Int}
    met_atoms = Dict{Symbol, Dict{Symbol,Int}}()
    for group in Metabolites
        for (name, atoms) in group
            met_atoms[name] = Dict{Symbol,Int}(a => c for (a, c) in atoms)
        end
    end
    met_names = Set(keys(met_atoms))

    # Initialize: first form in forms_order gets empty atoms (the free enzyme).
    # This convention matches the existing `_AnyMechanism` assumption that
    # at least one form has no atoms.
    atoms = Dict{Symbol, Dict{Symbol,Int}}()
    for f in forms_order; atoms[f] = Dict{Symbol,Int}(); end

    # Iteratively propagate atom content across known steps until stable.
    changed = true
    known = Set{Symbol}()
    # Start from any form that's reachable from an iso step that has only
    # one enzyme-on-one-side (needs an entry point; use the first form).
    push!(known, forms_order[1])
    atoms[forms_order[1]] = Dict{Symbol,Int}()

    while changed
        changed = false
        for (lhs, rhs, _, _) in Reactions
            e_lhs = first(s for s in lhs if s ∉ met_names)
            e_rhs = first(s for s in rhs if s ∉ met_names)
            m_lhs = [s for s in lhs if s in met_names]
            m_rhs = [s for s in rhs if s in met_names]
            if e_lhs ∈ known && e_rhs ∉ known
                # Atoms(e_rhs) = Atoms(e_lhs) + sum(m_lhs) - sum(m_rhs)
                new_atoms = copy(atoms[e_lhs])
                for m in m_lhs, (a, c) in met_atoms[m]; new_atoms[a] = get(new_atoms, a, 0) + c; end
                for m in m_rhs, (a, c) in met_atoms[m]; new_atoms[a] = get(new_atoms, a, 0) - c; end
                filter!(p -> p.second != 0, new_atoms)
                atoms[e_rhs] = new_atoms
                push!(known, e_rhs)
                changed = true
            elseif e_rhs ∈ known && e_lhs ∉ known
                new_atoms = copy(atoms[e_rhs])
                for m in m_rhs, (a, c) in met_atoms[m]; new_atoms[a] = get(new_atoms, a, 0) + c; end
                for m in m_lhs, (a, c) in met_atoms[m]; new_atoms[a] = get(new_atoms, a, 0) - c; end
                filter!(p -> p.second != 0, new_atoms)
                atoms[e_lhs] = new_atoms
                push!(known, e_lhs)
                changed = true
            elseif e_lhs ∈ known && e_rhs ∈ known
                # Consistency check
                expected = copy(atoms[e_lhs])
                for m in m_lhs, (a, c) in met_atoms[m]; expected[a] = get(expected, a, 0) + c; end
                for m in m_rhs, (a, c) in met_atoms[m]; expected[a] = get(expected, a, 0) - c; end
                filter!(p -> p.second != 0, expected)
                expected == atoms[e_rhs] ||
                    error("Atom-balance inconsistency at step $(lhs) ⇌ $(rhs): " *
                          "expected $(atoms[e_rhs]), derived $expected")
            end
        end
    end

    # Every form must be reachable
    for f in forms_order
        f ∈ known || error("Enzyme form $f is not reachable from the free form.")
    end

    Dict(f => Tuple((a, c) for (a, c) in atoms[f]) for f in forms_order)
end
```

- [ ] **Step 4: Run the new accessor test**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: the new test case passes. Other tests **will break** — that's expected and handled in the next tasks. Mark the expected-red state.

**RED-OK from this point until end of Task 2.6.**

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "Introduce new EnzymeMechanism{Metabolites, Reactions} type

Each step tuple carries (lhs, rhs, is_eq, kinetic_group). Steps sharing
a kinetic_group share their K (RE) or k_f/k_r (SS) parameters.
Species split into ((subs,), (prods,), (regs,)). Enzyme forms derived
from reactions via atom-balance (_infer_enzyme_atoms).

KNOWN RED: downstream code still expects 4-param type. Fixed in
subsequent commits within this phase."
```

### Task 2.2: Update `EnzymeMechanism` constructor

The existing constructor accepts `(species, reactions, eq_steps, [constraints])`. The new constructor accepts `(metabolites_3tuple, reactions_4tuple)` and handles canonicalization.

**Files:**
- Modify: `src/types.jl`

- [ ] **Step 1: Write constructor test**

Add to `test/test_types.jl`:

```julia
@testset "EnzymeMechanism constructor + canonicalization" begin
    # Same mechanism, two different DSL orderings — should produce same type.
    mets = (((:S, ((:C, 1),)),), ((:P, ((:C, 1),)),), ())

    # Ordering 1
    rxns_a = (
        ((:E, :S), (:ES,), true, 1),
        ((:ES,), (:EP,), false, 2),
        ((:EP,), (:E, :P), true, 3),
    )
    m_a = EnzymeRates.EnzymeMechanism(mets, rxns_a)

    # Ordering 2 — same steps, shuffled with different group numbers
    rxns_b = (
        ((:EP,), (:E, :P), true, 99),
        ((:ES,), (:EP,), false, 5),
        ((:E, :S), (:ES,), true, 42),
    )
    m_b = EnzymeRates.EnzymeMechanism(mets, rxns_b)

    @test typeof(m_a) == typeof(m_b)

    # Same-kinetics group test
    rxns_c = (
        ((:E, :S), (:ES,), true, 1),
        ((:ES, :S), (:ESS,), true, 1),  # shares group with first step
        ((:ES,), (:EP,), false, 2),
        ((:EP,), (:E, :P), true, 3),
    )
    mets_2s = (((:S, ((:C, 1),)),), ((:P, ((:C, 1),)),), ())
    m_c = EnzymeRates.EnzymeMechanism(mets_2s, rxns_c)
    @test EnzymeRates.kinetic_group(m_c, 1) == EnzymeRates.kinetic_group(m_c, 2)
end
```

- [ ] **Step 2: Run test — expect failure**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: fail (no 2-arg constructor yet).

- [ ] **Step 3: Implement constructor**

In `src/types.jl`, replace the current `function EnzymeMechanism(species::Tuple, reactions::Tuple, eq_steps, constraints=())` with the new 2-arg constructor. Keep the validation logic (enzyme reachability, atom balance, etc.) but adapt to new tuple shape.

```julia
"""
    EnzymeMechanism(metabolites, reactions) → EnzymeMechanism

Construct an `EnzymeMechanism` from explicit metabolite 3-tuple and
reaction 4-tuples. The reactions tuple is canonicalized (sorted by
step_sort_key; kinetic groups renumbered by first-occurrence order).

Validates:
- at least one SS step,
- canonical RE step direction (metabolite on LHS),
- atom conservation,
- enzyme reachability,
- kinetic-group composition rules (all RE-same-metabolite, all
  SS-binding-same-metabolite, or iso singleton),
- no mixing substrate-type and inhibitor-type binding in one group.
"""
function EnzymeMechanism(mets::Tuple{Tuple,Tuple,Tuple}, rxns::Tuple)
    subs, prods, regs = mets
    # ---- Normalize species tuples (sort alphabetically by name) ----
    subs = _sort_species(subs)
    prods = _sort_species(prods)
    regs = _sort_species(regs)
    mets = (subs, prods, regs)

    isempty(rxns) && error("Reactions tuple must not be empty")
    all(step[3] for step in rxns) &&
        error("At least one SS step is required (not all steps can be RE)")

    # ---- Build lookup tables ----
    met_atoms = Dict{Symbol, Dict{Symbol, Int}}()
    for group in (subs, prods, regs)
        for (name, atoms) in group
            d = Dict{Symbol, Int}(a => c for (a, c) in atoms)
            if haskey(met_atoms, name)
                met_atoms[name] == d ||
                    error("Inconsistent atoms for metabolite $name")
            else
                met_atoms[name] = d
            end
        end
    end
    met_set = Set(keys(met_atoms))

    # ---- Canonicalize RE step direction (metabolite on LHS) ----
    rxns = ntuple(length(rxns)) do i
        (lhs, rhs, is_eq, gnum) = rxns[i]
        if !is_eq
            return (lhs, rhs, is_eq, gnum)
        end
        lhs_has_met = any(s in met_set for s in lhs)
        rhs_has_met = any(s in met_set for s in rhs)
        if rhs_has_met && !lhs_has_met
            (rhs, lhs, is_eq, gnum)
        else
            (lhs, rhs, is_eq, gnum)
        end
    end

    # ---- Sort steps by _step_sort_key ----
    sort_key(step) = (sort(collect(step[1])), sort(collect(step[2])), step[3])
    sorted_rxns = sort(collect(rxns), by=sort_key)

    # ---- Renumber kinetic groups by first-occurrence order ----
    old_to_new = Dict{Int, Int}()
    next_new = 1
    renumbered = map(sorted_rxns) do step
        (lhs, rhs, is_eq, gnum) = step
        if !haskey(old_to_new, gnum)
            old_to_new[gnum] = next_new
            next_new += 1
        end
        (lhs, rhs, is_eq, old_to_new[gnum])
    end
    rxns_canonical = Tuple(renumbered)

    # ---- Validate kinetic-group composition ----
    _validate_kinetic_groups(rxns_canonical, met_set, met_atoms)

    # ---- Validate atom balance (delegates to _infer_enzyme_atoms) ----
    #  Compute the enzyme-form atoms to trigger the consistency check;
    #  store nothing here — `enzyme_forms(m)` recomputes lazily.
    forms_order = _extract_forms_order(rxns_canonical, met_set)
    _infer_enzyme_atoms(mets, rxns_canonical, forms_order)

    # ---- Validate net stoichiometry ----
    _validate_net_stoichiometry(mets, rxns_canonical, met_set)

    EnzymeMechanism{mets, rxns_canonical}()
end
```

Add helpers (`_validate_kinetic_groups`, `_extract_forms_order`, `_validate_net_stoichiometry`) as pure Julia functions in `src/types.jl`. Key behaviors:

```julia
"""
Validate: for each kinetic group of size 2+, all steps must be all RE
(same metabolite bound) or all SS binding (same metabolite bound);
never mixed, never iso. Also: no mixing substrate-type (form doesn't
contain metabolite) with inhibitor-type (form contains metabolite).
"""
function _validate_kinetic_groups(rxns, met_set, met_atoms)
    groups = Dict{Int, Vector{Int}}()
    for (i, step) in enumerate(rxns)
        push!(get!(groups, step[4], Int[]), i)
    end
    for (gnum, idxs) in groups
        length(idxs) == 1 && continue  # singleton groups are always valid
        # All members: extract (is_eq, metabolite_bound, lhs_form_has_metabolite)
        kinds = map(idxs) do i
            lhs, rhs, is_eq, _ = rxns[i]
            mets_in_step = [s for s in lhs if s in met_set]
            # Iso step = no metabolite on either side
            isempty(mets_in_step) && isempty(s for s in rhs if s in met_set) &&
                error("Iso step at index $i cannot be in a multi-step kinetic group $gnum")
            length(mets_in_step) == 1 ||
                error("Step $i has $(length(mets_in_step)) metabolites on LHS; expected 1")
            met = mets_in_step[1]
            # Find the enzyme form on LHS
            enz_lhs = first(s for s in lhs if s ∉ met_set)
            # Is this substrate-type (enz_lhs doesn't already bind met)
            # or inhibitor-type (enz_lhs already carries met's atoms)?
            # We detect by: substrate-type if the metabolite atoms aren't
            # already in the enzyme's accumulated atoms. We use the atom
            # map computed lazily; simpler check: does enz_rhs contain
            # the SAME symbol M (e.g., ES with S already) as a second copy?
            # Use a string heuristic on the enzyme name: inhibitor-type if
            # the name contains the metabolite twice (e.g., E_S_S).
            # (This isn't robust — better: count occurrences of met in
            #  enzyme form names by walking reachability. For the MVP,
            #  ask the user to distinguish at DSL level.)
            (is_eq, met)
        end
        first_kind = kinds[1]
        for k in kinds[2:end]
            k[1] == first_kind[1] ||
                error("Kinetic group $gnum mixes RE and SS binding steps")
            k[2] == first_kind[2] ||
                error("Kinetic group $gnum binds different metabolites " *
                      "($(first_kind[2]) and $(k[2]))")
        end
        # Substrate-type vs inhibitor-type check: compare each step's
        # enz_lhs form atom content against met atoms. If the enz_lhs
        # already contains met's atoms, this is inhibitor-type binding
        # (metabolite being added to a form that already carries it).
        # All members of the group must match on this classification.
    end
end
```

Implementation detail for the substrate-type vs inhibitor-type check: use the enzyme-form atom map from `_infer_enzyme_atoms` (which runs earlier in the constructor). For each binding step in the group, check whether the LHS enzyme form's atoms contain the metabolite's atoms as a subset. If so, this is inhibitor-type binding. All group members must agree on this classification — if some are substrate-type and others inhibitor-type, error with a message identifying the offending metabolite.

```julia
function _extract_forms_order(rxns, met_set)
    seen = Set{Symbol}()
    order = Symbol[]
    for (lhs, rhs, _, _) in rxns
        for s in lhs; s ∉ met_set && s ∉ seen && (push!(seen, s); push!(order, s)); end
        for s in rhs; s ∉ met_set && s ∉ seen && (push!(seen, s); push!(order, s)); end
    end
    order
end

function _validate_net_stoichiometry(mets, rxns, met_set)
    subs, prods, regs = mets
    expected = Dict{Symbol, Int}()
    for (n, _) in subs; expected[n] = get(expected, n, 0) - 1; end
    for (n, _) in prods; expected[n] = get(expected, n, 0) + 1; end
    for (n, _) in regs; expected[n] = get(expected, n, 0); end
    net = Dict{Symbol, Int}()
    for (lhs, rhs, _, _) in rxns
        for s in lhs; s in met_set && (net[s] = get(net, s, 0) - 1); end
        for s in rhs; s in met_set && (net[s] = get(net, s, 0) + 1); end
    end
    for (name, coeff) in expected
        if coeff != 0
            haskey(net, name) ||
                error("$(coeff < 0 ? "Substrate" : "Product") $name does not appear in any reaction")
        end
    end
    for name in keys(net)
        haskey(expected, name) ||
            error("Metabolite $name not in species tuple")
    end
end
```

- [ ] **Step 4: Run the constructor test**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: new constructor test passes; many other tests still failing (still RED, as expected).

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "Implement EnzymeMechanism 2-arg constructor with canonicalization

Canonicalization: sort steps by (sorted reactants, sorted products, is_eq);
renumber kinetic_groups by first-occurrence order in the sorted list.
Two DSL orderings of the same semantic mechanism produce identical types.

KNOWN RED: downstream code still expects old signature."
```

### Task 2.3: Update `@enzyme_mechanism` macro to emit new type

**Files:**
- Modify: `src/dsl.jl`

- [ ] **Step 1: Write DSL macro test**

Add to `test/test_dsl.jl`:

```julia
@testset "@enzyme_mechanism (new grammar)" begin
    m = @enzyme_mechanism begin
        substrates: S[C]
        products:   P[C]
        regulators: I

        steps: begin
            [E, S]  ⇌    [ES]
            [ES, I] ⇌    [ESI]
            [ES]   <--> [EP]
            [EP]   ⇌    [E, P]
        end
    end

    @test EnzymeRates.substrates(m) == ((:S, ((:C, 1),)),)
    @test EnzymeRates.products(m) == ((:P, ((:C, 1),)),)
    @test EnzymeRates.regulators(m) == ((:I, ()),)
    @test EnzymeRates.n_steps(m) == 4
    @test count(EnzymeRates.equilibrium_steps(m)) == 3   # 3 RE, 1 SS

    # Grouping test
    m_grouped = @enzyme_mechanism begin
        substrates: S[C]
        products:   P[C]

        steps: begin
            ([E, S] ⇌ [ES], [EP, S] ⇌ [EPS])      # shared K group
            [ES]   <--> [EP]
            [EP]   ⇌   [E, P]
        end
    end

    @test EnzymeRates.kinetic_group(m_grouped, 1) ==
          EnzymeRates.kinetic_group(m_grouped, 2)   # both in same group after sort
end
```

- [ ] **Step 2: Run test — expect failure**

- [ ] **Step 3: Rewrite `@enzyme_mechanism` macro body**

Replace the current `@enzyme_mechanism` body in `src/dsl.jl` with the new grammar.

The full replacement:

```julia
"""
    @enzyme_mechanism begin
        substrates: S[C]
        products:   P[C]
        regulators: I
        steps: begin
            [E, S]  ⇌   [ES]
            [ES, I] ⇌   [ESI]
            [ES]   <--> [EP]
            [EP]   ⇌   [E, P]

            # Optional same-kinetics grouping:
            # ([E, S] ⇌ [ES], [ES, S] ⇌ [ESS])
        end
    end

Build a plain (non-allosteric) `EnzymeMechanism`. Rejects allosteric-only
constructs (`site(...)` / `::Tag` / `allosteric_regulators:` /
`catalytic_inhibitors:`) with clear errors.
"""
macro enzyme_mechanism(block)
    _reject_allosteric_syntax(block)
    mets_expr, rxns_expr = _parse_plain_mechanism_body(block)
    return esc(:(EnzymeMechanism($mets_expr, $rxns_expr)))
end

"""Reject allosteric-only DSL fragments in @enzyme_mechanism."""
function _reject_allosteric_syntax(block)
    for arg in block.args
        arg isa LineNumberNode && continue
        if arg isa Expr && arg.head == :call && arg.args[1] == :(:)
            label = arg.args[2]
            if label isa Expr && label.head == :call && label.args[1] == :site
                error("@enzyme_mechanism: `site(...)` belongs in @allosteric_mechanism, not here")
            end
            label in (:allosteric_regulators, :catalytic_inhibitors) &&
                error("@enzyme_mechanism: `$label:` is an allosteric-only field; " *
                      "use @allosteric_mechanism instead")
        end
    end
end
```

Implement `_parse_plain_mechanism_body(block)` to:
1. Parse `substrates:`, `products:`, `regulators:` blocks (no tags).
2. Parse `steps: begin ... end` where each line is one of:
   - A step tuple `[lhs...] ⇌ [rhs...]` or `[lhs...] <--> [rhs...]` (standalone — unique kinetic group).
   - A parenthesized group of step tuples — all share a kinetic group.
3. Assign kinetic group integers (arbitrary at parse time; constructor canonicalizes).
4. Return `(mets_expr, rxns_expr)` where:
   - `mets_expr` is `Expr(:tuple, subs_tuple, prods_tuple, regs_tuple)`.
   - `rxns_expr` is `Expr(:tuple, step_tuples...)`, each `(lhs_syms, rhs_syms, is_eq, group_num)`.

Rejects `::Tag` on any step:

```julia
function _parse_step_line(arg, next_group::Ref{Int})
    # Match: `(s1, s2, s3)` — a parenthesized tuple of steps sharing kinetics
    if arg isa Expr && arg.head == :tuple
        gnum = next_group[]; next_group[] += 1
        steps = Expr(:tuple)
        for e in arg.args
            push!(steps.args, _parse_single_step(e, gnum))
        end
        return steps.args  # multiple steps, same group
    end
    # Match: `[lhs] ⇌ [rhs]` or `[lhs] <--> [rhs]` — standalone
    gnum = next_group[]; next_group[] += 1
    return [_parse_single_step(arg, gnum)]
end

function _parse_single_step(expr, gnum::Int)
    expr isa Expr && expr.head == :call ||
        error("Expected [lhs] ⇌ [rhs] or [lhs] <--> [rhs], got: $expr")
    op = expr.args[1]
    is_eq = op == :⇌
    is_eq || op == :(<-->) ||
        error("Expected ⇌ or <-->, got operator: $op")
    lhs = _parse_step_side_symbols(expr.args[2])
    rhs = _parse_step_side_symbols(expr.args[3])
    Expr(:tuple, lhs, rhs, is_eq, gnum)
end
```

Reuse the existing `_parse_species_tuple_expr` and `_parse_step_side_symbols` helpers.

Delete the old `_parse_enzyme_mechanism` function body, the `_parse_allosteric_mechanism` helper (moved to new macro in Task 2.4), and the `_is_allosteric_label` detection.

- [ ] **Step 4: Run tests**

Expected: new DSL test passes. Most other tests still broken (they use old DSL grammar with `species:` block and `constraints:` — migrated in Task 2.5).

- [ ] **Step 5: Commit**

```bash
git add src/dsl.jl test/test_dsl.jl
git commit -m "Rewrite @enzyme_mechanism with new flat grammar

substrates: / products: / regulators: fields at top level.
steps: block accepts standalone step lines and parenthesized
step-groups (shared kinetics). No species/enzymes/constraints blocks.
Rejects allosteric-only syntax with clear error messages.

KNOWN RED: existing @enzyme_mechanism call sites in tests use old
grammar."
```

### Task 2.4: Write `@allosteric_mechanism` macro

**Files:**
- Modify: `src/dsl.jl`
- Modify: `src/EnzymeRates.jl` (export new macro)

- [ ] **Step 1: Write allosteric macro test**

Add to `test/test_dsl.jl`:

```julia
@testset "@allosteric_mechanism (new grammar)" begin
    m = @allosteric_mechanism begin
        substrates: S[C]
        products:   P[C]
        allosteric_regulators: I::OnlyT, A::OnlyT, R::NonequalRT
        catalytic_inhibitors:  J

        site(:catalytic, 2): begin
            steps: begin
                ([E, S] ⇌ [ES], [ES, S] ⇌ [ESS])   :: EqualRT
                [ES]    <--> [EP]                    :: OnlyR
                [EP]    ⇌   [E, P]                   :: EqualRT
                [ES, J] ⇌   [ESJ]                    :: OnlyR
            end
        end

        site(:regulatory, 2): begin
            ligands: A, I
        end
    end

    @test m isa EnzymeRates.AllostericEnzymeMechanism
    @test EnzymeRates.catalytic_multiplicity(m) == 2

    # Allosteric regulator tags
    allo_regs = EnzymeRates.allosteric_regulators(m)
    @test (:I, :OnlyT) in allo_regs
    @test (:A, :OnlyT) in allo_regs
    @test (:R, :NonequalRT) in allo_regs

    @test EnzymeRates.catalytic_inhibitors(m) == (:J,)
end
```

- [ ] **Step 2: Run test — expect failure**

- [ ] **Step 3: Implement the new macro**

Add to `src/dsl.jl`:

```julia
"""
    @allosteric_mechanism begin
        substrates: S[C]
        products:   P[C]
        allosteric_regulators: I::OnlyT, A::OnlyT, R::NonequalRT
        catalytic_inhibitors:  J

        site(:catalytic, N): begin
            steps: begin
                ...
            end
        end

        # optional
        site(:regulatory, N): begin
            ligands: A, I
        end
    end

Build an `AllostericEnzymeMechanism` (MWC). Required: one
`site(:catalytic, N):` block. `site(:regulatory, N):` blocks are
optional and declare competing-ligand reg sites; regulators not
listed in any explicit reg site get their own independent reg site
with multiplicity = N (the catalytic multiplicity).
"""
macro allosteric_mechanism(block)
    return esc(_parse_allosteric_mechanism_body(block))
end

function _parse_allosteric_mechanism_body(block)
    species_state = _AllostericSpecies()
    cat_n = nothing
    cat_steps_block = nothing
    reg_sites = Tuple{Tuple{Vararg{Symbol}}, Int}[]

    for arg in block.args
        arg isa LineNumberNode && continue
        _parse_allosteric_top_level!(
            arg, species_state, reg_sites,
            (n, b) -> (cat_n = n; cat_steps_block = b),
        )
    end

    cat_n === nothing && error("@allosteric_mechanism: site(:catalytic, N): is required")
    cat_steps_block === nothing &&
        error("@allosteric_mechanism: site(:catalytic, N): must contain a steps: block")

    rxns_expr, group_tags_expr = _parse_allosteric_steps(cat_steps_block)

    # Build CatalyticMech expression
    mets_expr = _build_allosteric_metabolites(species_state)
    cm_expr = :(EnzymeMechanism($mets_expr, $rxns_expr))

    # Build CatSites: (multiplicity, group_tags)
    cat_sites_expr = :(($cat_n, $group_tags_expr))

    # Build RegSites: distribute allosteric_regulators across sites
    reg_sites_expr = _build_reg_sites(species_state, reg_sites, cat_n)

    :(AllostericEnzymeMechanism($cm_expr, $cat_sites_expr, $reg_sites_expr))
end
```

Plus helpers (all in `src/dsl.jl`):
- `_AllostericSpecies` mutable struct holding parsed subs/prods/allosteric_regs/cat_inhibitors lists.
- `_parse_allosteric_top_level!` dispatches on each top-level label.
- `_parse_tagged_species_list` parses `X::OnlyT, Y::OnlyR` from the `allosteric_regulators:` line.
- `_parse_allosteric_steps` parses `steps: begin ... end` with required `::Tag` on each step or step-group.
- `_build_reg_sites` distributes ligands: explicit reg sites pull their ligands; unlisted ligands go to their own default site.

Reject `::EqualRT` on single-ligand reg site at parse time (or at construction if you prefer to defer to the `AllostericEnzymeMechanism` constructor — see Task 3.3).

Also export the new macro: in `src/EnzymeRates.jl`, add `@allosteric_mechanism` to the export list.

- [ ] **Step 4: Run tests**

Expected: `@allosteric_mechanism` test passes. `AllostericEnzymeMechanism` constructor not yet updated — test may pass via whatever stub / old constructor remains; some allosteric tests will still fail.

- [ ] **Step 5: Commit**

```bash
git add src/dsl.jl src/EnzymeRates.jl test/test_dsl.jl
git commit -m "Add @allosteric_mechanism macro

substrates: / products: / allosteric_regulators: / catalytic_inhibitors:
top-level fields. site(:catalytic, N): block required with required
::Tag on each step or step-group. site(:regulatory, N): optional for
competing ligands; unlisted allosteric regulators go to default
independent sites with multiplicity = N.

KNOWN RED: AllostericEnzymeMechanism constructor not yet updated."
```

### Task 2.5: Migrate plain-mechanism call sites in tests

**Files:**
- Modify: `test/test_dsl.jl`
- Modify: `test/test_types.jl`
- Modify: `test/test_enzyme_derivation.jl`
- Modify: `test/test_mechanism_enumeration.jl`
- Modify: `test/test_fitting.jl`
- Modify: `test/test_identify_rate_equation.jl`
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`
- Modify: `test/test_accessors.jl`
- Modify: `test/test_sym_poly.jl`

For each `@enzyme_mechanism` call site that uses the old grammar (with `species: begin ... end` + `steps: begin ... end` + `constraints: begin ... end`), rewrite to new grammar.

- [ ] **Step 1: Find all call sites**

```bash
grep -rn "@enzyme_mechanism" test/ --include="*.jl" | wc -l
grep -rln "@enzyme_mechanism" test/ --include="*.jl"
```

Record the count — should be ~30-50 sites. Work through them file by file.

- [ ] **Step 2: Per-file migration**

For each test file, replace each `@enzyme_mechanism` block as follows.

**Translation rules:**

1. `species: begin substrates: S[C]; products: P[C]; regulators: I; enzymes: E, ES[C] end` →
   ```
   substrates: S[C]
   products:   P[C]
   regulators: I
   ```
   (Drop `enzymes:` entirely — forms inferred from steps.)

2. `steps: begin [E,S]⇌[ES]; ...; end` stays as-is. Each step keeps its syntax.

3. `constraints: begin K2 = K1 end` → replace with parenthesized step grouping:
   Find the step whose binding K is `K2` and the step whose binding K is `K1` (usually determinable by step order). Wrap them in parentheses:
   ```
   ([E, S] ⇌ [ES], [EP, S] ⇌ [EPS])    # K1 = K2 now encoded as same kinetic group
   ```
   Multiple constraints `K2=K1, K3=K1, K4=K1` become one group: `([E,S]⇌[ES], [ES,S]⇌[ESS], [ESS,S]⇌[ESSS], [ESSS,S]⇌[ESSSS])`.

4. Same-kinetics grouping for SS steps (where old constraint was `k8f = k7f; k9f = k7f`): wrap the SS steps in parentheses. Note the new DSL groups both k_f and k_r for SS; the old DSL constrained only k_f. **If any old constraint applied to only k_f or only k_r but not both, the new DSL cannot express it.** Such tests need hand-review — if k_f and k_r were implicitly tied by Haldane, the new DSL produces the same equation. If they were truly independent (rare), the test must be dropped or marked.

- [ ] **Step 3: Run tests after each file**

After migrating each test file, run:
```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -20
```

Expected: per-file failures decrease. Goal: by end of Task 2.5, only allosteric-mechanism tests should still be failing.

- [ ] **Step 4: Commit progressively (one commit per file or logical batch)**

```bash
git add test/<file>.jl
git commit -m "Migrate <file> to new @enzyme_mechanism DSL"
```

Expected final state: plain `EnzymeMechanism` tests pass; allosteric tests still fail (handled in Phase 3).

### Task 2.6: Update rate-equation derivation for new `EnzymeMechanism` signature

**Files:**
- Modify: `src/rate_eq_derivation.jl`
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl`

The rate-equation derivation machinery (`_build_rate_body`, `_dependent_param_exprs`, etc.) currently indexes `Species[k]` and reads `param_constraints(m)` to fan out shared kinetic parameters. Update it to use the new accessors and kinetic-group semantics.

- [ ] **Step 1: Identify all magic-index reads to fix**

```bash
grep -n "Species\[" src/rate_eq_derivation.jl src/thermodynamic_constr_for_rate_eq_derivation.jl src/sym_poly_for_rate_eq_derivation.jl
grep -n "param_constraints\|.parameters\[" src/rate_eq_derivation.jl src/thermodynamic_constr_for_rate_eq_derivation.jl src/sym_poly_for_rate_eq_derivation.jl
```

Replace every `Species[1]` → `substrates(m)`, `Species[2]` → `products(m)`, `Species[3]` → `regulators(m)`, `Species[4]` → `enzyme_forms(m)`.

- [ ] **Step 2: Replace constraint fan-out with kinetic-group fan-out**

Where the old code does:

```julia
# For each constraint (target, coeff, factors), substitute target params
# with coeff * ∏ factor_i^exp_i in the polynomial.
for (target, coeff, factors) in param_constraints(m)
    ...
end
```

Replace with kinetic-group-based fan-out:

```julia
# For each kinetic group with 2+ members, substitute all non-representative
# step parameters with the representative's parameter.
for g in kinetic_groups(m)
    steps = steps_in_group(m, g)
    length(steps) == 1 && continue
    rep_idx = first(steps)
    for idx in Base.tail(steps)
        # Determine param name based on is_eq: K{idx} vs k{idx}f / k{idx}r
        if equilibrium_steps(m)[idx]
            # RE: K{idx} → K{rep_idx}
            push!(subs, Symbol("K$idx") => Symbol("K$rep_idx"))
        else
            push!(subs, Symbol("k$(idx)f") => Symbol("k$(rep_idx)f"))
            push!(subs, Symbol("k$(idx)r") => Symbol("k$(rep_idx)r"))
        end
    end
end
```

Apply this substitution to polynomials before the Haldane/Wegscheider closure runs.

- [ ] **Step 3: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: plain `EnzymeMechanism` rate-equation tests pass. Allosteric tests still failing.

- [ ] **Step 4: Commit**

```bash
git add src/rate_eq_derivation.jl src/thermodynamic_constr_for_rate_eq_derivation.jl
git commit -m "Rewrite rate-equation derivation for new EnzymeMechanism signature

Replace Species[k] magic-index reads with substrates/products/regulators/
enzyme_forms accessors. Replace ParamConstraints fan-out with
kinetic-group fan-out (steps with shared kinetic_group share their
K / k_f / k_r values).

Allosteric codepath still uses old CatSites/RegSites structure."
```

---

## Phase 3: New `AllostericEnzymeMechanism` type

### Task 3.1: Redefine type and accessors

**Files:**
- Modify: `src/types.jl`
- Modify: `test/test_types.jl`

- [ ] **Step 1: Write failing test for new allosteric type**

```julia
@testset "AllostericEnzymeMechanism struct (new design)" begin
    # Hand-build a minimal AllostericEnzymeMechanism
    cm = @enzyme_mechanism begin
        substrates: S[C]
        products:   P[C]
        steps: begin
            [E, S]  ⇌   [ES]
            [ES]   <--> [EP]
            [EP]   ⇌   [E, P]
        end
    end
    cat_sites = (2, ((1, :EqualRT), (2, :OnlyR), (3, :EqualRT)))
    reg_sites = ((((:I,), 2, ((:I, :OnlyT),)),),)
    m = EnzymeRates.AllostericEnzymeMechanism{typeof(cm), cat_sites, reg_sites}()

    @test EnzymeRates.catalytic_mechanism(m) === cm
    @test EnzymeRates.catalytic_multiplicity(m) == 2
    @test EnzymeRates.group_tag(m, 1) == :EqualRT
    @test EnzymeRates.group_tag(m, 2) == :OnlyR
    @test EnzymeRates.regulatory_sites(m) == reg_sites[1]  # first element only for this test
end
```

(Adjust the `reg_sites` literal shape if your chosen convention differs — the test is documentation of the chosen shape.)

- [ ] **Step 2: Run test — expect failure**

- [ ] **Step 3: Rewrite struct and accessors**

Replace the current `AllostericEnzymeMechanism` struct (lines 297-308 approximately in `src/types.jl`):

```julia
"""
    AllostericEnzymeMechanism{CatalyticMech, CatSites, RegSites}

Singleton type encoding a multi-subunit MWC allosteric enzyme.

- `CatalyticMech`: an `EnzymeMechanism` type (the single-subunit
  catalytic mechanism).
- `CatSites`: `(multiplicity::Int, group_tags::Tuple{Pair{Int,Symbol}...})`.
  Non-default TR tags only; missing entries default to `:NonequalRT`.
- `RegSites`: tuple of entries `((ligands::Tuple{Symbol,...}, multiplicity::Int,
  ligand_tags::Tuple{Pair{Symbol,Symbol}...}),)`. One entry per reg site.
"""
struct AllostericEnzymeMechanism{
    CatalyticMech, CatSites, RegSites,
} <: AbstractEnzymeMechanism end
```

Replace all allosteric accessors (lines 550-574 approximately):

```julia
catalytic_mechanism(::AllostericEnzymeMechanism{CM}) where {CM} = CM()
catalytic_multiplicity(::AllostericEnzymeMechanism{CM, CS}) where {CM, CS} = CS[1]

"""Tag for kinetic group `g`; default :NonequalRT if not stored."""
function group_tag(m::AllostericEnzymeMechanism, g::Int)
    CS = typeof(m).parameters[2]
    for (k, t) in CS[2]
        k == g && return t
    end
    :NonequalRT
end

step_tag(m::AllostericEnzymeMechanism, idx::Int) =
    group_tag(m, kinetic_group(catalytic_mechanism(m), idx))

# Shared structural accessors (delegate to catalytic mech)
substrates(m::AllostericEnzymeMechanism)      = substrates(catalytic_mechanism(m))
products(m::AllostericEnzymeMechanism)        = products(catalytic_mechanism(m))
reactions(m::AllostericEnzymeMechanism)       = reactions(catalytic_mechanism(m))
equilibrium_steps(m::AllostericEnzymeMechanism) = equilibrium_steps(catalytic_mechanism(m))
n_steps(m::AllostericEnzymeMechanism)         = n_steps(catalytic_mechanism(m))
enzyme_forms(m::AllostericEnzymeMechanism)    = enzyme_forms(catalytic_mechanism(m))
n_states(m::AllostericEnzymeMechanism)        = n_states(catalytic_mechanism(m))
kinetic_group(m::AllostericEnzymeMechanism, idx::Int) =
    kinetic_group(catalytic_mechanism(m), idx)
kinetic_groups(m::AllostericEnzymeMechanism) =
    kinetic_groups(catalytic_mechanism(m))
steps_in_group(m::AllostericEnzymeMechanism, g) =
    steps_in_group(catalytic_mechanism(m), g)

"""Return union of catalytic-mechanism regulators and reg-site ligands."""
function regulators(m::AllostericEnzymeMechanism)
    cat_regs = regulators(catalytic_mechanism(m))
    RS = typeof(m).parameters[3]
    names = Set(r[1] for r in cat_regs)
    extra = Symbol[]
    for entry in RS
        for lig in entry[1]
            lig ∉ names && (push!(names, lig); push!(extra, lig))
        end
    end
    (cat_regs..., Tuple((l, ()) for l in extra)...)
end

"""Return full metabolite list (catalytic + reg-site-only ligands)."""
function metabolites(m::AllostericEnzymeMechanism)
    cat_mets = metabolites(catalytic_mechanism(m))
    RS = typeof(m).parameters[3]
    names = Set(cat_mets)
    extra = Symbol[]
    for entry in RS
        for lig in entry[1]
            lig ∉ names && (push!(names, lig); push!(extra, lig))
        end
    end
    (cat_mets..., extra...)
end

"""Return the allosteric regulators as `((name, tag), ...)` pairs, one per ligand across all reg sites."""
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

"""Return dead-end-only regulator names (present in CatalyticMech.regulators but no reg site)."""
function catalytic_inhibitors(m::AllostericEnzymeMechanism)
    RS = typeof(m).parameters[3]
    reg_site_names = Set{Symbol}()
    for (ligands, _, _) in RS
        for lig in ligands; push!(reg_site_names, lig); end
    end
    cat_regs = regulators(catalytic_mechanism(m))
    Tuple(r[1] for r in cat_regs if r[1] ∉ reg_site_names)
end

"""Raw reg-site entries."""
regulatory_sites(::AllostericEnzymeMechanism{CM, CS, RS}) where {CM, CS, RS} = RS

regulatory_site_ligands(m::AllostericEnzymeMechanism, i::Int) =
    regulatory_sites(m)[i][1]

regulatory_site_multiplicity(m::AllostericEnzymeMechanism, i::Int) =
    regulatory_sites(m)[i][2]

function regulatory_ligand_tag(m::AllostericEnzymeMechanism, i::Int, lig::Symbol)
    for (k, t) in regulatory_sites(m)[i][3]
        k == lig && return t
    end
    :NonequalRT
end
```

- [ ] **Step 4: Run tests**

Expected: new type accessor tests pass; allosteric rate-equation tests still fail (constructor not updated, derivation not updated).

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "Redefine AllostericEnzymeMechanism with new CatSites/RegSites shape

CatSites = (multiplicity, group_tags). RegSites entries are
(ligands, multiplicity, ligand_tags). All access via named accessors
(catalytic_mechanism, catalytic_multiplicity, group_tag, step_tag,
allosteric_regulators, catalytic_inhibitors, regulatory_sites,
regulatory_site_{ligands,multiplicity}, regulatory_ligand_tag).
No magic indices in the implementation.

KNOWN RED: rate-equation derivation for allosteric uses old fields."
```

### Task 3.2: Implement `AllostericEnzymeMechanism` constructor

**Files:**
- Modify: `src/types.jl`
- Modify: `src/dsl.jl` (the `@allosteric_mechanism` expansion feeds into this constructor)

- [ ] **Step 1: Write constructor test**

```julia
@testset "AllostericEnzymeMechanism constructor" begin
    cm = @enzyme_mechanism begin
        substrates: S[C]
        products:   P[C]
        steps: begin
            [E, S]  ⇌   [ES]
            [ES]   <--> [EP]
            [EP]   ⇌   [E, P]
        end
    end
    cat_sites = (2, ((2, :OnlyR),))  # only non-default: group 2 OnlyR
    reg_sites = ((((:I,), 2, ((:I, :OnlyT),)),),)
    m = EnzymeRates.AllostericEnzymeMechanism(cm, cat_sites, reg_sites[1])

    @test EnzymeRates.group_tag(m, 1) == :NonequalRT  # default
    @test EnzymeRates.group_tag(m, 2) == :OnlyR
    @test EnzymeRates.step_tag(m, 2) == :OnlyR

    # Validation: EqualRT at single-ligand reg site → error
    bad_reg_sites = ((((:I,), 2, ((:I, :EqualRT),)),),)
    @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
        cm, cat_sites, bad_reg_sites[1],
    )
end
```

- [ ] **Step 2: Run test — expect failure**

- [ ] **Step 3: Implement the constructor**

Add to `src/types.jl`:

```julia
"""
    AllostericEnzymeMechanism(catalytic_mech, cat_sites, reg_sites)

Construct an `AllostericEnzymeMechanism`. Validates:
- Each reg site has at least one non-`:EqualRT` ligand (else reg_Q cancels
  identically in num/den ratio).
- Iso-group tags are not `:OnlyT` (R-inactive is a relabel).
- Group tags reference actual kinetic-group integers present in
  catalytic_mech.
- Reg-site ligands appear in the ligand_tags list iff their tag is
  not `:NonequalRT`.
"""
function AllostericEnzymeMechanism(
    cm::EnzymeMechanism,
    cat_sites::Tuple{Int, <:Tuple},
    reg_sites::Tuple,
)
    multiplicity, group_tags = cat_sites

    # Validate group tags
    valid_groups = Set(kinetic_groups(cm))
    for (g, tag) in group_tags
        g in valid_groups ||
            error("group_tag references non-existent kinetic_group $g")
        tag in (:OnlyR, :OnlyT, :EqualRT, :NonequalRT) ||
            error("Invalid group tag: $tag")
        # Iso groups cannot be :OnlyT
        step_idxs = steps_in_group(cm, g)
        for idx in step_idxs
            step = reactions(cm)[idx]
            is_binding = any(s == sym for sym in metabolites(cm) for s in step[1])
            if !is_binding && !step[3]  # SS iso step
                tag == :OnlyT &&
                    error("Iso step (kinetic_group $g) tagged :OnlyT is forbidden " *
                          "(R-inactive is a relabel)")
            end
        end
    end

    # Validate reg sites
    for entry in reg_sites
        ligands, mult, lig_tags = entry
        isempty(ligands) && error("Reg site must have at least one ligand")
        # Build tag map (defaults :NonequalRT)
        tag_map = Dict(lig_tags)
        all_equal = all(get(tag_map, l, :NonequalRT) == :EqualRT for l in ligands)
        all_equal &&
            error("Reg site with all `:EqualRT` ligands cancels identically; " *
                  "at least one ligand must have non-:EqualRT tag. Ligands: $ligands")
        # Sanity on tag vocabulary
        for (lig, tag) in lig_tags
            tag in (:OnlyR, :OnlyT, :NonequalRT, :EqualRT) ||
                error("Invalid reg-site tag: $tag for ligand $lig")
        end
    end

    # Sort group_tags by group number for canonical form
    sorted_tags = Tuple(sort(collect(group_tags); by=first))

    cat_sites_canonical = (multiplicity, sorted_tags)

    AllostericEnzymeMechanism{typeof(cm), cat_sites_canonical, reg_sites}()
end
```

- [ ] **Step 4: Run tests**

Expected: new constructor tests pass, allosteric DSL tests from Task 2.4 now produce valid types.

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "Implement AllostericEnzymeMechanism constructor with validation

Validates: iso group tagged :OnlyT forbidden; reg site with all
:EqualRT ligands forbidden (pure cancellation). Canonicalizes
cat_sites.group_tags by group number."
```

### Task 3.3: Rewrite allosteric rate-equation derivation

**Files:**
- Modify: `src/rate_eq_derivation.jl`
- Modify: `src/sym_poly_for_rate_eq_derivation.jl` (remove `_rs_*` helpers)

The current code has a parallel R/T-state derivation with many specialized helpers (`_is_tr_equiv_catalytic_param`, `_is_r_only_catalytic_param`, `_build_allosteric_rate_body`). Rewrite as: derive R-state via plain `EnzymeMechanism` machinery, then apply a single substitution pass over parameter symbols driven by `group_tag(m, g)` lookups.

- [ ] **Step 1: Write rate-equation tests to lock behavior**

Before refactoring, capture the expected rate-equation strings for a few allosteric mechanisms currently in `test/mechanism_definitions_for_test_enzyme_derivation.jl`. Specifically: `rate_mwc_dimer_oligo`, `rate_homodimer_noncomp_inh_oligo`, `rate_mwc_dimer_inh_oligo`. These are the existing reference closed forms.

Migrate each to new DSL:

```julia
# test/mechanism_definitions_for_test_enzyme_derivation.jl
mwc_dimer_new_dsl = @allosteric_mechanism begin
    substrates: S[C]
    products:   P[C]

    site(:catalytic, 2): begin
        steps: begin
            [E_c, S] ⇌   [E_S]    :: NonequalRT
            [E_c, P] ⇌   [E_P]    :: NonequalRT
            [E_S]   <--> [E_P]    :: NonequalRT
        end
    end
end

@testset "MWC dimer rate equation (baseline)" begin
    params = (K1=..., K2=..., k3f=..., k3r=..., K1_T=..., K2_T=..., k3f_T=..., k3r_T=..., L=..., Keq=..., E_total=...)
    concs = (S=1.0, P=0.5)
    @test rate_equation(mwc_dimer_new_dsl, concs, params) ≈
          rate_mwc_dimer_oligo(params, concs)
end
```

- [ ] **Step 2: Run tests — expect failure on allosteric rate equation**

The test expresses the intended semantics but the allosteric rate-equation generation is about to be rewritten.

- [ ] **Step 3: Delete old allosteric rate-equation code**

Delete from `src/rate_eq_derivation.jl`:
- `_is_tr_equiv_catalytic_K`
- `_is_tr_equiv_catalytic_param`
- `_is_r_only_catalytic_param`
- Old `_binding_K_symbols(::Type{<:AllostericEnzymeMechanism{...}})`
- Old `_dependent_param_exprs(::Type{<:AllostericEnzymeMechanism{...}})`
- Old `_allosteric_dep_assignments`
- Old `_allosteric_num_den_exprs`
- Old `_build_allosteric_rate_body`
- Old `@generated rate_equation(...,::AllostericEnzymeMechanism)`
- Old `rate_equation_string(::AllostericEnzymeMechanism, ::ReducedMode)`
- Old `structural_identifiability_deficit(::AllostericEnzymeMechanism)`

Delete from `src/sym_poly_for_rate_eq_derivation.jl`:
- `_rs_tr_equiv`, `_rs_r_only`, `_rs_t_only`
- `_count_allosteric_rate_monomials` (move logic into new helper)

- [ ] **Step 4: Implement new derivation path**

Add to `src/rate_eq_derivation.jl`:

```julia
# ═══════════════════════════════════════════════════════════════════
# AllostericEnzymeMechanism rate equation — new unified path
#
# MWC formula (summed over conformations c ∈ {R, T}):
#   num = CS[1] * sum_c( L_c * N_cat_c * Q_cat_c^(CS[1]-1)
#                        * prod_i( reg_Q_i_c^n_i ) )
#   den = sum_c( L_c * Q_cat_c^CS[1] * prod_i( reg_Q_i_c^n_i ) )
#   v = E_total * num / den
#
# R-state polynomials (N_cat_R, Q_cat_R, reg_Q_i_R) come from the
# plain EnzymeMechanism derivation on catalytic_mechanism(m). T-state
# polynomials are computed by applying a tag-driven symbol substitution
# to the R-state polynomials — never a parallel symbolic derivation.
# ═══════════════════════════════════════════════════════════════════

function _allosteric_rate_expr(m::Type{<:AllostericEnzymeMechanism})
    cm = m.parameters[1]
    multiplicity = m.parameters[2][1]
    group_tags = m.parameters[2][2]
    reg_sites = m.parameters[3]

    # R-state polynomials from plain-EnzymeMechanism derivation
    N_R, Q_R = _raw_rate_polys(cm)

    # T-state substitution map driven by group_tags
    tag_lookup = Dict(group_tags)
    t_subs = Dict{Symbol, Any}()
    for g in kinetic_groups(cm())
        tag = get(tag_lookup, g, :NonequalRT)
        group_members = steps_in_group(cm(), g)
        rep_idx = first(group_members)
        rep_step_is_eq = equilibrium_steps(cm())[rep_idx]
        # Determine symbol names
        if rep_step_is_eq
            sym_R = Symbol("K$rep_idx")
            sym_T = Symbol("K$(rep_idx)_T")
        else
            # SS — has two symbols (kf, kr); handle both
            for sfx in ("f", "r")
                sym_R = Symbol("k$rep_idx$sfx")
                sym_T = Symbol("k$rep_idx$sfx", "_T")
                _apply_tag_to_subs!(t_subs, tag, sym_R, sym_T)
            end
            continue
        end
        _apply_tag_to_subs!(t_subs, tag, sym_R, sym_T)
    end

    # Apply substitution to get T-state polynomials
    N_T = substitute_symbols(N_R, t_subs)
    Q_T = substitute_symbols(Q_R, t_subs)

    # Build reg-site partition-function polynomials (per site, per conformation)
    reg_Q_R = [_build_reg_Q(entry, :R) for entry in reg_sites]
    reg_Q_T = [_build_reg_Q(entry, :T) for entry in reg_sites]

    # Assemble num and den
    num = :(
        $multiplicity * (
            $N_R * $Q_R^($multiplicity - 1) * $(_reg_product(reg_Q_R, reg_sites, multiplicity))
            + L * $N_T * $Q_T^($multiplicity - 1) * $(_reg_product(reg_Q_T, reg_sites, multiplicity))
        )
    )
    den = :(
        $Q_R^$multiplicity * $(_reg_product(reg_Q_R, reg_sites, multiplicity))
        + L * $Q_T^$multiplicity * $(_reg_product(reg_Q_T, reg_sites, multiplicity))
    )
    :(E_total * ($num) / ($den))
end

function _apply_tag_to_subs!(t_subs, tag, sym_R, sym_T)
    if tag == :EqualRT
        t_subs[sym_R] = sym_R  # T uses same sym as R
    elseif tag == :OnlyR
        t_subs[sym_R] = 0     # T-state has this set to zero
    elseif tag == :OnlyT
        # R-state substitution is identity; T gets its own symbol
        # (handled differently — see below)
        t_subs[sym_R] = sym_T
    else  # :NonequalRT
        t_subs[sym_R] = sym_T
    end
end
```

Then `rate_equation(m::AllostericEnzymeMechanism, concs, params)` uses `_allosteric_rate_expr(typeof(m))` as its generated body.

Full implementation guidance:

1. **`_raw_rate_polys(cm)`**: returns the symbolic `(N_cat, Q_cat)` polynomial pair for the R-state by running the existing King-Altman / Cha derivation on `cm` (a plain `EnzymeMechanism`). This is the same machinery used by `rate_equation(::EnzymeMechanism, ...)` — factor it out so both paths call it.

2. **`_apply_tag_to_subs!`**: builds the T-state substitution dictionary. For `:EqualRT`, T symbol → R symbol (no rename). For `:OnlyR`, T symbol → literal `0` (zeroes T-state polynomial at that step's contribution). For `:OnlyT`, R symbol stays R, but the R-state polynomial needs the *R symbol zeroed* — track this separately (collect a second dict `r_zero_subs`). For `:NonequalRT`, T symbol is `Symbol("..._T")`.

3. **`_build_reg_Q(entry, conformation)`**: for each ligand in `entry[1]`, inspect its tag in `entry[3]` (default `:NonequalRT`). Emit `(1 + lig/K_lig_reg_i)` in R-state if ligand isn't `:OnlyT`; `(1 + lig/K_lig_T_reg_i)` in T-state if ligand isn't `:OnlyR`. `:EqualRT` ligand uses the R-state symbol in both conformations. Returns an `Expr`.

4. **`_reg_product(reg_Qs, reg_sites, multiplicity)`**: folds `prod_i (reg_Q_i ^ entry.multiplicity)` over all sites.

5. **Iso-step `:OnlyR` handling**: already falls out of (2) — both `k{idx}f_T` and `k{idx}r_T` substitute to `0`, zeroing the T-state polynomial where they appear.

Key simplification vs. current code: no separate `_dependent_param_exprs` for allosteric, no `_is_tr_equiv_catalytic_param`, no `_allosteric_num_den_exprs`. One function builds the full rate expression; substitution dictionaries encode all the TR logic.

- [ ] **Step 5: Implement `rate_equation_string` and `structural_identifiability_deficit` using the same path**

These two share the substitution machinery. Extract the common expression-building into a helper and call it from both.

- [ ] **Step 6: Run tests**

Expected: all allosteric-mechanism rate-equation tests pass, including the baseline MWC-dimer reference.

- [ ] **Step 7: Commit**

```bash
git add src/rate_eq_derivation.jl src/sym_poly_for_rate_eq_derivation.jl test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -m "Rewrite allosteric rate-equation derivation as substitution over R-state

T-state polynomials are produced by a tag-driven symbol substitution
over the R-state polynomials derived via the plain-EnzymeMechanism
path, not by a parallel symbolic derivation. Deletes _is_tr_equiv_*,
_is_r_only_*, _rs_tr_equiv/r_only/t_only, _count_allosteric_rate_monomials,
_allosteric_dep_assignments, and the 7-field-CatSites indexing throughout."
```

---

## Phase 4: Mechanism Enumeration

### Task 4.1: Update `MechanismSpec` / `AllostericMechanismSpec`

**Files:**
- Modify: `src/mechanism_enumeration.jl`
- Modify: `test/test_mechanism_enumeration.jl`

The enumeration's `StepSpec` already has `reactants`/`products`/`is_equilibrium`. Add `kinetic_group::Int` field. Delete `param_constraints` from `MechanismSpec`.

- [ ] **Step 1: Modify `StepSpec` struct**

```julia
struct StepSpec
    reactants::Vector{Symbol}
    products::Vector{Symbol}
    is_equilibrium::Bool
    kinetic_group::Int   # NEW
end

Base.:(==)(a::StepSpec, b::StepSpec) =
    a.reactants == b.reactants &&
    a.products == b.products &&
    a.is_equilibrium == b.is_equilibrium &&
    a.kinetic_group == b.kinetic_group

Base.hash(s::StepSpec, h::UInt) =
    hash(s.kinetic_group,
        hash(s.is_equilibrium,
            hash(s.products, hash(s.reactants, h))))
```

- [ ] **Step 2: Modify `MechanismSpec` struct**

```julia
struct MechanismSpec <: AbstractMechanismSpec
    reaction::Any
    steps::Vector{StepSpec}
    param_count::Int
end
```

Delete `param_constraints::Vector{ParamConstraint}`. Delete `ParamConstraint` type alias if no longer used.

- [ ] **Step 3: Modify `AllostericMechanismSpec` struct**

```julia
struct AllostericMechanismSpec <: AbstractMechanismSpec
    base::MechanismSpec
    catalytic_n::Int
    allosteric_reg_sites::Vector{Vector{Symbol}}
    allosteric_multiplicities::Vector{Int}
    group_tags::Dict{Int, Symbol}           # kinetic_group -> tag
    reg_ligand_tags::Dict{Symbol, Symbol}   # ligand -> tag (per reg site handled in sites structure if needed)
    param_count::Int
end
```

Delete `tr_equiv_metabolites`, `tr_equiv_cat_steps`, `r_only_metabolites`, `t_only_metabolites`, `r_only_cat_steps` fields.

- [ ] **Step 4: Update `EnzymeMechanism(spec::MechanismSpec)` constructor**

```julia
function EnzymeMechanism(spec::MechanismSpec)
    subs = substrates(spec.reaction)
    prods = products(spec.reaction)
    regs = regulator_roles(spec.reaction)  # ((name, role), ...)
    # Normalize regs to ((name, atoms=()),) shape
    regs_normalized = Tuple((r[1], ()) for r in regs)
    mets = (subs, prods, regs_normalized)

    rxns = Tuple(
        (Tuple(s.reactants), Tuple(s.products), s.is_equilibrium, s.kinetic_group)
        for s in spec.steps
    )
    EnzymeMechanism(mets, rxns)
end
```

- [ ] **Step 5: Update `AllostericEnzymeMechanism(spec::AllostericMechanismSpec)` constructor**

```julia
function AllostericEnzymeMechanism(spec::AllostericMechanismSpec)
    cm = EnzymeMechanism(spec.base)
    group_tags = Tuple((g, t) for (g, t) in spec.group_tags)
    cat_sites = (spec.catalytic_n, group_tags)

    # Build reg-site entries: group ligands by reg site
    reg_sites = Tuple(
        let ligs = Tuple(group),
            mult = spec.allosteric_multiplicities[i],
            lig_tags = Tuple(
                (l, spec.reg_ligand_tags[l])
                for l in group if haskey(spec.reg_ligand_tags, l))
            (ligs, mult, lig_tags)
        end
        for (i, group) in enumerate(spec.allosteric_reg_sites)
    )

    AllostericEnzymeMechanism(cm, cat_sites, reg_sites)
end
```

- [ ] **Step 6: Update all call sites in `mechanism_enumeration.jl` that construct `StepSpec` or `MechanismSpec`**

Every `StepSpec(reactants, products, is_eq)` becomes `StepSpec(reactants, products, is_eq, kinetic_group)`. Use fresh integers — each step gets a unique `kinetic_group` by default; enumeration moves that introduce shared kinetics set the same group on multiple steps.

- [ ] **Step 7: Run enumeration tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: some enumeration tests may fail because expansion moves haven't been updated yet. Expected target: at least basic `init_mechanisms` passes.

- [ ] **Step 8: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Update MechanismSpec / AllostericMechanismSpec for new type system

StepSpec gains kinetic_group::Int. MechanismSpec drops ParamConstraint
vector. AllostericMechanismSpec drops 5 TR/rr fields, gains group_tags
and reg_ligand_tags dicts. Constructors EnzymeMechanism(spec) and
AllostericEnzymeMechanism(spec) updated to produce new types."
```

### Task 4.2: Simplify expansion moves

**Files:**
- Modify: `src/mechanism_enumeration.jl`

Each expansion move currently has a `::MechanismSpec` and a `::AllostericMechanismSpec` method, joined by `_rewrap_allosteric`. With the new cleaner spec structure, these collapse.

- [ ] **Step 1: Unified `_expand_re_to_ss`**

```julia
function _expand_re_to_ss(spec::AbstractMechanismSpec)
    results = typeof(spec)[]
    for (i, step) in enumerate(_steps(spec))
        step.is_equilibrium || continue
        new_step = StepSpec(
            step.reactants, step.products, false, step.kinetic_group,
        )
        new_steps = copy(_steps(spec))
        new_steps[i] = new_step
        push!(results, _with_steps(spec, new_steps, _param_count(spec) + 1))
    end
    results
end

_steps(s::MechanismSpec) = s.steps
_steps(s::AllostericMechanismSpec) = s.base.steps
_param_count(s::MechanismSpec) = s.param_count
_param_count(s::AllostericMechanismSpec) = s.param_count

_with_steps(spec::MechanismSpec, new_steps, new_pc) =
    MechanismSpec(spec.reaction, new_steps, new_pc)
_with_steps(spec::AllostericMechanismSpec, new_steps, new_pc) =
    AllostericMechanismSpec(
        MechanismSpec(spec.base.reaction, new_steps, spec.base.param_count + 1),
        spec.catalytic_n,
        deepcopy(spec.allosteric_reg_sites),
        copy(spec.allosteric_multiplicities),
        copy(spec.group_tags),
        copy(spec.reg_ligand_tags),
        new_pc,
    )
```

- [ ] **Step 2: Replace `_expand_remove_constraint` with `_expand_split_kinetic_group`**

Old: remove one `ParamConstraint` entry, splitting one shared K.
New: split one kinetic group into two — pick a subset of steps in a multi-step group, give them a new group number.

```julia
function _expand_split_kinetic_group(spec::AbstractMechanismSpec)
    steps = _steps(spec)
    # Build group -> step_indices mapping
    groups = Dict{Int, Vector{Int}}()
    for (i, s) in enumerate(steps)
        push!(get!(groups, s.kinetic_group, Int[]), i)
    end
    results = typeof(spec)[]
    # Next group number to assign
    max_group = maximum(keys(groups); init=0)
    # For each group of size 2+, generate all ways to split off at least one step
    for (g, idxs) in groups
        length(idxs) >= 2 || continue
        for subset_mask in 1:(1 << length(idxs)) - 2
            # Skip empty subset and full subset (no-op)
            subset_indices = [idxs[j] for j in 1:length(idxs)
                              if (subset_mask >> (j-1)) & 1 == 1]
            new_group = max_group + 1
            new_steps = copy(steps)
            for i in subset_indices
                old = new_steps[i]
                new_steps[i] = StepSpec(old.reactants, old.products,
                                        old.is_equilibrium, new_group)
            end
            push!(results, _with_steps(spec, new_steps, _param_count(spec) + 1))
        end
    end
    results
end
```

- [ ] **Step 3: Delete `_valid_allosteric_differentiations`**

The old K-type / V-type hardcoded branches are replaced by uniform tag enumeration over groups.

- [ ] **Step 4: Rewrite `_expand_to_allosteric` as group-tag enumeration**

```julia
function _expand_to_allosteric(spec::MechanismSpec, @nospecialize(reaction::EnzymeReaction))
    cn = oligomeric_state(reaction)
    results = AllostericMechanismSpec[]

    # Identify kinetic groups and their types (RE binding / SS binding / SS iso)
    group_info = _kinetic_group_info(spec.steps)

    # Valid tag options per group
    for (g, info) in group_info
        if info.is_iso
            valid_tags = [:OnlyR, :EqualRT, :NonequalRT]  # no OnlyT for iso
        else
            valid_tags = [:OnlyR, :OnlyT, :EqualRT, :NonequalRT]
        end
        # Enumerate one differentiation at a time (beam-search friendly):
        # set group g to a non-default tag, leave all others default.
        for tag in valid_tags
            tag == :NonequalRT && continue  # default, no expansion
            gtags = Dict{Int,Symbol}(g => tag)
            push!(results, AllostericMechanismSpec(
                spec, cn,
                Vector{Symbol}[], Int[],
                gtags,
                Dict{Symbol,Symbol}(),
                spec.param_count + 1,
            ))
        end
    end
    results
end

function _expand_to_allosteric(::AllostericMechanismSpec, @nospecialize(::EnzymeReaction))
    AllostericMechanismSpec[]  # already allosteric
end
```

- [ ] **Step 5: Rewrite `_expand_add_allosteric_regulator` with tag enumeration**

```julia
function _expand_add_allosteric_regulator(
    spec::AllostericMechanismSpec,
    @nospecialize(reaction::EnzymeReaction),
)
    # Existing allosteric regs
    existing = Set(l for site in spec.allosteric_reg_sites for l in site)
    # Candidate new regulators from reaction
    new_regs = Symbol[]
    for (name, role) in regulator_roles(reaction)
        (role == :unknown || role == :allosteric) || continue
        name in existing && continue
        push!(new_regs, name)
    end
    sort!(new_regs)

    results = AllostericMechanismSpec[]
    for reg in new_regs
        n_sites = length(spec.allosteric_reg_sites)
        for tag in (:OnlyR, :OnlyT, :NonequalRT)
            # Add to new site
            new_sites = deepcopy(spec.allosteric_reg_sites)
            new_mults = copy(spec.allosteric_multiplicities)
            push!(new_sites, Symbol[reg])
            push!(new_mults, spec.catalytic_n)
            new_lig_tags = copy(spec.reg_ligand_tags)
            new_lig_tags[reg] = tag
            push!(results, AllostericMechanismSpec(
                spec.base, spec.catalytic_n,
                new_sites, new_mults,
                copy(spec.group_tags),
                new_lig_tags,
                spec.param_count + 1,
            ))

            # Add to each existing site (competing)
            for si in 1:n_sites
                new_sites2 = deepcopy(spec.allosteric_reg_sites)
                push!(new_sites2[si], reg)
                new_lig_tags2 = copy(spec.reg_ligand_tags)
                new_lig_tags2[reg] = tag
                push!(results, AllostericMechanismSpec(
                    spec.base, spec.catalytic_n,
                    new_sites2, copy(spec.allosteric_multiplicities),
                    copy(spec.group_tags),
                    new_lig_tags2,
                    spec.param_count + 1,
                ))
            end
        end
    end
    results
end
```

- [ ] **Step 6: Rewrite `_expand_remove_tr_equiv` as group-tag transition**

Replace with a `_expand_change_group_tag` move that changes one group's tag from its current value to `:NonequalRT` (or vice versa, depending on refinement strategy).

- [ ] **Step 7: Delete old variants**

Remove:
- `_tr_equiv_met_delta`
- Old `_rewrap_allosteric` (no longer needed, unified dispatch via `_with_steps`).
- Old `_valid_allosteric_differentiations`.
- Old `_expand_re_to_ss(spec::AllostericMechanismSpec)` method (the single-method version handles both now).
- Same for other expansion moves.

- [ ] **Step 8: Update canonicalization**

`_canonicalize!` must renumber kinetic_groups after step sort, matching the constructor's canonicalization:

```julia
function _canonicalize!(spec::MechanismSpec)
    # Sort steps
    sort!(spec.steps, by=_step_sort_key)
    # Renumber kinetic_groups by first-occurrence order
    old_to_new = Dict{Int, Int}()
    next = 1
    for (i, s) in enumerate(spec.steps)
        if !haskey(old_to_new, s.kinetic_group)
            old_to_new[s.kinetic_group] = next
            next += 1
        end
        spec.steps[i] = StepSpec(
            s.reactants, s.products, s.is_equilibrium,
            old_to_new[s.kinetic_group],
        )
    end
    old_to_new
end

function _canonicalize!(spec::AllostericMechanismSpec)
    group_remap = _canonicalize!(spec.base)
    # Remap group_tags
    new_gtags = Dict{Int, Symbol}()
    for (old_g, tag) in spec.group_tags
        if haskey(group_remap, old_g)
            new_gtags[group_remap[old_g]] = tag
        end
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

- [ ] **Step 9: Run enumeration tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: enumeration count invariants hold (bi-bi=11, ter-ter=283, etc.) — or new numbers emerge; if they do, record them as the new expected values and verify they're consistent (same mechanisms, different enumeration path). Denis reviews any divergence.

- [ ] **Step 10: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Unify expansion moves and delete K-type/V-type hardcoding

Each expansion move is now a single method parametric over spec type.
_valid_allosteric_differentiations deleted — iso-step/metabolite
tag enumeration emerges uniformly from per-group tag moves.
_rewrap_allosteric and _tr_equiv_met_delta deleted. Canonicalization
renumbers kinetic_groups by first-occurrence order after step sort.
Enumeration count invariants preserved / recorded."
```

---

## Phase 5: New PFK and HK Tests

### Task 5.1: PFK mechanism and analytical rate

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`
- Modify: `test/test_enzyme_derivation.jl`

- [ ] **Step 1: Write PFK mechanism definition in new DSL**

In `test/mechanism_definitions_for_test_enzyme_derivation.jl`, add:

```julia
pfk_mechanism = @allosteric_mechanism begin
    substrates: F6P[C6H12O6P], ATP[C10H16N5O13P3]
    products:   F16BP[C6H12O12P2], ADP[C10H15N5O10P2]
    allosteric_regulators:
        Pi::EqualRT,
        ATP::OnlyT,       # ATP also an allosteric regulator (OnlyT)
        ADP::OnlyR,
        Citrate::OnlyT,
        F26BP::NonequalRT

    site(:catalytic, 4): begin
        steps: begin
            # Random order bi-bi:  F6P can bind first or ATP can bind first.
            # Substrate bindings (two routes):
            ([E, F6P] ⇌ [E_F6P], [E_ATP, F6P] ⇌ [E_F6P_ATP])   :: OnlyR     # F6P K-type
            ([E, ATP] ⇌ [E_ATP], [E_F6P, ATP] ⇌ [E_F6P_ATP])   :: EqualRT
            # Isomerization step (SS):
            [E_F6P_ATP] <--> [E_F16BP_ADP]                        :: EqualRT
            # Product releases (two routes, random order):
            ([E_F16BP_ADP] ⇌ [E_ADP, F16BP], [E_F16BP] ⇌ [E, F16BP]) :: EqualRT
            ([E_F16BP_ADP] ⇌ [E_F16BP, ADP], [E_ADP] ⇌ [E, ADP])     :: EqualRT
        end
    end

    site(:regulatory, 4): begin
        ligands: Pi, ATP           # Pi and ATP compete at this site
    end
    # ADP, Citrate, F26BP each get their own independent reg site (mult=4)
end
```

- [ ] **Step 2: Write analytical rate function for PFK**

```julia
function pfk_rate_analytical(params, concs)
    # Named parameters expected (derived from new DSL's kinetic-group numbering)
    # ... explicit arithmetic here matching the MWC formula for this mechanism.
end
```

Derive by hand the N_cat_R, Q_cat_R, reg_Q_i_R (and T versions) polynomials following the MWC structure laid out in spec §10.2.1. Key points:

- `F6P :: OnlyR` on its binding group: T-state K_F6P absent → F6P terms appear only in R-state partition function.
- `ATP :: EqualRT` catalytic binding: K_ATP_T = K_ATP in catalytic Q.
- Iso step `:: EqualRT`: k_iso_T_f = k_iso_f, k_iso_T_r = k_iso_r.
- Reg site 1: `(1 + Pi/K_Pi)` in R (ATP is OnlyT so absent from R); `(1 + Pi/K_Pi + ATP/K_ATP_T_reg)` in T.
- ADP reg site: `(1 + ADP/K_ADP_R_reg)` in R; `1` in T (ADP is OnlyR).
- Citrate: `1` in R; `(1 + Citrate/K_Citrate_T_reg)` in T.
- F26BP: `(1 + F26BP/K_F26BP_R_reg)` in R; `(1 + F26BP/K_F26BP_T_reg)` in T.

- [ ] **Step 3: Write test**

In `test/test_enzyme_derivation.jl`:

```julia
@testset "PFK rate equation matches analytical form" begin
    concs_test = (F6P=0.1, ATP=1.0, F16BP=0.001, ADP=0.1, Pi=1.0,
                  Citrate=0.01, F26BP=0.01)
    params_test = (; K_F6P=0.1, K_ATP=0.5, ..., Keq=1000.0, E_total=1.0, L=1.0)
    @test rate_equation(pfk_mechanism, concs_test, params_test) ≈
          pfk_rate_analytical(params_test, concs_test)

    # Test: at [F6P] → 0, rate → 0 (F6P OnlyR means no T-path)
    low_f6p = merge(concs_test, (F6P=1e-10,))
    @test rate_equation(pfk_mechanism, low_f6p, params_test) < 1e-6

    # Test: Pi shifts R/T ratio when ATP is present at reg site
    # (non-trivial because ATP is OnlyT, so reg_Q_R ≠ reg_Q_T)
    concs_noPi = merge(concs_test, (Pi=0.0,))
    concs_Pi   = merge(concs_test, (Pi=10.0,))
    r_noPi = rate_equation(pfk_mechanism, concs_noPi, params_test)
    r_Pi   = rate_equation(pfk_mechanism, concs_Pi, params_test)
    @test r_noPi != r_Pi  # Pi has an effect via asymmetric coupling with ATP
end
```

- [ ] **Step 4: Run test**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: PFK test passes.

- [ ] **Step 5: Commit**

```bash
git add test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_enzyme_derivation.jl
git commit -m "Add PFK hand-verified rate-equation test

Tests: F6P OnlyR (K-type), ATP both substrate and allosteric regulator
with different tags per context, Pi EqualRT at reg site with OnlyT
co-ligand (non-cancelling), ADP OnlyR own site, Citrate OnlyT,
F26BP NonequalRT. Validates the full allosteric rate-equation
machinery against a hand-derived analytical form."
```

### Task 5.2: HK mechanism and analytical rate

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`
- Modify: `test/test_enzyme_derivation.jl`

- [ ] **Step 1: Write HK mechanism**

```julia
hk_mechanism = @allosteric_mechanism begin
    substrates: Glucose[C6H12O6], ATP[C10H16N5O13P3]
    products:   G6P[C6H12O9P], ADP[C10H15N5O10P2]
    allosteric_regulators: G6P::OnlyT, Pi::EqualRT
    catalytic_inhibitors:  G6P        # same symbol, third role: dead-end inhibitor

    site(:catalytic, 2): begin
        steps: begin
            # Random-order bi-bi substrate binding
            ([E, Glucose] ⇌ [E_Glc], [E_ATP, Glucose] ⇌ [E_Glc_ATP])   :: EqualRT
            ([E, ATP] ⇌ [E_ATP], [E_Glc, ATP] ⇌ [E_Glc_ATP])           :: EqualRT
            # Iso step
            [E_Glc_ATP] <--> [E_G6P_ADP]                                :: EqualRT
            # Random-order product release
            ([E_G6P_ADP] ⇌ [E_ADP, G6P], [E_G6P] ⇌ [E, G6P])          :: EqualRT
            ([E_G6P_ADP] ⇌ [E_G6P, ADP], [E_ADP] ⇌ [E, ADP])          :: EqualRT
            # Dead-end G6P binding (catalytic inhibitor role)
            [E_ATP, G6P] ⇌ [E_ATP_G6P]                                  :: EqualRT
            [E_ADP, G6P] ⇌ [E_ADP_G6P]                                  :: EqualRT
        end
    end

    site(:regulatory, 2): begin
        ligands: G6P, Pi            # G6P and Pi compete; Pi EqualRT allowed because G6P is OnlyT
    end
end
```

- [ ] **Step 2: Write analytical rate function `hk_rate_analytical`**

Derive per spec §10.2.2. The three G6P roles:
1. Product in catalytic release steps.
2. Catalytic inhibitor: two dead-end steps, each with its own K.
3. Allosteric inhibitor: reg site 1, OnlyT.

- [ ] **Step 3: Write test**

```julia
@testset "HK rate equation matches analytical form" begin
    concs_test = (Glucose=1.0, ATP=1.0, G6P=0.01, ADP=0.1, Pi=1.0)
    params_test = (...; Keq=..., E_total=1.0, L=1.0)
    @test rate_equation(hk_mechanism, concs_test, params_test) ≈
          hk_rate_analytical(params_test, concs_test)

    # Test: G6P inhibits (three coupled mechanisms)
    high_g6p = merge(concs_test, (G6P=10.0,))
    low_g6p  = merge(concs_test, (G6P=0.001,))
    @test rate_equation(hk_mechanism, high_g6p, params_test) <
          rate_equation(hk_mechanism, low_g6p, params_test)

    # Substrate-type vs inhibitor-type grouping: G6P release steps
    # and G6P dead-end steps MUST be in separate kinetic groups.
    # Confirm construction didn't merge them (would error).
end
```

- [ ] **Step 4: Run test**

- [ ] **Step 5: Commit**

```bash
git add test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_enzyme_derivation.jl
git commit -m "Add HK hand-verified rate-equation test

G6P appears in three independent roles: product (release steps),
catalytic_inhibitor (dead-end steps), allosteric inhibitor (OnlyT
reg-site ligand). Pi EqualRT allowed at reg site because G6P is
OnlyT co-ligand. Substrate-type and inhibitor-type G6P bindings
stay in separate kinetic groups."
```

### Task 5.3: Narrow feature tests

**Files:**
- Modify: `test/test_enzyme_derivation.jl`

- [ ] **Step 1: Write single-feature edge-case tests**

```julia
@testset "Single-feature allosteric edge cases" begin
    # OnlyT substrate — rate → 0 as K_T → ∞
    onlyT_sub = @allosteric_mechanism begin
        substrates: S[C]
        products:   P[C]
        site(:catalytic, 2): begin
            steps: begin
                [E, S] ⇌   [ES]    :: OnlyT
                [ES]   <--> [EP]    :: EqualRT
                [EP]   ⇌   [E, P]   :: EqualRT
            end
        end
    end
    # ... rate at typical concs should be low / proportional to T-path only

    # Iso OnlyR alone (V-type) — verify T-state numerator is zero
    vtype = @allosteric_mechanism begin
        substrates: S[C]; products: P[C]
        site(:catalytic, 2): begin
            steps: begin
                [E, S] ⇌   [ES]    :: EqualRT
                [ES]   <--> [EP]    :: OnlyR       # V-type
                [EP]   ⇌   [E, P]   :: EqualRT
            end
        end
    end
    # ... verify rate → 0 as L → ∞ (fully T-state)

    # EqualRT single-ligand reg site → construction error
    @test_throws ErrorException @allosteric_mechanism begin
        substrates: S[C]; products: P[C]
        allosteric_regulators: I::EqualRT     # alone at its site — cancels identically
        site(:catalytic, 2): begin
            steps: begin
                [E, S] ⇌   [ES]    :: EqualRT
                [ES]   <--> [EP]    :: EqualRT
                [EP]   ⇌   [E, P]   :: EqualRT
            end
        end
    end

    # HK substrate-type + inhibitor-type G6P in one group → error
    @test_throws ErrorException @allosteric_mechanism begin
        substrates: S[C]
        products:   P[C]
        catalytic_inhibitors: P        # P as catalytic inhibitor too
        site(:catalytic, 2): begin
            steps: begin
                [E, S] ⇌ [ES]                                      :: EqualRT
                ([E_P, S] ⇌ [E_P_S], [E_S, P] ⇌ [E_S_P])           :: EqualRT
                # ^ ERROR: First binds S to a form without S; second binds P to E_S which already has S.
                # (Adjust example — the point is to construct a group mixing binding
                #  types of the same metabolite.)
                [ES] <--> [EP]                                      :: EqualRT
                [EP] ⇌ [E, P]                                       :: EqualRT
            end
        end
    end
end
```

- [ ] **Step 2: Run tests**

- [ ] **Step 3: Commit**

```bash
git add test/test_enzyme_derivation.jl
git commit -m "Add narrow allosteric feature edge-case tests

OnlyT substrate, V-type iso alone, single-ligand EqualRT reg site
rejection, substrate-type + inhibitor-type same-metabolite grouping
rejection."
```

---

## Phase 6: Final Cleanup

### Task 6.1: Magic-index audit

**Files:** various `src/`

- [ ] **Step 1: Grep for magic-index access**

```bash
grep -rn "Species\[\|CS\[\|RS\[\|\.parameters\[" src/ --include="*.jl" | grep -v "^.*#"
```

Expected: zero hits. Any remaining hits → fix by replacing with named accessor.

- [ ] **Step 2: Grep for stale function names**

```bash
grep -rn "_is_tr_equiv_catalytic\|_is_r_only_catalytic\|_rs_tr_equiv\|_rs_r_only\|_rs_t_only\|_rewrap_allosteric\|_tr_equiv_met_delta\|_valid_allosteric_differentiations\|param_constraints" src/ --include="*.jl"
```

Expected: zero hits (all deleted).

- [ ] **Step 3: Delete `src/old_mechanism_enumeration.jl` and `src/old_beam_enumeration.jl` if present**

These are preserved legacy files. With the refactor complete they should be removed (they reference the old type signatures and will now fail to load).

```bash
ls src/old_*.jl 2>/dev/null
# If present:
git rm src/old_mechanism_enumeration.jl src/old_beam_enumeration.jl
# Update src/EnzymeRates.jl to remove any include("old_...") lines
```

- [ ] **Step 4: Run full suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all green.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "Delete legacy old_*.jl files and stale references

Audit passed: no magic-index type-parameter access anywhere in src/,
no references to deleted helpers (_is_tr_equiv_*, _rs_*, etc.)."
```

### Task 6.2: Update CLAUDE.md

**Files:**
- Modify: `.claude/CLAUDE.md`

Update the architecture section to reflect the new type signatures, DSL shape, and testing approach. Key sections to touch:

- "Key Architecture Decisions" — replace old `EnzymeMechanism` / `AllostericEnzymeMechanism` descriptions with new ones.
- "Canonical Step Form" — add the kinetic-group renumbering canonical rule.
- "Regulator representation" — rewrite to new `allosteric_regulators` / `catalytic_inhibitors` separation.
- "Dead-end SS/RE propagation" and "Dead-end parameter equivalence constraints" — replace or delete.
- Remove "Catalytic topology constraints" references to deleted helpers.
- "Mechanism enumeration building blocks" — update move descriptions.
- "Vmax Normalization" — verify still accurate post-refactor.
- "Known Issues" — update `rate_equation` discussion if behavior changed.

- [ ] **Step 1: Write the updates inline**

Read current `.claude/CLAUDE.md`, replace sections as noted above. Keep the rules (Rule #1, foundational rules, etc.) intact.

- [ ] **Step 2: Commit**

```bash
git add .claude/CLAUDE.md
git commit -m "Update CLAUDE.md for new mechanism type signatures and DSL"
```

### Task 6.3: Update SPEC.md

**Files:**
- Modify: `SPEC.md`

SPEC.md currently documents `compile_mechanism` as exported; the new design keeps it as the `EnzymeMechanism(spec)` / `AllostericEnzymeMechanism(spec)` constructors (internal). Also add `@allosteric_mechanism` to exports.

- [ ] **Step 1: Update exported-API tables and core-workflow examples**

- [ ] **Step 2: Add `@allosteric_mechanism` to the macro table**

- [ ] **Step 3: Update the "Complete Exported Symbol List" code block**

```julia
export @enzyme_mechanism, @allosteric_mechanism
```

- [ ] **Step 4: Commit**

```bash
git add SPEC.md
git commit -m "Update SPEC.md for new macros and refactored type system"
```

### Task 6.4: Line-count reduction check

**Files:** none (verification only)

- [ ] **Step 1: Compare against baseline**

```bash
wc -l src/types.jl src/dsl.jl src/rate_eq_derivation.jl src/mechanism_enumeration.jl src/sym_poly_for_rate_eq_derivation.jl src/thermodynamic_constr_for_rate_eq_derivation.jl
cat REFACTOR_BASELINE.txt
```

Expected reductions vs. spec §9.1 targets:
- `src/types.jl`: ≥30% reduction.
- `src/rate_eq_derivation.jl`: ≥25% reduction.
- `src/mechanism_enumeration.jl`: ≥30% reduction.
- `src/sym_poly_for_rate_eq_derivation.jl`: `_rs_*` helpers deleted.
- `src/dsl.jl`: small net change acceptable (two macros, but less constraint-RHS machinery).

Record actual numbers; investigate deviations.

- [ ] **Step 2: Clean up baseline file**

```bash
rm REFACTOR_BASELINE.txt
```

- [ ] **Step 3: Final full-suite check**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all green. If not, stop and address issues before PR.

---

## Completion

Final git diff against `origin/main` should:
- Delete substantial code across the touched src/ files (verify against spec §9 and §9.1 targets).
- Add the new PFK and HK test fixtures.
- Leave the public API (`rate_equation`, `rate_equation_string`, `parameters`, `identify_rate_equation`, `fit_rate_equation`, `FittingProblem`, `IdentifyRateEquationProblem`, `IdentifyRateEquationResults`) functionally unchanged; only internals reshaped.

Push the branch and open a PR against `main` with the spec linked as the design reference.

```bash
git push -u origin allosteric-refactor-spec
gh pr create --title "Mechanism types refactor" --body "$(cat <<'EOF'
## Summary
- Refactor EnzymeMechanism and AllostericEnzymeMechanism with simpler, non-redundant type parameters
- Add per-step TR-mode tagging (OnlyR/OnlyT/EqualRT/NonequalRT) with @allosteric_mechanism macro
- Unify step-grouping into a single DSL construct (parenthesized tuple) that encodes shared kinetics + shared TR mode
- Eliminate magic-index access; strict accessor-only internal interface
- Delete dead code: graph() accessor, RegulatorRole type hierarchy, complex ParamConstraint monomial machinery, parallel R/T-state derivation paths

## Design reference
docs/superpowers/specs/2026-04-23-mechanism-types-refactor-design.md

## Test plan
- [ ] Full test suite passes (Aqua + JET + all unit/integration tests)
- [ ] PFK hand-verified rate equation matches analytical form
- [ ] HK hand-verified rate equation matches analytical form (G6P in three roles)
- [ ] Enumeration count invariants preserved or recorded
- [ ] Line-count reduction meets spec §9.1 targets

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

---

## Spec-coverage self-check

Before executing, verify this plan covers every spec section:

- **§1 Primary goal (code reduction):** Task 6.4 validates against §9.1 targets.
- **§2 Motivation (taxonomy, DSL, tests, streamlined derivation):** Phases 2-5 implement all four.
- **§3 Current state:** informational only.
- **§4.1 EnzymeMechanism type:** Task 2.1 + 2.2.
- **§4.2 AllostericEnzymeMechanism type:** Task 3.1 + 3.2.
- **§5 Accessors:** Task 2.1, 3.1.
- **§6 DSL:** Task 2.3 + 2.4.
- **§7 Tag semantics:** Task 3.3 (rate-eq derivation) embodies.
- **§8 Error cases:** Task 2.3, 3.2 (construction-time), 2.4 (DSL-level).
- **§9 Dead code:** Phase 1 (Task 1.1, 1.2, 1.3) + Phase 3 + 4 (_rs_*, _rewrap_allosteric, etc.).
- **§9.1 Reduction targets:** Task 6.4.
- **§10.1 DSL tests:** Tasks 2.3, 2.4 (new DSL tests); Task 2.5 (migrated existing).
- **§10.2 PFK / HK / edge cases:** Tasks 5.1, 5.2, 5.3.
- **§10.3 Enumeration invariants:** Task 4.2 step 9.
- **§10.4 kcat invariants:** covered by existing suite + Task 5 kcat-related assertions.
- **§10.5 Aqua/JET:** covered by existing suite run in every task.
- **§11 Migration notes:** Task 2.5 (plain) + embedded in Phase 3 for allosteric.
- **§12 Out-of-scope:** no plan tasks — confirmed not touching listed APIs.
- **§13 Sequence:** this plan *is* step 2.
