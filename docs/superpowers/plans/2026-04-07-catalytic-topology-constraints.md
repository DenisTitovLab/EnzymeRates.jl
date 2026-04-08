# Catalytic Topology Constraints Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the path combining bug and add biochemical constraints to `_catalytic_topologies` per the design spec at `docs/superpowers/specs/2026-04-06-catalytic-topology-constraints-design.md`.

**Architecture:** Modify `_catalytic_topologies` in `src/mechanism_enumeration.jl` to: (1) replace arbitrary subset combining with weak-ordering enumeration, (2) allow empty-residual ping-pong, (3) cap isomerization at 3×3, (4) require substrate participation for isomerization, (5) use product-only iso form names, (6) support multi-product release. Fix `_bound_metabolites_at_forms` for partial isomerization. Fix DSL multi-atom parsing. Tests use two real enzymes — pyruvate carboxylase and pyruvate dehydrogenase.

**Review findings addressed:**
- DSL fix uses per-arg atom parsing (not string join) to avoid `CO2,H` → `Co2H` bug
- `_bound_metabolites_at_forms` fixed for partial iso (ping-pong)
- C8 product-only naming verified safe (no collisions within a topology)
- C7 test strengthened to check all Estar-prefixed iso sources, not just bare `Estar`
- `_steps_for_weak_ordering` processes substrate and product orderings independently

**Tech Stack:** Julia, EnzymeRates.jl internals

**Spec:** `docs/superpowers/specs/2026-04-06-catalytic-topology-constraints-design.md`

---

### Task 1: Fix DSL multi-atom parsing

The `@enzyme_reaction` DSL currently drops all atoms after the first in
bracket syntax (`A[C,N]` → only captures `C`). Fix this so all atoms are
captured. This is needed for pyruvate carboxylase and pyruvate dehydrogenase
test reactions.

**Files:**
- Modify: `src/dsl.jl:23-34`
- Test: `test/test_dsl.jl`

- [ ] **Step 1: Write failing test for multi-atom DSL**

Add to the `@enzyme_reaction` testset in `test/test_dsl.jl`:

```julia
@testset "multi-atom metabolites" begin
    rxn = @enzyme_reaction begin
        substrates: A[C2H3], B[N,P]
        products: P[C2,N], Q[H3,P]
    end
    subs = EnzymeRates.substrates(rxn)
    @test subs[1] == (:A, ((:C, 2), (:H, 3)))
    @test subs[2] == (:B, ((:N, 1), (:P, 1)))
    prods = EnzymeRates.products(rxn)
    @test prods[1] == (:P, ((:C, 2), (:N, 1)))
    @test prods[2] == (:Q, ((:H, 3), (:P, 1)))
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -20`
Expected: FAIL — `subs[2]` will be `(:B, ((:N, 1),))` (missing P)

- [ ] **Step 3: Fix `_parse_species_tuple_expr` in `src/dsl.jl`**

Replace lines 23-34 with:

```julia
function _parse_species_tuple_expr(expr)
    if expr isa Symbol
        atoms = Expr(:tuple)
        return Expr(:tuple, QuoteNode(expr), atoms)
    elseif expr isa Expr && expr.head == :ref
        name = expr.args[1]
        # Parse each ref arg individually and merge
        # (handles A[C,N] where Julia parses as
        # ref with args [:A, :C, :N]).
        # Each arg is parsed as a chemical formula
        # separately to avoid ambiguity (e.g.,
        # A[CO2,H] must not become "CO2H" → cobalt).
        atoms = Expr(:tuple)
        for arg in expr.args[2:end]
            parsed = _parse_chemical_formula(string(arg))
            for atom in parsed.args
                push!(atoms.args, atom)
            end
        end
        return Expr(:tuple, QuoteNode(name), atoms)
    else
        error("Cannot parse species definition: $expr")
    end
end
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass including new multi-atom test

- [ ] **Step 5: Commit**

```bash
git add src/dsl.jl test/test_dsl.jl
git commit -m "Fix DSL to capture all atoms in multi-atom metabolites

A[C,N] was only capturing the first atom (C). Now parses each
bracket arg individually to avoid ambiguity (e.g., A[CO2,H]
must not be parsed as cobalt)."
```

---

### Task 2: Define test reactions for pyruvate carboxylase and pyruvate dehydrogenase

Create shared test reaction definitions that will be used by all subsequent
tasks. These use multi-atom metabolites (requires Task 1).

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (add after existing reaction
  definitions, around line 78)

- [ ] **Step 1: Add reaction definitions**

Add after the existing `ter_ter_rxn` and `ter_bi_rxn` definitions:

```julia
# Pyruvate carboxylase: Pyr + HCO3 + ATP = OAA + ADP + Pi
# Mechanism: ATP+HCO3 → ADP+Pi+CO2_residual,
#            then Pyr+CO2 → OAA
const pyruvate_carboxylase_rxn = @enzyme_reaction begin
    substrates: Pyr[C3H3O3], HCO3[HCO3], ATP[C10H16N5O13P3]
    products: OAA[C4H3O5], ADP[C10H15N5O10P2], Pi[H2PO4]
end

# Pyruvate dehydrogenase: Pyr + NAD + CoA = AcCoA + NADH + CO2
# Mechanism: Pyr → CO2+residual, CoA+residual → AcCoA+residual,
#            NAD+residual → NADH
const pyruvate_dehydrogenase_rxn = @enzyme_reaction begin
    substrates: Pyr[C3H3O3], NAD[C21H28N7O14P2], CoA[C21H36N7O16P3S]
    products: AcCoA[C23H38N7O17P3S], NADH[C21H29N7O14P2], CO2[CO2]
end
```

- [ ] **Step 2: Add atom balance sanity tests**

```julia
@testset "test reaction atom balance" begin
    for rxn in [pyruvate_carboxylase_rxn,
                pyruvate_dehydrogenase_rxn]
        sub_atoms = Dict{Symbol,Int}()
        for (_, atoms) in EnzymeRates.substrates(rxn)
            for (a, c) in atoms
                sub_atoms[a] = get(sub_atoms, a, 0) + c
            end
        end
        prod_atoms = Dict{Symbol,Int}()
        for (_, atoms) in EnzymeRates.products(rxn)
            for (a, c) in atoms
                prod_atoms[a] = get(prod_atoms, a, 0) + c
            end
        end
        @test sub_atoms == prod_atoms
    end
end
```

- [ ] **Step 3: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All pass

- [ ] **Step 4: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Add pyruvate carboxylase and dehydrogenase test reactions"
```

---

### Task 3: Implement weak-ordering path combining (C1 — bug fix)

Replace the arbitrary 2^n subset combining with weak-ordering enumeration.
This is the highest-priority fix.

**Files:**
- Modify: `src/mechanism_enumeration.jl` (lines 509-628 in
  `_catalytic_topologies`)
- Test: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing test for weak-ordering combining**

Add a new testset in the mechanism enumeration tests:

```julia
@testset "weak-ordering combining" begin
    # For bi-bi: 2 subs × 2 prods, all sequential
    # Weak orderings of 2 items = 3 (AB, BA, random)
    # Total: 3 × 3 = 9 topologies (unchanged from before)
    bi_bi_rxn_test = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    topos = EnzymeRates._catalytic_topologies(
        bi_bi_rxn_test)
    @test length(topos) == 9

    # For ter-ter: 3 subs × 3 prods, all sequential
    # Weak orderings of 3 items = 13
    # Total: 13 × 13 = 169 topologies
    topos_tt = EnzymeRates._catalytic_topologies(
        ter_ter_rxn)
    @test length(topos_tt) == 169
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep "weak-ordering"`
Expected: FAIL — ter-ter produces 3969 instead of 169

- [ ] **Step 3: Implement weak-ordering combining**

Replace lines 509-628 of `_catalytic_topologies` in
`src/mechanism_enumeration.jl`. The new logic:

1. Deduplicate paths by step content (unchanged)
2. Group paths by isomerization pattern (replaces pp/seq partition)
3. For each group, identify the metabolites that vary in ordering within
   each half-reaction
4. Enumerate weak orderings of those metabolites
5. Build topology step sets from each weak ordering

```julia
    # --- Deduplicate paths by step content (as sets) ---
    StepKey = Tuple{Vector{Symbol}, Vector{Symbol}}
    _step_key(s::StepSpec) = (
        sort(s.reactants), sort(s.products)
    )::StepKey
    unique_paths = Vector{Vector{StepSpec}}()
    seen_path_keys = Set{Set{StepKey}}()
    for path in all_paths
        key = Set(_step_key(s) for s in path)
        key ∈ seen_path_keys && continue
        push!(seen_path_keys, key)
        push!(unique_paths, path)
    end

    # --- Group paths by isomerization pattern ---
    met_names = Set{Symbol}(sub_names)
    union!(met_names, prod_names)

    function _iso_pattern(path)
        Set(
            (_step_key(s) for s in path
             if length(s.reactants) == 1 &&
                length(s.products) == 1 &&
                s.reactants[1] ∉ met_names &&
                s.products[1] ∉ met_names)
        )
    end

    iso_groups = Dict{
        Set{StepKey}, Vector{Vector{StepSpec}}
    }()
    for path in unique_paths
        pat = _iso_pattern(path)
        paths_vec = get!(iso_groups, pat,
            Vector{Vector{StepSpec}}())
        push!(paths_vec, path)
    end

    # --- Enumerate weak orderings within each group ---
    # A weak ordering partitions n items into ordered
    # levels. For each group, identify which metabolites
    # vary in ordering, then enumerate all weak orderings
    # and build the corresponding topology step sets.

    function _weak_orderings(items::Vector{T}) where T
        n = length(items)
        n == 0 && return [Vector{Vector{T}}()]
        n == 1 && return [[items]]
        result = Vector{Vector{Vector{T}}}()
        # Enumerate all partitions into ordered levels
        # by iterating over all non-empty subset
        # sequences
        _wo_recurse!(result, Vector{Vector{T}}(), items)
        result
    end

    function _wo_recurse!(
        result, prefix, remaining::Vector{T},
    ) where T
        if isempty(remaining)
            push!(result, copy(prefix))
            return
        end
        # Each level is a non-empty subset of remaining
        for mask in 1:(2^length(remaining) - 1)
            level = T[]
            rest = T[]
            for (i, item) in enumerate(remaining)
                if (mask >> (i - 1)) & 1 == 1
                    push!(level, item)
                else
                    push!(rest, item)
                end
            end
            push!(prefix, sort(level))
            _wo_recurse!(result, prefix, rest)
            pop!(prefix)
        end
    end

    function _steps_for_ordering(
        all_group_steps, ordering, met_set,
    )
        # Given a weak ordering (levels of metabolites
        # from ONE side — either substrates or products)
        # and the full set of steps, select the binding
        # steps for those metabolites consistent with
        # this ordering.
        #
        # A metabolite M at level k can bind to any form
        # whose existing metabolites (from met_set) are
        # all from levels 1..k-1 or the same level k.
        selected = Set{StepKey}()
        accessible = Set{Symbol}()
        for level in ordering
            for met in level
                for sk in all_group_steps
                    lhs_mets = [
                        s for s in sk[1] if s ∈ met_names
                    ]
                    length(lhs_mets) == 1 || continue
                    lhs_mets[1] == met || continue
                    form = [
                        s for s in sk[1] if s ∉ met_names
                    ]
                    length(form) == 1 || continue
                    # Extract metabolites from met_set
                    # that are in this form name
                    form_met_parts = [
                        Symbol(p) for p in split(
                            replace(string(form[1]),
                                    "Estar" => "E"),
                            "_")
                        if Symbol(p) ∈ met_set
                    ]
                    if all(m ∈ accessible
                           for m in form_met_parts)
                        push!(selected, sk)
                    end
                end
            end
            union!(accessible, level)
        end
        selected
    end

    # --- Build topologies ---
    result = MechanismSpec[]
    for (iso_pat, group_paths) in iso_groups
        # Collect all steps from all paths in this group
        all_group_steps = Set{StepKey}()
        step_dict = Dict{StepKey, StepSpec}()
        for path in group_paths
            for s in path
                sk = _step_key(s)
                push!(all_group_steps, sk)
                step_dict[sk] = s
            end
        end

        # Identify metabolites in binding steps, split
        # by half-reaction (substrate-side vs product-side
        # based on form prefix E vs Estar)
        sub_binding_mets = Set{Symbol}()
        prod_binding_mets = Set{Symbol}()
        for sk in all_group_steps
            lhs_mets = [
                s for s in sk[1] if s ∈ met_names
            ]
            length(lhs_mets) == 1 || continue
            met = lhs_mets[1]
            form = [s for s in sk[1] if s ∉ met_names]
            length(form) == 1 || continue
            if met ∈ Set(sub_names)
                push!(sub_binding_mets, met)
            else
                push!(prod_binding_mets, met)
            end
        end

        sub_orderings = _weak_orderings(
            sort(collect(sub_binding_mets)))
        prod_orderings = _weak_orderings(
            sort(collect(prod_binding_mets)))

        # Always include all isomerization steps
        iso_keys = Set{StepKey}()
        for sk in all_group_steps
            lhs_mets = [s for s in sk[1] if s ∈ met_names]
            rhs_mets = [s for s in sk[2] if s ∈ met_names]
            if isempty(lhs_mets) && isempty(rhs_mets)
                push!(iso_keys, sk)
            end
        end

        sub_met_set = Set(sub_names)
        prod_met_set = Set(prod_names)

        seen_topos = Set{Set{StepKey}}()
        for sub_ord in sub_orderings
            for prod_ord in prod_orderings
                # Process substrate and product orderings
                # INDEPENDENTLY — each only checks its
                # own metabolite set for accessibility
                sub_keys = _steps_for_ordering(
                    all_group_steps, sub_ord, sub_met_set,
                )
                prod_keys = _steps_for_ordering(
                    all_group_steps, prod_ord, prod_met_set,
                )
                topo_keys = union(iso_keys, sub_keys,
                    prod_keys)
                topo_keys ∈ seen_topos && continue
                push!(seen_topos, topo_keys)

                steps = [step_dict[sk] for sk in topo_keys]
                sort!(steps; by=s -> (
                    length(s.reactants) == 1 ? 1 : 0,
                    join(sort(s.reactants), "_")
                ))

                # Default RE/SS: first iso is SS, rest RE
                iso_idx = findfirst(
                    s -> length(s.reactants) == 1, steps
                )
                tagged = [
                    StepSpec(
                        s.reactants, s.products,
                        i != iso_idx
                    )
                    for (i, s) in enumerate(steps)
                ]

                # Compute param_count
                form_names = Set{Symbol}()
                for s in tagged
                    union!(form_names, s.reactants)
                    union!(form_names, s.products)
                end
                setdiff!(form_names, met_names)
                n_forms = length(form_names)
                n_steps = length(tagged)
                n_cycles = n_steps - n_forms + 1
                n_re = count(s -> s.is_equilibrium, tagged)
                n_ss = n_steps - n_re
                n_thermo = n_cycles
                param_count = n_re + 2 * n_ss -
                    n_thermo + 2

                push!(result, MechanismSpec(
                    reaction, tagged,
                    ParamConstraint[], param_count
                ))
            end
        end
    end
    result
```

- [ ] **Step 4: Update existing ter-ter test**

The test at line 522 expects 3969. Update it:

```julia
@test length(topos) == 169
```

Also update the dead-end test bounds (line 531) if needed — the topology
structure has changed so bounds may differ.

- [ ] **Step 5: Run tests to verify they pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All pass, including bi-bi (9) and ter-ter (169)

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Fix path combining to use weak orderings

Replace arbitrary 2^n subset combining with weak-ordering
enumeration (Fubini numbers). For ter-ter: 3969 → 169
topologies. For bi-bi: unchanged at 9."
```

---

### Task 4: Allow empty-residual ping-pong (C4)

Remove the `isempty(residual) && continue` check that prevents
conformational ping-pong without covalently attached atoms.

**Files:**
- Modify: `src/mechanism_enumeration.jl` (line 343 in backtracking)
- Test: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing test**

```julia
@testset "empty-residual ping-pong" begin
    # Simple ter-ter: A[C], B[N], D[X]
    # With empty residual allowed, ping-pong paths exist
    # (e.g., E_A → Estar_P with empty residual)
    topos = EnzymeRates._catalytic_topologies(
        ter_ter_rxn)
    # Check that at least some topologies have Estar
    # forms (ping-pong)
    has_estar = any(topos) do spec
        any(spec.steps) do s
            any(sym -> startswith(string(sym), "Estar"),
                Iterators.flatten(
                    (s.reactants, s.products)))
        end
    end
    @test has_estar
end
```

- [ ] **Step 2: Run test — verify it fails**

Expected: FAIL — current code has 0 ping-pong paths for simple ter-ter

- [ ] **Step 3: Remove the empty-residual check**

In `src/mechanism_enumeration.jl`, in the backtracking function
(around line 343), remove or comment out:

```julia
# Remove this line:
isempty(residual) && continue
```

- [ ] **Step 4: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All pass. The ter-ter count will change from 169 (need to update
the test from Task 3 — the new count includes both sequential and
ping-pong topologies).

- [ ] **Step 5: Update topology count tests**

After allowing empty-residual ping-pong, new paths appear. The exact count
depends on which other constraints are active. Update the `@test` assertions
to match the actual count. For ter-ter with weak orderings + empty residual
(but no other new constraints yet): verify the count experimentally and
update.

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Allow empty-residual ping-pong (conformational change)

Enzymes can be in a different conformation (Estar) after
isomerization even without covalently attached atoms."
```

---

### Task 5: Add isomerization constraints (C5, C6, C7, C8)

Add the remaining backtracking constraints: max bound metabolites,
isomerization size limit (≤3 per side), substrate participation required,
and product-only iso form names.

**Files:**
- Modify: `src/mechanism_enumeration.jl` (backtracking function and
  `_form_name` calls)
- Test: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing test for isomerization constraints**

```julia
@testset "isomerization constraints" begin
    topos = EnzymeRates._catalytic_topologies(
        ter_ter_rxn)

    # C5: max bound metabolites = max(3,3) = 3
    met_names = Set([:A, :B, :D, :P, :Q, :R])
    for spec in topos
        for s in spec.steps
            for sym_list in (s.reactants, s.products)
                for sym in sym_list
                    str = replace(
                        string(sym), "Estar" => "E"
                    )
                    parts = split(str, "_")
                    n_mets = count(
                        p -> Symbol(p) ∈ met_names,
                        parts,
                    )
                    @test n_mets <= 3
                end
            end
        end
    end

    # C7: Isomerization requires substrate participation.
    # Every iso source form must contain at least one
    # substrate name. Forms with only products or bare
    # Estar cannot isomerize.
    sub_names_for_c7 = Set([:A, :B, :D])
    for spec in topos
        for s in spec.steps
            if length(s.reactants) == 1 &&
                    length(s.products) == 1
                src = string(s.reactants[1])
                src_parts = Symbol.(split(
                    replace(src, "Estar" => "E"), "_"
                ))
                has_sub = any(
                    p -> p ∈ sub_names_for_c7, src_parts
                )
                @test has_sub "Iso source $src has no substrate"
            end
        end
    end

    # C8: iso forms should not contain substrate names
    # on the product side of an isomerization
    sub_names_set = Set([:A, :B, :D])
    for spec in topos
        for s in spec.steps
            if length(s.reactants) == 1 &&
                    length(s.products) == 1
                dst = string(s.products[1])
                if startswith(dst, "Estar")
                    dst_parts = Symbol.(
                        split(dst[7:end], "_"))
                    for p in dst_parts
                        @test p ∉ sub_names_set "Iso product form $dst contains substrate $p"
                    end
                end
            end
        end
    end
end
```

- [ ] **Step 2: Run test — verify it fails**

Expected: FAIL — current iso forms include substrate names (e.g.,
`Estar_A_P`)

- [ ] **Step 3: Implement constraints in backtracking**

Modify the backtracking in `_catalytic_topologies`:

**C5 — max bound metabolites**: Add `max_bound = max(length(sub_names), length(prod_names))` before the backtracking function. In all substrate binding branches, check `length(on_enzyme_subs) < max_bound` before binding more.

**C6 — isomerization size limit**: At isomerization points, compute `n_prods_eff = k + (isempty(residual) ? 0 : 1)` and check `length(on_enzyme_subs) <= 3 && n_prods_eff <= 3`.

**C7 — substrate participation**: In Case 3 (no subs bound, has residual), remove the final-isomerize branch. Only allow binding remaining substrates.

**C8 — product-only iso forms**: Change all `_form_name(on_enzyme_subs, prods, true)` calls to `_form_name(Symbol[], prods, true)` in ping-pong isomerization branches.

**C9 — multi-product release**: Add a branch that enumerates product subsets of size 1..k at isomerization using `combinations()`, checking atom conservation for each subset. Release products one at a time after iso.

These changes are substantial — implement them incrementally, running tests after each sub-change.

- [ ] **Step 4: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add isomerization constraints C5-C9

- Max bound metabolites: max(n_subs, n_prods)
- Iso size limit: n_subs <= 3 AND n_prods_eff <= 3
- Estar requires substrate for isomerization
- Iso forms contain only products
- Multi-product release at isomerization"
```

---

### Task 6: Fix `_bound_metabolites_at_forms` for partial isomerization

The current implementation (lines 666-681) assumes isomerization converts
ALL substrates to ALL products (`setdiff(bound[from], sub_names) ∪
prod_names`). For ping-pong (partial isomerization), `E_A_B → Estar_P`
only produces P — but the function computes `{P, Q, R}`. Fix to determine
bound metabolites from binding/release edges connected to each form.

**Files:**
- Modify: `src/mechanism_enumeration.jl` (function
  `_bound_metabolites_at_forms`, around line 639-686)
- Test: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing test**

```julia
@testset "bound metabolites for ping-pong" begin
    # Build a simple ping-pong topology manually and
    # check that bound metabolites are correct
    topos = EnzymeRates._catalytic_topologies(
        ter_ter_rxn)
    # Find a ping-pong topology (has Estar forms)
    pp_topo = nothing
    for spec in topos
        has_estar = any(spec.steps) do s
            any(sym -> startswith(string(sym), "Estar"),
                Iterators.flatten(
                    (s.reactants, s.products)))
        end
        if has_estar
            pp_topo = spec
            break
        end
    end
    @test pp_topo !== nothing

    bound = EnzymeRates._bound_metabolites_at_forms(
        pp_topo, ter_ter_rxn)

    # For any form, the number of bound metabolites
    # should not exceed max(n_subs, n_prods) = 3
    for (form, mets) in bound
        @test length(mets) <= 3 "Form $form has $(length(mets)) bound: $mets"
    end

    # Estar_P should have only P bound, not {P, Q, R}
    if haskey(bound, :Estar_P)
        @test :P ∈ bound[:Estar_P]
        @test :Q ∉ bound[:Estar_P]
        @test :R ∉ bound[:Estar_P]
    end
end
```

- [ ] **Step 2: Run test — verify it fails**

Expected: FAIL — `bound[:Estar_P]` contains `{P, Q, R}` instead of `{P}`

- [ ] **Step 3: Fix `_bound_metabolites_at_forms`**

Replace the isomerization handler. Instead of assuming all subs → all prods,
determine bound metabolites at each form from the binding/release steps
connected to it:

```julia
# For each form, find which metabolites have binding
# steps TO this form. A metabolite M is bound at form F
# if there exists a step [F, M] ⇌ [F'] or [F', M] ⇌ [F]
# (i.e., M can bind to or unbind from F).
function _bound_metabolites_at_forms(
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReaction),
)
    subs = substrates(reaction)
    prods = products(reaction)
    met_names = Set{Symbol}()
    for (name, _) in subs
        push!(met_names, name)
    end
    for (name, _) in prods
        push!(met_names, name)
    end

    # Collect all form names
    form_names = Set{Symbol}()
    for s in spec.steps
        for sym in s.reactants
            sym ∉ met_names && push!(form_names, sym)
        end
        for sym in s.products
            sym ∉ met_names && push!(form_names, sym)
        end
    end

    # For each form, find metabolites with binding edges
    bound = Dict{Symbol, Set{Symbol}}()
    for form in form_names
        bound[form] = Set{Symbol}()
    end
    for s in spec.steps
        # Binding step: [form, met] ⇌ [form']
        for (side_a, side_b) in [
            (s.reactants, s.products),
            (s.products, s.reactants),
        ]
            if length(side_a) == 2
                forms_in_a = [
                    x for x in side_a if x ∉ met_names
                ]
                mets_in_a = [
                    x for x in side_a if x ∈ met_names
                ]
                if length(forms_in_a) == 1 &&
                        length(mets_in_a) == 1
                    # met binds to form to produce side_b
                    product_form = side_b[1]
                    push!(bound[product_form],
                        mets_in_a[1])
                    # Also: met is bound at product_form
                    # AND all metabolites already at
                    # forms_in_a[1] are also at
                    # product_form
                end
            end
        end
    end

    # Propagate: if form F2 is produced by binding M to
    # F1, then bound[F2] = bound[F1] ∪ {M}
    changed = true
    while changed
        changed = false
        for s in spec.steps
            for (side_a, side_b) in [
                (s.reactants, s.products),
                (s.products, s.reactants),
            ]
                if length(side_a) == 2 &&
                        length(side_b) == 1
                    mets_in_a = [
                        x for x in side_a if x ∈ met_names
                    ]
                    forms_in_a = [
                        x for x in side_a if x ∉ met_names
                    ]
                    length(mets_in_a) == 1 || continue
                    length(forms_in_a) == 1 || continue
                    product_form = side_b[1]
                    source_form = forms_in_a[1]
                    expected = union(
                        bound[source_form],
                        Set([mets_in_a[1]]),
                    )
                    if !issubset(expected,
                                 bound[product_form])
                        union!(bound[product_form],
                            expected)
                        changed = true
                    end
                end
            end
        end
    end
    bound
end
```

- [ ] **Step 4: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All pass

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Fix _bound_metabolites_at_forms for partial isomerization

The function now determines bound metabolites from binding/
release edges rather than assuming all subs → all prods at
isomerization. Fixes incorrect dead-end expansion for
ping-pong mechanisms."
```

---

### Task 7: Verify pyruvate carboxylase mechanism is generated (was Task 6)

Verify that the known pyruvate carboxylase mechanism appears among the
generated catalytic topologies: ATP+HCO₃ bind, isomerize releasing
ADP+Pi (CO₂ residual), then Pyr binds to Estar and converts to OAA.

**Files:**
- Test: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write test**

```julia
@testset "pyruvate carboxylase mechanism" begin
    topos = EnzymeRates._catalytic_topologies(
        pyruvate_carboxylase_rxn)

    # The known mechanism has 2 isomerization steps:
    # 1. E_ATP_HCO3 → Estar_ADP_Pi
    #    (ATP+HCO3 convert, ADP+Pi produced, CO2 residual)
    # 2. Estar_Pyr → E_OAA
    #    (Pyr+CO2_residual → OAA, cycle complete)
    #
    # Note: Estar_ADP_Pi uses product-only naming (C8).
    # The form names depend on sort order.

    # Find topology containing both iso steps
    target_iso_1 = Set([
        (sort([:E_ATP_HCO3]), sort([:Estar_ADP_Pi])),
        (sort([:E_HCO3_ATP]), sort([:Estar_ADP_Pi])),
    ])
    # Try both possible form name orderings
    found = false
    for spec in topos
        iso_steps = Set{Tuple{Vector{Symbol},Vector{Symbol}}}()
        for s in spec.steps
            if length(s.reactants) == 1 &&
                    length(s.products) == 1
                push!(iso_steps,
                    (sort(s.reactants), sort(s.products)))
            end
        end
        # Check for the PC mechanism iso pattern
        has_atp_hco3_iso = any(iso_steps) do (r, p)
            r == [Symbol("E_ATP_HCO3")] &&
                p == [Symbol("Estar_ADP_Pi")]
        end
        has_pyr_iso = any(iso_steps) do (r, p)
            r == [Symbol("Estar_Pyr")] &&
                p == [Symbol("E_OAA")]
        end
        if has_atp_hco3_iso && has_pyr_iso
            found = true
            break
        end
    end
    @test found "Pyruvate carboxylase mechanism not found"

    # Verified count (scripts/verify_counts.py):
    # 312 = 169 seq + 99 bi-uni + 36 uni-bi + 8 hexa-uni
    @test length(topos) == 312
    seq_count = count(topos) do spec
        n_iso = count(spec.steps) do s
            length(s.reactants) == 1 && length(s.products) == 1
        end
        n_iso == 1
    end
    pp_count = length(topos) - seq_count
    @test seq_count == 169  # 13 × 13 weak orderings
    @test pp_count == 143   # 99 + 36 + 8
end
```

- [ ] **Step 2: Run test**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep "pyruvate carboxylase"`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Add test verifying pyruvate carboxylase mechanism is generated"
```

---

### Task 8: Verify pyruvate dehydrogenase mechanism is generated

The PDH mechanism is a hexa-uni ping-pong: Pyr→CO₂ (residual C2H3O),
then CoA→AcCoA (residual H), then NAD→NADH (no residual). This tests
the empty-residual ping-pong (last step has empty residual).

**Files:**
- Test: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write test**

```julia
@testset "pyruvate dehydrogenase mechanism" begin
    topos = EnzymeRates._catalytic_topologies(
        pyruvate_dehydrogenase_rxn)

    # The known mechanism has 3 iso steps (hexa-uni):
    # 1. E_Pyr → Estar_CO2
    #    (Pyr → CO2, residual C2H3O)
    # 2. Estar_CoA → Estar_AcCoA
    #    (CoA+residual → AcCoA, residual H)
    # 3. Estar_NAD → E_NADH
    #    (NAD+residual → NADH, no residual)
    found = false
    for spec in topos
        iso_steps = Set{Tuple{Vector{Symbol},Vector{Symbol}}}()
        for s in spec.steps
            if length(s.reactants) == 1 &&
                    length(s.products) == 1
                push!(iso_steps,
                    (sort(s.reactants), sort(s.products)))
            end
        end
        has_pyr = any(iso_steps) do (r, p)
            r == [:E_Pyr] && p == [:Estar_CO2]
        end
        has_coa = any(iso_steps) do (r, p)
            r == [:Estar_CoA] && p == [:Estar_AcCoA]
        end
        has_nad = any(iso_steps) do (r, p)
            r == [:Estar_NAD] && p == [:E_NADH]
        end
        if has_pyr && has_coa && has_nad
            found = true
            break
        end
    end
    @test found "Pyruvate dehydrogenase mechanism not found"

    # Verified count (scripts/verify_counts.py):
    # 334 = 169 seq + 117 bi-uni + 36 uni-bi + 12 hexa-uni
    @test length(topos) == 334
    seq_count = count(topos) do spec
        n_iso = count(spec.steps) do s
            length(s.reactants) == 1 && length(s.products) == 1
        end
        n_iso == 1
    end
    pp_count = length(topos) - seq_count
    @test seq_count == 169  # 13 × 13 weak orderings
    @test pp_count == 165   # 117 + 36 + 12
end
```

- [ ] **Step 2: Run test**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep "pyruvate dehydrogenase"`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Add test verifying pyruvate dehydrogenase mechanism is generated"
```

---

### Task 9: Verify topology counts and update existing tests

Update all existing tests that assert specific topology counts to match
the new constrained behavior.

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Determine actual counts**

Run the topology enumeration for all test reactions and record counts:

```julia
julia --project -e '
include("test/test_mechanism_enumeration.jl")
for (name, rxn) in [
    ("bi-bi", bi_bi_rxn),
    ("ter-ter", ter_ter_rxn),
    ("ter-bi", ter_bi_rxn),
    ("PC", pyruvate_carboxylase_rxn),
    ("PDH", pyruvate_dehydrogenase_rxn),
]
    topos = EnzymeRates._catalytic_topologies(rxn)
    println("$name: $(length(topos)) topologies")
end
'
```

- [ ] **Step 2: Update all count assertions**

Update `@test length(topos) == N` assertions. Independently verified
counts (scripts/verify_counts.py):

| Reaction | Total | Sequential | Ping-pong |
|---|---|---|---|
| bi-bi | 9 | 9 | 0 |
| ter-ter (A[C],B[N],D[X]) | 283 | 169 | 114 (81+27+6) |
| pyruvate carboxylase | 312 | 169 | 143 (99+36+8) |
| pyruvate dehydrogenase | 334 | 169 | 165 (117+36+12) |

- [ ] **Step 3: Update dead-end expansion tests**

The dead-end expansion tests (line 524-540) reference specific topologies.
After the constraint changes, the topology structure is different. Update
the dead-end tests to use representative topologies from the new set and
verify bounds.

- [ ] **Step 4: Run full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass

- [ ] **Step 5: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Update topology count assertions for new constraints

bi-bi: 9 (unchanged)
ter-ter: 283 (was 3969)
Pyruvate carboxylase: 312
Pyruvate dehydrogenase: 334"
```

---

### Task 10: Clean up brainstorming scripts

Remove the estimation scripts created during brainstorming that are not
part of the package.

**Files:**
- Delete: `scripts/estimate_constraint_impact.jl`
- Delete: `scripts/estimate_with_empty_residual.jl`
- Delete: `scripts/show_terter_paths.jl`
- Delete: `scripts/estimate_init_mechanisms.jl`
- Delete: `scripts/relaxed_constraint.jl`
- Delete: `scripts/multi_product_release.jl`
- Delete: `scripts/final_estimate.jl`
- Delete: `scripts/iso_pattern_analysis.jl`
- Delete: `scripts/count_combined.jl`
- Delete: `scripts/bystander_count_only.jl`
- Delete: `scripts/bystander_quick.jl`
- Delete: `scripts/estimate_bystander_impact.jl`
- Delete: `scripts/debug_pingpong.jl`
- Delete: `scripts/show_49_topologies.jl`
- Delete: `scripts/pyruvate_carboxylase.jl`

- [ ] **Step 1: Remove scripts**

```bash
rm -f scripts/estimate_*.jl scripts/show_*.jl \
      scripts/relaxed_constraint.jl \
      scripts/multi_product_release.jl \
      scripts/final_estimate.jl \
      scripts/iso_pattern_analysis.jl \
      scripts/count_combined.jl \
      scripts/bystander_*.jl \
      scripts/debug_pingpong.jl \
      scripts/pyruvate_carboxylase.jl
```

- [ ] **Step 2: Verify no scripts are referenced**

```bash
grep -r "scripts/" test/ src/ || echo "No references"
```

- [ ] **Step 3: Commit**

```bash
git add -u scripts/
git commit -m "Remove brainstorming estimation scripts"
```
