# Finish-Refactor Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the remaining gaps on the concrete-types refactor (restore dropped validation, fix three enumeration/doc bugs, dedup load-bearing code, reorganize three migration-era test files, sweep stale comments) so the branch is PR-ready.

**Architecture:** Four ordered clusters. **A** restores `EnzymeReaction`/`ReactantAtoms` constructor validation and fixes the regulator/multiplicity/README bugs (TDD, may surface fixture breakage — handle first). **B** is behavior-preserving dedup. **C** moves three test files into their ~1:1 src-correspondence homes. **D** sweeps stale comments/docs. Full `Pkg.test()` is green after each cluster.

**Tech Stack:** Julia 1.9+; package `EnzymeRates`. Test with `julia --project -e 'using Pkg; Pkg.test()'`. Single-test-file runs are NOT supported standalone (files depend on shared fixtures included by `runtests.jl`); use the full suite or a targeted `julia --project -e` snippet that `using EnzymeRates` + reproduces the case.

**Design doc:** `docs/superpowers/specs/2026-05-31-finish-refactor-fixes-design.md`

**Reference baselines (verified on branch `refactor-to-concrete-types-instead-of-symbols`):**
- Full suite passes: 27073 pass / 1 broken (deliberate, out of scope) / 0 fail, ~9m34s.
- Export count is 18 (do NOT change).
- Perf gate: `rate_equation` must stay `allocs == 0` and `t < 100e-9`.

**Memory note:** a prior agent's full run was OOM-killed (~2.1 GB RSS) at the compile-budget section; the reviewer's completed fine. Watch RSS; if a run dies, it's likely environment memory pressure, not a regression — re-run before assuming a code defect.

---

## How to run one mechanism's rate_equation by hand (used in several verify steps)

```bash
julia --project -e '
using EnzymeRates
const ER = EnzymeRates
# build whatever reaction/mechanism the step needs, then call ER.<fn>(...)
'
```

The shared test fixtures (`MECHANISM_TEST_SPECS`, etc.) only exist inside the `Pkg.test()`
process via `test/mechanism_definitions_for_test_enzyme_derivation.jl`. When a step says
"run the full suite," it means `julia --project -e 'using Pkg; Pkg.test()'`.

---

# CLUSTER A — Bug fixes (PR blockers). Do first; TDD.

### Task A1: `ReactantAtoms` constructor validates atoms

**Files:**
- Modify: `src/types.jl:224-231` (the `ReactantAtoms` struct + inner ctor)
- Test: `test/test_types.jl` (the `EnzymeReaction (new concrete)` testset region, ~line 834 — note this label is renamed in Task D1; add tests near the reaction-construction tests)

- [ ] **Step 1: Write the failing tests**

Add to `test/test_types.jl` inside the reaction-construction testset (search for `@testset` containing `ReactantAtoms` or `EnzymeReaction`; if none, add a new `@testset "ReactantAtoms validation" begin ... end` near the other type tests):

```julia
@testset "ReactantAtoms validation" begin
    # Mandatory atoms: empty atom list is rejected.
    @test_throws ErrorException EnzymeRates.ReactantAtoms(
        EnzymeRates.Substrate(:S), Pair{Symbol,Int}[])
    # Positive counts: zero / negative rejected.
    @test_throws ErrorException EnzymeRates.ReactantAtoms(
        EnzymeRates.Substrate(:S), [:C => 0])
    @test_throws ErrorException EnzymeRates.ReactantAtoms(
        EnzymeRates.Substrate(:S), [:C => -1])
    # Bool is not a valid count (true === 1 in Julia, must be excluded).
    @test_throws ErrorException EnzymeRates.ReactantAtoms(
        EnzymeRates.Substrate(:S), Pair{Symbol,Int}[:C => true])
    # Valid construction still works.
    ra = EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 6, :H => 12])
    @test EnzymeRates.atoms(ra) == [:C => 6, :H => 12]
end
```

- [ ] **Step 2: Run to verify failure**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: the new `@test_throws` lines FAIL (no error currently thrown), valid construction PASSES.

- [ ] **Step 3: Add validation to the constructor**

Replace `src/types.jl:224-231` (`ReactantAtoms` struct) inner constructor body. Current:

```julia
    function ReactantAtoms(metabolite::Reactant,
                           atoms::Vector{<:Pair{Symbol,<:Integer}})
        new(metabolite, sort(Vector{Pair{Symbol,Int}}(atoms); by=first))
    end
```

New:

```julia
    function ReactantAtoms(metabolite::Reactant,
                           atoms::Vector{<:Pair{Symbol,<:Integer}})
        isempty(atoms) && error(
            "ReactantAtoms: $(name(metabolite)) has no declared atoms; atoms " *
            "are mandatory (use `[C…]` bracket syntax in @enzyme_reaction).")
        for (elem, count) in atoms
            count isa Integer && !(count isa Bool) && count > 0 ||
                error("ReactantAtoms: $(name(metabolite)) has non-positive " *
                      "atom count for element $elem ($count); counts must be " *
                      "positive integers (not Bool).")
        end
        new(metabolite, sort(Vector{Pair{Symbol,Int}}(atoms); by=first))
    end
```

- [ ] **Step 4: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: the A1 testset PASSES. **If other fixtures now fail** because they built `ReactantAtoms` with empty atoms, that is expected breakage — note which, they get repaired in Task A2/A3 or, if a fixture is a deliberately atom-less toy, STOP and ask Denis.

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "fix: ReactantAtoms constructor validates mandatory positive atoms"
```

---

### Task A2: `EnzymeReaction` constructor validates structure + name identity

**Files:**
- Modify: `src/types.jl:271-282` (the `EnzymeReaction` inner ctor)
- Test: `test/test_types.jl` (reaction-construction testset)

This restores `main`'s non-empty / duplicate-name / atom-balance checks AND adds the
cross-category name-uniqueness that makes finding #2's regulator/substrate collision
impossible by construction.

- [ ] **Step 1: Write the failing tests**

Add to `test/test_types.jl`:

```julia
@testset "EnzymeReaction validation" begin
    S = EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1])
    P = EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])
    P2 = EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])
    Punbal = EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 2])
    noregs = EnzymeRates.RegulatorMults[]
    regS = EnzymeRates.RegulatorMults(EnzymeRates.CompetitiveInhibitor(:S), [1])
    regA = EnzymeRates.RegulatorMults(EnzymeRates.CompetitiveInhibitor(:A), [1])
    regA2 = EnzymeRates.RegulatorMults(EnzymeRates.CompetitiveInhibitor(:A), [1])

    # Empty substrate set rejected.
    @test_throws ErrorException EnzymeRates.EnzymeReaction(
        EnzymeRates.ReactantAtoms[], [P], Int[1])
    # Empty product set rejected.
    @test_throws ErrorException EnzymeRates.EnzymeReaction(
        [S], EnzymeRates.ReactantAtoms[], Int[1])
    # Duplicate product names rejected.
    @test_throws ErrorException EnzymeRates.EnzymeReaction(
        [S], [P, P2], noregs, Int[1])
    # Atom imbalance rejected (S has 1 C, Punbal has 2 C).
    @test_throws ErrorException EnzymeRates.EnzymeReaction(
        [S], [Punbal], noregs, Int[1])
    # Duplicate regulator names rejected.
    @test_throws ErrorException EnzymeRates.EnzymeReaction(
        [S], [P], [regA, regA2], Int[1])
    # Cross-category name collision rejected (regulator named :S == substrate :S).
    @test_throws ErrorException EnzymeRates.EnzymeReaction(
        [S], [P], [regS], Int[1])
    # Valid balanced reaction with a distinct-named regulator still constructs.
    rxn = EnzymeRates.EnzymeReaction([S], [P], [regA], Int[1])
    @test EnzymeRates.name.(EnzymeRates.substrates(rxn)) == [:S]
end
```

- [ ] **Step 2: Run to verify failure**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: the new `@test_throws` lines FAIL (no error currently); valid construction PASSES.

- [ ] **Step 3: Add validation to the constructor**

Replace the `EnzymeReaction` inner constructor body at `src/types.jl:271-282`. Current:

```julia
    function EnzymeReaction(reactants::Vector{ReactantAtoms},
                            regulators::Vector{RegulatorMults},
                            allowed_catalytic_multiplicities::Vector{Int})
        sorted_reactants = sort(reactants; by = ra -> name(metabolite(ra)))
        sorted_regulators = sort(regulators; by = rm -> name(regulator(rm)))
        sorted_mults = sort(unique(allowed_catalytic_multiplicities))
        all(m -> m ≥ 1, sorted_mults) ||
            error("EnzymeReaction: allowed_catalytic_multiplicities must " *
                  "all be ≥ 1, got $allowed_catalytic_multiplicities")
        new(sorted_reactants, sorted_regulators, sorted_mults)
    end
```

New:

```julia
    function EnzymeReaction(reactants::Vector{ReactantAtoms},
                            regulators::Vector{RegulatorMults},
                            allowed_catalytic_multiplicities::Vector{Int})
        sorted_reactants = sort(reactants; by = ra -> name(metabolite(ra)))
        sorted_regulators = sort(regulators; by = rm -> name(regulator(rm)))
        sorted_mults = sort(unique(allowed_catalytic_multiplicities))
        all(m -> m ≥ 1, sorted_mults) ||
            error("EnzymeReaction: allowed_catalytic_multiplicities must " *
                  "all be ≥ 1, got $allowed_catalytic_multiplicities")

        subs = Reactant[metabolite(ra) for ra in sorted_reactants
                        if metabolite(ra) isa Substrate]
        prods = Reactant[metabolite(ra) for ra in sorted_reactants
                         if metabolite(ra) isa Product]
        isempty(subs)  && error("EnzymeReaction: substrates must not be empty")
        isempty(prods) && error("EnzymeReaction: products must not be empty")

        sub_names  = Symbol[name(m) for m in subs]
        prod_names = Symbol[name(m) for m in prods]
        reg_names  = Symbol[name(regulator(rm)) for rm in sorted_regulators]
        length(sub_names)  == length(Set(sub_names))  ||
            error("EnzymeReaction: duplicate substrate names")
        length(prod_names) == length(Set(prod_names)) ||
            error("EnzymeReaction: duplicate product names")
        length(reg_names)  == length(Set(reg_names))  ||
            error("EnzymeReaction: duplicate regulator names")

        # Names are identities at the mechanism level (concs.X drives one
        # species), so a regulator may not share a substrate/product name.
        met_names = Set{Symbol}(sub_names) ∪ Set{Symbol}(prod_names)
        for rn in reg_names
            rn in met_names && error(
                "EnzymeReaction: regulator name $rn collides with a " *
                "substrate/product name; names must be unique across categories")
        end

        # Atom mass-balance: per-element sum over substrate vs product reactants.
        atom_sum = function (mets)
            d = Dict{Symbol,Int}()
            for ra in sorted_reactants
                metabolite(ra) in mets || continue
                for (elem, c) in atoms(ra)
                    d[elem] = get(d, elem, 0) + c
                end
            end
            d
        end
        sub_atoms  = atom_sum(Set(subs))
        prod_atoms = atom_sum(Set(prods))
        for elem in union(keys(sub_atoms), keys(prod_atoms))
            s_c = get(sub_atoms, elem, 0)
            p_c = get(prod_atoms, elem, 0)
            s_c == p_c || error(
                "EnzymeReaction: atom imbalance — element $elem appears " *
                "$s_c time(s) on substrate side and $p_c on product side. " *
                "Declared atoms must balance.")
        end

        new(sorted_reactants, sorted_regulators, sorted_mults)
    end
```

- [ ] **Step 4: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: A2 testset PASSES. Watch for fixture breakage from the new checks — repair
unbalanced/duplicate fixtures, or STOP and ask Denis if a fixture is a deliberately
unbalanced toy. The atom-balance check is the most likely to surface issues.

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "fix: EnzymeReaction validates non-empty sets, unique names, atom balance"
```

---

### Task A3: `_add_competitive_inhibitor` is idempotent on existing regulators

**Files:**
- Modify: `src/mechanism_enumeration.jl:1230-1235`
- Test: `test/test_mechanism_enumeration.jl` (near other `_add_competitive_inhibitor` / dead-end tests)

- [ ] **Step 1: Write the failing test**

Add to `test/test_mechanism_enumeration.jl`:

```julia
@testset "_add_competitive_inhibitor idempotent on existing regulator" begin
    S = EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1])
    P = EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])
    regA = EnzymeRates.RegulatorMults(EnzymeRates.CompetitiveInhibitor(:A), [1])
    rxn = EnzymeRates.EnzymeReaction([S], [P], [regA], Int[1])
    # Adding :A again must NOT create a duplicate regulator (ctor would reject
    # the dup anyway; the helper must be idempotent instead of throwing).
    rxn2 = EnzymeRates._add_competitive_inhibitor(rxn, :A)
    @test length(EnzymeRates.regulators(rxn2)) == 1
end
```

- [ ] **Step 2: Run to verify failure**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — currently pushes a duplicate, and after A2 the `EnzymeReaction` ctor throws on the duplicate regulator name.

- [ ] **Step 3: Make the helper idempotent**

Replace `src/mechanism_enumeration.jl:1230-1235`. Current:

```julia
function _add_competitive_inhibitor(rxn::EnzymeReaction, reg_name::Symbol)
    new_regs = copy(regulators(rxn))
    push!(new_regs, RegulatorMults(CompetitiveInhibitor(reg_name), Int[1]))
    EnzymeReaction(copy(reactants(rxn)), new_regs,
                   copy(allowed_catalytic_multiplicities(rxn)))
end
```

New:

```julia
function _add_competitive_inhibitor(rxn::EnzymeReaction, reg_name::Symbol)
    any(rm -> name(regulator(rm)) == reg_name, regulators(rxn)) && return rxn
    new_regs = copy(regulators(rxn))
    push!(new_regs, RegulatorMults(CompetitiveInhibitor(reg_name), Int[1]))
    EnzymeReaction(copy(reactants(rxn)), new_regs,
                   copy(allowed_catalytic_multiplicities(rxn)))
end
```

- [ ] **Step 4: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: A3 testset PASSES; no enumeration counts regress.

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "fix: _add_competitive_inhibitor is idempotent for declared regulators"
```

---

### Task A4: Enumerate catalytic multiplicity in `_expand_to_allosteric`

**Files:**
- Modify: `src/mechanism_enumeration.jl:1445-1462`
- Test: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write the failing test**

Add to `test/test_mechanism_enumeration.jl`:

```julia
@testset "_expand_to_allosteric enumerates all allowed multiplicities" begin
    S = EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1])
    P = EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])
    # Multi-valued allowed multiplicities.
    rxn = EnzymeRates.EnzymeReaction([S], [P], EnzymeRates.RegulatorMults[], Int[2, 4])
    inits = EnzymeRates.init_mechanisms(rxn)
    m = first(inits)
    allo = EnzymeRates._expand_to_allosteric(m, rxn)
    mults = Set(EnzymeRates.catalytic_multiplicity(am) for am in allo)
    @test mults == Set([2, 4])

    # Single-valued case is unchanged (regression guard).
    rxn1 = EnzymeRates.EnzymeReaction([S], [P], EnzymeRates.RegulatorMults[], Int[2])
    allo1 = EnzymeRates._expand_to_allosteric(first(EnzymeRates.init_mechanisms(rxn1)), rxn1)
    @test all(EnzymeRates.catalytic_multiplicity(am) == 2 for am in allo1)
end
```

- [ ] **Step 2: Run to verify failure**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `only(...)` throws `ArgumentError` for the `[2, 4]` reaction.

- [ ] **Step 3: Loop over multiplicities**

Replace `src/mechanism_enumeration.jl:1445-1462`. Current:

```julia
function _expand_to_allosteric(m::Mechanism, rxn::EnzymeReaction)
    cn = only(allowed_catalytic_multiplicities(rxn))
    n_g = length(steps(m))
    base_tags = Symbol[:EqualAI for _ in 1:n_g]
    empty_sites = RegulatorySite[]
    results = AllostericMechanism[]
    push!(results, AllostericMechanism(
        reaction(m), copy(steps(m)), copy(base_tags),
        cn, copy(empty_sites)))
    for g in 1:n_g
        new_tags = copy(base_tags)
        new_tags[g] = :OnlyA
        push!(results, AllostericMechanism(
            reaction(m), copy(steps(m)), new_tags,
            cn, copy(empty_sites)))
    end
    results
end
```

New:

```julia
function _expand_to_allosteric(m::Mechanism, rxn::EnzymeReaction)
    n_g = length(steps(m))
    base_tags = Symbol[:EqualAI for _ in 1:n_g]
    empty_sites = RegulatorySite[]
    results = AllostericMechanism[]
    for cn in allowed_catalytic_multiplicities(rxn)
        push!(results, AllostericMechanism(
            reaction(m), copy(steps(m)), copy(base_tags),
            cn, copy(empty_sites)))
        for g in 1:n_g
            new_tags = copy(base_tags)
            new_tags[g] = :OnlyA
            push!(results, AllostericMechanism(
                reaction(m), copy(steps(m)), new_tags,
                cn, copy(empty_sites)))
        end
    end
    results
end
```

Also update the docstring above the function (`mechanism_enumeration.jl:1437-1443`) — replace
"inherits `rxn`'s oligomeric state as `catalytic_multiplicity`" with "emits the variant set
for each value in `rxn`'s `allowed_catalytic_multiplicities`."

- [ ] **Step 4: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: A4 testset PASSES; existing single-valued enumeration counts unchanged.

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "fix: _expand_to_allosteric enumerates all allowed catalytic multiplicities"
```

---

### Task A5: README — switch the recovery example to `allosteric_regulators: A`

**Files:**
- Modify: `README.md:138-150`
- Test: `test/test_readme_runs.jl` (already runs README code blocks; no new test)

- [ ] **Step 1: Edit the example reaction**

In `README.md`, replace the reaction block:

```julia
rxn = @enzyme_reaction begin
    substrates: S[C]
    products:   P[C]
    competitive_inhibitors: A
    oligomeric_state: 2
end
```

with:

```julia
rxn = @enzyme_reaction begin
    substrates: S[C]
    products:   P[C]
    allosteric_regulators: A
    oligomeric_state: 2
end
```

- [ ] **Step 2: Fix the surrounding prose**

Replace the paragraph that begins `` `competitive_inhibitors: A` declares `A` as a catalytic-site binder; ... `` (README.md ~149-151) with:

```markdown
`allosteric_regulators: A` declares `A` as an allosteric effector; the search
enumerates MWC allosteric variants (which conformations `A` binds and which
kinetic groups are state-dependent) and selects among them by cross-validation
score. (To instead model `A` as a dead-end catalytic-site binder, declare it
with `competitive_inhibitors: A` and the search enumerates dead-end variants.)
```

- [ ] **Step 3: Verify the README still runs**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: `test_readme_runs.jl` PASSES (the `# README-SKIP-IN-TEST` block is still skipped; the reaction-construction line runs).

- [ ] **Step 4: Commit**

```bash
git add README.md
git commit -m "docs: README recovery example uses allosteric_regulators (matches enumeration)"
```

---

# CLUSTER B — Dedup / cleanup (behavior-preserving). Existing suite is the guard.

### Task B1: `metabolites()` lifts via `Mechanism(em)` instead of hand-decoding Sig

**Files:**
- Modify: `src/types.jl:1055-1075` (the `@generated metabolites(::EnzymeMechanism{Sig})`)
- Test: existing suite (the accessor is exercised everywhere); add a focused equality test in `test/test_types.jl`.

- [ ] **Step 1: Write a pin test (current behavior is correct; lock it before refactor)**

Add to `test/test_types.jl`:

```julia
@testset "metabolites() lift equals declaration order" begin
    m = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    @test EnzymeRates.metabolites(m) == (:S, :P)
end
```

- [ ] **Step 2: Run to verify it passes NOW (pin, not red)**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS (this captures current behavior so the refactor can't change it).

- [ ] **Step 3: Rewrite the generated body to lift**

Replace the `@generated function metabolites(::EnzymeMechanism{Sig}) where {Sig}` body
(`src/types.jl:1055-1075`). Current body indexes `Sig[1][1]` / `Sig[1][2]`. New body:

```julia
@generated function metabolites(::EnzymeMechanism{Sig}) where {Sig}
    m = Mechanism(EnzymeMechanism{Sig}())
    rxn = reaction(m)
    names = Symbol[]
    seen = Set{Symbol}()
    for s in substrates(rxn)
        nm = name(s); nm ∉ seen && (push!(seen, nm); push!(names, nm))
    end
    for p in products(rxn)
        nm = name(p); nm ∉ seen && (push!(seen, nm); push!(names, nm))
    end
    for rm in regulators(rxn)
        nm = name(regulator(rm)); nm ∉ seen && (push!(seen, nm); push!(names, nm))
    end
    return Tuple(names)
end
```

Keep the existing docstring (it explains why this stays `@generated` for the hot path) but
remove any sentence describing the old `Sig[1][1]` tuple-indexing if present.

- [ ] **Step 4: Run to verify pass + no `Sig[` left in src**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite PASSES (B1 pin test + all existing).
Run: `grep -rn 'Sig\[' src/`
Expected: zero matches (the only hand-decode is gone).

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "refactor: metabolites() lifts via Mechanism(em); drop Sig hand-decode"
```

---

### Task B2: Extract `_canonicalize_iso_groups` shared by both constructors

**Files:**
- Modify: `src/types.jl` — add helper before the `Mechanism` struct (~line 375); use it in `Mechanism` ctor (`385-390`) and `AllostericMechanism` ctor (`445-452`).
- Test: existing suite (canonicalization is covered by enumeration + naming goldens).

- [ ] **Step 1: Add the helper**

Insert before `# Mechanism: groups elementary steps...` (just above `struct Mechanism`, ~line 375):

```julia
# Canonicalize iso-step storage direction (RE + SS) to physical-forward for
# every group. Shared by the `Mechanism` and `AllostericMechanism`
# constructors so the Canonical Step Form invariant cannot drift between them.
function _canonicalize_iso_groups(reaction::EnzymeReaction,
                                  groups::Vector{Vector{Step}})
    subs  = Set{Symbol}(name(s) for s in substrates(reaction))
    prods = Set{Symbol}(name(s) for s in products(reaction))
    flat0 = Step[s for group in groups for s in group]
    binding_steps = filter(is_binding, flat0)
    [[_canonical_iso_direction(s, subs, prods, binding_steps)
      for s in group] for group in groups]
end
```

- [ ] **Step 2: Use it in the `Mechanism` constructor**

Replace `src/types.jl:385-390` (inside `Mechanism` inner ctor):

```julia
        subs  = Set{Symbol}(name(s) for s in substrates(reaction))
        prods = Set{Symbol}(name(s) for s in products(reaction))
        flat0 = Step[s for group in steps for s in group]
        binding_steps = filter(is_binding, flat0)
        steps = [[_canonical_iso_direction(s, subs, prods, binding_steps)
                  for s in group] for group in steps]
```

with:

```julia
        steps = _canonicalize_iso_groups(reaction, steps)
```

- [ ] **Step 3: Use it in the `AllostericMechanism` constructor**

Replace `src/types.jl:445-452` (the canonicalization block inside `AllostericMechanism`
inner ctor — the comment line `# Canonicalize iso-step storage direction...` through the
`cat_steps = [[...]]` assignment, stopping BEFORE the `# Detect Kreg name collision` block
at line 453):

```julia
        # Canonicalize iso-step storage direction (RE + SS) to physical-
        # forward, mirroring the non-allosteric `Mechanism` constructor.
        subs  = Set{Symbol}(name(s) for s in substrates(reaction))
        prods = Set{Symbol}(name(s) for s in products(reaction))
        flat0 = Step[s for group in cat_steps for s in group]
        binding_steps = filter(is_binding, flat0)
        cat_steps = [[_canonical_iso_direction(s, subs, prods, binding_steps)
                      for s in group] for group in cat_steps]
```

with:

```julia
        cat_steps = _canonicalize_iso_groups(reaction, cat_steps)
```

- [ ] **Step 4: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite PASSES (canonicalization behavior identical; goldens unchanged).

- [ ] **Step 5: Commit**

```bash
git add src/types.jl
git commit -m "refactor: extract _canonicalize_iso_groups shared by both Mechanism ctors"
```

---

### Task B3: Route Step→Parameter emit through one helper

**Files:**
- Move: `_emit_cat_params_for_rep` from `src/rate_eq_derivation.jl:1016-1024` to `src/types.jl` (beside the Parameter family, after `_enumerate_parameters_full`, ~line 1451).
- Modify: `src/types.jl:1376-1394` (`_onlyA_parameters_for_sym`), `1398-1422` (`_all_params_for_sym`), `1433-1451` (`_enumerate_parameters_full`) to call it.
- Test: existing suite (parameter enumeration goldens).

> **Note:** the design also mentioned the two rename-pair builders (`_I_rename_parameters`, `_A_rename_parameters` in `rate_eq_derivation.jl`). Those build *pair* maps (A→I, None→A), a different shape than the emit list. Refactoring them onto the helper requires a variant; do it only if it stays clearly simpler. This task covers the three `types.jl` walkers (the unambiguous win). The rename-pair builders are left as-is unless the executor confirms a clean shared form.

- [ ] **Step 1: Move the helper to types.jl**

Delete `_emit_cat_params_for_rep` from `src/rate_eq_derivation.jl:1016-1024` (keep its
docstring with it). Add to `src/types.jl` immediately after `_enumerate_parameters_full`
(end of file region, ~line 1451):

```julia
"""
Emit the Parameter(s) governing a single kinetic-group representative step
`rep` at allosteric `state`. The 4-way switch on `is_equilibrium(rep)` ×
`is_binding(rep)` is the structural truth of how a Step maps to Parameters;
every parameter-enumeration walker routes through here. Returns 1 element for
RE steps (`Kd`/`Kiso`) and 2 for SS steps (`Kon`+`Koff` / `Kfor`+`Krev`).
"""
function _emit_cat_params_for_rep(rep::Step, state::Symbol)
    if is_equilibrium(rep)
        return Parameter[is_binding(rep) ? Kd(rep, state) : Kiso(rep, state)]
    end
    if is_binding(rep)
        return Parameter[Kon(rep, state), Koff(rep, state)]
    end
    Parameter[Kfor(rep, state), Krev(rep, state)]
end
```

- [ ] **Step 2: Route `_enumerate_parameters_full` through it**

Replace the per-group emit body in `_enumerate_parameters_full` (`src/types.jl:1436-1449`):

```julia
    for group in steps(m)
        rep = _group_rep(group, fes)
        if is_equilibrium(rep)
            push!(out, is_binding(rep) ? Kd(rep, :None) : Kiso(rep, :None))
        else
            if is_binding(rep)
                push!(out, Kon(rep, :None))
                push!(out, Koff(rep, :None))
            else
                push!(out, Kfor(rep, :None))
                push!(out, Krev(rep, :None))
            end
        end
    end
```

with:

```julia
    for group in steps(m)
        append!(out, _emit_cat_params_for_rep(_group_rep(group, fes), :None))
    end
```

- [ ] **Step 3: Route `_onlyA_parameters_for_sym` through it**

Replace the per-group emit body in `_onlyA_parameters_for_sym` (`src/types.jl:1379-1393`):

```julia
    for (g, group) in enumerate(steps(am))
        st = cat_allo_state(am, g) === :EqualAI ? :EqualAI : :A
        rep = _group_rep(group, fes)
        if is_equilibrium(rep)
            push!(out, is_binding(rep) ? Kd(rep, st) : Kiso(rep, st))
        else
            if is_binding(rep)
                push!(out, Kon(rep, st))
                push!(out, Koff(rep, st))
            else
                push!(out, Kfor(rep, st))
                push!(out, Krev(rep, st))
            end
        end
    end
```

with:

```julia
    for (g, group) in enumerate(steps(am))
        st = cat_allo_state(am, g) === :EqualAI ? :EqualAI : :A
        append!(out, _emit_cat_params_for_rep(_group_rep(group, fes), st))
    end
```

- [ ] **Step 4: Route `_all_params_for_sym` through it (catalytic part only)**

In `_all_params_for_sym` (`src/types.jl:1398-1423`), replace the per-group emit body
(the `if is_equilibrium(rep) ... end` block inside the `for (g, group)` loop, lines
1404-1414) with:

```julia
        append!(out, _emit_cat_params_for_rep(_group_rep(group, fes), :I))
```

Leave the `cat_allo_state(am, g) === :OnlyA && continue` guard and the trailing reg-site
loop (`for site in regulatory_sites(am)`) unchanged.

- [ ] **Step 5: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite PASSES — `parameters`, `fitted_params`, `rate_equation_string` goldens
unchanged (the helper produces byte-identical Parameter lists). The chokepoint test must
still pass (helper added to types.jl, before the `name(p,m)` definitions — verify include
order: `_emit_cat_params_for_rep` only constructs Parameters, doesn't render names, so order
is fine).

- [ ] **Step 6: Commit**

```bash
git add src/types.jl src/rate_eq_derivation.jl
git commit -m "refactor: route Step->Parameter emit through shared _emit_cat_params_for_rep"
```

---

### Task B4: Collapse 4 metabolite `_to_sig` encoders into 1

**Files:**
- Modify: `src/types.jl:556-559`
- Test: existing suite (Sig round-trip is exercised by every compiled mechanism); add a round-trip pin.

- [ ] **Step 1: Write a Sig round-trip pin test**

Add to `test/test_types.jl`:

```julia
@testset "_to_sig metabolite encoding round-trips" begin
    for M in (EnzymeRates.Substrate, EnzymeRates.Product,
              EnzymeRates.AllostericRegulator, EnzymeRates.CompetitiveInhibitor)
        sig = EnzymeRates._to_sig(M(:X))
        @test EnzymeRates._metabolite_from_sig(sig) == M(:X)
    end
end
```

- [ ] **Step 2: Run to verify it passes NOW (pin)**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS.

- [ ] **Step 3: Collapse the 4 methods**

Replace `src/types.jl:556-559`:

```julia
_to_sig(s::Substrate)            = (:Substrate, name(s))
_to_sig(p::Product)              = (:Product, name(p))
_to_sig(r::AllostericRegulator)  = (:AllostericRegulator, name(r))
_to_sig(c::CompetitiveInhibitor) = (:CompetitiveInhibitor, name(c))
```

with:

```julia
# One encoder for every Metabolite leaf: (TypeTag, name). The tag Symbol is
# `nameof(typeof(m))`, identical to the four hand-written tags it replaces, so
# the Sig layout is unchanged and `_metabolite_from_sig` still decodes it.
_to_sig(m::Metabolite) = (nameof(typeof(m)), name(m))
```

- [ ] **Step 4: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite PASSES — `nameof(typeof(Substrate(:S)))` is `:Substrate`, so emitted
Sigs are byte-identical to before; the round-trip pin and all compiled-mechanism tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "refactor: single _to_sig(::Metabolite) encoder via nameof type tag"
```

---

### Task B5: Extract `_reject_opaque_bound_forms`; fix wrong-macro error

**Files:**
- Add helper + use in `src/dsl.jl:569-585` (`@enzyme_mechanism` path, the `let ... end` block ending just before `return _build_mechanism_expr(...)`) and `src/dsl.jl:1164-1180` (`@allosteric_mechanism` path, the byte-identical `let ... end` block ending just before `cm_expr = _build_mechanism_expr(...)`).
- Test: `test/test_dsl.jl` (opaque-rejection tests already exist; add an allosteric-path error-message assertion).

> **Confirmed:** both blocks use the local variable `side_terms_per_step` and call only `_is_conformation_shape`; nothing else from the surrounding macro scope. The extraction is clean. The current error string in BOTH blocks says `@enzyme_mechanism` (lines 580 and 1175) — the allosteric one is the wrong-macro bug this task fixes.

- [ ] **Step 1: Write the failing test (wrong macro name)**

Add to `test/test_dsl.jl` (near existing opaque-rejection tests):

```julia
@testset "@allosteric_mechanism opaque rejection names itself" begin
    err = try
        @eval @allosteric_mechanism begin
            substrates: S
            products:   P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                ES <--> EP
            end
        end
        nothing
    catch e
        e
    end
    @test err isa LoadError || err isa ErrorException
    msg = err isa LoadError ? sprint(showerror, err.error) : sprint(showerror, err)
    @test occursin("@allosteric_mechanism", msg)
    @test !occursin("@enzyme_mechanism", msg)
end
```

- [ ] **Step 2: Run to verify failure**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — the allosteric path currently emits `@enzyme_mechanism` in the message.

- [ ] **Step 3: Add the shared helper**

Add to `src/dsl.jl` (near the other parse helpers, before the first macro body that uses it):

```julia
# Reject opaque bound-form bare-enzyme names. A bare-enzyme term `:X` is
# acceptable iff `:X` is a call-form head seen in this steps block (`E` in
# `E(S)`) or matches the conformation shape (`:E`, `:Estar`, `:E_c`).
# Multi-capital (`:ES`) and underscore-then-uppercase (`:E_S`) names are
# opaque and rejected in favor of decomposed call notation. `macro_name`
# names the invoking macro so the error points at the right docs.
function _reject_opaque_bound_forms(side_terms_per_step, macro_name::String)
    call_heads = Set{Symbol}()
    for (_, lhs, rhs, _) in side_terms_per_step
        for t in (lhs..., rhs...)
            t.kind === :call && push!(call_heads, t.conformation)
        end
    end
    for (_, lhs, rhs, _) in side_terms_per_step
        for t in (lhs..., rhs...)
            t.kind === :bare_enzyme || continue
            (t.sym in call_heads || _is_conformation_shape(t.sym)) && continue
            error("$macro_name: `$(t.sym)` looks like an opaque bound-form " *
                  "name; write it as decomposed call notation, e.g. `E(S)` " *
                  "or `E(A, B)`.")
        end
    end
end
```

- [ ] **Step 4: Call it from the `@enzyme_mechanism` path**

Replace the comment block + `let ... end` opaque-rejection block at `src/dsl.jl:563-585`
(from the `# Reject opaque bound-form...` comment through the closing `end` of the `let`,
i.e. the block immediately before `return _build_mechanism_expr(...)`) with:

```julia
    _reject_opaque_bound_forms(side_terms_per_step, "@enzyme_mechanism")
```

- [ ] **Step 5: Call it from the `@allosteric_mechanism` path**

Replace the corresponding comment + `let ... end` block at `src/dsl.jl:1157-1180` (the
byte-identical second copy, immediately before `cm_expr = _build_mechanism_expr(...)`) with:

```julia
    _reject_opaque_bound_forms(side_terms_per_step, "@allosteric_mechanism")
```

- [ ] **Step 6: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: B5 testset PASSES; existing opaque-rejection tests (the `@enzyme_mechanism` path)
still PASS.

- [ ] **Step 7: Commit**

```bash
git add src/dsl.jl test/test_dsl.jl
git commit -m "refactor: share opaque-bound-form rejection; fix @allosteric_mechanism error name"
```

---

### Task B6 + C1: Move chokepoint test into test_types.jl and broaden its classifier

**Files:**
- Modify: `test/test_chokepoint.jl:50` (broaden classifier), then move the file's content into `test/test_types.jl`, delete `test/test_chokepoint.jl`, remove its `include` from `test/runtests.jl:` (the `include("test_chokepoint.jl")` line).
- Test: the moved testset is itself the test.

- [ ] **Step 1: Broaden the classifier in place**

In `test/test_chokepoint.jl`, replace line 50:

```julia
    return occursin(r"Parameter|::K[a-z]", arg_str)
```

with:

```julia
    return occursin(
        r"Parameter|::(Kd|Kiso|Kon|Koff|Kfor|Krev|Kreg|Keq|Etot|Lallo)\b",
        arg_str)
```

- [ ] **Step 2: Run to verify the broadened guard still passes**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: chokepoint testset PASSES (the `name(::Etot)`/`name(::Lallo)` bodies are now
recognized as renderer bodies; no false violation).

- [ ] **Step 3: Move the content into test_types.jl**

Append the entire body of `test/test_chokepoint.jl` (the helper functions
`_sig_fn_name`, `_sig_first_arg_str`, `_is_chokepoint_def`, `_symbol_call_pattern`,
`_walk_violations!`, the `const _CHOKEPOINT_PREFIX`, and the final `@testset`) to the end of
`test/test_types.jl`. Drop the `using Test` / `using EnzymeRates` lines (the host file
already has them). Keep the two `# ABOUTME:` lines as a section comment above the moved block.

- [ ] **Step 4: Delete the file and its include**

```bash
git rm test/test_chokepoint.jl
```

In `test/runtests.jl`, delete the line `    include("test_chokepoint.jl")`.

- [ ] **Step 5: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite PASSES; chokepoint assertions now run inside `test_types.jl`. Total test
count unchanged (modulo testset nesting).

- [ ] **Step 6: Commit**

```bash
git add test/test_types.jl test/runtests.jl
git commit -m "refactor: move+broaden chokepoint guard into test_types.jl (1:1 src<->test)"
```

---

### Task B7: Harden the `rate_equation` perf gate

**Files:**
- Modify: `test/test_rate_eq_derivation.jl:332-337` (the `test_rate_equation_performance` helper). The call site / assertions at `746-748` stay unchanged (`allocs, t = test_rate_equation_performance(...)`; `@test allocs == 0`; `@test t < 100e-9`).

- [ ] **Step 1: Inspect current code**

Current helper (`test/test_rate_eq_derivation.jl:332-337`):

```julia
function test_rate_equation_performance(m, params, concs)
    rate_equation(m, concs, params) # warmup/compile
    allocs = @allocated rate_equation(m, concs, params)
    t = @elapsed for _ in 1:10_000; rate_equation(m, concs, params); end
    return allocs, t / 10_000
end
```

It returns the **mean** (`t / 10_000`) over a loop whose result is discarded — so the
compiler may elide the calls, and the mean is GC/scheduling-inflated. Caller asserts
`allocs == 0` and `t < 100e-9` at lines 746-748.

- [ ] **Step 2: Replace mean-over-discard with sink + min-of-batches**

Replace the helper body (lines 333-336, keeping the function signature and the `return`
tuple shape `(allocs, t)` so the call site at 746-748 is unchanged):

```julia
function test_rate_equation_performance(m, params, concs)
    rate_equation(m, concs, params) # warmup/compile
    allocs = @allocated rate_equation(m, concs, params)
    # Minimum over several batches defeats GC/scheduling inflation a single
    # mean suffers; summing results into `acc` (then returning it via the
    # `isfinite` guard's side-effect-free observation) prevents the compiler
    # from eliding the calls as dead code.
    best = Inf
    acc = 0.0
    for _ in 1:5
        acc = 0.0
        t = @elapsed for _ in 1:10_000
            acc += rate_equation(m, concs, params)
        end
        best = min(best, t / 10_000)
    end
    isfinite(acc) || error("rate_equation produced non-finite result")
    return allocs, best
end
```

The caller (lines 746-748) is unchanged: `allocs, t = test_rate_equation_performance(m,
params, concs)`; `@test allocs == 0`; `@test t < 100e-9`. `t` is now the min-of-batches.

- [ ] **Step 3: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: the performance testset PASSES for every mechanism in `MECHANISM_TEST_SPECS`
(min-of-batches is ≤ the old mean, so anything that passed before still passes).

Optional elision sanity check — confirm `acc` is genuinely consumed (the `isfinite(acc) ||
error(...)` line forces the loop's result to be observed, so the optimizer cannot drop the
accumulation). No `@code_typed` needed; the `error` branch makes `acc` load-bearing.

- [ ] **Step 4: Commit**

```bash
git add test/test_rate_eq_derivation.jl
git commit -m "test: harden rate_equation perf gate (sink result, min-of-batches)"
```

---

### Task B8 + C3: Move canonical-hash code to identify_rate_equation.jl, test alongside

**Files:**
- Move: the canonical-rate-eq-hash block from `src/mechanism_enumeration.jl:1765-2009` to `src/identify_rate_equation.jl` (its sole consumer).
- Move: `test/test_canonical_hash_partition.jl` content into `test/test_identify_rate_equation.jl`; delete file; remove include.
- Update ABOUTME lines in both src files.

- [ ] **Step 1: Identify the exact block to move**

Run: `grep -n '_canonical_rate_eq_hash\|canonical' src/mechanism_enumeration.jl | head`
Confirm the contiguous block (functions feeding `_canonical_rate_eq_hash`) spans roughly
`1765-2009`. Read it fully before moving:
Run: `julia --project -e 'true'` (no-op; just read the file region with your editor/Read).

- [ ] **Step 2: Move the block**

Cut the canonical-hash functions from `src/mechanism_enumeration.jl` and paste them into
`src/identify_rate_equation.jl` (place near `_project_cached_params`, the consumer). Do not
change the code itself. If any moved function references a helper defined only in
`mechanism_enumeration.jl`, leave that helper where it is (resolution is at call time, both
files are included in the module) — only move the canonical-hash cluster.

- [ ] **Step 3: Update ABOUTME comments**

In `src/mechanism_enumeration.jl`, if the ABOUTME mentions canonical hashing, remove that
clause. In `src/identify_rate_equation.jl`, extend the ABOUTME to mention canonical
rate-equation hashing for fit reuse.

- [ ] **Step 4: Run to verify the move is behavior-neutral**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite PASSES (pure relocation; `test_canonical_hash_partition.jl` still
green at its current location).

- [ ] **Step 5: Move the test**

Append the body of `test/test_canonical_hash_partition.jl` (drop its `using` lines) to the
end of `test/test_identify_rate_equation.jl`, keeping the two ABOUTME lines as a section
comment.

```bash
git rm test/test_canonical_hash_partition.jl
```

In `test/runtests.jl`, delete the line `    include("test_canonical_hash_partition.jl")`.

- [ ] **Step 6: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite PASSES; partition golden (uni_uni→1, bi_bi→21) unchanged.

- [ ] **Step 7: Commit**

```bash
git add src/mechanism_enumeration.jl src/identify_rate_equation.jl test/test_identify_rate_equation.jl test/runtests.jl
git commit -m "refactor: relocate canonical-hash code+test to the identify layer (1:1 src<->test)"
```

---

### Task B9: Reuse allo-state validity consts in the singleton constructor

**Files:**
- Modify: `src/types.jl:417` (existing `_VALID_CAT_ALLO_STATES`), add `_VALID_REG_ALLO_STATES`; use both in the singleton `AllostericEnzymeMechanism` ctor (cat tags at `770-772`, reg tags at `787-789`) and the `RegulatorySite` ctor (`110-114`).
- Test: existing suite (validation tests already cover bad tags).

- [ ] **Step 1: Add the reg-state const**

Below the existing `const _VALID_CAT_ALLO_STATES = (:OnlyA, :EqualAI, :NonequalAI)`
(`src/types.jl:417`), add:

```julia
const _VALID_REG_ALLO_STATES = (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)
```

- [ ] **Step 2: Use it in `RegulatorySite` ctor**

In the `RegulatorySite` inner ctor (`src/types.jl:110-114`), replace:

```julia
        for st in allo_states
            st in (:OnlyA, :OnlyI, :EqualAI, :NonequalAI) ||
                error("RegulatorySite: allo state $st must be one of " *
                      ":OnlyA, :OnlyI, :EqualAI, :NonequalAI")
        end
```

with:

```julia
        for st in allo_states
            st in _VALID_REG_ALLO_STATES ||
                error("RegulatorySite: allo state $st must be one of " *
                      "$_VALID_REG_ALLO_STATES")
        end
```

> Note: `_VALID_REG_ALLO_STATES` is defined at line ~418 but `RegulatorySite` is at ~99.
> Consts are resolved at the `RegulatorySite` *call* time (runtime), not at definition parse
> time, so the forward reference is fine. Verify by running the suite in Step 5.

- [ ] **Step 3: Use the cat const in the singleton ctor**

In `AllostericEnzymeMechanism(cm, cat_sites, reg_sites)` (`src/types.jl:770-772`), replace:

```julia
        st in (:OnlyA, :EqualAI, :NonequalAI) ||
            error("Catalytic kinetic group $g has unknown allo state $st; " *
                  "must be one of (:OnlyA, :EqualAI, :NonequalAI)")
```

with:

```julia
        st in _VALID_CAT_ALLO_STATES ||
            error("Catalytic kinetic group $g has unknown allo state $st; " *
                  "must be one of $_VALID_CAT_ALLO_STATES")
```

Keep the separate `st === :OnlyI && error(...)` check above it unchanged (it gives a
specific message and is not redundant with the membership test).

- [ ] **Step 4: Use the reg const in the singleton ctor**

In the same function's reg-site loop (`src/types.jl:787-789`), replace:

```julia
            st in (:OnlyA, :OnlyI, :EqualAI, :NonequalAI) ||
                error("Reg site $i, ligand $(ligands[k]): unknown allo state $st")
```

with:

```julia
            st in _VALID_REG_ALLO_STATES ||
                error("Reg site $i, ligand $(ligands[k]): unknown allo state $st")
```

- [ ] **Step 5: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite PASSES; existing bad-tag rejection tests still throw (now via the
shared consts). The forward const reference resolves at runtime.

- [ ] **Step 6: Commit**

```bash
git add src/types.jl
git commit -m "refactor: share allo-state validity consts across all three validators"
```

---

### Task C2: Move dep-set-invariance test into test_rate_eq_derivation.jl

**Files:**
- Move: `test/test_dep_set_invariance.jl` content into `test/test_rate_eq_derivation.jl`; delete file; remove include from `test/runtests.jl`.

- [ ] **Step 1: Append the content**

Append the body of `test/test_dep_set_invariance.jl` (the `_dep_struct_key`,
`_dep_struct_key_set` helpers and the `@testset "dependent-param choice invariant to
group-rep"`) to the end of `test/test_rate_eq_derivation.jl`. Keep the three ABOUTME lines as
a section comment. The host file already runs inside the suite with `MECHANISM_TEST_SPECS`
in scope (both are included after `mechanism_definitions_for_test_enzyme_derivation.jl`), so
no new imports are needed — verify in Step 3.

- [ ] **Step 2: Delete file + include**

```bash
git rm test/test_dep_set_invariance.jl
```

In `test/runtests.jl`, delete the line `    include("test_dep_set_invariance.jl")`.

- [ ] **Step 3: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite PASSES; the invariance testset runs inside `test_rate_eq_derivation.jl`.

- [ ] **Step 4: Commit**

```bash
git add test/test_rate_eq_derivation.jl test/runtests.jl
git commit -m "refactor: move dep-set-invariance test into test_rate_eq_derivation.jl"
```

---

# CLUSTER D — Stale comment / doc sweep

### Task D1: Remove temporal labels from tests

**Files:**
- Modify: `test/test_dsl.jl:79`, `test/test_types.jl:2, 99, 296, 299, 834`

- [ ] **Step 1: Re-locate the exact strings (line numbers may have shifted after B/C edits)**

Run: `grep -rniE "new (design|concrete|grammar)" test/test_dsl.jl test/test_types.jl`
Expected: ~6 hits.

- [ ] **Step 2: Rewrite each to an evergreen description**

Apply these replacements (match on the testset string, not line number):
- `@testset "@enzyme_mechanism (new grammar)"` → `@testset "@enzyme_mechanism decomposed-Species grammar"`
- `@testset "EnzymeMechanism struct + accessors (new design)"` → `@testset "EnzymeMechanism struct + accessors"`
- `@testset "AllostericEnzymeMechanism (new design)"` → `@testset "AllostericEnzymeMechanism struct + accessors"`
- `@testset "EnzymeReaction (new concrete)"` → `@testset "EnzymeReaction struct + accessors"`
- For the two prose comments at `test_types.jl:296,299` ("were dropped — the new design accepts both", "design's terms; in the new design, the kinetic_group integer"): rewrite to describe current behavior without "new design", e.g. "Two reactions differing only in kinetic_group numbering compare equal; the kinetic_group integer is positional."

- [ ] **Step 3: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite PASSES (testset renames don't change behavior).

- [ ] **Step 4: Commit**

```bash
git add test/test_dsl.jl test/test_types.jl
git commit -m "docs: remove temporal labels from test names (evergreen descriptions)"
```

---

### Task D2: Fix stale source comments + CLAUDE.md

**Files:**
- Modify: `src/types.jl` (binding-canonicalization comments near `:129`, `:319`; indexed `k1f/k1r` doc near `:864` — **verify current line numbers**)
- Modify: `.claude/CLAUDE.md` (multiplicity "not enumerated" line; public-surface wording; guidance near `:245`)

- [ ] **Step 1: Re-locate the stale src comments**

Run: `grep -nE "k1f|k1r|k_1|index" src/types.jl`
Run: `grep -n "canonicaliz" src/types.jl | head`
Read each hit and confirm whether it describes current behavior. Only edit comments that are
**provably false** about the current code (per CLAUDE.md: never remove a comment unless it's
actively false). If the `:864` region (the `EnzymeMechanism{Metabolites,Reactions}` doc on
`main`) no longer exists on this branch, there is nothing to fix there — note it and move on.

- [ ] **Step 2: Fix the binding-canonicalization comments if stale**

For any comment that describes binding canonicalization in terms that no longer match the
current `Step` constructor (which canonicalizes the bound metabolite onto `to_species`),
rewrite to match. Example pattern — if a comment says "metabolite on the from side", correct
to "metabolite on the to side". Make the smallest change that makes the comment true.

- [ ] **Step 3: Update CLAUDE.md multiplicity line**

Run: `grep -n "not enumerated\|oligomeric_state\|allowed_catalytic_multiplicities" .claude/CLAUDE.md`
Replace the clause stating multiplicity is "not enumerated" with text reflecting Task A4:
`allowed_catalytic_multiplicities` IS enumerated by `_expand_to_allosteric` (one allosteric
variant set per allowed value); `oligomeric_state: N` is the single-value shorthand.

- [ ] **Step 4: Soften the public-surface wording**

Run: `grep -n "public mechanism-construction surface\|18 exported" .claude/CLAUDE.md`
Reword so `Mechanism` / `AllostericMechanism` / `init_mechanisms` are described as accessed
via `EnzymeRates.X` (internal but usable), NOT as exported public names. Leave the "18
exported public names" statement intact (it is correct).

- [ ] **Step 5: Review the `:245` guidance**

Read `.claude/CLAUDE.md` around the line about EqualAI dependent params / synth-dep. Only
edit if it's factually wrong about current code; otherwise leave it.

- [ ] **Step 6: Run to verify pass (docs don't affect tests, but confirm nothing broke)**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: full suite PASSES.

- [ ] **Step 7: Commit**

```bash
git add src/types.jl .claude/CLAUDE.md
git commit -m "docs: correct stale comments; CLAUDE.md multiplicity + public-surface wording"
```

---

# FINAL VERIFICATION

- [ ] **Full suite green**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: 27073+ pass, 1 broken (the deliberate Case-B item), 0 fail. New tests added in
A1-A5, B1, B4, B5 increase the pass count. Perf gate and chokepoint guard still green.

- [ ] **No Sig hand-decode remains**

Run: `grep -rn 'Sig\[' src/`
Expected: zero matches.

- [ ] **No stray parameter-symbol literals**

Run: `grep -rnE 'Symbol\("[KkVL]' src/`
Expected: zero matches (chokepoint intact).

- [ ] **Export count still 18**

Run: `grep -E '^export ' src/EnzymeRates.jl | sed 's/^export //' | tr ',' '\n' | sed 's/[[:space:]]//g' | grep -c .`
Expected: `18`.

- [ ] **Test files reduced to ~1:1**

Run: `ls test/test_chokepoint.jl test/test_dep_set_invariance.jl test/test_canonical_hash_partition.jl 2>&1`
Expected: all three "No such file" (moved into host files).

---

## Self-review notes (author)

- **Spec coverage:** A1/A2 (validation #1+#2), A3 (#2 regulator dup), A4 (#3 multiplicity),
  A5 (#4 README); B1 (metabolites bypass), B2 (iso-canon dup), B3 (emit-ladder), B4 (Sig
  metabolite collapse — the "free win" half of the Sig decision), B5 (DSL dup+wrong macro),
  B6 (chokepoint guard) folded into C1, B7 (perf gate), B8 (canonical-hash placement) folded
  into C3, B9 (allo-state consts); C1/C2/C3 (test reorg); D1/D2 (#5 stale comments + CLAUDE.md
  multiplicity + public surface). B10 (Wegscheider move) intentionally omitted per design.
  Sig full-rewrite intentionally omitted per design.
- **Known approximations:** D2 line numbers (`types.jl:864` etc.) come from the second
  agent's report and may be stale; each D2 step re-greps before editing and only touches
  provably-false comments.
- **Ordering rationale:** A first (validation surfaces fixture breakage early). B6→C1 and
  B8→C3 are paired (broaden/move code, then move test). C2 independent. D last.
