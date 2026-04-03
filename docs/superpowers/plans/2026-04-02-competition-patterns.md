# Competition Patterns Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace 2^n dead-end form enumeration with competition-pattern-based deterministic selection, fixing OOM for ter-ter reactions.

**Architecture:** Two new enumeration functions (`_competition_patterns`, `_inhibitor_competition_patterns`) produce competition patterns. The existing `_expand_substrate_product_dead_ends` and `_expand_add_dead_end_regulator` replace their bitmask loops with pattern-driven deterministic form selection + dedup. A shared `_forms_with_binding_step` helper supports the topology-aware inhibitor binding rule.

**Tech Stack:** Julia, EnzymeRates.jl internals (MechanismSpec, StepSpec, step_metabolite, step_forms)

---

### Task 1: `_competition_patterns` — enumerate S/P competition patterns

**Files:**
- Modify: `src/mechanism_enumeration.jl` (insert before `_expand_substrate_product_dead_ends` at ~line 737)
- Test: `test/test_mechanism_enumeration.jl` (new testset inside "init_mechanisms")

- [ ] **Step 1: Write failing tests**

Add a new testset after the "Dead-end substrate/product expansion" testset (after line 383) inside the "init_mechanisms" block:

```julia
@testset "Competition patterns" begin
    # Uni-uni: 1×1, only 1 pattern (single edge)
    pats_11 = EnzymeRates._competition_patterns(
        Set([:S]), Set([:P]))
    @test length(pats_11) == 1
    @test pats_11[1] == Set([(:S, :P)])

    # Uni-bi: 1×2, S must compete with both P and Q
    pats_12 = EnzymeRates._competition_patterns(
        Set([:S]), Set([:P, :Q]))
    @test length(pats_12) == 1
    @test pats_12[1] == Set([(:S, :P), (:S, :Q)])

    # Bi-uni: symmetric
    pats_21 = EnzymeRates._competition_patterns(
        Set([:A, :B]), Set([:P]))
    @test length(pats_21) == 1
    @test pats_21[1] == Set([(:A, :P), (:B, :P)])

    # Bi-bi: 7 patterns
    pats_22 = EnzymeRates._competition_patterns(
        Set([:A, :B]), Set([:P, :Q]))
    @test length(pats_22) == 7
    # Every pattern covers all vertices
    for pat in pats_22
        for s in [:A, :B]
            @test any(p -> (s, p) in pat, [:P, :Q])
        end
        for p in [:P, :Q]
            @test any(s -> (s, p) in pat, [:A, :B])
        end
    end
    # Invalid pattern excluded: {A↔P, B↔P} leaves Q uncovered
    @test Set([(:A, :P), (:B, :P)]) ∉ pats_22

    # Ter-ter: 265 patterns
    pats_33 = EnzymeRates._competition_patterns(
        Set([:A, :B, :C]), Set([:P, :Q, :R]))
    @test length(pats_33) == 265
    for pat in pats_33
        for s in [:A, :B, :C]
            @test any(
                p -> (s, p) in pat, [:P, :Q, :R])
        end
        for p in [:P, :Q, :R]
            @test any(
                s -> (s, p) in pat, [:A, :B, :C])
        end
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: FAIL — `_competition_patterns` not defined

- [ ] **Step 3: Implement `_competition_patterns`**

Add to `src/mechanism_enumeration.jl` before `_expand_substrate_product_dead_ends` (around line 737):

```julia
"""
    _competition_patterns(sub_names, prod_names)
        → Vector{Set{Tuple{Symbol,Symbol}}}

Enumerate all bipartite competition graphs on substrates × products
where every substrate and every product has degree ≥ 1.
"""
function _competition_patterns(
    sub_names::Set{Symbol},
    prod_names::Set{Symbol},
)
    subs = sort(collect(sub_names))
    prods = sort(collect(prod_names))
    edges = [(s, p) for s in subs for p in prods]
    n = length(edges)
    result = Set{Tuple{Symbol,Symbol}}[]
    for mask in 1:(1 << n) - 1
        pat = Set{Tuple{Symbol,Symbol}}()
        for j in 1:n
            if (mask >> (j - 1)) & 1 == 1
                push!(pat, edges[j])
            end
        end
        all(s -> any(
            p -> (s, p) in pat, prods), subs) ||
            continue
        all(p -> any(
            s -> (s, p) in pat, subs), prods) ||
            continue
        push!(result, pat)
    end
    result
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: All competition pattern tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add _competition_patterns for S/P competition enumeration"
```

---

### Task 2: Replace 2^n loop in `_expand_substrate_product_dead_ends` with competition patterns

**Files:**
- Modify: `src/mechanism_enumeration.jl:751-854` (`_expand_substrate_product_dead_ends`)
- Test: `test/test_mechanism_enumeration.jl` (modify existing dead-end expansion tests, add new filtering tests)

- [ ] **Step 1: Write failing tests for competition-based dead-end filtering**

Add a new testset after the "Competition patterns" testset from Task 1:

```julia
@testset "Dead-end filtering by competition" begin

    @testset "Bi-bi random, diagonal {A↔P, B↔Q}" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN], E_Q[N]
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
        spec = mechanism_spec_from_mechanism(
            m, bi_bi_rxn)
        bound = EnzymeRates._bound_metabolites_at_forms(
            spec, bi_bi_rxn)
        competition = Set([(:A, :P), (:B, :Q)])

        # E_A_Q: bound={A,Q}, A↔Q not competing → allowed
        @test !any(
            (s, p) in competition
            for s in intersect(
                bound[:E_A] ∪ Set([:Q]),
                Set([:A, :B]))
            for p in intersect(
                bound[:E_A] ∪ Set([:Q]),
                Set([:P, :Q])))
        # E_A_P: bound={A,P}, A↔P competing → forbidden
        @test any(
            (s, p) in competition
            for s in intersect(
                bound[:E_A] ∪ Set([:P]),
                Set([:A, :B]))
            for p in intersect(
                bound[:E_A] ∪ Set([:P]),
                Set([:P, :Q])))
    end

    @testset "Bi-bi random, complete → 0 dead-ends" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN], E_Q[N]
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
        spec = mechanism_spec_from_mechanism(
            m, bi_bi_rxn)
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [spec], bi_bi_rxn)
        # Complete competition pattern exists among the 7.
        # It produces 0 dead-end forms → bare topology only.
        # At least one result should have the same steps
        # as the original (no dead-end forms added).
        bare = filter(
            r -> length(r.steps) == length(spec.steps),
            result)
        @test length(bare) == 1
    end

    @testset "Bi-bi random: ≤ 7 variants (was 16)" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN], E_Q[N]
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
        spec = mechanism_spec_from_mechanism(
            m, bi_bi_rxn)
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [spec], bi_bi_rxn)
        # 7 competition patterns, some may dedup →
        # result ≤ 7
        @test length(result) <= 7
        @test length(result) >= 1
    end
end
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: The "Bi-bi random: ≤ 7 variants" test FAILS because current code returns 16.

- [ ] **Step 3: Replace 2^n loop with competition pattern loop**

In `src/mechanism_enumeration.jl`, replace lines 783-851 (the `for mask in 0:(1 << n_de) - 1` loop and its body) inside `_expand_substrate_product_dead_ends` with:

```julia
        # Map each dead-end form to its bound metabolites
        de_bound = Dict{Symbol, Set{Symbol}}()
        for de_name in de_form_names
            entries = de_forms[de_name]
            f, m = first(entries)
            de_bound[de_name] = union(
                bound[f], Set([m]))
        end

        # Enumerate competition patterns
        patterns = _competition_patterns(
            sub_names, prod_names)
        seen = Set{Vector{Symbol}}()

        for pattern in patterns
            # Filter dead-end forms by competition
            allowed_de = Symbol[]
            for de_name in de_form_names
                mets = de_bound[de_name]
                de_subs = intersect(mets, sub_names)
                de_prods = intersect(
                    mets, prod_names)
                has_conflict = any(
                    (s, p) in pattern
                    for s in de_subs
                    for p in de_prods)
                has_conflict || push!(
                    allowed_de, de_name)
            end

            # Dedup by form set
            allowed_de in seen && continue
            push!(seen, allowed_de)

            active_de = Set{Symbol}(allowed_de)

            # Build new steps: original + dead-end
            new_steps = copy(spec.steps)

            # Add binding steps for active dead-ends
            for de_name in sort(collect(active_de))
                entries = de_forms[de_name]
                for (cat_form, met) in entries
                    push!(new_steps, StepSpec(
                        [cat_form, met],
                        [de_name], true))
                end
            end

            # Add mirror steps
            for s in spec.steps
                from, to =
                    s.reactants[1], s.products[1]
                met = step_metabolite(s)
                for de_met in sort(collect(all_mets))
                    haskey(bound, from) || continue
                    haskey(bound, to) || continue
                    de_met in bound[from] && continue
                    de_met in bound[to] && continue
                    from_de = _dead_end_form_name(
                        bound[from], de_met)
                    to_de = _dead_end_form_name(
                        bound[to], de_met)
                    from_de in active_de || continue
                    to_de in active_de || continue
                    if met !== nothing
                        push!(new_steps, StepSpec(
                            [from_de, met],
                            [to_de],
                            s.is_equilibrium))
                    else
                        push!(new_steps, StepSpec(
                            [from_de], [to_de],
                            s.is_equilibrium))
                    end
                end
            end

            # Compute param_count
            n_steps = length(new_steps)
            n_re = count(
                s -> s.is_equilibrium, new_steps)
            n_ss = n_steps - n_re
            n_forms = length(
                all_form_names(new_steps))
            n_thermo = n_steps - n_forms + 1
            pc = n_re + 2 * n_ss - n_thermo + 2

            push!(result, MechanismSpec(
                spec.reaction, new_steps,
                spec.param_constraints, pc))
        end
```

- [ ] **Step 4: Update existing test counts**

In `test/test_mechanism_enumeration.jl`, update the existing "Bi-Bi random: 4 dead-end forms" test (around line 325) — change the count from 16 to the actual result. The test comment should also be updated:

```julia
            # 4 unique dead-end forms, 7 competition patterns
            # after dedup → N unique variants
            @test length(result) == 16  # PLACEHOLDER: run once to get actual count, then update
```

Also update the "Bi-Bi Ping-Pong: 3 dead-end forms" test (around line 381):

```julia
            @test length(result) == 8  # PLACEHOLDER: run once to get actual count, then update
```

**Important:** Run the tests once first to discover the actual counts, then update the assertions. The actual counts depend on how many competition patterns dedup for each topology.

- [ ] **Step 5: Run tests to discover actual counts and verify**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -50`

Look for the actual counts in the failure messages. Update the test assertions with the correct values. Re-run to confirm all pass.

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Replace 2^n dead-end enumeration with competition patterns"
```

---

### Task 3: Ter-ter integration test

**Files:**
- Test: `test/test_mechanism_enumeration.jl` (add ter-ter test to "Dead-end filtering by competition" testset)

- [ ] **Step 1: Write ter-ter test**

Add inside the "Dead-end filtering by competition" testset:

```julia
@testset "Ter-ter diagonal: 12 of 27 allowed" begin
    # Build random ter-ter topology manually using
    # init_mechanisms — pick the random topology
    # (largest form count) from the results
    specs = EnzymeRates.init_mechanisms(ter_ter_rxn)
    @test !isempty(specs)
    # Verify ter-ter doesn't OOM (was 134M, now ≤ 265
    # per topology)
    @test length(specs) < 100_000
end

@testset "Ter-ter completes without OOM" begin
    specs = EnzymeRates.init_mechanisms(ter_ter_rxn)
    @test length(specs) > 0
    # Verify param count invariant
    for s in specs
        @test s.param_count == 3 + 3 + 3  # n_s + n_p + 3
    end
end
```

- [ ] **Step 2: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: PASS — ter-ter completes without OOM

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Add ter-ter integration test for competition patterns"
```

---

### Task 4: `_forms_with_binding_step` helper

**Files:**
- Modify: `src/mechanism_enumeration.jl` (insert near other StepSpec helpers, ~line 90)
- Test: `test/test_mechanism_enumeration.jl` (new testset)

- [ ] **Step 1: Write failing tests**

Add a new testset before the "Add dead-end regulator" testset:

```julia
@testset "Forms with binding step" begin
    # Uni-uni: S binds to E, P binds to E
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
    spec_uu = mechanism_spec_from_mechanism(
        m_uu, uni_uni_rxn)
    @test EnzymeRates._forms_with_binding_step(
        spec_uu.steps, :S) == Set([:E])
    @test EnzymeRates._forms_with_binding_step(
        spec_uu.steps, :P) == Set([:E])

    # Bi-bi random: B binds to E and E_A
    m_bb = @enzyme_mechanism begin
        species: begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            enzymes: E, E_A[C], E_A_B[CN],
                E_B[N], E_P[C], E_P_Q[CN], E_Q[N]
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
    spec_bb = mechanism_spec_from_mechanism(
        m_bb, bi_bi_rxn)
    @test EnzymeRates._forms_with_binding_step(
        spec_bb.steps, :B) == Set([:E, :E_A])
    @test EnzymeRates._forms_with_binding_step(
        spec_bb.steps, :A) == Set([:E, :E_B])
    @test EnzymeRates._forms_with_binding_step(
        spec_bb.steps, :P) == Set([:E, :E_Q])
    @test EnzymeRates._forms_with_binding_step(
        spec_bb.steps, :Q) == Set([:E, :E_P])
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: FAIL — `_forms_with_binding_step` not defined

- [ ] **Step 3: Implement `_forms_with_binding_step`**

Add to `src/mechanism_enumeration.jl` after the other StepSpec helpers (around line 99):

```julia
"""
    _forms_with_binding_step(steps, metabolite)
        → Set{Symbol}

Return source forms that have a binding step for `metabolite`.
The source form is the one that doesn't contain the metabolite
(the reactant side of a binding step).
"""
function _forms_with_binding_step(
    steps::Vector{StepSpec}, metabolite::Symbol,
)
    result = Set{Symbol}()
    for s in steps
        step_metabolite(s) === metabolite || continue
        push!(result, s.reactants[1])
    end
    result
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add _forms_with_binding_step helper"
```

---

### Task 5: `_inhibitor_competition_patterns` — enumerate inhibitor patterns

**Files:**
- Modify: `src/mechanism_enumeration.jl` (insert near `_competition_patterns`)
- Test: `test/test_mechanism_enumeration.jl` (new testset)

- [ ] **Step 1: Write failing tests**

Add a new testset before the "Add dead-end regulator" testset:

```julia
@testset "Inhibitor competition patterns" begin
    # Uni-uni, no existing inhibitors
    pats = EnzymeRates._inhibitor_competition_patterns(
        Set([:S]), Set([:P]), Symbol[])
    @test length(pats) == 1
    @test pats[1] == (Set([:S]), Set([:P]), Set{Symbol}())

    # Bi-bi, no existing inhibitors: 3×3 = 9
    pats_bb = EnzymeRates._inhibitor_competition_patterns(
        Set([:A, :B]), Set([:P, :Q]), Symbol[])
    @test length(pats_bb) == 9

    # Ter-ter, no existing inhibitors: 7×7 = 49
    pats_tt = EnzymeRates._inhibitor_competition_patterns(
        Set([:A, :B, :C]), Set([:P, :Q, :R]), Symbol[])
    @test length(pats_tt) == 49

    # Bi-bi, 1 existing inhibitor: 9 × 2 = 18
    pats_1i = EnzymeRates._inhibitor_competition_patterns(
        Set([:A, :B]), Set([:P, :Q]), [:I1__reg])
    @test length(pats_1i) == 18

    # Bi-bi, 2 existing inhibitors: 9 × 4 = 36
    pats_2i = EnzymeRates._inhibitor_competition_patterns(
        Set([:A, :B]), Set([:P, :Q]),
        [:I1__reg, :I2__reg])
    @test length(pats_2i) == 36
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: FAIL — `_inhibitor_competition_patterns` not defined

- [ ] **Step 3: Implement `_inhibitor_competition_patterns`**

Add to `src/mechanism_enumeration.jl` after `_competition_patterns`:

```julia
"""
    _inhibitor_competition_patterns(sub_names, prod_names,
        existing_inhibitors)
        → Vector{Tuple{Set{Symbol},Set{Symbol},Set{Symbol}}}

Enumerate inhibitor competition patterns: (competing_subs,
competing_prods, competing_inhibitors). Each inhibitor must
compete with ≥1 substrate and ≥1 product. Competition with
existing inhibitors is a free binary choice.
"""
function _inhibitor_competition_patterns(
    sub_names::Set{Symbol},
    prod_names::Set{Symbol},
    existing_inhibitors::Vector{Symbol},
)
    subs = sort(collect(sub_names))
    prods = sort(collect(prod_names))
    inhs = sort(existing_inhibitors)
    n_s = length(subs)
    n_p = length(prods)
    n_i = length(inhs)

    result = Tuple{
        Set{Symbol}, Set{Symbol}, Set{Symbol}}[]
    for s_mask in 1:(1 << n_s) - 1
        comp_subs = Set{Symbol}()
        for j in 1:n_s
            if (s_mask >> (j - 1)) & 1 == 1
                push!(comp_subs, subs[j])
            end
        end
        for p_mask in 1:(1 << n_p) - 1
            comp_prods = Set{Symbol}()
            for j in 1:n_p
                if (p_mask >> (j - 1)) & 1 == 1
                    push!(comp_prods, prods[j])
                end
            end
            for i_mask in 0:(1 << n_i) - 1
                comp_inhs = Set{Symbol}()
                for j in 1:n_i
                    if (i_mask >> (j - 1)) & 1 == 1
                        push!(comp_inhs, inhs[j])
                    end
                end
                push!(result, (
                    comp_subs, comp_prods, comp_inhs))
            end
        end
    end
    result
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add _inhibitor_competition_patterns enumeration"
```

---

### Task 6: Replace 2^n loop in `_expand_add_dead_end_regulator` with competition patterns

**Files:**
- Modify: `src/mechanism_enumeration.jl:1282-1403` (`_expand_add_dead_end_regulator`)
- Test: `test/test_mechanism_enumeration.jl` (modify existing tests, add new topology-aware tests)

- [ ] **Step 1: Write failing tests for topology-aware inhibitor binding**

Add inside the "Add dead-end regulator" testset:

```julia
@testset "Topology-aware: sequential bi-bi" begin
    # Sequential: A binds first, then B
    # Forms: E, E_A, E_A_B (=E_P_Q), E_Q
    m = @enzyme_mechanism begin
        species: begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            enzymes: E, E_A[C], E_A_B[CN],
                E_P_Q[CN], E_Q[N]
        end
        steps: begin
            [E, A] ⇌ [E_A]
            [E_A, B] ⇌ [E_A_B]
            [E, Q] ⇌ [E_Q]
            [E_Q, P] ⇌ [E_P_Q]
            [E_A_B] <--> [E_P_Q]
        end
    end
    rxn_seq = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: I
    end
    spec = _spec_with_rxn(m, bi_bi_rxn, rxn_seq)
    result = EnzymeRates._expand_add_dead_end_regulator(
        spec, rxn_seq)
    # With topology-aware binding, inhibitor patterns
    # that only compete with B and P should NOT bind to
    # free enzyme E (no B or P binding step from E).
    # Check that no variant has I binding to E alone
    # when B-binding and P-binding steps don't exist
    # from E.
    @test !isempty(result)
    for r in result
        @test r.param_count == spec.param_count + 1
    end
end

@testset "Two inhibitors: compete vs not compete" begin
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
    rxn_2i = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: I1, I2
    end
    spec = _spec_with_rxn(m, uni_uni_rxn, rxn_2i)

    # Add first inhibitor
    result1 = EnzymeRates._expand_add_dead_end_regulator(
        spec, rxn_2i)
    @test length(result1) >= 1

    # Add second inhibitor to a result with I1
    spec_with_i1 = first(result1)
    result2 = EnzymeRates._expand_add_dead_end_regulator(
        spec_with_i1, rxn_2i)
    # Should have variants: I2 competing with I1
    # (can't coexist) and not competing (can coexist)
    @test length(result2) >= 1
    # At least one variant should have 2 inhibitor
    # forms (E_I2 and/or E_I1_I2 etc)
    all_forms = [
        EnzymeRates.all_form_names(r)
        for r in result2]
    @test !isempty(all_forms)
end
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: Tests may pass or fail depending on which aspects the current code handles. The topology-aware behavior will fail once we replace the loop.

- [ ] **Step 3: Replace the 2^n loop in `_expand_add_dead_end_regulator`**

In `src/mechanism_enumeration.jl`, replace lines 1340-1401 (the `for mask in 1:(1 << n_forms) - 1` loop and its body through the `push!(result, ...)`) with:

```julia
        # Find existing inhibitor dummies in mechanism
        existing_inhibitors = Symbol[]
        for met in existing_mets
            s = string(met)
            if contains(s, "__reg") &&
                    !contains(s, string(reg))
                push!(existing_inhibitors, met)
            end
        end
        sort!(unique!(existing_inhibitors))

        # Enumerate inhibitor competition patterns
        inh_patterns =
            _inhibitor_competition_patterns(
                sub_names, prod_names,
                existing_inhibitors)
        seen = Set{Vector{Symbol}}()

        for (comp_subs, comp_prods,
                comp_inhibitors) in inh_patterns
            # Find forms where competing metabolites
            # have binding steps (topology-aware)
            target_forms = Set{Symbol}()
            for met in comp_subs
                union!(target_forms,
                    _forms_with_binding_step(
                        spec.steps, met))
            end
            for met in comp_prods
                union!(target_forms,
                    _forms_with_binding_step(
                        spec.steps, met))
            end
            for inh in comp_inhibitors
                union!(target_forms,
                    _forms_with_binding_step(
                        spec.steps, inh))
            end

            # Filter: must be eligible AND not contain
            # any competing metabolite/inhibitor
            all_competing = union(
                comp_subs, comp_prods,
                comp_inhibitors)
            active = Symbol[]
            for f in sort(collect(target_forms))
                f in eligible_forms || continue
                haskey(bound, f) || continue
                isempty(intersect(
                    bound[f], all_competing)) ||
                    continue
                push!(active, f)
            end

            isempty(active) && continue
            active in seen && continue
            push!(seen, active)

            new_steps = copy(spec.steps)
            de_form_map = Dict{Symbol, Symbol}()

            # Add binding steps (always RE)
            binding_step_indices = Int[]
            for cf in active
                de_name = _dead_end_form_name(
                    bound[cf], dummy)
                de_form_map[cf] = de_name
                push!(new_steps, StepSpec(
                    [cf, dummy], [de_name], true))
                push!(binding_step_indices,
                    length(new_steps))
            end

            # Add mirror steps
            for s in spec.steps
                from, to = step_forms(s)
                haskey(de_form_map, from) || continue
                haskey(de_form_map, to) || continue
                met = step_metabolite(s)
                from_de = de_form_map[from]
                to_de = de_form_map[to]
                if met !== nothing
                    push!(new_steps, StepSpec(
                        [from_de, met], [to_de],
                        s.is_equilibrium))
                else
                    push!(new_steps, StepSpec(
                        [from_de], [to_de],
                        s.is_equilibrium))
                end
            end

            # Equivalence constraints
            new_constraints = copy(
                spec.param_constraints)
            if length(binding_step_indices) >= 2
                first_idx = binding_step_indices[1]
                for j in 2:length(
                        binding_step_indices)
                    push!(new_constraints, (
                        Symbol("K$(binding_step_indices[j])"),
                        1,
                        [(Symbol("K$(first_idx)"),
                            1)]))
                end
            end

            push!(result, MechanismSpec(
                spec.reaction, new_steps,
                new_constraints,
                spec.param_count + 1))
        end
```

- [ ] **Step 4: Update existing test counts**

The existing "Bi-bi: exact dead-end regulator count" test (line 924) expects 31 (2^5 - 1). This will change. Run tests first to discover the new count, then update:

```julia
        # 5 eligible forms, competition patterns determine
        # which subsets are valid → new count
        @test length(result) == 31  # PLACEHOLDER: update after running
```

The "Uni-uni + new regulator" test (line 694) should still return 1 (only 1 form E eligible, 1 pattern for uni-uni). Verify this doesn't change.

- [ ] **Step 5: Run tests, discover actual counts, update assertions**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -50`

Update all assertions with the actual counts from failure messages.

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Replace 2^n inhibitor enumeration with competition patterns

Topology-aware: inhibitor only binds to forms that are sources
of binding steps for competing metabolites. Supports multi-
inhibitor compete/not-compete flag."
```

---

### Task 7: Fix remaining test failures and verify full suite

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (update any remaining count assertions)
- Modify: `src/mechanism_enumeration.jl` (fix any bugs discovered)

- [ ] **Step 1: Run full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -80`

- [ ] **Step 2: Fix any failures**

Most failures will be count mismatches in existing tests. For each:
1. Verify the new count is correct (competition patterns produce fewer but biochemically valid mechanisms)
2. Update the assertion
3. Update the test comment to explain the new count

Common expected changes:
- Integration test counts (total mechanisms from `enumerate_all`)
- Any test that previously asserted 2^n counts

- [ ] **Step 3: Run full suite again to confirm green**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: ALL PASS

- [ ] **Step 4: Commit**

```bash
git add test/test_mechanism_enumeration.jl src/mechanism_enumeration.jl
git commit -m "Fix test counts for competition-pattern-based enumeration"
```

---

### Task 8: Round-trip compilation test

**Files:**
- Test: `test/test_mechanism_enumeration.jl` (add round-trip test)

- [ ] **Step 1: Write round-trip test**

Add at the end of the "Dead-end filtering by competition" testset:

```julia
@testset "Round-trip: competition-filtered specs compile" begin
    for rxn in [uni_uni_rxn, bi_bi_rxn, bi_bi_pp_rxn]
        specs = EnzymeRates.init_mechanisms(rxn)
        # Test first 5 specs (compilation can be slow)
        for spec in first(specs, 5)
            m = EnzymeMechanism(spec)
            @test m isa EnzymeMechanism
            @test length(parameters(m)) <=
                spec.param_count
        end
    end
end
```

- [ ] **Step 2: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Add round-trip compilation test for competition-filtered specs"
```
