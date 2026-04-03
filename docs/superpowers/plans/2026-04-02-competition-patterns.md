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

    # Shared bi-bi random mechanism for multiple tests
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

    @testset "Bi-bi random: 7 variants (was 16)" begin
        # 4 dead-end forms × 7 competition patterns.
        # Each pattern yields a distinct dead-end set:
        #   {A↔P,B↔Q}: {E_A_Q,E_B_P}
        #   {A↔Q,B↔P}: {E_A_P,E_B_Q}
        #   {A↔P,A↔Q,B↔P}: {E_B_Q}
        #   {A↔P,A↔Q,B↔Q}: {E_B_P}
        #   {A↔P,B↔P,B↔Q}: {E_A_Q}
        #   {A↔Q,B↔P,B↔Q}: {E_A_P}
        #   {A↔P,A↔Q,B↔P,B↔Q}: {} (no dead-ends)
        # All 7 sets are distinct → 7 variants after dedup
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [spec_bb], bi_bi_rxn)
        @test length(result) == 7
    end

    @testset "Bi-bi random: complete competition → bare topology" begin
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [spec_bb], bi_bi_rxn)
        # Complete pattern {A↔P,A↔Q,B↔P,B↔Q} forbids
        # all dead-end forms → 1 variant has no dead-end
        # steps (same step count as original)
        bare = filter(
            r -> length(r.steps) == length(spec_bb.steps),
            result)
        @test length(bare) == 1
    end

    @testset "Bi-bi random: diagonal has exactly 2 dead-end forms" begin
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [spec_bb], bi_bi_rxn)
        # Diagonal patterns {A↔P,B↔Q} and {A↔Q,B↔P}
        # each allow exactly 2 dead-end forms. Find the
        # variant with 2 dead-end binding steps.
        # Original has 9 steps. Each dead-end form adds
        # 2 binding steps (from 2 catalytic parents) +
        # mirror steps. With 2 forms: 9 + 4 binding +
        # mirror steps.
        two_de = filter(result) do r
            de_forms = setdiff(
                EnzymeRates.all_form_names(r),
                EnzymeRates.all_form_names(spec_bb))
            length(de_forms) == 2
        end
        @test length(two_de) == 2  # diagonal + anti-diagonal
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

In `test/test_mechanism_enumeration.jl`, update the existing "Bi-Bi random: 4 dead-end forms" test (around line 325). Change the count and comment:

```julia
            # 4 unique dead-end forms, 7 competition patterns,
            # all 7 produce distinct dead-end sets → 7 variants
            @test length(result) == 7
```

Also update the "Bi-Bi Ping-Pong: 3 dead-end forms" test (around line 381). Ping-pong has 3 dead-end forms (E_A_P, E_A_Q, E_B_Q) but E_B_P doesn't exist, so B↔P edges have no effect, causing patterns 4&7 and 1&5 to collapse:

```julia
            # 3 dead-end forms, 7 competition patterns,
            # 5 unique dead-end sets after dedup → 5 variants
            @test length(result) == 5
```

- [ ] **Step 5: Run tests to verify all pass**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: All tests PASS

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
@testset "Ter-ter completes without OOM" begin
    # Previously OOM: 27 dead-end forms → 2^27 = 134M
    # Now: 265 competition patterns per topology, each
    # deterministic, with dedup
    specs = EnzymeRates.init_mechanisms(ter_ter_rxn)
    @test length(specs) > 0
    @test length(specs) < 100_000
    # Param count invariant: n_s + n_p + 3 = 9
    for s in specs
        @test s.param_count == 9
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
    # Sequential: A first, B second, P released, Q last
    # Forms: E, E_A, E_A_B (=E_P_Q), E_Q
    # Binding step sources:
    #   A: from E    B: from E_A
    #   P: from E_Q  Q: from E
    # Eligible forms: E, E_A, E_Q
    # 9 inhibitor patterns → 4 unique form sets:
    #   {E,E_Q}: ({A},{P}), ({A,B},{P})
    #   {E,E_A}: ({B},{Q}), ({B},{P,Q})
    #   {E_A,E_Q}: ({B},{P})
    #   {E}: ({A},{Q}), ({A},{P,Q}), ({A,B},{Q}),
    #         ({A,B},{P,Q})
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
    @test length(result) == 4
    for r in result
        @test r.param_count == spec.param_count + 1
    end
    # Verify E is NOT a target when I competes only
    # with B and P: {B,P} binding sources = {E_A,E_Q},
    # so the variant binding {E_A,E_Q} must exist
    form_sets = [sort(collect(setdiff(
        EnzymeRates.all_form_names(r),
        EnzymeRates.all_form_names(spec))))
        for r in result]
    # At least one variant has dead-end forms from E_A
    # and E_Q but NOT E
    @test any(fs -> all(
        f -> !startswith(string(f), "E_I") ||
            !contains(string(f), "_A_") &&
            !contains(string(f), "_Q_"),
        fs) == false, form_sets)
end

@testset "Bi-bi random: 9 inhibitor variants" begin
    # Random bi-bi with inhibitor: 9 patterns, each
    # produces a distinct form set → 9 variants
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
    bb_with_reg = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: I
    end
    spec = _spec_with_rxn(m, bi_bi_rxn, bb_with_reg)
    result = EnzymeRates._expand_add_dead_end_regulator(
        spec, bb_with_reg)
    # 9 inhibitor patterns (3 sub subsets × 3 prod
    # subsets), all distinct form sets → 9 variants
    @test length(result) == 9
end

@testset "Two inhibitors: compete vs not compete" begin
    # Use bi-bi random so I1 has mirror steps.
    # I1 competes with ({A},{P}) → binds to {E,E_B,E_Q}.
    # Mirror steps: E_I1+B→E_B_I1, E_I1+Q→E_Q_I1.
    # When adding I2: compete-with-I1 excludes I1 forms,
    # not-compete-with-I1 can use I1 forms via mirrors.
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
    rxn_2i = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: I1, I2
    end
    spec = _spec_with_rxn(m, bi_bi_rxn, rxn_2i)
    # Add I1 first
    result1 = EnzymeRates._expand_add_dead_end_regulator(
        spec, rxn_2i)
    @test length(result1) == 9
    # Pick a variant where I1 binds to multiple forms
    # (so mirror steps exist)
    multi_form = filter(result1) do r
        inh_forms = filter(
            f -> contains(string(f), "I1__reg"),
            collect(EnzymeRates.all_form_names(r)))
        length(inh_forms) >= 2
    end
    @test !isempty(multi_form)
    spec_i1 = first(multi_form)
    # Add I2
    result2 = EnzymeRates._expand_add_dead_end_regulator(
        spec_i1, rxn_2i)
    # With 1 existing inhibitor: 9 × 2 = 18 patterns.
    # After dedup some collapse → verify > 9 (the extra
    # variants come from compete/not-compete with I1)
    @test length(result2) >= 9
    # Verify compete vs not-compete produces different
    # form sets: not-competing variants can have forms
    # containing I1__reg, competing variants cannot
    has_i1_coexist = any(result2) do r
        any(f -> contains(string(f), "I1__reg") &&
                 contains(string(f), "I2__reg"),
            collect(EnzymeRates.all_form_names(r)))
    end
    @test has_i1_coexist  # at least one non-competing variant
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

The existing "Bi-bi: exact dead-end regulator count" test (line 924) expects 31 (2^5 - 1). With competition patterns, bi-bi random has 9 inhibitor patterns, all producing distinct form sets → 9 variants:

```julia
        # 5 eligible forms, 9 inhibitor competition patterns
        # (3 sub subsets × 3 prod subsets), all distinct → 9
        @test length(result) == 9
```

The "Uni-uni + new regulator" test (line 694) stays 1 (only 1 form E eligible, 1 pattern for uni-uni).

- [ ] **Step 5: Run tests to verify all pass**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30`
Expected: All tests PASS

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
