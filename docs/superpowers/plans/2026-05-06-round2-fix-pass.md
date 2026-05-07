# Round 2 Fix Pass — Implementation Plan

**Goal:** Address the highest-leverage findings from round 2 of the competing-subagent review (~32 findings combined). Skip stylistic-only items; focus on structural patterns and missing direct coverage.

**Architecture:** Test-only changes. ~7 commits to `test/test_mechanism_enumeration.jl`. No `src/` changes planned; surface bugs per protocol if any.

**Reference:** Round 2 review findings (in conversation history). Original plan: `docs/superpowers/plans/2026-05-06-fix-test-review-findings.md`.

## Common verification

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -3
```

Baseline: 26255. Each task may add or remove tests; final expected ~26200-26350.

If a NEW assertion fails: surface as DONE_WITH_CONCERNS or BLOCKED. Per round 1's experience, plan-stated counts can be wrong — re-derive from source.

---

## Task 1: Remove tautological `result isa Vector{...}` checks from Task 1's fix

**Findings:** R2 #1 (Round 2). The `result isa Vector{AllostericMechanismSpec}` checks added in commit `9e860c3` are tautological — Julia's signature guarantees the return type.

**Files:** `test/test_mechanism_enumeration.jl`

- [ ] **Step 1.1:** Locate the three `@test result isa Vector{AllostericMechanismSpec}` occurrences. They were added in commit `9e860c3` to the negative-fallback testsets:
  - `_expand_to_allosteric` AllostericMechanismSpec → empty
  - `_expand_change_allo_state` MechanismSpec → empty
  - `_expand_add_allosteric_regulator` MechanismSpec → empty

```bash
grep -n "result isa Vector{AllostericMechanismSpec}" test/test_mechanism_enumeration.jl
```

- [ ] **Step 1.2:** Replace each with a behavioral check. The simplest meaningful behavioral check: feed the empty result through `_push_to_dict!` (which `expand_mechanisms` does) and verify it doesn't add anything to the dict. OR: assert `eltype(result) === AllostericMechanismSpec`. The eltype check is also signature-guaranteed BUT is at least more honest about what's being tested.

  **Decision:** delete the `isa` lines entirely. The `isempty(result)` already covers the meaningful check. Document with a one-line comment that the type is signature-guaranteed.

- [ ] **Step 1.3:** Run tests, verify green, commit.

```
test: drop tautological `isa Vector{T}` checks from negative-fallback tests

The Round 1 fix at 9e860c3 added type-check assertions to three
negative-fallback tests. Round 2 review surfaced these as tautological:
the source's fallback methods return concrete-typed empty literals
(`AllostericMechanismSpec[]`), so Julia infers the return type by
signature alone. The `isa Vector{AllostericMechanismSpec}` assertion
cannot fail without breaking Julia's type system.

Drop the three lines. The `isempty(result)` assertion covers the
meaningful check; the type contract is enforced by the function
signature itself.
```

---

## Task 2: Replace vacuous "inter-move overlap" test; tighten "different mechanisms preserved"

**Findings:** R1 #3, R2 #10 (inter-move overlap is vacuous — uni-uni doesn't have inter-move overlap); R1 #2, R2 #11 (different-mechanisms-preserved is bound-only).

**Files:** `test/test_mechanism_enumeration.jl`

- [ ] **Step 2.1:** Locate the two testsets:

```bash
grep -n "@testset \"Inter-move overlap collapses via dedup\!\"\|@testset \"Different mechanisms preserved\"" test/test_mechanism_enumeration.jl
```

- [ ] **Step 2.2:** Rewrite "Inter-move overlap" with a real bi-bi seed where two different moves CAN produce the same target spec.

The simplest inter-move overlap: take a bi-bi init spec, apply RE→SS on group 1 → spec_A. Apply split + RE→SS on group 1 (split to fresh group, then convert original group to SS) → spec_B. If spec_A and spec_B compile to the same EnzymeMechanism (which would happen if split-then-RE→SS produces a structural equivalent of just RE→SS for that group), dedup should collapse them.

Actually, simpler: just verify the **count** behavior. Run `expand_mechanisms` and `dedup!`, save the pre-dedup and post-dedup counts, assert that for AT LEAST ONE param-count bucket, post-dedup < pre-dedup. This proves dedup actually fired.

```julia
@testset "Inter-move overlap: dedup actually fires" begin
    # Run expand_mechanisms on a bi-bi init seed, then dedup!. Assert
    # that at least one param-count bucket has post-dedup count < pre-dedup
    # count — proving dedup actually collapsed something.
    init_specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
    expanded = EnzymeRates.expand_mechanisms(init_specs, bi_bi_rxn)
    pre_dedup_counts = Dict(pc => length(specs) for (pc, specs) in expanded)
    EnzymeRates.dedup!(expanded)
    post_dedup_counts = Dict(pc => length(specs) for (pc, specs) in expanded)
    # At least one bucket should have shrunk.
    @test any(pre_dedup_counts[pc] > post_dedup_counts[pc]
              for pc in keys(post_dedup_counts))
    # Within each surviving bucket, all specs compile to distinct mechanisms.
    for (_pc, specs) in expanded
        compiled = Set(EnzymeRates.compile_mechanism(s) for s in specs)
        @test length(compiled) == length(specs)
    end
end
```

- [ ] **Step 2.3:** Tighten "Different mechanisms preserved" to compare unique compiled mechanisms count.

Locate the testset (around line 3863). Replace the bound-only assertion:

```julia
# Before:
@test length(cache[pc]) >= 1
@test length(cache[pc]) <= length(specs)

# After:
unique_compiled = Set(EnzymeRates.compile_mechanism(s) for s in specs)
EnzymeRates.dedup!(cache)
post_dedup_compiled = Set(EnzymeRates.compile_mechanism(s) for s in cache[pc])
# Dedup should preserve every distinct compiled mechanism that was input.
@test post_dedup_compiled == unique_compiled
# And the surviving spec count equals the unique-compiled count
# (no duplicates left, but no over-collapse either).
@test length(cache[pc]) == length(post_dedup_compiled)
```

- [ ] **Step 2.4:** Run tests, verify green, commit.

```
test(dedup): replace vacuous inter-move overlap test; tighten preservation

Round 2 review found two vacuous tests:
- "Inter-move overlap collapses via dedup!" used a uni-uni seed where
  no inter-move overlap exists, so the post-dedup `length(compiled) ==
  length(specs)` invariant was tautological by construction. Replace
  with a bi-bi seed and assert that AT LEAST ONE bucket shrinks
  post-dedup (proving dedup actually fires) AND that the surviving
  specs compile to distinct mechanisms.
- "Different mechanisms preserved" had bound-only assertions
  (`>= 1 && <= length(specs)`). Tighten to verify dedup preserves
  every distinct compiled mechanism AND collapses to exactly the
  unique-compiled count.
```

---

## Task 3: Add exact counts to asymmetric `_competition_patterns` testsets

**Findings:** R1 #1. The 2×3 and 3×2 cases lack count assertions.

**Files:** `test/test_mechanism_enumeration.jl`

- [ ] **Step 3.1:** Determine the actual counts.

`_competition_patterns(subs, prods)` enumerates bipartite covers — every substrate paired with at least one product, every product paired with at least one substrate, edges chosen subject to that. Known counts:
- 1×1 = 1, 1×2 = 1, 2×1 = 1, 2×2 = 7, 3×3 = 265

For 2×3: the formula is the number of bipartite covers on K(2,3). By inclusion-exclusion or direct enumeration. 

The simplest empirical approach: run the function and observe. Per the independent-derivation rule, this IS reading the code's output, but the function is itself a closed-form enumeration — observing the output and verifying it doesn't change is acceptable IF the count is small enough to spot-check by hand. For 2×3, the count is small.

Alternative: cite the known result for bipartite covers and assert that. The number of bipartite covers on K(m,n) is given by:

```
sum((-1)^i * C(m, i) * (2^(n-i) - 1)^... )
```

Actually the formula for "number of bipartite graphs on K(m,n) with no isolated vertex" is non-trivial. Just compute and verify by reading the source — `_competition_patterns(2 subs, 3 prods)` enumerates pairs; the count by symmetry equals 3×2 = same.

**Pragmatic path:** spot-check by running the function once during plan-write. The implementer should:
1. Read `_competition_patterns` source.
2. Run `_competition_patterns(Set([:A, :B]), Set([:P, :Q, :R]))` once, count.
3. Independently sanity-check by listing the patterns by hand (this is feasible for 2×3).
4. Add `@test length(pats) == <derived value>`.

If derivation matches running output, commit. If not, surface as bug.

- [ ] **Step 3.2:** Apply, run, commit.

```
test: add exact-count assertions to asymmetric _competition_patterns

Round 1 fix added 2×3 and 3×2 testsets with only `!isempty + per-vertex
coverage` checks. Tighten to exact counts derived by spot-check
enumeration of bipartite covers on K(2,3) and K(3,2). By symmetry,
both should produce the same count.
```

---

## Task 4: Add direct unit tests for `_canonicalize!` and `_dedup_key`

**Findings:** R1 #4, R2 #3, R2 #4 (non-contiguous IDs). Both helpers have no direct unit tests.

**Files:** `test/test_mechanism_enumeration.jl`

- [ ] **Step 4.1:** Add a `@testset "_canonicalize!"` block in Section 5 (composition).

```julia
@testset "_canonicalize!" begin
    # Test 1: kinetic-group renumbering preserves equivalence.
    # Build two specs differing only by kinetic_group integer values; after
    # _canonicalize!, both should have the same renumbering.
    spec_a = MechanismSpec(uni_uni_rxn,
        [StepSpec([:E, :S], [:E_S], true, 5),
         StepSpec([:E, :P], [:E_P], true, 7),
         StepSpec([:E_S], [:E_P], false, 12)],
        3)
    spec_b = MechanismSpec(uni_uni_rxn,
        [StepSpec([:E, :S], [:E_S], true, 1),
         StepSpec([:E, :P], [:E_P], true, 2),
         StepSpec([:E_S], [:E_P], false, 3)],
        3)
    EnzymeRates._canonicalize!(spec_a)
    EnzymeRates._canonicalize!(spec_b)
    @test EnzymeRates._dedup_key(spec_a) == EnzymeRates._dedup_key(spec_b)
    # Both should have IDs {1, 2, 3} after canonicalization
    @test Set(s.kinetic_group for s in spec_a.steps) == Set([1, 2, 3])

    # Test 2: AllostericMechanismSpec — site permutation with DISTINCT
    # multiplicities. The canonicalizer must permute multiplicities
    # alongside reg_sites; an off-by-one would collapse spec_c and spec_d
    # incorrectly.
    base = MechanismSpec(uni_uni_rxn,
        [StepSpec([:E, :S], [:E_S], true, 1),
         StepSpec([:E, :P], [:E_P], true, 2),
         StepSpec([:E_S], [:E_P], false, 3)],
        3)
    spec_c = AllostericMechanismSpec(
        base, 2,
        Vector{Symbol}[[:R1], [:R2]], Int[2, 4],
        Dict(1 => :EqualRT, 2 => :EqualRT, 3 => :EqualRT),
        Dict(:R1 => :OnlyR, :R2 => :OnlyT),
        4)
    spec_d = AllostericMechanismSpec(
        base, 2,
        Vector{Symbol}[[:R2], [:R1]], Int[4, 2],   # sites swapped, multiplicities follow
        Dict(1 => :EqualRT, 2 => :EqualRT, 3 => :EqualRT),
        Dict(:R1 => :OnlyR, :R2 => :OnlyT),
        4)
    EnzymeRates._canonicalize!(spec_c)
    EnzymeRates._canonicalize!(spec_d)
    @test EnzymeRates._dedup_key(spec_c) == EnzymeRates._dedup_key(spec_d)
    # After canonicalization, sites should be sorted alphabetically;
    # multiplicities should follow.
    @test spec_c.allosteric_reg_sites == [[:R1], [:R2]]
    @test spec_c.allosteric_multiplicities == [2, 4]
end
```

- [ ] **Step 4.2:** Add a `@testset "_dedup_key"` block in Section 5.

```julia
@testset "_dedup_key" begin
    # Same content → same key.
    base = MechanismSpec(uni_uni_rxn,
        [StepSpec([:E, :S], [:E_S], true, 1),
         StepSpec([:E, :P], [:E_P], true, 2),
         StepSpec([:E_S], [:E_P], false, 3)],
        3)
    @test EnzymeRates._dedup_key(base) == EnzymeRates._dedup_key(base)

    # AllostericMechanismSpec: differing multiplicities → different keys.
    spec1 = AllostericMechanismSpec(
        base, 2,
        Vector{Symbol}[[:R]], Int[2],
        Dict(1 => :EqualRT, 2 => :EqualRT, 3 => :EqualRT),
        Dict(:R => :OnlyR),
        4)
    spec2 = AllostericMechanismSpec(
        base, 2,
        Vector{Symbol}[[:R]], Int[4],   # different multiplicity
        Dict(1 => :EqualRT, 2 => :EqualRT, 3 => :EqualRT),
        Dict(:R => :OnlyR),
        4)
    @test EnzymeRates._dedup_key(spec1) != EnzymeRates._dedup_key(spec2)
end
```

- [ ] **Step 4.3:** Run tests, verify green, commit.

```
test(composition): add direct unit tests for _canonicalize! and _dedup_key

Both helpers were tested only indirectly via dedup! — a bug in either
would only surface through dedup's observable behavior. Adds:
- _canonicalize! kinetic-group renumbering: build two specs with
  different ID sets ({5,7,12} vs {1,2,3}); both canonicalize to the
  same dedup key.
- _canonicalize! AllostericMechanismSpec site permutation with DISTINCT
  multiplicities: an off-by-one in the permutation would mismatch
  multiplicities to sites; verified the canonical (sorted) order has
  the right multiplicity per site.
- _dedup_key sensitivity: two specs differing only by allosteric
  multiplicities produce different keys.
```

---

## Task 5: Add missing 7-checklist items to fix-pass tests

**Findings:** R2 #9 (bi-bi PP RE→SS missing item 4), R2 #15 (SS :NonequalRT split missing item 5).

**Files:** `test/test_mechanism_enumeration.jl`

- [ ] **Step 5.1:** Locate the bi-bi ping-pong RE→SS testset (added in commit `b345692`):

```bash
grep -n "MechanismSpec — bi-bi ping-pong: 5 RE groups" test/test_mechanism_enumeration.jl
```

Add an item-4 property check between items 3 and 5:

```julia
# 4. property-style: each variant flips exactly one RE step from
# is_equilibrium=true to false. The flipped step's kinetic_group
# differs across the 5 variants (one per RE group).
flipped_groups = Int[]
for r in result
    flipped = [s_new.kinetic_group
        for (s_old, s_new) in zip(spec.steps, r.steps)
        if s_old.is_equilibrium && !s_new.is_equilibrium]
    @test length(flipped) == 1
    push!(flipped_groups, only(flipped))
end
@test length(unique(flipped_groups)) == 5
```

- [ ] **Step 5.2:** Locate the SS :NonequalRT split testset (added in commit `c5a709c`):

```bash
grep -n "AllostericMechanismSpec — SS multi-step :NonequalRT split" test/test_mechanism_enumeration.jl
```

Add the missing item-5 preservation block (matching the EqualRT split sibling):

```julia
# 5. preservation
for r in result
    @test r.catalytic_n == spec.catalytic_n
    @test r.allosteric_reg_sites == spec.allosteric_reg_sites
    @test r.allosteric_multiplicities == spec.allosteric_multiplicities
    @test r.reg_ligand_tags == spec.reg_ligand_tags
    @test r.base.reaction === spec.base.reaction
end
```

- [ ] **Step 5.3:** Run tests, verify green, commit.

```
test: add missing 7-checklist items to fix-pass testsets

- bi-bi ping-pong RE→SS (commit b345692) had items 1, 2, 3, 5 — add
  item 4: each variant flips exactly one RE step; flipped kinetic_groups
  cover 5 distinct values across the 5 variants.
- SS :NonequalRT split (commit c5a709c) had no item 5 preservation
  block. Add the standard one (catalytic_n, reg_sites, multiplicities,
  reg_ligand_tags, base.reaction) matching the sibling :EqualRT
  split's preservation.
```

---

## Task 6: Add edge case tests

**Findings:** R2 #6 (empty inputs), R2 #4 (non-contiguous kinetic_group dedup), R2 #8 (substrate-as-product overlap), R1 #10 (Estar dead-end forms in comments but never asserted).

**Files:** `test/test_mechanism_enumeration.jl`

- [ ] **Step 6.1:** Add empty-input tests for `expand_mechanisms` and `dedup!`. Place inside the existing `expand_mechanisms` and `Dedup` testsets respectively:

```julia
@testset "Empty input" begin
    # expand_mechanisms with empty input returns empty Dict.
    @test isempty(EnzymeRates.expand_mechanisms(
        EnzymeRates.AbstractMechanismSpec[], uni_uni_rxn))
    # dedup! on empty cache is a no-op (and cleans up empty buckets).
    cache = Dict{Int, Vector{EnzymeRates.AbstractMechanismSpec}}()
    EnzymeRates.dedup!(cache)
    @test isempty(cache)
end
```

- [ ] **Step 6.2:** Add a non-contiguous-ID dedup test. Place inside the existing `Dedup` testset:

```julia
@testset "Non-contiguous kinetic_group IDs collapse via canonicalization" begin
    # Two specs with the same step shape but different kinetic_group
    # integer values (different subsets of the integers) should collapse
    # to one after dedup!. This exercises _canonicalize!'s renumbering
    # path with non-contiguous input IDs.
    spec_a = MechanismSpec(uni_uni_rxn,
        [StepSpec([:E, :S], [:E_S], true, 5),
         StepSpec([:E, :P], [:E_P], true, 9),
         StepSpec([:E_S], [:E_P], false, 14)],
        3)
    spec_b = MechanismSpec(uni_uni_rxn,
        [StepSpec([:E, :S], [:E_S], true, 2),
         StepSpec([:E, :P], [:E_P], true, 6),
         StepSpec([:E_S], [:E_P], false, 11)],
        3)
    cache = Dict(3 => EnzymeRates.AbstractMechanismSpec[spec_a, spec_b])
    EnzymeRates.dedup!(cache)
    @test length(cache[3]) == 1
end
```

- [ ] **Step 6.3:** Add a substrate-as-product overlap test (racemase / self-isomerase shape). Place in section 2 init testset or as a standalone testset:

```julia
@testset "Substrate-as-product overlap (racemase shape)" begin
    # A reaction where the same metabolite name might appear as both a
    # substrate and a product is rejected at the EnzymeReaction level
    # (no metabolite can be in both lists by name). But a substrate
    # racemase has DIFFERENT name (e.g., L-Ala → D-Ala). Test that
    # such a reaction's init mechanisms compile correctly.
    rxn = @enzyme_reaction begin
        substrates: L_Ala[CHN]
        products: D_Ala[CHN]
    end
    specs = EnzymeRates.init_mechanisms(rxn)
    @test !isempty(specs)
    for spec in first(specs, 3)
        m = EnzymeMechanism(spec)
        @test m isa EnzymeMechanism
    end
end
```

(If the prompt's intent was a metabolite literally appearing in both lists, that's caught at the EnzymeReaction constructor — already tested in Task 2's error-path tests. The racemase shape above tests the closely-related case that DOES need to work.)

- [ ] **Step 6.4:** Add an Estar-dead-end-form assertion to "Bi-Bi PP: 5 dead-end forms → 7 variants". Locate around line 1098. Add inside the existing testset:

```julia
# Assert that some result variants contain Estar-prefixed dead-end forms
# (proving the Estar branch of _dead_end_form_name is reached).
seed_forms = EnzymeRates.all_form_names(spec)
new_estar_forms = Set{Symbol}()
for r in result
    new_forms = setdiff(EnzymeRates.all_form_names(r), seed_forms)
    for f in new_forms
        startswith(string(f), "Estar_") && push!(new_estar_forms, f)
    end
end
@test !isempty(new_estar_forms)   # At least one Estar dead-end form was emitted.
```

- [ ] **Step 6.5:** Run tests, verify green, commit.

```
test: add edge-case coverage (empty inputs, non-contiguous IDs, Estar)

- Empty inputs to expand_mechanisms and dedup!.
- Non-contiguous kinetic_group IDs ({5,9,14} vs {2,6,11}) round-trip
  via dedup!'s canonicalization to a single equivalent spec.
- Substrate-as-product overlap (racemase L-Ala → D-Ala): init_mechanisms
  produces compileable specs.
- Estar dead-end forms in "Bi-Bi PP: 5 dead-end forms → 7 variants":
  assert at least one result variant contains an Estar-prefixed
  dead-end form (proving the Estar branch of _dead_end_form_name fires).
```

---

## Task 7: Floor activation and helper edge-case tests

**Findings:** R2 #7 (floor activation never proven), R1 #6 (existing_de exclusion branch never tested), R1 #7 (new-site multiplicity = catalytic_n contract).

**Files:** `test/test_mechanism_enumeration.jl`

- [ ] **Step 7.1:** Add `_apply_equivalence_grouping` floor-activation test.

The floor `pc = max(formula, n_subs + n_prods + 1)` fires when mirror cycles cause the formula to underestimate. Find a specific spec where formula < floor and verify pc == floor.

The simplest case: bi-bi with a dead-end inhibitor at multiple forms creates mirror cycles. After `_apply_equivalence_grouping`, the formula counts those mirror-cycle constraints incorrectly.

```julia
@testset "_apply_equivalence_grouping floor activates for mirror cycles" begin
    # Per CLAUDE.md and src/mechanism_enumeration.jl:1418-1438, the floor
    # fires when mirror cycles cause the formula to underestimate.
    # init_mechanisms applies _apply_equivalence_grouping; we verify
    # at least one bi-bi-with-inhibitor init spec has the formula
    # value < floor (= n_subs + n_prods + 1 = 5).
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: I
    end
    specs = EnzymeRates.init_mechanisms(rxn)
    floor_pc = 2 + 2 + 1   # n_subs + n_prods + 1
    # All init specs must have pc >= floor.
    for spec in specs
        @test spec.n_fit_params_estimate >= floor_pc
    end
    # At least one spec should have pc EQUAL to floor (proving the floor
    # was the binding constraint, not the formula).
    @test any(spec.n_fit_params_estimate == floor_pc for spec in specs)
end
```

If this assertion fails (no spec with pc == floor), the floor never actually fires — which would itself be a finding worth surfacing. Per protocol: surface, don't dodge.

- [ ] **Step 7.2:** Add `existing_de` exclusion-branch test for `_expand_add_allosteric_regulator`.

```julia
@testset "existing_de exclusion prevents adding bound dead-end as allo regulator" begin
    # SEED: uni-uni with rxn declaring I as :unknown role (could be either
    # dead-end or allosteric). First add I as a dead-end; then convert to
    # allosteric; then call _expand_add_allosteric_regulator. The move
    # should exclude :I (already bound as dead-end) → only fires on other
    # eligible regs. With only :I declared, result is empty.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        regulators: I    # :unknown role
        oligomeric_state: 2
    end
    init_specs = EnzymeRates.init_mechanisms(rxn)
    seed_spec = first(init_specs)
    # Add I as dead-end
    de_specs = EnzymeRates._expand_add_dead_end_regulator(seed_spec, rxn)
    @test !isempty(de_specs)
    plain_with_i = first(de_specs)
    # Convert to allosteric
    allo_specs = EnzymeRates._expand_to_allosteric(plain_with_i, rxn)
    spec = first(allo_specs)
    # Now call _expand_add_allosteric_regulator — :I is in existing_de,
    # so it's excluded. With no other declared regs, result is empty.
    result = EnzymeRates._expand_add_allosteric_regulator(spec, rxn)
    @test isempty(result)
end
```

- [ ] **Step 7.3:** Add new-site-multiplicity assertion to existing `_expand_add_allosteric_regulator` testsets.

Locate the "Allosteric uni-uni + first allo regulator R: 3 variants" testset. Add at the end of its preservation block:

```julia
# new-site multiplicity contract: when site_idx=0 (new site), the
# move sets multiplicity to spec.catalytic_n.
for r in result
    @test r.allosteric_multiplicities[end] == spec.catalytic_n
end
```

- [ ] **Step 7.4:** Run tests, verify green, commit.

```
test: add floor activation, existing_de exclusion, new-site multiplicity

- _apply_equivalence_grouping floor: bi-bi with dead-end inhibitor I has
  at least one init spec where pc == floor (proving the floor was the
  binding constraint, not the formula). If this test fails, the floor
  never fires — surface that as a finding.
- existing_de exclusion in _expand_add_allosteric_regulator: build a
  spec with I bound as dead-end, then convert to allosteric; calling
  the move with rxn declaring only :I yields empty (excluded).
- new-site multiplicity contract: when the move creates a new
  regulatory site (site_idx=0), multiplicity equals spec.catalytic_n.
  Asserted in the existing first-regulator testset's preservation block.
```

---

## Self-Review

After all 7 tasks:

- [ ] Test count grew or stayed roughly stable (within ±50 of 26255).
- [ ] No new vacuous tests introduced.
- [ ] All 7 commits land green.
- [ ] Round 2 findings #1, #2, #3, #4, #6, #9, #10, #11, #15, R1 #1, R1 #2, R1 #3, R1 #4, R1 #5, R1 #6, R1 #7, R1 #10 are addressed.

Findings deferred (not addressed in this pass):
- R2 #2 (23+ tautological `compile_mechanism(r) isa <Type>` checks) — structural pattern; would touch every move test. Documented for future cleanup.
- R2 #5 (multi-product release shape) — would require constructing a specific topology; lower-leverage given coverage via integration.
- R2 #12 (`_assert_spec_invariants` only at 2 callsites) — structural; would touch 30+ tests.
- R2 #13 (Δ=4 claim is loose) — accepted as upper-bound documentation.
- R2 #14 (helper masks formula bugs) — documented but the current pattern is intentional.
- R2 #16 (sample-size masking) — accepted given combinatorial cost.
- R1 #8, R1 #11, R1 #12, R1 #13, R1 #14, R1 #15 — minor inconsistencies, documentation, sample-size, etc.
