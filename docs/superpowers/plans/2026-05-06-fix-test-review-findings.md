# Fix Test Review Findings — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Address all 25 findings from the post-rewrite competing-subagent review of `test/test_mechanism_enumeration.jl`. Findings span vacuous/tautological tests, missing constructor/helper error-path coverage, independent-derivation rule violations, coverage gaps for tag combinations and mechanism shapes, and quick cleanups.

**Architecture:** Test-only changes. ~12 commits to `test/test_mechanism_enumeration.jl`, each addressing one cohesive group of related findings. No `src/` changes are planned; if a test surfaces a bug, follow the §6 protocol from the rewrite spec.

**Reference:**
- Review findings synthesis: see the conversation that produced this plan, or grep `git log --oneline` for review commits.
- Rewrite spec: `docs/superpowers/specs/2026-05-04-mechanism-enumeration-test-rewrite-design.md`
- Rewrite plan: `docs/superpowers/plans/2026-05-04-mechanism-enumeration-test-rewrite.md`

**Tech Stack:** Julia 1.x, Test stdlib, EnzymeRates.jl pipeline.

---

## Common verification

After every task's edits:

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -3
```

Expected: green. Baseline before this plan: 25593 tests passing. Each task adds tests; final count expected to grow by ~150-200 across all tasks. If a count regresses, the task absorbed/removed a vacuous test — verify the removal was intentional.

If a NEW assertion fails: surface as DONE_WITH_CONCERNS or BLOCKED per the bug-handling protocol. Do NOT weaken assertions to dodge bugs.

## Important conventions reminder

- Storage: dense Dict for `AllostericMechanismSpec`. Constructor validates density.
- Helpers: `mechanism_spec_from_mechanism_and_rxn`, `allosteric_spec_from_mechanism_and_rxn`. Both validate consistency internally.
- Use `EnzymeRates.compile_mechanism(...)` (not exported).
- Independent-derivation rule: every count/delta has a comment deriving the prediction from the seed and move semantics, not from observed output.
- File sections: 0 (infra), 1 (support fns), 2 (init/compile), 3 (base moves), 4 (allosteric moves), 5 (composition), 6 (integration), plus end-of-file orphan blocks.

---

## Task 1: Eliminate tautological and vacuous-pass tests

**Files:** `test/test_mechanism_enumeration.jl`

**Findings addressed:** #1 (exclude_regs vacuous), #2 (`r.base === spec || r.base == spec` tautology), #3 (method-dispatch fallback tautologies).

- [ ] **Step 1.1: Fix `exclude_regs` testset to actually exercise the kwarg's gating**

Locate the testset `@testset "exclude_regs kwarg suppresses regulator addition (negative)"` (around line 2081).

Current bug: the seed has no `:I` bound, and the rxn declares `:I`. The move would naturally yield empty regardless of `exclude_regs` because `:I` is not yet eligible. Need to make `:I` eligible BUT excluded.

Replace the testset body with a two-step setup:

```julia
@testset "exclude_regs kwarg overrides eligible regulators" begin
    # SEED: uni-uni catalytic with rxn declaring two dead-end inhibitors I, J.
    # Both are eligible to add (no overlap with bound metabolites).
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
        dead_end_inhibitors: I, J
    end
    spec = mechanism_spec_from_mechanism_and_rxn(m_seed, rxn)

    # Without exclude_regs: 2 eligible regs × 1 form pattern → 2 variants.
    # Property: both I and J appear across the result set.
    baseline = EnzymeRates._expand_add_dead_end_regulator(spec, rxn)
    @test length(baseline) == 2
    has_i_baseline = any(any(contains(string(sym), "I__reg")
                             for st in s.steps
                             for sym in Iterators.flatten((st.reactants, st.products)))
                         for s in baseline)
    has_j_baseline = any(any(contains(string(sym), "J__reg")
                             for st in s.steps
                             for sym in Iterators.flatten((st.reactants, st.products)))
                         for s in baseline)
    @test has_i_baseline && has_j_baseline

    # With exclude_regs=Set([:I]): only J is eligible → 1 variant, only J in result.
    # Derivation: exclude_regs filters from eligible_regs at the head of the move,
    # so :I is dropped before any form-pattern enumeration.
    excluded = EnzymeRates._expand_add_dead_end_regulator(spec, rxn; exclude_regs=Set([:I]))
    @test length(excluded) == 1
    has_i_excluded = any(any(contains(string(sym), "I__reg")
                             for st in s.steps
                             for sym in Iterators.flatten((st.reactants, st.products)))
                         for s in excluded)
    @test !has_i_excluded
    has_j_excluded = any(any(contains(string(sym), "J__reg")
                             for st in s.steps
                             for sym in Iterators.flatten((st.reactants, st.products)))
                         for s in excluded)
    @test has_j_excluded

    # With exclude_regs=Set([:I, :J]): no eligible regs → empty.
    @test isempty(EnzymeRates._expand_add_dead_end_regulator(
        spec, rxn; exclude_regs=Set([:I, :J])))
end
```

- [ ] **Step 1.2: Replace `r.base === spec || r.base == spec` tautology**

Locate the assertion in `_expand_to_allosteric` testset's preservation block (around line 2659):

```julia
# 5. preservation
for r in result
    @test r.base === spec || r.base == spec
end
```

Replace with explicit field-level checks:

```julia
# 5. preservation: base spec's reaction and step content unchanged
# (the move only attaches allosteric tags; the catalytic mechanism is
# carried through unmodified).
for r in result
    @test r.base.reaction === spec.reaction
    @test r.base.steps == spec.steps
    @test r.base.n_fit_params_estimate == spec.n_fit_params_estimate
end
```

- [ ] **Step 1.3: Strengthen method-dispatch fallback tests**

Three method-dispatch fallback tests trivially pass via Julia dispatch alone:
- `_expand_to_allosteric(::AllostericMechanismSpec, ...)` returns empty.
- `_expand_change_allo_state(::MechanismSpec, ...)` returns empty.
- `_expand_add_allosteric_regulator(::MechanismSpec, ...)` returns empty.

Locate each (search for `→ empty (negative)` testsets in sections 4a, 4b, 4c). For each, ADD a stronger assertion that exercises the dispatch BEHIND a meaningful round-trip:

For the `_expand_to_allosteric` already-allosteric case, add:

```julia
# Verify the result is empty AND that calling expand_mechanisms (which
# wraps all moves) on this spec doesn't add any to-allosteric variants
# from this spec — only the OTHER moves' outputs.
result_via_dispatch = EnzymeRates._expand_to_allosteric(spec, uni_uni_allo)
@test isempty(result_via_dispatch)
@test result_via_dispatch isa Vector{AllostericMechanismSpec}  # right element type
```

For the `_expand_change_allo_state(::MechanismSpec, ...)` case, similarly assert the return type:

```julia
result = EnzymeRates._expand_change_allo_state(spec, uni_uni_allo)
@test isempty(result)
@test result isa Vector{AllostericMechanismSpec}  # not Vector{MechanismSpec}
```

For `_expand_add_allosteric_regulator(::MechanismSpec, ...)`:

```julia
result = EnzymeRates._expand_add_allosteric_regulator(spec, uni_uni_allo_reg)
@test isempty(result)
@test result isa Vector{AllostericMechanismSpec}
```

The element-type assertion is the key strengthening: a regression that returned `Vector{MechanismSpec}` (e.g., from a copy-paste typo) would silently pass the old `isempty` test but fail this one.

- [ ] **Step 1.4: Run tests and commit**

Expected: green. Test count should increase by ~6-8.

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test: eliminate tautological and vacuous-pass tests

- Rewrite "exclude_regs kwarg" testset substantively. Old test had no
  :I bound, so exclude_regs=Set([:I]) yielded empty regardless of the
  kwarg — the test passed for the wrong reason. New test compares
  baseline (no kwarg) to filtered (kwarg) and verifies :I appears in
  baseline but is excluded by the kwarg.
- Replace `r.base === spec || r.base == spec` tautology in
  _expand_to_allosteric preservation with explicit field-level checks.
  The disjunction always reduced to === ∨ === since == falls back
  to === for structs.
- Strengthen the three method-dispatch fallback tests
  (_expand_to_allosteric on allosteric input, _expand_change_allo_state
  on plain spec, _expand_add_allosteric_regulator on plain spec) with
  return-type assertions. Previously these passed by Julia dispatch
  alone — a regression returning the wrong element type would have
  been invisible.

Surfaced by post-rewrite review (findings #1, #2, #3).
EOF
)"
```

---

## Task 2: Add error-path tests for constructor and helpers

**Files:** `test/test_mechanism_enumeration.jl`

**Findings addressed:** #4 (AllostericMechanismSpec constructor + helper consistency-validation error paths untested).

- [ ] **Step 2.1: Add error-path testset for `AllostericMechanismSpec` constructor**

Insert RIGHT BEFORE the existing `@testset "allosteric_spec_from_mechanism_and_rxn round-trip"` (in section 0, near line 110):

```julia
@testset "AllostericMechanismSpec constructor density validation" begin
    # The constructor rejects sparse Dicts: every kinetic group used in
    # base.steps must have a group_tags entry; every ligand listed in
    # allosteric_reg_sites must have a reg_ligand_tags entry. This guards
    # the dense-storage invariant that both spec and compiled mechanism
    # use throughout the pipeline.
    base = MechanismSpec(uni_uni_rxn,
        [StepSpec([:E, :P], [:E_P], true, 1),
         StepSpec([:E, :S], [:E_S], true, 2),
         StepSpec([:E_S], [:E_P], false, 3)],
        3)

    # Missing group_tags entry (group 3 omitted).
    @test_throws ErrorException AllostericMechanismSpec(
        base, 2,
        Vector{Symbol}[], Int[],
        Dict(1 => :EqualRT, 2 => :EqualRT),  # group 3 missing
        Dict{Symbol, Symbol}(),
        4)

    # Missing reg_ligand_tags entry (ligand R declared in reg_sites but
    # not in reg_ligand_tags).
    @test_throws ErrorException AllostericMechanismSpec(
        base, 2,
        Vector{Symbol}[[:R]], Int[2],
        Dict(1 => :EqualRT, 2 => :EqualRT, 3 => :EqualRT),
        Dict{Symbol, Symbol}(),  # :R missing
        4)

    # Both Dicts complete → constructor succeeds.
    valid = AllostericMechanismSpec(
        base, 2,
        Vector{Symbol}[[:R]], Int[2],
        Dict(1 => :EqualRT, 2 => :EqualRT, 3 => :EqualRT),
        Dict(:R => :EqualRT),
        4)
    @test valid isa AllostericMechanismSpec
end
```

- [ ] **Step 2.2: Add error-path testset for the round-trip helpers**

Insert RIGHT BEFORE the existing `@testset "allosteric_spec_from_mechanism_and_rxn round-trip"`:

```julia
@testset "spec-from-mechanism helpers reject inconsistent inputs" begin
    # Both helpers validate that mechanism and reaction agree on
    # substrates/products/regulators (and oligomeric_state for the
    # allosteric variant). Mismatches throw ErrorException at the helper
    # level, before any spec is constructed.

    m_uu = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ E_P
            E + S ⇌ E_S
            E_S <--> E_P
        end
    end

    # Substrate name mismatch
    rxn_wrong_sub = @enzyme_reaction begin
        substrates: X[C]   # not :S
        products: P[C]
    end
    @test_throws ErrorException mechanism_spec_from_mechanism_and_rxn(m_uu, rxn_wrong_sub)

    # Product name mismatch
    rxn_wrong_prod = @enzyme_reaction begin
        substrates: S[C]
        products: Y[C]   # not :P
    end
    @test_throws ErrorException mechanism_spec_from_mechanism_and_rxn(m_uu, rxn_wrong_prod)

    # Regulator subset rule: m_with_I has :I bound; rxn lacks :I → reject.
    m_with_I = @enzyme_mechanism begin
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
    rxn_no_I = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    @test_throws ErrorException mechanism_spec_from_mechanism_and_rxn(m_with_I, rxn_no_I)

    # Allosteric helper: oligomeric_state mismatch
    m_allo_2 = @allosteric_mechanism begin
        substrates: S
        products: P
        site(:catalytic, 2): begin
            steps: begin
                E + P ⇌ E_P    :: EqualRT
                E + S ⇌ E_S    :: EqualRT
                E_S <--> E_P   :: EqualRT
            end
        end
    end
    rxn_oligo_4 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        oligomeric_state: 4   # mismatch — m_allo_2 has catalytic_n=2
    end
    @test_throws ErrorException allosteric_spec_from_mechanism_and_rxn(m_allo_2, rxn_oligo_4)
end
```

- [ ] **Step 2.3: Run tests and commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test: add error-path tests for constructor and round-trip helpers

The AllostericMechanismSpec constructor errors when group_tags or
reg_ligand_tags are not dense for the spec's content. The two round-trip
helpers (mechanism_spec_from_mechanism_and_rxn and
allosteric_spec_from_mechanism_and_rxn) error on substrate/product/
regulator/oligomeric_state mismatches between mechanism and reaction.
None of these error branches had tests — a regression that dropped a
validation check would have produced malformed specs silently.

Adds @test_throws cases for each error path plus a positive control
verifying the constructor accepts dense input.

Surfaced by post-rewrite review (finding #4).
EOF
)"
```

---

## Task 3: Fix precision violations and misleading comments/names

**Files:** `test/test_mechanism_enumeration.jl`

**Findings addressed:** #5 (Δ-multiset precision loss in "Two regulators with site options"), #6 (comments codify wrong number contradicting assertion in Bi-Bi random and Bi-Bi PP testsets), #7 ("Bi-Bi PP: 3 dead-end forms" testset name lies).

- [ ] **Step 3.1: Tighten "Two regulators with site options" Δ-multiset assertion**

Locate the testset around line 2855-2885 in `_expand_add_allosteric_regulator`. Find the assertion that currently uses `1 in deltas` / `2 in deltas`. Replace with:

```julia
# 2. Δ params: 7 variants total.
#   Non-:EqualRT branch (3 tags × 2 sites = 6):
#     :OnlyR/:OnlyT at new site → +1 each (cheap-tag delta vs :EqualRT base = 0,
#       plus +1 for the new K_R binding param) = 4 variants × +1
#     :NonequalRT at new site or existing site → +2 each (cost(:NonequalRT)=2,
#       plus +1 for new K) = 2 variants × +2
#   :EqualRT-at-existing-site branch (1 variant; only fires because R1 is
#   non-:EqualRT in the seed): +1 (cost(:EqualRT)=1, plus +1 for new K) = 1 × +1
# Sorted multiset: [1, 1, 1, 1, 1, 1, 2]
deltas = sort([r.n_fit_params_estimate -
               spec.n_fit_params_estimate for r in result])
@test deltas == [1, 1, 1, 1, 1, 1, 2]
```

If running tests reveals that the actual delta multiset is different (e.g., the derivation is wrong), STOP and surface as DONE_WITH_CONCERNS — don't change the assertion to match the observed output. Re-derive from the move's source.

- [ ] **Step 3.2: Fix `Bi-Bi random: 4 dead-end forms` comment**

Locate the testset (around line 877-907 in `_expand_substrate_product_dead_ends`). The current comment derives `2^4 = 16 variants` but the assertion is `length(result) == 7`. The 16 number is the unfiltered upper bound; the actual 7 comes from competition-pattern dedup.

Replace the comment block with derivation that matches the assertion:

```julia
# 4 unique mixed-substrate-product forms across competition patterns.
# Competition patterns for bi-bi (2 subs × 2 prods): 7 patterns
# (the count from _competition_patterns(2, 2)). Each pattern produces
# a distinct dead-end-form set:
#   {A↔P, B↔Q}: forbids E_A_P, E_B_Q → emits {E_A_Q, E_B_P}
#   {A↔Q, B↔P}: forbids E_A_Q, E_B_P → emits {E_A_P, E_B_Q}
#   ... (one set per pattern, all distinct)
#   {A↔P, A↔Q, B↔P, B↔Q}: forbids all → emits {} (bare topology)
# All 7 sets are distinct → 7 variants after dedup.
```

- [ ] **Step 3.3: Fix `Bi-Bi PP: 3 dead-end forms` testset name AND comment**

Locate the testset (around line 931 in `_expand_substrate_product_dead_ends`). Body derives 5 dead-end forms (E_A_P, E_A_Q, E_B_Q from E-side + Estar_B_P, Estar_B_Q from Estar-side) producing 7 competition-filtered variants.

Rename:

```julia
@testset "Bi-Bi Ping-Pong: 5 dead-end forms → 7 variants" begin
    # 5 dead-end forms total (E-side: E_A_P, E_A_Q, E_B_Q; Estar-side:
    # Estar_B_P, Estar_B_Q). 7 competition patterns; each yields a
    # distinct dead-end-form set after dedup → 7 variants.
    ...
end
```

Update the inline comment to derive 7 from the competition-pattern formula, not from `2^5 = 32`.

- [ ] **Step 3.4: Run tests and commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test: tighten precision violations and fix misleading comments

- Replace the loose `1 in deltas, 2 in deltas` assertion in
  "Two regulators with site options" with the exact-multiset
  `deltas == [1, 1, 1, 1, 1, 1, 2]`. The previous assertion would
  have passed for grossly wrong outputs (5 ones + 2 twos, etc.) —
  violating the spec's independent-derivation rule.
- Fix "Bi-Bi random: 4 dead-end forms" comment that derived
  `2^4 = 16 variants` but asserted length == 7. Rewrite the comment
  to derive 7 from competition-pattern dedup.
- Rename "Bi-Bi PP: 3 dead-end forms" → "5 dead-end forms → 7 variants"
  to match the body derivation.

Surfaced by post-rewrite review (findings #5, #6, #7).
EOF
)"
```

---

## Task 4: Add `:OnlyT` input seeds across allosteric moves

**Files:** `test/test_mechanism_enumeration.jl`

**Findings addressed:** #8 (`:OnlyT` never tested as input).

- [ ] **Step 4.1: Add `:OnlyT` group seed to `_expand_re_to_ss`**

Locate the `_expand_re_to_ss` testset (around line 1568). After the `:OnlyR group: Δ=+1` sub-testset, add:

```julia
@testset "AllostericMechanismSpec — :OnlyT group: Δ=+1" begin
    # SEED: uni-uni allosteric with one binding group :OnlyT (active only
    # in T-state; absent in R-state). :OnlyT is forbidden for catalytic-
    # iso-only groups but allowed for binding groups. After RE→SS, the
    # converted group's tag is preserved as :OnlyT.
    # Δ derivation: _re_to_ss_delta returns 1 for any non-:NonequalRT tag;
    # :OnlyT is cheap → Δ = +1.
    m_seed = @allosteric_mechanism begin
        substrates: S
        products: P
        site(:catalytic, 2): begin
            steps: begin
                E + P ⇌ E_P       :: EqualRT
                E + S ⇌ E_S       :: OnlyT
                E_S <--> E_P      :: EqualRT
            end
        end
    end
    spec = allosteric_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo)

    result = EnzymeRates._expand_re_to_ss(spec)

    # 1. count: 2 RE binding groups → 2 variants.
    @test length(result) == 2

    # 2. Δ params: +1 each (cheap-tag rule).
    for r in result
        @test r.n_fit_params_estimate ==
            spec.n_fit_params_estimate + 1
    end

    # 3. compilability
    for r in result
        @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
    end

    # 4. property-style: tag preservation. The flipped group keeps :OnlyT.
    for r in result
        @test r.group_tags == spec.group_tags
    end

    # 5. preservation
    for r in result
        @test r.catalytic_n == spec.catalytic_n
        @test r.allosteric_reg_sites == spec.allosteric_reg_sites
        @test r.allosteric_multiplicities == spec.allosteric_multiplicities
        @test r.reg_ligand_tags == spec.reg_ligand_tags
        @test r.base.reaction === spec.base.reaction
    end
end
```

- [ ] **Step 4.2: Add `:OnlyT` regulator-ligand seed to `_expand_change_allo_state`**

Locate the `_expand_change_allo_state` testset (around line 2995). After the existing seeds, add:

```julia
@testset ":OnlyT regulator-ligand relaxation" begin
    # SEED: uni-uni allosteric with one regulator R tagged :OnlyT.
    # _expand_change_allo_state should produce variants for each
    # non-:NonequalRT entry, including the :OnlyT ligand. Δ for the
    # ligand-relaxation variant: cost(:NonequalRT) - cost(:OnlyT) = +1.
    m_seed = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R::OnlyT
        site(:catalytic, 2): begin
            steps: begin
                E + P ⇌ E_P       :: EqualRT
                E + S ⇌ E_S       :: EqualRT
                E_S <--> E_P      :: EqualRT
            end
        end
    end
    spec = allosteric_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo_reg)

    result = EnzymeRates._expand_change_allo_state(spec, uni_uni_allo_reg)

    # 1. count: 3 group_tags entries (all :EqualRT) + 1 reg_ligand_tags
    # entry (:OnlyT) → 4 variants.
    @test length(result) == 4

    # 2. Δ params: 2 RE-binding-group EqualRT relaxations (+1 each),
    # 1 SS-iso EqualRT relaxation (+2), 1 reg-ligand :OnlyT relaxation
    # (cost(:NonequalRT) - cost(:OnlyT) = 2 - 1 = +1). Sorted: [1, 1, 1, 2].
    deltas = sort([r.n_fit_params_estimate -
                   spec.n_fit_params_estimate for r in result])
    @test deltas == [1, 1, 1, 2]

    # 3. compilability
    for r in result
        @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
    end

    # 4. property-style: exactly one ligand-relaxation variant has the
    # R tag flipped to :NonequalRT.
    n_r_relaxed = count(r -> r.reg_ligand_tags[:R] == :NonequalRT, result)
    @test n_r_relaxed == 1
end
```

- [ ] **Step 4.3: Add `:OnlyT` regulator at existing reg site to `_expand_add_allosteric_regulator`**

Locate the `_expand_add_allosteric_regulator` testset. After the existing seeds, add:

```julia
@testset "Adding :EqualRT R2 at site with :OnlyT R1" begin
    # SEED: allosteric uni-uni with R1::OnlyT already present at site 1.
    # Adding R2 should enumerate non-:EqualRT tags (×2 sites: new + existing) +
    # :EqualRT at existing site (because R1 is non-:EqualRT, the
    # :EqualRT-at-existing branch fires). Total: 3×2 + 1 = 7 variants.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R1, R2
        oligomeric_state: 2
    end
    m_seed = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R1::OnlyT
        site(:catalytic, 2): begin
            steps: begin
                E + P ⇌ E_P       :: EqualRT
                E + S ⇌ E_S       :: EqualRT
                E_S <--> E_P      :: EqualRT
            end
        end
    end
    spec = allosteric_spec_from_mechanism_and_rxn(m_seed, rxn)

    result = EnzymeRates._expand_add_allosteric_regulator(spec, rxn)

    # 1. count: 3 non-:EqualRT tags × 2 site options + 1 :EqualRT-at-existing
    # = 7 variants.
    @test length(result) == 7

    # 2. Δ params: same multiset as the analogous :OnlyR seed
    # ([1,1,1,1,1,1,2] — see Two regulators with site options).
    deltas = sort([r.n_fit_params_estimate -
                   spec.n_fit_params_estimate for r in result])
    @test deltas == [1, 1, 1, 1, 1, 1, 2]

    # 3. compilability
    for r in result
        @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
    end

    # 4. property-style: at least one variant has R2 :EqualRT at site 1
    # (the EqualRT-at-existing branch, gated on R1 being non-:EqualRT).
    has_eq_at_site1 = any(result) do r
        :R2 in r.allosteric_reg_sites[1] && r.reg_ligand_tags[:R2] == :EqualRT
    end
    @test has_eq_at_site1
end
```

- [ ] **Step 4.4: Run tests and commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test: add :OnlyT input seeds across allosteric moves

The :OnlyT tag was used only as expected output in round-trip tests but
never as input to expansion moves. Added three sub-testsets:

- _expand_re_to_ss with one :OnlyT binding group (cheap-tag delta +1).
- _expand_change_allo_state with one :OnlyT regulator-ligand
  (verifies ligand-relaxation delta +1).
- _expand_add_allosteric_regulator with :OnlyT R1 at existing site
  (verifies the :EqualRT-at-existing-site branch fires when R1 is
  non-:EqualRT).

Surfaced by post-rewrite review (finding #8).
EOF
)"
```

---

## Task 5: Add ping-pong (Estar) seeds across expansion moves

**Files:** `test/test_mechanism_enumeration.jl`

**Findings addressed:** #9 (ping-pong seeds missing from `_expand_re_to_ss`, `_expand_to_allosteric`, allosteric moves).

- [ ] **Step 5.1: Add ping-pong seed to `_expand_re_to_ss`**

After the existing bi-bi-multi-step testset in `_expand_re_to_ss`:

```julia
@testset "MechanismSpec — bi-bi ping-pong: 3 RE groups → 3 variants" begin
    # SEED: bi-bi ping-pong topology with Estar (residual) form.
    # 3 singleton RE groups (A binding, Q release on E-side; B binding
    # via Estar; A→Estar iso; Estar→E iso). Iso steps are SS.
    m_seed = @enzyme_mechanism begin
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
    spec = mechanism_spec_from_mechanism_and_rxn(m_seed, bi_bi_pp_rxn)

    result = EnzymeRates._expand_re_to_ss(spec)

    # 1. count: RE groups = E+A, Estar+B, E+Q, Estar+P, Estar_B → E_Q
    # (5 groups). The iso E_A↔Estar_A_P is SS so excluded. Among the
    # remaining 5, count how many are RE in the seed: E+A, Estar+B,
    # E+Q, Estar+P are RE; Estar_B → E_Q is RE; iso is SS. So 5 RE
    # groups → 5 variants. (Note: derivation depends on exact group
    # numbering in the macro — verify via the seed's compiled reactions
    # tuple if it differs.)
    @test length(result) == 5

    # 2. Δ params: +1 each (plain MechanismSpec).
    for r in result
        @test r.n_fit_params_estimate ==
            spec.n_fit_params_estimate + 1
    end

    # 3. compilability
    for r in result
        @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
    end

    # 5. preservation
    for r in result
        @test r.reaction === spec.reaction
    end
end
```

If the count assertion fails, surface as DONE_WITH_CONCERNS — do NOT modify the count to match observed output. Re-derive from the seed's actual kinetic_group structure.

- [ ] **Step 5.2: Add ping-pong seed to `_expand_to_allosteric`**

After the bi-bi sequential testset in `_expand_to_allosteric`:

```julia
@testset "Bi-bi ping-pong: 6 groups → 7 variants" begin
    # SEED: bi-bi ping-pong topology mapped to allosteric.
    # 6 kinetic groups (5 RE binding + 1 iso). Plus the second iso step
    # makes 7 groups total. Move emits 1 baseline + 6 :OnlyR variants.
    m_seed = @enzyme_mechanism begin
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
    bi_bi_pp_allo_rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        oligomeric_state: 2
    end
    spec = mechanism_spec_from_mechanism_and_rxn(m_seed, bi_bi_pp_allo_rxn)
    result = EnzymeRates._expand_to_allosteric(spec, bi_bi_pp_allo_rxn)

    # 1. count: 6 kinetic groups (5 binding + 1 iso) → 7 variants
    # (1 baseline + 6 :OnlyR per group). Verify the actual group count
    # from the seed before accepting this number.
    n_groups = length(unique(s.kinetic_group for s in spec.steps))
    @test length(result) == n_groups + 1

    # 2. Δ params: +1 (just L) per variant.
    for r in result
        @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
    end

    # 3. compilability
    for r in result
        @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
    end
end
```

- [ ] **Step 5.3: Run tests and commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test: add ping-pong (Estar) seeds across expansion moves

Spec §3 listed ping-pong as a seed for both _expand_re_to_ss and
_expand_to_allosteric. Only _expand_add_dead_end_regulator had a
ping-pong seed; the others were untested with Estar-form mechanisms.

Adds bi-bi ping-pong seeds for both moves. The Estar form-name handling
in _bound_metabolites_at_forms and downstream code is now exercised by
move-level tests, not just topology-generation tests.

Surfaced by post-rewrite review (finding #9).
EOF
)"
```

---

## Task 6: Complete `_expand_add_allosteric_regulator` seed battery

**Files:** `test/test_mechanism_enumeration.jl`

**Findings addressed:** #10 (3 missing planned seeds: two-reg-different-sites, product-as-allo-reg, all-regs-added → empty).

- [ ] **Step 6.1: Add two-regulators-at-different-sites seed**

Inside `_expand_add_allosteric_regulator` testset:

```julia
@testset "Two regulators at different sites" begin
    # SEED: allosteric uni-uni with R1 already present at site 1; we add R2.
    # R2 can go to: a NEW site (site_idx=0) OR R1's existing site (site_idx=1).
    # The cross-site placement (site_idx ≥ 2 in the source) requires 2+
    # existing sites; this seed has only 1 existing, so we focus on the
    # site_idx=0 vs site_idx=1 distinction.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R1, R2
        oligomeric_state: 2
    end
    m_seed = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R1::OnlyR
        site(:catalytic, 2): begin
            steps: begin
                E + P ⇌ E_P       :: EqualRT
                E + S ⇌ E_S       :: EqualRT
                E_S <--> E_P      :: EqualRT
            end
        end
    end
    spec = allosteric_spec_from_mechanism_and_rxn(m_seed, rxn)

    result = EnzymeRates._expand_add_allosteric_regulator(spec, rxn)

    # 1. count: 3 non-:EqualRT tags × 2 site options + 1 :EqualRT-at-existing
    # = 7 variants (same as the existing "Two regulators with site options"
    # but verifying via two-different-sites placement specifically).
    @test length(result) == 7

    # 4. property-style: separate the new-site and existing-site placements.
    new_site_variants = filter(r -> length(r.allosteric_reg_sites) == 2, result)
    existing_site_variants = filter(r -> length(r.allosteric_reg_sites) == 1, result)
    @test length(new_site_variants) == 3   # 3 tags × new site
    @test length(existing_site_variants) == 4  # 3 tags + 1 :EqualRT
end
```

- [ ] **Step 6.2: Add product-as-allosteric-regulator overlap seed**

```julia
@testset "Product-as-allosteric-regulator overlap" begin
    # SEED: uni-uni allosteric where product P is ALSO declared as an
    # allosteric regulator. Adding :P as allo regulator should produce
    # 3 tag variants × 1 site option = 3 variants. Verifies the move
    # treats name-overlapping ligand and product as independent.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: P
        oligomeric_state: 2
    end
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
    spec = allosteric_spec_from_mechanism_and_rxn(m_seed, rxn)

    result = EnzymeRates._expand_add_allosteric_regulator(spec, rxn)

    # 1. count: 3 non-:EqualRT tags × 1 new site = 3 variants.
    @test length(result) == 3

    # 2. Δ params: [1, 1, 2] (sorted: 2 cheap + 1 :NonequalRT).
    deltas = sort([r.n_fit_params_estimate -
                   spec.n_fit_params_estimate for r in result])
    @test deltas == [1, 1, 2]

    # 3. compilability
    for r in result
        @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
    end

    # 4. property: :P appears in allosteric_reg_sites
    for r in result
        @test any(:P in site for site in r.allosteric_reg_sites)
    end
end
```

- [ ] **Step 6.3: Add all-regs-already-added → empty seed**

```julia
@testset "All declared regs already present → empty (negative)" begin
    # SEED: allosteric uni-uni with R already added. The reaction declares
    # only R as a regulator. eligible_regs computation excludes R because
    # it's in existing_allo → result is empty.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R
        oligomeric_state: 2
    end
    m_seed = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R::OnlyR
        site(:catalytic, 2): begin
            steps: begin
                E + P ⇌ E_P       :: EqualRT
                E + S ⇌ E_S       :: EqualRT
                E_S <--> E_P      :: EqualRT
            end
        end
    end
    spec = allosteric_spec_from_mechanism_and_rxn(m_seed, rxn)
    @test isempty(EnzymeRates._expand_add_allosteric_regulator(spec, rxn))
end
```

- [ ] **Step 6.4: Run tests and commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test: complete _expand_add_allosteric_regulator seed battery

Spec §3 listed seeds: one-reg, two-reg-same-site, two-reg-different-sites,
existing :EqualRT site, substrate-as-allo-reg, product-as-allo-reg,
non-allosteric-negative, all-regs-added-negative. Three were missing:

- Two regulators at different sites (verifies site_idx=0 vs site_idx=1
  placement separately).
- Product-as-allosteric-regulator overlap (verifies name-overlap between
  product and allosteric ligand is handled correctly).
- All declared regs already added → empty (verifies the existing_allo
  filter at the head of the move).

Surfaced by post-rewrite review (finding #10).
EOF
)"
```

---

## Task 7: Add SS-allosteric seed for `_expand_split_kinetic_group`

**Files:** `test/test_mechanism_enumeration.jl`

**Findings addressed:** #11 (SS × :NonequalRT path Δ=4 in `_split_group_delta` never exercised).

- [ ] **Step 7.1: Add SS-allosteric split seed**

Inside `_expand_split_kinetic_group`, after the existing allosteric-tag-inheritance testset:

```julia
@testset "AllostericMechanismSpec — SS multi-step :NonequalRT split: Δ=+4" begin
    # SEED: bi-bi allosteric where one multi-step group is BOTH SS AND
    # :NonequalRT. _split_group_delta returns 4 for this case (factor 2
    # for SS × factor 2 for :NonequalRT R/T-state pair).
    m_seed = @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        site(:catalytic, 2): begin
            steps: begin
                (E + A <--> E_A, E_B + A <--> E_A_B)        :: NonequalRT
                (E + B ⇌ E_B, E_A + B ⇌ E_A_B)             :: EqualRT
                (E + P ⇌ E_P, E_Q + P ⇌ E_P_Q)             :: EqualRT
                (E + Q ⇌ E_Q, E_P + Q ⇌ E_P_Q)             :: EqualRT
                E_A_B <--> E_P_Q                            :: EqualRT
            end
        end
    end
    bi_bi_allo_rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        oligomeric_state: 2
    end
    spec = allosteric_spec_from_mechanism_and_rxn(m_seed, bi_bi_allo_rxn)

    result = EnzymeRates._expand_split_kinetic_group(spec)

    # 1. count: 4 multi-step groups (A-binding SS×2 :NonequalRT,
    # B-binding RE×2 :EqualRT, P-binding RE×2 :EqualRT, Q-binding RE×2 :EqualRT).
    # 4 × 2 members = 8 variants.
    @test length(result) == 8

    # 2. Δ params:
    # - SS × :NonequalRT split: factor 2 (SS) × factor 2 (NonequalRT) = +4
    #   → 2 variants × +4
    # - RE × :EqualRT split: factor 1 × factor 1 = +1
    #   → 6 variants × +1
    deltas = sort([r.n_fit_params_estimate -
                   spec.n_fit_params_estimate for r in result])
    @test deltas == [1, 1, 1, 1, 1, 1, 4, 4]

    # 3. compilability
    for r in result
        @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
    end

    # 4. tag inheritance: split's new group inherits parent's tag.
    pre_groups = Set(s.kinetic_group for s in spec.base.steps)
    for r in result
        post_groups = Set(s.kinetic_group for s in r.base.steps)
        new_g = only(setdiff(post_groups, pre_groups))
        # find the parent group
        pre_counts = Dict(g => count(s -> s.kinetic_group == g, spec.base.steps)
                          for g in pre_groups)
        post_counts = Dict(g => count(s -> s.kinetic_group == g, r.base.steps)
                           for g in pre_groups)
        old_g = only(g for g in pre_groups
                     if post_counts[g] < pre_counts[g])
        @test r.group_tags[new_g] == spec.group_tags[old_g]
    end
end
```

- [ ] **Step 7.2: Run tests and commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test(split_kinetic_group): add SS multi-step :NonequalRT seed (Δ=+4)

The _split_group_delta function returns 4 for SS × :NonequalRT splits
(factor 2 × factor 2). No previous test exercised this path. Added a
bi-bi allosteric seed where one multi-step group is both SS and
:NonequalRT; verifies the +4 delta along with tag inheritance.

Surfaced by post-rewrite review (finding #11).
EOF
)"
```

---

## Task 8: Tighten composition-layer property assertions

**Files:** `test/test_mechanism_enumeration.jl`

**Findings addressed:** #12 (`expand_mechanisms` sub-testsets use loose `!isempty/any` predicates).

- [ ] **Step 8.1: Strengthen "Allosteric expansion included" testset**

Locate the testset (around line 3205). Replace the loose `has_allo` assertion with a count-based check. The seed is uni-uni allo, which produces n_groups+1=4 allosteric variants from `_expand_to_allosteric` alone. Add:

```julia
# Before: just `@test has_allo`. Now derive the count.
allo_count = sum(count(s -> s isa AllostericMechanismSpec, ss)
                 for (_, ss) in result)
# Derivation: _expand_to_allosteric on uni-uni init produces n_groups+1=4
# allosteric variants. Other moves don't produce allosteric output from
# a plain MechanismSpec input. So at least 4 allosteric variants exist.
@test allo_count >= 4
```

- [ ] **Step 8.2: Strengthen "Allosteric rewrap preserves structure" testset**

Locate the testset (around line 3245). Currently asserts only `has_rewrapped`. Add a structural check that the rewrapped specs preserve the input's `catalytic_n` and `allosteric_reg_sites`:

```julia
allo_results = filter(s -> s isa AllostericMechanismSpec,
                     vcat([ss for (_, ss) in result]...))
@test !isempty(allo_results)

# Property: every rewrapped allosteric result has the same catalytic_n
# and a base.reaction matching the input's. The base.steps may differ
# (a base move may have changed them) but the allosteric-side metadata
# is preserved.
for r in allo_results
    @test r.catalytic_n == spec.catalytic_n
    @test r.base.reaction === spec.base.reaction
end
```

- [ ] **Step 8.3: Strengthen "Multiple levels populated" testset**

Locate the testset (around line 3398). Currently asserts `length(results) >= 2`. Strengthen to check that param-count keys are CONSECUTIVE (or document why gaps are allowed):

```julia
results = enumerate_all(uni_uni_rxn; max_params=8)
@test length(results) >= 2

# Param-count buckets should be a contiguous range from min to max
# (every level reachable by the +1 / +2 expansion moves should appear).
pcs = sort(collect(keys(results)))
@test all(pcs[i+1] - pcs[i] in (1, 2) for i in 1:length(pcs)-1)
```

- [ ] **Step 8.4: Run tests and commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test(expand_mechanisms): tighten loose !isempty/any property assertions

Three sub-testsets in the composition layer used predicates so loose
that they passed if even one element existed:
- "Allosteric expansion included" → now asserts at least n_groups+1=4
  allosteric variants exist (the count from _expand_to_allosteric).
- "Allosteric rewrap preserves structure" → now asserts catalytic_n
  and base.reaction are preserved across rewrap.
- "Multiple levels populated" → now asserts param-count keys form a
  contiguous range with gaps of 1 or 2 (the expansion-move deltas).

Surfaced by post-rewrite review (finding #12).
EOF
)"
```

---

## Task 9: Add coverage for under-tested mechanism shapes

**Files:** `test/test_mechanism_enumeration.jl`

**Findings addressed:** #13 (multi-product iso step never seeded), #14 (ter-ter and pyruvate reactions never seed expansion moves), #20 (non-default site multiplicities), #21 (RE→SS substrate-overlap missing cross-type), #22 (`_expand_change_allo_state` lacks multi-reg-ligands seed).

- [ ] **Step 9.1: Add ter-ter sequential seed for `_expand_re_to_ss`**

```julia
@testset "MechanismSpec — ter-ter sequential" begin
    # SEED: ter-ter sequential ordered. 6 RE binding steps + 1 SS iso = 7 groups.
    # _expand_re_to_ss fires per RE group → 6 variants.
    m_seed = @enzyme_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            E + A ⇌ E_A
            E_A + B ⇌ E_A_B
            E_A_B + D ⇌ E_A_B_D
            E + R ⇌ E_R
            E_R + Q ⇌ E_Q_R
            E_Q_R + P ⇌ E_P_Q_R
            E_A_B_D <--> E_P_Q_R
        end
    end
    spec = mechanism_spec_from_mechanism_and_rxn(m_seed, ter_ter_rxn)
    result = EnzymeRates._expand_re_to_ss(spec)
    @test length(result) == 6
    for r in result
        @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
        @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
    end
end
```

- [ ] **Step 9.2: Add multi-product iso seed**

Multi-product iso is when a single step like `E_S <--> E + P + Q` releases multiple products. This requires building a topology that has such a step. Use the `init_mechanisms(uni_bi_rxn)` output that may produce them, OR construct via `@enzyme_mechanism` with a step like `E_S <--> E_P_Q` followed by individual product release steps.

Actually, looking at the topology generator: `_release_products!` releases products one at a time, so a multi-product iso form has subsequent release steps. The "iso step" itself is `E_S → E_P_Q` (one product complex). The release of P from E_P_Q is a separate step. So this is already covered by the existing seeds.

If on closer inspection multi-product release IS a distinct shape worth testing: add a uni-bi seed via `init_mechanisms(uni_bi_rxn)` for a topology with multi-product release, document the structure, and call `_expand_re_to_ss`/`_expand_split_kinetic_group` on it.

If after investigation this is found to already be covered by existing seeds, add a comment to that effect in the relevant testset and skip this step.

- [ ] **Step 9.3: Add non-default-multiplicity allosteric seed**

```julia
@testset "Non-default site multiplicity (catalytic_n=4, reg site multiplicity=2)" begin
    # SEED: catalytic 4-mer with regulator at multiplicity-2 site (less than catalytic_n).
    # _expand_change_allo_state and other allo moves should preserve site
    # multiplicities independently of catalytic_n.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R
        oligomeric_state: 4
    end
    m_seed = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R::OnlyR
        site(:catalytic, 4): begin
            steps: begin
                E + P ⇌ E_P       :: EqualRT
                E + S ⇌ E_S       :: EqualRT
                E_S <--> E_P      :: EqualRT
            end
        end
        site(:regulatory, 2): begin
            ligands: R
        end
    end
    spec = allosteric_spec_from_mechanism_and_rxn(m_seed, rxn)
    @test spec.catalytic_n == 4
    @test spec.allosteric_multiplicities == [2]

    # _expand_change_allo_state should preserve multiplicities.
    result = EnzymeRates._expand_change_allo_state(spec, rxn)
    @test !isempty(result)
    for r in result
        @test r.catalytic_n == 4
        @test r.allosteric_multiplicities == [2]
    end
end
```

- [ ] **Step 9.4: Add cross-type substrate-overlap seed for `_expand_re_to_ss`**

```julia
@testset "AllostericMechanismSpec — substrate-as-dead-end-I overlap" begin
    # Allosteric counterpart of the substrate-as-I overlap test; verifies
    # that the move correctly handles the overlap when the spec is
    # AllostericMechanismSpec (cross-type).
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: S
        oligomeric_state: 2
    end
    init_specs = EnzymeRates.init_mechanisms(rxn)
    seed_spec = first(init_specs)
    de_specs = EnzymeRates._expand_add_dead_end_regulator(seed_spec, rxn)
    @test !isempty(de_specs)
    plain_spec = first(de_specs)
    # Convert to allosteric
    allo_specs = EnzymeRates._expand_to_allosteric(plain_spec, rxn)
    @test !isempty(allo_specs)
    spec = first(allo_specs)

    result = EnzymeRates._expand_re_to_ss(spec)
    # Same count derivation as the plain-spec overlap test: 3 RE groups
    # (substrate-S binding, product-P binding, dead-end-S__reg binding) → 3 variants.
    @test length(result) == 3

    for r in result
        @test r isa AllostericMechanismSpec
        @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
        @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
    end
end
```

- [ ] **Step 9.5: Add multi-regulator-ligand seed for `_expand_change_allo_state`**

```julia
@testset "Multiple regulator ligands at independent tags" begin
    # SEED: allosteric uni-uni with two regulators R1::OnlyR, R2::OnlyT.
    # Move should produce a relaxation variant for EACH non-:NonequalRT
    # entry: 3 group-tag (all :EqualRT) + 2 reg-ligand-tag (R1, R2) = 5 variants.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R1, R2
        oligomeric_state: 2
    end
    m_seed = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R1::OnlyR, R2::OnlyT
        site(:catalytic, 2): begin
            steps: begin
                E + P ⇌ E_P       :: EqualRT
                E + S ⇌ E_S       :: EqualRT
                E_S <--> E_P      :: EqualRT
            end
        end
    end
    spec = allosteric_spec_from_mechanism_and_rxn(m_seed, rxn)
    result = EnzymeRates._expand_change_allo_state(spec, rxn)

    # 1. count: 3 group_tags + 2 reg_ligand_tags = 5 variants.
    @test length(result) == 5

    # 2. Δ params: 3 group relaxations ([1, 1, 2]) + 2 ligand relaxations
    # (R1 :OnlyR → :NonequalRT = +1; R2 :OnlyT → :NonequalRT = +1) → [1, 1, 1, 1, 2].
    deltas = sort([r.n_fit_params_estimate -
                   spec.n_fit_params_estimate for r in result])
    @test deltas == [1, 1, 1, 1, 2]

    # 3. property: exactly one variant has each ligand independently relaxed.
    n_r1_relaxed_only = count(r -> r.reg_ligand_tags[:R1] == :NonequalRT &&
                                    r.reg_ligand_tags[:R2] == :OnlyT, result)
    n_r2_relaxed_only = count(r -> r.reg_ligand_tags[:R1] == :OnlyR &&
                                    r.reg_ligand_tags[:R2] == :NonequalRT, result)
    @test n_r1_relaxed_only == 1
    @test n_r2_relaxed_only == 1
end
```

- [ ] **Step 9.6: Run tests and commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test: add coverage for under-tested mechanism shapes

Five gaps identified by post-rewrite review:
- Ter-ter sequential as input to _expand_re_to_ss (was only tested in
  topology generation).
- Non-default site multiplicities (catalytic_n=4 with reg site
  multiplicity=2) verified across allo moves.
- Cross-type allosteric seed for _expand_re_to_ss substrate-as-dead-end-I
  overlap (was previously plain-only).
- Multi-regulator-ligand seed for _expand_change_allo_state with
  independent R1::OnlyR + R2::OnlyT relaxations.
- Multi-product iso shape investigation (covered by existing seeds; see
  comment in relevant testset).

Surfaced by post-rewrite review (findings #13, #14, #20, #21, #22).
EOF
)"
```

---

## Task 10: Strengthen `init_mechanisms` tests

**Files:** `test/test_mechanism_enumeration.jl`

**Findings addressed:** #15 (cap=30 may skip exotic specs), #16 (same-metabolite-RE-bindings test passes vacuously when loop body empty), #18 (`_assert_spec_invariants` too weak for MechanismSpec), #19 (`compile_mechanism` not invoked by name).

- [ ] **Step 10.1: Add sentinel counter to "Same-metabolite RE bindings share kinetic_group"**

Locate the testset (around line 1245). After the existing loop, add an assertion that the inner loop body executed at least once:

```julia
# Sentinel: the test would silently pass if no spec had any metabolite
# with 2+ RE binding steps. Verify that at least one such case was
# actually checked.
n_assertions_fired = 0
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
        n_assertions_fired += 1
    end
end
@test n_assertions_fired >= 1   # at least one multi-binding case existed
```

(Replace the existing loop body — don't duplicate the inner assertions.)

- [ ] **Step 10.2: Strengthen `_assert_spec_invariants` for `MechanismSpec`**

Locate the helper definition (around line 287). Add structural invariants:

```julia
function _assert_spec_invariants(spec::MechanismSpec)
    @test spec.n_fit_params_estimate >= 0
    # Structural invariants every valid MechanismSpec should satisfy:
    for s in spec.steps
        @test !isempty(s.reactants)
        @test !isempty(s.products)
        @test s.kinetic_group >= 1
        # The "from form" (first reactant) should differ from the
        # "to form" (first product) — a step without form change is degenerate.
        @test s.reactants[1] != s.products[1]
    end
end
```

- [ ] **Step 10.3: Add `compile_mechanism` dispatch test**

In the existing `compile_mechanism round-trip` testset (around line 1100), add explicit dispatch tests:

```julia
# Dispatch: compile_mechanism dispatches correctly to the right type-specific
# constructor. Verify by comparing to the explicit constructor.
@test EnzymeRates.compile_mechanism(spec_uu) === EnzymeMechanism(spec_uu)
@test EnzymeRates.compile_mechanism(spec_seq) === EnzymeMechanism(spec_seq)

# Also verify dispatch for AllostericMechanismSpec via a small allosteric
# seed (this overlaps with the dedicated allosteric round-trip testset
# but is part of the dispatch contract).
m_allo_dispatch = @allosteric_mechanism begin
    substrates: S
    products: P
    site(:catalytic, 2): begin
        steps: begin
            E + P ⇌ E_P    :: EqualRT
            E + S ⇌ E_S    :: EqualRT
            E_S <--> E_P   :: EqualRT
        end
    end
end
spec_allo_dispatch = allosteric_spec_from_mechanism_and_rxn(m_allo_dispatch, uni_uni_allo)
@test EnzymeRates.compile_mechanism(spec_allo_dispatch) ===
    AllostericEnzymeMechanism(spec_allo_dispatch)
```

- [ ] **Step 10.4: Run tests and commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test: strengthen init_mechanisms tests against silent passes

- Add sentinel counter to "Same-metabolite RE bindings share
  kinetic_group" — previously the test passed if no spec had any
  metabolite with 2+ RE bindings (loop body skipped). Sentinel asserts
  the assertion fired at least once.
- Strengthen `_assert_spec_invariants` for MechanismSpec from a single
  `n_fit_params_estimate >= 0` check (effectively impossible to violate)
  to structural invariants on every step (non-empty reactants/products,
  kinetic_group >= 1, distinct from/to forms).
- Add explicit `compile_mechanism === Constructor(spec)` dispatch
  assertions so a regression that broke the dispatch path would surface
  immediately rather than via downstream comparison failure.

Surfaced by post-rewrite review (findings #16, #18, #19).

Note: cap=30 in the upper-bound testset (#15) is unchanged; the cap
exists to bound @generated compile cost, and reducing it would risk
timeout. Documenting as accepted.
EOF
)"
```

---

## Task 11: Add `_competition_patterns` boundary cases

**Files:** `test/test_mechanism_enumeration.jl`

**Findings addressed:** #17 (boundary cases for `_competition_patterns` and `_inhibitor_competition_patterns` missing).

- [ ] **Step 11.1: Add 2×3, 3×2 cases for `_competition_patterns`**

Locate the testset (around line 686). Add:

```julia
@testset "Asymmetric: 2 × 3" begin
    # 2 substrates × 3 products. Counts are derivable from the inclusion-
    # exclusion formula; the actual count is covered by the source's
    # generation logic. Just verify it produces non-empty results and
    # all patterns cover all metabolites.
    pats = EnzymeRates._competition_patterns(Set([:A, :B]), Set([:P, :Q, :R]))
    @test !isempty(pats)
    for pat in pats
        for s in [:A, :B]
            @test any((s, p) in pat for p in [:P, :Q, :R])
        end
        for p in [:P, :Q, :R]
            @test any((s, p) in pat for s in [:A, :B])
        end
    end
end

@testset "Asymmetric: 3 × 2" begin
    pats = EnzymeRates._competition_patterns(Set([:A, :B, :D]), Set([:P, :Q]))
    @test !isempty(pats)
    for pat in pats
        for s in [:A, :B, :D]
            @test any((s, p) in pat for p in [:P, :Q])
        end
        for p in [:P, :Q]
            @test any((s, p) in pat for s in [:A, :B, :D])
        end
    end
end
```

- [ ] **Step 11.2: Add 3-existing-inhibitor case for `_inhibitor_competition_patterns`**

Locate the testset (around line 732):

```julia
@testset "Bi-bi with 3 existing inhibitors: 9 × 8 = 72" begin
    # 3 existing inhibitors → 2^3 = 8 inhibitor-competition combinations.
    # Combined with 9 base patterns → 9 × 8 = 72 variants.
    pats = EnzymeRates._inhibitor_competition_patterns(
        Set([:A, :B]), Set([:P, :Q]),
        [:I1__reg, :I2__reg, :I3__reg])
    @test length(pats) == 72
end
```

- [ ] **Step 11.3: Run tests and commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test: add boundary cases for _competition_patterns and _inhibitor_*

_competition_patterns was tested for 1×1, 1×2, 2×1, 2×2, 3×3 only.
Added 2×3 and 3×2 cases (asymmetric multi). _inhibitor_competition_patterns
was tested for 0/1/2 existing inhibitors; added the 3-inhibitor case
(72 variants).

Surfaced by post-rewrite review (finding #17).
EOF
)"
```

---

## Task 12: Quick cleanups

**Files:** `test/test_mechanism_enumeration.jl`

**Findings addressed:** #23 (RE→SS NonequalRT preservation incomplete), #24 (residual `first(init_mechanisms)` at line 3075), #25 (inter-move dedup interaction untested).

- [ ] **Step 12.1: Add missing preservation assertions to RE→SS `:NonequalRT` testset**

Locate the testset (around line 1689). The :EqualRT and :OnlyR siblings have:

```julia
@test r.allosteric_multiplicities == spec.allosteric_multiplicities
@test r.reg_ligand_tags == spec.reg_ligand_tags
@test r.base.reaction === spec.base.reaction
```

Add the same three assertions to the `:NonequalRT` testset's preservation block.

- [ ] **Step 12.2: Replace `first(init_mechanisms)` at line 3075**

Locate the line. The "MechanismSpec → empty" testset for `_expand_change_allo_state` uses `init_mechanisms(uni_uni_allo) |> first`. Replace with a literal seed:

```julia
@testset "MechanismSpec → empty (negative)" begin
    m_seed = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ E_P
            E + S ⇌ E_S
            E_S <--> E_P
        end
    end
    spec = mechanism_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo)
    result = EnzymeRates._expand_change_allo_state(spec, uni_uni_allo)
    @test isempty(result)
    @test result isa Vector{AllostericMechanismSpec}
end
```

- [ ] **Step 12.3: Add inter-move dedup interaction test**

In the `dedup!` testset (around line 3122), add a sub-testset that exercises inter-move overlap:

```julia
@testset "Inter-move overlap collapses via dedup!" begin
    # Two different expansion paths from the same seed can produce the
    # same target spec. Verify dedup! collapses such duplicates after
    # canonicalization.
    m_seed = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ E_P
            E + S ⇌ E_S
            E_S <--> E_P
        end
    end
    spec = mechanism_spec_from_mechanism_and_rxn(m_seed, uni_uni_rxn)

    # Path 1: RE→SS on the S-binding group.
    re_to_ss_results = EnzymeRates._expand_re_to_ss(spec)
    # Path 2: split the iso group (no-op since it's singleton, so empty)
    # plus another move that produces a structurally identical result.
    # In practice, expand_mechanisms calls all moves on the same spec
    # and dedup! collapses inter-move overlaps. Verify by running
    # expand_mechanisms and checking that the result count is no larger
    # than the union of unique compiled mechanisms.
    expanded = EnzymeRates.expand_mechanisms([spec], uni_uni_rxn)
    EnzymeRates.dedup!(expanded)
    for (pc, specs) in expanded
        compiled = Set(EnzymeRates.compile_mechanism(s) for s in specs)
        # After dedup, compiled count should equal spec count
        # (no duplicates surviving).
        @test length(compiled) == length(specs)
    end
end
```

- [ ] **Step 12.4: Run tests and commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
test: quick cleanups from post-rewrite review

- Add missing allosteric_multiplicities, reg_ligand_tags, and
  base.reaction preservation assertions to the RE→SS :NonequalRT
  testset (the :EqualRT and :OnlyR siblings had them).
- Replace residual `init_mechanisms(...) |> first` at the
  _expand_change_allo_state "MechanismSpec → empty" testset with a
  literal @enzyme_mechanism seed.
- Add inter-move dedup-interaction test verifying that expand_mechanisms
  + dedup! produces no compiled-mechanism duplicates within a param-count
  bucket.

Surfaced by post-rewrite review (findings #23, #24, #25).
EOF
)"
```

---

## Self-Review Checklist

After all 12 tasks complete:

- [ ] **All 25 findings addressed:** spot-check that each finding number is referenced in at least one commit message.

- [ ] **Test count grew, didn't regress:** baseline 25593; after all tasks expect ~25750-25800 (depending on how many new sub-testsets land).

- [ ] **No vacuous tests introduced:** new tests should have meaningful assertions, not `@test !isempty(result)` as their only check.

- [ ] **Anti-pattern grep clean:** `grep -nE "init_mechanisms\(.*\)[[:space:]]*\|>[[:space:]]*first|first\(EnzymeRates\.init_mechanisms\(|first\(filter\(.*init_mechanisms" test/test_mechanism_enumeration.jl` returns at most the integration-test sites where it's intentional.

- [ ] **Final test run green:** `julia --project -e 'using Pkg; Pkg.test()'` passes.

- [ ] **Commit count check:** ~12 commits between this plan's first task and the final task. Each is reviewable independently.
