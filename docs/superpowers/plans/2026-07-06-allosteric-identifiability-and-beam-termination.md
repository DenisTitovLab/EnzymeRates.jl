# Allosteric identifiability + beam termination — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix three defects from the v0.1.5 LDH run — beam non-termination (split/canonicalize self-loop), a Haldane-dependent parameter leaking into `fitted_params`, and enumeration of empirically-indistinguishable allosteric mechanisms.

**Architecture:** All three are fixed at the source. Bug 1: the split move canonicalizes its own output and drops no-op self-loops, which makes the PASS-1 canonicalization redundant (removed) and relocates the canonicalization helpers next to the split move. Bug 2: the allosteric independent-parameter set becomes a uniform complement of the dependent set. Bug 3: `_expand_to_allosteric` emits only mechanisms whose conformational split the data can resolve.

**Tech Stack:** Julia 1.12, package `EnzymeRates` (dev at `/home/denis.linux/.julia/dev/EnzymeRates`). Tests via `julia --project -e 'using Pkg; Pkg.test()'` or a single file via `include`.

**Design spec:** `docs/superpowers/specs/2026-07-06-allosteric-identifiability-and-beam-termination-design.md` (read it first).

## Global Constraints

- 92-character line limit, 4-space indent (`.claude/CLAUDE.md`).
- `rate_equation` must stay allocation-free and sub-120 ns/call (`test/test_rate_eq_derivation.jl::test_rate_equation_performance`). None of these changes touch the rate-equation hot path, but the perf test must stay green.
- Canonical Step Form is load-bearing: `fitted_params` order and `eq_hash` depend on step/group order. Bug 2 must preserve `merged_indep` ordering.
- Parameter-name rendering must flow through the `name(p, m)` chokepoint (AST-walker guard, `test/test_types.jl`).
- TDD: write the failing test, run it red, implement minimally, run it green, commit. Run the full suite before each commit that changes behavior.
- Run the full suite: `cd /home/denis.linux/.julia/dev/EnzymeRates && julia --project -e 'using Pkg; Pkg.test()'`.
- All new code files start with two `# ABOUTME:` lines. Match surrounding style; no temporal/"new"/"old" comments.

---

## Bug 1 — beam non-termination

### Task 1.1: Relocate canonicalization helpers to the enumeration file

Pure move — no behavior change. `_canonical_mechanism` and `_merge_tied_kinetic_groups` currently live in `src/rate_eq_derivation.jl`; move them beside `_expand_split_kinetic_group` in `src/mechanism_enumeration.jl`. They keep calling `_build_wegscheider_rename_map` (stays in `rate_eq_derivation.jl`; same module, no import needed).

**Files:**
- Modify: `src/rate_eq_derivation.jl` (remove `_merge_tied_kinetic_groups` at `:1877` and `:1923`, `_canonical_mechanism` at `:1982` and `:1984`, plus their docstrings)
- Modify: `src/mechanism_enumeration.jl` (paste them near `_expand_split_kinetic_group`, ~`:1249`)

**Interfaces:**
- Produces: `_canonical_mechanism(m::Mechanism)::Mechanism`, `_canonical_mechanism(am::AllostericMechanism)::AllostericMechanism` (unchanged signatures/behavior).

- [ ] **Step 1: Confirm the exact function spans.** Run:
  `grep -n 'function _merge_tied_kinetic_groups\|^_canonical_mechanism\|function _canonical_mechanism' src/rate_eq_derivation.jl`
  Read each function body plus its leading docstring so you move the whole unit.
- [ ] **Step 2: Cut the four function definitions (both `_merge_tied_kinetic_groups` methods, both `_canonical_mechanism` methods) with their docstrings out of `src/rate_eq_derivation.jl`.**
- [ ] **Step 3: Paste them into `src/mechanism_enumeration.jl` immediately after the two `_expand_split_kinetic_group` methods.** Keep the code identical.
- [ ] **Step 4: Verify no other definition or stale reference remains.** Run:
  `grep -rn 'function _canonical_mechanism\|function _merge_tied_kinetic_groups' src/` — expect the two functions only in `mechanism_enumeration.jl`.
- [ ] **Step 5: Run the full suite** (`julia --project -e 'using Pkg; Pkg.test()'`). Expected: all green — a pure move changes nothing.
- [ ] **Step 6: Commit.**
  `git add src/rate_eq_derivation.jl src/mechanism_enumeration.jl && git commit -m "Move canonicalization helpers next to the split move"`

### Task 1.2: Split move canonicalizes + drops self-loops; add the canonicality invariant

The split move must return canonical, self-loop-free output. Add a test that every `expand_mechanisms` child is canonical, watch it fail (the split move currently emits non-canonical no-op splits), then fix the split move.

**Files:**
- Test: `test/test_mechanism_enumeration.jl` (add a testset)
- Modify: `src/mechanism_enumeration.jl` (`_expand_split_kinetic_group`, both methods, `:1249`/`:1260`)

**Interfaces:**
- Consumes: `_canonical_mechanism` (Task 1.1), `expand_mechanisms(parents::Vector, rxn)`, `init_mechanisms(rxn)`.

- [ ] **Step 1: Write the failing invariant test.** Add to `test/test_mechanism_enumeration.jl`:

```julia
@testset "expand_mechanisms output is canonical" begin
    rxn = @enzyme_reaction begin
        substrates: NADH[C21H29N7O14P2], Pyruvate[C3H4O3]
        products: Lactate[C3H6O3], NAD[C21H27N7O14P2]
        oligomeric_state: 4
    end
    MECH = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}
    base = collect(EnzymeRates.init_mechanisms(rxn))
    children = EnzymeRates.expand_mechanisms(MECH[base...], rxn)
    noncanon = [c for c in children if EnzymeRates._canonical_mechanism(c) != c]
    @test isempty(noncanon)
end
```

- [ ] **Step 2: Run it, expect FAIL** (thousands of non-canonical split children).
  Run: `julia --project -e 'using Pkg; Pkg.activate("."); include("test/test_mechanism_enumeration.jl")'` (or the whole suite). Expected: this testset fails.
- [ ] **Step 3: Fix both `_expand_split_kinetic_group` methods** to canonicalize each split and drop self-loops. Replace the `Mechanism` method:

```julia
function _expand_split_kinetic_group(m::Mechanism)
    results = Mechanism[]
    mc = _canonical_mechanism(m)
    for g in kinetic_groups(m)
        length(steps(m)[g]) >= 2 || continue
        for split_idx in 1:length(steps(m)[g])
            child = _canonical_mechanism(
                _with_steps(m, _split_one_step(steps(m), g, split_idx)))
            child == mc || push!(results, child)
        end
    end
    results
end
```

  And the `AllostericMechanism` method:

```julia
function _expand_split_kinetic_group(am::AllostericMechanism)
    results = AllostericMechanism[]
    mc = _canonical_mechanism(am)
    for g in kinetic_groups(am)
        length(steps(am)[g]) >= 2 || continue
        for split_idx in 1:length(steps(am)[g])
            new_groups = _split_one_step(steps(am), g, split_idx)
            new_states = vcat(cat_allo_states(am), [cat_allo_states(am)[g]])
            child = _canonical_mechanism(
                _with_steps_and_cat_states(am, new_groups, new_states))
            child == mc || push!(results, child)
        end
    end
    results
end
```

- [ ] **Step 4: Run the invariant test, expect PASS.**
- [ ] **Step 5: Run the full suite.** Expected: green. (If a pre-existing enumeration-count test now sees fewer duplicate mechanisms, that is the intended effect; adjust only counts that were asserting the now-removed no-op splits, and note it in the commit.)
- [ ] **Step 6: Commit.**
  `git add -u && git commit -m "Split move drops Wegscheider-tied no-op splits (fixes beam self-loop)"`

### Task 1.3: Remove the now-redundant PASS-1 canonicalization

With every `expand_mechanisms` child canonical and `init_mechanisms` canonical, PASS 1 no longer needs to canonicalize.

**Files:**
- Modify: `src/identify_rate_equation.jl` (`_process_batch` PASS 1, ~`:514`)

- [ ] **Step 1: Read `_process_batch` PASS 1.** The line is `m = _canonical_mechanism(m0)`, with `(mech = m, orig = m0, …)` below.
- [ ] **Step 2: Replace `m = _canonical_mechanism(m0)` so PASS 1 uses `m0` directly.** Set `mech = m0` in the record (drop the separate `m`; `orig = m0` stays). Keep everything else — `eq_hash`, the memo, failure routing — untouched.
- [ ] **Step 3: Confirm the only remaining caller of `_canonical_mechanism` is the split move.** Run:
  `grep -rn '_canonical_mechanism' src/ | grep -v 'function _canonical_mechanism'` — expect only the two `_expand_split_kinetic_group` methods.
- [ ] **Step 4: Run the full suite.** Expected: green.
- [ ] **Step 5: Add a wide-beam termination regression** to `test/test_identify_rate_equation.jl`. Use synthetic data (quality irrelevant), a wide beam, and a bounded `max_param_count`; assert it returns rather than hangs. Model it on the existing smoke test but with `min_beam_width` well above 1:

```julia
@testset "wide beam terminates (no split self-loop)" begin
    rxn = @enzyme_reaction begin
        substrates: NADH[C21H29N7O14P2], Pyruvate[C3H4O3]
        products: Lactate[C3H6O3], NAD[C21H27N7O14P2]
        oligomeric_state: 4
    end
    n = 30; rnd() = 0.1 + 9.9 * rand()
    data = (group=["G$(mod(i,5))" for i in 1:n], Rate=[rnd() for _ in 1:n],
            NADH=[rnd() for _ in 1:n], Pyruvate=[rnd() for _ in 1:n],
            Lactate=[rnd() for _ in 1:n], NAD=[rnd() for _ in 1:n])
    prob = IdentifyRateEquationProblem(rxn, data; Keq=20000.0, scale_k_to_kcat=1.0)
    results = identify_rate_equation(prob; optimizer=cmaes_opt,
        min_beam_width=20, loss_rel_threshold=1.3, loss_abs_threshold=0.001,
        loss_parsimony_threshold=1.01, max_param_count=7,
        n_restarts=1, maxtime=0.1, save_dir=mktempdir())
    @test results.best !== nothing   # returned at all == terminated
end
```

  (`cmaes_opt` is already defined at the top of that file.)
- [ ] **Step 6: Run that test, expect PASS** (terminates). Sanity: it would hang on the pre-fix code.
- [ ] **Step 7: Commit.**
  `git add -u && git commit -m "Drop redundant PASS-1 canonicalization; add wide-beam termination test"`

---

## Bug 2 — Haldane-dependent parameter in `fitted_params`

### Task 2.1: Uniform `∉ keys(dep)` filter, driven by a promoted LDH fixture

The current spec suite is clean (0 dep∩indep overlaps across all 39 specs), so the fix must be TDD'd by first promoting a leak-triggering LDH mechanism into the suite. Fixture-add, fix, and golden regen are coupled, so this is one task.

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (add leak-triggering LDH mechanisms to `MECHANISM_TEST_SPECS`)
- Modify: `test/test_rate_eq_derivation.jl` (add the invariant testset)
- Modify: `src/rate_eq_derivation.jl` (`_dependent_param_exprs(::Type{AllostericEnzymeMechanism})`, ~`:1555`–`:1605`)
- Modify: `test/reference/allosteric_golden_reference.txt` (regenerate)

**Interfaces:**
- Consumes: `_dependent_param_exprs(M)::Tuple{dep::Dict, indep::Tuple}` for a compiled allosteric mechanism type.

- [ ] **Step 1: Read `LDH_ISTATE_FAILURE_MECHS`** (`test/test_rate_eq_derivation.jl`, ~`:1114`) and the caseB mechanisms there. Identify the ones that trigger the leak (shared-`:EqualAI` catalytic reverse-rate with a dead I-cycle — the "caseB_reverse"/"caseB_binding" shapes). Read `test/mechanism_definitions_for_test_enzyme_derivation.jl` to learn the exact `MECHANISM_TEST_SPECS` entry struct/fields (name, mechanism, expected counts).
- [ ] **Step 2: Add those mechanisms as new `MECHANISM_TEST_SPECS` entries**, matching the existing entry fields. Give them descriptive names (e.g. `m_ldh_caseB_reverse`).
- [ ] **Step 3: Write the failing invariant test** in `test/test_rate_eq_derivation.jl`:

```julia
@testset "indep ∩ keys(dep) == ∅ for all specs" begin
    for spec in MECHANISM_TEST_SPECS
        M = typeof(EnzymeRates.compile_mechanism(spec.mechanism))
        dep, indep = EnzymeRates._dependent_param_exprs(M)
        @test isempty(intersect(Set(keys(dep)), Set(indep)))
    end
end
```

  (Use the actual mechanism field name from Step 1 in place of `spec.mechanism`.)
- [ ] **Step 4: Run it, expect FAIL** on the newly-added LDH mechanisms (the shared reverse rate is in both `dep` and `indep`). The `test_allosteric_golden` test will also fail (new mechanisms not yet in the golden) — expected; Step 8 regenerates it.
- [ ] **Step 5: Fix `_dependent_param_exprs(::Type{AllostericEnzymeMechanism})`.** Read `:1555`–`:1605`. `dep` is fully assembled by `:1601`; `merged_indep` is built at `:1603` as `(indep_A..., indep_I_list..., reg_params_a..., reg_params_i_indep..., :L)`. Replace the return of `merged_indep` with an ordered complement: `Tuple(p for p in merged_indep if p ∉ keys(dep))`. You may also simplify `indep_I_list`'s filter (`:1583`) by removing its now-redundant conditions, but keep `p ∈ S_I` and `p ∉ collapse_targets` (membership gates, not dep-exclusion). The uniform filter is what closes the leak.
- [ ] **Step 6: Run the invariant test, expect PASS.**
- [ ] **Step 7: Run the full suite except the golden.** Confirm every test except `test_allosteric_golden` is green. Pre-existing specs' `indep` must be unchanged (they had 0 overlaps, so the complement leaves them identical) — if a pre-existing spec's non-golden test changes, STOP and inspect.
- [ ] **Step 8: Regenerate the golden.** Grep `test/` for how `allosteric_golden_reference.txt` is produced (a regen flag/script). Regenerate, then DIFF against the committed version: pre-existing mechanisms must be byte-identical; the only additions are the new LDH mechanisms, whose `fitted_params` must NOT contain the reverse-rate symbol that also appears in a Haldane line. If a pre-existing mechanism's golden changed, STOP.
- [ ] **Step 9: Run the full suite, expect green.**
- [ ] **Step 10: Commit.**
  `git add -u && git commit -m "Allosteric fitted_params = ordered complement of dependent set; add LDH i-state specs"`

---

## Bug 3 — enumeration emits indistinguishable allosteric mechanisms

### Task 3.1: Stop emitting the all-`:EqualAI` baseline

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_expand_to_allosteric(m::Mechanism, rxn)`, `:1554`)
- Modify: `test/test_mechanism_enumeration.jl` (the assertions that require the all-`:EqualAI` mechanism)

- [ ] **Step 1: Write a failing test** that `_expand_to_allosteric` emits no all-`:EqualAI` mechanism:

```julia
@testset "_expand_to_allosteric emits no all-EqualAI baseline" begin
    rxn = @enzyme_reaction begin
        substrates: NADH[C21H29N7O14P2], Pyruvate[C3H4O3]
        products: Lactate[C3H6O3], NAD[C21H27N7O14P2]
        oligomeric_state: 4
    end
    base = first(EnzymeRates.init_mechanisms(rxn))
    allo = EnzymeRates._expand_to_allosteric(base, rxn)
    @test all(am -> !all(==(:EqualAI), EnzymeRates.cat_allo_states(am)), allo)
end
```

- [ ] **Step 2: Run it, expect FAIL** (the baseline at `:1560` is all-`:EqualAI`).
- [ ] **Step 3: Remove the all-`:EqualAI` baseline emission** — delete the `push!(results, ...base_tags...)` at `:1560`-`:1562`, keeping the per-group-`:OnlyA` loop.
- [ ] **Step 4: Run the new test, expect PASS.**
- [ ] **Step 5: Update the enumeration tests that assert the baseline exists** — `test/test_mechanism_enumeration.jl:2929`, `:2968` (they require an all-`:EqualAI` mechanism). Re-read them; change them to assert the baseline is NOT produced, and update the per-multiplicity counts. Update the Δ-count assertions at `:2837`/`:2892` (`[1,1,2,2,2,2]`) — the two Δ=1 entries were the degenerate baseline + catalysis-OnlyA-all-EqualAI; recompute the expected deltas for the reduced emission set.
- [ ] **Step 6: Run the full suite, expect green.**
- [ ] **Step 7: Commit.**
  `git add -u && git commit -m "Enumeration: drop degenerate all-EqualAI allosteric baseline"`

### Task 3.2: Catalysis-`:OnlyA` only as a V-type (paired with a regulator)

The bare catalytic-step-`:OnlyA` variant (all bindings `:EqualAI`) is degenerate (`L` folds into `kcat`); it is identifiable only with an `:OnlyA`/`:OnlyI` regulator. Emit the combined V-type instead of the bare one, and only when the reaction declares regulators.

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_expand_to_allosteric`, and study `_expand_add_allosteric_regulator` at `:1608` for how a regulator site is attached)
- Test: `test/test_mechanism_enumeration.jl`

**Interfaces:**
- Consumes: `regulators(rxn)`, `RegulatorySite`, `AllostericRegulator`, the `AllostericMechanism` constructor `AllostericMechanism(reaction, steps, cat_allo_states, cat_mult, regulatory_sites)` (`src/types.jl:597`).

- [ ] **Step 1: Determine which kinetic group is the catalytic step.** Read `kinetic_groups`/`rep_step`/`bound_metabolite` usage in `test/test_identify_rate_equation.jl:28-35` — a group whose representative step binds no `Reactant` (i.e., an isomerization/conversion step) is catalytic; binding groups bind a substrate/product. Write a helper predicate or reuse an existing one.
- [ ] **Step 2: Write a failing test** capturing the rule for a reaction WITH a declared regulator and one WITHOUT:

```julia
@testset "catalysis-OnlyA is V-type only (needs a regulator)" begin
    # regulator-free: no bare catalysis-OnlyA emitted
    rxn0 = @enzyme_reaction begin
        substrates: S[C]; products: P[C]
    end
    base0 = first(EnzymeRates.init_mechanisms(rxn0))
    allo0 = EnzymeRates._expand_to_allosteric(base0, rxn0)
    cat_g(am) = findfirst(g -> EnzymeRates.bound_metabolite(
        EnzymeRates.rep_step(am, g)) === nothing, 1:length(EnzymeRates.steps(am)))
    for am in allo0
        g = cat_g(am)
        if g !== nothing && EnzymeRates.cat_allo_states(am)[g] == :OnlyA
            @test !isempty(EnzymeRates.regulatory_sites(am))
        end
    end
    # regulator-declared: a V-type (catalysis-OnlyA + regulator) IS reachable
    rxn1 = @enzyme_reaction begin
        substrates: S[C]; products: P[C]
        competitive_inhibitors: R
    end
    base1 = first(EnzymeRates.init_mechanisms(rxn1))
    allo1 = EnzymeRates._expand_to_allosteric(base1, rxn1)
    @test any(am -> begin g = cat_g(am);
        g !== nothing && EnzymeRates.cat_allo_states(am)[g] == :OnlyA &&
        !isempty(EnzymeRates.regulatory_sites(am)) end, allo1)
end
```

  (Verify `bound_metabolite` of a catalytic step's rep returns `nothing`/non-`Reactant`; adjust the predicate to the actual API if needed. Verify `rep_step`/`bound_metabolite` names — the identify test uses `bound_metabolite(rep)` and `Reactant`.)
- [ ] **Step 3: Run it, expect FAIL.**
- [ ] **Step 4: Implement the rule in `_expand_to_allosteric(m::Mechanism, rxn)`:** for each group flipped to `:OnlyA`, if that group is a **binding** step, emit it bare (K-type, as today); if it is the **catalytic** step, emit it only combined with each declared regulator ligand as an `:OnlyA` or `:OnlyI` regulatory site (V-type), and emit nothing for it when `regulators(rxn)` is empty. Reuse `_expand_add_allosteric_regulator`'s site-construction pattern.
- [ ] **Step 5: Run the test, expect PASS.**
- [ ] **Step 6: Run the full suite.** Update any enumeration-count tests affected by the catalysis-OnlyA change; confirm distinguishable K-type/regulator paths are unchanged.
- [ ] **Step 7: Commit.**
  `git add -u && git commit -m "Enumeration: catalysis-OnlyA emitted only as V-type (with a regulator)"`

### Task 3.3: Confirm change-allo-state degeneracy handling and final sweep

- [ ] **Step 1: Re-read `test/test_allosteric_collapse.jl`** and `_expand_change_allo_state`. Confirm the single-binding-`:NonequalAI` `K_I = K_A` collapse (Δ=0) is still covered and consistent with the distinguishability principle. Add a one-line comment cross-referencing the principle if helpful; add a test only if a gap exists.
- [ ] **Step 2: Run the full suite** (`julia --project -e 'using Pkg; Pkg.test()'`). Expected: fully green.
- [ ] **Step 3: Commit** any test/comment additions.
  `git add -u && git commit -m "Confirm change-allo-state degeneracy consistent with distinguishability rule"`

---

## Final verification

- [ ] Run the full suite once more, clean: `cd /home/denis.linux/.julia/dev/EnzymeRates && julia --project -e 'using Pkg; Pkg.test()'`. All green.
- [ ] Re-run the delta-0 termination sanity from the spec if desired (`canon(child) != parent` holds for all split output; the invariant test from Task 1.2 already enforces it).
- [ ] Leave the branch `fix-allosteric-identifiability-and-beam-termination` for Denis to review; do not merge.
