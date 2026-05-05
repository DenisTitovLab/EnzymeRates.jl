# Mechanism enumeration test rewrite — design

## 1. Overview

Target file: `test/test_mechanism_enumeration.jl` (~2451 lines).

What this rewrite delivers:

- Every spec-consuming testset seeds from a literal `@enzyme_mechanism` or
  `@allosteric_mechanism` block. No `init_mechanisms |> first`, no
  `first(filter(specs) do …)`, no implicit "whatever variant comes out first"
  picks.
- Every per-move testset follows a fixed 7-item checklist (count, Δ params,
  compilability, structural change, preservation, negative case, plus
  cross-type for the polymorphic moves).
- The file is reorganized into pipeline-execution order: support functions →
  initialization → base-spec moves → allosteric moves → composition →
  integration.
- Equivalence-style structural assertions (`Set(compile_mechanism.(results))
  == Set(<expected mechanism literals>)`) for moves with ≤6 outputs OR where
  existing comments enumerate the expected results. Property-style for larger
  fan-outs.

What this rewrite does NOT do:

- Add new public API.
- Change `src/` (except for bugs uncovered during the rewrite — see §6).
- Add new pipeline behavior.

The only test-infrastructure addition is one helper, `allosteric_spec_from_mechanism_and_rxn`,
living alongside the existing `mechanism_spec_from_mechanism_and_rxn` at the top of the
test file. No new files.

## 2. Test infrastructure

Two helpers at the top of `test_mechanism_enumeration.jl`:

**Existing — keep (renamed):**

```julia
mechanism_spec_from_mechanism_and_rxn(m::EnzymeMechanism, rxn) → MechanismSpec
```

Builds a `MechanismSpec` from a compiled `EnzymeMechanism` (built via
`@enzyme_mechanism`) and a reaction. The helper takes BOTH inputs because
`EnzymeMechanism` carries no reaction-level metadata (atoms, regulator
declarations) — the reaction is a separate dual input. The helper validates
internally that they're consistent: substrates and products names must match
exactly; the mechanism's regulators must be a subset of the reaction's
declared regulators (so the test can seed a mechanism that doesn't yet bind
every declared regulator). Throws a descriptive error otherwise.

Per-call-site `EnzymeMechanism(spec) === m_seed` round-trip assertions are
NOT required — the helper validates internally. Round-trip equality is
asserted in a single dedicated testset (Task 1).

**New — add:**

```julia
allosteric_spec_from_mechanism_and_rxn(m::AllostericEnzymeMechanism, rxn)
    → AllostericMechanismSpec
```

Symmetric helper. Same internal-consistency rules plus
`oligomeric_state(rxn) == catalytic_multiplicity(m)`. Implementation:

1. Build base `MechanismSpec` from `catalytic_mechanism(m)` via the existing
   helper.
2. Extract `catalytic_n` from `catalytic_multiplicity(m)`.
3. Extract `group_tags::Dict{Int,Symbol}` by iterating kinetic groups (1 to
   n_groups) and reading `cat_allo_state(m, g)` for each. All entries —
   including `:NonequalRT` — are stored explicitly (dense storage).
4. Extract reg sites by iterating `regulatory_site_ligands(m, i)` /
   `regulatory_site_multiplicity(m, i)` / `reg_allo_state(m, i, lig)`. Every
   ligand tag is stored explicitly (dense storage).
5. Set `n_fit_params_estimate = length(fitted_params(m))`.

Storage is dense in both `AllostericMechanismSpec` and the compiled
`AllostericEnzymeMechanism`, so the round-trip is straightforward
pass-through with no filtering.

**Note on `n_fit_params_estimate` semantics:** the helper sets the field to
`length(fitted_params(m))` — the *exact* fitted-param count for the
compiled mechanism. By contrast, `init_mechanisms` sets it to an
*upper-bound estimate* (via `_apply_equivalence_grouping`) that can be
strictly greater than the actual count when mirror cycles exist. This
asymmetry is benign for the rewrite because all per-move delta assertions
(`r.n_fit_params_estimate == spec.n_fit_params_estimate + delta`) are
baseline-independent — every move computes its delta from tag/RE-vs-SS
properties of the affected group, not from the baseline value. The
upper-bound invariant `n_fit_params_estimate >= length(fitted_params(m))`
is exercised separately under §4 (`init_mechanisms` testset, init-seeded);
helper-seeded specs satisfy it trivially via equality.

**Round-trip validation pattern (Task 1's dedicated testset only):**

```julia
m_allo = @allosteric_mechanism begin … end
spec = allosteric_spec_from_mechanism_and_rxn(m_allo, rxn)
@test AllostericEnzymeMechanism(spec) === m_allo  # round-trip lossless
```

The `===` check catches drift between the macro, the spec, and the compiler.
This assertion is NOT replicated at every call site — the helper's internal
consistency check covers the dual-input mismatch case, and the round-trip
equality lives in a single dedicated testset where it IS the thing being
tested.

## 3. Standard test checklist

Every per-move testset follows this template. Items 1–6 are mandatory; item 7
applies only to the three polymorphic moves (`_expand_re_to_ss`,
`_expand_split_kinetic_group`, `_expand_add_dead_end_regulator`).

**The independent-derivation rule (non-negotiable).** Every numerical or
structural prediction in a test — the count `N` (item 1), the
`EXPECTED_DELTA` (item 2), the contents of the equivalence set (item 4) —
MUST be accompanied by a comment that derives the predicted value from
**independent reasoning about the seed mechanism and the move's
specification**, not from observing what the code produces. Just running the
code, seeing 9 outputs, and writing `@test length(result) == 9` codifies
whatever the code does — including any bugs — and turns the test into a
regression-snapshot rather than a correctness check.

Concretely, every count and delta needs a comment of the form
"<seed-property-1>; <seed-property-2>; therefore <prediction>". Example for
uni-uni RE→SS: `# 3 singleton groups: S-binding (RE), P-binding (RE),`
`# iso (SS). RE→SS fires on each all-RE group atomically. 2 RE groups → 2`
`# variants. Δ = +1 per variant (1 RE param K becomes (kf, kr) = 2 SS`
`# params, net +1 over the kinetic-group count).` For equivalence-style
tests, the listed expected mechanism literals are the derivation — but each
seed-and-move pair still needs a one-line comment explaining why the list
contains exactly those mechanisms and no others (e.g., "S-group can flip
to SS, P-group can flip to SS, iso group is already SS so excluded").

What this rule rules out:

- "code produces 9 outputs, so I'll write `@test length(result) == 9`".
- Listing expected mechanism literals by copy-pasting compile output.
- Comments that describe what the code does ("returns 2 variants") rather
  than why those variants are the right answer ("RE→SS fires per all-RE
  group; uni-uni has 2 such groups").

```julia
@testset "_expand_<move>" begin

    @testset "<move>: <seed description>" begin
        # SEED — literal @enzyme_mechanism or @allosteric_mechanism
        m_seed = @enzyme_mechanism begin … end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, rxn)
        # Helper validates seed/rxn consistency internally; no per-call-site
        # round-trip assertion needed.

        # MOVE
        result = EnzymeRates._expand_<move>(spec, …)

        # 1. count — derivation REQUIRED, independent of code output
        # <seed-property-1>; <seed-property-2>; therefore N variants because <reason>.
        @test length(result) == N

        # 2. Δ params — derivation REQUIRED, independent of code output
        # <reason this move adds EXPECTED_DELTA params: e.g., "RE→SS adds
        # 1 net param per kinetic group: 1 K → (kf, kr)">.
        for r in result
            @test r.n_fit_params_estimate ==
                spec.n_fit_params_estimate + EXPECTED_DELTA
        end

        # 3. compilability — REQUIRED only when item 4 is property-style.
        # When item 4 is equivalence-style, the
        # `Set(compile_mechanism(r) for r in result) == expected` call
        # already invokes compile_mechanism per `r`, so a separate item-3
        # loop just doubles the @generated compile cost.
        if !equivalence_style
            for r in result
                @test compile_mechanism(r) isa
                    Union{EnzymeMechanism, AllostericEnzymeMechanism}
            end
        end

        # 4. structural change
        if N <= 6 || expected_listed_inline
            # EQUIVALENCE-STYLE — expected mechanisms via DSL macros.
            # Derivation comment REQUIRED: explain why this list is
            # complete and exclusive given the seed and move spec.
            # E.g., "S-group can flip to SS (variant 1), P-group can flip
            # to SS (variant 2), iso group is already SS so excluded".
            expected = Set([
                @enzyme_mechanism begin … end,    # variant 1: <derivation>
                @enzyme_mechanism begin … end,    # variant 2: <derivation>
                …
            ])
            @test Set(compile_mechanism(r) for r in result) == expected
        else
            # PROPERTY-STYLE — assert what changed.
            # Derivation comment REQUIRED on the property itself.
            for r in result
                @test <property of the move>
            end
        end

        # 5. preservation
        for r in result
            @test r.reaction === spec.reaction
            # for allosteric: tags on non-target groups equal; reg sites equal
            # for plain: kinetic_group of untouched steps equal
        end
    end

    # 6. negative cases — each its own @testset
    @testset "<move>: <negative seed> → empty" begin
        m_no_op = @enzyme_mechanism begin … end
        spec = mechanism_spec_from_mechanism_and_rxn(m_no_op, rxn)
        @test isempty(EnzymeRates._expand_<move>(spec, …))
    end

    # 7. cross-type — polymorphic moves only
    @testset "<move>: AllostericMechanismSpec — <case>" begin
        m_seed = @allosteric_mechanism begin … end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, rxn)
        # full checklist again, plus tag-inheritance assertions
    end
end
```

**Cross-type tag-inheritance assertions:**

- `_expand_re_to_ss` on allosteric: converted group keeps its tag; Δ = +1 for
  cheap tags, +2 for `:NonequalRT`.
- `_expand_split_kinetic_group` on allosteric: new group inherits parent's
  tag.
- `_expand_add_dead_end_regulator` on allosteric: new binding-step group is
  auto-tagged `:EqualRT`; allosteric ligands are excluded from
  `eligible_regs`.

**Seed battery per move** (refined when each commit is written):

| Move | Plain seeds | Allosteric seeds | Overlap seeds | Negative seeds |
|---|---|---|---|---|
| `_expand_re_to_ss` | uni-uni, bi-bi (random), ping-pong | `:EqualRT`/`:OnlyR`/`:NonequalRT` group flavors | substrate-as-dead-end-I | all-SS |
| `_expand_split_kinetic_group` | bi-bi size-2 groups, bi-bi after RE→SS | bi-bi allo `:EqualRT`/`:NonequalRT` parent | — | all-singleton groups |
| `_expand_add_dead_end_regulator` | uni-uni+I, bi-bi+I, ping-pong+I, sequential bi-bi+I, two-inh chain | uni-uni allo+I, mixed allo+dead-end | substrate-as-I, product-as-I | no regs, allo-only reg |
| `_expand_to_allosteric` | uni-uni, bi-bi, ping-pong | n/a | — | already-allosteric |
| `_expand_add_allosteric_regulator` | n/a | one-reg, two-reg-same-site, two-reg-different-sites, existing `:EqualRT` site | substrate-as-allo-reg, product-as-allo-reg | non-allosteric, all regs added |
| `_expand_change_allo_state` | n/a | one tagged, multiple tagged, fully relaxed | substrate-as-allo-reg (its ligand tag removable) | non-allosteric, fully `:NonequalRT` |

**Stoichiometry-2 (out of scope):** the `EnzymeReaction` constructor (`src/types.jl:46-49`) errors on duplicate substrate or product names, and the macros pass through to that constructor. Reactions like `2 ADP ↔ ATP + AMP` cannot be expressed in the current API. Adding stoichiometric-coefficient support is a separate scope item; this rewrite tests only the reaction shapes the existing API supports.

**Naming convention:**

```julia
@testset "_expand_re_to_ss" begin
    @testset "MechanismSpec — uni-uni: 2 RE binding groups → 2 variants" begin … end
    @testset "MechanismSpec — all SS → empty" begin … end
    @testset "AllostericMechanismSpec — :EqualRT group: Δ=+1" begin … end
    @testset "AllostericMechanismSpec — :NonequalRT group: Δ=+2" begin … end
end
```

## 4. File structure (pipeline order)

```
test_mechanism_enumeration.jl
│
├── ─── 0. Test infrastructure ──────────────────────────────────
│   ├── mechanism_spec_from_mechanism_and_rxn(m, rxn)         (existing)
│   ├── allosteric_spec_from_mechanism_and_rxn(m, rxn)        (NEW)
│   └── shared @enzyme_reaction defs
│
├── ─── 1. Support functions (no spec input) ───────────────────
│   ├── _catalytic_topologies
│   ├── _competition_patterns
│   ├── _inhibitor_competition_patterns
│   ├── _forms_with_binding_step
│   ├── _substrate_product_dead_end_opportunities
│   └── _expand_substrate_product_dead_ends
│
├── ─── 2. Initialization ──────────────────────────────────────
│   ├── compile_mechanism / EnzymeMechanism constructors
│   └── init_mechanisms
│
├── ─── 3. Base-spec expansion moves (polymorphic) ─────────────
│   ├── _expand_re_to_ss
│   ├── _expand_split_kinetic_group
│   └── _expand_add_dead_end_regulator
│
├── ─── 4. Allosteric expansion moves ──────────────────────────
│   ├── _expand_to_allosteric
│   ├── _expand_add_allosteric_regulator
│   └── _expand_change_allo_state
│
├── ─── 5. Composition ─────────────────────────────────────────
│   ├── dedup!
│   └── expand_mechanisms
│
└── ─── 6. Integration ─────────────────────────────────────────
    └── enumerate_all
```

Sections 1 → 6 monotonically reflect pipeline depth: support fns are called
by init; init is consumed by base moves; allosteric moves run after base
moves; composition wraps moves; integration covers end-to-end.

**Rewrite intensity per section:**

- Sections **3, 4** — full rewrite (literal seeds + 7-item checklist + seed
  battery).
- Section **2** — moderate rewrite (literal seeds where they replace
  `init_mechanisms |> first`).
- Sections **1, 5, 6** — light touch (reorganize into pipeline order;
  assertions mostly preserved).

**What's absorbed / moved out:**

| Existing testset (line) | Disposition |
|---|---|
| "AllostericEnzymeMechanism TR equivalence" (217) | Move to `test_types.jl` — accessor test, not enumeration |
| "test reaction atom balance" (198) | Move to `test_dsl.jl` — `@enzyme_reaction` macro test |
| "Tagged groups exclude T-state params" (2120) | Move to `test_rate_eq_derivation.jl` — `parameters(m)`/canonicalizer test |
| "Metabolite overlap: substrate as dead-end inhibitor" (2203) | Absorbed into per-move overlap seeds |
| "Metabolite overlap: substrate as allosteric regulator" (2233) | Absorbed into per-move overlap seeds |
| "Base-level moves on allosteric specs" (2279) | Absorbed into cross-type sub-testsets within polymorphic-move testsets |
| "C6 iso size limit blocks 4x4" (2408) | Folded into `_catalytic_topologies` testset |
| "init_mechanisms drops unbound regulators" (2428) | Folded into `init_mechanisms` testset |

The "move out of file" rows are flagged for confirmation per-row before
moving; they're not core to the brittleness fix.

## 5. Execution plan

**Pre-work (one commit):**

- Add `allosteric_spec_from_mechanism_and_rxn(m, rxn)` helper at the top of
  the file with a round-trip-validation testset:
  - 3-4 cases via `@allosteric_mechanism`: K-type uni-uni, K-type bi-bi, with
    two reg sites, with `:NonequalRT` regulator. Each asserts
    `AllostericEnzymeMechanism(allosteric_spec_from_mechanism_and_rxn(m, rxn)) === m`.

**Sequenced rewrite (one commit per section):**

1. **Section 1 — support functions**: pure reorganization. No assertion
   changes. Light commit.

2. **Section 2 — init/compile**: rewrite `init_mechanisms` testset to use
   literal seeds; add a dedicated `compile_mechanism` round-trip testset.

3. **Section 3a — `_expand_re_to_ss`**: full rewrite. Acts as the
   **template-validation commit** — if the test pattern has problems, find
   them here on the simplest move.

4. **Section 3b — `_expand_split_kinetic_group`**: same template applied.

5. **Section 3c — `_expand_add_dead_end_regulator`**: same template; absorbs
   "Regulator dummy naming stability" testset and the substrate-as-dead-end
   overlap testset.

6. **Section 4a — `_expand_to_allosteric`**: rewrite with literal seeds.

7. **Section 4b — `_expand_add_allosteric_regulator`**: rewrite; absorbs
   substrate-as-allosteric-regulator overlap testset.

8. **Section 4c — `_expand_change_allo_state`**: rewrite.

9. **Section 5 — composition**: replace `init_mechanisms |> first` patterns
   in `expand_mechanisms` tests with literal seeds; preserve `dedup!` tests
   as-is.

10. **Section 6 — integration**: keep as-is or with minimal cleanup.

11. **Final cleanup commit (optional, per-row confirmation):** move three
    testsets out of this file:
    - `AllostericEnzymeMechanism TR equivalence` → `test_types.jl`
    - `test reaction atom balance` → `test_dsl.jl`
    - `Tagged groups exclude T-state params` → `test_rate_eq_derivation.jl`

**Per-commit verification protocol:**

- Run `julia --project -e 'using Pkg; Pkg.test()'` to confirm green before
  commit.
- Show diff of which existing testsets the new commit absorbs/replaces, so
  coverage regressions are visible at review.
- Each commit message names the move being rewritten and the pre-rewrite
  testsets it absorbs.

**Risk register:**

- **Compilation cost**: ~30–50 literal `@(allosteric_)mechanism` seeds across
  the file. Each triggers `@generated` derivation. Likely tolerable but worth
  measuring after section 3a — if compile time blows up, fall back to
  property-style assertions for some seeds.
- **`==` on compiled mechanisms**: equivalence-style tests rely on
  `Set(compile_mechanism.(results)) == Set(expected)`. Compiled mechanisms
  are singleton types so `===` works directly; `==` falls back to `===` by
  default. Already implicitly used by existing round-trip assertions.
- **Coverage drop**: per-commit diff review is the safety net.

## 6. Bug-handling protocol

When a new test fails, the test does NOT get modified to match buggy
behavior. The test is paused, the bug is surfaced, the production code in
`src/` gets fixed, and the test runs green.

**Protocol on failure:**

1. **Stop the rewrite of the current move.** Don't keep adding tests on top
   of an unaddressed failure.

2. **Diagnose via the systematic-debugging discipline in CLAUDE.md** (Phase
   1: read the error; Phase 2: compare working vs broken; Phase 3: form one
   hypothesis, test minimally).

3. **Surface the failure with a structured report:**
   - The test that fails (file:line + assertion)
   - The seed mechanism (literal `@(allosteric_)mechanism` block)
   - Expected output vs actual output (counts, Δ params, equivalence-set
     diff)
   - Hypothesis: is the expectation wrong, or the code wrong? Why?
   - If unambiguous: state which side is wrong and propose the fix.
   - If unclear: stop and ask.

4. **Decide who's wrong, by category:**
   - **Test expectation wrong** → fix the test. Document why I got it wrong
     in the commit message.
   - **Code wrong** → fix `src/`, NOT the test.
   - **Routine/clear fix** → just fix it, separate commit.
   - **Architectural / unclear** → stop and discuss with Denis before
     changing code.

5. **What's forbidden, regardless of pressure to keep moving:**
   - `@test_broken` to make a known-failing assertion stop blocking.
   - `@test_skip` to disable a test.
   - Weakening an assertion to dodge a count mismatch.
   - Removing the assertion entirely and not replacing it.
   - Commenting the test out.
   - Changing the seed to one that happens to avoid the bug.
   - Any "TODO: revisit later" without a tracked issue.

   If any of these is genuinely the right answer, that's an architectural
   decision per rule 4 — pause and discuss, don't unilaterally do it.

**Commit ordering when a bug is found:**

```
A. fix: <one-line bug summary> (uncovered by <move> test rewrite)
   - src/ change
   - existing tests still pass
   - new test that exercises the bug NOT yet committed
B. test: rewrite <move> tests with literal seeds + checklist
   - test/ changes only
   - now green because of A
```

Both commits land green at HEAD. Reverting B leaves A as a standalone bug
fix. Reverting A only would re-fail B's new tests.

**Bug ledger:** as bugs are uncovered during the rewrite, maintain a running
list in conversation so we can scan for patterns across moves — e.g., if
three moves all forget to update `group_tags`, that's a systemic issue worth
a single broader fix.

**The independent-derivation rule (§3) is what makes bugs visible.** If
the count `N` is reasoned about from the seed and move spec independently
of the code output, a discrepancy between expected and actual surfaces a
bug. If the count is just "what the code returned when I ran it", the
test silently codifies whatever behavior — including bugs — exists today.
This rule is the precondition for the bug-handling protocol to work; the
two are inseparable.

**Most likely bug surfaces:**

- **Equivalence-set mismatches** — code produces extras / drops a variant.
  The new equivalence-style assertions are far stronger than the existing
  count-only checks.
- **Metabolite-overlap seeds in the polymorphic moves** — existing overlap
  tests only cover `_expand_add_dead_end_regulator` and
  `_expand_add_allosteric_regulator`. Splitting / RE→SS /
  `_expand_change_allo_state` aren't currently exercised with overlapping
  names.
- **Cross-type tag-inheritance in moves 4–6** — split tag inheritance is
  partial today; RE→SS tag preservation isn't asserted at all.
