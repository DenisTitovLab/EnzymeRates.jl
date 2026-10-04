# Merged Seeds and Fused Steps Implementation Plan (sub-project C)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Simplify B, classify steps by stoichiometry, let every move handle fused and
Theorell–Chance steps, stop flipping qualifying-chain flanks, add merged and Theorell–Chance
seeds, and expand degenerate seeds without fitting them.

**Architecture:** All enumeration logic lives in `src/mechanism_enumeration.jl`; the step model
in `src/types.jl`; the `:OnlyA` check in `src/thermodynamic_constr_for_rate_eq_derivation.jl`;
the beam in `src/identify_rate_equation.jl`. New seeds are built in `init_mechanisms` from
today's seeds by structural operations (MERGE, ELIM) and a minimal-set search that reuses
`_minimal_gaining_sets`; every test is structural, before any `Mechanism` is built.

**Tech Stack:** Julia 1.12, the package's own `@enzyme_mechanism` / `@allosteric_mechanism` /
`@enzyme_reaction` DSL, `Test`, `TestEnv`.

**Spec:** `docs/superpowers/specs/2026-10-03-merged-seeds-and-fused-steps-design.md` (read it
first; section numbers below refer to it). The read-only review of B that Task 1–4 implement is
`docs/superpowers/specs/2026-10-03-b-simplification-review.md` (line numbers at b340822).

## Global Constraints

- Branch `step-explicit-stoichiometry`. No version bump (stays 0.8.0). Never push.
- One `julia` process at a time (7.7 GB RAM, no swap). The VS Code language server is fine.
  Check with `pgrep -fa 'bin/julia' | grep -v language`. Give long runs `--heap-size-hint=2500M`.
- Full suite (about 15 minutes; it must pass before every commit):
  `julia --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])' > /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/f13ff7d8-8f39-476b-b9d6-fead260f40c0/scratchpad/suite-<task>.log 2>&1`
  started in the background, then polled in bounded foreground loops
  (`for i in $(seq 1 28); do grep -q "Testing EnzymeRates tests passed\|Test Summary\|ERROR" LOG && break; sleep 20; done`).
  A subagent's background job dies when it hands back: poll to the end yourself.
- Focused run of the enumeration file (about 4 minutes; it skips Aqua and JET):
  `julia --project --heap-size-hint=2500M -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_mechanism_enumeration.jl")'`.
  A helper from another test file (`random_reduced_params`) is undefined there; an
  `UndefVarError` for it is an artifact. The fastest red/green loop is one testset copied into a
  scratch file under the scratchpad with that prelude.
- CLAUDE.md rules bind every task: TDD; 92-character lines, 4-space indent; every new file starts
  with two `# ABOUTME:` lines; no comment narrates former behaviour or says "new", "now",
  "improved", "refactored"; never delete a test because it fails (rewrite it to the rule);
  `git add` named files only, never `-A`; commit messages end with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` and
  `Claude-Session: https://claude.ai/code/session_015DG9FSn9ZSbZhzeSDkqqnS`.
- Enumeration-test rules (CLAUDE.md): every fixture written inline in its testset with the
  macros; a move test asserts `length(children) == n` and `Set(children) == Set(expected)` with
  every expected child written out; chemistry written the enumerator's way; file-level helpers
  prefixed `_testhelper_`. Every expected child is hand-derived in a comment and cross-checked
  with `_testhelper_identifiable_rank` (`test/test_mechanism_enumeration.jl:92`). If a
  derivation disagrees with this plan's expectation, STOP and report; do not "fix" the test to
  match the code.
- `first(init_mechanisms(rxn))` must stay today's first seed for every reaction: the new seeds
  are appended after today's seeds, in today's order (many tests rely on it).
- The `rate_equation` performance test (`test/test_rate_eq_derivation.jl`,
  `test_rate_equation_performance`) must stay green; nothing here should touch it.
- Restricting the flip rule to `Mechanism` and seeds to `init_mechanisms` is deliberate (spec
  Non-goals); do not extend either to `AllostericMechanism`.

## Review Focus

- Ter-ter seeding: `init_mechanisms(ter_ter_rxn)` must finish in seconds, not minutes, with the
  variants added (Task 10 pins `< 60 s`).
- Required regulators: `seed_mechanisms` on bi-bi with a required competitive inhibitor must
  still return seeds, merged ones among them (Task 10).
- Unequal substrate and product counts: uni-bi gains exactly one merged seed per undecorated
  ordered topology (n + m − 2 = 1) (Task 10).
- Multi-subunit allosteric children of a merged seed: `_expand_to_allosteric` at
  `oligomeric_state: 2` constructs every child without error (Task 6).
- Result rows of children of expand-only seeds carry no parent, and `_rows_to_dataframe`
  accepts them (Task 11).

---

## Wave 0: simplify B (spec section 0)

### Task 1: One home for each of B's arguments; stale text; `ss_mirror` rank pin

**Files:**
- Modify: `src/mechanism_enumeration.jl` (docstrings only, the line ranges below are at b340822)
- Modify: `test/test_mechanism_enumeration.jl:7417` (one added assertion), and the phrase test
  near 10477 if the `seed_mechanisms` message changes
- Read first: `docs/superpowers/specs/2026-10-03-b-simplification-review.md` items 2 and 14,
  and handoff `docs/superpowers/specs/2026-10-02-c-handoff.md` item 4

**Interfaces:** none change. Output of every function is identical.

- [ ] **Step 1: Add the `ss_mirror` rank pin (handoff cleanup 2).** After
  `@test isempty(ER._redundant_copy_groups(ss_mirror))` (line 7417) add:

```julia
    @test fitted(ss_mirror) == 8 && _testhelper_identifiable_rank(ss_mirror) == 6
```

  (`fitted` is defined at the top of that testset.) Run that testset in a scratch file; expect
  PASS (8/6 is in B's spec table).

- [ ] **Step 2: Deduplicate the docstrings.** Keep each argument in one home and cut it elsewhere,
  exactly per the review's item 2 table:
  - "a complex that duplicates only a copy's complex is no twin": keep in `_productive_twin`;
    delete from `_assert_emission_rules`, the split docstring, `_redundant_copy_groups`, the
    dead-end docstring;
  - "all-twin with a failing gauge is kept; the gauge is the proof": keep in
    `_redundant_copy_groups`;
  - "a placement never makes an older copy redundant": keep in `_dead_end_child`;
  - Theorem 2's finite merge K_h + K\*: keep in `_gauge_rescaling`;
  - the twin definition: keep in `_productive_twin`;
  - "a zero-flux group's two constants enter only as their ratio": keep in
    `_flux_carrying_groups`;
  - per-site competition: keep in `_expand_add_dead_end_regulator`; reduce
    `_forms_with_binding_step_native`'s docstring to its contract.
  Delete the evidence clause "(four such copies at depth 2 of R6 have full rank)", the label
  "(Case 3's is)", and the narration of the tolerated two-copies class in the split docstring.
  Shrink `_assert_emission_rules`'s docstring to: the two rules, why a parent must obey them
  (the flip tests only its flipped groups), and which moves filter their children.

- [ ] **Step 3: Fix the stale and inaccurate text.**
  - `expand_mechanisms` docstring: replace "every competitive-inhibitor group creates a new
    complex" with "no competitive-inhibitor group is redundant (`_redundant_copy_groups`)".
  - `_gauge_rescaling` (handoff cleanup 3): the ρ of a shared class is one symbol; say
    "twins formed by RE bindings of the copied ligand in one kinetic group h share one factor
    (Theorem 2's finite merge gives it as ρ = K_h/(K_h + K\*))".
  - The dead-end docstring (handoff cleanup 4): a `Mechanism` candidate takes the
    `_gauge_rescaling` path in `_dead_end_child`; an allosteric child runs
    `_redundant_copy_groups`.

- [ ] **Step 4: Shorten the `seed_mechanisms` error and docstring.** The error keeps the
  regulator names, one sentence of why ("a substrate or product declared as a competitive
  inhibitor binds only where its complex is not redundant; a uni-uni seed has no such site"),
  and the two keywords. Remove the docstring sentence that repeats the message. Update the
  phrase test near `test/test_mechanism_enumeration.jl:10477` to the new wording (keep it
  asserting the regulator name and both keywords).

- [ ] **Step 5: Full suite, then commit.**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Give each of the copy and flux rules' arguments one home"  # plus the two trailer lines
```

### Task 2: Pure refactors of B's code

**Files:**
- Modify: `src/mechanism_enumeration.jl`
- Read first: the review's items 3, 6, 7, 8, 9, 10, 12 and handoff cleanup 1

**Interfaces:**
- Produces: `_with_equilibrium(s::Step, flag::Bool)::Step` in `src/mechanism_enumeration.jl`
  (Tasks 9 and 10 use it); one method
  `_expand_add_dead_end_regulator(m::Union{Mechanism, AllostericMechanism}, rxn; exclude_regs)`
  replacing the two wrappers and `_expand_add_dead_end_regulator_native` (Task 7 edits it).
- Output of every function is identical; no test changes except call sites of renamed internals.

- [ ] **Step 1: `_with_equilibrium`.** Add near `_flip_group_to_ss`:

```julia
"""`s` with its rapid-equilibrium flag set to `flag`."""
_with_equilibrium(s::Step, flag::Bool) =
    Step(from_species(s), to_species(s), consumed(s), released(s), flag)
```

  Use it at the four rebuild sites (init's SS-tagging loop, `_flip_group_to_ss`,
  `_all_steady_state`, `_revert_zero_flux_parts`). `_flip_group_to_ss` becomes
  `[gi == g ? _with_equilibrium.(gr, false) : gr for (gi, gr) in enumerate(groups)]`;
  `_all_steady_state(groups) = [_with_equilibrium.(group, false) for group in groups]`.

- [ ] **Step 2: One dead-end method.** Merge the two `_expand_add_dead_end_regulator` methods and
  `_expand_add_dead_end_regulator_native` into one method on
  `Union{Mechanism, AllostericMechanism}` that computes
  `additional_excluded = m isa AllostericMechanism ? Set(name(l) for site in regulatory_sites(m) for l in ligands(site)) : Set{Symbol}()`.
  Keep the docstring of the public-facing name; grep `src/` and `test/` for `_native` callers and
  update them.

- [ ] **Step 3: One RE-segment preamble helper.** Add in the enumeration file:

```julia
"""`_re_segment_extras(groups)` with each form's index and segment: `(species, segments,
extras, index, segment_of)`, where `index[sp]` is the form's position in `species` and
`segment_of[i]` the segment holding form `i`."""
function _indexed_re_segments(groups::Vector{Vector{Step}})
    species, segments, extras = _re_segment_extras(groups)
    index = Dict(sp => i for (i, sp) in enumerate(species))
    segment_of = zeros(Int, length(species))
    for (k, members) in enumerate(segments), i in members
        segment_of[i] = k
    end
    species, segments, extras, index, segment_of
end
```

  and use it in `_flux_carrying_steps`, `_productive_twin` and `_hyperbolic_catalysis` in place
  of their six-line preambles. Do not change `_re_segment_extras`.

- [ ] **Step 4: Inline `_copy_complex` and `_composition`.** `_copy_complex` has one caller (the
  dead-end kernel); `_composition` is used only in `_productive_twin` (make it a local closure
  there). Rewrite the test line that calls `_composition` (near
  `test/test_mechanism_enumeration.jl:7075`) to assert the same fact through `_productive_twin`.

- [ ] **Step 5: Tidy `_split_gain_test`.** One loop over both parts of each selected unit (set
  the kind of every step; give part 2's steps the new id). Optionally make
  `_revert_zero_flux_parts` return `(parts, reverted::Bool)` instead of relying on
  `parts !== bp` identity, and update its one caller.

- [ ] **Step 6: `_redundant_copy_groups(am)`.** Remove the dead clauses
  `tags[g] !== :OnlyA &&` in `binds_inactive` and `h != g &&` in the "binds nothing" Dict
  (review item 6 proves both dead), and build the inactive view with `_reachable_from_free`:
  keep a step of a non-`:OnlyA` group iff both its forms' names are in
  `_reachable_from_free(non_onlya_groups)` (review item 12). This drops the
  `_state_mechanism(am, :I)` round trip and the `_forward_sides` string keys.

- [ ] **Step 7: Build `_productive_twin(steps(m))` once per dead-end parent (handoff cleanup
  1).** The kernel computes it and `_copy_twin_test(m)`; for a `Mechanism` the latter is the
  former tested against `nothing`. Compute the twin function once and derive the predicate from
  it (Task 3 deletes `_copy_twin_test` entirely; here only stop computing the same thing twice).

- [ ] **Step 8: Full suite, then commit** (one commit per step is fine; at least one per
  logical item). Message example: "Build one RE-segment index in each B predicate".

### Task 3: One "every complex has a twin" predicate

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_copy_twin_test`, `_twin_only`, `_all_twin`,
  `_gauge_rescaling`, the split, `_redundant_copy_groups`, `_dead_end_child`, the dead-end
  kernel)
- Modify: `test/test_mechanism_enumeration.jl:7180-7221, 7275, 7318`
- Modify: `docs/src/developer.md:104`
- Read first: the review's item 1 (it argues output identity branch by branch)

**Interfaces:**
- Produces: `_all_twin(part::Vector{Step}, twin)::Bool` (part-level; `twin` is a
  `_productive_twin(groups)` closure); `_gauge_rescaling(groups, g, twin, label)` returns
  `nothing` unless `_all_twin(groups[g], twin)`.
- Removes: `_copy_twin_test`, `_twin_only`.

- [ ] **Step 1: Rewrite the tests that call the deleted helpers, first.** In the testset that
  holds `onlya`, `nonequal`, `all_onlya` (lines 7180–7221), replace each
  `_copy_twin_test(am)(form(am, :E), Sinh, tag)` assertion by the move-level fact it encodes:
  build the child with the copy placed and assert `_redundant_copy_groups(child)`:

```julia
    # An :EqualAI copy of S at E: new in the inactive state when S binds :OnlyA there,
    # so not redundant; redundant when S binding is :NonequalAI (a twin in both states).
    place(am, tag) = ER.AllostericMechanism(ER.reaction(am),
        vcat(ER.steps(am), [[ER.Step(form(am, :E),
            ER.Species(ER.Metabolite[Sinh], :E), ER.Metabolite[Sinh], ER.Metabolite[], true)]]),
        vcat(ER.cat_allo_states(am), [tag]), ER.catalytic_multiplicity(am),
        ER.RegulatorySite[])
    @test isempty(ER._redundant_copy_groups(place(onlya, :EqualAI)))
    @test !isempty(ER._redundant_copy_groups(place(onlya, :OnlyA)))
    @test !isempty(ER._redundant_copy_groups(place(nonequal, :EqualAI)))
    @test isempty(ER._redundant_copy_groups(place(all_onlya, :EqualAI)))
    @test !isempty(ER._redundant_copy_groups(place(all_onlya, :OnlyA)))
```

  (`place` must use the reaction extended with the inhibitor if the constructor requires it:
  use `ER._add_competitive_inhibitor(ER.reaction(am), :S)` when `Sinh` is `S`'s copy.) Each
  verdict must equal the deleted assertion's (true twin ⇔ redundant here, since each of these
  copy groups is a single step whose gauge is consistent; derive that in a comment). Add one
  rank pin for the `:OnlyA`-copy case:
  `@test fitted(c) == _testhelper_identifiable_rank(c) + 1` for `c = place(onlya, :OnlyA)`.
  Replace lines 7275 and 7318 with
  `@test ER._all_twin(ER.steps(case3)[g3], ER._productive_twin(ER.steps(case3)))` (and the same
  for `h1`, `gh`).

- [ ] **Step 2: Run the rewritten testsets against the old code**; expect PASS (they assert
  behaviour, not helpers) except the two `_all_twin` lines (the part-level signature does not
  exist yet): expect `MethodError`.

- [ ] **Step 3: Implement.** Replace `_all_twin(groups, g, twin)` by

```julia
"""Whether `part`, steps of one competitive-inhibitor group, binds the inhibitor only where
its complex has a productive twin (`twin`, a `_productive_twin` of one state's graph)."""
function _all_twin(part::Vector{Step}, twin)
    ligand = bound_metabolite(first(part))
    ligand isa CompetitiveInhibitor &&
        all(s -> twin(from_species(s), ligand) !== nothing, part)
end
```

  make `_gauge_rescaling` begin with `_all_twin(groups[g], twin) || return nothing`, collapse
  every `_all_twin(...) && _gauge_rescaling(...) !== nothing` pair to the second half, change
  the split's `redundant_part` to `any(part -> _all_twin(part, twin), bp) && begin … end` with
  `twin = _productive_twin(groups)`, drop `_dead_end_child`'s `all_twin` argument (the
  `Mechanism` method: `_gauge_rescaling(groups, length(groups), twin, identity) !== nothing &&
  return nothing`; the allosteric method: `_all_twin(groups[end], twin) &&
  !isempty(_redundant_copy_groups(child))`), compute `_productive_twin(steps(m))` once in the
  kernel, and delete `_copy_twin_test` and `_twin_only`.

- [ ] **Step 4: Focused run of the enumeration file; then `docs/src/developer.md:104`** stops
  naming `_copy_twin_test` (name `_redundant_copy_groups` instead).

- [ ] **Step 5: Full suite; confirm the population pins at
  `test/test_mechanism_enumeration.jl:9700` and `:9716` are unchanged; commit** ("Test every
  copy complex for a twin with one predicate").

### Task 4: Delete the split's unit pre-filter; rerun B's regression on the opening wave

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`redundant_part` and its use in
  `_expand_split_kinetic_group`; the split docstring's sentence about it)
- Modify: `docs/superpowers/specs/2026-09-30-exact-filters-design.md` §4 (the "Copy groups"
  bullet)
- Create: `docs/superpowers/specs/2026-10-03-opening-wave-regression.md`
- Read first: the review's item 4; B's record `2026-09-30-exact-filters-regression.md` (its
  method)

**Interfaces:** none. The split's children can only grow, each satisfying both rules.

- [ ] **Step 1: Delete `redundant_part`** and the `redundant_part(g, bp) && continue` line; keep
  the final `filter!(c -> isempty(_redundant_copy_groups(c)), …)`. Remove the split docstring's
  sentence that a redundant-part bipartition "is not a unit"; say instead that children with a
  redundant copy group are not emitted.

- [ ] **Step 2: Correct B's spec §4.** Replace "Every superset of such a unit keeps the part, and
  splitting other groups only relaxes the part's gauge. Excluding the unit loses only children
  whose part is provably a phantom." with: "Splitting a group that forms the part's twins can
  break the gauge (H1 with its two A groups merged into one is redundant; splitting that A group
  makes it identifiable), so such a bipartition stays a unit, and the child filter alone removes
  children with a redundant copy group."

- [ ] **Step 3: Regression on the opening wave.** Write a scratch script (in the scratchpad) that
  enumerates, on a detached worktree at b340822 (`git worktree add --detach ../er-b340822
  b340822`) and on the current tip, `init_mechanisms` plus two levels of the three catalytic
  moves on R1–R6 (the six reactions of the findings document, §"How the findings were
  verified"; R6 level 2 in full) and the allosteric sample of A's record, keyed by steps with
  explicit lists (`_forward_sides` of every step, per group, sorted). Assert: the two populations
  are identical except for mechanisms only on the tip; every tip-only mechanism is a split child
  satisfying both emission rules; report their count and, for each, fitted and rank
  (`_testhelper_identifiable_rank` logic copied into the script). Run the two sides one after
  the other (one Julia process at a time).

- [ ] **Step 4: Write the record** `docs/superpowers/specs/2026-10-03-opening-wave-regression.md`
  (two ABOUTME-free markdown headers are fine; it is not code): method, per-reaction per-level
  counts on both sides, the tip-only mechanisms with fitted/rank, and timings of R6 levels 0–2
  on both sides. If any tip-only mechanism has fitted > rank, list it; do not change code.

- [ ] **Step 5: Full suite; re-pin a population count only if the record shows a change, with a
  comment that states the count and what makes it (no narration of the pre-filter); commit**
  ("Let the split's child filter alone enforce the copy rule"), then remove the worktree
  (`git worktree remove ../er-b340822`).

---

## Wave 1: steps classified by stoichiometry (spec section 1)

### Task 5: A binding is a step that takes up exactly one metabolite and gives off none

**Files:**
- Modify: `src/types.jl:162-217` (the `Step` docstring and constructor, `bound_metabolite`,
  `is_binding`), `src/types.jl:643-673` (`_step_kind`, `_assert_uniform_groups` docstring)
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl:323` (`_count_kind` docstring if
  it says "pure binding")
- Test: `test/test_types.jl`, `test/test_rate_eq_derivation.jl` (fused fixtures near 1186–1360)

**Interfaces:**
- Produces: `bound_metabolite(s)` = `only(consumed(s))` when `length(consumed(s)) == 1 &&
  isempty(released(s))`, else `nothing`; `is_binding(s) = bound_metabolite(s) !== nothing`;
  `_is_chemistry(s::Step)::Bool` in `src/types.jl`; `_step_kind(s) = (consumed(s),
  released(s))`.
- Every derivation `is_binding` test then admits fused bindings: `_re_weight_ratio`
  (`src/rate_eq_derivation.jl:231`), parameter types (`src/rate_eq_derivation.jl:1145`,
  `src/thermodynamic_constr_for_rate_eq_derivation.jl:40`), binding-K sets
  (`src/rate_eq_derivation.jl:121-126`, `:1218`), `_count_kind`. No code change there; verify
  each still reads "takes up exactly one metabolite" correctly.

- [ ] **Step 1: Failing tests** in `test/test_types.jl`:

```julia
@testset "a step that takes up one metabolite and gives off none binds it" begin
    ER = EnzymeRates
    sp(mets...) = ER.Species(ER.Metabolite[mets...], :E)
    A, B, P, Q = ER.Substrate(:A), ER.Substrate(:B), ER.Product(:P), ER.Product(:Q)
    plain = ER.Step(sp(), sp(A), ER.Metabolite[A], ER.Metabolite[], true)
    @test ER.bound_metabolite(plain) == A && !ER._is_chemistry(plain)
    fused = ER.Step(sp(A), sp(P, Q), ER.Metabolite[B], ER.Metabolite[], true)
    @test ER.bound_metabolite(fused) == B && ER._is_chemistry(fused)
    # A step that only gives off one metabolite is stored as the binding it reverses.
    release = ER.Step(sp(A, B), sp(Q), ER.Metabolite[], ER.Metabolite[P], false)
    @test ER.from_species(release) == sp(Q) && ER.to_species(release) == sp(A, B)
    @test ER.consumed(release) == ER.Metabolite[P] && ER._is_chemistry(release)
    tc = ER.Step(sp(A), sp(Q), ER.Metabolite[B], ER.Metabolite[P], false)
    @test ER.bound_metabolite(tc) === nothing && ER._is_chemistry(tc)
    iso = ER.Step(sp(A, B), sp(P, Q), ER.Metabolite[], ER.Metabolite[], false)
    @test ER.bound_metabolite(iso) === nothing && ER._is_chemistry(iso)
end

@testset "a fused binding shares a kinetic group with a plain binding of its metabolite" begin
    # Merged decorated ordered bi-bi, variant {A, Pˣ}: the B group holds the fused
    # E(A) + B → E(P, Q) and the dead-end E(Q) + B ⇌ E(B, Q). Fitted: A (2) + B (1)
    # + P (2) + Q (1) − 1 Haldane = 5.
    em = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))
            E(Q) + P <--> E(P, Q)
            E + Q ⇌ E(Q)
        end
    end
    @test length(EnzymeRates.fitted_params(em)) == 5
end
```

  and in `test/test_rate_eq_derivation.jl` add the same mechanism to `_testhelper_fused_cases`
  (so the mass-action oracle checks the mixed group). Run: expect the first testset to fail at
  `_is_chemistry` (undefined) and the fused `bound_metabolite` (`nothing` today), the second
  to fail in the constructor (`_assert_uniform_groups`: "differ in kind").

- [ ] **Step 2: Implement** in `src/types.jl`:

```julia
function bound_metabolite(s::Step)
    length(s.consumed) == 1 && isempty(s.released) ? only(s.consumed) : nothing
end
is_binding(s::Step) = bound_metabolite(s) !== nothing

"""Whether `s` changes the enzyme beyond the metabolite it binds: an isomerization, a fused
binding (`E(A) + B → E(P, Q)`) or a Theorell–Chance step; every step but a plain binding, whose
`to_species` is its `from_species` with the metabolite added. The allosteric moves and the
`:OnlyA` check read it to tell catalysis from binding."""
_is_chemistry(s::Step) =
    !(is_binding(s) && _binds_ligand(from_species(s), to_species(s), bound_metabolite(s)))
```

  In the constructor replace
  `if isempty(c) && length(r) == 1 && _binds_ligand(to_species, from_species, only(r))` with
  `if isempty(c) && length(r) == 1`. Set `_step_kind(s) = (consumed(s), released(s))` and
  rewrite its docstring and `_assert_uniform_groups`' docstring: a group's steps take up and give
  off the same metabolites. Rewrite the `Step` and `bound_metabolite` docstrings: a binding takes
  up exactly one metabolite and gives off none, whatever happens to the enzyme's composition; a
  plain binding only adds it.

- [ ] **Step 3: Run `test/test_types.jl` and `test/test_rate_eq_derivation.jl` focused**; expect
  the new tests to pass. Fix every test that pins the old semantics by rewriting it to the new
  rule (grep both files and `test/test_dsl.jl` for `bound_metabolite`, `is_binding`, fused-step
  name pins such as `K_E_S_to_EP`, and orientation tests near
  `test/test_rate_eq_derivation.jl:1337-1360`). An RE fused binding's constant is now a
  dissociation constant named in the release direction (`K_EP_to_E_S` for `E + S ⇌ E(P)`): the
  mass-action oracle must still pass on every fused case, and the orientation test must still
  find written-direction independence.

- [ ] **Step 4: Regression of names on non-fused mechanisms.** In a scratch script, render
  `rate_equation_string` for every `MECHANISM_TEST_SPECS` mechanism and for
  `init_mechanisms` of `bi_bi_pp_rxn` plus one level of `_expand_re_to_ss` and
  `_expand_split_kinetic_group`, on a detached worktree at the task's base commit and on the
  working tree; assert the strings are identical for every mechanism without a fused or
  Theorell–Chance step, and list those that differ (only fused fixtures may). Paste the counts
  into the commit message body.

- [ ] **Step 5: Full suite, commit** ("Classify a binding by its stoichiometry").

---

## Wave 2: the moves on fused steps (spec section 2)

### Task 6: Accept fused parents; the allosteric moves read `_is_chemistry`

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_assert_chemistry_is_iso` and its call in
  `expand_mechanisms`; `_assert_mechanism_invariants(::Mechanism)`; `_expand_to_allosteric`;
  `_partial_onlya_catalysis`; `_expand_change_allo_state`)
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl:568-580`
  (`_onlya_haldane_violation`)
- Test: `test/test_mechanism_enumeration.jl` (near 9740–9750, which asserts
  `_assert_chemistry_is_iso`), `test/test_types.jl:1853`, and the test that asserts refusal of a
  fused `:OnlyA` ping-pong (grep "Tag the cycle's chemical step" in `test/`)

**Interfaces:**
- Consumes: `_is_chemistry` (Task 5).
- A group is a chemistry group iff `any(_is_chemistry, group)` (a group may hold a fused binding
  and a plain binding of the same metabolite; the whole group then counts as chemistry).

- [ ] **Step 1: Failing tests.** Replace the two tests that assert `_assert_chemistry_is_iso`
  rejections (`test/test_types.jl:1853`, `test/test_mechanism_enumeration.jl:9746-9749`) with
  this move test (exact children):

```julia
@testset "expand_mechanisms on a merged uni-uni" begin
    # E + S → E(P) fused (SS), E + P ⇌ E(P) (RE). Flip: the P group cuts E from E(P),
    # carries flux, no bottomless segment: one child, both steps SS. Split: single-step
    # groups, nothing. Dead end: no regulator. To-allosteric at multiplicity 1: the
    # chemistry group (the fused step) is :OnlyA, the binding subsets range over {P}:
    # one K-type child; no regulator, so no V-type.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    withrxn(em) = EnzymeRates.Mechanism(rxn, EnzymeRates.steps(EnzymeRates.Mechanism(em)))
    m = withrxn(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(P)
            E + P ⇌ E(P)
        end
    end)
    flipped = withrxn(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(P)
            E + P <--> E(P)
        end
    end)
    k_type = EnzymeRates.AllostericMechanism(rxn, EnzymeRates.steps(m),
        [:OnlyA, :OnlyA], 1, EnzymeRates.RegulatorySite[])
    kids = EnzymeRates.expand_mechanisms([m], rxn)
    @test length(kids) == 2
    @test Set(kids) == Set([flipped, k_type])
end
```

  (Order the `[:OnlyA, :OnlyA]` tags after the constructor's canonical group order; build
  `k_type` with tags in the order of `EnzymeRates.steps(m)`.) Add an `_expand_change_allo_state`
  exact-children test on `k_type`: relaxing P gives one child (P `:NonequalAI`); relaxing the
  chemistry group gives a partial `:OnlyA` catalysis and is dropped; so exactly one child.
  Add a Theorell–Chance K-type test: `E + A ⇌ E(A)`, `E(A) + B <--> E(Q) + P`,
  `E + Q ⇌ E(Q)` on a bi-bi reaction at multiplicity 1: chemistry = the TC group; subsets of
  {A, Q}: three children, each with the TC group `:OnlyA` (the Haldane check drops the
  `:OnlyA` chemistry group, so no cycle row remains and nothing is rejected; derive this).
  Rewrite the refusal test of a fused `:OnlyA` ping-pong into an acceptance test
  (`AllostericMechanism(...)` constructs).

  Review Focus: add
  `@test all(c -> c isa EnzymeRates.AllostericMechanism, EnzymeRates._expand_to_allosteric(m2, rxn2))`
  for the merged ordered variant {A, Pˣ} (`E + A <--> E(A)`, `E(A) + B ⇌ E(P, Q)`,
  `E(Q) + P <--> E(P, Q)`, `E + Q ⇌ E(Q)`) on a bi-bi reaction with `oligomeric_state: 2`,
  and assert the exact count of multiplicity-2 children (derive it: chemistry group :OnlyA,
  binding subsets of {A, P, Q} = 7, filtered by `_hyperbolic_catalysis`).

- [ ] **Step 2: Run; expect failures** (`expand_mechanisms` throws "folds chemistry";
  `_expand_to_allosteric` treats the fused group as a binding).

- [ ] **Step 3: Implement.** Delete `_assert_chemistry_is_iso` and its call. Remove the
  "binding or isomerization" loop and docstring line from `_assert_mechanism_invariants`. In
  `_expand_to_allosteric`:
  `chem = [g for g in 1:n_g if any(_is_chemistry, steps(m)[g])]`,
  `bind = [g for g in 1:n_g if !(g in chem)]` (rename the locals from `iso`). In
  `_partial_onlya_catalysis`:
  `live = any(any(_is_chemistry, cat_steps[g]) && cat_allo_states[g] !== :OnlyA for g in eachindex(cat_steps))`.
  In `_expand_change_allo_state`: the chemistry groups are
  `[g for g in eachindex(cs) if any(_is_chemistry, cs[g])]`, and the binding loop skips them.
  In `_onlya_haldane_violation`:
  `keep = [g for g in eachindex(cat_steps) if !(cat_allo_states[g] === :OnlyA && any(_is_chemistry, cat_steps[g]))]`
  and the `:OnlyA` bindings are the groups with `cat_allo_states[g] === :OnlyA &&
  !any(_is_chemistry, cat_steps[g])`. Update every docstring that says "isomerization" for
  chemistry in these functions (they describe chemistry steps; say "chemistry step").

- [ ] **Step 4: Focused runs** of the enumeration file and `test/test_types.jl`; then the full
  suite. Every allosteric exact-children test must still pass unchanged (today's mechanisms have
  no fused step, so `any(_is_chemistry, group)` equals `is_iso(rep)` for them).

- [ ] **Step 5: Commit** ("Let the moves expand fused and Theorell–Chance parents").

### Task 7: Dead-end sites are the forms where the competing metabolite is free

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_forms_with_binding_step_native`, its two call
  sites in the dead-end method)
- Test: `test/test_mechanism_enumeration.jl` (new testset near the dead-end testsets, ~2669)

**Interfaces:**
- Produces: `_forms_where_free(m, role::Type{<:Metabolite}, met_name::Symbol)::Set{Symbol}`
  replacing `_forms_with_binding_step_native` (same role semantics: `Reactant` for productive
  competition, `CompetitiveInhibitor` for an existing inhibitor).

- [ ] **Step 1: Failing test** (exact children, hand-derived; the expected list below is the
  derivation, verify it):

```julia
@testset "Mechanism — dead-end sites on a Theorell–Chance step" begin
    # Ordered Theorell–Chance, all SS: E + A → E(A), E(A) + B → E(Q) + P, E + Q → E(Q).
    # Sites: A is free at E, B at E(A) (taken up by the TC step), P at E(Q) (given off
    # there), Q at E. A foreign inhibitor I competes with one or both substrates and one
    # or both products; a site is skipped when it already holds a competing reactant.
    # {A}×{P}: E, E(Q). {A}×{Q}: E. {A}×{P,Q}: E (E(Q) holds Q). {B}×{P}: E(A), E(Q).
    # {B}×{Q}: E, E(A). {B}×{P,Q}: E, E(A) (E(Q) holds Q). {A,B}×{P}: E, E(Q) (E(A)
    # holds A). {A,B}×{Q}: E. {A,B}×{P,Q}: E. Distinct placements: {E}, {E, E(Q)},
    # {E(A), E(Q)}, {E, E(A)}. Mirrors join two sites: E–E(Q) by Q's binding, E–E(A) by
    # A's, E(A)–E(Q) by the TC step, each in its parent's group, at its parent's flag.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: I
    end
    withrxn(em) = EnzymeRates.Mechanism(
        EnzymeRates._add_competitive_inhibitor(rxn, :I),
        EnzymeRates.steps(EnzymeRates.Mechanism(em)))
    m = EnzymeRates.Mechanism(rxn, EnzymeRates.steps(EnzymeRates.Mechanism(
        @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(Q) + P
                E + Q <--> E(Q)
            end
        end)))
    at_E = withrxn(@enzyme_mechanism begin
        substrates: A, B; products: P, Q; regulators: I
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(Q) + P
            E + Q <--> E(Q)
            E + I ⇌ E(I)
        end
    end)
    at_E_EQ = withrxn(@enzyme_mechanism begin
        substrates: A, B; products: P, Q; regulators: I
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(Q) + P
            (E + Q <--> E(Q), E(I) + Q <--> E(I, Q))
            (E + I ⇌ E(I), E(Q) + I ⇌ E(I, Q))
        end
    end)
    at_EA_EQ = withrxn(@enzyme_mechanism begin
        substrates: A, B; products: P, Q; regulators: I
        steps: begin
            E + A <--> E(A)
            (E(A) + B <--> E(Q) + P, E(A, I) + B <--> E(I, Q) + P)
            E + Q <--> E(Q)
            (E(A) + I ⇌ E(A, I), E(Q) + I ⇌ E(I, Q))
        end
    end)
    at_E_EA = withrxn(@enzyme_mechanism begin
        substrates: A, B; products: P, Q; regulators: I
        steps: begin
            (E + A <--> E(A), E(I) + A <--> E(A, I))
            E(A) + B <--> E(Q) + P
            E + Q <--> E(Q)
            (E + I ⇌ E(I), E(A) + I ⇌ E(A, I))
        end
    end)
    kids = EnzymeRates._expand_add_dead_end_regulator(m, rxn)
    @test length(kids) == 4
    @test Set(kids) == Set([at_E, at_E_EQ, at_EA_EQ, at_E_EA])
end
```

  (Check the DSL accepts `; regulators: I` on one line as used here; if not, use separate
  lines as elsewhere in the file. A foreign inhibitor has no twins, so the copy rule keeps
  every placement; rank-check each child.) Run: expect 1 child today ({E} only: B and P have no
  binding step).

- [ ] **Step 2: Implement** `_forms_where_free`:

```julia
"""Names of the forms where the metabolite named `met_name`, in the role `role`, is free on
some step: the step's `from_species` when it takes the metabolite up, its `to_species` when it
gives it off. `Reactant` asks where a substrate or product binds productively;
`CompetitiveInhibitor` where an inhibitor, a reactant's copy included, binds its dead-end
site. The dead-end move targets these forms for competition with the metabolite."""
function _forms_where_free(m::Union{Mechanism, AllostericMechanism},
                           role::Type{<:Metabolite}, met_name::Symbol)
    forms = Set{Symbol}()
    for group in steps(m), s in group,
        (form, mets) in ((from_species(s), consumed(s)), (to_species(s), released(s)))
        any(x -> x isa role && name(x) == met_name, mets) && push!(forms, name(form))
    end
    forms
end
```

  and replace both calls. Delete `_forms_with_binding_step_native`.

- [ ] **Step 3: Focused run of every dead-end testset, then the full suite** (today's mechanisms
  have no Theorell–Chance step, so every existing expectation holds). Commit ("Place dead-end
  inhibitors where the competing metabolite is free").

### Task 8: The flip never flips a flank of a qualifying chain with steady-state chemistry

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_expand_re_to_ss` and its docstring; a new
  `_chain_flank_groups`)
- Test: `test/test_mechanism_enumeration.jl` (new testsets near the flip testsets; every
  existing exact-children flip test whose parent has a qualifying chain; the population pin
  at 9700)

**Interfaces:**
- Produces: `_chain_flank_groups(m::Mechanism)::Set{Int}` and
  `_chain_flank_groups(::AllostericMechanism) = Set{Int}()`.

- [ ] **Step 1: Failing tests** (exact children):

```julia
@testset "Mechanism — the flip skips the flanks of a qualifying chain" begin
    # Ordered RE seed: E + A ⇌ E(A), E(A) + B ⇌ E(A, B), E(A, B) → E(P, Q) (SS),
    # E(Q) + P ⇌ E(P, Q), E + Q ⇌ E(Q). E(A, B) and E(P, Q) each have one other step, a
    # binding into them, and every group holds one step: a qualifying chain (RSR). Its
    # flanks, the B and P groups, never flip: a flipped flank gives its parent's family
    # plus a phantom (chain lemma). The A and Q flips each cut a segment and carry flux.
    seed = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    a_flip = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    q_flip = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q <--> E(Q)
        end
    end)
    kids = EnzymeRates._expand_re_to_ss(seed)
    @test length(kids) == 2
    @test Set(kids) == Set([a_flip, q_flip])
end
```

  plus: a decorated seed whose B group is shared (`(E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))`):
  no qualifying chain, so the A, B, P and Q flips all appear (derive each gain; exact set of
  four); and an allosteric K-type variant of the undecorated seed at multiplicity 2 whose flank
  flip is still a child (exact children; derive them with `_hyperbolic_catalysis`; rank-check
  the flank child: the spec measured 80 of 248 such flips raise the rank).

- [ ] **Step 2: Run; expect 4 children today** for the undecorated seed.

- [ ] **Step 3: Implement:**

```julia
"""
Kinetic groups of `m` that are a flank of a qualifying chain whose isomerization is steady
state: an isomerization X1 → X2 alone in its group, where X1 and X2 each have exactly one
other step, each a binding into that form and alone in its group. By the chain lemma the
mechanism sees such a chain only through its flux and enzyme content, four numbers that the
form with both flanks at rapid equilibrium already covers, so flipping a flank adds parameters
and no rate law. An allosteric mechanism has none: at more than one catalytic subunit a
flank's flip can be visible through the two conformations.
"""
function _chain_flank_groups(m::Mechanism)
    groups = steps(m)
    group_of = Dict(s => g for (g, group) in enumerate(groups) for s in group)
    at = Dict{Species, Vector{Step}}()
    for group in groups, s in group, sp in (from_species(s), to_species(s))
        push!(get!(at, sp, Step[]), s)
    end
    flank(s0, x) = begin
        rest = filter(!=(s0), at[x])
        length(rest) == 1 || return nothing
        s = only(rest)
        is_binding(s) && to_species(s) == x && length(groups[group_of[s]]) == 1 ? s : nothing
    end
    out = Set{Int}()
    for group in groups
        length(group) == 1 && is_iso(only(group)) && !is_equilibrium(only(group)) || continue
        s0 = only(group)
        f1, f2 = flank(s0, from_species(s0)), flank(s0, to_species(s0))
        f1 === nothing || f2 === nothing || union!(out, (group_of[f1], group_of[f2]))
    end
    out
end
_chain_flank_groups(::AllostericMechanism) = Set{Int}()
```

  In `_expand_re_to_ss`, compute `flanks = _chain_flank_groups(m)` and add `!(g in flanks) &&`
  to the unit condition. Add one docstring sentence: a flank of a qualifying chain whose
  isomerization is steady state is no unit (`_chain_flank_groups`).

- [ ] **Step 4: Run the focused file.** Every failing exact-children flip test whose parent has
  a qualifying chain is re-derived by hand (remove the flank-flip children; rank-check that each
  removed child had its parent's rank). Re-pin the R4 population at
  `test/test_mechanism_enumeration.jl:9700`: the spec measured `[62, 349, 1134]` with this rule;
  re-pin R6 at 9716 to the measured count. State each pin as the count and what makes it, not as
  a change.

- [ ] **Step 5: Full suite, commit** ("Never flip the flank of a qualifying chain").

---

## Wave 3: the new seeds (spec sections 3 and 4)

### Task 9: MERGE, ELIM, validity and non-degeneracy

**Files:**
- Modify: `src/mechanism_enumeration.jl` (new functions next to the flux predicate; refactor
  `_flux_carrying_steps` to use `_reactant_signs`)
- Test: `test/test_mechanism_enumeration.jl` (new testsets; a `_testhelper_degenerate` probe at
  the top next to `_testhelper_identifiable_rank`)

**Interfaces (produced, used by Tasks 10 and 11):**
- `_reactant_signs(rxn)::Dict{Symbol, Int}`: +1 per substrate name, −1 per product name.
- `_merge_isomerization(groups::Vector{Vector{Step}}, s0::Step)::Vector{Vector{Step}}`.
- `_eliminate_form(groups::Vector{Vector{Step}}, x::Species)::Union{Vector{Vector{Step}}, Nothing}`.
- `_re_turnover_cycle(groups, rxn)::Bool`.
- `_has_vmax(groups, rxn, side::Type)::Bool` with `side` ∈ `(Substrate, Product)`.
- `_chemistry_equilibrates_both_sides(groups, rxn)::Bool`.
- `_degenerate(m::Union{Mechanism, AllostericMechanism})::Bool`.

- [ ] **Step 1: The probe helper** at the top of `test/test_mechanism_enumeration.jl`:

```julia
"""Numeric degeneracy probe (test oracle only): true when the rate at zero products does not
need some substrate (ratio of the rate at 10⁻⁸ of it to the rate at 1 stays above 10⁻³), or
grows without bound as the substrates scale (10⁹ against 10⁶), or the mirror for products."""
function _testhelper_degenerate(m; ndraws = 2)
    em = EnzymeRates.compile_mechanism(m)
    fp = collect(EnzymeRates.fitted_params(em))
    cm = m isa EnzymeRates.Mechanism ? m : EnzymeRates._state_mechanism(m, :A)
    subs = [EnzymeRates.name(s) for s in EnzymeRates.substrates(EnzymeRates.reaction(cm))]
    prods = [EnzymeRates.name(p) for p in EnzymeRates.products(EnzymeRates.reaction(cm))]
    mets = sort!(collect(EnzymeRates._concentration_symbols(cm)))
    rng = MersenneTwister(7)
    for _ in 1:ndraws
        θ = exp.(randn(rng, length(fp)))
        p = NamedTuple{(fp..., :Keq, :E_total)}((θ..., exp(randn(rng)), 1.0))
        v(d) = rate_equation(em, NamedTuple{Tuple(mets)}(Tuple(get(d, x, 0.0) for x in mets)), p)
        for (side, other) in ((subs, prods), (prods, subs))
            base = v(Dict(x => 1.0 for x in side))
            for x in side
                r = v(Dict(y => (y == x ? 1e-8 : 1.0) for y in side)) / base
                (isfinite(r) && abs(r) < 1e-3) || return true
            end
            big = v(Dict(x => 1e9 for x in side)) / v(Dict(x => 1e6 for x in side))
            big > 10 && return true
        end
    end
    false
end
```

  (Regulators and copies stay at 0.)

- [ ] **Step 2: Failing tests** (all on hand-written groups, no init involved). Write the
  undecorated ordered seed and its merge inline:

```julia
@testset "MERGE, ELIM and the seed tests on the ordered bi-bi seed" begin
    ER = EnzymeRates
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    groups(em) = ER.steps(ER.Mechanism(rxn, ER.steps(ER.Mechanism(em))))
    seed = groups(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    iso = only(s for g in seed for s in g if ER.is_iso(s))
    merged = ER._merge_isomerization(seed, iso)
    # The isomerization's product side E(P, Q) takes over E(A, B)'s steps.
    @test Set(s for g in merged for s in g) == Set(s for g in groups(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end) for s in g)
    @test ER._re_turnover_cycle(merged, rxn)        # every step RE: infinite rate
    @test !ER._re_turnover_cycle(seed, rxn)
    variant(ss) = [[ER._with_equilibrium(s, !(ER.bound_metabolite(s) !== nothing &&
                    ER.name(ER.bound_metabolite(s)) in ss)) for s in g] for g in merged]
    clean(gs) = !ER._re_turnover_cycle(gs, rxn) && ER._has_vmax(gs, rxn, ER.Substrate) &&
        ER._has_vmax(gs, rxn, ER.Product) && !ER._chemistry_equilibrates_both_sides(gs, rxn)
    # Track 4 table 4.3: {A, Pˣ}, {Bˣ, Q}, {Bˣ, Pˣ} clean; {A, Q} fails C (B is not
    # needed); {A, Bˣ} has no forward Vmax; {Pˣ, Q} no reverse Vmax; singletons none.
    @test clean(variant([:A, :P])) && clean(variant([:B, :Q])) && clean(variant([:B, :P]))
    @test ER._chemistry_equilibrates_both_sides(variant([:A, :Q]), rxn)
    @test !ER._has_vmax(variant([:A, :B]), rxn, ER.Substrate)
    @test !ER._has_vmax(variant([:P, :Q]), rxn, ER.Product)
    @test all(x -> !clean(variant([x])), (:A, :B, :P, :Q))
    # ELIM of the merged complex: one step E(A) + B → E(Q) + P; all three SS is clean.
    x = ER.to_species(iso)
    tc = ER._eliminate_form(merged, x)
    is_tc(s) = !isempty(ER.consumed(s)) && !isempty(ER.released(s))
    binds(name) = s -> ER.bound_metabolite(s) !== nothing &&
                       ER.name(ER.bound_metabolite(s)) == name
    group_of(pred) = only(g for (g, grp) in enumerate(tc) if any(pred, grp))
    gA, gQ, gTC = group_of(binds(:A)), group_of(binds(:Q)), group_of(is_tc)
    @test length(tc) == 3 && ER.consumed(only(tc[gTC])) == ER.Metabolite[ER.Substrate(:B)] &&
        ER.released(only(tc[gTC])) == ER.Metabolite[ER.Product(:P)]
    ss(which) = [g in which ? ER._with_equilibrium.(grp, false) : grp
                 for (g, grp) in enumerate(tc)]
    @test clean(ss([gA, gQ, gTC]))
    @test !clean(ss([gA, gQ])) && !clean(ss([gA, gTC])) && !clean(ss([gQ, gTC]))
    # The seed itself is not degenerate; the probe agrees on every verdict above.
    @test !ER._degenerate(ER.Mechanism(rxn, seed))
end
```

  Add a second
  testset on the undecorated ping-pong seed (`bi_bi_pp_rxn`; write it inline: `E + A ⇌ E(A)`,
  `E(A) <--> E(P; residual = A - P)`, `E(; residual = A - P) + P ⇌ E(P; residual = A - P)`,
  `E(; residual = A - P) + B ⇌ E(B; residual = A - P)`,
  `E(B; residual = A - P) ⇌ E(Q)`, `E + Q ⇌ E(Q)`): `_degenerate` is true (C fails on the node
  {E(B; res), EQ}) and the probe agrees; flipping B or Q gives a non-degenerate mechanism and
  flipping A or P a degenerate one (spec section 3, the four C examples; probe each). Add a case
  where `_eliminate_form` returns `nothing`: the uni-uni merged complex, whose two steps lead to
  the same form E (a Theorell–Chance step E + S → E + P would join E to itself).

- [ ] **Step 3: Run; expect `UndefVarError`s.**

- [ ] **Step 4: Implement:**

```julia
"""Net uptake sign of each reactant name of `rxn`: +1 for a substrate, −1 for a product."""
function _reactant_signs(rxn::EnzymeReaction)
    rho = Dict{Symbol, Int}()
    for s in substrates(rxn); rho[name(s)] = get(rho, name(s), 0) + 1; end
    for p in products(rxn);   rho[name(p)] = get(rho, name(p), 0) - 1; end
    rho
end
```

  (use it in `_flux_carrying_steps` in place of its inline `rho`), and

```julia
"""
`groups` with the isomerization `s0`, alone in its group, merged onto its product side: `s0`
and its group are removed, and every other step at `from_species(s0)` moves onto
`to_species(s0)`. The last substrate's binding into the substrate side becomes a fused binding
into the product side; the releases stay plain.
"""
function _merge_isomerization(groups::Vector{Vector{Step}}, s0::Step)
    x1, x2 = from_species(s0), to_species(s0)
    move(sp) = sp == x1 ? x2 : sp
    [Step[Step(move(from_species(s)), move(to_species(s)), consumed(s), released(s),
               is_equilibrium(s)) for s in group]
     for group in groups if group != [s0]]
end

"""
`groups` with the form `x` eliminated, or `nothing` when that is not possible: `x` must have
exactly two steps, both bindings into it, Y + L → x and Z + R → x with Y ≠ Z. They become one
Theorell–Chance step Y + L → Z + R, rapid equilibrium, in a group of its own; their groups
keep their other steps.
"""
function _eliminate_form(groups::Vector{Vector{Step}}, x::Species)
    at_x = [s for group in groups for s in group if x in (from_species(s), to_species(s))]
    length(at_x) == 2 && all(s -> is_binding(s) && to_species(s) == x, at_x) ||
        return nothing
    entry, exit = sort(at_x; by = s -> !(bound_metabolite(s) isa Substrate))
    from_species(entry) == from_species(exit) && return nothing
    fused = Step(from_species(entry), from_species(exit), consumed(entry), consumed(exit), true)
    kept = [filter(s -> !(s in at_x), group) for group in groups]
    push!(filter!(!isempty, kept), [fused])
end

"""Whether a cycle of rapid-equilibrium steps of `groups` runs net turnover, which makes the
rate infinite. Each RE step is an edge weighted by the uptake of the reaction's substrates
minus its products (`_reactant_signs`); a turnover cycle has weight n·Σρ², never zero, and a
block of the RE graph holds one iff it is unbalanced (`_block_balanced`)."""
function _re_turnover_cycle(groups::Vector{Vector{Step}}, rxn::EnzymeReaction)
    rho = _reactant_signs(rxn)
    idx = Dict{Species, Int}()
    vertex(sp) = get!(idx, sp, length(idx) + 1)
    edges = Tuple{Int, Int}[]; weights = Int[]
    for group in groups, s in group
        is_equilibrium(s) || continue
        push!(edges, (vertex(from_species(s)), vertex(to_species(s))))
        push!(weights, sum(get(rho, name(m), 0) for m in consumed(s); init = 0) -
                       sum(get(rho, name(m), 0) for m in released(s); init = 0))
    end
    isempty(edges) && return false
    block = _edge_blocks(length(idx), edges)
    any(b -> !_block_balanced(edges, weights, findall(==(b), block)), unique(block))
end

"""Whether the rate has a maximum as the metabolites of `side` (`Substrate` or `Product`)
grow: turning rapid equilibrium every steady-state step that takes up or gives off one of them
leaves no rapid-equilibrium turnover cycle (track 4's condition V)."""
function _has_vmax(groups::Vector{Vector{Step}}, rxn::EnzymeReaction, side::Type)
    names = Set(name(x) for x in (side === Substrate ? substrates(rxn) : products(rxn)))
    touches(s) = any(m -> name(m) in names, Iterators.flatten((consumed(s), released(s))))
    !_re_turnover_cycle([[touches(s) ? _with_equilibrium(s, true) : s for s in group]
                         for group in groups], rxn)
end

"""
Whether some chemistry node of `groups` is left both by a rapid-equilibrium step that releases
a substrate and by one that releases a product (track 4's condition C). A chemistry node is a
merged complex (a form entered by a fused binding of a substrate) with every form joined to it
by RE isomerizations, a set of forms joined to each other by RE isomerizations, or an RE
Theorell–Chance step. An RE step leaves a node through a node form F and releases M when
traversing it away from F gives M off: F is its `to_species` and M is consumed, or F is its
`from_species` and M is released. Then the chemistry sits in rapid equilibrium with both sides
and some reactant is not needed at zero products or substrates.
"""
function _chemistry_equilibrates_both_sides(groups::Vector{Vector{Step}},
                                            rxn::EnzymeReaction)
    subs = Set(name(s) for s in substrates(rxn))
    prods = Set(name(p) for p in products(rxn))
    re = [s for group in groups for s in group if is_equilibrium(s)]
    re_iso = filter(is_iso, re)
    function node(start)
        out = Set{Species}(start); frontier = collect(start)
        while !isempty(frontier)
            f = pop!(frontier)
            for s in re_iso, (a, b) in ((from_species(s), to_species(s)),
                                        (to_species(s), from_species(s)))
                a == f && !(b in out) && (push!(out, b); push!(frontier, b))
            end
        end
        out
    end
    releases(n, names) = any(re) do t
        (to_species(t) in n && !(from_species(t) in n) &&
            any(m -> name(m) in names, consumed(t))) ||
        (from_species(t) in n && !(to_species(t) in n) &&
            any(m -> name(m) in names, released(t)))
    end
    merged = [to_species(s) for group in groups for s in group
              if is_binding(s) && _is_chemistry(s) && name(bound_metabolite(s)) in subs]
    nodes = vcat([node([x]) for x in merged],
                 [node([from_species(s), to_species(s)]) for s in re_iso])
    any(n -> releases(n, subs) && releases(n, prods), nodes) && return true
    any(re) do t
        !isempty(consumed(t)) && !isempty(released(t)) &&
            (any(m -> name(m) in subs, consumed(t)) && any(m -> name(m) in prods, released(t)) ||
             any(m -> name(m) in prods, consumed(t)) && any(m -> name(m) in subs, released(t)))
    end
end

"""Whether `m` is degenerate by track 4's structural conditions: no maximal rate in one
direction (`_has_vmax`), or chemistry in rapid equilibrium with both sides
(`_chemistry_equilibrates_both_sides`). Read on `steps(m)`, an allosteric mechanism's
active-state graph. The beam expands a degenerate seed without fitting it."""
function _degenerate(m::Union{Mechanism, AllostericMechanism})
    groups, rxn = steps(m), reaction(m)
    !(_has_vmax(groups, rxn, Substrate) && _has_vmax(groups, rxn, Product)) ||
        _chemistry_equilibrates_both_sides(groups, rxn)
end
```

  Check `_block_balanced`'s signature (`edges, weights, es::Vector{Int}`) and that `_edge_blocks`
  numbers blocks `1:nblocks`. A merged base's first `_re_turnover_cycle` call short-circuits
  before `_flux_carrying_groups` (which assumes no RE turnover cycle).

- [ ] **Step 5: Run; every assertion of Step 2 passes; cross-check every hand case with
  `_testhelper_degenerate` in the same testsets** (`@test _testhelper_degenerate(Mechanism(rxn,
  gs)) == !clean(gs)` for each constructible `gs`; a merged base with an RE turnover cycle is
  not constructible, skip it).

- [ ] **Step 6: Full suite, commit** ("Test seed candidates for turnover, maximal rates and
  equilibrated chemistry").

### Task 10: `init_mechanisms` emits the merged and Theorell–Chance variants

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`init_mechanisms` and its docstring; new
  `_seed_variants`)
- Modify: `test/fixtures/phase2_init_golden.txt` (regenerated deliberately)
- Test: `test/test_mechanism_enumeration.jl` (new testsets; pins at 1411, 9700, 9716, the
  golden test near 10200), `test/test_identify_rate_equation.jl:736-739` (the bi-bi distinct
  equation count)

**Interfaces:**
- Consumes: Task 9's functions, `_minimal_gaining_sets`, `_flip_group_to_ss`,
  `_all_steady_state`, `_flux_carrying_groups`, `_bottomless_re_segment`, `_with_equilibrium`.
- Produces: `_seed_variants(m::Mechanism)::Vector{Mechanism}`; `init_mechanisms(rxn)` returns
  today's seeds, in today's order, followed by the variants not among them.

- [ ] **Step 1: Failing tests, exact variant sets** (each expected variant hand-derived with the
  Task 9 predicates in its comment, then checked: fitted = rank by the oracle and
  `!_testhelper_degenerate`):
  - **Uni-uni seed** (`E + S ⇌ E(S)`, `E(S) <--> E(P)`, `E + P ⇌ E(P)` on `uni_uni_rxn`):
    `_seed_variants` is empty (the {S, P} merge is the lumping twin, singletons have no Vmax,
    and ELIM would join E to itself).
  - **Undecorated ordered seed** (Task 9's `seed` on a bi-bi reaction with atoms): exactly
    three variants, {A, Pˣ} (`E + A <--> E(A)`, `E(A) + B ⇌ E(P, Q)`, `E(Q) + P <--> E(P, Q)`,
    `E + Q ⇌ E(Q)`), {Bˣ, Q} (`E + A ⇌ E(A)`, `E(A) + B <--> E(P, Q)`, `E(Q) + P ⇌ E(P, Q)`,
    `E + Q <--> E(Q)`), and its Theorell–Chance seed (`E + A <--> E(A)`,
    `E(A) + B <--> E(Q) + P`, `E + Q <--> E(Q)`); each 5 fitted, rank 5.
  - **Decorated ordered seed** with the B group shared (`(E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))`):
    exactly four, {Bˣ, Pˣ} (no twin: the B group is shared), {A, Pˣ}, {Bˣ, Q}, and the
    Theorell–Chance seed whose B group keeps `E(Q) + B ⇌ E(B, Q)` alone (6 fitted); derive and
    rank each.
  - **Undecorated ping-pong seed** (Task 9's ping-pong seed on `bi_bi_pp_rxn`): exactly four,
    PP-AQ ({A, Q} SS), PP-BP ({B, P} SS), and the two half-Theorell–Chance seeds (X1 eliminated
    with B and Q SS; X2 eliminated with A and P SS); each 5/5.
  - **Undecorated ordered/random seed** (A first, random products: `E + A ⇌ E(A)`,
    `E(A) + B ⇌ E(A, B)`, `E(A, B) <--> E(P, Q)`, `(E(Q) + P ⇌ E(P, Q), E + P ⇌ E(P))`,
    `(E(P) + Q ⇌ E(P, Q), E + Q ⇌ E(Q))`): the track's expectation is {Bˣ, Pˣ}, {Bˣ, Qˣ} at 5
    and {A, Pˣ, Qˣ} at 6; derive each candidate with V and C by hand and STOP if the
    derivation disagrees.
  Then the aggregate tests (each states the count it pins and what makes it):
  - `init_mechanisms(bi_bi_pp_rxn)`: its first 62 are today's seeds (each holds an
    isomerization), every later one holds none, and the total is the measured count (spec
    estimate about 264); `init_mechanisms(bi_bi_rxn)` likewise from 55.
  - Review Focus: `init_mechanisms(uni_bi_rxn)` gains exactly one merged variant per undecorated
    ordered topology (n + m − 2 = 1); pin the total.
  - Review Focus: `@test (@elapsed EnzymeRates.init_mechanisms(ter_ter_rxn)) < 60` after a warm
    call; report the count. If it exceeds 60 s, STOP and report to the controller.
  - Review Focus: `seed_mechanisms` on a bi-bi reaction with `dead_end_inhibitors: I` (required)
    returns seeds, at least one of them holding a fused binding.

- [ ] **Step 2: Run; expect failures** (`_seed_variants` undefined; counts unchanged).

- [ ] **Step 3: Implement:**

```julia
"""
The merged and Theorell–Chance variants of the seed `m`. Merging every isomerization of `m`
onto its product side (`_merge_isomerization`) and setting every step at rapid equilibrium
gives the merged base; eliminating one merged complex with two steps (`_eliminate_form`) gives
a Theorell–Chance base. For each base, every inclusion-minimal set of its groups whose flip to
steady state gives a valid, flux-carrying, non-degenerate candidate is a variant: no
rapid-equilibrium turnover cycle (`_re_turnover_cycle`), no bottomless segment, every
steady-state group carrying flux, a maximal rate both ways (`_has_vmax`) and no chemistry in
equilibrium with both sides (`_chemistry_equilibrates_both_sides`). A merged variant whose every
merged complex has both its steps steady state, each alone in its group, is skipped: the
unmerged form with rapid-equilibrium flanks has its family at the same count.
"""
function _seed_variants(m::Mechanism)
    rxn = reaction(m)
    isos = [s for group in steps(m) for s in group if is_iso(s)]
    merged = foldl(_merge_isomerization, isos; init = steps(m))
    merged = [_with_equilibrium.(group, true) for group in merged]
    complexes = [to_species(s) for s in isos]
    bases = Tuple{Vector{Vector{Step}}, Bool}[(merged, true)]
    for x in complexes
        b = _eliminate_form(merged, x)
        b === nothing || push!(bases, (b, false))
    end
    admissible(gs) = !_re_turnover_cycle(gs, rxn) && _bottomless_re_segment(gs) === nothing &&
        (flux = _flux_carrying_groups(gs, rxn);
         all(g -> is_equilibrium(first(gs[g])) || flux[g], eachindex(gs))) &&
        _has_vmax(gs, rxn, Substrate) && _has_vmax(gs, rxn, Product) &&
        !_chemistry_equilibrates_both_sides(gs, rxn)
    lumping_twin(gs) = all(complexes) do x
        at_x = [(s, group) for group in gs for s in group
                if x in (from_species(s), to_species(s))]
        length(at_x) == 2 && all(((s, group),) -> !is_equilibrium(s) && length(group) == 1, at_x)
    end
    variants = Mechanism[]
    for (base, is_merged) in bases
        flux = _flux_carrying_groups(_all_steady_state(base), rxn)
        units = [g for g in eachindex(base) if flux[g]]
        flipped(sel) = foldl((gs, u) -> _flip_group_to_ss(gs, units[u]), sel; init = base)
        for sel in _minimal_gaining_sets(sel -> admissible(flipped(sel)), _ -> 1:length(units))
            gs = flipped(sel)
            is_merged && lumping_twin(gs) && continue
            push!(variants, Mechanism(rxn, gs))
        end
    end
    variants
end
```

  and in `init_mechanisms`, after building `mechs` exactly as today:

```julia
    variants = Mechanism[]
    for m in unique(mechs)
        append!(variants, _seed_variants(m))
    end
    vcat(mechs, filter!(v -> !(v in mechs), unique!(variants)))
```

  `mechs` is returned unchanged as the prefix (today's output exactly, duplicates included if
  any), so every `first(init_mechanisms(rxn))` and every caller's `unique!` behave as before.
  Run `_assert_atom_conserving` on each variant (MERGE and ELIM conserve atoms; assert it).
  Rewrite the `init_mechanisms` docstring: today's seeds first, then their variants; mixed
  parameter counts.

- [ ] **Step 4: Run the new testsets; then regenerate the golden fixture deliberately** in a
  scratch script (same key function as the golden test), write it to
  `test/fixtures/phase2_init_golden.txt`, and add to the golden test an assertion that the 55
  keys of today's seeds are a subset (hard-code nothing: compute today's keys from the
  isomerization-holding prefix). Re-pin `test/test_mechanism_enumeration.jl:1411` (55 → the
  `bi_bi_rxn` total), the population pins at 9700 and 9716 (measure), and the distinct-equation
  count at `test/test_identify_rate_equation.jl:739` (measure: compile every bi-bi seed; the
  comment states the count and that the 8 ordered/random and 2 Theorell–Chance duplicate
  families are structurally distinct and so render distinct equations).

- [ ] **Step 5: Full suite (expect it to run longer: the identify test compiles every bi-bi
  seed), commit** ("Seed the enumeration with merged and Theorell–Chance variants").

### Task 11: The beam expands degenerate seeds without fitting them

**Files:**
- Modify: `src/identify_rate_equation.jl` (`_beam_search` base tier; a new `_base_tier`)
- Test: `test/test_identify_rate_equation.jl`

**Interfaces:**
- Consumes: `_degenerate` (Task 9), `_expand_parent`.
- Produces: `_base_tier(mechs::Vector, rxn)::Tuple{Vector{Union{Mechanism, AllostericMechanism}}, Vector{FitFailure}}`.

- [ ] **Step 1: Failing test:**

```julia
@testset "_base_tier expands degenerate seeds instead of fitting them" begin
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    seeds = unique!(collect(EnzymeRates.init_mechanisms(rxn)))
    degenerate = filter(EnzymeRates._degenerate, seeds)
    @test length(degenerate) == 7                  # the seven ping-pong seeds
    base, failures = EnzymeRates._base_tier(seeds, rxn)
    @test isempty(failures)
    @test !any(m -> m in base, degenerate)
    @test all(m -> m in base, filter(!EnzymeRates._degenerate, seeds))
    kids = EnzymeRates.expand_mechanisms(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[degenerate...], rxn)
    @test all(k -> k in base, kids)
    @test length(base) == length(unique!(vcat(filter(!EnzymeRates._degenerate, seeds), kids)))
end
```

  Review Focus: add a test that `_rows_to_dataframe` accepts a row whose `parent_n_params` is
  `missing` and whose mechanism is a child of an expand-only seed (build the row with the same
  fields `_fit_batch` writes for a base-tier entry; grep `parent_n_params` in
  `src/identify_rate_equation.jl:280-300`).

- [ ] **Step 2: Run; expect `UndefVarError: _base_tier`.**

- [ ] **Step 3: Implement:**

```julia
"""
The mechanisms the beam fits first: every seed of `mechs` that is not degenerate
(`_degenerate`), and the children of every one that is. A degenerate seed's law ignores a
substrate or never saturates, so it is not worth a fit, but it is the only parent of
mechanisms that are not degenerate (the ping-pong seeds' children). Expansion errors are
returned as failures, one per seed, as `_expand_parent` records them.
"""
function _base_tier(mechs::Vector, rxn::EnzymeReaction)
    base = Union{Mechanism, AllostericMechanism}[m for m in mechs if !_degenerate(m)]
    failures = FitFailure[]
    for m in mechs
        _degenerate(m) || continue
        kids, failure = _expand_parent(m, rxn)
        failure === nothing || push!(failures, failure)
        append!(base, kids)
    end
    unique!(base), failures
end
```

  In `_beam_search`, replace `base = (…) ? unique!(collect(init_mechanisms(…))) :
  unique!(collect(seed_mechanisms(…)))` by computing the seeds the same way and then
  `base, base_expand_failures = _base_tier(seeds, prob.reaction)`, and append
  `base_expand_failures` to `base_failures` after `_fit_batch` (so they are written as failure
  rows in `initial_mechanisms.csv`). Update the "Base tier" docstring/comment and the
  `_required_regulators` docstring sentence that describes the seed.

- [ ] **Step 4: Run the identify test file focused**
  (`julia --project --heap-size-hint=2500M -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_identify_rate_equation.jl")'`);
  every beam test must still pass (they run uni-uni, which has no degenerate seed). Full suite;
  commit ("Expand degenerate seeds without fitting them").

---

## Wave 4: record and documentation

### Task 12: Regression record and timing for C

**Files:**
- Create: `docs/superpowers/specs/2026-10-03-merged-seeds-regression.md`
- Scripts in the scratchpad only (describe them in the record).

- [ ] **Step 1:** On a detached worktree at the opening wave's last commit (Task 4) and on the
  tip, run in turn (one Julia process at a time):
  - `rate_equation_string` for `MECHANISM_TEST_SPECS` and R1–R6 to depth 2 (flip, split,
    dead-end), keyed by `_forward_sides` per group: identical for every mechanism without a
    fused or Theorell–Chance step; fused fixtures equal under the name map (evaluate both
    strings at 5 random points with a small expression interpreter, not by compiling).
  - `init_mechanisms` counts per reaction: R1–R6, `uni_bi_rxn`, `ter_ter_rxn`; for every new
    seed: derives, `_testhelper_check_against_mass_action` logic on a sample of 30, numeric probe
    non-degenerate, fitted and rank (oracle).
  - V and C against the probe on every candidate the seed search tests for R4 and R6 (log each
    `admissible` verdict): a single disagreement stops the task; report it.
  - The flip rule: every mechanism the tip lacks relative to the worktree, at R4 levels 0–3, has
    its RSR twin (both flanks RE) in the tip's population; phantom counts before and after at
    levels 0–2.
  - Per-level counts for R1–R6 to depth 2 (R6 level 2 sampled with `MersenneTwister(20261003)`,
    4,000 mechanisms), the level-3 phantom fraction on R4, the allosteric sample of A's record.
  - Timing: R6 levels 0–2 enumeration, per emitted mechanism, tip against worktree (limit +20%);
    `init_mechanisms` and `seed_mechanisms` (with `dead_end_inhibitors: I` required) on
    `ter_ter_rxn`.
- [ ] **Step 2:** Write the record (method, tables, every disagreement or limit exceeded named).
  If any limit is exceeded or any check fails, STOP and report before committing.
- [ ] **Step 3:** Commit the record ("Record the regression with merged seeds and the flip rule").

### Task 13: Documentation and CLAUDE.md

**Files:**
- Modify: `docs/src/identify/enumeration_engine.md`, `docs/src/deriving/kinetic_groups.md:22`,
  `docs/src/deriving/re_vs_ss.md`, `docs/src/identify/model_selection.md`,
  `docs/src/identify/combinatorics.md` (only if it quotes seed counts), `docs/src/developer.md`,
  `.claude/CLAUDE.md` (enumeration-test rule 3),
  `docs/superpowers/specs/2026-09-28-identifiable-enumeration-findings-and-plan.md` (section
  "C"), docstrings of `init_mechanisms` and `expand_mechanisms` if Tasks 10/6 left anything
  stale.
- Read first: the spec's "Documentation" section, which lists every change and gives the
  approved CLAUDE.md wording verbatim. Apply the elements-of-style rules (active voice,
  positive form, omit needless words); Denis asks for them on all documentation.

- [ ] **Step 1:** `enumeration_engine.md`: the seed construction (merged and Theorell–Chance
  variants, V and C, skipped lumping twins, mixed counts), the flip rule (and why allosteric
  flank flips stay), dead-end sites on Theorell–Chance steps; rewrite the paragraphs that call
  chemistry the isomerization step (around lines 36–42 and 313–319) around `_is_chemistry`;
  reflect Task 3/4 (one twin predicate; the split's child filter).
- [ ] **Step 2:** `kinetic_groups.md:22`: a group's steps take up and give off the same
  metabolites and share one flag. `re_vs_ss.md`: an RE label on a flank of a steady-state chain
  is an apparent (Michaelis-type) constant; same rate laws as the all-SS chain, fewer
  parameters. `model_selection.md`: degenerate seeds expanded, not fitted.
  `developer.md`: stoichiometric step kinds, `_is_chemistry`, seed construction.
- [ ] **Step 3:** CLAUDE.md rule 3: replace the rule with the spec's approved wording, verbatim.
- [ ] **Step 4:** The findings document's section C: one paragraph pointing to the spec and
  listing the decisions (stoichiometric bindings; flip rule only; expand-only degenerate seeds;
  duplicates kept; B simplified; gauge ruling confirmed).
- [ ] **Step 5:** Build the docs if the docs environment is available
  (`julia --project=docs docs/make.jl`), else check every `@ref` and code reference by grep;
  full suite (doctests run in it if configured); commit ("Document merged seeds, fused steps
  and the flip rule").
