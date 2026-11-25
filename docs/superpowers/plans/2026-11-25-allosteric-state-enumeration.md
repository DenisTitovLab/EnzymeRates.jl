# Allosteric state enumeration: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development
> (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps
> use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the enumerator build every allosteric version of a mechanism that the
package accepts, with the four additions the design lists.

**Architecture:** `_expand_to_allosteric` emits every set of `:OnlyA` kinetic groups;
`_expand_change_allo_state` tags minimal gaining sets `:NonequalAI`. One predicate,
`_enumerator_returns`, says whether the enumerator returns a candidate, and the two
places that return candidates apply it. All changes are inside `src/Enumeration/`.

**Tech stack:** Julia 1.12 (tests), Julia 1.13 for `test/test_type_stability.jl`.

**Spec:** `docs/superpowers/specs/2026-11-24-allosteric-state-enumeration-design.md`. Read
it whole before a task; it defines the states, the functions and the terms "version",
"minimal gaining set", "irreducible" and "kind (a), (b), (c)". Those kind labels are the
spec's; do not put them in code or docstrings.

## Global Constraints

- **Where to work.** Only in `/home/denis.linux/claude-scratchpad/EnzymeRates-allosteric-states`
  (branch `allosteric-state-enumeration`). Never read-write
  `/home/denis.linux/.julia/dev/EnzymeRates`: another session uses it.
- **One Julia process at a time** on this machine (7.7 GB, no swap). Run `pgrep -a julia`
  first and wait if one is running.
- **Focused test run** (from the worktree root):
  `julia --project=. --heap-size-hint=2500M -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/testhelpers.jl"); include("test/<path>.jl")'`.
  `test/Enumeration/test_expansion_moves.jl` takes several minutes. Read all output; it
  must be pristine.
- **Project rules.** Follow `.claude/CLAUDE.md`: 92-character lines, docstrings that say
  what a caller needs in the code's own names, comments only for what the code cannot say,
  no "new/old/legacy" in names or comments, every helper under `test/` prefixed
  `_testhelper_`.
- **Enumeration tests** (`.claude/CLAUDE.md`, "Enumeration-engine tests"): write the
  mechanism inline with `@enzyme_mechanism`, hand it over as `_testhelper_candidate(m)`,
  assert `length(children) == n` and `Set(children) == Set(expected)` with every expected
  child written out, and name forms as the enumerator does.
- **Expected children come from the spec, not from the move's output.** To learn which
  sets gain on a fixture, call `Derivation.fitted_parameter_count` and
  `Enumeration._enumerator_keeps` on versions you write out. Never paste what the move
  returned.
- **Changing an existing test.** A testset's expected set or pinned count changes only
  when the spec changes the behaviour it pins. The task report lists every such change:
  file, testset name, old value, new value, and the spec sentence behind it. A testset
  that pins a removed rule is deleted and listed. Any other failing assertion is reported
  to the controller and not edited.
- **TDD.** Write each new test first and watch it fail for the stated reason.
- **No public name changes.** `Enumeration` calls only public names of `Derivation`.
- **Developer page tables.** After a change to which names one file of `src/Enumeration/`
  takes from another, run `julia --project=. test/update_module_map.jl`, then the focused
  run of `test/test_boundaries.jl`, and commit the regenerated `docs/src/developer.md`
  with the task.
- **Commits** end with these two lines:
  `Co-Authored-By: Claude Fable 5.2 <noreply@anthropic.com>` and
  `Claude-Session: https://claude.ai/code/session_015DwkaAKoGZ5ULMU7mDCAhF`.
  Never skip a pre-commit hook. Never push.

## Review Focus

1. **A reaction with no competitive inhibitor.** No child of any move may hold an `:OnlyI`
   kinetic group. Pinned in Task 7.
2. **A parent that is not hyperbolic.** With allowed oligomeric states `(1, 2)` it gets
   children at one subunit only; with `2` alone it gets none and no error. Pinned by the
   existing testsets at `test_expansion_moves.jl` 5202, 5224, 5259, which must keep
   passing in Task 5.
3. **A parent with nothing left to change.** `_expand_change_allo_state` on a parent whose
   every group is `:NonequalAI` returns nothing, and on a ten-group parent it returns
   within the time bound of Task 4.
4. **A regulator with a declared effect.** With `R::Activator` required, no starting
   mechanism holds `R` tagged `:OnlyI`. Pinned in Task 6.
5. **Determinism.** Two calls of `candidate_mechanisms(reaction, parent)` return equal
   vectors, in the same order. Pinned in Task 5.

---

### Task 1: Starting mechanisms from the base tier

**Files:**
- Modify: `src/Enumeration/starting_mechanisms.jl` (`seed_mechanisms` at 554,
  `candidate_mechanisms(reaction)` at 621, docstrings)
- Test: `test/Enumeration/test_starting_mechanisms.jl`

**Interfaces:**
- Produces: `seed_mechanisms(rxn, required_allo, required_comp)` with its signature
  unchanged; its search starts from `_base_tier(unique(init_mechanisms(rxn)), rxn)`.
  `candidate_mechanisms(reaction)` calls `_base_tier` only when nothing is required and no
  longer calls it on the result of `seed_mechanisms`.

- [ ] **Step 1: Write the failing tests** in `test_starting_mechanisms.jl`:
  - `"seed_mechanisms starts from the base tier"`: for
    `A[CX] + B[N] → P[C] + Q[NX]` with `allosteric_regulators: R` and
    `oligomeric_state: (1, 2)`, `length(candidate_mechanisms(rxn)) == 7126` (an aggregate
    pin; say so in a comment), and no returned candidate satisfies
    `Enumeration._degenerate`.
  - In the same testset, for the same reaction with `competitive_inhibitors: I` in place
    of the regulator: `length == 1548`.
- [ ] **Step 2: Run the file; the first count is 6,428, so the test fails.**
- [ ] **Step 3: Implement.** Update the docstrings of `seed_mechanisms` and `_base_tier`.
- [ ] **Step 4: Resolve the existing testsets** under the test-change rule. Expected to
  need work: `"seed_mechanisms wave-parallel equivalence"` (its serial search must start
  from the base tier too), `"seed_mechanisms drops a seed the guard rejects"` (count
  1,638; its three named mechanisms `half1`, `half2`, `both`: if `both ∈ seeds` or either
  `∉` no longer holds, stop and report to the controller), the `_base_tier` testsets at
  2220–2409 and `test_expansion_moves.jl` 12366 and 12500 if they pass allosteric
  mechanisms that only the old order produced.
- [ ] **Step 5: Run** `test/Enumeration/test_starting_mechanisms.jl`,
  `test/Enumeration/test_expansion_moves.jl`, `test/Mechanisms/test_mechanism_type.jl` and
  `test/test_oracles.jl`. Re-pin population counts that changed, and list them.
- [ ] **Step 6: Commit.**

### Task 2: The hyperbolic rule where candidates are returned

**Files:**
- Modify: `src/Enumeration/candidate_mechanism_type.jl` (after `_enumerator_keeps`, 666),
  `src/Enumeration/expansion_moves.jl` (65, 218, 910, 1026 and the three docstrings),
  `src/Enumeration/starting_mechanisms.jl` (586)
- Test: `test/Enumeration/test_candidate_mechanism_type.jl`,
  `test/Enumeration/test_expansion_moves.jl`

**Interfaces:**
- Consumes: Task 1 (`_base_tier` now runs `_expand_re_to_ss` on plain mechanisms only).
- Produces: `_enumerator_returns(c::CandidateMechanism)::Bool`, true when
  `_hyperbolic_where_required(c)` and `_enumerator_keeps(c.mechanism)`. Tasks 4 and 5 call
  it. `candidate_mechanisms(reaction, parent)` and `seed_mechanisms` filter with it where
  they call `_enumerator_keeps` today. `_expand_re_to_ss`, `_expand_to_allosteric` and
  `_expand_change_allo_state` no longer call `_hyperbolic_where_required`.

- [ ] **Step 1: Write the failing test** `"_enumerator_returns"` in
  `test_candidate_mechanism_type.jl`: on the written-out child `k_pe` of the testset at
  `test_expansion_moves.jl` 4506 (active conformation hyperbolic, inactive not, two
  subunits) it is false; on the same mechanism at one subunit it equals
  `_enumerator_keeps`; on a plain mechanism it equals `_enumerator_keeps`.
- [ ] **Step 2: Run; it fails with `UndefVarError`.**
- [ ] **Step 3: Implement** the predicate and move the filter. `_expand_to_allosteric`
  keeps its lines 883–886 (a parent that fails `_hyperbolic_catalysis` gets oligomeric
  state 1 only).
- [ ] **Step 4: Carry the testsets over** (spec, "Tests", "Hyperbolic rule"):
  - Flip: 8829 and 8907 assert through `_testhelper_children(rxn, [parent])` instead of
    `_testhelper_move(_expand_re_to_ss, …)`; 502's "every child is hyperbolic" likewise.
  - `_expand_to_allosteric`: 4228 and 4506, and `_expand_change_allo_state`: 6542. Each
    asserts that `_hyperbolic_where_required` is false on the written-out child and that
    `_testhelper_children` does not return it. Their expected sets of the bare move are
    rewritten in Tasks 4 and 5; here only remove the assumption that the move filters.
  - `test_candidate_mechanism_type.jl` 1651–1729 line 1719: must still hold, since the
    population comes through `candidate_mechanisms`. If it fails, report.
- [ ] **Step 5: Run the two test files, `test_starting_mechanisms.jl`, and
  `test/test_boundaries.jl` after regenerating the developer tables. Commit.**

### Task 3: Chain flanks in allosteric mechanisms

**Files:**
- Modify: `src/Enumeration/expansion_moves.jl` (`_chain_flank_groups` 221–260, docstring
  of `_expand_re_to_ss`)
- Test: `test/Enumeration/test_expansion_moves.jl`

**Interfaces:**
- Produces: `_chain_flank_groups(am)` returns the flanks of every qualifying chain,
  whether or not `am` is allosteric.

- [ ] **Step 1: Write the failing test** `"an allosteric mechanism's chain flanks are not
  flipped"`: the bi-bi mechanism of the testset at line 502 (five groups, a steady-state
  isomerization between `EA`-type forms with one binding each side), at one subunit and at
  two. `_chain_flank_groups` returns the two flank groups, and
  `_testhelper_move(_expand_re_to_ss, parent)` has no child whose flank step is at steady
  state. Write the expected children out.
- [ ] **Step 2: Run; it fails because the function returns an empty set.**
- [ ] **Step 3: Delete line 233 and rewrite the docstring's last sentence:** the rule
  covers allosteric mechanisms, where it is a modeling choice and not a proof.
- [ ] **Step 4: Give the testsets at 502, 793, 1056, 1147 and 12500 their new expected
  sets** (each loses the children that flip a flank; 502 also loses
  `@test isempty(_chain_flank_groups(am))`). The comment at
  `test/test_identify_rate_equation.jl` 1945–1981 cites the old behaviour: correct the
  comment only.
- [ ] **Step 5: Run the file. Commit.**

### Task 4: `_expand_change_allo_state`

**Files:**
- Modify: `src/Enumeration/expansion_moves.jl` (976–1027)
- Test: `test/Enumeration/test_expansion_moves.jl`

**Interfaces:**
- Consumes: `_enumerator_returns(c)` (Task 2), `_minimal_gaining_sets(gains, partners)`
  (line 426), `_redundant_copy_groups(c)`, `Derivation.fitted_parameter_count(m)`.
- Produces: `_expand_change_allo_state(c)` with its signature unchanged:
  - one child per minimal set of kinetic groups not yet `:NonequalAI` whose child, with
    those groups tagged `:NonequalAI`, passes `_enumerator_returns`, has no redundant
    inhibitor, and has a higher fitted count than `c`. A set that fails is extended by any
    other group not yet `:NonequalAI` (`partners = _ -> those groups`);
  - then one child per regulatory ligand not yet `:NonequalAI`, kept only when its fitted
    count exceeds `c`'s.

- [ ] **Step 1: Write the failing tests**, each with exact children, on the six-step
  ping-pong with both chemistry steps at steady state (spec, "A worked case") at one
  subunit:
  - `"change_allo_state tags an :OnlyA half-reaction :NonequalAI alone"`: parent with the
    first half-reaction `:OnlyA`. Among the children is the one with that half
    `:NonequalAI` and nothing else changed.
  - `"change_allo_state tags a pair of bindings together"`: parent with A and P `:OnlyA`.
    One child tags both `:NonequalAI`; no child tags A alone.
  - `"change_allo_state never tags one binding alone"`: parent with the first
    half-reaction `:OnlyA` and every binding `:EqualAI`. No child changes exactly one
    binding group.
  - `"change_allo_state on a parent with every group :NonequalAI"`: no child.
  - `"change_allo_state keeps a regulatory ligand's child that adds a constant"`: the
    uni-uni mechanism with `R::OnlyA` at a site; the ligand's `:NonequalAI` child is
    present and has one more fitted constant.
  - `"change_allo_state returns in time on a ten-group parent"`: a ter-ter ordered
    mechanism with one binding `:OnlyA`; `@elapsed` of the move's second call is under
    10 s. If it is not, stop and report the time and the number of sets tried.
- [ ] **Step 2: Run; the pair test and the half-reaction test fail.**
- [ ] **Step 3: Implement.** Order the tests inside `gains` by cost.
- [ ] **Step 4: Carry over the 18 testsets of the move** (inventory: 1627, 1747, 1845,
  1936, 3679, 6274, 6317, 6363, 6424, 6453, 6526, 6542, 9769, 10418, 10955, 11661, 11691,
  11815). Those that pin "one binding group at a time" or "all chemistry groups together"
  pin removed rules: delete and list them, or give them the new expected set when the
  fixture has at most four groups. 6317 pins a ligand child with a fitted-count change of
  0: under the new rule that child is not emitted; re-pin and list.
- [ ] **Step 5: Run the file. Commit.**

### Task 5: `_expand_to_allosteric`

**Files:**
- Modify: `src/Enumeration/expansion_moves.jl` (790–911: delete `_valid_tagging`,
  `_stops_inactive_turnover`, `_tagging_bases`; rewrite the move),
  `src/Enumeration/starting_mechanisms.jl` (`_is_seed_node` 608–619 and the docstring of
  `seed_mechanisms`)
- Test: `test/Enumeration/test_expansion_moves.jl`,
  `test/Enumeration/test_starting_mechanisms.jl`

**Interfaces:**
- Consumes: `_enumerator_returns`, `_minimal_gaining_sets`, `_chemistry_groups(c)`,
  `_candidate(c, states, oligomeric_state, sites)`.
- Produces: `_expand_to_allosteric(c, rxn)` with its signature unchanged. For each allowed
  oligomeric state `cn` (1 only when `c` fails `_hyperbolic_catalysis`), in this order:
  1. one child per nonempty set of kinetic groups tagged `:OnlyA`, the rest `:EqualAI`, no
     regulatory site, sets in `_subsets` order;
  2. for each such child whose `:OnlyA` groups are all chemistry groups, and each
     allosteric regulator of `rxn`, two children with the regulator at a new site of
     multiplicity `cn`, tagged `:OnlyA` and `:OnlyI`;
  3. for each minimal gaining set (a set gains when its child, tagged `:NonequalAI` with
     the rest `:EqualAI`, passes `_enumerator_returns` and has a fitted count above
     `fitted_parameter_count(c.mechanism) + 1`), the child tagged `:NonequalAI`, emitted
     only when the same set tagged `:OnlyA` at `cn` fails `_enumerator_returns`.

  The move does not call `_enumerator_keeps` on the children of item 1. It still errors
  on an allosteric `c`.
- Produces: `_is_seed_node` no longer rejects a candidate for a `:NonequalAI` kinetic
  group; it still rejects a `:NonequalAI` regulatory ligand.

- [ ] **Step 1: Write the failing tests**, exact children of the bare move:
  - `"to_allosteric emits every :OnlyA set of a uni-uni mechanism"`: the uni-uni mechanism
    of the testset at 4707, no regulator, one subunit. Seven children from item 1, plus
    the item 3 children you derive with `fitted_parameter_count` and `_enumerator_keeps`.
  - `"to_allosteric on a three-step ping-pong"`: three groups, seven item 1 children, and
    item 3 as derived.
  - `"to_allosteric pairs a chemistry-only :OnlyA set with each regulator"`: the uni-uni
    mechanism with `uni_uni_allo_reg`; the chemistry-only child appears three times (bare,
    `R::OnlyA`, `R::OnlyI`); no child with an `:OnlyA` binding group has a site; no child
    with a site has a `:NonequalAI` group.
  - `"to_allosteric emits a :NonequalAI set only where its :OnlyA set is refused"`, on the
    six-step ping-pong of the worked case, in addition to the exact-children testsets: the
    children tagged `:NonequalAI` are exactly the pairs A with B, Q with P, B with P.
  - `"candidate_mechanisms is deterministic"`: two calls on one plain bi-bi parent return
    equal vectors.
- [ ] **Step 2: Run; they fail on the missing children.**
- [ ] **Step 3: Implement the move and the `_is_seed_node` change; delete the three
  helpers.** Rewrite both docstrings in the code's names, without "base", "K-type" or
  "V-type".
- [ ] **Step 4: Carry over the testsets that call the move directly** (inventory, section
  A, 33 rows) under the spec's two carry-over rules. Delete and list those that pin a
  result of `_tagging_bases` or the every-chemistry set (10193, 10313, 11028, 11119,
  11381 and the base assertions inside 3975). Move those about another aspect (oligomeric
  states 4763, regulator pairing 4875 and 7226, the parent gate 5202/5224/5259, fused
  chemistry groups 9871–10107, Theorell–Chance 9804) to a fixture of at most four groups
  or one oligomeric state. Counts such as "15 K-types" (3923), 31 (5284) and "69 = 6
  flips + 63" (9565) become `2^g − 1` plus the derived item 3 children; keep 9565's time
  budget and report the measured time.
- [ ] **Step 5: Add to `test_starting_mechanisms.jl`** the testset
  `"starting mechanisms hold one half-reaction :OnlyA and a :NonequalAI pair"`: with `R`
  required at one subunit, the six-step ping-pong with both chemistry steps at steady
  state, first half-reaction `:OnlyA`, `R::OnlyA`, is a starting mechanism, and so is the
  same mechanism with A and B `:NonequalAI`, `R::OnlyA`.
- [ ] **Step 6: Run `test_expansion_moves.jl` and `test_starting_mechanisms.jl`.** Re-pin
  `expand_bytes` (line 831, today `< 17_200_000`) at 3.5% above the measured value and
  report bytes before and after; the maintainer approved this re-pin. Re-pin
  `length(children) == 21` (786–832) and the seed counts. Regenerate the developer tables.
  Commit.

### Task 6: Populations built from the moves

**Files:**
- Test: `test/test_oracles.jl`, `test/Mechanisms/test_mechanism_type.jl`,
  `test/Enumeration/test_candidate_mechanism_type.jl`,
  `test/Derivation/test_rate_equation_correctness_guard.jl`,
  `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Run each file and collect the failures.**
- [ ] **Step 2: `test/test_oracles.jl`.** The testset at 1160 calls `_tagging_bases`,
  `_valid_tagging` and `_stops_inactive_turnover`: rewrite it against the children of
  `_expand_to_allosteric` that `_enumerator_keeps` keeps (admitted, distinct,
  oracle-checked), and re-pin its counts. Re-pin the counts at 1324, 1599 and 1622. An
  oracle check that now fails on a child is a finding: stop and report the mechanism.
- [ ] **Step 3: The other files.** Re-pin population counts (for example
  `test_mechanism_type.jl` 2975 and 3020). An assertion about a property, such as
  `test_candidate_mechanism_type.jl` 1719, that fails is a finding: stop and report.
- [ ] **Step 4: Add the testset `"a declared activator is never :OnlyI in a starting
  mechanism"`** to `test_starting_mechanisms.jl`: uni-uni reaction with `R::Activator`
  required; no starting mechanism has `R` tagged `:OnlyI`.
- [ ] **Step 5: Run all five files and `test_starting_mechanisms.jl`. Commit, with the
  list of re-pins in the message body.**

### Task 7: A competitive inhibitor's group can be `:OnlyI`

**Files:**
- Modify: `src/Enumeration/expansion_moves.jl` (`_expand_to_allosteric`,
  `_assert_emission_rules` docstring)
- Test: `test/Enumeration/test_expansion_moves.jl`,
  `test/Derivation/test_rate_equation_derivation.jl`

**Interfaces:**
- Consumes: `_bound_copy(c, k)` (names the inhibitor a binding step takes up, or
  `nothing`), `Mechanisms.steps_running_in_allosteric_conformations(m)` (fields `active`,
  `inactive`, one flag per step), `_redundant_copy_groups(c)`.
- Produces: `_inhibitors_bind(c::CandidateMechanism)::Bool` in `expansion_moves.jl`: true
  when every kinetic group that binds a competitive inhibitor has a step that runs in the
  active or the inactive conformation. Task 8 calls it.
- Produces: in item 1 of `_expand_to_allosteric`, a group that binds a competitive
  inhibitor takes `:EqualAI`, `:OnlyA` or `:OnlyI`, so `2^(g−k) · 3^k − 1` children for
  `k` such groups. Every child of the move (items 1 to 3) is dropped when `_inhibitors_bind`
  is false or `_redundant_copy_groups` is not empty. Item 2 skips a child with an `:OnlyI`
  group. In the search of item 3, and in the search of `_expand_change_allo_state`
  (Task 4 already tests `_redundant_copy_groups` there), a set whose child would be
  dropped this way does not gain and is extended.

- [ ] **Step 1: Write the failing tests**, exact children:
  - `"to_allosteric gives an inhibitor's group three states"`: uni-uni with `I` bound at
    free `E` (four groups), one subunit: `2^3 · 3 − 1 = 23` item 1 children less those
    dropped; among them the child with only the inhibitor's group `:OnlyI`.
  - `"to_allosteric drops an :OnlyI inhibitor that binds nowhere"`: a mechanism whose only
    inhibitor site is a form the `:OnlyA` set removes from the inactive conformation; no
    child has that `:OnlyA` set with the inhibitor `:OnlyI`, and the children with the
    inhibitor `:EqualAI` and `:OnlyA` are present.
  - `"to_allosteric drops a redundant :OnlyI reactant copy"`: derive a fixture from the
    spec's measurement (a reactant bound as an inhibitor); if none with at most five
    groups exists, report instead of weakening.
  - `"no child holds :OnlyI without a competitive inhibitor"`: the item 1 children of the
    uni-uni mechanism without an inhibitor hold no `:OnlyI` group.
  - `"change_allo_state tags an :OnlyI inhibitor :NonequalAI"`: one more fitted constant.
  - In `test/Derivation/`: one mechanism with an `:OnlyI` inhibitor, written out, compared
    with `_testhelper_reference_qssa` at random parameters, unless the file already has
    one (search for `:: OnlyI` on a step that binds an inhibitor).
- [ ] **Step 2: Run; they fail on the missing `:OnlyI` children.**
- [ ] **Step 3: Implement.** In the docstring of `_assert_emission_rules`, rewrite the
  sentence "No seed or move tags a catalytic group `:OnlyI`" and the sentence that names
  the moves that filter their children.
- [ ] **Step 4: Run `test_expansion_moves.jl`, the derivation file and
  `test/test_oracles.jl` (the testset at 1324 asserts "none tags catalysis :OnlyI": its
  assertion now reads "only an inhibitor's group is `:OnlyI`"; list the change). Commit.**

### Task 8: `_expand_add_dead_end_regulator` on an allosteric parent

**Files:**
- Modify: `src/Enumeration/expansion_moves.jl` (`_dead_end_child` 704–723 and its
  docstring, `_expand_add_dead_end_regulator` 782–787 and its docstring)
- Test: `test/Enumeration/test_expansion_moves.jl`

**Interfaces:**
- Consumes: `_inhibitors_bind(c)` (Task 7).
- Produces: `_dead_end_children(m, metabolites, forms, steps, reg)::Vector{CandidateMechanism}`
  replaces `_dead_end_child`. Plain `m`: the one child, or none when its inhibitor is
  redundant, as today. Allosteric `m`: the children with the new group `:EqualAI`,
  `:OnlyA` and `:OnlyI`, in that order, each kept when `_inhibitors_bind` holds and
  `_redundant_copy_groups` is empty.

- [ ] **Step 1: Write the failing tests**, exact children:
  - `"the dead-end move gives its inhibitor three states on an allosteric parent"`: the
    parent of the testset at 3637 (one placement): three children.
  - `"the dead-end move drops an :OnlyI inhibitor the inactive conformation cannot bind"`:
    a parent whose `:OnlyA` groups remove the placement's forms from the inactive
    conformation: the `:EqualAI` and `:OnlyA` children only.
  - `"an inhibitor added before or after to_allosteric gives the same mechanisms"`: one
    plain uni-uni mechanism, the placement at free `E`, the `:OnlyA` set {chemistry}: the
    three children of the dead-end move on that allosteric child equal the three item 1
    children of the mechanism with the inhibitor that have that `:OnlyA` set.
  - `"the dead-end move adds a second inhibitor beside an :OnlyI one"`: exact children.
- [ ] **Step 2: Run; they fail with one child where three are expected.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Give the testsets at 3582, 3637, 3788, 3857 and 5847 their new expected
  sets** (two more children per placement where the spec keeps them).
- [ ] **Step 5: Run `test_expansion_moves.jl`, `test_starting_mechanisms.jl`,
  `test/test_oracles.jl`, `test/Mechanisms/test_mechanism_type.jl`; re-pin populations and
  list them. Regenerate the developer tables. Commit.**

### Task 9: Reachability and measurements

**Files:**
- Test: `test/Enumeration/test_expansion_moves.jl`

- [ ] **Step 1: Write the aggregate testset** `"to_allosteric and change_allo_state reach
  every irreducible version"`, with a comment that it is an aggregate check. Helper
  `_testhelper_irreducible_versions(c)` tries every state for every group of a plain
  candidate at one subunit and returns the versions that are irreducible by the spec's
  definition (use `_enumerator_returns`, `_inhibitors_bind`, `_redundant_copy_groups` and
  `fitted_parameter_count`; it never calls the two moves). Helper
  `_testhelper_reached_versions(c, rxn)` closes `_expand_to_allosteric` and
  `_expand_change_allo_state` over their returned children. Assert, per parent:
  irreducible ⊆ reached; every child has more fitted constants than its parent;
  `_enumerator_keeps` on a version equals `_enumerator_keeps` with its `:NonequalAI`
  groups made `:EqualAI`. Parents: every seventh of the 278 base-tier mechanisms of the
  bi-bi reaction; every 150th of the 1,548 with `competitive_inhibitors: I`; five with a
  reactant bound as an inhibitor (`competitive_inhibitors: I, A`) of at most seven groups.
  Keep the testset under two minutes; thin the samples if needed and say how.
- [ ] **Step 2: Add the gain check on two regulators**: for five starting mechanisms of
  the reaction with `allosteric_regulators: R, S`, and their children under
  `_expand_merge_regulatory_sites` (the move that merges two regulatory sites), every
  child of `_expand_change_allo_state` has more fitted constants than its parent.
- [ ] **Step 3: Run. A version that is irreducible and not reached is a finding: stop and
  report the parent and the version.**
- [ ] **Step 4: Add a derivation test** to `test/Derivation/test_rate_equation_derivation.jl`:
  the six-step ping-pong with A and B `:NonequalAI` and both half-reactions `:EqualAI`,
  written out, compared with `_testhelper_reference_qssa` at random parameters and
  concentrations (`.claude/CLAUDE.md`, "Derivation tests").
- [ ] **Step 5: Measure and report** (not tests): the number of starting mechanisms for
  the bi-bi reaction with `R` required at `(1, 2)`, with `I` required, and with both; the
  time of `candidate_mechanisms(rxn, parent)` on ten ter-ter allosteric parents; and a
  recovery check with `identify_rate_equation` on data generated from the six-step
  ping-pong with A and B `:NonequalAI`.
- [ ] **Step 6: Commit.**

### Task 10: Docs, roadmap and version

**Files:**
- Modify: `docs/src/identify/enumeration_engine.md`, `docs/src/deriving/re_vs_ss.md`,
  `docs/src/developer.md`, `docs/src/deriving/mwc_allostery.md`, `ROADMAP.md`,
  `docs/src/roadmap.md`, `Project.toml`

- [ ] **Step 1: Engine page.** Rewrite moves 4 and 6 and the passages the spec lists under
  "Guards and docs" (lines 198–212, 303–372, 389–390, 439–442 and "Starting from the
  required regulators"), in the code's names and for a reader who has not read the code.
- [ ] **Step 2: The other pages** named in the spec's "Guards and docs" (developer page
  697, 745–746, 781–782; MWC page 150; `re_vs_ss.md`).
- [ ] **Step 3: `ROADMAP.md`.** Apply every bullet of the spec's two `ROADMAP.md` items,
  with the numbers Task 9 measured. Write the new Known issues for items 1, 2, 4, 5 and 6
  of "Known costs" and the new Idea.
- [ ] **Step 4: `Project.toml`** version `0.11.1`.
- [ ] **Step 5: Build the docs** (`julia --project=docs docs/make.jl`) and run
  `test/test_boundaries.jl`. Commit.

### Task 11: Whole suite and cleanup

- [ ] **Step 1: Run the full suite** in the background:
  `timeout 1200 julia --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`.
- [ ] **Step 2: Run** `test/test_type_stability.jl` under `julia +1.13`.
- [ ] **Step 3: Resolve failures** under the test-change rule. A changed selection in
  `test/test_identify_rate_equation.jl` is a finding for the maintainer.
- [ ] **Step 4: Re-read the changed source files** for dead code, stale comments and
  helpers the change made unnecessary.
- [ ] **Step 5: `git rm` the spec and this plan. Commit.**
