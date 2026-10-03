# Merged and Theorell–Chance seeds, fused steps in the moves (sub-project C)

Date: 2026-10-03. Status: design agreed with Denis section by section on 2026-10-02 and
2026-10-03. Part of the A + B + C plan in `2026-09-28-identifiable-enumeration-findings-and-plan.md`
(plan section "C", points 3–7, 9, 10 and 12); the handoff `2026-10-02-c-handoff.md` gives the state
after A and B. Built on branch `step-explicit-stoichiometry` after B, for the one squashed pull
request that carries A, B and C. No version bump (the branch is at 0.8.0, one over `main`).

## Why

The enumerator still misses textbook families at their true parameter counts, and it still emits
phantoms of the chain class.

- **Missing families.** Merging a seed's central complex and making two of its binding groups
  steady state gives 114 new bi-bi families at 5 parameters, the seed count, none of them reached
  today at any level (findings point 9, EXACT). Among them are the full steady-state ordered
  forward law V·A·B/(Kia·Kb + Kb·A + Ka·B + A·B) and the classic ping-pong law
  V·A·B/(Kb·A + Ka·B + A·B). Theorell–Chance adds 6 more families at 5, and random-order seeds
  first become clean at 6 parameters.
- **Degenerate ping-pong seeds.** All 7 of today's ping-pong seeds are degenerate: at zero
  products their rate does not depend on B (PROVED). They are also the only route to 83
  identifiable ping-pong mechanisms at levels 1–2 (measured below), so the beam must keep them as
  parents without fitting them as models.
- **Chain phantoms.** Flipping a flank of a qualifying chain (RE binding, SS isomerization, RE
  release) adds one or two parameters and no rate law (chain lemma, findings points 3–5). On bi-bi
  at depth ≤ 2 such mechanisms carry 154 of the 266 phantoms left after B.
- **Fused steps.** The moves recognize chemistry only as an isomerization and reject any other
  parent (`_assert_chemistry_is_iso`). Merged and Theorell–Chance seeds need the moves to accept
  them.

## Goals

1. `init_mechanisms` emits, beside today's seeds, every merged and Theorell–Chance variant that is
   inclusion-minimal, valid, flux-carrying and non-degenerate (section 3).
2. The beam expands degenerate seeds without fitting them (section 4).
3. The flip never flips a flank of a qualifying chain whose isomerization is steady state, in a
   plain `Mechanism` (section 2).
4. Every move accepts and handles fused and Theorell–Chance steps (sections 1 and 2).
5. B's code is simplified first, with identical output except where the split pre-filter's
   deletion adds rule-conforming children (section 0).

## Non-goals

- **A canonicalization pass.** No move output is rewritten to a chain's class representative:
  no merge of rapid-equilibrium-chemistry chains (classes C, D and RSTOR), no unmerge of a merged
  complex whose two steps are steady state, no rewrite of a chain that a split makes private. The
  flip rule (goal 3) gets 105 of the 154 chain-related phantoms at depth ≤ 2 for a fraction of the
  code; the remaining non-canonical mechanisms (124 on bi-bi to depth 3) stay as tolerated phantoms
  and duplicates.
- **Chain rules for allosteric mechanisms.** At catalytic multiplicity 2 a chain flank's flip
  raises the rank in 80 of 248 cases (measured below), so the flip rule applies to `Mechanism`
  only.
- **Duplicate families among the new seeds.** 8 ordered/random pairs and 2 Theorell–Chance pairs
  have equal families but different structures; both members stay (10 extra fits per bi-bi run).
- **Expand-only beyond seeds.** A degenerate mechanism below the seeds is fitted and selected as
  today; 4 identifiable level-2 ping-pong mechanisms have only degenerate level-1 parents.
- **D and E** of the findings (MERGE and ELIM as moves anywhere; rank on the beam's axis).
- **The `_block_balanced` rewrite** proposed by the review of B.
- **Uni-uni with S required as an inhibitor** still errors (B's behaviour): its only clean merged
  variant is the lumping twin of the seed.
- **Allosteric variants of a merged mechanism** tag its fused binding together with its chemistry;
  they cannot tag that binding on its own.

## Decisions (Denis, 2026-10-02 and 2026-10-03)

- **A binding is defined by its stoichiometry.** A step that takes up exactly one metabolite and
  gives off none binds it, whether or not the enzyme's composition changes beyond that metabolite.
  The merged complex is the product-bound form.
- **Flip rule only**, no canonicalization pass. Denis first chose RSR as the representative of
  the FULL chain class; with no pass to apply it, that choice survives as the seeds' lumping-twin
  skip.
- **Degenerate seeds are expanded without being fitted**, seeds only.
- **Keep both members of duplicate families** among the new seeds; the decorated Theorell–Chance
  seeds (6 or 7 parameters) are in.
- **Simplify B first**: every item of the review except the `_block_balanced` rewrite. The unit
  tests of the deleted twin helpers are rewritten, not dropped; the split pre-filter is deleted.
- **Copy gauge: one shared scale factor** on all of a copy's complexes, never their twins'
  factors (confirms B's controller ruling; handoff item 2).
- **CLAUDE.md enumeration-test rule 3** is reworded (section 6).
- The branch-finishing menu and E remain open (handoff items 1 and 3).

## Terms

- **Plain binding**: a binding whose `to_species` is its `from_species` with the metabolite added
  (`_binds_ligand`). **Fused binding**: a binding that is not plain, such as EA + B → E(P, Q).
- **Chemistry step** (`_is_chemistry`): every step that is not a plain binding: isomerizations,
  fused bindings and Theorell–Chance steps.
- **Theorell–Chance step**: a step that takes up one metabolite and gives off another,
  EA + B → EQ + P.
- **MERGE** of an isomerization X1 → X2 that is alone in its group: delete it and its group, and
  move every other step at X1 onto X2 (the product side). Steps keep their metabolites, flags and
  groups.
- **ELIM** of a form X with exactly two steps, Y + L → X and Z + R → X: replace both by one step
  Y + L → Z + R in a new group; the two old groups keep their other steps.
- **Qualifying chain**: an isomerization X1 → X2 alone in its group, where X1 and X2 each have
  exactly one other step, each a binding into that form (the form is the binding's `to_species`)
  and each alone in its group. Its **flanks** are those two bindings; its **type** lists the three
  flags, R for rapid equilibrium and S for steady state, in the order entry, isomerization, exit.
- **Lumping twin**: a merged variant whose every merged complex has both its steps steady state,
  each alone in its group. It has the family of the unmerged RSR form at the same count (lumping
  theorem, track 4 §5.1).
- **Degenerate**: the rate at zero products does not depend on some substrate, or has no maximum
  as the substrates grow; or the mirror statement for products (track 4, flags e and f).

## Design

### 0. Opening wave: simplify B

The items of `2026-10-03-b-simplification-review.md` (the read-only review of B, line numbers at
b340822), each a commit of its own, before any other change:

- **Docstrings.** Each argument gets one home: "a complex that duplicates only a copy's complex
  is no twin" (keep `_productive_twin`), "all-twin with a failing gauge is kept" (keep
  `_redundant_copy_groups`), "a placement never makes an older copy redundant" (keep
  `_dead_end_child`), Theorem 2's finite merge (keep `_gauge_rescaling`), the twin definition
  (keep `_productive_twin`), "a zero-flux group's two constants enter only as their ratio" (keep
  `_flux_carrying_groups`), per-site competition (keep `_expand_add_dead_end_regulator`). Remove
  measurements and spec case labels from docstrings. Fix the stale `expand_mechanisms` text, which
  still states the withdrawn "creates a new complex" rule. Handoff cleanups 3 (`_gauge_rescaling`
  calls the finite ρ the shared factor) and 4 (a docstring names `_dead_end_child` where a
  `Mechanism` takes the `_gauge_rescaling` path). Shorten the `seed_mechanisms` error and docstring,
  with the test phrase that pins it.
- **Pure refactors**, identical output:
  - one `_expand_add_dead_end_regulator` method for both mechanism types, absorbing
    `_expand_add_dead_end_regulator_native`;
  - one helper for the RE-segment index preamble that `_flux_carrying_steps`, `_productive_twin`
    and `_hyperbolic_catalysis` repeat (in the enumeration file; `_re_segment_extras` itself is
    untouched, since every `Mechanism` construction calls it);
  - `_with_equilibrium(s, flag)` for the four step rebuilds, making `_flip_group_to_ss` and
    `_all_steady_state` one-liners;
  - `_copy_complex` and `_composition` inlined;
  - `_split_gain_test` as one loop over both parts;
  - the two dead clauses of `_redundant_copy_groups(am)` removed, and its inactive view built from
    `_reachable_from_free` instead of a `Mechanism` round trip;
  - `_productive_twin(steps(m))` built once per dead-end parent (handoff cleanup 1).
- **One twin predicate.** Keep `_all_twin(part, twin)`; `_gauge_rescaling` checks it first; delete
  `_copy_twin_test` and `_twin_only`. The allosteric dead-end child runs the active-state twin test
  as a guard and `_redundant_copy_groups` on the built child, which gives the inactive-state verdict
  of the deleted function in every branch (the review's item 1 argues each). The unit tests that
  call the deleted helpers (`test/test_mechanism_enumeration.jl` 7194–7221, 7275, 7318) are
  rewritten onto `_redundant_copy_groups` and `_all_twin`, with one new `:OnlyA`-copy fixture and
  a rank pin. `docs/src/developer.md:104` stops naming `_copy_twin_test`.
- **Pin `ss_mirror`'s rank** at 8 fitted, rank 6 (handoff cleanup 2).
- **Delete the split's unit pre-filter** (`redundant_part`). B's spec §4 claims that splitting
  other groups only relaxes a copy part's gauge; it can break it instead. In H1 with its two A
  groups merged into one, both twins share the A group's class and the copy is redundant; splitting
  that A group gives the twins two classes and the shared B group reads two ratios, so the gauge
  fails. Without the pre-filter the split can only add children, each satisfying the copy rule,
  and only when the copy-group unit fails its gain test alone. The child post-filter stays. B's
  spec §4 gets the corrected sentence.
- Kept as they are: the flip's all-steady-state unit filter (proved and essential), the split's
  bottomless guard, and the `Mechanism` pre-construction gauge in the dead-end move (its removal
  would take R6 enumeration from +18% to about +23% over the commit before B).

B's regression harness reruns on the opening wave alone. Every population is identical except for
children the pre-filter deletion adds; the record counts and ranks each.

### 1. Steps are classified by stoichiometry

- **`bound_metabolite(s)`** returns M when `consumed(s) == [M]` and `released(s)` is empty, and
  `nothing` otherwise; `is_binding` follows it. The composition check (`_binds_ligand`) leaves
  both.
- **The `Step` constructor** stores a step that takes up nothing and gives off exactly one
  metabolite as the binding it reverses, whatever its forms' compositions; today it does so only
  for a plain release.
- **`_is_chemistry(s)`** is true for every step that is not a plain binding:
  `!(is_binding(s) && _binds_ligand(from_species(s), to_species(s), bound_metabolite(s)))`. Only
  the allosteric moves and the `:OnlyA` check read it (section 2).
- **Kinetic groups.** `_step_kind(s)` becomes the pair `(consumed(s), released(s))`, the key
  `_apply_equivalence_grouping` already groups by. A group holds steps of one kind and one flag,
  so EA + B → E(P, Q) and EQ + B ⇌ EBQ can share the B group.
- **Derivation.** Every `is_binding` test in the derivation now admits fused bindings: the RE
  weight [M]/K (`_re_weight_ratio`), the parameter types (`Kd`, `Kon`, `Koff`), the binding-K sets
  of `_assemble_constraints` and `_build_wegscheider_rename_map`, and `_count_kind`. A fused
  binding's RE constant is therefore a dissociation constant named in the release direction
  (`K_EPQ_to_EA_B`), and a group holding a fused and a plain binding shares one consistent
  constant. Nothing else in the derivation changes.
- **Effect on existing mechanisms.** None for a mechanism without fused steps: rate strings and
  names stay identical. In a hand-written fused mechanism an RE fused binding's constant is
  inverted and renamed, and a fused release is stored reversed; its rate law is unchanged.
- **The merged complex** is the product side of the merged isomerization. A merge therefore moves
  every step at the substrate-side form onto it: the last substrate's binding becomes fused, and
  the releases stay plain.

### 2. The moves on fused steps, and the flip rule

- **Removed:** `_assert_chemistry_is_iso` and its call in `expand_mechanisms`; the "binding or
  isomerization" check of `_assert_mechanism_invariants(::Mechanism)`.
- **The allosteric moves read `_is_chemistry` instead of `is_iso`:** the chemistry and binding
  groups of `_expand_to_allosteric` (K-type sets every chemistry group `:OnlyA` and ranges its
  subsets over the other groups; V-type sets the chemistry `:OnlyA`), `_partial_onlya_catalysis`,
  the chemistry groups that `_expand_change_allo_state` relaxes together, and
  `_onlya_haldane_violation`, which drops `:OnlyA` chemistry groups from its cycle graph and treats
  only `:OnlyA` plain bindings as `:OnlyA` bindings. A hand-written ping-pong whose fused release
  is tagged `:OnlyA` is then accepted.
- **Dead-end sites.** Competition with M targets every form where M is free on some step of the
  role's kind: `from_species` when M is consumed, `to_species` when M is released.
  `_forms_with_binding_step_native` becomes this rule. For plain and fused bindings it gives
  today's sites, since a binding is stored with M consumed; a Theorell–Chance step adds its own
  (EA for B and EQ for P in EA + B → EQ + P, without which the ordered Theorell–Chance seed has no
  site for B or P). The capacity test and the exclusion of forms that already hold M read
  compositions, unchanged.
- **The split is unchanged.** A fused binding's context form is the form its metabolite binds to
  (`_context_form`).
- **The flip rule.** In a `Mechanism`, a group is not a flip unit when it is a flank of a
  qualifying chain whose isomerization is steady state. Such a flank's flip always gains alone (it
  cuts the chain form off its segment) and carries flux (the isomerization, its only neighbour, is
  a private steady-state group, so rule 1 makes it carry flux), so the minimal-set search never
  extends it, and excluding the unit removes exactly the single-flank children. The chain lemma
  makes each removed child's family that of its parent. An `AllostericMechanism` keeps flipping
  flanks.
- **Unchanged:** the split, the flip's flux and bottomless tests, the copy rule,
  `_hyperbolic_catalysis` and `_assert_emission_rules`.

### 3. The new seeds

`init_mechanisms` keeps every seed it builds today and adds, for each one, its merged and
Theorell–Chance variants, all as plain `Mechanism`s:

1. **Merged base.** MERGE every chemistry isomerization of the seed: one for a sequential seed,
   two for a ping-pong seed. Every step of the base is rapid equilibrium, so its rate is infinite.
2. **Theorell–Chance bases.** In the merged base, ELIM each merged complex that has exactly two
   steps, one complex at a time. The new step is rapid equilibrium, like the rest of the base.
3. **Search.** For each base, the units are its groups that carry flux in its all-steady-state
   copy (as in the flip). `_minimal_gaining_sets` finds every inclusion-minimal set of units whose
   flip to steady state gives a candidate that is valid, carries flux and is non-degenerate. Each
   such set is a variant.
4. **Skip lumping twins** among the merged variants. In bi-bi these are the merged {B, P} variants
   of the 8 seeds whose two groups next to the central complex are private.
5. **Deduplicate** today's seeds and the variants structurally (`unique!`).

The tests on a candidate run on its groups, before any `Mechanism` is built:

- **Valid:** no all-RE turnover cycle, and no bottomless RE segment (`_bottomless_re_segment`, the
  constructor's rule). `_re_turnover_cycle(groups, rxn)` tests the first: give each RE step the
  weight Σ_x ρ_x·σ[x] of section 1 of B's spec (ρ = +1 for a substrate, −1 for a product), and an
  RE component holds a turnover cycle iff it is unbalanced (`_block_balanced`).
- **Carries flux:** every steady-state group carries flux (`_flux_carrying_groups`).
- **Non-degenerate**, track 4 §7's conditions V and C. With the flux rule in force, V and C alone
  gave the same verdicts as V, C and I on about 5,600 flux-valid mechanisms. V's known false
  rejections need copies of substrates or products, which seeds never carry.
  - **V (a maximal rate exists both ways).** Turn rapid equilibrium every steady-state step that
    takes up or gives off a substrate; `_re_turnover_cycle` must stay false. The same holds with
    products.
  - **C (chemistry does not equilibrate both sides).** A chemistry node is one of: a merged
    complex (a form entered by a fused binding of a substrate) with every form joined to it by RE
    isomerizations; a set of forms joined to each other by RE isomerizations; an RE
    Theorell–Chance step. An RE step t leaves a node through a node form F and releases M when
    traversing t away from F gives M off: F is `to_species(t)` and M is consumed, or F is
    `from_species(t)` and M is released. An RE Theorell–Chance node leaves through its own step in
    both directions. C fails when one node is left both by an RE step that releases a substrate
    and by an RE step that releases a product. Examples:
    - merged {A, Q} fails: E(P, Q) releases B through its RE entry and P through its RE exit;
    - merged {A, P} passes: P's release is steady state;
    - every ping-pong seed fails: the node {E(B; res), EQ} of the RE isomerization releases B and Q
      by RE steps;
    - the SRR ping-pong half passes: its B binding is steady state;
    - an RE Theorell–Chance step fails by itself.

    A merge can also turn a dead-end binding at the substrate-side form into a fused binding of a
    product (a decorated ping-pong seed's EA + P becomes E(P; res) + P → E(A, P)); the node
    definition counts only fused bindings of substrates, so that form is no chemistry node.
- **`_degenerate(m)`** (`!(V && C)` on `steps(m)`, the active-state graph for an allosteric
  mechanism) serves section 4.

Output and counts:

- Every reaction gets these variants. An undecorated ordered topology with n substrates and m
  products gains n + m − 2 merged seeds (findings point 9, uni-uni to ter-ter) plus its
  Theorell–Chance seeds.
- `init_mechanisms` returns mixed parameter counts. The beam already fits every base mechanism in
  its base tier and expands each when its sweep reaches that count.
- Expected for bi-bi, from the findings; the regression record confirms or corrects it:

  | seeds | count |
  |---|---|
  | today's seeds (the 7 ping-pong ones expanded without fitting) | 62 |
  | two-group merged, clean (114 families plus the 8 ordered/random pairs) | 122 |
  | Theorell–Chance: 8 ordered at 5, 12 decorated ordered at 6–7, 4 ping-pong halves at 5 | 24 |
  | three-group merged at 6: 28 random/random, 28 ordered/random | 56 |
  | total | about 264 |

- `seed_mechanisms` builds on `init_mechanisms`, so required regulators get the new seeds. Uni-uni
  gains none (its only clean merged variant is the lumping twin), so uni-uni with S required as an
  inhibitor still errors.

### 4. The beam expands degenerate seeds without fitting them

- `_beam_search` splits its base set (from `init_mechanisms`, or `seed_mechanisms` when regulators
  are required) with `_degenerate`. Non-degenerate seeds are fitted as today. Degenerate seeds are
  expanded at once (`expand_mechanisms`), and their children join the base tier, fitted beside the
  other seeds.
- `init_mechanisms` and `seed_mechanisms` return the same sets as without the split.
  Expand-only seeds spend none of the per-count expansion budget. Their children's rows carry no
  parent, as seeds' rows do.
- In bi-bi the expand-only seeds are the 7 ping-pong seeds. Their copy and allosteric children
  stay degenerate: a copy or a conformation adds no chemistry, so C still fails and B is still not
  needed.

### 5. Measurements behind the decisions

All on the branch at b340822 (after B), with the rank oracle `_testhelper_identifiable_rank` and a
numeric degeneracy probe (rate ratios at 10⁻⁸ against 1 for each substrate and product, and at
10⁹ against 10⁶ for the maximal rates). The probe flags all 7 ping-pong seeds and reproduces the
track's 83 identifiable descendants. Scripts in the session scratchpad: `deadend_chain_check.jl`,
`pingpong_seed_check.jl`, `allo_chain_check2.jl`, `canon_worth.jl` and its log, `pp_paths.jl`.

- **The dead-end move never targets a chain form.** R6 (bi-bi with A, B, P, Q also inhibitors),
  levels 0–1: 1,711 parents, 234 with a qualifying chain (70 ping-pong), 44,312 dead-end children,
  0 adding a step at a chain form. R5 (one foreign inhibitor), levels 0–2: 4,527 parents, 587 with
  a chain (187 ping-pong), 9,963 children, 0. The move does target ping-pong isomerization ends that
  carry an abortive binding (4,084 and 1,000 placements); those ends have three steps and form no
  qualifying chain.
- **The degenerate ping-pong seeds' descendants** (flip and split, depth 2): 202, of which 128 are
  non-degenerate and 83 identifiable. Making both chemistry steps steady state in the seeds instead
  (Denis's alternative) reaches only 32 of the 83, even one level deeper. Expanding the seeds
  without fitting reaches all 83; 4 of the 67 at level 2 have only degenerate level-1 parents
  (level 1: 47 mechanisms, 26 degenerate, 16 identifiable; level 2: 155, 48 degenerate, 67
  identifiable).
- **Chain flank flips in allosteric variants** of the four undecorated ordered seeds: plain
  `Mechanism` 0 of 8 raise the rank; K-type at multiplicity 1, 0 of 120 (L itself is a phantom
  there); K-type and `:NonequalAI`-relaxed flanks at multiplicity 2, 80 of 248 (16 of 64 with a
  `:NonequalAI` flank). Every multiplicity-2 parent is identifiable.
- **The flip rule** (bi-bi R4, flip and split, levels 0–3: 62, 369, 1,200, 2,795 mechanisms; with
  the rule 62, 349, 1,134, 2,677). Mechanisms with a non-canonical qualifying chain: 2, 28, 102,
  196. The rule removes 204 and adds none, and every removed mechanism keeps its RSR twin in the
  population. Phantoms at levels 0–2: 266 in all, 154 in mechanisms with a non-canonical chain, 105
  in mechanisms the rule removes. Level 3, a sample of 150 non-canonical mechanisms: 185 phantoms,
  119 of them in the 88 the rule removes. 124 non-canonical mechanisms remain to depth 3 (chains
  made private by a split, and the ping-pong rapid-equilibrium-chemistry chains).

## Test changes

Every change is made test-first. New tests follow CLAUDE.md's enumeration-test rules: fixtures
inline in the testset, exact child sets and counts, every expected child hand-derived in a comment
and checked with the rank oracle, `_testhelper_` prefix for file-level helpers.

- **Step kinds:** `bound_metabolite` of a fused binding (B), of a fused release (stored reversed)
  and of a Theorell–Chance step (`nothing`); `_is_chemistry` on each kind; a group holding a fused
  and a plain binding of B is accepted; a new mass-action oracle fixture with such a mixed group
  (the merged decorated ordered seed, B group {EA + B → E(P, Q), EQ + B ⇌ EBQ}); the renamed fused
  fixtures in `test/test_rate_eq_derivation.jl`.
- **Moves:** exact allosteric children of a merged and of a Theorell–Chance mechanism; the
  `:OnlyA` fused ping-pong accepted by `_onlya_haldane_violation`; exact dead-end children of the
  ordered Theorell–Chance seed (its B and P sites); the split of a mixed group; the flip of a merged
  mechanism.
- **Flip rule:** the undecorated ordered seed flips no flank; a decorated seed whose flank group is
  shared still flips it; an allosteric variant at multiplicity 2 still flips its flank.
- **Seeds:** unit tests of MERGE, ELIM, `_re_turnover_cycle`, V and C on hand cases (merged {A, Q},
  {A, P}, a ping-pong seed, an RE Theorell–Chance step); exact variant sets, each child with
  fitted = rank and non-degenerate by a `_testhelper_degenerate` probe, for uni-uni (none), the
  undecorated ordered seed ({A, Pˣ}, {Bˣ, Q} and its Theorell–Chance seed; twin skipped), one
  decorated ordered seed, one ordered/random seed, and the undecorated ping-pong seed (PP-AQ,
  PP-BP and two half-Theorell–Chance seeds).
- **Pins**, each stated as a counterfactual: `init_mechanisms` counts for bi-bi and ter-ter; the
  depth-2 population pins (today 62/369/1,200 and 62/1,649/31,730).
- **Beam:** a degenerate seed is absent from the initial CSV and its children are present.
- **Rewritten, not deleted:** the twin-helper unit tests (section 0); every exact-children test
  whose parent has a chain flank the flip now skips; tests that assert `expand_mechanisms`
  rejects fused parents become tests that it expands them.

## Verification

1. **Regression record**, B's method: a detached worktree at the opening wave's tip against C's
   tip, keyed by steps with explicit lists.
   - Rate strings and names identical for every mechanism without fused steps
     (`MECHANISM_TEST_SPECS`, R1–R6 to depth 2); fused fixtures give equal laws under the name map.
   - `init_mechanisms` counts per reaction (R1–R6, uni-bi, ter-ter), bi-bi expected near 264. Every
     new seed derives, passes the mass-action oracle on a sample, is non-degenerate by the probe,
     and is reported with fitted and rank.
   - V and C against the probe on every candidate the seed search tests. A single disagreement
     stops the work for investigation.
   - The flip rule: every removed mechanism has its RSR twin present; phantom counts before and
     after.
   - Per-level counts for R1–R6 to depth 2 (R6 sampled), the level-3 phantom fraction, and the
     allosteric sample.
2. **Timing:** enumeration time per emitted mechanism on R6 levels 0–2 within +20% of the opening
   wave's tip (the total grows with the seeds); `init_mechanisms` and `seed_mechanisms` times on
   ter-ter reported, and ter-ter seeding over one minute is brought to Denis before going on; the
   ter-ter worst-case split test and the `rate_equation` performance test stay green.
3. **Full suite** before every commit.

## Risks

- **Mixed groups in the derivation** are new. The mass-action oracle fixture and the regression
  sample of new seeds guard them.
- **Fitting cost.** The bi-bi base tier grows from 62 to about 264 fits, by design.
- **Ter-ter seeding cost** is unmeasured; the timing check bounds it.
- **V and C** are checked against a numeric probe, not proved for every topology; the regression
  record compares them on every candidate.
- **Hand-derived expectations.** Each is derived in its test comment and checked with the rank
  oracle; five of B's plan expectations were wrong, and derivation caught them.

## Documentation

- `docs/src/identify/enumeration_engine.md`: the seed construction (merged and Theorell–Chance
  variants, the non-degeneracy conditions, skipped lumping twins), the flip rule, the paragraphs
  that call chemistry the isomerization step (around lines 36–42 and 313–319) rewritten around
  `_is_chemistry`, dead-end sites on Theorell–Chance steps, and the opening wave's changes to B's
  sections.
- `docs/src/deriving/kinetic_groups.md:22`: a group's steps share their consumed and released
  metabolites and one flag.
- `docs/src/deriving/re_vs_ss.md`: an RE label on a flank of a steady-state chain is an apparent
  constant, a Michaelis-type combination of the chain's rate constants, not a dissociation
  constant; it gives the same rate laws as the all-SS chain with fewer parameters.
- `docs/src/identify/model_selection.md`: degenerate seeds are expanded without being fitted.
- `docs/src/identify/combinatorics.md`: seed counts, wherever it quotes them.
- `docs/src/developer.md`: stoichiometric step kinds, `_is_chemistry`, seed construction.
- Docstrings of `init_mechanisms` (it no longer returns only the minimum count) and
  `expand_mechanisms`.
- The findings document's section C points here and records the decisions; B's spec §4 gets its
  corrected sentence.
- CLAUDE.md, enumeration-test rule 3, approved wording:

  > 3. **Write chemistry the enumerator's way.** Unmerged chemistry is an isomerization to the
  > product-bound form followed by release steps (`E(A) <--> E(P; residual = A - P)` then
  > `E(P; residual = A - P) ⇌ E(; residual = A - P) + P`). A merged complex is the product-bound
  > form: the last substrate binds into it in one step (`E(A) + B <--> E(P, Q)`) and the products
  > leave by plain releases. A Theorell–Chance step takes up the substrate and gives off the
  > product in one step (`E(A) + B <--> E(Q) + P`). Never write a merged complex on the substrate
  > side (`E(A, B) <--> E(Q) + P`). The enumerator never emits that form, so a fixture written that
  > way is a different mechanism from the one the move produces.
