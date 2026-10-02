# Exact filters: regression check

Date: 2026-10-02. This checks the exact enumeration filters (design
`2026-09-30-exact-filters-design.md`, plan `../plans/2026-09-30-exact-filters.md`): the flux
predicate, the flip's extend semantics, the split's revert, the copy rule in its one reading
(section 2 of the design, "One reading") and the dead-end move's per-site competition
(section 5). Every mechanism the filters remove must be one a rule names; every kept and added
mechanism must satisfy both emission rules; every kept mechanism must derive as before; added
mechanisms are reported with their phantom fraction.

## Method

A scratch script (not committed) enumerates the populations of one checkout and writes them as
keys; a second mode compares two key sets mechanism by mechanism. The baseline ran on a2a02b1,
before the filters; the comparison ran on d0706a7, the branch commit that holds them. Both
enumerated independently, one Julia process at a time; the a2a02b1 run is the snapshot of the
earlier comparisons, made with the current `Manifest.toml`.

The copy rule has one reading. A copy complex is a twin only of a productive complex, a form
bound to no competitive inhibitor, judged in every conformational state where the copy binds
(`_copy_twin_test`). Every conformation holds its free enzyme, a form with no bound metabolite
and no residual (the forms `_reachable_from_free` starts from): a free enzyme that the inactive
step graph lacks has no complex there to duplicate, so a copy bound to it is new in that state.

Competition in the dead-end move is decided per site. A substrate or product declared as a
competitive inhibitor binds as a copy at a dead-end site of its own, so the move targets, for
competition with a reactant, the forms where that reactant binds productively, and for
competition with an inhibitor the forms where that inhibitor binds; it excludes a form that
holds a competing reactant productively or a competing inhibitor, and its capacity test counts
productive bindings only. a2a02b1 decided competition by name, a reactant and its copy counting
as one competitor.

Populations: `init_mechanisms` plus two levels of the flip (`_expand_re_to_ss`), split
(`_expand_split_kinetic_group`) and dead-end (`_expand_add_dead_end_regulator`) moves on
reactions R1–R6 of the findings document, deduplicated across levels; ALLO, the allosteric
children of the R4 seeds (`_expand_to_allosteric`) with their flip and split children; and
ALLO6, which carries copies into allosteric mechanisms. R1 is uni-uni, R2 uni-uni with one
dead-end inhibitor, R3 uni-uni with S and P as inhibitors, R4 the bi-bi reaction whose atoms
admit ping-pong, R5 R4 with one dead-end inhibitor, R6 R4 with A, B, P and Q as inhibitors.
ALLO6 level 1 is 150 of the allosteric children of the R6 seeds, drawn with
`MersenneTwister(20260930)` from the children sorted by key; level 2 is their dead-end,
tag-relaxation (`_expand_change_allo_state`), flip and split children; level 3 is the tag
relaxations of the level-2 mechanisms that bind a copy, with each one's level-2 parents
recorded. A mechanism is keyed by the sig strings of its reaction and steps, group by group,
with an allosteric mechanism's states and multiplicity appended.

For each removed mechanism the script records the first rule that names it: a steady-state
group that carries no flux (`_flux_carrying_groups`); else a copy group whose every site is a
twin (`_duplicate_copy_groups`); else a mechanism the new code reaches by placing a copy on its
parent without that copy, the parent being in the new population (other order); else
UNEXPLAINED. For a zero-flux removal it checks that the mechanism's rapid-equilibrium twin is in
the new population, and for the removed R4 zero-flux mechanisms that the numerical rank of
∂v/∂log θ (finite differences, h = 1e-5, 60 points, three draws, singular values above 1e-7
of the largest) is at least 1 and fitted parameters minus rank at least the number of
zero-flux groups. It runs `_assert_emission_rules` on every kept and added mechanism and records
each population's enumeration time (one run per checkout).

For the added mechanisms it counts them per set and level, checks both rules on each, and draws
up to 150 per set and level (sorted keys, `MersenneTwister(20260930)`). Each drawn mechanism is
ranked against every parent the new population holds, the mechanism without one of its copies,
and the script checks that the dead-end move reaches it from such a parent. A mechanism adds a
phantom when its fitted count minus rank exceeds its parent's.

It also compares the new keys with the previous snapshot (87f49a7, which decided competition by
name under the one reading of the copy rule).

Level-3 sample (`level3_sample.jl`, new checkout only): on R6, 300 level-2 parents drawn with
`MersenneTwister(20260930)`, their split and dead-end children, the children with a copy group
whose every site duplicates some form when every form counts as a twin source but not when only
productive forms do, and the rank oracle on 150 of them (a seeded draw) against their parents.

The derivation comparison draws, per set, up to 40 keys (all of R1–R3) from the sorted
intersection of the two key sets with `MersenneTwister(20260930)`, derives each on both
checkouts (`fitted_params` and the Reduced `rate_equation_string`) and compares them. a2a02b1
ran as a detached worktree with the current `Manifest.toml`, removed afterwards.

## Results

| Set | Level | Old | New | Kept | Removed | Added |
|---|---|---|---|---|---|---|
| R1 | 0 | 1 | 1 | 1 | 0 | 0 |
| R1 | 1 | 2 | 2 | 2 | 0 | 0 |
| R1 | 2 | 1 | 1 | 1 | 0 | 0 |
| R2 | 0 | 1 | 1 | 1 | 0 | 0 |
| R2 | 1 | 3 | 3 | 3 | 0 | 0 |
| R2 | 2 | 3 | 3 | 3 | 0 | 0 |
| R3 | 0 | 1 | 1 | 1 | 0 | 0 |
| R3 | 1 | 4 | 2 | 2 | 2 | 0 |
| R3 | 2 | 6 | 1 | 1 | 5 | 0 |
| R4 | 0 | 62 | 62 | 62 | 0 | 0 |
| R4 | 1 | 369 | 369 | 369 | 0 | 0 |
| R4 | 2 | 1388 | 1200 | 1200 | 188 | 0 |
| R5 | 0 | 62 | 62 | 62 | 0 | 0 |
| R5 | 1 | 719 | 719 | 719 | 0 | 0 |
| R5 | 2 | 3934 | 3746 | 3746 | 188 | 0 |
| R6 | 0 | 62 | 62 | 62 | 0 | 0 |
| R6 | 1 | 1769 | 1409 | 1409 | 360 | 0 |
| R6 | 2 | 28304 | 22552 | 17858 | 10446 | 4694 |
| ALLO | 1 | 930 | 930 | 930 | 0 | 0 |
| ALLO | 2 | 5471 | 5471 | 5471 | 0 | 0 |
| ALLO6 | 1 | 150 | 150 | 150 | 0 | 0 |
| ALLO6 | 2 | 4750 | 4352 | 4352 | 398 | 0 |
| ALLO6 | 3 | 15756 | 13462 | 13462 | 2294 | 0 |

Every new population except R6 level 2 is a subset of the old. The test suite pins 62, 369 and
1,200 for R4 and 62, 1,409 and 22,552 for R6.

## Removed mechanisms

Removed: 13,881.

| Set | Zero-flux | Twin-only copy group | Other order | Unexplained |
|---|---|---|---|---|
| R1 | 0 | 0 | 0 | 0 |
| R2 | 0 | 0 | 0 | 0 |
| R3 | 0 | 7 | 0 | 0 |
| R4 | 188 | 0 | 0 | 0 |
| R5 | 188 | 0 | 0 | 0 |
| R6 | 188 | 10618 | 0 | 0 |
| ALLO | 0 | 0 | 0 | 0 |
| ALLO6 | 0 | 2681 | 11 | 0 |

- Zero-flux (564): the rapid-equilibrium twin of each is in the new population, 564 of 564.
  The 188 removed R4 mechanisms all have rank at least 1 and fitted minus rank at least their
  zero-flux group count (fitted 8 to 9, rank 6 to 8); 0 failures.
- Twin-only copy groups (13,306): every site of a copy group duplicates a productive form in
  every state where the copy binds. R3 loses all 7 dead-end children: a uni-uni mechanism has no
  new complex to give.
- Other order (11, ALLO6 level 3): mechanisms whose level-2 parents were all removed. For each,
  the mechanism without the copy is in the new level 2, and the dead-end move on it returns the
  mechanism: a level-3 path ALLO6 does not enumerate.
- ALLO6 by level: level 2 loses 398, all with a `_duplicate_copy_groups` hit; level 3 loses
  2,294, 2,283 with one. Of the level-3 removals, 417 have a level-2 parent kept in the new
  population: the tag-relaxation filter dropped them, all 417 with a hit. The first eight by key
  add 1, 0, 1, 1, 1, 1, 1 and 1 phantoms over the kept parent. All 417 were ranked on the
  previous snapshot, whose ALLO6 population this one matches mechanism for mechanism: 365 add
  one phantom, 9 add two and 43 add none.
- R1, R2 and ALLO lose nothing. None is unexplained.

R6 level 2 keeps 368 two-copy mechanisms that the previous snapshot could not reach. In each,
one copy's group is valid only through a site at the other copy's form. In the first of them, A*
binds at {E, E(B), E(Q), E(Q*)} and Q* at {E, E(P), E(A*)}. Q* at E duplicates E(Q) and at E(P)
duplicates E(P, Q); its complex at E(A*), E(A*, Q*), shares its composition only with E(A*, Q),
a copy's complex, so the Q* group is not twin-only. The parent without Q* is in the new level 1,
and on it the dead-end move places Q* at {E, E(P), E(A*)} while Q* competes with A, B and Q:
E(A*) carries no productive A. The mirror of E + A* ⇌ E(A*) onto E(Q*) + A* ⇌ E(A*, Q*) joins
the A* group. The dead-end move reaches all 368 from a parent in the new level 1. Against that
parent, 144 add one identifiable constant and 224 add a phantom (ranked on the previous
snapshot, with the same parents).

Against the previous snapshot (87f49a7), this one removes none and adds 5,062, all in R6 level
2; 368 of them are the mechanisms above, which the baseline holds as well.

## Added mechanisms

Added against a2a02b1: 4,694, all in R6 level 2, none violating a rule. Each binds two copies,
the second placed where competition by name refused it.

| Second copy binds | Count |
|---|---|
| at a form that carries the first copy | 4,054 |
| at a form that holds the first copy's reactant productively | 612 |
| at neither | 28 |

In the two inspected of the 28, the second copy competes with a reactant whose copy the parent
binds; competition by name also targets the copy's sites, so these site sets were out of reach.

The rank oracle on a draw of 150:

| Ranked | Reached from a parent in the new population | Adds a phantom over every such parent | Over some, not all | Over none | Carries ≥ 1 phantom |
|---|---|---|---|---|---|
| 150 | 150 | 3 | 0 | 147 | 4 |

The 150 have 260 parents in the new population; 256 of the comparisons add no phantom and 4 add
one. Each of the three that add a phantom binds two copies whose groups are twin-only when
copies' complexes count, the class the rule tolerates (two copies that are twins only of each
other). The fourth phantom is inherited from the parent. 3 of 150 (2.0%) add a phantom.

## Kept-mechanism rules

`_assert_emission_rules` holds on every kept and added mechanism:

| Set | Checked | Violating |
|---|---|---|
| R1 | 4 | 0 |
| R2 | 7 | 0 |
| R3 | 4 | 0 |
| R4 | 1631 | 0 |
| R5 | 4527 | 0 |
| R6 | 24023 | 0 |
| ALLO | 6401 | 0 |
| ALLO6 | 17964 | 0 |

The rules do not cover a `:NonequalAI` copy group, one constant per state, that binds only at
twin sites in the active state, or only at twin sites among its inactive-state sites: 628 of
the 17,964 ALLO6 mechanisms hold one. Fitted and rank for the first eight by key: 9 and 9, 8
and 6, 8 and 6, 8 and 8, 9 and 9, 8 and 7, 8 and 7, 8 and 8.

## Derivation of kept mechanisms

| Set | Compared | Differing | Errors |
|---|---|---|---|
| R1 | 4 | 0 | 0 |
| R2 | 7 | 0 | 0 |
| R3 | 4 | 0 | 0 |
| R4 | 40 | 0 | 0 |
| R5 | 40 | 0 | 0 |
| R6 | 40 | 0 | 0 |
| ALLO | 40 | 0 | 0 |
| ALLO6 | 40 | 0 | 0 |

215 mechanisms, 0 differences, 0 errors.

## Level-3 phantom fraction

R6 levels 0–2 hold 62, 1,409 and 22,552; 2,118 level-2 mechanisms hold a copy group that is
twin-only when complexes of copies count. The 300 drawn level-2 parents have 984 split and 8,959
dead-end children. 2,494 of these 9,943 children hold a copy group whose every site duplicates
some form when every form counts and not when only productive forms count (split 134, dead-end
2,360); 817 of the 2,494 have a parent that already holds such a group.

| Move | Parent already holds one | Ranked | Adds a phantom over parent | Does not | Carries ≥ 1 phantom |
|---|---|---|---|---|---|
| split | no | 3 | 1 | 2 | 1 |
| split | yes | 10 | 2 | 8 | 3 |
| dead-end | no | 100 | 8 | 92 | 8 |
| dead-end | yes | 37 | 0 | 37 | 5 |
| all | — | 150 | 11 | 139 | 17 |

11 of the 150 ranked children (7.3%) add a phantom over their parent, and one has one phantom
fewer than its parent; scaled to all children of the drawn parents, about 1.84% of level-3
children are phantoms that the productive-only reading admits. The 8 ping-pong level-2 split
children (R5 2, R6 6) are each 8 fitted and rank 8 over a parent of 7 and 7: none adds a
phantom.

## Enumeration time

| Set | Old (s) | New (s) | Change |
|---|---|---|---|
| R1 | 0.89 | 0.77 | -14% |
| R2 | 0.056 | 0.058 | 3% |
| R3 | 0.005 | 0.0003 | -94% |
| R4 | 0.44 | 0.48 | 8% |
| R5 | 1.29 | 1.40 | 8% |
| R6 | 7.84 | 6.40 | -18% |
| ALLO | 53.06 | 53.24 | 0% |
| ALLO6 | 100.40 | 91.38 | -9% |

One run per checkout. R6 enumerates 18% faster; no population is more than 8% slower.

## Conclusion

Against a2a02b1 the filters remove 13,881 mechanisms: 564 with a steady-state group that
carries no flux, whose rapid-equilibrium twin stays; 13,306 with a copy group whose every
complex duplicates a productive form; and 11 ALLO6 mechanisms the new code reaches in another
order. None is unexplained. Per-site competition adds 4,694 R6 level-2 mechanisms with two
copies, each satisfying both rules; 3 of a draw of 150 add a phantom over their parents, each
two copies that are twins only of each other. Every kept mechanism satisfies both rules, and the
215 sampled derivations are unchanged. At level 3, about 1.84% of the children of a sample of
R6 level-2 parents are phantoms that the productive-only reading admits.
