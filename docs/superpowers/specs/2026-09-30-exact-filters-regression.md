# Exact filters: regression check

Date: 2026-10-01. This checks the exact enumeration filters (design
`2026-09-30-exact-filters-design.md`, plan `../plans/2026-09-30-exact-filters.md`): the flux
predicate, the flip's extend semantics, the twin predicate, the dead-end move's new-complex
rule, the split's revert and its copy-unit filter (judged against all forms), the
mechanism-level check against copy-free forms and the tag-relaxation filter. Every mechanism
the filters remove must be one a rule names; every kept and added mechanism must satisfy both
emission rules; every kept mechanism must derive as before.

## Method

A scratch script (not committed) enumerates the populations of one checkout and writes them
as keys; a second mode compares two key sets mechanism by mechanism. The baseline ran on
a2a02b1, before the filters; the comparison ran on 0cb540f, which holds them. Both enumerated
independently, one Julia process at a time. The a2a02b1 checkout ran with the current
`Manifest.toml`, so both runs used the same dependency versions.

Populations: `init_mechanisms` plus two levels of the flip (`_expand_re_to_ss`), split
(`_expand_split_kinetic_group`) and dead-end (`_expand_add_dead_end_regulator`) moves on
reactions R1–R6 of the findings document, deduplicated across levels; ALLO, the allosteric
children of the R4 seeds (`_expand_to_allosteric`) with their flip and split children; and
ALLO6, which carries copies into allosteric mechanisms. R1 is uni-uni, R2 uni-uni with one
dead-end inhibitor, R3 uni-uni with S and P as inhibitors, R4 the bi-bi reaction whose atoms
admit ping-pong, R5 R4 with one dead-end inhibitor, R6 R4 with A, B, P and Q as inhibitors.

ALLO6 has three levels. Level 1 is 150 of the allosteric children of the R6 seeds, drawn
with `MersenneTwister(20260930)` from the children sorted by key, so both checkouts draw the
same ones. Level 2 is their children under the dead-end, tag-relaxation
(`_expand_change_allo_state`), flip and split moves. Level 3 is the tag relaxations of the
level-2 mechanisms that bind a copy; the R6 seeds bind none, so only level 3 exercises the
tag-relaxation filter on a copy. Each level-3 child's level-2 parents are recorded.

A mechanism is keyed by its reaction's sig string and the sig strings of its steps
(`_to_sig`), group by group, with an allosteric mechanism's conformational states and
catalytic multiplicity appended; it is rebuilt by parsing each piece (`_step_from_sig`,
`_reaction_from_sig`). Both checkouts produce the same key. Caching the per-step strings
made keying 1,000 R6 mechanisms take 0.018 s instead of 88.7 s, and each snapshot took about
7 minutes.

For each removed mechanism the script records the first rule that names it:

- a steady-state group that carries no flux (`_flux_carrying_groups`);
- else a copy group whose every site is a twin among copy-free forms
  (`_duplicate_copy_groups`, the mechanism-level check);
- else a copy group whose every site is a twin among all forms of the mechanism without that
  copy, other copies' complexes included: the dead-end move's placement rule, which judges
  a new copy against every form of its parent;
- else a copy group whose every site is a twin among all forms of the mechanism, the copy's
  own complexes at its other sites included (`_twin_only` with `_copy_twin_test(m)`): the
  split's part filter;
- else, for an allosteric mechanism, a valid mechanism that the new code reaches by placing
  the copy on its relaxed parent (`_expand_add_dead_end_regulator` on the mechanism without
  that copy returns it);
- else UNEXPLAINED.

For a zero-flux removal it checks that the mechanism with those groups at rapid equilibrium
(its RE twin) is in the new population, and for the removed R4 zero-flux mechanisms that the
numerical rank of ∂v/∂log θ (finite differences with step `h = 1e-5` over 60 points, three
draws, singular values above 1e-7 of the largest) is at least 1 and fitted parameters minus
the rank is at least the number of zero-flux groups. For ALLO6 it counts, per level, the
removals with a `_duplicate_copy_groups` hit under the new code; the level-3 removals with a
level-2 parent kept in the new population, which the tag-relaxation filter dropped; and,
with the same rank, the phantoms each level-2 placement-rule removal adds over the mechanism
without its copy, and each of the first eight filtered level-3 removals by key adds over its
kept parent.

It runs `_assert_emission_rules` on every kept and every added mechanism of every set. It
also records the enumeration time of each population (one run per checkout).

It also compares the new keys with the snapshot of the previous comparison, which ran on
9215223, where the split judged its parts against copy-free forms. For each mechanism that
snapshot holds and this one does not, it records the two twin predicates, the fitted
parameters and the rank.

The derivation comparison is a separate step. It draws, per set, up to 40 keys (all of them
for R1–R3) from the sorted intersection of the two key sets with `MersenneTwister(20260930)`,
then derives each on both checkouts (`fitted_params` joined by spaces and the Reduced
`rate_equation_string`) and compares the two.

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
| R5 | 2 | 3934 | 3744 | 3744 | 190 | 0 |
| R6 | 0 | 62 | 62 | 62 | 0 | 0 |
| R6 | 1 | 1769 | 1409 | 1409 | 360 | 0 |
| R6 | 2 | 28304 | 16980 | 16980 | 11324 | 0 |
| ALLO | 1 | 930 | 930 | 930 | 0 | 0 |
| ALLO | 2 | 5471 | 5471 | 5471 | 0 | 0 |
| ALLO6 | 1 | 150 | 150 | 150 | 0 | 0 |
| ALLO6 | 2 | 4750 | 4245 | 4245 | 505 | 0 |
| ALLO6 | 3 | 15756 | 13088 | 13088 | 2668 | 0 |

Every new population is a subset of the old. The test suite pins 62, 369 and 1,200 for R4
and 62, 1,409 and 16,980 for R6.

## Removed mechanisms

Removed: 15,242, each named by a rule; none is unexplained.

| Set | Zero-flux | Twin-only (copy-free) | Placement rule | Twin-only part (all forms) | Other order | Unexplained |
|---|---|---|---|---|---|---|
| R1 | 0 | 0 | 0 | 0 | 0 | 0 |
| R2 | 0 | 0 | 0 | 0 | 0 | 0 |
| R3 | 0 | 7 | 0 | 0 | 0 | 0 |
| R4 | 188 | 0 | 0 | 0 | 0 | 0 |
| R5 | 188 | 0 | 0 | 2 | 0 | 0 |
| R6 | 188 | 10618 | 872 | 6 | 0 | 0 |
| ALLO | 0 | 0 | 0 | 0 | 0 | 0 |
| ALLO6 | 0 | 2681 | 346 | 0 | 146 | 0 |

- Zero-flux removals (564): the RE twin of each is in the new population, 564 of 564. For
  the 188 removed R4 mechanisms the rank is at least 1 and fitted parameters minus the rank
  is at least the number of zero-flux groups in all 188; 0 failures. Fitted parameters range
  over 8 to 9 and the rank over 6 to 8.
- Twin-only copy groups (13,306): dead-end children whose copy binds only where its complex
  duplicates a copy-free form. R3 loses all 7 of its dead-end children, as a uni-uni
  mechanism has no new complex to give.
- Placement rule (1,218): dead-end children that the copy-free check passes and the
  placement rule rejects. In R6 (872) the copy duplicates a form of its parent at every
  site, and at some site the only such forms are other copies' complexes. In ALLO6 (346: 107
  at level 2, 239 at level 3), whose level-1 parents bind no copy, the two checks read the
  inactive state differently. Measured against the mechanism without the copy, 62 of the 107
  level-2 ones add no phantom and 45 add one. In the first of them, the parent's inactive
  conformation keeps no step, so the placement rule finds free E absent from that state and
  judges the copy by the active state alone, while in the child the copy binding reaches the
  inactive state and the copy is new there.
- Twin-only parts against all forms (8: R5 2, R6 6): split children of a ping-pong mechanism
  whose second chemistry step E(B; res) ⇌ E(Q) is at rapid equilibrium, with the copy (I in R5,
  B or Q in R6) bound at E(Q) and at E(B; res) and the mirrored isomerization between its two
  complexes in a kinetic group of its own. The split separates the two sites, and the part
  holding one of them binds only at twin sites; at that site the only twin is the copy's own
  complex in the other part, with the same segment and offsets. Each has 8 fitted parameters and
  rank 8: identifiable mechanisms that the all-forms reading does not emit, the cost stated when
  Denis chose this reading (2026-10-01). A scratch probe judged the same parts with every form a
  twin source except the complexes of the copy being split: all eight then pass, while the
  depth-4 phantom of the test suite (A* split down to E(Q) beside Q* at E(A)) and the B* split
  of its neighboring test are still rejected.
- R5's only inhibitor is foreign. No copy-free form binds it, and the placement rule drops
  every form that does, so a foreign inhibitor's complex is never a twin under either check:
  188 of R5's 190 removals are zero-flux. The split's all-forms reading can match the
  inhibitor's complex to its own complex at another site, which removes the other 2. These
  figures are measured per set, not derived from the totals.
- Other order (146, ALLO6 level 3): valid mechanisms whose level-2 parents were all removed;
  the new code reaches each by placing the copy after the relaxation, a path ALLO6 does not
  enumerate.
- ALLO6 by level: level 2 loses 505, of which 398 have a `_duplicate_copy_groups` hit; level
  3 loses 2,668, of which 2,283 have one. Of the level-3 removals, 310 have a level-2 parent
  kept in the new population: the tag-relaxation filter dropped them, and all 310 have a
  `_duplicate_copy_groups` hit. Over the first eight of them by key, each relaxation adds one
  phantom over its kept parent. The other 2,358 lost every parent.
- R1, R2 and ALLO lose nothing.

Against the previous comparison's snapshot (9215223), this one loses exactly the eight
twin-only parts above, two at R5 level 2 and six at R6 level 2, and gains none; every other
set and level holds the same keys. Each of the eight has a copy group that is twin-only
against all forms and none against copy-free forms.

## Added mechanisms

None: 0 added.

## Kept-mechanism rules

`_assert_emission_rules` holds on every kept and added mechanism:

| Set | Checked | Violating |
|---|---|---|
| R1 | 4 | 0 |
| R2 | 7 | 0 |
| R3 | 4 | 0 |
| R4 | 1631 | 0 |
| R5 | 4525 | 0 |
| R6 | 18451 | 0 |
| ALLO | 6401 | 0 |
| ALLO6 | 17483 | 0 |

The rules do not cover one class the script also counted: a `:NonequalAI` copy group, which
has one constant per conformational state, that binds only at twin sites in the active state,
or only at twin sites among the sites it binds in the inactive state. 521 of the 17,483 ALLO6
mechanisms hold one. Fitted parameters and rank for the first eight by key: 9 and 9, 8 and 6,
8 and 6, 8 and 8, 9 and 9, 8 and 7, 8 and 7, 8 and 8. The parents were not ranked, so these
do not attribute a deficit to the copy; the test suite's tag-relaxation case attributes one
(relaxing the kept copy of S to `:NonequalAI` gives 6 fitted parameters and rank 5, its parent
5 and 5).

## Derivation of kept mechanisms

The sample is drawn from the sorted intersection of the kept keys, seed 20260930, up to 40
per set, all of R1–R3. Fitted names and Reduced string, old against new:

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

## Enumeration time

| Set | Old (s) | New (s) | Change |
|---|---|---|---|
| R1 | 0.89 | 0.89 | 0% |
| R2 | 0.056 | 0.069 | 24% |
| R3 | 0.005 | 0.0003 | -93% |
| R4 | 0.44 | 0.52 | 17% |
| R5 | 1.29 | 1.48 | 15% |
| R6 | 7.84 | 5.33 | -32% |
| ALLO | 53.06 | 55.00 | 4% |
| ALLO6 | 100.40 | 94.62 | -6% |

Each time is one run per checkout; there were no repeats. R6 enumerates 32% faster and ALLO6
6% faster. R4 is 17% and R5 15% slower. The previous comparison ran the same populations on
code that differs from this one's only in the split's twin sources and two refactors that
change no behavior, and measured 0.44 s and 1.36 s. R4 binds no copy, so the split's twin
test never reaches a site there, and its rise is variation between single runs. R2 is 24%
slower on a population that takes 0.07 s.

## Conclusion

The filters remove 15,242 mechanisms across R1–R6, ALLO and ALLO6, each named by a rule: 564
with a steady-state group that carries no flux, whose rapid-equilibrium twin stays in the
population; 13,306 whose copy duplicates a copy-free form at every site; 1,218 the dead-end
placement rule rejects; 8 split children with a part that binds only at twins among all
forms; and 146 valid ALLO6 mechanisms the new code reaches in another order. The
tag-relaxation filter drops 310 ALLO6 relaxations, each with a twin-only copy group. The
filters add none, every kept mechanism satisfies both emission rules, and no kept mechanism's
derivation changes in the 215-mechanism sample. Three questions remain open for Denis: the
107 ALLO6 placements that the placement rule and the child-level check judge differently (62
of them identifiable), the `:NonequalAI` copy class the rules do not cover, and the eight
identifiable split children the all-forms reading drops, which a reading that excludes only
the split copy's own complexes would keep.
