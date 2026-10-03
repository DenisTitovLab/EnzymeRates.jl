# Exact filters: regression check

Date: 2026-10-02. This checks the exact enumeration filters (design
`2026-09-30-exact-filters-design.md`, plan `../plans/2026-09-30-exact-filters.md`): the flux
predicate, the flip's extend semantics, the split's revert, the copy rule (section 2 of the
design: a copy group is rejected only when every complex has a productive twin and the dwell
gauge of Theorem 2 is consistent) and the dead-end move's per-site competition (section 5).
Every mechanism the filters remove must be one a rule names; every kept and added mechanism must
satisfy both emission rules; every kept mechanism must derive as before; added mechanisms are
reported with their phantom fraction; the rejections are checked against the rank oracle.

## Method

A scratch script (not committed) enumerates the populations of one checkout and writes them as
keys; a second mode compares two key sets mechanism by mechanism. The baseline ran on a2a02b1,
before the filters; the comparison ran on the tree committed as c3481a4, the branch commit that
holds the gauge rule (the commit adds only a type assertion and the wording of an error message
to the enumeration code that ran). Both enumerated independently, one Julia process at a time;
the a2a02b1 run is the snapshot of the earlier comparisons, made with the current
`Manifest.toml`. The previous branch snapshot, ea0ada5, decided copies by the all-twin rule
(every complex a twin, no gauge).

The copy rule. A copy complex's twin is a productive complex, a form bound to no competitive
inhibitor, in the complex's rapid-equilibrium segment with the complex's offsets, so the two
weights are proportional (`_productive_twin`). A copy group is redundant when every complex has
a twin in every conformational state where the copy binds and the dwell gauge is consistent
over the states together, an `:EqualAI` group taking one rescaling in both
(`_redundant_copy_groups`). Every conformation holds its free enzyme.

Competition in the dead-end move is decided per site: for competition with a reactant the move
targets the forms where that reactant binds productively, for competition with an inhibitor
the forms where that inhibitor binds; it excludes a form that holds a competing reactant
productively or a competing inhibitor, and its capacity test counts productive bindings only.
a2a02b1 decided competition by name.

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
group that carries no flux (`_flux_carrying_groups`); else a redundant copy group
(`_redundant_copy_groups`); else a mechanism the new code reaches by placing a copy on its
parent without that copy, the parent being in the new population (other order); else a
mechanism reachable only through a redundant intermediate: for each copy, the mechanism without
it, or with its split groups joined, has a redundant copy group, or (ALLO6 level 3) every
level-2 parent the baseline recorded has one; else UNEXPLAINED. Every mechanism of the last
named category is ranked. For a zero-flux removal the script checks that the mechanism's
rapid-equilibrium twin is in the new population, and for the removed R4 zero-flux mechanisms
that the numerical rank of ∂v/∂log θ (finite differences, h = 1e-5, 60 points, three draws,
singular values above 1e-7 of the largest) is at least 1 and fitted parameters minus rank at
least the number of zero-flux groups. It runs `_assert_emission_rules` on every kept and added
mechanism and records each population's enumeration time.

Ranked mechanisms are compared with their parents in the new population, the mechanism without
one of its copies; a mechanism adds a phantom when its fitted count minus rank exceeds its
parent's. Two oracle checks run on R6 levels 1–2 with `MersenneTwister(20260930)`: (a) 100
kept mechanisms with a copy group whose every complex duplicates a productive form under the
all-twin rule, which ea0ada5 rejected, are ranked; (b) 100 mechanisms removed for a redundant
copy group must all have fitted minus rank at least 1, in R6 and in ALLO6.

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
| R3 | 2 | 6 | 3 | 3 | 3 | 0 |
| R4 | 0 | 62 | 62 | 62 | 0 | 0 |
| R4 | 1 | 369 | 369 | 369 | 0 | 0 |
| R4 | 2 | 1388 | 1200 | 1200 | 188 | 0 |
| R5 | 0 | 62 | 62 | 62 | 0 | 0 |
| R5 | 1 | 719 | 719 | 719 | 0 | 0 |
| R5 | 2 | 3934 | 3746 | 3746 | 188 | 0 |
| R6 | 0 | 62 | 62 | 62 | 0 | 0 |
| R6 | 1 | 1769 | 1649 | 1649 | 120 | 0 |
| R6 | 2 | 28304 | 31870 | 25772 | 2532 | 6098 |
| ALLO | 1 | 930 | 930 | 930 | 0 | 0 |
| ALLO | 2 | 5471 | 5471 | 5471 | 0 | 0 |
| ALLO6 | 1 | 150 | 150 | 150 | 0 | 0 |
| ALLO6 | 2 | 4750 | 4609 | 4609 | 141 | 0 |
| ALLO6 | 3 | 15756 | 14980 | 14980 | 776 | 0 |

Every new population except R6 level 2 is a subset of the old. The test suite pins 62, 369 and
1,200 for R4 and 62, 1,649 and 31,870 for R6.

## Removed mechanisms

Removed: 3,950.

| Set | Zero-flux | Redundant copy group | Other order | Reachable only through a redundant intermediate | Unexplained |
|---|---|---|---|---|---|
| R1 | 0 | 0 | 0 | 0 | 0 |
| R2 | 0 | 0 | 0 | 0 | 0 |
| R3 | 0 | 5 | 0 | 0 | 0 |
| R4 | 188 | 0 | 0 | 0 | 0 |
| R5 | 188 | 0 | 0 | 0 | 0 |
| R6 | 188 | 2286 | 0 | 178 | 0 |
| ALLO | 0 | 0 | 0 | 0 | 0 |
| ALLO6 | 0 | 776 | 20 | 121 | 0 |

- Zero-flux (564): the rapid-equilibrium twin of each is in the new population, 564 of 564.
  The 188 removed R4 mechanisms all have rank at least 1 and fitted minus rank at least their
  zero-flux group count (fitted 8 to 9, rank 6 to 8); 0 failures.
- Redundant copy groups (3,067): every complex has a productive twin in every state where the
  copy binds and the gauge is consistent. R3 loses five copies at E of uni-uni mechanisms
  whose copied ligand binds at rapid equilibrium.
- Other order (20, ALLO6 level 3): mechanisms whose level-2 parents were all removed; for
  each, the mechanism without the copy is in the new level 2, and the dead-end move on it
  returns the mechanism.
- Reachable only through a redundant intermediate (299), the rule admits each but no move
  reaches it. R6, 178: 132 split a copy group whose joined group is redundant, and 46 bind two
  copies each redundant without the other. ALLO6, 121: level-3 relaxations of a copy to
  `:NonequalAI`, whose every recorded level-2 parent holds the copy at `:EqualAI` as a
  redundant group; the per-state constants break the ties that made the parent redundant.
  Ranked all 299: every one has fitted minus rank at least 1; against the mechanism without
  its copies, 199 carry one phantom and 100 two.
- ALLO6 by level: level 2 loses 141, all with a redundant group; level 3 loses 776, 635 with
  one. Of the level-3 removals, 115 have a level-2 parent kept in the new population: the
  tag-relaxation filter dropped them, all 115 with a redundant group; the first eight by key
  add one phantom each over the kept parent.
- R1, R2 and ALLO lose nothing. None is unexplained.

Against the previous snapshot (ea0ada5), this one removes none and adds 11,335: R3 level 2 2,
R6 level 1 240, R6 level 2 9,318, ALLO6 level 2 257, ALLO6 level 3 1,518. 9,931 of them are in
the baseline as well.

## Added mechanisms

Added against a2a02b1: 6,098, all in R6 level 2, none violating a rule. Each binds two copies,
5,110 with a copy group bound at a form that carries the other copy: per-site competition
reaches them, and the gauge returns those whose placement it no longer rejects. The rank oracle
on a draw of 150:

| Ranked | Reached from a parent in the new population | Adds a phantom over every such parent | Over some, not all | Over none | Carries ≥ 1 phantom |
|---|---|---|---|---|---|
| 150 | 150 | 2 | 0 | 148 | 2 |

The 150 have 284 parents in the new population; 281 of the comparisons add no phantom and 3
add one.

## Oracle checks

- (a) R6 mechanisms kept with a copy group whose every complex duplicates a productive form
  under the all-twin rule: 240 at level 1 and 9,214 at level 2. Of 100 ranked, **90 have full
  rank** and 10 do not; 8 of the 10 add a phantom over their parents (the copy's gauge fails,
  yet no other constant pins it), 2 inherit it. The test is conservative: it rejects only on a
  proof.
- (b) Removed for a redundant copy group: of 100 ranked R6 mechanisms (of 2,286), **100 have
  fitted minus rank at least 1**, and in all 100 the copy adds at least one phantom over the
  mechanism without it; ALLO6, 100 of 776 ranked, the same 100 and 100. None of the ranked R6
  rejections has a twin outside Theorem 2's reach.
- Twins of the same composition formed by a steady-state binding are not twins. Taken as twins
  they rejected 20 R6 mechanisms, 4 of them identifiable (an SS binding of A, P, Q or B with its
  copy at E and a random-order product or substrate side: 7 fitted, rank 7). The 20 return (4
  full rank, 16 with a phantom), and so do two R3 copies at E of a uni-uni mechanism with S or
  P bound at steady state, each with a phantom.
- Ties across conformations: 23 ALLO6 mechanisms in the new population would be redundant if
  each state were judged alone; ranked all, each copy adds a phantom over the mechanism without
  it. The 121 unreachable relaxations above are the other mechanisms the ties change.

## Phantom intake

Estimated from draws of 150 per stratum (`MersenneTwister(20260930)`), mechanisms added against
ea0ada5 that add a phantom over every parent in the new population:

| Stratum | Added | Ranked | Add a phantom | Estimate |
|---|---|---|---|---|
| R6 level 1 | 240 | 150 | 2 | about 3 |
| R6 level 2 | 9,318 | 150 | 23 | about 1,430 (15%) |
| ALLO6 level 2 | 257 | 150 | 18 | about 31 |
| ALLO6 level 3 | 1,518 | 150 | 65 | about 660 |

The per-site wave's level-2 intake, also an estimate: 224 of the 368 mechanisms it restored add
a phantom (measured) and about 2% of its other 4,694 additions do (sampled), about 6% of its
5,062 additions.

## Kept-mechanism rules

`_assert_emission_rules` holds on every kept and added mechanism:

| Set | Checked | Violating |
|---|---|---|
| R1 | 4 | 0 |
| R2 | 7 | 0 |
| R3 | 6 | 0 |
| R4 | 1631 | 0 |
| R5 | 4527 | 0 |
| R6 | 33581 | 0 |
| ALLO | 6401 | 0 |
| ALLO6 | 19739 | 0 |

The rules do not cover a `:NonequalAI` copy group, one constant per state, that binds only at
twin sites in the active state, or only at twin sites among its inactive-state sites: 885 of
the 19,739 ALLO6 mechanisms hold one. Fitted and rank for the first eight by key: 9 and 9, 8
and 6, 8 and 6, 8 and 8, 9 and 9, 8 and 7, 8 and 7, 8 and 8.

## Derivation of kept mechanisms

| Set | Compared | Differing | Errors |
|---|---|---|---|
| R1 | 4 | 0 | 0 |
| R2 | 7 | 0 | 0 |
| R3 | 6 | 0 | 0 |
| R4 | 40 | 0 | 0 |
| R5 | 40 | 0 | 0 |
| R6 | 40 | 0 | 0 |
| ALLO | 40 | 0 | 0 |
| ALLO6 | 40 | 0 | 0 |

217 mechanisms, 0 differences, 0 errors.

## Level-3 phantom fraction

R6 levels 0–2 hold 62, 1,649 and 31,870. 11,332 level-2 mechanisms hold a copy group that is
twin-only when complexes of copies count, 8,854 a copy group whose every complex has a
productive twin and whose gauge fails. The 300 drawn level-2 parents have 1,290 split and
10,855 dead-end children. 6,609 of these 12,145 children are admitted only by a relaxation of
the all-twin rule over all forms (split 698, dead-end 5,911): a copy group twin-only when
copies' complexes count, or all-twin with a failing gauge (4,341); 3,951 of the 6,609 have a
parent that already holds such a group.

| Move | Parent already holds one | Ranked | Adds a phantom over parent | Does not | Carries ≥ 1 phantom |
|---|---|---|---|---|---|
| split | no | 6 | 0 | 6 | 0 |
| split | yes | 9 | 3 | 6 | 4 |
| dead-end | no | 46 | 5 | 41 | 5 |
| dead-end | yes | 89 | 6 | 83 | 10 |
| all | — | 150 | 14 | 136 | 19 |

14 of the 150 ranked children (9.3%) add a phantom over their parent, and three have one
phantom fewer; 9 of the 14 hold an all-twin copy group with a failing gauge. Scaled to all
children of the drawn parents, about 5.1% of level-3 children are phantoms that a relaxation
admits. Level-2 split children so admitted: R5 2 and R6 1,180 (1,174 all-twin with a failing
gauge); of 20 ranked R6 ones, none adds a phantom.

## Seeds with required copies

`seed_mechanisms(R6, ∅, {A, Q})`: 3,569 seeds in 40 s, none violating a rule; 993 hold a copy
group whose every complex has a productive twin and whose gauge fails. Of 20 ranked, 18 have
full rank and 2 are 7 fitted, rank 6. (ea0ada5: 2,550 seeds.)

## Enumeration time

| Set | Old (s) | New (s) | Change |
|---|---|---|---|
| R1 | 0.89 | 0.80 | -10% |
| R2 | 0.06 | 0.06 | 0% |
| R3 | 0.00 | 0.27 | — |
| R4 | 0.44 | 0.48 | 8% |
| R5 | 1.29 | 1.32 | 3% |
| R6 | 7.84 | 9.58 | 22% |
| ALLO | 53.06 | 55.23 | 4% |
| ALLO6 | 100.40 | 104.05 | 4% |

One run per checkout. R6 warm (best of three runs after a warm-up): a2a02b1 7.58 s, ea0ada5
6.33 s, the gauge rule 9.86 s, a rise of 30% over a2a02b1 for 11% more mechanisms; per
enumerated mechanism 0.252, 0.264 and 0.294 ms. The copy rule's own checks cost about 0.25 s
in the split (0.14 s for 980 one-unit children of the part filter, 0.10 s for the children
filter); the rest follows the larger population (16% more level-1 parents than ea0ada5, 43%
more dead-end children). R3's 0.27 s is first-use compilation.

## Conclusion

Against a2a02b1 the filters remove 3,950 mechanisms: 564 with a steady-state group that carries
no flux, whose rapid-equilibrium twin stays; 3,067 with a redundant copy group; 20 ALLO6
mechanisms the new code reaches in another order; and 299 that the rule admits but no move
reaches, each with at least one phantom. None is unexplained. All 100 ranked rejections in R6
and in ALLO6 have a phantom; 90 of 100 ranked kept all-twin mechanisms have full rank. Every
kept mechanism satisfies both rules, and the 217 sampled derivations are unchanged. The gauge
returns the shared-group family Denis chose to keep: R6 level 2 grows from 22,552 to 31,870,
about 15% of its additions adding a phantom (estimate), and at level 3 about 5.1% of the
children of a sample of level-2 parents are phantoms that a relaxation admits.
