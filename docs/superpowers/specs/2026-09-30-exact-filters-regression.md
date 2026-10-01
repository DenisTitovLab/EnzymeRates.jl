# Exact filters: regression check

Date: 2026-09-30. This checks the exact enumeration filters (design
`2026-09-30-exact-filters-design.md`, plan `../plans/2026-09-30-exact-filters.md`): the flux
predicate, the flip's extend semantics, the twin predicate, the dead-end move's new-complex
rule and the split's revert and copy-unit filter. Every mechanism the filters remove must be
one the rules name; every kept mechanism must derive as before; nothing may be added.

## Method

A scratch script (not committed) enumerates the populations of one checkout and writes them
as keys; a second mode compares two key sets mechanism by mechanism. The baseline ran on
a2a02b1, before the filters; the comparison ran on 33f4a3d, which holds the filters. Both
enumerated independently, one Julia process at a time.

Populations: `init_mechanisms` plus two levels of the flip (`_expand_re_to_ss`), split
(`_expand_split_kinetic_group`) and dead-end (`_expand_add_dead_end_regulator`) moves on
reactions R1–R6 of the findings document, deduplicated across levels; and ALLO, the
allosteric children of the R4 seeds (`_expand_to_allosteric`) with their flip and split
children. R1 is uni-uni, R2 uni-uni with one dead-end inhibitor, R3 uni-uni with S and P as
inhibitors, R4 the bi-bi reaction whose atoms admit ping-pong, R5 R4 with one dead-end
inhibitor, R6 R4 with A, B, P and Q as inhibitors.

A mechanism is keyed by `_sig_of` (an allosteric one also by its conformational states and
catalytic multiplicity), a key both versions of the code produce. For each removed mechanism
the script checks:

- the groups at steady state that carry no flux (`_flux_carrying_groups`), or, if none, the
  copy groups whose twin test fails (`_duplicate_copy_groups`); a removed mechanism with
  neither is UNEXPLAINED;
- for a zero-flux removal, that the mechanism with those groups at rapid equilibrium (its RE
  twin) is in the new population (TWIN MISSING otherwise);
- for the removed R4 zero-flux mechanisms, that fitted parameters minus the numerical rank of
  ∂v/∂log θ (finite differences over 60 points, three draws, singular values above 1e-7 of
  the largest) is at least the number of zero-flux groups (RANK otherwise).

For each added mechanism it runs `_assert_emission_rules`. It also records the enumeration
time of each population (one run per checkout). The compare mode ran on the new checkout, so
it uses the new predicates; the snapshot mode uses only functions both checkouts have. The
rank is a finite-difference one with step `h = 1e-5`, and a removed R4 mechanism passes only
if its rank is at least 1 and fitted minus rank is at least its number of zero-flux groups.
The 11,505 copy removals were classified by the predicate under test
(`_duplicate_copy_groups`); no independent rank check covers them.

The derivation comparison is a separate step. It draws, per set, up to 40 keys (all of them
for R1–R3) from the sorted intersection of the two key sets with `MersenneTwister(20260930)`,
then derives each on both checkouts (`fitted_params` joined by spaces and the Reduced
`rate_equation_string`) and compares the two. The scripts lived in the session scratchpad
and are not committed. The a2a02b1 checkout ran with the current `Manifest.toml`, so both
runs used the same dependency versions.

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

Every new population is a subset of the old. R4's levels are 62, 369 and 1,388 before and
62, 369 and 1,200 after; the test suite pins 62, 369 and 1,200 for R4 and 1,409 for R6's
first level.

## Removed mechanisms

Removed: 12,069. Of these, 564 have a steady-state group that carries no flux and 11,505
have a copy group whose twin test fails; none is unexplained.

- Zero-flux removals (564): the RE twin of each is in the new population, 564 of 564. R4
  loses 188. The findings predict 148 flip and 40 split children; the script does not tell
  flip from split children. R5 loses 190, all zero-flux (derived from the totals), 2 more
  than R4, and the other 186 zero-flux removals are R6's.
- Twin-only copy groups (11,505): the dead-end children whose inhibitor copy of a substrate or
  product binds a complex a form already holds. R3 loses all of its dead-end children, as a
  uni-uni mechanism has no new complex to give; R6 loses 20% of its level 1 (360 of 1,769) and 40% of its level 2
  (11,324 of 28,304).
- Rank: for the 188 removed R4 zero-flux mechanisms, the rank is at least 1 and fitted
  parameters minus the rank is at least the number of zero-flux groups in all 188; 0
  failures. Fitted parameters range over 8 to 9 and the rank over 6 to 8. The removed
  parameters are ones the rate equation cannot see.
- R1, R2 and ALLO lose nothing, and the allosteric populations are unchanged.

## Added mechanisms

None: 0 added, so 0 violate an emission rule. The extended flip pairs the plan expects appear
only at level 3, and these populations do not reach level 3.

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

175 mechanisms, 0 differences, 0 errors.

## Enumeration time

| Set | Old (s) | New (s) | Change |
|---|---|---|---|
| R1 | 0.78 | 0.82 | 6% |
| R2 | 0.06 | 0.07 | 8% |
| R3 | 0.01 | 0.00 | -94% |
| R4 | 0.46 | 0.48 | 5% |
| R5 | 1.17 | 1.31 | 12% |
| R6 | 8.15 | 5.38 | -34% |
| ALLO | 71.69 | 56.06 | -22% |

Each time is one run per checkout; there were no repeats. R6, the largest of R1–R6,
enumerates 34% faster; R5 is 12% slower and R1, R2 and R4 are 5% to 8% slower, on populations
that take between 0.06 and 1.3 s. ALLO is 22% faster.

## Conclusion

The filters remove 12,069 mechanisms, each named by a rule: 564 with a steady-state group
that carries no flux, whose rapid-equilibrium twin stays in the population, and 11,505 whose
dead-end copy adds no new complex. They add none and change no kept mechanism's derivation in
the 175-mechanism sample. Enumeration is 34% faster on R6 and up to 12% slower on the
sub-two-second sets.

The snapshot script's own derivation sample overlapped only 24 mechanisms between the runs,
because both runs shared one random stream and R3's population differs; the separate
intersection sample above replaced it.
