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

For each added mechanism it runs `_assert_emission_rules`. For a random sample of the
populations (100 per reaction, 50 for ALLO, fixed seed) it records `fitted_params` and the
Reduced `rate_equation_string`, and compares them where both runs sampled the same mechanism.
It also records the enumeration time of each population. The compare mode ran on the new
checkout, so it uses the new predicates; the snapshot mode uses only functions both
checkouts have.

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

Every new population is a subset of the old. R4's levels, 62, 369 and 1,388 before and 62, 369
and 1,200 after, are the counts the test suite pins.

## Removed mechanisms

Removed: 12,069. Of these, 564 have a steady-state group that carries no flux and 11,505
have a copy group whose twin test fails; none is unexplained.

- Zero-flux removals (564): the RE twin of each is in the new population, 564 of 564. R4
  loses 188 (the 148 zero-flux flip children and the 40 zero-flux split children that the
  test's comment names), R5 loses 190 with the same structures, and the rest are R6's.
- Twin-only copy groups (11,505): the dead-end children whose inhibitor copy of a substrate or
  product binds a complex a form already holds. R3 loses all of its dead-end children, as a
  uni-uni mechanism has no new complex to give; R6 loses the bulk of its levels 1 and 2.
- Rank: for the 188 removed R4 zero-flux mechanisms, fitted parameters minus the numerical
  rank is at least the number of zero-flux groups in all 188; 0 failures. The removed
  parameters are ones the rate equation cannot see.
- R1, R2 and ALLO lose nothing, and the allosteric populations are unchanged.

## Added mechanisms

None: 0 added, so 0 violate an emission rule. The extended flip pairs the plan expects appear
only at level 3, which these populations do not reach.

## Derivation of kept mechanisms

The two runs drew their samples from different populations, so only 24 sampled mechanisms
occur in both. For those 24, the fitted names and the Reduced string are identical in both
runs; 0 differ. The comparison is small, but the filters change which mechanisms exist, not
how one derives, and the kept mechanisms are unchanged as keys.

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

R6, the largest population, enumerates 34% faster; no set is more than 12% slower, and those
differences are within run-to-run variation on populations that take about a second.

## Conclusion

The filters remove 12,069 mechanisms, each named by a rule: 564 with a steady-state group
that carries no flux, whose rapid-equilibrium twin stays in the population, and 11,505 whose
dead-end copy adds no new complex. They add none, change no kept mechanism's derivation in
the sample, and cost no enumeration time.

The script's first compare run stopped at the read of the sampled strings, because a Reduced
string spans several lines; the reader was fixed to join them and the comparison rerun on the
same snapshots.
