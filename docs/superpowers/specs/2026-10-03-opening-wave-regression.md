# Opening wave: regression check

Date: 2026-10-03. This checks the opening wave of sub-project C (section 0 of
`2026-10-03-merged-seeds-and-fused-steps-design.md`, Tasks 1–4 of
`../plans/2026-10-03-merged-seeds-and-fused-steps.md`) against b340822, B's tip. Tasks 1–3
simplify B's code and must change no output. Task 4 deletes the split's unit pre-filter
(`redundant_part`), which can only add split children, each satisfying both emission rules.

## Method

Scratch scripts (not committed) follow B's record, `2026-09-30-exact-filters-regression.md`:

- `wave_regress.jl snap` enumerates the populations of one checkout, writes one key per
  mechanism, and runs `_assert_emission_rules` on every mechanism. A key holds the sig strings
  (`repr(_to_sig(s))`) of the reaction and of every step, group by group in stored order, with an
  allosteric mechanism's tags and multiplicity appended, so it does not depend on names.
- `wave_regress.jl compare` compares two snapshots per set and level. For a mechanism found only
  on the tip it would re-enumerate the set, list the moves of previous-level parents that emit
  it, run `_assert_emission_rules`, and rank it (the test suite's
  `_testhelper_identifiable_rank`).
- `split_diff.jl` runs the split with and without the pre-filter on every mechanism of a
  population as a parent and compares the child sets. Its split without the pre-filter is the
  checkout's own `_expand_split_kinetic_group` with the `redundant_part` line removed. With R6
  level-2 parents, this covers the split children at depth 3.
- `timing.jl` times R6 levels 0–2 warm: three runs after a warm-up, per checkout.

Populations, as in B's record: `init_mechanisms` plus two levels of the flip, split and dead-end
moves on R1–R6, deduplicated across levels (R6 level 2 in full); ALLO, the allosteric children of
the R4 seeds with their flip and split children; ALLO6, 150 allosteric children of the R6 seeds
drawn with `MersenneTwister(20260930)` from the children sorted by key, their dead-end,
tag-relaxation, flip and split children, and the tag relaxations of those that bind a copy.

Three snapshots ran, one Julia process at a time: b340822 (a detached worktree with the current
`Manifest.toml`, removed afterwards), 689f44d (Tasks 1–3), and the tree of Task 4.

## Populations

| Set | Level | b340822 | 689f44d | Task 4 | Only b340822 | Only Task 4 |
|---|---|---|---|---|---|---|
| R1 | 0, 1, 2 | 1, 2, 1 | 1, 2, 1 | 1, 2, 1 | 0 | 0 |
| R2 | 0, 1, 2 | 1, 3, 3 | 1, 3, 3 | 1, 3, 3 | 0 | 0 |
| R3 | 0, 1, 2 | 1, 2, 3 | 1, 2, 3 | 1, 2, 3 | 0 | 0 |
| R4 | 0, 1, 2 | 62, 369, 1200 | 62, 369, 1200 | 62, 369, 1200 | 0 | 0 |
| R5 | 0, 1, 2 | 62, 719, 3746 | 62, 719, 3746 | 62, 719, 3746 | 0 | 0 |
| R6 | 0, 1, 2 | 62, 1649, 31730 | 62, 1649, 31730 | 62, 1649, 31730 | 0 | 0 |
| ALLO | 1, 2 | 930, 5471 | 930, 5471 | 930, 5471 | 0 | 0 |
| ALLO6 | 1, 2, 3 | 150, 4609, 14980 | 150, 4609, 14980 | 150, 4609, 14980 | 0 | 0 |

The three key files are byte-identical, 65,756 lines in the same order. No mechanism is only on
the tip, so none is ranked. `_assert_emission_rules` holds on all 65,756 mechanisms of each
snapshot. The counts equal B's record, and the test suite's pins (62, 369 and 1,200 for R4; 62,
1,649 and 31,730 for R6) stand.

## The split with and without the pre-filter

| Set | Parents | Levels | Differ |
|---|---|---|---|
| R1–R3 | 17 | 0–2 | 0 |
| R4 | 1,631 | 0–2 | 0 |
| R5 | 4,527 | 0–2 | 0 |
| R6 | 33,441 | 0–2 | 0 |
| ALLO | 6,401 | 1–2 | 0 |
| ALLO6 | 4,759 | 1–2 | 0 |

None of these 50,776 parents has a split child that the pre-filter changes, down to depth 3.
The runs used 689f44d, except ALLO6's parents 2,001 to 4,759, which used b340822 (identical
populations). ALLO6's 14,980 level-3 parents were not run: the allosteric split costs about
0.24 s per parent there. The case the deletion admits needs a copy-group unit that fails its
gain test alone and a gaining superset whose other split breaks the part's gauge; no parent in
these populations has one, so no unit test reproduces it.

The corrected spec sentence holds: H1 with its two A groups merged into one has a redundant copy
group (6 fitted, rank 5), and its split children include H1 (7 fitted, rank 7, no redundant
group).

## Timing

R6 levels 0–2 (33,441 mechanisms), warm, best of three, two rounds alternating the checkouts:

| Checkout | Round 1 (s) | Round 2 (s) | Best (s) | Per mechanism (ms) |
|---|---|---|---|---|
| b340822 | 8.74 | 8.96 | 8.74 | 0.261 |
| Task 4 | 8.77 | 9.10 | 8.77 | 0.262 |

The tip is 0.3% slower than b340822, inside the spread between rounds. One run per snapshot
(first-use compilation in R1):

| Set | b340822 (s) | 689f44d (s) | Task 4 (s) |
|---|---|---|---|
| R1 | 2.79 | 2.68 | 2.82 |
| R2 | 0.05 | 0.05 | 0.06 |
| R3 | 0.38 | 0.25 | 0.27 |
| R4 | 0.45 | 0.56 | 0.55 |
| R5 | 1.22 | 1.22 | 1.27 |
| R6 | 9.18 | 9.10 | 9.11 |
| ALLO | 53.87 | 53.76 | 54.59 |
| ALLO6 | 93.77 | 92.62 | 93.23 |

The split alone, summed over the differential's parents, with and without the pre-filter: R5
3.93 s and 3.95 s; R6 281.4 s and 268.5 s; ALLO 422.7 s and 427.2 s; ALLO6 1,182.0 s and
1,171.3 s. Dropping the pre-filter saves the child it built per all-twin part and costs the
larger sets the gain test then tries; R6's split is 5% faster without it, the allosteric sets
within 1%.

## Conclusion

The opening wave changes no mechanism in R1–R6, ALLO or ALLO6 against b340822, and every
mechanism satisfies both emission rules. The split's pre-filter removed no rule-conforming child
in these populations or at depth 3, so its deletion is a correction of B's argument, not of its
output. R6 levels 0–2 enumerate in the same time as at b340822.
