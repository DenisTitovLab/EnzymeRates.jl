# Merged seeds and the flip rule: regression check

Date: 2026-10-04. This checks the functional changes of sub-project C (design
`2026-10-03-merged-seeds-and-fused-steps-design.md`, Tasks 5–11 of
`../plans/2026-10-03-merged-seeds-and-fused-steps.md`) against 7857b6d, the opening wave's tip
(`2026-10-03-opening-wave-regression.md`). The tip is d7c64d4. Between them: bindings classified
by stoichiometry (19cdc12, bae52fb), the moves on fused and Theorell–Chance parents (d6c52e3),
dead-end sites where the metabolite is free (12acf7a), the flip rule (51e232b), the seed
predicates (7c26f63), the merged and Theorell–Chance seeds (27c829b, 0469ab6), and degenerate
seeds expanded without fitting (d7c64d4).

Every check passes. No compared mechanism without a fused or Theorell–Chance step changes its
rate string, fitted names or rank; the populations grown from today's seeds lose exactly the
flank flips of qualifying chains and their descendants; every new seed derives, fits its rank
and is non-degenerate; V and C agree with the numeric probe on every candidate the seed search
tests; and R6 levels 0–2 enumerate in 29% less time per mechanism. Ter-ter seeding with a
required inhibitor is estimated at about 15 minutes, over the ten minutes set for running it
and over the design's one-minute line for ter-ter seeding, which the base already exceeds.

## Method

Scratch scripts (not committed) follow the opening wave's record:

- `lib.jl`: the reactions, the populations, sig-string keys, and copies of the test file's rank
  oracle (`_testhelper_identifiable_rank`) and degeneracy probe (`_testhelper_degenerate`). A
  key holds `repr(_to_sig(s))` for the reaction and every step, group by group in stored order,
  with an allosteric mechanism's tags and multiplicity appended; it is finer than
  `_forward_sides` per group and does not depend on names. `oracle.jl` copies the mass-action
  oracle of `test/test_rate_eq_derivation.jl` (consistent point, full linear steady state).
- `snap.jl` enumerates the populations of one checkout, writes one key per mechanism, and runs
  `_assert_emission_rules` on every mechanism; `compare.jl` compares two snapshots per set and
  level and draws the rate-string samples.
- `flipcheck.jl` (base checkout) takes every mechanism found only on the base, finds each
  previous-level parent and move that emits it, and builds its RSR twin: the mechanism with
  both flanks of every steady-state-isomerization qualifying chain at rapid equilibrium. Its
  chain test is a copy of the tip's `_chain_flank_groups`.
- `derive.jl` derives a list of mechanisms on one checkout and writes fitted names, rank, probe
  verdict and the Reduced `rate_equation_string`; `cmpstrings.jl` compares two such lists. Where
  two strings differ it evaluates both, with a small expression evaluator (no compilation), at
  five consistent points built as the mass-action oracle builds them on the tip, mapping a
  base name absent on the tip to the reciprocal of its reversed constant (K_X_to_Y = 1/K_Y_to_X).
- `newseeds.jl` and `massaction_all.jl` check every merged and Theorell–Chance seed.
- `vc.jl` runs an instrumented copy of `_seed_variants` that logs every candidate its
  `admissible` tests, with the screen's verdict and its bottomless and flux tests, then
  evaluates the documented predicates (`_re_turnover_cycle`, `_has_vmax`,
  `_chemistry_equilibrates_both_sides`) and the probe on each candidate.
- `initseed.jl`, `terter_order.jl`, `terter_init_timing.jl`, `timing.jl`, `timing_gc.jl` and
  `timing_gcstats.jl` count and time `init_mechanisms` and `seed_mechanisms` and time R6.

Populations are those of B's record: `init_mechanisms` plus two levels of the flip
(`_expand_re_to_ss`), split (`_expand_split_kinetic_group`) and dead-end
(`_expand_add_dead_end_regulator`) moves on R1–R6, deduplicated across levels, with R4 also to
level 3; ALLO, the allosteric children of the R4 init mechanisms with their flip and split
children; and ALLO6, 150 allosteric children of the R6 seeds drawn with
`MersenneTwister(20260930)`, their dead-end, tag-relaxation, flip and split children, and the
tag relaxations of those that bind a copy. R1 is uni-uni, R2 uni-uni with one dead-end
inhibitor, R3 uni-uni with S and P as inhibitors, R4 the bi-bi reaction whose atoms admit
ping-pong (`bi_bi_pp_rxn`), R5 R4 with one dead-end inhibitor, R6 R4 with A, B, P and Q as
inhibitors. On the tip each population is grown twice: from the seeds (the init mechanisms
that hold an isomerization, the base's init list) and from all init mechanisms. Samples are
drawn with `MersenneTwister(20261003)` from keys sorted as strings.

The base ran as a detached worktree with the current `Manifest.toml`, removed afterwards. One
Julia process ran at a time, with `--heap-size-hint=2500M` unless a timing says otherwise.

## Populations

| Set | Levels | 7857b6d | Tip, from the seeds | Only on the base | Tip, all init mechanisms | Grown from the new seeds |
|---|---|---|---|---|---|---|
| R1 | 0, 1, 2 | 1, 2, 1 | 1, 0, 0 | 0, 2, 1 | 1, 0, 0 | — |
| R2 | 0, 1, 2 | 1, 3, 3 | 1, 1, 0 | 0, 2, 3 | 1, 1, 0 | — |
| R3 | 0, 1, 2 | 1, 2, 3 | 1, 0, 0 | 0, 2, 3 | 1, 0, 0 | — |
| R4 | 0, 1, 2, 3 | 62, 369, 1,200, 2,795 | 62, 349, 1,134, 2,677 | 0, 20, 66, 118 | 264, 1,018, 2,371, 4,437 | 202, 669, 1,237, 1,760 |
| R5 | 0, 1, 2 | 62, 719, 3,746 | 62, 699, 3,600 | 0, 20, 146 | 264, 2,488, 10,120 | 202, 1,789, 6,520 |
| R6 | 0, 1, 2 | 62, 1,649, 31,730 | 62, 1,629, 31,382 | 0, 20, 348 | 264, 6,738, 128,493 | 202, 5,109, 97,111 |
| ALLO | 1, 2 | 930, 5,471 | 930, 5,471 | 0 | 1,954, 8,676 | 1,024, 3,205 |
| ALLO6 | 1, 2, 3 | 150, 4,609, 14,980 | 150, 4,609, 14,980 | 0 | not built | — |

- Grown from the seeds, every tip population is a subset of the base's, level by level, and
  adds no mechanism. ALLO and ALLO6 are identical. The removals are the flip rule's (next
  section). Uni-uni loses both of its flips: each binding is a flank of the seed's chain.
- Grown from all init mechanisms, each tip population is the seeds' population plus a disjoint
  part grown from the new seeds, level by level. Every mechanism of that part holds a fused or
  Theorell–Chance step and no isomerization; no mechanism grown from the seeds holds either.
  The moves neither add nor remove an isomerization.
- R1–R3 gain nothing: uni-uni has no new seed.
- R6 level 2 was enumerated in full (26 s on the tip), not sampled.
- `_assert_emission_rules` holds on every mechanism of every set: 120,625 keys on the base and
  240,089 on the tip, counted per set, 0 violations.
- The test suite's pins stand: R4 62, 349, 1,134 and 264, 1,018, 2,371; R6 62, 1,629, 31,382 and
  264, 6,738, 128,493.

## The flip rule

R4 levels 0–3, grown from the seeds: the base holds 204 mechanisms the tip lacks, 20, 66 and 118
at levels 1, 2 and 3, and the tip holds none the base lacks.

- Every emission of every removed mechanism is a flip of a flank of a qualifying chain of a
  parent both checkouts hold, or any move of a removed parent: 166 flank flips (20, 56 and 90
  by level) and 38 descendants of removed mechanisms (0, 10 and 28). None is unexplained.
- Each of the 204 holds a qualifying chain with a steady-state isomerization and at least one
  steady-state flank, and its RSR twin is in the tip's population, 204 of 204.

The other sets, by the same test:

| Set | Only on the base | Flank flips | Descendants | RSR twin present |
|---|---|---|---|---|
| R1 | 3 | 2 | 1 | 3 |
| R2 | 5 | 4 | 1 | 5 |
| R3 | 5 | 2 | 3 | 3 |
| R5 | 166 | 156 | 10 | 166 |
| R6 | 368 | 356 | 12 | 366 |

Four RSR twins are absent, two in R3 and two in R6. Each of the four is a dead-end child of a
flank flip whose copy, with the flanks at rapid equilibrium, forms a redundant copy group, so
the copy rule rejects the twin on both checkouts. The R3 pair are uni-uni mechanisms at 5
fitted, rank 3 (their twins 4 fitted, rank 3; the seed 3 and 3). The R6 pair are children of a
degenerate ping-pong seed at 8 fitted, rank 5 (their twins 7 and 5), the seed's rank; the probe
flags both degenerate. All 13 uni-uni removals carry 1 or 2 phantoms (fitted 4 or 5, rank 3 or
4).

Phantoms (fitted minus rank, by the oracle) at R4 levels 0–2, every mechanism ranked on its
own checkout:

| Level | 7857b6d: mechanisms, phantoms, with one or more | Tip: mechanisms, phantoms, with one or more | Removed: mechanisms, phantoms |
|---|---|---|---|
| 0 | 62, 2, 2 | 62, 2, 2 | 0, 0 |
| 1 | 369, 35, 33 | 349, 13, 13 | 20, 22 |
| 2 | 1,200, 229, 195 | 1,134, 146, 129 | 66, 83 |
| All | 1,631, 266, 230 | 1,545, 161, 144 | 86, 105 |

Every removed mechanism carries at least one phantom. The 1,545 mechanisms on both checkouts
have the same fitted names, rank, probe verdict and rate string on each.

Level 3, a draw of 150 per stratum:

| Stratum | Mechanisms | Ranked | Phantoms | With one or more |
|---|---|---|---|---|
| Base, level 3 | 2,795 | 150 | 30 | 24 (16.0%) |
| Base, level 3, removed by the rule | 118 | 118 | 163 | 118 (100%) |
| Tip, level 3 grown from the seeds | 2,677 | 150 | 27 | 20 (13.3%) |
| Tip, level 3 grown from the new seeds | 1,760 | 150 | 0 | 0 |

Scaled, about 447 of the base's 2,795 level-3 mechanisms carry a phantom and about 357 of the
tip's 4,437 (8.0%). The new seeds' descendants carry none in draws of 150 at levels 1 and 2
either, and none of those 450 is degenerate by the probe.

## Rate strings

| Set | Compared | Identical | Different string, equal law |
|---|---|---|---|
| `MECHANISM_TEST_SPECS` | 43 | 41 | 2 |
| Fused and Theorell–Chance fixtures of `test_rate_eq_derivation.jl` | 24 | 18 | 6 |
| R1–R3, every mechanism on both | 4 | 4 | 0 |
| R4, every mechanism on both | 1,545 | 1,545 | 0 |
| R5, drawn from 4,361 on both | 400 | 400 | 0 |
| R6, drawn from 33,073 on both | 400 | 400 | 0 |
| ALLO, drawn from 6,401 on both | 100 | 100 | 0 |
| ALLO6, drawn from 19,739 on both | 100 | 100 | 0 |

Every enumerated mechanism compared also has the same fitted names, rank and probe verdict on
both checkouts. The eight different strings all belong to mechanisms with a fused step; each
keeps its law and changes only its parameterization:

- Uni-Uni and fixtures 1 and 4 (`E(S) <--> E + P`, a steady-state fused release now stored as
  the binding it reverses): the base fits `k_ES_to_E_P`; the tip fits `k_ES_to_E_S` and derives
  `k_ES_to_E_P` by Haldane. RE Uni-Uni and fixture 2: `k_EA_to_E_P` and `k_ES_to_E_P` become
  `k_E_P_to_EA` and `k_E_P_to_ES`.
- Fixture 7 (`E(P) ⇌ E + S`, a rapid-equilibrium fused binding): `K_E_S_to_EP` becomes its
  dissociation constant `K_EP_to_E_S`, the one name the map inverts.
- Fixtures 14 and 16 (bi-bi with steady-state fused releases from E(A, B)): `k_EAB_to_EP_Q` and
  `k_EAB_to_EQ_P` become `k_EP_Q_to_EAB` and `k_EQ_P_to_EAB`.

At the five consistent points the two strings of each agree to a relative 7e-76 (BigFloat): the
laws are equal. The 25th fixture, a fused and a plain binding of B in one group, builds only on
the tip; the base's constructor rejects the mixed group.

## The new seeds

`init_mechanisms` on both checkouts:

| Reaction | 7857b6d | Tip | Merged and Theorell–Chance |
|---|---|---|---|
| R1, R2, R3 (uni-uni) | 1 | 1 | 0 |
| R4, R5, R6 (`bi_bi_pp_rxn`) | 62 | 264 | 202 |
| `bi_bi_rxn` | 55 | 239 | 184 |
| `uni_bi_rxn` | 3 | 7 | 4 |
| `ter_ter_rxn` | 35,665 | 250,855 | 215,190 |

On every reaction the tip's list begins with the base's list in the base's order, and every
later mechanism holds no isomerization. R4, R5 and R6 give the same steps in the same order:
`init_mechanisms` ignores regulators. Bi-bi meets the design's estimate of about 264 exactly.

| Kind | R4 | `bi_bi_rxn` | `uni_bi_rxn` | Fitted = rank |
|---|---|---|---|---|
| Merged, two steady-state groups | 122 | 108 | 4 | 5 (uni-bi 4) |
| Merged, three steady-state groups | 56 | 56 | 0 | 6 |
| Theorell–Chance | 8, 8, 4 | 8, 8, 4 | 0 | 5, 6, 7 |
| Half Theorell–Chance (ping-pong) | 4 | 0 | 0 | 5 |

All 390 new seeds derive, fit exactly their rank, are non-degenerate by the probe and by
`_degenerate`, and satisfy both emission rules. The mass-action oracle on a draw of 30 (15 of
R4, 15 of `bi_bi_rxn`) gives a largest relative error of 6.6e-15, and on all 390 of 8.2e-15,
against the test's tolerance of 1e-8.

## V and C against the probe

| Reaction | Seeds | Candidates tested | Admissible | Lumping twins | Variants | V and C agree with the probe |
|---|---|---|---|---|---|---|
| R4 | 62 | 844 | 210 | 8 | 202 | 844 |
| R6 | 62 | 844 | 210 | 8 | 202 | 844 |
| `bi_bi_rxn` | 55 | 746 | 192 | 8 | 184 | 746 |
| `uni_bi_rxn` | 3 | 24 | 6 | 2 | 4 | 24 |

- Every candidate the search tests has no rapid-equilibrium turnover cycle, no bottomless
  segment and only flux-carrying steady-state groups, so the probe runs on all of them, and
  each is admissible exactly when the probe finds it non-degenerate. There is no disagreement.
- The screen's verdict equals the documented predicates' on every candidate, and the
  instrumented copy returns `_seed_variants`' variants on every seed. R6's log is R4's,
  candidate by candidate.
- R4's 634 rejections: C alone 142, V for substrates alone 120, V for products alone 120, V for
  substrates and C 114, V for products and C 114, both V 20, all three 4. `bi_bi_rxn`'s 554: 124,
  109, 109, 96, 96, 20 and 0.

## Timing

R6 levels 0–2, warm, three runs after a warm-up, two rounds alternating the checkouts:

| Checkout | Mechanisms | Round 1 (s) | Round 2 (s) | Best (s) | Per mechanism (ms) |
|---|---|---|---|---|---|
| 7857b6d | 33,441 | 8.90, 8.78, 8.83 | 8.85, 8.70, 8.79 | 8.70 | 0.260 |
| d7c64d4 | 135,495 | 25.24, 39.50, 42.21 | 25.01, 39.62, 42.23 | 25.01 | 0.185 |

Per mechanism the tip is 29% faster; in total it takes 2.9 times as long for 4.05 times the
mechanisms. The tip's slow second and third runs are garbage collection under the heap hint. A
third round that recorded collection time took 25.01, 39.57 and 42.87 s, of which 5.0, 19.4 and
22.6 s collecting. Without the hint the tip's three runs take 25.03, 24.90 and 24.72 s
(collection 4.2–4.6 s) and the base's 8.95, 8.85 and 8.83 s (1.3–1.5 s). The slowest hinted run,
42.87 s or 0.316 ms per mechanism, is 22% over the base's best. The limit is set on the best
run, as in B's record, and the tip's best is well inside it.

`init_mechanisms(ter_ter_rxn)`, one process each, the first result released before the warm
call: the base 16.8 s cold and 12.9 s warm; the tip 48.9 s cold and 43.2 s warm, with a 2.0 GB
peak. The same with and without the heap hint.

`seed_mechanisms` on R5 with I required: the base 350 seeds, 2.88 s on the first call and 0.04 to
0.12 s warm; the tip 1,470 seeds, 2.95 s and 0.10 to 0.13 s.

Ter-ter with a required dead-end inhibitor I, estimated rather than run. A copy of
`seed_mechanisms`' search that takes its first frontier as an argument ran from 500 init
mechanisms drawn with `MersenneTwister(20261003)`. The counts scale linearly, since removing
the I steps from a seed recovers its init mechanism, so no seed is reached from two:

| Checkout | Init mechanisms | Seeds from 500 | Time (s) | Scaled seeds | Scaled time |
|---|---|---|---|---|---|
| 7857b6d | 35,665 | 11,332 | 2.9 | about 808,000 | about 3.4 min |
| d7c64d4 | 250,855 | 10,220 | 1.8 | about 5,130,000 | about 15 min |

The tip's estimate exceeds the ten minutes set for running the search, so it was not run. The
design brings ter-ter seeding over one minute to Denis: `init_mechanisms` stays inside that
line, and seeding with a required inhibitor exceeds it on both checkouts. At the size of the
init mechanisms (Task 10 measured 534 MB live for 250,855), five million seeds need about
11 GB, more than this machine holds; the base's 808,000 need about 1.7 GB.

## Limitations

- Rate strings of R5 and R6 were compared on draws of 400 of 4,361 and 33,073 shared mechanisms,
  ALLO and ALLO6 on draws of 100; R1–R4 and the fixtures in full.
- Mechanisms grown from the new seeds exist only on the tip. All 390 seeds were derived and
  checked; their descendants only on draws (450 ranked at R4 levels 1–3; 100 of the 4,229 tip-only
  ALLO mechanisms, all of which derive, 26 with one phantom, against 49 of the 100 shared ALLO
  mechanisms drawn).
- A known MWC derivation defect (the five `@test_broken` gates of Task 9) gives an allosteric
  mechanism whose free enzyme shares a rapid-equilibrium segment with a residual form, among
  them K-type children of the ping-pong seeds, a law whose L term misses a factor. The defect
  predates C and both checkouts share it, so it is no regression; the 200 allosteric strings
  drawn match. The probe reads such laws as they are derived, and the V and C check runs on plain
  mechanisms only.
- ALLO6 was not built from all of the tip's init mechanisms: its draw of 150 would differ from
  the base's.
- Ter-ter seeding with a required inhibitor is an estimate on both checkouts.
- The test suite was not run for this record; the last full run, on d7c64d4, passed.

## Conclusion

C changes nothing it was not meant to change. Grown from today's seeds, the tip keeps every
mechanism but the flip rule's removals and adds none, and every kept mechanism compared derives
exactly as at 7857b6d: all 1,545 R4 mechanisms of levels 0–2 and the draws of R5, R6, ALLO and
ALLO6. The flip rule removes 204 R4 mechanisms to level 3, each a
flank flip or its descendant with its RSR twin kept, and adds none; at levels 0–2 it removes 105
of 266 phantoms. The 202 merged and Theorell–Chance seeds of R4 (184 of bi-bi, 4 of uni-bi) all
fit their rank and pass the mass-action oracle, and their descendants carry no phantom in the
draws. V and C agree with the probe on all 1,614 candidates of R4, bi-bi and uni-bi. R6 levels
0–2 enumerate 29% faster per mechanism. Ter-ter seeding with a required inhibitor would reach
about five million seeds in about 15 minutes, against 808,000 at 7857b6d; that cost goes to
Denis with the ter-ter base-tier concerns of Tasks 10 and 11.
