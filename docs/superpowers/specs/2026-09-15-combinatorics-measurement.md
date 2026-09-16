# Enumeration combinatorics measurement

Date: 2026-09-15
Branch: `refactor-re-to-ss-steps-to-only-keep-mechanisms-with-more-params`, commit `254ecb4` (the two commits before `ca05783` touch only docs)
Script: `docs/superpowers/specs/2026-09-15-combinatorics-measurement.jl`
(run with `julia --project docs/superpowers/specs/2026-09-15-combinatorics-measurement.jl`)

Reaction:

```julia
rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    allosteric_regulators: R(1)
    competitive_inhibitors: I
end
```

## Method notes and deviations from the brief

- The measurement writes the regulator as `R(1)`, the single-subunit case; a
  bare `R` inherits the default `allowed_catalytic_multiplicities: (1,)` and
  gives the same reaction.
- Per-move counts are **raw move output**. `expand_mechanisms` unions the seven
  moves and then applies `_filter_by_reg_type`, so its child count per parent is
  at most the row sums in the tables below.
- `seed_mechanisms(rxn, Set([:R]), Set([:I]))` returns 19,904 seeds and
  `_expand_split_kinetic_group` costs ~0.5 s per parent, so sweeping all seven
  moves over the whole seed set would take about three hours. Item 3 therefore
  reports scope `all` (all 19,904 seeds) for the five cheap moves and scope
  `sample` (every 40th seed in canonical printed order, 498 parents) for
  `_expand_split_kinetic_group` and `_expand_change_allo_state`. Item 4's
  depth-1 parent pool is built from the same 498 sampled seeds.
- Mechanisms are printed in `@allosteric_mechanism` style: `(step, step)` for a
  kinetic group, `⇌`/`<-->` for RE/SS, `E(A, B)` for forms, `::Inh` on a bound
  competitive inhibitor, the group's allosteric tag after `::`, and one
  `regulatory_site(...)` line per site. Each printed mechanism is preceded by its
  step, group, RE-segment and independent-parameter counts
  (`_re_segment_count`, `_independent_param_count`).
- Item 5's "RE segmentation" key is the partition of the mechanism's enzyme
  forms into RE-connected segments (from `_compute_re_groups`, the same
  partition `_re_segment_count` counts), rendered as a sorted tuple of sorted
  form-name tuples. This is the complete segmentation; the brief's step-string
  variant is a projection of it.
- Each closure is capped at 20,000 visited nodes **and** 120 s of wall clock; a
  truncated closure is marked `(cap hit)` or `(time hit)` in the table. The split
  closure of a largest seed does not finish inside that budget, so item 5 also
  probes the three seeds with the FEWEST steps, where both closures complete.
- The whole script runs in 1005 s (16.7 min); items 1-4 are identical across
  two runs apart from the elapsed-time lines, so the numbers are deterministic.

## Why four of the seven moves emit nothing here

Four moves have no work left to do on a fully-required seed set with one
allosteric regulator and one competitive inhibitor, which is why their rows are
zero across all 19,904 seeds:

- `_expand_add_dead_end_regulator`: the only dead-end-eligible regulator is `I`,
  which every seed already binds (it is required); `R` is excluded on an
  `AllostericMechanism` because it belongs on a regulatory site.
- `_expand_to_allosteric`: every seed is already an `AllostericMechanism`, and
  the move is a no-op on that type.
- `_expand_add_allosteric_regulator`: `R` is the only allosteric regulator and
  every seed already binds it.
- `_expand_merge_regulatory_sites`: merging needs two sites; every seed has one.

So at depth 1 and depth 2 the branching factor is carried entirely by
`_expand_re_to_ss`, `_expand_split_kinetic_group` and
`_expand_change_allo_state`.

## Open question for Denis

The split closure of an 18-step seed is unresolved: it exceeds 62 nodes and
grows at roughly 0.5 nodes/s (`_expand_split_kinetic_group` costs ~2 s on a
parent that size), so 120 s buys only 62 nodes. At the other end, a 6-step seed
has a split closure of 1 — every kinetic group is a single step, so no group
admits a context bipartition. Nothing in between was measured. If the
combinatorics page needs the real figure for a rich seed, it needs a dedicated
run of tens of minutes on that one closure.

---

## 1. `init_mechanisms`

| quantity | value |
|---|---|
| catalytic topologies | 9 |
| dead-end (competition) patterns | 7 |
| `init_mechanisms(rxn)` | 55 |
| seeds per topology: min | 5 |
| seeds per topology: median | 7.0 |
| seeds per topology: max | 7 |
| sum over topologies == `init_mechanisms` | true |

Topology with the FEWEST init mechanisms (5):
      E + A ⇌ E(A)
      E + P ⇌ E(P)
      E(A) + B ⇌ E(A, B)
      E(P) + Q ⇌ E(P, Q)
      E(A, B) <--> E(P, Q)

Topology with the MOST init mechanisms (7):
      E + B ⇌ E(B)
      E + P ⇌ E(P)
      E + Q ⇌ E(Q)
      E(B) + A ⇌ E(A, B)
      E(P) + Q ⇌ E(P, Q)
      E(Q) + P ⇌ E(P, Q)
      E(A, B) <--> E(P, Q)

4.0 s elapsed

## 2. `seed_mechanisms(rxn, Set([:R]), Set([:I]))`

| quantity | value |
|---|---|
| seeds | 19904 |
| allosteric seeds | 19904 |
| plain `Mechanism` seeds | 0 |
| steps per seed: min / median / max | 6 / 12.0 / 18 |

144.9 s elapsed

## 3. Children per parent, over the seed set of item 2

Counts are raw per-move output; `expand_mechanisms` applies `_filter_by_reg_type`
to the union, so its child count is ≤ the row sums below. Scope `all` = all 19904
seeds; scope `sample` = every 40-th seed in canonical order (498 parents).

| move | scope | parents | parents with ≥1 child | children | min | median | max |
|---|---|---|---|---|---|---|---|
| `_expand_re_to_ss` | all | 19904 | 19904 | 79616 | 4 | 4.0 | 4 |
| `_expand_split_kinetic_group` | sample | 498 | 493 | 1646 | 0 | 3.0 | 7 |
| `_expand_add_dead_end_regulator` | all | 19904 | 0 | 0 | 0 | 0.0 | 0 |
| `_expand_to_allosteric` | all | 19904 | 0 | 0 | 0 | 0.0 | 0 |
| `_expand_add_allosteric_regulator` | all | 19904 | 0 | 0 | 0 | 0.0 | 0 |
| `_expand_change_allo_state` | sample | 498 | 498 | 3020 | 6 | 6.0 | 7 |
| `_expand_merge_regulatory_sites` | all | 19904 | 0 | 0 | 0 | 0.0 | 0 |

**_expand_re_to_ss** — parent at the minimum (4 children):
    AllostericMechanism, 11 steps, 6 groups, 1 RE segments, 8 indep params
      (E + A ⇌ E(A), E(B) + A ⇌ E(A, B)) :: EqualAI
      (E + B ⇌ E(B), E(A) + B ⇌ E(A, B)) :: EqualAI
      (E + I ⇌ E(I::Inh), E(P) + I ⇌ E(I::Inh, P)) :: EqualAI
      (E + P ⇌ E(P), E(I::Inh) + P ⇌ E(I::Inh, P)) :: EqualAI
      (E(A) + Q ⇌ E(A, Q), E(P) + Q ⇌ E(P, Q)) :: EqualAI
      E(A, B) <--> E(P, Q) :: OnlyA
      catalytic_multiplicity: 1
      regulatory_site(multiplicity = 1): R::OnlyA
**_expand_re_to_ss** — parent at the maximum (4 children):
    AllostericMechanism, 11 steps, 6 groups, 1 RE segments, 8 indep params
      (E + A ⇌ E(A), E(B) + A ⇌ E(A, B)) :: EqualAI
      (E + B ⇌ E(B), E(A) + B ⇌ E(A, B)) :: EqualAI
      (E + I ⇌ E(I::Inh), E(P) + I ⇌ E(I::Inh, P)) :: EqualAI
      (E + P ⇌ E(P), E(I::Inh) + P ⇌ E(I::Inh, P)) :: EqualAI
      (E(A) + Q ⇌ E(A, Q), E(P) + Q ⇌ E(P, Q)) :: EqualAI
      E(A, B) <--> E(P, Q) :: OnlyA
      catalytic_multiplicity: 1
      regulatory_site(multiplicity = 1): R::OnlyA

**_expand_split_kinetic_group** — parent at the minimum (0 children):
    AllostericMechanism, 7 steps, 6 groups, 1 RE segments, 9 indep params
      E + A ⇌ E(A) :: OnlyA
      E + I ⇌ E(I::Inh) :: EqualAI
      E + Q ⇌ E(Q) :: EqualAI
      E(A) + B ⇌ E(A, B) :: OnlyA
      (E(A) + P ⇌ E(A, P), E(Q) + P ⇌ E(P, Q)) :: EqualAI
      E(A, B) <--> E(P, Q) :: OnlyA
      catalytic_multiplicity: 1
      regulatory_site(multiplicity = 1): R::OnlyA
**_expand_split_kinetic_group** — parent at the maximum (7 children):
    AllostericMechanism, 18 steps, 6 groups, 1 RE segments, 8 indep params
      (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(I::Inh) + A ⇌ E(A, I::Inh), E(P) + A ⇌ E(A, P)) :: OnlyA
      (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q)) :: OnlyA
      (E + I ⇌ E(I::Inh), E(A) + I ⇌ E(A, I::Inh), E(Q) + I ⇌ E(I::Inh, Q)) :: OnlyA
      (E + P ⇌ E(P), E(A) + P ⇌ E(A, P), E(Q) + P ⇌ E(P, Q)) :: EqualAI
      (E + Q ⇌ E(Q), E(B) + Q ⇌ E(B, Q), E(I::Inh) + Q ⇌ E(I::Inh, Q), E(P) + Q ⇌ E(P, Q)) :: EqualAI
      E(A, B) <--> E(P, Q) :: OnlyA
      catalytic_multiplicity: 1
      regulatory_site(multiplicity = 1): R::OnlyA

**_expand_add_dead_end_regulator** emits no children on any parent in scope.

**_expand_to_allosteric** emits no children on any parent in scope.

**_expand_add_allosteric_regulator** emits no children on any parent in scope.

**_expand_change_allo_state** — parent at the minimum (6 children):
    AllostericMechanism, 10 steps, 6 groups, 1 RE segments, 8 indep params
      (E + A ⇌ E(A), E(B) + A ⇌ E(A, B)) :: EqualAI
      (E + B ⇌ E(B), E(A) + B ⇌ E(A, B)) :: EqualAI
      (E + I ⇌ E(I::Inh), E(P) + I ⇌ E(I::Inh, P)) :: OnlyA
      (E + P ⇌ E(P), E(I::Inh) + P ⇌ E(I::Inh, P)) :: EqualAI
      E(A, B) <--> E(P, Q) :: OnlyA
      E(P) + Q ⇌ E(P, Q) :: EqualAI
      catalytic_multiplicity: 1
      regulatory_site(multiplicity = 1): R::OnlyA
**_expand_change_allo_state** — parent at the maximum (7 children):
    AllostericMechanism, 11 steps, 6 groups, 1 RE segments, 8 indep params
      (E + A ⇌ E(A), E(B) + A ⇌ E(A, B)) :: EqualAI
      (E + B ⇌ E(B), E(A) + B ⇌ E(A, B)) :: EqualAI
      (E + I ⇌ E(I::Inh), E(P) + I ⇌ E(I::Inh, P)) :: EqualAI
      (E + P ⇌ E(P), E(I::Inh) + P ⇌ E(I::Inh, P)) :: EqualAI
      (E(A) + Q ⇌ E(A, Q), E(P) + Q ⇌ E(P, Q)) :: EqualAI
      E(A, B) <--> E(P, Q) :: OnlyA
      catalytic_multiplicity: 1
      regulatory_site(multiplicity = 1): R::OnlyA

**_expand_merge_regulatory_sites** emits no children on any parent in scope.

535.3 s elapsed

## 4. Children per parent, one level down

depth-1 children (all moves over the 498 sampled seeds, deduped): 6574
parents: every 33-th in canonical order — 200 mechanisms

| move | parents with ≥1 child | children | min | median | max |
|---|---|---|---|---|---|
| `_expand_re_to_ss` | 200 | 767 | 3 | 4.0 | 8 |
| `_expand_split_kinetic_group` | 197 | 654 | 0 | 3.0 | 8 |
| `_expand_add_dead_end_regulator` | 0 | 0 | 0 | 0.0 | 0 |
| `_expand_to_allosteric` | 0 | 0 | 0 | 0.0 | 0 |
| `_expand_add_allosteric_regulator` | 0 | 0 | 0 | 0.0 | 0 |
| `_expand_change_allo_state` | 200 | 1171 | 5 | 6.0 | 9 |
| `_expand_merge_regulatory_sites` | 0 | 0 | 0 | 0.0 | 0 |

**_expand_re_to_ss** — parent at the minimum (3 children):
    AllostericMechanism, 11 steps, 6 groups, 2 RE segments, 9 indep params
      (E + A <--> E(A), E(B) + A <--> E(A, B)) :: EqualAI
      (E + B ⇌ E(B), E(A) + B ⇌ E(A, B)) :: EqualAI
      (E + I ⇌ E(I::Inh), E(P) + I ⇌ E(I::Inh, P)) :: EqualAI
      (E + P ⇌ E(P), E(I::Inh) + P ⇌ E(I::Inh, P)) :: EqualAI
      (E(A) + Q ⇌ E(A, Q), E(P) + Q ⇌ E(P, Q)) :: EqualAI
      E(A, B) <--> E(P, Q) :: OnlyA
      catalytic_multiplicity: 1
      regulatory_site(multiplicity = 1): R::OnlyA
**_expand_re_to_ss** — parent at the maximum (8 children):
    AllostericMechanism, 16 steps, 8 groups, 1 RE segments, 9 indep params
      (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(Q) + A ⇌ E(A, Q)) :: EqualAI
      (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(I::Inh) + B ⇌ E(B, I::Inh)) :: EqualAI
      (E + I ⇌ E(I::Inh), E(B) + I ⇌ E(B, I::Inh), E(P) + I ⇌ E(I::Inh, P)) :: EqualAI
      (E + P ⇌ E(P), E(I::Inh) + P ⇌ E(I::Inh, P)) :: EqualAI
      (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q)) :: EqualAI
      E(A, B) <--> E(P, Q) :: OnlyA
      E(P) + Q ⇌ E(P, Q) :: EqualAI
      E(Q) + P ⇌ E(P, Q) :: EqualAI
      catalytic_multiplicity: 1
      regulatory_site(multiplicity = 1): R::OnlyA

**_expand_split_kinetic_group** — parent at the minimum (0 children):
    AllostericMechanism, 7 steps, 7 groups, 1 RE segments, 9 indep params
      E + A ⇌ E(A) :: EqualAI
      E + I ⇌ E(I::Inh) :: EqualAI
      E + P ⇌ E(P) :: EqualAI
      E(A) + B ⇌ E(A, B) :: OnlyA
      E(A) + Q ⇌ E(A, Q) :: EqualAI
      E(A, B) <--> E(P, Q) :: OnlyA
      E(P) + Q ⇌ E(P, Q) :: EqualAI
      catalytic_multiplicity: 1
      regulatory_site(multiplicity = 1): R::OnlyA
**_expand_split_kinetic_group** — parent at the maximum (8 children):
    AllostericMechanism, 15 steps, 6 groups, 2 RE segments, 9 indep params
      (E + A <--> E(A), E(B) + A <--> E(A, B), E(I::Inh) + A <--> E(A, I::Inh), E(Q) + A <--> E(A, Q)) :: OnlyA
      (E + B ⇌ E(B), E(A) + B ⇌ E(A, B)) :: EqualAI
      (E + I ⇌ E(I::Inh), E(A) + I ⇌ E(A, I::Inh), E(Q) + I ⇌ E(I::Inh, Q)) :: EqualAI
      (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q), E(I::Inh) + Q ⇌ E(I::Inh, Q)) :: EqualAI
      E(A, B) <--> E(P, Q) :: OnlyA
      (E(B) + P ⇌ E(B, P), E(Q) + P ⇌ E(P, Q)) :: EqualAI
      catalytic_multiplicity: 1
      regulatory_site(multiplicity = 1): R::OnlyA

**_expand_add_dead_end_regulator** emits no children on any parent in scope.

**_expand_to_allosteric** emits no children on any parent in scope.

**_expand_add_allosteric_regulator** emits no children on any parent in scope.

**_expand_change_allo_state** — parent at the minimum (5 children):
    AllostericMechanism, 11 steps, 6 groups, 1 RE segments, 9 indep params
      (E + A ⇌ E(A), E(B) + A ⇌ E(A, B)) :: EqualAI
      (E + B ⇌ E(B), E(A) + B ⇌ E(A, B)) :: EqualAI
      (E + I ⇌ E(I::Inh), E(Q) + I ⇌ E(I::Inh, Q)) :: OnlyA
      (E + Q ⇌ E(Q), E(I::Inh) + Q ⇌ E(I::Inh, Q)) :: OnlyA
      E(A, B) <--> E(P, Q) :: OnlyA
      (E(B) + P ⇌ E(B, P), E(Q) + P ⇌ E(P, Q)) :: EqualAI
      catalytic_multiplicity: 1
      regulatory_site(multiplicity = 1): R::NonequalAI
**_expand_change_allo_state** — parent at the maximum (9 children):
    AllostericMechanism, 16 steps, 8 groups, 1 RE segments, 9 indep params
      (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(Q) + A ⇌ E(A, Q)) :: EqualAI
      (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(I::Inh) + B ⇌ E(B, I::Inh)) :: EqualAI
      (E + I ⇌ E(I::Inh), E(B) + I ⇌ E(B, I::Inh), E(P) + I ⇌ E(I::Inh, P)) :: EqualAI
      (E + P ⇌ E(P), E(I::Inh) + P ⇌ E(I::Inh, P)) :: EqualAI
      (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q)) :: EqualAI
      E(A, B) <--> E(P, Q) :: OnlyA
      E(P) + Q ⇌ E(P, Q) :: EqualAI
      E(Q) + P ⇌ E(P, Q) :: EqualAI
      catalytic_multiplicity: 1
      regulatory_site(multiplicity = 1): R::OnlyA

**_expand_merge_regulatory_sites** emits no children on any parent in scope.

640.6 s elapsed

## 5. Flip and split closures

Closure caps: 20000 visited nodes or 120.0 s per closure.

| seed | steps | groups | RE segments | indep params | flip closure | distinct RE segmentations | split closure |
|---|---|---|---|---|---|---|---|
| most-1 | 18 | 6 | 1 | 8 | 16 | 16 | 62 (time hit) |
| most-2 | 18 | 6 | 1 | 8 | 16 | 16 | 61 (time hit) |
| most-3 | 18 | 6 | 1 | 8 | 16 | 16 | 59 (time hit) |
| fewest-1 | 6 | 6 | 1 | 8 | 16 | 16 | 1 |
| fewest-2 | 6 | 6 | 1 | 8 | 16 | 16 | 1 |
| fewest-3 | 6 | 6 | 1 | 8 | 16 | 16 | 1 |

seed most-1:
    AllostericMechanism, 18 steps, 6 groups, 1 RE segments, 8 indep params
      (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(I::Inh) + A ⇌ E(A, I::Inh), E(P) + A ⇌ E(A, P)) :: EqualAI
      (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q)) :: EqualAI
      (E + I ⇌ E(I::Inh), E(A) + I ⇌ E(A, I::Inh), E(P) + I ⇌ E(I::Inh, P)) :: EqualAI
      (E + P ⇌ E(P), E(A) + P ⇌ E(A, P), E(I::Inh) + P ⇌ E(I::Inh, P), E(Q) + P ⇌ E(P, Q)) :: EqualAI
      (E + Q ⇌ E(Q), E(B) + Q ⇌ E(B, Q), E(P) + Q ⇌ E(P, Q)) :: EqualAI
      E(A, B) <--> E(P, Q) :: OnlyA
      catalytic_multiplicity: 1
      regulatory_site(multiplicity = 1): R::OnlyA

seed most-2:
    AllostericMechanism, 18 steps, 6 groups, 1 RE segments, 8 indep params
      (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(I::Inh) + A ⇌ E(A, I::Inh), E(P) + A ⇌ E(A, P)) :: EqualAI
      (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q)) :: EqualAI
      (E + I ⇌ E(I::Inh), E(A) + I ⇌ E(A, I::Inh), E(P) + I ⇌ E(I::Inh, P)) :: EqualAI
      (E + P ⇌ E(P), E(A) + P ⇌ E(A, P), E(I::Inh) + P ⇌ E(I::Inh, P), E(Q) + P ⇌ E(P, Q)) :: EqualAI
      (E + Q ⇌ E(Q), E(B) + Q ⇌ E(B, Q), E(P) + Q ⇌ E(P, Q)) :: EqualAI
      E(A, B) <--> E(P, Q) :: OnlyA
      catalytic_multiplicity: 1
      regulatory_site(multiplicity = 1): R::OnlyI

seed most-3:
    AllostericMechanism, 18 steps, 6 groups, 1 RE segments, 8 indep params
      (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(I::Inh) + A ⇌ E(A, I::Inh), E(P) + A ⇌ E(A, P)) :: EqualAI
      (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q)) :: EqualAI
      (E + I ⇌ E(I::Inh), E(A) + I ⇌ E(A, I::Inh), E(P) + I ⇌ E(I::Inh, P)) :: EqualAI
      (E + P ⇌ E(P), E(A) + P ⇌ E(A, P), E(I::Inh) + P ⇌ E(I::Inh, P), E(Q) + P ⇌ E(P, Q)) :: EqualAI
      (E + Q ⇌ E(Q), E(B) + Q ⇌ E(B, Q), E(P) + Q ⇌ E(P, Q)) :: OnlyA
      E(A, B) <--> E(P, Q) :: OnlyA
      catalytic_multiplicity: 1
      regulatory_site(multiplicity = 1): R::OnlyA

seed fewest-1:
    AllostericMechanism, 6 steps, 6 groups, 1 RE segments, 8 indep params
      E + B ⇌ E(B) :: OnlyA
      E + I ⇌ E(I::Inh) :: OnlyA
      E + Q ⇌ E(Q) :: OnlyA
      E(A, B) <--> E(P, Q) :: OnlyA
      E(B) + A ⇌ E(A, B) :: OnlyA
      E(Q) + P ⇌ E(P, Q) :: OnlyA
      catalytic_multiplicity: 1
      regulatory_site(multiplicity = 1): R::OnlyI

seed fewest-2:
    AllostericMechanism, 6 steps, 6 groups, 1 RE segments, 8 indep params
      E + B ⇌ E(B) :: OnlyA
      E + I ⇌ E(I::Inh) :: OnlyA
      E + Q ⇌ E(Q) :: OnlyA
      E(A, B) <--> E(P, Q) :: OnlyA
      E(B) + A ⇌ E(A, B) :: OnlyA
      E(Q) + P ⇌ E(P, Q) :: OnlyA
      catalytic_multiplicity: 1
      regulatory_site(multiplicity = 1): R::OnlyA

seed fewest-3:
    AllostericMechanism, 6 steps, 6 groups, 1 RE segments, 8 indep params
      E + B ⇌ E(B) :: OnlyA
      E + I ⇌ E(I::Inh) :: OnlyA
      E + Q ⇌ E(Q) :: OnlyA
      E(A, B) <--> E(P, Q) :: OnlyA
      E(B) + A ⇌ E(A, B) :: OnlyA
      E(Q) + P ⇌ E(P, Q) :: EqualAI
      catalytic_multiplicity: 1
      regulatory_site(multiplicity = 1): R::OnlyI

1004.5 s elapsed

DONE
