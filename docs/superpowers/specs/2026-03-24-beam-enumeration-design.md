# Beam-Based Mechanism Enumeration Design

## Problem

The current `enumerate_mechanisms` pipeline materializes all mechanism variants
eagerly, causing OOM for reactions like bi-bi + 2 regulators. It needs to handle
ter-ter + 3 regulators. Additionally, the future `identify_rate_equation` requires
level-by-level iteration grouped by `param_count` for beam search with
cross-validation.

## Solution Overview

Replace `enumerate_mechanisms` with a beam search that expands mechanisms
level-by-level (ascending `param_count`). Only one level is materialized at a
time. The current eager pipeline is preserved as `old_enumerate_mechanisms` for
correctness verification.

## File and Naming Changes

- `enumerate_mechanisms` → `old_enumerate_mechanisms` (stays in
  `mechanism_enumeration.jl`)
- `test/test_mechanism_enumeration.jl` →
  `test/old_test_mechanism_enumeration.jl`
- New file: `src/beam_enumeration.jl` — public `enumerate_mechanisms`,
  `expand_mechanisms_same_param_count`, `expand_mechanisms_by_one_param`,
  `expand_mechanisms_by_two_params`
- New file: `test/test_beam_enumeration.jl`
- `_catalytic_topologies` stays shared (used by both old and new)

## Core Data Flow

```
catalytic_topologies (all, grouped by param_count)
    │
    ▼
Level N_min: catalytic seeds with smallest param_count
    │ expand_mechanisms_same_param_count → +0 variants (to fixed point)
    │ deduplicate
    │ expand_mechanisms_by_one_param → candidates at +1
    │ expand_mechanisms_by_two_params → candidates at +2 → cached[N_min+2]
    │
    ▼
Level N_min+1: expanded(+1) + catalytic seeds at this count + cached[N_min+1]
    │ expand_mechanisms_same_param_count → +0 variants (to fixed point)
    │ deduplicate
    │ expand_mechanisms_by_one_param → candidates at +1
    │ expand_mechanisms_by_two_params → candidates at +2 → cached[N_min+3]
    │
    ▼
  ... continues until max_param_count or no new mechanisms ...
```

- Catalytic topologies are computed eagerly upfront (small set even for
  ter-ter). Each topology has exactly 1 SS step and all other steps RE,
  giving its minimum param_count. Different topologies may have different
  minimum param_counts, so they enter the beam at their respective levels.
- Only one level is materialized at a time — previous level is discarded.
- A cache dictionary maps future `param_count` values to vectors of pending
  specs. When expanding level N, +2 candidates are inserted at cache[N+2].
  When starting level M, specs from cache[M] are merged with the +1
  expansion results from level M-1 and any catalytic seeds at M.
- The three expansion functions are pure:
  `Vector{AbstractMechanismSpec} → Vector{AbstractMechanismSpec}`.
- Deduplication happens within each level via concentration fingerprints +
  constraint descriptors (same as current `_deduplicate`). Within-level
  dedup is sufficient because two expansion paths that arrive at the same
  mechanism always arrive at the same param_count (the moves are additive
  and param_count is a property of the mechanism structure, not the path).
- Termination: expansion stops when all three expansion functions plus cache
  plus catalytic seeds produce an empty set after dedup, or when
  `param_count` exceeds `max_param_count`.
- For `identify_rate_equation` (future): the loop wraps this with fitting +
  top-X% filtering between levels. For plain `enumerate_mechanisms` (X=100%),
  every mechanism at each level passes through to expansion.

## Expansion Moves

Three separate functions handle expansion at different param_count deltas.

### `expand_mechanisms_by_one_param` — +1 Parameter Moves

1. **RE→SS**: Convert one RE step to SS. One candidate per RE step. Always
   exactly +1 param because the topology (forms, steps, cycles) is unchanged,
   so `n_thermo` is unchanged; RE contributes 1 param (K), SS contributes 2
   (kf, kr), net +1. Since each flip is independent, repeated application
   covers all 2^(n-1) RE/SS assignments reachable from the single-SS seed
   (the catalytic topology's isomerization step is always SS).

2. **Remove equivalence constraint**: Drop one parameter constraint (e.g., let
   K_dead_end ≠ K_catalytic). One candidate per active constraint. Always +1
   because each constraint eliminates exactly one parameter.

3. **Add dead-end binding**: For each metabolite or dead-end regulator not yet
   bound at a given form, generate all dead-end binding configurations that
   result in exactly +1 net parameter. Binding targets include both catalytic
   forms AND existing dead-end forms (multi-level dead-end), subject to the
   binding capacity limit (max(n_substrates, n_products) bound entities per
   form). The regulator binds to 1, 2, 3, ... eligible forms with maximal
   equivalence constraints:
   - All R-binding steps share one K_R
   - All substrate/product-binding steps on dead-end forms share K with their
     catalytic counterpart
   - Only configurations where the actual computed `param_count` of the
     resulting mechanism equals `source.param_count + 1` are kept
   Coverage argument: every dead-end configuration can be built by (a) starting
   with the maximally-constrained version (all K's shared), which is generated
   here, then (b) relaxing constraints one at a time via move 2 at subsequent
   levels. The order of relaxation doesn't matter — each relaxation is an
   independent +1 move.

4. **Remove TR equivalence** (allosteric only): Make one metabolite's K_T ≠ K_R.
   One candidate per metabolite currently in the TR-equivalent set.

### `expand_mechanisms_by_two_params` — +2 Parameter Moves

5. **Add allosteric regulation**: Convert a base mechanism to allosteric. Adds L
   (conformational equilibrium) + one K_T≠K_R. All remaining metabolites start
   TR-equivalent. Generates:
   - All variants of which single metabolite has K_T≠K_R
   - All `catalytic_n` values (1 to a configurable max, typically the number
     of subunits)
   - All regulator site partitions (groupings of allosteric regulators into
     binding sites with Bell-number enumeration)
   Only candidates where the actual computed `param_count` equals
   `source.param_count + 2` are kept. Note: `catalytic_n` and site partitions
   affect rate equation structure but not always param_count, so multiple
   structural variants may coexist at the +2 level.

### Regulator Role Handling

Regulators with known roles (`:dead_end`, `:allosteric`) are handled by their
respective moves (move 3 for dead-end, move 5 for allosteric). Regulators with
`:unknown` role generate candidates for BOTH roles independently: each unknown
regulator produces dead-end candidates (move 3) AND allosteric candidates
(move 5). This replaces the current pipeline's 2^n_unknown bitmask partitioning
with per-regulator, per-move exploration.

### `expand_mechanisms_same_param_count` — +0 Parameter Moves

Adds dead-end configurations that result in +0 net parameter change.

Concrete example: ordered bi-bi mechanism with forms E, EA, EAB, EPQ, E_Q.
Adding dead-end complex E_Q + A ⇌ E_QA where K_A is constrained to equal the
catalytic K_A (from E + A ⇌ EA). This adds 1 RE step (+1 param) and 1
equivalence constraint (-1 param), net = 0. The rate equation changes (E_QA
sequesters enzyme) but param_count does not.

Other +0 cases:
- Regulator R binding to 2 forms with both equivalence constraints (same K_R
  across forms AND substrate K matches catalytic): new cycle creates 1 thermo
  constraint that, combined with equivalence constraints, cancels the new
  parameters.
- Multi-level dead-end binding: adding a metabolite to an existing dead-end
  form where the K is shared with an existing step (e.g., S2 binding to
  E_S1_P1 with K_S2 = K_S2_catalytic).

Iterates to a fixed point (adding one +0 complex may enable another). Runs at
every level before +1/+2 expansion, because topology changes at each level
(e.g., new enzyme forms from adding a regulator) may enable new +0 moves.
Multi-level dead-end binding respects the binding capacity limit:
max(n_substrates, n_products) bound entities per form (counting substrates,
products, and dead-end regulators but not allosteric regulators).

## Interface

```julia
# New public enumerate_mechanisms — lazy, level-by-level
enumerate_mechanisms(
    reaction::EnzymeReaction;
    max_param_count=nothing,
) → MechanismIterator

# Three separate expansion functions (pure, stateless)
expand_mechanisms_same_param_count(
    specs::Vector{<:AbstractMechanismSpec},
    reaction::EnzymeReaction,
) → Vector{<:AbstractMechanismSpec}  # +0 dead-end additions

expand_mechanisms_by_one_param(
    specs::Vector{<:AbstractMechanismSpec},
    reaction::EnzymeReaction,
) → Vector{<:AbstractMechanismSpec}  # RE→SS, remove constraint,
                                     # add dead-end (+1), remove TR equiv

expand_mechanisms_by_two_params(
    specs::Vector{<:AbstractMechanismSpec},
    reaction::EnzymeReaction,
) → Vector{<:AbstractMechanismSpec}  # add allosteric regulation
```

`enumerate_mechanisms` returns `MechanismIterator` for backward compatibility.
Internally it calls all three expansion functions at each level.

`max_param_count` caps how far the expansion goes. Without it, it expands until
no new mechanisms are generated.

## Memory and Performance

Each `MechanismSpec` is ~2-5 KB. Per-level memory usage:

| Per-level count | Memory       | Fitting time (5s each) |
|----------------:|-------------:|----------------------:|
|           10,000 |   20-50 MB  |              14 hours |
|          100,000 |  200-500 MB |               6 days  |
|        1,000,000 | 2-5 GB      |              58 days  |

Only one level is held in memory at a time. Peak memory is ~2× one level
(current level + expansion candidates). For `identify_rate_equation` with
X < 100%, per-level count is bounded by beam_width × moves_per_mechanism.

A `max_mechanisms_per_level` safety cap (deferred) will warn/error if a level
exceeds a threshold.

## Dead-End Binding Rules

### Multi-Level Dead-End Binding

Dead-end forms can themselves be targets for further metabolite/regulator
binding. For example, in a ter-ter reaction: E_S1 (catalytic) → E_S1_P1
(dead-end) → E_S1_P1_S2 (dead-end binding to dead-end form).

This is a change from the current `old_enumerate_mechanisms` pipeline, which
only creates single-level dead-end complexes (binding to catalytic forms only).

### Binding Capacity Limit

Total bound entities per enzyme form ≤ max(n_substrates, n_products).

Entities that count toward the limit:
- Substrates
- Products
- Dead-end regulators (bind at catalytic site)

Entities that do NOT count:
- Allosteric regulators (bind at non-catalytic site)

Rationale: the catalytic site has limited physical space. For a ter-ter
reaction (limit = 3), forms like E_S1_S2_P1 (3 bound) are valid but
E_S1_S2_P1_P2 (4 bound) exceeds the limit.

### Maximal Constraints on Addition

When adding a dead-end regulator, the initial move is maximally constrained:

- All R-binding steps share one K_R across all enzyme forms
- All substrate/product-binding steps on dead-end forms share K with their
  catalytic counterpart

This minimizes the parameter count increase. Each constraint can be relaxed
individually at later levels via "remove equivalence constraint" moves.

## Thermodynamic Constraint Handling

Thermodynamic constraints (Haldane/Wegscheider) are handled by computing
`param_count` from the actual mechanism structure, not by assuming fixed deltas
per move type. Each expansion move generates complete `MechanismSpec` objects,
computes the actual `param_count`, and filters by target `param_count`.

The `param_count` computation follows the codebase convention: a base count
from `n_re + 2*n_ss - n_thermo + 2` (where `n_thermo = n_steps - n_forms + 1`
is the number of independent cycles, and `+2` accounts for E_total and Keq),
with equivalence constraints applied as a separate delta reduction. For SS
equivalence groups, each constrained pair eliminates 2 parameters (kf and kr);
for RE groups, each constrained pair eliminates 1 parameter (K).

## Testing Strategy

File: `test/test_beam_enumeration.jl`

1. **Catalytic topologies** — copied from current tests, same expected counts
2. **RE→SS expansion** — known mechanism, verify each candidate flips one step,
   param_count = original + 1
3. **Remove equivalence constraint** — mechanism with constraints, verify each
   drops one, param_count = original + 1
4. **Add dead-end regulator** — verify all configurations yielding exactly +1
   param (1 form, 2 forms with shared K, etc.), all with maximal equivalence
   constraints
5. **Add allosteric** — verify all variants of which single metabolite has
   K_T≠K_R, param_count = original + 2
6. **Remove TR equivalence** — allosteric mechanism with TR-equivalent
   metabolites, verify one metabolite becomes non-equivalent, param_count =
   original + 1
7. **expand_same_param_count** — verify +0 dead-end additions, iterate to fixed
   point
8. **Expansion integration** — all three expansion functions together produce
   the expected candidates, correctly bucketed by param_count delta
9. **Level-by-level equivalence** — for small reactions, verify new
   `enumerate_mechanisms` produces same final set as `old_enumerate_mechanisms`
10. **Deduplication within levels** — equivalent candidates from different
    expansion paths collapse to one
11. **param_count vs parameters()** — for a sample of compiled mechanisms at
    each level, verify `spec.param_count == length(parameters(compile_mechanism(spec)))`
12. **Multi-level dead-end binding** — for a bi-bi or ter-ter reaction, verify
    that dead-end forms can receive additional binding (e.g., E_S1_P1 + S2),
    respecting the binding capacity limit
13. **Binding capacity limit** — verify forms exceeding
    max(n_substrates, n_products) bound entities are not generated

## Scope (First PR)

**In scope:**
- Rename old pipeline and tests
- `expand_mechanisms_same_param_count` (+0 moves)
- `expand_mechanisms_by_one_param` (+1 moves, forward direction only)
- `expand_mechanisms_by_two_params` (+2 moves, forward direction only)
- Multi-level dead-end binding with binding capacity limit (new capability
  beyond `old_enumerate_mechanisms`)
- New `enumerate_mechanisms` using the beam loop with X=100%
- Tests for all of the above
- Equivalence test against `old_enumerate_mechanisms` for small reactions
  (excluding multi-level dead-end mechanisms, which are new)

**Deferred:**
- Reverse direction (max→min)
- `identify_rate_equation` integration
- `max_mechanisms_per_level` safety cap
- Removal of `old_enumerate_mechanisms` (after equivalence tests pass for all
  reactions up to bi-bi + 2 regulators)
