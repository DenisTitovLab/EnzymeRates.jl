# Mechanism Enumeration Test Redesign

## Problem

The current mechanism enumeration tests have several issues:

1. **Incorrect expected counts**: Some hardcoded values appear reverse-engineered from
   code output rather than independently calculated. Examples:
   - RE/SS expansion for Uni-Uni: expected 3 (2^2-1, only toggles binding edges) but
     should be 7 (2^3-1, all steps toggleable including isomerization)
   - Dead-end for Bi-Bi (no reg): expected 1 (passthrough) but should enumerate
     substrate/product dead-end complexes (e.g., EAQ, EBP)

2. **Missing pipeline stage**: Substrate/product dead-end complexes (enzyme forms where
   non-clashing substrates and products co-occupy the enzyme) are not generated. This
   requires a new stage between RE/SS expansion and regulator dead-end expansion.

3. **Insufficient isolation**: `StageExpansionTestSpec` tests each stage independently on
   a single base mechanism. Deduplication, equivalence constraints, and skip behavior
   are never tested with carefully constructed multi-mechanism inputs.

4. **No round-trip verification**: `compile_mechanism(spec)` is not systematically
   verified against the original `@enzyme_mechanism` definitions used to build test
   inputs.

## Design

### Scope

- **Test redesign only** — rewrite tests with correct hand-calculated expectations
- Existing pipeline bugs will produce failing tests marked `@test_broken`
- Add tests for future substrate/product dead-end stage (all `@test_broken`)
- No pipeline code changes

### File structure

```
test/
  mechanism_enumeration_test_specs.jl  — reaction definitions, helper functions
  test_mechanism_enumeration.jl        — rewritten, organized by stage:

    @testset "Stage 1: Catalytic topologies"
    @testset "Stage 2: RE/SS expansion"
    @testset "Stage 2.5: Substrate/product dead-end expansion"
    @testset "Stage 3: Regulator dead-end expansion"
    @testset "Stage 4: Equivalence constraints"
    @testset "Stage 5: Deduplication"
    @testset "Stage 6: Allosteric expansion"
    @testset "Stage 7: TR equivalence"
    @testset "Stage 8: Allosteric deduplication"
    @testset "End-to-end pipeline"
```

### Test methodology

Each stage testset:
- Defines test mechanisms inline using `@enzyme_mechanism` macro
- Converts to `MechanismSpec` via `mechanism_spec_from_mechanism`
- Verifies round-trip: `compile_mechanism(spec)` matches original mechanism
- Runs single stage function on hand-built inputs
- Checks counts against hand-calculated expectations
- Checks structural properties of outputs
- Tests passthrough/skip behavior where applicable
- Tests edge cases

### Test reactions

```julia
# No regulators
uni_uni             # S[C] → P[C]
uni_bi              # S[AB] → P[A] + Q[B]
bi_bi               # A[C] + B[N] → P[C] + Q[N]
bi_bi_ping_pong     # A[CX] + B[N] → P[C] + Q[NX]
ter_ter             # A[C] + B[N] + D[X] → P[C] + Q[N] + R[X]
# ter_ter_pp skipped — 4 ping-pong intermediates (E_X, E_C, E_N, E_Y),
# combinatorics too complex for initial test suite. TODO: add later.

# Dead-end inhibitor variants
uni_uni_dead_end_I
uni_bi_dead_end_I
bi_bi_dead_end_I
bi_bi_ping_pong_dead_end_I
uni_uni_dead_end_I_J          # 2 dead-end inhibitors

# Allosteric regulator variants
uni_uni_allosteric_R
uni_bi_allosteric_R
bi_bi_ping_pong_allosteric_R
uni_bi_allosteric_R_cn2       # catalytic_n=2

# Mixed
bi_bi_dead_end_I_allosteric_R
bi_bi_allosteric_R1_R2         # 2 allosteric regulators (for Stage 7 TR equiv)

# For end-to-end (regulator partitioning)
uni_uni_reg_unknown
uni_bi_reg_unknown
bi_bi_ping_pong_reg_unknown
```

## Stage-by-stage test specification

### Stage 1: Catalytic topologies

| Reaction | Expected topologies | Reasoning |
|----------|-------------------|-----------|
| uni_uni | 1 | Single cycle: E→ES→EP→E |
| uni_bi | 3 | 2 sequential + 1 random product release |
| bi_bi | 9 | 2 seq-bind/rand-release + 2 rand-bind/seq-release + 4 seq-bind/seq-release + 1 rand-bind/rand-release |
| bi_bi_ping_pong | 10 | 9 bi-bi + 1 ping-pong |
| ter_ter | 49 | 1 random + 6 rand-bind/seq-release + 6 seq-bind/rand-release + 36 seq/seq |

Verification approach:
- uni_uni and uni_bi: define ALL expected mechanisms by hand via `@enzyme_mechanism`
- bi_bi: define 3-4 representative mechanisms (one per category) + verify total count
- bi_bi_ping_pong: same as bi_bi + define the ping-pong mechanism
- ter_ter: verify count only (49 too many to define by hand)

Properties checked:
- Every topology forms a valid catalytic cycle
- `n_catalytic_edges == length(edges)`
- At least 1 SS step per topology
- Round-trip: `compile_mechanism(mechanism_spec_from_mechanism(m, rxn)) == m`

### Stage 2: RE/SS expansion

**Core formula**: For a mechanism with `n` total steps, `2^n - 1` variants (all combos
minus all-RE). Any step (binding or isomerization) can be SS.

| Input | n_steps | Expected | Formula |
|-------|---------|----------|---------|
| uni_uni (3 steps) | 3 | 7 | 2^3 - 1 |
| uni_bi topology (4 steps) | 4 | 15 | 2^4 - 1 |
| bi_bi topology (5 steps) | 5 | 31 | 2^5 - 1 |

**`@test_broken`**: Current code gives `2^(n_binding_edges) - 1` because it doesn't
toggle isomerization steps. Expected values use the correct `2^n - 1` formula.

**`max_re_groups` filtering**: Ter-Ter mechanism where the all-SS assignment creates
> 7 RE groups (each form is its own group). Verify `count < 2^n - 1`. Specifically
verify the all-SS assignment is excluded.

Properties:
- Every output has at least 1 SS step
- Every output has same edges as input (only `equilibrium_steps` changes)
- `n_catalytic_edges` unchanged

### Stage 2.5: Substrate/product dead-end expansion (future stage)

**All tests `@test_broken`** — stage does not exist yet.

Dead-end forms are all enzyme forms from `enumerate_enzyme_forms` that are NOT part of
the catalytic cycle for a given topology. This includes forms with atom-overlapping
metabolites (e.g., EAP in Bi-Bi where A and P share atom C) — we follow
`enumerate_enzyme_forms` as source of truth and let fitting filter biochemically
implausible mechanisms.

Expansion is per-topology: off-cycle forms depend on which forms the specific catalytic
topology uses.

| Input | Reaction | Off-cycle forms | Expected per input |
|-------|----------|----------------|-------------------|
| uni_uni topology | uni_uni | 0 (all forms on cycle) | 1 (passthrough) |
| bi_bi topology | bi_bi | 4 (EAP, EBP, EAQ, EBQ) | 2^4 = 16 |
| bi_bi_ping_pong topology | bi_bi_ping_pong | count from enumerate_enzyme_forms minus on-cycle | 2^n |
| ter_ter topology | ter_ter | count from enumerate_enzyme_forms minus on-cycle | 2^n |

Properties:
- Passthrough when no off-cycle forms exist
- New edges connect dead-end forms to catalytic forms only
- Dead-end edges inherit RE/SS from catalytic counterpart
- `n_catalytic_edges` unchanged

### Stage 3: Regulator dead-end expansion

**Formula**: `(2^n_regulators)^n_eligible_forms` per input spec — each eligible enzyme
form independently decides which subset of regulators bind to it.

**Inhibitor binding rule**: An enzyme form is **fully occupied** when it has ALL
substrates bound OR ALL products bound. Dead-end inhibitor I cannot bind to fully
occupied forms. All other forms (including substrate/product dead-end forms from
stage 2.5 that don't have all substrates or all products) are eligible.

Each substrate, product, and inhibitor can only bind to an enzyme once (unless the
inhibitor is listed multiple times, indicating multiple binding sites).

**Eligible forms by reaction type (catalytic cycle forms only, before stage 2.5):**

| Reaction | Fully occupied (ineligible) | Eligible | n_eligible |
|----------|---------------------------|----------|------------|
| Uni-Uni | ES (all subs), EP (all prods) | E | 1 |
| Uni-Bi | E_S (all subs), E_P_Q (all prods) | E, E_P, E_Q | 3 |
| Bi-Bi | E_A_B (all subs), E_P_Q (all prods) | E, E_A, E_B, E_P, E_Q | 5 |
| Ter-Ter | E_A_B_D, E_P_Q_R | E, all partial forms | many |

**After stage 2.5** (substrate/product dead-end forms):
- Bi-Bi E_A_Q: not all subs, not all prods → eligible
- Ter-Ter E_A_Q: eligible (partial occupancy)
- Ter-Ter E_A_B_Q: all subs bound → NOT eligible

| Input | Reaction | Regulators | n_eligible | Expected |
|-------|----------|-----------|------------|----------|
| uni_uni topology | uni_uni (no regs) | 0 | 1 | 1 (passthrough) |
| uni_uni topology | uni_uni_dead_end_I | 1 | 1 | (2^1)^1 = 2 |
| uni_bi topology | uni_bi_dead_end_I | 1 | 3 | (2^1)^3 = 8 |
| bi_bi topology | bi_bi_dead_end_I | 1 | 5 | (2^1)^5 = 32 |
| uni_uni topology | uni_uni_dead_end_I_J | 2 | 1 | (2^2)^1 = 4 |
| uni_uni topology | uni_uni_allosteric_R | allosteric only | — | 1 (passthrough) |

**`@test_broken`**: Current code allows inhibitors to bind to ALL catalytic forms
regardless of site occupancy. Tests with correct eligible form counts will fail.

Properties:
- Allosteric regulators are skipped (passthrough)
- New edges are all RE (regulator-binding)
- `n_catalytic_edges` unchanged
- Inhibitors only bind to forms that are not fully occupied

### Stage 4: Equivalence constraints

Test cases:
- **No equiv groups** (S and P bind different atoms): passthrough (count unchanged)
- **1 equiv group of size 2** (inhibitor I binds E and ES — same metabolite, same role):
  2 variants (unconstrained + K_I_E = K_I_ES)
- **Multiple equiv groups**: multiplicative expansion
- **Same metabolite as substrate AND dead-end regulator**: verify only substrate-binding
  edges are grouped together and only regulator-binding edges are grouped together.
  Substrate and regulator bindings of the same metabolite MUST NOT be put in the same
  equivalence group.
- **Substrate/product dead-end equivalence** (future, `@test_broken`): equivalence
  constraints must apply to substrate/product dead-end complexes as well

Properties:
- Output count ≥ input count
- Constrained variants have lower `param_count` than unconstrained
- Constraints reference valid edge indices

### Stage 5: Deduplication

**Hand-verified cases**:

1. **Uni-Uni**: 7 RE/SS variants → dedup to **1**. All variants have the same
   concentration fingerprint {1, [S], [P]} because the 3-form cycle is structurally
   symmetric. Verified by comparing `rate_equation_string` output — all produce
   equations with the same denominator monomial structure.

2. **Bi-Bi ordered** (single topology, 5 single-SS variants): all 5 have distinct
   fingerprints. Verified via `rate_equation_string`:

   | SS step | Denominator monomials |
   |---------|----------------------|
   | 1 (A bind) | {1, Q, PQ, PQ/B} |
   | 2 (B bind) | {1, A, Q, PQ} |
   | 3 (iso) | {1, A, AB, Q, PQ} |
   | 4 (P bind) | {1, A, AB, Q} |
   | 5 (Q bind) | {1, A, AB, AB/P} |

   No deduplication among single-SS variants.

**Additional test cases**:
- Feed 2 identical specs → verify dedup removes one
- Feed 2 specs with same fingerprint but different `param_count` → verify the lower
  `param_count` one is kept
- Multi-SS Bi-Bi dedup counts: compute programmatically as regression values

### Stage 6: Allosteric expansion

| Input | Allosteric regs | catalytic_n | Expected | Notes |
|-------|----------------|-------------|----------|-------|
| 1 spec, no regs | — | 1 | 1 | passthrough |
| 1 spec, 1 reg R | [R] | 1 | 1 | m=1 only |
| 1 spec, 1 reg R | [R] | 2 | 2 | m=1, m=2 |
| 1 spec, 1 reg R | [R] | 3 | 3 | m=1, m=2, m=3 |
| 1 spec, 2 regs R1, R2 | [R1, R2] | 1 | 2 | same site vs different sites |
| 1 spec, 2 regs R1, R2 | [R1, R2] | 2 | 6 | 2 partitions: separate sites (2×2=4) + same site (2) |
| 1 spec with dead-end I edges | — | 1 | 1 | dead-end passthrough |
| 1 spec with dead-end I + allosteric R | [R] | 1 | 1 | only R expands, I passes through |

Properties:
- Every output has `catalytic_n` matching input
- Every output has non-empty `allosteric_reg_sites` when regulators present
- `param_count` increases relative to base mechanism

### Stage 7: TR equivalence

**Formula**: `2^n` variants where n = number of metabolites with T-state params
(RE-binding metabolites + regulator ligands).

| Input | Metabolites with T-state params | Expected |
|-------|--------------------------------|----------|
| Uni-Uni + allosteric R | 3 (S, P, R) | 2^3 = 8 |
| Uni-Bi + allosteric R | 4 (S, P, Q, R) | 2^4 = 16 |
| Bi-Bi + allosteric R1, R2 | 6 (A, B, P, Q, R1, R2) | 2^6 = 64 |

Properties:
- Variants with non-empty `tr_equiv_metabolites` have fewer parameters
- The all-empty variant (no TR equiv) has most parameters
- The all-equivalent variant has fewest parameters

### Stage 8: Allosteric deduplication

T↔R mirror: two allosteric specs are mirrors when swapping T and R conformations
produces the same mechanism. For example, if metabolite X has K_T=a, K_R=b in one
spec and K_T=b, K_R=a in another, these are mirrors — the labels T and R are
arbitrary. The `tr_equiv_metabolites` field determines which metabolites have K_T=K_R
(making them invariant to T↔R swap).

Test cases:
- Hand-craft two specs that are T↔R mirrors (same base mechanism, complementary
  `tr_equiv_metabolites` sets) → verify they dedup to 1
- Two specs with different site assignments → verify both survive
- Two specs with different base mechanisms → verify both survive (even if
  `tr_equiv_metabolites` match)
- Verify dedup keeps the spec with fewer parameters when mirrors have different counts

### End-to-end pipeline

Keep existing `EnumerationTestSpec` structure with these reactions:
- uni_uni (no regs): simplest end-to-end validation
- uni_uni + 1 unknown reg: tests regulator partitioning (2^1 = 2 assignments)
- uni_bi (no regs): multiple topologies
- uni_bi + 1 unknown reg: combines multiple topologies with partitioning
- bi_bi (no regs): largest non-regulated case
- bi_bi_ping_pong (no regs): includes ping-pong topology

Expected counts will be updated after pipeline bugs are fixed. Currently-wrong counts
marked `@test_broken` with comments explaining the discrepancy. Counts that are
correct under current (buggy) behavior but wrong under correct behavior will be
replaced with correct values and marked `@test_broken`.

## Known bugs to be exposed

1. **RE/SS expansion**: Only toggles binding edges, not isomerization steps.
   Fix: toggle all steps. Expected: `2^n_total_steps - 1`.

2. **Missing stage 2.5**: Substrate/product dead-end complexes not enumerated.
   Fix: add new pipeline stage between stages 2 and 3.

3. **Regulator dead-end expansion**: Current code allows inhibitors to bind to ALL
   catalytic forms. Should only allow binding to forms that are not fully occupied
   (don't have all substrates or all products bound). Uncompetitive and noncompetitive
   inhibition patterns are special cases of allosteric control, not dead-end inhibition.

4. **Expected counts in existing tests**: Many appear reverse-engineered from code
   output rather than independently calculated. All will be replaced with
   hand-calculated values.
