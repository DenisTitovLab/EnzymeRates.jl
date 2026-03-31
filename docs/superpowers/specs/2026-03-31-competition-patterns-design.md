# Competition Patterns for Dead-End Form Enumeration

## Problem

`init_mechanisms` enumerates all 2^n subsets of dead-end substrate/product forms. For ter-ter (3+3 metabolites), the random topology has 27 dead-end forms → 2^27 = 134M variants per topology, causing OOM.

The root cause: the current code allows any substrate+product combination in dead-end forms, even when they compete for the same binding site. Biochemically, if A and P occupy overlapping active-site regions, they cannot coexist in any enzyme complex.

Similarly, `_expand_add_dead_end_regulator` enumerates all 2^n subsets of eligible forms for each inhibitor. Competition patterns should constrain which forms an inhibitor can bind to.

## Design

### Competition Pattern Definition

A **competition pattern** is a bipartite graph on substrates × products. An edge (S, P) means S and P compete for overlapping active-site regions and cannot coexist in any enzyme complex.

**Constraints:**
- Every substrate has degree ≥ 1 (competes with at least one product)
- Every product has degree ≥ 1 (competed by at least one substrate)

**Counts** (by inclusion-exclusion):
- Uni-uni (1×1): 1 pattern
- Bi-bi (2×2): 7 patterns
- Ter-ter (3×3): 265 patterns
- General formula: Σ over subsets S⊆subs, T⊆prods of (-1)^(|S|+|T|) × C(n_s,|S|) × C(n_p,|T|) × 2^((n_s-|S|)(n_p-|T|))

### Substrate/Product Dead-End Forms (init_mechanisms)

**Current behavior:** For each topology, enumerate all possible dead-end forms, then iterate over 2^n subsets.

**New behavior:** For each topology, enumerate competition patterns. For each pattern, the dead-end forms are **deterministic** — a form exists iff:
1. Its parent catalytic form exists in the topology
2. The binding satisfies existing rules (mixed binding, not all subs or all prods, etc.)
3. No two metabolites in the form are connected by a competition edge

No 2^n subset enumeration. Each (topology, competition pattern) pair yields exactly one mechanism.

**Dedup:** Different competition patterns may yield the same dead-end form set for a given topology. Dedup at the form-set level (cheap set comparison) before building MechanismSpec, avoiding duplicate mechanism construction.

### Dead-End Inhibitor Competition (expand_mechanisms)

When adding a dead-end inhibitor I via `_expand_add_dead_end_regulator`:

**Current behavior:** Enumerate all non-empty subsets of eligible enzyme forms (2^n).

**New behavior:** Enumerate inhibitor competition patterns. Each pattern determines which forms I can bind to (those not containing any metabolite that competes with I).

**Inhibitor competition pattern constraints:**
- I competes with ≥ 1 substrate (independent-site regulators are allosteric, not dead-end)
- I competes with ≥ 1 product
- At least one (S, P) pair exists where I competes with both S and P, AND S↔P is an edge in the substrate/product competition pattern

**Inhibitor pattern count:** For each substrate/product competition pattern, enumerate subsets of (substrate ∪ product) metabolites that I competes with, subject to the three constraints above. For ter-ter diagonal pattern: 37 inhibitor patterns (down from 49 without the third constraint).

**Determinism:** For a given inhibitor competition pattern, the set of forms I can bind to is deterministic — I binds to every eligible form that doesn't contain any metabolite competing with I. No 2^n form-subset enumeration.

**Dedup:** Same strategy as substrate/product dead-ends — dedup by the set of forms I binds to before building MechanismSpec.

### Where Competition Patterns Live

The substrate/product competition pattern is established in `init_mechanisms` and must be stored on `MechanismSpec` so that `_expand_add_dead_end_regulator` can reference it when generating inhibitor patterns.

**MechanismSpec change:** Add a `competition::Set{Tuple{Symbol,Symbol}}` field — the set of (substrate, product) competition edges. `AllostericMechanismSpec` inherits this from its `.base`.

All expansion moves (`_expand_re_to_ss`, `_expand_remove_constraint`, etc.) propagate the competition pattern unchanged.

## Source Changes

### 1. New function: `_competition_patterns(sub_names, prod_names)`

Enumerates all valid substrate/product competition patterns for a reaction.

```
Input: sub_names::Set{Symbol}, prod_names::Set{Symbol}
Output: Vector{Set{Tuple{Symbol,Symbol}}}
```

Iterates over all subsets of the n_s × n_p edge set, keeps those where every substrate and every product has degree ≥ 1.

### 2. New function: `_inhibitor_competition_patterns(sub_names, prod_names, competition)`

Enumerates all valid inhibitor competition patterns given a substrate/product competition pattern.

```
Input: sub_names, prod_names, competition::Set{Tuple{Symbol,Symbol}}
Output: Vector{Tuple{Set{Symbol}, Set{Symbol}}}  # (competing_subs, competing_prods)
```

Constraints:
- competing_subs non-empty
- competing_prods non-empty
- ∃ (s,p) with s ∈ competing_subs, p ∈ competing_prods, (s,p) ∈ competition

### 3. New function: `_dead_end_allowed(de_bound, competition, sub_names, prod_names)`

Returns true iff no (s, p) pair in the dead-end form's bound metabolites is a competition edge.

### 4. Modified: `MechanismSpec`

Add `competition::Set{Tuple{Symbol,Symbol}}` field. Provide an outer constructor with default empty set so existing call sites don't need changing:

```julia
MechanismSpec(reaction, steps, constraints, pc) =
    MechanismSpec(reaction, steps, constraints, pc,
        Set{Tuple{Symbol,Symbol}}())
```

~20 existing call sites remain unchanged. Only `_expand_substrate_product_dead_ends` and `init_mechanisms` pass the competition field explicitly.

### 5. Modified: `_expand_substrate_product_dead_ends`

Replace the `for mask in 0:(1 << n_de) - 1` loop with:

```
patterns = _competition_patterns(sub_names, prod_names)
seen = Set{Vector{Symbol}}()
for pattern in patterns
    filter dead-end forms by competition
    dedup by form set
    build mechanism with competition field set
end
```

### 6. Modified: `_expand_add_dead_end_regulator`

Replace the `for mask in 1:(1 << n_forms) - 1` loop with:

```
inh_patterns = _inhibitor_competition_patterns(
    sub_names, prod_names, spec.competition)
seen = Set{Vector{Symbol}}()
for (comp_subs, comp_prods) in inh_patterns
    filter eligible forms: exclude those containing
        any metabolite in comp_subs ∪ comp_prods
    dedup by form set
    build mechanism
end
```

### 7. Propagation in expansion moves

All functions that create new MechanismSpec from existing ones copy the `competition` field:
- `_expand_re_to_ss`
- `_expand_remove_constraint`
- `_expand_to_allosteric`
- `_expand_add_allosteric_regulator`
- `_expand_remove_tr_equiv`
- `_rewrap_allosteric`

## Test Plan

### Diagnostic tests for `_competition_patterns`

1. **Uni-uni:** 1 pattern (the single edge)
2. **Bi-bi:** 7 patterns; verify each has all 4 vertices covered
3. **Ter-ter:** 265 patterns; verify each has all 6 vertices covered
4. **Uni-bi:** verify count (should be 1 — the single substrate must compete with both products)
5. **Bi-uni:** verify count (symmetric)

### Diagnostic tests for substrate/product dead-end filtering

6. **Bi-bi random topology, diagonal competition {A↔P, B↔Q}:**
   - E_A_Q allowed (A,Q don't compete) ✓
   - E_B_P allowed (B,P don't compete) ✓
   - E_A_P forbidden (A↔P compete) ✗
   - E_B_Q forbidden (B↔Q compete) ✗
   - Exactly 2 dead-end forms

7. **Bi-bi random topology, complete competition {A↔P, A↔Q, B↔P, B↔Q}:**
   - All dead-end forms forbidden (every pair competes)
   - 0 dead-end forms

8. **Bi-bi {A↔P, B↔P} is not a valid pattern:**
   - Q has degree 0 (no substrate competes with Q)
   - Verify this subset is excluded from `_competition_patterns` output

9. **Ter-ter random topology, diagonal {A↔P, B↔Q, C↔R}:**
   - 1S+1P: E_A_Q, E_A_R, E_B_P, E_B_R, E_C_P, E_C_Q allowed (6); E_A_P, E_B_Q, E_C_R forbidden (3)
   - 2S+1P: E_AB_R, E_AC_Q, E_BC_P allowed (3); rest forbidden (6) — each forbidden form contains a competing pair
   - 1S+2P: E_A_QR, E_B_PR, E_C_PQ allowed (3); rest forbidden (6)
   - Total: **12 allowed** out of 27

### Diagnostic tests for `_inhibitor_competition_patterns`

10. **Bi-bi diagonal {A↔P, B↔Q}:**
    - Verify I competing with {A,P} is valid (A↔P is a competition edge) ✓
    - Verify I competing with {A,Q} is invalid (A↔Q not a competition edge, only pair) ✗
    - Verify I competing with {A,B,P,Q} is valid ✓
    - Count total valid patterns

11. **Ter-ter diagonal {A↔P, B↔Q, C↔R}:**
    - Count: 37 valid inhibitor patterns
    - Verify specific valid/invalid patterns

### Diagnostic tests for inhibitor dead-end form selection

12. **Uni-uni + inhibitor competing with {S, P}:**
    - I can only bind to E (form without S or P)
    - 1 dead-end form: E_I

13. **Bi-bi random + inhibitor competing with {A, P}:**
    - I can bind to forms without A and without P: E, E_B, E_Q, E_B_Q (if they exist as catalytic forms)
    - Verify exact forms

### Integration tests

14. **init_mechanisms bi-bi:** verify total mechanism count matches topologies × unique-dead-end-sets (should match or be less than current 16 per random topology)

15. **init_mechanisms ter-ter:** verify completes without OOM, count mechanisms

16. **Round-trip:** for a sample of generated mechanisms, verify `EnzymeMechanism(spec)` compiles and `length(parameters(m)) <= spec.param_count`

### Existing test changes

Current tests that check exact counts will change:
- **Bi-bi random dead-end count:** currently 16 (2^4 subsets of 4 dead-end forms) → will become ≤ 7 (one per competition pattern, after dedup)
- **Bi-bi Ping-Pong dead-end count:** currently 8 (2^3) → will change similarly
- **Dead-end regulator form counts:** currently 2^n subsets → will change to competition-pattern-determined counts
- **Integration test total counts:** will change throughout

These count changes are expected and correct — the new counts reflect biochemically valid mechanisms only.
