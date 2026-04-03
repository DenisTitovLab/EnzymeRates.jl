# Competition Patterns for Dead-End Form Enumeration

## Problem

`init_mechanisms` enumerates all 2^n subsets of dead-end substrate/product forms. For ter-ter (3+3 metabolites), the random topology has 27 dead-end forms → 2^27 = 134M variants per topology, causing OOM.

The root cause: the current code allows any substrate+product combination in dead-end forms, even when they compete for the same binding site. Biochemically, if A and P occupy overlapping active-site regions, they cannot coexist in any enzyme complex.

Similarly, `_expand_add_dead_end_regulator` enumerates all 2^n subsets of eligible forms for each inhibitor, which also needs constraining.

## Design

### Competition Pattern Definition

A **competition pattern** is a bipartite graph on substrates × products. An edge (S, P) means S and P compete for overlapping active-site regions and cannot coexist in any enzyme complex.

#### Substrate/product competition constraints

1. Every substrate competes with ≥ 1 product (degree ≥ 1)
2. Every product is competed by ≥ 1 substrate (degree ≥ 1)
3. A dead-end form must bind ≥ 1 substrate AND ≥ 1 product (mixed binding only — binding to free enzyme is not a dead-end)
4. A dead-end form cannot bind all substrates or all products (those are catalytic forms)
5. No two metabolites in a dead-end form can be connected by a competition edge (competing metabolites cannot coexist)
6. For a given (topology, competition pattern), the set of dead-end forms is **deterministic**: a form exists iff its parent catalytic form exists, it satisfies constraints 3-5, and it's not already a catalytic form

#### Dead-end inhibitor competition constraints

A dead-end inhibitor I binds to the catalytic site (not allosteric). Its competition pattern against substrates/products determines which enzyme forms it can bind to:

1. I competes with ≥ 1 substrate (independent-site regulators are allosteric, not dead-end)
2. I competes with ≥ 1 product (catalytic site overlaps with both substrate and product regions)
3. I can only bind to catalytic form F if F is the **source of a binding step** for at least one metabolite that I competes with. I blocks a specific metabolite's binding event, so it can only bind where that event exists in the topology.
4. For a given inhibitor competition pattern and topology, the set of forms I binds to is **deterministic**: take the union of source forms across all competing metabolites' binding steps.
5. An eligible form has neither all substrates nor all products bound (same rule as substrate/product dead-ends)

**Example — sequential bi-bi (A→B→catalysis→P→Q), I competes with {B, P}:**
- B-binding step: E_A + B → E_AB. Source form: E_A
- P-release step: E_PQ → E_Q + P. Source form: E_Q
- I binds to: {E_A, E_Q} — NOT free enzyme E (no B or P binding step from E)

**Example — random bi-bi, I competes with {B, P}:**
- B-binding steps: E + B → E_B, E_A + B → E_AB. Sources: {E, E_A}
- P-binding steps: E + P → E_P, E_Q + P → E_PQ. Sources: {E, E_Q}
- I binds to: {E, E_A, E_Q} — includes E because B and P binding steps exist from E

**Pattern counts** (by inclusion-exclusion formula):

| Reaction type | S/P competition patterns | Inhibitor patterns per S/P pattern |
|---------------|------------------------|------------------------------------|
| Uni-uni (1×1) | 1 | (2^1 - 1)(2^1 - 1) = 1 |
| Uni-bi (1×2) | 1 | (2^1 - 1)(2^2 - 1) = 3 |
| Bi-bi (2×2) | 7 | (2^2 - 1)(2^2 - 1) = 9 |
| Ter-ter (3×3) | 265 | (2^3 - 1)(2^3 - 1) = 49 |

S/P formula: Σ_{S⊆subs, T⊆prods} (-1)^(|S|+|T|) C(n_s,|S|) C(n_p,|T|) 2^((n_s-|S|)(n_p-|T|))

Inhibitor formula: (2^n_s - 1)(2^n_p - 1) — non-empty substrate subset × non-empty product subset.

### Substrate/Product Dead-End Forms (init_mechanisms)

**Current behavior:** For each topology, enumerate all possible dead-end forms, then iterate over 2^n subsets.

**New behavior:** For each topology, enumerate S/P competition patterns. For each pattern, dead-end forms are deterministic per constraint 6 above. No 2^n subset enumeration.

**Dedup:** Different competition patterns may yield the same dead-end form set for a given topology. Dedup at the form-set level (cheap set comparison) before building MechanismSpec, avoiding duplicate mechanism construction.

### Dead-End Inhibitor Competition (expand_mechanisms)

**Current behavior:** Enumerate all non-empty subsets of eligible enzyme forms (2^n).

**New behavior:** Enumerate inhibitor competition patterns. For each pattern, find all catalytic forms that are sources of binding steps for competing metabolites (using `step_metabolite(s)` and `step_forms(s)`). This set is deterministic per topology.

**Dedup:** Same strategy — dedup by the set of forms I binds to before building MechanismSpec.

### No MechanismSpec Changes

Competition patterns are used at two enumeration stages and need not be stored:
- `_expand_substrate_product_dead_ends`: enumerates S/P patterns from the reaction
- `_expand_add_dead_end_regulator`: enumerates inhibitor patterns from the reaction

Both derive competition patterns from the `EnzymeReaction`, not from the mechanism. The resulting MechanismSpec already encodes competition through which steps/forms exist.

### Code Reuse

S/P and inhibitor dead-end expansion differ in how they select forms (competition-edge filtering vs step-source lookup) but share:

1. **Competition pattern enumeration**: both enumerate subsets with coverage constraints. `_competition_patterns` (bipartite edges with min degree) and `_inhibitor_competition_patterns` (non-empty sub/prod subsets) are similar loops over bitmasks. Consider a shared `_nonempty_covered_subsets` helper.

2. **Dedup-then-build**: both dedup by form set before constructing MechanismSpec. The build loop (binding steps + mirror steps) is structurally identical — extract shared step-building logic.

## Source Changes

### 1. New function: `_competition_patterns(sub_names, prod_names)`

Enumerates all valid substrate/product competition patterns.

```
Input: sub_names::Set{Symbol}, prod_names::Set{Symbol}
Output: Vector{Set{Tuple{Symbol,Symbol}}}
```

Iterates over all subsets of the n_s × n_p edge set, keeps those where every substrate and every product has degree ≥ 1.

### 2. New function: `_inhibitor_competition_patterns(sub_names, prod_names)`

Enumerates all valid inhibitor competition patterns.

```
Input: sub_names::Set{Symbol}, prod_names::Set{Symbol}
Output: Vector{Tuple{Set{Symbol}, Set{Symbol}}}  # (competing_subs, competing_prods)
```

Returns all (non-empty substrate subset, non-empty product subset) pairs. Count: (2^n_s - 1)(2^n_p - 1).

### 3. New function: `_forms_with_binding_step(steps, metabolite)`

Returns the set of source forms that have a binding step for `metabolite`. Uses `step_metabolite(s)` and `step_forms(s)` to find steps involving the metabolite, returns the form that doesn't contain it (the source side).

Used by `_expand_add_dead_end_regulator` to determine where an inhibitor can bind.

### 4. Modified: `_expand_substrate_product_dead_ends`

Replace the `for mask in 0:(1 << n_de) - 1` loop with:

```
patterns = _competition_patterns(sub_names, prod_names)
seen = Set{Vector{Symbol}}()
for pattern in patterns
    competing_pairs = pattern
    allowed_de = [de for de in de_form_names
        if no (s,p) pair in bound[de] is in competing_pairs]
    allowed_de in seen && continue
    push!(seen, allowed_de)
    build mechanism from allowed_de
end
```

### 5. Modified: `_expand_add_dead_end_regulator`

Replace the `for mask in 1:(1 << n_forms) - 1` loop with:

```
inh_patterns = _inhibitor_competition_patterns(sub_names, prod_names)
seen = Set{Vector{Symbol}}()
for (comp_subs, comp_prods) in inh_patterns
    # Find forms where competing metabolites have binding steps
    target_forms = Set{Symbol}()
    for met in union(comp_subs, comp_prods)
        union!(target_forms, _forms_with_binding_step(spec.steps, met))
    end
    # Filter to eligible forms only
    allowed = sort([f for f in target_forms if f in eligible_forms])
    isempty(allowed) && continue
    allowed in seen && continue
    push!(seen, allowed)
    build mechanism with inhibitor binding to allowed
end
```

## Test Plan

Tests follow TDD: write failing test with expected count first, then implement.

### Diagnostic tests for `_competition_patterns`

1. **Uni-uni:** 1 pattern (the single edge S↔P)
2. **Bi-bi:** 7 patterns; verify each has all 4 vertices with degree ≥ 1
3. **Ter-ter:** 265 patterns; verify each has all 6 vertices with degree ≥ 1
4. **Uni-bi (1S, 2P):** 1 pattern (S must compete with both products)
5. **Bi-uni (2S, 1P):** 1 pattern (symmetric)

### Diagnostic tests for substrate/product dead-end filtering

6. **Bi-bi random topology, diagonal competition {A↔P, B↔Q}:**
   - E_A_Q allowed (A,Q don't compete) ✓
   - E_B_P allowed (B,P don't compete) ✓
   - E_A_P forbidden (A↔P compete) ✗
   - E_B_Q forbidden (B↔Q compete) ✗
   - Exactly 2 dead-end forms

7. **Bi-bi random topology, complete competition {A↔P, A↔Q, B↔P, B↔Q}:**
   - All dead-end forms forbidden (every S/P pair competes)
   - 0 dead-end forms

8. **Bi-bi {A↔P, B↔P} is not a valid pattern:**
   - Q has degree 0 (no substrate competes with Q)
   - Verify excluded from `_competition_patterns` output

9. **Ter-ter random topology, diagonal {A↔P, B↔Q, C↔R}:**
   - 1S+1P: 6 allowed (E_A_Q, E_A_R, E_B_P, E_B_R, E_C_P, E_C_Q), 3 forbidden
   - 2S+1P: 3 allowed (E_AB_R, E_AC_Q, E_BC_P), 6 forbidden
   - 1S+2P: 3 allowed (E_A_QR, E_B_PR, E_C_PQ), 6 forbidden
   - Total: **12 allowed** out of 27

### Diagnostic tests for `_inhibitor_competition_patterns`

10. **Uni-uni:** 1 pattern — I competes with {S} and {P}
11. **Bi-bi:** 9 patterns — (2^2-1)(2^2-1) = 3×3
12. **Ter-ter:** 49 patterns — (2^3-1)(2^3-1) = 7×7

### Diagnostic tests for inhibitor dead-end form selection

13. **Uni-uni + inhibitor competing with {S, P}:**
    - S-binding step from E, P-binding step from E
    - I binds to E → 1 dead-end form: E_I

14. **Bi-bi sequential (A→B, P→Q release) + I competes with {B, P}:**
    - B-binding step: E_A + B → E_AB. Source: E_A
    - P-release step: E_PQ → E_Q + P. Source: E_Q
    - I binds to {E_A, E_Q} — NOT E (no B or P step from E)
    - 2 dead-end forms: E_A_I, E_Q_I

15. **Bi-bi random + inhibitor competing with {A, P}:**
    - A-binding steps from: E, E_B, E_Q, ... (all forms without A)
    - P-binding steps from: E, E_A, E_B, ... (all forms without P)  
    - Union of source forms, filtered to eligible
    - Verify exact forms

16. **Bi-bi random + inhibitor competing with {A, B, P, Q}:**
    - A/B/P/Q binding steps all include E as source
    - I binds to E (and possibly others) → verify exact set

### Integration tests with expected counts

17. **init_mechanisms bi-bi:** compute expected total = Σ over topologies of (unique dead-end sets per topology across 7 competition patterns). Verify exact count.

18. **init_mechanisms ter-ter:** verify completes without OOM. Verify count is reasonable (265 patterns × n_topologies, minus dedup).

19. **expand_mechanisms with dead-end regulator:** for a specific bi-bi mechanism + regulator, verify inhibitor variant count matches 9 inhibitor patterns minus dedup.

20. **Round-trip:** for a sample of generated mechanisms, verify `EnzymeMechanism(spec)` compiles and `length(parameters(m)) <= spec.param_count`

### Existing test changes

Current tests that check exact counts will change:
- **Bi-bi random dead-end count:** currently 16 (2^4) → ≤ 7 (competition patterns after dedup)
- **Bi-bi Ping-Pong dead-end count:** currently 8 (2^3) → will change similarly
- **Dead-end regulator form counts:** currently 2^n subsets → competition-pattern-determined counts
- **Integration test total counts:** will change throughout

These count changes are expected — the new counts reflect biochemically valid mechanisms only.
