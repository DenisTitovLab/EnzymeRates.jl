# Investigation: cross-hash duplication in EnzymeRates `eq_hash`

## Problem

`eq_hash` is missing algebraically-equivalent rate equations: distinct
`EnzymeMechanism` structs with different `eq_hash` values converge to the same
loss because their printed rate equations are mathematically identical despite
being textually different. Within-hash multiplicity (multiple structs sharing
one `eq_hash`) is by design — preserved for divergence at higher `n_params` —
and is NOT in scope here.

Empirically (LDH, `params_estimate_{5,6,7,8}.csv`):

| n | distinct eq_hash | distinct true rate eq | cross-hash dup % |
|---|------------------|-----------------------|-------------------|
| 5 |  13 |  11 | 15% |
| 6 |  62 |  24 | 61% |
| 7 | 493 | 156 | 68% |
| 8 |1390 | 467 | 66% |

"distinct true rate eq" = number of distinct loss values to ~10 sig figs
(astronomically unlikely to collide accidentally).

The four n=7 mechanisms at `loss = 0.012722696603…` (rows 22, 27, 31, 36 of
`cv_results.csv`, eq_hashes `831e36af`, `9c7141ac`, `89f33d51`, `b362dd75`) are
the canonical example. All four produce the same algebraic rate law and fit to
the same 7 parameters with values agreeing to 5–6 sig figs.

## Root cause: three independent textual sources

`identify_rate_equation.jl:151` `_canonicalize_rate_eq_with_map` already
alpha-renames parameter symbols to `p_1, p_2, …` in first-appearance order
before hashing, so the hash is insensitive to parameter *names*. The remaining
sensitivity is to the *structural form* of the printed expression. Three
independent sources of structural variation produce cross-hash duplicates:

**Source A — factored vs. factored (different factorings of the same polynomial).**
`_rate_v_line` (`rate_eq_derivation.jl:750`) emits the numerator via
`_factored_sigma_to_expr` and the denominator via `_denom_terms_to_expr`. Two
mechanisms can naturally fall out of the generator with different factorings
of the same polynomial. Among the four n=7 mechanisms: Mech 2/4's leading
denominator factor is `(1 + NAD/K2) * (1 + Pyruvate/K4)`; Mech 1/3's is
`(1 + Pyruvate/K4 + (NAD/K2) * (1 + Pyruvate/K8))`. These expand to the same
sum of monomials.

**Source B — `1/(1/X)` alias artifact (Ka↔Kd inversion).**
`_apply_kd_inversion` (`thermodynamic_constr_for_rate_eq_derivation.jl:295`)
double-wraps with `inv_fn = K → 1/K` when both the dependent symbol and the
RHS are binding-K's, producing constraint lines like `K8 = 1/(1/K4)` and
`K12 = 1/(1/K2)`. These dominate the K-LHS constraint forms in the data: 1507
of 2457 K-LHS constraints (61%) are the trivial `K = 1/(1/X)` form.

**Source C — parameter-tying redundancy (split-with-tie ≡ pre-merged).**
For each mechanism with a "tied K" relationship `K8 = K4` (modeler-chosen
shared affinity across two binding events), the enumerator also produces the
parametrically-equivalent "pre-merged" mechanism that has only `K4` and no
`K8`. These have identical rate equations but differ structurally in whether
a `K8` symbol exists at all in the body. Source C is what survives after
fixing A and B — see "residual" example in the table below.

Effect on the four n=7 mechanisms:

| Pair | Differ in | Source |
|------|-----------|--------|
| Mech 2 ↔ Mech 1 | Different factorings of same polynomial | A |
| Mech 1 ↔ Mech 3 | Both define `K8`; only Mech 1 uses it in `v` | B |
| Mech 2 ↔ Mech 4 | Mech 4 has dead alias `K12 = 1/(1/K2)` | B |
| Mech 1 ↔ Mech 2 (post-A,B) | Mech 1 still has `K8 = K4` line; Mech 2 has no `K8` | C |

All three need addressing for a full collapse.

## Proposed fix: three-step refactor that simplifies the codebase

Rather than patching the canonicalizer to mask each artifact, eliminate them
at the source. Each step is independently shippable; together they should
collapse `eq_hash` to match the distinct-loss count, and the codebase shrinks
meaningfully (~94 occurrences of inversion + factoring machinery across
`rate_eq_derivation.jl`, `sym_poly_for_rate_eq_derivation.jl`, and
`thermodynamic_constr_for_rate_eq_derivation.jl`).

### Step 1 — K_d-by-construction at the polynomial level

Eliminate the `Ka` (association) representation throughout. Build the
symbolic polynomials directly in terms of dissociation constants `K_d`. Removes
Source B entirely: with no inversion, no `1/(1/X)` text can be produced.

What goes away:
- `_apply_kd_inversion`
  (`thermodynamic_constr_for_rate_eq_derivation.jl:295`) and all its callers.
- `_binding_K_symbols`, `inverted_params`, the `inv_set` parameter threaded
  through `to_rate_expr`, `_factored_sigma_to_expr`, `_denom_terms_to_expr`,
  `_poly_to_expr`.
- The `inv_fn = K → 1/K` callbacks in
  `thermodynamic_constr_for_rate_eq_derivation.jl:312, 348`.
- The `K → 1/K` mental layer between internal storage and printed form.

What changes: polynomial construction uses `K_d` directly, so monomials
contain `1/K_d` rather than `K_a`. Rate-equation evaluation is mathematically
identical; the printed form matches the internal form.

`eq_hash` improves incidentally; the primary win is code clarity.

### Step 2 — always-expanded emission for non-allosteric

Drop algebraic factoring of the printed numerator and denominator. Use the
existing expansion + `_poly_to_expr` infrastructure. Allosteric enzymes
retain factoring because it happens at derivation and is structurally
required (R-state and T-state factor as `(1 + L*…)`); the change is scoped
to non-allosteric.

What goes away (or becomes allosteric-only):
- `_factor_poly` and `_try_algebraic_factor_sigma`
  (`rate_eq_derivation.jl:405`, `rate_eq_derivation.jl:210`).
- `_factored_sigma_to_expr` and `_denom_terms_to_expr`
  (`sym_poly_for_rate_eq_derivation.jl:392, 416`).
- `FactoredPoly`, `FactoredSigma`, `DenomTerm` types — or scoped to
  the allosteric path.
- The `check_benefit` heuristic and `_estimate_expanded_term_count`.

`_rate_v_line` (`rate_eq_derivation.jl:750`) becomes ~4 lines: expand num,
expand den, build canonical exprs via `_poly_to_expr`, format. Removes
Source A.

`MAX_RATE_EQUATION_TERMS = 5000` (`sym_poly_for_rate_eq_derivation.jl:7`)
already caps the *expanded* term count — current code estimates expansion
size before deciding to factor. Switching to always-expanded does not
introduce a new failure mode; mechanisms that compile today continue to.
Sanity-check that LDH/PGD/etc. at the largest enumerated `n_params` stay
under the cap.

### Step 3 — mechanism-level canonicalization that pre-merges tied K's

When the enumerator would produce a mechanism with a `K_a = K_b` (single-
symbol equality) constraint between two binding-K's, pre-merge those
elementary-step groups instead of emitting a tie. After this, no trivial
`K_a = K_b` line ever appears — only genuine cycle-balance constraints
(multi-symbol RHS) and Haldane closures for `k_r`. Removes Source C.

What changes: the enumerator's group-partition logic. When two groups would
be tied to the same K, return the merged partition rather than the split-
with-tie partition. The merged form has fewer groups and the same number of
fitted parameters.

This is the same change as eliminating "category 3" from
`investigation_eq_hash_duplication`'s constraint-classification table:
parameter-tying-by-modeler-assumption is no longer enumerated as a separate
mechanism from its pre-merged equivalent. Categories 1 (cycle-balance) and
2 (allosteric T/R mirroring) remain — they're inherent to chemistry/
allostery and don't suffer the same redundancy.

After this step, the LDH `cv_results.csv` and `params_estimate_*.csv` files
would also be much shorter — fewer "different mechanisms" to track at every
`n_params` level.

## Validation

After each step, run on `params_estimate_7.csv` (and 5/6/8 for completeness):

1. Recompute `eq_hash` for all rows using the post-step canonicalizer.
2. Group rows by loss to 10 sig figs and check that within each group all
   rows share the new `eq_hash`.
3. Expected progression toward `distinct eq_hash == distinct true rate eq`:
   - After Step 1 (K_d-by-construction): partial collapse — Source-B-only
     pairs (e.g., Mech 1 ↔ Mech 3, Mech 2 ↔ Mech 4) merge. Source A and C
     pairs remain split.
   - After Step 2 (always-expanded): further collapse — Source-A pairs
     (e.g., Mech 1 ↔ Mech 2) merge.
   - After Step 3 (mechanism canonicalization): full collapse —
     `eq_hash` count should match distinct-loss count
     (493 → 156 at n=7, 1390 → 467 at n=8, etc.).

If any residual cross-hash duplicates remain after Step 3, dump the canonical
strings for two non-collapsing mechanisms and diff them to identify the next
quirk.

## Fallback: localized canonicalizer patches

If any of the three steps proves too invasive to land at once, a textual
patch in `_canonicalize_rate_eq_with_map` (`identify_rate_equation.jl:151`)
can simulate the collapse for hashing only, leaving the runtime rate
equation untouched:

- For Source B: `body = replace(body, r"\(\s*1\s*/\s*\(\s*1\s*/\s*(\w+)\s*\)\s*\)" => s"\1")`.
- For Source C: detect single-symbol equality lines `X = Y` (both bare
  symbols) and substitute `X → Y` in the body, then drop the line.
- For Source A: requires expanded-form polynomial canonicalization, which
  is essentially Step 2 narrowed to the hash path.

These fallbacks deliver the `eq_hash` correctness without the code
simplification benefit. Use them only as a stop-gap.

## Out of scope

- Within-hash multiplicity (multiple `mechanism_type` structs sharing one
  `eq_hash`). By design.
- The "extra parameter pinned at 1.0" pattern. Consequence of within-hash
  multiplicity at the current `n_params` level; resolves at higher `n_params`
  when the parameter actually enters the rate equation.
- Changing the optimizer pipeline, beam search, or model-selection logic.
