# Catalytic Topology Enumeration Constraints — Design Spec

## Problem

The current `_catalytic_topologies` backtracking produces 3,969 catalytic
topologies for a simple ter-ter reaction (A[C]+B[N]+D[X] → P[C]+Q[N]+R[X]),
which expands to ~1M mechanisms after dead-end expansion — causing OOM.

Two independent problems cause this:
1. **Incorrect path combining**: the code produces 3,969 = 63 × 63
   topologies by taking all 2^6-1=63 non-empty subsets of 6 substrate
   binding paths × 63 product release subsets. Most subsets don't correspond
   to valid mechanism types. For example, subset {ABC, BCA} allows both A
   and B to bind to free E (implying random access at that level), but then
   E_A allows only B to bind next while E_B allows only C — an asymmetric
   behavior at structurally equivalent forms. The correct count using weak
   orderings (Fubini numbers) is 13 substrate variants × 13 product
   variants = 169.
2. **Missing ping-pong mechanisms**: the current backtracking requires
   non-empty residual for ping-pong and only supports releasing 1 product
   or all products at isomerization. This misses biochemically important
   mechanisms like pyruvate carboxylase (ATP+HCO₃ → ADP+Pi+CO₂_residual,
   then Pyr+CO₂ → OAA).

## Goal

Fix the path combining logic to use weak orderings instead of arbitrary
subsets, and extend the backtracking to generate ping-pong mechanisms with
empty residual and multi-product release at isomerization.

## Definitions

- **Step**: A single elementary event in the mechanism — either metabolite
  binding/release (`[E, S] ⇌ [ES]`) or isomerization (`[ES] ⇌ [EP]`).
- **Binding/release step**: A step where one metabolite binds to or
  dissociates from an enzyme form. Exactly one metabolite on one side.
- **Isomerization step**: A step where one enzyme form converts to another
  with no metabolite binding/release. Single form on each side.
- **Enzyme form**: A node in the mechanism graph representing a distinct
  state of the enzyme (e.g., `E`, `E_A`, `Estar`, `Estar_P`).
- **Estar**: Enzyme in a conformationally changed state after a ping-pong
  isomerization. Distinct from `E` even if no atoms are covalently attached.
- **Residual**: Atoms remaining on the enzyme after an isomerization releases
  some products. Counted as a product entity for constraint purposes.
- **Catalytic path**: A single complete catalytic cycle
  (`E → ... → E`) — one specific ordering of binding, isomerization,
  and release events.
- **Catalytic topology**: A mechanism graph formed by combining one or more
  catalytic paths that share the same isomerization steps. Contains all
  binding/release orderings from the constituent paths.

## Constraints

### C1. Path combining by weak orderings (bug fix — implement first)

Catalytic paths that share identical isomerization steps (same set of
`form_A → form_B` conversions) but differ in binding/release ordering are
combined into topologies using **weak orderings** (also called Fubini
numbers or ordered Bell numbers).

A weak ordering of n metabolites partitions them into ordered priority
levels. All metabolites in the same level can bind/release in random order;
levels are strictly ordered. For n metabolites:

- n=1: 1 ordering (trivial)
- n=2: 3 orderings (AB, BA, A/B random)
- n=3: 13 orderings (6 fully ordered + 6 partially random + 1 fully random)

For 3 substrates, the 13 orderings are:
1. 6 fully ordered: ABC, ACB, BAC, BCA, CAB, CBA
2. 6 partially random: A{B/C}, {B/C}A, B{A/C}, {A/C}B, C{A/B}, {A/B}C
3. 1 fully random: {A/B/C}

**Rationale**: If two substrates A and B both bind to the same enzyme form
(e.g., free enzyme E), they occupy the same priority level — their binding
sites are both accessible in that conformation. After either one binds, the
remaining substrates' accessibility must be consistent: E_A and E_B must
offer the same set of next binders. The current code's arbitrary subset
combining violates this. For example, subset {ABC, BCA} allows A and B to
both bind to E, but then E_A allows only B next while E_B allows only C
next — an asymmetry that implies structurally equivalent forms (E_A, E_B)
have different binding site accessibility, which is biochemically
implausible.

Each half-reaction has its own set of metabolites to order (substrates
binding, products releasing). The weak orderings for each half-reaction are
independent, so the total topologies per isomerization pattern =
product of weak ordering counts across half-reactions.

Paths are combined regardless of whether individual paths use `Estar` forms
(ping-pong) or not — the ping-pong/sequential partition from the current
code is superseded by this isomerization-based grouping.

*Bug fix: current code takes all 2^n-1 subsets of n paths, producing
biochemically meaningless ordering variants (e.g., 63 × 63 = 3,969 instead
of the correct 13 × 13 = 169 for 3-substrate × 3-product sequential
mechanisms).*

### C2. One metabolite per step

Each binding/release step involves exactly one metabolite on one side
of the reaction. No step has 2+ metabolites binding or releasing
simultaneously.

*Unchanged from current code. Enforced in `types.jl`.*

### C3. Atom conservation

Isomerization is only valid when the accumulated substrate atoms on the
enzyme contain (as a subset) the atoms of the products being released.
Validated via `_can_pingpong` checks.

*Unchanged from current code.*

### C4. Empty-residual ping-pong allowed

Ping-pong isomerization is allowed even when the residual is empty (no atoms
remain on the enzyme after conversion). The enzyme can be in a
conformationally changed state (`Estar`) without covalently attached atoms.

**Rationale**: Enzymes using ATP hydrolysis often undergo conformational
changes that don't leave covalently attached atoms but do change the enzyme's
binding properties.

*Changed from current code, which requires non-empty residual
(`isempty(residual) && continue`).*

### C5. Maximum bound metabolites

At most `max(n_substrates, n_products)` metabolites may be simultaneously
bound to the enzyme at any point in the mechanism.

**Rationale**: Catalytic site(s) have limited capacity. For a ter-ter
reaction, at most 3 metabolites can be bound simultaneously.

*New constraint. Current code has no binding limit.*

### C6. Isomerization size limit

At any isomerization step:

```
n_substrates_reacting ≤ 3  AND  n_products_effective ≤ 3
```

where `n_products_effective = n_released_products + (1 if residual is
non-empty)`. Residual counts as a product because it represents converted
substrate material that hasn't been released yet.

This is a hard cap at 3 per side. For bi-bi reactions this is a no-op
(already ≤ 2). For ter-ter, 3→3 sequential isomerization
(`E_A_B_D → E_P_Q_R`) is allowed — with the weak ordering fix (C1), this
produces a manageable 13 × 13 = 169 topologies. For hypothetical quad-quad
reactions, this forces at least one ping-pong step (4→4 blocked).

**Allowed examples**:
- 3 subs → 3 prods (no residual) = 3→3 ✓ (sequential ter-ter)
- 2 subs → 1 prod + residual = 2→2 ✓ (pyruvate carboxylase: ATP+HCO₃ →
  ADP+Pi+CO₂_residual, where CO₂_residual is the residual)
- 2 subs → 2 prods (no residual) = 2→2 ✓
- 1 sub → 1 prod + residual = 1→2 ✓
- 1 sub → 1 prod (no residual) = 1→1 ✓

**Forbidden**: either side > 3 (e.g., 4 substrates reacting simultaneously).

Note: 3→2+residual (=3→3) and 3→1+residual (=3→2) are allowed by this
constraint, but paths using them dead-end if no remaining substrates exist
to bind (C7 prevents Estar from isomerizing without substrate input). In
practice, the only valid 3-sub iso for ter-ter is the full 3→3 sequential
conversion (no residual).

*New constraint. Current code allows unlimited isomerization size.*

### C7. Isomerization requires substrate participation

Every isomerization step must have at least one substrate bound to the
enzyme. An `Estar` form with only residual (no substrates bound) cannot
isomerize — it must first bind a substrate.

If `Estar` has residual atoms and no remaining substrates to bind, the
catalytic path is a dead end (invalid path, pruned during backtracking).

**Rationale**: Residual atoms on the enzyme cannot spontaneously rearrange
into products without chemical input from a substrate. The isomerization
is a chemical reaction, not a spontaneous decomposition.

*New constraint. Current code allows `Estar → E_products` without substrate
input.*

### C8. Isomerization forms contain only products

After an isomerization step, the resulting enzyme form name includes only
the products ready for release — not the consumed substrates. For example,
if `E_A_B` isomerizes releasing product P, the resulting form is `Estar_P`,
not `Estar_A_B_P`.

**Rationale**: Substrates are consumed in the isomerization. The enzyme
form after isomerization holds products (ready for release) and possibly
residual atoms (which will participate in a future isomerization with
another substrate). Mixed substrate-product forms are biochemically
implausible in catalytic mechanisms.

*Changed from current code, which includes both substrate and product names
in isomerization form names (e.g., `Estar_A_P`).*

### C9. Multi-product release at isomerization

A single isomerization step may produce multiple products on the enzyme
(subject to C6). The products are then released one at a time in subsequent
binding/release steps (per C2).

For example, `E_ATP_HCO3 → Estar_ADP_Pi` is an isomerization that produces
ADP and Pi simultaneously on the enzyme. ADP and Pi are then released in
separate steps: `Estar_Pi + ADP ⇌ Estar_ADP_Pi` and
`Estar + Pi ⇌ Estar_Pi`.

*New capability. Current code only supports releasing 1 product per
isomerization (ping-pong) or all remaining products at once (final
isomerization).*

## Impact

| Reaction | Current (buggy) | New rules |
|---|---|---|
| Simple ter-ter | 3,969 | **283** (169 seq + 114 pp) |
| Pyruvate carboxylase | ~4,000 | **312** (169 seq + 143 pp) |
| Pyruvate dehydrogenase | ~4,000 | **334** (169 seq + 165 pp) |

The current code produces 3,969 = 63 × 63 due to the arbitrary subset
combining bug (C1). With the weak ordering fix alone, the sequential
3→3 isomerization produces 13 × 13 = 169 topologies. The new constraints
add ping-pong mechanisms on top.

Breakdown for simple ter-ter (283 topologies):
- 1 sequential iso pattern (E_A_B_D → E_P_Q_R):
  13 × 13 = 169 topologies (weak orderings of 3 substrates × 3 products)
- 12 bi-uni iso patterns × 9 topologies = 108
  (9 = 3 substrate orderings × 3 product orderings per half-reaction,
  where each half-reaction has 2 metabolites → 3 weak orderings each)
- 6 hexa-uni iso patterns × 1 topology = 6
  (each half-reaction has 1 metabolite → 1 weak ordering)

These three groups don't overlap — they have different numbers of
isomerization steps (1, 2, and 3 respectively).

## Constraints removed from current code

- **Arbitrary subset combining** (bug fix, C1 supersedes): current lines
  546-580 taking all 2^n-1 subsets of paths and unioning step sets. Produces
  biochemically meaningless ordering variants. Replaced by weak ordering
  enumeration in C1.
- **Ping-pong/sequential partition**: current lines 527-541 partitioning
  paths into ping-pong and sequential groups for separate combining.
  Superseded by C1 — combining is based on shared isomerization steps, not
  on presence of Estar forms.
- **Non-empty residual requirement** (C4 supersedes): current line 343
  `isempty(residual) && continue` is removed.
- **Ping-pong requires remaining substrates**: current line 333
  `!isempty(remaining_subs)` check on ping-pong branch. Superseded by C7 —
  paths that ping-pong with no remaining substrates dead-end naturally.
- **Substrate names in isomerization forms**: current `_form_name` calls
  passing `on_enzyme_subs` to isomerization form construction. Superseded
  by C8.

## Bystander mechanisms — not needed at init level

A bystander is a substrate/product that binds to the enzyme during a
half-reaction where it does not participate in the isomerization, creating
parallel catalytic paths.

Analysis showed that when bystander binding uses the same kinetic constants
as catalytic binding (K_D_bystander = K_D_catalytic) and the parallel
isomerization has the same rate constants, the bystander terms cancel in
the King-Altman derivation — the rate equation is identical to the base
mechanism. Bystander only produces a different rate equation when kinetic
constants are allowed to differ from the base.

Therefore bystander expansion is NOT a separate pipeline stage. If needed,
it would be an `expand_mechanisms` move that adds bystander binding with
relaxed constraints (+N parameters). This is deferred — the catalytic
topology constraints in this spec are sufficient for the initial
enumeration.
