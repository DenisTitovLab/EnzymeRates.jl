# Solve-then-limit MWC allosteric derivation

## Motivation

The allosteric MWC derivation builds each conformation's rate polynomial by
**deleting** graph elements for `:OnlyA` tags before the thermodynamic
constraints are solved. Deleting before solving discards the constraints that
tie the deleted step to the rest of the cycle, and every allosteric derivation
bug this year traces to it: the L-term leak, the `_kcat_forward` "no components"
crash, the trap where a deleted inactive form is re-populated by reverse
catalysis (giving a metabolite-bearing free-enzyme weight and a wrong equation),
and the thermodynamic inconsistency that PR #70 now rejects. Each fix has patched
a symptom; the `d_free` cross-weighting machinery exists *only* to paper over the
graph fragmentation the deletion causes.

The single-conformation King–Altman derivation is mature and unaffected. This
design reuses it for both conformations and combines them, so consistency holds
by construction.

## The approach (validated by prototype)

Replace graph-deletion with **solve then limit**:

1. **Derive both conformations in full.** The active (A) and inactive (I)
   conformations are complete King–Altman networks with the same topology. Each
   step carries its own parameter, except `:EqualAI` steps where the I parameter
   equals the A parameter (one shared symbol). `:NonequalAI` and `:OnlyA` steps
   get separate `K_A`/`K_I` (or `k_A`/`k_I`).
2. **Solve the thermodynamic constraints on the full coupled system.** Each
   conformation obeys its own Haldane/Wegscheider relations; both share one
   `Keq`.
3. **Normalize each conformation to its own free-enzyme weight,** then combine:
   `v = E_total · (N_A + L·N_I) / (D_A + L·D_I)` at `n = 1`, and
   `(N_A·D_A^{n-1} + L·N_I·D_I^{n-1}) / (D_A^n + L·D_I^n)` for `n` protomers,
   where `N_state`/`D_state` are the single-conformation numerator/denominator
   **normalized so the free-enzyme term is 1**.
4. **Apply `:OnlyA` as a limit on the combined equation** — `K_I → ∞` drops an
   `:OnlyA` RE binding's terms; `k_I → 0` drops an `:OnlyA` SS step's flux. The
   limit is taken *after* the constraint solve, so the constraint has already
   fixed the surviving relationships.

### Prototype evidence (n = 1 uni-uni, adversarially attacked)

- **Correctness:** matches the mass-action ground truth on all 6 constructable
  tag combos (≤ 4.8e-7); `v = 0` exact at equilibrium. n = 2 Family A matches an
  explicit 2-protomer network (3.1e-7).
- **Simplicity:** in the metabolite-bearing-`D` regime the cross-weighting was
  built for (LDH i-state, `D_A ≠ D_I`), the **plain** combine matches to 9.6e-9
  while a **cross-weighted** denominator is 8–37% off. Cross-weighting is not
  merely unnecessary — it is *wrong* there. The `d_free` machinery goes away.
- **Guard boundary:** the one break the attack found (an `:OnlyA` product binding
  with active inactive-catalysis, where a pendant `EP_I` node holds mass the
  limit drops) is **structurally unreachable** — PR #70's `_onlya_haldane_violation`
  rejects exactly those mechanisms at construction. The correctness boundary of
  solve-then-limit is *precisely* the merged guard.

## Architecture

All changes are in `src/rate_eq_derivation.jl` (the `@generated` allosteric path).
The single-conformation King–Altman (`_raw_symbolic_rate_polys` and its callees)
and PR #70's constructor guard (`_onlya_haldane_violation`) are reused unchanged.

**Removed:**

- `_state_allo_mechanism`'s graph deletion (`rate_eq_derivation.jl:1214`) — the
  I-state is no longer a mutilated subgraph.
- The `d_free` cross-weighting in `_allosteric_num_den_exprs` (`:966-978`,
  `:1691+`) and `_kcat_forward`, and its helpers `_is_metabolite_free_monomial`,
  `_invert_monomial`, and the three-way rendering branch.
- `_state_rate_polys`'s third return value `d_free` (`:381`, `:400`).

**Added / changed:**

- **Per-conformation derivation.** `_state_rate_polys` derives each
  conformation's full King–Altman `(N, D)` with conformation-scoped parameters
  (A-names, I-names, shared `:EqualAI` names), on the **undeleted** graph.
- **Per-state free-enzyme normalization** (the load-bearing new step, symbolic —
  see below). Normalize each `(N, D)` so the free-enzyme spanning-tree weight is
  `1` before combining.
- **The combine.** `_allosteric_num_den_exprs` becomes the plain MWC weighting of
  the two normalized `(N, D)` pairs, raised to `n`.
- **The `:OnlyA` limit.** A pass over the combined polynomial that sends
  `:OnlyA` `K_I → ∞` / `k_I → 0` and drops vanishing monomials.
- **The constraint solve.** `_combined_state_dependent_exprs` already stacks the
  A- and I-state constraint rows; it stays, now operating on the full (undeleted)
  two-state system. `:OnlyA` and `:NonequalAI` I-parameters are genuine unknowns
  the solve resolves, rather than deletions.
- **`_kcat_forward`.** Recomputed from the normalized combine — the saturating
  limit of the new equation — without cross-weighting.

## The load-bearing piece: symbolic per-state normalization

The prototype proved the *target* (plain combine of per-state-normalized `N`/`D`
reproduces the ground truth) but produced `N`/`D` by numerically solving each
conformation and dividing by `[E_free]`. The design must do this **symbolically**:
derive each conformation's King–Altman and divide by its free-enzyme spanning-tree
weight `D[g_free]` before combining. The current code already surfaces
`D[g_free]` (that is what `d_free` was); the change is to **normalize each state
by it up front** rather than cross-weight the two states in coupled form
afterward. Working out the exact symbolic operation — and confirming it leaves
`rate_equation` allocation-free — is the primary implementation task and is
deferred to the plan.

## The guard

**Depends on PR #70.** The correctness boundary of this derivation is PR #70's
`_onlya_haldane_violation`, so implementation must land on top of it — branch off
`main` once #70 merges (or off `onlya-haldane-validator` meanwhile). The design is
written against a `main` that contains #70.

No new guard. `_onlya_haldane_violation`, in the `AllostericMechanism`
constructor, is the correctness boundary and runs upstream of every derivation.
The collapse-detection floated during brainstorming is **narrower** (it misses
the pendant-node trap) and is not built. The one requirement: the new derivation
must remain downstream of the constructor check, so no invalid mechanism reaches
it. If any future path constructs a derivation input without going through the
constructor, it must carry this gate.

## `n > 1` and the `N_I ≠ 0` cross term

The `^n` combine is validated only for Family A (`N_I = 0`, inactive numerator
vanishes) at `n = 2`. The cross term `L·N_I·D_I^{n-1}` — an active inactive
conformation (`:NonequalAI` catalysis) at `n ≥ 2` — is untested, and `n ≥ 3` is
untested entirely. The plan must validate the `N_I ≠ 0` branch against a
multi-protomer ground truth before the derivation is trusted for it.

## Migration

Every allosteric rate equation rederives (to the correct, cross-weight-free
form), so `test/reference/allosteric_golden_reference.txt` regenerates. The
non-allosteric path is untouched and must stay byte-identical. Fitted-parameter
counts should be unchanged (the `:OnlyA` limit drops the same parameters the old
deletion did); any change is a finding to explain, not to accept.

## Testing

`test/allosteric_ground_truth.jl` — the n = 1 two-conformation mass-action solver
that caught every prior bug — is the acceptance gate. Each derived equation must
match it (rtol 1e-4) and give `v = 0` at the equilibrium metabolite ratio, across:
the existing gates (uni/multi-`:OnlyA`, metabolite-bearing-`D`, `:NonequalAI`),
the TR mechanism, and a **new** `N_I ≠ 0` multi-protomer gate for the `^n` cross
term. The `rate_equation` performance contract (0 allocations, < 120 ns) must
hold — verified by the existing perf gate.

## Risks and open items

- **Symbolic normalization** is proven correct as a target but not yet built; it
  is the crux and could surface complications (deferred to the plan).
- **`N_I ≠ 0` at `n ≥ 2`** is unvalidated.
- **Scope:** this rewrites a perf-critical `@generated` subsystem that has been
  iterated on repeatedly. The mitigation is the ground-truth harness as a
  hard gate and PR #70's guard as a fixed correctness boundary.
- **What this does not change:** the enumeration moves, the DSL, the
  single-conformation derivation, and the constructor guard all stay.

## Evidence

Prototype scripts under the session scratchpad: `solve_then_limit.jl`,
`n2_spotcheck.jl`, `attack.jl`, `attack_claimB_strong.jl` (the decisive
metabolite-bearing-`D` no-cross-weighting result), `brief_confirm_trap.jl`,
`brief_enum_check.jl` (constructor rejects the trap), `brief_all_constructable.jl`
(correct on exactly the constructable set).
