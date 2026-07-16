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
conformation and dividing by `[E_free]`. The design must do this **symbolically**,
and the symbolic operation splits by whether the free-enzyme spanning-tree weight
`D[g_free]` carries a metabolite:

- **`D[g_free]` metabolite-free** (a constant or pure rate-constant expression —
  every rapid-equilibrium-binding mechanism). Trivial: divide `N`/`D` by it, a
  common factor that `L` and the rate constants absorb. Deferring it changes only
  `L`'s numerical value, not `v`.
- **`D[g_free]` metabolite-bearing** (a steady-state substrate binding, as in the
  LDH i-state; ping-pong bi-bi). Dividing would put a concentration in the
  denominator (forbidden). Instead **factor** the denominator into a dimensionless
  polynomial — terms like `1`, `S/Km`, with non-negative exponents and `Km` a
  function of the `K`s and `k`s — an operation in the spirit of `_kcat_forward`
  (which already infers a saturating-limit `kcat`). This is the delicate part.
  **Open question (being prototyped): whether skipping this factoring keeps `v`
  correct.** For a metabolite-bearing weight the A and I states normalize by
  different concentration-dependent factors, so `L` (a constant) may not absorb
  the difference — meaning `v`, not just `L`, could be wrong. The prototype
  validated the *normalized* combine, never the un-normalized one, so this is
  unresolved. If `v` does stay correct, the factoring is a pure second-PR
  refinement (clean `L`); if not, it is required in the first PR for LDH-class
  mechanisms.

**Sequencing (pending the ping-pong result).** If deferral preserves `v`: PR 1
delivers the correct rate equation with the metabolite-free normalization only,
and PR 2 makes every denominator dimensionless (the clean-`L` refinement, package
-wide). If deferral does not preserve `v` for metabolite-bearing weights, the
factoring lands in PR 1 for those mechanisms. `ping-pong bi-bi` is the decisive
test case either way.

## The guard: limit the solved constraints

The consistency check falls out of solve-then-limit for free. After step 2 we
hold the full two-conformation thermodynamic relations — each an equality among
log-parameters (a Haldane or Wegscheider row). **Apply the `:OnlyA` limits to
those solved relations and require each to stay a consistent equality.** A valid
mechanism's `K_I → ∞` / `k_I → 0` exponents cancel and the relation reduces to a
finite equality (often `0 = 0`); an impossible one produces a contradiction — a
lone `∞`, or `K = ∞` forced onto an `:EqualAI`/shared step. Reject on any
contradiction.

This is **complete**, unlike PR #70's per-row sign test: it operates on the
actual solved constraints (all cycles at once), so it catches the ter-and-up
multi-cycle inconsistencies the sign heuristic admits. It is the Stiemke
feasibility condition expressed as "the infinite exponents balance." The same
rule catches a `:NonequalAI` that thermodynamics forces to `K_I = K_A`.

**Relationship to PR #70.** #70's `_onlya_haldane_violation` (the sign heuristic,
in the `AllostericMechanism` constructor) stays as a cheap early reject during
enumeration — it never wrongly rejects a valid mechanism, only admits some benign
ter+ cases. This design's limit-the-constraints check is the *complete*
correctness guarantee, run at derivation time. Whether the constructor keeps the
cheap heuristic or adopts the full check is a design sub-decision; the derivation
must carry the complete check regardless. Implementation lands on top of #70
(branch off `main` once #70 merges, or off `onlya-haldane-validator` meanwhile).

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
that caught every prior bug — is the acceptance gate, with **one change to its
methodology**: it currently encodes `:OnlyA` by *deleting* the inactive edge, the
same move this derivation abandons. Validating a limit-derivation against a
graph-deletion oracle is circular on exactly the cases that matter (the trap).
Rebuild the oracle to encode `:OnlyA` as the **limit** — the full two-conformation
network with `K_I` large / `k_I` small — so it is a genuinely independent check.
(The prototype's adversary already did this to surface the trap.)

Each derived equation must match the oracle (rtol 1e-4) and give `v = 0` at the
equilibrium metabolite ratio, across: the existing gates (uni/multi-`:OnlyA`,
metabolite-bearing-`D`, `:NonequalAI`), the TR mechanism, ping-pong bi-bi (the
normalization test case), and a **new** `N_I ≠ 0` multi-protomer gate for the
`^n` cross term. The `rate_equation` performance contract (0 allocations,
< 120 ns) must hold — verified by the existing perf gate.

## Code deletion is a primary goal

Simpler here means less code. Solve-then-limit exists to *remove* machinery, and
the diff should be net-negative in `src/rate_eq_derivation.jl`. Targeted for
deletion:

- `_state_allo_mechanism`'s graph pruning (`:1214`).
- the `d_free` cross-weighting in `_allosteric_num_den_exprs` (`:966-978`,
  `:1691+`) and in `_kcat_forward`.
- its helpers `_is_metabolite_free_monomial`, `_invert_monomial`, and the
  three-way rendering branch.
- `d_free` surfacing from `_state_rate_polys` (`:331`, `:381`, `:400`).

A refactor that adds more than it removes has missed the point; the LOC delta is
a review metric.

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
