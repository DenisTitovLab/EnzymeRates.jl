# Solve-then-limit MWC allosteric derivation

> **Status: DECLINED (2026-07-16).** Superseded by
> `docs/superpowers/specs/2026-07-16-mwc-derivation-targeted-fixes-design.md`,
> which reproduces every claim below and records what measured true. Kept as a
> record. In particular: the normalization rewrite (Task 5) is an algebraic
> no-op; the `:OnlyA` limit equals graph deletion on every constructable
> mechanism; the guard gap is real but already documented in
> `_onlya_haldane_violation`'s docstring, and graph deletion is what makes it
> benign; the three-way `d_free` branch is required by `:NonequalAI` and does
> not go away; and the `^n` cross term is already correct at n=2 and n=3. The
> "12-43%" figure measured 0.08-86% (median 28%, n=400), and the dividing line
> is metabolite-in-`D`, not the cross-weight branch. Task 4's named witness
> (a lone `:OnlyA` edge) is provably impossible — the minimum is 7 of 12 cube
> edges — and the guard rule it proposes admits the witness class anyway,
> because the solve pins `K_I` finite rather than infinite.

## Motivation

**This is a simplicity and robustness rewrite, not a bug fix.** The shipped
allosteric derivation — post PR #70 (which rejects the thermodynamically
impossible mechanisms) — is *correct* for every valid mechanism tested, including
the `:NonequalAI` ping-pong case the ground-truth harness had left deferred
(shipped `rate_equation` matches an independent limit oracle to 3e-16). The goal
here is not to fix a wrong equation.

The goal is to remove a *structural* fragility. The derivation builds each
conformation's rate polynomial by **deleting** graph elements for `:OnlyA` tags
*before* the thermodynamic constraints are solved. Deleting before solving can
silently drop a constraint tying the deleted step to the rest of the cycle — and
that is the mechanism behind every allosteric derivation bug this year (the L-term
leak, the `_kcat_forward` "no components" crash, the trap, the thermodynamic
inconsistency). Each was patched after it appeared. Solving the constraints on the
full system *before* taking the `:OnlyA` limit makes the whole class
unreachable: consistency holds by construction, so a novel tag mix or topology is
far less likely to produce a new derivation bug. That structural guarantee — not
a smaller diff — is the return.

The single-conformation King–Altman derivation is mature and reused unchanged.

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

### Prototype evidence (n = 1 uni-uni + ping-pong + LDH, adversarially attacked)

- **Correctness:** matches the mass-action ground truth on all 6 constructable
  tag combos (≤ 4.8e-7); `v = 0` exact at equilibrium. n = 2 Family A matches an
  explicit 2-protomer network (3.1e-7).
- **Normalization is required, not optional, and is the same operation the shipped
  code already does correctly.** The `:OnlyA` limit is clean only on a
  *dimensionless* denominator, so each conformation must be normalized to its
  free-enzyme weight *before* the limit. For a metabolite-bearing weight (LDH
  steady-state substrate binding; ping-pong), skipping normalization makes `v`
  concentration-dependently wrong by 12–43% — not a relabeling of `L`. Per-state
  normalization is algebraically **identical** to the shipped `d_free`
  cross-weighting (multiply each state by the other's weight vs divide each by its
  own); the rewrite re-expresses the same correct operation in decoupled form, it
  does not delete it. (An earlier prototype's "cross-weight 8–37% off" compared the
  normalized combine to a *bogus* formula, not the shipped cross-weight — corrected.)
- **Guard boundary:** the one break the attack found (an `:OnlyA` product binding
  with active inactive-catalysis, where a pendant `EP_I` node holds mass the
  limit drops) is **structurally unreachable** — PR #70's `_onlya_haldane_violation`
  rejects exactly those mechanisms at construction. The correctness boundary of
  solve-then-limit is *precisely* the merged guard.
- **`:OnlyA` SS limit convention:** for a valid `S:OnlyA` mechanism only
  `k_I → 0` (inactive cycle dead, `N_I = 0`) is thermodynamically legal — the
  other conventions violate `v = 0` at equilibrium. In that regime the limit
  *coincides* with graph deletion, which is why the shipped LDH equation is
  already correct. The design's payoff is ping-pong / `:NonequalAI` (no deletion
  involved), trap rejection, and structural robustness — not fixing LDH.

## Architecture

All changes are in `src/rate_eq_derivation.jl` (the `@generated` allosteric path).
The single-conformation King–Altman (`_raw_symbolic_rate_polys` and its callees)
and PR #70's constructor guard (`_onlya_haldane_violation`) are reused unchanged.

**Removed:**

- `_state_allo_mechanism`'s graph deletion (`rate_eq_derivation.jl:1214`) — the
  I-state is no longer a mutilated subgraph.
- the graph-fragmentation branch selection in `_allosteric_num_den_exprs`
  (`:966-978`, `:1691+`) — the `d_free_A == d_free_I` special case and the
  three-way rendering branch — since the undeleted I-state no longer fragments
  from `:OnlyA`.

**Kept (re-expressed), not removed:**

- the free-enzyme normalization itself. It is correct (identical to the shipped
  cross-weighting) and becomes decoupled per-state normalization. `d_free` is
  still computed per state; `_is_metabolite_free_monomial` / `_invert_monomial`
  may remain as the metabolite-free fast path.

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
  denominator (forbidden), so normalization is done in polynomial form — either
  the shipped cross-weighting (`den = d_free_I·Q_A + L·d_free_A·Q_I`, algebraically
  the per-state-normalized combine) or an equivalent `_kcat_forward`-style
  factoring of the denominator into a dimensionless polynomial (terms like `1`,
  `S/Km`, `Km` a function of the `K`s and `k`s). **Settled by prototype: this
  normalization is required — it cannot be deferred.** Skipping it makes `v`
  concentration-dependently wrong by 12–43% for these mechanisms (the required
  per-state factor is concentration-dependent, so no constant `L` absorbs it). The
  metabolite-bearing set includes LDH, the primary target, and ping-pong.

**Sequencing.** PR 1 ships the normalization for every mechanism (it is correct
today in cross-weight form; the rewrite carries it forward in decoupled per-state
form). Only the cosmetic "make every denominator dimensionless with a clean `L`"
refinement — which changes `L`'s numeric value but not `v`, and only for
metabolite-free mechanisms — is a legitimate PR 2. There is no correctness content
to defer.

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
that caught every prior bug — is the acceptance gate, strengthened by encoding
`:OnlyA` as the **limit** (the full two-conformation network with `K_I` large /
`k_I` small) rather than by deleting the inactive edge. For a *valid* `S:OnlyA`
mechanism the two coincide (the legal SS limit `k_I → 0` equals deletion), so this
is not a correctness fix to the existing gates; it matters because (a) it lets the
oracle exercise the trap and `:NonequalAI` cases where deletion and the limit
diverge, and (b) it keeps the oracle methodologically aligned with the derivation
it checks. (The prototype's adversary used exactly this form to surface the trap.)

Each derived equation must match the oracle (rtol 1e-4) and give `v = 0` at the
equilibrium metabolite ratio, across: the existing gates (uni/multi-`:OnlyA`,
metabolite-bearing-`D`, `:NonequalAI`), the TR mechanism, ping-pong bi-bi (the
normalization test case), and a **new** `N_I ≠ 0` multi-protomer gate for the
`^n` cross term. The `rate_equation` performance contract (0 allocations,
< 120 ns) must hold — verified by the existing perf gate.

## What is removed, and what is kept

Simpler here means less code, but expect a *moderate* net reduction, not a
dramatic one — the normalization is correct machinery that stays (re-expressed),
not buggy machinery deleted.

**Removed** (the structural-fragility surface):

- `_state_allo_mechanism`'s graph pruning (`:1214`) and the delete-then-solve
  ordering it serves.
- the graph-fragmentation handling the deletion forced — the `else`/cross-weight
  *branch selection* in `_allosteric_num_den_exprs` (`:966-978`, `:1691+`), the
  `d_free_A == d_free_I` special-casing, and the three-way rendering branch.
- PR #70's per-row sign guard, superseded by the complete limit-the-constraints
  check (or demoted to a cheap enumeration pre-filter).

**Kept, re-expressed** (correct, not deletable):

- the free-enzyme normalization. Today it is coupled cross-weighting
  (`d_free_I·Q_A + L·d_free_A·Q_I`); the rewrite makes it decoupled per-state
  normalization. Same math, clearer form. `_is_metabolite_free_monomial` /
  `_invert_monomial` may survive as the metabolite-free fast path.

The LOC delta is a review signal, but the primary metric is that the removed code
is exactly the constraint-losing surface — the diff should read as "delete
delete-then-solve, keep the correct normalization."

## Risks and open items

- **Rewriting a currently-correct subsystem.** This is the central risk and it is
  honest to name it: the shipped derivation is correct for every tested valid
  mechanism, so the rewrite trades a real (if bounded) regression risk for
  structural robustness and clarity. The mitigation is the ground-truth harness as
  a hard gate — every mechanism must match it — plus a byte-identical
  non-allosteric path and the unchanged perf gate. If the rewrite cannot beat the
  shipped equations on the harness, it is not worth landing.
- **Symbolic normalization** — the recipe is settled (per-state normalization,
  algebraically the shipped cross-weight), but carrying it into decoupled symbolic
  form while keeping `rate_equation` allocation-free is the main implementation
  task and may surface complications (the plan's first target).
- **`N_I ≠ 0` at `n ≥ 2`** (active inactive-catalysis, multiple protomers) is
  unvalidated; needs a multi-protomer gate.
- **What this does not change:** the enumeration moves, the DSL, the
  single-conformation derivation, and the constructor guard all stay.

## Evidence

Prototype scripts under the session scratchpad: `solve_then_limit.jl`,
`n2_spotcheck.jl`, `attack.jl`, `attack_claimB_strong.jl` (the decisive
metabolite-bearing-`D` no-cross-weighting result), `brief_confirm_trap.jl`,
`brief_enum_check.jl` (constructor rejects the trap), `brief_all_constructable.jl`
(correct on exactly the constructable set).
