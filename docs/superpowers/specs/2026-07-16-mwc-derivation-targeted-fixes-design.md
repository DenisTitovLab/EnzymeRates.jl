# Targeted MWC derivation fixes (solve-then-limit declined)

**Date:** 2026-07-16
**Status:** Design, approved. Supersedes
`2026-07-16-mwc-solve-then-limit-derivation-design.md` and its plan, which are
**declined** — see "Why the refactor is declined" and the measurements below.

## Summary

The ping-pong `:OnlyA` bug documented in
`docs/superpowers/findings/2026-07-16-pingpong-onlya-kcat-bug.md` (112 of 691 PFK-P
base-tier mechanisms) is real and reproduced. It is the *only* surviving correctness
argument for the solve-then-limit rewrite. Both of its failure classes have local
fixes. This design ships those fixes plus three independently-valuable changes, and
does not rewrite the derivation.

## Why the refactor is declined

Every load-bearing claim in the solve-then-limit design was reproduced
independently. Results:

| Claim | Measured |
|---|---|
| Per-state normalization is the "load-bearing crux" (Task 5) | **Algebraic no-op** — identical to the shipped cross-weight, term for term; 0 to 5.1e-16 over 7 mechanisms across all branches, including `CatN=4` and regulator multiplicity 4 |
| Delete-before-solve caused "every allosteric bug this year" | **1-for-4 by name, 1-for-2 by root cause.** The L-term leak and the `_kcat_forward` crash were normalization/basis defects; the refactor *keeps* the normalization |
| `:OnlyA` limit ≢ deletion (so the rewrite changes something) | **Limit ≡ deletion** on every constructable `:OnlyA` mechanism; converges as 1/X to 1.6e-13. The rewrite is behaviour-preserving there |
| Guard boundary is "precisely" the merged guard | **False from ter up.** 0 over-rejection, 1074 under-rejection over 17,814 tag assignments |
| The complete guard is a new discovery | **Already documented** in `_onlya_haldane_violation`'s docstring (`thermodynamic_constr_for_rate_eq_derivation.jl:421-429`), *including* the argument that graph deletion makes the gap benign |
| Removes the three-way `d_free` branch | **False.** Gate4 (`:NonequalAI`, no `:OnlyA`, no deletion) requires cross-weight. The branch survives; the net diff is likely an addition |
| Fitted-param counts unchanged | **False** — 4→6, 8→11, 9→12. Reconciling them requires an indep filter, i.e. deletion re-implemented one level up |
| `:OnlyA` limit "drops vanishing monomials" | **Not a monomial filter.** The solve makes some `k_I` *dependent*, so it is an order-of-vanishing analysis over the dep DAG |
| `^n` cross term at `N_I ≠ 0` may be broken | **Correct.** Exact at n=2 *and* n=3 against an independent concerted-MWC oracle (3.1e-7, converging as 1/FAST) |

Two further findings:

- **The proposed guard rule is probably mis-specified.** The spec's rule is "apply the
  `:OnlyA` limits to the solved relations and require the infinite exponents to
  balance." But for the witness class the solve *pins* `K_I` to a finite dependent
  expression (ε-exponent `w = 0`), so there is no infinite exponent and the rule sees
  a clean equality and admits it — missing exactly the cases it advertises catching.
- **The perf contract was never measured** by any part of the investigation, and Task 5
  rewrites the compiled expression. It is the plan's own hard gate.

Because `:OnlyA` ⟹ `:OnlyA` catalysis ⟹ a dead inactive cycle (measured: 54/54 tag
combos over the uni-uni skeleton, 0 exceptions), the `:OnlyA` limit only ever kills an
already-dead cycle. Solve-then-limit therefore cannot change any shipped `:OnlyA`
equation. Its remaining value was the ping-pong bug, which this design fixes locally.

The declined spec and plan stay in the repo as a record. They are not deleted.

## The bug, reproduced

`docs/superpowers/prototypes/2026-07-16-mwc-solve-then-limit/` holds the original
prototypes. The reproduction below is fresh (scratchpad `pp_repro.jl`), on the PFK-P
shape, all four mechanisms guard-**admitted** (they are valid — the guard is right):

```
E + ATP ⇌ E(ATP)                                                        :: EqualAI
E(ATP) <--> E(F16BP; residual = ATP - F16BP)                            :: tag3
E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP)  :: EqualAI
E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP)      :: tag6
E(F6P; residual = ATP - F16BP) ⇌ E(ADP)                                 :: tag4
E + ADP ⇌ E(ADP)                                                        :: EqualAI
```

| case | tags | `d_free_I` | `rate_equation` | `_kcat_forward` |
|---|---|---|---|---|
| control | all `:EqualAI` | `1` | 0.2231 | ok |
| V-system | tag3=`:OnlyA` | `1` | 0.0970 | ok |
| **err1** | tag3,tag6=`:OnlyA` | **`0`** | **NaN** | **crash** |
| **err2** | tag4,tag6=`:OnlyA` | `F16BP·k/K` | 0.0498 (finite) | **crash** |

**Root cause (err1).** `_reachable_from_free` (`rate_eq_derivation.jl:1173-1193`) seeds
reachability from *every* empty-`bound` form. A ping-pong covalent intermediate has an
empty `bound` but a non-empty `residual`, so it seeds itself as a free root. The
`:OnlyA`-deleted I-graph then keeps a covalent island that free `E` cannot reach, no
spanning tree rooted at free `E` exists, and `D[g_free] = 0`. The normalization divides
by it → `NaN` structurally, at every concentration.

The seeding is deliberate (its docstring says so) but wrong for the I-state. Under
formulation 1 **only free enzyme flips**, so the inactive conformation is entered *only*
via free `E`; a component free `E` cannot reach holds no inactive mass and must be
stranded. `_reachable_from_free` is used only by `_state_allo_mechanism`'s I-state
pruning, so tightening the seed cannot affect the active state.

**Root cause (err2).** `d_free_I` carries a product. The cross-weight multiplies the
A-numerator by `d_free_I^n`, so every saturating-substrate group key acquires a product
factor; `_kcat_forward` evaluates at products = 0 and filters product-bearing keys, so
`a_keys` empties and it throws. `rate_equation` is unaffected — the factor cancels
between numerator and denominator, which the measured finite 0.0498 confirms. This is a
**group-key** defect in one function, not a derivation defect. `_kcat_forward`'s own
comment (`:960-965`) already states the normalization "is a common factor of the
saturating-limit ratio, so it leaves kcat's value unchanged" — it is applied there only
to match `rate_equation`'s branch.

## Deliverables

### 1. err1 — tighten the I-state free-root seed

Seed `_reachable_from_free` from forms with an empty `bound` **and** an empty
`residual`. One predicate.

**Measured on this change (spike, reverted):** err1 `d_free_I` 0→1, `rate_equation`
NaN→0.0966, `_kcat_forward` crash→ok; controls unchanged; **1803/1803 derivation tests
and all 12 ground-truth gates pass**, golden reference byte-identical, perf gate green.

The docstring must be rewritten: it currently states the opposite rationale.

### 2. err2 — make `_kcat_forward` group on un-normalized polynomials

The `d_free` factor is a common factor of each conformation's saturating-limit ratio,
so it cannot change kcat's value — but it pollutes the group key. Compute the kcat
group keys without it.

**Open sub-decision, resolved at implementation time by the tests:** the normalization
does *not* cancel across the `L`-weighted A/I combine, only within a conformation. If
grouping on un-normalized polys breaks the cross-conformation pattern match, the
fallback is to strip only the metabolite part of the common factor from the key while
keeping the value path unchanged. The acceptance test is the same either way: err2's
`_kcat_forward` returns finite, and `kcat consistent with rate_equation` stays green.

If neither approach yields a clean local fix, **stop and report** — that is the one
finding that would legitimately reopen the rewrite, and it must not be worked around.

### 3. Complete the `:OnlyA` guard (Stiemke)

Replace the per-row sign test in `_onlya_haldane_violation` with exact ε-feasibility:
reject unless some `w > 0` satisfies `M·w = 0`, where `M` is the constraint matrix
restricted to the ε-normalized `:OnlyA` columns. Exact rational Fourier–Motzkin; a
working reference implementation exists (scratchpad `dirB_lib.jl`).

**Measured:** 0 over-rejection over 17,814 tag assignments (uni-uni, ordered/random
bi-uni, ping-pong, random-order bi-bi, ter-uni cube); agrees with the sign test on all
1042 assignments up to bi-bi; 1074 under-rejections at ter. Enumeration-reachable:
**136 of 13,005** mechanisms at BFS depth 4 from `init_mechanisms(A+B+C→P)`, first
appearing at depth 3 with only **6 groups**.

**Honest value.** These 1074 are exact I-state duplicates of feasible assignments, so
today they derive *working* equations. This buys deduplication and semantic honesty —
not a wrong number fixed. It is cheap (constructor-only, zero `rate_equation` runtime
impact) and costs no valid mechanism, which is why it ships.

Constructor-only: no `rate_equation` codegen change, so the perf contract is untouched.

### 4. Oracle gates

- **n=2/n=3 concerted MWC, `N_I ≠ 0`.** Port the validated oracle. This closes the
  largest validation hole in the derivation: n≥2 has **no** mass-action oracle for any
  family today. Self-validating (L=0 reduction, all-`:EqualAI` L-independence, `v=0` at
  equilibrium, and n=1 reduction to the existing free-flip reference — the last is what
  proves it is formulation 1 and not formulation 2). Discriminating: catches
  `Q_I^n`-for-`Q_I^{n-1}` at 96%, `L^n` at 10%, a dropped cross term at 21%.
- **Ping-pong value gate.** Gate at **1e-15**, not 1e-4: the fast-flip limit is *exact*
  for this topology (the free-enzyme flip is the only cut between the two conformation
  subnetworks, so it carries zero net flux at any FAST), and shipped achieves 4.25e-16.
  Delete the DEFERRED comment at `allosteric_ground_truth.jl:596-600`.
- **err1/err2 regression gates**, from deliverables 1 and 2.

**Formulation-1 warning, to be recorded in the harness:** a multi-protomer oracle that
lets *ligated* oligomers flip is formulation 2 and disagrees by 0.1–3% for live
`:NonequalAI` — a real number that is not a bug. The three usual self-checks do not
distinguish the two formulations; only the n=1 reduction to `biuni_nonequalAI_freeflip_flux`
does.

### 5. Record the decision

Mark the solve-then-limit spec and plan declined, with a pointer to this document, so
the misattributions are not re-litigated. Correct in place: the "12–43%" figure (measured
0.08–86%, median 28%, n=400; the dividing line is metabolite-in-`D`, not cross-weight),
and Task 4's impossible witness (a lone `:OnlyA` edge provably cannot be one — the
minimum is 7 of 12 cube edges).

## Architecture

No architectural change. Deliverable 1 touches one predicate in
`_reachable_from_free`; 2 touches `_kcat_forward`'s group-key construction; 3 replaces
the body of `_onlya_haldane_violation`; 4 is test-only; 5 is docs-only. The
single-conformation King–Altman, the combine, the normalization, and the `@generated`
codegen are untouched.

## Testing

Every deliverable is TDD'd: failing test, confirm the failure, minimal fix, confirm
green.

Non-negotiable gates, verified directly (not via a subagent's claim):

- `rate_equation` **0 allocations, < 120 ns** (`test_rate_equation_performance`).
  Deliverables 1-3 change no codegen, so this must stay green *unchanged*.
- The **golden reference** (`test/reference/allosteric_golden_reference.txt`) must stay
  **byte-identical**. Deliverable 1 only alters I-state pruning for mechanisms whose
  I-graph has an unreachable residual island — no current spec has one (measured: 1803/1803
  pass, golden green). If any golden block moves, that is a finding to explain, not to
  regenerate.
- All 12 `allosteric_ground_truth.jl` gates.
- Full `Pkg.test()` before the branch is declared done.

## Risks

- **Deliverable 2 is the only genuinely open one.** Its resolution is gated on tests,
  and its failure mode is "stop and report", not "work around".
- **Deliverable 3 changes what the enumerator admits.** 0 over-rejection is measured, so
  no valid mechanism is lost, but 136-per-13,005 fewer mechanisms enter the search.
  Enumeration-count tests may move; each move must be explained by the guard, not
  accepted blindly. `test_types.jl:1779` asserts `uni(:OnlyA, :NonequalAI, :EqualAI)` is
  a violation — flagged as a possible pre-existing over-rejection (a `:NonequalAI` iso
  group's free `k_I` ratio can absorb the imbalance, but the guard's `keep` filter drops
  only `:OnlyA` iso groups). If the exact test flips it, that is a finding for Denis, not
  a test to edit.
- **Fourier–Motzkin blowup** on large mechanisms. The reference implementation errors at
  >4000 rows. It ran the full 16,384-assignment ter-uni sweep inside budget, but the
  constructor is on the enumeration hot path; if a real mechanism blows up, fall back to
  the sign test for that mechanism (sound, just incomplete) rather than erroring.
