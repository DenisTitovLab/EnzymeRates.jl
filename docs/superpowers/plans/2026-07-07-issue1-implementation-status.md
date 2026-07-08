# ABOUTME: Status of the allosteric combined-constraint-solve implementation —
# ABOUTME: what shipped, what works uncommitted, and the one remaining derivation gap.

# Issue-1 implementation status (branch `allosteric-combined-constraint-solve`)

Autonomous execution of `docs/superpowers/plans/2026-07-07-allosteric-combined-constraint-solve.md`.
Stopped at the plan's explicit stop-condition: a genuine derivation gap that needs a design
decision from Denis before proceeding. **No wrong derivation was committed.**

## Committed (green base, in git log)

- **Spec + plan** (`docs/superpowers/specs|plans/2026-07-07-allosteric-combined-constraint-solve*`).
- **Issues 2-4 findings** (`docs/superpowers/specs/2026-07-07-ldh-issues-2-4-findings.md`).
- **Phase 0** (`6d734af`, `94fa644`): three failing reproducers captured as fixtures + a
  `_dep_graph_is_sound` regression + a `rate_equation`-callable regression. These are
  **RED on this branch until the fix lands** — the intended TDD state.
- **Phase 1** (`ee7292b`): behaviour-preserving factoring of `_dependent_param_exprs_kernel`
  into `_assemble_constraints` + `_solve_dependent_set`. Verified: Enzyme Derivation Tests
  1802/1802, split_resolution + allosteric_collapse green. **Safe, green, keep this.**

## Uncommitted in the working tree (Phase 2 — works, one gap)

`git diff` on `src/rate_eq_derivation.jl` and `test/mechanism_definitions_for_test_enzyme_derivation.jl`:

- `_combined_state_dependent_exprs`: stacks the A-state and I-state constraint rows over one
  combined tagged-column space and solves once via `_solve_dependent_set`.
- `_dependent_param_exprs(::AllostericEnzymeMechanism)`: combined solve + the verbatim
  regulator/`L` handling.
- `_build_dep_assignments`: consumes the combined dep, splits A/I by `_i_state_symbol_set`.
- `rate_equation_string(::AllostericEnzymeMechanism, ::ReducedMode)` (line ~1807): fixed to
  render the combined-solve assignments (it was rendering the A-state-native `dep_A`, which
  diverged from the combined partition and left `koff_A_Lactate_E` undefined).
- Two LDH i-state spec counts re-baselined 7→8 (see "what works").

### What works (verified)

- **All three reproducers fixed:** graph-sound 15/15, `rate_equation` callable 3/3
  (`koff_Pyruvate_ENAD`, `K_NAD_ELactate`, `kon_I_NAD_E` all return finite values).
- **Enzyme Derivation Tests 1803/1803** — the QSSA oracle, equilibrium (`v=0` at `Q=Keq`),
  analytical-rate, and 0-alloc/120 ns performance all green on every allosteric spec.
- **Honorable splits are correct:** "LDH i-state NonequalAI 5-group" now has 8 fitted params,
  and the identifiability rank is 8 — the old count of 7 was *under*-parameterised (the S_I
  gate dropped a genuinely free DOF). This is the case Denis's "prefer `K_I` fitted" rule is
  about, and the combined solve gets it right.

### The gap (why I stopped)

The combined solve **loses the honorable-vs-forbidden split distinction** that
`_split_resolution` computed. For a *forbidden* `:NonequalAI` split — e.g. the single S-binding
in the `test_allosteric_collapse` uni-uni case — the correct behaviour is `K_I_S_E = K_A_S_E`
(collapse; the split has no thermodynamic freedom). The combined solve instead keeps `K_I_S_E`
**free**, re-introducing exactly the kind of non-identifiable phantom the #58/#61 work removed.

Root cause: my `[A-rows; I-rows]` formulation relies on the difference `I-row − A-row` to
supply the affinity-split constraint `C_A·δ = 0` (δ = `log K_I − log K_A`). That works when the
I-state cycle is **live**, but when it is **dead** (pruned / no steady-state cut) there is no
I-row to subtract, so the split constraint is missing and the split stays free.

Confirmed: the erroring "Allosteric edge cases" test and the collapse assertions in
`test/test_allosteric_collapse.jl` are the symptom; `_dependent_param_exprs` returns `K_I_S_E`
in `indep` where it must be a dependent mirror.

## The fix (elegant, but needs your call)

Add the affinity-split constraints as **explicit rows** in the combined matrix: for each
A-cycle, a row `C_A·δ = 0` over the `:NonequalAI` affinity columns (`+coeff` on `K_I`, `-coeff`
on `K_A`; for a steady-state binding the affinity is `koff/kon`, so the row decomposes as
`+koff_I −kon_I −koff_A +kon_I` scaled by the cycle coefficient). Then:

- The single RREF **automatically** classifies honorable vs forbidden: an honorable split lies
  in the nullspace of the split constraints and stays free; a forbidden split is pinned and
  collapses — **no separate `_split_resolution` classification needed.**
- This unifies everything into one solve and still deletes `_split_resolution`'s *solve*
  (the RREF partition), reusing only its constraint *construction* (the `C_A`-on-`:NonequalAI`
  matrix, `rate_eq_derivation.jl:1304-1321`, plus the affinity/speed decomposition at
  `:1264-1282`).

**Decision needed from you:** confirm that forbidden `:NonequalAI` splits should collapse
(`K_I = K_A`), i.e. your "prefer `K_I` fitted" rule applies only to *honorable* splits. The
identifiability data supports this: honorable splits are fully identifiable (keep them fitted);
forbidden splits are non-identifiable phantoms (collapse them). If you agree, the fix above is
mechanical; if you actually want all `:NonequalAI` splits free, then the current uncommitted
Phase 2 is already correct and only the re-baselines (edge-case params + `test_allosteric_collapse`
assertions + golden) need updating.

## Update — box/split constraints attempted (uncommitted, promising but not exact)

Per Denis's decision, I implemented the box/split constraint rows in the combined solve
(all uncommitted in `src/rate_eq_derivation.jl` and `src/thermodynamic_constr…jl`):

- `_split_constraint_rows(am, col_index, n, indep_A)`: for each A-cycle and each *collapsible*
  `:NonequalAI` group (all its A-state rate constants independent), a row `C_A·δα = 0`,
  `δα = log(effK_I/effK_A)` (Kd for RE, koff/kon for SS), RHS 0. Stacked into the combined
  matrix; I-column pivot priority bumped (+1) so a pinned split collapses `K_I` onto `K_A`.
- `_solve_dependent_set`: `best_pri` init changed from `-1` to `typemin(Int)` so a
  binding-only constraint (which a box row is) can pivot a binding column — catalytic columns
  still dominate when present, so non-allosteric output is unchanged (verified indirectly).

**Verified working:** all three reproducers fixed (callable, graph-sound); the forbidden-split
edge case now collapses `K_I_S_E = K_A_S_E` with detailed balance restored (`|v| ≈ 2e-15` at
equilibrium, was 3.7); the LDH i-state specs match HEAD (5-group → 7). **Zero oracle,
equilibrium, or callable failures across the whole derivation file.**

**Not exact:** 12 `test_constraint_counting` assertions fail — ~6 allosteric specs (e.g. `m_all`:
indep 12→11) get a *different partition* than the tested `_split_resolution`. The equations stay
thermodynamically consistent, but the box construction does not reproduce `_split_resolution`'s
collapsible/SS-affinity-decomposition/multi-split-nullspace logic exactly, and the `best_pri`
change reaches beyond the box rows. Making it exact is essentially re-implementing
`_split_resolution`.

**Two ways to finish (Denis's call):**
1. **Refine the box rows** to match `_split_resolution` precisely — port its `collapsible`
   affinity/speed decomposition (`_split_resolution` :1264-1301) and confine the `best_pri`
   relaxation to the box rows only. Highest-fidelity to the "one solve" goal; most fiddly.
2. **Hybrid (recommended for correctness/risk):** keep the combined A/I solve for the catalytic
   partition, then reuse the *tested* `_collapse_mirror_exprs(am)` to move each forbidden split's
   `K_I` from `indep` to `dep`. Exact by construction (matches HEAD), low-risk; keeps
   `_split_resolution`/`_collapse_mirror_exprs` (less net deletion, but they compute a genuine
   thing — the box collapse — not merge patchwork). Revert the `best_pri` change and box rows.

## Update 2 — box refuted, hybrid attempted, both blocked on split reconciliation

**Empirical double-check (Denis's ask): is box more correct than `_split_resolution`? No.**
For the six specs where box's partition deviates, BOTH box and `_split_resolution` (HEAD) hold
detailed balance, but the identifiability rank shows HEAD is closer to the true dimension on all
six: on four specs HEAD is `count = rank` (fully identifiable) while box is one lower — box collapses
a *genuinely identifiable* differential-affinity/cooperativity DOF. So `_split_resolution` is more
correct; box over-collapses honorable splits (it pins any split with a cycle incidence, missing the
multiplicity/nullspace structure `_split_resolution` computes). Evidence: `scratchpad/rank6.jl`,
`thoro_eq.jl`.

**Hybrid attempted (combined solve + `_collapse_mirror_exprs`), UNCOMMITTED in the tree — INCOMPLETE,
DO NOT COMMIT AS-IS.** It fixes the MWC honorable splits (match HEAD) and gets the derivation file to
1800/1802, but it is **thermodynamically wrong for reproducers [2] and [3]**: they are callable yet
violate detailed balance (|v| ≈ 0.14, 0.08 at `Q = Keq`; `scratchpad/repro_eq.jl`). Root cause: for a
forbidden `:NonequalAI` split whose **I-cycle is live**, the combined solve's I-row already pins one
side of the split (e.g. `koff_A = f(koff_I)`), while `_collapse_mirror_exprs` pins the other
(`koff_I = f(koff_A)`) — the two disagree. Applying the mirror is circular ([3]) or mis-ordered ([2]);
skipping it (the `_expr_references_any` decoupling filter now in the tree) leaves the forbidden split
free → detailed balance broken. Neither is correct.

**The real blocker:** correctly partitioning allosteric affinity splits (honorable-free vs
forbidden-collapsed, over multiplicity-weighted nullspaces) is exactly `_split_resolution`'s job, and it
is entangled with the old per-state merge. Re-deriving it (box) over-collapses; bolting it onto the
combined solve (hybrid) conflicts with the I-row pinning.

**Most promising next direction (untried):** make the combined solve pin the *I-side* of a live-cycle
split (bump I-column pivot priority in `_combined_state_dependent_exprs` — NOT the box rows, NOT the
`best_pri` change), so its choice agrees with the mirror direction (`K_I`/`koff_I` dependent); then
`_collapse_mirror_exprs` only needs to handle *dead*-I-cycle forbidden splits, and a topological sort
in `_build_dep_assignments` covers the remaining mirror→dep ordering. If that reconciles [2]/[3]'s
detailed balance while keeping MWC honorable splits free, the hybrid is complete.

**Recommendation:** given how entangled the split logic is, the lowest-risk path to a *correct* Issue-1
fix may be to keep the tested `_split_resolution`/`_collapse_mirror_exprs` verbatim and fix ONLY the
three merge/S_I defects surgically (the original targeted-guard option from the spec, which a subagent
verified fixes 20/23), rather than the full combined-solve rewrite. The rewrite is elegant but the
affinity-split reconciliation is the hard 20%.

## Remaining tasks after the gap is resolved

- Task 2.3: delete the now-dead chain — `_i_state_referenced_syms`, `_collapse_mirror_exprs`,
  `_split_resolution`/`SplitResolution`, `_state_dependent_exprs`, `_state_raw_dependent_exprs`
  (keep `_rref_partition`/`_integer_nullspace` — used by `_thermodynamic_constraints`; keep the
  `_state_mechanism`/`_state_step_params`/`_state_all_params`/`_state_rate_polys` family).
- Phase 3.1: delete `test/test_split_resolution.jl` + its `runtests.jl` include; update the
  mirror-text assertions in `test/test_allosteric_collapse.jl` (keep every `veq` equilibrium
  assertion).
- Phase 3.2: regenerate `test/reference/allosteric_golden_reference.txt`.
- Phase 3.3: full suite green; confirm the net `src/` diff is strongly negative in logic lines.

## How to resume

Everything above is uncommitted in the working tree. Read this file, run `git diff`, decide the
forbidden-split semantics, then either add the split-constraint rows (recommended) or update the
re-baselines. The reproducer regressions (`test/reference/allosteric_undefvar_reproducers.jl`) are
the fast red→green signal; `julia --project=<testenv>` derivation runs use the test env at
`…/scratchpad/testenv` (has `OrdinaryDiffEqFIRK`).
