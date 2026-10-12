# The July solve-then-limit MWC design and the :EqualAI thread prototype

Repository root: `/home/denis.linux/.julia/dev/EnzymeRates`. All paths below are relative to it unless absolute.

**What I read:**
- Every file added by b067aa0, 651a821, 4504f84, 9c95f12, 6f0011e and 938c502, read with `git show`.
- The thread prototype and its report in `.superpowers/harness/record/`. The `.patch` file is byte-identical to the patch embedded in the report.
- `test/allosteric_ground_truth.jl`.
- `src/rate_eq_derivation.jl:200-298` and `:540-1168`.
- `src/thermodynamic_constr_for_rate_eq_derivation.jl`.
- `ROADMAP.md`, plus commits 9ffba17 (#71) and 37f91d6.

I ran no Julia. The prototype outputs were never committed, so every number below is quoted from the spec, the findings or the reports. Anything marked **Inference** is my own reasoning and has not been run.

## Bottom line

1. **The July design had no scaling parameter and no leading-order extraction.** Its "limit" is a pass that "drops vanishing monomials" after per-state normalization (spec :108-109; plan :172). The requirement to take a leading-order term, and the 0/0 warning, first appear in the ROADMAP commit 37f91d6 (`ROADMAP.md:54-56`).
2. **The prototypes are numeric only, and the headline "6/6 uni-uni combos match" compares deletion with deletion.**
   - Both the derivation stand-in (`solve_then_limit.jl:76-78`) and the oracle (`brief_all_constructable.jl:25-31`) encode `:OnlyA` by removing edges.
   - True limits (BIG up to 1e12, BigFloat) appear only in the LDH, uni-uni and ping-pong normalization scripts.
3. **Round 2 retracted "delete the d_free cross-weighting".** The shipped cross-weight is per-state normalization, and it is required. Even so, `ROADMAP.md:60-62` still says the design "would delete the `d_free` cross-weighting", and it cites only b067aa0 and 651a821, not 4504f84.
4. **The design was declined.** The 9ffba17 (#71) message says: "The declined solve-then-limit rewrite and superseded design docs were removed." #71 shipped targeted fixes instead. The plan's limit-encoded oracle (Task 1, steps 1-3) was never built. The ping-pong value gate and the n = 1-3 oligomer gate from the same task did land, in #71.
5. **The thread prototype was never merged.** It fixes the inactive state's second name for `:EqualAI` constants: mechanisms with undeclared fitted names drop from 898 to 0. On 7.10 Denis chose the single-equation derivation instead (37f91d6).
6. **Biggest risks for the redesign:** deletion and the limit disagree; the limit depends on which constants are held fixed; the guard and the ROADMAP treat a steady-state (SS) `:OnlyA` binding differently; and ping-pong has two free forms. See section 7.

## 1. The design (math)

**Model: formulation 1.** Only free enzyme changes conformation. Each conformation runs its own cycle, and L = [E_I]/[E_A] for free enzyme. The per-site rate is (spec :38-42):

`v = E_t·(N_A·D_A^{n−1} + L·N_I·D_I^{n−1}) / (D_A^n + L·D_I^n)`

Here N and D are each conformation's King–Altman numerator and partition, divided by that conformation's own free-enzyme weight `d_free`. Regulator factors W multiply each branch (`rate_eq_derivation.jl:1144-1150`; `docs/src/deriving/mwc_allostery.md:116-120`).

**The four steps (spec :28-46):**
1. Derive both conformations on the full, undeleted graph. Every constant gets its own symbol, except that `:EqualAI` constants are shared.
2. Solve Haldane and Wegscheider on the stacked A+I system with one shared Keq. The code already stacks the two systems (`rate_eq_derivation.jl:963-975`; "cross-state ties emerge as I-row − A-row").
3. Normalize each state and combine them plainly.
4. Apply `:OnlyA` as a limit to the combined equation: `K_I→∞` for a rapid-equilibrium (RE) binding, `k_I→0` for an SS step, then drop vanishing monomials. In plan :172 this reads "drop every monomial carrying an `:OnlyA` `K_I`-inverse (RE) or `:OnlyA` `k_I` (SS)". That is exactly the substitution the ROADMAP now warns against.

**Partner limits** have no algorithm in the design.
- They are meant to be implicit: the solve runs first, so the dependent constants are monomials in the independent ones and inherit the limit. The spec says "the constraint has already fixed the surviving relationships" (:45-46).
- The proposed guard applies the limits to the solved relations. It rejects when the infinite exponents do not balance, for example "K = ∞ forced onto an `:EqualAI`/shared step". The spec calls this "the Stiemke feasibility condition" (:148-163).
- That check now exists in the constructor (from #71), at `thermodynamic_constr_for_rate_eq_derivation.jl:402-418` and `:474-509`:
  - Each cycle row of the cycle basis C gets the ε-exponents of the `:OnlyA` groups, with `K_I = K_A/ε`.
  - It runs a per-row sign test, then a strictly-positive nullspace test.
  - It drops `:OnlyA` chemistry groups from the check graph, because their free `k_I` ratio absorbs any imbalance (:443-448).
  - It scales SS bindings exactly like RE bindings, so `kon/koff→0` (:465-470).
- The cycle basis is the integer nullspace of the form×step incidence. Each cycle is classified as Haldane (net change c·ν_net) or Wegscheider (c = 0) (`thermodynamic_constr_for_rate_eq_derivation.jl:171-221`).

**0/0** is handled only by ordering: normalize each conformation before the limit, because the limit "is clean only on a dimensionless denominator" (spec :53-58).
- This fails when `d_free` itself vanishes in the limit (section 7).
- The shipped code avoids 0/0 by pruning instead. `_inactive_groups` keeps the forms reachable from residual-free unbound forms over all non-`:OnlyA` steps (`rate_eq_derivation.jl:779-813`). An empty inactive graph returns (0, 1, 1) (:235-237).

**SS limit convention** (spec :69-74; `ldh_koff_ambiguity.jl`, `eq_consistency.jl`). For an SS `:OnlyA` substrate binding, only "both k→0" (koff_I→0 together with kon_I) keeps v = 0 at equilibrium for every L. Setting koff_I = koff_A, or koff_I→∞, violates the second law. When the vanishing step is not a cut, this limit equals deletion.

**L and normalization** (spec :117-146):
- **Metabolite-free `d_free`:** divide by it. Skipping the division changes only L's value.
- **Metabolite-bearing `d_free`** (an SS substrate binding, or ping-pong): dividing would put a concentration in a denominator. The polynomial form is the cross-weight `den = d_I^n·Q_A^n·W_A + L·d_A^n·Q_I^n·W_I`, which is exactly what ships (`_mwc_state_polys` :1102-1115, `_mwc_cross_weight` :1058-1059, `_num_den_exprs` :1152-1166).
- Skipping the normalization gives errors of 12–43%, depending on concentration.
- A dimensionless "clean L" form (kcat/Km-style factoring) was deferred to a second PR as cosmetic.

## 2. Spec review rounds

**Round 1 (651a821, titled "Spec review round 1 (Denis)"):**
- **Guard:** replace "reuse the #70 sign check" with applying the limits to the solved full-system constraints. This is complete for ter-substrate and up; #70's sign test stays as an enumeration pre-filter.
- **Oracle:** encode `:OnlyA` as the K→∞ / k→0 limit of the full network, because a deletion-encoded oracle is "circular for the trap".
- **Code deletion** became a primary goal, with a delete-list and a LOC metric.
- **Normalization** was split into a metabolite-free case (trivial, deferrable) and a metabolite-bearing case, with an open question about whether deferring it keeps v correct.

**Round 2 (4504f84).** Its body says "Prototype corrections (verified this session)"; it does not attribute these changes to Denis.
- The shipped derivation is correct for every valid mechanism tested, including the `:NonequalAI` ping-pong (3e-16). The motivation was reframed as robustness and clarity, "not a bug fix".
- "Cross-weighting is 8–37% off" was a strawman: `attack.jl:82` compared the combine against `N_A/(D_A(1+L))`. Normalization stays, re-expressed per state.
- Normalization is required for metabolite-bearing weights.
- For a valid S:`:OnlyA` mechanism the SS limit coincides with deletion, so the deletion oracle was not circular there.
- The deletion goal was downgraded to "moderate".

**Text that round 2 left stale:**
- The Migration section still says "cross-weight-free form" (spec :184-185).
- The Evidence section still calls `attack_claimB_strong.jl` "the decisive … no-cross-weighting result" (:257-258).
- The README for the prototypes (:14) and `ROADMAP.md:61-62` repeat the retracted claim.

## 3. Plan (9c95f12): eight tasks, gate first

1. **Oracle.**
   - Limit encoding for the uni, multi and metabolite-D gates (`K_I = 1e10×K`, SS `k_I = 0`).
   - Port the ping-pong `oracle_Eflip` as a value gate.
   - Add an n = 2 gate with a live inactive numerator (N_I ≠ 0).
2. **Inactive state on the full graph.** Test: the I-state form count equals the A-state form count.
3. **Stacked solve on the full graph.** Test: the I reverse rate references the I forward rate and Keq; `_dep_graph_is_sound`.
4. **Limit-the-constraints guard.** It must catch a ter-uni case the sign test admits.
5. **The crux.** Per-state normalization (divide on the fast path, cross-weight when metabolite-bearing), plain combine, then the limit: "iterate against the oracle … If a case cannot be made to match, STOP".
6. Recompute `_kcat_forward`.
7. Delete the graph-deletion code and report the LOC change.
8. Regenerate the goldens; fitted counts must not change.

Throughout: the perf gate runs on every task, and the non-allosteric path must stay byte-identical.

## 4. Prototype results

| Script | Compares | Result (as quoted) | Caveat |
|---|---|---|---|
| `solve_then_limit.jl`, `brief_all_constructable.jl` | v_new vs a full-network oracle, 6 constructable uni-uni combos | ≤ 4.8e-7; v = 0 at equilibrium (spec :50-52) | Both sides delete edges; FAST = 1e7 sets the floor |
| `n2_spotcheck.jl` | Family A combine with exponent n vs an explicit 2-protomer network | 3.1e-7 | Only N_I = 0. P-bound T forms also flip, so formulation 1 vs 2 is not tested |
| `attack_claimB_strong.jl` | Per-state normalized combine vs free-flip oracle; bi-uni `:NonequalAI` with metabolite-bearing D | 9.6e-9 | Normalization done numerically (free E from a linear solve) |
| `attack.jl` | "cross-weight" | 8–37% off | Bogus formula; retracted |
| `brief_confirm_trap.jl`, `brief_enum_check.jl` | Trap (S:EqualAI, cat:NonequalAI, P:OnlyA) | v_new diverges; EP_I holds mass; the constructor rejects all four combos that v_new gets wrong | The oracle deletes P binding in I |
| `redesign_norm_vs_skip.jl`, `redesign_pingpong.jl`, `pp_eqcheck.jl`, `final_check.jl` | Normalized vs raw combine against LDH limit oracles (BIG = 1e4…1e12) and ping-pong E-flip oracles | Raw: 12–43% wrong; normalized: matches | BigFloat numerics |
| `ldh_koff_ambiguity.jl`, `eq_consistency.jl`, `limit_vs_deletion.jl` | SS `:OnlyA` convention | Only the "both k→0" convention keeps v = 0 at equilibrium | The I-state keeps A's catalysis (`eq_consistency.jl:12-22`) |
| `compare_pingpong.jl` | Shipped `rate_equation`, ping-pong `:NonequalAI` | 3e-16 | Now a gate at rtol 1e-10 (`allosteric_ground_truth.jl:523-611`) |
| `RC_gt.jl`, `RC_probe.jl`, `RC_dfree.jl` | PFK-P ping-pong `:OnlyA` | v = 0 at equilibrium (1.7e-16, 7.6e-17); `d_free_I` = 0 or product-bearing | `RC_gt` strands the covalent island by construction (:5-9) |

**The ping-pong `:OnlyA` `_kcat_forward` bug** (findings doc, 6f0011e and 938c502):
- **Scope:** 112 of 691 PFK-P base-tier mechanisms. All were ping-pong, and all were accepted by the guard: 14 V-system, 14 balanced K-system, 84 one-sided plus chemistry; 56 of the 84 return NaN.
- **Root cause:** ping-pong has two empty-bound forms, and `_reachable_from_free` seeded both as free enzyme.
  - **err1:** `d_free_I = 0`, so `rate_equation` returns NaN.
  - **err2:** `d_free_I` carries a product. The cross-weight injects that product into every saturating key, `_kcat_keys` ends up empty, and it throws "no kcat components" (now `rate_eq_derivation.jl:616-622`).
- **Options listed:** a sentinel value, a derivation-time degeneracy check, or the redesign.
- **What #71 shipped:**
  - Free-enzyme roots now come only from forms with no bound metabolite and no residual (`rate_eq_derivation.jl:767-796`).
  - An err1 gate (`allosteric_ground_truth.jl:749-784`).
  - A "dead-inactive" enumeration rule, so the crash class can no longer be enumerated.
- **Untested:** no gate covers a hand-written err2 shape (only chem2 `:OnlyA`). Whether it still throws is unverified.

## 5. Open problems, and why it was not executed

**Open at the time:**
- Symbolic decoupled normalization (spec :245-248).
- The semantics of the limit pass.
- A live inactive numerator (N_I ≠ 0) at n ≥ 2.
- Where the guard should run (plan :217).

**Why it was not executed:**
- Round 2 removed the correctness case for the rewrite.
- The spec names the central risk itself: "Rewriting a currently-correct subsystem … If the rewrite cannot beat the shipped equations on the harness, it is not worth landing" (:238-244).
- The next day, #71 did four things: fixed err1, made the guard complete with a Stiemke test, added the dead-inactive enumeration rule, and added the n = 1-3 and ping-pong oracles.
- The message then records the decline (9ffba17). The reason for declining is not written down; I infer it was that no correctness payoff remained.

**Why it was revived (37f91d6, 8.10).** The motivation is now parameter-count correctness. Two known issues are both caused by deriving the inactive state on its own pruned graph with its own Wegscheider rename (`ROADMAP.md:10-34`):
- 766 mechanisms give a shared `:EqualAI` constant a second fitted name.
- 256 fold a tied pair in opposite directions, and 512 more fold in one state only.

## 6. The :EqualAI thread prototype

**What it changes** (src +106/−85, mostly docstrings):
- `_state_parts` anchors each inactive-state Parameter on the step as `am` stores it (`s in stored ? s : _reversed(s)`). As a result, `name(p, am)` renders the group representative of the whole mechanism.
- A new keyword `naming` (default `mech`) is threaded through `_raw_symbolic_rate_polys`, `_build_wegscheider_rename_map`, `_assemble_constraints` and `_param_columns`.
  - The graph still comes from the pruned `cm`: forms, segments, `i_free`, and the free-enzyme set for pivot priority.
  - Names come from `am`.
- A new `_anchored_steps(step_params)`. `_thermodynamic_constraints(rxn, flat)` and `_compute_alpha(mech, flat, …)` now read step direction from the anchored steps, which fixes the `F ⇌ E` orientation flip.
- `_raw_symbolic_rate_polys` loses its `rename_map` positional argument.

**New tests:**
- T1 and T2: mass-action gates, +127 lines in the ground-truth file.
- T3, T3b and T4: in the derivation tests.
- T5: an aggregate pin requiring all 977 `_expand_to_allosteric` children of the bi-bi seeds to fit only declared constants.
- T7: a split-move test (one child; parameter count 6 → 7).

**Measured:**
- **Before/after the fix:** at HEAD, 7 tests fail and 2 error, exactly as predicted; with the fix, 41/41 pass.
- **Undeclared symbols** in inactive-state polynomials or fitted lists: 898 mechanisms → 0. Of the 15,422 mechanisms common to both runs:
  - 14,524 are unchanged.
  - Class (a), 766 mechanisms: drop 1-3 fitted parameters (578 drop one, 172 drop two, 16 drop three; 970 names in total, 626 `K_` and 344 `k_`).
  - Class (b), 132 mechanisms: only a `:NonequalAI` inactive name changes.
- **Numeric equivalence:** 2.0e-16 on 60 sampled mechanisms when each phantom name takes its whole-mechanism value. In 4 of 40, the phantom was the inactive copy of a Haldane-dependent reverse rate, which HEAD fitted freely.
- **Direction gate:** reading direction from the pruned graph fails the ground truth 5/5 (0.0663 vs 0.1135), while `fitted_params` still passes.
- **Suites:** the full suite passes 32,797/32,797 in 6m05s, and the perf gate held.
- **Harness:** 3 of 15 files differ. Two equation records moved; +72 bi-bi allosteric children; +75 ping-pong depth-2 children. In the ping-pong case, 139 parents that had phantoms gained 147 children, because the split gain test had been comparing against an inflated parent count.

**What it leaves alone:**
- The opposite-fold hazard (256 + 512; rank equals the fitted count, so these are identifiable but larger models).
- The `i_free` ambiguity.
- The fact that names and directions now come from two step lists.

**Status:** not on main (no `_anchored_steps` there). `ROADMAP.md:23-24` calls it "a narrower fix, prototyped and measured".

## 7. Lessons, reusable pieces, risks

### Reusable

- **The mass-action oracle.**
  - `_testhelper_mwc_ground_truth_flux` (`allosteric_ground_truth.jl:16-27`).
  - The n-protomer oracle with its formulation-1 pin (:797-880).
  - The exact ping-pong gate, with FAST-invariance as its justification (:523-611).
  - The RE segment holding both free E and F (:652-694), the inert-inactive gate (:696-747) and the dead-inactive ping-pong gate (:981-1077).
  - What is missing is any limit-encoded `:OnlyA` network. The current `:OnlyA` gates still delete edges and add extra EP flips (:37-50).
- **A ten-line numeric spanning-tree oracle:** `ND(forms, edges, cat_edges, free)` (`redesign_norm_vs_skip.jl:14-20`). It computes each form's Volkenstein–Goldstein D_r as a principal minor, plus `d_free`, in BigFloat. With the BIG sweep (:56-59) it is a ready differential oracle for a symbolic leading-order engine.
- **Constraint machinery:** the cycle space (`thermodynamic_constr_for_rate_eq_derivation.jl:171-221`); the ε-exponent Stiemke test (:402-509); and the stacked A+I solve with inactive-before-active pivot priority (`rate_eq_derivation.jl:953-975`).
  - A design that only relabels constants still needs that stacked solve to find forced A = I ties (the "collapse" example in `mwc_allostery.md:157-185`).
- **Direction-invariant numerator weights ω_e** (`rate_eq_derivation.jl:276-294`). This already satisfies goal 9 for the numerator.
- **From the thread prototype:** tests T1, T2 and T5, and the "undeclared fitted names" scan as a metric for "one symbol per constant".

### Risks

1. **Deletion differs from "both k→0" when the vanishing step is a cut** (Inference).
   - Toy: for `E ⇌ X` with both rates scaled by ε, `[X]/[E] = k_f/k_r` for every ε > 0. Deletion gives 0.
   - The leading-order term is therefore the vanishing step's own ratio. That ratio belongs to no cycle, so it would stay in the rate law.
   - Family A, with SS S binding and chemistry both `:OnlyA` (the LDH gate, `allosteric_ground_truth.jl:217-273`), cuts off {E(S), E(S,B)} in exactly this way.
   - Deletion is recovered only if association vanishes faster than dissociation, i.e. K→∞, which is the guard's semantics. The July claim "coincides with deletion" was tested only where the step was not a cut.
2. **An RE binding with K→∞ is not deletion** (Inference).
   - With `K_I_P→∞` and `k_I_r` Haldane-dependent, `k_I_r·[EP_I] = k_I·[E_I]·P/(K_S·Keq)`. `K_I_P` cancels exactly, EP_I holds no mass, and the inactive cycle stays live.
   - So for the July "trap" there are three different answers: the prototype gives N_I = 0, the deletion oracle traps mass in EP_I, and the true limit keeps a live cycle. The "pendant EP_I mass" is a deletion artifact.
   - A 15 July memory note says the same: "Deleting the step manufactures a sink the real limit never has" (`~/.claude/projects/-home-denis-linux--julia-dev-EnzymeRates/memory/archive/project_kcat_no_components_onlya_tags.md`).
   - The guard rejects this combination by convention, because the partner limit would land on a `:NonequalAI` constant.
3. **Partner limits depend on which constants are held fixed** (Inference). For S:`:OnlyA` with `:NonequalAI` catalysis, fixing `k_I` gives N_I → 0, while fixing `k_I/K_I_S` keeps the inactive cycle live. Pivot priority currently chooses the dependent constant, so it also chooses the limit. That choice needs an explicit rule.
4. **The guard and the ROADMAP disagree on SS `:OnlyA` bindings.**
   - The guard treats an SS binding as K→∞ and rejects it next to non-`:OnlyA` chemistry.
   - `ROADMAP.md:51-52` ("both `k → 0` for a steady-state step") and the July consistency test (`eq_consistency.jl`, whose inactive state keeps A's catalysis) admit it.
   - Combined with risk 1, my suggestion is K→∞ for every binding and both k→0 only for chemistry. That is a decision for Denis.
5. **Dual free forms.** Today residual-bearing forms are never roots (`rate_eq_derivation.jl:767-796`), and the first residual-free form in species order carries L (:267).
   - (Inference) For err1 under a true K→∞ limit, the covalent island sits on the free side of the vanishing F6P binding. Its inactive weight then grows like [ADP]·K_I/[F6P] instead of vanishing, and the limit is non-uniform as ADP → 0.
   - A per-state check that the order of `d_free` equals the order of D is the precise form of the findings' option 2.
6. **Normalization.** Division is a Laurent operation only when `d_free` is a metabolite-free monomial (`rate_eq_derivation.jl:1069-1077`). Taking the limit by substitution on the cross-weighted polynomial gives 0/0 whenever `d_free → 0`. A dead inactive cycle (handled today by the branch at :1157-1163) becomes N_I = O(ε) under limits.
7. **Naming and orientation from a subgraph.** Taking names or directions from the pruned graph silently produces wrong rate laws. The thread prototype's direction gate shows this.
8. **Formulation.** Per-state normalization is formulation 1. A per-form flip changes `:NonequalAI` laws by 0.1–3% (`allosteric_ground_truth.jl:786-793`), and the raw combine misses ping-pong by 0.95–93% (:531-535).
9. **kcat.** `_kcat_forward` is a saturation limit implemented by grouping monomials (`rate_eq_derivation.jl:557-726`), and it broke on the same normalization. One leading-order engine could serve both the `:OnlyA` limits and kcat.

### Toy prototypes to run first

- (a) Symbolic leading-order extraction on the six uni-uni combos, checked against a limit-encoded oracle with a BIG sweep.
- (b) A Family A cut toy, to decide the SS binding convention.
- (c) err1 under a true K→∞ limit.
- (d) Relabel-only versus the stacked solve on the "collapse" example.

### Files

- `/home/denis.linux/.julia/dev/EnzymeRates/.superpowers/harness/record/equalai-thread-prototype.patch`
- `/home/denis.linux/.julia/dev/EnzymeRates/.superpowers/harness/record/equalai-prototype-report.md`
- `/home/denis.linux/.julia/dev/EnzymeRates/test/allosteric_ground_truth.jl`
- `/home/denis.linux/.julia/dev/EnzymeRates/src/rate_eq_derivation.jl`
- `/home/denis.linux/.julia/dev/EnzymeRates/src/thermodynamic_constr_for_rate_eq_derivation.jl`
- `/home/denis.linux/.julia/dev/EnzymeRates/ROADMAP.md`

The July spec, plan, findings and prototype scripts exist only on branch `mwc-solve-then-limit` (read them with `git show`):
- `docs/superpowers/specs/2026-07-16-mwc-solve-then-limit-derivation-design.md`
- `docs/superpowers/plans/2026-07-16-mwc-solve-then-limit-derivation.md`
- `docs/superpowers/findings/2026-07-16-pingpong-onlya-kcat-bug.md`
- `docs/superpowers/prototypes/2026-07-16-mwc-solve-then-limit/`