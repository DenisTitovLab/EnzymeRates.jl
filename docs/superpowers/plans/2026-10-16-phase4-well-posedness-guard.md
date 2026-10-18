# Phase 4: well-posedness guard and admission — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reject every mechanism whose rate law is not the unique fast-RE limit of mass-action
kinetics where data sit, with a message that names the step, and drop such children in
enumeration; bring the flip at free E, the `:OnlyI` rule on modifier bindings and the tagging
move that keeps allosteric ping-pong reachable.

**Architecture:** One guard in `src/admission.jl` reads phase 3's reaction graph through masks
and flag overrides, never by constructing a mechanism. It builds two networks of one catalytic
subunit (the law network and the ε-network, flip at free E), enumerates the spec's faces, and
runs, in order: class 1 (no turnover cycle), class 2 (RE turnover on the ε-network), class 4's
interior admission, then per face class 3 (procedure A with a mod 2^61−1 limit test) and class 4's
face condition (B v3). It returns an isbits `_Verdict`; `_check_well_posed` throws its message.
Class 5 runs on the derived law. Test oracles in the shared fixtures file check every verdict.

**Tech Stack:** Julia 1.10–1.13; EnzymeRates internals after phase 3 (graph `Mechanism`, read API
of S3 §7); `Rational{BigInt}`, 256/512-bit `BigFloat` and `UInt64` arithmetic mod 2^61−1; Test,
TestEnv, Aqua, JET; Documenter with doctests.

**Spec:** `docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md` at 7205ba4
(source copy: `$S/graph_mechanism_wt/docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md`),
sections 7, 8, 13, 14, the phase-4 row of 15, 16, 17, 18 and Appendix A (B1, B2, B8, E2, E10).

Inputs every task reads (`$S` below):
- RES = `$S/phase45/research/phase4_research.md`: §1 rules in pseudo-code (the authority over
  research scripts), §2 rulings, §3 test plan, §4 task order, §5 Task 0 results.
- RUL = `$S/phase45/RULINGS.md` "Phase 4": P4-1 to P4-31. The amendments P4-8b, P4-9a/b/c,
  P4-11a, P4-12a, P4-19a, P4-24a, P4-29a/b override the research text they amend.
- S3D = `$S/phase3/proto/S3_design.md` (§1 types, §7 read API, §9 boundary function, §10 S3g-S3i);
  S3R = `$S/phase3/proto/S3_report.md` (S3d-S3f); R# = `$S/phase3/RULINGS.md` (R11, R18, R19, R27, R52).
- Research code under `$S/phase45/research/`; some scripts predate RES §1 (RECOVERY.md). Where a
  script differs from RES §1, follow RES §1 and say so in the task's commit body.
- Plan probes and their logs: `$S/phase45/plan/probes/` (the evidence before Tasks 4, 8
  and 11 cites them). They built on the S3f snapshot of phase 3's prototype, which no
  longer exists; their logs stay as evidence. The P4-32 refresh pass's probes at 08c4043,
  with their logs: `$S/phase45/refresh/`.
- Code cites name functions and testsets as they stand at 08c4043, the tip of the phase-3
  PR (#93) on graph-mechanism, where the P4-32 refresh pass checked them.

## Task index

Tasks run in order; each ends with the full suite green and pristine, and a commit. The last
column is the row of RES §4 that each task carries out.

| Task | Title | RES §4 |
|---|---|---|
| 0 | Confirmatory research (done 2026-10-14) | 0 |
| 1 | Test oracles and their self-tests | 1 |
| 2 | Networks, faces and classes 1-2, not wired | 2 |
| 3 | Class 3: procedure A, the limit test and its diagnosis, not wired | 3 |
| 4 | Admission interior | 4 |
| 5 | Admission's face condition B v3 | 5 |
| 6 | Class 5, the support test on the derived law, not wired | 6 |
| 7 | Trap specs become rejection tests (before the wiring, whose macros would throw on them) | 11 |
| 8 | Wire the guard | 7 |
| 9 | Flip at free E | 8 |
| 10 | `:OnlyI` on modifier-binding groups | 9 |
| 11 | The tagging move's `:pumped` drops through the oracles (P4-9c; outside the suite) | new |
| 12 | Tagging move (P4-1, P4-19) and change-tag through admission | 10 |
| 13 | Population tests | 12 |
| 14 | Docs, ROADMAP, the R18 scan and the draft PR | 13 |

Two departures from RES §4's order, which P4-31 says the plan follows: Task 7 (row 11) runs
before the wiring, whose macros would throw on the trap specs, and Task 11 is new (P4-9c).
Open for Denis A lists both.

The interface decisions that span tasks follow the File Structure. "Open for Denis", at the
end, lists the rulings the PR body names and every question the plan leaves open.

## Global Constraints

- Version stays 0.10.0 (`Project.toml` untouched); R19.
- Julia compat floor 1.10; CI runs Julia 1.13 with the newest SciML; on 1.13 a closure over a
  `@nospecialize` argument is untyped, so use `Base.Fix2` there.
- `rate_equation` stays allocation-free and under 120 ns per call for every
  `MECHANISM_TEST_SPECS` mechanism (`_testhelper_test_rate_equation_performance`: `allocs == 0`,
  `t < 120e-9`); the Expr-shape and flat-string pins in `test/test_rate_eq_derivation.jl` hold.
- Compile budgets in `test/test_compile_budget.jl` hold (init trace-compile ≤ 200, first
  `rate_equation` ≤ 6 s, Sig encoding ≤ 250 instances and 100 KB, fit compile ≤ 250 and 300 KB).
  Task 8, which puts the guard into `init_mechanisms`, `_minimal_flips` and the derivation
  entry, measures them before its full run; a breach stops the task for Denis, and no
  budget is raised unasked.
- R18: no function a caller loops over returns a `Union` of immutable reference-holding
  parametric structs; no loop reads a field whose type differs across a `Union`'s members. The
  guard returns the isbits `_Verdict`. Task 14 reruns the remark scan (`--pkgimages=no`) on
  derivation and identify workloads on Julia 1.12 and 1.13.
- Memory: one Julia process at a time (7.7 GB, no swap); the full suite runs with
  `--heap-size-hint=2500M` and needs about 4 GB free; each full run reports maxrss, and a rise
  over the branch point's figure is explained in the commit body; research runs stay ≤ 3 GB in
  chunks ≤ 20 min.
- Suite time: phase 4 adds at most 30 s of population tests to phase 2's 31 s (P4-28); ter-ter
  `init_mechanisms` stays under 150 s ("init_mechanisms on ter-ter within 150 s").
- R52's three figures at the branch point come from phase 3's final report (S3i's gates:
  full-suite maxrss at most 3.1 GB, ter-ter init at most 31 s, `test_mechanism_enumeration.jl`
  focused at most 145 s). Tasks 8 and 12 report all three beside the branch point's.
- Runner. Every command runs from the worktree `$S/phase45/guard_wt` (branch graph-guard, cut
  from graph-mechanism after the phase-3 PR), one Julia process at a time, through the session
  runner. Each shell call starts fresh, so every command block sets `S` first:
  `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh <log> <timeout_s> --project=. -e '<cmd>'`.
  Focused `<cmd>` (timeout 1800): `using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<file>.jl")`.
  Full `<cmd>` (timeout 3600): `using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])`.
  Logs go to `$S/phase45/logs/P4_T<n>_<what>.log`. A run longer than a few minutes goes to the
  background, and its log is polled in a bounded loop:
  `timeout 580 bash -c 'until grep -q "^exit " LOG; do sleep 20; done'; tail -40 LOG`.
- Measured figures. A `<t>`, `<kb>`, `<N>` or similar slot in an expected output or a commit
  body stands for the figure that step measures; the step names its source, and the commit
  carries the number. `CHILDREN` in Task 12's P6 testset is such a figure, read from Task 11's
  data before the testset is written.
- Commits: plain imperative summary, optional body, ending with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` and
  `Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE`. Never skip a hook.
- Draft PR only (graph-guard into graph-mechanism, stacked on the phase-3 PR); never merge. The PR
  body ends with the "Generated with Claude Code" line and no session link; `gh pr edit` fails on
  the Projects-classic deprecation, so edit through `gh api ... -X PATCH`.
- CLAUDE.md applies to every task: TDD (failing test first, watched failing), `_testhelper_`
  prefix on every test helper, no mocks, pristine output, ABOUTME headers on new files, 92-column
  lines, the three enumeration-test rules (inline macro fixture per testset, exact children with
  count and `Set`, chemistry written the enumerator's way), ROADMAP updated in the commit that
  fixes, adds or retires an entry.
- Fixtures are inline DSL as RES §3.2 writes them, with `::Inh` written `::Mod` and
  modifier-only names declared under `modifiers:`. A changed pin follows a value check (an oracle
  or an exact law) and its commit says why it changed.
- Messages follow RES §1.11: `<Mechanism|AllostericMechanism>: <cause naming each step by its
  DSL line and the face>. <fix>.`; a fix is printed only where a test shows its variant passes
  the guard (P4-22). Prose follows elements-of-style.

## Review Focus

Most likely to bite a user first; each line names the task whose tests pin it.

1. An admitted law evaluated where one substrate and every product are zero (a plotting grid):
   `rate_equation` returns NaN for the 357 census-type laws P4-12a lists, the ROADMAP known
   issue says so, and a fit never meets that corner (`FittingProblem` refuses a zero-rate row;
   the loss scores NaN as it scores 0). Test: Task 2.
2. A user pastes a step from an error message back into the DSL: every step line a message
   prints reads back as that edge, for `::Mod` copies, residual forms, fused and
   Theorell–Chance steps. Test: Task 2.
3. A mechanism phase 3 wrote (an identify CSV row or `show` text) that phase 4 rejects, such as
   an allosteric ping-pong with both halves `:OnlyA`: the parser raises the guard's message naming
   the step, not a parse error. Test: Task 8.
4. A hand-written `:OnlyA` group that mixes a plain and a fused binding of one metabolite whose
   limit exists (spec §7 "Mixed groups stay allowed"; the ROADMAP issue Task 12 retires): the
   guard admits it and the ε-path converges. Test: Task 4.
5. Renaming metabolites so the table order and the first substrate change: the guard returns the
   same fault, conformation and zeroed names (spec §14 "renames"; RES §3.3 P7 covers step order,
   direction and regrouping only). Test: Task 13.

---

## File Structure

Source (the guard runs where RES §1.10 places it):

| File | Change | Responsibility |
|---|---|---|
| `src/admission.jl` | modify (most code) | The guard: `_Fault`, `_LimitVerdict`, `_Verdict`, `_Network`, faces, classes 1-5, admission (interior and face), the limit test mod 2^61−1, messages (`_guard_message`, naming forms by `_form_text` and steps by `_step_text`), `_check_well_posed`, `_is_well_posed`, `_faces_pass`, `_assert_support`, beside `_derivation_input` (P4-23); `_assert_bottomed_segments` deleted (Task 8), `_assert_one_flip_segment` and `_free_forms_flip_together` deleted (Task 9) |
| `src/types.jl` | modify | Sub-graph masks next to `_inactive_mask` (`_excludes`, `_law_mask`; `_inactive_mask` and `_active_mask` on `_law_mask`, Tasks 9-10); `_outside_onlya` deleted (Task 9); the `:OnlyI` rule `_check_onlyi_groups`, with `_VALID_CAT_ALLO_STATES` and `_VALID_REG_ALLO_STATES` deleted for `_TAG_CODES` (Task 10); `_chemistry_groups(m)` (Task 12) |
| `src/thermo.jl` | modify | The joint solve's `:OnlyI` mirror (`_active_law`, `_in_active_law`, `_has_inactive_columns`, `_class`, `_ordered_columns`, `_kept_members`, `_dependent_param_exprs`; Task 10); the old `:OnlyA` guard deleted (`_onlya_haldane_violation`, `_onlya_violation_message` in Task 8; `_onlya_offenders`, `_onlya_cone_offenders`, `_has_strict_positive_combination` in Task 12) |
| `src/names.jl` | modify | `_edge_constants`' docstring (Task 8); `_constant_sides` for `:OnlyI` copies and classes (Task 10) |
| `src/rate_eq_derivation.jl` | modify | Class 5 inside `_num_den_exprs` (both methods, `_assert_support`) and `_compute_alpha`'s RE-cycle error removed (Task 8); `d_free` at free E (`_raw_symbolic_rate_polys`, Task 9); `_state_mask(am, :A)` on `_active_mask` (Task 10) |
| `src/mechanism_enumeration.jl` | modify | Class 3 in `_minimal_flips`' predicate on plain and allosteric mechanisms; exits drop (`expand_mechanisms`, `init_mechanisms` and `seed_mechanisms`); `_split_gain_test`'s pre-check and balancers removed (Task 8); `_redundant_copy_groups` per-state masks and `_degenerate` on the active graph, through an edge mask that `_has_vmax`, `_chemistry_equilibrates_both_sides` and `_merged_complexes` gain (Task 10); `_admits`, `_stops_inactive_turnover`, `_tagging_bases`, `_expand_to_allosteric` and change-tag through admission (Task 12) |
| `src/identify_rate_equation.jl` | modify | `_base_tier` drops flip children the guard rejects (Task 8) |
| `src/dsl.jl` | modify | Macro docstrings: the rejected classes and where each error surfaces (Task 8, R27); the `@allosteric_mechanism` docstring's `:OnlyI` rule and flip at free E (Task 14) |

Tests:

| File | Change | Responsibility |
|---|---|---|
| `test/test_admission.jl` | create (Task 2) | Verdict, message and oracle-agreement tests of the guard and admission on inline fixtures (RES §3.2); unit tests of networks, faces and classes 1-5; the tag sets `_testhelper_tag_sets`; the macros, the parser and the copy rule under the guard |
| `test/runtests.jl` | modify | Includes `test_admission.jl` after `test_oracles.jl` |
| `test/mechanism_definitions_for_test_enzyme_derivation.jl` | modify | Oracles of RES §3.1 (P4-24, P4-24a), `@_testhelper_unguarded` (sharing `_testhelper_mechanism_call` with `@_testhelper_with_raw`), the fault-per-face helpers, extended `_testhelper_consistent_point`, `_testhelper_law_value`, `_testhelper_single_subunit`, `_testhelper_graph`, `_testhelper_declared_dimension` and `_testhelper_detailed_balance`; flip form 1 (`_testhelper_steady_state_rate`'s `flip`); `_testhelper_check_law_against_limit`; `_testhelper_onlyi_fixtures`; `_testhelper_existence_agrees` and `_testhelper_plain_moves`; the two trap specs replaced by their variants (Task 7) |
| `test/test_oracles.jl` | modify | Oracle self-tests (RES §3.1); class 2 against the moves' RE-turnover test (Task 2); the iso families (Task 9); `:OnlyI` laws (Task 10); population tests P1-P7 (RES §3.3; P6 in Task 12); the regrouping test's pre-check rewritten (Task 8) |
| `test/test_types.jl` | modify | Bottomless, one-flip, S1, R1-R4, `:OnlyA` validator, ping-pong fused-release and `:OnlyI` testsets converted (RES §3.4); `_has_strict_positive_combination`'s test removed (Task 12) |
| `test/test_dsl.jl` | modify | Rejected fixtures converted (Task 8); catalytic `:OnlyI` messages and modifier `:OnlyI` accepted (Task 10) |
| `test/test_rate_eq_derivation.jl` | modify | "All-RE catalytic cycle raises" asserts on the macro (Task 8) |
| `test/allosteric_ground_truth.jl` | modify | The three ping-pong testsets become rejection tests with half-2 variants (P4-29b) |
| `test/test_mechanism_enumeration.jl` | modify | Exit drops and class 3 in the flip move (Task 8); split pre-check subjects (P4-29a); ter-ter count; `_testhelper_identifiable_rank` on unchecked laws; P4-1 re-pins and tagging fixtures (Task 12); `_testhelper_plain_moves` moved out (Task 13) |
| `test/test_identify_rate_equation.jl` | modify | `_base_tier` drops a rejected flip child (Task 8, updated in Task 9); `_degenerate` reads the active law on an `:OnlyI` fixture (Task 10) |
| `test/test_reaction_graph.jl` | modify | The census allosteric children's count the exit thins (Task 8); comments that name the flip form and `_testhelper_inactive_reach` seeded from free E (Task 9); the `:OnlyI` message pins (Task 10); the census population pins P4-1 moves (Task 12) |
| `test/reference/allosteric_golden_reference.txt` | modify | The LDH trap spec's variant replaces its two lines (Task 7) |
| `test/fixtures/canonical_text_golden.txt` | modify | Both trap specs' canonical text (Task 7) |

Docs and records:

| File | Change | Responsibility |
|---|---|---|
| `docs/src/deriving/mwc_allostery.md` | modify | Flip at E, limits, admission, rejected classes, `:OnlyI` rule, L's meaning for a free form sorting before E (P4-21) |
| `docs/src/identify/enumeration_engine.md` | modify | The tagging move (P4-1), change-tag, exit drops |
| `docs/src/developer.md` | modify | The guard's architecture; columns outside the law with `:OnlyI`; the split pre-check gone |
| `docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md` | modify | Every passage that asks for an exact face condition: section 7 "Admission" item 2 (lines 392-394), section 8 "Phase 4's first task" (478), the phase-4 row of section 15 (748) and section 17 row 2 (813) now say "sound; conservative on R-B1 and R-B2, listed" (Task 5, P4-9a) |
| `ROADMAP.md` | modify | Entries added, rewritten and retired in the task that causes each (CLAUDE.md): mixed corners (Task 2), the Freidlin–Wentzell Idea (Task 5) and its measured count (Task 11), the `BigInt` and mixed-group entries (Tasks 8, 12), the inactive-limit Planned entry retired (Task 9), the both-halves entry (Task 12), the rewrite entry (Task 14) |
| `docs/src/roadmap.md` | modify | The both-halves ping-pong limitation (Task 8, P4-3) |

Scratch, outside the repository:

| Path | Responsibility |
|---|---|
| `$S/phase45/probes/P4_T<n>_*.jl` | Cost probes (Tasks 2-6), Task 7's golden rewrite, Task 12's value check and split parents, Task 13's timing |
| `$S/phase45/scripts/P4_terter.jl`, `P4_bibi_kids.jl`, `P4_census_explain.jl` | Task 8's ter-ter timing, the value checks of its re-pins, and the census explanation |
| `$S/phase45/pumped/{pumped_check.jl, run_queue.sh, summarize.py, research_parent.jl}`, `data/` | Task 11's `:pumped` check (P4-9c) and Task 12's `move` stage |
| `$S/phase45/census/explain_tagging.py` | Task 12's census explanation |
| `$S/phase3/census/env_guard`, `census/P4_base`, `P4_T8`, `P4_T9`, `P4_T10`, `P4_T12` | Census environment and corpora (`census/README.md`; adapter `neutral_graph.jl`); Task 11 reads `census/base6` |
| `$S/phase45/gcscan/{mkenv,derive,identify}.jl`, `env112`, `env113` | Task 14's R18 remark scan |
| `$S/phase45/pr_body.md` | The draft PR's body (Task 14) |
| `$S/phase45/logs/` | Every run's log |

---

## Interface decisions that span tasks

Each task's Interfaces block gives exact signatures; these are the choices a later task relies
on, with the reason for each.

- One method on `_AnyMech` (`Union{Mechanism, AllostericMechanism}`) each: `_check_face`,
  `_face_limit_check`, `_faces_pass`, `_faces_verdict`, `_guard_message` (Tasks 2-3),
  `_check_well_posed` (Task 8) and, from Task 9, `_is_well_posed`. Tasks 4 and 6 insert one
  line each into `_guard_message` (the admission and class-5 branches); Tasks 4 and 5 replace
  `_guard(am)` and its docstring in full. Task 5's `_guard(am)` runs the face loop itself
  (`_check_face`, then `_admission_face`), so `_faces_verdict` and `_faces_pass` stay class 3
  only. `_support_message(who, m, v)` is one method on the catalytic graph.
- `BigFloat` keywords default to `BigFloat`s (`big(1.0)`, `big(10.0)^10`): Julia checks a
  keyword's declared type on its default too, and `big(1)` is a `BigInt`.
  `_testhelper_network` takes `sigma` and `theta` as one value or one per edge, which
  P4-24a's corners need (σ per blocked edge, θ per vanishing rapid-equilibrium binding).
  `_testhelper_steady_state_rate`'s `Λ = big(10)^30` is untyped and stays.
- `_testhelper_respeed` and `_testhelper_respeed_verdicts` judge one catalytic subunit: an
  oligomer is read at n = 1, as P4-24a and Task 0 judge it. No test makes an RE-speed verdict
  at n > 1, so P4-24a's wide set for one has no caller and the oracle leaves it out (YAGNI).
  `_testhelper_protomer_rate` reads `_testhelper_protomer_network` and takes `speeds` and `Λ`.
- `_testhelper_epath_verdict` returns `:diverges` too (a fast cycle of the ε-network turns
  over: RES §3.2's "RE-first diverges" for the B1 rows). Island tests accept either
  non-converging verdict, `:direction_dependent` or `:wrong_limit`: the probes read two
  orders, the ported oracle five.
- `_testhelper_law_value` and `_testhelper_single_subunit` take concrete mechanisms and read
  the law without the guard (Task 1), so the oracles judge laws the guard rejects (the
  compiled `_testhelper_single_subunit` lifts and calls the concrete one);
  `_testhelper_detailed_balance` (Task 10) and TME's `_testhelper_identifiable_rank`
  (Task 8) read the unchecked law the same way. Class 5 still throws inside `_num_den_exprs`
  (Task 8), so a fixture whose law is 0/0 on a face is judged with `law = :limit`.
- `_face_structure` returns `segments` (each segment's nodes, breadth-first), which the limit
  test and the messages read. Task 5 moves Task 3's Tarjan pass into
  `_strong_components(next, nodes)`, shared by `_bottom_sccs` and T0+;
  `_potential_on_sccs(net, keep, x)` takes `x`, the first substrate's `change` row
  (`_turnover_name`), since `_Network` does not record it. `_bridging_edges` returns an
  inclusion-minimal set of vanishing steps on a cheapest tree; the verdict's edge is its first.
- Messages name a form by `_form_text` (free E prints `E`) and a step by `_step_text`,
  phase 3's renderers in `src/names.jl`, whose step text `show` prints. Tagged chemistry
  is read in the direction that runs toward products (`_runs_toward_products`), so the
  island message names half 1 of a both-halves ping-pong in either writing. A face reads
  "at zero A and B, ".
- Task 8: `_minimal_flips(admissible, m::_AnyMech, eligible)` takes the mechanism itself, and
  `_expand_re_to_ss` passes `m`, so class 3 reads an allosteric parent's law network.
  `_split_gain_test` returns `gains` alone; the balancer partners leave
  `_expand_split_kinetic_group`. `_compute_alpha(m, flat, segments, names)` loses `conc_set`.
  `_assert_support(m, branches)` is class 5's throw, called by both `_num_den_exprs` methods
  with the polynomials they hold. Until Task 9, `_is_well_posed(am)` also runs the one-flip
  check (`_free_forms_flip_together`).
- Task 9 deletes `_testhelper_flip_form`; `_testhelper_steady_state_rate(...; flip = 1)` takes
  the flip form. `_testhelper_check_law_against_limit` builds the oracle's own network (own
  pruning), while phase 3's `_testhelper_check_law_against_oracle` reads the package's masks.
- Task 10 keeps `_check_tags(cat_allo_states, G, multiplicity)`, which the raw path calls
  before a graph exists; the `:OnlyI` rule is `_check_onlyi_groups`, run in the inner
  constructor. `_ordered_columns` returns `alaw`; `_active_law` and `_in_active_law` serve it
  and `_dependent_param_exprs`. `_testhelper_onlyi_fixtures()` is a function so that the
  fixtures file loads before the constructor accepts `:OnlyI`. `_has_vmax`,
  `_chemistry_equilibrates_both_sides` and `_merged_complexes` gain a last argument `on`,
  the edges that count (every edge by default, so their other callers and tests keep their
  calls), and `_degenerate`
  passes `_active_mask(m)`. Admission's turnover check (`_NO_ACTIVE_TURNOVER`) lands in
  Task 10, the first task whose constructor accepts a catalytic `:OnlyI` group.
- Task 12: `_expand_change_allo_state(am)` keeps one argument (the move reads nothing from the
  reaction); `_chemistry_groups(m)` loses its `group` and count arguments once
  `_onlya_offenders`, their only other caller, goes.
- Placement (CLAUDE.md): the tests of `expand_mechanisms` and the moves live in
  `test/test_mechanism_enumeration.jl`, so Task 8's "enumeration drops the children the guard
  rejects" (renamed by Task 12) and "class 3 inside the flip move emits the superset flips"
  go there; Task 10's copy-rule testset reads `_testhelper_onlyi_fixtures`, which that file's
  rule 1 forbids, so it goes to `test/test_admission.jl`. Task 2's class-2-against-the-moves
  testset reads `_testhelper_population_sample` and goes to `test/test_oracles.jl`. Task 13
  moves `_testhelper_plain_moves` to the fixtures file, where a focused run of
  `test_oracles.jl` sees it.
- Testset names are string literals that must fit 92 columns, so these are worded short:
  "the RE-speed oracle sees a trap's pattern and an unfed closed class", "the
  ε-path oracle rejects islands, and R2 only at P = 0", "the law and ε
  networks: a node per form and copy, two transitions per edge", "on a plain mechanism class
  2 equals the moves' RE-turnover test", "R2's limit at zero product depends on how fast its
  isomerization vanishes", and the population tests P4, P5, P6 and P7.

Smaller helpers Task 1 adds beside its oracles (Task 10 also reads `_testhelper_excludes`):
`_testhelper_mechanism_call`, `_testhelper_excludes`, `_testhelper_flip_ratio`, `_testhelper_gth`, `_testhelper_offsets`,
`_testhelper_relgap`, `_testhelper_package_law`, `_testhelper_limit_law`,
`_testhelper_excluding`, `_testhelper_component_potentials`,
`_testhelper_standard_directions`, `_testhelper_corner_directions`,
`_testhelper_epsilon_span`, `_testhelper_path_limit`, `_testhelper_feasible`.

---

## Tasks

### Task 0: Confirmatory research — DONE

**Files:** none in the repository. Reports: `$S/phase45/research/task0/{a,b,c,defg}.md`; RES §5.

**Interfaces:**
- Consumes: the research code under `$S/phase45/research/` on neutral step lists.
- Produces: no code. It hands later tasks the changes P4-31 lists: Task 1 ports the oracle
  with P4-24a's settings (every corner when 2^k ≤ 128, else 128 drawn; κ ∈ {RE first, 8, 2,
  ½, ε first}; oligomers judged at n = 1; the wide set at n > 1 has no caller); Tasks 5 and
  10 test the joint-R1 check (P4-8b); Task 11 checks the 5,342 `:pumped` drops (P4-9c);
  RES §3.4's three ill-posed variants take RES §5.2 item 1's replacements (P4-29b).

**Gate (passed):** RES §5.1. (a) 0 misses under corners and five orders, false alarms only in
R-B1 (332) and R-B2 (21); (b) 88,434 children, all admitted, none duplicated, 0 oracle misses
on 54,760 face rows, each of the 30 regulated ping-pong parents keeps 2-4 V-types; (c) the n = 1
verdict holds at n = 2; (d) 0 misses on uni-bi; (e) 357 admitted mixed-corner laws (P4-12a);
(f) both trap variants pass, three §3.4 variants replaced; (g) 130 of 130 ter-ter islands
direction-dependent. No miss and no new class, so P4-8 stays closed.

- [x] **Step 1: Run the confirmatory research.** Done 2026-10-14, about 4.6 h of Julia;
  results in RES §5, ruled by P4-31. Part (b)'s sampled ε-path stands.

---

### Task 1: Test oracles and their self-tests

The oracles of RES §3.1, ported from `research/oracle/` and `vi/vor.jl` onto the graph read
API, with P4-24a's settings. They never call the guard's code (P4-24). Where a research
script differs from RES §1 or the rulings, this task follows the rulings: both potential kinds
rise by integers (P4-25; `oracle/oracles.jl` jitters the second), κ is P4-24a's set
(`vi/vor.jl` runs seven orders), and corners are 128 when drawn (Task 0(a)'s first sweep
drew 32, RES §5.2 item 4).

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (new section "Phase-4
  oracles" after `_testhelper_population`; `@_testhelper_with_raw`, which shares the new
  `_testhelper_mechanism_call`; `_testhelper_consistent_point`,
  `_testhelper_steady_state_rate`, `_testhelper_law_value`, `_testhelper_single_subunit`)
- Modify: `test/test_oracles.jl` (ABOUTME lines; eleven testsets appended at the end)

**Interfaces:**
- Consumes (phase 3): `EnzymeRates._raw_mechanism(block, allosteric::Bool)`,
  `EnzymeRates.Mechanism(raw)`, `EnzymeRates.AllostericMechanism(raw)` and the composition
  constructor `AllostericMechanism(catalytic, cat_allo_states, multiplicity, sites)`;
  `catalytic`, `m.edges`, `m.forms`, `m.table`, `m.group`, `cat_allo_state`,
  `_is_plain_binding`, `entries`, `regulatory_sites`, `metabolites`, `_reg_K`,
  `_dependent_param_exprs(m::_AnyMech)`, `_flat_edges`; `_testhelper_graph(m)`,
  `_testhelper_cycle_basis(g)`, `_testhelper_lift`, `_testhelper_steady_state_rate`,
  `_testhelper_consistent_point`.
- Produces (fixtures file):
  ```julia
  _testhelper_mechanism_call(call, who::String)
      # -> (allosteric::Bool, raw::Union{RawMechanism, RawAllostericMechanism})
  macro _testhelper_unguarded(call)   # -> EnzymeRates.Mechanism | EnzymeRates.AllostericMechanism
  _testhelper_faces(m)::Vector{Set{Symbol}}          # RES §1.2 order; faces[1] == Set{Symbol}()
  struct _testhelper_Link
      i::Int; j::Int; f::BigFloat; r::BigFloat; w::Int; fast::Bool
  end
  struct _testhelper_Network
      nodes::Vector{Tuple{Int, Int}}
      links::Vector{_testhelper_Link}
  end
  _testhelper_exclusion(am, c::Int)::NamedTuple{(:component, :arcs, :failure, :lost)}
  _testhelper_draw_speeds(m, rng)::Matrix{BigFloat}                     # edges × 2
  _testhelper_network(m, point, concs::AbstractDict; speeds::Matrix{BigFloat},
                      mode::Symbol = :law, ε::BigFloat = big(1.0),
                      potentials::Vector{Vector{BigFloat}} = Vector{BigFloat}[],
                      sigma::Union{BigFloat, Vector{BigFloat}} = big(1.0),
                      theta::Union{BigFloat, Vector{BigFloat}} = big(1.0) / 2
                      )::_testhelper_Network                          # mode ∈ (:law, :vanish)
  _testhelper_closed_classes(n::Int, arcs::Vector{Tuple{Int, Int}})::Vector{Vector{Int}}
  _testhelper_steady_rate(nw::_testhelper_Network; Λ::BigFloat)::NamedTuple{(:v, :nclasses)}
  _testhelper_re_limit(nw::_testhelper_Network)::BigFloat   # NaN non-unique, Inf divergent
  _testhelper_respeed(m, point, concs::AbstractDict, face::Set{Symbol}; rng,
                      law::Symbol = :package, Λs = (big(10.0)^10, big(10.0)^20),
                      npat::Int = 2)::Symbol
      # ∈ (:agree, :pattern, :nonunique, :diverges, :unconverged, :law_nan, :mismatch)
  _testhelper_respeed_verdicts(m; rng, npoints::Int = 1,
                               law::Symbol = :package)::Vector{Symbol}   # one per face
  _testhelper_corners(k::Int, rng)::Vector{Vector{Int}}
  _testhelper_epath_verdict(m; rng, face::Set{Symbol} = Set{Symbol}(),
                            corners::Bool = false, precision::Int = 256)::Symbol
      # ∈ (:no_path, :nonunique, :diverges, :del_nonunique, :converges,
      #    :direction_dependent, :wrong_limit)
  _testhelper_in_epath_window(m; rng)::Bool
  _testhelper_certified(am, c::Int)::Bool
  _testhelper_lp_existence(am, c::Int)::Bool
  _testhelper_uniuni_grid()::Vector{EnzymeRates.AllostericMechanism}       # 108
  _testhelper_protomer_network(m, point, concs::AbstractDict; n::Int,
                               speeds::Matrix{BigFloat})::_testhelper_Network
  _testhelper_protomer_rate(m, point, concs::AbstractDict; n::Int,
                            speeds::Matrix{BigFloat} = fill(big(1.0), size(point.K)),
                            Λ::BigFloat = big(10.0)^30
                            )::BigFloat
  _testhelper_consistent_point(m, rng)   # adds K::Matrix{BigFloat}, kf::Matrix{BigFloat},
                                         # edge × conformation, every edge
  _testhelper_law_value(m, params, concs)  # also takes a concrete mechanism, unchecked
  ```

**Gate:** the eleven testsets below pass; RES §3.1's self-tests pass with P4-24a's settings,
except the two RES §3.1 moves to later tasks (the 512-bit `bridge_REcycle2` run needs an
`:OnlyI` group, Task 10; the iso families at n = 1-3 need the flip at E, Task 9); the new
testsets add at most 4 s to `test/test_oracles.jl`; the full suite is green and pristine.

- [ ] **Step 1: Record the branch point.** Every shell call starts fresh, so each command
  sets `S` before it reads it. If the worktree does not exist yet, cut it from
  graph-mechanism (after the phase-3 PR) and instantiate it first:
  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; mkdir -p $S/phase45/logs $S/phase45/probes && git -C /home/denis.linux/.julia/dev/EnzymeRates worktree add $S/phase45/guard_wt -b graph-guard graph-mechanism
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T1_inst.log 1800 --project=. -e 'using Pkg; Pkg.instantiate()'
  ```
  Then:

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; mkdir -p $S/phase45/logs $S/phase45/probes && cd $S/phase45/guard_wt && git status --short && git log --oneline -1
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T1_base.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'
  ```
  Run the suite in the background and poll. Expected: a clean tree on graph-guard; the log
  ends `exit 0 elapsed <n>s` after `Testing EnzymeRates tests passed`, no warning after
  "Testing Running tests". Record the pass count N0, the wall time and `maxrss_kb` M0 from
  the log's last lines: every later full run compares with them.

- [ ] **Step 2: Write the failing tests for the unguarded macro and the oracle's faces.**

  In `test/test_oracles.jl`, replace the two ABOUTME lines with:
  ```julia
  # ABOUTME: Tests of the test helpers' mechanism reader and of the exact, RE-speed, ε-path,
  # ABOUTME: LP and n-protomer oracles, and of the laws over regenerated populations.
  ```
  Append at the end of the file:
  ```julia
  @testset "the unguarded macro builds what the guard refuses" begin
      # R2 (spec section 7), which the macro refuses: at zero product its inactive limit
      # depends on how fast the isomerization vanishes. The unguarded macro builds it, as
      # the constructors do from the walker's raw value.
      r2 = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)      :: OnlyA
              E(S) <--> E(P)    :: EqualAI
              F + P ⇌ E(P)      :: EqualAI
              F <--> E          :: OnlyA
          end
      end
      @test r2 isa EnzymeRates.AllostericMechanism
      @test r2 == EnzymeRates.AllostericMechanism(EnzymeRates._raw_mechanism(:(begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)      :: OnlyA
              E(S) <--> E(P)    :: EqualAI
              F + P ⇌ E(P)      :: EqualAI
              F <--> E          :: OnlyA
          end
      end), true))
      @test_throws ErrorException @allosteric_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)      :: OnlyA
              E(S) <--> E(P)    :: EqualAI
              F + P ⇌ E(P)      :: EqualAI
              F <--> E          :: OnlyA
          end
      end
      mm = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E + P
          end
      end
      @test mm == EnzymeRates.Mechanism(@enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E + P
          end
      end)
  end

  @testset "the oracle's faces follow the spec's list" begin
      bibi = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) + B ⇌ E(A, B)
              E(A, B) <--> E(P, Q)
              E(P, Q) ⇌ E(Q) + P
              E(Q) ⇌ E + Q
          end
      end
      @test _testhelper_faces(bibi) == [Set{Symbol}(), Set([:P, :Q]), Set([:A, :B]),
                                        Set([:A]), Set([:B]), Set([:P]), Set([:Q])]
      inhibited = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          modifiers: I
          steps: begin
              E + A ⇌ E(A)
              E(A) + B ⇌ E(A, B)
              E(A, B) <--> E(P, Q)
              E(P, Q) ⇌ E(Q) + P
              E(Q) ⇌ E + Q
              (E + I ⇌ E(I), E(A) + I ⇌ E(A, I))
          end
      end
      faces = _testhelper_faces(inhibited)
      @test length(faces) == 14
      @test faces[8:14] == [union(z, Set([:I])) for z in faces[1:7]]
      terter = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B, C
          products: P, Q, R
          steps: begin
              E + A ⇌ E(A)
              E(A) + B ⇌ E(A, B)
              E(A, B) + C ⇌ E(A, B, C)
              E(A, B, C) <--> E(P, Q, R)
              E(P, Q, R) ⇌ E(Q, R) + P
              E(Q, R) ⇌ E(R) + Q
              E(R) ⇌ E + R
          end
      end
      @test length(_testhelper_faces(terter)) == 9
      # P's modifier copy carries P's name: it is zeroed with P, and P adds no modifier face.
      copy_P = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          modifiers: P
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E(P)
              E(P) ⇌ E + P
              E(S) + P::Mod ⇌ E(P::Mod, S)
          end
      end
      @test _testhelper_faces(copy_P) == [Set{Symbol}(), Set([:P]), Set([:S])]
      # A regulatory-site ligand is never zeroed.
      regulated = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          allosteric_regulators: R
          steps: begin
              E + S ⇌ E(S)      :: NonequalAI
              E(S) <--> E(P)    :: EqualAI
              E(P) ⇌ E + P      :: EqualAI
          end
          regulatory_site(1): R::OnlyA
      end
      @test _testhelper_faces(regulated) == [Set{Symbol}(), Set([:P]), Set([:S])]
  end
  ```

- [ ] **Step 3: Run them and watch them fail.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T1_red1.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_oracles.jl")'
  ```
  Expected: the file's earlier testsets pass, then the include fails with `LoadError:
  UndefVarError: `@_testhelper_unguarded` not defined`, raised when the first new testset is
  expanded.

- [ ] **Step 4: Add the unguarded macro and the oracle's faces.**

  In `test/mechanism_definitions_for_test_enzyme_derivation.jl`, the new macro and
  `@_testhelper_with_raw` check their call the same way, so they share one helper. Directly
  above `@_testhelper_with_raw`'s docstring insert:
  ```julia
  """
      _testhelper_mechanism_call(call, who) -> (allosteric, raw)

  Whether `call`, an `@enzyme_mechanism` or `@allosteric_mechanism` call, is the allosteric
  one, and the raw value the package's DSL walker reads from its block; an error naming the
  test macro `who` for any other expression.
  """
  function _testhelper_mechanism_call(call, who::String)
      head = call isa Expr && call.head === :macrocall ? call.args[1] : nothing
      head in (Symbol("@enzyme_mechanism"), Symbol("@allosteric_mechanism")) ||
          error("$who wraps an @enzyme_mechanism or @allosteric_mechanism call; got $call")
      allosteric = head === Symbol("@allosteric_mechanism")
      allosteric, EnzymeRates._raw_mechanism(call.args[end], allosteric)
  end
  ```
  and replace the body of `macro _testhelper_with_raw(call)` (its docstring stays) with:
  ```julia
  macro _testhelper_with_raw(call)
      _, raw = _testhelper_mechanism_call(call, "@_testhelper_with_raw")
      esc(:(($call, $raw)))
  end
  ```
  The census's `corpus.jl` evaluates the file's statements in order, so it sees the helper
  before the specs expand the macro.

  Then, after the function `_testhelper_population` and before the line
  `# ── Mechanism test specifications ──`, insert:
  ```julia
  # ── Phase-4 oracles ─────────────────────────────────────────────────────────
  #
  # The well-posedness guard's oracles (spec section 14): mass-action networks of one
  # catalytic subunit, read from the graph read API and never from the guard's code, whose
  # conformations flip at free E (form 1); the RE-speed oracle on the guard's faces; the
  # ε-path oracle of admission, along integer potentials or the corners of the admissible
  # cone and five orders of the two limits; the exact LP of admission's existence test; and
  # the concerted n-protomer oracle.

  """
      @_testhelper_unguarded call

  The concrete mechanism that `call`, an `@enzyme_mechanism` or `@allosteric_mechanism`
  call, writes: the constructors on the DSL walker's raw value, without the well-posedness
  check the macros run, so a test can hold a mechanism the macros refuse.
  """
  macro _testhelper_unguarded(call)
      allosteric, raw = _testhelper_mechanism_call(call, "@_testhelper_unguarded")
      esc(allosteric ? :(EnzymeRates.AllostericMechanism($raw)) :
                       :(EnzymeRates.Mechanism($raw)))
  end

  """
  The faces of `m` (concrete or compiled) that the guard checks (spec section 8), each as
  the set of names it zeroes, in the guard's order: the interior, every product zero, every
  substrate zero, each substrate and then each product alone on a side of two or more; then
  each of these with every modifier-only name zero too. A reactant's `::Mod` copy carries
  its reactant's name, so it is zeroed with it; a regulatory-site ligand never is.
  """
  function _testhelper_faces(m)
      cm = EnzymeRates.catalytic(_testhelper_lift(m))
      mods = Set(x.name for x in cm.table
                 if x.modifier && !(x.name in cm.substrates || x.name in cm.products))
      base = [Set{Symbol}(), Set(cm.products), Set(cm.substrates)]
      length(cm.substrates) >= 2 && append!(base, [Set([x]) for x in cm.substrates])
      length(cm.products) >= 2 && append!(base, [Set([x]) for x in cm.products])
      isempty(mods) ? base : [base; [union(z, mods) for z in base]]
  end
  ```

- [ ] **Step 5: Run them and watch them pass.** The Step 3 command, log `P4_T1_green1.log`.
  Expected: `exit 0`; "the unguarded macro builds what the guard refuses" passes 4 tests and
  "the oracle's faces follow the spec's list" 6, every earlier testset of the file still
  passes, and the output holds no warning.

- [ ] **Step 6: Write the failing tests of the networks, the solvers and the RE-speed oracle.**

  Append to `test/test_oracles.jl`:
  ```julia
  @testset "the RE-speed oracle's error falls as 1/Λ on Michaelis–Menten" begin
      mm = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E + P
          end
      end
      setprecision(BigFloat, 256) do
          rng = Random.MersenneTwister(11)
          point = _testhelper_consistent_point(mm, rng)
          concs = Dict(:S => big"0.7", :P => big"1.3")
          speeds = _testhelper_draw_speeds(mm, rng)
          net = _testhelper_network(mm, point, concs; speeds)
          law = _testhelper_package_law(mm, point, concs)
          # The steady state approaches the law as c/Λ, so the gap times Λ is constant.
          gap(Λ) = (_testhelper_steady_rate(net; Λ).v - law) * Λ
          @test isapprox(gap(big(10.0)^10), gap(big(10.0)^20); rtol = 1e-3)
          @test abs(_testhelper_re_limit(net) - law) <= big(10)^-60 * abs(law)
          @test _testhelper_respeed_verdicts(mm; rng) == [:agree, :agree, :agree]
      end
  end

  @testset "the RE-speed oracle sees a trap's pattern and an unfed closed class" begin
      # LDH-like: at zero products the forms behind E(P) + Q ⇌ E(P, Q) trap enzyme, and
      # at Q = 0 the segment of free E splits; both limits move with the RE speeds.
      ldh = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) + B <--> E(A, B)
              E(A, B) <--> E(P, Q)
              E(P) + Q ⇌ E(P, Q)
              E(Q) + P ⇌ E(P, Q)
              E + P <--> E(P)
              E + Q ⇌ E(Q)
          end
      end
      @test _testhelper_respeed_verdicts(ldh; rng = Random.MersenneTwister(1)) ==
            [:agree, :pattern, :agree, :agree, :agree, :agree, :pattern]
      # At zero products F(P, P) is a closed set the enzyme never leaves.
      fuzz5 = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E ⇌ F
              E + I ⇌ F(I)
              (E + P <--> E(P), E(P) + P <--> E(P, P))
              E + S <--> E(P)
              E + S ⇌ F(S)
              E(P, P) + S <--> F(P, P) + P
          end
      end
      @test _testhelper_respeed_verdicts(fuzz5; rng = Random.MersenneTwister(2))[2] ===
            :nonunique
  end

  @testset "the RE-speed oracle sees RE turnover beside a steady-state route diverge" begin
      beside = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) ⇌ E(P)
              E + P ⇌ E(P)
              E(S) <--> F(P)
              F(P) <--> E + P
          end
      end
      @test _testhelper_respeed_verdicts(beside; rng = Random.MersenneTwister(3))[1] ===
            :diverges
  end
  ```

- [ ] **Step 7: Run them and watch them fail.** The Step 3 command, log `P4_T1_red2.log`.
  Expected: the three new testsets each report `Error During Test` with `UndefVarError:
  `_testhelper_draw_speeds` not defined` (Michaelis–Menten, whose first new call it is) and
  `UndefVarError: `_testhelper_respeed_verdicts` not defined` (the other two); every other
  testset passes.

- [ ] **Step 8: Add the networks, the solvers and the RE-speed oracle.**

  (a) In `_testhelper_steady_state_rate`, replace the block that starts
  `if g.catalytic_multiplicity !== nothing` (its comment, the local
  `_testhelper_site_factor` and the flip `push!`) with the version below; the comment's
  explanation moves into `_testhelper_flip_ratio`'s docstring:
  ```julia
          if g.catalytic_multiplicity !== nothing
              flip = _testhelper_flip_form(g)
              push!(edges, (id(:A, flip), id(:I, flip),
                            _testhelper_flip_ratio(mech, point, concs), big(1), 0))
          end
  ```

  (b) In `_testhelper_consistent_point`, replace the returned `(; vals, Keq, active,
  inactive, L)` (its last line) with:
  ```julia
      # Every edge's equilibrium constant in both conformations and, for a steady-state
      # edge, its forward rate: the law's own where the law holds the edge, an `:EqualAI`
      # edge's active one, and otherwise one draw per kinetic group and conformation (each
      # member of an `:OnlyA` group in the inactive conformation its own). The phase-4
      # oracles read these for the edges the inactive law deletes or prunes.
      # (`forward` holds the forward rates: a local named `kf` would be the one the closure
      # `_testhelper_conformation_steps` assigns.)
      cm = EnzymeRates.catalytic(mech)
      flat = EnzymeRates._flat_edges(cm)
      K = ones(BigFloat, length(cm.edges), 2)
      forward = ones(BigFloat, length(cm.edges), 2)
      rates = Dict{Tuple{Int, Int}, BigFloat}()
      for c in 1:nconf, s in g.steps
          names = c == 1 || tag(s.group) === :EqualAI ? s.names : s.inactive_names
          s.rapid_equilibrium || isempty(names) || (rates[(s.group, c)] = vals[names[1]])
      end
      for c in 1:nconf, (j, s) in enumerate(g.steps)
          k = flat[j]
          K[k, c] = exp(u[(c - 1) * nf + s.from] + sum(gm, s.consumed; init = big(0)) -
                        u[(c - 1) * nf + s.to] - sum(gm, s.released; init = big(0)))
          s.rapid_equilibrium && continue
          key = c == 2 && tag(s.group) === :OnlyA ? (-k, 2) : (s.group, c)
          forward[k, c] = get!(() -> exp(4 * big(rand(rng)) - 2), rates, key)
      end
      allo || (K[:, 2] .= K[:, 1]; forward[:, 2] .= forward[:, 1])
      (; vals, Keq, active, inactive, L, K, kf = forward)
  ```
  In its docstring, change the signature line to
  `    _testhelper_consistent_point(m, rng) -> (; vals, Keq, active, inactive, L, K, kf)` and
  append, after the paragraph on `active` and `inactive`:
  ```julia
  `K[k, c]` is edge `k`'s equilibrium constant in its stored direction in conformation `c`
  and `kf[k, c]`, for a steady-state edge, its forward rate, for every edge of the
  catalytic graph, the ones the inactive law deletes or prunes included; a plain
  mechanism repeats column 1 in column 2. The law's constants are among them.
  ```

  (c) Change the first lines of `_testhelper_law_value`'s body and its docstring:
  ```julia
  """
  The rate law of `m` per enzyme, evaluated in 256-bit BigFloat from the polynomials
  `rate_equation` compiles (`_num_den_exprs`) at the parameter values `params` (by name: the
  fitted constants, `Keq` and, for an allosteric mechanism, `L` and the regulator
  constants) and the concentrations `concs`. A compiled `m` goes through its derivation
  entry, which runs the boundary check; a concrete `m` is read as it is, so the oracles can
  evaluate the law of a mechanism the guard rejects. The dependent constants are computed
  from `_dependent_param_exprs`, as the compiled body computes them.
  """
  function _testhelper_law_value(m, params::AbstractDict, concs::AbstractDict)
      setprecision(BigFloat, 256) do
          mech = m isa EnzymeRates.AbstractEnzymeMechanism ?
                 EnzymeRates._derivation_input(typeof(m)) : m
  ```
  (the rest of the body is unchanged).

  (d) Replace the two methods of `_testhelper_single_subunit` with the four below, and
  widen its docstring's first line to "`m`, compiled or concrete, with one catalytic
  subunit". The compiled allosteric method lifts the mechanism and calls the concrete one;
  at the branch point it lifts the mechanism and calls the composition constructor itself.
  ```julia
  _testhelper_single_subunit(m::EnzymeRates.EnzymeMechanism) = m
  _testhelper_single_subunit(m::EnzymeRates.AllostericEnzymeMechanism) =
      EnzymeRates.AllostericEnzymeMechanism(
          _testhelper_single_subunit(EnzymeRates.AllostericMechanism(m)))
  _testhelper_single_subunit(m::EnzymeRates.Mechanism) = m
  _testhelper_single_subunit(am::EnzymeRates.AllostericMechanism) =
      EnzymeRates.AllostericMechanism(am.catalytic, copy(am.cat_allo_states), 1,
                                      am.regulatory_sites)
  ```

  (e) At the end of the "Phase-4 oracles" section, add:
  ```julia
  """Whether conformation `c` (1 active, 2 inactive) of `m` excludes edge `k` of its
  catalytic graph: an `:OnlyA` group's edge in the inactive conformation, an `:OnlyI`
  group's in the active one. A plain mechanism excludes nothing."""
  function _testhelper_excludes(m, k::Int, c::Int)
      m isa EnzymeRates.AllostericMechanism || return false
      tag = EnzymeRates.cat_allo_state(m, m.catalytic.group[k])
      (c == 2 && tag === :OnlyA) || (c == 1 && tag === :OnlyI)
  end

  """One reversible step of an oracle network between nodes `i` and `j`: forward rate `f`
  and reverse rate `r`, concentrations included, `w` copies of the first substrate taken up
  from `i` to `j`; a `fast` link runs at Λ times its rates."""
  struct _testhelper_Link
      i::Int
      j::Int
      f::BigFloat
      r::BigFloat
      w::Int
      fast::Bool
  end

  """A mass-action network of the oracles: its nodes, each `(form, conformation)` (`(state,
  conformation)` for the n-protomer oracle's product states), and its links."""
  struct _testhelper_Network
      nodes::Vector{Tuple{Int, Int}}
      links::Vector{_testhelper_Link}
  end

  """The equilibrium ratio of the flip at free E, inactive over active, of the allosteric
  mechanism `mech` at `concs`: `L` times the ratio of the two conformations'
  regulatory-site factors. Regulators bind each conformation at rapid equilibrium whatever
  its catalytic form, so summed over regulator occupancy the flip takes this ratio."""
  function _testhelper_flip_ratio(mech, point, concs::AbstractDict)
      function _testhelper_site_factor(inactive)
          prod((1 + sum((concs[lig] / point.vals[K]
                         for (lig, tag) in zip(site.ligands, site.allo_states)
                         for K in (EnzymeRates._reg_K(mech, lig, tag, inactive),)
                         if K !== nothing); init = big(0)))^site.multiplicity
               for site in EnzymeRates.regulatory_sites(mech); init = big(1))
      end
      point.L * _testhelper_site_factor(true) / _testhelper_site_factor(false)
  end

  """Speed factors of the rapid-equilibrium edges, one per edge of `m`'s catalytic graph and
  conformation, log-uniform over `10^±3`: a speed pattern of the RE-speed oracle."""
  _testhelper_draw_speeds(m, rng) =
      BigFloat[big(10)^(3 * (2 * big(rand(rng)) - 1))
               for _ in EnzymeRates.catalytic(_testhelper_lift(m)).edges, _ in 1:2]

  """
      _testhelper_network(m, point, concs; speeds, mode = :law, ε = 1, potentials, sigma,
                          theta) -> _testhelper_Network

  The mass-action network of `m`'s catalytic graph, one catalytic subunit, at the constants
  of `point` (`_testhelper_consistent_point`) and the concentrations `concs`: one copy of
  the forms for a plain mechanism, two for an allosteric one, joined by a fast flip between
  the copies of free E at `_testhelper_flip_ratio`. A rapid-equilibrium edge of
  conformation `c` runs fast at speed `speeds[k, c]`, a steady-state one at
  `point.kf[k, c]`.

  `mode = :law` is the law's network: each conformation keeps the edges it does not
  exclude, on the forms they join to free E (deletion and pruning). `mode = :vanish` is the
  ε-network: every form and edge, each excluded edge vanishing as ε → 0 along a direction,
  with `potentials[c][f]` the potential of form `f`'s component in conformation `c` and
  `d = potentials[c][to] - potentials[c][from]`. A plain binding's equilibrium constant
  scales by `ε^d`: a rapid-equilibrium one stays fast, its on rate scaled by `ε^(θ d)` and
  its off rate by `ε^(-(1 - θ) d)`; a steady-state one keeps its off rate. Any other
  excluded edge is slow, its rates scaled by `ε^(σ + max(±d, 0))`. `sigma` (σ) and `theta`
  (θ) are one value for every edge or one per edge.
  """
  function _testhelper_network(m, point, concs::AbstractDict; speeds::Matrix{BigFloat},
                               mode::Symbol = :law, ε::BigFloat = big(1.0),
                               potentials::Vector{Vector{BigFloat}} = Vector{BigFloat}[],
                               sigma::Union{BigFloat, Vector{BigFloat}} = big(1.0),
                               theta::Union{BigFloat, Vector{BigFloat}} = big(1.0) / 2)
      mech = _testhelper_lift(m)
      cm = EnzymeRates.catalytic(mech)
      nf = length(cm.forms)
      nc = mech isa EnzymeRates.AllostericMechanism ? 2 : 1
      x1 = first(cm.substrates)
      met_names(c) = [cm.table[i].name for i in EnzymeRates.entries(c)]
      keep = [trues(nf) for _ in 1:nc]
      if mode === :law
          for c in 1:nc
              reached = falses(nf)
              reached[1] = true
              grew = true
              while grew
                  grew = false
                  for (k, e) in enumerate(cm.edges)
                      !_testhelper_excludes(mech, k, c) &&
                          reached[e.from] != reached[e.to] || continue
                      reached[e.from] = reached[e.to] = true
                      grew = true
                  end
              end
              keep[c] = reached
          end
      end
      nodes = [(f, c) for c in 1:nc for f in 1:nf if keep[c][f]]
      id = Dict(node => n for (n, node) in enumerate(nodes))
      links = _testhelper_Link[]
      for c in 1:nc, (k, e) in enumerate(cm.edges)
          keep[c][e.from] && keep[c][e.to] || continue
          excluded = _testhelper_excludes(mech, k, c)
          mode === :law && excluded && continue
          K, kf = point.K[k, c], point.kf[k, c]
          fast = e.rapid_equilibrium
          f, r = fast ? (speeds[k, c] * sqrt(K), speeds[k, c] / sqrt(K)) : (kf, kf / K)
          if excluded
              d = potentials[c][e.to] - potentials[c][e.from]
              if EnzymeRates._is_plain_binding(cm, e)
                  θ = theta isa BigFloat ? theta : theta[k]
                  f, r = fast ? (f * ε^(θ * d), r * ε^(-(1 - θ) * d)) : (f * ε^d, r)
              else
                  σ = sigma isa BigFloat ? sigma : sigma[k]
                  f, r, fast = f * ε^(σ + max(d, 0)), r * ε^(σ + max(-d, 0)), false
              end
          end
          w = count(==(x1), met_names(e.consumed)) - count(==(x1), met_names(e.released))
          conc_in = prod((concs[x] for x in met_names(e.consumed)); init = big(1))
          conc_out = prod((concs[x] for x in met_names(e.released)); init = big(1))
          push!(links, _testhelper_Link(id[(e.from, c)], id[(e.to, c)], f * conc_in,
                                        r * conc_out, w, fast))
      end
      if nc == 2
          ratio = _testhelper_flip_ratio(mech, point, concs)
          push!(links, _testhelper_Link(id[(1, 1)], id[(1, 2)], sqrt(ratio),
                                        1 / sqrt(ratio), 0, true))
      end
      _testhelper_Network(nodes, links)
  end

  """The closed classes (strongly connected, no transition leaving) of the digraph on nodes
  `1:n` with transitions `arcs`, by Tarjan's algorithm, each as its sorted nodes."""
  function _testhelper_closed_classes(n::Int, arcs::Vector{Tuple{Int, Int}})
      next = [Int[] for _ in 1:n]
      foreach(((a, b),) -> push!(next[a], b), arcs)
      index, low, comp = zeros(Int, n), zeros(Int, n), zeros(Int, n)
      onstack, stack = falses(n), Int[]
      counter, ncomp = 0, 0
      for s in 1:n
          index[s] == 0 || continue
          counter += 1
          index[s] = low[s] = counter
          push!(stack, s)
          onstack[s] = true
          work = [(s, 1)]
          while !isempty(work)
              u, i = work[end]
              if i <= length(next[u])
                  work[end] = (u, i + 1)
                  v = next[u][i]
                  if index[v] == 0
                      counter += 1
                      index[v] = low[v] = counter
                      push!(stack, v)
                      onstack[v] = true
                      push!(work, (v, 1))
                  elseif onstack[v]
                      low[u] = min(low[u], index[v])
                  end
              else
                  pop!(work)
                  isempty(work) || (low[work[end][1]] = min(low[work[end][1]], low[u]))
                  if low[u] == index[u]
                      ncomp += 1
                      while true
                          v = pop!(stack)
                          onstack[v] = false
                          comp[v] = ncomp
                          v == u && break
                      end
                  end
              end
          end
      end
      leaves = falses(ncomp)
      for (a, b) in arcs
          comp[a] != comp[b] && (leaves[comp[a]] = true)
      end
      [findall(==(k), comp) for k in 1:ncomp if !leaves[k]]
  end

  """The stationary distribution of the irreducible chain with off-diagonal rates `R[i, j]`
  (`i → j`), by Grassmann–Taksar–Heyman state reduction: sums of non-negative terms only, so
  every entry keeps its relative precision whatever the spread of the rates."""
  function _testhelper_gth(R::Matrix{BigFloat})
      n = size(R, 1)
      A = copy(R)
      for k in n:-1:2
          s = sum(A[k, j] for j in 1:(k - 1))
          iszero(s) && error("_testhelper_gth: state $k has no exit to states 1:$(k - 1)")
          for i in 1:(k - 1)
              iszero(A[i, k]) && continue
              A[i, k] /= s
              for j in 1:(k - 1)
                  i == j || iszero(A[k, j]) || (A[i, j] += A[i, k] * A[k, j])
              end
          end
      end
      p = zeros(BigFloat, n)
      p[1] = 1
      for j in 2:n
          p[j] = sum((p[i] * A[i, j] for i in 1:(j - 1) if !iszero(A[i, j])); init = big(0))
      end
      p ./ sum(p)
  end

  """Each node's potential of the first substrate along `links` (node `j`'s is node `i`'s
  plus the link's uptake `w`), component by component on nodes `1:n`; `nothing` when a
  cycle of the links takes the substrate up."""
  function _testhelper_offsets(n::Int, links)
      h = fill(typemin(Int), n)
      adjacent = [Tuple{Int, Int}[] for _ in 1:n]
      for l in links
          push!(adjacent[l.i], (l.j, l.w))
          push!(adjacent[l.j], (l.i, -l.w))
      end
      for s in 1:n
          h[s] == typemin(Int) || continue
          h[s] = 0
          todo = [s]
          while !isempty(todo)
              u = pop!(todo)
              for (v, w) in adjacent[u]
                  if h[v] == typemin(Int)
                      h[v] = h[u] + w
                      push!(todo, v)
                  elseif h[v] != h[u] + w
                      return nothing
                  end
              end
          end
      end
      h
  end

  """
      _testhelper_steady_rate(nw; Λ) -> (; v, nclasses)

  The steady-state turnover per enzyme of `nw` (net uptake of the first substrate) with
  every fast link at `Λ` times its rates, in the current BigFloat precision; `v` is NaN when
  the chain has more than one closed class (`nclasses` of them). It weights each link's flux
  by its uptake plus the fall of the first substrate's potential along the fast links that
  run both ways, so those links, whose fluxes nearly cancel at large Λ, drop out; when they
  close a cycle that takes the substrate up, every link counts with its own uptake.
  """
  function _testhelper_steady_rate(nw::_testhelper_Network; Λ::BigFloat)
      n = length(nw.nodes)
      rates(l) = l.fast ? (Λ * l.f, Λ * l.r) : (l.f, l.r)
      arcs = Tuple{Int, Int}[]
      for l in nw.links
          l.f > 0 && push!(arcs, (l.i, l.j))
          l.r > 0 && push!(arcs, (l.j, l.i))
      end
      classes = _testhelper_closed_classes(n, arcs)
      length(classes) == 1 || return (v = big(NaN), nclasses = length(classes))
      C = only(classes)
      at = zeros(Int, n)
      at[C] .= eachindex(C)
      inside = [l for l in nw.links if at[l.i] > 0 && at[l.j] > 0]
      R = zeros(BigFloat, length(C), length(C))
      for l in inside
          rf, rr = rates(l)
          R[at[l.i], at[l.j]] += rf
          R[at[l.j], at[l.i]] += rr
      end
      p = length(C) == 1 ? BigFloat[1] : _testhelper_gth(R)
      two_way(l) = l.fast && l.f > 0 && l.r > 0
      h = _testhelper_offsets(n, filter(two_way, inside))
      v = big(0.0)
      for l in inside
          rf, rr = rates(l)
          k = h === nothing ? l.w : two_way(l) ? 0 : l.w + h[l.i] - h[l.j]
          k == 0 || (v += k * (rf * p[at[l.i]] - rr * p[at[l.j]]))
      end
      (v = v, nclasses = 1)
  end

  """
  The rapid-equilibrium limit (Λ = ∞) of `nw`, exactly, at an interior point: the forms of
  each segment of fast links take their equilibrium weights, the segments form the chain of
  the slow links, and the turnover sums the slow links as `_testhelper_steady_rate` does.
  `Inf` when the fast links close a cycle that takes the first substrate up, NaN when the
  slow chain has more than one closed class.
  """
  function _testhelper_re_limit(nw::_testhelper_Network)
      n = length(nw.nodes)
      fast = [l for l in nw.links if l.fast]
      all(l -> l.f > 0 && l.r > 0, fast) ||
          error("_testhelper_re_limit: a fast link has a zero rate; the limit is taken at " *
                "interior points")
      h = _testhelper_offsets(n, fast)
      h === nothing && return big(Inf)
      segment, weight = zeros(Int, n), zeros(BigFloat, n)
      adjacent = [Tuple{Int, BigFloat}[] for _ in 1:n]
      for l in fast
          push!(adjacent[l.i], (l.j, l.f / l.r))
          push!(adjacent[l.j], (l.i, l.r / l.f))
      end
      ns = 0
      for s in 1:n
          segment[s] == 0 || continue
          ns += 1
          segment[s], weight[s] = ns, big(1)
          todo = [s]
          while !isempty(todo)
              u = pop!(todo)
              for (v, q) in adjacent[u]
                  segment[v] == 0 || continue
                  segment[v], weight[v] = ns, weight[u] * q
                  push!(todo, v)
              end
          end
      end
      W = zeros(BigFloat, ns)
      foreach(i -> W[segment[i]] += weight[i], 1:n)
      R, arcs = zeros(BigFloat, ns, ns), Tuple{Int, Int}[]
      for l in nw.links
          l.fast && continue
          a, b = segment[l.i], segment[l.j]
          a == b && continue
          rf, rr = l.f * weight[l.i] / W[a], l.r * weight[l.j] / W[b]
          rf > 0 && (R[a, b] += rf; push!(arcs, (a, b)))
          rr > 0 && (R[b, a] += rr; push!(arcs, (b, a)))
      end
      classes = _testhelper_closed_classes(ns, arcs)
      length(classes) == 1 || return big(NaN)
      C = only(classes)
      P = zeros(BigFloat, ns)
      P[C] .= length(C) == 1 ? BigFloat[1] : _testhelper_gth(R[C, C])
      p = [P[segment[i]] * weight[i] / W[segment[i]] for i in 1:n]
      v = big(0.0)
      for l in nw.links
          l.fast && continue
          k = l.w + h[l.i] - h[l.j]
          k == 0 || (v += k * (l.f * p[l.i] - l.r * p[l.j]))
      end
      v
  end

  """The relative gap between `x` and `y`, with an absolute floor `least`."""
  _testhelper_relgap(x, y; least = big(10)^-14) = abs(x - y) / max(abs(x), abs(y), least)

  """The package's rate law of the concrete mechanism `mech` at the constants of `point`
  (its fitted values and `Keq`) and the concentrations `concs`, read without the boundary
  check (`_testhelper_law_value`); NaN where it has no value."""
  function _testhelper_package_law(mech, point, concs::AbstractDict)
      params = Dict{Symbol, BigFloat}(k => point.vals[k]
                                      for k in EnzymeRates._dependent_param_exprs(mech)[2])
      params[:Keq] = point.Keq
      _testhelper_law_value(mech, params, concs)
  end

  """
  The law's value on the face of `mech` where the names `face` vanish, computed without the
  package: the rapid-equilibrium limit of the law's own network (`_testhelper_re_limit`) at
  the interior points `face = δ·u`, δ = 1e-40, 1e-80, 1e-160, along two random directions
  `u`; NaN when it grows without bound, has not settled by δ = 1e-160 or depends on `u` (the
  law's 0/0 does not cancel). With `face` empty, the limit at `concs`.
  """
  function _testhelper_limit_law(mech, point, concs::AbstractDict, face::Set{Symbol}; rng)
      speeds = _testhelper_draw_speeds(mech, rng)
      value(c) = _testhelper_re_limit(_testhelper_network(mech, point, c; speeds))
      inner = value(concs)
      isempty(face) && return inner
      tiny = big(10)^-30 * max(abs(inner), big(10)^-40)
      limits = BigFloat[]
      for _ in 1:2
          u = Dict(z => big(10)^(2 * big(rand(rng)) - 1) for z in face)
          vals = [value(merge(concs, Dict(z => δ * u[z] for z in face)))
                  for δ in (big(10)^-40, big(10)^-80, big(10)^-160)]
          all(isfinite, vals) || return big(NaN)
          abs(vals[3]) > 1e10 * max(abs(vals[2]), tiny) && return big(NaN)
          _testhelper_relgap(vals[2], vals[3]; least = tiny) > big(10)^-25 && return big(NaN)
          push!(limits, vals[3])
      end
      _testhelper_relgap(limits[1], limits[2]; least = tiny) > big(10)^-20 ? big(NaN) :
          limits[1]
  end

  """
      _testhelper_respeed(m, point, concs, face; rng, law = :package, Λs, npat = 2) -> Symbol

  The RE-speed oracle (spec section 14) on the face of `m` where the names `face` vanish:
  the law's network (`_testhelper_network`, one subunit) at each Λ of `Λs` under `npat`
  speed patterns of 10^±3, its rapid-equilibrium edges at Λ times their pattern's factor.
  The first verdict that applies:
  - `:nonunique`: the face chain has more than one closed class;
  - `:diverges`: |v| grows more than 1e3-fold from the first Λ to the last (RE turnover);
  - `:unconverged`: the change from the first Λ to the last, scaled by their ratio, exceeds
    1e-10 relative, so the last Λ is too small to judge 1e-8;
  - `:pattern`: the values at the last Λ differ between patterns by more than 1e-8;
  - `:law_nan`: the law has no value on the face;
  - `:mismatch`: the limit differs from the law by more than 1e-8 relative;
  - `:agree`.
  The law is the package's (`law = :package`, read unchecked) or the oracle's own
  (`law = :limit`, `_testhelper_limit_law`).
  """
  function _testhelper_respeed(m, point, concs::AbstractDict, face::Set{Symbol}; rng,
                               law::Symbol = :package, Λs = (big(10.0)^10, big(10.0)^20),
                               npat::Int = 2)
      mech = _testhelper_lift(m)
      c0 = merge(concs, Dict(z => big(0) for z in face))
      patterns = [_testhelper_draw_speeds(mech, rng) for _ in 1:npat]
      solved = [_testhelper_steady_rate(_testhelper_network(mech, point, c0; speeds); Λ)
                for speeds in patterns, Λ in Λs]
      any(x -> x.nclasses > 1, solved) && return :nonunique
      v = [x.v for x in solved]
      top = size(v, 2)
      any(i -> abs(v[i, top]) > 1e3 * max(abs(v[i, 1]), big(10)^-30), axes(v, 1)) &&
          return :diverges
      any(i -> abs(v[i, 1] - v[i, top]) * Λs[1] / Λs[end] >
               1e-10 * max(abs(v[i, top]), big(10)^-14), axes(v, 1)) && return :unconverged
      any(i -> _testhelper_relgap(v[i, top], v[1, top]) > 1e-8, axes(v, 1)) &&
          return :pattern
      vlaw = law === :package ? _testhelper_package_law(mech, point, c0) :
             _testhelper_limit_law(mech, point, concs, face; rng)
      isfinite(vlaw) || return :law_nan
      _testhelper_relgap(v[1, top], vlaw) > 1e-8 ? :mismatch : :agree
  end

  """
  The RE-speed oracle's verdict on each face of `m` (`_testhelper_faces`), the worst over
  `npoints` points (`_testhelper_consistent_point`, concentrations log-uniform over
  10^±1) in the order nonunique, diverges, unconverged, pattern, law_nan, mismatch, agree.
  An allosteric mechanism is judged at one catalytic subunit (P4-24a,
  `_testhelper_single_subunit`).
  """
  function _testhelper_respeed_verdicts(m; rng, npoints::Int = 1, law::Symbol = :package)
      mech = _testhelper_single_subunit(_testhelper_lift(m))
      faces = _testhelper_faces(mech)
      order = (:nonunique, :diverges, :unconverged, :pattern, :law_nan, :mismatch, :agree)
      worst = fill(:agree, length(faces))
      setprecision(BigFloat, 256) do
          for _ in 1:npoints
              point = _testhelper_consistent_point(mech, rng)
              concs = Dict{Symbol, BigFloat}(x => big(10)^(2 * big(rand(rng)) - 1)
                                             for x in EnzymeRates.metabolites(mech))
              for (i, face) in enumerate(faces)
                  v = _testhelper_respeed(mech, point, concs, face; rng, law)
                  findfirst(==(v), order) < findfirst(==(worst[i]), order) &&
                      (worst[i] = v)
              end
          end
      end
      worst
  end
  ```

- [ ] **Step 9: Run them and watch them pass.** The Step 3 command, log `P4_T1_green2.log`.
  Expected: `exit 0`; the three new testsets pass (3, 2 and 1 tests); the earlier testsets of
  the file, among them "Exact steady-state oracle" and "Laws over a regenerated population",
  still pass (they read the refactored `_testhelper_steady_state_rate` and the extended
  point); no warning.

- [ ] **Step 10: Write the failing tests of the ε-path, the LP, the corners and the
  n-protomer oracle.**

  Append to `test/test_oracles.jl`:
  ```julia
  @testset "the ε-path oracle rejects islands, and R2 only at P = 0" begin
      # R3 and R4 strand F behind tagged chemistry in the inactive conformation, an island:
      # its weight follows the potential the path gives it. The research ε-path read two
      # orders; this oracle reads five, so an island may read direction-dependent or a
      # wrong limit.
      r3 = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)      :: EqualAI
              E(S) <--> E(P)    :: EqualAI
              E + P ⇌ E(P)      :: EqualAI
              E <--> F          :: OnlyA
              F + I ⇌ F(I)      :: EqualAI
          end
      end
      r4 = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)      :: EqualAI
              E(S) <--> E(P)    :: NonequalAI
              E + P ⇌ E(P)      :: EqualAI
              F ⇌ E             :: OnlyA
          end
      end
      @test _testhelper_epath_verdict(r3; rng = Random.MersenneTwister(31)) in
            (:direction_dependent, :wrong_limit)
      @test _testhelper_epath_verdict(r4; rng = Random.MersenneTwister(32)) in
            (:direction_dependent, :wrong_limit)
      r2 = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)      :: OnlyA
              E(S) <--> E(P)    :: EqualAI
              F + P ⇌ E(P)      :: EqualAI
              F <--> E          :: OnlyA
          end
      end
      @test _testhelper_epath_verdict(r2; rng = Random.MersenneTwister(33)) === :converges
      @test _testhelper_epath_verdict(r2; rng = Random.MersenneTwister(34),
                                      face = Set([:P])) === :direction_dependent
      @test _testhelper_epath_verdict(r2; rng = Random.MersenneTwister(35),
                                      face = Set([:S])) === :converges
      # An exclusion that no direction can reach: no path.
      no_path = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)      :: OnlyA
              E(S) <--> E(P)    :: EqualAI
              E + P ⇌ E(P)      :: EqualAI
          end
      end
      @test _testhelper_epath_verdict(no_path; rng = Random.MersenneTwister(36)) === :no_path
  end

  @testset "an admitted law converges inside the ε window" begin
      # The inactive conformation's chemistry is blocked: the error falls as ε^σ.
      chem = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)      :: EqualAI
              E(S) <--> E(P)    :: OnlyA
              E(P) ⇌ E + P      :: EqualAI
          end
      end
      @test _testhelper_in_epath_window(chem; rng = Random.MersenneTwister(41))
      @test _testhelper_epath_verdict(chem; rng = Random.MersenneTwister(42)) === :converges
  end

  @testset "the test's certified direction satisfies every cycle row" begin
      grid = _testhelper_uniuni_grid()
      @test length(grid) == 108
      certified = 0
      for am in grid
          ex = _testhelper_exclusion(am, 2)
          @test _testhelper_certified(am, 2) == (ex.failure === :none)
          certified += ex.failure === :none
      end
      # 44 tag sets with an :OnlyA group have a direction, and the 32 without one.
      @test certified == 76
  end

  @testset "the exact LP agrees with the test's graph test on the uni-uni grid" begin
      grid = _testhelper_uniuni_grid()
      @test count(am -> _testhelper_lp_existence(am, 2) ==
                        (_testhelper_exclusion(am, 2).failure === :none), grid) == 108
      @test count(am -> _testhelper_lp_existence(am, 2), grid) == 76
  end

  @testset "the corner set holds every corner up to 128" begin
      rng = Random.MersenneTwister(51)
      seven = _testhelper_corners(7, rng)
      @test length(seven) == 128
      @test Set(seven) == Set(digits(x; base = 2, pad = 7) for x in 0:127)
      eight = _testhelper_corners(8, rng)
      @test length(eight) == 128 && length(Set(eight)) == 128
      @test all(c -> length(c) == 8 && all(in((0, 1)), c), eight)
      @test _testhelper_corners(0, rng) == [Int[]]
  end

  @testset "the n-protomer oracle equals the one-subunit solver at n = 1" begin
      am = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          allosteric_regulators: R
          steps: begin
              E + S ⇌ E(S)      :: NonequalAI
              E(S) <--> E(P)    :: OnlyA
              E(P) ⇌ E + P      :: EqualAI
          end
          regulatory_site(2): R::NonequalAI
      end
      setprecision(BigFloat, 256) do
          rng = Random.MersenneTwister(61)
          point = _testhelper_consistent_point(am, rng)
          concs = Dict{Symbol, BigFloat}(x => big(10)^(2 * big(rand(rng)) - 1)
                                         for x in EnzymeRates.metabolites(am))
          exact = _testhelper_steady_state_rate(am, point, concs)
          @test abs(_testhelper_protomer_rate(am, point, concs; n = 1) - exact) <=
                big(10)^-20 * abs(exact)
          # Two subunits, each at one of a conformation's three law forms: nine states per
          # conformation.
          two = _testhelper_protomer_network(am, point, concs; n = 2,
                                             speeds = _testhelper_draw_speeds(am, rng))
          @test length(two.nodes) == 3^2 + 3^2
      end
  end
  ```

- [ ] **Step 11: Run them and watch them fail.** The Step 3 command, log `P4_T1_red3.log`.
  Expected: the six new testsets report `Error During Test` with `UndefVarError` for
  `_testhelper_epath_verdict`, `_testhelper_in_epath_window`, `_testhelper_uniuni_grid`,
  `_testhelper_corners` and `_testhelper_protomer_rate`; every other testset passes.

- [ ] **Step 12: Add the ε-path oracle, the window, the LP, the corners, the grid and the
  n-protomer oracle.** Append to the "Phase-4 oracles" section:
  ```julia
  """
      _testhelper_exclusion(am, c) -> (; component, arcs, failure, lost)

  The exclusion structure of conformation `c` of `am` (spec section 7): `component[f]`, the
  component of form `f` under the edges `c` keeps, numbered by first form; `arcs`,
  `(unbound side's component, bound side's, edge)` for each plain binding `c` excludes;
  `failure`, `:inside` when an arc stays inside one component, `:cycle` when the arcs close
  a directed cycle, else `:none` (a path to the limit exists); `lost`, the components free
  E's does not reach along the arcs (islands).
  """
  function _testhelper_exclusion(am, c::Int)
      cm = EnzymeRates.catalytic(am)
      nf = length(cm.forms)
      parent = collect(1:nf)
      root(i) = parent[i] == i ? i : root(parent[i])
      for (k, e) in enumerate(cm.edges)
          _testhelper_excludes(am, k, c) || (parent[root(e.from)] = root(e.to))
      end
      roots = [root(f) for f in 1:nf]
      label = Dict(r => q for (q, r) in enumerate(unique(roots)))
      component = [label[r] for r in roots]
      n = length(label)
      arcs = Tuple{Int, Int, Int}[(component[e.from], component[e.to], k)
                                  for (k, e) in enumerate(cm.edges)
                                  if _testhelper_excludes(am, k, c) &&
                                     EnzymeRates._is_plain_binding(cm, e)]
      failure = any(a -> a[1] == a[2], arcs) ? :inside : :none
      if failure === :none
          indegree = zeros(Int, n)
          foreach(a -> indegree[a[2]] += 1, arcs)
          todo = findall(iszero, indegree)
          while !isempty(todo)
              u = pop!(todo)
              for a in arcs
                  a[1] == u && (indegree[a[2]] -= 1) == 0 && push!(todo, a[2])
              end
          end
          any(>(0), indegree) && (failure = :cycle)
      end
      reached = falses(n)
      reached[component[1]] = true
      grew = true
      while grew
          grew = false
          for a in arcs
              reached[a[1]] && !reached[a[2]] && (reached[a[2]] = true; grew = true)
          end
      end
      (; component, arcs, failure, lost = findall(!, reached))
  end

  """The excluding conformations of `am` with their exclusion structures, `[(c, ex)]`."""
  function _testhelper_excluding(am)
      edges = eachindex(EnzymeRates.catalytic(am).edges)
      [(c, _testhelper_exclusion(am, c)) for c in 1:2
       if any(k -> _testhelper_excludes(am, k, c), edges)]
  end

  """Potentials of the components of `ex` (`_testhelper_exclusion` with no failure): each
  starts at `start[x]`, each arc raises its bound side to at least its unbound side's plus
  `rise[q]`, and free E's component is shifted back to 0."""
  function _testhelper_component_potentials(ex, start::Vector{BigFloat},
                                            rise::Vector{BigFloat})
      φ = copy(start)
      for _ in eachindex(φ), (q, a) in enumerate(ex.arcs)
          φ[a[2]] = max(φ[a[2]], φ[a[1]] + rise[q])
      end
      φ .- φ[ex.component[1]]
  end

  """
  The ε-path's six standard directions for `am` (research section 3.1, P4-25): integer
  potentials, arcs rising by 1 and every component other than free E's starting at +1 or
  at -1, each with one speed exponent σ ∈ {½, 1, 2} for every blocked edge and the even
  split θ = ½ of each vanishing rapid-equilibrium binding. Each is `(; potentials, sigma,
  theta)`, `potentials[c][f]` the potential of form `f` in conformation `c` (empty where `c`
  excludes nothing).
  """
  function _testhelper_standard_directions(am)
      ne = length(EnzymeRates.catalytic(am).edges)
      directions = []
      for offset in (1, -1), σ in (big(0.5), big(1.0), big(2.0))
          potentials = [BigFloat[] for _ in 1:2]
          for (c, ex) in _testhelper_excluding(am)
              start = [x == ex.component[1] ? big(0.0) : BigFloat(offset)
                       for x in 1:maximum(ex.component)]
              φ = _testhelper_component_potentials(ex, start, ones(BigFloat, length(ex.arcs)))
              potentials[c] = [φ[x] for x in ex.component]
          end
          push!(directions, (; potentials, sigma = fill(σ, ne), theta = fill(big(0.5), ne)))
      end
      directions
  end

  """The corners of a box of `k` binary choices, each a vector of `k` zeros and ones: all
  `2^k` when there are at most 128, else 128 distinct ones drawn from `rng` (P4-24a)."""
  function _testhelper_corners(k::Int, rng)
      k == 0 && return [Int[]]
      k <= 7 && return [digits(x; base = 2, pad = k) for x in 0:(2^k - 1)]
      drawn = Set{Vector{Int}}()
      while length(drawn) < 128
          push!(drawn, rand(rng, 0:1, k))
      end
      sort!(collect(drawn))
  end

  """
  The ε-path's corner directions for `am` (P4-24a): the corners of the admissible cone,
  every binary choice made one way or the other. Every arc rises by 1 or 20, every component
  other than free E's that no arc enters starts at -5 or 5, every blocked edge's σ is 1/8 or
  16 and every vanishing rapid-equilibrium binding's θ is 0 or 1; all corners when there
  are at most 128, else 128 of them (`_testhelper_corners`).
  """
  function _testhelper_corner_directions(am, rng)
      cm = EnzymeRates.catalytic(am)
      excluding = _testhelper_excluding(am)
      slots = Tuple{Symbol, Int, Int}[]
      for (c, ex) in excluding
          append!(slots, [(:rise, c, q) for q in eachindex(ex.arcs)])
          entered = Set(a[2] for a in ex.arcs)
          append!(slots, [(:start, c, x) for x in 1:maximum(ex.component)
                          if x != ex.component[1] && !(x in entered)])
      end
      for (k, e) in enumerate(cm.edges)
          any(c -> _testhelper_excludes(am, k, c), 1:2) || continue
          if !EnzymeRates._is_plain_binding(cm, e)
              push!(slots, (:sigma, k, 0))
          elseif e.rapid_equilibrium
              push!(slots, (:theta, k, 0))
          end
      end
      map(_testhelper_corners(length(slots), rng)) do bits
          high(slot) = (q = findfirst(==(slot), slots); q !== nothing && bits[q] == 1)
          potentials = [BigFloat[] for _ in 1:2]
          for (c, ex) in excluding
              entered = Set(a[2] for a in ex.arcs)
              start = [x == ex.component[1] ? big(0.0) : x in entered ? big(-Inf) :
                       high((:start, c, x)) ? big(5.0) : big(-5.0)
                       for x in 1:maximum(ex.component)]
              rise = [high((:rise, c, q)) ? big(20.0) : big(1.0) for q in eachindex(ex.arcs)]
              φ = _testhelper_component_potentials(ex, start, rise)
              potentials[c] = [φ[x] for x in ex.component]
          end
          edges = eachindex(cm.edges)
          (; potentials,
           sigma = [high((:sigma, k, 0)) ? big(16.0) : big(0.125) for k in edges],
           theta = [high((:theta, k, 0)) ? big(1.0) : big(0.0) for k in edges])
      end
  end

  """How far the ε-network of `am` along direction `d` spreads its rates' ε-exponents: the
  sum over the excluded edges of the gap between their two exponents. The RE-first order
  lets Λ outgrow it."""
  function _testhelper_epsilon_span(am, d)
      cm = EnzymeRates.catalytic(am)
      span = big(0.0)
      for c in 1:2, (k, e) in enumerate(cm.edges)
          _testhelper_excludes(am, k, c) || continue
          rise = abs(d.potentials[c][e.to] - d.potentials[c][e.from])
          span += EnzymeRates._is_plain_binding(cm, e) ? rise : 2 * d.sigma[k] + rise
      end
      span
  end

  """
  The limit of the turnover of the ε-network `network(ε)` along Λ = ε^-κ as ε → 0, ε =
  10^-E with E doubling: `(; v, converged)`. κ is `:re_first` (Λ outgrows every ε-scale,
  `2 span + 10`), `:eps_first` (Λ grows slower than any rate vanishes, `1/(8(nodes + 1))`)
  or a number. Converged once Λ ≥ 1e40 and the last two steps move by less than 1e-12 and
  1e-6 of `scale`; infinite when the rate keeps growing by orders of magnitude.
  """
  function _testhelper_path_limit(network, κ, span::BigFloat, nnodes::Int, scale::BigFloat)
      k = κ === :re_first ? 2 * span + 10 :
          κ === :eps_first ? big(1) / (8 * (nnodes + 1)) : BigFloat(κ)
      E, Emax = κ === :eps_first ? (4, 2^18) : (2, 2^14)
      vals = BigFloat[]
      while E <= Emax
          ε = big(10.0)^(-E)
          push!(vals, _testhelper_steady_rate(network(ε); Λ = ε^(-k)).v)
          n = length(vals)
          if n >= 3
              d1, d2 = abs(vals[n] - vals[n - 1]), abs(vals[n - 1] - vals[n - 2])
              k * E >= 40 && d1 <= big(10)^-12 * scale && d2 <= big(10)^-6 * scale &&
                  return (v = vals[n], converged = true)
              abs(vals[n]) > 1e6 * max(abs(vals[n - 1]), scale) &&
                  abs(vals[n - 1]) > 1e3 * max(abs(vals[n - 2]), scale) &&
                  return (v = big(Inf), converged = true)
          end
          E *= 2
      end
      (v = vals[end], converged = false)
  end

  """
      _testhelper_epath_verdict(m; rng, face = Set{Symbol}(), corners = false,
                                precision = 256) -> Symbol

  The ε-path oracle (spec section 14; P4-7, P4-24a) on the face of `m` where the names
  `face` vanish, one catalytic subunit: along each direction (the six standard ones, or the
  corners with `corners`) and each order κ ∈ (RE first, 8, 2, ½, ε first) of the two
  limits, the turnover of the ε-network against the deletion law's value there
  (`_testhelper_limit_law`). The first verdict that applies:
  - `:no_path`: an excluding conformation has no path to the limit (`_testhelper_exclusion`);
  - `:nonunique`: an ε-network has more than one closed class on the face;
  - `:diverges`: some path's turnover grows without bound (a fast cycle of the ε-network
    turns over);
  - `:del_nonunique`: the deletion law has no value on the face;
  - `:converges`: every path tends to it, to 1e-8 of the rate's scale;
  - `:wrong_limit`: every path settles, none on it, all on one other value;
  - `:direction_dependent`: otherwise (some paths tend to it and some do not, or some never
    settle).
  """
  function _testhelper_epath_verdict(m; rng, face::Set{Symbol} = Set{Symbol}(),
                                     corners::Bool = false, precision::Int = 256)
      setprecision(BigFloat, precision) do
          am = _testhelper_single_subunit(_testhelper_lift(m))
          any(((_, ex),) -> ex.failure !== :none, _testhelper_excluding(am)) &&
              return :no_path
          point = _testhelper_consistent_point(am, rng)
          concs = Dict{Symbol, BigFloat}(x => big(10)^(2 * big(rand(rng)) - 1)
                                         for x in EnzymeRates.metabolites(am))
          c0 = merge(concs, Dict(z => big(0) for z in face))
          vdel = _testhelper_limit_law(am, point, concs, face; rng)
          vin = _testhelper_limit_law(am, point, concs, Set{Symbol}(); rng)
          scale = max(isfinite(vdel) ? abs(vdel) : big(0), isfinite(vin) ? abs(vin) : big(1),
                      big(10)^-30)
          directions = corners ? _testhelper_corner_directions(am, rng) :
                       _testhelper_standard_directions(am)
          speeds = _testhelper_draw_speeds(am, rng)
          outcomes, ends = Symbol[], BigFloat[]
          for d in directions
              function _testhelper_eps_network(ε)
                  _testhelper_network(am, point, c0; speeds, mode = :vanish, ε,
                                      potentials = d.potentials, sigma = d.sigma,
                                      theta = d.theta)
              end
              nw = _testhelper_eps_network(big(10)^-2)
              _testhelper_steady_rate(nw; Λ = big(1.0)).nclasses > 1 && return :nonunique
              span = _testhelper_epsilon_span(am, d)
              for κ in (:re_first, 8, 2, 1 // 2, :eps_first)
                  r = _testhelper_path_limit(_testhelper_eps_network, κ, span,
                                             length(nw.nodes), scale)
                  if !r.converged
                      push!(outcomes, :unsettled)
                  elseif !isfinite(r.v)
                      push!(outcomes, :infinite)
                  elseif abs(r.v - vdel) <= big(10)^-8 * scale
                      push!(outcomes, :tends)
                  else
                      push!(outcomes, :elsewhere)
                      push!(ends, r.v)
                  end
              end
          end
          :infinite in outcomes && return :diverges
          isfinite(vdel) || return :del_nonunique
          all(==(:tends), outcomes) && return :converges
          all(==(:elsewhere), outcomes) &&
              all(x -> abs(x - ends[1]) <= big(10)^-6 * scale, ends) && return :wrong_limit
          :direction_dependent
      end
  end

  """
  Whether every standard direction (`_testhelper_standard_directions`) of `m`'s ε-path at an
  interior point approaches the deletion law at a power of ε (P4-25): over ε = 1e-9 … 1e-18
  the error relative to the rate is `C·ε^p` within 10% at every decade, `p` a positive
  multiple of ½, with exact rapid-equilibrium limits (`_testhelper_re_limit`) in 512-bit
  BigFloat; an error below 1e-60 throughout passes.
  """
  function _testhelper_in_epath_window(m; rng)
      setprecision(BigFloat, 512) do
          am = _testhelper_single_subunit(_testhelper_lift(m))
          point = _testhelper_consistent_point(am, rng)
          concs = Dict{Symbol, BigFloat}(x => big(10)^(2 * big(rand(rng)) - 1)
                                         for x in EnzymeRates.metabolites(am))
          speeds = _testhelper_draw_speeds(am, rng)
          vdel = _testhelper_re_limit(_testhelper_network(am, point, concs; speeds))
          scale = max(abs(vdel), big(10)^-30)
          eps = [big(10.0)^-j for j in 9:18]
          all(_testhelper_standard_directions(am)) do d
              errors = [abs(_testhelper_re_limit(_testhelper_network(
                            am, point, concs; speeds, mode = :vanish, ε,
                            potentials = d.potentials, sigma = d.sigma,
                            theta = d.theta)) - vdel) / scale for ε in eps]
              all(<(big(10)^-60), errors) && return true
              any(iszero, errors) && return false
              slopes = sort([log10(errors[j] / errors[j + 1]) for j in 1:(length(eps) - 1)])
              p = slopes[(length(slopes) + 1) ÷ 2]
              q = round(2 * p) / 2
              q >= 1 // 2 && abs(p - q) < 0.05 || return false
              C = errors[end] / eps[end]^q
              all(j -> abs(errors[j] / (C * eps[j]^q) - 1) < 0.1, eachindex(eps))
          end
      end
  end

  """Whether `{x ≥ 0 : M x = b}` is nonempty: phase 1 of the simplex method with Bland's
  rule, in exact rational arithmetic."""
  function _testhelper_feasible(M::Matrix{Rational{BigInt}}, b::Vector{Rational{BigInt}})
      m, n = size(M)
      flip = [x < 0 ? -1 : 1 for x in b]
      T = hcat(flip .* M, Matrix{Rational{BigInt}}(LinearAlgebra.I, m, m), flip .* b)
      basis = collect((n + 1):(n + m))
      cost = [zeros(Rational{BigInt}, n); ones(Rational{BigInt}, m)]
      while true
          reduced(j) = cost[j] - sum(cost[basis[i]] * T[i, j] for i in 1:m)
          j = findfirst(j -> reduced(j) < 0, 1:(n + m))
          j === nothing && break
          rows = [i for i in 1:m if T[i, j] > 0]
          ratio = minimum(i -> T[i, end] / T[i, j], rows)
          i = argmin(i -> basis[i], [i for i in rows if T[i, end] / T[i, j] == ratio])
          T[i, :] ./= T[i, j]
          for k in 1:m
              k == i || iszero(T[k, j]) || (T[k, :] .-= T[k, j] .* T[i, :])
          end
          basis[i] = j
      end
      iszero(sum(cost[basis[i]] * T[i, end] for i in 1:m))
  end

  """
  Admission's existence test as an exact LP (spec section 7): whether conformation `c` of
  `am` has a direction of vanishing exponents for its excluded constants that keeps every
  cycle row of the catalytic graph (`_testhelper_cycle_basis`) with every other constant
  fixed. A plain binding's equilibrium constant must fall (exponent at least 1, after
  scaling); any other excluded edge's exponent is free.
  """
  function _testhelper_lp_existence(am, c::Int)
      cm = EnzymeRates.catalytic(am)
      excluded = [k for k in eachindex(cm.edges) if _testhelper_excludes(am, k, c)]
      plain = [k for k in excluded if EnzymeRates._is_plain_binding(cm, cm.edges[k])]
      isempty(plain) && return true
      rows = _testhelper_cycle_basis((forms = cm.forms, steps = cm.edges))
      isempty(rows) && return true
      Z = Rational{BigInt}.(permutedims(reduce(hcat, rows)))
      AP, AF = Z[:, plain], Z[:, setdiff(excluded, plain)]
      _testhelper_feasible([AP AF -AF], -AP * ones(Rational{BigInt}, length(plain)))
  end

  """
  Whether the test's own certified direction for conformation `c` of `am` exists and keeps
  every cycle row (P4-26): the component potentials of `_testhelper_exclusion` (free E's at
  0, every other source at 1, arcs rising by 1) give each excluded edge the exponent of its
  bound side over its unbound side; each arc's must be at least 1 and every row of
  `_testhelper_cycle_basis` must sum to 0.
  """
  function _testhelper_certified(am, c::Int)
      ex = _testhelper_exclusion(am, c)
      ex.failure === :none || return false
      cm = EnzymeRates.catalytic(am)
      start = [x == ex.component[1] ? big(0.0) : big(1.0) for x in 1:maximum(ex.component)]
      φ = _testhelper_component_potentials(ex, start, ones(BigFloat, length(ex.arcs)))
      a = [_testhelper_excludes(am, k, c) ? φ[ex.component[e.to]] - φ[ex.component[e.from]] :
           big(0.0) for (k, e) in enumerate(cm.edges)]
      rows = _testhelper_cycle_basis((forms = cm.forms, steps = cm.edges))
      all(arc -> a[arc[3]] >= 1, ex.arcs) && all(z -> iszero(sum(z .* a)), rows)
  end

  """The uni-uni grid of admission's existence test: `E + S ⇌ E(S)` or `<-->`,
  `E(S) <--> E(P)` and `E + P ⇌ E(P)` or `<-->`, each step with each of `:OnlyA`,
  `:EqualAI` and `:NonequalAI`, 4 × 27 = 108 allosteric mechanisms, built unguarded."""
  function _testhelper_uniuni_grid()
      tags = (:OnlyA, :EqualAI, :NonequalAI)
      [EnzymeRates.AllostericMechanism(EnzymeRates._raw_mechanism(Meta.parse("""
           begin
               substrates: S
               products: P
               steps: begin
                   E + S $s E(S) :: $tS
                   E(S) <--> E(P) :: $tC
                   E + P $p E(P) :: $tP
               end
           end"""), true))
       for s in ("⇌", "<-->") for p in ("⇌", "<-->")
       for tS in tags for tC in tags for tP in tags]
  end

  """
  The concerted network of `n` catalytic subunits (research section 3.1, P4-24a): each
  subunit moves through its conformation's law network (`_testhelper_network`) on its own,
  and the oligomer flips only with every subunit at free E, at the single-subunit flip
  ratio. Its nodes are `(state, conformation)`, a state numbering one assignment of the law
  network's nodes to the subunits.
  """
  function _testhelper_protomer_network(m, point, concs::AbstractDict; n::Int,
                                        speeds::Matrix{BigFloat})
      single = _testhelper_network(m, point, concs; speeds)
      conformation(u) = single.nodes[u][2]
      states = Vector{Int}[]
      for c in sort!(unique(last.(single.nodes)))
          members = [u for u in eachindex(single.nodes) if conformation(u) == c]
          append!(states, [collect(t) for t in Iterators.product(ntuple(_ -> members, n)...)])
      end
      index = Dict(s => q for (q, s) in enumerate(states))
      links = _testhelper_Link[]
      for (q, s) in enumerate(states), p in 1:n, l in single.links
          conformation(l.i) == conformation(l.j) && l.i == s[p] || continue
          t = copy(s)
          t[p] = l.j
          push!(links, _testhelper_Link(q, index[t], l.f, l.r, l.w, l.fast))
      end
      for l in single.links
          conformation(l.i) == conformation(l.j) && continue
          push!(links, _testhelper_Link(index[fill(l.i, n)], index[fill(l.j, n)], l.f, l.r,
                                        0, true))
      end
      _testhelper_Network([(q, conformation(s[1])) for (q, s) in enumerate(states)], links)
  end

  """The steady-state turnover per active site of the concerted `n`-subunit network of `m`
  (`_testhelper_protomer_network`), its rapid-equilibrium edges at `Λ` times `speeds`."""
  function _testhelper_protomer_rate(m, point, concs::AbstractDict; n::Int,
                                     speeds::Matrix{BigFloat} = fill(big(1.0), size(point.K)),
                                     Λ::BigFloat = big(10.0)^30)
      network = _testhelper_protomer_network(m, point, concs; n, speeds)
      _testhelper_steady_rate(network; Λ).v / n
  end
  ```

- [ ] **Step 13: Run them and watch them pass.** The Step 3 command, log `P4_T1_green3.log`.
  Expected: `exit 0`; the six testsets pass (6, 2, 110, 2, 5 and 2 tests); no warning. Note the
  times of the eleven new testsets in the summary; their sum stays at most 4 s, else profile
  the ε-path calls (the largest part) and say so before the commit. The refresh pass
  (P4-32) ran Task 1's code at 08c4043 (`$S/phase45/refresh/map-head-t0-3_task1.py`, log
  `map-head-t0-3_6.log`): the eleven testsets passed with these counts, the grid's 76 of
  108 tag sets with a direction and the certified-direction testset's 110 tests among them;
  this run repeats that check on the branch. A pin that fails is checked against the
  research LP (`slice/allo/g1_lp.jl`) on that grid member before any change, and is never
  re-pinned to the code's output.

- [ ] **Step 14: Run the full suite.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T1_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'
  ```
  In the background; poll. Expected: `exit 0`, `Testing EnzymeRates tests passed`, the pass
  count N0 plus the new tests (about 140), Aqua and JET pass, no warning; `maxrss_kb` within
  5% of M0 (the oracles hold no population; a larger rise is explained in the commit body).

- [ ] **Step 15: Commit.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && git add test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_oracles.jl && git commit -F - <<'EOF'
  Add the RE-speed, ε-path, LP and n-protomer test oracles

  The guard's oracles, read from the graph read API: networks of one catalytic subunit
  flipping at free E, the RE-speed verdicts on every face, the ε-path along integer
  potentials or the cone's corners in five orders of the limits, the exact LP of
  admission's existence test and the concerted n-protomer network. Both potential kinds
  rise by integers (P4-25) and κ is P4-24a's set, where the research scripts jittered one
  kind and ran seven orders.

  Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
  EOF
  ```
  Append the full run's pass count, time and maxrss to the body if maxrss rose more than 5%.

---

### Task 2: Networks, faces and classes 1-2, not wired

The guard's verdict type, its two networks (RES §1.1), the face list (§1.2), classes 1 and 2
(§1.3), step lines and the class 1-2 messages (§1.11). Nothing calls `_guard` yet; Task 8
wires it.

**Files:**
- Modify: `src/admission.jl` (new sections "The guard: verdicts, networks, faces, classes 1
  and 2" and "Guard messages", after `_derivation_input`)
- Modify: `src/types.jl` (sub-graph section, after `_inactive_mask`: `_excludes`,
  `_law_mask`)
- Modify: `ROADMAP.md` (Known issues: the mixed-corner entry, P4-12a)
- Create: `test/test_admission.jl`; Modify: `test/runtests.jl`, `test/test_oracles.jl` (one
  testset, beside `_testhelper_population_sample`)

**Interfaces:**
- Consumes: Task 1's `@_testhelper_unguarded`, `_testhelper_faces`, `_testhelper_respeed`,
  `_testhelper_respeed_verdicts`, `_testhelper_network`, `_testhelper_steady_rate`,
  `_testhelper_consistent_point`, `_testhelper_draw_speeds`; phase 3's `_flags`,
  `_group_members`, `_name_index`, `_is_plain_binding`, `_is_reactant`, `_form_text` and
  `_step_text` (names.jl),
  `entries`, `copies`, `_raw(m)`, `_raw_step` and `_step_parts` (dsl.jl), `modifiers`,
  `allosteric_regulators`, `_has_unbalanced_cycle` (mechanism_enumeration.jl),
  `_re_turnover_cycle(m, fast)` (mechanism_enumeration.jl); `_testhelper_edge`,
  `_testhelper_lift`, `_testhelper_random_independent_params`,
  `_testhelper_population_sample`.
- Produces:
  ```julia
  @enum _Fault::UInt8 _NO_FAULT _NO_CYCLE _RE_TURNOVER          # Tasks 3-6 append values
  @enum _LimitVerdict::UInt8 _LIMIT_UNTESTED _LIMIT_EQUAL _LIMIT_DIFFER _LIMIT_SPEED _LIMIT_LAW_NAN
  struct _Verdict
      fault::_Fault; limit::_LimitVerdict; face::Int; conformation::Int; edge::Int
  end
  const _WELL_POSED = _Verdict(_NO_FAULT, _LIMIT_UNTESTED, 0, 0, 0)
  struct _Network
      nforms::Int; ncopies::Int; present::BitVector; from::Vector{Int}; to::Vector{Int}
      edge::Vector{Int}; conformation::Vector{Int}; forward::BitVector; fast::BitVector
      vanishing::BitVector; change::Matrix{Int}; takes::Vector{UInt64}
  end   # transitions 2j - 1 (stored direction) and 2j (reverse) form a pair
  _node(net::_Network, f::Int, c::Int)::Int
  _excludes(am::AllostericMechanism, k::Int, c::Int)::Bool                  # types.jl
  _law_mask(am::AllostericMechanism, c::Int)::BitVector                     # types.jl
  _law_network(m::Mechanism, re::BitVector = _flags(m))::_Network
  _law_network(am::AllostericMechanism, re::BitVector = _flags(am.catalytic))::_Network
  _eps_network(m::Mechanism, re::BitVector = _flags(m))::_Network
  _eps_network(am::AllostericMechanism, re::BitVector = _flags(am.catalytic))::_Network
  _faces(m::Mechanism)::Vector{UInt64}                 # faces[1] == 0, the interior
  _face_names(m::Mechanism, face::UInt64)::Vector{Symbol}
  _turnover_name(m::Mechanism)::Int
  _has_turnover_cycle(m::Mechanism, mask::BitVector)::Bool
  _re_turnover_transitions(net::_Network)::Vector{Int}  # in cycle order; empty when none
  _re_turnover_verdict(net::_Network)::_Verdict         # class 2's verdict on `net`
  _any_step_cures(mech::Union{Mechanism, AllostericMechanism}, net::_Network,
                  cycle::Vector{Int})::Bool               # class 2's fix holds
  _guard(m::Mechanism)::_Verdict
  _guard(am::AllostericMechanism)::_Verdict
  _guard_message(mech::Union{Mechanism, AllostericMechanism}, v::_Verdict)::String
  ```

**Gate:** RES §3.2's class 1-2 rows: `no_turnover_path` `_NO_CYCLE` with no fix;
`re_turnover_beside_ss`, `re_turnover_flip2_453` and both mechanisms of TRD "All-RE catalytic
cycle raises" `_RE_TURNOVER`, the RE-speed oracle `:diverges` at the interior of each, each
fix's variant well-posed; the four fitted classes pass. `_faces` equals `_testhelper_faces` on
all 42 specs; bi-bi 7, bi-bi with I 14, ter-ter 9. Class 2 on `_eps_network` equals
`_re_turnover_cycle` on the plain specs and the plain population members, with RE overrides
that make it fire. Every step line of every spec reads back as its edge. Cost of classes 1-2
reported per mechanism (prototype 5.7 µs). Full suite green.

- [ ] **Step 0: Check the names this task reads at the branch point.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && grep -n "function _re_turnover_cycle\|function _has_unbalanced_cycle" src/mechanism_enumeration.jl && grep -n "function _form_text\|function _step_text\|function _raw_step\|function _step_parts" src/names.jl src/dsl.jl && grep -rn '::Mod"' src/ | head
  ```
  Expected: `function _has_unbalanced_cycle(from, to, weight, on, parent, offset)` and
  `function _re_turnover_cycle(m::Mechanism, fast::BitVector)` (it returns a `Bool`) in
  `src/mechanism_enumeration.jl`; `function _form_text(m::Mechanism, f::Int)` and
  `function _step_text(m::Mechanism, k::Int)` in `src/names.jl`;
  `_raw_steps(steps_block, names, macro_name, allosteric::Bool)`,
  `_step_parts(expr, macro_name)` and `_raw_step(op, lhs, rhs, names, macro_name)` in
  `src/dsl.jl`; the third grep finds `_met_text` in `src/names.jl`, which writes a
  reactant's copy `P::Mod` for `_form_text` and `_step_text`, and the decodability check's
  spellings in `src/types.jl`. The messages name a form by `_form_text`, which prints a free
  form bare (`E`), and a step by `_step_text`, whose line `show` prints and the macros read
  back as that edge. If a signature differs, call it as the branch point has it in Step 1's
  tests and say so in the commit body.

- [ ] **Step 1: Write the failing tests.**

  In `test/runtests.jl`, after `include("test_oracles.jl")`, add
  `    include("test_admission.jl")`.

  Create `test/test_admission.jl`:
  ```julia
  # ABOUTME: Tests of the well-posedness guard on inline fixtures: its networks, faces and
  # ABOUTME: classes, their verdicts and messages, each checked against the test oracles.
  using Random

  @testset "the guard's faces follow the spec's list and the oracle's" begin
      for spec in MECHANISM_TEST_SPECS
          cm = EnzymeRates.catalytic(_testhelper_lift(spec.mechanism))
          faces = EnzymeRates._faces(cm)
          @test faces[1] == 0
          @test [Set(EnzymeRates._face_names(cm, f)) for f in faces] == _testhelper_faces(cm)
          # No face zeroes a substrate with a product: the guard lists no mixed corner.
          @test !any(f -> (n = EnzymeRates._face_names(cm, f);
                           any(in(cm.substrates), n) && any(in(cm.products), n)), faces)
      end
      bibi = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) + B ⇌ E(A, B)
              E(A, B) <--> E(P, Q)
              E(P, Q) ⇌ E(Q) + P
              E(Q) ⇌ E + Q
          end
      end
      inhibited = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          modifiers: I
          steps: begin
              E + A ⇌ E(A)
              E(A) + B ⇌ E(A, B)
              E(A, B) <--> E(P, Q)
              E(P, Q) ⇌ E(Q) + P
              E(Q) ⇌ E + Q
              (E + I ⇌ E(I), E(A) + I ⇌ E(A, I))
          end
      end
      terter = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B, C
          products: P, Q, R
          steps: begin
              E + A ⇌ E(A)
              E(A) + B ⇌ E(A, B)
              E(A, B) + C ⇌ E(A, B, C)
              E(A, B, C) <--> E(P, Q, R)
              E(P, Q, R) ⇌ E(Q, R) + P
              E(Q, R) ⇌ E(R) + Q
              E(R) ⇌ E + R
          end
      end
      @test length.(EnzymeRates._faces.((bibi, inhibited, terter))) == (7, 14, 9)
      # P's modifier copy shares P's name and bit: it is zeroed with P, and P adds no
      # modifier face.
      copy_P = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          modifiers: P
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E(P)
              E(P) ⇌ E + P
              E(S) + P::Mod ⇌ E(P::Mod, S)
          end
      end
      @test EnzymeRates._faces(copy_P) == UInt64[0, 1, 2]
      @test EnzymeRates._face_names(copy_P, UInt64(1)) == [:P]
  end

  @testset "a law the guard admits can be 0/0 at a mixed corner, which no fit meets" begin
      # Census mechanism d7a09031d5dd2fc6 (ter-bi, census/base6/mechs.tsv). At B = P = Q = 0
      # the rapid-equilibrium segment {E(B, D), E(B, D, P), E(D, P)} has no form free of B
      # and P, so the law is 0/0 there; the guard's faces never zero one substrate with
      # every product (spec section 8, P4-12a). At A = P = Q = 0, another such corner, the
      # law is 0 (both values measured on phase 3, `probes/mixed_corner_d7a0.log`).
      am = @allosteric_mechanism begin
          substrates: A, B, D
          products: P, Q
          steps: begin
              (E + B ⇌ E(B), E(D, P) + B ⇌ E(B, D, P), E(P) + B ⇌ E(B, P)) :: OnlyA
              (E + P ⇌ E(P), E(B) + P ⇌ E(B, P), E(B, D) + P ⇌ E(B, D, P),
               E(Q) + P ⇌ E(P, Q)) :: EqualAI
              (E + Q <--> E(Q), E(P) + Q <--> E(P, Q)) :: EqualAI
              (E(A, B) + D <--> E(P, Q), E(B) + D <--> E(B, D),
               E(B, P) + D <--> E(B, D, P), E(P) + D <--> E(D, P)) :: OnlyA
              (E(B) + A <--> E(A, B), E(B, D) + A <--> E(P, Q)) :: OnlyA
              E(Q) + A ⇌ E(A, Q) :: EqualAI
          end
      end
      @test EnzymeRates._guard(EnzymeRates.AllostericMechanism(am)) ===
            EnzymeRates._WELL_POSED
      params = _testhelper_random_independent_params(am; rng = Random.MersenneTwister(1))
      @test isnan(rate_equation(am, (A = 1.0, B = 0.0, D = 1.0, P = 0.0, Q = 0.0), params))
      @test rate_equation(am, (A = 0.0, B = 1.0, D = 1.0, P = 0.0, Q = 0.0), params) == 0.0
      # A fit never meets the corner: a zero rate is refused, and the loss scores a NaN
      # prediction as it scores the 0 a well-posed law predicts there.
      _testhelper_table(rate, A, B) =
          (group = [1, 1], Rate = [1.0, rate], A = [1.0, A], B = [1.0, B], D = [1.0, 1.0],
           P = [0.1, 0.0], Q = [0.1, 0.0])
      @test_throws(ErrorException("Zero rate at row 2: log(0) is undefined"),
                   FittingProblem(am, _testhelper_table(0.0, 1.0, 0.0); Keq = 1.0))
      x = [log(params[p]) for p in EnzymeRates.fitted_params(am)]
      at_nan = FittingProblem(am, _testhelper_table(0.5, 1.0, 0.0); Keq = 1.0)
      at_zero = FittingProblem(am, _testhelper_table(0.5, 0.0, 1.0); Keq = 1.0)
      @test EnzymeRates.loss!(x, at_nan) == EnzymeRates.loss!(x, at_zero)
  end

  @testset "the law and ε networks: a node per form and copy, two transitions per edge" begin
      r2 = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)      :: OnlyA
              E(S) <--> E(P)    :: EqualAI
              F + P ⇌ E(P)      :: EqualAI
              F <--> E          :: OnlyA
          end
      end
      bind, iso = _testhelper_edge(r2, "E + S ⇌ ES"), _testhelper_edge(r2, "F <--> E")
      @test EnzymeRates._excludes(r2, bind, 2) && !EnzymeRates._excludes(r2, bind, 1)
      @test EnzymeRates._law_mask(r2, 1) == trues(4)
      # The inactive law deletes both tagged steps, which strands every form but E.
      @test EnzymeRates._law_mask(r2, 2) == falses(4)
      law = EnzymeRates._law_network(r2)
      @test (law.nforms, law.ncopies, length(law.from)) == (4, 2, 10)
      @test law.present == [true, true, true, true, true, false, false, false]
      @test law.edge[9:10] == [0, 0] && law.conformation[9:10] == [0, 0]
      @test (law.from[9], law.to[9]) ==
            (EnzymeRates._node(law, 1, 1), EnzymeRates._node(law, 1, 2))
      @test law.fast[9] && law.fast[10] && law.forward[9] && !law.forward[10]
      @test !any(law.vanishing)
      eps = EnzymeRates._eps_network(r2)
      @test (length(eps.from), eps.present) == (18, trues(8))
      _testhelper_transition(net, k, c, forward) =
          only(t for t in eachindex(net.from)
               if net.edge[t] == k && net.conformation[t] == c && net.forward[t] == forward)
      uptake, release = (_testhelper_transition(eps, bind, 2, d) for d in (true, false))
      blocked = [_testhelper_transition(eps, iso, 2, d) for d in (true, false)]
      # The tagged rapid-equilibrium binding stays fast while its affinity vanishes: its
      # uptake vanishes, its release keeps its rate.
      @test eps.fast[uptake] && eps.fast[release]
      @test eps.vanishing[uptake] && !eps.vanishing[release]
      # The tagged isomerization is blocked: slow, vanishing both ways.
      @test !any(eps.fast[blocked]) && all(eps.vanishing[blocked])
      @test !any(eps.vanishing[t] for t in eachindex(eps.from) if eps.conformation[t] == 1)
      # Names index the table [P, S]: the binding takes up one S, its release gives it off.
      active_uptake = _testhelper_transition(eps, bind, 1, true)
      active_release = _testhelper_transition(eps, bind, 1, false)
      @test eps.change[:, active_uptake] == [0, 1]
      @test (eps.takes[active_uptake], eps.takes[active_release]) == (UInt64(2), UInt64(0))
      mm = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E + P
          end
      end
      law_mm, eps_mm = EnzymeRates._law_network(mm), EnzymeRates._eps_network(mm)
      @test (law_mm.nforms, law_mm.ncopies, length(law_mm.from)) == (2, 1, 4)
      @test all(f -> getfield(law_mm, f) == getfield(eps_mm, f),
                fieldnames(EnzymeRates._Network))
      s_bind = _testhelper_edge(mm, "E + S ⇌ ES")
      @test [law_mm.fast[t] for t in eachindex(law_mm.from) if law_mm.edge[t] == s_bind] ==
            [true, true]
      @test count(law_mm.fast) == 2
  end

  @testset "class 1 rejects a mechanism whose cycles never turn the reaction over" begin
      no_turnover = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E(P)
              E(P) <--> F + P
          end
      end
      @test !EnzymeRates._has_turnover_cycle(no_turnover, trues(3))
      v = EnzymeRates._guard(no_turnover)
      @test v === EnzymeRates._Verdict(EnzymeRates._NO_CYCLE, EnzymeRates._LIMIT_UNTESTED,
                                       0, 0, 0)
      @test EnzymeRates._guard_message(no_turnover, v) ==
            "Mechanism: no cycle of steps turns substrates into products, so the rate is " *
            "zero everywhere."
      # The oracle's steady state carries no flux.
      setprecision(BigFloat, 256) do
          rng = Random.MersenneTwister(1)
          point = _testhelper_consistent_point(no_turnover, rng)
          net = _testhelper_network(no_turnover, point, Dict(:S => big"0.7", :P => big"1.3");
                                    speeds = _testhelper_draw_speeds(no_turnover, rng))
          @test abs(_testhelper_steady_rate(net; Λ = big(10.0)^20).v) < big(10)^-60
      end
  end

  @testset "class 2 rejects a rapid-equilibrium cycle that turns the reaction over" begin
      beside = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) ⇌ E(P)
              E + P ⇌ E(P)
              E(S) <--> F(P)
              F(P) <--> E + P
          end
      end
      v = EnzymeRates._guard(beside)
      @test (v.fault, v.conformation) == (EnzymeRates._RE_TURNOVER, 1)
      @test EnzymeRates._guard_message(beside, v) ==
            "Mechanism: the rapid-equilibrium steps E + P ⇌ E(P), E + S ⇌ E(S) and " *
            "E(S) ⇌ E(P) close a cycle that turns the reaction over, so the mechanism has " *
            "no finite rate. Make one step of the cycle steady-state."
      beside_allo = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)      :: EqualAI
              E(S) ⇌ E(P)       :: EqualAI
              E + P ⇌ E(P)      :: EqualAI
              E(S) <--> F(P)    :: EqualAI
              F(P) <--> E + P   :: EqualAI
          end
      end
      allo_message = EnzymeRates._guard_message(beside_allo, EnzymeRates._guard(beside_allo))
      @test startswith(allo_message,
                       "AllostericMechanism: the rapid-equilibrium steps E + P ⇌ E(P), ")
      flip2 = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E + P ⇌ E(P)
              E + Q ⇌ E(Q)
              E(A) + B ⇌ E(P, Q)
              E(P) + Q ⇌ E(P, Q)
              E(Q) + P <--> E(P, Q)
          end
      end
      # test_rate_eq_derivation.jl "All-RE catalytic cycle raises": its two mechanisms.
      allre = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S1, S2
          products: P1, P2
          steps: begin
              E + S1 ⇌ E(S1)
              E + S2 ⇌ E(S2)
              E(S1) + S2 <--> E(S1, S2)
              E(S2) + S1 ⇌ E(S1, S2)
              E(S1, S2) ⇌ E(P1, P2)
              E(P1, P2) ⇌ E(P1) + P2
              E(P1, P2) <--> E(P2) + P1
              E(P1) ⇌ E + P1
              E(P2) ⇌ E + P2
          end
      end
      sib = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          modifiers: R
          steps: begin
              E + S ⇌ E(S)
              E(S) ⇌ E(P)
              E + P <--> E(P)
              E(R) + S ⇌ E(S, R)
              E(S, R) ⇌ E(P, R)
              E(R) + P ⇌ E(P, R)
              (E + R ⇌ E(R), E(S) + R ⇌ E(S, R), E(P) + R ⇌ E(P, R))
          end
      end
      fix = " Make one step of the cycle steady-state."
      for (i, m) in enumerate((beside, flip2, allre, sib))
          net = EnzymeRates._eps_network(m)
          cycle = EnzymeRates._re_turnover_transitions(net)
          @test !isempty(cycle) && all(t -> net.fast[t], cycle)
          v = EnzymeRates._guard(m)
          @test v.fault === EnzymeRates._RE_TURNOVER
          msg = EnzymeRates._guard_message(m, v)
          @test all(t -> occursin(EnzymeRates._step_text(m, net.edge[t]), msg), cycle)
          # The fix is printed only where any one step of the cycle made steady-state cures
          # every fast turnover cycle; on sib some of its steps leave a second fast cycle.
          @test endswith(msg, "so the mechanism has no finite rate." *
                              (m === sib ? "" : fix))
          # The RE-speed oracle at the interior. On sib's faces with R zero no fast cycle
          # holds enzyme, so the oracle would compare with the package law there, and a
          # mechanism class 2 rejects has no law.
          verdict = setprecision(BigFloat, 256) do
              rng = Random.MersenneTwister(i)
              point = _testhelper_consistent_point(m, rng)
              concs = Dict{Symbol, BigFloat}(x => big(10)^(2 * big(rand(rng)) - 1)
                                             for x in EnzymeRates.metabolites(m))
              _testhelper_respeed(m, point, concs, Set{Symbol}(); rng)
          end
          @test verdict === :diverges
      end
      # The fix: one step of the cycle steady-state. Each variant is well-posed.
      fixed = (
          @_testhelper_unguarded(@enzyme_mechanism begin
              substrates: S
              products: P
              steps: begin
                  E + S ⇌ E(S)
                  E(S) <--> E(P)
                  E + P ⇌ E(P)
                  E(S) <--> F(P)
                  F(P) <--> E + P
              end
          end),
          @_testhelper_unguarded(@enzyme_mechanism begin
              substrates: A, B
              products: P, Q
              steps: begin
                  E + A ⇌ E(A)
                  E + P ⇌ E(P)
                  E + Q ⇌ E(Q)
                  E(A) + B <--> E(P, Q)
                  E(P) + Q ⇌ E(P, Q)
                  E(Q) + P <--> E(P, Q)
              end
          end),
          @_testhelper_unguarded(@enzyme_mechanism begin
              substrates: S1, S2
              products: P1, P2
              steps: begin
                  E + S1 ⇌ E(S1)
                  E + S2 ⇌ E(S2)
                  E(S1) + S2 <--> E(S1, S2)
                  E(S2) + S1 ⇌ E(S1, S2)
                  E(S1, S2) <--> E(P1, P2)
                  E(P1, P2) ⇌ E(P1) + P2
                  E(P1, P2) <--> E(P2) + P1
                  E(P1) ⇌ E + P1
                  E(P2) ⇌ E + P2
              end
          end))
      @test all(m -> EnzymeRates._guard(m) === EnzymeRates._WELL_POSED, fixed)
      # sib turns over through two fast cycles, and some steps of the reported one leave
      # the other when made steady-state, so its message offers no fix: with both
      # chemistry steps steady-state no fast cycle turns over, but at R = 0 that variant
      # splits the segment of free E, which class 3 rejects.
      sib_fixed = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          modifiers: R
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E(P)
              E + P <--> E(P)
              E(R) + S ⇌ E(S, R)
              E(S, R) <--> E(P, R)
              E(R) + P ⇌ E(P, R)
              (E + R ⇌ E(R), E(S) + R ⇌ E(S, R), E(P) + R ⇌ E(P, R))
          end
      end
      @test isempty(EnzymeRates._re_turnover_transitions(EnzymeRates._eps_network(sib_fixed)))
  end

  @testset "the fitted classes pass classes 1 and 2" begin
      michaelis_menten = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E + P
          end
      end
      abortive_complex = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) + B <--> E(A, B)
              E(A, B) <--> E(P, Q)
              E(Q) + P <--> E(P, Q)
              E + Q <--> E(Q)
              E(Q) + A ⇌ E(A, Q)
          end
      end
      product_activated = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + P ⇌ E(P)
              E(P) + S ⇌ E(P, S)
              E(P, S) <--> E(P, P)
              E(P) + P ⇌ E(P, P)
          end
      end
      ss_flux_free_dead_end = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E + P
              E + I <--> E(I)
              E(S) <--> F(S)
          end
      end
      for m in (michaelis_menten, abortive_complex, product_activated, ss_flux_free_dead_end)
          @test EnzymeRates._guard(m) === EnzymeRates._WELL_POSED
      end
  end

  @testset "each step a message names reads back as that step" begin
      pingpong = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) <--> E(P; residual = A - P)
              E(P; residual = A - P) ⇌ E(; residual = A - P) + P
              E(; residual = A - P) + B ⇌ E(B; residual = A - P)
              E(B; residual = A - P) <--> E(Q)
              E(Q) ⇌ E + Q
          end
      end
      theorell_chance = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) + B <--> E(Q) + P
              E(Q) ⇌ E + Q
          end
      end
      merged = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) + B <--> E(P, Q)
              E(P, Q) ⇌ E(Q) + P
              E(Q) ⇌ E + Q
          end
      end
      copies = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          modifiers: I, P
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E(P)
              E(P) ⇌ E + P
              (E + I ⇌ E(I), E(S) + I ⇌ E(I, S))
              E(S) + P::Mod ⇌ E(P::Mod, S)
          end
      end
      mechs = [[_testhelper_lift(s.mechanism) for s in MECHANISM_TEST_SPECS];
               [pingpong, theorell_chance, merged, copies]]
      for mech in mechs
          cm = EnzymeRates.catalytic(mech)
          names = (substrates = cm.substrates, products = cm.products,
                   modifiers = EnzymeRates.modifiers(cm),
                   regulators = mech isa EnzymeRates.AllostericMechanism ?
                                EnzymeRates.allosteric_regulators(mech) : Symbol[])
          raw = EnzymeRates._raw(cm)
          for k in eachindex(cm.edges)
              op, lhs, rhs = EnzymeRates._step_parts(
                  Meta.parse(EnzymeRates._step_text(cm, k)), "@enzyme_mechanism")
              @test EnzymeRates._raw_step(op, lhs, rhs, names, "@enzyme_mechanism") ==
                    raw.steps[k]
          end
      end
      line(m, text) = EnzymeRates._step_text(m, _testhelper_edge(m, text))
      @test line(pingpong, "E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP") ==
            "E(; residual = A - P) + B ⇌ E(B; residual = A - P)"
      @test line(theorell_chance, "EA + B <--> EQ + P") == "E(A) + B <--> E(Q) + P"
      @test line(merged, "EA + B <--> EPQ") == "E(A) + B <--> E(P, Q)"
      @test line(copies, "ES + Pmod ⇌ EPmodS") == "E(S) + P::Mod ⇌ E(P::Mod, S)"
      @test line(copies, "E + I ⇌ EI") == "E + I ⇌ E(I)"
  end
  ```

  Append to `test/test_oracles.jl` (it reads `_testhelper_population_sample`, defined in that
  file):
  ```julia
  @testset "on a plain mechanism class 2 equals the moves' RE-turnover test" begin
      # A plain mechanism's ε-network is its law network, so class 2 must flag exactly the
      # cycles the moves' test flags: on each plain spec and population member as written,
      # with one more kinetic group at rapid equilibrium, and with every step at it.
      plain = [m for m in [[_testhelper_lift(s.mechanism) for s in MECHANISM_TEST_SPECS];
                           _testhelper_population_sample] if m isa EnzymeRates.Mechanism]
      found = Bool[]
      for m in plain
          flags = EnzymeRates._flags(m)
          for re in [[flags]; [flags .| (m.group .== g) for g in 1:EnzymeRates.ngroups(m)];
                     [trues(length(m.edges))]]
              net = EnzymeRates._eps_network(m, re)
              turnover = !isempty(EnzymeRates._re_turnover_transitions(net))
              @test turnover == EnzymeRates._re_turnover_cycle(m, re)
              push!(found, turnover)
          end
      end
      @test any(found) && !all(found)
  end
  ```

- [ ] **Step 2: Run them and watch them fail.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T2_red.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_admission.jl")'
  ```
  Expected: every testset of `test_admission.jl` but "each step a message names reads back
  as that step" reports `Error During Test` with `UndefVarError: `_faces` not defined in
  `EnzymeRates`` (or `_guard`, `_law_network`, `_excludes`, `_has_turnover_cycle`); that one
  passes already, since it reads phase 3's `_step_text`, the line the messages print. Then
  the same command with
  `include("test/test_oracles.jl")` (log `P4_T2_red_oracles.log`): the new testset errors
  with `UndefVarError: `_eps_network` not defined in `EnzymeRates``; the others pass.

- [ ] **Step 3: Implement.**

  (a) `src/types.jl`, after `_inactive_mask`:
  ```julia
  """Whether conformation `c` (1 active, 2 inactive) of `am` excludes edge `k`: an `:OnlyA`
  group's edge in the inactive conformation, an `:OnlyI` group's in the active one. There a
  plain binding loses its affinity and any other edge is blocked (spec section 7)."""
  function _excludes(am::AllostericMechanism, k::Int, c::Int)
      tag = am.cat_allo_states[am.catalytic.group[k]]
      (c == 2 && tag === :OnlyA) || (c == 1 && tag === :OnlyI)
  end

  """
  The edges of conformation `c`'s law (1 active, 2 inactive) in `am`'s catalytic graph:
  those `c` does not exclude (`_excludes`) between forms they join to free E, form 1. The
  law deletes the excluded edges and prunes the forms this strands; admission decides when
  that is the limit of the vanishing constants (spec section 7).
  """
  function _law_mask(am::AllostericMechanism, c::Int)
      m = am.catalytic
      reached = falses(length(m.forms))
      reached[1] = true
      grew = true
      while grew
          grew = false
          for (k, e) in enumerate(m.edges)
              !_excludes(am, k, c) && reached[e.from] != reached[e.to] || continue
              reached[e.from] = reached[e.to] = true
              grew = true
          end
      end
      BitVector([!_excludes(am, k, c) && reached[e.from] for (k, e) in enumerate(m.edges)])
  end
  ```

  (b) `src/names.jl` is unchanged: the messages name a form by its `_form_text` and a
  step by its `_step_text`, the line `show` prints.

  (c) `src/admission.jl`, after `_derivation_input`:
  ```julia
  # ─── The guard: verdicts, networks, faces, classes 1 and 2 ─────────────

  """
  What the guard finds wrong with a mechanism (spec section 8): nothing (`_NO_FAULT`); no
  cycle turns the reaction over (`_NO_CYCLE`, class 1); a cycle of rapid-equilibrium steps
  turns it over (`_RE_TURNOVER`, class 2).
  """
  @enum _Fault::UInt8 begin
      _NO_FAULT
      _NO_CYCLE
      _RE_TURNOVER
  end

  """
  The limit test's verdict on a face (class 3): not run (`_LIMIT_UNTESTED`); the fast-RE
  limit of the steady state equals the law (`_LIMIT_EQUAL`), differs from it
  (`_LIMIT_DIFFER`) or moves with the relative speeds of the rapid-equilibrium steps
  (`_LIMIT_SPEED`); the law has no value there (`_LIMIT_LAW_NAN`).
  """
  @enum _LimitVerdict::UInt8 begin
      _LIMIT_UNTESTED
      _LIMIT_EQUAL
      _LIMIT_DIFFER
      _LIMIT_SPEED
      _LIMIT_LAW_NAN
  end

  """
  The guard's verdict on a mechanism, an isbits value, so a loop over verdicts keeps its
  types concrete: its `fault`, the limit test's verdict where it decided (`limit`), the face
  (an index into `_faces`; 0 when the fault belongs to no face), the conformation it lies in
  (1 active, 2 inactive, 0 none) and the edge its message leads with (0 none).
  """
  struct _Verdict
      fault::_Fault
      limit::_LimitVerdict
      face::Int
      conformation::Int
      edge::Int
  end

  """The verdict on a mechanism that passes every check."""
  const _WELL_POSED = _Verdict(_NO_FAULT, _LIMIT_UNTESTED, 0, 0, 0)

  """
  A network the guard reads (research section 1.1): `ncopies` copies of a mechanism's
  `nforms` forms, node `(c - 1) * nforms + f` for form `f` in copy `c` (`_node`), `present`
  marking the nodes it holds. Transitions come in pairs: `2j - 1` runs an edge (`edge`; 0
  for the flip at free E) in its stored direction in copy `conformation` (0 for the flip),
  and `2j` runs it in reverse; `forward` marks the first of each pair. `fast` marks a
  rapid-equilibrium transition and `vanishing` one whose rate goes to zero in the
  ε-network. `change[n, t]` counts the copies of metabolite name `n` (`_name_index`) that
  transition `t` takes up less those it gives off; bit `n - 1` of `takes[t]` is set when it
  takes up name `n`.
  """
  struct _Network
      nforms::Int
      ncopies::Int
      present::BitVector
      from::Vector{Int}
      to::Vector{Int}
      edge::Vector{Int}
      conformation::Vector{Int}
      forward::BitVector
      fast::BitVector
      vanishing::BitVector
      change::Matrix{Int}
      takes::Vector{UInt64}
  end

  """The node of form `f` in copy `c` of `net`."""
  _node(net::_Network, f::Int, c::Int) = (c - 1) * net.nforms + f

  """
  The network over the `present` nodes holding both directions of edge `k` of `m` in copy
  `c` wherever `use[k, c]` holds, fast where `fast[k, c]` does, the stored direction
  vanishing where `vanish_forward[k, c]` does and the reverse where `vanish_reverse[k, c]`
  does; with two copies, the fast flip between the copies of free E, which takes up and
  gives off nothing.
  """
  function _network(m::Mechanism, present::BitVector, use::BitMatrix, fast::BitMatrix,
                    vanish_forward::BitMatrix, vanish_reverse::BitMatrix)
      nf, ncopies = length(m.forms), size(use, 2)
      name_of = _name_index(m)
      nt = 2 * count(use) + (ncopies == 2 ? 2 : 0)
      from, to, edge, conformation = zeros(Int, nt), zeros(Int, nt), zeros(Int, nt),
                                     zeros(Int, nt)
      forward, isfast, vanishing = falses(nt), falses(nt), falses(nt)
      change, takes = zeros(Int, maximum(name_of; init = 0), nt), zeros(UInt64, nt)
      t = 0
      for c in 1:ncopies, (k, e) in enumerate(m.edges)
          use[k, c] || continue
          a, b = (c - 1) * nf + e.from, (c - 1) * nf + e.to
          for (u, v, gain, loss, gone) in
                  ((a, b, e.consumed, e.released, vanish_forward[k, c]),
                   (b, a, e.released, e.consumed, vanish_reverse[k, c]))
              t += 1
              from[t], to[t], edge[t], conformation[t] = u, v, k, c
              forward[t], isfast[t], vanishing[t] = u == a, fast[k, c], gone
              for i in entries(gain)
                  change[name_of[i], t] += 1
                  takes[t] |= UInt64(1) << (name_of[i] - 1)
              end
              for i in entries(loss)
                  change[name_of[i], t] -= 1
              end
          end
      end
      if ncopies == 2
          from[t + 1], to[t + 1], forward[t + 1], isfast[t + 1] = 1, nf + 1, true, true
          from[t + 2], to[t + 2], isfast[t + 2] = nf + 1, 1, true
      end
      _Network(nf, ncopies, present, from, to, edge, conformation, forward, isfast,
               vanishing, change, takes)
  end

  """
  The law network L (research section 1.1) of `m`, its rapid-equilibrium flags overridden by
  `re`: a plain mechanism's forms and edges as they are; an allosteric mechanism's two
  copies, copy `c` holding the edges of its law (`_law_mask`) and the forms they touch, free
  E always, joined by the flip at free E.
  """
  function _law_network(m::Mechanism, re::BitVector = _flags(m))
      n = length(m.edges)
      fast = falses(n, 1)
      fast[:, 1] .= re
      _network(m, trues(length(m.forms)), trues(n, 1), fast, falses(n, 1), falses(n, 1))
  end
  function _law_network(am::AllostericMechanism, re::BitVector = _flags(am.catalytic))
      m = am.catalytic
      n, nf = length(m.edges), length(m.forms)
      use, fast, present = falses(n, 2), falses(n, 2), falses(2 * nf)
      for c in 1:2
          mask = _law_mask(am, c)
          present[(c - 1) * nf + 1] = true
          for (k, e) in enumerate(m.edges)
              mask[k] || continue
              use[k, c], fast[k, c] = true, re[k]
              present[(c - 1) * nf + e.from] = present[(c - 1) * nf + e.to] = true
          end
      end
      _network(m, present, use, fast, falses(n, 2), falses(n, 2))
  end

  """
  The ε-network N (research section 1.1) of `m`, its rapid-equilibrium flags overridden by
  `re`: every form and edge in each copy. In the copy that excludes an edge (`_excludes`), a
  rapid-equilibrium plain binding stays fast while its uptake vanishes and its release keeps
  its rate; any other excluded edge is slow and vanishes both ways. A plain mechanism's
  ε-network is its law network.
  """
  _eps_network(m::Mechanism, re::BitVector = _flags(m)) = _law_network(m, re)
  function _eps_network(am::AllostericMechanism, re::BitVector = _flags(am.catalytic))
      m = am.catalytic
      n = length(m.edges)
      fast, vanish_forward, vanish_reverse = falses(n, 2), falses(n, 2), falses(n, 2)
      for c in 1:2, (k, e) in enumerate(m.edges)
          excluded, plain = _excludes(am, k, c), _is_plain_binding(m, e)
          fast[k, c] = re[k] && (!excluded || plain)
          vanish_forward[k, c] = excluded
          vanish_reverse[k, c] = excluded && !plain
      end
      _network(m, trues(2 * length(m.forms)), trues(n, 2), fast, vanish_forward,
               vanish_reverse)
  end

  """
  The faces the guard checks (spec section 8) of `m`, each a mask of the metabolite names it
  zeroes (bit `n - 1` for name `n` of `_name_index`), in order: the interior (0), every
  product, every substrate, each substrate and then each product alone on a side of two or
  more; then each of these with every modifier-only name too. A reactant's `::Mod` copy
  shares its reactant's name and so is zeroed with it; a regulatory-site ligand never is.
  No face zeroes a substrate together with a product.
  """
  function _faces(m::Mechanism)
      name_of = _name_index(m)
      bit(x::Symbol) = UInt64(1) << (name_of[searchsortedfirst(m.table, Met(x, false))] - 1)
      side(xs::Vector{Symbol}) = reduce(|, (bit(x) for x in xs); init = UInt64(0))
      base = UInt64[0, side(m.products), side(m.substrates)]
      length(m.substrates) >= 2 && append!(base, bit.(m.substrates))
      length(m.products) >= 2 && append!(base, bit.(m.products))
      mods = UInt64(0)
      for (i, x) in enumerate(m.table)
          x.modifier && !_is_reactant(m, x.name) && (mods |= UInt64(1) << (name_of[i] - 1))
      end
      mods == 0 ? base : [base; base .| mods]
  end

  """The metabolite names that face `face` of `m` zeroes (`_faces`), in table order."""
  function _face_names(m::Mechanism, face::UInt64)
      name_of = _name_index(m)
      Symbol[x.name for (i, x) in enumerate(m.table)
             if (i == 1 || name_of[i] != name_of[i - 1]) &&
                (face >> (name_of[i] - 1)) & 1 == 1]
  end

  """The metabolite-name index (`_name_index`) of `m`'s first substrate, whose net uptake
  is the rate."""
  _turnover_name(m::Mechanism) =
      _name_index(m)[searchsortedfirst(m.table, Met(m.substrates[1], false))]

  """
  Whether a cycle of the edges `mask` keeps turns the reaction over: no potential of the
  forms rises along each edge by its uptake of the first substrate (a modifier copy counting
  as it), which a weighted union-find finds (`_has_unbalanced_cycle`). Class 1 rejects a
  mechanism with none; its rate is zero everywhere.
  """
  function _has_turnover_cycle(m::Mechanism, mask::BitVector)
      name_of, x1 = _name_index(m), _turnover_name(m)
      uptake(c::Counts) =
          sum((copies(c, i) for i in eachindex(m.table) if name_of[i] == x1); init = 0)
      n = length(m.forms)
      _has_unbalanced_cycle([e.from for e in m.edges], [e.to for e in m.edges],
                            [uptake(e.consumed) - uptake(e.released) for e in m.edges],
                            mask, zeros(Int, n), zeros(Int, n))
  end

  """The tree transitions from node `x` up to its root in the forest `parent` (each node's
  tree transition, 0 at a root)."""
  function _tree_path(net::_Network, parent::Vector{Int}, x::Int)
      path = Int[]
      while parent[x] != 0
          t = parent[x]
          push!(path, t)
          x = net.from[t] == x ? net.to[t] : net.from[t]
      end
      path
  end

  """The transitions of the cycle that the non-tree transition `t` between nodes `u` and
  `v` closes in the forest `parent`, in cycle order from `u`."""
  function _tree_cycle(net::_Network, parent::Vector{Int}, u::Int, v::Int, t::Int)
      up, down = _tree_path(net, parent, u), _tree_path(net, parent, v)
      while !isempty(up) && !isempty(down) && up[end] == down[end]
          pop!(up)
          pop!(down)
      end
      [up; reverse(down); t]
  end

  """
  The transitions of a cycle of fast transitions of `net` that takes up some metabolite name
  net (class 2: the rate is infinite), in cycle order; empty when none does. Every name's
  offsets along a breadth-first forest of the fast transitions find the first transition
  they disagree on, which closes the cycle with its tree path. Class 2 runs once, at the
  interior, since a fast cycle that runs on a face runs at the interior too.
  """
  function _re_turnover_transitions(net::_Network)
      nn = net.nforms * net.ncopies
      adjacent = [Int[] for _ in 1:nn]
      for t in eachindex(net.from)
          net.fast[t] && net.forward[t] || continue
          push!(adjacent[net.from[t]], t)
          push!(adjacent[net.to[t]], t)
      end
      offset = zeros(Int, size(net.change, 1), nn)
      parent = zeros(Int, nn)
      seen = falses(nn)
      for root in 1:nn
          net.present[root] && !seen[root] || continue
          seen[root] = true
          queue = [root]
          head = 0
          while head < length(queue)
              head += 1
              u = queue[head]
              for t in adjacent[u]
                  t == parent[u] && continue
                  v, dir = net.from[t] == u ? (net.to[t], 1) : (net.from[t], -1)
                  if !seen[v]
                      seen[v] = true
                      parent[v] = t
                      push!(queue, v)
                      for n in axes(offset, 1)
                          offset[n, v] = offset[n, u] + dir * net.change[n, t]
                      end
                  else
                      for n in axes(offset, 1)
                          offset[n, v] == offset[n, u] + dir * net.change[n, t] ||
                              return _tree_cycle(net, parent, u, v, t)
                      end
                  end
              end
          end
      end
      Int[]
  end

  """
      _guard(m) -> _Verdict

  The well-posedness guard (spec section 8) on `m`, plain or allosteric: class 1 (a cycle
  turns the reaction over, `_has_turnover_cycle`), then class 2 (no cycle of fast
  transitions does, `_re_turnover_transitions`, on the ε-network, where an excluded
  rapid-equilibrium binding stays fast while its affinity vanishes). `_WELL_POSED` when
  every check passes. It reads the reaction graph through masks and builds nothing.
  """
  function _guard(m::Mechanism)
      _has_turnover_cycle(m, trues(length(m.edges))) ||
          return _Verdict(_NO_CYCLE, _LIMIT_UNTESTED, 0, 0, 0)
      _re_turnover_verdict(_eps_network(m))
  end
  function _guard(am::AllostericMechanism)
      m = am.catalytic
      _has_turnover_cycle(m, trues(length(m.edges))) ||
          return _Verdict(_NO_CYCLE, _LIMIT_UNTESTED, 0, 0, 0)
      _re_turnover_verdict(_eps_network(am))
  end

  """Class 2's verdict on the ε-network `net`: `_RE_TURNOVER` in the conformation of the
  first fast cycle that turns the reaction over, led by its first edge; else
  `_WELL_POSED`."""
  function _re_turnover_verdict(net::_Network)
      cycle = _re_turnover_transitions(net)
      isempty(cycle) && return _WELL_POSED
      _Verdict(_RE_TURNOVER, _LIMIT_UNTESTED, 0, net.conformation[cycle[1]],
               net.edge[cycle[1]])
  end

  # ─── Guard messages ────────────────────────────────────────────────────

  """
      _guard_message(mech, v) -> String

  The error for the guard's verdict `v` on `mech` (research section 1.11):
  `<Mechanism|AllostericMechanism>: <cause, naming each step by its line>. <fix>.` A fix is
  printed only where a test shows that its variant passes the guard.
  """
  function _guard_message(mech::_AnyMech, v::_Verdict)
      who = mech isa AllostericMechanism ? "AllostericMechanism" : "Mechanism"
      v.fault === _NO_CYCLE &&
          return "$who: no cycle of steps turns substrates into products, so the rate is " *
                 "zero everywhere."
      v.fault === _RE_TURNOVER && return _re_turnover_message(who, mech)
      error("_guard_message: verdict $v holds no fault")
  end

  """Whether making any one step of `cycle`, transitions of `mech`'s ε-network `net`,
  steady-state with its kinetic group leaves no fast cycle that turns the reaction over: then
  class 2's fix holds whichever step of the cycle is picked (P4-22)."""
  function _any_step_cures(mech::_AnyMech, net::_Network, cycle::Vector{Int})
      m = catalytic(mech)
      members = _group_members(m)
      for t in cycle
          re = _flags(m)
          re[members[m.group[net.edge[t]]]] .= false
          isempty(_re_turnover_transitions(_eps_network(mech, re))) || return false
      end
      true
  end

  """Class 2's message: the steps of the first fast cycle that turns the reaction over,
  each tagged binding among them marked with the conformation where its affinity
  vanishes, and the fix where any one of those steps steady-state cures it
  (`_any_step_cures`)."""
  function _re_turnover_message(who::String, mech::_AnyMech)
      m = catalytic(mech)
      net = _eps_network(mech)
      cycle = _re_turnover_transitions(net)
      c = net.conformation[cycle[1]]
      state = c == 2 ? "inactive" : "active"
      lines = String[]
      for t in cycle
          k = net.edge[t]
          line = _step_text(m, k)
          mech isa AllostericMechanism && _excludes(mech, k, c) &&
              (line *= " (whose affinity vanishes in the $state conformation)")
          push!(lines, line)
      end
      place = mech isa AllostericMechanism && c == 2 ? " in the inactive conformation" : ""
      msg = "$who: the rapid-equilibrium steps $(join(lines, ", ", " and ")) close a " *
            "cycle$place that turns the reaction over, so the mechanism has no finite rate."
      _any_step_cures(mech, net, cycle) ? msg * " Make one step of the cycle steady-state." :
                                          msg
  end
  ```

- [ ] **Step 4: Run them and watch them pass.** The Step 2 commands, logs
  `P4_T2_green.log` and `P4_T2_green_oracles.log`. Expected: `exit 0` for both; the seven
  testsets of `test_admission.jl` and the new testset of `test_oracles.jl` pass; no warning.
  The refresh pass (P4-32) ran Task 2's testsets at 08c4043 with this plan's code
  (`$S/phase45/refresh/map-head-t0-3_task23.py`, logs `map-head-t0-3_7.log` and
  `map-head-t0-3_8.log`), and each passed, Task 2's own pins among them (Open item 1): the
  RE-speed oracle's `:diverges` on `beside`, `flip2`, `allre` and `sib`, R2's network
  counts, edge order and stored directions, and the class-2 messages and fixes; this run
  repeats that check on the branch. A pin that fails is checked against its oracle before
  any change, never re-pinned to the code.

- [ ] **Step 5: Add the ROADMAP known issue (P4-12a).** In `ROADMAP.md`, under "## Known
  issues", after the entry "**A rapid-equilibrium isomerization of a covalent intermediate
  can be unidentifiable.**", add:
  ```markdown
  - **A law the guard admits can be 0/0 where one substrate and every product are zero.**
    The guard checks the faces of the rewrite spec's section 8 (every product zero, every
    substrate zero, each single reactant zero, each also with every modifier zero), never a
    corner that zeroes one substrate with every product. Of 65,423 census mechanisms, 357
    pass the guard yet have a law that is 0/0 at such a corner: 293 plain ter-bi, 38 plain
    ping-pong and 26 allosteric ter-bi mechanisms; there 264 of these laws have no limit and
    93 a removable 0/0, census mechanism `d7a09031d5dd2fc6` at B = P = Q = 0 among them.
    Another 12 allosteric ping-pongs that pass the guard have a non-unique steady state at
    such a corner. `rate_equation` returns NaN there, where a well-posed law returns 0. Fits
    never meet it: `FittingProblem` refuses a zero rate, and the loss scores a NaN prediction
    as a sign mismatch, the penalty a 0 prediction gets (`test/test_admission.jl`). Adding
    these corners to the faces would drop at least these 369 mechanisms, 93 of whose laws
    have a limit there.
  ```

- [ ] **Step 6: Measure the cost of classes 1-2.** Write
  `$S/phase45/probes/P4_T2_cost.jl` (outside the repository):
  ```julia
  # ABOUTME: Times the guard's classes 1 and 2 per mechanism over the shared test specs, the
  # ABOUTME: cost Task 2 of the phase-4 plan reports (prototype 5.7 µs).
  using Test, EnzymeRates, LinearAlgebra, Random
  # `include` resolves a relative path against this file's folder; the run's working
  # directory is the worktree.
  include(joinpath(pwd(), "test", "mechanism_definitions_for_test_enzyme_derivation.jl"))
  mechs = [_testhelper_lift(s.mechanism) for s in MECHANISM_TEST_SPECS]
  foreach(EnzymeRates._guard, mechs)
  t = @elapsed for _ in 1:200, m in mechs
      EnzymeRates._guard(m)
  end
  println("classes 1-2: ", round(1e6 * t / (200 * length(mechs)); digits = 2),
          " µs per mechanism over ", length(mechs), " specs")
  ```
  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T2_cost.log 900 --project=. -e 'using TestEnv; TestEnv.activate(); include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/probes/P4_T2_cost.jl")'
  ```
  Expected: `exit 0` and one line `classes 1-2: <x> µs per mechanism over 42 specs`. Record
  `<x>` for the commit body.

- [ ] **Step 7: Run the full suite.** The Task 1 Step 14 command, log `P4_T2_full.log`.
  Expected: `exit 0`, the previous count plus this task's tests, Aqua and JET pass (JET reads
  the new `src` code: no possible error reported), no warning, maxrss within 5% of M0.

- [ ] **Step 8: Commit.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && git add src/admission.jl src/types.jl ROADMAP.md test/test_admission.jl test/runtests.jl test/test_oracles.jl && git commit -F - <<'EOF'
  Check classes 1 and 2 on the guard's law and ε networks

  The guard's verdict, its two networks (copies per conformation, the flip at free E), the
  spec's faces, and classes 1 and 2: no turnover cycle, and no fast cycle that turns the
  reaction over on the ε-network, where a tagged rapid-equilibrium binding stays fast while
  its affinity vanishes. Messages name each step by the line the DSL reads back as it.
  Nothing calls the guard yet. Classes 1-2 cost <x> µs per mechanism over the 42 specs.

  Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
  EOF
  ```
  with `<x>` from Step 6.

---

### Task 3: Class 3 — procedure A, the limit test and its diagnosis, not wired

Procedure A on the law network (RES §1.4) on every face: a second closed class, a bottomless
segment, Case I and A2 from the structure, else the limit test mod 2^61 − 1 (P4-11, P4-11a)
with the stage it stopped at as the fault. Code from RES §1.4 and `rules/lib.jl`
(`face_structure`, `faceA`, `residual_test` and its helpers) with `rules/modp.jl`, not
`graphcomb/alg.jl` (about 65× slower). Where they differ, RES §1.11 decides: a trap entered
by a slow transition names the rapid-equilibrium step that joins it to the rest of its
segment, where `rules/lib.jl`'s `trap_feeders` names the slow step. The residue stream is a
splitmix64 generator in `src`, since `Random` is a test-only dependency.

**Files:**
- Modify: `src/admission.jl` (`_Fault` gains six values; `_guard` runs the faces; new sections
  "Class 3: procedure A on every face" and "Class 3: the limit test mod 2^61 − 1"; the class-3
  messages)
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (`_testhelper_face_verdicts`,
  `_testhelper_face_faults`, `_testhelper_agrees`, `_testhelper_names_its_step`, at the end
  of "Phase-4 oracles")
- Test: `test/test_admission.jl`

**Interfaces:**
- Consumes: Task 2's `_Network`, `_law_network`, `_eps_network`, `_faces`, `_face_names`,
  `_Verdict`, `_WELL_POSED`, `_guard`, `_guard_message`, `_turnover_name`,
  `_re_turnover_transitions`, `_re_turnover_verdict`; Task 1's
  `_testhelper_respeed_verdicts`, `_testhelper_faces`, `@_testhelper_unguarded`; phase 3's
  `_step_text`, `_form_text` (names.jl), `hash(m)` (content hash, the same in every
  process), `_group_members`, `ngroups`,
  `_name_index`, `_integer_nullspace` (thermo.jl), `catalytic`.
- Produces:
  ```julia
  # _Fault gains, in this order: _NONUNIQUE _BOTTOMLESS _TRAP _SPLIT _LAW_MULTI _MULTILEVEL
  _closed_classes(net::_Network, keep::BitVector)::Vector{Vector{Int}}
  _face_structure(net::_Network, face::UInt64)
      # NamedTuple: alive::BitVector, closed::Vector{Vector{Int}}, segment::Vector{Int},
      # segments::Vector{Vector{Int}}, lowest::BitVector, bottomless::Int,
      # fast_classes::Vector{Vector{Int}}, class_of::Vector{Int}, trap::BitVector,
      # aggregated::Vector{Vector{Int}}, law_classes::Vector{Vector{Int}}
  _face_limit_check(mech::Union{Mechanism, AllostericMechanism}, net::_Network, s,
                    face::UInt64, seed::UInt64)::_LimitVerdict
  _check_face(mech::Union{Mechanism, AllostericMechanism}, net::_Network,
              faces::Vector{UInt64}, i::Int)::_Verdict      # seed hash(i, hash(mech))
  _faces_pass(mech::Union{Mechanism, AllostericMechanism}, re::BitVector)::Bool
  _faces_verdict(mech::Union{Mechanism, AllostericMechanism}, net::_Network)::_Verdict
  _bottom_sccs(next::Vector{Vector{Int}}, nodes::BitVector)::Vector{Vector{Int}}
  _copy_of(net::_Network, u::Int)::Int   # the copy (conformation) of node u
  _form_of(net::_Network, u::Int)::Int   # the form of node u
  _testhelper_face_verdicts(m)::Vector{EnzymeRates._Verdict}     # fixtures file
  _testhelper_face_faults(m)::Vector{EnzymeRates._Fault}         # fixtures file
  _testhelper_agrees(fault::EnzymeRates._Fault, verdict::Symbol)::Bool
  _testhelper_names_its_step(m, v::EnzymeRates._Verdict)::Bool  # fixtures file
  ```
  `_guard` runs `_check_face` on each face after class 2 and returns the first fault.

**Gate:** RES §3.2's plain rows except the two trap specs (Task 7), each with its verdict,
named step and oracle verdict: the fitted classes, `R2_active_conformation`, C3d, C3e, C3f and
the lumpable splits accept every face (`:agree`), `lumpable_split_172` at products = 0 and
Q = 0 and `deadend_split` at P = 0 through the limit test (`_LIMIT_EQUAL`); `ldh_like_trap`
products = 0 `_TRAP` naming `E(P) + Q ⇌ E(P, Q)`, Q = 0 `_SPLIT`; `trap_flip_child9` products
= 0 `_TRAP` naming `E(Q) + P ⇌ E(P, Q)`, P = 0 `_SPLIT`; `split_flip_child1` products = 0 and
Q = 0 `_SPLIT`; `bottomless_entered` and `unfed_bottomless` products = 0 `_BOTTOMLESS`;
`ss_dead_end_at_face` products = 0 `_LIMIT_LAW_NAN`; `fuzz5`, `deadend_closed_unweighted`,
`deadend_closed_weighted`, `B_detour` `_NONUNIQUE` on the products-zero faces;
`unlumped_split_172` products = 0 and Q = 0 `_LIMIT_SPEED`; `lumpable_split_732` Q = 0
`_LIMIT_DIFFER`; `unfed_trap` Q = 0 `_LIMIT_SPEED`, products = 0 accepted by A2;
`mixed_inhibitor_ss_bind` rejected on its three modifier faces; `C3b_multilevel` products = 0
and `C3c_trap_two_transients` products = 0 and Q = 0 `_LIMIT_SPEED`. The guard agrees with
the RE-speed oracle (`_testhelper_agrees`) on every face of every fixture, except
`unfed_bottomless` at products = 0, rejected on purpose (P4-10). None of these fixtures has a
second free form outside free E's rapid-equilibrium segment in an allosteric mechanism, so
every row reads the package law (`law = :package`) and none needs Task 9. Adds at most 3 s;
cost per face reported (prototype 6.6 µs, target ≤ 3 µs, RES §4 risk 2).

- [ ] **Step 1: Write the failing tests.** Add to the "Phase-4 oracles" section of
  `test/mechanism_definitions_for_test_enzyme_derivation.jl`:
  ```julia
  """The guard's class-3 verdict (`_check_face`) on each face of `m`, in `_faces` order."""
  function _testhelper_face_verdicts(m)
      mech = _testhelper_lift(m)
      net = EnzymeRates._law_network(mech)
      faces = EnzymeRates._faces(EnzymeRates.catalytic(mech))
      [EnzymeRates._check_face(mech, net, faces, i) for i in eachindex(faces)]
  end

  """The guard's class-3 fault on each face of `m` (`_testhelper_face_verdicts`)."""
  _testhelper_face_faults(m) = [v.fault for v in _testhelper_face_verdicts(m)]

  """Whether the guard's fault on a face agrees with the RE-speed oracle's verdict there: an
  accepted face exactly when the oracle agrees, a rejected one exactly when it reads a speed
  pattern, a second closed class, divergence, no law value or a mismatch."""
  _testhelper_agrees(fault::EnzymeRates._Fault, verdict::Symbol) =
      fault === EnzymeRates._NO_FAULT ? verdict === :agree :
          verdict in (:pattern, :nonunique, :diverges, :law_nan, :mismatch)

  """Whether the guard's message for the class-3 verdict `v` on `m` names the face and, when
  `v` leads with an edge, that edge's step line."""
  function _testhelper_names_its_step(m, v::EnzymeRates._Verdict)
      mech = _testhelper_lift(m)
      msg = EnzymeRates._guard_message(mech, v)
      line = v.edge == 0 ? "" : EnzymeRates._step_text(EnzymeRates.catalytic(mech), v.edge)
      occursin(": at zero ", msg) && occursin(line, msg)
  end
  ```

  Append to `test/test_admission.jl`:
  ```julia
  @testset "class 3 accepts the fitted classes and lumpable splits on every face" begin
      michaelis_menten = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E + P
          end
      end
      abortive_complex = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) + B <--> E(A, B)
              E(A, B) <--> E(P, Q)
              E(Q) + P <--> E(P, Q)
              E + Q <--> E(Q)
              E(Q) + A ⇌ E(A, Q)
          end
      end
      product_activated = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + P ⇌ E(P)
              E(P) + S ⇌ E(P, S)
              E(P, S) <--> E(P, P)
              E(P) + P ⇌ E(P, P)
          end
      end
      ss_flux_free_dead_end = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E + P
              E + I <--> E(I)
              E(S) <--> F(S)
          end
      end
      # R2's active conformation: a second free form, and no trap (spec section 8).
      r2_active = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E(P)
              F + P ⇌ E(P)
              F <--> E
          end
      end
      tied_isoforms = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              (E + S <--> E(S), F + S <--> F(S))
              (E(S) <--> E + P, F(S) <--> F + P)
              (E + P ⇌ G(P), F + P ⇌ G(P))
              E <--> F
          end
      end
      untied_isoforms = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S <--> E(S)
              F + S <--> F(S)
              (E(S) <--> E + P, F(S) <--> F + P)
              (E + P ⇌ G(P), F + P ⇌ G(P))
              E <--> F
          end
      end
      copy_trap = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          modifiers: P
          steps: begin
              E + S <--> E(S)
              E(S) <--> E(P)
              E(P) ⇌ E + P
              E(P) + P::Mod ⇌ E(P, P::Mod)
              E(P, P::Mod) <--> E(P::Mod) + P
              E + P::Mod <--> E(P::Mod)
          end
      end
      lumpable_172 = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              (E + A <--> E(A), E(B) + A <--> E(P, Q))
              (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))
              (E + Q ⇌ E(Q), E(B) + Q ⇌ E(B, Q))
              (E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q))
          end
      end
      lumpable_532 = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          modifiers: I
          steps: begin
              (E(B) + A ⇌ E(P, Q), E(Q) + A ⇌ E(A, Q))
              (E + B <--> E(B), E(A) + B <--> E(P, Q), E(I) + B <--> E(B, I))
              (E + I ⇌ E(I), E(B) + I ⇌ E(B, I), E(Q) + I ⇌ E(I, Q))
              (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q), E(I) + Q ⇌ E(I, Q))
              E(Q) + P <--> E(P, Q)
          end
      end
      deadend_split = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E + P
              E + P ⇌ E(P)
              F + P ⇌ E(P)
              F <--> E
          end
      end
      accepted = (michaelis_menten, abortive_complex, product_activated,
                  ss_flux_free_dead_end, r2_active, tied_isoforms, untied_isoforms, copy_trap,
                  lumpable_172, lumpable_532, deadend_split)
      for (i, m) in enumerate(accepted)
          @test all(==(EnzymeRates._NO_FAULT), _testhelper_face_faults(m))
          @test all(==(:agree),
                    _testhelper_respeed_verdicts(m; rng = Random.MersenneTwister(i)))
          @test EnzymeRates._guard(m) === EnzymeRates._WELL_POSED
      end
      # Where the structure leaves a split, the limit test finds the law exact.
      limits(m) = [v.limit for v in _testhelper_face_verdicts(m)]
      @test limits(lumpable_172)[[2, 7]] ==
            [EnzymeRates._LIMIT_EQUAL, EnzymeRates._LIMIT_EQUAL]
      @test EnzymeRates._LIMIT_EQUAL in limits(lumpable_532)
      @test limits(deadend_split)[2] === EnzymeRates._LIMIT_EQUAL
  end

  @testset "class 3 rejects a trap and names the step that feeds it" begin
      ok, trap, split = EnzymeRates._NO_FAULT, EnzymeRates._TRAP, EnzymeRates._SPLIT
      ldh = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) + B <--> E(A, B)
              E(A, B) <--> E(P, Q)
              E(P) + Q ⇌ E(P, Q)
              E(Q) + P ⇌ E(P, Q)
              E + P <--> E(P)
              E + Q ⇌ E(Q)
          end
      end
      faults = _testhelper_face_faults(ldh)
      @test faults == [ok, trap, ok, ok, ok, ok, split]
      v = EnzymeRates._guard(ldh)
      @test (v.fault, v.limit, v.face, v.conformation, v.edge) ==
            (trap, EnzymeRates._LIMIT_SPEED, 2, 1, _testhelper_edge(ldh, "EP + Q ⇌ EPQ"))
      @test EnzymeRates._guard_message(ldh, v) ==
            "Mechanism: at zero P and Q, the rapid-equilibrium step E(P) + Q ⇌ E(P, Q) " *
            "fills the forms {E(P)}, which no rapid-equilibrium step running there " *
            "empties, so the rate law misweights them. Make that step's kinetic group " *
            "steady-state."
      verdicts = _testhelper_respeed_verdicts(ldh; rng = Random.MersenneTwister(1))
      @test verdicts[[2, 7]] == [:pattern, :pattern]
      @test all(_testhelper_agrees.(faults, verdicts))
      # The fix: the feeding step steady-state.
      ldh_fed = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) + B <--> E(A, B)
              E(A, B) <--> E(P, Q)
              E(P) + Q <--> E(P, Q)
              E(Q) + P ⇌ E(P, Q)
              E + P <--> E(P)
              E + Q ⇌ E(Q)
          end
      end
      @test EnzymeRates._guard(ldh_fed) === EnzymeRates._WELL_POSED
      @test all(==(:agree),
                _testhelper_respeed_verdicts(ldh_fed; rng = Random.MersenneTwister(2)))
      child9 = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A <--> E(A)
              E + P ⇌ E(P)
              E + Q <--> E(Q)
              E(A) + B ⇌ E(A, B)
              E(A) + Q <--> E(A, Q)
              E(P) + B ⇌ E(B, P)
              E(P) + Q ⇌ E(P, Q)
              E(Q) + A <--> E(A, Q)
              E(Q) + P ⇌ E(P, Q)
              E(A, B) <--> E(P, Q)
          end
      end
      faults9 = _testhelper_face_faults(child9)
      @test faults9 == [ok, trap, ok, ok, ok, split, ok]
      v9 = _testhelper_face_verdicts(child9)[2]
      @test v9.edge == _testhelper_edge(child9, "EQ + P ⇌ EPQ")
      @test occursin("the rapid-equilibrium step E(Q) + P ⇌ E(P, Q) fills the forms",
                     EnzymeRates._guard_message(child9, v9))
      verdicts9 = _testhelper_respeed_verdicts(child9; rng = Random.MersenneTwister(3))
      @test verdicts9[[2, 6]] == [:pattern, :pattern]
      @test all(_testhelper_agrees.(faults9, verdicts9))
      child9_fed = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A <--> E(A)
              E + P ⇌ E(P)
              E + Q <--> E(Q)
              E(A) + B ⇌ E(A, B)
              E(A) + Q <--> E(A, Q)
              E(P) + B ⇌ E(B, P)
              E(P) + Q ⇌ E(P, Q)
              E(Q) + A <--> E(A, Q)
              E(Q) + P <--> E(P, Q)
              E(A, B) <--> E(P, Q)
          end
      end
      @test EnzymeRates._guard(child9_fed) === EnzymeRates._WELL_POSED
      @test all(==(:agree),
                _testhelper_respeed_verdicts(child9_fed; rng = Random.MersenneTwister(4)))
  end

  @testset "class 3 rejects a split, with no fix" begin
      ok, split = EnzymeRates._NO_FAULT, EnzymeRates._SPLIT
      split_child1 = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A <--> E(A)
              E + B <--> E(B)
              E + Q ⇌ E(Q)
              E(A) + B <--> E(P, Q)
              E(B) + Q ⇌ E(B, Q)
              E(B) + A <--> E(P, Q)
              E(Q) + B ⇌ E(B, Q)
              E(Q) + P <--> E(P, Q)
          end
      end
      faults = _testhelper_face_faults(split_child1)
      @test faults == [ok, split, ok, ok, ok, ok, split]
      at_Q = _testhelper_face_verdicts(split_child1)[7]
      msg = EnzymeRates._guard_message(split_child1, at_Q)
      @test startswith(msg, "Mechanism: at zero Q, the forms ")
      @test endswith(msg, "so the rate law cannot weight one against the other.")
      @test !occursin("Make", msg)
      verdicts = _testhelper_respeed_verdicts(split_child1; rng = Random.MersenneTwister(5))
      @test verdicts[[2, 7]] == [:pattern, :pattern]
      @test all(_testhelper_agrees.(faults, verdicts))
      # LDH-like at Q = 0: the lowest forms of free E's segment fall into two classes.
      ldh = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) + B <--> E(A, B)
              E(A, B) <--> E(P, Q)
              E(P) + Q ⇌ E(P, Q)
              E(Q) + P ⇌ E(P, Q)
              E + P <--> E(P)
              E + Q ⇌ E(Q)
          end
      end
      @test EnzymeRates._guard_message(ldh, _testhelper_face_verdicts(ldh)[7]) ==
            "Mechanism: at zero Q, the forms {E, E(A)} and {E(P)} of one rapid-equilibrium " *
            "segment are joined only by steps that take up a zeroed metabolite " *
            "(E + Q ⇌ E(Q) and E(P) + Q ⇌ E(P, Q)), so the rate law cannot weight one " *
            "against the other."
  end

  @testset "class 3 rejects a bottomless segment" begin
      ok, bottomless = EnzymeRates._NO_FAULT, EnzymeRates._BOTTOMLESS
      entered = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) + B ⇌ E(A, B)
              E(A, B) <--> E(P, Q)
              E(P) + Q ⇌ E(P, Q)
              E(Q) + P ⇌ E(P, Q)
              E + P <--> E(P)
              E + Q <--> E(Q)
          end
      end
      faults = _testhelper_face_faults(entered)
      @test faults == [ok, bottomless, ok, ok, ok, ok, ok]
      v = EnzymeRates._guard(entered)
      @test (v.fault, v.face) == (bottomless, 2)
      @test EnzymeRates._guard_message(entered, v) ==
            "Mechanism: at zero P and Q, the rapid-equilibrium segment " *
            "{E(P), E(P, Q), E(Q)} has no form free of the zeroed metabolites, so the rate " *
            "is 0/0 there. Make one of the segment's binding steps steady-state."
      verdicts = _testhelper_respeed_verdicts(entered; rng = Random.MersenneTwister(6))
      @test verdicts[2] === :pattern
      @test all(_testhelper_agrees.(faults, verdicts))
      # The fix: one of the segment's binding steps steady-state.
      entered_fixed = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) + B ⇌ E(A, B)
              E(A, B) <--> E(P, Q)
              E(P) + Q ⇌ E(P, Q)
              E(Q) + P <--> E(P, Q)
              E + P <--> E(P)
              E + Q <--> E(Q)
          end
      end
      @test EnzymeRates._guard(entered_fixed) === EnzymeRates._WELL_POSED
      @test all(==(:agree),
                _testhelper_respeed_verdicts(entered_fixed; rng = Random.MersenneTwister(7)))
      # Catalysis never feeds this segment, so the face's steady state has the limit the
      # oracle's own law finds; the package law is 0/0 there all the same, and the guard
      # rejects the face on purpose (P4-10).
      unfed = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) + B <--> E(A, B)
              E(A, B) <--> E + P + Q
              E + P <--> E(P)
              E(P) + Q ⇌ E(P, Q)
              E(Q) + P ⇌ E(P, Q)
              E + Q <--> E(Q)
          end
      end
      @test _testhelper_face_faults(unfed)[2] === bottomless
      @test _testhelper_respeed_verdicts(unfed; rng = Random.MersenneTwister(8),
                                         law = :limit)[2] === :agree
      @test _testhelper_respeed_verdicts(unfed; rng = Random.MersenneTwister(9))[2] ===
            :law_nan
  end

  @testset "class 3 rejects a second closed class on a face" begin
      fuzz5 = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E ⇌ F
              E + I ⇌ F(I)
              (E + P <--> E(P), E(P) + P <--> E(P, P))
              E + S <--> E(P)
              E + S ⇌ F(S)
              E(P, P) + S <--> F(P, P) + P
          end
      end
      # At P = 0 the dead end F(P, S) is closed and unreachable from E.
      unweighted = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E + P
              E + P <--> E(P)
              E(P) + P <--> E(P, P)
              E(P, P) + S <--> F(P, S) + P
          end
      end
      weighted = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E + P
              E + P <--> E(P)
              E(P) + S <--> F(S) + P
          end
      end
      detour = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: X
          steps: begin
              E + S ⇌ E(S)              :: EqualAI
              E(S) <--> E + P           :: EqualAI
              E + P ⇌ E(P)              :: EqualAI
              E(P) + P ⇌ E(P, P)        :: EqualAI
              F(P) + P <--> E(P, P)     :: EqualAI
              E + X ⇌ E(X)              :: OnlyA
              F(P) + X <--> F(P, X)     :: OnlyA
              E(X) + S ⇌ E(S, X)        :: EqualAI
              E(S, X) <--> F(P, X)      :: EqualAI
          end
      end
      cases = ((fuzz5, [2, 5]), (unweighted, [2]), (weighted, [2]), (detour, [2, 5]))
      for (i, (m, zero_products)) in enumerate(cases)
          faults = _testhelper_face_faults(m)
          @test all(==(EnzymeRates._NONUNIQUE), faults[zero_products])
          verdicts = _testhelper_respeed_verdicts(m; rng = Random.MersenneTwister(10 + i))
          @test all(==(:nonunique), verdicts[zero_products])
          @test all(_testhelper_agrees.(faults, verdicts))
          @test all(v -> _testhelper_names_its_step(m, v),
                    _testhelper_face_verdicts(m)[zero_products])
      end
      # The message names the closed set and the steps out of it, with no fix.
      v = _testhelper_face_verdicts(weighted)[2]
      leaving = only(k for (k, e) in enumerate(weighted.edges)
                     if weighted.forms[e.from].conformation !=
                        weighted.forms[e.to].conformation)
      @test v.edge == leaving
      @test EnzymeRates._guard_message(weighted, v) ==
            "Mechanism: at zero P, the forms {F(S)} keep enzyme away from free E: every " *
            "step that leaves them takes up a zeroed metabolite " *
            "($(EnzymeRates._step_text(weighted, leaving))), so the steady state is not " *
            "unique."
  end

  @testset "the limit test decides what structure cannot" begin
      unlumped = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A <--> E(A)
              E(B) + A <--> E(P, Q)
              (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))
              (E + Q ⇌ E(Q), E(B) + Q ⇌ E(B, Q))
              (E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q))
          end
      end
      lumpable_732 = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              (E + A <--> E(A), E(B) + A <--> E(A, B), E(P) + A <--> E(A, P))
              (E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
              (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q))
              (E + Q ⇌ E(Q), E(B) + Q ⇌ E(B, Q), E(P) + Q ⇌ E(P, Q))
              E(A, B) <--> E(P, Q)
          end
      end
      unfed_trap = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) + B <--> E(A, B)
              E(A, B) <--> E(Q) + P
              E + Q ⇌ E(Q)
              E(Q) + P ⇌ E(P, Q)
              E(P) + Q ⇌ E(P, Q)
              E + P <--> E(P)
          end
      end
      dead_end = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P, Q
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E(P, Q)
              E(Q) + P ⇌ E(P, Q)
              E + Q <--> E(Q)
              E(P) + Q <--> E(P, Q)
          end
      end
      multilevel = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E(P)
              E + P ⇌ E(P)
              E(P) <--> F + P
              F + S ⇌ F(S)
              F(S) <--> F(P)
              F + P ⇌ F(P)
              F(P) <--> E + P
          end
      end
      two_transients = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) + B <--> E(A, B)
              E(A, B) <--> E(P, Q)
              E(P, Q) ⇌ G(P, Q)
              E(Q) + P ⇌ G(P, Q)
              E(P) + Q ⇌ G(P, Q)
              E + P <--> E(P)
              E + Q ⇌ E(Q)
          end
      end
      rows = ((unlumped, [2, 7], EnzymeRates._LIMIT_SPEED, :pattern),
              (lumpable_732, [7], EnzymeRates._LIMIT_DIFFER, :mismatch),
              (unfed_trap, [7], EnzymeRates._LIMIT_SPEED, :pattern),
              (dead_end, [2], EnzymeRates._LIMIT_LAW_NAN, :law_nan),
              (multilevel, [2], EnzymeRates._LIMIT_SPEED, :pattern),
              (two_transients, [2, 7], EnzymeRates._LIMIT_SPEED, :pattern))
      for (i, (m, faces, limit, oracle)) in enumerate(rows)
          verdicts = _testhelper_face_verdicts(m)
          @test all(f -> verdicts[f].fault !== EnzymeRates._NO_FAULT &&
                         verdicts[f].limit === limit, faces)
          oracles = _testhelper_respeed_verdicts(m; rng = Random.MersenneTwister(20 + i))
          @test all(==(oracle), oracles[faces])
          @test all(_testhelper_agrees.([v.fault for v in verdicts], oracles))
          @test all(f -> _testhelper_names_its_step(m, verdicts[f]), faces)
      end
      # unfed_trap at zero products: A2 accepts, since no transition feeds the trap there,
      # and the limit test does not run.
      @test _testhelper_face_verdicts(unfed_trap)[2] ===
            EnzymeRates._Verdict(EnzymeRates._NO_FAULT, EnzymeRates._LIMIT_UNTESTED, 2, 0, 0)
  end

  @testset "the modifier faces catch a steady-state inhibitor binding" begin
      # At I = 0 the rapid-equilibrium square through E(I) and E(I, S) no longer joins E
      # and E(S), whose ratio the steady-state binding of S then sets.
      mixed = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S <--> E(S)
              E(S) <--> E + P
              (E + I ⇌ E(I), E(S) + I ⇌ E(I, S))
              E(I) + S ⇌ E(I, S)
          end
      end
      faults = _testhelper_face_faults(mixed)
      @test all(==(EnzymeRates._NO_FAULT), faults[1:3])
      @test all(!=(EnzymeRates._NO_FAULT), faults[4:6])
      verdicts = _testhelper_respeed_verdicts(mixed; rng = Random.MersenneTwister(31))
      @test verdicts[4:6] == [:mismatch, :mismatch, :mismatch]
      @test all(_testhelper_agrees.(faults, verdicts))
      @test all(v -> _testhelper_names_its_step(mixed, v),
                _testhelper_face_verdicts(mixed)[4:6])
  end

  @testset "the limit test draws one verdict for every writing of a mechanism" begin
      written = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              (E + A <--> E(A), E(B) + A <--> E(A, B), E(P) + A <--> E(A, P))
              (E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
              (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q))
              (E + Q ⇌ E(Q), E(B) + Q ⇌ E(B, Q), E(P) + Q ⇌ E(P, Q))
              E(A, B) <--> E(P, Q)
          end
      end
      # The same mechanism with its groups and steps in reverse order, each binding
      # written as its release.
      rewritten = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E(A, B) <--> E(P, Q)
              (E(P, Q) ⇌ E(P) + Q, E(B, Q) ⇌ E(B) + Q, E(Q) ⇌ E + Q)
              (E(P, Q) <--> E(Q) + P, E(A, P) <--> E(A) + P, E(P) <--> E + P)
              (E(B, Q) ⇌ E(Q) + B, E(A, B) ⇌ E(A) + B)
              (E(A, P) <--> E(P) + A, E(A, B) <--> E(B) + A, E(A) <--> E + A)
          end
      end
      @test written == rewritten
      @test _testhelper_face_verdicts(written) == _testhelper_face_verdicts(rewritten)
      @test _testhelper_face_verdicts(written)[7].limit === EnzymeRates._LIMIT_DIFFER
  end

  @testset "the flip move's predicate reads class 3 through the rapid-equilibrium flags" begin
      ldh = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) + B <--> E(A, B)
              E(A, B) <--> E(P, Q)
              E(P) + Q ⇌ E(P, Q)
              E(Q) + P ⇌ E(P, Q)
              E + P <--> E(P)
              E + Q ⇌ E(Q)
          end
      end
      flags = EnzymeRates._flags(ldh)
      @test !EnzymeRates._faces_pass(ldh, flags)
      # The feeding step steady-state through the flags: the graph of the trap test's fix.
      fed = copy(flags)
      fed[_testhelper_edge(ldh, "EP + Q ⇌ EPQ")] = false
      @test EnzymeRates._faces_pass(ldh, fed)
      # Every step at rapid equilibrium: a fast cycle turns the reaction over.
      @test !EnzymeRates._faces_pass(ldh, trues(length(ldh.edges)))
  end
  ```

- [ ] **Step 2: Run them and watch them fail.** The Task 2 Step 2 command for
  `test/test_admission.jl`, log `P4_T3_red.log`. Expected: the nine new testsets report
  `Error During Test` with `UndefVarError: `_check_face` not defined in `EnzymeRates``
  (`_testhelper_face_verdicts` calls it), or with the first fault each names
  (`_TRAP`, `_SPLIT`, `_BOTTOMLESS`) or `_faces_pass` undefined; Task 2's seven testsets
  still pass.

- [ ] **Step 3: Implement procedure A.** In `src/admission.jl`:

  (a) Replace the `_Fault` enum and its docstring with:
  ```julia
  """
  What the guard finds wrong with a mechanism (spec section 8): nothing (`_NO_FAULT`); no
  cycle turns the reaction over (`_NO_CYCLE`, class 1); a cycle of rapid-equilibrium steps
  turns it over (`_RE_TURNOVER`, class 2); on a face (class 3), a second closed class
  (`_NONUNIQUE`), a rapid-equilibrium segment with no lowest form (`_BOTTOMLESS`), or a law
  the limit test rejects where procedure A stopped at a trap fed on the face (`_TRAP`), a
  segment whose lowest forms split (`_SPLIT`), a law segment chain with several closed
  classes (`_LAW_MULTI`) or an aggregated chain with several (`_MULTILEVEL`).
  """
  @enum _Fault::UInt8 begin
      _NO_FAULT
      _NO_CYCLE
      _RE_TURNOVER
      _NONUNIQUE
      _BOTTOMLESS
      _TRAP
      _SPLIT
      _LAW_MULTI
      _MULTILEVEL
  end
  ```

  (b) Replace the docstring and the two methods of `_guard` with:
  ```julia
  """
      _guard(m) -> _Verdict

  The well-posedness guard (spec section 8) on `m`, plain or allosteric: class 1 (a cycle
  turns the reaction over, `_has_turnover_cycle`); class 2 (no cycle of fast transitions
  does, `_re_turnover_transitions`, on the ε-network, where an excluded rapid-equilibrium
  binding stays fast while its affinity vanishes); then class 3 on every face
  (`_check_face`, on the law network), the first fault winning. `_WELL_POSED` when every
  check passes. It reads the reaction graph through masks and builds nothing.
  """
  function _guard(m::Mechanism)
      _has_turnover_cycle(m, trues(length(m.edges))) ||
          return _Verdict(_NO_CYCLE, _LIMIT_UNTESTED, 0, 0, 0)
      net = _law_network(m)
      v = _re_turnover_verdict(net)
      v === _WELL_POSED || return v
      _faces_verdict(m, net)
  end
  function _guard(am::AllostericMechanism)
      m = am.catalytic
      _has_turnover_cycle(m, trues(length(m.edges))) ||
          return _Verdict(_NO_CYCLE, _LIMIT_UNTESTED, 0, 0, 0)
      v = _re_turnover_verdict(_eps_network(am))
      v === _WELL_POSED || return v
      _faces_verdict(am, _law_network(am))
  end
  ```

  (c) After `_re_turnover_verdict`, add:
  ```julia
  # ─── Class 3: procedure A on every face ────────────────────────────────

  """The closed classes (bottom strongly connected components) of the digraph with successor
  lists `next`, over the nodes `nodes` marks, by Tarjan's algorithm, each as its sorted
  nodes; a node with no successor is a class of its own."""
  function _bottom_sccs(next::Vector{Vector{Int}}, nodes::BitVector)
      n = length(next)
      index, low, comp = zeros(Int, n), zeros(Int, n), zeros(Int, n)
      onstack, stack, work = falses(n), Int[], Tuple{Int, Int}[]
      counter = ncomp = 0
      for s in 1:n
          nodes[s] && index[s] == 0 || continue
          counter += 1
          index[s] = low[s] = counter
          push!(stack, s)
          onstack[s] = true
          push!(work, (s, 1))
          while !isempty(work)
              u, i = work[end]
              if i <= length(next[u])
                  work[end] = (u, i + 1)
                  v = next[u][i]
                  if index[v] == 0
                      counter += 1
                      index[v] = low[v] = counter
                      push!(stack, v)
                      onstack[v] = true
                      push!(work, (v, 1))
                  elseif onstack[v]
                      low[u] = min(low[u], index[v])
                  end
              else
                  pop!(work)
                  isempty(work) || (low[work[end][1]] = min(low[work[end][1]], low[u]))
                  if low[u] == index[u]
                      ncomp += 1
                      while true
                          v = pop!(stack)
                          onstack[v] = false
                          comp[v] = ncomp
                          v == u && break
                      end
                  end
              end
          end
      end
      leaves = falses(ncomp)
      for u in 1:n, v in next[u]
          nodes[u] && comp[u] != comp[v] && (leaves[comp[u]] = true)
      end
      classes = [Int[] for _ in 1:ncomp]
      for u in 1:n
          nodes[u] && push!(classes[comp[u]], u)
      end
      [classes[k] for k in 1:ncomp if !leaves[k]]
  end

  """The closed classes of the transitions of `net` that `keep` marks, over its present
  nodes."""
  function _closed_classes(net::_Network, keep::BitVector)
      next = [Int[] for _ in 1:(net.nforms * net.ncopies)]
      for t in eachindex(net.from)
          keep[t] && push!(next[net.from[t]], net.to[t])
      end
      _bottom_sccs(next, net.present)
  end

  """The fast classes (`class_of`) that the alive fast transitions (`fast_next`) reach from
  node `w`, its own included."""
  function _classes_reached(fast_next::Vector{Vector{Int}}, class_of::Vector{Int}, w::Int)
      seen = falses(length(fast_next))
      seen[w] = true
      todo, found = [w], Int[]
      while !isempty(todo)
          x = pop!(todo)
          class_of[x] > 0 && !(class_of[x] in found) && push!(found, class_of[x])
          for y in fast_next[x]
              seen[y] || (seen[y] = true; push!(todo, y))
          end
      end
      found
  end

  """
      _face_structure(net, face) -> NamedTuple

  The structure procedure A (research section 1.4) reads on face `face` of the law network
  `net`: `alive`, the transitions that take up no zeroed name; `closed`, the closed classes
  of the alive ones; the rapid-equilibrium segments (`segment` per node, `segments`
  breadth-first), each node's uptake of every zeroed name from its segment's first node,
  and `lowest`, the nodes at the segment's componentwise least uptake, whose weight stays
  positive on the face; `bottomless`, the first segment with no lowest node (0 when none);
  the fast classes (closed classes of the alive fast transitions, a node with no alive fast
  exit its own), `class_of` each node's (0 for a node that drains elsewhere) and `trap`,
  the classes of no lowest node; `aggregated`, the closed classes of the chain over fast
  classes that an alive slow transition and then alive fast paths make; and `law_classes`,
  the closed classes of the law's segment chain, its alive slow transitions out of lowest
  nodes.
  """
  function _face_structure(net::_Network, face::UInt64)
      nn, nt = net.nforms * net.ncopies, length(net.from)
      alive = BitVector([net.takes[t] & face == 0 for t in 1:nt])
      closed = _closed_classes(net, alive)
      zeroed = [n for n in axes(net.change, 1) if (face >> (n - 1)) & 1 == 1]
      adjacent = [Int[] for _ in 1:nn]
      for t in 1:nt
          net.fast[t] && push!(adjacent[net.from[t]], t)
      end
      segment, segments = zeros(Int, nn), Vector{Int}[]
      uptake = zeros(Int, length(zeroed), nn)
      for root in 1:nn
          net.present[root] && segment[root] == 0 || continue
          members = [root]
          push!(segments, members)
          segment[root] = length(segments)
          head = 0
          while head < length(members)
              head += 1
              u = members[head]
              for t in adjacent[u]
                  v = net.to[t]
                  segment[v] == 0 || continue
                  segment[v] = length(segments)
                  push!(members, v)
                  for (j, n) in enumerate(zeroed)
                      uptake[j, v] = uptake[j, u] + net.change[n, t]
                  end
              end
          end
      end
      lowest, bottomless = falses(nn), 0
      for (g, members) in enumerate(segments)
          least = [minimum(uptake[j, u] for u in members) for j in eachindex(zeroed)]
          for u in members
              lowest[u] = all(j -> uptake[j, u] == least[j], eachindex(zeroed))
          end
          bottomless == 0 && !any(u -> lowest[u], members) && (bottomless = g)
      end
      fast_alive = alive .& net.fast
      fast_classes = _closed_classes(net, fast_alive)
      class_of = zeros(Int, nn)
      for (k, C) in enumerate(fast_classes), u in C
          class_of[u] = k
      end
      all(C -> all(u -> lowest[u] == lowest[C[1]], C), fast_classes) ||
          error("_face_structure: a fast class holds lowest and other forms")
      trap = BitVector([!lowest[C[1]] for C in fast_classes])
      fast_next = [Int[] for _ in 1:nn]
      for t in 1:nt
          fast_alive[t] && push!(fast_next[net.from[t]], net.to[t])
      end
      aggregated_next = [Int[] for _ in fast_classes]
      law_next = [Int[] for _ in segments]
      for t in 1:nt
          alive[t] && !net.fast[t] || continue
          u, v = net.from[t], net.to[t]
          if class_of[u] > 0
              for k in _classes_reached(fast_next, class_of, v)
                  k == class_of[u] || push!(aggregated_next[class_of[u]], k)
              end
          end
          lowest[u] && segment[u] != segment[v] && push!(law_next[segment[u]], segment[v])
      end
      aggregated = _bottom_sccs(aggregated_next, trues(length(fast_classes)))
      law_classes = _bottom_sccs(law_next, trues(length(segments)))
      (; alive, closed, segment, segments, lowest, bottomless, fast_classes, class_of, trap,
       aggregated, law_classes)
  end

  """
  Where procedure A (research section 1.4) settles a face from its structure `s`
  (`_face_structure`): a second closed class (`_NONUNIQUE`); a segment with no lowest form
  (`_BOTTOMLESS`); `_NO_FAULT` by Case I (every segment holds one fast class, the law's
  segment chain one closed class) or A2 (the aggregated chain has one closed class, holding
  no trap, and each segment it meets holds one lowest class, the law's chain one closed
  class); else the stage at which the limit test decides: the law's chain has several
  closed classes (`_LAW_MULTI`), the aggregated chain has (`_MULTILEVEL`), its closed class
  holds a trap (`_TRAP`), or a segment it meets splits (`_SPLIT`).
  """
  function _face_stage(s)
      length(s.closed) > 1 && return _NONUNIQUE
      s.bottomless > 0 && return _BOTTOMLESS
      classes, lowest = zeros(Int, length(s.segments)), zeros(Int, length(s.segments))
      for (k, C) in enumerate(s.fast_classes)
          g = s.segment[C[1]]
          classes[g] += 1
          s.trap[k] || (lowest[g] += 1)
      end
      one_law = length(s.law_classes) == 1
      one_law && all(==(1), classes) && return _NO_FAULT
      if one_law && length(s.aggregated) == 1
          R = s.aggregated[1]
          !any(k -> s.trap[k], R) &&
              all(k -> lowest[s.segment[s.fast_classes[k][1]]] == 1, R) && return _NO_FAULT
      end
      !one_law ? _LAW_MULTI : length(s.aggregated) > 1 ? _MULTILEVEL :
          any(k -> s.trap[k], s.aggregated[1]) ? _TRAP : _SPLIT
  end

  """The copy (1 active, 2 inactive) of node `u` of `net`."""
  _copy_of(net::_Network, u::Int) = (u - 1) ÷ net.nforms + 1
  """The form of node `u` of `net`."""
  _form_of(net::_Network, u::Int) = (u - 1) % net.nforms + 1

  """The closed class of face structure `s` that holds no active free E (node 1), or its
  first when each holds none."""
  _away_class(s) = s.closed[something(findfirst(C -> !insorted(1, C), s.closed), 1)]

  """The edges whose transitions leave the nodes `C` and take up a zeroed name on the face
  of `s`, in transition order, each once."""
  _exits(net::_Network, s, C::Vector{Int}) =
      unique!([net.edge[t] for t in eachindex(net.from)
               if insorted(net.from[t], C) && !insorted(net.to[t], C) && !s.alive[t] &&
                  net.edge[t] != 0])

  """
  The trap class 3 found fed on a face (stage `_TRAP`) and the rapid-equilibrium edge that
  feeds it (research section 1.11), as `(trap, edge)`, the trap an index into
  `s.fast_classes`. For an alive slow transition out of a class of the aggregated chain's
  closed class, the edge is the last fast transition of an alive fast path from its target
  into a trap of that class; when the slow transition lands in the trap, or no transition
  feeds one (the closed class is a trap alone), it is the rapid-equilibrium edge that joins
  the trap to the rest of its segment, whose way out takes up a zeroed name.
  """
  function _trap_feeder(net::_Network, s)
      nn = net.nforms * net.ncopies
      in_R = falses(length(s.fast_classes))
      in_R[s.aggregated[1]] .= true
      trap_of(u) = (q = s.class_of[u]; q > 0 && in_R[q] && s.trap[q] ? q : 0)
      function joining(q)
          for t in eachindex(net.from)
              net.fast[t] && s.class_of[net.from[t]] == q && s.class_of[net.to[t]] != q &&
                  return (q, net.edge[t])
          end
          (q, 0)
      end
      fast_out = [Int[] for _ in 1:nn]
      for t in eachindex(net.from)
          s.alive[t] && net.fast[t] && push!(fast_out[net.from[t]], t)
      end
      for t in eachindex(net.from)
          k = s.class_of[net.from[t]]
          s.alive[t] && !net.fast[t] && k > 0 && in_R[k] || continue
          v = net.to[t]
          landed = trap_of(v)
          landed > 0 && return joining(landed)
          seen = falses(nn)
          seen[v] = true
          todo = [v]
          while !isempty(todo)
              x = pop!(todo)
              for t2 in fast_out[x]
                  y = net.to[t2]
                  fed = trap_of(y)
                  fed > 0 && s.class_of[x] != fed && return (fed, net.edge[t2])
                  seen[y] || (seen[y] = true; push!(todo, y))
              end
          end
      end
      joining(s.aggregated[1][findfirst(k -> s.trap[k], s.aggregated[1])])
  end

  """
  The split class 3 found on a face (stage `_SPLIT`): the first rapid-equilibrium segment
  that the aggregated chain's closed class meets and that holds more than one lowest fast
  class, those classes (`parts`), and the segment's rapid-equilibrium edges whose
  transitions take up a zeroed name (`joins`), the only edges that could join them.
  """
  function _split_parts(net::_Network, s)
      lowest_classes = [Int[] for _ in s.segments]
      for (k, C) in enumerate(s.fast_classes)
          s.trap[k] || push!(lowest_classes[s.segment[C[1]]], k)
      end
      segment_of(k) = s.segment[s.fast_classes[k][1]]
      R = s.aggregated[1]
      first_split = findfirst(k -> length(lowest_classes[segment_of(k)]) != 1, R)
      segment = segment_of(R[something(first_split, 1)])
      joins = unique!([net.edge[t] for t in eachindex(net.from)
                       if net.fast[t] && !s.alive[t] && net.edge[t] != 0 &&
                          s.segment[net.from[t]] == segment])
      (; segment, parts = [s.fast_classes[k] for k in lowest_classes[segment]], joins)
  end

  """The kinetic-group tags the limit test ties: none for a plain mechanism."""
  _tags(::Mechanism) = Symbol[]
  _tags(am::AllostericMechanism) = am.cat_allo_states

  """
      _check_face(mech, net, faces, i) -> _Verdict

  Class 3 (spec section 8, research section 1.4) on face `faces[i]` of `mech`, whose law
  network is `net`: procedure A rejects a second closed class or a bottomless segment and
  accepts a face whose structure fixes the limit (`_face_stage`); otherwise the limit test
  (`_face_limit_check`, seeded by the content hash and the face) decides, and a rejection
  carries the stage it stopped at, the trap's feeding edge and the split's first joining
  edge. A face the limit test accepts reports `_LIMIT_EQUAL`.
  """
  function _check_face(mech::_AnyMech, net::_Network, faces::Vector{UInt64}, i::Int)
      s = _face_structure(net, faces[i])
      stage = _face_stage(s)
      stage === _NO_FAULT && return _Verdict(_NO_FAULT, _LIMIT_UNTESTED, i, 0, 0)
      if stage === _NONUNIQUE
          C = _away_class(s)
          exits = _exits(net, s, C)
          return _Verdict(_NONUNIQUE, _LIMIT_UNTESTED, i, _copy_of(net, C[1]),
                          isempty(exits) ? 0 : exits[1])
      end
      stage === _BOTTOMLESS &&
          return _Verdict(_BOTTOMLESS, _LIMIT_UNTESTED, i,
                          _copy_of(net, s.segments[s.bottomless][1]), 0)
      limit = _face_limit_check(mech, net, s, faces[i], hash(i, hash(mech)))
      limit === _LIMIT_EQUAL && return _Verdict(_NO_FAULT, _LIMIT_EQUAL, i, 0, 0)
      if stage === _TRAP
          trap, edge = _trap_feeder(net, s)
          return _Verdict(_TRAP, limit, i, _copy_of(net, s.fast_classes[trap][1]), edge)
      elseif stage === _SPLIT
          split = _split_parts(net, s)
          return _Verdict(_SPLIT, limit, i, _copy_of(net, s.segments[split.segment][1]),
                          isempty(split.joins) ? 0 : split.joins[1])
      end
      _Verdict(stage, limit, i, 0, 0)
  end

  """The first class-3 fault of `mech` on its law network `net`, face by face in the order
  of `_faces`; `_WELL_POSED` when every face passes."""
  function _faces_verdict(mech::_AnyMech, net::_Network)
      faces = _faces(catalytic(mech))
      for i in eachindex(faces)
          v = _check_face(mech, net, faces, i)
          v.fault === _NO_FAULT || return v
      end
      _WELL_POSED
  end

  """
  Whether `mech`, its rapid-equilibrium flags overridden by `re`, passes class 2 and then
  class 3 on every face: the flip move asks it of a candidate (`_minimal_flips`), which it
  reads through the override, building nothing. Class 2 runs first because a candidate's
  flags may close a fast cycle its parent's did not.
  """
  function _faces_pass(mech::_AnyMech, re::BitVector)
      isempty(_re_turnover_transitions(_eps_network(mech, re))) || return false
      _faces_verdict(mech, _law_network(mech, re)) === _WELL_POSED
  end
  ```

- [ ] **Step 4: Implement the limit test.** After `_faces_pass`, add:
  ```julia
  # ─── Class 3: the limit test mod 2^61 − 1 ───────────────────────────────

  """The prime 2^61 − 1: the limit test computes every rate, weight and limit modulo it."""
  const _P61 = (Int64(1) << 61) - 1
  """A limit that grows without bound, outside the residues 0…p − 1."""
  const _INFINITE = Int64(-1)
  """A law that has no value on a face, outside the residues 0…p − 1."""
  const _UNDEFINED = Int64(-2)

  _mulp(a::Int64, b::Int64) = Int64(mod(widemul(a, b), _P61))
  _addp(a::Int64, b::Int64) = (s = a + b; s >= _P61 ? s - _P61 : s)
  _subp(a::Int64, b::Int64) = (s = a - b; s < 0 ? s + _P61 : s)
  _modp(x::Integer) = Int64(mod(x, _P61))

  """`a^k` mod p; a negative `k` raises the inverse."""
  function _powp(a::Int64, k::Integer)
      k < 0 && return _powp(_invp(a), -k)
      r = Int64(1)
      while k > 0
          isodd(k) && (r = _mulp(r, a))
          a = _mulp(a, a)
          k >>= 1
      end
      r
  end

  """The inverse of `a` mod p."""
  _invp(a::Int64) = (a == 0 && throw(DivideError()); _powp(a, _P61 - 2))

  """A deterministic stream of residues in 2…p − 2 (splitmix64), seeded by a mechanism's
  content hash and a face, so every process draws the same points."""
  mutable struct _ResidueStream
      state::UInt64
  end

  """The next residue of `r`."""
  function _residue!(r::_ResidueStream)
      r.state += 0x9e3779b97f4a7c15
      z = r.state
      z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
      z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
      Int64(2 + (z ⊻ (z >> 31)) % UInt64(_P61 - 3))
  end

  """`A \\ B` mod p for square `A`; an empty matrix when `A` is singular."""
  function _solvep(A::Matrix{Int64}, B::Matrix{Int64})
      n = size(A, 1)
      M = hcat(A, B)
      for j in 1:n
          k = findfirst(i -> M[i, j] != 0, j:n)
          k === nothing && return Matrix{Int64}(undef, 0, 0)
          k += j - 1
          k != j && (M[[j, k], :] = M[[k, j], :])
          pivot = _invp(M[j, j])
          for c in j:size(M, 2)
              M[j, c] = _mulp(M[j, c], pivot)
          end
          for i in 1:n
              (i == j || M[i, j] == 0) && continue
              f = M[i, j]
              for c in j:size(M, 2)
                  M[i, c] = _subp(M[i, c], _mulp(f, M[j, c]))
              end
          end
      end
      M[:, (n + 1):end]
  end

  """The determinant of `A` mod p."""
  function _detp(A::Matrix{Int64})
      n = size(A, 1)
      M = copy(A)
      d = Int64(1)
      for j in 1:n
          k = findfirst(i -> M[i, j] != 0, j:n)
          k === nothing && return Int64(0)
          k += j - 1
          if k != j
              M[[j, k], :] = M[[k, j], :]
              d = _subp(Int64(0), d)
          end
          d = _mulp(d, M[j, j])
          pivot = _invp(M[j, j])
          for i in (j + 1):n
              M[i, j] == 0 && continue
              f = _mulp(M[i, j], pivot)
              for c in j:n
                  M[i, c] = _subp(M[i, c], _mulp(f, M[j, c]))
              end
          end
      end
      d
  end

  """The stationary distribution mod p of the chain with off-diagonal rates `R[i, j]`
  (`i → j`) and one closed class covering every state; empty when singular."""
  function _stationaryp(R::Matrix{Int64})
      n = size(R, 1)
      n == 1 && return Int64[1]
      A = zeros(Int64, n, n)
      for i in 1:n, j in 1:n
          i == j && continue
          A[j, i] = _addp(A[j, i], R[i, j])
          A[i, i] = _subp(A[i, i], R[i, j])
      end
      A[n, :] .= 1
      b = zeros(Int64, n, 1)
      b[n] = 1
      x = _solvep(A, b)
      isempty(x) ? Int64[] : vec(x)
  end

  # Polynomials mod p: coefficient vectors, constant term first.
  _ptrimp!(a::Vector{Int64}) = (while length(a) > 1 && a[end] == 0; pop!(a); end; a)
  _degp(a::Vector{Int64}) =
      (b = _ptrimp!(copy(a)); length(b) == 1 && b[1] == 0 ? -1 : length(b) - 1)
  _paddp(a::Vector{Int64}, b::Vector{Int64}) =
      _ptrimp!([_addp(get(a, i, Int64(0)), get(b, i, Int64(0)))
                for i in 1:max(length(a), length(b))])
  _psubp(a::Vector{Int64}, b::Vector{Int64}) =
      _ptrimp!([_subp(get(a, i, Int64(0)), get(b, i, Int64(0)))
                for i in 1:max(length(a), length(b))])
  function _pmulp(a::Vector{Int64}, b::Vector{Int64})
      c = zeros(Int64, length(a) + length(b) - 1)
      for i in eachindex(a), j in eachindex(b)
          c[i + j - 1] = _addp(c[i + j - 1], _mulp(a[i], b[j]))
      end
      _ptrimp!(c)
  end
  function _pdivremp(a::Vector{Int64}, b::Vector{Int64})
      a, b = _ptrimp!(copy(a)), _ptrimp!(copy(b))
      db = _degp(b)
      db < 0 && throw(DivideError())
      q = zeros(Int64, max(_degp(a) - db + 1, 1))
      lead = _invp(b[end])
      while _degp(a) >= db
          k = _degp(a) - db
          c = _mulp(a[end], lead)
          q[k + 1] = c
          for i in 0:db
              a[i + k + 1] = _subp(a[i + k + 1], _mulp(c, b[i + 1]))
          end
          _ptrimp!(a)
      end
      _ptrimp!(q), a
  end
  _pevalp(a::Vector{Int64}, x::Int64) =
      foldr((c, acc) -> _addp(c, _mulp(acc, x)), a; init = Int64(0))
  """The lowest-order term of `a` as `(order, coefficient)`; `(typemax(Int), 0)` for 0."""
  _plowp(a::Vector{Int64}) =
      (k = findfirst(!=(0), a); k === nothing ? (typemax(Int), Int64(0)) : (k - 1, a[k]))

  """Newton interpolation mod p through the points `(x, y)`."""
  function _interpp(x::Vector{Int64}, y::Vector{Int64})
      n = length(x)
      c = copy(y)
      for j in 2:n, i in n:-1:j
          c[i] = _mulp(_subp(c[i], c[i - 1]), _invp(_subp(x[i], x[i - j + 1])))
      end
      p = Int64[c[n]]
      for i in (n - 1):-1:1
          p = _paddp(_pmulp(p, Int64[_subp(Int64(0), x[i]), 1]), Int64[c[i]])
      end
      p
  end

  """`(N, D)` with degrees at most `dmax` and `N/D` through the points: the extended
  Euclidean algorithm on the interpolant (it needs at least `2 dmax + 1` points)."""
  function _ratrecp(x::Vector{Int64}, y::Vector{Int64}, dmax::Int)
      M = Int64[1]
      for xi in x
          M = _pmulp(M, Int64[_subp(Int64(0), xi), 1])
      end
      r0, r1 = M, _interpp(x, y)
      t0, t1 = Int64[0], Int64[1]
      while _degp(r1) > dmax
          q, r = _pdivremp(r0, r1)
          r0, r1 = r1, r
          t0, t1 = t1, _psubp(t0, _pmulp(q, t1))
      end
      r1, t1
  end

  """The rational function of one variable that `f` evaluates mod p, reconstructed from
  random points with a degree bound doubled until two more points agree; empty vectors past
  degree 1024."""
  function _reconstructp(f, rng::_ResidueStream)
      d = 8
      while d <= 1024
          xs = [_residue!(rng) for _ in 1:(2d + 3)]
          ys = [f(x) for x in xs]
          N, D = _ratrecp(xs[1:(end - 2)], ys[1:(end - 2)], d)
          all(k -> (dv = _pevalp(D, xs[k]); dv != 0 &&
                    _mulp(_pevalp(N, xs[k]), _invp(dv)) == ys[k]),
              (length(xs) - 1):length(xs)) && return N, D
          d *= 2
      end
      Int64[], Int64[]
  end

  """
  A random point mod p of `m`'s model (research section 1.4) for `ncopies` conformations on
  face `face`: log-free energies of every form in each conformation and of every
  metabolite name, an integer null-space basis of the kinetic-group ties (members equal in
  each conformation that keeps them, an `:EqualAI` member equal across) raised to random
  residues; from them each edge's equilibrium constant `K[k, c]` and `L`; a steady-state
  speed per kinetic group and conformation (one for both of an `:EqualAI` group), a
  rapid-equilibrium speed per edge and conformation, the flip's speed; concentrations
  `conc` per name, zero on the face, and `concw`, one on the face, which the weights'
  leading coefficients read.
  """
  function _model_point(m::Mechanism, tags::Vector{Symbol}, ncopies::Int, face::UInt64,
                        rng::_ResidueStream)
      nf, name_of = length(m.forms), _name_index(m)
      nnames = maximum(name_of; init = 0)
      nv = ncopies * nf + nnames
      excluded(k, c) = !isempty(tags) &&
          ((c == 2 && tags[m.group[k]] === :OnlyA) || (c == 1 && tags[m.group[k]] === :OnlyI))
      function drop!(row, k, c, dir)
          e = m.edges[k]
          row[(c - 1) * nf + e.from] += dir
          row[(c - 1) * nf + e.to] -= dir
          foreach(i -> row[ncopies * nf + name_of[i]] += dir, entries(e.consumed))
          foreach(i -> row[ncopies * nf + name_of[i]] -= dir, entries(e.released))
          row
      end
      ties = Vector{Int}[]
      for members in _group_members(m), c in 1:ncopies
          lead = first(members)
          excluded(lead, c) && continue
          for k in members[2:end]
              push!(ties, drop!(drop!(zeros(Int, nv), k, c, 1), lead, c, -1))
          end
          if c == 2 && tags[m.group[lead]] === :EqualAI
              for k in members
                  push!(ties, drop!(drop!(zeros(Int, nv), k, 2, 1), k, 1, -1))
              end
          end
      end
      B = _integer_nullspace(isempty(ties) ? zeros(Int, 0, nv) :
                             permutedims(reduce(hcat, ties)))
      t = [_residue!(rng) for _ in axes(B, 2)]
      x = [reduce(_mulp, (_powp(t[j], B[i, j]) for j in axes(B, 2) if B[i, j] != 0);
                  init = Int64(1)) for i in 1:nv]
      K = ones(Int64, length(m.edges), ncopies)
      for (k, e) in enumerate(m.edges), c in 1:ncopies
          num, den = x[(c - 1) * nf + e.from], x[(c - 1) * nf + e.to]
          for i in entries(e.consumed)
              num = _mulp(num, x[ncopies * nf + name_of[i]])
          end
          for i in entries(e.released)
              den = _mulp(den, x[ncopies * nf + name_of[i]])
          end
          K[k, c] = _mulp(num, _invp(den))
      end
      ss = [_residue!(rng) for _ in 1:ngroups(m), _ in 1:ncopies]
      for g in 1:ngroups(m)
          ncopies == 2 && tags[g] === :EqualAI && (ss[g, 2] = ss[g, 1])
      end
      re = [_residue!(rng) for _ in m.edges, _ in 1:ncopies]
      flip = _residue!(rng)
      L = ncopies == 2 ? _mulp(x[1], _invp(x[nf + 1])) : Int64(1)
      conc = [(face >> (n - 1)) & 1 == 1 ? Int64(0) : _residue!(rng) for n in 1:nnames]
      concw = [(face >> (n - 1)) & 1 == 1 ? Int64(1) : conc[n] for n in 1:nnames]
      (; K, ss, re, flip, L, conc, concw)
  end

  """`pt` with new rapid-equilibrium and flip speeds: another speed pattern."""
  _respeed_point(pt, rng::_ResidueStream) =
      merge(pt, (re = [_residue!(rng) for _ in axes(pt.re, 1), _ in axes(pt.re, 2)],
                 flip = _residue!(rng)))

  """The rate mod p of every transition of `net` at the concentrations `conc` (per name):
  its speed times `K` in the stored direction, its speed alone in reverse, times each
  metabolite it takes up, a fast one times `Λ` unless `Λ` is 0."""
  function _rates_p(m::Mechanism, net::_Network, pt, conc::Vector{Int64}, Λ::Int64)
      name_of = _name_index(m)
      r = zeros(Int64, length(net.from))
      for t in eachindex(net.from)
          k, c = net.edge[t], net.conformation[t]
          if k == 0
              r[t] = net.forward[t] ? _mulp(pt.flip, pt.L) : pt.flip
          else
              speed = net.fast[t] ? pt.re[k, c] : pt.ss[m.group[k], c]
              r[t] = net.forward[t] ? _mulp(speed, pt.K[k, c]) : speed
              e = m.edges[k]
              for i in entries(net.forward[t] ? e.consumed : e.released)
                  r[t] = _mulp(r[t], conc[name_of[i]])
              end
          end
          net.fast[t] && Λ != 0 && (r[t] = _mulp(r[t], Λ))
      end
      r
  end

  """Each node's potential of the turnover name `x1` along the fast transitions of `net`,
  from its segment's first node; class 2 has ruled out a fast cycle that takes it up."""
  function _fast_potential(net::_Network, x1::Int)
      nn = net.nforms * net.ncopies
      h, seen = zeros(Int, nn), falses(nn)
      next = [Int[] for _ in 1:nn]
      for t in eachindex(net.from)
          net.fast[t] && push!(next[net.from[t]], t)
      end
      for root in 1:nn
          net.present[root] && !seen[root] || continue
          seen[root] = true
          todo = [root]
          while !isempty(todo)
              u = pop!(todo)
              for t in next[u]
                  v = net.to[t]
                  seen[v] && continue
                  seen[v] = true
                  h[v] = h[u] + net.change[x1, t]
                  push!(todo, v)
              end
          end
      end
      h
  end

  """ρ(u) mod p: over the slow transitions out of node `u`, their rate `r` times their
  uptake of `x1` plus the fall of its potential `h` along them."""
  function _slow_rho(net::_Network, r::Vector{Int64}, h::Vector{Int}, x1::Int)
      rho = zeros(Int64, net.nforms * net.ncopies)
      for t in eachindex(net.from)
          (net.fast[t] || r[t] == 0) && continue
          u = net.from[t]
          rho[u] = _addp(rho[u],
                         _mulp(r[t], _modp(net.change[x1, t] + h[u] - h[net.to[t]])))
      end
      rho
  end

  """Equilibrium weights mod p within each part of `parts` over the fast transitions inside
  it, from its first node, at rates `r`."""
  function _eq_weights_p(net::_Network, parts::Vector{Vector{Int}}, r::Vector{Int64})
      nn = net.nforms * net.ncopies
      part = zeros(Int, nn)
      for (q, P) in enumerate(parts), u in P
          part[u] = q
      end
      next = [Int[] for _ in 1:nn]
      for t in eachindex(net.from)
          u = net.from[t]
          net.fast[t] && part[u] > 0 && part[u] == part[net.to[t]] && push!(next[u], t)
      end
      w = zeros(Int64, nn)
      for P in parts
          w[P[1]] = 1
          done = falses(nn)
          done[P[1]] = true
          todo = [P[1]]
          while !isempty(todo)
              u = pop!(todo)
              for t in next[u]
                  v = net.to[t]
                  done[v] && continue
                  back = isodd(t) ? t + 1 : t - 1
                  r[t] != 0 && r[back] != 0 ||
                      error("_face_limit_check: a fast transition inside a part runs one way")
                  w[v] = _mulp(w[u], _mulp(r[t], _invp(r[back])))
                  done[v] = true
                  push!(todo, v)
              end
          end
          all(u -> done[u], P) || error("_face_limit_check: a part is not fast-connected")
      end
      w
  end

  """The face chain's limit mod p as every rapid-equilibrium speed grows, by first-order
  aggregation over the fast classes of the aggregated chain's one closed class."""
  function _limit_by_aggregation(m::Mechanism, net::_Network, s, pt, h::Vector{Int},
                                 x1::Int)
      nn = net.nforms * net.ncopies
      r = _rates_p(m, net, pt, pt.conc, Int64(0))
      R = s.aggregated[1]
      place = zeros(Int, length(s.fast_classes))
      place[R] .= eachindex(R)
      nu = _eq_weights_p(net, [s.fast_classes[k] for k in R],
                         _rates_p(m, net, pt, pt.concw, Int64(0)))
      for k in R
          C = s.fast_classes[k]
          inverse = _invp(reduce(_addp, nu[C]; init = Int64(0)))
          nu[C] .= _mulp.(nu[C], inverse)
      end
      transient = [u for u in 1:nn if net.present[u] && s.class_of[u] == 0]
      slot = zeros(Int, nn)
      slot[transient] .= eachindex(transient)
      X = zeros(Int64, length(transient), length(R))
      if !isempty(transient)
          M = zeros(Int64, length(transient), length(transient))
          rhs = zeros(Int64, length(transient), length(R))
          for t in eachindex(net.from)
              u, v = net.from[t], net.to[t]
              net.fast[t] && s.alive[t] && r[t] != 0 && slot[u] > 0 || continue
              M[slot[u], slot[u]] = _addp(M[slot[u], slot[u]], r[t])
              if slot[v] > 0
                  M[slot[u], slot[v]] = _subp(M[slot[u], slot[v]], r[t])
              elseif s.class_of[v] > 0 && place[s.class_of[v]] > 0
                  q = place[s.class_of[v]]
                  rhs[slot[u], q] = _addp(rhs[slot[u], q], r[t])
              end
          end
          X = _solvep(M, rhs)
          isempty(X) && error("_face_limit_check: singular absorption system")
      end
      Q = zeros(Int64, length(R), length(R))
      for t in eachindex(net.from)
          u, v = net.from[t], net.to[t]
          !net.fast[t] && s.alive[t] && r[t] != 0 && s.class_of[u] > 0 &&
              place[s.class_of[u]] > 0 || continue
          a, f = place[s.class_of[u]], _mulp(nu[u], r[t])
          if s.class_of[v] > 0
              b = place[s.class_of[v]]
              b == 0 && error("_face_limit_check: the aggregated class is not closed")
              b != a && (Q[a, b] = _addp(Q[a, b], f))
          else
              for b in eachindex(R)
                  b == a || (Q[a, b] = _addp(Q[a, b], _mulp(f, X[slot[v], b])))
              end
          end
      end
      mu = _stationaryp(Q)
      isempty(mu) && error("_face_limit_check: singular aggregated chain")
      rho = _slow_rho(net, r, h, x1)
      v = Int64(0)
      for (q, k) in enumerate(R), u in s.fast_classes[k]
          v = _addp(v, _mulp(_mulp(mu[q], nu[u]), rho[u]))
      end
      v
  end

  """The face chain's limit mod p as every rapid-equilibrium speed grows, by interpolating
  the rate's numerator and its tree-sum denominator in Λ (degree below the closed class's
  size) on the face chain's closed class; `_INFINITE` when the rate grows with Λ."""
  function _limit_by_interpolation(m::Mechanism, net::_Network, s, pt, h::Vector{Int},
                                   x1::Int, rng::_ResidueStream)
      C = s.closed[1]
      n = length(C)
      pos = zeros(Int, net.nforms * net.ncopies)
      pos[C] .= 1:n
      xs, Ds, Ns = Int64[], Int64[], Int64[]
      for _ in 1:(n + 1)
          Λ = _residue!(rng)
          r = _rates_p(m, net, pt, pt.conc, Λ)
          Q = zeros(Int64, n, n)
          for t in eachindex(net.from)
              a, b = pos[net.from[t]], pos[net.to[t]]
              s.alive[t] && r[t] != 0 && a > 0 && b > 0 && a != b || continue
              Q[a, b] = _addp(Q[a, b], r[t])
              Q[a, a] = _subp(Q[a, a], r[t])
          end
          A = permutedims(Q)
          A[n, :] .= 1
          rhs = zeros(Int64, n, 1)
          rhs[n] = 1
          x = _solvep(A, rhs)
          isempty(x) && error("_face_limit_check: singular at a random Λ")
          rho = _slow_rho(net, r, h, x1)
          v = reduce(_addp, (_mulp(x[pos[u], 1], rho[u]) for u in C); init = Int64(0))
          D = _detp(A)
          push!(xs, Λ)
          push!(Ds, D)
          push!(Ns, _mulp(D, v))
      end
      cD, cN = _interpp(xs[1:n], Ds[1:n]), _interpp(xs[1:n], Ns[1:n])
      _pevalp(cD, xs[end]) == Ds[end] && _pevalp(cN, xs[end]) == Ns[end] ||
          error("_face_limit_check: the degree bound failed")
      dD, dN = _degp(cD), _degp(cN)
      dN > dD && return _INFINITE
      dN < dD && return Int64(0)
      _mulp(cN[dD + 1], _invp(cD[dD + 1]))
  end

  """The law's value mod p on the face, from the law's segment chain with its one closed
  class: each segment's weight on its lowest forms, its alive slow transitions out of
  them."""
  function _law_on_face(m::Mechanism, net::_Network, s, pt, h::Vector{Int}, x1::Int)
      r = _rates_p(m, net, pt, pt.conc, Int64(0))
      parts = [sort(S; by = u -> !s.lowest[u]) for S in s.segments]
      w = _eq_weights_p(net, parts, _rates_p(m, net, pt, pt.concw, Int64(0)))
      share = zeros(Int64, length(w))
      for S in s.segments
          low = [u for u in S if s.lowest[u]]
          inverse = _invp(reduce(_addp, w[low]; init = Int64(0)))
          share[low] .= _mulp.(w[low], inverse)
      end
      Lc = s.law_classes[1]
      at = zeros(Int, length(s.segments))
      at[Lc] .= eachindex(Lc)
      Q = zeros(Int64, length(Lc), length(Lc))
      for t in eachindex(net.from)
          u = net.from[t]
          sa, sb = s.segment[u], s.segment[net.to[t]]
          !net.fast[t] && s.alive[t] && s.lowest[u] && at[sa] > 0 && sa != sb || continue
          at[sb] == 0 && error("_face_limit_check: the law's closed class is not closed")
          Q[at[sa], at[sb]] = _addp(Q[at[sa], at[sb]], _mulp(share[u], r[t]))
      end
      mu = _stationaryp(Q)
      isempty(mu) && error("_face_limit_check: singular law segment chain")
      rho = _slow_rho(net, r, h, x1)
      v = Int64(0)
      for (q, g) in enumerate(Lc), u in s.segments[g]
          s.lowest[u] && (v = _addp(v, _mulp(_mulp(mu[q], share[u]), rho[u])))
      end
      v
  end

  """The law mod p at concentrations `conc` with every name present: segment weights and
  the slow segment chain over every segment."""
  function _law_inside(m::Mechanism, net::_Network, s, pt, conc::Vector{Int64},
                       h::Vector{Int}, x1::Int)
      r = _rates_p(m, net, pt, conc, Int64(0))
      w = _eq_weights_p(net, s.segments, r)
      W = [_invp(reduce(_addp, w[S]; init = Int64(0))) for S in s.segments]
      Q = zeros(Int64, length(s.segments), length(s.segments))
      for t in eachindex(net.from)
          net.fast[t] && continue
          a, b = s.segment[net.from[t]], s.segment[net.to[t]]
          a == b && continue
          Q[a, b] = _addp(Q[a, b], _mulp(_mulp(w[net.from[t]], W[a]), r[t]))
      end
      mu = _stationaryp(Q)
      isempty(mu) && error("_face_limit_check: singular segment chain")
      rho = _slow_rho(net, r, h, x1)
      v = Int64(0)
      for (g, S) in enumerate(s.segments), u in S
          v = _addp(v, _mulp(_mulp(mu[g], _mulp(w[u], W[g])), rho[u]))
      end
      v
  end

  """The law's value mod p on the face when its segment chain has several closed classes:
  its limit as the zeroed concentrations go to zero along two random directions, by
  rational reconstruction in their scale; `_UNDEFINED` when a direction meets a pole or the
  two differ."""
  function _law_by_reconstruction(m::Mechanism, net::_Network, s, pt, h::Vector{Int},
                                  x1::Int, face::UInt64, rng::_ResidueStream)
      found = Int64[]
      for _ in 1:2
          direction = [_residue!(rng) for _ in pt.conc]
          law_at(δ) = _law_inside(m, net, s, pt,
                                  [(face >> (n - 1)) & 1 == 1 ? _mulp(direction[n], δ) :
                                   pt.conc[n] for n in eachindex(pt.conc)], h, x1)
          num, den = _reconstructp(law_at, rng)
          isempty(den) && error("_face_limit_check: no reconstruction under the degree cap")
          (a, lead_num), (b, lead_den) = _plowp(num), _plowp(den)
          push!(found, a > b ? Int64(0) : a < b ? _INFINITE :
                       _mulp(lead_num, _invp(lead_den)))
      end
      found[1] == _INFINITE || found[1] != found[2] ? _UNDEFINED : found[1]
  end

  """
      _face_limit_check(mech, net, s, face, seed) -> _LimitVerdict

  The limit test of class 3 (research section 1.4) on face `face` of `mech`, whose law
  network is `net` and face structure `s`: at two random points of the model mod 2^61 − 1,
  drawn from `seed`, the fast-RE limit of the face chain (by aggregation when the aggregated
  chain has one closed class, else by interpolation in Λ) against the law's face value
  (from the law's segment chain when it has one closed class, else by reconstruction as the
  zeroed concentrations vanish). On a difference a third limit at new speeds tells
  `_LIMIT_DIFFER` from `_LIMIT_SPEED`. A false `_LIMIT_EQUAL` has probability at most twice
  the degree over p.
  """
  function _face_limit_check(mech::_AnyMech, net::_Network, s, face::UInt64, seed::UInt64)
      m, tags = catalytic(mech), _tags(mech)
      rng = _ResidueStream(seed)
      x1 = _turnover_name(m)
      h = _fast_potential(net, x1)
      limit(pt) = length(s.aggregated) == 1 ? _limit_by_aggregation(m, net, s, pt, h, x1) :
                  _limit_by_interpolation(m, net, s, pt, h, x1, rng)
      for _ in 1:2
          pt = _model_point(m, tags, net.ncopies, face, rng)
          law = length(s.law_classes) == 1 ? _law_on_face(m, net, s, pt, h, x1) :
                _law_by_reconstruction(m, net, s, pt, h, x1, face, rng)
          law == _UNDEFINED && return _LIMIT_LAW_NAN
          v = limit(pt)
          # The face chain's fast transitions are the law network's, a subset of the
          # ε-network's, where class 2, which runs first, finds no turnover cycle.
          v == _INFINITE && error("_face_limit_check: the limit grows without bound")
          v == law && continue
          return limit(_respeed_point(pt, rng)) == v ? _LIMIT_DIFFER : _LIMIT_SPEED
      end
      _LIMIT_EQUAL
  end
  ```

- [ ] **Step 5: Add the class-3 messages.** In `_guard_message`, replace its last line
  `error("_guard_message: verdict $v holds no fault")` with:
  ```julia
      v.fault === _NO_FAULT && error("_guard_message: verdict $v holds no fault")
      _face_message(who, mech, v)
  ```
  and after `_re_turnover_message` add:
  ```julia
  """Where a class-3 fault of `mech` lies: nothing for a plain mechanism; the conformation
  `c` of an allosteric one."""
  _place(::Mechanism, ::Int) = ""
  _place(::AllostericMechanism, c::Int) =
      c == 2 ? " in the inactive conformation" : c == 1 ? " in the active conformation" : ""

  """
  The message of a class-3 verdict `v` on `mech` (research section 1.11): the face, the
  forms and steps at fault, and a fix only for a trap and a bottomless segment, whose
  variants the tests show passing.
  """
  function _face_message(who::String, mech::_AnyMech, v::_Verdict)
      m = catalytic(mech)
      net = _law_network(mech)
      face = _faces(m)[v.face]
      s = _face_structure(net, face)
      at = "at zero " * join(String.(_face_names(m, face)), ", ", " and ") *
           _place(mech, v.conformation)
      forms(nodes) = "{" * join([_form_text(m, _form_of(net, u)) for u in nodes], ", ") * "}"
      lines(edges) = join([_step_text(m, k) for k in edges], ", ", " and ")
      if v.fault === _NONUNIQUE
          C = _away_class(s)
          return "$who: $at, the forms $(forms(C)) keep enzyme away from free E: every " *
                 "step that leaves them takes up a zeroed metabolite " *
                 "($(lines(_exits(net, s, C)))), so the steady state is not unique."
      elseif v.fault === _BOTTOMLESS
          return "$who: $at, the rapid-equilibrium segment " *
                 "$(forms(s.segments[s.bottomless])) has no form free of the zeroed " *
                 "metabolites, so the rate is 0/0 there. Make one of the segment's binding " *
                 "steps steady-state."
      elseif v.fault === _TRAP
          trap, edge = _trap_feeder(net, s)
          return "$who: $at, the rapid-equilibrium step $(_step_text(m, edge)) fills the " *
                 "forms $(forms(s.fast_classes[trap])), which no rapid-equilibrium step " *
                 "running there empties, so the rate law misweights them. Make that " *
                 "step's kinetic group steady-state."
      elseif v.fault === _SPLIT
          split = _split_parts(net, s)
          return "$who: $at, the forms $(join(forms.(split.parts), " and ")) of one " *
                 "rapid-equilibrium segment are joined only by steps that take up a " *
                 "zeroed metabolite ($(lines(split.joins))), so the rate law cannot weight " *
                 "one against the other."
      elseif v.fault === _LAW_MULTI || v.fault === _MULTILEVEL
          why = v.limit === _LIMIT_LAW_NAN ?
                    "the rate law has no value (0/0), while the steady state does" :
                v.limit === _LIMIT_SPEED ?
                    "the steady state depends on how fast the rapid-equilibrium steps " *
                    "run, which the rate law does not carry" :
                    "the steady state's fast-equilibrium limit differs from the rate law"
          return "$who: $at, $why."
      end
      error("_face_message: no message for $(v.fault)")
  end
  ```

- [ ] **Step 6: Run the tests and watch them pass.** The Task 2 Step 2 command for
  `test/test_admission.jl`, log `P4_T3_green.log`. Expected: `exit 0`; all sixteen testsets
  pass, Task 2's seven among them (the fix variants of the class-2 testset now pass class 3
  as well); no warning. A failure of `_testhelper_agrees` on a face is a guard-oracle
  disagreement: stop and run `superpowers:systematic-debugging` on that face (compare
  `_face_structure` with `rules/lib.jl`'s `face_structure` on the same network) before any
  change; never loosen the assertion. A message pin that differs while the verdicts agree is
  checked against RES §1.11's wording and the structure, then re-pinned with the reason in
  the commit body. The refresh pass (P4-32) ran Task 3's testsets at 08c4043 with this
  plan's code (`$S/phase45/refresh/map-head-t0-3_task23.py`, log `map-head-t0-3_7.log`),
  and each passed, its unmeasured pins among them (Open item 1): the per-face fault
  vectors of `ldh_like_trap`, `trap_flip_child9`, `split_flip_child1` and
  `bottomless_entered`, `_LIMIT_SPEED` for the LDH trap, `lumpable_split_532`'s faces,
  C3b's and C3c's stages, the face order, and the bottomless, non-unique and split messages
  with the split message's form order; this run repeats that check on the branch, and the
  oracle-agreement assertion on each face confirms any pin before anything changes.

- [ ] **Step 7: Measure the cost of class 3 per face.** Write
  `$S/phase45/probes/P4_T3_cost.jl`:
  ```julia
  # ABOUTME: Times the guard's class 3 per face over the shared test specs, the cost Task 3
  # ABOUTME: of the phase-4 plan reports (prototype 6.6 µs, target at most 3 µs).
  using Test, EnzymeRates, LinearAlgebra, Random
  # `include` resolves a relative path against this file's folder; the run's working
  # directory is the worktree.
  include(joinpath(pwd(), "test", "mechanism_definitions_for_test_enzyme_derivation.jl"))
  mechs = [_testhelper_lift(s.mechanism) for s in MECHANISM_TEST_SPECS]
  cases = [(m, EnzymeRates._law_network(m), EnzymeRates._faces(EnzymeRates.catalytic(m)))
           for m in mechs]
  check(m, net, faces) = [EnzymeRates._check_face(m, net, faces, i) for i in eachindex(faces)]
  verdicts = [check(c...) for c in cases]
  nfaces = sum(length, verdicts)
  tested = count(v -> v.limit !== EnzymeRates._LIMIT_UNTESTED, Iterators.flatten(verdicts))
  structural = [(m, net, faces) for (m, net, faces) in cases
                if all(v -> v.limit === EnzymeRates._LIMIT_UNTESTED, check(m, net, faces))]
  t = @elapsed for _ in 1:50, (m, net, faces) in structural
      check(m, net, faces)
  end
  n = sum(length(f) for (_, _, f) in structural)
  println("class 3: ", round(1e6 * t / (50 * n); digits = 2), " µs per face over ", n,
          " faces the structure decides; ", tested, " of ", nfaces,
          " spec faces reach the limit test")
  ```
  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T3_cost.log 900 --project=. -e 'using TestEnv; TestEnv.activate(); include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/probes/P4_T3_cost.jl")'
  ```
  Expected: `exit 0` and one line `class 3: <y> µs per face over <n> faces the structure
  decides; <z> of <N> spec faces reach the limit test`, `<z>` small (the two trap specs among
  them). Record the line for the commit body. Above 6.6 µs, note in the commit body that
  Task 8's ter-ter gate (150 s) is at risk (RES §4 risk 2) and where the time goes
  (`@profile` on the same loop).

- [ ] **Step 8: Run the full suite.** The Task 1 Step 14 command, log `P4_T3_full.log`.
  Expected: `exit 0`, the previous count plus this task's tests, Aqua and JET pass, no
  warning, maxrss within 5% of M0; `test/test_admission.jl` takes at most 3 s more than after
  Task 2 (its summary line).

- [ ] **Step 9: Commit.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && git add src/admission.jl test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_admission.jl && git commit -F - <<'EOF'
  Check class 3 on every face with procedure A and a mod-p limit test

  On each face of the spec's list, procedure A rejects a second closed class or a segment
  with no lowest form, accepts a face whose structure fixes the limit, and otherwise lets a
  test mod 2^61 - 1, seeded by the content hash, compare the fast-RE limit with the law.
  Messages name the trap's feeding step and the split's joining steps. A trap entered by a
  slow step names the rapid-equilibrium step that joins it to its segment (research
  section 1.11), where the research script named the slow step. Not wired yet.
  Class 3 costs <y> µs per face; <z> of <N> spec faces reach the limit test.

  Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
  EOF
  ```
  with the figures from Step 7.

---

### Evidence for Tasks 4-7

The expected verdicts below were measured, not inferred. `$S/phase45/plan/probes/p2_verdicts.jl`
ran research §1's rules (`rules/lib.jl`: classes 1-3, interior admission, B v1, plus B v3's
ε-network count and joint R1 from `vi/fix_b1_lib.jl` and `verify_claim/joint_r1.jl`) and both
research oracles on every fixture of these tasks (log `$S/phase45/plan/probes/p2_verdicts.log`,
`p2_verdicts2.log`; 17 s). `p2_counts.jl` (log `p2_counts2.log`) counted the trap variants'
constraints and timed their `rate_equation` on phase 2. `p2_fz101.jl` regenerated
`fz101_10232` and `fz101_1964` from `vi/fuzz.jl` (seed 101, Julia 1.12.5; `fz101_1964` matches
VI:119-126 exactly). Findings that shape the tasks:

- The ε-path oracle does not see class-3 failures: it converges on the LDH spec's trap face.
  Class-3 rejections are therefore checked against the RE-speed oracle, admission's against the
  ε-path (`_testhelper_admission_agreement`, Task 5).
- RES §1.11's `:arc_cycle` fix ("untag one of them") fails on the smallest arc cycle (an
  ordered bi-uni with both bindings `:OnlyA`): untagging either binding leaves the other inside
  one component. The message offers instead the two fixes the probe confirms: untag them all,
  or tag the steps that join them.
- R2funnel's bridging binding made steady-state (R2funnel_SS) is still a bridge, while
  ADbridge's and bridge_deadend's become well-posed. The `:bridge` message therefore offers the
  steady-state fix only after checking, on that face, that the variant passes.
- Under B v3, S1 fails admission's face condition too (`_EPS_NONUNIQUE` at substrates = 0),
  not only class 3; the guard reports class 3, which runs first on each face.
- AD6 cannot be built after phase 3: its two edges between E(S) and E(P) are one reaction in
  two groups (`_check_each_reaction_once`). Research's AD2a_iso_chemNonequal is R1, AD2e is R2
  and AD2c is R3 written again; each is tested once.
- The variants' counts: spec 14's 7 independent parameters (Haldane 1, Wegscheider 1), LDH's 10
  (Haldane 1, Wegscheider 2). `rate_equation`: spec 14's variant 9.1 ns, LDH's 65.4 ns (the
  spec 20.1 ns), 0 allocations, on phase 2 locally.
- Every fixture these tasks write builds at 08c4043, 75 unguarded calls and the 10 AD13
  sets, and stores each edge as the S3f snapshot did: the refresh probe
  `$S/phase45/refresh/map-head-t0-3_probe.jl` (log `map-head-t0-3_2.log`) repeated
  `$S/phase45/plan/probes/p2check_build.jl` (log `p2check_build.log`), which prints each
  edge's stored text, and printed the same 46 stored texts. A step that takes up or gives
  off one metabolite is stored as its uptake, so the ping-pong's half-reactions are stored
  `E(; residual = A - P) + P <--> E(A)` and `E + Q <--> E(B; residual = A - P)`. The island
  message therefore reads each tagged chemistry step in the direction that runs toward
  products (`_runs_toward_products`): it names half 1 of a both-halves ping-pong in either
  writing, the half whose untagging leaves a well-posed mechanism for AGT:663 too (RES §5.2
  item 1, P4-29b). By stored direction it would name half 2 here and half 1 in AGT:663.
- `$S/phase45/plan/probes/p2check_rules.jl` (log `p2check_rules.log`) loads this plan's code
  into the S3f snapshot (Task 2's networks, faces and messages, Task 3's closed classes,
  Tasks 4-6) and checks, without class 3 or the oracles, every admission verdict, named
  edge, admission fault per face and message string Tasks 4-6 pin, the hand-built joint-R1
  network and class 5 on the specs: 461 tests pass. The two trap specs keep a monomial on
  every face too. Rerun at 08c4043 on the assembled plan
  (`$S/phase45/refresh/map-head-t0-3_rules.jl`, log `map-head-t0-3_3.log`), the same 461
  tests pass with the same output. `p2check_fixes.jl` (log `p2check_fixes.log`, research
  rules and both
  oracles) shows the no-limit message's two fixes well-posed on AD1 SS and AD3c as well as on
  the uni-uni fixture: every face accepts and the ε-path converges. `nolimit_fixes.jl` (log
  `nolimit_fixes.log`) shows the same for AD11's and the ter cube's two fixes.

---

### Task 4: Admission interior

Admission's existence and uniqueness checks per excluding conformation (RES §1.5 ADMIT, spec
§7 "Admission" items 1-2), wired into `_guard(am)` after class 2. The turnover check (item
3) needs a catalytic `:OnlyI` group, which no constructor accepts before Task 10, so Task 10
adds it after watching its test fail.

**Files:**
- Modify: `src/admission.jl` (`@enum _Fault`; new interior functions and messages; `_guard(am)`;
  `_guard_message`)
- Test: `test/test_admission.jl`

**Interfaces:**
- Consumes (Task 2): `_excludes(am::AllostericMechanism, k::Int, c::Int)::Bool`,
  `_law_mask(am::AllostericMechanism, c::Int)::BitVector`,
  `_has_turnover_cycle(m::Mechanism, mask::BitVector)::Bool`, `_Verdict`, `_WELL_POSED`,
  `_LIMIT_UNTESTED`, `_re_turnover_verdict`, `_eps_network(am)`, `_law_network(am)`,
  `_guard(am)`, `_guard_message(mech::_AnyMech, v)`. (Task 3) `_faces_verdict(mech, net)`,
  the `_Fault` enum and its docstring. (Task 1) `@_testhelper_unguarded` (used as
  `@_testhelper_unguarded @allosteric_mechanism begin ... end`),
  `_testhelper_uniuni_grid()`, `_testhelper_lp_existence(am, c)`, `_testhelper_exclusion(am,
  c)`, `_testhelper_epath_verdict(m; rng, face, corners, precision)`. (Phase 3) the read API:
  `_is_plain_binding`, `_is_chemistry`, `_role_masks(m)`, `total`, `_form_name`,
  `_step_text(m::Mechanism, k::Int)`, `_form_text(m::Mechanism, f::Int)` (names.jl);
  `_testhelper_edge(m, text)`, whose text writes a reactant's modifier copy as its label
  (`Smod`, `_entry_label`).
- Produces:
  ```julia
  # _Fault gains, after _MULTILEVEL: _NO_LIMIT _ARC_CYCLE _ISLAND
  const _ADMISSION_FAULTS::NTuple{3, _Fault}    # six from Task 5, seven from Task 10
  _exclusion_components(am::AllostericMechanism, c::Int)::Vector{Int}  # per form: its component's smallest form
  _arcs(am::AllostericMechanism, c::Int)::Vector{Int}                  # plain bindings c excludes
  _first_excluded(am::AllostericMechanism, c::Int)::Int                # 0 when c excludes nothing
  _arc_cycle(m::Mechanism, comp::Vector{Int}, arcs::Vector{Int})::Vector{Int}   # cycle order; empty
  _reached_components(m::Mechanism, comp::Vector{Int}, arcs::Vector{Int})::BitVector
  _runs_toward_products(m::Mechanism, e::Edge)::Bool  # as stored, from the substrate-richer side
  _island_edge(am::AllostericMechanism, c::Int, comp::Vector{Int}, reached::BitVector)::Int
  _kept_path(am::AllostericMechanism, c::Int, a::Int, b::Int)::Vector{Int}      # edges, a to b
  _admission_interior(am::AllostericMechanism)::_Verdict   # face 0; conformation and edge set
  _admission_message(am::AllostericMechanism, v::_Verdict)::String
  const _testhelper_tag_sets   # test_admission.jl: Dict{String, AllostericMechanism}
  ```
  `_guard(am)` runs `_admission_interior` after class 2 and before the faces; `_guard_message`
  hands an allosteric admission fault to `_admission_message`.

**Gate:** the graph test equals `_testhelper_lp_existence` on 108 of 108 grid mechanisms, and
its fault equals `_testhelper_exclusion`'s on each; AD1 SS, AD1 RE, AD11, AD3c, the ter cube and
uni_S_OnlyA are `_NO_LIMIT` (ε-path `:no_path`); biuni_arc_cycle is `_ARC_CYCLE`; R3, R4, AD2b,
pp_V_both_halves, pp_K_Bbind_half1 are `_ISLAND` (ε-path `:direction_dependent` or
`:wrong_limit`); every interior-admitted row of RES §3.2 passes; every printed fix's variant
passes `_guard` (for the no-limit fix: the uni-uni, AD1, AD3c, AD11 and the ter cube, each
variant checked by the research rules and both oracles, `nolimit_fixes.log` and
`p2check_fixes.log`); mixed_group_limit is admitted and converges (Review Focus 4);
test_admission.jl grows by at most 1 s; the interior's cost per mechanism is reported
(prototype 2.5 µs).

- [ ] **Step 1: Write the failing tests**

Append to `test/test_admission.jl`:

```julia
# ── Tag sets admission's interior accepts ───────────────────────────────────

"""
The tag sets of research section 3.2 that admission's interior accepts, by name, each built
without the guard (`@_testhelper_unguarded`): those that pass every face, and those a face
rejects. The research's covalent ping-pong form `F` is written
`E(; residual = A - P)`, the same graph; AD2a_iso_chemNonequal is R1, AD2e is R2 and AD2c is
R3 written again, and AD6 holds one reaction in two groups, which the constructor rejects.
"""
const _testhelper_tag_sets = Dict{String, EnzymeRates.AllostericMechanism}()

_testhelper_tag_sets["R1"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> E(P)      :: NonequalAI
        F + P ⇌ E(P)        :: EqualAI
        F <--> E            :: NonequalAI
    end
end
_testhelper_tag_sets["A5"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> E(P)      :: NonequalAI
        E + P ⇌ E(P)        :: EqualAI
        F ⇌ E               :: NonequalAI
    end
end
_testhelper_tag_sets["AD2a"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> E(P)      :: OnlyA
        F + P ⇌ E(P)        :: EqualAI
        F <--> E            :: EqualAI
    end
end
_testhelper_tag_sets["AD2d"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> E(P)      :: EqualAI
        F + P ⇌ E(P)        :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["AD3a"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> E(P)      :: NonequalAI
        E + P ⇌ E(P)        :: EqualAI
        E + I ⇌ E(I)        :: OnlyA
    end
end
_testhelper_tag_sets["AD3b"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: S
    steps: begin
        E + S ⇌ E(S)                                       :: EqualAI
        E(S) <--> E(P)                                     :: EqualAI
        E + P ⇌ E(P)                                       :: EqualAI
        (E + S::Mod ⇌ E(S::Mod), E(S) + S::Mod ⇌ E(S, S::Mod)) :: OnlyA
        E(S::Mod) + S ⇌ E(S, S::Mod)                       :: EqualAI
    end
end
_testhelper_tag_sets["AD3d"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: S
    steps: begin
        E + S ⇌ E(S)                    :: EqualAI
        E(S) <--> E(P)                  :: EqualAI
        E + P ⇌ E(P)                    :: EqualAI
        E + S::Mod <--> E(S::Mod)       :: OnlyA
        E(S) + S::Mod <--> E(S, S::Mod) :: OnlyA
        E(S::Mod) + S ⇌ E(S, S::Mod)    :: EqualAI
    end
end
_testhelper_tag_sets["AD4"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P, Q
    steps: begin
        E + S ⇌ E(S)            :: EqualAI
        E(S) <--> E(P, Q)       :: EqualAI
        E(P, Q) <--> E + P + Q  :: OnlyA
    end
end
_testhelper_tag_sets["AD4b"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P, Q
    steps: begin
        E + S ⇌ E(S)            :: EqualAI
        E(S) <--> E(P, Q)       :: OnlyA
        E(P, Q) ⇌ E + P + Q     :: EqualAI
    end
end
_testhelper_tag_sets["AD5"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E <--> F                        :: EqualAI
        (E + S ⇌ E(S), F + S ⇌ F(S))    :: OnlyA
        E(S) <--> F(S)                  :: EqualAI
        F(S) <--> F(P)                  :: OnlyA
        F + P ⇌ F(P)                    :: EqualAI
    end
end
_testhelper_tag_sets["AD7a"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A
    products: P
    modifiers: I
    steps: begin
        E + A <--> E(A)                 :: EqualAI
        (E + I ⇌ E(I), F + I ⇌ E(I))    :: OnlyA
        E + P <--> E(P)                 :: EqualAI
        E ⇌ F                           :: EqualAI
        E(A) <--> E(P)                  :: EqualAI
    end
end
_testhelper_tag_sets["AD7b"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A
    products: P
    modifiers: I
    steps: begin
        E + A <--> E(A)                 :: EqualAI
        (E + I ⇌ E(I), F + I ⇌ E(I))    :: EqualAI
        E + P <--> E(P)                 :: EqualAI
        E ⇌ F                           :: OnlyA
        E(A) <--> E(P)                  :: EqualAI
    end
end
_testhelper_tag_sets["AD7c"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A
    products: P
    modifiers: I
    steps: begin
        E + A <--> E(A)                 :: EqualAI
        (E + I ⇌ E(I), F + I ⇌ E(I))    :: EqualAI
        E + P <--> E(P)                 :: EqualAI
        E ⇌ F                           :: EqualAI
        E(A) <--> E(P)                  :: OnlyA
    end
end
_testhelper_tag_sets["AD8"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ E(A)            :: EqualAI
        E(A) + B <--> E(Q) + P  :: OnlyA
        E + Q ⇌ E(Q)            :: EqualAI
        E(Q) + B ⇌ E(B, Q)      :: EqualAI
    end
end
_testhelper_tag_sets["AD9"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ E(A)            :: EqualAI
        E(A) + B ⇌ E(A, B)      :: EqualAI
        E(A, B) ⇌ E(P, Q)       :: OnlyA
        E(Q) + P <--> E(P, Q)   :: EqualAI
        E + Q <--> E(Q)         :: EqualAI
        E(A) + Q <--> E(A, Q)   :: OnlyA
    end
end
_testhelper_tag_sets["AD10"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)         :: OnlyA
        E + B ⇌ E(B)            :: EqualAI
        E(A) + B ⇌ E(A, B)      :: EqualAI
        E(B) + A <--> E(A, B)   :: OnlyA
        E(A, B) <--> E(P, Q)    :: OnlyA
        E(Q) + P ⇌ E(P, Q)      :: EqualAI
        E + Q ⇌ E(Q)            :: EqualAI
    end
end
# AD13: a second free form joined to E at rapid equilibrium, each label sorting before or
# after E, with the chemistry :NonequalAI or :OnlyA.
for label in (:D, :Ca, :Ea, :F, :Z), chem in (:NonequalAI, :OnlyA)
    block = quote
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: EqualAI
            E(S) <--> E(P)      :: $chem
            E + P ⇌ E(P)        :: EqualAI
            $label ⇌ E          :: NonequalAI
        end
    end
    _testhelper_tag_sets["AD13_$(chem)_$label"] =
        EnzymeRates.AllostericMechanism(EnzymeRates._raw_mechanism(block, true))
end
_testhelper_tag_sets["AD14"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ E(S)                            :: EqualAI
        E + P ⇌ E(P)                            :: EqualAI
        (E(S) <--> E(P), E(S, I) <--> E(P, I))  :: OnlyA
        E(I) + S ⇌ E(S, I)                      :: OnlyA
        E(I) + P ⇌ E(P, I)                      :: EqualAI
        E + I ⇌ E(I)                            :: EqualAI
    end
end
_testhelper_tag_sets["AD14b"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ E(S)                            :: EqualAI
        E + P ⇌ E(P)                            :: EqualAI
        (E(S) <--> E(P), E(S, I) <--> E(P, I))  :: OnlyA
        E(I) + S <--> E(S, I)                   :: OnlyA
        E(I) + P ⇌ E(P, I)                      :: EqualAI
        E + I ⇌ E(I)                            :: EqualAI
    end
end
_testhelper_tag_sets["AD14c"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ E(S)                            :: EqualAI
        E + P ⇌ E(P)                            :: EqualAI
        (E(S) <--> E(P), E(S, I) <--> E(P, I))  :: OnlyA
        E(I) + S ⇌ E(S, I)                      :: EqualAI
        E(I) + P ⇌ E(P, I)                      :: EqualAI
        E + I ⇌ E(I)                            :: EqualAI
    end
end
_testhelper_tag_sets["AD15"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))  :: OnlyA
        E(A) + B ⇌ E(A, B)                  :: OnlyA
        E + B ⇌ E(B)                        :: EqualAI
        E(A, B) <--> E(P, Q)                :: OnlyA
        E(Q) + P ⇌ E(P, Q)                  :: EqualAI
        E + Q ⇌ E(Q)                        :: EqualAI
    end
end
_testhelper_tag_sets["R2P_sealed"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> E + P     :: EqualAI
        E + P ⇌ E(P)        :: OnlyA
        F + P ⇌ E(P)        :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["B_uphill"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: X
    steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> E(P)      :: EqualAI
        E + P ⇌ E(P)        :: EqualAI
        E + X ⇌ E(X)        :: OnlyA
        F + X <--> E(X)     :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["TC1"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ E(A)            :: EqualAI
        E(A) + B ⇌ E(Q) + P     :: OnlyA
        E(Q) <--> E + Q         :: EqualAI
        E(A) + B <--> E(A, B)   :: EqualAI
        E(A, B) <--> E(P, Q)    :: EqualAI
        E(P, Q) ⇌ E(Q) + P      :: EqualAI
    end
end
_testhelper_tag_sets["pp_V_half1"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)                                       :: EqualAI
        E(A) <--> E(; residual = A - P) + P                   :: OnlyA
        E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
        E(B; residual = A - P) <--> E + Q                     :: EqualAI
    end
end
_testhelper_tag_sets["pp_V_half2"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)                                       :: EqualAI
        E(A) <--> E(; residual = A - P) + P                   :: EqualAI
        E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
        E(B; residual = A - P) <--> E + Q                     :: OnlyA
    end
end
_testhelper_tag_sets["pp_K_Abind_half1"] =
    @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)                                       :: OnlyA
        E(A) <--> E(; residual = A - P) + P                   :: OnlyA
        E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
        E(B; residual = A - P) <--> E + Q                     :: EqualAI
    end
end
_testhelper_tag_sets["PP1"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ E(A)                                          :: EqualAI
        E(A) <--> E(; residual = A - P) + P                   :: OnlyA
        E(; residual = A - P) + B ⇌ E(B; residual = A - P)    :: EqualAI
        E(B; residual = A - P) <--> E + Q                     :: EqualAI
    end
end
_testhelper_tag_sets["PP2"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    modifiers: I
    steps: begin
        E + A <--> E(A)                                       :: EqualAI
        E(A) <--> E(; residual = A - P) + P                   :: EqualAI
        E(; residual = A - P) + B ⇌ E(B; residual = A - P)    :: EqualAI
        E(B; residual = A - P) <--> E + Q                     :: OnlyA
        E(; residual = A - P) + I ⇌ E(I; residual = A - P)    :: OnlyA
    end
end
# Census ter-bi d75d2415361c4dae (`phase3/census/S3f/mechs.tsv`): its faces with P zero pass
# by T0+ alone.
_testhelper_tag_sets["terbi_d75d2415361c4dae"] =
    @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B, D
    products: P, Q
    steps: begin
        (E + A ⇌ E(A), E(B, P) + A ⇌ E(A, B, P), E(P) + A ⇌ E(A, P))       :: OnlyA
        (E + P ⇌ E(P), E(A) + P ⇌ E(A, P), E(A, B) + P ⇌ E(A, B, P),
         E(Q) + P ⇌ E(P, Q))                                               :: EqualAI
        (E + Q <--> E(Q), E(P) + Q <--> E(P, Q))                           :: OnlyA
        (E(A) + B <--> E(A, B), E(A, P) + B <--> E(A, B, P),
         E(P) + B <--> E(B, P))                                            :: EqualAI
        (E(A, B) + D <--> E(P, Q), E(Q) + D <--> E(D, Q))                  :: OnlyA
    end
end
# An :OnlyA group mixing B's fused binding into E(P, Q) with its plain binding at E(Q),
# which no untagged step joins to E(B, Q): its limit exists (spec section 7, "Mixed groups
# stay allowed").
_testhelper_tag_sets["mixed_group_limit"] =
    @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ E(A)                                    :: EqualAI
        (E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q))  :: OnlyA
        E(Q) + P ⇌ E(P, Q)                              :: EqualAI
        E + Q ⇌ E(Q)                                    :: EqualAI
    end
end
_testhelper_tag_sets["R2"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)        :: OnlyA
        E(S) <--> E(P)      :: EqualAI
        F + P ⇌ E(P)        :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["BS1"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> E(P)      :: EqualAI
        E(P) ⇌ E + P        :: EqualAI
        E <--> F            :: OnlyA
        E + I ⇌ F(I)        :: OnlyA
        F + I ⇌ F(I)        :: EqualAI
        F + S ⇌ F(S)        :: EqualAI
        F(S) <--> F(P)      :: EqualAI
        F(P) ⇌ F + P        :: EqualAI
    end
end
_testhelper_tag_sets["BS2"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)        :: OnlyA
        E(S) <--> E(P)      :: EqualAI
        F + P ⇌ E(P)        :: EqualAI
        F ⇌ E               :: OnlyA
    end
end
_testhelper_tag_sets["BS3"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S <--> E(S)     :: OnlyA
        E(S) <--> E(P)      :: EqualAI
        F + P ⇌ E(P)        :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["BS4"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ E(S)            :: EqualAI
        E(S) <--> E + P         :: EqualAI
        E + I ⇌ E(I)            :: OnlyA
        E(I) + S <--> E(I, S)   :: EqualAI
        E(I, S) <--> F(I) + P   :: EqualAI
        F(I) ⇌ F + I            :: EqualAI
        F <--> E                :: OnlyA
    end
end
_testhelper_tag_sets["ADbridge"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S <--> E(S)                     :: EqualAI
        E(S) <--> E + P                     :: EqualAI
        (E + I ⇌ E(I), E(S) + I ⇌ E(I, S))  :: OnlyA
        E(I) + S ⇌ E(I, S)                  :: EqualAI
    end
end
_testhelper_tag_sets["ADbridge_noturnover"] =
    @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S <--> E(S)                     :: EqualAI
        E(S) <--> E + P                     :: OnlyA
        (E + I ⇌ E(I), E(S) + I ⇌ E(I, S))  :: OnlyA
        E(I) + S ⇌ E(I, S)                  :: EqualAI
    end
end
_testhelper_tag_sets["bridge_drain"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> F + P     :: EqualAI
        F <--> E            :: OnlyA
        E + I ⇌ E(I)        :: OnlyA
        F + I ⇌ F(I)        :: OnlyA
        F(I) ⇌ E(I)         :: EqualAI
    end
end
_testhelper_tag_sets["bridge_deadend"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I, J
    steps: begin
        E + S <--> E(S)         :: EqualAI
        E(S) <--> E + P         :: EqualAI
        E + J <--> E(J)         :: EqualAI
        E + I ⇌ E(I)            :: OnlyA
        E(J) + I ⇌ E(I, J)      :: OnlyA
        E(I) + J ⇌ E(I, J)      :: EqualAI
    end
end
_testhelper_tag_sets["R2funnel"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)        :: OnlyA
        E(S) <--> E(P)      :: EqualAI
        F + P ⇌ E(P)        :: EqualAI
        F <--> E            :: OnlyA
        F + S ⇌ F(S)        :: EqualAI
        F(S) <--> F(P)      :: EqualAI
        F + P ⇌ F(P)        :: EqualAI
        F(P) ⇌ E(P)         :: EqualAI
    end
end
_testhelper_tag_sets["R2funnel_SS"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S <--> E(S)     :: OnlyA
        E(S) <--> E(P)      :: EqualAI
        F + P ⇌ E(P)        :: EqualAI
        F <--> E            :: OnlyA
        F + S ⇌ F(S)        :: EqualAI
        F(S) <--> F(P)      :: EqualAI
        F + P ⇌ F(P)        :: EqualAI
        F(P) ⇌ E(P)         :: EqualAI
    end
end
_testhelper_tag_sets["S1"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P
    modifiers: A
    steps: begin
        E + A <--> E(A)                 :: OnlyA
        E(A) + B <--> E(A, B)           :: EqualAI
        E(A, B) <--> E(P)               :: OnlyA
        E(P) <--> E + P                 :: EqualAI
        E + A::Mod ⇌ E(A::Mod)          :: OnlyA
        E(A::Mod) + B ⇌ E(A::Mod, B)    :: EqualAI
        E(B) + A::Mod ⇌ E(A::Mod, B)    :: EqualAI
        F + B <--> E(B)                 :: EqualAI
    end
end
# Fuzz mechanism fz101_10232 (`vi/fuzz.jl`, seed 101).
_testhelper_tag_sets["fz101_10232"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S <--> E(S)                             :: OnlyA
        E(S) <--> E(P)                              :: OnlyA
        E(P) <--> F + P                             :: EqualAI
        F <--> E                                    :: OnlyA
        (E(S) + I <--> E(I, S), F + I <--> F(I))    :: OnlyA
        E(P) + S ⇌ E(S) + P                         :: EqualAI
    end
end

@testset "the graph test equals the exact LP on the uni-uni grid" begin
    # Existence holds exactly when no arc stays inside a component and the arcs close no
    # cycle (spec section 7); the test oracle's exact simplex decides the same question
    # from the cycle rows, and its own exclusion structure names the same failure.
    grid = _testhelper_uniuni_grid()
    @test length(grid) == 108
    for am in grid
        v = EnzymeRates._admission_interior(am)
        exists = v.fault ∉ (EnzymeRates._NO_LIMIT, EnzymeRates._ARC_CYCLE)
        @test exists == _testhelper_lp_existence(am, 2)
        ex = _testhelper_exclusion(am, 2)
        @test v.fault === (ex.failure === :inside ? EnzymeRates._NO_LIMIT :
                           ex.failure === :cycle ? EnzymeRates._ARC_CYCLE :
                           isempty(ex.lost) ? EnzymeRates._NO_FAULT : EnzymeRates._ISLAND)
    end
end

@testset "a tagged binding that untagged steps join on both sides has no limit" begin
    rng = Random.MersenneTwister(31)
    NL, U = EnzymeRates._NO_LIMIT, EnzymeRates._LIMIT_UNTESTED
    # The uni-uni S binding :OnlyA: the chemistry and the P binding join E to E(S) in the
    # inactive conformation, so the S affinity cannot vanish there.
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: OnlyA
            E(S) <--> E(P)      :: EqualAI
            E + P ⇌ E(P)        :: EqualAI
        end
    end
    m = am.catalytic
    bind, release, chem = (_testhelper_edge(am, t) for t in
                           ("E + S ⇌ ES", "E + P ⇌ EP", "ES <--> EP"))
    v = EnzymeRates._guard(am)
    @test v == EnzymeRates._Verdict(NL, U, 0, 2, bind)
    @test EnzymeRates._admission_interior(am) == v
    @test _testhelper_epath_verdict(am; rng) === :no_path
    line(k) = EnzymeRates._step_text(m, k)
    @test EnzymeRates._guard_message(am, v) ==
          "AllostericMechanism: in the inactive conformation, steps outside :OnlyA " *
          "groups still join the two sides of the :OnlyA binding $(line(bind)), so " *
          "thermodynamics keeps its affinity finite and the binding has no limit " *
          "there. Tag its group :EqualAI or :NonequalAI, or tag :OnlyA the steps that " *
          "join its sides: $(line(release)); $(line(chem))."
    # Both fixes pass the guard: the binding untagged, and the steps that join its sides
    # tagged.
    untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: EqualAI
            E(S) <--> E(P)      :: EqualAI
            E + P ⇌ E(P)        :: EqualAI
        end
    end
    joins_tagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: OnlyA
            E(S) <--> E(P)      :: OnlyA
            E + P ⇌ E(P)        :: OnlyA
        end
    end
    @test EnzymeRates._guard(untagged) == EnzymeRates._WELL_POSED
    @test EnzymeRates._guard(joins_tagged) == EnzymeRates._WELL_POSED
    # AD1 and AD11: a mixed :OnlyA group whose plain member's two sides untagged steps
    # still join (spec section 7: AD1 and AD11 are rejected). AD3c: a ::Mod copy's binding.
    ad1_ss = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                    :: EqualAI
            (E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q))  :: OnlyA
            E(Q) + P ⇌ E(P, Q)                              :: EqualAI
            E + Q ⇌ E(Q)                                    :: EqualAI
            E + B ⇌ E(B)                                    :: EqualAI
            E(B) + Q ⇌ E(B, Q)                              :: EqualAI
        end
    end
    ad1_re = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                :: EqualAI
            (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))    :: OnlyA
            E(Q) + P <--> E(P, Q)                       :: EqualAI
            E + Q ⇌ E(Q)                                :: EqualAI
            E + B ⇌ E(B)                                :: EqualAI
            E(B) + Q ⇌ E(B, Q)                          :: EqualAI
        end
    end
    ad11 = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                            :: EqualAI
            (E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q),
             E(P) + B <--> E(B, P))                                 :: OnlyA
            E(Q) + P ⇌ E(P, Q)                                      :: EqualAI
            E + Q ⇌ E(Q)                                            :: EqualAI
            E + P ⇌ E(P)                                            :: EqualAI
            E + B ⇌ E(B)                                            :: EqualAI
            E(B) + Q ⇌ E(B, Q)                                      :: EqualAI
        end
    end
    ad3c = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: S
        steps: begin
            E + S ⇌ E(S)                    :: EqualAI
            E(S) <--> E(P)                  :: OnlyA
            E + P ⇌ E(P)                    :: EqualAI
            E + S::Mod ⇌ E(S::Mod)          :: OnlyA
            E(S) + S::Mod ⇌ E(S, S::Mod)    :: EqualAI
            E(S::Mod) + S ⇌ E(S, S::Mod)    :: EqualAI
        end
    end
    for (x, text) in ((ad1_ss, "EQ + B <--> EBQ"), (ad1_re, "EQ + B ⇌ EBQ"),
                      (ad11, "EQ + B <--> EBQ"), (ad3c, "E + Smod ⇌ ESmod"))
        @test EnzymeRates._guard(x) ==
              EnzymeRates._Verdict(NL, U, 0, 2, _testhelper_edge(x, text))
        @test _testhelper_exclusion(x, 2).failure === :inside
        @test _testhelper_epath_verdict(x; rng) === :no_path
    end
    # On AD1 and AD3c the message's steps that join the sides run from the binding's
    # unbound side to its bound side, and both fixes pass there too.
    for (x, path) in ((ad1_ss, ("E + Q ⇌ EQ", "E + B ⇌ EB", "EB + Q ⇌ EBQ")),
                      (ad3c, ("E + S ⇌ ES", "ES + Smod ⇌ ESSmod", "ESmod + S ⇌ ESSmod")))
        joins = join((EnzymeRates._step_text(x.catalytic, _testhelper_edge(x, t))
                      for t in path), "; ")
        @test endswith(EnzymeRates._guard_message(x, EnzymeRates._guard(x)),
                       "tag :OnlyA the steps that join its sides: $joins.")
    end
    ad1_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                    :: EqualAI
            (E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q))  :: EqualAI
            E(Q) + P ⇌ E(P, Q)                              :: EqualAI
            E + Q ⇌ E(Q)                                    :: EqualAI
            E + B ⇌ E(B)                                    :: EqualAI
            E(B) + Q ⇌ E(B, Q)                              :: EqualAI
        end
    end
    ad1_joins_tagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                    :: EqualAI
            (E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q))  :: OnlyA
            E(Q) + P ⇌ E(P, Q)                              :: EqualAI
            E + Q ⇌ E(Q)                                    :: OnlyA
            E + B ⇌ E(B)                                    :: OnlyA
            E(B) + Q ⇌ E(B, Q)                              :: OnlyA
        end
    end
    ad3c_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: S
        steps: begin
            E + S ⇌ E(S)                    :: EqualAI
            E(S) <--> E(P)                  :: OnlyA
            E + P ⇌ E(P)                    :: EqualAI
            E + S::Mod ⇌ E(S::Mod)          :: EqualAI
            E(S) + S::Mod ⇌ E(S, S::Mod)    :: EqualAI
            E(S::Mod) + S ⇌ E(S, S::Mod)    :: EqualAI
        end
    end
    ad3c_joins_tagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: S
        steps: begin
            E + S ⇌ E(S)                    :: OnlyA
            E(S) <--> E(P)                  :: OnlyA
            E + P ⇌ E(P)                    :: EqualAI
            E + S::Mod ⇌ E(S::Mod)          :: OnlyA
            E(S) + S::Mod ⇌ E(S, S::Mod)    :: OnlyA
            E(S::Mod) + S ⇌ E(S, S::Mod)    :: OnlyA
        end
    end
    # AD11's message names the steps AD1's does (E + Q ⇌ E(Q); E + B ⇌ E(B);
    # E(B) + Q ⇌ E(B, Q)); both its fixes pass every face, the oracles agreeing
    # (`nolimit_fixes.log`).
    ad11_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                            :: EqualAI
            (E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q),
             E(P) + B <--> E(B, P))                                 :: EqualAI
            E(Q) + P ⇌ E(P, Q)                                      :: EqualAI
            E + Q ⇌ E(Q)                                            :: EqualAI
            E + P ⇌ E(P)                                            :: EqualAI
            E + B ⇌ E(B)                                            :: EqualAI
            E(B) + Q ⇌ E(B, Q)                                      :: EqualAI
        end
    end
    ad11_joins_tagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                            :: EqualAI
            (E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q),
             E(P) + B <--> E(B, P))                                 :: OnlyA
            E(Q) + P ⇌ E(P, Q)                                      :: EqualAI
            E + Q ⇌ E(Q)                                            :: OnlyA
            E + P ⇌ E(P)                                            :: EqualAI
            E + B ⇌ E(B)                                            :: OnlyA
            E(B) + Q ⇌ E(B, Q)                                      :: OnlyA
        end
    end
    joins = join((EnzymeRates._step_text(ad11.catalytic, _testhelper_edge(ad11, t))
                  for t in ("E + Q ⇌ EQ", "E + B ⇌ EB", "EB + Q ⇌ EBQ")), "; ")
    @test endswith(EnzymeRates._guard_message(ad11, EnzymeRates._guard(ad11)),
                   "tag :OnlyA the steps that join its sides: $joins.")
    for y in (ad1_untagged, ad1_joins_tagged, ad3c_untagged, ad3c_joins_tagged,
              ad11_untagged, ad11_joins_tagged)
        @test EnzymeRates._guard(y) == EnzymeRates._WELL_POSED
    end
end

@testset "tagged bindings that close a directed cycle have no limit" begin
    rng = Random.MersenneTwister(37)
    U = EnzymeRates._LIMIT_UNTESTED
    # The ter cube of seven :OnlyA edges: the untagged E(A) + C, E(C) + B and E(B) + A
    # bindings join E(A) to E(A, B), so its :OnlyA binding of B stays inside one component.
    cube = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, C
        products: P
        steps: begin
            E + A ⇌ E(A)             :: OnlyA
            E + B ⇌ E(B)             :: OnlyA
            E + C ⇌ E(C)             :: OnlyA
            E(A) + B ⇌ E(A, B)       :: OnlyA
            E(A) + C ⇌ E(A, C)       :: EqualAI
            E(B) + A ⇌ E(A, B)       :: EqualAI
            E(B) + C ⇌ E(B, C)       :: EqualAI
            E(C) + A ⇌ E(A, C)       :: EqualAI
            E(C) + B ⇌ E(B, C)       :: EqualAI
            E(A, B) + C ⇌ E(A, B, C) :: OnlyA
            E(A, C) + B ⇌ E(A, B, C) :: OnlyA
            E(B, C) + A ⇌ E(A, B, C) :: OnlyA
            E(A, B, C) <--> E(P)     :: OnlyA
            E + P ⇌ E(P)             :: EqualAI
        end
    end
    @test _testhelper_exclusion(cube, 2).failure === :inside
    @test EnzymeRates._guard(cube) ==
          EnzymeRates._Verdict(EnzymeRates._NO_LIMIT, U, 0, 2,
                               _testhelper_edge(cube, "EA + B ⇌ EAB"))
    @test _testhelper_epath_verdict(cube; rng) === :no_path
    # Both fixes the message offers pass every face, the oracles agreeing
    # (`nolimit_fixes.log`): the B binding at E(A) untagged, and the five untagged bindings
    # that join its sides (E(A) + C, E(C) + A, E(C) + B, E(B) + C, E(B) + A) tagged.
    cube_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, C
        products: P
        steps: begin
            E + A ⇌ E(A)             :: OnlyA
            E + B ⇌ E(B)             :: OnlyA
            E + C ⇌ E(C)             :: OnlyA
            E(A) + B ⇌ E(A, B)       :: EqualAI
            E(A) + C ⇌ E(A, C)       :: EqualAI
            E(B) + A ⇌ E(A, B)       :: EqualAI
            E(B) + C ⇌ E(B, C)       :: EqualAI
            E(C) + A ⇌ E(A, C)       :: EqualAI
            E(C) + B ⇌ E(B, C)       :: EqualAI
            E(A, B) + C ⇌ E(A, B, C) :: OnlyA
            E(A, C) + B ⇌ E(A, B, C) :: OnlyA
            E(B, C) + A ⇌ E(A, B, C) :: OnlyA
            E(A, B, C) <--> E(P)     :: OnlyA
            E + P ⇌ E(P)             :: EqualAI
        end
    end
    cube_joins_tagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, C
        products: P
        steps: begin
            E + A ⇌ E(A)             :: OnlyA
            E + B ⇌ E(B)             :: OnlyA
            E + C ⇌ E(C)             :: OnlyA
            E(A) + B ⇌ E(A, B)       :: OnlyA
            E(A) + C ⇌ E(A, C)       :: OnlyA
            E(B) + A ⇌ E(A, B)       :: OnlyA
            E(B) + C ⇌ E(B, C)       :: OnlyA
            E(C) + A ⇌ E(A, C)       :: OnlyA
            E(C) + B ⇌ E(B, C)       :: OnlyA
            E(A, B) + C ⇌ E(A, B, C) :: OnlyA
            E(A, C) + B ⇌ E(A, B, C) :: OnlyA
            E(B, C) + A ⇌ E(A, B, C) :: OnlyA
            E(A, B, C) <--> E(P)     :: OnlyA
            E + P ⇌ E(P)             :: EqualAI
        end
    end
    @test EnzymeRates._guard(cube_untagged) == EnzymeRates._WELL_POSED
    @test EnzymeRates._guard(cube_joins_tagged) == EnzymeRates._WELL_POSED
    # Both bindings of the ordered bi-uni :OnlyA: the chemistry and the P binding join
    # E(A, B) back to E, so the arcs E → E(A) → {E, E(A, B), E(P)} close a cycle and the
    # two affinities cannot vanish together.
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A ⇌ E(A)            :: OnlyA
            E(A) + B ⇌ E(A, B)      :: OnlyA
            E(A, B) <--> E(P)       :: EqualAI
            E + P ⇌ E(P)            :: EqualAI
        end
    end
    m = am.catalytic
    a, b, chem, release = (_testhelper_edge(am, t) for t in
                           ("E + A ⇌ EA", "EA + B ⇌ EAB", "EAB <--> EP", "E + P ⇌ EP"))
    @test _testhelper_exclusion(am, 2).failure === :cycle
    @test EnzymeRates._exclusion_components(am, 2) == [1, 2, 1, 1]  # E, E(A), E(A, B), E(P)
    @test EnzymeRates._exclusion_components(am, 1) == [1, 1, 1, 1]
    v = EnzymeRates._guard(am)
    @test v == EnzymeRates._Verdict(EnzymeRates._ARC_CYCLE, U, 0, 2, b)
    @test _testhelper_epath_verdict(am; rng) === :no_path
    line(k) = EnzymeRates._step_text(m, k)
    @test EnzymeRates._guard_message(am, v) ==
          "AllostericMechanism: the :OnlyA bindings $(line(b)); $(line(a)) close a cycle " *
          "in the inactive conformation: steps outside :OnlyA groups join each one's " *
          "bound side to the next one's unbound side, so thermodynamics forbids their " *
          "affinities to vanish together. Tag them :EqualAI or :NonequalAI, or tag " *
          ":OnlyA the steps that join them: $(line(chem)); $(line(release))."
    # Both fixes pass the guard; untagging one binding alone leaves the other inside a
    # component, so the message does not offer it.
    untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A ⇌ E(A)            :: EqualAI
            E(A) + B ⇌ E(A, B)      :: EqualAI
            E(A, B) <--> E(P)       :: EqualAI
            E + P ⇌ E(P)            :: EqualAI
        end
    end
    joins_tagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A ⇌ E(A)            :: OnlyA
            E(A) + B ⇌ E(A, B)      :: OnlyA
            E(A, B) <--> E(P)       :: OnlyA
            E + P ⇌ E(P)            :: OnlyA
        end
    end
    one_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A ⇌ E(A)            :: OnlyA
            E(A) + B ⇌ E(A, B)      :: EqualAI
            E(A, B) <--> E(P)       :: EqualAI
            E + P ⇌ E(P)            :: EqualAI
        end
    end
    @test EnzymeRates._guard(untagged) == EnzymeRates._WELL_POSED
    @test EnzymeRates._guard(joins_tagged) == EnzymeRates._WELL_POSED
    @test EnzymeRates._guard(one_untagged).fault === EnzymeRates._NO_LIMIT
end

@testset "forms reached only through tagged chemistry are an island" begin
    rng = Random.MersenneTwister(41)
    IS, U = EnzymeRates._ISLAND, EnzymeRates._LIMIT_UNTESTED
    _testhelper_edge_between(x, a, b) = only(k for (k, e) in enumerate(x.catalytic.edges)
        if Set(EnzymeRates._form_name.(Ref(x.catalytic), (e.from, e.to))) == Set((a, b)))
    _testhelper_island_verdict(x) =
        _testhelper_epath_verdict(x; rng) in (:direction_dependent, :wrong_limit)
    # Both half-reactions of the ping-pong :OnlyA: the covalent forms join free E only
    # through chemistry. The message names half 1, which runs toward products into the
    # island (it is stored as the uptake of P), and its untagging leaves half 2 :OnlyA, a
    # well-posed variant (P4-29b).
    both = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)                                       :: EqualAI
            E(A) <--> E(; residual = A - P) + P                   :: OnlyA
            E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
            E(B; residual = A - P) <--> E + Q                     :: OnlyA
        end
    end
    half1 = _testhelper_edge(both, "E_res_plusA_minusP + P <--> EA")
    @test !EnzymeRates._runs_toward_products(both.catalytic, both.catalytic.edges[half1])
    v = EnzymeRates._guard(both)
    @test v == EnzymeRates._Verdict(IS, U, 0, 2, half1)
    @test !isempty(_testhelper_exclusion(both, 2).lost)
    @test _testhelper_island_verdict(both)
    line = EnzymeRates._step_text(both.catalytic, half1)
    @test EnzymeRates._guard_message(both, v) ==
          "AllostericMechanism: in the inactive conformation the forms " *
          "E(; residual = A - P), E(B; residual = A - P) join free E only through " *
          ":OnlyA steps that are not bindings toward them, such as $line, so their " *
          "weight depends on how fast those steps vanish. Tag the group of $line " *
          ":EqualAI or :NonequalAI."
    @test EnzymeRates._guard(_testhelper_tag_sets["pp_V_half2"]) == EnzymeRates._WELL_POSED
    # R3 and R4: a second free form behind a tagged isomerization; untagging it passes.
    r3 = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)        :: EqualAI
            E(S) <--> E(P)      :: EqualAI
            E + P ⇌ E(P)        :: EqualAI
            E <--> F            :: OnlyA
            F + I ⇌ F(I)        :: EqualAI
        end
    end
    r3_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)        :: EqualAI
            E(S) <--> E(P)      :: EqualAI
            E + P ⇌ E(P)        :: EqualAI
            E <--> F            :: EqualAI
            F + I ⇌ F(I)        :: EqualAI
        end
    end
    r4 = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: EqualAI
            E(S) <--> E(P)      :: NonequalAI
            E + P ⇌ E(P)        :: EqualAI
            F ⇌ E               :: OnlyA
        end
    end
    r4_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: EqualAI
            E(S) <--> E(P)      :: NonequalAI
            E + P ⇌ E(P)        :: EqualAI
            F ⇌ E               :: EqualAI
        end
    end
    for (x, fixed, forms) in ((r3, r3_untagged, "F, F(I)"), (r4, r4_untagged, "F"))
        iso = _testhelper_edge_between(x, "E", "F")
        v = EnzymeRates._guard(x)
        @test v == EnzymeRates._Verdict(IS, U, 0, 2, iso)
        @test _testhelper_island_verdict(x)
        msg = EnzymeRates._guard_message(x, v)
        line = EnzymeRates._step_text(x.catalytic, iso)
        @test occursin("the forms $forms join free E", msg)
        @test endswith(msg, "Tag the group of $line :EqualAI or :NonequalAI.")
        @test EnzymeRates._guard(fixed) == EnzymeRates._WELL_POSED
    end
    # AD2b and the ping-pong with B's binding and half 1 :OnlyA: a tagged binding borders
    # the island too, so untagging the named step leaves that binding inside a component,
    # and the message offers no fix.
    ad2b = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: EqualAI
            E(S) <--> E(P)      :: EqualAI
            F + P ⇌ E(P)        :: OnlyA
            F <--> E            :: OnlyA
        end
    end
    bbind = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)                                       :: EqualAI
            E(A) <--> E(; residual = A - P) + P                   :: OnlyA
            E(; residual = A - P) + B <--> E(B; residual = A - P) :: OnlyA
            E(B; residual = A - P) <--> E + Q                     :: EqualAI
        end
    end
    for (x, named) in ((ad2b, _testhelper_edge_between(ad2b, "E", "F")),
                       (bbind, _testhelper_edge(bbind, "E_res_plusA_minusP + P <--> EA")))
        v = EnzymeRates._guard(x)
        @test v == EnzymeRates._Verdict(IS, U, 0, 2, named)
        @test _testhelper_island_verdict(x)
        @test !occursin("Tag the group", EnzymeRates._guard_message(x, v))
    end
end

@testset "the admitted tag sets pass admission's interior" begin
    @testset "$name" for name in sort!(collect(keys(_testhelper_tag_sets)))
        @test EnzymeRates._admission_interior(_testhelper_tag_sets[name]) ==
              EnzymeRates._WELL_POSED
    end
end

@testset "a mixed :OnlyA group whose limit exists is admitted" begin
    # B's fused binding into E(P, Q) and its plain binding at E(Q) share one :OnlyA group;
    # in the inactive conformation no untagged step joins E(Q) to E(B, Q), so the plain
    # binding's affinity can vanish, and the limit is deletion plus pruning.
    am = _testhelper_tag_sets["mixed_group_limit"]
    @test _testhelper_exclusion(am, 2).failure === :none
    @test _testhelper_lp_existence(am, 2)
    @test EnzymeRates._admission_interior(am) == EnzymeRates._WELL_POSED
    @test EnzymeRates._guard(am) == EnzymeRates._WELL_POSED
    @test _testhelper_epath_verdict(am; rng = Random.MersenneTwister(43)) === :converges
end
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T4_red.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_admission.jl")'`
Expected: each new testset errors: ``UndefVarError: `_admission_interior` not defined in
`EnzymeRates` `` in the grid, admitted-tag-set and mixed-group testsets, ``UndefVarError:
`_NO_LIMIT` not defined`` in the no-limit and arc-cycle testsets and ``UndefVarError:
`_ISLAND` not defined`` in the island testset; Tasks 2-3's testsets still pass.

- [ ] **Step 3: Add the admission faults**

In `src/admission.jl`, replace the docstring and the `@enum _Fault` declaration of Task 3 with:

```julia
"""
What the guard finds wrong with a mechanism (spec section 8): nothing (`_NO_FAULT`); no
cycle turns the reaction over (`_NO_CYCLE`, class 1); a cycle of rapid-equilibrium steps
turns it over (`_RE_TURNOVER`, class 2); on a face (class 3), a second closed class
(`_NONUNIQUE`), a rapid-equilibrium segment with no lowest form (`_BOTTOMLESS`), or a law
the limit test rejects where procedure A stopped at a trap fed on the face (`_TRAP`), a
segment whose lowest forms split (`_SPLIT`), a law segment chain with several closed
classes (`_LAW_MULTI`) or an aggregated chain with several (`_MULTILEVEL`); in admission's
interior (class 4, spec section 7), a tagged binding whose two sides the conformation's
other steps join (`_NO_LIMIT`), tagged bindings that close a directed cycle
(`_ARC_CYCLE`), or forms free E's component does not reach along the tagged bindings
(`_ISLAND`).
"""
@enum _Fault::UInt8 begin
    _NO_FAULT
    _NO_CYCLE
    _RE_TURNOVER
    _NONUNIQUE
    _BOTTOMLESS
    _TRAP
    _SPLIT
    _LAW_MULTI
    _MULTILEVEL
    _NO_LIMIT
    _ARC_CYCLE
    _ISLAND
end
```

Directly after the declaration add:

```julia
"""The faults admission reports (spec section 7); `_admission_message` writes their text."""
const _ADMISSION_FAULTS = (_NO_LIMIT, _ARC_CYCLE, _ISLAND)
```

- [ ] **Step 4: Write the interior checks**

In `src/admission.jl`, insert before the section line `# ─── Guard messages` (Task 2), so
that it follows Task 3's limit test:

```julia
# ─── Admission: the interior (spec section 7, research section 1.5) ────

"""
Each form's component under the edges of `am`'s catalytic graph that conformation `c`
(1 active, 2 inactive) keeps (`_excludes`), every form included, labelled by its smallest
form index. Admission reads each plain binding `c` excludes as an arc from its unbound
side's component to its bound side's (`_arcs`).
"""
function _exclusion_components(am::AllostericMechanism, c::Int)
    m = am.catalytic
    parent = collect(eachindex(m.forms))
    top(i) = (while parent[i] != i; i = parent[i]; end; i)
    for (k, e) in enumerate(m.edges)
        _excludes(am, k, c) && continue
        a, b = top(e.from), top(e.to)
        a == b || (parent[max(a, b)] = min(a, b))
    end
    [top(f) for f in eachindex(m.forms)]
end

"""The plain bindings of `am`'s catalytic graph that conformation `c` excludes, in edge
order: admission's arcs, each stored as its uptake, unbound side first."""
_arcs(am::AllostericMechanism, c::Int) =
    [k for (k, e) in enumerate(am.catalytic.edges)
     if _excludes(am, k, c) && _is_plain_binding(am.catalytic, e)]

"""The first edge of `am`'s catalytic graph that conformation `c` excludes; 0 when it
excludes none."""
function _first_excluded(am::AllostericMechanism, c::Int)
    for k in eachindex(am.catalytic.edges)
        _excludes(am, k, c) && return k
    end
    0
end

"""
The arcs of one directed cycle that `arcs` close over the components `comp` (each arc an
edge of `m` read from its unbound side's component to its bound side's), in cycle order;
empty when they close none. Kahn's algorithm removes every component no cycle feeds; a
walk back along the remaining arcs then lands on a cycle.
"""
function _arc_cycle(m::Mechanism, comp::Vector{Int}, arcs::Vector{Int})
    tail(k) = comp[m.edges[k].from]
    head(k) = comp[m.edges[k].to]
    indeg = zeros(Int, length(comp))
    for k in arcs
        indeg[head(k)] += 1
    end
    todo = [r for r in unique(comp) if indeg[r] == 0]
    while !isempty(todo)
        u = pop!(todo)
        for k in arcs
            tail(k) == u && (indeg[head(k)] -= 1) == 0 && push!(todo, head(k))
        end
    end
    start = 0
    for k in arcs
        indeg[tail(k)] > 0 && indeg[head(k)] > 0 && (start = k; break)
    end
    start == 0 && return Int[]
    function back(x)
        for k in arcs
            head(k) == x && indeg[tail(k)] > 0 && return k
        end
        0
    end
    x = head(start)
    for _ in eachindex(comp)
        x = tail(back(x))
    end
    cycle = [back(x)]
    while tail(cycle[end]) != x
        push!(cycle, back(tail(cycle[end])))
    end
    reverse!(cycle)
end

"""The components (`comp` labels) that free E's reaches along `arcs`, edges of `m` read
from their unbound side's component to their bound side's."""
function _reached_components(m::Mechanism, comp::Vector{Int}, arcs::Vector{Int})
    reached = falses(length(comp))
    reached[comp[1]] = true
    grew = true
    while grew
        grew = false
        for k in arcs
            a, b = comp[m.edges[k].from], comp[m.edges[k].to]
            reached[a] && !reached[b] && (reached[b] = true; grew = true)
        end
    end
    reached
end

"""
Whether edge `e` of `m` runs, as stored, toward products: from the side that holds more
substrates, bound or free, then fewer products, the order canonical form gives a step that
exchanges several metabolites (`_orientation`). A step that takes up or gives off one
metabolite is stored as its uptake, so a half-reaction that releases a product,
`E(A) <--> E(; residual = A - P) + P`, runs toward products against its stored direction.
A step whose two sides tie, such as an isomerization of free forms, runs as stored.
"""
function _runs_toward_products(m::Mechanism, e::Edge)
    substrate_mask, product_mask = _role_masks(m)
    score(f::Int, free::Counts) =
        (total(m.forms[f].bound & substrate_mask) + total(free & substrate_mask),
         -total(m.forms[f].bound & product_mask) - total(free & product_mask))
    score(e.from, e.consumed) >= score(e.to, e.released)
end

"""
The edge an island's message names, 0 when free E's component reaches every component
(`reached`): of the edges conformation `c` excludes that join a reached component to an
unreached one, the first chemistry edge (`_is_chemistry`) that, run toward products
(`_runs_toward_products`), enters the unreached one; else the first chemistry edge; else
the first edge.
"""
function _island_edge(am::AllostericMechanism, c::Int, comp::Vector{Int},
                      reached::BitVector)
    m = am.catalytic
    best, rank = 0, 4
    for (k, e) in enumerate(m.edges)
        _excludes(am, k, c) || continue
        a, b = reached[comp[e.from]], reached[comp[e.to]]
        a == b && continue
        r = !_is_chemistry(m, e) ? 3 : a == _runs_toward_products(m, e) ? 1 : 2
        if r < rank
            best, rank = k, r
        end
    end
    best
end

"""
    _admission_interior(am) -> _Verdict

Admission's interior checks (spec section 7), for each conformation `c` whose tags exclude
an edge (`_excludes`), on the components of the forms its other edges join
(`_exclusion_components`) with its tagged plain bindings as arcs (`_arcs`):
- existence: no arc stays inside one component (`_NO_LIMIT`) and the arcs close no
  directed cycle (`_ARC_CYCLE`); equivalent to the exact cone test, since with every other
  constant fixed a direction's equilibrium-ratio exponents are differences of component
  potentials and only the arcs need a sign;
- uniqueness: free E's component reaches every component along the arcs, so every island
  lies behind a tagged plain binding on its bound side and the limit is deletion plus
  pruning (`_ISLAND`).
The verdict carries the conformation and the edge its message leads with (face 0).
"""
function _admission_interior(am::AllostericMechanism)
    m = am.catalytic
    for c in 1:2
        _first_excluded(am, c) == 0 && continue
        comp = _exclusion_components(am, c)
        arcs = _arcs(am, c)
        for k in arcs
            comp[m.edges[k].from] == comp[m.edges[k].to] &&
                return _Verdict(_NO_LIMIT, _LIMIT_UNTESTED, 0, c, k)
        end
        cycle = _arc_cycle(m, comp, arcs)
        isempty(cycle) || return _Verdict(_ARC_CYCLE, _LIMIT_UNTESTED, 0, c, cycle[1])
        k = _island_edge(am, c, comp, _reached_components(m, comp, arcs))
        k == 0 || return _Verdict(_ISLAND, _LIMIT_UNTESTED, 0, c, k)
    end
    _WELL_POSED
end

"""The edges of a shortest path from form `a` to form `b` over the edges of `am`'s
catalytic graph that conformation `c` keeps, in path order; empty when `a == b`. Callers
pass two forms those edges join."""
function _kept_path(am::AllostericMechanism, c::Int, a::Int, b::Int)
    m = am.catalytic
    via = zeros(Int, length(m.forms))
    seen = falses(length(m.forms))
    seen[a] = true
    queue = [a]
    head = 0
    while !seen[b]
        head += 1
        u = queue[head]
        for (k, e) in enumerate(m.edges)
            _excludes(am, k, c) && continue
            v = e.from == u ? e.to : e.to == u ? e.from : 0
            (v == 0 || seen[v]) && continue
            seen[v] = true
            via[v] = k
            push!(queue, v)
        end
    end
    path = Int[]
    f = b
    while f != a
        k = via[f]
        pushfirst!(path, k)
        f = m.edges[k].from == f ? m.edges[k].to : m.edges[k].from
    end
    path
end

"""
The message of admission's verdict `v` on `am`: the cause, naming each step by its DSL
line (`_step_text`), and a fix where a test shows its variant passes the guard (P4-22).
A fix that tags steps is offered for the inactive conformation only, since a step tagged
`:OnlyI` must bind modifiers alone.
"""
function _admission_message(am::AllostericMechanism, v::_Verdict)
    m = am.catalytic
    c, k = v.conformation, v.edge
    conf = c == 1 ? "active" : "inactive"
    tag = c == 1 ? ":OnlyI" : ":OnlyA"
    line(j) = _step_text(m, j)
    lines(js) = join((_step_text(m, j) for j in js), "; ")
    if v.fault === _NO_LIMIT
        msg = "AllostericMechanism: in the $conf conformation, steps outside $tag " *
              "groups still join the two sides of the $tag binding $(line(k)), so " *
              "thermodynamics keeps its affinity finite and the binding has no limit " *
              "there."
        c == 2 || return msg
        e = m.edges[k]
        return msg * " Tag its group :EqualAI or :NonequalAI, or tag :OnlyA the steps " *
               "that join its sides: $(lines(_kept_path(am, c, e.from, e.to)))."
    elseif v.fault === _ARC_CYCLE
        cycle = _arc_cycle(m, _exclusion_components(am, c), _arcs(am, c))
        msg = "AllostericMechanism: the $tag bindings $(lines(cycle)) close a cycle in " *
              "the $conf conformation: steps outside $tag groups join each one's bound " *
              "side to the next one's unbound side, so thermodynamics forbids their " *
              "affinities to vanish together."
        c == 2 || return msg
        joins = Int[]
        for (i, j) in enumerate(cycle)
            next = cycle[mod1(i + 1, length(cycle))]
            append!(joins, _kept_path(am, c, m.edges[j].to, m.edges[next].from))
        end
        return msg * " Tag them :EqualAI or :NonequalAI, or tag :OnlyA the steps that " *
               "join them: $(lines(unique!(joins)))."
    elseif v.fault === _ISLAND
        comp = _exclusion_components(am, c)
        reached = _reached_components(m, comp, _arcs(am, c))
        e = m.edges[k]
        lost = reached[comp[e.from]] ? comp[e.to] : comp[e.from]
        # Untagging the named step joins every lost form to free E's component only when
        # one component is lost and the tagged steps it borders are all chemistry.
        fixable = all(f -> reached[comp[f]] || comp[f] == lost, eachindex(m.forms))
        for (j, d) in enumerate(m.edges)
            _excludes(am, j, c) && reached[comp[d.from]] != reached[comp[d.to]] &&
                !_is_chemistry(m, d) && (fixable = false)
        end
        msg = "AllostericMechanism: in the $conf conformation the forms " *
              join([_form_text(m, f) for f in eachindex(m.forms) if comp[f] == lost],
                   ", ") *
              " join free E only through $tag steps that are not bindings toward them, " *
              "such as $(line(k)), so their weight depends on how fast those steps vanish."
        return c == 2 && fixable ?
               msg * " Tag the group of $(line(k)) :EqualAI or :NonequalAI." : msg
    else
        error("_admission_message: $(v.fault) is no admission fault")
    end
end
```

- [ ] **Step 5: Wire the interior into the guard and its message**

In `src/admission.jl`, replace the docstring of `_guard` and its allosteric method (Task 3's)
with the following; the plain method `_guard(m::Mechanism)` between them stays as it is:

```julia
"""
    _guard(m) -> _Verdict

The well-posedness guard (spec section 8) on `m`, plain or allosteric: class 1 (a cycle
turns the reaction over, `_has_turnover_cycle`); class 2 (no cycle of fast transitions
does, `_re_turnover_transitions`, on the ε-network, where an excluded rapid-equilibrium
binding stays fast while its affinity vanishes); for an allosteric mechanism, admission's
interior (`_admission_interior`); then class 3 on every face (`_check_face`, on the law
network), the first fault winning. `_WELL_POSED` when every check passes. It reads the
reaction graph through masks and builds nothing.
"""
```

```julia
function _guard(am::AllostericMechanism)
    m = am.catalytic
    _has_turnover_cycle(m, trues(length(m.edges))) ||
        return _Verdict(_NO_CYCLE, _LIMIT_UNTESTED, 0, 0, 0)
    v = _re_turnover_verdict(_eps_network(am))
    v === _WELL_POSED || return v
    v = _admission_interior(am)
    v === _WELL_POSED || return v
    _faces_verdict(am, _law_network(am))
end
```

In `_guard_message(mech::_AnyMech, v::_Verdict)` (Tasks 2-3), directly after the line
`    v.fault === _RE_TURNOVER && return _re_turnover_message(who, mech)`, insert:

```julia
    mech isa AllostericMechanism && v.fault in _ADMISSION_FAULTS &&
        return _admission_message(mech, v)
```

- [ ] **Step 6: Run the tests to verify they pass**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T4_green.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_admission.jl")'`
Expected: every testset of test_admission.jl passes, no `Error`, `Fail` or warning line. Compare
the file's Test Summary time with `P4_T3_green.log`: at most 1 s more.

- [ ] **Step 7: Report the interior's cost**

Create `$S/phase45/probes/P4_T4_cost.jl`:

```julia
# ABOUTME: Times admission's interior checks per mechanism over the uni-uni grid, for the
# ABOUTME: cost report of Task 4; run from the guard worktree under TestEnv.
using Test, EnzymeRates, LinearAlgebra, Random
include(joinpath(pwd(), "test", "mechanism_definitions_for_test_enzyme_derivation.jl"))
grid = _testhelper_uniuni_grid()
foreach(EnzymeRates._admission_interior, grid)
t = minimum(@elapsed(foreach(EnzymeRates._admission_interior, grid)) for _ in 1:50)
println("admission interior: ", round(1e6 * t / length(grid); digits = 2),
        " µs per mechanism over ", length(grid))
```

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T4_cost.log 900 --project=. -e 'using TestEnv; TestEnv.activate(); include("'$S'/phase45/probes/P4_T4_cost.jl")'`
Expected: one line `admission interior: <t> µs per mechanism over 108`, `<t>` of the order of the
prototype's 2.5 µs; the commit body reports it.

- [ ] **Step 8: Run the full suite**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T4_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`
Expected: `exit 0`, `Testing EnzymeRates tests passed`, no `Error`/`Fail`/warning lines, Aqua and
JET green; note the `maxrss_kb` line and explain any rise over the branch point's figure in the
commit body.

- [ ] **Step 9: Commit**

Write the measured figures in place of `<t>` (Step 7) and `<kb>` (Step 8).

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
git add src/admission.jl test/test_admission.jl
git commit -m "$(cat <<'EOF'
Admit tag sets by the existence and island tests

Admission's interior runs after class 2: a tagged plain binding inside one
component or arcs closing a cycle has no limit, and forms that join free E
only through tagged chemistry are an island. The graph test equals the exact
LP on the 108 uni-uni grid mechanisms. The arc-cycle message offers the two
fixes a test confirms; research's "untag one of them" fails on the ordered
bi-uni. The island message names the tagged chemistry that runs toward
products into the island, since a half-reaction that releases a product is
stored as its uptake. Interior cost: <t> µs per mechanism. Full-suite maxrss:
<kb> KB.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
)"
```

---

### Task 5: Admission's face condition B v3

The per-face half of admission (RES §1.5 ADMIT_FACE): the ε-network's closed classes, stranded
basins, R1 or T0+ per excluding conformation, and joint R1 (P4-8, P4-8b). Wired into `_guard(am)`
after `_check_face` on each face.

**Files:**
- Modify: `src/admission.jl`
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (section "Phase-4
  oracles": `_testhelper_class3_faults`, `_testhelper_admission_faults`,
  `_testhelper_admission_agrees`, `_testhelper_admission_agreement`)
- Modify: `ROADMAP.md` (Ideas: the conservative R-B1 and R-B2 rejections, P4-9)
- Modify: `docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md` (section 7
  "Admission" item 2, section 8 "Phase 4's first task", the phase-4 row of section 15,
  section 17 row 2: every passage that asks for an exact face condition; P4-9a)
- Test: `test/test_admission.jl`

**Interfaces:**
- Consumes (Task 2): `_Network` (fields `nforms, ncopies, present, from, to, edge,
  conformation, forward, fast, vanishing, change, takes`; `from`/`to` are node indices
  `(c - 1) * nforms + f`, transitions `2j - 1` and `2j` a pair), `_node(net, f, c)`,
  `_eps_network(am, re)`, `_law_network(am, re)`, `_faces(m)::Vector{UInt64}`,
  `_face_names(m, face)`, `_turnover_name(m)::Int`, `_re_turnover_verdict`, `_guard(am)`;
  (Task 3) `_closed_classes(net::_Network, keep::BitVector)::Vector{Vector{Int}}`,
  `_bottom_sccs(next, nodes)` (rewritten here),
  `_copy_of(net, u)`, `_form_of(net, u)`, `_check_face(mech, net, faces, i)`,
  `_testhelper_agrees(fault, verdict)`, `_testhelper_face_faults(m)`; (Task 4)
  `_ADMISSION_FAULTS`, `_admission_interior`, `_admission_message`, `_testhelper_tag_sets`;
  (Task 1) `_testhelper_faces(m)`, `_testhelper_respeed_verdicts(m; rng, npoints, law)`,
  `_testhelper_epath_verdict(m; rng, face, corners, precision)`; (phase 3)
  `_flags`, `_group_members`, `_form_text(m, f)` and `_step_text(m, k)` (names.jl).
- Produces:
  ```julia
  # _Fault gains, after _ISLAND: _EPS_NONUNIQUE _PUMPED _BRIDGE
  _reach(net::_Network, keep::BitVector, start::BitVector, backward::Bool)::BitVector
  _face_sets(net::_Network, law_nodes::BitVector, face::UInt64)
      # -> (; alive::BitVector, fed::BitVector, n0::BitVector, stranded::BitVector)
  _copy_transitions(net::_Network, alive::BitVector, nodes::BitVector, c::Int)
      # -> (; van::BitVector, held::BitVector)
  _min_arborescence(n::Int, arcs::Vector{NTuple{3, Int}}, root::Int)::Int   # -1 when none
  _fewest_slow(net::_Network, nodes::BitVector, keep::BitVector)::Int       # MINSLOW; -1 none
  _strong_components(next::Vector{Vector{Int}}, nodes::BitVector)::Vector{Int}  # Tarjan; 0 unvisited
  _potential_on_sccs(net::_Network, keep::BitVector, x::Int)::Bool          # T0+_c
  _bridges_jointly(net::_Network, fed::BitVector, base::BitVector, vanishing::BitVector)::Bool
  _bridging_edges(net::_Network, nodes::BitVector, base::BitVector, van::BitVector)::Vector{Int}
  _unfed_class(closed::Vector{Vector{Int}}, fed::BitVector)::Vector{Int}
  _face_phrase(m::Mechanism, face::UInt64)::String
  _face_condition(net::_Network, law_nodes::BitVector, face::UInt64, x::Int)::_Verdict  # face 0
  _admission_face(am::AllostericMechanism, eps::_Network, law::_Network,
                  faces::Vector{UInt64}, i::Int)::_Verdict
  _admission(am::AllostericMechanism)::_Verdict   # interior, then every face; no class 3
  # fixtures file:
  const _testhelper_class3_faults::NTuple{6, EnzymeRates._Fault}
  _testhelper_admission_faults(am)::Vector{EnzymeRates._Fault}       # class 3, then B v3
  _testhelper_admission_agrees(fault::EnzymeRates._Fault, verdict::Symbol)::Bool
  _testhelper_admission_agreement(am; rng)::Vector{Bool}             # one per face
  ```
  `_guard(am)` runs `_admission_face` on each face after `_check_face` accepts it.

**Gate:** each `:OnlyA` allosteric row equals its oracle (`_testhelper_admission_agreement`):
the 39 tag sets outside `_testhelper_face_rejected` pass every face, R1, A5, AD2a, AD2d, AD3a,
AD3b, AD3d, AD4, AD4b, AD5, AD7a, AD7c, AD8, AD9, AD10, AD13 ×10, AD14, AD14b, AD14c, AD15,
R2P_sealed, B_uphill, TC1, pp_V_half1, pp_V_half2, pp_K_Abind_half1, PP1, PP2, terbi
d75d2415361c4dae (T0+ at products = 0 and P = 0) and mixed_group_limit among them; R2 products = 0
`_PUMPED` by the face condition alone (ε-path `:direction_dependent`); BS2, BS3, BS4 products = 0
`_PUMPED`; ADbridge, bridge_drain `_BRIDGE` at the interior, products = 0 and substrates = 0,
class 3 on 3 and 2 modifier faces; R2funnel, R2funnel_SS products = 0 `_BRIDGE`;
ADbridge_noturnover admitted at its first three faces by T0+; S1 class 3 `_NONUNIQUE` at
substrates = 0, A = 0, B = 0; AD7b `_NONUNIQUE` on its three modifier faces. Known conservative
(P4-9a), checked with corners: bridge_deadend `_BRIDGE` ×3 (R-B1), BS1 `_PUMPED` at products = 0
and substrates = 0 with I zero (R-B2), fz101_10232 `_PUMPED` at substrates = 0 (R-B2), the oracle
`:converges` on each; fz101_10232 at products = 0 `_PUMPED` with a corner direction failing. The
hand-built joint-R1 network is `_BRIDGE` (P4-8b; Task 10 adds a mechanism that realizes the
pattern). test_admission.jl grows by at most 3 s
without the corner rows and at most 8 s with them; the face condition's cost is reported
(prototype 149 µs per mechanism over all faces).

- [ ] **Step 1: Write the oracle-side helpers**

Append to the section "Phase-4 oracles" of
`test/mechanism_definitions_for_test_enzyme_derivation.jl`:

```julia
"""The faults class 3 reports on a face (`_check_face`): the faces the RE-speed oracle
judges."""
const _testhelper_class3_faults = (EnzymeRates._NONUNIQUE, EnzymeRates._BOTTOMLESS,
    EnzymeRates._TRAP, EnzymeRates._SPLIT, EnzymeRates._LAW_MULTI, EnzymeRates._MULTILEVEL)

"""
The guard's fault on each face of the allosteric mechanism `am`, in the order of `_faces`
(and `_testhelper_faces`), after its interior checks pass: class 3's where `_check_face`
rejects the face, else that of admission's face condition (`_admission_face`).
"""
function _testhelper_admission_faults(am)
    eps, law = EnzymeRates._eps_network(am), EnzymeRates._law_network(am)
    faces = EnzymeRates._faces(am.catalytic)
    map(eachindex(faces)) do i
        v = EnzymeRates._check_face(am, law, faces, i)
        v.fault === EnzymeRates._NO_FAULT ?
            EnzymeRates._admission_face(am, eps, law, faces, i).fault : v.fault
    end
end

"""Whether admission's `fault` on a face agrees with the ε-path oracle's `verdict` there:
an accepted face converges, a rejected one does not."""
_testhelper_admission_agrees(fault::EnzymeRates._Fault, verdict::Symbol) =
    fault === EnzymeRates._NO_FAULT ? verdict === :converges :
    verdict in (:no_path, :del_nonunique, :nonunique, :direction_dependent, :wrong_limit)

"""
Whether the guard's fault on each face of the allosteric mechanism `am`
(`_testhelper_admission_faults`) agrees with the oracle that judges it: a class-3 fault
with the RE-speed oracle (`_testhelper_agrees`), an admission fault with the ε-path oracle
(`_testhelper_admission_agrees`), and an accepted face with both. The RE-speed oracle
compares with its own deletion limit (`law = :limit`): the package cannot derive the law
of a mechanism the guard rejects. One `Bool` per face.
"""
function _testhelper_admission_agreement(am; rng)
    faults = _testhelper_admission_faults(am)
    faces = _testhelper_faces(am)
    respeed = _testhelper_respeed_verdicts(am; rng, law = :limit)
    map(eachindex(faults)) do i
        f = faults[i]
        f in _testhelper_class3_faults && return _testhelper_agrees(f, respeed[i])
        epath = _testhelper_admission_agrees(f, _testhelper_epath_verdict(am; rng,
                                                                          face = faces[i]))
        f === EnzymeRates._NO_FAULT ? epath && respeed[i] === :agree : epath
    end
end
```

- [ ] **Step 2: Write the failing tests**

Append to `test/test_admission.jl`:

```julia
"""The tag sets of `_testhelper_tag_sets` that a face rejects; every other one passes every
face."""
const _testhelper_face_rejected = Set(["R2", "BS1", "BS2", "BS3", "BS4", "ADbridge",
    "ADbridge_noturnover", "bridge_drain", "bridge_deadend", "R2funnel", "R2funnel_SS",
    "S1", "AD7b", "fz101_10232"])

@testset "R2's limit at zero product depends on how fast its isomerization vanishes" begin
    rng = Random.MersenneTwister(43)
    NO, PU = EnzymeRates._NO_FAULT, EnzymeRates._PUMPED
    # Classes 1-3 and admission's interior accept R2, so the face condition alone rejects
    # it (spec section 7): at P zero, F's exits vanish (F <--> E) or take up P (F + P ⇌
    # E(P)), and E(P) feeds it by a release whose reverse is dead.
    am = _testhelper_tag_sets["R2"]
    m = am.catalytic
    @test _testhelper_face_faults(am) == fill(NO, 3)
    @test EnzymeRates._admission_interior(am) == EnzymeRates._WELL_POSED
    @test _testhelper_admission_faults(am) == [NO, PU, NO]
    @test _testhelper_epath_verdict(am; rng, face = _testhelper_faces(am)[2]) ===
          :direction_dependent
    @test all(_testhelper_admission_agreement(am; rng))
    k, iso = _testhelper_edge(am, "F + P ⇌ EP"), _testhelper_edge(am, "F <--> E")
    v = EnzymeRates._guard(am)
    @test v == EnzymeRates._Verdict(PU, EnzymeRates._LIMIT_UNTESTED, 2, 2, k)
    line(j) = EnzymeRates._step_text(m, j)
    @test EnzymeRates._guard_message(am, v) ==
          "AllostericMechanism: at zero P, $(line(k)) feeds the forms F of the inactive " *
          "conformation, which leave only through the vanishing steps $(line(iso)) while " *
          "nothing runs back the way it came, so their weight depends on how fast those " *
          "steps vanish."
end

@testset "admitted tag sets converge on every face" begin
    rng = Random.MersenneTwister(47)
    names = sort!([n for n in keys(_testhelper_tag_sets) if n ∉ _testhelper_face_rejected])
    @test length(names) == 39
    @testset "$name" for name in names
        am = _testhelper_tag_sets[name]
        @test EnzymeRates._guard(am) == EnzymeRates._WELL_POSED
        @test all(==(EnzymeRates._NO_FAULT), _testhelper_admission_faults(am))
        @test all(_testhelper_admission_agreement(am; rng))
    end
end

@testset "a sealed basin drains, and a pumped one does not" begin
    rng = Random.MersenneTwister(53)
    NO, PU = EnzymeRates._NO_FAULT, EnzymeRates._PUMPED
    function _testhelper_stranded(am, i)
        face = EnzymeRates._faces(am.catalytic)[i]
        s = EnzymeRates._face_sets(EnzymeRates._eps_network(am),
                                   EnzymeRates._law_network(am).present, face)
        any(s.stranded)
    end
    # A stranded basin entered only reversibly from the law network is sealed: its weight
    # vanishes with the tagged constants, and the face passes. R2P_sealed at P zero,
    # B_uphill and BS4 with X or I zero.
    for (name, sealed) in (("R2P_sealed", [2]), ("B_uphill", [4, 5, 6]), ("BS4", [4, 5, 6]))
        am = _testhelper_tag_sets[name]
        @test all(i -> _testhelper_stranded(am, i), sealed)
        @test all(==(NO), _testhelper_admission_faults(am)[sealed])
    end
    # A basin fed by a step whose reverse is dead is pumped: BS2, BS3 and BS4 at P zero.
    @test _testhelper_admission_faults(_testhelper_tag_sets["BS2"]) == [NO, PU, NO]
    @test _testhelper_admission_faults(_testhelper_tag_sets["BS3"]) == [NO, PU, NO]
    @test _testhelper_admission_faults(_testhelper_tag_sets["BS4"]) ==
          [NO, PU, NO, NO, NO, NO]
    for name in ("R2P_sealed", "B_uphill", "BS2", "BS3", "BS4")
        @test all(_testhelper_admission_agreement(_testhelper_tag_sets[name]; rng))
    end
end

@testset "a tagged binding that bridges two basins is rejected" begin
    rng = Random.MersenneTwister(59)
    NO, BR = EnzymeRates._NO_FAULT, EnzymeRates._BRIDGE
    SP, TR = EnzymeRates._SPLIT, EnzymeRates._TRAP
    # ADbridge: the :OnlyA bindings of I join E and E(S) of the inactive conformation, which
    # turns over through its steady-state S binding. With them the fed forms' cheapest
    # in-tree needs one slow step fewer (R1 fails), and the turnover is no potential (T0+
    # fails). With I zero the bridge is dead, and class 3 finds the inhibitor square split.
    am = _testhelper_tag_sets["ADbridge"]
    m = am.catalytic
    @test _testhelper_admission_faults(am) == [BR, BR, BR, SP, SP, TR]
    @test all(_testhelper_admission_agreement(am; rng))
    k = _testhelper_edge(am, "ES + I ⇌ EIS")
    v = EnzymeRates._guard(am)
    @test v == EnzymeRates._Verdict(BR, EnzymeRates._LIMIT_UNTESTED, 1, 2, k)
    line = EnzymeRates._step_text(m, k)
    @test EnzymeRates._guard_message(am, v) ==
          "AllostericMechanism: the :OnlyA step $line bridges forms of the inactive " *
          "conformation that its other steps keep apart, so the limit depends on how " *
          "fast the tagged constants vanish. Make the group of $line steady-state."
    fixed = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S <--> E(S)                         :: EqualAI
            E(S) <--> E + P                         :: EqualAI
            (E + I <--> E(I), E(S) + I <--> E(I, S)) :: OnlyA
            E(I) + S ⇌ E(I, S)                      :: EqualAI
        end
    end
    @test EnzymeRates._guard(fixed) == EnzymeRates._WELL_POSED
    # The same bridge with the chemistry :OnlyA turns nothing over in the inactive
    # conformation: R1 still fails, but the first substrate's uptake is a potential on
    # every strongly connected part, so T0+ admits the first three faces.
    noturn = _testhelper_tag_sets["ADbridge_noturnover"]
    @test _testhelper_admission_faults(noturn) == [NO, NO, NO, SP, SP, TR]
    @test all(_testhelper_admission_agreement(noturn; rng))
    for (x, potential) in ((am, false), (noturn, true))
        eps, law = EnzymeRates._eps_network(x), EnzymeRates._law_network(x)
        row = EnzymeRates._turnover_name(x.catalytic)
        for face in EnzymeRates._faces(x.catalytic)[1:3]
            s = EnzymeRates._face_sets(eps, law.present, face)
            nodes = s.fed .& .!s.stranded
            (; van, held) = EnzymeRates._copy_transitions(eps, s.alive, nodes, 2)
            @test EnzymeRates._bridges_jointly(eps, nodes, s.n0, van)
            @test EnzymeRates._potential_on_sccs(eps, held, row) == potential
        end
    end
    # bridge_drain: :OnlyA bindings of I let the inactive F return to E.
    @test _testhelper_admission_faults(_testhelper_tag_sets["bridge_drain"]) ==
          [BR, BR, BR, SP, NO, SP]
    @test all(_testhelper_admission_agreement(_testhelper_tag_sets["bridge_drain"]; rng))
    # R2funnel: at P zero the trap {F, F(S)} drains only through a slow step out of a fast
    # transient; its bridge stays with the S binding steady-state (R2funnel_SS), so the
    # message offers no fix.
    for name in ("R2funnel", "R2funnel_SS")
        x = _testhelper_tag_sets[name]
        @test _testhelper_admission_faults(x) == [NO, BR, NO]
        @test all(_testhelper_admission_agreement(x; rng))
        v = EnzymeRates._guard(x)
        @test (v.fault, v.face, v.conformation) === (BR, 2, 2)
        @test !occursin("Make the group", EnzymeRates._guard_message(x, v))
    end
end

@testset "the census ter-bi mechanism passes by a potential on every component" begin
    rng = Random.MersenneTwister(61)
    am = _testhelper_tag_sets["terbi_d75d2415361c4dae"]
    m = am.catalytic
    @test EnzymeRates._guard(am) == EnzymeRates._WELL_POSED
    @test all(_testhelper_admission_agreement(am; rng))
    # With P zero its :OnlyA bindings shorten the inactive conformation's in-tree (R1
    # fails), yet the first substrate's uptake is a potential on each strongly connected
    # part of its kept and fast transitions (T0+): the census's case for T0+.
    eps, law = EnzymeRates._eps_network(am), EnzymeRates._law_network(am)
    row = EnzymeRates._turnover_name(m)
    faces = EnzymeRates._faces(m)
    by_potential = map(faces) do face
        s = EnzymeRates._face_sets(eps, law.present, face)
        nodes = s.fed .& .!s.stranded
        (; van, held) = EnzymeRates._copy_transitions(eps, s.alive, nodes, 2)
        EnzymeRates._bridges_jointly(eps, nodes, s.n0, van) &&
            EnzymeRates._potential_on_sccs(eps, held, row)
    end
    @test [Set(EnzymeRates._face_names(m, f)) for f in faces[by_potential]] ==
          [Set([:P, :Q]), Set([:P])]
end

@testset "conservative rejections the oracle admits stay rejected" begin
    rng = Random.MersenneTwister(67)
    NO, BR, PU = EnzymeRates._NO_FAULT, EnzymeRates._BRIDGE, EnzymeRates._PUMPED
    _testhelper_corner_verdict(am, i) =
        _testhelper_epath_verdict(am; rng, face = _testhelper_faces(am)[i], corners = true)
    # R-B1, an equilibrium bridge (P4-9a): the bridge joins E and the dead end E(J), which
    # deletion already holds at its equilibrium ratio, so the limit is deletion's, but no
    # potential proves it. The steady-state fix passes and the message offers it.
    am = _testhelper_tag_sets["bridge_deadend"]
    @test _testhelper_admission_faults(am) == [BR, BR, BR, NO, NO, NO]
    @test all(i -> _testhelper_corner_verdict(am, i) === :converges, 1:3)
    @test all(_testhelper_admission_agreement(am; rng)[4:6])
    k = _testhelper_edge(am, "EJ + I ⇌ EIJ")
    v = EnzymeRates._guard(am)
    @test v == EnzymeRates._Verdict(BR, EnzymeRates._LIMIT_UNTESTED, 1, 2, k)
    @test endswith(EnzymeRates._guard_message(am, v),
                   "Make the group of $(EnzymeRates._step_text(am.catalytic, k)) " *
                   "steady-state.")
    fixed = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S <--> E(S)         :: EqualAI
            E(S) <--> E + P         :: EqualAI
            E + J <--> E(J)         :: EqualAI
            E + I ⇌ E(I)            :: OnlyA
            E(J) + I <--> E(I, J)   :: OnlyA
            E(I) + J ⇌ E(I, J)      :: EqualAI
        end
    end
    @test EnzymeRates._guard(fixed) == EnzymeRates._WELL_POSED
    # R-B2, a stranded basin fed only through vanishing sources: BS1 with I zero and P or
    # S zero, fz101_10232 with S zero. Every corner direction converges.
    bs1 = _testhelper_tag_sets["BS1"]
    @test _testhelper_admission_faults(bs1) == [NO, NO, NO, NO, PU, PU]
    @test all(i -> _testhelper_corner_verdict(bs1, i) === :converges, 5:6)
    @test all(_testhelper_admission_agreement(bs1; rng)[1:4])
    fz = _testhelper_tag_sets["fz101_10232"]
    @test _testhelper_admission_faults(fz) == [NO, PU, PU, NO, PU, PU]
    @test _testhelper_corner_verdict(fz, 3) === :converges
    # At P zero the rejection is right: F's only alive exit, the blocked F <--> E, can
    # vanish faster than its inflow, which one corner direction shows (VI F3).
    @test _testhelper_corner_verdict(fz, 2) !== :converges
    @test all(_testhelper_admission_agreement(fz; rng)[[1, 4]])
end

@testset "S1 and AD7b pass admission's interior and fail class 3" begin
    rng = Random.MersenneTwister(71)
    NO, NU = EnzymeRates._NO_FAULT, EnzymeRates._NONUNIQUE
    # S1's stray free form F binds B at steady state: with A or B zero the law network has
    # a second closed class. The ε-network's closed-class count sees the same split, so
    # admission alone rejects those faces too; the guard reports class 3, which runs first.
    s1 = _testhelper_tag_sets["S1"]
    @test EnzymeRates._admission_interior(s1) == EnzymeRates._WELL_POSED
    @test _testhelper_admission_faults(s1) == [NO, NO, NU, NU, NU]
    @test EnzymeRates._guard(s1).fault === NU
    @test all(_testhelper_admission_agreement(s1; rng))
    v = EnzymeRates._admission(s1)
    @test (v.fault, v.face) === (EnzymeRates._EPS_NONUNIQUE, 3)
    @test startswith(EnzymeRates._guard_message(s1, v),
                     "AllostericMechanism: at zero A and B, the forms ")
    # AD7b: with I zero no spanning in-tree reaches the fed forms, so the face condition
    # passes and class 3 decides.
    ad7b = _testhelper_tag_sets["AD7b"]
    @test EnzymeRates._admission(ad7b) == EnzymeRates._WELL_POSED
    @test _testhelper_admission_faults(ad7b) == [NO, NO, NO, NU, NU, NU]
    @test all(_testhelper_admission_agreement(ad7b; rng))
end

@testset "joint R1 catches a tree that only both conformations' bridges shorten" begin
    # VC X4's in-tree costs, written in the guard's own type (P4-8b): the smallest network
    # that shows the pattern, for the in-tree arithmetic.
    # Copy A holds E and X (nodes 1, 2), copy I holds E, Y and Z (nodes 4, 5, 6), and the
    # flip joins the two copies of E. Alone, copy A's vanishing E → X and copy I's
    # vanishing Y → E each leave the cheapest spanning in-tree at one slow transition;
    # together they make it free, rooted at X. Transitions in order: E→X slow, X→E slow,
    # E→X fast and vanishing (copy A); E→Y fast, Y→Z slow, Z→Y fast, Z→E slow, Y→E fast and
    # vanishing (copy I); the flip. Only the transitions these costs need are present, so
    # most lack the reverse partner a mechanism's would have; the face condition reads the
    # partners only on a stranded basin, which this network has none of.
    net = EnzymeRates._Network(3, 2, BitVector([1, 1, 0, 1, 1, 1]),
        [1, 2, 1, 4, 5, 6, 6, 5, 1, 4],                  # from
        [2, 1, 2, 5, 6, 5, 4, 4, 4, 1],                  # to
        [1, 1, 2, 3, 4, 5, 6, 7, 0, 0],                  # edge
        [1, 1, 1, 2, 2, 2, 2, 2, 0, 0],                  # conformation
        BitVector([1, 0, 1, 1, 1, 1, 1, 1, 1, 0]),       # forward
        BitVector([0, 0, 1, 1, 0, 1, 0, 1, 1, 1]),       # fast
        BitVector([0, 0, 1, 0, 0, 0, 0, 1, 0, 0]),       # vanishing
        zeros(Int, 1, 10), zeros(UInt64, 10))
    s = EnzymeRates._face_sets(net, net.present, UInt64(0))
    nodes = s.fed .& .!s.stranded
    @test nodes == net.present
    @test EnzymeRates._fewest_slow(net, nodes, s.n0) == 1
    for c in 1:2
        (; van) = EnzymeRates._copy_transitions(net, s.alive, nodes, c)
        @test count(van) == 1
        @test !EnzymeRates._bridges_jointly(net, nodes, s.n0, van)
    end
    @test EnzymeRates._fewest_slow(net, nodes, s.alive) == 0
    @test EnzymeRates._face_condition(net, net.present, UInt64(0), 1) ==
          EnzymeRates._Verdict(EnzymeRates._BRIDGE, EnzymeRates._LIMIT_UNTESTED, 0, 0, 0)
    # The in-tree test reports -1 when no in-tree spans the nodes: nodes 2 and 6 share no
    # transition.
    @test EnzymeRates._fewest_slow(net, BitVector([0, 1, 0, 0, 0, 1]), trues(10)) == -1
    # T0+'s potential: the copy-I cycle E → Y → Z → E taking up one first substrate is no
    # potential; without its first transition the cycle is broken and one exists.
    change = zeros(Int, 1, 10)
    change[1, 4] = 1
    turning = EnzymeRates._Network(net.nforms, net.ncopies, net.present, net.from, net.to,
        net.edge, net.conformation, net.forward, net.fast, net.vanishing, change, net.takes)
    @test EnzymeRates._potential_on_sccs(net, trues(10), 1)
    @test !EnzymeRates._potential_on_sccs(turning, trues(10), 1)
    broken = BitVector([1, 1, 1, 0, 1, 1, 1, 1, 1, 1])
    @test EnzymeRates._potential_on_sccs(turning, broken, 1)
end
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T5_red.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_admission.jl")'`
Expected: Task 4's testsets pass; each new one errors, with ``UndefVarError: `_admission_face`
not defined in `EnzymeRates` `` (through `_testhelper_admission_faults`), ``UndefVarError:
`_PUMPED` not defined``, ``UndefVarError: `_BRIDGE` not defined`` or ``UndefVarError:
`_face_sets` not defined``.

- [ ] **Step 4: Add the face faults**

In the docstring of `_Fault` (Task 4), replace its last lines

```julia
(`_ARC_CYCLE`), or forms free E's component does not reach along the tagged bindings
(`_ISLAND`).
"""
```

with

```julia
(`_ARC_CYCLE`), or forms free E's component does not reach along the tagged bindings
(`_ISLAND`); on a face, in admission's face condition, a second closed class of the
ε-network (`_EPS_NONUNIQUE`), a stranded basin a step feeds that does not run back from the
law's forms (`_PUMPED`), or tagged steps that shorten the fed forms' cheapest spanning
in-tree where no potential absorbs their turnover (`_BRIDGE`).
"""
```

then add these three lines to the `@enum _Fault` block, after `    _ISLAND`:

```julia
    _EPS_NONUNIQUE
    _PUMPED
    _BRIDGE
```

and widen `_ADMISSION_FAULTS` and its docstring:

```julia
"""The faults admission reports (spec section 7); `_admission_message` writes their text."""
const _ADMISSION_FAULTS = (_NO_LIMIT, _ARC_CYCLE, _ISLAND, _EPS_NONUNIQUE, _PUMPED,
                           _BRIDGE)
```

- [ ] **Step 5: Write the face condition**

In `src/admission.jl`, replace Task 3's `_bottom_sccs` (its docstring and function) with
these two functions, which move its Tarjan pass unchanged into `_strong_components` for T0+
to share:

```julia
"""Each node's strongly connected component in the digraph with successor lists `next`,
searched from the nodes `nodes` marks, numbered in the order Tarjan's algorithm closes
them; 0 for a node no search reaches."""
function _strong_components(next::Vector{Vector{Int}}, nodes::BitVector)
    n = length(next)
    index, low, comp = zeros(Int, n), zeros(Int, n), zeros(Int, n)
    onstack, stack, work = falses(n), Int[], Tuple{Int, Int}[]
    counter = ncomp = 0
    for s in 1:n
        nodes[s] && index[s] == 0 || continue
        counter += 1
        index[s] = low[s] = counter
        push!(stack, s)
        onstack[s] = true
        push!(work, (s, 1))
        while !isempty(work)
            u, i = work[end]
            if i <= length(next[u])
                work[end] = (u, i + 1)
                v = next[u][i]
                if index[v] == 0
                    counter += 1
                    index[v] = low[v] = counter
                    push!(stack, v)
                    onstack[v] = true
                    push!(work, (v, 1))
                elseif onstack[v]
                    low[u] = min(low[u], index[v])
                end
            else
                pop!(work)
                isempty(work) || (low[work[end][1]] = min(low[work[end][1]], low[u]))
                if low[u] == index[u]
                    ncomp += 1
                    while true
                        v = pop!(stack)
                        onstack[v] = false
                        comp[v] = ncomp
                        v == u && break
                    end
                end
            end
        end
    end
    comp
end

"""The closed classes (bottom strongly connected components, `_strong_components`) of the
digraph with successor lists `next`, over the nodes `nodes` marks, each as its sorted
nodes; a node with no successor is a class of its own."""
function _bottom_sccs(next::Vector{Vector{Int}}, nodes::BitVector)
    comp = _strong_components(next, nodes)
    ncomp = maximum(comp; init = 0)
    leaves = falses(ncomp)
    for u in eachindex(next), v in next[u]
        nodes[u] && comp[u] != comp[v] && (leaves[comp[u]] = true)
    end
    classes = [Int[] for _ in 1:ncomp]
    for u in eachindex(next)
        nodes[u] && push!(classes[comp[u]], u)
    end
    [classes[k] for k in 1:ncomp if !leaves[k]]
end
```

Then append to `src/admission.jl`, after Task 4's `_admission_message`:

```julia
# ─── Admission: the face condition B v3 (research section 1.5) ─────────

"""The nodes of `net` that the nodes `start` reach along the transitions `keep` flags,
`start` included; against the transitions' direction when `backward`."""
function _reach(net::_Network, keep::BitVector, start::BitVector, backward::Bool)
    seen = copy(start)
    todo = findall(seen)
    while !isempty(todo)
        u = pop!(todo)
        for t in eachindex(net.from)
            keep[t] || continue
            a, b = backward ? (net.to[t], net.from[t]) : (net.from[t], net.to[t])
            a == u && !seen[b] && (seen[b] = true; push!(todo, b))
        end
    end
    seen
end

"""
    _face_sets(net, law_nodes, face) -> (; alive, fed, n0, stranded)

What admission's face condition reads on face `face` of the ε-network `net`, whose law
network holds the nodes `law_nodes`: the transitions alive there (taking up no zeroed
name), the nodes free E's active copy reaches along them, the alive transitions that do
not vanish, and the stranded basin: fed nodes outside the law network that reach it along
no chain of alive transitions that do not vanish.
"""
function _face_sets(net::_Network, law_nodes::BitVector, face::UInt64)
    alive = BitVector([(net.takes[t] & face) == 0 for t in eachindex(net.takes)])
    n0 = alive .& .!net.vanishing
    start = falses(length(net.present))
    start[_node(net, 1, 1)] = true
    fed = _reach(net, alive, start, false)
    drains = _reach(net, n0, law_nodes .& net.present, true)
    (; alive, fed, n0, stranded = fed .& .!law_nodes .& .!drains)
end

"""Copy `c`'s alive transitions with both ends in `nodes`: `van`, those that vanish with
their tags, and `held`, those T0+ reads (not vanishing, or fast)."""
function _copy_transitions(net::_Network, alive::BitVector, nodes::BitVector, c::Int)
    van, held = falses(length(net.from)), falses(length(net.from))
    for t in eachindex(net.from)
        alive[t] && net.conformation[t] == c && nodes[net.from[t]] && nodes[net.to[t]] ||
            continue
        van[t] = net.vanishing[t]
        held[t] = !net.vanishing[t] || net.fast[t]
    end
    (; van, held)
end

"""
The cost of a cheapest spanning out-arborescence of the nodes `1:n` rooted at `root` over
`arcs`, each `(from, to, cost)`, by Chu–Liu/Edmonds: take each node's cheapest entering
arc, contract the cycles they close and repeat; -1 when some node has no entering arc.
"""
function _min_arborescence(n::Int, arcs::Vector{NTuple{3, Int}}, root::Int)
    INF = typemax(Int) ÷ 4
    cost = 0
    while true
        inw, pre = fill(INF, n), zeros(Int, n)
        for (u, v, w) in arcs
            u != v && w < inw[v] && (inw[v] = w; pre[v] = u)
        end
        for v in 1:n
            v != root && inw[v] == INF && return -1
        end
        inw[root] = 0
        id, vis = zeros(Int, n), zeros(Int, n)
        cnt = 0
        for v in 1:n
            cost += inw[v]
            x = v
            while vis[x] != v && id[x] == 0 && x != root
                vis[x] = v
                x = pre[x]
            end
            if x != root && id[x] == 0
                cnt += 1
                y = pre[x]
                while y != x
                    id[y] = cnt
                    y = pre[y]
                end
                id[x] = cnt
            end
        end
        cnt == 0 && return cost
        for v in 1:n
            id[v] == 0 && (cnt += 1; id[v] = cnt)
        end
        next = NTuple{3, Int}[]
        for (u, v, w) in arcs
            id[u] != id[v] && push!(next, (id[u], id[v], w - inw[v]))
        end
        arcs, n, root = next, cnt, id[root]
    end
end

"""
    _fewest_slow(net, nodes, keep) -> Int

MINSLOW (research section 1.5): the fewest slow transitions in a spanning in-tree of the
nodes `nodes` over the transitions `keep` flags (cost 1 slow, 0 fast), any root; -1 when no
in-tree spans them. A virtual root joined to every node at a cost above any tree's turns
the in-tree into Chu–Liu/Edmonds's out-arborescence (`_min_arborescence`) on reversed
transitions; a cost of two such joins means no single root serves.
"""
function _fewest_slow(net::_Network, nodes::BitVector, keep::BitVector)
    at = zeros(Int, length(nodes))
    n = 0
    for u in eachindex(nodes)
        nodes[u] && (n += 1; at[u] = n)
    end
    arcs = NTuple{3, Int}[]
    for t in eachindex(net.from)
        a, b = at[net.from[t]], at[net.to[t]]
        keep[t] && a != 0 && b != 0 && push!(arcs, (b, a, net.fast[t] ? 0 : 1))
    end
    for u in 1:n
        push!(arcs, (n + 1, u, n + 1))
    end
    c = _min_arborescence(n + 1, arcs, n + 1)
    c < 0 || c >= 2 * (n + 1) ? -1 : c - (n + 1)
end

"""
    _potential_on_sccs(net, keep, x) -> Bool

T0+ (research section 1.5): whether the uptake of name `x` (a `change` row,
`_turnover_name`) is a potential on every strongly connected part
(`_strong_components`) of the transitions `keep` flags, so no cycle of them turns the
reaction over.
"""
function _potential_on_sccs(net::_Network, keep::BitVector, x::Int)
    n = length(net.present)
    next = [Int[] for _ in 1:n]
    for t in eachindex(net.from)
        keep[t] && push!(next[net.from[t]], net.to[t])
    end
    comp = _strong_components(next, net.present)
    adj = [Tuple{Int, Int}[] for _ in 1:n]
    for t in eachindex(net.from)
        u, v = net.from[t], net.to[t]
        keep[t] && comp[u] == comp[v] || continue
        push!(adj[u], (v, net.change[x, t]))
        push!(adj[v], (u, -net.change[x, t]))
    end
    h, seen = zeros(Int, n), falses(n)
    for s in 1:n
        seen[s] && continue
        seen[s] = true
        todo = [s]
        while !isempty(todo)
            u = pop!(todo)
            for (v, w) in adj[u]
                if !seen[v]
                    seen[v] = true
                    h[v] = h[u] + w
                    push!(todo, v)
                elseif h[v] != h[u] + w
                    return false
                end
            end
        end
    end
    true
end

"""Whether the transitions `vanishing` shorten the cheapest spanning in-tree of `fed` over
`base` (`_fewest_slow`): R1's failure for one conformation's vanishing transitions, joint
R1's for both (research section 1.5, X4)."""
function _bridges_jointly(net::_Network, fed::BitVector, base::BitVector,
                          vanishing::BitVector)
    without = _fewest_slow(net, fed, base)
    with = _fewest_slow(net, fed, base .| vanishing)
    with >= 0 && (without < 0 || with < without)
end

"""
The edges of the vanishing transitions on a cheapest spanning in-tree of `nodes` over
`base` and `van` (research section 1.5, `:bridge`), in transition order: an
inclusion-minimal set of the transitions `van` with which `base` still gives that tree's
cost (`_fewest_slow`), so every such tree uses each of them.
"""
function _bridging_edges(net::_Network, nodes::BitVector, base::BitVector, van::BitVector)
    best = _fewest_slow(net, nodes, base .| van)
    kept = copy(van)
    for t in eachindex(van)
        kept[t] || continue
        kept[t] = false
        _fewest_slow(net, nodes, base .| kept) == best || (kept[t] = true)
    end
    unique!([net.edge[t] for t in findall(kept)])
end

"""Of two or more closed classes, the first that no `fed` node lies in, else the second:
the class an `_EPS_NONUNIQUE` verdict names."""
function _unfed_class(closed::Vector{Vector{Int}}, fed::BitVector)
    for cls in closed
        any(view(fed, cls)) || return cls
    end
    closed[2]
end

"""How a message names face `face` of `m`, as the class-3 messages do: nothing for the
interior, else `at zero ` and its zeroed names (`_face_names`), then a comma."""
_face_phrase(m::Mechanism, face::UInt64) =
    face == 0 ? "" : "at zero " * join(String.(_face_names(m, face)), ", ", " and ") * ", "

"""
    _face_condition(net, law_nodes, face, x) -> _Verdict

Admission's face condition, B v3 (research section 1.5), on face `face` of the ε-network
`net` of an allosteric mechanism, whose law network holds the nodes `law_nodes`; `x` is
the `change` row of its first substrate (`_turnover_name`). In order:
- `_EPS_NONUNIQUE`: the alive transitions leave two closed classes, so at every ε the
  steady state depends on where the enzyme starts;
- `_PUMPED`: a stranded basin (`_face_sets`) is entered by, or holds, an alive transition
  whose reverse is dead, or is entered from outside it and the law network, so its weight
  depends on how fast its exits vanish; a basin entered only reversibly from the law
  network is sealed: its weight vanishes, and it leaves the fed nodes;
- `_BRIDGE`: in some conformation the vanishing transitions shorten the cheapest spanning
  in-tree of the fed nodes (R1 fails, `_bridges_jointly`) and the first substrate's uptake
  is no potential on the strongly connected parts of its alive transitions that are not
  vanishing or are fast (T0+ fails, `_potential_on_sccs`), led by the first of its
  vanishing steps on a cheapest tree (`_bridging_edges`); or both conformations pass R1 and
  their vanishing transitions together shorten it (joint R1, no conformation or step).
A face whose fed nodes no in-tree spans passes: class 3 decides it. The verdict leaves the
face index 0 (`_admission_face` sets it).
"""
function _face_condition(net::_Network, law_nodes::BitVector, face::UInt64, x::Int)
    s = _face_sets(net, law_nodes, face)
    closed = _closed_classes(net, s.alive)
    length(closed) > 1 && return _Verdict(_EPS_NONUNIQUE, _LIMIT_UNTESTED, 0,
                                          _copy_of(net, _unfed_class(closed, s.fed)[1]), 0)
    for t in eachindex(net.from)
        u, w = net.from[t], net.to[t]
        s.alive[t] && s.fed[u] && s.stranded[w] || continue
        s.alive[isodd(t) ? t + 1 : t - 1] && (s.stranded[u] || law_nodes[u]) && continue
        return _Verdict(_PUMPED, _LIMIT_UNTESTED, 0, net.conformation[t], net.edge[t])
    end
    nodes = s.fed .& .!s.stranded
    _fewest_slow(net, nodes, s.n0) < 0 && return _WELL_POSED
    passed = 0
    for c in 1:net.ncopies
        (; van, held) = _copy_transitions(net, s.alive, nodes, c)
        any(van) || continue
        if !_bridges_jointly(net, nodes, s.n0, van)
            passed += 1
        elseif !_potential_on_sccs(net, held, x)
            return _Verdict(_BRIDGE, _LIMIT_UNTESTED, 0, c,
                            first(_bridging_edges(net, nodes, s.n0, van)))
        end
    end
    passed == 2 && _bridges_jointly(net, nodes, s.n0, s.alive .& net.vanishing) &&
        return _Verdict(_BRIDGE, _LIMIT_UNTESTED, 0, 0, 0)
    _WELL_POSED
end

"""
    _admission_face(am, eps, law, faces, i) -> _Verdict

Admission's face condition (`_face_condition`) on face `faces[i]` of `am`, from its
ε-network `eps` and law network `law`; a rejection carries the face index `i`.
"""
function _admission_face(am::AllostericMechanism, eps::_Network, law::_Network,
                         faces::Vector{UInt64}, i::Int)
    v = _face_condition(eps, law.present, faces[i], _turnover_name(am.catalytic))
    v.fault === _NO_FAULT ? v : _Verdict(v.fault, v.limit, i, v.conformation, v.edge)
end

"""
    _admission(am) -> _Verdict

Admission alone (spec section 7): the interior checks (`_admission_interior`), then the
face condition on every face (`_admission_face`), class 3 not run, as the tagging and
change-tag moves filter their children (research section 1.9, ok(t)).
"""
function _admission(am::AllostericMechanism)
    v = _admission_interior(am)
    v.fault === _NO_FAULT || return v
    eps, law = _eps_network(am), _law_network(am)
    faces = _faces(am.catalytic)
    for i in eachindex(faces)
        v = _admission_face(am, eps, law, faces, i)
        v.fault === _NO_FAULT || return v
    end
    _WELL_POSED
end
```

- [ ] **Step 6: Write the face messages**

In `_admission_message` (Task 4), insert these branches before its final `else`:

```julia
    elseif v.fault in (_EPS_NONUNIQUE, _PUMPED, _BRIDGE)
        eps, law = _eps_network(am), _law_network(am)
        faces = _faces(m)
        face = faces[v.face]
        phrase = _face_phrase(m, face)
        s = _face_sets(eps, law.present, face)
        if v.fault === _EPS_NONUNIQUE
            unfed = _unfed_class(_closed_classes(eps, s.alive), s.fed)
            forms = [_form_text(m, _form_of(eps, u)) *
                     (_copy_of(eps, u) == 1 ? " (active)" : " (inactive)") for u in unfed]
            return "AllostericMechanism: $(phrase)the forms $(join(forms, ", ")) " *
                   "hold the enzyme they receive and are reached only through tagged " *
                   "steps, so the steady state depends on where the enzyme starts."
        elseif v.fault === _PUMPED
            basin = [u for u in findall(s.stranded) if _copy_of(eps, u) == c]
            exits = unique!([eps.edge[t] for t in eachindex(eps.from)
                             if s.alive[t] && eps.vanishing[t] && eps.from[t] in basin &&
                                !s.stranded[eps.to[t]]])
            out = isempty(exits) ?
                "which no step running there leaves, so the enzyme collects in them" :
                "which leave only through the vanishing steps $(lines(exits)) while " *
                "nothing runs back the way it came, so their weight depends on how fast " *
                "those steps vanish"
            return "AllostericMechanism: $(phrase)$(line(k)) feeds the forms " *
                   "$(join([_form_text(m, _form_of(eps, u)) for u in basin], ", ")) of " *
                   "the $conf conformation, $out."
        elseif c == 0
            return "AllostericMechanism: $(phrase)the tagged steps of both " *
                   "conformations together bridge forms that each conformation's other " *
                   "steps keep apart, so the limit depends on how fast those steps vanish."
        end
        nodes = s.fed .& .!s.stranded
        (; van) = _copy_transitions(eps, s.alive, nodes, c)
        ks = _bridging_edges(eps, nodes, s.n0, van)
        bridge = length(ks) == 1 ? "step $(line(k)) bridges" : "steps $(lines(ks)) bridge"
        msg = "AllostericMechanism: $(phrase)the $tag $bridge forms of the $conf " *
              "conformation that its other steps keep apart, so the limit depends on how " *
              "fast the tagged constants vanish."
        # Offer the steady-state fix only where it passes this face: R2funnel keeps its
        # bridge with its tagged S binding steady-state (R2funnel_SS).
        re = _flags(m)
        re[_group_members(m)[m.group[k]]] .= false
        fixed = _admission_face(am, _eps_network(am, re), _law_network(am, re), faces,
                                v.face)
        return fixed.fault === _NO_FAULT ?
               msg * " Make the group of $(line(k)) steady-state." : msg
```

- [ ] **Step 7: Wire the face condition into the guard**

In `src/admission.jl`, replace the docstring of `_guard` (Task 4's) and its allosteric method
with the following; the plain method `_guard(m::Mechanism)` between them stays as it is. The
allosteric method now runs the faces itself, since `_faces_verdict` also serves
`_faces_pass`, which the flip move reads for class 3 alone:

```julia
"""
    _guard(m) -> _Verdict

The well-posedness guard (spec section 8) on `m`, plain or allosteric: class 1 (a cycle
turns the reaction over, `_has_turnover_cycle`); class 2 (no cycle of fast transitions
does, `_re_turnover_transitions`, on the ε-network, where an excluded rapid-equilibrium
binding stays fast while its affinity vanishes); for an allosteric mechanism, admission's
interior (`_admission_interior`); then on every face class 3 (`_check_face`, on the law
network) and, for an allosteric mechanism, admission's face condition (`_admission_face`,
on the ε-network), the first fault winning. `_WELL_POSED` when every check passes. It
reads the reaction graph through masks and builds nothing.
"""
```

```julia
function _guard(am::AllostericMechanism)
    m = am.catalytic
    _has_turnover_cycle(m, trues(length(m.edges))) ||
        return _Verdict(_NO_CYCLE, _LIMIT_UNTESTED, 0, 0, 0)
    eps = _eps_network(am)
    v = _re_turnover_verdict(eps)
    v === _WELL_POSED || return v
    v = _admission_interior(am)
    v === _WELL_POSED || return v
    law, faces = _law_network(am), _faces(m)
    for i in eachindex(faces)
        v = _check_face(am, law, faces, i)
        v.fault === _NO_FAULT || return v
        v = _admission_face(am, eps, law, faces, i)
        v.fault === _NO_FAULT || return v
    end
    _WELL_POSED
end
```

- [ ] **Step 8: Run the tests to verify they pass**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T5_green.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_admission.jl")'`
Expected: every testset passes, pristine. Compare the file's Test Summary time with
`P4_T4_green.log`: at most 3 s more without "conservative rejections the oracle admits stay
rejected", at most 8 s more with it (its time shows in the summary's nested line).

- [ ] **Step 9: Report the face condition's cost**

Create `$S/phase45/probes/P4_T5_cost.jl`:

```julia
# ABOUTME: Times admission, interior and every face, per mechanism over the uni-uni grid,
# ABOUTME: for the cost report of Task 5; run from the guard worktree under TestEnv.
using Test, EnzymeRates, LinearAlgebra, Random
include(joinpath(pwd(), "test", "mechanism_definitions_for_test_enzyme_derivation.jl"))
grid = _testhelper_uniuni_grid()
foreach(EnzymeRates._admission, grid)
t = minimum(@elapsed(foreach(EnzymeRates._admission, grid)) for _ in 1:20)
println("admission, interior and every face: ", round(1e6 * t / length(grid); digits = 1),
        " µs per mechanism over ", length(grid))
```

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T5_cost.log 900 --project=. -e 'using TestEnv; TestEnv.activate(); include("'$S'/phase45/probes/P4_T5_cost.jl")'`
Expected: one line `admission, interior and every face: <t> µs per mechanism over 108`, below the
prototype's 149 µs; the commit body reports it.

- [ ] **Step 10: Record the conservative rejections (P4-9, P4-9a)**

In `ROADMAP.md`, append this entry as the last bullet of `## Ideas`, before `## Decided against`
(one blank line before the heading, as the other entries have):

```markdown
- **Admit the equilibrium bridges and the vanishing-fed basins that admission rejects.**
  Admission's face condition (`_face_condition`, `src/admission.jl`) is sound but not
  exact: it rejects two rare classes whose limit exists. In an equilibrium bridge (R-B1) a
  tagged binding joins forms its conformation's other steps keep apart, yet those forms
  already sit at equilibrium, so the limit is deletion's (`bridge_deadend` in
  `test/test_admission.jl`). In a stranded basin fed only through vanishing sources (R-B2)
  the basin's weight vanishes along every path, but no structural test proves it (`BS1`).
  The phase-4 research found R-B1 on 332 faces and R-B2 on 21, all in fuzzed or hand-built
  mechanisms, none among the census mechanisms or their `:OnlyI` variants. No exact
  narrower structural rule is known: the one candidate tried, T0++, admits faces the
  ε-path oracle rejects. The exact test compares Freidlin–Wentzell exponents over the
  admissible directions: the basin vanishes when, for every rise of the tagged constants,
  every per-member rate and every order of the two limits, its cheapest in-tree costs more
  than the law network's; each comparison is one linear-programming feasibility test, not
  yet prototyped.
```

In `docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md`, P4-9a amends every
passage that asks for an exact face condition. In section 7 "Admission", item 2
("Uniqueness"), replace the lines

```markdown
   rejects cases such as R2 below. That face condition is not designed yet (the slice
   admits R2); phase 4 designs it and ships it only when it matches the RE-speed and ε-path
   oracles.
```

with

```markdown
   rejects cases such as R2 below. Phase 4's face condition (B v3) rejects R2 and is sound
   against the ε-path oracle; conservative on R-B1 and R-B2, listed: it also rejects
   equilibrium bridges (R-B1) and stranded basins fed only through vanishing steps (R-B2)
   whose limit exists, and ROADMAP lists both.
```

In section 8 "Phase 4's first task", replace the line

```markdown
Make the face check exact against the RE-speed oracle. The probes found false alarms (6 on
```

with

```markdown
Make the face check exact against the RE-speed oracle. Admission's face condition is
sound; conservative on R-B1 and R-B2, listed: it rejects equilibrium bridges (R-B1) and
stranded basins fed only through vanishing steps (R-B2) that the ε-path oracle admits, and
ROADMAP lists both. The probes found false alarms (6 on
```

in section 15's table, in the phase-4 row, replace the last cell

```markdown
| face check and admission's face condition exact first |
```

with

```markdown
| face check exact first; admission's face condition sound, conservative on R-B1 and R-B2, listed |
```

and in section 17, replace the row

```markdown
| The face check or admission's face condition rejects valid mechanisms or misses invalid ones | Phase 4 starts by making both exact against the RE-speed and ε-path oracles |
```

with

```markdown
| The face check or admission's face condition rejects valid mechanisms or misses invalid ones | Phase 4 makes the face check exact against the RE-speed oracle and admission's face condition sound against the ε-path oracle; conservative on R-B1 and R-B2, listed in ROADMAP |
```

- [ ] **Step 11: Run the full suite**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T5_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`
Expected: `exit 0`, `Testing EnzymeRates tests passed`, pristine, Aqua and JET green; report
`maxrss_kb`.

- [ ] **Step 12: Commit**

Write the measured figures in place of `<t>` (Step 9) and `<kb>` (Step 11).

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
git add src/admission.jl test/test_admission.jl \
    test/mechanism_definitions_for_test_enzyme_derivation.jl ROADMAP.md \
    docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md
git commit -m "$(cat <<'EOF'
Check admission's face condition on every face

After class 3 accepts a face, B v3 rejects an epsilon-network with two closed
classes, a stranded basin that is pumped rather than sealed, and a tagged step
that bridges forms without a potential to absorb its turnover; joint R1 is
tested on a network in the guard's own type, since a mechanism that needs it
holds a catalytic :OnlyI group, which the constructor still refuses.
Every :OnlyA row agrees with its oracle; the R-B1 and R-B2 rejections are
conservative and listed in ROADMAP and the spec. The bridge message names the
tagged steps every cheapest tree uses and offers the steady-state fix only
where it passes the face; the pumped message names the basin's vanishing
exits. T0+ shares class 3's Tarjan pass, moved into _strong_components.
Cost: <t> us per mechanism. Full-suite maxrss: <kb> KB.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
)"
```

---

### Task 6: Class 5, the support test on the derived law, not wired

On every face, some monomial of the derived denominator in lowest terms must hold none of the
face's names (RES §1.6; spec §8 item 5). Task 8 calls it inside `_num_den_exprs`.

**Files:**
- Modify: `src/admission.jl`
- Test: `test/test_admission.jl`

**Interfaces:**
- Consumes: Task 2's `_faces`, `_face_names`, `_Verdict`, `_guard_message(mech::_AnyMech, v)`;
  Task 4's `_Fault` docstring; Task 5's `_face_phrase`; phase 3's
  `_raw_symbolic_rate_polys(m::Mechanism)` (returns `(num, den, d_free)` in lowest terms),
  `_mwc_state_polys(am, mets::Set{Symbol}, kept)` (returns `(num_A, den_A, num_I, den_I,
  d_A, d_I)`), `_kept_members`, `metabolites`, `catalytic`, `POLY`, `poly_one`, `_mono`;
  the tests read `_testhelper_lift` and `MECHANISM_TEST_SPECS`.
- Produces:
  ```julia
  # _Fault gains, after _BRIDGE: _DENOMINATOR_VANISHES
  _law_branches(m::Mechanism)::Vector{Tuple{POLY, POLY}}              # [(den, poly_one())]
  _law_branches(am::AllostericMechanism)::Vector{Tuple{POLY, POLY}}   # [(den_A, d_I), (den_I, d_A)]
  _keeps_monomial(p::POLY, z::Vector{Symbol})::Bool
  _vanishing_face(branches::Vector{Tuple{POLY, POLY}}, faces::Vector{Vector{Symbol}})::Int
  _support_verdict(m::Mechanism, branches::Vector{Tuple{POLY, POLY}})::_Verdict
  _support_verdict(am::AllostericMechanism, branches::Vector{Tuple{POLY, POLY}})::_Verdict
  _support_message(who::String, m::Mechanism, v::_Verdict)::String  # m the catalytic graph
  ```
  `_guard_message` hands `_DENOMINATOR_VANISHES` to `_support_message`.

**Gate:** `ss_dead_end_at_face`'s law fails at products = 0 (face 2); every spec of
`MECHANISM_TEST_SPECS` but the two trap specs keeps a monomial on every face; cost per law
reported (research 1.4-2.3 µs); test_admission.jl grows by at most 1 s.

- [ ] **Step 1: Write the failing tests**

Append to `test/test_admission.jl`:

```julia
@testset "class 5 finds the face where a denominator loses every monomial" begin
    mono(pairs...) = EnzymeRates._mono(pairs...)
    unit = EnzymeRates.poly_one()
    faces = [Symbol[], [:P, :Q], [:S]]
    # K·S + P·Q: zeroing P and Q leaves K·S, zeroing S leaves P·Q.
    den = EnzymeRates.POLY(mono(:K => 1, :S => 1) => 1, mono(:P => 1, :Q => 1) => 1)
    @test EnzymeRates._vanishing_face([(den, unit)], faces) == 0
    # P·S + Q·S: zeroing P and Q leaves no term.
    den = EnzymeRates.POLY(mono(:P => 1, :S => 1) => 1, mono(:Q => 1, :S => 1) => 1)
    @test EnzymeRates._vanishing_face([(den, unit)], faces) == 2
    # Two conformations, each denominator times the other's free-enzyme weight: S·1 and
    # (1 + S)·P. Zeroing S keeps 1·P, zeroing P keeps S·1, zeroing both keeps nothing.
    den_A = EnzymeRates.POLY(mono(:S => 1) => 1)
    den_I = EnzymeRates.POLY(mono() => 1, mono(:S => 1) => 1)
    d_A = EnzymeRates.POLY(mono(:P => 1) => 1)
    branches = [(den_A, unit), (den_I, d_A)]
    @test EnzymeRates._vanishing_face(branches, [Symbol[], [:S], [:P]]) == 0
    @test EnzymeRates._vanishing_face(branches, [Symbol[], [:S, :P]]) == 2
end

@testset "the law of a steady-state dead end vanishes at zero products" begin
    # E(Q) + P ⇌ E(P, Q) puts every form downstream of the chemistry behind a product: each
    # term of the denominator holds P or Q, so the law is 0/0 with both zero (face 2).
    m = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q <--> E(Q)
            E(P) + Q <--> E(P, Q)
        end
    end
    v = EnzymeRates._support_verdict(m, EnzymeRates._law_branches(m))
    @test v == EnzymeRates._Verdict(EnzymeRates._DENOMINATOR_VANISHES,
                                    EnzymeRates._LIMIT_UNTESTED, 2, 0, 0)
    @test EnzymeRates._guard_message(m, v) ==
          "Mechanism: at zero P and Q, every term of the rate law's denominator " *
          "vanishes, so the rate is 0/0 there."
    am = _testhelper_tag_sets["AD3a"]
    @test length(EnzymeRates._law_branches(am)) == 2
    w = EnzymeRates._Verdict(EnzymeRates._DENOMINATOR_VANISHES,
                             EnzymeRates._LIMIT_UNTESTED, 2, 0, 0)
    @test EnzymeRates._guard_message(am, w) ==
          "AllostericMechanism: at zero P, every term of the rate law's denominator " *
          "vanishes, so the rate is 0/0 there."
end

@testset "every spec's denominator keeps a monomial on every face" begin
    traps = ("Numerator: redundant SS-bind", "LDH i-state NonequalAI 6-group")
    @testset "$(spec.name)" for spec in filter(s -> s.name ∉ traps, MECHANISM_TEST_SPECS)
        x = _testhelper_lift(spec.mechanism)
        @test EnzymeRates._support_verdict(x, EnzymeRates._law_branches(x)) ==
              EnzymeRates._WELL_POSED
    end
end
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T6_red.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_admission.jl")'`
Expected: Tasks 4-5's testsets pass; the three new testsets error with ``UndefVarError:
`_vanishing_face` not defined in `EnzymeRates` `` and ``UndefVarError: `_law_branches` not
defined``.

- [ ] **Step 3: Write class 5**

In the docstring of `_Fault` (Task 5), replace its last lines

```julia
in-tree where no potential absorbs their turnover (`_BRIDGE`).
"""
```

with

```julia
in-tree where no potential absorbs their turnover (`_BRIDGE`); on a face of the derived
law (class 5), a denominator that loses every term (`_DENOMINATOR_VANISHES`).
"""
```

and add the line `    _DENOMINATOR_VANISHES` to the `@enum _Fault` block, after `    _BRIDGE`. In
`src/admission.jl`, insert before the section line `# ─── Guard messages` (Task 2), after
Task 5's `_admission`:

```julia
# ─── Class 5: the support test on the derived law (spec section 8) ─────

"""
    _law_branches(m) -> Vector{Tuple{POLY, POLY}}

The derived denominator of `m`'s rate law as branches of positive terms, so that it keeps
a monomial free of a face's names exactly when some branch's two factors both do: a plain
mechanism's denominator in lowest terms (`_raw_symbolic_rate_polys`) times 1; an allosteric
one's two conformations, each denominator times the other conformation's free-enzyme
weight (`_mwc_state_polys`). Powers of a branch and the regulatory sites' factors, whose 1
never vanishes, keep the same faces.
"""
_law_branches(m::Mechanism) = [(_raw_symbolic_rate_polys(m)[2], poly_one())]
function _law_branches(am::AllostericMechanism)
    _, den_A, _, den_I, d_A, d_I =
        _mwc_state_polys(am, Set{Symbol}(metabolites(am.catalytic)), _kept_members(am))
    [(den_A, d_I), (den_I, d_A)]
end

"""Whether some monomial of `p` holds none of the names `z` to a positive power."""
function _keeps_monomial(p::POLY, z::Vector{Symbol})
    for mono in keys(p)
        any(((s, e),) -> e > 0 && s in z, mono) || return true
    end
    false
end

"""The first face, each given by its zeroed names, on which no branch of `branches`
(`_law_branches`) keeps a monomial in both factors; 0 when every face keeps one."""
function _vanishing_face(branches::Vector{Tuple{POLY, POLY}},
                         faces::Vector{Vector{Symbol}})
    for (i, z) in enumerate(faces)
        any(b -> _keeps_monomial(b[1], z) && _keeps_monomial(b[2], z), branches) ||
            return i
    end
    0
end

"""
    _support_verdict(m, branches) -> _Verdict

Class 5 (spec section 8): `_DENOMINATOR_VANISHES` on the first face of `m` (`_faces`)
where the denominator `branches` (`_law_branches`) loses every monomial, so the law is 0/0
there; it catches the steady-state dead ends class 3's structure misses.
"""
function _support_verdict(m::Mechanism, branches::Vector{Tuple{POLY, POLY}})
    i = _vanishing_face(branches, [_face_names(m, f) for f in _faces(m)])
    i == 0 ? _WELL_POSED : _Verdict(_DENOMINATOR_VANISHES, _LIMIT_UNTESTED, i, 0, 0)
end
_support_verdict(am::AllostericMechanism, branches::Vector{Tuple{POLY, POLY}}) =
    _support_verdict(am.catalytic, branches)

"""The message of class 5's verdict `v` on a mechanism that `who` names and whose catalytic
graph is `m`: the face where every term of the law's denominator vanishes."""
_support_message(who::String, m::Mechanism, v::_Verdict) =
    "$who: $(_face_phrase(m, _faces(m)[v.face]))every term of the rate law's " *
    "denominator vanishes, so the rate is 0/0 there."
```

In `_guard_message(mech::_AnyMech, v::_Verdict)`, directly after the two lines Task 4 inserted
(`mech isa AllostericMechanism && v.fault in _ADMISSION_FAULTS && ...`), insert:

```julia
    v.fault === _DENOMINATOR_VANISHES && return _support_message(who, catalytic(mech), v)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T6_green.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_admission.jl")'`
Expected: every testset passes, pristine; 40 spec testsets in "every spec's denominator keeps a
monomial on every face"; at most 1 s more than `P4_T5_green.log`.

- [ ] **Step 5: Report class 5's cost**

Create `$S/phase45/probes/P4_T6_cost.jl`:

```julia
# ABOUTME: Times class 5's support test per law over the specs, against deriving their
# ABOUTME: denominators, for the cost report of Task 6; run under TestEnv from the worktree.
using Test, EnzymeRates, LinearAlgebra, Random
include(joinpath(pwd(), "test", "mechanism_definitions_for_test_enzyme_derivation.jl"))
xs = [_testhelper_lift(s.mechanism) for s in MECHANISM_TEST_SPECS]
laws = [(x, EnzymeRates._law_branches(x)) for x in xs]
foreach(((x, b),) -> EnzymeRates._support_verdict(x, b), laws)
t = minimum(@elapsed(foreach(((x, b),) -> EnzymeRates._support_verdict(x, b), laws))
            for _ in 1:50)
d = minimum(@elapsed(foreach(EnzymeRates._law_branches, xs)) for _ in 1:3)
println("class 5: ", round(1e6 * t / length(laws); digits = 2), " µs per law over ",
        length(laws), "; ", round(100 * t / d; digits = 2), "% of deriving the branches")
```

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T6_cost.log 900 --project=. -e 'using TestEnv; TestEnv.activate(); include("'$S'/phase45/probes/P4_T6_cost.jl")'`
Expected: one line `class 5: <t> µs per law over 42; <p>% of deriving the branches`, `<t>` near
research's 1.4-2.3 µs and `<p>` under 1%.

- [ ] **Step 6: Run the full suite**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T6_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`
Expected: `exit 0`, tests passed, pristine, Aqua and JET green; report `maxrss_kb`.

- [ ] **Step 7: Commit**

Write the measured figures in place of `<t>`, `<p>` (Step 5) and `<kb>` (Step 6).

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
git add src/admission.jl test/test_admission.jl
git commit -m "$(cat <<'EOF'
Find the face where a derived denominator vanishes

Class 5 reads the derived law: on each face some monomial of the denominator,
or for an allosteric law of one conformation's denominator and the other's
free-enzyme weight, must hold none of the face's names. The steady-state dead
end vanishes at zero products; every other spec keeps a monomial on every face.
Cost: <t> us per law, <p>% of the derivation. Full-suite maxrss: <kb> KB.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
)"
```

---

### Task 7: Trap specs become rejection tests

The wired macros (Task 8) would throw on the two trap specs when the fixtures file loads, so
each moves to an unguarded rejection test and its well-posed variant takes its place in
`MECHANISM_TEST_SPECS` (P4-27; spec §8 "Phase 4's first task").

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (specs "Numerator:
  redundant SS-bind" and "LDH i-state NonequalAI 6-group" become their variants)
- Modify: `test/reference/allosteric_golden_reference.txt` (the LDH variant's two lines)
- Modify: `test/fixtures/canonical_text_golden.txt` (both specs' entries)
- Test: `test/test_admission.jl` (two rejection testsets; Task 6's spec loop loses its
  exclusion)

**Interfaces:**
- Consumes: Task 2's `_guard`, `_guard_message`; phase 3's `_step_text`; Task 3's
  `_testhelper_face_faults`, `_TRAP`, `_SPLIT`, `_LIMIT_SPEED`, the trap message naming the
  feeding step; Task 5's
  `_testhelper_admission_faults`, `_testhelper_admission_agreement`; Task 6's `_law_branches`,
  `_support_verdict`; Task 1's `_testhelper_respeed_verdicts(m; rng, law)` (a compiled
  mechanism for `law = :package`), `@_testhelper_unguarded`; phase 3's
  `_testhelper_check_law_against_oracle(m)`, `_testhelper_edge`, `_testhelper_lift`.
- Produces: specs "Numerator: redundant SS-bind, steady-state feeder" (`E(S1) + S2 <--> E(S1,
  S2)`; `expected_n_independent_params=7`) and "LDH i-state NonequalAI 6-group, steady-state
  NAD group" (its `:EqualAI` NAD group steady-state, 6 groups; `expected_n_independent_params=10`),
  in the trap specs' places.

**Gate:** values before pins: each variant's law equals the exact oracle to 1e-20, every face
passes the guard and agrees with both oracles; spec 14 rejected at substrates = 0 (`_TRAP`
naming `E(S1) + S2 ⇌ E(S1, S2)`) and S2 = 0 (`_SPLIT`), LDH at products = 0 (`_TRAP` naming
`E(Lactate) + NAD ⇌ E(Lactate, NAD)`) and NAD = 0 (`_SPLIT`), the RE-speed oracle `:pattern` on
each; the golden reference changes only in the LDH variant's two lines; the canonical golden
changes only in the two specs' entries; `rate_equation` stays allocation-free and under 120 ns
for both variants (measured on phase 2: 9.1 ns and 65.4 ns) and the string pins hold;
`length(MECHANISM_TEST_SPECS) == 42`.

- [ ] **Step 1: Write the failing tests**

Append to `test/test_admission.jl`:

```julia
@testset "the redundant steady-state binding spec traps substrate-free enzyme" begin
    rng = Random.MersenneTwister(73)
    NO, TR, SP = EnzymeRates._NO_FAULT, EnzymeRates._TRAP, EnzymeRates._SPLIT
    # With S1 and S2 zero, E(S1, S2) feeds E(S1) by the fast release of S2, and E(S1)
    # leaves only by its slow release of S1: a trap, whose limit depends on the RE speeds.
    # With S2 zero, the lowest forms are joined only through E(S1, S2): a split.
    m = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S1, S2
        products: P
        steps: begin
            E + S1 <--> E(S1)
            E + S2 ⇌ E(S2)
            E(S1) + S2 ⇌ E(S1, S2)
            E(S2) + S1 ⇌ E(S1, S2)
            E(S1, S2) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    @test _testhelper_face_faults(m) == [NO, NO, TR, NO, SP]
    @test _testhelper_respeed_verdicts(m; rng, law = :limit) ==
          [:agree, :agree, :pattern, :agree, :pattern]
    feed = _testhelper_edge(m, "ES1 + S2 ⇌ ES1S2")
    v = EnzymeRates._guard(m)
    @test (v.fault, v.limit, v.face, v.edge) ===
          (TR, EnzymeRates._LIMIT_SPEED, 3, feed)
    @test occursin(EnzymeRates._step_text(m, feed), EnzymeRates._guard_message(m, v))
    # The variant makes the feeding step steady-state: every face passes, and its law is
    # the exact steady state, before any pin.
    spec = only(s for s in MECHANISM_TEST_SPECS
                if s.name == "Numerator: redundant SS-bind, steady-state feeder")
    x = _testhelper_lift(spec.mechanism)
    _testhelper_check_law_against_oracle(spec.mechanism)
    @test _testhelper_face_faults(x) == fill(NO, 5)
    @test _testhelper_respeed_verdicts(spec.mechanism; rng, law = :package) ==
          fill(:agree, 5)
    @test EnzymeRates._guard(x) == EnzymeRates._WELL_POSED
    @test EnzymeRates._support_verdict(x, EnzymeRates._law_branches(x)) ==
          EnzymeRates._WELL_POSED
    @test length(EnzymeRates.fitted_params(spec.mechanism)) == 7
end

@testset "the LDH inactive-state spec traps lactate-bound enzyme" begin
    rng = Random.MersenneTwister(79)
    NO, TR, SP = EnzymeRates._NO_FAULT, EnzymeRates._TRAP, EnzymeRates._SPLIT
    # With Lactate and NAD zero, E(Lactate, NAD) feeds E(Lactate) by the fast release of
    # NAD, and E(Lactate) leaves only slowly: a trap. With NAD zero, a split.
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: NADH, Pyruvate
        products: Lactate, NAD
        catalytic_multiplicity: 4
        steps: begin
            (E + Lactate <--> E(Lactate),
             E(NADH) + Lactate <--> E(Lactate, NADH)) :: NonequalAI
            (E + NAD ⇌ E(NAD),
             E(Lactate) + NAD ⇌ E(Lactate, NAD),
             E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate)) :: EqualAI
            (E + NADH <--> E(NADH),
             E(Lactate) + NADH <--> E(Lactate, NADH),
             E(Pyruvate) + NADH <--> E(NADH, Pyruvate)) :: OnlyA
            (E + Pyruvate ⇌ E(Pyruvate),
             E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate),
             E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate)) :: EqualAI
            E(NAD) + Lactate ⇌ E(Lactate, NAD) :: EqualAI
            E(NADH, Pyruvate) <--> E(Lactate, NAD) :: OnlyA
        end
    end
    @test _testhelper_admission_faults(am) == [NO, TR, NO, NO, NO, NO, SP]
    @test _testhelper_respeed_verdicts(am; rng, law = :limit) ==
          [:agree, :pattern, :agree, :agree, :agree, :agree, :pattern]
    feed = _testhelper_edge(am, "ELactate + NAD ⇌ ELactateNAD")
    v = EnzymeRates._guard(am)
    @test (v.fault, v.limit, v.face, v.edge) ===
          (TR, EnzymeRates._LIMIT_SPEED, 2, feed)
    @test occursin(EnzymeRates._step_text(am.catalytic, feed),
                   EnzymeRates._guard_message(am, v))
    # The variant makes the :EqualAI NAD group steady-state, keeping its six groups: every
    # face passes the guard and agrees with both oracles, and its law is exact.
    spec = only(s for s in MECHANISM_TEST_SPECS
                if s.name == "LDH i-state NonequalAI 6-group, steady-state NAD group")
    x = _testhelper_lift(spec.mechanism)
    _testhelper_check_law_against_oracle(spec.mechanism)
    @test _testhelper_admission_faults(x) == fill(NO, 7)
    @test all(_testhelper_admission_agreement(x; rng))
    @test _testhelper_respeed_verdicts(spec.mechanism; rng, law = :package) ==
          fill(:agree, 7)
    @test EnzymeRates._guard(x) == EnzymeRates._WELL_POSED
    @test EnzymeRates._support_verdict(x, EnzymeRates._law_branches(x)) ==
          EnzymeRates._WELL_POSED
    @test length(EnzymeRates.fitted_params(spec.mechanism)) == 10
end
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T7_red.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_admission.jl")'`
Expected: the trap fixtures' assertions pass; each new testset then errors at
`only(s for s in MECHANISM_TEST_SPECS ...)` with `ArgumentError: Collection is empty, must
contain exactly 1 element`; every other testset passes.

- [ ] **Step 3: Replace the trap specs by their variants**

In `test/mechanism_definitions_for_test_enzyme_derivation.jl`, in spec 11d, replace

```julia
# 11d. Single-RE-segment with a redundant SS binding branch (one branch
#      RE, one SS). ODE is ground truth.
```

with

```julia
# 11d. Two binding orders into one complex: S1 then S2 at steady state, S2 then S1 at
#      rapid equilibrium. ODE is ground truth.
```

in its steps replace `            E(S1) + S2 ⇌ E(S1, S2)` with `            E(S1) + S2 <--> E(S1, S2)`,
and replace

```julia
        name="Numerator: redundant SS-bind",
```

with

```julia
        name="Numerator: redundant SS-bind, steady-state feeder",
```

and, in the same spec, `expected_n_independent_params=6,` with
`expected_n_independent_params=7,` (`expected_n_haldane_constraints=1` and
`expected_n_wegscheider_constraints=1` stay; measured on phase 2, `p2_counts2.log`).

In the spec "LDH i-state NonequalAI 6-group", replace the NAD group

```julia
            (E + NAD ⇌ E(NAD),
             E(Lactate) + NAD ⇌ E(Lactate, NAD),
             E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate)) :: EqualAI
```

(the group of that spec just above its `(E + NADH <--> E(NADH),` group) with

```julia
            (E + NAD <--> E(NAD),
             E(Lactate) + NAD <--> E(Lactate, NAD),
             E(Pyruvate) + NAD <--> E(NAD, Pyruvate)) :: EqualAI
```

and replace

```julia
        name="LDH i-state NonequalAI 6-group",
        mechanism=m,
        raw=raw,
        metabolite_names=[:NADH, :Pyruvate, :Lactate, :NAD],
        expected_n_states=9, expected_n_steps=13,
        expected_n_haldane_constraints=1,
        expected_n_wegscheider_constraints=2, expected_n_independent_params=9,
```

with

```julia
        name="LDH i-state NonequalAI 6-group, steady-state NAD group",
        mechanism=m,
        raw=raw,
        metabolite_names=[:NADH, :Pyruvate, :Lactate, :NAD],
        expected_n_states=9, expected_n_steps=13,
        expected_n_haldane_constraints=1,
        expected_n_wegscheider_constraints=2, expected_n_independent_params=10,
```

With both variants in place, in Task 6's testset "every spec's denominator keeps a monomial
on every face" (`test/test_admission.jl`) replace

```julia
    traps = ("Numerator: redundant SS-bind", "LDH i-state NonequalAI 6-group")
    @testset "$(spec.name)" for spec in filter(s -> s.name ∉ traps, MECHANISM_TEST_SPECS)
```

with

```julia
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
```

so that it checks all 42 specs.

- [ ] **Step 4: Regenerate the LDH variant's golden lines only**

Create `$S/phase45/probes/P4_T7_golden.jl`:

```julia
# ABOUTME: Rewrites the two lines of test/reference/allosteric_golden_reference.txt that pin
# ABOUTME: the LDH trap spec with its variant's name and rate-equation string (Task 7).
using Test, EnzymeRates, LinearAlgebra, Random
include(joinpath(pwd(), "test", "mechanism_definitions_for_test_enzyme_derivation.jl"))
path = joinpath(pwd(), "test", "reference", "allosteric_golden_reference.txt")
lines = readlines(path)
i = only(findall(==("### LDH i-state NonequalAI 6-group"), lines))
spec = only(s for s in MECHANISM_TEST_SPECS
            if s.name == "LDH i-state NonequalAI 6-group, steady-state NAD group")
lines[i] = "### " * spec.name
lines[i + 1] = "RATE_EQUATION_STRING " *
               replace(EnzymeRates.rate_equation_string(spec.mechanism), "\n" => "\\n")
write(path, join(lines, "\n") * "\n")
println("rewrote lines ", i, " and ", i + 1)
```

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T7_golden.log 900 --project=. -e 'using TestEnv; TestEnv.activate(); include("'$S'/phase45/probes/P4_T7_golden.jl")'`
then `git diff --stat test/reference/allosteric_golden_reference.txt`.
Expected: `rewrote lines 17 and 18` (or wherever the LDH header sits), and the diff changes
exactly those two lines; the new header line reads `### LDH i-state NonequalAI 6-group,
steady-state NAD group`, and the new string's parameter header lists 10 fitted parameters, the
last `L`, then `Keq, E_total`.

- [ ] **Step 5: Update the canonical golden fixture**

In `test/fixtures/canonical_text_golden.txt` (read by `test/test_reaction_graph.jl`),
replace `# Numerator: redundant SS-bind` with
`# Numerator: redundant SS-bind, steady-state feeder` and, in that entry,
`4: ES1 + S2 ⇌ ES1S2 | group 4` with `4: ES1 + S2 <--> ES1S2 | group 4`. Replace
`# LDH i-state NonequalAI 6-group` (the entry whose next lines start `substrates: NADH,
Pyruvate` and `forms: E, ELactate, ELactateNAD, ELactateNADH, ENAD, ENADPyruvate, ...`) with
`# LDH i-state NonequalAI 6-group, steady-state NAD group`, and in that entry:

```
2: E + NAD ⇌ ENAD | group 2 :: EqualAI
5: ELactate + NAD ⇌ ELactateNAD | group 2 :: EqualAI
12: EPyruvate + NAD ⇌ ENADPyruvate | group 2 :: EqualAI
```

with

```
2: E + NAD <--> ENAD | group 2 :: EqualAI
5: ELactate + NAD <--> ELactateNAD | group 2 :: EqualAI
12: EPyruvate + NAD <--> ENADPyruvate | group 2 :: EqualAI
```

The edges keep their order and orientation: a binding is stored as its uptake whatever its
flag, and the flag only breaks ties between edges joining the same two forms.

- [ ] **Step 6: Run the tests to verify they pass**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T7_green.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_admission.jl"); include("test/test_rate_eq_derivation.jl"); include("test/test_allosteric_golden.jl"); include("test/test_reaction_graph.jl")'`
Expected: all four files pass, pristine; in test_rate_eq_derivation.jl the per-spec testsets of
both variants pass (structure, constraints with 7 and 10 independent parameters, the
mass-action oracle, `rate_equation` with `allocs == 0` and `t < 120e-9`, the string checks),
and "Enzyme Derivation Tests" still counts 42 specs.

- [ ] **Step 7: Run the full suite**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T7_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`
Expected: `exit 0`, tests passed, pristine, Aqua and JET green; report `maxrss_kb`.

- [ ] **Step 8: Commit**

Write the measured `maxrss_kb` of Step 7 in place of `<kb>`.

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
git add test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_admission.jl \
    test/reference/allosteric_golden_reference.txt test/fixtures/canonical_text_golden.txt
git commit -m "$(cat <<'EOF'
Replace the two trap specs by well-posed variants

"Numerator: redundant SS-bind" traps substrate-free enzyme and "LDH i-state
NonequalAI 6-group" lactate-bound enzyme, so the guard rejects both. Each
becomes a rejection test, and its variant takes its place: the feeding
binding, or the LDH spec's :EqualAI NAD group, steady-state. Before the pins
changed, both variants matched the exact steady state to 1e-20 and passed
every face with both oracles agreeing; they fit 7 and 10 parameters
(Haldane 1, Wegscheider 1 and 2). The golden reference changes in the LDH
variant's two lines only. Full-suite maxrss: <kb> KB.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
)"
```

---

### Evidence for Tasks 8-10

The expected verdicts were measured, not inferred, where the inputs did not already give
them; every probe and log is in `$S/phase45/plan/probes/`.
- `p3_verdicts.jl` ran research §1's rules (`rules/lib.jl`, with B v3's additions as
  `p2_verdicts.jl` wrote them) and both research oracles, face by face, on every fixture of these
  tasks that RES §3.2 and Task 0 do not cover (`p3_verdicts.log`; the tied-exchange J
  variants `p3_verdicts_x.log`; Task 10's activator route `p3_verdicts_t10.log` and the
  `:NonequalAI` fixes of the rejected `:OnlyI` shapes `p3_verdicts_t10b.log`).
- `p3_scan.jl` ran the same rules on every mechanism macro block of the S3g test files
  that phase 3's check accepts (`p3_scan.log`: 40 of 850 rejected, each listed by file,
  line and testset); Task 8 Step 3 converts each.
- `p3_class5.jl` ran class 5 (Task 6) on the derived laws of every fixture of Tasks 1-7
  (`p3_class5.log`) and of the S3g test files (`p3_class5_tests.log`): the laws Task 8
  stops deriving.
- `p3_moves.jl` ran the flip move with class 3 inside its predicate against the move with
  the check at the exits only (R4.5 (ii), `extras/pkg_moves.jl`) on the census bi-bi and
  ping-pong parents, which the 2026-10-14 outage had lost (`p3_moves_bibi.log`,
  `p3_moves_pingpong.log`: one parent differs in each, and `p3_parents_*.log` find both
  parents admitted); `p3_split.jl` counted the split testset's parents the guard admits
  (`p3_split.log`).
- `p3_tied.jl`, `p3_ssmirror.jl` and `p3_degenerate.jl` ran phase 3's package on the
  variants Tasks 8 and 10 substitute: the tied exchanges' orientation and names, the copy
  rule's steady-state mirror and its variants, and the degeneracy conditions on each
  conformation's law of Task 10's activator route.
- `p3c_sink.jl` built the three half-2 ping-pong variants Task 8 substitutes on the S3f
  snapshot (`--project=$S/phase45/plan/probes/env_p2`): each has a rate of exactly 0.0 with
  both products at zero and a positive peak turnover, and the dead variant's fitted names
  map onto its plain twin's (`p3c_sink.log`).
- `joint_r1_mech.jl` (log `joint_r1_mech.log`, rows `joint_r1_mech.tsv`) ran the research
  rules, B v3 and both research oracles (two points, 128 corners, five orders) on Task 10's
  joint-R1 mechanism (research name jointR1_a): fewest slow transitions (none vanishing,
  the active conformation's, the inactive one's, both) = (4, 4, 4, 3) at zero products and
  at P = 0, where B v3 rejects by joint R1 alone; every other check and face passes; on
  those two faces the RE-speed oracle agrees with class 3 and the ε-path fails (research
  RE-first direction-dependent; `vi/vor.jl` converges to 0.0056 against the deletion law's
  0), so the rejection is right. `onlyi_active_faults.jl` gave Task 10's three
  active-conformation faults: `G1` inside and `G1` cycle with the ε-path `nopath`, and an
  island whose ε-path reaches a wrong limit.
- The probes ran research code, not this plan's: their ε-path reads two orders, Task 1's
  oracle five (P4-24a). Where a research verdict mixes the orders (an island going
  elsewhere RE-first and converging ε-first), the tests accept either of the oracle's two
  non-converging verdicts, `:direction_dependent` or `:wrong_limit`, as Task 4's island
  tests do.

---

### Task 8: Wire the guard (RES §4 row 7)

The macros, the parser and every derivation entry throw the guard's message; the enumerator's
exits and identify's base tier drop what it rejects; class 3 runs inside the flip move's
predicate; class 5 runs where the law is derived. The old `:OnlyA` guard, the split's
pre-check and balancers, and the derive-time RE-cycle error go. The one-flip check stays until
Task 9.

**Files:**
- Modify: `src/admission.jl` (`_check_well_posed`, `_is_well_posed`, `_free_forms_flip_together`,
  `_assert_support`; `_assert_bottomed_segments` deleted; `_admission_message`'s `_ISLAND`
  branch)
- Modify: `src/rate_eq_derivation.jl` (`_num_den_exprs`, both methods; `_compute_alpha` and its
  call; the `_state_rate_polys` docstring)
- Modify: `src/mechanism_enumeration.jl` (`_minimal_flips`, `_expand_re_to_ss`,
  `_expand_split_kinetic_group`, `_split_gain_test` (both); `expand_mechanisms`,
  `init_mechanisms`, `seed_mechanisms`)
- Modify: `src/identify_rate_equation.jl` (`_base_tier`)
- Modify: `src/thermo.jl` (delete `_onlya_haldane_violation(am)` and `_onlya_violation_message`;
  `_onlya_offenders`' docstring)
- Modify: `src/names.jl` (`_edge_constants`' docstring)
- Modify: `src/dsl.jl` (both macro docstrings: what the macros reject and where class 5
  surfaces)
- Modify: `ROADMAP.md` (rewrite "`BigInt` arithmetic inflates the GC's live-byte count" and
  "The `:OnlyA` guard misjudges groups that mix a binding and a chemistry step"),
  `docs/src/roadmap.md` (P4-3)
- Test: `test/test_admission.jl`, `test/test_types.jl`, `test/test_dsl.jl`,
  `test/test_rate_eq_derivation.jl`, `test/allosteric_ground_truth.jl`,
  `test/test_oracles.jl`, `test/test_mechanism_enumeration.jl`,
  `test/test_identify_rate_equation.jl`, `test/test_reaction_graph.jl` (Step 7)
- Create (scratch): `$S/phase45/scripts/P4_terter.jl`, `$S/phase45/scripts/P4_bibi_kids.jl`,
  `$S/phase45/scripts/P4_census_explain.jl`;
  census `census/env_guard`, `census/P4_base`, `census/P4_T8`

**Interfaces:**
- Consumes: Tasks 2-6's `_guard(m)`, `_guard_message(mech, v)`, `_faces_pass(mech, re)`,
  `_support_verdict(m, branches)`, `_admission(am)`, `_admission_message`, `_Verdict`,
  `_WELL_POSED`, the `_Fault` values; Task 1's `@_testhelper_unguarded`,
  `_testhelper_respeed_verdicts(m; rng, npoints, law)`, `_testhelper_epath_verdict(m; rng,
  face, corners, precision)`, `_testhelper_single_subunit(m)` (concrete methods),
  `_testhelper_law_value(m, params, concs)`; phase 3's `_minimal_flips`, `_expand_re_to_ss`,
  `_split_gain_test`, `_expand_split_kinetic_group`, `_expand_to_allosteric`,
  `_expand_change_allo_state`, `expand_mechanisms`, `init_mechanisms`, `seed_mechanisms`,
  `_base_tier`, the parser (`_parse_mechanism`), `show`, `_mwc_state_polys`,
  `_raw_symbolic_rate_polys`, `poly_one`, `_step_text(m, k)` (names.jl).
- Produces:
  ```julia
  _check_well_posed(m::Union{Mechanism, AllostericMechanism})::Nothing  # _guard_message
  _is_well_posed(m::Mechanism)::Bool                     # _guard(m).fault === _NO_FAULT
  _is_well_posed(am::AllostericMechanism)::Bool          # and the one flip, until Task 9
  _free_forms_flip_together(am::AllostericMechanism)::Bool   # Task 9 deletes it
  _assert_support(m::Union{Mechanism, AllostericMechanism},
                  branches::Vector{Tuple{POLY, POLY}})::Nothing   # class 5's throw
  _minimal_flips(admissible, m::Union{Mechanism, AllostericMechanism},
                 eligible = _ -> true)::Vector{BitVector}
  _split_gain_test(m::Mechanism, units, flags_of)                 # -> gains
  _split_gain_test(am::AllostericMechanism, units, flags_of)      # -> gains
  ```
  Kept for Task 12: `_onlya_offenders`, `_onlya_cone_offenders`,
  `_has_strict_positive_combination` (the tagging and change-tag filters still read them).

**Gate:**
- The new and converted testsets below pass, each rejection with its fault, face, named step
  and oracle verdict.
- Census `census/P4_T8` against `census/P4_base` (this task's starting tree, the branch point
  plus Tasks 1-7's unwired code and spec variants): compare.jl reports differences only in
  child sets, seeds and the law sample; `P4_census_explain.jl` rebuilds every child, seed and
  law only in `P4_base` and finds a guard fault on each, and finds every child only in
  `P4_T8` a guard-passing flip of its parent. Expected by class (R4.5, R4.7): ping-pong
  tagging children `_ISLAND` (pingpong, prefix, terbi parents), mixed `:OnlyA` groups
  `_NO_LIMIT` (bibi level 2, 64 in the census), class-3 rejections of split and change-tag
  children, the pingpong seeds that tag both half-reactions `_ISLAND`; init sets equal; split
  children equal but for exit drops; about one parent per reaction swaps a class-3 flip for
  superset flips.
- Ter-ter `init_mechanisms` under 150 s, its time before and after this task measured and
  written in the commit body; its count re-pinned to the measured value after the census
  explains the drop (RES: 249,559 at the exit with the check at the exits only, 1,296 class-3
  drops; the check inside the seed variants' flips may add superset variants).
- Bi-bi init golden and the init tallies (TME "bi-bi init_mechanisms structural golden", the
  tally testsets) unchanged; TID "_base_tier expands degenerate seeds instead of fitting them",
  "_base_tier replaces a degenerate seed by the flips that cure it" and "_base_tier records a
  degenerate seed's expansion error" unchanged.
- `rate_equation` gates (allocations 0, under 120 ns) and the compile budgets hold, the
  budgets measured in Step 6 before the full run; the full suite passes, pristine, with Aqua
  and JET. R52's three figures (full-suite maxrss, ter-ter init, `test_mechanism_enumeration.jl`
  focused) go in the commit body beside the branch point's.
- Class 5's throw (`_assert_support`) has no public path: at every derivation entry the
  guard runs first, and its class 3 rejects every known 0/0 law (P4-15). Task 6's tests of
  the internal function are its whole coverage; identify's pass 1 would record it as a
  `FitFailure` through the try/catch that TID "derivation guard (filter >
  MAX_RATE_EQUATION_TERMS) → FitFailure, not crash" exercises.

- [ ] **Step 0: Check the names this task reads at the branch point.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  grep -n "^function _minimal_flips\|^function _split_gain_test\|^function _expand_split_kinetic_group\|^function _expand_re_to_ss\|^function expand_mechanisms\|^function init_mechanisms\|^function seed_mechanisms\|^function _expand_to_allosteric\|^function _expand_change_allo_state" src/mechanism_enumeration.jl
  grep -n "^function _base_tier" src/identify_rate_equation.jl
  grep -n "_onlya_haldane_violation\|_onlya_violation_message\|_assert_bottomed_segments\|_bottomless_segment(\|_check_well_posed" src/*.jl
  grep -n "_raw_mechanism(" src/*.jl
  ```
  Expected (as at 08c4043): one `_minimal_flips(admissible, m::Mechanism,
  eligible = _ -> true)` (the seed variants call it on a plain base);
  `_split_gain_test(m::Mechanism, units, flags_of)` and `(am::AllostericMechanism, units,
  flags_of)`, the latter reading `_onlya_offenders`; `_bottomless_segment(` called by
  `_assert_bottomed_segments`, `_minimal_flips` and `_expand_split_kinetic_group` only;
  `_onlya_haldane_violation(am)` and `_onlya_violation_message` called only by
  `_check_well_posed`; `_check_well_posed` called by `_well_posed` (admission.jl), the three
  exits (mechanism_enumeration.jl) and `_base_tier` (its dsl.jl hits are docstrings: the
  macros and the parser, `_parse_mechanism`, reach it through `_checked_mechanism`, which
  calls `_well_posed`); `_raw_mechanism(` called by the macros and by the parser. Step 2's
  testset "the parser rejects a mechanism the guard rejects, naming its step" calls it as
  `EnzymeRates._parse_mechanism(text)`. If a signature differs, call it as the branch point
  defines it and say so in the commit body.

- [ ] **Step 1: Record the starting tree's census and ter-ter time.**

  Create `$S/phase45/scripts/P4_terter.jl`:
  ```julia
  # ABOUTME: Times init_mechanisms on the ter-ter reaction, counts its mechanisms and those
  # ABOUTME: the guard rejects; run on the tree before and after the guard is wired.
  using EnzymeRates
  small = @enzyme_reaction begin
      substrates: S[C]
      products: P[C]
  end
  rxn = @enzyme_reaction begin
      substrates: A[C], B[N], D[X]
      products: P[C], Q[N], R[X]
  end
  EnzymeRates.init_mechanisms(small)
  t = @elapsed ms = EnzymeRates.init_mechanisms(rxn)
  rejected = count(m -> EnzymeRates._guard(m).fault !== EnzymeRates._NO_FAULT, ms)
  println("ter-ter init: ", length(ms), " mechanisms in ", round(t; digits = 1),
          " s; the guard rejects ", rejected)
  ```
  and `$S/phase45/scripts/P4_bibi_kids.jl`, the value check of TME's bi-bi raw-input
  aggregate (Step 7, kind 4):
  ```julia
  # ABOUTME: Counts the bi-bi seeds' first-round children, all and allosteric, and those the
  # ABOUTME: guard rejects; run before and after the guard is wired, for TME's bi-bi pins.
  using EnzymeRates
  allosteric(m) = m isa EnzymeRates.AllostericMechanism
  rxn = @enzyme_reaction begin
      substrates: A[C], B[N]
      products: P[C], Q[N]
  end
  kids = EnzymeRates.expand_mechanisms(EnzymeRates.init_mechanisms(rxn), rxn)
  bad = [k for k in kids if EnzymeRates._guard(k).fault !== EnzymeRates._NO_FAULT]
  println("bi-bi children: ", length(kids), " (", count(allosteric, kids),
          " allosteric); the guard rejects ", length(bad), " (", count(allosteric, bad),
          " allosteric)")
  ```
  Then, from the worktree at the Task 7 commit, with nothing edited yet:
  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && git status --short
  test -d $S/phase3/census/env_guard || $S/phase3/jl.sh $S/phase45/logs/P4_T8_mkenv.log 1800 --heap-size-hint=2500M $S/phase3/census/mkenv.jl $S/phase45/guard_wt $S/phase3/census/env_guard
  $S/phase3/jl.sh $S/phase45/logs/P4_T8_census_base.log 1800 --project=$S/phase3/census/env_guard --heap-size-hint=2500M $S/phase3/census/corpus.jl $S/phase45/guard_wt $S/phase3/census/P4_base $S/phase3/census/rename_phase3.toml
  $S/phase3/jl.sh $S/phase45/logs/P4_T8_terter_base.log 900 --project=. --heap-size-hint=2500M $S/phase45/scripts/P4_terter.jl
  $S/phase3/jl.sh $S/phase45/logs/P4_T8_bibi_kids_base.log 900 --project=. --heap-size-hint=2500M $S/phase45/scripts/P4_bibi_kids.jl
  ```
  Expected: a clean tree; the corpus log ends `exit 0` after about 110 s (phase 3's final
  census, `census/logs/replay_TFW.log`, took 112 s) and `$S/phase3/census/P4_base` holds
  `mechs.tsv`, `reactions.tsv`, `init.tsv`, `parents.tsv`, `children.tsv`, `seeds.tsv`,
  `laws.tsv`, `meta.toml`; the ter-ter log reads `ter-ter init: 250855 mechanisms
  in <t0> s; the guard rejects <R0>`, `<t0>` near phase 3's final ter-ter figure (R52's
  gate at S3i: at most 31 s) and `<R0>` near 1,296 (RES §1.10: the class-3 drops of the
  exits). Record `<t0>` and `<R0>`: Step 8's re-pin reads them. The bi-bi log reads
  `bi-bi children: <n0> (<a0> allosteric); the guard rejects <K0> (<A0> allosteric)`, with
  `<n0>` and `<a0>` the counts TME's raw-input testset pins at the branch point (2,831 and
  1,833 at 08c4043, as at S3g); record all four for Step 7.

- [ ] **Step 2: Write the new tests.** Append to `test/test_admission.jl`:
  ```julia
  # ── The guard wired: macros, exits, flip move and parser ──────────────────────

  @testset "the macros reject what the guard rejects when the expanded code runs" begin
      # One mechanism per class: no cycle turns the reaction over (class 1); a cycle of
      # rapid-equilibrium steps does (class 2); a trap at zero products (class 3); an
      # :OnlyA binding whose two sides untagged steps join (class 4). Each `@test_throws`
      # holds the macro call itself, not an `eval`: had a macro thrown at expansion, this
      # file would not load (R27).
      no_cycle = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E(P)
              E(P) <--> F + P
          end
      end
      re_cycle = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) ⇌ E(P)
              E + P ⇌ E(P)
              E(S) <--> F(P)
              F(P) <--> E + P
          end
      end
      trap = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) + B <--> E(A, B)
              E(A, B) <--> E(P, Q)
              E(P) + Q ⇌ E(P, Q)
              E(Q) + P ⇌ E(P, Q)
              E + P <--> E(P)
              E + Q ⇌ E(Q)
          end
      end
      no_limit = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)      :: OnlyA
              E(S) <--> E(P)    :: EqualAI
              E + P ⇌ E(P)      :: EqualAI
          end
      end
      faults = [EnzymeRates._guard(m).fault for m in (no_cycle, re_cycle, trap, no_limit)]
      @test faults == [EnzymeRates._NO_CYCLE, EnzymeRates._RE_TURNOVER, EnzymeRates._TRAP,
                       EnzymeRates._NO_LIMIT]
      thrown(m) = ErrorException(EnzymeRates._guard_message(m, EnzymeRates._guard(m)))
      @test_throws thrown(no_cycle) @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E(P)
              E(P) <--> F + P
          end
      end
      @test_throws thrown(re_cycle) @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) ⇌ E(P)
              E + P ⇌ E(P)
              E(S) <--> F(P)
              F(P) <--> E + P
          end
      end
      @test_throws thrown(trap) @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) + B <--> E(A, B)
              E(A, B) <--> E(P, Q)
              E(P) + Q ⇌ E(P, Q)
              E(Q) + P ⇌ E(P, Q)
              E + P <--> E(P)
              E + Q ⇌ E(Q)
          end
      end
      @test_throws thrown(no_limit) @allosteric_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)      :: OnlyA
              E(S) <--> E(P)    :: EqualAI
              E + P ⇌ E(P)      :: EqualAI
          end
      end
  end

  @testset "the parser rejects a mechanism the guard rejects, naming its step" begin
      # Review Focus 3: phase 3 wrote allosteric ping-pongs with both half-reactions :OnlyA
      # (an identify CSV row, `show` text). Read back, such text raises the guard's message,
      # not a parse error, and the message names the half-reaction whose untagging fixes
      # it (P4-29b).
      both = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)                                          :: EqualAI
              E(A) <--> E(P; residual = A - P)                      :: OnlyA
              E(; residual = A - P) + P ⇌ E(P; residual = A - P)    :: EqualAI
              E(; residual = A - P) + B ⇌ E(B; residual = A - P)    :: EqualAI
              E(B; residual = A - P) ⇌ E(Q)                         :: OnlyA
              E + Q ⇌ E(Q)                                          :: EqualAI
          end
      end
      v = EnzymeRates._guard(both)
      @test v.fault === EnzymeRates._ISLAND
      msg = EnzymeRates._guard_message(both, v)
      @test occursin("Tag the group of E(A) <--> E(P; residual = A - P) :EqualAI or " *
                     ":NonequalAI.", msg)
      @test occursin("is not supported yet (see the Roadmap page)", msg)
      @test_throws ErrorException(msg) EnzymeRates._parse_mechanism(sprint(show, both))
  end
  ```

  Append to `test/test_mechanism_enumeration.jl` (CLAUDE.md: the tests of
  `expand_mechanisms` and the moves live there, under its three rules):
  ```julia
  @testset "enumeration drops the children the guard rejects" begin
      # A plain ping-pong (census 3fdbe6e1dc001bb7) whose A binding is fused into the
      # product-bound covalent form. Its one rapid-equilibrium group flips to steady state:
      # one plain child. Today's tagging move tags both half-reactions' fused bindings :OnlyA
      # in each of its five children, so in the inactive conformation the covalent forms join
      # free E only through tagged chemistry: an island (`_ISLAND`), whose limit depends on
      # how fast the tagged constants vanish. The exit drops the five and keeps the flip.
      rxn = @enzyme_reaction begin
          substrates: A[CX], B[N]
          products: P[C], Q[NX]
          allosteric_regulators: R
          oligomeric_state: 2
          shared_catalytic_site: (A, P)
      end
      parent = EnzymeRates.Mechanism(@enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(P; residual = A - P)
              E(; residual = A - P) + P <--> E(P; residual = A - P)
              E(; residual = A - P) + B <--> E(Q)
              E + Q <--> E(Q)
          end
      end)
      flip = EnzymeRates.Mechanism(@enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A <--> E(P; residual = A - P)
              E(; residual = A - P) + P <--> E(P; residual = A - P)
              E(; residual = A - P) + B <--> E(Q)
              E + Q <--> E(Q)
          end
      end)
      # Untyped: inside `T[...]` the parser reads `begin` as an index keyword.
      tagged = [
          @_testhelper_unguarded(@allosteric_mechanism begin
              substrates: A, B
              products: P, Q
              catalytic_multiplicity: 2
              steps: begin
                  E + A ⇌ E(P; residual = A - P)                        :: OnlyA
                  E(; residual = A - P) + P <--> E(P; residual = A - P) :: EqualAI
                  E(; residual = A - P) + B <--> E(Q)                   :: OnlyA
                  E + Q <--> E(Q)                                       :: OnlyA
              end
          end),
          @_testhelper_unguarded(@allosteric_mechanism begin
              substrates: A, B
              products: P, Q
              catalytic_multiplicity: 2
              steps: begin
                  E + A ⇌ E(P; residual = A - P)                        :: OnlyA
                  E(; residual = A - P) + P <--> E(P; residual = A - P) :: OnlyA
                  E(; residual = A - P) + B <--> E(Q)                   :: OnlyA
                  E + Q <--> E(Q)                                       :: EqualAI
              end
          end),
          @_testhelper_unguarded(@allosteric_mechanism begin
              substrates: A, B
              products: P, Q
              catalytic_multiplicity: 2
              steps: begin
                  E + A ⇌ E(P; residual = A - P)                        :: OnlyA
                  E(; residual = A - P) + P <--> E(P; residual = A - P) :: OnlyA
                  E(; residual = A - P) + B <--> E(Q)                   :: OnlyA
                  E + Q <--> E(Q)                                       :: OnlyA
              end
          end),
          @_testhelper_unguarded(@allosteric_mechanism begin
              substrates: A, B
              products: P, Q
              allosteric_regulators: R
              catalytic_multiplicity: 2
              steps: begin
                  E + A ⇌ E(P; residual = A - P)                        :: OnlyA
                  E(; residual = A - P) + P <--> E(P; residual = A - P) :: EqualAI
                  E(; residual = A - P) + B <--> E(Q)                   :: OnlyA
                  E + Q <--> E(Q)                                       :: EqualAI
              end
              regulatory_site(2): R::OnlyA
          end),
          @_testhelper_unguarded(@allosteric_mechanism begin
              substrates: A, B
              products: P, Q
              allosteric_regulators: R
              catalytic_multiplicity: 2
              steps: begin
                  E + A ⇌ E(P; residual = A - P)                        :: OnlyA
                  E(; residual = A - P) + P <--> E(P; residual = A - P) :: EqualAI
                  E(; residual = A - P) + B <--> E(Q)                   :: OnlyA
                  E + Q <--> E(Q)                                       :: EqualAI
              end
              regulatory_site(2): R::OnlyI
          end)]
      moved = EnzymeRates._expand_to_allosteric(parent, rxn)
      @test length(moved) == 5
      @test Set(moved) == Set(tagged)
      @test all(t -> EnzymeRates._guard(t).fault === EnzymeRates._ISLAND, tagged)
      # The ε-path oracle, on one subunit (P4-24a), finds no single limit of the island:
      # direction-dependent or a wrong limit, as test_admission.jl's island tests read it.
      @test _testhelper_epath_verdict(_testhelper_single_subunit(tagged[4]);
                                      rng = Random.MersenneTwister(81)) in
            (:direction_dependent, :wrong_limit)
      @test EnzymeRates._expand_re_to_ss(parent) == [flip]
      kids = EnzymeRates.expand_mechanisms([parent], rxn)
      @test length(kids) == 1
      @test Set(kids) == Set([flip])
  end

  @testset "class 3 inside the flip move emits the superset flips" begin
      # Random-order merged bi-bi with the abortive complex E(A, P) (census
      # 0988ab459f638946, the bi-bi parent of R4.5). Flipping the A group alone gains a
      # segment but leaves a face whose steady state moves with the relative speeds of the
      # rapid-equilibrium steps, a class-3 rejection (`_LIMIT_SPEED`) the RE-speed oracle
      # reads as a pattern. Inside the flip move's predicate (`_faces_pass`) that set fails,
      # so the move extends it: the A group with each of the other three rapid-equilibrium
      # groups. The other two minimal sets pass and stay: five children where the check at
      # the exit alone kept two.
      parent = EnzymeRates.Mechanism(@enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              (E + A ⇌ E(A), E(B) + A ⇌ E(P, Q))
              (E + B <--> E(B), E(A) + B <--> E(P, Q))
              E + P ⇌ E(P)
              E(A) + P ⇌ E(A, P)
              E(P) + A ⇌ E(A, P)
              E(P) + Q <--> E(P, Q)
          end
      end)
      a_alone = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              (E + A <--> E(A), E(B) + A <--> E(P, Q))
              (E + B <--> E(B), E(A) + B <--> E(P, Q))
              E + P ⇌ E(P)
              E(A) + P ⇌ E(A, P)
              E(P) + A ⇌ E(A, P)
              E(P) + Q <--> E(P, Q)
          end
      end
      v = EnzymeRates._guard(a_alone)
      @test v.fault !== EnzymeRates._NO_FAULT && v.limit === EnzymeRates._LIMIT_SPEED
      @test _testhelper_respeed_verdicts(a_alone; rng = Random.MersenneTwister(91),
                                         law = :limit)[v.face] === :pattern
      flags = EnzymeRates._flags(parent)
      flags[[_testhelper_edge(parent, "E + A ⇌ EA"),
             _testhelper_edge(parent, "EB + A ⇌ EPQ")]] .= false
      @test !EnzymeRates._faces_pass(parent, flags)
      expected = EnzymeRates.Mechanism.([
          (@enzyme_mechanism begin      # A and P@E
              substrates: A, B
              products: P, Q
              steps: begin
                  (E + A <--> E(A), E(B) + A <--> E(P, Q))
                  (E + B <--> E(B), E(A) + B <--> E(P, Q))
                  E + P <--> E(P)
                  E(A) + P ⇌ E(A, P)
                  E(P) + A ⇌ E(A, P)
                  E(P) + Q <--> E(P, Q)
              end
          end),
          (@enzyme_mechanism begin      # A and P@E(A)
              substrates: A, B
              products: P, Q
              steps: begin
                  (E + A <--> E(A), E(B) + A <--> E(P, Q))
                  (E + B <--> E(B), E(A) + B <--> E(P, Q))
                  E + P ⇌ E(P)
                  E(A) + P <--> E(A, P)
                  E(P) + A ⇌ E(A, P)
                  E(P) + Q <--> E(P, Q)
              end
          end),
          (@enzyme_mechanism begin      # A and A@E(P)
              substrates: A, B
              products: P, Q
              steps: begin
                  (E + A <--> E(A), E(B) + A <--> E(P, Q))
                  (E + B <--> E(B), E(A) + B <--> E(P, Q))
                  E + P ⇌ E(P)
                  E(A) + P ⇌ E(A, P)
                  E(P) + A <--> E(A, P)
                  E(P) + Q <--> E(P, Q)
              end
          end),
          (@enzyme_mechanism begin      # P@E and P@E(A)
              substrates: A, B
              products: P, Q
              steps: begin
                  (E + A ⇌ E(A), E(B) + A ⇌ E(P, Q))
                  (E + B <--> E(B), E(A) + B <--> E(P, Q))
                  E + P <--> E(P)
                  E(A) + P <--> E(A, P)
                  E(P) + A ⇌ E(A, P)
                  E(P) + Q <--> E(P, Q)
              end
          end),
          (@enzyme_mechanism begin      # P@E and A@E(P)
              substrates: A, B
              products: P, Q
              steps: begin
                  (E + A ⇌ E(A), E(B) + A ⇌ E(P, Q))
                  (E + B <--> E(B), E(A) + B <--> E(P, Q))
                  E + P <--> E(P)
                  E(A) + P ⇌ E(A, P)
                  E(P) + A <--> E(A, P)
                  E(P) + Q <--> E(P, Q)
              end
          end)])
      children = EnzymeRates._expand_re_to_ss(parent)
      @test length(children) == 5
      @test Set(children) == Set(expected)
      @test a_alone ∉ children
  end
  ```


- [ ] **Step 3: Convert the tests whose verdicts or messages change.** Each conversion keeps
  its subject; the probe log `p3_verdicts.log` holds the measured verdict of every fixture
  below (research rules and both research oracles).

  (a) `test/test_types.jl`, testset "the macro rejects a rapid-equilibrium segment with no
  bottom form": replace its body (from its first comment through the closing `end` of the
  testset; the raw-input part, `inh_raw` and `inh_m`, goes as well) with:
  ```julia
      @testset "the macro rejects a rapid-equilibrium segment with no bottom form" begin
          # Random product release from E(P, Q) at rapid equilibrium, but release from E(P)
          # and E(Q) at steady state: the RE segment {E(P), E(P, Q), E(Q)} has weights
          # P·Q : Kq·P : Kp·Q, all zero at P = Q = 0. The RE approximation carries no
          # parameter for how E(P, Q) splits between E(P) and E(Q) there, so the rate is
          # 0/0 at zero products: class 3 rejects that face (`_BOTTOMLESS`), and the RE-speed
          # oracle sees the steady state there move with the release speeds.
          random_release = @_testhelper_unguarded @enzyme_mechanism begin
              substrates: S
              products: P, Q
              steps: begin
                  E + S ⇌ E(S)
                  E + P <--> E(P)
                  E + Q <--> E(Q)
                  E(P) + Q ⇌ E(P, Q)
                  E(Q) + P ⇌ E(P, Q)
                  E(S) <--> E(P, Q)
              end
          end
          v = ER._guard(random_release)
          @test (v.fault, v.face) == (ER._BOTTOMLESS, 2)
          msg = ER._guard_message(random_release, v)
          @test occursin("has no form free of the zeroed metabolites", msg)
          @test _testhelper_respeed_verdicts(random_release; rng = Random.MersenneTwister(1),
                                             law = :limit)[2] === :pattern
          @test_throws ErrorException(msg) @enzyme_mechanism begin
              substrates: S
              products: P, Q
              steps: begin
                  E + S ⇌ E(S)
                  E + P <--> E(P)
                  E + Q <--> E(Q)
                  E(P) + Q ⇌ E(P, Q)
                  E(Q) + P ⇌ E(P, Q)
                  E(S) <--> E(P, Q)
              end
          end

          # Mirror on the substrate side: the segment {E(A), E(B), E(A, B)} is 0/0 at
          # A = B = 0 (face 3).
          mirror = @_testhelper_unguarded @enzyme_mechanism begin
              substrates: A, B
              products: P
              steps: begin
                  E + A <--> E(A)
                  E + B <--> E(B)
                  E(A) + B ⇌ E(A, B)
                  E(B) + A ⇌ E(A, B)
                  E + P ⇌ E(P)
                  E(A, B) <--> E(P)
              end
          end
          v = ER._guard(mirror)
          @test (v.fault, v.face) == (ER._BOTTOMLESS, 3)
          @test _testhelper_respeed_verdicts(mirror; rng = Random.MersenneTwister(2),
                                             law = :limit)[3] === :pattern
          @test_throws ErrorException(ER._guard_message(mirror, v)) @enzyme_mechanism begin
              substrates: A, B
              products: P
              steps: begin
                  E + A <--> E(A)
                  E + B <--> E(B)
                  E(A) + B ⇌ E(A, B)
                  E(B) + A ⇌ E(A, B)
                  E + P ⇌ E(P)
                  E(A, B) <--> E(P)
              end
          end

          # The same segment with A's modifier copy binding E(B) at rapid equilibrium:
          # E(A::Mod, B) joins it, but every form still carries A or B (the copy's
          # concentration is A's), so the segment stays bottomless, whichever order the
          # groups are written in.
          copy_last = @_testhelper_unguarded @enzyme_mechanism begin
              substrates: A, B
              products: P
              modifiers: A
              steps: begin
                  E + A <--> E(A)
                  E + B <--> E(B)
                  E(A) + B ⇌ E(A, B)
                  E(B) + A ⇌ E(A, B)
                  E(B) + A::Mod ⇌ E(A::Mod, B)
                  E + P ⇌ E(P)
                  E(A, B) <--> E(P)
              end
          end
          copy_first = @_testhelper_unguarded @enzyme_mechanism begin
              substrates: A, B
              products: P
              modifiers: A
              steps: begin
                  E(B) + A::Mod ⇌ E(A::Mod, B)
                  E(A, B) <--> E(P)
                  E + P ⇌ E(P)
                  E(B) + A ⇌ E(A, B)
                  E(A) + B ⇌ E(A, B)
                  E + B <--> E(B)
                  E + A <--> E(A)
              end
          end
          v = ER._guard(copy_last)
          @test (v.fault, v.face) == (ER._BOTTOMLESS, 3)
          @test ER._guard(copy_first) == v
          @test_throws ErrorException(ER._guard_message(copy_last, v)) @enzyme_mechanism begin
              substrates: A, B
              products: P
              modifiers: A
              steps: begin
                  E + A <--> E(A)
                  E + B <--> E(B)
                  E(A) + B ⇌ E(A, B)
                  E(B) + A ⇌ E(A, B)
                  E(B) + A::Mod ⇌ E(A::Mod, B)
                  E + P ⇌ E(P)
                  E(A, B) <--> E(P)
              end
          end

          # Ordered release: the segment {E(P), E(P, Q)} has weights Kq : Q, so E(P) is its
          # bottom form and the rate is defined at zero products; every face agrees with
          # the RE-speed oracle.
          m_ordered = @enzyme_mechanism begin
              substrates: S
              products: P, Q
              steps: begin
                  E + S ⇌ E(S)
                  E + P <--> E(P)
                  E(P) + Q ⇌ E(P, Q)
                  E(S) <--> E(P, Q)
              end
          end
          @test m_ordered isa EnzymeMechanism
          @test all(==(:agree),
                    _testhelper_respeed_verdicts(m_ordered; rng = Random.MersenneTwister(3)))

          # Mixed abortive complex: the segment {E(S), E(P), E(S, P)} has weights
          # K·S : K′·P : S·P, zero only at S = P = 0, a corner outside the guard's faces
          # where no turnover is possible anyway. Accepted.
          m_abortive = @enzyme_mechanism begin
              substrates: S
              products: P
              steps: begin
                  E + S <--> E(S)
                  E + P <--> E(P)
                  E(S) + P ⇌ E(S, P)
                  E(P) + S ⇌ E(S, P)
                  E(S) <--> E(P)
              end
          end
          @test m_abortive isa EnzymeMechanism
          @test all(==(:agree),
                    _testhelper_respeed_verdicts(m_abortive; rng = Random.MersenneTwister(4)))
      end
  ```
  The oracle reads its own deletion limit (`law = :limit`) for the three rejected mechanisms:
  their derived denominators lose every monomial on the bottomless face, so class 5 refuses
  to derive their laws (`_assert_support`).

  (b) `test/test_types.jl`, testset "the inactive conformation's segments need a bottom form
  too": replace its body with:
  ```julia
      @testset "the inactive conformation's segments need a bottom form too" begin
          # The :OnlyA P group holds the fused binding E(Q) + P ⇌ Estar(S) beside the plain
          # binding E + P ⇌ E(P). In the active conformation that binding joins E(P),
          # E(P, Q) and E(Q) to free E at rapid equilibrium, and E is the segment's bottom
          # form. The inactive conformation lacks it and reaches those forms only through
          # the steady-state chemistry and the steady-state Q binding at E, so its segment
          # {E(P), E(P, Q), E(Q)} holds a product in every form: class 3 finds it bottomless
          # at zero products. The same untagged steps join the two sides of the tagged
          # binding E + P ⇌ E(P), so its affinity cannot vanish there: admission rejects the
          # tag set first (`_NO_LIMIT`), and the ε-path oracle finds no path.
          am = @_testhelper_unguarded @allosteric_mechanism begin
              substrates: S
              products: P, Q
              steps: begin
                  E + S ⇌ E(S)                               :: EqualAI
                  E(S) <--> E(P, Q)                          :: EqualAI
                  E(P, Q) ⇌ E(P) + Q                         :: EqualAI
                  E(P, Q) ⇌ E(Q) + P                         :: EqualAI
                  (E + P ⇌ E(P), E(Q) + P ⇌ Estar(S))        :: OnlyA
                  E + Q <--> E(Q)                            :: EqualAI
              end
          end
          @test _testhelper_face_faults(am)[2] === ER._BOTTOMLESS
          v = ER._guard(am)
          @test v.fault === ER._NO_LIMIT
          @test ER._step_text(am.catalytic, v.edge) == "E + P ⇌ E(P)"
          @test _testhelper_epath_verdict(am; rng = Random.MersenneTwister(5)) === :no_path
          @test_throws ErrorException(ER._guard_message(am, v)) @allosteric_mechanism begin
              substrates: S
              products: P, Q
              steps: begin
                  E + S ⇌ E(S)                               :: EqualAI
                  E(S) <--> E(P, Q)                          :: EqualAI
                  E(P, Q) ⇌ E(P) + Q                         :: EqualAI
                  E(P, Q) ⇌ E(Q) + P                         :: EqualAI
                  (E + P ⇌ E(P), E(Q) + P ⇌ Estar(S))        :: OnlyA
                  E + Q <--> E(Q)                            :: EqualAI
              end
          end
      end
  ```

  (c) `test/test_types.jl` (its helpers read `_onlya_offenders`),
  testset "OnlyA Haldane validator": the moves' `:OnlyA` test (`_onlya_offenders`, until
  Task 12) and admission give every tag set one verdict. In `_testhelper_uni`, replace its
  last two statements (`offenders = ...` and the returned list) with:
  ```julia
          offenders = ER._onlya_offenders(am.catalytic, am.catalytic.group, tags)
          # Admission (`_admission`) gives the same verdict: a tag set the moves' test
          # rejects has no limit (`_NO_LIMIT`).
          v = ER._admission(ER.AllostericMechanism(am.catalytic, tags, 2,
                                                   ER.RegulatorySite[]))
          @test v.fault === (isempty(offenders) ? ER._NO_FAULT : ER._NO_LIMIT)
          [only(firsts[g].consumed) for g in offenders]
  ```
  in `_testhelper_both_bindings_onlya`, replace its last statement with:
  ```julia
          offenders = ER._onlya_offenders(am.catalytic, am.catalytic.group, tags)
          v = ER._admission(ER.AllostericMechanism(am.catalytic, tags, 2,
                                                   ER.RegulatorySite[]))
          @test v.fault === (isempty(offenders) ? ER._NO_FAULT : ER._NO_LIMIT)
          offenders
  ```
  and replace `_testhelper_verdict(t) = ER._onlya_offenders(...)` with:
  ```julia
          function _testhelper_verdict(t)
              offenders = ER._onlya_offenders(bu_am.catalytic, bu_am.catalytic.group, t)
              v = ER._admission(ER.AllostericMechanism(bu_am.catalytic, t, 2,
                                                       ER.RegulatorySite[]))
              @test v.fault === (isempty(offenders) ? ER._NO_FAULT : ER._NO_LIMIT)
              offenders
          end
  ```
  Replace the `lone_onlya` message and its `@test_throws` with:
  ```julia
          # ...and the macro rejects that tagging: untagged steps join the two sides of the
          # :OnlyA binding of A, so its affinity has no limit (`_NO_LIMIT`).
          lone = @_testhelper_unguarded @allosteric_mechanism begin
              substrates: A, B
              products:   P
              catalytic_multiplicity: 2
              steps: begin
                  E + A ⇌ E(A)          :: OnlyA
                  E + B ⇌ E(B)          :: EqualAI
                  E(A) + B ⇌ E(A, B)    :: EqualAI
                  E(B) + A ⇌ E(A, B)    :: EqualAI
                  E(A, B) <--> E(P)     :: EqualAI
                  E(P) ⇌ E + P          :: EqualAI
              end
          end
          v = ER._guard(lone)
          @test v.fault === ER._NO_LIMIT
          @test ER._step_text(lone.catalytic, v.edge) == "E + A ⇌ E(A)"
          @test_throws ErrorException(ER._guard_message(lone, v)) @allosteric_mechanism begin
              substrates: A, B
              products:   P
              catalytic_multiplicity: 2
              steps: begin
                  E + A ⇌ E(A)          :: OnlyA
                  E + B ⇌ E(B)          :: EqualAI
                  E(A) + B ⇌ E(A, B)    :: EqualAI
                  E(B) + A ⇌ E(A, B)    :: EqualAI
                  E(A, B) <--> E(P)     :: EqualAI
                  E(P) ⇌ E + P          :: EqualAI
              end
          end
  ```
  The verdicts the probe measured: uni (S, chemistry, P) tags (OnlyA, EqualAI, EqualAI),
  (EqualAI, EqualAI, OnlyA) and (OnlyA, NonequalAI, EqualAI) `G1`; the other four and both
  RE/SS mixes admitted on every face; the bi-uni square's one-edge, two-edge and
  chemistry-plus-edge tag sets `G1`, the three-edge set admitted.

  (d) `test/test_types.jl`, testset ":OnlyA guard rejects a multi-cycle ter-substrate
  inconsistency": replace `@test_throws ErrorException @allosteric_mechanism begin ... end`
  with:
  ```julia
      # Admission finds untagged steps that join the two sides of an :OnlyA binding
      # (`_NO_LIMIT`), as the moves' Stiemke test finds the coupled system infeasible.
      cube = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B, C
          products:   P
          catalytic_multiplicity: 2
          steps: begin
              E + A ⇌ E(A)             :: OnlyA
              E + B ⇌ E(B)             :: OnlyA
              E + C ⇌ E(C)             :: OnlyA
              E(A) + B ⇌ E(A, B)       :: OnlyA
              E(A) + C ⇌ E(A, C)       :: EqualAI
              E(B) + A ⇌ E(A, B)       :: EqualAI
              E(B) + C ⇌ E(B, C)       :: EqualAI
              E(C) + A ⇌ E(A, C)       :: EqualAI
              E(C) + B ⇌ E(B, C)       :: EqualAI
              E(A, B) + C ⇌ E(A, B, C) :: OnlyA
              E(A, C) + B ⇌ E(A, B, C) :: OnlyA
              E(B, C) + A ⇌ E(A, B, C) :: OnlyA
              E(A, B, C) <--> E(P)     :: OnlyA
              E + P ⇌ E(P)             :: EqualAI
          end
      end
      @test !isempty(ER._onlya_offenders(cube.catalytic, cube.catalytic.group,
                                         ER.cat_allo_states(cube)))
      v = ER._guard(cube)
      @test v.fault === ER._NO_LIMIT
      @test_throws ErrorException(ER._guard_message(cube, v)) @allosteric_mechanism begin
          substrates: A, B, C
          products:   P
          catalytic_multiplicity: 2
          steps: begin
              E + A ⇌ E(A)             :: OnlyA
              E + B ⇌ E(B)             :: OnlyA
              E + C ⇌ E(C)             :: OnlyA
              E(A) + B ⇌ E(A, B)       :: OnlyA
              E(A) + C ⇌ E(A, C)       :: EqualAI
              E(B) + A ⇌ E(A, B)       :: EqualAI
              E(B) + C ⇌ E(B, C)       :: EqualAI
              E(C) + A ⇌ E(A, C)       :: EqualAI
              E(C) + B ⇌ E(B, C)       :: EqualAI
              E(A, B) + C ⇌ E(A, B, C) :: OnlyA
              E(A, C) + B ⇌ E(A, B, C) :: OnlyA
              E(B, C) + A ⇌ E(A, B, C) :: OnlyA
              E(A, B, C) <--> E(P)     :: OnlyA
              E + P ⇌ E(P)             :: EqualAI
          end
      end
  ```
  In ":OnlyA guard admits feasible ter-substrate mechanisms", the macro now runs admission,
  not the Stiemke test, so the testset reads that test directly: after
  `@test cube isa ER.AllostericEnzymeMechanism` add
  ```julia
      # The moves' :OnlyA test decides it in its Stiemke stage (feasible, w = 1).
      cube_m = ER.AllostericMechanism(cube)
      @test isempty(ER._onlya_offenders(cube_m.catalytic, cube_m.catalytic.group,
                                        ER.cat_allo_states(cube_m)))
  ```
  and in the comment above the testset, "The rejection testset above pins the guard's reject
  direction; this one pins its ACCEPT direction, which nothing else covers. The guard runs"
  becomes "The rejection testset above pins the reject direction of the moves' :OnlyA test
  (`_onlya_offenders`); this one pins its ACCEPT direction, which nothing else covers. The
  test runs". (Probe: the feasible cube and the ordered ter mechanism are admitted on every
  face, the rejected cube `G1`.)

  (e) `test/test_types.jl`: replace the testset ":OnlyA guard admits a ping-pong whose fused
  releases are :OnlyA" with:
  ```julia
  @testset "a ping-pong with both fused releases :OnlyA is an island; one half passes" begin
      # E(A) → E(; residual = A - P) + P runs the first half-reaction and releases P in
      # one step, E(B; residual = A - P) → E + Q the second and releases Q: chemistry
      # (`_is_chemistry`), though each gives off one metabolite and takes up none. Tagged
      # :OnlyA, both leave the moves' :OnlyA test, which reads no cycle through chemistry.
      # Admission reads each edge: in the inactive conformation the covalent forms join
      # free E only through the tagged chemistry, an island (`_ISLAND`) whose weight depends
      # on how fast the tagged constants vanish. Either half alone :OnlyA is admitted, and
      # both oracles agree on every face.
      both = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              E + A <--> E(A)                                       :: EqualAI
              E(A) <--> E(; residual = A - P) + P                   :: OnlyA
              E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
              E(B; residual = A - P) <--> E + Q                     :: OnlyA
          end
      end
      @test isempty(ER._onlya_offenders(both.catalytic, both.catalytic.group,
                                        ER.cat_allo_states(both)))
      v = ER._guard(both)
      @test v.fault === ER._ISLAND
      @test _testhelper_epath_verdict(both; rng = Random.MersenneTwister(6)) in
            (:direction_dependent, :wrong_limit)
      @test_throws ErrorException(ER._guard_message(both, v)) @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              E + A <--> E(A)                                       :: EqualAI
              E(A) <--> E(; residual = A - P) + P                   :: OnlyA
              E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
              E(B; residual = A - P) <--> E + Q                     :: OnlyA
          end
      end
      half1 = @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              E + A <--> E(A)                                       :: EqualAI
              E(A) <--> E(; residual = A - P) + P                   :: OnlyA
              E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
              E(B; residual = A - P) <--> E + Q                     :: EqualAI
          end
      end
      half2 = @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              E + A <--> E(A)                                       :: EqualAI
              E(A) <--> E(; residual = A - P) + P                   :: EqualAI
              E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
              E(B; residual = A - P) <--> E + Q                     :: OnlyA
          end
      end
      # Each face: the RE-speed oracle (class 3) and the ε-path (admission) agree with the
      # guard, which accepts it.
      for (i, half) in enumerate((half1, half2))
          @test all(_testhelper_admission_agreement(ER.AllostericMechanism(half);
                                                    rng = Random.MersenneTwister(6 + i)))
      end
  end
  ```

  (f) `test/test_types.jl`, testset ":OnlyA guard reads cycles past an abortive segment it
  strands": after `@test am isa ER.AllostericEnzymeMechanism` add
  ```julia
      # Admission reads each edge: the abortive bindings reach {E(A, P), E(A, Q),
      # E(A, P, Q)} from the product forms on their unbound side only, so the inactive
      # conformation holds no enzyme there and every face passes (the macro accepted it).
      # The moves' :OnlyA test, which reads no cycle through the chemistry group, agrees.
      am_m = ER.AllostericMechanism(am)
      @test isempty(ER._onlya_offenders(am_m.catalytic, am_m.catalytic.group,
                                        ER.cat_allo_states(am_m)))
  ```


  (g) `test/test_rate_eq_derivation.jl`: replace the testset "All-RE catalytic cycle raises"
  with:
  ```julia
  @testset "All-RE catalytic cycle raises" begin
      # A cycle of rapid-equilibrium steps that turns the reaction over gives the mechanism
      # no finite rate. The guard's class 2 finds it, and the macro throws when its expanded
      # code runs (R27), naming the cycle's steps.
      msg = "close a cycle that turns the reaction over, so the mechanism has no finite rate"
      # Binding stage mixed (S2 SS, S1 RE), release stage mixed (P1 SS, P2 RE), chemistry
      # RE: a complete all-RE catalytic cycle.
      @test_throws msg @enzyme_mechanism begin
          substrates: S1, S2
          products: P1, P2
          steps: begin
              E + S1 ⇌ E(S1)
              E + S2 ⇌ E(S2)
              E(S1) + S2 <--> E(S1, S2)
              E(S2) + S1 ⇌ E(S1, S2)
              E(S1, S2) ⇌ E(P1, P2)
              E(P1, P2) ⇌ E(P1) + P2
              E(P1, P2) <--> E(P2) + P1
              E(P1) ⇌ E + P1
              E(P2) ⇌ E + P2
          end
      end
      # The activator route E(R) + S ⇌ E(S, R) ⇌ E(P, R) ⇌ E(R) + P is an all-RE catalytic
      # cycle.
      @test_throws msg @enzyme_mechanism begin
          substrates: S
          products: P
          modifiers: R
          steps: begin
              E + S ⇌ E(S)
              E(S) ⇌ E(P)
              E + P <--> E(P)
              E(R) + S ⇌ E(S, R)
              E(S, R) ⇌ E(P, R)
              E(R) + P ⇌ E(P, R)
              (E + R ⇌ E(R),
               E(S) + R ⇌ E(S, R),
               E(P) + R ⇌ E(P, R))
          end
      end
  end
  ```

  (h) `test/allosteric_ground_truth.jl` (P4-29a, P4-29b; Task 0 (f) confirmed each variant
  with both oracles). Each half-2 variant below holds its inactive enzyme in the covalent
  forms when both products are zero: half 2 cannot run there, and half 1 runs back only
  with the first product. Its rate is then exactly 0, the mass-action limit too (both
  oracles agree on that face; `$S/phase45/plan/probes/p3c_sink.log` measured 0.0 for all
  three variants on the S3f snapshot, and `$S/phase45/refresh/map-t8-10_1.log` again at
  08c4043, each with a positive peak turnover).

  Replace the comment block above the testset "ping-pong MWC derivation with free E and F
  in one RE segment" (from `# The two ping-pong value gates above` through `(0.28977 at
  these parameters).`) with
  ```julia
  # The two ping-pong value gates above write every step steady state. This gate checks,
  # against the exact two-conformation steady state, a mechanism whose active
  # rapid-equilibrium segment holds both free E and F. Clearing that segment's polynomial
  # gives free E the weight B (F sits at Q/B relative to E), so the inactive term must
  # carry the same factor beside L.
  ```
  and the testset with:
  ```julia
  @testset "ping-pong MWC derivation with free E and F in one RE segment" begin
      # Both half-reactions :OnlyA strand the covalent forms between the halves in the
      # inactive conformation: they join free E only through tagged chemistry, an island
      # (`_ISLAND`) whose weight depends on how fast the tagged constants vanish, so the
      # macro rejects it and the ε-path oracle finds the limit direction-dependent. With the
      # second half alone :OnlyA, the inactive conformation runs the first half both ways,
      # its active rapid-equilibrium segment still holds free E and F = E(; residual = A - P),
      # and the law equals the exact two-conformation steady state. (Half 1 alone fails
      # class 3 at zero substrates and at B = 0.)
      both = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              E + A ⇌ E(A)                                          :: EqualAI
              E(A) <--> E(P; residual = A - P)                      :: OnlyA
              E(; residual = A - P) + P ⇌ E(P; residual = A - P)    :: EqualAI
              E(; residual = A - P) + B ⇌ E(B; residual = A - P)    :: EqualAI
              E(B; residual = A - P) ⇌ E(Q)                         :: OnlyA
              E + Q ⇌ E(Q)                                          :: EqualAI
          end
      end
      v = ER._guard(both)
      @test v.fault === ER._ISLAND
      @test _testhelper_epath_verdict(both; rng = MersenneTwister(663)) in
            (:direction_dependent, :wrong_limit)
      @test_throws ErrorException(ER._guard_message(both, v)) @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              E + A ⇌ E(A)                                          :: EqualAI
              E(A) <--> E(P; residual = A - P)                      :: OnlyA
              E(; residual = A - P) + P ⇌ E(P; residual = A - P)    :: EqualAI
              E(; residual = A - P) + B ⇌ E(B; residual = A - P)    :: EqualAI
              E(B; residual = A - P) ⇌ E(Q)                         :: OnlyA
              E + Q ⇌ E(Q)                                          :: EqualAI
          end
      end
      half2 = @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              E + A ⇌ E(A)                                          :: EqualAI
              E(A) <--> E(P; residual = A - P)                      :: EqualAI
              E(; residual = A - P) + P ⇌ E(P; residual = A - P)    :: EqualAI
              E(; residual = A - P) + B ⇌ E(B; residual = A - P)    :: EqualAI
              E(B; residual = A - P) ⇌ E(Q)                         :: OnlyA
              E + Q ⇌ E(Q)                                          :: EqualAI
          end
      end
      _testhelper_check_law_against_oracle(half2)
  end
  ```
  The harness `_testhelper_pingpong_re_chemistry_flux` and its self-validation testset stay;
  the commit body says that no gate reads it now.

  In the testset "ping-pong :OnlyA I-state keeps a reachable free-enzyme root", replace the
  comment line "The mechanism below is accepted by the `:OnlyA` thermodynamic guard — it is
  valid, and the derivation must handle it." with "With the F6P binding and the first
  half-reaction :OnlyA, the covalent forms E(; residual) and E(F16BP; residual) join free E
  in the inactive conformation only through tagged chemistry: an island (`_ISLAND`), which
  the macro rejects. With the F6P binding and the second half :OnlyA, the inactive
  conformation keeps its free-enzyme root and every face passes; the derivation must handle
  it." and replace the body from `err1 = @allosteric_mechanism begin` through
  `am = ER.AllostericMechanism(err1)` with:
  ```julia
      err1 = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: ATP, F6P
          products: ADP, F16BP
          catalytic_multiplicity: 1
          steps: begin
              E + ATP ⇌ E(ATP)                                                      :: EqualAI
              E(ATP) <--> E(F16BP; residual = ATP - F16BP)                          :: OnlyA
              E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP):: EqualAI
              E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP)    :: OnlyA
              E(F6P; residual = ATP - F16BP) ⇌ E(ADP)                               :: EqualAI
              E + ADP ⇌ E(ADP)                                                      :: EqualAI
          end
      end
      v = ER._guard(err1)
      @test v.fault === ER._ISLAND
      @test _testhelper_epath_verdict(err1; rng = MersenneTwister(761)) in
            (:direction_dependent, :wrong_limit)
      @test_throws ErrorException(ER._guard_message(err1, v)) @allosteric_mechanism begin
          substrates: ATP, F6P
          products: ADP, F16BP
          catalytic_multiplicity: 1
          steps: begin
              E + ATP ⇌ E(ATP)                                                      :: EqualAI
              E(ATP) <--> E(F16BP; residual = ATP - F16BP)                          :: OnlyA
              E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP):: EqualAI
              E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP)    :: OnlyA
              E(F6P; residual = ATP - F16BP) ⇌ E(ADP)                               :: EqualAI
              E + ADP ⇌ E(ADP)                                                      :: EqualAI
          end
      end
      half2 = @allosteric_mechanism begin
          substrates: ATP, F6P
          products: ADP, F16BP
          catalytic_multiplicity: 1
          steps: begin
              E + ATP ⇌ E(ATP)                                                      :: EqualAI
              E(ATP) <--> E(F16BP; residual = ATP - F16BP)                          :: EqualAI
              E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP):: EqualAI
              E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP)    :: OnlyA
              E(F6P; residual = ATP - F16BP) ⇌ E(ADP)                               :: OnlyA
              E + ADP ⇌ E(ADP)                                                      :: EqualAI
          end
      end
      _testhelper_check_law_against_oracle(half2)
      am = ER.AllostericMechanism(half2)
  ```
  and in the rest of that testset read `half2` where it read `err1` (the `fitted_params`,
  `rate_equation` and `_peak_turnover` lines), except the two lines
  ```julia
      # No covalent sink: with both products at zero the enzyme still turns over.
      @test real(ER.rate_equation(err1, merge(concs, (ADP = 0.0, F16BP = 0.0)), prm)) > 0
  ```
  which become
  ```julia
      # With both products at zero the inactive conformation holds its enzyme in the
      # covalent forms: half 2 cannot run there, and half 1 runs back only with F16BP. The
      # rate is then 0, as the mass-action limit is.
      @test real(ER.rate_equation(half2, merge(concs, (ADP = 0.0, F16BP = 0.0)), prm)) == 0
  ```

  Replace the comment block above the testset "dead-inactive ping-pong :OnlyA derivation"
  (from `# ── The gate: dead-inactive ping-pong :OnlyA` through `# a positive peak turnover.`)
  with
  ```julia
  # ── The gate: a dead-inactive ping-pong is rejected; its half-2 variant derives ──────
  # With the F6P binding and both chemical steps :OnlyA the inactive conformation runs
  # neither half-reaction, which waits for iso mechanisms (ROADMAP). The half-2 variant
  # derives: a finite rate as a product goes to zero, a rate of 0 with both products at
  # zero (its inactive conformation holds the enzyme in the covalent forms there), a
  # positive peak turnover, and the non-allosteric twin's rate at L = 0.
  ```
  and the testset with:
  ```julia
  @testset "a dead-inactive ping-pong is rejected; its half-2 variant is its L = 0 twin" begin
      # With the F6P binding and both half-reactions :OnlyA, the inactive conformation runs
      # neither half: the covalent forms between the halves join free E only through tagged
      # chemistry (`_ISLAND`), and the message points at the roadmap, which lists the flip at
      # the covalent form that would hold them (ROADMAP "Let an allosteric ping-pong stop
      # both half-reactions in the inactive conformation"). With the second half alone
      # :OnlyA (half 1 alone fails class 3, as it does for the ping-pong above), the
      # mechanism is admitted.
      dead = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: ATP, F6P
          products: ADP, F16BP
          catalytic_multiplicity: 1
          steps: begin
              E + ATP ⇌ E(ATP)                                                      :: EqualAI
              E(ATP) <--> E(F16BP; residual = ATP - F16BP)                          :: OnlyA
              E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP):: EqualAI
              E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP)    :: OnlyA
              E(F6P; residual = ATP - F16BP) ⇌ E(ADP)                               :: OnlyA
              E + ADP ⇌ E(ADP)                                                      :: EqualAI
          end
      end
      v = ER._guard(dead)
      @test v.fault === ER._ISLAND
      msg = ER._guard_message(dead, v)
      @test occursin("runs neither half-reaction of a ping-pong is not supported yet " *
                     "(see the Roadmap page)", msg)
      @test _testhelper_epath_verdict(dead; rng = MersenneTwister(997)) in
            (:direction_dependent, :wrong_limit)
      @test_throws ErrorException(msg) @allosteric_mechanism begin
          substrates: ATP, F6P
          products: ADP, F16BP
          catalytic_multiplicity: 1
          steps: begin
              E + ATP ⇌ E(ATP)                                                      :: EqualAI
              E(ATP) <--> E(F16BP; residual = ATP - F16BP)                          :: OnlyA
              E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP):: EqualAI
              E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP)    :: OnlyA
              E(F6P; residual = ATP - F16BP) ⇌ E(ADP)                               :: OnlyA
              E + ADP ⇌ E(ADP)                                                      :: EqualAI
          end
      end
      half2 = @allosteric_mechanism begin
          substrates: ATP, F6P
          products: ADP, F16BP
          catalytic_multiplicity: 1
          steps: begin
              E + ATP ⇌ E(ATP)                                                      :: EqualAI
              E(ATP) <--> E(F16BP; residual = ATP - F16BP)                          :: EqualAI
              E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP):: EqualAI
              E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP)    :: EqualAI
              E(F6P; residual = ATP - F16BP) ⇌ E(ADP)                               :: OnlyA
              E + ADP ⇌ E(ADP)                                                      :: EqualAI
          end
      end
      afp = ER.fitted_params(half2)
      # Non-allosteric twin: the same six steps with no tags. At L = 0 the inactive
      # conformation is unpopulated, so the allosteric rate must reduce to the twin's, an
      # independent re-derivation through the plain King–Altman path. An :OnlyA constant
      # carries the A_ tag its twin lacks.
      nonallo = @enzyme_mechanism begin
          substrates: ATP, F6P
          products: ADP, F16BP
          steps: begin
              E + ATP ⇌ E(ATP)
              E(ATP) <--> E(F16BP; residual = ATP - F16BP)
              E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP)
              E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP)
              E(F6P; residual = ATP - F16BP) ⇌ E(ADP)
              E + ADP ⇌ E(ADP)
          end
      end
      nfp = ER.fitted_params(nonallo)
      twin_of(s) = Symbol(replace(String(s), r"^([Kk])_A_" => s"\1_"))
      @test Set(twin_of.(filter(!=(:L), collect(afp)))) == Set(nfp)
      rng = MersenneTwister(20260717)
      for _ in 1:6
          vals = Dict(s => 0.5 + 2rand(rng) for s in nfp)
          Keq = 2.0 + 2rand(rng)
          ATP = 0.5+2rand(rng); F6P = 0.5+2rand(rng)
          ADP = 0.5+2rand(rng); F16BP = 0.5+2rand(rng)
          concs = (ATP=ATP, F6P=F6P, ADP=ADP, F16BP=F16BP)
          prm = NamedTuple{(afp..., :Keq, :E_total)}(
              ((s === :L ? 0.7 : vals[twin_of(s)] for s in afp)..., Keq, 1.0))
          @test isfinite(real(ER.rate_equation(half2, concs, prm)))
          @test isfinite(real(ER.rate_equation(half2,
              (ATP=ATP, F6P=F6P, ADP=ADP, F16BP=1e-8), prm)))
          # With both products at zero the inactive conformation holds its enzyme in the
          # covalent forms (half 2 cannot run, half 1 runs back only with F16BP): rate 0.
          @test real(ER.rate_equation(half2, merge(concs, (ADP=0.0, F16BP=0.0)), prm)) == 0
          @test ER._peak_turnover(half2, prm, (1e-12, 1e3)).value > 0
          # v = 0 at the equilibrium metabolite ratio ADP·F16BP = Keq·ATP·F6P.
          @test isapprox(real(ER.rate_equation(half2,
              (ATP=ATP, F6P=F6P, ADP=Keq*ATP*F6P/F16BP, F16BP=F16BP), prm)),
              0.0; atol=1e-8)
          # L = 0 cross-check: the allosteric rate reduces to the non-allosteric twin.
          pn = NamedTuple{(nfp..., :Keq, :E_total)}(((vals[s] for s in nfp)..., Keq, 1.0))
          p0 = NamedTuple{(afp..., :Keq, :E_total)}(
              ((s === :L ? 0.0 : vals[twin_of(s)] for s in afp)..., Keq, 1.0))
          @test isapprox(real(ER.rate_equation(half2, concs, p0)),
                         real(ER.rate_equation(nonallo, concs, pn)); rtol=1e-4)
      end
  end
  ```
  The dead mechanism's `d_free_I == 1` assertion and its `afp` pin go with its law
  (P4-29a/b, listed in the commit body), and so does its "no covalent sink" check: the
  half-2 variant's rate is 0 with both products at zero, and the commit body says so.

  (i) `test/test_oracles.jl`, testset "a same-model regrouping keeps every fitted name":
  replace the comment's last sentence ("An allosteric split is built only when the boundary
  check's :OnlyA test accepts it.") with "Only a split the exit keeps counts: one the guard
  accepts (`_is_well_posed`).", delete the `if m isa EnzymeRates.AllostericMechanism ... end`
  block that calls `_onlya_offenders`, and after
  `child = EnzymeRates._apply_bipartitions(m, [(g, bp)])` add
  `EnzymeRates._is_well_posed(child) || continue`.


  (j) `test/test_mechanism_enumeration.jl` (P4-29a: the split pre-check's
  fixtures become exit tests).

  Replace the testset "AllostericMechanism — a split the boundary check would reject gains
  nothing" with:
  ```julia
  @testset "AllostericMechanism — an :OnlyA binding with no limit leaves no child" begin
      # Merged random-order bi-bi with the abortive complex E(A, Q). The steady-state A group
      # binds A at E and, fused, at E(B) into E(P, Q), all :OnlyA. In the inactive
      # conformation the :EqualAI steps E + Q <--> E(Q), E(Q) + A <--> E(A, Q) and
      # E(A) + Q <--> E(A, Q) still join E to E(A), the two sides of the :OnlyA binding
      # E + A <--> E(A): thermodynamics keeps its affinity finite there, so it has no limit
      # (`_NO_LIMIT`), and the ε-path oracle finds no path. Every move keeps that binding
      # :OnlyA and those steps untagged: a split's parts inherit their group's tag, and the
      # change-tag move relaxes binding groups only, while the A group holds chemistry. So
      # every child has no limit either, and the exit drops them all.
      am = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              (E + A <--> E(A), E(B) + A <--> E(P, Q))                         :: OnlyA
              (E + B ⇌ E(B), E(A) + B ⇌ E(P, Q))                               :: OnlyA
              (E + P <--> E(P), E(Q) + P <--> E(P, Q))                         :: OnlyA
              (E + Q <--> E(Q), E(A) + Q <--> E(A, Q), E(P) + Q <--> E(P, Q))  :: EqualAI
              E(Q) + A <--> E(A, Q)                                            :: EqualAI
          end
      end
      v = EnzymeRates._guard(am)
      @test v.fault === EnzymeRates._NO_LIMIT
      @test EnzymeRates._step_text(am.catalytic, v.edge) == "E + A <--> E(A)"
      @test _testhelper_epath_verdict(am; rng = Random.MersenneTwister(2922)) === :no_path
      msg = EnzymeRates._guard_message(am, v)
      @test_throws ErrorException(msg) @allosteric_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              (E + A <--> E(A), E(B) + A <--> E(P, Q))                         :: OnlyA
              (E + B ⇌ E(B), E(A) + B ⇌ E(P, Q))                               :: OnlyA
              (E + P <--> E(P), E(Q) + P <--> E(P, Q))                         :: OnlyA
              (E + Q <--> E(Q), E(A) + Q <--> E(A, Q), E(P) + Q <--> E(P, Q))  :: EqualAI
              E(Q) + A <--> E(A, Q)                                            :: EqualAI
          end
      end
      moved = [EnzymeRates._expand_re_to_ss(am); EnzymeRates._expand_split_kinetic_group(am);
               EnzymeRates._expand_change_allo_state(am)]
      @test !isempty(moved)
      @test all(c -> EnzymeRates._guard(c).fault === EnzymeRates._NO_LIMIT, moved)
      kids = EnzymeRates.expand_mechanisms([am], bi_bi_rxn)
      @test length(kids) == 0
      @test Set(kids) == Set()
  end
  ```
  Replace the testset "AllostericMechanism — a rejected split is extended by every :OnlyA
  group" with:
  ```julia
  @testset "AllostericMechanism — beside an inhibitor, :OnlyA A has no limit; no child" begin
      # Merged random-order bi-bi with the competitive inhibitor I and the abortive
      # complexes E(A, P) and E(B, Q); the A, B, I and Q groups :OnlyA. In the inactive
      # conformation the :EqualAI steps E + P <--> E(P), E(P) + A <--> E(A, P) and
      # E(A) + P <--> E(A, P) join E to E(A), the two sides of the :OnlyA binding
      # E + A <--> E(A), so it has no limit (`_NO_LIMIT`). Every child keeps it so: the exit
      # drops them all.
      am = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B
          products: P, Q
          modifiers: I
          steps: begin
              (E + A <--> E(A), E(B) + A <--> E(P, Q), E(I) + A <--> E(A, I))   :: OnlyA
              (E + B ⇌ E(B), E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))            :: OnlyA
              (E + I ⇌ E(I), E(A) + I ⇌ E(A, I))                                :: OnlyA
              (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q))   :: EqualAI
              (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(P) + Q <--> E(P, Q))   :: OnlyA
              E(P) + A <--> E(A, P)                                             :: EqualAI
          end
      end
      v = EnzymeRates._guard(am)
      @test v.fault === EnzymeRates._NO_LIMIT
      @test EnzymeRates._step_text(am.catalytic, v.edge) == "E + A <--> E(A)"
      msg = EnzymeRates._guard_message(am, v)
      @test_throws ErrorException(msg) @allosteric_mechanism begin
          substrates: A, B
          products: P, Q
          modifiers: I
          steps: begin
              (E + A <--> E(A), E(B) + A <--> E(P, Q), E(I) + A <--> E(A, I))   :: OnlyA
              (E + B ⇌ E(B), E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))            :: OnlyA
              (E + I ⇌ E(I), E(A) + I ⇌ E(A, I))                                :: OnlyA
              (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q))   :: EqualAI
              (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(P) + Q <--> E(P, Q))   :: OnlyA
              E(P) + A <--> E(A, P)                                             :: EqualAI
          end
      end
      rxn = @enzyme_reaction begin
          substrates: A[C], B[N]
          products: P[C], Q[N]
          competitive_inhibitors: I
      end
      moved = [EnzymeRates._expand_re_to_ss(am); EnzymeRates._expand_split_kinetic_group(am);
               EnzymeRates._expand_change_allo_state(am)]
      @test !isempty(moved)
      @test all(c -> EnzymeRates._guard(c).fault === EnzymeRates._NO_LIMIT, moved)
      kids = EnzymeRates.expand_mechanisms([am], rxn)
      @test length(kids) == 0
      @test Set(kids) == Set()
  end
  ```
  Replace the testset "AllostericMechanism — the split's :OnlyA test agrees with
  _check_well_posed" with:
  ```julia
  @testset "AllostericMechanism — a split keeps its parent's admission verdict" begin
      # An aggregate check over the bi-bi seed set (the exception to writing every fixture
      # inline): each split child of a seed lifted to its allosteric variants, one parent in
      # ten. Both parts of a split group inherit its tag, and admission's interior reads
      # each edge's tag, so every split the move makes (zero-flux parts at rapid
      # equilibrium) has its parent's interior verdict (P4-14): the move needs no pre-check,
      # and the exit drops the splits of a rejected parent.
      parents = unique(EnzymeRates.AllostericMechanism[
          a for m in EnzymeRates.init_mechanisms(bi_bi_rxn)
          for c in EnzymeRates._expand_split_kinetic_group(m)
          for a in EnzymeRates._expand_to_allosteric(c, bi_bi_rxn)])[1:10:end]
      accepted = rejected = 0
      for am in parents
          cm, states = am.catalytic, EnzymeRates.cat_allo_states(am)
          verdict = EnzymeRates._admission_interior(am).fault
          flux = EnzymeRates._flux_carrying_edges(cm, EnzymeRates._flags(cm))
          units = [(g, bp) for g in 1:EnzymeRates.ngroups(cm)
                   for bp in EnzymeRates._context_bipartitions(cm, g)]
          reverted = [EnzymeRates._revert_zero_flux_parts(cm, bp, flux) for (_, bp) in units]
          for sel in [[[u] for u in eachindex(units)];
                      [[u, v] for u in eachindex(units) for v in u+1:length(units)
                       if units[u][1] != units[v][1]]]
              re = EnzymeRates._flags(cm)
              foreach(u -> re[reverted[u]] .= true, sel)
              child = EnzymeRates._apply_bipartitions(am, units[sel], re)
              @test EnzymeRates._admission_interior(child).fault === verdict
          end
          verdict === EnzymeRates._NO_FAULT ? (accepted += 1) : (rejected += 1)
      end
      # Aggregate pin over the sampled parents: both verdicts occur, so the agreement
      # means something.
      @test (accepted, rejected) == (716, 4)
      @test all(am -> EnzymeRates._expand_split_kinetic_group(am) isa Vector, parents)
  end
  ```
  The pin is `p3_split.log`'s count (720 parents, 716 admitted, 4 rejected by the
  research rules' interior admission on phase 2's moves). If the red run (Step 4) counts
  another population on the branch's moves, first check that each rejected parent's
  `_admission_interior` fault is one `_testhelper_admission_faults` and the ε-path agree
  with (`_testhelper_admission_agreement`), then pin the counts the run prints and give
  both counts in the commit body.

  In the testset "C on the ping-pong seed and its single flips", replace the
  conformational variant `am` (both chemistry groups :OnlyA, an island since Task 8) with
  the half-2 V-type (P4-29b), and the comment above it with:
  ```julia
      # A conformational variant is read on its active-state graph, the seed's: still
      # degenerate. The variant tags the second half-reaction :OnlyA; tagging both would
      # strand the covalent forms in the inactive conformation (`_ISLAND`).
      am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              E + A ⇌ E(A)                                          :: EqualAI
              E(A) <--> E(P; residual = A - P)                      :: EqualAI
              E(; residual = A - P) + P ⇌ E(P; residual = A - P)    :: EqualAI
              E(; residual = A - P) + B ⇌ E(B; residual = A - P)    :: EqualAI
              E(B; residual = A - P) ⇌ E(Q)                         :: OnlyA
              E + Q ⇌ E(Q)                                          :: EqualAI
          end
      end)
  ```
  and replace the closing comment "The probe reads the derived allosteric law. Free E and
  E(; residual = A - P) lie in one rapid-equilibrium segment of the active state; the law
  divides each conformation by its free enzyme's weight, so at zero products it does not
  depend on B (the gate in test/allosteric_ground_truth.jl checks it against mass action)."
  with "The probe reads the derived allosteric law, which vanishes at zero products: the
  inactive conformation holds its enzyme in the covalent forms there, since half 2 cannot
  run and half 1 runs back only with P, so the probe reads no rate to compare and calls the
  law degenerate too (the gate in test/allosteric_ground_truth.jl checks the half-2 V-type's
  law against mass action)."

  In the testset "Mechanism — flip past a bottomless RE segment", the children do not change
  (the probe: each of the four passes class 3 on every face and agrees with the RE-speed
  oracle; the pair {P@E, Q@E} fails class 3 at zero products, `_BOTTOMLESS`). Replace "would
  leave {E(P), E(Q), E(P, Q)} with no bottom form — a mechanism the boundary check rejects —
  so that pair counts as failed and is not emitted." with "would leave {E(P), E(Q), E(P, Q)}
  with no bottom form, which class 3 rejects at zero products inside the flip move's
  predicate (`_faces_pass`), so that pair counts as failed and is not emitted.", and after
  `@test Set(children) == Set(expected)` add:
  ```julia
      for (i, c) in enumerate(children)
          @test all(==(:agree),
                    _testhelper_respeed_verdicts(c; rng = Random.MersenneTwister(1302 + i)))
      end
      # The pair {P@E, Q@E} fails class 3 at zero products inside the move's predicate.
      pq = [_testhelper_edge(parent, "E + P ⇌ EP"), _testhelper_edge(parent, "E + Q ⇌ EQ")]
      pair = EnzymeRates._flags(parent)
      pair[pq] .= false
      @test !EnzymeRates._faces_pass(parent, pair)
  ```

  In the testset "Mechanism — Bi-bi ping-pong: binding-group :OnlyA variants only", the
  fifteen children are islands (both chemistry groups :OnlyA) that the exit drops; Task 12's
  move replaces them. Their laws are not derived any more, so the parameter count comes from
  the joint solve, which runs no guard: replace
  `@test _testhelper_param_deltas(m, result) == fill(1, 15)` with
  ```julia
      @test [EnzymeRates._independent_param_count(r) -
             EnzymeRates._independent_param_count(m) for r in result] == fill(1, 15)
      # Every child strands the covalent forms in the inactive conformation, so the exit
      # drops it (`_ISLAND`).
      @test all(r -> EnzymeRates._guard(r).fault === EnzymeRates._ISLAND, result)
  ```
  and in the comment above it "All 15 variants are Δ=1." stays. The next comment's
  "(the deltas above compile every result)" is false once the deltas come from the joint
  solve: "3. invariants (the deltas above compile every result), and no result is the
  all-:EqualAI baseline." becomes "3. invariants, and no result is the all-:EqualAI
  baseline."

  In the testsets "change_allo_state drops partial-catalysis relaxations" and
  "change_allo_state drops mixed-iso V-type relaxations (no :OnlyA binding)", the parents
  (both iso steps :OnlyA) are islands, which the macro rejects: build each with
  `@_testhelper_unguarded @allosteric_mechanism begin ... end` (same block) and use the
  result as `am`: write `am = @_testhelper_unguarded @allosteric_mechanism begin` in place
  of `dead = @allosteric_mechanism begin` (`vtype = @allosteric_mechanism begin` in the
  second) and delete the line `am = EnzymeRates.AllostericMechanism(dead)`
  (`AllostericMechanism(vtype)`). In the second, the
  end-to-end rate check derives each child's law: put
  `EnzymeRates._is_well_posed(k) || continue` before `cm = EnzymeRates.compile_mechanism(k)`,
  and after the loop add
  ```julia
      # The binding relaxations keep both iso steps :OnlyA, islands the exit drops; the
      # chemistry relaxation is admitted and turns over.
      @test count(EnzymeRates._is_well_posed, kids) >= 1
      @test all(k -> EnzymeRates._is_well_posed(k) ||
                     EnzymeRates._guard(k).fault === EnzymeRates._ISLAND, kids)
  ```
  Task 12 re-pins these move-level tests on P4-1's sets.

  (k) `test/test_identify_rate_equation.jl`: replace the testset "_base_tier
  records a flip child whose rate law is not defined" with:
  ```julia
  @testset "_base_tier drops a flip child whose rate law the guard rejects" begin
      # Every step is rapid equilibrium, so the catalytic cycle E, E(S), E(P), F turns over
      # at rapid equilibrium: the seed is degenerate and expanded. Its free forms E and F
      # share one rapid-equilibrium segment through F ⇌ E. The six flips of two steps of the
      # ring cut it; the three that cut it at F ⇌ E leave E and F in two segments, which the
      # one-flip check rejects, so the base tier drops them and keeps the other three (each
      # passes every face of the guard). No failure is recorded.
      ER = EnzymeRates
      rxn = @enzyme_reaction begin
          substrates: S[C]
          products: P[C]
      end
      seed = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)      :: EqualAI
              E(S) ⇌ E(P)       :: EqualAI
              E(P) ⇌ F + P      :: EqualAI
              F ⇌ E             :: EqualAI
          end
      end
      # Untyped: inside `T[...]` the parser reads `begin` as an index keyword.
      kept = [
          ER.AllostericMechanism(@allosteric_mechanism begin
              substrates: S
              products: P
              steps: begin
                  E + S <--> E(S)   :: EqualAI
                  E(S) <--> E(P)    :: EqualAI
                  E(P) ⇌ F + P      :: EqualAI
                  F ⇌ E             :: EqualAI
              end
          end),
          ER.AllostericMechanism(@allosteric_mechanism begin
              substrates: S
              products: P
              steps: begin
                  E + S <--> E(S)   :: EqualAI
                  E(S) ⇌ E(P)       :: EqualAI
                  E(P) <--> F + P   :: EqualAI
                  F ⇌ E             :: EqualAI
              end
          end),
          ER.AllostericMechanism(@allosteric_mechanism begin
              substrates: S
              products: P
              steps: begin
                  E + S ⇌ E(S)      :: EqualAI
                  E(S) <--> E(P)    :: EqualAI
                  E(P) <--> F + P   :: EqualAI
                  F ⇌ E             :: EqualAI
              end
          end)]
      @test ER._degenerate(seed)
      @test length(ER._expand_re_to_ss(seed)) == 6
      @test count(!ER._is_well_posed, ER._expand_re_to_ss(seed)) == 3
      base, failures = ER._base_tier([seed], rxn)
      @test isempty(failures)
      @test length(base) == 3
      @test Set(base) == Set(kept)
  end
  ```
  (The seed itself is rejected, class 2, so it is built unguarded; `_base_tier` reads the
  seeds it is given. The probe finds all six flips well-posed under the guard's classes;
  Task 9 deletes the one-flip check and re-pins this testset to all six.)


  (l) Task 3's class-3 testsets read the package law of fixtures whose derived denominator
  loses every monomial on a face; from this task the derivation refuses those laws
  (`_assert_support`), so the oracle reads its own deletion limit there. In
  `test/test_admission.jl`:
  - "class 3 rejects a bottomless segment": `_testhelper_respeed_verdicts(entered; rng =
    Random.MersenneTwister(6))` becomes `_testhelper_respeed_verdicts(entered; rng =
    Random.MersenneTwister(6), law = :limit)`, and the last assertion
    (`_testhelper_respeed_verdicts(unfed; rng = Random.MersenneTwister(9))[2] === :law_nan`)
    becomes
    ```julia
        # The package law is 0/0 at zero products: class 5 refuses to derive it.
        @test EnzymeRates._support_verdict(unfed, EnzymeRates._law_branches(unfed)).face == 2
        @test_throws "every term of the rate law's denominator vanishes" (
            EnzymeRates._num_den_exprs(unfed))
    ```
  - "the limit test decides what structure cannot": the `dead_end` row's oracle verdict
    `:law_nan` comes from the package law, which class 5 now refuses. The oracle's own
    deletion limit gives the same verdict there: approaching zero products, E(P) gathers
    enzyme in proportion to P/Q, so the face value depends on the direction and the limit
    law is NaN. Keep the row, and in the loop replace
    `oracles = _testhelper_respeed_verdicts(m; rng = Random.MersenneTwister(20 + i))` with
    ```julia
            # dead_end's package law is 0/0 at zero products, which class 5 refuses to
            # derive; the oracle's own deletion limit has no value there either.
            oracles = _testhelper_respeed_verdicts(m; rng = Random.MersenneTwister(20 + i),
                                                   law = m === dead_end ? :limit : :package)
    ```
    and add after the loop:
    ```julia
        @test EnzymeRates._support_verdict(dead_end,
                                           EnzymeRates._law_branches(dead_end)).face == 2
    ```
  - "class 3 rejects a second closed class on a face": `detour`'s law loses every
    denominator monomial at zero products and modifiers (face 5). In the loop over `cases`,
    replace `verdicts = _testhelper_respeed_verdicts(m; rng = Random.MersenneTwister(10 +
    i))` with
    ```julia
            # detour's package law is 0/0 at face 5, which class 5 refuses to derive.
            verdicts = _testhelper_respeed_verdicts(m; rng = Random.MersenneTwister(10 + i),
                                                    law = m === detour ? :limit : :package)
    ```
    The face chain's second closed class decides face 5 before any law is read, so the
    verdicts stay `:nonunique` there and `:agree` elsewhere.
  - "forms reached only through tagged chemistry are an island" (Task 4): `both`
    runs neither half-reaction in the inactive conformation, so its message gains the
    roadmap sentence of Step 5 (a). In its pinned message replace the last line
    `":EqualAI or :NonequalAI."` with
    ```julia
              ":EqualAI or :NonequalAI. An inactive conformation that runs neither " *
              "half-reaction of a ping-pong is not supported yet (see the Roadmap page)."
    ```
    The other `_ISLAND` messages Task 4 pins (R3, R4, AD2b, the B-binding ping-pong) keep a
    chemistry step live in the inactive conformation and do not change.
  - `p3_class5.log` lists no other fixture of Tasks 1-7 whose package law a test
    reads: its other rows are the second-free-form tag sets of `test/test_admission.jl`
    (R2, R3, BS1, BS4, S1) and R2 in "the ε-path oracle rejects islands, and R2 only at
    P = 0", which the tests judge by the oracles' own deletion law
    (`law = :limit`, the ε-path) only. If Step 7 finds another, its oracle call takes
    `law = :limit`, or reads class 5's verdict where it asserted `:law_nan`.

  (m) The other fixtures `p3_scan.log` lists, each a testset whose subject is not the law
  the guard rejects. The probe gave every verdict below, and `p3_tied.log`,
  `p3_ssmirror.log` and `p3_class5_tests.log` the variants' orientation, names, counts and
  class-5 verdicts (class 5 holds for every fixture here).
  - `test/test_dsl.jl`, testset "DSL: call-form notation parses each form's label and bound
    metabolites": its first fixture runs all three
    steps at rapid equilibrium, a fast turnover cycle. Replace
    `m = @enzyme_mechanism begin` with `mech = @_testhelper_unguarded @enzyme_mechanism begin`
    (the block unchanged), and the line `mech = EnzymeRates.Mechanism(m)` with
    ```julia
            # Every step at rapid equilibrium closes a fast cycle that turns the reaction
            # over (class 2); the parse is the subject, so the constructor builds it
            # unchecked.
            @test EnzymeRates._guard(mech).fault === EnzymeRates._RE_TURNOVER
    ```
  - `test/test_types.jl`, testsets "a kinetic group ties one metabolite exchange at two
    conformations" and "a tied exchange turns to take up its lead's metabolites": at zero
    modifiers E(J) and F(J) are dead ends, so the face has two closed classes (`_NONUNIQUE`,
    face 4). In every fixture of both testsets except the `@test_throws "make different
    conformational changes"` one, add after the line `(E + I ⇌ E(I), F + I ⇌ F(I))` the
    line `(E + J ⇌ E(J), F + J ⇌ F(J))` (in the `@allosteric_mechanism` fixtures
    `(E + J ⇌ E(J), F + J ⇌ F(J))                      :: EqualAI`); in the two `mixed`
    fixtures of the second testset add `E + J ⇌ E(J)` instead (J binds only E there). The
    probe checked every variant: admitted on every face (`p3_verdicts_x.log`), each
    exchange group oriented as the testsets assert and each constant named as they assert
    (`p3_tied.log`). Append to each testset's first comment: "J binds E and F as well:
    without those bindings E(J) and F(J) are dead ends at zero modifiers, a second closed
    class the guard rejects." and at the end of the first testset add
    ```julia
        # Without the J bindings the guard rejects the mechanism at zero modifiers.
        dead_ends = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: A
            products: P
            modifiers: I, J
            steps: begin
                E + A <--> E(A)
                E(A) <--> E(P)
                E(P) <--> F + P
                F <--> E
                (E + I ⇌ E(I), F + I ⇌ F(I))
                (E(I) + J <--> E(J) + I, F(I) + J <--> F(J) + I)
            end
        end
        v = EnzymeRates._guard(dead_ends)
        @test (v.fault, v.face) == (EnzymeRates._NONUNIQUE, 4)
        @test _testhelper_respeed_verdicts(dead_ends; rng = Random.MersenneTwister(81),
                                           law = :limit)[4] === :nonunique
    ```
  - `test/test_mechanism_enumeration.jl`, testset "_flux_carrying_groups": the
    fixture `loop` fails class 3 at zero products. Replace
    `loop = EnzymeRates.Mechanism(@enzyme_mechanism begin` and its closing `end)` with
    `loop = @_testhelper_unguarded @enzyme_mechanism begin` and `end`, and add after it
    ```julia
        # At zero products its segment's steady state depends on how fast the
        # rapid-equilibrium steps run (class 3); the flux test reads the graph alone.
        v = EnzymeRates._guard(loop)
        @test v.face == 2 && v.limit === EnzymeRates._LIMIT_SPEED
    ```
  - Same file, testset "MERGE, ELIM and the seed tests on the ordered bi-bi seed": the two
    expected results, the merged groups and the Theorell–Chance groups, run every step at
    rapid equilibrium (the testset asserts `_re_turnover_cycle` on both). Each is written
    `ER.Mechanism(@enzyme_mechanism begin ... end)`: write
    `@_testhelper_unguarded(@enzyme_mechanism begin` in place of
    `ER.Mechanism(@enzyme_mechanism begin` (the block and its closing `end)` unchanged; the
    macro returns the concrete mechanism, for which `Mechanism` has no method), and add
    above the first:
    `# The expected groups turn the reaction over at rapid equilibrium, which the guard
    rejects; the constructor builds them unchecked.` Testset "ELIM refuses the uni-uni
    merged complex": the same change to its expected `E + S ⇌ E(P)`, `E + P ⇌ E(P)`
    groups, with the same comment. Testset "_seed_variants: the ordered/random seed": the
    same change to its expected merged base, every step at rapid equilibrium
    (`EnzymeRates.Mechanism(@enzyme_mechanism begin` there), with the same comment.
  - Same file, testset "_expand_re_to_ss: a split metabolite flips as a pair, never
    singly": the first absent single-A flip (`E + A <--> E(A)` alone) fails class 3 at zero
    substrates. Replace `for single in (` and that tuple member
    `EnzymeRates.Mechanism(@enzyme_mechanism begin ... end),` with
    ```julia
        # The first also fails class 3 at zero substrates: its segment's steady state
        # depends on how fast the rapid-equilibrium steps run.
        single_A = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(B) + A ⇌ E(A, B)
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
                E(A, B) <--> E(P, Q)
            end
        end
        v = EnzymeRates._guard(single_A)
        @test v.face == 3 && v.limit === EnzymeRates._LIMIT_SPEED
        for single in (
            single_A,
    ```
    The parent, the four expected children and the second single flip pass the guard
    (the probe flagged none of them), so `kids` keeps its four.
  - Same file, testset "_redundant_copy_groups: every complex a twin and a consistent dwell
    gauge": `ss_mirror` fails class 3 at zero substrates; its subject, the copy rule on its
    graph, stands. Replace `ss_mirror = ER.Mechanism(@enzyme_mechanism begin` and its
    closing `end)` with `ss_mirror = @_testhelper_unguarded @enzyme_mechanism begin` and
    `end`, and the assertion
    `@test _testhelper_fitted(ss_mirror) == 8 &&` with its continuation line
    `_testhelper_identifiable_rank(ss_mirror) == 6`
    with
    ```julia
        # The guard rejects ss_mirror at zero substrates, where its segment's steady state
        # depends on how fast its rapid-equilibrium steps run; its law is still counted
        # and ranked, unchecked.
        v = EnzymeRates._guard(ss_mirror)
        @test v.face == 3 && v.limit === EnzymeRates._LIMIT_SPEED
        @test ER._independent_param_count(ss_mirror) == 8 &&
              _testhelper_identifiable_rank(ss_mirror) == 6
    ```
    and give `_testhelper_identifiable_rank` (top of the file) the unchecked path: replace
    its docstring and its first three lines (`em = ...`, `fp = ...`, `mets = ...`) with
    ```julia
    """Finite-difference rank of ∂v/∂log θ over the fitted parameters (test oracle only). The
    law of a mechanism the guard rejects is read unchecked from its polynomials
    (`_num_den_exprs`) and evaluated in BigFloat."""
    function _testhelper_identifiable_rank(m; npts = 60, ndraws = 3, h = 1e-5)
        mech = _testhelper_lift(m)
        dep, fitted = EnzymeRates._dependent_param_exprs(mech)
        rate = if EnzymeRates._is_well_posed(mech)
            em = EnzymeRates.compile_mechanism(mech)
            (c, p) -> rate_equation(em, c, p)
        else
            num, den = EnzymeRates._num_den_exprs(mech)
            function (c, p)
                env = Dict{Symbol, BigFloat}(k => BigFloat(v)
                                             for (k, v) in pairs(merge(p, c)))
                for (k, ex) in dep
                    env[k] = _testhelper_eval_expr(ex, env)
                end
                Float64(_testhelper_eval_expr(num, env) / _testhelper_eval_expr(den, env))
            end
        end
        fp = collect(fitted)
        mets = sort(collect(EnzymeRates.metabolites(mech)))
    ```
    and in its loop `J[i, j] += sgn * rate_equation(em, c, p) / (2h)` with
    `J[i, j] += sgn * rate(c, p) / (2h)`. A well-posed mechanism takes the compiled path as
    before, so every other rank pin reads the same law.

- [ ] **Step 4: Run the new and converted tests and watch them fail.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  for f in test_admission test_types test_dsl test_rate_eq_derivation allosteric_ground_truth test_oracles test_identify_rate_equation; do $S/phase3/jl.sh $S/phase45/logs/P4_T8_red_$f.log 1800 --project=. -e "using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include(\"test/mechanism_definitions_for_test_enzyme_derivation.jl\"); include(\"test/$f.jl\")"; done
  ```
  then the same for `test_mechanism_enumeration` in the background (about 150 s), log
  `P4_T8_red_test_mechanism_enumeration.log`. `_is_well_posed` does not exist until Step 5,
  so every test that calls it, directly or through `_testhelper_identifiable_rank`, errors
  with ``UndefVarError: `_is_well_posed` not defined``. Expected failures, and nothing else:
  - test_admission.jl: the two new testsets fail (`No exception thrown` for each
    `@test_throws` on a macro and on the parser); in "class 3 rejects a bottomless segment"
    the class-5 lines fail (`_num_den_exprs` throws nothing yet); the island testset's
    `both` pin fails (its message lacks the roadmap sentence).
  - test_types.jl: every converted `@test_throws ErrorException(...)` fails with a different
    message (the old bottomless and `:OnlyA` texts) or `No exception thrown` (the
    ping-pong); the `_admission` assertions inside the validator pass.
  - test_rate_eq_derivation.jl: "All-RE catalytic cycle raises" fails `No exception thrown`
    (the macro builds; the error came at derivation).
  - allosteric_ground_truth.jl: the three converted testsets fail at their `@test_throws`;
    their half-2 checks pass (phase 3 derives the same laws).
  - test_identify_rate_equation.jl: "_base_tier drops a flip child whose rate law the guard
    rejects" errors at `_is_well_posed`.
  - test_mechanism_enumeration.jl: the two new testsets fail (`expand_mechanisms` returns
    six children where one is expected; `_expand_re_to_ss` returns the three children the
    exit-only check gives, where five are expected); the two split exit testsets fail
    (`kids` holds the children the old exit kept, and the macros throw no guard message);
    every testset that
    reads `_testhelper_identifiable_rank` and "change_allo_state drops mixed-iso V-type
    relaxations" error at `_is_well_posed`; the aggregate split testset, "C on the
    ping-pong seed", "flip past a bottomless RE segment" and the ping-pong tagging testset
    pass.
  - test_oracles.jl: "a same-model regrouping keeps every fitted name" errors at
    `_is_well_posed`; every other testset passes.
  - test_dsl.jl and the other testsets (m) converts: pass. Each conversion holds under
    phase 3's check and under the guard, so these runs show only that the edits are sound;
    the full suite after the implementation (Step 7) shows the guard needs nothing more
    there.


- [ ] **Step 5: Wire the guard into the boundary function and add class 5's throw.**

  (a) `src/admission.jl`: replace the docstring and both methods of `_check_well_posed` with:
  ```julia
  """
      _check_well_posed(m)

  Error unless `m`'s rate law is the unique limit of mass-action kinetics, as its
  rapid-equilibrium steps speed up and its tagged constants vanish, on every face where data
  sit (spec section 8): the guard (`_guard`) finds no fault, and the error is its message
  (`_guard_message`), which names each step at fault by its DSL line. An allosteric
  mechanism must also hold every free form in one rapid-equilibrium segment of the edges
  outside `:OnlyA` groups (`_assert_one_flip_segment`). It reads the reaction graph through
  masks and builds nothing. The mechanism macros, the parser and every derivation entry
  (`_derivation_input`) call it; the enumerator's exits and identify's base tier keep what
  `_is_well_posed` accepts; the constructors and the Sig decoder call neither, so code that
  builds a mechanism internally never trips a check meant for the mechanisms it returns.
  Class 5 reads the derived law and runs where the law is derived (`_assert_support`).
  """
  function _check_well_posed(m::_AnyMech)
      m isa AllostericMechanism && _assert_one_flip_segment(m)
      v = _guard(m)
      v.fault === _NO_FAULT || error(_guard_message(m, v))
      nothing
  end

  """Whether `_check_well_posed(m)` passes: the enumerator's exits and identify's base tier
  keep a child exactly when it does."""
  _is_well_posed(m::Mechanism) = _guard(m).fault === _NO_FAULT
  _is_well_posed(am::AllostericMechanism) =
      _free_forms_flip_together(am) && _guard(am).fault === _NO_FAULT
  ```
  Delete `_assert_bottomed_segments` with its docstring (class 3's `_BOTTOMLESS` replaces it on
  every face; `_bottomless_segment` stays for the split move). Replace the body of
  `_assert_one_flip_segment(am)` (its docstring stays) with:
  ```julia
  function _assert_one_flip_segment(am::AllostericMechanism)
      _free_forms_flip_together(am) && return nothing
      m = am.catalytic
      free = [f for f in eachindex(m.forms) if _is_free(m.forms[f])]
      error("AllostericMechanism: the free forms " *
            "$(join([_form_name(m, f) for f in free], ", ")) " *
            "are not joined by rapid-equilibrium steps outside :OnlyA groups. The " *
            "conformations interconvert at one free form, so the other free forms must be " *
            "in rapid equilibrium with it: join them by ⇌ steps outside :OnlyA groups, or " *
            "keep one free form.")
  end
  ```
  and add, directly before it:
  ```julia
  """Whether every free form (`_is_free`) of `am`'s catalytic graph lies in one
  rapid-equilibrium segment of the edges outside `:OnlyA` groups
  (`_assert_one_flip_segment`)."""
  function _free_forms_flip_together(am::AllostericMechanism)
      m = am.catalytic
      free = [f for f in eachindex(m.forms) if _is_free(m.forms[f])]
      length(free) <= 1 && return true
      segment_of = _segments(m, trues(length(m.edges)),
                             _flags(m) .& _outside_onlya(am)).segment_of
      all(f -> segment_of[f] == segment_of[first(free)], free)
  end
  ```
  After `_support_verdict` (Task 6), add:
  ```julia
  """
      _assert_support(m, branches)

  Error with class 5's message (`_guard_message`) when the derived denominator of `m`, given
  as `branches` (`_law_branches`), loses every monomial on some face (`_support_verdict`):
  the law is 0/0 there. Both `_num_den_exprs` methods call it with the polynomials they hold.
  """
  function _assert_support(m::_AnyMech, branches::Vector{Tuple{POLY, POLY}})
      v = _support_verdict(m, branches)
      v.fault === _NO_FAULT || error(_guard_message(m, v))
      nothing
  end
  ```
  In `_admission_message` (Task 4), replace the `_ISLAND` branch's last statement
  ```julia
          return c == 2 && fixable ?
                 msg * " Tag the group of $(line(k)) :EqualAI or :NonequalAI." : msg
  ```
  with
  ```julia
          msg = c == 2 && fixable ?
                msg * " Tag the group of $(line(k)) :EqualAI or :NonequalAI." : msg
          # A ping-pong whose inactive conformation runs neither half-reaction strands the
          # covalent forms between the halves (P4-2); the roadmap lists what would hold them.
          neither = c == 2 && all(j -> !_is_chemistry(m, m.edges[j]) || _excludes(am, j, 2),
                                  eachindex(m.edges))
          covalent = any(f -> comp[f] == lost && _has_residual(m.forms[f]),
                         eachindex(m.forms))
          return neither && covalent ?
                 msg * " An inactive conformation that runs neither half-reaction of a " *
                       "ping-pong is not supported yet (see the Roadmap page)." : msg
  ```

  (b) `src/rate_eq_derivation.jl`: in `_compute_alpha`, replace the docstring's last three
  sentences ("`segments` are the sub-graph's RE segments (`_segments`), `names[k]` the
  parameter Symbols of edge `k` (the K first), and `conc_set` the concentration symbols.
  Raises when the RE edges close a catalytic cycle, which gives the mechanism no finite
  rate.") with
  ```julia
  `segments` are the sub-graph's RE segments (`_segments`) and `names[k]` the parameter
  Symbols of edge `k` (the K first). The guard's class 2 rejects every mechanism whose RE
  edges close a cycle that turns the reaction over (`_re_turnover_transitions`), so the
  weights reproduce each RE edge's own ratio.
  ```
  its signature with
  ```julia
  function _compute_alpha(m::Mechanism, flat::Vector{Int},
                          segments::Vector{Vector{Int}}, names::Vector{Vector{Symbol}})
  ```
  and delete the block from the comment `# Every RE edge must reproduce its own concentration
  ratio from the weights; a` through the `for k in flat ... end` loop that follows it, so the
  function ends `alpha` after the segment loop. In `_raw_symbolic_rate_polys`, the call becomes
  `alpha = _compute_alpha(cm, flat, segments.segments, names)`.

  Replace `_num_den_exprs(m::Mechanism; horner = false)` and its docstring with:
  ```julia
  """The numerator and denominator Exprs of the rate equation of `m`, from its raw
  symbolic rate polys: flat sums for the printed equation, or in Horner form
  (`_poly_to_horner`) for the generated code when `horner` is true. Errors when the
  denominator vanishes on a face of the guard (`_assert_support`)."""
  function _num_den_exprs(m::Mechanism; horner = false)
      names = _edge_names(m, :None)
      num, den, _ = _raw_symbolic_rate_polys(m, trues(length(m.edges)), names)
      _assert_support(m, [(den, poly_one())])
      horner && return _poly_to_horner(num), _poly_to_horner(den)
      param_syms = Set{Symbol}(n for ns in names for n in ns)
      _poly_to_expr(num, param_syms), _poly_to_expr(den, param_syms)
  end
  ```
  In `_num_den_exprs(am::AllostericMechanism; horner = false)`, directly after the statement
  that assigns `num_A_poly, den_A_poly, num_i_poly, den_i_poly, d_A, d_I`, add
  ```julia
      # Each conformation's denominator times the other's free-enzyme weight
      # (`_law_branches`): every term is positive, so no face can cancel them.
      _assert_support(am, [(den_A_poly, d_I), (den_i_poly, d_A)])
  ```
  and append to its docstring: "Errors when the denominator vanishes on a face of the guard
  (`_assert_support`)." In `_state_rate_polys`' docstring, replace "Every derivation entry
  first runs the boundary check (`_check_well_posed`), which rejects a rapid-equilibrium
  segment with no bottom form in either graph." with "Every derivation entry first runs the
  boundary check (`_check_well_posed`), whose class 3 rejects, among other faults, a face
  where a segment of either conformation's law network has no bottom form."

  (c) `src/mechanism_enumeration.jl`: replace `_minimal_flips` and its docstring
  with:
  ```julia
  """
      _minimal_flips(admissible, m, eligible = _ -> true) -> Vector{BitVector}

  The rapid-equilibrium flags of the catalytic graph of `m`, plain or allosteric, after
  every inclusion-minimal flip to steady state that passes, one flag per edge for each flip
  set. A unit is a rapid-equilibrium kinetic group that carries flux with every edge steady
  state (`_flux_carrying_groups`) and passes `eligible`. A candidate turns its units' edges
  steady state and keeps every other flag. It passes when `admissible` holds on its flags,
  every flipped group carries flux in it, and the guard's classes 2 and 3 pass on every face
  (`_faces_pass`, which reads `m`'s law network under the candidate's flags). The sets come
  from `_minimal_gaining_sets`, so a failing set is extended and a passing one is not.
  Nothing is built: a candidate is its flags.
  """
  function _minimal_flips(admissible, m::_AnyMech, eligible = _ -> true)
      cm = catalytic(m)
      flags = _flags(cm)
      members = _group_members(cm)
      flux = _flux_carrying_groups(cm, falses(length(cm.edges)))
      units = [g for g in eachindex(members)
               if all(k -> flags[k], members[g]) && flux[g] && eligible(g)]
      function flipped(sel)
          re = copy(flags)
          for u in units[sel], k in members[u]
              re[k] = false
          end
          re
      end
      passes(sel) = (re = flipped(sel);
                     admissible(re) && all(_flux_carrying_groups(cm, re)[units[sel]]) &&
                     _faces_pass(m, re))
      [flipped(sel) for sel in _minimal_gaining_sets(passes, _ -> 1:length(units))]
  end
  ```
  In `_expand_re_to_ss`, the call becomes
  `flips = _minimal_flips(re -> _re_segment_count(cm, re) > base, m, eligible)`, and in its
  docstring the sentence "A flip set that leaves a rapid-equilibrium segment with no bottom
  form (`_bottomless_segment`, which the boundary check rejects) also counts as failed, so the
  minimal-set search extends it instead of emitting it;" becomes "A flip set whose child the
  guard's class 3 rejects on some face (`_faces_pass`) also counts as failed, so the
  minimal-set search extends it instead of emitting it;". The seed variants call
  `_minimal_flips` on a plain base and need no change.

  Replace both `_split_gain_test` methods and their docstrings with:
  ```julia
  """The split move's gain test for a plain mechanism: `gains(sel)` is true when applying
  the selected units, with the flags `flags_of(sel)`, raises the independent-parameter count
  above the parent's. Each edge is counted under its group in the child
  (`_bipartitioned_groups`) and the kind its flag in the child gives it."""
  function _split_gain_test(m::Mechanism, units, flags_of)
      counter = _partition_independent_count(m)
      flat = _flat_edges(m)
      kinds(re) = Symbol[_count_kind(m.edges[k], re[k]) for k in flat]
      base = counter(m.group[flat], kinds(_flags(m)))
      function gains(sel)
          group, _ = _bipartitioned_groups(m, units[sel])
          counter(group[flat], kinds(flags_of(sel))) > base
      end
      gains
  end

  """
  The split move's gain test for an allosteric mechanism, which counts each candidate on the
  built child by the joint solve (`_independent_param_count`); it derives no rate
  polynomial. Both parts of a split group inherit its tag and admission's interior reads
  each edge's tag, so a split never changes the interior verdict (`_admission_interior`);
  the faces can change through the new group ties and the parts it reverts to rapid
  equilibrium, and the exit (`expand_mechanisms`) drops such a child.
  """
  function _split_gain_test(am::AllostericMechanism, units, flags_of)
      base = _independent_param_count(am)
      function gains(sel)
          group, origin = _bipartitioned_groups(am.catalytic, units[sel])
          child = AllostericMechanism(_raw(am, flags_of(sel), group,
                                           cat_allo_states(am)[origin]))
          _independent_param_count(child) > base
      end
      gains
  end
  ```
  In `_expand_split_kinetic_group`, `gain, balancers = _split_gain_test(m, units, flags_of)`
  becomes `gain = _split_gain_test(m, units, flags_of)`, and `partners` loses the balancers:
  ```julia
      partners(sel) = begin
          isempty(sel) && return 1:length(units)
          used = Set(units[u][1] for u in sel)
          touched = reduce(union, (segments[units[u][1]] for u in sel))
          [u for u in 1:length(units)
           if !(units[u][1] in used) && !isempty(intersect(segments[units[u][1]], touched))]
      end
  ```
  In its docstring, the passage "and by the joint solve on the built child for an
  `AllostericMechanism`, where a candidate the boundary check would reject (`_onlya_offenders`)
  is not built and gains nothing (`_split_gain_test`). A candidate the boundary check would
  reject also takes as partners the units of every `:OnlyA` group: only a split of an `:OnlyA`
  group changes that check, and the split that balances the candidate's cycle may share no
  rapid-equilibrium segment with it." becomes "and by the joint solve on the built child for
  an `AllostericMechanism` (`_split_gain_test`)."

  In `expand_mechanisms`, replace the loop
  ```julia
      for child in result
          _assert_atom_conserving(child, rxn)
          _check_well_posed(child)
      end
  ```
  with
  ```julia
      foreach(child -> _assert_atom_conserving(child, rxn), result)
      filter!(_is_well_posed, result)
  ```
  and in its docstring "Every child that passes the regulator-type filter must conserve atoms
  (`_assert_atom_conserving`) and give a defined rate law (`_check_well_posed`); either check
  errors on the first child that fails it." becomes "Every child that passes the
  regulator-type filter must conserve atoms (`_assert_atom_conserving`), which errors on the
  first child that does not; a child whose rate law the guard rejects (`_is_well_posed`) is
  dropped."

  In `init_mechanisms`, `foreach(_check_well_posed, out)` becomes
  `filter!(_is_well_posed, out)`, and its docstring's "Errors when the reaction has no
  catalytic topology, and when a mechanism it returns gives no defined rate law
  (`_check_well_posed`)." becomes "Errors when the reaction has no catalytic topology; a
  mechanism whose rate law the guard rejects (`_is_well_posed`) is dropped."

  In `seed_mechanisms`, `foreach(_check_well_posed, seeds)` becomes
  `filter!(_is_well_posed, seeds)`, and its docstring's "Errors when no mechanism binds every
  required regulator, and when a seed it returns gives no defined rate law
  (`_check_well_posed`)." becomes "Errors when no mechanism binds every required regulator;
  a seed whose rate law the guard rejects (`_is_well_posed`) is dropped. The search still
  expands such a node: only the seeds it returns are checked."

  (d) `src/identify_rate_equation.jl`, in `_base_tier`: replace
  ```julia
              children = _expand_re_to_ss(m)
              foreach(_check_well_posed, children)
  ```
  with
  ```julia
              children = filter!(_is_well_posed, _expand_re_to_ss(m))
  ```
  and in its docstring "a flip child whose rate law is not defined (`_check_well_posed`) is
  one, and then none of that seed's flips enters the tier." becomes "a flip child whose rate
  law the guard rejects (`_is_well_posed`) is dropped, as `expand_mechanisms` drops one."

  (e) `src/thermo.jl`: delete `_onlya_violation_message` and `_onlya_haldane_violation(am)`
  with their docstrings. In `_onlya_offenders`' docstring, "The boundary check
  (`_check_well_posed`) names them in its error (`_onlya_violation_message`); the moves read
  only whether any group offends." becomes "The tagging and change-tag moves read only whether
  any group offends; admission (`_admission`) decides what the guard accepts.", and "which the
  boundary check rejects (`_bottomless_segment`)" becomes "which class 3 rejects
  (`_check_face`)".

  (f) `src/names.jl`: in `_edge_constants`' docstring, "`_edge_names` applies it to every edge
  of one conformation, with that conformation's tag, `_column_parameter` picks a joint-solve
  column's, and `_onlya_violation_message` applies it to a group's first edge to label the
  group." becomes "`_edge_names` applies it to every edge of one conformation, with that
  conformation's tag, and `_column_parameter` picks a joint-solve column's."

  (g) `src/dsl.jl`: in the `@enzyme_mechanism` docstring, before the paragraph that begins
  "Below, the inhibitor I binds E and E(S) with one constant", insert:
  ```markdown
  When the macro's code runs, it rejects a mechanism whose rate law is not defined where data
  sit, with an error that names the steps at fault: no cycle of steps turns substrates into
  products; a cycle of rapid-equilibrium steps does, so the rate has no finite value; or, with
  every product, every substrate, one reactant or every modifier absent, the steady state is
  not unique or is not the limit the rate law gives as the rapid-equilibrium steps speed up. A
  law whose denominator vanishes there (0/0) is rejected where it is first derived, by
  [`rate_equation`](@ref) or [`rate_equation_string`](@ref).
  ```
  In the `@allosteric_mechanism` docstring, after the paragraph that begins "A tag says how a
  statement's constants", insert:
  ```markdown
  A tagged step's constants vanish in the conformation that lacks it, and the macro rejects a
  tag set whose limit does not exist or is not unique: an `:OnlyA` binding whose two sides
  untagged steps still join, forms that join free E only through tagged chemistry, or a face
  where the limit depends on how fast the tagged constants vanish. The checks of
  [`@enzyme_mechanism`](@ref) apply as well.
  ```


- [ ] **Step 6: Run the converted files and watch them pass.** The Step 4 commands with
  logs `P4_T8_green_<file>.log`. Expected: `exit 0` for each, every testset passing, no
  warning, no `Error During Test`. Record `test_mechanism_enumeration.jl`'s time from its
  Test Summary (R52: at most 145 s at the branch point) for the commit body.

  Then measure the compile budgets, which the guard's new code in `init_mechanisms`,
  `_minimal_flips` and the derivation entry could break:
  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T8_budget.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_compile_budget.jl")'
  grep -E "init_mechanisms trace-compile|first-call wall-clock|compile reuse" $S/phase45/logs/P4_T8_budget.log
  ```
  Expected: `exit 0`, every budget testset passing; the three `@info` lines give the init
  trace-compile count (budget 200, baseline 57-87), the first `rate_equation` call (budget
  6 s) and the bi-bi→uni-uni warm/cold ratio (budget 0.01). Record the three for the commit
  body. A breach is Step 7's kind 6.

- [ ] **Step 7: Run the full suite and convert what else the wired guard rejects.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T8_full1.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'
  ```
  In the background; poll. `$S/phase45/plan/probes/p3_scan.log` lists every mechanism macro
  block of the S3g test files that phase 3's check accepts and the phase-4 rules reject, by
  file, line and testset; each row names a fixture this step converts. Rerun at 08c4043
  (`$S/phase45/refresh/map-t8-10_2.log`) it finds the same rows, the DSL testset renamed,
  plus TME "_seed_variants: the ordered/random seed" (Step 3 (m)). Read with the mechanism
  texts and quoted blocks the tests hold, and by the stage that added each
  (`$S/phase45/refresh/settle_6.log`), the blocks S3h, S3i and S4 added hold no other
  fixture that phase 3's check accepts and the rules reject. Beyond those rows the wired
  guard changes three testsets:
  - TRG "a graph's DSL text reads back as the graph": `length(allosteric) == 5 + 934 + 5`
    (kind 4, below).
  - TT "every derivation entry rejects a bottomless rapid-equilibrium segment" (kind 2): it
    pins phase 3's message at the macro's code and every derivation entry. Its mechanism
    is Step 3 (a)'s `random_release` (`$S/phase45/refresh/settle_4.log`), whose verdict,
    oracle and fix Step 3 (a) pins: build `m = EnzymeRates.Mechanism(raw)` before the first
    `@test_throws`, pin `(v.fault, v.face) == (EnzymeRates._BOTTOMLESS, 2)` for
    `v = EnzymeRates._guard(m)`, and throw `ErrorException(EnzymeRates._guard_message(m,
    v))` in place of the literal `msg`.
  - TME "expand_mechanisms checks that every child's rate law is defined": the exit drops
    the children it rejected. Replace the `@test_throws` with the exact children it keeps,
    `length(kids) == 2` and `Set(kids) == Set(expected)` for
    `kids = EnzymeRates.expand_mechanisms([parent], rxn)`, the two flips
    `_expand_re_to_ss(parent)` gives (`E(P) + Q <--> E(P, Q)` or `E(Q) + P <--> E(P, Q)`
    steady-state), each written out as `EnzymeRates.Mechanism(@enzyme_mechanism begin ...
    end)` (CLAUDE.md; the macro returns the compiled mechanism and the exit concrete ones,
    as TME writes expected children); the rules
    admit both and reject the 3 dead-end and 19 tagging children, bottomless at zero
    products (`$S/phase45/refresh/settle_4.log`). In its comment, "and `expand_mechanisms`
    rejects them (`_check_well_posed`) instead of returning them" becomes "and
    `expand_mechanisms` drops them".
  TD "a text the macros reject is rejected with their message" keeps passing: its
  bottomless text is `random_release`, and it reads the macro's error at run time and asks
  only for "rapid-equilibrium segment", which the `_BOTTOMLESS` message holds. Every
  failure is one of these kinds, handled so:
  1. A fixture built by a macro that the guard now rejects (the error is the guard's
     message, raised in a macro's expanded code) and that is the subject of a move-level pin
     P4-1 changes (a tagging or change-tag child or parent of a ping-pong, as in the
     testsets of Step 3 (j)): build it with `@_testhelper_unguarded` (same block), add
     `@test EnzymeRates._guard(x).fault === <its fault>` beside it, and leave the pin to
     Task 12; list the testset in the commit body.
  2. Any other fixture the guard now rejects: make the testset a rejection test (the
     unguarded build, its verdict and named step, `@test_throws
     ErrorException(EnzymeRates._guard_message(x, v)) <the macro call>`, and the oracle that
     judges its class: RE-speed for classes 1-3, the ε-path for class 4), plus an admitted
     variant that keeps the testset's subject: the fix the message prints, RES §3.4's or
     §5.2's variant where one is listed, or the probe row's admitted neighbour; check the
     variant with `_testhelper_respeed_verdicts` (all `:agree`) or
     `_testhelper_epath_verdict` (`:converges`) before keeping it.
  3. A law derivation (`fitted_params`, `rate_equation`, `compile_mechanism` and a call,
     `_testhelper_check_law_against_oracle`) of a mechanism the guard rejects: the subject is
     either the move's output (count it with `_independent_param_count`, which runs no guard,
     as in Step 3 (j)) or the law itself (kind 2).
  4. An aggregate count over a population the exit now thins (the population sample's
     `count(Allosteric) >= 30`, the regrouping counts `n_plain >= 20` and `n_allo >= 50`,
     census or init counts): read the new count from the log; if a `>=` bound still holds,
     leave it; if it fails, stop and ask Denis (never lower a bound unasked). An `==` pin
     moves only after a value check attributes the whole difference to the guard:
     - "every bi-bi init mechanism and child is what its raw input builds" (`<n0>`
       children, `<a0>` allosteric at the branch point): rerun Step 1's `P4_bibi_kids.jl`
       on the wired tree (log `P4_T8_bibi_kids.log`). Its children must number
       `<n0> − <K0>` (`<a0> − <A0>` allosteric) plus the superset flips the check inside
       `_minimal_flips` adds (each a guard-passing flip of its parent, as Step 10's
       `is_flip` reads one), and the guard must reject none of them; then re-pin both
       counts and give `<K0>`, `<A0>` and the flips in the commit body;
     - an init tally or the bi-bi init golden: the Gate expects them unchanged, so a change
       is kind 5 unless Step 10's census explanation (it needs only Step 5's code; run it
       first) attributes it to a guard fault;
     - TRG "a graph's DSL text reads back as the graph", `length(allosteric) == 5 + 934 +
       5` (the allosteric children of census reactions 1, 3 and 4): the research rules
       reject 90 of reaction 3's 934, ping-pong tagging children (`island`), and none of
       the other ten, nor any of the 5508 + 5505 of TRG's other population pin
       (`$S/phase45/refresh/settle_5.log`). On the Task 7 tree, count the reaction-3
       children `_guard` rejects, rebuilt as `_testhelper_census_allosteric` builds them;
       each must be `_ISLAND`, and the wired count must be 934 less that count; then re-pin
       `5 + <n3> + 5` and give the count in the commit body;
     - any other `==` pin: stop and ask Denis.
     A difference the check leaves unexplained is kind 5. Ter-ter's count is Step 8's.
  5. A disagreement between the guard and an oracle, or a fixture rejected for a class the
     probe did not report: stop and run `superpowers:systematic-debugging`; never loosen an
     assertion to pass.
  6. A compile budget that fails (`test/test_compile_budget.jl`, measured first in Step 6):
     find which specializations are new (the trace-compile file of a bi-bi
     `init_mechanisms` run before and after this task's src change, filtered to
     `EnzymeRates.`), then stop and ask Denis; never raise a budget unasked (CLAUDE.md).
  Rerun the full suite after the conversions (log `P4_T8_full2.log`) until it is green.
  Expected at the end: `exit 0`, Aqua and JET pass, no warning; the pass count, the time and
  `maxrss_kb` recorded for the commit body.

- [ ] **Step 8: Re-pin the ter-ter count.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T8_terter.log 900 --project=. --heap-size-hint=2500M $S/phase45/scripts/P4_terter.jl
  ```
  Expected: `ter-ter init: <N> mechanisms in <t1> s; the guard rejects 0`, `<t1>` under 150
  and `<N>` near 249,559 (RES §1.10: 1,296 class-3 drops at the exit with the check at the
  exits only; the check inside the seed variants' flips may swap a drop for superset
  variants). The value check of the re-pin: `<N>` is at least `250855 − <R0>` (Step 1), and
  the excess is the superset variants the inside check adds, each passing the guard. Over
  150 s: stop, profile `_faces_pass` inside `_minimal_flips` (RES §4 risk 2: the residual's
  share), and report to Denis before any change. Then, in the TME testset "init_mechanisms on
  ter-ter within 150 s", replace `@test length(ms) == 250855` with
  `@test length(ms) == <N>` and the comment's first sentence with "Aggregate pin over the
  whole ter-ter seed set: its seeds and their merged and Theorell–Chance variants that the
  guard admits (<D> fewer than the bottomless test alone kept)." where `<D>` is
  250,855 − `<N>` (class 3 drops variants at the exit and, inside the flips, may swap one
  for its supersets); run the testset's file
  focused (log `P4_T8_tme_terter.log`, background) and expect it to pass.


- [ ] **Step 9: Update ROADMAP.md and the user-facing roadmap.**

  In `ROADMAP.md`, replace the entry "**`BigInt` arithmetic inflates the GC's live-byte
  count.**" (its whole text) with:
  ```markdown
  - **`BigInt` arithmetic inflates the GC's live-byte count.**
    Three functions that compute in `BigInt` raise `Base.gc_live_bytes` by bytes that no
    live object holds: `_onlya_offenders` (`src/thermo.jl`; a `Rational{BigInt}` null
    space), the tagging and change-tag moves' `:OnlyA` test, by about 15 KiB per call; the
    counter that `_split_gain_test` uses for a plain mechanism
    (`_partition_independent_count`; a `Rational{BigInt}` row reduction) by 18 KiB; and
    `_eq_complexity` (a `BigInt` determinant, skipped for a single rapid-equilibrium
    segment) by 15 KiB. A 6×6 determinant outside the package does the same in `BigInt`
    and `Rational{BigInt}`, not in `Rational{Int}`. The moves run the `:OnlyA` test on
    every tag set they try: 100,000 calls add 1.46 GiB and 12 full sweeps. The joint solve
    adds nothing, on a plain mechanism or on a concrete allosteric one; it runs in
    `Rational{Int}` on 10,226 enumerated bi-bi and ter-bi mechanisms without overflow. The
    count stays up after a full collection: 100,000 calls that start 1.47 GiB high make 52
    full sweeps in 10.0 s, against 12 in 6.5 s from a clean start. The count holds no
    memory and raises no peak: a process whose count reaches 2.8 GiB peaks at 1.2 GB RSS,
    and the enumeration tests' count after a full collection stays near 1.8 GiB whether
    their peak RSS is 2.2 or 3.5 GB. The higher peak came from mechanisms that kept vectors
    at the capacity their builds grew to; `_canonicalize` stores each at its length, and
    `test/test_reaction_graph.jl` pins it. Filtering the moves through admission replaces
    the `:OnlyA` test, not the other two. A loop that calls one of these and prints
    `Base.gc_live_bytes()` and `Base.gc_num().full_sweep` shows it. The figures come from a
    random-order bi-bi of nine steady-state steps, plain and as a two-subunit allosteric
    mechanism with `:OnlyA` on both bindings of A and on the chemical step and `:EqualAI`
    elsewhere (Julia 1.12.5, aarch64, `--heap-size-hint=2500M`).
  ```
  and replace the entry "**The `:OnlyA` guard misjudges groups that mix a binding and a
  chemistry step.**" with:
  ```markdown
  - **The tagging and change-tag moves misjudge `:OnlyA` groups that mix a binding and a
    chemistry step.** Their filter, `_onlya_offenders` (`src/thermo.jl`), skips every
    `:OnlyA` group that holds a chemistry step, including one that also holds a plain
    binding of the same metabolite (enumeration builds these by merging and tags every
    allosteric one `:OnlyA`), and gives the bindings of one group a single rate of
    divergence. The moves therefore emit tag sets that admission rejects, such as a mixed
    `:OnlyA` group whose plain member sits on a square of `:EqualAI` bindings (the exit
    drops them), and miss tag sets that admission admits, such as a split `:OnlyA` binding
    group whose limit exists. Filtering the moves through admission replaces the test.
  ```
  In `docs/src/roadmap.md`, under "## Support more mechanisms", after the paragraph on iso
  mechanisms, add:
  ```markdown
  **Allosteric ping-pong that stops both half-reactions.** An allosteric ping-pong may stop
  one half-reaction in the inactive conformation, as phosphorylation stops pyruvate
  dehydrogenase's E1 and leaves E2 and E3 running. Stopping both strands the covalent
  intermediate between the halves: its weight, and with it the rate, then depends on how fast
  each half stops, so `@allosteric_mechanism` rejects such a mechanism. Supporting it needs
  the conformations to interconvert at the covalent intermediate as well, which will come
  with iso mechanisms.
  ```


- [ ] **Step 10: Build the wired tree's census and explain every difference.**

  Create `$S/phase45/scripts/P4_census_explain.jl`:
  ```julia
  # ABOUTME: Explains a census difference between two corpora of the guard worktree: every
  # ABOUTME: child, seed or law only in the first is rebuilt and tallied by the guard's fault.
  # Usage: julia --project=<census env of the second tree> P4_census_explain.jl <A> <B>
  import EnzymeRates
  const ER = EnzymeRates
  const A_DIR, B_DIR = ARGS[1], ARGS[2]

  rows(path) = [split(l, '\t') for l in eachline(path) if !isempty(l)]
  oriented = Dict{String, String}()
  for dir in (A_DIR, B_DIR), r in rows(joinpath(dir, "mechs.tsv"))
      oriented[r[1]] = r[3]
  end
  # The names a reaction's text lists after `key:` (`role_names`, not `names`, which
  # Base exports).
  role_names(text, key) = (m = match(Regex("$key:\\s*([^\\n\\\\]*)"), text);
                           m === nothing ? Symbol[] :
                           [Symbol(strip(split(x, "[")[1])) for x in split(m[1], ",")])
  roles = Dict(String(r[1]) => (role_names(r[end], "substrates"),
                                role_names(r[end], "products"))
               for r in rows(joinpath(A_DIR, "reactions.tsv")))

  """The mechanism of census oriented text `t` (README "Neutral text"), built by the
  constructors, for the reaction with substrates `subs` and products `prods`."""
  function rebuild(t::AbstractString, subs, prods)
      reactant(x) = x in subs || x in prods
      function met(l)
          x = Symbol(l)
          reactant(x) && return ER.Met(x, false)
          endswith(l, "mod") && reactant(Symbol(l[1:end-3])) &&
              return ER.Met(Symbol(l[1:end-3]), true)
          ER.Met(x, true)
      end
      function form(s)
          i = findfirst('(', s)
          conf, inner = Symbol(s[1:i-1]), s[i+1:end-1]
          parts = split(inner, "; residual=")
          bound = isempty(parts[1]) ? ER.Met[] :
                  ER.Met[met(lab) for lab in split(parts[1], ",")]
          res = length(parts) > 1 ? split(parts[2], "_") : SubString{String}[]
          added = Symbol[Symbol(r[5:end]) for r in res if startswith(r, "plus")]
          subtracted = Symbol[Symbol(r[6:end]) for r in res if startswith(r, "minus")]
          ER.RawForm(conf, bound, added, subtracted)
      end
      allo = startswith(t, "allosteric ")
      n, body, sitetext = 1, t[length("plain ") + 1:end], "none"
      if allo
          mm = match(r"^allosteric n=(\d+) (.*) sites (.*)$", t)
          n, body, sitetext = parse(Int, mm[1]), mm[2], mm[3]
      end
      steps, group, tags = ER.RawStep[], Int[], Symbol[]
      for (g, gt) in enumerate(split(body, " / "))
          mm = match(r"^\{(.*)\}(?::(\w+))?$", gt)
          push!(tags, mm[2] === nothing ? :EqualAI : Symbol(mm[2]))
          for st in split(mm[1], " | ")
              l, r = split(st[1:end-3], " -> ")
              ls, rs = split(l, " + "), split(r, " + ")
              push!(steps, ER.RawStep(form(ls[1]), form(rs[1]),
                                      ER.Met[met(lab) for lab in ls[2:end]],
                                      ER.Met[met(lab) for lab in rs[2:end]],
                                      endswith(st, " RE")))
              push!(group, g)
          end
      end
      mods = sort!(unique(Symbol[x.name for s in steps for x in [s.consumed; s.released]
                                 if x.modifier]))
      raw = ER.RawMechanism(sort(subs), sort(prods), mods, steps, group)
      allo || return ER.Mechanism(raw)
      sites = ER.RawSite[]
      if sitetext != "none"
          for st in split(sitetext, " / ")
              mm = match(r"^site\((\d+): (.*)\)$", st)
              ligs = [split(x, "::") for x in split(mm[2], ", ")]
              push!(sites, ER.RawSite(Symbol.(first.(ligs)), parse(Int, mm[1]),
                                      Symbol.(last.(ligs))))
          end
      end
      ER.AllostericMechanism(ER.RawAllostericMechanism(raw, tags, n, sites))
  end

  fault(rx, id) = (m = rebuild(oriented[id], roles[rx]...); ER._guard(m).fault)
  """Whether child `c` is a flip of parent `p`: the same edges in the same groups, every
  rapid-equilibrium flag of `c` one of `p`'s."""
  function is_flip(rx, c, p)
      a, b = (ER.catalytic(rebuild(oriented[x], roles[rx]...)) for x in (c, p))
      key(m, k) = (m.forms[m.edges[k].from], m.forms[m.edges[k].to], m.edges[k].consumed,
                   m.edges[k].released)
      flags(m) = Dict(key(m, k) => e.rapid_equilibrium for (k, e) in enumerate(m.edges))
      groups(m) = Set(Set(key(m, k) for k in ks) for ks in ER._group_members(m))
      fa, fb = flags(a), flags(b)
      keys(fa) == keys(fb) && groups(a) == groups(b) && all(k -> fa[k] <= fb[k], keys(fa))
  end
  tally!(d, k) = (d[k] = get(d, k, 0) + 1)

  ids(dir) = Set(r[1] for r in rows(joinpath(dir, "mechs.tsv")))
  idsB = ids(B_DIR)
  children(dir) = (d = Dict{Tuple{String, String}, Set{String}}();
                   foreach(r -> push!(get!(Set{String}, d, (r[1], r[3])), r[5]),
                           rows(joinpath(dir, "children.tsv"))); d)
  ca, cb = children(A_DIR), children(B_DIR)
  childB = Set(c for (_, cs) in cb for c in cs)
  unexplained = 0
  # An item only in A is explained when the guard rejects it, or when B still holds it but
  # a quota picked other mechanisms (a superset flip can take a smaller id's place).
  explained(rx, id) = fault(rx, id) !== ER._NO_FAULT || id in idsB
  for rx in sort(collect(keys(roles)))
      drops, flips, differ = Dict{Any, Int}(), 0, 0
      for key in sort(collect(keys(ca)))
          key[1] == rx && haskey(cb, key) || continue
          a, b = ca[key], cb[key]
          a == b && continue
          differ += 1
          for c in setdiff(a, b)
              f = fault(rx, c)
              tally!(drops, f)
              f === ER._NO_FAULT && (global unexplained += 1;
                                     println("UNEXPLAINED drop $rx $(key[2]) -> $c"))
          end
          for c in setdiff(b, a)
              fault(rx, c) === ER._NO_FAULT && is_flip(rx, c, key[2]) ? (flips += 1) :
                  (global unexplained += 1; println("UNEXPLAINED child $rx $(key[2]) -> $c"))
          end
      end
      # laws.tsv keys a row by law and names its reaction in column 2; the others in 1.
      col(dir, file, j) = Set(r[j] for r in rows(joinpath(dir, file))
                              if r[file == "laws.tsv" ? 2 : 1] == rx)
      out = Dict{String, Any}()
      for (label, file, j) in (("seeds", "seeds.tsv", 3), ("init", "init.tsv", 3),
                               ("parents", "parents.tsv", 3), ("laws", "laws.tsv", 3))
          only_a = setdiff(col(A_DIR, file, j), col(B_DIR, file, j))
          only_b = setdiff(col(B_DIR, file, j), col(A_DIR, file, j))
          t = Dict{Any, Int}()
          for id in only_a
              tally!(t, fault(rx, id))
              explained(rx, id) || (global unexplained += 1;
                                    println("UNEXPLAINED $label only in A: $rx $id"))
          end
          for id in only_b
              fault(rx, id) === ER._NO_FAULT && continue
              global unexplained += 1
              println("UNEXPLAINED $label only in B: $rx $id")
          end
          label in ("seeds", "init") && !isempty(only_b) &&
              (global unexplained += length(only_b);
               println("UNEXPLAINED new $label: $rx $(collect(only_b))"))
          out[label] = (only_a = t, only_b = length(only_b))
      end
      println(rx, ": ", differ, " parents differ; child drops ", drops, "; superset flips ",
              flips, "; ", out)
  end
  println("unexplained: ", unexplained)
  ```
  Then:
  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  $S/phase3/jl.sh $S/phase45/logs/P4_T8_census.log 1800 --project=$S/phase3/census/env_guard --heap-size-hint=2500M $S/phase3/census/corpus.jl $S/phase45/guard_wt $S/phase3/census/P4_T8 $S/phase3/census/rename_phase3.toml
  $S/phase3/jl.sh $S/phase45/logs/P4_T8_census_cmp.log 600 --project=$S/phase3/census/env --heap-size-hint=2500M $S/phase3/census/compare.jl $S/phase3/census/P4_base $S/phase3/census/P4_T8 $S/phase3/census/rename_phase3.toml
  $S/phase3/jl.sh $S/phase45/logs/P4_T8_census_explain.log 900 --project=$S/phase3/census/env_guard --heap-size-hint=2500M $S/phase45/scripts/P4_census_explain.jl $S/phase3/census/P4_base $S/phase3/census/P4_T8
  ```
  Expected: the corpus `exit 0`; compare.jl exits 1 with `RESULT FAIL` whose differences are
  child sets (pingpong, prefix and terbi parents with ping-pong tagging children, bibi level-2
  parents with mixed `:OnlyA` groups), the pingpong seeds and the law sample, and no
  stored-direction, init, value or law-string difference on a law both corpora sample; the
  explain log prints per reaction a drop tally whose faults are `_ISLAND`, `_NO_LIMIT`,
  `_SPLIT` and `_TRAP` (no `_NO_FAULT`), superset flips at about one parent per reaction, no
  init drop, no new seed, and ends `unexplained: 0`. A positive count is a failure: run
  `superpowers:systematic-debugging` on the first unexplained line before going on. Copy the
  tallies into the commit body.


- [ ] **Step 11: Run the full suite.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T8_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'
  ```
  In the background; poll. Expected: `exit 0`, `Testing EnzymeRates tests passed`, Aqua and
  JET pass, no warning after "Testing Running tests"; `rate_equation` performance gates and
  the compile budgets pass (they run in the suite); `maxrss_kb` within 5% of Task 1's M0, or
  the rise explained in the commit body.

- [ ] **Step 12: Commit.** First `git add`, by name, every other test file Step 7 converted
  (`git status --short` lists them; never `git add -A`), and write the measured figures and
  lists in place of the `<...>` below. Then:

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && git add src/admission.jl src/rate_eq_derivation.jl src/mechanism_enumeration.jl src/identify_rate_equation.jl src/thermo.jl src/names.jl src/dsl.jl ROADMAP.md docs/src/roadmap.md test/test_admission.jl test/test_types.jl test/test_dsl.jl test/test_rate_eq_derivation.jl test/allosteric_ground_truth.jl test/test_oracles.jl test/test_mechanism_enumeration.jl test/test_identify_rate_equation.jl && git status --short && git commit -F - <<'EOF2'
  Reject what the guard rejects at the macros, the parser and the derivation

  The macros, the parser and every derivation entry throw the guard's message; the
  enumerator's exits and identify's base tier drop what it rejects; class 3 runs inside the
  flip move's predicate in place of the bottomless test; class 5 runs where the law is
  derived. The :OnlyA Haldane guard, the split's pre-check and balancers, and the
  derive-time RE-cycle error go. The one-flip check stays until the flip moves to free E.

  Census against the starting tree: <the explain log's per-reaction tallies>; unexplained 0.
  Ter-ter init: <t0> s before, <t1> s after, <N> mechanisms (<R0> rejected by the guard
  before). test_mechanism_enumeration.jl focused: <te> s. Compile budgets: init
  trace-compile <n>, first rate_equation <s> s, bi-bi to uni-uni reuse <r>. Full suite:
  <pass count>, <time>, maxrss <kb> KB (branch point <kb0> KB).

  Converted to rejection tests with admitted variants: <list of testsets>. Move-level pins
  that P4-1 changes, built unguarded until Task 12: <list>. Re-pinned after a value check:
  <each pin, its old and new value, and the check (Step 7 kind 4, Step 8)>. Subjects removed for Denis's
  approval (P4-29, P4-29a/b): the split pre-check and balancer branch; the dead-inactive
  ping-pong's d_free_I == 1 and its fitted-name pin; the two "no covalent sink" checks,
  since the half-2 variants hold their inactive enzyme in the covalent forms at zero
  products (rate 0 there, the mass-action limit too); the RE-chemistry ground-truth
  harness now gates nothing.

  Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
  EOF2
  ```


---

### Task 9: Flip at free E (RES §4 row 8)

The conformations interconvert at free E, form 1, alone: the inactive conformation keeps the
forms its untagged edges join to E, and `d_free` is E's own weight. Admission's face condition
(Task 5) now catches what the one-flip check guarded (R2, AM risk 1), so the one-flip check
goes. Phase 3 already roots E's segment at form 1 (`_segment_root` ranks it first), so a
second free form in rapid equilibrium with E keeps its law, with L meaning [E_I]/[E_A]
(P4-21).

**Files:**
- Modify: `src/types.jl` (`_inactive_mask`; `_outside_onlya` deleted when Step 0 finds no
  other caller)
- Modify: `src/rate_eq_derivation.jl` (`_raw_symbolic_rate_polys`: `d_free` at form 1)
- Modify: `src/admission.jl` (delete `_assert_one_flip_segment` and
  `_free_forms_flip_together`; `_check_well_posed`, `_is_well_posed`)
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`
  (`_testhelper_flip_form` deleted; `_testhelper_steady_state_rate`'s `flip` keyword; new
  `_testhelper_check_law_against_limit`)
- Modify: `ROADMAP.md` (retire Planned "Take `:OnlyA` steps to their limit in the inactive
  conformation")
- Test: `test/test_oracles.jl`, `test/test_types.jl`, `test/test_admission.jl`,
  `test/test_identify_rate_equation.jl`; comments and `_testhelper_inactive_reach` in
  `test/test_reaction_graph.jl`

**Interfaces:**
- Consumes: Task 2's `_law_mask(am, c)`; Task 1's `_testhelper_protomer_rate(m, point,
  concs; n)`, `_testhelper_consistent_point` (with `K`, `kf`), `_testhelper_law_value`,
  `_testhelper_network`, `_testhelper_re_limit`, `_testhelper_draw_speeds`,
  `_testhelper_single_subunit` (concrete methods), `_testhelper_respeed_verdicts`,
  `_testhelper_epath_verdict`; Task 3's `_testhelper_face_faults`; Task 4's
  `_testhelper_tag_sets` (test_admission.jl); Task 5's `_testhelper_admission_agreement`;
  phase 3's `_testhelper_check_law_against_oracle`.
- Produces:
  ```julia
  _inactive_mask(am::AllostericMechanism)::BitVector   # _law_mask(am, 2): pruned from form 1
  _check_well_posed(m::Union{Mechanism, AllostericMechanism})::Nothing  # the guard alone
  _is_well_posed(m::Union{Mechanism, AllostericMechanism})::Bool        # one method
  _testhelper_steady_state_rate(m, point, concs::AbstractDict; Λ = big(10)^30,
                                flip::Int = 1)::BigFloat                # fixtures file
  _testhelper_check_law_against_limit(m; n_points::Int = 3, seed::Int = 1)  # fixtures file:
      # the law against the RE limit of the oracle's own network, which prunes by itself
  ```
  `_testhelper_check_law_against_oracle` reads the package's masks for its network
  (`point.inactive`, from `_inactive_mask`), so it cannot see a pruning error: the second
  free forms that this task's pruning strands (R2P_sealed, B_uphill) and Task 10's `:OnlyI`
  fixtures are checked against `_testhelper_network`, whose pruning is its own (P4-24).

**Gate:** the iso families (`round4/flipE/scripts/exact.jl`: IsoUU steady-state and RE
bindings × 8 tag sets, IsoBiBi × 9, `D ⇌ E` × 2) admitted by the macros and exact at n = 1
against the two-conformation solver (1e-20), three of them at n = 2 and one at n = 3 against
the concerted protomer oracle (1e-20); a flip at another free form gives another law; the
package law at L equals the oracle flipping at D at L·[D_I]/[D_A]·[E_A]/[E_I] (P4-21); R1, its
D-label twin, A5, AD2a, AD2d, AD5, AD7a, AD7c, R2P_sealed, B_uphill and AD13 × 10 admitted with
law checks at 1e-20 against the oracle's own pruning (the Task 4-5 rows with a second free form,
which used the oracles' own law until now); S1 rejected by class 3 at three faces with both
variants admitted by both oracles; R2 `_PUMPED` with AD2d admitted; R3, R4
`_ISLAND` with their P4-1 bases admitted; TID's base tier keeps all six flips; census `P4_T9`
equals `P4_T8` (compare.jl `RESULT PASS`, zero differences); full suite green. The probe
measured every verdict (`p3_verdicts.log`: the iso families, D cases, R1-R4, A5, S1 and its
variants, the bases pass or fail exactly as listed).

- [ ] **Step 0: Check the names this task reads.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  grep -n "_outside_onlya\|_assert_one_flip_segment\|_free_forms_flip_together\|_testhelper_flip_form" src/*.jl test/*.jl
  grep -n "i_free\|d_free" src/rate_eq_derivation.jl
  grep -n "interconvert at one free form\|interconvert at a free form\|must be in rapid equilibrium with\|reach from a free form\|from any free form" src/*.jl test/*.jl
  ```
  Expected: `_outside_onlya` called by `_inactive_mask` and `_free_forms_flip_together` only
  (delete it in Step 3; if another caller exists, keep it and say so in the commit body);
  `_assert_one_flip_segment` and `_free_forms_flip_together` in admission.jl and the
  one-flip testsets only; `_testhelper_flip_form` in `_testhelper_steady_state_rate` only;
  `i_free` in `_raw_symbolic_rate_polys` only. The last grep lists the comments and
  docstrings that state the one-flip rule or prune from any free form; Step 3 (e) rewrites
  each that this task makes false.

- [ ] **Step 1: Write the failing tests.**

  (a) Append to `test/test_oracles.jl`:
  ```julia
  """
  The iso-mechanism families of the flip-at-free-E research (`round4/flipE/scripts/exact.jl`)
  as `@allosteric_mechanism` text: Iso Uni Uni (`:uu_ss`, steady-state bindings; `:uu_re`,
  rapid-equilibrium bindings) and iso ordered Bi Bi (`:bibi`), the tags in step order, and
  `extra` steps appended.
  """
  function _testhelper_iso_text(family::Symbol, tags; n::Int = 1, extra::String = "")
      steps = family === :bibi ?
          ["E + A ⇌ E(A)", "E(A) + B <--> E(A, B)", "E(A, B) <--> Eprime(P, Q)",
           "Eprime(Q) + P <--> Eprime(P, Q)", "Eprime + Q <--> Eprime(Q)", "Eprime <--> E"] :
          family === :uu_ss ?
          ["E + S <--> E(S)", "E(S) <--> Eprime(P)", "Eprime + P <--> Eprime(P)",
           "Eprime <--> E"] :
          ["E + S ⇌ E(S)", "E(S) <--> Eprime(P)", "Eprime + P ⇌ Eprime(P)",
           "Eprime <--> E"]
      reactants = family === :bibi ? "substrates: A, B\nproducts: P, Q" :
                                     "substrates: S\nproducts: P"
      lines = join(["$s :: $t" for (s, t) in zip(steps, tags)], "\n")
      "@allosteric_mechanism begin\n$reactants\ncatalytic_multiplicity: $n\n" *
          "steps: begin\n$lines\n$extra\nend\nend"
  end

  @testset "the conformations flip at free E" begin
      E_, N_, O_ = :EqualAI, :NonequalAI, :OnlyA
      uu_tags = ((E_, E_, E_, N_), (E_, N_, E_, E_), (E_, O_, E_, E_), (N_, E_, E_, E_),
                 (E_, E_, E_, O_), (E_, O_, E_, N_), (N_, N_, N_, N_), (O_, O_, E_, E_))
      bibi_tags = ((E_, E_, E_, E_, E_, N_), (E_, E_, N_, E_, E_, E_),
                   (E_, E_, O_, E_, E_, E_), (E_, N_, E_, E_, E_, E_),
                   (E_, E_, E_, E_, E_, O_), (N_, E_, E_, E_, E_, E_),
                   (E_, E_, E_, N_, E_, N_), (O_, E_, O_, E_, E_, E_),
                   (N_, N_, N_, N_, N_, N_))
      texts = [[_testhelper_iso_text(:uu_ss, t) for t in uu_tags];
               [_testhelper_iso_text(:uu_re, t) for t in uu_tags];
               [_testhelper_iso_text(:bibi, t) for t in bibi_tags];
               _testhelper_iso_text(:uu_ss, (E_, E_, E_, E_); extra = "D ⇌ E :: NonequalAI");
               _testhelper_iso_text(:uu_ss, (E_, N_, E_, E_); extra = "D ⇌ E :: NonequalAI")]
      # Each family member at one subunit equals the exact two-conformation steady state,
      # the conformations flipping at free E only. The macro builds every one: a second
      # free form needs no rapid-equilibrium link to E.
      for (i, text) in enumerate(texts)
          @testset "family member $i" begin
              _testhelper_check_law_against_oracle(eval(Meta.parse(text)); seed = i)
          end
      end
      # More subunits: the concerted oracle, in which only the fully unliganded oligomer
      # flips, at free E.
      for (i, (family, tags, n)) in enumerate(((:uu_ss, (E_, E_, E_, N_), 2),
                                               (:uu_ss, (E_, O_, E_, N_), 2),
                                               (:bibi, (E_, E_, E_, N_, E_, N_), 2),
                                               (:uu_ss, (E_, E_, E_, N_), 3)))
          m = eval(Meta.parse(_testhelper_iso_text(family, tags; n)))
          am = EnzymeRates.AllostericMechanism(m)
          fitted = EnzymeRates._dependent_param_exprs(typeof(m))[2]
          rng = Random.MersenneTwister(100 + i)
          point = setprecision(() -> _testhelper_consistent_point(am, rng), BigFloat, 256)
          concs = Dict{Symbol, BigFloat}(x => exp(3 * big(rand(rng)) - big(1.5))
                                         for x in EnzymeRates.metabolites(am))
          params = Dict{Symbol, BigFloat}(k => point.vals[k] for k in fitted)
          params[:Keq] = point.Keq
          exact = setprecision(() -> _testhelper_protomer_rate(am, point, concs; n),
                               BigFloat, 256)
          @test abs(_testhelper_law_value(m, params, concs) - exact) <=
                big(10)^-20 * abs(exact)
      end
      # A flip at the other free form gives another law: the law tells the flip form.
      m = eval(Meta.parse(_testhelper_iso_text(:uu_ss, (E_, E_, E_, N_))))
      am = EnzymeRates.AllostericMechanism(m)
      eprime = only(f for f in 2:length(am.catalytic.forms)
                    if EnzymeRates._is_free(am.catalytic.forms[f]))
      fitted = EnzymeRates._dependent_param_exprs(typeof(m))[2]
      rng = Random.MersenneTwister(201)
      for _ in 1:3
          point = setprecision(() -> _testhelper_consistent_point(am, rng), BigFloat, 256)
          concs = Dict{Symbol, BigFloat}(x => exp(3 * big(rand(rng)) - big(1.5))
                                         for x in EnzymeRates.metabolites(am))
          params = Dict{Symbol, BigFloat}(k => point.vals[k] for k in fitted)
          params[:Keq] = point.Keq
          law = _testhelper_law_value(m, params, concs)
          at_e = _testhelper_steady_state_rate(am, point, concs)
          at_eprime = _testhelper_steady_state_rate(am, point, concs; flip = eprime)
          @test abs(law - at_e) <= big(10)^-20 * abs(at_e)
          @test abs(law - at_eprime) > big(10)^-4 * abs(at_e)
      end
      # A free form D in rapid equilibrium with E (D's label sorts before E's) changes no
      # law: flipping at D with ratio L_D is flipping at E with ratio L_D·r_I/r_A, where
      # r_c = [D_c]/[E_c] (P4-21). So L means [E_I]/[E_A].
      m = eval(Meta.parse(_testhelper_iso_text(:uu_ss, (E_, E_, E_, E_);
                                               extra = "D ⇌ E :: NonequalAI")))
      am = EnzymeRates.AllostericMechanism(m)
      cm = am.catalytic
      d = only(f for f in 2:length(cm.forms)
               if EnzymeRates._is_free(cm.forms[f]) &&
                  cm.conformations[cm.forms[f].conformation] === :D)
      kd = only(j for (j, e) in enumerate(cm.edges) if Set((e.from, e.to)) == Set((1, d)))
      fitted = EnzymeRates._dependent_param_exprs(typeof(m))[2]
      rng = Random.MersenneTwister(202)
      for _ in 1:3
          point = setprecision(() -> _testhelper_consistent_point(am, rng), BigFloat, 256)
          concs = Dict{Symbol, BigFloat}(x => exp(3 * big(rand(rng)) - big(1.5))
                                         for x in EnzymeRates.metabolites(am))
          params = Dict{Symbol, BigFloat}(k => point.vals[k] for k in fitted)
          params[:Keq] = point.Keq
          ratio(c) = cm.edges[kd].from == d ? 1 / point.K[kd, c] : point.K[kd, c]
          at_d = _testhelper_steady_state_rate(am, merge(point, (L = point.L * ratio(2) /
                                                                 ratio(1),)),
                                               concs; flip = d)
          @test abs(_testhelper_law_value(m, params, concs) - at_d) <=
                big(10)^-20 * abs(at_d)
      end
  end
  ```

  (b) Append to `test/test_admission.jl`:
  ```julia
  @testset "a second free form keeps its law exact once the conformations flip at E" begin
      # Every admitted tag set of `_testhelper_tag_sets` that holds a free form besides E:
      # R1, AD2a, AD2d, AD5, R2P_sealed and B_uphill join it to E by a steady-state or
      # tagged step, A5, AD7a, AD7c and the AD13 family in rapid equilibrium (its label
      # sorting before or after E). The conformations flip at E alone, so every one
      # derives through its derivation entry, and each law equals the exact limit of the
      # oracle's own network, which prunes from free E by itself (R2P_sealed and B_uphill
      # lose their second free form in the inactive conformation); Tasks 4 and 5 judged
      # these by the oracles' own deletion law.
      tag_sets = [["R1", "AD2a", "AD2d", "AD5", "R2P_sealed", "B_uphill", "A5", "AD7a",
                   "AD7c"];
                  ["AD13_$(chem)_$label" for chem in (:NonequalAI, :OnlyA)
                   for label in (:D, :Ca, :Ea, :F, :Z)]]
      for (i, name) in enumerate(tag_sets)
          @testset "$name" begin
              _testhelper_check_law_against_limit(
                  EnzymeRates.compile_mechanism(_testhelper_tag_sets[name]); seed = i)
          end
      end
  end
  ```

  (b′) At the end of the "Phase-4 oracles" section of
  `test/mechanism_definitions_for_test_enzyme_derivation.jl` add:
  ```julia
  """
  Check that the rate law of `m`, with one catalytic subunit, equals to 1e-20 the exact
  rapid-equilibrium limit of the oracle's own law network (`_testhelper_network`,
  `_testhelper_re_limit`), which deletes each conformation's excluded steps and prunes from
  free E by itself, at `n_points` consistent points (`_testhelper_consistent_point`). A
  compiled `m` is read through its derivation entry, which runs the guard; a concrete `m` is
  read without it (`_testhelper_law_value`).
  """
  function _testhelper_check_law_against_limit(m; n_points::Int = 3, seed::Int = 1)
      sub = _testhelper_single_subunit(m)
      am = _testhelper_lift(sub)
      fitted = EnzymeRates._dependent_param_exprs(
          sub isa EnzymeRates.AbstractEnzymeMechanism ? typeof(sub) : sub)[2]
      rng = Random.MersenneTwister(seed)
      setprecision(BigFloat, 256) do
          for _ in 1:n_points
              point = _testhelper_consistent_point(am, rng)
              concs = Dict{Symbol, BigFloat}(x => exp(3 * big(rand(rng)) - big(1.5))
                                             for x in EnzymeRates.metabolites(am))
              params = Dict{Symbol, BigFloat}(k => point.vals[k] for k in fitted)
              params[:Keq] = point.Keq
              speeds = _testhelper_draw_speeds(am, rng)
              exact = _testhelper_re_limit(_testhelper_network(am, point, concs; speeds))
              @test abs(_testhelper_law_value(sub, params, concs) - exact) <=
                    big(10)^-20 * abs(exact)
          end
      end
  end
  ```

  (c) `test/test_types.jl`: replace the testset "a stray free form cannot leave the inactive
  state's segment bottomless" with:
  ```julia
      @testset "a stray free form cannot leave the inactive state's segment bottomless" begin
          # Every step that takes up the substrate A is :OnlyA, and the conformations flip at
          # free E, so the inactive conformation keeps E and E(P) alone: a second free form
          # F, which binds B and joins E only through F + B <--> E(B), holds no inactive
          # enzyme. With every substrate absent, or A, or B, the forms F, E(B), E(A::Mod)
          # and E(A::Mod, B) hold the enzyme they receive and never return it to E: a second
          # closed class, which class 3 rejects (`_NONUNIQUE`) and the RE-speed oracle sees.
          s1 = @_testhelper_unguarded @allosteric_mechanism begin
              substrates: A, B
              products: P
              modifiers: A
              steps: begin
                  E + A <--> E(A)                   :: OnlyA
                  E(A) + B <--> E(A, B)             :: EqualAI
                  E(A, B) <--> E(P)                 :: OnlyA
                  E(P) <--> E + P                   :: EqualAI
                  E + A::Mod ⇌ E(A::Mod)            :: OnlyA
                  E(A::Mod) + B ⇌ E(A::Mod, B)      :: EqualAI
                  E(B) + A::Mod ⇌ E(A::Mod, B)      :: EqualAI
                  F + B <--> E(B)                   :: EqualAI
              end
          end
          @test _testhelper_face_faults(s1)[3:5] == fill(ER._NONUNIQUE, 3)
          v = ER._guard(s1)
          @test (v.fault, v.face) == (ER._NONUNIQUE, 3)
          @test _testhelper_respeed_verdicts(s1; rng = Random.MersenneTwister(12),
                                             law = :limit)[3:5] == fill(:nonunique, 3)
          @test_throws ErrorException(ER._guard_message(s1, v)) @allosteric_mechanism begin
              substrates: A, B
              products: P
              modifiers: A
              steps: begin
                  E + A <--> E(A)                   :: OnlyA
                  E(A) + B <--> E(A, B)             :: EqualAI
                  E(A, B) <--> E(P)                 :: OnlyA
                  E(P) <--> E + P                   :: EqualAI
                  E + A::Mod ⇌ E(A::Mod)            :: OnlyA
                  E(A::Mod) + B ⇌ E(A::Mod, B)      :: EqualAI
                  E(B) + A::Mod ⇌ E(A::Mod, B)      :: EqualAI
                  F + B <--> E(B)                   :: EqualAI
              end
          end
          # Two variants keep F and pass every face (RES §5.2 item 1): F <--> E tagged :OnlyA;
          # or F <--> E untagged with the modifier binding of A untagged, so that no :OnlyA
          # binding has its two sides joined by untagged steps.
          f_onlya = @allosteric_mechanism begin
              substrates: A, B
              products: P
              modifiers: A
              steps: begin
                  E + A <--> E(A)                   :: OnlyA
                  E(A) + B <--> E(A, B)             :: EqualAI
                  E(A, B) <--> E(P)                 :: OnlyA
                  E(P) <--> E + P                   :: EqualAI
                  E + A::Mod ⇌ E(A::Mod)            :: OnlyA
                  E(A::Mod) + B ⇌ E(A::Mod, B)      :: EqualAI
                  E(B) + A::Mod ⇌ E(A::Mod, B)      :: EqualAI
                  F + B <--> E(B)                   :: EqualAI
                  F <--> E                          :: OnlyA
              end
          end
          copy_untagged = @allosteric_mechanism begin
              substrates: A, B
              products: P
              modifiers: A
              steps: begin
                  E + A <--> E(A)                   :: OnlyA
                  E(A) + B <--> E(A, B)             :: EqualAI
                  E(A, B) <--> E(P)                 :: OnlyA
                  E(P) <--> E + P                   :: EqualAI
                  E + A::Mod ⇌ E(A::Mod)            :: EqualAI
                  E(A::Mod) + B ⇌ E(A::Mod, B)      :: EqualAI
                  E(B) + A::Mod ⇌ E(A::Mod, B)      :: EqualAI
                  F + B <--> E(B)                   :: EqualAI
                  F <--> E                          :: EqualAI
              end
          end
          # Both oracles agree on every face: the RE-speed oracle (class 3) and the ε-path
          # (admission of the :OnlyA tags).
          for (i, variant) in enumerate((f_onlya, copy_untagged))
              @test all(_testhelper_admission_agreement(
                  ER.AllostericMechanism(variant); rng = Random.MersenneTwister(12 + i)))
          end
      end
  ```
  Replace the testset "allosteric free forms share one rapid-equilibrium segment" with:
  ```julia
      @testset "a second free form is admitted unless a tag strands it" begin
          # The conformations flip at free E alone. A second free form joined to E by a
          # steady-state step is an ordinary form, whatever its label (R1, and its twin
          # labelled D, which sorts before E). Where tags strand it, admission decides:
          # every step at E :OnlyA leaves F fed at zero product with a dead way back
          # (`_PUMPED`); a tagged isomerization alone joining F to E makes F an island
          # (`_ISLAND`). Each rejection has an admitted variant: the S binding untagged
          # (AD2d), or the catalysis tagged in place of the isomerization (P4-1's base).
          r1 = @allosteric_mechanism begin
              substrates: S
              products: P
              catalytic_multiplicity: 2
              steps: begin
                  E + S ⇌ E(S)      :: EqualAI
                  E(S) <--> E(P)    :: NonequalAI
                  F + P ⇌ E(P)      :: EqualAI
                  F <--> E          :: NonequalAI
              end
          end
          r1_d = @allosteric_mechanism begin
              substrates: S
              products: P
              catalytic_multiplicity: 2
              steps: begin
                  E + S ⇌ E(S)      :: EqualAI
                  E(S) <--> E(P)    :: NonequalAI
                  D + P ⇌ E(P)      :: EqualAI
                  D <--> E          :: NonequalAI
              end
          end
          # R1's oracle is the RE-speed one (no tag excludes a step, P4-29b); the law also
          # equals the exact limit at one subunit, where L means [E_I]/[E_A] for D too.
          for (i, m) in enumerate((r1, r1_d))
              @test all(==(:agree), _testhelper_respeed_verdicts(
                  m; rng = Random.MersenneTwister(20 + i)))
              _testhelper_check_law_against_limit(m; seed = 20 + i)
          end
          r2 = @_testhelper_unguarded @allosteric_mechanism begin
              substrates: S
              products: P
              catalytic_multiplicity: 2
              steps: begin
                  E + S ⇌ E(S)      :: OnlyA
                  E(S) <--> E(P)    :: EqualAI
                  F + P ⇌ E(P)      :: EqualAI
                  F <--> E          :: OnlyA
              end
          end
          v = ER._guard(r2)
          @test (v.fault, v.face) == (ER._PUMPED, 2)
          @test _testhelper_epath_verdict(r2; rng = Random.MersenneTwister(23),
                                          face = Set([:P])) === :direction_dependent
          @test_throws ErrorException(ER._guard_message(r2, v)) @allosteric_mechanism begin
              substrates: S
              products: P
              catalytic_multiplicity: 2
              steps: begin
                  E + S ⇌ E(S)      :: OnlyA
                  E(S) <--> E(P)    :: EqualAI
                  F + P ⇌ E(P)      :: EqualAI
                  F <--> E          :: OnlyA
              end
          end
          ad2d = @allosteric_mechanism begin
              substrates: S
              products: P
              catalytic_multiplicity: 2
              steps: begin
                  E + S ⇌ E(S)      :: EqualAI
                  E(S) <--> E(P)    :: EqualAI
                  F + P ⇌ E(P)      :: EqualAI
                  F <--> E          :: OnlyA
              end
          end
          @test _testhelper_epath_verdict(ER.AllostericMechanism(ad2d);
                                          rng = Random.MersenneTwister(24)) === :converges
          r3 = @_testhelper_unguarded @allosteric_mechanism begin
              substrates: S
              products: P
              modifiers: I
              catalytic_multiplicity: 2
              steps: begin
                  E + S ⇌ E(S)      :: EqualAI
                  E(S) <--> E(P)    :: EqualAI
                  E + P ⇌ E(P)      :: EqualAI
                  E <--> F          :: OnlyA
                  F + I ⇌ F(I)      :: EqualAI
              end
          end
          r4 = @_testhelper_unguarded @allosteric_mechanism begin
              substrates: S
              products: P
              catalytic_multiplicity: 2
              steps: begin
                  E + S ⇌ E(S)      :: EqualAI
                  E(S) <--> E(P)    :: NonequalAI
                  E + P ⇌ E(P)      :: EqualAI
                  F ⇌ E             :: OnlyA
              end
          end
          for (i, m) in enumerate((r3, r4))
              v = ER._guard(m)
              @test v.fault === ER._ISLAND
              @test _testhelper_epath_verdict(m; rng = Random.MersenneTwister(25 + i)) in
                    (:direction_dependent, :wrong_limit)
          end
          msg = ER._guard_message(r3, ER._guard(r3))
          @test_throws ErrorException(msg) @allosteric_mechanism begin
              substrates: S
              products: P
              modifiers: I
              catalytic_multiplicity: 2
              steps: begin
                  E + S ⇌ E(S)      :: EqualAI
                  E(S) <--> E(P)    :: EqualAI
                  E + P ⇌ E(P)      :: EqualAI
                  E <--> F          :: OnlyA
                  F + I ⇌ F(I)      :: EqualAI
              end
          end
          msg = ER._guard_message(r4, ER._guard(r4))
          @test_throws ErrorException(msg) @allosteric_mechanism begin
              substrates: S
              products: P
              catalytic_multiplicity: 2
              steps: begin
                  E + S ⇌ E(S)      :: EqualAI
                  E(S) <--> E(P)    :: NonequalAI
                  E + P ⇌ E(P)      :: EqualAI
                  F ⇌ E             :: OnlyA
              end
          end
          r3_base = @allosteric_mechanism begin
              substrates: S
              products: P
              modifiers: I
              catalytic_multiplicity: 2
              steps: begin
                  E + S ⇌ E(S)      :: EqualAI
                  E(S) <--> E(P)    :: OnlyA
                  E + P ⇌ E(P)      :: EqualAI
                  E <--> F          :: EqualAI
                  F + I ⇌ F(I)      :: EqualAI
              end
          end
          r4_base = @allosteric_mechanism begin
              substrates: S
              products: P
              catalytic_multiplicity: 2
              steps: begin
                  E + S ⇌ E(S)      :: EqualAI
                  E(S) <--> E(P)    :: OnlyA
                  E + P ⇌ E(P)      :: EqualAI
                  F ⇌ E             :: EqualAI
              end
          end
          for (i, m) in enumerate((r3_base, r4_base))
              @test _testhelper_epath_verdict(
                  ER.AllostericMechanism(m); rng = Random.MersenneTwister(27 + i)) ===
                  :converges
          end
          # F in E's rapid-equilibrium segment, untagged, is an ordinary form too (A5).
          @test @allosteric_mechanism(begin
              substrates: S
              products: P
              catalytic_multiplicity: 2
              steps: begin
                  E + S ⇌ E(S)      :: EqualAI
                  E(S) <--> E(P)    :: NonequalAI
                  E + P ⇌ E(P)      :: EqualAI
                  F ⇌ E             :: NonequalAI
              end
          end) isa AllostericEnzymeMechanism
      end
  ```
  The one-flip message pins go with the check (P4-29a); each fixture keeps a verdict test.

  (d) `test/test_identify_rate_equation.jl`: in "_base_tier drops a flip child whose rate law
  the guard rejects", the three flips that cut the ring at F ⇌ E are well-posed now that the
  conformations flip at E (the probe: every face passes). Rename the testset "_base_tier keeps
  every flip child the guard admits", replace its comment's last three sentences with "The six
  flips of two steps of the ring cut it, and each passes every face of the guard: the
  conformations flip at E, so F needs no rapid-equilibrium link to it. No failure is
  recorded.", and add the three flips to `kept`: replace the vector's closing `end)]` with
  `end),` and follow it with these three entries, the last one closing the vector:
  ```julia
          ER.AllostericMechanism(@allosteric_mechanism begin
              substrates: S
              products: P
              steps: begin
                  E + S <--> E(S)   :: EqualAI
                  E(S) ⇌ E(P)       :: EqualAI
                  E(P) ⇌ F + P      :: EqualAI
                  F <--> E          :: EqualAI
              end
          end),
          ER.AllostericMechanism(@allosteric_mechanism begin
              substrates: S
              products: P
              steps: begin
                  E + S ⇌ E(S)      :: EqualAI
                  E(S) <--> E(P)    :: EqualAI
                  E(P) ⇌ F + P      :: EqualAI
                  F <--> E          :: EqualAI
              end
          end),
          ER.AllostericMechanism(@allosteric_mechanism begin
              substrates: S
              products: P
              steps: begin
                  E + S ⇌ E(S)      :: EqualAI
                  E(S) ⇌ E(P)       :: EqualAI
                  E(P) <--> F + P   :: EqualAI
                  F <--> E          :: EqualAI
              end
          end)]
  ```
  and replace `@test count(!ER._is_well_posed, ER._expand_re_to_ss(seed)) == 3` with
  `@test all(ER._is_well_posed, ER._expand_re_to_ss(seed))` and `@test length(base) == 3`
  with `@test length(base) == 6`.

- [ ] **Step 2: Run them and watch them fail.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  for f in test_oracles test_admission test_types test_identify_rate_equation; do $S/phase3/jl.sh $S/phase45/logs/P4_T9_red_$f.log 1800 --project=. -e "using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include(\"test/mechanism_definitions_for_test_enzyme_derivation.jl\"); include(\"test/$f.jl\")"; done
  ```
  Expected: in "the conformations flip at free E" every family member's testset errors, the
  macro throwing the one-flip message ("are not joined by rapid-equilibrium steps outside
  :OnlyA groups"), and the outer testset then errors at its first n = 2 mechanism with the
  same message; "a second free form keeps its law exact ..." errors for R1, AD2a, AD2d,
  AD5, R2P_sealed and B_uphill with the one-flip message (their derivation entry runs it)
  and passes for A5, AD7a, AD7c and the ten AD13 members, whose free forms share E's
  rapid-equilibrium segment; in test_types.jl the S1 testset fails at its `@test_throws`
  (the macro throws the one-flip message) and errors at `f_onlya` (the same message), and
  the free-form testset errors at `r1` (one-flip message); the base-tier testset fails
  (three flips are not well-posed, `length(base) == 3`). Nothing else fails.

- [ ] **Step 3: Flip at free E.**

  (a) `src/types.jl`: replace `_inactive_mask` and its docstring with:
  ```julia
  """
  The edges of `am`'s inactive conformation: those outside `:OnlyA` groups between forms
  they join to free E (`_law_mask`). The conformations interconvert at free E alone, so a
  form those edges do not join to it holds no inactive enzyme, a second free form included.
  The inactive law (`_state_rate_polys`), the joint solve (`_inactive_law`) and the copy rule
  (`_redundant_copy_groups`) read the inactive conformation from it, each edge in `am`'s
  stored direction.
  """
  _inactive_mask(am::AllostericMechanism) = _law_mask(am, 2)
  ```
  and delete `_outside_onlya` with its docstring (Step 0 found no other caller).

  (b) `src/rate_eq_derivation.jl`, in `_raw_symbolic_rate_polys`: replace
  ```julia
      # A graph with no unbound, residual-free form (e.g. a hand-written mechanism whose
      # every form carries a residual) has no free enzyme to weight, so `d_free` is 1.
      i_free = findfirst(f -> seg[f] != 0 && _is_free(cm.forms[f]), eachindex(cm.forms))
      d_free = i_free === nothing ? poly_one() : D[seg[i_free]]
  ```
  with
  ```julia
      # The conformations interconvert at free E, form 1, which every state's edges reach
      # (`_state_mask`). It roots its segment (`_segment_root` ranks it first), so its own
      # alpha is 1 and its weight is its segment's tree weight.
      d_free = D[seg[1]]
  ```
  and in its docstring replace "Also returns `d_free`, the weight in the returned denominator
  of the free resting enzyme (the first form with empty `bound` and empty `residual`):" with
  "Also returns `d_free`, the weight in the returned denominator of free E, form 1, where the
  conformations interconvert:".

  (c) `src/admission.jl`: delete `_assert_one_flip_segment` and `_free_forms_flip_together`
  with their docstrings; in `_check_well_posed` delete the line
  `m isa AllostericMechanism && _assert_one_flip_segment(m)` and, in its docstring, the
  sentence "An allosteric mechanism must also hold every free form in one rapid-equilibrium
  segment of the edges outside `:OnlyA` groups (`_assert_one_flip_segment`)."; replace both
  `_is_well_posed` methods with
  ```julia
  _is_well_posed(m::_AnyMech) = _guard(m).fault === _NO_FAULT
  ```

  (d) `test/mechanism_definitions_for_test_enzyme_derivation.jl`: delete
  `_testhelper_flip_form` with its docstring; give `_testhelper_steady_state_rate` the
  keyword `flip::Int = 1` (signature `_testhelper_steady_state_rate(m, point,
  concs::AbstractDict; Λ = big(10)^30, flip::Int = 1)`), delete the line
  `flip = _testhelper_flip_form(g)` in its body (the push that follows reads the keyword),
  and in its docstring replace "one flip between the two copies of the flip form
  (`_testhelper_flip_form`)" with "one flip between the two copies of form `flip`, free E
  unless a test asks for another".

  (e) Rewrite each comment and docstring Step 0's last grep found that states the one-flip
  rule or prunes the inactive conformation from any free form; keep what stays true. In
  `test/test_types.jl`, testset "a kinetic group ties one metabolite exchange at two
  conformations", the sentences "The allosteric twin joins F to E at rapid equilibrium: its
  conformations interconvert at one free form, so a second free form must be in rapid
  equilibrium with it." become "The allosteric twin joins F to E at rapid equilibrium." In
  `test/test_reaction_graph.jl`, the comment "The inactive conformation keeps the edges
  outside :OnlyA groups between forms that such edges reach from a free form" ends "from
  free E" instead. The same file's `_testhelper_inactive_reach`, the oracle
  `_testhelper_subgraphs_agree` checks `_inactive_mask` against, seeds its reach from every
  free form (`g.free[s.from]`, `g.free[s.to]`): seed it from free E alone (`s.from == 1`,
  `s.to == 1`; the test reader lists free E first), and its docstring's "reach from a free
  form (no metabolite, no residual)" ends "reach from free E" instead. List every rewritten
  comment in the commit body.

- [ ] **Step 4: Run them and watch them pass.** The Step 2 commands, logs
  `P4_T9_green_<file>.log`. Expected: `exit 0` for each; no warning. Then the
  derivation-heavy files that read the inactive law, focused: `test_rate_eq_derivation`,
  `allosteric_ground_truth`, `test_joint_solve`, `test_reaction_graph` (logs
  `P4_T9_green_<file>.log`); expected `exit 0`, no pin moved (the census shows why: every
  enumerated and spec mechanism holds one free form, so its masks and law are unchanged).

- [ ] **Step 5: Retire the ROADMAP entry.** In `ROADMAP.md`, under "## Planned", delete the
  entry "**Take `:OnlyA` steps to their limit in the inactive conformation.**" (its whole
  text): the guard takes every tagged constant to its limit or rejects the tag set, and the
  conformations flip at free E.

- [ ] **Step 6: Census against Task 8's.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  $S/phase3/jl.sh $S/phase45/logs/P4_T9_census.log 1800 --project=$S/phase3/census/env_guard --heap-size-hint=2500M $S/phase3/census/corpus.jl $S/phase45/guard_wt $S/phase3/census/P4_T9 $S/phase3/census/rename_phase3.toml
  $S/phase3/jl.sh $S/phase45/logs/P4_T9_census_cmp.log 600 --project=$S/phase3/census/env --heap-size-hint=2500M $S/phase3/census/compare.jl $S/phase3/census/P4_T8 $S/phase3/census/P4_T9 $S/phase3/census/rename_phase3.toml
  ```
  Expected: `RESULT PASS` with zero differences (verify_flipE: 9,308 enumerated mechanisms
  hold one free form). A difference stops the task: run `superpowers:systematic-debugging`
  on the first one.

- [ ] **Step 7: Run the full suite.** As Task 8 Step 11, log `P4_T9_full.log`. Expected:
  `exit 0`, pristine, Aqua and JET pass, the performance and compile gates hold.

- [ ] **Step 8: Commit.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && git add src/types.jl src/rate_eq_derivation.jl src/admission.jl ROADMAP.md test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_oracles.jl test/test_admission.jl test/test_types.jl test/test_identify_rate_equation.jl test/test_reaction_graph.jl && git commit -F - <<'EOF2'
  Flip the conformations at free E alone

  The inactive conformation keeps the forms its untagged steps join to free E, and the
  law weights E itself; a second free form is an ordinary form. Admission's face condition
  covers what the one-flip check guarded, so the check goes. The iso families are exact at
  one to three subunits; a flip elsewhere gives another law, and a free form in rapid
  equilibrium with E only rescales L, which means [E_I]/[E_A]. Census equals Task 8's.

  Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
  EOF2
  ```

---

### Task 10: `:OnlyI` on modifier-binding groups (RES §4 row 9)

A catalytic group may be `:OnlyI` when its steps take up or give off modifiers alone (R4.4,
P4-16, P4-17): the constructor checks it; the active law drops the group and prunes from
free E (`_active_mask`, the mirror of `_inactive_mask`); the joint solve gives the group
inactive columns tied by group rows and each excluded active copy a class of its own; the
keep-in-law clause keeps an `:EqualAI` constant that only the inactive law reads; the copy
rule builds each conformation from its own mask; `_degenerate` reads the active law (AM A8).
Admission's `:OnlyI` branches (Tasks 4-5) get their tests here.

**Files:**
- Modify: `src/types.jl` (`_VALID_CAT_ALLO_STATES` and `_VALID_REG_ALLO_STATES` deleted, the
  checks reading `_TAG_CODES`; `_check_tags`' message; new `_check_onlyi_groups` and its
  call in the inner constructor; `_active_mask`)
- Modify: `src/admission.jl` (admission's turnover check: `_NO_ACTIVE_TURNOVER`, its line in
  `_admission_interior`, its message, `_ADMISSION_FAULTS`)
- Modify: `src/thermo.jl` (`_active_law`, `_in_active_law`, `_has_inactive_columns`,
  `_class`, `_ordered_columns`, `_kept_members`, `_dependent_param_exprs`; `_priority_key`'s
  docstring)
- Modify: `src/names.jl` (`_constant_sides`)
- Modify: `src/rate_eq_derivation.jl` (`_state_mask`; the `_state_rate_polys` and
  `_mwc_state_polys` docstrings)
- Modify: `src/mechanism_enumeration.jl` (`_redundant_copy_groups`; `_has_vmax`,
  `_chemistry_equilibrates_both_sides`, `_merged_complexes` and `_degenerate`)
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (`_testhelper_graph`'s
  `active`, `_testhelper_consistent_point`, `_testhelper_declared_dimension`,
  `_testhelper_detailed_balance`, new `_testhelper_onlyi_fixtures`)
- Test: `test/test_oracles.jl`, `test/test_admission.jl`, `test/test_types.jl`,
  `test/test_dsl.jl`, `test/test_identify_rate_equation.jl`,
  `test/test_mechanism_enumeration.jl`, `test/test_reaction_graph.jl`

**Interfaces:**
- Consumes: Task 2's `_law_mask`, `_excludes`, `_eps_network`; Task 1's
  oracles (`_testhelper_steady_state_rate`, `_testhelper_epath_verdict`,
  `_testhelper_respeed_verdicts`, `_testhelper_single_subunit` (concrete),
  `_testhelper_law_value` (concrete, unchecked), `_testhelper_excludes`,
  `@_testhelper_unguarded`); Task 3's `_testhelper_face_faults`; Task 9's
  `_testhelper_check_law_against_limit`; Tasks 4-5's admission (`_admission_interior`,
  `_ADMISSION_FAULTS`, `_NO_LIMIT`, `_ARC_CYCLE`, `_ISLAND`, `_EPS_NONUNIQUE`, `_BRIDGE`,
  `_admission_message`, `_face_sets`, `_face_condition`, `_fewest_slow`,
  `_copy_transitions`, `_turnover_name`, `_has_turnover_cycle`, `_first_excluded`) and
  Task 5's `_testhelper_admission_faults`, `_testhelper_admission_agreement`; phase 3's
  `_modifier_mask`, `_gauge_rescaling`, `_productive_twin`, `_has_vmax(m, re, names)`,
  `_chemistry_equilibrates_both_sides(m, re)`, `_merged_complexes(m, subs)`,
  `_degenerate(m)`, `_step_text(m, k)` (names.jl).
- Produces:
  ```julia
  _check_tags(cat_allo_states::Vector{Symbol}, G::Int, catalytic_multiplicity::Int)::Nothing
  _check_onlyi_groups(catalytic::Mechanism, cat_allo_states::Vector{Symbol})::Nothing
  _active_mask(am::AllostericMechanism)::BitVector        # _law_mask(am, 1)
  _active_law(m::Union{Mechanism, AllostericMechanism})::BitVector
  _in_active_law(k::Int, group::Vector{Int}, tags::Vector{Symbol}, law::BitVector,
                 alaw::BitVector)::Bool
  _ordered_columns(m)::NamedTuple{(:columns, :flat, :tags, :law, :alaw)}
  # _Fault gains, after _BRIDGE: _NO_ACTIVE_TURNOVER (admission's turnover check)
  const _ADMISSION_FAULTS::NTuple{7, _Fault}
  _has_vmax(m::Mechanism, re::BitVector, names::Vector{Symbol},
            on::BitVector = trues(length(m.edges)))::Bool
  _chemistry_equilibrates_both_sides(m::Mechanism, re::BitVector,
                                     on::BitVector = trues(length(m.edges)))::Bool
  _merged_complexes(m::Mechanism, subs::Counts,
                    on::BitVector = trues(length(m.edges)))::Vector{Int}
  _degenerate(m::Union{Mechanism, AllostericMechanism})::Bool      # the active law
  _testhelper_onlyi_fixtures()::NamedTuple   # fixtures file: R4.4's 13, built unguarded
  _testhelper_declared_dimension(m)::Int     # extended to :OnlyI classes
  _testhelper_detailed_balance(m, rng)::Vector{BigFloat}   # reads the concrete law
  ```
  `_VALID_CAT_ALLO_STATES` and `_VALID_REG_ALLO_STATES` are deleted: with `:OnlyI` allowed
  both would equal `_TAG_CODES`.

**Gate:** R4.4's 13 fixtures give a law equal to the exact limit of the oracle's own
deletion network to 1e-20 (`_testhelper_check_law_against_limit`, whose pruning is the
oracle's, not `_active_mask`), detailed balance to 1e-30 and a declared dimension equal to the
fitted count, 13 of 13. The constructor rejects R4.4's 3 other shapes (a reactant binding; an
isomerization of an inhibitor complex; a Theorell–Chance step that takes up a modifier and
gives off a product), naming the step, and builds each `:NonequalAI` fix, which the guard
admits. (R4.4 also wrote the reactant binding in one group with its modifier copy's
binding; after phase 3 the kinetic-group rule refuses that group whatever its tag, before
the `:OnlyI` rule runs, so no test reaches the rule with it.) The probe's verdicts hold
(`p3_verdicts.log`, `p3_verdicts_t10.log`, `p3_verdicts_t10b.log`): the other 11 fixtures
admitted; copy_group_bridge (RES §3.2's OnlyI_bridge) and B2b `_BRIDGE` at the interior in
the active conformation, the ε-path not converging, the steady-state fix admitted;
essential_activator `_NO_ACTIVE_TURNOVER`, both fixes admitted; B1a, B1b, B1d, B1e, B1f,
bridge_REcycle and bridge_REcycle2 `_RE_TURNOVER` on the ε-network, class 3 and the RE-speed
oracle agreeing on every face of each, and the ε-path not converging on B1a (bridge_REcycle2's
512-bit run covers its zero-modifier face); B1c admitted;
bridge_REcycle_SSbind class 3 `_LIMIT_SPEED` at zero products (RE-speed `:pattern`), its
other faces agreeing; fz101_1964 `_EPS_NONUNIQUE` at zero products with and without
modifiers (ε-path `:nonunique`), class 3 passing every face; B2a admitted by admission on
every face and rejected by class 3 on its three zero-modifier faces (`_LIMIT_DIFFER`,
RE-speed `:mismatch`); a 512-bit ε-path of bridge_REcycle2 at zero modifiers converges. No
move emits a catalytic `:OnlyI` group and the tagging move keeps `:OnlyI` site tags;
`_degenerate` reads the active law (AM A8), pinned on an admitted activator route that the
whole graph would call degenerate (`p3_degenerate.log`), and the edge mask restricts the
maximal-rate test, which folding the mask into the flags would not; the copy rule flags none
of the 13 (R4.4). P4-8b: a ping-pong with half 2 `:OnlyA` and the M bindings at E(A) and at
the covalent form `:OnlyI` realizes VC X4's pattern (`joint_r1_mech.log`, research name
jointR1_a): every check passes but B v3's joint R1 at zero products and at P = 0, where each
conformation's vanishing transitions alone leave the cheapest in-tree at 4 slow transitions
and both together cut it to 3, and the ε-path does not converge there, so the rejection is
correct; the RE-speed oracle agrees with class 3 on all 14 faces and the ε-path converges on
the other 12. Admission's interior faults in the active conformation, each new to the tests:
an `:OnlyI` binding with no limit and an `:OnlyI` arc cycle (ε-path `:no_path`) and an
`:OnlyI` island (ε-path wrong limit) (`onlyi_active_faults.log`), each message without a
fix. Census `P4_T10` equals `P4_T9`. The new testsets add at most 2 s.

- [ ] **Step 0: Check the names this task reads.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  grep -n "_VALID_CAT_ALLO_STATES\|_VALID_REG_ALLO_STATES\|const _TAG_CODES\|function _check_tags\|_check_tags(" src/*.jl
  grep -n "^const .*RegulatorySite(\|^const .*AllostericMechanism(" src/*.jl
  grep -n "function _degenerate\|function _has_vmax\|function _chemistry_equilibrates_both_sides" -A12 src/*.jl
  grep -n "const _ADMISSION_FAULTS\|^    _BRIDGE$\|^    _DENOMINATOR_VANISHES$\|return _Verdict(_ISLAND" src/admission.jl
  grep -n "function _redundant_copy_groups" -A40 src/mechanism_enumeration.jl | grep -n "_inactive_mask\|twin_active\|trues"
  grep -n "_state_mask\|_inactive_law\|_has_inactive_columns\|_kept_members(m" src/*.jl
  grep -n "in_law = allo\|inactive = Int\[\|key = c == 2" test/mechanism_definitions_for_test_enzyme_derivation.jl
  grep -rn "OnlyI is rejected\|catalytic group\" eval" test/*.jl
  ```
  Expected: `_VALID_CAT_ALLO_STATES` read by `_check_tags` only, `_VALID_REG_ALLO_STATES` by
  the `RegulatorySite` constructor only, `_TAG_CODES` defined once in types.jl; no
  top-level constant builds a `RegulatorySite` or an `AllostericMechanism` (so the checks
  may read `_TAG_CODES`, defined later in the file, at run time); `_check_tags` called by
  the inner constructor and by `AllostericMechanism(raw)`. `_degenerate` reads the
  catalytic graph with every edge kept, through `_has_vmax(m::Mechanism, re::BitVector,
  names::Vector{Symbol})` and `_chemistry_equilibrates_both_sides(m::Mechanism,
  re::BitVector)`, each `re` the rapid-equilibrium flag of every edge; the chemistry test
  reads its merged complexes through `_merged_complexes(m, subs)`, over every edge. Step 4
  (f) gives all three an edge mask. `_redundant_copy_groups` builds its
  active state from `trues(length(cm.edges))` and its inactive state from
  `_inactive_mask`; the fixtures file
  has `_testhelper_graph`'s `in_law = allo ? EnzymeRates._inactive_mask(mech) : ...`, its
  `inactive = Int[...]` field and Task 1's `key = c == 2 && tag(s.group) === :OnlyA ? ...`
  in `_testhelper_consistent_point`; the `:OnlyI` rejection texts in test_dsl.jl (three
  places), the old `_check_tags` message in test_types.jl (`_testhelper_invalid_tag`, which
  "AllostericMechanism validation errors" reads for `:OnlyI` and `:Bogus`) and in
  test_reaction_graph.jl (two pins), and the `ErrorException` tests in test_types.jl (two
  places, found by "Catalytic group :OnlyI" and "R-state-active convention"); in
  `src/admission.jl`,
  `_ADMISSION_FAULTS` as Task 5 widened it, the enum's `_BRIDGE` line followed by
  `_DENOMINATOR_VANISHES` (Task 6), and `_admission_interior`'s island return, which Step
  4 (g) follows with the turnover check.

- [ ] **Step 1: Teach the oracles' helpers the `:OnlyI` conformation.**
  `test/mechanism_definitions_for_test_enzyme_derivation.jl`:

  (a) In `_testhelper_graph`, after the line `in_law = allo ? ...`, add
  ```julia
      # The active law keeps the steps outside :OnlyI groups on forms they join to free E.
      in_active = allo ? EnzymeRates._law_mask(mech, 1) : trues(length(cm.edges))
  ```
  in the returned tuple, after `inactive = Int[j for (j, k) in enumerate(edges) if in_law[k]],`,
  add `active = Int[j for (j, k) in enumerate(edges) if in_active[k]],`, and in its
  docstring replace "`inactive`, which steps the inactive law keeps;" with "`inactive` and
  `active`, which steps the inactive and the active law keep (every step of a plain
  mechanism is active);". The helper reads `_law_mask` (Task 2), not Step 4's
  `_active_mask`, so every other test runs while Step 4 is unwritten.

  (b) In `_testhelper_consistent_point`: replace the tie loop's body with
  ```julia
      for (j, s) in enumerate(g.steps)
          j1 = findfirst(t -> t.group == s.group, g.steps)
          s1 = g.steps[j1]
          # Each conformation's copy of an excluded step is a class of its own: no
          # conformation-1 tie for :OnlyI, no conformation-2 tie for :OnlyA.
          j == j1 || tag(s.group) === :OnlyI ||
              push!(ties, _testhelper_drop(1, s) - _testhelper_drop(1, s1))
          allo || continue
          tag(s.group) in (:NonequalAI, :OnlyI) && j != j1 &&
              push!(ties, _testhelper_drop(2, s) - _testhelper_drop(2, s1))
          tag(s.group) === :EqualAI &&
              push!(ties, _testhelper_drop(2, s) - _testhelper_drop(1, s))
      end
  ```
  in `_testhelper_conformation_steps`, replace `inactive && j ∉ g.inactive && continue`
  with `j ∉ (inactive ? g.inactive : g.active) && continue`; in Task 1's block that fills
  `K` and `kf`, replace
  `s.rapid_equilibrium || isempty(names) || (rates[(s.group, c)] = vals[names[1]])` with
  ```julia
          s.rapid_equilibrium || isempty(names) || !haskey(vals, names[1]) ||
              (rates[(s.group, c)] = vals[names[1]])
  ```
  and `key = c == 2 && tag(s.group) === :OnlyA ? (-k, 2) : (s.group, c)` with
  `key = _testhelper_excludes(mech, k, c) ? (-k, c) : (s.group, c)`; in the docstring,
  replace "(each member of an `:OnlyA` group in the inactive conformation its own)" with
  "(each member of a group in the conformation that excludes it its own: an `:OnlyA`
  group's in the inactive conformation, an `:OnlyI` group's in the active one)" and
  "`active[j]` (and `inactive[j]` for the inactive law's steps)" with "`active[j]` for the
  active law's steps (and `inactive[j]` for the inactive law's)".

  (c) Replace `_testhelper_declared_dimension`'s docstring and its inner
  `_testhelper_classes` with
  ```julia
  """
  The dimension of the model `m` declares, counted without the package's thermodynamic
  solve: the classes of constants in its law (a kinetic group's constant, slot by slot, in
  each conformation whose law keeps a step of the group; an `:EqualAI` group shares one
  class between the two), less the rank the cycle rows of both conformations impose on
  them, plus each regulator constant and `L`. A cycle row also reads the classes outside
  the law (the copy of a step in the conformation that excludes it, an `:OnlyA` step's
  inactive copy or an `:OnlyI` step's active one, each a class of its own, and the classes
  on forms a conformation prunes), which absorb what they can: the law keeps
  `n_law − (rank(rows) − rank(rows over the classes outside the law))` dimensions.
  """
  ```
  ```julia
      function _testhelper_classes(j, inactive)
          s = g.steps[j]
          t = tag(s.group)
          slots = s.rapid_equilibrium ? (1,) : (1, 2)
          t === (inactive ? :OnlyA : :OnlyI) &&
              return [_testhelper_class((:copy, j, k, inactive), false) for k in slots]
          kept = j in (inactive ? g.inactive : g.active)
          key = !inactive || t === :EqualAI ? :A : :I
          [_testhelper_class((key, s.group, k), kept) for k in slots]
      end
  ```
  An own class's key now carries `:copy` and the conformation, so it can no longer equal a
  group class's key, as `(:I, j, k)` could when step `j`'s number was a group's.

  (d) Replace `_testhelper_detailed_balance`'s first lines, from `em = ...` to
  `params[:Keq] = ...`, with
  ```julia
          mech = _testhelper_lift(m)
          dep, fitted = EnzymeRates._dependent_param_exprs(mech)
          params = Dict{Symbol, BigFloat}(k => exp(4 * big(rand(rng)) - 2) for k in fitted)
          params[:Keq] = exp(4 * big(rand(rng)) - 2)
  ```
  delete its later line `dep, _ = EnzymeRates._dependent_param_exprs(typeof(em))`, read
  `EnzymeRates.metabolites(mech)` and `_testhelper_law_value(mech, ...)` where it read
  `em`, and append to its docstring: "`m`, compiled or concrete, is read as the concrete
  mechanism, without the guard (`_testhelper_law_value`), so a law the guard rejects is
  checked too." The draws come in the same order, so every existing pin reads the same
  values.

  (e) At the end of the "Phase-4 oracles" section add:
  ```julia
  """
  R4.4's thirteen mechanisms whose `:OnlyI` groups bind modifiers alone, by name, built
  without the guard (`@_testhelper_unguarded`): a competitive and an uncompetitive
  inhibitor, each with rapid-equilibrium and with steady-state steps; a substrate's and a
  product's modifier copy (the product's twice); a copy group; an inhibitor that bridges
  it; an essential and a nonessential activator; an `:OnlyI` binding beside an `:OnlyA`
  one; and an ordered bi-bi with a regulator. The guard rejects `copy_group_bridge`
  (`_BRIDGE`) and `essential_activator` (`_NO_ACTIVE_TURNOVER`); their laws still exist.
  """
  function _testhelper_onlyi_fixtures()
      competitive_RE = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)            :: EqualAI
              E(S) <--> E + P         :: EqualAI
              E + I ⇌ E(I)            :: OnlyI
          end
      end
      competitive_SS = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)            :: EqualAI
              E(S) <--> E + P         :: NonequalAI
              E + I <--> E(I)         :: OnlyI
          end
      end
      uncompetitive_RE = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)            :: EqualAI
              E(S) <--> E + P         :: EqualAI
              E(S) + I ⇌ E(S, I)      :: OnlyI
          end
      end
      uncompetitive_SS = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S <--> E(S)         :: EqualAI
              E(S) <--> E + P         :: EqualAI
              E(S) + I <--> E(S, I)   :: OnlyI
          end
      end
      mod_copy_substrate = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: S
          steps: begin
              E + S ⇌ E(S)                    :: EqualAI
              E(S) <--> E + P                 :: EqualAI
              E(S) + S::Mod ⇌ E(S, S::Mod)    :: OnlyI
          end
      end
      mod_copy_product_twin = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: P
          steps: begin
              E + S ⇌ E(S)            :: EqualAI
              E(S) <--> E(P)          :: EqualAI
              E + P ⇌ E(P)            :: EqualAI
              E + P::Mod ⇌ E(P::Mod)  :: OnlyI
          end
      end
      mod_copy_product_twin_nonequal = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: P
          steps: begin
              E + S ⇌ E(S)            :: EqualAI
              E(S) <--> E(P)          :: EqualAI
              E + P ⇌ E(P)            :: NonequalAI
              E + P::Mod ⇌ E(P::Mod)  :: OnlyI
          end
      end
      copy_group_mixed = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)                          :: EqualAI
              E(S) <--> E + P                       :: EqualAI
              (E + I ⇌ E(I), E(S) + I ⇌ E(S, I))    :: OnlyI
          end
      end
      copy_group_bridge = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S <--> E(S)                       :: EqualAI
              E(S) <--> E + P                       :: EqualAI
              (E + I ⇌ E(I), E(S) + I ⇌ E(S, I))    :: OnlyI
              E(I) + S ⇌ E(S, I)                    :: EqualAI
          end
      end
      essential_activator = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: X
          steps: begin
              E + X ⇌ E(X)            :: OnlyI
              E(X) + S ⇌ E(S, X)      :: EqualAI
              E(S, X) <--> E(X) + P   :: EqualAI
          end
      end
      nonessential_activator = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: X
          steps: begin
              E + S ⇌ E(S)            :: EqualAI
              E(S) <--> E + P         :: EqualAI
              E + X ⇌ E(X)            :: OnlyI
              E(X) + S ⇌ E(S, X)      :: EqualAI
              E(S, X) <--> E(X) + P   :: NonequalAI
          end
      end
      with_OnlyA = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I, J
          steps: begin
              E + S ⇌ E(S)            :: EqualAI
              E(S) <--> E + P         :: NonequalAI
              E + I ⇌ E(I)            :: OnlyI
              E(S) + J ⇌ E(S, J)      :: OnlyA
          end
      end
      bibi_regulator = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B
          products: P, Q
          modifiers: I
          allosteric_regulators: R
          steps: begin
              E + A ⇌ E(A)            :: EqualAI
              E(A) + B ⇌ E(A, B)      :: EqualAI
              E(A, B) <--> E(P, Q)    :: NonequalAI
              E(Q) + P ⇌ E(P, Q)      :: EqualAI
              E + Q <--> E(Q)         :: EqualAI
              E(A) + I ⇌ E(A, I)      :: OnlyI
          end
          regulatory_site(1): R::NonequalAI
      end
      (; competitive_RE, competitive_SS, uncompetitive_RE, uncompetitive_SS,
       mod_copy_substrate, mod_copy_product_twin, mod_copy_product_twin_nonequal,
       copy_group_mixed, copy_group_bridge, essential_activator, nonessential_activator,
       with_OnlyA, bibi_regulator)
  end
  ```
  A function, not a constant: the fixtures file is included before every test file, and
  until Step 4 the constructor rejects these tags, which would stop every file at its
  include.

- [ ] **Step 2: Write the failing tests.**

  (a) Append to `test/test_oracles.jl`:
  ```julia
  @testset "an :OnlyI law is exact, balanced and counts its declared dimension" begin
      # R4.4's thirteen mechanisms whose :OnlyI groups bind modifiers alone. Each law, read
      # unchecked from the concrete mechanism, equals the exact limit of the oracle's own
      # deletion network (each conformation without the steps it excludes and the forms
      # that strands, pruned by the oracle itself; the conformations flip at free E),
      # satisfies detailed balance, and fits as many constants as the model declares. The
      # guard rejects two of them; the identities hold for their laws as well.
      for (i, (label, am)) in enumerate(pairs(_testhelper_onlyi_fixtures()))
          @testset "$label" begin
              _testhelper_check_law_against_limit(am; seed = 400 + i)
              rng = Random.MersenneTwister(450 + i)
              @test all(<(big(10)^-30), _testhelper_detailed_balance(am, rng))
              @test length(EnzymeRates._dependent_param_exprs(am)[2]) ==
                    _testhelper_declared_dimension(am)
          end
      end
  end

  @testset "a 512-bit ε-path of bridge_REcycle2 converges at zero modifiers" begin
      # RES §3.1's bridge_REcycle2: at zero modifiers the I bindings vanish with I, and
      # every path of the ε-network tends to the law. At 256 bits the paths lose the limit
      # to rounding (RES §3.1), so the check runs at 512. Class 2 rejects the mechanism on
      # the interior ("an RE binding whose affinity vanishes still closes a fast cycle").
      m = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)                      :: EqualAI
              E(S) ⇌ F + P                      :: OnlyA
              F <--> E                          :: EqualAI
              (E + I ⇌ E(I), F + I ⇌ F(I))      :: OnlyI
              F(I) ⇌ E(I)                       :: OnlyA
          end
      end
      @test _testhelper_epath_verdict(m; rng = Random.MersenneTwister(91),
                                      face = Set([:I]), precision = 512) === :converges
  end
  ```

  (b) Append to `test/test_admission.jl`:
  ```julia
  # ── :OnlyI on modifier-binding groups ─────────────────────────────────────────

  @testset "an :OnlyI modifier binding mirrors :OnlyA in the active conformation" begin
      fx = _testhelper_onlyi_fixtures()
      # R4.4's fixtures but the inhibitor bridge and the essential activator pass the
      # guard (`p3_verdicts.log`), and the oracles agree on every face of four of them:
      # steady-state and rapid-equilibrium steps, a substrate's modifier copy and an
      # activator route.
      admitted = [k for k in keys(fx) if !(k in (:copy_group_bridge, :essential_activator))]
      @test length(admitted) == 11
      for k in admitted
          @test (k, EnzymeRates._guard(fx[k]).fault) == (k, EnzymeRates._NO_FAULT)
      end
      for (i, k) in enumerate((:competitive_SS, :uncompetitive_RE, :mod_copy_substrate,
                               :nonessential_activator))
          agree = _testhelper_admission_agreement(
              fx[k]; rng = Random.MersenneTwister(500 + i))
          @test (k, all(agree)) == (k, true)
      end
      # The competitive inhibitor tagged :OnlyI is the :OnlyA one with the conformations
      # swapped: each conformation's law keeps what the other's drops.
      onlya = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)            :: EqualAI
              E(S) <--> E + P         :: EqualAI
              E + I ⇌ E(I)            :: OnlyA
          end
      end
      onlyi = fx.competitive_RE
      @test EnzymeRates._active_mask(onlyi) == EnzymeRates._inactive_mask(onlya)
      @test EnzymeRates._inactive_mask(onlyi) == EnzymeRates._active_mask(onlya)
      @test !all(EnzymeRates._active_mask(onlyi))
  end

  @testset "an RE binding whose affinity vanishes still closes a fast cycle" begin
      # The B1 set (P4-6): an :OnlyI binding at rapid equilibrium loses its affinity in
      # the active conformation, yet on the ε-network it stays fast and closes, with the
      # other rapid-equilibrium steps, a cycle that turns the reaction over: class 2
      # rejects on the interior. The law networks hold no such cycle, so class 3 and the
      # RE-speed oracle agree on every face, while the ε-path does not converge (RES §3.2
      # reads it diverging; checked on B1a).
      b1a = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)            :: EqualAI
              E(S) <--> E(P)          :: EqualAI
              E + P ⇌ E(P)            :: EqualAI
              E + I ⇌ E(I)            :: OnlyI
              E(I) + S ⇌ E(I, S)      :: EqualAI
              E(I, S) ⇌ E(I, P)       :: OnlyA
              E(P) + I ⇌ E(I, P)      :: OnlyI
          end
      end
      b1b = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)            :: OnlyA
              E + P ⇌ E(P)            :: OnlyA
              E(S) <--> E(P)          :: EqualAI
              E(S) + I ⇌ E(I, S)      :: OnlyI
              E(I, S) ⇌ E(I, P)       :: EqualAI
              E(P) + I ⇌ E(I, P)      :: OnlyI
          end
      end
      b1d = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)                          :: EqualAI
              E(S) <--> E(P)                        :: EqualAI
              E + P ⇌ E(P)                          :: EqualAI
              (E + I ⇌ E(I), E(P) + I ⇌ E(I, P))    :: OnlyI
              E(I) + S ⇌ E(I, S)                    :: EqualAI
              E(I, S) ⇌ E(I, P)                     :: OnlyA
          end
      end
      b1e = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: P
          steps: begin
              E + S ⇌ E(S)                    :: EqualAI
              E(S) <--> E(P)                  :: EqualAI
              E + P ⇌ E(P)                    :: EqualAI
              E + P::Mod ⇌ E(P::Mod)          :: OnlyI
              E(P::Mod) + S ⇌ E(P::Mod, S)    :: EqualAI
              E(P::Mod, S) ⇌ E(P::Mod, P)     :: OnlyA
              E(P) + P::Mod ⇌ E(P::Mod, P)    :: OnlyI
          end
      end
      b1f = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B
          products: P, Q
          modifiers: I
          steps: begin
              E + A ⇌ E(A)                    :: EqualAI
              E(A) + B ⇌ E(A, B)              :: EqualAI
              E(A, B) <--> E(P, Q)            :: EqualAI
              E(Q) + P ⇌ E(P, Q)              :: EqualAI
              E + Q ⇌ E(Q)                    :: EqualAI
              E(A, B) + I ⇌ E(A, B, I)        :: OnlyI
              E(A, B, I) ⇌ E(I, P, Q)         :: OnlyA
              E(P, Q) + I ⇌ E(I, P, Q)        :: OnlyI
          end
      end
      bridge_REcycle = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)                      :: EqualAI
              E(S) ⇌ F + P                      :: EqualAI
              F <--> E                          :: EqualAI
              (E + I ⇌ E(I), F + I ⇌ F(I))      :: OnlyI
              F(I) ⇌ E(I)                       :: OnlyA
          end
      end
      bridge_REcycle2 = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)                      :: EqualAI
              E(S) ⇌ F + P                      :: OnlyA
              F <--> E                          :: EqualAI
              (E + I ⇌ E(I), F + I ⇌ F(I))      :: OnlyI
              F(I) ⇌ E(I)                       :: OnlyA
          end
      end
      for (i, (label, m)) in enumerate(pairs((; b1a, b1b, b1d, b1e, b1f, bridge_REcycle,
                                                 bridge_REcycle2)))
          v = EnzymeRates._guard(m)
          @test (label, v.fault, v.face) == (label, EnzymeRates._RE_TURNOVER, 0)
          msg = EnzymeRates._guard_message(m, v)
          @test (label, occursin("(whose affinity vanishes in the ", msg) &&
                        occursin("turns the reaction over", msg)) == (label, true)
          # The law networks agree with the RE-speed oracle face by face: class 3 accepts
          # every face but bridge_REcycle's zero products, where both see the steady state
          # move with the speeds.
          respeed = _testhelper_respeed_verdicts(m; rng = Random.MersenneTwister(60 + i),
                                                 law = :limit)
          @test (label, all(_testhelper_agrees.(_testhelper_face_faults(m), respeed))) ==
                (label, true)
      end
      @test all(==(:agree), _testhelper_respeed_verdicts(
          b1a; rng = Random.MersenneTwister(61), law = :limit))
      @test _testhelper_epath_verdict(b1a; rng = Random.MersenneTwister(62)) in
            (:diverges, :direction_dependent, :wrong_limit)
      # B1c: the inhibitor complex's chemistry at steady state breaks the fast cycle.
      b1c = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)            :: EqualAI
              E(S) <--> E(P)          :: EqualAI
              E + P ⇌ E(P)            :: EqualAI
              E + I ⇌ E(I)            :: OnlyI
              E(I) + S ⇌ E(I, S)      :: EqualAI
              E(I, S) <--> E(I, P)    :: OnlyA
              E(P) + I ⇌ E(I, P)      :: OnlyI
          end
      end)
      @test all(_testhelper_admission_agreement(b1c; rng = Random.MersenneTwister(63)))
      # With bridge_REcycle's I bindings at steady state no fast cycle is left; class 3
      # then rejects zero products, where the steady state depends on how fast the
      # rapid-equilibrium steps run, and passes the other faces.
      ssbind = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)                          :: EqualAI
              E(S) ⇌ F + P                          :: EqualAI
              F <--> E                              :: EqualAI
              (E + I <--> E(I), F + I <--> F(I))    :: OnlyI
              F(I) ⇌ E(I)                           :: OnlyA
          end
      end
      v = EnzymeRates._guard(ssbind)
      @test v.face == 2 && v.limit === EnzymeRates._LIMIT_SPEED
      @test _testhelper_respeed_verdicts(ssbind; rng = Random.MersenneTwister(64),
                                         law = :limit) ==
            [:agree, :pattern, :agree, :agree, :agree, :agree]
  end

  @testset "the ε-network's second closed class rejects the face" begin
      # fz101_1964 (P4-5): at zero products the :OnlyA steps vanish in the inactive
      # conformation and the :OnlyI binding in the active one, and the alive transitions
      # of the ε-network leave two closed classes, so the steady state depends on where
      # the enzyme starts; the same at zero products and modifiers. Class 3 passes every
      # face: the law networks keep one class each.
      fz = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)            :: OnlyA
              E(S) <--> E + P         :: OnlyA
              E(S) + I ⇌ E(I, S)      :: OnlyI
              E(I, S) ⇌ E(I) + P      :: NonequalAI
          end
      end
      v = EnzymeRates._guard(fz)
      @test (v.fault, v.face) == (EnzymeRates._EPS_NONUNIQUE, 2)
      NO, NU = EnzymeRates._NO_FAULT, EnzymeRates._EPS_NONUNIQUE
      @test _testhelper_admission_faults(fz) == [NO, NU, NO, NO, NU, NO]
      @test all(==(NO), _testhelper_face_faults(fz))
      @test _testhelper_epath_verdict(fz; rng = Random.MersenneTwister(65),
                                      face = Set([:P])) === :nonunique
      @test _testhelper_epath_verdict(fz; rng = Random.MersenneTwister(66),
                                      face = Set([:P, :I])) === :nonunique
  end

  @testset "an :OnlyI bridge is rejected" begin
      # copy_group_bridge (R4.4) and B2b: the :OnlyI bindings of I at E and E(S) bridge,
      # in the active conformation, forms its other steps keep apart, so the limit depends
      # on how fast the vanishing affinities fall; the ε-path settles on another value.
      # With the bridging group at steady state the guard admits the mechanism, the fix
      # the message prints. At zero modifiers class 3 also rejects copy_group_bridge.
      bridge = _testhelper_onlyi_fixtures().copy_group_bridge
      b2b = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I, J
          steps: begin
              E + S <--> E(S)                       :: EqualAI
              E(S) <--> E + P                       :: NonequalAI
              (E + I ⇌ E(I), E(S) + I ⇌ E(I, S))    :: OnlyI
              E(I) + S ⇌ E(I, S)                    :: EqualAI
              E + J ⇌ E(J)                          :: OnlyA
          end
      end
      for m in (bridge, b2b)
          v = EnzymeRates._guard(m)
          @test (v.fault, v.face, v.conformation) == (EnzymeRates._BRIDGE, 1, 1)
          @test _testhelper_admission_faults(m)[1:3] == fill(EnzymeRates._BRIDGE, 3)
          @test all(!=(EnzymeRates._NO_FAULT), _testhelper_face_faults(m)[4:6])
      end
      msg = EnzymeRates._guard_message(bridge, EnzymeRates._guard(bridge))
      @test startswith(msg, "AllostericMechanism: the :OnlyI step ") &&
            occursin(" bridges forms of the active conformation that its other steps " *
                     "keep apart, so the limit depends on how fast the tagged constants " *
                     "vanish. Make the group of ", msg) && endswith(msg, " steady-state.")
      # The research ε-path went elsewhere in the RE-first order and converged ε-first,
      # so the oracle's five orders read it wrong-limit or direction-dependent.
      @test _testhelper_epath_verdict(bridge; rng = Random.MersenneTwister(71)) in
            (:wrong_limit, :direction_dependent)
      fixed = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S <--> E(S)                           :: EqualAI
              E(S) <--> E + P                           :: EqualAI
              (E + I <--> E(I), E(S) + I <--> E(S, I))  :: OnlyI
              E(I) + S ⇌ E(S, I)                        :: EqualAI
          end
      end)
      @test all(_testhelper_admission_agreement(fixed; rng = Random.MersenneTwister(72)))
      # B2a bridges in both conformations, yet admission passes every face; class 3
      # rejects the three faces with zero modifiers, where the limit differs from the law.
      b2a = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I, J
          steps: begin
              E + S <--> E(S)                       :: EqualAI
              E(S) <--> E + P                       :: EqualAI
              (E + I ⇌ E(I), E(S) + I ⇌ E(I, S))    :: OnlyA
              E(I) + S ⇌ E(I, S)                    :: EqualAI
              (E + J ⇌ E(J), E(S) + J ⇌ E(J, S))    :: OnlyI
              E(J) + S ⇌ E(J, S)                    :: EqualAI
          end
      end
      v = EnzymeRates._guard(b2a)
      @test v.face == 4 && v.limit === EnzymeRates._LIMIT_DIFFER
      @test _testhelper_admission_faults(b2a)[1:3] == fill(EnzymeRates._NO_FAULT, 3)
      @test _testhelper_respeed_verdicts(b2a; rng = Random.MersenneTwister(73),
                                         law = :limit) ==
            [:agree, :agree, :agree, :mismatch, :mismatch, :mismatch]
  end

  @testset "an essential activator tagged :OnlyI stops catalysis" begin
      # X binds only the inactive conformation, and every turnover passes through E(X):
      # the active conformation, which alone the convention lets catalyze, turns nothing
      # over. Tagging the X binding :EqualAI or :NonequalAI, the message's fix, passes.
      ea = _testhelper_onlyi_fixtures().essential_activator
      v = EnzymeRates._guard(ea)
      @test (v.fault, v.face, v.conformation) == (EnzymeRates._NO_ACTIVE_TURNOVER, 0, 1)
      @test EnzymeRates._guard_message(ea, v) ==
            "AllostericMechanism: without its :OnlyI steps E + X ⇌ E(X) the active " *
            "conformation turns no substrate into product, so the enzyme never turns " *
            "the reaction over. Tag their groups :EqualAI or :NonequalAI."
      for (i, tag) in enumerate((:EqualAI, :NonequalAI))
          fixed = EnzymeRates.AllostericMechanism(eval(:(@allosteric_mechanism begin
              substrates: S
              products: P
              modifiers: X
              steps: begin
                  E + X ⇌ E(X)            :: $tag
                  E(X) + S ⇌ E(S, X)      :: EqualAI
                  E(S, X) <--> E(X) + P   :: EqualAI
              end
          end)))
          @test all(_testhelper_admission_agreement(fixed;
                                                    rng = Random.MersenneTwister(80 + i)))
      end
  end

  @testset "the active conformation's interior faults name the :OnlyI step" begin
      # The mirror of Task 4's :OnlyA faults, in conformation 1. The messages offer no fix:
      # a step tagged :OnlyI binds modifiers alone, and no test shows an untagged variant
      # well-posed. The research rules and both oracles gave each verdict
      # (`onlyi_active_faults.log`).
      rng = Random.MersenneTwister(91)
      U = EnzymeRates._LIMIT_UNTESTED
      # I binds E and E(S); without E + I ⇌ E(I), the untagged I binding at E(S) and S
      # binding at E(I) still join E to E(I), so the I affinity cannot vanish.
      no_limit = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)            :: EqualAI
              E(S) <--> E + P         :: EqualAI
              E + I ⇌ E(I)            :: OnlyI
              E(I) + S ⇌ E(I, S)      :: EqualAI
              E(S) + I ⇌ E(I, S)      :: EqualAI
          end
      end
      k = _testhelper_edge(no_limit, "E + I ⇌ EI")
      v = EnzymeRates._guard(no_limit)
      @test v == EnzymeRates._Verdict(EnzymeRates._NO_LIMIT, U, 0, 1, k)
      @test EnzymeRates._guard_message(no_limit, v) ==
            "AllostericMechanism: in the active conformation, steps outside :OnlyI " *
            "groups still join the two sides of the :OnlyI binding " *
            "$(EnzymeRates._step_text(no_limit.catalytic, k)), so thermodynamics keeps " *
            "its affinity finite and the binding has no limit there."
      @test _testhelper_epath_verdict(no_limit; rng) === :no_path
      # I binds E and J binds E(I), both :OnlyI; the untagged J binding at E and I binding
      # at E(J) join E(I, J) back to E, so the two arcs close a cycle.
      arc_cycle = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I, J
          steps: begin
              E + S ⇌ E(S)            :: EqualAI
              E(S) <--> E + P         :: EqualAI
              E + I ⇌ E(I)            :: OnlyI
              E(I) + J ⇌ E(I, J)      :: OnlyI
              E + J ⇌ E(J)            :: EqualAI
              E(J) + I ⇌ E(I, J)      :: EqualAI
          end
      end
      arcs = [_testhelper_edge(arc_cycle, t) for t in ("E + I ⇌ EI", "EI + J ⇌ EIJ")]
      v = EnzymeRates._guard(arc_cycle)
      @test (v.fault, v.face, v.conformation) == (EnzymeRates._ARC_CYCLE, 0, 1)
      @test v.edge in arcs
      msg = EnzymeRates._guard_message(arc_cycle, v)
      @test startswith(msg, "AllostericMechanism: the :OnlyI bindings ")
      @test all(j -> occursin(EnzymeRates._step_text(arc_cycle.catalytic, j), msg), arcs)
      @test endswith(msg, " close a cycle in the active conformation: steps outside " *
                          ":OnlyI groups join each one's bound side to the next one's " *
                          "unbound side, so thermodynamics forbids their affinities to " *
                          "vanish together.")
      @test _testhelper_epath_verdict(arc_cycle; rng) === :no_path
      # E(I) is reached only as E(I, J) gives off J, and in the active conformation J's
      # affinity for E(I) vanishes: an island on the :OnlyI binding's unbound side, which
      # collects the active enzyme.
      island = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I, J
          steps: begin
              E + S ⇌ E(S)            :: EqualAI
              E(S) <--> E + P         :: EqualAI
              E + J ⇌ E(J)            :: EqualAI
              E(J) + I ⇌ E(I, J)      :: EqualAI
              E(I) + J ⇌ E(I, J)      :: OnlyI
          end
      end
      k = _testhelper_edge(island, "EI + J ⇌ EIJ")
      v = EnzymeRates._guard(island)
      @test v == EnzymeRates._Verdict(EnzymeRates._ISLAND, U, 0, 1, k)
      @test EnzymeRates._guard_message(island, v) ==
            "AllostericMechanism: in the active conformation the forms E(I) join free E " *
            "only through :OnlyI steps that are not bindings toward them, such as " *
            "$(EnzymeRates._step_text(island.catalytic, k)), so their weight depends on " *
            "how fast those steps vanish."
      @test _testhelper_epath_verdict(island; rng) in (:direction_dependent, :wrong_limit)
  end

  @testset "joint R1 rejects a mechanism whose bridges cut the tree only together" begin
      # VC X4's pattern on a mechanism (P4-8b). Half 2 :OnlyA strands the inactive
      # conformation's covalent forms at zero P, so no in-tree of that conformation reaches
      # its free E without half 2's vanishing steps. M binds E(A) and the covalent form
      # :OnlyI, and M-bound half 1 runs at rapid equilibrium in the active conformation,
      # where the vanishing M binding at E(A) lets free E's fast class leave, one way, for
      # the M-bound covalent form. Each conformation's vanishing transitions alone leave the
      # cheapest spanning in-tree at 4 slow transitions; together they cut it to 3. Every
      # other check passes, and the ε-path does not converge on those faces, so the
      # rejection is right (research rules and both oracles, `joint_r1_mech.log`).
      am = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B
          products: P, Q
          modifiers: M
          steps: begin
              E + A ⇌ E(A)                                             :: EqualAI
              E(A) <--> E(; residual = A - P) + P                      :: EqualAI
              E(; residual = A - P) + B ⇌ E(B; residual = A - P)       :: EqualAI
              E(B; residual = A - P) <--> E + Q                        :: OnlyA
              E(A) + M ⇌ E(A, M)                                       :: OnlyI
              E(; residual = A - P) + M <--> E(M; residual = A - P)    :: OnlyI
              E(A, M) ⇌ E(M; residual = A - P) + P                     :: OnlyA
          end
      end
      NO, BR, U = EnzymeRates._NO_FAULT, EnzymeRates._BRIDGE, EnzymeRates._LIMIT_UNTESTED
      v = EnzymeRates._guard(am)
      @test v == EnzymeRates._Verdict(BR, U, 2, 0, 0)
      @test _testhelper_admission_faults(am) == [[NO, BR, NO, NO, NO, BR, NO]; fill(NO, 7)]
      eps, law = EnzymeRates._eps_network(am), EnzymeRates._law_network(am)
      faces = EnzymeRates._faces(am.catalytic)
      for i in (2, 6)   # zero P and Q, zero P
          s = EnzymeRates._face_sets(eps, law.present, faces[i])
          nodes = s.fed .& .!s.stranded
          @test EnzymeRates._fewest_slow(eps, nodes, s.n0) == 4
          for c in 1:2
              (; van) = EnzymeRates._copy_transitions(eps, s.alive, nodes, c)
              @test EnzymeRates._fewest_slow(eps, nodes, s.n0 .| van) == 4
          end
          @test EnzymeRates._fewest_slow(eps, nodes, s.alive) == 3
          @test EnzymeRates._face_condition(eps, law.present, faces[i],
                                            EnzymeRates._turnover_name(am.catalytic)) ==
                EnzymeRates._Verdict(BR, U, 0, 0, 0)
      end
      @test EnzymeRates._guard_message(am, v) ==
            "AllostericMechanism: at zero P and Q, the tagged steps of both conformations " *
            "together bridge forms that each conformation's other steps keep apart, so " *
            "the limit depends on how fast those steps vanish."
      @test all(_testhelper_admission_agreement(am; rng = Random.MersenneTwister(95)))
  end
  ```
  In Task 5's testset "joint R1 catches a tree that only both conformations' bridges
  shorten", replace the comment line `# that shows the pattern, for the in-tree arithmetic.`
  with
  ```julia
      # that shows the pattern, for the in-tree arithmetic; "joint R1 rejects a mechanism
      # whose bridges cut the tree only together" tests the check on a mechanism.
  ```

  (c) `test/test_types.jl`. In testset "AllostericEnzymeMechanism lift + DSL", replace the
  comment `# Catalytic group :OnlyI → error` and `@test_throws ErrorException
  @allosteric_mechanism begin` (the block unchanged) with
  ```julia
          # A catalytic group may be :OnlyI only when its steps bind modifiers alone.
          @test_throws ("AllostericMechanism: E + S ⇌ E(S) is tagged :OnlyI but does more " *
                        "than bind a modifier.") @allosteric_mechanism begin
  ```
  and add after that block
  ```julia
          onlyi = @allosteric_mechanism begin
              substrates: S
              products:   P
              modifiers: I
              steps: begin
                  E + S ⇌ E(S)      :: EqualAI
                  E(S) <--> E + P   :: EqualAI
                  E + I ⇌ E(I)      :: OnlyI
              end
          end
          @test :OnlyI in ER.cat_allo_states(ER.AllostericMechanism(onlyi))
  ```
  In testset "AllostericMechanism validation errors" (its fixture is raw uni-uni steps,
  `_testhelper_tagged_uniuni`), replace the comment
  `# :OnlyI for catalytic group is rejected (R-state-active convention)` and its
  `@test_throws` with an assertion on the composition constructor over a graph:
  ```julia
          # :OnlyI fits only a group that binds modifiers alone: the active conformation
          # is the one that catalyzes.
          uni = ER.Mechanism(@enzyme_mechanism begin
              substrates: S
              products: P
              steps: begin
                  E + S ⇌ E(S)
                  E(S) <--> E(P)
                  E(P) ⇌ E + P
              end
          end)
          onlyi_tags = fill(:EqualAI, ER.ngroups(uni))
          onlyi_tags[uni.group[_testhelper_edge(uni, "E + S ⇌ ES")]] = :OnlyI
          @test_throws ErrorException(
              "AllostericMechanism: E + S ⇌ E(S) is tagged :OnlyI but does more than " *
              "bind a modifier. Only the active conformation catalyzes, so :OnlyI fits " *
              "only a group whose steps take up or give off modifiers alone. Tag its " *
              "group :NonequalAI.") ER.AllostericMechanism(uni, onlyi_tags, 1,
                                                            ER.RegulatorySite[])
  ```
  The fix it prints passes: test_dsl.jl's rule testset below builds each `:NonequalAI`
  variant through the macro, which runs the guard. The testset's helper
  `_testhelper_invalid_tag(tag)`, which its `:Bogus` check still reads, pins `_check_tags`'
  old message: replace its text with `"AllostericMechanism: catalytic group 1 has " *
  "invalid allo state $tag (must be one of (:OnlyA, :OnlyI, :EqualAI, :NonequalAI))"`. In
  `test/test_reaction_graph.jl` two tests pin the old message for an `:OnlyI` catalytic
  group, which now meets `_check_onlyi_groups`: in "an allosteric graph sorts its sites and
  checks its tags and ligands" (group 2 of `_testhelper_graph_uniuni()` is `E + S ⇌ E(S)`)
  and in "the constructors reject what no graph represents, each with its message" (written
  label 2 is `E(S) <--> E(P)`), replace the pinned text with
  `"AllostericMechanism: <step> is tagged :OnlyI but does more than bind a modifier. " *
  "Only the active conformation catalyzes, so :OnlyI fits only a group whose steps " *
  "take up or give off modifiers alone. Tag its group :NonequalAI."`, `<step>` being
  `E + S ⇌ E(S)` and `E(S) <--> E(P)`.

  (d) `test/test_dsl.jl`. In testset "@allosteric_mechanism (parsing & validation)",
  replace the comment `# Reject :OnlyI on a catalytic step (V-type allostery not supported)`
  and `@test_throws ":OnlyI is rejected for catalytic groups" eval(` with
  ```julia
          # Reject :OnlyI on a catalytic step: only the active conformation catalyzes.
          @test_throws ("AllostericMechanism: E(F6P) <--> E(F16BP) is tagged :OnlyI but " *
                        "does more than bind a modifier.") eval(
  ```
  In "@allosteric_mechanism checks names and sites; constructors check states", replace
  `for state in (:Bogus, :OnlyI)` with `for state in (:Bogus,)` and the comment above it
  with `# An unknown catalytic state is reported by AllostericMechanism.`, and add after
  the loop
  ```julia
          # So is :OnlyI on a step that does more than bind a modifier.
          @test_throws ("AllostericMechanism: E + P ⇌ E(P) is tagged :OnlyI but does " *
                        "more than bind a modifier.") eval(:(@allosteric_mechanism begin
              substrates: S
              products:   P
              steps: begin
                  E + S ⇌ E(S)      :: EqualAI
                  E(S) <--> E(P)    :: EqualAI
                  E(P) ⇌ E + P      :: OnlyI
              end
          end))
  ```
  In "constructor errors throw when the expanded code runs", replace
  `@test_throws ":OnlyI is rejected for catalytic groups" @allosteric_mechanism begin`
  with
  ```julia
          @test_throws ("AllostericMechanism: E(S) <--> E(P) is tagged :OnlyI but does " *
                        "more than bind a modifier.") @allosteric_mechanism begin
  ```
  Then append to the same file's top-level `@testset "DSL"`, before its closing `end`:
  ```julia
      @testset "the :OnlyI rule rejects every step that does more than bind a modifier" begin
          # R4.4's three other shapes: a reactant binding, an isomerization of an inhibitor
          # complex, and a Theorell–Chance step that takes up a modifier and gives off a
          # product. Each fix, the group tagged :NonequalAI, passes the guard
          # (`p3_verdicts_t10b.log`).
          lead = "is tagged :OnlyI but does more than bind a modifier."
          for (text, line) in (
                  ("E + S ⇌ E(S) :: OnlyI\nE(S) <--> E + P :: EqualAI\n" *
                   "E + I ⇌ E(I) :: EqualAI", "E + S ⇌ E(S)"),
                  ("E + S ⇌ E(S) :: EqualAI\nE(S) <--> E + P :: EqualAI\n" *
                   "E + I ⇌ E(I) :: EqualAI\nE(I) <--> F(I) :: OnlyI", "E(I) <--> F(I)"),
                  ("E + S ⇌ E(S) :: EqualAI\nE(S) <--> E + P :: EqualAI\n" *
                   "E + I ⇌ E(I) :: EqualAI\nE(S) + I <--> E(I) + P :: OnlyI",
                   "E(S) + I <--> E(I) + P"))
              _testhelper_block(steps) = Meta.parse(
                  "@allosteric_mechanism begin\nsubstrates: S\nproducts: P\n" *
                  "modifiers: I\nsteps: begin\n$steps\nend\nend")
              @test_throws "AllostericMechanism: $line $lead" eval(_testhelper_block(text))
              @test eval(_testhelper_block(replace(text, ":: OnlyI" => ":: NonequalAI"))) isa
                    EnzymeRates.AllostericEnzymeMechanism
          end
      end
  ```
  Each pinned line is `_step_text`'s text (phase 3) for the edge the rule rejects, in its
  stored direction: a binding toward its complex. Where the run of Step 5 shows the same
  edge written the other way round (an isomerization's or a Theorell–Chance step's stored
  direction comes from the orientation tiers), pin the message's text after checking it
  equals `_step_text` on that edge, and say so in the commit body; a message that names
  another edge is a defect: run `superpowers:systematic-debugging`.

  (e) `test/test_identify_rate_equation.jl`, after the `_base_tier` testsets, which read
  `_degenerate` (the file has no testset of `_degenerate` alone; the last is Task 9's
  "_base_tier keeps every flip child the guard admits"):
  ```julia
  @testset "_degenerate reads the active conformation's law" begin
      # An activator X that binds only the inactive conformation opens a route whose
      # chemistry only the inactive law holds. The active law, three steady-state steps,
      # has a maximal rate both ways (AM A8); read over every step, the route's product
      # release turns rapid equilibrium and closes a turnover cycle with the X complex's
      # substrate binding, so the whole graph has no reverse maximal rate
      # (`p3_degenerate.log`). The guard admits the mechanism (`p3_verdicts_t10.log`).
      route = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: X
          steps: begin
              E + S <--> E(S)         :: EqualAI
              E(S) <--> E(P)          :: EqualAI
              E(P) <--> E + P         :: EqualAI
              E + X ⇌ E(X)            :: OnlyI
              E(X) + S ⇌ E(S, X)      :: EqualAI
              E(S, X) <--> E(X) + P   :: EqualAI
          end
      end)
      @test !EnzymeRates._degenerate(route)
      @test EnzymeRates._degenerate(EnzymeRates.Mechanism(@enzyme_mechanism begin
          substrates: S
          products: P
          modifiers: X
          steps: begin
              E + S <--> E(S)
              E(S) <--> E(P)
              E(P) <--> E + P
              E + X ⇌ E(X)
              E(X) + S ⇌ E(S, X)
              E(S, X) <--> E(X) + P
          end
      end))
      # The maximal-rate test turns rapid equilibrium every step that takes up or gives
      # off a substrate, so masking only the flags would bring the excluded steps back:
      # here the X route's two Theorell–Chance steps, each taking up S and giving off P,
      # would close a turnover cycle. Read through the edge mask, the active law keeps
      # only its three steady-state steps, which have a maximal rate both ways.
      tc_route = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: X
          steps: begin
              E + S <--> E(S)           :: EqualAI
              E(S) <--> E(P)            :: EqualAI
              E(P) <--> E + P           :: EqualAI
              E + X ⇌ E(X)              :: OnlyI
              E(X) + S <--> F(X) + P    :: EqualAI
              F(X) + S <--> E(X) + P    :: EqualAI
          end
      end
      cm, on = tc_route.catalytic, EnzymeRates._active_mask(tc_route)
      re = EnzymeRates._flags(cm)
      @test !EnzymeRates._has_vmax(cm, re .& on, cm.substrates)
      @test EnzymeRates._has_vmax(cm, re, cm.substrates, on)
      @test !EnzymeRates._degenerate(tc_route)
  end
  ```

  (f) Append to `test/test_mechanism_enumeration.jl`:
  ```julia
  @testset "no move tags a catalytic group :OnlyI" begin
      # Aggregate check over the whole uni-uni seed set with an allosteric regulator and a
      # competitive inhibitor, and its first expansion round: no seed or child tags a
      # catalytic group :OnlyI (only a hand-written mechanism may), while the tagging move
      # still gives regulatory sites the tag.
      rxn = @enzyme_reaction begin
          substrates: S[C]
          products: P[C]
          allosteric_regulators: R
          competitive_inhibitors: I
          oligomeric_state: 2
      end
      seeds = EnzymeRates.init_mechanisms(rxn)
      children = EnzymeRates.expand_mechanisms(seeds, rxn)
      allo = [m for m in [seeds; children] if m isa EnzymeRates.AllostericMechanism]
      @test !isempty(allo)
      @test all(m -> :OnlyI ∉ EnzymeRates.cat_allo_states(m), allo)
      @test any(m -> any(s -> :OnlyI in EnzymeRates.allo_states(s),
                         EnzymeRates.regulatory_sites(m)), allo)
  end
  ```
  and append to `test/test_admission.jl`, after the `:OnlyI` testsets of (b) (it reads
  `_testhelper_onlyi_fixtures`, which the enumeration file's rule 1 forbids there: a
  fixture pulled from a helper cannot be reviewed inside the testset):
  ```julia
  @testset "_redundant_copy_groups builds each conformation from its own law" begin
      # R4.4's thirteen: the copy rule reads the active conformation without its :OnlyI
      # groups and the inactive one without its :OnlyA groups; it flags none of them,
      # as R4.4's prototype found.
      for (label, am) in pairs(_testhelper_onlyi_fixtures())
          @test (label, EnzymeRates._redundant_copy_groups(am)) == (label, Int[])
      end
  end
  ```
  The copy-rule testset pins a regression, not a difference: R4.4 found no fixture whose
  verdict the per-conformation masks change, and none is known (Open questions, 6).

- [ ] **Step 3: Run them and watch them fail.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  for f in test_oracles test_admission test_types test_dsl test_identify_rate_equation test_reaction_graph; do $S/phase3/jl.sh $S/phase45/logs/P4_T10_red_$f.log 1800 --project=. -e "using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include(\"test/mechanism_definitions_for_test_enzyme_derivation.jl\"); include(\"test/$f.jl\")"; done
  ```
  then `test_mechanism_enumeration` the same way in the background, log
  `P4_T10_red_test_mechanism_enumeration.log`. Expected, and nothing else:
  - test_oracles.jl and test_admission.jl: every new testset errors at its first
    `:OnlyI` fixture with phase 3's message ("catalytic group <g> has invalid allo state
    OnlyI ... :OnlyI is rejected for catalytic groups").
  - test_types.jl, test_dsl.jl and test_reaction_graph.jl: each converted `@test_throws`
    and each `@test_throws` of the new rule testset fails (phase 3's message differs); the
    `onlyi` acceptance errors with phase 3's message; the rule testset's `:NonequalAI` fixes
    pass.
  - test_identify_rate_equation.jl: "_degenerate reads the active conformation's law"
    errors at `route`.
  - test_mechanism_enumeration.jl: "no move tags a catalytic group :OnlyI" passes (no move
    tags one today; it guards the move code against starting). test_admission.jl's
    "_redundant_copy_groups builds each conformation from its own law" errors with the
    other `:OnlyI` testsets.
  - Every other testset passes: Step 1's helper changes are no-ops without `:OnlyI`, so
    the law, balance and dimension pins of the other files read the same values.

- [ ] **Step 4: Admit `:OnlyI` on modifier-binding groups.**

  (a) `src/types.jl`. Replace the comment and the two constants
  ```julia
  # Each catalytic kinetic group carries an allosteric-state tag (`:OnlyA`, `:EqualAI` or
  # `:NonequalAI`; `:OnlyI` is rejected by the active-state convention); a site ligand
  # may also be `:OnlyI`.
  const _VALID_CAT_ALLO_STATES = (:OnlyA, :EqualAI, :NonequalAI)
  const _VALID_REG_ALLO_STATES = (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)
  ```
  with nothing: both checks now accept the four tags, which `_TAG_CODES` already lists.
  In the `RegulatorySite` constructor read `_TAG_CODES` where it read
  `_VALID_REG_ALLO_STATES` (twice; its message keeps the same text, since the tuple
  prints the same). Replace `_TAG_CODES`' docstring with "The allosteric tags, coded by
  position: catalytic groups and site ligands share them, and the constructors accept no
  other (`_check_tags`, `RegulatorySite`)." Replace `_check_tags`' docstring and its third
  check with
  ```julia
  """Error unless `cat_allo_states` holds one valid tag (`_TAG_CODES`) for each of `G`
  kinetic groups and the catalytic multiplicity is at least 1."""
  ```
  ```julia
      bad = findfirst(∉(_TAG_CODES), cat_allo_states)
      bad === nothing ||
          error("AllostericMechanism: catalytic group $bad has invalid " *
                "allo state $(cat_allo_states[bad]) (must be one of $_TAG_CODES)")
  ```
  After `_check_tags` add
  ```julia
  """Error when a kinetic group of `catalytic` tagged `:OnlyI` (`cat_allo_states`) holds a
  step that does more than take up or give off modifiers (`_modifier_mask`): only the
  active conformation catalyzes, so a group it lacks binds modifiers alone (spec
  section 7)."""
  function _check_onlyi_groups(catalytic::Mechanism, cat_allo_states::Vector{Symbol})
      modifiers = _modifier_mask(catalytic)
      for (k, e) in enumerate(catalytic.edges)
          cat_allo_states[catalytic.group[k]] === :OnlyI || continue
          exchanged = e.consumed | e.released
          exchanged != 0 && exchanged & ~modifiers == 0 && continue
          error("AllostericMechanism: $(_step_text(catalytic, k)) is tagged :OnlyI but " *
                "does more than bind a modifier. Only the active conformation " *
                "catalyzes, so :OnlyI fits only a group whose steps take up or give off " *
                "modifiers alone. Tag its group :NonequalAI.")
      end
      nothing
  end
  ```
  and in the inner constructor of `AllostericMechanism`, after the line
  `_check_tags(cat_allo_states, ngroups(catalytic), catalytic_multiplicity)`, add
  `_check_onlyi_groups(catalytic, cat_allo_states)`. After `_inactive_mask` add
  ```julia
  """
  The edges of `am`'s active conformation: those outside `:OnlyI` groups between forms they
  join to free E (`_law_mask`), the mirror of `_inactive_mask`. The active law
  (`_state_rate_polys`), the joint solve (`_active_law`), the copy rule
  (`_redundant_copy_groups`) and the degeneracy test (`_degenerate`) read the active
  conformation from it, each edge in `am`'s stored direction.
  """
  _active_mask(am::AllostericMechanism) = _law_mask(am, 1)
  ```

  (b) `src/thermo.jl`. Replace `_inactive_law`'s two methods and docstring, and
  `_has_inactive_columns` and `_class` with their docstrings, by
  ```julia
  """Whether each edge lies in the inactive conformation's law (`_inactive_mask`); false
  for every edge of a plain mechanism."""
  _inactive_law(m::Mechanism) = falses(length(m.edges))
  _inactive_law(am::AllostericMechanism) = _inactive_mask(am)

  """Whether each edge lies in the active conformation's law (`_active_mask`); true for
  every edge of a plain mechanism."""
  _active_law(m::Mechanism) = trues(length(m.edges))
  _active_law(am::AllostericMechanism) = _active_mask(am)

  """Whether the active column of edge `k`, in a group `group[k]` tagged `tags`, enters the
  law: its edge lies in the active conformation's law (`alaw`), or it is an `:EqualAI`
  edge of the inactive conformation's law (`law`), whose one constant the inactive law
  reads."""
  _in_active_law(k::Int, group::Vector{Int}, tags::Vector{Symbol}, law::BitVector,
                 alaw::BitVector) = alaw[k] || (tags[group[k]] === :EqualAI && law[k])

  """Whether the steps of a group with allosteric tag `tag` take their own constants in the
  inactive conformation: those of a `:NonequalAI`, `:OnlyA` or `:OnlyI` group do; an
  `:EqualAI` group's read the active ones."""
  _has_inactive_columns(tag::Symbol) = tag in (:NonequalAI, :OnlyA, :OnlyI)

  """
  The kinetic class of column `c`, whose edge lies in group `g` with tag `tag`: the group's
  constant in the column's slot and conformation, except that the copy of an edge in the
  conformation that excludes it (an `:OnlyA` edge's inactive copy, an `:OnlyI` edge's
  active one) is a class of its own.
  """
  _class(c::_Column, g::Int, tag::Symbol) =
      c.inactive && tag === :OnlyA || !c.inactive && tag === :OnlyI ?
          (c.edge, c.slot, c.inactive, true) : (g, c.slot, c.inactive, false)
  ```
  In `_priority_key`'s docstring replace "(2: an inactive copy outside the law, 1: another
  inactive copy, 0: active)" with "(2: a copy outside its conformation's law, 1: an
  inactive copy in the law, 0: an active copy in the law)". In `_ordered_columns`: in the
  docstring's signature line write `-> (; columns, flat, tags, law, alaw)`, replace "of a
  `:NonequalAI` or `:OnlyA` group" with "of a `:NonequalAI`, `:OnlyA` or `:OnlyI` group"
  and "and `law[k]` whether edge `k` lies in the inactive conformation's law
  (`_inactive_law`)" with "`law[k]` whether edge `k` lies in the inactive conformation's
  law (`_inactive_law`) and `alaw[k]` whether its active column enters the law (in the
  active conformation's law, `_active_law`, or read by the inactive one,
  `_in_active_law`)"; in the body, after `law = _inactive_law(m)` add
  `alaw = _active_law(m)`, replace
  `block(c) = !c.inactive ? 0 : law[c.edge] ? 1 : 2` with
  ```julia
      block(c) = c.inactive ? (law[c.edge] ? 1 : 2) :
                 _in_active_law(c.edge, cm.group, tags, law, alaw) ? 0 : 2
  ```
  and the last line with
  `(; columns = columns[sortperm(priorities; rev = true)], flat, tags, law, alaw)`.
  Replace `_kept_members`' docstring sentence "A `:NonequalAI` group has an inactive
  class, kept by a member in the inactive law when it has one; an `:EqualAI` group's edges
  share their active class, and each `:OnlyA` edge's inactive copy is a class of its own,
  so neither appears." with "A `:NonequalAI` or `:OnlyI` group has an inactive class, kept
  by a member in the inactive law when it has one; an `:EqualAI` group's edges share their
  active class; the copy of an edge in the conformation that excludes it (an `:OnlyA`
  edge's inactive copy, an `:OnlyI` edge's active one) is a class of its own, so neither
  of these appears.", and its loop body with
  ```julia
          g = group[c.edge]
          tag = cols.tags[g]
          (c.inactive ? !(tag in (:NonequalAI, :OnlyI)) : tag === :OnlyI) && continue
          kept[(g, c.slot, c.inactive)] = c.edge
  ```
  In `_dependent_param_exprs` replace
  ```julia
      kept_in_law(i) = (c = sys.columns[i];
                        (!c.inactive || sys.law[c.edge]) && last_of[class(c)] == i)
  ```
  with
  ```julia
      kept_in_law(i) = (c = sys.columns[i];
                        (c.inactive ? sys.law[c.edge] :
                         _in_active_law(c.edge, group, sys.tags, sys.law, cols.alaw)) &&
                        last_of[class(c)] == i)
  ```
  and in its docstring "(an `:OnlyA` edge's inactive copy, an inactive copy on a pruned
  form)" with "(an `:OnlyA` edge's inactive copy, an `:OnlyI` edge's active copy, a copy
  on a form its conformation prunes)".

  (c) `src/names.jl`, `_constant_sides`: replace
  ```julia
      state === :I && tag === :OnlyA && return _sides(cm, k)
      inactive = state === :I && tag === :NonequalAI
  ```
  with
  ```julia
      (state === :I && tag === :OnlyA || state === :A && tag === :OnlyI) &&
          return _sides(cm, k)
      inactive = state === :I && tag in (:NonequalAI, :OnlyI)
  ```
  and in its docstring "the inactive copy of an `:OnlyA` edge is a class of its own, named
  after its own edge" with "the copy of an edge in the conformation that excludes it (an
  `:OnlyA` edge's inactive copy, an `:OnlyI` edge's active one) is a class of its own,
  named after its own edge".

  (d) `src/rate_eq_derivation.jl`: replace `_state_mask`'s docstring and methods with
  ```julia
  """The edges of `m`'s law in conformation `state`: every edge for `:None`; for `:A` the
  active conformation's (`_active_mask`), for `:I` the inactive one's
  (`_inactive_mask`)."""
  _state_mask(m::Mechanism, ::Symbol) = trues(length(m.edges))
  _state_mask(am::AllostericMechanism, state::Symbol) =
      state === :I ? _inactive_mask(am) : _active_mask(am)
  ```
  in `_state_rate_polys`' docstring replace "the active conformation on the full graph,
  the inactive one on its pruned graph (`_inactive_mask`)," with "each conformation on
  its own pruned graph (`_state_mask`: the active one without its `:OnlyI` groups, the
  inactive one without its `:OnlyA` groups, each without the forms that strands)," and in
  `_mwc_state_polys`' docstring "The I-state polys come from the inactive conformation's
  pruned graph (`_inactive_mask` drops `:OnlyA` groups and every form they disconnect
  from free E)." with "Each state's polys come from its conformation's pruned graph
  (`_state_mask`: `_active_mask` drops the `:OnlyI` groups and `_inactive_mask` the
  `:OnlyA` ones, each with every form they disconnect from free E)."

  (e) `src/mechanism_enumeration.jl`, `_redundant_copy_groups`: replace the
  lines from `twin_active = _productive_twin(cm, trues(length(cm.edges)))` to the end of
  the function with
  ```julia
      if m isa Mechanism
          twin = _productive_twin(cm, trues(length(cm.edges)))
          return [g for g in eachindex(members)
                  if _gauge_rescaling(cm, members, g, twin, _one_state_class) !== nothing]
      end
      tags = cat_allo_states(m)
      masks = (_active_mask(m), _inactive_mask(m))
      kept = [[filter(k -> mask[k], ks) for ks in members] for mask in masks]
      twins = [_productive_twin(cm, mask) for mask in masks]
      unscaled(c) = Dict(h => (_NO_RESCALING, _NO_RESCALING)
                         for h in eachindex(kept[c]) if !isempty(kept[c][h]))
      filter(collect(eachindex(members))) do g
          all(c -> isempty(kept[c][g]), 1:2) && return false
          one_number(h) = h > 0 && tags[h] === :EqualAI && tags[g] === :EqualAI
          label(state) = (h, t) -> _GaugeClass(h, t, one_number(h) ? 0 : state)
          rescaling = [isempty(kept[c][g]) ? unscaled(c) :
                       _gauge_rescaling(cm, kept[c], g, twins[c], label(c)) for c in 1:2]
          any(isnothing, rescaling) && return false
          all(h -> tags[h] !== :EqualAI || !haskey(rescaling[2], h) ||
                   rescaling[1][h] == rescaling[2][h],
              keys(rescaling[1]))
      end
  ```
  and the docstring's paragraph that begins "An allosteric copy binds the active state
  always" up to "every factor of that state is 1." with "An allosteric copy binds the
  active state unless its tag is `:OnlyI` and the inactive state unless it is `:OnlyA`.
  Each state's graph is its law's (`_active_mask`, `_inactive_mask`): the edges its tags
  keep between forms they join to free E. The copy binds in a state only at the sites
  that graph keeps; in a state where it binds nothing, every factor of that state is 1.
  A copy that binds in neither state is not redundant." (the rest of the paragraph, from
  "The two states form one system.", stays). A mechanism without `:OnlyI` takes every
  edge for its active state as before, so its verdict does not move.

  (f) Same file, `_degenerate`: read the active conformation's law (AM A8). At 08c4043
  it reads the catalytic graph with its own rapid-equilibrium
  flags over every edge. Folding the active mask into those flags would not restrict the
  maximal-rate test, which turns rapid equilibrium every edge that takes up or gives off
  one of the names and so brings excluded edges back; the helpers take the mask as an edge
  set instead. Replace `_has_vmax` and its docstring with
  ```julia
  """Whether the rate of graph `m`, edge `k` at rapid equilibrium where `re[k]` holds, has a
  maximum as the metabolites `names` (its substrates or its products, sorted) grow: turning
  rapid equilibrium every steady-state edge that takes up or gives off one of them leaves no
  rapid-equilibrium turnover cycle (condition V). Only the edges `on` flags count, so an
  allosteric mechanism's test can read its active law (`_degenerate`)."""
  function _has_vmax(m::Mechanism, re::BitVector, names::Vector{Symbol},
                     on::BitVector = trues(length(m.edges)))
      side = _names_mask(m, names)
      touches = BitVector([(e.consumed | e.released) & side != 0 for e in m.edges])
      !_re_turnover_cycle(m, (re .| touches) .& on)
  end
  ```
  In `_chemistry_equilibrates_both_sides`, append to the docstring "Only the edges `on`
  flags count (`_degenerate`)."; replace its signature and first line
  ```julia
  function _chemistry_equilibrates_both_sides(m::Mechanism, re::BitVector)
      subs, prods = _names_mask(m, m.substrates), _names_mask(m, m.products)
  ```
  with
  ```julia
  function _chemistry_equilibrates_both_sides(m::Mechanism, re::BitVector,
                                              on::BitVector = trues(length(m.edges)))
      re = re .& on
      subs, prods = _names_mask(m, m.substrates), _names_mask(m, m.products)
  ```
  and its merged complexes: replace the line
  ```julia
      nodes = [s.segments; [[x] for x in _merged_complexes(m, subs) if s.segment_of[x] == 0]]
  ```
  with
  ```julia
      nodes = [s.segments; [[x] for x in _merged_complexes(m, subs, on)
                            if s.segment_of[x] == 0]]
  ```
  and `_merged_complexes` with
  ```julia
  _merged_complexes(m::Mechanism, subs::Counts, on::BitVector = trues(length(m.edges))) =
      unique(e.to for (k, e) in enumerate(m.edges)
             if on[k] && _is_binding(e) && e.consumed & subs != 0 && _is_chemistry(m, e))
  ```
  appending to its docstring "Only the edges `on` flags count." Every other read of an edge
  there is gated by `re[k]`, which now holds only the kept edges. Then replace
  `_degenerate` and its docstring with
  ```julia
  """Whether `m` is degenerate by the structural conditions V and C: no maximal rate in one
  direction (`_has_vmax`), or chemistry in rapid equilibrium with both sides
  (`_chemistry_equilibrates_both_sides`). Read on the catalytic graph with its own flags,
  over the active conformation's law: every edge of a plain mechanism, and the edges
  `_active_mask` keeps of an allosteric one. The beam expands a degenerate seed without
  fitting it (`_base_tier`)."""
  function _degenerate(m::Union{Mechanism, AllostericMechanism})
      cm = catalytic(m)
      re = _flags(cm)
      on = m isa AllostericMechanism ? _active_mask(m) : trues(length(cm.edges))
      !(_has_vmax(cm, re, cm.substrates, on) && _has_vmax(cm, re, cm.products, on)) ||
          _chemistry_equilibrates_both_sides(cm, re, on)
  end
  ```
  The `on` argument defaults to every edge, so the helpers' other callers (the seed moves,
  `_seed_candidate_screen` and `_seed_variants`, which call `_merged_complexes`, and TME's
  tests of `_has_vmax` and `_chemistry_equilibrates_both_sides`) keep their calls. A
  mechanism without `:OnlyI` keeps every edge in its active law, so its verdict does not
  move.

  (g) `src/admission.jl`: admission's turnover check (spec section 7, "Admission" item 3;
  RES §1.5's last interior line), the subject of Step 2 (b)'s essential-activator testset,
  watched failing in Step 3. In the docstring of `_Fault`, replace
  ```julia
  (`_ARC_CYCLE`), or forms free E's component does not reach along the tagged bindings
  (`_ISLAND`); on a face, in admission's face condition, a second closed class of the
  ```
  with
  ```julia
  (`_ARC_CYCLE`), forms free E's component does not reach along the tagged bindings
  (`_ISLAND`), or an active conformation that no longer turns the reaction over
  (`_NO_ACTIVE_TURNOVER`); on a face, in admission's face condition, a second closed class
  of the
  ```
  add the line `    _NO_ACTIVE_TURNOVER` to the `@enum _Fault` block after `    _BRIDGE`
  (before `    _DENOMINATOR_VANISHES`), and replace `_ADMISSION_FAULTS` with
  ```julia
  const _ADMISSION_FAULTS = (_NO_LIMIT, _ARC_CYCLE, _ISLAND, _NO_ACTIVE_TURNOVER,
                             _EPS_NONUNIQUE, _PUMPED, _BRIDGE)
  ```
  In `_admission_interior`'s docstring, replace
  ```julia
    pruning (`_ISLAND`).
  The verdict carries the conformation and the edge its message leads with (face 0).
  ```
  with
  ```julia
    pruning (`_ISLAND`);
  - turnover: when tags exclude edges of the active conformation, its law network still
    holds a cycle that turns the reaction over (`_NO_ACTIVE_TURNOVER`).
  The verdict carries the conformation and the edge its message leads with (face 0).
  ```
  and its last lines
  ```julia
          k == 0 || return _Verdict(_ISLAND, _LIMIT_UNTESTED, 0, c, k)
      end
      _WELL_POSED
  end
  ```
  with
  ```julia
          k == 0 || return _Verdict(_ISLAND, _LIMIT_UNTESTED, 0, c, k)
      end
      k = _first_excluded(am, 1)
      k != 0 && !_has_turnover_cycle(m, _law_mask(am, 1)) &&
          return _Verdict(_NO_ACTIVE_TURNOVER, _LIMIT_UNTESTED, 0, 1, k)
      _WELL_POSED
  end
  ```
  In `_admission_message`, insert before its final `else`:
  ```julia
      elseif v.fault === _NO_ACTIVE_TURNOVER
          ks = [j for j in eachindex(m.edges) if _excludes(am, j, 1)]
          return "AllostericMechanism: without its :OnlyI steps $(lines(ks)) the active " *
                 "conformation turns no substrate into product, so the enzyme never turns " *
                 "the reaction over. Tag their groups :EqualAI or :NonequalAI."
  ```

  Leave `src/dsl.jl`'s `@allosteric_mechanism` docstring and
  `docs/src/deriving/mwc_allostery.md`, which still say a catalytic step is never
  `:OnlyI`: Task 14 rewrites both paragraphs and expects their present text.

- [ ] **Step 5: Run them and watch them pass.** The Step 3 commands, logs
  `P4_T10_green_<file>.log`, then `test_joint_solve`, `test_rate_eq_derivation`,
  `allosteric_ground_truth` and `test_reaction_graph` the same way (logs
  `P4_T10_green_<file>.log`), which read the joint solve, the names and the state masks
  this task changed. Expected: `exit 0` for each, no warning, no pin moved. Record each
  file's testset time from Step 3's and this step's logs: the new testsets add at most
  2 s in all. Over that, report the times to Denis before any change; never drop a
  fixture to meet it. A law that misses the oracle by more than 1e-20, a detailed-balance
  residual over 1e-30 or a dimension that differs from the fitted count is a defect: run
  `superpowers:systematic-debugging` on that fixture (compare the joint solve's classes,
  `_kept_members`, with the helper's `_testhelper_classes`), never loosen the assertion.

- [ ] **Step 6: Census against Task 9's.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  $S/phase3/jl.sh $S/phase45/logs/P4_T10_census.log 1800 --project=$S/phase3/census/env_guard --heap-size-hint=2500M $S/phase3/census/corpus.jl $S/phase45/guard_wt $S/phase3/census/P4_T10 $S/phase3/census/rename_phase3.toml
  $S/phase3/jl.sh $S/phase45/logs/P4_T10_census_cmp.log 600 --project=$S/phase3/census/env --heap-size-hint=2500M $S/phase3/census/compare.jl $S/phase3/census/P4_T9 $S/phase3/census/P4_T10 $S/phase3/census/rename_phase3.toml
  ```
  Expected: `RESULT PASS` with zero differences: no enumerated or spec mechanism has a
  catalytic `:OnlyI` group, and without one every mask, column, class and name above is
  unchanged. A difference stops the task: run `superpowers:systematic-debugging` on the
  first one.

- [ ] **Step 7: Run the full suite.** As Task 8 Step 11, log `P4_T10_full.log`. Expected:
  `exit 0`, pristine, Aqua and JET pass, the `rate_equation` performance gates (0
  allocations, under 120 ns) and the compile budgets hold.

- [ ] **Step 8: Commit.** Write Step 5's time in place of `<t>` below, then:

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && git add src/types.jl src/admission.jl src/thermo.jl src/names.jl src/rate_eq_derivation.jl src/mechanism_enumeration.jl test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_oracles.jl test/test_admission.jl test/test_types.jl test/test_dsl.jl test/test_identify_rate_equation.jl test/test_mechanism_enumeration.jl test/test_reaction_graph.jl && git commit -F - <<'EOF2'
  Let a catalytic group bind modifiers in the inactive conformation only

  A catalytic group may be :OnlyI when its steps take up or give off modifiers alone;
  the constructor rejects any other step with the fix. The active law drops such a group
  and the forms it strands, the joint solve gives the group inactive constants and each
  dropped active copy a class of its own, and an :EqualAI constant the inactive law
  alone reads stays in the law. The copy rule builds each conformation from its own law,
  and the degeneracy test reads the active one through an edge mask, which the
  maximal-rate and chemistry tests now take. Admission gains its turnover check, which
  only an :OnlyI group can fail. Both tag checks now read _TAG_CODES.

  R4.4's thirteen mechanisms: the law equals the deletion network's exact steady state
  to 1e-20, detailed balance holds to 1e-30, and the declared dimension equals the fitted
  count. Admission's :OnlyI verdicts are tested: the B1 set's fast cycles on the
  ε-network, the bridges, the essential activator, fz101_1964's second closed class, the
  interior faults in the active conformation, and a mechanism that only joint R1
  rejects (P4-8b), whose epsilon-path does not converge on the two faces it rejects.
  Census equals Task 9's. New testsets: <t> s.

  Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
  EOF2
  ```

---

### Evidence for Tasks 11-14

The expected children and verdicts of Task 12's fixtures were computed before the plan, with
Task 0(b)'s own implementation of RES §1.9 (`$S/phase45/research/task0/b/p41v3.jl`: B v3
admission, the research guard with the mod-p residual) and the research oracles (`vi/vor.jl`
at 256 bits, RE-speed on every face, the ε-path along 32 corners and five orders). Scripts
and logs: `$S/phase45/plan/probes/p4/` (`tagging.jl`, `tagging2.jl`, `tagging3.jl`,
`tagging4.jl`, `scan.py`, with `.log` beside each; `gen/fixtures.py` wrote the fixture
blocks). Findings that shape the tasks:

- `census/base6` reproduces Task 0(b)'s parents: pingpong 894, bibi 987, terbi 4,885 init
  seeds plus 31 level-1 parents (`scan.py`). Exactly 12 bi-bi parents fail existence with
  every chemistry group `:OnlyA`; their interior-passing K-types (6, 12, 24 per parent) equal
  Task 0(b)'s P4-1 counts. 757 ter-bi parents are ping-pong; P4-1 leaves 98 childless, 22 of
  them without a binding group (`tagging2.log`), as Task 0(b) found.
- The tagging-count testsets keep their children: uni-uni 3, 7, 9; uni-bi 7; 31; TC 3;
  merged 7 and 7; random-order 3; bi-bi sequential 15 (`tagging.log`).
- The six-step bi-bi ping-pong gets 6 K-types, 3 per half-reaction base. Class 3 rejects the
  three on the first half at zero substrates and zero B (stage `split`, limit `differ`), so
  three survive the exit; the oracles agree (`mismatch` on those faces, `agree` and
  `converges` elsewhere and on the kept children; `tagging.log`, `tagging3.log`).
- On that ping-pong at one subunit, the one-half V-type's chemistry relaxation (both halves
  `:NonequalAI`) passes admission and fails class 3 (`split`, `speed`): its second
  half-reaction is rapid-equilibrium chemistry, so the covalent forms of both conformations
  share one segment and split at zero substrates (`tagging.log`, change-tag of `FP2_base`,
  and `tagging3.log`, `FP2 chem relaxed`: FP2 is this ping-pong written with ATP, F16BP, F6P
  and ADP for A, P, B and Q). On the three-step ping-pong 8fe71fded09a27ec the same
  relaxation passes (`tagging4.log`). Change-tag reaches a fully live ping-pong only where
  class 3 allows.
- Task 8's exit testset parent, census ping-pong 3fdbe6e1dc001bb7 with regulator R at
  n = 2, gets six children from RES §1.9's rules, all passing the research guard: a K-type
  per half-reaction base (the first half with P's binding, the second with Q's) and two
  V-types per base. The other seven K-type candidates fail admission (two `:pumped`, five
  islands). On bi-bi 1291abb91f537c6b's K-type with Q's binding `:OnlyA`, change-tag's Q
  relaxation fails existence (G1) while its P and B-at-E(Q) relaxations pass the guard
  (`tagging5.jl`, `tagging5.log`, written while this plan was checked).
- Admission gives every tag set of `test_types.jl`'s `:OnlyA` testsets the verdict
  `_onlya_offenders` gave: `_NO_LIMIT` where it named an offender, `_NO_FAULT` where it
  found none (`tagging3.log`).
- The pingpong census reaction's 18 ping-pong seeds get two bases (14 seeds) or one (4);
  each base tags one half-reaction's groups (`tagging2.log`).
- Bi-bi parent 1291abb91f537c6b's dropped K-type and its every-chemistry set fail existence
  `inside` (`_NO_LIMIT`); ter-bi parent 0dc399e7dab6e4f5's K-types fail `:pumped` on its
  three bases and stay an island on the every-chemistry set (`tagging2.log`).

---

### Task 11: The tagging move's `:pumped` drops through the oracles (P4-9c)

A measurement outside the suite. Task 0(b) found that B v3 drops 5,342 tagging candidates as
`:pumped` (RES §5.2 item 3) and no oracle checked them; the 2026-10-14 outage took its candidate
lists (P4-9c). This task rebuilds RES §1.9's candidates from the package's own admission
(Tasks 4-5) on Task 0(b)'s parents and runs both test oracles (Task 1) on every `:pumped` one.
Each candidate the oracles find well-posed is an R-B2 false alarm in enumeration (P4-9b): the
ROADMAP Idea of Task 5 cites the count, and the PR body lists the parents (Task 14).

The package gets P4-1's move only in Task 12, so the script builds the candidates itself, as
RES §1.9 states them, and records which ones the move must emit. Task 12 then runs the
package's move on the same parents and compares (stage `move`), which tests the move against
an implementation written apart from it.

**Files:**
- Create (scratch, outside the repository): `$S/phase45/pumped/pumped_check.jl` (stages
  `candidates`, `oracles`, `move`), `$S/phase45/pumped/run_queue.sh` (resumable chunks),
  `$S/phase45/pumped/summarize.py` (merge and report), `$S/phase45/pumped/research_parent.jl`
  (Task 0(b)'s research code on one parent, for Step 6)
- Data: `$S/phase45/pumped/data/{parents,candidates,children,oracles}.tsv`, merged from the
  chunk files `data/<stage>_<rxn>_<lo>.tsv`; `data/false_alarms.tsv`
- Logs: `$S/phase45/logs/P4_T11_<stage>_<rxn>_<lo>.log`, `P4_T11_report.log`,
  `P4_T11_full.log`
- Modify: `ROADMAP.md` (Task 5's Idea "Admit the equilibrium bridges and the vanishing-fed
  basins that admission rejects" gains the measured count, P4-9b)

**Interfaces:**
- Consumes (Tasks 2-5, 8): `EnzymeRates._admission(am)::_Verdict`,
  `_admission_face(am, eps, law, faces, i)::_Verdict`, `_eps_network(am)`, `_law_network(am)`,
  `_faces(m)::Vector{UInt64}`, `_face_names(m, face)`, `_law_mask(am, c)::BitVector`,
  `_has_turnover_cycle(m, mask)::Bool`, `_guard(am)::_Verdict`, `_is_well_posed(m)::Bool`,
  `_Fault` values `_NO_FAULT`, `_PUMPED`. (Phase 3) `_chemistry_groups(m)::BitVector` and
  `_hyperbolic_catalysis(m)`; `_subsets`, `allowed_catalytic_multiplicities`,
  `regulators`, `regulator`, `name`, `AllostericRegulator`, `Met`, `RawForm`, `RawStep`,
  `RawMechanism`, `Mechanism(raw)`, the composition constructor
  `AllostericMechanism(m, tags, n, sites)`, `RegulatorySite(ligands, n, tags)`,
  `_has_residual`. (Task 1, fixtures file) `_testhelper_faces(m)`,
  `_testhelper_respeed_verdicts(m; rng, law = :limit)`,
  `_testhelper_epath_verdict(m; rng, face, corners = true)`. (Census) `census_common.jl`'s
  `read_rows`, `mech_id`, `census_text`; `neutral_graph.jl`'s `census_nmech`; the corpus
  `census/base6` (the five reactions of `census/base`, which the outage took; base6 holds the
  same parents: pingpong 894, bibi 987, terbi 4,885 init + 31 level-1, checked 2026-10-16).
- Produces: `candidates.tsv` (one row per candidate of RES §1.9: the V-type set `vt`, every
  stopping chemistry set when `vt` fails, every K-type on `kbases`), `children.tsv` (the
  children RES §1.9 emits, with the exit verdict of each network), `parents.tsv`,
  `oracles.tsv`, `false_alarms.tsv`; the count of R-B2 false alarms and their parents for the
  PR body; the `move` stage that Task 12 runs.

**Gate:** the `:pumped` drops reproduce Task 0(b)'s 5,342: 4,808 K-type (pingpong 352 at 52
parents, ter-bi 4,456 at 667) and 534 base candidates (pingpong 38, ter-bi 496, tried and
superset sets together, as Task 0's `verify_b/spot.jl` part 2 counted them); every difference
is explained (a parent the guard now rejects, or a package verdict that differs from the
research code's, which Step 6 traces). The emitted children reproduce P4-19a: 88,434 children
on 88,174 networks (pingpong 4,878 = 4,358 K + 520 V; bi-bi 15,479; ter-bi 68,077); class 3
drops 820 networks at the exit (pingpong 56 = 46 K + 10 V; ter-bi 764); bi-bi's 12 parents
whose `vt` fails existence get 198 K-types on `vt`; each of pingpong's 30 ping-pong parents
that allow n = 2 keeps 2-4 V-types past the exit; 98 of ter-bi's 757 ping-pong parents end
childless (22 without a binding group, 76 losing every K-type to `_ISLAND` or `_PUMPED`).
Every `:pumped` candidate runs the RE-speed oracle on every face and the ε-path at 128
corners on each `:pumped` face; each well-posed one is listed in `false_alarms.tsv`. About
1.5 h of Julia in chunks of at most 20 min and 3 GB, resumable by chunk. The commit carries
only the ROADMAP change.

- [ ] **Step 1: Confirm the branch state.** Task 10 is committed and its full suite passed.

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && git status --short && git log --oneline -1
  grep -c . $S/phase3/census/base6/mechs.tsv
  ```
  Expected: a clean tree whose last commit is Task 10's; `81498` (base6's mechanisms).

- [ ] **Step 2: Write the script.** Create `$S/phase45/pumped/pumped_check.jl`:

  ```julia
  # ABOUTME: Task 11 of the phase-4 plan: builds research section 1.9's tagging candidates
  # ABOUTME: from package admission and checks the ones it drops as :pumped with the oracles.
  #
  # Run from the guard worktree through phase3/jl.sh, one stage at a time:
  #   candidates <rxn> <lo> <hi> <budget_s>   parents lo:hi of <rxn> (bibi, pingpong, terbi)
  #   oracles <lo> <hi> <budget_s>            :pumped rows lo:hi of data/candidates.tsv
  #   move <rxn> <lo> <hi> <budget_s>         the package's tagging move (Task 12) on parents
  # A stage writes data/<stage>_<rxn>_<lo>.tsv (oracles: data/oracles_all_<lo>.tsv) and
  # prints "NEXT <i>", the first index it left undone. The budget counts from the process
  # start, loading included, so a chunk ends inside phase3/jl.sh's timeout.
  const T_LOAD = time()
  const S = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
  const CENSUS = joinpath(S, "phase3", "census")
  const BASE = joinpath(CENSUS, "base6")
  const DATA = joinpath(@__DIR__, "data")
  using TestEnv; TestEnv.activate()
  using Test, EnzymeRates, LinearAlgebra, Random, Printf
  const ER = EnzymeRates
  include(joinpath(pwd(), "test", "mechanism_definitions_for_test_enzyme_derivation.jl"))
  include(joinpath(CENSUS, "census_common.jl"))   # read_rows, mech_id, census_text
  include(joinpath(CENSUS, "neutral_graph.jl"))   # census_nmech of a reaction graph

  const MECHS = Dict(r[1] => (neutral = r[2], oriented = r[3])
                     for r in read_rows(joinpath(BASE, "mechs.tsv")))
  const RXN_TEXT = Dict(r[1] => r[5] for r in read_rows(joinpath(BASE, "reactions.tsv")))

  """The substrates and the products an `@enzyme_reaction` text declares, each sorted."""
  function reactants(text)
      names(key) = sort([Symbol(strip(split(x, "[")[1]))
                         for x in split(match(Regex("\\s$key:\\s*(.*)"), text)[1], ",")])
      (names("substrates"), names("products"))
  end

  """
  The plain parents Task 0(b) ran (research section 5.1), as `(level, id)` sorted: every
  plain census parent of `which` (bibi levels 0-2, pingpong 0-1), and for terbi every
  `init_mechanisms` seed of the census too.
  """
  function parent_list(which)
      level = Dict{String, Int}()
      for r in read_rows(joinpath(BASE, "parents.tsv"))
          r[1] == which && startswith(MECHS[r[3]].neutral, "plain") || continue
          level[r[3]] = min(get(level, r[3], typemax(Int)), parse(Int, r[2]))
      end
      if which == "terbi"
          for r in read_rows(joinpath(BASE, "init.tsv"))
              r[1] == which && (level[r[3]] = 0)
          end
      end
      sort([(l, id) for (id, l) in level])
  end

  """
  The raw mechanism of a plain census oriented text (`census/README.md`, "Neutral text"):
  each step `from-side -> to-side` in its stored direction, a label `Xmod` for a reactant
  `X` bound as a modifier, `residual=plusA_minusP` for a residual form.
  """
  function census_raw(text, subs::Vector{Symbol}, prods::Vector{Symbol})
      startswith(text, "plain ") || error("census_raw: not a plain mechanism: $text")
      reactant(x::Symbol) = x in subs || x in prods
      function met(label::AbstractString)
          stem = Symbol(chopsuffix(label, "mod"))
          endswith(label, "mod") && reactant(stem) && return ER.Met(stem, true)
          x = Symbol(label)
          ER.Met(x, !reactant(x))
      end
      function form(t::AbstractString)
          i = findfirst('(', t)
          parts = split(t[(i + 1):(end - 1)], "; residual=")
          bound = ER.Met[met(x) for x in split(parts[1], ",") if !isempty(x)]
          runs = length(parts) > 1 ? split(parts[2], "_") : SubString{String}[]
          ER.RawForm(Symbol(t[1:(i - 1)]), bound,
                     Symbol[Symbol(r[5:end]) for r in runs if startswith(r, "plus")],
                     Symbol[Symbol(r[6:end]) for r in runs if startswith(r, "minus")])
      end
      steps, group = ER.RawStep[], Int[]
      for (g, block) in enumerate(split(text[(length("plain ") + 1):end], " / "))
          for st in split(block[2:(end - 1)], " | ")
              l, r = split(st[1:(end - 3)], " -> ")
              ls, rs = split(l, " + "), split(r, " + ")
              push!(steps, ER.RawStep(form(ls[1]), form(rs[1]),
                                      ER.Met[met(x) for x in ls[2:end]],
                                      ER.Met[met(x) for x in rs[2:end]], endswith(st, " RE")))
              push!(group, g)
          end
      end
      mods = sort!(unique!([x.name for s in steps
                            for x in [s.consumed; s.released; s.from.bound; s.to.bound]
                            if x.modifier]))
      ER.RawMechanism(copy(subs), copy(prods), mods, steps, group)
  end

  """The tags of `n` groups: those in `onlya` `:OnlyA`, every other one `:EqualAI`."""
  function tagset(n::Int, onlya)
      t = fill(:EqualAI, n)
      t[collect(Int, onlya)] .= :OnlyA
      t
  end
  allo(m, tags) = ER.AllostericMechanism(m, tags, 1, ER.RegulatorySite[])
  """Whether `tags` stop every catalytic cycle of the inactive conformation (stops, §1.9)."""
  stops(m, tags) = !ER._has_turnover_cycle(m, ER._law_mask(allo(m, tags), 2))
  """The nonempty subsets of `v` by size, each size in `_subsets`' order."""
  by_size(v) = [s for k in 1:length(v) for s in ER._subsets(v)[2:end] if length(s) == k]
  tagtext(t) = join(string.(t), ",")

  """
  Research section 1.9's candidates of the plain mechanism `m`, as rows `(kind, src, tags,
  v)` with `v` the package's admission verdict (`_admission`, class 3 not run): the V-type
  set `vt` (every chemistry group `:OnlyA`); when `vt` fails, every set of chemistry groups
  that stops the inactive conformation's cycles, by size (`src` "tried" when no admitted set
  lies inside it, else "superset"); then the K-types, every nonempty set of binding groups
  added to each base and, when `vt` fails, to `vt` too (`src` names the base). Also the
  bases (the chemistry sets the V-types use) and whether `vt` passed.
  """
  function candidates(m)
      chem = ER._chemistry_groups(m)
      n, cg, bg = length(chem), findall(chem), findall(!, chem)
      rows = NamedTuple[]
      isempty(cg) && return (rows = rows, bases = Vector{Int}[], vt_ok = false)
      add!(kind, src, tags) = push!(rows, (kind = kind, src = src, tags = tags,
                                           v = ER._admission(allo(m, tags))))
      add!("vt", "-", tagset(n, cg))
      vt_ok = rows[1].v.fault === ER._NO_FAULT
      bases = vt_ok ? [cg] : Vector{Int}[]
      if !vt_ok
          for sel in by_size(cg)
              t = tagset(n, sel)
              stops(m, t) || continue
              tried = !any(b -> b ⊆ sel, bases)
              add!("base", tried ? "tried" : "superset", t)
              tried && rows[end].v.fault === ER._NO_FAULT && push!(bases, sel)
          end
      end
      kbases = vt_ok ? [("today", cg)] :
               [[("base:" * join(b, ","), b) for b in bases]; ("vt", cg)]
      for (src, b) in kbases, sel in by_size(bg)
          add!("K", src, tagset(n, [b; sel]))
      end
      (rows = rows, bases = bases, vt_ok = vt_ok)
  end

  """Whether the move emits candidate row `r` (a passing K-type, or a passing tried base)."""
  emitted(r) = r.v.fault === ER._NO_FAULT && (r.kind == "K" || r.kind == "vt" ||
                                              (r.kind == "base" && r.src == "tried"))

  """The faces (indices into `_faces`) on which admission's face condition finds a stranded
  basin that it cannot seal (`_PUMPED`)."""
  function pumped_faces(am)
      eps, law = ER._eps_network(am), ER._law_network(am)
      faces = ER._faces(am.catalytic)
      [i for i in eachindex(faces)
       if ER._admission_face(am, eps, law, faces, i).fault === ER._PUMPED]
  end

  """The reaction `which` and what its parents need: substrates, products, the declared
  allosteric regulators."""
  function reaction(which)
      rxn = Core.eval(Main, Meta.parse(RXN_TEXT[which]))
      subs, prods = reactants(RXN_TEXT[which])
      regs = [ER.name(ER.regulator(rm)) for rm in ER.regulators(rxn)
              if ER.regulator(rm) isa ER.AllostericRegulator]
      (rxn = rxn, subs = subs, prods = prods, regs = regs)
  end

  """Parent `id` of reaction `R` as the package's graph, checked to round-trip its id."""
  function census_parent(R, id)
      m = ER.Mechanism(census_raw(MECHS[id].oriented, R.subs, R.prods))
      mech_id(census_text(m)) == id || error("parent $id does not round-trip its census id")
      m
  end

  multiplicities(R, m) = [cn for cn in ER.allowed_catalytic_multiplicities(R.rxn)
                          if cn == 1 || ER._hyperbolic_catalysis(m)]

  function run_candidates(which, lo, hi, budget)
      R = reaction(which)
      plist = parent_list(which)
      io = Dict(k => open(joinpath(DATA, "$(k)_$(which)_$(lo).tsv"), "w")
                for k in ("parents", "candidates", "children"))
      t0, next = time(), min(hi, length(plist)) + 1
      for i in lo:min(hi, length(plist))
          time() - T_LOAD > budget && (next = i; break)
          level, id = plist[i]
          m = census_parent(R, id)
          cns = multiplicities(R, m)
          c = candidates(m)
          exits = Dict{Vector{Symbol}, String}()
          verdict(t) = get!(() -> string(ER._guard(allo(m, t)).fault), exits, t)
          for r in c.rows
              pf = r.v.fault === ER._PUMPED ? join(pumped_faces(allo(m, r.tags)), ",") : "-"
              names = r.v.face == 0 ? "-" :
                      join(ER._face_names(m, ER._faces(m)[r.v.face]), ",")
              # The child's census id at one subunit and no site keys the row for the
              # research code's cross-check (research_parent.jl, Step 6).
              println(io["candidates"], join([which, id, i, r.kind, r.src, tagtext(r.tags),
                  string(r.v.fault), r.v.face, names, r.v.conformation, emitted(r), pf,
                  mech_id(census_text(allo(m, r.tags)))], '\t'))
          end
          ktags = [r.tags for r in c.rows if r.kind == "K" && emitted(r)]
          vtags = [tagset(length(ER._chemistry_groups(m)), b) for b in c.bases]
          nK = nV = drops = 0
          for cn in cns
              for t in ktags
                  println(io["children"],
                          join([which, id, "K", tagtext(t), cn, "-", verdict(t)], '\t'))
                  nK += 1
                  drops += verdict(t) != "_NO_FAULT"
              end
              for t in vtags, reg in R.regs, tag in (:OnlyA, :OnlyI)
                  println(io["children"], join([which, id, "V", tagtext(t), cn, "$reg::$tag",
                                                verdict(t)], '\t'))
                  nV += 1
                  drops += verdict(t) != "_NO_FAULT"
              end
          end
          chem = ER._chemistry_groups(m)
          kind = !any(chem) ? "nochem" : c.vt_ok ? "today" : "minimal"
          vt_fault = isempty(c.rows) ? "-" : string(c.rows[1].v.fault)
          println(io["parents"], join([which, i, level, id, any(ER._has_residual, m.forms),
              count(chem), count(!, chem), vt_fault, kind, length(c.bases),
              join(cns, ","), nK, nV, drops, ER._is_well_posed(m)], '\t'))
          foreach(flush, values(io))
          i % 200 == 0 && (@printf("[%.0f s] %s %d\n", time() - t0, which, i); flush(stdout))
      end
      foreach(close, values(io))
      @printf("DONE candidates %s %d:%d in %.0f s\n", which, lo, next - 1, time() - t0)
      println("NEXT ", next)
  end

  # A `vt` row's tags equal those of the `base` row of the full chemistry set, which
  # `candidates` also writes whenever `vt` fails, so the oracles skip `vt` rows.
  function run_oracles(lo, hi, budget)
      rows = [split(l, '\t') for l in eachline(joinpath(DATA, "candidates.tsv"))]
      sel = [i for (i, r) in enumerate(rows) if r[7] == "_PUMPED" && r[4] != "vt"]
      Rs = Dict(w => reaction(w) for w in ("pingpong", "terbi", "bibi"))
      io = open(joinpath(DATA, "oracles_all_$(lo).tsv"), "w")
      t0, next = time(), min(hi, length(sel)) + 1
      for j in lo:min(hi, length(sel))
          time() - T_LOAD > budget && (next = j; break)
          i = sel[j]
          which, id, _, kind, src, tags = rows[i][1:6]
          am = allo(census_parent(Rs[which], id), Symbol.(split(tags, ",")))
          rng = Random.Xoshiro(hash((id, tags)))
          faces = _testhelper_faces(am)
          respeed = _testhelper_respeed_verdicts(am; rng, law = :limit)
          pumped = parse.(Int, split(rows[i][12], ","))
          for (f, face) in enumerate(faces)
              epath = f in pumped ?
                      _testhelper_epath_verdict(am; rng, face, corners = true) : :-
              println(io, join([i, which, id, kind, src, tags, f,
                                join(sort!(collect(face)), ","), respeed[f], epath], '\t'))
          end
          flush(io)
          j % 25 == 0 &&
              (@printf("[%.0f s] %d/%d\n", time() - t0, j, length(sel)); flush(stdout))
      end
      close(io)
      @printf("DONE oracles %d:%d in %.0f s\n", lo, next - 1, time() - t0)
      println("NEXT ", next)
  end

  """A child of the tagging move as `children.tsv` writes it: tags, multiplicity, site."""
  childkey(am) = (tagtext(am.cat_allo_states), string(am.catalytic_multiplicity),
                  isempty(am.regulatory_sites) ? "-" :
                  string(only(am.regulatory_sites).ligands[1], "::",
                         only(am.regulatory_sites).allo_states[1]))

  function run_move(which, lo, hi, budget)
      R = reaction(which)
      plist = parent_list(which)
      expected = Dict{String, Vector{NTuple{3, String}}}()
      for l in eachline(joinpath(DATA, "children.tsv"))
          r = split(l, '\t')
          r[1] == which || continue
          push!(get!(expected, r[2], NTuple{3, String}[]), (r[4], r[5], r[6]))
      end
      io = open(joinpath(DATA, "move_$(which)_$(lo).tsv"), "w")
      t0, next = time(), min(hi, length(plist)) + 1
      for i in lo:min(hi, length(plist))
          time() - T_LOAD > budget && (next = i; break)
          id = plist[i][2]
          kids = childkey.(ER._expand_to_allosteric(census_parent(R, id), R.rxn))
          want = get(expected, id, NTuple{3, String}[])
          println(io, join([which, id, length(want), length(kids), Set(kids) == Set(want),
                            allunique(kids), join(sort(string.(setdiff(want, kids))), ";"),
                            join(sort(string.(setdiff(kids, want))), ";")], '\t'))
          flush(io)
      end
      close(io)
      @printf("DONE move %s %d:%d in %.0f s\n", which, lo, next - 1, time() - t0)
      println("NEXT ", next)
  end

  @printf("loaded in %.0f s\n", time() - T_LOAD)
  mkpath(DATA)
  if ARGS[1] == "candidates"
      run_candidates(ARGS[2], parse(Int, ARGS[3]), parse(Int, ARGS[4]),
                     parse(Float64, ARGS[5]))
  elseif ARGS[1] == "oracles"
      run_oracles(parse(Int, ARGS[2]), parse(Int, ARGS[3]), parse(Float64, ARGS[4]))
  elseif ARGS[1] == "move"
      run_move(ARGS[2], parse(Int, ARGS[3]), parse(Int, ARGS[4]), parse(Float64, ARGS[5]))
  else
      error("unknown stage $(ARGS[1])")
  end
  ```

  The candidate stage follows RES §1.9 where Task 0's `task0/b/p41v3.jl` (a script that
  predates the outage) differs in order only: it walks subsets in lexicographic order within a
  size, this script in `_subsets`' order, which is the order Task 12's `_tagging_bases` reads.
  Inclusion-minimal sets do not depend on the order within a size. The base rows include the
  full chemistry set, as Task 0's count of 534 base drops did (`verify_b/spot.jl` part 2), so
  a failing `vt` appears twice, as `vt` and as a base. `census_common.jl` loads `SHA` and
  `TOML` and this script `Printf`: stdlibs, which load from `@stdlib` under TestEnv's
  environment. `chopsuffix` needs Julia 1.8; the worktree runs 1.12.

- [ ] **Step 3: Write the queue and the summary.** Create `$S/phase45/pumped/run_queue.sh`:

  ```bash
  #!/bin/bash
  # ABOUTME: Runs one stage of pumped_check.jl over an index range in resumable chunks of at
  # ABOUTME: most 20 min (budget 1080 s), one Julia at a time through phase3/jl.sh.
  # Usage: run_queue.sh <candidates|move> <rxn> <lo> <hi>  |  run_queue.sh oracles - <lo> <hi>
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
  STAGE=$1; RXN=$2; LO=$3; HI=$4
  cd $S/phase45/guard_wt || exit 1
  lo=$LO
  while [ "$lo" -le "$HI" ]; do
    PFX=$([ "$STAGE" = move ] && echo P4_T12 || echo P4_T11)
    LOG=$S/phase45/logs/${PFX}_${STAGE}_${RXN}_${lo}.log
    if [ "$STAGE" = oracles ]; then
      $S/phase3/jl.sh $LOG 1200 --project=. --heap-size-hint=2500M \
        $S/phase45/pumped/pumped_check.jl oracles $lo $HI 1080
    else
      $S/phase3/jl.sh $LOG 1200 --project=. --heap-size-hint=2500M \
        $S/phase45/pumped/pumped_check.jl $STAGE $RXN $lo $HI 1080
    fi
    next=$(grep -o "^NEXT [0-9]*" $LOG | awk '{print $2}')
    if [ -z "$next" ] || [ "$next" -le "$lo" ]; then echo "QUEUE STOP: no progress in $LOG"; exit 1; fi
    lo=$next
  done
  echo "QUEUE DONE $STAGE $RXN $LO:$HI"
  ```

  and `$S/phase45/pumped/summarize.py`:

  ```python
  # ABOUTME: Task 11 of the phase-4 plan: merges pumped_check.jl's chunk files and reports the
  # ABOUTME: :pumped drops, the emitted children and the oracles' verdicts against Task 0(b).
  # Usage: python3 -I summarize.py merge | report
  import collections, glob, os, sys
  D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

  def merge():
      for stage in ('parents', 'candidates', 'children', 'move'):
          files = sorted(glob.glob(os.path.join(D, f'{stage}_*_*.tsv')),
                         key=lambda f: (os.path.basename(f).split('_')[1],
                                        int(os.path.basename(f).split('_')[-1][:-4])))
          if not files:
              continue
          with open(os.path.join(D, f'{stage}.tsv'), 'w') as out:
              for f in files:
                  out.write(open(f).read())
          print(f'{stage}.tsv: {sum(1 for _ in open(os.path.join(D, stage + ".tsv")))} rows '
                f'from {len(files)} files')
      files = sorted(glob.glob(os.path.join(D, 'oracles_all_*.tsv')),
                     key=lambda f: int(f.split('_')[-1][:-4]))
      if files:
          with open(os.path.join(D, 'oracles.tsv'), 'w') as out:
              for f in files:
                  out.write(open(f).read())

  def rows(name):
      p = os.path.join(D, name)
      return [l.rstrip('\n').split('\t') for l in open(p)] if os.path.exists(p) else []

  def report():
      par, cand, kids = rows('parents.tsv'), rows('candidates.tsv'), rows('children.tsv')
      pr = print
      for rxn in ('pingpong', 'bibi', 'terbi'):
          P = [p for p in par if p[0] == rxn]
          C = [c for c in cand if c[0] == rxn]
          K = [k for k in kids if k[0] == rxn]
          pr(f'== {rxn}: {len(P)} parents, {sum(p[14] == "false" for p in P)} the guard rejects')
          for kind in ('vt', 'K', 'base'):
              pumped = [c for c in C if c[3] == kind and c[6] == '_PUMPED']
              pr(f'  {kind} :pumped {len(pumped)} at {len({c[1] for c in pumped})} parents'
                 + (f' (tried {sum(c[4] == "tried" for c in pumped)}, '
                    f'superset {sum(c[4] == "superset" for c in pumped)})' if kind == 'base' else ''))
          nets = {(k[1], k[3]) for k in K}
          pr(f'  children {len(K)} (K {sum(k[2] == "K" for k in K)}, V {sum(k[2] == "V" for k in K)}) '
             f'on {len(nets)} networks')
          dropped = {(k[1], k[3], k[2]) for k in K if k[6] != '_NO_FAULT'}
          pr(f'  exit drops: {len(dropped)} networks (K {sum(d[2] == "K" for d in dropped)}, '
             f'V {sum(d[2] == "V" for d in dropped)}); faults '
             f'{dict(collections.Counter(k[6] for k in K if k[6] != "_NO_FAULT"))}')
          if rxn == 'bibi':
              fail = {p[3] for p in P if p[7] not in ('_NO_FAULT', '-')}
              on_vt = {(c[1], c[5]) for c in C if c[3] == 'K' and c[4] == 'vt'}
              onvt = [k for k in K if k[2] == 'K' and (k[1], k[3]) in on_vt]
              pr(f'  parents whose vt fails: {len(fail)} '
                 f'({dict(collections.Counter(p[7] for p in P if p[3] in fail))}); '
                 f'K-types on vt {len(onvt)} at {len({k[1] for k in onvt})} parents')
          pp = [p for p in P if p[4] == 'true']
          if rxn == 'pingpong':
              open_ = [p for p in pp if '2' in p[10].split(',')]
              keep = collections.Counter(sum(1 for k in K if k[1] == p[3] and k[2] == 'V'
                                             and k[6] == '_NO_FAULT') for p in open_)
              pr(f'  ping-pong parents allowing n = 2: {len(open_)}; V-types past the exit '
                 f'per parent: {dict(sorted(keep.items()))}')
          if rxn == 'terbi':
              childless = [p for p in pp if int(p[11]) + int(p[12]) == 0]
              pr(f'  ping-pong parents {len(pp)}; childless {len(childless)} '
                 f'(no binding group {sum(p[6] == "0" for p in childless)})')
      O = rows('oracles.tsv')
      if not O:
          return
      by = collections.defaultdict(list)
      for o in O:
          by[int(o[0])].append(o)
      verdicts, wp = collections.Counter(), []
      for i, os_ in by.items():
          bad_re = sorted({o[8] for o in os_ if o[8] != 'agree'})
          bad_ep = sorted({o[9] for o in os_ if o[9] not in ('-', 'converges')})
          v = ('class3:' + ','.join(bad_re)) if bad_re else \
              ('ill_posed:' + ','.join(bad_ep)) if bad_ep else 'well_posed'
          verdicts[(os_[0][1], os_[0][3], v.split(':')[0])] += 1
          if v == 'well_posed':
              wp.append(os_[0])
      npumped = sum(1 for c in cand if c[6] == '_PUMPED' and c[3] != 'vt')
      pr(f'== oracles: {len(by)} of {npumped} :pumped candidates checked')
      for k, n in sorted(verdicts.items()):
          pr(f'  {k[0]} {k[1]} {k[2]}: {n}')
      with open(os.path.join(D, 'false_alarms.tsv'), 'w') as out:
          for o in wp:
              out.write('\t'.join([o[0], o[1], o[2], o[3], o[4], o[5]]) + '\n')
      pr(f'== R-B2 false alarms: {len(wp)} at {len({o[2] for o in wp})} parents '
         f'(data/false_alarms.tsv)')
      M = rows('move.tsv')
      if M:
          pr(f'== move: {len(M)} parents, {sum(m[4] == "true" for m in M)} with equal child sets, '
             f'{sum(m[5] == "false" for m in M)} with a repeated child')
          for m in M:
              if m[4] != 'true':
                  pr('  DIFF', m[0], m[1], 'expected-only', m[6], 'package-only', m[7])

  {'merge': merge, 'report': report}[sys.argv[1]]()
  ```

  `chmod +x $S/phase45/pumped/run_queue.sh` and `mkdir -p $S/phase45/pumped/data $S/phase45/logs`.

- [ ] **Step 4: Smoke-run the candidate stage.** Five parents of each reaction:

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  for r in pingpong bibi terbi; do $S/phase3/jl.sh $S/phase45/logs/P4_T11_smoke_$r.log 1200 --project=. --heap-size-hint=2500M $S/phase45/pumped/pumped_check.jl candidates $r 1 5 1080; done
  tail -3 $S/phase45/logs/P4_T11_smoke_*.log
  ```
  Expected: each log ends `DONE candidates <r> 1:5 in <t> s`, `NEXT 6`, `exit 0`; no
  "does not round-trip". `data/parents_<r>_1.tsv` has 5 rows, the last field `true` on each
  unless the guard rejects that parent (Step 6 explains such a parent). Delete the smoke
  files before Step 5:
  `rm $S/phase45/pumped/data/*_1.tsv`.

- [ ] **Step 5: Run the candidate stage on every parent.** In the background:

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  Q=$S/phase45/pumped/run_queue.sh
  ($Q candidates pingpong 1 894 && $Q candidates bibi 1 987 && $Q candidates terbi 1 4916) > $S/phase45/logs/P4_T11_candidates_queue.log 2>&1
  ```
  Poll `P4_T11_candidates_queue.log` until it ends `QUEUE DONE candidates terbi 1:4916`
  (the three queues take about 15 min of Julia). Then:

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  python3 -I $S/phase45/pumped/summarize.py merge && python3 -I $S/phase45/pumped/summarize.py report | tee $S/phase45/logs/P4_T11_report.log
  ```
  Expected report lines (Task 0(b), RES §5.1 and §5.3; a `vt :pumped` row is also counted
  among the base rows, as Task 0 counted it, and a ping-pong parent's `vt` fails as an
  island, not as `:pumped`):
  ```
  == pingpong: 894 parents, 0 the guard rejects
    vt :pumped <v> at <v> parents
    K :pumped 352 at 52 parents
    base :pumped 38 at <p> parents (tried <a>, superset <b>)
    children 4878 (K 4358, V 520) on 4618 networks
    exit drops: 56 networks (K 46, V 10); faults {...}
    ping-pong parents allowing n = 2: 30; V-types past the exit per parent: {2: .., 3: .., 4: ..}
  == bibi: 987 parents, 0 the guard rejects
    vt :pumped <v> at <v> parents
    K :pumped 0 at 0 parents
    base :pumped 0 at 0 parents (tried 0, superset 0)
    children 15479 (K 15479, V 0) on 15479 networks
    exit drops: 0 networks (K 0, V 0); faults {}
    parents whose vt fails: 12 ({'_NO_LIMIT': 12}); K-types on vt 198 at 12 parents
  == terbi: 4916 parents, 0 the guard rejects
    vt :pumped <v> at <v> parents
    K :pumped 4456 at 667 parents
    base :pumped 496 at <p> parents (tried <a>, superset <b>)
    children 68077 (K 68077, V 0) on 68077 networks
    exit drops: 764 networks (K 764, V 0); faults {...}
    ping-pong parents 757; childless 98 (no binding group 22)
  ```
  The pingpong V-types per parent lie in {2, 3, 4} and their counts sum to 30. Every
  `faults` entry names a class-3 fault (`_SPLIT` with the limit test `differ`, RES §5.1).

- [ ] **Step 6: Explain every difference from Task 0(b).** For each count that differs:
  - a parent the guard rejects (the parents' last field `false`): its candidates do not count
    as the package's, since enumeration never tags it; give its fault
    (`EnzymeRates._guard(m)`) in the report;
  - otherwise, take the first parent whose counts differ and run Task 0(b)'s research code
    on it with `$S/phase45/pumped/research_parent.jl` (below; `task0/b/move.jl` itself
    rebuilds its parents from `census/base`, which the outage took). It prints one line per
    candidate, `kind  child-id  verdict  [tried|superset]`, the child id being the census id
    of the candidate at one subunit and no site, which `candidates.tsv` holds in its last
    column, so the two lists join on it. The research verdict names the rule: `ok`, `G1`
    (`_NO_LIMIT` or `_ARC_CYCLE`), `island` (`_ISLAND`), `face_<x>` for a face fault
    (`face_pumped` is `_PUMPED`). A candidate whose package verdict differs from RES §1.5's
    rule is a defect in Tasks 4-5: STOP and debug it (superpowers:systematic-debugging)
    before going on; a difference the research script makes against RES §1 (some predate
    the outage) is recorded as such.

  `$S/phase45/pumped/research_parent.jl`:

  ```julia
  # ABOUTME: Task 11's cross-check: Task 0(b)'s research code (p41v3.jl, B v3 admission) on
  # ABOUTME: one census parent, each candidate keyed by its child's census id.
  # Usage, in the census environment: research_parent.jl <bibi|pingpong|terbi> <parent id>
  const S = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
  const RS = joinpath(S, "phase45", "research")
  const BASE = joinpath(S, "phase3", "census", "base6")
  using Printf, Random
  include(joinpath(RS, "vi", "compare.jl"))           # modules VOr and R, 256 bits
  include(joinpath(RS, "task0", "a", "bv3.jl"))       # bv3_face, class2_eps
  include(joinpath(RS, "task0", "b", "p41v3.jl"))     # admit_v3, rtext

  const WHICH, ID = ARGS[1], ARGS[2]
  const ROLES = R.reaction_roles(only(r for r in R.census_rows(joinpath(BASE,
                                                                        "reactions.tsv"))
                                      if r[1] == WHICH)[5])
  const NET = R.census_net(only(r[3] for r in R.census_rows(joinpath(BASE, "mechs.tsv"))
                                if r[1] == ID), ROLES; name = ID)
  child(t) = R.child_net(NET, t, 1, R.Site[])
  cid(t) = R.mech_id(rtext(child(t)))
  verdict(t) = (v = admit_v3(child(t)); v[1] ? :ok : v[2])

  G = R.ngroups(NET)
  cg, bg = findall(R.chem_groups(NET)), findall(!, R.chem_groups(NET))
  vt = R.tagset(G, cg)
  vt_ok = verdict(vt) === :ok
  println("vt\t", cid(vt), "\t", verdict(vt))
  bases = vt_ok ? [cg] : Vector{Int}[]
  for k in 1:length(cg), sel in sort(R.subsets_of(cg, k))
      vt_ok && break
      t = R.tagset(G, sel)
      R.copy_turns_over(child(t), 2) && continue
      v, tried = verdict(t), !any(b -> issubset(b, sel), bases)
      tried && v === :ok && push!(bases, sel)
      println("base\t", cid(t), "\t", v, "\t", tried ? "tried" : "superset")
  end
  kbases = vt_ok ? bases : [bases; [cg]]
  for b in kbases, k in 1:length(bg), sel in sort(R.subsets_of(bg, k))
      t = R.tagset(G, [b; sel])
      println("K\t", cid(t), "\t", verdict(t))
  end
  ```
  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad
  $S/phase3/jl.sh $S/phase45/logs/P4_T11_research_<id>.log 900 --project=$S/phase3/census/env --heap-size-hint=2500M $S/phase45/pumped/research_parent.jl <rxn> <id>
  awk -F'\t' -v id=<id> '$2 == id { print $4 "\t" $13 "\t" $7 "\t" $5 }' $S/phase45/pumped/data/candidates.tsv
  ```
  Record the explanation in `$S/phase45/logs/P4_T11_report.log` under `== differences`.

- [ ] **Step 7: Run the oracles on every `:pumped` candidate.** The row count `N` is the
  sum of the K and base `:pumped` lines of Step 5, 5,342 when Task 0(b)'s counts reproduce;
  the oracles skip the `vt` rows, each a repeat of a base row. In the background:

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  N=<the sum above>; $S/phase45/pumped/run_queue.sh oracles - 1 $N > $S/phase45/logs/P4_T11_oracles_queue.log 2>&1
  ```
  About 1.5 h of Julia in five or six chunks of 18 min; poll the queue log in bounded loops
  until `QUEUE DONE oracles - 1:<N>`.
  A chunk that stops (`QUEUE STOP`) is rerun from its `NEXT`. Peak RSS per chunk under 3 GB
  (`maxrss_kb` in each chunk log).

- [ ] **Step 8: Report the oracles' verdicts.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  python3 -I $S/phase45/pumped/summarize.py merge && python3 -I $S/phase45/pumped/summarize.py report | tee $S/phase45/logs/P4_T11_report.log
  ```
  Expected: `== oracles: N of N :pumped candidates checked`, one line per (reaction,
  kind, verdict) with `class3`, `ill_posed` and `well_posed` counts, and
  `== R-B2 false alarms: <F> at <G> parents`. A `class3` candidate fails the RE-speed oracle
  on some face (class 3 would drop it at the exit anyway); an `ill_posed` one fails the
  ε-path on a `:pumped` face, which confirms the rejection; a `well_posed` one agrees on every
  face and converges on every `:pumped` face along all 128 corners and five orders: an R-B2
  false alarm. `F` is the count the ROADMAP and the PR body cite.

- [ ] **Step 9: Cite the count in ROADMAP.** In `ROADMAP.md`, in the Ideas entry
  "**Admit the equilibrium bridges and the vanishing-fed basins that admission rejects.**",
  replace

  ```markdown
    mechanisms, none among the census mechanisms or their `:OnlyI` variants. No exact
  ```
  with, when `F` is 0,
  ```markdown
    mechanisms, none among the census mechanisms or their `:OnlyI` variants. Tag sets that
    stop the inactive conformation's catalysis one half-reaction at a time, built on the
    census ping-pong and ter-bi parents and filtered by admission, lose 5,342 candidates
    to such basins, and the RE-speed and ε-path oracles confirm every one. No exact
  ```
  and when `F` is positive, with (writing `F` and `G` as numbers)
  ```markdown
    mechanisms, none among the census mechanisms or their `:OnlyI` variants. Tag sets that
    stop the inactive conformation's catalysis one half-reaction at a time, built on the
    census ping-pong and ter-bi parents and filtered by admission, lose 5,342 candidates
    to such basins; the RE-speed and ε-path oracles find F of them, at G parents,
    well-posed: R-B2 false alarms. No exact
  ```
  Write `N` for 5,342 if the two differ. Keep the entry's lines under 92 columns. The
  wording describes the measurement, not the package's move, which Task 12 brings; Task 12
  Step 11 names the move once it filters this way.

- [ ] **Step 10: Run the full suite.** `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T11_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'` in the background. Expected: `exit 0`, the pass count of Task 10's
  run, Aqua and JET pass, no warning; `maxrss_kb` within 5% of Task 10's.

- [ ] **Step 11: Commit.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && git add ROADMAP.md && git commit -F - <<'EOF'
  Count the tagging candidates admission drops as pumped that are well-posed

  The oracles checked every candidate tag set that admission's face condition drops as a
  stranded, unsealed basin on the census ping-pong and ter-bi parents (research section
  1.9 candidates built from the package's admission): <N> candidates, <F> of them
  well-posed (R-B2 false alarms). The ROADMAP idea for an exact face condition cites the
  count. The check runs outside the suite ($S/phase45/pumped); its candidate lists match
  Task 0(b)'s counts <or: differ as follows, each explained>.

  Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
  EOF
  ```

---

### Task 12: Tagging move (P4-1, P4-19) and change-tag through admission (RES §4 row 10)

The tagging move builds on the bases of RES §1.9: today's V-type set when admission accepts it,
otherwise each inclusion-minimal set of chemistry groups that stops the inactive conformation's
catalysis and passes admission (reading (a), P4-19). K-types build on every base and also on
today's set when admission rejects it; V-types on every base. Change-tag keeps a relaxation when
admission accepts it (P4-20). `_onlya_offenders`, `_onlya_cone_offenders` and
`_has_strict_positive_combination` lose their last callers and go; `_rational_nullspace` stays
for `_integer_nullspace`.

**Files:**
- Modify: `src/mechanism_enumeration.jl` (new `_admits`, `_stops_inactive_turnover`,
  `_tagging_bases`; `_expand_to_allosteric` and its docstring; `_expand_change_allo_state`'s
  filter and docstring)
- Modify: `src/thermo.jl` (delete `_onlya_offenders`, `_onlya_cone_offenders`,
  `_has_strict_positive_combination` and their docstrings)
- Modify: `src/types.jl` (`_chemistry_groups(m)` loses its `group` and count arguments, whose
  only other caller was `_onlya_offenders`; `_is_chemistry`'s docstring)
- Modify: `ROADMAP.md` (retire the Known issue Task 8 renamed "The tagging and change-tag
  moves misjudge `:OnlyA` groups that mix a binding and a chemistry step."; rewrite the
  `BigInt` entry Task 8 rewrote, which names `_onlya_offenders`; the both-halves entry
  states what the tagging move does; a Known issue on the six-step bi-bi ping-pong's
  change-tag, open item 45)
- Test: `test/test_mechanism_enumeration.jl` (Task 8's exit testset re-pinned),
  `test/test_types.jl`, `test/test_oracles.jl`, `test/test_reaction_graph.jl` (the census
  population pins P4-1 moves, re-pinned in Step 10)
- Scratch: `$S/phase45/probes/P4_T12_value_check.jl`, `P4_T12_split_parents.jl`,
  `$S/phase45/census/explain_tagging.py`; logs `$S/phase45/logs/P4_T12_*.log`; census
  `census/P4_T12`

**Interfaces:**
- Consumes: Task 5's `_admission(am)::_Verdict`; Task 2's `_law_mask(am, c)`,
  `_has_turnover_cycle(m, mask)`; Tasks 2-5's `_guard`, `_faces`, `_face_names`,
  `_LIMIT_DIFFER`, `_LIMIT_SPEED`, `_NO_LIMIT`, `_ISLAND`, `_PUMPED`; Task 8's
  `_is_well_posed`; Task 11's
  `data/children.tsv`, `data/parents.tsv` and the `move` stage of `pumped_check.jl`;
  `_hyperbolic_catalysis`, `_subsets`, `_chemistry_groups`, `_groups_with`, `ngroups`,
  `_with(am; states)`; `allowed_catalytic_multiplicities`, `regulators`,
  `regulator`, `name`, `AllostericRegulator`, `RegulatorySite`, `_is_binding`,
  `_bound_entry`; Task 1's `@_testhelper_unguarded`, `_testhelper_faces`,
  `_testhelper_respeed_verdicts`, `_testhelper_epath_verdict`; Task 3's
  `_testhelper_face_faults`, `_testhelper_agrees`; Task 5's `_testhelper_admission_agrees`;
  Task 8's TME testset "enumeration drops the children the guard rejects" (its
  Step 2), TME testset "AllostericMechanism — a split keeps its parent's admission verdict"
  (Step 3(j)) and TT `_onlya_offenders` assertions (Step 3(c)-(f)); the census tools
  (`corpus.jl`, `compare.jl`, the environment `census/env_guard` and corpus `census/P4_T10`
  of Task 10).
- Produces:
  ```julia
  _admits(m::Mechanism, tags::Vector{Symbol})::Bool                   # ok(t), RES §1.9
  _stops_inactive_turnover(m::Mechanism, tags::Vector{Symbol})::Bool  # stops(t)
  _tagging_bases(m::Mechanism)::Vector{Vector{Symbol}}   # [vt] when ok(vt), else the
                                                         # inclusion-minimal t_S (reading a)
  _expand_to_allosteric(m::Mechanism, rxn::EnzymeReaction)::Vector{AllostericMechanism}
  _expand_change_allo_state(am::AllostericMechanism)::Vector{AllostericMechanism}
  _chemistry_groups(m::Mechanism)::BitVector
  const _testhelper_pingpong_census_rxn   # test_oracles.jl: P6's census ping-pong reaction
  ```

**Gate:** the package move reproduces Task 11's `children.tsv` on every one of its parents
(`move` stage: equal child sets on all 6,797) and so Task 0(b)'s counts (P4-19a): 88,434
children on 88,174 networks (pingpong 4,878 = 4,358 K + 520 V; bi-bi 15,479; ter-bi 68,077);
class 3 drops 820 networks at the exit (pingpong 56, ter-bi 764); K-types on today's V-type set
give 198 bi-bi children at 12 parents; each of pingpong's 30 regulated ping-pong parents allowing
n = 2 keeps 2-4 V-types; 98 of 757 ter-bi ping-pong parents end childless (22 without a
binding group, 76 losing every K-type to `_ISLAND` or `_PUMPED`). Census `P4_T12` against
`P4_T10`: every difference classified as the tagging move's or change-tag's, none unexplained.
Every new child of the re-pinned fixtures passes the value check (both oracles). RES §3.3 P6
passes. `test_mechanism_enumeration.jl`, `test_types.jl` and `test_oracles.jl` together grow by
at most 3 s. Denis confirms P4-19's counts in the PR (Task 14).

- [ ] **Step 1: Write the failing tests in `test/test_mechanism_enumeration.jl`.**

  (a) Replace the testset "Mechanism — Bi-bi ping-pong: binding-group :OnlyA variants only"
  (section `_expand_to_allosteric`) with:

  ```julia
  @testset "Mechanism — Bi-bi ping-pong: binding-group :OnlyA variants only" begin
      # SEED: bi-bi ping-pong, 6 kinetic groups: 4 bindings (A and Q at free E, B and P at
      # the covalent form E(; residual = A - P)) and the two half-reactions E(A) <-->
      # E(P; residual = A - P) and E(B; residual = A - P) ⇌ E(Q). Tagging both halves
      # :OnlyA strands the covalent forms behind chemistry in the inactive conformation,
      # an island admission rejects, so the move builds on one base per half: that half
      # :OnlyA, the other running in both conformations (P4-1). No regulator is declared,
      # so no V-type: the children are K-types, a nonempty set of bindings :OnlyA on a
      # base, each kept when admission accepts it. On the first half's base it keeps the A
      # and P bindings, while B strands the covalent forms and Q pumps a basin; on the
      # second half's base it keeps B and Q, while A pumps a basin and P strands forms.
      em_seed = @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(; residual = A - P) + B ⇌ E(B; residual = A - P)
              E + Q ⇌ E(Q)
              E(; residual = A - P) + P ⇌ E(P; residual = A - P)
              E(A) <--> E(P; residual = A - P)
              E(B; residual = A - P) ⇌ E(Q)
          end
      end
      m = EnzymeRates.Mechanism(em_seed)
      _testhelper_assert_mechanism_invariants(m)
      bi_bi_pp_allo_rxn = @enzyme_reaction begin
          substrates: A[CX], B[N]
          products: P[C], Q[NX]
          oligomeric_state: 2
      end
      a_half1 = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 2
          steps: begin
              E + A ⇌ E(A)                                          :: OnlyA
              E(; residual = A - P) + B ⇌ E(B; residual = A - P)    :: EqualAI
              E + Q ⇌ E(Q)                                          :: EqualAI
              E(; residual = A - P) + P ⇌ E(P; residual = A - P)    :: EqualAI
              E(A) <--> E(P; residual = A - P)                      :: OnlyA
              E(B; residual = A - P) ⇌ E(Q)                         :: EqualAI
          end
      end
      p_half1 = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 2
          steps: begin
              E + A ⇌ E(A)                                          :: EqualAI
              E(; residual = A - P) + B ⇌ E(B; residual = A - P)    :: EqualAI
              E + Q ⇌ E(Q)                                          :: EqualAI
              E(; residual = A - P) + P ⇌ E(P; residual = A - P)    :: OnlyA
              E(A) <--> E(P; residual = A - P)                      :: OnlyA
              E(B; residual = A - P) ⇌ E(Q)                         :: EqualAI
          end
      end
      ap_half1 = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 2
          steps: begin
              E + A ⇌ E(A)                                          :: OnlyA
              E(; residual = A - P) + B ⇌ E(B; residual = A - P)    :: EqualAI
              E + Q ⇌ E(Q)                                          :: EqualAI
              E(; residual = A - P) + P ⇌ E(P; residual = A - P)    :: OnlyA
              E(A) <--> E(P; residual = A - P)                      :: OnlyA
              E(B; residual = A - P) ⇌ E(Q)                         :: EqualAI
          end
      end
      b_half2 = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 2
          steps: begin
              E + A ⇌ E(A)                                          :: EqualAI
              E(; residual = A - P) + B ⇌ E(B; residual = A - P)    :: OnlyA
              E + Q ⇌ E(Q)                                          :: EqualAI
              E(; residual = A - P) + P ⇌ E(P; residual = A - P)    :: EqualAI
              E(A) <--> E(P; residual = A - P)                      :: EqualAI
              E(B; residual = A - P) ⇌ E(Q)                         :: OnlyA
          end
      end
      q_half2 = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 2
          steps: begin
              E + A ⇌ E(A)                                          :: EqualAI
              E(; residual = A - P) + B ⇌ E(B; residual = A - P)    :: EqualAI
              E + Q ⇌ E(Q)                                          :: OnlyA
              E(; residual = A - P) + P ⇌ E(P; residual = A - P)    :: EqualAI
              E(A) <--> E(P; residual = A - P)                      :: EqualAI
              E(B; residual = A - P) ⇌ E(Q)                         :: OnlyA
          end
      end
      bq_half2 = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 2
          steps: begin
              E + A ⇌ E(A)                                          :: EqualAI
              E(; residual = A - P) + B ⇌ E(B; residual = A - P)    :: OnlyA
              E + Q ⇌ E(Q)                                          :: OnlyA
              E(; residual = A - P) + P ⇌ E(P; residual = A - P)    :: EqualAI
              E(A) <--> E(P; residual = A - P)                      :: EqualAI
              E(B; residual = A - P) ⇌ E(Q)                         :: OnlyA
          end
      end

      result = EnzymeRates._expand_to_allosteric(m, bi_bi_pp_allo_rxn)
      @test length(result) == 6
      @test Set(result) == Set([a_half1, p_half1, ap_half1, b_half2, q_half2, bq_half2])
      # On the first half's base class 3 rejects the faces with every substrate zero and
      # with B zero: a split whose limit the speeds do not move but which differs from the
      # law (limit test `differ`), so expand_mechanisms drops those three children.
      for c in (a_half1, p_half1, ap_half1)
          v = EnzymeRates._guard(c)
          names = EnzymeRates._face_names(m, EnzymeRates._faces(m)[v.face])
          @test (v.fault, v.limit, names) ==
                (EnzymeRates._SPLIT, EnzymeRates._LIMIT_DIFFER, [:A, :B])
      end
      kept = [b_half2, q_half2, bq_half2]
      @test all(EnzymeRates._is_well_posed, kept)
      @test Set(c for c in EnzymeRates.expand_mechanisms([m], bi_bi_pp_allo_rxn)
                if c isa EnzymeRates.AllostericMechanism) == Set(kept)
      # Δ params: each kept child adds L only; every tag is :EqualAI or :OnlyA.
      @test _testhelper_param_deltas(m, kept) == fill(1, 3)
      foreach(_testhelper_assert_mechanism_invariants, kept)
  end
  ```

  (b) Replace the testset ":OnlyA to_allosteric emits dead-inactive combos only" with this
  testset, which also tests that the tagging move gives a ping-pong one V-type base per
  half:

  ```julia
  @testset "the tagging move gives a ping-pong one V-type base per half" begin
      # The three-step ping-pong of census pingpong seed 8fe71fded09a27ec: free E takes up
      # A and gives off P in one Theorell–Chance step, leaving the covalent form E(;
      # residual = A - P), which takes up B into E(Q) (fused); Q's release is the one
      # binding group. Tagging both half-reactions :OnlyA strands the covalent form, an
      # island, so each half alone is a base (P4-1). The V-types pair each base with the
      # declared regulator R. The one K-type tags Q's binding on the second half's base; on
      # the first half's base Q's binding pumps a basin, and on the island it stays one.
      m = EnzymeRates.Mechanism(@enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A <--> E(; residual = A - P) + P
              E + Q <--> E(Q)
              E(; residual = A - P) + B <--> E(Q)
          end
      end)
      rxn = @enzyme_reaction begin
          substrates: A[CX], B[N]
          products: P[C], Q[NX]
          allosteric_regulators: R
      end
      half1 = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              E + A <--> E(; residual = A - P) + P :: OnlyA
              E + Q <--> E(Q)                      :: EqualAI
              E(; residual = A - P) + B <--> E(Q)  :: EqualAI
          end
      end
      half2 = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              E + A <--> E(; residual = A - P) + P :: EqualAI
              E + Q <--> E(Q)                      :: EqualAI
              E(; residual = A - P) + B <--> E(Q)  :: OnlyA
          end
      end
      k_q_half2 = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              E + A <--> E(; residual = A - P) + P :: EqualAI
              E + Q <--> E(Q)                      :: OnlyA
              E(; residual = A - P) + B <--> E(Q)  :: OnlyA
          end
      end
      v_half1_a = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          allosteric_regulators: R
          steps: begin
              E + A <--> E(; residual = A - P) + P :: OnlyA
              E + Q <--> E(Q)                      :: EqualAI
              E(; residual = A - P) + B <--> E(Q)  :: EqualAI
          end
          regulatory_site(1): R::OnlyA
      end
      v_half1_i = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          allosteric_regulators: R
          steps: begin
              E + A <--> E(; residual = A - P) + P :: OnlyA
              E + Q <--> E(Q)                      :: EqualAI
              E(; residual = A - P) + B <--> E(Q)  :: EqualAI
          end
          regulatory_site(1): R::OnlyI
      end
      v_half2_a = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          allosteric_regulators: R
          steps: begin
              E + A <--> E(; residual = A - P) + P :: EqualAI
              E + Q <--> E(Q)                      :: EqualAI
              E(; residual = A - P) + B <--> E(Q)  :: OnlyA
          end
          regulatory_site(1): R::OnlyA
      end
      v_half2_i = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          allosteric_regulators: R
          steps: begin
              E + A <--> E(; residual = A - P) + P :: EqualAI
              E + Q <--> E(Q)                      :: EqualAI
              E(; residual = A - P) + B <--> E(Q)  :: OnlyA
          end
          regulatory_site(1): R::OnlyI
      end
      k_q_half1 = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              E + A <--> E(; residual = A - P) + P :: OnlyA
              E + Q <--> E(Q)                      :: OnlyA
              E(; residual = A - P) + B <--> E(Q)  :: EqualAI
          end
      end
      @test Set(EnzymeRates.AllostericMechanism(m, t, 1, EnzymeRates.RegulatorySite[])
                for t in EnzymeRates._tagging_bases(m)) == Set([half1, half2])
      kids = EnzymeRates._expand_to_allosteric(m, rxn)
      @test length(kids) == 5
      @test Set(kids) == Set([k_q_half2, v_half1_a, v_half1_i, v_half2_a, v_half2_i])
      @test all(EnzymeRates._is_well_posed, kids)
      @test EnzymeRates._admission(k_q_half1).fault === EnzymeRates._PUMPED
      both = [:OnlyA, :OnlyA, :OnlyA]   # every group :OnlyA (the K-type on vt)
      @test EnzymeRates._admission(EnzymeRates.AllostericMechanism(m, both, 1,
          EnzymeRates.RegulatorySite[])).fault === EnzymeRates._ISLAND
      @test Set(c for c in EnzymeRates.expand_mechanisms([m], rxn)
                if c isa EnzymeRates.AllostericMechanism) == Set(kids)
  end
  ```

  (c) Replace the testset "change_allo_state drops partial-catalysis relaxations" with:

  ```julia
  @testset "change_allo_state relaxes a ping-pong K-type's binding, not its chemistry" begin
      # The tagging move's K-type of the three-step ping-pong (census pingpong seed
      # 8fe71fded09a27ec): the second half-reaction and Q's binding :OnlyA, the first half
      # running in both conformations. Q's binding relaxes, and admission accepts it; the
      # chemistry does not relax while a binding is :OnlyA. One child.
      k_type = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              E + A <--> E(; residual = A - P) + P :: EqualAI
              E + Q <--> E(Q)                      :: OnlyA
              E(; residual = A - P) + B <--> E(Q)  :: OnlyA
          end
      end
      q_relaxed = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              E + A <--> E(; residual = A - P) + P :: EqualAI
              E + Q <--> E(Q)                      :: NonequalAI
              E(; residual = A - P) + B <--> E(Q)  :: OnlyA
          end
      end
      kids = EnzymeRates._expand_change_allo_state(k_type)
      @test length(kids) == 1
      @test Set(kids) == Set([q_relaxed])
      @test EnzymeRates._is_well_posed(k_type) && EnzymeRates._is_well_posed(q_relaxed)
  end
  ```

  (d) Replace the testset "change_allo_state drops mixed-iso V-type relaxations (no :OnlyA
  binding)" with:

  ```julia
  @testset "change_allo_state relaxes a one-half V-type, leaving class 3 to the exit" begin
      # The tagging move's V-type on the second half-reaction of the six-step bi-bi
      # ping-pong (testset "Mechanism — Bi-bi ping-pong: …" above), at one subunit: that
      # half :OnlyA, every binding :EqualAI, R at an :OnlyA site. Each binding relaxes
      # alone, the two half-reactions relax together, and the site ligand relaxes;
      # admission accepts every catalytic relaxation. With both halves :NonequalAI the
      # conformations run
      # different chemistry on one rapid-equilibrium segment, and with every substrate zero
      # or B zero the covalent forms hold the enzyme in both, joined only through free E: a
      # split whose limit depends on the speeds. Class 3 rejects it, so expand_mechanisms
      # drops that child.
      rxn = @enzyme_reaction begin
          substrates: A[CX], B[N]
          products: P[C], Q[NX]
          allosteric_regulators: R
      end
      v_type = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          allosteric_regulators: R
          steps: begin
              E + A ⇌ E(A)                                       :: EqualAI
              E(; residual = A - P) + B ⇌ E(B; residual = A - P) :: EqualAI
              E + Q ⇌ E(Q)                                       :: EqualAI
              E(; residual = A - P) + P ⇌ E(P; residual = A - P) :: EqualAI
              E(A) <--> E(P; residual = A - P)                   :: EqualAI
              E(B; residual = A - P) ⇌ E(Q)                      :: OnlyA
          end
          regulatory_site(1): R::OnlyA
      end
      a_relaxed = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          allosteric_regulators: R
          steps: begin
              E + A ⇌ E(A)                                       :: NonequalAI
              E(; residual = A - P) + B ⇌ E(B; residual = A - P) :: EqualAI
              E + Q ⇌ E(Q)                                       :: EqualAI
              E(; residual = A - P) + P ⇌ E(P; residual = A - P) :: EqualAI
              E(A) <--> E(P; residual = A - P)                   :: EqualAI
              E(B; residual = A - P) ⇌ E(Q)                      :: OnlyA
          end
          regulatory_site(1): R::OnlyA
      end
      b_relaxed = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          allosteric_regulators: R
          steps: begin
              E + A ⇌ E(A)                                       :: EqualAI
              E(; residual = A - P) + B ⇌ E(B; residual = A - P) :: NonequalAI
              E + Q ⇌ E(Q)                                       :: EqualAI
              E(; residual = A - P) + P ⇌ E(P; residual = A - P) :: EqualAI
              E(A) <--> E(P; residual = A - P)                   :: EqualAI
              E(B; residual = A - P) ⇌ E(Q)                      :: OnlyA
          end
          regulatory_site(1): R::OnlyA
      end
      q_relaxed = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          allosteric_regulators: R
          steps: begin
              E + A ⇌ E(A)                                       :: EqualAI
              E(; residual = A - P) + B ⇌ E(B; residual = A - P) :: EqualAI
              E + Q ⇌ E(Q)                                       :: NonequalAI
              E(; residual = A - P) + P ⇌ E(P; residual = A - P) :: EqualAI
              E(A) <--> E(P; residual = A - P)                   :: EqualAI
              E(B; residual = A - P) ⇌ E(Q)                      :: OnlyA
          end
          regulatory_site(1): R::OnlyA
      end
      p_relaxed = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          allosteric_regulators: R
          steps: begin
              E + A ⇌ E(A)                                       :: EqualAI
              E(; residual = A - P) + B ⇌ E(B; residual = A - P) :: EqualAI
              E + Q ⇌ E(Q)                                       :: EqualAI
              E(; residual = A - P) + P ⇌ E(P; residual = A - P) :: NonequalAI
              E(A) <--> E(P; residual = A - P)                   :: EqualAI
              E(B; residual = A - P) ⇌ E(Q)                      :: OnlyA
          end
          regulatory_site(1): R::OnlyA
      end
      chemistry_relaxed = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          allosteric_regulators: R
          steps: begin
              E + A ⇌ E(A)                                       :: EqualAI
              E(; residual = A - P) + B ⇌ E(B; residual = A - P) :: EqualAI
              E + Q ⇌ E(Q)                                       :: EqualAI
              E(; residual = A - P) + P ⇌ E(P; residual = A - P) :: EqualAI
              E(A) <--> E(P; residual = A - P)                   :: NonequalAI
              E(B; residual = A - P) ⇌ E(Q)                      :: NonequalAI
          end
          regulatory_site(1): R::OnlyA
      end
      site_relaxed = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          allosteric_regulators: R
          steps: begin
              E + A ⇌ E(A)                                       :: EqualAI
              E(; residual = A - P) + B ⇌ E(B; residual = A - P) :: EqualAI
              E + Q ⇌ E(Q)                                       :: EqualAI
              E(; residual = A - P) + P ⇌ E(P; residual = A - P) :: EqualAI
              E(A) <--> E(P; residual = A - P)                   :: EqualAI
              E(B; residual = A - P) ⇌ E(Q)                      :: OnlyA
          end
          regulatory_site(1): R::NonequalAI
      end
      kids = EnzymeRates._expand_change_allo_state(v_type)
      @test length(kids) == 6
      @test Set(kids) == Set([a_relaxed, b_relaxed, q_relaxed, p_relaxed,
                              chemistry_relaxed, site_relaxed])
      v = EnzymeRates._guard(chemistry_relaxed)
      cm = chemistry_relaxed.catalytic
      names = EnzymeRates._face_names(cm, EnzymeRates._faces(cm)[v.face])
      @test (v.fault, v.limit, names) ==
            (EnzymeRates._SPLIT, EnzymeRates._LIMIT_SPEED, [:A, :B])
      kept = [a_relaxed, b_relaxed, q_relaxed, p_relaxed, site_relaxed]
      @test all(EnzymeRates._is_well_posed, kept)
      exits = EnzymeRates.expand_mechanisms([v_type], rxn)
      @test chemistry_relaxed ∉ exits && all(in(exits), kept)
      # Every kept child turns the reaction over once products are present. At zero
      # products its rate is 0: the inactive conformation's first half-reaction fills
      # the covalent forms, which only P's binding drains while the second half is
      # stopped.
      for k in kept
          cm = EnzymeRates.compile_mechanism(k)
          fp = EnzymeRates.fitted_params(cm)
          prm = NamedTuple{(fp..., :Keq, :E_total)}(((1.3 for _ in fp)..., 3.0, 1.0))
          @test EnzymeRates.rate_equation(cm, (A = 1.3, B = 1.3, P = 0.0, Q = 0.0,
                                               R = 1.0), prm) == 0
          @test EnzymeRates.rate_equation(cm, (A = 1.3, B = 1.3, P = 0.1, Q = 0.1,
                                               R = 1.0), prm) > 0
      end
  end
  ```

  (e) Add, after the testset "_expand_to_allosteric on a random-order merged bi-bi led by
  plain bindings":

  ```julia
  @testset "K-types build on today's V-type set when admission rejects it" begin
      # Census bi-bi parent 1291abb91f537c6b (census/base6), one of the 12 bi-bi parents
      # whose V-type set fails admission (P4-19a). A's group mixes its plain bindings at E
      # and E(P) with its fused binding into E(P, Q), and B's mixes its plain binding at E
      # with its fused one, so both are chemistry. Tagged :OnlyA with every binding group
      # :EqualAI, B's plain binding at E keeps both sides joined, through Q's bindings and
      # B's binding at E(Q), so its affinity has no path to zero. No smaller set of
      # chemistry groups stops the inactive conformation's catalysis, so there is no base,
      # and the K-types build on that set itself: a tagged Q binding or B binding at E(Q)
      # cuts the join and restores the limit. Of the seven candidates admission keeps six;
      # P's bindings alone leave B's binding at E joined.
      m = EnzymeRates.Mechanism(@enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              (E + A ⇌ E(A), E(P) + A ⇌ E(A, P), E(B) + A ⇌ E(P, Q))
              (E + B <--> E(B), E(A) + B <--> E(P, Q))
              (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q))
              (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(P) + Q <--> E(P, Q))
              E(Q) + B <--> E(B, Q)
          end
      end)
      rxn = @enzyme_reaction begin
          substrates: A[C], B[N]
          products: P[C], Q[N]
      end
      q_tagged = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              (E + A ⇌ E(A), E(P) + A ⇌ E(A, P), E(B) + A ⇌ E(P, Q))          :: OnlyA
              (E + B <--> E(B), E(A) + B <--> E(P, Q))                        :: OnlyA
              (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q)) :: EqualAI
              (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(P) + Q <--> E(P, Q)) :: OnlyA
              E(Q) + B <--> E(B, Q)                                           :: EqualAI
          end
      end
      b_at_eq_tagged = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              (E + A ⇌ E(A), E(P) + A ⇌ E(A, P), E(B) + A ⇌ E(P, Q))          :: OnlyA
              (E + B <--> E(B), E(A) + B <--> E(P, Q))                        :: OnlyA
              (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q)) :: EqualAI
              (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(P) + Q <--> E(P, Q)) :: EqualAI
              E(Q) + B <--> E(B, Q)                                           :: OnlyA
          end
      end
      p_q_tagged = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              (E + A ⇌ E(A), E(P) + A ⇌ E(A, P), E(B) + A ⇌ E(P, Q))          :: OnlyA
              (E + B <--> E(B), E(A) + B <--> E(P, Q))                        :: OnlyA
              (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q)) :: OnlyA
              (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(P) + Q <--> E(P, Q)) :: OnlyA
              E(Q) + B <--> E(B, Q)                                           :: EqualAI
          end
      end
      p_b_at_eq_tagged = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              (E + A ⇌ E(A), E(P) + A ⇌ E(A, P), E(B) + A ⇌ E(P, Q))          :: OnlyA
              (E + B <--> E(B), E(A) + B <--> E(P, Q))                        :: OnlyA
              (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q)) :: OnlyA
              (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(P) + Q <--> E(P, Q)) :: EqualAI
              E(Q) + B <--> E(B, Q)                                           :: OnlyA
          end
      end
      q_b_at_eq_tagged = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              (E + A ⇌ E(A), E(P) + A ⇌ E(A, P), E(B) + A ⇌ E(P, Q))          :: OnlyA
              (E + B <--> E(B), E(A) + B <--> E(P, Q))                        :: OnlyA
              (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q)) :: EqualAI
              (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(P) + Q <--> E(P, Q)) :: OnlyA
              E(Q) + B <--> E(B, Q)                                           :: OnlyA
          end
      end
      all_tagged = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              (E + A ⇌ E(A), E(P) + A ⇌ E(A, P), E(B) + A ⇌ E(P, Q))          :: OnlyA
              (E + B <--> E(B), E(A) + B <--> E(P, Q))                        :: OnlyA
              (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q)) :: OnlyA
              (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(P) + Q <--> E(P, Q)) :: OnlyA
              E(Q) + B <--> E(B, Q)                                           :: OnlyA
          end
      end
      p_tagged = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              (E + A ⇌ E(A), E(P) + A ⇌ E(A, P), E(B) + A ⇌ E(P, Q))          :: OnlyA
              (E + B <--> E(B), E(A) + B <--> E(P, Q))                        :: OnlyA
              (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q)) :: OnlyA
              (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(P) + Q <--> E(P, Q)) :: EqualAI
              E(Q) + B <--> E(B, Q)                                           :: EqualAI
          end
      end
      _testhelper_admission_fault(t) = EnzymeRates._admission(
          EnzymeRates.AllostericMechanism(m, t, 1, EnzymeRates.RegulatorySite[])).fault
      vt = [c ? :OnlyA : :EqualAI for c in EnzymeRates._chemistry_groups(m)]
      @test _testhelper_admission_fault(vt) === EnzymeRates._NO_LIMIT
      @test isempty(EnzymeRates._tagging_bases(m))
      kids = EnzymeRates._expand_to_allosteric(m, rxn)
      @test length(kids) == 6
      @test Set(kids) == Set([q_tagged, b_at_eq_tagged, p_q_tagged, p_b_at_eq_tagged,
                              q_b_at_eq_tagged, all_tagged])
      @test all(EnzymeRates._is_well_posed, kids)
      @test EnzymeRates._admission(p_tagged).fault === EnzymeRates._NO_LIMIT
  end
  ```

  (f) Add, after the testset of (e):

  ```julia
  @testset "a ter-bi ping-pong without a regulator can end with no child" begin
      # Census ter-bi parent 0dc399e7dab6e4f5 (census/base6), one of the 98 ter-bi
      # ping-pong parents left childless (P4-19a): its covalent form E(; residual = B + D -
      # Q) takes up A and gives off P. Four groups hold chemistry, those of B, D, A and Q,
      # each mixing a plain binding with a fused one; P's binding at E is the only binding
      # group. Tagging every chemistry group :OnlyA strands the covalent forms, an island,
      # so the bases are the smallest stopping sets admission accepts: B's and D's groups,
      # B's and Q's, D's and Q's. Each K-type adds P's binding to a base or to the
      # every-chemistry set, and admission drops all four: on a base a stranded basin is
      # pumped, on the every-chemistry set the island stays. The reaction declares no
      # regulator, so there is no V-type either: no child, by design. Each base passes the
      # guard, so a declared regulator would give each two V-types.
      m = EnzymeRates.Mechanism(@enzyme_mechanism begin
          substrates: A, B, D
          products: P, Q
          steps: begin
              (E + B <--> E(B), E(D) + B <--> E(Q; residual = B + D - Q))
              (E + D ⇌ E(D), E(B) + D ⇌ E(Q; residual = B + D - Q))
              E + P <--> E(P)
              (E(; residual = B + D - Q) + A ⇌ E(P),
               E(Q; residual = B + D - Q) + A ⇌ E(A, Q; residual = B + D - Q))
              (E(; residual = B + D - Q) + Q <--> E(Q; residual = B + D - Q),
               E(P) + Q <--> E(A, Q; residual = B + D - Q))
          end
      end)
      rxn = @enzyme_reaction begin
          substrates: A[C], B[N], D[X]
          products: P[CN], Q[X]
      end
      base_bd = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B, D ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              (E + B <--> E(B), E(D) + B <--> E(Q; residual = B + D - Q))      :: OnlyA
              (E + D ⇌ E(D), E(B) + D ⇌ E(Q; residual = B + D - Q))            :: OnlyA
              E + P <--> E(P)                                                  :: EqualAI
              (E(; residual = B + D - Q) + A ⇌ E(P),
               E(Q; residual = B + D - Q) + A ⇌ E(A, Q; residual = B + D - Q)) :: EqualAI
              (E(; residual = B + D - Q) + Q <--> E(Q; residual = B + D - Q),
               E(P) + Q <--> E(A, Q; residual = B + D - Q))                    :: EqualAI
          end
      end
      base_bq = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B, D ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              (E + B <--> E(B), E(D) + B <--> E(Q; residual = B + D - Q))      :: OnlyA
              (E + D ⇌ E(D), E(B) + D ⇌ E(Q; residual = B + D - Q))            :: EqualAI
              E + P <--> E(P)                                                  :: EqualAI
              (E(; residual = B + D - Q) + A ⇌ E(P),
               E(Q; residual = B + D - Q) + A ⇌ E(A, Q; residual = B + D - Q)) :: EqualAI
              (E(; residual = B + D - Q) + Q <--> E(Q; residual = B + D - Q),
               E(P) + Q <--> E(A, Q; residual = B + D - Q))                    :: OnlyA
          end
      end
      base_dq = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B, D ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              (E + B <--> E(B), E(D) + B <--> E(Q; residual = B + D - Q))      :: EqualAI
              (E + D ⇌ E(D), E(B) + D ⇌ E(Q; residual = B + D - Q))            :: OnlyA
              E + P <--> E(P)                                                  :: EqualAI
              (E(; residual = B + D - Q) + A ⇌ E(P),
               E(Q; residual = B + D - Q) + A ⇌ E(A, Q; residual = B + D - Q)) :: EqualAI
              (E(; residual = B + D - Q) + Q <--> E(Q; residual = B + D - Q),
               E(P) + Q <--> E(A, Q; residual = B + D - Q))                    :: OnlyA
          end
      end
      tags = EnzymeRates._tagging_bases(m)
      bases = [EnzymeRates.AllostericMechanism(m, t, 1, EnzymeRates.RegulatorySite[])
               for t in tags]
      @test length(bases) == 3
      @test Set(bases) == Set([base_bd, base_bq, base_dq])
      @test all(EnzymeRates._is_well_posed, bases)
      chem = EnzymeRates._chemistry_groups(m)
      p = only(findall(!, chem))
      vt = [c ? :OnlyA : :EqualAI for c in chem]
      with_p(t) = (u = copy(t); u[p] = :OnlyA; u)
      _testhelper_admission_fault(t) = EnzymeRates._admission(
          EnzymeRates.AllostericMechanism(m, t, 1, EnzymeRates.RegulatorySite[])).fault
      @test _testhelper_admission_fault(vt) === EnzymeRates._ISLAND
      @test all(t -> _testhelper_admission_fault(with_p(t)) === EnzymeRates._PUMPED, tags)
      @test _testhelper_admission_fault(with_p(vt)) === EnzymeRates._ISLAND
      @test isempty(EnzymeRates._expand_to_allosteric(m, rxn))
  end
  ```

  (g) Add, after the testset of (f), the test of change-tag's filter: today's filter keeps a
  relaxation that admission refuses.

  ```julia
  @testset "change_allo_state keeps only the relaxations admission accepts" begin
      # A K-type of census bi-bi parent 1291abb91f537c6b (testset "K-types build on
      # today's V-type set when admission rejects it"): A's and B's mixed groups and Q's
      # bindings :OnlyA. Relaxing Q's bindings to :NonequalAI gives back the every-chemistry
      # set, whose plain binding of B at E has its two sides joined in the inactive
      # conformation, through Q's bindings and B's binding at E(Q) (`_NO_LIMIT`), so the
      # move drops it. Relaxing P's bindings or B's binding at E(Q) leaves Q's cut, and
      # admission accepts both. The chemistry does not relax while a binding is :OnlyA.
      # Two children.
      q_tagged = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              (E + A ⇌ E(A), E(P) + A ⇌ E(A, P), E(B) + A ⇌ E(P, Q))          :: OnlyA
              (E + B <--> E(B), E(A) + B <--> E(P, Q))                        :: OnlyA
              (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q)) :: EqualAI
              (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(P) + Q <--> E(P, Q)) :: OnlyA
              E(Q) + B <--> E(B, Q)                                           :: EqualAI
          end
      end
      p_relaxed = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              (E + A ⇌ E(A), E(P) + A ⇌ E(A, P), E(B) + A ⇌ E(P, Q))          :: OnlyA
              (E + B <--> E(B), E(A) + B <--> E(P, Q))                        :: OnlyA
              (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q)) :: NonequalAI
              (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(P) + Q <--> E(P, Q)) :: OnlyA
              E(Q) + B <--> E(B, Q)                                           :: EqualAI
          end
      end
      b_at_eq_relaxed = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              (E + A ⇌ E(A), E(P) + A ⇌ E(A, P), E(B) + A ⇌ E(P, Q))          :: OnlyA
              (E + B <--> E(B), E(A) + B <--> E(P, Q))                        :: OnlyA
              (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q)) :: EqualAI
              (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(P) + Q <--> E(P, Q)) :: OnlyA
              E(Q) + B <--> E(B, Q)                                           :: NonequalAI
          end
      end
      q_relaxed = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
          steps: begin
              (E + A ⇌ E(A), E(P) + A ⇌ E(A, P), E(B) + A ⇌ E(P, Q))          :: OnlyA
              (E + B <--> E(B), E(A) + B <--> E(P, Q))                        :: OnlyA
              (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q)) :: EqualAI
              (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(P) + Q <--> E(P, Q)) :: NonequalAI
              E(Q) + B <--> E(B, Q)                                           :: EqualAI
          end
      end
      kids = EnzymeRates._expand_change_allo_state(q_tagged)
      @test length(kids) == 2
      @test Set(kids) == Set([p_relaxed, b_at_eq_relaxed])
      @test EnzymeRates._admission(q_relaxed).fault === EnzymeRates._NO_LIMIT
      @test all(EnzymeRates._is_well_posed, kids)
  end
  ```
  (Plan probe `$S/phase45/plan/probes/p4/tagging5.log`: the research rules admit the P and
  B-at-E(Q) relaxations, both pass the research guard, and the Q relaxation fails existence.
  Step 8's value check runs both oracles on the two children.)

  (h) Comments that name the deleted test (`_onlya_offenders`) as the moves' filter. Each
  `old` is phase 3's text, each found at 08c4043; skip one that Task 8 already rewrote.

  | Testset | Old | New |
  |---|---|---|
  | "AllostericMechanism — no relaxation makes a kept copy redundant" | `# Relaxing the chemistry while S binds the active state only is rejected by` / `` # `_onlya_offenders` first, since the Haldane row then has one `:OnlyA` `` / `# column, and partial catalysis would reject it too: the inactive state would` / `# catalyze in part. The P binding and the copy relax. The relaxed copy, with one` | `# The move never relaxes the chemistry while S binds the active state only, and` / `# admission would refuse it: the chemistry and the P binding would join E to E(S)` / `# in the inactive state, leaving S's affinity no path to zero. The P binding and` / `# the copy relax. The relaxed copy, with one` |
  | "Mechanism — non-hyperbolic catalytic scheme: multiplicity 1 only" | `# the chemistry :OnlyA: 2^3 - 1 = 7 K-type children, none rejected by` / ``# `_onlya_offenders`.`` | `# the chemistry :OnlyA: 2^3 - 1 = 7 K-type children, all admitted (`_admits`).` |
  | "Mechanism — substrate declared as an inhibitor: promotion proceeds" | `# 2^5 - 1 = 31 candidate K-type children; _onlya_offenders drops` / `# none of them. The reaction declares no allosteric regulator, so no` | `# 2^5 - 1 = 31 candidate K-type children; admission drops` / `# none of them. The reaction declares no allosteric regulator, so no` |
  | "_expand_change_allo_state on a merged uni-uni K-type" | `# P's binding is :OnlyA: that would restore inactive catalysis beside it, a partial` / `# :OnlyA catalysis whose one-sided :OnlyA binding also leaves the Haldane cycle` / `# unsatisfiable once the chemistry is back in the check graph` / ``# (`_onlya_offenders`). One child.`` | `# P's binding is :OnlyA: restored, it would join E to E(P) in the inactive` / `# conformation beside the :OnlyA P binding, which then has no path to its limit` / `# (admission's existence test). One child.` |
  | "_expand_to_allosteric on an ordered Theorell–Chance mechanism" | ``# {A}, {Q} and {A, Q}. `_onlya_offenders` drops the :OnlyA chemistry group`` / `# from its cycle graph; the A and Q bindings left form the tree E(A) – E – E(Q),` / `# so no cycle row remains and no subset is refused. No regulator, so no V-type:` | `# {A}, {Q} and {A, Q}. In the inactive conformation the :OnlyA chemistry leaves` / `# the A and Q bindings as the tree E(A) – E – E(Q), so every subset's bindings` / `# point away from E and admission refuses none. No regulator, so no V-type:` |
  | "_expand_to_allosteric on a merged ordered bi-bi at multiplicity 2" | ``# `_onlya_offenders` refuses none. The scheme is hyperbolic, so`` | `# admission refuses none. The scheme is hyperbolic, so` |
  | "_expand_to_allosteric on a merged bi-bi whose B group mixes fused and plain" | `# and Q groups: 2^3 − 1 = 7. Without the :OnlyA B group the check graph keeps the` / `# tree E(A) – E – E(Q) – E(P, Q) (E(B, Q) leaves with B's group), so` / ``# `_onlya_offenders` refuses none. No regulator, so no V-type: seven`` | `# and Q groups: 2^3 − 1 = 7. Without the :OnlyA B group the inactive conformation` / `# keeps the tree E(A) – E – E(Q) – E(P, Q), and E(B, Q) lies behind B's plain` / `# binding at E(Q), so admission refuses none. No regulator, so no V-type: seven` |
  | "_expand_to_allosteric on a random-order merged bi-bi led by plain bindings" | `# the :OnlyA A and B groups the check graph keeps the path E – E(Q) – E(P, Q), so` / ``# `_onlya_offenders` refuses none. The reaction allows one catalytic subunit`` | `# the :OnlyA A and B groups the inactive conformation keeps the path E – E(Q) –` / `# E(P, Q), and E(A) and E(B) lie behind A's and B's bindings at E, so admission` / `# refuses none. The reaction allows one catalytic subunit` |
  | "to_allosteric shares :EqualAI downstream constants across conformations" | ``# `_expand_to_allosteric` emits dead-inactive `:OnlyA`-binding combos (all`` / `# chemical steps `:OnlyA`), so the case this guards is the child that sets` | ``# Among its K-types `_expand_to_allosteric` always emits dead-inactive`` / `# `:OnlyA`-binding combos (all chemical steps `:OnlyA`), so the case this guards is` / `# the child that sets` |
  | "_expand_change_allo_state drops the chemical step's relaxation" (and its section line) | `# ─── _expand_change_allo_state Haldane filter ──────────────────────────` ; `# A regulatory ligand's tag is not an argument to the Haldane check —` | `# ─── _expand_change_allo_state admission filter ────────────────────────` ; `# A regulatory ligand's tag does not enter admission —` |
  | "change_allo_state relaxes the bindings of a hand-built partial parent" | `# enumeration never builds one, since every mechanism it emits is fully dead or fully` / `# live in the inactive conformation; this one is written by hand. Each binding` / `# relaxation leaves three :OnlyA bindings of both signs on the catalytic cycle, so` / `# the Haldane check accepts all four. The chemistry relaxation is not tried while a` | `# enumeration builds such a parent only as a K-type on a ping-pong's half-reaction` / `# base; this one is written by hand. Each binding relaxation leaves three :OnlyA` / `# bindings whose arrows still point away from E and close no cycle, so admission` / `# accepts all four. The chemistry relaxation is not tried while a` |

  Each `/` separates lines of one comment; keep each line under 92 columns. No code
  changes.

- [ ] **Step 2: Write the failing tests in `test/test_types.jl` and `test/test_oracles.jl`.**

  (a) In the testset "rational nullspace + Stiemke feasibility helpers", delete the inner
  testset `@testset "_has_strict_positive_combination" begin … end` (its six assertions and
  their comments) and rename the outer testset to `"the rational nullspace helper"`; the
  inner testset `"_rational_nullspace"` and its assertions stay (P4-29a: the helper goes,
  so its test goes).

  (b) In the testset "OnlyA Haldane validator", Task 8 made three helpers
  (`_testhelper_uni`, `_testhelper_both_bindings_onlya`, and `_testhelper_verdict` in
  "random-order binding square") run `ER._onlya_offenders` and assert that admission agrees.
  Their verdicts become admission's. Rename the testset "admission rejects a one-sided :OnlyA
  binding and admits balanced ones", put this helper first in its body:

  ```julia
      # Admission's fault for `am`'s graph tagged `tags`, and the metabolite of the
      # binding it names (`nothing` when it names none).
      function _testhelper_admission_reads(am, tags)
          cm = am.catalytic
          v = ER._admission(ER.AllostericMechanism(cm, tags, 1, ER.RegulatorySite[]))
          e = v.edge == 0 ? nothing : cm.edges[v.edge]
          (v.fault, e === nothing || !ER._is_binding(e) ? nothing :
                    cm.table[ER._bound_entry(e)].name)
      end
      NO, NL = ER._NO_FAULT, ER._NO_LIMIT
  ```
  and make these replacements (each helper's old text is Task 8's, from its Step 3(c);
  each assertion's is phase 3's, as at 08c4043, which Task 8 left):

  - in `_testhelper_uni`, replace Task 8's closing statements
    ```julia
            offenders = ER._onlya_offenders(am.catalytic, am.catalytic.group, tags)
            # Admission (`_admission`) gives the same verdict: a tag set the moves' test
            # rejects has no limit (`_NO_LIMIT`).
            v = ER._admission(ER.AllostericMechanism(am.catalytic, tags, 2,
                                                     ER.RegulatorySite[]))
            @test v.fault === (isempty(offenders) ? ER._NO_FAULT : ER._NO_LIMIT)
            [only(firsts[g].consumed) for g in offenders]
    ```
    with
    ```julia
            _testhelper_admission_reads(am, tags)
    ```
    and the testset's opening comment
    ```julia
        # Uni-uni S -> P. Tags: (S binding, chemical step, P binding). The names of the
        # metabolites the offending :OnlyA groups bind.
    ```
    with
    ```julia
        # Uni-uni S -> P. Tags: (S binding, chemical step, P binding). Admission's fault,
        # and the metabolite of the binding it names.
    ```
  - replace the seven uni-uni assertions (from `# no :OnlyA anywhere -> valid` to the
    `:NonequalAI catalysis` one) with
    ```julia
        # no :OnlyA anywhere -> admitted
        @test _testhelper_uni(:EqualAI, :EqualAI, :EqualAI) == (NO, nothing)
        # :OnlyA on the substrate only, catalysis :EqualAI: catalysis and the P binding
        # join E to E(S) in the inactive conformation -> no limit, by the S binding
        @test _testhelper_uni(:OnlyA, :EqualAI, :EqualAI) == (NL, :S)
        # :OnlyA on the product only, catalysis :EqualAI -> no limit, by the P binding
        @test _testhelper_uni(:EqualAI, :EqualAI, :OnlyA) == (NL, :P)
        # :OnlyA on the substrate and the catalysis: E(S) lies behind the binding
        @test _testhelper_uni(:OnlyA, :OnlyA, :EqualAI) == (NO, nothing)
        # balanced: :OnlyA on both bindings, catalysis :EqualAI -> admitted
        @test _testhelper_uni(:OnlyA, :EqualAI, :OnlyA) == (NO, nothing)
        # V-system: :OnlyA chemical step only -> admitted
        @test _testhelper_uni(:EqualAI, :OnlyA, :EqualAI) == (NO, nothing)
        # :NonequalAI catalysis joins E to E(S) too -> same verdict
        @test _testhelper_uni(:OnlyA, :NonequalAI, :EqualAI) == (NL, :S)
    ```
  - in `_testhelper_both_bindings_onlya`, replace Task 8's closing statements
    ```julia
            offenders = ER._onlya_offenders(am.catalytic, am.catalytic.group, tags)
            v = ER._admission(ER.AllostericMechanism(am.catalytic, tags, 2,
                                                     ER.RegulatorySite[]))
            @test v.fault === (isempty(offenders) ? ER._NO_FAULT : ER._NO_LIMIT)
            offenders
    ```
    with
    ```julia
            _testhelper_admission_reads(am, tags)
    ```
    and its two assertions
    `@test isempty(_testhelper_both_bindings_onlya(mixed_uni_re_ss))` and
    `@test isempty(_testhelper_both_bindings_onlya(mixed_uni_ss_re))` with
    ```julia
        @test _testhelper_both_bindings_onlya(mixed_uni_re_ss) == (NO, nothing)
        @test _testhelper_both_bindings_onlya(mixed_uni_ss_re) == (NO, nothing)
    ```
    and the comment above `mixed_uni_re_ss` ("Whether a binding is rapid-equilibrium or
    steady-state does not change which affinities diverge, …") with
    ```julia
        # Whether a binding is rapid-equilibrium or steady-state does not change
        # admission's verdict: with both bindings :OnlyA, E(S) and E(P) lie behind them.
        # Both orientations of the mix are admitted.
    ```
  - in "random-order binding square", replace Task 8's
    ```julia
            function _testhelper_verdict(t)
                offenders = ER._onlya_offenders(bu_am.catalytic, bu_am.catalytic.group, t)
                v = ER._admission(ER.AllostericMechanism(bu_am.catalytic, t, 2,
                                                         ER.RegulatorySite[]))
                @test v.fault === (isempty(offenders) ? ER._NO_FAULT : ER._NO_LIMIT)
                offenders
            end
    ```
    with
    ```julia
            _testhelper_verdict(t) = _testhelper_admission_reads(bu_am, t)
    ```
    the assertion pair
    `@test isempty(_testhelper_verdict(_testhelper_tags()))` and
    `@test _testhelper_verdict(_testhelper_tags((:E, :A))) == [_testhelper_group(:E, :A)]`
    with
    ```julia
            # the square alone trips nothing: no :OnlyA anywhere -> admitted
            @test _testhelper_verdict(_testhelper_tags()) == (NO, nothing)
            # one square edge :OnlyA: the square's other three bindings join its sides
            @test _testhelper_verdict(_testhelper_tags((:E, :A))) == (NL, :A)
    ```
    (the `lone` block that follows, which replaced `lone_onlya` in Task 8, stays), and the
    three assertions after it with
    ```julia
            # both A-side bindings :OnlyA: the B bindings, the catalysis and the P binding
            # still join E to E(A) -> no limit
            @test _testhelper_verdict(_testhelper_tags((:E, :A), (:EB, :A)))[1] === NL
            # ...and :OnlyA on P cuts the last join -> admitted
            @test _testhelper_verdict(
                _testhelper_tags((:E, :A), (:EB, :A), (:E, :P))) == (NO, nothing)
            # :OnlyA on the chemical step leaves the square, whose lone :OnlyA edge still
            # has its sides joined -> no limit, by that edge
            @test _testhelper_verdict(_testhelper_tags((:EAB, nothing), (:E, :A))) ==
                  (NL, :A)
    ```
    Delete the helper `_testhelper_group`, whose last reader was the old assertion, and
    the comments that read the square through "Wegscheider row", "Haldane row" or
    "VIOLATION": the block comment above `biuni` becomes
    ```julia
        # A random-order binding square: E binds A and B, and each of E(A) and E(B) binds
        # the other. A lone :OnlyA edge has both sides joined by the other three, so
        # admission finds no limit; tagging P's binding too separates them.
    ```

  (c) In Task 8's testset "a ping-pong with both fused releases :OnlyA is an island; one half
  passes" (which replaced ":OnlyA guard admits a ping-pong whose fused releases are
  :OnlyA"), delete
  ```julia
      @test isempty(ER._onlya_offenders(both.catalytic, both.catalytic.group,
                                        ER.cat_allo_states(both)))
  ```
  and in its first comment replace the lines
  ```julia
      # (`_is_chemistry`), though each gives off one metabolite and takes up none. Tagged
      # :OnlyA, both leave the moves' :OnlyA test, which reads no cycle through chemistry.
  ```
  with
  ```julia
      # (`_is_chemistry`), though each gives off one metabolite and takes up none.
  ```
  Task 8's admission and oracle assertions stay.

  (d) The three testsets that follow describe the deleted sign and Stiemke stages, and Task 8
  gave each an `_onlya_offenders` assertion. In ":OnlyA guard rejects a multi-cycle
  ter-substrate inconsistency" delete Task 8's
  ```julia
      @test !isempty(ER._onlya_offenders(cube.catalytic, cube.catalytic.group,
                                         ER.cat_allo_states(cube)))
  ```
  and replace Task 8's comment
  ```julia
      # Admission finds untagged steps that join the two sides of an :OnlyA binding
      # (`_NO_LIMIT`), as the moves' Stiemke test finds the coupled system infeasible.
  ```
  with
  ```julia
      # Admission finds untagged steps that join the two sides of an :OnlyA binding
      # (`_NO_LIMIT`).
  ```
  In ":OnlyA guard admits feasible ter-substrate mechanisms" delete Task 8's
  ```julia
      # The moves' :OnlyA test decides it in its Stiemke stage (feasible, w = 1).
      cube_m = ER.AllostericMechanism(cube)
      @test isempty(ER._onlya_offenders(cube_m.catalytic, cube_m.catalytic.group,
                                        ER.cat_allo_states(cube_m)))
  ```
  In ":OnlyA guard reads cycles past an abortive segment it strands" delete Task 8's
  ```julia
      am_m = ER.AllostericMechanism(am)
      @test isempty(ER._onlya_offenders(am_m.catalytic, am_m.catalytic.group,
                                        ER.cat_allo_states(am_m)))
  ```
  and the last sentence of Task 8's comment above it ("The moves' :OnlyA test, which reads no
  cycle through the chemistry group, agrees."), and replace the testset's first comment
  (phase 3's, from "The A group mixes the fused binding" to "lies on none of them.") with
  ```julia
      # The A group mixes the fused binding into E(P, Q) with A's abortive bindings at
      # E(P) and E(Q), so it is chemistry. Tagged :OnlyA, it strands
      # {E(A, P), E(A, Q), E(A, P, Q)} in the inactive conformation: a rapid-equilibrium
      # segment whose every form holds a product, which a mechanism may not have.
  ```
  Then replace the comment block above the second testset, from "# The rejection testset
  above pins" to "(M is 5x12, nullity 7, feasible = true).", whichever wording Task 8 gave
  its first sentence, with
  ```julia
  # The rejection testset above pins admission's reject direction on a ter-substrate
  # binding cube; this one pins its accept direction there. An over-rejection does not
  # fail loudly: it silently shrinks the enumeration's search space by dropping valid
  # mechanisms before they are fitted.
  ```
  its inner comment on `cube` with
  ```julia
      # A full random-order A/B/C binding lattice, every substrate binding and the
      # chemical step :OnlyA: in the inactive conformation every form stands alone but E
      # and E(P), the bindings point away from free E and close no cycle, so admission
      # accepts the cube.
  ```
  the comment on `ordered` with
  ```julia
      # Ordered ter with an :OnlyA B binding and an :OnlyA chemical step: E(A, B) and
      # E(A, B, C) lie behind the B binding, on its bound side, so this is admitted.
  ```
  and the first testset's comment (phase 3's, from "A full random-order ter binding cube"
  to "perpetual motion.") with
  ```julia
      # A full random-order ter binding cube with seven :OnlyA edges and the cross
      # bindings :EqualAI. In the inactive conformation the untagged bindings join forms
      # across the cube, so a tagged binding's two sides share a component: no path takes
      # its affinity to zero, and admission rejects the cube (`_NO_LIMIT`).
  ```

  (e) Append to `test/test_oracles.jl`:

  ```julia
  # The census ping-pong declaration (census/README.md, "Reactions"): the
  # ping-pong-capable bi-bi reaction with the allosteric regulator R, two catalytic
  # subunits and A and P on one shared site.
  const _testhelper_pingpong_census_rxn = @enzyme_reaction begin
      substrates: A[CX], B[N]
      products: P[C], Q[NX]
      allosteric_regulators: R
      oligomeric_state: 2
      shared_catalytic_site: (A, P)
  end

  @testset "ping-pong seeds' tagging children are admitted, distinct and oracle-checked" begin
      # Research section 3.3, P6, an aggregate pin over the seed set: the 18 ping-pong seeds
      # (seeds with a covalent form) of the census ping-pong reaction. Each seed's bases
      # stop the inactive conformation's catalysis, pass admission, and tag disjoint sets
      # of chemistry groups, one half-reaction each: 14 seeds get two bases and 4 get one
      # (the research rules on census pingpong level 0). Every child is admitted and
      # distinct; on every face its guard verdict agrees with the RE-speed oracle, and a
      # child the guard admits converges along the ε-path at the interior. Research section
      # 3.3 asks for the ε-path on every face; the suite's time budget (P4-28) keeps it to
      # the interior. The children number CHILDREN, as a separate implementation of the
      # rules counts them.
      rxn = _testhelper_pingpong_census_rxn
      rng = Random.MersenneTwister(406)
      seeds = [m for m in EnzymeRates.init_mechanisms(rxn)
               if any(EnzymeRates._has_residual, m.forms)]
      @test length(seeds) == 18
      nbases, n_kids = Dict{Int, Int}(), 0
      for m in seeds
          tags = EnzymeRates._tagging_bases(m)
          nbases[length(tags)] = get(nbases, length(tags), 0) + 1
          onlya = [Set(findall(==(:OnlyA), t)) for t in tags]
          @test all(isdisjoint(onlya[i], onlya[j]) for i in eachindex(onlya)
                    for j in (i + 1):length(onlya))
          @test all(t -> EnzymeRates._admits(m, t) &&
                         EnzymeRates._stops_inactive_turnover(m, t), tags)
          kids = EnzymeRates._expand_to_allosteric(m, rxn)
          n_kids += length(kids)
          @test allunique(kids)
          for c in kids
              @test EnzymeRates._admits(c.catalytic, c.cat_allo_states)
              faults = _testhelper_face_faults(c)
              admitted = all(==(EnzymeRates._NO_FAULT), faults)
              verdicts = _testhelper_respeed_verdicts(c; rng,
                                                      law = admitted ? :package : :limit)
              @test all(_testhelper_agrees.(faults, verdicts))
              admitted && @test _testhelper_epath_verdict(c; rng) === :converges
          end
      end
      @test nbases == Dict(1 => 4, 2 => 14)
      @test n_kids == CHILDREN
  end
  ```
  `CHILDREN` is a number: before writing the testset, count the children Task 11's
  candidate stage recorded for the census pingpong level-0 ping-pong parents (the same 18
  seeds):

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  awk -F'\t' 'FNR == NR { if ($1 == "pingpong" && $3 == 0 && $5 == "true") s[$4] = 1; next }
              $1 == "pingpong" && ($2 in s) { n++ } END { print n }' \
      $S/phase45/pumped/data/parents.tsv $S/phase45/pumped/data/children.tsv
  ```
  and write that number for `CHILDREN` in the comment and the last assertion. (Task 11's
  parents file has 15 columns: reaction, index, level, id, ping-pong, …; its children file
  has the parent id in column 2.)

- [ ] **Step 3: Run the tests to verify they fail.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  $S/phase3/jl.sh $S/phase45/logs/P4_T12_red.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_mechanism_enumeration.jl")'
  $S/phase3/jl.sh $S/phase45/logs/P4_T12_red_oracles.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_oracles.jl")'
  ```
  Expected in `P4_T12_red.log`: "Mechanism — Bi-bi ping-pong: binding-group :OnlyA variants
  only" fails `length(result) == 6` (Evaluated: 15 == 6, today's K-types on both halves);
  "the tagging move gives a ping-pong one V-type base per half", "K-types build on today's
  V-type set when admission rejects it" and "a ter-bi ping-pong without a regulator can end
  with no child" error with `UndefVarError: `_tagging_bases` not defined in `EnzymeRates``;
  "change_allo_state keeps only the relaxations admission accepts" fails `length(kids) == 2`
  (Evaluated: 3 == 2: today's filter keeps the Q relaxation) and `all(_is_well_posed,
  kids)`. The two change-tag testsets (c) and (d) pass already: today's filter keeps the same
  relaxations of these parents, and Task 8's exit drops the class-3 failure; they are
  conversions of their subjects to P4-1's parents. In `P4_T12_red_oracles.log` the P6
  testset errors with the same `UndefVarError`. `test_types.jl`'s converted assertions read
  `_admission`, which exists, so they pass before the deletion of Step 6.

- [ ] **Step 4: Implement the move.** In `src/mechanism_enumeration.jl`, insert
  before `_expand_to_allosteric`'s docstring:

  ```julia
  """
      _admits(m::Mechanism, tags::Vector{Symbol}) -> Bool

  Whether admission (`_admission`) accepts `m` with the tags `tags`, one per kinetic group: a
  path to the limit of the tagged constants exists, every form free E reaches only through
  tagged steps lies behind a tagged plain binding, the active conformation turns over, and
  the face condition holds on every face. It judges one catalytic subunit with no regulatory
  site, which admission never reads. Class 3 is not run: the guard at the enumerator's exits
  drops a child that fails it (research section 1.9, `ok`).
  """
  _admits(m::Mechanism, tags::Vector{Symbol}) =
      _admission(AllostericMechanism(m, tags, 1, RegulatorySite[])).fault === _NO_FAULT

  """
      _stops_inactive_turnover(m::Mechanism, tags::Vector{Symbol}) -> Bool

  Whether the tags `tags` leave the inactive conformation of `m` no catalytic cycle: no
  cycle of the edges its law keeps (`_law_mask`) turns the reaction over
  (`_has_turnover_cycle`).
  """
  _stops_inactive_turnover(m::Mechanism, tags::Vector{Symbol}) =
      !_has_turnover_cycle(m, _law_mask(AllostericMechanism(m, tags, 1, RegulatorySite[]), 2))

  """
      _tagging_bases(m::Mechanism) -> Vector{Vector{Symbol}}

  The tag sets the tagging move (`_expand_to_allosteric`) builds on, one tag per kinetic
  group of `m` (spec section 7, "Allosteric ping-pong"). The first candidate tags every
  chemistry group (`_chemistry_groups`) `:OnlyA` and every binding group `:EqualAI`; when
  admission accepts it (`_admits`) it is the only base. Otherwise the bases are the
  smallest sets of chemistry groups, by inclusion, that stop every catalytic cycle of the
  inactive conformation (`_stops_inactive_turnover`) and pass admission, each tagged
  `:OnlyA` with every other group `:EqualAI`. A ping-pong mechanism thus gets one base per
  half-reaction, the other half running in both conformations. The bases come by size,
  each size in `_subsets` order.
  """
  function _tagging_bases(m::Mechanism)
      chem = _chemistry_groups(m)
      tagged(onlya) = [g in onlya ? :OnlyA : :EqualAI for g in eachindex(chem)]
      groups = findall(chem)
      _admits(m, tagged(groups)) && return [tagged(groups)]
      sets = _subsets(groups)
      found = Vector{Int}[]
      for k in 1:length(groups), s in sets
          length(s) == k && !any(b -> b ⊆ s, found) || continue
          t = tagged(s)
          _stops_inactive_turnover(m, t) && _admits(m, t) && push!(found, s)
      end
      [tagged(s) for s in found]
  end
  ```

  Replace `_expand_to_allosteric`'s docstring and body with:

  ```julia
  """
      _expand_to_allosteric(m::Mechanism, rxn::EnzymeReaction)
          → Vector{AllostericMechanism}

  Convert a non-allosteric `Mechanism` into allosteric variants, each built on a base of
  `_tagging_bases`: the tag set with every chemistry group `:OnlyA` and every binding group
  `:EqualAI` when admission accepts it, else each smallest set of chemistry groups that stops
  the inactive conformation's catalysis and passes admission. A chemistry group is one
  holding a chemistry step (`_is_chemistry`): an isomerization, a fused binding or a
  Theorell–Chance step; every other group is a binding group. On a ping-pong mechanism,
  tagging both half-reactions strands the covalent forms between them, so each half alone
  is a base and the other half runs in both conformations (spec section 7, "Allosteric
  ping-pong"). Two kinds whose conformational constant `L` never shows are not enumerated:
  the all-`:EqualAI` baseline and a V-type variant without a regulator.

    * The all-`:EqualAI` baseline is never emitted — the two conformations
      are identical, `L` cancels, and the mechanism is indistinguishable
      from `m`.
    * K-type: a nonempty set of binding groups set `:OnlyA` on top of a base, and on top
      of the every-chemistry set when admission rejects that set, since a tagged binding
      can restore its limit; each is emitted bare and kept when admission accepts it
      (`_admits`). Over more than one catalytic subunit the conformational equilibrium
      enters the rate with the multiplicity's power. Over one subunit the inactive
      conformation is a dead-end branch of the free enzyme, and `L` is a phantom in some
      K-type variants and not in others. In every K-type variant of the rapid-equilibrium
      ordered bi-bi seed a rescaling of the active state's constants absorbs it, whatever
      that conformation binds. In the merged ordered bi-bi whose B group mixes a fused and
      a plain binding it is a phantom only in the variants whose inactive conformation
      binds nothing: no rescaling absorbs an A or Q binding kept there. Binding nothing is
      not enough: in a random-order merged bi-bi, whose A group takes A up at E and at
      E(B) under one rate constant, `L` shows in every K-type variant.
    * V-type: a base, with no `:OnlyA` binding, paired with a declared allosteric
      regulator at a new site, one variant per `(regulator, tag)` with
      `tag ∈ {:OnlyA, :OnlyI}`. On the every-chemistry base the inactive state binds
      substrate and product as the active state does and runs no chemistry, so without a
      regulator `L` folds entirely into `kcat` (`v = kcat/(1+L)·shape`) and is not
      observable; every base is emitted only with a regulator. A reaction with no
      declared allosteric regulators emits no V-type.

  Admission does not run class 3 of the guard: a child that fails a face is dropped where
  `expand_mechanisms` returns it. A ping-pong parent with no declared regulator can end
  with no child, when every K-type strands a form or a basin.

  For each value in `rxn`'s `allowed_catalytic_multiplicities`, the move emits
  the K-types, then the V-types, at that `catalytic_multiplicity`. Every variant
  composes `m`'s graph, which it shares, with its tags, multiplicity and sites. The
  variants are distinct by construction: each differs from the others in its tags, its
  regulatory site or its multiplicity.

  A parent whose catalytic scheme fails `_hyperbolic_catalysis` (random-order
  steady-state binding, or a substrate that traps a steady-state intermediate in
  an abortive complex, whose own equation carries concentration powers) emits
  variants at multiplicity 1 only: above one subunit the conformational
  equilibrium would add a second source of powers
  (`_requires_hyperbolic_catalysis`).
  """
  function _expand_to_allosteric(m::Mechanism, rxn::EnzymeReaction)
      hyperbolic = _hyperbolic_catalysis(m)
      cns = [cn for cn in allowed_catalytic_multiplicities(rxn) if cn == 1 || hyperbolic]
      isempty(cns) && return AllostericMechanism[]
      chem = _chemistry_groups(m)
      vtags = [c ? :OnlyA : :EqualAI for c in chem]
      regs = [name(regulator(rm)) for rm in regulators(rxn)
              if regulator(rm) isa AllostericRegulator]
      bases = _tagging_bases(m)
      # K-type: a nonempty set of binding groups :OnlyA on top of a base, and on top of the
      # every-chemistry set when admission rejects it; admission drops the rest.
      kbases = vtags in bases ? bases : [bases; [vtags]]
      ktags = Vector{Symbol}[]
      for base in kbases, sel in _subsets(findall(!, chem))[2:end]
          tags = copy(base)
          tags[sel] .= :OnlyA
          _admits(m, tags) && push!(ktags, tags)
      end
      results = AllostericMechanism[]
      for cn in cns
          for tags in ktags
              push!(results, AllostericMechanism(m, tags, cn, RegulatorySite[]))
          end
          # V-type: a base paired with a declared regulator, whose binding reveals L.
          any(chem) || continue
          for base in bases, reg in regs, tag in (:OnlyA, :OnlyI)
              push!(results,
                    AllostericMechanism(m, base, cn, [RegulatorySite([reg], cn, [tag])]))
          end
      end
      results
  end
  ```

  In `_expand_change_allo_state`, replace its docstring with:

  ```julia
  """
      _expand_change_allo_state(am::AllostericMechanism)
          → Vector{AllostericMechanism}

  Relax a "constrained" tag (`:EqualAI`, `:OnlyA`, `:OnlyI`) to `:NonequalAI`. A binding
  catalytic group and a regulatory ligand each relax individually — one variant per group
  not already `:NonequalAI`. The chemistry groups (those holding a chemistry step,
  `_is_chemistry`) relax **together**, one variant setting every non-`:NonequalAI`
  chemistry group to `:NonequalAI` at once, so a fully live inactive conformation is one
  move from every V-type, the tagging move's one-half bases included. The base catalytic
  steps, multiplicity, and untouched tags are preserved.

  The chemistry relaxation is tried only while no binding group is `:OnlyA`. Beside an
  `:OnlyA` binding the restored chemistry can join the binding's two sides in the inactive
  conformation and leave its affinity no path to zero; where the inactive conformation
  binds nothing (every binding `:OnlyA`), the skipped variant is rate-equivalent to the
  form emitted directly, so no observable hypothesis is lost.

  A catalytic relaxation is kept when admission accepts it (`_admits`): relaxing an
  `:OnlyA` binding joins its two sides, which can leave another tagged binding inside one
  component or close a cycle of tagged bindings. A regulatory ligand's tag does not enter
  admission — a regulator site completes no catalytic cycle — so that branch needs no
  filter. Admission does not run class 3 of the guard: a child that fails a face is dropped
  where `expand_mechanisms` returns it.

  A child with a redundant competitive-inhibitor group (`_redundant_copy_groups`)
  is dropped. The per-state test reads the tags this move changes: relaxing an
  `:OnlyA` binding to `:NonequalAI` brings its complex into the inactive state,
  where a copy that was new there can become its twin with a consistent gauge, and
  the copy's constant then shows in neither state.
  """
  ```
  and in its body replace the comment block
  ```julia
      # Binding catalytic groups relax individually. Chemistry groups relax together.
      # Inactive catalysis is all-or-nothing, so a fully-`:NonequalAI` catalytic inactive
      # conformation is unreachable by relaxing chemistry groups one at a time: each
      # mixed intermediate is partial. One variant sets every non-`:NonequalAI` chemistry
      # group to `:NonequalAI` at once, and only while no binding group is `:OnlyA`, which
      # would also leave the inactive catalysis partial.
  ```
  with
  ```julia
      # Binding catalytic groups relax individually; the chemistry groups relax together,
      # in one variant, and only while no binding group is :OnlyA.
  ```
  and the filter line
  ```julia
          isempty(_onlya_offenders(am.catalytic, am.catalytic.group, new_states)) || continue
  ```
  with
  ```julia
          _admits(am.catalytic, new_states) || continue
  ```

- [ ] **Step 5: Run the tests to verify they pass.** The two commands of Step 3, logs
  `P4_T12_green.log` and `P4_T12_green_oracles.log`. Expected: `exit 0`; every testset
  passes, the converted and new ones included, and the five tagging-count testsets RES §3.4
  lists ("V-type param counts + competitive-inhibitor exclusion" 7 and 9, "Mechanism —
  substrate declared as an inhibitor: promotion proceeds" 31, "… on an ordered
  Theorell–Chance mechanism" 3, "… on a merged ordered bi-bi at multiplicity 2" 7, "… on a
  merged bi-bi whose B group mixes fused and plain" 7, "… on a random-order merged bi-bi led
  by plain bindings" 3, and "Mechanism — uni-uni: binding-group :OnlyA variants only" 3,
  "Mechanism — enumerates all allowed multiplicities" 6, "Mechanism — Bi-bi sequential" 15,
  "Mechanism — non-hyperbolic catalytic scheme: multiplicity 1 only" 7) pass unchanged: the
  research rules give each the same children (plan probe `$S/phase45/plan/probes/p4/
  tagging.log`), since admission accepts their every-chemistry sets and every K-type
  `_onlya_offenders` accepted. A count that changes is not re-pinned before Step 8 has
  checked its children with the oracles; STOP and debug it if the oracles disagree. In
  `test/test_mechanism_enumeration.jl`, Task 8's testset "enumeration drops the children
  the guard rejects" now fails: its ping-pong parent's tagging children are P4-1's. Step 8 re-pins it
  after the oracles check them.

  Two aggregates of `test_mechanism_enumeration.jl` read the move's output:
  - Task 8's "AllostericMechanism — a split keeps its parent's admission verdict" draws its
    parents from `_expand_to_allosteric`, which now emits only tag sets admission accepts,
    so `(accepted, rejected) == (716, 4)` fails with `rejected == 0` and the reject
    direction goes untested. Build the parents from each split child's every-chemistry set
    and its K-types unfiltered: replace
    ```julia
          parents = unique(EnzymeRates.AllostericMechanism[
              a for m in EnzymeRates.init_mechanisms(bi_bi_rxn)
              for c in EnzymeRates._expand_split_kinetic_group(m)
              for a in EnzymeRates._expand_to_allosteric(c, bi_bi_rxn)])[1:10:end]
    ```
    with
    ```julia
          # Each split child's every-chemistry tag set with each nonempty set of its
          # binding groups :OnlyA, unfiltered: the tagging move keeps only the tag sets
          # admission accepts (P4-1), and the reject direction needs parents it refuses.
          function _testhelper_unfiltered_tags(c)
              chem = EnzymeRates._chemistry_groups(c)
              vt = [x ? :OnlyA : :EqualAI for x in chem]
              [(t = copy(vt); t[s] .= :OnlyA; t)
               for s in EnzymeRates._subsets(findall(!, chem))[2:end]]
          end
          parents = unique(EnzymeRates.AllostericMechanism[
              EnzymeRates.AllostericMechanism(c, t, 1, EnzymeRates.RegulatorySite[])
              for m in EnzymeRates.init_mechanisms(bi_bi_rxn)
              for c in EnzymeRates._expand_split_kinetic_group(m)
              for t in _testhelper_unfiltered_tags(c)])[1:10:end]
    ```
    and in its comment "each split child of a seed lifted to its allosteric variants, one
    parent in ten" with "each split child of a seed tagged as below, one parent in ten".
    Rerun the focused file: the testset reports the new `(accepted, rejected)`. Before
    pinning it, check every rejected parent against the ε-path oracle, as Task 8 did for
    its four, with `$S/phase45/probes/P4_T12_split_parents.jl`:
    ```julia
    # ABOUTME: Task 12's check of the split-verdict aggregate's rejected parents: the
    # ABOUTME: ε-path oracle agrees with admission's interior fault on each one.
    using Test, EnzymeRates, LinearAlgebra, Random
    include(joinpath(pwd(), "test", "mechanism_definitions_for_test_enzyme_derivation.jl"))
    const bi_bi_rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    function unfiltered_tags(c)
        chem = EnzymeRates._chemistry_groups(c)
        vt = [x ? :OnlyA : :EqualAI for x in chem]
        [(t = copy(vt); t[s] .= :OnlyA; t)
         for s in EnzymeRates._subsets(findall(!, chem))[2:end]]
    end
    parents = unique(EnzymeRates.AllostericMechanism[
        EnzymeRates.AllostericMechanism(c, t, 1, EnzymeRates.RegulatorySite[])
        for m in EnzymeRates.init_mechanisms(bi_bi_rxn)
        for c in EnzymeRates._expand_split_kinetic_group(m)
        for t in unfiltered_tags(c)])[1:10:end]
    rng = Random.MersenneTwister(4)
    rejected = [(am, EnzymeRates._admission_interior(am).fault) for am in parents
                if EnzymeRates._admission_interior(am).fault !== EnzymeRates._NO_FAULT]
    bad = count(r -> !_testhelper_admission_agrees(r[2], _testhelper_epath_verdict(r[1]; rng)),
                rejected)
    println("parents ", length(parents), ", rejected ", length(rejected), ", ε-path disagrees ",
            bad)
    ```
    ```bash
    S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T12_split_parents.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/probes/P4_T12_split_parents.jl")'
    ```
    Expected: `ε-path disagrees 0`, and `parents` and `rejected` equal the testset's
    counts. Then pin `(accepted, rejected)` and give the old and new counts in the commit
    body.
  - "every bi-bi init mechanism and child is what its raw input builds" pins
    `length(kids)` and the allosteric count over `expand_mechanisms` of bi_bi_rxn's 239
    seeds (2,831 and 1,833 at 08c4043, or as Task 8 re-pinned them). P4-1 adds the K-types
    that `_onlya_offenders` refused and admission accepts (Task 0(b): 176 bi-bi children
    at 37 parents over levels 0-2), so a count may rise. If it fails, the extra children
    are those of the level-0 bi-bi parents in Step 9's move check, which compares them with
    Task 11's independently built candidates; once Step 9 passes, re-pin both counts to
    the values the failure prints and say why in the commit body.

- [ ] **Step 6: Delete the old `:OnlyA` test.** In `src/thermo.jl`, delete
  `_has_strict_positive_combination`, `_onlya_offenders` and `_onlya_cone_offenders` with
  their docstrings (and `_onlya_haldane_violation` and `_onlya_violation_message` if Task 8
  left either), and rename the section line
  `# ─── Cycle-space helpers of the :OnlyA guard and the split move's counter ─────` to
  `# ─── Cycle-space helpers of the split move's counter ───────────────────────`.
  `_rref_partition`, `_rational_nullspace` and `_integer_nullspace` stay.

  In `src/types.jl`, replace `_chemistry_groups` and the last sentence of
  `_is_chemistry`'s docstring:

  ```julia
  """Whether each kinetic group of `m` holds a chemistry edge (`_is_chemistry`). The
  allosteric moves tag a chemistry group apart from a binding group."""
  _chemistry_groups(m::Mechanism) =
      _groups_with(k -> _is_chemistry(m, m.edges[k]), m.group, ngroups(m))
  ```
  (at 08c4043, `_chemistry_groups(m, group = m.group, n = maximum(group; init = 0))`
  takes the two arguments for `_onlya_offenders` alone and already builds the vector
  with `_groups_with`), and in `_is_chemistry`'s
  docstring replace the lines
  ```julia
  with several metabolites on a side. The joint solve's priority (`_priority_key`) and the
  `:OnlyA` test (`_onlya_offenders`) read it to tell catalysis from binding."""
  ```
  with
  ```julia
  with several metabolites on a side. The joint solve's priority (`_priority_key`) and the
  allosteric moves (`_chemistry_groups`) read it to tell catalysis from binding."""
  ```

  Then check that no code or test names the deleted functions:

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && grep -rn "_onlya_offenders\|_onlya_cone_offenders\|_has_strict_positive_combination\|_onlya_haldane_violation\|_onlya_violation_message\|Stiemke" src test
  ```
  Expected: no output. A hit in a comment is rewritten to say what the code does now (the
  table of Step 1(h) gives the enumeration test's). `ROADMAP.md` still names
  `_onlya_offenders` until Step 11, which reruns the grep over it and `docs/src`.

- [ ] **Step 7: Run the affected files.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  for f in test_types test_mechanism_enumeration test_oracles test_reaction_graph; do $S/phase3/jl.sh $S/phase45/logs/P4_T12_${f}.log 1800 --project=. -e "using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include(\"test/mechanism_definitions_for_test_enzyme_derivation.jl\"); include(\"test/$f.jl\")"; done
  grep -h "Test Summary\|^\(test_\|EnzymeRates\)" $S/phase45/logs/P4_T12_test_*.log | head
  ```
  Expected: `test_types.jl` and `test_oracles.jl` pass every testset; in
  `test_mechanism_enumeration.jl` every testset passes but two that later steps re-pin
  after their value checks: Task 8's "enumeration drops the children the guard rejects"
  (Step 8; its parent's tagging children are now P4-1's), and, if P4-1 changed its counts,
  "every bi-bi init mechanism and child is what its raw input builds" (Step 9). In
  `test_reaction_graph.jl` every testset passes but, if P4-1 changed their counts, the two
  that pin the census populations, which Step 10 re-pins: "a graph's canonical form,
  pinned as text, read back from its Sig and raw" (`length(population) == 5508 + 5505`,
  whose second term counts the children of census bi-bi's init mechanisms, among them
  the K-types P4-1 adds) and "a graph's DSL text reads back as the graph"
  (`length(allosteric) == 5 + <n3> + 5` as Task 8 re-pinned it, whose middle term counts
  census reaction 3's allosteric children: P4-1 builds the tagging children of its 18
  ping-pong seeds on P6's bases, whose old ones Task 8's exit dropped as islands). Any other
  failure is a defect: STOP and debug it. No warning. Compare the total times of
  `test_types.jl`, `test_mechanism_enumeration.jl` and `test_oracles.jl` with Task 11's
  full-suite log: together they grow by at most 3 s (P6 costs about 2 s; the
  deleted Stiemke helper test and the old 15-child pins give some back). Record
  `test_mechanism_enumeration.jl`'s time and its ter-ter testset's ("init_mechanisms on
  ter-ter within 150 s", a top-level testset with its own Test Summary line) for the
  commit body, beside R52's branch-point figures.

- [ ] **Step 8: Check the new children's values (before the commit).** Create
  `$S/phase45/probes/P4_T12_value_check.jl`:

  ```julia
  # ABOUTME: Task 12's value check: both test oracles on every child the tagging and
  # ABOUTME: change-tag moves give the fixtures whose pins Task 12 writes.
  using Test, EnzymeRates, LinearAlgebra, Random
  # A relative path would resolve against this file's directory, not the worktree's.
  include(joinpath(pwd(), "test", "mechanism_definitions_for_test_enzyme_derivation.jl"))
  const ER = EnzymeRates
  # Task 8's exit-test parent (TME), census ping-pong 3fdbe6e1dc001bb7.
  pp4 = ER.Mechanism(@enzyme_mechanism begin
      substrates: A, B
      products: P, Q
      steps: begin
          E + A ⇌ E(P; residual = A - P)
          E(; residual = A - P) + P <--> E(P; residual = A - P)
          E(; residual = A - P) + B <--> E(Q)
          E + Q <--> E(Q)
      end
  end)
  rxn_pp4 = @enzyme_reaction begin
      substrates: A[CX], B[N]
      products: P[C], Q[NX]
      allosteric_regulators: R
      oligomeric_state: 2
      shared_catalytic_site: (A, P)
  end
  pp = ER.Mechanism(@enzyme_mechanism begin
      substrates: A, B
      products: P, Q
      steps: begin
          E + A ⇌ E(A)
          E(; residual = A - P) + B ⇌ E(B; residual = A - P)
          E + Q ⇌ E(Q)
          E(; residual = A - P) + P ⇌ E(P; residual = A - P)
          E(A) <--> E(P; residual = A - P)
          E(B; residual = A - P) ⇌ E(Q)
      end
  end)
  tc = ER.Mechanism(@enzyme_mechanism begin
      substrates: A, B
      products: P, Q
      steps: begin
          E + A <--> E(; residual = A - P) + P
          E + Q <--> E(Q)
          E(; residual = A - P) + B <--> E(Q)
      end
  end)
  bibi = ER.Mechanism(@enzyme_mechanism begin
      substrates: A, B
      products: P, Q
      steps: begin
          (E + A ⇌ E(A), E(P) + A ⇌ E(A, P), E(B) + A ⇌ E(P, Q))
          (E + B <--> E(B), E(A) + B <--> E(P, Q))
          (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q))
          (E + Q <--> E(Q), E(B) + Q <--> E(B, Q), E(P) + Q <--> E(P, Q))
          E(Q) + B <--> E(B, Q)
      end
  end)
  terbi = ER.Mechanism(@enzyme_mechanism begin
      substrates: A, B, D
      products: P, Q
      steps: begin
          (E + B <--> E(B), E(D) + B <--> E(Q; residual = B + D - Q))
          (E + D ⇌ E(D), E(B) + D ⇌ E(Q; residual = B + D - Q))
          E + P <--> E(P)
          (E(; residual = B + D - Q) + A ⇌ E(P),
           E(Q; residual = B + D - Q) + A ⇌ E(A, Q; residual = B + D - Q))
          (E(; residual = B + D - Q) + Q <--> E(Q; residual = B + D - Q),
           E(P) + Q <--> E(A, Q; residual = B + D - Q))
      end
  end)
  rxn_pp2 = @enzyme_reaction begin
      substrates: A[CX], B[N]
      products: P[C], Q[NX]
      oligomeric_state: 2
  end
  rxn_ppR = @enzyme_reaction begin
      substrates: A[CX], B[N]
      products: P[C], Q[NX]
      allosteric_regulators: R
  end
  rxn_bibi = @enzyme_reaction begin
      substrates: A[C], B[N]
      products: P[C], Q[N]
  end
  rxn_terbiR = @enzyme_reaction begin
      substrates: A[C], B[N], D[X]
      products: P[CN], Q[X]
      allosteric_regulators: R
  end
  children = Tuple{String, ER.AllostericMechanism}[]
  for (name, m, rxn) in (("bi-bi ping-pong, n = 2", pp, rxn_pp2),
                         ("bi-bi ping-pong, R", pp, rxn_ppR),
                         ("three-step ping-pong, R", tc, rxn_ppR),
                         ("four-step ping-pong 3fdbe6e1dc001bb7, R", pp4, rxn_pp4),
                         ("bi-bi 1291abb91f537c6b", bibi, rxn_bibi),
                         ("ter-bi 0dc399e7dab6e4f5, R", terbi, rxn_terbiR))
      for c in ER._expand_to_allosteric(m, rxn)
          push!(children, (name, c))
          ER._is_well_posed(c) || continue
          append!(children,
                  [(name * ", change-tag", d) for d in ER._expand_change_allo_state(c)])
      end
  end
  bad = 0
  for (name, c) in children
      rng = Random.Xoshiro(hash(c))
      faults = _testhelper_face_faults(c)
      v = ER._guard(c)
      admitted = v.fault === ER._NO_FAULT
      respeed = _testhelper_respeed_verdicts(c; rng, law = admitted ? :package : :limit)
      epath = admitted ?
              [_testhelper_epath_verdict(c; rng, face) for face in _testhelper_faces(c)] :
              Symbol[]
      ok = all(_testhelper_agrees.(faults, respeed)) && all(==(:converges), epath)
      ok || (global bad += 1)
      println(ok ? "ok  " : "BAD ", name, " | ", join(string.(c.cat_allo_states), ","),
              " n=", c.catalytic_multiplicity, " sites=", length(c.regulatory_sites),
              " | guard ", v.fault, " face ", v.face, " | RE-speed ", join(respeed, ","),
              " | ε-path ", join(epath, ","))
  end
  println("VALUE CHECK: ", length(children), " children, ", bad, " disagree")
  ```
  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T12_value_check.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/probes/P4_T12_value_check.jl")'
  tail -3 $S/phase45/logs/P4_T12_value_check.log
  ```
  Expected: `VALUE CHECK: <n> children, 0 disagree`. Every child the guard admits agrees with
  the RE-speed oracle on every face and converges along the ε-path on every face; each
  child class 3 rejects (the bi-bi ping-pong's first-half children, and the change-tag
  child with both halves `:NonequalAI`) reads `mismatch` or `pattern` from the RE-speed
  oracle on its rejected faces. The research oracles gave the same verdicts on these
  children before the plan (`$S/phase45/plan/probes/p4/tagging.log`, `tagging4.log`;
  `tagging5.log` for the four-step ping-pong's six, by the research guard). A `BAD` line is
  a defect: STOP and debug it before the commit.

  Then re-pin Task 8's testset "enumeration drops the children the guard rejects"
  (`test/test_mechanism_enumeration.jl`). Its parent is the script's four-step ping-pong, whose tagging
  children are now P4-1's six: a K-type and two V-types per half-reaction base, each passing
  the guard (the script's lines for "four-step ping-pong 3fdbe6e1dc001bb7, R" read `guard
  _NO_FAULT`), so the exit drops none of them and the testset's subject becomes what the
  exit emits. Rename it "enumeration emits only children the guard admits" and replace its
  comment (Task 8's, from "# A plain ping-pong (census 3fdbe6e1dc001bb7)" to "The exit drops
  the five and keeps the flip.") with

  ```julia
      # A plain ping-pong (census 3fdbe6e1dc001bb7) whose A binding is fused into the
      # product-bound covalent form. Its one rapid-equilibrium group flips to steady state:
      # one plain child. Tagging both half-reactions' fused bindings :OnlyA, as every tag
      # set in `tagged` does, leaves the covalent forms joined to free E in the inactive
      # conformation only through tagged chemistry: an island (`_ISLAND`), whose limit
      # depends on how fast the tagged constants vanish. The tagging move builds on one
      # half-reaction per base instead (P4-1) and keeps the six children admission
      # accepts; each passes the guard, so the exit keeps them and the flip.
  ```
  After the vector `tagged` insert

  ```julia
      # Untyped: inside `T[...]` the parser reads `begin` as an index keyword.
      admitted = [
          @_testhelper_unguarded(@allosteric_mechanism begin
              substrates: A, B
              products: P, Q
              catalytic_multiplicity: 2
              steps: begin
                  E + A ⇌ E(P; residual = A - P)                        :: OnlyA
                  E(; residual = A - P) + P <--> E(P; residual = A - P) :: OnlyA
                  E(; residual = A - P) + B <--> E(Q)                   :: EqualAI
                  E + Q <--> E(Q)                                       :: EqualAI
              end
          end),
          @_testhelper_unguarded(@allosteric_mechanism begin
              substrates: A, B
              products: P, Q
              catalytic_multiplicity: 2
              steps: begin
                  E + A ⇌ E(P; residual = A - P)                        :: EqualAI
                  E(; residual = A - P) + P <--> E(P; residual = A - P) :: EqualAI
                  E(; residual = A - P) + B <--> E(Q)                   :: OnlyA
                  E + Q <--> E(Q)                                       :: OnlyA
              end
          end),
          @_testhelper_unguarded(@allosteric_mechanism begin
              substrates: A, B
              products: P, Q
              allosteric_regulators: R
              catalytic_multiplicity: 2
              steps: begin
                  E + A ⇌ E(P; residual = A - P)                        :: OnlyA
                  E(; residual = A - P) + P <--> E(P; residual = A - P) :: EqualAI
                  E(; residual = A - P) + B <--> E(Q)                   :: EqualAI
                  E + Q <--> E(Q)                                       :: EqualAI
              end
              regulatory_site(2): R::OnlyA
          end),
          @_testhelper_unguarded(@allosteric_mechanism begin
              substrates: A, B
              products: P, Q
              allosteric_regulators: R
              catalytic_multiplicity: 2
              steps: begin
                  E + A ⇌ E(P; residual = A - P)                        :: OnlyA
                  E(; residual = A - P) + P <--> E(P; residual = A - P) :: EqualAI
                  E(; residual = A - P) + B <--> E(Q)                   :: EqualAI
                  E + Q <--> E(Q)                                       :: EqualAI
              end
              regulatory_site(2): R::OnlyI
          end),
          @_testhelper_unguarded(@allosteric_mechanism begin
              substrates: A, B
              products: P, Q
              allosteric_regulators: R
              catalytic_multiplicity: 2
              steps: begin
                  E + A ⇌ E(P; residual = A - P)                        :: EqualAI
                  E(; residual = A - P) + P <--> E(P; residual = A - P) :: EqualAI
                  E(; residual = A - P) + B <--> E(Q)                   :: OnlyA
                  E + Q <--> E(Q)                                       :: EqualAI
              end
              regulatory_site(2): R::OnlyA
          end),
          @_testhelper_unguarded(@allosteric_mechanism begin
              substrates: A, B
              products: P, Q
              allosteric_regulators: R
              catalytic_multiplicity: 2
              steps: begin
                  E + A ⇌ E(P; residual = A - P)                        :: EqualAI
                  E(; residual = A - P) + P <--> E(P; residual = A - P) :: EqualAI
                  E(; residual = A - P) + B <--> E(Q)                   :: OnlyA
                  E + Q <--> E(Q)                                       :: EqualAI
              end
              regulatory_site(2): R::OnlyI
          end)]
  ```
  replace

  ```julia
      moved = EnzymeRates._expand_to_allosteric(parent, rxn)
      @test length(moved) == 5
      @test Set(moved) == Set(tagged)
  ```
  with

  ```julia
      moved = EnzymeRates._expand_to_allosteric(parent, rxn)
      @test length(moved) == 6
      @test Set(moved) == Set(admitted)
      @test all(EnzymeRates._is_well_posed, admitted)
  ```
  and

  ```julia
      kids = EnzymeRates.expand_mechanisms([parent], rxn)
      @test length(kids) == 1
      @test Set(kids) == Set([flip])
  ```
  with

  ```julia
      kids = EnzymeRates.expand_mechanisms([parent], rxn)
      @test length(kids) == 7
      @test Set(kids) == Set([flip; admitted])
  ```
  The `tagged` vector, its `_ISLAND` assertion and the ε-path check on `tagged[4]` stay:
  they pin why the move no longer emits today's tag sets. If the script's lines for this
  parent show other children, STOP: the move and the research rules disagree. Rerun the
  focused `test_mechanism_enumeration.jl` (log `P4_T12_tme.log`): expected `exit 0` and every
  testset passing but, if Step 7 found it failing, the bi-bi raw-input aggregate, which
  Step 9 re-pins.

- [ ] **Step 9: Check the move on Task 0(b)'s parents.** Task 11's queue writes the `move`
  stage's logs as `P4_T12_move_<rxn>_<lo>.log`. In the background:

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  Q=$S/phase45/pumped/run_queue.sh
  ($Q move pingpong 1 894 && $Q move bibi 1 987 && $Q move terbi 1 4916) > $S/phase45/logs/P4_T12_move_queue.log 2>&1
  python3 -I $S/phase45/pumped/summarize.py merge && python3 -I $S/phase45/pumped/summarize.py report | tee $S/phase45/logs/P4_T12_move_report.log
  ```
  Expected, after `QUEUE DONE move terbi 1:4916`: `== move: 6797 parents, 6797 with equal
  child sets, 0 with a repeated child` and no `DIFF` line. With Task 11's report this
  reproduces P4-19a's counts (Gate). A `DIFF` line names the children only one side gives;
  the move and the script implement RES §1.9 apart, so a difference is a defect in one of
  them: STOP and debug it. Once the report shows no `DIFF`, re-pin the bi-bi raw-input
  aggregate if Step 7 found it failing (Step 5: both counts to the values its failure
  prints, the reason in the commit body), and rerun the focused
  `test_mechanism_enumeration.jl` (log `P4_T12_tme2.log`): expected `exit 0`, every testset
  passes.

- [ ] **Step 10: Run the census against Task 10's.** Create
  `$S/phase45/census/explain_tagging.py`:

  ```python
  # ABOUTME: Task 12 of the phase-4 plan: classifies each difference between two census
  # ABOUTME: corpora as the tagging move's (P4-1), change-tag's (P4-20), a sample shift, or none.
  # Usage: python3 -I explain_tagging.py <corpus A> <corpus B>
  import collections, re, sys

  def rows(d, f):
      return [l.rstrip('\n').split('\t') for l in open(f'{d}/{f}')]

  def corpus(d):
      text = {r[0]: r[1] for r in rows(d, 'mechs.tsv')}
      kids = collections.defaultdict(set)
      for r in rows(d, 'children.tsv'):
          kids[(r[0], int(r[1]), r[2])].add(r[4])
      parents = {(r[0], int(r[1]), r[2]) for r in rows(d, 'parents.tsv')}
      seeds = collections.defaultdict(set)
      for r in rows(d, 'seeds.tsv'):
          seeds[r[0]].add(r[2])
      laws = {r[0]: r[1] for r in rows(d, 'laws.tsv')}
      return dict(text=text, kids=kids, parents=parents, seeds=seeds, laws=laws)

  A, B = corpus(sys.argv[1]), corpus(sys.argv[2])
  text = {**A['text'], **B['text']}
  allo = lambda i: text[i].startswith('allosteric')

  def groups(i):
      t = re.sub(r' sites .*$', '', re.sub(r'^(plain|allosteric n=\d+) ', '', text[i]))
      return frozenset(re.sub(r'\}:\w+$', '}', g) for g in t.split(' / '))

  def why_child(p, c):
      if allo(c) and groups(p) == groups(c):
          return 'change-tag (P4-20)' if allo(p) else 'tagging (P4-1)'
      return None

  def pool(C, rxn, level):
      return set().union(set(), *[v for (r, l, _), v in C['kids'].items()
                                  if r == rxn and l == level - 1])

  tally, bad = collections.Counter(), []
  for k in sorted(A['parents'] & B['parents']):
      for c in A['kids'].get(k, set()) ^ B['kids'].get(k, set()):
          why = why_child(k[2], c)
          tally[why or 'unexplained child'] += 1
          why or bad.append(f'child {c} of {k}')
  for k in sorted(A['parents'] ^ B['parents']):
      shifted = k[1] > 0 and pool(A, k[0], k[1]) != pool(B, k[0], k[1])
      tally['parent sample shift' if shifted else 'unexplained parent'] += 1
      shifted or bad.append(f'parent {k}')
  for rxn in sorted(set(A['seeds']) | set(B['seeds'])):
      for s in A['seeds'].get(rxn, set()) ^ B['seeds'].get(rxn, set()):
          tally['seed closure (P4-1)' if allo(s) else 'unexplained seed'] += 1
          allo(s) or bad.append(f'seed {s} of {rxn}')
  named = lambda C: set(C['text'])
  for key in sorted(set(A['laws']) ^ set(B['laws'])):
      shift = named(A) != named(B)
      tally['law sample shift' if shift else 'unexplained law'] += 1
      shift or bad.append(f'law {key}')
  for k, n in sorted(tally.items()):
      print(f'{k}: {n}')
  print(f'UNEXPLAINED: {len(bad)}')
  for b in bad[:20]:
      print('  ', b)
  ```
  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  $S/phase3/jl.sh $S/phase45/logs/P4_T12_census.log 1800 --project=$S/phase3/census/env_guard --heap-size-hint=2500M $S/phase3/census/corpus.jl $S/phase45/guard_wt $S/phase3/census/P4_T12 $S/phase3/census/rename_phase3.toml
  $S/phase3/jl.sh $S/phase45/logs/P4_T12_census_cmp.log 600 --project=$S/phase3/census/env --heap-size-hint=2500M $S/phase3/census/compare.jl $S/phase3/census/P4_T10 $S/phase3/census/P4_T12 $S/phase3/census/rename_phase3.toml
  python3 -I $S/phase45/census/explain_tagging.py $S/phase3/census/P4_T10 $S/phase3/census/P4_T12 | tee $S/phase45/logs/P4_T12_census_explain.log
  ```
  (The arguments Task 8 used for `P4_T8` and Task 10 for `P4_T10`, the rename map included,
  which compare.jl reads to key both corpora's values alike.) Expected: `compare.jl` reports `RESULT: FAIL` with only child-set, parent-sample,
  seed and law-sample lines (no stored-direction, value or law-string difference on a common
  law: the move changes no derivation); `explain_tagging.py` prints counts under `tagging
  (P4-1)`, `parent sample shift`, `seed closure (P4-1)`, `law sample shift` and possibly
  `change-tag (P4-20)`, and `UNEXPLAINED: 0`. The tagging lines fall where P4-1 changes
  children: the ping-pong parents of pingpong, terbi and prefix, and the bi-bi parents whose
  every-chemistry set fails existence.

  The census's level 0 is the whole init set of uniuni, bibi, pingpong and twosites, TRG's
  census reactions 1-4 (`corpus.jl`'s `REACTIONS` quotas take each), so its child sets
  hold every child TRG's census population pins count. Once it prints `UNEXPLAINED: 0`,
  re-pin each such pin Step 7 found failing to the count its failure prints (old and new
  counts in the commit body), and rerun the focused `test_reaction_graph.jl` (log
  `P4_T12_trg.log`): expected `exit 0`, every testset passes.

- [ ] **Step 11: ROADMAP.** In `ROADMAP.md`:
  - in the Ideas entry "**Admit the equilibrium bridges and the vanishing-fed basins that
    admission rejects.**", name the move that now filters this way: replace the lines Task
    11 wrote

    ```markdown
      mechanisms, none among the census mechanisms or their `:OnlyI` variants. Tag sets that
      stop the inactive conformation's catalysis one half-reaction at a time, built on the
      census ping-pong and ter-bi parents and filtered by admission, lose 5,342 candidates
    ```
    (with Task 11's count in place of 5,342 if it wrote another) with

    ```markdown
      mechanisms, none among the census mechanisms or their `:OnlyI` variants. The tagging
      move (`_tagging_bases`) stops the inactive conformation's catalysis one
      half-reaction at a time and filters by admission; on the census ping-pong and ter-bi
      parents it loses 5,342 candidates
    ```
  - delete the Known issue that Task 8 renamed "**The tagging and change-tag moves misjudge
    `:OnlyA` groups that mix a binding and a chemistry step.**" (phase 3's title "**The
    `:OnlyA` guard misjudges groups that mix a binding and a chemistry step.**"), its
    whole bullet: the moves now filter by admission, which reads each member of a mixed
    group by its own kind, and Task 4's "a mixed :OnlyA group whose limit exists is
    admitted" tests it;
  - replace the Known issue "**`BigInt` arithmetic inflates the GC's live-byte count.**",
    its whole bullet as Task 8 rewrote it, with

    ```markdown
    - **`BigInt` arithmetic inflates the GC's live-byte count.**
      Two functions that compute in `BigInt` raise `Base.gc_live_bytes` by bytes that no
      live object holds: the counter that `_split_gain_test` uses for a plain mechanism
      (`_partition_independent_count`; a `Rational{BigInt}` row reduction), by 18 KiB per
      call, and `_eq_complexity` (a `BigInt` determinant, skipped for a single
      rapid-equilibrium segment), by 15 KiB. A 6×6 determinant outside the package does
      the same in `BigInt` and `Rational{BigInt}`, not in `Rational{Int}`. The joint solve
      adds nothing, on a plain mechanism or on a concrete allosteric one; it runs in
      `Rational{Int}` on 10,226 enumerated bi-bi and ter-bi mechanisms without overflow.
      The count holds no memory and raises no peak: a process whose count reaches
      2.8 GiB peaks at 1.2 GB RSS, and the enumeration tests' count after a full
      collection stays near 1.8 GiB whether their peak RSS is 2.2 or 3.5 GB. The higher
      peak came from mechanisms that kept vectors at the capacity their builds grew to;
      `_canonicalize` stores each at its length, and `test/test_reaction_graph.jl` pins
      it. A loop that calls one of these and prints `Base.gc_live_bytes()` and
      `Base.gc_num().full_sweep` shows it. The figures come from a random-order bi-bi of
      nine steady-state steps, plain and as a two-subunit allosteric mechanism with
      `:OnlyA` on both bindings of A and on the chemical step and `:EqualAI` elsewhere
      (Julia 1.12.5, aarch64, `--heap-size-hint=2500M`).
    ```
    (`_onlya_offenders`, the third function and the source of the 100,000-run and
    full-sweep figures, is gone);
  - in the Planned bullet "**Let an allosteric ping-pong stop both half-reactions in the
    inactive conformation.**", replace the lines

    ```markdown
      at one catalytic subunit), so admission rejects it. The admission PR keeps allosteric
      ping-pong reachable by stopping one half per V-type child, which covers PDH, whose
      phosphorylation stops E1 and leaves E2 and E3 running. Stopping both halves needs the
    ```
    with

    ```markdown
      at one catalytic subunit), so admission rejects it. The tagging move keeps allosteric
      ping-pong reachable by stopping one half per V-type base (`_tagging_bases`,
      `src/mechanism_enumeration.jl`), which covers PDH, whose phosphorylation stops E1 and
      leaves E2 and E3 running. Stopping both halves needs the
    ```
  - at the end of `## Known issues`, add (open item 45; testset "change_allo_state relaxes
    a one-half V-type, leaving class 3 to the exit" of Step 1 (d) pins it):

    ```markdown
    - **Change-tag cannot give the six-step bi-bi ping-pong a live inactive conformation.**
      The tagging move gives the bi-bi ping-pong of six steps, whose second half-reaction
      is the rapid-equilibrium step `E(B; residual = A - P) ⇌ E(Q)`, one V-type per
      half-reaction, each stopping that half in the inactive conformation. On the V-type
      that stops the second half, change-tag's relaxation of both half-reactions to
      `:NonequalAI` passes admission, but the guard rejects it: at zero substrates, and at
      zero B, the covalent forms of both conformations lie on one rapid-equilibrium
      segment, joined only through free E, so the rate's limit depends on how fast each
      conformation's steps run (`_SPLIT`, `_LIMIT_SPEED`). `expand_mechanisms` drops that
      child, so the search never gives this ping-pong an inactive conformation that runs
      both half-reactions. On the three-step ping-pong, whose steps are all steady-state,
      the same relaxation passes.
    ```

  Then check that nothing outside Task 14's developer page names the deleted functions:

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && grep -rn "_onlya_offenders\|_onlya_cone_offenders\|_has_strict_positive_combination\|_onlya_haldane_violation\|_onlya_violation_message\|Stiemke" src test docs/src ROADMAP.md
  ```
  Expected: hits in `docs/src/developer.md` only, which Task 14 Step 5 rewrites; any other
  hit is rewritten now to say what the code does.

- [ ] **Step 12: Run the full suite.** `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T12_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'` in the background. Expected: `exit 0`, Aqua and JET pass (JET
  sees no call to a deleted function), no warning; `maxrss_kb` within 5% of Task 11's.

- [ ] **Step 13: Commit.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && git add src/mechanism_enumeration.jl src/thermo.jl src/types.jl ROADMAP.md test/test_mechanism_enumeration.jl test/test_types.jl test/test_oracles.jl test/test_reaction_graph.jl && git commit -F - <<'EOF'
  Build allosteric variants on one half-reaction per base when admission refuses both

  The tagging move keeps today's V-type set when admission accepts it; otherwise it
  builds on each smallest set of chemistry groups that stops the inactive
  conformation's catalysis and passes admission, so a ping-pong gets one base per
  half-reaction and stays allosteric (P4-1). K-types build on every base, and on today's
  set when admission rejects it (P4-19); change-tag keeps a relaxation admission
  accepts (P4-20). On Task 0(b)'s 6,797 parents the move equals the candidates built
  apart from it in Task 11: <88,434> children on <88,174> networks, <820> dropped by
  class 3 at the exit. The old :OnlyA test (_onlya_offenders and its Stiemke helper)
  loses its last callers and goes, with its helper test (P4-29a). Pins changed after the
  oracles checked every new child: the bi-bi ping-pong's K-types (15 to 6, three dropped
  at the exit), the ping-pong tagging and change-tag testsets rebuilt on P4-1's parents,
  Task 8's exit testset (today's five islands give way to P4-1's six admitted children),
  and the split-verdict aggregate, whose parents are now unfiltered tag sets so that both
  verdicts occur (<716, 4> to <a, r>, every rejected parent checked by the epsilon-path
  oracle)<; the bi-bi raw-input aggregate, <2831/1833> to <n/k>, by P4-1's added K-types><;
  test_reaction_graph.jl's census population pins, <5505> to <m> children of census bi-bi
  and <n3> to <m3> allosteric children of census reaction 3, by P4-1's children>.
  Census against Task 10: every difference is the move's or a sample shift it causes.
  ROADMAP gains a known issue: change-tag cannot give the six-step bi-bi ping-pong a live
  inactive conformation, since class 3 rejects its relaxation of both half-reactions.
  test_mechanism_enumeration.jl focused: <te> s, its ter-ter testset <tt> s. Full-suite
  maxrss: <kb> KB.

  Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
  EOF
  ```

---

### Task 13: Population tests (RES §4 row 12)

The sampled population tests of RES §3.3 (P6 is Task 12's), each regenerated inside the suite
from `init_mechanisms` and `expand_mechanisms` (spec §14) and judged by Task 1's oracles. All
go in `test/test_oracles.jl`, which holds the sampled population; two shared helpers move to
the fixtures file so a focused run of either test file sees them.

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (section "Phase-4
  oracles": `_testhelper_existence_agrees`; `_testhelper_plain_moves`, moved from
  `test/test_mechanism_enumeration.jl`)
- Modify: `test/test_mechanism_enumeration.jl` (the definition of `_testhelper_plain_moves`
  goes; its callers stay)
- Modify: `test/test_admission.jl` (Task 4's grid testset reads `_testhelper_existence_agrees`)
- Test: `test/test_oracles.jl` (six testsets, three helpers, the reaction constant, and a
  `let` block that builds P4's and P5's population and frees it after them)

**Interfaces:**
- Consumes: Task 2's `_guard`, `_faces`, `_face_names`; Task 4's `_admission_interior`,
  `_NO_LIMIT`, `_ARC_CYCLE`, `_ISLAND`; Task 8's `_is_well_posed`; Task 1's
  `_testhelper_faces`, `_testhelper_respeed_verdicts`, `_testhelper_epath_verdict`,
  `_testhelper_in_epath_window`, `_testhelper_lp_existence`, `_testhelper_exclusion`,
  `_testhelper_uniuni_grid`, `@_testhelper_unguarded`; Task 3's `_testhelper_face_faults`,
  `_testhelper_agrees`; phase 3's `_raw`, `Met`, `RawForm`, `RawStep`, `RawMechanism`,
  `RawSite`, `RawAllostericMechanism`, `metabolites`, `_context_bipartitions(m, g)` and
  `_apply_bipartitions(m, selection)`, `_independent_param_count`;
  `test_oracles.jl`'s `_testhelper_population_sample`, `_testhelper_population_one_round`,
  `_testhelper_content_text`.
- Produces (nothing later tasks use):
  ```julia
  _testhelper_existence_agrees(am)::Bool               # fixtures file
  _testhelper_plain_moves(m, rxn)::Vector{Union{Mechanism, AllostericMechanism}}  # moved
  _testhelper_renamed(m, rename::Dict{Symbol, Symbol})  # test_oracles.jl; unguarded
  _testhelper_reversing_rename(m)::Dict{Symbol, Symbol} # test_oracles.jl
  _testhelper_rewritten(m)                              # test_oracles.jl; steps reversed
  const _testhelper_inhibitor_rxn                       # test_oracles.jl: P4, P5
  # the inhibitor population: a let block around P4 and P5 in test_oracles.jl
  ```

**Gate:** P1-P5 and P7 pass, with the rename check of Review Focus 5; the six testsets,
Task 12's P6 and the build of the inhibitor population take at most 30 s together (P4-28,
measured in Step 5); the full suite is green with Aqua and JET, pristine.

- [ ] **Step 1: Move the shared helpers.** In `test/test_mechanism_enumeration.jl`, delete
  the comment and definition

  ```julia
  # The children of the four moves `expand_mechanisms` runs on a plain `Mechanism` `m`, in
  # its order, before its regulator-type filter.
  _testhelper_plain_moves(m, rxn) = append!(
      Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[],
      EnzymeRates._expand_re_to_ss(m), EnzymeRates._expand_split_kinetic_group(m),
      EnzymeRates._expand_add_dead_end_regulator(m, rxn),
      EnzymeRates._expand_to_allosteric(m, rxn))
  ```
  (test/test_mechanism_enumeration.jl:202-208 at 08c4043) and append to the
  section "Phase-4 oracles" of `test/mechanism_definitions_for_test_enzyme_derivation.jl`:

  ```julia
  # The children of the four moves `expand_mechanisms` runs on a plain `Mechanism` `m`, in
  # its order, before its regulator-type filter and before it drops the children the
  # guard rejects.
  _testhelper_plain_moves(m, rxn) = append!(
      Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[],
      EnzymeRates._expand_re_to_ss(m), EnzymeRates._expand_split_kinetic_group(m),
      EnzymeRates._expand_add_dead_end_regulator(m, rxn),
      EnzymeRates._expand_to_allosteric(m, rxn))

  """
  Whether admission's existence test (`_admission_interior`: no tagged plain binding inside
  one component of the inactive conformation's other steps, no cycle of such bindings)
  agrees with the exact LP (`_testhelper_lp_existence`) on `am`, and its fault with the
  test's own exclusion structure (`_testhelper_exclusion`): inside, cycle, an island, or
  none.
  """
  function _testhelper_existence_agrees(am)
      v = EnzymeRates._admission_interior(am)
      ex = _testhelper_exclusion(am, 2)
      fault = ex.failure === :inside ? EnzymeRates._NO_LIMIT :
              ex.failure === :cycle ? EnzymeRates._ARC_CYCLE :
              isempty(ex.lost) ? EnzymeRates._NO_FAULT : EnzymeRates._ISLAND
      exists = v.fault ∉ (EnzymeRates._NO_LIMIT, EnzymeRates._ARC_CYCLE)
      exists == _testhelper_lp_existence(am, 2) && v.fault === fault
  end
  ```

  In `test/test_admission.jl`, in the testset "the graph test equals the exact LP on the
  uni-uni grid" (Task 4), replace the loop

  ```julia
      for am in grid
          v = EnzymeRates._admission_interior(am)
          exists = v.fault ∉ (EnzymeRates._NO_LIMIT, EnzymeRates._ARC_CYCLE)
          @test exists == _testhelper_lp_existence(am, 2)
          ex = _testhelper_exclusion(am, 2)
          @test v.fault === (ex.failure === :inside ? EnzymeRates._NO_LIMIT :
                             ex.failure === :cycle ? EnzymeRates._ARC_CYCLE :
                             isempty(ex.lost) ? EnzymeRates._NO_FAULT : EnzymeRates._ISLAND)
      end
  ```
  with
  ```julia
      for am in grid
          @test _testhelper_existence_agrees(am)
      end
  ```

- [ ] **Step 2: Write the population tests.** Append to `test/test_oracles.jl`:

  ```julia
  # ── Phase-4 population tests (research section 3.3) ─────────────────────────

  # The census bi-bi declaration with A and I as competitive inhibitors (census/README.md,
  # "Reactions"), whose seeds' children are the population of research section 3.3's P4
  # and P5 (built in a let block around those two testsets, below).
  const _testhelper_inhibitor_rxn = @enzyme_reaction begin
      substrates: A[C], B[N]
      products: P[C], Q[N]
      competitive_inhibitors: A, I
  end

  """`m` with every metabolite renamed by `rename`, built by the constructors from its raw
  form (`_raw`), without the guard."""
  function _testhelper_renamed(m, rename::Dict{Symbol, Symbol})
      r(x::Symbol) = get(rename, x, x)
      met(x::EnzymeRates.Met) = EnzymeRates.Met(r(x.name), x.modifier)
      _testhelper_raw_form(f::EnzymeRates.RawForm) = EnzymeRates.RawForm(f.conformation,
          EnzymeRates.Met[met(x) for x in f.bound], Symbol[r(x) for x in f.added],
          Symbol[r(x) for x in f.subtracted])
      _testhelper_raw_step(s::EnzymeRates.RawStep) = EnzymeRates.RawStep(
          _testhelper_raw_form(s.from), _testhelper_raw_form(s.to),
          EnzymeRates.Met[met(x) for x in s.consumed],
          EnzymeRates.Met[met(x) for x in s.released], s.rapid_equilibrium)
      _testhelper_raw_graph(c::EnzymeRates.RawMechanism) = EnzymeRates.RawMechanism(
          Symbol[r(x) for x in c.substrates], Symbol[r(x) for x in c.products],
          Symbol[r(x) for x in c.modifiers],
          EnzymeRates.RawStep[_testhelper_raw_step(s) for s in c.steps], copy(c.group))
      raw = EnzymeRates._raw(m)
      m isa EnzymeRates.Mechanism && return EnzymeRates.Mechanism(_testhelper_raw_graph(raw))
      EnzymeRates.AllostericMechanism(EnzymeRates.RawAllostericMechanism(
          _testhelper_raw_graph(raw.catalytic), copy(raw.cat_allo_states),
          raw.catalytic_multiplicity,
          [EnzymeRates.RawSite(Symbol[r(x) for x in s.ligands], s.multiplicity,
                               copy(s.allo_states)) for s in raw.sites]))
  end

  """The renaming of `m`'s metabolites that reverses their name order, so the table order and
  the first substrate change: the first name in order becomes `Z`, the next `Y`, and so on,
  skipping the reserved `N` and `L`."""
  function _testhelper_reversing_rename(m)
      names = sort(EnzymeRates.metabolites(m))
      Dict(zip(names, [:Z, :Y, :X, :W, :U, :T, :R, :M, :H, :G][1:length(names)]))
  end

  """`m` built by the constructors from its raw form with its steps in reverse order, each
  written the other way round, without the guard: canonical form gives `m` back."""
  function _testhelper_rewritten(m)
      _testhelper_turn(s::EnzymeRates.RawStep) =
          EnzymeRates.RawStep(s.to, s.from, s.released, s.consumed, s.rapid_equilibrium)
      _testhelper_raw_graph(c::EnzymeRates.RawMechanism) = EnzymeRates.RawMechanism(
          c.substrates, c.products, c.modifiers, reverse(_testhelper_turn.(c.steps)),
          reverse(c.group))
      raw = EnzymeRates._raw(m)
      m isa EnzymeRates.Mechanism && return EnzymeRates.Mechanism(_testhelper_raw_graph(raw))
      EnzymeRates.AllostericMechanism(EnzymeRates.RawAllostericMechanism(
          _testhelper_raw_graph(raw.catalytic), raw.cat_allo_states,
          raw.catalytic_multiplicity, raw.sites))
  end

  @testset "the guard equals the RE-speed oracle on every face of the population sample" begin
      # Research section 3.3, P1: the guard admits each of the 120 sampled mechanisms, as
      # the exits they passed require, and on every face its class-3 verdict agrees with
      # the RE-speed oracle at one point; an allosteric member is judged at one catalytic
      # subunit (P4-24a).
      rng = Random.MersenneTwister(401)
      for m in _testhelper_population_sample
          @test EnzymeRates._is_well_posed(m)
          faults = _testhelper_face_faults(m)
          @test all(_testhelper_agrees.(faults, _testhelper_respeed_verdicts(m; rng)))
      end
  end

  @testset "the population's allosteric members converge along the ε-path on every face" begin
      # Research section 3.3, P2: along the standard directions and every order of the two
      # limits, each allosteric member of the sample tends to its deletion law on every
      # face, and at the interior its error falls inside the window of P4-25.
      rng = Random.MersenneTwister(402)
      allo = [m for m in _testhelper_population_sample
              if m isa EnzymeRates.AllostericMechanism]
      @test length(allo) >= 30
      for am in allo
          @test all(face -> _testhelper_epath_verdict(am; rng, face) === :converges,
                    _testhelper_faces(am))
          @test _testhelper_in_epath_window(am; rng)
      end
  end

  @testset "the graph test equals the exact LP over the grid and the population" begin
      # Research section 3.3, P3: admission's existence test agrees with the exact LP, and
      # its fault with the test's own exclusion structure, on the 108 uni-uni grid
      # mechanisms and on every allosteric member of the one-round population.
      allo = [m for m in _testhelper_population_one_round
              if m isa EnzymeRates.AllostericMechanism]
      @test length(allo) >= 30
      for am in [_testhelper_uniuni_grid(); allo]
          @test _testhelper_existence_agrees(am)
      end
  end

  # P4 and P5 read one population, built in a let block around them so that the collector
  # frees it after them: the children the moves give the inhibitor declaration's seeds
  # before expand_mechanisms drops the ones the guard rejects.
  let children = [c for m in EnzymeRates.init_mechanisms(_testhelper_inhibitor_rxn)
                  for c in _testhelper_plain_moves(m, _testhelper_inhibitor_rxn)]
  @testset "bi-bi inhibitor children pass the guard and no move tags catalysis :OnlyI" begin
      # Research section 3.3, P4: every child the moves give the census bi-bi seeds with
      # inhibitors A and I passes the guard, so expand_mechanisms drops none of them, and no
      # allosteric child tags a catalytic group :OnlyI (spec section 13). The census counts
      # 5,505 such children at phase 3 (census/base6, bibi level 0); the tagging move of
      # P4-1 changes a few.
      kids = children
      @test length(kids) >= 5400
      @test count(!EnzymeRates._is_well_posed, kids) == 0
      allo = [c for c in kids if c isa EnzymeRates.AllostericMechanism]
      @test !isempty(allo)
      @test all(c -> :OnlyI ∉ EnzymeRates.cat_allo_states(c), allo)
  end

  @testset "the guard equals the RE-speed oracle on inhibitor children's modifier faces" begin
      # Research section 3.3, P5: 60 of the children above that bind I, spaced evenly
      # through their content order; on every face that zeroes I the guard's class-3
      # verdict agrees with the RE-speed oracle at one point.
      rng = Random.MersenneTwister(405)
      binds_i = [c for c in children
                 if any(x -> x.name === :I, EnzymeRates.catalytic(c).table)]
      @test length(binds_i) >= 60
      sorted = binds_i[sortperm(_testhelper_content_text.(binds_i))]
      for c in sorted[round.(Int, range(1, length(sorted); length = 60))]
          modifier_faces = findall(face -> :I in face, _testhelper_faces(c))
          @test !isempty(modifier_faces)
          faults = _testhelper_face_faults(c)
          verdicts = _testhelper_respeed_verdicts(c; rng)
          @test all(f -> _testhelper_agrees(faults[f], verdicts[f]), modifier_faces)
      end
  end
  end

  @testset "the guard's verdict survives reordering, reversal, regrouping and renames" begin
      # Research section 3.3, P7, and Review Focus 5. The sample, plus four mechanisms the
      # guard rejects: a lumpable split rejected at Q zero only (limit test `differ`), an
      # unfed trap rejected at Q zero only (`speed`), R2's pumped inactive basin at zero
      # product, and S1's second closed class at zero substrates. Reversing the step order
      # and every step's direction gives the same mechanism. A same-model regrouping (a split
      # the cycles tie back, the count unchanged) keeps the fault. A renaming that reverses
      # the name order moves the faces' indices and the first substrate, and keeps the fault,
      # the limit verdict, the conformation and the zeroed names.
      lumpable_732 = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              (E + A <--> E(A), E(B) + A <--> E(A, B), E(P) + A <--> E(A, P))
              (E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
              (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q))
              (E + Q ⇌ E(Q), E(B) + Q ⇌ E(B, Q), E(P) + Q ⇌ E(P, Q))
              E(A, B) <--> E(P, Q)
          end
      end
      unfed_trap = @_testhelper_unguarded @enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          steps: begin
              E + A ⇌ E(A)
              E(A) + B <--> E(A, B)
              E(A, B) <--> E(Q) + P
              E + Q ⇌ E(Q)
              E(Q) + P ⇌ E(P, Q)
              E(P) + Q ⇌ E(P, Q)
              E + P <--> E(P)
          end
      end
      r2 = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)      :: OnlyA
              E(S) <--> E(P)    :: EqualAI
              F + P ⇌ E(P)      :: EqualAI
              F <--> E          :: OnlyA
          end
      end
      s1 = @_testhelper_unguarded @allosteric_mechanism begin
          substrates: A, B
          products: P
          modifiers: A
          steps: begin
              E + A <--> E(A)                 :: OnlyA
              E(A) + B <--> E(A, B)           :: EqualAI
              E(A, B) <--> E(P)               :: OnlyA
              E(P) <--> E + P                 :: EqualAI
              E + A::Mod ⇌ E(A::Mod)          :: OnlyA
              E(A::Mod) + B ⇌ E(A::Mod, B)    :: EqualAI
              E(B) + A::Mod ⇌ E(A::Mod, B)    :: EqualAI
              F + B <--> E(B)                 :: EqualAI
          end
      end
      rejected = [lumpable_732, unfed_trap, r2, s1]
      verdicts = EnzymeRates._guard.(rejected)
      @test all(v -> v.fault !== EnzymeRates._NO_FAULT, verdicts)
      @test [(v.limit, v.face) for v in verdicts[1:2]] ==
            [(EnzymeRates._LIMIT_DIFFER, 7), (EnzymeRates._LIMIT_SPEED, 7)]
      @test (verdicts[3].fault, verdicts[3].face, verdicts[3].conformation) ==
            (EnzymeRates._PUMPED, 2, 2)
      @test (verdicts[4].fault, verdicts[4].face) == (EnzymeRates._NONUNIQUE, 3)
      _testhelper_zeroed(m, v) = v.face == 0 ? Set{Symbol}() :
          Set(EnzymeRates._face_names(EnzymeRates.catalytic(m),
                                      EnzymeRates._faces(EnzymeRates.catalytic(m))[v.face]))
      for m in [_testhelper_population_sample; rejected]
          v = EnzymeRates._guard(m)
          rewritten = _testhelper_rewritten(m)
          @test rewritten == m && EnzymeRates._guard(rewritten) == v
          rename = _testhelper_reversing_rename(m)
          renamed = _testhelper_renamed(m, rename)
          w = EnzymeRates._guard(renamed)
          @test (w.fault, w.limit, w.conformation) == (v.fault, v.limit, v.conformation)
          @test _testhelper_zeroed(renamed, w) ==
                Set(rename[x] for x in _testhelper_zeroed(m, v))
      end
      # The renaming moves Q's face: it is the last face of the lumpable split, and its
      # renamed image W sorts before X, P's image, so the fault moves to face 6.
      rename_732 = _testhelper_reversing_rename(lumpable_732)
      @test EnzymeRates._guard(_testhelper_renamed(lumpable_732, rename_732)).face == 6
      n_regrouped = 0
      for m in _testhelper_population_sample
          fault, n = EnzymeRates._guard(m).fault, EnzymeRates._independent_param_count(m)
          cm = EnzymeRates.catalytic(m)
          for g in 1:EnzymeRates.ngroups(cm), bp in EnzymeRates._context_bipartitions(cm, g)
              child = EnzymeRates._apply_bipartitions(m, [(g, bp)])
              EnzymeRates._independent_param_count(child) == n || continue
              n_regrouped += 1
              @test EnzymeRates._guard(child).fault === fault
          end
      end
      @test n_regrouped >= 20
  end
  ```

  The verdict pins of the four rejected fixtures
  come from the gates of Tasks 3 and 5 (RES §3.2): `lumpable_split_732` Q = 0
  `_LIMIT_DIFFER`, `unfed_trap` Q = 0 `_LIMIT_SPEED`, R2 products = 0 `_PUMPED` in the
  inactive conformation, S1 substrates = 0 `_NONUNIQUE`. Bi-bi faces run interior, products,
  substrates, A, B, P, Q, so Q = 0 is face 7; R2's run interior, products, substrates. The
  two plain rows pin the limit test's verdict and the face here; Step 6 adds the class-3
  fault's name, read from the guard whose verdicts Task 3 confirmed, and the renamed twin
  must repeat the fault and the limit (its face moves to 6). `n_regrouped >= 20`: TO's "a
  same-model regrouping keeps every fitted name" asserts at least 20 plain same-model
  regroupings over this same sample, and this loop counts those too.

- [ ] **Step 3: Run the tests to verify they fail where they should.**

  `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T13_red.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_oracles.jl")'`.
  The guard, admission, the oracles and the move all exist, so these tests are new
  coverage of finished code, not of new code; the run shows whether the code already meets
  them. Expected: every new testset passes. A failing assertion is a defect in an earlier
  task: STOP and debug it (superpowers:systematic-debugging) with the oracle the assertion
  reads, rather than changing the test. To confirm the tests can fail, run once with the
  P7 assertion `@test _testhelper_zeroed(renamed, w) == Set(rename[x] for x in
  _testhelper_zeroed(m, v))` changed to compare against `_testhelper_zeroed(m, v)` (no rename): expected failures on the rejected fixtures
  whose zeroed names the renaming changes (`lumpable_732`, `unfed_trap`, `r2`, `s1`), then
  restore the line.

- [ ] **Step 4: Run test_admission.jl** (the refactored grid testset):
  `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T13_admission.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_admission.jl")'`.
  Expected: `exit 0`, every testset passes; "the graph test equals the exact LP on the
  uni-uni grid" passes 108 tests plus its count.

- [ ] **Step 5: Measure the population tests' time (P4-28).** Create
  `$S/phase45/probes/P4_T13_time.jl`:

  ```julia
  # ABOUTME: Times the population that Task 13's P4 and P5 build (the moves' children of the
  # ABOUTME: bi-bi seeds with inhibitors), for the P4-28 budget; run under TestEnv.
  using Test, EnzymeRates, LinearAlgebra, Random
  # A relative path would resolve against this file's directory, not the worktree's.
  include(joinpath(pwd(), "test", "mechanism_definitions_for_test_enzyme_derivation.jl"))
  rxn = @enzyme_reaction begin
      substrates: A[C], B[N]
      products: P[C], Q[N]
      competitive_inhibitors: A, I
  end
  build() = [c for m in EnzymeRates.init_mechanisms(rxn)
             for c in _testhelper_plain_moves(m, rxn)]
  kids = build()
  t = @elapsed build()
  println("inhibitor children: ", length(kids), " built in ", round(t; digits = 1), " s")
  ```
  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T13_time.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/probes/P4_T13_time.jl")'
  grep -E "inhibitor children|^(the guard|the population|the graph test|bi-bi inhibitor|ping-pong seeds)" $S/phase45/logs/P4_T13_time.log $S/phase45/logs/P4_T13_red.log
  ```
  Expected: `inhibitor children: <n> built in <t> s` with `n` at least 5,400. Add `t` (the
  second, warm build) to the `Time` column of the Test Summary lines of the seven population
  testsets in `P4_T13_red.log` (the six above and Task 12's "ping-pong seeds' tagging
  children are admitted, distinct and oracle-checked"). Expected: the sum is at most 30 s. If
  it is more, shrink in this order and rerun Steps 3 and 5: P2's faces to the interior and
  zero products (`_testhelper_faces(am)[1:2]`), P5's sample from 60 to 40, P4's seeds to
  every second one in `init_mechanisms` order, P6's seeds (Task 12) to every second one,
  keeping the seed count and base-count pins of the seeds it keeps; record the shrink in
  the commit body.

- [ ] **Step 6: Pin P4's count and P7's fault names (open items 46 and 47).**

  (a) P4: CLAUDE.md enumeration rule 1 allows an aggregate regression pin over a whole seed
  set whose comment says so. In the testset "bi-bi inhibitor children pass the guard and no
  move tags catalysis :OnlyI", replace `@test length(kids) >= 5400` with
  `@test length(kids) == <n>`, where `<n>` is the count Step 5's last run printed
  (`inhibitor children: <n>`; if Step 5 shrank P4's seeds, `build` in `P4_T13_time.jl`
  reads the same seeds and is rerun first), and replace the comment's
  last sentence, "The census counts 5,505 such children at phase 3 (census/base6, bibi
  level 0); the tagging move of P4-1 changes a few.", with "An aggregate regression pin
  over the declaration's whole seed set (CLAUDE.md enumeration rule 1's exception): 5,505
  children at phase 3 (census/base6, bibi level 0), <n> with phase 4's moves." If Step 5
  shrank P4's seeds, write "every second seed" for "the declaration's whole seed set".

  (b) P7: Task 3's testset "the limit test decides what structure cannot" confirmed the
  guard's limit verdicts on `lumpable_split_732` and `unfed_trap` against the RE-speed
  oracle; read their class-3 fault names from the guard. Create
  `$S/phase45/probes/P4_T13_faults.jl`:

  ```julia
  # ABOUTME: Prints the guard's verdict on P7's two plain rejections, the lumpable split and
  # ABOUTME: the unfed trap, whose class-3 fault names Task 13 pins; run under TestEnv.
  using Test, EnzymeRates, LinearAlgebra, Random
  # A relative path would resolve against this file's directory, not the worktree's.
  include(joinpath(pwd(), "test", "mechanism_definitions_for_test_enzyme_derivation.jl"))
  lumpable_732 = @_testhelper_unguarded @enzyme_mechanism begin
      substrates: A, B
      products: P, Q
      steps: begin
          (E + A <--> E(A), E(B) + A <--> E(A, B), E(P) + A <--> E(A, P))
          (E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
          (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q))
          (E + Q ⇌ E(Q), E(B) + Q ⇌ E(B, Q), E(P) + Q ⇌ E(P, Q))
          E(A, B) <--> E(P, Q)
      end
  end
  unfed_trap = @_testhelper_unguarded @enzyme_mechanism begin
      substrates: A, B
      products: P, Q
      steps: begin
          E + A ⇌ E(A)
          E(A) + B <--> E(A, B)
          E(A, B) <--> E(Q) + P
          E + Q ⇌ E(Q)
          E(Q) + P ⇌ E(P, Q)
          E(P) + Q ⇌ E(P, Q)
          E + P <--> E(P)
      end
  end
  for (name, m) in (("lumpable_split_732", lumpable_732), ("unfed_trap", unfed_trap))
      v = EnzymeRates._guard(m)
      println(name, ": ", v.fault, " ", v.limit, " face ", v.face)
  end
  ```
  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T13_faults.log 900 --project=. -e 'using TestEnv; TestEnv.activate(); include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/phase45/probes/P4_T13_faults.jl")'
  ```
  Expected: `lumpable_split_732: _SPLIT _LIMIT_DIFFER face 7` and
  `unfed_trap: _SPLIT _LIMIT_SPEED face 7`, the verdicts the refresh pass (P4-32) read from
  Tasks 1-3's code at 08c4043, where every testset of Task 3 passed
  (`$S/phase45/refresh/apply_probe47.jl`, log `apply_47.log`). A limit or face other than
  P7's pins, or another fault name, is checked against `_face_structure` at Q = 0 and the
  message before anything is pinned. Then replace P7's assertion
  ```julia
      @test [(v.limit, v.face) for v in verdicts[1:2]] ==
            [(EnzymeRates._LIMIT_DIFFER, 7), (EnzymeRates._LIMIT_SPEED, 7)]
  ```
  with the names the run printed:
  ```julia
      @test [(v.fault, v.limit, v.face) for v in verdicts[1:2]] ==
            [(EnzymeRates._SPLIT, EnzymeRates._LIMIT_DIFFER, 7),
             (EnzymeRates._SPLIT, EnzymeRates._LIMIT_SPEED, 7)]
  ```

  Then rerun Step 3's command, log `P4_T13_pinned.log`. Expected: `exit 0`, every new
  testset passes.

- [ ] **Step 7: Run the full suite.** `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T13_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'` in the background. Expected: `exit 0`, Task 12's pass count plus
  the new tests, Aqua and JET pass, no warning; `maxrss_kb` within 5% of Task 12's (the
  inhibitor population, about 5,500 mechanisms, lives only while P4 and P5 run; a larger
  rise is explained in the commit body). The commit body gives both figures, Task 12's
  before and this run's after.

- [ ] **Step 8: Commit.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && git add test/test_oracles.jl test/test_admission.jl test/test_mechanism_enumeration.jl test/mechanism_definitions_for_test_enzyme_derivation.jl && git commit -F - <<'EOF'
  Check the guard against the oracles over regenerated populations

  Over the sampled ping-pong-capable bi-bi population the guard's face verdicts agree
  with the RE-speed oracle and every allosteric member converges along the epsilon-path;
  the existence test equals the exact LP over the uni-uni grid and the population; every
  round-one child of the bi-bi seeds with inhibitors passes the guard and none tags
  catalysis :OnlyI, and their modifier faces agree with the oracle; the guard's verdict
  holds under step order, direction, same-model regrouping and a renaming that reorders
  the table and the first substrate. P4 pins the population's size, <n> children (5,505
  at phase 3), as an aggregate regression pin, and P7 pins both plain rejections as
  splits. Phase 4's population tests take <t> s (P4-28 allows 30 s). Full-suite maxrss:
  <kb0> KB before, <kb> KB after.

  Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
  EOF
  ```

---

### Task 14: Docs, ROADMAP, the R18 scan and the draft PR (RES §4 row 13)

The pages phase 4 makes false (spec §16: `deriving/mwc_allostery.md` and the `:OnlyI` rule in
the macro docstrings), the identify and developer pages that describe the moves and the guard,
the ROADMAP entry that counts the remaining PRs, GC rule R18's remark scan, and the draft PR.
Phase 3's S5 rewrote these pages for the graph types after the snapshot this plan was written
on, so each edit names its paragraph by heading and content: replace the paragraph that says
what the old text says, wherever S5 left it, with the new text given here. The refresh
pass (P4-32) located each paragraph at 08c4043; the old texts below quote it there.

**Files:**
- Modify: `docs/src/deriving/mwc_allostery.md`, `src/dsl.jl` (the `@allosteric_mechanism`
  docstring), `docs/src/identify/enumeration_engine.md`, `docs/src/developer.md`
- Modify: `ROADMAP.md` (Planned "Rewrite the derivation around the reaction graph": one PR
  remains)
- Create (scratch): `$S/phase45/gcscan/mkenv.jl`, `derive.jl`, `identify.jl`; environments
  `$S/phase45/gcscan/env112`, `env113`; `$S/phase45/pr_body.md`
- Logs: `$S/phase45/logs/P4_T14_{docs,full,remarks_*}.log`

**Interfaces:**
- Consumes: every earlier task's behavior; Task 8's census and ter-ter figures, Task 11's
  report (`P4_T11_report.log`), Task 12's move check and census explanation, Task 13's
  population time; `$S/phase3/gcprobe/hoistscan.py` (the remark scanner, final after the
  outage per `RECOVERY.md`; it reads only Julia's JuliaLICM remarks).
- Produces: the draft PR graph-guard into graph-mechanism.

**Gate:** `docs/make.jl` builds with its doctests and its check that every export has a
docstring, no warning; the R18 remark scan of both workloads on Julia 1.12 and 1.13 reports 0
hoisted allocations with a bits store into a reference field and no hoist at a phase-4
function; the full suite is green and pristine; the PR body lists every ruling P4-1 to P4-31
with the six Denis may overturn first (P4-8, P4-9, P4-11, P4-12, P4-25, P4-29), the verdicts
that changed, P4-29's removed subjects for approval, R-B1 and R-B2 (P4-9a), Task 11's count,
the mod-p residual as the largest part he could delete (P4-11a, about 300 lines), P4-19a's
counts, and version 0.10.0 (R19).

- [ ] **Step 1: Read the pages as phase 3 left them.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && grep -n "OnlyI\|second free form\|N_I = 0\|rejects the inconsistent" docs/src/deriving/mwc_allostery.md
  grep -n "Promote a non-allosteric\|Relax an\|Every child is checked\|bottomless" docs/src/identify/enumeration_engine.md
  grep -n "Columns outside\|_onlya\|bottomless\|_is_chemistry\|check'" docs/src/developer.md
  grep -n "A catalytic statement takes" src/dsl.jl
  ```
  Expected: each paragraph Steps 2-5 replace is found once. No earlier task edits
  `docs/src/developer.md`: at 08c4043 it never names `_onlya_haldane_violation`, and it
  names `_onlya_offenders` twice, in the sentence on `_is_chemistry`'s readers (Step 5
  (d)) and in the split paragraph (Step 5 (c)).

- [ ] **Step 2: `docs/src/deriving/mwc_allostery.md`.**

  (a) Section "Setting up an MWC mechanism": in the paragraph after the tag table, replace the
  sentence that says a catalytic step can never be `:OnlyI` ("A catalytic step can never be
  `:OnlyI` — the active state must always be able to catalyze — so tagging one that way is a
  hard constructor error.") with:

  ```markdown
  A catalytic group takes `:OnlyI` only when its steps bind a modifier and nothing else (a
  name under `modifiers:` or a reactant's `::Mod` copy), such as an inhibitor that binds only
  the inactive conformation. The active conformation must still catalyze, so the
  constructor rejects `:OnlyI` on a step that takes up or gives off a substrate or product.
  ```

  (b) Section "Derivation of the MWC rate equation", first paragraph: replace the two
  sentences on a second free form ("A hand-written mechanism may name a second free form,
  … which would flip at different forms in the two conformations.") with:

  ```markdown
  The conformations interconvert only at free E, the form `E` with nothing bound, and `L` is
  its inactive-to-active ratio. A mechanism may hold other free forms, such as a
  conformation `F` with `F <--> E`; they change conformation only by way of E. When another
  free form sits at rapid equilibrium with E, its own inactive-to-active ratio is `L` times
  the ratio of the inactive to the active equilibrium constant of the step that joins it to
  E, so a fit reports `L` at E whatever the free forms are named.
  ```

  (c) Same section, the paragraph on `D`: replace its last two sentences ("When a catalytic
  group is `:OnlyA`, the inactive cycle cannot close, so the package sets `N_I = 0`: the
  inactive branch still contributes to the denominator as enzyme mass but carries no forward
  flux.") with:

  ```markdown
  When the `:OnlyA` groups stop every catalytic cycle of the inactive conformation, as the
  enumerated tag sets do, its numerator `N_I` is 0: the inactive branch still contributes
  to the denominator as enzyme mass but carries no forward flux.
  ```

  (d) Insert, before the heading `## Thermodynamic constraints of MWC equations`:

  ```markdown
  ## Tags as limits, and admission

  An `:OnlyA` tag says what the inactive conformation lacks, as the limit of ordinary
  constants; `:OnlyI` says the same of the active conformation:

  - a plain binding loses its affinity: its dissociation constant grows without bound, or, for
    a steady-state binding, its binding rate constant falls to zero with its release rate
    fixed;
  - any other step (an isomerization, a fused binding, a Theorell–Chance step) is blocked:
    both rate constants fall to zero with their ratio held by the thermodynamic cycles, and a
    rapid-equilibrium step keeps the constant the cycles set while its speed vanishes;
  - each step of a tagged group vanishes at its own rate, and a group may mix a plain and a
    fused binding of one metabolite.

  The derivation reads the inactive conformation's law by deleting the tagged steps and the
  forms that the deletion cuts off from free E. **Admission** decides when that is the limit,
  whatever the rates at which the tagged constants vanish:

  1. **Existence.** Join the forms by the conformation's untagged steps, and draw each tagged
     plain binding as an arrow from its unbound side's component to its bound side's. A path
     to the limit exists exactly when no arrow stays inside one component and the arrows
     close no cycle.
  2. **Uniqueness.** Every form that free E reaches only through tagged steps lies behind a
     tagged plain binding, on its bound side, and the limit is the same on every face where
     data sit (below).
  3. **Turnover.** With `:OnlyI` groups, the active conformation still turns the reaction
     over.

  A tag set that fails is rejected with an error that names the step. Tagging `E + S ⇌ E(S)`
  `:OnlyA` while the catalytic step and the product binding stay `:EqualAI` fails existence:
  in the inactive conformation those steps still join E to E(S), so no path takes that
  affinity to zero. Admission is sound but not exact. It also rejects two rare shapes whose
  limit exists, equilibrium bridges and stranded basins fed only through vanishing steps,
  which the repository's `ROADMAP.md` lists.

  ### What the package rejects

  A mechanism is well defined when its rate law is the unique fast-equilibrium limit of
  mass-action kinetics wherever data sit: at interior concentrations and on every **face**,
  where every product is zero, every substrate, each reactant alone on a side of two or more,
  or any of these with every modifier-only name zero. The macros throw an error that names
  the step and the face; `init_mechanisms` and `expand_mechanisms` drop the mechanism. The
  package rejects:

  1. no catalytic cycle: no cycle of steps turns substrates into products, so the rate is
     zero;
  2. a cycle of rapid-equilibrium steps that turns the reaction over, so the rate is
     infinite;
  3. a face where the law is not the unique limit: a rapid-equilibrium segment with no form
     free of the face's metabolites, a segment whose lowest forms are joined only through
     other forms, a trap that the face feeds and its fast exits cannot drain, or a second
     closed set of forms;
  4. tags with no limit (admission, above);
  5. a denominator that vanishes on a face, found when the law is first derived.

  The faces leave out the corners that zero one substrate with every product. A law the
  package accepts can be 0/0 there, where a well-defined law is 0, so `rate_equation`
  returns `NaN` at such a point; a fit never meets one, since a measured rate there is zero.
  ```

  (e) Section "Thermodynamic constraints of MWC equations": replace the paragraph that begins
  "Often a K-system uses exclusive `:OnlyA` binding" and ends "`@allosteric_mechanism`
  rejects the inconsistent combination." with:

  ```markdown
  Often a K-system uses exclusive `:OnlyA` binding: the substrate binds the active state
  alone, and the inactive cycle is pruned (`N_I = 0`). Its catalytic step **must** be
  `:OnlyA` too. This is a thermodynamic requirement, not a convention: with `:EqualAI` or
  `:NonequalAI` catalysis and an `:EqualAI` product binding, the inactive conformation still
  joins E to E(S) through catalysis and product release, so no path takes the substrate
  affinity to zero while every cycle keeps `Keq`. Admission's existence test rejects the
  combination.
  ```

  The paragraph that follows ("The exception is a balanced K-system. …") stays: with both
  bindings `:OnlyA` the arrows leave E's component and close no cycle, so the limit exists.

- [ ] **Step 3: The `@allosteric_mechanism` docstring (`src/dsl.jl`).** Replace

  ```julia
  one only. A catalytic statement takes `EqualAI`, `NonequalAI` or `OnlyA`; a ligand takes
  any of the four.
  ```
  with
  ```julia
  one only. A catalytic statement takes `EqualAI`, `NonequalAI` or `OnlyA`, or `OnlyI` when
  its steps bind only modifiers (`modifiers:` names or `::Mod` copies); a ligand takes any of
  the four. The conformations interconvert only at free E.
  ```
  (the first line of the old text continues the sentence "… `OnlyA` acts in the active
  conformation only and `OnlyI` in the inactive"; keep that sentence's start as phase 3 and
  Task 8 left it).

- [ ] **Step 4: `docs/src/identify/enumeration_engine.md`.**

  (a) Section `EnzymeRates.init_mechanisms`, the bullet **Valid.**: replace it with

  ```markdown
  - **Valid.** No cycle of rapid-equilibrium steps runs the net reaction
    (`_re_turnover_cycle`), and the guard's face check accepts every face (`_faces_pass`,
    read on the candidate's flags).
  ```

  (b) Section "The seven expansion moves", first paragraph: replace its last sentence ("Every
  child is checked to conserve atoms before it is returned.") with

  ```markdown
  Every child is checked to conserve atoms, and a child the well-posedness guard rejects
  (see [Mechanisms with allosteric regulators](@ref)) is dropped.
  ```

  (c) Section "### 4. Promote a non-allosteric to allosteric mechanism": replace its first
  paragraph (from "Converts a `Mechanism` to an `AllostericMechanism` variant set." to "No-op
  on an already allosteric input.") with

  ```markdown
  Converts a `Mechanism` to `AllostericMechanism` variants, each built on a **base** tag set
  (`_tagging_bases`). A chemistry group holds an isomerization, a fused binding or a
  Theorell–Chance step (`_is_chemistry`); every other group is a binding group. The first
  candidate base tags every chemistry group `:OnlyA` and every binding group `:EqualAI`: an
  inactive conformation that binds ligands but runs no chemistry. When admission accepts it
  (see [Mechanisms with allosteric regulators](@ref)), it is the only base. Otherwise the
  bases are the smallest sets of chemistry groups, tagged `:OnlyA` with every other group
  `:EqualAI`, that stop every catalytic cycle of the inactive conformation and pass
  admission. A ping-pong mechanism is the common case: tagging both half-reactions strands
  the covalent forms between them, whose weight then depends on how fast each half stops,
  so it gets one base per half-reaction, the other half running in both conformations, as
  phosphorylation stops PDH's E1 and leaves E2 and E3 running. Stopping both halves needs the
  conformations to flip at the covalent form too, with a constant of its own, which waits for
  iso mechanisms.

  On its bases the engine emits, per multiplicity:

  - **K-types:** every nonempty set of binding groups tagged `:OnlyA` on top of a base, and,
    when admission rejects the every-chemistry set, on top of that set too, since a tagged
    binding can restore the limit (as in a bi-bi whose group mixes a plain binding with a
    fused one); each is emitted bare and kept when admission accepts it;
  - **V-types:** each base paired with a declared allosteric regulator at a new site, one
    variant per `(regulator, tag)` with `tag ∈ {:OnlyA, :OnlyI}`, where the regulator makes
    `L` identifiable.

  The all-`:EqualAI` baseline is never emitted (the conformations would be identical and `L`
  would cancel). Admission does not run the guard's face check, so a child that fails it is
  dropped where `expand_mechanisms` returns it, as for every move. A ping-pong parent with no
  declared regulator can end with no child, when every K-type strands a form or a basin.
  Enumeration runs over `allowed_catalytic_multiplicities`. No-op on an already allosteric
  input.
  ```

  (d) Section "### 6. Relax an `:EqualAI` group to independent A/I states": append to its
  first paragraph

  ```markdown
  The chemistry groups relax together, and only while no binding group is `:OnlyA`. A
  catalytic relaxation is kept when admission accepts it (`_admits`); a child the guard's face
  check rejects is dropped where `expand_mechanisms` returns it.
  ```

- [ ] **Step 5: `docs/src/developer.md`.**

  (a) Insert, before the heading `## Enumeration engine architecture`:

  ```markdown
  ## Well-posedness guard

  `src/admission.jl` decides whether a mechanism's rate law is the unique fast-equilibrium
  limit of mass-action kinetics wherever data sit. `_guard(m)` reads the reaction graph
  through masks and flag overrides and never constructs a mechanism, so a move can judge a
  candidate on its flags (`_faces_pass(m, re)`). It returns an isbits `_Verdict` (fault,
  limit-test verdict, face, conformation, edge), never a `Union` that a loop reads, the rule
  the GC probe set for hot code. `_check_well_posed(m)` throws `_guard_message(m, v)`, which
  names each step by its DSL line (`_step_text`) and the face; `_is_well_posed(m)` is the
  enumerator's drop test.

  The guard builds two networks of one catalytic subunit whose conformations flip at free E
  (`_Network`): the **law network**, each conformation's kept steps on the forms free E
  reaches through them (`_law_mask`), and the **ε-network**, every step of both
  conformations, each excluded step marked vanishing and an excluded rapid-equilibrium plain
  binding kept fast (`_eps_network`). The faces (`_faces`) are bit masks over the metabolite
  names. In order it checks class 1, a cycle that turns the reaction over
  (`_has_turnover_cycle`); class 2, no fast cycle of the ε-network that takes up a name net;
  admission's interior (`_admission_interior`: existence by components and arcs, islands,
  the active conformation's turnover); then on each face class 3 (`_check_face`: closed
  classes, bottomless segments, the aggregated chains, and a limit test modulo 2^61 − 1 where
  structure decides nothing) and admission's face condition (`_admission_face`: the
  ε-network's closed classes, stranded basins, and the fewest-slow-transitions in-tree test
  per excluding conformation and for both together). Class 5, a denominator that keeps no
  monomial on a face, runs on the derived law in `_num_den_exprs`. The macros, the parser
  and the derivation entry throw; `init_mechanisms`, `seed_mechanisms`, `expand_mechanisms`
  and identify's `_base_tier` drop; the flip move runs class 3 on each candidate's flags; the
  tagging and change-tag moves filter by admission alone (`_admits`) and leave class 3 to the
  exit.
  ```

  (b) Section "The joint thermodynamic solve": replace the sentence "Columns outside the
  inactive conformation's law (an `:OnlyA` edge's inactive copy, an edge on a form the
  inactive conformation cannot reach) come first, so no kept dependent reads them." with

  ```markdown
  Columns outside a conformation's law (an `:OnlyA` edge's inactive copy, an `:OnlyI` edge's
  active copy, an edge on a form that conformation cannot reach) come first, so no kept
  dependent reads them.
  ```

  (c) Section "Enumeration engine architecture", the paragraph on the split: replace the
  sentences "A candidate whose `:OnlyA` binding leaves a thermodynamic cycle
  unsatisfiable, which the boundary check rejects, is never built: … since the split that
  balances its cycle may be steady state." with

  ```markdown
  A split never changes admission, since both parts keep the group's tag, and a child the
  guard rejects is dropped where `expand_mechanisms` returns it. The search extends a failed
  set by the splits of the groups that share its rapid-equilibrium segments.
  ```

  (d) Same section: in the sentence that lists the readers of `_is_chemistry`, replace
  "the `:OnlyA` test (`_onlya_offenders`)" with "admission"; in the paragraph on
  `init_mechanisms`' variants, replace "leaves no bottomless segment" with "passes the
  guard's face check".

  (e) Same section, the paragraph on what the constructors check: replace its last two
  sentences ("Whether a mechanism's rate law is defined is the boundary check's question
  … the constructors and the Sig decoder do not.", which name the one-flip-segment and
  `:OnlyA` cycle conditions that Tasks 8 and 9 delete) with

  ```markdown
  Whether a mechanism's rate law is defined is the well-posedness guard's question
  (`src/admission.jl`; see [Well-posedness guard](@ref)). The macros, the parser and every
  derivation entry throw on a mechanism it rejects (`_check_well_posed`), and the
  enumerator's exits and identify's base tier drop one (`_is_well_posed`); the
  constructors and the Sig decoder do not run it.
  ```

  (f) Same section, the paragraph on `_hyperbolic_catalysis`: replace "which the boundary
  check's `_bottomless_segment` also reads" with "which the split move's
  `_bottomless_segment` also reads" (after Tasks 8 and 9 the split move is its only
  caller).

- [ ] **Step 6: ROADMAP.** In `ROADMAP.md`, under `## Planned`, replace the entry
  "**Rewrite the derivation around the reaction graph.**" (its whole bullet, as phase 3 left
  it) with

  ```markdown
  - **Rewrite the derivation around the reaction graph.**
    `docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md` records the design,
    decided between 2026-10-09 and 2026-10-12. One PR remains: the engine, the numerator,
    normalization and printing, a rename-proof dedup key, identify's key cache and beam,
    the L-phantom flag, the code-generation locals and the speed gates' re-calibration
    (its section 15, phase 5). The table in its section 16 lists the entries it retires.
    Priority: high.
  ```

- [ ] **Step 7: Build the docs.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T14_docs.log 3600 --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path = pwd())); Pkg.instantiate(); include("docs/make.jl")'
  grep -n "Error\|Warning\|doctest" $S/phase45/logs/P4_T14_docs.log | head -20
  git status --short docs/Project.toml docs/Manifest.toml
  ```
  Expected: `exit 0`; no `Error` and no `Warning` from Documenter except deploydocs' note
  that it skips deployment outside CI; every doctest passes (the `@allosteric_mechanism` and
  `@enzyme_mechanism` doctests included); no missing-docstring report. If `Pkg.develop`
  changed `docs/Project.toml` or created `docs/Manifest.toml`, restore them
  (`git checkout docs/Project.toml`, `rm -f docs/Manifest.toml`) so the commit carries only
  the page edits.

- [ ] **Step 8: Run the full suite.** `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && $S/phase3/jl.sh $S/phase45/logs/P4_T14_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'` in the background. Expected: `exit 0`, Task 13's pass count, Aqua
  and JET pass, no warning; `maxrss_kb` within 5% of Task 13's.

- [ ] **Step 9: Commit the docs.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && git add docs/src/deriving/mwc_allostery.md docs/src/identify/enumeration_engine.md docs/src/developer.md src/dsl.jl ROADMAP.md && git commit -F - <<'EOF'
  Document the flip at free E, the tags as limits, admission and the guard

  The MWC page states the flip at free E and what L means there, the limits the tags
  stand for, admission's three checks and the five rejected classes; the macro docstring
  gives the :OnlyI rule for modifier bindings; the enumeration page describes the tagging
  move's bases (one per half-reaction on a ping-pong), the exit drops and change-tag's
  admission filter; the developer page describes the guard. ROADMAP's rewrite entry now
  counts one remaining PR.

  Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
  EOF
  ```

- [ ] **Step 10: Write the R18 workloads.** Create `$S/phase45/gcscan/mkenv.jl`:

  ```julia
  # ABOUTME: Creates an environment that develops the guard worktree by path and holds the
  # ABOUTME: optimizer the identify workload needs, so remark scans never touch its manifest.
  using Pkg
  Pkg.offline(true)
  Pkg.activate(ARGS[2])
  Pkg.develop(path = ARGS[1])
  Pkg.add(["Random", "OptimizationCMAEvolutionStrategy"])
  Pkg.instantiate()
  using EnzymeRates
  println("ok ", pkgversion(EnzymeRates), " from ", pathof(EnzymeRates))
  ```

  `$S/phase45/gcscan/derive.jl`:

  ```julia
  # ABOUTME: R18 remark-scan workload: enumerates a ping-pong bi-bi population with an
  # ABOUTME: allosteric regulator (the guard at every exit) and derives laws and messages.
  using EnzymeRates
  const rxn = @enzyme_reaction begin
      substrates: A[CX], B[N]
      products: P[C], Q[NX]
      allosteric_regulators: R
      oligomeric_state: 2
  end
  seeds = EnzymeRates.init_mechanisms(rxn)
  kids = EnzymeRates.expand_mechanisms(seeds, rxn)
  pop = [seeds; kids]
  total = 0
  for m in pop[1:min(end, 150)]
      global total += length(EnzymeRates.parameters(m)) +
                      length(EnzymeRates.rate_equation_string(m))
  end
  rejected = try
      @allosteric_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)      :: OnlyA
              E(S) <--> E(P)    :: EqualAI
              F + P ⇌ E(P)      :: EqualAI
              F <--> E          :: OnlyA
          end
      end
      "admitted"
  catch err
      sprint(showerror, err)
  end
  println("population ", length(pop), "; checksum ", total, "; R2: ", first(rejected, 40))
  ```

  and `$S/phase45/gcscan/identify.jl`:

  ```julia
  # ABOUTME: R18 remark-scan workload: a small identify_rate_equation search on a uni-uni
  # ABOUTME: reaction with an inhibitor, the suite's smoke settings, three rounds.
  using EnzymeRates, Random, OptimizationCMAEvolutionStrategy
  rxn = @enzyme_reaction begin
      substrates: S[C]
      products: P[C]
      competitive_inhibitors: R
  end
  truth = @allosteric_mechanism begin
      substrates: S
      products: P
      allosteric_regulators: R
      steps: begin
          E + S ⇌ E(S)      :: EqualAI
          E(S) <--> E(P)    :: OnlyA
          E + P ⇌ E(P)      :: EqualAI
      end
      regulatory_site(1): R::OnlyI
  end
  fp = EnzymeRates.fitted_params(truth)
  values = ((1.0 + 0.1i for i in eachindex(fp))..., 2.0, 1.0)
  params = NamedTuple{(fp..., :Keq, :E_total)}(values)
  Random.seed!(42)
  rows = [(g, 0.1 + 9.9rand(), 0.1 + 9.9rand(), 0.1 + 9.9rand()) for g in 1:5 for _ in 1:10]
  data = (group = ["G$(r[1])" for r in rows], S = [r[2] for r in rows],
          P = [r[3] for r in rows], R = [r[4] for r in rows],
          Rate = [rate_equation(truth, (S = r[2], P = r[3], R = r[4]), params) for r in rows])
  prob = IdentifyRateEquationProblem(rxn, data; Keq = 2.0)
  for round in 1:3
      res = redirect_stdout(devnull) do
          identify_rate_equation(prob; min_beam_width = 1, loss_rel_threshold = 1.0,
              loss_abs_threshold = 0.0, max_param_count = 8, n_cv_candidates = 1,
              save_dir = mktempdir(), optimizer = CMAEvolutionStrategyOpt(),
              n_restarts = 1, maxtime = 1.0, show_progress = false)
      end
      println("round ", round, ": ", size(res.cv_results, 1), " cv rows")
      GC.gc(true)
  end
  ```

  The two workloads follow phase 3's GC probe (`$S/phase3/gcprobe/pkg_stress.jl`,
  `pkg_identify.jl`), rewritten on the graph API: the identify workload's truth comes from the
  macro, its parameter names from `fitted_params`, so no name is spelled by hand. The derive
  workload adds the macro path on a mechanism the guard rejects (R2), so the message code is
  compiled under the remark flags too.

- [ ] **Step 11: Run the remark scan on Julia 1.12 and 1.13.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt
  G=$S/phase45/gcscan; cd $S/phase45/guard_wt
  $S/phase3/jl.sh $S/phase45/logs/P4_T14_env112.log 1800 $G/mkenv.jl $PWD $G/env112
  $S/phase3/jl.sh $S/phase45/logs/P4_T14_env113.log 1800 +1.13 $G/mkenv.jl $PWD $G/env113
  for v in 112 113; do for w in derive identify; do
    ch=$([ $v = 113 ] && echo +1.13)
    JULIA_LLVM_ARGS="-pass-remarks=julia-licm -pass-remarks-analysis=julia-licm" \
      $S/phase3/jl.sh $S/phase45/logs/P4_T14_remarks_${w}_${v}.log 3600 $ch --project=$G/env$v --pkgimages=no --heap-size-hint=2500M $G/$w.jl
  done; done
  python3 -I $S/phase3/gcprobe/hoistscan.py $S/phase45/logs/P4_T14_remarks_*_1*.log
  ```
  (Each scan compiles the package and its dependencies with the remarks on, a few minutes
  apiece; run them in the background one at a time.) Expected: each workload log ends
  `exit 0` (`population <n>; …; R2: AllostericMechanism: …` and three `round` lines); the
  scanner prints, for each of the four logs, `escape analyses: <n>  hoisted allocations: 0
  hoisted with a bits store into a reference field: 0`, as phase 3's final scan found for
  its derivation and identify workloads at 67eaea5 and at phase 2, on Julia 1.12 and 1.13
  (`$S/phase3/replay/logs/T63_remarks_scan.txt`; the GC probe's package scan at 4f9a34a,
  `$S/phase3/gcprobe/report.txt` section 4, found 0 too). A `hoist` or
  `DANGER` line at a phase-4 function (`src/admission.jl`, `_tagging_bases`, `_admits`) is a
  defect: STOP and remove the `Union` the loop reads (R18) before the PR. A hoist outside
  phase-4 code is compared with the same scan on the branch point
  (`git worktree add --detach $G/base_wt $(git merge-base HEAD graph-mechanism)`, envs built
  the same way) and reported in the PR body.

- [ ] **Step 12: Write the PR body.** Create `$S/phase45/pr_body.md` with the text below,
  each `<…>` replaced by the figure its source names:

  ```markdown
  ## Summary

  Phase 4 of the reaction-graph rewrite (spec `docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md`,
  sections 7, 8, 13, 14). A mechanism is fitted only when its rate law is the unique
  fast-equilibrium limit of mass-action kinetics wherever data sit; the rest is rejected
  with an error naming the step and the face, or dropped in enumeration.

  - **Guard** (`src/admission.jl`): class 1 (no catalytic cycle), class 2 (RE turnover, on
    the ε-network), class 3 per face (closed classes, bottomless segments, aggregated chains,
    a mod 2^61 − 1 limit test), class 5 at derivation (a denominator that vanishes on a face).
    Isbits `_Verdict`, no `Union` read in a loop (R18).
  - **Admission** (class 4): existence by components and arcs (equal to an exact LP on the
    108 uni-uni grid mechanisms and the population), islands behind tagged bindings, active
    turnover, and the face condition B v3 (ε-network closed classes, stranded basins, R1 or
    T0+ per conformation, joint R1).
  - **Flip at free E**; the one-flip-segment check is gone.
  - **`:OnlyI` on modifier-binding groups**, mirrored in the active conformation.
  - **Enumeration:** every exit drops what the guard rejects; class 3 runs inside the flip
    move; the tagging move stops one half-reaction per base on a ping-pong (P4-1);
    change-tag filters by admission; the split's allosteric pre-check is gone.
  - **Tests:** the RE-speed, ε-path, LP and n-protomer oracles in the shared fixtures file;
    `test/test_admission.jl`; population tests P1-P7.
  - Version stays 0.10.0 (R19). Draft, stacked on the phase-3 PR; not for merge.

  ## For Denis first

  The six rulings he may overturn first:
  - **P4-8** B v3 is the face condition; the neutral-bridge exemption (B v2b, about 90 lines)
    is not shipped.
  - **P4-9 / P4-9a** R-B1 (equilibrium bridges) and R-B2 (stranded basins fed only through
    vanishing sources) ship as conservative rejections: the face condition is sound, not
    exact; spec sections 7, 8, 15 and 17 now say "sound; conservative on R-B1 and R-B2,
    listed".
    Task 0(a) found R-B1 on 332 faces and R-B2 on 21, none in the census or the `:OnlyI`
    variants. In enumeration (P4-9b/c): the tagging move drops <N> candidates as `:pumped`;
    the oracles find <F> of them well-posed (R-B2 false alarms), at <G> parents
    <: parent ids>. ROADMAP's Idea "Admit the equilibrium bridges …" cites this.
  - **P4-11 / P4-11a** the mod-p residual (about 300 lines) ships; it is the largest part
    that could go: deleting it leaves enumeration unchanged (every enumerated residual face
    rejects: census 77/77, ter-ter 1,296/1,296) and refuses hand-written lumpable splits.
  - **P4-12 / P4-12a** no mixed corners: 357 of 65,423 census mechanisms pass the guard yet
    are 0/0 where one substrate and every product are zero (264 with no limit, 93
    removable); 12 allosteric ping-pongs are non-unique there. Fits are unaffected. ROADMAP
    known issue added.
  - **P4-25** ε-path window: error ∝ ε^p, p a multiple of ½, within 10% per decade over
    1e-9…1e-18 (the spec's 1e-3…1e-12 fails 11.5% of interior directions).
  - **P4-29 / P4-29a / P4-29b** removed test subjects, for approval: the split's allosteric
    pre-check and balancer branch (TME "a split the boundary check would reject gains
    nothing", "a rejected split is extended by every :OnlyA group", "the split's :OnlyA test
    agrees with _check_well_posed": the first two parents are now `_NO_LIMIT` rejection
    tests and exit tests with exact surviving children; the third testset, an aggregate over
    the bi-bi seed set with no fixture of its own, became a test of P4-14's claim, that a
    split keeps its parent's interior verdict, rather than an exit test, which P4-29a asked
    for); `_has_strict_positive_combination`'s helper
    test (the helper is gone); the dead-inactive ping-pong derivation (AGT, now a rejection
    test pointing at ROADMAP's both-halves entry, plus half 2); the one-flip-segment message
    tests (their fixtures keep verdict tests).

  ## Rulings P4-1 to P4-31

  - P4-1 ping-pong stays allosteric: one half-reaction per V-type base; P4-2 both halves wait
    for iso mechanisms (ROADMAP); P4-3 `docs/src/roadmap.md` gains the limitation.
  - P4-4 class 3 is procedure A on the combined law network, flip at E; P4-5 a second closed
    class rejects (law and ε networks); P4-6 class 2 runs on the ε-network, once; P4-7
    admission needs the deletion limit along every joint path; P4-8 (above); P4-8b joint R1
    is tested on a mechanism that realizes VC X4's pattern (below) and, for its in-tree
    arithmetic, on a hand-built `_Network`.
  - P4-9, 9a, 9b, 9c (above); P4-10 a bottomless face is rejected, fed or not; P4-11, 11a
    (above); P4-12, 12a (above); P4-13 class 3 inside `_minimal_flips`, every exit drops;
    P4-14 the split loses its pre-check; P4-15 class 5 on the derived law; P4-16 the `:OnlyI`
    rule lives in the constructor; P4-17 `:OnlyI` mirrors `:OnlyA` (13 fixtures exact to
    1e-20); P4-18 turnover is structural.
  - P4-19 / P4-19a the tagging move's details and counts (below); P4-20 change-tag filters
    by admission, its all-or-nothing docstring rewritten; P4-21 flip at E, L at E documented;
    P4-22 a fix is printed only where a test shows its variant passes; P4-23 the guard's API;
    P4-24 / P4-24a the oracles in the fixtures file with Task 0's settings (every corner up
    to 128, κ ∈ {RE first, 8, 2, ½, ε first}); P4-25 (above); P4-26 the certified direction
    is the test's own; P4-27 the two trap specs become rejection tests with variants; P4-28
    population tests add <t> s (≤ 30 s); P4-29, 29a, 29b (above); P4-30 / P4-30a Task 0 ran
    before the plan; P4-31 Task 0's results stand, its sampled ε-path accepted.

  ## Where the plan departs from the rulings

  - **P4-31** follows research section 4's order but for two changes: the trap specs'
    task runs before the wiring, whose macros would throw on them, and the `:pumped`
    check is a task of its own (P4-9c).
  - **P4-9c** asked for the `:pumped` check on the package's own tagging move, before the
    move exists. The check ran the package's admission on research section 1.9's
    candidates, rebuilt apart from the move, and the move, once written, gives the same
    children on all <6,797> parents, which also tests it.
  - **P4-8b**: a mechanism realizes VC X4's pattern (a ping-pong with half 2 `:OnlyA` and
    the M bindings at E(A) and the covalent form `:OnlyI`): joint R1 alone rejects it at
    zero products and at P = 0, and the ε-path does not converge there. It needs a
    catalytic `:OnlyI` group, so its test came with `:OnlyI`; the guard's own network type
    tests the in-tree arithmetic before then, on VC X4's literal costs.
  - **P4-24a**: oligomers are judged at one subunit; the wide set for an RE-speed verdict
    at n > 1 has no caller, so the oracle leaves it out.
  - **P6** (research section 3.3) runs the ε-path at the interior only, within P4-28's
    budget; the RE-speed oracle covers every face.
  - **P4-22** prints a fix only where a test shows it working, so some messages leave out
    part of research section 1.11 (open questions).
  - **P4-29a**: one split testset became a verdict aggregate, not an exit test (above).

  ## Verdicts that changed

  - Rejected now, admitted before: an allosteric ping-pong with both half-reactions `:OnlyA`
    (an island; the TT, AGT and TME fixtures become rejection tests with a half-2 variant);
    AD1 and AD11 (`:no_limit`); the parents of TME's split pre-check testsets (`:no_limit`);
    the two trap specs (class 3 trap and split; their variants replace them in
    `MECHANISM_TEST_SPECS`).
  - Admitted now, rejected before: R1 and the other mechanisms whose second free form lies
    outside E's rapid-equilibrium segment (the flip at E); AD15.
  - Rejected before and now, by another check: R2 (`:pumped`; before, the one-flip check),
    R3 and R4 (`:island`), S1 (class 3, a second closed class), the bottomless-segment
    testsets (class 3).
  - "All-RE catalytic cycle raises" when the macro's code runs (class 2), not at derivation.

  ## Tagging move (P4-1, P4-19a)

  On Task 0(b)'s parents the package move gives <88,434> children on <88,174> networks
  (pingpong <4,878> = <4,358> K + <520> V; bi-bi <15,479>; ter-bi <68,077>), equal to the
  candidates built independently in Task 11 on every parent; class 3 drops <820> networks at
  the exit (pingpong <56>, ter-bi <764>); K-types on today's V-type set give <198> bi-bi
  children at <12> parents; each of pingpong's <30> regulated ping-pong parents allowing
  n = 2 keeps 2-4 V-types; <98> of <757> ter-bi ping-pong parents end with no child (<22>
  without a binding group, <76> losing every K-type to `:island` or `:pumped`), by design:
  today's children there are all islands, and a declared regulator would give each a
  V-type. Change-tag relaxes both half-reactions of a one-half V-type together; on the
  six-step bi-bi ping-pong, whose second half-reaction is rapid-equilibrium chemistry, that
  child passes admission and fails class 3 at zero substrates (a speed-dependent split), so
  the exit drops it: a fully live inactive conformation is reachable only where class 3
  allows. Please confirm these counts.

  ## Open questions

  <each item of the plan's "Open for Denis" part B that the runs left open, one line each,
  with its task and the choice it asks for>

  ## Measurements

  - Census against the branch point (Task 8) and through Tasks 9, 10, 12: every difference
    explained by class (<Task 8's and Task 12's census summaries>).
  - Ter-ter `init_mechanisms`: <s> s (budget 150 s; branch point <s0> s); count <n>.
  - `test_mechanism_enumeration.jl` focused: <te> s (branch point <te0> s).
  - Compile budgets: init trace-compile <n>, first `rate_equation` <s> s, bi-bi to uni-uni
    reuse <r> (budgets 200, 6 s, 0.01).
  - Full suite: <count> tests, <time>, maxrss <kb> KB (branch point <kb0> KB).
  - R18 remark scan, derivation and identify workloads, Julia 1.12 and 1.13: 0 hoisted
    allocations with a bits store into a reference field.

  🤖 Generated with [Claude Code](https://claude.com/claude-code)
  ```

- [ ] **Step 13: Push and open the draft PR.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/guard_wt && git push -u origin graph-guard
  gh pr create --draft --base graph-mechanism --head graph-guard \
    --title "Reject mechanisms whose rate law is not the fast-equilibrium limit (phase 4)" \
    --body-file $S/phase45/pr_body.md
  ```
  Expected: the PR URL. The body ends with the "Generated with Claude Code" line and holds
  no session link. To change the body later, use
  `gh api repos/DenisTitovLab/EnzymeRates.jl/pulls/<n> -X PATCH -F body=@$S/phase45/pr_body.md`
  (`gh pr edit` fails on the Projects-classic deprecation). Never merge.

---

## Open for Denis

Task 14's PR body names everything here. Part A lists the rulings; part B lists each question
the plan leaves open, by task. "Settled" marks an item the plan as assembled answers; it stays
listed so Denis can overturn the answer.

### A. Rulings the PR body names

The six Denis may overturn first (RES §4 "For Denis", with RUL's amendments):

1. **P4-8**: admission's face condition is B v3; the neutral-bridge exemption (B v2b, about 90
   lines) does not ship. P4-8b's first branch applies: a mechanism realizes VC X4's pattern
   (`joint_r1_mech.log`), joint R1 alone rejects it at zero products and at P = 0, and the
   ε-path confirms both rejections. Its test needs a catalytic `:OnlyI` group, so it lands
   in Task 10, not Task 5 as P4-8b says; Task 5 tests the in-tree arithmetic on a `_Network`
   with VC X4's literal costs. The check stays and is tested on a mechanism.
2. **P4-9**: R-B1 (equilibrium bridges) and R-B2 (stranded basins fed only through vanishing
   sources) ship as conservative rejections, and the exact Freidlin–Wentzell fallback becomes
   a ROADMAP Idea (Task 5). P4-9a: every spec passage that asked for an exact face
   condition (section 7 "Admission" item 2, section 8, the phase-4 row of section 15,
   section 17 row 2) now says "sound; conservative on R-B1 and R-B2, listed" (Task 5).
   P4-9b/c: Task 11 counts the `:pumped` drops the oracles find well-posed (R-B2 false
   alarms in enumeration); the Idea cites the count and the PR body lists their parents.
   P4-9c asked for the check before the tagging move (RES §4 row 10) yet "on the package's
   own tagging move", which cannot both hold: Task 11 runs the package's admission on a
   rebuild of RES §1.9's candidates written apart from the move, and Task 12 compares the
   move with that rebuild on every parent (stage `move`), which tests the move as well.
3. **P4-11**: the mod 2^61−1 residual ships. P4-11a: at about 300 lines it is the largest part
   Denis could delete; enumeration would not change, and hand-written lumpable splits would
   be refused.
4. **P4-12**: faces follow the spec, no mixed corners. P4-12a: 357 of 65,423 census mechanisms
   pass the guard yet are 0/0 where one substrate and every product are zero (264 with no
   limit, 93 removable, d7a09031d5dd2fc6 among them), and 12 allosteric ping-pongs are
   non-unique there; fits are unaffected; ROADMAP's known issue says so (Task 2).
5. **P4-25**: the ε-path window is error ∝ ε^p, p a multiple of ½, within 10% per decade over
   ε = 1e-9…1e-18, replacing spec section 14's 1e-3…1e-12.
6. **P4-29**: test subjects that go with their code, for approval: the split's pre-check and
   balancers (TME:2853, 2936, 3074), `_has_strict_positive_combination`'s test, AGT:997's
   dead-inactive ping-pong derivation, the one-flip-segment message tests. P4-29a keeps each
   fixture as a rejection, exit or verdict test; P4-29b takes the half-2 variants. P4-29a
   asks for exit tests of the split pre-check's fixtures: two become such tests, and the
   third testset, which has no fixture, becomes a verdict aggregate (item 38).

Also for Denis, **P4-31**, which says the plan follows RES §4's order: Task 7 (RES row 11,
the trap specs) runs before the wiring (row 7), whose macros would throw on the trap
specs, and Task 11 is new (P4-9c).

The PR body also lists: P4-1 to P4-3 (Denis's own rulings on allosteric ping-pong; P4-3's
user-facing roadmap entry lands in Task 8); P4-4 to P4-31 as adopted; P4-19a's counts, which
Denis confirms (Tasks 11-12); P4-24a's oracle settings; P4-31 (Task 0 done; part (b)'s sampled
ε-path stands and the 16 h full run is not ordered); the verdicts that changed; version 0.10.0
(R19).

### B. Open questions, by task

Tasks 1-3:

1. The oracle verdict vectors, limit verdicts (`_LIMIT_SPEED`, `_LIMIT_DIFFER`,
   `_LIMIT_EQUAL`) and message texts pinned in Tasks 1-3 were derived by hand, not run: RES
   §3.2's source logs (RD/hand.tsv, rules/vi) were lost in the 2026-10-14 outage, and the
   Julia lock was busy. Hand-derived: the ldh trap and split structure, feeder edge, messages
   and face order; the full per-face fault vectors of ldh_like_trap, trap_flip_child9,
   split_flip_child1 and bottomless_entered, and `_LIMIT_SPEED` for the ldh trap (procedure
   A by hand); that beside, flip2 and allre diverge on every face; sib's failure at R = 0;
   R2's network counts, edge order and stored directions; the bottomless, non-unique and
   split messages; the grid's 76/108. Not derived: lumpable_532's 14 faces, C3b's and C3c's
   stages, and whether d7a09031d5dd2fc6 passes classes 3-5. Measured since:
   d7a09031d5dd2fc6's law is NaN at B = P = Q = 0 and 0 at A = P = Q = 0 on the S3f
   snapshot (`mixed_corner_d7a0.log`) and at 08c4043, and the refresh pass (P4-32) ran
   Tasks 1-3's code at 08c4043, where every testset of the three tasks passed, these pins
   and d7a09031d5dd2fc6's classes 1-3 among them (`$S/phase45/refresh/map-head-t0-3_6.log`,
   `_7.log`, `_8.log`). The green steps of Tasks 1-3 list the rest as their
   first checks; each testset's oracle-agreement assertion is the value check, and a
   disagreement goes to systematic debugging (Task 3 Step 6), never a re-pin.
2. Settled: the BigInt keyword defaults. Every task uses `BigFloat` defaults (interface
   decisions).
3. Settled (P4-33): at 08c4043 `_testhelper_graph` (fixtures file :103-159) keeps every
   field the extended `_testhelper_consistent_point` and `_testhelper_steady_state_rate`
   read (`free`; `steps[j].names`, `inactive_names`, `group`, `from`, `to`, `consumed`,
   `released`; `inactive`, `tags`, `sites`), and the `in_law = allo ? ...` line (:131) and
   `inactive` field (:155) that Task 10 Step 1 (a) extends
   (`$S/phase45/refresh/settle_1.log`). Task 1's code runs on it
   (`$S/phase45/refresh/map-head-t0-3_6.log`: 159 pass).
4. `_faces_pass(mech, re)` runs class 2 on the ε-network before class 3: more than class 3
   alone, harmless, and it shields `_face_structure` from a fast turnover cycle. Settled: an
   allosteric parent's candidates reach it through `_minimal_flips(admissible, m::_AnyMech)`
   (Task 8).
5. Settled, with a choice for Denis: class 5 inside `_num_den_exprs` makes the package refuse
   0/0 laws, so the oracle rows of those fixtures read `law = :limit` (Task 8 Step 3 (l):
   entered, unfed, dead end, detour; `ss_dead_end_at_face`'s `:law_nan` row,
   unfed_bottomless's check). The alternative is class 5 in the guard alone, which would
   leave `_num_den_exprs` total but let a hand-built `rate_equation` call compile a 0/0 law.
6. The fix variants the full guard is asserted to admit were reasoned, not run: class 2's
   one step steady-state (re_turnover_beside_ss, re_turnover_flip2_453, m_allre), the trap
   feeder steady-state (ldh_like_trap, trap_flip_child9), and the bottomless fix
   `E(Q) + P <--> E(P, Q)`. Task 0 oracle-checked only the two trap specs' variants. If one
   fails a later class, P4-22 drops or rewords that fix.
7. Settled: the class-2 fix ("Make one step of the cycle steady-state") is printed only
   when any one step of the reported cycle, made steady-state with its group, leaves no fast
   turnover cycle (`_any_step_cures`, Task 2). TRD's m_sib, whose two fast cycles survive
   some single steps (and whose both-chemistry variant fails class 3 at R = 0), gets the
   cause alone. The check reads class 2 only, not the whole guard.
8. Settled (P4-33): S4 gave `show` a per-step renderer, `_step_text(m, k)` (names.jl:65),
   which `_print_steps` prints (dsl.jl:913-916). In the stored direction it writes each side
   as the form (`_form_text`, names.jl:41: a free form bare, `E`), then the metabolites the
   step takes up or gives off there in table order (`_met_text`, names.jl:34: `P::Mod` for a
   reactant's copy), the sides joined by `⇌` or `<-->`: the plan's `_step_line`, with
   `_met_text` its `_entry_text` and `_form_text` its `_form_line` (0 of 826 edges differ,
   all read back: `$S/phase45/refresh/map-head-t0-3_2.log`,
   `$S/phase45/refresh/settle_1.log`). The messages call `_form_text` and `_step_text`
   directly, so Task 2 writes none of the three and leaves `src/names.jl` alone; the pinned
   strings hold. Task 2 Step 0 allowed `_step_line` as a forward to the renderer; Denis may
   prefer the two one-line forwards.
9. Settled: P4-24a's RE-speed wide set at n > 1 has no caller (every oligomer is judged at
   n = 1, as Task 0 judged them), so `_testhelper_respeed` leaves the `protomers` branch
   out (YAGNI). The n-protomer oracle stays for Task 9's iso families.
10. `_faces_pass` seeds the limit test with the parent's content hash under any RE override,
    so a candidate and the child it becomes can draw different residue points (chance of a
    different verdict at most 2·deg/p). Accept this, or seed with `hash(re, hash(m))`?
11. The class-3 port allocates per face. If Task 3's probe reports much more than the
    prototype's 6.6 µs per face (target at most 3 µs), Task 8's 150 s ter-ter gate is at
    risk, and the structural part needs an allocation-free rewrite before the wiring.

Tasks 4-7:

12. The verdicts Tasks 4-7 pin (`_testhelper_epath_verdict`, `_testhelper_respeed_verdicts`,
    `_testhelper_admission_agreement`) come from the research oracles (probe p2_verdicts);
    the ported oracles first run in Task 1. Which of `:direction_dependent` and
    `:wrong_limit` the five-order ε-path returns for each island is unmeasured; the tests
    accept either.
13. Suite-time gates are unmeasured: Task 4 at most 1 s, Task 5 at most 3 s (8 s with
    corners), Task 6 at most 1 s. Task 5's "admitted tag sets converge" runs both oracles on
    39 mechanisms (about 190 faces) and may exceed 3 s; the plan gives no fallback.
14. The `_EPS_NONUNIQUE` message follows RES §1.11 ("reached only through tagged steps"). For
    S1's unfed active F that is inaccurate, but the guard reports class 3 first there; only
    `_admission` shows it.
15. The no-limit fix's two variants are now tested on AD11 and the ter cube too
    (`nolimit_fixes.log`, Task 4), and the arc-cycle fix's on the ordered bi-uni. The
    messages still print them for every shape of their class, reading P4-22 per class;
    checking each variant at the message would build a re-tagged mechanism, which the guard
    never does. Accept the per-class reading?
16. Settled by running the gates: JET on `_guard_message`'s narrowing, the R18 scan and the
    class-3 parts of every gate could not be checked before Tasks 1-3 existed.
17. Settled: `_testhelper_tag_sets`' docstring named `_testhelper_face_rejected`, which only
    Task 5 defines; the docstring no longer names it.
18. Settled: Tasks 4-6 edit Task 2-3 code by anchors (`_guard(am)`, `_guard_message`'s
    `_RE_TURNOVER` line, the `@enum _Fault` block and its docstring, `_bottom_sccs`) that
    Tasks 2-3 write as quoted; Task 5 reuses Task 3's `_bottom_sccs` through
    `_strong_components`.
19. Settled: the Task 1 contracts Tasks 4-7 assumed hold as Task 1 writes them:
    `@_testhelper_unguarded @allosteric_mechanism begin ... end` builds unguarded;
    `_testhelper_respeed_verdicts` lifts a compiled mechanism for `law = :package` and
    judges an oligomer at one protomer; `_testhelper_epath_verdict` returns
    `:converges`, `:no_path`, `:direction_dependent`, `:wrong_limit`, `:del_nonunique`,
    `:nonunique` or `:diverges`.
20. Settled: admission's turnover check (`_NO_ACTIVE_TURNOVER`, its message and its place in
    `_ADMISSION_FAULTS`) lands in Task 10, the first task whose constructor accepts a
    catalytic `:OnlyI` group, after its essential-activator test is watched failing (TDD).
21. RES §1.11's `:arc_cycle` fix ("untag one of them") fails on the smallest arc cycle (the
    ordered bi-uni with both bindings `:OnlyA`). The message offers the two fixes the probe
    confirms: untag them all, or tag the steps that join them. Confirm the wording (P4-22).
22. Settled: two more spec passages claimed exactness, section 7 "Admission" item 2 (lines
    392-394: "ships it only when it matches the RE-speed and ε-path oracles") and section 15's
    phase-4 row (748). P4-9a names lines 478 and 813; Task 5 amends all four the same way,
    so the spec stays consistent. Denis may prefer the two extra edits left out.
23. AD6 cannot be built after phase 3 (one reaction in two groups,
    `_check_each_reaction_once`), so the plan drops it and tests the duplicates
    AD2a_iso_chemNonequal (R1), AD2e (R2) and AD2c (R3) once each. Approve dropping AD6.
24. fz101_10232's faces 5 and 6 (P or S zero, with I zero) are pinned `_PUMPED` without an
    oracle verdict: the two-potential ε-path converges there and corners were not run. R-B2
    false alarms or correct rejections? Task 11's measurement may settle it.
25. The LDH variant's `rate_equation` takes 65.4 ns locally on phase 2, against the 120 ns
    gate and the original spec's 20.1 ns, so it has less CI margin than most specs. If CI
    timing fails, P4-27's other variant (only the feeding step split out, 7 groups) is cheaper
    but changes the ruling.

Tasks 8-10:

26. Settled (P4-33): at 08c4043 `_has_vmax(m::Mechanism, re::BitVector,
    names::Vector{Symbol})` (mechanism_enumeration.jl:749) and
    `_chemistry_equilibrates_both_sides(m::Mechanism, re::BitVector)` (:781) take flags over
    every edge and no mask, `_has_vmax`'s text as at proto3 2395401. The chemistry test's
    text differs: S3i (1182f84) moved its merged complexes into `_merged_complexes(m, subs)`
    (:758), which `_seed_candidate_screen` (:834) and `_seed_variants` (:903) share, and no
    `_fused_substrate_binding` remains. Task 10 gives all three an edge mask `on`,
    defaulting to every edge (Step 4 (f)), since folding the mask into the flags would not
    restrict the maximal-rate test, and a test pins the difference.
27. Settled (P4-33): the parser is `EnzymeRates._parse_mechanism(text::AbstractString)`
    (dsl.jl:935), which reads the text with the macros' walker and checks it through
    `_checked_mechanism`, so the wired guard's message reaches it. `show` of a concrete
    `AllostericMechanism` (dsl.jl:872) prints the `@allosteric_mechanism` text, tags and
    sites included, and reads back as the mechanism: Task 8's both-halves ping-pong `both`,
    built unguarded, gives `_parse_mechanism(sprint(show, both)) == both` at 08c4043, and
    all 11 allosteric specs read back (`$S/phase45/refresh/settle_1.log`); TRG "a graph's
    DSL text reads back as the graph" pins the same for 944 census allosteric mechanisms.
28. What ran before execution: `p2check_rules.jl` loaded Task 2's code, Task 3's
    `_closed_classes` and the code of Tasks 4-6 into the S3f snapshot and checked, without
    class 3 or the oracles, every admission verdict, named edge, per-face admission fault
    and message Tasks 4-6 pin (461 tests); the research rules, both research oracles and the
    S3f probes gave the expected verdicts. What did not run: Task 1's oracles, all of Task 3
    (procedure A, the limit test mod p, the class-3 messages), Task 2's own testsets, and
    Tasks 7-14's code. Pins that depend on that code (copy_group_bridge's singular "step …
    bridges" message, `_Verdict.conformation` values, the 716/4 split pin measured on phase
    2's moves) rest on the red and green runs; Tasks 1-3's green steps list theirs.
29. A modeling fact: P4-29b's half-2 variants, and P4-1's half-2 V-types in Task 12, predict
    a rate of exactly 0 whenever every product is zero and L > 0: the inactive conformation
    holds the enzyme in the covalent forms, since half 2 cannot run there and half 1 runs back
    only with the first product. The mass-action limit is 0 too, so the guard admits them
    correctly; the AGT "no covalent sink" checks become `== 0` pins (`p3c_sink.log`). A fit to
    initial-rate data penalizes such a variant but does not refuse it.
30. `seed_mechanisms` now filters seeds after its "no mechanism binds every required
    regulator" check; if the exit drops every seed it returns an empty list without an error.
    Not addressed; Task 12's tagging makes it unlikely.
31. Suite-time budgets are unmeasured: Task 10's at most 2 s covers RE-speed checks on all 7
    B1 fixtures; Task 9 runs 19 law checks; S1 and the TT halves get per-face admission
    agreement. Task 10 Step 5 reports an overrun rather than drop fixtures.
32. `_degenerate` reads the active law only (AM A8). On Task 10's admitted activator route the
    active law has a maximal rate both ways while the inactive law's product release closes a
    turnover cycle (`p3_degenerate.log`), so a law the inactive conformation dominates at high
    product never saturates. Only hand-written mechanisms reach this. Read both conformations'
    laws (degenerate if either is)? Its census effect is unmeasured.
33. The `_ISLAND` message's roadmap pointer for a ping-pong whose inactive conformation runs
    neither half-reaction ("is not supported yet (see the Roadmap page)", P4-29a): its
    wording and target (the user-facing roadmap) are a guess.
34. After Task 8 the RE-chemistry ground-truth harness in `test/allosteric_ground_truth.jl`
    gates nothing: its subjects (the RE-cycle error, the dead-inactive ping-pong's
    `d_free_I == 1`) go (P4-29, P4-29a/b). Keep it or delete it?
35. The copy rule's per-conformation masks (Task 10) are pinned by regression only: no known
    fixture gets another verdict from them. A copy group whose sites the active law prunes
    would tell them apart.
36. `ss_mirror`'s rank pin reads the unchecked law (TME's `_testhelper_identifiable_rank`).
    The alternative moves the pin to an admitted variant (`E(A) + Q <--> E(A, Q)`
    steady-state: fitted 9, rank 8), which may not keep the testset's subject (the gauge fails
    because the mirror is steady-state).
37. The split aggregate pin (716, 4) was measured with phase 2's moves and the research rules;
    Task 8 re-measures it before pinning if the population differs, and Task 12 re-measures it
    on unfiltered tag sets.
38. P4-29a turns the split pre-check's fixtures into exit tests with exact surviving
    children. Two of the three testsets hold a fixture and become such tests (no child
    survives); the third is an aggregate over the bi-bi seed set with no fixture of its own,
    so Task 8 turns it into a test of P4-14's claim (a split keeps its parent's interior
    verdict). The PR body names this departure from P4-29a for approval.

Tasks 11-14:

39. Settled: Task 12's edits quote Task 8's text as Task 8 writes it (checked while assembling
    the plan), and Task 8's exit testset lives in `test/test_mechanism_enumeration.jl`, where
    Task 12 re-pins it. The census names agree across tasks: Tasks 8, 9, 10 and 12 use
    `census/env_guard`, and Task 12 compares `census/P4_T12` with Task 10's `census/P4_T10`.
40. The bi-bi raw-input aggregate ("every bi-bi init mechanism and child is what its raw input
    builds") may change under P4-1: Task 0(b) reports 176 bi-bi gains over levels 0-2 but not
    by level, so the plan gives a measure-then-re-pin procedure, not a number.
41. The split-verdict aggregate's new (accepted, rejected) counts can only be measured once
    the code exists; the ε-path check precedes the re-pin.
42. Task 12 (a)'s `_testhelper_param_deltas(m, kept) == fill(1, 3)` rests on reasoning (half 1
    tied with `:EqualAI`, half 2 excluded, so only L is added); no probe measured it.
43. Task 12 (d)'s verdicts (`_SPLIT`/`_LIMIT_SPEED` at A, B; the four binding relaxations
    well-posed) were probed only on FP2, which is isomorphic to the six-step ping-pong. Step 8's
    value check confirms them on the package.
44. Task 11 assumes `using SHA, TOML, Printf` resolve from @stdlib under `TestEnv.activate()`;
    Step 4's smoke run would expose a failure in seconds.
45. Settled (P4-33): yes. Change-tag cannot give the six-step bi-bi ping-pong a fully live
    inactive conformation: on its one-half V-type, relaxing both halves to `:NonequalAI`
    passes admission and fails class 3 at zero substrates (a speed-dependent split),
    because its second half-reaction is rapid-equilibrium chemistry. Task 12 Step 11 adds
    the ROADMAP Known issue "**Change-tag cannot give the six-step bi-bi ping-pong a live
    inactive conformation.**", which Step 13 commits; the PR body raises it too.
46. Settled (P4-33): yes. P4 pinned `length(children) >= 5400` (5,505 at phase 3, changed
    by P4-1 and the flips' class-3 check). Task 13 Step 6 (a) pins it exactly to the count
    Step 5 measures, `@test length(kids) == <n>`, and the testset's comment calls it an
    aggregate regression pin over the seed set (CLAUDE.md enumeration rule 1's exception).
    At 08c4043 the population is 5,505 (`$S/phase45/refresh/map-t11-14_4.log`).
47. Settled (P4-33): yes. Task 13 Step 6 (b) reads the class-3 fault of
    lumpable_split_732 and unfed_trap from the guard, whose verdicts on them Task 3
    confirmed, and P7 pins `(_SPLIT, _LIMIT_DIFFER, 7)` and `(_SPLIT, _LIMIT_SPEED, 7)`:
    the refresh pass read these from Tasks 1-3's code at 08c4043, where every Task 3
    testset passed (`$S/phase45/refresh/apply_47.log`, `map-head-t0-3_7.log`).
48. The R18 scan's "no new hoist" compares with phase 3's final figure (0 hoisted
    allocations in its derivation and identify workloads at 67eaea5,
    `T63_remarks_scan.txt`, and in the GC probe's package scan at 4f9a34a); no scan of the
    branch point itself is in the inputs, so Task 14 scans the branch point only if a hoist
    appears.
49. Settled (P4-33): the refresh pass located each passage Task 14 edits at 08c4043, after
    S5's rewrite: Step 1's expectations and the old texts of Steps 2 and 5 quote S5's words
    (mwc_allostery.md's K-system paragraph; developer.md's "edge" columns sentence, split
    paragraph and `_is_chemistry` readers), and Step 5 gains (e) and (f) for the two
    developer.md passages S5 wrote that phase 4 makes false (the boundary-check paragraph
    and "the boundary check's `_bottomless_segment`").
50. P2 (the ε-path on every face of at least 30 allosteric members) may break P4-28's 30 s
    budget. Task 13 gives a shrink order (P2 to two faces, P5 to 40 children, P4 to every
    second seed, P6 to every second seed); Denis may prefer another.
51. Settled (P4-33): of the 113 mechanism blocks, texts and quoted blocks that S3h, S3i, S4
    and later commits added to the tests or moved under a new testset title
    (`$S/phase45/refresh/settle_6.log`), the phase-4 rules reject two that phase 3's check
    accepts (it could not read test_dsl.jl's read-back fixture with Greek names, which
    Tasks 1-3's guard admits: `$S/phase45/refresh/check-fidelity-1_1.log`), both converted
    in Step 3 (m): the DSL fixture already in `p3_scan.log`, whose
    testset S3i renamed, and the all-rapid-equilibrium merged base of TME "_seed_variants:
    the ordered/random seed" (S3i), converted with MERGE/ELIM. The exit also thins one pin
    S4 added, TRG's `length(allosteric) == 5 + 934 + 5` (the rules reject 90 ping-pong
    tagging children, `$S/phase45/refresh/settle_5.log`): Step 7 kind 4 re-pins it after a
    value check. Two S3g testsets the macro scan cannot see pin phase 3's bottomless
    message, TT "every derivation entry rejects a bottomless rapid-equilibrium segment" and
    TME "expand_mechanisms checks that every child's rate law is defined"
    (`$S/phase45/refresh/settle_4.log`); Step 7 names their conversions.

Across tasks:

52. The messages leave out these parts of RES §1.11, since no test shows the variant
    working (P4-22): `:nonunique`'s fix "give those forms a step that runs on the face";
    `:island`'s fix "tag a binding that leads to them" (the message offers only untagging
    the chemistry, and only where that alone joins the island to free E); `:bridge`'s
    "the forms they join" (the message names the steps); every fix for an `:OnlyI` step in
    the active conformation (`:no_limit`, `:arc_cycle`, `:island`); the constructor's
    `:OnlyI` fix "or bind a modifier"; and `:arc_cycle`'s "untag one of them" (item 21).
    Restore any of them with a test of its variant?
53. P6 (Task 12) runs the ε-path at the interior only, where RES §3.3 asks for both oracles
    on every face of each child; the RE-speed oracle covers every face. The narrowing keeps
    P4-28's budget. Run the ε-path on every face of a sampled subset instead?
