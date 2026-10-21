# Phase 5: law engine, unit rule, printer, key and identify cache — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Derive every rate law with one exact polynomial engine on phase 3's graph and
phase 4's masks, print it in the independent parameters with dimensionless denominator
terms, key it by a rename-proof canonical form, and let identify fit one law per key, count
keys in the beam and cross-validation, and pass fits on to one-subunit L phantoms.

**Architecture:** `_solve` returns the joint solve as one struct, `_Solve`, which the
engine, `all_parameters`, code generation and the printer read. `poly.jl` holds sparse
polynomials with checked `Int64` coefficients. `engine.jl` derives one conformation's law,
`_Law`, on an edge mask: rapid-equilibrium forest, spanning trees of the segment graph, the
edge-sum numerator, the solve substituted, the Haldane factor divided out exactly. It also
derives the MWC law, `_AlloLaw`, on `_active_mask` and `_inactive_mask`, under a cap of 5,000
on V×τ summed over both conformations. `codegen.jl` emits `rate_equation` from the law,
`key.jl` labels it canonically and certifies L phantoms, and `normalize.jl` and `print.jl`
print it. identify fits per key. The old engine (`POLY`, `sym_det`, `_num_den_exprs`, the
MWC `Expr` assembly) stays beside the new code until Task 13, so each task switches one
consumer; `eq_hash` goes in Task 8.

**Tech Stack:** Julia 1.10–1.13; EnzymeRates internals after phase 4 (graph `Mechanism`,
`_segments`, `_active_mask`/`_inactive_mask`, `_check_well_posed`); checked `Int64` and
`Rational{Int}`; `SHA` (already a dependency); `@generated` functions; Test, TestEnv, Aqua,
JET; CSV and DataFrames; Documenter with doctests.

**Spec:** `docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md` (read at
f9361cf in `$S/phase45/law_wt`; the fix wave left it unchanged),
sections 9, 10 and 11, the phase-5 parts of 13 and 14, the phase-5 row of 15, sections 16,
17 and 18, and the Appendix A cards (B6, B9, C1, C2, C4, C6, C10, E1, E8) and Appendix B
entries they cite.

`$S` = `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad`
(a symlink to `/home/denis.linux/claude-scratchpad/EnzymeRates-12f42995`).

Inputs every task reads:
- RES = `$S/phase45/research5/phase5_research.md`: §1 the designs in pseudo-code (the
  authority over research scripts), §2 the rulings, §3 the test plan, §4 the task order.
  Its reports (LAW, UNIT, PH, KC, VI, VC) sit beside it in `$S/phase45/research5/`.
- RUL = `$S/phase45/RULINGS.md` "Phase 5". P5-0 adopts P5-1 to P5-44 with these
  exceptions: P5-10a (N stays the spec's graded basis; no Cleland-first rule), P5-31a (the
  spec's pre-filter stays and the certificate runs inside it), P5-43a (research logs stay
  in the scratchpad), P5-45 (this plan; branch graph-law; PR base main once #94 merged,
  P5-48) and P5-46 (this refresh; graph-law cut at f9361cf; Task 0 merges graph-guard's
  head 991355f, which changes no `src` file, into it for Task 0's baseline, as P4-37(b)
  orders).
- SLC = `$S/slice/src`, the slice code the ports follow; R# = `$S/phase3/RULINGS.md` (R16,
  R18, R22, R37, R52).
- Code cites name functions and testsets as they stand at f9361cf, where graph-law is cut;
  graph-guard's head 991355f (draft PR #94), which Task 0 merges into graph-law, changes
  no `src` file after it. The refresh pass (P5-45, P5-46) re-anchored them there from the
  snapshot `$S/phase45/plan5/snap_a9e0699`; a figure measured on that snapshot keeps its
  commit.

## Task index

Tasks run in order. Each ends with the full suite green and pristine, then a commit. Task
numbers equal the rows of RES §4; Task 12b is the PASS1 half of RES row 12, split from the
printer so that a PASS1 compile surprise cannot hold back the printer's commit.

| Task | Title |
|---|---|
| 0 | Baseline on graph-guard's final head |
| 1 | Structured solve output (P5-4) |
| 2 | Exact polynomials, `poly.jl` (P5-1) |
| 3 | Plain engine, not wired; `_eq_complexity` by checked Bareiss (P5-3, P5-5) |
| 4 | MWC law on both masks; summed cap and filter (P5-8, P5-26) |
| 5 | Test helpers and class 5 on the new law (P5-6) |
| 6 | Code generation on the new law; the gates (P5-21, P5-22, P5-25, P5-39, P5-7) |
| 7 | The key (P5-18, P5-19, P5-21) |
| 8 | identify's key cache (P5-20, P5-27, P5-28) |
| 9 | Beam per key (P5-29) |
| 10 | L phantoms (P5-30 to P5-36, P5-31a) |
| 11 | The unit rule (P5-9 to P5-13, P5-10a) |
| 12 | The printer and coefficient form (P5-14 to P5-17) |
| 12b | PASS1 on the concrete law (P5-18) |
| 13 | Delete the old engine and printer |
| 14 | The size gate's constant on CI (P5-24) |
| 15 | Docs, ROADMAP, PR body, final GC gate, Julia 1.10 |

## Global Constraints

- Version stays 0.10.0 (`Project.toml` untouched; P5-42, R19).
- Julia compat floor 1.10 (`julia = "1.10"`); CI runs Julia 1.13 with the newest SciML;
  local runs use 1.12.5 (`julia`), 1.13.1 (`julia +1.13`) and 1.10.12 (`julia +1.10`). On
  1.13 a closure over a `@nospecialize` argument is untyped, so use `Base.Fix2` there.
- `rate_equation` allocates nothing and runs under 120 ns per call for every
  `MECHANISM_TEST_SPECS` mechanism (`_testhelper_test_rate_equation_performance`:
  `allocs == 0`, `t < 120e-9`). If the port allocates or misses a bound, the task stops
  for Denis.
- Size-scaled gate (P5-22 to P5-25): `t < max(120e-9, c * D)` on the six gate laws
  (`_testhelper_largest_enumerated_mechanisms`, `_testhelper_large_ter_bi_mechanisms`).
  - D is the law's denominator terms after substitution; an MWC law counts D_A + D_I, and a
    flip-only D_I counts 1. The pins [396, 396, 396] and [589, 799, 858] stay.
  - c (`_testhelper_seconds_per_denominator_term`) stays 0.8e-9 until Task 14 sets it by
    phase 1's rule: two CI runs on three runners, twice the worst ns per term, rounded up
    to one significant figure (RES predicts 1.0).
  - The helper times through a `@noinline` barrier (P5-25), once Task 6 confirms hoisting
    on package code.
  - Roots emit `sqrt(x)` for q = 2 and `x^(one(x)/q)` for q > 2.
- Expr-shape gates hold over every local of the generated body: binary `+` and `*` only,
  and greedy Horner form.
- Compile budgets (`test/test_compile_budget.jl`) hold: init trace-compile ≤ 200; first
  `rate_equation` ≤ 6.0 s; Sig encoding ≤ 250 instances and 100 KB; fit compile ≤ 250
  instances and 300 KB; 0 show instances in PASS1; from Task 12b, 0 method instances in
  PASS1 for a fresh mechanism. No budget is raised without Denis.
- R18: no function that callers loop over returns a `Union` of immutable,
  reference-holding parametric structs, or a `Tuple` whose element types differ by
  mechanism type. No loop reads a field whose type differs across a `Union`'s members.
- P5-1, P5-7: `_Poly`, `_Law` and `_AlloLaw` are `mutable struct`s with `const` fields and
  no type parameter. No phase-5 struct is an immutable wrapper of one reference.
- GC gate (P5-7): a typed workload outside the suite, in `synth/scan_f.jl`'s pattern, run
  on the package with `--pkgimages=no` and LICM remarks, read by
  `$S/phase3/gcprobe/hoistscan.py`, on Julia 1.12.5 and 1.13.1.
  - The workload runs the six merging loops (`merge_vec`, `merge_nothing`,
    `fresh_or_nothing`, `fresh_or_vec`, `fresh_active_or_poly`, `fresh_poly_or_vec`), with
    full GCs inside, plus a derivation workload and an identify workload.
  - It passes with 0 hoisted allocations with a bits store into a reference field.
  - Positive controls must show the hazard: the same loops on an immutable one-field
    wrapper defined in the script, and `$S/phase3/gcprobe/repro_upstream.jl`.
  - It runs in Tasks 2 (`_Poly` alone), 6 and 15.
- Arithmetic (P5-1, P5-2, P5-3): checked `Int64` only, with no `BigInt` retry.
  - An `OverflowError` while `_law` derives a law, plain or MWC, escapes as an internal
    error that names the mechanism. Every coefficient that grows does so inside the plain
    law `_law(m, s, mask)`, which catches it; the MWC assembly in `_law(am, s)` only
    divides contents (`_cross_weights`) and adds unit terms (`_site_weight`).
  - Code generation's products (`fI·DA`, `xA^n`), the normalizer's small rationals and
    the key's lowest terms let a raw `OverflowError` escape: under the cap they stay below
    τ_I·F·τ_A ≈ 4e8 (P5-2), so one would mean the cap failed. This narrows RES §1.3's
    sentence to the law and goes to Denis (Open for Denis, "Deviations from the research").
  - `_eq_complexity` returns `Int`; it takes the `BigInt` determinant only after an
    overflow and saturates at `typemax(Int)`.
- Cap (P5-3, P5-8): `MAX_RATE_EQUATION_TERMS = 5000` applies to V×τ(active) +
  V×τ(inactive); a flip-only inactive conformation counts 0. The check runs at the law's
  entry, before any polynomial product. `RateEquationTooLargeError` keeps its type and
  text.
- Filter (P5-26): `eq_complexity_filter` compares the same sum, with default 337.
- Masks (P5-8): every MWC use (law, `_eq_complexity`, filter, codegen, key parts) reads
  `_active_mask(am)` and `_inactive_mask(am)`, never `trues`.
- Key (P5-18 to P5-20): concentrations are labelled `c_<name>`; the key holds all 64 hex
  digits of the SHA-256. At 100 leaves the candidate's key is
  `bytes2hex(sha256("mechanism|" * mechanism_id))`, and identify's log names it.
  `mechanism_id`, `parent_id` and `fit_source_id` stay 16 hex digits (R16).
- Unit rule (P5-9 to P5-13, P5-10a): constant term; else the graded basis with fractional
  λ; else the ratio form. There is no Cleland-first rule.
- L phantoms (P5-30, P5-31a): a variant is a phantom when it has one catalytic subunit,
  `!any(_inactive_mask(am))`, every site ligand `:OnlyA`, and a certificate.
- Canonical Step Form stays load-bearing: no change to how `Step`, `Mechanism` or
  `AllostericMechanism` canonicalize.
- Parameter naming chokepoint: every constant name comes from `name(p, m)` or the solve's
  `names`. New code writes no `Symbol("K…")`, `Symbol("k…")`, `Symbol("V…")` or
  `Symbol("L…")` literal (the AST-walker testset in `test/test_types.jl`).
- Memory: one Julia process at a time (7.7 GB, no swap).
  - The full suite runs with `--heap-size-hint=2500M` and needs about 4 GB free.
  - Every full run reports maxrss; a rise over Task 0's figure is explained in the commit
    body.
  - Research and probe runs stay ≤ 3 GB, in chunks of at most 20 min.
- Suite time: Q1–Q6 add at most 20 s (P5-41) and Q7's subprocess about 10 s (P5-21).
  Ter-ter `init_mechanisms` stays under 150 s. Tasks 6, 12 and 13 report R52's three
  figures (full-suite maxrss, ter-ter init, `test_mechanism_enumeration.jl` focused)
  beside Task 0's.
- Runner. Every command runs from the worktree `$S/phase45/law_wt` (branch graph-law, cut
  at f9361cf; Task 0 merges graph-guard's head 991355f) through the session runner, one
  Julia process at a time. Each shell call starts fresh, so every block sets `S` first:
  `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh <log> <timeout_s> --project=. -e '<cmd>'`.
  - Focused `<cmd>` (timeout 1800): `using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<file>.jl")`.
  - Full `<cmd>` (timeout 3600): `using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])`.
  - Another Julia version goes first: `$S/phase3/jl.sh <log> <t> +1.13 --project=. -e '<cmd>'`.
  - Logs go to `$S/phase45/logs/P5_T<n>_<what>.log`.
  - A run longer than a few minutes goes to the background and is polled in a bounded
    loop: `timeout 580 bash -c 'until grep -q "^exit " LOG; do sleep 20; done'; tail -40 LOG`.
  - A focused run stops at the first top-level testset that fails or errors: the testset
    throws when it ends. A red run that expects more than one top-level testset to fail
    wraps the include in one outer testset, as `runtests.jl` does, so every testset runs:
    the focused `<cmd>` with `include("test/<file>.jl")` written
    `@testset "red" begin include("test/<file>.jl") end`. Its log ends with the outer
    testset's `Some tests did not pass` and `exit 1`. The rerun that watches the tests
    pass drops the wrapper.
- Measured figures. A `<t>`, `<kb>`, `<N>` or similar slot stands for the figure its step
  measures; the commit carries the number.
- Commits: a plain imperative summary like the repo's ("Reject the conformation labels A
  and I"), body optional, ending with
  `Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>` and
  `Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE`. Never skip a
  hook.
- Draft PR only, graph-law into main (P5-48); never merge. The controller pushes, opens
  the PR and runs CI (Task 14). The PR body ends with the "🤖 Generated with
  [Claude Code](https://claude.com/claude-code)" line and no session link. Edit it with
  `gh api repos/{owner}/{repo}/pulls/<n> -X PATCH -F body=@$S/phase45/pr_body5.md`, since
  `gh pr edit` fails.
- CLAUDE.md applies to every task:
  - TDD, with each failing test watched failing;
  - the `_testhelper_` prefix on every test helper, no mocks, pristine output;
  - ABOUTME headers on new files, 92-column lines;
  - the three enumeration-test rules (an inline macro fixture per testset; exact children
    with count and `Set`; chemistry written the enumerator's way);
  - ROADMAP updated in the commit that fixes, adds or retires an entry.
- Fixtures are inline DSL. A changed pin follows a value check (an oracle or the exact
  law), and its commit says why it changed (spec section 14; RES §3.4). A pin that fails
  its value check is a bug, not a re-pin. No test is deleted. A test whose subject code
  goes is converted, or listed for Denis in the PR body (Task 15).
- Prose (docs, docstrings, comments, commit and PR text) follows elements-of-style.

## Review Focus

These failure modes are the most likely to bite a user. The spec implies each one, and no
RES test exercises it. Each line names the task that adds its test.

1. A hand-written `:OnlyI` inhibitor (phase 4's new tag) printed by `rate_equation_string`
   in full and coefficient form. The printed law should evaluate to the derived law and to
   the exact limit, with no term of the `:OnlyI` binding in `D_A` and that term in `D_I`.
   The slice never derived `:OnlyI` (VC M3), and no RES test prints it. Test: Task 12.
2. A user reads the one-row `best_equation.csv` with `CSV.read(path, DataFrame)`. `key`
   should read back as its 64-hex `String`. Read with the documented
   `types = Dict(:mechanism_id => String, :parent_id => String, :fit_source_id => String)`,
   the three ids should come back unchanged, including an id of digits and `e` only, such
   as `"1234567e89012345"`. Test: Task 8; Task 10 adds `L_phantom`.
3. The 1-SE rule selects an L phantom whose source is absent from the CV pool. Its
   `best_equation.csv` row carries `L = 0.0`, `L_phantom = true` and the source's
   `fit_source_id`. `rate_equation(results.best, concs, params)` with that row's
   params should equal the source's rate at every point, finite. `rescale_parameter_values`
   should accept the params. Test: Task 10.
4. The fit of a phantom's source key throws, or returns a non-finite loss. When it throws,
   the phantom's row should be a failure row that carries the source's error; when the
   loss is not finite, the phantom's row carries that loss and retcode, as the source's
   does. Neither should enter the frontier or the CV pool, and the run should go on.
   Test: Task 10.
5. PASS1 runs on Distributed workers, so `key` and the CSV's `rate_equation` text come from
   other processes. A second process should give byte-identical keys and printed text, in
   full and coefficient form, for the 42 specs and the six gate laws. Test: Task 12b extends
   Task 6's subprocess testset.

---

## File Structure

Source:

| File | Change | Responsibility |
|---|---|---|
| `src/poly.jl` | create (T2) | Exact sparse polynomials over integer-indexed variables: `_Vars` (with `==` and `hash`), `_Exp`, `_Mono`, `_Poly` (mutable, checked `Int`), `_divexact`, `_subst`, lowest-terms shifts, `_split_conc`, `_content` |
| `src/engine.jl` | create (T3, T4); modify (T13) | The law of one conformation on a mask (`_Law`, `_bare_law`, `_law(m, s, mask)`), every edge (`_law(m, s)`, `_law(m)`, T4), the MWC law on both masks (`_AlloLaw`, `_flip_only_law`, `_flip_only`, `_cross_weights`, `_site_weight`), `_spanning_trees`, `_tree_count`, `_eq_complexity` and the cap (`RateEquationTooLargeError`, `_assert_derivable`, moved from `rate_eq_derivation.jl` in T3; `MAX_RATE_EQUATION_TERMS` moved here in T13) |
| `src/codegen.jl` | create (T6) | `_Program`, `_local!`, `_factor!`, `_column_values!`, `_state!`, `_program`, `_loc`, `_power`, greedy `_horner`, `_emit`, `_body`, `_rate_body`, the generated `rate_equation` and its concrete method (moved from `rate_eq_derivation.jl`) |
| `src/key.jl` | create (T7, T10) | `_Canonical` (with `leaves`), `_KEY_LEAF_CAP`, `_key_parts`, `_ranks`, `_refine`, `_encode`, `_KeySearch`, `_search!`, `_canonical`; the L-phantom certificate `_phantom_weights` and `_is_phantom` |
| `src/normalize.jl` | create (T11) | The unit rule: `_Norm`, `_term_key`, `_normalizer` per conformation |
| `src/print.jl` | create (T12) | `rate_equation_string` (moved from `rate_eq_derivation.jl`) with the `coefficients` keyword, `_rate_text`, `_law_text` (plain and MWC), the text helpers (`_power_text`, `_ratio_text`, `_term_text`, `_sum_text`, `_conc_label`), `_Term`, `_apart!`, `_printed`, `_group_sum`, `_label_times`, `_coefficient_sum`, `_v_line`, `_reference_lines`, `_header_lines`, `_conformation`, `_site_line!` |
| `src/thermo.jl` | modify (T1, T13) | `_Dependent`, `_Solve`, `_solve`; `_dependent_param_exprs` renders `_Solve`; `build_power_expr` moves here in T13 (exact exponents from T12) |
| `src/rate_eq_derivation.jl` | modify (T1, T3, T4, T6, T8, T12, T13) | Keeps the parameter API (`parameters`, `fitted_params`, `all_parameters` on `_Solve`) and the per-conformation names (`_edge_names`, `_state_mask`, `_kreg_params`, `_reg_K`, `_ss_rate_constant_names`). Loses `_evaluate_dependent` (T1), the cap code (T3), the generated `rate_equation`, `_destructuring_expr`, `_dep_assignments` and the `horner` keyword (T6), the printer's entry points (T12) and the old engine with `_equation_text` (T13). `_num_den_exprs(am)` runs the summed cap first (T4); `_equation_text`'s docstring drops `eq_hash` (T8). ABOUTME rewritten in T6 and T13. |
| `src/sym_poly_for_rate_eq_derivation.jl` | modify (T6, T12); delete (T13) | `_poly_to_horner` goes (T6); `build_power_expr` writes `Rational{Int}` exponents (T12) |
| `src/admission.jl` | modify (T5, T13) | Class 5 on `_Law`/`_AlloLaw` (`_law_branches`, `_keeps_monomial`, `_vanishing_face`, `_support_verdict(m, vars, branches)`, `_assert_support(m, L)`); the `POLY` methods deleted in T13 |
| `src/identify_rate_equation.jl` | modify (T4, T7, T8, T9, T10, T12b) | Filter docstring (T4); `using SHA: sha256` moves to `EnzymeRates.jl` (T7); `_Fit`, `_Candidate`, `_fit_order`, `_row_params`, `BatchEntry` (`key`, `fit_key`), memo by key, CSV `key` and `fit_source_id`, `_offer_cv!` per key, `_beam_search` returns `(pool, memo)`, folds per fit key, `_rate_eq_dedup_key` deleted (T8); `_select_count!` per key (T9); phantoms and `L_phantom` (T10); PASS1 on the concrete law, behind the guard and without compiling, its result `records` (T12b) |
| `src/mechanism_enumeration.jl` | modify (T10) | `_expand_to_allosteric`'s docstring: K-type L phantoms stay in the frontier and inherit a fit |
| `src/types.jl` | modify (T13) | The `_inactive_mask` and `_active_mask` docstrings name `_law` |
| `src/EnzymeRates.jl` | modify (T2, T3, T6, T7, T11, T12, T13) | Include order: `types`, `names`, `dsl`, (`sym_poly…` until T13), `rate_eq_derivation`, `thermo`, `poly`, `engine`, `normalize`, `print`, `key`, `codegen`, `admission`, `fitting`, `mechanism_enumeration`, `identify_rate_equation`; `using SHA: sha256` (T7); the four `precompile` lines at the end stay |

Tests:

| File | Change | Responsibility |
|---|---|---|
| `test/test_poly.jl` | create (T2) | Unit tests of `poly.jl` |
| `test/test_law.jl` | create (T3, T4) | The plain and MWC laws against the exact oracles; `_spanning_trees`, `_tree_count`, `_eq_complexity`; a conformation's law on its mask; cap, overflow and RE-turnover refusals; `_cross_weights`; the summed cap |
| `test/test_key.jl` | create (T7, T10) | Key toys, attacks, leaf cap, pinned digests; certificate unit tests |
| `test/test_print.jl` | create (T11, T12) | Unit-rule fixtures (F1–F4, MM, Segel, the power case, MWC per rule) and the 42 specs' units; printer and coefficient-form fixtures; the `:OnlyI` printed law |
| `test/runtests.jl` | modify (T2, T3, T7, T11) | Includes `test_poly.jl`, `test_law.jl`, `test_key.jl`, `test_print.jl` in that order after `test_rate_eq_derivation.jl` |
| `test/mechanism_definitions_for_test_enzyme_derivation.jl` | modify (T3–T7, T10–T13) | The oracle helpers of RES §3.1 (see Interface decisions), BigFloat inputs through the body among them (T6); the six gate laws (T6); the rank oracle (T10); the 11 factored pairs (T12) |
| `test/test_joint_solve.jl` | modify (T1) | `_Solve` structure tests |
| `test/test_accessors.jl` | modify (T1) | `_evaluate_dependent`'s testset becomes an `all_parameters` testset |
| `test/test_rate_eq_derivation.jl` | modify (T4, T6, T12, T13, T14) | Pins 5–9, 11–12, 15–17 of RES §3.4; the gates and the helper's barrier; the shuffled-NamedTuple, cube-root and zero-concentration tests; BigFloat inputs through the body on 15 specs (T6); c (T14) |
| `test/test_oracles.jl` | modify (T4, T5, T7, T10, T11, T12) | Population tests Q1–Q3, Q5–Q6 on `_testhelper_population_sample`; the sample's same-model regroupings, computed once (T7); the half-integer dependent's `Rational{Int}` exponent (T12); FW-23's testset "an admitted :OnlyI mechanism compiles to its law"; the ABOUTME FW-39 wrote |
| `test/test_admission.jl` | modify (T5) | Class 5 on `_Poly`; `_num_den_exprs` calls rebased; FW-07's seven verdict calls read `_testhelper_support_verdict` |
| `test/allosteric_ground_truth.jl` | modify (T13) | Internal pins on `_AlloLaw` fields |
| `test/test_allosteric_golden.jl`, `test/reference/allosteric_golden_reference.txt` | modify (T12) | Regenerated golden, both printed forms |
| `test/test_allosteric_collapse.jl` | run (T12; modify only if a line moves) | Mirror lines as reference lines |
| `test/test_mechanism_enumeration.jl` | modify (T4, T5, T7, T10, T13) | Term-cap testsets (pin 13); `_testhelper_identifiable_rank` on the new law (moved to the fixtures file in T10); dedup-key testsets on the key; phantom verdicts on the K-type fixtures; internal pins |
| `test/test_identify_rate_equation.jl` | modify (T4, T7, T8, T9, T10, T12b) | The filter fixture (T4); Q4 (T7); pin 18; same-key fit sharing; beam per key; phantom rows; the CSV round trip; PASS1's guard (T12b) |
| `test/test_compile_budget.jl` | modify (T6, T7, T8, T12b) | Q7's subprocess testset (bodies, keys, texts); the PASS1 workload (pin 19); `PASS1_INSTANCES` (T12b) |

Docs and records:

| File | Change | Responsibility |
|---|---|---|
| `docs/src/deriving/textbooks.md` | modify (T12 doctest, T15 prose) | Normalization, coefficient form |
| `docs/src/deriving/ping_pong.md` | modify (T15) | Normalization (graded N), coefficient form |
| `docs/src/deriving/mwc_allostery.md` | modify (T15) | Per-conformation names `Q_s`, `D_s`, `F_s`, `W_s`, `N_s`; the summed cap |
| `docs/src/deriving/iso.md`, `docs/src/deriving/cha_king_altman.md` | modify (T15) | A printed term's description; the segment root at a segment's first form, the cofactor sentence, the ratio-form caveat |
| `docs/src/identify/tutorial.md`, `model_selection.md`, `enumeration_engine.md` | modify (T15) | Key, inherited fits, L phantoms, CSV columns, `types = String` |
| `docs/src/developer.md` | modify (T6 cites, T14 c, T15) | Architecture of phase 5; gates with the new c |
| `README.md`, `docs/src/getting_started.md` | modify (T15) | Printed laws |
| `ROADMAP.md` | modify (T3, T4, T8, T10, T13, T15) | Each entry in the commit that settles it (see Tasks) |
| `docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md` | modify (T2, T10) | P5-1's amendment of section 9 "Engine"; P5-30's of section 13 "One-subunit L phantoms"; both listed for Denis |
| `.claude/CLAUDE.md` | modify (T14) | The size gate's constant (P5-24; ◆ Denis) |

Scratch, outside the repository:

| Path | Responsibility |
|---|---|
| `$S/phase45/p5_base/` | Task 0's baseline files and `figures.md` |
| `$S/phase45/scripts/P5_sets.jl`, `P5_baseline.jl`, `P5_compare.jl`, `P5_identify_run.jl` | The baseline's mechanism sets; baseline dump; comparison of a later tree against it (read through Task 6 Step 1, while the gate laws live in TR); the seeded identify run (T0, T9) |
| `$S/phase45/scripts/P5_r55_recount.jl`, `P5_T9_check.jl`, `P5_phantom_rerun.jl`, `P5_T12_pins.jl`, `P5_T12_golden.jl`, `P5_census_nostr.jl`, `P5_T14_c.jl`, `P5_mkenv110.jl` | R5.5 on P4-1 tagging (T4); the identify run's check (T9); R5.2 on P4-1's ping-pong children (T10); Task 12's texts with their value checks, the golden file, the census without strings; c from CI logs (T14); the Julia 1.10 environment `$S/phase45/env110` (T15) |
| `$S/phase45/gcscan5/{mkenv.jl, scan_law.jl, derive.jl, derive_control.jl, identify.jl}`, `env112`, `env113` | The GC gate (T2, T6, T15) |
| `$S/phase45/probes/P5_T<n>_*.jl` | Probes: head renderings (T1), live bytes (T3), cap fixture (T4), hoisting (T6), key cost (T7), same-key values (T8), phantom verdicts (T10), old and new denominators (T13) |
| `$S/phase3/census/env_law`, `P5_base`, `P5_T6`, `P5_T12`, `P5_T13` | Census environment and corpora (`census/README.md`; adapter `neutral_graph.jl`); `compare.jl` runs from `census/env` with `rename_phase3.toml`, as phase 4 ran it |
| `$S/phase45/pr_body5.md` | The draft PR's body (T14, T15) |
| `$S/phase45/logs/` | Every run's log |

---

## Interface decisions that span tasks

- Names. Internal types and functions carry a leading `_`, as `_Column` and `_Verdict` do.
  `_eq_complexity` is the one complexity function: Task 3 rewrites it, Task 4 sums the
  masks, and there is no alias (RES §1.11's `complexity`).
- Laws. `_law(m, s, mask; inactive, concentrations)` derives one conformation (Task 3);
  `_law(m::Mechanism, s)` every edge and `_law(am, s)` the MWC law (Task 4), so
  `_law(x, _solve(x))` reads both mechanism types; `_law(x)` solves first. `_rate_body`
  (Task 6), PASS1 (Task 8) and the printer (Task 12) call `_law(m, s)`.
- Variables. Every `_Law` derived from one `_Solve` numbers its `_Vars` alike:
  concentrations (`concentrations`, default `metabolites(m)`; `metabolites(am)` for both
  conformations of an MWC law), then `s.names[i]` for `i in s.independents`, then `:Keq`.
  So the two conformations' polynomials multiply directly. `_AlloLaw.vars` appends the
  regulator constants (`s.fitted` order) and `:L`. Segment 1 holds free E in both masks
  (`_segments` roots segments in form order), so `L.trees[1] * L.shift` is free E's weight.
- Equality. `_Vars`, `_Poly`, `_Law` (Task 3) and `_AlloLaw` (Task 4) define `==` and
  `hash` over their fields (RES §1.3; a mutable struct's default `==` is identity). Task 3's
  rewrite tests compare laws with `==`, and Task 4 checks that two derivations of one MWC
  law are `==` with equal hashes.
- Arithmetic. `_law` turns an `OverflowError` into an internal error that names the
  mechanism (Task 3; the MWC law's growing arithmetic runs inside its two plain laws).
  Code generation, `_normalizer` and the key let one escape, which P5-2's bound rules out
  (Global Constraints). `_phantom_weights` returns no certificate, so the variant is
  fitted.
- Class 5 runs where a law is derived for use, never inside `_law`: the generated body
  (T6), PASS1 (T8; from T12b after `_check_well_posed`) and the printer (T12). The oracles
  can then read the law of a mechanism the guard rejects (P4-24).
- The old engine stays live until Task 13. From Task 4, `_num_den_exprs(am)` calls
  `_assert_derivable(am)` first, so every public path refuses a law over the summed cap.
  Task 6 deletes the `horner` keyword of both `_num_den_exprs` methods and
  `_poly_to_horner`; `eq_hash` and `_rate_eq_dedup_key` go in Task 8.
- R18 care. Call `_solve`, `_law` and `_assert_support` separately. No function returns
  `(s, L)` as a tuple, whose type would differ between plain and MWC mechanisms.
  `Union{_Law, _AlloLaw}` arguments read only `vars` and `H`, which share one type.
- Mechanism and law together. A function that derives `L = _law(m, s)` and passes `m` and
  `L` to methods defined per pair (`_program`, `_law_text`) has one method per mechanism
  type (`_rate_body`, Task 6; `_rate_text`, Task 12), so `L` is concrete in each. On a
  `Union` argument inference splits `m` and `L` apart, and JET reports the two mismatched
  pairs (`$S/phase45/plan5/probes/review_exec/jet_union.log`). PASS1, whose `m` is an
  element of a `Union` vector, narrows both as Task 10's phantom line does: `m isa
  AllostericMechanism ? f(m, s, L::_AlloLaw) : f(m, s, L::_Law)` (Task 12b).
- The key. `_Canonical(order, text, key, leaves, capped)` (Task 7); `leaves` is the
  labellings compared, `_KEY_LEAF_CAP` when capped. `using SHA: sha256` sits in
  `src/EnzymeRates.jl` from Task 7, so tests read `EnzymeRates.sha256`.
- identify's records (Tasks 8–12):
  - `_Fit(values, loss, retcode, peak_at_range_edge, source, source_id)`: `values` at
    canonical positions, then `L` for an MWC law.
  - `_Candidate(mech, fitted, mechanism, mechanism_id, eq_text, canon)`, Task 10 appending
    `phantom::Bool` and `canon0::_Canonical`: one concrete PASS1 record (R18), not a
    `NamedTuple`. PASS2's `pmap` returns a `_Fit` or the error's text.
  - `BatchEntry(mech, n_params, loss, retcode, key, fit_key, row)`; `fit_key` is `key`, or
    for an L phantom its L→0 key (Task 10).
  - A row keeps its own `fitted_params` order (`_row_params`). `fit_inherited` means the
    row is not its fit's source, and `L_phantom` is a required row field from Task 10.
  - `_select_count!`'s key argument is `law_keys` (Task 9). PASS1's result is `compiled`
    in Tasks 8–12 and `records` from Task 12b.
- Test helpers (shared fixtures file unless noted):
  - `_testhelper_poly_value(p::EnzymeRates._Poly, vars::EnzymeRates._Vars, env::AbstractDict{Symbol})::BigFloat` (T3); the two-argument `POLY` method goes in T13
  - `_testhelper_evaluate_law(L::EnzymeRates._Law, env::AbstractDict{Symbol})::BigFloat` (T3); the `_AlloLaw` method (T4). `env` holds the concentrations, independents, `Keq`, regulator constants and `L`.
  - `_testhelper_laplacian_tree_count(m::EnzymeRates.Mechanism, mask::BitVector)::BigInt`: V×τ from the test's own segment Laplacian (T3)
  - `_testhelper_detailed_balance(m, rng)::Vector{BigFloat}`, rebased on `_law` (T4)
  - `_testhelper_law_value(m, params::AbstractDict, concs::AbstractDict)::BigFloat`, rebased on `_law` (T5)
  - `_testhelper_support_verdict(mech)::EnzymeRates._Verdict` replaces FW-07's `_testhelper_law_branches(mech)` (T5)
  - `_testhelper_largest_enumerated_mechanisms`, `_testhelper_large_ter_bi_mechanisms`, moved unchanged from TR (T6)
  - `_testhelper_check_law_against_oracle(m; n_points = 3, seed = 1, bigfloat = false, rtol = big(10)^-20)`: with `bigfloat`, BigFloat inputs through the compiled body at `rtol` (T6, RES §3.1)
  - TR (T6): `_testhelper_rate_body(m)::Expr`, `_testhelper_body_locals(m)::Vector{Pair{Symbol, Any}}`, `_testhelper_generated_rate_return(m)`, `_testhelper_denominator_terms(m)::Int`, `_testhelper_rate_barrier` (when the hoisting probe finds hoists); test_compile_budget.jl: `_testhelper_subprocess_output(script::String)::String`
  - `_testhelper_same_law_under_map(L1, L2, order1::Vector{Symbol}, order2::Vector{Symbol}; rng, swap::Bool = false)::Bool`, `_testhelper_renamed_law(L, rename::Dict{Symbol, Symbol})`, `_testhelper_renamed_vars(vars, rename)::EnzymeRates._Vars` (T7)
  - `_testhelper_identifiable_rank(m; npts = 60, ndraws = 3, h = 1e-5)::Int`, deriving the law once per call (T5), moved from TME (T10)
  - `_testhelper_phantom_map_holds(am; rng)::Bool` (T10)
  - `_testhelper_name_units(m)::Dict{Symbol, Tuple{Rational{Int}, Rational{Int}}}`, `_testhelper_term_units(units::AbstractDict, vars::EnzymeRates._Vars, mo::EnzymeRates._Mono)::Tuple{Rational{Int}, Rational{Int}}` (the (concentration, time) exponents), `_testhelper_check_units(m, L)::Nothing`, `_testhelper_norm_value(nz, D, vars, env)::BigFloat` (T11)
  - `_testhelper_eval_rate_lines(s::AbstractString, params, concs)::Dict{Symbol, BigFloat}` and `_testhelper_eval_rate_string(s::AbstractString, params, concs)::BigFloat`, `params` and `concs` NamedTuples or dictionaries, in 256-bit, the latter moved from TR (T12)
  - `_testhelper_print_artifacts(s::AbstractString)::Vector{String}`, `_testhelper_printed_law_error(m; points::Int = 2, seed::Int = 31)::BigFloat` (T12)
  - `_testhelper_steady_state_rate` gains `occupancy::Int = 0` (T13)
- Placement. Population tests live in `test/test_oracles.jl`, beside
  `_testhelper_population_sample`. Enumeration-move assertions stay in
  `test/test_mechanism_enumeration.jl` under CLAUDE.md's rules. Gates stay in
  `test/test_rate_eq_derivation.jl`, which CLAUDE.md names. The subprocess test lives in
  `test/test_compile_budget.jl`, which already spawns processes. A helper two files read
  lives in the fixtures file, since a focused run sees only it and the file it names.

---

## Tasks

### Task 0: Baseline on graph-guard's final head

This task records today's engine and printer before any change, so later tasks can
compare against them outside the suite. It changes no repository file.

**Files:**
- None in the repository. The worktree `$S/phase45/law_wt` holds branch graph-law, cut
  at f9361cf (P5-46); Step 1 merges graph-guard's head 991355f, which changes no `src`
  file after f9361cf, into it (P4-37(b)).
- Create (scratch): `$S/phase45/scripts/P5_sets.jl`, `P5_baseline.jl`, `P5_compare.jl`,
  `P5_identify_run.jl`; the directory `$S/phase45/p5_base/` and its `figures.md`; the
  census environment `$S/phase3/census/env_law` and corpus `$S/phase3/census/P5_base`.
- Logs: `$S/phase45/logs/P5_T0_{inst,full,terter,tme,baseline,selfcheck_<kind>,mkenv,census,identify}.log`.

**Interfaces:**
- Consumes:
  - fixtures file: `MECHANISM_TEST_SPECS`, `_testhelper_lift`, `_testhelper_single_subunit`,
    `_testhelper_law_inputs`, `_testhelper_law_value`, `_testhelper_onlyi_fixtures`,
    `_testhelper_plain_moves`, `_testhelper_population`;
  - `test/test_oracles.jl`: `_testhelper_pingpong_allo_rxn`, `_testhelper_content_text`,
    `_testhelper_population_one_round`, `_testhelper_population_sample`,
    `_testhelper_inhibitor_rxn`;
  - `test/test_rate_eq_derivation.jl`: `_testhelper_largest_enumerated_mechanisms`,
    `_testhelper_large_ter_bi_mechanisms`,
    `_testhelper_random_independent_params_concs`,
    `_testhelper_test_rate_equation_performance`;
  - src: `_dependent_param_exprs`, `parameters`, `all_parameters`,
    `rate_equation_string`, `_rate_eq_dedup_key`, `_eq_complexity`, `_active_mask`,
    `_inactive_mask`, `_raw_symbolic_rate_polys`, `_mechanism_id`, `identify_rate_equation`;
  - `$S/phase45/scripts/P4_terter.jl`; census `mkenv.jl`, `corpus.jl`.
- Produces:
  ```julia
  # P5_sets.jl (included after the fixtures file, from the worktree root)
  p5_defined_name(ex)::Union{Symbol, Nothing}
  p5_eval_definitions(path::AbstractString, wanted::Vector{Symbol})::Nothing
  p5_id(m)::String                                   # _mechanism_id(sprint(show, m))
  p5_sets()::Vector{Tuple{String, String, Any}}      # (set, label, concrete mechanism)
  p5_counts(sets)::Vector{Pair{String, Int}}
  p5_try(f)::String                                  # escape_string of f() or "error: …"
  ```
  - `julia --project=. P5_baseline.jl <out_dir>` writes, for the sets `spec` (42),
    `sample` (120), `bibi_init` (239), `onlyi` (13), `gate` (6) and `inhibitor` (60),
    seven tab-separated files with a header row:
    - `laws.tsv`: `set id point n params concs value`. There are three points per
      mechanism, seeded by its id. `n` is 1 (one catalytic subunit) and, for an
      allosteric mechanism with more, its own multiplicity. `params` and `concs` are
      `name=value` pairs joined by `;`, and every number is a 256-bit decimal.
    - `dependents.tsv`: `set id dep indep parameters p_values all_keys all_values`. `dep`
      is the escaped `repr` of `_dependent_param_exprs(m)[1]` sorted by name. `p_values`
      holds the Float64 inputs, in `(indep..., :Keq)` order, at point 1. `all_keys` lists
      `all_parameters`' keys joined by `,`, or the escaped error text, and `all_values`
      gives its values' `repr`s.
    - `strings.tsv` (`set id string`): the escaped `rate_equation_string`, or
      `error: …`.
    - `eqhash.tsv` (`id eq_hash`): the 239 bi-bi seeds.
    - `complexity.tsv` (`set id all active inactive`): `_eq_complexity(m)`; for an
      allosteric mechanism also per mask.
    - `dterms.tsv` (`label id denominator_terms`): the six gate laws.
    - `timing.tsv` (`kind label allocs ns_per_call`): 42 specs and 6 gate laws.
  - `julia --project=. P5_compare.jl <base_dir> <kind>` rebuilds the sets, compares the
    current tree with the baseline and exits 1 past the kind's tolerance:
    - `laws`: 1e-20 relative, through `EnzymeRates._law` and `_testhelper_evaluate_law`
      when the tree defines both for the mechanism, else through `_testhelper_law_value`;
      it prints the path per set;
    - `dependents`: byte equality of `dep`, `indep` and `parameters`, and `all_parameters`
      to 1e-14 relative;
    - `strings`, `complexity`, `dterms`: byte equality. `dterms` counts
      `_nterms(_law(m).D)` once `_law` exists.
  - `julia --project=. P5_identify_run.jl <out_dir>` writes the run's CSVs and
    `summary.txt` (seconds, fits run, the selected `n_params` and text, the true text).
  - `$S/phase45/p5_base/figures.md` holds Task 0's figures. Tasks 1, 2, 3, 6, 12 and 13
    compare against it.

**Gate:**
- The full suite on the head is green and pristine.
- `p5_base/` holds the seven files, and the set counts read 42, 120, 239, 13, 6 and 60.
- `dterms.tsv` reads 396, 396, 396, 589, 799, 858, and `eqhash.tsv` holds 239 classes.
- Each kind of `P5_compare.jl` passes against the baseline on the unchanged tree.
- `census/P5_base` builds, and `meta.toml` has empty `sig_mismatches` and `conflations`.
- The identify run completes and writes `summary.txt`.
- Timings: the suite about 7 min, the baseline about 2 min, the self-checks about
  1 min each; the census (about 2 min) and the identify run (about 2 min: 102 s of
  search on f9361cf) run in the background.

**Testsets added:** none. **Pins changed:** none.

- [ ] **Step 1: Bring the worktree to graph-guard's final head.**

  The worktree is cut (P5-46): branch graph-law at f9361cf, with graph-guard's
  `Manifest.toml`, and at most the refreshed plan committed on top. graph-guard has since
  moved to 991355f (draft PR #94, whose base is main since #93 merged): a98c071 raises the
  constructor benchmark's gate in `test/test_reaction_graph.jl` to 20 µs, which CI failed
  at 10 µs on PRs #93 and #94; 3543abe and 7de259f name "LDH i-state NonequalAI 6-group,
  steady-state NAD group" (about 55 ns) as the slowest spec under the 120 ns gate, which
  the ubuntu and windows runners miss, in a test comment, `.claude/CLAUDE.md` and a
  `ROADMAP.md` known issue; 991355f merges main (#93's squash, a98c071's tree) and changes
  no file. None changes a `src` file, so the refresh's f9361cf figures
  hold on it. Merge it into graph-law with a merge commit, as P4-37(b) orders, so that
  Task 14's CI runs carry the 20 µs gate and its CLAUDE.md edit applies to graph-guard's
  text:
  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && git status --short && [ "$(git rev-parse --short graph-guard)" = 991355f ] && git rev-parse --short graph-guard && git diff --stat f9361cf graph-guard -- src && git merge --no-ff graph-guard -m "Merge graph-guard: gate the constructor benchmark at 20 µs and name the slowest spec" -m $'Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>\nClaude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE' && git log --oneline graph-guard..HEAD && git diff --stat graph-guard HEAD -- src test && cmp Manifest.toml $S/phase45/guard_wt/Manifest.toml && echo manifest-equal
  ```
  Expected: no status lines; `991355f`; no change under `src` from f9361cf to graph-guard;
  `Merge made by the 'ort' strategy.` and the four files graph-guard changed after f9361cf
  (`.claude/CLAUDE.md`, `ROADMAP.md`, `test/test_rate_eq_derivation.jl`,
  `test/test_reaction_graph.jl`); the merge commit and at most the plan's commit; no
  change under `src` or `test` from graph-guard to HEAD; `manifest-equal`. If graph-guard
  has moved past 991355f, the chain stops before the merge: stop and report.
  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; mkdir -p $S/phase45/logs $S/phase45/probes $S/phase45/scripts $S/phase45/p5_base
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T0_inst.log 1800 --project=. -e 'using Pkg; Pkg.instantiate()' && git status --short && git log --oneline -1
  ```
  Expected: the log ends `exit 0`, and the status is empty (`Manifest.toml` is ignored).

- [ ] **Step 2: Run the full suite on the head.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T0_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'
  ```
  Run it in the background, then poll:
  `timeout 580 bash -c 'until grep -q "^exit " LOG; do sleep 20; done'; tail -40 LOG`
  with `LOG=$S/phase45/logs/P5_T0_full.log`.

  Expected:
  - `Testing EnzymeRates tests passed`, then `exit 0`;
  - no warning or error after "Testing Running tests", other than the compile-budget and
    size-gate `[ Info:` lines.

  The fix wave's last full run (`$S/phase45/logs/P4_FW_full.log`, at 342ca14; f9361cf
  changes one test comment) passed 45,895 tests in 7m01.6s, maxrss 2,884,192 KB. Record
  the `EnzymeRates.jl` pass count N0, the wall time and `maxrss_kb` M0. Every later
  full run compares with them.

- [ ] **Step 3: Time ter-ter init and the enumeration file.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T0_terter.log 900 --project=. --heap-size-hint=2500M $S/phase45/scripts/P4_terter.jl
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T0_tme.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_mechanism_enumeration.jl")'
  ```
  Expected:
  - `ter-ter init: 251107 mechanisms in <t> s; the guard rejects 0`. FW-09 measured
    29.5 s fresh; the gate is 150 s.
  - The enumeration log passes and ends `exit 0 elapsed <s>s` (about 150 s; CLAUDE.md
    says about 2.5 min).

  Record both figures; they are two of R52's three figures.

- [ ] **Step 4: Write the set builder.** Create `$S/phase45/scripts/P5_sets.jl`:

  ```julia
  # ABOUTME: The phase-5 baseline's mechanism sets, rebuilt from the worktree's test files: the
  # ABOUTME: specs, population sample, bi-bi seeds, :OnlyI fixtures, gate laws, inhibitor children.
  # Include after the fixtures file, from the worktree root, with the test environment active.

  """The name a top-level expression of a test file defines (`const x = …`, `function f(…)`,
  `f(…) = …`, either under a docstring), or `nothing`."""
  function p5_defined_name(ex)
      ex isa Expr || return nothing
      ex.head === :macrocall && ex.args[1] isa GlobalRef &&
          ex.args[1].name === Symbol("@doc") && return p5_defined_name(ex.args[end])
      ex.head === :const && return p5_defined_name(ex.args[1])
      ex.head === :(=) && ex.args[1] isa Symbol && return ex.args[1]
      ex.head in (:(=), :function) && ex.args[1] isa Expr && ex.args[1].head === :call &&
          return ex.args[1].args[1]
      nothing
  end

  """Evaluate in `Main`, in file order, the top-level definitions of the test file `path`
  whose names are in `wanted`; error when one is missing."""
  function p5_eval_definitions(path::AbstractString, wanted::Vector{Symbol})
      found = Symbol[]
      for ex in Meta.parseall(read(path, String); filename = path).args
          x = p5_defined_name(ex)
          x in wanted || continue
          Core.eval(Main, ex)
          push!(found, x)
      end
      missing_names = setdiff(wanted, found)
      isempty(missing_names) || error("$path defines none of $missing_names")
  end

  p5_eval_definitions("test/test_oracles.jl",
                      [:_testhelper_pingpong_allo_rxn, :_testhelper_content_text,
                       :_testhelper_population_one_round, :_testhelper_population_sample,
                       :_testhelper_inhibitor_rxn])
  p5_eval_definitions("test/test_rate_eq_derivation.jl",
                      [:_testhelper_largest_enumerated_mechanisms,
                       :_testhelper_large_ter_bi_mechanisms,
                       :_testhelper_random_independent_params_concs,
                       :_testhelper_test_rate_equation_performance])

  """The concrete mechanism's id in identify's CSVs (`_mechanism_id` of its text)."""
  p5_id(m) = EnzymeRates._mechanism_id(sprint(show, m))

  """
  The baseline's sets as `(set, label, mechanism)`, each mechanism concrete: `spec` (the 42
  specs, by name), `sample` (test_oracles.jl's 120-member population sample), `bibi_init`
  (the 239 census bi-bi seeds), `onlyi` (the 13 `_testhelper_onlyi_fixtures`), `gate` (the
  six laws of the size gate) and `inhibitor` (the 60 I-binding inhibitor children of
  test_oracles.jl, picked as that file picks them).
  """
  function p5_sets()
      out = Tuple{String, String, Any}[]
      for s in MECHANISM_TEST_SPECS
          push!(out, ("spec", s.name, _testhelper_lift(s.mechanism)))
      end
      for (i, m) in enumerate(_testhelper_population_sample)
          push!(out, ("sample", string(i), m))
      end
      bibi = @enzyme_reaction begin
          substrates: A[C], B[N]
          products: P[C], Q[N]
      end
      for (i, m) in enumerate(EnzymeRates.init_mechanisms(bibi))
          push!(out, ("bibi_init", string(i), m))
      end
      for (label, am) in pairs(_testhelper_onlyi_fixtures())
          push!(out, ("onlyi", string(label), am))
      end
      for (i, m) in enumerate([_testhelper_largest_enumerated_mechanisms;
                               _testhelper_large_ter_bi_mechanisms])
          push!(out, ("gate", string(i), _testhelper_lift(m)))
      end
      children = [c for m in EnzymeRates.init_mechanisms(_testhelper_inhibitor_rxn)
                  for c in _testhelper_plain_moves(m, _testhelper_inhibitor_rxn)]
      binds_i = unique!([c for c in children
                         if any(x -> x.name === :I, EnzymeRates.catalytic(c).table)])
      sorted = binds_i[sortperm(_testhelper_content_text.(binds_i))]
      for (i, c) in enumerate(sorted[round.(Int, range(1, length(sorted); length = 60))])
          push!(out, ("inhibitor", string(i), c))
      end
      out
  end

  """The counts per set, in set order."""
  p5_counts(sets) = [s => count(x -> x[1] == s, sets)
                     for s in ("spec", "sample", "bibi_init", "onlyi", "gate", "inhibitor")]

  """`f()`'s value as text, or `error: <its message>` when it throws; tabs and newlines
  escaped."""
  p5_try(f) = try
      escape_string(f())
  catch err
      escape_string("error: " * first(sprint(showerror, err), 300))
  end
  ```

  The inhibitor pick repeats test_oracles.jl's testset "the guard equals the RE-speed
  oracle on inhibitor children's modifier faces". If FW changed that selection, copy the
  file's version.

- [ ] **Step 5: Write the baseline script.** Create `$S/phase45/scripts/P5_baseline.jl`:

  ```julia
  # ABOUTME: Writes the phase-5 baseline of the worktree it runs from: law values at seeded points,
  # ABOUTME: solve renderings, printed laws, eq_hash classes, V×τ, gate-law sizes and call timings.
  # Usage, from the worktree root: julia --project=. P5_baseline.jl <out_dir>
  using TestEnv; TestEnv.activate()
  using Test, EnzymeRates, LinearAlgebra, Random
  include(joinpath(pwd(), "test", "mechanism_definitions_for_test_enzyme_derivation.jl"))
  include(joinpath(@__DIR__, "P5_sets.jl"))
  const OUT = ARGS[1]
  mkpath(OUT)
  const ER = EnzymeRates
  t0 = time()
  sets = p5_sets()
  println("sets: ", p5_counts(sets), " in ", round(time() - t0; digits = 1), " s")
  flush(stdout)

  """`name=value` pairs of the dictionary `d`, sorted by name, values as `string` writes
  them (256-bit decimals for BigFloat)."""
  pairs_text(d) = join(("$k=$(string(v))" for (k, v) in sort!(collect(d); by = first)), ";")

  open(joinpath(OUT, "laws.tsv"), "w") do io
      println(io, "set\tid\tpoint\tn\tparams\tconcs\tvalue")
      for (set, _, m) in sets
          id = p5_id(m)
          sub = _testhelper_single_subunit(m)
          fitted = collect(ER._dependent_param_exprs(sub)[2])
          rng = Random.MersenneTwister(parse(UInt64, id; base = 16))
          n = m isa ER.AllostericMechanism ? ER.catalytic_multiplicity(m) : 1
          for point in 1:3
              _, concs, params = _testhelper_law_inputs(sub, fitted, rng)
              for (k, mech) in (n > 1 ? ((1, sub), (n, m)) : ((1, sub),))
                  v = _testhelper_law_value(mech, params, concs)
                  println(io, join((set, id, point, k, pairs_text(params), pairs_text(concs),
                                    string(v)), '\t'))
              end
          end
      end
  end
  println("laws.tsv in ", round(time() - t0; digits = 1), " s"); flush(stdout)

  open(joinpath(OUT, "dependents.tsv"), "w") do io
      println(io, "set\tid\tdep\tindep\tparameters\tp_values\tall_keys\tall_values")
      for (set, _, m) in sets
          id = p5_id(m)
          dep, indep = ER._dependent_param_exprs(m)
          _, _, params = _testhelper_law_inputs(_testhelper_single_subunit(m), collect(indep),
                                                Random.MersenneTwister(parse(UInt64, id;
                                                                             base = 16)))
          p = NamedTuple{(indep..., :Keq)}(Tuple(Float64(params[k]) for k in (indep..., :Keq)))
          allp = try
              all_parameters(m, p)
          catch err
              "error: " * first(sprint(showerror, err), 300)
          end
          all_keys, all_values = allp isa NamedTuple ?
              (join(keys(allp), ","), join(repr.(values(allp)), ";")) : (escape_string(allp), "")
          println(io, join((set, id, escape_string(repr(sort!(collect(dep); by = first))),
                            repr(indep), p5_try(() -> repr(ER.parameters(m))),
                            join(repr.(values(p)), ";"), all_keys, all_values), '\t'))
      end
  end
  println("dependents.tsv in ", round(time() - t0; digits = 1), " s"); flush(stdout)

  open(joinpath(OUT, "strings.tsv"), "w") do io
      println(io, "set\tid\tstring")
      for (set, _, m) in sets
          println(io, join((set, p5_id(m), p5_try(() -> rate_equation_string(m))), '\t'))
      end
  end
  println("strings.tsv in ", round(time() - t0; digits = 1), " s"); flush(stdout)

  open(joinpath(OUT, "eqhash.tsv"), "w") do io
      println(io, "id\teq_hash")
      for (set, _, m) in sets
          set == "bibi_init" || continue
          println(io, p5_id(m), '\t',
                  repr(ER._rate_eq_dedup_key(rate_equation_string(m))))
      end
  end

  open(joinpath(OUT, "complexity.tsv"), "w") do io
      println(io, "set\tid\tall\tactive\tinactive")
      for (set, _, m) in sets
          allo = m isa ER.AllostericMechanism
          cm = ER.catalytic(m)
          println(io, join((set, p5_id(m), string(ER._eq_complexity(m)),
                            allo ? string(ER._eq_complexity(cm, ER._active_mask(m))) : "",
                            allo ? string(ER._eq_complexity(cm, ER._inactive_mask(m))) : ""),
                           '\t'))
      end
  end

  open(joinpath(OUT, "dterms.tsv"), "w") do io
      println(io, "label\tid\tdenominator_terms")
      for (set, label, m) in sets
          set == "gate" || continue
          println(io, label, '\t', p5_id(m), '\t', length(ER._raw_symbolic_rate_polys(m)[2]))
      end
  end
  println("eqhash, complexity, dterms in ", round(time() - t0; digits = 1), " s")
  flush(stdout)

  open(joinpath(OUT, "timing.tsv"), "w") do io
      println(io, "kind\tlabel\tallocs\tns_per_call")
      timed = [[("spec", s.name, s.mechanism, s.metabolite_names) for s in MECHANISM_TEST_SPECS];
               [("gate", string(i), em, collect(metabolites(em)))
                for (i, em) in enumerate([_testhelper_largest_enumerated_mechanisms;
                                          _testhelper_large_ter_bi_mechanisms])]]
      for (kind, label, em, mets) in timed
          params, concs, _ = _testhelper_random_independent_params_concs(
              em, mets; rng = Random.MersenneTwister(42))
          allocs, t = _testhelper_test_rate_equation_performance(em, params, concs)
          println(io, join((kind, label, allocs, round(1e9 * t; digits = 2)), '\t'))
      end
  end
  println("timing.tsv; done in ", round(time() - t0; digits = 1), " s")
  ```

  `_testhelper_law_value` reads a concrete mechanism without the guard. So the two
  `:OnlyI` fixtures the guard rejects (`copy_group_bridge`, `essential_activator`) get
  law values, and their `parameters`, `all_parameters` and strings record the guard's
  message.

- [ ] **Step 6: Write the comparison script.** Create `$S/phase45/scripts/P5_compare.jl`:

  ```julia
  # ABOUTME: Compares the worktree it runs from with the phase-5 baseline (P5_baseline.jl): law values,
  # ABOUTME: solve renderings, printed laws, V×τ or gate-law sizes; exit 1 past the kind's tolerance.
  # Usage, from the worktree root: julia --project=. P5_compare.jl <base_dir> <kind>
  # kind: laws (1e-20 relative), dependents (repr bytes; all_parameters 1e-14 relative),
  # strings, complexity, dterms (bytes).
  using TestEnv; TestEnv.activate()
  using Test, EnzymeRates, LinearAlgebra, Random
  include(joinpath(pwd(), "test", "mechanism_definitions_for_test_enzyme_derivation.jl"))
  include(joinpath(@__DIR__, "P5_sets.jl"))
  const BASE, KIND = ARGS[1], ARGS[2]
  const ER = EnzymeRates
  const MECH = Dict((set, p5_id(m)) => m for (set, _, m) in p5_sets())

  """The rows of the baseline file `name` as vectors of fields, header dropped."""
  rows(name) = [split(l, '\t'; keepempty = true)
                for l in Iterators.drop(eachline(joinpath(BASE, name)), 1)]

  """The `name=value` pairs of a laws.tsv field, as 256-bit BigFloats."""
  values_of(field) = Dict{Symbol, BigFloat}(
      Symbol(k) => setprecision(() -> parse(BigFloat, v), BigFloat, 256)
      for (k, v) in (split(kv, '='; limit = 2) for kv in split(field, ';') if !isempty(kv)))

  relerr(a, b) = a == b ? zero(float(a)) : abs(a - b) / max(abs(a), abs(b))

  """The law value of `mech` at `params` and `concs`, and how it was read: through `_law`
  and `_testhelper_evaluate_law` when the tree defines both for `mech` (the law cached under
  `key`), else through `_testhelper_law_value`."""
  function law_value(mech, params, concs, cache, key)
      if isdefined(ER, :_law) && hasmethod(ER._law, Tuple{typeof(mech)}) &&
         isdefined(Main, :_testhelper_evaluate_law)
          L = get!(() -> ER._law(mech), cache, key)
          applicable(_testhelper_evaluate_law, L, Dict{Symbol, BigFloat}()) &&
              return _testhelper_evaluate_law(L, merge(params, concs)), "_law"
      end
      _testhelper_law_value(mech, params, concs), "_num_den_exprs"
  end

  function compare_laws()
      worst = Dict{String, BigFloat}()
      path = Dict{String, Set{String}}()
      cache = Dict{Any, Any}()
      for r in rows("laws.tsv")
          set, id, n = r[1], r[2], parse(Int, r[4])
          m = MECH[(set, id)]
          mech = n == 1 ? _testhelper_single_subunit(m) : m
          v, how = law_value(mech, values_of(r[5]), values_of(r[6]), cache, (set, id, n))
          e = relerr(v, setprecision(() -> parse(BigFloat, r[7]), BigFloat, 256))
          worst[set] = max(get(worst, set, big(0.0)), e)
          push!(get!(path, set, Set{String}()), how)
      end
      for set in sort!(collect(keys(worst)))
          println("laws ", set, ": worst relative error ", Float64(worst[set]), " through ",
                  join(sort!(collect(path[set])), " and "))
      end
      all(<=(big(10)^-20), values(worst))
  end

  function compare_dependents()
      ndiff, worst = 0, 0.0
      for r in rows("dependents.tsv")
          set, id = r[1], r[2]
          m = MECH[(set, id)]
          dep, indep = ER._dependent_param_exprs(m)
          if escape_string(repr(sort!(collect(dep); by = first))) != r[3] ||
             repr(indep) != r[4] || p5_try(() -> repr(ER.parameters(m))) != r[5]
              ndiff += 1
              ndiff <= 5 && println("dependents differ: ", set, " ", id)
              continue
          end
          p = NamedTuple{(indep..., :Keq)}(Tuple(parse.(Float64, split(r[6], ';'))))
          now = try
              all_parameters(m, p)
          catch err
              "error: " * first(sprint(showerror, err), 300)
          end
          if now isa NamedTuple && join(keys(now), ",") == r[7]
              worst = max(worst, maximum(relerr.(collect(values(now)),
                                                 parse.(Float64, split(r[8], ';')))))
          elseif !(now isa String && escape_string(now) == r[7])
              ndiff += 1
              ndiff <= 5 && println("all_parameters differs: ", set, " ", id)
          end
      end
      println("dependents: ", ndiff, " rows differ in bytes; all_parameters worst relative ",
              "error ", worst)
      ndiff == 0 && worst <= 1e-14
  end

  """The fields the current tree gives for baseline row `r` of `kind`."""
  function fields_now(kind, r)
      if kind == "strings"
          return [p5_try(() -> rate_equation_string(MECH[(r[1], r[2])]))]
      elseif kind == "complexity"
          m = MECH[(r[1], r[2])]
          allo = m isa ER.AllostericMechanism
          cm = ER.catalytic(m)
          return [string(ER._eq_complexity(m)),
                  allo ? string(ER._eq_complexity(cm, ER._active_mask(m))) : "",
                  allo ? string(ER._eq_complexity(cm, ER._inactive_mask(m))) : ""]
      end
      m = MECH[("gate", r[2])]
      [string(isdefined(ER, :_law) ? ER._nterms(ER._law(m).D) :
                                     length(ER._raw_symbolic_rate_polys(m)[2]))]
  end

  function compare_bytes(kind)
      file = Dict("strings" => "strings.tsv", "complexity" => "complexity.tsv",
                  "dterms" => "dterms.tsv")[kind]
      ndiff, nrows = 0, 0
      for r in rows(file)
          nrows += 1
          now = fields_now(kind, r)
          now == r[3:end] && continue
          ndiff += 1
          ndiff <= 5 && println(kind, " differ: ", r[1], " ", r[2], ": ", now, " vs ",
                                r[3:end])
      end
      println(kind, ": ", ndiff, " of ", nrows, " rows differ")
      ndiff == 0
  end

  ok = KIND == "laws" ? compare_laws() :
       KIND == "dependents" ? compare_dependents() :
       KIND in ("strings", "complexity", "dterms") ? compare_bytes(KIND) :
       error("kind must be laws, dependents, strings, complexity or dterms; got $KIND")
  exit(ok ? 0 : 1)
  ```

- [ ] **Step 7: Write the identify run.** Create `$S/phase45/scripts/P5_identify_run.jl`:

  ```julia
  # ABOUTME: One seeded identify_rate_equation search on synthetic data from a known ordered bi-bi
  # ABOUTME: mechanism; writes the run's CSVs and a summary of fits run and the selected law.
  # Usage, from the worktree root: julia --project=. P5_identify_run.jl <out_dir>
  using TestEnv; TestEnv.activate()
  using EnzymeRates, Random, CSV, DataFrames, OptimizationCMAEvolutionStrategy
  const OUT = ARGS[1]
  isdir(OUT) && !isempty(readdir(OUT)) && error("$OUT is not empty")
  rxn = @enzyme_reaction begin
      substrates: A[C], B[N]
      products: P[C], Q[N]
  end
  # Ordered bi-bi: A and B bind at rapid equilibrium, the chemistry and both releases are
  # steady-state.
  truth = @enzyme_mechanism begin
      substrates: A, B
      products: P, Q
      steps: begin
          E + A ⇌ E(A)
          E(A) + B ⇌ E(A, B)
          E(A, B) <--> E(P, Q)
          E(P, Q) <--> E(Q) + P
          E(Q) <--> E + Q
      end
  end
  fp = EnzymeRates.fitted_params(truth)
  params = NamedTuple{(fp..., :Keq, :E_total)}((ntuple(i -> 0.5 + 0.5i, length(fp))..., 2.0, 1.0))
  Random.seed!(1)
  draw() = exp(log(0.1) + rand() * log(100.0))
  # Each product a fraction of its substrate drawn from [0.1, 1], independently, so the
  # product terms of the steady-state releases show in the rates.
  frac() = 0.1 + 0.9 * rand()
  rows = map(1:48) do i                     # four groups of 12 points
      a, b = draw(), draw()
      (cld(i, 12), a, b, frac() * a, frac() * b)
  end
  data = (group = ["G$(r[1])" for r in rows], A = [r[2] for r in rows],
          B = [r[3] for r in rows], P = [r[4] for r in rows], Q = [r[5] for r in rows],
          Rate = [rate_equation(truth, (A = r[2], B = r[3], P = r[4], Q = r[5]), params)
                  for r in rows])
  all(>(0), data.Rate) || error("the data are not forward-favorable")
  prob = IdentifyRateEquationProblem(rxn, data; Keq = 2.0)
  t = @elapsed res = identify_rate_equation(prob; optimizer = CMAEvolutionStrategyOpt(),
                                            min_beam_width = 10, n_restarts = 4, maxtime = 5,
                                            save_dir = OUT, show_progress = false)
  searched = [f for f in readdir(OUT) if f == "initial_mechanisms.csv" ||
                                         startswith(f, "equation_search_iteration_")]
  fits = sum(nrow(filter(r -> r.fit_inherited === false,
                         CSV.read(joinpath(OUT, f), DataFrame))) for f in searched)
  best = CSV.read(joinpath(OUT, "best_equation.csv"), DataFrame)
  open(joinpath(OUT, "summary.txt"), "w") do io
      println(io, "seconds: ", round(t; digits = 1))
      println(io, "fits run (rows with fit_inherited == false): ", fits)
      println(io, "selected n_params: ", only(best.n_params))
      println(io, "selected mechanism:\n", only(best.mechanism))
      println(io, "true mechanism:\n", sprint(show, EnzymeRates.Mechanism(truth)))
  end
  print(read(joinpath(OUT, "summary.txt"), String))
  ```

  Each product is 10–100% of its substrate, so the product terms of the truth's
  steady-state releases, which tell it from laws with the same forward family, show in
  the rates; with products at 1% they would not, and the 1-SE rule would prefer a law with
  fewer constants. The mass-action ratio `P·Q/(A·B)` is at most 1, below `Keq = 2`, so
  every rate is positive; the script checks this before it fits.

- [ ] **Step 8: Write the baseline and check the comparison against it.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T0_baseline.log 2400 --project=. --heap-size-hint=2500M $S/phase45/scripts/P5_baseline.jl $S/phase45/p5_base
  ```
  Expected log lines, the first:
  `sets: ["spec" => 42, "sample" => 120, "bibi_init" => 239, "onlyi" => 13, "gate" => 6, "inhibitor" => 60] in <t> s`,
  then one line per file, the last `timing.tsv; done in <t> s`, and `exit 0`. Then:
  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; B=$S/phase45/p5_base; tail -n +2 $B/dterms.tsv | cut -f3 | paste -sd,; tail -n +2 $B/eqhash.tsv | cut -f2 | sort -u | wc -l; tail -n +2 $B/laws.tsv | cut -f1 | sort | uniq -c; awk -F'\t' 'NR > 1 && $3 != 0' $B/timing.tsv | wc -l
  ```
  Expected:
  - `396,396,396,589,799,858`;
  - `239`;
  - laws.tsv row counts `159 spec`, `519 sample`, `717 bibi_init`, `39 onlyi`, `18 gate`,
    `180 inhibitor` (three rows per mechanism, and three more for each allosteric one with
    more than one catalytic subunit: 11 specs and 53 sampled; f9361cf gives these, as
    a9e0699 did);
  - `0` timing rows with allocations.

  A dry run on the a9e0699 snapshot (`$S/phase45/plan5/probes/p1_baseline_dry.sh`) wrote
  the seven files in 95 s. They held 2 `error:` strings, those of `copy_group_bridge` and
  `essential_activator`, and the slowest call was 115.3 ns, on the fifth gate law. The
  refresh's dry run at f9361cf (`$S/phase45/refresh5/map-head-t0-3_run.sh`) wrote them in
  75 s, with the same counts and the same two `error:` strings in the fix wave's wording,
  and its slowest call was 119.0 ns, on the fifth gate law.

  Then check that each comparison passes on the unchanged tree:
  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && for k in laws dependents strings complexity dterms; do $S/phase3/jl.sh $S/phase45/logs/P5_T0_selfcheck_$k.log 1800 --project=. --heap-size-hint=2500M $S/phase45/scripts/P5_compare.jl $S/phase45/p5_base $k; done
  ```
  Expected: each log ends `exit 0`.
  - `laws`: six lines, each `worst relative error 0.0 through _num_den_exprs`;
  - `dependents`: `dependents: 0 rows differ in bytes; all_parameters worst relative
    error 0.0`;
  - `strings: 0 of 480 rows differ`, `complexity: 0 of 480 rows differ`,
    `dterms: 0 of 6 rows differ`.

  `dependents.tsv`, `strings.tsv` and `complexity.tsv` hold one row per mechanism, 480
  each. The a9e0699 dry run printed the `laws`,
  `dependents` and `dterms` lines above, in 30, 55 and 25 s; at f9361cf all five kinds
  printed theirs, in 26, 55, 26, 21 and 21 s. Run the loop in the background and poll.

- [ ] **Step 9: Build the census corpus.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; test -d $S/phase3/census/env_law || $S/phase3/jl.sh $S/phase45/logs/P5_T0_mkenv.log 1800 --heap-size-hint=2500M $S/phase3/census/mkenv.jl $S/phase45/law_wt $S/phase3/census/env_law
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; $S/phase3/jl.sh $S/phase45/logs/P5_T0_census.log 1800 --project=$S/phase3/census/env_law --heap-size-hint=2500M $S/phase3/census/corpus.jl $S/phase45/law_wt $S/phase3/census/P5_base $S/phase3/census/rename_phase3.toml
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; grep -E '^(sig_mismatches|bibi|pingpong|prefix|terbi|twosites|uniuni) = \[' $S/phase3/census/P5_base/meta.toml
  ```
  Expected:
  - the corpus log ends `exit 0` after about 2 min (P4_T12's census took about 110 s);
  - `P5_base` holds `mechs.tsv`, `reactions.tsv`, `init.tsv`, `parents.tsv`,
    `children.tsv`, `seeds.tsv`, `laws.tsv` and `meta.toml`;
  - the grep prints seven lines, each ending `= []`.

- [ ] **Step 10: Run the identify search.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T0_identify.log 5400 --project=. --heap-size-hint=2500M $S/phase45/scripts/P5_identify_run.jl $S/phase45/p5_T0_identify
  ```
  Run it in the background and poll the log in bounded loops. Expected: the log ends
  `exit 0` after the summary:
  - `seconds: <t>`;
  - `fits run (rows with fit_inherited == false): <F0>`;
  - `selected n_params: <n0>`;
  - the selected and the true mechanism's texts.

  The selection is recorded, not required to equal the truth. Task 9 compares its own run
  with these figures: when this run selects the truth, Task 9's must too; when it does
  not, Task 9's must select the truth, the law of this run's selection or a near-exact
  law (`cv_score` below 1e-6) with at most its `n_params`, and the summary names this
  run's selection for that check (`best_equation.csv`'s `mechanism`). A settle run on
  f9361cf (`$S/phase45/refresh5/settle_1.log`) gave `seconds: 102.4`, 599 fits and
  `selected n_params: 7`, a one-subunit L phantom, not the truth. If the
  run passes 90 min, stop and report: the controller decides whether to lower
  `n_restarts` in this script and in Task 9's run alike.

- [ ] **Step 11: Write the figures.** Create `$S/phase45/p5_base/figures.md` from the logs:

  ```markdown
  # Phase-5 baseline (graph-guard <head sha>, Julia <version>)

  - Full suite: <N0> tests, <wall> s, maxrss <M0> KB (logs/P5_T0_full.log)
  - ter-ter init_mechanisms: <t> s, 251107 mechanisms (logs/P5_T0_terter.log)
  - test_mechanism_enumeration.jl focused: <s> s (logs/P5_T0_tme.log)
  - Compile budgets: the `[ Info:` lines of P5_T0_full.log for init trace-compile,
    first rate_equation, Sig encoding, fit compile and show instances, copied verbatim
  - Gate laws, local: the six `[ Info: rate_equation: … ns per term` lines, copied verbatim
  - timing.tsv: slowest spec <name> <ns> ns; gate laws <ns list> ns
  - Identify run: <seconds> s, <F0> fits, selected n_params <n0>; selected = truth: <yes|no>
  - Census: census/P5_base, 0 conflations, 0 sig mismatches
  ```

  Fill every slot from its log; leave no `<…>`.

- [ ] **Step 12: No commit.** Task 0 changes no repository file. `git status --short` in
  `law_wt` prints nothing.

---

### Task 1: Structured solve output (P5-4)

`_solve` returns the reduction phase 2 already runs as one struct, `_Solve`, and
`_dependent_param_exprs` renders it. `all_parameters` reads it directly, so the `Expr`
interpreter `_evaluate_dependent` goes. No output changes.

**Files:**
- Modify: `src/thermo.jl` (`_Dependent`, `_Solve`, `_solve`, `_dependent_param_exprs`
  on `_Solve` and on a mechanism, `_independent_param_count`)
- Modify: `src/rate_eq_derivation.jl` (`all_parameters` reads `_solve`; delete
  `_evaluate_dependent`)
- Test: `test/test_joint_solve.jl` (four testsets appended), `test/test_accessors.jl`
  (one testset converted)
- Create (scratch): `$S/phase45/probes/P5_T1_dependents.jl`
- Logs: `P5_T1_{dependents,solve,accessors,compare,full}.log`

**Interfaces:**
- Consumes: `_ordered_columns`, `_joint_system`, `_reduce_joint_system`, `_kept_members`,
  `_column_parameter`, `_class` (on `_excluded_by`), `_in_active_law`,
  `_has_inactive_columns`, `_kreg_params`, `name(p, m, kept)`, `build_power_expr`,
  `_derivation_input`, `_well_posed`.
- Produces:
  ```julia
  struct _Dependent
      col::Int                                  # the pivot column
      keq::Rational{Int}                        # exponent of Keq
      factors::Vector{Pair{Int, Rational{Int}}} # free column => exponent, by column
      label::Symbol                             # :kinetic_group, :haldane or :wegscheider
  end
  struct _Solve
      columns::Vector{_Column}       # every bare constant's column, in elimination order
      names::Vector{Symbol}          # column i's name (name(_column_parameter(m, c), m, kept))
      column::Dict{_Column, Int}     # each bare constant => its column; an :EqualAI edge's
                                     #   inactive constant => its active one
      independents::Vector{Int}      # the law's free columns, sorted by name
      dependents::Vector{_Dependent} # every pivot column of the law, by column
      fitted::Vector{Symbol}         # fitted_params order: independents' names, then the
                                     #   regulator constants (active, then inactive), then :L
  end
  _solve(m::Union{Mechanism, AllostericMechanism})::_Solve
  _solve(@nospecialize(M::Type{<:AbstractEnzymeMechanism}))::_Solve   # _derivation_input
  _dependent_param_exprs(s::_Solve)::Tuple{Dict{Symbol, Union{Symbol, Expr, Float64}}, Tuple}
  _dependent_param_exprs(m::Union{Mechanism, AllostericMechanism})    # of _solve(m)
  _dependent_param_exprs(@nospecialize(M::Type{<:AbstractEnzymeMechanism}))
  _independent_param_count(m::Union{Mechanism, AllostericMechanism})::Int  # length(fitted)
  all_parameters(m, params::NamedTuple)::NamedTuple   # signature and docstring unchanged
  ```
  A pivot of the law that is not its class's kept member is `:kinetic_group`. Another
  pivot is `:haldane` when `keq != 0`, else `:wegscheider` (RES §1.2).
  `_dependent_param_exprs` renders every dependent except the `:kinetic_group` ones. The
  reduction's `Rational{BigInt}` retry stays in `_reduce_joint_system`; `_solve` narrows
  each value to `Rational{Int}`, and a value outside `Int` throws `InexactError`.

**Gate:**
- `P5_compare.jl $S/phase45/p5_base dependents`: `_dependent_param_exprs` and
  `parameters` are byte-identical on all 480 baseline mechanisms, and `all_parameters`
  agrees to 1e-14 relative (RES: worst 6.2e-16; the a9e0699 probe: within 1e-14 on the
  42 specs; the refresh probe at f9361cf: worst 4.3e-16 on all 480).
- The half-integer fixture's dependent has `keq == -1//2`.
- The full suite is green; +0 s.

**Testsets added (`test/test_joint_solve.jl`):**
- "_solve renders a half-integer dependent with keq -1//2" (the mechanism
  of TR "rate_equation_string runs a half-integer dependent as written", inline);
- "_solve labels kinetic-group, Haldane and Wegscheider dependents" (TR's `folded`
  fixture: a tied group, a binding square and a Haldane cycle);
- "_solve maps an :EqualAI edge's inactive constant to its active column";
- "_solve lists the fitted names in parameter order" (the 42 specs).

**Pins changed:** test_accessors "_evaluate_dependent reads the solve's power products
only" becomes "all_parameters computes rational powers of Keq from the solve". Value
check: the head's 0.25 for `Keq ^ -0.5 * a / a ^ 2` at `a = 2`, `Keq = 4` is reproduced
through `all_parameters` on a fixture with a `Keq^-1/2` dependent. Step 1 reads the value
on the head. The `sqrt(Keq)` assertion's subject goes, and the PR body lists it for Denis.

- [ ] **Step 1: Check the pinned renderings on the head.** The new testsets pin how the
  head renders two fixtures, so first read them off the head. Create
  `$S/phase45/probes/P5_T1_dependents.jl`:

  ```julia
  # ABOUTME: Prints the head's dependent renderings of Task 1's two solve fixtures and the
  # ABOUTME: head's all_parameters on the half-integer one: the values Task 1's pins must keep.
  using EnzymeRates
  folded = @enzyme_mechanism begin
      substrates: A, B
      products: P
      steps: begin
          (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
          E + B ⇌ E(B)
          E(A) + B ⇌ E(A, B)
          E(A, B) <--> E(P)
          E(P) ⇌ E + P
      end
  end
  half = @enzyme_mechanism begin
      substrates: S
      products: P
      steps: begin
          (E(S) <--> F(S), E <--> F)
          F(S) <--> E + P
          F + S <--> E(S)
      end
  end
  for (label, m) in (("folded", folded), ("half", half))
      dep, indep = EnzymeRates._dependent_param_exprs(EnzymeRates.Mechanism(m))
      println(label, " dep ", repr(sort!(collect(dep); by = first)))
      println(label, " indep ", indep)
  end
  p = (k_E_to_F = 1.0, k_E_P_to_FS = 1.0, k_FS_to_E_P = 1.0, k_ES_to_F_S = 4.0,
       k_F_S_to_ES = 1.0, Keq = 4.0)
  println("half k_F_to_E ", all_parameters(half, p).k_F_to_E, " ",
          all_parameters(half, merge(p, (k_F_S_to_ES = 16.0,))).k_F_to_E)
  ```
  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T1_dependents.log 600 --project=. $S/phase45/probes/P5_T1_dependents.jl && head -5 $S/phase45/logs/P5_T1_dependents.log
  ```
  Expected (as on a9e0699, `probes/p1_solve_facts.log`, and at f9361cf,
  `$S/phase45/refresh5/map-head-t0-3_7.log`):
  ```
  folded dep Pair{Symbol, Union{Float64, Expr, Symbol}}[:K_EAB_to_EA_B => :K_EB_to_E_B, :k_EP_to_EAB => :((1 / Keq) * (1 / K_EA_to_E_A) * (1 / K_EB_to_E_B) * K_EP_to_E_P * k_EAB_to_EP)]
  folded indep (:K_EA_to_E_A, :K_EB_to_E_B, :K_EP_to_E_P, :k_EAB_to_EP)
  half dep Pair{Symbol, Union{Float64, Expr, Symbol}}[:k_F_to_E => :(Keq ^ -0.5 * k_ES_to_F_S ^ -0.5 * k_E_P_to_FS ^ -0.5 * k_E_to_F * k_FS_to_E_P ^ 0.5 * k_F_S_to_ES ^ 0.5)]
  half indep (:k_ES_to_F_S, :k_E_P_to_FS, :k_E_to_F, :k_FS_to_E_P, :k_F_S_to_ES)
  half k_F_to_E 0.25 1.0
  ```
  If the head prints otherwise, the head's text is the pin: write Step 2's literals from
  it.

- [ ] **Step 2: Write the failing tests.** Append to the end of `test/test_joint_solve.jl`,
  after the last `end`:

  ```julia

  @testset "_solve renders a half-integer dependent with keq -1//2" begin
      # The conformational square ties two isomerizations of one group, so the Haldane
      # relation fixes k_F_to_E through square roots: Keq = (k_E_to_F / k_F_to_E)² ·
      # (k_F_S_to_ES / k_ES_to_F_S) · (k_FS_to_E_P / k_E_P_to_FS).
      half = EnzymeRates.Mechanism(@enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              (E(S) <--> F(S), E <--> F)
              F(S) <--> E + P
              F + S <--> E(S)
          end
      end)
      s = EnzymeRates._solve(half)
      @test issorted([d.col for d in s.dependents])
      @test sort([d.label for d in s.dependents]) ==
            [:haldane, :kinetic_group, :kinetic_group]
      d = only(d for d in s.dependents if d.label === :haldane)
      @test s.names[d.col] === :k_F_to_E && d.keq == -1//2
      @test Dict(s.names[f] => e for (f, e) in d.factors) ==
            Dict(:k_E_to_F => 1//1, :k_E_P_to_FS => -1//2, :k_FS_to_E_P => 1//2,
                 :k_ES_to_F_S => -1//2, :k_F_S_to_ES => 1//2)
      @test EnzymeRates._dependent_param_exprs(s)[1] ==
            Dict(:k_F_to_E => :((Keq ^ -0.5) * (k_ES_to_F_S ^ -0.5) * (k_E_P_to_FS ^ -0.5) *
                                k_E_to_F * (k_FS_to_E_P ^ 0.5) * (k_F_S_to_ES ^ 0.5)))
  end

  @testset "_solve labels kinetic-group, Haldane and Wegscheider dependents" begin
      # A binds E and E(B) under one group (a kinetic group: its occupied member's
      # constant equals the kept one); B's two bindings close a binding square with A's,
      # which ties them (Wegscheider); the cycle through the catalytic step fixes its
      # reverse rate (Haldane).
      folded = EnzymeRates.Mechanism(@enzyme_mechanism begin
          substrates: A, B
          products: P
          steps: begin
              (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
              E + B ⇌ E(B)
              E(A) + B ⇌ E(A, B)
              E(A, B) <--> E(P)
              E(P) ⇌ E + P
          end
      end)
      s = EnzymeRates._solve(folded)
      described = [(s.names[d.col], d.label, d.keq, [s.names[f] => e for (f, e) in d.factors])
                   for d in s.dependents]
      @test described == [
          (:k_EP_to_EAB, :haldane, -1//1, [:k_EAB_to_EP => 1//1, :K_EP_to_E_P => 1//1,
                                           :K_EB_to_E_B => -1//1, :K_EA_to_E_A => -1//1]),
          (:K_EA_to_E_A, :kinetic_group, 0//1, [:K_EA_to_E_A => 1//1]),
          (:K_EAB_to_EA_B, :wegscheider, 0//1, [:K_EB_to_E_B => 1//1])]
      @test s.fitted == [:K_EA_to_E_A, :K_EB_to_E_B, :K_EP_to_E_P, :k_EAB_to_EP]
      @test s.names[s.independents] == s.fitted
      @test EnzymeRates._dependent_param_exprs(s) ==
            (Dict(:K_EAB_to_EA_B => :K_EB_to_E_B,
                  :k_EP_to_EAB => :((1 / Keq) * (1 / K_EA_to_E_A) * (1 / K_EB_to_E_B) *
                                    K_EP_to_E_P * k_EAB_to_EP)),
             (:K_EA_to_E_A, :K_EB_to_E_B, :K_EP_to_E_P, :k_EAB_to_EP))
  end

  @testset "_solve maps an :EqualAI edge's inactive constant to its active column" begin
      am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)        :: EqualAI
              E(S) <--> E + P     :: NonequalAI
          end
      end)
      s = EnzymeRates._solve(am)
      C = EnzymeRates._Column
      binding = _testhelper_edge(am, "E + S ⇌ ES")
      release = _testhelper_edge(am, "E + P <--> ES")
      @test length(s.columns) == 5
      @test s.column[C(binding, 1, true)] == s.column[C(binding, 1, false)]
      @test s.columns[s.column[C(binding, 1, false)]] == C(binding, 1, false)
      for slot in 1:2, inactive in (false, true)
          @test s.columns[s.column[C(release, slot, inactive)]] == C(release, slot, inactive)
      end
      @test s.fitted == [:K_ES_to_E_S, :k_A_ES_to_E_P, :k_I_ES_to_E_P, :L]
  end

  @testset "_solve lists the fitted names in parameter order" begin
      for spec in MECHANISM_TEST_SPECS
          @test Tuple(EnzymeRates._solve(_testhelper_lift(spec.mechanism)).fitted) ==
                EnzymeRates.fitted_params(spec.mechanism)
      end
  end
  ```

  The constructor stores `E(S) <--> E + P` as `E + P <--> E(S)`: a step with one
  metabolite takes it up. That is why the release edge's text reads `E + P <--> ES`.

  In `test/test_accessors.jl`, replace the last testset:
  ```julia
  @testset "_evaluate_dependent reads the solve's power products only" begin
      p = (a = 2.0, Keq = 4.0)
      @test EnzymeRates._evaluate_dependent(:(Keq ^ -0.5 * a / a ^ 2), p) == 0.25
      @test_throws ErrorException EnzymeRates._evaluate_dependent(:(sqrt(Keq)), p)
  end
  ```
  with:
  ```julia
  @testset "all_parameters computes rational powers of Keq from the solve" begin
      # The dependent of this conformational square is k_F_to_E = Keq^-1/2 ·
      # k_ES_to_F_S^-1/2 · k_E_P_to_FS^-1/2 · k_E_to_F · k_FS_to_E_P^1/2 · k_F_S_to_ES^1/2
      # (test_joint_solve.jl): 1/4 at Keq = k_ES_to_F_S = 4 and every other constant 1,
      # and 1 once k_F_S_to_ES is 16. BigFloat inputs give the exact value.
      half = @enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              (E(S) <--> F(S), E <--> F)
              F(S) <--> E + P
              F + S <--> E(S)
          end
      end
      p = (k_E_to_F = 1.0, k_E_P_to_FS = 1.0, k_FS_to_E_P = 1.0, k_ES_to_F_S = 4.0,
           k_F_S_to_ES = 1.0, Keq = 4.0)
      @test all_parameters(half, p).k_F_to_E == 0.25
      @test all_parameters(half, merge(p, (k_F_S_to_ES = 16.0,))).k_F_to_E == 1.0
      @test all_parameters(half, map(big, p)).k_F_to_E == big(0.25)
  end
  ```
  The `sqrt(Keq)` assertion's subject, `_evaluate_dependent`, goes in this task. The PR
  body lists it for Denis (Task 15).

- [ ] **Step 3: Watch the solve tests fail.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T1_solve.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_joint_solve.jl")'
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T1_accessors.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_accessors.jl")'
  ```
  Expected:
  - `P5_T1_solve.log`: "Joint thermodynamic solve" passes. The first new testset errors
    with ``UndefVarError: `_solve` not defined in `EnzymeRates` ``. A top-level testset
    that errors throws when it ends, so the include stops there with `Some tests did not
    pass`, the other three new testsets do not run, and the log ends `exit 1`.
  - `P5_T1_accessors.log` passes. The converted testset is the value check: it holds on
    the head and must hold after the change.

- [ ] **Step 4: Write `_solve` and its renderer.** In `src/thermo.jl`, replace everything
  from the line
  `# ─── Dependents ────────────────────────────────────────────────────────────`
  through the line
  `_independent_param_count(m::_AnyMech) = length(_dependent_param_exprs(m)[2])`
  (the long `_dependent_param_exprs(m)` docstring and function, the `Type` method and
  `_independent_param_count`) with:

  ```julia
  # ─── The solve ─────────────────────────────────────────────────────────────

  """
  A dependent constant of the joint solve: pivot column `col` equals `Keq^keq · Π x_f^e`
  over the free columns `f => e` of `factors`, in column order. `label` says why it is
  dependent: `:kinetic_group` for a column that is not its class's kept member, which it
  equals; otherwise `:haldane` when the relation reads `Keq` and `:wegscheider` when it does
  not.
  """
  struct _Dependent
      col::Int
      keq::Rational{Int}
      factors::Vector{Pair{Int, Rational{Int}}}
      label::Symbol
  end

  """
  The joint solve of a mechanism (`_solve`). `columns` holds a column (`_Column`) per bare
  constant, in elimination order, and `names` each column's name, that of its class's kept
  member (`_kept_members`). `column` maps each bare constant to its column; an `:EqualAI`
  edge's inactive constant reads its active column. `independents` lists the law's free
  columns, sorted by name (`isless` on the strings, so `K_…` precedes `k_…`), and
  `dependents` every pivot column of the law, in column order. `fitted` holds the fitted
  constants in `fitted_params` order: the independents' names, then an allosteric
  mechanism's regulator constants (`_kreg_params`: active, then inactive, site by site)
  and `L`.
  """
  struct _Solve
      columns::Vector{_Column}
      names::Vector{Symbol}
      column::Dict{_Column, Int}
      independents::Vector{Int}
      dependents::Vector{_Dependent}
      fitted::Vector{Symbol}
  end

  """
      _solve(m) -> _Solve

  The joint solve of `m` (`_joint_system`, `_reduce_joint_system`), its exact rationals
  narrowed to `Rational{Int}`; a value outside `Int` throws. Every class of constants of the
  law is independent or dependent under the name of its kept member (`_kept_members`).
  Columns outside the law (an `:OnlyA` edge's inactive copy, an `:OnlyI` edge's active copy,
  a copy on a form its conformation prunes) are eliminated first, so no dependent of the
  law reads them. The `Type` method decodes a compiled mechanism and checks that its rate
  law is defined (`_derivation_input`).

  Errors when two constants of the law, regulator constants included, render one name,
  which the law would read as one constant: with a conformation labelled `A`, the uptake
  constant of the `:EqualAI` step `A(Eb) <--> A + Eb` renders `k_A_Eb_to_AEb`, the name of
  the active constant of a `:NonequalAI` step `E(b) <--> A(Eb)`. The macros and the
  constructors reject the labels `A` and `I` (`_check_conformation_label`). A metabolite
  name can still imitate a residual's: with a modifier named `minusQ`, the side
  `E(; residual = A - P) + minusQ` and the form `E(; residual = A - P - Q)` both render
  `E_res_plusA_minusP_minusQ`, and the macros accept a step between them, so this check can
  fire for a mechanism they build.
  """
  function _solve(m::_AnyMech)
      cols = _ordered_columns(m)
      sys = _joint_system(m, cols)
      M, pivots = _reduce_joint_system(sys.A, sys.rhs)
      n = length(sys.columns)
      kept = _kept_members(m, cols)
      names = [name(_column_parameter(m, c), m, kept) for c in sys.columns]
      group = catalytic(m).group
      class(c) = _class(c, group[c.edge], sys.tags[group[c.edge]])
      last_of = Dict(class(c) => i for (i, c) in enumerate(sys.columns))
      in_law(c) = c.inactive ? sys.law[c.edge] :
                  _in_active_law(c.edge, group, sys.tags, sys.law, cols.alaw)
      kept_in_law(i) = in_law(sys.columns[i]) && last_of[class(sys.columns[i])] == i
      reg = m isa Mechanism ? Symbol[] :
            Symbol[name(p, m) for p in [_kreg_params(m, :A); _kreg_params(m, :I)]]
      rendered = [names[filter(kept_in_law, 1:n)]; reg]
      allunique(rendered) || error(
          "Constants of this mechanism's rate equation share a name (" *
          join(unique(x for x in rendered if count(==(x), rendered) > 1), ", ") *
          "); rename a metabolite or a conformation label so that each constant has a " *
          "name of its own")
      free = setdiff(1:n, pivots)
      dependents = _Dependent[]
      for (r, p) in enumerate(pivots)
          in_law(sys.columns[p]) || continue
          keq = Rational{Int}(M[r, n + 1])
          factors = Pair{Int, Rational{Int}}[f => Rational{Int}(-M[r, f])
                                             for f in free if !iszero(M[r, f])]
          label = !kept_in_law(p) ? :kinetic_group : iszero(keq) ? :wegscheider : :haldane
          push!(dependents, _Dependent(p, keq, factors, label))
      end
      independents = sort!(filter(kept_in_law, free); by = i -> string(names[i]))
      column = Dict(c => i for (i, c) in enumerate(sys.columns))
      if m isa AllostericMechanism
          for (i, c) in enumerate(sys.columns)
              c.inactive || _has_inactive_columns(sys.tags[group[c.edge]]) ||
                  (column[_Column(c.edge, c.slot, true)] = i)
          end
      end
      fitted = names[independents]
      m isa AllostericMechanism && append!(fitted, reg, (:L,))
      _Solve(sys.columns, names, column, independents, dependents, fitted)
  end

  _solve(@nospecialize(M::Type{<:AbstractEnzymeMechanism})) = _solve(_derivation_input(M))

  """
      _dependent_param_exprs(s::_Solve) -> (dep, indep)

  The solve `s` as expressions: `dep` maps each dependent constant that is its class's kept
  member to `Keq^q · Π independent^e` with rational exponents (`build_power_expr`; the empty
  product is `1.0`); `indep` is `Tuple(s.fitted)`. A kinetic group's other members equal the
  kept one and add no entry. The methods on a mechanism render `_solve(m)`.
  """
  function _dependent_param_exprs(s::_Solve)
      dep = Dict{Symbol, Union{Symbol, Expr, Float64}}()
      for d in s.dependents
          d.label === :kinetic_group && continue
          dep[s.names[d.col]] =
              build_power_expr(d.keq, [(s.names[f], e) for (f, e) in d.factors])
      end
      dep, Tuple(s.fitted)
  end
  _dependent_param_exprs(m::_AnyMech) = _dependent_param_exprs(_solve(m))
  _dependent_param_exprs(@nospecialize(M::Type{<:AbstractEnzymeMechanism})) =
      _dependent_param_exprs(_solve(M))

  """Number of independent (fitted) constants of a concrete mechanism, from the joint solve,
  without compiling the mechanism."""
  _independent_param_count(m::_AnyMech) = length(_solve(m).fitted)
  ```

  The closures `class`, `in_law` and `kept_in_law` capture locals that are assigned once,
  so nothing is boxed. `M` is `Matrix{Rational{Int}}` or, after the overflow retry,
  `Matrix{Rational{BigInt}}`; `Rational{Int}(x)` is the identity on the first and narrows
  the second.

- [ ] **Step 5: Read the solve in `all_parameters`.** In `src/rate_eq_derivation.jl`,
  replace:
  ```julia
  function all_parameters(m::Union{AbstractEnzymeMechanism, Mechanism, AllostericMechanism},
                          params::NamedTuple)
      dep, indep = _dependent_param_exprs(m isa AbstractEnzymeMechanism ? typeof(m) :
                                          _well_posed(m))
      absent = [k for k in (indep..., :Keq) if !haskey(params, k)]
      isempty(absent) || error("all_parameters: params must hold every name of " *
                               "parameters(m) except E_total; it lacks " *
                               join(absent, ", "))
      names = sort!(collect(keys(dep)))
      merge(params, NamedTuple{Tuple(names)}(
          Tuple(_evaluate_dependent(dep[k], params) for k in names)))
  end

  """The value at `params` of a dependent constant's expression (`build_power_expr`): a
  power product of `Keq` and independent constants, or the number it is pinned to."""
  _evaluate_dependent(x::Symbol, params) = params[x]
  _evaluate_dependent(x::Real, params) = x
  function _evaluate_dependent(e::Expr, params)
      op, args = e.args[1], [_evaluate_dependent(a, params) for a in e.args[2:end]]
      op === :* && return prod(args)
      op === :/ && return args[1] / args[2]
      op === :^ && return args[1]^args[2]
      error("_evaluate_dependent: no rule for the operator $op")
  end
  ```
  with:
  ```julia
  function all_parameters(m::Union{AbstractEnzymeMechanism, Mechanism, AllostericMechanism},
                          params::NamedTuple)
      s = _solve(m isa AbstractEnzymeMechanism ? typeof(m) : _well_posed(m))
      absent = [k for k in (s.fitted..., :Keq) if !haskey(params, k)]
      isempty(absent) || error("all_parameters: params must hold every name of " *
                               "parameters(m) except E_total; it lacks " *
                               join(absent, ", "))
      # A negative integer exponent inverts the positive power; a fractional one is taken
      # in the value's own precision, so a BigFloat value stays exact.
      pow(x, e) = !isone(denominator(e)) ?
                      x^(oftype(float(x), numerator(e)) / denominator(e)) :
                  numerator(e) >= 0 ? x^numerator(e) : inv(x^-numerator(e))
      value(d) = prod((pow(params[s.names[f]], e) for (f, e) in d.factors);
                      init = pow(params.Keq, d.keq))
      deps = sort!(filter(d -> d.label !== :kinetic_group, s.dependents);
                   by = d -> s.names[d.col])
      merge(params, NamedTuple{Tuple(s.names[d.col] for d in deps)}(
          Tuple(value(d) for d in deps)))
  end
  ```
  The docstring above `all_parameters` stays as it is. For a `Float64` value the powers
  are those the old interpreter computed (`inv(y)` is `1.0 / y`), and an integer value
  still takes a negative power.

- [ ] **Step 6: Watch the tests pass.** Rerun both commands of Step 3. Expected: both logs
  pass, with the four new testsets and the converted one green.

- [ ] **Step 7: Compare the solve with the baseline.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T1_compare.log 1800 --project=. --heap-size-hint=2500M $S/phase45/scripts/P5_compare.jl $S/phase45/p5_base dependents
  ```
  Expected:
  `dependents: 0 rows differ in bytes; all_parameters worst relative error <e>`, with `<e>`
  at most 1e-14, then `exit 0`. A byte difference is a bug in Step 4, not a re-pin.

- [ ] **Step 8: Run the full suite** (the Step 2 command of Task 0, log
  `$S/phase45/logs/P5_T1_full.log`, in the background, polled). Expected: green and
  pristine. The pass count is N0 plus this task's new tests (60: 5, 4, 8 and 42
  across the four testsets, and one more in the converted one). Wall time and maxrss are
  within noise of Task 0's.

- [ ] **Step 9: Commit.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && git add src/thermo.jl src/rate_eq_derivation.jl test/test_joint_solve.jl test/test_accessors.jl && git commit -F - <<'EOF'
  Return the joint solve as one struct, _Solve

  _solve returns the reduction as columns, names, the law's independents and
  dependents, each dependent with its Keq exponent and its factors as
  Rational{Int}, and the fitted names. _dependent_param_exprs renders it,
  byte-identically on the 480 mechanisms of the phase-5 baseline, and
  all_parameters reads it directly, so the Expr interpreter
  _evaluate_dependent goes (its test of an unknown operator is listed for
  Denis). The engine, code generation and the printer read the same struct
  in later tasks.

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
  EOF
  ```

---

### Task 2: Exact polynomials, `poly.jl` (P5-1)

**Files:**
- Create: `src/poly.jl`; modify `src/EnzymeRates.jl` (include it after `thermo.jl`)
- Create: `test/test_poly.jl`; modify `test/runtests.jl`
- Modify: the spec, section 9 "Engine", third bullet (P5-1 ◆, listed for Denis)
- Create (scratch): `$S/phase45/gcscan5/mkenv.jl`, `$S/phase45/gcscan5/scan_law.jl`
  (workload `poly`); environments `$S/phase45/gcscan5/env112`, `env113`
- Logs: `P5_T2_{poly,env112,env113,scan112_poly,scan112_control,scan112_repro,scan113_poly,scan113_control,scan113_repro,hoist,full}.log`

**Interfaces:**
- Consumes: none from earlier tasks.
- Produces:
  ```julia
  struct _Vars; names::Vector{Symbol}; conc::BitVector; index::Dict{Symbol, Int}; end
  _Vars()::_Vars
  Base.:(==)(a::_Vars, b::_Vars); Base.hash(t::_Vars, h::UInt)   # names and conc
  _var!(t::_Vars, s::Symbol; conc::Bool = false)::Int
  const _Exp = Rational{Int}
  const _Mono = Vector{Pair{Int, _Exp}}          # sorted by variable, no zero exponent
  _mmul(a::_Mono, b::_Mono)::_Mono
  _mpow(a::_Mono, e::Union{Integer, _Exp})::_Mono
  _monomial(ps::Pair{Int}...)::_Mono             # a variable may repeat
  _mexp(a::_Mono, v::Int)::_Exp
  mutable struct _Poly; const terms::Dict{_Mono, Int}; end
  _Poly()::_Poly
  _term(m::_Mono, c::Integer = 1)::_Poly
  Base.one(::Type{_Poly}); Base.zero(::Type{_Poly}); Base.iszero(p::_Poly)
  Base.:(==)(p::_Poly, q::_Poly); Base.hash(p::_Poly, h::UInt)   # over terms
  _nterms(p::_Poly)::Int
  _addto!(p::_Poly, m::_Mono, c::Integer)::_Poly                 # checked add
  Base.:+(p::_Poly, q::_Poly); Base.:-(p::_Poly); Base.:-(p::_Poly, q::_Poly)
  Base.:*(p::_Poly, q::_Poly); Base.:*(p::_Poly, m::_Mono)       # checked mul
  _scale(p::_Poly, c::Integer)::_Poly
  _subst(f, p::_Poly)::_Poly        # f(v::Int)::Union{Nothing, _Mono}; nothing keeps v
  _divexact(p::_Poly, d::_Poly)::_Poly   # two-term d; ArgumentError when d does not divide
  _lowest_shift(ps, counted = _ -> true)::_Mono
  _conc_shift(ps, t::_Vars)::_Mono
  _split_conc(p::_Poly, t::_Vars)::Dict{_Mono, _Poly}   # concentration monomial => coefficient
  _content(p::_Poly)::Tuple{Int, _Mono, _Poly}          # integer, monomial, primitive part
  ```
  These port `SLC/poly.jl`, plus `content` from `SLC/print.jl`, without the type
  parameter. The arithmetic goes through `Base.checked_add`, `Base.checked_sub` and
  `Base.checked_mul`. `_content`'s third part is `one(_Poly)` where the slice returned
  `nothing`, so the result has one concrete type (R18).
- Scratch: `julia --project=$S/phase45/gcscan5/env<v> --pkgimages=no scan_law.jl poly [control]`
  runs six loops three times; the loops merge freshly built polynomials, held in toy
  law-like structs, with `Nothing`, `Vector` or a holder. It prints
  `round <r>: same (…), full gc ok`. With `control`, the polynomial is the script's
  immutable one-field `ImmPoly`. Task 6 adds the workloads `derive` and `identify`.

**Gate:**
- `test_poly.jl` is green, and a forced overflow throws `OverflowError`.
- The GC scan (P5-7) on 1.12.5 and 1.13.1:
  - `scan_law.jl poly` reports 0 hoisted allocations with a bits store into a reference
    field;
  - `scan_law.jl poly control` reports at least 1 and aborts the collector;
  - `repro_upstream.jl` reports at least 1 and aborts;
  - each scan takes a few minutes.
- The full suite is green.

**Testsets added (`test/test_poly.jl`):** "monomials multiply, raise and read exponents
exactly", "polynomials add, subtract and multiply exactly", "an Int64 overflow throws
OverflowError", "_divexact divides out the Haldane factor", "_subst replaces variables by
power products", "the shifts bring polynomials to lowest terms in the concentrations",
"_content splits off the integer and monomial content". **Pins changed:** none.

- [ ] **Step 1: Write the failing test.** Create `test/test_poly.jl`:

  ```julia
  # ABOUTME: Tests of the exact sparse polynomials (src/poly.jl): monomials, checked
  # ABOUTME: arithmetic, division by the Haldane factor, substitution, the shifts and content.

  @testset "monomials multiply, raise and read exponents exactly" begin
      P = EnzymeRates
      @test P._monomial(2 => 1, 1 => 1//2, 2 => 2) == [1 => 1//2, 2 => 3//1]
      @test P._monomial(1 => 1, 1 => -1) == P._Mono()
      @test P._monomial(1 => 0, 3 => 2) == [3 => 2//1]
      @test P._mmul(P._monomial(1 => 1//2), P._monomial(1 => -1//2, 3 => 2)) == [3 => 2//1]
      @test P._mpow(P._monomial(1 => 1//2, 2 => -1), 3) == [1 => 3//2, 2 => -3//1]
      @test P._mpow(P._monomial(1 => 2), 1//2) == [1 => 1//1]
      @test P._mpow(P._monomial(1 => 2), 0) == P._Mono()
      @test P._mexp(P._monomial(1 => 1//2), 1) == 1//2
      @test P._mexp(P._monomial(1 => 1//2), 2) == 0
  end

  @testset "polynomials add, subtract and multiply exactly" begin
      P = EnzymeRates
      x, y = P._term(P._monomial(1 => 1)), P._term(P._monomial(2 => 1))
      p = x + y
      @test (p * p).terms == Dict(P._monomial(1 => 2) => 1, P._monomial(1 => 1, 2 => 1) => 2,
                                  P._monomial(2 => 2) => 1)
      @test P._nterms(p * p) == 3
      @test (x - y) * (x + y) == x * x - y * y
      @test p - y == x && -x == P._scale(x, -1)
      @test iszero(p - p) && iszero(zero(P._Poly)) && !iszero(one(P._Poly))
      @test p * P._monomial(1 => -1) == one(P._Poly) + P._term(P._monomial(1 => -1, 2 => 1))
      @test P._scale(p, 0) == zero(P._Poly) && P._term(P._Mono(), 0) == zero(P._Poly)
      @test P._addto!(P._term(P._monomial(1 => 1), 2), P._monomial(1 => 1), -2) ==
            zero(P._Poly)
      # Equal by terms, hashed by terms, and distinct objects.
      @test x + y == y + x && hash(x + y) == hash(y + x) && x + y !== y + x
      @test ismutabletype(P._Poly) && isconst(P._Poly, :terms)
      @test fieldtypes(P._Poly) == (Dict{P._Mono, Int},)
  end

  @testset "an Int64 overflow throws OverflowError" begin
      P = EnzymeRates
      m = P._monomial(1 => 1)
      @test_throws OverflowError P._term(m, typemax(Int)) + P._term(m, 1)
      @test_throws OverflowError P._term(m, typemax(Int) ÷ 2 + 1) * P._term(P._Mono(), 2)
      @test_throws OverflowError P._scale(P._term(m, typemax(Int)), 2)
  end

  @testset "_divexact divides out the Haldane factor" begin
      P = EnzymeRates
      # Variables: 1 = A, 2 = P, 3 = Keq.
      mono(ps...) = P._term(P._monomial(ps...))
      H = mono(1 => 1) - mono(2 => 1, 3 => -1)
      q = mono(1 => 1, 3 => 1//2) + P._scale(mono(2 => 1, 3 => -1//2), 2) + one(P._Poly)
      @test P._divexact(H * q, H) == q
      # A substrate with two copies: H = A² − P/Keq.
      H2 = mono(1 => 2) - mono(2 => 1, 3 => -1)
      q2 = mono(1 => 1) + mono(2 => 1) + mono(1 => 1, 2 => 1, 3 => 1//2)
      @test P._divexact(H2 * q2, H2) == q2
      @test iszero(P._divexact(zero(P._Poly), H))
      @test_throws ArgumentError P._divexact(H * q + one(P._Poly), H)
      @test_throws ArgumentError P._divexact(q, q)
  end

  @testset "_subst replaces variables by power products" begin
      P = EnzymeRates
      mono(ps...) = P._term(P._monomial(ps...))
      p = mono(1 => 1, 2 => 2) + P._scale(mono(3 => 1), 3)
      f(v) = v == 2 ? P._monomial(4 => 1//2, 5 => -1) : nothing
      @test P._subst(f, p) == mono(1 => 1, 4 => 1, 5 => -2) + P._scale(mono(3 => 1), 3)
      # Terms that substitution makes alike merge, and cancel to zero.
      g(v) = v == 2 ? P._monomial(1 => 1) : nothing
      @test P._subst(g, mono(1 => 1) + mono(2 => 1)) == P._scale(mono(1 => 1), 2)
      @test iszero(P._subst(g, mono(1 => 1) - mono(2 => 1)))
  end

  @testset "the shifts bring polynomials to lowest terms in the concentrations" begin
      P = EnzymeRates
      t = P._Vars()
      a, b, k = P._var!(t, :A; conc = true), P._var!(t, :B; conc = true), P._var!(t, :k)
      @test (a, b, k) == (1, 2, 3) && P._var!(t, :A; conc = true) == 1
      @test t.names == [:A, :B, :k] && t.conc == [true, true, false]
      @test_throws ErrorException P._var!(t, :k; conc = true)
      @test t == P._Vars([:A, :B, :k], BitVector([1, 1, 0]), Dict(:A => 1, :B => 2, :k => 3))
      mono(ps...) = P._term(P._monomial(ps...))
      Q = mono(a => 1, b => 1, k => 1) + mono(a => 2, k => -1)
      D = mono(a => 1) + mono(a => 1, b => 1)
      # A is in every term, lowest 1; B is missing from some, lowest 0; k reaches −1.
      @test P._lowest_shift((Q, D)) == P._monomial(a => -1, k => 1)
      @test P._conc_shift((Q, D), t) == P._monomial(a => -1)
      D2 = mono(a => 1, k => 1) + mono(a => 1, b => 1) + P._scale(mono(k => 1//2), 2) +
           mono(b => 1, k => 1)
      @test P._split_conc(D2, t) == Dict(P._monomial(a => 1) => mono(k => 1),
                                        P._monomial(a => 1, b => 1) => one(P._Poly),
                                        P._Mono() => P._scale(mono(k => 1//2), 2),
                                        P._monomial(b => 1) => mono(k => 1))
  end

  @testset "_content splits off the integer and monomial content" begin
      P = EnzymeRates
      mono(ps...) = P._term(P._monomial(ps...))
      p = P._scale(mono(1 => 2, 3 => 1), 6) + P._scale(mono(1 => 1, 3 => 1//2), 4)
      @test P._content(p) ==
            (2, P._monomial(1 => 1, 3 => 1//2),
             P._scale(mono(1 => 1, 3 => 1//2), 3) + P._scale(one(P._Poly), 2))
      @test P._content(P._scale(mono(1 => 1), -2) + P._scale(one(P._Poly), -4)) ==
            (2, P._Mono(), -mono(1 => 1) - P._scale(one(P._Poly), 2))
      @test P._content(P._scale(mono(1 => -1, 2 => 1), 5)) ==
            (5, P._monomial(1 => -1, 2 => 1), one(P._Poly))
  end
  ```

  In `test/runtests.jl`, replace
  ```julia
      include("test_rate_eq_derivation.jl")
      include("test_oracles.jl")
  ```
  with
  ```julia
      include("test_rate_eq_derivation.jl")
      include("test_poly.jl")
      include("test_oracles.jl")
  ```

- [ ] **Step 2: Watch it fail.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T2_poly.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_poly.jl")'
  ```
  Expected: the first testset, "monomials multiply, raise and read exponents exactly",
  records nine errors, each an `UndefVarError` naming a `poly.jl` function (`_monomial`,
  `_mmul`, `_mpow` or `_mexp`) not defined in `EnzymeRates`.
  A top-level testset that errors throws when it ends, so the include stops with `Some
  tests did not pass: 0 passed, 0 failed, 9 errored, 0 broken`, the other six testsets do
  not run, and the log ends `exit 1`.

- [ ] **Step 3: Write `src/poly.jl`.**

  ```julia
  # ABOUTME: Exact sparse polynomials over integer-indexed variables: concentrations with
  # ABOUTME: integer exponents, constants with rational exponents, checked Int64 coefficients.

  """
  A variable table: variable `i` is `names[i]`, a concentration when `conc[i]`; `index` maps
  each name to its variable. Two tables are equal when they name the same variables in the
  same order.
  """
  struct _Vars
      names::Vector{Symbol}
      conc::BitVector
      index::Dict{Symbol, Int}
  end
  _Vars() = _Vars(Symbol[], falses(0), Dict{Symbol, Int}())
  Base.:(==)(a::_Vars, b::_Vars) = a.names == b.names && a.conc == b.conc
  Base.hash(t::_Vars, h::UInt) = hash(t.conc, hash(t.names, hash(:_Vars, h)))

  """The index of variable `s` in `t`, added when new; `conc` marks a concentration. Errors
  when `s` is already a variable of the other kind."""
  function _var!(t::_Vars, s::Symbol; conc::Bool = false)
      i = get(t.index, s, 0)
      i == 0 || return t.conc[i] == conc ? i :
          error("internal error: $s is a concentration in one use and a constant in another")
      push!(t.names, s)
      push!(t.conc, conc)
      t.index[s] = length(t.names)
  end

  """An exponent: an integer for a concentration, a rational for a constant."""
  const _Exp = Rational{Int}

  """A monomial: `variable => exponent` pairs sorted by variable, no zero exponent."""
  const _Mono = Vector{Pair{Int, _Exp}}

  """The product of two monomials (a merge of the sorted pairs)."""
  function _mmul(a::_Mono, b::_Mono)
      r = _Mono()
      i, j = 1, 1
      while i <= length(a) || j <= length(b)
          if j > length(b) || (i <= length(a) && a[i].first < b[j].first)
              push!(r, a[i])
              i += 1
          elseif i > length(a) || b[j].first < a[i].first
              push!(r, b[j])
              j += 1
          else
              e = a[i].second + b[j].second
              iszero(e) || push!(r, a[i].first => e)
              i += 1
              j += 1
          end
      end
      r
  end

  """`a` to the power `e`."""
  _mpow(a::_Mono, e::Union{Integer, _Exp}) =
      iszero(e) ? _Mono() : _Mono([v => x * e for (v, x) in a])

  """The monomial of the `variable => exponent` pairs `ps`; a variable may repeat."""
  _monomial(ps::Pair{Int}...) =
      foldl(_mmul, (_Mono([v => _Exp(e)]) for (v, e) in ps if e != 0); init = _Mono())

  """The exponent of variable `v` in `a`."""
  function _mexp(a::_Mono, v::Int)
      i = findfirst(p -> p.first == v, a)
      i === nothing ? zero(_Exp) : a[i].second
  end

  """
  A polynomial: monomial => nonzero coefficient. Every sum and product is checked, so a
  coefficient that leaves `Int64` throws `OverflowError` instead of wrapping. It is mutable,
  with a `const` field, so a loop that merges it with another type into a `Union` holds a
  reference, never an inline box the compiler may hoist. Two polynomials are equal when
  their terms are.
  """
  mutable struct _Poly
      const terms::Dict{_Mono, Int}
  end
  _Poly() = _Poly(Dict{_Mono, Int}())

  """Add `c·m` into `p`, dropping a zero sum; returns `p`."""
  function _addto!(p::_Poly, m::_Mono, c::Integer)
      v = Base.checked_add(get(p.terms, m, 0), Int(c))
      iszero(v) ? delete!(p.terms, m) : (p.terms[m] = v)
      p
  end

  """The one-term polynomial `c·m`."""
  _term(m::_Mono, c::Integer = 1) = _addto!(_Poly(), m, c)

  Base.one(::Type{_Poly}) = _term(_Mono())
  Base.zero(::Type{_Poly}) = _Poly()
  Base.iszero(p::_Poly) = isempty(p.terms)
  Base.:(==)(p::_Poly, q::_Poly) = p.terms == q.terms
  Base.hash(p::_Poly, h::UInt) = hash(p.terms, hash(:_Poly, h))

  """The number of terms of `p`."""
  _nterms(p::_Poly) = length(p.terms)

  function Base.:+(p::_Poly, q::_Poly)
      r = _Poly(copy(p.terms))
      for (m, c) in q.terms
          _addto!(r, m, c)
      end
      r
  end
  Base.:-(p::_Poly) = _scale(p, -1)
  Base.:-(p::_Poly, q::_Poly) = p + (-q)

  function Base.:*(p::_Poly, q::_Poly)
      r = _Poly()
      for (m1, c1) in p.terms, (m2, c2) in q.terms
          _addto!(r, _mmul(m1, m2), Base.checked_mul(c1, c2))
      end
      r
  end
  # Multiplying by a monomial maps distinct monomials to distinct ones, so no terms merge.
  Base.:*(p::_Poly, m::_Mono) =
      _Poly(Dict{_Mono, Int}(_mmul(k, m) => c for (k, c) in p.terms))

  """`c·p` for an integer `c`."""
  _scale(p::_Poly, c::Integer) = iszero(c) ? _Poly() :
      _Poly(Dict{_Mono, Int}(m => Base.checked_mul(x, Int(c)) for (m, x) in p.terms))

  """
  `p` with each variable `v` for which `f(v)` returns a monomial (a power product) replaced
  by it, `v^e` by `f(v)^e`; `f(v) === nothing` keeps `v`.
  """
  function _subst(f, p::_Poly)
      r = _Poly()
      for (m, c) in p.terms
          out = _Mono()
          for (v, e) in m
              s = f(v)
              out = _mmul(out, s === nothing ? _Mono([v => e]) : _mpow(s, e))
          end
          _addto!(r, out, c)
      end
      r
  end

  """
      _divexact(p, d) -> p / d

  The exact quotient of `p` by a two-term `d` (the Haldane factor); an `ArgumentError` when
  `d` does not divide `p`. With `d = c_u·u + c_w·w` and `t = w/u`, the terms of `p/u` fall
  into chains `key·t^k` over integer `k`; dividing a chain by `c_u + c_w·t` is the recurrence
  `b_k = (a_k − c_w·b_{k−1}) / c_u`, which must end at zero.
  """
  function _divexact(p::_Poly, d::_Poly)
      _nterms(d) == 2 || throw(ArgumentError("_divexact divides by a two-term polynomial"))
      (u, cu), (w, cw) = collect(d.terms)
      t = _mmul(w, _mpow(u, -1))
      v, et = first(t)
      chains = Dict{_Mono, Vector{Pair{Int, Int}}}()
      for (m, c) in p.terms
          m = _mmul(m, _mpow(u, -1))
          k = floor(Int, _mexp(m, v) / et)
          push!(get!(chains, _mmul(m, _mpow(t, -k)), Pair{Int, Int}[]), k => c)
      end
      q = _Poly()
      for (key, chain) in chains
          sort!(chain; by = first)
          b, j = 0, 1
          for k in chain[1].first:chain[end].first
              a = 0
              if chain[j].first == k
                  a = chain[j].second
                  j += 1
              end
              r = Base.checked_sub(a, Base.checked_mul(cw, b))
              rem(r, cu) == 0 ||
                  throw(ArgumentError("_divexact: the quotient is not integral"))
              b = div(r, cu)
              k == chain[end].first || iszero(b) || _addto!(q, _mmul(key, _mpow(t, k)), b)
          end
          iszero(b) || throw(ArgumentError("_divexact: the divisor does not divide"))
      end
      q
  end

  """
      _lowest_shift(ps, counted = every variable) -> _Mono

  The monomial in the `counted` variables that, multiplied into every polynomial of `ps`,
  brings each such variable's lowest exponent over all their terms to zero, a term without
  the variable counting as exponent 0: the shift to lowest terms.
  """
  function _lowest_shift(ps, counted = _ -> true)
      low, seen, n = Dict{Int, _Exp}(), Dict{Int, Int}(), 0
      for p in ps, m in keys(p.terms)
          n += 1
          for (v, e) in m
              counted(v) || continue
              low[v] = min(get(low, v, e), e)
              seen[v] = get(seen, v, 0) + 1
          end
      end
      _monomial((v => -(seen[v] < n ? min(e, 0) : e) for (v, e) in low)...)
  end

  """The shift to lowest terms in the concentrations of `t` alone (`_lowest_shift`)."""
  _conc_shift(ps, t::_Vars) = _lowest_shift(ps, v -> t.conc[v])

  """`p` as `concentration monomial => its coefficient`, a polynomial in the constants."""
  function _split_conc(p::_Poly, t::_Vars)
      out = Dict{_Mono, _Poly}()
      for (m, c) in p.terms
          _addto!(get!(_Poly, out, filter(q -> t.conc[q.first], m)),
                  filter(q -> !t.conc[q.first], m), c)
      end
      out
  end

  """
  The nonzero polynomial `p` as `c·mu·rest`: its integer content `c` (the coefficients'
  gcd, positive, or a one-term polynomial's own coefficient), its monomial content `mu`
  (each variable's lowest exponent over the terms, a term without the variable counting as
  exponent 0) and the rest, `one(_Poly)` for a one-term polynomial.
  """
  function _content(p::_Poly)
      c = foldl(gcd, values(p.terms))
      mu = _monomial((v => minimum(mo -> _mexp(mo, v), keys(p.terms))
                      for v in unique(first(t) for mo in keys(p.terms) for t in mo))...)
      inv_mu = _mpow(mu, -1)
      c, mu, _Poly(Dict{_Mono, Int}(_mmul(mo, inv_mu) => div(x, c) for (mo, x) in p.terms))
  end
  ```

  In `src/EnzymeRates.jl`, replace
  ```julia
  include("thermo.jl")
  include("admission.jl")
  ```
  with
  ```julia
  include("thermo.jl")
  include("poly.jl")
  include("admission.jl")
  ```
  Phase 4's `precompile` lines at the end of the module stay where FW-10 put them.

- [ ] **Step 4: Watch it pass.** Rerun Step 2's command. Expected: seven testsets, every
  test passing (41 on a9e0699 with this code), `exit 0`.

- [ ] **Step 5: Amend the spec (P5-1).** In
  `docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md`, section 9,
  subsection "Engine (C1)", replace the bullet
  ```markdown
  - Exact polynomial arithmetic: checked `Int64`, retried in `BigInt` on overflow.
  ```
  with
  ```markdown
  - Exact polynomial arithmetic: checked `Int64`. Under the cap no coefficient leaves
    `Int64` (proof in the PR), so an overflow is an internal error.
  ```
  The PR body lists this amendment for Denis (Task 15, ◆ P5-1).

- [ ] **Step 6: Write the GC gate's environment script and workload.** Create
  `$S/phase45/gcscan5/mkenv.jl`:

  ```julia
  # ABOUTME: Creates an environment that develops a worktree by path and holds the optimizer
  # ABOUTME: the identify workload needs, so the phase-5 GC scans never touch its manifest.
  using Pkg
  Pkg.offline(true)
  Pkg.activate(ARGS[2])
  Pkg.develop(path = ARGS[1])
  Pkg.add(["Random", "OptimizationCMAEvolutionStrategy"])
  Pkg.instantiate()
  using EnzymeRates
  println("ok ", pkgversion(EnzymeRates), " from ", pathof(EnzymeRates))
  ```

  and `$S/phase45/gcscan5/scan_law.jl`:

  ```julia
  # ABOUTME: The phase-5 GC gate's typed workloads (P5-7), run with LICM remarks, --pkgimages=no:
  # ABOUTME: `poly` merges fresh polynomials with other types; `poly control` on the crash layout.
  # Usage: julia --project=<env> --pkgimages=no scan_law.jl poly [control]
  using EnzymeRates
  const ER = EnzymeRates
  const WORK = ARGS[1]
  const CONTROL = length(ARGS) > 1 && ARGS[2] == "control"

  # The control layout, which phase 3's GC probe crashed: an immutable wrapper of one
  # reference, the Dict of a polynomial.
  struct ImmPoly
      terms::Dict{ER._Mono, Int}
  end
  wrap(::Type{ER._Poly}, p::ER._Poly) = p
  wrap(::Type{ImmPoly}, p::ER._Poly) = ImmPoly(p.terms)
  nterms(p::ER._Poly) = ER._nterms(p)
  nterms(p::ImmPoly) = length(p.terms)

  # Law-like holders over a polynomial type P, as the slice's Law and AlloLaw held their
  # polynomials: a loop over their Union reads fields whose type is P or a holder.
  struct ToyLaw{P}
      H::P
      Q::P
      D::P
  end
  struct ToyAllo{P}
      H::P
      active::ToyLaw{P}
      sites::Vector{Tuple{P, P, Int}}
  end

  """A fresh law of polynomials over 4 variables, built by the package's checked arithmetic;
  an odd `i` gives a plain law, an even one an allosteric holder."""
  @noinline function fresh(::Type{P}, i::Int) where {P}
      x, y, z = (ER._term(ER._monomial(v => 1)) for v in 1:3)
      one_ = one(ER._Poly)
      H = x - y * ER._monomial(4 => -1)
      D = (one_ + x) * (one_ + y + z) * (one_ + ER._scale(x, i % 5 + 1))
      Q = ER._divexact(H * D, H)
      L = ToyLaw{P}(wrap(P, H), wrap(P, Q), wrap(P, D))
      isodd(i) ? L : ToyAllo{P}(wrap(P, H), L, [(wrap(P, D), wrap(P, Q), i % 3 + 1)])
  end

  function merge_vec(ls)
      t = 0
      for L in ls
          y = L isa ToyLaw ? L.H : L.sites                   # Union{P, Vector}
          t += y isa Vector ? length(y) : nterms(y)
      end
      t
  end
  function merge_nothing(ls)
      t = 0
      for (i, L) in enumerate(ls)
          z = isodd(i) ? L.H : nothing                       # Union{Nothing, P}
          t += z === nothing ? 0 : nterms(z)
      end
      t
  end
  function fresh_or_nothing(::Type{P}, n::Int, full::Bool) where {P}
      t = 0
      for i in 1:n
          z = isodd(i) ? fresh(P, i) : nothing               # Union{Nothing, holders}
          t += z === nothing ? 0 : z isa ToyLaw ? nterms(z.D) : nterms(z.active.D)
          full && GC.gc(true)
      end
      t
  end
  function fresh_or_vec(::Type{P}, n::Int, full::Bool) where {P}
      t = 0
      for i in 1:n
          L = fresh(P, i)
          y = L isa ToyLaw ? L : L.sites                     # Union{ToyLaw, Vector}
          t += y isa ToyLaw ? nterms(y.Q) : length(y)
          full && GC.gc(true)
      end
      t
  end
  function fresh_active_or_poly(::Type{P}, n::Int, full::Bool) where {P}
      t = 0
      for i in 1:n
          L = fresh(P, i)
          y = L isa ToyAllo ? L.active : L.H                 # Union{ToyLaw, P}
          t += y isa ToyLaw ? nterms(y.D) : nterms(y)
          full && GC.gc(true)
      end
      t
  end
  function fresh_poly_or_vec(::Type{P}, n::Int, full::Bool) where {P}
      t = 0
      for i in 1:n
          L = fresh(P, i)
          y = L isa ToyLaw ? L.H : L.sites                   # Union{P, Vector}
          t += y isa Vector ? length(y) : nterms(y)
          full && GC.gc(true)
      end
      t
  end

  function poly_rounds(::Type{P}) where {P}
      ref = nothing
      for rd in 1:3
          ls = Union{ToyLaw{P}, ToyAllo{P}}[fresh(P, i) for i in 1:400]
          sig = (merge_vec(ls), merge_nothing(ls),
                 fresh_or_nothing(P, 400, false), fresh_or_vec(P, 400, false),
                 fresh_active_or_poly(P, 400, false), fresh_poly_or_vec(P, 400, false),
                 fresh_or_nothing(P, 60, true), fresh_or_vec(P, 60, true),
                 fresh_active_or_poly(P, 60, true), fresh_poly_or_vec(P, 60, true))
          ref = something(ref, sig)
          GC.gc(true)
          println("round ", rd, ": ", sig == ref ? "same" : "DIFFERENT", " ", sig,
                  ", full gc ok")
          flush(stdout)
      end
  end

  t0 = time()
  if WORK == "poly"
      P = CONTROL ? ImmPoly : ER._Poly
      println("julia ", VERSION, ", workload poly on ", P, "; ismutabletype ", ismutabletype(P))
      poly_rounds(P)
  else
      error("unknown workload $WORK")
  end
  println("done in ", round(time() - t0; digits = 1), " s")
  ```

  The six loops follow `$S/phase45/research5/synth/scan_f.jl`. Its law is replaced by toy
  holders, since `_Law` arrives in Task 3. As in research variants E (mutable `Poly`) and C
  (immutable `Poly`), the holders are immutable and only the polynomial's layout differs.
  Task 6 extends this file with the package's laws.

- [ ] **Step 7: Run the GC scan on Julia 1.12.5 and 1.13.1.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; G=$S/phase45/gcscan5; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T2_env112.log 1800 $G/mkenv.jl $PWD $G/env112 && $S/phase3/jl.sh $S/phase45/logs/P5_T2_env113.log 1800 +1.13 $G/mkenv.jl $PWD $G/env113
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; G=$S/phase45/gcscan5; L=$S/phase45/logs; cd $S/phase45/law_wt
  for v in 112 113; do ch=$([ $v = 113 ] && echo +1.13)
    for w in poly control; do arg=$([ $w = control ] && echo "poly control" || echo poly)
      JULIA_LLVM_ARGS="-pass-remarks=julia-licm -pass-remarks-analysis=julia-licm" \
        $S/phase3/jl.sh $L/P5_T2_scan${v}_$w.log 1800 $ch --project=$G/env$v --pkgimages=no --heap-size-hint=2000M $G/scan_law.jl $arg
    done
    JULIA_LLVM_ARGS="-pass-remarks=julia-licm -pass-remarks-analysis=julia-licm" \
      $S/phase3/jl.sh $L/P5_T2_scan${v}_repro.log 600 $ch --startup-file=no --pkgimages=no $S/phase3/gcprobe/repro_upstream.jl
  done
  python3 -I $S/phase3/gcprobe/hoistscan.py $L/P5_T2_scan11*.log > $L/P5_T2_hoist.log 2>&1; cat $L/P5_T2_hoist.log
  ```
  Each scan recompiles the package and its dependencies with the remarks on, so run the
  loop in the background and poll. Expected:
  - Each environment log prints `ok 0.10.0 from <path>/law_wt/src/EnzymeRates.jl` and
    ends `exit 0`.
  - Each `P5_T2_scan<v>_poly.log` prints three lines
    `round <r>: same (600, 400, 1800, 2000, 2200, 600, 270, 300, 330, 90), full gc ok`,
    then `done in <t> s` and `exit 0`.
  - Each `_control.log` and each `_repro.log` prints `GC error (probable corruption)` and
    `signal 6`, and ends `exit 134`: the hazard itself.
  - The scanner prints `hoisted allocations: 0  hoisted with a bits store into a
    reference field: 0` for each `_poly` log.
  - For each `_control` and `_repro` log it prints a nonzero bits-store figure and a
    `DANGER` line.

  A stand-in run on Julia 1.12.5 (`$S/phase45/plan5/probes/p1_scan_check.sh`) used this
  workload with these polynomials in a stand-in module. It found 0 of 44 for `_Poly`, and
  10 of 50 with a bits store for the control, all at `range.jl:921`, the inlined loop of
  the merging functions. The control then aborted the collector. On 1.13.1 the same run
  found 0 of 55 and 10 of 65, at `range.jl:928`, and the control aborted
  (`$S/phase45/refresh5/settle_7.log`, `settle_8.log`).

  If a `_poly` log shows a `DANGER` line, stop: `_Poly`'s layout must change before any
  law holds it (R18; P5-1). If a `_control` log shows none, the workload cannot tell the
  layouts apart. Stop and report the logs, since the gate needs a control that fires.

- [ ] **Step 8: Run the full suite** (log `$S/phase45/logs/P5_T2_full.log`, background,
  polled). Expected: green and pristine; the pass count is Task 1's plus 41; wall time and
  maxrss within noise of Task 0's.

- [ ] **Step 9: Commit.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && git add src/poly.jl src/EnzymeRates.jl test/test_poly.jl test/runtests.jl docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md && git commit -F - <<'EOF'
  Add exact sparse polynomials in checked Int64

  poly.jl holds the polynomials the law engine works in: monomials with
  rational exponents over a variable table, coefficients in checked Int64,
  exact division by the Haldane factor, substitution, the shifts to lowest
  terms and the content. _Poly is mutable with a const field. The GC scan of
  six loops that merge fresh polynomials with other types finds no hoisted
  allocation with a bits store on Julia 1.12.5 and 1.13.1. The same loops on
  an immutable one-field wrapper hoist <n112> and <n113> and abort the
  collector. The spec's engine bullet drops the BigInt retry (P5-1).

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
  EOF
  ```
  Replace `<n112>` and `<n113>` with the control's counts from `P5_T2_hoist.log`.

---

### Task 3: Plain engine, not wired; `_eq_complexity` by checked Bareiss (P5-3, P5-5)

`engine.jl` derives one conformation's law, `_Law`, on an edge mask. No consumer reads it
yet; `rate_equation`, the printer and identify still read the old engine.
`_eq_complexity` moves here and counts in checked `Int`.

**Files:**
- Create: `src/engine.jl`; modify `src/EnzymeRates.jl`
- Modify: `src/rate_eq_derivation.jl` (move `_eq_complexity`, `RateEquationTooLargeError`
  with its `showerror`, and `_assert_derivable` to `engine.jl`)
- Modify: `ROADMAP.md` ("`BigInt` arithmetic inflates the GC's live-byte count": drop the
  `_eq_complexity` clause, after the probe measures 0 KiB)
- Create: `test/test_law.jl`; modify `test/runtests.jl`; modify the fixtures file
  (`_testhelper_poly_value` on `_Poly`, `_testhelper_evaluate_law`,
  `_testhelper_laplacian_tree_count`)
- Create (scratch): `$S/phase45/probes/P5_T3_live_bytes.jl`
- Logs: `P5_T3_{live_before,law,compare_laws,compare_complexity,compare_dterms,live_after,full}.log`

**Interfaces:**
- Consumes: Task 1's `_solve`, `_Solve`, `_Dependent`; Task 2's poly API and `_Vars`
  equality; `_segments(m, mask, _flags(m))`, `_flat_edges`, `_is_binding`, `copies`,
  `metabolites`, `substrates`, `products`, `m.table[i].name` (a `::Mod` copy reads as its
  reactant), `det` (`LinearAlgebra`, already imported); in the tests, `catalytic`,
  `_active_mask` and `_inactive_mask` (`_law_mask`, on `_exclusion_components`).
- Produces:
  ```julia
  mutable struct _Law          # v = E_total·H·Q/D in one conformation
      const vars::_Vars
      const H::_Poly           # ΠS − ΠP/Keq
      const Q::_Poly           # positive coefficients
      const D::_Poly           # shift·Σ_g sigma[g]·trees[g]
      const sigma::Vector{_Poly}   # each segment's weight sum
      const trees::Vector{_Poly}   # each segment's spanning-tree sum
      const shift::_Mono           # concentration monomial to lowest terms, in Q and D
  end
  Base.:(==)(a::_Law, b::_Law); Base.hash(L::_Law, h::UInt)      # over the fields
  _spanning_trees(n::Int, pairs::Vector{Tuple{Int, Int}})::Vector{Vector{Int}}
  _bare_law(m::Mechanism, mask::BitVector, concentrations::Vector{Symbol})::
      @NamedTuple{vars::_Vars, cols::Vector{Tuple{Int, Int}}, num::_Poly, den::_Poly,
                  sigma::Vector{_Poly}, trees::Vector{_Poly}}
  _law(m::Mechanism, s::_Solve, mask::BitVector; inactive::Bool = false,
       concentrations::Vector{Symbol} = metabolites(m))::_Law
  _law(m::Mechanism)::_Law                    # _solve(m), every edge
  _tree_count(lap::Matrix{Int})::Int          # spanning trees of a Laplacian: checked
                                              #   Bareiss; BigInt det after overflow,
                                              #   saturated at typemax(Int)
  _eq_complexity(m::Mechanism, mask::BitVector)::Int   # V × τ, on _tree_count
  _eq_complexity(m::Union{Mechanism, AllostericMechanism})::Int  # whole catalytic graph
                                                                 #   until Task 4
  struct RateEquationTooLargeError <: Exception end    # moved, text unchanged
  _assert_derivable(m::Mechanism, mask::BitVector)::Nothing     # moved
  ```
  - `_law` runs `_assert_derivable(m, mask)` first.
  - It reads the segments from `_segments`, so each is rooted at its lowest form and free
    E (form 1) roots its own. The forms join the tree in `_segments`' order, each by the
    first kept RE edge in edge (content) order from a form placed before it, as spec
    section 9 asks ("a content-ordered tree").
  - Each RE edge must close its cycle with no net uptake of any metabolite (P5-5).
    Otherwise `_law` errors "internal error: a rapid-equilibrium cycle of this mechanism
    turns the reaction over (class 2)".
  - It then divides by H with `_divexact` and shifts to lowest terms.
  - An `OverflowError` becomes "internal error: the rate law of this mechanism leaves
    Int64 arithmetic, which the term cap rules out:" and the mechanism's text.
  - It ports `SLC/engine.jl`'s `bare_law` and `law` (RES §1.4) onto the graph types.
- Test helpers (fixtures file, new section "Phase-5 law oracles"):
  ```julia
  _testhelper_poly_value(p::EnzymeRates._Poly, vars::EnzymeRates._Vars,
                         env::AbstractDict{Symbol})::BigFloat
  _testhelper_evaluate_law(L::EnzymeRates._Law, env::AbstractDict{Symbol})::BigFloat  # H·Q/D
  _testhelper_laplacian_tree_count(m::EnzymeRates.Mechanism, mask::BitVector)::BigInt
  ```

**Gate:**
- `P5_compare.jl $S/phase45/p5_base laws`: every row read through `_law` (each plain
  mechanism of the baseline) equals Task 0's value within 1e-20 relative.
- `P5_compare.jl … complexity`: the new `_eq_complexity` equals the head's on all 480
  baseline mechanisms, and on each allosteric one's active and inactive masks (RES §4's
  "`complexity` = `_eq_complexity`"). `P5_compare.jl … dterms`: the six gate laws'
  `_nterms(_law(m).D)` read 396, 396, 396, 589, 799, 858, the head's counts (P5-23).
- `test_law.jl` is green:
  - the oracle agrees within 1e-20 on every plain spec;
  - every Q coefficient is positive, H divides, and (Q, D) is in lowest terms;
  - detailed balance holds;
  - rewrites and a same-model regrouping give `==` laws;
  - an `:OnlyI` fixture's two conformations, each on its own mask, give the laws of the
    steps each keeps, over one variable table (P5-8);
  - `_eq_complexity` equals the test's own Laplacian count on every plain spec (31);
  - the engine refuses RE turnover, a law over the cap, and an overflow, each with its
    message.
- The probe finds 0.00 KiB of live bytes per `_eq_complexity` call after a full
  collection (RES: 0.39 before).
- The full suite is green; +3 s (the a9e0699 run of `test_law.jl` took 2.8 s).

**Testsets added (`test/test_law.jl`):** "_spanning_trees enumerates each spanning tree
once", "_tree_count is the BigInt determinant, saturated at typemax(Int)",
"_eq_complexity is segments times spanning trees", "the plain law equals the exact steady
state on every plain spec", "a plain law is H times a positive polynomial in lowest terms",
"rewriting or regrouping a mechanism changes no plain law", "plain laws satisfy detailed
balance", "the engine refuses a rapid-equilibrium cycle that turns the reaction over",
"_law refuses a law over the term cap", and two more: "a
conformation's law on its mask is the law of the steps it keeps" (the `mask`, `inactive`
and `concentrations` paths `_law` adds in this task, on P5-8's `:OnlyI` fixture
`competitive_RE`, written inline) and "an Int64 overflow in a law is an internal error
naming the mechanism".

**Pins changed:** none. TR "_eq_complexity (V×τ term-count estimate)" keeps its values; only
the return type changes, from `BigInt` to `Int`. TR "Rate equation too large error" keeps
its type and text (RES §3.4 row 14).

- [ ] **Step 1: Measure `_eq_complexity`'s live bytes on the head.** Create
  `$S/phase45/probes/P5_T3_live_bytes.jl`:

  ```julia
  # ABOUTME: Measures the live bytes a full collection leaves per _eq_complexity call, on the
  # ABOUTME: ROADMAP's plain random-order bi-bi of nine steady-state steps (V×τ = 336).
  using EnzymeRates
  m = EnzymeRates.Mechanism(@enzyme_mechanism begin
      substrates: A, B
      products: P, Q
      steps: begin
          E + A <--> E(A)
          E + B <--> E(B)
          E(A) + B <--> E(A, B)
          E(B) + A <--> E(A, B)
          E(A, B) <--> E(P, Q)
          E(P, Q) <--> E(Q) + P
          E(P, Q) <--> E(P) + Q
          E(P) <--> E + P
          E(Q) <--> E + Q
      end
  end)
  mask = trues(length(m.edges))
  function churn(m, mask, n)
      t = 0
      for _ in 1:n
          t += EnzymeRates._eq_complexity(m, mask)
      end
      t
  end
  churn(m, mask, 200)
  GC.gc(true); GC.gc(true)
  before = Base.gc_live_bytes()
  total = churn(m, mask, 20_000)
  GC.gc(true); GC.gc(true)
  after = Base.gc_live_bytes()
  println("_eq_complexity = ", EnzymeRates._eq_complexity(m, mask), " (",
          typeof(EnzymeRates._eq_complexity(m, mask)), "); live bytes per call after a ",
          "full collection: ", round((after - before) / 20_000 / 1024; digits = 2), " KiB")
  ```
  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T3_live_before.log 600 --project=. $S/phase45/probes/P5_T3_live_bytes.jl && head -1 $S/phase45/logs/P5_T3_live_before.log
  ```
  Expected:
  `_eq_complexity = 336 (BigInt); live bytes per call after a full collection: <x> KiB`,
  with `<x>` clearly above 0 (RES measured 0.39 KiB; the ROADMAP gives 15 KiB on its own
  loop). Record `<x>`.

- [ ] **Step 2: Add the law oracles to the fixtures file.** In
  `test/mechanism_definitions_for_test_enzyme_derivation.jl`, insert before the line
  `# ── Mechanism test specifications ───────────────────────────────────────────`:

  ```julia
  # ── Phase-5 law oracles ─────────────────────────────────────────────────────

  """The value of `p`, a polynomial over the variables `vars`, at the values `env` (by name),
  in 256-bit BigFloat; every power, a fractional one included, is correctly rounded to that
  precision."""
  function _testhelper_poly_value(p::EnzymeRates._Poly, vars::EnzymeRates._Vars,
                                  env::AbstractDict{Symbol})
      setprecision(BigFloat, 256) do
          pow(x, e) = x^(BigFloat(numerator(e)) / denominator(e))
          sum((c * prod((pow(BigFloat(env[vars.names[v]]), e) for (v, e) in mo);
                        init = one(BigFloat))
               for (mo, c) in p.terms); init = zero(BigFloat))
      end
  end

  """The rate per enzyme of the law `L` (`EnzymeRates._law`), `H·Q/D`, at the values `env`
  (by name: the concentrations, the independent constants and `Keq`), in 256-bit BigFloat."""
  _testhelper_evaluate_law(L::EnzymeRates._Law, env::AbstractDict{Symbol}) =
      setprecision(BigFloat, 256) do
          _testhelper_poly_value(L.H, L.vars, env) *
              _testhelper_poly_value(L.Q, L.vars, env) /
              _testhelper_poly_value(L.D, L.vars, env)
      end

  """
  V×τ of the edges of the plain mechanism `m` that `mask` keeps, counted by the test from its
  own reading of the graph (`_testhelper_re_segments`): the segments times the spanning
  trees of the multigraph whose edges are the kept steady-state steps between segments,
  parallel steps counted apart, by the matrix-tree theorem in `BigInt`. One segment counts
  1, an empty mask 0.
  """
  function _testhelper_laplacian_tree_count(m::EnzymeRates.Mechanism, mask::BitVector)
      g = _testhelper_graph(m)
      kept = BitVector([mask[k] for k in EnzymeRates._flat_edges(m)])
      segments, _ = _testhelper_re_segments(g, kept)
      V = length(segments)
      V <= 1 && return big(V)
      seg = zeros(Int, length(g.forms))
      for (i, s) in enumerate(segments), f in s
          seg[f] = i
      end
      lap = zeros(BigInt, V, V)
      for (j, s) in enumerate(g.steps)
          kept[j] && !s.rapid_equilibrium || continue
          a, b = seg[s.from], seg[s.to]
          a == b && continue
          lap[a, a] += 1; lap[b, b] += 1; lap[a, b] -= 1; lap[b, a] -= 1
      end
      V * det(lap[2:end, 2:end])
  end
  ```

  `_testhelper_graph` lists the steps in `_flat_edges` order, so `kept[j]` is the mask of
  step `j`'s edge.

- [ ] **Step 3: Write the failing test.** Create `test/test_law.jl`:

  ```julia
  # ABOUTME: Tests of the rate-law engine (src/engine.jl): spanning trees, the tree count,
  # ABOUTME: V×τ, and the plain law against the exact steady state, its form and refusals.

  @testset "_spanning_trees enumerates each spanning tree once" begin
      trees(n, pairs) = Set(Set.(EnzymeRates._spanning_trees(n, pairs)))
      @test trees(3, [(1, 2), (2, 3), (1, 3)]) == Set([Set([1, 2]), Set([1, 3]), Set([2, 3])])
      @test length(EnzymeRates._spanning_trees(4, [(1, 2), (2, 3), (3, 4), (1, 4)])) == 4
      @test EnzymeRates._spanning_trees(2, [(1, 2), (1, 2)]) == [[1], [2]]
      @test isempty(EnzymeRates._spanning_trees(3, [(1, 2)]))
      @test EnzymeRates._spanning_trees(1, Tuple{Int, Int}[]) == [Int[]]
  end

  @testset "_tree_count is the BigInt determinant, saturated at typemax(Int)" begin
      complete(n) = [i == j ? n - 1 : -1 for i in 1:n, j in 1:n]
      @test EnzymeRates._tree_count(complete(5)) == 125
      @test EnzymeRates._tree_count(complete(20)) == typemax(Int)
      @test EnzymeRates._tree_count(zeros(Int, 1, 1)) == 1
      @test EnzymeRates._tree_count([1 -1 0; -1 1 0; 0 0 0]) == 0
      rng = Random.MersenneTwister(31)
      for _ in 1:50
          n = rand(rng, 2:14)
          lap = zeros(Int, n, n)
          for a in 1:n, b in (a + 1):n
              rand(rng) < 0.5 || continue
              k = rand(rng, 1:3)
              lap[a, a] += k; lap[b, b] += k; lap[a, b] -= k; lap[b, a] -= k
          end
          exact = det(BigInt.(lap[2:end, 2:end]))
          @test EnzymeRates._tree_count(lap) ==
                (exact > typemax(Int) ? typemax(Int) : Int(exact))
      end
  end

  @testset "_eq_complexity is segments times spanning trees" begin
      n_plain = 0
      for spec in MECHANISM_TEST_SPECS
          m = _testhelper_lift(spec.mechanism)
          m isa EnzymeRates.Mechanism || continue
          n_plain += 1
          mask = trues(length(m.edges))
          @test EnzymeRates._eq_complexity(m, mask) ===
                Int(_testhelper_laplacian_tree_count(m, mask))
      end
      @test n_plain == 31
      random_bibi = EnzymeRates.Mechanism(@enzyme_mechanism begin
          substrates: S1, S2
          products: P1, P2
          steps: begin
              E + S1 <--> E(S1)
              E + S2 <--> E(S2)
              E(S1) + S2 <--> E(S1, S2)
              E(S2) + S1 <--> E(S1, S2)
              E(S1, S2) <--> E(P1, P2)
              E(P1, P2) <--> E(P1) + P2
              E(P1, P2) <--> E(P2) + P1
              E(P1) <--> E + P1
              E(P2) <--> E + P2
          end
      end)
      @test EnzymeRates._eq_complexity(random_bibi) === 336
      one_segment = EnzymeRates.Mechanism(@enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E + P
          end
      end)
      @test EnzymeRates._eq_complexity(one_segment) === 1
  end

  @testset "the plain law equals the exact steady state on every plain spec" begin
      for (i, spec) in enumerate(MECHANISM_TEST_SPECS)
          m = _testhelper_lift(spec.mechanism)
          m isa EnzymeRates.Mechanism || continue
          @testset "$(spec.name)" begin
              L = EnzymeRates._law(m)
              rng = Random.MersenneTwister(500 + i)
              for _ in 1:3
                  point, concs, params = _testhelper_law_inputs(m,
                                                                EnzymeRates._solve(m).fitted,
                                                                rng)
                  oracle = _testhelper_steady_state_rate(m, point, concs)
                  @test abs(_testhelper_evaluate_law(L, merge(params, concs)) - oracle) <=
                        big(10)^-20 * abs(oracle)
              end
          end
      end
  end

  @testset "a plain law is H times a positive polynomial in lowest terms" begin
      rng = Random.MersenneTwister(7)
      for spec in MECHANISM_TEST_SPECS
          m = _testhelper_lift(spec.mechanism)
          m isa EnzymeRates.Mechanism || continue
          @testset "$(spec.name)" begin
              L = EnzymeRates._law(m)
              @test L.vars.names == [EnzymeRates.metabolites(m);
                                     EnzymeRates._solve(m).fitted; :Keq]
              env = Dict{Symbol, BigFloat}(x => exp(2 * big(rand(rng)) - 1)
                                           for x in L.vars.names)
              H = prod(env[x] for x in EnzymeRates.substrates(m)) -
                  prod(env[x] for x in EnzymeRates.products(m)) / env[:Keq]
              @test isapprox(_testhelper_poly_value(L.H, L.vars, env), H;
                             rtol = big(10)^-60)
              @test !iszero(L.Q) && all(>(0), values(L.Q.terms))
              @test all(>(0), values(L.D.terms))
              @test isempty(EnzymeRates._conc_shift((L.Q, L.D), L.vars))
              @test L.D == sum(L.sigma .* L.trees; init = zero(EnzymeRates._Poly)) * L.shift
          end
      end
  end

  @testset "a conformation's law on its mask is the law of the steps it keeps" begin
      # An inhibitor that binds only the inactive conformation (:OnlyI). The active mask
      # drops its binding and the form E(I), keeping the Michaelis–Menten steps, whose
      # constants both conformations share (:EqualAI); the inactive mask keeps all three
      # steps, the binding under its inactive constant. Both laws number their variables
      # alike: the mechanism's concentrations, the solve's independents, then Keq.
      am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)            :: EqualAI
              E(S) <--> E + P         :: EqualAI
              E + I ⇌ E(I)            :: OnlyI
          end
      end)
      mm = EnzymeRates.Mechanism(@enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E + P
          end
      end)
      competitive = EnzymeRates.Mechanism(@enzyme_mechanism begin
          substrates: S
          products: P
          modifiers: I
          steps: begin
              E + S ⇌ E(S)
              E(S) <--> E + P
              E + I ⇌ E(I)
          end
      end)
      s = EnzymeRates._solve(am)
      cm, mets = EnzymeRates.catalytic(am), EnzymeRates.metabolites(am)
      active = EnzymeRates._law(cm, s, EnzymeRates._active_mask(am); concentrations = mets)
      inactive = EnzymeRates._law(cm, s, EnzymeRates._inactive_mask(am); inactive = true,
                                  concentrations = mets)
      @test active.vars == inactive.vars
      @test active.vars.names == [mets; s.names[s.independents]; :Keq]
      rng = Random.MersenneTwister(71)
      for _ in 1:3
          env = Dict{Symbol, BigFloat}(x => exp(2 * big(rand(rng)) - 1)
                                       for x in active.vars.names)
          env[:K_EI_to_E_I] = env[:K_I_EI_to_E_I]
          @test isapprox(_testhelper_evaluate_law(active, env),
                         _testhelper_evaluate_law(EnzymeRates._law(mm), env);
                         rtol = big(10)^-60)
          @test isapprox(_testhelper_evaluate_law(inactive, env),
                         _testhelper_evaluate_law(EnzymeRates._law(competitive), env);
                         rtol = big(10)^-60)
      end
  end

  @testset "rewriting or regrouping a mechanism changes no plain law" begin
      # Written backwards: every step the other way round, in the reverse order, which
      # canonical form gives back.
      for spec in MECHANISM_TEST_SPECS
          m = _testhelper_lift(spec.mechanism)
          m isa EnzymeRates.Mechanism || continue
          raw = EnzymeRates._raw(m)
          back = EnzymeRates.Mechanism(EnzymeRates.RawMechanism(raw.substrates, raw.products,
              raw.modifiers, reverse([EnzymeRates.RawStep(s.to, s.from, s.released,
                                                          s.consumed, s.rapid_equilibrium)
                                      for s in raw.steps]), reverse(raw.group)))
          @test EnzymeRates._law(back) == EnzymeRates._law(m)
      end
      # A same-model regrouping: the binding square ties B's two dissociation constants,
      # so tying them as one group as well leaves the model and its fitted names alone.
      folded = EnzymeRates.Mechanism(@enzyme_mechanism begin
          substrates: A, B
          products: P
          steps: begin
              (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
              E + B ⇌ E(B)
              E(A) + B ⇌ E(A, B)
              E(A, B) <--> E(P)
              E(P) ⇌ E + P
          end
      end)
      tied = EnzymeRates.Mechanism(@enzyme_mechanism begin
          substrates: A, B
          products: P
          steps: begin
              (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
              (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
              E(A, B) <--> E(P)
              E(P) ⇌ E + P
          end
      end)
      @test folded != tied
      @test EnzymeRates._solve(folded).fitted == EnzymeRates._solve(tied).fitted
      @test EnzymeRates._law(folded) == EnzymeRates._law(tied)
  end

  @testset "plain laws satisfy detailed balance" begin
      for (i, spec) in enumerate(MECHANISM_TEST_SPECS)
          m = _testhelper_lift(spec.mechanism)
          m isa EnzymeRates.Mechanism || continue
          L = EnzymeRates._law(m)
          rng = Random.MersenneTwister(600 + i)
          env = Dict{Symbol, BigFloat}(x => exp(4 * big(rand(rng)) - 2)
                                       for x in L.vars.names)
          p1 = first(EnzymeRates.products(m))
          env[p1] = env[:Keq] * prod(env[x] for x in EnzymeRates.substrates(m)) /
                    prod((env[x] for x in EnzymeRates.products(m) if x != p1);
                         init = big(1))
          off = copy(env)
          off[p1] /= 10
          @test abs(_testhelper_evaluate_law(L, env)) <=
                big(10)^-30 * abs(_testhelper_evaluate_law(L, off))
      end
  end

  @testset "the engine refuses a rapid-equilibrium cycle that turns the reaction over" begin
      # E + S ⇌ E(S) ⇌ E(P) ⇌ E + P closes a rapid-equilibrium cycle that turns S into P
      # beside a steady-state route; the macro rejects it (class 2), and the engine,
      # reading it unchecked, refuses to derive a law for it.
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
      msg = "internal error: a rapid-equilibrium cycle of this mechanism turns the " *
            "reaction over (class 2)"
      @test_throws ErrorException(msg) EnzymeRates._law(beside)
      err = @test_throws ErrorException @enzyme_mechanism begin
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
      @test startswith(err.value.msg,
                       "Mechanism: the rapid-equilibrium steps E + P ⇌ E(P), " *
                       "E + S ⇌ E(S) and E(S) ⇌ E(P) close a cycle")
  end

  @testset "_law refuses a law over the term cap" begin
      # The mechanism of "Rate equation too large error" (test_rate_eq_derivation.jl):
      # 11 forms, 16 steps, V×τ about 29,000.
      m_manual = EnzymeRates.Mechanism(@enzyme_mechanism begin
          substrates: A, B
          products: P, Q
          modifiers: R1
          steps: begin
              E + A <--> E(A)
              E(A) <--> F(P)
              F(P) <--> F + P
              F + B <--> F(B)
              F(B) <--> E(Q)
              E(Q) <--> E + Q
              E + R1 <--> E(R1)
              E(A) + R1 <--> E(A, R1)
              F(P) + R1 <--> F(P, R1)
              F + R1 <--> F(R1)
              F(B) + R1 <--> F(B, R1)
              E(R1) + A <--> E(A, R1)
              E(A, R1) <--> F(P, R1)
              F(P, R1) <--> F(R1) + P
              F(R1) + B <--> F(B, R1)
              F(B, R1) <--> E(R1) + Q
          end
      end)
      @test EnzymeRates._eq_complexity(m_manual) > EnzymeRates.MAX_RATE_EQUATION_TERMS
      err = @test_throws EnzymeRates.RateEquationTooLargeError EnzymeRates._law(m_manual)
      @test sprint(showerror, err.value) ==
            "Rate equation for this mechanism has more than 5000 polynomial terms (limit: " *
            "5000). Equations this large take a very long time to compile and are " *
            "unlikely to be practically useful for parameter fitting."
  end

  @testset "an Int64 overflow in a law is an internal error naming the mechanism" begin
      # Under the term cap no coefficient or exponent leaves Int64. A doctored solve whose
      # dependent reads one constant twice, each time to the power typemax(Int), forces an
      # exponent past it.
      mm = EnzymeRates.Mechanism(@enzyme_mechanism begin
          substrates: S
          products: P
          steps: begin
              E + S <--> E(S)
              E(S) <--> E + P
          end
      end)
      s = EnzymeRates._solve(mm)
      d = only(s.dependents)
      f = first(d.factors).first
      big_power = EnzymeRates._Dependent(d.col, d.keq,
                                         [f => typemax(Int) // 1, f => typemax(Int) // 1],
                                         d.label)
      bad = EnzymeRates._Solve(s.columns, s.names, s.column, s.independents, [big_power],
                               s.fitted)
      err = @test_throws ErrorException EnzymeRates._law(mm, bad, trues(length(mm.edges)))
      @test startswith(err.value.msg, "internal error: the rate law of this mechanism " *
                                      "leaves Int64 arithmetic")
      @test endswith(err.value.msg, sprint(show, mm))
  end
  ```

  The overflow testset builds a `_Solve` by hand: under the cap no real mechanism
  overflows (LAW §1.4), so a doctored dependent is the one way to reach the error. In
  `test/runtests.jl`, replace
  ```julia
      include("test_poly.jl")
      include("test_oracles.jl")
  ```
  with
  ```julia
      include("test_poly.jl")
      include("test_law.jl")
      include("test_oracles.jl")
  ```

- [ ] **Step 4: Watch it fail.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T3_law.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_law.jl")'
  ```
  Expected: the run stops while including the fixtures file:
  ``LoadError: UndefVarError: `_Law` not defined in `EnzymeRates` ``, from
  `_testhelper_evaluate_law`'s signature. Nothing in the tree defines the law yet.

- [ ] **Step 5: Write `src/engine.jl`.**

  ```julia
  # ABOUTME: The rate law of one conformation on an edge mask: rapid-equilibrium segments,
  # ABOUTME: spanning trees, the numerator, the solve and H divided out; V×τ and the term cap.

  """
  The rate law of one conformation, `v = E_total·H·Q/D`, over `vars`: the concentrations,
  then the solve's independent constants (`_Solve.independents` order), then `:Keq`.
  `H = ΠS − ΠP/Keq`, and `Q` has positive coefficients. `D = shift·Σ_g sigma[g]·trees[g]`:
  `sigma[g]` sums the weights of rapid-equilibrium segment `g`'s forms, `trees[g]` is its
  spanning-tree sum, and `shift`, the concentration monomial that brings `Q` and `D` to
  lowest terms, is already in `Q` and `D`. Two laws are equal when their fields are.
  """
  mutable struct _Law
      const vars::_Vars
      const H::_Poly
      const Q::_Poly
      const D::_Poly
      const sigma::Vector{_Poly}
      const trees::Vector{_Poly}
      const shift::_Mono
  end
  Base.:(==)(a::_Law, b::_Law) =
      a.vars == b.vars && a.H == b.H && a.Q == b.Q && a.D == b.D && a.sigma == b.sigma &&
      a.trees == b.trees && a.shift == b.shift
  Base.hash(L::_Law, h::UInt) = hash((L.vars, L.H, L.Q, L.D), hash(:_Law, h))

  """
      _spanning_trees(n, pairs) -> trees

  The spanning trees of the graph on the nodes `1:n` with the edges `pairs`, each as a list
  of positions in `pairs`: each edge is included or left out in turn, and a union-find
  rejects a cycle. A disconnected graph has none; a graph of at most one node has one, with
  no edge.
  """
  function _spanning_trees(n::Int, pairs::Vector{Tuple{Int, Int}})
      n <= 1 && return [Int[]]
      trees = Vector{Int}[]
      function grow(i, chosen, root)
          if length(chosen) == n - 1
              push!(trees, copy(chosen))
              return
          end
          length(pairs) - i + 1 < n - 1 - length(chosen) && return
          top(j) = (while root[j] != j; j = root[j]; end; j)
          a, b = top(pairs[i][1]), top(pairs[i][2])
          if a != b
              r = copy(root)
              r[a] = b
              push!(chosen, i)
              grow(i + 1, chosen, r)
              pop!(chosen)
          end
          grow(i + 1, chosen, root)
      end
      grow(1, Int[], collect(1:n))
      trees
  end

  """The product of the weights `W[(child, parent)]` of the tree `t` (positions in `pairs`)
  oriented toward node `g`."""
  function _toward(t::Vector{Int}, g::Int, pairs::Vector{Tuple{Int, Int}},
                   W::Dict{Tuple{Int, Int}, _Poly})
      p, seen, todo = one(_Poly), Set([g]), [g]
      while !isempty(todo)
          q = pop!(todo)
          for j in t
              u, v = pairs[j]
              c = u == q && !(v in seen) ? v : v == q && !(u in seen) ? u : 0
              c == 0 && continue
              p *= W[(c, q)]
              push!(seen, c)
              push!(todo, c)
          end
      end
      p
  end

  """
      _bare_law(m, mask, concentrations) -> (; vars, cols, num, den, sigma, trees)

  The steady-state law of the edges of `m` that `mask` keeps (at least one), with one
  variable per bare constant: the concentrations `concentrations`, then variable
  `length(concentrations) + j` for slot `cols[j][2]` of edge `cols[j][1]`, in edge order.
  `num/den` is the net consumption of `substrates(m)[1]` per enzyme. Each rapid-equilibrium
  segment (`_segments`) is a node, rooted at its lowest form, so free E (form 1) roots its
  own; a form weighs its root's weight, 1, times the equilibrium ratio of each edge on its
  path from the root. The forms join the tree in `_segments`' order, each by the first kept
  rapid-equilibrium edge in edge (content) order from a form placed before it.
  `den = Σ_g sigma[g]·trees[g]`: `sigma[g]` sums segment `g`'s weights and `trees[g]` its
  spanning trees of the segment graph oriented toward it, each tree edge weighing the
  steady-state flux out of its source. The numerator sums each steady-state edge's net flux
  times the copies of the first substrate it moves, its forms' weights included.

  A rapid-equilibrium edge off the tree must close a cycle with no net uptake of any
  metabolite. Class 2 of the guard rejects every mechanism where one does not, so this is
  an internal error.
  """
  function _bare_law(m::Mechanism, mask::BitVector, concentrations::Vector{Symbol})
      vars = _Vars()
      foreach(x -> _var!(vars, x; conc = true), concentrations)
      cols = Tuple{Int, Int}[]
      at = Dict{Tuple{Int, Int}, Int}()
      for (k, e) in enumerate(m.edges), slot in 1:(e.rapid_equilibrium ? 1 : 2)
          mask[k] || continue
          push!(cols, (k, slot))
          at[(k, slot)] = _var!(vars, Symbol("c", k, "_", slot))
      end
      conc(c::Counts) = _monomial((vars.index[m.table[i].name] => copies(c, i)
                                   for i in eachindex(m.table))...)
      bare(k, slot) = _monomial(at[(k, slot)] => 1)
      # The equilibrium ratio w(to)/w(from) of rapid-equilibrium edge k: [M]/K for a binding
      # (K a dissociation constant), K·Π[consumed]/Π[released] for any other edge.
      ratio(k) = _mmul(_mpow(bare(k, 1), _is_binding(m.edges[k]) ? -1 : 1),
                       _mmul(conc(m.edges[k].consumed), _mpow(conc(m.edges[k].released), -1)))
      s = _segments(m, mask, _flags(m))
      re = [k for k in eachindex(m.edges) if mask[k] && m.edges[k].rapid_equilibrium]
      alpha = fill(_Mono(), length(m.forms))
      placed = falses(length(m.forms))
      for segment in s.segments
          placed[segment[1]] = true
          for v in segment[2:end]
              k = first(k for k in re if (m.edges[k].to == v && placed[m.edges[k].from]) ||
                                         (m.edges[k].from == v && placed[m.edges[k].to]))
              e = m.edges[k]
              alpha[v] = e.to == v ? _mmul(alpha[e.from], ratio(k)) :
                                     _mmul(alpha[e.to], _mpow(ratio(k), -1))
              placed[v] = true
          end
      end
      for k in re
          e = m.edges[k]
          any(q -> vars.conc[q.first],
              _mmul(_mmul(alpha[e.from], ratio(k)), _mpow(alpha[e.to], -1))) && error(
              "internal error: a rapid-equilibrium cycle of this mechanism turns the " *
              "reaction over (class 2)")
      end
      pairs, W = Tuple{Int, Int}[], Dict{Tuple{Int, Int}, _Poly}()
      flux = Tuple{Int, Int, Int, _Mono, _Mono}[]
      x = vars.index[first(substrates(m))]
      for (k, e) in enumerate(m.edges)
          mask[k] && !e.rapid_equilibrium || continue
          f = _mmul(_mmul(bare(k, 1), conc(e.consumed)), alpha[e.from])
          r = _mmul(_mmul(bare(k, 2), conc(e.released)), alpha[e.to])
          a, b = s.segment_of[e.from], s.segment_of[e.to]
          push!(flux, (Int(_mexp(_mmul(f, _mpow(r, -1)), x)), a, b, f, r))
          a == b && continue
          minmax(a, b) in pairs || push!(pairs, minmax(a, b))
          _addto!(get!(_Poly, W, (a, b)), f, 1)
          _addto!(get!(_Poly, W, (b, a)), r, 1)
      end
      G = length(s.segments)
      spanning = _spanning_trees(G, pairs)
      isempty(spanning) && error("internal error: the steady-state steps of this mechanism " *
                                 "do not join its rapid-equilibrium segments")
      trees = [sum((_toward(t, g, pairs, W) for t in spanning); init = zero(_Poly))
               for g in 1:G]
      sigma = [zero(_Poly) for _ in 1:G]
      for (g, segment) in enumerate(s.segments), f in segment
          _addto!(sigma[g], alpha[f], 1)
      end
      num = zero(_Poly)
      for (w, a, b, f, r) in flux
          w == 0 || (num += _scale(trees[a] * f - trees[b] * r, w))
      end
      (; vars, cols, num, den = sum(sigma .* trees; init = zero(_Poly)), sigma, trees)
  end

  """
      _law(m, s, mask; inactive = false, concentrations = metabolites(m)) -> _Law

  The rate law of the edges of `m` that `mask` keeps, in the inactive conformation's
  constants when `inactive`, over the independents of the solve `s`: the bare law
  (`_bare_law`) with each bare constant replaced by its column's value (an independent, or a
  dependent's power product of independents and `Keq`), the numerator divided exactly by `H`,
  then both brought to lowest terms in the concentrations. Every law derived from one solve
  numbers its variables alike: `concentrations`, then `s.independents`, then `Keq`.

  Runs `_assert_derivable` first. Checked `Int64` arithmetic throughout: under the term cap
  no coefficient leaves it, so an `OverflowError` becomes an internal error that names the
  mechanism. Class 5 (`_assert_support`) is left to the callers that use the law.
  """
  function _law(m::Mechanism, s::_Solve, mask::BitVector; inactive::Bool = false,
                concentrations::Vector{Symbol} = metabolites(m))
      _assert_derivable(m, mask)
      try
          b = _bare_law(m, mask, concentrations)
          vars = _Vars()
          foreach(x -> _var!(vars, x; conc = true), concentrations)
          value = Dict{Int, _Mono}(i => _monomial(_var!(vars, s.names[i]) => 1)
                                   for i in s.independents)
          keq = _var!(vars, :Keq)
          for d in s.dependents
              value[d.col] = foldl(_mmul, (_mpow(value[f], e) for (f, e) in d.factors);
                                   init = _monomial(keq => d.keq))
          end
          nc = length(concentrations)
          sub(p) = _subst(v -> v <= nc ? nothing :
                               value[s.column[_Column(b.cols[v - nc]..., inactive)]], p)
          side(xs) = _monomial((vars.index[x] => 1 for x in xs)...)
          H = _term(side(substrates(m))) -
              _term(_mmul(side(products(m)), _monomial(keq => -1)))
          Q, D = _divexact(sub(b.num), H), sub(b.den)
          shift = _conc_shift((Q, D), vars)
          _Law(vars, H, Q * shift, D * shift, map(sub, b.sigma), map(sub, b.trees), shift)
      catch err
          err isa OverflowError || rethrow()
          error("internal error: the rate law of this mechanism leaves Int64 arithmetic, " *
                "which the term cap rules out:\n" * sprint(show, m))
      end
  end

  _law(m::Mechanism) = _law(m, _solve(m), trues(length(m.edges)))

  """
  The spanning trees of the multigraph whose Laplacian is `lap`: the determinant of `lap`
  without its first row and column (Kirchhoff's matrix-tree theorem), by fraction-free
  (Bareiss) elimination in checked `Int`. After an `OverflowError` it takes the `BigInt`
  determinant, saturated at `typemax(Int)`, which is above every cap. A graph of at most one
  node has one tree; a disconnected graph has none.
  """
  function _tree_count(lap::Matrix{Int})
      n = size(lap, 1) - 1
      n <= 0 && return 1
      M = lap[2:end, 2:end]
      sgn, prev = 1, 1
      try
          for k in 1:(n - 1)
              if M[k, k] == 0
                  p = findfirst(i -> M[i, k] != 0, (k + 1):n)
                  p === nothing && return 0
                  p += k
                  for j in 1:n
                      M[k, j], M[p, j] = M[p, j], M[k, j]
                  end
                  sgn = -sgn
              end
              for i in (k + 1):n, j in (k + 1):n
                  M[i, j] = div(Base.checked_sub(Base.checked_mul(M[i, j], M[k, k]),
                                                 Base.checked_mul(M[i, k], M[k, j])), prev)
              end
              prev = M[k, k]
          end
          sgn * M[n, n]
      catch err
          err isa OverflowError || rethrow()
          Int(min(det(BigInt.(lap[2:end, 2:end])), big(typemax(Int))))
      end
  end

  """
  Estimate the term count of the law of the edges of `m` that `mask` keeps, V×τ, from its
  segment graph, without deriving the law: `V` is the number of rapid-equilibrium segments
  (`_segments`), `τ` the spanning trees of the segment graph whose edges are the kept
  steady-state edges between segments, parallel ones counted apart (`_tree_count`). V×τ is
  the number of spanning-tree products in the law's denominator before kinetic grouping, so
  it bounds fit cost. It guards the derivation (`_assert_derivable` aborts when V×τ exceeds
  `MAX_RATE_EQUATION_TERMS`, before any polynomial is built) and backs the
  `eq_complexity_filter` keyword of [`identify_rate_equation`](@ref). Saturates at
  `typemax(Int)`. The method for a mechanism counts every edge of its catalytic graph,
  which for an allosteric mechanism with `:OnlyI` groups holds more than its active
  conformation's law.
  """
  function _eq_complexity(m::Mechanism, mask::BitVector)
      s = _segments(m, mask, _flags(m))
      G = length(s.segments)
      G <= 1 && return G
      lap = zeros(Int, G, G)
      for (k, e) in enumerate(m.edges)
          mask[k] && !e.rapid_equilibrium || continue
          g1, g2 = s.segment_of[e.from], s.segment_of[e.to]
          g1 == g2 && continue
          lap[g1, g1] += 1; lap[g2, g2] += 1; lap[g1, g2] -= 1; lap[g2, g1] -= 1
      end
      τ = _tree_count(lap)
      τ > typemax(Int) ÷ G ? typemax(Int) : G * τ
  end
  _eq_complexity(m::_AnyMech) =
      _eq_complexity(catalytic(m), trues(length(catalytic(m).edges)))

  """
  The error `_assert_derivable` raises for a mechanism whose rate equation would have more
  than `MAX_RATE_EQUATION_TERMS` terms.
  """
  struct RateEquationTooLargeError <: Exception end

  Base.showerror(io::IO, ::RateEquationTooLargeError) = print(io,
      "Rate equation for this mechanism has more than $MAX_RATE_EQUATION_TERMS " *
      "polynomial terms (limit: $MAX_RATE_EQUATION_TERMS). Equations this large " *
      "take a very long time to compile and are unlikely to be practically " *
      "useful for parameter fitting.")

  """
  Abort the derivation of the law of the edges of `m` that `mask` keeps when its term count
  (V×τ, `_eq_complexity`) exceeds `MAX_RATE_EQUATION_TERMS`: a quick count on the segment
  graph that throws a `RateEquationTooLargeError` before any polynomial is built. Every
  derivation entry calls it (`_law`, `_raw_symbolic_rate_polys`).
  """
  function _assert_derivable(m::Mechanism, mask::BitVector)
      _eq_complexity(m, mask) > MAX_RATE_EQUATION_TERMS && throw(RateEquationTooLargeError())
      return nothing
  end
  ```

  Notes for the port:
  - The bare variables are named `c<edge>_<slot>`. A metabolite name holds no `_`, so no
    concentration shares a name, and the chokepoint test sees no `K`, `k`, `V` or `L`
    literal.
  - The RE assertion compares only the concentration part of each cycle's ratio. The
    constants' part holds after the solve is substituted (thermodynamics), not before.
  - The closures capture locals assigned once.

- [ ] **Step 6: Remove the moved code from `src/rate_eq_derivation.jl`.** Delete
  everything from the docstring that begins `Estimate a mechanism's King–Altman
  denominator term count, V×τ, from its` (the `"""` line above it included) through the
  end of `_assert_derivable` and the blank line after it:

  ```julia
  function _assert_derivable(m::Mechanism, mask::BitVector)
      _eq_complexity(m, mask) > MAX_RATE_EQUATION_TERMS && throw(RateEquationTooLargeError())
      return nothing
  end
  ```

  That removes the `_eq_complexity` docstring and both methods, `RateEquationTooLargeError`,
  its `showerror` and `_assert_derivable`'s docstring and body. One blank line then
  separates `_edge_names`' `end` from
  `# ─── Expr generation from POLY ──────────────────────────────`. In `src/EnzymeRates.jl`,
  replace
  ```julia
  include("poly.jl")
  include("admission.jl")
  ```
  with
  ```julia
  include("poly.jl")
  include("engine.jl")
  include("admission.jl")
  ```
  `_raw_symbolic_rate_polys` still calls `_assert_derivable`. Julia resolves the call when
  it runs, so loading `engine.jl` after `rate_eq_derivation.jl` is fine.
  `MAX_RATE_EQUATION_TERMS` stays in `sym_poly_for_rate_eq_derivation.jl` until Task 13.

- [ ] **Step 7: Watch the engine tests pass.** Rerun Step 4's command. Expected: all 11
  testsets pass, 454 tests (as on a9e0699 with this code), in about 3 s, then `exit 0`.

- [ ] **Step 8: Compare the plain laws, V×τ and the gate laws' sizes with the baseline.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && for k in laws complexity dterms; do $S/phase3/jl.sh $S/phase45/logs/P5_T3_compare_$k.log 1800 --project=. --heap-size-hint=2500M $S/phase45/scripts/P5_compare.jl $S/phase45/p5_base $k; done
  ```
  Run the loop in the background and poll. Expected, each log ending `exit 0`:
  - `P5_T3_compare_laws.log`:
    - six `laws <set>` lines;
    - `bibi_init`, `gate` and `inhibitor` read `through _law` only (all plain at f9361cf);
    - `spec` (31 plain, 11 allosteric) and `sample` (67 plain, 53 allosteric) read
      `through _law and _num_den_exprs`, their allosteric members still on the old engine
      until Task 4;
    - `onlyi` reads `through _num_den_exprs`;
    - every worst relative error is at most 1e-20.
  - `P5_T3_compare_complexity.log`: `complexity: 0 of 480 rows differ`. Task 4 changes
    the allosteric rows' first field on purpose; until then every V×τ is the head's.
  - `P5_T3_compare_dterms.log`: `dterms: 0 of 6 rows differ`, now counted on the new law.

  A value past 1e-20 on a plain row, or any `complexity` or `dterms` difference, is a bug
  in Step 5, not a re-pin.

- [ ] **Step 9: Measure the live bytes after the change.** Rerun Step 1 with log
  `$S/phase45/logs/P5_T3_live_after.log`. Expected:
  `_eq_complexity = 336 (Int64); live bytes per call after a full collection: 0.0 KiB`
  (`-0.0` also passes; RES measured 0.00).

- [ ] **Step 10: Update the ROADMAP entry.** In `ROADMAP.md`, entry "**`BigInt` arithmetic
  inflates the GC's live-byte count.**", replace
  ```markdown
    Three functions that compute in `BigInt` raise `Base.gc_live_bytes` by bytes that no live
    object holds: the counter that `_split_gain_test` uses for a plain mechanism
    (`_partition_independent_count`; a `Rational{BigInt}` row reduction), by 18 KiB per call;
    `_eq_complexity` (a `BigInt` determinant, skipped for a single rapid-equilibrium segment),
    by 15 KiB; and class 3's model point (`_model_point`, `src/admission.jl`; an integer null
    space in `Rational{BigInt}`), by 101 KiB on the six-step bi-bi ping-pong at two
    conformations. A 6×6 determinant outside the package does the same in `BigInt` and
  ```
  with
  ```markdown
    Two functions that compute in `BigInt` raise `Base.gc_live_bytes` by bytes that no live
    object holds: the counter that `_split_gain_test` uses for a plain mechanism
    (`_partition_independent_count`; a `Rational{BigInt}` row reduction), by 18 KiB per call,
    and class 3's model point (`_model_point`, `src/admission.jl`; an integer null space in
    `Rational{BigInt}`), by 101 KiB on the six-step bi-bi ping-pong at two conformations. A
    6×6 determinant outside the package does the same in `BigInt` and
  ```
  and replace
  ```markdown
    `Base.gc_live_bytes()` and `Base.gc_num().full_sweep` shows it. The first two figures come
    from a random-order bi-bi of nine steady-state steps, plain and as a two-subunit
    allosteric mechanism with `:OnlyA` on both bindings of A and on the chemical step and
    `:EqualAI` elsewhere (Julia 1.12.5, aarch64, `--heap-size-hint=2500M`).
  ```
  with
  ```markdown
    `Base.gc_live_bytes()` and `Base.gc_num().full_sweep` shows it. The first figure comes
    from a plain random-order bi-bi of nine steady-state steps (Julia 1.12.5, aarch64,
    `--heap-size-hint=2500M`).
  ```
  The phase-4 fix wave does not touch this entry (FW-37 edits three others), so the old
  text matches as quoted.

  The known issue "`eq_complexity_filter` counts an allosteric mechanism's catalytic graph
  once." (FW-37's title; Task 4 retires it) names `_eq_complexity`'s file.
  Wherever `ROADMAP.md` reads
  ```markdown
  `_eq_complexity` (`src/rate_eq_derivation.jl`)
  ```
  write
  ```markdown
  `_eq_complexity` (`src/engine.jl`)
  ```
  so that entry stays true until Task 4 retires it.

- [ ] **Step 11: Run the full suite** (log `$S/phase45/logs/P5_T3_full.log`, background,
  polled). Expected:
  - green and pristine;
  - the pass count is Task 2's plus `test_law.jl`'s 454;
  - TR "_eq_complexity (V×τ term-count estimate)" and "Rate equation too large error"
    pass unchanged;
  - test_identify_rate_equation.jl's filter and guard testsets pass, now reading an
    `Int`;
  - wall time within a few seconds of Task 2's (the budget is +3 s).

  A maxrss rise over Task 0's figure is explained in the commit body.

- [ ] **Step 12: Commit.**

  ```bash
  S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && git add src/engine.jl src/EnzymeRates.jl src/rate_eq_derivation.jl test/test_law.jl test/runtests.jl test/mechanism_definitions_for_test_enzyme_derivation.jl ROADMAP.md && git commit -F - <<'EOF'
  Derive the plain rate law with the exact engine, not yet wired

  engine.jl derives one conformation's law on an edge mask: the
  rapid-equilibrium forest from _segments, spanning trees of the segment
  graph, the edge-sum numerator, the solve substituted and the Haldane
  factor divided out exactly, all in checked Int64. On every plain
  mechanism of the phase-5 baseline it equals the old engine's law to
  1e-20, and on the plain specs the exact steady state. Each conformation of
  an :OnlyI mechanism, on its own mask, gives the law of the steps it keeps.
  The engine refuses a rapid-equilibrium cycle that turns the reaction
  over, which class 2 rejects first, and turns an overflow into an internal
  error that names the mechanism. _eq_complexity moves here with the cap
  and counts in checked Int by Bareiss, the BigInt determinant only after
  an overflow. Its counts equal the old ones on the 480 baseline
  mechanisms, and it leaves <x_after> KiB of live bytes per call where it
  left <x_before>, so the ROADMAP's live-byte entry drops it. Nothing reads
  the new law yet.

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
  EOF
  ```
  Replace `<x_before>` and `<x_after>` with Steps 1 and 9's figures.

---

### Task 4: MWC law on both masks; summed cap and filter (P5-8, P5-26)

`_law(am, s)` derives each conformation's law on its own mask, over the same variables, and
combines them as RES §1.5 writes. The cap and identify's filter both read
`_eq_complexity(am)`, the active count plus the inactive one. The old MWC path refuses the
same laws, because `_num_den_exprs(am)` calls `_assert_derivable(am)` first.

**Files:**
- Modify: `src/engine.jl` (`_AlloLaw` with `==` and `hash`, `_flip_only_law`,
  `_flip_only`, `_cross_weights`, `_site_weight`, `_law(am, s)`, `_law(am)`, `_law(m, s)`,
  `_eq_complexity(m)`, `_eq_complexity(am)`, `_assert_derivable(am)`)
- Modify: `src/rate_eq_derivation.jl` (`_num_den_exprs(am)` calls `_assert_derivable(am)`
  first)
- Modify: `src/identify_rate_equation.jl` (the `eq_complexity_filter` keyword's docstring)
- Modify: `ROADMAP.md` (retire "`eq_complexity_filter` counts an allosteric
  mechanism's catalytic graph once")
- Modify: the fixtures file (`_testhelper_evaluate_law(::_AlloLaw, env)`;
  `_testhelper_detailed_balance` rebased on `_law`)
- Test: `test/test_law.jl`, `test/test_oracles.jl` (Q1), `test/test_rate_eq_derivation.jl`
  (pin 12), `test/test_mechanism_enumeration.jl` (pin 13),
  `test/test_identify_rate_equation.jl`
- Create (scratch): `$S/phase45/probes/P5_T4_cap_fixture.jl`,
  `$S/phase45/scripts/P5_r55_recount.jl`
- Logs: `$S/phase45/logs/P5_T4_{tid_red,law_red,law,cap_red,cap,tr,tme,tid,oracles,compare,recount,full}.log`

**Interfaces:**
- Consumes: Task 3's `_law(m::Mechanism, s::_Solve, mask::BitVector; inactive,
  concentrations)::_Law`, `_Law` (`vars`, `H`, `Q`, `D`, `sigma`, `trees`, `shift`),
  `_eq_complexity(m::Mechanism, mask::BitVector)::Int`,
  `_assert_derivable(m::Mechanism, mask)`, `RateEquationTooLargeError`; Task 2's `_Vars`,
  `_var!`, `_Poly`, `_term`, `_monomial`, `_mmul`, `_mpow`, `_mexp`, `_content`, `+`, `*`;
  Task 1's `_solve`, `_Solve` (`independents`, `fitted`); `_active_mask`, `_inactive_mask`
  (on `_exclusion_components`), `catalytic`, `catalytic_multiplicity`,
  `regulatory_sites`, `ligands`, `allo_states`, `multiplicity`, `metabolites(am)`,
  `_reg_K`. Tests: Task 3's `_testhelper_poly_value(p, vars, env)`,
  `_testhelper_evaluate_law(::_Law, env)`, `_testhelper_laplacian_tree_count`; phase 4's
  `_testhelper_law_inputs`, `_testhelper_steady_state_rate`, `_testhelper_limit_law`,
  `_testhelper_protomer_rate`, `_testhelper_onlyi_fixtures`, `_testhelper_declared_dimension`,
  `_testhelper_single_subunit`, `_testhelper_lift`, `_testhelper_graph`.
- Produces:
  ```julia
  mutable struct _AlloLaw
      const vars::_Vars        # active.vars, then the regulator constants and :L
      const H::_Poly
      const active::_Law
      const inactive::_Law     # flip-only: Q = 0, D = 1
      const fA::_Poly          # free E's weights over their common factor
      const fI::_Poly
      const sites::Vector{Tuple{_Poly, _Poly, Int}}   # (W_A, W_I, multiplicity)
      const n::Int
  end
  Base.:(==)(a::_AlloLaw, b::_AlloLaw); Base.hash(A::_AlloLaw, h::UInt)   # over the fields
  _flip_only_law(LA::_Law)::_Law
  _flip_only(L::_Law)::Bool
  _cross_weights(fA::_Poly, fI::_Poly)::Tuple{_Poly, _Poly}
  _site_weight(am::AllostericMechanism, vars::_Vars, site::RegulatorySite,
               inactive::Bool)::_Poly
  _law(m::Mechanism, s::_Solve)::_Law             # every edge
  _law(am::AllostericMechanism, s::_Solve)::_AlloLaw
  _law(am::AllostericMechanism)::_AlloLaw         # _law(am, _solve(am))
  _eq_complexity(m::Mechanism)::Int
  _eq_complexity(am::AllostericMechanism)::Int    # active + (any(imask) ? inactive : 0)
  _assert_derivable(am::AllostericMechanism)::Nothing
  # tests (fixtures file)
  _testhelper_evaluate_law(A::EnzymeRates._AlloLaw, env::AbstractDict{Symbol})::BigFloat
  _testhelper_detailed_balance(m, rng)::Vector{BigFloat}   # on _law
  ```

**Gate:**
- `test_law.jl` is green: every allosteric spec (11) equals the one-subunit oracle to 1e-20
  at 3 points; the 1–3-protomer fixtures equal `_testhelper_protomer_rate` to 1e-20; the 13
  `:OnlyI` laws equal `_testhelper_limit_law` to 1e-20 (`essential_activator` and
  `copy_group_bridge` included), and `D_A` holds none of their `:OnlyI` constants; detailed
  balance holds per conformation and combined; `length(s.fitted)` equals
  `_testhelper_declared_dimension`.
- The cap fixture counts 4,320 + 4,320 = 8,640; `_law(am)` and `rate_equation_string(am)`
  throw `RateEquationTooLargeError` with today's text, the second at once, where the head's
  path does not refuse it (Step 8).
- The filter fixture counts 336 + 16 = 352: skipped at 337; at 700 it passes the filter
  and meets the parameter cap. The test is watched failing first (Step 2).
- Pins 12 and 13 move after their value checks.
- `P5_compare.jl … laws` passes on every allosteric baseline mechanism to 1e-20.
- `P5_r55_recount.jl` reports the candidates newly over 337 per round; about 10 min.
- The full suite is green and pristine; +3 s.

- [ ] **Step 1: Write the failing filter test**

Append to `test/test_identify_rate_equation.jl`, after the testset "derivation guard
(filter > MAX_RATE_EQUATION_TERMS) → FitFailure, not crash":

```julia
@testset "eq_complexity_filter counts both conformations" begin
    # The 336-term random-order bi-bi at one catalytic subunit with its chemistry and its
    # product bindings at free E :OnlyA, a K-type the tagging move emits: the inactive
    # conformation keeps the substrate square E, E(S1), E(S1, S2), E(S2), four segments and
    # four spanning trees, so the law counts 336 + 16 = 352. The default filter skips it,
    # where the active count alone would let it through. A filter of 700 lets it through
    # to the parameter cap, which skips its 16 fitted constants without a fit.
    rxn = @enzyme_reaction begin
        substrates: S1[C2H4], S2[C2H2]
        products:   P1[C2H2], P2[C2H4]
    end
    data = DataFrame(
        S1 = [1.0, 2.0, 1.0, 2.0], S2 = [1.0, 1.0, 2.0, 2.0],
        P1 = [0.1, 0.2, 0.1, 0.2], P2 = [0.1, 0.1, 0.2, 0.2],
        Rate = [0.5, 0.8, 0.9, 1.1], group = [1, 1, 2, 2])
    prob = IdentifyRateEquationProblem(rxn, data; Keq=2.0)
    ktype = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S1, S2
        products: P1, P2
        catalytic_multiplicity: 1
        steps: begin
            E + P1 <--> E(P1) :: OnlyA
            E + P2 <--> E(P2) :: OnlyA
            E + S1 <--> E(S1) :: EqualAI
            E + S2 <--> E(S2) :: EqualAI
            E(P1) + P2 <--> E(P1, P2) :: EqualAI
            E(P2) + P1 <--> E(P1, P2) :: EqualAI
            E(S1) + S2 <--> E(S1, S2) :: EqualAI
            E(S1, S2) <--> E(P1, P2) :: OnlyA
            E(S2) + S1 <--> E(S1, S2) :: EqualAI
        end
    end)
    cm = ktype.catalytic
    @test _testhelper_laplacian_tree_count(cm, EnzymeRates._active_mask(ktype)) == 336
    @test _testhelper_laplacian_tree_count(cm, EnzymeRates._inactive_mask(ktype)) == 16
    @test EnzymeRates._eq_complexity(ktype) == 352
    @test EnzymeRates._independent_param_count(ktype) == 16
    batch = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[ktype]
    _, _, n_param_skip, n_cx_skip = EnzymeRates._process_batch(
        batch, prob; optimizer=CMAEvolutionStrategyOpt(),
        max_param_count=15, eq_complexity_filter=337, n_restarts=1, maxtime=1.0)
    @test (n_cx_skip, n_param_skip) == (1, 0)
    _, _, n_param_wide, n_cx_wide = EnzymeRates._process_batch(
        batch, prob; optimizer=CMAEvolutionStrategyOpt(),
        max_param_count=15, eq_complexity_filter=700, n_restarts=1, maxtime=1.0)
    @test (n_cx_wide, n_param_wide) == (0, 1)
end
```

The parameter cap of 15 stops the mechanism after the filter in both calls, so the test
fits nothing (`_process_batch` checks the filter, then the fitted count).

- [ ] **Step 2: Run the filter test to verify it fails**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T4_tid_red.log 3000 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); @testset "red" begin include("test/test_identify_rate_equation.jl") end'`, polled with the Runner's bounded loop (Global Constraints).
Expected: every other testset passes. "eq_complexity_filter counts both conformations"
fails twice: `_eq_complexity(ktype) == 352` evaluates `336 == 352` (Task 3 counts the whole
catalytic graph), and the first batch gives `(0, 1) == (1, 0)` (the filter lets 336
through). The two Laplacian counts, the fitted count and the second batch pass. The
outer `red` testset then reports both failures, and the log ends `exit 1`.

- [ ] **Step 3: Write the test helpers**

In `test/mechanism_definitions_for_test_enzyme_derivation.jl`, directly after Task 3's
`_testhelper_evaluate_law(L::EnzymeRates._Law, env::AbstractDict{Symbol})`, add:

```julia
"""
The MWC law `A` per enzyme at the values `env` (by name: the concentrations, the
independents, `Keq`, the regulator constants and `L`), in 256-bit BigFloat:
`H·(fI^n·QA·DA^(n−1)·WA + L·fA^n·QI·DI^(n−1)·WI) / (fI^n·DA^n·WA + L·fA^n·DI^n·WI)`, each
site's binding polynomial raised to its multiplicity in `WA` and `WI`.
"""
function _testhelper_evaluate_law(A::EnzymeRates._AlloLaw, env::AbstractDict{Symbol})
    setprecision(BigFloat, 256) do
        value(p) = _testhelper_poly_value(p, A.vars, env)
        n, L = A.n, BigFloat(env[:L])
        WA = prod((value(w)^k for (w, _, k) in A.sites); init = big(1))
        WI = prod((value(w)^k for (_, w, k) in A.sites); init = big(1))
        fA, fI = value(A.fA), value(A.fI)
        QA, DA = value(A.active.Q), value(A.active.D)
        QI, DI = value(A.inactive.Q), value(A.inactive.D)
        value(A.H) * (fI^n * QA * DA^(n - 1) * WA + L * fA^n * QI * DI^(n - 1) * WI) /
            (fI^n * DA^n * WA + L * fA^n * DI^n * WI)
    end
end
```

Replace the whole of `_testhelper_detailed_balance` (docstring and function) by:

```julia
"""
    _testhelper_detailed_balance(m, rng) -> Vector{BigFloat}

Detailed-balance residuals of `m`'s derived law (`EnzymeRates._law`) at a random parameter
point, in 256-bit BigFloat: the law's value at equilibrium concentrations (the first product
scaled so that Π[P]/Π[S] = Keq) over its value with that product ten times lower, then, for
an allosteric mechanism, the same ratio of `H·Q` for each conformation whose `Q` is not
zero. Each is zero in exact arithmetic when the law satisfies every Haldane and Wegscheider
relation. `m`, compiled or concrete, is read as the concrete mechanism, without the guard,
so a law the guard rejects is checked too.
"""
function _testhelper_detailed_balance(m, rng)
    setprecision(BigFloat, 256) do
        mech = _testhelper_lift(m)
        s = EnzymeRates._solve(mech)
        L = EnzymeRates._law(mech, s)
        params = Dict{Symbol, BigFloat}(k => exp(4 * big(rand(rng)) - 2) for k in s.fitted)
        params[:Keq] = exp(4 * big(rand(rng)) - 2)
        g = _testhelper_graph(mech)
        subs, prods = g.substrates, g.products
        concs = Dict{Symbol, BigFloat}(x => exp(2 * big(rand(rng)) - 1)
                                       for x in EnzymeRates.metabolites(mech))
        p1 = first(prods)
        concs[p1] = params[:Keq] * prod(concs[x] for x in subs) /
                    prod((concs[x] for x in prods if x != p1); init = big(1))
        off = copy(concs)
        off[p1] /= 10
        at, below = merge(params, concs), merge(params, off)
        ratios = BigFloat[abs(_testhelper_evaluate_law(L, at)) /
                          abs(_testhelper_evaluate_law(L, below))]
        L isa EnzymeRates._AlloLaw || return ratios
        for C in (L.active, L.inactive)
            iszero(C.Q) && continue
            h(env) = _testhelper_poly_value(C.H, C.vars, env)
            q(env) = _testhelper_poly_value(C.Q, C.vars, env)
            push!(ratios, abs(h(at) * q(at)) / abs(h(below) * q(below)))
        end
        ratios
    end
end
```

The draws keep the head's order (fitted constants in `fitted_params` order, `Keq`, then the
concentrations), so every caller sees the head's points.

- [ ] **Step 4: Write the failing law tests**

Append to `test/test_law.jl`:

```julia
# ── The MWC law ─────────────────────────────────────────────────────────────

@testset "the MWC law equals the one-subunit oracle on every allosteric spec" begin
    # At one catalytic subunit the MWC law is the exact steady state of one protomer that
    # flips at free E: 3 points per spec, 1e-20 relative in 256-bit BigFloat. Each
    # conformation's Q has positive coefficients, the law vanishes at equilibrium in each
    # conformation and combined, and the fitted names count the declared dimension. A
    # second derivation gives an equal law with an equal hash.
    specs = [s for s in MECHANISM_TEST_SPECS
             if s.mechanism isa EnzymeRates.AllostericEnzymeMechanism]
    @test length(specs) == 11
    for (i, spec) in enumerate(specs)
        @testset "$(spec.name)" begin
            am = _testhelper_lift(_testhelper_single_subunit(spec.mechanism))
            s = EnzymeRates._solve(am)
            A = EnzymeRates._law(am, s)
            rng = Random.MersenneTwister(i)
            for _ in 1:3
                point, concs, params = _testhelper_law_inputs(am, s.fitted, rng)
                oracle = _testhelper_steady_state_rate(am, point, concs)
                law = _testhelper_evaluate_law(A, merge(params, concs))
                @test abs(law - oracle) <= big(10)^-20 * abs(oracle)
            end
            @test all(>(0), values(A.active.Q.terms))
            @test all(>(0), values(A.inactive.Q.terms))
            @test all(<(big(10)^-30),
                      _testhelper_detailed_balance(am, Random.MersenneTwister(100 + i)))
            @test length(s.fitted) == _testhelper_declared_dimension(am)
            again = EnzymeRates._law(am, s)
            @test again == A && hash(again) == hash(A)
        end
    end
end

@testset "the MWC law equals the n-protomer oracle at 1–3 protomers" begin
    # The bi-uni of the :NonequalAI ground-truth gate: A binds at steady state, B and P at
    # rapid equilibrium, and the chemistry differs between the conformations, so the
    # inactive numerator is live. The oracle solves the concerted network of n protomers,
    # each moving through its conformation's forms, all flipping together at free E.
    mechs = (
        @allosteric_mechanism(begin
            substrates: A, B
            products: P
            catalytic_multiplicity: 1
            steps: begin
                E + A <--> E(A)        :: EqualAI
                E(A) + B ⇌ E(A, B)     :: EqualAI
                E(A, B) <--> E(P)      :: NonequalAI
                E + P ⇌ E(P)           :: EqualAI
            end
        end),
        @allosteric_mechanism(begin
            substrates: A, B
            products: P
            catalytic_multiplicity: 2
            steps: begin
                E + A <--> E(A)        :: EqualAI
                E(A) + B ⇌ E(A, B)     :: EqualAI
                E(A, B) <--> E(P)      :: NonequalAI
                E + P ⇌ E(P)           :: EqualAI
            end
        end),
        @allosteric_mechanism(begin
            substrates: A, B
            products: P
            catalytic_multiplicity: 3
            steps: begin
                E + A <--> E(A)        :: EqualAI
                E(A) + B ⇌ E(A, B)     :: EqualAI
                E(A, B) <--> E(P)      :: NonequalAI
                E + P ⇌ E(P)           :: EqualAI
            end
        end))
    for (n, em) in enumerate(mechs)
        am = _testhelper_lift(em)
        s = EnzymeRates._solve(am)
        A = EnzymeRates._law(am, s)
        @test A.n == n
        rng = Random.MersenneTwister(30 + n)
        for _ in 1:2
            point, concs, params = _testhelper_law_inputs(am, s.fitted, rng)
            oracle = setprecision(() -> _testhelper_protomer_rate(am, point, concs; n),
                                  BigFloat, 256)
            law = _testhelper_evaluate_law(A, merge(params, concs))
            @test abs(law - oracle) <= big(10)^-20 * abs(oracle)
        end
    end
end

@testset "an :OnlyI law is the exact limit of its vanishing constants" begin
    # The thirteen mechanisms whose :OnlyI groups bind modifiers alone, read without the
    # guard: each law equals the rapid-equilibrium limit of the oracle's own deletion
    # network (`_testhelper_limit_law`) at 3 points, and the active conformation's
    # denominator holds no constant of an :OnlyI group, since `_active_mask` drops them.
    for (i, (label, am)) in enumerate(pairs(_testhelper_onlyi_fixtures()))
        @testset "$label" begin
            sub = _testhelper_single_subunit(am)
            s = EnzymeRates._solve(sub)
            A = EnzymeRates._law(sub, s)
            rng = Random.MersenneTwister(400 + i)
            for _ in 1:3
                point, concs, params = _testhelper_law_inputs(sub, s.fitted, rng)
                exact = setprecision(BigFloat, 256) do
                    _testhelper_limit_law(sub, point, concs, Set{Symbol}(); rng)
                end
                law = _testhelper_evaluate_law(A, merge(params, concs))
                @test abs(law - exact) <= big(10)^-20 * abs(exact)
            end
            cm = EnzymeRates.catalytic(sub)
            onlyi = Set(x for (k, ns) in enumerate(EnzymeRates._edge_names(sub, :I))
                        if EnzymeRates.cat_allo_state(sub, cm.group[k]) === :OnlyI
                        for x in ns)
            held = Set(A.vars.names[v] for mo in keys(A.active.D.terms) for (v, _) in mo)
            @test !isempty(onlyi)
            @test isempty(intersect(onlyi, held))
        end
    end
end

@testset "a flip-only inactive conformation is Q = 0, D = 1" begin
    # Both bindings :OnlyA: the inactive conformation's free E binds nothing, so that
    # conformation is free E alone. It turns nothing over and weighs 1, so the law is
    # H·QA/(DA + L), and it adds nothing to the term count.
    inert = @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: EqualAI
            E + P ⇌ E(P)      :: OnlyA
        end
    end
    am = _testhelper_lift(inert)
    @test !any(EnzymeRates._inactive_mask(am))
    s = EnzymeRates._solve(am)
    A = EnzymeRates._law(am, s)
    @test EnzymeRates._flip_only(A.inactive) && !EnzymeRates._flip_only(A.active)
    @test iszero(A.inactive.Q) && A.inactive.D == one(EnzymeRates._Poly)
    @test A.inactive.trees == [one(EnzymeRates._Poly)] && isempty(A.inactive.shift)
    @test A.fA == one(EnzymeRates._Poly) && A.fI == one(EnzymeRates._Poly)
    @test EnzymeRates._eq_complexity(am) == 1
    rng = Random.MersenneTwister(77)
    for _ in 1:3
        point, concs, params = _testhelper_law_inputs(am, s.fitted, rng)
        oracle = _testhelper_steady_state_rate(am, point, concs)
        law = _testhelper_evaluate_law(A, merge(params, concs))
        @test abs(law - oracle) <= big(10)^-20 * abs(oracle)
    end
end

@testset "_cross_weights divides out the common integer and monomial content" begin
    vars = EnzymeRates._Vars()
    B = EnzymeRates._var!(vars, :B; conc = true)
    k = EnzymeRates._var!(vars, :k_ES_to_EP)
    q = EnzymeRates._var!(vars, :K_EQ_to_E_Q)
    t(c, ps...) = EnzymeRates._term(EnzymeRates._monomial(ps...), c)
    # Two monomials, 2·k·B and 4·B², share the content 2·B.
    @test EnzymeRates._cross_weights(t(2, k => 1, B => 1), t(4, B => 2)) ==
          (t(1, k => 1), t(2, B => 1))
    # Equal weights both become 1, whatever their polynomial part.
    w = t(1, k => 1) + t(1, B => 1)
    @test EnzymeRates._cross_weights(w, w) ==
          (one(EnzymeRates._Poly), one(EnzymeRates._Poly))
    # A two-term weight keeps its primitive part; only the shared B comes out.
    @test EnzymeRates._cross_weights(t(1, q => 1, B => 1) * w, t(1, B => 1)) ==
          (t(1, q => 1) * w, one(EnzymeRates._Poly))
end

@testset "the cap counts both conformations" begin
    # A random-order bi-bi with the abortive complexes E(A, Q) and E(B, P), every step
    # steady state: 9 forms and 480 spanning trees, 4,320 products. With its chemistry
    # :NonequalAI at one catalytic subunit both conformations keep every step, so the law
    # counts 4,320 + 4,320 = 8,640: over the cap of 5,000, though each conformation is
    # under it. The engine refuses it before any polynomial product, with
    # RateEquationTooLargeError's message.
    em = @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)          :: EqualAI
            E + B <--> E(B)          :: EqualAI
            E(A) + B <--> E(A, B)    :: EqualAI
            E(B) + A <--> E(A, B)    :: EqualAI
            E(A, B) <--> E(P, Q)     :: NonequalAI
            E(P, Q) <--> E(P) + Q    :: EqualAI
            E(P, Q) <--> E(Q) + P    :: EqualAI
            E(P) <--> E + P          :: EqualAI
            E(Q) <--> E + Q          :: EqualAI
            E(A) + Q <--> E(A, Q)    :: EqualAI
            E(Q) + A <--> E(A, Q)    :: EqualAI
            E(B) + P <--> E(B, P)    :: EqualAI
            E(P) + B <--> E(B, P)    :: EqualAI
        end
    end
    am = _testhelper_lift(em)
    cm = am.catalytic
    amask, imask = EnzymeRates._active_mask(am), EnzymeRates._inactive_mask(am)
    @test _testhelper_laplacian_tree_count(cm, amask) == 4320
    @test _testhelper_laplacian_tree_count(cm, imask) == 4320
    @test EnzymeRates._eq_complexity(cm, amask) <= EnzymeRates.MAX_RATE_EQUATION_TERMS
    @test EnzymeRates._eq_complexity(am) == 8640
    @test_throws EnzymeRates.RateEquationTooLargeError EnzymeRates._law(am)
    @test_throws "has more than 5000 polynomial terms" EnzymeRates._law(am)
end
```

- [ ] **Step 5: Run the law tests to verify they fail**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T4_law_red.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_law.jl")'`
Expected: the run stops while including the fixtures file, at Step 3's method signature:
``LoadError: UndefVarError: `_AlloLaw` not defined in `EnzymeRates` ``, as Task 3's first
run stopped at `_Law`. Nothing defines the MWC law yet, so every test of this step fails
to load; Step 2 already watched the filter test fail on its assertions.

- [ ] **Step 6: Write the MWC law, the summed count and the cap**

In `src/engine.jl` make three changes.

First, replace Task 3's `_law(m::Mechanism) = _law(m, _solve(m), trues(length(m.edges)))`
with

```julia
"""The law of the edges of `m` over the independents of `s`, every edge kept."""
_law(m::Mechanism, s::_Solve) = _law(m, s, trues(length(m.edges)))
_law(m::Mechanism) = _law(m, _solve(m))
```

and insert after it:

```julia
# ─── The MWC law ─────────────────────────────────────────────────────────

"""
An MWC rate law over `vars`: `active.vars` (the concentrations of `metabolites(am)`, the
independents of both conformations and `:Keq`), then the regulator constants (in
`_Solve.fitted` order) and `:L`. With `QA, DA` from `active` and `QI, DI` from `inactive`
(0 and 1 when the inactive conformation is free E alone, `_flip_only_law`), and
`WA = Π_s sites[s][1]^sites[s][3]`, `WI = Π_s sites[s][2]^sites[s][3]`:

    v = E_total·H·(fI^n·QA·DA^(n−1)·WA + L·fA^n·QI·DI^(n−1)·WI) /
                  (fI^n·DA^n·WA + L·fA^n·DI^n·WI)

`fA` and `fI` are free E's weights in the two laws over their common factor
(`_cross_weights`); `n` is the catalytic multiplicity.
"""
mutable struct _AlloLaw
    const vars::_Vars
    const H::_Poly
    const active::_Law
    const inactive::_Law
    const fA::_Poly
    const fI::_Poly
    const sites::Vector{Tuple{_Poly, _Poly, Int}}
    const n::Int
end
Base.:(==)(a::_AlloLaw, b::_AlloLaw) =
    a.vars == b.vars && a.H == b.H && a.active == b.active && a.inactive == b.inactive &&
    a.fA == b.fA && a.fI == b.fI && a.sites == b.sites && a.n == b.n
Base.hash(A::_AlloLaw, h::UInt) =
    hash((A.vars, A.active, A.inactive, A.n), hash(:_AlloLaw, h))

"""The law of an inactive conformation that is free E alone, over `LA`'s variables: it
turns nothing over (`Q = 0`), weighs 1 (`D = 1`) and has no segment."""
_flip_only_law(LA::_Law) =
    _Law(LA.vars, LA.H, zero(_Poly), one(_Poly), _Poly[], [one(_Poly)], _Mono())

"""Whether `L` is the law of a conformation that is free E alone (`_flip_only_law`)."""
_flip_only(L::_Law) = isempty(L.sigma)

"""
    _cross_weights(fA, fI) -> (fA′, fI′)

The free-E weights `fA` and `fI` of the two conformations over their common factor: each
divided by the greatest common divisor of their integer contents and by the monomial of
the lower exponent of each variable in their monomial contents (`_content`). When their
primitive parts are equal those cancel too, leaving monomials, 1 for equal weights.
"""
function _cross_weights(fA::_Poly, fI::_Poly)
    (cA, uA, rA), (cI, uI, rI) = _content(fA), _content(fI)
    same = rA == rI
    g = gcd(cA, cI)
    u = _monomial((v => min(_mexp(uA, v), _mexp(uI, v))
                   for v in unique!([first.(uA); first.(uI)]))...)
    weight(c, mu, r) = (p = _term(_mmul(mu, _mpow(u, -1)), div(c, g)); same ? p : p * r)
    weight(cA, uA, rA), weight(cI, uI, rI)
end

"""
The binding polynomial of regulatory site `site` of `am` in the active (`inactive =
false`) or inactive conformation, over `vars`: 1 plus each ligand's concentration over its
dissociation constant (`_reg_K`), an `:OnlyI` ligand left out of the active conformation
and an `:OnlyA` one out of the inactive.
"""
function _site_weight(am::AllostericMechanism, vars::_Vars, site::RegulatorySite,
                      inactive::Bool)
    W = one(_Poly)
    for (lig, tag) in zip(ligands(site), allo_states(site))
        K = _reg_K(am, lig, tag, inactive)
        K === nothing && continue
        W = W + _term(_monomial(vars.index[lig] => 1, vars.index[K] => -1))
    end
    W
end

"""
    _law(am, s = _solve(am)) -> _AlloLaw

The MWC law of `am` over the independents of `s`: the active conformation's law on
`_active_mask(am)`, the inactive one's on `_inactive_mask(am)` (`_flip_only_law` when that
mask is empty), both over the concentrations of `metabolites(am)`, so their polynomials
multiply directly; free E's weight in each over their common factor; and each regulatory
site's binding polynomials. Free E roots segment 1, so its weight is that segment's tree
sum times the law's shift. A law over the summed cap is refused before any polynomial
product (`_assert_derivable`).
"""
function _law(am::AllostericMechanism, s::_Solve)
    _assert_derivable(am)
    m, imask = catalytic(am), _inactive_mask(am)
    concentrations = metabolites(am)
    LA = _law(m, s, _active_mask(am); concentrations)
    LI = any(imask) ? _law(m, s, imask; inactive = true, concentrations) :
         _flip_only_law(LA)
    vars = deepcopy(LA.vars)
    for x in view(s.fitted, (length(s.independents) + 1):(length(s.fitted) - 1))
        _var!(vars, x)
    end
    _var!(vars, :L)
    fA, fI = _cross_weights(LA.trees[1] * LA.shift, LI.trees[1] * LI.shift)
    sites = Tuple{_Poly, _Poly, Int}[(_site_weight(am, vars, site, false),
                                      _site_weight(am, vars, site, true),
                                      multiplicity(site))
                                     for site in regulatory_sites(am)]
    _AlloLaw(vars, LA.H, LA, LI, fA, fI, sites, catalytic_multiplicity(am))
end
_law(am::AllostericMechanism) = _law(am, _solve(am))
```

Second, delete Task 3's method `_eq_complexity(m::_AnyMech)`, which counts the whole
catalytic graph,

```julia
_eq_complexity(m::_AnyMech) =
    _eq_complexity(catalytic(m), trues(length(catalytic(m).edges)))
```

and the last sentence of the docstring above `_eq_complexity(m::Mechanism,
mask::BitVector)`, "The method for a mechanism counts every edge of its catalytic graph,
which for an allosteric mechanism with `:OnlyI` groups holds more than its active
conformation's law.", and put in the deleted method's place:

```julia
"""
    _eq_complexity(m) -> Int

The size of `m`'s rate law in spanning-tree products, V×τ (`_eq_complexity(m, mask)`):
over every edge of a plain mechanism; for an allosteric one, the active law's count plus
the inactive law's, 0 when the inactive conformation is free E alone, saturated at
`typemax(Int)`. The cap (`_assert_derivable`) and identify's `eq_complexity_filter` read
it.
"""
_eq_complexity(m::Mechanism) = _eq_complexity(m, trues(length(m.edges)))
function _eq_complexity(am::AllostericMechanism)
    m, imask = catalytic(am), _inactive_mask(am)
    a = _eq_complexity(m, _active_mask(am))
    i = any(imask) ? _eq_complexity(m, imask) : 0
    a > typemax(Int) - i ? typemax(Int) : a + i
end
```

Third, append after Task 3's `_assert_derivable(m::Mechanism, mask::BitVector)`:

```julia
"""Abort the derivation of `am`'s MWC law when its two conformations together exceed
`MAX_RATE_EQUATION_TERMS` spanning-tree products (`_eq_complexity(am)`): a numeric check
that throws a `RateEquationTooLargeError` before any polynomial product."""
function _assert_derivable(am::AllostericMechanism)
    _eq_complexity(am) > MAX_RATE_EQUATION_TERMS && throw(RateEquationTooLargeError())
    return nothing
end
```

- [ ] **Step 7: Run the law tests to verify they pass**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T4_law.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_law.jl")'`
Expected: every testset passes, pristine; 11 spec testsets and 13 `:OnlyI` testsets. A
mismatch past 1e-20 is a bug in the engine or the combination, not a tolerance to widen:
stop and debug it (systematic-debugging).

- [ ] **Step 8: Watch the public path run past the cap**

Create `$S/phase45/probes/P5_T4_cap_fixture.jl`:

```julia
# ABOUTME: Task 4's public-path check of the summed cap: rate_equation_string on the fixture whose
# ABOUTME: conformations count 4,320 each; the head's path derives on, the fixed one throws at once.
using EnzymeRates
am = @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)          :: EqualAI
        E + B <--> E(B)          :: EqualAI
        E(A) + B <--> E(A, B)    :: EqualAI
        E(B) + A <--> E(A, B)    :: EqualAI
        E(A, B) <--> E(P, Q)     :: NonequalAI
        E(P, Q) <--> E(P) + Q    :: EqualAI
        E(P, Q) <--> E(Q) + P    :: EqualAI
        E(P) <--> E + P          :: EqualAI
        E(Q) <--> E + Q          :: EqualAI
        E(A) + Q <--> E(A, Q)    :: EqualAI
        E(Q) + A <--> E(A, Q)    :: EqualAI
        E(B) + P <--> E(B, P)    :: EqualAI
        E(P) + B <--> E(B, P)    :: EqualAI
    end
end
t = @elapsed r = try
    rate_equation_string(am)
    "derived"
catch e
    sprint(showerror, e)
end
println("rate_equation_string after ", round(t; digits = 1), " s: ", first(r, 100))
```

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T4_cap_red.log 300 --project=. $S/phase45/probes/P5_T4_cap_fixture.jl`
Expected: no `RateEquationTooLargeError`, since the old MWC path caps each conformation
alone (4,320 ≤ 5,000): `rate_equation_string after <t> s: derived` and `exit 0` (f9361cf:
4.6 s, `$S/phase45/refresh5/settle_6.log`). Record `<t>` for the commit. A
`RateEquationTooLargeError` here means the head already refuses the fixture: stop and
find out why.

- [ ] **Step 9: Cap the old MWC path on the sum**

In `src/rate_eq_derivation.jl`, in `_num_den_exprs(am::AllostericMechanism; horner =
false)`, make the first line of the body

```julia
    _assert_derivable(am)
```

and after its docstring's last sentence, "Errors when the denominator vanishes on a face of
the guard (`_assert_support`).", add the sentence "It first refuses a law whose two
conformations together exceed `MAX_RATE_EQUATION_TERMS` spanning-tree products
(`_assert_derivable`)."

Append to the testset "the cap counts both conformations" in `test/test_law.jl`, before its
`end`:

```julia
    @test_throws EnzymeRates.RateEquationTooLargeError rate_equation_string(em)
```

Run the probe again: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T4_cap.log 300 --project=. $S/phase45/probes/P5_T4_cap_fixture.jl`
Expected: `rate_equation_string after <t> s: Rate equation for this mechanism has more than
5000 polynomial terms (limit: 5000). Equations this l` and `exit 0`; `<t>` is seconds,
the guard's and the compiler's, with no derivation in it. Then rerun Step 7's command into
`P5_T4_law.log`: every testset passes.

- [ ] **Step 10: Re-pin 12 after its value check**

In `test/test_rate_eq_derivation.jl`, testset "_eq_complexity (V×τ term-count estimate)",
replace

```julia
    # Allosteric: computed on the catalytic core, so it equals the catalytic V×τ
    # (an all-SS uni-uni: G = 3 segments, τ = 3 spanning trees ⇒ 9).
```

with

```julia
    # Allosteric: both conformations count. An all-SS uni-uni has G = 3 segments and
    # τ = 3 spanning trees, 9, and its inactive conformation keeps every step: 9 more.
```

and replace

```julia
    @test EnzymeRates._eq_complexity(EnzymeRates.AllostericMechanism(allo)) == 9
```

with

```julia
    am = EnzymeRates.AllostericMechanism(allo)
    @test _testhelper_laplacian_tree_count(am.catalytic, EnzymeRates._active_mask(am)) == 9
    @test _testhelper_laplacian_tree_count(am.catalytic,
                                           EnzymeRates._inactive_mask(am)) == 9
    @test EnzymeRates._eq_complexity(am) == 18
```

The value check is the two Laplacian counts, the test's own (probe: 9 and 9).

- [ ] **Step 11: Re-pin 13 after its value check**

The split move counts by the joint solve and derives no law (`_expand_split_kinetic_group`
reads neither `_eq_complexity` nor the cap), so both testsets keep their exact children. Each
gains the per-mask counts the summed cap reads. In `test/test_mechanism_enumeration.jl`,
testset "AllostericMechanism — the split count of a law over the term cap", after

```julia
    @test EnzymeRates._eq_complexity(am) >
          EnzymeRates.MAX_RATE_EQUATION_TERMS
```

insert

```julia
    # Both conformations count: each inactive conformation keeps only E + Q <--> E(Q), two
    # segments joined by one tree, so the sum is the active law's count plus 2.
    for x in (am, p_by_a, p_by_b, a_by_p)
        @test _testhelper_laplacian_tree_count(x.catalytic,
                                               EnzymeRates._inactive_mask(x)) == 2
        @test EnzymeRates._eq_complexity(x) ==
              _testhelper_laplacian_tree_count(x.catalytic,
                                               EnzymeRates._active_mask(x)) + 2
    end
```

In testset "AllostericMechanism — a parent over the term cap splits by its count", after

```julia
    @test EnzymeRates._eq_complexity(am) > EnzymeRates.MAX_RATE_EQUATION_TERMS
```

insert

```julia
    # Both conformations count: the inactive conformation keeps only E + Q <--> E(Q).
    for x in (am, a_by_p, p_by_b)
        @test _testhelper_laplacian_tree_count(x.catalytic,
                                               EnzymeRates._inactive_mask(x)) == 2
        @test EnzymeRates._eq_complexity(x) ==
              _testhelper_laplacian_tree_count(x.catalytic,
                                               EnzymeRates._active_mask(x)) + 2
    end
```

- [ ] **Step 12: Write the filter docstring and retire the ROADMAP entry**

In `src/identify_rate_equation.jl`, in the `identify_rate_equation` docstring, replace the
`eq_complexity_filter` item

```
- `eq_complexity_filter::Int = 337`: skip (before fitting, before any
  derivation) any mechanism whose rate equation is more complex than this —
  roughly the number of terms in the rate equation's denominator: V×τ, the
  number of RE-segments times the spanning-tree count of the catalytic segment
  graph, i.e. the products the compiled equation evaluates on every call. The
  default 337 is one above a fully steady-state random-order bi-bi (V×τ = 336),
  so such a mechanism passes and anything more complex is skipped — equations
  denser than that are impractical to fit and can blow up derivation/codegen.
  Computed from the mechanism graph alone.
```

with

```
- `eq_complexity_filter::Int = 337`: skip (before fitting, before any
  derivation) any mechanism whose rate equation is more complex than this —
  roughly the number of terms in the rate equation's denominator: V×τ, the
  number of RE-segments times the spanning-tree count of the segment graph,
  i.e. the products the compiled equation evaluates on every call. An
  allosteric mechanism counts both conformations, V×τ of the active one plus
  V×τ of the inactive one. The default 337 is one above a fully steady-state
  random-order bi-bi (V×τ = 336), so such a mechanism passes and anything more
  complex is skipped — equations denser than that are impractical to fit and
  can blow up derivation/codegen. A one-subunit variant of that bi-bi whose
  inactive conformation keeps its substrate bindings counts 336 + 16 = 352 and
  is skipped. Computed from the mechanism graph alone.
```

In `ROADMAP.md`, under "Known issues", delete the bullet whose bold title is
"**`eq_complexity_filter` counts an allosteric mechanism's catalytic graph once.**" (its
text as Task 3 left it, citing `src/engine.jl`), through the line before the next bullet.

- [ ] **Step 13: Convert the population law test into Q1**

In `test/test_oracles.jl`, in the testset "Laws over a regenerated population", delete the
inner testset

```julia
    @testset "every sampled law matches the exact oracle" begin
        for (i, m) in enumerate(_testhelper_population_sample)
            _testhelper_check_law_against_oracle(EnzymeRates.compile_mechanism(m);
                                                 n_points = 2, seed = i, compiled = false)
        end
    end
```

and replace the testset "Detailed balance over a regenerated population" (from its
`@testset` line through its `end`) by Q1, which runs both checks with their seeds:

```julia
@testset "the derived laws of the population sample are exact, positive and balanced" begin
    # Q1 (spec section 14). Each sampled law, derived from the concrete mechanism at one
    # catalytic subunit, equals the exact steady state at two points, holds positive
    # coefficients in each conformation's Q, and vanishes at equilibrium, per conformation
    # and combined.
    @test count(m -> m isa EnzymeRates.AllostericMechanism,
                _testhelper_population_sample) >= 30
    for (i, m) in enumerate(_testhelper_population_sample)
        @testset "population mechanism $i" begin
            mech = _testhelper_single_subunit(m)
            s = EnzymeRates._solve(mech)
            L = EnzymeRates._law(mech, s)
            rng = Random.MersenneTwister(i)
            for _ in 1:2
                point, concs, params = _testhelper_law_inputs(mech, s.fitted, rng)
                oracle = _testhelper_steady_state_rate(mech, point, concs)
                law = _testhelper_evaluate_law(L, merge(params, concs))
                @test abs(law - oracle) <= big(10)^-20 * abs(oracle)
            end
            Qs = L isa EnzymeRates._Law ? [L.Q] : [L.active.Q, L.inactive.Q]
            @test all(Q -> all(>(0), values(Q.terms)), Qs)
            @test all(<(big(10)^-30),
                      _testhelper_detailed_balance(m, Random.MersenneTwister(i)))
        end
    end
end
```

The oracle points are the deleted inner testset's (seed `i`, fitted constants in
`fitted_params` order, one subunit) and the balance draws the deleted population testset's
(seed `i`, the mechanism at its own multiplicity).

- [ ] **Step 14: Run the touched test files**

Run (background), the four focused runs one after the other:
`S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && for f in test_rate_eq_derivation:tr test_identify_rate_equation:tid test_oracles:oracles test_mechanism_enumeration:tme; do $S/phase3/jl.sh $S/phase45/logs/P5_T4_${f#*:}.log 3000 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/'${f%%:*}'.jl")'; done`
and poll the logs in run order, each with the Runner's bounded loop:
`S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; for k in tr tid oracles tme; do L=$S/phase45/logs/P5_T4_$k.log; timeout 580 bash -c "until grep -q '^exit ' $L 2>/dev/null; do sleep 20; done"; tail -30 $L; done`
(rerun it while a log has no `exit` line yet).
Expected: every file passes, pristine. In `P5_T4_tid.log`, "eq_complexity_filter counts
both conformations" passes. In `P5_T4_oracles.log`, "Detailed balance of every spec, per
conformation and combined" still counts `n_inactive_checked == 5` and Q1 runs 120
population testsets; in `P5_T4_tme.log` both term-cap testsets keep their children.

- [ ] **Step 15: Compare the laws with Task 0's baseline**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T4_compare.log 1800 --project=. --heap-size-hint=2500M $S/phase45/scripts/P5_compare.jl $S/phase45/p5_base laws`
(the script activates the test environment itself, as Task 0 runs it).
Expected: six `laws <set>` lines, each `through _law`, the allosteric rows of the specs, the
population sample and the 13 `:OnlyI` fixtures now among them; every worst relative error
is at most 1e-20; `exit 0`.

- [ ] **Step 16: Recount R5.5 on P4-1's tagging**

Create `$S/phase45/scripts/P5_r55_recount.jl`:

```julia
# ABOUTME: R5.5 on P4-1's tagging: identify's seed tier and two expansion rounds of the regulated
# ABOUTME: bi-bi and ping-pong at one and two subunits; counts what the summed filter newly drops.
using EnzymeRates, Random, Printf
const ER = EnzymeRates
const FILTER = 337
const MechU = Union{ER.Mechanism, ER.AllostericMechanism}
const bibi2 = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    allosteric_regulators: R
    oligomeric_state: 2
end
const pingpong2 = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    allosteric_regulators: R
    oligomeric_state: 2
    shared_catalytic_site: (A, P)
end
const bibi1 = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    allosteric_regulators: R
end
const pingpong1 = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    allosteric_regulators: R
    shared_catalytic_site: (A, P)
end
"""The head's filter count: V×τ of the whole catalytic graph."""
whole(m) = ER._eq_complexity(ER.catalytic(m), trues(length(ER.catalytic(m).edges)))
function report(label, round, ms, t0)
    newly = count(m -> whole(m) <= FILTER < ER._eq_complexity(m), ms)
    @printf("%s round %d: %d candidates, %d allosteric, %d newly over %d under the sum, largest sum %d; %.0f s\n",
            label, round, length(ms), count(m -> m isa ER.AllostericMechanism, ms), newly,
            FILTER, maximum(ER._eq_complexity, ms; init = 0), time() - t0)
    flush(stdout)
end
function recount(label, rxn; rounds = 2, quota = 400)
    t0 = time()
    rng = Xoshiro(20261019)
    seeds = unique!(collect(ER.seed_mechanisms(rxn, Set([:R]), Set{Symbol}())))
    base, _ = ER._base_tier(seeds, rxn)
    report(label, 0, base, t0)
    seen = Set{UInt64}(hash.(base))
    pool = filter(m -> ER._eq_complexity(m) <= FILTER, base)
    for round in 1:rounds
        parents = length(pool) <= quota ? pool :
                  pool[sort(randperm(rng, length(pool))[1:quota])]
        kids = MechU[]
        for p in parents, c in ER.expand_mechanisms([p], rxn)
            h = hash(c)
            h in seen && continue
            push!(seen, h)
            push!(kids, c)
        end
        report(label, round, kids, t0)
        pool = filter(m -> ER._eq_complexity(m) <= FILTER, kids)
    end
end
recount("bi-bi + R, n = 2", bibi2)
recount("ping-pong + R, n = 2", pingpong2)
recount("bi-bi + R, n = 1", bibi1)
recount("ping-pong + R, n = 1", pingpong1)
```

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T4_recount.log 1800 --project=. --heap-size-hint=2500M $S/phase45/scripts/P5_r55_recount.jl`, then `timeout 580 bash -c 'until grep -q "^exit " '$S'/phase45/logs/P5_T4_recount.log; do sleep 20; done'; tail -20 $S/phase45/logs/P5_T4_recount.log`
Expected: twelve report lines and `exit 0`, about 10 min, maxrss under 3 GB. KC §1 found 0
newly over at n = 2 and 0.05–0.5% per round at n = 1; record the counts for the PR body
(VC m12). A count far above KC's is not a test failure: write it down and say so in the
commit.

- [ ] **Step 17: Run the full suite**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T4_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`, polled in bounded loops.
Expected: `exit 0`, every test passing, pristine, Aqua and JET green; about 3 s over Task 3's
time; report `maxrss_kb` against Task 0's.

- [ ] **Step 18: Commit**

Write the measured figures in place of `<…>`.

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt
git add src/engine.jl src/rate_eq_derivation.jl src/identify_rate_equation.jl ROADMAP.md \
    test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_law.jl \
    test/test_oracles.jl test/test_rate_eq_derivation.jl \
    test/test_mechanism_enumeration.jl test/test_identify_rate_equation.jl
git commit -m "$(cat <<'EOF'
Derive the MWC law on both conformations and cap them together

Each conformation's law comes from its own mask over one variable table;
free E's weights cross over their common factor. The term cap and
eq_complexity_filter count the active law plus the inactive one, and the
old MWC path refuses the same laws; before, it did not refuse the cap
fixture (<derived after <t> s, or killed at 300 s>). Pin 12 counts 9 + 9
on the test's own Laplacian; pin 13's children are unchanged, since the
split move derives no law, and gain their per-mask counts. The
detailed-balance helper and Q1 read the new law. R5.5 on P4-1's tagging:
<counts per round>. Full-suite maxrss: <kb> KB (Task 0: <kb0> KB).

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
)"
```

---
### Task 5: Test helpers and class 5 on the new law (P5-6)

Class 5 reads the new law: a plain law's `D`, or an MWC law's two branches `fI·DA` and
`fA·DI`, must keep a monomial free of each face's concentrations. The law oracle and the
identifiable rank read `_law` too, so the RE-speed and ε-path oracles (`law = :package`)
judge the new law from here on. The `POLY` methods stay for the old engine until Task 13.

**Files:**
- Modify: `src/admission.jl` (`_law_branches(::_Law)`, `_law_branches(::_AlloLaw)`,
  `_keeps_monomial(::_Poly, ::Vector{Int})`, `_vanishing_face` on `_Poly` branches and
  variable faces, `_support_verdict(m, vars, branches)`, `_assert_support(m, L)`)
- Modify: the fixtures file (`_testhelper_law_value` on `_law`; FW-07's
  `_testhelper_law_branches` becomes `_testhelper_support_verdict`)
- Modify: `test/test_mechanism_enumeration.jl` (`_testhelper_identifiable_rank` reads the
  new law, derived once per rank, through `_testhelper_evaluate_law`)
- Test: `test/test_admission.jl`, `test/test_oracles.jl` (run only)
- Logs: `$S/phase45/logs/P5_T5_{admission_red,admission,oracles,tme,full}.log`

**Interfaces:**
- Consumes: Tasks 3–4's `_law`, `_Law`, `_AlloLaw`; Task 2's `_Vars`, `_Poly`; phase 4's
  `_faces`, `_face_names`, `_Verdict`, `_WELL_POSED`, `_NO_FAULT`, `_DENOMINATOR_VANISHES`,
  `_LIMIT_UNTESTED`, `_guard_message`, `_support_message`; Task 4's
  `_testhelper_evaluate_law`.
- Produces:
  ```julia
  _law_branches(L::_Law)::Vector{Tuple{_Poly, _Poly}}       # [(L.D, one(_Poly))]
  _law_branches(A::_AlloLaw)::Vector{Tuple{_Poly, _Poly}}   # [(A.active.D, A.fI),
                                                            #  (A.inactive.D, A.fA)]
  _keeps_monomial(p::_Poly, z::Vector{Int})::Bool           # z: zeroed concentration vars
  _vanishing_face(branches::Vector{Tuple{_Poly, _Poly}}, faces::Vector{Vector{Int}})::Int
  _support_verdict(m::Mechanism, vars::_Vars,
                   branches::Vector{Tuple{_Poly, _Poly}})::_Verdict
  _assert_support(m::Union{Mechanism, AllostericMechanism},
                  L::Union{_Law, _AlloLaw})::Nothing        # class 5's throw
  # tests (fixtures file)
  _testhelper_law_value(m, params::AbstractDict, concs::AbstractDict)::BigFloat  # on _law
  _testhelper_support_verdict(mech)::EnzymeRates._Verdict   # class 5 on _law(mech)
  ```
  `_assert_support` reads only `L.vars`, which `_Law` and `_AlloLaw` type alike, and
  dispatches `_law_branches` (R18).

**Gate:**
- Phase 4's verdicts are unchanged: `test_admission.jl` and `test_oracles.jl` pass with every
  count as before (the RE-speed and ε-path oracles with `law = :package` now read the new
  law on every face they check).
- Class 5 rejects exactly `ss_dead_end_at_face` (plain and with its `:OnlyA` release),
  `bottomless_entered` and `unfed` ("class 3 rejects a bottomless segment") among its hand
  cases, each at face 2, and passes all 42 specs.
- TME's `_testhelper_identifiable_rank` gives the same ranks on every fixture that calls it.
- The full suite is green and pristine; +0 s.

- [ ] **Step 1: Write the class-5 tests on the new law**

In the fixtures file, replace FW-07's `_testhelper_law_branches(mech)` (docstring and
method) by:

```julia
"""Class 5's verdict (`EnzymeRates._support_verdict`) on the derived law of `mech`
(`EnzymeRates._law`), compiled or concrete, read without the guard."""
function _testhelper_support_verdict(mech)
    x = _testhelper_lift(mech)
    L = EnzymeRates._law(x)
    EnzymeRates._support_verdict(EnzymeRates.catalytic(x), L.vars,
                                 EnzymeRates._law_branches(L))
end
```

In `test/test_admission.jl`, replace each of FW-07's seven calls
`EnzymeRates._support_verdict(<x>, _testhelper_law_branches(<x>))` by
`_testhelper_support_verdict(<x>)`: in "class 3 rejects a bottomless segment" (`unfed`),
"the limit test decides what structure cannot" (`dead_end`, written over two lines), "the
law of a steady-state dead end vanishes at zero products" (`m`, `am`), "every spec's
denominator keeps a monomial on every face" (`x`), "a redundant steady-state binding traps
S1-bound enzyme at zero substrates" (`x`) and "the LDH i-state 6-group mechanism traps active
lactate-bound enzyme" (`x`). Then `grep -n "_testhelper_law_branches" test/` prints nothing.

In "class 3 rejects a bottomless segment", replace

```julia
    # The package law is 0/0 at zero products: class 5 refuses to derive it.
```

with

```julia
    # The package law is 0/0 at zero products: class 5 refuses it.
```

and

```julia
    @test_throws "every term of the rate law's denominator vanishes" (
        EnzymeRates._num_den_exprs(unfed))
```

with

```julia
    @test_throws "every term of the rate law's denominator vanishes" (
        EnzymeRates._assert_support(unfed, EnzymeRates._law(unfed)))
```

In "the law of a steady-state dead end vanishes at zero products", replace

```julia
    @test_throws(ErrorException(EnzymeRates._guard_message(am, w)),
                 EnzymeRates._num_den_exprs(am))
```

with

```julia
    @test_throws(ErrorException(EnzymeRates._guard_message(am, w)),
                 EnzymeRates._assert_support(am, EnzymeRates._law(am)))
```

Replace the whole testset "class 5 finds the face where a denominator loses every monomial"
by:

```julia
@testset "class 5 finds the face where a denominator loses every monomial" begin
    # Concentrations S, P and Q and a constant K; a face is the list of the concentration
    # variables it zeroes.
    vars = EnzymeRates._Vars()
    S, P, Q = (EnzymeRates._var!(vars, x; conc = true) for x in (:S, :P, :Q))
    K = EnzymeRates._var!(vars, :K)
    t(ps...) = EnzymeRates._term(EnzymeRates._monomial(ps...))
    unit = one(EnzymeRates._Poly)
    faces = [Int[], [P, Q], [S]]
    # K·S + P·Q: zeroing P and Q leaves K·S, zeroing S leaves P·Q.
    den = t(K => 1, S => 1) + t(P => 1, Q => 1)
    @test EnzymeRates._vanishing_face([(den, unit)], faces) == 0
    # P·S + Q·S: zeroing P and Q leaves no term.
    den = t(P => 1, S => 1) + t(Q => 1, S => 1)
    @test EnzymeRates._vanishing_face([(den, unit)], faces) == 2
    # Two conformations, each denominator times the other's free-enzyme weight: S·1 and
    # (1 + S)·P. Zeroing S keeps 1·P, zeroing P keeps S·1, zeroing both keeps nothing.
    den_A = t(S => 1)
    den_I = unit + t(S => 1)
    d_A = t(P => 1)
    branches = [(den_A, unit), (den_I, d_A)]
    @test EnzymeRates._vanishing_face(branches, [Int[], [S], [P]]) == 0
    @test EnzymeRates._vanishing_face(branches, [Int[], [S, P]]) == 2
    # Each polynomial by its variable names, term by term.
    _testhelper_named(p) =
        Dict(sort!([vars.names[v] => Int(e) for (v, e) in mo]; by = first) => c
             for (mo, c) in p.terms)
    @test _testhelper_named(t(K => 1, S => 1) + t(P => 1, Q => 1)) ==
          Dict([:K => 1, :S => 1] => 1, [:P => 1, :Q => 1] => 1)
    @test _testhelper_named(den_I) == Dict(Pair{Symbol, Int}[] => 1, [:S => 1] => 1)
end
```

After the testset "every spec's denominator keeps a monomial on every face", add:

```julia
@testset "class 5 reads the derived law's branches" begin
    # A plain law pairs its denominator with 1; an MWC law pairs each conformation's
    # denominator with the other conformation's free-E weight. The steady-state dead end
    # is 0/0 at zero products (face 2), plain and with E(P) + Q <--> E(P, Q) :OnlyA, and so
    # is bottomless_entered's law; Michaelis–Menten keeps a term on every face.
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
    L = EnzymeRates._law(dead_end)
    @test EnzymeRates._law_branches(L) == [(L.D, one(EnzymeRates._Poly))]
    v = EnzymeRates._support_verdict(dead_end, L.vars, EnzymeRates._law_branches(L))
    @test v == EnzymeRates._Verdict(EnzymeRates._DENOMINATOR_VANISHES,
                                    EnzymeRates._LIMIT_UNTESTED, 2, 0, 0)
    @test EnzymeRates._face_names(dead_end, EnzymeRates._faces(dead_end)[v.face]) ==
          [:P, :Q]
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ E(S)          :: EqualAI
            E(S) <--> E(P, Q)     :: EqualAI
            E(Q) + P ⇌ E(P, Q)    :: EqualAI
            E + Q <--> E(Q)       :: EqualAI
            E(P) + Q <--> E(P, Q) :: OnlyA
        end
    end
    A = EnzymeRates._law(am)
    @test EnzymeRates._law_branches(A) == [(A.active.D, A.fI), (A.inactive.D, A.fA)]
    w = EnzymeRates._support_verdict(am.catalytic, A.vars, EnzymeRates._law_branches(A))
    @test (w.fault, w.face) == (EnzymeRates._DENOMINATOR_VANISHES, 2)
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
    @test _testhelper_support_verdict(entered).face == 2
    mm = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E + P
        end
    end
    @test _testhelper_support_verdict(mm) === EnzymeRates._WELL_POSED
    @test EnzymeRates._assert_support(mm, EnzymeRates._law(mm)) === nothing
end
```

- [ ] **Step 2: Run the admission tests to verify they fail**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T5_admission_red.log 2400 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); @testset "red" begin include("test/test_admission.jl") end'`
Expected: the other testsets pass. Every testset that reads `_testhelper_support_verdict`,
`_law_branches` of a law, `_assert_support(…, _law(…))` or `_vanishing_face` on `_Poly`
branches fails or errors with a `MethodError`: no method matching
`_law_branches(::EnzymeRates._Law)` or `_law_branches(::EnzymeRates._AlloLaw)`,
`_assert_support(::EnzymeRates.Mechanism, ::EnzymeRates._Law)` (and its `AllostericMechanism`,
`_AlloLaw` form; inside `@test_throws` this shows as a wrong exception), or
`_vanishing_face(::Vector{Tuple{EnzymeRates._Poly, EnzymeRates._Poly}}, ::Vector{Vector{Int64}})`.
The outer `red` testset reports them, and the log ends `exit 1`.

- [ ] **Step 3: Write class 5 on the new law**

In `src/admission.jl`, insert after `_assert_support(m::_AnyMech, branches::Vector{Tuple{POLY,
POLY}})` and before `_support_message`:

```julia
"""
    _law_branches(L) -> Vector{Tuple{_Poly, _Poly}}

The derived law `L`'s denominator as branches of positive terms, so that it keeps a
monomial free of a face's concentrations exactly when some branch's two factors both do:
a plain law's `D` times 1; an MWC law's two conformations, each `D` times the other
conformation's free-E weight (`fI·DA` and `fA·DI`). Powers of a branch and the regulatory
sites' factors, whose 1 never vanishes, keep the same faces.
"""
_law_branches(L::_Law) = [(L.D, one(_Poly))]
_law_branches(A::_AlloLaw) = [(A.active.D, A.fI), (A.inactive.D, A.fA)]

"""Whether some monomial of `p` holds none of the variables `z` to a positive power."""
function _keeps_monomial(p::_Poly, z::Vector{Int})
    for mono in keys(p.terms)
        any(((v, e),) -> e > 0 && v in z, mono) || return true
    end
    false
end

"""The first face, each given by its zeroed concentration variables, on which no branch
of `branches` (`_law_branches`) keeps a monomial in both factors; 0 when every face keeps
one."""
function _vanishing_face(branches::Vector{Tuple{_Poly, _Poly}},
                         faces::Vector{Vector{Int}})
    for (i, z) in enumerate(faces)
        any(b -> _keeps_monomial(b[1], z) && _keeps_monomial(b[2], z), branches) ||
            return i
    end
    0
end

"""
    _support_verdict(m, vars, branches) -> _Verdict

Class 5 (spec section 8) on a law over `vars`: `_DENOMINATOR_VANISHES` on the first face
of the catalytic graph `m` (`_faces`) where the denominator `branches` (`_law_branches`)
loses every monomial, so the law is 0/0 there.
"""
function _support_verdict(m::Mechanism, vars::_Vars,
                          branches::Vector{Tuple{_Poly, _Poly}})
    faces = [Int[vars.index[x] for x in _face_names(m, f)] for f in _faces(m)]
    i = _vanishing_face(branches, faces)
    i == 0 ? _WELL_POSED : _Verdict(_DENOMINATOR_VANISHES, _LIMIT_UNTESTED, i, 0, 0)
end

"""
    _assert_support(m, L)

Error with class 5's message (`_guard_message`) when the derived law `L` of `m` is 0/0 on
some face (`_support_verdict`). `_law` does not call it, so a test can read the law of a
mechanism the guard rejects.
"""
function _assert_support(m::_AnyMech, L::Union{_Law, _AlloLaw})
    v = _support_verdict(catalytic(m), L.vars, _law_branches(L))
    v.fault === _NO_FAULT || error(_guard_message(m, v))
    nothing
end
```

- [ ] **Step 4: Run the admission tests to verify they pass**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T5_admission.log 2400 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_admission.jl")'`
Expected: every testset passes, pristine; 42 spec testsets in "every spec's denominator
keeps a monomial on every face"; every verdict assertion as in `P5_T5_admission_red.log`'s
passing testsets.

- [ ] **Step 5: Rebase the law oracle and the identifiable rank on the new law**

In the fixtures file, replace the whole of `_testhelper_law_value` (docstring and function)
by:

```julia
"""
The rate law of `m` per enzyme, evaluated in 256-bit BigFloat from its derived law
(`EnzymeRates._law`) at the parameter values `params` (by name: the fitted constants, `Keq`
and, for an allosteric mechanism, `L` and the regulator constants) and the concentrations
`concs`. A compiled `m` goes through its derivation entry, which runs the boundary check; a
concrete `m` is read as it is, so the oracles can evaluate the law of a mechanism the guard
rejects. Class 5 does not run here (`_assert_support`): a law that is 0/0 on a face gives
NaN there.
"""
function _testhelper_law_value(m, params::AbstractDict, concs::AbstractDict)
    setprecision(BigFloat, 256) do
        mech = m isa EnzymeRates.AbstractEnzymeMechanism ?
               EnzymeRates._derivation_input(typeof(m)) : m
        env = Dict{Symbol, BigFloat}(k => BigFloat(v) for (k, v) in params)
        for (k, v) in concs
            env[k] = BigFloat(v)
        end
        _testhelper_evaluate_law(EnzymeRates._law(mech), env)
    end
end
```

In `test/test_mechanism_enumeration.jl`, replace the docstring and the first lines of
`_testhelper_identifiable_rank`,

```julia
"""Finite-difference rank of ∂v/∂log θ over the fitted parameters (test oracle only). The
law of a mechanism the guard rejects is read unchecked from its polynomials
(`_num_den_exprs`) and evaluated in BigFloat; class 5 still refuses a law that is 0/0 on a
face (`_num_den_exprs` throws its message)."""
function _testhelper_identifiable_rank(m; npts = 60, ndraws = 3, h = 1e-5)
    mech = _testhelper_lift(m)
    dep, fitted = EnzymeRates._dependent_param_exprs(mech)
    rate = if EnzymeRates._is_well_posed(mech)
        em = EnzymeRates.compile_mechanism(mech)
        (c, p) -> rate_equation(em, c, p)
    else
        num, den = EnzymeRates._num_den_exprs(mech)
        function _testhelper_unchecked_rate(c, p)
            env = Dict{Symbol, BigFloat}(k => BigFloat(v)
                                         for (k, v) in pairs(merge(p, c)))
            for (k, ex) in dep
                env[k] = _testhelper_eval_expr(ex, env)
            end
            Float64(_testhelper_eval_expr(num, env) / _testhelper_eval_expr(den, env))
        end
    end
```

with

```julia
"""Finite-difference rank of ∂v/∂log θ over the fitted parameters (test oracle only). The
law of a mechanism the guard rejects is derived once, unchecked (`EnzymeRates._law`), and
evaluated in BigFloat (`_testhelper_evaluate_law`), as `_testhelper_law_value` reads it."""
function _testhelper_identifiable_rank(m; npts = 60, ndraws = 3, h = 1e-5)
    mech = _testhelper_lift(m)
    fitted = EnzymeRates._dependent_param_exprs(mech)[2]
    rate = if EnzymeRates._is_well_posed(mech)
        em = EnzymeRates.compile_mechanism(mech)
        (c, p) -> rate_equation(em, c, p)
    else
        L = EnzymeRates._law(mech)
        function _testhelper_unchecked_rate(c, p)
            env = Dict{Symbol, BigFloat}(k => BigFloat(v)
                                         for (k, v) in pairs(merge(p, c)))
            Float64(_testhelper_evaluate_law(L, env))
        end
    end
```

The rest of the helper is unchanged; the `E_total = 1` in `p` is an unused entry of the
law's environment. The law is derived once per call of the helper, as the head derived its
polynomials once: a draw evaluates it 2·60 times per fitted constant, so deriving it per
evaluation would cost thousands of derivations per rank.

- [ ] **Step 6: Run the oracle and enumeration tests**

Run (background), one after the other: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && for f in test_oracles:oracles test_mechanism_enumeration:tme; do $S/phase3/jl.sh $S/phase45/logs/P5_T5_${f#*:}.log 3000 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/'${f%%:*}'.jl")'; done`, polling each log with the Runner's bounded loop (Global Constraints).
Expected: both pass, pristine, with the pass counts of `P5_T4_oracles.log` and
`P5_T4_tme.log`. "an :OnlyI law is exact, balanced and counts its declared
dimension" and FW-23's "an admitted :OnlyI mechanism compiles to its law" (the old
engine's compiled `rate_equation` against the law, to 1e-12) now read the new law at
the same tolerances; "the guard equals the RE-speed
oracle on every face of the population sample" and the inhibitor-children faces agree with
`law = :package` on the new law; every TME rank assertion holds unchanged. A count that
moves is a bug in the rebase, not a re-pin: stop and find it.

- [ ] **Step 7: Run the full suite**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T5_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`, polled.
Expected: `exit 0`, every test passing, pristine; the time within a few seconds of Task 4's;
report `maxrss_kb`.

- [ ] **Step 8: Commit**

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt
git add src/admission.jl test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/test_admission.jl test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Run class 5 on the derived law and read it in the law oracle

Class 5 reads a plain law's denominator, or an MWC law's two conformations
each times the other's free-E weight, and keeps phase 4's verdicts: the
steady-state dead ends and bottomless_entered are 0/0 at zero products,
every spec keeps a term on every face. The law oracle, the RE-speed and
eps-path oracles and the identifiable rank read the new law with unchanged
counts. The POLY methods serve the old engine until it goes; their failure
path is reached only by laws the guard already rejects. Full-suite maxrss:
<kb> KB.

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
)"
```

---

### Task 6: Code generation on the new law; the gates (P5-21, P5-22, P5-25, P5-39, P5-7)

The generated `rate_equation` emits the derived law as straight-line locals: the
dependents first, in solve-column order, then each segment's weight and tree sums, `Q`,
`H` and the MWC terms, all in greedy Horner form, with each input read by name into a
positional local. Expr-shape tests walk every local. The speed gates, the compile budgets,
a second process and the GC gate check the result.

**Files:**
- Create: `src/codegen.jl`; modify `src/EnzymeRates.jl` (include after `engine.jl`)
- Modify: `src/rate_eq_derivation.jl` (move the `@generated rate_equation`, its docstring and
  the concrete-mechanism method to `codegen.jl`; delete the 3-argument `_rate_body`,
  `_destructuring_expr` and `_dep_assignments`; drop the `horner` keyword of both
  `_num_den_exprs` methods; ABOUTME; `_assert_no_name_clash`'s docstring)
- Modify: `src/sym_poly_for_rate_eq_derivation.jl` (delete `_poly_to_horner`)
- Modify: `docs/src/developer.md` (the two cites of the moved code; Task 15 rewrites the page)
- Modify: the fixtures file (the six gate laws move here from
  `test/test_rate_eq_derivation.jl`, unchanged; `_testhelper_check_law_against_oracle`
  gains `bigfloat`, RES §3.1)
- Test: `test/test_rate_eq_derivation.jl` (helpers, pin 16 first, then the new testsets, the
  gate's counter, the helper's barrier), `test/test_compile_budget.jl` (the subprocess
  testset)
- Create (scratch): `$S/phase45/probes/P5_T6_hoist.jl`;
  `$S/phase45/gcscan5/{derive.jl,derive_control.jl,identify.jl}`; census `census/P5_T6`
- Logs: `$S/phase45/logs/P5_T6_{dterms,tr_red,budget_red,tr,aqua_jet,hoist,tr_barrier,budget,scan112_derive,scan112_control,scan112_identify,scan113_derive,scan113_control,scan113_identify,repro112,repro113,hoistscan,census,census_cmp,full,terter,tme}.log`

**Interfaces:**
- Consumes: Task 1's `_solve`, `_Solve` (`names`, `column`, `independents`, `dependents`
  by column), `_Dependent` (`col`, `keq`, `factors`), `_Column`; Task 2's `_Vars`, `_var!`,
  `_Poly`, `_term`, `_monomial`, `_mmul`, `_mpow`, `_subst`, `_conc_shift`, `_nterms`, `+`,
  `*`; Task 3's `_bare_law(m, mask, concentrations)` (`vars`, `cols`, `sigma`, `trees`),
  `_law`, `_Law`; Task 4's `_AlloLaw`, `_flip_only`, `_law(m, s)`; Task 5's
  `_assert_support`; `_derivation_input`, `_active_mask`, `_inactive_mask`, `catalytic`.
- Produces:
  ```julia
  struct _Program
      vars::_Vars                          # inputs, then locals
      defs::Vector{Pair{Int, _Poly}}       # local v => polynomial over earlier variables
      num::_Poly
      den::_Poly
  end
  _local!(vars::_Vars, defs::Vector{Pair{Int, _Poly}}, q::_Poly)::Int
  _factor!(vars::_Vars, defs::Vector{Pair{Int, _Poly}}, p::_Poly)::_Poly
  _column_values!(vars::_Vars, defs::Vector{Pair{Int, _Poly}}, s::_Solve)::Dict{Int, _Mono}
  _state!(vars::_Vars, defs::Vector{Pair{Int, _Poly}}, m::Mechanism, s::_Solve,
          L::_Law, mask::BitVector, inactive::Bool)::Tuple{_Poly, _Poly}   # (Q, D)
  _program(m::Mechanism, s::_Solve, L::_Law)::_Program
  _program(am::AllostericMechanism, s::_Solve, A::_AlloLaw)::_Program
  _loc(v::Int, d::_Exp)::Symbol            # x<v>, x<v>m<q>, x<v>p<q>
  _power(x::Symbol, d::_Exp)               # Symbol or Expr: inv, sqrt, x^(one(x)/q)
  _horner(ts::Vector{Tuple{Int, Dict{Symbol, Int}}})       # Int, Symbol or Expr
  _emit(p::_Poly, atoms)                                   # Int, Symbol or Expr
  _body(p::_Program)::Expr
  _rate_body(m::Mechanism)::Expr                 # one method per mechanism type, so the
  _rate_body(am::AllostericMechanism)::Expr      #   law's type follows the mechanism's
  @generated rate_equation(m::M_type, concs::NamedTuple, params::NamedTuple)
      where {M_type <: AbstractEnzymeMechanism}            # _rate_body(_derivation_input(M_type))
  rate_equation(m::Union{Mechanism, AllostericMechanism}, concs, params)   # moved
  # tests
  _testhelper_rate_body(m)::Expr                            # TR
  _testhelper_body_locals(m)::Vector{Pair{Symbol, Any}}     # TR
  _testhelper_generated_rate_return(m)                      # TR: the returned Expr
  _testhelper_denominator_terms(m)::Int                     # TR: D after substitution
  _testhelper_bigfloat_body_specs::NTuple{4, String}        # TR: BigFloat-body plain specs
  _testhelper_check_law_against_oracle(m; n_points = 3, seed = 1,
                                       bigfloat = false, rtol = big(10)^-20)  # fixtures
  _testhelper_subprocess_output(script::String)::String     # test_compile_budget.jl
  ```
  The dependent locals come first, in `s.dependents` order, so a body's numbering never
  follows a hash (P5-21); a dependent equal to one symbol reads it directly. Roots are
  `sqrt` for q = 2 and `x^(one(x)/q)` for q > 2 (P5-22).

**Gate:**
- The Expr-shape tests (pin 16) are rewritten first and watched failing on the head's body,
  then pass on the new one.
- 0 allocations for every tested law; every spec under 120 ns (slice: at most 13.1 ns;
  package before: 41.6 ns); the six gate laws under `max(120 ns, 0.8 ns × D)` with the D
  pins unchanged (slice: bi-bi 42.7 ns, ter-bi 68.9 ns). A miss or an allocation stops the
  task for Denis.
- BigFloat inputs through the compiled body equal the oracle to 1e-20 on the 11
  allosteric specs, four plain ones and the cube-root fixture (RES §3.1).
- `test_aqua_jet.jl` passes in a focused run before the full suite: JET finds no
  `_rate_body` call without a method (Interface decisions, "Mechanism and law together").
- The hoisting probe reports how many specs the head's helper times under 0.5 ns (RES: 4 of
  42). Some → the helper times through a `@noinline` barrier (P5-25); none → it stays, and
  the PR body says so.
- The compile budgets hold, measured before the full run; none is raised.
- The subprocess bodies equal the parent's on the 42 specs and six gate laws.
- The GC gate (`derive`, `identify`) shows 0 hoisted allocations with a bits store into a
  reference field on 1.12.5 and 1.13.1; `derive_control` and `repro_upstream.jl` show the
  hazard on both; about 4 min per version.
- `census/P5_T6` against `census/P5_base`: `RESULT: PASS`, or a `FAIL` whose every
  difference is a law that `P5_base` derives and `P5_T6` refuses with
  `RateEquationTooLargeError`, Task 4's summed cap (Step 11's rule); those laws are named
  in the commit and the PR body. Any other difference is a bug.
- R52's three figures beside Task 0's; the full suite green and pristine; +15 s (the
  subprocess about 10 s, the BigFloat bodies about 5 s).

- [ ] **Step 1: Check the gate's term count against Task 0**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T6_dterms.log 1800 --project=. --heap-size-hint=2500M $S/phase45/scripts/P5_compare.jl $S/phase45/p5_base dterms`
(the script activates the test environment itself, as Task 0 runs it).
Expected: `dterms: 0 of 6 rows differ`: the six laws' D after substitution
(`_nterms(_law(m).D)`, which the script reads now that the tree defines `_law`) read 396,
396, 396, 589, 799 and 858, byte-equal to `dterms.tsv`; `exit 0`. This is pin 15's value
check: the counter switches in Step 6.

- [ ] **Step 2: Move the gate laws and rewrite the Expr-shape tests (pin 16) on the head's body**

Move the two constants `_testhelper_largest_enumerated_mechanisms` and
`_testhelper_large_ter_bi_mechanisms`, with the comment above each, from
`test/test_rate_eq_derivation.jl` to the end of the fixtures file, unchanged. The comment
above `_testhelper_seconds_per_denominator_term` stays in `test/test_rate_eq_derivation.jl`
and points at the laws it no longer follows: replace its first two lines,

```julia
# Seconds per call allowed per denominator term: each of the six laws above must run in
# under the larger of 120 ns and this constant times its denominator terms. It covers
```

with

```julia
# Seconds per call allowed per denominator term: each of the six gate laws (the fixtures
# file) must run in under the larger of 120 ns and this constant times its denominator
# terms. It covers
```

In `test/test_rate_eq_derivation.jl`, replace `_testhelper_generated_rate_return(m)` (the
function before the testset "rate_equation polynomial body uses 2-arg +/* calls") by:

```julia
"""The body of the generated `rate_equation` of `m`, as the generator builds it."""
function _testhelper_rate_body(m)
    mech = ER._derivation_input(typeof(m))
    dep, indep = ER._dependent_param_exprs(mech)
    ER._rate_body(mech, (indep..., :Keq, :E_total), dep)
end

"""Each local the generated `rate_equation` body of `m` binds, in order, with its
right-hand side."""
_testhelper_body_locals(m) =
    Pair{Symbol, Any}[ex.args[1] => ex.args[2] for ex in _testhelper_rate_body(m).args
                      if ex isa Expr && ex.head === :(=) && ex.args[1] isa Symbol]

"""The expression the generated `rate_equation` body of `m` returns."""
_testhelper_generated_rate_return(m) =
    only(ex.args[1] for ex in _testhelper_rate_body(m).args
         if ex isa Expr && ex.head === :return)
```

This first `_testhelper_rate_body` reads the head's body (3-argument `_rate_body`), so the
tests below are watched failing on it; Step 5 switches it.

In the testset "rate_equation polynomial body uses 2-arg +/* calls", replace

```julia
    _testhelper_walk!(_testhelper_generated_rate_return(spec.mechanism))
    @test isempty(bad)
```

with

```julia
    # Every local the body binds, and what it returns, for a plain law and an MWC law.
    mwc = only(s for s in MECHANISM_TEST_SPECS
               if s.name == "MWC Tetramer Random Bi-Bi RE + Two Allosteric Sites")
    for m in (spec.mechanism, mwc.mechanism)
        foreach(_testhelper_walk!, last.(_testhelper_body_locals(m)))
        _testhelper_walk!(_testhelper_generated_rate_return(m))
    end
    @test isempty(bad)
```

Replace the testset "the generated rate_equation body is in Horner form" by:

```julia
@testset "the generated rate_equation body is in Horner form" begin
    # Horner form multiplies by each factor once per branch instead of once per term,
    # which keeps the 396-term gate laws at tens of ns per call. Greedy Horner takes a
    # factor out of every summand it divides, so no local divides two summands of one sum,
    # in any local's definition or in what the body returns.
    for m in (_testhelper_largest_enumerated_mechanisms[2],
              only(s for s in MECHANISM_TEST_SPECS
                   if s.name == "MWC Tetramer Random Bi-Bi RE + Two Allosteric Sites"
                  ).mechanism)
        locals = _testhelper_body_locals(m)
        names = first.(locals)
        for ex in [last.(locals); _testhelper_generated_rate_return(m)]
            @test isempty(_testhelper_unfactored_sums(ex, names))
        end
    end
end
```

- [ ] **Step 3: Write the new code-generation tests**

In `test/test_rate_eq_derivation.jl`, after the testset "the generated rate_equation body is
in Horner form", add:

```julia
@testset "rate_equation reads concs and params by name in any order" begin
    # The body reads each input by name (`getfield`), so reordered NamedTuples, a data row
    # with Rate and group fields, and all_parameters' output with the dependents give the
    # same value without allocating. A metabolite named x1 reads as a concentration beside
    # the body's positional local x1.
    m = only(s for s in MECHANISM_TEST_SPECS if s.name == "Random-order Bi-Bi").mechanism
    params, concs, all_params = _testhelper_random_independent_params_concs(
        m, collect(metabolites(m)); rng = Random.MersenneTwister(91))
    v = rate_equation(m, concs, params)
    flipped(nt) = NamedTuple{reverse(keys(nt))}(reverse(values(nt)))
    row = merge((Rate = 1.0, group = 2), flipped(concs))
    for (c, p) in ((flipped(concs), params), (concs, flipped(params)),
                   (row, all_params), (row, flipped(all_params)))
        @test rate_equation(m, c, p) == v
        allocs, _ = _testhelper_test_rate_equation_performance(m, p, c)
        @test allocs == 0
    end
    x1 = @enzyme_mechanism begin
        substrates: x1
        products: P
        steps: begin
            E + x1 ⇌ E(x1)
            E(x1) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    _testhelper_check_law_against_oracle(x1)
    params, concs, _ = _testhelper_random_independent_params_concs(
        x1, collect(metabolites(x1)); rng = Random.MersenneTwister(92))
    allocs, _ = _testhelper_test_rate_equation_performance(x1, params, concs)
    @test allocs == 0
end

@testset "a cube-root dependent emits x^(one(x)/3) and is exact in BigFloat" begin
    # One isomerization group, crossed three times by the only cycle: its reverse rate is
    # the cube root of a power product (Keq^(-1/3)·…). The body computes each root once as
    # x^(one(x)/3), exact in BigFloat, where a Float64 exponent 1/3 would err by 1e-17.
    m = @enzyme_mechanism begin
        substrates: S, B
        products: P
        steps: begin
            (E <--> F, E(S) <--> F(S), E(P) <--> F(P))
            F + S <--> E(S)
            F(S) + B <--> E(P)
            F(P) <--> E + P
        end
    end
    mech = _testhelper_lift(m)
    @test any(d -> denominator(d.keq) == 3, EnzymeRates._solve(mech).dependents)
    roots = [(name, rhs) for (name, rhs) in _testhelper_body_locals(m)
             if occursin(r"^x\d+[mp]3$", String(name))]
    @test !isempty(roots)
    for (name, rhs) in roots
        x = Symbol(match(r"^(x\d+)", String(name)).captures[1])
        root = :($x^(one($x) / 3))
        @test rhs == (endswith(String(name), "p3") ? root : :(inv($root)))
    end
    fitted = EnzymeRates.fitted_params(m)
    rng = Random.MersenneTwister(93)
    setprecision(BigFloat, 256) do
        for _ in 1:3
            point, concs, params = _testhelper_law_inputs(mech, fitted, rng)
            oracle = _testhelper_steady_state_rate(mech, point, concs)
            p = NamedTuple{(fitted..., :Keq, :E_total)}(
                (getindex.(Ref(params), fitted)..., params[:Keq], big(1)))
            c = NamedTuple{metabolites(m)}(Tuple(concs[x] for x in metabolites(m)))
            @test abs(rate_equation(m, c, p) - oracle) <= big(10)^-20 * abs(oracle)
        end
    end
end

@testset "dependent locals are numbered by solve column" begin
    # The program defines a local for each dependent that is a power product of two or
    # more symbols before any other local, in the solve's column order, so a body's
    # numbering never follows a hash (spec section 2, criterion 3); a dependent equal to
    # one symbol takes no local. The random-order bi-bi has a Haldane dependent and a
    # Wegscheider one, each a product of several constants.
    m = _testhelper_lift(only(s for s in MECHANISM_TEST_SPECS
                              if s.name == "Random-order Bi-Bi").mechanism)
    s = EnzymeRates._solve(m)
    L = EnzymeRates._law(m, s)
    p = EnzymeRates._program(m, s, L)
    @test issorted([d.col for d in s.dependents])
    function _testhelper_power(d, s, L)
        EnzymeRates._monomial(L.vars.index[:Keq] => d.keq,
                              (L.vars.index[s.names[f]] => e for (f, e) in d.factors)...)
    end
    products = unique!([_testhelper_power(d, s, L) for d in s.dependents
                        if length(_testhelper_power(d, s, L)) >= 2])
    @test length(s.dependents) == length(products) == 2
    @test [q for (_, q) in p.defs[1:length(products)]] == EnzymeRates._term.(products)
    # Task 1's labels fixture: a kinetic-group dependent and a Wegscheider one that each
    # equal one constant, read directly, and a Haldane product, the program's first local.
    # Every later local holds a sum of two or more terms.
    folded = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
            E + B ⇌ E(B)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
    sf = EnzymeRates._solve(folded)
    Lf = EnzymeRates._law(folded, sf)
    pf = EnzymeRates._program(folded, sf, Lf)
    @test [d.label for d in sf.dependents] == [:haldane, :kinetic_group, :wegscheider]
    haldane = only(d for d in sf.dependents if d.label === :haldane)
    @test first(pf.defs).second == EnzymeRates._term(_testhelper_power(haldane, sf, Lf))
    @test count(d -> EnzymeRates._nterms(d.second) == 1, pf.defs) == 1
end

@testset "the generated body is finite at zero concentrations" begin
    # Every term of the program holds each concentration to a nonnegative power, so with
    # one concentration at zero the body returns the law's value there, which the law in
    # lowest terms gives exactly; "Zero-metabolite finiteness" checks finiteness alone.
    for spec in MECHANISM_TEST_SPECS
        @testset "$(spec.name)" begin
            m = spec.mechanism
            L = EnzymeRates._law(_testhelper_lift(m))
            rng = Random.MersenneTwister(778 + hash(spec.name) % 1000)
            params = _testhelper_random_independent_params(m; rng)
            mets = collect(metabolites(m))
            for zeroed in mets
                concs = NamedTuple{Tuple(mets)}(
                    Tuple(n == zeroed ? 0.0 : 0.5 + rand(rng) for n in mets))
                v = rate_equation(m, concs, params)
                env = Dict{Symbol, BigFloat}(k => BigFloat(x)
                                             for (k, x) in pairs(merge(params, concs)))
                limit = params.E_total * _testhelper_evaluate_law(L, env)
                @test isfinite(v) && isapprox(v, Float64(limit); rtol = 1e-8)
            end
        end
    end
end
```

In the fixtures file, in `_testhelper_check_law_against_oracle` (RES §3.1: BigFloat
inputs through the compiled body), replace the end of its docstring and its signature,

```julia
its dependents from them, so a wrong dependent shows as a mismatch. The compiled
`rate_equation` is checked too, in Float64 to 1e-8. An allosteric mechanism is checked
with one catalytic subunit.
"""
function _testhelper_check_law_against_oracle(m; n_points = 3, seed = 1,
                                              rtol = big(10)^-20)
```

with

```julia
its dependents from them, so a wrong dependent shows as a mismatch. The compiled
`rate_equation` is checked too, in Float64 to 1e-8, and with `bigfloat` also with 256-bit
BigFloat inputs, to `rtol`. An allosteric mechanism is checked with one catalytic subunit.
"""
function _testhelper_check_law_against_oracle(m; n_points = 3, seed = 1,
                                              bigfloat = false, rtol = big(10)^-20)
```

and after its line
`        @test isapprox(rate_equation(em, c64, p64), Float64(oracle); rtol = 1e-8)` add

```julia
        bigfloat || continue
        pbig = NamedTuple{(fitted..., :Keq, :E_total)}(
            (getindex.(Ref(params), fitted)..., point.Keq, big(1)))
        cbig = NamedTuple{mets}(Tuple(concs[x] for x in mets))
        @test abs(rate_equation(em, cbig, pbig) - oracle) <= rtol * abs(oracle)
```

In `test/test_rate_eq_derivation.jl`, replace `_testhelper_test_mass_action_oracle` (its
docstring and function) with:

```julia
# The plain specs whose compiled body also reads BigFloat inputs in the mass-action check,
# beside every allosteric spec: an isomerization, a ping-pong, a ter-bi and a law with
# both a Haldane and a Wegscheider dependent.
const _testhelper_bigfloat_body_specs = ("Segel Iso Uni Uni", "Segel Ping Pong Bi Bi",
                                         "Segel Ordered Ter Bi", "Random-order Bi-Bi")

"""The spec's rate law against the exact steady state of its mass-action network
(`_testhelper_check_law_against_oracle`); an allosteric spec runs with one catalytic
subunit, where the MWC law is the exact steady state. The allosteric specs and
`_testhelper_bigfloat_body_specs` also run 256-bit BigFloat inputs through the compiled
body, to 1e-20."""
function _testhelper_test_mass_action_oracle(spec::_testhelper_MechanismTestSpec)
    @testset "Mass-action oracle" begin
        bigfloat = spec.mechanism isa EnzymeRates.AllostericEnzymeMechanism ||
                   spec.name in _testhelper_bigfloat_body_specs
        _testhelper_check_law_against_oracle(spec.mechanism; bigfloat)
    end
end
```

Each BigFloat check compiles one more `rate_equation` method per spec (15 in all), about
5 s; Step 7 records the TR run's time against Task 5's. The other callers of the helper
keep Float64 alone (RES §3.1 asks for BigFloat everywhere; the 15 specs and the cube-root
fixture bound the compile cost, a deviation listed for Denis).

In the testset "the reserved names cover every non-metabolite name of a law", replace

```julia
        mech = EnzymeRates._derivation_input(typeof(m))
        dep, indep = EnzymeRates._dependent_param_exprs(mech)
        body = EnzymeRates._rate_body(mech, (indep..., :Keq, :E_total), dep)
        printed = Meta.parseall(rate_equation_string(m))
        others = setdiff(_testhelper_names_read(body) ∪ _testhelper_names_read(printed),
                         metabolites(m))
        @test all(n -> occursin('_', String(n)) || n in reserved, others)
        # The law destructures params and concs and reads Keq; the printed one assigns v.
```

with

```julia
        printed = Meta.parseall(rate_equation_string(m))
        others = setdiff(_testhelper_names_read(printed), metabolites(m))
        @test all(n -> occursin('_', String(n)) || n in reserved, others)
        # The printed law destructures params and concs, reads Keq and assigns v.
```

and fix its comment above the testset, "Besides the metabolites, the generated
rate_equation body and the printed law bind or read only such names", to read "Besides the
metabolites, the printed law binds or reads only such names; the generated body binds only
positional locals (the next testset)". After that testset add:

```julia
@testset "the generated body binds positional locals and reads inputs by getfield" begin
    # Each local is x<i>, its reciprocal x<i>m<q> or its root x<i>p<q>, and each input is
    # read as getfield(concs, :X) of a metabolite or getfield(params, :k) of a parameter,
    # so no metabolite name can shadow a constant or a local.
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        body = _testhelper_rate_body(m)
        locals = _testhelper_body_locals(m)
        @test all(p -> occursin(r"^x\d+([mp]\d+)?$", String(p.first)), locals)
        @test _testhelper_names_read(body) ⊆ union(Set(first.(locals)), [:concs, :params])
        reads = Expr[]
        function _testhelper_getfields!(e)
            e isa Expr || return
            e.head === :call && e.args[1] === :getfield && push!(reads, e)
            foreach(_testhelper_getfields!, e.args)
        end
        _testhelper_getfields!(body)
        @test all(r -> r.args[2] in (:concs, :params) && r.args[3] isa QuoteNode, reads)
        @test Set(r.args[3].value for r in reads if r.args[2] === :concs) ⊆
              Set(metabolites(m))
        @test Set(r.args[3].value for r in reads if r.args[2] === :params) ⊆
              Set(parameters(m))
    end
end
```

In `test/test_compile_budget.jl`, replace `_testhelper_measure_labeled_subprocess` with the
two helpers

```julia
# Runs `script` in a fresh Julia subprocess on the test environment, so the script can
# load test dependencies, and returns what it printed, or "" when it fails.
function _testhelper_subprocess_output(script::String)
    julia_exe = Base.julia_cmd().exec[1]
    out_buf = IOBuffer()
    try
        run(pipeline(Cmd([julia_exe, "--project=$(Base.active_project())", "-e", script]);
                     stdout=out_buf, stderr=devnull); wait=true)
    catch e
        @warn "subprocess failed: $e"
        return ""
    end
    String(take!(out_buf))
end

# Runs `script` as `_testhelper_subprocess_output` does; the script is expected to print
# `<label>:<float>` for each label in `labels`. Returns a Vector{Float64} parallel to
# `labels` (NaN for any label not found or on subprocess failure).
_testhelper_measure_labeled_subprocess(script::String, labels::Vector{String}) =
    _testhelper_parse_labeled(_testhelper_subprocess_output(script), labels)
```

and append after the testset "compile-budget":

```julia
@testset "a second process derives the same rate_equation bodies" begin
    # Q7 (spec section 2, criterion 3): a fresh process reads each mechanism back from its
    # text and derives the same generated body, compared by SHA-256, for the 42 specs and
    # the six gate laws.
    mechs = [_testhelper_lift(s.mechanism) for s in MECHANISM_TEST_SPECS]
    append!(mechs, _testhelper_lift.(_testhelper_largest_enumerated_mechanisms))
    append!(mechs, _testhelper_lift.(_testhelper_large_ter_bi_mechanisms))
    @test length(mechs) == 48
    file = tempname()
    write(file, join([sprint(show, m) for m in mechs], "\n----\n"))
    script = """
        using EnzymeRates
        const ER = EnzymeRates
        for (i, text) in enumerate(split(read($(repr(file)), String), "\\n----\\n"))
            m = ER._parse_mechanism(String(text))
            println("BODY ", i, " ", bytes2hex(ER.sha256(string(ER._rate_body(m)))))
        end
        """
    lines = filter(startswith("BODY "), split(_testhelper_subprocess_output(script), '\n'))
    rm(file; force = true)
    @test length(lines) == 48
    for (i, m) in enumerate(mechs)
        mine = bytes2hex(EnzymeRates.sha256(string(EnzymeRates._rate_body(m))))
        @test "BODY $i $mine" in lines
    end
end
```

- [ ] **Step 4: Run the tests to verify they fail**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T6_tr_red.log 3600 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); @testset "red" begin include("test/test_rate_eq_derivation.jl") end' && $S/phase3/jl.sh $S/phase45/logs/P5_T6_budget_red.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_compile_budget.jl")'`, polled.
Expected in `P5_T6_tr_red.log`:
- "rate_equation polynomial body uses 2-arg +/* calls" fails: the head's dependent
  assignments are n-ary products (`Expr(:call, :*, terms...)` from `build_power_expr`);
- the cube-root testset fails: no `x<i>p3` local, and the head's Float64 exponent
  `-0.3333333333333333` misses the BigFloat oracle by about 1e-17;
- "dependent locals are numbered by solve column" errors with
  ``UndefVarError: `_program` not defined in `EnzymeRates` ``;
- "the generated body binds positional locals …" fails: the head binds `k_…` and
  metabolite names;
- the Horner testset, "the reserved names …", "rate_equation reads concs and params by name
  in any order" and "the generated body is finite at zero concentrations" pass: the head's
  body already has these properties, and the tests hold the new body to them;
- the 15 "Mass-action oracle" testsets with BigFloat inputs pass on the head's body when
  the spec's dependents carry integer or half-integer exponents, which the head writes as
  exact Float64 literals; a spec with a third-root exponent misses by about 1e-17, as the
  cube-root testset does. Record which specs miss, if any;
- everything else passes, and the outer `red` testset ends the log with `exit 1`.

In `P5_T6_budget_red.log` the compile-budget testsets pass and the
second-process testset fails `length(lines) == 48` (the subprocess hits `MethodError: no
method matching _rate_body(::EnzymeRates.Mechanism)`, and the helper warns that it failed),
then errors on the same `MethodError` in the parent's first `_rate_body(m)`.

- [ ] **Step 5: Write the code generator**

Create `src/codegen.jl`:

```julia
# ABOUTME: Code generation: the generated rate_equation emits the derived law as
# ABOUTME: straight-line locals in greedy Horner form, reading its inputs by name.

"""
A rate law as straight-line code, `v = E_total·num/den`. `vars` lists the law's inputs
(concentrations and parameters, read by name), then its locals; `defs[i] = v => p`
defines local `v` as the polynomial `p` over the inputs and earlier locals. Every term of
`num` and `den` holds each concentration to a nonnegative power.
"""
struct _Program
    vars::_Vars
    defs::Vector{Pair{Int, _Poly}}
    num::_Poly
    den::_Poly
end

"""The local of `vars` that `defs` defines as `q`, defined anew unless one already is."""
function _local!(vars::_Vars, defs::Vector{Pair{Int, _Poly}}, q::_Poly)
    i = findfirst(d -> d.second == q, defs)
    i === nothing || return defs[i].first
    push!(defs, _var!(vars, Symbol("#", length(defs) + 1)) => q)
    defs[end].first
end

"""`p` as one term: `p` itself when it has at most one, else the monomial of its terms'
lowest concentration exponents times a local holding the rest, which is then in lowest
terms and so finite at zero concentrations."""
function _factor!(vars::_Vars, defs::Vector{Pair{Int, _Poly}}, p::_Poly)
    _nterms(p) <= 1 && return p
    s = _conc_shift((p,), vars)
    _term(_mmul(_monomial(_local!(vars, defs, p * s) => 1), _mpow(s, -1)))
end

"""Each solve column's value as a monomial over `vars`: an independent's own variable; a
dependent's power product of `Keq` and independents, read directly when it is one symbol,
else from a local (`_local!`). The locals are defined in `s.dependents` order, by column,
so a body's numbering never follows a hash."""
function _column_values!(vars::_Vars, defs::Vector{Pair{Int, _Poly}}, s::_Solve)
    value = Dict{Int, _Mono}(c => _monomial(vars.index[s.names[c]] => 1)
                             for c in s.independents)
    for d in s.dependents
        r = _monomial(vars.index[:Keq] => d.keq,
                      (vars.index[s.names[f]] => e for (f, e) in d.factors)...)
        value[d.col] = length(r) <= 1 ? r : _monomial(_local!(vars, defs, _term(r)) => 1)
    end
    value
end

"""
    _state!(vars, defs, m, s, L, mask, inactive) -> (Q, D)

The terms over `vars` of the law `L` of the edges `mask` of `m`, in the inactive
conformation's constants when `inactive`, with locals added to `defs`: each bare constant
reads its column's value (`_column_values!`), each segment's weight sum `σ_g` and tree sum
`T_g` and `Q` get a local (`_factor!`). Returns `Q` as one term and
`D = shift·Σ_g σ_g·T_g`.
"""
function _state!(vars::_Vars, defs::Vector{Pair{Int, _Poly}}, m::Mechanism, s::_Solve,
                 L::_Law, mask::BitVector, inactive::Bool)
    value = _column_values!(vars, defs, s)
    b = _bare_law(m, mask, L.vars.names[L.vars.conc])
    nc = count(L.vars.conc)
    bare(v) = v <= nc ? nothing : value[s.column[_Column(b.cols[v - nc]..., inactive)]]
    shift = _term(L.shift)
    D = zero(_Poly)
    for (σ, T) in zip(b.sigma, b.trees)
        D = D + shift * _factor!(vars, defs, _subst(bare, σ)) *
                _factor!(vars, defs, _subst(bare, T))
    end
    _factor!(vars, defs, L.Q), D
end

"""The program of a plain law: `num = H·Q`, `den = D` (`_state!`)."""
function _program(m::Mechanism, s::_Solve, L::_Law)
    vars, defs = deepcopy(L.vars), Pair{Int, _Poly}[]
    Q, D = _state!(vars, defs, m, s, L, trues(length(m.edges)), false)
    _Program(vars, defs, _factor!(vars, defs, L.H) * Q, D)
end

"""
The program of an MWC law: each conformation's terms (`_state!`; 0 and 1 for an inactive
conformation that is free E alone), locals for the free-E weights `fA`, `fI`, the
conformation terms `xA = fI·DA` and `xI = fA·DI` and each site's two binding polynomials;
`num = H·(fI·QA·xA^(n−1)·WA + L·fA·QI·xI^(n−1)·WI)`, `den = xA^n·WA + L·xI^n·WI`.
"""
function _program(am::AllostericMechanism, s::_Solve, A::_AlloLaw)
    m, n = catalytic(am), A.n
    vars, defs = deepcopy(A.vars), Pair{Int, _Poly}[]
    QA, DA = _state!(vars, defs, m, s, A.active, _active_mask(am), false)
    QI, DI = _flip_only(A.inactive) ? (zero(_Poly), one(_Poly)) :
             _state!(vars, defs, m, s, A.inactive, _inactive_mask(am), true)
    f(p) = _factor!(vars, defs, p)
    pw(p, k) = foldl(*, (p for _ in 1:k); init = one(_Poly))
    fA, fI = f(A.fA), f(A.fI)
    xA, xI = f(fI * DA), f(fA * DI)
    WA = foldl(*, (pw(f(w), k) for (w, _, k) in A.sites); init = one(_Poly))
    WI = foldl(*, (pw(f(w), k) for (_, w, k) in A.sites); init = one(_Poly))
    L = _term(_monomial(vars.index[:L] => 1))
    num = f(A.H) * (fI * QA * pw(xA, n - 1) * WA + L * fA * QI * pw(xI, n - 1) * WI)
    _Program(vars, defs, num, pw(xA, n) * WA + L * pw(xI, n) * WI)
end

"""The local holding variable `v` to the unit power `d` (`1`, `-1`, `1//q` or `-1//q`):
`x<v>`, its reciprocal `x<v>m1`, its root `x<v>p<q>` or that root's reciprocal
`x<v>m<q>`."""
_loc(v::Int, d::_Exp) =
    d == 1 ? Symbol(:x, v) : Symbol(:x, v, d < 0 ? :m : :p, denominator(d))

"""`x` to the unit power `d`: a reciprocal (`inv`), a square root, or for `q > 2` the root
`x^(one(x)/q)`, exact in any float type."""
function _power(x::Symbol, d::_Exp)
    q = denominator(d)
    base = q == 1 ? x : q == 2 ? :(sqrt($x)) : :($x^(one($x) / $q))
    d < 0 ? :(inv($base)) : base
end

"""Greedy Horner form of the terms `ts`, each `(coefficient, local => power)`: take out
the local most terms hold (ties by name), then recurse on the terms with it and on the
rest, with binary `+` and `*` only."""
function _horner(ts::Vector{Tuple{Int, Dict{Symbol, Int}}})
    c0 = sum((c for (c, a) in ts if isempty(a)); init = 0)
    rest = Tuple{Int, Dict{Symbol, Int}}[(c, a) for (c, a) in ts if !isempty(a)]
    isempty(rest) && return c0
    held = Dict{Symbol, Int}()
    for (_, a) in rest, x in keys(a)
        held[x] = get(held, x, 0) + 1
    end
    best = first(sort!(collect(keys(held)); by = x -> (-held[x], String(x))))
    lower(a) = (b = copy(a); b[best] == 1 ? delete!(b, best) : (b[best] -= 1); b)
    inner = _horner(Tuple{Int, Dict{Symbol, Int}}[(c, lower(a))
                                                  for (c, a) in rest if haskey(a, best)])
    ex = inner == 1 ? best : :($best * $inner)
    others = Tuple{Int, Dict{Symbol, Int}}[(c, a) for (c, a) in rest if !haskey(a, best)]
    iszero(c0) || push!(others, (c0, Dict{Symbol, Int}()))
    isempty(others) ? ex : :($ex + $(_horner(others)))
end

"""`p` as code, `atoms(mono)` giving each monomial's locals and their powers: the Horner
form of its positive terms minus that of its negative ones."""
function _emit(p::_Poly, atoms)
    ts = Tuple{Int, Dict{Symbol, Int}}[(c, atoms(mo))
                                       for (mo, c) in sort!(collect(p.terms); by = first)]
    pos = Tuple{Int, Dict{Symbol, Int}}[t for t in ts if t[1] > 0]
    neg = Tuple{Int, Dict{Symbol, Int}}[(-c, a) for (c, a) in ts if c < 0]
    isempty(neg) ? _horner(pos) : isempty(pos) ? :(-$(_horner(neg))) :
        :($(_horner(pos)) - $(_horner(neg)))
end

"""
    _body(p) -> Expr

The body of `rate_equation` for the program `p`: each variable it reads, in variable
order, bound to a positional local, an input read by name (`getfield(concs, :X)`,
`getfield(params, :k)`) or a local in greedy Horner form, followed by the reciprocals and
roots the body reads (`_power`); then `E_total·num/den`. `x^e` reads `x` to the integer
part of `e` and `x^(±1/q)` to the rest (`x^(3/2) = x·sqrt(x)`).
"""
function _body(p::_Program)
    used = Set{Tuple{Int, _Exp}}()
    function atoms(mo::_Mono)
        a = Dict{Symbol, Int}()
        for (v, e) in mo
            n = trunc(Int, e)
            r = e - n
            for (d, k) in ((sign(n) // 1, abs(n)),
                           (sign(r) // denominator(r), numerator(abs(r))))
                k == 0 && continue
                push!(used, (v, d))
                a[_loc(v, d)] = k
            end
        end
        a
    end
    def = Dict{Int, Any}(v => _emit(q, atoms) for (v, q) in p.defs)
    num, den = _emit(p.num, atoms), _emit(p.den, atoms)
    lines = Expr[]
    for v in sort!(unique!([first(u) for u in used]))
        x = get(def, v) do
            :(getfield($(p.vars.conc[v] ? :concs : :params), $(QuoteNode(p.vars.names[v]))))
        end
        push!(lines, :($(_loc(v, 1 // 1)) = $x))
        for d in sort!([d for (w, d) in used if w == v && d != 1])
            push!(lines, :($(_loc(v, d)) = $(_power(_loc(v, 1 // 1), d))))
        end
    end
    Expr(:block, lines..., :(return getfield(params, :E_total) * $num / $den))
end

"""The body of the generated `rate_equation` of `m`: its solve, its law, class 5's check
of that law (`_assert_support`), then the program's code (`_body`). One method per
mechanism type, so the law's type follows the mechanism's and `_program` dispatches on a
matching pair."""
function _rate_body(m::Mechanism)
    s = _solve(m)
    L = _law(m, s)
    _assert_support(m, L)
    _body(_program(m, s, L))
end

function _rate_body(am::AllostericMechanism)
    s = _solve(am)
    A = _law(am, s)
    _assert_support(am, A)
    _body(_program(am, s, A))
end

# ─── rate_equation ─────────────────────────────────────────────────────────────
```

Then move, from `src/rate_eq_derivation.jl` to the end of `src/codegen.jl`, the section from
the docstring `"""\n    rate_equation(m, concs, params)` through the concrete-mechanism
method `rate_equation(m::Union{Mechanism, AllostericMechanism}, concs, params) =
rate_equation(compile_mechanism(m), concs, params)` with its comment, and change it there to:

```julia
"""
    rate_equation(m, concs, params)

Return the net reaction rate of mechanism `m` at the metabolite concentrations
`concs` and parameters `params` (both `NamedTuple`s). The rate equation itself
is derived from the mechanism's elementary steps by the King–Altman/Cha method:
rapid-equilibrium steps collapse to binding constants, steady-state steps are
assembled into the King–Altman determinant, and the two combine into the
quasi-steady-state flux through the whole mechanism — so one call returns the
overall turnover, not the rate of a single step.

`params` holds the independent constants (`parameters(m)`: the fitted constants, `Keq`
and `E_total`); the rate law computes each dependent constant from them by its
Haldane or Wegscheider relation. Both are read by name, so their fields may come in any
order, and fields the law does not read are ignored.

The derivation runs once, at compile time: the body is generated as straight-line
arithmetic with no allocations, loops, or matrix operations. Use
[`rate_equation_string`](@ref) to inspect the derived equation symbolically.
"""
@generated function rate_equation(m::M_type, concs::NamedTuple,
        params::NamedTuple) where {M_type <: AbstractEnzymeMechanism}
    # The generator decodes the mechanism and checks its rate law once
    # (`_derivation_input`).
    _rate_body(_derivation_input(M_type))
end

# Concrete-mechanism convenience: lift to the singleton (allocates, not the hot path).
rate_equation(m::Union{Mechanism, AllostericMechanism}, concs, params) =
    rate_equation(compile_mechanism(m), concs, params)
```

In `src/rate_eq_derivation.jl`:
- delete `_destructuring_expr`, `_dep_assignments` and the 3-argument `_rate_body`, each with
  its docstring, and the now-empty section line `# ─── rate_equation ───…`;
- in `_num_den_exprs(m::Mechanism; horner = false)` make the signature
  `_num_den_exprs(m::Mechanism)`, delete the line
  `horner && return _poly_to_horner(num), _poly_to_horner(den)`, and write its docstring
  as "The numerator and denominator Exprs of the rate equation of `m`, from its raw symbolic
  rate polys, as flat sums for the printed equation. Errors when the denominator vanishes
  on a face of the guard (`_assert_support`).";
- in `_num_den_exprs(am::AllostericMechanism; horner = false)` make the signature
  `_num_den_exprs(am::AllostericMechanism)`, make `to_expr(p) = _poly_to_expr(p,
  cat_params)`, and in its docstring replace "written with flat sums for the printed
  equation, or in Horner form (`_poly_to_horner`) for the generated code when `horner` is
  true." with "written with flat sums for the printed equation.";
- in `_assert_no_name_clash`'s docstring replace "The rate equation binds the parameters"
  with "The printed rate equation binds the parameters";
- replace the two ABOUTME lines with
  ```julia
  # ABOUTME: King-Altman/Cha rate polynomials of the printed rate equation, the parameter API
  # ABOUTME: and the per-conformation constant names; MWC assembly of the printed law.
  ```

In `src/sym_poly_for_rate_eq_derivation.jl`, delete `_poly_to_horner` with its docstring.

In `src/EnzymeRates.jl`, add `include("codegen.jl")` directly after `include("engine.jl")`.

In `docs/src/developer.md`, "Derivation architecture": replace "in `@generated` functions
(`src/rate_eq_derivation.jl`)" with "in `@generated` functions (`src/codegen.jl`)", and
"`rate_equation(m, conc, params)` as one arithmetic expression, written in Horner\nform
(`_poly_to_horner`), that must be" with "`rate_equation(m, conc, params)` as straight-line
arithmetic, written in Horner\nform (`_horner`), that must be".

In `test/test_rate_eq_derivation.jl`, switch `_testhelper_rate_body` to the new body:

```julia
"""The body of the generated `rate_equation` of `m`, as the generator builds it."""
_testhelper_rate_body(m) = ER._rate_body(ER._derivation_input(typeof(m)))
```

Run `grep -rn "_poly_to_horner\|_dep_assignments\|_destructuring_expr\|horner = \|horner=" src test docs/src`: it prints nothing.

- [ ] **Step 6: Count the gate's terms on the new law (pin 15)**

In `test/test_rate_eq_derivation.jl`, replace

```julia
"""The number of terms in the denominator of `m`'s derived rate law."""
_testhelper_denominator_terms(m) =
    length(EnzymeRates._raw_symbolic_rate_polys(typeof(m))[2])
```

with

```julia
"""The number of terms in the denominator of `m`'s derived rate law after the solve is
substituted: an MWC law counts both conformations, an inactive conformation that is free E
alone counting 1."""
function _testhelper_denominator_terms(m)
    L = EnzymeRates._law(_testhelper_lift(m))
    L isa EnzymeRates._Law && return EnzymeRates._nterms(L.D)
    EnzymeRates._nterms(L.active.D) + EnzymeRates._nterms(L.inactive.D)
end
```

The pins `[396, 396, 396]` and `[589, 799, 858]` stay (Step 1's value check) and
`_testhelper_seconds_per_denominator_term` stays `0.8e-9` until Task 14.

- [ ] **Step 7: Run the derivation tests to verify they pass**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T6_tr.log 3600 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_rate_eq_derivation.jl")'`, polled.
Expected: every testset passes, pristine: every "Performance" testset `allocs == 0` and
`t < 120e-9`; the size gate prints six `[ Info: rate_equation: <t> ns per call, <D>
denominator terms, <c> ns per term` lines, with D 396, 396, 396, 589, 799 and 858, and its
`allocs == 0` and `t < max(120e-9, 0.8e-9 * D)` tests pass (slice: 42.7 ns bi-bi, 68.9 ns
ter-bi). Record the six ns-per-term readings and the run's `elapsed` for the commit. An
allocation or a missed bound stops the task: report the law, the reading and the body to
Denis before any other change (CLAUDE.md).

Then run JET before the full suite, which is the first run that includes it: the focused
command with `include("test/test_aqua_jet.jl")`, log `P5_T6_aqua_jet.log`, timeout 1800.
Expected: the "Aqua" and "JET" testsets pass, `exit 0`. A JET report of a call with no
matching method (a `(2/4 union split)` line) is a bug in this task's code: give the caller
one method per mechanism type (Interface decisions, "Mechanism and law together"), never
filter the report.

- [ ] **Step 8: Measure hoisting in the gate's helper**

Create `$S/phase45/probes/P5_T6_hoist.jl`:

```julia
# ABOUTME: Task 6's hoisting probe: times each spec's rate_equation as the gate's helper does,
# ABOUTME: inline in its loop, and through a @noinline barrier, on the package's generated code.
using Test, EnzymeRates, LinearAlgebra, Random, Printf
include(joinpath(pwd(), "test", "mechanism_definitions_for_test_enzyme_derivation.jl"))
@noinline barrier(m, concs, params) = rate_equation(m, concs, params)
"""Seconds per call of `f`, the minimum over 5 batches of 10,000, as the helper times."""
function per_call(f::F, m, params, concs) where {F}
    best, acc = Inf, 0.0
    for _ in 1:5
        acc = 0.0
        t = @elapsed for _ in 1:10_000
            acc += f(m, concs, params)
        end
        best = min(best, t / 10_000)
    end
    isfinite(acc) || error("non-finite rate")
    best
end
rows = map(MECHANISM_TEST_SPECS) do spec
    m = spec.mechanism
    rng = Random.MersenneTwister(42)
    params = _testhelper_random_independent_params(m; rng)
    mets = collect(metabolites(m))
    concs = NamedTuple{Tuple(mets)}(Tuple(0.1 + 9.9 * rand(rng) for _ in mets))
    rate_equation(m, concs, params)
    barrier(m, concs, params)
    (spec.name, per_call(rate_equation, m, params, concs), per_call(barrier, m, params, concs))
end
for (name, inline, through) in rows
    @printf("%-62s inline %7.2f ns   barrier %7.2f ns\n", name, 1e9 * inline, 1e9 * through)
end
println("under 0.5 ns: inline ", count(r -> r[2] < 0.5e-9, rows), " of ", length(rows),
        ", barrier ", count(r -> r[3] < 0.5e-9, rows))
```

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T6_hoist.log 900 --project=. -e 'using TestEnv; TestEnv.activate(); include("'$S'/phase45/probes/P5_T6_hoist.jl")'`
Expected: 42 lines and the summary `under 0.5 ns: inline <h> of 42, barrier 0`. RES predicts
`<h>` = 4 on slice bodies. Record the two slowest specs and their readings (the barrier
column; the inline one when `<h>` = 0) for the commit and for Task 14's CLAUDE.md edit,
whose headroom sentence names 7de259f's figures (about 55 ns on "LDH i-state NonequalAI
6-group, steady-state NAD group", about 30 ns on PFK-1 and the MWC tetramer).

If `<h>` > 0, the helper times calls the compiler hoisted out of its loop: make it time
through a barrier (P5-25). In `test/test_rate_eq_derivation.jl` replace
`_testhelper_test_rate_equation_performance` (docstring and function) with:

```julia
"""The call the gate times: `rate_equation` behind a function boundary the compiler may not
inline, so it cannot hoist the call out of the timing loop."""
@noinline _testhelper_rate_barrier(m, concs, params) = rate_equation(m, concs, params)

"""
Test that `rate_equation` is non-allocating and fast for the given mechanism.
Must be a standalone function to avoid @testset closure boxing. Each call goes through
`_testhelper_rate_barrier`, so every one of the 10,000 calls per batch runs.
"""
function _testhelper_test_rate_equation_performance(m, params, concs)
    _testhelper_rate_barrier(m, concs, params) # warmup/compile
    allocs = @allocated _testhelper_rate_barrier(m, concs, params)
    # Minimum over several batches defeats the GC/scheduling inflation a
    # single mean suffers; accumulating into `acc` (and observing it via the
    # finite-result check) prevents the optimizer from eliding the calls.
    best = Inf
    acc = 0.0
    for _ in 1:5
        acc = 0.0
        t = @elapsed for _ in 1:10_000
            acc += _testhelper_rate_barrier(m, concs, params)
        end
        best = min(best, t / 10_000)
    end
    isfinite(acc) || error("rate_equation produced a non-finite result")
    return allocs, best
end
```

and rerun Step 7's command into `$S/phase45/logs/P5_T6_tr_barrier.log`: every gate passes,
the readings 1–2 ns above Step 7's. CLAUDE.md names this helper, so the PR body lists the
change for Denis (◆ P5-25) with the probe's counts. If `<h>` = 0, keep the helper and note
for the PR body that the package's bodies were not hoisted.

- [ ] **Step 9: Measure the compile budgets**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T6_budget.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_compile_budget.jl")'`
Expected: every testset passes, "a second process derives the same rate_equation bodies"
with 48 bodies; the `[ Info` lines read init trace-compile ≤ 200, first `rate_equation` under
6.0 s, Sig ≤ 250 instances and 100 KB, fit ≤ 250 instances and 300 KB, 0 show instances.
Record each reading beside Task 0's (`p5_base/figures.md`).

If the first `rate_equation` reading is over the budget, precompile the generator: append
to the four precompile lines at the end of `src/EnzymeRates.jl` (`_guard` and
`_faces_pass` on each mechanism type)

```julia
precompile(_rate_body, (Mechanism,))
precompile(_rate_body, (AllostericMechanism,))
```

and measure again. A reading still over a budget stops the task for Denis: no budget is
raised.

- [ ] **Step 10: Run the GC gate on the package**

Create `$S/phase45/gcscan5/derive.jl`:

```julia
# ABOUTME: R18 remark-scan workload of the derivation: laws, programs and bodies of a regulated
# ABOUTME: ping-pong population, with scan_f.jl's six merging loops, full collections inside.
using EnzymeRates
const ER = EnzymeRates
const MechU = Union{ER.Mechanism, ER.AllostericMechanism}
const LawU = Union{ER._Law, ER._AlloLaw}
const rxn = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    allosteric_regulators: R
    oligomeric_state: 2
end
seeds = MechU[m for m in ER.seed_mechanisms(rxn, Set([:R]), Set{Symbol}())]
grown = MechU[seeds[1:min(end, 60)]; ER.expand_mechanisms(seeds[1:min(end, 20)], rxn)]
pop = MechU[m for m in grown if ER._eq_complexity(m) <= 337][1:min(end, 300)]
println("julia ", VERSION, ": ", length(pop), " mechanisms, ",
        count(m -> m isa ER.AllostericMechanism, pop), " allosteric; mutable _Poly/_Law/_AlloLaw: ",
        ismutabletype(ER._Poly), "/", ismutabletype(ER._Law), "/", ismutabletype(ER._AlloLaw))
flush(stdout)
function derive_typed(ms::Vector{MechU})
    tot = zeros(Int, 3)
    for m in ms
        s = ER._solve(m)
        L = ER._law(m, s)
        ER._assert_support(m, L)
        nD = L isa ER._Law ? ER._nterms(L.D) : ER._nterms(L.active.D) + ER._nterms(L.inactive.D)
        tot .+= (1, nD, length(ER._body(ER._program(m, s, L)).args))
    end
    tot
end
function merge_vec(ls::Vector{LawU})
    t = 0
    for L in ls
        y = L isa ER._Law ? L.H : L.sites                 # Union{_Poly, Vector}
        t += y isa ER._Poly ? ER._nterms(y) : length(y)
    end
    t
end
function merge_nothing(ls::Vector{LawU})
    t = 0
    for (i, L) in enumerate(ls)
        z = isodd(i) ? L.H : nothing                      # Union{Nothing, _Poly}
        t += z === nothing ? 0 : ER._nterms(z)
    end
    t
end
function fresh_or_nothing(ms::Vector{MechU}, full::Bool)
    t = 0
    for (i, m) in enumerate(ms)
        z = isodd(i) ? ER._law(m) : nothing               # Union{Nothing, _Law, _AlloLaw}
        t += z === nothing ? 0 : z isa ER._Law ? ER._nterms(z.D) : ER._nterms(z.active.D)
        full && GC.gc(true)
    end
    t
end
function fresh_or_vec(ms::Vector{MechU}, full::Bool)
    t = 0
    for m in ms
        L = ER._law(m)
        y = L isa ER._Law ? L : L.sites                   # Union{_Law, Vector}
        t += y isa ER._Law ? ER._nterms(y.Q) : length(y)
        full && GC.gc(true)
    end
    t
end
function fresh_active_or_poly(ms::Vector{MechU}, full::Bool)
    t = 0
    for m in ms
        L = ER._law(m)
        y = L isa ER._AlloLaw ? L.active : L.H            # Union{_Law, _Poly}
        t += y isa ER._Law ? ER._nterms(y.D) : ER._nterms(y)
        full && GC.gc(true)
    end
    t
end
function fresh_poly_or_vec(ms::Vector{MechU}, full::Bool)
    t = 0
    for m in ms
        L = ER._law(m)
        y = L isa ER._Law ? L.H : L.sites                 # Union{_Poly, Vector}
        t += y isa ER._Poly ? ER._nterms(y) : length(y)
        full && GC.gc(true)
    end
    t
end
ref = nothing
for round in 1:3
    ls = LawU[ER._law(m) for m in pop]
    sig = (derive_typed(pop), merge_vec(ls), merge_nothing(ls),
           fresh_or_nothing(pop, false), fresh_or_vec(pop, false),
           fresh_active_or_poly(pop, false), fresh_poly_or_vec(pop, false),
           fresh_or_nothing(pop[1:60], true), fresh_or_vec(pop[1:60], true),
           fresh_active_or_poly(pop[1:60], true), fresh_poly_or_vec(pop[1:60], true))
    global ref = ref === nothing ? sig : ref
    GC.gc(true)
    println("round ", round, ": ", sig == ref ? "same" : "DIFFERENT", " ", sig,
            ", full gc ok")
    flush(stdout)
end
```

Create `$S/phase45/gcscan5/derive_control.jl`, the positive control (P5-7): derive.jl's six
merging loops over the same laws, whose polynomials ride in immutable one-reference
wrappers defined here, the layout that crashed (LAW §1.2):

```julia
# ABOUTME: R18 remark-scan positive control: derive.jl's six merging loops over the package's laws
# ABOUTME: copied into immutable one-field wrappers defined here, the layout that crashed (LAW §1.2).
using EnzymeRates
const ER = EnzymeRates
const MechU = Union{ER.Mechanism, ER.AllostericMechanism}
const rxn = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    allosteric_regulators: R
    oligomeric_state: 2
end
struct CPoly
    terms::Dict{ER._Mono, Int}
end
struct CLaw
    H::CPoly
    Q::CPoly
    D::CPoly
end
struct CAllo
    H::CPoly
    active::CLaw
    sites::Vector{Tuple{CPoly, CPoly, Int}}
end
const WrapU = Union{CLaw, CAllo}
cpoly(p::ER._Poly) = CPoly(p.terms)
claw(L::ER._Law) = CLaw(cpoly(L.H), cpoly(L.Q), cpoly(L.D))
"""`m`'s law, freshly derived, in the wrappers: a `CLaw` or a `CAllo`."""
@noinline function wrapped(m)
    L = ER._law(m)
    L isa ER._Law && return claw(L)
    CAllo(cpoly(L.H), claw(L.active),
          Tuple{CPoly, CPoly, Int}[(cpoly(a), cpoly(b), k) for (a, b, k) in L.sites])
end
nt(p::CPoly) = length(p.terms)
seeds = MechU[m for m in ER.seed_mechanisms(rxn, Set([:R]), Set{Symbol}())]
grown = MechU[seeds[1:min(end, 60)]; ER.expand_mechanisms(seeds[1:min(end, 20)], rxn)]
pop = MechU[m for m in grown if ER._eq_complexity(m) <= 337][1:min(end, 300)]
println("julia ", VERSION, ": ", length(pop), " mechanisms; immutable CPoly/CLaw/CAllo: ",
        !ismutabletype(CPoly), "/", !ismutabletype(CLaw), "/", !ismutabletype(CAllo))
flush(stdout)
function merge_vec(ls::Vector{WrapU})
    t = 0
    for L in ls
        y = L isa CLaw ? L.H : L.sites                    # Union{CPoly, Vector}
        t += y isa CPoly ? nt(y) : length(y)
    end
    t
end
function merge_nothing(ls::Vector{WrapU})
    t = 0
    for (i, L) in enumerate(ls)
        z = isodd(i) ? L.H : nothing                      # Union{Nothing, CPoly}
        t += z === nothing ? 0 : nt(z)
    end
    t
end
function fresh_or_nothing(ms::Vector{MechU}, full::Bool)
    t = 0
    for (i, m) in enumerate(ms)
        z = isodd(i) ? wrapped(m) : nothing               # Union{Nothing, CLaw, CAllo}
        t += z === nothing ? 0 : z isa CLaw ? nt(z.D) : nt(z.active.D)
        full && GC.gc(true)
    end
    t
end
function fresh_or_vec(ms::Vector{MechU}, full::Bool)
    t = 0
    for m in ms
        L = wrapped(m)
        y = L isa CLaw ? L : L.sites                      # Union{CLaw, Vector}
        t += y isa CLaw ? nt(y.Q) : length(y)
        full && GC.gc(true)
    end
    t
end
function fresh_active_or_poly(ms::Vector{MechU}, full::Bool)
    t = 0
    for m in ms
        L = wrapped(m)
        y = L isa CAllo ? L.active : L.H                  # Union{CLaw, CPoly}
        t += y isa CLaw ? nt(y.D) : nt(y)
        full && GC.gc(true)
    end
    t
end
function fresh_poly_or_vec(ms::Vector{MechU}, full::Bool)
    t = 0
    for m in ms
        L = wrapped(m)
        y = L isa CLaw ? L.H : L.sites                    # Union{CPoly, Vector}
        t += y isa CPoly ? nt(y) : length(y)
        full && GC.gc(true)
    end
    t
end
ref = nothing
for round in 1:3
    ls = WrapU[wrapped(m) for m in pop]
    sig = (merge_vec(ls), merge_nothing(ls),
           fresh_or_nothing(pop, false), fresh_or_vec(pop, false),
           fresh_active_or_poly(pop, false), fresh_poly_or_vec(pop, false),
           fresh_or_nothing(pop[1:60], true), fresh_or_vec(pop[1:60], true),
           fresh_active_or_poly(pop[1:60], true), fresh_poly_or_vec(pop[1:60], true))
    global ref = ref === nothing ? sig : ref
    GC.gc(true)
    println("round ", round, ": ", sig == ref ? "same" : "DIFFERENT", " ", sig,
            ", full gc ok")
    flush(stdout)
end
```

The control may abort in a collection, as research variant C did; its remarks are in the
log by then, since the loops compile before they run.

Create `$S/phase45/gcscan5/identify.jl`:

```julia
# ABOUTME: R18 remark-scan workload: a small identify_rate_equation search on a uni-uni
# ABOUTME: reaction with a regulator, the suite's smoke settings, three rounds.
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
vals = ((1.0 + 0.1i for i in eachindex(fp))..., 2.0, 1.0)   # not `values`, a Base name
params = NamedTuple{(fp..., :Keq, :E_total)}(vals)
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

Make `scan_law.jl` (Task 2, which names its workload in `ARGS[1]`) run these files: replace
its usage line,

```julia
# Usage: julia --project=<env> --pkgimages=no scan_law.jl poly [control]
```

with

```julia
# Usage: julia --project=<env> --pkgimages=no scan_law.jl <workload>, the workload one of
# `poly`, `poly control`, `derive`, `derive_control` and `identify` (the last three are
# the files of those names beside this one).
ARGS[1] in ("derive", "derive_control", "identify") &&
    (include(joinpath(@__DIR__, ARGS[1] * ".jl")); exit())
```

Run (background), one Julia at a time:

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; G=$S/phase45/gcscan5; LG=$S/phase45/logs
export JULIA_LLVM_ARGS="-pass-remarks=julia-licm -pass-remarks-analysis=julia-licm"
for w in derive derive_control identify; do
  n=$w; [ $w = derive_control ] && n=control
  $S/phase3/jl.sh $LG/P5_T6_scan112_$n.log 1200 --project=$G/env112 --startup-file=no --pkgimages=no --heap-size-hint=2500M $G/scan_law.jl $w
  $S/phase3/jl.sh $LG/P5_T6_scan113_$n.log 1200 +1.13 --project=$G/env113 --startup-file=no --pkgimages=no --heap-size-hint=2500M $G/scan_law.jl $w
done
$S/phase3/jl.sh $LG/P5_T6_repro112.log 300 --startup-file=no --pkgimages=no $S/phase3/gcprobe/repro_upstream.jl
$S/phase3/jl.sh $LG/P5_T6_repro113.log 300 +1.13 --startup-file=no --pkgimages=no $S/phase3/gcprobe/repro_upstream.jl
unset JULIA_LLVM_ARGS
python3 -I $S/phase3/gcprobe/hoistscan.py $LG/P5_T6_scan11*.log $LG/P5_T6_repro11*.log > $LG/P5_T6_hoistscan.log 2>&1
cat $LG/P5_T6_hoistscan.log
```

Expected: `derive` and `identify` on both versions print three rounds `same … full gc ok`
and report `hoisted with a bits store into a reference field: 0`; `derive_control` and
`repro_upstream.jl` report at least one `DANGER` on both versions. A package DANGER line
stops the task: find the loop, restore a concrete type, rescan. A control without DANGER
means the scan cannot see the hazard: stop and report it.

- [ ] **Step 11: Compare the census with Task 0's**

Run (background), with the arguments Task 0 gave `census/P5_base` (the rename map
included), so both corpora are keyed alike, and compare.jl from `census/env` as phase 4
ran it:
`S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; $S/phase3/jl.sh $S/phase45/logs/P5_T6_census.log 1800 --project=$S/phase3/census/env_law --heap-size-hint=2500M $S/phase3/census/corpus.jl $S/phase45/law_wt $S/phase3/census/P5_T6 $S/phase3/census/rename_phase3.toml && $S/phase3/jl.sh $S/phase45/logs/P5_T6_census_cmp.log 600 --project=$S/phase3/census/env --heap-size-hint=2500M $S/phase3/census/compare.jl $S/phase3/census/P5_base $S/phase3/census/P5_T6 $S/phase3/census/rename_phase3.toml`, polled.
Expected: `RESULT: PASS` in `P5_T6_census_cmp.log`: the same child sets, `rate_equation`
values within compare.jl's `TOL` (1e-9), and fitted names and printed strings equal, since
the printer is unchanged.

`P5_base` was built before Task 4's summed cap, so a law whose two conformations sum past
5,000 is `ok` there and an error here, which compare.jl reports as a law difference. The
gate accepts a `FAIL` only when every line of its list is a `<group> laws: <n>
differences, 0 strings the rename cannot map exactly` line and the counts `<n>` add up to
the lines this prints:

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; C=$S/phase3/census; awk -F'\t' 'NR == FNR { st[$1] = $5; next } st[$1] == "ok" && $5 == "error" { print $1 "\t" substr($8, 1, 80) }' $C/P5_base/laws.tsv $C/P5_T6/laws.tsv | tee $S/phase45/logs/P5_T6_census_refused.log
```

each line's text starting `Rate equation for this mechanism has more than 5000 polynomial
terms` (RateEquationTooLargeError's message). Name those laws (their keys) in the commit;
Task 12's census meets the same list, and the PR body carries it. Any other difference is
a bug.

- [ ] **Step 12: Run the full suite and R52's figures**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T6_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`, polled; then
the ter-ter timing with Task 0's script, so the figures compare:
`S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T6_terter.log 900 --project=. --heap-size-hint=2500M $S/phase45/scripts/P4_terter.jl`
and the TME focused run:
`S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T6_tme.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_mechanism_enumeration.jl")'`.
Expected: `exit 0`, every test passing, pristine, Aqua and JET green, about 15 s over Task
5's time (the subprocess about 10 s, the BigFloat bodies about 5 s); `ter-ter init: 251107 mechanisms in <t> s; the guard rejects 0`,
well under 150; the TME log's `exit 0 elapsed <s>s`. Report the three figures (full-suite
`maxrss_kb`, ter-ter seconds, TME seconds) beside Task 0's (`p5_base/figures.md`).

- [ ] **Step 13: Commit**

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt
git add src/codegen.jl src/EnzymeRates.jl src/rate_eq_derivation.jl \
    src/sym_poly_for_rate_eq_derivation.jl docs/src/developer.md \
    test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/test_rate_eq_derivation.jl test/test_compile_budget.jl
git commit -m "$(cat <<'EOF'
Generate rate_equation from the derived law

The body binds the dependents first, in solve-column order, then each
segment's weight and tree sums, Q, H and the MWC terms, all in greedy Horner
form, and reads every input by name into a positional local; roots are sqrt
or x^(one(x)/q). The Expr-shape tests now walk every local, and BigFloat
inputs run through the body of the allosteric specs and four plain ones.
_rate_body has one method per mechanism type, so JET sees each law's type.
Readings: specs max <t> ns (slowest <spec 1>, <spec 2>); gate laws
<ns-per-term list> ns per term, all under 0.8 ns x D, 0 allocations.
Hoisting probe: <h> of 42 specs under 0.5 ns inline<, so the gate's helper
times through a @noinline barrier>. Compile budgets: <readings>. GC gate: 0
hoisted on 1.12.5 and 1.13.1, controls positive. Census: <PASS, or FAIL
only on the n laws the summed cap refuses, with their keys>. R52: maxrss
<kb> KB, ter-ter init <s> s, TME <s> s (Task 0: <figures>).

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
)"
```

---

### Task 7: The key (P5-18, P5-19, P5-21)

The key labels a law canonically: each unit to lowest terms, fitted constants coloured by
refinement on the term–constant incidence and individualized while colours tie,
concentrations labelled by name, and the text's full SHA-256. At `_KEY_LEAF_CAP` labellings
the search stops and the law gets its mechanism's own key. identify does not read the key
until Task 8; here the dedup-key tests, the 239-class pin and the population tests move to
it.

**Files:**
- Create: `src/key.jl`; modify `src/EnzymeRates.jl` (include `key.jl` before `codegen.jl`;
  `using SHA: sha256` moves here from `src/identify_rate_equation.jl`)
- Modify: `src/identify_rate_equation.jl` (delete its `using SHA: sha256`)
- Create: `test/test_key.jl`; modify `test/runtests.jl` (include `test_key.jl` after
  `test_law.jl`)
- Modify: the fixtures file (`_testhelper_same_law_under_map`, `_testhelper_renamed_law`)
- Test: `test/test_oracles.jl` (Q2, Q3), `test/test_identify_rate_equation.jl` (Q4),
  `test/test_mechanism_enumeration.jl` ("Rate-equation dedup key"),
  `test/test_compile_budget.jl` (keys in the subprocess)
- Create (scratch): `$S/phase45/probes/P5_T7_keycost.jl`
- Logs: `$S/phase45/logs/P5_T7_{key_red,oracles_red,budget_red,key,pins112,pins113,oracles,tid,tme,budget,keycost,full}.log`

**Interfaces:**
- Consumes: Task 2's `_Vars`, `_Mono`, `_Exp`, `_Poly`, `_lowest_shift`, `_mmul`, `_scale`,
  `_term`, `_monomial`; Tasks 3–4's `_Law`, `_AlloLaw`, `_law`; `SHA.sha256`; tests: Task 3's
  `_testhelper_evaluate_law`, test_oracles.jl's `_testhelper_written_backwards`,
  `_testhelper_content_text`, `_testhelper_population_sample`,
  `_testhelper_population_one_round`, `_testhelper_inhibitor_rxn`'s `children`.
- Produces:
  ```julia
  struct _Canonical
      order::Vector{Symbol}  # fitted constants in canonical order, constants no term
                             #   holds last by name; L is fixed and not listed
      text::String           # header|P=…|absent=…, then sorted term lines
      key::String            # bytes2hex(sha256(text)), 64 hex digits
      leaves::Int            # the labellings the search compared (≤ _KEY_LEAF_CAP)
      capped::Bool           # the cap was reached: key is the per-mechanism key
  end
  const _KEY_LEAF_CAP = 100
  _key_parts(L::_Law)::Tuple{Vector{_Poly}, Vector{Vector{Int}}, String}     # [H, Q, D]
  _key_parts(A::_AlloLaw)::Tuple{Vector{_Poly}, Vector{Vector{Int}}, String}
      # [H, QA, DA, fA, QI, DI, fI, (W_A, W_I)…]; units [[1], [2,3,4], [5,6,7], one per
      #   site]; header "mwc|n=…|sites=…"
  _ranks(xs::AbstractVector)::Vector{Int}
  _refine(pcol::Vector{Int}, tcol::Vector{Int}, inc::Vector{Vector{Tuple{Int, _Exp}}},
          occ::Vector{Vector{Tuple{Int, _Exp}}})::Tuple{Vector{Int}, Vector{Int}}
  _encode(pcol::Vector{Int}, inc::Vector{Vector{Tuple{Int, _Exp}}},
          labels::Vector{String})::String
  mutable struct _KeySearch; text::String; order::Vector{Int}; leaves::Int; end
  _search!(st::_KeySearch, pcol::Vector{Int}, tcol::Vector{Int}, inc, occ,
           labels::Vector{String})::_KeySearch
  _canonical(L::Union{_Law, _AlloLaw}, mechanism_id::AbstractString)::_Canonical
  # tests (fixtures file)
  _testhelper_same_law_under_map(L1, L2, order1::Vector{Symbol}, order2::Vector{Symbol};
                                 rng, swap::Bool = false)::Bool
  _testhelper_renamed_law(L, rename::Dict{Symbol, Symbol})   # _Law or _AlloLaw
  _testhelper_renamed_vars(vars::EnzymeRates._Vars,
                           rename::Dict{Symbol, Symbol})::EnzymeRates._Vars
  ```
  `_canonical` ports `SLC/key.jl` `canonical` with KC §6's change (`klib.jl`
  `canonical_named`): concentrations are labelled `c_<name>`, each term's fixed factors
  sorted. Fixed variables are concentrations, `Keq` and `L`. The search keeps its best
  labelling in a `_KeySearch`, so no closure captures a variable it reassigns. At
  `_KEY_LEAF_CAP` leaves `text = "mechanism|" * mechanism_id`, `key` its SHA-256, `order`
  the variable order and `capped = true`.

**Gate:**
- `test_key.jl` is green; the five pinned keys are equal on 1.12.5 and 1.13.1.
- Q2, Q3 and Q4 pass; the dedup-key testsets keep their three verdicts.
- The subprocess keys equal the parent's on the 42 specs and six gate laws.
- The key-cost probe on the 239 bi-bi init laws and 200 stratified ter-ter init laws
  (V×τ ≤ 5,000) finds 0 leaf-cap hits; medians and maxima go to the commit and the PR body
  (slice: 351 ms at most).
- The full suite is green and pristine; +4 s.

- [ ] **Step 1: Write the shared test helpers**

Append to the fixtures file:

```julia
"""
    _testhelper_same_law_under_map(L1, L2, order1, order2; rng, swap = false) -> Bool

Whether the laws `L1` and `L2` (`_Law` or `_AlloLaw`) agree to 1e-20 relative at 3
random points when constant `order2[i]` of `L2` takes the value of `order1[i]` of `L1`, the
canonical positions of `_canonical`, and each concentration, `Keq` and `L` take one value
in both. With `swap`, `order2`'s first two constants trade places, a map that must fail
when the law tells them apart.
"""
function _testhelper_same_law_under_map(L1, L2, order1::Vector{Symbol},
                                        order2::Vector{Symbol}; rng, swap::Bool = false)
    length(order1) == length(order2) || return false
    o2 = copy(order2)
    swap && length(o2) >= 2 && ((o2[1], o2[2]) = (o2[2], o2[1]))
    concs = union(L1.vars.names[L1.vars.conc], L2.vars.names[L2.vars.conc])
    setprecision(BigFloat, 256) do
        all(1:3) do _
            env1 = Dict{Symbol, BigFloat}(x => exp(4 * big(rand(rng)) - 2) for x in order1)
            env2 = Dict{Symbol, BigFloat}(o2[i] => env1[order1[i]] for i in eachindex(o2))
            for x in [concs; :Keq; :L]
                env1[x] = env2[x] = exp(3 * big(rand(rng)) - big(1.5))
            end
            v1 = _testhelper_evaluate_law(L1, env1)
            v2 = _testhelper_evaluate_law(L2, env2)
            abs(v1 - v2) <= big(10)^-20 * abs(v1)
        end
    end
end

"""The variable table `vars` with each name `x` renamed `get(rename, x, x)`."""
function _testhelper_renamed_vars(vars::EnzymeRates._Vars, rename::Dict{Symbol, Symbol})
    names = Symbol[get(rename, x, x) for x in vars.names]
    EnzymeRates._Vars(names, copy(vars.conc), Dict(x => i for (i, x) in enumerate(names)))
end

"""The law `L` with its variables renamed by `rename` (`_testhelper_renamed_vars`) and its
polynomials unchanged: a renaming of its constants when `rename` permutes them."""
_testhelper_renamed_law(L::EnzymeRates._Law, rename::Dict{Symbol, Symbol}) =
    EnzymeRates._Law(_testhelper_renamed_vars(L.vars, rename), L.H, L.Q, L.D, L.sigma,
                     L.trees, L.shift)
_testhelper_renamed_law(A::EnzymeRates._AlloLaw, rename::Dict{Symbol, Symbol}) =
    EnzymeRates._AlloLaw(_testhelper_renamed_vars(A.vars, rename), A.H,
                         _testhelper_renamed_law(A.active, rename),
                         _testhelper_renamed_law(A.inactive, rename), A.fA, A.fI, A.sites,
                         A.n)
```

- [ ] **Step 2: Write the failing key tests**

Create `test/test_key.jl`:

```julia
# ABOUTME: Tests of the dedup key: canonical labelling of a law's fitted constants, lowest
# ABOUTME: terms, the leaf cap, which metabolite binds, and the pinned keys of five specs.

# ── Toy laws over named constants ───────────────────────────────────────────

"""The polynomial Σ c·Π x^e of `terms`, each `c => [x => e, …]`, its variables looked up
by name in `vars`."""
_testhelper_poly(vars::EnzymeRates._Vars, terms) =
    foldl(+, (EnzymeRates._term(EnzymeRates._monomial(
                  (vars.index[x] => e for (x, e) in xs)...), c) for (c, xs) in terms);
          init = zero(EnzymeRates._Poly))

"""A plain law over the concentrations `concs`, the constants `constants` and `:Keq`,
with `H`, `Q` and `D` written as `_testhelper_poly` terms; the key reads no other field."""
function _testhelper_toy_law(concs, constants, H, Q, D)
    vars = EnzymeRates._Vars()
    foreach(x -> EnzymeRates._var!(vars, x; conc = true), concs)
    foreach(x -> EnzymeRates._var!(vars, x), constants)
    EnzymeRates._var!(vars, :Keq)
    p(ts) = _testhelper_poly(vars, ts)
    EnzymeRates._Law(vars, p(H), p(Q), p(D), [p(D)], [one(EnzymeRates._Poly)],
                     EnzymeRates._Mono())
end

"""Reversible Michaelis–Menten over the constants `k1`, `k2`, `k3`, `k4` (and `extra`,
which no term holds): `H = A − P/Keq`, `Q = k1·k2`, `D = k3 + k2 + k1·A + k4·P`."""
_testhelper_mm_law(k1, k2, k3, k4; extra = Symbol[]) = _testhelper_toy_law(
    [:A, :P], [k1, k2, k3, k4, extra...],
    [1 => [:A => 1], -1 => [:P => 1, :Keq => -1]], [1 => [k1 => 1, k2 => 1]],
    [1 => [k3 => 1], 1 => [k2 => 1], 1 => [k1 => 1, :A => 1], 1 => [k4 => 1, :P => 1]])

"""A law whose constants come in two symmetric pairs: `a`, `b` only as `a·b` and `c`, `d`
only as `c·d`."""
_testhelper_pairs_law(a, b, c, d) = _testhelper_toy_law(
    [:A, :P], [a, b, c, d], [1 => [:A => 1], -1 => [:P => 1, :Keq => -1]],
    [1 => Pair{Symbol, Int}[]],
    [1 => Pair{Symbol, Int}[], 1 => [a => 1, b => 1, :A => 1],
     1 => [c => 1, d => 1, :P => 1]])

"""A law whose denominator holds one term `x·y` per edge `(x, y)` of a graph on the six
constants p1…p6: every constant has degree 2, so refinement alone cannot split them."""
_testhelper_graph_law(edges) = _testhelper_toy_law(
    [:A, :P], [:p1, :p2, :p3, :p4, :p5, :p6],
    [1 => [:A => 1], -1 => [:P => 1, :Keq => -1]], [1 => Pair{Symbol, Int}[]],
    [1 => [x => 1, y => 1] for (x, y) in edges])

"""A law over the concentration S and nine constants a1…a9 that multiply S alike in
`D = 1 + Σ S·a_i`: 9! labellings tie."""
_testhelper_nine_law() = _testhelper_toy_law(
    [:S], [Symbol(:a, i) for i in 1:9], [1 => [:S => 1]], [1 => Pair{Symbol, Int}[]],
    [1 => Pair{Symbol, Int}[]; [1 => [Symbol(:a, i) => 1, :S => 1] for i in 1:9]])

const _testhelper_cycle6 = [(:p1, :p2), (:p2, :p3), (:p3, :p4), (:p4, :p5), (:p5, :p6),
                            (:p6, :p1)]

@testset "renaming fitted constants keeps the key and gives the induced order" begin
    a = EnzymeRates._canonical(_testhelper_mm_law(:a, :b, :c, :d), "toy")
    b = EnzymeRates._canonical(_testhelper_mm_law(:w, :x, :y, :z), "toy")
    c = EnzymeRates._canonical(_testhelper_mm_law(:z, :y, :x, :w), "toy")
    @test a.key == b.key == c.key
    @test length(a.key) == 64 && all(in("0123456789abcdef"), a.key)
    @test Dict(zip(a.order, c.order)) == Dict(:a => :z, :b => :y, :c => :x, :d => :w)
    @test a.leaves == 1 && !a.capped
    # Random renamings of three toys, each a permutation of the constant names.
    rng = Random.MersenneTwister(94)
    for L in (_testhelper_mm_law(:a, :b, :c, :d), _testhelper_pairs_law(:a, :b, :c, :d),
              _testhelper_graph_law(_testhelper_cycle6))
        names = filter(!=(:Keq), L.vars.names[.!L.vars.conc])
        cL = EnzymeRates._canonical(L, "toy")
        for _ in 1:3
            R = _testhelper_renamed_law(L, Dict(zip(names, shuffle(rng, names))))
            cR = EnzymeRates._canonical(R, "toy")
            @test cR.key == cL.key
            @test _testhelper_same_law_under_map(L, R, cL.order, cR.order; rng)
        end
    end
    # Swapping two constants' roles in the denominator changes the function and the key.
    swapped = _testhelper_toy_law([:A, :P], [:a, :b, :c, :d],
        [1 => [:A => 1], -1 => [:P => 1, :Keq => -1]], [1 => [:a => 1, :b => 1]],
        [1 => [:c => 1], 1 => [:b => 1], 1 => [:b => 1, :A => 1], 1 => [:d => 1, :P => 1]])
    @test EnzymeRates._canonical(swapped, "toy").key != a.key
end

@testset "the key brings each unit to lowest terms" begin
    # Q and D times 3·a⁻²·A, and H times 2: each unit's lowest-terms form, and so the key,
    # stay the same.
    base = _testhelper_mm_law(:a, :b, :c, :d)
    v = base.vars
    m = EnzymeRates._monomial(v.index[:a] => -2, v.index[:A] => 1)
    scaled = EnzymeRates._Law(v, EnzymeRates._scale(base.H, 2),
                              EnzymeRates._scale(base.Q * m, 3),
                              EnzymeRates._scale(base.D * m, 3), base.sigma, base.trees,
                              base.shift)
    @test EnzymeRates._canonical(scaled, "toy").key ==
          EnzymeRates._canonical(base, "toy").key
end

@testset "symmetric constants take 4 leaves and one key" begin
    # a and b tie, and c and d: two branches for each pair, four labellings in all.
    s = _testhelper_pairs_law(:a, :b, :c, :d)
    cs = EnzymeRates._canonical(s, "toy")
    @test cs.leaves == 4 && !cs.capped
    s2 = _testhelper_renamed_law(s, Dict(:a => :d, :b => :a, :c => :c, :d => :b))
    c2 = EnzymeRates._canonical(s2, "toy")
    @test c2.key == cs.key
    @test _testhelper_same_law_under_map(s, s2, cs.order, c2.order;
                                         rng = Random.MersenneTwister(95))
end

@testset "a 6-cycle and two triangles get different keys" begin
    # Both graphs are 2-regular, so refinement alone ties every constant; individualizing
    # tells them apart, and a relabelled 6-cycle keeps its key.
    cycle = _testhelper_graph_law(_testhelper_cycle6)
    triangles = _testhelper_graph_law([(:p1, :p2), (:p2, :p3), (:p3, :p1), (:p4, :p5),
                                       (:p5, :p6), (:p6, :p4)])
    k6 = EnzymeRates._canonical(cycle, "toy")
    k33 = EnzymeRates._canonical(triangles, "toy")
    @test k6.key != k33.key
    @test k6.leaves > 1 && k33.leaves > 1 && !k6.capped && !k33.capped
    names = [:p1, :p2, :p3, :p4, :p5, :p6]
    rename = Dict(zip(names, shuffle(Random.MersenneTwister(96), names)))
    @test EnzymeRates._canonical(_testhelper_renamed_law(cycle, rename), "toy").key ==
          k6.key
end

@testset "an extra declared constant changes the key" begin
    plain = EnzymeRates._canonical(_testhelper_mm_law(:a, :b, :c, :d), "toy")
    extra = EnzymeRates._canonical(_testhelper_mm_law(:a, :b, :c, :d; extra = [:e]), "toy")
    @test extra.key != plain.key
    @test extra.order[end] === :e
    @test occursin("|P=5|absent=1", extra.text) && occursin("|P=4|absent=0", plain.text)
end

@testset "which metabolite binds changes the key" begin
    # Concentrations are labelled by name, so models that differ only in which inhibitor
    # or regulator binds, or where an inhibitor binds, get different keys.
    key(m) = EnzymeRates._canonical(EnzymeRates._law(_testhelper_lift(m)), "attack").key
    uni_I = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
            E + I ⇌ E(I)
        end
    end
    uni_J = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: J
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
            E + J ⇌ E(J)
        end
    end
    unc_I = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
            E(S) + I ⇌ E(S, I)
        end
    end
    allo_R = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)       :: NonequalAI
            E(S) <--> E + P    :: NonequalAI
        end
        regulatory_site(1): R::OnlyA
    end
    allo_T = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: T
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)       :: NonequalAI
            E(S) <--> E + P    :: NonequalAI
        end
        regulatory_site(1): T::OnlyA
    end
    @test key(uni_I) != key(uni_J)
    @test key(uni_I) != key(unc_I)
    @test key(allo_R) != key(allo_T)
end

@testset "at 100 leaves a law gets its mechanism's key" begin
    # Nine constants that tie: the search stops at the cap, and the law takes the key of
    # its mechanism's id, so twins of it are fitted apart.
    id = "1234567e89012345"
    c = EnzymeRates._canonical(_testhelper_nine_law(), id)
    @test c.capped && c.leaves == EnzymeRates._KEY_LEAF_CAP == 100
    @test c.text == "mechanism|" * id
    @test c.key == bytes2hex(EnzymeRates.sha256("mechanism|" * id))
    @test c.order == [Symbol(:a, i) for i in 1:9]
end
```

In `test/runtests.jl`, add `    include("test_key.jl")` directly after
`    include("test_law.jl")`.

In `test/test_oracles.jl`, add this helper at top level, just before the testset "Laws over
a regenerated population":

```julia
"""The regroupings of `m` that keep its model: each split of one group by a context
bipartition (`_context_bipartitions`) that the guard accepts and that keeps the
independent count."""
function _testhelper_same_model_regroupings(m)
    n = EnzymeRates._independent_param_count(m)
    cm = EnzymeRates.catalytic(m)
    kept = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[]
    for g in 1:EnzymeRates.ngroups(cm), bp in EnzymeRates._context_bipartitions(cm, g)
        child = EnzymeRates._apply_bipartitions(m, [(g, bp)])
        EnzymeRates._is_well_posed(child) &&
            EnzymeRates._independent_param_count(child) == n && push!(kept, child)
    end
    kept
end

# The same-model regroupings of each mechanism of `_testhelper_population_sample`, in its
# order, computed once for the two testsets that read them (the guard runs on every child).
const _testhelper_sample_regroupings =
    map(_testhelper_same_model_regroupings, _testhelper_population_sample)
```

In "a same-model regrouping keeps every fitted name", replace

```julia
        for m in _testhelper_population_sample
            names = Set(EnzymeRates._dependent_param_exprs(m)[2])
            n = EnzymeRates._independent_param_count(m)
            cm = EnzymeRates.catalytic(m)
            for g in 1:EnzymeRates.ngroups(cm),
                bp in EnzymeRates._context_bipartitions(cm, g)
                child = EnzymeRates._apply_bipartitions(m, [(g, bp)])
                EnzymeRates._is_well_posed(child) || continue
                EnzymeRates._independent_param_count(child) == n || continue
                m isa EnzymeRates.Mechanism ? (n_plain += 1) : (n_allo += 1)
                @test Set(EnzymeRates._dependent_param_exprs(child)[2]) == names
            end
        end
```

with

```julia
        for (m, regrouped) in zip(_testhelper_population_sample,
                                  _testhelper_sample_regroupings)
            names = Set(EnzymeRates._dependent_param_exprs(m)[2])
            for child in regrouped
                m isa EnzymeRates.Mechanism ? (n_plain += 1) : (n_allo += 1)
                @test Set(EnzymeRates._dependent_param_exprs(child)[2]) == names
            end
        end
```

After the testset "the guard's verdict survives reordering, reversal, regrouping and
renames" (it follows `_testhelper_written_backwards`, which Q2 reads), add:

```julia
@testset "rewrites, regroupings and constant renames keep every sampled key" begin
    # Q2 (spec section 14): writing a mechanism backwards, regrouping it without changing
    # its model, and renaming its fitted constants change no law and no key, and the
    # canonical order maps one labelling onto the other. A same-model regrouping keeps
    # every fitted name too ("a same-model regrouping keeps every fitted name").
    rng = Random.MersenneTwister(612)
    n_regrouped = 0
    for (m, regrouped) in zip(_testhelper_population_sample, _testhelper_sample_regroupings)
        L = EnzymeRates._law(m)
        c = EnzymeRates._canonical(L, "q2")
        @test !c.capped
        back = EnzymeRates._law(_testhelper_written_backwards(m))
        @test EnzymeRates._canonical(back, "q2").key == c.key
        for child in regrouped
            n_regrouped += 1
            Lc = EnzymeRates._law(child)
            cc = EnzymeRates._canonical(Lc, "q2")
            @test cc.key == c.key
            @test _testhelper_same_law_under_map(L, Lc, c.order, cc.order; rng)
        end
        constants = filter(x -> x ∉ (:Keq, :L), L.vars.names[.!L.vars.conc])
        rename = Dict(zip(constants, shuffle(rng, constants)))
        R = _testhelper_renamed_law(L, rename)
        r = EnzymeRates._canonical(R, "q2")
        @test r.key == c.key
        @test Set(r.order) == Set(rename[x] for x in c.order)
        @test _testhelper_same_law_under_map(L, R, c.order, r.order; rng)
    end
    @test n_regrouped >= 70
end
```

In `test/test_oracles.jl`, just before the line `let children = [c for m in …` (at top
level: a function defined in the `let` block would be local to it), add

```julia
"""60 children of `children` that bind I, spaced evenly through their content order
(`_testhelper_content_text`); when fewer than 60 bind I, some repeat."""
function _testhelper_i_binding_sample(children)
    binds_i = unique!([c for c in children
                       if any(x -> x.name === :I, EnzymeRates.catalytic(c).table)])
    sorted = binds_i[sortperm(_testhelper_content_text.(binds_i))]
    sorted[round.(Int, range(1, length(sorted); length = 60))]
end
```

In "the guard equals the RE-speed oracle on inhibitor children's modifier faces",
replace

```julia
    binds_i = unique!([c for c in children
                       if any(x -> x.name === :I, EnzymeRates.catalytic(c).table)])
    @test length(binds_i) >= 60
    sorted = binds_i[sortperm(_testhelper_content_text.(binds_i))]
    for c in sorted[round.(Int, range(1, length(sorted); length = 60))]
```

with

```julia
    sample = _testhelper_i_binding_sample(children)
    # At least 60 children bind I, so the 60 spaced picks are distinct.
    @test allunique(sample)
    for c in sample
```

In the comment above the `let` line, replace "The next two testsets read one population"
with "The next three testsets read one population", since Q3 joins them; and add after
that testset, inside the `let` block:

```julia
@testset "laws that share a key are one law under the canonical map" begin
    # Q3 (spec section 14): the sample, the 60 I-binding inhibitor children, and every
    # tagging child of ten catalytic graphs of the one-round population, whose siblings
    # often share a law (tag sets that differ only on groups the inactive conformation
    # never reaches, or only in names). Members of one key evaluate alike with constants
    # mapped by canonical position, and swapping two of them breaks the map.
    by_graph = Dict{EnzymeRates.Mechanism, Vector{EnzymeRates.AllostericMechanism}}()
    for m in _testhelper_population_one_round
        m isa EnzymeRates.AllostericMechanism || continue
        push!(get!(by_graph, m.catalytic, EnzymeRates.AllostericMechanism[]), m)
    end
    graphs = sort!(collect(keys(by_graph)); by = _testhelper_content_text)
    siblings = reduce(vcat, [by_graph[g] for g in
                             graphs[round.(Int, range(1, length(graphs); length = 10))]])
    # A mechanism listed twice would form a merged group with itself.
    pool = unique!([_testhelper_population_sample; _testhelper_i_binding_sample(children);
                    siblings])
    keyed = Dict{String, Vector{Any}}()
    for m in pool
        L = EnzymeRates._law(m)
        c = EnzymeRates._canonical(L, EnzymeRates._mechanism_id(sprint(show, m)))
        @test !c.capped
        push!(get!(keyed, c.key, Any[]), (L, c))
    end
    merged = [v for v in values(keyed) if length(v) > 1]
    @test !isempty(merged)
    rng = Random.MersenneTwister(613)
    for group in merged
        (L1, c1) = group[1]
        for (L2, c2) in group[2:end]
            @test _testhelper_same_law_under_map(L1, L2, c1.order, c2.order; rng)
        end
        length(c1.order) >= 2 || continue
        (L2, c2) = group[2]
        @test !_testhelper_same_law_under_map(L1, L2, c1.order, c2.order; rng, swap = true)
    end
end
```

In `test/test_identify_rate_equation.jl`, replace the testset "rate-eq dedup-key partition
stability" by:

```julia
@testset "the key partitions the bi-bi init laws into 239 classes" begin
    # Q4. The 239 bi-bi init mechanisms (55 seeds and their 184 merged and Theorell–Chance
    # variants) have 239 laws no renaming of fitted constants makes equal, so the key
    # gives 239 classes, as the text key did on these laws. On a fixed sample the same law
    # gives the same key twice.
    reaction = @enzyme_reaction(begin
        substrates: A[C], B[N]
        products:   P[C], Q[N]
    end)
    buckets = Dict{String, Vector{Int}}()
    for (i, m) in enumerate(EnzymeRates.init_mechanisms(reaction))
        id = EnzymeRates._mechanism_id(sprint(show, m))
        c = EnzymeRates._canonical(EnzymeRates._law(m), id)
        @test !c.capped
        push!(get!(buckets, c.key, Int[]), i)
        i % 25 == 1 && @test EnzymeRates._canonical(EnzymeRates._law(m), id).key == c.key
    end
    @test length(buckets) == 239
end
```

In `test/test_mechanism_enumeration.jl`, testset "Rate-equation dedup key": in "Distinct
mechanisms produce distinct dedup keys" replace

```julia
        @test EnzymeRates._rate_eq_dedup_key(rate_equation_string(m_ordered)) !=
              EnzymeRates._rate_eq_dedup_key(rate_equation_string(m_random))
```

with

```julia
        key(m) = EnzymeRates._canonical(EnzymeRates._law(_testhelper_lift(m)), "tme").key
        @test key(m_ordered) != key(m_random)
```

and in "LDH Pattern-A: graph-distinct mechanisms with equivalent v share a dedup key"
replace

```julia
        key_a = EnzymeRates._rate_eq_dedup_key(rate_equation_string(ldh_m_a))
        @test key_a == EnzymeRates._rate_eq_dedup_key(rate_equation_string(ldh_m_b))
```

with

```julia
        law_a = EnzymeRates._law(_testhelper_lift(ldh_m_a))
        law_b = EnzymeRates._law(_testhelper_lift(ldh_m_b))
        canon_a = EnzymeRates._canonical(law_a, "a")
        canon_b = EnzymeRates._canonical(law_b, "b")
        key_a = canon_a.key
        @test key_a == canon_b.key
        @test _testhelper_same_law_under_map(law_a, law_b, canon_a.order, canon_b.order;
                                             rng = Random.MersenneTwister(13))
```

and

```julia
        @test key_a != EnzymeRates._rate_eq_dedup_key(rate_equation_string(ldh_m_c))
```

with

```julia
        @test key_a != EnzymeRates._canonical(EnzymeRates._law(_testhelper_lift(ldh_m_c)),
                                              "c").key
```

In `test/test_compile_budget.jl`, rename the testset "a second process derives the same
rate_equation bodies" to "a second process derives the same rate_equation bodies and keys",
and in it replace

```julia
            println("BODY ", i, " ", bytes2hex(ER.sha256(string(ER._rate_body(m)))))
```

with

```julia
            println("BODY ", i, " ", bytes2hex(ER.sha256(string(ER._rate_body(m)))))
            println("KEY ", i, " ",
                    ER._canonical(ER._law(m), ER._mechanism_id(String(text))).key)
```

and

```julia
    lines = filter(startswith("BODY "), split(_testhelper_subprocess_output(script), '\n'))
    rm(file; force = true)
    @test length(lines) == 48
    for (i, m) in enumerate(mechs)
        mine = bytes2hex(EnzymeRates.sha256(string(EnzymeRates._rate_body(m))))
        @test "BODY $i $mine" in lines
    end
```

with

```julia
    lines = split(_testhelper_subprocess_output(script), '\n')
    rm(file; force = true)
    @test count(startswith("BODY "), lines) == 48 && count(startswith("KEY "), lines) == 48
    for (i, m) in enumerate(mechs)
        mine = bytes2hex(EnzymeRates.sha256(string(EnzymeRates._rate_body(m))))
        @test "BODY $i $mine" in lines
        key = EnzymeRates._canonical(EnzymeRates._law(m),
                                     EnzymeRates._mechanism_id(sprint(show, m))).key
        @test "KEY $i $key" in lines
    end
```

and its comment

```julia
    # Q7 (spec section 2, criterion 3): a fresh process reads each mechanism back from its
    # text and derives the same generated body, compared by SHA-256, for the 42 specs and
    # the six gate laws.
```

with

```julia
    # Q7 (spec section 2, criterion 3): a fresh process reads each mechanism back from its
    # text and derives the same generated body, compared by SHA-256, and the same key, for
    # the 42 specs and the six gate laws.
```

- [ ] **Step 3: Run the key tests to verify they fail**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T7_key_red.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); @testset "red" begin include("test/test_key.jl") end'`
Expected: every testset errors with ``UndefVarError: `_canonical` not defined in
`EnzymeRates` `` (and ``UndefVarError: `_KEY_LEAF_CAP` ``); the outer `red` testset
reports every one, and the log ends `exit 1`.

Then watch the new population tests and the subprocess's key lines fail too (background,
one after the other, polled):
`S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && for f in test_oracles:oracles_red test_compile_budget:budget_red; do $S/phase3/jl.sh $S/phase45/logs/P5_T7_${f#*:}.log 3000 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); @testset "red" begin include("test/'${f%%:*}'.jl") end'; done`
Expected: in `P5_T7_oracles_red.log`, Q2 ("rewrites, regroupings and constant renames keep
every sampled key") and Q3 ("laws that share a key are one law under the canonical map")
error with the same `UndefVarError`, and every other testset passes, the RE-speed testset
on `_testhelper_i_binding_sample` included. In `P5_T7_budget_red.log`, "a second process
derives the same rate_equation bodies and keys" fails its line counts (the subprocess
stops at `_canonical`, and the helper warns that it failed), then errors at the parent's
`_canonical`; the compile-budget testsets pass. The TID and TME testsets of Step 2 are
re-pins, checked against their values in Step 7.

- [ ] **Step 4: Write the key**

Create `src/key.jl`:

```julia
# ABOUTME: The dedup key: a rate law brought to lowest terms, its fitted constants labelled
# ABOUTME: canonically by colour refinement and individualization, and that text's SHA-256.

"""
A rate law's canonical form (spec section 10): `order`, its fitted constants in canonical
order, the constants no term holds last by name (`L`, a fixed variable, is not listed);
`text`, the canonical text; `key`, the text's SHA-256 in 64 hex digits; `leaves`, the
labellings the search compared; `capped`, whether the search stopped at `_KEY_LEAF_CAP`,
in which case `text` is the mechanism's own (`"mechanism|"` and its id) and `order` is the
variable order. Below the cap, two laws get one key exactly when one is the other with its
fitted constants renamed, each brought to lowest terms unit by unit; concentrations count
by metabolite name.
"""
struct _Canonical
    order::Vector{Symbol}
    text::String
    key::String
    leaves::Int
    capped::Bool
end

"""The most labellings the key's search compares before the law takes its mechanism's
key. Enumerated laws need one; the cap bounds the search on laws with many symmetric
constants."""
const _KEY_LEAF_CAP = 100

"""
    _key_parts(L) -> (components, units, header)

The polynomials of law `L` in role order, the units of components that rescale together
without changing the rate, and a header naming the law's shape: `H`, then `(Q, D)` for a
plain law; for an MWC law `H`, each conformation's `(Q, D, free-E weight)` and each
regulatory site's binding polynomials in both conformations, the header naming the
catalytic and the site multiplicities.
"""
_key_parts(L::_Law) = (_Poly[L.H, L.Q, L.D], [[1], [2, 3]], "plain")

function _key_parts(A::_AlloLaw)
    comps = _Poly[A.H, A.active.Q, A.active.D, A.fA, A.inactive.Q, A.inactive.D, A.fI]
    units = [[1], [2, 3, 4], [5, 6, 7]]
    for (WA, WI, _) in A.sites
        push!(comps, WA, WI)
        push!(units, [length(comps) - 1, length(comps)])
    end
    comps, units, string("mwc|n=", A.n, "|sites=", join((s[3] for s in A.sites), ","))
end

"""The rank of each element of `xs` among its distinct values, 1 for the smallest."""
_ranks(xs::AbstractVector) = (u = sort!(unique(xs)); [searchsortedfirst(u, x) for x in xs])

"""Colour refinement of the constants and the terms from the colours `pcol` and `tcol`,
over the incidences `inc[t] = [(constant, exponent)]` and `occ[p] = [(term, exponent)]`:
the coarsest stable ordered partition below the given one."""
function _refine(pcol::Vector{Int}, tcol::Vector{Int},
                 inc::Vector{Vector{Tuple{Int, _Exp}}},
                 occ::Vector{Vector{Tuple{Int, _Exp}}})
    pcol, tcol = _ranks(pcol), _ranks(tcol)
    while true
        t2 = _ranks([(tcol[t], sort!([(pcol[p], e) for (p, e) in inc[t]]))
                     for t in eachindex(inc)])
        p2 = _ranks([(pcol[p], sort!([(t2[t], e) for (t, e) in occ[p]]))
                     for p in eachindex(occ)])
        maximum(p2; init = 0) == maximum(pcol; init = 0) &&
            maximum(t2; init = 0) == maximum(tcol; init = 0) && return p2, t2
        pcol, tcol = p2, t2
    end
end

"""The terms as sorted lines of text, constant `p` written `x<pcol[p]>`."""
function _encode(pcol::Vector{Int}, inc::Vector{Vector{Tuple{Int, _Exp}}},
                 labels::Vector{String})
    factors(t) = join(("x$(pcol[p])^$e" for (p, e) in sort(inc[t]; by = x -> pcol[x[1]])),
                      "*")
    join(sort!([string(labels[t], "|", factors(t)) for t in eachindex(inc)]), "\n")
end

"""The key search's state: the smallest text found, the constants' order in it, and the
labellings compared."""
mutable struct _KeySearch
    text::String
    order::Vector{Int}
    leaves::Int
end

"""Refine the colours `pcol`, `tcol`; at a labelling, keep it in `st` when its text is the
smallest so far; else individualize each constant of the first tied class in turn and
search on, until `_KEY_LEAF_CAP` labellings."""
function _search!(st::_KeySearch, pcol::Vector{Int}, tcol::Vector{Int}, inc, occ,
                  labels::Vector{String})
    st.leaves >= _KEY_LEAF_CAP && return st
    pcol, tcol = _refine(pcol, tcol, inc, occ)
    if allunique(pcol)
        st.leaves += 1
        text = _encode(pcol, inc, labels)
        if st.leaves == 1 || text < st.text
            st.text, st.order = text, sortperm(pcol)
        end
        return st
    end
    c = minimum(x for x in pcol if count(==(x), pcol) > 1)
    for p in findall(==(c), pcol)
        q = 2 .* pcol
        q[p] -= 1
        _search!(st, q, tcol, inc, occ, labels)
    end
    st
end

"""
    _canonical(L, mechanism_id) -> _Canonical

The canonical form of the rate law `L` (a `_Law` or an `_AlloLaw`). Each unit of
`_key_parts` is multiplied by the monomial that brings every variable to lowest terms
(`_lowest_shift`) and divided by its integer content. A term is labelled by its
component, coefficient and fixed factors (concentrations as `c_<name>`, `Keq`, `L`),
sorted. The fitted constants, every other variable a term holds, are coloured by
refinement on the term–constant incidence; while colours tie, each constant of the first
tied class is individualized in turn, keeping the smallest text (`_search!`). At
`_KEY_LEAF_CAP` labellings the law takes the key of `"mechanism|" * mechanism_id`.
"""
function _canonical(L::Union{_Law, _AlloLaw}, mechanism_id::AbstractString)
    comps, units, header = _key_parts(L)
    vars = L.vars
    fixed(v) = vars.conc[v] || vars.names[v] === :Keq || vars.names[v] === :L
    terms = Tuple{Int, _Mono, Int}[]
    for u in units
        s = _lowest_shift(comps[u])
        g = foldl(gcd, (c for i in u for c in values(comps[i].terms)); init = 0)
        for i in u, (mo, c) in comps[i].terms
            push!(terms, (i, _mmul(mo, s), div(c, g)))
        end
    end
    params = [v for v in eachindex(vars.names) if !fixed(v)]
    held = Set(v for (_, mo, _) in terms for (v, _) in mo)
    present = filter(in(held), params)
    at = Dict(v => j for (j, v) in enumerate(present))
    label(v) = vars.conc[v] ? "c_" * String(vars.names[v]) : String(vars.names[v])
    labels = [string(i, "|", c, "|",
                     join(sort!([string(label(v), "^", e) for (v, e) in mo if fixed(v)]),
                          "*"))
              for (i, mo, c) in terms]
    inc = [Tuple{Int, _Exp}[(at[v], e) for (v, e) in mo if !fixed(v)]
           for (_, mo, _) in terms]
    occ = [Tuple{Int, _Exp}[] for _ in present]
    for (t, ps) in enumerate(inc), (p, e) in ps
        push!(occ[p], (t, e))
    end
    st = _search!(_KeySearch("", Int[], 0), ones(Int, length(present)), _ranks(labels),
                  inc, occ, labels)
    if st.leaves >= _KEY_LEAF_CAP
        text = "mechanism|" * mechanism_id
        return _Canonical(vars.names[params], text, bytes2hex(sha256(text)), st.leaves,
                          true)
    end
    absent = sort!(vars.names[setdiff(params, present)])
    text = string(header, "|P=", length(params), "|absent=", length(absent), "\n", st.text)
    _Canonical([vars.names[present[st.order]]; absent], text, bytes2hex(sha256(text)),
               st.leaves, false)
end
```

In `src/EnzymeRates.jl`, add `using SHA: sha256` after `import Unicode`, and
`include("key.jl")` directly before `include("codegen.jl")`. In
`src/identify_rate_equation.jl`, delete the line `using SHA: sha256`.

- [ ] **Step 5: Run the key tests to verify they pass**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T7_key.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_key.jl")'`
Expected: every testset passes, pristine; the 6-cycle and triangle toys under the cap
(about 12 and 72 labellings).

- [ ] **Step 6: Read and pin five spec keys**

Run on each version:
`S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && for v in 112: 113:+1.13; do $S/phase3/jl.sh $S/phase45/logs/P5_T7_pins${v%%:*}.log 1200 ${v#*:} --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); for n in ("Uni-Uni", "Segel Ordered Bi Bi", "Random-order Bi-Bi", "Segel Ping Pong Bi Bi", "MWC Dimer [AllostericEnzymeMechanism]"); m = _testhelper_lift(only(s for s in MECHANISM_TEST_SPECS if s.name == n).mechanism); println("PIN ", n, " => ", EnzymeRates._canonical(EnzymeRates._law(m), "pin").key); end'; done; diff <(grep "^PIN" $S/phase45/logs/P5_T7_pins112.log) <(grep "^PIN" $S/phase45/logs/P5_T7_pins113.log) && echo "keys agree"`
Expected: five `PIN … => <64 hex>` lines in each log and `keys agree`. A difference between
the versions is a bug (KC: keys agree across 1.10, 1.12 and 1.13): find it before pinning.

Append to `test/test_key.jl`, with the five keys the logs printed in place of the
`<…>` slots:

```julia
@testset "five spec keys are pinned" begin
    # A key is the same in every process and Julia version (spec section 2, criterion 3).
    # These digests were read once and checked on Julia 1.12.5 and 1.13.1.
    pinned = ("Uni-Uni" => "<Uni-Uni key, P5_T7_pins112.log>",
              "Segel Ordered Bi Bi" => "<Segel Ordered Bi Bi key>",
              "Random-order Bi-Bi" => "<Random-order Bi-Bi key>",
              "Segel Ping Pong Bi Bi" => "<Segel Ping Pong Bi Bi key>",
              "MWC Dimer [AllostericEnzymeMechanism]" => "<MWC Dimer key>")
    for (name, key) in pinned
        spec = only(s for s in MECHANISM_TEST_SPECS if s.name == name)
        m = _testhelper_lift(spec.mechanism)
        @test EnzymeRates._canonical(EnzymeRates._law(m), "pin").key == key
    end
end
```

Rerun Step 5's command: every testset passes.

- [ ] **Step 7: Run the files that read the key**

Run (background), one after the other: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && for f in test_oracles:oracles test_identify_rate_equation:tid test_mechanism_enumeration:tme test_compile_budget:budget; do $S/phase3/jl.sh $S/phase45/logs/P5_T7_${f#*:}.log 3000 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/'${f%%:*}'.jl")'; done`, polling each log with the Runner's bounded loop (Global Constraints).
Expected: all four pass, pristine. Q2 counts at least 70 regroupings; Q3's `merged` holds at
least one group (record how many groups and members for the commit); Q4 gives 239 classes,
the value Task 0's baseline holds (`tail -n +2 $S/phase45/p5_base/eqhash.tsv | cut -f2 |
sort -u | wc -l` prints 239; the first line is the header) and K1's; the dedup-key testsets keep their three verdicts, the LDH pair one law
under the canonical map; the subprocess prints 48 bodies and 48 keys equal to the parent's,
and the compile-budget readings stay within their budgets. A Q3 with no merged group would
test nothing: stop and report the population, since KM found same-law siblings in 13–27% of
fits.

- [ ] **Step 8: Measure the key's cost**

Create `$S/phase45/probes/P5_T7_keycost.jl`:

```julia
# ABOUTME: Times the key (_canonical) on the 239 bi-bi init laws and 200 ter-ter init laws, 40 per
# ABOUTME: V×τ bin up to the cap, beside each law's derivation, and counts leaf-cap hits.
using EnzymeRates, Random, Statistics, Printf
const ER = EnzymeRates
const bibi = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end
const terter = @enzyme_reaction begin
    substrates: A[C], B[N], D[X]
    products: P[C], Q[N], R[X]
end
function timed(ms)
    map(ms) do m
        tl = @elapsed L = ER._law(m)
        id = ER._mechanism_id(sprint(show, m))
        ER._canonical(L, id)
        tk = @elapsed c = ER._canonical(L, id)
        (law_ms = 1e3 * tl, key_ms = 1e3 * tk, leaves = c.leaves, capped = c.capped,
         terms = sum(ER._nterms, ER._key_parts(L)[1]))
    end
end
function report(label, rows)
    ks = [r.key_ms for r in rows]
    @printf("%s: %d laws; key ms median %.2f, p90 %.2f, max %.1f; law ms median %.2f; leaves max %d; leaf-cap hits %d; key-part terms max %d\n",
            label, length(rows), median(ks), quantile(ks, 0.9), maximum(ks),
            median(r.law_ms for r in rows), maximum(r.leaves for r in rows),
            count(r -> r.capped, rows), maximum(r.terms for r in rows))
    flush(stdout)
end
report("bi-bi init", timed(ER.init_mechanisms(bibi)))
rng = Xoshiro(20261019)
tt = ER.init_mechanisms(terter)
cx = ER._eq_complexity.(tt)
for (lo, hi) in ((0, 50), (50, 337), (337, 1000), (1000, 2000), (2000, 5000))
    bin = findall(x -> lo < x <= hi, cx)
    pick = tt[bin[sort(randperm(rng, length(bin))[1:min(40, length(bin))])]]
    report("ter-ter init, V×τ in ($lo, $hi]", timed(pick))
end
```

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T7_keycost.log 1800 --project=. --heap-size-hint=2500M $S/phase45/probes/P5_T7_keycost.jl`, polled.
Expected: six report lines, every one with `leaf-cap hits 0` and `leaves max 1`; the largest
bin's key maximum near the slice's 351 ms; `exit 0`, maxrss under 3 GB. A leaf-cap hit on a
real law is not a test failure but a finding: name the mechanism in the commit and the PR
body.

- [ ] **Step 9: Run the full suite**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T7_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`, polled.
Expected: `exit 0`, every test passing, pristine, Aqua and JET green; about 4 s over Task
6's time; report `maxrss_kb`.

- [ ] **Step 10: Commit**

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt
git add src/key.jl src/EnzymeRates.jl src/identify_rate_equation.jl test/test_key.jl \
    test/runtests.jl test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/test_oracles.jl test/test_identify_rate_equation.jl \
    test/test_mechanism_enumeration.jl test/test_compile_budget.jl
git commit -m "$(cat <<'EOF'
Key each rate law by its canonical form

Each unit of the law goes to lowest terms, the fitted constants are labelled
by colour refinement and individualization, concentrations by name, and the
key is the text's full SHA-256; at 100 labellings a law takes its
mechanism's own key. Rewrites, same-model regroupings and constant renames
keep every sampled key; laws that share a key are one law under the
canonical map (<groups> groups, <members> members); the bi-bi init laws
give 239 classes, as before. Five spec keys are pinned, equal on Julia
1.12.5 and 1.13.1, and a second process gives the same keys. Key cost:
bi-bi median <ms> ms; ter-ter by bin <medians>, max <ms> ms; 0 leaf-cap
hits. Full-suite maxrss: <kb> KB.

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
)"
```

---

### Task 8: identify's key cache (P5-20, P5-27, P5-28)

PASS1 derives each mechanism's law from the concrete mechanism and labels it canonically.
The memo holds one fit per key at canonical positions, and every row reads it through its
own canonical order, so a row carries its own parameter names. The CSV gains `key` and
`fit_source_id`; `eq_hash` and `_rate_eq_dedup_key` go. LOOCV fits each fit key's source
once per fold. The CSV's `rate_equation` text is still the old printer's until Task 12,
and `n_params` still comes from the compiled mechanism.

**Files:**
- Modify: `src/identify_rate_equation.jl` (ABOUTME; delete `_rate_eq_dedup_key`; `_Fit`,
  `_Candidate`, `_fit_order`, `_row_params`, `BatchEntry`; `_rows_to_dataframe`;
  `_failure_row`; `_prefit_summary`'s docstring; `_process_batch`; `_ingest!`'s
  docstring; `_offer_cv!`; `_beam_search` (memo, return); `_cv_model_selection`;
  `identify_rate_equation` (call sites; `n_cv_candidates` docstring))
- Modify: `src/rate_eq_derivation.jl` (`_equation_text`'s docstring, which names the
  `eq_hash` dedup)
- Modify: `ROADMAP.md` (retire "Splitting a tied kinetic group changes `eq_hash` without
  changing the model"; reword "Profile the first v0.9.0 HPC runs" and "Seed each identify
  fit from its mechanism's content")
- Test: `test/test_identify_rate_equation.jl`, `test/test_compile_budget.jl` (the PASS1
  workload)
- Create (scratch): `$S/phase45/probes/P5_T8_values.jl`
- Logs: `$S/phase45/logs/P5_T8_{values,tid_red,tid,budget,full}.log`

**Interfaces:**
- Consumes: Task 1's `_solve`, `_Solve.fitted`; Tasks 3–4's `_law(m::Mechanism,
  s::_Solve)`, `_law(am::AllostericMechanism, s::_Solve)`, `_law(m)`; Task 5's
  `_assert_support(m, L)`; Task 7's `_canonical(L, mechanism_id)::_Canonical` (fields
  `order`, `text`, `key`, `capped`, `leaves`) and `sha256` in the `EnzymeRates` namespace;
  `_mechanism_id`, `_eq_complexity`, `fit_rate_equation`, `FittingProblem`,
  `compile_mechanism`, `_cv_fold_loss`, `_select_best_row`. Tests: Task 7's
  `_testhelper_same_law_under_map(L1, L2, order1, order2; rng, swap = false)`.
- Produces:
  ```julia
  struct _Fit
      values::Vector{Float64}   # fitted constants at canonical positions, then L (MWC)
      loss::Float64
      retcode::Symbol
      peak_at_range_edge::Bool
      source::Union{Mechanism, AllostericMechanism}   # the mechanism that was fitted
      source_id::String                               # its _mechanism_id
  end
  struct _Candidate                                   # the PASS1 record
      mech::Union{Mechanism, AllostericMechanism}
      fitted::Vector{Symbol}                          # fitted_params order
      mechanism::String
      mechanism_id::String
      eq_text::String
      canon::_Canonical
  end
  _fit_order(c::_Candidate)::Vector{Symbol}           # c.canon.order, then :L for MWC
  _row_params(c::_Candidate, order::Vector{Symbol}, values::Vector{Float64})::NamedTuple
  struct BatchEntry
      mech::Union{Mechanism, AllostericMechanism}
      n_params::Int
      loss::Float64
      retcode::Symbol
      key::String
      fit_key::String
      row::NamedTuple
  end
  _process_batch(mechs, prob::IdentifyRateEquationProblem; optimizer, max_param_count,
                 eq_complexity_filter::Int = typemax(Int),
                 memo::Dict{String, _Fit} = Dict{String, _Fit}(),
                 seen::Set{UInt64} = Set{UInt64}(), parent_of::AbstractDict = Dict(),
                 log = msg -> nothing, kwargs...)
      # -> (entries, failures, n_param_skip, n_cx_skip, n_seen_skip)
  _offer_cv!(pool::Vector{BatchEntry}, e::BatchEntry, n::Int)::Vector{BatchEntry}
  _beam_search(prob; …)::Tuple{Vector{BatchEntry}, Dict{String, _Fit}}
  _cv_model_selection(pool::Vector{BatchEntry}, memo::Dict{String, _Fit},
                      prob::IdentifyRateEquationProblem; optimizer,
                      se_threshold::Float64, save_dir, show_progress,
                      kwargs...)::IdentifyRateEquationResults
  ```
  - CSV columns, in order: `n_params`, `parent_n_params`, `loss`, `mechanism`,
    `mechanism_id`, `parent_id`, `rate_equation`, `retcode`, `error`, `key`,
    `fit_inherited`, `fit_source_id`, `peak_at_range_edge`, then the parameters. A failure
    row has `missing` in `key`, `fit_inherited` and `fit_source_id`.
  - A capped key adds the line `leaf cap reached for <mechanism_id>; fitted under its own
    key` to the batch's pre-fit summary.
  - Test helpers (TID): `_testhelper_SpreadStubOpt`, `_testhelper_unibi_pair()`; the
    testset-local `_testhelper_pool_entry`.

**Gate:**
- `P5_T8_values.jl` prints `one key true` and `same law true` for each of the four
  same-key pairs before any pin moves, and `swapped map false` unless the commit records
  that the two swapped constants play symmetric roles (Step 1).
- `test_identify_rate_equation.jl` is green, pins 18 and 19 moved after their value checks.
- The uni-bi pair fits once: the stub counts 1 solve, `fit_inherited` is true on exactly
  one row, the two rows' parameter names differ, and their rates agree to 1e-12.
- `fit_inherited` means the row is not its fit's source: "fit-dedup by key" re-pins the
  source's own row in a later batch from `true` to `false`.
- The compile budgets hold, read from `P5_T8_budget.log` before the full run.
- The full suite is green and pristine; about +8 s (the uni-bi pair, the pinwheel, whose
  key search alone takes about 3 s, and the LDH pair compile and fit in `_process_batch`);
  maxrss reported against Task 0's.

- [ ] **Step 1: Check the values the re-pins rest on**

Four tests change their dedup assertion from `eq_hash` to the key. Each pair must have one
key and be one law under the canonical map, and a swapped map must fail, before its pin
moves (spec section 14). On a9e0699 the dedup twins, the LDH renaming pair and the LDH
allosteric split/merge pair have laws equal field by field, names included, and the
uni-bi pair differs only in one constant's name (`$S/phase45/plan5/probes/p3_probe2.log`,
section K; `p3_probe.log`, section D). With Task 7's key the twins and the uni-bi pair
print `one key true; same law true; swapped map false` (`p3c_probe.log`, section C).

Create `$S/phase45/probes/P5_T8_values.jl`:

```julia
# ABOUTME: Task 8's value checks before its re-pins: each same-key pair of identify's tests
# ABOUTME: has one key and is one law under the canonical map, and a swapped map fails.
# Usage, from the worktree root: julia --project=. P5_T8_values.jl
using TestEnv; TestEnv.activate()
using Test, EnzymeRates, LinearAlgebra, Random
include(joinpath(pkgdir(EnzymeRates), "test",
                 "mechanism_definitions_for_test_enzyme_derivation.jl"))
const ER = EnzymeRates
canon(m) = ER._canonical(ER._law(m), ER._mechanism_id(sprint(show, m)))
twin_a = ER.AllostericMechanism(@allosteric_mechanism begin
    substrates: S
    products: P
    catalytic_multiplicity: 2
    steps: begin
        E + S ⇌ E(S)       :: OnlyA
        E(S) <--> E(P)     :: NonequalAI
        E(P) ⇌ E + P       :: OnlyA
    end
end)
twin_b = ER.AllostericMechanism(@allosteric_mechanism begin
    substrates: S
    products: P
    catalytic_multiplicity: 2
    steps: begin
        E + S ⇌ E(S)       :: OnlyA
        E(S) <--> E(P)     :: OnlyA
        E(P) ⇌ E + P       :: OnlyA
    end
end)
ldh_merged = ER.Mechanism(@enzyme_mechanism begin
    substrates: NADH, Pyruvate
    products:   Lactate, NAD
    steps: begin
        (E + Lactate ⇌ E(Lactate), E(NAD) + Lactate ⇌ E(Lactate, NAD),
         E(NADH) + Lactate ⇌ E(Lactate, NADH))
        (E + NAD ⇌ E(NAD), E(Lactate) + NAD ⇌ E(Lactate, NAD))
        (E + NADH ⇌ E(NADH), E(Lactate) + NADH ⇌ E(Lactate, NADH))
        E + Pyruvate ⇌ E(Pyruvate)
        (E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate),
         E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate))
        E(NADH, Pyruvate) <--> E(Lactate, NAD)
        E(Pyruvate) + NADH ⇌ E(NADH, Pyruvate)
        E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate)
    end
end)
ldh_split = ER.Mechanism(@enzyme_mechanism begin
    substrates: NADH, Pyruvate
    products:   Lactate, NAD
    steps: begin
        (E + Lactate ⇌ E(Lactate), E(NAD) + Lactate ⇌ E(Lactate, NAD))
        (E + NAD ⇌ E(NAD), E(Lactate) + NAD ⇌ E(Lactate, NAD))
        (E + NADH ⇌ E(NADH), E(Lactate) + NADH ⇌ E(Lactate, NADH))
        E + Pyruvate ⇌ E(Pyruvate)
        (E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate),
         E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate))
        E(NADH) + Lactate ⇌ E(Lactate, NADH)
        E(NADH, Pyruvate) <--> E(Lactate, NAD)
        E(Pyruvate) + NADH ⇌ E(NADH, Pyruvate)
        E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate)
    end
end)
ldh_allo_split = ER.AllostericMechanism(@allosteric_mechanism begin
    substrates: NADH, Pyruvate
    products:   Lactate, NAD
    catalytic_multiplicity: 4
    steps: begin
        (E + Lactate ⇌ E(Lactate), E(NAD) + Lactate ⇌ E(Lactate, NAD))   :: EqualAI
        (E + NAD ⇌ E(NAD), E(Lactate) + NAD ⇌ E(Lactate, NAD),
         E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate))                           :: EqualAI
        (E + NADH <--> E(NADH), E(Lactate) + NADH <--> E(Lactate, NADH),
         E(Pyruvate) + NADH <--> E(NADH, Pyruvate))                      :: OnlyA
        E + Pyruvate ⇌ E(Pyruvate)                                       :: NonequalAI
        E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate)                             :: EqualAI
        E(NADH) + Lactate ⇌ E(Lactate, NADH)                             :: EqualAI
        E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate)                           :: EqualAI
        E(NADH, Pyruvate) <--> E(Lactate, NAD)                           :: OnlyA
    end
end)
ldh_allo_merged = ER.AllostericMechanism(@allosteric_mechanism begin
    substrates: NADH, Pyruvate
    products:   Lactate, NAD
    catalytic_multiplicity: 4
    steps: begin
        (E + Lactate ⇌ E(Lactate), E(NAD) + Lactate ⇌ E(Lactate, NAD),
         E(NADH) + Lactate ⇌ E(Lactate, NADH))                           :: EqualAI
        (E + NAD ⇌ E(NAD), E(Lactate) + NAD ⇌ E(Lactate, NAD),
         E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate))                           :: EqualAI
        (E + NADH <--> E(NADH), E(Lactate) + NADH <--> E(Lactate, NADH),
         E(Pyruvate) + NADH <--> E(NADH, Pyruvate))                      :: OnlyA
        E + Pyruvate ⇌ E(Pyruvate)                                       :: NonequalAI
        E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate)                             :: EqualAI
        E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate)                           :: EqualAI
        E(NADH, Pyruvate) <--> E(Lactate, NAD)                           :: OnlyA
    end
end)
unibi_equal = ER.AllostericMechanism(@allosteric_mechanism begin
    substrates: A
    products: P, Q
    catalytic_multiplicity: 2
    steps: begin
        E + A <--> E(A)          :: OnlyA
        E(A) <--> E(P, Q)        :: OnlyA
        E(Q) + P ⇌ E(P, Q)       :: EqualAI
        E + Q ⇌ E(Q)             :: OnlyA
    end
end)
unibi_only = ER.AllostericMechanism(@allosteric_mechanism begin
    substrates: A
    products: P, Q
    catalytic_multiplicity: 2
    steps: begin
        E + A <--> E(A)          :: OnlyA
        E(A) <--> E(P, Q)        :: OnlyA
        E(Q) + P ⇌ E(P, Q)       :: OnlyA
        E + Q ⇌ E(Q)             :: OnlyA
    end
end)
for (label, a, b) in (("dedup twins", twin_a, twin_b),
                      ("LDH renaming pair", ldh_merged, ldh_split),
                      ("LDH allosteric split/merge pair", ldh_allo_split, ldh_allo_merged),
                      ("uni-bi pair", unibi_equal, unibi_only))
    ca, cb = canon(a), canon(b)
    La, Lb = ER._law(a), ER._law(b)
    same = _testhelper_same_law_under_map(La, Lb, ca.order, cb.order;
                                          rng = MersenneTwister(1))
    swapped = _testhelper_same_law_under_map(La, Lb, ca.order, cb.order;
                                             rng = MersenneTwister(1), swap = true)
    println(label, ": one key ", ca.key == cb.key, "; same law ", same,
            "; swapped map ", swapped)
end
```

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T8_values.log 1800 --project=. $S/phase45/probes/P5_T8_values.jl`
Expected: four lines, each `one key true; same law true; swapped map false`, and
`exit 0`. A pair whose laws are equal but whose keys differ is a key bug (Task 7): stop
and debug it (systematic-debugging). A pair whose laws differ is a wrong premise in its
test: stop and report it to the controller. A `swapped map true` means the two swapped
constants play symmetric roles in that law, so the control says nothing there: record it
in the commit; it is no failure.

- [ ] **Step 2: Write the failing tests**

All edits are in `test/test_identify_rate_equation.jl` unless a step names another file.
Each "replace" quotes the text as f9361cf holds it.

(a) Inside the testset "identify_rate_equation", directly after the two lines that build
`best_row` (`best_row = results.cv_results[` and
`EnzymeRates._select_best_row(results.cv_results), :]`), add:

```julia
    # The rate-law key of a mechanism's text: the key PASS1 gives it.
    _testhelper_key_of(text) =
        EnzymeRates._canonical(EnzymeRates._law(EnzymeRates._parse_mechanism(text)),
                               EnzymeRates._mechanism_id(text)).key
```

(b) Testset "_rows_to_dataframe": replace the first rows literal's line

```julia
            eq_hash = "0123456789abcdef",
```

with

```julia
            key = repeat("0123456789abcdef", 4),
            fit_source_id = "0123456789abcdef",
```

and replace

```julia
        @test names(df) == ["n_params", "parent_n_params", "loss", "mechanism",
                            "mechanism_id", "parent_id", "rate_equation", "retcode",
                            "error", "eq_hash", "fit_inherited", "peak_at_range_edge",
                            "a", "b"]
```

with

```julia
        @test names(df) == ["n_params", "parent_n_params", "loss", "mechanism",
                            "mechanism_id", "parent_id", "rate_equation", "retcode",
                            "error", "key", "fit_inherited", "fit_source_id",
                            "peak_at_range_edge", "a", "b"]
        @test df.key == [repeat("0123456789abcdef", 4)]
        @test df.fit_source_id == ["0123456789abcdef"]
```

Replace `@test "eq_hash" in names(df)` with `@test "key" in names(df)`. In
`rows_with_failure`, replace the fitted row's `eq_hash = "0123456789abcdef",` with
`key = "0123456789abcdef", fit_source_id = "0000000000000001",` and the failure row's
`eq_hash = missing,` with `key = missing, fit_source_id = missing,`. Replace
`@test ismissing(df_with_failure.eq_hash[2])` with

```julia
        @test ismissing(df_with_failure.key[2])
        @test ismissing(df_with_failure.fit_source_id[2])
```

In the last row literal (`… params = (;), eq_hash = "0", fit_inherited = false)`),
replace `eq_hash = "0",` with `key = "0", fit_source_id = "0",`.

(c) Testset "a failure row holds its mechanism's text and id": before its `end`, add

```julia
        @test ismissing(row.key) && ismissing(row.fit_inherited) &&
              ismissing(row.fit_source_id)
```

(d) Replace the whole testset "_rate_eq_dedup_key" (its subject goes; listed for Denis)
with:

```julia
    @testset "the key of a law is the SHA-256 of its canonical text" begin
        # identify fits one law per key (spec section 10). The key reads the law, not its
        # printed text: it is the 64-hex SHA-256 of the canonical text, the same on every
        # call, and another law gets another key.
        m = _testhelper_uni_seed
        id = EnzymeRates._mechanism_id(sprint(show, m))
        c = EnzymeRates._canonical(EnzymeRates._law(m), id)
        again = EnzymeRates._canonical(EnzymeRates._law(m), id)
        @test (again.key, again.text, again.order) == (c.key, c.text, c.order)
        @test c.key == bytes2hex(EnzymeRates.sha256(c.text))
        @test occursin(r"^[0-9a-f]{64}$", c.key)
        ss_release = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) <--> E + P
            end
        end)
        other = EnzymeRates._mechanism_id(sprint(show, ss_release))
        @test EnzymeRates._canonical(EnzymeRates._law(ss_release), other).key != c.key
    end
```

(e) Testset "results structure": replace

```julia
        @test "eq_hash" in names(
            results.cv_results)
        # LOOCV candidate dedup invariant: within each n_params
        # bucket, each eq_hash should appear at most once (the cv
        # pool `_offer_cv!` keeps one slot per eq_hash). This
        # catches a regression where duplicates would enter LOOCV
        # and waste compute / bias the per-bucket "best".
        for gdf in groupby(
                results.cv_results, :n_params)
            @test allunique(gdf.eq_hash)
        end
```

with

```julia
        @test "key" in names(results.cv_results)
        @test "fit_source_id" in names(results.cv_results)
        # LOOCV candidate dedup invariant: within each n_params bucket, each key appears
        # at most once (the cv pool `_offer_cv!` keeps one slot per key). This catches a
        # regression where duplicates would enter LOOCV and waste compute / bias the
        # per-bucket "best". Each row's key is its mechanism's.
        for gdf in groupby(results.cv_results, :n_params)
            @test allunique(gdf.key)
        end
        @test all(r -> r.key == _testhelper_key_of(r.mechanism),
                  eachrow(results.cv_results))
```

(f) Testset "CSV output: every file's columns, texts and ids": replace

```julia
        _testhelper_read_csv(f) = CSV.read(joinpath(save_dir, f), DataFrame;
            types = Dict(:mechanism_id => String, :parent_id => String))
```

with

```julia
        _testhelper_read_csv(f) = CSV.read(joinpath(save_dir, f), DataFrame;
            types = Dict(:mechanism_id => String, :parent_id => String,
                         :fit_source_id => String))
```

replace `@test "eq_hash" in names(init_df)` with `@test "key" in names(init_df)`, and
replace

```julia
            @test "eq_hash" in names(df_file)
            @test all(length.(string.(skipmissing(df_file.eq_hash))) .== 16)
```

with

```julia
            @test "key" in names(df_file) && "fit_source_id" in names(df_file)
            @test all(k -> occursin(r"^[0-9a-f]{64}$", k), skipmissing(df_file.key))
            # Each fitted row's key is its mechanism's, its fit's source is a 16-hex id,
            # and the row inherits its fit exactly when the source is another mechanism.
            for r in eachrow(df_file)
                ismissing(r.key) && continue
                @test r.key == _testhelper_key_of(r.mechanism)
                @test occursin(r"^[0-9a-f]{16}$", r.fit_source_id)
                @test r.fit_inherited == (r.fit_source_id != r.mechanism_id)
            end
```

(g) Testset "loocv_results.csv and best_equation.csv saved": replace

```julia
        best_hash = string(
            EnzymeRates._rate_eq_dedup_key(rate_equation_string(results.best)),
            base=16, pad=16)
        @test string(bestf.eq_hash[1]) == best_hash
```

with

```julia
        # The key reads back as text, and it is the selected mechanism's key.
        @test bestf.key[1] isa AbstractString
        @test bestf.key[1] == _testhelper_key_of(best_row.mechanism)
```

(h) Testset "csv writer": in the rows literal replace `eq_hash = "abc",` with
`key = "abc", fit_source_id = "00000000000000aa",`, and replace
`@test "eq_hash" in names(df)` with `@test "key" in names(df) && df.key == ["abc"]`.

(i) Testset "all base fits fail: failure CSV written, then raises": replace the comment
line `# are all failures (non-missing `error`, missing `eq_hash`).` with
`# are all failures (non-missing `error`, missing `key`).` and replace
`@test all(ismissing.(fail_df.eq_hash))` with `@test all(ismissing.(fail_df.key))`.

(j) Testset "_progress": replace

```julia
    row = (n_params=3, loss=0.5, mechanism="M", rate_equation="v",
           retcode="Success", error=missing, params=(K = 1.0,), eq_hash="abc",
           fit_inherited=false)
    e_succ = EnzymeRates.BatchEntry(mech, 3, 0.5, :Success, hash(:a), row)
    e_mt   = EnzymeRates.BatchEntry(mech, 3, 0.9, :MaxTime, hash(:b), row)
```

with

```julia
    row = (n_params=3, loss=0.5, mechanism="M", rate_equation="v",
           retcode="Success", error=missing, params=(K = 1.0,), key="abc",
           fit_inherited=false)
    e_succ = EnzymeRates.BatchEntry(mech, 3, 0.5, :Success, "a", "a", row)
    e_mt   = EnzymeRates.BatchEntry(mech, 3, 0.9, :MaxTime, "b", "b", row)
```

(k) Testset "_process_batch": replace
`@test all(e -> occursin(r"^[0-9a-f]{16}$", e.row.eq_hash), entries)` with

```julia
    @test all(e -> occursin(r"^[0-9a-f]{64}$", e.key) && e.key == e.row.key == e.fit_key,
              entries)
```

and, in this testset and in "_process_batch failures report the ORIGINAL mechanism",
replace every `memo=Dict{UInt64,NamedTuple}()` and `memo=Dict{UInt64, NamedTuple}()` with
`memo=Dict{String,EnzymeRates._Fit}()` (four places).

(l) Testset "_ingest! and cv pool": replace

```julia
    _testhelper_mk(n, loss, h) = EnzymeRates.BatchEntry(
        _testhelper_uni_seed, n, loss, :Success, hash(h),
        (n_params=n, loss=loss, mechanism="M",
         rate_equation="v", retcode="Success", error=missing,
         params=(K = 1.0,), eq_hash=string(hash(h),base=16,pad=16),
         fit_inherited=false))
```

with

```julia
    _testhelper_mk(n, loss, h) = EnzymeRates.BatchEntry(
        _testhelper_uni_seed, n, loss, :Success, string(h), string(h),
        (n_params=n, loss=loss, mechanism="M",
         rate_equation="v", retcode="Success", error=missing,
         params=(K = 1.0,), key=string(h), fit_inherited=false))
```

Replace

```julia
    @test length(cv_pool[5]) == 2              # bounded, distinct eq_hash
    # the kept :a entry is the lower-loss one (2.0, not 3.0);
    # BatchEntry.eq_hash is a UInt64, so compare against hash(:a), not hex:
    a = only(filter(e -> e.eq_hash == hash(:a), cv_pool[5]))
```

with

```julia
    @test length(cv_pool[5]) == 2              # bounded, distinct keys
    # the kept :a entry is the lower-loss one (2.0, not 3.0)
    a = only(filter(e -> e.key == "a", cv_pool[5]))
```

replace `@test [e.eq_hash for e in frontier_nf[5]] == [hash(:h)]` with
`@test [e.key for e in frontier_nf[5]] == ["h"]`, and
`@test [e.eq_hash for e in cv_pool_nf[5]] == [hash(:h)]` with
`@test [e.key for e in cv_pool_nf[5]] == ["h"]`. Replace

```julia
    # `_offer_cv!` keeps at most one entry per eq_hash: a repeat hash updates its own
    # slot to the lower loss, never consuming a second.
```

with

```julia
    # `_offer_cv!` keeps at most one entry per key: a repeat key updates its own slot to
    # the lower loss, never consuming a second.
```

and in the lines after it replace `allunique([e.eq_hash for e in pool])` with
`allunique([e.key for e in pool])` and `filter(e -> e.eq_hash == hash(:a), pool)` with
`filter(e -> e.key == "a", pool)`.

(m) Testsets "_expand_parents records an expansion error instead of aborting" and
"_expand_parents matches each parent's expand_mechanisms": replace
`EnzymeRates.BatchEntry(m, 0, 0.0, :Success, UInt64(0), (mechanism="M",))` with
`EnzymeRates.BatchEntry(m, 0, 0.0, :Success, "", "", (mechanism="M",))` and
`EnzymeRates.BatchEntry(m, 3, 0.0, :Success, hash(m), (mechanism="M",))` with
`EnzymeRates.BatchEntry(m, 3, 0.0, :Success, "", "", (mechanism="M",))`. The testset
"n_cv_candidates below 1 is refused before any fit" builds no `BatchEntry` and needs no
edit.

(n) After the `_testhelper_dedup_twins` function, replace the whole testset "fit-dedup by
eq_hash in _process_batch" with:

```julia
@testset "fit-dedup by key in _process_batch" begin
    m1, m2 = _testhelper_dedup_twins()
    em1 = EnzymeRates.compile_mechanism(m1)
    em2 = EnzymeRates.compile_mechanism(m2)
    c1, c2 = (EnzymeRates._canonical(EnzymeRates._law(m),
                                     EnzymeRates._mechanism_id(sprint(show, m)))
              for m in (m1, m2))
    # Preconditions: distinct structure and compiled types, one key and one law under
    # the canonical map, and the same fitted names.
    @test m1 != m2
    @test typeof(em1) != typeof(em2)
    @test c1.key == c2.key
    @test _testhelper_same_law_under_map(EnzymeRates._law(m1), EnzymeRates._law(m2),
                                         c1.order, c2.order;
                                         rng = Random.MersenneTwister(1))
    @test EnzymeRates.fitted_params(em1) == EnzymeRates.fitted_params(em2)

    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            S = [1.0, 2.0, 1.0, 2.0], P = [0.1, 0.2, 0.1, 0.2])
    # The range stops below the stub fit's binding constants (5), so the fit still rises at
    # the range's upper edge and carries the edge flag.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        oligomeric_state: 2
    end
    prob = IdentifyRateEquationProblem(rxn, data; Keq=2.0, peak_range=(0.1, 1.0))
    pair = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m1, m2]

    memo = Dict{String, EnzymeRates._Fit}()
    opt = _testhelper_CountingStubOpt(; uval = log(5.0))
    # Each logged line is stored with the solve count at the moment it was logged.
    batch_log = Tuple{String,Int}[]
    entries, failures = EnzymeRates._process_batch(pair, prob;
        optimizer=opt, max_param_count=20, eq_complexity_filter=337, n_restarts=1,
        maxtime=1.0, memo, log = msg -> push!(batch_log, (msg, opt.count)))

    # The shared law is fit exactly ONCE (n_restarts=1 → one solve).
    @test opt.count == 1
    # The pre-fit summary counts one new fit and one inherited row, and is logged
    # before the fit runs.
    @test batch_log == [("1 new fits + 1 inherited + 0 skipped (already fit) + " *
                         "0 skipped (>20 params) + 0 skipped (>337 complexity)", 0)]
    @test length(entries) == 2
    @test isempty(failures)
    # loss + retcode belong to the law's fit → one key → identical.
    @test entries[1].loss == entries[2].loss
    @test entries[1].retcode === entries[2].retcode === :Success
    # Representative fit first (false); the duplicate is inherited (true), and both
    # name the representative as the fit's source.
    @test [e.row.fit_inherited for e in entries] == [false, true]
    @test [e.row.fit_source_id for e in entries] == fill(entries[1].row.mechanism_id, 2)

    # The memo stores the fit at canonical positions; every row reads it through its own
    # canonical order, so both rows carry the fit's values under their own names.
    fit = memo[c1.key]
    @test fit.source == m1
    fkeys = EnzymeRates.fitted_params(em1)
    for (e, c) in zip(entries, (c1, c2))
        @test keys(e.row.params) == fkeys
        @test Float64[e.row.params[x] for x in [c.order; :L]] == fit.values
        # The edge flag belongs to the law's fit, so each row copies it.
        @test e.row.peak_at_range_edge === fit.peak_at_range_edge
    end
    @test fit.peak_at_range_edge === true
    @test entries[1].row.params == entries[2].row.params
    # The peak scale: the copied params are the rescaled fit, not the raw
    # 5.0 the stub optimizer returned.
    @test !all(v -> v ≈ 5.0, entries[1].row.params)

    # Cross-batch memo hit: a later batch with the same key refits NOTHING.
    single = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m1]
    reuse_log = String[]
    reused, _ = EnzymeRates._process_batch(single, prob;
        optimizer=opt, max_param_count=20, eq_complexity_filter=337, n_restarts=1,
        maxtime=1.0, memo, log = msg -> push!(reuse_log, msg))
    @test opt.count == 1                                # no new solve
    # The mechanism that was fitted reads its own fit in a later batch: not inherited.
    @test [e.row.fit_inherited for e in reused] == [false]
    @test reuse_log == ["0 new fits + 1 inherited + 0 skipped (already fit) + " *
                        "0 skipped (>20 params) + 0 skipped (>337 complexity)"]

    # A representative whose fit throws fails ALL the mechanisms with its key
    # (all-or-nothing per law).
    opt_bad = _testhelper_CountingStubOpt(; throwit=true)
    bad_entries, bad_failures = EnzymeRates._process_batch(pair, prob;
        optimizer=opt_bad, max_param_count=20, n_restarts=1, maxtime=1.0,
        memo = Dict{String, EnzymeRates._Fit}())
    @test isempty(bad_entries)
    @test length(bad_failures) == 2
    @test all(f -> f isa EnzymeRates.FitFailure, bad_failures)
    @test [f.mech for f in bad_failures] == [m1, m2]   # the mechanisms as handed in
end

# A stub optimizer that counts `solve` calls and returns the log-space point
# u_i = log(i + 1), so each fitted constant gets a value of its own and a row read through
# the wrong names gives another rate.
mutable struct _testhelper_SpreadStubOpt
    count::Int
end
_testhelper_SpreadStubOpt() = _testhelper_SpreadStubOpt(0)
Optimization.allowsbounds(::_testhelper_SpreadStubOpt) = true
function Optimization.SciMLBase.__solve(
        prob::Optimization.OptimizationProblem, opt::_testhelper_SpreadStubOpt; kwargs...)
    opt.count += 1
    u = [log(i + 1.0) for i in eachindex(prob.u0)]
    cache = DefaultOptimizationCache(prob.f, prob.p)
    build_solution(cache, opt, u, prob.f(u, prob.p); retcode = ReturnCode.Success)
end

# Two allosteric uni-bi mechanisms with one law and different fitted names (KM §4): A and
# Q bind only the active conformation (:OnlyA), so the inactive one keeps free E alone, and
# P's binding at E(Q), which that conformation never reaches, is :EqualAI in the first and
# :OnlyA in the second. The first names its constant without a conformation
# (`K_EPQ_to_EQ_P`), the second as the active conformation's (`K_A_EPQ_to_EQ_P`), so the
# two sort it to different fitted positions. Returns the two and a problem on their
# reaction.
function _testhelper_unibi_pair()
    equal = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A
        products: P, Q
        catalytic_multiplicity: 2
        steps: begin
            E + A <--> E(A)          :: OnlyA
            E(A) <--> E(P, Q)        :: OnlyA
            E(Q) + P ⇌ E(P, Q)       :: EqualAI
            E + Q ⇌ E(Q)             :: OnlyA
        end
    end)
    only_active = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A
        products: P, Q
        catalytic_multiplicity: 2
        steps: begin
            E + A <--> E(A)          :: OnlyA
            E(A) <--> E(P, Q)        :: OnlyA
            E(Q) + P ⇌ E(P, Q)       :: OnlyA
            E + Q ⇌ E(Q)             :: OnlyA
        end
    end)
    rxn = @enzyme_reaction begin
        substrates: A[CN]
        products: P[C], Q[N]
        oligomeric_state: 2
    end
    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            A = [1.0, 2.0, 1.0, 2.0], P = [0.1, 0.2, 0.1, 0.2], Q = [0.3, 0.3, 0.4, 0.4])
    equal, only_active, IdentifyRateEquationProblem(rxn, data; Keq = 2.0)
end

@testset "two same-key mechanisms read one fit under their own names" begin
    # Spec section 10: identify fits one representative per key and stores the fit at
    # canonical positions; the other mechanism reads it through its own labelling, so its
    # row carries its own names and both rows compute one rate.
    a, b, prob = _testhelper_unibi_pair()
    ca, cb = (EnzymeRates._canonical(EnzymeRates._law(m),
                                     EnzymeRates._mechanism_id(sprint(show, m)))
              for m in (a, b))
    fa, fb = EnzymeRates.fitted_params(a), EnzymeRates.fitted_params(b)
    @test ca.key == cb.key
    @test Set(fa) != Set(fb)
    # The canonical map sends a constant to another fitted position, so a fit copied by
    # position would hand the second row wrong values.
    @test [findfirst(==(x), fa) for x in ca.order] !=
          [findfirst(==(x), fb) for x in cb.order]
    opt = _testhelper_SpreadStubOpt()
    entries, failures = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[a, b], prob;
        optimizer = opt, max_param_count = 20, n_restarts = 1, maxtime = 1.0)
    @test opt.count == 1
    @test isempty(failures) && length(entries) == 2
    @test count(e -> e.row.fit_inherited, entries) == 1
    @test [e.key for e in entries] == [ca.key, ca.key]
    @test keys(entries[1].row.params) == fa && keys(entries[2].row.params) == fb
    rng = Random.MersenneTwister(10)
    for _ in 1:3
        concs = (A = exp(randn(rng)), P = exp(randn(rng)), Q = exp(randn(rng)))
        rates = [rate_equation(EnzymeRates.compile_mechanism(e.mech), concs,
                               merge(e.row.params, (Keq = 2.0, E_total = 1.0)))
                 for e in entries]
        @test isapprox(rates[1], rates[2]; rtol = 1e-12)
    end
end

@testset "a fit is stored by canonical position and read through each row's order" begin
    # The memo's values sit at the canonical positions of the representative's law (its
    # `_Canonical.order`, then L); each row's own canonical order reads them back.
    a, b, prob = _testhelper_unibi_pair()
    memo = Dict{String, EnzymeRates._Fit}()
    entries, _ = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[a, b], prob;
        optimizer = _testhelper_SpreadStubOpt(), max_param_count = 20, memo,
        n_restarts = 1, maxtime = 1.0)
    fit = only(values(memo))
    @test fit.source == a && fit.source_id == entries[1].row.mechanism_id
    @test fit.loss == entries[1].loss == entries[2].loss
    for e in entries
        c = EnzymeRates._canonical(EnzymeRates._law(e.mech), e.row.mechanism_id)
        @test Float64[e.row.params[x] for x in [c.order; :L]] == fit.values
    end
end

@testset "_offer_cv! keeps one slot per key" begin
    _testhelper_entry(key, loss) = EnzymeRates.BatchEntry(_testhelper_uni_seed, 3, loss,
                                                          :Success, key, key, (key = key,))
    pool = EnzymeRates.BatchEntry[]
    EnzymeRates._offer_cv!(pool, _testhelper_entry("a", 1.0), 2)
    EnzymeRates._offer_cv!(pool, _testhelper_entry("b", 2.0), 2)
    # A full pool: a repeat key with a lower loss takes its own slot, not the worst one.
    EnzymeRates._offer_cv!(pool, _testhelper_entry("a", 0.5), 2)
    @test [(e.key, e.loss) for e in pool] == [("a", 0.5), ("b", 2.0)]
    # A repeat key with a higher loss changes nothing; a new key evicts the worst.
    EnzymeRates._offer_cv!(pool, _testhelper_entry("a", 3.0), 2)
    EnzymeRates._offer_cv!(pool, _testhelper_entry("c", 1.5), 2)
    @test [(e.key, e.loss) for e in pool] == [("a", 0.5), ("c", 1.5)]
end

@testset "LOOCV fits each fit key's source once per fold" begin
    # Spec section 10, P5-28: rows that read one fit share their fold losses, so LOOCV fits
    # that fit's source once per fold. Handed to model selection beside a second law, the
    # dedup twins (one key) cost two fits per fold, not three, and share their fold columns.
    m1, m2 = _testhelper_dedup_twins()
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        oligomeric_state: 2
    end
    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            S = [1.0, 2.0, 1.0, 2.0], P = [0.1, 0.2, 0.1, 0.2])
    prob = IdentifyRateEquationProblem(rxn, data; Keq = 2.0)
    memo = Dict{String, EnzymeRates._Fit}()
    entries, _ = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[
            m1, m2, _testhelper_uni_seed], prob;
        optimizer = _testhelper_CountingStubOpt(; uval = log(5.0)), max_param_count = 20,
        memo, n_restarts = 1, maxtime = 1.0)
    @test length(unique(e.fit_key for e in entries)) == 2
    opt = _testhelper_CountingStubOpt(; uval = log(5.0))
    res = EnzymeRates._cv_model_selection(entries, memo, prob; optimizer = opt,
        se_threshold = 1.0, save_dir = mktempdir(), show_progress = false,
        n_restarts = 1, maxtime = 1.0)
    groups = unique(prob.data.group)
    @test opt.count == 2 * length(groups)
    folds(r) = [r[Symbol("cv_fold_$g")] for g in groups]
    twins = [r for r in eachrow(res.cv_results) if r.key == entries[1].key]
    @test length(twins) == 2 && folds(twins[1]) == folds(twins[2])
end

@testset "best_equation.csv reads back with its key and ids" begin
    # Review Focus 2. In a one-row CSV, CSV.jl guesses a column's type from one value: a
    # 16-hex id of digits and one `e` would read as a number, so the documented read names
    # the three id columns `String`; the 64-hex key reads back as text without help.
    prob = _testhelper_uni_prob(NamedTuple)
    memo = Dict{String, EnzymeRates._Fit}()
    entries, _ = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[_testhelper_uni_seed],
        prob; optimizer = _testhelper_CountingStubOpt(; uval = 0.0), max_param_count = 20,
        memo, n_restarts = 1, maxtime = 1.0)
    e = only(entries)
    ids = (mechanism_id = "1234567e89012345", parent_id = "0123456789012345",
           fit_source_id = "1234567e89012345")
    pool = [EnzymeRates.BatchEntry(e.mech, e.n_params, e.loss, e.retcode, e.key, e.fit_key,
                                   merge(e.row, ids))]
    mktempdir() do dir
        EnzymeRates._cv_model_selection(pool, memo, prob;
            optimizer = _testhelper_CountingStubOpt(; uval = 0.0), se_threshold = 1.0,
            save_dir = dir, show_progress = false, n_restarts = 1, maxtime = 1.0)
        path = joinpath(dir, "best_equation.csv")
        plain = CSV.read(path, DataFrame)
        @test plain.key[1] isa AbstractString
        @test plain.key[1] == e.key && occursin(r"^[0-9a-f]{64}$", plain.key[1])
        typed = CSV.read(path, DataFrame; types = Dict(:mechanism_id => String,
                                                       :parent_id => String,
                                                       :fit_source_id => String))
        @test (typed.mechanism_id[1], typed.parent_id[1], typed.fit_source_id[1]) ==
              Tuple(ids)
    end
end

@testset "a law at the leaf cap is fitted under its mechanism's key" begin
    # Five two-step cycles E → X → E, each taking up A and giving off P twice, alike up to
    # the label of X: every permutation of the cycles maps the law onto itself, so a
    # canonical labelling would compare 5! = 120 labellings and the search stops at the
    # leaf cap (P5-19). The mechanism is fitted under a key of its own, the SHA-256 of
    # "mechanism|" and its mechanism_id, and the pre-fit summary names it. (The key search
    # takes about 3 s here, so it runs once, inside the batch.)
    rxn = @enzyme_reaction begin
        substrates: A[C]
        products: P[C]
    end
    pinwheel = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
            E + A <--> G + P
            G + A <--> E + P
            E + A <--> H + P
            H + A <--> E + P
            E + A <--> J + P
            J + A <--> E + P
            E + A <--> M + P
            M + A <--> E + P
        end
    end)
    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            A = [1.0, 2.0, 3.0, 4.0], P = [0.1, 0.2, 0.3, 0.4])
    prob = IdentifyRateEquationProblem(rxn, data; Keq = 10.0)
    id = EnzymeRates._mechanism_id(sprint(show, pinwheel))
    batch_log = String[]
    entries, failures = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[pinwheel], prob;
        optimizer = _testhelper_CountingStubOpt(), max_param_count = 20,
        eq_complexity_filter = 337, n_restarts = 1, maxtime = 1.0,
        log = msg -> push!(batch_log, msg))
    @test isempty(failures)
    @test only(entries).key == only(entries).row.key ==
          bytes2hex(EnzymeRates.sha256("mechanism|" * id))
    @test batch_log == ["1 new fits + 0 inherited + 0 skipped (already fit) + " *
                        "0 skipped (>20 params) + 0 skipped (>337 complexity)\n" *
                        "  leaf cap reached for $id; fitted under its own key"]
end
```

(o) Replace the whole testset "renaming-dup pair: same independent count, different
eq_hash" (with its three-line comment above it, "# LDH renaming-dup pair (same graph,
tied kinetic-group split): merged form" through "# form's extra binding K is tied straight
back by a Wegscheider cycle.") with:

```julia
# LDH renaming-dup pair (same graph, tied kinetic-group split): merged form (8 groups) vs
# split form (9 groups). The split form's extra binding K is tied straight back by a
# Wegscheider cycle, so the two forms have one law: one key, one fit.
@testset "renaming-dup pair: one key, one fit" begin
    # Merged, 8 kinetic groups.
    m1 = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: NADH, Pyruvate
        products:   Lactate, NAD
        steps: begin
            (E + Lactate ⇌ E(Lactate), E(NAD) + Lactate ⇌ E(Lactate, NAD),
             E(NADH) + Lactate ⇌ E(Lactate, NADH))
            (E + NAD ⇌ E(NAD), E(Lactate) + NAD ⇌ E(Lactate, NAD))
            (E + NADH ⇌ E(NADH), E(Lactate) + NADH ⇌ E(Lactate, NADH))
            E + Pyruvate ⇌ E(Pyruvate)
            (E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate),
             E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate))
            E(NADH, Pyruvate) <--> E(Lactate, NAD)
            E(Pyruvate) + NADH ⇌ E(NADH, Pyruvate)
            E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate)
        end
    end)
    # Split, 9 groups — same graph, Wegscheider-tied.
    m2 = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: NADH, Pyruvate
        products:   Lactate, NAD
        steps: begin
            (E + Lactate ⇌ E(Lactate), E(NAD) + Lactate ⇌ E(Lactate, NAD))
            (E + NAD ⇌ E(NAD), E(Lactate) + NAD ⇌ E(Lactate, NAD))
            (E + NADH ⇌ E(NADH), E(Lactate) + NADH ⇌ E(Lactate, NADH))
            E + Pyruvate ⇌ E(Pyruvate)
            (E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate),
             E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate))
            E(NADH) + Lactate ⇌ E(Lactate, NADH)
            E(NADH, Pyruvate) <--> E(Lactate, NAD)
            E(Pyruvate) + NADH ⇌ E(NADH, Pyruvate)
            E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate)
        end
    end)
    @test m1 != m2
    @test EnzymeRates._independent_param_count(m1) ==
          EnzymeRates._independent_param_count(m2)
    c1, c2 = (EnzymeRates._canonical(EnzymeRates._law(m),
                                     EnzymeRates._mechanism_id(sprint(show, m)))
              for m in (m1, m2))
    @test c1.key == c2.key
    @test _testhelper_same_law_under_map(EnzymeRates._law(m1), EnzymeRates._law(m2),
                                         c1.order, c2.order;
                                         rng = Random.MersenneTwister(2))
    # identify fits the pair once; the second form inherits the first's fit.
    rxn = @enzyme_reaction begin
        substrates: NADH[N1H1], Pyruvate[C1]
        products:   Lactate[C1H1], NAD[N1]
    end
    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            NADH = [1.0, 2.0, 1.0, 2.0], Pyruvate = [0.5, 0.5, 1.0, 1.0],
            Lactate = [0.1, 0.2, 0.1, 0.2], NAD = [0.3, 0.3, 0.4, 0.4])
    prob = IdentifyRateEquationProblem(rxn, data; Keq = 2.0)
    opt = _testhelper_CountingStubOpt(; uval = 0.0)
    entries, failures = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m1, m2], prob;
        optimizer = opt, max_param_count = 20, n_restarts = 1, maxtime = 1.0)
    @test opt.count == 1 && isempty(failures)
    @test [e.row.fit_inherited for e in entries] == [false, true]
    @test [e.key for e in entries] == [c1.key, c1.key]
end
```

(p) Testset "allosteric split/merge pair: the split form is the merged form's model":
replace the comment above it,

```julia
# LDH ALLOSTERIC split/merge pair: split form (8 cat groups) vs merged (7).
# Different eq_hash, though the cycles tie the split form's extra constant back, so
# both are one model.
```

with

```julia
# LDH ALLOSTERIC split/merge pair: split form (8 cat groups) vs merged (7). The cycles tie
# the split form's extra constant back, so both are one law, with one key.
```

and replace

```julia
    _testhelper_key(m) = EnzymeRates._rate_eq_dedup_key(
        rate_equation_string(EnzymeRates.compile_mechanism(m)))
    @test am1 != am2
    # the two render different equations
    @test _testhelper_key(am1) != _testhelper_key(am2)
```

with

```julia
    _testhelper_canon(m) = EnzymeRates._canonical(
        EnzymeRates._law(m), EnzymeRates._mechanism_id(sprint(show, m)))
    @test am1 != am2
    # one law under the canonical map, so one key
    c1, c2 = _testhelper_canon(am1), _testhelper_canon(am2)
    @test c1.key == c2.key
    @test _testhelper_same_law_under_map(EnzymeRates._law(am1), EnzymeRates._law(am2),
                                         c1.order, c2.order;
                                         rng = Random.MersenneTwister(3))
```

(q) Testset "_cv_model_selection: flatten reproduces serial LOOCV per candidate": replace
its comment's first sentence "Two candidates with distinct equations." with "Three
candidates with distinct keys, each the source of its own fit.", and replace

```julia
    function _testhelper_mkrow(m, loss)
        em = EnzymeRates.compile_mechanism(m)
        fkeys = EnzymeRates.fitted_params(em)
        (n_params=length(fkeys), loss=loss, mechanism=sprint(show, m),
         mechanism_id=EnzymeRates._mechanism_id(sprint(show, m)),
         rate_equation="v", retcode="Success", error=missing,
         params=NamedTuple{fkeys}(ntuple(_ -> 1.0, length(fkeys))),
         eq_hash=string(EnzymeRates._rate_eq_dedup_key(rate_equation_string(em)),
                        base=16, pad=16),
         fit_inherited=false)
    end
    # `cands` and the rows are parallel and go in against (n_params, loss) order; they
    # must come back sorted by it: the 6-parameter row has the lowest loss but sorts last.
    df = EnzymeRates._rows_to_dataframe([_testhelper_mkrow(m6, 0.1),
        _testhelper_mkrow(m1, 0.5), _testhelper_mkrow(m3, 0.2)])
    save_dir = mktempdir()
    stub() = _testhelper_CountingStubOpt(; uval=log(5.0))
    res = EnzymeRates._cv_model_selection(cands, df, prob;
        optimizer=stub(), se_threshold=1.0, save_dir, show_progress=false,
        n_restarts=1, maxtime=1.0)
    @test nrow(res.cv_results) == 3
    @test df.n_params == [6, 5, 5]
```

with

```julia
    # The pool as `_beam_search` hands it over: each candidate with its row and key, and
    # the memo holding each key's fit with the candidate as its source.
    memo = Dict{String, EnzymeRates._Fit}()
    function _testhelper_pool_entry(m, loss)
        em = EnzymeRates.compile_mechanism(m)
        fkeys = EnzymeRates.fitted_params(em)
        text = sprint(show, m)
        id = EnzymeRates._mechanism_id(text)
        key = EnzymeRates._canonical(EnzymeRates._law(m), id).key
        memo[key] = EnzymeRates._Fit(ones(length(fkeys)), loss, :Success, false, m, id)
        row = (n_params=length(fkeys), loss=loss, mechanism=text, mechanism_id=id,
               rate_equation="v", retcode="Success", error=missing,
               params=NamedTuple{fkeys}(ntuple(_ -> 1.0, length(fkeys))), key,
               fit_inherited=false, fit_source_id=id)
        EnzymeRates.BatchEntry(m, length(fkeys), loss, :Success, key, key, row)
    end
    # The pool goes in against (n_params, loss) order and must come back sorted by it:
    # the 6-parameter candidate has the lowest loss but sorts last.
    pool = [_testhelper_pool_entry(m6, 0.1), _testhelper_pool_entry(m1, 0.5),
            _testhelper_pool_entry(m3, 0.2)]
    save_dir = mktempdir()
    stub() = _testhelper_CountingStubOpt(; uval=log(5.0))
    res = EnzymeRates._cv_model_selection(pool, memo, prob;
        optimizer=stub(), se_threshold=1.0, save_dir, show_progress=false,
        n_restarts=1, maxtime=1.0)
    @test nrow(res.cv_results) == 3
    @test [e.n_params for e in pool] == [6, 5, 5]
```

Delete the now unused line `cands = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m6, m1, m3]`,
and in the loop replace `m = only(c for c in cands if sprint(show, c) == r.mechanism)`
with `m = only(e.mech for e in pool if sprint(show, e.mech) == r.mechanism)`. The serial
loop's `@test folds(r) == serial` is the value check: each candidate is the source of its
own fit, so its fold losses are its serial ones exactly.

(r) In `test/test_compile_budget.jl`, in `fit_script`'s `_testhelper_measure`, replace

```julia
            eq_text = rate_equation_string(em)
            ER._rate_eq_dedup_key(eq_text)
            ER._mechanism_id(sprint(show, m))
```

with

```julia
            rate_equation_string(em)
            id = ER._mechanism_id(sprint(show, m))
            s = ER._solve(m)
            law = ER._law(m, s)
            ER._assert_support(m, law)
            ER._canonical(law, id)
```

and in the comment above `_testhelper_measure` nothing changes: PASS1 still compiles,
fits and prints here (Task 12b removes the compile).

(s) Rename the section line above `_testhelper_CountingStubOpt`,

```julia
# ── §2 fit-dedup by eq_hash ──────────────────────────────────────────────────
```

to

```julia
# ── §2 fit-dedup by key ──────────────────────────────────────────────────────
```

(the same 79 columns), so no test line names `eq_hash` (Step 8's grep, Task 13's).

- [ ] **Step 3: Run the identify tests to verify they fail**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T8_tid_red.log 3600 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); @testset "red" begin include("test/test_identify_rate_equation.jl") end'`,
then `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; L=$S/phase45/logs/P5_T8_tid_red.log; timeout 580 bash -c "until grep -q '^exit ' $L; do sleep 20; done"; grep -E "Error|Fail|Test Summary" $L | head -60`
(repeat the poll while the log is open).
Expected: the converted and new testsets fail or error, among them
``UndefVarError: `_Fit` not defined in `EnzymeRates` ``, a `MethodError` for the
seven-argument `BatchEntry`, ``FieldError: type NamedTuple has no field `key` `` (the old
rows), and the
leaf-cap testset (`BatchEntry` has no field `key`); the testsets this step left alone
still pass, and the outer `red` testset ends the log with `exit 1`.

- [ ] **Step 4: Write the key cache**

In `src/identify_rate_equation.jl`:

Replace the three ABOUTME lines with

```julia
# ABOUTME: Beam-search pipeline to identify the best rate equation.
# ABOUTME: Enumerates mechanisms, fits one representative per rate-law key
# ABOUTME: (`_canonical`), and selects among the fits by cross-validation.
```

Delete `_rate_eq_dedup_key` with its docstring (from `"""` "Equation-identity key: …"
through the function's `end`).

Replace the `BatchEntry` docstring and struct (`"""One fitted mechanism: its own params +
retcode + eq_hash + the CSV row."""` through its `end`) with:

```julia
"""
A fit as identify stores it, once per rate-law key (`_canonical`): `values` holds the
fitted constants at their canonical positions (`_Canonical.order`), then `L` for an MWC
law, so every mechanism with that key reads them through its own canonical order;
`loss`, `retcode` and `peak_at_range_edge` are `fit_rate_equation`'s; `source` is the
mechanism that was fitted and `source_id` its `mechanism_id`.
"""
struct _Fit
    values::Vector{Float64}
    loss::Float64
    retcode::Symbol
    peak_at_range_edge::Bool
    source::Union{Mechanism, AllostericMechanism}
    source_id::String
end

"""
A mechanism PASS1 of `_process_batch` derived: its fitted constants in `fitted_params`
order, its text and the text's id (`_mechanism_id`), its printed law, and its law's
canonical form (`canon`), under whose key it is fitted.
"""
struct _Candidate
    mech::Union{Mechanism, AllostericMechanism}
    fitted::Vector{Symbol}
    mechanism::String
    mechanism_id::String
    eq_text::String
    canon::_Canonical
end

"""The order in which candidate `c` reads its fit's values (`_Fit`): its canonical order,
then `L` for an MWC law."""
_fit_order(c::_Candidate) =
    c.mech isa AllostericMechanism ? [c.canon.order; :L] : c.canon.order

"""
The parameters of candidate `c` as its row holds them, in its `fitted_params` order: a
fit's `values`, stored at canonical positions, read through `order`, the candidate's own
order for those positions.
"""
function _row_params(c::_Candidate, order::Vector{Symbol}, values::Vector{Float64})
    at = Dict(zip(order, values))
    NamedTuple{Tuple(c.fitted)}(Tuple(at[x] for x in c.fitted))
end

"""One fitted mechanism: its parameter count, its fit's loss and retcode, its rate-law key
(`_canonical`), the key of the fit it reads (`fit_key`), and its CSV row."""
struct BatchEntry
    mech::Union{Mechanism, AllostericMechanism}
    n_params::Int
    loss::Float64
    retcode::Symbol
    key::String
    fit_key::String
    row::NamedTuple
end
```

In `_rows_to_dataframe`, replace the line `eq_hash = [r.eq_hash for r in rows],` and the
`fit_inherited` line after it with

```julia
        key = Union{Missing,String}[r.key for r in rows],
        fit_inherited = Union{Missing,Bool}[r.fit_inherited for r in rows],
        fit_source_id = Union{Missing,String}[r.fit_source_id for r in rows],
```

In `_failure_row`, replace

```julia
     eq_hash = missing,
     fit_inherited = missing,
```

with

```julia
     key = missing,
     fit_inherited = missing,
     fit_source_id = missing,
```

In `_prefit_summary`'s docstring, replace "`new fits` (distinct new equations being fit
this batch)" with "`new fits` (distinct new laws, one per key, being fit this batch)".

Replace the whole `_process_batch` (docstring and function) with:

```julia
"""
Fit one batch of mechanisms. Returns `(entries, failures, n_param_skip, n_cx_skip,
n_seen_skip)`.

PASS-1: skip mechanisms already produced, in an earlier batch or earlier in this one (a
mechanism whose structural `hash` is already in `seen` is dropped before compiling;
`seen` gets every structure on first sight regardless of outcome). For each fresh
mechanism, in parallel (`pmap`): skip it over `eq_complexity_filter` or
`max_param_count`; else compile it and print its law, take its text and the text's id
(`_mechanism_id`), derive its law from the concrete mechanism (`_solve`, `_law`), check
the law on the guard's faces (`_assert_support`) and label it canonically
(`_canonical`). A law whose labelling reaches the leaf cap gets a key of its own, and the
pre-fit summary names its mechanism. Then pick one fit representative per key not
already in `memo`, and pass the pre-fit summary (`_prefit_summary`) to `log`. `mechs` is
already structurally deduped by the caller (`unique!`).

PASS-2: fit ONE representative per key not already in `memo`, in parallel across all
workers (`pmap`), and store the fit in `memo` at canonical positions (`_Fit`). Every
mechanism with that key reads it through its own canonical order, so its row carries its
own parameter names: same law ⟹ same fitted values and same rescaling, so no refit is
needed. `memo` persists across iterations (threaded from `_beam_search`), so a law fit in
an earlier iteration is never refit. A representative whose fit throws fails ALL the
mechanisms with its key (all-or-nothing per law). A row's `fit_inherited` is `false`
when its mechanism is its fit's source and `true` otherwise, and `fit_source_id` names
the source. `parent_of` maps a child to its parent's `BatchEntry`, which fills the row's
parent columns. `entries::Vector{BatchEntry}` are the mechanisms with a usable fit (each
keeping its own row, `retcode` and key); `failures::Vector{FitFailure}` are mechanisms
that threw at compile, derivation or fit — captured WITH the exception text, never
silently swallowed.
"""
function _process_batch(
    mechs, prob::IdentifyRateEquationProblem;
    optimizer, max_param_count, eq_complexity_filter::Int = typemax(Int),
    memo::Dict{String,_Fit}=Dict{String,_Fit}(),
    seen::Set{UInt64}=Set{UInt64}(),
    parent_of::AbstractDict = Dict(), log = msg -> nothing, kwargs...
)
    # Skip structures already produced — expand each once. Added to `seen` on first
    # sight regardless of outcome (fit / cap / error).
    fresh = empty(mechs)
    for m in mechs
        h = hash(m)
        h in seen || (push!(seen, h); push!(fresh, m))
    end

    # PASS 1 (workers): complexity filter (checked first, before any derivation),
    # compile, parameter cap, printed law, then the concrete mechanism's law, its check on
    # the guard's faces and its canonical form. `:complexity_skip` and `:param_skip` are
    # skips; `FitFailure` = threw; else a `_Candidate`.
    compiled = pmap(fresh) do m
        try
            _eq_complexity(m) > eq_complexity_filter && return :complexity_skip
            em = compile_mechanism(m)
            length(fitted_params(em)) > max_param_count && return :param_skip
            eq_text = rate_equation_string(em)
            text = sprint(show, m)
            id = _mechanism_id(text)
            s = _solve(m)
            L = _law(m, s)
            _assert_support(m, L)
            _Candidate(m, s.fitted, text, id, eq_text, _canonical(L, id))
        catch e
            FitFailure(m, _exc_string(e))
        end
    end

    # Pick one representative per key not already fit.
    rep_idx = Dict{String,Int}()
    for (i, c) in enumerate(compiled)
        c isa _Candidate && !haskey(memo, c.canon.key) && get!(rep_idx, c.canon.key, i)
    end
    reps = _Candidate[compiled[i] for i in values(rep_idx)]
    n_param_skip = count(==(:param_skip), compiled)
    n_cx_skip = count(==(:complexity_skip), compiled)
    n_seen_skip = length(mechs) - length(fresh)
    capped = ["leaf cap reached for $(c.mechanism_id); fitted under its own key"
              for c in compiled if c isa _Candidate && c.canon.capped]
    log(join([_prefit_summary(length(reps),
                  count(c -> c isa _Candidate, compiled) - length(reps),
                  n_param_skip, n_cx_skip, n_seen_skip;
                  max_param_count, eq_complexity_filter); capped], "\n  "))

    # PASS 2: fit them in PARALLEL (`pmap`) — fitting dominates cost, so it must run
    # across all workers. Each representative gives its fit at canonical positions, or
    # the text of the error its fit threw.
    rep_fits = pmap(reps) do r
        try
            fp = FittingProblem(compile_mechanism(r.mech), prob.data;
                Keq=prob.Keq, absolute=prob.absolute, peak_range=prob.peak_range)
            fit = fit_rate_equation(fp, optimizer; kwargs...)
            _Fit(Float64[fit.params[x] for x in _fit_order(r)], fit.loss, fit.retcode,
                 fit.peak_at_range_edge, r.mech, r.mechanism_id)
        catch e
            _exc_string(e)
        end
    end
    fit_error = Dict{String,String}()   # key → its representative's fit threw
    for (r, f) in zip(reps, rep_fits)
        f isa _Fit ? (memo[r.canon.key] = f) : (fit_error[r.canon.key] = f::String)
    end

    # Build a row per derived mechanism from its key's fit, read through its own order.
    entries  = BatchEntry[]
    failures = FitFailure[]
    for c in compiled
        c isa FitFailure && (push!(failures, c); continue)
        c isa _Candidate || continue                       # a skip
        key = c.canon.key
        if haskey(fit_error, key)                          # representative fit threw
            push!(failures, FitFailure(c.mech, fit_error[key]))
            continue
        end
        fit = memo[key]
        parent = get(parent_of, c.mech, nothing)
        row = (
            n_params = length(c.fitted),
            parent_n_params = parent === nothing ? missing : parent.n_params,
            loss = fit.loss,
            mechanism = c.mechanism,
            mechanism_id = c.mechanism_id,
            parent_id = parent === nothing ? missing : parent.row.mechanism_id,
            rate_equation = c.eq_text,
            retcode = string(fit.retcode),
            error = missing,
            params = _row_params(c, _fit_order(c), fit.values),
            key,
            fit_inherited = fit.source_id != c.mechanism_id,
            fit_source_id = fit.source_id,
            peak_at_range_edge = fit.peak_at_range_edge,
        )
        push!(entries, BatchEntry(c.mech, length(c.fitted), fit.loss, fit.retcode, key,
                                  key, row))
    end
    (entries, failures, n_param_skip, n_cx_skip, n_seen_skip)
end
```

In `_ingest!`'s docstring, replace

```
joins the `frontier` (the unexpanded work queue — ALL structurally-distinct
mechanisms, no eq-dedup); `best_loss_by_count` tracks the per-count running
min (the beam-cutoff reference); `cv_pool` keeps the top `n_cv_candidates`
DISTINCT equations (by `eq_hash`, lowest loss each) per param count. An entry whose
```

with

```
joins the `frontier` (the unexpanded work queue — ALL structurally-distinct
mechanisms, no key dedup); `best_loss_by_count` tracks the per-count running
min (the beam-cutoff reference); `cv_pool` keeps the top `n_cv_candidates`
entries with DISTINCT rate-law keys (lowest loss each) per param count. An entry whose
```

Replace `_offer_cv!`'s docstring and its first line,

```julia
"""
Keep `pool` at the top `n` distinct-`eq_hash` entries by loss. A repeat
`eq_hash` only ever updates its own slot (to the lower loss); it never
consumes a second slot.
"""
function _offer_cv!(pool::Vector{BatchEntry}, e::BatchEntry, n::Int)
    idx = findfirst(p -> p.eq_hash == e.eq_hash, pool)
```

with

```julia
"""
Keep `pool` at the top `n` entries with distinct rate-law keys, by loss. A repeat key
only ever updates its own slot (to the lower loss); it never consumes a second slot.
"""
function _offer_cv!(pool::Vector{BatchEntry}, e::BatchEntry, n::Int)
    idx = findfirst(p -> p.key == e.key, pool)
```

In `_beam_search`, replace

```julia
    # Fits as `fit_rate_equation` returns them (scaled to peak turnover 1 in relative
    # mode), keyed by `eq_hash`, shared across iterations so each distinct equation is
    # fit exactly once over the whole search.
    memo = Dict{UInt64,NamedTuple}()
```

with

```julia
    # Fits by rate-law key, at canonical positions (`_Fit`; scaled to peak turnover 1 in
    # relative mode), shared across iterations so each distinct law is fit exactly once
    # over the whole search. LOOCV refits each key's source on its folds.
    memo = Dict{String,_Fit}()
```

replace

```julia
    base_batch === nothing && return (
        Union{Mechanism, AllostericMechanism}[], _rows_to_dataframe(NamedTuple[]))
```

with `base_batch === nothing && return (BatchEntry[], memo)`, and replace

```julia
    pool_entries = BatchEntry[e for v in values(cv_pool) for e in v]
    mechs = Union{Mechanism, AllostericMechanism}[
        e.mech for e in pool_entries]
    df = _rows_to_dataframe([e.row for e in pool_entries])
    return mechs, df
```

with

```julia
    return BatchEntry[e for v in values(cv_pool) for e in v], memo
```

In `identify_rate_equation`, replace `mechanisms, df = _beam_search(prob;` with
`pool, memo = _beam_search(prob;` and

```julia
    result = _cv_model_selection(
        mechanisms, df, prob;
```

with

```julia
    result = _cv_model_selection(
        pool, memo, prob;
```

and in its docstring replace

```
- `n_cv_candidates::Int = 5`: LOOCV top N
  **unique-rate-equation** candidates per param count (at least 1)
```

with

```
- `n_cv_candidates::Int = 5`: LOOCV top N candidates with **distinct rate-law
  keys** per param count (at least 1)
```

Replace `_cv_model_selection` from its `function` line through the line
`best_mechanism = compile_mechanism(candidate_mechs[best_row_idx])` with:

```julia
"""
    _cv_model_selection(pool, memo, prob; optimizer, se_threshold, save_dir,
                        show_progress, kwargs...) -> IdentifyRateEquationResults

Cross-validate the candidates of `pool` (the cv pool `_beam_search` returns) by leaving
one group of `prob.data` out at a time, select one by the 1-SE rule (`_select_best_row`)
and write `loocv_results.csv` and `best_equation.csv` to `save_dir`. A candidate's fold
losses come from fitting, on each fold, the source of the fit it reads
(`memo[e.fit_key].source`), so candidates that read one fit share their fold losses and
each fit's folds run once.
"""
function _cv_model_selection(
    pool::Vector{BatchEntry}, memo::Dict{String,_Fit},
    prob::IdentifyRateEquationProblem;
    optimizer, se_threshold::Float64,
    save_dir, show_progress,
    kwargs...
)
    isempty(pool) && error(
        "No mechanisms were successfully " *
        "fitted during beam search")

    # LOOCV candidates: the cv pool, ordered by `n_params`, then `loss`. `_offer_cv!`
    # keeps the top `n_cv_candidates` entries with DISTINCT rate-law keys (lowest `loss`
    # each) per `n_params` bucket, so one law never takes two LOOCV slots at a count and
    # folds are never wasted on, or biased by, its duplicates.
    candidates = pool[sortperm([(e.n_params, e.loss) for e in pool])]
    cv_df = _rows_to_dataframe([e.row for e in candidates])

    _progress(save_dir, show_progress,
        "Cross-validating $(length(candidates)) candidate equations (LOOCV)…")
    # The fold losses of each fit key, from fitting its source once per fold, flattened to
    # a (fit key, fold) grid so all folds run across every worker. `pmap` keeps the grid
    # order (the fold varies fastest), so column j of `by_fit` holds fit key j's fold
    # losses in `groups` order, and column ci of `scores` candidate ci's.
    groups = unique(prob.data.group)
    fit_keys = unique(e.fit_key for e in candidates)
    fold_losses = pmap([(memo[k].source, g) for k in fit_keys for g in groups]) do (m, g)
        _cv_fold_loss(compile_mechanism(m), prob, g; optimizer, kwargs...)
    end
    by_fit = reshape(Float64.(fold_losses), length(groups), :)
    scores = by_fit[:, [findfirst(==(e.fit_key), fit_keys) for e in candidates]]

    cv_df.cv_score = [mean(c) for c in eachcol(scores)]
    cv_df.cv_score_se = [std(c) / sqrt(length(groups)) for c in eachcol(scores)]

    best_row_idx = _select_best_row(cv_df; se_threshold)
    best_mechanism = compile_mechanism(candidates[best_row_idx].mech)
```

and, further down, replace `"(eq_hash=$(cv_df.eq_hash[best_row_idx])), " *` with
`"(key=$(cv_df.key[best_row_idx])), " *`.

In `src/rate_eq_derivation.jl`, in `_equation_text`'s docstring, replace

```
`# Wegscheider constraints:`. Each section is sorted by name, which is load-bearing: the
eq_hash dedup of rate-equivalent mechanisms compares these strings, so their line order
must not depend on the order the solve emits its dependents in. It first runs
```

with

```
`# Wegscheider constraints:`. Each section is sorted by name, so a law's text does not
depend on the order the solve emits its dependents in. It first runs
```

- [ ] **Step 5: Run the identify tests to verify they pass**

Run (background) Step 3's command with log `$S/phase45/logs/P5_T8_tid.log`, polled the
same way.
Expected: every testset passes, pristine, `exit 0`. Check the counts the gate names in
the log's summary lines: "two same-key mechanisms read one fit under their own names"
and "a law at the leaf cap is fitted under its mechanism's key" pass all their tests.

- [ ] **Step 6: Measure the compile budgets**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T8_budget.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_compile_budget.jl")'`
Expected: `exit 0`, every budget testset passes; record the `[ Info: fit compile per
fresh mechanism: method instances …` and `[ Info: show instances per fresh mechanism: …`
lines against Task 0's `figures.md`. The key computes in code that does not specialize
on the mechanism, so the per-mechanism instance counts rise by little. A budget that
fails is not raised: stop and report the readings to the controller (CLAUDE.md).

- [ ] **Step 7: Update the ROADMAP**

In `ROADMAP.md`, under "Known issues", delete the bullet that starts
"- **Splitting a tied kinetic group changes `eq_hash` without changing the model.**"
through its last line, "fits both. The rename-proof key under Planned removes this."

Under "Planned", in "**Profile the first v0.9.0 HPC runs.**", replace

```
  the final parsimony cutoff, children and survivors per move, same-function equations
  that `eq_hash` misses (`_rate_eq_dedup_key`, `src/identify_rate_equation.jl`),
```

with

```
  the final parsimony cutoff, children and survivors per move, same-function laws that
  the key misses, equal only after a reparametrization (`_canonical`, `src/key.jl`),
```

Under "Ideas", in "**Seed each identify fit from its mechanism's content.**", replace

```
  (`src/identify_rate_equation.jl`) fits a batch's representatives in the iteration order
  of a `Dict` keyed by `eq_hash`, filled in the order the moves emit the mechanisms. So a
  change of emission order or of a law's printed text moves the starting points of every
```

with

```
  (`src/identify_rate_equation.jl`) fits a batch's representatives in the iteration order
  of a `Dict` keyed by the rate-law key (`_canonical`), filled in the order the moves emit
  the mechanisms. So a change of emission order or of a law's key moves the starting
  points of every
```

and rewrap that paragraph to 92 columns without changing its other words.

- [ ] **Step 8: Run the full suite**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T8_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`,
polled in bounded loops (`timeout 580 bash -c 'until grep -q "^exit " LOG; do sleep 20; done'; tail -40 LOG`).
Expected: `exit 0`, every test passing, pristine, Aqua and JET green; about 8 s over Task
7's time; report `maxrss_kb` against Task 0's. `grep -rn "eq_hash\|_rate_eq_dedup_key" src test`
prints nothing (the docs still name `eq_hash` until Task 15 rewrites the identify pages).

- [ ] **Step 9: Commit**

Write the measured figures in place of `<…>`.

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt
git add src/identify_rate_equation.jl src/rate_eq_derivation.jl ROADMAP.md \
    test/test_identify_rate_equation.jl test/test_compile_budget.jl
git commit -m "$(cat <<'EOF'
Fit one law per key and read each fit through its own names

PASS1 derives the law of the concrete mechanism and labels it canonically.
The memo holds one fit per key at canonical positions, and every row reads
it through its own order, so a same-key mechanism keeps its own parameter
names. The CSV gains key and fit_source_id; eq_hash goes. The cv pool keeps
one slot per key, and LOOCV fits each fit key's source once per fold. A law
at the leaf cap is fitted under its mechanism's key, and the log names it.

The dedup pins moved after their value checks: each same-key pair has one
key and one law under the canonical map, and a swapped map fails. A row's
fit_inherited means the row is not its fit's source, so the source read
again in a later batch is not inherited. The _rate_eq_dedup_key testset's
subject is gone; it became the key's determinism test.

Fit compile per fresh mechanism: <instances> instances, <kb> KB (budgets
250, 300 KB). Full-suite maxrss: <kb> KB (Task 0: <kb0> KB).

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
)"
```

---
### Task 9: Beam per key (P5-29)

The beam ranks rate-law keys, not mechanisms. The mechanisms of one key share one fit, so
a key takes one rank and one unit of the cumulative budget, and duplicates cannot crowd
out a distinct key. Every member of a selected key expands, and a key keeps its verdict
for the rest of the search: a late member of a selected key expands without budget, a
late member of a rejected key stays out.

**Files:**
- Modify: `src/identify_rate_equation.jl` (`_select_count!`; `_beam_search` keeps the
  `selected_keys` and `rejected_keys` sets and passes the keys; the `seen` comment;
  `identify_rate_equation`'s `min_beam_width` and "Beam selection" docstring text)
- Test: `test/test_identify_rate_equation.jl`
- Create (scratch): `$S/phase45/scripts/P5_T9_check.jl`
- Run (scratch): `$S/phase45/scripts/P5_identify_run.jl $S/phase45/p5_T9_identify` (Task 0)
- Logs: `$S/phase45/logs/P5_T9_{tid_red,tid,identify,check,full}.log`

**Interfaces:**
- Consumes: Task 8's `BatchEntry.key`, `_beam_search`, `memo`; Task 7's `_canonical`;
  `_parse_mechanism`; Task 0's `P5_identify_run.jl`, `$S/phase45/p5_T0_identify/summary.txt`
  and that run's `best_equation.csv`.
- Produces:
  ```julia
  _select_count!(expanded::Dict{Int, Int}, best_loss_by_count::Dict{Int, Float64},
                 c::Int, losses::AbstractVector{<:Real},
                 law_keys::AbstractVector{String};
                 loss_rel_threshold::Float64, loss_abs_threshold::Float64,
                 loss_parsimony_threshold::Float64, min_beam_width::Int,
                 selected::Set{String}, rejected::Set{String})::Vector{Int}
  ```
  It returns, in input order, the indices of every mechanism whose key is selected,
  now or in an earlier sweep. `expanded[c]` counts the keys selected at count `c`.

**Gate:**
- `test_identify_rate_equation.jl` is green; the three converted `_select_count!`
  testsets keep every expected index.
- The identify run selects a mechanism whose key equals the true mechanism's key, or,
  when Task 0's run did not select the truth, the key of Task 0's selection (Task 0
  Step 10) or a near-exact law (`cv_score` below 1e-6) with at most Task 0's selected
  `n_params`, since on these noise-free data the choice among near-exact laws follows
  the optimizer's random stream (Open for Denis, item 34); its `loocv_results.csv` holds
  one row per key per count; the fits run, the selected `n_params` and `cv_score`, and
  whether each run selected the truth are reported beside Task 0's `summary.txt`. Any
  other selection is investigated, not re-pinned.
- The full suite is green and pristine; about +1 s; the identify run takes about 2 min
  (Task 0's run on f9361cf: 102 s).

- [ ] **Step 1: Write the failing tests**

In `test/test_identify_rate_equation.jl`, testset "_select_count!: thresholds, floor, best
loss, parsimony", replace the helper

```julia
    _testhelper_selected(losses, best; rel, add = 0.0, width = 1,
                         others = Dict{Int,Float64}(), parsimony = 1.0) =
        EnzymeRates._select_count!(Dict{Int,Int}(), merge(Dict(5 => best), others), 5,
            losses; loss_rel_threshold = rel, loss_abs_threshold = add,
            loss_parsimony_threshold = parsimony, min_beam_width = width)
```

with

```julia
    # Each mechanism its own key, so keys rank as the mechanisms did.
    _testhelper_selected(losses, best; rel, add = 0.0, width = 1,
                         others = Dict{Int,Float64}(), parsimony = 1.0) =
        EnzymeRates._select_count!(Dict{Int,Int}(), merge(Dict(5 => best), others), 5,
            losses, string.(eachindex(losses)); loss_rel_threshold = rel,
            loss_abs_threshold = add, loss_parsimony_threshold = parsimony,
            min_beam_width = width, selected = Set{String}(), rejected = Set{String}())
```

In testset "_select_count! cumulative per-count floor", replace

```julia
    kw = (loss_rel_threshold=1.0, loss_abs_threshold=0.0,
          loss_parsimony_threshold=1.0, min_beam_width=3)
```

with

```julia
    # Each mechanism its own key; the key sets persist across the sweeps, as in the beam.
    kw = (loss_rel_threshold=1.0, loss_abs_threshold=0.0,
          loss_parsimony_threshold=1.0, min_beam_width=3,
          selected=Set{String}(), rejected=Set{String}())
```

and add a key vector to each of its three calls:
`[1.0, 2.0, 3.0, 4.0, 5.0]; kw...)` becomes `[1.0, 2.0, 3.0, 4.0, 5.0],
["a", "b", "c", "d", "e"]; kw...)`, `[10.0, 11.0, 12.0]; kw...)` becomes `[10.0, 11.0,
12.0], ["f", "g", "h"]; kw...)`, and `[1.0, 20.0]; kw...)` becomes `[1.0, 20.0],
["i", "j"]; kw...)`.

In testset "§1 parsimony cutoff = threshold * min over all counts < c", replace

```julia
    _testhelper_selected(best_loss_by_count, c, losses) = EnzymeRates._select_count!(
        Dict{Int,Int}(), best_loss_by_count, c, losses; loss_rel_threshold=10.0,
        loss_abs_threshold=0.0, loss_parsimony_threshold=1.01, min_beam_width=0)
```

with

```julia
    _testhelper_selected(best_loss_by_count, c, losses) = EnzymeRates._select_count!(
        Dict{Int,Int}(), best_loss_by_count, c, losses, string.(eachindex(losses));
        loss_rel_threshold=10.0, loss_abs_threshold=0.0, loss_parsimony_threshold=1.01,
        min_beam_width=0, selected=Set{String}(), rejected=Set{String}())
```

These three testsets keep every expected index: with one key per mechanism the key ranks
are the mechanism ranks, so their pinned selections are the head's own values.

After the testset "§1 parsimony cutoff = threshold * min over all counts < c", add:

```julia
@testset "_select_count! charges the budget once per key" begin
    # Spec section 10: the beam counts each key once. Five mechanisms on two keys, the
    # three of key a sharing their fit's loss. With a budget of two, ranking mechanisms
    # would spend it on two of a's three and leave b out; ranking keys spends one unit on
    # each key and expands all five mechanisms.
    expanded = Dict{Int,Int}()
    selected, rejected = Set{String}(), Set{String}()
    sel = EnzymeRates._select_count!(expanded, Dict(5 => 1.0), 5,
        [1.0, 1.0, 1.0, 2.0, 2.0], ["a", "a", "a", "b", "b"];
        loss_rel_threshold = 1.0, loss_abs_threshold = 0.0, loss_parsimony_threshold = 1.0,
        min_beam_width = 2, selected, rejected)
    @test sel == [1, 2, 3, 4, 5]
    @test expanded[5] == 2
    @test selected == Set(["a", "b"]) && isempty(rejected)
end

@testset "a late member of a selected key expands without budget" begin
    # Key a is selected in the first sweep at loss 2.0, then the count's best. A later
    # sweep brings a structurally new member of a after another key lowered the best to
    # 1.0: the cutoff (1.0) and the spent budget would refuse a mechanism at 2.0, but a's
    # verdict stands, so the member expands and the budget stays spent on a alone.
    expanded = Dict{Int,Int}()
    selected, rejected = Set{String}(), Set{String}()
    kw = (loss_rel_threshold = 1.0, loss_abs_threshold = 0.0,
          loss_parsimony_threshold = 1.0, min_beam_width = 1, selected, rejected)
    @test EnzymeRates._select_count!(expanded, Dict(5 => 2.0), 5, [2.0], ["a"];
                                     kw...) == [1]
    @test expanded[5] == 1
    @test EnzymeRates._select_count!(expanded, Dict(5 => 1.0), 5, [2.0, 3.0],
                                     ["a", "c"]; kw...) == [1]
    @test expanded[5] == 1
    @test selected == Set(["a"]) && rejected == Set(["c"])
end

@testset "a late member of a rejected key stays out" begin
    # Key b loses on rank in the first sweep and is rejected. A later member of b is not
    # ranked again and stays out, while a later member of the selected key a expands.
    expanded = Dict{Int,Int}()
    selected, rejected = Set{String}(), Set{String}()
    kw = (loss_rel_threshold = 1.0, loss_abs_threshold = 0.0,
          loss_parsimony_threshold = 1.0, min_beam_width = 1, selected, rejected)
    @test EnzymeRates._select_count!(expanded, Dict(5 => 1.0), 5, [1.0, 1.5],
                                     ["a", "b"]; kw...) == [1]
    @test rejected == Set(["b"])
    @test EnzymeRates._select_count!(expanded, Dict(5 => 1.0), 5, [1.5, 1.0],
                                     ["b", "a"]; kw...) == [2]
    @test expanded[5] == 1
    @test selected == Set(["a"]) && rejected == Set(["b"])
end
```

- [ ] **Step 2: Run the identify tests to verify they fail**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T9_tid_red.log 3600 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); @testset "red" begin include("test/test_identify_rate_equation.jl") end'`,
then `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; L=$S/phase45/logs/P5_T9_tid_red.log; timeout 580 bash -c "until grep -q '^exit ' $L; do sleep 20; done"; grep -E "Error|Fail|Test Summary" $L | head -40`
(repeat while the log is open).
Expected: the six `_select_count!` testsets error with `MethodError: no method matching
_select_count!(::Dict{Int64, Int64}, ::Dict{Int64, Float64}, ::Int64, ::Vector{Float64},
::Vector{String}; …)`; everything else passes, and the outer `red` testset ends the log
with `exit 1`.

- [ ] **Step 3: Write the selection per key**

In `src/identify_rate_equation.jl`, replace `_select_count!` (docstring and function)
with:

```julia
"""
Select this sweep's parents among the count-`c` mechanisms with `losses` and rate-law
keys `law_keys`, and advance the cumulative floor budget. The beam ranks keys, not
mechanisms: the mechanisms of one key share one fit (`_process_batch`), so a key enters
at its lowest loss, takes one rank and spends one unit of budget, and duplicates cannot
crowd out a distinct key. A key this sweep meets for the first time is selected if
either:
  • its loss ≤ cutoff, where
    cutoff = min(loss_rel_threshold * best_loss_by_count[c] + loss_abs_threshold,
                 loss_parsimony_threshold * best(<c))
    and best(<c) is the best loss over ALL counts strictly below c (not just c-1): an
    added parameter must beat the best simpler model of any size. The parsimony term
    only tightens the cutoff, and is dropped while no smaller count has been fit,
  • OR its rank (1-indexed by ascending loss among this sweep's new keys) is within the
    *cumulative* floor budget. `expanded[c]` counts the keys at count `c` the whole
    search has selected so far; the width floor may add at most
    `min_beam_width - expanded[c]` more. Once the budget is spent, only the loss cutoff
    admits at that count.
A new key joins `selected` or `rejected` and keeps that verdict for the rest of the
search: a later member of a selected key is returned without spending budget, and a
later member of a rejected key is not. Returns the indices, in input order, of every
mechanism whose key is selected.

Every loss is finite: a failed or non-converging fit (`Inf`, `NaN`) should not
propagate to the next level, so `_ingest!` keeps it out of the frontier.
"""
function _select_count!(
    expanded::Dict{Int,Int}, best_loss_by_count::Dict{Int,Float64}, c::Int,
    losses::AbstractVector{<:Real}, law_keys::AbstractVector{String};
    loss_rel_threshold::Float64, loss_abs_threshold::Float64,
    loss_parsimony_threshold::Float64, min_beam_width::Int,
    selected::Set{String}, rejected::Set{String},
)
    cutoff = loss_rel_threshold * best_loss_by_count[c] + loss_abs_threshold
    smaller = [l for (k, l) in best_loss_by_count if k < c]
    isempty(smaller) ||
        (cutoff = min(cutoff, loss_parsimony_threshold * minimum(smaller)))
    budget = max(0, min_beam_width - get(expanded, c, 0))
    # The keys this sweep meets for the first time, in input order, each at its lowest
    # loss.
    fresh, lowest = String[], Dict{String,Float64}()
    for (l, k) in zip(losses, law_keys)
        (k in selected || k in rejected) && continue
        haskey(lowest, k) || push!(fresh, k)
        lowest[k] = min(get(lowest, k, Inf), l)
    end
    for (rank, j) in enumerate(sortperm([lowest[k] for k in fresh]))
        k = fresh[j]
        push!(lowest[k] <= cutoff || rank <= budget ? selected : rejected, k)
    end
    expanded[c] = get(expanded, c, 0) + count(in(selected), fresh)
    # Indices in original (input) order, so callers do not rely on the by-loss sort
    # order, a side-effect of the rank computation rather than part of the contract.
    findall(in(selected), law_keys)
end
```

In `_beam_search`, replace

```julia
    # Mechanisms expanded so far per parameter count — the cumulative floor
    # budget. Spent once over the whole search, never re-granted per sweep.
    expanded_by_count = Dict{Int, Int}()
```

with

```julia
    # Keys selected so far per parameter count — the cumulative floor budget. Spent once
    # over the whole search, never re-granted per sweep.
    expanded_by_count = Dict{Int, Int}()
    # Each key's verdict (`_select_count!`), kept for the whole search: a late member of
    # a selected key expands without budget, a late member of a rejected key stays out.
    selected_keys, rejected_keys = Set{String}(), Set{String}()
```

replace the `seen` comment

```julia
    # Structures already produced — each is expanded at most once (termination).
    # This preserves the selected model: expansion is deterministic, so a
    # re-produced structure yields identical children (no new reachable
    # mechanism), and its FIRST production is its best beam-selection chance —
    # the width-floor budget (`expanded_by_count`) only shrinks and the loss
    # cutoff (`best_loss_by_count` running min) only tightens, so anything the
    # old code selected on a later re-production it would already have selected
    # on the first. If `_select_count!` ever gains a non-monotone budget, this
    # invariant breaks.
```

with

```julia
    # Structures already produced — each is expanded at most once (termination).
    # This preserves the selected model: expansion is deterministic, so a
    # re-produced structure yields identical children (no new reachable
    # mechanism), and its FIRST production is its best beam-selection chance —
    # the width-floor budget (`expanded_by_count`) only shrinks, the loss cutoff
    # (`best_loss_by_count` running min) only tightens, and a key keeps its first
    # verdict (`selected_keys`, `rejected_keys`), so anything a later
    # re-production would select was selected on the first. If `_select_count!`
    # ever gains a non-monotone budget, this invariant breaks.
```

and replace

```julia
            sel = _select_count!(expanded_by_count, best_loss_by_count, c,
                [e.loss for e in entries_at_count];
                loss_rel_threshold, loss_abs_threshold, loss_parsimony_threshold,
                min_beam_width)
```

with

```julia
            sel = _select_count!(expanded_by_count, best_loss_by_count, c,
                [e.loss for e in entries_at_count], [e.key for e in entries_at_count];
                loss_rel_threshold, loss_abs_threshold, loss_parsimony_threshold,
                min_beam_width, selected = selected_keys, rejected = rejected_keys)
```

In `identify_rate_equation`'s docstring, replace

```
- `min_beam_width::Int = 50`: minimum mechanisms
  to keep per param-count tier
```

with

```
- `min_beam_width::Int = 50`: minimum rate-law keys to expand per
  param-count tier; every mechanism of a selected key expands
```

and in its "# Beam selection" section replace

```
A mechanism at parameter count `n` qualifies for the next-level
beam if either:
```

with

```
The beam selects rate-law keys: the mechanisms of one key share one fit,
so a key counts once in the ranks and the budget, and every mechanism of a
selected key expands. A key at parameter count `n` qualifies for the
next-level beam if either:
```

and

```
- OR it falls under the `min_beam_width` budget: a cumulative
  per-count allowance that expands at least `min_beam_width`
  mechanisms at each count over the whole search, then stops —
  spent once, not re-granted each time the count is revisited.
```

with

```
- OR it falls under the `min_beam_width` budget: a cumulative
  per-count allowance that expands at least `min_beam_width`
  keys at each count over the whole search, then stops —
  spent once, not re-granted each time the count is revisited.
  A key keeps its verdict: a mechanism found later with a key
  already selected expands, one with a key already rejected does not.
```

- [ ] **Step 4: Run the identify tests to verify they pass**

Run (background) Step 2's command with log `$S/phase45/logs/P5_T9_tid.log`.
Expected: every testset passes, pristine, `exit 0`.

- [ ] **Step 5: Run the seeded identify search**

Create `$S/phase45/scripts/P5_T9_check.jl`:

```julia
# ABOUTME: Checks a seeded identify run (P5_identify_run.jl) against its true mechanism: the
# ABOUTME: selected law's key, one LOOCV row per key and count, and the fits run.
# Usage, from the worktree root: julia --project=. P5_T9_check.jl <run_dir> <base_summary>
using TestEnv; TestEnv.activate()
using EnzymeRates, CSV, DataFrames
const ER = EnzymeRates
const DIR, BASE = ARGS[1], ARGS[2]                  # BASE: Task 0's summary.txt
# The true mechanism of P5_identify_run.jl, written out the same way.
truth = ER.Mechanism(@enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ E(A)
        E(A) + B ⇌ E(A, B)
        E(A, B) <--> E(P, Q)
        E(P, Q) <--> E(Q) + P
        E(Q) <--> E + Q
    end
end)
truth_key = ER._canonical(ER._law(truth), ER._mechanism_id(sprint(show, truth))).key
types = Dict(:mechanism_id => String, :parent_id => String, :fit_source_id => String)
load(f) = CSV.read(joinpath(DIR, f), DataFrame; types)
best, cv = load("best_equation.csv"), load("loocv_results.csv")
searched = [f for f in readdir(DIR) if f == "initial_mechanisms.csv" ||
                                       startswith(f, "equation_search_iteration_")]
fits = sum(nrow(filter(r -> r.fit_inherited === false, load(f))) for f in searched)
# Task 0's selection, read back from its text and keyed by this tree (the head's CSV has
# no fit_source_id column, so it is read without `types`).
base_csv = joinpath(dirname(BASE), "best_equation.csv")
base_text = only(CSV.read(base_csv, DataFrame).mechanism)
base_mech = ER._parse_mechanism(String(base_text))
base_key = ER._canonical(ER._law(base_mech), ER._mechanism_id(sprint(show, base_mech))).key
println("selected key equals the true key: ", only(best.key) == truth_key)
println("Task 0 selected the truth: ", base_key == truth_key)
println("selected key equals Task 0's selection: ", only(best.key) == base_key)
println("selected n_params: ", only(best.n_params), " (true: ",
        length(ER._solve(truth).fitted), ")")
println("loocv rows: ", nrow(cv), "; one row per key per count: ",
        allunique(zip(cv.n_params, cv.key)))
println("fits run: ", fits)
println("Task 0: ", join(filter(l -> startswith(l, "fits run") ||
                                     startswith(l, "selected n_params"),
                                readlines(BASE)), "; "))
```

Run (background, about 2 min): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T9_identify.log 5400 --project=. --heap-size-hint=2500M $S/phase45/scripts/P5_identify_run.jl $S/phase45/p5_T9_identify`,
polled in bounded loops; then
`S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T9_check.log 1800 --project=. $S/phase45/scripts/P5_T9_check.jl $S/phase45/p5_T9_identify $S/phase45/p5_T0_identify/summary.txt`
Expected, in `P5_T9_check.log`: `selected key equals the true key: true`, or, when the
line `Task 0 selected the truth: false` precedes it, either `true` or `selected key equals
Task 0's selection: true`, or a `selected n_params` at most `<n0>` with the selected row's
`cv_score` (`best_equation.csv`) below 1e-6; `one row per key per count: true`, `fits
run: <F9>` and Task 0's line `fits run (rows with fit_inherited == false): <F0>; selected
n_params: <n0>`, then `exit 0`. <F9> can fall or rise against
<F0>: one fit per key drops same-law refits, while the beam per key expands every member
of a selected key and lets distinct keys in that duplicates crowded out; the commit gives
both. Any other selection is investigated, not re-pinned: compare the two laws'
`cv_score` rows in `loocv_results.csv` and the per-round CSVs with Task 0's, and
report to the controller before the commit. Whether each run selected the truth goes to
the commit and the PR body's Measurements.

- [ ] **Step 6: Run the full suite**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T9_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`,
polled in bounded loops.
Expected: `exit 0`, every test passing, pristine, Aqua and JET green; about 1 s over Task
8's time; report `maxrss_kb` against Task 0's.

- [ ] **Step 7: Commit**

Write the measured figures in place of `<…>`.

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "$(cat <<'EOF'
Rank and charge the beam per key

The mechanisms of one key share one fit, so the beam ranks keys at their
lowest loss and charges the cumulative budget once per key: duplicates can
no longer crowd out a distinct key. Each key keeps its verdict for the
search, so a late member of a selected key expands without budget and a
late member of a rejected key stays out, and the cutoff and budget stay
monotone. The three _select_count! pins pass one key per mechanism and keep
every expected index.

Seeded identify run (P5_identify_run.jl): selected key equals the true
key: <yes or no> (Task 0 selected the truth: <yes or no>; equal to Task 0's
selection: <yes or no>); n_params <n9>, cv_score <cv9> (Task 0: <n0>);
<F9> fits (Task 0: <F0>); LOOCV holds one row per key and count.
Full-suite maxrss: <kb> KB (Task 0: <kb0> KB).

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
)"
```

---
### Task 10: L phantoms (P5-30 to P5-36, P5-31a)

A one-subunit variant whose inactive conformation binds nothing, and whose law has a
rescaling certificate, has the family of its law at L = 0 (RES §1.10). identify keeps it in
the frontier and fits none of its constants: it reads the fit of whoever holds that law's
key, with L = 0, and takes that source's fold losses in cross-validation. A phantom whose
law at L = 0 has no fit is fitted as an ordinary candidate (P5-32's b2). The pre-filter
stays (P5-31a).

**Files:**
- Modify: `src/key.jl` (ABOUTME; append `_phantom_weights`, `_is_phantom`)
- Modify: `src/identify_rate_equation.jl` (`_Candidate` gains `phantom` and `canon0`;
  PASS1 marks phantoms; representatives; rows with `fit_key = canon0.key` and
  `L_phantom`; `BatchEntry`'s and `_process_batch`'s docstrings; `_rows_to_dataframe`;
  `_failure_row`; `identify_rate_equation`'s "Model selection" docstring)
- Modify: `src/mechanism_enumeration.jl` (`_expand_to_allosteric`'s K-type paragraph)
- Modify: `docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md`, section 13
  "One-subunit L phantoms" (P5-30, P5-32) and section 10's L-phantom bullet (P5-32); ◆,
  listed for Denis in Task 15
- Modify: `ROADMAP.md` (retire the Planned L-phantom entry; rewrite the Idea "Count
  identifiable parameters with an exact modular rank"; add the Idea "Apply the chain-flank
  exclusion at one catalytic subunit" (P5-37, ◆); reword the PGK clause of "Profile the
  first v0.9.0 HPC runs")
- Modify: the fixtures file (`_testhelper_identifiable_rank` moved from TME;
  `_testhelper_phantom_map_holds` new)
- Test: `test/test_key.jl`, `test/test_mechanism_enumeration.jl` (five K-type fixture
  testsets; the rank helper's definition removed), `test/test_oracles.jl` (Q5, Q6),
  `test/test_identify_rate_equation.jl`
- Create (scratch): `$S/phase45/probes/P5_T10_verdicts.jl`,
  `$S/phase45/scripts/P5_phantom_rerun.jl`
- Logs: `$S/phase45/logs/P5_T10_{red,key,verdicts,tme,oracles,tid_red,tid,rerun,full}.log`

**Interfaces:**
- Consumes: Task 4's `_AlloLaw` (`vars`, `H`, `active`, `inactive`, `fA`, `fI`, `sites`,
  `n`), `_flip_only_law`; Task 7's `_canonical`; Task 8's `_Fit`, `_Candidate`,
  `_fit_order`, `_row_params`, `BatchEntry`, `_process_batch`, `_cv_model_selection`;
  Task 9's `_select_count!`; `_rref!` (`thermo.jl`), `_inactive_mask`, `_active_mask`,
  `catalytic_multiplicity`, `regulatory_sites`, `allo_states`. Tests: Task 5's
  `_testhelper_identifiable_rank` (as Task 5 left it); `_testhelper_law_inputs`,
  `_testhelper_steady_state_rate`, `_testhelper_is_binding`; test_oracles.jl's
  `_testhelper_content_text`; TID's `_testhelper_SpreadStubOpt`,
  `_testhelper_CountingStubOpt`.
- Produces:
  ```julia
  _phantom_weights(A::_AlloLaw)::Vector{Rational{Int}}  # w in A.vars order over the
                                                        #   fitted constants but L;
                                                        #   empty when no certificate
  _is_phantom(am::AllostericMechanism, A::_AlloLaw)::Bool
      # catalytic_multiplicity(am) == 1 && !any(_inactive_mask(am)) &&
      # every site ligand :OnlyA && !isempty(_phantom_weights(A))
  struct _Candidate                         # Task 8's, plus:
      …
      phantom::Bool                         # inherits the fit of its law at L = 0
      canon0::_Canonical                    # that law's canonical form; canon otherwise
  end
  # tests (fixtures file)
  _testhelper_identifiable_rank(m; npts = 60, ndraws = 3, h = 1e-5)::Int   # moved
  _testhelper_phantom_map_holds(am; rng)::Bool
  ```
  - CSV gains `L_phantom` after `fit_source_id`: `true` on a phantom's row, `false` on
    every other fitted row, `missing` on a failure row.
  - A phantom row: params `canon0.order => memo[canon0.key].values` plus `L => 0.0`, in
    its own `fitted_params` order; `loss`, `retcode` and `peak_at_range_edge` from the
    source; `n_params` its own; `fit_inherited = true`; `fit_source_id` the source's id.
  - `BatchEntry.fit_key = canon0.key`.
  - Test helpers (TID): `_testhelper_phantom_pair()`, `_testhelper_InfiniteStubOpt`;
    test_oracles.jl: `_testhelper_phantom_bibi_rxn`, `_testhelper_phantom_pingpong_rxn`,
    `_testhelper_spaced(ms, n)`, `_testhelper_tagging_children(parents, rxn)`; test_key.jl:
    the testset-local `_testhelper_exponent_toy_law(p, q)`.

**Gate:**
- `P5_T10_verdicts.jl` re-derives the K-type fixtures' verdicts under P5-31a and checks
  each flagged one against the rank oracle, `_testhelper_phantom_map_holds` and its
  catalytic mechanism's key before any verdict is pinned (Step 7 comes after it); its
  table matches Step 6's.
- Q5 passes (a9e0699 probe: 200 variants, 21 flagged, every map holds, four ranks
  equal); Q6 passes (172 children, none with an empty inactive conformation or a
  certificate).
- `P5_phantom_rerun.jl` prints `PASS`: every flagged ping-pong variant of round 1 and of
  a 200-parent round 2 has key(active law) = key(`_law(catalytic(am))`), np = np₀ + 1, a
  holding map and equal ranks. Its counts of flagged and certificate-only variants go to
  the commit for the PR body (P5-36), and so do the median and maximum time of
  `_is_phantom` per variant (spec section 17: "Phase 5 measures it").
- `test_key.jl`, TME, `test_oracles.jl` and TID are green; Q5 and Q6 add at most 10 s to
  `test_oracles.jl` (P5-41's 20 s covers Q1–Q6).
- The full suite is green and pristine; about +12 s; maxrss reported against Task 0's.

Facts this task rests on, measured on a9e0699 with Tasks 1–4's law code and this task's
certificate (`$S/phase45/plan5/probes/p3_probe.log`, `p3_probe2.log`, `p3c_probe.log`):
- The merged uni-uni K-type has w = [1, −1] over `K_A_EP_to_E_P`, `k_A_E_S_to_EP`.
- The Theorell–Chance `only_a` and `only_q` have certificates whose maps hold, yet keep a
  binding at free E in the inactive conformation: P5-31a's pre-filter fits them.
- With an `:OnlyA` site ligand the pre-filter passes, but no certificate exists: the
  regulator reveals L. This holds for every such site: each term of `X` (free of the
  regulator, since `W_I = 1`) has its multiple by `R/K_R` in `Y − X`, and the two rows ask
  `e·w = s + 1` and `e·w − w(K_R) = s` while the terms of `N` force `w(K_R) = 0`. So a
  phantom has no regulator constant, and its active law holds every fitted name but `L`.
- Every flagged fixture's active law, and #3268's, has its catalytic mechanism's key, so
  the merged uni-uni phantom of identify's tests reads its parent's fit.
- Census #3268 (the chain-flank flip) is admitted, flagged and its map holds.
- No ping-pong tagging child of 30 ping-pong parents at one subunit binds nothing in the
  inactive conformation or has a certificate, so ping-pong phantom coverage stays zero
  outside the rerun (P5-36).

- [ ] **Step 1: Move the rank oracle and add the map oracle**

Cut the docstring and function `_testhelper_identifiable_rank` from
`test/test_mechanism_enumeration.jl` (the block Task 5 left there, from its `"""` line
through its `end`) and paste it, unchanged, into
`test/mechanism_definitions_for_test_enzyme_derivation.jl` directly after the function
`_testhelper_declared_dimension`. test_oracles.jl runs before TME in the suite and reads it
from Q5.

After it, add:

```julia
"""
    _testhelper_phantom_map_holds(am; rng) -> Bool

Whether the one-subunit variant `am` reproduces the family of its law at `L = 0` through the
rescaling its certificate names (`EnzymeRates._phantom_weights`): at 3 consistent points
(`_testhelper_law_inputs`) for each `L` in {1/10, 10}, the exact steady-state rate of `am`
equals that of `catalytic(am)` with every fitted constant θ scaled by (1 + L)^w, to 1e-20
relative in 256-bit BigFloat (`_testhelper_steady_state_rate` on both sides). A dependent
constant scales by the power product the solve gives it, and each step's equilibrium
constant and forward rate are rebuilt from its scaled constants. False when `am` has no
certificate. Errors unless the active conformation keeps every edge, since only then is
the law at `L = 0` the law of `catalytic(am)`.
"""
function _testhelper_phantom_map_holds(am; rng)
    setprecision(BigFloat, 256) do
        s = EnzymeRates._solve(am)
        A = EnzymeRates._law(am, s)
        w = EnzymeRates._phantom_weights(A)
        isempty(w) && return false
        all(EnzymeRates._active_mask(am)) || error("_testhelper_phantom_map_holds: " *
            "the active conformation drops an edge, so the law at L = 0 is not the law " *
            "of catalytic(am)")
        θ = [A.vars.names[v] for v in eachindex(A.vars.names)
             if !A.vars.conc[v] && !(A.vars.names[v] in (:Keq, :L))]
        μ = Dict{Symbol, Rational{Int}}(zip(θ, w))
        for d in s.dependents
            μ[s.names[d.col]] = sum((e * μ[s.names[f]] for (f, e) in d.factors);
                                    init = 0 // 1)
        end
        holds = true
        for L in (big(1) / 10, big(10)), _ in 1:3
            point, concs, _ = _testhelper_law_inputs(am, s.fitted, rng)
            v = _testhelper_steady_state_rate(am, merge(point, (; L)), concs)
            scaled(x) = point.vals[x] * (1 + L)^BigFloat(μ[x])
            active = map(point.active) do (st, _, kf)
                kf === nothing || return (st, scaled(st.names[1]) / scaled(st.names[2]),
                                          scaled(st.names[1]))
                K = scaled(only(st.names))
                (st, _testhelper_is_binding(st) ? 1 / K : K, nothing)
            end
            v0 = _testhelper_steady_state_rate(EnzymeRates.catalytic(am),
                                               (; active, inactive = []), concs)
            holds &= abs(v - v0) <= big(10)^-20 * abs(v0)
        end
        holds
    end
end
```

- [ ] **Step 2: Write the failing certificate tests**

Append to `test/test_key.jl`:

```julia
# ── One-subunit L phantoms ─────────────────────────────────────────────────

@testset "a certificate maps the phantom onto its L→0 law" begin
    # The merged uni-uni K-type at one catalytic subunit: E + S <--> E(P) and
    # E + P ⇌ E(P), both :OnlyA. Its inactive conformation is free E alone, so the law is
    # v = N/(Y + L·X) with X a term of Y. The certificate scales P's dissociation constant
    # by (1 + L) and the uptake rate constant by 1/(1 + L): at every L the variant is its
    # law at L = 0, the law of its catalytic mechanism, at those constants (spec
    # section 13).
    am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(P)    :: OnlyA
            E + P ⇌ E(P)       :: OnlyA
        end
    end)
    A = EnzymeRates._law(am)
    θ = [A.vars.names[v] for v in eachindex(A.vars.names)
         if !A.vars.conc[v] && !(A.vars.names[v] in (:Keq, :L))]
    @test θ == [:K_A_EP_to_E_P, :k_A_E_S_to_EP]
    @test EnzymeRates._phantom_weights(A) == [1 // 1, -1 // 1]
    @test EnzymeRates._is_phantom(am, A)
    @test _testhelper_phantom_map_holds(am; rng = Random.MersenneTwister(21))
end

@testset "no phantom without the pre-filter's conditions" begin
    # P5-31a: the spec's pre-filter stays, and the certificate runs inside it. The merged
    # uni-uni K-type at two subunits, with P's binding :NonequalAI (the inactive
    # conformation binds P), and with an :OnlyI site ligand fails the pre-filter and has
    # no certificate. With an :OnlyA site ligand it passes the pre-filter, but the
    # regulator reveals L and no certificate exists. The Theorell–Chance K-type only_a
    # has a certificate, yet its inactive conformation binds Q at free E, so it is fitted.
    two_subunits = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        steps: begin
            E + S <--> E(P)    :: OnlyA
            E + P ⇌ E(P)       :: OnlyA
        end
    end)
    p_relaxed = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 1
        steps: begin
            E + S <--> E(P)    :: OnlyA
            E + P ⇌ E(P)       :: NonequalAI
        end
    end)
    onlyi_site = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R
        steps: begin
            E + S <--> E(P)    :: OnlyA
            E + P ⇌ E(P)       :: OnlyA
        end
        regulatory_site(1): R::OnlyI
    end)
    onlya_site = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R
        steps: begin
            E + S <--> E(P)    :: OnlyA
            E + P ⇌ E(P)       :: OnlyA
        end
        regulatory_site(1): R::OnlyA
    end)
    only_a = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        steps: begin
            E + A ⇌ E(A)                :: OnlyA
            E(A) + B <--> E(Q) + P      :: OnlyA
            E + Q ⇌ E(Q)                :: EqualAI
        end
    end)
    for am in (two_subunits, p_relaxed, onlyi_site, onlya_site, only_a)
        @test !EnzymeRates._is_phantom(am, EnzymeRates._law(am))
    end
    for am in (two_subunits, p_relaxed, onlyi_site, onlya_site)
        @test isempty(EnzymeRates._phantom_weights(EnzymeRates._law(am)))
    end
    @test !any(EnzymeRates._inactive_mask(onlya_site))
    @test !isempty(EnzymeRates._phantom_weights(EnzymeRates._law(only_a)))
    @test any(EnzymeRates._inactive_mask(only_a))
end

@testset "a certificate's rows overflow to no certificate" begin
    # A one-subunit toy law v = A·k1^p / (1 + k1^p·A + k1·k2^q·A + L), its inactive
    # conformation free E alone. Its certificate weighs k1 by −1/p and k2 by
    # (1 − p)/(p·q): −1/2 and −1/6 for p = 2, q = 3. For p = 2^40 + 15 and q = 2^40 + 3 the
    # second weight's denominator leaves Int64, so the row reduction overflows and the
    # variant gets no certificate: identify fits it.
    function _testhelper_exponent_toy_law(p, q)
        vars = EnzymeRates._Vars()
        a = EnzymeRates._var!(vars, :A; conc = true)
        k1, k2 = EnzymeRates._var!(vars, :k_one), EnzymeRates._var!(vars, :k_two)
        EnzymeRates._var!(vars, :Keq)
        EnzymeRates._var!(vars, :L)
        t(ps...) = EnzymeRates._term(EnzymeRates._monomial(ps...))
        D = one(EnzymeRates._Poly) + t(k1 => p, a => 1) + t(k1 => 1, k2 => q, a => 1)
        LA = EnzymeRates._Law(vars, t(a => 1), t(k1 => p), D, [D],
                              [one(EnzymeRates._Poly)], EnzymeRates._Mono())
        EnzymeRates._AlloLaw(vars, t(a => 1), LA, EnzymeRates._flip_only_law(LA),
                             one(EnzymeRates._Poly), one(EnzymeRates._Poly),
                             Tuple{EnzymeRates._Poly, EnzymeRates._Poly, Int}[], 1)
    end
    @test EnzymeRates._phantom_weights(_testhelper_exponent_toy_law(2, 3)) ==
          [-1 // 2, -1 // 6]
    p, q = 2^40 + 15, 2^40 + 3
    @test denominator(big(1 - p) // (big(p) * q)) > typemax(Int)
    @test isempty(EnzymeRates._phantom_weights(_testhelper_exponent_toy_law(p, q)))
end
```

- [ ] **Step 3: Run the certificate tests to verify they fail**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T10_red.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); @testset "red" begin include("test/test_key.jl") end'`
Expected: Task 7's testsets pass; the three new ones error with
``UndefVarError: `_phantom_weights` not defined in `EnzymeRates` `` (or `_is_phantom`);
the outer `red` testset ends the log with `exit 1`.

- [ ] **Step 4: Write the certificate**

In `src/key.jl`, replace the two ABOUTME lines with

```julia
# ABOUTME: The canonical key of a rate law, by colour refinement and individualization,
# ABOUTME: and the rescaling certificate that flags a one-subunit L phantom.
```

and append:

```julia
# ─── One-subunit L phantoms ─────────────────────────────────────────────────

"""
    _phantom_weights(A) -> w

The rescaling certificate of the MWC law `A` (spec section 13, "One-subunit L phantoms").
At one catalytic subunit, with an inactive conformation that turns nothing over
(`QI = 0`), the law is `v = N/(Y + L·X)` with `N = H·fI·QA·WA`, `Y = fI·DA·WA` and
`X = fA·DI·WI`. When every term of `X` is a term of `Y` with the same coefficient,
`Y + L·X = (Y − X) + (1 + L)·X`; weights `w` over the fitted constants other than `L`, and
a shift `s`, with `e·w = s` on the constant exponents `e` of every term of `N` and of
`Y − X` and `e·w = s + 1` on those of every term of `X`, then give
`v(θ, L) = v₀(θ ⊙ (1 + L)^w)` for all positive `θ` and `L`, `v₀ = N/Y` being the law at
`L = 0`: the two laws have one family. Returns `w` in `A.vars` order (the independents,
then the regulator constants), or an empty vector when no weights exist, when `A` has
another shape, or when the arithmetic leaves `Int64`, so that such a variant is fitted.
"""
function _phantom_weights(A::_AlloLaw)
    none = Rational{Int}[]
    A.n == 1 && iszero(A.inactive.Q) || return none
    try
        WA = foldl(*, (w for (w, _, k) in A.sites for _ in 1:k); init = one(_Poly))
        WI = foldl(*, (w for (_, w, k) in A.sites for _ in 1:k); init = one(_Poly))
        N = A.H * A.fI * A.active.Q * WA
        Y = A.fI * A.active.D * WA
        X = A.fA * A.inactive.D * WI
        all(((m, c),) -> get(Y.terms, m, 0) == c, X.terms) || return none
        vars = A.vars
        θ = [v for v in eachindex(vars.names)
             if !vars.conc[v] && vars.names[v] !== :Keq && vars.names[v] !== :L]
        at = Dict(v => j for (j, v) in enumerate(θ))
        n = length(θ)
        # One row per distinct exponent vector: [e; -1] · [w; s] = rhs.
        rows = Dict{Vector{Rational{Int}}, Rational{Int}}()
        function row!(p::_Poly, rhs::Rational{Int})
            for m in keys(p.terms)
                e = zeros(Rational{Int}, n + 1)
                for (v, x) in m
                    j = get(at, v, 0)
                    j == 0 || (e[j] += x)
                end
                e[end] = -1
                get!(rows, e, rhs) == rhs || return false
            end
            true
        end
        row!(N, 0 // 1) && row!(Y - X, 0 // 1) && row!(X, 1 // 1) || return none
        M = zeros(Rational{Int}, length(rows), n + 2)
        for (i, (e, rhs)) in enumerate(rows)
            M[i, 1:(n + 1)] = e
            M[i, n + 2] = rhs
        end
        pivots = _rref!(M, n + 1)
        all(iszero, view(M, (length(pivots) + 1):size(M, 1), n + 2)) || return none
        w = zeros(Rational{Int}, n + 1)
        for (r, c) in enumerate(pivots)
            w[c] = M[r, n + 2]
        end
        w[1:n]
    catch err
        err isa OverflowError || rethrow()
        none
    end
end

"""
    _is_phantom(am, A) -> Bool

Whether the allosteric mechanism `am`, whose MWC law is `A`, is a one-subunit L phantom
(spec section 13): one catalytic subunit, an inactive conformation that binds nothing and
runs no step (`_inactive_mask` empty), only `:OnlyA` site ligands, and a rescaling
certificate (`_phantom_weights`), so its family is the family of its law at `L = 0`.
"""
_is_phantom(am::AllostericMechanism, A::_AlloLaw) =
    catalytic_multiplicity(am) == 1 && !any(_inactive_mask(am)) &&
    all(site -> all(==(:OnlyA), allo_states(site)), regulatory_sites(am)) &&
    !isempty(_phantom_weights(A))
```

The closure `row!` captures `rows`, `at` and `n`, each assigned once, so nothing is boxed.
The certificate's products stay far below `Int64` under the cap (P5-2: `τ_I·F·τ_A` about
4e8).

- [ ] **Step 5: Run the certificate tests to verify they pass**

Rerun Step 3's command with log `$S/phase45/logs/P5_T10_key.log`.
Expected: every testset passes, pristine, `exit 0`.

- [ ] **Step 6: Check the fixture verdicts before pinning them**

Create `$S/phase45/probes/P5_T10_verdicts.jl`:

```julia
# ABOUTME: Task 10's value check before the K-type verdict pins: each TME fixture's
# ABOUTME: pre-filter, certificate, map oracle, and L→0 key and rank against its parent's.
# Usage, from the worktree root: julia --project=. P5_T10_verdicts.jl
using TestEnv; TestEnv.activate()
using Test, EnzymeRates, LinearAlgebra, Random, Printf
include(joinpath(pkgdir(EnzymeRates), "test",
                 "mechanism_definitions_for_test_enzyme_derivation.jl"))
const ER = EnzymeRates
prefilter(am) = ER.catalytic_multiplicity(am) == 1 && !any(ER._inactive_mask(am)) &&
                all(s -> all(==(:OnlyA), ER.allo_states(s)), ER.regulatory_sites(am))
function verdict(label, am)
    A = ER._law(am)
    cert = !isempty(ER._phantom_weights(A))
    flagged = ER._is_phantom(am, A)
    holds = cert ? _testhelper_phantom_map_holds(am; rng = MersenneTwister(7)) : missing
    # A flagged variant reads the fit of its active law's key; the merged uni-uni pair of
    # identify's tests needs that key to be its catalytic mechanism's.
    key0 = flagged ? ER._canonical(A.active, "variant").key ==
                     ER._canonical(ER._law(ER.catalytic(am)), "parent").key : missing
    ranks = flagged ? (_testhelper_identifiable_rank(am),
                       _testhelper_identifiable_rank(ER.catalytic(am))) : missing
    @printf("%-44s pre-filter %-5s certificate %-5s phantom %-5s map %-7s key %-7s ",
            label, prefilter(am), cert, flagged, holds, key0)
    println("ranks ", ranks)
end
bibi = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end
bibi2 = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    oligomeric_state: 2
end
uni = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
end
parents = [
    ("merged uni-uni", uni, ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(P)
            E + P ⇌ E(P)
        end
    end)),
    ("Theorell–Chance", bibi, ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(Q) + P
            E + Q ⇌ E(Q)
        end
    end)),
    ("mixed B group", bibi, ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))
            E(Q) + P <--> E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)),
    ("random-order merged", bibi, ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> E(A), E(B) + A <--> E(P, Q))
            (E + B ⇌ E(B), E(A) + B ⇌ E(P, Q))
            E + Q ⇌ E(Q)
            E(Q) + P <--> E(P, Q)
        end
    end)),
    ("merged ordered, n = 2", bibi2, ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B ⇌ E(P, Q)
            E(Q) + P <--> E(P, Q)
            E + Q ⇌ E(Q)
        end
    end))]
for (label, rxn, m) in parents
    println("== ", label, ", parent rank ", _testhelper_identifiable_rank(m))
    for k in ER._expand_to_allosteric(m, rxn)
        verdict(join(string.(k.cat_allo_states), ","), k)
    end
end
p_relaxed = ER.AllostericMechanism(@allosteric_mechanism begin
    substrates: S
    products: P
    catalytic_multiplicity: 1
    steps: begin
        E + S <--> E(P)    :: OnlyA
        E + P ⇌ E(P)       :: NonequalAI
    end
end)
verdict("p_relaxed", p_relaxed)
```

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T10_verdicts.log 1800 --project=. $S/phase45/probes/P5_T10_verdicts.jl`
Expected (the a9e0699 probes' verdicts, `p3_probe.log` section A and `p3c_probe.log`
section B; the tags are each child's `cat_allo_states`, in its canonical group order):
- merged uni-uni (parent rank 2): its one child `OnlyA,OnlyA` flagged, map true, key true,
  ranks (2, 2).
- Theorell–Chance (parent rank 3): `OnlyA,OnlyA,OnlyA` (a_and_q) flagged, map true, key
  true, ranks (3, 3); `OnlyA,EqualAI,OnlyA` (only_a) and `EqualAI,OnlyA,OnlyA` (only_q)
  pre-filter false, certificate true, map true, not flagged.
- mixed B group (parent rank 5): the two children with one `EqualAI` group at P (a_and_q)
  and with all four `OnlyA` (all_three) flagged, map true, key true, ranks (5, 5); the
  other five pre-filter false, certificate false.
- random-order merged (parent rank 5): two children pass the pre-filter (only_q, P's
  group `EqualAI`, and p_and_q) with certificate false; only_p pre-filter false; none
  flagged.
- merged ordered, n = 2: seven children, all pre-filter false, certificate false.
- `p_relaxed`: pre-filter false, certificate false.
A verdict that differs is not pinned: compare it with the rank pair and the map oracle,
and stop and report to the controller (spec section 14; RES §3.4).

- [ ] **Step 7: Pin the verdicts and add Q5 and Q6**

Step 6's table is the value check these pins rest on (spec section 14): write them only
once every verdict there matches. In `test/test_mechanism_enumeration.jl`:

In testset "expand_mechanisms on a merged uni-uni", replace

```julia
    @test _testhelper_fitted(k_type) == 3 && _testhelper_identifiable_rank(k_type) == 2
end
```

with

```julia
    @test _testhelper_fitted(k_type) == 3 && _testhelper_identifiable_rank(k_type) == 2
    # The K-type child is an L phantom (spec section 13): one catalytic subunit, an
    # inactive conformation that is free E alone, and a rescaling certificate that maps
    # it onto the law of its catalytic mechanism.
    @test EnzymeRates._is_phantom(k_type, EnzymeRates._law(k_type))
    @test _testhelper_phantom_map_holds(k_type; rng = Random.MersenneTwister(1))
end
```

In testset "_expand_to_allosteric on an ordered Theorell–Chance mechanism", replace

```julia
    for k in (only_a, a_and_q)
        @test _testhelper_identifiable_rank(k) == 3
    end
end
```

with

```julia
    for k in (only_a, a_and_q)
        @test _testhelper_identifiable_rank(k) == 3
    end
    # Only a_and_q is an L phantom. only_a and only_q keep a binding at free E in the
    # inactive conformation, so the pre-filter fits them, though each has a certificate
    # (P5-31a).
    @test EnzymeRates._is_phantom(a_and_q, EnzymeRates._law(a_and_q))
    @test _testhelper_phantom_map_holds(a_and_q; rng = Random.MersenneTwister(2))
    for k in (only_a, only_q)
        @test !EnzymeRates._is_phantom(k, EnzymeRates._law(k))
        @test !isempty(EnzymeRates._phantom_weights(EnzymeRates._law(k)))
    end
end
```

In testset "_expand_to_allosteric on a merged bi-bi whose B group mixes fused and plain",
replace

```julia
    for k in (a_and_q, all_three)
        @test _testhelper_fitted(k) == 6 && _testhelper_identifiable_rank(k) == 5
    end
end
```

with

```julia
    for k in (a_and_q, all_three)
        @test _testhelper_fitted(k) == 6 && _testhelper_identifiable_rank(k) == 5
    end
    # The two whose inactive conformation binds nothing are L phantoms, each mapped onto
    # the parent's law by its certificate; the five that keep a binding there are fitted.
    for k in (a_and_q, all_three)
        @test EnzymeRates._is_phantom(k, EnzymeRates._law(k))
        @test _testhelper_phantom_map_holds(k; rng = Random.MersenneTwister(3))
    end
    for k in (only_a, only_p, only_q, a_and_p, p_and_q)
        @test !EnzymeRates._is_phantom(k, EnzymeRates._law(k))
    end
end
```

In testset "_expand_to_allosteric on a random-order merged bi-bi led by plain bindings",
replace

```julia
    for k in (only_q, only_p, p_and_q)
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 6
    end
end
```

with

```julia
    for k in (only_q, only_p, p_and_q)
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 6
    end
    # None is an L phantom. only_q and p_and_q bind nothing in the inactive conformation
    # and pass the pre-filter, but L raises their rank and no certificate exists; only_p
    # keeps Q's binding there.
    for k in (only_q, p_and_q)
        @test !any(EnzymeRates._inactive_mask(k))
        @test isempty(EnzymeRates._phantom_weights(EnzymeRates._law(k)))
    end
    @test !any(k -> EnzymeRates._is_phantom(k, EnzymeRates._law(k)),
               (only_q, only_p, p_and_q))
end
```

In testset "_expand_to_allosteric on a merged ordered bi-bi at multiplicity 2", replace

```julia
    @test _testhelper_fitted(m2) == _testhelper_identifiable_rank(m2) == 5
    for k in kids
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 6
    end
end
```

with

```julia
    @test _testhelper_fitted(m2) == _testhelper_identifiable_rank(m2) == 5
    for k in kids
        @test _testhelper_fitted(k) == _testhelper_identifiable_rank(k) == 6
    end
    # At two subunits no variant is an L phantom.
    @test !any(k -> EnzymeRates._is_phantom(k, EnzymeRates._law(k)), kids)
end
```

Append to `test/test_oracles.jl`, at top level after its last testset. (Task 7's Q3 sits
inside the `let` block of the inhibitor children, where these `const`s would be refused
and these functions would be local, so the new code does not go next to it.)

```julia
# ── L phantoms over regenerated populations ─────────────────────────────────

# The census bi-bi reaction at one catalytic subunit, whose seeds' tagging children are
# Q5's population.
const _testhelper_phantom_bibi_rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end

# The census ping-pong declaration at one catalytic subunit, whose ping-pong parents'
# tagging children are Q6's population.
const _testhelper_phantom_pingpong_rxn = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    allosteric_regulators: R
    shared_catalytic_site: (A, P)
end

"""`n` of the mechanisms `ms`, spaced evenly through their content order
(`_testhelper_content_text`), or all of them when there are at most `n`."""
function _testhelper_spaced(ms, n)
    sorted = ms[sortperm(_testhelper_content_text.(ms))]
    length(sorted) <= n ? sorted :
        sorted[round.(Int, range(1, length(sorted); length = n))]
end

"""The well-posed tagging children (`_expand_to_allosteric`) of `parents` over `rxn`,
without duplicates, in content order (`_testhelper_content_text`)."""
function _testhelper_tagging_children(parents, rxn)
    kids = unique!([c for p in parents for c in EnzymeRates._expand_to_allosteric(p, rxn)
                    if EnzymeRates._is_well_posed(c)])
    kids[sortperm(_testhelper_content_text.(kids))]
end

@testset "every flagged one-subunit bi-bi variant is its L→0 law's family" begin
    # Q5 (spec section 13; P5-30, P5-31a): the tagging children of 40 bi-bi seeds spaced
    # through their content order, 200 of them spaced the same way. Each variant
    # `_is_phantom` flags reproduces its law at L = 0 under the certificate's rescaling, so
    # its family is that law's and no flagged variant gains rank; the rank oracle confirms
    # it on the first four. The counts are aggregate regression pins over the
    # declaration's seed set (CLAUDE.md enumeration rule 1's exception).
    rxn = _testhelper_phantom_bibi_rxn
    parents = _testhelper_spaced(collect(EnzymeRates.init_mechanisms(rxn)), 40)
    kids = _testhelper_spaced(_testhelper_tagging_children(parents, rxn), 200)
    @test length(kids) == 200
    flagged = [am for am in kids if EnzymeRates._is_phantom(am, EnzymeRates._law(am))]
    @test length(flagged) == 21
    rng = Random.MersenneTwister(505)
    @test all(am -> _testhelper_phantom_map_holds(am; rng), flagged)
    for am in first(flagged, 4)
        @test _testhelper_identifiable_rank(am) ==
              _testhelper_identifiable_rank(EnzymeRates.catalytic(am))
    end
end

@testset "no one-subunit ping-pong tagging child is an L phantom" begin
    # Q6 (P5-36): the tagging children at one catalytic subunit of 30 ping-pong parents,
    # the seeds with a covalent form and their well-posed flips with one, spaced through
    # their content order. Every child keeps a step in the inactive conformation, so the
    # pre-filter fits it, and none has a certificate, so dropping the pre-filter (P5-31)
    # would flag none either. The counts are aggregate regression pins over the
    # declaration's seed set (CLAUDE.md enumeration rule 1's exception).
    rxn = _testhelper_phantom_pingpong_rxn
    covalent(m) = any(EnzymeRates._has_residual, m.forms)
    seeds = filter(covalent, collect(EnzymeRates.init_mechanisms(rxn)))
    flips = [c for m in seeds for c in EnzymeRates._expand_re_to_ss(m)
             if EnzymeRates._is_well_posed(c) && covalent(c)]
    parents = _testhelper_spaced(unique!([seeds; flips]), 30)
    @test length(parents) == 30
    kids = _testhelper_tagging_children(parents, rxn)
    @test length(kids) == 172
    @test all(am -> any(EnzymeRates._inactive_mask(am)), kids)
    @test all(am -> isempty(EnzymeRates._phantom_weights(EnzymeRates._law(am))), kids)
end
```

- [ ] **Step 8: Run the TME and population tests**

Run (background) the TME focused run, then test_oracles.jl:
`S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && for f in test_mechanism_enumeration:tme test_oracles:oracles; do $S/phase3/jl.sh $S/phase45/logs/P5_T10_${f#*:}.log 3000 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/'${f%%:*}'.jl")'; done`,
polled in turn: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; L=$S/phase45/logs/P5_T10_tme.log; timeout 580 bash -c "until grep -q '^exit ' $L; do sleep 20; done"; tail -30 $L`.
Expected: both pass, pristine; Q5 flags 21 of its 200 variants and Q6 reads 172 children,
the a9e0699 probe's counts (f9361cf builds the same 200 variants and 172 children). A
count that differs is checked before it is pinned: every flagged variant passes the map
oracle and the first four keep their catalytic mechanism's rank; report it to the
controller. In `P5_T10_oracles.log` read the `Test Summary` times of
the Q5 and Q6 testsets: together at most 10 s (a9e0699 probe: Q5 6.5 s including its
first-call compile, Q6 under 2 s once warm; f9361cf, the same code in a cold fresh
process with the head's compile-based rank oracle, which Task 5 replaces: Q5 12.8 s, Q6
0.3 s, warm 0.2 s each). The cap reads the suite's warm `Test Summary` time. Over 10 s, stop and report: P5-41's 20 s for
Q1–Q6 is Denis's budget.

- [ ] **Step 9: Write the failing identify tests**

In `test/test_identify_rate_equation.jl`:

(a) Testset "_rows_to_dataframe": in the first rows literal add `L_phantom = false,`
after `fit_source_id = "0123456789abcdef",`; in the expected names replace
`"fit_source_id", "peak_at_range_edge", "a", "b"]` with
`"fit_source_id", "L_phantom", "peak_at_range_edge", "a", "b"]`, and after the line
`@test df.fit_source_id == ["0123456789abcdef"]` add `@test df.L_phantom == [false]`. In
`rows_with_failure` add `L_phantom = false,` to the fitted row and `L_phantom = missing,`
to the failure row, and after `@test ismissing(df_with_failure.fit_source_id[2])` add
`@test ismissing(df_with_failure.L_phantom[2])`. In the last row literal add
`L_phantom = false,` after `fit_source_id = "0",`.

(b) Testset "a failure row holds its mechanism's text and id": replace

```julia
        @test ismissing(row.key) && ismissing(row.fit_inherited) &&
              ismissing(row.fit_source_id)
```

with

```julia
        @test ismissing(row.key) && ismissing(row.fit_inherited) &&
              ismissing(row.fit_source_id) && ismissing(row.L_phantom)
```

(c) Testset "CSV output: every file's columns, texts and ids": replace

```julia
            @test "key" in names(df_file) && "fit_source_id" in names(df_file)
```

with

```julia
            @test "key" in names(df_file) && "fit_source_id" in names(df_file) &&
                  "L_phantom" in names(df_file)
```

and, in the loop that follows, replace

```julia
            for r in eachrow(df_file)
                ismissing(r.key) && continue
```

with

```julia
            for r in eachrow(df_file)
                # A failure row has no flag; a flagged row is an L phantom that reads
                # another mechanism's fit with L = 0.
                if ismissing(r.key)
                    @test ismissing(r.L_phantom)
                    continue
                end
                if r.L_phantom
                    m = EnzymeRates._parse_mechanism(r.mechanism)
                    @test EnzymeRates._is_phantom(m, EnzymeRates._law(m))
                    @test r.L == 0.0 && r.fit_inherited
                end
```

(d) Testset "csv writer": add `L_phantom = false,` to the rows literal after
`fit_source_id = "00000000000000aa",`.

(e) Testset "_cv_model_selection: flatten reproduces serial LOOCV per candidate": in
`_testhelper_pool_entry`'s row replace `fit_inherited=false, fit_source_id=id)` with
`fit_inherited=false, fit_source_id=id, L_phantom=false)`.

(f) Testset "best_equation.csv reads back with its key and ids": after
`@test plain.key[1] == e.key && occursin(r"^[0-9a-f]{64}$", plain.key[1])` add
`@test plain.L_phantom[1] === false`.

(g) After the testset "a law at the leaf cap is fitted under its mechanism's key", add:

```julia
# A stub optimizer whose every restart reports an infinite objective, as a fit with no
# finite loss does; it counts `solve` calls.
mutable struct _testhelper_InfiniteStubOpt
    count::Int
end
Optimization.allowsbounds(::_testhelper_InfiniteStubOpt) = true
function Optimization.SciMLBase.__solve(
        prob::Optimization.OptimizationProblem, opt::_testhelper_InfiniteStubOpt; kwargs...)
    opt.count += 1
    u = zeros(length(prob.u0))
    build_solution(DefaultOptimizationCache(prob.f, prob.p), opt, u, Inf;
                   retcode = ReturnCode.Failure)
end

# The merged uni-uni and its K-type child at one catalytic subunit
# (test_mechanism_enumeration.jl, "expand_mechanisms on a merged uni-uni"): the child's
# inactive conformation is free E alone and its certificate maps it onto the parent's
# law, so it is an L phantom whose law at L = 0 has the parent's key.
function _testhelper_phantom_pair()
    parent = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(P)
            E + P ⇌ E(P)
        end
    end)
    phantom = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(P)    :: OnlyA
            E + P ⇌ E(P)       :: OnlyA
        end
    end)
    parent, phantom
end

@testset "a phantom inherits the fit of whoever holds its L→0 key" begin
    # Spec section 13, P5-32: the phantom is not fitted. It reads the fit of its law at
    # L = 0 from the parent, which holds that key and is fitted in the same batch, with
    # L = 0, the parent's loss and its own parameter count; its row names the parent as
    # the fit's source, and the two rates agree.
    parent, phantom = _testhelper_phantom_pair()
    prob = _testhelper_uni_prob(NamedTuple)
    @test EnzymeRates._is_phantom(phantom, EnzymeRates._law(phantom))
    opt = _testhelper_SpreadStubOpt()
    batch_log = String[]
    entries, failures = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[phantom, parent],
        prob; optimizer = opt, max_param_count = 20, eq_complexity_filter = 337,
        n_restarts = 1, maxtime = 1.0, log = msg -> push!(batch_log, msg))
    @test opt.count == 1 && isempty(failures)
    @test batch_log == ["1 new fits + 1 inherited + 0 skipped (already fit) + " *
                        "0 skipped (>20 params) + 0 skipped (>337 complexity)"]
    p = only(filter(e -> e.mech == phantom, entries))
    f = only(filter(e -> e.mech == parent, entries))
    @test p.row.L_phantom && !f.row.L_phantom
    @test p.row.params.L == 0.0
    @test p.loss == f.loss && p.retcode === f.retcode
    @test (p.n_params, f.n_params) == (3, 2)
    @test p.fit_key == f.key != p.key
    @test p.row.fit_inherited && !f.row.fit_inherited
    @test p.row.fit_source_id == f.row.mechanism_id == f.row.fit_source_id
    @test p.row.peak_at_range_edge == f.row.peak_at_range_edge
    rng = Random.MersenneTwister(13)
    for _ in 1:3
        concs = (S = exp(randn(rng)), P = exp(randn(rng)))
        _testhelper_rate(e) =
            rate_equation(EnzymeRates.compile_mechanism(e.mech), concs,
                          merge(e.row.params, (Keq = 10.0, E_total = 1.0)))
        @test isapprox(_testhelper_rate(p), _testhelper_rate(f); rtol = 1e-12)
    end
end

@testset "a phantom whose L→0 key has no fit is fitted unflagged" begin
    # P5-32's b2: census #3268, the RE→SS flip of B's binding in the K-type of an ordered
    # bi-bi, is an L phantom whose law at L = 0 no mechanism of the batch or the memo
    # holds (the plain flip move never emits it: `_chain_flank_groups` keeps B's group
    # at rapid equilibrium). It is fitted as an ordinary candidate under its own key,
    # unflagged, its L fitted.
    m3268 = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)              :: OnlyA
            E + Q ⇌ E(Q)              :: OnlyA
            E(A) + B <--> E(A, B)     :: EqualAI
            E(A, B) <--> E(P, Q)      :: OnlyA
            E(Q) + P ⇌ E(P, Q)        :: EqualAI
        end
    end)
    @test EnzymeRates._is_phantom(m3268, EnzymeRates._law(m3268))
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            A = [1.0, 2.0, 1.0, 2.0], B = [0.5, 0.5, 1.0, 1.0],
            P = [0.1, 0.2, 0.1, 0.2], Q = [0.3, 0.3, 0.4, 0.4])
    prob = IdentifyRateEquationProblem(rxn, data; Keq = 2.0)
    opt = _testhelper_SpreadStubOpt()
    entries, failures = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m3268], prob;
        optimizer = opt, max_param_count = 20, n_restarts = 1, maxtime = 1.0)
    e = only(entries)
    @test opt.count == 1 && isempty(failures)
    @test !e.row.L_phantom && !e.row.fit_inherited
    @test e.fit_key == e.key
    @test e.n_params == 7 && e.row.params.L != 0.0
end

@testset "a phantom reads its source's fold losses" begin
    # P5-34: equal families give equal fold optima, so LOOCV fits the source once per
    # fold and the phantom takes those fold losses. Beside its source the phantom, one
    # parameter larger, loses the 1-SE tie.
    parent, phantom = _testhelper_phantom_pair()
    prob = _testhelper_uni_prob(NamedTuple)
    memo = Dict{String, EnzymeRates._Fit}()
    entries, _ = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[parent, phantom],
        prob; optimizer = _testhelper_SpreadStubOpt(), max_param_count = 20, memo,
        n_restarts = 1, maxtime = 1.0)
    opt = _testhelper_CountingStubOpt(; uval = 0.0)
    res = EnzymeRates._cv_model_selection(entries, memo, prob; optimizer = opt,
        se_threshold = 1.0, save_dir = mktempdir(), show_progress = false,
        n_restarts = 1, maxtime = 1.0)
    groups = unique(prob.data.group)
    @test opt.count == length(groups)
    folds(r) = [r[Symbol("cv_fold_$g")] for g in groups]
    @test res.cv_results.n_params == [2, 3]
    @test folds(res.cv_results[1, :]) == folds(res.cv_results[2, :])
    @test res.best === EnzymeRates.compile_mechanism(parent)
end

@testset "a phantom chosen by the 1-SE rule reports L = 0 and evaluates" begin
    # Review Focus 3: the pool holds the phantom without its source, so the 1-SE rule
    # selects it. best_equation.csv reports L = 0, the flag and the source's id, and its
    # parameters, read back from the CSV, give the source's rates through the phantom's
    # own rate equation; rescale_parameter_values accepts them.
    parent, phantom = _testhelper_phantom_pair()
    prob = _testhelper_uni_prob(NamedTuple)
    memo = Dict{String, EnzymeRates._Fit}()
    entries, _ = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[parent, phantom],
        prob; optimizer = _testhelper_SpreadStubOpt(), max_param_count = 20, memo,
        n_restarts = 1, maxtime = 1.0)
    p = only(filter(e -> e.mech == phantom, entries))
    f = only(filter(e -> e.mech == parent, entries))
    mktempdir() do dir
        res = EnzymeRates._cv_model_selection([p], memo, prob;
            optimizer = _testhelper_CountingStubOpt(; uval = 0.0), se_threshold = 1.0,
            save_dir = dir, show_progress = false, n_restarts = 1, maxtime = 1.0)
        @test res.best === EnzymeRates.compile_mechanism(phantom)
        best = CSV.read(joinpath(dir, "best_equation.csv"), DataFrame;
            types = Dict(:mechanism_id => String, :parent_id => String,
                         :fit_source_id => String))
        @test best.L[1] == 0.0 && best.L_phantom[1] === true
        @test best.fit_source_id[1] == f.row.mechanism_id
        fp = EnzymeRates.fitted_params(res.best)
        params = merge(NamedTuple{fp}(Tuple(Float64(best[1, x]) for x in fp)),
                       (Keq = prob.Keq, E_total = 1.0))
        source = merge(f.row.params, (Keq = prob.Keq, E_total = 1.0))
        rng = Random.MersenneTwister(14)
        for _ in 1:3
            concs = (S = exp(randn(rng)), P = exp(randn(rng)))
            v = rate_equation(res.best, concs, params)
            v0 = rate_equation(EnzymeRates.compile_mechanism(parent), concs, source)
            @test isfinite(v) && isapprox(v, v0; rtol = 1e-12)
        end
        scaled, spread = rescale_parameter_values(res.best, params)
        @test spread === nothing && scaled.L == 0.0
    end
end

@testset "a phantom whose source fit fails is a failure row" begin
    # Review Focus 4. When the source's fit throws, the phantom fails with the source's
    # error and the batch returns. When the source's fit reaches no finite loss, the
    # phantom's row carries the source's loss and retcode, and `_ingest!` keeps both out of
    # the frontier and the cv pool.
    parent, phantom = _testhelper_phantom_pair()
    prob = _testhelper_uni_prob(NamedTuple)
    batch = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[parent, phantom]
    entries, failures = EnzymeRates._process_batch(batch, prob;
        optimizer = _testhelper_CountingStubOpt(; throwit = true), max_param_count = 20,
        n_restarts = 1, maxtime = 1.0)
    @test isempty(entries)
    @test [f.mech for f in failures] == [parent, phantom]
    @test failures[2].error == failures[1].error
    @test occursin("stub solver forced failure", failures[2].error)
    @test ismissing(EnzymeRates._failure_row(failures[2]).L_phantom)
    opt = _testhelper_InfiniteStubOpt(0)
    entries, failures = EnzymeRates._process_batch(batch, prob; optimizer = opt,
        max_param_count = 20, n_restarts = 1, maxtime = 1.0)
    @test opt.count == 1 && isempty(failures)
    @test all(e -> e.loss == Inf && e.retcode === :NoFiniteLoss, entries)
    @test only(filter(e -> e.mech == phantom, entries)).row.L_phantom
    frontier = Dict{Int, Vector{EnzymeRates.BatchEntry}}()
    pool = Dict{Int, Vector{EnzymeRates.BatchEntry}}()
    best = Dict{Int, Float64}()
    @test isempty(EnzymeRates._ingest!(frontier, pool, best, entries; n_cv_candidates = 5))
    @test isempty(frontier) && isempty(pool) && isempty(best)
end
```

- [ ] **Step 10: Run the identify tests to verify they fail**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T10_tid_red.log 3600 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); @testset "red" begin include("test/test_identify_rate_equation.jl") end'`,
polled as before.
Expected: the edited and new testsets fail or error, the rest pass. The rows have no
`L_phantom` yet (``FieldError: type NamedTuple has no field `L_phantom` ``, and no
`L_phantom` column), and
the phantom is fitted as a representative of its own key, so "a phantom inherits …" fails
`opt.count == 1 && isempty(failures)` and its log line reads `2 new fits + 0 inherited`,
"a phantom reads its source's fold losses" counts twice the folds, and the selected
phantom's `best.L` is its own fitted value, not 0. The outer `red` testset ends the log
with `exit 1`.

- [ ] **Step 11: Let phantoms inherit**

In `src/identify_rate_equation.jl`:

Replace Task 8's `_Candidate` docstring and struct with:

```julia
"""
A mechanism PASS1 of `_process_batch` derived: its fitted constants in `fitted_params`
order, its text and the text's id (`_mechanism_id`), its printed law, its law's canonical
form (`canon`), whether it is an L phantom that reads another law's fit (`_is_phantom`),
and `canon0`, the canonical form of the law whose fit it reads: for a phantom its active
law, its law at `L = 0`; `canon` otherwise.
"""
struct _Candidate
    mech::Union{Mechanism, AllostericMechanism}
    fitted::Vector{Symbol}
    mechanism::String
    mechanism_id::String
    eq_text::String
    canon::_Canonical
    phantom::Bool
    canon0::_Canonical
end
```

Replace `BatchEntry`'s docstring with

```julia
"""One fitted mechanism: its parameter count, its fit's loss and retcode, its rate-law key
(`_canonical`), the key of the fit it reads (`fit_key`: its own, or for an L phantom the
key of its law at `L = 0`), and its CSV row."""
```

In `_rows_to_dataframe`, after the line
`fit_source_id = Union{Missing,String}[r.fit_source_id for r in rows],` add

```julia
        L_phantom = Union{Missing,Bool}[r.L_phantom for r in rows],
```

In `_failure_row`, after `fit_source_id = missing,` add `L_phantom = missing,`.

In `_process_batch`, replace the PASS1 line

```julia
            _Candidate(m, s.fitted, text, id, eq_text, _canonical(L, id))
```

with

```julia
            canon = _canonical(L, id)
            # An L phantom reads the fit of its law at L = 0, its active law, when that law
            # holds every fitted constant but L. The certificate rules out a phantom with
            # regulator constants; the check keeps a row from reading names its fit lacks.
            phantom = m isa AllostericMechanism && _is_phantom(m, L::_AlloLaw)
            canon0 = phantom ? _canonical((L::_AlloLaw).active, id) : canon
            if phantom && Set(canon0.order) != Set(s.fitted[1:(end - 1)])
                phantom, canon0 = false, canon
            end
            _Candidate(m, s.fitted, text, id, eq_text, canon, phantom, canon0)
```

replace

```julia
    # Pick one representative per key not already fit.
    rep_idx = Dict{String,Int}()
    for (i, c) in enumerate(compiled)
        c isa _Candidate && !haskey(memo, c.canon.key) && get!(rep_idx, c.canon.key, i)
    end
```

with

```julia
    # Pick one representative per key not already fit; an L phantom is none.
    rep_idx = Dict{String,Int}()
    for (i, c) in enumerate(compiled)
        c isa _Candidate && !c.phantom && !haskey(memo, c.canon.key) &&
            get!(rep_idx, c.canon.key, i)
    end
    # A phantom inherits the fit of whoever holds its L → 0 key, in `memo` or among this
    # batch's representatives; any other phantom becomes an ordinary candidate (P5-32).
    for (i, c) in enumerate(compiled)
        c isa _Candidate && c.phantom || continue
        (haskey(memo, c.canon0.key) || haskey(rep_idx, c.canon0.key)) && continue
        compiled[i] = _Candidate(c.mech, c.fitted, c.mechanism, c.mechanism_id,
                                 c.eq_text, c.canon, false, c.canon)
        haskey(memo, c.canon.key) || get!(rep_idx, c.canon.key, i)
    end
```

and replace the row loop's body from `key = c.canon.key` through the `push!(entries, …)`
call with:

```julia
        fit_key = c.canon0.key                    # its own key, or a phantom's L → 0 key
        if haskey(fit_error, fit_key)             # the fit's source threw
            push!(failures, FitFailure(c.mech, fit_error[fit_key]))
            continue
        end
        fit = memo[fit_key]
        params = c.phantom ? _row_params(c, [c.canon0.order; :L], [fit.values; 0.0]) :
                             _row_params(c, _fit_order(c), fit.values)
        parent = get(parent_of, c.mech, nothing)
        row = (
            n_params = length(c.fitted),
            parent_n_params = parent === nothing ? missing : parent.n_params,
            loss = fit.loss,
            mechanism = c.mechanism,
            mechanism_id = c.mechanism_id,
            parent_id = parent === nothing ? missing : parent.row.mechanism_id,
            rate_equation = c.eq_text,
            retcode = string(fit.retcode),
            error = missing,
            params,
            key = c.canon.key,
            fit_inherited = fit.source_id != c.mechanism_id,
            fit_source_id = fit.source_id,
            L_phantom = c.phantom,
            peak_at_range_edge = fit.peak_at_range_edge,
        )
        push!(entries, BatchEntry(c.mech, length(c.fitted), fit.loss, fit.retcode,
                                  c.canon.key, fit_key, row))
```

In `_process_batch`'s docstring, after the PASS-2 paragraph's sentence "A representative
whose fit throws fails ALL the mechanisms with its key (all-or-nothing per law).", add:

```
An L phantom (`_is_phantom`: one catalytic subunit, an inactive conformation that binds
nothing, and a rescaling certificate) is never a representative. It reads the fit of its
law at `L = 0`, its active law's key (`fit_key`), from `memo` or from this batch's
representatives, through its own order and with `L = 0`, and its row sets `L_phantom`. A
phantom whose law at `L = 0` has no fit there is fitted as an ordinary candidate; one
whose source's fit throws fails with the source's error.
```

In `identify_rate_equation`'s docstring, at the end of the "# Model selection (LOOCV)"
section (after "…to within the best one's fold-to-fold standard error."), add:

```
A candidate's fold losses come from fitting, on each fold, the mechanism whose fit it
reads. An L phantom (`L_phantom` in the CSVs: a one-subunit variant whose `L` the data
cannot determine) takes its source's fold losses, so beside its source it loses the
tie by its extra parameter, and selected without it, it reports `L = 0`.
```

- [ ] **Step 12: Run the identify tests to verify they pass**

Run (background) Step 10's command with log `$S/phase45/logs/P5_T10_tid.log`.
Expected: every testset passes, pristine, `exit 0`.

- [ ] **Step 13: Update the enumeration docstring, the spec and the ROADMAP**

In `src/mechanism_enumeration.jl`, in `_expand_to_allosteric`'s docstring, after the K-type
bullet's sentence "…in a random-order merged bi-bi, whose A group
takes A up at E and at E(B) under one rate constant, `L` shows in every K-type
variant.", add:

```
    identify keeps every K-type variant in the frontier, since its children can carry
    laws its parent's children lack. A variant that `_is_phantom` (`src/key.jl`) flags
    is not fitted when another mechanism holds the fit of its law at `L = 0`: it reads
    that fit, with `L = 0`.
```

(indented as the bullet's other continuation lines), rewrapping that bullet to 92
columns.

In `docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md` (P5-30 and P5-32,
both ◆ for Denis), section 13, replace the whole bullet that starts "**One-subunit L
phantoms** stay in the frontier, inherit the parent's fit and are" and ends "A variant
whose inactive conformation binds something is fitted.", that is

```
- **One-subunit L phantoms** stay in the frontier, inherit the parent's fit and are
  flagged (phase 5). A variant is an L phantom when its inactive conformation binds nothing
  and does not catalyze, and an exact rank test in `src/` (exact arithmetic, such as
  modular) shows that adding L does not raise the identifiable rank over the parent's. As
  L → 0 the variant's law becomes the parent's exactly, so its reported fit is the
  parent's constants under the variant's names with L = 0, and the parent's loss. In the
  beam and in cross-validation it counts its own `length(fitted_params)` (B9).
  Cross-validation does not refit it either: on every fold its best fit equals the
  parent's, so it takes the parent's fold losses. Its CSV row records the fit source's
  `mechanism_id` and the L-phantom flag (section 5). A variant whose inactive conformation
  binds something is fitted.
```

with

```
- **One-subunit L phantoms** stay in the frontier, inherit the fit of the mechanism that
  holds the key of their law at L = 0, their source (often the parent), and are flagged
  (phase 5). A variant is an L phantom when its inactive conformation binds nothing and
  does not catalyze and, inside that pre-filter, an exact rescaling certificate in `src/`
  (`key.jl`) shows that its family is the family of its law at L = 0; identifiable ranks
  are test oracles. As L → 0 the variant's law becomes the source's, up to the names of
  its constants, so its reported fit is the source's constants under the variant's names
  with L = 0, and the source's loss. In the beam and in cross-validation it counts its
  own `length(fitted_params)` (B9). Cross-validation does not refit it either: on every
  fold its best fit equals the source's, so it takes the source's fold losses. Its CSV row
  records the fit source's `mechanism_id` and the L-phantom flag (section 5). A phantom
  whose law at L = 0 no mechanism of the search holds is fitted as an ordinary candidate,
  unflagged, and a variant whose inactive conformation binds something is fitted.
```

In section 10, "identify", replace the bullet

```
- One-subunit L phantoms (section 13) inherit the parent's fit and are flagged; the parent
  has another key, so section 13 defines how the fit maps onto the variant.
```

with

```
- One-subunit L phantoms (section 13) inherit the fit of the mechanism that holds the key
  of their law at L = 0 (often the parent) and are flagged; that mechanism has another
  key, so section 13 defines how the fit maps onto the variant.
```

Run `grep -n "parent's fit\|parent's loss\|parent's fold" docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md`:
it prints nothing.

In `ROADMAP.md`:
- Under "Planned", delete the bullet "**Flag one-subunit K-type variants whose L is a
  phantom instead of fitting them.**" through its last line,
  "`docs/src/identify/enumeration_engine.md` documents the phantom. Priority: medium."
- Under "Planned", in "**Profile the first v0.9.0 HPC runs.**", replace
  ```
    zero-gain `change_allo_state` children, and how often PGK meets the one-subunit L
    phantom. These numbers rank the open search-speed and deduplication ideas.
  ```
  with
  ```
    zero-gain `change_allo_state` children, and how many PGK rows identify flags
    `L_phantom`. These numbers rank the open search-speed and deduplication ideas.
  ```
- Under "Ideas", replace the whole bullet "**Count identifiable parameters with an exact
  modular rank.**" with:
  ```
  - **Count identifiable parameters with an exact modular rank.**
    Several phantom classes make the beam count a model above its true dimension: chain
    flanks kept after a later split, copy-copy twin pairs, `:NonequalAI` copies and the
    iso dwell-time gauge. The one-subunit L has an exact test of its own: identify flags
    an L phantom whose inactive conformation binds nothing by a rescaling certificate
    (`_phantom_weights`, `src/key.jl`) and lets it read the fit of its law at L = 0. The
    rank of ∂v/∂log θ modulo a prime at random points of the derived law (`_law`,
    `src/engine.jl`) is exact up to a negligible failure probability, took 3.5 ms per
    mechanism in a prototype, and has a finite-difference oracle in
    `_testhelper_identifiable_rank`
    (`test/mechanism_definitions_for_test_enzyme_derivation.jl`). Run it over several
    saved runs to find the other phantom classes; performance may keep it out of `src/`,
    where it would otherwise replace the beam's parameter count. It must never filter
    children, since a phantom child can carry families its parent lacks.
  ```
- Under "Ideas", after that bullet, add:
  ```
  - **Apply the chain-flank exclusion at one catalytic subunit.**
    The flip move keeps a chain's flanks at rapid equilibrium only in plain mechanisms
    (`_chain_flank_groups`, `src/mechanism_enumeration.jl`): an allosteric mechanism has
    none, since above one subunit a flank's flip can show through the two
    conformations. At one subunit an L phantom's flank flip stays invisible, so the move
    emits a phantom whose law at L = 0 no plain mechanism of the search holds, and
    identify fits it as an ordinary candidate (census #3268; the phase-5 research found
    18 such bi-bi variants on 8 keys in a sampled round 2). Applying the exclusion to
    one-subunit variants would drop these fits; first check that no other law becomes
    unreachable.
  ```

- [ ] **Step 14: Rerun R5.2 on P4-1's ping-pong tagging children**

Create `$S/phase45/scripts/P5_phantom_rerun.jl`:

```julia
# ABOUTME: R5.2 on P4-1's ping-pong tagging children at one catalytic subunit: each flagged
# ABOUTME: L phantom against its law at L = 0 (key, count, map, rank); counts and costs.
# Usage, from the worktree root: julia --project=. P5_phantom_rerun.jl
using TestEnv; TestEnv.activate()
using Test, EnzymeRates, LinearAlgebra, Random, Printf
include(joinpath(pkgdir(EnzymeRates), "test",
                 "mechanism_definitions_for_test_enzyme_derivation.jl"))
const ER = EnzymeRates
const MechU = Union{ER.Mechanism, ER.AllostericMechanism}
const RXN = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    allosteric_regulators: R
    shared_catalytic_site: (A, P)
end
covalent(m) = any(ER._has_residual, ER.catalytic(m).forms)
key_of(L, m) = ER._canonical(L, ER._mechanism_id(sprint(show, m))).key
"""The median and the maximum of `ts`, in µs."""
spread(ts) = (1e6 * sort(ts)[cld(length(ts), 2)], 1e6 * maximum(ts))
"""Check every one-subunit allosteric variant of `kids` with a covalent form: count the
flagged phantoms and the certificate-only variants (a certificate, but the pre-filter fits
them), test each flagged one against its law at L = 0, and time the detector on each
variant: `_is_phantom`, as PASS1 calls it, and the certificate `_phantom_weights`."""
function check(label, kids, rng)
    t0 = time()
    variants = [am for am in kids if am isa ER.AllostericMechanism &&
                ER.catalytic_multiplicity(am) == 1 && covalent(am)]
    # Compile both before timing them.
    isempty(variants) || ER._is_phantom(first(variants), ER._law(first(variants)))
    isempty(variants) || ER._phantom_weights(ER._law(first(variants)))
    flagged, certonly, bad = 0, 0, String[]
    t_detect, t_cert = Float64[], Float64[]
    for am in variants
        s = ER._solve(am)
        A = ER._law(am, s)
        push!(t_detect, @elapsed ER._is_phantom(am, A))
        push!(t_cert, @elapsed w = ER._phantom_weights(A))
        isempty(w) && continue
        if !ER._is_phantom(am, A)
            certonly += 1
            continue
        end
        flagged += 1
        cm = ER.catalytic(am)
        ok = key_of(A.active, am) == key_of(ER._law(cm), cm) &&
             length(s.fitted) == length(ER._solve(cm).fitted) + 1 &&
             _testhelper_phantom_map_holds(am; rng) &&
             _testhelper_identifiable_rank(am) == _testhelper_identifiable_rank(cm)
        ok || push!(bad, sprint(show, am))
    end
    @printf("%s: %d one-subunit covalent variants, %d flagged, %d certificate-only, ",
            label, length(variants), flagged, certonly)
    @printf("%d failing; %.0f s\n", length(bad), time() - t0)
    if !isempty(variants)
        @printf("%s: _is_phantom median %.1f µs, max %.1f µs; ", label, spread(t_detect)...)
        @printf("_phantom_weights median %.1f µs, max %.1f µs\n", spread(t_cert)...)
    end
    foreach(println, bad)
    flush(stdout)
    isempty(bad)
end
seeds = unique!(collect(ER.init_mechanisms(RXN)))
base, _ = ER._base_tier(seeds, RXN)
round1 = unique!(ER.expand_mechanisms(MechU[m for m in base if covalent(m)], RXN))
ok1 = check("round 1", round1, MersenneTwister(1))
pool = MechU[m for m in round1 if covalent(m)]
parents2 = length(pool) <= 200 ? pool :
           pool[sort(randperm(Xoshiro(60), length(pool))[1:200])]
round2 = unique!(ER.expand_mechanisms(parents2, RXN))
ok2 = check("round 2, $(length(parents2)) sampled parents", round2, MersenneTwister(2))
println(ok1 && ok2 ? "PASS" : "FAIL")
```

Run (background, about 20 min): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T10_rerun.log 2400 --project=. --heap-size-hint=2500M $S/phase45/scripts/P5_phantom_rerun.jl`,
polled in bounded loops.
Expected: two count lines, each followed by its timing line, `PASS`, `exit 0`, maxrss
under 3 GB. Q6 found no flagged ping-pong child among the tagging children of 30 parents,
so round 1 may flag none; round 2 adds the other moves' children. Record the count lines
for the PR body (P5-36) and the timing lines for its Measurements (spec section 17; the
research measured about 25 µs per certificate on slice code). A `FAIL`
names each failing variant: stop and debug it (systematic-debugging) before the commit.

- [ ] **Step 15: Run the full suite**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T10_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`,
polled in bounded loops.
Expected: `exit 0`, every test passing, pristine, Aqua and JET green; about 12 s over Task
9's time (Q5 and Q6 about 8 s, the identify phantom tests the rest); report `maxrss_kb`
against Task 0's.

- [ ] **Step 16: Commit**

Write the measured figures in place of `<…>`.

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt
git add src/key.jl src/identify_rate_equation.jl src/mechanism_enumeration.jl \
    docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md ROADMAP.md \
    test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_key.jl \
    test/test_mechanism_enumeration.jl test/test_oracles.jl \
    test/test_identify_rate_equation.jl
git commit -m "$(cat <<'EOF'
Flag one-subunit L phantoms and let them inherit a fit

A one-subunit variant whose inactive conformation binds nothing, and whose
law has a rescaling certificate, has the family of its law at L = 0. It
stays in the frontier, reads the fit of whoever holds that law's key with
L = 0, and takes the source's fold losses; one without such a source is
fitted unflagged. The CSV gains L_phantom. The spec's pre-filter stays, and
the certificate replaces the exact rank of spec section 13, which now says
so; sections 10 and 13 name the fit's source as the holder of the phantom's
L = 0 key.

The K-type verdicts were re-derived under the pre-filter and checked
against the rank oracle and the map oracle before pinning. Q5: <n5> of 200
bi-bi variants flagged, every map holds. Q6: no ping-pong tagging child of
30 parents binds nothing in the inactive conformation or has a
certificate. R5.2 rerun on ping-pong: <round 1 line>; <round 2 line>.
_is_phantom per variant: median <us> µs, max <us> µs (_phantom_weights:
median <us> µs, max <us> µs). Full-suite maxrss: <kb> KB (Task 0: <kb0> KB).

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
)"
```

---
### Task 11: The unit rule (P5-9 to P5-13, P5-10a)

Every denominator term of a law shares one unit. The unit rule picks, per conformation, a
normalizer N with that unit, so each term divided by N is dimensionless: the constant
term; else weights on a graded basis of terms, fractional ones allowed; else the first
term, monomial included (the ratio form). There is no Cleland-first rule (P5-10a). This
task computes N; Task 12 prints with it. A units oracle read from the graph checks every
law of the population sample and the 42 specs.

**Files:**
- Create: `src/normalize.jl`; modify `src/EnzymeRates.jl` (include it after `engine.jl`)
- Create: `test/test_print.jl`; modify `test/runtests.jl` (include it after
  `test_key.jl`)
- Modify: the fixtures file (`_testhelper_name_units`, `_testhelper_term_units`,
  `_testhelper_check_units`, `_testhelper_norm_value`)
- Test: `test/test_oracles.jl` (Q1's units)
- Logs: `$S/phase45/logs/P5_T11_{print_red,print,oracles,full}.log`

**Interfaces:**
- Consumes: Task 2's `_Mono`, `_Exp`, `_mexp`, `_monomial`, `_split_conc`, `_Vars`; Tasks
  3–4's `_law`, `_Law`, `_AlloLaw`; `_rref!(M::Matrix{<:Rational}, n::Int)` (`thermo.jl`).
  Tests: Task 3's `_testhelper_poly_value(p, vars, env)`; Task 4's
  `_testhelper_evaluate_law(::_AlloLaw, env)`; `_testhelper_law_inputs`,
  `_testhelper_graph`, `_testhelper_is_binding`, `_testhelper_lift`,
  `_testhelper_single_subunit`, `MECHANISM_TEST_SPECS`.
- Produces:
  ```julia
  struct _Norm
      rule::Symbol                       # :constant, :graded or :ratio
      lambda::Vector{Pair{_Mono, _Exp}}  # concentration monomial of D => weight, in basis
                                         #   order, zero weights left out
      shift::_Mono                       # the ratio form's monomial; empty otherwise
  end
  _term_key(mo::_Mono, vars::_Vars, substrates::Vector{Int})
      # (degree, distinct non-substrates, -exponent vector over the concentrations)
  _normalizer(D::_Poly, vars::_Vars, substrates::Vector{Int})::_Norm
  # tests (fixtures file)
  _testhelper_name_units(m)::Dict{Symbol, Tuple{Rational{Int}, Rational{Int}}}
  _testhelper_term_units(units::AbstractDict, vars::EnzymeRates._Vars,
                         mo::EnzymeRates._Mono)::Tuple{Rational{Int}, Rational{Int}}
  _testhelper_check_units(m, L)::Nothing           # @test inside
  _testhelper_norm_value(nz, D, vars, env)::BigFloat
  ```
  An MWC law calls `_normalizer` once per conformation, on `A.active.D` over
  `A.active.vars` and on `A.inactive.D` over `A.inactive.vars`; both number their
  concentrations alike, so one `substrates` vector serves both.
  - Test helpers in `test_print.jl`: `_testhelper_normalizer(m)`,
    `_testhelper_conc_mono(L, ps...)`.

**Gate:**
- Every fixture takes its rule and its exact weights.
- Q1's units hold over the 120 sample laws and the 42 specs: every term of each
  conformation's D shares one unit, the normalizer has it, and `H·Q/D` is a rate per
  time.
- The full suite is green and pristine; about +2 s.

Facts this task rests on (`$S/phase45/plan5/probes/p3_probe.log`, sections E–G, on
a9e0699 with Tasks 1–4's law code and this task's code): each fixture below took the rule
and weights its testset pins; the guard admits F1–F4 and the power case; each MWC
fixture's normalized combination equals its law to 2e-77; the units check passed on the
42 specs (111 tests, 1.9 s) and, in `p3c_probe.log` section D, on the 120 sampled laws at
one subunit (346 tests, 1.2 s). Two figures of RES §3.2 concern the printed N, not the
normalizer: "N is the product of two sums" (Segel ter bi) and "graded N = Cleland's
squared" (the power case) hold for N's polynomial factor once Task 12's printer moves the
integer and monomial content into the coefficients (the `N =` lines of both in
`$S/phase45/research5/unit/logs/examples.log`). On the full normalizers the probe finds
ratios 1.41 and 0.129, so this task pins the weights and leaves those two figures to
Task 12.

- [ ] **Step 1: Add the units oracle and the normalizer's value to the fixtures file**

In `test/mechanism_definitions_for_test_enzyme_derivation.jl`, after Task 4's
`_testhelper_evaluate_law(A::EnzymeRates._AlloLaw, env)`, add:

```julia
"""
    _testhelper_name_units(m) -> Dict{Symbol, Tuple{Rational{Int}, Rational{Int}}}

The units, as (concentration, time) exponents, of every name a rate law of `m` (concrete or
compiled) reads, from the graph alone (`_testhelper_graph`), never from the law: each
concentration and `E_total` (1, 0); `Keq` (product copies − substrate copies, 0); a
steady-state rate constant (−copies its step takes up in its direction, −1); a
rapid-equilibrium binding's dissociation constant (1, 0); any other rapid-equilibrium
constant (copies released − copies consumed, 0); a regulator's dissociation constant
(1, 0); `L` (0, 0).
"""
function _testhelper_name_units(m)
    mech = _testhelper_lift(m)
    g = _testhelper_graph(mech)
    u = Dict{Symbol, Tuple{Rational{Int}, Rational{Int}}}(
        x => (1 // 1, 0 // 1) for x in EnzymeRates.metabolites(mech))
    u[:E_total] = (1 // 1, 0 // 1)
    u[:Keq] = (length(g.products) - length(g.substrates), 0 // 1)
    for st in g.steps, names in (st.names, st.inactive_names), (slot, x) in enumerate(names)
        u[x] = !st.rapid_equilibrium ?
                   (-(slot == 1 ? length(st.consumed) : length(st.released)), -1 // 1) :
               _testhelper_is_binding(st) ? (1 // 1, 0 // 1) :
                   (length(st.released) - length(st.consumed), 0 // 1)
    end
    if mech isa EnzymeRates.AllostericMechanism
        for state in (:A, :I), p in EnzymeRates._kreg_params(mech, state)
            u[EnzymeRates.name(p, mech)] = (1 // 1, 0 // 1)
        end
        u[:L] = (0 // 1, 0 // 1)
    end
    u
end

"""The units of the monomial `mo` over `vars`, from the name units `units`
(`_testhelper_name_units`): the sum of each variable's units times its exponent."""
_testhelper_term_units(units::AbstractDict, vars::EnzymeRates._Vars,
                       mo::EnzymeRates._Mono) =
    foldl(((c, t), (v, e)) -> (c + e * units[vars.names[v]][1],
                               t + e * units[vars.names[v]][2]), mo;
          init = (0 // 1, 0 // 1))

"""
    _testhelper_check_units(m, L)

Test the units of `L`, the derived law of the concrete mechanism `m` (a `_Law`, or an
`_AlloLaw` for an allosteric `m`), against the graph's units (`_testhelper_name_units`).
In each conformation every term of `D` has one unit, the unit rule's normalizer
(`EnzymeRates._normalizer`) has that unit too, so each denominator term over `N` is
dimensionless, and a numerator that is not zero makes `H·Q/D` a rate per time. In an MWC
law `fI·D_A` and `fA·D_I` share a unit, so the two conformations' terms add.
"""
function _testhelper_check_units(m, L)
    units = _testhelper_name_units(m)
    vars = L isa EnzymeRates._Law ? L.vars : L.active.vars
    subs = Int[vars.index[x] for x in unique(EnzymeRates.substrates(m))]
    unit(p) = only(unique(_testhelper_term_units(units, vars, mo) for mo in keys(p.terms)))
    plus(a, b) = (a[1] + b[1], a[2] + b[2])
    for C in (L isa EnzymeRates._Law ? (L,) : (L.active, L.inactive))
        uD = unit(C.D)
        nz = EnzymeRates._normalizer(C.D, vars, subs)
        Dc = EnzymeRates._split_conc(C.D, vars)
        uN = foldl((acc, (mo, λ)) -> plus(acc, λ .* unit(Dc[mo])), nz.lambda;
                   init = _testhelper_term_units(units, vars, nz.shift))
        @test uN == uD
        iszero(C.Q) || @test plus(unit(L.H), unit(C.Q)) == plus(uD, (0 // 1, -1 // 1))
    end
    L isa EnzymeRates._AlloLaw &&
        @test plus(unit(L.fI), unit(L.active.D)) == plus(unit(L.fA), unit(L.inactive.D))
    nothing
end

"""The value at `env` (by name), in BigFloat, of the normalizer `nz` of the denominator
`D` over `vars`: Π c(m)^λ_m over `nz.lambda`, `c(m)` the coefficient of the concentration
monomial `m` in `D`, times the monomial `nz.shift`."""
function _testhelper_norm_value(nz, D, vars, env)
    Dc = EnzymeRates._split_conc(D, vars)
    x(v) = BigFloat(env[vars.names[v]])
    prod((_testhelper_poly_value(Dc[mo], vars, env)^BigFloat(λ) for (mo, λ) in nz.lambda);
         init = big(1)) * prod((x(v)^BigFloat(e) for (v, e) in nz.shift); init = big(1))
end
```

`unit` reads `only`, so a denominator or coefficient whose terms disagree in units errors
there, naming the law; `_testhelper_name_units` reads the step tables of the graph, whose
constant names are the solve's (`name(p, m, kept)`).

- [ ] **Step 2: Write the failing tests**

Create `test/test_print.jl`:

```julia
# ABOUTME: Tests of what a printed rate law rests on: the unit rule's normalizer of each
# ABOUTME: conformation's denominator, and the units of every derived law from the graph.

# The unit rule's normalizer of the plain mechanism `m`'s law, and the law.
function _testhelper_normalizer(m)
    mech = _testhelper_lift(m)
    L = EnzymeRates._law(mech)
    subs = Int[L.vars.index[x] for x in unique(EnzymeRates.substrates(mech))]
    EnzymeRates._normalizer(L.D, L.vars, subs), L
end

# The concentration monomial of law `L` with the exponents `x => e`.
_testhelper_conc_mono(L, ps::Pair{Symbol}...) =
    EnzymeRates._monomial((L.vars.index[x] => e for (x, e) in ps)...)

# ── The unit rule (spec section 9) ───────────────────────────────────────────

@testset "Michaelis–Menten takes its constant term" begin
    # Rule 1: the denominator has a term free of concentrations, the free enzyme's, and
    # N is its coefficient, so every term reads as a form's abundance relative to free E.
    mm = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E + P
        end
    end
    nz, L = _testhelper_normalizer(mm)
    @test nz.rule === :constant
    @test nz.lambda == [EnzymeRates._Mono() => 1 // 1]
    @test isempty(nz.shift)
    _testhelper_check_units(_testhelper_lift(mm), L)
end

@testset "Segel Ping Pong Bi Bi takes the graded basis A, B, A·B" begin
    # Rule 2: no constant term. A, B, P and Q (degree 1), then A·B, are affinely
    # independent; the weights A:1, B:1, A·B:−1 solve Σλ = 1 and Σλ·a = 0, and P and Q
    # weigh 0. Here the graded basis is Cleland's N: ΠS/X over each substrate X, over ΠS.
    pp = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) <--> F + P
            F + B <--> F(B)
            F(B) <--> E + Q
        end
    end
    nz, L = _testhelper_normalizer(pp)
    m(ps...) = _testhelper_conc_mono(L, ps...)
    @test nz.rule === :graded
    @test nz.lambda == [m(:A => 1) => 1 // 1, m(:B => 1) => 1 // 1,
                        m(:A => 1, :B => 1) => -1 // 1]
    @test isempty(nz.shift)
    _testhelper_check_units(_testhelper_lift(pp), L)
end

@testset "Segel Bi Uni Uni Uni PP Ter Bi takes the graded basis with C" begin
    # The degree-1 term C enters the basis before Cleland's monomials A·B, A·C, B·C, so
    # the weights are C:2, A·B:1, A·C:−1, B·C:−1. P5-10a keeps this graded basis: N is
    # not Cleland's N = c(B·C)·c(A·C)·c(A·B)/c(A·B·C)².
    terbi = @enzyme_mechanism begin
        substrates: A, B, C
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> F + P
            F + C <--> F(C)
            F(C) <--> E + Q
        end
    end
    nz, L = _testhelper_normalizer(terbi)
    m(ps...) = _testhelper_conc_mono(L, ps...)
    @test nz.rule === :graded
    @test nz.lambda == [m(:C => 1) => 2 // 1, m(:A => 1, :B => 1) => 1 // 1,
                        m(:A => 1, :C => 1) => -1 // 1, m(:B => 1, :C => 1) => -1 // 1]
    rng = Random.MersenneTwister(31)
    env = Dict{Symbol, BigFloat}(x => exp(2 * big(rand(rng)) - 1) for x in L.vars.names)
    Dc = EnzymeRates._split_conc(L.D, L.vars)
    c(ps...) = _testhelper_poly_value(Dc[m(ps...)], L.vars, env)
    cleland = c(:B => 1, :C => 1) * c(:A => 1, :C => 1) * c(:A => 1, :B => 1) /
              c(:A => 1, :B => 1, :C => 1)^2
    @test !isapprox(_testhelper_norm_value(nz, L.D, L.vars, env), cleland; rtol = 1e-6)
    _testhelper_check_units(_testhelper_lift(terbi), L)
end

@testset "the ping-pong power case puts A² before A·B" begin
    # Census ping-pong mechanism 6e46f55b0b68e5d9 (UNIT §2). A² and A·B are both degree
    # 2 with no non-substrate; the earlier name's larger exponent comes first (P5-11), so
    # A² joins A, B, P and Q in the basis, and the weights are A:2, A²:−1.
    power = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> E(P; residual = A - P), E(Q) + A <--> E(A, Q))
            (E + Q <--> E(Q), E(P; residual = A - P) + Q <--> E(A, Q))
            (E(; residual = A - P) + B <--> E(Q),
             E(P; residual = A - P) + B <--> E(B, P; residual = A - P))
            (E(; residual = A - P) + P ⇌ E(P; residual = A - P),
             E(Q) + P ⇌ E(B, P; residual = A - P))
        end
    end
    nz, L = _testhelper_normalizer(power)
    m(ps...) = _testhelper_conc_mono(L, ps...)
    @test nz.rule === :graded
    @test nz.lambda == [m(:A => 1) => 2 // 1, m(:A => 2) => -1 // 1]
    _testhelper_check_units(_testhelper_lift(power), L)
end

@testset "F1 divides by its first term, A" begin
    # Rule 3: two half-reactions between the free forms E and F, each A → P. Every term
    # holds A or P, and no weights on A and P solve the units equation, so the law reads A
    # and P only through their ratio: N is the first term, its monomial A included.
    f1 = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
        end
    end
    nz, L = _testhelper_normalizer(f1)
    a = _testhelper_conc_mono(L, :A => 1)
    @test nz.rule === :ratio
    @test nz.lambda == [a => 1 // 1] && nz.shift == a
    _testhelper_check_units(_testhelper_lift(f1), L)
end

@testset "F2 takes fractional weights" begin
    # F1 with F binding two A at rapid equilibrium: the terms A, P, A³ and A²·P. The basis
    # A, P, A³ solves the units equation with A:3/2, A³:−1/2 (P5-12).
    f2 = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
            F + A + A ⇌ F(A, A)
        end
    end
    nz, L = _testhelper_normalizer(f2)
    m(ps...) = _testhelper_conc_mono(L, ps...)
    @test nz.rule === :graded
    @test nz.lambda == [m(:A => 1) => 3 // 2, m(:A => 3) => -1 // 2]
    _testhelper_check_units(_testhelper_lift(f2), L)
end

@testset "F3, a rapid-equilibrium Theorell–Chance step, divides by B" begin
    # Every term holds B or P (the rapid-equilibrium step E(A) + B ⇌ E(Q) + P ties them),
    # and the units equation has no solution on the basis B, P, A·B, B·Q, so N is the
    # first term, B with its coefficient (P5-13).
    f3 = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B ⇌ E(Q) + P
            E(Q) <--> E + Q
        end
    end
    nz, L = _testhelper_normalizer(f3)
    b = _testhelper_conc_mono(L, :B => 1)
    @test nz.rule === :ratio
    @test nz.lambda == [b => 1 // 1] && nz.shift == b
    _testhelper_check_units(_testhelper_lift(f3), L)
end

@testset "F4 takes integer weights A:2, A²:−1" begin
    # F1 with F binding A twice in two rapid-equilibrium steps: the terms A, P, A², A·P,
    # A³ and A²·P. The basis A, P, A² solves the units equation with A:2, A²:−1.
    f4 = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
            F + A ⇌ F(A)
            F(A) + A ⇌ F(A, A)
        end
    end
    nz, L = _testhelper_normalizer(f4)
    m(ps...) = _testhelper_conc_mono(L, ps...)
    @test nz.rule === :graded
    @test nz.lambda == [m(:A => 1) => 2 // 1, m(:A => 2) => -1 // 1]
    _testhelper_check_units(_testhelper_lift(f4), L)
end

@testset "each MWC conformation takes its own rule" begin
    # Spec section 9: an MWC law normalizes each conformation's polynomials by that
    # conformation's rule, dividing its Q, D and free-E weight alike, so the MWC
    # combination does not change. A two-subunit ping-pong whose chemistry differs between
    # the conformations takes Segel's graded basis in both. F3 with A's binding and the
    # Theorell–Chance step :OnlyA takes the ratio form in the active conformation and the
    # constant term in the inactive one, which keeps free E and E(Q).
    pp2 = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        catalytic_multiplicity: 2
        steps: begin
            E + A <--> E(A)       :: EqualAI
            E(A) <--> F + P       :: NonequalAI
            F + B <--> F(B)       :: EqualAI
            F(B) <--> E + Q       :: EqualAI
        end
    end)
    tc = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)            :: OnlyA
            E(A) + B ⇌ E(Q) + P        :: OnlyA
            E(Q) <--> E + Q            :: EqualAI
        end
    end)
    for (i, am) in enumerate((pp2, tc))
        s = EnzymeRates._solve(am)
        A = EnzymeRates._law(am, s)
        subs = Int[A.vars.index[x] for x in unique(EnzymeRates.substrates(am))]
        nA = EnzymeRates._normalizer(A.active.D, A.active.vars, subs)
        nI = EnzymeRates._normalizer(A.inactive.D, A.inactive.vars, subs)
        m(ps...) = _testhelper_conc_mono(A.active, ps...)
        if i == 1
            segel = [m(:A => 1) => 1 // 1, m(:B => 1) => 1 // 1,
                     m(:A => 1, :B => 1) => -1 // 1]
            @test (nA.rule, nI.rule) == (:graded, :graded)
            @test nA.lambda == nI.lambda == segel
        else
            @test (nA.rule, nI.rule) == (:ratio, :constant)
            @test nA.lambda == [m(:B => 1) => 1 // 1] && nA.shift == m(:B => 1)
            @test nI.lambda == [EnzymeRates._Mono() => 1 // 1]
        end
        _testhelper_check_units(am, A)
        rng = Random.MersenneTwister(32 + i)
        for _ in 1:3
            _, concs, params = _testhelper_law_inputs(am, s.fitted, rng)
            env = merge(params, concs)
            value(p) = _testhelper_poly_value(p, A.vars, env)
            NA = _testhelper_norm_value(nA, A.active.D, A.active.vars, env)
            NI = _testhelper_norm_value(nI, A.inactive.D, A.inactive.vars, env)
            n, l = A.n, env[:L]
            qa, da, fa = value(A.active.Q) / NA, value(A.active.D) / NA, value(A.fA) / NA
            qi, di, fi = (value(A.inactive.Q) / NI, value(A.inactive.D) / NI,
                          value(A.fI) / NI)
            num = fi^n * qa * da^(n - 1) + l * fa^n * qi * di^(n - 1)
            normalized = value(A.H) * num / (fi^n * da^n + l * fa^n * di^n)
            law = _testhelper_evaluate_law(A, env)
            @test abs(normalized - law) <= big(10)^-20 * abs(law)
        end
    end
end

@testset "every spec's law has one unit per denominator and a rate per time" begin
    # Q1's units on the 42 specs, each at one catalytic subunit (spec section 14).
    for spec in MECHANISM_TEST_SPECS
        @testset "$(spec.name)" begin
            am = _testhelper_lift(_testhelper_single_subunit(spec.mechanism))
            _testhelper_check_units(am, EnzymeRates._law(am, EnzymeRates._solve(am)))
        end
    end
end
```

In `test/runtests.jl`, after the line `include("test_key.jl")` (Task 7), add
`    include("test_print.jl")`.

In `test/test_oracles.jl`, in Task 4's testset "the derived laws of the population sample
are exact, positive and balanced", replace

```julia
    # Q1 (spec section 14). Each sampled law, derived from the concrete mechanism at one
    # catalytic subunit, equals the exact steady state at two points, holds positive
    # coefficients in each conformation's Q, and vanishes at equilibrium, per conformation
    # and combined.
```

with

```julia
    # Q1 (spec section 14). Each sampled law, derived from the concrete mechanism at one
    # catalytic subunit, equals the exact steady state at two points, holds positive
    # coefficients in each conformation's Q, vanishes at equilibrium, per conformation
    # and combined, and has one unit per denominator, which its normalizer shares
    # (`_testhelper_check_units`).
```

and, in its inner testset "population mechanism $i", replace

```julia
            @test all(<(big(10)^-30),
                      _testhelper_detailed_balance(m, Random.MersenneTwister(i)))
        end
```

with

```julia
            @test all(<(big(10)^-30),
                      _testhelper_detailed_balance(m, Random.MersenneTwister(i)))
            _testhelper_check_units(mech, L)
        end
```

(`mech` is the one-subunit mechanism and `L` its law, both bound earlier in that testset.)

- [ ] **Step 3: Run the printing tests to verify they fail**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T11_print_red.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); @testset "red" begin include("test/test_print.jl") end'`
Expected: every testset errors with ``UndefVarError: `_normalizer` not defined in
`EnzymeRates` ``; the outer `red` testset ends the log with `exit 1`.

- [ ] **Step 4: Write the unit rule**

Create `src/normalize.jl`:

```julia
# ABOUTME: The unit rule of a rate law: the normalizer that makes every denominator term
# ABOUTME: dimensionless, by the constant term, a graded basis of terms or the ratio form.

"""
The normalizer of one conformation's denominator (`_normalizer`): `N = Π_m c(m)^λ_m ·
x^shift` over the pairs `m => λ_m` of `lambda`, `c(m)` being the coefficient of the
concentration monomial `m` in the denominator. `rule` names the case that chose it:
`:constant`, the constant term; `:graded`, weights on the first affinely independent
terms; `:ratio`, the first term, whose concentration monomial `shift` stays in `N`.
`shift` is empty under the other two rules.
"""
struct _Norm
    rule::Symbol
    lambda::Vector{Pair{_Mono, _Exp}}
    shift::_Mono
end

"""
    _term_key(mo, vars, substrates) -> key

The unit rule's order on the concentration monomials of a denominator (spec section 9):
total degree first, then the number of distinct concentrations that are not substrates
(`substrates` holds the substrates' variables), then the exponent vector over the
concentrations of `vars`, in `vars` order, with an earlier name's larger exponent first,
so `A²` comes before `A·B`.
"""
function _term_key(mo::_Mono, vars::_Vars, substrates::Vector{Int})
    conc = findall(vars.conc)
    ex = _Exp[_mexp(mo, v) for v in conc]
    (sum(ex), count(i -> !iszero(ex[i]) && !(conc[i] in substrates), eachindex(conc)), -ex)
end

"""
    _normalizer(D, vars, substrates) -> _Norm

The normalizer of the denominator `D` over `vars` (spec section 9, "Normalization: the
unit rule"). Every term of `D` has one unit, and every term divided by `N` is
dimensionless. The first rule that applies gives `N`:

1. `:constant`: `D` has a term free of concentrations, and `N` is its coefficient.
2. `:graded`: taking the concentration monomials `m` of `D` in `_term_key` order, the
   first ones whose vectors `[1; a_m]` (`a_m` the exponent vector) are affinely
   independent form a basis; when weights on the basis satisfy `Σ λ_m = 1` and
   `Σ λ_m·a_m = 0`, `N = Π c(m)^λ_m`. The weights are exact and may be fractional.
3. `:ratio`: otherwise the law reads some concentrations only through their ratio, and
   `N` is the first term in that order, its concentration monomial included.

Exact `Rational{Int}` arithmetic; an `OverflowError` escapes.
"""
function _normalizer(D::_Poly, vars::_Vars, substrates::Vector{Int})
    Dc = _split_conc(D, vars)
    haskey(Dc, _Mono()) && return _Norm(:constant, Pair{_Mono, _Exp}[_Mono() => 1], _Mono())
    conc = findall(vars.conc)
    monos = sort!(collect(keys(Dc)); by = mo -> _term_key(mo, vars, substrates))
    column(mo) = _Exp[1; [_mexp(mo, v) for v in conc]]
    basis = _Mono[]
    for mo in monos
        M = reduce(hcat, [column(b) for b in [basis; [mo]]])
        length(_rref!(M, length(basis) + 1)) > length(basis) && push!(basis, mo)
    end
    k = length(basis)
    M = [reduce(hcat, [column(b) for b in basis]) _Exp[1; zeros(Int, length(conc))]]
    _rref!(M, k)
    all(iszero, view(M, (k + 1):size(M, 1), k + 1)) ||
        return _Norm(:ratio, Pair{_Mono, _Exp}[first(monos) => 1], first(monos))
    _Norm(:graded, Pair{_Mono, _Exp}[b => M[i, k + 1] for (i, b) in enumerate(basis)
                                     if !iszero(M[i, k + 1])], _Mono())
end
```

The basis vectors are linearly independent, so `_rref!` puts their pivots on the first
`k` rows in basis order, and row `i`'s last entry is basis term `i`'s weight; a nonzero
entry below them means no weights exist. `_term_key` returns one concrete tuple type and
`_normalizer` one concrete struct (R18). In `src/EnzymeRates.jl`, add
`include("normalize.jl")` on the line after `include("engine.jl")`.

- [ ] **Step 5: Run the printing and population tests to verify they pass**

Rerun Step 3's command with log `$S/phase45/logs/P5_T11_print.log`. Expected: every
testset passes, pristine, `exit 0`; the 42 spec testsets among them.

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T11_oracles.log 3000 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_oracles.jl")'`
(background, polled).
Expected: every testset passes, pristine; Q1 runs its 120 population testsets with the
units check, about 1 s more than in Task 10's log.

- [ ] **Step 6: Run the full suite**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T11_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`,
polled in bounded loops.
Expected: `exit 0`, every test passing, pristine, Aqua and JET green; about 2 s over Task
10's time; report `maxrss_kb` against Task 0's.

- [ ] **Step 7: Commit**

Write the measured figures in place of `<…>`.

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt
git add src/normalize.jl src/EnzymeRates.jl test/test_print.jl test/runtests.jl \
    test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_oracles.jl
git commit -m "$(cat <<'EOF'
Normalize each conformation's denominator by the unit rule

The normalizer is the constant term, else weights on the first affinely
independent terms in degree order, fractional ones allowed, else the first
term with its monomial: the ratio form. There is no Cleland-first rule. An
MWC law normalizes each conformation by its own rule. A units oracle read
from the graph checks every term of the population sample's laws and the
42 specs: one unit per denominator, shared by its normalizer, and H·Q/D a
rate per time. Nothing prints differently yet. Full-suite maxrss: <kb> KB
(Task 0: <kb0> KB).

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
)"
```

---

### Task 12: The printer and coefficient form (P5-14 to P5-17)

`rate_equation_string` prints the law in the independent constants over the unit rule's
normalizer, its reference lines in ratio-of-products form with exact exponents, and, with
`coefficients = true`, the coefficient form. Value checks land first (test_print.jl), then
the switch, then each re-pin. PASS1 still prints through the compiled mechanism here; Task
12b moves it onto the concrete law.

**Files:**
- Create: `src/print.jl`; modify `src/EnzymeRates.jl` (include after `normalize.jl`)
- Modify: `src/rate_eq_derivation.jl` (delete `rate_equation_string`'s docstring, its
  `function rate_equation_string end`, its concrete method and its compiled method;
  `_equation_text` and `_assert_no_name_clash` stay until Task 13)
- Modify: `src/sym_poly_for_rate_eq_derivation.jl` (`build_power_expr`: a non-integer
  exponent is an exact `Rational{Int}` literal)
- Modify: `docs/src/deriving/textbooks.md` (the doctest's output; the prose is Task 15's)
- Modify: the fixtures file (the printed-law oracles; the 11 `expected_factored_num` /
  `expected_factored_denom` pairs and their field comment)
- Test: `test/test_print.jl` (four testsets appended), `test/test_rate_eq_derivation.jl`,
  `test/test_oracles.jl` (Q1; the half-integer dependent), `test/test_allosteric_golden.jl`,
  `test/reference/allosteric_golden_reference.txt`; run only:
  `test/test_allosteric_collapse.jl`, `test/test_mechanism_enumeration.jl`,
  `test/test_identify_rate_equation.jl`, `test/test_types.jl`, `test/test_dsl.jl`,
  `test/test_aqua_jet.jl`
- Create (scratch): `$S/phase45/scripts/P5_T12_pins.jl`, `P5_T12_golden.jl`,
  `P5_census_nostr.jl`; census `census/P5_T12`
- Logs: `$S/phase45/logs/P5_T12_{print_red,tr_red,print,pins,tr,golden_gen,golden,oracles,collapse,tme,tid,types,dsl,aqua_jet,docs,census,census_nostr,census_cmp,full,terter,tme_focused}.log`

**Interfaces:**
- Consumes: Task 1's `_solve`, `_Solve` (`names`, `dependents`, `fitted`), `_Dependent`
  (`col`, `keq`, `factors`, `label`); Task 2's `_Poly`, `_Mono`, `_Exp`, `_term`,
  `_monomial`, `_mmul`, `_mpow`, `_split_conc`, `_content`, `_nterms`; Tasks 3–4's `_law`,
  `_Law` (`vars`, `H`, `Q`, `D`, `trees`, `shift`), `_AlloLaw` (`vars`, `H`, `active`,
  `inactive`, `sites`, `n`), `_flip_only`; Task 5's `_assert_support`; Task 11's
  `_normalizer(D, vars, substrates)::_Norm` (`lambda`, `shift`) and
  `_term_key(mo, vars, substrates)`; `_check_well_posed`, `_derivation_input`,
  `metabolites`, `substrates`; tests: Task 3's `_testhelper_poly_value`,
  `_testhelper_evaluate_law`, Task 5's `_testhelper_law_value`, phase 4's
  `_testhelper_law_inputs`, `_testhelper_limit_law`, `_testhelper_single_subunit`,
  `_testhelper_lift`.
- Produces:
  ```julia
  _power_text(x, e::_Exp)::String              # x, x^2, sqrt(x), x^(3/2)
  _ratio_text(up::Vector{String}, dn::Vector{String})::String
  _term_text(c, mo::_Mono, names::Vector{Symbol}; rest::_Poly = one(_Poly),
             over::Bool = false)::String       # ratio of products
  _sum_text(p::_Poly, names::Vector{Symbol})::String
  _conc_label(mo::_Mono, names::Vector{Symbol})::String   # "0", "A²_B", "P_B⁻¹"
  struct _Term; conc::_Mono; c::Rational{Int}; mo::_Mono; rest::_Poly; over::Bool; end
  _printed(L::_Law, nz::_Norm, polys::Vector{_Poly}, subs::Vector{Int})::
      @NamedTuple{N::String, terms::Vector{Vector{_Term}}}
  _coefficient_sum(prefix::String, ts::Vector{_Term},
                   names::Vector{Symbol})::Tuple{String, Vector{String}}
  _v_line(q::AbstractString, h::AbstractString, d::AbstractString)::String
  _reference_lines(s::_Solve, vars::_Vars)::Vector{String}
  _header_lines(m, s::_Solve, vars::_Vars)::Vector{String}
  _law_text(m::Mechanism, s::_Solve, L::_Law; coefficients::Bool = false)::String
  _law_text(am::AllostericMechanism, s::_Solve, A::_AlloLaw;
            coefficients::Bool = false)::String
  rate_equation_string(m::Union{Mechanism, AllostericMechanism};
                       coefficients::Bool = false)::String    # runs _check_well_posed
  rate_equation_string(@nospecialize(m::AbstractEnzymeMechanism);
                       coefficients::Bool = false)::String    # runs _derivation_input
  _rate_text(m::Mechanism, coefficients::Bool)::String          # one method per
  _rate_text(am::AllostericMechanism, coefficients::Bool)::String   # mechanism type
  build_power_expr(keq_exp::Rational, factors)   # a fractional exponent: Rational{Int}
  # tests (fixtures file)
  _testhelper_eval_rate_lines(s::AbstractString, params, concs)::Dict{Symbol, BigFloat}
  _testhelper_eval_rate_string(s::AbstractString, params, concs)::BigFloat
  _testhelper_print_artifacts(s::AbstractString)::Vector{String}
  _testhelper_printed_law_error(m; points::Int = 2, seed::Int = 31)::BigFloat
  # tests (test_rate_eq_derivation.jl)
  _testhelper_extract_factors(v_line::AbstractString)::NTuple{3, String}   # (Q, H, D)
  _testhelper_reference_lines_match(m)::Bool
  _testhelper_sum_depth(t::AbstractString)::Int
  ```
  The full law reads, in order: the `params` line (`s.fitted`, `Keq`, `E_total`), the
  `concs` line (`metabolites(m)`), the Wegscheider and then the Haldane reference lines,
  `N = …` when the normalizer has a factor that is not a monomial, and
  `v = E_total * (Q) * (H) / (D)`. The MWC law prints, after the header, `N_c`, `Q_c`,
  `D_c` per conformation, `F_A`/`F_I` when they differ, `W_A`/`W_I`, then `v`.

**Gate:**
- test_print.jl's four new testsets pass before any pin moves: both forms of every spec
  equal the derived law to 1e-20 with no artifact (pin 1's value check), and the
  `:OnlyI` laws equal `_testhelper_limit_law` to 1e-20.
- `P5_T12_pins.jl` prints the texts quoted in Step 7, each with its value check passing.
- Q1 passes on the sample: both printed forms equal the law to 1e-20, with no artifact.
- The docs build runs the moved docstring's doctests and textbooks.md's.
- Pin 17's `Rational{Int}` check fails on the head (Step 3's TR run) and passes after
  Step 5 (Step 11).
- `test_aqua_jet.jl` passes in Step 11's focused run, before the full suite: JET finds no
  `_rate_text` call without a method (Interface decisions, "Mechanism and law together").
- `census/P5_T12` against `census/P5_base`: with the printed laws blanked, compare.jl
  `RESULT: PASS` (child sets, values ≤ 1e-9, fitted names), or a `FAIL` on exactly the
  summed-cap refusals of Task 6's `P5_T6_census_refused.log`; the printed text differs on
  every law both derive.
- R52's three figures beside Task 0's; the full suite green and pristine; about +8 s,
  most of it the 42-spec printed-law check (the commit reports the measured figure).

**Testsets added:**
- test_print.jl "the printed law and its coefficient form equal the law on every spec".
- test_print.jl "coefficient names read powers and ratios, and N prints the normalizer
  once", with the printed N of F4, Segel's ter bi and the power case, which Task 11 leaves
  to the printer.
- test_print.jl "an :OnlyI inhibitor prints a law equal to its limit" (Review Focus 1).
- test_print.jl "the MWC lines name each conformation's polynomials".

**Pins changed (each after its value check, RES §3.4):** 1, 2, 3, 4, 5 (value checks
added; the pinned lines hold), 6, 7, 8, 9, 17, R37's printed half, and test_oracles.jl's
half-integer dependent, whose exponent now prints `-1//2` (Step 8); 10 and 20 run
unchanged (Step 11). Pin 19 and TID's CSV-text check move in Task 12b.

- [ ] **Step 1: Add the printed-law oracles to the fixtures file**

In `test/mechanism_definitions_for_test_enzyme_derivation.jl`, at the end of the
"Phase-5 law oracles" section (Tasks 3, 4 and 11's helpers), directly before the line
`# ── Mechanism test specifications ───────────────────────────────────────────`, insert:

```julia
# ── Phase-5 printed-law oracles ─────────────────────────────────────────────

"""The exact value of an exponent a printed law writes: an integer or a ratio `(3/2)`."""
_testhelper_exponent(x::Integer) = x // 1
function _testhelper_exponent(ex::Expr)
    ex.head === :call && ex.args[1] === :/ && length(ex.args) == 3 ||
        error("unexpected exponent $ex in a printed law")
    _testhelper_exponent(ex.args[2]) // _testhelper_exponent(ex.args[3])
end

"""The value of the expression `ex` of a printed law at the names `env`, in the current
BigFloat precision: `+`, `-`, `*`, `/`, `sqrt` and `^` with an exact exponent."""
function _testhelper_printed_value(ex, env::AbstractDict{Symbol, BigFloat})
    ex isa Symbol && return env[ex]
    ex isa Number && return BigFloat(ex)
    ex isa Expr && ex.head === :call || error("unexpected $ex in a printed law")
    op, args = ex.args[1], ex.args[2:end]
    op === :^ && return _testhelper_printed_value(args[1], env)^
                        BigFloat(_testhelper_exponent(args[2]))
    a = [_testhelper_printed_value(x, env) for x in args]
    op === :+ ? sum(a) : op === :* ? prod(a) : op === :/ ? a[1] / a[2] :
    op === :- ? (length(a) == 1 ? -a[1] : a[1] - a[2]) : op === :sqrt ? sqrt(a[1]) :
    error("unexpected operator $op in a printed law")
end

"""The names an expression of a printed law reads, a called function's name aside."""
_testhelper_printed_names(ex::Symbol) = Set([ex])
_testhelper_printed_names(ex::Expr) =
    union(Set{Symbol}(), (_testhelper_printed_names(a)
                          for a in (ex.head === :call ? ex.args[2:end] : ex.args))...)
_testhelper_printed_names(ex) = Set{Symbol}()

"""
    _testhelper_eval_rate_lines(s, params, concs) -> Dict{Symbol, BigFloat}

Every name a printed law `s` (`rate_equation_string`, in full or in coefficients) binds,
evaluated in 256-bit BigFloat without `eval`, at `params` and `concs` (NamedTuples or
dictionaries by name), whose names are bound first. A destructuring line binds its names
from them; every other line runs once each name it reads is bound, so the coefficient form,
which prints `v` before its legend, runs too. Comment lines are skipped.
"""
function _testhelper_eval_rate_lines(s::AbstractString, params, concs)
    setprecision(BigFloat, 256) do
        given = Dict{Symbol, BigFloat}(k => BigFloat(x) for (k, x) in pairs(params))
        for (k, x) in pairs(concs)
            given[k] = BigFloat(x)
        end
        env = copy(given)
        todo = Expr[Meta.parse(l) for l in split(s, '\n')
                    if !isempty(strip(l)) && !startswith(lstrip(l), "#")]
        while !isempty(todo)
            n = length(todo)
            filter!(todo) do ex
                lhs, rhs = ex.args
                if lhs isa Expr
                    foreach(x -> env[x] = given[x], lhs.args[1].args)
                    return false
                end
                all(x -> haskey(env, x), _testhelper_printed_names(rhs)) || return true
                env[lhs] = _testhelper_printed_value(rhs, env)
                false
            end
            length(todo) == n && error("a printed law reads names it never binds: $todo")
        end
        env
    end
end

"""The rate `v` a printed law `s` gives at `params` and `concs`, in 256-bit BigFloat
(`_testhelper_eval_rate_lines`)."""
_testhelper_eval_rate_string(s::AbstractString, params, concs) =
    _testhelper_eval_rate_lines(s, params, concs)[:v]

const _testhelper_artifact_patterns = (r"\^1(?![\d/.)])", r"(?<![\w.)])1 \^",
                                       r"\* 1(?![\d./])", r"L \* 1(?!\d)", r"\(1\)",
                                       r"\^\(1/1\)")

"""The printing artifacts in the printed law `s`: an exponent 1 (`^1`, `^(1/1)`), 1 raised
to a power (`1 ^`), a factor 1 (`* 1`, `L * 1`) or a parenthesized 1."""
_testhelper_print_artifacts(s::AbstractString) =
    String[m.match for r in _testhelper_artifact_patterns for m in eachmatch(r, s)]

"""
The worst relative error, over `points` consistent points (`_testhelper_law_inputs`, seed
`seed`), of `m`'s printed law and its coefficient form (`rate_equation_string`), read back
in 256-bit BigFloat (`_testhelper_eval_rate_string`), against its derived law
(`EnzymeRates._law`). `m` is compiled or concrete and passes the guard; an allosteric
mechanism's point is drawn at one catalytic subunit, and the law is read at its own
multiplicity.
"""
function _testhelper_printed_law_error(m; points::Int = 2, seed::Int = 31)
    mech = _testhelper_lift(m)
    L = EnzymeRates._law(mech)
    fitted = collect(EnzymeRates._solve(mech).fitted)
    texts = (rate_equation_string(m), rate_equation_string(m; coefficients = true))
    rng = Random.MersenneTwister(seed)
    worst = big(0.0)
    for _ in 1:points
        _, concs, params = _testhelper_law_inputs(_testhelper_single_subunit(mech), fitted,
                                                  rng)
        law = _testhelper_evaluate_law(L, merge(params, concs))
        p = merge(params, Dict(:E_total => big(1)))
        for text in texts
            worst = max(worst, abs(_testhelper_eval_rate_string(text, p, concs) - law) /
                               abs(law))
        end
    end
    worst
end
```

In `test/test_rate_eq_derivation.jl`, delete the section that defined the old evaluator:

```julia
# ── Rate equation string evaluation helper ──────────────────────────────────

# Runs the rendered code as written: its destructuring lines read `params` and
# `concs`, its constraint lines define the dependent parameters, and its `v` line
# is the value.
function _testhelper_eval_rate_string(s, params, concs)
    code = "let params = $params, concs = $concs\n" *
           s * "\nend"
    eval(Meta.parse(code))
end
```

Its three callers in that file (the per-spec "Rate Equation String", "rate_equation_string
prints residual names that run as written" and "rate_equation_string runs a half-integer
dependent as written") read the fixtures file's evaluator unchanged: they compare its
BigFloat to `rate_equation(m, concs, params)` with `isapprox`.

- [ ] **Step 2: Write the failing printer tests and pin 17**

In `test/test_print.jl`, replace Task 11's two ABOUTME lines with

```julia
# ABOUTME: Tests of the unit rule (src/normalize.jl) and the printed laws (src/print.jl):
# ABOUTME: each rule's normalizer, both printed forms against the law, names, MWC lines.
```

and append at the end of the file:

```julia
@testset "the printed law and its coefficient form equal the law on every spec" begin
    # Both forms, read back in 256-bit BigFloat, equal the derived law to 1e-20 at a
    # consistent point of every spec, the concrete and the compiled mechanism print one
    # text, and neither form holds a printing artifact.
    for spec in MECHANISM_TEST_SPECS
        @testset "$(spec.name)" begin
            m = spec.mechanism
            full = rate_equation_string(m)
            @test rate_equation_string(_testhelper_lift(m)) == full
            @test isempty(_testhelper_print_artifacts(full))
            @test isempty(_testhelper_print_artifacts(
                rate_equation_string(m; coefficients = true)))
            @test _testhelper_printed_law_error(m; points = 1) <= big(10)^-20
        end
    end
end

@testset "coefficient names read powers and ratios, and N prints the normalizer once" begin
    # F1 reads A and P only through P/A, so the ratio form divides by A: the coefficients
    # of P/A take ⁻¹, and the polynomial factor of the divisor prints as N.
    f1 = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
        end
    end
    c1 = split(rate_equation_string(f1; coefficients = true), '\n')
    @test c1[1] == "v = E_total * (n_0 + n_P_A⁻¹ * P / A) * (A - P / Keq) / " *
                   "(1 + d_P_A⁻¹ * P / A)"
    @test c1[2] == "# The coefficients are functions of the fitted constants, not " *
                   "parameters:"
    @test c1[3:end] == [
        "N = k_E_A_to_F_P + k_F_A_to_E_P",
        "n_0 = 2 * k_E_A_to_F_P * k_F_A_to_E_P / N",
        "n_P_A⁻¹ = 2 * k_E_A_to_F_P * k_F_A_to_E_P / (Keq * N)",
        "d_P_A⁻¹ = (k_E_A_to_F_P * k_F_A_to_E_P + k_F_P_to_E_A^2 * Keq^2) / " *
        "(k_F_P_to_E_A * Keq^2 * N)"]
    # F3's rapid-equilibrium Theorell–Chance step makes every denominator term carry B or
    # P: the ratio form divides by B, whose coefficient is one constant, so no N prints.
    f3 = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B ⇌ E(Q) + P
            E(Q) <--> E + Q
        end
    end
    c3 = split(rate_equation_string(f3; coefficients = true), '\n')
    @test c3[1] == "v = E_total * (n_B⁻¹ / B) * (A * B - P * Q / Keq) / (1 + " *
                   "d_P_B⁻¹ * P / B + d_A * A + d_A_P_B⁻¹ * A * P / B + d_Q * Q + " *
                   "d_P_Q_B⁻¹ * P * Q / B)"
    @test "n_B⁻¹ = k_E_A_to_EA" in c3
    @test !any(startswith("N = "), c3)
    # F2's weights are fractional: N is a power 3/2, and its coefficients hold roots.
    f2 = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
            F + A + A ⇌ F(A, A)
        end
    end
    @test occursin("\nN = (k_E_A_to_F_P + k_F_A_to_E_P)^(3/2)\n", rate_equation_string(f2))
    @test "n_A = 2 * sqrt(K_F_A_A_to_FAA) * k_E_A_to_F_P^(3/2) * k_F_A_to_E_P / N" in
          split(rate_equation_string(f2; coefficients = true), '\n')
    # The random-order bi-bi's denominator holds A^2 * B, named by its superscript.
    ro = only(s for s in MECHANISM_TEST_SPECS if s.name == "Random-order Bi-Bi").mechanism
    @test occursin(" + d_A²_B * A^2 * B + ", rate_equation_string(ro; coefficients = true))
    # The unit-rule testsets above pin the weights of F4, Segel's ter bi and the ping-pong
    # power case (UNIT §2); N prints the normalizer's polynomial factor alone, its integer
    # and monomial content moved into the coefficients. F4's is a sum squared. Segel's
    # graded N is the product of two sums, where Cleland's would divide it by a third sum
    # squared (P5-10a). The power case's is one sum squared, Cleland's N squared.
    f4 = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
            F + A ⇌ F(A)
            F(A) + A ⇌ F(A, A)
        end
    end
    @test occursin("\nN = (k_E_A_to_F_P + k_F_A_to_E_P)^2\n", rate_equation_string(f4))
    terbi = only(s for s in MECHANISM_TEST_SPECS
                 if s.name == "Segel Bi Uni Uni Uni PP Ter Bi").mechanism
    power = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> E(P; residual = A - P), E(Q) + A <--> E(A, Q))
            (E + Q <--> E(Q), E(P; residual = A - P) + Q <--> E(A, Q))
            (E(; residual = A - P) + B <--> E(Q),
             E(P; residual = A - P) + B <--> E(B, P; residual = A - P))
            (E(; residual = A - P) + P ⇌ E(P; residual = A - P),
             E(Q) + P ⇌ E(B, P; residual = A - P))
        end
    end
    n_line(m) = only(filter(startswith("N = "), split(rate_equation_string(m), '\n')))
    @test occursin(r"^N = \([^()]+\) \* \([^()]+\)$", n_line(terbi))
    @test occursin(r"^N = \([^()]+\)\^2$", n_line(power))
    for m in (f1, f2, f3, f4, terbi, power)
        @test _testhelper_printed_law_error(m) <= big(10)^-20
        @test isempty(_testhelper_print_artifacts(rate_equation_string(m))) &&
              isempty(_testhelper_print_artifacts(
                  rate_equation_string(m; coefficients = true)))
    end
end

@testset "an :OnlyI inhibitor prints a law equal to its limit" begin
    # Review Focus 1. I binds only the inactive conformation, at rapid equilibrium and at
    # steady state. Both printed forms equal the rapid-equilibrium limit of the oracle's
    # own network (`_testhelper_limit_law`) to 1e-20 at two points; the active
    # conformation's denominator holds no constant of I's binding, the inactive one's does.
    competitive_RE = @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)            :: EqualAI
            E(S) <--> E + P         :: EqualAI
            E + I ⇌ E(I)            :: OnlyI
        end
    end
    competitive_SS = @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)            :: EqualAI
            E(S) <--> E + P         :: NonequalAI
            E + I <--> E(I)         :: OnlyI
        end
    end
    cases = ((competitive_RE, "D_I = 1 + S / K_ES_to_E_S + I / K_I_EI_to_E_I",
              "dI_I = 1 / K_I_EI_to_E_I"),
             (competitive_SS,
              "D_I = 1 + S / K_ES_to_E_S + I * k_I_E_I_to_EI / k_I_EI_to_E_I",
              "dI_I = k_I_E_I_to_EI / k_I_EI_to_E_I"))
    for (i, (m, d_i, d_i_coefficient)) in enumerate(cases)
        full = rate_equation_string(m)
        coef = rate_equation_string(m; coefficients = true)
        @test "D_A = 1 + S / K_ES_to_E_S" in split(full, '\n')
        @test d_i in split(full, '\n')
        @test last(split(full, '\n')) ==
              "v = E_total * (Q_A + L * Q_I) * (S - P / Keq) / (D_A + L * D_I)"
        @test "D_A = 1 + dA_S * S" in split(coef, '\n')
        @test d_i_coefficient in split(coef, '\n')
        am = _testhelper_lift(m)
        cm = EnzymeRates.catalytic(am)
        onlyi = [String(x) for (k, ns) in enumerate(EnzymeRates._edge_names(am, :I))
                 if EnzymeRates.cat_allo_state(am, cm.group[k]) === :OnlyI for x in ns]
        d_a = only(l for l in split(full, '\n') if startswith(l, "D_A = "))
        @test !isempty(onlyi)
        @test !any(x -> occursin(x, d_a), onlyi)
        @test any(x -> occursin(x, d_i), onlyi)
        s = EnzymeRates._solve(am)
        rng = Random.MersenneTwister(500 + i)
        for _ in 1:2
            point, concs, params = _testhelper_law_inputs(am, s.fitted, rng)
            exact = setprecision(BigFloat, 256) do
                _testhelper_limit_law(am, point, concs, Set{Symbol}(); rng)
            end
            p = merge(params, Dict(:E_total => big(1)))
            for text in (full, coef)
                @test abs(_testhelper_eval_rate_string(text, p, concs) - exact) <=
                      big(10)^-20 * abs(exact)
            end
        end
        @test _testhelper_printed_law_error(m) <= big(10)^-20
    end
end

@testset "the MWC lines name each conformation's polynomials" begin
    # A two-site homodimer prints Q, D and W for each conformation and no F: both
    # free-enzyme weights are 1. test_rate_eq_derivation.jl pins its text.
    two_site = @allosteric_mechanism begin
        substrates: S
        products:   P
        allosteric_regulators: R
        catalytic_multiplicity: 2
        steps: begin
            E + P ⇌ E(P)     :: NonequalAI
            E + S ⇌ E(S)     :: NonequalAI
            E(S) <--> E(P)   :: NonequalAI
        end
        regulatory_site(2): R::NonequalAI
    end
    # A ping-pong whose second half runs in the active conformation only: neither
    # conformation's denominator has a constant term, so each prints its own normalizer,
    # and free E weighs F6P in one and F16BP in the other, so F_A and F_I print.
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
    # Both bindings :OnlyA: the inactive conformation is free E alone, so the
    # denominator's inactive term is L with no factor.
    inert = @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: EqualAI
            E + P ⇌ E(P)      :: OnlyA
        end
    end
    heads(t) = [l[1:3] for l in split(t, '\n') if occursin(r"^[QDFWN]_[AI] = ", l)]
    @test heads(rate_equation_string(two_site)) ==
          ["Q_A", "D_A", "Q_I", "D_I", "W_A", "W_I"]
    t2 = rate_equation_string(half2)
    @test heads(t2) == ["N_A", "Q_A", "D_A", "N_I", "D_I", "F_A", "F_I"]
    @test "F_A = F6P / (K_A_EF6P_res_plusATP_minusF16BP_to_E_res_plusATP_minusF16BP_F6P " *
          "* N_A)" in split(t2, '\n')
    @test last(split(t2, '\n')) == "v = E_total * (F_I * Q_A) * (ATP * F6P - ADP * " *
                                    "F16BP / Keq) / (F_I * D_A + L * F_A * D_I)"
    ti = rate_equation_string(inert)
    @test heads(ti) == ["Q_A", "D_A"]
    @test last(split(ti, '\n')) == "v = E_total * (Q_A * D_A) * (S - P / Keq) / (D_A^2 + L)"
    for m in (two_site, half2, inert)
        @test _testhelper_printed_law_error(m) <= big(10)^-20
        @test isempty(_testhelper_print_artifacts(rate_equation_string(m))) &&
              isempty(_testhelper_print_artifacts(
                  rate_equation_string(m; coefficients = true)))
    end
end
```

In `test/test_rate_eq_derivation.jl` (pin 17, written before Step 5 changes
`build_power_expr`), in the testset "build_power_expr return types", after
`@test bpe(R(1), [(:k1f, R(1))]) isa Union{Float64, Symbol, Expr}` add:

```julia

        # A fractional exponent is an exact Rational{Int} literal, so the expression
        # evaluates exactly in BigFloat.
        @test bpe(R(-1//2), [(:k1f, R(1//2))]) ==
              Expr(:call, :*, Expr(:call, :^, :Keq, -1//2), Expr(:call, :^, :k1f, 1//2))
        @test bpe(R(-1//2), [(:k1f, R(1//2))]).args[2].args[3] isa Rational{Int}
        half = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                (E(S) <--> F(S), E <--> F)
                F(S) <--> E + P
                F + S <--> E(S)
            end
        end
        dep, fitted = EnzymeRates._dependent_param_exprs(typeof(half))
        setprecision(BigFloat, 256) do
            p = NamedTuple{(fitted..., :Keq)}(
                Tuple(big(1) + big(i) / 7 for i in 1:(length(fitted) + 1)))
            env = Dict{Symbol, BigFloat}(k => v for (k, v) in pairs(p))
            @test abs(_testhelper_eval_expr(dep[:k_F_to_E], env) -
                      all_parameters(half, p).k_F_to_E) <=
                  big(10)^-60 * all_parameters(half, p).k_F_to_E
        end
```

- [ ] **Step 3: Run the new tests to verify they fail**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T12_print_red.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); @testset "red" begin include("test/test_print.jl") end'`
Expected: Task 11's unit-rule testsets pass; inside the outer `red` testset each of the
four new testsets errors with
`MethodError: no method matching rate_equation_string(…; coefficients::Bool)` (the head's
method takes no keyword), so every spec's subtestset errors; `exit 1`.

Then run (background, polled):
`S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T12_tr_red.log 3600 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); @testset "red" begin include("test/test_rate_eq_derivation.jl") end'`
Expected: "build_power_expr return types" fails at the new `isa Rational{Int}` assertion
alone (the head emits the Float64 `-0.5`; `==` compares the literals by value, and 0.5 is
exact in binary, so the `==` and the BigFloat check pass); "rate_equation_string runs a
half-integer dependent as written" errors in `_testhelper_exponent`, since Step 1's
evaluator reads exponents exactly and has no method for the head's printed `^ -0.5`
(Step 4's printer writes `sqrt`); every other testset passes, so Step 1's evaluator reads
the head's printed laws; `exit 1`.

- [ ] **Step 4: Write the printer**

Create `src/print.jl`:

````julia
# ABOUTME: Printed rate laws: the law in the independent constants, with reference lines and
# ABOUTME: the normalizer's factor, its coefficient form, and MWC laws per conformation.

# ─── Text of powers, terms and sums ─────────────────────────────────────────

"""`x^e` as text for a positive exponent `e`: `x`, `x^2`, `sqrt(x)`, `x^(3/2)`."""
_power_text(x, e::_Exp) = e == 1 ? string(x) : e == 1 // 2 ? "sqrt($x)" :
    isinteger(e) ? "$x^$(numerator(e))" : "$x^($(numerator(e))/$(denominator(e)))"

"""Factors as a ratio, `a * b / (c * d)`; `1` when no factor is above the bar."""
_ratio_text(up::Vector{String}, dn::Vector{String}) =
    string(isempty(up) ? "1" : join(up, " * "),
           isempty(dn) ? "" :
           " / " * (length(dn) == 1 ? only(dn) : "(" * join(dn, " * ") * ")"))

"""
    _term_text(c, mo, names; rest = one(_Poly), over = false) -> String

The positive term `c·mo` as a ratio of products over the variable names `names`, times the
polynomial `rest` and over the normalizer's printed factor `N` when `over`:
`2 * A * k1 * (k2 + k3) / (k4 * N)`, and `1` for the empty term.
"""
function _term_text(c, mo::_Mono, names::Vector{Symbol}; rest::_Poly = one(_Poly),
                    over::Bool = false)
    up = String[_power_text(names[v], e) for (v, e) in mo if e > 0]
    dn = String[_power_text(names[v], -e) for (v, e) in mo if e < 0]
    numerator(c) == 1 || pushfirst!(up, string(numerator(c)))
    denominator(c) == 1 || pushfirst!(dn, string(denominator(c)))
    rest == one(_Poly) || push!(up, "(" * _sum_text(rest, names) * ")")
    over && push!(dn, "N")
    _ratio_text(up, dn)
end

"""The nonzero polynomial `p` over the variable names `names` as a sum, `a + 2 * b - c`,
its terms in monomial order."""
function _sum_text(p::_Poly, names::Vector{Symbol})
    s = join(string(c < 0 ? " - " : " + ", _term_text(abs(c), mo, names))
             for (mo, c) in sort!(collect(p.terms); by = first))
    (startswith(s, " + ") ? "" : "-") * chop(s; head = 3, tail = 0)
end

const _SUPERSCRIPTS = Dict(zip("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻"))

"""
A coefficient's name for the concentration monomial `mo` over the names `names`: `0` for
none, else each concentration with its exponent as a superscript, positive exponents first
and each part in variable order, joined by `_`: `A²_B`, `P_B⁻¹`.
"""
function _conc_label(mo::_Mono, names::Vector{Symbol})
    isempty(mo) && return "0"
    sup(e) = e == 1 ? "" : String([_SUPERSCRIPTS[d] for d in string(numerator(e))])
    join((string(names[v], sup(e))
          for (v, e) in sort(mo; by = p -> (p.second < 0, p.first))), "_")
end

# ─── A law's polynomials over its normalizer ─────────────────────────────────

"""
One printed term of a law's polynomial: the concentration monomial `conc` times its
coefficient `c·mo·rest`, over the normalizer's printed factor `N` when `over`. `c` and `mo`
have the normalizer's integer and monomial contents divided out; `rest` is the
coefficient's polynomial factor, `one(_Poly)` when there is none.
"""
struct _Term
    conc::_Mono
    c::Rational{Int}
    mo::_Mono
    rest::_Poly
    over::Bool
end

"""Multiply the factors printed apart, `apart` (polynomial and exponent), by `p^e`: an
equal polynomial adds its exponent, and a factor whose exponent sums to 0 goes."""
function _apart!(apart::Vector{Tuple{_Poly, _Exp}}, p::_Poly, e::_Exp)
    i = findfirst(x -> x[1] == p, apart)
    i === nothing ? push!(apart, (p, e)) : (apart[i] = (p, apart[i][2] + e))
    filter!(x -> !iszero(x[2]), apart)
end

"""
    _printed(L, nz, polys, subs) -> (; N, terms)

The polynomials `polys` of the conformation law `L` divided by its normalizer `nz`
(`_normalizer`, `subs` the substrates' variables). The normalizer's integer and monomial
contents divide every coefficient, and its other factors print apart: `N` is their text,
empty when there is none. `terms[j]` holds a `_Term` for each concentration monomial of
`polys[j]`, in the unit rule's order (`_term_key`), divided by `nz.shift`. A coefficient
whose polynomial factor is `N`'s one factor drops it and is not over `N`.
"""
function _printed(L::_Law, nz::_Norm, polys::Vector{_Poly}, subs::Vector{Int})
    D, names = _split_conc(L.D, L.vars), L.vars.names
    k, mu, apart = 1 // 1, _Mono(), Tuple{_Poly, _Exp}[]
    for (mo, e) in nz.lambda
        c, u, rest = _content(D[mo])
        if isinteger(e)
            k *= (c // 1)^numerator(e)
        elseif c != 1
            _apart!(apart, _term(_Mono(), c), e)
        end
        mu = _mmul(mu, _mpow(u, e))
        rest == one(_Poly) || _apart!(apart, rest, e)
    end
    paren(p) = "(" * _sum_text(p, names) * ")"
    N = isempty(apart) ? "" :
        length(apart) == 1 && apart[1][2] == 1 ? _sum_text(apart[1][1], names) :
        _ratio_text(String[_power_text(paren(p), e) for (p, e) in apart if e > 0],
                    String[_power_text(paren(p), -e) for (p, e) in apart if e < 0])
    whole = length(apart) == 1 && apart[1][2] == 1 ? apart[1][1] : zero(_Poly)
    # `term` captures only names assigned once, so none is boxed (the loop reassigns `k`).
    scale, inv_mu, inv_shift = k, _mpow(mu, -1), _mpow(nz.shift, -1)
    function term(conc, q)
        c, u, rest = _content(q)
        w = rest == whole
        _Term(_mmul(conc, inv_shift), (c // 1) / scale, _mmul(u, inv_mu),
              w ? one(_Poly) : rest, !isempty(apart) && !w)
    end
    terms = Vector{_Term}[
        _Term[term(conc, q) for (conc, q) in sort!(collect(_split_conc(p, L.vars));
                                                   by = x -> _term_key(x[1], L.vars, subs))]
        for p in polys]
    (; N, terms)
end

"""The terms `ts` (`_printed`) as a sum, each concentration monomial times its
coefficient."""
_group_sum(ts::Vector{_Term}, names::Vector{Symbol}) =
    join((_term_text(t.c, _mmul(t.conc, t.mo), names; rest = t.rest, over = t.over)
          for t in ts), " + ")

"""`label` times the concentration monomial `conc`: `d_A * A`, `d_P_B⁻¹ * P / B`,
`n_B⁻¹ / B`."""
function _label_times(label::String, conc::_Mono, names::Vector{Symbol})
    t = _term_text(1, conc, names)
    startswith(t, "1 / ") ? label * " / " * t[5:end] : label * " * " * t
end

"""
    _coefficient_sum(prefix, ts, names) -> (sum, legend)

The terms `ts` (`_printed`) as a sum of named coefficients, `prefix_A_B * A * B`, with
`prefix_0` alone for the constant monomial and `1` when its coefficient is 1, and one
legend line per named coefficient giving its value in the fitted constants.
"""
function _coefficient_sum(prefix::String, ts::Vector{_Term}, names::Vector{Symbol})
    unit(t) = isempty(t.conc) && t.c == 1 && isempty(t.mo) && t.rest == one(_Poly) &&
              !t.over
    label(t) = string(prefix, "_", _conc_label(t.conc, names))
    sum = join((unit(t) ? "1" : isempty(t.conc) ? label(t) :
                _label_times(label(t), t.conc, names) for t in ts), " + ")
    legend = String[string(label(t), " = ",
                           _term_text(t.c, t.mo, names; rest = t.rest, over = t.over))
                    for t in ts if !unit(t)]
    sum, legend
end

"""The `v` line, `v = E_total * (Q) * (H) / (D)`: `Q` bare when it is one name and left out
when it is 1."""
function _v_line(q::AbstractString, h::AbstractString, d::AbstractString)
    lead = q == "1" ? "" : occursin(r"^[^ ()*/+^-]+$", q) ? q * " * " : "(" * q * ") * "
    "v = E_total * $(lead)($h) / ($d)"
end

const _COEFFICIENTS_COMMENT =
    "# The coefficients are functions of the fitted constants, not parameters:"

# ─── Reference lines and header ──────────────────────────────────────────────

"""
The reference lines of the solve `s` over the variables `vars`: a `# Wegscheider
constraints:` section, then a `# Haldane constraints:` section, each dependent constant by
name as a ratio of products of `Keq` and the independents with exact exponents (`sqrt`,
`^(3/2)`), `= 1` for a constant that thermodynamics pins to 1. A section without a line is
left out; a kinetic group's other members equal its kept constant and print no line.
"""
function _reference_lines(s::_Solve, vars::_Vars)
    lines = String[]
    for (label, title) in ((:wegscheider, "Wegscheider"), (:haldane, "Haldane"))
        ds = sort!(filter(d -> d.label === label, s.dependents); by = d -> s.names[d.col])
        isempty(ds) || push!(lines, "# $title constraints:")
        for d in ds
            mo = _monomial((vars.index[s.names[f]] => e for (f, e) in d.factors)...,
                           vars.index[:Keq] => d.keq)
            push!(lines, "$(s.names[d.col]) = $(_term_text(1, mo, vars.names))")
        end
    end
    lines
end

"""The `params` and `concs` destructuring lines of `m`'s printed law, then the reference
lines (`_reference_lines`)."""
_header_lines(m::Union{Mechanism, AllostericMechanism}, s::_Solve, vars::_Vars) =
    String["(; $(join((s.fitted..., :Keq, :E_total), ", "))) = params";
           "(; $(join(metabolites(m), ", "))) = concs";
           _reference_lines(s, vars)]

# ─── Plain and MWC laws ──────────────────────────────────────────────────────

"""
    _law_text(m, s, L; coefficients = false) -> String

The printed law `L` of `m` over the solve `s` (`rate_equation_string`). In full: the header
(`_header_lines`), `N = …` when the normalizer has a factor that is not a monomial, and the
`v` line (`_v_line`) with `Q` and `D` over the normalizer. With `coefficients`, the `v`
line in named coefficients, then `N` and one legend line per coefficient.
"""
function _law_text(m::Mechanism, s::_Solve, L::_Law; coefficients::Bool = false)
    subs = [L.vars.index[x] for x in substrates(m)]
    P = _printed(L, _normalizer(L.D, L.vars, subs), [L.Q, L.D], subs)
    names, (q, d) = L.vars.names, P.terms
    h = _sum_text(L.H, names)
    if coefficients
        (qs, lq) = _coefficient_sum("n", q, names)
        (ds, ld) = _coefficient_sum("d", d, names)
        lines = [_v_line(qs, h, ds), _COEFFICIENTS_COMMENT]
        isempty(P.N) || push!(lines, "N = $(P.N)")
        return join([lines; lq; ld], "\n")
    end
    lines = _header_lines(m, s, L.vars)
    isempty(P.N) || push!(lines, "N = $(P.N)")
    push!(lines, _v_line(_group_sum(q, names), h, _group_sum(d, names)))
    join(lines, "\n")
end

"""
The conformation law `L` of an MWC law, tagged `tag` (`"A"` or `"I"`), printed over its own
normalizer: the terms (`_printed`) of its `Q`, its `D` and its free-E weight (the tree sum
of free E's segment times the law's shift), each over the normalizer's printed factor,
which reads as the name `N_<tag>` appended to the variable names.
"""
function _conformation(L::_Law, tag::String, subs::Vector{Int})
    names = [L.vars.names; Symbol("N_", tag)]
    n_var = _monomial(length(names) => -1)
    over(ts) = _Term[_Term(t.conc, t.c, t.over ? _mmul(t.mo, n_var) : t.mo, t.rest, false)
                     for t in ts]
    P = _printed(L, _normalizer(L.D, L.vars, subs), [L.Q, L.D, L.trees[1] * L.shift], subs)
    (; tag, names, N = P.N, q = over(P.terms[1]), d = over(P.terms[2]),
     f = over(P.terms[3]))
end

"""Push `W_<tag> = (1 + …)^k * …`, the product of the regulatory sites' binding
polynomials `ws` (polynomial and multiplicity) in one conformation, onto `defs`; a site
whose polynomial is 1 is left out, and so is the line when every site is."""
function _site_line!(defs::Vector{String}, tag::String, ws::Vector{Tuple{_Poly, Int}},
                     names::Vector{Symbol})
    parts = [string("(", _sum_text(w, names), ")", k == 1 ? "" : "^$k")
             for (w, k) in ws if _nterms(w) > 1]
    isempty(parts) || push!(defs, "W_$tag = " * join(parts, " * "))
    defs
end

"""
The printed MWC law `A` of `am` over the solve `s`. Each conformation `c` (`A`, then `I`
unless it is free E alone) defines `Q_c` (unless it turns nothing over) and `D_c`, over its
own normalizer `N_c`, which prints as a line when it has a factor that is not a monomial;
`F_A` and `F_I`, the free-E weights over the normalizers, print only when they differ and
each only when it is not 1; `W_A` and `W_I` multiply the regulatory sites' factors. The `v`
line is `E_total·H·(F_I^n·Q_A·D_A^(n−1)·W_A + L·F_A^n·Q_I·D_I^(n−1)·W_I) /
(F_I^n·D_A^n·W_A + L·F_A^n·D_I^n·W_I)` without its factors of 1. In full, the header comes
first and the definitions before `v`; with `coefficients`, `v` comes first, the sums read
in coefficients named by conformation (`nA_0`, `dA_S`, `fI_B`), and the normalizers and the
legend follow.
"""
function _law_text(am::AllostericMechanism, s::_Solve, A::_AlloLaw;
                   coefficients::Bool = false)
    subs = [A.vars.index[x] for x in substrates(am)]
    confs = [_conformation(A.active, "A", subs)]
    _flip_only(A.inactive) || push!(confs, _conformation(A.inactive, "I", subs))
    defs, legend = String[], String[]
    function define!(name, ts, c)
        if coefficients
            prefix = string(name[1] == 'Q' ? 'n' : lowercase(name[1]), c.tag)
            text, lines = _coefficient_sum(prefix, ts, c.names)
            push!(defs, "$name = $text")
            append!(legend, lines)
        else
            push!(defs, "$name = $(_group_sum(ts, c.names))")
        end
    end
    for c in confs
        isempty(c.N) || push!(coefficients ? legend : defs, "N_$(c.tag) = $(c.N)")
        isempty(c.q) || define!("Q_$(c.tag)", c.q, c)
        define!("D_$(c.tag)", c.d, c)
    end
    flips = [_group_sum(c.f, c.names) for c in confs]
    if flips[1] != (length(flips) == 2 ? flips[2] : "1")
        for (c, t) in zip(confs, flips)
            t == "1" || define!("F_$(c.tag)", c.f, c)
        end
    end
    _site_line!(defs, "A", Tuple{_Poly, Int}[(w, k) for (w, _, k) in A.sites], A.vars.names)
    _site_line!(defs, "I", Tuple{_Poly, Int}[(w, k) for (_, w, k) in A.sites], A.vars.names)
    has(x) = any(l -> startswith(l, x * " = "), defs)
    pw(x, e) = !has(x) || e == 0 ? "" : e == 1 ? x : "$x^$e"
    product(fs...) = (f = filter(!isempty, String[fs...]);
                      isempty(f) ? "1" : join(f, " * "))
    branch(c, o, e) = product(pw("F_$o", A.n), pw("Q_$c", e), pw("D_$c", A.n - e),
                              pw("W_$c", 1))
    num = String[]
    has("Q_A") && push!(num, branch("A", "I", 1))
    has("Q_I") && push!(num, "L * " * branch("I", "A", 1))
    rest = branch("I", "A", 0)
    den = branch("A", "I", 0) * " + L" * (rest == "1" ? "" : " * " * rest)
    v = _v_line(isempty(num) ? "0" : join(num, " + "), _sum_text(A.H, A.vars.names), den)
    coefficients ? join([v; defs; _COEFFICIENTS_COMMENT; legend], "\n") :
                   join([_header_lines(am, s, A.vars); defs; v], "\n")
end

# ─── rate_equation_string ────────────────────────────────────────────────────

"""
    rate_equation_string(m; coefficients = false) -> String

Return the rate law of mechanism `m` as a multi-line `String` (it returns the text; it does
not print). Pass a concrete `Mechanism` / `AllostericMechanism` or its compiled
[`EnzymeMechanism`](@ref) singleton.

The text runs as Julia code once `params` and `concs` are bound. It destructures
[`parameters`](@ref)`(m)` from `params` and the metabolites from `concs`, lists each
dependent constant under `# Wegscheider constraints:` or `# Haldane constraints:` as a
reference line in `Keq` and the independent constants, and ends with
`v = E_total * (Q) * (H) / (D)` in the independent constants. `H = Π[S] − Π[P]/Keq`
vanishes at equilibrium; `Q` and `D` have positive terms. Numerator and denominator are
divided by a normalizer that makes every term of `D` dimensionless: `D`'s constant term
when it has one, so `D` reads `1 + S/K…`; otherwise a product of powers of `D`'s
coefficients; and when no such product has `D`'s unit, `D`'s first term, so terms such as
`P / B` appear. A factor of the normalizer that is not a monomial prints once as `N = …`.
An allosteric mechanism prints each conformation's numerator `Q_A`, `Q_I`, denominator
`D_A`, `D_I`, normalizer `N_A`, `N_I`, free-enzyme weight `F_A`, `F_I` (when they differ)
and regulatory-site factor `W_A`, `W_I`, then the MWC combination.

With `coefficients = true` the law reads in named coefficients,
`v = E_total * n_0 * (H) / (1 + d_S * S + …)`, each named by the concentrations it
multiplies (`d_A²_B` for `A^2 * B`, and `d_P_B⁻¹` for `P / B` in a law that reads them only
through their ratio), followed by one line per coefficient giving its value in the fitted
constants.

```jldoctest
julia> using EnzymeRates

julia> m = @enzyme_mechanism begin
           substrates: S
           products: P
           steps: begin
               E + S ⇌ E(S)
               E(S) <--> E(P)
               E(P) ⇌ E + P
           end
       end;

julia> print(rate_equation_string(m))
(; K_EP_to_E_P, K_ES_to_E_S, k_ES_to_EP, Keq, E_total) = params
(; S, P) = concs
# Haldane constraints:
k_EP_to_ES = K_EP_to_E_P * k_ES_to_EP / (K_ES_to_E_S * Keq)
v = E_total * (k_ES_to_EP / K_ES_to_E_S) * (S - P / Keq) / (1 + S / K_ES_to_E_S + P / K_EP_to_E_P)

julia> print(rate_equation_string(m; coefficients = true))
v = E_total * n_0 * (S - P / Keq) / (1 + d_S * S + d_P * P)
# The coefficients are functions of the fitted constants, not parameters:
n_0 = k_ES_to_EP / K_ES_to_E_S
d_S = 1 / K_ES_to_E_S
d_P = 1 / K_EP_to_E_P
```
"""
function rate_equation_string(m::Union{Mechanism, AllostericMechanism};
                              coefficients::Bool = false)
    _check_well_posed(m)
    _rate_text(m, coefficients)
end

rate_equation_string(@nospecialize(m::AbstractEnzymeMechanism);
                     coefficients::Bool = false) =
    _rate_text(_derivation_input(typeof(m)), coefficients)

"""The printed law of the well-posed mechanism `m` (`_law_text`), once class 5 passes on
its derived law (`_assert_support`). One method per mechanism type, so the law's type
follows the mechanism's and `_law_text` dispatches on a matching pair."""
function _rate_text(m::Mechanism, coefficients::Bool)
    s = _solve(m)
    L = _law(m, s)
    _assert_support(m, L)
    _law_text(m, s, L; coefficients)
end

function _rate_text(am::AllostericMechanism, coefficients::Bool)
    s = _solve(am)
    A = _law(am, s)
    _assert_support(am, A)
    _law_text(am, s, A; coefficients)
end
````

Notes:
- R18: `_printed` and `_conformation` return one concrete `NamedTuple` type each, and
  `_Term` holds three references and two bits fields. `_rate_text` has one method per
  mechanism type, and `_law` dispatches on that type, so `L` is concrete in each method
  and the call to `_law_text` matches a method. One method on the `Union` would leave `L`
  a `Union` unlinked to `m`, and JET would report the two mismatched pairs (Interface
  decisions, "Mechanism and law together").
- The closures (`term`, `paren`, `define!`, `has`, `pw`, `product`, `branch`, `over`)
  capture only locals assigned once, so nothing is boxed: `term` reads `scale`, not the
  loop's `k`. Every comprehension that can be empty (a conformation that turns nothing
  over, a mechanism without sites) names its element type, so Julia 1.10 builds the same
  vectors (P5-38).
- The parameter-naming chokepoint test sees no `Symbol("K…")`, `Symbol("k…")`,
  `Symbol("V…")` or `Symbol("L…")` literal: the one `Symbol` call builds `N_A` or `N_I`.

In `src/EnzymeRates.jl`, add `include("print.jl")` directly after `include("normalize.jl")`.

In `src/rate_eq_derivation.jl`, delete from the line `"""` above
`    rate_equation_string(m) -> String` through the two lines

```julia
rate_equation_string(m::Union{Mechanism, AllostericMechanism}) =
    rate_equation_string(compile_mechanism(m))
```

(the old docstring, `function rate_equation_string end`, the comment line above the
concrete method and the method), and delete the compiled method below `_equation_text`:

```julia
function rate_equation_string(@nospecialize(m::AbstractEnzymeMechanism))
    mech = _derivation_input(typeof(m))
    dep, indep = _dependent_param_exprs(mech)
    _equation_text(mech, (indep..., :Keq, :E_total), dep)
end
```

`_equation_text`, its docstring and `_assert_no_name_clash` stay one task with no caller:
Task 13 deletes the old engine and printer as one unit, under one grep over every old
name, and deleting the printer's entry here would orphan `_expr_to_string` and its string
helpers one at a time. The commit says so (CLAUDE.md's rule against unused code; listed
for Denis, item 43). Run `grep -n "rate_equation_string" src/rate_eq_derivation.jl`: only
the
`all_parameters` docstring's `[`rate_equation_string`](@ref)` and `_equation_text`'s
docstring remain.

- [ ] **Step 5: Exact exponents in `build_power_expr`**

In `src/sym_poly_for_rate_eq_derivation.jl`, in `build_power_expr`, replace

```julia
        !isinteger(exp) ? :($sym ^ $(Float64(exp))) :
```

with

```julia
        !isinteger(exp) ? :($sym ^ $(Rational{Int}(exp))) :
```

and in its docstring replace "The empty product is the\nFloat64 `1.0`: a constant that
thermodynamics pins to 1." with "A fractional exponent is an\nexact `Rational{Int}` literal.
The empty product is the Float64 `1.0`: a constant that\nthermodynamics pins to 1." The
renderer `_dependent_param_exprs` is the only caller; no src code evaluates the `Expr`.
Task 1's test_joint_solve.jl testset "_solve renders a half-integer dependent with keq
-1//2" compares its `Expr` with `==`, which reads number literals by value (`-0.5 ==
-1//2`), so it passes unchanged; the fixtures file's `_testhelper_eval_expr` reads a
`Rational` literal as `BigFloat(x)`. One test reads the `Expr`'s text: test_oracles.jl's
"a dependent with a half-integer exponent matches the oracle" asserts
`occursin("^ -0.5", string(ex))`, and the `Rational{Int}` literal prints `Keq ^ -1//2`
(`$S/phase45/plan5/probes/review_exec/rat_expr.log`). Step 8 re-pins it after Step 7's
value check.

- [ ] **Step 6: Run the printer tests to verify they pass**

Rerun Step 3's first command (test_print.jl) with log `P5_T12_print.log`; TR, pin 17
among it, runs green in Step 11, once Step 8 has moved its printed-law pins.
Expected: every testset passes, pristine, `exit 0`; the 42 spec subtestsets of "the
printed law and its coefficient form equal the law on every spec" pass (pin 1's and pin
2's value checks: every spec's printed laws, the 11 allosteric ones among them, equal the
derived law to 1e-20, which Tasks 3–4 tie to the oracle at 1e-20). The four new testsets
take about 10 s with compilation (snapshot probe, at two points per spec: 221 tests in
11.8 s); record the time for the commit. The `N` checks on Segel's ter bi and the power
case match a shape, not a probed text (UNIT §2 printed them in the slice's names,
`$S/phase45/research5/unit/logs/examples.log`, lines 21 and 36): one that fails while its
value check passes is read against those lines before any change, and a difference
P5-10a's graded basis does not explain is a bug.

- [ ] **Step 7: Read the new texts and their value checks**

Create `$S/phase45/scripts/P5_T12_pins.jl`:

```julia
# ABOUTME: Prints the texts Task 12 pins (doctest, factored pairs, reference lines, the
# ABOUTME: allosteric fixture, the flat-sum depths) with each pin's value check, once.
using TestEnv; TestEnv.activate()
using Test, EnzymeRates, LinearAlgebra, Random
const ER = EnzymeRates
include(joinpath(pwd(), "test", "mechanism_definitions_for_test_enzyme_derivation.jl"))
check(name, m) = println("value check ", name, ": ", _testhelper_printed_law_error(m))
mm = @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E(P) ⇌ E + P
    end
end
println(rate_equation_string(mm)); println(rate_equation_string(mm; coefficients = true))
check("doctest", mm)
for spec in MECHANISM_TEST_SPECS
    spec.expected_factored_num === nothing && continue
    v = last(split(rate_equation_string(spec.mechanism), "\n"))
    s = String(chopprefix(v, "v = E_total * "))
    println(spec.name, "\n  ", s)
    check(spec.name, spec.mechanism)
end
fixtures = [
    "half" => @enzyme_mechanism(begin substrates: S; products: P; steps: begin
        (E(S) <--> F(S), E <--> F)
        F(S) <--> E + P
        F + S <--> E(S) end end),
    "folded" => @enzyme_mechanism(begin substrates: A, B; products: P; steps: begin
        (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
        E + B ⇌ E(B)
        E(A) + B ⇌ E(A, B)
        E(A, B) <--> E(P)
        E(P) ⇌ E + P end end),
    "one" => @enzyme_mechanism(begin substrates: A; products: P; modifiers: I
        steps: begin
        E + A <--> E(A)
        E(A) <--> E(P)
        E(P) <--> E + P
        E ⇌ F
        (E + I ⇌ E(I), F + I ⇌ E(I)) end end)]
for (name, m) in fixtures
    println("== ", name); println(rate_equation_string(m)); check(name, m)
    fitted = ER.fitted_params(m)
    rng = Random.MersenneTwister(4)
    p = NamedTuple{(fitted..., :Keq, :E_total)}(
        (Tuple(0.5 + rand(rng) for _ in fitted)..., 2.0, 1.0))
    c = NamedTuple{metabolites(m)}(Tuple(0.5 + rand(rng) for _ in metabolites(m)))
    env = _testhelper_eval_rate_lines(rate_equation_string(m), p, c)
    ap = all_parameters(m, p)
    println("reference lines vs all_parameters: ",
            [(k, Float64(env[k] / ap[k] - 1)) for k in keys(ap) if !(k in keys(p))])
    # The dependents' Exprs, Rational exponents included, read back by value.
    envp = Dict{Symbol, BigFloat}(k => BigFloat(x) for (k, x) in pairs(p))
    dep = ER._dependent_param_exprs(typeof(m))[1]
    println("dependent Exprs vs all_parameters: ",
            [(k, Float64(_testhelper_eval_expr(ex, envp) / ap[k] - 1)) for (k, ex) in dep])
    name == "half" && println("half dep: ", join([string(ex) for (_, ex) in dep], "; "))
end
m_allo = @allosteric_mechanism begin
    substrates: S
    products:   P
    allosteric_regulators: R
    catalytic_multiplicity: 2
    steps: begin
        E + P ⇌ E(P)     :: NonequalAI
        E + S ⇌ E(S)     :: NonequalAI
        E(S) <--> E(P)   :: NonequalAI
    end
    regulatory_site(2): R::NonequalAI
end
println("== m_allo"); println(rate_equation_string(m_allo))
println(rate_equation_string(m_allo; coefficients = true)); check("m_allo", m_allo)
function depth(t)
    d, w = 0, -1
    for ch in t
        d += (ch == '(') - (ch == ')')
        ch in ('+', '-') && (w = max(w, d))
    end
    w
end
function factors(v)
    s = String(chopprefix(v, "v = E_total * "))
    parts, d, from, i = String[], 0, 1, 1
    while i <= lastindex(s)
        d += (s[i] == '(') - (s[i] == ')')
        if d == 0 && (startswith(SubString(s, i), " * (") ||
                      startswith(SubString(s, i), " / ("))
            push!(parts, s[from:prevind(s, i)])
            i += 3
            from = i
            continue
        end
        i = nextind(s, i)
    end
    push!(parts, s[from:end])
    [p[2:(end - 1)] for p in parts]
end
ro = only(s for s in MECHANISM_TEST_SPECS if s.name == "Random-order Bi-Bi").mechanism
println("random-order sum depths (Q, H, D): ",
        depth.(factors(last(split(rate_equation_string(ro), '\n')))))
check("random-order", ro)
println("bpe: ", repr(ER.build_power_expr(-1 // 2, [(:k1f, 1 // 2)])))
```

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T12_pins.log 1800 --project=. $S/phase45/scripts/P5_T12_pins.jl`
Expected: every `value check` line at most 1e-20 (the snapshot probe: at most 6e-77); every
reference line and every dependent `Expr` within 1e-14 of `all_parameters`; a `half dep:`
line whose `k_F_to_E` reads `Keq ^ -1//2` and the other exponents as `±1//2` literals;
`random-order sum depths (Q, H, D): [1, 0, 1]`; a `bpe:` line whose two exponents are the `Rational{Int}` literals `-1//2` and
`1//2`; and these texts, which the next steps pin:
- the doctest mechanism: the docstring's two outputs (Step 4) verbatim;
- the factored pairs, each printed after `v = E_total * ` as `(Q) * (S - P / Keq) / (D)`
  with the Q and D of Step 8;
- `half`'s Haldane line
  `k_F_to_E = k_E_to_F * sqrt(k_FS_to_E_P) * sqrt(k_F_S_to_ES) / (sqrt(k_ES_to_F_S) * sqrt(k_E_P_to_FS) * sqrt(Keq))`;
- `folded`'s lines `# Wegscheider constraints:` and `K_EAB_to_EA_B = K_EB_to_E_B`;
- `one`'s lines `# Wegscheider constraints:` and `K_E_to_F = 1`;
- `m_allo`: the texts of Step 8's pin 3.

A text that differs from these while its value check passes is read against the rules
before any re-pin: term order follows `_term_key` (P5-11), factors follow the variable
order (concentrations in `metabolites(m)` order, then the independents by name, then
`Keq`). A difference those rules explain is a mismatch with Task 11's `_term_key` or with
Task 3's variable order: stop and reconcile them. A value check above 1e-20 is a bug,
not a re-pin.

- [ ] **Step 8: Re-pin the printed laws (pins 1, 3, 5–9, R37, the half-integer exponent)**

*Pin 1, the factored pairs.* In the fixtures file, replace the field comment

```julia
    # Factored form validation (optional — for mechanisms with known factored forms)
    # Expected numerator/denominator from rate_equation_string output
```

with

```julia
    # Factored form validation (optional — for mechanisms with known factored forms):
    # the printed v line's Q (expected_factored_num) and D (expected_factored_denom)
```

and replace each spec's pair (the specs are the eleven whose pairs are set; the old
strings are the ones at f9361cf, unchanged since a9e0699):

| Spec | `expected_factored_num` (Q) | `expected_factored_denom` (D) |
|---|---|---|
| Competitive Inhibitor | `k_ES_to_EP / K_ES_to_E_S` | `1 + S / K_ES_to_E_S + P / K_EP_to_E_P + R / K_ER_to_E_R` |
| Non-competitive Inhibitor | `k_ES_to_EP / K_ES_to_E_S` | `1 + S / K_ES_to_E_S + P / K_EP_to_E_P + R / K_ER_to_E_R + S * R / (K_ER_to_E_R * K_ES_to_E_S)` |
| Uncompetitive Inhibitor | `k_ES_to_EP / K_ES_to_E_S` | `1 + S / K_ES_to_E_S + P / K_EP_to_E_P + S * R / (K_ERS_to_ES_R * K_ES_to_E_S)` |
| Essential Activator | `R * k_ERS_to_EPR / (K_ERS_to_ER_S * K_ER_to_E_R)` | `1 + R / K_ER_to_E_R + S * R / (K_ERS_to_ER_S * K_ER_to_E_R) + P * R / (K_EPR_to_ER_P * K_ER_to_E_R)` |
| Non-essential Activator | `k_ES_to_EP / K_ES_to_E_S + R * k_ERS_to_EPR / (K_ER_to_E_R * K_ES_to_E_S)` | `1 + S / K_ES_to_E_S + P / K_EP_to_E_P + R / K_ER_to_E_R + S * R / (K_ER_to_E_R * K_ES_to_E_S) + P * R / (K_EP_to_E_P * K_ER_to_E_R)` |
| Activator + Competitive Inhibitor | `k_ES_to_EP / K_ES_to_E_S + A * k_EAS_to_EAP / (K_EA_to_E_A * K_ES_to_E_S)` | `1 + S / K_ES_to_E_S + P / K_EP_to_E_P + A / K_EA_to_E_A + I / K_EI_to_E_I + S * A / (K_EA_to_E_A * K_ES_to_E_S) + P * A / (K_EA_to_E_A * K_EP_to_E_P)` |
| Two Competitive Inhibitors | `k_ES_to_EP / K_ES_to_E_S` | `1 + S / K_ES_to_E_S + P / K_EP_to_E_P + I1 / K_EI1_to_E_I1 + I2 / K_EI2_to_E_I2` |
| Two Non-competitive Inhibitors | `k_ES_to_EP / K_ES_to_E_S` | `1 + S / K_ES_to_E_S + P / K_EP_to_E_P + I1 / K_EI1_to_E_I1 + I2 / K_EI2_to_E_I2 + S * I1 / (K_EI1_to_E_I1 * K_ES_to_E_S) + S * I2 / (K_EI2_to_E_I2 * K_ES_to_E_S) + P * I1 / (K_EI1_to_E_I1 * K_EP_to_E_P) + P * I2 / (K_EI2_to_E_I2 * K_EP_to_E_P) + I1 * I2 / (K_EI1_to_E_I1 * K_EI2_to_E_I2) + S * I1 * I2 / (K_EI1_to_E_I1 * K_EI2_to_E_I2 * K_ES_to_E_S) + P * I1 * I2 / (K_EI1_to_E_I1 * K_EI2_to_E_I2 * K_EP_to_E_P)` |
| Non-competitive + Competitive Inhibitor | `k_ES_to_EP / K_ES_to_E_S` | `1 + S / K_ES_to_E_S + P / K_EP_to_E_P + I1 / K_EI1_to_E_I1 + I2 / K_EI2_to_E_I2 + S * I1 / (K_EI1_to_E_I1 * K_ES_to_E_S) + P * I1 / (K_EI1_to_E_I1 * K_EP_to_E_P)` |
| Uncompetitive + Competitive Inhibitor | `k_ES_to_EP / K_ES_to_E_S` | `1 + S / K_ES_to_E_S + P / K_EP_to_E_P + I2 / K_EI2_to_E_I2 + S * I1 / (K_EI1S_to_ES_I1 * K_ES_to_E_S)` |
| Two Same-site Inhibitors | `k_ES_to_EP / K_ES_to_E_S` | `1 + S / K_ES_to_E_S + P / K_EP_to_E_P + I1 / K_EI1_to_E_I1 + I2 / K_EI2_to_E_I2 + S * I1 / (K_EI1_to_E_I1 * K_ES_to_E_S) + S * I2 / (K_EI2_to_E_I2 * K_ES_to_E_S) + P * I1 / (K_EI1_to_E_I1 * K_EP_to_E_P) + P * I2 / (K_EI2_to_E_I2 * K_EP_to_E_P)` |

Write each as the file writes strings now: a value over 92 columns splits into `"…" *`
continuation lines, as the old values do. For example, Competitive Inhibitor's pair
becomes

```julia
        expected_factored_num=
        "k_ES_to_EP / K_ES_to_E_S",
        expected_factored_denom=
        "1 + S / K_ES_to_E_S + P / K_EP_to_E_P + R / K_ER_to_E_R",
```

and Non-essential Activator's

```julia
        expected_factored_num=
        "k_ES_to_EP / K_ES_to_E_S + R * k_ERS_to_EPR / (K_ER_to_E_R * K_ES_to_E_S)",
        expected_factored_denom=
        "1 + S / K_ES_to_E_S + P / K_EP_to_E_P + R / K_ER_to_E_R" *
        " + S * R / (K_ER_to_E_R * K_ES_to_E_S) + P * R / (K_EP_to_E_P * K_ER_to_E_R)",
```

The comments above each spec ("Denom factors as …", the textbook `v = …` lines) describe
the textbook law, which is unchanged, and stay.

In `test/test_rate_eq_derivation.jl`, replace `_testhelper_test_rate_equation_string`
and the extractor below it (from `function _testhelper_test_rate_equation_string` through
the `end` of `_testhelper_extract_num_denom`, its docstring included) with:

```julia
function _testhelper_test_rate_equation_string(spec::_testhelper_MechanismTestSpec)
    m = spec.mechanism
    met_names = spec.metabolite_names
    @testset "Rate Equation String" begin
        s = rate_equation_string(m)
        @test !occursin("params.", s)   # No field-access prefixes
        @test !occursin("concs.", s)    # No field-access prefixes
        @test !occursin("+ -", s)       # No malformed signs
        @test !occursin("- -", s)       # No malformed signs
        @test startswith(last(split(s, "\n")), "v = E_total * ")  # Proper format
        @test occursin(") / (", s)      # Has denominator
        @test occursin("= params", s)   # Has params destructuring
        @test occursin("= concs", s)    # Has concs destructuring
        @test isempty(_testhelper_print_artifacts(s))

        # Numerical equivalence test
        rng = Random.MersenneTwister(9000 + hash(spec.name) % 1000)
        @test all(1:10) do _
            new_params, concs, _ =
                _testhelper_random_independent_params_concs(
                    m, met_names; rng=rng)
            isapprox(
                rate_equation(m, concs, new_params),
                _testhelper_eval_rate_string(s, new_params, concs);
                rtol=1e-10)
        end

        # The eleven pinned specs turn S into P, so H reads S - P / Keq.
        if spec.expected_factored_num !== nothing ||
           spec.expected_factored_denom !== nothing
            q, h, d = _testhelper_extract_factors(last(split(s, "\n")))
            @test h == "S - P / Keq"
            spec.expected_factored_num === nothing || @test q == spec.expected_factored_num
            spec.expected_factored_denom === nothing ||
                @test d == spec.expected_factored_denom
        end
    end
end

"""
The three factors of a printed `v` line, `v = E_total * (Q) * (H) / (D)`, as written
inside their parentheses: `Q` (a bare name as it is, `"1"` when the line leaves it out),
`H` and `D`. Each factor is a parenthesized group or a bare name, so the split runs at the
` * (` and ` / (` outside every parenthesis.
"""
function _testhelper_extract_factors(v_line::AbstractString)
    s = String(chopprefix(v_line, "v = E_total * "))
    parts, depth, from, i = String[], 0, firstindex(s), firstindex(s)
    while i <= lastindex(s)
        depth += (s[i] == '(') - (s[i] == ')')
        if depth == 0 && (startswith(SubString(s, i), " * (") ||
                          startswith(SubString(s, i), " / ("))
            push!(parts, s[from:prevind(s, i)])
            i += 3
            from = i
            continue
        end
        i = nextind(s, i)
    end
    push!(parts, s[from:end])
    inner(t) = startswith(t, "(") && endswith(t, ")") ? t[2:prevind(t, lastindex(t))] : t
    length(parts) == 2 ? ("1", inner(parts[1]), inner(parts[2])) :
                         (inner(parts[1]), inner(parts[2]), inner(parts[3]))
end
```

*Pin 3, "rate_equation_string allosteric byte-identical fixture".* Replace the
assignment `expected = raw"""…"""` and `@test actual == expected` with:

```julia
    expected = raw"""(; K_A_EP_to_E_P, K_A_ES_to_E_S, K_I_EP_to_E_P, K_I_ES_to_E_S, k_A_ES_to_EP, k_I_ES_to_EP, K_A_Rreg, K_I_Rreg, L, Keq, E_total) = params
(; S, P, R) = concs
# Haldane constraints:
k_A_EP_to_ES = K_A_EP_to_E_P * k_A_ES_to_EP / (K_A_ES_to_E_S * Keq)
k_I_EP_to_ES = K_I_EP_to_E_P * k_I_ES_to_EP / (K_I_ES_to_E_S * Keq)
Q_A = k_A_ES_to_EP / K_A_ES_to_E_S
D_A = 1 + S / K_A_ES_to_E_S + P / K_A_EP_to_E_P
Q_I = k_I_ES_to_EP / K_I_ES_to_E_S
D_I = 1 + S / K_I_ES_to_E_S + P / K_I_EP_to_E_P
W_A = (1 + R / K_A_Rreg)^2
W_I = (1 + R / K_I_Rreg)^2
v = E_total * (Q_A * D_A * W_A + L * Q_I * D_I * W_I) * (S - P / Keq) / (D_A^2 * W_A + L * D_I^2 * W_I)"""
    @test actual == expected
    expected_coefficients = raw"""v = E_total * (Q_A * D_A * W_A + L * Q_I * D_I * W_I) * (S - P / Keq) / (D_A^2 * W_A + L * D_I^2 * W_I)
Q_A = nA_0
D_A = 1 + dA_S * S + dA_P * P
Q_I = nI_0
D_I = 1 + dI_S * S + dI_P * P
W_A = (1 + R / K_A_Rreg)^2
W_I = (1 + R / K_I_Rreg)^2
# The coefficients are functions of the fitted constants, not parameters:
nA_0 = k_A_ES_to_EP / K_A_ES_to_E_S
dA_S = 1 / K_A_ES_to_E_S
dA_P = 1 / K_A_EP_to_E_P
nI_0 = k_I_ES_to_EP / K_I_ES_to_E_S
dI_S = 1 / K_I_ES_to_E_S
dI_P = 1 / K_I_EP_to_E_P"""
    @test rate_equation_string(m_allo; coefficients = true) == expected_coefficients
    # Both texts equal the derived law (value check of this pin).
    @test _testhelper_printed_law_error(m_allo) <= big(10)^-20
```

*Pin 5, "rate_equation_string prints every tie as a constraint line".* The three pinned
lines print as before (`K_EAB_to_EA_B = K_EB_to_E_B`, `K_ES_to_EstarS = K_E_to_Estar`,
`K_E_S_to_Estar_P = Keq`): a one-factor ratio of products is the factor. Add the value
check: after each of the three `_testhelper_check_law_against_oracle(…)` lines of the
testset, add `@test _testhelper_reference_lines_match(folded)` (and `iso_tie`, `keq_tie`
in their places), and insert before the testset:

```julia
"""Whether each reference line of `m`'s printed law, read back in 256-bit BigFloat
(`_testhelper_eval_rate_lines`), equals the dependent constant `all_parameters` computes,
to 1e-12 relative, at a random point; false when the law prints no reference line."""
function _testhelper_reference_lines_match(m)
    rng = Random.MersenneTwister(4)
    fitted = EnzymeRates.fitted_params(m)
    p = NamedTuple{(fitted..., :Keq, :E_total)}(
        (Tuple(0.5 + rand(rng) for _ in fitted)..., 2.0, 1.0))
    c = NamedTuple{metabolites(m)}(Tuple(0.5 + rand(rng) for _ in metabolites(m)))
    env = _testhelper_eval_rate_lines(rate_equation_string(m), p, c)
    ap = all_parameters(m, p)
    dependents = [k for k in keys(ap) if !haskey(p, k)]
    !isempty(dependents) &&
        all(k -> abs(env[k] - ap[k]) <= 1e-12 * abs(ap[k]), dependents)
end
```

*Pin 6, "rate_equation_string runs a half-integer dependent as written".* Replace
`@test occursin("0.5", s)` with

```julia
    @test occursin("\nk_F_to_E = k_E_to_F * sqrt(k_FS_to_E_P) * sqrt(k_F_S_to_ES) / " *
                   "(sqrt(k_ES_to_F_S) * sqrt(k_E_P_to_FS) * sqrt(Keq))\n", s)
    @test _testhelper_reference_lines_match(half)
```

and the comment's "through square roots" sentence stays.

*Pin 7, "rate_equation_string prints flat +/* sums (no nested parens)".* Rename the testset
"rate_equation_string prints Q, H and D as flat sums" and replace its body with:

```julia
    # Each of Q, H and D prints as one sum of terms, one per concentration monomial. A
    # coefficient that is a polynomial prints in parentheses inside its term, and none
    # holds another sum, so no sum sits deeper than one parenthesis. Random-order Bi-Bi is
    # steady state throughout, so its coefficients are polynomials; H holds none.
    spec = only(s for s in MECHANISM_TEST_SPECS
                if s.name == "Random-order Bi-Bi")
    s = last(split(rate_equation_string(spec.mechanism), '\n'))
    @test startswith(s, "v = E_total * (")
    q, h, d = _testhelper_extract_factors(s)
    @test (_testhelper_sum_depth(q), _testhelper_sum_depth(h), _testhelper_sum_depth(d)) ==
          (1, 0, 1)
    @test h == "A * B - P * Q / Keq"
    @test _testhelper_printed_law_error(spec.mechanism) <= big(10)^-20
```

and insert before the testset:

```julia
"""The deepest parenthesis level of the printed text `t` that holds a `+` or a `-`; -1
when it holds neither."""
function _testhelper_sum_depth(t::AbstractString)
    depth, worst = 0, -1
    for c in t
        depth += (c == '(') - (c == ')')
        c in ('+', '-') && (worst = max(worst, depth))
    end
    worst
end
```

*Pin 8, "a constant thermodynamics pins to 1 is the dependent constant 1.0".* Rename it "a
constant thermodynamics pins to 1 prints as 1", replace

```julia
        @test occursin("# Wegscheider constraints:\nK_E_to_F = 1.0\n",
                       rate_equation_string(em))
```

with

```julia
        @test occursin("# Wegscheider constraints:\nK_E_to_F = 1\n",
                       rate_equation_string(em))
```

and after the line `params = NamedTuple{syms}(ntuple(i -> 0.5 + i / 10, length(syms)))`
add `@test all_parameters(em, params).K_E_to_F == 1` (its value check).

*Pin 9, the §5a helper.* Replace `_testhelper_undefined_rhs_symbols`'s regex line

```julia
    ident = r"[A-Za-z_][A-Za-z0-9_]*"
```

with

```julia
    # A name may hold superscripts (a coefficient's); a called function (`sqrt`) is no name.
    ident = r"[A-Za-z_][A-Za-z0-9_⁰¹²³⁴⁵⁶⁷⁸⁹⁻]*(?![A-Za-z0-9_⁰¹²³⁴⁵⁶⁷⁸⁹⁻(])"
```

The testset "§5a I-state parameters are all defined (regression)" passes on the three LDH
specs (snapshot probe: no undefined name); its MWC laws define `N_A`, `Q_A`, `D_A`, `W_A`
and the `I` lines before `v` reads them. Value check: Step 6's printed-law testset, where
the three LDH laws, run line by line, equal the derived law to 1e-20.

*Pin 17, "Degenerate constraint handling" → "build_power_expr return types",* landed
in Step 2 and was watched failing in Step 3; its value check is the BigFloat agreement
there.

*R37's printed half, "the reserved names cover every non-metabolite name of a law".* In
its loop body (as Task 6 left it), replace

```julia
        printed = Meta.parseall(rate_equation_string(m))
        others = setdiff(_testhelper_names_read(printed), metabolites(m))
        @test all(n -> occursin('_', String(n)) || n in reserved, others)
        # The printed law destructures params and concs, reads Keq and assigns v.
        @test (:params, :concs, :Keq) ⊆ others
        @test :v in _testhelper_names_read(printed)
```

with

```julia
        printed = Meta.parseall(rate_equation_string(m))
        coefficients = Meta.parseall(rate_equation_string(m; coefficients = true))
        for p in (printed, coefficients)
            others = setdiff(_testhelper_names_read(p), metabolites(m))
            # N, N_A, Q_A, D_A, F_A, W_A, n_0, d_A²_B, nA_0 and the rest hold `_` or are N.
            @test all(n -> occursin('_', String(n)) || n in reserved, others)
            @test :v in _testhelper_names_read(p)
        end
        # The full law destructures params and concs and reads Keq.
        @test (:params, :concs, :Keq) ⊆ setdiff(_testhelper_names_read(printed),
                                               metabolites(m))
```

(snapshot probe: no name outside the reserved set in either form of the 42 specs).

*Pin 4, the doctests.* The `rate_equation_string` docstring's doctest is Step 4's. In
`docs/src/deriving/textbooks.md`, in the block after
`julia> print(rate_equation_string(m))`, replace the two lines

```
k_EP_to_ES = (1 / Keq) * K_EP_to_E_P * (1 / K_ES_to_E_S) * k_ES_to_EP
v = E_total * (k_ES_to_EP * S / K_ES_to_E_S - k_EP_to_ES * P / K_EP_to_E_P) / (1 + P / K_EP_to_E_P + S / K_ES_to_E_S)
```

with

```
k_EP_to_ES = K_EP_to_E_P * k_ES_to_EP / (K_ES_to_E_S * Keq)
v = E_total * (k_ES_to_EP / K_ES_to_E_S) * (S - P / Keq) / (1 + S / K_ES_to_E_S + P / K_EP_to_E_P)
```

Value check: the doctest mechanism's `value check doctest` line of Step 7. Task 15 rewrites
the page's prose.

*The half-integer dependent's exponent.* In `test/test_oracles.jl`, testset "Exact
steady-state oracle" → "a dependent with a half-integer exponent matches the oracle",
replace

```julia
        # The conformational square ties two isomerizations of one group, so the
        # Haldane relation fixes k_F_to_E through square roots (Keq ^ -0.5 · …).
```

with

```julia
        # The conformational square ties two isomerizations of one group, so the
        # Haldane relation fixes k_F_to_E through square roots (Keq ^ -1//2 · …).
```

and

```julia
        @test any(ex -> occursin("^ -0.5", string(ex)), values(dep))
```

with

```julia
        @test any(ex -> occursin("Keq ^ -1//2", string(ex)), values(dep))
```

Value check: Step 7's `dependent Exprs vs all_parameters` line for `half` (the same
mechanism) within 1e-14 and its `half dep:` line; test_joint_solve.jl's "_solve renders a
half-integer dependent with keq -1//2" still compares the `Expr` by value. Only the
literal's spelling moves, from `-0.5` to `-1//2`.

- [ ] **Step 9: Regenerate the golden reference (pin 2)**

In `test/test_allosteric_golden.jl`, replace the ABOUTME lines with

```julia
# ABOUTME: Byte-identical golden reference for the allosteric specs' printed laws, in full
# ABOUTME: and in coefficients, and for their parameter lists.
```

in `_testhelper_allosteric_golden_lines`, after the `push!` of the `RATE_EQUATION_STRING`
line, add

```julia
        push!(lines, "RATE_EQUATION_COEFFICIENTS " *
                     replace(EnzymeRates.rate_equation_string(m; coefficients = true),
                             "\n" => "\\n"))
```

and in the testset's first comment replace the lines

```julia
    # subunit, against the exact steady state of the spec's mass-action network, and the
    # per-spec string check (`_testhelper_test_rate_equation_string`) ties the printed
    # string to the compiled law.
```

with

```julia
    # subunit, against the exact steady state of the spec's mass-action network; the
    # per-spec string check (`_testhelper_test_rate_equation_string`) ties the printed
    # string to the compiled law, and test_print.jl ties both printed forms to the
    # derived law.
```

Create `$S/phase45/scripts/P5_T12_golden.jl`:

```julia
# ABOUTME: Rewrites test/reference/allosteric_golden_reference.txt from the current printer,
# ABOUTME: after test_print.jl has checked every allosteric spec's printed laws (Task 12).
using TestEnv; TestEnv.activate()
using Test, EnzymeRates, LinearAlgebra, Random
include(joinpath(pwd(), "test", "mechanism_definitions_for_test_enzyme_derivation.jl"))
lines = String[]
for spec in MECHANISM_TEST_SPECS
    spec.mechanism isa EnzymeRates.AllostericEnzymeMechanism || continue
    m = spec.mechanism
    push!(lines, "### " * spec.name)
    push!(lines, "RATE_EQUATION_STRING " *
                 replace(EnzymeRates.rate_equation_string(m), "\n" => "\\n"))
    push!(lines, "RATE_EQUATION_COEFFICIENTS " *
                 replace(EnzymeRates.rate_equation_string(m; coefficients = true),
                         "\n" => "\\n"))
end
write("test/reference/allosteric_golden_reference.txt", join(lines, "\n") * "\n")
println(length(lines), " lines written")
```

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T12_golden_gen.log 1800 --project=. $S/phase45/scripts/P5_T12_golden.jl && $S/phase3/jl.sh $S/phase45/logs/P5_T12_golden.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_allosteric_golden.jl")'`
Expected: `33 lines written`; the golden testset passes, `exit 0`. Then
`git diff --stat test/reference/allosteric_golden_reference.txt` shows the file rewritten,
and `grep -c "^RATE_EQUATION_COEFFICIENTS " test/reference/allosteric_golden_reference.txt`
prints `11`. The value check is Step 6's (every allosteric spec's two printed forms equal
the law to 1e-20). MWC Dimer's line, for one, reads
`RATE_EQUATION_STRING (; K_A_EP_to_E_P, K_A_ES_to_E_S, K_I_EP_to_E_P, K_I_ES_to_E_S, k_A_ES_to_EP, k_I_ES_to_EP, L, Keq, E_total) = params\n(; S, P) = concs\n# Haldane constraints:\nk_A_EP_to_ES = K_A_EP_to_E_P * k_A_ES_to_EP / (K_A_ES_to_E_S * Keq)\nk_I_EP_to_ES = K_I_EP_to_E_P * k_I_ES_to_EP / (K_I_ES_to_E_S * Keq)\nQ_A = k_A_ES_to_EP / K_A_ES_to_E_S\nD_A = 1 + S / K_A_ES_to_E_S + P / K_A_EP_to_E_P\nQ_I = k_I_ES_to_EP / K_I_ES_to_E_S\nD_I = 1 + S / K_I_ES_to_E_S + P / K_I_EP_to_E_P\nv = E_total * (Q_A * D_A + L * Q_I * D_I) * (S - P / Keq) / (D_A^2 + L * D_I^2)`
(snapshot probe).

- [ ] **Step 10: Extend Q1 with the printed forms**

In `test/test_oracles.jl`, in Q1 ("the derived laws of the population sample are exact,
positive and balanced"), at the end of the testset `"population mechanism $i"` (after the
detailed-balance line and Task 11's units check), add:

```julia
            # The printed law and its coefficient form of the law above, read back in
            # 256-bit BigFloat, equal it at a consistent point, with no printing artifact.
            full = EnzymeRates._law_text(mech, s, L)
            coef = EnzymeRates._law_text(mech, s, L; coefficients = true)
            @test isempty(_testhelper_print_artifacts(full)) &&
                  isempty(_testhelper_print_artifacts(coef))
            _, pconcs, pparams = _testhelper_law_inputs(mech, s.fitted,
                                                        Random.MersenneTwister(1000 + i))
            law = _testhelper_evaluate_law(L, merge(pparams, pconcs))
            p = merge(pparams, Dict(:E_total => big(1)))
            @test all(t -> abs(_testhelper_eval_rate_string(t, p, pconcs) - law) <=
                           big(10)^-20 * abs(law), (full, coef))
```

and extend the testset's opening comment with the sentence "Both printed forms of each law
equal it at a further point, with no printing artifact." The draws use their own seed, so
the oracle and balance points are those of Task 4.

- [ ] **Step 11: Run the files that read the printer**

Run (background), one after another, each with the focused command (log
`$S/phase45/logs/P5_T12_<what>.log`): `test_rate_eq_derivation` (`tr`), `test_oracles`
(`oracles`), `test_allosteric_collapse` (`collapse`), `test_mechanism_enumeration`
(`tme`), `test_identify_rate_equation` (`tid`), `test_types` (`types`), `test_dsl` (`dsl`),
`test_aqua_jet` (`aqua_jet`):

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && for f in test_rate_eq_derivation:tr test_oracles:oracles test_allosteric_collapse:collapse test_mechanism_enumeration:tme test_identify_rate_equation:tid test_types:types test_dsl:dsl test_aqua_jet:aqua_jet; do $S/phase3/jl.sh $S/phase45/logs/P5_T12_${f#*:}.log 3600 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/'${f%%:*}'.jl")'; done
```

polled with `L=$S/phase45/logs/P5_T12_aqua_jet.log; timeout 580 bash -c "until grep -q '^exit ' $L; do sleep 20; done"; tail -30 $L` (repeat while open).
Expected: every file passes, pristine. In particular:
- TR: the per-spec "Rate Equation String" testsets (pin 1), the byte fixture (pin 3), the
  tie lines (pin 5), the half-integer line (pin 6), the flat sums (pin 7), `K_E_to_F = 1`
  (pin 8), §5a (pin 9), `build_power_expr` (pin 17), the reserved names (R37), "reversing
  written steps changes nothing" (equal texts across rewrites), "Allosteric edge cases"
  (pin 10: `S / K_I_ES_to_E_S` and `P / K_I_EP_to_E_P` sit in `D_I`, snapshot probe),
  "Rate equation too large error" and "constants that would share a name are rejected"
  (their errors come from `_law` and `_solve`, before printing);
- collapse (pin 10): `K_I_ES_to_E_S = K_A_ES_to_E_S` and `k_I_ES_to_E_S = …` print as
  Wegscheider reference lines, so the five string assertions hold unchanged;
- TME (pin 20): "parent and children are rate-equation-distinct" (one law, one text) and
  "rate_equation_string emits section labels" (both section headers) unchanged;
- oracles: "a dependent with a half-integer exponent matches the oracle" finds
  `Keq ^ -1//2` (Step 8) and its oracle checks pass; Q1's printed forms (Step 10);
- TID: every testset keeps its counts; the CSV's `rate_equation` text is the new
  printer's, still printed through the compiled mechanism;
- test_types: the exchange names (`k_EI_J_to_EJ_I`, …) print in the `params` line or a
  reference line, and "Wegscheider" is absent where no Wegscheider dependent exists;
- test_dsl: `(; S, P) = concs` and no `Pmod = concs`;
- aqua_jet: the "Aqua" and "JET" testsets pass. A JET report of a call with no matching
  method is a bug in this task's code: give the caller one method per mechanism type,
  never filter the report.
A pin among these that fails is checked against `_testhelper_printed_law_error` (≤ 1e-20)
before it moves; if it moves, its commit line says why.

- [ ] **Step 12: Build the docs**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T12_docs.log 3600 --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path = pwd())); Pkg.instantiate(); include("docs/make.jl")'`
Expected: the build finishes with no `Error` line and one `Warning`, `deploydocs`'
"Documenter could not auto-detect the building environment. Skipping deployment.": the
docstring's two doctests and textbooks.md's pass, every export has a docstring, and the
`@example` blocks that print laws run.

- [ ] **Step 13: Compare the census with Task 0's**

Create `$S/phase45/scripts/P5_census_nostr.jl`:

```julia
# ABOUTME: Copies two census corpora with every derived law's printed text replaced by "-",
# ABOUTME: so compare.jl checks all but the strings; counts the laws whose text differs.
a, b, a_out, b_out = ARGS
"""Copy corpus `src` to `dst` with each `ok` row's law text blanked; the texts by key."""
function blank(src, dst)
    mkpath(dst)
    for f in readdir(src)
        f == "laws.tsv" || cp(joinpath(src, f), joinpath(dst, f); force = true)
    end
    texts = Dict{String, String}()
    open(joinpath(dst, "laws.tsv"), "w") do io
        for l in eachline(joinpath(src, "laws.tsv"))
            fs = String.(split(l, '\t'))
            if length(fs) >= 8 && fs[5] == "ok"
                texts[fs[1]] = fs[8]
                fs[8] = "-"
            end
            println(io, join(fs, '\t'))
        end
    end
    texts
end
ta, tb = blank(a, a_out), blank(b, b_out)
both = intersect(keys(ta), keys(tb))
println("laws ok in both: ", length(both), "; printed text differs: ",
        count(k -> ta[k] != tb[k], both))
```

Run (background):

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase3 && ./jl.sh $S/phase45/logs/P5_T12_census.log 3600 --project=census/env_law --heap-size-hint=2500M census/corpus.jl $S/phase45/law_wt census/P5_T12 census/rename_phase3.toml && ./jl.sh $S/phase45/logs/P5_T12_census_nostr.log 600 $S/phase45/scripts/P5_census_nostr.jl census/P5_base census/P5_T12 census/P5_base_nostr census/P5_T12_nostr && ./jl.sh $S/phase45/logs/P5_T12_census_cmp.log 900 --project=census/env --heap-size-hint=2500M census/compare.jl census/P5_base_nostr census/P5_T12_nostr census/rename_phase3.toml
```

Expected: the corpus log ends `exit 0`; `laws ok in both: <n>; printed text differs: <n>`
with the two figures equal (every law prints anew); `RESULT: PASS` in
`P5_T12_census_cmp.log` (same init, seed and child sets, the same fitted names, values
within 1e-9), or, when Task 6 recorded summed-cap refusals, a `FAIL` whose only lines are
`<group> laws` lines whose counts add up to the lines of `P5_T6_census_refused.log`. Any
other `FAIL` is a bug in this task: the printer changes no value, name or child.

- [ ] **Step 14: Run the full suite and R52's figures**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T12_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`, polled; then
`S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T12_terter.log 900 --project=. --heap-size-hint=2500M $S/phase45/scripts/P4_terter.jl`
and the TME focused run, `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T12_tme_focused.log 3000 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_mechanism_enumeration.jl")'`.
Expected: `exit 0`, every test passing, pristine, Aqua and JET green, about 8 s over Task
11's time; `ter-ter init: 251107 mechanisms in <t> s; the guard rejects 0`, the count of
Tasks 0 and 6, well under 150 s; the TME run's `elapsed`. Report the three figures
(full-suite `maxrss_kb`, ter-ter seconds, TME seconds) beside Task 0's.

- [ ] **Step 15: Commit**

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt
git add src/print.jl src/EnzymeRates.jl src/rate_eq_derivation.jl \
    src/sym_poly_for_rate_eq_derivation.jl docs/src/deriving/textbooks.md \
    test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_print.jl \
    test/test_rate_eq_derivation.jl test/test_oracles.jl test/test_allosteric_golden.jl \
    test/reference/allosteric_golden_reference.txt
git commit -m "$(cat <<'EOF'
Print the law in its independent constants and in coefficients

rate_equation_string prints v = E_total * (Q) * (H) / (D) over the unit
rule's normalizer, the reference lines as ratios of products with exact
exponents, N when the normalizer has a polynomial factor, and per
conformation Q, D, N, F and W lines for an MWC law. coefficients = true
prints the law in named coefficients with one legend line each. Both forms
of every spec and of the population sample equal the derived law to 1e-20
in 256-bit, with no artifact; the :OnlyI laws equal their limit. The pins
moved after those checks: the eleven factored pairs (Q and D now), the
allosteric fixture and golden file, the doctests, sqrt for the half-integer
dependent, = 1 for a pinned constant, flat sums of depth one,
build_power_expr's Rational exponent, and the half-integer dependent's -1//2
in test_oracles. _rate_text has one method per mechanism type, so JET sees
each law's type. _equation_text and _assert_no_name_clash have no caller
now; the next commit deletes them with the rest of the old engine and
printer, as one unit. Census: values, names and children equal; <n> laws
print anew. R52: maxrss <kb> KB, ter-ter init <s> s, TME <s> s (Task 0:
<figures>).

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
)"
```

---
### Task 12b: PASS1 on the concrete law (P5-18)

identify's PASS1 runs the guard, then derives, keys and prints the law from the concrete
mechanism, so it compiles nothing; PASS2 and LOOCV still compile each mechanism they fit.
A second process prints the same laws in both forms (Review Focus 5).

**Files:**
- Modify: `src/identify_rate_equation.jl` (PASS1 on the concrete law, its result
  `records`; the `_process_batch` docstring)
- Test: `test/test_identify_rate_equation.jl` (the CSV-text check; "PASS1 runs the guard
  before it derives a law"), `test/test_compile_budget.jl` (PASS1's workload and
  `PASS1_INSTANCES`; the second-process testset's printed forms); run only:
  `test/test_aqua_jet.jl`
- Logs: `$S/phase45/logs/P5_T12b_{budget,tid,aqua_jet,full}.log`; `P5_T12b_pass1_trace.jl`
  only when Step 2 finds a PASS1 instance

**Interfaces:**
- Consumes: Task 12's `_law_text(m, s, L; coefficients)` (one method per mechanism type)
  and `rate_equation_string(m; coefficients)`; Task 7's `_canonical`; Task 8's PASS1
  record `_Candidate` (`mech`, `fitted`, `mechanism`, `mechanism_id`, `eq_text`, `canon`)
  and Task 10's `phantom`, `canon0`; Task 5's `_assert_support`; `_check_well_posed`,
  `_solve`, `_law`, `_mechanism_id`. Tests: TID's `_testhelper_CountingStubOpt` and
  `@_testhelper_unguarded`; test_compile_budget.jl's `fit_script`,
  `_testhelper_all_instances`, `_testhelper_instances`,
  `_testhelper_measure_labeled_subprocess` and Task 6's `_testhelper_subprocess_output`.
- Produces: no new src interface. PASS1's result is `records`, a `Vector` of
  `_Candidate`, `FitFailure`, `:complexity_skip` and `:param_skip`. `fit_script`'s
  `_testhelper_measure(m)` returns `(instances, native_bytes, compile_s, show_instances,
  pass1_instances)`.

**Gate:**
- PASS1 runs the guard, then derives, keys and prints from the concrete mechanism, and
  compiles nothing for a fresh mechanism (`PASS1_INSTANCES_k` = 0); the compile budgets
  are read before the full run, and none is raised. A PASS1 reading above 0 stops the task
  (Step 2).
- The second-process testset counts 48 `BODY`, `KEY`, `FULL` and `COEF` lines, each equal
  to the parent's.
- TID is green: "PASS1 runs the guard before it derives a law" records the guard's class-2
  message, and "loocv_results.csv and best_equation.csv saved" reads the selected
  mechanism's printed law back from the CSV.
- `test_aqua_jet.jl` passes before the full run: PASS1's `_law_text` call narrows `m` and
  `L` together (Interface decisions, "Mechanism and law together").
- The full suite is green and pristine; about +1 s.

**Testsets added:**
- test_identify_rate_equation.jl "PASS1 runs the guard before it derives a law".
- test_compile_budget.jl "fit compile: PASS1 compiles nothing for a fresh mechanism".
- test_compile_budget.jl: the second-process testset also compares both printed forms
  (Review Focus 5).

**Pins changed (each after its value check, RES §3.4):** 19, PASS1's workload in
`fit_script` (value check: the budgets read in Step 2 before the full run). TID
"loocv_results.csv and best_equation.csv saved" gains a check that the CSV's text is
`rate_equation_string` of the selected mechanism.

- [ ] **Step 1: PASS1 on the concrete law**

In `src/identify_rate_equation.jl`, `_process_batch`, PASS1 (the `pmap(fresh) do m` block
as Task 10 left it) stops compiling. Replace the lines

```julia
            em = compile_mechanism(m)
            length(fitted_params(em)) > max_param_count && return :param_skip
            eq_text = rate_equation_string(em)
            text = sprint(show, m)
            id = _mechanism_id(text)
            s = _solve(m)
            L = _law(m, s)
            _assert_support(m, L)
```

with

```julia
            _check_well_posed(m)
            s = _solve(m)
            length(s.fitted) > max_param_count && return :param_skip
            L = _law(m, s)
            _assert_support(m, L)
            # `m` and `L` narrow together, so the call matches a `_law_text` method.
            eq_text = m isa AllostericMechanism ? _law_text(m, s, L::_AlloLaw) :
                      _law_text(m, s, L::_Law)
            text = sprint(show, m)
            id = _mechanism_id(text)
```

`m` is an element of a `Union` vector, so the narrowing line keeps `m` and `L` paired, as
Task 10's phantom line does (Interface decisions, "Mechanism and law together"); a bare
`_law_text(m, s, L)` would leave two mismatched pairs for JET.

PASS1 now derives the law itself, so it is a derivation entry and runs the guard first,
as RES §1.1 and R11 ask of every one; the compiled path ran it inside `fitted_params(em)`
(`_derivation_input`), so the guard keeps its place before the parameter cap and a
mechanism it rejects still fails in PASS1 with the guard's message. The enumerator's exits
already drop such mechanisms, so on identify's own candidates the call never throws.

The rest of the block (Task 10's `canon`, `phantom`, `canon0` lines and the
`_Candidate(m, s.fitted, text, id, eq_text, canon, phantom, canon0)` record) stays. Replace
the comment above the block,

```julia
    # PASS 1 (workers): complexity filter (checked first, before any derivation),
    # compile, parameter cap, printed law, then the concrete mechanism's law, its check on
    # the guard's faces and its canonical form. `:complexity_skip` and `:param_skip` are
    # skips; `FitFailure` = threw; else a `_Candidate`.
```

with

```julia
    # PASS 1 (workers), from the concrete mechanism and compiling nothing: the complexity
    # filter (checked first, before any derivation), the guard, the solve and the
    # parameter cap, the law, its check on the guard's faces, the printed law and the
    # canonical form. `:complexity_skip` and `:param_skip` are skips; `FitFailure` = threw;
    # else a `_Candidate`.
```

and rename the block's result `compiled` to `records`, in PASS1 and at every later use in
`_process_batch` (the representative pick, Task 10's phantom loop that reassigns
`compiled[i]`, the skip counts, the leaf-cap log lines, the pre-fit summary and the row
loop): afterwards `grep -nw "compiled" src/identify_rate_equation.jl` prints prose lines
only (no variable `compiled` is left).

In `_process_batch`'s docstring, replace "(a mechanism whose structural `hash` is already in
`seen` is dropped before compiling;" with "(a mechanism whose structural `hash` is already
in `seen` is dropped before it is derived;"; replace "For each fresh mechanism, in
parallel (`pmap`): skip it over `eq_complexity_filter` or `max_param_count`; else compile
it and print its law, take its text and the text's id (`_mechanism_id`), derive its law
from the concrete mechanism (`_solve`, `_law`), check the law on the guard's faces
(`_assert_support`) and label it canonically (`_canonical`)." with "For each fresh
mechanism, in parallel (`pmap`) and from the concrete mechanism alone, without compiling
it: skip it over `eq_complexity_filter`; check it with the guard (`_check_well_posed`);
skip it over `max_param_count` by its solve's fitted constants (`_solve`); else derive its
law (`_law`), check the law on the guard's faces
(`_assert_support`), print it (`_law_text`), take its text and the text's id
(`_mechanism_id`) and label the law canonically (`_canonical`)."; and replace "that threw
at compile, derivation or fit" with "that threw at derivation or fit". `compile_mechanism`
now runs only in PASS2 and LOOCV: `grep -n "compile_mechanism"
src/identify_rate_equation.jl` lists PASS2's `FittingProblem`, `_cv_model_selection`'s
fold fits and its `best_mechanism`.

In `test/test_identify_rate_equation.jl`, testset "loocv_results.csv and best_equation.csv
saved", after `@test nrow(bestf) == 1` add:

```julia
        # The CSV's law is the selected mechanism's printed law: PASS1 prints it from the
        # concrete mechanism, rate_equation_string from the compiled one.
        @test bestf.rate_equation[1] == rate_equation_string(results.best)
```

and after the testset "_process_batch failures report the ORIGINAL mechanism" add:

```julia
@testset "PASS1 runs the guard before it derives a law" begin
    # PASS1 derives the law from the concrete mechanism, so it runs the guard first, as
    # every derivation entry does (R11). Rapid-equilibrium steps that turn S into P beside
    # a steady-state route fail class 2: PASS1 records the guard's message, not the
    # engine's internal error, and fits nothing.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
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
    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            S = [1.0, 2.0, 1.0, 2.0], P = [0.1, 0.2, 0.1, 0.2])
    prob = IdentifyRateEquationProblem(rxn, data; Keq = 2.0)
    opt = _testhelper_CountingStubOpt()
    entries, failures = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[beside], prob;
        optimizer = opt, max_param_count = 20, n_restarts = 1, maxtime = 1.0)
    @test isempty(entries) && opt.count == 0
    @test only(failures).mech == beside
    @test occursin("E + S ⇌ E(S) and E(S) ⇌ E(P) close a cycle", only(failures).error)
    @test !occursin("internal error", only(failures).error)
end
```

The head passes it too, since its PASS1 reaches the guard through `fitted_params(em)`; it
pins that PASS1 keeps the guard once it stops compiling (without Step 1's
`_check_well_posed(m)` line it records the engine's "internal error: a rapid-equilibrium
cycle … (class 2)"). The fixture is Task 3's (test_law.jl "the engine refuses a
rapid-equilibrium cycle that turns the reaction over"), inline.

- [ ] **Step 2: Measure PASS1's compile and compare printed laws across processes (pin 19,
  Review Focus 5)**

In `test/test_compile_budget.jl`, in `fit_script`, replace the whole function
`_testhelper_measure(m)` and the comment above it with:

```julia
        # The method instances and native code bytes that PASS1, PASS2 and one LOOCV
        # fold of m add, their compile time, the show instances that PASS1 and then
        # printing the compiled mechanism add, and the method instances PASS1 alone adds.
        function _testhelper_measure(m)
            n0, b0 = _testhelper_all_instances(), Base.jit_total_bytes()
            t0 = Base.cumulative_compile_time_ns()[1]
            s0 = _testhelper_instances(methods(show))
            ER._eq_complexity(m)
            ER._check_well_posed(m)
            s = ER._solve(m)
            length(s.fitted)
            L = ER._law(m, s)
            ER._assert_support(m, L)
            text = sprint(show, m)
            ER._canonical(L, ER._mechanism_id(text))
            ER._law_text(m, s, L)
            pass1 = _testhelper_all_instances() - n0
            s1 = _testhelper_instances(methods(show))
            fp = FittingProblem(ER.compile_mechanism(m), PROB.data;
                                Keq = PROB.Keq, absolute = PROB.absolute,
                                peak_range = PROB.peak_range)
            fit_rate_equation(fp, OPT; FIT_KW...)
            ER._cv_fold_loss(ER.compile_mechanism(m), PROB, "G1";
                             optimizer = OPT, FIT_KW...)
            fit = (_testhelper_all_instances() - n0, Base.jit_total_bytes() - b0,
                   (Base.cumulative_compile_time_ns()[1] - t0) / 1e9)
            s2 = _testhelper_instances(methods(show))
            sprint(show, ER.compile_mechanism(m))
            (fit..., s1 - s0 + _testhelper_instances(methods(show)) - s2, pass1)
        end
```

and replace

```julia
        for (k, (n, b, t, s)) in enumerate(map(_testhelper_measure, picks[2:end]))
            println("FIT_INSTANCES_", k, ":", n, " FIT_NATIVE_", k, ":", b,
                    " FIT_COMPILE_", k, ":", t, " PASS1_SHOW_", k, ":", s)
        end
```

with

```julia
        for (k, (n, b, t, s, p)) in enumerate(map(_testhelper_measure, picks[2:end]))
            println("FIT_INSTANCES_", k, ":", n, " FIT_NATIVE_", k, ":", b,
                    " FIT_COMPILE_", k, ":", t, " PASS1_SHOW_", k, ":", s,
                    " PASS1_INSTANCES_", k, ":", p)
        end
```

In the comment above the scripts, in the "fit" bullet, replace the parenthetical after
"PASS1" with "(complexity, the guard, the solve and its fitted count, the law and class 5,
the mechanism's text, id and key, and the printed law, all on the concrete mechanism)", and add
after the bullet's sentence on show instances: "PASS1 compiles nothing for a fresh
mechanism, since it reads only the concrete one."

In the parent, replace

```julia
        fit_labels = [string(metric, "_", k)
                      for metric in ("FIT_INSTANCES", "FIT_NATIVE", "FIT_COMPILE",
                                     "PASS1_SHOW")
                      for k in 1:4]
        fit_vals = _testhelper_measure_labeled_subprocess(fit_script,
            [fit_labels; "FIT_SHAPES"; "FIT_NAMES"; "SIG_INSTANCES"; "SIG_NATIVE"])
        # One column per metric, one row per fresh mechanism.
        instances, native, compile_s, pass1_show = eachcol(reshape(fit_vals[1:16], 4, 4))
        n_shapes, n_name_tuples, sig_instances, sig_native = fit_vals[17:20]
```

with

```julia
        fit_labels = [string(metric, "_", k)
                      for metric in ("FIT_INSTANCES", "FIT_NATIVE", "FIT_COMPILE",
                                     "PASS1_SHOW", "PASS1_INSTANCES")
                      for k in 1:4]
        fit_vals = _testhelper_measure_labeled_subprocess(fit_script,
            [fit_labels; "FIT_SHAPES"; "FIT_NAMES"; "SIG_INSTANCES"; "SIG_NATIVE"])
        # One column per metric, one row per fresh mechanism.
        instances, native, compile_s, pass1_show, pass1_instances =
            eachcol(reshape(fit_vals[1:20], 4, 5))
        n_shapes, n_name_tuples, sig_instances, sig_native = fit_vals[21:24]
```

and after the testset "fit compile: PASS1 and a compiled mechanism's text compile no show"
add, inside the same `if isdefined(Core, :methodtable)` branch:

```julia
        # PASS1 reads the concrete mechanism alone (its solve, law, key and printed law),
        # so once the first mechanism has compiled it, a fresh mechanism compiles nothing.
        @testset "fit compile: PASS1 compiles nothing for a fresh mechanism" begin
            @info "PASS1 method instances per fresh mechanism: $pass1_instances (budget: 0)"
            @test maximum(pass1_instances) == 0
        end
```

The subprocess testset (Task 7's "a second process derives the same rate_equation bodies
and keys"): rename it "a second process derives the same rate_equation bodies, keys and
printed laws"; replace its comment, as Task 7 left it,

```julia
    # Q7 (spec section 2, criterion 3): a fresh process reads each mechanism back from its
    # text and derives the same generated body, compared by SHA-256, and the same key, for
    # the 42 specs and the six gate laws.
```

with

```julia
    # Q7 (spec section 2, criterion 3): a fresh process reads each mechanism back from its
    # text and derives the same generated body, compared by SHA-256, the same key and the
    # same printed law in both forms, for the 42 specs and the six gate laws (Review
    # Focus 5: PASS1 derives and prints on the workers).
```

and in its script, after the two `println("KEY ", …)` lines, add

```julia
            println("FULL ", i, " ", bytes2hex(ER.sha256(ER.rate_equation_string(m))))
            println("COEF ", i, " ",
                    bytes2hex(ER.sha256(ER.rate_equation_string(m; coefficients = true))))
```

and in the parent replace

```julia
    @test count(startswith("BODY "), lines) == 48 && count(startswith("KEY "), lines) == 48
```

with

```julia
    @test all(t -> count(startswith(t * " "), lines) == 48, ("BODY", "KEY", "FULL", "COEF"))
```

and add at the end of its loop body:

```julia
        full = bytes2hex(EnzymeRates.sha256(rate_equation_string(m)))
        coef = bytes2hex(EnzymeRates.sha256(rate_equation_string(m; coefficients = true)))
        @test "FULL $i $full" in lines
        @test "COEF $i $coef" in lines
```

Printing the six gate laws costs about 0.1 s per process (snapshot probe: the 858-term
ter-bi law, 28 ms in full and 24 ms in coefficients).

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T12b_budget.log 1800 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_compile_budget.jl")'`
Expected: every testset passes; the `[ Info` lines read the trace count (≤ 200), the
first-call wall clock (≤ 6.0 s), the Sig budgets (≤ 250, ≤ 100 KB), the fit compile per
fresh mechanism (≤ 250 instances, ≤ 300 KB; record them beside Task 8's readings: PASS2
still compiles each fresh type), `show instances per fresh mechanism: [0.0, 0.0, 0.0,
0.0]` and `PASS1 method instances per fresh mechanism: [0.0, 0.0, 0.0, 0.0]` (the parsed
readings are `Float64`); the second-process testset counts 48 of each line. Record every
reading for the commit. A PASS1 reading above 0 stops the task: rerun the fit script with
`--trace-compile=$S/phase45/logs/P5_T12b_pass1_trace.jl` on the four picks and find the
methods compiled for a fresh mechanism. A method specialized on something that differs
by mechanism is a bug: remove the specialization and rerun. A method compiled because a
fresh pick takes a branch the warm-up pick does not (no mechanism-dependent type in its
signature) is not one: stop and report the trace to Denis. Never raise the budget.

- [ ] **Step 3: Run identify's tests and JET**

Run (background), one after the other, the focused command on
`test_identify_rate_equation` (`tid`) and `test_aqua_jet` (`aqua_jet`):

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && for f in test_identify_rate_equation:tid test_aqua_jet:aqua_jet; do $S/phase3/jl.sh $S/phase45/logs/P5_T12b_${f#*:}.log 3600 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/'${f%%:*}'.jl")'; done
```

polled with `L=$S/phase45/logs/P5_T12b_aqua_jet.log; timeout 580 bash -c "until grep -q '^exit ' $L; do sleep 20; done"; tail -30 $L` (repeat while open).
Expected: both pass, pristine. TID: "PASS1 runs the guard before it derives a law" records
the guard's class-2 message and fits nothing; "loocv_results.csv and best_equation.csv
saved" reads the selected mechanism's printed law back from the CSV; every other testset
keeps its counts. aqua_jet: the "Aqua" and "JET" testsets pass; a JET report of a call
with no matching method is a bug in Step 1's code, never a report to filter.

- [ ] **Step 4: Run the full suite**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T12b_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`, polled.
Expected: `exit 0`, every test passing, pristine, Aqua and JET green, about 1 s over Task
12's time; report `maxrss_kb` against Task 0's.

- [ ] **Step 5: Commit**

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl \
    test/test_compile_budget.jl
git commit -m "$(cat <<'EOF'
Derive, key and print PASS1's law from the concrete mechanism

PASS1 runs the guard, then derives the law, checks it on the guard's faces,
prints it and keys it from the concrete mechanism, so it compiles nothing:
PASS1 method instances per fresh mechanism <readings>. PASS2 and LOOCV
still compile each mechanism they fit. The CSV's rate_equation is the
selected mechanism's rate_equation_string, and a second process prints the
same laws in both forms. Budgets: <readings>. Full-suite maxrss: <kb> KB
(Task 0: <kb0> KB).

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
)"
```

---
### Task 13: Delete the old engine and printer

Every consumer reads the new law from Task 12b on. This task converts the last tests that
read the old engine's polynomials, checks them against the old engine while it still
exists, then deletes it.

**Files:**
- Delete: `src/sym_poly_for_rate_eq_derivation.jl`. `MAX_RATE_EQUATION_TERMS` moves to
  `src/engine.jl` and `build_power_expr` to `src/thermo.jl`; `MONO`, `POLY`, `_mono`,
  `poly_zero`, `poly_one`, `poly_sym`, `_poly_addop`, `poly_add`, `poly_sub`, `poly_mul`,
  `_mono_mul`, `_reduce_conc_lowest_terms`, `sym_det`, `_poly_to_expr`, `_nest_binary`,
  `_op_prec`, `_str_paren`, `_expr_to_string` and `_mentions` go (`_poly_to_horner` went
  in Task 6).
- Modify: `src/rate_eq_derivation.jl`: delete `_segment_root`, `_re_weight_ratio`,
  `_compute_alpha`, `_ss_contrib`, `_raw_symbolic_rate_polys` (three methods),
  `_num_den_exprs` (both), `_assert_no_name_clash`, `_equation_text`,
  `_state_rate_polys`, `_reg_site_expr`, `_power_expr`, `_mwc_combine`,
  `_mwc_cross_weight`, `_invert_monomial`, `_is_metabolite_free_monomial`,
  `_mwc_state_polys`, the MWC comment block and the emptied section lines; rewrite the
  ABOUTME.
- Modify: `src/admission.jl` (the `POLY` methods of `_law_branches`,
  `_keeps_monomial`, `_vanishing_face`, `_support_verdict` and `_assert_support`, with their
  docstrings); `src/types.jl` (the `_inactive_mask` and `_active_mask` docstrings);
  `src/engine.jl` (`MAX_RATE_EQUATION_TERMS`; `_assert_derivable`'s and `_AlloLaw`'s
  docstrings);
  `src/thermo.jl` (`build_power_expr`); `src/EnzymeRates.jl` (the include line)
- Modify: `ROADMAP.md` (any pointer to the deleted file or functions)
- Test: `test/test_rate_eq_derivation.jl`, `test/allosteric_ground_truth.jl`,
  `test/test_mechanism_enumeration.jl`, `test/test_admission.jl` (its two
  `_num_den_exprs` throw assertions go with their subject), the fixtures file
  (`_testhelper_steady_state_rate`'s `occupancy`; the two-argument `_testhelper_poly_value`
  on `POLY` goes; comments naming the old engine)
- Create (scratch): `$S/phase45/probes/P5_T13_den_check.jl`; census `census/P5_T13`
- Logs: `$S/phase45/logs/P5_T13_{den_check,tr_conv,agt_conv,tme_conv,grep,tr,agt,tme,admission,law,print,full,census,census_cmp,terter,tme_focused}.log`

**Interfaces:**
- Consumes: every earlier task's code; tests: Task 3's `_testhelper_poly_value`,
  phase 4's `_testhelper_law_inputs`, `_testhelper_steady_state_rate`.
- Produces: no new src interface. `_testhelper_steady_state_rate(m, point, concs; Λ, flip,
  occupancy::Int = 0)`: with `occupancy > 0`, the steady-state fraction of the enzyme in
  active form `occupancy` (a `BigFloat`). After the task,
  `grep -rnE "\bPOLY\b|\bMONO\b|\b_mono(\(|_mul\b)|sym_det|_num_den_exprs|_raw_symbolic_rate_polys|_state_rate_polys|_mwc_(combine|cross_weight|state_polys)|_invert_monomial|_is_metabolite_free_monomial|_segment_root|_re_weight_ratio|_compute_alpha|_ss_contrib|_assert_no_name_clash|_equation_text|_expr_to_string|_poly_to_(expr|horner)|_nest_binary|_op_prec|_str_paren|_poly_addop|\b_power_expr\b|_reg_site_expr|_reduce_conc_lowest_terms|\bpoly_(one|zero|sym|add|sub|mul)\b|\b_mentions\b|_rate_eq_dedup_key|eq_hash|sym_poly_for_rate_eq_derivation" src test ROADMAP.md`
  prints nothing.

**Gate:**
- Before any deletion, `P5_T13_den_check.jl` finds the old and the new denominators carry
  the same concentration exponents on every mechanism it reads, and the converted tests pass
  against the old engine's presence (Steps 1–3).
- The grep above prints nothing.
- The full suite is green and pristine, with Aqua and JET.
- `census/P5_T13` against `census/P5_T12`: compare.jl `RESULT: PASS`, strings included.
- R52's figures beside Task 0's.
- The deleted subjects go to the PR body's list for Denis (Task 15).

**Testsets added:** none.

**Pins changed (11, and two more):**
- TR "the free-enzyme weight is free E's weight in the reduced denominator": on `_Law`,
  `L.trees[1] * L.shift` is `B`'s term. Value check, added to the testset: the weight over
  `D` is free E's steady-state share of the enzyme (`occupancy = 1`), to 1e-20.
- TR "rendering helpers": its subjects go; Task 4's `_cross_weights` testset carries its
  three cases. The testset is removed and listed for Denis. Value check: Task 4's MWC laws,
  built with those weights, equal the oracle (test_law.jl).
- TR "a segment's root ties go to the canonical form order" becomes "each
  rapid-equilibrium segment is rooted at its first form": `_segment_root`'s rule goes with
  the old engine, and the engine roots each segment at its first form (`_segments`).
  Listed for Denis. Value check: the fixture's law equals the oracle
  (`_testhelper_check_law_against_oracle`).
- TR `_testhelper_classify_dep_expr` reads `_testhelper_mentions`, since `_mentions` goes.
- allosteric_ground_truth "OnlyA MWC derivation matches mass-action ground truth": the
  free-E weights read `A.active` and `A.inactive`, and `A.fA == A.fI == one(_Poly)`. Value
  check: the testset's own ground-truth check passes.
- allosteric_ground_truth "ping-pong :OnlyA I-state keeps a reachable free-enzyme root":
  the inactive flip weight from `_law(am).inactive`. Value check: the same testset's oracle
  check (`_testhelper_check_law_against_oracle(half2)`).
- TME "zero inactive-state numerator keeps the :NonequalAI binding constants":
  `iszero(_law(am).inactive.Q)`.
- TME "_hyperbolic_catalysis" and "_hyperbolic_catalysis matches the derived
  denominator": read `_law(m).D`, `_law(am).active.D` and `.inactive.D`. Value check:
  Step 1's comparison, then the converted testsets' verdicts with the old engine still
  present (Step 3).
- TR "Allosteric edge cases": the comment that names `_num_den_exprs`'s flat path is
  rewritten.

- [ ] **Step 1: Compare the old and the new denominators while both exist**

Create `$S/phase45/probes/P5_T13_den_check.jl`:

```julia
# ABOUTME: Before Task 13 deletes the old engine: the concentration exponents of every
# ABOUTME: denominator and free-E weight, old against new, per conformation; 0 differences.
using TestEnv; TestEnv.activate()
using Test, EnzymeRates, LinearAlgebra, Random
const ER = EnzymeRates
include(joinpath(pwd(), "test", "mechanism_definitions_for_test_enzyme_derivation.jl"))

"""The concentration-exponent vectors of an old `POLY` over the names `mets`."""
old_exps(p, mets) = Set([Int(get(Dict(mono), x, 0)) for x in mets] for mono in keys(p))
"""The concentration-exponent vectors of a new `_Poly` over the names `mets`."""
new_exps(p, vars, mets) = Set([Int(numerator(ER._mexp(mono, vars.index[x]))) for x in mets]
                               for mono in keys(p.terms))

"""Old against new: the denominators and free-E weights per conformation of `m`."""
function differs(m)
    mets = collect(ER.metabolites(m))
    if m isa ER.Mechanism
        _, den, d_free = ER._raw_symbolic_rate_polys(m)
        L = ER._law(m)
        return old_exps(den, mets) != new_exps(L.D, L.vars, mets) ||
               old_exps(d_free, mets) != new_exps(L.trees[1] * L.shift, L.vars, mets)
    end
    A = ER._law(m)
    any(((:A, A.active), (:I, A.inactive))) do (state, C)
        _, den, d_free = ER._state_rate_polys(m, state)
        old_exps(den, mets) != new_exps(C.D, C.vars, mets) ||
            old_exps(d_free, mets) != new_exps(C.trees[1] * C.shift, C.vars, mets)
    end
end

bibi = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    allosteric_regulators: R
    oligomeric_state: 2
end
pingpong = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    oligomeric_state: 2
end
sets = [
    "specs" => [_testhelper_lift(s.mechanism) for s in MECHANISM_TEST_SPECS],
    "bi-bi init" => collect(ER.init_mechanisms(@enzyme_reaction(begin
        substrates: A[C], B[N]; products: P[C], Q[N] end))),
    "bi-bi + R level 1" => unique!(ER.expand_mechanisms(
        collect(ER.init_mechanisms(bibi))[1:20], bibi)),
    "ping-pong level 1" => unique!(ER.expand_mechanisms(
        collect(ER.init_mechanisms(pingpong))[1:20], pingpong)),
]
for (name, ms) in sets
    ms = filter(m -> ER._eq_complexity(m) <= 337, ms)
    bad = filter(differs, ms)
    println(name, ": ", length(ms), " mechanisms (",
            count(m -> m isa ER.AllostericMechanism, ms), " allosteric), ",
            length(bad), " differ")
    foreach(m -> println("  ", sprint(show, m)), first(bad, 3))
end
```

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T13_den_check.log 1800 --project=. --heap-size-hint=2500M $S/phase45/probes/P5_T13_den_check.jl`
Expected: four lines, each with `0 differ`: `specs: 42 mechanisms (11 allosteric)`,
`bi-bi init: 239 mechanisms (0 allosteric)`, and the two level-1 lines with allosteric
children (snapshot probe, on Tasks 1–4's engine: 439 with 340 allosteric, and 399 with 300
allosteric; f9361cf's populations give the same counts, and the zeros are the gate).
Both engines bring their polynomials to lowest terms in the concentrations by the
same shift, so the exponent sets agree; a difference is a bug in the new law to find before
any deletion.

- [ ] **Step 2: Convert the tests that read the old engine**

In the fixtures file, `_testhelper_steady_state_rate`: change the signature line to

```julia
function _testhelper_steady_state_rate(m, point, concs::AbstractDict; Λ = big(10)^30,
                                       flip::Int = 1, occupancy::Int = 0)
```

replace

```julia
        x = M \ rhs
        sum((ω * (rf * x[a] - rr * x[b]) for (a, b, rf, rr, ω) in edges); init = big(0))
```

with

```julia
        x = M \ rhs
        occupancy == 0 || return x[node[(:A, occupancy)]]
        sum((ω * (rf * x[a] - rr * x[b]) for (a, b, rf, rr, ω) in edges); init = big(0))
```

and in its docstring replace the signature line
`    _testhelper_steady_state_rate(m, point, concs; Λ = big(10)^30, flip = 1) -> BigFloat`
with
`    _testhelper_steady_state_rate(m, point, concs; Λ = big(10)^30, flip = 1, occupancy = 0)`
and append to the docstring's first paragraph the sentence "With `occupancy`, it returns
instead the steady-state fraction of the enzyme in active form `occupancy`."

In the same file, delete the two-argument `_testhelper_poly_value(p, env)` with its
docstring ("The value of the polynomial `p` at the values `env` (by symbol), in
BigFloat."): it sums the old engine's `POLY` dictionaries, and since Task 4 no caller reads
it. `grep -rn "_testhelper_poly_value(" test/` then lists only three-argument calls on
`_Poly`.

In `test/test_rate_eq_derivation.jl`:

(a) Replace the testset "a segment's root ties go to the canonical form order" with:

```julia
@testset "each rapid-equilibrium segment is rooted at its first form" begin
    # The engine weighs each form of a rapid-equilibrium segment against the segment's
    # first form in the canonical form order (`_segments`), so free E, form 1, roots its
    # own segment. A's modifier copy comes before Ab there: the table sorts A's copy just
    # after A, before Ab. The law in lowest terms does not depend on the root.
    em = @enzyme_mechanism begin
        substrates: A, Ab
        products: P
        modifiers: A
        steps: begin
            E + A ⇌ E(A)
            E(A) + Ab ⇌ E(A, Ab)
            E(A, Ab) <--> E(P)
            E(P) ⇌ E + P
            E + A::Mod ⇌ E(A::Mod)
            E + Ab ⇌ E(Ab)
            E(Ab) + A ⇌ E(A, Ab)
            E(A) <--> E(; residual = A - P) + P
        end
    end
    m = ER.Mechanism(em)
    form(text) = only(f for f in eachindex(m.forms) if ER._form_name(m, f) == text)
    a_copy, ab, e = form("EAmod"), form("EAb"), form("E")
    @test e == 1 && a_copy < ab
    s = ER._segments(m, trues(length(m.edges)), ER._flags(m))
    @test all(seg -> seg[1] == minimum(seg), s.segments)
    @test s.segments[s.segment_of[e]][1] == e
    @test s.segment_of[a_copy] == s.segment_of[ab] == s.segment_of[e]
    covalent = form("E_res_plusA_minusP")
    @test s.segments[s.segment_of[covalent]] == [covalent]
    _testhelper_check_law_against_oracle(em)
end
```

(b) In the testset "the free-enzyme weight is free E's weight in the reduced denominator",
replace the comment's last clause "; the MWC combination divides each conformation by this
weight." with "; the MWC law reads this weight as free E's (`_cross_weights`).", and
replace

```julia
    _, _, d_free = EnzymeRates._raw_symbolic_rate_polys(typeof(em))
    @test d_free == EnzymeRates.poly_sym(:B)
```

with

```julia
    m = EnzymeRates.Mechanism(em)
    s = EnzymeRates._solve(m)
    L = EnzymeRates._law(m, s)
    @test L.trees[1] * L.shift ==
          EnzymeRates._term(EnzymeRates._monomial(L.vars.index[:B] => 1))
    # Free E's weight over the denominator is its steady-state share of the enzyme.
    point, concs, params = _testhelper_law_inputs(m, s.fitted, Random.MersenneTwister(8))
    env = merge(params, concs)
    share = setprecision(BigFloat, 256) do
        _testhelper_poly_value(L.trees[1] * L.shift, L.vars, env) /
            _testhelper_poly_value(L.D, L.vars, env)
    end
    @test abs(share - _testhelper_steady_state_rate(m, point, concs; occupancy = 1)) <=
          big(10)^-20 * share
```

(c) In `_testhelper_classify_dep_expr`, replace `EnzymeRates._mentions(expr, :Keq)` with
`_testhelper_mentions(expr, :Keq)`, and insert above its docstring:

```julia
"""Whether the expression `ex` mentions the symbol `s`."""
_testhelper_mentions(ex, s::Symbol) =
    ex === s || (ex isa Expr && any(a -> _testhelper_mentions(a, s), ex.args))
```

(d) In "Allosteric edge cases", replace the comment

```julia
    # Regression: T-state binding K's must be in Kd convention even when
    # `:OnlyA` and `:NonequalAI` catalytic groups coexist. Without the fix,
    # the flat-poly path in _num_den_exprs renders T-state K's
    # as `K_T * met` (Ka) instead of `met / K_T` (Kd), silently producing
    # wrong rates whenever a mechanism mixes these two tags. Regression
    # for src/rate_eq_derivation.jl:1395-1396.
```

with

```julia
    # The inactive conformation's binding constants are dissociation constants (`S / K`)
    # also when `:OnlyA` and `:NonequalAI` catalytic groups coexist; read as association
    # constants (`K * S`) they would halve the rate below.
```

(e) Delete the testset "rendering helpers" (from its `@testset` line through its `end`).

In `test/allosteric_ground_truth.jl`:

(f) In "OnlyA MWC derivation matches mass-action ground truth", replace "so its
free-enzyme graph is a single segment and both states surface D = 1." with "so its
free-enzyme graph is a single segment and both conformations weigh free E 1.", and replace

```julia
    am = ER.AllostericMechanism(onlyA)
    _, _, dA = ER._state_rate_polys(am, :A)
    _, _, dI = ER._state_rate_polys(am, :I)
    @test dA == ER.poly_one()                       # single active segment
    @test dI == ER.poly_one()                       # single inactive segment
```

with

```julia
    am = ER.AllostericMechanism(onlyA)
    A = ER._law(am)
    @test A.active.trees[1] * A.active.shift == one(ER._Poly)       # single active segment
    @test A.inactive.trees[1] * A.inactive.shift == one(ER._Poly)   # one inactive segment
    @test A.fA == A.fI == one(ER._Poly)
```

(g) In "ping-pong :OnlyA I-state keeps a reachable free-enzyme root", replace

```julia
    _, _, d_free_I = ER._state_rate_polys(am, :I)
    @test !isempty(d_free_I)
```

with

```julia
    inactive = ER._law(am).inactive
    @test !iszero(inactive.trees[1] * inactive.shift)
```

In `test/test_mechanism_enumeration.jl`:

(h) In "zero inactive-state numerator keeps the :NonequalAI binding constants", replace

```julia
    @test isempty(first(
        EnzymeRates._state_rate_polys(EnzymeRates.AllostericMechanism(m), :I)))
```

with

```julia
    @test iszero(EnzymeRates._law(EnzymeRates.AllostericMechanism(m)).inactive.Q)
```

(i) In "_hyperbolic_catalysis", replace

```julia
    _testhelper_den_a_degree(k) = maximum(get(Dict(mono), :A, 0) for mono in keys(
        EnzymeRates._raw_symbolic_rate_polys(k)[2]))
```

with

```julia
    _testhelper_den_a_degree(k) = (L = EnzymeRates._law(k);
        maximum(numerator(EnzymeRates._mexp(mono, L.vars.index[:A]))
                for mono in keys(L.D.terms)))
```

and

```julia
    _testhelper_inactive_degree(am, x) = maximum(get(Dict(mono), x, 0) for mono in keys(
        EnzymeRates._state_rate_polys(am, :I)[2]))
```

with

```julia
    _testhelper_inactive_degree(am, x) = (C = EnzymeRates._law(am).inactive;
        maximum(numerator(EnzymeRates._mexp(mono, C.vars.index[x]))
                for mono in keys(C.D.terms)))
```

(j) In "_hyperbolic_catalysis matches the derived denominator", replace

```julia
    hyperbolic(p, mets) = all(e <= 1 for mono in keys(p) for (s, e) in mono if s in mets)
```

with (two lines, so the helper takes the `_testhelper_` prefix)

```julia
    _testhelper_hyperbolic(C, mets) = all(e <= 1 for mono in keys(C.D.terms)
                                          for (v, e) in mono if C.vars.names[v] in mets)
```

replace

```julia
            hyperbolic(EnzymeRates._state_rate_polys(m, state)[2], mets)
```

with

```julia
            L = EnzymeRates._law(m)
            _testhelper_hyperbolic(state === :A ? L.active : L.inactive, mets)
```

and replace

```julia
                _, den, _ = EnzymeRates._raw_symbolic_rate_polys(m)
                @test structural == hyperbolic(den, mets)
```

with

```julia
                @test structural == _testhelper_hyperbolic(EnzymeRates._law(m), mets)
```

`_state_mask(m, state)`, the cache's key, is `_active_mask(am)` for `:A` and
`_inactive_mask(am)` for `:I`, the masks `_law` derives each conformation on, so the cache
still holds one derivation per conformation graph.

(k) Run `grep -n "_num_den_exprs\|_raw_symbolic\|_state_rate_polys\|_mwc_state_polys\|POLY" test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_admission.jl test/test_oracles.jl`
and rewrite each comment or docstring it prints to name the new law (`EnzymeRates._law`)
in place of the old function; earlier tasks rewrote the code that called them, so only
prose is left. The exception is test_admission.jl's two `@test_throws` on
`EnzymeRates._num_den_exprs(unfed)` and `EnzymeRates._num_den_exprs(am)`, kept by Task 5's
ruling on item 16: delete both; their testsets keep the `_assert_support(…, _law(…))`
assertions.

- [ ] **Step 3: Run the converted tests against the old engine's presence**

Run (background), one after another:

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && for f in test_rate_eq_derivation:tr_conv allosteric_ground_truth:agt_conv test_mechanism_enumeration:tme_conv; do $S/phase3/jl.sh $S/phase45/logs/P5_T13_${f#*:}.log 3600 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/'${f%%:*}'.jl")'; done
```

polled. Expected: all three pass, pristine. "_hyperbolic_catalysis matches the derived
denominator" keeps `n_checked > 100` and `0 < n_nonhyperbolic < n_checked`; the
`_hyperbolic_catalysis` degrees `[1, 2, 2, 3]` and the inactive degrees 2 and ≤ 1 hold as
before; the free-E share agrees to 1e-20. These are the pins' value checks: the old engine
still exists, so a verdict that moved here is a conversion bug.

- [ ] **Step 4: Delete the old engine**

(a) In `src/engine.jl`, insert before the docstring of `RateEquationTooLargeError`:

```julia
"""
The term cap: the most spanning-tree products (V×τ, `_eq_complexity`) a derived rate law
may have, both conformations of an MWC law together. A larger law takes too long to
compile and is unlikely to be useful for fitting.
"""
const MAX_RATE_EQUATION_TERMS = 5000
```

in `_assert_derivable(m::Mechanism, mask)`'s docstring replace "Every derivation entry
calls it (`_law`, `_raw_symbolic_rate_polys`)." with "Every derivation entry calls it
(`_law`).", and append to `_AlloLaw`'s docstring the paragraph that keeps the facts of the
MWC comment block (d) deletes:

```
The law is per active site: `E_total` is the active-site concentration, so `n` enters
only as an exponent (of `DA`, `DI` and the free-E weights), never as a leading factor, and
each regulatory site's factor enters numerator and denominator at the site's own
multiplicity, whatever `n` is.
```

(b) Move `build_power_expr`, docstring, code and comment unchanged, from
`src/sym_poly_for_rate_eq_derivation.jl` to `src/thermo.jl`, after the
`_dependent_param_exprs` methods. As Task 12 left it, it reads:

```julia
"""
    build_power_expr(keq_exp::Rational, factors)

Build a power-product expression for thermodynamic-constraint
substitution of the form `Keq^keq_exp * prod(k_i^exp_i)` from
`factors`, an iterable of `(k_i, exp_i)` pairs. Used by the
constraint solver to materialize the Keq×rate-constant power product
that replaces a dependent rate constant. A fractional exponent is an
exact `Rational{Int}` literal. The empty product is the Float64 `1.0`: a constant that
thermodynamics pins to 1.
"""
function build_power_expr(keq_exp::Rational, factors)
    _pa(sym, exp) = exp == 1 ? sym : exp == -1 ? :(1 / $sym) :
        !isinteger(exp) ? :($sym ^ $(Rational{Int}(exp))) :
        exp > 0 ? :($sym ^ $(Int(exp))) : :(1 / $sym ^ $(Int(-exp)))
    terms = Any[]
    keq_exp != 0 && push!(terms, _pa(:Keq, keq_exp))
    # Sort factors by name so the rendered power product is content-canonical
    # (independent of the order factors were collected during elimination):
    # two mechanisms with the same constraint produce the identical RHS string.
    for (sym, exp) in sort(collect(factors); by = x -> string(first(x)))
        exp != 0 && push!(terms, _pa(sym, exp))
    end
    isempty(terms) ? 1.0 : length(terms) == 1 ? only(terms) : Expr(:call, :*, terms...)
end
```

(c) `git rm src/sym_poly_for_rate_eq_derivation.jl`, and in `src/EnzymeRates.jl` delete the
line `include("sym_poly_for_rate_eq_derivation.jl")`.

(d) In `src/rate_eq_derivation.jl`, delete each of these with its docstring and the comment
lines that belong to it: `_segment_root`, `_re_weight_ratio`, `_compute_alpha`,
`_ss_contrib`, the three `_raw_symbolic_rate_polys` methods, `_num_den_exprs(m::Mechanism)`,
`_assert_no_name_clash`, `_equation_text`, `_state_rate_polys`, `_reg_site_expr`,
`_power_expr`, `_mwc_combine`, `_mwc_cross_weight`, `_invert_monomial`,
`_is_metabolite_free_monomial`, `_mwc_state_polys` and `_num_den_exprs(am::AllostericMechanism)`;
the MWC comment block (from `# ═══…` "AllostericEnzymeMechanism rate equations (MWC)" to
its closing `# ═══…` line); and the section lines left with nothing under them (`RE Group
Helpers`, `Expr generation from POLY`, `String Representation`). Rename the section line
`# ─── Raw Rate Equation Derivation (Unified Cha / King-Altman) ───` to
`# ─── Edge masks and constant names ───` (`_state_mask` and `_edge_names` stay under it),
rename `# ─── Rate body building helpers ───` to `# ─── Regulator constants ───`, and move
`_kreg_params` (docstring and method) from below `_ss_rate_constant_names` to just below
that line, beside `_reg_K`; pad each renamed line with `─` to the width of the others.
Replace the two ABOUTME lines with

```julia
# ABOUTME: The parameter API (parameters, fitted_params, all_parameters) and the constant
# ABOUTME: names of each conformation's law: edge masks, edge names and regulator constants.
```

The file keeps `parameters`, `fitted_params` and their concrete methods, `all_parameters`,
`_state_mask`, `_edge_names`, `_ss_rate_constant_names`, `_kreg_params` and `_reg_K`.
Check:
`grep -nE "^(function |[A-Za-z_]+\(.*\) =|@generated)" src/rate_eq_derivation.jl` lists
those functions' definitions and nothing else.

(e) In `src/admission.jl`, section "Class 5: the support test on the derived law", delete
the `POLY` methods with their docstrings: the docstring headed
`_law_branches(den) -> Vector{Tuple{POLY, POLY}}` and
`_law_branches(den_A, den_I, d_A, d_I) -> Vector{Tuple{POLY, POLY}}` with
`_law_branches(den::POLY)` and `_law_branches(den_A::POLY, den_I::POLY, d_A::POLY,
d_I::POLY)`; `_keeps_monomial(p::POLY, z::Vector{Symbol})`;
`_vanishing_face(branches::Vector{Tuple{POLY, POLY}}, faces::Vector{Vector{Symbol}})`;
both `_support_verdict(…, branches::Vector{Tuple{POLY, POLY}})` methods; and
`_assert_support(m::_AnyMech, branches::Vector{Tuple{POLY, POLY}})`. Task 5's methods on
`_Poly`, `_Law` and `_AlloLaw`, each with its own docstring, stay.

(f) In `src/types.jl`, in the docstrings of `_inactive_mask` and `_active_mask`, replace
`` (`_state_rate_polys`) `` with `` (`_law`) ``.

(g) In `ROADMAP.md`, run
`grep -n "sym_poly_for_rate_eq_derivation\|_num_den_exprs\|_raw_symbolic\|_state_rate_polys" ROADMAP.md`
and rewrite each mention: "the rate polynomials (`src/sym_poly_for_rate_eq_derivation.jl`)"
(Idea "Count identifiable parameters with an exact modular rank", if Task 10's rewrite kept
it) becomes "the derived law (`_law`, `src/engine.jl`)"; any other mention names `_law` the
same way.

- [ ] **Step 5: Check that nothing reads the old engine**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && grep -rnE "\bPOLY\b|\bMONO\b|\b_mono(\(|_mul\b)|sym_det|_num_den_exprs|_raw_symbolic_rate_polys|_state_rate_polys|_mwc_(combine|cross_weight|state_polys)|_invert_monomial|_is_metabolite_free_monomial|_segment_root|_re_weight_ratio|_compute_alpha|_ss_contrib|_assert_no_name_clash|_equation_text|_expr_to_string|_poly_to_(expr|horner)|_nest_binary|_op_prec|_str_paren|_poly_addop|\b_power_expr\b|_reg_site_expr|_reduce_conc_lowest_terms|\bpoly_(one|zero|sym|add|sub|mul)\b|\b_mentions\b|_rate_eq_dedup_key|eq_hash|sym_poly_for_rate_eq_derivation" src test ROADMAP.md > $S/phase45/logs/P5_T13_grep.log; echo "lines: $(wc -l < $S/phase45/logs/P5_T13_grep.log)"`
Expected: `lines: 0`. A line printed is a reader the earlier steps missed: convert it, then
rerun.

- [ ] **Step 6: Run the touched test files**

Run (background), one after another, with the focused command: `test_rate_eq_derivation`
(`tr`), `allosteric_ground_truth` (`agt`), `test_mechanism_enumeration` (`tme`),
`test_admission` (`admission`), `test_law` (`law`), `test_print` (`print`), logs
`$S/phase45/logs/P5_T13_<what>.log`, polled as in Task 12 Step 11.
Expected: every file passes, pristine, with Step 3's counts.

- [ ] **Step 7: Run the full suite, the census and R52's figures**

Run (background), one after another:

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T13_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'
cd $S/phase3 && ./jl.sh $S/phase45/logs/P5_T13_census.log 3600 --project=census/env_law --heap-size-hint=2500M census/corpus.jl $S/phase45/law_wt census/P5_T13 census/rename_phase3.toml && ./jl.sh $S/phase45/logs/P5_T13_census_cmp.log 900 --project=census/env --heap-size-hint=2500M census/compare.jl census/P5_T12 census/P5_T13 census/rename_phase3.toml
```

then the ter-ter and TME figures as in Task 12 Step 14 (logs `P5_T13_terter.log`,
`P5_T13_tme_focused.log`).
Expected: the full suite `exit 0`, every test passing, pristine, Aqua and JET green, its
test count Task 12's minus the removed "rendering helpers" tests (8 at f9361cf) and the
two `_num_den_exprs` throw assertions plus the converted testsets' new assertions;
`RESULT: PASS` with strings exact (`strings exact <n>/<n>`): deleting the old engine
changes no child, value, name or text; the three R52 figures beside Task 0's.

- [ ] **Step 8: Commit**

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt
git add src/EnzymeRates.jl src/engine.jl src/thermo.jl src/rate_eq_derivation.jl \
    src/admission.jl src/types.jl ROADMAP.md \
    test/mechanism_definitions_for_test_enzyme_derivation.jl \
    test/test_rate_eq_derivation.jl test/allosteric_ground_truth.jl \
    test/test_mechanism_enumeration.jl test/test_admission.jl test/test_oracles.jl
git status --short
git commit -m "$(cat <<'EOF'
Delete the old derivation engine and printer

POLY, sym_det, the raw symbolic rate polynomials, the MWC Expr assembly and
the old printer go, and the test suite's POLY evaluator with them;
MAX_RATE_EQUATION_TERMS moves to engine.jl, build_power_expr to thermo.jl,
and the MWC comment block's per-active-site facts to _AlloLaw's docstring. Before the deletion the old and new
denominators and free-E weights carried the same concentration exponents
on the specs, the bi-bi seeds and two level-1 populations, and the
converted tests passed beside the old engine. The free-E weight pin now
also checks the weight against free E's steady-state share. "rendering
helpers" goes with its subjects (Task 4's _cross_weights tests carry its
cases); the segment-root test pins the engine's rooting at a segment's
first form. Census: equal to Task 12's, strings included. R52: maxrss
<kb> KB, ter-ter init <s> s, TME <s> s (Task 0: <figures>).

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
)"
```

`git status --short` before the commit shows `D  src/sym_poly_for_rate_eq_derivation.jl`
staged by Step 4's `git rm`, the files above, and nothing untracked. A file that Step 5's
conversions touched beyond this list (a docstring in `src/codegen.jl`, say) shows as
` M`: add it by name before the commit.

---
### Task 14: The size gate's constant on CI (P5-24)

c, the size gate's nanoseconds per denominator term, follows phase 1's rule on the PR's
CI: two runs on the three runners, twice the worst reading over the six gate laws, rounded
up to one significant figure. Until now it stays 0.8 ns; the research predicts 1.0.

**Files:**
- Modify: `test/test_rate_eq_derivation.jl` (`_testhelper_seconds_per_denominator_term`
  and the comment above it; the headroom comment above `@test t < 120e-9` in
  `_testhelper_test_performance`)
- Modify: `.claude/CLAUDE.md` (the "`rate_equation` runtime perf is non-negotiable"
  paragraph: `0.8 ns × denominator terms`, `0.8e-9 * n`, and the headroom sentence's
  slowest-spec reading, from Task 6's measurement; P5-24 ◆, first among the CLAUDE.md
  edits listed for Denis)
- Modify: `docs/src/developer.md` (`0.8 ns × denominator terms`)
- Modify: `ROADMAP.md` (the known issue "**Ubuntu and windows CI miss `rate_equation`'s
  120 ns gate on one spec.**", which 3543abe added and 7de259f reworded for phase 5's code
  generation to settle)
- Create (scratch): `$S/phase45/pr_body5.md` (provisional), `$S/phase45/scripts/P5_T14_c.jl`
- Logs: `$S/phase45/logs/P5_T14_{ci1_raw,ci1,ci2_raw,ci2,c,tr,ci3_raw,ci3}.log`

**Interfaces:**
- Consumes: Task 6's gate testset "rate_equation time per call grows at most in proportion
  to the law's size" and its `@info` line `rate_equation: <t> ns per call, <n> denominator
  terms, <x> ns per term`; the CI workflow (`.github/workflows/CI.yml`: `pull_request` and
  `workflow_dispatch`; runners ubuntu-latest x64, windows-latest x64, macOS-latest aarch64,
  Julia 1).
- Produces: c (`<c>` below: the figure Step 5 computes) and the draft PR.

**Gate:**
- The controller pushes graph-law and opens the draft PR (base main, P5-48) with the
  provisional body; CI runs twice, green on all three runners both times.
- c is twice the worst ns per term over the two runs, rounded up to one significant figure.
  A c above 0.8 is said in the commit, and the CLAUDE.md edit and the looser bound go to
  Denis in the PR body; a c at or below 0.8 tightens the gate and goes to him the same way.
- A third CI run, on the new c, is green.
- A spec that misses 120 ns on a runner is rerun once, not loosened (CLAUDE.md). The same
  spec missing again on that runner, as the LDH variant P4-37 kept does on ubuntu and
  windows unless Task 6 cut it, stops the task for Denis with the readings; Task 6 Step 8's
  local reading of that law predicts it (at about 50 ns or more, CI misses again).

**Testsets added:** none.

**Pins changed:** TR "rate_equation time per call grows at most in proportion to the law's
size": c. Value check: every CI reading of the two runs sits under the new bound (c is
twice the worst).

- [ ] **Step 1: Write the provisional PR body**

Create `$S/phase45/pr_body5.md`:

````markdown
## Summary

Phase 5 of the reaction-graph rewrite (spec `docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md`,
sections 9, 10, 11 and the phase-5 parts of 13 and 14). One exact engine derives every rate
law on phase 3's graph and phase 4's masks; `rate_equation` is generated from it, the
printer prints it in the independent constants, and identify fits one mechanism per law.

- **Engine** (`src/poly.jl`, `src/engine.jl`): the law of one conformation on an edge mask
  in checked `Int64`; the MWC law on `_active_mask` and `_inactive_mask`; the 5,000 cap on
  both conformations together.
- **Code generation** (`src/codegen.jl`): locals for the dependents, the segment sums, `Q`,
  `H` and the MWC terms; inputs read by name.
- **Unit rule and printer** (`src/normalize.jl`, `src/print.jl`): dimensionless
  denominator terms; `rate_equation_string(m; coefficients = true)`.
- **Key and identify** (`src/key.jl`, `src/identify_rate_equation.jl`): one fit per key,
  rows under their own names, the beam and LOOCV per key, L phantoms that inherit a fit;
  CSV columns `key`, `fit_inherited`, `fit_source_id`, `L_phantom`.
- **Deleted:** the old engine and printer (`POLY`, `sym_det`, `_num_den_exprs`, the MWC
  `Expr` assembly, `eq_hash`).
- Version stays 0.10.0. Draft, stacked on the phase-4 PR; not for merge.

**CI calibration in progress (P5-24):** this run and one rerun set the size gate's constant
c from the gate laws' ns per term on the three runners; a third run checks the new c.

The full body (rulings, questions for Denis, measurements) follows in Task 15.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
````

- [ ] **Step 2 (controller): Push graph-law and open the draft PR**

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && git status --short && git push -u origin graph-law && gh pr create --draft --base main --head graph-law --title "Phase 5 of the reaction-graph rewrite: one law engine, its printer, its key" --body-file $S/phase45/pr_body5.md
```

Expected: a clean tree; the push succeeds; `gh pr create` prints the PR's URL. Record its
number `<n>` for Step 7 and Task 15.

- [ ] **Step 3 (controller): Read the first CI run**

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt
SHA=$(git rev-parse HEAD)
timeout 600 bash -c "until gh run list --branch graph-law --workflow CI.yml --commit $SHA --json databaseId --jq '.[0].databaseId' | grep -q '[0-9]'; do sleep 20; done"
RUN=$(gh run list --branch graph-law --workflow CI.yml --commit $SHA --limit 1 --json databaseId --jq '.[0].databaseId'); echo "run $RUN"
gh run watch $RUN --exit-status --interval 60 > /dev/null; echo "exit $?"
gh run view $RUN --log > $S/phase45/logs/P5_T14_ci1_raw.log
grep -P '^Julia [^\t]*\t.*rate_equation: [0-9.]+ ns per call' $S/phase45/logs/P5_T14_ci1_raw.log > $S/phase45/logs/P5_T14_ci1.log; wc -l < $S/phase45/logs/P5_T14_ci1.log
```

Expected: `exit 0` (every job green, the documentation job included); 18 lines in
`P5_T14_ci1.log`, six per runner. `gh run watch` blocks; run the block in the background
and poll for `exit ` in its output. A job that fails on a 120-ns spec bound is rerun once
with `gh run rerun $RUN --failed`, not loosened; the same spec missing again stops the task
for Denis (Gate); any other failure stops the task.

- [ ] **Step 4 (controller): Rerun and read the second CI run**

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt
SHA=$(git rev-parse HEAD)
RUN=$(gh run list --branch graph-law --workflow CI.yml --commit $SHA --limit 1 --json databaseId --jq '.[0].databaseId'); echo "run $RUN"
gh run rerun $RUN
timeout 300 bash -c "until [ \"\$(gh run view $RUN --json status --jq .status)\" != completed ]; do sleep 10; done"
gh run watch $RUN --exit-status --interval 60 > /dev/null; echo "exit $?"
gh run view $RUN --log > $S/phase45/logs/P5_T14_ci2_raw.log
grep -P '^Julia [^\t]*\t.*rate_equation: [0-9.]+ ns per call' $S/phase45/logs/P5_T14_ci2_raw.log > $S/phase45/logs/P5_T14_ci2.log; wc -l < $S/phase45/logs/P5_T14_ci2.log
```

(`RUN` is Step 3's run, found again by the head commit: each shell call starts fresh.)
Expected: `exit 0` and 18 lines. The wait loop lets the rerun start before `gh run watch`
reads the run, and `gh run view` reads the latest attempt of a rerun run.

- [ ] **Step 5: Compute c**

Create `$S/phase45/scripts/P5_T14_c.jl`:

```julia
# ABOUTME: The size gate's constant from CI logs: each runner's worst ns per term over the
# ABOUTME: six gate laws (bi-bi, ter-bi), and c, twice the worst rounded up to one figure.
worst = Dict{Tuple{String, String}, Float64}()
for file in ARGS, line in eachline(file)
    m = match(r"^(Julia [^\t]*)\t.*rate_equation: [0-9.]+ ns per call, (\d+) denominator " *
              r"terms, ([0-9.]+) ns per term", line)
    m === nothing && continue
    family = parse(Int, m[2]) == 396 ? "bi-bi" : "ter-bi"
    worst[(m[1], family)] = max(get(worst, (m[1], family), 0.0), parse(Float64, m[3]))
end
for k in sort!(collect(keys(worst)))
    println(k[1], " ", k[2], ": ", worst[k], " ns per term")
end
x = 2 * maximum(values(worst))
e = floor(Int, log10(x))
c = round(ceil(x / 10.0^e) * 10.0^e; sigdigits = 1)
println("worst ", maximum(values(worst)), " ns per term; c = ", c, " ns per term")
```

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; $S/phase3/jl.sh $S/phase45/logs/P5_T14_c.log 300 --startup-file=no $S/phase45/scripts/P5_T14_c.jl $S/phase45/logs/P5_T14_ci1.log $S/phase45/logs/P5_T14_ci2.log; cat $S/phase45/logs/P5_T14_c.log`
Expected: six lines, one per runner and family,
then `worst <w> ns per term; c = <c> ns per term`. The research predicts c = 1.0 (local
worst 0.108 times the worst CI/local ratio 4.58, doubled).

- [ ] **Step 6: Set c**

In `test/test_rate_eq_derivation.jl`, replace the comment and constant as Task 6 Step 2
left them,

```julia
# Seconds per call allowed per denominator term: each of the six gate laws (the fixtures
# file) must run in under the larger of 120 ns and this constant times its denominator
# terms. It covers
# the worst CI readings on each runner, in ns per term, ter-bi / bi-bi: ubuntu x64
# 0.305 / 0.313, windows x64 0.346 / 0.459, macOS aarch64 0.248 / 0.732.
const _testhelper_seconds_per_denominator_term = 0.8e-9
```

with (each `<…>` the figure of `P5_T14_c.log`)

```julia
# Seconds per call allowed per denominator term: each of the six gate laws (the fixtures
# file) must run in under the larger of 120 ns and this constant times its denominator
# terms. It is twice the worst CI reading over two runs, rounded up to one significant
# figure; the worst readings on each runner, in ns per term, ter-bi / bi-bi: ubuntu x64
# <a> / <b>, windows x64 <c1> / <d>, macOS aarch64 <e> / <f>.
const _testhelper_seconds_per_denominator_term = <c>e-9
```

In `.claude/CLAUDE.md`, in the paragraph under "### `rate_equation` runtime perf is
non-negotiable", replace `max(120 ns, 0.8 ns × denominator terms)` with
`max(120 ns, <c> ns × denominator terms)` and `max(120e-9, 0.8e-9 * n)` with
`max(120e-9, <c>e-9 * n)`. In the same paragraph replace the stale headroom figures,

```markdown
Locally the 120-ns bound leaves about 2× headroom over the slowest `MECHANISM_TEST_SPECS` law ("LDH i-state NonequalAI 6-group, steady-state NAD group", about 55 ns) and about 4× over the next two (about 30 ns: the tetramer MWC laws PFK-1 and "MWC Tetramer Random Bi-Bi RE + Two Allosteric Sites"); shared CI runners run these 2-7× slower: the ubuntu and windows runners read about 130 ns on the LDH law and miss the bound (a ROADMAP known issue), and a job that misses it on another law is rerun rather than the bound loosened.
```

with the two slowest specs Task 6 Step 8 recorded (`P5_T6_hoist.log`; `<k>` is 120
divided by the slower reading, rounded down)

```markdown
Locally the 120-ns bound leaves about <k>× headroom over the slowest `MECHANISM_TEST_SPECS` laws (about <t> ns: <spec 1> and <spec 2>); shared CI runners can run a law 4-7× slower, so rerun a job that misses the bound rather than loosening it.
```

In `test/test_rate_eq_derivation.jl`, `_testhelper_test_performance`, rewrite the comment
above `@test t < 120e-9` (it names the LDH law at about 55 ns and 2× headroom, P4-37 (a))
with the same figures, in its own words: the headroom over the two slowest specs Task 6
Step 8 recorded, that shared CI runners can run a law 4-7× slower, and that such a job is
rerun rather than the bound loosened; keep its last sentence, "The 0-alloc bound above is
the strict part."

In `docs/src/developer.md`, "Derivation architecture", replace
`0.8 ns × denominator terms` with `<c> ns × denominator terms`.

In `ROADMAP.md`, under "Known issues", the bullet "**Ubuntu and windows CI miss
`rate_equation`'s 120 ns gate on one spec.**" (through its last text line, keeping the
blank line before "## Planned"): when
Task 6 Step 8's slowest spec reads about 30 ns or less (the 4× headroom the entry wants
back) and neither CI run of Steps 3–4 missed 120 ns on a spec, delete it; otherwise
rewrite it in the present tense with Task 6 Step 8's slowest spec and reading and the
worst CI reading on it. The commit body says which.

Run `grep -rn "0\.8 ns\|0\.8e-9" --include=*.md --include=*.jl . | grep -v docs/superpowers`:
it prints nothing (when c is 0.8, it prints the three lines that keep it), and
`grep -c "about 30 ns" .claude/CLAUDE.md` prints `0`.

- [ ] **Step 7: Check locally and commit**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T14_tr.log 3600 --project=. -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_rate_eq_derivation.jl")'`
Expected: passes, pristine; the six `[ Info: rate_equation: …` lines sit under the new
bound.

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt
git add test/test_rate_eq_derivation.jl .claude/CLAUDE.md docs/src/developer.md ROADMAP.md
git commit -m "$(cat <<'EOF'
Set the size gate's constant from two CI runs on the new code

c is twice the worst ns per denominator term over two CI runs on the
three runners, rounded up to one significant figure: <c> ns (worst <w>,
<runner>). It replaces 0.8 ns, set on the old code, in the gate and in
CLAUDE.md<, which loosens the bound: Denis decides in the PR>.

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
)"
```

The bracketed clause stays only when c is above 0.8; when it is below, write "which
tightens the bound: Denis decides in the PR" in its place. Add to the body: "CLAUDE.md's
headroom sentence reads the new code's slowest specs, <spec 1> and <spec 2>, at about <t>
ns." When c is 0.8, the bound stays: commit the comment's new readings and CLAUDE.md's
headroom sentence, with the summary "Record the size gate's CI readings on the new code"
and the body's second sentence "c stays 0.8 ns."

- [ ] **Step 8 (controller): Push and read the third CI run**

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && git push
SHA=$(git rev-parse HEAD)
timeout 600 bash -c "until gh run list --branch graph-law --workflow CI.yml --commit $SHA --json databaseId --jq '.[0].databaseId' | grep -q '[0-9]'; do sleep 20; done"
RUN3=$(gh run list --branch graph-law --workflow CI.yml --commit $SHA --limit 1 --json databaseId --jq '.[0].databaseId'); echo "run $RUN3"
gh run watch $RUN3 --exit-status --interval 60 > /dev/null; echo "exit $?"
gh run view $RUN3 --log > $S/phase45/logs/P5_T14_ci3_raw.log
grep -P '^Julia [^\t]*\t.*rate_equation: [0-9.]+ ns per call' $S/phase45/logs/P5_T14_ci3_raw.log > $S/phase45/logs/P5_T14_ci3.log; wc -l < $S/phase45/logs/P5_T14_ci3.log
```

Expected: `exit 0` and 18 lines, each reading under `max(120, <c> × n)` ns. Task 15 quotes
the three runs' readings.

---

### Task 15: Docs, ROADMAP, PR body, final GC gate, Julia 1.10

The pages phase 5 makes false (spec section 16: `deriving/textbooks.md`,
`deriving/ping_pong.md`, the identify pages and `developer.md`, plus every page that shows a
printed law or an API this phase changes), the ROADMAP entries it settles, the GC gate on
the final head, the Julia 1.10 run, and the final PR body. Earlier phases and the fix wave
rewrote these pages, so each edit names its paragraph by heading and content, as it stands
at f9361cf.

**Files:**
- Modify: `docs/src/deriving/textbooks.md`, `docs/src/deriving/ping_pong.md`,
  `docs/src/deriving/mwc_allostery.md`, `docs/src/deriving/iso.md`,
  `docs/src/deriving/cha_king_altman.md`, `docs/src/deriving/thermodynamic_constraints.md`,
  `docs/src/identify/tutorial.md`, `docs/src/identify/model_selection.md`,
  `docs/src/identify/enumeration_engine.md`, `docs/src/developer.md`,
  `README.md`, `docs/src/getting_started.md`
- Modify: `ROADMAP.md`
- Create (scratch): `$S/phase45/pr_body5.md` (final); `$S/phase45/scripts/P5_mkenv110.jl`
  and the environment `$S/phase45/env110`
- Logs: `$S/phase45/logs/P5_T15_{docs,full,scan112_*,scan113_*,repro112,repro113,hoistscan,v110_env,v110_poly,v110_law,v110_key,v110_print,v110_tr,v110_oracles,v110_tid}.log`

**Interfaces:**
- Consumes: Task 0's `p5_base/figures.md`; Task 4's recount; Task 6's gate readings,
  hoisting probe and scans; Task 7's key cost; Task 9's identify run; Task 10's rerun and
  counts; Task 12's and 13's census and R52 figures; Task 12b's compile budgets; Task 14's
  c and CI readings.
- Produces: the final PR body, which the controller pushes
  (`gh api … -X PATCH`).

**Gate:**
- `docs/make.jl` builds with its doctests and its check that every export has a docstring,
  with no warning but `deploydocs`' "Skipping deployment."
- The GC gate on the final head reports 0 hoisted allocations with a bits store into a
  reference field on Julia 1.12.5 and 1.13.1, with both positive controls firing.
- Julia 1.10.12 (P5-38), in an environment of its own (`$S/phase45/env110`):
  `test_poly.jl`, `test_law.jl`, `test_key.jl`, `test_print.jl`,
  `test_rate_eq_derivation.jl`, `test_oracles.jl` and `test_identify_rate_equation.jl`
  pass in focused runs, one at a time.
- The full suite is green and pristine.
- The PR body lists every item of this plan's "Open for Denis".

**Testsets added:** none. **Pins changed:** none.

- [ ] **Step 1: Locate the paragraphs**

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt
grep -n "The string has four sections\|The final .v = E_total" docs/src/deriving/textbooks.md
grep -n "signature of a ping-pong\|residual-bearing forms" docs/src/deriving/ping_pong.md
grep -n "For a mechanism with catalytic multiplicity\|In the printed equation above" docs/src/deriving/mwc_allostery.md
grep -n "S·P. cross" docs/src/deriving/iso.md
grep -n "fewest bound\|cofactor determinants\|division-free" docs/src/deriving/cha_king_altman.md
grep -n "(1 / Keq)\|K_E_to_F = 1.0" docs/src/deriving/thermodynamic_constraints.md
grep -n "partition function,\|names its mechanism by its text\|roughly the" docs/src/identify/tutorial.md
grep -n "Width floor\|min_beam_width. | .50\|complexity bound roughly\|eq_hash\|fit_inherited\|mechanism_id. |" docs/src/identify/model_selection.md
grep -n "phantom all the same\|whose equation is the parent's\|On its bases the engine emits" docs/src/identify/enumeration_engine.md
grep -n "keeps flat sums\|_dependent_param_exprs. (.src/thermo.jl.)\|engine derives each conformation\|_num_den_exprs\|eq_hash" docs/src/developer.md
grep -n "rate_equation_string" README.md docs/src/getting_started.md
```

Expected: each paragraph Steps 2–9 replace is found once. The pages wrap their prose, so
a sentence a step quotes can break across lines (the grep anchors on one line of it):
read each hit's neighbourhood. A paragraph the fix wave reworded is replaced wherever it
now stands; its new text is the one given here.

- [ ] **Step 2: `docs/src/deriving/textbooks.md`**

Replace from the line `The string has four sections:` through the list item
`4. The final \`v = E_total * (num) / (den)\` line.` with:

````markdown
The string reads, in order:

1. a `params` destructuring line: the independent fitted constants, then `Keq` and
   `E_total`;
2. a `concs` destructuring line: the concentrations;
3. the reference lines under `# Haldane constraints:` (after `# Wegscheider constraints:`
   when a cycle of bindings ties constants together): each dependent constant in `Keq` and
   the independent constants, for reference, since the law reads only the independent
   ones;
4. when the normalizer has a factor that is not a single term, the line `N = …`
   ([Normalization](@ref));
5. the law, `v = E_total * (Q) * (H) / (D)`. `H = S - P / Keq` vanishes at equilibrium,
   and `Q` and `D` have positive terms.

## Normalization

Every term of `D` is dimensionless: the package divides numerator and denominator by a
normalizer with `D`'s unit, chosen by the first rule that applies.

1. **Constant term.** When `D` has a term free of concentrations, divide by its
   coefficient. Each term is then an enzyme form's abundance relative to free E, as in
   `1 + S / K_ES_to_E_S + P / K_EP_to_E_P` above.
2. **Virtual constant term.** Otherwise take `D`'s terms in a fixed order (lowest total
   degree first, then fewest metabolites other than substrates, then metabolite order) and
   weigh the first affinely independent ones so that the product of their coefficients,
   each to its weight, has `D`'s unit. A ping-pong bi-bi gets Cleland's normalizer this way
   ([Ping-pong mechanisms](@ref)); fractional weights print as roots.
3. **Ratio form.** When no such weights exist, the law reads some concentrations only
   through their ratio, as a rapid-equilibrium Theorell–Chance step `E(A) + B ⇌ E(Q) + P`
   makes every term carry `B` or `P`. Divide by the first term in the same order; terms
   such as `P / B` are then dimensionless.

A factor of the normalizer that is not a single term prints once, as the line `N = …`, and
the terms it divides read `/ N`.

## The coefficient form

`coefficients = true` prints the same law with one name per coefficient, then a legend
that gives each coefficient in the fitted constants:

```jldoctest textbook
julia> print(rate_equation_string(m; coefficients = true))
v = E_total * n_0 * (S - P / Keq) / (1 + d_S * S + d_P * P)
# The coefficients are functions of the fitted constants, not parameters:
n_0 = k_ES_to_EP / K_ES_to_E_S
d_S = 1 / K_ES_to_E_S
d_P = 1 / K_EP_to_E_P
```

A coefficient is named by the concentrations it multiplies: `d_A_B` for `A * B`, `d_A²_B`
for `A^2 * B`, and, in the ratio form, `d_P_B⁻¹` for `P / B`; `n_0` is the numerator's
constant term. The fit still runs over the fitted constants.
````

- [ ] **Step 3: `docs/src/deriving/ping_pong.md`**

Replace the last paragraph, from "The derived equation has the signature of a ping-pong
mechanism:" through "(with `B` additionally bound).", with:

````markdown
The derived equation has the signature of a ping-pong mechanism: its denominator has **no
constant term**. A sequential mechanism always carries a `1` in its denominator for the
free enzyme; a ping-pong enzyme never sits idle as free `E` during turnover (it is always
either loaded with substrate or carrying the covalent residual), so that term is absent.

Without a constant term the printed law divides by a virtual one ([Normalization](@ref)).
The first affinely independent terms of the fixed order are those in `A`, `B`, `P`, `Q`
and `A * B`, and the weights that give the normalizer the denominator's unit fall on `A`,
`B` and `A * B` alone: the coefficients of `A` and `B` over that of `A * B`, Cleland's
normalizer for a ping-pong bi-bi. Its part that is not a single term
prints once as `N = …`, and every denominator term is dimensionless. In other mechanisms,
ter-reactant ones above all, the fixed order can pick terms other than Cleland's, which
changes how `N` reads but never the law's value.

The coefficient form names each coefficient by the concentrations it multiplies:

```@example pingpong
print(rate_equation_string(m; coefficients = true))
```

The residual-bearing forms take descriptive names such as `E_res_plusA_minusP` (the
modified free enzyme) and `EB_res_plusA_minusP` (with `B` additionally bound).
````

- [ ] **Step 4: `docs/src/deriving/mwc_allostery.md`**

In section "Derivation of the MWC rate equation", replace from the paragraph that begins
"For a mechanism with catalytic multiplicity `cat_n`, the package normalizes each
conformation" through the paragraph that ends "with `(1 + I / K_I_Ireg) ^ 2`." (the formula
block, the paragraph on `N_A`, `Q_A`, `D_A`, the paragraph on the free-enzyme weight `D`
with phase 4's sentence on `N_I = 0`, and the paragraph on the printed equation) with:

````markdown
For a mechanism with catalytic multiplicity `n`, the package derives each conformation's
law on the graph it keeps (the active conformation without its `:OnlyI` groups, the
inactive one without its `:OnlyA` groups, each without the forms that cuts off from free
E) and combines the two:

```
v = E_total * (H) * (F_I^n * Q_A * D_A^(n-1) * W_A + L * F_A^n * Q_I * D_I^(n-1) * W_I) /
                    (F_I^n * D_A^n * W_A + L * F_A^n * D_I^n * W_I)
```

`H = Π[S] − Π[P]/Keq` is the Haldane factor both conformations share. `Q_A` is the active
conformation's numerator over `H` and `D_A` its denominator, the King–Altman/Cha
polynomials of its catalytic graph; `F_A` is free E's weight in `D_A`, the spanning-tree
sum of free E's segment times the concentration monomial that brings the law to lowest
terms. `Q_I`, `D_I` and `F_I` are the inactive conformation's, with `:EqualAI` groups
keeping the shared constants and `:NonequalAI` groups their own `I` names. `W_A` and `W_I`
multiply the regulatory-site factors `(1 + X/K)^m`, and `L` weighs the inactive branch.
Each branch carries the other conformation's free-E weight, so both stand on free E, where
the conformations interconvert; that keeps a bare rate constant out of the `L` term.

The printed law lists these polynomials by name. Each conformation's `Q`, `D` and `F` are
divided by that conformation's own normalizer, `N_A` or `N_I` ([Normalization](@ref)),
which the combination does not see, so every term of `D_A` and `D_I` is dimensionless.
`F_A` and `F_I` print only when they differ. They are 1 whenever the catalytic graph is
one rapid-equilibrium segment without a covalent form, the common case, and the law is then
the textbook MWC form. When the `:OnlyA` groups stop every catalytic cycle of the inactive
conformation, as the enumerated tag sets do, the law has no `Q_I`: the inactive branch adds
enzyme to the denominator and no flux. When the inactive conformation is free E alone, the
law has no `D_I` either, and its branch of the denominator is `L` times `F_A^n` and `W_I`,
each left out when it is 1. With `coefficients = true` each polynomial reads in named
coefficients (`dA_S` for the coefficient of `S` in `D_A`, `nI_0` for the constant term of
`Q_I`).

In the printed equation above `n = 2`, so each `D` is squared in the denominator, the
activator's site gives `W_A = (1 + A / K_A_Areg)^2`, and the inhibitor's
`W_I = (1 + I / K_I_Ireg)^2`.

The term cap counts both conformations: a law whose two conformations hold more than 5,000
spanning-tree products together is refused, and identify's `eq_complexity_filter` compares
the same sum.
````

In section "Thermodynamic constraints of MWC equations", in the paragraph "Often a K-system
uses exclusive `:OnlyA` binding: … the inactive cycle is pruned (`N_I = 0`).", replace
"(`N_I = 0`)" with "(no `Q_I`)".

- [ ] **Step 5: `docs/src/deriving/iso.md` and `docs/src/deriving/cha_king_altman.md`**

(a) In `iso.md`, replace "an `S·P` cross term in the denominator (the `… P * S /
(K_ES_to_E_S * K_EprimeP_to_Eprime_P)` terms above)." with "an `S·P` cross term in the
denominator, the term of `D` above that holds both `S` and `P`."

(b) In `cha_king_altman.md`, section "Overview of the algorithm", the segment-root
paragraph describes the rule of `_segment_root`, which Task 13 deleted. Replace it, from
"**Within each RE segment**, the package picks the form with the fewest bound" through
"`1 + S / K_ES_to_E_S + P / K_EP_to_E_P` denominator in the final equation.", with:

```markdown
**Within each RE segment**, the first form in the mechanism's form order is the
**segment root**, so free E, the first form, roots its own segment. Every other form in
the segment gets an alpha factor, its abundance relative to the root. The law in lowest
terms does not depend on which form roots a segment. The readable
`1 + S / K_ES_to_E_S + P / K_EP_to_E_P` denominator comes from the printer, which divides
the law by the denominator's constant term, free E's weight ([Normalization](@ref)).
```

Replace the next paragraph, "**Across segments**, the SS steps form an inter-segment rate
matrix. The Cha numerator and denominator come from cofactor determinants of that matrix,
with each diagonal entry weighted by the sigma sum of alpha factors for its segment.",
with:

```markdown
**Across segments**, the SS steps form an inter-segment rate matrix. The denominator
sums, over segments, each segment's sigma (the sum of its forms' alpha factors) times its
cofactor of that matrix, the sum of the spanning trees of the segment graph directed
toward it. The numerator sums each SS step's forward flux times the trees toward its
source segment, less its reverse flux times the trees toward its target, weighted by the
copies of the first substrate the step moves.
```

In section "Keeping the equation finite at zero concentration", after the paragraph that
begins "The derived rate equation is **division-free in the concentrations**", add: "The
printed law of a mechanism in the ratio form ([Normalization](@ref)) writes terms such as
`P / B`, but the code `rate_equation` runs is the law in lowest terms, free of such
divisions."

(c) In `docs/src/deriving/thermodynamic_constraints.md`, section "Haldane and Wegscheider
relations", replace the code line
`k_EP_to_ES = (1 / Keq) * K_EP_to_E_P * (1 / K_ES_to_E_S) * k_ES_to_EP` with the reference
line the printer writes, `k_EP_to_ES = K_EP_to_E_P * k_ES_to_EP / (K_ES_to_E_S * Keq)` (the
docstring's doctest, Task 12 Step 4), and "its line reads `K_E_to_F = 1.0`" with "its line
reads `K_E_to_F = 1`" (Task 12's pin 8).

- [ ] **Step 6: The identify pages**

(a) `docs/src/identify/tutorial.md`, section "Read the result": replace the paragraph from
"On noiseless data the search recovers the generating MWC mechanism." through "from `Keq`
and the independent constants." with:

```markdown
On noiseless data the search recovers the generating MWC mechanism. Its law ends
`v = E_total * (Q_A * D_A * W_A) * (S - P / Keq) / (D_A^2 * W_A + L * D_I^2)`: each
conformation's denominator is squared because the enzyme is a dimer, the activator's site
`W_A = (1 + A / K_A_Areg)^2` weighs the active conformation alone, and there is no `Q_I`,
since the inactive conformation does not catalyze. The Haldane relation, printed as a
reference line, fixes the dependent reverse rate `k_A_EP_to_ES` from `Keq` and the
independent constants.
```

After the docs build (Step 10), check this `v` line against the page the build renders:
`grep -o "v = E_total \* (Q_A[^<]*" docs/build/identify/tutorial/index.html` prints the line
quoted above. If the rendered line differs, quote the rendered one and fix the sentence
that reads it.

In the same section, replace "Every row names its mechanism by its text (`mechanism`) and a
short id (`mechanism_id`), and a child names its parent's id (`parent_id`)." with:

```markdown
Every row names its mechanism by its text (`mechanism`) and a short id (`mechanism_id`),
and a child names its parent's id (`parent_id`). Each row also carries its law's key
(`key`), which mechanisms whose laws agree up to renaming the fitted constants share: the
search fits one mechanism per key, and the others read that fit under their own names
(`fit_inherited`, `fit_source_id`). Read the CSVs with the ids as text,
`CSV.read(file, DataFrame; types = Dict(:mechanism_id => String, :parent_id => String,
:fit_source_id => String))`; [Best mechanism selection](@ref) lists every column.
```

In section "Run the search", replace "and `eq_complexity_filter` caps it by rate-equation
complexity — roughly the number of terms in the equation's denominator," with "and
`eq_complexity_filter` caps it by rate-equation complexity — roughly the number of terms in
the equation's denominator, both conformations together for an allosteric mechanism,".

(b) `docs/src/identify/model_selection.md`:
- Before the paragraph that begins "The search walks parameter counts in ascending
  order.", insert the paragraph:

  ```markdown
  The beam ranks laws, not mechanisms. Mechanisms whose laws agree up to renaming the
  fitted constants share a key, count once at the key's lowest loss, and share one fit,
  which each reads under its own names. Every member of a key the beam selects expands; a
  key the beam passes over stays out, members met later included.
  ```

- In the bullet "**Width floor**", replace "until `min_beam_width` of them have been
  expanded" with "until `min_beam_width` distinct keys have been expanded".
- In the keyword table, replace the `min_beam_width` row's meaning "Cumulative floor:
  expand at least this many mechanisms per parameter count over the whole search." with
  "Cumulative floor: expand at least this many distinct laws (keys) per parameter count
  over the whole search."
- Replace "The complexity bound roughly counts the terms in a rate equation's denominator;
  its default (`337`) passes a fully steady-state random-order bi-bi and skips anything
  denser." with "The complexity bound roughly counts the terms in a rate equation's
  denominator, both conformations together for an allosteric mechanism; its default
  (`337`) passes a fully steady-state random-order bi-bi and skips anything denser, and a
  one-subunit allosteric variant of a law near it can exceed it."
- In the `cv_results` column table: in the `mechanism_id` row replace
  `types = Dict(:mechanism_id => String, :parent_id => String)` with
  `types = Dict(:mechanism_id => String, :parent_id => String, :fit_source_id => String)`;
  replace the `rate_equation` row's description with "The printed law
  (`rate_equation_string`), in full."; replace the `eq_hash` and `fit_inherited` rows with:

  ```markdown
  | `key` | The SHA-256 (64 hex digits) of the law's canonical form: its terms with the fitted constants labelled canonically and the concentrations by name. Mechanisms share a key exactly when their laws agree up to renaming the fitted constants. A law whose canonical search reaches 100 labellings takes its mechanism's own key, the SHA-256 of `mechanism|` and its `mechanism_id`. `missing` in a failure row. |
  | `fit_inherited` | `true` when the row reads another mechanism's fit: one with the same key, under the row's own names, or, for an L phantom, the fit of the law it reduces to; `false` for the mechanism actually fit. |
  | `fit_source_id` | `mechanism_id` of the mechanism whose fit the row reads, the row's own when it was fit. |
  | `L_phantom` | `true` for a one-subunit allosteric variant whose inactive conformation binds nothing (free E alone, every regulator `:OnlyA`) and whose `L` an exact rescaling of its active constants absorbs: the data cannot determine `L`, so the row reads the fit of its active conformation's law, with `L = 0`, and keeps its own `n_params`. |
  ```

- Replace "Distinctness is by `eq_hash`." with "Distinctness is by `key`: LOOCV fits each
  key's source once per fold, and every candidate with that key, an L phantom of it
  included, reads those fold losses."

(c) `docs/src/identify/enumeration_engine.md`:
- In the allosteric section, after the list that follows "On its bases the engine emits,
  per multiplicity:", insert the sentence "Tag sets whose laws agree up to renaming the
  constants share one key, so identify fits the law once." at the start of the paragraph
  that follows the list.
- In the paragraph on K-type variants over one subunit, replace its last sentence ("The
  engine emits the variants where `L` is a phantom all the same, so a fit can carry an `L`
  the data cannot determine.") with: "The engine emits the variants where
  `L` is a phantom all the same, since their children can carry laws the parent's children
  lack. A variant whose inactive conformation binds nothing and whose `L` an exact
  rescaling of its active constants absorbs is an L phantom: when a mechanism with its
  active conformation's law is fitted, identify does not fit the phantom but reads that
  fit with `L = 0` and flags the row (`L_phantom` in the CSVs). Any other variant, a
  phantom whose active law nobody fits included, is fitted like any other candidate."
- Replace "A few children whose equation is the parent's up to renaming the constants
  survive; they cost one fit each and never win selection." with "A few children whose law
  is the parent's up to renaming the constants survive; they share the parent's key, so
  identify fits their law once, and they never win selection."

- [ ] **Step 7: `docs/src/developer.md`**

(a) "Derivation architecture", second paragraph: replace "The printed equation of
`rate_equation_string` keeps flat sums." with "`rate_equation_string` (`src/print.jl`)
prints the law rather than the code: each denominator term over the unit rule's normalizer,
dimensionless."

(b) "The joint thermodynamic solve": replace "`_dependent_param_exprs` (`src/thermo.jl`)
builds one integer system per mechanism" with "`_solve` (`src/thermo.jl`) builds one
integer system per mechanism"; after "Pivot columns are dependents." insert "`_solve`
returns them with the independents and the fitted names as a `_Solve`, which the engine,
`all_parameters`, code generation and the printer read; `_dependent_param_exprs` renders it
as expressions."; and replace the sentence "The engine derives each conformation's
polynomials from the mechanism, a mask of the edges its law keeps (`_state_mask`; the
inactive conformation's is `_inactive_mask(am)`) and the names of their constants
(`_edge_names`)." with "The engine (`_law`, `src/engine.jl`) derives each conformation's
law on a mask of the edges it keeps (`_active_mask(am)`, `_inactive_mask(am)`; every edge
of a plain mechanism)."

(c) Insert, after that section and before "## Well-posedness guard":

```markdown
### The law, its code, its text and its key

`src/poly.jl` holds sparse polynomials over integer-indexed variables (`_Poly`), with
rational exponents for constants and checked `Int64` coefficients: under the term cap no
coefficient leaves `Int64`, so an overflow is an internal error. `_law(m, s, mask)`
(`src/engine.jl`) derives one conformation's law, `_Law`, as `v = E_total·H·Q/D`. It
contracts each rapid-equilibrium segment to a node rooted at its first form (free E roots
its own), sums the spanning trees of the segment graph, builds the numerator from each
steady-state edge's net flux, substitutes the solve, divides out `H = Π[S] − Π[P]/Keq`
exactly and brings `Q` and `D` to lowest terms in the concentrations. An MWC law
(`_AlloLaw`) holds both conformations' laws, their free-E weights over a common factor and
the regulatory sites' polynomials. `_eq_complexity` counts V×τ, both conformations of an
MWC law together, by a checked Bareiss determinant; the cap (`_assert_derivable`) and
identify's `eq_complexity_filter` read it before any polynomial is built. `_Poly`, `_Law`
and `_AlloLaw` are mutable structs with `const` fields, never an immutable wrapper of one
reference: Julia 1.12 and 1.13 can hoist such a wrapper's box out of a loop that merges it
with another type, without a write barrier.

`_rate_body` (`src/codegen.jl`) emits the law as straight-line locals: the dependents in
solve-column order, each segment's weight and tree sums, `Q`, `H` and the MWC terms, all in
greedy Horner form with binary `+` and `*`. Each input is read by name into a positional
local, so a `NamedTuple` in any field order works, and reciprocals and roots (`sqrt`, or
`x^(one(x)/q)`) are computed once.

`rate_equation_string` divides the law by the unit rule's normalizer (`_normalizer`,
`src/normalize.jl`: the constant term, else a virtual one on the first affinely
independent terms, else the ratio form), so every denominator term is dimensionless, and
with `coefficients = true` prints it in named coefficients.

`_canonical` (`src/key.jl`) labels a law canonically: each part to lowest terms, the fitted
constants coloured by refinement on the term–constant incidence and individualized while
colours tie, the concentrations by name. The SHA-256 of its text is the law's key, by which
identify fits one mechanism per law and reads that fit through each mechanism's own names
(the canonical order). At 100 labellings the search stops, and the law takes its
mechanism's own key. The L-phantom certificate (`_phantom_weights`) solves exactly for a
rescaling of the active constants that absorbs `L` in a one-subunit variant whose inactive
conformation binds nothing (free E alone, every regulator `:OnlyA`); such a variant reads
the fit of its active law when a mechanism with that law is fitted.
```

(d) "Well-posedness guard": replace "Class 5, a denominator that keeps no monomial on a
face, runs on the derived law in `_num_den_exprs`." with "Class 5, a denominator that
keeps no monomial on a face, runs on the derived law (`_assert_support`) wherever a law is
derived for use: the generated `rate_equation`, identify's PASS1 and
`rate_equation_string`."

(e) "Enumeration engine architecture": replace "Duplicate equations that survive these
proofs are collapsed at compile time by `eq_hash`. A numerical identifiability rank exists
only in the test suite, as an oracle for the proofs; nothing in `src/` estimates
identifiability numerically." with "Duplicate laws that survive these proofs share one fit:
identify fits one mechanism per key (`_canonical`, `src/key.jl`). A numerical
identifiability rank exists only in the test suite, as an oracle for the proofs and the
L-phantom certificate; nothing in `src/` estimates identifiability numerically."

- [ ] **Step 8: `README.md` and `docs/src/getting_started.md`**

In `README.md`, after the line

```julia
print(rate_equation_string(m))             # the symbolic rate equation
```

add

```julia
print(rate_equation_string(m; coefficients = true))  # the law in named coefficients
```

In `docs/src/getting_started.md`, section "Derive the rate equation", after the block

````markdown
```@example getting-started
print(rate_equation_string(m))
```
````

add

````markdown
```@example getting-started
print(rate_equation_string(m; coefficients = true))
```
````

and at the start of the paragraph that follows, insert "`rate_equation_string` prints the
law in its independent constants, each denominator term relative to free enzyme;
`coefficients = true` prints it in named coefficients, each given by a legend line in the
fitted constants."

- [ ] **Step 9: `ROADMAP.md`**

(a) In "## Planned", delete the entry "**Rewrite the derivation around the reaction
graph.**" (from its `- **` line through its "Priority: high." line): this PR is its last
phase.

(b) Replace the entry "**Drop `change_allo_state` relaxations that free no parameter, if
measurement finds any.**" (its bullet through "Priority: medium-low.") with:

```markdown
- **Drop `change_allo_state` relaxations that free no parameter, if measurement finds any.**
  `_expand_change_allo_state` (`src/mechanism_enumeration.jl`) has no gain test, unlike
  the split and RE→SS moves, and in a past LDH run it produced 20 of 33 zero-gain edges.
  A child whose law equals one already fitted, up to renaming its constants, shares that
  fit through the key (`_canonical`, `src/key.jl`) and costs no fit; a zero-gain child
  with a law of its own still costs one. First count, over the bi-bi allosteric depth-2
  population, the children whose `_independent_param_count` does not exceed their
  parent's and whose key is new, and confirm them with `_testhelper_identifiable_rank`. A
  naive filter can drop a different model with the same count, so check reachability
  before filtering. Priority: medium-low.
```

(c) Consistency pass over the entries Tasks 3, 4, 8, 10 and 13 changed:
`grep -n "eq_hash\|_rate_eq_dedup_key\|sym_poly\|_num_den_exprs\|POLY\|[Pp]hase 5\|reaction-graph rewrite\|the rewrite" ROADMAP.md`
prints no line that describes phase 5 as future work or names a deleted function; rewrite
any it prints in the present tense against the code as it stands. `docs/src/roadmap.md`
lists no graph-rewrite direction (`grep -n "rewrite\|eq_hash" docs/src/roadmap.md` prints
nothing), so it stays.

- [ ] **Step 10: Build the docs**

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T15_docs.log 3600 --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path = pwd())); Pkg.instantiate(); include("docs/make.jl")'`
Expected: no `Error` line, and no `Warning` but `deploydocs`' "Skipping deployment."; the
doctests of `rate_equation_string` and textbooks.md pass; `[Normalization](@ref)` resolves
to textbooks.md's heading (the only `## Normalization` heading under `docs/src`); then
Step 6(a)'s check of the rendered tutorial line.

- [ ] **Step 11: Run the full suite and commit**

Run (background): `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && $S/phase3/jl.sh $S/phase45/logs/P5_T15_full.log 3600 --project=. -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`, polled.
Expected: `exit 0`, every test passing, pristine; record the time, test count and
`maxrss_kb` beside Task 0's.

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt
git add docs/src/deriving/textbooks.md docs/src/deriving/ping_pong.md \
    docs/src/deriving/mwc_allostery.md docs/src/deriving/iso.md \
    docs/src/deriving/cha_king_altman.md docs/src/deriving/thermodynamic_constraints.md \
    docs/src/identify/tutorial.md docs/src/identify/model_selection.md \
    docs/src/identify/enumeration_engine.md docs/src/developer.md README.md \
    docs/src/getting_started.md ROADMAP.md
git commit -m "$(cat <<'EOF'
Document the printed law, the key and the inherited fits

The derivation pages describe the unit rule's three cases, the N line and
the coefficient form, Cleland's normalizer as the ping-pong case of the
virtual constant term, the MWC law's per-conformation Q, D, F, W and N,
and the engine's segment root at a segment's first form; the identify
pages the key, the beam and LOOCV per key, the inherited fits and L
phantoms, the CSV columns key, fit_inherited, fit_source_id and L_phantom,
and reading the ids as text; the developer page the solve
struct, the engine, the code generation, the printer, the key and the
certificate. ROADMAP retires the graph rewrite, whose last phase this is,
and notes that same-law change_allo_state children now share a fit.

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE
EOF
)"
```

- [ ] **Step 12: Run the GC gate on the final head**

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; G=$S/phase45/gcscan5; LG=$S/phase45/logs; cd $S/phase45/law_wt
export JULIA_LLVM_ARGS="-pass-remarks=julia-licm -pass-remarks-analysis=julia-licm"
for v in 112 113; do ch=$([ $v = 113 ] && echo +1.13)
  for w in poly "poly control" derive derive_control identify; do n=$(echo $w | tr ' ' '_')
    $S/phase3/jl.sh $LG/P5_T15_scan${v}_$n.log 1800 $ch --project=$G/env$v --startup-file=no --pkgimages=no --heap-size-hint=2500M $G/scan_law.jl $w
  done
  $S/phase3/jl.sh $LG/P5_T15_repro$v.log 600 $ch --startup-file=no --pkgimages=no $S/phase3/gcprobe/repro_upstream.jl
done
unset JULIA_LLVM_ARGS
python3 -I $S/phase3/gcprobe/hoistscan.py $LG/P5_T15_scan11*.log $LG/P5_T15_repro11*.log > $LG/P5_T15_hoistscan.log 2>&1
cat $LG/P5_T15_hoistscan.log
```

Run it in the background and poll; each scan recompiles the package with the remarks on
(about 2 min per version per workload). Expected, on both versions: `poly`, `derive` and
`identify` print three rounds `same … full gc ok` and report `hoisted with a bits store
into a reference field: 0`; `poly_control`, `derive_control` and `repro` each report at
least one `DANGER`, the hazard the gate must see. A package `DANGER` line stops the task:
find the loop, give it a concrete type, rescan. A control without `DANGER` means the scan
cannot see the hazard: stop and report the logs.

- [ ] **Step 13: Run the phase-5 test files on Julia 1.10.12 (P5-38)**

`law_wt`'s Manifest.toml was resolved on Julia 1.12.5 with the newest SciML, so 1.10 gets
an environment of its own, built once as `gcscan5/mkenv.jl` builds the scan environments.
Create `$S/phase45/scripts/P5_mkenv110.jl`:

```julia
# ABOUTME: Creates a Julia 1.10 environment that develops the phase-5 worktree by path and
# ABOUTME: holds TestEnv, then checks that the package and its test environment load.
# Usage: julia +1.10 P5_mkenv110.jl <worktree> <env_dir>
using Pkg
Pkg.activate(ARGS[2])
Pkg.develop(path = ARGS[1])
Pkg.add("TestEnv")
Pkg.instantiate()
using TestEnv
TestEnv.activate("EnzymeRates")
using Test, EnzymeRates
println("ok ", VERSION, " ", pkgversion(EnzymeRates), " from ", pathof(EnzymeRates))
```

Run: `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; $S/phase3/jl.sh $S/phase45/logs/P5_T15_v110_env.log 1800 +1.10 --heap-size-hint=2500M $S/phase45/scripts/P5_mkenv110.jl $S/phase45/law_wt $S/phase45/env110`
Expected: `ok 1.10.12 0.10.0 from …/law_wt/src/EnzymeRates.jl`, `exit 0`. Pkg resolves the
test dependencies for 1.10 and may fetch versions the 1.12 depot lacks. A resolution or
precompile failure here is the environment's, not P5-38's: stop and report the log; the
PR body then lists the 1.10 run as not done, for Denis.

Then, one at a time, in the background, the phase-5 test files, `test_oracles.jl` (Q1–Q6)
and `test_identify_rate_equation.jl` (the key cache, the beam per key, the phantoms) among
them:

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && for f in test_poly:poly test_law:law test_key:key test_print:print test_rate_eq_derivation:tr test_oracles:oracles test_identify_rate_equation:tid; do $S/phase3/jl.sh $S/phase45/logs/P5_T15_v110_${f#*:}.log 3600 +1.10 --project=$S/phase45/env110 --heap-size-hint=2500M -e 'using TestEnv; TestEnv.activate("EnzymeRates"); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/'${f%%:*}'.jl")'; done
```

polled one log at a time. Expected: each `exit 0` with every test passing; the unicode
coefficient names (`d_A²_B`, `n_B⁻¹`) parse on 1.10 (test_print.jl); the `rate_equation`
gates of test_rate_eq_derivation.jl hold (0 allocations, under 120 ns, under `<c>` ns per
term). The full suite does not run on 1.10 here: it can peak above 4 GB (ROADMAP known
issue). A file that fails only on 1.10 stops the task: report the failure (P5-38 names the
slice's untyped empty `map`, which explicit element types remove).

- [ ] **Step 14: Write the final PR body**

Replace `$S/phase45/pr_body5.md` with the text below, each `<…>` the figure of the log or
file named beside it (Task 0's `p5_base/figures.md`, the task logs named in each line).
The open-questions slot takes, one line each, every item of this plan's "Open for Denis"
→ "Unresolved items from the checkers" that the tasks did not settle, with its task. The
body's other sections carry the rest of "Open for Denis": the rulings, the spec and
CLAUDE.md edits, the deviations from the research, the deleted test subjects and the suite
time. `<T6 census result>` is `PASS`, or the count and keys of the laws the summed cap
refuses (`P5_T6_census_refused.log`), the only differences.

````markdown
## Summary

Phase 5 of the reaction-graph rewrite (spec `docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md`,
sections 9, 10, 11 and the phase-5 parts of 13 and 14). One exact engine derives every rate
law on phase 3's graph and phase 4's masks; `rate_equation` is generated from it,
`rate_equation_string` prints it, and identify fits one mechanism per law.

- **Solve** (`src/thermo.jl`): `_solve` returns the joint solve as `_Solve`;
  `all_parameters` reads it, and `_evaluate_dependent` is gone.
- **Engine** (`src/poly.jl`, `src/engine.jl`): one conformation's law on an edge mask in
  checked `Int64`, `H` divided out exactly, lowest terms; the MWC law on `_active_mask` and
  `_inactive_mask`; V×τ by a checked Bareiss determinant; the 5,000 cap and the
  complexity filter on both conformations together.
- **Code generation** (`src/codegen.jl`): locals for the dependents in solve-column order,
  the segment sums, `Q`, `H` and the MWC terms, in greedy Horner form; inputs read by name;
  roots computed once. 0 allocations; specs at most <t_specs> ns; the six gate laws under
  `max(120 ns, <c> ns × D)` (Task 14).
- **Unit rule and printer** (`src/normalize.jl`, `src/print.jl`): every denominator term
  dimensionless (constant term, else the graded virtual one, else the ratio form);
  reference lines as ratios of products with exact exponents; `N = …`; MWC `Q_s`, `D_s`,
  `F_s`, `W_s`, `N_s`; `rate_equation_string(m; coefficients = true)`.
- **Key** (`src/key.jl`): a rename-proof canonical labelling, its SHA-256 in 64 hex digits;
  at 100 labellings a mechanism's own key.
- **identify**: one fit per key, read by each row under its own names; `_offer_cv!`, the
  beam and the LOOCV folds per key; L phantoms (certificate inside the spec's pre-filter)
  inherit the fit of whoever holds their L→0 key, else they are fitted unflagged; PASS1
  derives, keys and prints the concrete law and compiles nothing. CSV: `key`,
  `fit_inherited`, `fit_source_id`, `L_phantom`; `eq_hash` is gone.
- **Deleted:** the old engine and printer (`POLY`, `sym_det`, `_raw_symbolic_rate_polys`,
  `_num_den_exprs`, `_state_rate_polys`, the MWC `Expr` assembly, `_segment_root`,
  `_equation_text`, `_expr_to_string`), `_rate_eq_dedup_key` and `eq_hash`,
  `_evaluate_dependent`.
- Version stays 0.10.0 (R19, P5-42). Draft, stacked on the phase-4 PR (#<phase-4 PR>);
  not for merge.

## For Denis first

Two questions the rulings leave open, then the four calls the research proposed:
- **P5-10a: restore Cleland-first?** The printed `N` follows the spec's graded basis. It
  differs from Cleland's `N` in 3 of 5 eligible specs and in 257 of 292 eligible ter-ter
  conformations (13.8% of ter-ter laws); values never change; Cleland-first is 4 lines.
  Saying yes re-pins the printed laws where the two differ.
- **P5-31a: drop the L-phantom pre-filter?** The certificate runs inside the spec's
  pre-filter (one subunit, inactive conformation free E alone, every site ligand
  `:OnlyA`). Without the pre-filter, 17–21% of bi-bi fits go instead of 2.5–3.2%, each
  with an exact proof; this branch flags <flagged> variants in Task 10's rerun, and the
  certificate alone would flag <certificate-only> more (P5_T10_rerun.log).
- **P5-14** the coefficient form sits behind `rate_equation_string(m; coefficients =
  true)`; the default and identify's CSV keep the full law (the two are about as long).
- **P5-20** the `key` column holds all 64 hex digits; `mechanism_id`, `parent_id` and
  `fit_source_id` stay 16, read as text with the documented `types`.
- **P5-26** the complexity filter compares both conformations' V×τ, summed, with 337.
  Recount on P4-1's tagging (P5_T4_recount.log): <per-round figures>.
- **P5-32** a phantom inherits the fit of whoever holds its L→0 key; when nobody does, it
  is fitted as an ordinary candidate, unflagged.

## Spec and CLAUDE.md edits for approval

- **P5-24 (CLAUDE.md)** the size gate's constant: 0.8 → <c> ns per denominator term,
  twice the worst of two CI runs (worst <w> ns per term on <runner>), rounded up to one
  figure; a third run on the new value passed (P5_T14_c.log, P5_T14_ci3.log). The same
  paragraph's headroom sentence now reads the slowest specs, <spec 1> and <spec 2>, at
  about <t> ns (<k>× under 120 ns; P5_T6_hoist.log), where it read about 55 ns and 2×.
- **P5-1 (spec §9)** "checked `Int64`, retried in `BigInt` on overflow" → "checked `Int64`;
  overflow is an internal error". Proof note: a coefficient is at most max(F·τ, 2τ·Σ|w_k|)
  with τ ≤ 5,000; the largest seen is 36 over 2,642 laws; code generation's `fI·DA`, `fA·DI`
  and the certificate's products stay below τ_I·F·τ_A ≈ 4e8, and the cap is checked before
  any product (P5-2).
- **P5-25 (CLAUDE.md names the helper)** the gate's helper times `rate_equation` through a
  `@noinline` barrier: on package code the old helper timed <h> of 42 specs under 0.5 ns
  (P5_T6_hoist.log)<; with none, the helper stayed as it was>.
- **P5-30 (spec §13)** "an exact rank test in `src/`" → "an exact rescaling certificate in
  `src/` (`key.jl`), inside the pre-filter"; ranks are test oracles.
- **P5-32 (spec §10, §13)** a phantom inherits "the parent's fit" → "the fit of the
  mechanism that holds the key of its law at L = 0 (often the parent)", its loss and
  fold losses that source's; §13 adds that a phantom with no such source is fitted
  unflagged.
- **P5-37 (ROADMAP)** the Idea "Apply the chain-flank exclusion at one catalytic subunit",
  which would remove the never-fitted chain-flank phantoms.
- **P5-41** phase 5's population tests add <s5> s (Q1–Q6) and the second-process check
  <s7> s; with phases 2 and 4, about <total> s against the spec's one minute.
- **P5-43a** may the research logs and data (`/home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase45/research5/`,
  `phase3/gcprobe/`) go to the reference branch? Nothing is pushed without your OK.

## Deviations from the research, for approval

- **Overflow (P5-1, RES §1.3)** `_law` names the mechanism on an `OverflowError`, plain and
  MWC (the MWC law's growing arithmetic runs inside its two plain laws); code generation's
  products, the normalizer and the key let a raw `OverflowError` escape, which the cap
  rules out (products below τ_I·F·τ_A ≈ 4e8, P5-2).
- **BigFloat through the body (RES §3.1)** on the 11 allosteric specs, four plain ones and
  the cube-root fixture, not at every oracle call: each check compiles another
  `rate_equation` method (<s6> s for the 15).
- **Task 9's identify gate (RES §4 row 9)** asks for the truth's key only when Task 0's run
  selected the truth; otherwise Task 0's selection, or a near-exact law with at most its
  constants, also passes, since on these noise-free data the choice among near-exact laws
  follows the optimizer's random stream. Task 0 selected the truth: <yes or no>; Task 9:
  <yes or no>; constants <n0> and <n9>. The data put each product at 10–100% of its
  substrate, so the releases' product terms show.
- **RES §4 row 12** runs as Tasks 12 (the printer) and 12b (PASS1 on the concrete law), so
  PASS1's compile gate cannot hold back the printer.

## Rulings P5-0 to P5-46

- P5-0 P5-1 to P5-44 adopted as the research writes them, except P5-10a, P5-31a, P5-43a and
  P5-45.
- P5-1 ◆ checked `Int64`, no `BigInt` retry; `_Poly`, `_Law`, `_AlloLaw` mutable with
  `const` fields; a flip-only inactive conformation is a law with Q = 0, D = 1.
- P5-2 the overflow note covers code generation's and the certificate's products; the cap
  is checked before any product.
- P5-3 `_eq_complexity` returns `Int` (checked Bareiss; the `BigInt` determinant after an
  overflow, saturated); an MWC law counts active plus inactive.
- P5-4 `_solve` returns `_Solve`; `_dependent_param_exprs` renders it.
- P5-5 RE turnover stays class 2; the engine asserts zero net uptake on every non-tree RE
  edge.
- P5-6 class 5 reads the new law.
- P5-7 the GC gate is a typed workload outside the suite with positive controls, on 1.12.5
  and 1.13.1.
- P5-8 every MWC use reads `_active_mask` and `_inactive_mask`.
- P5-9 the unit rule per conformation.
- P5-10 ★ Cleland-first proposed; P5-10a: declined until Denis answers.
- P5-11 term order: degree, then non-substrates, then the exponent vector, `A²` before
  `A·B`.
- P5-12 fractional λ; roots print as `sqrt` and `^(p/q)`.
- P5-13 the ratio form (`d_P_B⁻¹`, `n_B⁻¹`; no `n_0` without an origin term).
- P5-14 ★ the coefficient form behind the keyword.
- P5-15 MWC names `Q_s`, `D_s`, `F_s`, `W_s`, `N_s`; `H` never prints as a name.
- P5-16 `v = E_total * (Q) * (H) / (D)`.
- P5-17 reference lines as ratios of products with exact exponents; `= 1` for a pinned
  constant.
- P5-18 the key from the concrete law in PASS1, concentrations by name.
- P5-19 leaf cap 100; the mechanism's own key at the cap, logged.
- P5-20 ★ the 64-hex key; 16-hex ids.
- P5-21 a second process checks bodies, keys and printed laws; five spec keys pinned;
  dependent locals by solve column.
- P5-22 roots: `sqrt` for q = 2, `x^(one(x)/q)` above.
- P5-23 the size gate counts `D` after substitution (MWC: `D_A + D_I`).
- P5-24 ◆ c by phase 1's rule on CI.
- P5-25 ◆ the gate's helper behind a `@noinline` barrier.
- P5-26 ★ the summed complexity filter at 337.
- P5-27 fits stored by canonical position, read through each row's order; `fit_source_id`;
  the source's `peak_at_range_edge`.
- P5-28 `_offer_cv!` per key; LOOCV per fit key.
- P5-29 the beam per key; selected and rejected keys persist across sweeps.
- P5-30 ◆ the certificate in `key.jl`; ranks are oracles.
- P5-31 no pre-filter proposed; P5-31a: the pre-filter stays until Denis answers.
- P5-32 ★ a phantom inherits whoever holds its L→0 key; otherwise fitted unflagged (b2).
- P5-33 a phantom's names come from its active law's canonical order; L = 0.
- P5-34 a phantom takes its source's fold losses.
- P5-35 the oracle rank in absolute mode.
- P5-36 distinct counts reported; the ping-pong rerun below.
- P5-37 ◆ the chain-flank Idea.
- P5-38 Julia 1.10 run; the site list has an explicit element type.
- P5-39 shuffled NamedTuples and a metabolite named `x1` through the body.
- P5-40 the ROADMAP entries this phase settles.
- P5-41 ◆ population tests within 20 s.
- P5-42 version 0.10.0.
- P5-43 ◆ the research archive; P5-43a: deferred.
- P5-44 the evidence limits below.
- P5-45 this plan's method; branch graph-law on graph-guard; P5-48 its PR targets main.
- P5-46 the refresh on f9361cf, where graph-law is cut, and Task 0's baseline; Task 0
  brings graph-law to graph-guard's head 991355f (draft PR #94), which changes no `src`
  file.

## Deleted test subjects, for approval

- test_accessors: `_evaluate_dependent`'s `sqrt(Keq)` error assertion (Task 1); the
  testset now checks `all_parameters`'s rational powers.
- TID "_rate_eq_dedup_key" (Task 8): converted to the key's determinism.
- TR "rendering helpers" (Task 13): `_invert_monomial`, `_is_metabolite_free_monomial`,
  `_mwc_cross_weight`; test_law.jl's `_cross_weights` testset carries its three cases.
- TR "a segment's root ties go to the canonical form order" (Task 13): `_segment_root`'s
  rule (fewest bound metabolites, then no residual) goes; the testset now pins the engine's
  rooting at a segment's first form.
- The old engine and printer, listed in the summary, and with them the `POLY` class-5
  failure path, whose last direct tests, test_admission.jl's two `_num_den_exprs` throw
  assertions, go with it in Task 13.

## L phantoms (P5-36)

The research's distinct counts: the spec's rule flags 738 bi-bi variants (188 in round 1,
550 in round 2), the certificate 3,102. Ping-pong coverage was zero there, since ping-pong
round 1 duplicated bi-bi round 1. Task 10 reran R5.2 on P4-1's ping-pong tagging children
(n = 1, round 1 and a seeded round-2 sample, P5_T10_rerun.log): <flagged> flagged,
<certificate-only> certificate-only, every flagged variant with key(active law) =
key(L→0 law), np = np₀ + 1, the map oracle holding and equal ranks.

## Evidence limits (P5-44)

- R5.1's speed gate was judged by a same-run ratio, since the 70.7 µs population was lost.
- The research's identify GC workload fitted pass 1 only, with a stand-in loss; this
  branch's gate runs a real identify workload (Tasks 6 and 15).
- Containment was sampled (43 of 926) and numeric.
- Ter-ter gave 1,865 laws of the 2,000 asked.

## Open questions from the plan

<open plan items>

## Measurements

- Gates (Task 14): CI ns per term, ter-bi / bi-bi, worst of three runs: ubuntu x64
  <a> / <b>, windows x64 <c1> / <d>, macOS aarch64 <e> / <f>; c = <c>; specs at most
  <t_specs> ns locally; 0 allocations.
- Compile budgets (P5_T12b_budget.log): init trace <n>; first `rate_equation` <s> s; Sig
  <i> instances, <kb> KB; fit per fresh mechanism <instances> instances, <kb> KB; PASS1
  per fresh mechanism 0 instances, 0 show instances.
- GC gate (P5_T15_hoistscan.log): 0 hoisted with a bits store on 1.12.5 and 1.13.1 for
  `_Poly`, the derivation and identify; the controls fire.
- Key cost (P5_T7_keycost.log): median <ms> ms, max <ms> ms over the 239 bi-bi and 200
  ter-ter laws; 0 leaf-cap hits.
- L-phantom detector (P5_T10_rerun.log; spec section 17): `_is_phantom` median <us> µs,
  max <us> µs per one-subunit variant, `_phantom_weights` median <us> µs, max <us> µs,
  over the ping-pong rerun's rounds 1 and 2 (the research: about 25 µs on slice code).
- R5.5 recount (P5_T4_recount.log): <figures>.
- identify run (Task 9 against Task 0): fits run <n9> against <n0>; selected the truth:
  Task 0 <yes or no>, Task 9 <yes or no>.
- Census: P5_T6 <T6 census result>; P5_T12
  children, values and names equal, <n> laws print anew; P5_T13 PASS.
- R52: full-suite maxrss, ter-ter init and TME focused, Task 0 → 6 → 12 → 13: <figures>.
- Full suite: <t0> → <t15> s, <n0> → <n15> tests.
- Suite time against RES §4's estimates: Task 3 +<s> s (RES +2), Task 6 +<s> s (+10; the
  BigFloat bodies <s6> s of it), Task 8 +<s> s (+2), Task 10 +<s> s (+7), Task 12 +<s> s
  (+2); P5-41's 20 s covers Q1–Q6 only.
- Julia 1.10.12, own environment: test_poly, test_law, test_key, test_print,
  test_rate_eq_derivation, test_oracles and test_identify_rate_equation pass.

## Version

0.10.0, unchanged (R19, P5-42).

🤖 Generated with [Claude Code](https://claude.com/claude-code)
````

Check: `grep -c "claude.ai/code/session" $S/phase45/pr_body5.md` prints `0`, and
`grep -nE "<[A-Za-z0-9_ -]+>" $S/phase45/pr_body5.md` prints no unfilled slot.

- [ ] **Step 15 (controller): Push and update the PR body**

```bash
S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad; cd $S/phase45/law_wt && git push && gh api repos/{owner}/{repo}/pulls/<n> -X PATCH -F body=@$S/phase45/pr_body5.md --jq '.html_url'
```

(`<n>` is Task 14's PR number.) Expected: the push succeeds and the PR's URL prints; the
PR page shows the new body. CI on the pushed head is green; never merge.

---

## Open for Denis

The PR body (Task 15 Step 14) names every item below. The rulings come from RES §4 "For
Denis" and RULINGS.md "Phase 5"; the rest are what the plan's checkers left unresolved.

### Rulings, ★ first

The two questions the rulings leave open lead:
- ★ **P5-10a** (on P5-10, Cleland-first): restore it? N stays the spec's graded basis until
  Denis answers, since he chose "lowest-degree terms + fixed order" on 2026-10-11. N differs
  from Cleland's in 3 of 5 eligible specs and 257 of 292 eligible ter-ter conformations
  (13.8% of ter-ter laws); values never change; Cleland-first is 4 lines. A yes re-pins the
  printed laws where the two differ.
- **P5-31a** (on P5-31, no pre-filter): drop the pre-filter? It stays
  (`!any(_inactive_mask(am))`, every site ligand `:OnlyA`), and the certificate runs inside
  it. Dropping it removes 17–21% of bi-bi fits instead of 2.5–3.2%, each with a proof; until
  then about 16% more fits run; the change is one line.

The other ★ calls:
- ★ **P5-14** the coefficient form sits behind `rate_equation_string(m; coefficients =
  true)`; the default and identify's CSV keep the full law (Task 12).
- ★ **P5-20** the `key` column holds all 64 hex digits; `mechanism_id`, `parent_id` and
  `fit_source_id` stay 16, read as text with the documented `types` (Task 8).
- ★ **P5-26** the complexity filter compares V×τ(active) + V×τ(inactive) with 337; Task 4's
  recount gives the newly filtered candidates per round.
- ★ **P5-32** a phantom inherits the fit of whoever holds its L→0 key; when nobody does, it
  is fitted as an ordinary candidate, unflagged (b2; Task 10).

◆ Spec, CLAUDE.md and ROADMAP edits:
- **P5-1** spec §9 "Engine": "checked `Int64`, retried in `BigInt`" becomes "checked
  `Int64`; overflow is an internal error", with the proof note of P5-2 (Task 2).
- **P5-24** CLAUDE.md's size-gate constant c, set from two CI runs, and the same
  paragraph's headroom sentence, which named the slowest specs at about 55 and 30 ns
  (Task 14, from Task 6's readings); the ROADMAP known issue "Ubuntu and windows CI miss
  `rate_equation`'s 120 ns gate on one spec" is <deleted | rewritten> (Task 14 says
  which).
- **P5-25** the gate's helper times through a `@noinline` barrier; CLAUDE.md names the
  helper (Task 6, when the hoisting probe finds hoists).
- **P5-30** spec §13: "an exact rank test in `src/`" becomes "an exact rescaling certificate
  in `src/` (`key.jl`), inside the pre-filter"; ranks are test oracles (Task 10).
- **P5-32** spec §10 and §13: a phantom inherits "the fit of the mechanism that holds the
  key of its law at L = 0 (often the parent)", not "the parent's fit", with that source's
  loss and fold losses; §13 adds that a phantom with no source is fitted unflagged
  (Task 10).
- **P5-31** see P5-31a.
- **P5-37** the ROADMAP Idea "Apply the chain-flank exclusion at one catalytic subunit"
  (Task 10).
- **P5-41** phase 5's population tests stay within 20 s; with phases 2 and 4 the total is
  about 80 s against the spec's one minute.
- **P5-43**, deferred by **P5-43a**: may the research logs and data
  (`/home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase45/research5/`,
  `phase3/gcprobe/`) go to the reference branch? Nothing is pushed without his OK.
- **P5-45** this plan's method: written against the read-only snapshot a9e0699 while
  phase 4's fix wave ran (skeleton, four writers with checkers, assembly, reviews, one fix
  pass); the refresh pass of P5-46 re-anchored every cite on f9361cf, where graph-law was
  cut, and Task 0 brought graph-law to graph-guard's head 991355f (draft PR #94), which
  changes no `src` file; branch graph-law, draft PR on main (P5-48), version 0.10.0.

### Deviations from the research

- **Overflow** (P5-1, RES §1.3, Global Constraints "Arithmetic"): `_law`, plain and MWC,
  turns an `OverflowError` into an internal error naming the mechanism; code generation,
  the normalizer and the key let a raw one escape, which P5-2's bound (τ_I·F·τ_A ≈ 4e8)
  rules out under the cap.
- **BigFloat through the body** (RES §3.1): Task 6 runs it on the 11 allosteric specs, four
  plain ones and the cube-root fixture, not at every oracle call site (about 80), since
  each check compiles another `rate_equation` method.
- **Task 9's identify gate** (RES §4 row 9): the truth's key is required only when Task 0's
  run selected the truth; otherwise Task 0's selection, or a near-exact law (`cv_score`
  below 1e-6) with at most its constants, also passes, and both runs' truth checks are
  reported. The synthetic data put products at 10–100% of their substrates, yet carry no
  noise: on f9361cf the run misses the truth, and the choice among near-exact laws
  follows the optimizer's random stream (item 34).
- **RES §4 row 12** runs as Tasks 12 and 12b.

### Deleted test subjects

- test_accessors: `_evaluate_dependent`'s `sqrt(Keq)` error assertion (Task 1); the testset
  now checks `all_parameters`'s rational powers.
- TID "_rate_eq_dedup_key" (Task 8), converted to the key's determinism test.
- TR "rendering helpers" (Task 13); test_law.jl's `_cross_weights` testset carries its three
  cases.
- TR "a segment's root ties go to the canonical form order" (Task 13), converted to "each
  rapid-equilibrium segment is rooted at its first form".
- The old engine and printer (Task 13), and with them the `POLY` class-5 failure path, whose
  last direct test goes in Task 5 (item 16 below).

### Suite time

Five tasks run over RES §4's estimates: Task 3 +3 s (RES +2), Task 6 +15 s (+10: the
BigFloat bodies add about 5 s), Task 8 +8 s (+2), Task 10 +12 s (+7), Task 12 +8 s (+2).
Each commit reports the measured figure, and the PR body's Measurements carry them; P5-41's
20 s covers Q1–Q6 only.

### Unresolved items from the checkers

"Resolved at assembly", "Resolved at review" and "Settled (P5-46)" (the refresh on
f9361cf) mark an item this file settles; every other item is open.

Tasks 0–3:
1. P5-8 says phase 4's 13 `:OnlyI` fixtures gate Tasks 3–4; RES §4 puts them in Task 4. Task 3 tests one (`competitive_RE`) on both masks against plain laws.
   All 13 per conformation would need a one-conformation oracle; the only one today is the
   old engine's `_state_rate_polys`, which Task 13 deletes. The full set stays in Task 4.
2. Task 4's red run expected six new top-level test_law.jl testsets to fail, though a
   focused run stops at the first. Resolved at assembly: that run now stops at the fixtures
   file's `LoadError`, and the Global Constraints' red-run rule fixes the same defect in the
   red runs of Tasks 4–12.
3. Task 3 Step 10's path edit to the `eq_complexity_filter` ROADMAP entry waited on FW-37.
   Settled (P5-46): FW-37 retitled it "`eq_complexity_filter` counts an allosteric
   mechanism's catalytic graph once." (ROADMAP.md:24 at f9361cf) and kept the
   `_eq_complexity` (`src/rate_eq_derivation.jl`) cite that the step rewrites (:25).
4. Task 0's identify run is unmeasured (estimated 10 min; 5,400 s timeout; at 90 min the
   controller lowers `n_restarts` for Tasks 0 and 9 alike). `P5_compare.jl`'s `laws`,
   `dependents` and `dterms` kinds were dry-run and `complexity` re-run by the checker;
   `strings` never was. Settled (P5-46): on f9361cf (Julia 1.12.5) the run as Task 0
   writes it takes 102 s (120 s with loading), runs 599 fits and selects a 7-constant
   one-subunit L phantom (`E(A) + B <--> E(P, Q)` merged, every step `:OnlyA`); no
   fitted mechanism has the truth's steps (`refresh5/settle_1.log`, `settle_14.log`).
   All five comparisons, `strings` included, pass on f9361cf's baseline
   (`refresh5/map-head-t0-3_2.log` to `_6.log`). Item 34 has what this means for Task 9.
5. The RE spanning tree is not the slice's `re_forest` BFS: each form joins by the first
   kept RE edge in edge order from a form already placed. The substituted law is the same,
   since the solve makes every RE cycle's constants multiply to 1, and the oracle,
   regrouping and baseline tests confirm it; it is not a line-for-line port. Denis may
   prefer `_segments`' BFS parents.
6. RES §3.2 and the first draft of Task 1 said `keq == 1//2`. Resolved at assembly: Task 1
   pins `keq == -1//2` (probe `p1_solve_facts.log`).
7. ROADMAP live-byte entry: Task 3 attributes the remaining 18 KiB figure to the plain
   random-order bi-bi. If it came from the allosteric form, Task 3 Step 10's last sentence
   needs rewording. Settled (P5-46): the figure is `_partition_independent_count`'s,
   whose one method takes a `Mechanism` (`src/thermo.jl:545`), so it is the plain form's
   and Step 10's sentence stands.
8. Task 3's overflow testset builds a malformed `_Solve` by hand, since no real law
   overflows under the cap. Does this count as a real-path test, not a mock, under
   CLAUDE.md?
9. The GC gate's toy control fired (10 hoists, GC abort) only in a stand-in run on Julia
   1.12.5; 1.13.1 is unverified, and Task 2 stops if the control does not fire there.
   Settled (P5-46): on 1.13.1 the same stand-in (`p1_scan_check.jl`) finds 10 of 65 with
   a bits store, at `range.jl:928`, and aborts the collector (exit 134); `_Poly` gives 0
   of 55 and three `same` rounds (`refresh5/settle_7.log`, `_8.log`, `_9_hoist.log`).
10. The MWC law and code generation's products (P5-2) do not turn an overflow into the
    internal error; only the plain `_law` wraps its arithmetic. Resolved at review: the
    Global Constraint now names the law (the MWC law's growing arithmetic runs inside its
    two plain laws) and cites P5-2's bound for code generation, the normalizer and the
    key; the narrowing is under "Deviations from the research".
11. test_law.jl pins `n_plain == 31`. Settled (P5-46): the fix wave changed only a
    comment in `MECHANISM_TEST_SPECS`, which at f9361cf holds 42 specs, 11 allosteric, so
    31 plain.

Tasks 4–7:
12. RES §3.1 asks `_testhelper_check_law_against_oracle` to feed BigFloat inputs through
    the compiled body at 1e-20. No task does it; Task 6 checks it on the cube-root fixture
    only. Adding it compiles a BigFloat `rate_equation` per checked mechanism (about 80
    call sites), likely tens of seconds against P5-41's 20 s. Resolved at review: Task 6
    adds a `bigfloat` keyword and runs it on 15 specs (about 5 s); the sample is under
    "Deviations from the research".
13. Q3 (Task 7) adds every tagging child of ten catalytic graphs of the one-round
    population. That pool's size and time were never measured (a graph can have about 99
    K-types); Task 7 Step 9's "+4 s" is the only guard. Q3 asserts `!isempty(merged)`,
    also unmeasured on the package: with no merged key the task stops and reports.
    Settled (P5-46): on f9361cf, with the plan's extracted law and key code, the ten
    graphs (of 87) hold 101 siblings (2 to 15 each); the pool has 275 mechanisms, 259
    keys, none capped, and 12 merged groups of 28; the testset passes (304) in 0.5 s
    after the pool's first laws and keys take 1.7 s (`refresh5/settle_2.log`).
14. Task 5 routes the RE-speed and ε-path oracles (`law = :package`) and the rank of
    guard-rejected fixtures through `_law`, which throws on RE turnover where the old
    engine returned a finite law. No checked call site reaches such a mechanism, but not
    all 98 rank call sites were checked; Step 6's "a count that moves is a bug" is the
    guard. Settled (P5-46): a recorder on f9361cf's TME saw 171 rank calls on 162
    mechanisms (39 allosteric); one is guard-rejected (`_SPLIT`) and none has a turnover
    rapid-equilibrium cycle on a law mask, so `_law` throws at none
    (`refresh5/settle_5.log`).
15. test_law.jl's ":OnlyI law is the exact limit" (Task 4) overlaps test_oracles.jl's "an
    :OnlyI law is exact, balanced and counts its declared dimension" (phase 4, before the
    fix wave) at the same seeds, `400 + i`; FW-23's testset draws at `600 + i`. The plan
    keeps both, a few seconds' duplication.
16. Task 5 moves the class-5 throw assertions to `_assert_support(m, _law(m))`. The `POLY`
    failure path stays in src until Task 13 with no test (no public path reaches it; the
    guard rejects every hand case class 5 refuses). Acceptable as a backstop until Task 13,
    or keep both assertions and change Task 13's deletion list? Settled (Task 5 ruling):
    both assertions stay beside the new ones until Task 13 deletes them with
    `_num_den_exprs`.
17. Task 6's identify.jl GC workload runs real CMA fits (three rounds, `--pkgimages=no`),
    untimed; it may approach its 1,200 s timeout. Settled (P5-46): phase 4's Task 14 ran
    the same workload with `--pkgimages=no` and the remarks in 41 s on 1.12 and on 1.13
    (`logs/P4_T14_remarks_identify_112.log`, `_113.log`); `src/identify_rate_equation.jl`
    and `src/fitting.jl` have not changed since.
18. Task 4 Step 8 cannot say whether the head's MWC path derives the cap fixture or times
    out; is a 300 s timeout (`exit 137`) acceptable red evidence for the summed cap?
    Settled (P5-46): it derives; on f9361cf the probe prints `rate_equation_string after
    4.6 s: derived` and exits 0 (`refresh5/settle_6.log`), which Step 8 now expects.
19. Tasks 4 and 6 assume `trees[1]` is free E's segment in both masks and that
    `_bare_law`'s first variables are the concentrations in `L.vars` order. Resolved at
    assembly: `_bare_law` adds the concentrations first, and f9361cf's `_segments` roots
    segments in form order, so segment 1 holds form 1.
20. The P5-25 barrier depends on Task 6's probe finding sub-0.5 ns specs, and the
    `_rate_body` precompile lines on the first-call budget. Should either be unconditional?
    Settled (Task 6, T6-M5): the `_rate_body` precompile lines are unconditional; the
    first call reads 0.11 s (Task 6: 2.43 s) for package precompile 10.5 s (5.9 s); the
    barrier stayed conditional, and the probe found 4 of 42 on 1.12.5 and 1.13.1.
21. FW-07's call text at the seven test_admission sites. Settled (P5-46): each reads
    `EnzymeRates._support_verdict(<x>, _testhelper_law_branches(<x>))` (test_admission.jl
    :1080, :1269-1270, :3137, :3158, :3171, :3213, :3278 at f9361cf), as Task 5 quotes it.
22. `P5_compare.jl` must read allosteric rows through `_law(am)`, and `census/P5_T6` must be
    built with `census/P5_base`'s arguments. Resolved at assembly: `law_value` dispatches
    through `hasmethod` once Task 4 lands, and Task 6 Step 11 passes the rename map.

Tasks 8–11:
23. RES §3.2's printed-N figures (Segel ter bi, the power case) were left by Task 11 to the
    printer. Resolved at assembly: Task 12's test_print.jl testset pins F1's and F4's `N`
    exactly and the other two by shape.
24. Part 2 changed while Part 3's check ran. Resolved at assembly: Task 11's Q1 anchors and
    Task 8's use of `_Canonical` (with `leaves`), `_canonical` and `EnzymeRates.sha256`
    match the final Tasks 4 and 7.
25. Q5's and Q6's counts (21 flagged of 200; 172 children) and the TME verdicts were
    measured on a9e0699 with extracted Part 1–2 code; FW-09 and FW-33 could move them.
    Task 10 Step 8 rechecks before pinning; the refresh pass must re-measure them.
    Settled (P5-46): on f9361cf the same probe (`p3_probe.jl`) prints every verdict and
    count a9e0699 gave (TME's K-type fixtures, #3268, Q5's 200, 21 and 50, Q6's 172, the
    units), and Task 10 Step 7's Q5 and Q6 testsets pass there with their pins
    (`refresh5/settle_12.log`, `settle_2.log`).
26. Q6's warm time is unmeasured (8.7 s cold, compile included); the 10 s cap on Q5 + Q6
    rests on Task 10 Step 8's stop rule. With Q1–Q4 (about 10 s), the phase-5 population
    tests sit near P5-41's 20 s. Settled (P5-46): on f9361cf, with the extracted code and
    the head's rank oracle in a fresh process, Q5 takes 12.8 s (its first calls compile
    the certificate, the laws and the ranks' rate equations) and Q6 0.3 s after it; warm,
    0.2 s each (`refresh5/settle_2.log`). Q5 run cold alone exceeds the 10 s cap.
27. TDD against value-check-first: the TME verdict pins, Q5, Q6 and Q1's units run only
    green, since each follows a value check that needs the code; only the certificate's and
    the unit rule's own tests are watched failing. A cheap red run of TME and test_oracles
    may be wanted.
28. No Julia run checked Task 8's `_process_batch` and `_cv_model_selection` or Task 9's
    `_select_count!` as code; they were traced by hand against the snapshot and RES §1.9
    (`src/identify_rate_equation.jl` is the same at f9361cf), and their tests are the
    first execution.
29. P5-36: no one-subunit ping-pong tagging child of 30 parents binds nothing in the
    inactive conformation or has a certificate (0 of 172). Q6 pins that property; only the
    rerun looks further. Does Denis want a round-2 ping-pong population in the suite?
30. Spec §10 and §13 still say a phantom "inherits the parent's fit"; under P5-32 the source
    is whoever holds its L→0 key. Resolved at review: Task 10 rewrites §13's bullet and
    §10's, listed with P5-30 (◆ P5-32).
31. `model_selection.md` and `developer.md` name `eq_hash` from Task 8 until Task 15
    rewrites them. Acceptable for that interval?
32. PASS1's check that a phantom's names match its L→0 order cannot fire, since the
    certificate rules out regulator constants; it stays, untested.
33. Task 8's value script relies on Task 7's `swap` control; when the two swapped constants
    play symmetric roles the control says nothing, and Step 1 records it.
34. Task 9's gate needs the selected key to equal the true mechanism's; neither Task 0 nor
    the research ran this search, so another selection stops Task 9 for investigation.
    Resolved at review: the data now put products at 10–100% of their substrates, and
    when Task 0's run misses the truth, Task 9 may select Task 0's law instead (under
    "Deviations from the research"). Refresh (P5-46), for Denis: on f9361cf the run
    misses the truth and selects a 7-constant L phantom; reseeding the optimizer alone
    selects another 7-constant law (one tag differs; another key, the same L→0 law). On
    these noise-free data 47 of the 56 LOOCV candidates of each run fit to the noise
    floor (scores 2e-13 to 1.1e-10, 6 to 17 constants; the next scores 1.9e-3), so which
    of them is selected follows the optimizer's random stream. Task 9's gate therefore
    also accepts a near-exact law (`cv_score` below 1e-6) with at most Task 0's selected
    `n_params`, and investigates any other selection. Noisy data, or a truth the search
    reaches, would restore a key-level gate.

Tasks 12–15:
35. `PASS1_INSTANCES = 0` was never measured. A fresh pick can compile a branch the warm-up
    pick skips; Task 12b Step 2 then stops for Denis rather than raising the budget. Should
    the gate allow a small constant?
36. Segel ter bi's and the power case's printed `N` are checked by shape, taken from UNIT
    §2's slice-name output; nobody printed them in package names. F4's exact `N` is the
    writer's probe text.
37. The checker's edits to Tasks 12–15 were parse-checked only (`p4c_parse.jl`): the TID
    guard testset, the `N` pins, the `isa Rational{Int}` assertion, the typed
    comprehensions, `scale`, `_testhelper_hyperbolic`. The guard message comes from Task
    3's pin. Settled (P5-46): the fix wave kept class 2's message, and at f9361cf PASS1
    records "Mechanism: the rapid-equilibrium steps E + P ⇌ E(P), E + S ⇌ E(S) and
    E(S) ⇌ E(P) close a cycle …" for the TID fixture, which holds the asserted phrase.
38. The printer probe ran Parts 1–3 as stand-ins taken before they last changed (diffed:
    docstrings and one `pow` simplification; the normalize stand-in matches Task 11's
    `_normalizer`). The pinned texts (the 11 factored pairs, the allosteric fixture, the
    doctests) should be reprinted with `P5_T12_pins.jl` on the real Task 11 code before
    Task 12 runs.
39. PASS1's `_check_well_posed` call goes beyond RES §1.9's pseudo-code, on RES §1.1 and
    R11 (every derivation entry runs the guard) and to keep today's PASS1 behaviour.
40. Not strict TDD: Q1's printed forms, TID's CSV-text check, the `PASS1_INSTANCES`
    testset and the TID guard testset pass on the head or follow the code they check; the
    printer and pin 17 are watched failing.
41. `P5_census_nostr.jl` reads `laws.tsv` column 8 and compare.jl's "strings exact"; both
    were verified in the census source, neither run. Settled (P5-46): both ran on
    f9361cf's corpus `census/P4_FW` against itself: `laws ok in both: 542; printed text
    differs: 0`, then `RESULT: PASS` with every `strings exact` full
    (`refresh5/settle_10.log`, `_11.log`).
42. Should the coefficient form also print the `params` and `concs` lines? The plan keeps
    the research's layout: `v`, the comment, `N`, the legend.
43. Task 12 leaves `_equation_text` and `_assert_no_name_clash` as dead code until Task 13,
    which deletes the old engine whole. Delete them at once (CLAUDE.md: remove unused
    code)? The plan keeps them one task, and Task 12's commit says why: Task 13 deletes
    the old engine and printer as one unit under one grep, and deleting them first would
    orphan `_expr_to_string` and its helpers piecemeal.
44. Task 12 adds about 8 s against RES §4's 2 s; accept, or sample the 42-spec printed-law
    check?
45. Task 15's tutorial `v` line is quoted from reasoning; Step 6 checks it against the
    rendered page and follows the page.
46. Task 14: a c below 0.8 still edits CLAUDE.md and goes to Denis. Apply a tighter c
    without the CLAUDE.md change?
47. R52's figures and the Tasks 12–13 census run on graph-guard's final head, whose `src`
    is f9361cf's (991355f changes no `src` file). Settled (P5-46): f9361cf gives Task 13's
    check the snapshot's populations, 439 and 399 level-1 mechanisms (340 and 300
    allosteric) under the summed filter; only the zero-difference result is a gate.

---

## Appendix: the probes behind the plan

Each part of this plan was written and checked against the a9e0699 snapshot. These are the
runs its figures come from, as each part recorded them.

Tasks 0–3 (`$S/phase45/plan5/probes/`):

- `p1_engine_check.jl`, on the code before its final edits:
  - `_dependent_param_exprs` through `_solve` was byte-identical on the 42 specs, and
    `all_parameters` agreed to 1e-14;
  - 41 `test_poly.jl` tests passed;
  - 446 `test_law.jl` tests passed in 2.8 s, the oracle within 1e-20.
- `p1_plan_check.jl` ran the code blocks extracted from Tasks 0–3 (`p1_extract.py`):
  - the solve, `all_parameters`, `poly.jl` and `engine.jl` were evaluated into the
    package;
  - every new testset passed, with a9e0699's `test_accessors.jl` and
    `test_joint_solve.jl`: 1159 tests.
- `p1_baseline_dry.sh` ran Task 0's scripts on the snapshot; Task 0 quotes its figures.
- `c1_check.jl` (the checker's) ran those code blocks again after the checker's
  fixes: 557 tests passed, `test_law.jl`'s 454 among them; the new `_eq_complexity`
  equalled the head's on the 480 dry-run baseline rows, and the six gate laws' new D sizes
  read 396, 396, 396, 589, 799, 858.
- `c1_red_solve.jl` and `c1_red_poly.jl` ran the new tests on the unchanged snapshot: the
  converted accessor testset passed (3 tests), and the failures read as Task 1 Step 3 and
  Task 2 Step 2 say.

Tasks 4–7:

Measured on the snapshot (`$S/phase45/plan5/probes/p2_*.log`):
the cap fixture (random-order bi-bi plus E(A, Q) and E(B, P), chemistry `:NonequalAI`,
n = 1) is admitted with V×τ 4,320 per conformation; pin 12's fixture counts 9 + 9; pin 13's
parent counts 8,832 + 2; the inert uni-uni counts 1 + 0; Task 6's cube-root fixture gives
`k_F_to_E = Keq ^ -0.333… · …` (`p2_cube.log`; it is the slice's `cube_root` with a second
substrate B bound in `F(S) + B <--> E(P)`, since the slice's own fixture solves to integer
exponents on the snapshot, `p2_fixtures.log`); the key-attack fixtures and a metabolite
named `x1` pass the macros; the plain and MWC steady-state dead ends of test_admission are
rejected by the guard (`_LAW_MULTI`, `_ISLAND`), so no public path reaches class 5 on them.

Tasks 8–11:

Three probes ran Tasks 8–11's code on the a9e0699 snapshot with Tasks 1–3's extracted
solve, polynomials and engine and Task 4's code (`$S/phase45/plan5/probes/p3_probe.jl`,
`p3_probe2.jl`, their logs; 418 s and 22 s, under 1.4 GB); the third, `p3c_probe.jl`, adds
Task 7's key (283 s, 1.3 GB):
- the certificate and `_is_phantom` give the TME K-type verdicts that match TME's rank
  pins, with every map holding, and census #3268 is admitted and flagged;
- the overflow toy: p = 2, q = 3 gives w = [−1/2, −1/6]; p = 2^40 + 15, q = 2^40 + 3
  overflows to no certificate;
- the unit rule gives each Task 11 fixture its pinned rule and weights; the guard
  admits F1–F4 and the power case; the units check passes on the 42 specs;
- Q5: 280 tagging children of 40 bi-bi seeds, 200 picked, 21 flagged, 50 through the
  pre-filter, every map holding, four rank pairs equal, about 6.5 s;
- Q6: 172 children of 30 ping-pong parents, 76 K-types, none flagged, none with a
  certificate;
- the dedup twins, the LDH renaming pair and the LDH allosteric split/merge pair have laws
  equal field by field; the uni-bi pair's fitted names differ in one constant;
- the pinwheel has V×τ 192 and 15 fitted constants, and its law is invariant under every
  permutation of its cycles; its key search stops at the cap, 100 leaves, in about 4 s,
  and gives the per-mechanism key;
- the dedup twins and the uni-bi pair each have one key, one law under the canonical map,
  and a swapped map that fails; the uni-bi pair's canonical positions sit at different
  fitted positions, the twins' at the same;
- every flagged TME K-type fixture, the two certificate-only Theorell–Chance variants and
  census #3268 have an active law whose key is their catalytic mechanism's, and their
  canonical names are the fitted names but L;
- Q1's units check passes on the 120 sampled laws at one subunit (346 tests, 1.2 s).

Tasks 12–15:

Measured on the snapshot (probe `$S/phase45/plan5/probes/p4_print_check.jl`, log
`p4_print_check.log`): Task 12's `src/print.jl`, loaded over Tasks 1–3's `_solve`,
`poly.jl` and `engine.jl`, Tasks 4–5's MWC law and class 5, and a stand-in of Task 11's
interface (`p4_code/normalize_standin.jl`, the research's `spec_norm`), printed every
text Task 12 quotes. Both forms of all 42 specs equal
the derived law to 6e-77 in 256-bit BigFloat, with no artifact; the reference lines equal
`all_parameters` to 1.4e-16; the LDH i-state laws define every name they read; every
non-metabolite name of both forms holds `_` or is reserved; F1–F4 print as UNIT §2 has
them; the 13 `:OnlyI` fixtures pass the guard except `copy_group_bridge` and
`essential_activator`; the ter-bi gate law prints in 28 ms (59,645 characters) and its
coefficient form in 24 ms. A second probe (`p4_plan_check.jl`, log `p4_plan_check.log`)
ran Tasks 12–14's test code as the plan writes it (extracted by `p4_extract2.py`): Task 12's
four test_print.jl testsets, the factored-pin table, pins 3, 5–8, 17 and R37, the Q1
addition on the 42 specs, Task 13's converted segment-root and free-E-weight testsets
(`occupancy` applied) and Task 14's c script: 892 tests, all passing; and Task 13's
`P5_T13_den_check.jl` found 0 differences on its four sets. Texts the D term order decides
follow `_term_key` (P5-11); a term-order difference in a Task 12 pin is a mismatch with
Task 11's `_term_key`: compare it with P5-11 before any re-pin.
The edits made after those probes (PASS1's guard and its TID testset, now Task 12b; the
`N` pins and pin 17's `isa` check in Task 12; Task 13's helper deletion, docstring move and
`_testhelper_hyperbolic`; Task 14's run lookup; Task 15's `cha_king_altman.md` edits; the
typed comprehensions and the single-assigned `scale` of `src/print.jl`) were parsed
(`p4c_parse.jl`, log `p4c_parse.log`: every block parses but the quoted one-line
fragments and templates) and not run. The printed `N` of F4 is the probe's text, and those
of Segel's ter bi and the power case are checked by shape (Task 12 Step 6).
