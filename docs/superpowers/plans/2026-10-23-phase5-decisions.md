# Phase 5 decisions Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Carry out Denis's 2026-10-23 answers to draft PR #95's open questions on branch
graph-law, so #95 lands phase 5 complete.

**Architecture:** Eight small tasks on top of 4427dc91, each one decision or a tight group:
the printer's normalizer and layouts (Tasks 1-3), the L-phantom rule (Task 4), the beam's
key sets (Task 5), the parameter-name API (Task 6), two test-only changes (Task 7), and the
spec and CLAUDE.md text (Task 8). Evidence runs and the research push follow as controller
steps (Task 9). No task touches code generation, so `rate_equation`'s body is unchanged.

**Tech Stack:** Julia 1.12 (local), 1.13 (CI); Test, TestEnv; Documenter.

**Spec:** `docs/superpowers/specs/2026-10-12-reaction-graph-rewrite-design.md`; the decisions
are `/home/denis.linux/claude-scratchpad/EnzymeRates-12f42995/phase45/denis_answers5.md`
(one line each, quoted below by their ids).

## Global Constraints

- CLAUDE.md binds: TDD for every behaviour change (red run watched failing for the right
  reason, then green); no mocks; pristine output; ABOUTME headers; `_testhelper_` prefix on
  every test helper; 92 columns; no hand whitespace changes outside edited lines; evergreen
  names and comments (no task, phase, ruling or decision ids in `src/` or `test/`);
  ROADMAP.md updated in the commit that changes an entry.
- `rate_equation` stays allocation-free, under 120 ns on every spec and under
  `max(120 ns, 0.8 ns × terms)` on the gate laws. No task here may change `_rate_body`'s
  output; a change that would is BLOCKED for Denis.
- Never raise a compile budget. PASS1 stays at 0 method instances per fresh mechanism
  (item 35). The first-call budget falls to 1.5 s (DN-01).
- Julia only through `$S/phase3/jl.sh <log> <timeout> <args>`, one process at a time
  (7.7 GB, no swap); focused runs through TestEnv; the full suite with
  `--heap-size-hint=2500M`. `S=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad`.
- A pin a run contradicts is checked against an oracle or the exact law before it is
  re-pinned; each re-pin names its evidence in the commit body.
- Commit trailers, exactly:
  `Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>` and
  `Claude-Session: https://claude.ai/code/session_01PJZFTkVDnVn4EyCYVUNmWE`.
- Version stays 0.10.0.

## Review Focus

1. An allosteric conformation with no constant term (an allosteric ping-pong): Cleland's
   rule applies per conformation, with that conformation's own terms (Task 1 test).
2. A one-substrate law with no constant term (fixtures F2, F4): Cleland's rule would need
   ΠS/A = 1, the constant term, which is absent, so the graded basis still applies; F2's
   and F4's existing testsets must stay green unchanged (Task 1 Step 4).
3. The coefficient form of an MWC law with regulatory sites and `N_A`/`N_I` lines must
   evaluate top to bottom to the law (Task 2 test over all allosteric specs).
4. A certificate-only phantom with regulator constants meets PASS1's name check; it must be
   demoted, not mislabelled (Task 4, item 32).
5. The beam's buckets are a `Dict` keyed by count; selection must not depend on its
   iteration order (Task 5 test with two counts in one sweep).

## File Structure

| File | Change |
|---|---|
| `src/normalize.jl` | Cleland's rule between the constant term and the graded basis |
| `src/print.jl` | coefficient form in evaluable order; docstrings |
| `src/identify_rate_equation.jl` | PASS1 prints the coefficient form; `_select_count!` loses `rejected`; `_sweep_parents!` |
| `src/key.jl` | `_is_phantom(A)` is the certificate alone |
| `src/rate_eq_derivation.jl`, `src/thermo.jl` | `parameters`/`fitted_params` from `_solve`; the renderers leave `src` |
| `test/mechanism_definitions_for_test_enzyme_derivation.jl` | `_testhelper_dependent_param_exprs`, `_testhelper_power_expr`; `_testhelper_eval_expr` goes |
| tests named per task | new and updated testsets |
| `docs/src/...`, `ROADMAP.md`, `.claude/CLAUDE.md`, the spec | text per task |

---

### Task 1: Cleland's normalizer before the graded basis (P5-10a, item 36)

**Files:**
- Modify: `src/normalize.jl` (`_Norm`'s docstring, `_normalizer`)
- Modify: `src/print.jl` (`rate_equation_string`'s rule list), `docs/src/deriving/textbooks.md` (the unit rule's list)
- Test: `test/test_print.jl`

**Interfaces:**
- Produces: `_normalizer(D, vars, substrates) -> _Norm` with a fourth rule symbol,
  `:cleland`, tried after `:constant` and before `:graded`.

- [ ] **Step 1: Write the failing tests** in `test/test_print.jl`, after "F4 takes integer
  weights A:2, A²:−1":
  - `@testset "a ter-bi ping-pong takes Cleland's normalizer"`: for spec
    "Segel Bi Uni Uni Uni PP Ter Bi", `nz.rule === :cleland` and
    `nz.lambda == [BC => 1//1, AC => 1//1, AB => 1//1, ABC => -2//1]` (monomials through
    `_testhelper_conc_mono`), `nz.shift == EnzymeRates._Mono()`; `_testhelper_check_units`;
    and, at 256 bits on random positive constants, the new `N` equals the graded `N`
    divided by `(k_FC_to_E_Q + k_EAB_to_F_P)^2` (the graded `N` is
    `(k_EAB_to_EA_B + k_EAB_to_F_P) * (k_FC_to_E_Q + k_FC_to_F_C)`).
  - `@testset "a bi-bi ping-pong's Cleland normalizer is its graded one"`: "Segel Ping Pong
    Bi Bi" gives `:cleland`, and its `N` equals, at 256 bits, today's
    `(k_FB_to_E_Q + k_FB_to_F_B) * (k_EA_to_E_A + k_EA_to_F_P) / (k_EA_to_F_P + k_FB_to_E_Q)`.
  - `@testset "each MWC conformation takes Cleland's rule on its own terms"` (Review Focus
    1): an inline one-subunit allosteric bi-bi ping-pong (`:EqualAI` throughout); both
    conformations' `_normalizer` give `:cleland`; `_testhelper_printed_law_error(am) <=
    big(10)^-20`.
- [ ] **Step 2: Run them; they fail** (`:graded` where `:cleland` is expected).
  Log `$S/phase45/logs/P5D_T1_red.log`.
- [ ] **Step 3: Implement.** In `_normalizer`, after the constant-term return: with
  `S = _monomial((v => 1 for v in unique(substrates))...)` and the pairs
  `_mmul(S, _monomial(v => -1)) => 1` for each substrate variable `v`, then
  `S => 1 - k` (`k` substrate variables), return `_Norm(:cleland, pairs, _Mono())` when
  every monomial is a key of `Dc`. Update `_Norm`'s and `_normalizer`'s
  docstrings (the rule list, `:cleland`), `rate_equation_string`'s rule list, and
  textbooks.md's list ("Cleland's normalizer, when the denominator holds ΠS/X for every
  substrate X and ΠS: the substrate terms then read as products of Michaelis constants,
  `A·B·C/(K_a·K_b·K_c)`").
- [ ] **Step 4: Re-pin exactly (item 36).** Run test_print.jl, test_allosteric_golden.jl,
  test_oracles.jl and the doctests (docs build). For each printed `N` or golden law that
  changes, check its value against the law at 256 bits (`_testhelper_printed_law_error`)
  and then pin the new text exactly, including Segel ter bi's and the ping-pong power
  case's `N`, which are pinned by shape today. List every re-pin in the commit body.
- [ ] **Step 5: Full suite**, pristine; log `P5D_T1_full.log`. Commit:
  "Print Cleland's normalizer when the denominator holds its terms".

### Task 2: The coefficient form evaluates top to bottom (item 42)

**Files:**
- Modify: `src/print.jl` (both `_law_text` methods, their docstrings, `rate_equation_string`'s doctest)
- Modify: docs pages that show a coefficient form (find with `grep -rn "coefficients = true" docs/src`)
- Test: `test/test_print.jl`, `test/mechanism_definitions_for_test_enzyme_derivation.jl`

**Interfaces:**
- Consumes: Task 1's `_normalizer`.
- Produces: `_law_text(m, s, L; coefficients = true)` lines, in order: the header
  (`_header_lines`: params, concs, reference lines), `_COEFFICIENTS_COMMENT`, `N`/`N_c`
  lines, the coefficient definitions, for MWC the `Q_c`/`D_c`/`F_c`/`W_c` definitions in
  coefficients, and last the `v` line.

- [ ] **Step 1: Write the failing test** `@testset "both printed forms paste and evaluate
  as Julia"` (Review Focus 3): for every spec in `MECHANISM_TEST_SPECS` and the F1-F4
  fixtures, and for both `coefficients = false` and `true`, evaluate the text with
  `Base.include_string` in a fresh `Module()` in which `params` and `concs` are bound to
  NamedTuples of random positive Float64 values; the result equals
  `rate_equation(m, concs, params)` to `rtol = 1e-10`. Also: the first line ends with
  `) = params` and the last starts with `v = `. (`_testhelper_printed_law_error` already
  checks both forms' values in an order-free reader; this test checks the order.)
- [ ] **Step 2: Run; it fails** (the first line is `v = …`). Log `P5D_T2_red.log`.
- [ ] **Step 3: Implement** the order above in both `_law_text` methods; update their
  docstrings and the doctest at `rate_equation_string`.
- [ ] **Step 4: Run test_print.jl, test_allosteric_golden.jl, test_compile_budget.jl and
  the docs build**; re-pin coefficient-form texts (the second-process `COEF` hashes in
  test_compile_budget move; the values do not).
- [ ] **Step 5: Full suite**, pristine; commit "Print the coefficient form in the order it
  evaluates".

### Task 3: identify's CSV holds the coefficient form (P5-14)

**Files:**
- Modify: `src/identify_rate_equation.jl` (PASS1's `eq_text`, the `rate_equation` column in the docstring)
- Modify: `test/test_compile_budget.jl` (PASS1 replica), `docs/src/identify/model_selection.md` (the `rate_equation` row)
- Test: `test/test_identify_rate_equation.jl`

**Interfaces:**
- Consumes: Task 2's layout. `rate_equation_string(m)` keeps the full law as its default.

- [ ] **Step 1: Change the failing assertion first**: TID's
  `bestf.rate_equation[1] == rate_equation_string(…)` check reads
  `rate_equation_string(…; coefficients = true)`. Run; it fails. Log `P5D_T3_red.log`.
- [ ] **Step 2: Implement**: PASS1's `eq_text = _law_text(m, s, L; coefficients = true)`;
  the replica in test_compile_budget does the same. The `rate_equation` row reads "The
  printed law in named coefficients (`rate_equation_string(m; coefficients = true)`)".
- [ ] **Step 3: Run TID and test_compile_budget**: PASS1 instances stay `[0, 0, 0, 0]`, fit
  78. If PASS1 compiles, warm the coefficient path in the budget test's warm-up, never
  raise the budget.
- [ ] **Step 4: Full suite**, pristine; commit "Write each law's coefficient form to
  identify's CSVs".

### Task 4: The certificate alone makes an L phantom (P5-31a, item 32)

**Files:**
- Modify: `src/key.jl` (`_is_phantom`), `src/identify_rate_equation.jl` (PASS1's call, its comments and docstrings)
- Modify: `docs/src/identify/model_selection.md` (`L_phantom` row), `docs/src/identify/enumeration_engine.md`, `docs/src/developer.md`, `ROADMAP.md` (the phantom entry)
- Test: `test/test_key.jl`, `test/test_identify_rate_equation.jl`, pins in `test/test_oracles.jl`

**Interfaces:**
- Produces: `_is_phantom(A::_AlloLaw) -> Bool = !isempty(_phantom_weights(A))`. The
  certificate already requires one subunit and an inactive conformation that turns nothing
  over (`A.n == 1 && iszero(A.inactive.Q)`).

- [ ] **Step 1: Write the failing test** in test_key.jl, `@testset "an inactive conformation
  that binds substrate can be an L phantom"`: an inline one-subunit
  `@allosteric_mechanism`, uni-uni, S binding `:EqualAI`, catalysis `:OnlyA`: `_is_phantom`
  is true and the certificate's map reproduces the law (`v(θ, L) = v₀(θ ⊙ (1+L)^w)` at
  256 bits, as the existing certificate testset checks). Control in the same testset:
  S binding `:NonequalAI`, catalysis `:OnlyA`: false.
- [ ] **Step 2: Run; the first assertion fails.** Log `P5D_T4_red.log`.
- [ ] **Step 3: Implement** `_is_phantom(A)` and update its one `src` caller and every test
  caller (`git grep -w _is_phantom`).
- [ ] **Step 4: Item 32.** Over the bi-bi population with regulators (the research's
  `bibi_R`; regenerate with `_testhelper_population` on the reaction with a regulator),
  run PASS1's name check (`Set(canon0.order) != Set(s.fitted[1:(end - 1)])`) on every
  certificate-only phantom; log `P5D_T4_namecheck.log`. If one trips it, add it inline as
  `@testset "a phantom whose source lacks its constants is fitted itself"` asserting PASS1
  demotes it (`phantom == false` in its `_Candidate`). If none does, delete the check and
  its comment, and say so in the commit body (Denis's item 32 answer allows it).
- [ ] **Step 5: Re-pin the phantom counts** (test_oracles' population figures, TID's phantom
  testsets) against the certificate-alone figures; each new count is checked by the rank
  oracle on that population before it is pinned.
- [ ] **Step 6: Text**: `_is_phantom`'s docstring, identify's docstring (the `L_phantom`
  column and PASS-2), model_selection.md's `L_phantom` row ("a one-subunit allosteric variant
  whose inactive conformation turns nothing over and whose `L` an exact rescaling of its
  other constants absorbs …", keeping FW-32's timing clause), enumeration_engine.md,
  developer.md, ROADMAP.md's phantom entry.
- [ ] **Step 7: Full suite** and docs build; commit "Let the rescaling certificate alone
  decide an L phantom".

### Task 5: One key set and a beam-level test of it (DN-03, DN-04)

**Files:**
- Modify: `src/identify_rate_equation.jl` (`_select_count!`, `_beam_search`, new `_sweep_parents!`)
- Test: `test/test_identify_rate_equation.jl`

**Interfaces:**
- Produces: `_select_count!(expanded, best_loss_by_count, c, losses, law_keys; loss_rel_threshold, loss_abs_threshold, loss_parsimony_threshold, min_beam_width, selected::Set{String}) -> Vector{Int}`;
  `_sweep_parents!(frontier::Dict{Int, Vector{BatchEntry}}, target::Int, expanded::Dict{Int, Int}, best_loss_by_count::Dict{Int, Float64}, selected::Set{String}; loss_rel_threshold, loss_abs_threshold, loss_parsimony_threshold, min_beam_width) -> Vector{BatchEntry}`,
  which pops every bucket at a count up to `target` and returns the entries of selected
  keys, bucket by bucket; `_beam_search`'s sweep calls it.

- [ ] **Step 1: Write the beam-level test** `@testset "a sweep keeps each key's verdict"`
  (Review Focus 5). Hand-built `BatchEntry`s (any inline mechanism, `row = (;)`):
  - Sweep 1, target 3, count 3: `a` (loss 1.0), `b` (loss 5.0); thresholds rel 1.3, abs 0,
    parsimony `Inf`, `min_beam_width = 1`, `best_loss_by_count[3] = 1.0`. Returns `[a]`.
  - Sweep 2, target 3: bucket 3 holds `a2` (key of `a`, loss 1.0), `b2` (key of `b`, loss
    5.0) and `c` (new key, loss 0.5); bucket 2 holds `d` (loss 0.2) so two counts share the
    sweep; set `best_loss_by_count[3] = 0.5`, `[2] = 0.2`. Returns exactly the entries
    `a2`, `c` and `d` (as a `Set` of `objectid`s, so bucket order cannot matter): `a2` by
    `a`'s verdict, though its loss is above the new cutoff 0.65 and count 3's budget is
    spent; `c` under the cutoff; `d` within count 2's budget; `b2` judged again and
    refused.
- [ ] **Step 2: Implement `_sweep_parents!`** by moving the sweep's loop out of
  `_beam_search`; run the test green. Then watch it fail under each mutation, restoring
  after each (logs `P5D_T5_mut_rebuild.log`, `P5D_T5_mut_reverse.log`): (a) `selected`
  replaced by a fresh `Set{String}()` inside the call; (b) the keys passed to
  `_select_count!` taken from `reverse(entries_at_count)`.
- [ ] **Step 3: Drop `rejected`** from `_select_count!` and `_beam_search`: a key not in
  `selected` is judged again by the same rule. Update the five TID call sites, delete the
  four assertions on `rejected`'s contents, keep "a late member of a rejected key stays
  out" with its selection assertions. Update `_select_count!`'s docstring (a key keeps a
  selected verdict; a rejected key's late member is judged again and rejected, since the
  cutoff only tightens and a rejection by rank spends the budget), the `seen` comment and
  identify's docstring.
- [ ] **Step 4: Full suite**; commit "Keep one key set in the beam and test it across
  sweeps".

### Task 6: Parameter names straight from the solve (DN-09, DN-10)

**Files:**
- Modify: `src/rate_eq_derivation.jl` (`parameters`, `fitted_params`), `src/thermo.jl` (delete `_dependent_param_exprs` ×3 and `build_power_expr`)
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`, every test calling the moved functions, `docs/src/developer.md:89`
- Test: the existing tests

**Interfaces:**
- Produces: `parameters(::M) = (_solve(M).fitted..., :Keq, :E_total)` and
  `fitted_params(::M) = Tuple(_solve(M).fitted)` (both stay `@generated`);
  `_testhelper_dependent_param_exprs(x)` (the solve, a mechanism or its type) and
  `_testhelper_power_expr(keq_exp::Rational, factors)` with the bodies moved verbatim;
  `_testhelper_exponent(x::Rational) = x`.

- [ ] **Step 1: Move** the two functions into the fixtures under the new names, rename every
  test call (`git grep -w _dependent_param_exprs build_power_expr` must print nothing in
  `src/` and only the helpers' definitions under their new names in `test/`).
- [ ] **Step 2: Rewrite `parameters` and `fitted_params`**; delete the src definitions.
- [ ] **Step 3: DN-10**: TR's half-integer testset evaluates `dep[:k_F_to_E]` with
  `_testhelper_printed_value`; add `_testhelper_exponent(x::Rational) = x`; delete
  `_testhelper_eval_expr`.
- [ ] **Step 4: Run test_accessors, test_joint_solve, test_oracles,
  test_rate_eq_derivation, test_types, test_allosteric_collapse**, then the full suite
  (Aqua, JET, first call). Commit "Read parameter names from the solve and keep the
  dependent renderer in the tests".

### Task 7: Tighter first-call budget and fresh seeds (DN-01, item 15)

**Files:**
- Modify: `test/test_compile_budget.jl`, `test/test_oracles.jl`

- [ ] **Step 1: DN-01.** `RATE_EQUATION_WALLCLOCK_BUDGET_S = 1.5` with the comment
  "baseline 2026-10-23: 0.12 s locally, 0.20-0.37 s on CI"; add
  `@testset "the generator is in the package image"` with
  `@test precompile(EnzymeRates._rate_body, (EnzymeRates.Mechanism,))` and the
  `AllostericMechanism` twin. Watch the budget fail with both `precompile(_rate_body, …)`
  lines removed from `src/EnzymeRates.jl` (log `P5D_T7_mut.log`), then restore.
- [ ] **Step 2: Item 15.** test_oracles' "an :OnlyI law is exact, balanced and counts its
  declared dimension" draws its mechanisms from a seed base no other testset uses (check
  with `grep -n "400 +\|MersenneTwister(" test/test_law.jl test/test_oracles.jl`); its
  assertions stay.
- [ ] **Step 3: Full suite**; commit "Hold rate_equation's first call to 1.5 s and give the
  :OnlyI oracle its own seeds".

### Task 8: Spec, CLAUDE.md and docs text (DN-07, DN-08, P5-10a, P5-14, P5-31a, P5-32)

**Files:**
- Modify: the spec (sections 5, 9, 10, 13), `.claude/CLAUDE.md` (the `rate_equation` paragraph)

- [ ] **Step 1: Spec section 9**: rule 1's second sentence becomes the docs' (FW-31: every
  term dimensionless; in a rapid-equilibrium law whose only form that binds nothing is
  free E, each term is a form's abundance relative to free E); a new rule between 1 and 2,
  Cleland's normalizer when the denominator holds ΠS/X for every substrate X and ΠS; the
  `N = …` sentence names an integer to a fractional power (FW-12); Printing: the
  coefficient form's line order (Task 2).
- [ ] **Step 2: Spec sections 5, 10, 13**: the CSV's `rate_equation` holds the coefficient
  form; an L phantom is any variant the certificate admits (one subunit, an inactive
  conformation that turns nothing over), no pre-filter; a phantom reads its source's fit
  when that fit was made in an earlier iteration or is made in its own, and otherwise is
  fitted unflagged (b2).
- [ ] **Step 3: CLAUDE.md (DN-08)**: "about 4× headroom over the slowest
  `MECHANISM_TEST_SPECS` law ("LDH i-state NonequalAI 6-group, steady-state NAD group", about
  26 ns) and about 9× over the next ("Random-order Bi-Bi", about 13 ns)".
- [ ] **Step 4: Docs build** clean; commit "Bring the spec and CLAUDE.md in line with the
  decided rules".

### Task 9: Evidence and the research archive (controller, no `src` change)

- [ ] **Item 27 and 40 mutation runs.** For each of: the TME verdict pins, test_oracles'
  Q5, Q6 and Q1-units tests, Q1's printed forms, TID's CSV-text check, the
  `PASS1_INSTANCES` testset and the TID guard testset: apply the smallest code change that
  breaks what it guards, run the file focused, record the failure (`P5D_mut_<name>.log`),
  restore (`git diff --stat` empty). A table goes in the PR body.
- [ ] **Item 34 noisy identify run.** `$S/phase45/scripts/P5_identify_noisy.jl`: the
  seeded run of `P5_identify_run.jl` with each rate times `1 + 0.05·randn()` (seeded);
  check with `P5_T9_check.jl` that the selected key is the truth's. If it is not, report to
  Denis with the selected law, its `cv_score` and the truth's; never loosen the gate.
- [ ] **P5-43a.** Push `phase45/research5/` and `phase3/gcprobe/` to `graph-rewrite-reference`
  in a worktree of that branch, under the directory its earlier archives use; check sizes
  first and leave out files over 50 MB, listing them in the PR body.
