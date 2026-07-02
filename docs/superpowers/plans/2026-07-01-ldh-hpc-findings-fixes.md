# LDH HPC Findings — Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the four bugs + one efficiency issue the LDH HPC run exposed in `identify_rate_equation`, per `docs/superpowers/specs/2026-07-01-ldh-hpc-findings-fixes-design.md`.

**Architecture:** Two are allosteric rate-equation **derivation** bugs (`src/rate_eq_derivation.jl`, compile-time `@generated` codegen): §5a undefined inactive-state parameters, §5b an uncancelled num/den product factor. Three are **beam-search** changes (`src/identify_rate_equation.jl`): §1 parsimony reference, §3 progress message, §2 fit-dedup. §6 is a one-line report + a kwarg-forwarding test.

**Tech Stack:** Julia; tests via `Test` stdlib; `Pkg.test()` runner (`test/runtests.jl`). CMA-ES fitting via `OptimizationCMAEvolutionStrategy`.

## Global Constraints

- **`rate_equation` runtime contract:** allocation-free and sub-100 ns per call for every mechanism in `MECHANISM_TEST_SPECS`. Enforced by `test_rate_equation_performance` in `test/test_rate_eq_derivation.jl`. §5a/§5b are compile-time only and must not regress it.
- **Parameter-name chokepoint:** all `Parameter → Symbol` rendering flows through `name(p, m)`; the AST-walker test at `test/test_types.jl:1577-1644` fails the build on any stray `Symbol("K…")`/`k…`/`V…`/`L…` literal outside a parameter-name renderer.
- **Canonical Step Form** is load-bearing; do not reorder/relax step/group canonicalization.
- **Style:** 92-char lines, 4-space indent. Match surrounding code. Files start with two `# ABOUTME:` lines.
- **Every existing test stays green.** Tests that pin specific LDH winners may change once §5 admits shed candidates — update to correct behavior, never to paper over.
- **Reproduction assets** (already on disk): failing mechanism strings live in `docs/ldh_hpc_results/2026_06_27_results/*.csv`; helper repro scripts in the session scratchpad (`repro_numeric.jl`, `investigate_kcat2.jl`, `verify_fix.jl`) show how to reconstruct a mechanism from its `mechanism_type` string via `Core.eval(EnzymeRates, Meta.parse(s))` then `T()`.
- **Do NOT `git add` `docs/ldh_hpc_results/`** (2.4 GB, kept untracked via `.gitignore`).

---

### Task 0: Branch + gitignore setup

**Files:**
- Create/modify: `.gitignore`
- Commit: the approved spec + this plan.

- [ ] **Step 1:** Create the feature branch from the current `several-bug-fixes` tip:

Run: `git checkout -b ldh-hpc-findings-fixes`

- [ ] **Step 2:** Add the results dir to `.gitignore` (create the file if absent; append if present):

```
# Large local-only LDH HPC results (2.4 GB) — never committed
docs/ldh_hpc_results/
```

- [ ] **Step 3:** Verify the results dir is now ignored and NOT staged:

Run: `git status --short`
Expected: `docs/ldh_hpc_results/` does NOT appear; `.gitignore`, the spec, and the plan do.

- [ ] **Step 4:** Commit the design docs (explicit paths only — never `git add -A` here):

```bash
git add .gitignore docs/superpowers/specs/2026-07-01-ldh-hpc-findings-fixes-design.md docs/superpowers/plans/2026-07-01-ldh-hpc-findings-fixes.md
git commit -m "docs: LDH HPC findings spec + implementation plan; ignore results dir"
```

---

### Task 1: §5a — Fix undefined inactive-state parameters (`UndefVarError`)

**Files:**
- Modify: `src/rate_eq_derivation.jl` — add `_i_state_referenced_syms`; rewrite the I-symbol gating in `_build_dep_assignments` (`:1470-1573`) and `_dependent_param_exprs` (`:1321-1420`); reconcile `_synthesized_dep_i_names` (`:1263`).
- Test: `test/test_rate_eq_derivation.jl` (invariant + finite-rate regression); `test/mechanism_definitions_for_test_enzyme_derivation.jl` (fixtures).

**Interfaces:**
- Produces: `_i_state_referenced_syms(am::AllostericMechanism) -> Set{Symbol}` returning `S_I` = { I-state param Symbols appearing in `den_i_poly` } ∪ ({ in `num_i_poly` } when `!_i_state_dead`). Consumed by both `_build_dep_assignments` and `_dependent_param_exprs`.
- Where `num_i_poly`/`den_i_poly` are the existing I-state polynomials:
  `den_i = _rename_symbols(_zero_symbols_in_poly(den_A, _a_only_syms(am)), rename_I)` (and same for `num_A`), with `rename_I = _a_to_i_rename(am)` then `_add_case_b_renames!(rename_I, dep_A, am)`. These are already computed inside `_allosteric_num_den_exprs`; factor them into a shared helper so all sites agree.

- [ ] **Step 1: Reproduce and pin the failing set.** Run the existing repro to confirm the current crash and capture 4 mechanism strings, one per trigger path:

Run: `julia /tmp/.../scratchpad/repro_numeric.jl` (or re-collect from the CSVs). Confirm `UndefVarError` on `rate_equation` for a Case-B reverse (`k_I_ELactateNAD_to_ENADHPyruvate`), a Case-B binding (`K_I_Lactate_ENAD`), an `i_dead` phantom (`kon_I_*`), and a non-`i_dead` (`kon_I_NAD_EPyruvate`) mechanism.
Expected: all four throw `UndefVarError` today.

- [ ] **Step 2: Write the failing invariant + finite-rate test.** In `test/test_rate_eq_derivation.jl`, add a testset that, for each of the four reconstructed mechanisms (store their `mechanism_type` strings as test constants, reconstructed via `Core.eval(EnzymeRates, Meta.parse(s))()`):
  1. asserts every free parameter Symbol in the generated body is defined — parse `rate_equation_string(em)`, collect symbols referenced on the RHS of the `v = …` line, and assert each appears as an LHS (a destructured name in a `(; …) = params`/`= concs` line or an assignment LHS);
  2. asserts `rate_equation(em, concs, params)` returns a finite number at representative concentrations (all metabolites = 1.5) and also at products = 0.

```julia
@testset "§5a I-state parameters are all defined (regression)" begin
    for s in LDH_ISTATE_FAILURE_MECHS  # 4 Sig strings, one per trigger path
        em = Core.eval(EnzymeRates, Meta.parse(s))()
        pnames = EnzymeRates.fitted_params(em); mets = EnzymeRates.metabolites(em)
        params = merge(NamedTuple{pnames}(ntuple(_->1.3, length(pnames))),
                       (Keq=20000.0, E_total=1.0))
        concs  = NamedTuple{mets}(ntuple(_->1.5, length(mets)))
        @test isfinite(EnzymeRates.rate_equation(em, concs, params))  # no UndefVarError
        # DEFINED ⊇ REFERENCED on the rendered transcript:
        @test _every_rhs_symbol_is_defined(EnzymeRates.rate_equation_string(em))
    end
end
```

(Add the `_every_rhs_symbol_is_defined` test helper: split the string into lines; LHS names = names before `=` on destructure/assignment lines; RHS symbols = identifier tokens on the `v = …` and assignment RHS that are not metabolites/`Keq`/`E_total`/`L`; assert `RHS ⊆ LHS ∪ metabolites ∪ {Keq,E_total,L}`.)

- [ ] **Step 3: Run the test to confirm it fails.**

Run: `julia --project -e 'using Pkg; Pkg.test()'` (or load the package in a REPL and run the testset).
Expected: FAIL — `UndefVarError` on the numeric `rate_equation` call for the four mechanisms.

- [ ] **Step 4: Add the shared `S_I` helper.** In `src/rate_eq_derivation.jl`, factor the I-state polynomials out of `_allosteric_num_den_exprs` into a helper and add:

```julia
# ABOUTME line not needed (mid-file). Docstring:
"""I-state parameter Symbols actually referenced by the retained rate-equation
polynomials: `den_i_poly` always (Q_I is kept as enzyme mass), plus `num_i_poly`
when the I-state cycle is live. This is the single source of truth for which
I-state names get defined, replacing the group-structure prediction that drifted
out of sync with the polynomials."""
function _i_state_referenced_syms(am::AllostericMechanism)
    num_A, den_A = _raw_symbolic_rate_polys_allosteric(am)
    a_only = _a_only_syms(am)
    rename_I = _a_to_i_rename(am)
    dep_A, _ = _dependent_param_exprs_allosteric(am)
    _add_case_b_renames!(rename_I, dep_A, am)
    den_i = _rename_symbols(_zero_symbols_in_poly(den_A, a_only), rename_I)
    S = _poly_param_syms(den_i)
    if !_i_state_dead(am)
        num_i = _rename_symbols(_zero_symbols_in_poly(num_A, a_only), rename_I)
        union!(S, _poly_param_syms(num_i))
    end
    S
end
```

(Add `_poly_param_syms(p::POLY) = Set{Symbol}(s for mono in keys(p) for (s,_) in mono)` if no equivalent exists.)

- [ ] **Step 5: Rewire `_build_dep_assignments`.** Replace the group-structure `i_names_set` construction (`:1499-1529`, including the `i_dead && tag !== :NonequalAI` gate and the `if !i_dead` synth block) with `i_names_set = _i_state_referenced_syms(am)`. Leave the catalytic dep-emission loop (`:1561-1571`) and its `a_only → 0` zeroing unchanged. Leave the reg-mirror and EqualAI-indep-mirror blocks (`:1537-1555`) driven by structure as today (do NOT gate on `S_I`).

- [ ] **Step 6: Rewire `_dependent_param_exprs`.** Gate `dep_I` and `indep_I_list` on `S_I` **with RHS closure**: an I-name is destructured if it is in `S_I` OR appears in the substituted RHS (`substitute_params_expr(v, rename_I)`) of any dep whose I-LHS is in `S_I`. Collect the RHS-referenced independent I-names in one pass and union them into the destructure set. (Reference logic: session scratchpad `verify_fix.jl`, `proposed_defined_vs_referenced` — verified 0/40 dangle with closure, 10/40 without.)

- [ ] **Step 7: Reconcile `_synthesized_dep_i_names`.** It early-returns `Symbol[]` when any `:OnlyA` group is present (`:1265-1266`). Change it to return the synth-dep I-names that are in `S_I` (so `parameters(Full)` enumerates exactly the emitted names), keeping the chokepoint routing.

- [ ] **Step 8: Run the §5a test — expect PASS.**

Run: the testset from Step 2.
Expected: PASS — all four mechanisms derive finite, transcript complete.

- [ ] **Step 9: Run the full suite — expect green (esp. perf + chokepoint AST-walker).**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS, including `test_rate_equation_performance` (0-alloc/sub-100 ns) and `test/test_types.jl:1577-1644`.

- [ ] **Step 10: Add permanent `MECHANISM_TEST_SPECS` fixtures.** In `test/mechanism_definitions_for_test_enzyme_derivation.jl`, add one `@allosteric_mechanism_src` spec per trigger path, transcribing each failing topology (the step lists are captured in the session scratchpad `investigate_kcat2.jl` output style — use `AllostericMechanism(em)`'s `steps`/`cat_allo_state` to read the topology, then write the DSL). Give each the correct `expected_n_independent_params`. If a full `analytical_rate_fn` is impractical to hand-derive, still register the spec (it inherits the derive/parameter/perf battery) and rely on the finite-rate assertion from Step 2.

- [ ] **Step 11: Run the full suite again; commit.**

```bash
git add src/rate_eq_derivation.jl test/test_rate_eq_derivation.jl test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -m "fix: define all inactive-state parameters via reference-polynomial S_I (§5a)"
```

---

### Task 2: §5b — Fix `_kcat_forward` "no components" + `NaN` at products=0

**Files:**
- Modify: `src/rate_eq_derivation.jl` (numerator/denominator derivation: `_allosteric_num_den_exprs` `:1581`, its `make_num_term`/`make_den_term`, and the free-enzyme-reference / conc-GCD reduction it relies on — exact site found in Step 1); possibly the shared reduction used by both allosteric and non-allosteric paths.
- Test: `test/test_rate_eq_derivation.jl`; `test/mechanism_definitions_for_test_enzyme_derivation.jl`.

**Interfaces:**
- No new public interface. `_kcat_forward(::AllostericEnzymeMechanism, params)` must return the finite peak forward turnover (products=0 limit) instead of erroring; `rate_equation` must be finite at products=0.

- [ ] **Step 1: Root-cause the uncancelled common factor.** Reconstruct one of the 28 mechanisms (session scratchpad `investigate_kcat2.jl`; e.g. `_kcat_forward`-failing, all-`EqualAI` + `:OnlyA`). Confirm: forward rate → finite nonzero as products→0⁺ but `NaN` at exactly 0; the denominator has **no constant `1` term** (both num and den carry a common product factor). Locate where the free-enzyme-per-segment reference (division-free derivation, PR #48) leaves the shared factor for these topologies — inspect `_raw_symbolic_rate_polys_allosteric` and the conc-GCD reduction it applies. Document the exact site in the commit message.

- [ ] **Step 2: Write the failing test.** In `test/test_rate_eq_derivation.jl`, for ≥1 reconstructed mechanism:

```julia
@testset "§5b kcat + finite rate at products=0 (regression)" begin
    em = Core.eval(EnzymeRates, Meta.parse(LDH_KCAT_FAILURE_MECH))()
    pnames = EnzymeRates.fitted_params(em); mets = EnzymeRates.metabolites(em)
    params = merge(NamedTuple{pnames}(ntuple(_->1.3, length(pnames))),
                   (Keq=20000.0, E_total=1.0))
    kc = EnzymeRates._kcat_forward(em, params)
    @test isfinite(kc) && kc > 0                      # currently: throws "no components"
    # rate at products=0, saturating substrates == the finite forward limit
    concs0 = NamedTuple{mets}(map(m -> (m in (:Lactate,:NAD)) ? 0.0 : 1e6, mets))
    @test isfinite(EnzymeRates.rate_equation(em, concs0, params))  # currently: NaN
    # kcat matches the numerical grid-peak (products→0⁺ limit)
    concs_eps = NamedTuple{mets}(map(m -> (m in (:Lactate,:NAD)) ? 1e-9 : 1e6, mets))
    @test isapprox(EnzymeRates.rate_equation(em, concs_eps, params), kc; rtol=1e-3)
end
```

- [ ] **Step 3: Run the test — expect fail** (`_kcat_forward` throws "no components"; `rate_equation` at products=0 is `NaN`).

- [ ] **Step 4: Implement the cancellation** at the site found in Step 1: divide the shared product factor out of numerator and denominator (extend the conc-GCD reduction / fix the free-enzyme reference) so the denominator regains its constant term. Do NOT touch the `isempty(a_keys)` guard directly — with the factor cancelled, the substrate-only pattern reappears and `_kcat_forward` works unchanged. If a residual guard is still reachable for a genuinely degenerate topology, keep it as an `error` (not a silent 0).

- [ ] **Step 5: Run the §5b test — expect PASS.**

- [ ] **Step 6: Run the full suite — expect green** (esp. perf; the cancellation reduces term count, so it cannot slow `rate_equation`). Confirm no currently-working mechanism's rendered equation changed (spot-check a few `MECHANISM_TEST_SPECS` `rate_equation_string`s against their committed expectations).

- [ ] **Step 7: Add a `MECHANISM_TEST_SPECS` fixture** for one of the 28 with an `analytical_kcat_fn` (the finite grid-peak value), plus a suite-wide guard asserting no valid mechanism's denominator lacks a constant term.

- [ ] **Step 8: Commit.**

```bash
git add src/rate_eq_derivation.jl test/test_rate_eq_derivation.jl test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -m "fix: cancel shared num/den product factor so kcat + products=0 rate are finite (§5b)"
```

---

### Task 3: §1 — Parsimony cutoff references all counts `< c`

**Files:**
- Modify: `src/identify_rate_equation.jl` — add `_parsimony_cutoff`; replace the inline expression at `:579-581`.
- Test: `test/test_identify_rate_equation.jl`.

**Interfaces:**
- Produces: `_parsimony_cutoff(best_loss_by_count::Dict{Int,Float64}, c::Int, loss_parsimony_threshold::Float64) -> Union{Nothing,Float64}`.

- [ ] **Step 1: Write the failing test.** In `test/test_identify_rate_equation.jl`:

```julia
@testset "§1 _parsimony_cutoff = threshold * min over all counts < c" begin
    f = EnzymeRates._parsimony_cutoff
    @test f(Dict(5=>0.02), 5, 1.01) === nothing            # no count < c
    @test f(Dict(5=>0.02,6=>0.05,7=>0.03), 8, 1.01) ≈ 1.01*0.02   # min over <c, not c-1
    @test f(Dict(5=>0.02), 7, 1.01) ≈ 1.01*0.02            # count gap: c-1=6 absent
    @test f(Dict(5=>0.01,6=>0.04), 7, 1.01) ≈ 1.01*0.01    # non-monotone → true min
end
```

- [ ] **Step 2: Run — expect fail** (`_parsimony_cutoff` undefined).

- [ ] **Step 3: Implement** (place near `_select_beam`):

```julia
# Parsimony reference = threshold × best loss over ALL counts strictly below c
# (not just c-1): an added parameter must beat the best simpler model of any size.
# Returns nothing when no simpler tier has been fit yet.
function _parsimony_cutoff(best_loss_by_count::Dict{Int,Float64}, c::Int,
                           loss_parsimony_threshold::Float64)
    prev = [best_loss_by_count[k] for k in keys(best_loss_by_count) if k < c]
    isempty(prev) && return nothing
    loss_parsimony_threshold * minimum(prev)
end
```

- [ ] **Step 4: Wire it in `_beam_search`.** Replace `:579-581`:

```julia
parsimony_cutoff = _parsimony_cutoff(best_loss_by_count, c, loss_parsimony_threshold)
```

- [ ] **Step 5: Run tests — expect PASS.**

- [ ] **Step 6: Commit.**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "fix: parsimony cutoff references best loss over all counts < c (§1)"
```

---

### Task 4: §3 — Progress message reports observed child param range

**Files:**
- Modify: `src/identify_rate_equation.jl:598-606`.
- Test: `test/test_identify_rate_equation.jl`.

**Interfaces:**
- Produces: `_np_range_label(child_entries) -> String` (`"n/a"` / `"8"` / `"7-13"`).

- [ ] **Step 1: Write the failing test:**

```julia
@testset "§3 child n_params range label" begin
    g = EnzymeRates._np_range_label
    @test g(NamedTuple[]) == "n/a"
    @test g([(n_params=8,)]) == "8"
    @test g([(n_params=7,), (n_params=13,), (n_params=9,)]) == "7-13"
end
```

- [ ] **Step 2: Run — expect fail.**

- [ ] **Step 3: Implement** the helper and use it in the progress line, replacing `target n_params=$target`:

```julia
_np_range_label(entries) = isempty(entries) ? "n/a" :
    let lo = minimum(e -> e.n_params, entries), hi = maximum(e -> e.n_params, entries)
        lo == hi ? string(lo) : "$lo-$hi"
    end
```

```julia
_progress(save_dir, show_progress,
    "Iteration $iteration (child n_params $(_np_range_label(child_entries))): " *
    "$(length(parents)) parents → $(length(children)) children | " *
    _batch_summary(child_entries, child_failures))
```

- [ ] **Step 4: Run tests — expect PASS. Step 5: Commit.**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "refactor: progress line shows observed child n_params range, not the iteration counter (§3)"
```

---

### Task 5: §2 — Fit-dedup by `eq_hash` + `fit_inherited` column + LOOCV dedup guard

**Files:**
- Modify: `src/identify_rate_equation.jl` — split `_process_batch` (`:447-482`) into compile/fit passes with a `memo`; thread `memo` from `_beam_search`; add `fit_inherited` to the row NamedTuple and `_rows_to_dataframe` (`:267-296`); expose the raw pre-rescale fit.
- Modify: `src/fitting.jl` — let `fit_rate_equation` (or a helper) return the raw pre-rescale params so the caller can re-rescale per mechanism.
- Test: `test/test_identify_rate_equation.jl`, `test/test_fitting.jl`.

**Interfaces:**
- Produces: `_process_batch` accepting/returning a `memo::Dict{UInt64,<rawfit>}`; each result row gains `fit_inherited::Bool`.
- Consumes: `eq_hash` (already computed via `_rate_eq_dedup_key`), `rescale_parameter_values` (per-mechanism), `fitted_params`.

- [ ] **Step 1: Write the failing dedup test** with a counting stub optimizer. In `test/test_identify_rate_equation.jl`, construct two structurally-distinct mechanisms that render the same `eq_hash` (reuse a known twin pair from the run, reconstructed from Sig strings), run them through `_process_batch` with a stub optimizer that increments a counter per `solve`, and assert the fit ran once, both rows carry equal `loss`, and each row's params equal a standalone per-mechanism rescale of the shared raw fit; assert `fit_inherited == [false, true]`.

- [ ] **Step 2: Run — expect fail** (no dedup; `fit_inherited` absent; fit runs twice).

- [ ] **Step 3: Expose the raw fit** in `src/fitting.jl`: add `fit_rate_equation_raw` (or return the pre-rescale params alongside) so the dedup layer holds the raw optimum and can call `rescale_parameter_values(m, raw; scale_k_to_kcat)` per mechanism. Keep `fit_rate_equation`'s existing behavior for other callers.

- [ ] **Step 4: Implement two-pass `_process_batch` + memo.** PASS-1 (`pmap`): compile + cap-check + render → `(mech, n, em_type, eq_text, eq_hash)` or `FitFailure`. PASS-2 (master): for each unseen `eq_hash`, recompile the representative and fit once (raw), store in `memo`; then build a `BatchEntry` for every mechanism from `memo[eq_hash]`, applying that mechanism's own `rescale_parameter_values`, setting `fit_inherited` (false for the representative, true otherwise). Thread `memo` from `_beam_search` so it persists across iterations. A representative whose fit throws → all its duplicates become `FitFailure`.

- [ ] **Step 5: Add `fit_inherited` to the row + DataFrame.** Extend the row NamedTuple (`:462-473`) and `_rows_to_dataframe` (`:277-285`) with `fit_inherited::Bool` (failure rows: `missing`/`false`).

- [ ] **Step 6: Run the dedup test — expect PASS.**

- [ ] **Step 7: Add the LOOCV `eq_hash`-uniqueness guard test** (§4): assert `_offer_cv!` and the `_cv_model_selection` candidate picker keep ≤1 row per `eq_hash` per param count. Add a one-line comment at `_rate_eq_dedup_key` noting its exact-render-only scope.

- [ ] **Step 8: Run the full suite — expect green.** Confirm selection on a small fixture is unchanged (loss is `eq_hash`-invariant).

- [ ] **Step 9: Commit.**

```bash
git add src/identify_rate_equation.jl src/fitting.jl test/test_identify_rate_equation.jl test/test_fitting.jl
git commit -m "perf: fit each distinct rate equation once (eq_hash memo) + fit_inherited column (§2)"
```

---

### Task 6: §6 — `maxtime` kwarg-forwarding test + report note

**Files:**
- Test: `test/test_fitting.jl` or `test/test_identify_rate_equation.jl`.
- Modify (doc): docstring of `identify_rate_equation` if the `maxtime` guidance is not already explicit.

- [ ] **Step 1: Write a test** asserting `maxtime` is forwarded to `Optimization.solve` (a stub optimizer that records the `maxtime` it received; call `fit_rate_equation(fp, stub; maxtime=1.23)` and assert `1.23` was seen). Run — expect it to pass if forwarding already works (it does per `fitting.jl:237`); this is a regression guard.

- [ ] **Step 2: Commit.**

```bash
git add test/test_fitting.jl
git commit -m "test: guard maxtime forwarding to Optimization.solve (§6)"
```

---

## Self-Review Notes

- **Spec coverage:** §5a→Task 1; §5b→Task 2; §1→Task 3; §3→Task 4; §2 + §4 (LOOCV dedup + canonical deferral)→Task 5; §6→Task 6. Task 0 handles branch + gitignore.
- **§4 canonical key** is intentionally NOT a task (deferred by decision); only the existing-LOOCV-dedup guard test (Task 5 Step 7) is in scope.
- **§5b Step 1 is investigation-first** by necessity (exact cancellation site unpinned); the failing test in Step 2 is the concrete gate. If the root-cause fix proves intractable in one pass, STOP and report to Denis rather than shipping a wrong kcat — do not fall back to `return 0.0`.
- **Reviewer gates:** each Task ends at an independently testable, committable deliverable.
