# Allosteric D1 — State-Parameterized I-State Re-Derivation — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the allosteric I-state "rename-the-output" reconstruction with per-state "tag-the-input" re-derivation on the shared King–Altman kernel, then adopt per-monomer normalization (drop `CatN`, `E_total` = active-site concentration).

**Architecture:** Derive each MWC conformational state (A, I) by running the *existing* non-allosteric derivation pipeline on a per-state view of the catalytic mechanism whose `Parameter`s carry the state tag (`:A`/`:I`/`:EqualAI`) and whose `:OnlyA` steps are pruned for the I-state. `_state_tag(:EqualAI) == ""` makes EqualAI sharing automatic. A single local naming rule handles the one irreducible Case-B (a shared EqualAI dependent whose Haldane references a NonequalAI symbol). The MWC combine (`num_A + L·num_I` over `den_A + L·den_I`), reg-site partition factors, and `L` are unchanged.

**Tech Stack:** Julia 1.12; `@generated` rate-equation derivation; `Rational{BigInt}` symbolic polynomial backend; the package's own `MECHANISM_TEST_SPECS` harness.

**Spec:** `docs/superpowers/specs/2026-07-02-allosteric-state-parameterized-derivation-design.md` (read it fully before starting; §2, §3, §4, §4a, §5 are load-bearing).

## Global Constraints

- **0-allocation / sub-100 ns `rate_equation` contract is a HARD GATE** for every mechanism in `MECHANISM_TEST_SPECS`. Enforced by `test_rate_equation_performance` in `test/test_rate_eq_derivation.jl` (`allocs == 0`, `t < 100e-9`). Any task that regresses this must STOP and escalate, not work around it. (Spec §8.)
- **Parameter naming chokepoint guard:** all `Parameter → Symbol` rendering flows through `name(p, m)`; no stray `Symbol("K…")`/`k…`/`V…`/`L…` literals. Enforced by the AST-walker at `test/test_types.jl:1577-1644`.
- **Canonical Step Form is load-bearing:** do not relax step/group/direction canonicalization in constructors.
- **Byte-identical golden reference (Commit 1):** every allosteric mechanism that exists in `MECHANISM_TEST_SPECS` today must produce identical `rate_equation_string(m, Reduced)`, `parameters(m, Full)`, and `parameters(m, Reduced)` after the refactor. This is the primary safety net.
- **Test-suite ops** (memory `project-test-suite-ops`): the full suite (`julia --project -e 'using Pkg; Pkg.test()'`) is ~11 min and memory-heavy — run only ONE at a time; before re-running, `pgrep -af runtests.jl` and `kill -9` orphans. For fast iteration use the TestEnv focused-run recipe (below). A `@test` failure does NOT abort an `include`, so grep output for a nonzero `Fail`/`Error` column, not just for absence of `ERROR`.
- **92-char line limit, 4-space indent. ABOUTME: header on any new file.** Commit after every green task; never commit a red tree.

**TestEnv focused-run recipe** (fast iteration on one file):
```julia
using TestEnv; TestEnv.activate()
using EnzymeRates, Test, DataFrames, CSV, Random, Statistics, LinearAlgebra,
      OptimizationCMAEvolutionStrategy, OptimizationBBO
include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
include("test/test_allosteric_golden.jl")   # or the file under test
```
Run headless: `julia --project=. /path/to/focused_script.jl 2>&1 | tee out.txt`, then grep for `Fail`/`Error`.

---

## File Structure

- **Modify** `src/rate_eq_derivation.jl` — the allosteric derivation surface (~lines 1040–1738). New per-state helpers replace the rename/reconstruction functions.
- **Modify** `src/thermodynamic_constr_for_rate_eq_derivation.jl` — parameterize `_dependent_param_exprs_kernel` (and, if needed, `_raw_symbolic_rate_polys`) to accept caller-supplied state-tagged `step_params`. Non-allosteric callers unchanged (default args).
- **Create** `test/test_allosteric_golden.jl` — byte-identical golden reference test.
- **Create** `test/reference/allosteric_golden_reference.txt` — committed golden data (generated from current `main` before any refactor).
- **Modify** `test/mechanism_definitions_for_test_enzyme_derivation.jl` — add `:OnlyA`-catalytic fixtures and the degeneracy characterization fixture.
- **Modify** `test/runtests.jl` — include the new golden test file.

---

## COMMIT 1 — Structure-preserving re-derivation (keeps `CatN`)

### Task 1: Golden-reference safety net

**Files:**
- Create: `test/test_allosteric_golden.jl`
- Create: `test/reference/allosteric_golden_reference.txt`
- Modify: `test/runtests.jl` (add include)

**Interfaces:**
- Produces: `_allosteric_golden_lines()` returning `Vector{String}` — the canonical serialization used by both the generator and the comparison test.

- [ ] **Step 1: Write the golden test file.** Create `test/test_allosteric_golden.jl`:

```julia
# ABOUTME: Byte-identical golden reference for allosteric rate-equation strings
# ABOUTME: and parameter lists; guards the D1 state-parameterized re-derivation.
using Test

const _ALLO_GOLDEN_PATH =
    joinpath(@__DIR__, "reference", "allosteric_golden_reference.txt")

"""Canonical serialization of every allosteric spec's derivation output."""
function _allosteric_golden_lines()
    lines = String[]
    for spec in MECHANISM_TEST_SPECS
        spec.mechanism isa EnzymeRates.AllostericEnzymeMechanism || continue
        m = spec.mechanism
        push!(lines, "### " * spec.name)
        push!(lines, "REDUCED_STRING " *
              EnzymeRates.rate_equation_string(m, EnzymeRates.Reduced))
        push!(lines, "PARAMS_FULL " * string(parameters(m, EnzymeRates.Full)))
        push!(lines, "PARAMS_REDUCED " *
              string(parameters(m, EnzymeRates.Reduced)))
    end
    lines
end

@testset "allosteric golden reference (D1)" begin
    @test isfile(_ALLO_GOLDEN_PATH)
    current = _allosteric_golden_lines()
    reference = readlines(_ALLO_GOLDEN_PATH)
    @test length(current) == length(reference)
    for (c, r) in zip(current, reference)
        @test c == r
    end
end
```

- [ ] **Step 2: Generate the reference from current `main` behavior.** In a focused-run REPL/script (recipe above), after including the definitions file and the golden test file, run:

```julia
mkpath(joinpath(@__DIR__, "reference"))  # or the absolute test/reference path
open(_ALLO_GOLDEN_PATH, "w") do io
    for l in _allosteric_golden_lines()
        println(io, l)
    end
end
```

Do this **before** any derivation change so it captures today's exact output (including current `1^2`/`CatN` artifacts — those are preserved through Commit 1).

- [ ] **Step 3: Add the include to `test/runtests.jl`** after `test_rate_eq_derivation.jl` (line 13):

```julia
    include("test_allosteric_golden.jl")
```

- [ ] **Step 4: Run the golden test — must PASS on unmodified code.**

Focused run including `test_allosteric_golden.jl`. Expected: the `allosteric golden reference (D1)` testset passes (0 fails). This confirms the net matches current behavior.

- [ ] **Step 5: Commit.**

```bash
git add test/test_allosteric_golden.jl test/reference/allosteric_golden_reference.txt test/runtests.jl
git commit -m "test: golden reference for allosteric rate eqs (D1 safety net)"
```

---

### Task 2: State-tagged per-state derivation helpers (A-state wired, byte-identical)

**Files:**
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl` (parameterize the kernel)
- Modify: `src/rate_eq_derivation.jl` (new helpers; route A-state through them)
- Test: `test/test_allosteric_golden.jl` (unchanged; must stay green)

**Interfaces:**
- Produces:
  - `_state_step_params(am::AllostericMechanism, state::Symbol)` → the same shape as `_step_parameters(::Mechanism)` (a per-flat-step tuple of `Parameter`s) but with each catalytic group's `Parameter`s tagged: `:NonequalAI`→`state`, `:EqualAI`→`:EqualAI`, `:OnlyA`→`state`. For `state == :I`, steps belonging to `:OnlyA` groups are pruned.
  - `_state_mechanism(am::AllostericMechanism, state::Symbol)` → `Mechanism` — the catalytic mechanism for the state (full for `:A`; `:OnlyA` steps removed for `:I`).
  - `_state_rate_polys(am, state)` → `(num_poly, den_poly)` via `_raw_symbolic_rate_polys(_state_mechanism(am,state), _state_step_params(am,state), Dict{Symbol,Symbol}(), subs_syms, prods_syms)`.
- Consumes: `_raw_symbolic_rate_polys(mech, step_params, rename_map, subs, prods)` (5-arg, `rate_eq_derivation.jl:343`); `_dependent_param_exprs_kernel` (kernel to be parameterized).

- [ ] **Step 1: Parameterize the kernel by `step_params`.** In `thermodynamic_constr_for_rate_eq_derivation.jl`, change `_dependent_param_exprs_kernel(mech, rename)` (line 299) to accept a caller-supplied `step_params` (default = `_step_parameters(mech)`), so its internal `step_params = _step_parameters(mech)` (line 320) becomes the argument. All symbol rendering already flows through `name(p, mech)` (line 321), so state-tagged params produce state-tagged symbols. Non-allosteric callers pass nothing → unchanged behavior.

Add a keyword or positional overload — match the existing dispatch style at `:422-424`. Confirm the plain-`Mechanism` type-dispatch wrapper still compiles.

- [ ] **Step 2: Run the full suite — non-allosteric path unaffected.** Confirm `test_rate_eq_derivation.jl` (non-allosteric mechanisms) is still green after the signature change. (Focused run of that file.) Expected: no new fails.

- [ ] **Step 3: Add the per-state helpers** in `rate_eq_derivation.jl` near the existing allosteric block (`_A_rename_parameters` region, ~1129). Implement `_state_step_params`, `_state_mechanism`, `_state_rate_polys` per the Interfaces above. For `_state_step_params`, walk `steps(am)` with `_group_rep`/`_step_parameters` exactly as the current `_A_rename_parameters` walks groups, but *construct* the `Parameter`s with the state tag instead of building a rename map.

- [ ] **Step 4: Route the A-state through `_state_rate_polys(am, :A)`.** In `_allosteric_num_den_exprs` (`:1582`), replace the `_raw_symbolic_rate_polys_allosteric(am)` call (which derives at `:None` then renames to `:A`) with `_state_rate_polys(am, :A)`. Leave the I-state path (`_i_state_num_den_polys`) untouched for now.

- [ ] **Step 5: Run the golden test — MUST stay byte-identical.** Focused run. Expected: `allosteric golden reference (D1)` passes with 0 fails. If any spec's `REDUCED_STRING` differs, the A-state re-derivation is not reproducing the `:None`→`:A` rename — diff the failing spec and reconcile before proceeding. Do NOT regenerate the reference to make it pass.

- [ ] **Step 6: Run the perf gate.** Focused run of `test_rate_equation_performance` in `test_rate_eq_derivation.jl`. Expected: `allocs == 0`, `t < 100e-9` for all specs.

- [ ] **Step 7: Commit.**

```bash
git add src/thermodynamic_constr_for_rate_eq_derivation.jl src/rate_eq_derivation.jl
git commit -m "refactor: derive allosteric A-state via state-tagged step_params"
```

---

### Task 3: Re-derive the I-state on the pruned/tagged graph

**Files:**
- Modify: `src/rate_eq_derivation.jl`
- Test: `test/test_allosteric_golden.jl` (must stay green)

**Interfaces:**
- Produces:
  - `_state_dependent_exprs(am, state)` → `(dep_exprs::Dict{Symbol,Union{Symbol,Expr}}, indep::Vector{Symbol})` — the per-state Haldane/Wegscheider assignments in native state names, via the parameterized kernel on `_state_mechanism(am,state)` with `_state_step_params(am,state)`. For `state == :I`, applies the **one-rule Case-B naming** (below).
  - Case-B rule: after the I-run, any dependent whose LHS renders as a bare/`:EqualAI` symbol but whose RHS references a state-differing (`:I`-tagged NonequalAI) symbol takes its forced-`:I` name (reuse `_force_inactive`/`_dep_inactive_name`'s core so the AST guard is satisfied). Spec §4, §4a. Because Gaussian elimination expresses each dependent purely in terms of independents, this check is one non-transitive pass.

- [ ] **Step 1: Implement `_state_dependent_exprs`** using the parameterized kernel from Task 2. For `:I`, apply the Case-B renaming pass.

- [ ] **Step 2: Replace `_i_state_num_den_polys(am)`** usage in `_allosteric_num_den_exprs` (`:1601`) with `_state_rate_polys(am, :I)`. Keep the existing `_i_state_dead(m)` gate (`:1602`) — but note that with `:OnlyA` steps pruned from the I-graph, `_state_rate_polys(am, :I)` should now yield `N_I` that is *natively* zero when the pruned graph has no productive cycle (verify: the reaction-cut numerator is 0 when substrate and product sub-graphs are disconnected). Keep `_i_state_dead` as a guard/assertion for now; do not delete it in this task.

- [ ] **Step 3: Run the golden test — MUST stay byte-identical** for all existing specs (none of which has an `:OnlyA` *catalytic* group, so the I-graph equals the A-graph modulo tags; re-derivation must reproduce the current rename output exactly). Focused run. Expected: 0 fails. Diff and reconcile any mismatch — common cause is Case-B naming not matching today's `_synthesized_dep_i_names` output.

- [ ] **Step 4: Run the perf gate** (as Task 2 Step 6). Expected: green.

- [ ] **Step 5: Commit.**

```bash
git add src/rate_eq_derivation.jl
git commit -m "refactor: re-derive allosteric I-state natively via tagged graph"
```

---

### Task 4: Reroute dep-assignments, `_dependent_param_exprs(::Allosteric)`, `_kcat_forward`, `parameters(Full)`

**Files:**
- Modify: `src/rate_eq_derivation.jl`
- Test: `test/test_allosteric_golden.jl` + full suite

**Interfaces:**
- Consumes: `_state_rate_polys`, `_state_dependent_exprs` (Tasks 2–3).
- Produces: `_build_dep_assignments(M)`, `_dependent_param_exprs(::AllostericEnzymeMechanism)`, `_kcat_forward(::AllostericEnzymeMechanism)`, and the allosteric branch of `parameters(m, Full)` all sourced from the per-state helpers (A-run + I-run), not from rename reconstruction.

- [ ] **Step 1: Rebuild `_build_dep_assignments`** (`:1504`) from `_state_dependent_exprs(am, :A)` and `_state_dependent_exprs(am, :I)` — emit `a_assignments` from the A-run and `i_assignments` from the I-run directly (each already in native state names). The `:OnlyA`-referencing `= 0` case is subsumed: pruned steps produce no assignment.

- [ ] **Step 2: Rebuild `_dependent_param_exprs(::AllostericEnzymeMechanism)`** (`:1358`) as the union of the A-run and I-run independent sets plus reg-site `Kreg`s and `:L`, sorted for content-canonical order. Drop the Pass-1/Pass-2 rename closure.

- [ ] **Step 3: Reroute `_kcat_forward(::AllostericEnzymeMechanism)`** (`:832`) to consume `_state_rate_polys`/`_state_dependent_exprs` instead of re-inlining I-state poly construction (`:846-859`). The reg-corner assembly stays; only the num/den/dep source changes.

- [ ] **Step 4: Reroute the allosteric `parameters(m, Full)` branch** (`:58-86`) to enumerate from the per-state derivations (A-run independents+deps, I-run independents+deps, reg `Kreg`s, `L`), removing the `_synthesized_dep_i_names` + `filter!`/`splice!` reconciliation. The output ORDER must match the current golden `PARAMS_FULL`.

- [ ] **Step 5: Run the golden test — MUST stay byte-identical** (`REDUCED_STRING`, `PARAMS_FULL`, `PARAMS_REDUCED` for all existing specs). Focused run. Expected: 0 fails.

- [ ] **Step 6: Run the FULL suite** (single run; kill orphans first). This exercises `test_accessors.jl` (parameters), `test_rate_eq_derivation.jl` (kcat via `analytical_kcat_fn`), `test_types.jl` (naming guard), fitting, enumeration. Expected: 0 fails / 0 errors. Grep the output for the `Fail`/`Error` columns.

- [ ] **Step 7: Commit.**

```bash
git add src/rate_eq_derivation.jl
git commit -m "refactor: source allosteric deps/kcat/Full-params from per-state derivation"
```

---

### Task 5: Delete the now-dead reconstruction functions

**Files:**
- Modify: `src/rate_eq_derivation.jl`
- Test: full suite

**Interfaces:**
- Removes (only if grep-confirmed unreferenced): `_A_rename_parameters`, `_a_to_i_rename`, `_I_rename_parameters`, `_raw_symbolic_rate_polys_allosteric`, `_dependent_param_exprs_allosteric`, `_synthesized_dep_i_names`, `_add_case_b_renames!`, `_i_state_referenced_syms`, `_all_i_state_parameters` (if now unused). Keep the minimal `_force_inactive`/`_dep_inactive_name` naming helper used by the Case-B rule.
- **RETAINED (do NOT delete — updated after Task 3):** `_i_state_num_den_polys` and `_i_state_dead` are load-bearing for the dead-state reachability partition (spec §3), as is whatever A→I rename helper `_i_state_num_den_polys` still calls to name its surviving binding groups (likely keeps `_onlyA_parameters`/`_a_only_syms` and part of the rename machinery alive). The grep-then-delete-only-if-unreferenced procedure is self-correcting, so this is an expectation-setting note, not a rewrite: the net deletion is smaller than the header's original estimate.

- [ ] **Step 1: For each candidate function, grep the whole `src/` and `test/` for references.**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
for f in _A_rename_parameters _a_to_i_rename _I_rename_parameters \
  _raw_symbolic_rate_polys_allosteric _dependent_param_exprs_allosteric \
  _i_state_num_den_polys _synthesized_dep_i_names _add_case_b_renames! \
  _i_state_referenced_syms _all_i_state_parameters; do
  echo "== $f =="; grep -rn "$f" src/ test/;
done
```

- [ ] **Step 2: Delete each function whose only remaining reference is its own definition.** Do NOT delete anything still referenced (e.g. `_force_inactive` if the Case-B rule uses it, or `_i_state_dead` if kept as an assertion). Re-run the grep after deletion to confirm zero dangling references.

- [ ] **Step 3: Run the golden test + FULL suite.** Expected: 0 fails / 0 errors.

- [ ] **Step 4: Run the perf gate.** Expected: `allocs == 0`, `t < 100e-9`.

- [ ] **Step 5: Re-read the changed region of `rate_eq_derivation.jl`** for dead code and further simplification (per repo Code Style). Remove any now-orphaned local helpers.

- [ ] **Step 6: Commit.**

```bash
git add src/rate_eq_derivation.jl
git commit -m "refactor: delete allosteric I-state rename reconstruction (~250-400 LOC)"
```

---

### Task 6: `:OnlyA`-catalytic fixtures + degeneracy characterization fixture (point 3 + §4a)

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`
- Test: the new fixtures via `test_rate_eq_derivation.jl`'s spec-driven checks

**Interfaces:**
- Consumes: `@allosteric_mechanism_src` / `@enzyme_mechanism_src` fixture macros; `MechanismTestSpec`.

- [ ] **Step 1: Add Fixture A — dimer uni-uni, `:OnlyA` catalysis (I-state dead by re-derivation).** Mechanism: `E + S ⇌ E(S)` (binding, `:EqualAI`), `E(S) <--> E(P)` (SS catalysis, `:OnlyA`), `E(P) ⇌ E + P` (release, `:EqualAI`), `catalytic_multiplicity: 2`, no regulators. Add it via `@allosteric_mechanism_src`. Then print `parameters(m, Reduced)` (focused run) to read the exact canonical symbol names, and write `analytical_rate_fn` with the math below, destructuring the confirmed names:

```julia
# math (align symbol names to parameters(m, Reduced)):
#   k_A_rev derived from A-state Haldane; N_I = 0 (catalysis pruned in I).
function onlyA_uniuni_rate(params, concs)
    (; K_S_E, K_P_E, k_A_ES_to_EP, L, Keq, E_total) = params  # confirm names
    (; S, P) = concs
    k_A_EP_to_ES = (1/Keq) * K_P_E * (1/K_S_E) * k_A_ES_to_EP
    Q_A = 1 + S/K_S_E + P/K_P_E
    Q_I = 1 + S/K_S_E + P/K_P_E
    N_A = k_A_ES_to_EP * S/K_S_E - k_A_EP_to_ES * P/K_P_E
    return E_total * 2.0 * (N_A * Q_A) / (Q_A^2 + L * Q_I^2)
end
```
Set `analytical_kcat_fn = p -> 2 * p.k_A_ES_to_EP` (per-oligomer; Commit 2 drops the `2`). Set `run_ode_test=false`.

- [ ] **Step 2: Run the fixture's analytical check.** Focused run of `test_rate_eq_derivation.jl` restricted to the new spec (or the whole file). Expected: the derived `rate_equation` matches `analytical_rate_fn` within `reference_rtol`. This proves the re-derived broken-cycle I-state (`N_I = 0`) is correct without the old `_i_state_dead` patch.

- [ ] **Step 3: Add Fixture B — `:OnlyA` on a redundant binding path (I-state still functions).** Random-order bi-uni where one substrate's binding-first path is `:OnlyA` and the other order remains, so the I-state retains flux via the surviving path. Validate with **properties** (safer than a hand-formula): (a) equilibrium flux is zero — pick concentrations at the `Keq` ratio and assert the derived rate ≈ 0; (b) the `:OnlyA` step's rate constant symbol does NOT appear in `rate_equation_string(m, Reduced)`'s I-state half but DOES appear in the A-state half; (c) a numeric spot-check at one non-equilibrium point against an independent King–Altman hand-calculation. Write these as explicit `@test`s in a dedicated `@testset` in `test_rate_eq_derivation.jl`.

- [ ] **Step 4: Add the degeneracy characterization fixture (§4a).** The existing "Random-order Bi-Bi" mechanism made allosteric with ONE binding group `:NonequalAI`, the rest `:EqualAI`. Assert: (a) equilibrium flux is zero; (b) capture the current `parameters(m, Reduced)` / `rate_equation_string` output as the documented behavior. Add a comment: this is the phantom-parameter config that D3 (`2026-05-29-nonequalai-rank-validity.md`) will reject; D1 only characterizes it.

- [ ] **Step 5: Run the FULL suite.** Expected: 0 fails / 0 errors. (New specs also feed the golden test's `_allosteric_golden_lines` — regenerate the golden reference to INCLUDE the new specs, since they are legitimately new, and re-commit the reference. This is the ONE sanctioned reference regeneration in Commit 1, and only adds the new specs' lines; existing specs' lines must be unchanged — verify via `git diff` that only additions appear.)

- [ ] **Step 6: Commit.**

```bash
git add test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_rate_eq_derivation.jl test/reference/allosteric_golden_reference.txt
git commit -m "test: OnlyA-catalytic + degeneracy fixtures for allosteric re-derivation"
```

---

## COMMIT 2 — Per-monomer normalization (drops `CatN`)

### Task 7: Drop `CatN`, `E_total` = active-site concentration

**Files:**
- Modify: `src/rate_eq_derivation.jl` (num/den assembly)
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (analytical formulas + kcat)
- Modify: `test/reference/allosteric_golden_reference.txt` (regenerate)
- Test: full suite

**Interfaces:**
- Changes `_allosteric_num_den_exprs` (`:1582`) so the numerator drops the leading `CatN *` factor: numerator becomes `num_A + L·num_I` (dead-I: `num_A`), unchanged `Q^(CatN-1)` / `Q^CatN` powers.

- [ ] **Step 1: Write the failing test first.** Update the analytical formulas to per-active-site: for every allosteric spec with an `analytical_rate_fn`, remove the leading `CatN` multiplier (`2 *`, `4 *`) and change `analytical_kcat_fn` from `p -> CatN * p.k…` to `p -> p.k…`. (PK: `4.0 * (num_R + L*num_T)` → `(num_R + L*num_T)`, `analytical_kcat_fn = p -> p.k5f`; MWC Dimer: drop `2 *`; Fixture A: `2.0 *` → `1.0 *`, kcat `2*` → bare.) Run the affected specs — they now FAIL against the still-`CatN`-scaled derivation.

- [ ] **Step 2: Run to confirm failure.** Focused run. Expected: analytical-rate mismatch (factor of `CatN`) for the updated specs.

- [ ] **Step 3: Drop the `CatN *` prefactor** in `_allosteric_num_den_exprs` — remove it from both the dead-I branch (`:1636`) and the live branch (`:1639`), so the returned numerator has no leading multiplier. Leave the `Q^(CatN-1)`/`Q^CatN` binding-statistics powers untouched.

- [ ] **Step 4: Run the updated analytical checks — now PASS.** Focused run. Expected: derived rate matches the per-active-site formulas.

- [ ] **Step 5: Regenerate the golden reference** (the `REDUCED_STRING`s lose their leading `2 *`/`4 *`). Regenerate `test/reference/allosteric_golden_reference.txt` from the new output; `git diff` it to confirm the ONLY changes are the dropped leading coefficients (the `1^n` artifacts and everything else are unchanged — those are D4a, out of scope).

- [ ] **Step 6: Document `E_total`.** Update the doc-comment/docstring where `E_total`/`Etot` is defined (and any user-facing docs that describe it) to state it is the **active-site (protomer) concentration**, so allosteric `kcat` is per active site. Grep for existing `E_total` documentation and keep the wording consistent.

- [ ] **Step 7: Run the FULL suite + perf gate.** Expected: 0 fails / 0 errors; `allocs == 0`, `t < 100e-9`.

- [ ] **Step 8: Commit.**

```bash
git add src/rate_eq_derivation.jl test/mechanism_definitions_for_test_enzyme_derivation.jl test/reference/allosteric_golden_reference.txt
git commit -m "feat: per-monomer allosteric normalization (drop CatN, E_total=active site)"
```

---

## Self-Review (completed by plan author)

**Spec coverage:** §1 root cause → Tasks 2–5. §2 tag-the-input → Tasks 2–3. §3 EqualAI/NonequalAI/OnlyA/shared-Keq → Tasks 2–3, Fixture A/B. §4 Case-B one rule → Task 3. §4a degeneracy boundary → Task 6 Step 4 characterization fixture. §5 per-monomer → Task 7. §6 delete-vs-stay → Task 5. §7 two-commit sequencing → Commit 1 / Commit 2 split. §8 test strategy → Task 1 golden, Task 6 fixtures, perf gate in every task. §9 open-verification (Wegscheider, `parameters(Full)` order, kcat reg-corners) → Tasks 4, 6.

**Placeholder scan:** Fixture B and the degeneracy fixture use property + numeric-spot-check tests rather than a hand-derived closed form (deliberate — a wrong hand-derivation is worse than a property test). Fixture A gives the full closed form with an explicit "confirm names against `parameters(m, Reduced)`" step because exact rendered symbol names (`E(S)`→`ES` etc.) must be read from the code, not guessed.

**Type consistency:** `_state_step_params`/`_state_mechanism`/`_state_rate_polys`/`_state_dependent_exprs` names are used consistently across Tasks 2–5. The parameterized kernel keeps its existing return type `(dep_exprs, indep)`.

**Known risk:** the "byte-identical" criterion assumes re-derivation reproduces the current rename output's exact Expr shape (term order, binary-tree association). If a spec's string differs only by benign reordering, STOP and escalate to Denis rather than regenerating the golden reference — a reorder means the assembly order changed and must be understood, not masked.
