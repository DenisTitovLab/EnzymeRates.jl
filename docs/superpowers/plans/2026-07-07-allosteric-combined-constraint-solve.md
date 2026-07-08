# Allosteric Combined Constraint Solve — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the fragmented allosteric dependent-parameter derivation (two per-state solves + `_split_resolution` + collapse mirrors + ad-hoc merge) with one combined constraint solve, deleting far more logic than it adds.

**Architecture:** Assemble the A-state and I-state thermodynamic-constraint rows into one rational matrix over state-tagged parameter columns and run the existing priority-pivoting Gaussian elimination once. Cross-state ties emerge as `I-row − A-row` (the Keq cancels), so `_split_resolution`, `_collapse_mirror_exprs`, the `S_I` gate, the per-state merge, and `#61`'s filter all disappear. Both consumers — `_dependent_param_exprs(::AllostericEnzymeMechanism)` (param list) and `_build_dep_assignments` (body constraint lines) — read the same solve.

**Tech Stack:** Julia; `Rational{BigInt}` linear algebra; `@generated` rate-equation derivation.

**Spec:** `docs/superpowers/specs/2026-07-07-allosteric-combined-constraint-solve-design.md`.

## Global Constraints

- **Ground-truth gates are NEVER re-baselined.** The QSSA oracle (`test_reference_qssa`), the equilibrium test (`test_haldane_equilibrium`, `v=0` at `Q=Keq`, tol `1e-10`), the ODE oracle (`test_ode_steadystate`), and the hand-written `analytical_rate_fn` checks are derivation-independent. If any fails, the code is wrong — fix the code, never the expected value.
- **Only these get re-baselined:** the byte-identical golden `test/reference/allosteric_golden_reference.txt`, and the tests of deleted internals (`test/test_split_resolution.jl`, mirror-text asserts in `test/test_allosteric_collapse.jl`). Re-baseline a golden ONLY after the ground-truth gates for that mechanism are green.
- **Non-allosteric path stays byte-identical.** The kernel factoring must not change any non-allosteric output. The non-allosteric oracle/golden tests must stay green with no edits.
- **`rate_equation` performance contract is sacred:** `allocs == 0` and `t < 120e-9` for every `MECHANISM_TEST_SPECS` entry (`test_performance`). Derivation runs at compile time, so this should be untouched — but verify it after Phase 2.
- **Parameter-name chokepoint:** every `K`/`k`/`V`/`L` symbol must be built through a `name(::Parameter…)` renderer. The AST guard at `test/test_types.jl:1694-1707` fails the build on any stray `Symbol("K…")` literal in `src/`.
- **Measure of success = net logic removed** (not comments/docs). The combined solve must stay very well documented, especially *why* stacked A/I rows produce the cross-state ties.
- **Run a single derivation-test file** with the fixture loaded:
  `julia --project=. -e 'using Test, EnzymeRates, LinearAlgebra, Random, OrdinaryDiffEqFIRK; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_rate_eq_derivation.jl")'`
  Full suite: `julia --project=. -e 'using Pkg; Pkg.test()'`.
- **Commit frequently**; the branch is `allosteric-combined-constraint-solve`.

---

## Phase 0 — Capture the bug as failing tests

The enumerated mechanisms that fail in production escaped `MECHANISM_TEST_SPECS` (the suite is green). Add them as regressions that fail on HEAD before any fix.

### Task 0.1: Store the three failure reproducers as fixtures

**Files:**
- Create: `test/reference/allosteric_undefvar_reproducers.jl` (the three `AllostericEnzymeMechanism` type literals, one per defect family: `koff` EqualAI-shared, `K_I` NonequalAI, `kon_I` SS-speed).
- Source of truth: the saved reproducers in the session scratchpad (`err_mech_koff_Pyruvate_ENAD.txt`, `err_mech_K_NAD_ELactate.txt`, and a `kon_I_NAD_E` case — grep the run CSVs' error column for one).

**Interfaces:**
- Produces: `const ALLOSTERIC_UNDEFVAR_REPRODUCERS::Vector{Type}` — the three lifted `AllostericEnzymeMechanism{…}` types, ready to instantiate with `T()`.

- [ ] **Step 1:** Recover the three type strings. Two are saved:
  `cat /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/24f0cbd6-a595-4254-91cd-de42c5cf165b/scratchpad/err_mech_koff_Pyruvate_ENAD.txt` and `…/err_mech_K_NAD_ELactate.txt`. For the `kon_I` family, extract one from the run:
  ```bash
  cd /home/denis.linux/.julia/dev/EnzymeRates/docs/ldh_hpc_results/2026_07_07_results
  python3 -c "import csv,sys; csv.field_size_limit(sys.maxsize)
  [print(r[2]) or sys.exit() for r in csv.reader(open('equation_search_iteration_10.csv')) if 'kon_I_NAD_E' in r[5]]"
  ```
- [ ] **Step 2:** Write `test/reference/allosteric_undefvar_reproducers.jl` defining `const ALLOSTERIC_UNDEFVAR_REPRODUCERS = [ eval(Meta.parse(s)) for s in (…three strings…) ]`. Keep the raw strings inline (they are the fixture). Add the two `ABOUTME:` header lines.
- [ ] **Step 3:** Sanity-run: instantiate each and confirm it is an `AllostericEnzymeMechanism`.
  `julia --project=. -e 'using EnzymeRates; include("test/reference/allosteric_undefvar_reproducers.jl"); @assert all(T -> T() isa EnzymeRates.AllostericEnzymeMechanism, ALLOSTERIC_UNDEFVAR_REPRODUCERS); println("ok ", length(ALLOSTERIC_UNDEFVAR_REPRODUCERS))'`
  Expected: `ok 3`.
- [ ] **Step 4: Commit.** `git add test/reference/allosteric_undefvar_reproducers.jl && git commit -m "test: capture the three allosteric UndefVar reproducers as fixtures"`

### Task 0.2: Add the structural-invariant regression test (fails on HEAD)

**Files:**
- Modify: `test/test_rate_eq_derivation.jl` — extend the existing `indep ∩ keys(dep) == ∅` block (around :1171-1189) with an acyclicity + all-RHS-defined check, and add a testset over `ALLOSTERIC_UNDEFVAR_REPRODUCERS`.

**Interfaces:**
- Consumes: `EnzymeRates._dependent_param_exprs(T)` → `(dep::Dict{Symbol,Union{Symbol,Expr}}, indep::Tuple)`.
- Consumes: `ALLOSTERIC_UNDEFVAR_REPRODUCERS` (Task 0.1).

- [ ] **Step 1: Write the failing test.** Add a helper that verifies the dependent-assignment graph is a DAG with every RHS symbol defined, then assert it over the reproducers plus every allosteric spec:

```julia
# every dep RHS symbol must be a fitted param or a PRECEDING dep (acyclic), and
# indep ∩ keys(dep) == ∅. A single combined solve guarantees this; the buggy
# merge does not.
function _dep_graph_is_sound(T)
    dep, indep = EnzymeRates._dependent_param_exprs(T)
    indepset = Set(indep)
    isempty(intersect(indepset, keys(dep))) || return false
    rhs_syms(x) = x isa Symbol ? [x] :
        x isa Expr ? reduce(vcat, map(rhs_syms, x.args); init=Symbol[]) : Symbol[]
    # topological soundness: no dep may (transitively) depend on itself
    depset = Set(keys(dep))
    edges = Dict(k => Symbol[s for s in rhs_syms(v) if s in depset] for (k,v) in dep)
    state = Dict{Symbol,Int}()  # 0=unseen 1=onstack 2=done
    ok = true
    function dfs(n)
        get(state,n,0) == 2 && return
        get(state,n,0) == 1 && (ok = false; return)
        state[n] = 1
        for m in get(edges,n,Symbol[]); dfs(m); ok || return; end
        state[n] = 2
    end
    for k in keys(dep); dfs(k); ok || break; end
    # every RHS symbol is either fitted, a dep, or a global (Keq/E_total/L/metabolite)
    known = union(indepset, depset)
    for (_, v) in dep, s in rhs_syms(v)
        (s in known || s in (:Keq, :E_total, :L)) || (occursin(r"^[A-Za-z]", string(s)) && continue)
    end
    ok
end

@testset "allosteric dependent-param graph is sound" begin
    for T in ALLOSTERIC_UNDEFVAR_REPRODUCERS
        @test _dep_graph_is_sound(typeof(T()) == T ? T : typeof(T()))
    end
    for spec in MECHANISM_TEST_SPECS
        spec.mechanism isa EnzymeRates.AllostericEnzymeMechanism || continue
        @test _dep_graph_is_sound(typeof(spec.mechanism))
    end
end
```

Include `test/reference/allosteric_undefvar_reproducers.jl` near the top of `test_rate_eq_derivation.jl` (after `using`).

- [ ] **Step 2: Run it — expect FAIL on the reproducers.**
  `julia --project=. -e 'using Test, EnzymeRates, LinearAlgebra, Random, OrdinaryDiffEqFIRK; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_rate_eq_derivation.jl")' 2>&1 | grep -A3 "graph is sound"`
  Expected: the three reproducer `@test`s FAIL (cycle/undefined detected); the promoted specs may pass.
- [ ] **Step 3: Also add a "rate_equation actually runs" regression** for the reproducers (the user-visible symptom):

```julia
@testset "allosteric reproducers: rate_equation is callable" begin
    for T in ALLOSTERIC_UNDEFVAR_REPRODUCERS
        m = T(); em = typeof(m)
        p = Dict(k => 1.0 + 0.01 for k in EnzymeRates.fitted_params(em))
        p[:Keq] = 20000.0; p[:E_total] = 1.0
        concs = (NADH=1.0, Pyruvate=1.0, Lactate=1.0, NAD=1.0)
        @test_nowarn EnzymeRates.rate_equation(m, concs, p)  # currently throws UndefVarError
    end
end
```
  Run: expect FAIL (`UndefVarError`) on HEAD. (Adjust the `concs`/param plumbing to match the real `rate_equation` call signature — read `test_rate_eq_derivation.jl:387-404` for the exact call shape.)
- [ ] **Step 4: Commit the failing tests.** `git add -u && git commit -m "test: failing regressions for the three allosteric UndefVar defects"`

---

## Phase 1 — Factor the kernel (behavior-preserving)

Split the non-allosteric kernel into an assemble step and a solve step, with **zero** behavior change. The full suite must stay green with no re-baselining. This de-risks Phase 2 by giving the allosteric path a clean solver to call.

### Task 1.1: Extract `_solve_dependent_set` from `_dependent_param_exprs_kernel`

**Files:**
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl:315-433` (`_dependent_param_exprs_kernel`).

**Interfaces:**
- Produces: `_solve_dependent_set(A::Matrix{Rational{BigInt}}, rhs::Vector{Rational{BigInt}}, columns::Vector{Symbol}, priority::Vector{Int}) → (dep::Dict{Symbol,Union{Symbol,Expr}}, indep::Tuple)` — the priority-pivoting Gaussian elimination + dependent-expression build (current lines 388-433).
- Produces: `_assemble_constraints(mech, rename; step_params, all_params) → (A, rhs, columns, priority)` — the matrix build + priority scoring (current lines 324-386).

- [ ] **Step 1:** Read `_dependent_param_exprs_kernel` (thermodynamic_constr:315-433) in full. Identify the clean seam: everything through building `A`, `rhs`, `sym_col`/`all_params` (columns), and `priority` is *assembly* (lines ~321-386); everything from `pivot_entries` onward is *solve* (lines ~388-433).
- [ ] **Step 2:** Introduce `_solve_dependent_set(A, rhs, columns, priority)` containing the current lines 388-433 verbatim (renaming `all_params`→`columns`, `n_vars=length(columns)`). Introduce `_assemble_constraints(...)` returning `(A, rhs, columns, priority)`. Rewrite `_dependent_param_exprs_kernel` as: `A,rhs,cols,pri = _assemble_constraints(...); return _solve_dependent_set(A,rhs,cols,pri)`. Keep the early `nc == 0` return path.
- [ ] **Step 3: Run the derivation + allosteric-collapse + split-resolution files.** Expect ALL PASS (behavior-preserving):
  ```
  julia --project=. -e 'using Test, EnzymeRates, LinearAlgebra, Random, OrdinaryDiffEqFIRK; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_rate_eq_derivation.jl")'
  julia --project=. test/test_split_resolution.jl
  julia --project=. test/test_allosteric_collapse.jl
  ```
  Expected: green (the Phase-0 reproducer tests still fail — untouched). If any previously-green test changed, the factoring altered behavior — revert and re-seam.
- [ ] **Step 4: Commit.** `git add -u && git commit -m "refactor: split kernel into _assemble_constraints + _solve_dependent_set (no behavior change)"`

---

## Phase 2 — The combined allosteric solve

Rewrite `_dependent_param_exprs(::AllostericEnzymeMechanism)` to assemble one matrix over state-tagged columns and call `_solve_dependent_set` once. This is the correctness-critical task; the QSSA/equilibrium/analytical oracles are the gate.

### Task 2.1: Assemble stacked A/I constraints and solve once

**Files:**
- Modify: `src/rate_eq_derivation.jl:1555-1610` (`_dependent_param_exprs(::Type{AllostericEnzymeMechanism})`).
- Read for reuse: `_state_mechanism`, `_state_step_params`, `_state_all_params`, `_state_wegscheider_rename_map` (used by `_state_raw_dependent_exprs`, :1236-1245); `_thermodynamic_constraints`; the `name(::Parameter, am)` / state-tagged column naming; how `L` and regulator `Kreg`s enter `indep` (current :1590-1604).

**Interfaces:**
- Consumes: `_assemble_constraints`, `_solve_dependent_set` (Task 1.1).
- Produces (unchanged signature): `_dependent_param_exprs(::Type{AllostericEnzymeMechanism}) → (dep, indep)` consumed by `fitted_params`/`parameters` and `_build_dep_assignments`.

- [ ] **Step 1:** For each state `∈ (:A, :I)`, build `(A_s, rhs_s, cols_s, pri_s)` via `_assemble_constraints` on `_state_mechanism(am, state)` with that state's `step_params`/`all_params`/`rename`. Columns are already state-tagged through `name(p, am)` (bare for `:EqualAI`, `K_A_`/`K_I_` for split states). This reuses the exact per-state assembly the current code already trusts.
- [ ] **Step 2:** Merge the two column sets into one ordered `columns` (shared `:EqualAI` symbols appear once — they are literally the same Symbol from both states). Build a combined matrix by placing each state's rows over the shared column index; add regulator `Kreg` and `L` columns (independent, no rows). Concatenate priorities; apply the allosteric rule: non-`:EqualAI` state-specific params (`K_A`/`K_I`, `k_A`/`k_I`) get priority to STAY independent (i.e. the tied/cross-state quantities pivot). Read the current kernel's `priority` scoring (thermodynamic_constr:368-386) and mirror it per state; document the chosen scores.
- [ ] **Step 3:** Call `_solve_dependent_set(A, rhs, columns, priority)` once → `(dep, indep)`. Append `L` and the reg `Kreg`s to `indep` as today (:1590-1604) if not already columns. Return `(dep, indep)`. Delete the merge loop, the `S_I` handling, and the line-1610 filter from this function.
- [ ] **Step 4: Gate on the invariant test (Phase 0).** Run the "graph is sound" testset — the three reproducers must now PASS:
  `julia --project=. -e '…include test_rate_eq_derivation.jl…' 2>&1 | grep -A6 "graph is sound"`
  If a reproducer still fails, the assembly is wrong (missing a column identity or a cross-state row) — do not proceed.
- [ ] **Step 5: Commit** (body builder still uses old sources; some tests will be red until Task 2.2). `git add -u && git commit -m "wip: combined allosteric constraint solve for fitted_params partition"`

### Task 2.2: Point `_build_dep_assignments` at the combined solve

**Files:**
- Modify: `src/rate_eq_derivation.jl:1674-1720` (`_build_dep_assignments`).

**Interfaces:**
- Consumes: `_dependent_param_exprs(M_type)` (Task 2.1) as the single source of dependent expressions.
- Produces (unchanged signature): `_build_dep_assignments(M_type) → (a_assignments::Vector{Expr}, i_assignments::Vector{Expr})`.

- [ ] **Step 1:** Replace the three fragmented sources (`_state_dependent_exprs` ×2, `_collapse_mirror_exprs`, `_i_state_referenced_syms`) with the combined `dep` from `_dependent_param_exprs`. Topologically order `dep` (a dep's RHS references only earlier deps or independents — the solve already guarantees a valid order; sort by that order). Split into A-block and I-block by which body polynomial references each symbol, computed from the actual num/den Exprs (`_allosteric_num_den_exprs`, :1728) — this replaces the `_i_state_referenced_syms` gate. Emit `:EqualAI` reg mirrors as today (:1696-1703) unless the solve already covers them.
- [ ] **Step 2: Gate on ground truth.** Run the full derivation file:
  `julia --project=. -e '…include test_rate_eq_derivation.jl…'`
  The QSSA oracle (`test_reference_qssa`), equilibrium (`test_haldane_equilibrium`), analytical-rate, and performance tests for every allosteric spec must be GREEN. These are not re-baselined — green here is the proof the rewrite is correct. The golden test will be RED (expected — re-baseline in Phase 3).
- [ ] **Step 3:** If any oracle/equilibrium test is red, the equation is genuinely wrong — debug the assembly/priority (Task 2.1), do not proceed. If green, commit. `git add -u && git commit -m "refactor: body dep-assignments read the combined solve; drop S_I gate"`

### Task 2.3: Delete the subsumed machinery

**Files:**
- Modify/delete in `src/rate_eq_derivation.jl`: `_split_resolution`, `struct SplitResolution`, `_collapse_mirror_exprs`, `_i_state_referenced_syms` (and any now-unused per-state merge helpers). Confirm no remaining callers with grep.
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl` — remove `_split_resolution`-only helpers (e.g. `_integer_nullspace`, `_rref_partition`) IF nothing else uses them (grep first).

**Interfaces:** none produced; pure deletion.

- [ ] **Step 1:** `grep -rn "_split_resolution\|_collapse_mirror_exprs\|_i_state_referenced_syms\|SplitResolution\|_rref_partition\|_integer_nullspace" src/` — confirm the only remaining references are the definitions themselves.
- [ ] **Step 2:** Delete the dead functions/structs and their docstrings. Keep documentation on the *surviving* combined-solve path rich (add a docstring explaining the `I-row − A-row` cross-state emergence).
- [ ] **Step 3: Run the full suite.** `julia --project=. -e 'using Pkg; Pkg.test()'`. Everything green EXCEPT the allosteric golden (Phase 3) and the deleted-internal tests (`test_split_resolution.jl`, mirror-text asserts) — those are handled next.
- [ ] **Step 4: Commit.** `git add -u && git commit -m "refactor: delete split-resolution, collapse mirrors, S_I gate (subsumed by combined solve)"`

---

## Phase 3 — Re-baseline and cleanup

Update ONLY the string-golden and deleted-internal tests, after the ground-truth gates are green.

### Task 3.1: Remove/adapt the deleted-internal tests

**Files:**
- Delete: `test/test_split_resolution.jl` (unit-tests the deleted `_split_resolution`). Remove its `include` from `test/runtests.jl:9-21`.
- Modify: `test/test_allosteric_collapse.jl` — the thermodynamic assertions (`abs(veq) < 1e-8`) and `fitted_params`-membership checks stay; update the mirror-text substring asserts (e.g. `occursin("K_I_S_E=K_A_S_E", …)`) to the new dependent-assignment rendering. Keep every `veq` equilibrium assertion unchanged.

- [ ] **Step 1:** Delete `test/test_split_resolution.jl` and its `runtests.jl` include line. (The collapse behavior remains covered by `test_allosteric_collapse.jl`'s equilibrium + membership asserts — note this coverage argument in the commit message.)
- [ ] **Step 2:** Run `test/test_allosteric_collapse.jl`; for each failing mirror-text assert, confirm the *equilibrium* assert in the same block is green, then update the text to match the new rendering. Never weaken a `veq` tolerance.
- [ ] **Step 3: Commit.** `git add -u && git commit -m "test: drop split-resolution unit test; update collapse mirror-text to combined-solve rendering"`

### Task 3.2: Re-baseline the allosteric golden + ligand naming

**Files:**
- Modify: `test/reference/allosteric_golden_reference.txt` (regenerate).
- Possibly modify: the `name(::Parameter…)` renderer(s) in `src/types.jl` if switching allosteric params to ligand-based names requires it (only if the spec's naming change is in scope this PR; if the current names are already ligand-based, no change).

- [ ] **Step 1:** Decide naming: confirm whether current allosteric names are already ligand-based (`K_A_NADH_E`) or positional (`K1`). If already ligand-based, no renderer change is needed; if positional, that is a separate follow-up — flag and defer rather than expand scope.
- [ ] **Step 2:** Regenerate the golden (no built-in script): 
  `julia --project=. -e 'using Test, EnzymeRates, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_allosteric_golden.jl"); write("test/reference/allosteric_golden_reference.txt", join(_allosteric_golden_lines(), "\n") * "\n")'`
  Then **manually inspect the diff** — every changed equation must correspond to a mechanism whose oracle/equilibrium tests are green (Task 2.2). A changed golden with a red oracle is a bug, not a re-baseline.
- [ ] **Step 3:** Run `test/test_allosteric_golden.jl` → green.
- [ ] **Step 4: Commit.** `git add -u && git commit -m "test: re-baseline allosteric golden to combined-solve equations"`

### Task 3.3: Full-suite green + net-diff review

- [ ] **Step 1:** `julia --project=. -e 'using Pkg; Pkg.test()'` → all green.
- [ ] **Step 2:** Verify the performance contract explicitly (allocs==0, t<120e-9) is still green in the derivation file output.
- [ ] **Step 3:** `git diff main --stat` — confirm the change is strongly net-negative in `src/` logic lines (the success criterion). If not, investigate whether patchwork survived.
- [ ] **Step 4: Commit** any final cleanup. Leave the branch ready for PR (do NOT open the PR without Denis).

---

## Self-review notes

- **Spec coverage:** combined solve (2.1), both consumers share it (2.1+2.2), deletions (2.3), naming/re-baseline (3.2), tests incl. equilibrium + acyclicity + promoted failures (0.1/0.2/2.x), fractional-`1//2` is out of scope here (the solve renders it; accept/reject is Denis's call — noted, not implemented).
- **Ground-truth vs re-baseline** separation is enforced in Global Constraints and gated per task.
- **Stop conditions:** if Task 2.1/2.2 cannot make the QSSA/equilibrium oracles green, STOP after Phase 1 (committed, safe) and write a status doc — do not push a plausibly-wrong derivation.
