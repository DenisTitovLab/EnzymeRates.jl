# Allosteric combined solve + unified pivot priority — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the fragmented per-state allosteric dependent-parameter derivation with one combined constraint solve, and remove the shared `-1` pivot sentinel that silently drops binding-only constraint rows.

**Architecture:** Encode the Gaussian-elimination pivot priority as a lexicographic tuple `(is_I_state, type_priority)` with a "none-found" marker `(false, typemin(Int))` strictly below every real priority — so no column is ever un-pivotable, and I-state columns outrank A-state columns (a pinned cross-state split collapses onto the free A-side). The allosteric path stacks both states' constraint blocks over one tagged-column space and calls the shared solver once, deleting `_split_resolution`/`_collapse_mirror_exprs`/`S_I`/merge.

**Tech Stack:** Julia, `Rational{BigInt}` exact linear algebra, `@generated` rate-equation bodies. Design: `docs/superpowers/specs/2026-07-08-allosteric-combined-solve-and-unified-pivot-design.md`. Spike patch (validated ordering): `docs/superpowers/plans/2026-07-08-combined-solve-spike.patch`. Non-allo root-cause evidence: `docs/superpowers/specs/2026-07-08-wegscheider-pivot-nonallosteric.md`.

## Global Constraints

- **Ground-truth oracles are NEVER re-baselined.** Detailed balance `v = 0` at `Q = Keq` (tol `1e-10`), QSSA reference (`test_reference_qssa`), ODE steady state (`test_ode_steadystate`), analytical rate. A failure in any of these means the fix is wrong — STOP, do not rebaseline.
- **Performance contract (hard gate):** `rate_equation` allocation-free and sub-120 ns for every mechanism in `MECHANISM_TEST_SPECS` (`test_rate_equation_performance`). If a change regresses this, STOP and discuss with Denis.
- **No forbidden priority values.** After this change, no parameter may carry a priority that makes it un-pivotable; the only never-pivoted value is the solver's `(false, typemin(Int))` init marker, which is not a real parameter priority.
- **Naming:** allosteric parameters use real structural names (`K_A_ATP_E`); non-allosteric keep positional textbook names (`K1`). Naming is decided upstream of the solver via `name(p, m)`.
- **Rational split coefficients are allowed** (e.g. `K_I = K_A·(…)^(1//2)`).
- **Preserve the exact relative order of `_step_priority`** — the tuple wraps its Int output; `_step_priority` itself is unchanged (it is also the `argmin` naming representative in `_group_rep`).
- **Test env:** the derivation load is ~3 min; `Bash` default 2-min timeout kills it. Run Julia via `run_in_background`. `kill` (SIGTERM) is ignored by Julia — use `kill -9` + verify with `ps`.
- **Commit frequently**; do not use `git add -A` without a preceding `git status`.

---

### Task 1: Characterization tests (RED on HEAD)

Pin the CORRECT behavior before touching code. These are RED on HEAD (`c59f1d7`) and GREEN after the fix.

**Files:**
- Modify: `test/test_allosteric_collapse.jl` (add dead-I identifiability testset)
- Create: `test/test_pivot_priority_regression.jl` (non-allo Wegscheider tie + enumeration monotonicity)
- Modify: `test/runtests.jl` (include the new file)

**Interfaces:**
- Consumes: `EnzymeRates.fitted_params`, `rate_equation`, `compile_mechanism`, `metabolites`, `init_mechanisms`, `_expand_re_to_ss`, `_expand_split_kinetic_group`; the `uni`/`evalrate` helpers already in `test_allosteric_collapse.jl`.
- Produces: three regression testsets the later tasks turn green.

- [ ] **Step 1: Add the dead-I identifiability testset** to `test/test_allosteric_collapse.jl`, inside the `@testset "strict :EqualAI collapse"` block (it reuses the module's `uni`/`evalrate` helpers):

```julia
    @testset "dead-I NonequalAI binding -> K_I identifiable, NOT collapsed" begin
        # I state cannot turn over (OnlyA catalysis) but binds S with its own
        # affinity: K_A_S_E and K_I_S_E are BOTH identifiable (a dead-end E_I·S is
        # in no cycle, so nothing pins K_I to K_A). HEAD over-collapses this.
        am = uni([:NonequalAI, :OnlyA, :EqualAI])
        fp, v, veq = evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
        @test :K_I_S_E in fp                          # NOT collapsed
        v1 = evalrate(am; split=(:K_I_S_E, 1.3))[2]
        v2 = evalrate(am; split=(:K_I_S_E, 5.0))[2]
        @test !isapprox(v1, v2)                       # identifiable (moves the rate)
        s = replace(ER.rate_equation_string(am), " " => "")
        @test !occursin("K_I_S_E=K_A_S_E", s)         # no collapse mirror
    end
```

- [ ] **Step 2: Create `test/test_pivot_priority_regression.jl`** — the non-allo over-count + enumeration monotonicity regressions (adapted from `2026-07-08-wegscheider-pivot-nonallosteric.md`):

```julia
# ABOUTME: Regressions for the shared -1 pivot-priority sentinel — a Wegscheider row
# ABOUTME: whose only pivot is a free-enzyme binding K must still be pivoted, not dropped.
module PivotPriorityRegressionTests
using Test, EnzymeRates
const ER = EnzymeRates

@testset "pivot priority: no split move reduces the fitted-param count" begin
    rxn = ER.@enzyme_reaction begin
        substrates: NADH[C21H29N7O14P2], Pyruvate[C3H4O3]
        products: Lactate[C3H6O3], NAD[C21H27N7O14P2]
        oligomeric_state: 4
    end
    np(m) = length(ER.fitted_params(m))
    base = unique!(collect(ER.init_mechanisms(rxn)))
    seen = Set(base); pool = collect(base); frontier = copy(base)
    for _ in 1:2
        kids = ER.Mechanism[]
        for m in frontier
            append!(kids, ER._expand_re_to_ss(m))
            append!(kids, ER._expand_split_kinetic_group(m))
        end
        newk = ER.Mechanism[]
        for c in kids
            c in seen && continue
            push!(seen, c); push!(pool, c); push!(newk, c)
        end
        frontier = newk
    end
    # On a correct kernel, a split never reduces the identifiable param count.
    offenders = Tuple{Int,Int}[]
    for P in pool, c in ER._expand_split_kinetic_group(P)
        np(c) < np(P) && push!(offenders, (np(P), np(c)))
    end
    @test isempty(offenders)
end
end # module
```

- [ ] **Step 3: Register the new file** in `test/runtests.jl` (match the existing `include(...)` style; place near the other derivation tests).

- [ ] **Step 4: Run the new tests, confirm RED on HEAD.**

Run (background): `julia --project=. -e 'using Pkg; Pkg.test()'` — or scope to the two files.
Expected: the dead-I testset FAILS (`:K_I_S_E in fp` is false — HEAD collapses it); the monotonicity testset FAILS (`offenders` non-empty, dozens of `−1` edges). Capture the offender list and the dead-I `fitted_params` for the commit message.

- [ ] **Step 5: Commit.**

```bash
git add test/test_allosteric_collapse.jl test/test_pivot_priority_regression.jl test/runtests.jl
git commit -m "test: characterize dead-I identifiability + pivot-sentinel over-count (RED)"
```

---

### Task 2: Unified sentinel-free pivot priority (shared kernel)

**Files:**
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl` — `_assemble_constraints` (~314-371), `_solve_dependent_set` (~383-437)

**Interfaces:**
- Consumes: `_step_priority(step, free_enz_set)::Int` (UNCHANGED — still the `argmin` naming rep).
- Produces: `_assemble_constraints(mech, rename; step_params, all_params, is_i_state=false)` now returns `priority::Vector{Tuple{Bool,Int}}`; `_solve_dependent_set(A, rhs, columns, priority::AbstractVector{Tuple{Bool,Int}})` unchanged return `(dep_exprs::Dict{Symbol,Union{Symbol,Expr}}, indep::Tuple)`.

- [ ] **Step 1: Grep for every `_solve_dependent_set` / priority producer** so no caller is missed:

Run: `grep -rn "_solve_dependent_set\|_assemble_constraints\|priority" src/`
Expected: producers are `_assemble_constraints` and its two callers (`_dependent_param_exprs_kernel`, and the new combined solve in Task 3). Confirm no third path builds an `Int` priority vector and calls the solver directly.

- [ ] **Step 2: Change the priority to a lexicographic tuple in `_assemble_constraints`.** Add the `is_i_state` kwarg and build `Tuple{Bool,Int}` priorities:

Replace the signature default block and the `priority = zeros(Int, n_vars)` construction:

```julia
function _assemble_constraints(
    mech::Mechanism,
    rename::AbstractDict{Symbol, Symbol};
    step_params = _step_parameters(mech),
    all_params = _raw_param_symbols(mech),
    is_i_state::Bool = false,
)
```

and

```julia
    # Pivot priority: (is_I_state, type_priority). Lexicographic — an I-state column
    # outranks any A-state / non-allosteric column, so a cross-state affinity split
    # collapses onto the free A-side; within a state the `_step_priority` order holds.
    # No value is a never-pivot sentinel (that lives only in `_solve_dependent_set`).
    priority = fill((is_i_state, 0), n_vars)
    for j in 1:nsteps
        step = step_params[j][1].step
        base = _step_priority(step, free_enz_set)
        if is_equilibrium(flat[j][1])
            s = step_name(step_params[j][1])
            haskey(sym_col, s) && (priority[sym_col[s]] = (is_i_state, base))
        else
            for (offset, p) in enumerate(step_params[j])
                s = step_name(p)
                haskey(sym_col, s) &&
                    (priority[sym_col[s]] = (is_i_state, base + offset - 1))
            end
        end
    end
    return A, rhs, columns, priority
end
```

- [ ] **Step 3: De-sentinelize `_solve_dependent_set`.** Change the priority type and the `best_pri` init so a `-1` type-priority column is a valid last-resort pivot:

Signature: `priority::AbstractVector{Tuple{Bool,Int}}`. Init and comparison:

```julia
        best_col, best_pri = 0, (false, typemin(Int))
        for c in 1:n_vars
            c in pivot_col_set && continue
            wA[i, c] == 0 && continue
            if priority[c] > best_pri
                best_pri = priority[c]
                best_col = c
            end
        end
```

(The `best_col == 0` branch still drops only genuinely redundant `rhs == 0` rows — now reached only when *every* column of the row is already pivoted, i.e. a true `0 = 0`, never merely because the survivors are binding K's.)

- [ ] **Step 4: Run the non-allo regressions, confirm GREEN.**

Run (background): scope to `test/test_pivot_priority_regression.jl`.
Expected: monotonicity testset PASSES (`offenders` empty). Do NOT expect the full suite green yet — many non-allo `fitted_params`/golden baselines now shift (Task 5).

- [ ] **Step 5: Spot-check the LDH over-count fix** with a scratch script: build the 6-group non-allo LDH parent reached in enumeration (or find it in the Task-1 pool) and confirm `fitted_params` length is 6 (was 7) and the rate string shows `K_Lactate_E = K_Lactate_ENADH`.

- [ ] **Step 6: Commit.**

```bash
git add src/thermodynamic_constr_for_rate_eq_derivation.jl
git commit -m "fix: lexicographic-tuple pivot priority, remove -1 never-pivot sentinel"
```

---

### Task 3: Allosteric combined constraint solve

**Files:**
- Modify: `src/rate_eq_derivation.jl` — `_dependent_param_exprs(::AllostericEnzymeMechanism)` (~1544-1641), `_build_dep_assignments` (~1688-1754), `rate_equation_string(::AllostericEnzymeMechanism, ::ReducedMode)` (~1861-1907)
- Reference (do NOT re-derive): `docs/superpowers/plans/2026-07-08-combined-solve-spike.patch`

**Interfaces:**
- Consumes: `_assemble_constraints(...; is_i_state)` (Task 2), `_solve_dependent_set` (Task 2), `_state_mechanism`, `_state_step_params`, `_state_wegscheider_rename_map`, `_state_all_params`, `_expr_references_any`, `name(::Kreg, am)`, `_state_rate_polys`.
- Produces: `_combined_state_dependent_exprs(am)::(Dict,Tuple)`, `_i_state_symbol_set(am)::Set{Symbol}`; `_dependent_param_exprs(::Allo)` and `_build_dep_assignments` route through the single combined solve.

- [ ] **Step 1: Add `_combined_state_dependent_exprs`** just before `_dependent_param_exprs(::Allo)`. The tuple priority (Task 2) carries the I-above-A ordering, so this only stacks — no manual offset:

```julia
function _combined_state_dependent_exprs(am::AllostericMechanism)
    function state_system(state)
        cm = _state_mechanism(am, state)
        sp = _state_step_params(am, state)
        _assemble_constraints(cm, _state_wegscheider_rename_map(am, state);
                              step_params = sp, all_params = _state_all_params(cm, sp),
                              is_i_state = (state === :I))
    end
    A_A, rhs_A, cols_A, pri_A = state_system(:A)
    A_I, rhs_I, cols_I, pri_I = state_system(:I)

    # Union columns: A-state first, then any I-only column (a shared `:EqualAI` group
    # carries the same bare Symbol in both states and coincides — it keeps its A tag).
    columns = copy(cols_A)
    col_index = Dict(c => i for (i, c) in enumerate(columns))
    priority = copy(pri_A)
    for (j, c) in enumerate(cols_I)
        haskey(col_index, c) && continue
        push!(columns, c)
        col_index[c] = length(columns)
        push!(priority, pri_I[j])
    end

    # Stack the two per-state constraint blocks over the combined column space and
    # solve once. Cross-state ties emerge as `I-row − A-row` (the `log Keq` cancels).
    A = zeros(Rational{BigInt}, size(A_A, 1) + size(A_I, 1), length(columns))
    for i in axes(A_A, 1), (j, c) in enumerate(cols_A)
        A_A[i, j] == 0 || (A[i, col_index[c]] = A_A[i, j])
    end
    off = size(A_A, 1)
    for i in axes(A_I, 1), (j, c) in enumerate(cols_I)
        A_I[i, j] == 0 || (A[off + i, col_index[c]] = A_I[i, j])
    end
    rhs = vcat(rhs_A, rhs_I)
    return _solve_dependent_set(A, rhs, columns, priority)
end
```

- [ ] **Step 2: Replace the body of `_dependent_param_exprs(::Allo)`** with the combined solve plus regulator/`L` handling (NO mirror loop). Replace everything from `am = AllostericMechanism(...)` through the final `return`:

```julia
function _dependent_param_exprs(
    ::Type{AllostericEnzymeMechanism{CM,CS,RS}},
) where {CM,CS,RS}
    am = AllostericMechanism(AllostericEnzymeMechanism{CM,CS,RS}())
    dep, indep = _combined_state_dependent_exprs(am)

    # Regulator-site affinities complete no catalytic thermodynamic cycle, so they
    # are independent — except an `:EqualAI` regulator, whose I-name mirrors its
    # shared A-name. `L` (the conformational constant) is always independent.
    reg_params_a = Symbol[]
    reg_params_i = Symbol[]
    for site in regulatory_sites(am)
        for (lig, tag) in zip(ligands(site), allo_states(site))
            tag === :OnlyI || push!(reg_params_a, name(Kreg(site, lig, :A), am))
            if tag === :EqualAI
                dep[name(Kreg(site, lig, :I), am)] = name(Kreg(site, lig, :A), am)
            elseif tag === :NonequalAI || tag === :OnlyI
                push!(reg_params_i, name(Kreg(site, lig, :I), am))
            end
        end
    end
    return dep, Tuple(p for p in (indep..., reg_params_a..., reg_params_i..., :L)
                      if p ∉ keys(dep))
end
```

- [ ] **Step 3: Replace `_build_dep_assignments` and add `_i_state_symbol_set`** (split the single combined `dep` into A-block / I-block for emission):

```julia
function _i_state_symbol_set(am::AllostericMechanism)
    cols_A = _state_all_params(_state_mechanism(am, :A), _state_step_params(am, :A))
    cols_I = _state_all_params(_state_mechanism(am, :I), _state_step_params(am, :I))
    syms = Set{Symbol}(setdiff(cols_I, cols_A))
    for site in regulatory_sites(am)
        for (lig, tag) in zip(ligands(site), allo_states(site))
            tag === :OnlyA && continue
            push!(syms, name(Kreg(site, lig, :I), am))
        end
    end
    syms
end

function _build_dep_assignments(M_type::Type{<:AllostericEnzymeMechanism})
    am = AllostericMechanism(M_type())
    dep, _ = _dependent_param_exprs(M_type)
    # A-block first so an I-block `:EqualAI` regulator mirror (`K_I_reg = K_A_reg`)
    # finds its A-name defined. The combined solve expresses every dependent purely
    # in independent columns, so no dependent reads another — order within a block
    # is free.
    i_syms = _i_state_symbol_set(am)
    a_assignments = Expr[]
    i_assignments = Expr[]
    for (sym, rhs) in sort(collect(dep); by = first)
        push!(sym in i_syms ? i_assignments : a_assignments, Expr(:(=), sym, rhs))
    end
    return a_assignments, i_assignments
end
```

- [ ] **Step 4: Repoint `rate_equation_string(::Allo)`** at the combined assignments. Replace the block from `am = AllostericMechanism(m)` through the inactive-state loop (delete the native `dep_A` + `_partition_constraint_lines!` rendering; keep the `sort!`/`v_line`/assembly below):

```julia
    _, indep = _dependent_param_exprs(M)
    hw_params = (indep..., :Keq, :E_total)
    mets = metabolites(m)

    # Every dependent assignment comes from the single combined solve — the same set
    # the compiled body assigns — split into Wegscheider/Haldane by Keq-reference.
    keq_set = Set([:Keq])
    a_assignments, i_assignments = _build_dep_assignments(M)
    weg_lines, hal_lines = String[], String[]
    for a in (a_assignments..., i_assignments...)
        sym = a.args[1]
        expr = a.args[2]
        is_haldane = _expr_references_any(expr, keq_set)
        line = "$sym = $(_expr_to_string(expr))"
        push!(is_haldane ? hal_lines : weg_lines, line)
    end
```

- [ ] **Step 5: Verify the derivation is correct** with the spike scratch scripts (they still apply): run the 6 uni collapse cases + the 3 reproducers and confirm the spike results reproduce (case1 collapses `K_I_S_E=K_A_S_E`; dead-I keeps `K_I_S_E` free+identifiable; SS `koff_I` collapses/`kon_I` free; all detbal `|v| ≲ 1e-16`; reproducer param counts 12→10, 13→12). Copy `scratchpad/spikeB.jl`, `spikeC.jl` from the session scratchpad or re-derive from Task-1 helpers.

Expected: identical to the spike (the tuple encoding reproduces the validated ordering).

- [ ] **Step 6: Flip the two `@test_broken`** in `test/test_rate_eq_derivation.jl:1301-1305` (reproducers 2/3 detailed balance) from `@test_broken maxv < 1e-8` to `@test maxv < 1e-8` (the single `if i == 1 ... else ... end` now uses `@test` for all three):

```julia
        @test maxv < 1e-8
```

- [ ] **Step 7: Run the dead-I test + reproducer detbal + oracles, confirm GREEN.** Do NOT expect allosteric goldens/eq_hash green yet (Task 5).

- [ ] **Step 8: Commit.**

```bash
git add src/rate_eq_derivation.jl test/test_rate_eq_derivation.jl
git commit -m "feat: allosteric derivation via one combined constraint solve"
```

---

### Task 4: Delete the subsumed machinery

**Files:**
- Modify: `src/rate_eq_derivation.jl` (delete functions), possibly `src/types.jl` (`SplitResolution` struct)
- Delete: `test/test_split_resolution.jl`; Modify: `test/runtests.jl`

**Interfaces:** none produced; this removes now-dead code.

- [ ] **Step 1: For each candidate, grep for remaining usages** before deleting:

Run: `grep -rn "_split_resolution\|SplitResolution\|_collapse_mirror_exprs\|_i_state_referenced_syms\|_state_dependent_exprs\|_flat_expr_syms\|_partition_constraint_lines" src/ test/`
Expected: after Task 3, each appears only at its own definition (and `test_split_resolution.jl`). Any surviving `src/` usage means the rewrite missed a path — investigate, do not force-delete.

- [ ] **Step 2: Delete the dead functions** in `src/rate_eq_derivation.jl`: `_split_resolution`, `_collapse_mirror_exprs`, `_i_state_referenced_syms`, `_state_dependent_exprs`, `_flat_expr_syms`, `_partition_constraint_lines!`. Delete the `SplitResolution` struct (wherever defined). Keep everything the combined solve still calls (`_state_mechanism`/`_state_step_params`/`_state_all_params`/`_state_rate_polys`/`_state_wegscheider_rename_map`).

- [ ] **Step 3: Delete `test/test_split_resolution.jl`** and remove its `include` from `test/runtests.jl`.

- [ ] **Step 4: Run the oracle subset** (detailed balance, QSSA, dead-I, reproducers) — confirm still GREEN (deletion changed no behavior).

- [ ] **Step 5: Commit.**

```bash
git status
git add -u src/ test/ && git rm test/test_split_resolution.jl
git commit -m "refactor: delete _split_resolution/_collapse_mirror_exprs/S_I/native per-state solve"
```

---

### Task 5: Single rebaseline pass (REVIEWED CHECKPOINT)

This is the delicate work. The pivot fix and combined solve change many hardcoded baselines. Rebaseline them **once**, understanding each change; never rebaseline an oracle.

**Files (expected; confirm by running the suite):**
- `test/reference/allosteric_golden_reference.txt` (regenerate under real names)
- `test/test_types.jl`, `test/test_rate_eq_derivation.jl`, `test/test_mechanism_enumeration.jl`, `test/test_allosteric_golden.jl` (hardcoded `fitted_params`/counts/eq_hashes/enumeration totals)

**Interfaces:** none; behavior is frozen from Tasks 2-4.

- [ ] **Step 1: Run the full suite, capture every failure** (background; write to a file, not a pipe):

Run: `julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/.../fullsuite.out 2>&1`

- [ ] **Step 2: Triage each failure into one of two buckets.**
  - **Oracle failure** (detailed balance, QSSA, ODE, analytical, perf) → STOP. The fix is wrong; return to Task 2/3. Do NOT edit the oracle.
  - **Baseline failure** (a hardcoded `fitted_params` count, golden string, `eq_hash`, enumeration total, mirror-direction assertion) → candidate rebaseline. For each, state in one line WHY it changed (over-count fixed → count −1; split collapses `K_I` onto `K_A`; real names; enumeration pool shrank as inconsistent duplicates collapsed). If a change cannot be explained by the fix, STOP and investigate.

- [ ] **Step 3: Regenerate the allosteric golden reference** and review the diff line-by-line (it encodes the identifiable dimension AND the naming). Use the repo's golden-regeneration path (grep `allosteric_golden_reference` for the writer/`JULIA_REGEN`-style switch). Confirm: real names throughout, forbidden splits show `K_I_… = K_A_…`, honorable splits stay fitted, no `UndefVar`-shaped references.

- [ ] **Step 4: Update each explained baseline** to its corrected value. Group edits by test file; one commit per file with the WHY in the message.

- [ ] **Step 5: Re-run the full suite, confirm GREEN** (except the perf gate, Task 6).

- [ ] **Step 6: Commit** (per file, e.g.):

```bash
git add test/reference/allosteric_golden_reference.txt
git commit -m "test: rebaseline allosteric goldens (real names, true identifiable dim)"
```

---

### Task 6: Performance gate + final verification

**Files:** none (verification only).

- [ ] **Step 1: Run the performance test.**

Run: `test_rate_equation_performance` in `test/test_rate_eq_derivation.jl`.
Expected: `allocs == 0` and `t < 120e-9` for every mechanism. If regressed → STOP and discuss with Denis (do not weaken the bound).

- [ ] **Step 2: Full-suite green confirmation.**

Run (background): `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: all green.

- [ ] **Step 3: End-to-end detailed-balance sweep** over all allosteric `MECHANISM_TEST_SPECS` + the 3 reproducers: `|v| < 1e-8` at `Q = Keq` for 100 random parameter draws each. Confirms the invariant holds beyond the pinned tests.

- [ ] **Step 4: Confirm naming** — grep the regenerated allosteric goldens for any positional `K<digit>` name; expect none (allosteric is fully real-named). Non-allo goldens keep positional names.

- [ ] **Step 5: Final commit / branch ready for review.**

```bash
git status
git commit -am "chore: combined solve + unified pivot — suite green, perf gate held" || true
```

---

## Notes for the executor

- The core code is the spike patch (`2026-07-08-combined-solve-spike.patch`) re-expressed with the tuple priority (Tasks 2-3); it is validated ordering, not a fresh derivation.
- The rebaseline (Task 5) is data-driven: exact new values come from running the suite. The plan gives the triage rule and the review protocol, not pre-guessed numbers — inventing baseline values would be worse than reading them from a correct run.
- If any oracle or the perf gate fails at any point, STOP: that is signal the fix is wrong, not a baseline to update.
