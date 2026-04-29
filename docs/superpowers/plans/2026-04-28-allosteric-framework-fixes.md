# Allosteric Framework Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix 20 issues identified by parallel code reviews of the allosteric refactor (commits `29aec43..aef9345`). Two are critical correctness bugs (`_count_allosteric_rate_monomials` parallel drift; `_dependent_param_exprs` mixed-tag bug). The rest are dead-code cleanup, naming consistency, validation tightening, enumeration coverage gap, test categorization, and documentation updates.

**Architecture:** Framework fixes (Phase 1) make the rate equation thermodynamically consistent for all valid tag combinations and stop locking in buggy snapshot values. Test reform (Phase 2) replaces the lumped `expected_n_haldane` field with three categorized fields and drops the `expected_identifiability_deficit` field that measured a methodologically-flawed monomial-count heuristic. DSL/constructor improvements (Phase 3) move parse-time errors closer to the source. Enumeration (Phase 4) closes a coverage gap so m_all-shape mechanisms are reachable. Cleanups (Phase 5) finish the rename + comment hygiene.

**Tech Stack:** Julia, EnzymeRates.jl, type-parameter-based mechanism encoding.

---

## Files Touched

- `src/rate_eq_derivation.jl` — `_dependent_param_exprs`, `_count_allosteric_rate_monomials`, `_kcat_forward`, `_allosteric_num_den_exprs`, `_T_rename`, `_onlyT_syms` deletion, `Base.show`, `@assert` relaxation. Tasks 1, 2, 3, 8.
- `src/types.jl` — constructor validation, `Base.show`, `cat_allo_state`/`reg_allo_state` accessors. Task 8.
- `src/dsl.jl` — `_parse_steps_block_with_groups` parse-time error. Task 7.
- `src/mechanism_enumeration.jl` — `_expand_add_allosteric_regulator` `:EqualRT` extension; rename `_expand_change_group_tag` → `_expand_change_allo_state` and other `tag` → `allo_state` function names. Tasks 9, 10.
- `test/test_rate_eq_derivation.jl` — `test_constraint_counting` rewrite; new regression tests for Haldane-equilibrium of `:NonequalRT binding + :EqualRT catalysis` mechanism, empty-ligand error, non-consecutive kinetic_group error, `Base.show` round-trip. Remove redundant `:OnlyR sub + :OnlyT prod` test. Tasks 4, 6, 8, 12.
- `test/mechanism_definitions_for_test_enzyme_derivation.jl` — replace `expected_n_haldane` / `expected_n_wegscheider` with three categorized fields; drop `expected_identifiability_deficit`; revert PK to `:EqualRT` catalysis; update analytical formula. Tasks 4, 5, 6.
- `test/test_dsl.jl` — bare-step error message verification; rename `(smoke)` testset. Task 7, 11.
- `.claude/CLAUDE.md` — terminology and source-layout documentation updates. Task 11.

---

### Task 1: `_T_rename` dataflow extension + dep generation + skip dead T-mirrors

The framework has a thermodynamic-consistency bug: when catalysis is `:EqualRT` but a binding step is `:NonequalRT`, the T-state Haldane closure is lost. At chemical equilibrium, the T-state numerator stays non-zero — second-law violation.

The fix has THREE coordinated parts that must land together:

1. **Extend `_T_rename`**: it currently maps only `:NonequalRT` catalytic-group symbols (e.g., `K1 → K1_T`). Extend it to also include synthesized dep-symbol mappings (e.g., `k5r → k5r_T`) for any `:EqualRT`-tagged dep whose RHS expression references a `:NonequalRT` symbol. This is the load-bearing change — every consumer of `rename_T` (`_allosteric_num_den_exprs`, `_kcat_forward`, `_count_allosteric_rate_monomials`, `_build_dep_assignments`) automatically picks up the new mappings without further edits.

2. **Generate the synthesized dep entries** in `_dependent_param_exprs` so the rate equation body has an assignment for `k5r_T = k5f * K6 * K8 / (Keq * K1_T * K3)`.

3. **Emit the assignments** in `_build_dep_assignments` (parallel function used by `_build_allosteric_rate_body` and `rate_equation_string`).

Same function (`_dependent_param_exprs`) also skips dead `:EqualRT` mirror dep entries (`K1_T = K1`, `k5f_T = k5f`, etc.) when `_t_state_dead` is true — they're already elided from the rate equation body and producing them only inflates `length(dep_exprs)` for tests.

**Files:**
- Modify: `src/rate_eq_derivation.jl`:
  - `_T_rename(m)` (around lines 1162-1174) — extend with dataflow second pass
  - `_dependent_param_exprs(::Type{<:AllostericEnzymeMechanism})` (around lines 1219-1259) — generate Case B dep entries + skip dead mirrors
  - `_build_dep_assignments` (around lines 1325-1389) — emit Case B assignments
- Add: `test/test_rate_eq_derivation.jl` — regression test in "Allosteric edge cases"

- [ ] **Step 1: Extend `_T_rename` with the dataflow pass (also drops `:OnlyT` since Task 9 of the prior refactor made it impossible)**

Locate `_T_rename(m::AllostericEnzymeMechanism)` in `src/rate_eq_derivation.jl` (around line 1162). Replace its body with a two-pass version: first pass collects catalytic-group `:NonequalRT` params, second pass adds synthesized dep-symbol mappings for `:EqualRT`-tagged deps whose RHS references a renamed symbol.

The current source still has a `:OnlyT` branch in the first pass. Since the constructor rejects `:OnlyT` catalytic groups (Task 9 of the prior refactor, commit `66a6119`), that branch is dead code. The new code drops it. **This also satisfies what Task 3 Step 7 was going to do** — Task 3's docstring says "Remove `:OnlyT` from `_T_rename`," and Task 1 Step 1 here does it. Drop Task 3 Step 7 from execution (or treat it as already-applied when Task 3 reaches it).

Old:
```julia
function _T_rename(m::AllostericEnzymeMechanism)
    cm = catalytic_mechanism(m)
    rename = Dict{Symbol, Symbol}()
    for g in kinetic_groups(cm)
        tag = cat_allo_state(m, g)
        (tag == :NonequalRT || tag == :OnlyT) || continue
        rep = first(steps_in_group(cm, g))
        for s in _group_param_symbols(cm, rep)
            rename[s] = _rename_params_T(s)
        end
    end
    rename
end
```

New:
```julia
function _T_rename(m::AllostericEnzymeMechanism)
    cm = catalytic_mechanism(m)
    rename = Dict{Symbol, Symbol}()
    # First pass: catalytic-group params with :NonequalRT (independent K_R / K_T).
    # :OnlyT catalytic groups are rejected at construction, so we don't
    # need a branch for them.
    for g in kinetic_groups(cm)
        cat_allo_state(m, g) == :NonequalRT || continue
        rep = first(steps_in_group(cm, g))
        for s in _group_param_symbols(cm, rep)
            rename[s] = _rename_params_T(s)
        end
    end
    # Second pass: derived dep symbols whose RHS references any T-renamed
    # symbol need their own T-state name. After Gaussian elimination of
    # the constraint matrix, dep RHSs reference only independent params
    # (never other deps), so a single non-iterating pass suffices.
    # Without this pass, a Haldane closure with :EqualRT catalysis whose
    # formula references a :NonequalRT binding K would have its T-state
    # value undefined at runtime — N_T uses the R-state value, breaking
    # Haldane consistency at chemical equilibrium.
    dep_R_all, _ = _dependent_param_exprs(typeof(cm))
    renamed_set = Set{Symbol}(keys(rename))
    for (k, v) in dep_R_all
        haskey(rename, k) && continue
        _expr_references_any(v, renamed_set) || continue
        rename[k] = _rename_params_T(k)
    end
    rename
end
```

- [ ] **Step 2: Generate Case B dep entries + skip dead mirrors in `_dependent_param_exprs`**

Locate `_dependent_param_exprs(::Type{AllostericEnzymeMechanism{CM,CS,RS}})` (around line 1219). Find the T-state iteration block:

Old:
```julia
    dep_T = Dict{Symbol, Union{Symbol, Expr}}()
    indep_T_list = Symbol[]
    for (k, v) in dep_R_all
        _expr_references_any(v, r_only_syms) && continue
        t_k = get(rename_T, k, k)
        # `:EqualRT` (k unchanged) is already in dep_R; skip duplicate.
        t_k == k && continue
        dep_T[t_k] = substitute_params_expr(v, T_subs)
    end
    for p in indep_R_all
        p ∈ r_only_syms && continue
        if !haskey(rename_T, p)
            # `:EqualRT` independent: its T-state mirror equals the R-state
            # symbol. Add p_T = p as a dep, do not duplicate as indep.
            dep_T[_rename_params_T(p)] = p
        else
            # `:NonequalRT` and `:OnlyT` get a distinct T-state independent.
            push!(indep_T_list, _rename_params_T(p))
        end
    end
```

New:
```julia
    t_state_dead_flag = _t_state_dead(m)

    dep_T = Dict{Symbol, Union{Symbol, Expr}}()
    indep_T_list = Symbol[]

    # Generate T-state dep entries for every R-state dep that has a
    # T-state version per `rename_T`. After Step 1's extension, this
    # includes both Case A (dep symbol's catalytic group is :NonequalRT —
    # the symbol itself is in rename_T) and Case B (dep symbol is
    # :EqualRT-tagged but its RHS references a :NonequalRT symbol — the
    # extended rename_T still contains a synthesized mapping).
    for (k, v) in dep_R_all
        _expr_references_any(v, r_only_syms) && continue
        t_k = get(rename_T, k, nothing)
        t_k === nothing && continue
        dep_T[t_k] = substitute_params_expr(v, T_subs)
    end

    # When the T-state cycle is dead (any :OnlyR group), skip generating
    # :EqualRT mirror entries (K1_T = K1, k5f_T = k5f, etc.). They're
    # already elided from the rate equation body in
    # _build_allosteric_rate_body, so producing them here only inflates
    # length(dep_exprs).
    for p in indep_R_all
        p ∈ r_only_syms && continue
        if haskey(rename_T, p)
            push!(indep_T_list, _rename_params_T(p))
        elseif !t_state_dead_flag
            dep_T[_rename_params_T(p)] = p
        end
    end
```

- [ ] **Step 3: Apply the corresponding update in `_build_dep_assignments`**

Locate `_build_dep_assignments` in `src/rate_eq_derivation.jl` (around line 1325). Find the T-state iteration block:

Old:
```julia
    for (sym, expr_kd) in sorted_deps
        t_sym = get(rename_T, sym, sym)
        t_sym == sym && continue
        if _expr_references_any(expr_kd, r_only_syms)
            push!(t_assignments, Expr(:(=), t_sym, 0))
        else
            push!(t_assignments, Expr(:(=), t_sym,
                substitute_params_expr(expr_kd, T_subs)))
        end
    end
```

New:
```julia
    # Emit a T-state assignment for every dep that has a T-state name
    # in rename_T. Step 1's extension to _T_rename includes both
    # :NonequalRT-tagged dep symbols (Case A) and synthesized T-names
    # for :EqualRT-tagged derived deps whose RHS references a
    # :NonequalRT symbol (Case B). The unified lookup catches both.
    for (sym, expr_kd) in sorted_deps
        t_sym = get(rename_T, sym, nothing)
        t_sym === nothing && continue
        if _expr_references_any(expr_kd, r_only_syms)
            push!(t_assignments, Expr(:(=), t_sym, 0))
        else
            push!(t_assignments, Expr(:(=), t_sym,
                substitute_params_expr(expr_kd, T_subs)))
        end
    end
```

- [ ] **Step 4: Add a Haldane-equilibrium regression test in "Allosteric edge cases"**

Open `test/test_rate_eq_derivation.jl`. Locate the `@testset "Allosteric edge cases" begin` block. Add this regression test at the end of the testset (just before the closing `end` of the `@testset`):

```julia
    # Regression: :NonequalRT substrate + :EqualRT catalysis must produce
    # zero rate at chemical equilibrium. The framework derives a T-state
    # Haldane (k2r_T) from the :EqualRT k2f because the dep expression
    # for k2r references :NonequalRT K1, so _T_rename's dataflow pass
    # synthesizes a T-name for k2r and substitutes it into N_T.
    cm_mixed = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            [E, S] ⇌ [E_S]
            [E_S] <--> [E_P]
            [E, P] ⇌ [E_P]
        end
    end
    m_mixed = EnzymeRates.AllostericEnzymeMechanism(
        cm_mixed,
        (2, (:NonequalRT, :EqualRT, :EqualRT)),
        (((:I,), 2, (:NonequalRT,)),),
    )
    Keq_val = 5.0
    p_eq = (K1=0.3, k2f=8.0, K3=0.7,
            K1_T=2.5,
            K_I_reg1=1.0, K_I_T_reg1=4.0,
            L=2.0, Keq=Keq_val, E_total=1.0)
    # At chemical equilibrium: P = Keq · S
    S_eq = 1.5
    P_eq = Keq_val * S_eq
    rate_eq = rate_equation(m_mixed, (S=S_eq, P=P_eq, I=0.5), p_eq)
    @test isapprox(rate_eq, 0.0; atol=1e-10)
```

- [ ] **Step 5: Verify the regression test passes — MUST pass before proceeding**

This is the gate. The framework fix is correct iff the new `m_mixed` test asserts `rate_eq ≈ 0.0` at chemical equilibrium.

Run the regression in isolation by extracting just the `m_mixed` block into a one-off script:

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
julia --project -e '
using Test
using EnzymeRates
@testset "m_mixed Haldane gate" begin
    cm_mixed = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            [E, S] ⇌ [E_S]
            [E_S] <--> [E_P]
            [E, P] ⇌ [E_P]
        end
    end
    m_mixed = EnzymeRates.AllostericEnzymeMechanism(
        cm_mixed,
        (2, (:NonequalRT, :EqualRT, :EqualRT)),
        (((:I,), 2, (:NonequalRT,)),),
    )
    Keq_val = 5.0
    p_eq = (K1=0.3, k2f=8.0, K3=0.7,
            K1_T=2.5,
            K_I_reg1=1.0, K_I_T_reg1=4.0,
            L=2.0, Keq=Keq_val, E_total=1.0)
    S_eq = 1.5
    P_eq = Keq_val * S_eq
    rate_eq = rate_equation(m_mixed, (S=S_eq, P=P_eq, I=0.5), p_eq)
    @test isapprox(rate_eq, 0.0; atol=1e-10)
end
'
```

Expected output: `Test Summary: ... Pass 1 ... Total 1`.

**If the test fails, the framework fix is incomplete — STOP and report.** Don't proceed to Step 6 until this test passes. Other test failures across the full suite (mismatch in `expected_n_haldane` etc.) are expected at this stage and are fixed in Tasks 4-5.

- [ ] **Step 6: Run the full test suite (failures from spec-count mismatches are expected)**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: the new regression test passes; many existing constraint-counting tests fail because `expected_n_haldane` values now don't match the new dep counts. Do NOT update spec values here — they're re-measured in Task 5.

- [ ] **Step 7: Commit**

```bash
git add src/rate_eq_derivation.jl test/test_rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
Make _T_rename dataflow-aware + skip dead t_state_dead mirrors

Three coordinated changes restoring Haldane consistency for
:EqualRT-catalysis + :NonequalRT-binding mechanisms:

1. _T_rename gains a dataflow pass: any dep symbol whose RHS
   expression references a T-renamed symbol gets its own T-name
   added to the rename map. Iterates to a fixpoint to handle chained
   deps. This is the load-bearing change — every consumer of
   rename_T (_allosteric_num_den_exprs, _kcat_forward,
   _count_allosteric_rate_monomials, _build_dep_assignments)
   automatically rewrites the dep symbol in N_T / Q_T to its T-state
   counterpart, picking up the corresponding T-state Haldane value.

2. _dependent_param_exprs generates the synthesized dep entries
   (e.g., k5r_T = k5f·K6·K8/(Keq·K1_T·K3)) so the rate equation body
   has an assignment for the new T-symbol.

3. _build_dep_assignments emits the same assignments in
   t_assignments for the rate-body and string consumers.

Same function (_dependent_param_exprs) also skips :EqualRT mirror
entries (K1_T = K1, k5f_T = k5f, ...) when _t_state_dead is true —
they're already elided from the rate equation body, so producing
them here only inflates length(dep_exprs).

Regression test added in "Allosteric edge cases": a mechanism with
:NonequalRT substrate + :EqualRT catalysis gives rate=0 at chemical
equilibrium (Keq·S = P) — verifies the T-state Haldane closure.

Many existing tests still fail because spec values reflect the
old (buggy) mirror counts. Those will be re-measured in Task 5.
EOF
)"
```

---

### Task 2: Patch `_count_allosteric_rate_monomials` parallel-drift

`_count_allosteric_rate_monomials` builds the same allosteric rate-equation polynomial structure as `_allosteric_num_den_exprs`, but in flat `POLY` form (for monomial counting). It missed two changes:
1. Task 3's removal of the `n_reg == CatN || continue` filter from the numerator.
2. Task 5's elision of the `L*num_T` term when `_t_state_dead`.

Patch both. Document the parallel constraint via comments.

**Files:**
- Modify: `src/rate_eq_derivation.jl` (function `_count_allosteric_rate_monomials`, around lines 1570-1643)

- [ ] **Step 1: Drop the `n_reg == CatN` filter in `num_poly_for_conf`**

Locate the inner function `num_poly_for_conf` (around line 1607). Replace:

Old:
```julia
    function num_poly_for_conf(N_cat, Q_cat, reg_Qs, L_factor)
        n_term = poly_mul(N_cat, _poly_power(Q_cat, CatN - 1))
        for i in eachindex(RS)
            n_reg = regulatory_site_multiplicity(m, i)
            n_reg == CatN || continue
            n_term = poly_mul(n_term, _poly_power(reg_Qs[i], n_reg))
        end
        L_factor === nothing ? n_term : poly_mul(poly_sym(L_factor), n_term)
    end
```

New:
```julia
    # Numerator: per-state catalytic flux × Q_cat^(CatN-1) × all reg-site
    # factors at their multiplicity. MUST mirror `make_num_term` in
    # `_allosteric_num_den_exprs` — they build the same polynomial in
    # different representations. Drift between them produces
    # inconsistent monomial counts vs the rate equation.
    function num_poly_for_conf(N_cat, Q_cat, reg_Qs, L_factor)
        n_term = poly_mul(N_cat, _poly_power(Q_cat, CatN - 1))
        for i in eachindex(RS)
            n_reg = regulatory_site_multiplicity(m, i)
            n_term = poly_mul(n_term, _poly_power(reg_Qs[i], n_reg))
        end
        L_factor === nothing ? n_term : poly_mul(poly_sym(L_factor), n_term)
    end
```

- [ ] **Step 2: Skip the `L*num_T` term when `_t_state_dead`**

Locate the assembly block (around lines 1626-1633). Replace:

Old:
```julia
    full_num = poly_add(
        num_poly_for_conf(N_cat_R, Q_cat_R, reg_Q_R, nothing),
        num_poly_for_conf(N_cat_T, Q_cat_T, reg_Q_T, :L),
    )
    full_den = poly_add(
        den_poly_for_conf(Q_cat_R, reg_Q_R, nothing),
        den_poly_for_conf(Q_cat_T, reg_Q_T, :L),
    )
```

New:
```julia
    # Drop the L*N_T numerator branch when t_state_dead: it expands to
    # L * 0 * polynomial, contributing zero to monomial count. MUST
    # mirror `_allosteric_num_den_exprs` post-Task 5.
    if _t_state_dead(m)
        full_num = num_poly_for_conf(N_cat_R, Q_cat_R, reg_Q_R, nothing)
    else
        full_num = poly_add(
            num_poly_for_conf(N_cat_R, Q_cat_R, reg_Q_R, nothing),
            num_poly_for_conf(N_cat_T, Q_cat_T, reg_Q_T, :L),
        )
    end
    full_den = poly_add(
        den_poly_for_conf(Q_cat_R, reg_Q_R, nothing),
        den_poly_for_conf(Q_cat_T, reg_Q_T, :L),
    )
```

- [ ] **Step 3: Add a comment at `make_num_term` and `make_den_term` in `_allosteric_num_den_exprs`**

Locate `make_num_term` (in `_allosteric_num_den_exprs`, around line 1469). Add a comment ABOVE the function:

```julia
    # Numerator: N × Q_cat^(CatN-1) × all reg-site factors at multiplicity.
    # MUST mirror `num_poly_for_conf` in `_count_allosteric_rate_monomials`.
    function make_num_term(N, Q, reg_Qs)
```

Same for `make_den_term`:

```julia
    # Denominator: Q_cat^CatN × all reg-site factors at multiplicity.
    # MUST mirror `den_poly_for_conf` in `_count_allosteric_rate_monomials`.
    function make_den_term(Q, reg_Qs)
```

- [ ] **Step 4: The test suite still has many failures from Task 1's framework changes**

Don't run the suite yet — wait until Task 5 re-measures all spec values.

- [ ] **Step 5: Commit**

```bash
git add src/rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
Patch _count_allosteric_rate_monomials parallel-drift bugs

Two divergences from _allosteric_num_den_exprs after Tasks 3 and 5:
- num_poly_for_conf retained the n_reg == CatN filter Task 3 removed
- full_num always added the L*num_T branch even when t_state_dead

Both produced monomial counts inconsistent with the actual rate
equation, locking in incorrect structural_identifiability_deficit
values for any mechanism with mismatched-multiplicity reg sites or
t_state_dead.

Patched both. Added parallel-constraint comments at the four shared
sites (make_num_term/make_den_term ↔ num_poly_for_conf/den_poly_for_conf)
to remind future authors.
EOF
)"
```

---

### Task 3: Sweep dead `:OnlyT` catalytic scaffolding

After Task 9 (in the prior refactor), the constructor rejects any catalytic kinetic group with state `:OnlyT`. So `_onlyT_syms(m)` always returns an empty `Set{Symbol}()`. Multiple call sites handle a "non-empty t_only_syms" branch that is now unreachable. Remove the dead scaffolding.

(Reg-site `:OnlyT` ligands remain valid. Those are handled by `reg_allo_state` per (site_idx, ligand) — they don't go through `_onlyT_syms`.)

**Files:**
- Modify: `src/rate_eq_derivation.jl` — delete `_onlyT_syms`; remove `t_only_syms` locals and their fast-path branches.

- [ ] **Step 1: Locate and delete `_onlyT_syms`**

In `src/rate_eq_derivation.jl`, find the function definition (around line 1119-1128):

```julia
"""Catalytic-cycle parameter symbols zeroed in the R-state (`:OnlyT` groups)."""
function _onlyT_syms(m::AllostericEnzymeMechanism)
    cm = catalytic_mechanism(m)
    syms = Set{Symbol}()
    for g in kinetic_groups(cm)
        cat_allo_state(m, g) == :OnlyT || continue
        rep = first(steps_in_group(cm, g))
        for s in _group_param_symbols(cm, rep); push!(syms, s); end
    end
    syms
end
```

Delete it entirely.

- [ ] **Step 2: Remove `t_only_syms` from `_dependent_param_exprs`**

In `_dependent_param_exprs(::Type{<:AllostericEnzymeMechanism})` (around line 1219), find:

```julia
    t_only_syms = _onlyT_syms(m)
    r_only_syms = _onlyR_syms(m)
```

Replace with:

```julia
    r_only_syms = _onlyR_syms(m)
```

Then find and remove the line:
```julia
    indep_R = Symbol[p for p in indep_R_all if p ∉ t_only_syms]
```

Replace with:
```julia
    indep_R = collect(indep_R_all)
```

Then find:
```julia
    dep_R = Dict{Symbol, Union{Symbol, Expr}}()
    for (k, v) in dep_R_all
        _expr_references_any(v, t_only_syms) && continue
        dep_R[k] = v
    end
```

Replace with:
```julia
    dep_R = Dict{Symbol, Union{Symbol, Expr}}(dep_R_all)
```

- [ ] **Step 3: Remove `t_only_syms` from `_build_dep_assignments`**

In `_build_dep_assignments` (around line 1325), find:
```julia
    t_only_syms = _onlyT_syms(m)
    r_only_syms = _onlyR_syms(m)
```

Replace with:
```julia
    r_only_syms = _onlyR_syms(m)
```

Find:
```julia
    r_assignments = Expr[]
    for (sym, expr_kd) in sorted_deps
        if _expr_references_any(expr_kd, t_only_syms)
            push!(r_assignments, Expr(:(=), sym, 0))
        else
            push!(r_assignments, Expr(:(=), sym, expr_kd))
        end
    end
```

Replace with:
```julia
    r_assignments = Expr[Expr(:(=), sym, expr_kd) for (sym, expr_kd) in sorted_deps]
```

- [ ] **Step 4: Remove `t_only_syms` from `_kcat_forward(::AllostericEnzymeMechanism, params)`**

Around line 893, find:
```julia
    t_only_syms = _onlyT_syms(m)
    r_only_syms = _onlyR_syms(m)
```

Replace with:
```julia
    r_only_syms = _onlyR_syms(m)
```

Then find the conditional `if isempty(t_only_syms)` block and its else branch (the `_zero_symbols_in_poly(...)` calls); the `else` is now unreachable. Replace the whole conditional with the (formerly true-branch) code path. Specifically, locate and replace:

```julia
    num_R_poly = _zero_symbols_in_poly(_expand_factored_sigma(num_fs), t_only_syms)
    den_R_poly = _zero_symbols_in_poly(_expand_to_poly(denom_terms), t_only_syms)
```

With:
```julia
    num_R_poly = _expand_factored_sigma(num_fs)
    den_R_poly = _expand_to_poly(denom_terms)
```

- [ ] **Step 5: Remove `t_only_syms` from `_allosteric_num_den_exprs`**

Around line 1373, find:
```julia
    t_only_syms = _onlyT_syms(m)
    r_only_syms = _onlyR_syms(m)
```

Replace with:
```julia
    r_only_syms = _onlyR_syms(m)
```

Find the conditional `if isempty(t_only_syms)` block (around lines 1384-1392):

```julia
    if isempty(t_only_syms)
        N_R = _factored_sigma_to_expr(num_fs, cat_params, cat_mets, binding_Ks_r)
        Q_R = _denom_terms_to_expr(denom_terms, cat_params, cat_mets, binding_Ks_r)
    else
        num_r_poly = _zero_symbols_in_poly(_expand_factored_sigma(num_fs), t_only_syms)
        den_r_poly = _zero_symbols_in_poly(_expand_to_poly(denom_terms), t_only_syms)
        N_R = _poly_to_expr(num_r_poly, cat_params, cat_mets, binding_Ks_r)
        Q_R = _poly_to_expr(den_r_poly, cat_params, cat_mets, binding_Ks_r)
    end
```

Replace with (just the now-always-taken branch):
```julia
    N_R = _factored_sigma_to_expr(num_fs, cat_params, cat_mets, binding_Ks_r)
    Q_R = _denom_terms_to_expr(denom_terms, cat_params, cat_mets, binding_Ks_r)
```

- [ ] **Step 6: Remove `t_only_syms` from `_count_allosteric_rate_monomials`**

Around line 1580, find:
```julia
    t_only_syms = _onlyT_syms(m)
    r_only_syms = _onlyR_syms(m)
```

Replace with:
```julia
    r_only_syms = _onlyR_syms(m)
```

Find:
```julia
    N_cat_R = _zero_symbols_in_poly(N_cat_base, t_only_syms)
    Q_cat_R = _zero_symbols_in_poly(Q_cat_base, t_only_syms)
```

Replace with:
```julia
    N_cat_R = N_cat_base
    Q_cat_R = Q_cat_base
```

- [ ] **Step 7: Skip — `_T_rename` `:OnlyT` removal is handled by Task 1 Step 1**

Task 1 Step 1's "New" body for `_T_rename` already drops the `:OnlyT` branch (its first pass only matches `:NonequalRT`). When Task 3 reaches this step, verify the branch is already gone and skip the edit.

If for some reason Task 1 hasn't yet run, then the `:OnlyT` removal can be applied directly here — but the standard execution order is Task 1 first.

- [ ] **Step 8: Run the full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -20
```

Tests will still fail because `expected_n_haldane` values are stale — that's expected, will be fixed in Task 5. Verify that no NEW failures appear from the dead-code removal (i.e., the same set of failures should remain).

- [ ] **Step 9: Commit**

```bash
git add src/rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
Sweep dead :OnlyT catalytic scaffolding from rate equation derivation

After the constructor rejects :OnlyT catalytic groups (R-state-active
convention), _onlyT_syms always returned an empty Set. The framework
had ~30 lines of "if isempty(t_only_syms)" fast paths and zero-
substitution branches that were unreachable.

Deleted _onlyT_syms entirely. Removed t_only_syms locals and their
conditional branches in _dependent_param_exprs, _build_dep_assignments,
_kcat_forward, _allosteric_num_den_exprs, _count_allosteric_rate_monomials.
Removed the :OnlyT branch from _T_rename.

Reg-site :OnlyT ligands remain valid — those are handled per-site via
reg_allo_state, not through _onlyT_syms.
EOF
)"
```

---

### Task 4: Replace `expected_n_haldane`/`expected_n_wegscheider`/`expected_identifiability_deficit` with three categorized fields

`MechanismTestSpec` currently has three test fields tied to `_dependent_param_exprs` output:
- `expected_n_haldane::Int` — actually counts ALL dependent params (Haldanes + mirrors + Wegscheider)
- `expected_n_wegscheider::Int` — separate count
- `expected_identifiability_deficit::Int` — locks in a methodologically-flawed monomial-counting heuristic

Replace with three meaningfully-categorized fields based on RHS structure of each dep expression:
- `expected_n_haldane_constraints::Int` — RHS is an Expr that references `Keq` (true thermodynamic closure for a catalytic cycle)
- `expected_n_mirror_constraints::Int` — RHS is a single Symbol (allosteric `:EqualRT` rename like `K1_T = K1`)
- `expected_n_wegscheider_constraints::Int` — RHS is an Expr that does NOT reference `Keq` (multi-cycle futile-cycle closure)

Drop `expected_identifiability_deficit` entirely. The boolean `expected_is_identifiable` (deficit ≤ 0) stays.

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` — struct definition.
- Modify: `test/test_rate_eq_derivation.jl` — `test_constraint_counting`, `test_identifiability` functions.

- [ ] **Step 1: Update the `MechanismTestSpec` struct**

Open `test/mechanism_definitions_for_test_enzyme_derivation.jl`. Find the struct (around line 13). Replace these three field lines:

Old:
```julia
    expected_n_haldane::Int
    expected_n_wegscheider::Int
    expected_n_independent_params::Int
    expected_identifiability_deficit::Int
    expected_is_identifiable::Bool
```

New:
```julia
    expected_n_haldane_constraints::Int       # RHS references Keq (catalytic-cycle closure)
    expected_n_mirror_constraints::Int        # RHS is a single Symbol (allosteric :EqualRT rename)
    expected_n_wegscheider_constraints::Int   # RHS Expr without Keq (multi-cycle futile-cycle closure)
    expected_n_independent_params::Int
    expected_is_identifiable::Bool
```

Update any default values accordingly (the comments above each field stay close to their fields).

- [ ] **Step 2: Replace `test_constraint_counting` in `test/test_rate_eq_derivation.jl`**

Find the function `test_constraint_counting(spec::MechanismTestSpec)` (around line 439). Replace it with:

```julia
"""Classify a dep expression as Haldane (RHS references Keq), Mirror
(RHS is a single Symbol), or Wegscheider (RHS Expr without Keq)."""
function _classify_dep_expr(expr)
    if expr isa Symbol
        return :mirror
    elseif EnzymeRates._expr_references_any(expr, Set([:Keq]))
        return :haldane
    else
        return :wegscheider
    end
end

function test_constraint_counting(spec::MechanismTestSpec)
    m = spec.mechanism
    @testset "Constraints" begin
        dep_exprs, indep = EnzymeRates._dependent_param_exprs(typeof(m))
        n_haldane = 0
        n_mirror = 0
        n_wegscheider = 0
        for (_, expr) in dep_exprs
            cat = _classify_dep_expr(expr)
            if cat == :haldane
                n_haldane += 1
            elseif cat == :mirror
                n_mirror += 1
            else
                n_wegscheider += 1
            end
        end
        @test n_haldane == spec.expected_n_haldane_constraints
        @test n_mirror == spec.expected_n_mirror_constraints
        @test n_wegscheider == spec.expected_n_wegscheider_constraints
        @test length(indep) == spec.expected_n_independent_params
    end
end
```

- [ ] **Step 3: Update `test_identifiability` to drop the deficit numerical check**

Find `test_identifiability(spec::MechanismTestSpec)` (around line 449). Replace with:

```julia
function test_identifiability(spec::MechanismTestSpec)
    m = spec.mechanism
    @testset "Identifiability" begin
        # The deficit is computed via a monomial-counting heuristic that
        # over-counts identifiable degrees of freedom for factored
        # polynomials (e.g. (Q_R)^catN). Use it only for the boolean
        # is_identifiable check — the magnitude is not biophysically
        # meaningful.
        @test (structural_identifiability_deficit(m) <= 0) ==
              spec.expected_is_identifiable
    end
end
```

- [ ] **Step 4: Don't run tests yet**

The spec entries still use the old field names. They'll be migrated in Task 5.

- [ ] **Step 5: Commit (test framework only — specs follow in Task 5)**

```bash
git add test/test_rate_eq_derivation.jl test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -m "$(cat <<'EOF'
Refactor MechanismTestSpec constraint fields into three categories

Replace expected_n_haldane (which actually counted Haldanes + mirrors)
and expected_n_wegscheider (separately counted multi-cycle closures)
with three categorized fields:
- expected_n_haldane_constraints (RHS references Keq)
- expected_n_mirror_constraints (RHS is a Symbol; :EqualRT rename)
- expected_n_wegscheider_constraints (RHS Expr without Keq)

A future framework regression that swaps a Haldane for a mirror
(same total count, different category) is now caught.

Drop expected_identifiability_deficit. The numerical value relies on
a monomial-counting heuristic that over-counts identifiable degrees
of freedom for factored polynomials (Q^catN expansions in allosteric
mechanisms). Keep the expected_is_identifiable boolean (deficit ≤ 0).

This commit only updates the framework; spec entries are migrated in
the next commit.
EOF
)"
```

---

### Task 5: Re-measure constraint counts for all `MECHANISM_TEST_SPECS` entries

After Tasks 1-4, every spec needs:
- Three categorized constraint counts replacing `expected_n_haldane`/`expected_n_wegscheider`
- `expected_identifiability_deficit` removed
- Possibly updated `expected_n_independent_params` (Task 1's framework fix may add synthesized T-state Haldanes)

Use a measurement script to compute the values from each mechanism programmatically, then paste into spec definitions.

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` — every `push!(specs, MechanismTestSpec(...))` call.

- [ ] **Step 1: Write a measurement script**

Create a temporary file `/tmp/measure_specs.jl`:

```julia
using EnzymeRates
include("test/mechanism_definitions_for_test_enzyme_derivation.jl")

function classify_dep_expr(expr)
    if expr isa Symbol
        return :mirror
    elseif EnzymeRates._expr_references_any(expr, Set([:Keq]))
        return :haldane
    else
        return :wegscheider
    end
end

for spec in MECHANISM_TEST_SPECS
    m = spec.mechanism
    dep_exprs, indep = EnzymeRates._dependent_param_exprs(typeof(m))
    n_h, n_m, n_w = 0, 0, 0
    for (_, expr) in dep_exprs
        c = classify_dep_expr(expr)
        c == :haldane     && (n_h += 1)
        c == :mirror      && (n_m += 1)
        c == :wegscheider && (n_w += 1)
    end
    deficit = structural_identifiability_deficit(m)
    is_id = deficit <= 0
    println("$(spec.name): haldane=$n_h, mirror=$n_m, wegscheider=$n_w, " *
            "indep=$(length(indep)), deficit=$deficit, is_identifiable=$is_id")
end
```

Run it:
```bash
julia --project /tmp/measure_specs.jl
```

This prints all five values for each mechanism. Copy them.

Also: verify there are no `MechanismTestSpec(` constructors outside the main `build_mechanism_test_specs()` function:
```bash
grep -rn "MechanismTestSpec(" test/
```
If any matches appear outside `mechanism_definitions_for_test_enzyme_derivation.jl`, audit those too.

- [ ] **Step 2: Update each spec entry**

For every `push!(specs, MechanismTestSpec(...))` call in `test/mechanism_definitions_for_test_enzyme_derivation.jl`, replace:

Old (example):
```julia
            expected_n_haldane=N1,
            expected_n_wegscheider=N2,
            expected_n_independent_params=N3,
            expected_identifiability_deficit=N4,
            expected_is_identifiable=true,
```

With (measured values from Step 1):
```julia
            expected_n_haldane_constraints=<measured n_h>,
            expected_n_mirror_constraints=<measured n_m>,
            expected_n_wegscheider_constraints=<measured n_w>,
            expected_n_independent_params=<measured indep>,
            expected_is_identifiable=<measured is_identifiable>,
```

Use the boolean from the script's `is_identifiable` column directly — DO NOT guess. The script reports `deficit <= 0` per spec; copy that boolean into `expected_is_identifiable`. Specs whose `is_identifiable=false` should remain `false` even after the framework fix.

- [ ] **Step 3: Run the test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -10
```

Expected: all tests pass.

If any test fails, the measured value was wrong — re-run the measurement script for that spec.

- [ ] **Step 4: Delete the temporary measurement script**

```bash
rm /tmp/measure_specs.jl
```

- [ ] **Step 5: Commit**

```bash
git add test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -m "$(cat <<'EOF'
Re-measure constraint counts across MECHANISM_TEST_SPECS

After Tasks 1-4 (framework fixes + field categorization), every spec
needs updated counts for the three new categorized constraint fields
and an updated independent_params count. Values measured directly
from _dependent_param_exprs output.

expected_identifiability_deficit field removed throughout.
EOF
)"
```

---

### Task 6: Revert PK to `:EqualRT` catalysis

With Task 1's framework fix (synthesized T-state Haldane via dataflow), the previous workaround in PK (`:NonequalRT` catalysis) is no longer needed. `:EqualRT` catalysis is what the original plan intended — cleaner analytical formula, simpler kcat, fewer fitted params.

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` — PK mechanism + analytical formula + spec.

- [ ] **Step 1: Change PK catalysis state**

Locate the PK block (search for `name="PK"`). Find the catalysis SS step in the `@allosteric_mechanism` block:

Old:
```julia
                    [E_PEP_ADP] <--> [E_Pyr_ATP]               :: NonequalRT
```

New:
```julia
                    [E_PEP_ADP] <--> [E_Pyr_ATP]               :: EqualRT
```

- [ ] **Step 2: Simplify the analytical formula**

Replace BOTH the existing param-mapping comment block above `pk_rate_analytical` AND the function body. Find the existing comment block (which currently reads "Param mapping (kinetic-group representative-step convention): K1, K1_T : PEP binding (group 1, NonequalRT) ... k5f, k5f_T : catalysis SS (group 3, NonequalRT) ...") and replace the WHOLE region (comment + function) with the New code below. Failing to replace the comment leaves a stale block contradicting the new function body.

Old (the existing comment block + function — DO NOT just match the function, match the comment block too):

Old:
```julia
        function pk_rate_analytical(params, concs)
            (; K1, K1_T, K3, k5f, k5f_T, K6, K8,
               K_ATP_T_reg1, K_F16BP_reg2,
               L, Keq, Et) = params
            (; PEP, ADP, Pyruvate, ATP, F16BP) = concs

            k5r   = k5f   * K6 * K8 / (Keq * K1   * K3)
            k5r_T = k5f_T * K6 * K8 / (Keq * K1_T * K3)

            Q_cat_R = 1 + PEP/K1   + ADP/K3 + PEP*ADP/(K1   * K3) +
                      Pyruvate/K6  + ATP/K8 + Pyruvate*ATP/(K6 * K8)
            Q_cat_T = 1 + PEP/K1_T + ADP/K3 + PEP*ADP/(K1_T * K3) +
                      Pyruvate/K6  + ATP/K8 + Pyruvate*ATP/(K6 * K8)

            N_R = k5f   * PEP * ADP / (K1   * K3) - k5r   * Pyruvate * ATP / (K6 * K8)
            N_T = k5f_T * PEP * ADP / (K1_T * K3) - k5r_T * Pyruvate * ATP / (K6 * K8)

            ...
```

New:
```julia
        # Param mapping:
        #   K1, K1_T : PEP binding (group 1, NonequalRT)
        #   K3       : ADP binding (group 2, EqualRT)
        #   k5f      : catalysis SS forward rate (group 3, EqualRT)
        #   K6       : Pyruvate release (group 4, EqualRT)
        #   K8       : ATP release (group 5, EqualRT)
        #
        # k5r derives via R-state Haldane: k5r = k5f·K6·K8/(Keq·K1·K3).
        # The framework auto-synthesizes k5r_T because k5r's RHS
        # references K1 (a :NonequalRT symbol with T-rename K1_T):
        #   k5r_T = k5f·K6·K8/(Keq·K1_T·K3).
        # Both Haldanes share the forward k5f — at saturation, forward
        # kcat = catN·k5f (shared between R and T).
        function pk_rate_analytical(params, concs)
            (; K1, K1_T, K3, k5f, K6, K8,
               K_ATP_T_reg1, K_F16BP_reg2,
               L, Keq, Et) = params
            (; PEP, ADP, Pyruvate, ATP, F16BP) = concs

            k5r   = k5f * K6 * K8 / (Keq * K1   * K3)
            k5r_T = k5f * K6 * K8 / (Keq * K1_T * K3)

            Q_cat_R = 1 + PEP/K1   + ADP/K3 + PEP*ADP/(K1   * K3) +
                      Pyruvate/K6  + ATP/K8 + Pyruvate*ATP/(K6 * K8)
            Q_cat_T = 1 + PEP/K1_T + ADP/K3 + PEP*ADP/(K1_T * K3) +
                      Pyruvate/K6  + ATP/K8 + Pyruvate*ATP/(K6 * K8)

            N_R = k5f * PEP * ADP / (K1   * K3) - k5r   * Pyruvate * ATP / (K6 * K8)
            N_T = k5f * PEP * ADP / (K1_T * K3) - k5r_T * Pyruvate * ATP / (K6 * K8)

            Q_reg1_R = 1
            Q_reg1_T = 1 + ATP / K_ATP_T_reg1
            Q_reg2_R = 1 + F16BP / K_F16BP_reg2
            Q_reg2_T = 1

            num_R = N_R * Q_cat_R^3 * Q_reg1_R^2 * Q_reg2_R^4
            num_T = N_T * Q_cat_T^3 * Q_reg1_T^2 * Q_reg2_T^4
            den_R = Q_cat_R^4 * Q_reg1_R^2 * Q_reg2_R^4
            den_T = Q_cat_T^4 * Q_reg1_T^2 * Q_reg2_T^4

            return Et * 4.0 * (num_R + L * num_T) / (den_R + L * den_T)
        end
```

(Diff: `k5f_T` → `k5f` in three places: `k5r_T` derivation, `N_T` flux. Comment block added explaining the framework auto-derivation.)

- [ ] **Step 3: Restore `analytical_kcat_fn`**

In the same `push!(specs, MechanismTestSpec(...))` call, change:

Old:
```julia
            analytical_kcat_fn=nothing,
```

New:
```julia
            analytical_kcat_fn = p -> 4 * p.k5f,
```

- [ ] **Step 4: Re-measure PK's constraint counts and indep_params**

Run the same measurement approach from Task 5 Step 1, just for PK:

```bash
julia --project -e '
using EnzymeRates
include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
m = _spec_by_name("PK").mechanism
dep_exprs, indep = EnzymeRates._dependent_param_exprs(typeof(m))
n_h, n_m, n_w = 0, 0, 0
for (_, expr) in dep_exprs
    if expr isa Symbol
        n_m += 1
    elseif EnzymeRates._expr_references_any(expr, Set([:Keq]))
        n_h += 1
    else
        n_w += 1
    end
end
println("PK: haldane=$n_h, mirror=$n_m, wegscheider=$n_w, indep=$(length(indep))")
'
```

Update the spec values accordingly.

- [ ] **Step 5: Verify analytical match numerically**

```bash
julia --project -e '
using EnzymeRates
using Random
include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
m = _spec_by_name("PK").mechanism
spec = _spec_by_name("PK")
rng = Random.MersenneTwister(42)
indep = (:K1, :K1_T, :K3, :k5f, :K6, :K8, :K_ATP_T_reg1, :K_F16BP_reg2, :L)
keys_t = (indep..., :Keq, :E_total)
vals_t = Tuple(0.1 + 9.9 * rand(rng) for _ in keys_t)
params = NamedTuple{keys_t}(vals_t)
concs = (PEP=0.5, ADP=0.3, Pyruvate=0.1, ATP=0.2, F16BP=0.4)
p_an = merge(params, (Et=params.E_total,))
println("rate_equation = ", rate_equation(m, concs, params))
println("analytical    = ", spec.analytical_rate_fn(p_an, concs))
println("kcat (orig)   = ", EnzymeRates._kcat_forward(m, params))
println("analytical_kcat = ", spec.analytical_kcat_fn(params))
'
```

Expected: rate match to ~10 decimals, kcat match to ~10 decimals.

If the analytical formula doesn't match `rate_equation`, **STOP and report**.

- [ ] **Step 6: Run the full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -5
```

Expected: all tests pass.

- [ ] **Step 7: Commit**

```bash
git add test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -m "$(cat <<'EOF'
Revert PK to :EqualRT catalysis after framework Haldane fix

Task 1's dataflow-aware _dependent_param_exprs auto-derives k5r_T from
the :EqualRT k5f when its RHS references the :NonequalRT K1. The PK
workaround (:NonequalRT catalysis with separate k5f_T) is no longer
needed.

Changes:
- Catalysis SS step: :NonequalRT → :EqualRT
- Analytical formula: drop k5f_T binding; k5r and k5r_T both derive
  from shared k5f via per-state Haldanes
- Restored analytical_kcat_fn = 4 * k5f (closed form, kcat is the
  same in both states because catalysis is :EqualRT)
- One fewer fitted parameter

Mismatched-multiplicity reg sites (ATP::OnlyT mult 2, F16BP::OnlyR
mult 4) remain — PK still tests the symmetric all-reg-sites
contribution from Task 3 of the prior refactor.
EOF
)"
```

---

### Task 7: Move bare-step error inside `_parse_steps_block_with_groups`

The previous refactor's Task 8 left the bare-step error as a post-parse `setdiff` check that names integer kinetic-group IDs. Move it inside the parser where the offending step's full Expr and LineNumberNode are available, so the error message can name the actual step.

**Files:**
- Modify: `src/dsl.jl` — `_parse_steps_block_with_groups` and remove the post-parse check.

- [ ] **Step 1: Locate the parser and add the in-loop error**

Open `src/dsl.jl`. Find `_parse_steps_block_with_groups` (around line 309). The function iterates over step expressions and either extracts a tag or doesn't.

Find the branches around lines 320-355 that handle parenthesized step groups and single steps. For each branch where a step or step group has NO tag and `allow_tag` is `true`, add an explicit error.

Specifically, the pattern is:
```julia
        if arg.head == :tuple
            # Parenthesized step group
            if has_tag_annotation
                # ... extract tag, push to tags ...
            else
                # group without tag
                allow_tag &&
                    error("@allosteric_mechanism: parenthesized step group " *
                          "`$(arg)` is missing `:: <:OnlyR|:EqualRT|:NonequalRT>` " *
                          "annotation. Add `:: <state>` after the closing paren.")
                # ... existing untagged-group handling for plain mechanism ...
            end
        else
            # Single step
            tag = _peel_step_tag!(arg)
            if tag === nothing && allow_tag
                error("@allosteric_mechanism: step `$(arg)` is missing " *
                      "`:: <:OnlyR|:EqualRT|:NonequalRT>` annotation. Add " *
                      "`:: <state>` after the step expression.")
            end
            # ... existing single-step handling ...
        end
```

Use the actual variable names from the existing code (`arg`, `tag`, etc.). The exact location of the if/else matches what's in the file.

- [ ] **Step 2: Remove the post-parse `setdiff` check**

In `_parse_allosteric_mechanism_body` (around line 638-643), find:

```julia
    rxns_expr, group_tags = _parse_steps_block_with_groups(
        cat_steps_block; allow_tag=true,
    )
    tagged = Set{Int}(g for (g, _) in group_tags)
    all_groups = Set{Int}(step.args[4] for step in rxns_expr.args)
    untagged = sort(collect(setdiff(all_groups, tagged)))
    isempty(untagged) ||
        error("@allosteric_mechanism: catalytic step-group(s) $untagged " *
              "missing a ::Tag annotation")
```

Replace with:
```julia
    rxns_expr, group_tags = _parse_steps_block_with_groups(
        cat_steps_block; allow_tag=true,
    )
    # Bare-step rejection now happens inside _parse_steps_block_with_groups.
```

- [ ] **Step 3: Verify no test asserts on the old error message content**

The old post-parse error message named integer group IDs (e.g., `"step-group(s) [3, 5] missing"`). The new in-parser error names a single step's Expr. If any test asserts on the old message format, it'll break.

```bash
grep -rn "step-group\|missing.*annotation\|missing a ::Tag" test/
```

If any matches appear in `@test` assertions (not just comments), update those tests to match the new message OR replace with `@test_throws Exception` (matching error type only).

- [ ] **Step 4: Run the existing bare-step tests in `test_dsl.jl`**

The tests should still pass — they assert `@test_throws Exception`, and the new error is also an Exception. Verify:

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep -E "(@allosteric_mechanism|Test Summary|tests passed)"
```

- [ ] **Step 5: Commit**

```bash
git add src/dsl.jl
git commit -m "$(cat <<'EOF'
Move bare-step error to _parse_steps_block_with_groups (parse-time)

The previous post-parse setdiff check named integer kinetic-group IDs
in the error message, which forced users to manually trace IDs back
to source steps. Moving the check inside the parser means the error
names the actual offending step expression, with Julia line info
available for free.

Trade-off: parser fires on the first bare step it sees; users with
multiple bare-step groups get one error at a time, rather than a
list. In practice this is fine — fix one, recompile, see the next.
EOF
)"
```

---

### Task 8: Constructor robustness bundle

Three improvements to `AllostericEnzymeMechanism` construction and display:
1. Empty-ligand reg site triggers a misleading "all :EqualRT" error. Add an explicit `length(ligands) >= 1` check.
2. No check that `kinetic_group` numbers in the catalytic mechanism are 1..n consecutive. Latent OOB risk for `cat_allo_state(m, g) = states[g]`.
3. `Base.show` filters out `:NonequalRT` ligand annotations, hiding what the user explicitly wrote under dense storage.
4. `_kcat_forward(::AllostericEnzymeMechanism)` has `@assert length(r_keys) == 1` that may opaquely fail for future multi-saturation-pattern mechanisms. Replace with `max` over multiple components.

**Files:**
- Modify: `src/types.jl` — constructor; `Base.show`.
- Modify: `src/rate_eq_derivation.jl` — `_kcat_forward(::AllostericEnzymeMechanism)`.
- Add: `test/test_rate_eq_derivation.jl` or `test/test_types.jl` — `@test_throws` for empty ligand and non-consecutive group.

- [ ] **Step 1: Add empty-ligand check in constructor**

In `src/types.jl`, find the 3-arg constructor. Inside the reg-sites loop (around line 289), BEFORE the all-`:EqualRT` check:

Add:
```julia
        length(ligands) >= 1 ||
            error("Reg site $i: must have at least one ligand; got empty " *
                  "ligand tuple")
```

So the order becomes:
```julia
    for (i, entry) in enumerate(reg_sites)
        ligands, n_reg, reg_allo_states = entry
        ligands isa Tuple && all(l isa Symbol for l in ligands) ||
            error("Reg site $i: ligands must be a Tuple of Symbol")
        length(ligands) >= 1 ||
            error("Reg site $i: must have at least one ligand; got empty " *
                  "ligand tuple")
        n_reg isa Int && n_reg ≥ 1 ||
            error("Reg site $i: multiplicity must be a positive Int")
        # ... rest of validation
```

- [ ] **Step 2: Add kinetic_group consecutiveness check in constructor**

In the same constructor, after `n_groups` is computed (around line 277):

Old:
```julia
    n_groups = length(unique(kinetic_group(cm, i) for i in 1:n_steps(cm)))
    length(cat_allo_states) == n_groups ||
        error("cat_allo_states length $(length(cat_allo_states)) does not " *
              "match catalytic kinetic-group count $n_groups")
```

New:
```julia
    n_groups = length(unique(kinetic_group(cm, i) for i in 1:n_steps(cm)))
    # Validate kinetic_group numbers are 1..n_groups consecutive — the
    # cat_allo_states tuple is indexed by group number, so non-consecutive
    # numbering would cause OOB or wrong-state lookup at runtime.
    observed_groups = sort!(unique(kinetic_group(cm, i) for i in 1:n_steps(cm)))
    observed_groups == collect(1:n_groups) ||
        error("Catalytic mechanism kinetic_group numbers must be 1..n " *
              "consecutive; got $observed_groups")
    length(cat_allo_states) == n_groups ||
        error("cat_allo_states length $(length(cat_allo_states)) does not " *
              "match catalytic kinetic-group count $n_groups")
```

- [ ] **Step 3: Fix `Base.show` to display all states densely (catalytic + reg-site)**

Read the actual `Base.show(io::IO, m::AllostericEnzymeMechanism)` in `src/types.jl` (around line 416-437) before editing — the format is multi-line with `[` `]` brackets for the per-site tagged-ligand summary, NOT the parens / single-line format the previous plan version assumed.

The current format produces output like:
```
AllostericEnzymeMechanism (cat_n=2, 2 reg sites):
  catalytic: <full catalytic mechanism repr>
  reg site 1 (n=2): G6P, Pi [G6P::OnlyT]
  reg site 2 (n=2): ...
```

(The bracket suffix `[G6P::OnlyT]` only appears when there are non-`:NonequalRT` ligands — `:NonequalRT` ligands are silently hidden, which is the bug.)

Two changes:

**Change 3a**: remove the `:NonequalRT` filter on reg-site display (so `:NonequalRT` ligands are no longer hidden).

Find:
```julia
        site_tags = filter(p -> p[1] in site_ligs, collect(tagged))
        non_default = filter(p -> p[2] != :NonequalRT, site_tags)
        if !isempty(non_default)
            print(io, " [")
            print(io, join(("$(n)::$(t)" for (n, t) in non_default), ", "))
            print(io, "]")
        end
```

Replace with:
```julia
        site_state_pairs = filter(p -> p[1] in site_ligs, collect(tagged))
        if !isempty(site_state_pairs)
            print(io, " [")
            print(io, join(("$(n)::$(t)" for (n, t) in site_state_pairs), ", "))
            print(io, "]")
        end
```

**Change 3b**: add a `cat_allo_states` line to the multi-line output, between the header and the `catalytic:` line.

Find (the line starting "AllostericEnzymeMechanism" through the `print(io, "):\n  catalytic: ", cm)` line):
```julia
    print(io, "AllostericEnzymeMechanism (cat_n=", catalytic_multiplicity(m))
    rs = regulatory_sites(m)
    if !isempty(rs)
        print(io, ", ", length(rs), " reg sites")
    end
    print(io, "):\n  catalytic: ", cm)
```

Replace with:
```julia
    print(io, "AllostericEnzymeMechanism (cat_n=", catalytic_multiplicity(m))
    rs = regulatory_sites(m)
    if !isempty(rs)
        print(io, ", ", length(rs), " reg sites")
    end
    print(io, "):\n  cat_allo_states: [")
    n_groups = length(unique(kinetic_group(cm, i) for i in 1:n_steps(cm)))
    print(io, join((string(cat_allo_state(m, g)) for g in 1:n_groups), ", "))
    print(io, "]\n  catalytic: ", cm)
```

After this change, `print(m)` for HK produces (catalytic states now visible, all ligand states shown):
```
AllostericEnzymeMechanism (cat_n=2, 1 reg sites):
  cat_allo_states: [EqualRT, OnlyR, EqualRT, EqualRT, EqualRT, EqualRT]
  catalytic: <full catalytic mechanism repr>
  reg site 1 (n=2): G6P, Pi [G6P::OnlyT, Pi::EqualRT]
```

Add a regression test in `test/test_types.jl`:
```julia
@testset "Base.show displays all dense states" begin
    # Round-trip: every state in the mechanism appears in show output.
    m = pfk_mechanism  # or any mechanism with mixed states
    s = sprint(show, m)
    # Every catalytic state appears
    cm = catalytic_mechanism(m)
    n_groups = length(unique(EnzymeRates.kinetic_group(cm, i)
                             for i in 1:EnzymeRates.n_steps(cm)))
    for g in 1:n_groups
        @test occursin(string(EnzymeRates.cat_allo_state(m, g)), s)
    end
    # No :NonequalRT ligand silently hidden from reg-site display
    for (i, _) in enumerate(EnzymeRates.regulatory_sites(m))
        for lig in EnzymeRates.regulatory_site_ligands(m, i)
            state = EnzymeRates.reg_allo_state(m, i, lig)
            @test occursin("$lig::$state", s)
        end
    end
end
```

- [ ] **Step 4: Replace `@assert length(r_keys) == 1` in `_kcat_forward` with a clear error**

YAGNI applies — current mechanisms (post-Task 6 PK revert) all produce a single saturating-substrate kcat component. The multi-component refactor is a future-only feature. For this round, just upgrade the assertion to a descriptive error.

In `src/rate_eq_derivation.jl`, find the assertion (around line 910). Replace:

Old:
```julia
    r_keys = sort!([k for k in keys(num_R_groups) if haskey(den_R_groups, k)])
    @assert length(r_keys) == 1 "Catalytic mechanism should have exactly 1 kcat component"
```

New:
```julia
    r_keys = sort!([k for k in keys(num_R_groups) if haskey(den_R_groups, k)])
    isempty(r_keys) &&
        error("_kcat_forward: AllostericEnzymeMechanism produced no kcat " *
              "components — saturating-substrate pattern not found in numerator")
    length(r_keys) == 1 ||
        error("_kcat_forward: AllostericEnzymeMechanism with multiple " *
              "saturating-substrate kcat components ($(length(r_keys)) found) " *
              "is not currently supported")
```

Two clear failure modes (no components vs multiple components), each with an actionable message. The body that follows continues to assume `length(r_keys) == 1`.

- [ ] **Step 5: Add error tests**

In `test/test_rate_eq_derivation.jl`, "Allosteric edge cases" testset, add:

```julia
    # Empty ligand list at reg site → constructor error
    cm_simple = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            [E, S] ⇌ [ES]
            [ES] <--> [EP]
            [EP] ⇌ [E, P]
        end
    end
    @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
        cm_simple,
        (2, (:NonequalRT, :EqualRT, :EqualRT)),
        (((), 2, ()),),  # empty ligand tuple
    )
```

- [ ] **Step 6: Run the full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass.

- [ ] **Step 7: Commit**

```bash
git add src/types.jl src/rate_eq_derivation.jl test/test_rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
Tighten AllostericEnzymeMechanism construction and display

Four robustness improvements:
- Empty-ligand reg site now errors with a clear message (previously
  triggered the misleading "all :EqualRT" error via vacuous all/[])
- kinetic_group numbers in catalytic mechanism must be 1..n
  consecutive; non-consecutive numbering would cause OOB or wrong-state
  lookup since cat_allo_states is indexed by group number
- Base.show now displays ALL ligand states under dense storage,
  no longer hiding :NonequalRT (the previous "non-default" filter
  reflected sparse-storage semantics that no longer apply)
- _kcat_forward replaces the @assert length(r_keys) == 1 with a
  clear error message documenting the single-component limitation

Regression tests for the empty-ligand and non-consecutive-group cases
added in Allosteric edge cases.
EOF
)"
```

---

### Task 9: Add `:EqualRT` to `_expand_add_allosteric_regulator` for existing sites

The enumeration cannot reach m_all-shape mechanisms (mixed-state reg sites with at least one `:EqualRT` ligand) because `_expand_add_allosteric_regulator` enumerates only `(:OnlyR, :OnlyT, :NonequalRT)` for new ligands. Extend with `:EqualRT` — but only when adding to an EXISTING site that already has at least one non-`:EqualRT` ligand (a brand-new site with a single `:EqualRT` ligand would error in the constructor).

**Files:**
- Modify: `src/mechanism_enumeration.jl` — `_expand_add_allosteric_regulator`.
- Add: `test/test_mechanism_enumeration.jl` — test that m_all-shape is reachable from enumeration.

- [ ] **Step 1: Extend the loop**

Open `src/mechanism_enumeration.jl`. Find `_expand_add_allosteric_regulator` (around line 1909). Replace:

Old:
```julia
    results = AllostericMechanismSpec[]
    for reg in new_regs
        n_sites = length(spec.allosteric_reg_sites)
        for tag in (:OnlyR, :OnlyT, :NonequalRT)
            for site_idx in 0:n_sites
                # ... existing body
            end
        end
    end
    results
```

New:
```julia
    results = AllostericMechanismSpec[]
    for reg in new_regs
        n_sites = length(spec.allosteric_reg_sites)
        # Enumerate non-:EqualRT states for any site (new or existing)
        for tag in (:OnlyR, :OnlyT, :NonequalRT)
            for site_idx in 0:n_sites
                # ... existing body
            end
        end
        # Enumerate :EqualRT only for existing sites where at least one
        # ligand is already non-:EqualRT (single-ligand or all-:EqualRT
        # site cancels identically — constructor would reject).
        for site_idx in 1:n_sites
            existing_ligs = spec.allosteric_reg_sites[site_idx]
            any(get(spec.reg_ligand_tags, l, :NonequalRT) != :EqualRT
                for l in existing_ligs) || continue
            new_sites = deepcopy(spec.allosteric_reg_sites)
            new_mults = copy(spec.allosteric_multiplicities)
            new_lig_tags = copy(spec.reg_ligand_tags)
            push!(new_sites[site_idx], reg)
            new_lig_tags[reg] = :EqualRT
            delta_cost = _allo_lig_tag_delta(:EqualRT, :EqualRT) + 1
            push!(results, AllostericMechanismSpec(
                spec.base, spec.catalytic_n,
                new_sites, new_mults,
                copy(spec.group_tags), new_lig_tags,
                spec.param_count + delta_cost))
        end
    end
    results
```

(Reuse the existing inner-body structure; this just adds a parallel loop for `:EqualRT` ligands at existing sites.)

- [ ] **Step 2: Add a regression test**

Open `test/test_mechanism_enumeration.jl`. Find a relevant testset (e.g., near other `_expand_add_allosteric_regulator` tests). Add:

```julia
@testset "EqualRT ligand reachable at existing reg site" begin
    # Set up a spec with a single regulator at one site, non-:EqualRT.
    # Then expand to add a SECOND regulator at the SAME site as :EqualRT.
    # Verify a result spec exists with both ligands at site 1, the
    # second tagged :EqualRT.
    rxn = @enzyme_reaction begin
        substrates: S
        products:   P
        regulators: I, J
    end
    base = first(EnzymeRates.init_mechanisms(rxn))
    # Promote to allosteric with one regulator I::OnlyR at site 1
    allo_specs = EnzymeRates._expand_to_allosteric(base, rxn)
    seed = first(filter(s -> haskey(s.reg_ligand_tags, :I) &&
                              s.reg_ligand_tags[:I] == :OnlyR, allo_specs))
    expanded = EnzymeRates._expand_add_allosteric_regulator(seed, rxn)
    # Find a result where J is at site 1 with :EqualRT
    target = findfirst(expanded) do s
        get(s.reg_ligand_tags, :J, nothing) == :EqualRT &&
            (:J in s.allosteric_reg_sites[1])
    end
    @test target !== nothing
end
```

(Adjust to your codebase's actual API for constructing seed specs — the snippet may need touch-up.)

- [ ] **Step 3: Run the test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass.

- [ ] **Step 4: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "$(cat <<'EOF'
Enumerate :EqualRT ligands at existing reg sites

_expand_add_allosteric_regulator previously enumerated only
(:OnlyR, :OnlyT, :NonequalRT) for new ligands. m_all-shape mechanisms
(mixed-state reg site like R1::NonequalRT + R2::EqualRT) were
unreachable from auto-enumeration, only constructible by hand.

Extend with :EqualRT — but only when adding to an existing site that
already has at least one non-:EqualRT ligand (a brand-new site with a
single :EqualRT ligand would cancel identically and be rejected by
the constructor).

dedup! handles multi-path equivalents correctly via its existing
ligand-name sort within each site.
EOF
)"
```

---

### Task 10: Rename remaining `tag` → `allo_state` in function names

Task 7 of the prior refactor renamed the public-API accessors but missed several internal function names. Complete the rename for consistency.

Functions to rename:
- `_expand_change_group_tag` → `_expand_change_allo_state`
- `_allo_tag_delta` → `_allo_state_delta`
- `_allo_lig_tag_delta` → `_allo_lig_state_delta`
- `_format_tag_set` → `_format_state_set`
- `_ALLOSTERIC_REG_TAGS` → `_ALLOSTERIC_REG_STATES`

**Files:**
- Modify: `src/mechanism_enumeration.jl`, `src/dsl.jl`, `test/test_mechanism_enumeration.jl`, `.claude/CLAUDE.md` — every reference to the old names.

- [ ] **Step 1: Find all occurrences**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
grep -rn "_expand_change_group_tag\|_allo_tag_delta\|_allo_lig_tag_delta\|_format_tag_set\|_ALLOSTERIC_REG_TAGS" src/ test/ .claude/
```

Record the per-name match count. Each name should have at least one match. If any name has ZERO matches, the symbol no longer exists in the codebase — drop it from the rename list rather than apply a vacuous rename.

- [ ] **Step 2: Apply the renames (only those with matches)**

Use editor find-and-replace for each name confirmed to have matches:
- `_expand_change_group_tag` → `_expand_change_allo_state`
- `_allo_tag_delta` → `_allo_state_delta`
- `_allo_lig_tag_delta` → `_allo_lig_state_delta`
- `_format_tag_set` → `_format_state_set`
- `_ALLOSTERIC_REG_TAGS` → `_ALLOSTERIC_REG_STATES`

Verify post-rename:
```bash
grep -rn "_expand_change_group_tag\|_allo_tag_delta\|_allo_lig_tag_delta\|_format_tag_set\|_ALLOSTERIC_REG_TAGS" src/ test/ .claude/
```

Expected: zero matches.

- [ ] **Step 3: Update CLAUDE.md mentions**

Find references to the old function names in `.claude/CLAUDE.md` and update them similarly.

- [ ] **Step 4: Run the test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "$(cat <<'EOF'
Complete tag → allo_state rename in internal function names

Previous refactor renamed the public API (group_tag → cat_allo_state,
regulatory_ligand_tag → reg_allo_state, group_tags → cat_allo_states,
reg_ligand_tags → reg_allo_states) but left these internal names:
- _expand_change_group_tag → _expand_change_allo_state
- _allo_tag_delta → _allo_state_delta
- _allo_lig_tag_delta → _allo_lig_state_delta
- _format_tag_set → _format_state_set
- _ALLOSTERIC_REG_TAGS → _ALLOSTERIC_REG_STATES

Pure rename; behavior unchanged.
EOF
)"
```

---

### Task 11: Documentation/comment updates

Bundle of minor doc/comment fixes from the review:

- CLAUDE.md `:NonequalRT` "default for absent entries" wording → distinguish dense type-parameter storage from the internal sparse Dict in `AllostericMechanismSpec`.
- Outdated MWC formula header in `src/rate_eq_derivation.jl:1086-1093` (describes the old `n_reg < cat_n` rule).
- `vtype` test comment in `test_rate_eq_derivation.jl` says "T-state numerator literally zero" — undersells Task 5's elision.
- m_all spec comment doesn't enumerate the 12 `expected_n_independent_params`.
- PK spec comment doesn't reflect the post-Task-6 `:EqualRT` catalysis.
- `(smoke)` testset in `test/test_dsl.jl` now hosts substantive validation tests — rename to `(parsing & validation)`.
- `AllostericMechanismSpec` docstring still uses "tag" terminology — clarify that internal Dict storage uses sparse default-`:NonequalRT` semantics, distinct from the dense type-parameter representation.

**Files:**
- Modify: `.claude/CLAUDE.md`, `src/rate_eq_derivation.jl`, `src/mechanism_enumeration.jl`, `test/test_dsl.jl`, `test/test_rate_eq_derivation.jl`, `test/mechanism_definitions_for_test_enzyme_derivation.jl`.

- [ ] **Step 1: Update CLAUDE.md**

Open `.claude/CLAUDE.md`. Find the "Allosteric state taxonomy" section.

First, fix the misleading "default" wording. Find any line saying:
```
:NonequalRT (default) — independent R and T symbols
```
Replace with:
```
:NonequalRT — independent R and T symbols
```
(Drop the parenthetical "(default)". Under dense storage there is no default.)

Then add a clarifying paragraph distinguishing dense vs sparse storage:

```markdown
- The `AllostericEnzymeMechanism` type-parameter storage is **dense** — every catalytic
  kinetic group has an explicit entry in `cat_allo_states`, every regulator ligand has
  an explicit entry in `reg_allo_states`. There is no defaulting.
- The internal `AllostericMechanismSpec` (in `mechanism_enumeration.jl`) uses **sparse**
  Dict storage where absent entries are interpreted as `:NonequalRT` during the dense
  conversion in the `AllostericEnzymeMechanism(spec)` constructor. This is an
  enumeration internal — DSL and constructor users see only the dense form.
```

- [ ] **Step 2: Update outdated MWC formula header**

Open `src/rate_eq_derivation.jl`. Find the comment block around line 1086-1093 (the `_t_state_dead` docstring or nearby formula explanation). Find any reference to "Sites with `n_reg_i < cat_n` appear only in the denominator" or similar. Update to reflect the current behavior:

Old:
```
> Regulatory sites with `n_reg_i < cat_n` appear only in the denominator. Sites with `n_reg_i == cat_n` appear in both numerator and denominator.
```

New:
```
> Regulatory sites contribute to BOTH numerator and denominator at their multiplicity, regardless of whether `n_reg_i` matches `cat_n`.
```

- [ ] **Step 3: Update vtype test comment**

In `test/test_rate_eq_derivation.jl`, find the `vtype` test (around line 935-950). Update the comment to evergreen wording (no "post-Task" temporal references — CLAUDE.md forbids historical context in comments):

Old:
```
# T-state numerator literally zero: rate ∝ 1/(1+L) at large L
```

New:
```
# T-state numerator branch is elided when t_state_dead (any :OnlyR catalytic group);
# rate is E_total · catN · num_R / (Q_R^catN + L · Q_T^catN). At large L, the T-state
# enzyme mass dominates the denominator → rate ∝ 1/(1+L).
```

- [ ] **Step 4: Update m_all spec comment**

In `test/mechanism_definitions_for_test_enzyme_derivation.jl`, find the m_all spec. Add an enumeration of the independent params:

```
# Independent parameters (12): K1, K1_T, K3, k5f, k5f_T, K6, K6_T, K8,
# K_R1_reg1, K_R1_T_reg1, K_R2_reg1, L
```

- [ ] **Step 5: Update PK spec comment to reflect :EqualRT catalysis**

In the PK spec block (post-Task 6), update the descriptive comment to reflect `:EqualRT` catalysis:

```
# Catalysis SS step is :EqualRT; k5r and k5r_T both derive from the
# shared k5f via per-state Haldanes (R-state uses K1, T-state uses K1_T).
# Independent parameters (9): K1, K1_T, K3, k5f, K6, K8,
# K_ATP_T_reg1, K_F16BP_reg2, L
```

(Adjust counts to match Task 5's measured values.)

- [ ] **Step 6: Rename `(smoke)` testset in `test_dsl.jl`**

Open `test/test_dsl.jl`. Find:
```julia
    @testset "@allosteric_mechanism (smoke)" begin
```

Rename to:
```julia
    @testset "@allosteric_mechanism (parsing & validation)" begin
```

- [ ] **Step 7: Update `AllostericMechanismSpec` docstring**

In `src/mechanism_enumeration.jl`, find the `AllostericMechanismSpec` docstring (around line 46-52). Add a clarifying note:

```
> Note: this struct uses sparse Dict storage internally (where absent entries default
> to `:NonequalRT`). When converted via `AllostericEnzymeMechanism(spec)`, the type
> parameters become dense — every catalytic kinetic group has an explicit
> `cat_allo_states` entry, every regulator ligand has an explicit `reg_allo_states`
> entry.
```

- [ ] **Step 8: Run the test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -3
```

Expected: all tests pass.

- [ ] **Step 9: Commit**

```bash
git add -A
git commit -m "$(cat <<'EOF'
Documentation and comment updates from review

Bundled doc fixes:
- CLAUDE.md clarifies dense (type-parameter) vs sparse (enumeration
  internal Dict) allo-state storage
- Outdated MWC formula header reflects Task 3's all-reg-sites
  contribution (no more n_reg < cat_n filter)
- vtype test comment reflects Task 5's elision of L*N_T
- m_all and PK spec comments enumerate their independent parameters
- (smoke) testset renamed to (parsing & validation) since it now hosts
  substantive validation tests
- AllostericMechanismSpec docstring distinguishes its internal sparse
  Dict storage from the dense type-parameter representation
EOF
)"
```

---

### Task 12: Remove redundant `:OnlyR sub + :OnlyT prod` test

The constructor-error test for the `:OnlyR sub + :OnlyT prod` combination is structurally identical to the `:OnlyT prod` test alone — both fail at the constructor's `:OnlyT` rejection. The combined test adds no coverage.

**Files:**
- Modify: `test/test_rate_eq_derivation.jl`.

- [ ] **Step 1: Remove the redundant test**

In `test/test_rate_eq_derivation.jl`, find the comment `# :OnlyR substrate + :OnlyT product → constructor error` (around line 995). Delete the entire `@test_throws Exception eval(...)` block including the comment.

- [ ] **Step 2: Run the test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -3
```

Expected: all tests pass with one fewer test in the count.

- [ ] **Step 3: Commit**

```bash
git add test/test_rate_eq_derivation.jl
git commit -m "$(cat <<'EOF'
Remove redundant :OnlyR sub + :OnlyT prod constructor-error test

The test was structurally identical to the :OnlyT prod test alone —
both fire the same constructor :OnlyT rejection. The combined test
added no coverage. The independent rules (single-:OnlyT errors and
mixed-:OnlyR-:OnlyT errors are subsumed) are still tested via the
single-tag :OnlyT cases.
EOF
)"
```

---

## Self-Review

Spec coverage check:

| Issue | Task |
|---|---|
| #1 (Critical: `_count_allosteric_rate_monomials` parallel-drift) | 2 |
| #2 (Critical: `_dependent_param_exprs` mixed-tag bug) | 1 |
| #3 (Dead `:OnlyT` catalytic branches) | 3 |
| #4 (`Base.show` filters `:NonequalRT`) | 8 |
| #5 (Empty-ligand error message) | 8 |
| #6 (Dead T-mirrors when t_state_dead) | 1 |
| #7 (Bare-step error placement) | 7 |
| #8 (`_expand_add_allosteric_regulator` `:EqualRT` gap) | 9 |
| #9 ("tag" leftover function names) | 10 |
| #10 (kinetic_group consecutiveness check) | 8 |
| #11 (CLAUDE.md `:NonequalRT` "default" docstring) | 11 |
| #12 (`AllostericMechanismSpec` docstring) | 11 |
| #13 (Outdated MWC formula header) | 11 |
| #14 (Redundant test) | 12 |
| #15 (`expected_n_haldane` field name) | 4 |
| #16 (PK plan-vs-implementation comment) | 11 |
| #17 (m_all `expected_n_independent_params=12` comment) | 11 |
| #18 (`@assert` may break for PK) | 8 |
| #19 (`(smoke)` testset name) | 11 |
| #20 (No iso-only `:OnlyR` end-to-end test) | covered by existing `vtype` test (`[ES] <--> [EP] :: OnlyR` is iso-only) |

Type-consistency check:
- `cat_allo_state(m, g)` (singular accessor) — used in Task 1.
- `_expand_change_allo_state` (post-rename) — used consistently across Tasks 9, 10.
- `expected_n_haldane_constraints`, `expected_n_mirror_constraints`, `expected_n_wegscheider_constraints` — three categorized fields used in Tasks 4, 5, 6.

Note: Task 1 covers BOTH issue #2 (mixed-tag bug) AND issue #6 (skip dead T-mirrors when t_state_dead) — they live in the same function and are bundled.

Issue #20 (iso-only `:OnlyR`): the existing `vtype` test in `test_rate_eq_derivation.jl` already covers this. Its catalysis SS step is `[ES] <--> [EP] :: OnlyR` — both sides are bound enzyme forms, no metabolite. That IS an iso-only `:OnlyR` step. End-to-end coverage already exists.

Order constraints:
- Task 1 must precede Task 3 (Task 1 Step 1's "New" body for `_T_rename` removes the `:OnlyT` branch, satisfying what Task 3 Step 7 was going to do; Task 3 Step 7 then becomes a no-op or skip).
- Task 1 must precede Task 5 (Task 1 changes dep counts that Task 5 measures).
- Task 1 must precede Task 6 (Task 6 reverts PK based on Task 1's framework fix).
- Task 4 must precede Task 5 (struct fields change before bulk update).
- Task 9 must precede Task 10 (Task 9 references `_allo_lig_tag_delta`; Task 10 renames it to `_allo_lig_state_delta`. If Task 10 runs first, Task 9's prescribed code references a non-existent symbol).
- Task 1 Step 5 is the gate: the `m_mixed` regression test MUST pass before proceeding to subsequent tasks. If it fails, Task 1's framework fix is incomplete and downstream tasks (5, 6) will produce wrong results.
- Other tasks (7, 8, 11, 12) are independent and can run in any order.
