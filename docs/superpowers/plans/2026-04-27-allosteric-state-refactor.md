# Allosteric State Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refactor `AllostericEnzymeMechanism` to (a) require explicit dense `cat_allo_states`/`reg_allo_states` declarations replacing the sparse default-`:NonequalRT` tag dicts, (b) drop the `n_reg == CatN` filter so all reg sites contribute to the rate-equation numerator at their multiplicity, (c) broaden `_t_state_dead` to fire for any `:OnlyR` catalytic group (not just substrate-binding), (d) drop the dead `L * N_T` numerator branch when t_state_dead, (e) error on construction for any `:OnlyT` catalytic group (R-state-active convention), and (f) add three new mechanisms (PK, m_all, m_OnlyR_prod) with hand-derived analytical formulas to `MECHANISM_TEST_SPECS`.

**Architecture:** All allosteric states are dense and validated at construction. The R-state is the convention for "active". `:OnlyT` regulator ligands stay valid (regulators don't gate the catalytic cycle). Rate-equation simplifications are scoped to `_allosteric_num_den_exprs` and `_kcat_forward`. Field rename and DSL changes ripple through `src/types.jl`, `src/dsl.jl`, `src/mechanism_enumeration.jl`, `src/rate_eq_derivation.jl`, and 4 test files.

**Tech Stack:** Julia, EnzymeRates.jl, type-parameter-based mechanism encoding.

---

## Files Touched

- `src/types.jl` — `AllostericEnzymeMechanism` struct, accessors `group_tag` / `regulatory_ligand_tag`, constructor validation.
- `src/dsl.jl` — `@allosteric_mechanism` macro, allow_tag handling, regulator parsing.
- `src/mechanism_enumeration.jl` — `AllostericMechanismSpec` struct, `_expand_to_allosteric` (drop `:OnlyT` from catalytic enumeration), all expansion-move callers that touch `group_tags` / `reg_ligand_tags`.
- `src/rate_eq_derivation.jl` — `_t_state_dead` simplification, `make_num_term` filter removal, `_allosteric_num_den_exprs` drop-`L*N_T` branch, `_build_dep_assignments` skip dead T-assignments.
- `test/test_enzyme_derivation.jl` → renamed to `test/test_rate_eq_derivation.jl` — remove PFK/HK standalone testsets, convert `onlyT_sub` to error test, add new error tests, remove `*_broken` field references.
- `test/runtests.jl` — update include path post-rename.
- `test/mechanism_definitions_for_test_enzyme_derivation.jl` — remove `*_broken` fields, migrate all 6 `@allosteric_mechanism` definitions to dense state declarations, add `m_PK`, `m_all`, `m_OnlyR_prod`.
- `test/test_dsl.jl` — migrate 4 `@allosteric_mechanism` test definitions.
- `test/test_types.jl` — migrate 2 `@allosteric_mechanism` test definitions.
- `.claude/CLAUDE.md` — update terminology (`tag` → `allo state`).

---

### Task 1: Remove `factored_num_broken` and `factored_denom_broken` placeholder flags

These two `Bool=false` flag fields in `MechanismTestSpec` were added as escape hatches but every call site sets them to `false`. The `@test_broken` branch in `test_factored_form` is dead code.

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`
- Modify: `test/test_enzyme_derivation.jl` (lines 626–637)

- [ ] **Step 1: Remove the two field declarations**

In `test/mechanism_definitions_for_test_enzyme_derivation.jl`, locate the `MechanismTestSpec` struct (around line 51–54) and remove:

```julia
    factored_num_broken::Bool = false
    factored_denom_broken::Bool = false
```

So the struct ends with `expected_factored_denom::Union{String,Nothing} = nothing` as the last field (no comma after `nothing`).

- [ ] **Step 2: Remove all `factored_num_broken=false` / `factored_denom_broken=false` lines from spec definitions**

In `test/mechanism_definitions_for_test_enzyme_derivation.jl`, search for `factored_num_broken=false` and `factored_denom_broken=false`. Every occurrence (8 lines) is `=false`. Remove each line. Verify with `grep -n "factored_num_broken\|factored_denom_broken" test/mechanism_definitions_for_test_enzyme_derivation.jl`. Expected: zero matches after edit.

- [ ] **Step 3: Replace the conditional in `test_factored_form`**

In `test/test_enzyme_derivation.jl`, find the function `test_factored_form` (around line 612). Replace the inner block (lines 624–638) that branches on `spec.factored_num_broken` / `spec.factored_denom_broken`. Old code:

```julia
        if num_str !== nothing && denom_str !== nothing
            if has_num
                if spec.factored_num_broken
                    @test_broken num_str == spec.expected_factored_num
                else
                    @test num_str == spec.expected_factored_num
                end
            end
            if has_denom
                if spec.factored_denom_broken
                    @test_broken denom_str == spec.expected_factored_denom
                else
                    @test denom_str == spec.expected_factored_denom
                end
            end
        end
```

New code:

```julia
        if num_str !== nothing && denom_str !== nothing
            has_num && @test num_str == spec.expected_factored_num
            has_denom && @test denom_str == spec.expected_factored_denom
        end
```

- [ ] **Step 4: Run the full test suite**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_enzyme_derivation.jl
git commit -m "Remove vestigial *_broken flags from MechanismTestSpec"
```

---

### Task 2: Rename `test/test_enzyme_derivation.jl` → `test/test_rate_eq_derivation.jl`

Mirrors the source-file naming convention `src/rate_eq_derivation.jl` with a `test_` prefix.

**Files:**
- Rename: `test/test_enzyme_derivation.jl` → `test/test_rate_eq_derivation.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Rename the file via `git mv`**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
git mv test/test_enzyme_derivation.jl test/test_rate_eq_derivation.jl
```

- [ ] **Step 2: Update the include path in `test/runtests.jl`**

Locate the line `include("test_enzyme_derivation.jl")` (or similar) in `test/runtests.jl` and change to `include("test_rate_eq_derivation.jl")`.

- [ ] **Step 3: Run the full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass.

- [ ] **Step 4: Commit**

```bash
git add test/test_rate_eq_derivation.jl test/runtests.jl
git commit -m "Rename test/test_enzyme_derivation.jl to test/test_rate_eq_derivation.jl

Aligns with src/rate_eq_derivation.jl file naming."
```

---

### Task 3: Drop `n_reg == CatN` filter from `make_num_term`

Currently `make_num_term` filters reg sites whose multiplicity differs from `CatN`, so they appear only in the denominator. Per Denis: all reg sites contribute symmetrically. PFK and HK have all reg sites at `n_reg == CatN`, so this change has no effect on existing tests.

**Files:**
- Modify: `src/rate_eq_derivation.jl` (around line 1423)

- [ ] **Step 1: Remove the filter line**

Locate `make_num_term` inside `_allosteric_num_den_exprs` (around line 1423). Remove the line `n_reg == CatN || continue` so the loop body always runs:

Before:

```julia
    function make_num_term(N, Q, reg_Qs)
        factors = Any[N]
        CatN > 1 && push!(factors, _power_expr(Q, CatN - 1))
        for i in eachindex(RS)
            n_reg = regulatory_site_multiplicity(m, i)
            n_reg == CatN || continue
            push!(factors, _power_expr(reg_Qs[i], n_reg))
        end
        _nest_binary(:*, factors)
    end
```

After:

```julia
    function make_num_term(N, Q, reg_Qs)
        factors = Any[N]
        CatN > 1 && push!(factors, _power_expr(Q, CatN - 1))
        for i in eachindex(RS)
            push!(factors, _power_expr(reg_Qs[i],
                                       regulatory_site_multiplicity(m, i)))
        end
        _nest_binary(:*, factors)
    end
```

Now `make_num_term` mirrors `make_den_term` exactly (modulo the catalytic-Q exponent of `CatN-1` vs `CatN`).

- [ ] **Step 2: Run the full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass — PFK and HK already have `n_reg == CatN` for every reg site so the filter was a no-op for them.

- [ ] **Step 3: Commit**

```bash
git add src/rate_eq_derivation.jl
git commit -m "Drop n_reg == CatN filter from make_num_term

Symmetric with make_den_term — all reg sites contribute to the rate
equation numerator at their multiplicity. Mismatched-mult sites
(e.g., PK with mult 2 + mult 4 reg sites) now correctly carry the
mult-2 factor through the numerator."
```

---

### Task 4: Simplify `_t_state_dead` — fire on any `:OnlyR` catalytic group

The current implementation has a substrate-only filter that misses product-binding `:OnlyR` groups. Per Denis: any `:OnlyR` catalytic group breaks the T-state cycle (substrate binding, product binding, isomerization, or catalysis SS — all four cases). The Cha polynomial-zeroing alone produces a non-zero `N_T` for substrate- and product-binding cases at chemical equilibrium, which is non-physical. Forcing `N_T = 0` whenever any `:OnlyR` catalytic group exists fixes both.

**Files:**
- Modify: `src/rate_eq_derivation.jl` (around lines 1086–1107)

- [ ] **Step 1: Replace `_t_state_dead` body with the simplified version**

Locate the function `_t_state_dead(m::AllostericEnzymeMechanism)` in `src/rate_eq_derivation.jl`. Replace the entire body:

```julia
"""
The T-state catalytic cycle cannot close — and therefore both forward
and reverse net flux vanish — when any `:OnlyR` kinetic group is
present. The Cha polynomial-zeroing approach kills only one half of
the catalytic flux (forward for substrate-OnlyR, reverse for
product-OnlyR), leaving the other half non-zero at chemical
equilibrium. Forcing `N_T = 0` ensures Haldane consistency. Used by
both `rate_equation` (via `_allosteric_num_den_exprs`) and
`_kcat_forward`.
"""
function _t_state_dead(m::AllostericEnzymeMechanism)
    cm = catalytic_mechanism(m)
    any(group_tag(m, g) == :OnlyR for g in kinetic_groups(cm))
end
```

The previous metabolite-LHS check (`isempty(mets_lhs) || any(met in sub_set for met in mets_lhs)`) is removed. Catalysis and isomerization-`:OnlyR` already produce `N_T = 0` from polynomial zeroing, so forcing it explicitly is redundant but harmless. Substrate-binding-`:OnlyR` and product-binding-`:OnlyR` need the explicit forcing.

- [ ] **Step 2: Run the full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass — the simplification is a strict superset of the previous behavior; everything that was `t_state_dead = true` before still is.

- [ ] **Step 3: Commit**

```bash
git add src/rate_eq_derivation.jl
git commit -m "Simplify _t_state_dead: fire on any :OnlyR catalytic group

The previous metabolite-LHS filter only caught substrate-binding and
isomerization :OnlyR groups. Product-binding :OnlyR was missed,
producing non-physical reverse flux at chemical equilibrium via the
surviving k_cat·S/K_S term. Broaden detection to any :OnlyR group;
catalysis and isomerization cases were already at 0 from polynomial
zeroing so the explicit force is harmless there."
```

---

### Task 5: Drop `L * N_T` numerator term and dead T-assignments when t_state_dead

Currently `_allosteric_num_den_exprs` returns a numerator of `CatN * (num_R + L * num_T)` even when `_t_state_dead = true`. Since `N_T = 0` in that case, the term expands at runtime to `L * 0 * Q_T_polynomial * Q_reg_T^...` — all the multiplicative factors get computed and then discarded by floating-point multiplication. The display string also carries the dead branch verbatim. Drop it. Also skip the t_assignments (T-state Haldane and `:EqualRT` catalytic mirrors `K1_T = K1`, etc.) since none are referenced once the L*N_T term is gone.

**Files:**
- Modify: `src/rate_eq_derivation.jl` (around lines 1493 in `_allosteric_num_den_exprs`, and `_build_allosteric_rate_body`, and `rate_equation_string` around line 1497–1540)

- [ ] **Step 1: Modify `_allosteric_num_den_exprs` to skip the L*N_T term when t_state_dead**

Locate the final return at the bottom of `_allosteric_num_den_exprs` (around line 1493). Replace:

```julia
    num_R = make_num_term(N_R, Q_R, reg_Q_R)
    den_R = make_den_term(Q_R, reg_Q_R)
    num_T = make_num_term(N_T, Q_T, reg_Q_T)
    den_T = make_den_term(Q_T, reg_Q_T)

    :($(CatN) * ($(num_R) + L * $(num_T))), :($(den_R) + L * $(den_T))
```

with:

```julia
    num_R = make_num_term(N_R, Q_R, reg_Q_R)
    den_R = make_den_term(Q_R, reg_Q_R)
    den_T = make_den_term(Q_T, reg_Q_T)

    if _t_state_dead(m)
        # T-state cycle broken: N_T = 0, so drop the L*num_T term
        # entirely (skip dead numerator branch). Q_T still contributes
        # to denominator as enzyme mass.
        :($(CatN) * $(num_R)), :($(den_R) + L * $(den_T))
    else
        num_T = make_num_term(N_T, Q_T, reg_Q_T)
        :($(CatN) * ($(num_R) + L * $(num_T))), :($(den_R) + L * $(den_T))
    end
```

- [ ] **Step 2: Skip T-assignments in `_build_allosteric_rate_body` when t_state_dead**

Locate `_build_allosteric_rate_body` (around line 1497). Modify the t_assignments handling to mirror what `_kcat_forward` already does. Old code:

```julia
function _build_allosteric_rate_body(M_type::Type{<:AllostericEnzymeMechanism})
    full_num, full_den = _allosteric_num_den_exprs(M_type)
    rate_expr = :(E_total * ($full_num) / ($full_den))

    r_assignments, t_assignments = _build_dep_assignments(M_type, K -> :(inv($K)))

    _, indep = _dependent_param_exprs(M_type)
    hw_params = (indep..., :Keq, :E_total)
    mets = metabolites(M_type())

    Expr(:block,
        _destructuring_expr(hw_params, :params),
        _destructuring_expr(mets, :concs),
        r_assignments...,
        t_assignments...,
        :(return $rate_expr))
end
```

New code:

```julia
function _build_allosteric_rate_body(M_type::Type{<:AllostericEnzymeMechanism})
    full_num, full_den = _allosteric_num_den_exprs(M_type)
    rate_expr = :(E_total * ($full_num) / ($full_den))

    r_assignments, t_assignments_ = _build_dep_assignments(M_type, K -> :(inv($K)))
    # When the T-state cycle is broken, t_assignments (T-state Haldanes
    # and :EqualRT catalytic mirrors K_T = K) become dead code — they're
    # only referenced from the L*num_T branch, which is now elided.
    t_assignments = _t_state_dead(M_type()) ? Expr[] : t_assignments_

    _, indep = _dependent_param_exprs(M_type)
    hw_params = (indep..., :Keq, :E_total)
    mets = metabolites(M_type())

    Expr(:block,
        _destructuring_expr(hw_params, :params),
        _destructuring_expr(mets, :concs),
        r_assignments...,
        t_assignments...,
        :(return $rate_expr))
end
```

- [ ] **Step 3: Skip T-assignments in `rate_equation_string` when t_state_dead**

Locate `rate_equation_string(::AllostericEnzymeMechanism, ::ReducedMode)` (around line 1526). Find the local that captures `r_assignments, t_assignments = _build_dep_assignments(...)` and apply the same `_t_state_dead` skip:

Old code:

```julia
function rate_equation_string(
    ::AllostericEnzymeMechanism{CM,CS,RS}, ::ReducedMode,
) where {CM,CS,RS}
    M = AllostericEnzymeMechanism{CM,CS,RS}
    _, indep = _dependent_param_exprs(M)
    hw_params = (indep..., :Keq, :E_total)
    mets = metabolites(M())

    r_assignments, t_assignments = _build_dep_assignments(M, K -> :(1 / $K))
```

New code:

```julia
function rate_equation_string(
    ::AllostericEnzymeMechanism{CM,CS,RS}, ::ReducedMode,
) where {CM,CS,RS}
    M = AllostericEnzymeMechanism{CM,CS,RS}
    _, indep = _dependent_param_exprs(M)
    hw_params = (indep..., :Keq, :E_total)
    mets = metabolites(M())

    r_assignments, t_assignments_ = _build_dep_assignments(M, K -> :(1 / $K))
    t_assignments = _t_state_dead(M()) ? Expr[] : t_assignments_
```

- [ ] **Step 4: Run the full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass. PFK and HK rate equation strings will no longer contain `+ L * 0 * (...)` or the dead `K1_T = K1` mirrors.

- [ ] **Step 5: Verify the cleaner output for HK**

```bash
julia --project -e '
using EnzymeRates
include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
println(rate_equation_string(hk_mechanism))
' | head -3
```

Expected output starts with `(; K1, K4, k6f, K7, K10, K12, K_Pi_reg1, K_G6P_T_reg1, L, Keq, E_total) = params` and the `v = ...` line should NOT contain `L * 0 *`.

- [ ] **Step 6: Commit**

```bash
git add src/rate_eq_derivation.jl
git commit -m "Drop L*N_T numerator term and dead T-assignments when t_state_dead

When _t_state_dead = true, N_T is forced to 0; the L*num_T branch
expanded to L * 0 * Q_T_polynomial * Q_reg_T^... — wasteful at the
symbolic level (carried in the @generated rate equation body) and
ugly in rate_equation_string output. Skip the branch entirely.
Also skip t_assignments (T-state Haldane + :EqualRT catalytic mirrors
K_T = K) since they're only referenced from the now-elided branch."
```

---

### Task 6: Drop `:OnlyT` from catalytic-group enumeration in `_expand_to_allosteric`

`_expand_to_allosteric` enumerates `(:OnlyR, :OnlyT, :EqualRT)` for each catalytic kinetic group. With the upcoming constructor rule (Task 8), `:OnlyT` catalytic groups will error. The enumeration should not produce them.

**Files:**
- Modify: `src/mechanism_enumeration.jl` (around line 1811)

- [ ] **Step 1: Replace the tag tuple in `_expand_to_allosteric`**

Locate `_expand_to_allosteric` in `src/mechanism_enumeration.jl`. Find the inner loop:

```julia
        for tag in (:OnlyR, :OnlyT, :EqualRT)
            tag == :OnlyT && iso_only && continue
            new_tags = copy(base_tags)
```

Replace with:

```julia
        for tag in (:OnlyR, :EqualRT)
            new_tags = copy(base_tags)
```

The `iso_only` filter was about iso-only groups not being able to be `:OnlyT` (the relabel-ambiguity case). With `:OnlyT` removed entirely, that filter is moot.

Also update the docstring 4 lines above (around line 1788) — change `{:OnlyR, :OnlyT, :EqualRT}` to `{:OnlyR, :EqualRT}` and remove the "Iso-only groups skip `:OnlyT`" sentence.

- [ ] **Step 2: Run the mechanism enumeration tests**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep -E "test_mechanism_enumeration|Mechanism Enumeration|Test Summary"
```

Expected: all enumeration tests pass. Some test counts may decrease (since `:OnlyT` catalytic variants are no longer produced).

- [ ] **Step 3: Commit**

```bash
git add src/mechanism_enumeration.jl
git commit -m "Drop :OnlyT from catalytic-group enumeration in _expand_to_allosteric

Per the R-state-active convention, :OnlyT catalytic groups will error
on construction (Task 8). The enumeration should never produce them."
```

---

### Task 7: Switch `AllostericEnzymeMechanism` to dense `cat_allo_states` / `reg_allo_states`

The biggest task. Renames the type-parameter slots, requires every catalytic kinetic group and every regulatory ligand to have an explicit allo-state entry, and drops the "tag" terminology in favor of "allo state" throughout. Validation moves from the constructor's `tag in (:OnlyR, :OnlyT, :EqualRT, :NonequalRT)` checks to a dense-coverage check.

**Files:**
- Modify: `src/types.jl` — struct fields, accessors, constructor.
- Modify: `src/dsl.jl` — macro emits dense state tuples.
- Modify: `src/mechanism_enumeration.jl` — `AllostericMechanismSpec` field renames + dense storage.
- Modify: `src/rate_eq_derivation.jl` — call sites of `group_tag` and `regulatory_ligand_tag`.
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` — migrate 6 mechanisms.
- Modify: `test/test_dsl.jl` — migrate 4 mechanisms.
- Modify: `test/test_types.jl` — migrate 2 mechanisms.
- Modify: `.claude/CLAUDE.md` — terminology update.

**Naming summary:**
- Field: `group_tags::Tuple{Pair{Int,Symbol}...}` → `cat_allo_states::Tuple{Symbol...}` (dense, indexed by kinetic group number; entry `i` is the state for group `i`).
- Field: regulator ligand tags inside `RS` entries — change shape from `Tuple{Pair{Symbol,Symbol}...}` (sparse) to `Tuple{Symbol...}` (dense, parallel to `ligands` tuple).
- Accessor: `group_tag(m, g)` → `cat_allo_state(m, g)`.
- Accessor: `regulatory_ligand_tag(m, site_idx, lig)` → `reg_allo_state(m, site_idx, lig)`.

- [ ] **Step 1: Define new type-parameter shape and accessors in `src/types.jl`**

Locate the `AllostericEnzymeMechanism` struct (around line 240–280) and its docstring. Update both:

The `CatSites` parameter changes shape: `(multiplicity::Int, group_tags::Tuple{Pair{Int,Symbol}...})` becomes `(multiplicity::Int, cat_allo_states::Tuple{Symbol...})` where the tuple's length must equal the number of catalytic kinetic groups in the catalytic mechanism.

The `RS` regulator-site entries change from `(ligands, multiplicity, ligand_tags::Tuple{Pair{Symbol,Symbol}...})` to `(ligands, multiplicity, reg_allo_states::Tuple{Symbol...})` where the states tuple is parallel to `ligands` (same length, same order).

Replace the existing `AllostericEnzymeMechanism` docstring intro paragraph and field descriptions:

```julia
"""
    AllostericEnzymeMechanism{CatalyticMech, CatSites, RegSites}

Multi-subunit MWC allosteric enzyme. `CatalyticMech` is an
`EnzymeMechanism` describing one catalytic subunit's cycle.

# Type parameters
- `CatSites`: `(multiplicity::Int, cat_allo_states::Tuple{Symbol...})`.
  `cat_allo_states[g]` is the allosteric state of catalytic kinetic
  group `g` (1-indexed, dense — every group must have an entry).
  Allowed values: `:EqualRT`, `:NonequalRT`, `:OnlyR`. `:OnlyT`
  catalytic groups error during construction (R-state-active
  convention).
- `RegSites`: tuple of regulator-site entries
  `(ligands::Tuple{Symbol...}, multiplicity::Int,
   reg_allo_states::Tuple{Symbol...})` where `reg_allo_states` is
  parallel to `ligands`. Allowed values: all four states
  (`:EqualRT`, `:NonequalRT`, `:OnlyR`, `:OnlyT`).

Constructor validates:
- Catalytic state count matches kinetic-group count.
- Regulator state tuple length matches ligand tuple length at each site.
- No catalytic group has `:OnlyT` state.
- At least one ligand at each reg site is non-`:EqualRT` (single-
  ligand `:EqualRT` site cancels identically and conveys no allosteric
  effect; two-ligand all-`:EqualRT` site likewise).
"""
struct AllostericEnzymeMechanism{
    CatalyticMech<:EnzymeMechanism, CatSites, RegSites,
} end
```

Replace the `group_tag` accessor (around line 583):

```julia
"""Return the allosteric state of catalytic kinetic group `g`."""
function cat_allo_state(::AllostericEnzymeMechanism{CM, CS, RS}, g::Int) where {CM, CS, RS}
    _, states = CS
    return states[g]
end
```

Replace the `regulatory_ligand_tag` accessor (around line 660):

```julia
"""Return the allosteric state of regulator ligand `lig` at site `site_idx`."""
function reg_allo_state(
    ::AllostericEnzymeMechanism{CM, CS, RS}, site_idx::Int, lig::Symbol,
) where {CM, CS, RS}
    ligands, _, states = RS[site_idx]
    idx = findfirst(==(lig), ligands)
    idx === nothing && error("Ligand $lig not at regulatory site $site_idx")
    return states[idx]
end
```

Update the `regulatory_site_ligands` accessor — it doesn't change behavior but the underlying tuple shape is now confirmed parallel.

- [ ] **Step 2: Update the 3-arg constructor in `src/types.jl`**

Locate the 3-arg `AllostericEnzymeMechanism(cm, cat_sites, reg_sites)` constructor (around line 270–315). Rewrite the validation:

```julia
function AllostericEnzymeMechanism(
    cm::EnzymeMechanism, cat_sites::Tuple, reg_sites::Tuple,
)
    multiplicity, cat_allo_states = cat_sites
    multiplicity isa Int && multiplicity ≥ 1 ||
        error("Catalytic multiplicity must be a positive Int, got $multiplicity")

    n_groups = length(unique(kinetic_group(cm, i) for i in 1:n_steps(cm)))
    length(cat_allo_states) == n_groups ||
        error("cat_allo_states length $(length(cat_allo_states)) does not " *
              "match catalytic kinetic-group count $n_groups")
    for (g, st) in enumerate(cat_allo_states)
        st in (:OnlyR, :EqualRT, :NonequalRT) ||
            (st == :OnlyT &&
                error("Catalytic kinetic group $g has state :OnlyT; the " *
                      "R-state is the active state by convention. Relabel " *
                      "your mechanism so the active state is R (use :OnlyR " *
                      "instead).")) ||
            error("Catalytic kinetic group $g has unknown allo state $st; " *
                  "must be one of (:OnlyR, :EqualRT, :NonequalRT)")
    end

    # Reject mechanism if any catalytic group is :OnlyR AND any other is :OnlyT
    # — already excluded above, but explicit for clarity if rule changes later.

    for (i, entry) in enumerate(reg_sites)
        ligands, n_reg, reg_allo_states = entry
        ligands isa Tuple && all(l isa Symbol for l in ligands) ||
            error("Reg site $i: ligands must be a Tuple of Symbol")
        n_reg isa Int && n_reg ≥ 1 ||
            error("Reg site $i: multiplicity must be a positive Int")
        length(reg_allo_states) == length(ligands) ||
            error("Reg site $i: reg_allo_states length $(length(reg_allo_states)) " *
                  "does not match ligand count $(length(ligands))")
        for (k, st) in enumerate(reg_allo_states)
            st in (:OnlyR, :OnlyT, :EqualRT, :NonequalRT) ||
                error("Reg site $i, ligand $(ligands[k]): unknown allo state $st")
        end
        # All-:EqualRT site cancels identically — error
        all(st == :EqualRT for st in reg_allo_states) &&
            error("Reg site $i: all ligands are :EqualRT, which produces " *
                  "Q_reg_R == Q_reg_T — no allosteric effect. At least one " *
                  "ligand must be :OnlyR, :OnlyT, or :NonequalRT.")
    end

    AllostericEnzymeMechanism{typeof(cm), cat_sites, reg_sites}()
end
```

- [ ] **Step 3: Update DSL emission in `src/dsl.jl`**

Locate the `@allosteric_mechanism` macro implementation. The macro currently:
- Parses step groups, collecting `(group_id => tag)` pairs.
- Emits sparse-storage `cat_sites = (multiplicity, group_tags::Tuple{Pair...})`.
- Parses regulators into `(name => tag)` pairs.
- Emits sparse-storage `reg_sites = (ligands, multiplicity, ligand_tags::Tuple{Pair...})`.

Modify the emission to dense:
- After collecting `(g => state)` pairs from steps, look up each group `1:n_groups` in the dict and emit `Tuple{Symbol}` of states in group order. Error if any group is missing a state.
- After collecting `(lig => state)` pairs from regulators, emit a `Tuple{Symbol}` parallel to the `ligands` tuple. Error if any ligand is missing a state.

Pseudo-code (the actual symbol manipulation lives in `dsl.jl`'s macro logic):

```julia
# After parsing all step-group tags into a Dict{Int, Symbol}:
n_groups = length(group_dict)  # all kinetic groups
state_tuple = ntuple(g -> get(group_states, g) do
    error("@allosteric_mechanism: catalytic kinetic group $g has no allo " *
          "state annotation. Every step (or step group) must be tagged " *
          "with `:: <:OnlyR|:EqualRT|:NonequalRT>`.")
end, n_groups)
cat_sites_expr = :(($multiplicity, $state_tuple))

# For each reg site, after parsing ligand tags:
for site in reg_sites
    state_tuple_lig = ntuple(k -> get(site.lig_states, site.ligands[k]) do
        error("@allosteric_mechanism: regulator ligand $(site.ligands[k]) at " *
              "site $(site.idx) has no allo state annotation. Every regulator " *
              "must be `<lig>::<:OnlyR|:OnlyT|:EqualRT|:NonequalRT>`.")
    end, length(site.ligands))
    site_entry_expr = :((($(site.ligands...),), $(site.mult), $state_tuple_lig))
end
```

Find the existing implementation around `_parse_steps_block_with_groups` and `parse_allosteric_mechanism` calls in `dsl.jl` (search for `:NonequalRT` defaults and `Pair{Int, Symbol}`). Replace defaults with explicit-required errors and switch tuple shapes to flat Symbol tuples.

- [ ] **Step 4: Update `AllostericMechanismSpec` in `src/mechanism_enumeration.jl`**

The struct has fields `group_tags::Dict{Int, Symbol}` and `reg_ligand_tags::Dict{Symbol, Symbol}`. These are sparse internally for the enumeration, but the `AllostericEnzymeMechanism(spec)` constructor (around line 1258) converts to type parameters. Update that constructor to emit dense tuples per the new type-parameter shape.

Find the constructor (search for `AllostericMechanismSpec` and `AllostericEnzymeMechanism(spec)`):

```julia
function AllostericEnzymeMechanism(spec::AllostericMechanismSpec)
    cm = EnzymeMechanism(spec.base)

    n_groups = length(unique(s.kinetic_group for s in spec.base.steps))
    cat_states = ntuple(g -> get(spec.group_tags, g, :NonequalRT), n_groups)

    reg_sites = ntuple(length(spec.allosteric_reg_sites)) do i
        ligs = Tuple(spec.allosteric_reg_sites[i])
        mult = spec.allosteric_multiplicities[i]
        lig_states = ntuple(k -> get(spec.reg_ligand_tags, ligs[k], :NonequalRT),
                            length(ligs))
        (ligs, mult, lig_states)
    end

    AllostericEnzymeMechanism(cm, (spec.catalytic_n, cat_states), reg_sites)
end
```

The `Dict`-based sparse storage stays internal to the enumeration logic — that's fine, only the type-parameter conversion changes.

- [ ] **Step 5: Update call sites of `group_tag` and `regulatory_ligand_tag` in `src/rate_eq_derivation.jl`**

Run:

```bash
grep -n "group_tag\b\|regulatory_ligand_tag" src/rate_eq_derivation.jl
```

Replace each call:
- `group_tag(m, g)` → `cat_allo_state(m, g)`
- `regulatory_ligand_tag(m, site_idx, lig)` → `reg_allo_state(m, site_idx, lig)`

Sites to update include `_T_rename`, `_onlyT_syms`, `_onlyR_syms`, `_t_state_dead` (already updated in Task 4 — re-check), `_reg_site_expr`, `_build_dep_assignments`, `_kcat_forward`, `_allosteric_num_den_exprs`, `rate_equation_string`. Verify post-edit:

```bash
grep -c "group_tag\b\|regulatory_ligand_tag" src/rate_eq_derivation.jl
```

Expected: 0.

- [ ] **Step 6: Migrate the 6 `@allosteric_mechanism` definitions in `test/mechanism_definitions_for_test_enzyme_derivation.jl`**

Search for all `@allosteric_mechanism begin` blocks. For each, ensure every step (or step group) has an explicit `:: <state>` annotation, and every entry in `allosteric_regulators:` has an explicit `<lig>::<state>` annotation. The PFK and HK definitions already do this. Check the others (any single-substrate edge cases, etc.) and add explicit `:: NonequalRT` where the previous default was relied on.

Verify with `grep -n ":: " test/mechanism_definitions_for_test_enzyme_derivation.jl` — every step in every `@allosteric_mechanism` block should be annotated.

- [ ] **Step 7: Migrate `@allosteric_mechanism` blocks in `test/test_dsl.jl` and `test/test_types.jl`**

Same as Step 6. Search:

```bash
grep -n "@allosteric_mechanism" test/test_dsl.jl test/test_types.jl
```

For each occurrence, audit the steps and regulators, add explicit allo-state annotations.

- [ ] **Step 8: Update CLAUDE.md terminology**

In `.claude/CLAUDE.md`, replace `tag` with `allo state` (or `allosteric state`) where it refers to the per-group / per-ligand R/T-state classifier. Specifically:
- "Allosteric tag taxonomy" → "Allosteric state taxonomy"
- `group_tags` → `cat_allo_states`
- `reg_ligand_tags` → `reg_allo_states`
- `group_tag()` → `cat_allo_state()`
- `regulatory_ligand_tag()` → `reg_allo_state()`

Search `grep -n "tag\|Tag" .claude/CLAUDE.md` and update every occurrence in the AllostericEnzymeMechanism context.

- [ ] **Step 9: Run the full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass. Test count may shift slightly due to dropped `:OnlyT` catalytic enumeration variants from Task 6.

- [ ] **Step 10: Commit**

```bash
git add src/types.jl src/dsl.jl src/mechanism_enumeration.jl src/rate_eq_derivation.jl
git add test/mechanism_definitions_for_test_enzyme_derivation.jl test/test_dsl.jl test/test_types.jl
git add .claude/CLAUDE.md
git commit -m "Switch AllostericEnzymeMechanism to dense cat_allo_states / reg_allo_states

Renames type parameters and accessors:
  group_tags         → cat_allo_states     (Tuple{Symbol...} per kinetic group)
  reg_ligand_tags    → reg_allo_states     (Tuple{Symbol...} parallel to ligands)
  group_tag()        → cat_allo_state()
  regulatory_ligand_tag() → reg_allo_state()

Drops sparse-storage default-:NonequalRT shortcut. Every catalytic
kinetic group must have an explicit allo-state entry; same for every
regulator ligand. Constructor errors on :OnlyT catalytic groups
(R-state-active convention). Migrates all 12 mechanism definitions
across tests + CLAUDE.md terminology."
```

---

### Task 8: DSL: error on bare step or bare regulator (no `::AlloState`)

The previous DSL allowed steps and regulators without explicit allo-state annotation, defaulting to `:NonequalRT`. After Task 7, the runtime constructor errors on missing entries, but a friendlier error at parse-time is preferable. Add explicit checks in `@allosteric_mechanism` macro so the error message points at the offending step/regulator line.

**Files:**
- Modify: `src/dsl.jl`
- Add tests: `test/test_rate_eq_derivation.jl` (the renamed file from Task 2) — in the "Allosteric edge cases" testset

- [ ] **Step 1: Reject step without `::AlloState` in `_parse_steps_block_with_groups`**

In `src/dsl.jl`, locate `_parse_steps_block_with_groups`. The function is called with `allow_tag::Bool=false` for plain mechanisms and `allow_tag=true` for allosteric. Currently `allow_tag=true` accepts steps both with AND without tags (defaulting to NonequalRT). Change to: when `allow_tag=true`, every step must have a tag — bare steps error.

Find the branch around line 337 (`# Parenthesized-group-without-tag (plain)`) and the branch around line 343 (`# Single step (with or without tag)`). When `allow_tag=true`, both must error if no tag is found.

Pseudo:

```julia
elseif arg.head == :tuple
    if allow_tag
        error("@allosteric_mechanism: step group $arg has no allo state " *
              "annotation. Add `:: <:OnlyR|:EqualRT|:NonequalRT>` to the " *
              "parenthesized step group.")
    end
    # else original plain-mechanism path
elseif <single-step pattern>
    tag = _peel_step_tag!(arg)
    if tag === nothing && allow_tag
        error("@allosteric_mechanism: step `$original` has no allo state " *
              "annotation. Add `:: <:OnlyR|:EqualRT|:NonequalRT>` after the " *
              "step.")
    end
    # ...
```

- [ ] **Step 2: Reject regulator without `::AlloState`**

In `src/dsl.jl`, locate the regulator parsing for `allosteric_regulators:` (and `dead_end_inhibitors:` / `regulators:` are for non-allosteric `@enzyme_reaction`, leave untouched). For each entry in the allosteric-regulators list, require an `Expr(:(::), name, state)` form; bare `name` errors.

Pseudo:

```julia
for entry in allosteric_regulators_list
    if entry isa Symbol
        error("@allosteric_mechanism: regulator `$entry` has no allo state " *
              "annotation. Use `$entry::<:OnlyR|:OnlyT|:EqualRT|:NonequalRT>`.")
    end
    # ... handle Expr(:(::), name, state) ...
end
```

- [ ] **Step 3: Add three error tests in "Allosteric edge cases"**

Open `test/test_rate_eq_derivation.jl` and locate `@testset "Allosteric edge cases" begin`. Add these three tests at the end of the testset:

```julia
    # Bare step (no :: AlloState) → DSL error
    @test_throws Exception eval(:(@allosteric_mechanism begin
        substrates: S
        products:   P
        site(:catalytic, 2): begin
            steps: begin
                [E, S] ⇌ [ES]                # missing :: AlloState
                [ES] <--> [EP]  :: EqualRT
                [EP] ⇌ [E, P]   :: EqualRT
            end
        end
    end))

    # Bare regulator (no ::AlloState) → DSL error
    @test_throws Exception eval(:(@allosteric_mechanism begin
        substrates: S
        products:   P
        allosteric_regulators: I, J::OnlyT
        site(:catalytic, 2): begin
            steps: begin
                [E, S] ⇌ [ES]   :: EqualRT
                [ES] <--> [EP]  :: EqualRT
                [EP] ⇌ [E, P]   :: EqualRT
            end
        end
    end))

    # Bare step group (parenthesized, no :: AlloState) → DSL error
    @test_throws Exception eval(:(@allosteric_mechanism begin
        substrates: S, A
        products:   P
        site(:catalytic, 2): begin
            steps: begin
                ([E, S] ⇌ [ES], [E_A, S] ⇌ [ES_A])    # missing :: AlloState
                [ES_A] <--> [EP]   :: EqualRT
                [EP] ⇌ [E, P]      :: EqualRT
            end
        end
    end))
```

- [ ] **Step 4: Run the full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass; the three new `@test_throws` tests are now in the count.

- [ ] **Step 5: Commit**

```bash
git add src/dsl.jl test/test_rate_eq_derivation.jl
git commit -m "DSL: error on bare step or bare regulator (no allo state)

Allosteric mechanisms must declare an explicit allosteric state for
every catalytic step (or step group) and every regulator ligand.
Defaulting to :NonequalRT silently is no longer allowed — the error
now points at the offending line in the macro body."
```

---

### Task 9: Constructor error on `:OnlyT` catalytic group + tests

The constructor in Task 7 already rejects `:OnlyT` catalytic states. Add the corresponding `@test_throws` tests, and convert the existing `onlyT_sub` test (which previously verified rate-vanishing) into a constructor-error test.

**Files:**
- Modify: `test/test_rate_eq_derivation.jl` — convert/add tests in "Allosteric edge cases"

- [ ] **Step 1: Find and remove the existing `onlyT_sub` test**

In `test/test_rate_eq_derivation.jl`, locate the block starting with `# OnlyT substrate: S binds only in T-state, ...` and the `onlyT_sub = @allosteric_mechanism ...` definition. Delete the entire block including the `rate_strong = ...`, `rate_weak = ...`, and `@test rate_strong > 1.0` assertions through `@test rate_weak / rate_strong < 1e-5`.

- [ ] **Step 2: Add `:OnlyT` substrate constructor-error test**

Insert at the location where `onlyT_sub` was:

```julia
    # :OnlyT on a substrate-binding catalytic group → constructor error
    # (R-state convention: relabel so the active state is R).
    @test_throws Exception eval(:(@allosteric_mechanism begin
        substrates: S
        products:   P
        site(:catalytic, 2): begin
            steps: begin
                [E, S] ⇌ [ES]    :: OnlyT
                [ES] <--> [EP]   :: EqualRT
                [EP] ⇌ [E, P]    :: EqualRT
            end
        end
    end))
```

- [ ] **Step 3: Add `:OnlyT` product constructor-error test**

Add immediately after Step 2:

```julia
    # :OnlyT on a product-binding catalytic group → constructor error
    @test_throws Exception eval(:(@allosteric_mechanism begin
        substrates: S
        products:   P
        site(:catalytic, 2): begin
            steps: begin
                [E, S] ⇌ [ES]    :: EqualRT
                [ES] <--> [EP]   :: EqualRT
                [EP] ⇌ [E, P]    :: OnlyT
            end
        end
    end))
```

- [ ] **Step 4: Add `:OnlyT` catalysis (SS) constructor-error test**

```julia
    # :OnlyT on the catalysis SS step → constructor error
    @test_throws Exception eval(:(@allosteric_mechanism begin
        substrates: S
        products:   P
        site(:catalytic, 2): begin
            steps: begin
                [E, S] ⇌ [ES]    :: EqualRT
                [ES] <--> [EP]   :: OnlyT
                [EP] ⇌ [E, P]    :: EqualRT
            end
        end
    end))
```

- [ ] **Step 5: Add mixed `:OnlyR sub + :OnlyT prod` test (already errors via :OnlyT alone, but explicit)**

```julia
    # :OnlyR substrate + :OnlyT product → constructor error
    # (subsumed by single-:OnlyT error, but explicit confirms the rule.)
    @test_throws Exception eval(:(@allosteric_mechanism begin
        substrates: S
        products:   P
        site(:catalytic, 2): begin
            steps: begin
                [E, S] ⇌ [ES]    :: OnlyR
                [ES] <--> [EP]   :: EqualRT
                [EP] ⇌ [E, P]    :: OnlyT
            end
        end
    end))
```

- [ ] **Step 6: Run the full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass; 4 new `@test_throws` tests are in the count.

- [ ] **Step 7: Commit**

```bash
git add test/test_rate_eq_derivation.jl
git commit -m "Constructor error on :OnlyT catalytic group + tests

Convert the previous onlyT_sub edge case (which built a T-state-active
mechanism) into a constructor-error test. Add three more cases:
:OnlyT substrate, :OnlyT product, :OnlyT catalysis. Plus the
:OnlyR sub + :OnlyT prod combo for explicit coverage of the rule."
```

---

### Task 10: Remove standalone PFK and HK rate-equation testsets

`test_analytical_rate(spec)` already covers the rate-equation match via random params over many trials. The standalone testsets only add monotonicity sweeps, which are nice but not essential and don't fit the spec battery cleanly.

**Files:**
- Modify: `test/test_rate_eq_derivation.jl`

- [ ] **Step 1: Delete `@testset "PFK rate equation matches analytical form"`**

Locate the testset (search for `PFK rate equation matches analytical`). Delete from the `@testset` opening through the matching `end`. Approximately lines 922–972.

- [ ] **Step 2: Delete `@testset "HK rate equation matches analytical form"`**

Locate (search for `HK rate equation matches analytical`). Delete from `@testset` through the matching `end`. Approximately lines 974–1017.

- [ ] **Step 3: Run the full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass. Test count drops by ~12 (the assertions inside the deleted testsets).

- [ ] **Step 4: Commit**

```bash
git add test/test_rate_eq_derivation.jl
git commit -m "Remove standalone PFK and HK rate-equation testsets

The exact-conc/exact-param rate match is already covered by
test_analytical_rate(spec) with random params over many trials, and
analytical_kcat_fn covers saturation behavior. The biology-flavored
monotonicity sweeps in the standalone tests are nice but not load-
bearing — they don't fit the spec battery cleanly."
```

---

### Task 11: Add `m_PK` (pyruvate kinase) to MECHANISM_TEST_SPECS

Pyruvate kinase: PEP + ADP → Pyruvate + ATP. Tests `:NonequalRT` substrate (PEP) and mismatched-multiplicity reg sites (ATP at mult=2 :OnlyT, F16BP at mult=4 :OnlyR).

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`

- [ ] **Step 1: Add the PK mechanism + analytical formula to `build_mechanism_test_specs`**

Insert immediately after the HK block (around line 2095) inside `build_mechanism_test_specs()`:

```julia
    # ── Pyruvate kinase (PK) hand-verified mechanism ──────────────────────────
    # Reaction: PEP + ADP ⇌ Pyruvate + ATP, 4 catalytic subunits.
    # PEP binding is :NonequalRT (independent K_R and K_T) so the T-state
    # cycle is alive. Reg sites have MISMATCHED multiplicities:
    #   ATP::OnlyT at mult 2
    #   F16BP::OnlyR at mult 4 (matches catalytic mult)
    # This exercises the symmetric all-reg-sites contribution to both
    # numerator and denominator (after dropping the n_reg == CatN filter).
    let
        m = @allosteric_mechanism begin
            substrates: PEP, ADP
            products:   Pyruvate, ATP
            allosteric_regulators: ATP::OnlyT, F16BP::OnlyR

            site(:catalytic, 4): begin
                steps: begin
                    ([E, PEP] ⇌ [E_PEP],
                     [E_ADP, PEP] ⇌ [E_PEP_ADP])              :: NonequalRT
                    ([E, ADP] ⇌ [E_ADP],
                     [E_PEP, ADP] ⇌ [E_PEP_ADP])              :: EqualRT
                    [E_PEP_ADP] <--> [E_Pyr_ATP]               :: EqualRT
                    ([E_Pyr_ATP] ⇌ [E_ATP, Pyruvate],
                     [E_Pyr] ⇌ [E, Pyruvate])                  :: EqualRT
                    ([E_Pyr_ATP] ⇌ [E_Pyr, ATP],
                     [E_ATP] ⇌ [E, ATP])                       :: EqualRT
                end
            end

            site(:regulatory, 2): begin
                ligands: ATP::OnlyT
            end
            site(:regulatory, 4): begin
                ligands: F16BP::OnlyR
            end
        end

        # Param mapping (kinetic-group representative-step convention):
        #   K1, K1_T : PEP binding (group 1, NonequalRT)
        #   K3       : ADP binding (group 2, EqualRT)
        #   k5f      : catalysis SS (group 3, EqualRT)
        #   K6       : Pyruvate release (group 4, EqualRT)
        #   K8       : ATP release (group 5, EqualRT)
        function pk_rate_analytical(params, concs)
            (; K1, K1_T, K3, k5f, K6, K8,
               K_ATP_T_reg1, K_F16BP_reg2,
               L, Keq, Et) = params
            (; PEP, ADP, Pyruvate, ATP, F16BP) = concs

            k5r   = k5f * K6 * K8 / (Keq * K1   * K3)
            k5r_T = k5f * K6 * K8 / (Keq * K1_T * K3)

            Q_cat_R = 1 + PEP/K1 + ADP/K3 + PEP*ADP/(K1*K3) +
                      Pyruvate/K6 + ATP/K8 + Pyruvate*ATP/(K6*K8)
            Q_cat_T = 1 + PEP/K1_T + ADP/K3 + PEP*ADP/(K1_T*K3) +
                      Pyruvate/K6 + ATP/K8 + Pyruvate*ATP/(K6*K8)

            N_R = k5f * PEP * ADP / (K1   * K3) - k5r   * Pyruvate * ATP / (K6 * K8)
            N_T = k5f * PEP * ADP / (K1_T * K3) - k5r_T * Pyruvate * ATP / (K6 * K8)

            Q_reg1_R = 1                                     # ATP::OnlyT, no R term
            Q_reg1_T = 1 + ATP / K_ATP_T_reg1
            Q_reg2_R = 1 + F16BP / K_F16BP_reg2               # F16BP::OnlyR, no T term
            Q_reg2_T = 1

            num_R = N_R * Q_cat_R^3 * Q_reg1_R^2 * Q_reg2_R^4
            num_T = N_T * Q_cat_T^3 * Q_reg1_T^2 * Q_reg2_T^4
            den_R = Q_cat_R^4 * Q_reg1_R^2 * Q_reg2_R^4
            den_T = Q_cat_T^4 * Q_reg1_T^2 * Q_reg2_T^4

            return Et * 4.0 * (num_R + L * num_T) / (den_R + L * den_T)
        end

        push!(specs, MechanismTestSpec(
            name="PK",
            mechanism=m,
            metabolite_names=[:PEP, :ADP, :Pyruvate, :ATP, :F16BP],
            expected_n_states=7,           # E, E_PEP, E_ADP, E_PEP_ADP, E_Pyr_ATP, E_Pyr, E_ATP
            expected_n_steps=9,
            expected_n_metabolites=5,
            expected_n_haldane=2,           # k5r derived in R-state, k5r_T in T-state
            expected_n_wegscheider=0,
            # 4 EqualRT mirrors (K3, K6, K8 each + k5f) plus 1 NonequalRT
            # T-state Haldane dependent on K1_T → expected_n_haldane covers k5r and k5r_T.
            # Independent params: K1, K1_T, K3, k5f, K6, K8, K_ATP_T_reg1, K_F16BP_reg2, L = 9
            expected_n_independent_params=9,
            expected_identifiability_deficit=0,  # placeholder — verify and update
            expected_is_identifiable=true,
            run_ode_test=false,
            analytical_rate_fn=pk_rate_analytical,
            # kcat: at saturating substrates and zero products, A_R/B_R = A_T/B_T = k5f
            # because catalysis is :EqualRT, so all corner kcats reduce to catN · k5f.
            analytical_kcat_fn = p -> 4 * p.k5f,
            expected_factored_num=nothing,
            expected_factored_denom=nothing,
        ))
    end
```

The placeholder `expected_identifiability_deficit=0` will need correction (next step).

- [ ] **Step 2: Run the test suite once and capture the actual identifiability deficit and other counts**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep -E "PK|Test Summary|Evaluated"
```

The Constraints / Identifiability tests will fail on the placeholder values. Read the failure messages to extract:
- `expected_n_haldane` actual
- `expected_n_independent_params` actual
- `expected_identifiability_deficit` actual

Update the spec lines accordingly.

Also independently verify analytical_rate_fn matches `rate_equation`:

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
'
```

Expected: both numbers equal to ~10 decimal places.

- [ ] **Step 3: Run the full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass.

- [ ] **Step 4: Commit**

```bash
git add test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -m "Add pyruvate kinase (PK) mechanism to MECHANISM_TEST_SPECS

PEP::NonequalRT (substrate), ATP::OnlyT (mult 2 reg site),
F16BP::OnlyR (mult 4 reg site). Tests :NonequalRT substrate-binding
end-to-end and exercises the symmetric all-reg-sites numerator
contribution after dropping the n_reg == CatN filter (Task 3) —
both the mult-2 site (Q_reg1_T²) and mult-4 site (Q_reg2_R⁴) appear
in num and den. Hand-derived analytical rate matches @generated
rate_equation; analytical_kcat_fn = 4·k5f."
```

---

### Task 12: Add `m_all` to MECHANISM_TEST_SPECS

A simple BiBi mechanism testing `:NonequalRT` on substrate (S1) and product (P1), `:EqualRT` on the rest, plus a 2-ligand mixed-state reg site (R1::NonequalRT + R2::EqualRT).

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`

- [ ] **Step 1: Add the `m_all` mechanism + analytical formula**

Insert immediately after the PK block:

```julia
    # ── m_all: NonequalRT coverage on substrate + product ───────────────────
    # Two-substrate two-product reaction with explicit :NonequalRT on
    # S1 binding and P1 release, and a 2-ligand mixed-state reg site
    # (R1::NonequalRT + R2::EqualRT). Catalysis is :EqualRT to keep the
    # kcat formula simple (kcat = catN · k5f at every corner).
    let
        m = @allosteric_mechanism begin
            substrates: S1, S2
            products:   P1, P2
            allosteric_regulators: R1::NonequalRT, R2::EqualRT

            site(:catalytic, 2): begin
                steps: begin
                    ([E, S1] ⇌ [E_S1],
                     [E_S2, S1] ⇌ [E_S1_S2])      :: NonequalRT
                    ([E, S2] ⇌ [E_S2],
                     [E_S1, S2] ⇌ [E_S1_S2])      :: EqualRT
                    [E_S1_S2] <--> [E_P1_P2]      :: EqualRT
                    ([E_P1_P2] ⇌ [E_P2, P1],
                     [E_P1] ⇌ [E, P1])             :: NonequalRT
                    ([E_P1_P2] ⇌ [E_P1, P2],
                     [E_P2] ⇌ [E, P2])             :: EqualRT
                end
            end

            site(:regulatory, 2): begin
                ligands: R1::NonequalRT, R2::EqualRT
            end
        end

        # Param mapping:
        #   K1, K1_T : S1 binding (group 1, NonequalRT)
        #   K3       : S2 binding (group 2, EqualRT)
        #   k5f      : catalysis (group 3, EqualRT)
        #   K6, K6_T : P1 release (group 4, NonequalRT)
        #   K8       : P2 release (group 5, EqualRT)
        function m_all_rate_analytical(params, concs)
            (; K1, K1_T, K3, k5f, K6, K6_T, K8,
               K_R1_reg1, K_R1_T_reg1, K_R2_reg1,
               L, Keq, Et) = params
            (; S1, S2, P1, P2, R1, R2) = concs

            k5r   = k5f * K6   * K8 / (Keq * K1   * K3)
            k5r_T = k5f * K6_T * K8 / (Keq * K1_T * K3)

            Q_cat_R = 1 + S1/K1   + S2/K3 + S1*S2/(K1   * K3) +
                      P1/K6   + P2/K8 + P1*P2/(K6   * K8)
            Q_cat_T = 1 + S1/K1_T + S2/K3 + S1*S2/(K1_T * K3) +
                      P1/K6_T + P2/K8 + P1*P2/(K6_T * K8)

            N_R = k5f * S1 * S2 / (K1   * K3) - k5r   * P1 * P2 / (K6   * K8)
            N_T = k5f * S1 * S2 / (K1_T * K3) - k5r_T * P1 * P2 / (K6_T * K8)

            Q_reg1_R = 1 + R1/K_R1_reg1   + R2/K_R2_reg1
            Q_reg1_T = 1 + R1/K_R1_T_reg1 + R2/K_R2_reg1

            num = N_R * Q_cat_R   * Q_reg1_R^2 + L * N_T * Q_cat_T   * Q_reg1_T^2
            den = Q_cat_R^2 * Q_reg1_R^2       + L *      Q_cat_T^2 * Q_reg1_T^2

            return Et * 2.0 * num / den
        end

        push!(specs, MechanismTestSpec(
            name="m_all",
            mechanism=m,
            metabolite_names=[:S1, :S2, :P1, :P2, :R1, :R2],
            expected_n_states=7,             # E, E_S1, E_S2, E_S1_S2, E_P1_P2, E_P1, E_P2
            expected_n_steps=9,
            expected_n_metabolites=6,
            expected_n_haldane=2,             # k5r in R-state, k5r_T in T-state
            expected_n_wegscheider=0,
            expected_n_independent_params=10, # K1, K1_T, K3, k5f, K6, K6_T, K8, K_R1_reg1, K_R1_T_reg1, K_R2_reg1, L = 11 → minus dependents... verify
            expected_identifiability_deficit=0,  # placeholder
            expected_is_identifiable=true,
            run_ode_test=false,
            analytical_rate_fn=m_all_rate_analytical,
            analytical_kcat_fn = p -> 2 * p.k5f,  # cat is :EqualRT
            expected_factored_num=nothing,
            expected_factored_denom=nothing,
        ))
    end
```

- [ ] **Step 2: Run the suite, capture actual count values, update placeholders**

Same procedure as Task 11 Step 2. Update `expected_n_independent_params`, `expected_n_haldane`, `expected_identifiability_deficit` to match the actual measured values from the failed test output. Also independently verify `analytical_rate_fn` matches `rate_equation`.

- [ ] **Step 3: Run the full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass.

- [ ] **Step 4: Commit**

```bash
git add test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -m "Add m_all mechanism to MECHANISM_TEST_SPECS

Two-substrate two-product BiBi covering :NonequalRT on substrate (S1)
and product (P1), :EqualRT elsewhere. 2-ligand mixed-state reg site
(R1::NonequalRT + R2::EqualRT) — exercises the per-ligand state lookup
and validates that two-EqualRT-only sites would error (single :EqualRT
ligand here is balanced by R1::NonequalRT). kcat = catN · k5f."
```

---

### Task 13: Add `m_OnlyR_prod` to MECHANISM_TEST_SPECS

Single substrate, single product, catN=2, with product-binding `:OnlyR` (T-state cycle dead via product side, R-state catalyzes). Tests the broadened `_t_state_dead` from Task 4 in the product-binding case.

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`

- [ ] **Step 1: Add the `m_OnlyR_prod` mechanism + analytical formula**

Insert immediately after the m_all block:

```julia
    # ── m_OnlyR_prod: single product :OnlyR (T-state cycle dead) ────────────
    # Tests broadened _t_state_dead detection for product-binding :OnlyR
    # groups (Task 4). Cycle in T-state is broken at product release;
    # forward catalysis from E_S_T → E_P_T appears blocked because E_P_T
    # form is absent (K3 zeroed in T-state poly). With t_state_dead = true,
    # N_T is forced to 0 and the L*num_T branch drops.
    let
        m = @allosteric_mechanism begin
            substrates: S
            products:   P

            site(:catalytic, 2): begin
                steps: begin
                    [E, S] ⇌ [E_S]    :: EqualRT      # group 1, K1
                    [E_S] <--> [E_P]  :: EqualRT      # group 2, k2f catalysis
                    [E, P] ⇌ [E_P]    :: OnlyR        # group 3, K3 (P binding)
                end
            end
        end

        # Param mapping:
        #   K1   : S binding (group 1, EqualRT)
        #   k2f  : catalysis SS (group 2, EqualRT, k2r derived via Haldane)
        #   K3   : P release (group 3, OnlyR)
        function m_OnlyR_prod_rate_analytical(params, concs)
            (; K1, k2f, K3, L, Keq, Et) = params
            (; S, P) = concs

            k2r = k2f * K3 / (Keq * K1)

            Q_cat_R = 1 + S/K1 + P/K3
            Q_cat_T = 1 + S/K1                    # P/K3 monomial dropped (OnlyR group 3)

            N_R = k2f * S/K1 - k2r * P/K3
            # N_T = 0 forced (t_state_dead via group 3 :OnlyR)

            num = N_R * Q_cat_R                   # L*N_T*Q_cat_T term elided
            den = Q_cat_R^2 + L * Q_cat_T^2

            return Et * 2.0 * num / den
        end

        push!(specs, MechanismTestSpec(
            name="m_OnlyR_prod",
            mechanism=m,
            metabolite_names=[:S, :P],
            expected_n_states=3,                  # E, E_S, E_P
            expected_n_steps=3,
            expected_n_metabolites=2,
            expected_n_haldane=1,                 # k2r derived in R-state only (T-state Haldane dropped via t_state_dead)
            expected_n_wegscheider=0,
            expected_n_independent_params=4,      # K1, k2f, K3, L = 4
            expected_identifiability_deficit=0,   # placeholder
            expected_is_identifiable=true,
            run_ode_test=false,
            analytical_rate_fn=m_OnlyR_prod_rate_analytical,
            # kcat at saturating S, zero P:
            #   A_R = k2f/K1², B_R = 1/K1², B_T = (1/K1)² (T-state pattern same as R)
            #   kcat = catN · A_R / (B_R + L · B_T) = 2 · k2f / (1 + L)
            analytical_kcat_fn = p -> 2 * p.k2f / (1 + p.L),
            expected_factored_num=nothing,
            expected_factored_denom=nothing,
        ))
    end
```

- [ ] **Step 2: Run the suite, capture actual count values, update placeholders**

Same procedure as Task 11 Step 2. Update `expected_n_independent_params`, `expected_n_haldane`, `expected_identifiability_deficit` to match the actual measured values.

Also verify analytical formula numerically:

```bash
julia --project -e '
using EnzymeRates
using Random
include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
m = _spec_by_name("m_OnlyR_prod").mechanism
spec = _spec_by_name("m_OnlyR_prod")
rng = Random.MersenneTwister(42)
indep = (:K1, :k2f, :K3, :L)
keys_t = (indep..., :Keq, :E_total)
vals_t = Tuple(0.1 + 9.9 * rand(rng) for _ in keys_t)
params = NamedTuple{keys_t}(vals_t)
concs = (S=0.5, P=0.1)
p_an = merge(params, (Et=params.E_total,))
println("rate_equation = ", rate_equation(m, concs, params))
println("analytical    = ", spec.analytical_rate_fn(p_an, concs))
println("kcat (orig)   = ", EnzymeRates._kcat_forward(m, params))
println("kcat (analytical) = ", spec.analytical_kcat_fn(params))
'
```

Expected: rate_equation == analytical and kcat_orig == analytical_kcat to ~10 decimals.

- [ ] **Step 3: Run the full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass.

- [ ] **Step 4: Commit**

```bash
git add test/mechanism_definitions_for_test_enzyme_derivation.jl
git commit -m "Add m_OnlyR_prod mechanism to MECHANISM_TEST_SPECS

Single substrate, single product, catN=2, with product-binding :OnlyR
(group 3). Exercises the broadened _t_state_dead detection from Task 4
— previous filter limited :OnlyR detection to substrate-binding and
isomerization, missing this case. Analytical rate and kcat both
L-dependent: kcat = 2·k2f/(1+L) at saturating substrate.

This is the canonical 'product-binding only-R' test that motivated
broadening _t_state_dead beyond the substrate-only filter."
```

---

## Self-Review

Spec coverage check:
- ✅ Remove `*_broken` flags → Task 1
- ✅ Rename test file → Task 2
- ✅ Drop `n_reg == CatN` filter → Task 3
- ✅ Simplify `_t_state_dead` → Task 4
- ✅ Drop `L*N_T` term + dead t_assignments → Task 5
- ✅ Drop `:OnlyT` from enumeration → Task 6
- ✅ Dense `cat_allo_states` / `reg_allo_states` rename + drop "tag" terminology + require explicit state → Task 7
- ✅ DSL: error on bare step / bare regulator → Task 8
- ✅ Constructor: error on `:OnlyT` catalytic group → Task 7 (constructor) + Task 9 (tests)
- ✅ Remove PFK/HK standalone testsets → Task 10
- ✅ Add PK → Task 11
- ✅ Add m_all → Task 12
- ✅ Add m_OnlyR_prod → Task 13

Type-consistency check: `cat_allo_states` (plural, used everywhere), `reg_allo_states` (plural), `cat_allo_state(m, g)` (singular accessor), `reg_allo_state(m, site_idx, lig)` (singular accessor). Consistent across tasks.

No placeholders remain in the body (the `expected_identifiability_deficit=0` lines in Tasks 11–13 are explicit-placeholder lines that the implementer is instructed to replace with measured values via running the failing test).

Order constraints: Task 7 introduces the rename + dense storage, so Tasks 8–13 reference `cat_allo_state`/`reg_allo_state` and the new field semantics — those come after. Task 6 must come before Task 7 (the constructor in Task 7 will error on `:OnlyT` catalytic groups; the enumeration must stop producing them first to keep tests green).
