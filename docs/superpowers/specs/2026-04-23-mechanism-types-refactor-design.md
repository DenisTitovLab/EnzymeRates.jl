# Mechanism Types Refactor — Design

**Status:** Draft (revised after first reviewer round)
**Date:** 2026-04-23 (revised 2026-04-24)
**Scope:** `EnzymeMechanism`, `AllostericEnzymeMechanism`, `@enzyme_mechanism`, `@allosteric_mechanism`, and all code that reads their internals.

---

## 1. Primary Goal: Remove Code

**Less code == simpler code.** The refactor is a success only if a substantial fraction of the affected files is deleted. The current fragmented representations (7-field `CatSites` tuples accessed by magic indices `CS[3]`, `CS[5]`, `CS[6]`, `CS[7]`; 5-per-entry `RegSites`; `_apply_param_constraints` machinery for fully-general monomial constraints; per-species atom inference; parallel R-state/T-state derivation) encode information in multiple places and force conditional branches throughout. A unified internal representation + accessor-only interface should let us delete substantial code in:

- `src/types.jl` — magic-index accessors, atom inference, `_count_side`, `_infer_enzyme_atoms`, `RegulatorRole` hierarchy, `graph()`, `param_constraints` stub.
- `src/rate_eq_derivation.jl` — `_is_tr_equiv_catalytic_param`, `_is_r_only_catalytic_param`, parallel-path R/T control flow in `_dependent_param_exprs` and `_build_allosteric_rate_body`, the duplicated allosteric `_kcat_forward`.
- `src/mechanism_enumeration.jl` — `_valid_allosteric_differentiations` K-type/V-type branches, `_rewrap_allosteric`, `_tr_equiv_met_delta`, `_is_mirror_of`, `_constrained_step_indices`, the `::AllostericMechanismSpec` variants of every expansion move.
- `src/sym_poly_for_rate_eq_derivation.jl` — `_rs_tr_equiv` / `_rs_r_only` / `_rs_t_only` magic-index helpers, `_count_allosteric_rate_monomials`, `_apply_param_constraints` methods on POLY / FactoredSigma / FactoredPoly / DenomTerm.
- `src/dsl.jl` — `_walk_rhs!`, `_parse_constraint_rhs`, `_push_constraint!`, the heuristic `_is_allosteric_label`, the dual `_parse_enzyme_mechanism` / `_parse_allosteric_mechanism` flavors.

**If the implementation PR doesn't show large net code reduction across these files, the design is still carrying redundancy somewhere and should be revisited.**

---

## 2. Motivation

Beyond code reduction, this refactor addresses four concrete capability gaps:

1. **Rigorous per-step TR-binding-mode taxonomy.** Each kinetic group in an allosteric mechanism explicitly declares whether it operates in R only, T only, with equal R/T kinetics, or with independent R/T kinetics.
2. **DSL support for TR modes.** Users can specify `OnlyR` / `OnlyT` / `EqualRT` / `NonequalRT` directly via `@allosteric_mechanism`.
3. **Hand-verified tests.** PFK-1 and HK as primary tests; rate equations checked against analytical closed forms.
4. **Streamlined enumeration and derivation.** K-type / V-type special-case branches in enumeration collapse into uniform per-group tag moves. The rate-equation derivation collapses parallel R/T-state symbolic derivation into a shared raw-polynomial derivation followed by tag-driven symbol substitution.

---

## 3. Current State (Reference)

### 3.1 Types as they are today

```julia
struct EnzymeMechanism{
    Species,              # (substrates, products, regulators, enzyme_species) 4-tuple
    Reactions,            # ((lhs_syms, rhs_syms),) per step
    EquilibriumSteps,     # (Bool,) parallel to Reactions
    ParamConstraints,     # ((target_sym, coeff::Int, factors::Tuple),)  — general linear
} <: AbstractEnzymeMechanism end

struct AllostericEnzymeMechanism{
    Metabolites,          # tuple of Symbol (full metabolite names, incl reg-only)
    CatalyticMech,        # EnzymeMechanism type (embedded)
    CatSites,             # 7-tuple
    RegSites,             # tuple of 5-tuples
} <: AbstractEnzymeMechanism end
```

### 3.2 DSL as it is today

Two flavors of `@enzyme_mechanism`, switched via `_is_allosteric_label` heuristic. Allosteric flavor cannot express TR tags from the DSL — every user-facing mechanism defaults to fully-independent K_T / K_R.

---

## 4. Target Types

### 4.1 `EnzymeMechanism`

```julia
struct EnzymeMechanism{
    Metabolites,          # ((subs_names,), (prod_names,), (reg_names,))   — each Tuple{Symbol,...}
    Reactions,            # ((lhs_syms, rhs_syms, is_eq::Bool, kinetic_group::Int),)
} <: AbstractEnzymeMechanism end
```

Changes vs. current:

- **`Species` (4-tuple) → `Metabolites` (3-tuple of plain Symbol tuples).** No atom content stored anywhere on the mechanism. Enzyme forms inferred from steps as "any name in steps not in Metabolites." Substrate/product atoms remain at the `EnzymeReaction` level for catalytic-topology enumeration; they are never used at the mechanism level.
- **`Reactions` + `EquilibriumSteps` + `ParamConstraints` collapse into a single tuple of 4-tuples.** Each step carries its `kinetic_group::Int`. Steps with identical `kinetic_group` share kinetic parameters (one `K` for the whole RE group; one `k_f` and one `k_r` for the whole SS group). Steps with unique `kinetic_group` values have independent parameters. No `SameKineticsSteps` type parameter, no general `ParamConstraints` machinery.

At construction time, the `@enzyme_mechanism` / `@allosteric_mechanism` parser assigns `kinetic_group` integers to steps: all steps inside a parenthesized DSL group share the same `kinetic_group`; standalone steps receive unique `kinetic_group` values.

#### 4.1.1 Kinetic-group composition rules (constructor-enforced)

- A kinetic group of 2+ steps contains either all RE steps binding the same metabolite, or all SS binding steps of the same metabolite. Never a mix of RE and SS. Never an iso step.
- Iso steps (SS, no metabolite) always form **singleton groups**. They cannot share kinetic parameters with any other step — Haldane closure relates iso-step k's, not same-kinetics equivalence.

#### 4.1.2 Step canonicalization: NONE

The constructor does **not** sort steps or renumber kinetic groups. Step indices in `Reactions` match the order the user wrote them in the DSL. This preserves predictable parameter naming (`K1` is whatever step the user wrote first as a binding step in its own group, etc.) for hand-derived analytical formulas in tests.

Two semantically equivalent mechanisms written in different DSL orders therefore produce **different Julia types**. This is acceptable: type uniqueness across DSL rewrites is not a load-bearing requirement, and `_canonicalize!` in `mechanism_enumeration.jl`'s `dedup!` path handles canonical form for spec-level dedup independently (sorting steps + renumbering groups before hashing — happens within enumeration, not at type construction).

#### 4.1.3 Constructor validation

The `EnzymeMechanism(metabolites, reactions)` constructor validates:

1. At least one SS step (`is_eq == false`).
2. Each substrate / product / regulator listed appears in some step.
3. Every name in steps is either a listed metabolite or an enzyme form (anything not a listed metabolite is treated as an enzyme form).
4. Each step has exactly one enzyme form on each side.
5. Enzyme-form graph (built from steps) is weakly connected.
6. Kinetic-group composition rules (§4.1.1).
7. **Stoichiometric feasibility**: `r ∈ col(S)` where `S` is the full stoichiometry matrix (enzymes + metabolites) and `r` is the target vector (0 on enzyme rows, `-count(M in subs:)` on substrate rows, `+count(M in prods:)` on product rows, 0 on regulator rows). Verified via `rank(Rational.(S)) == rank(Rational.(hcat(S, r)))`. Catches typos: missing substrate consumption, accidental side-product, regulator with net consumption/production.

No atom-balance check, no enzyme-form atom inference, no substrate-type vs inhibitor-type classification. The user is trusted on chemistry beyond the stoichiometric feasibility check.

### 4.2 `AllostericEnzymeMechanism`

```julia
struct AllostericEnzymeMechanism{
    CatalyticMech,        # embedded EnzymeMechanism type
    CatSites,             # (multiplicity::Int, group_tags::Tuple{Pair{Int,Symbol},...})
    RegSites,             # (((ligands::Tuple{Symbol,...}, multiplicity::Int, ligand_tags::Tuple{Pair{Symbol,Symbol},...}),),)
} <: AbstractEnzymeMechanism end
```

`CatSites` (2 fields, was 7):
- `multiplicity::Int` — catalytic subunit count.
- `group_tags::Tuple{Pair{Int, Symbol}, ...}` — `(kinetic_group, tag)` pairs for groups with non-default TR tags. A group not listed has tag `:NonequalRT`.

`RegSites` per-entry (3 fields, was 5):
- `ligands::Tuple{Symbol, ...}` — ligand names at this reg site.
- `multiplicity::Int` — number of binding sites of this kind.
- `ligand_tags::Tuple{Pair{Symbol, Symbol}, ...}` — `(ligand_name, tag)` pairs.

Tag vocabulary:
- Group tags (catalytic): `:OnlyR`, `:OnlyT`, `:EqualRT`, `:NonequalRT`. Iso groups forbid `:OnlyT` (R-inactive is a relabel).
- Ligand tags (reg sites): `:OnlyR`, `:OnlyT`, `:EqualRT`, `:NonequalRT`. `:EqualRT` is allowed only at multi-ligand reg sites where at least one co-ligand has a non-`:EqualRT` tag (otherwise the ligand's contribution cancels identically).

`Metabolites` (full list) is **derived** by accessor:
```julia
metabolites(m::AllostericEnzymeMechanism) =
    unique union of catalytic-mechanism metabolites and reg-site ligands
```

No information is stored twice. `cat_metabolites` is read from `CatalyticMech`. `same_kinetics_steps` (kinetic groups) live on the `CatalyticMech`'s reactions tuple.

#### 4.2.1 Constructor validation

The `AllostericEnzymeMechanism(catalytic_mech, cat_sites, reg_sites)` constructor validates:

1. `cat_sites.group_tags` references only kinetic groups that exist in `catalytic_mech`.
2. Iso-only groups do not carry `:OnlyT` tags.
3. Each reg-site entry has at least one ligand.
4. No reg site has all `:EqualRT` ligands (cancellation identity); single-ligand `:EqualRT` reg site is forbidden.
5. Reg-site ligand tags are in the allowed vocabulary.

---

## 5. Accessor Interface

**All reads of type-parameter data go through named accessor functions. No `m.parameters[k]`, `CS[k]`, `RS[i][k]`, `Species[k]` lookups in the implementation.** Audit at merge time: grep for `\.parameters\[` / `CS\[` / `RS\[` / `Species\[` returns zero hits in `src/`.

### 5.1 Shared accessors (both `EnzymeMechanism` and `AllostericEnzymeMechanism`)

| Accessor | Description |
|---|---|
| `substrates(m)` | `Tuple{Symbol, ...}` — substrate names |
| `products(m)` | `Tuple{Symbol, ...}` — product names |
| `regulators(m)` | `Tuple{Symbol, ...}` — regulator names (dead-end + allosteric, as a flat list) |
| `metabolites(m)` | `Tuple{Symbol, ...}` — full metabolite list (substrates ∪ products ∪ regulators) |
| `reactions(m)` | Step tuple `((lhs_syms, rhs_syms, is_eq, kinetic_group), ...)` |
| `equilibrium_steps(m)` | `Tuple{Vararg{Bool}}` — extracted from `reactions(m)` |
| `n_steps(m)` | `length(reactions(m))` |
| `enzyme_forms(m)` | `Tuple{Symbol, ...}` — enzyme form names (derived from steps; not stored) |
| `n_states(m)` | `length(enzyme_forms(m))` |
| `kinetic_group(m, step_idx)` | `Int` — the kinetic group assigned to step `step_idx` |
| `kinetic_groups(m)` | Tuple of unique kinetic-group integers in the mechanism |
| `steps_in_group(m, group_num)` | Tuple of step indices in group `group_num` |
| `stoich_matrix(m)` | Full stoichiometry matrix; rows are `(enzyme_forms..., metabolites...)` in that order, columns are step indices. Positive = produced, negative = consumed. |
| `enzyme_row_range(m)` / `metabolite_row_range(m)` | Index ranges for slicing `stoich_matrix(m)` |

### 5.2 Allosteric-specific accessors

| Accessor | Description |
|---|---|
| `catalytic_mechanism(m)` | `EnzymeMechanism` singleton (the embedded `CatalyticMech`) |
| `catalytic_multiplicity(m)` | `Int` — subunit count |
| `group_tag(m, group_num)` | `:OnlyR`/`:OnlyT`/`:EqualRT`/`:NonequalRT`. Defaults to `:NonequalRT` if absent. `:OnlyT` is forbidden for iso groups. |
| `step_tag(m, step_idx)` | Convenience: `group_tag(m, kinetic_group(m, step_idx))`. |
| `allosteric_regulators(m)` | Tuple of `(name, tag)` pairs, one per ligand across all reg sites |
| `catalytic_inhibitors(m)` | Tuple of dead-end-only regulator names (regulators in `CatalyticMech.regulators` that don't appear in any reg site) |
| `regulatory_sites(m)` | Tuple of reg-site descriptors |
| `regulatory_site_ligands(m, i)` | Tuple of ligand names at site `i` |
| `regulatory_site_multiplicity(m, i)` | `Int` |
| `regulatory_ligand_tag(m, i, lig)` | `:OnlyR`/`:OnlyT`/`:EqualRT`/`:NonequalRT` |

Accessors return defaults (`:NonequalRT`) for absent entries — callers never need to distinguish "present but default" from "absent."

---

## 6. DSL

### 6.1 Plain `@enzyme_mechanism`

Tag-free at species level. No `site(...)` blocks, no `::Tag` annotations on steps, no `allosteric_regulators:` / `catalytic_inhibitors:` fields, **no atom syntax** (no `S[C]` brackets — the bracket form is rejected at parse time with a clear error directing users to `@enzyme_reaction` if they want atom declarations).

```julia
@enzyme_mechanism begin
    substrates: S
    products:   P
    regulators: I

    steps: begin
        ([E, S] ⇌ [ES], [EP, S] ⇌ [EPS])    # parenthesized → shared kinetics
        [ES, I] ⇌ [ESI]                      # dead-end inhibitor
        [ES]   <--> [EP]
        [EP]   ⇌   [E, P]
    end
end
```

- `substrates:`, `products:`, `regulators:` accept comma-separated bare symbol lists.
- No `enzymes:` block — enzyme form names are inferred as "any name in steps that is not a listed metabolite."
- No `constraints:` block — same-kinetics groups are expressed by wrapping the sharing steps in a parenthesized tuple. A standalone step is its own group of one. A parenthesized group of steps declares shared kinetics.

### 6.2 `@allosteric_mechanism`

```julia
@allosteric_mechanism begin
    substrates: F6P, ATP
    products:   F16BP, ADP
    allosteric_regulators: Pi::EqualRT, ATP::OnlyT, ADP::OnlyR, Citrate::OnlyT, F26BP::NonequalRT
    catalytic_inhibitors:  J        # optional; dead-end regulators bound only in catalytic steps

    site(:catalytic, 4): begin
        steps: begin
            ([E, F6P] ⇌ [E_F6P], [E_ATP, F6P] ⇌ [E_F6P_ATP])    :: OnlyR
            ([E, ATP] ⇌ [E_ATP], [E_F6P, ATP] ⇌ [E_F6P_ATP])    :: EqualRT
            [E_F6P_ATP]   <--> [E_F16BP_ADP]                      :: EqualRT
            ([E_F16BP_ADP] ⇌ [E_ADP, F16BP], [E_F16BP] ⇌ [E, F16BP]) :: EqualRT
            ([E_F16BP_ADP] ⇌ [E_F16BP, ADP], [E_ADP] ⇌ [E, ADP])     :: EqualRT
        end
    end

    site(:regulatory, 4): begin           # optional — only for competing ligands
        ligands: Pi, ATP                   # tags come from allosteric_regulators: declaration
    end
    # ADP, Citrate, F26BP are tagged but not in an explicit site → each gets its own
    # independent reg site with multiplicity = 4 (the catalytic multiplicity).
end
```

- `substrates:` / `products:` carry **no tags** and **no atoms**. Their T/R binding behavior is fully determined by the kinetic-group tag on each of their binding steps.
- `allosteric_regulators:` carries **required tag** `::Tag` where `Tag ∈ {OnlyR, OnlyT, EqualRT, NonequalRT}`. `EqualRT` is permitted only when a non-`EqualRT` co-ligand exists at the same reg site (validation at construction).
- `catalytic_inhibitors:` carries **no tag**. These dead-end regulators bind in catalytic steps; their TR mode is set by the kinetic-group tag of those steps.
- `site(:catalytic, N):` is **required** in `@allosteric_mechanism` and specifies the catalytic multiplicity. Contains only a `steps:` block.
- `site(:regulatory, N):` is **optional**. Used to declare allosteric regulators that compete at one site (shared partition-function factor). Non-listed allosteric regulators get their own independent reg site with multiplicity = catalytic multiplicity.
- **Every step or step-group has a required `:: Tag`.** Tag vocabulary `{OnlyR, OnlyT, EqualRT, NonequalRT}`. Iso steps (SS, no metabolite) forbid `OnlyT`.

### 6.3 Detection: two macros, no heuristics

`@enzyme_mechanism` and `@allosteric_mechanism` are two separate macros. No heuristic detection. Each rejects the other's exclusive syntax with clear errors.

---

## 7. Tag Semantics

### 7.1 Allosteric regulator tag (at reg sites)

| Tag | Meaning at reg site |
|---|---|
| `::OnlyR` | Ligand binds only in R state (no T-state K) |
| `::OnlyT` | Ligand binds only in T state (no R-state K) |
| `::EqualRT` | K_T = K_R. Allowed only when at least one co-ligand at the same site is non-`EqualRT`. |
| `::NonequalRT` | Independent K_T, K_R |

### 7.2 Kinetic-group tag (catalytic steps)

| Tag | RE binding group | SS binding group | SS iso group |
|---|---|---|---|
| `:: OnlyR` | K_T absent (T doesn't bind here) | k_T_f = k_T_r = 0 (T doesn't bind productively) | k_T_f = k_T_r = 0 (T doesn't catalyze) |
| `:: OnlyT` | K_R absent (R doesn't bind here) | k_R_f = k_R_r = 0 | **Forbidden** (R-inactive is a relabel) |
| `:: EqualRT` | K_T = K_R | k_T_f = k_f, k_T_r = k_r | k_T_f = k_f, k_T_r = k_r |
| `:: NonequalRT` | Independent K_T, K_R | Independent T, R k's | Independent T, R k's |

### 7.3 Tag is a property of the kinetic group, not of metabolites or steps individually

Substrates, products, and catalytic inhibitors carry no metabolite-level TR tags. TR tags live on kinetic groups. A standalone step is its own group of one. Different groups binding the same metabolite can carry different tags.

### 7.4 Reg-site and catalytic contexts are independent

The `allosteric_regulators:` tag governs reg-site binding only. If the same regulator also appears in catalytic steps (e.g., as a dead-end inhibitor), the catalytic kinetic-group tag governs that context independently. ATP in PFK is `::EqualRT` as a catalytic substrate binding and `::OnlyT` as an allosteric regulator — two contexts, two specifications.

---

## 8. Error Cases (DSL-level)

All raised at macro-expansion or constructor time with specific diagnostic messages.

### 8.1 Kinetic-group composition

- Group of 2+ steps contains different metabolites (K/k values not interchangeable).
- Group of 2+ steps contains an iso step (iso steps must be singletons; Haldane handles cycle closure).
- Group of 2+ steps mixes RE and SS binding (different parameter structures).

### 8.2 Tag validation

- Iso group tagged `:: OnlyT`.
- `allosteric_regulators:` ligand tagged `::EqualRT` at a single-ligand reg site, or where every co-ligand is also `::EqualRT` (cancellation identity).
- `catalytic_inhibitors:` entry carrying a tag.
- Step or step-group in `@allosteric_mechanism` without a `:: Tag`.

### 8.3 DSL syntax

- Atom bracket syntax (`S[C]`) at mechanism level: rejected with message directing user to `@enzyme_reaction` for atom declarations.
- `@enzyme_mechanism` block contains `site(...)` / `::Tag` / `allosteric_regulators:` / `catalytic_inhibitors:` — wrong macro.
- `@allosteric_mechanism` block missing `site(:catalytic, N):`.

### 8.4 Stoichiometric feasibility

- `rank(S) ≠ rank([S | r])` — mechanism cannot implement the declared net reaction with the declared substrate / product / regulator stoichiometry. Error message identifies the residual.

---

## 9. Dead Code / Unused Surface to Delete

Confirmed unused or replaceable in production code paths:

### 9.1 `src/types.jl`
- `graph()` accessor (only consumer is a non-allocation test).
- `RegulatorRole` abstract type and concrete `Allosteric` / `DeadEnd` / `UnconstrainedRegulator` subtypes (never dispatched on; symbols `:unknown`/`:dead_end`/`:allosteric` are used everywhere).
- `_count_side` per-step atom counter (no atoms stored at mechanism level).
- `param_constraints(::AllostericEnzymeMechanism) = ()` stub.
- 4-parameter `EnzymeMechanism{Species, Reactions, EquilibriumSteps, ParamConstraints}` and 4-parameter `AllostericEnzymeMechanism{Metabolites, ...}` — replaced by the 2-parameter and 3-parameter shapes in §4.

### 9.2 `src/dsl.jl`
- `_walk_rhs!` / `_parse_constraint_rhs` / `_push_constraint!` — constraint-DSL parser (no `constraints:` block in new DSL).
- `_is_allosteric_label` — heuristic detection (replaced by two distinct macros).
- The `metabolites:` / `site(...):` allosteric flavor of `_parse_enzyme_mechanism`.
- Bracket-atom parsing within mechanism DSL — atoms only appear in `@enzyme_reaction`.

### 9.3 `src/sym_poly_for_rate_eq_derivation.jl`
- `_rs_tr_equiv` / `_rs_r_only` / `_rs_t_only` magic-index helpers.
- `_count_allosteric_rate_monomials` — replaced by direct counting on the post-substitution expression tree.
- `_apply_param_constraints` methods on `POLY` / `FactoredSigma` / `FactoredPoly` / `DenomTerm` (and the recursive walker). Replaced by `_rename_symbols(poly, rename_map)` — simple symbol substitution driven by kinetic-group representatives.

### 9.4 `src/rate_eq_derivation.jl`
- `_is_tr_equiv_catalytic_K`, `_is_tr_equiv_catalytic_param`, `_is_r_only_catalytic_param` — tag-routing helpers replaced by direct accessor calls.
- `_dependent_param_exprs(::AllostericEnzymeMechanism)` parallel-path implementation (~100 lines). Replaced by a single function deriving R-state from the embedded plain-mechanism path and producing T-state via tag-driven symbol substitution at POLY level (zero `:OnlyT` syms for R, zero `:OnlyR` syms + rename `:NonequalRT` syms for T).
- `_allosteric_dep_assignments`, `_allosteric_num_den_exprs`, `_build_allosteric_rate_body` — collapsed into a single rate-expression builder.
- `_kcat_forward(::AllostericEnzymeMechanism)` (~150 lines of magic-index access) — rewritten using a shared `_kcat_from_poly(num, den, params)` helper that also serves the plain `EnzymeMechanism` path.

### 9.5 `src/mechanism_enumeration.jl`
- `_valid_allosteric_differentiations` K-type / V-type hardcoded branches — replaced by uniform per-group tag enumeration.
- `_rewrap_allosteric` — replaced by single methods parametric over spec type via `_steps(s)` / `_with_steps(s, ...)` accessors.
- `_tr_equiv_met_delta`, `_constrained_step_indices` — group structure makes these unnecessary.
- `_is_mirror_of` — mirror propagation is subsumed by kinetic-group atomicity (a group's RE→SS conversion is atomic; if `init_mechanisms` puts catalytic + dead-end-mirror steps in the same kinetic group, they convert together by construction).
- The `::AllostericMechanismSpec` variants of `_expand_re_to_ss`, `_expand_remove_constraint`, `_expand_add_dead_end_regulator`, `_expand_to_allosteric`, `_expand_add_allosteric_regulator`, `_expand_remove_tr_equiv` — collapsed into single methods.
- `compile_mechanism` remains as an internal dispatcher; **removed from the export list**. Internal callers (`identify_rate_equation.jl`, the enumeration pipeline) keep using it.

### 9.6 Documentation drift
- `.claude/CLAUDE.md` references `src/old_mechanism_enumeration.jl` / `src/old_beam_enumeration.jl` / `test/old_test_*.jl` — these files no longer exist; the references are stale and removed in the migration step.

---

## 10. Testing Strategy

### 10.1 Unit tests for the new DSL

- Tag-free plain `@enzyme_mechanism` mirrors existing semantics for migrated mechanisms.
- `@enzyme_mechanism` rejects atom bracket syntax (`S[C]`) — test the error message.
- `@enzyme_mechanism` rejects `site(...)` / `::Tag` / `allosteric_regulators:` / `catalytic_inhibitors:` — wrong-macro detection.
- `@allosteric_mechanism` rejects tag-free step-groups; rejects iso step `:: OnlyT`; rejects single-ligand `::EqualRT` reg sites; rejects all-`::EqualRT` reg site.
- Step-group validation: cross-metabolite group, iso-step in a 2+ group, RE+SS mix in one group.
- Stoichiometric feasibility: a deliberately-broken mechanism (substrate listed but never consumed) errors; the error message is informative.

### 10.2 Hand-verified rate-equation tests (new)

Two realistic mechanisms exercise the full feature set together. Each is compared against a hand-derived closed-form rate equation at multiple concentration / parameter points; agreement to floating-point tolerance is required.

#### 10.2.1 PFK-1 (phosphofructokinase-1)

Reaction: `F6P + ATP ⇌ F16BP + ADP`. Random-order bi-bi at catalytic site; all binding steps RE, iso step SS. Oligomeric state 4.

- **Catalytic-site TR modes**:
  - F6P binding groups: `:: OnlyR` (K-type allosteric — F6P absent from T).
  - ATP, F16BP, ADP binding groups: `:: EqualRT`.
  - Iso step group: `:: EqualRT` (Vmax equal in R and T).
- **Allosteric regulators**:
  - Pi: `::EqualRT`, competes with ATP at reg site 1 (allowed because ATP is `::OnlyT` at the same site — non-cancelling).
  - ATP: `::OnlyT` at reg site 1 (also a substrate; catalytic and reg-site contexts have independent tags).
  - ADP: `::OnlyR`, own reg site.
  - Citrate: `::OnlyT`, own reg site.
  - F26BP: `::NonequalRT`, own reg site.

Tests:
- `rate_equation_string(m)` matches hand-derived form.
- At saturating F6P / ATP and zero allosteric regulators, rate matches Vmax × E_total.
- At F6P → 0, rate → 0 (F6P OnlyR; no T-path).
- ADP, Citrate, ATP-as-regulator, F26BP each modulate rate via R/T equilibrium shift in the directions implied by their tags.

#### 10.2.2 HK (hexokinase)

Reaction: `Glucose + ATP ⇌ G6P + ADP`. Random-order bi-bi; all-`:: EqualRT` catalytic groups. Oligomeric state 2.

- **Catalytic inhibitor**: G6P. The same G6P symbol is in `products:` AND `catalytic_inhibitors:`. Dead-end binding steps `[E_ATP, G6P] ⇌ [E_ATP_G6P]` and `[E_ADP, G6P] ⇌ [E_ADP_G6P]` are tagged `:: EqualRT`.
- **Allosteric regulators**:
  - G6P (third role): `::OnlyT`, competes with Pi at reg site 1.
  - Pi: `::EqualRT`, competes with G6P at reg site 1 — allowed because G6P is `::OnlyT`.

Tests:
- `rate_equation_string(m)` matches hand-derived form. G6P appears in three independent contributions: catalytic-product release, catalytic-inhibitor dead-end, reg-site T-state binding.
- High G6P inhibits across all three mechanisms.
- Reg-site 1 partition functions: at [G6P] = 0, the Pi term cancels in the R/T ratio (both reg_Q's contain only `1 + Pi/K_Pi`); at [G6P] > 0, Pi competes with G6P for the T-state site and shifts kinetics non-trivially.

### 10.2.3 Single-feature edge-case tests

- OnlyT substrate (exotic but well-defined): rate → 0 as K_T → ∞.
- V-type only: every iso group `:: OnlyR`; verify T-state numerator is zero.
- Constructor errors:
  - `::EqualRT` at single-ligand reg site → error.
  - Iso group tagged `:: OnlyT` → error.
  - All-`::EqualRT` reg site → error.
  - Atom bracket in mechanism DSL → error.
  - `same_kinetics` group across different metabolites / RE+SS / iso-included → error.
  - Stoichiometrically infeasible mechanism (substrate listed but never consumed) → error.

### 10.3 Enumeration invariants

`init_mechanisms(reaction)` and `expand_mechanisms(specs, reaction)` count checks for bi-bi, ter-ter, pyruvate carboxylase, pyruvate dehydrogenase. Counts may shift due to:
- Per-group tag enumeration (replaces hard-coded K-type / V-type subsets).
- Atomic kinetic-group RE→SS conversion (replaces independent-mirror conversion).

The new counts are recorded as the post-refactor expected baseline. Denis reviews any divergence from the old counts.

### 10.4 kcat / rescale invariants

Existing kcat-analytical-formula and `rescale_parameter_values` tests apply unchanged for plain mechanisms. For allosteric mechanisms (PFK + HK), the new `_kcat_forward` (built from the shared `_kcat_from_poly` helper) is verified to:
- Match a hand-derived kcat formula on at least one test point per mechanism.
- Make `rescale_parameter_values(m, params; kcat=target)` produce `rate_equation(m, sat_concs, rescaled_params) ≈ E_total * target`.

### 10.5 Aqua / JET

Unchanged. Type stability checks cover the new accessors. Aqua's stale-deps check passes (after potential `Graphs` removal — see §10.6).

### 10.6 `Graphs.jl` dependency

If the rank-based stoichiometric feasibility check (§4.1.3) replaces the only consumer of `Graphs`, the `Graphs` dependency is removed from `Project.toml`. Checked during the migration audit.

---

## 11. Migration Notes

This is a **breaking change** to `EnzymeMechanism`, `AllostericEnzymeMechanism`, the `@enzyme_mechanism` macro grammar, and adds `@allosteric_mechanism` as a new macro.

- `compile_mechanism(spec::MechanismSpec)` and `compile_mechanism(spec::AllostericMechanismSpec)` remain as **internal** dispatchers used by `identify_rate_equation` and the enumeration pipeline. They are not exported. `compile_mechanism` is removed from the SPEC.md exports table and from `src/EnzymeRates.jl` exports. Users construct mechanisms via the DSL macros or via the typed constructors `EnzymeMechanism(metabolites, reactions)` / `AllostericEnzymeMechanism(catalytic_mech, cat_sites, reg_sites)`.
- Existing `@enzyme_mechanism` call sites that used the allosteric flavor (`metabolites:` + `site(:catalytic, N):`) migrate to `@allosteric_mechanism`.
- Existing `species: begin ... end` blocks split into top-level `substrates:` / `products:` / `regulators:` lines; `enzymes:` block deleted (forms inferred from steps).
- Existing `constraints: K2 = K1` blocks translate mechanically into parenthesized step-groups: the steps whose K's were equal go inside one parenthesized tuple in the `steps:` block.
- Existing atom syntax `S[C]` in mechanism-level DSL is dropped — atoms are declared in `@enzyme_reaction` only.
- The three flat-written-out homodimer test mechanisms (`MWC Dimer`, `Homodimer + Non-competitive Inhibitor`, `MWC Dimer + Independent Inhibitor` — non-`[AllostericEnzymeMechanism]` variants) are **deleted**. Their `[AllostericEnzymeMechanism]` siblings (using `site(:catalytic, 2):`) are the canonical encoding and stay.
- Any caller doing `.parameters[k]` indexing on mechanism types switches to named accessors.

All migrations are mechanical. Test files carry the bulk of the migration work.

---

## 12. Out of Scope

- Changes to `rate_equation`, `rate_equation_string`, `parameters`, `fit_rate_equation`, `identify_rate_equation`, `FittingProblem`, `IdentifyRateEquationProblem`, `IdentifyRateEquationResults` public signatures.
- Changes to the rate equation's mathematical semantics.
- Changes to the beam-search pipeline (though `mechanism_enumeration.jl` internals simplify).
- First-class `Vmax` parameter — V-type allostery remains expressible via uniform per-iso-step `:: OnlyR` tagging.
- Per-k_f / per-k_r granularity in TR-equivalence — not expressible. The old code's asymmetric `kNf = kMf` (without matching `kNr = kMr`) constraints appear only in the three flat homodimer fixtures, which are deleted (their `[AllostericEnzymeMechanism]` siblings are equivalent under MWC multiplicity).

---

## 13. Sequence

1. **User reviews this design** and requests edits.
2. Implementation plan (separate doc) executes the design.
3. Plan execution PRs should foreground the code-reduction metric; if a PR doesn't delete code, the design revisits before proceeding.
