# Mechanism Types Refactor — Design

**Status:** Draft — pending user review before implementation planning
**Date:** 2026-04-23
**Scope:** `EnzymeMechanism`, `AllostericEnzymeMechanism`, `@enzyme_mechanism`, `@allosteric_mechanism`, and all code that reads their internals.

---

## 1. Primary Goal: Remove Code

**Less code == simpler code.** The refactor is a success only if a substantial fraction of the affected files is deleted. The current fragmented representations (7-field `CatSites` tuples accessed by magic indices `CS[3]`, `CS[5]`, `CS[6]`, `CS[7]`; 5-per-entry `RegSites`; parallel `tr_equiv_metabolites` / `tr_equiv_cat_steps` / `r_only_metabolites` / `t_only_metabolites` / `r_only_cat_steps` fields; general `ParamConstraint` machinery used only for trivial equalities; duplicated R-state/T-state derivation paths) encode the same underlying information in multiple places and force conditional branches throughout the codebase.

A unified internal representation + accessor-only interface should let us delete, not just move, a large fraction of:

- `src/types.jl` (~574 lines) — magic-index accessors, duplicated structural accessors for the two types, `CatSites`/`RegSites` ad-hoc decoding.
- `src/rate_eq_derivation.jl` (~1591 lines) — parallel R-state/T-state symbol-routing logic, `_is_tr_equiv_catalytic_param` / `_is_r_only_catalytic_param` helpers that re-derive per-param-context what tags should declare directly.
- `src/mechanism_enumeration.jl` (~2430 lines) — `_valid_allosteric_differentiations` K-type/V-type hardcoded branches; `_tr_equiv_met_delta`, `_rewrap_allosteric`, `_canonicalize!` step-index remapping for `ParamConstraints`; separate `_expand_re_to_ss` / `_expand_remove_constraint` / `_expand_add_dead_end_regulator` / `_expand_to_allosteric` / `_expand_add_allosteric_regulator` / `_expand_remove_tr_equiv` functions that each handle both `MechanismSpec` and `AllostericMechanismSpec` variants with near-duplicate code.
- `src/sym_poly_for_rate_eq_derivation.jl` — `_rs_tr_equiv` / `_rs_r_only` / `_rs_t_only` magic-index helpers for `RegSites`.
- `src/dsl.jl` — the two DSL grammars (plain and allosteric) currently branch early via `_is_allosteric_label` detection; shared helpers are minimal.

Target outcome: one unified representation per type, zero magic-index access, a single rate-equation derivation path (the R-state path of `AllostericEnzymeMechanism` *is* plain-`EnzymeMechanism` derivation called directly), and a single DSL grammar with MWC-specific extensions layered on top.

**If the implementation PR doesn't show a net reduction in the touched files, the design is still carrying redundancy somewhere and should be revisited.**

---

## 2. Motivation

Beyond code reduction, this refactor addresses four concrete capability gaps:

1. **Rigorous per-step TR-binding-mode taxonomy.** Each step in an allosteric mechanism explicitly declares whether it operates in R only, T only, with equal R/T kinetics, or with independent R/T kinetics. Currently TR-mode information is fragmented across five parallel fields and inaccessible from the DSL.
2. **DSL support for TR modes.** Users can currently write `@enzyme_mechanism` for `AllostericEnzymeMechanism` but cannot specify `OnlyR` / `OnlyT` / `EqualRT` / `NonequalRT` — every user-facing mechanism defaults to fully-independent K_T / K_R, and refined mechanisms are only reachable via enumeration expansion.
3. **Tests for TR semantics.** New tests cover derived rate equations against hand-written closed forms for: OnlyR substrate, OnlyT substrate, OnlyR / OnlyT / NonequalRT regulator, EqualRT regulator at catalytic dead-end binding, substrate-also-regulator with different tags per context.
4. **Streamlined enumeration and derivation.** The mechanism enumeration pipeline's K-type/V-type special-case branches collapse into uniform per-step/per-species tag moves. The rate-equation derivation's T-state path becomes an expression substitution on the R-state path, not a duplicated symbolic derivation.

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
    CatSites,             # 7-tuple: (cat_mets, multiplicity, tr_equiv_mets,
                          #           tr_equiv_cat_steps, r_only_mets, t_only_mets, r_only_cat_steps)
    RegSites,             # tuple of 5-tuples: ((ligands,), multiplicity, tr_equiv_ligs,
                          #                       r_only_ligs, t_only_ligs)
} <: AbstractEnzymeMechanism end
```

### 3.2 DSL as it is today

Two flavors of `@enzyme_mechanism`, switched via `_is_allosteric_label` heuristic:
- Plain: `species: / steps: / constraints:` blocks. `enzymes:` explicitly lists forms. `constraints:` uses general monomial equality.
- Allosteric: `metabolites:` + `site(:catalytic, N): / site(:regulatory, N):` blocks. **No DSL for `OnlyR` / `OnlyT` / `EqualRT`** — all mechanisms default to `NonequalRT`.

### 3.3 Symptoms of the current design

- Magic-index reads sprinkled throughout (`CS[3]`, `RS[i][5]`, `Species[4]`, etc.).
- `_valid_allosteric_differentiations` generates K-type and V-type mechanisms as separate hardcoded branches.
- `_dependent_param_exprs(M::Type{<:AllostericEnzymeMechanism})` duplicates the R-state logic with `_T`-suffixed symbols, conditionally zeroing `r_only_syms` / `t_only_syms`.
- `_expand_re_to_ss`, `_expand_remove_constraint`, `_expand_add_dead_end_regulator`, etc. each have a `::MechanismSpec` and a `::AllostericMechanismSpec` method with near-duplicate bodies, joined by `_rewrap_allosteric`.

---

## 4. Target Types

### 4.1 `EnzymeMechanism`

```julia
struct EnzymeMechanism{
    Metabolites,          # ((substrates,), (products,), (regulators,)); each ((name, atoms),)
    Reactions,            # ((lhs_syms, rhs_syms, is_eq::Bool, kinetic_group::Int),)
} <: AbstractEnzymeMechanism end
```

Changes vs. current:
- `Species` (4-tuple) → `Metabolites` (3-tuple of substrates/products/regulators). Enzyme forms no longer stored in type parameters; derived from `Reactions` via atomic-balance bookkeeping (already exists).
- `Reactions` + `EquilibriumSteps` + `ParamConstraints` collapse into a single tuple of 4-tuples. Each step carries its `kinetic_group::Int`. Steps with identical `kinetic_group` share kinetic parameters (one `K` for the whole RE group, one `k_f` and one `k_r` for the whole SS group). Steps with unique `kinetic_group` values have independent parameters. No parallel `SameKineticsSteps` type parameter, no general-linear `ParamConstraints` machinery.

At construction time, the `@enzyme_mechanism` / `@allosteric_mechanism` parser assigns `kinetic_group` integers to steps: all steps inside a parenthesized DSL group share the same `kinetic_group`; standalone steps receive unique `kinetic_group` values. Canonicalization (sort by `_step_sort_key`) preserves group membership because `kinetic_group` travels on the step tuple itself — no separate index mapping needed.

### 4.2 `AllostericEnzymeMechanism`

```julia
struct AllostericEnzymeMechanism{
    CatalyticMech,        # embedded EnzymeMechanism type
    CatSites,             # (multiplicity, group_tags) — non-default-only storage
    RegSites,             # (((ligands,), multiplicity, ligand_tags),)  — one entry per reg site
} <: AbstractEnzymeMechanism end
```

`CatSites` fields (2, down from 7):
- `multiplicity::Int` — catalytic subunit count (was `CS[2]`).
- `group_tags::Tuple{Pair{Int, Symbol}, ...}` — `(kinetic_group, tag)` pairs for kinetic groups with non-default TR tags. A group with no pair has tag `:NonequalRT`. Keyed by `kinetic_group` (not `step_idx`), so a group of 3 substrate-binding steps needs one entry, not three.

No `species_tags` field: substrates, products, and catalytic inhibitors carry no metabolite-level TR tags (group tags govern). Allosteric-regulator tags live in `RegSites.ligand_tags`.

Tag vocabulary: `:OnlyR`, `:OnlyT`, `:EqualRT`, `:NonequalRT`. Only three valid for iso steps (`:OnlyT` forbidden — it is a mere relabeling of `:OnlyR` since the theory is asymmetric in R/T).

`RegSites` entry fields (3, down from 5):
- `ligands::Tuple{Symbol, ...}` — ligand names at this reg site.
- `multiplicity::Int` — number of binding sites of this kind.
- `ligand_tags::Tuple{Pair{Symbol, Symbol}, ...}` — `(ligand_name, tag)` pairs. Tag vocabulary at reg sites: `:OnlyR`, `:OnlyT`, `:NonequalRT` (no `:EqualRT` — cancels identically at reg sites).

`Metabolites` is **derived** (not stored):

```julia
metabolites(m::AllostericEnzymeMechanism) =
    unique union of catalytic-mechanism metabolites and reg-site ligands
```

No field carries information derivable from another field. `cat_metabolites` (was `CS[1]`) is read from `CatalyticMech`. `same_kinetics_steps` is on the `CatalyticMech`. Atom content for metabolites is on the `CatalyticMech` (for anything that binds catalytically) or omitted (pure reg-site ligands have no atoms stored).

---

## 5. Accessor Interface

**All reads of type-parameter data go through named accessor functions. No `m.parameters[k]` or `CS[k]` or `RS[i][k]` anywhere in the implementation.** Grep-for-magic-indices becomes part of the verification before merging.

### 5.1 Shared accessors (both types)

| Accessor | Description |
|---|---|
| `substrates(m)` | Tuple of `(name, atoms)` substrate entries |
| `products(m)` | Tuple of `(name, atoms)` product entries |
| `regulators(m)` | Tuple of `(name,)` or `(name, atoms)` regulator entries (inclusive — dead-end and allosteric) |
| `metabolites(m)` | Full metabolite list (union of the above) |
| `reactions(m)` | Step tuple `((lhs, rhs, is_eq, kinetic_group),)` |
| `equilibrium_steps(m)` | `Tuple{Vararg{Bool}}` — extracted from `reactions(m)` |
| `n_steps(m)` | `length(reactions(m))` |
| `enzyme_forms(m)` | Tuple of `(name, atoms)` enzyme forms (derived, not stored) |
| `n_states(m)` | `length(enzyme_forms(m))` |
| `kinetic_group(m, step_idx)` | `Int` — the kinetic group assigned to step `step_idx` |
| `kinetic_groups(m)` | Tuple of unique kinetic-group integers present in the mechanism |
| `steps_in_group(m, group_num)` | Tuple of step indices that share this kinetic group |
| `stoich_matrix(m)` | `Matrix{Int}` (mets × steps) |
| `graph(m)` | `(SimpleDiGraph, enzyme_forms_tuple)` |

### 5.2 Allosteric-specific accessors

| Accessor | Description |
|---|---|
| `catalytic_mechanism(m)` | `EnzymeMechanism` singleton (embedded) |
| `catalytic_multiplicity(m)` | `Int` — subunit count |
| `group_tag(m, group_num)` | `:OnlyR`/`:OnlyT`/`:EqualRT`/`:NonequalRT` — tag for kinetic group `group_num` (`:OnlyT` forbidden for iso groups). Defaults to `:NonequalRT` if not stored. |
| `step_tag(m, step_idx)` | Derived as `group_tag(m, kinetic_group(m, step_idx))`. |
| `allosteric_regulators(m)` | Tuple of `(name, reg_site_tag)` pairs |
| `catalytic_inhibitors(m)` | Tuple of dead-end inhibitor names |
| `regulatory_sites(m)` | Tuple of site descriptors |
| `regulatory_site_ligands(m, i)` | Tuple of ligand names at site `i` |
| `regulatory_site_multiplicity(m, i)` | `Int` |
| `regulatory_ligand_tag(m, i, lig)` | `:OnlyR`/`:OnlyT`/`:NonequalRT` |

Accessors return defaults (`:NonequalRT`) for absent entries in `step_tags` / `ligand_tags` — callers never need to distinguish "present but default" from "absent."

---

## 6. DSL

### 6.1 Plain `@enzyme_mechanism`

Tag-free. No `site(...)` blocks, no `::Tag` annotations, no `allosteric_regulators:` or `catalytic_inhibitors:` fields.

```julia
@enzyme_mechanism begin
    substrates: S[C]
    products:   P[C]
    regulators: I

    steps: begin
      ([E, S] ⇌ [ES], [ES, S] ⇌ [ESS], [EP, S] ⇌ [EPS])   # grouped → same kinetics (K shared)
      [ES, I] ⇌ [ESI]                                      # dead-end
      [ES]   <--> [EP]
      [EP]   ⇌   [E, P]
    end
end
```

- `substrates:`, `products:`, `regulators:` replace the current `species: begin ... end` block; atoms use bracket syntax (`S[C]`, `S[C6H12O6]`, `A[C,N]`).
- `enzymes:` block removed — enzyme form names are the species symbols appearing in steps that are neither substrates, products, nor regulators. Atom content inferred from atomic-balance (existing machinery).
- `constraints:` block removed. **Same-kinetics groups are expressed by wrapping the sharing steps in a parenthesized tuple** — no separate `same_kinetics:` directive, no step duplication. A standalone step is itself; a parenthesized group of steps declares shared kinetics (shared K for RE binding, shared k_f and k_r for SS binding).

### 6.2 `@allosteric_mechanism`

```julia
@allosteric_mechanism begin
    substrates:            S[C]
    products:              P[C]
    allosteric_regulators: I::OnlyT, A::OnlyT, R::NonequalRT
    catalytic_inhibitors:  J

    site(:catalytic, 2): begin
      steps: begin
        ([E, S] ⇌ [ES],
         [ES, S] ⇌ [ESS],
         [EP, S] ⇌ [EPS])         :: EqualRT     # 3 steps, shared kinetics, all EqualRT
        [ES]     <--> [EP]        :: OnlyR       # iso step, T inactive
        [EP]     ⇌   [E, P]       :: EqualRT
        [ES, J]  ⇌   [ESJ]        :: OnlyR       # dead-end J in R state only
      end
    end

    site(:regulatory, 2): begin                  # optional — only for competing ligands
      ligands: A, I
    end
    # R tagged in allosteric_regulators: but not in an explicit reg site →
    # its own independent site, multiplicity = catalytic multiplicity (N = 2 here)
end
```

- `substrates:` and `products:` carry **no tags** (no magic defaults). Their T/R binding behavior is fully determined by step-group tags.
- `allosteric_regulators:` carries **required tag** `::Tag` where `Tag ∈ {OnlyR, OnlyT, NonequalRT}`. No `EqualRT` (cancels at reg sites).
- `catalytic_inhibitors:` carries **no tag**. These regulators appear only in catalytic steps (dead-end binding) and their T/R mode is set by the step-group tag there.
- `site(:catalytic, N):` is **required** and specifies the catalytic multiplicity. Contains only a `steps:` block.
- `site(:regulatory, N):` is **optional**. Used only to declare allosteric regulators that *compete for the same site* (shared partition-function factor). Non-listed allosteric regulators are placed at their own independent site with multiplicity = catalytic multiplicity (the `N` from `site(:catalytic, N):`).
- **Every step or step-group has a required `:: Tag`.** Tag vocabulary: `{OnlyR, OnlyT, EqualRT, NonequalRT}`; iso steps (single, never grouped) forbid `OnlyT`. No site-level aggregate, no defaults.
- A parenthesized group of steps declares shared kinetics and a shared TR tag (both same-kinetics and same TR mode for the group — they can't differ, since grouped steps share K or k_f/k_r values outright, which would be inconsistent with different TR modes).

### 6.3 Detection: two macros, no heuristics

`@enzyme_mechanism` and `@allosteric_mechanism` are two separate macros. No detection/inference of which type to build from block contents. Each macro rejects blocks that don't match its grammar with clear errors.

---

## 7. Tag Semantics

### 7.1 Allosteric-regulator tag (at reg sites)

| Tag | Meaning at reg site |
|---|---|
| `::OnlyR` | Ligand binds only in R state (no T-state K) |
| `::OnlyT` | Ligand binds only in T state (no R-state K) |
| `::EqualRT` | K_T = K_R. Allowed **only** when at least one other ligand at the same reg site has a non-`EqualRT` tag. |
| `::NonequalRT` | Independent K_T, K_R |

**`EqualRT` at a reg site is not an identity.** At a multi-ligand site with mixed tags, an `EqualRT` ligand still affects the rate because its partition-function contribution combines additively with a non-`EqualRT` co-ligand whose R-state and T-state factors differ. Example: PFK reg site with Pi (`EqualRT`) and ATP (`OnlyT`) competing — reg_Q_R = 1 + Pi/K_Pi + ATP_R_missing_so_just_Pi, reg_Q_T = 1 + Pi/K_Pi + ATP/K_ATP_T. The Pi term couples with the ATP term asymmetrically, and Pi's concentration shifts the R/T ratio even though K_Pi is the same in both states.

`EqualRT` **is disallowed** at a reg site if every ligand at that site is `EqualRT` (pure cancellation), or if it is the only ligand at its site.

### 7.2 Kinetic-group tag (catalytic steps)

| Tag | On RE binding group | On SS binding group | On SS iso group |
|---|---|---|---|
| `:: OnlyR` | K_T absent (T state doesn't bind here) | k_T_f = k_T_r = 0 (T doesn't bind productively) | k_T_f = k_T_r = 0 (T doesn't catalyze) |
| `:: OnlyT` | K_R absent (R state doesn't bind here) | k_R_f = k_R_r = 0 | **Forbidden** (R-inactive is a relabel) |
| `:: EqualRT` | K_T = K_R | k_T_f = k_f, k_T_r = k_r | k_T_f = k_f, k_T_r = k_r |
| `:: NonequalRT` | Independent K_T, K_R | Independent T, R k's | Independent T, R k's |

### 7.3 Granularity: kinetic group, not metabolite or step

Substrates, products, and catalytic inhibitors carry no metabolite-level TR tags. TR tags live on kinetic groups (a standalone step forms its own group). Different groups binding the same metabolite can carry different tags (e.g., a group of 3 substrate-binding steps tagged `:: EqualRT` and a standalone substrate-binding step tagged `:: NonequalRT`). Within a group, all steps share *both* kinetics and TR mode — the two properties can't diverge because the DSL binds them in one declaration. Nonsensical configurations surface as unidentifiable-parameter warnings via `structural_identifiability_deficit`, not as DSL errors.

### 7.4 Tag contexts are independent

The `allosteric_regulators:` tag governs **reg-site binding only**. If the same regulator also appears in catalytic steps (dead-end or productive binding), the group tag there governs that context independently. A metabolite like ATP in PFK can be `::EqualRT` as a catalytic substrate binding and `::OnlyT` as an allosteric regulator — two contexts, two specifications.

---

## 8. Error Cases (DSL-level)

All raised at macro-expansion time with specific diagnostic messages:

### 8.1 Step-group validation (same-kinetics rules)

- Group contains steps that bind **different metabolites** — K/k values are not interchangeable across metabolites.
- Group contains an **iso step** — iso cycle closure is handled by Haldane, not manual equality.
- Group **mixes RE and SS** binding steps — different parameter structures (K vs k_f/k_r) cannot be equated.
- Group contains both **substrate-type** binding (metabolite M binding at a form that does not yet contain M) **and inhibitor-type** binding (M binding at a form that already contains M). These use different binding pockets and must have independent K values, even when the metabolite symbol is the same (Denis's edge case: M appears as both substrate and catalytic inhibitor).

### 8.2 Tag validation

- Iso step tagged `:: OnlyT`.
- `allosteric_regulators:` ligand tagged `::EqualRT` at a single-ligand reg site, or at a multi-ligand reg site where every ligand is `::EqualRT` (cancellation identity).
- `catalytic_inhibitors:` entry carrying a tag.
- Step or step-group in `@allosteric_mechanism` without a `:: Tag` (tags required on every step / group).

### 8.3 Macro scope

- `@enzyme_mechanism` block contains `site(...)` / `::Tag` / `allosteric_regulators:` / `catalytic_inhibitors:` — wrong macro.
- `@allosteric_mechanism` block missing `site(:catalytic, N):`.
- Atomic-balance inconsistency across steps (e.g., ping-pong mis-declaration).

---

## 9. Dead Code / Unused Surface to Delete Outright

Audited during brainstorm (2026-04-23/24), all confirmed unused in production code paths:

- **`graph(::EnzymeMechanism)` accessor** (`src/types.jl`). Defined, allocates a `SimpleDiGraph`, tested once in `test_accessors.jl` for non-allocation. Not called anywhere in `src/` computation. Delete the accessor, the `@generated` body, and the graph-allocation test; the `Graphs` import may become unnecessary.
- **`RegulatorRole` type hierarchy** (`src/types.jl`): `abstract type RegulatorRole`, `struct Allosteric <: RegulatorRole end`, `struct DeadEnd <: RegulatorRole end`, `struct UnconstrainedRegulator <: RegulatorRole end`. Every usage in `src/` and `test/` compares against the symbols `:unknown` / `:dead_end` / `:allosteric` directly; the type hierarchy is never dispatched on. Delete.
- **Complex `ParamConstraint` monomial form and the whole `ParamConstraints` type parameter.** User-facing `constraints:` with coefficients or multi-symbol products (`k3r = 2 * k1r`, `k3r = k1f * k2f / k2r`) is exercised only by DSL-parser generality tests in `test/test_dsl.jl`. No test mechanism or use case consumes the monomial form; every real constraint is a simple equality. The 4th type parameter `ParamConstraints` on `EnzymeMechanism`, the `ParamConstraint = Tuple{Symbol, Int, Vector{...}}` definition, the monomial-RHS parser (`_walk_rhs!`, `_parse_constraint_rhs`, `_push_constraint!`), the `coeff::Int, factors::Tuple` carrier, and `_remap_constraint_sym` all disappear. Constraint semantics are carried by the `kinetic_group::Int` stored on each step tuple.
- **`param_constraints(::AllostericEnzymeMechanism) = ()`** stub (`src/types.jl:571`). Returns empty tuple because `AllostericEnzymeMechanism` never carried constraints. Removed along with the `param_constraints` accessor in general (callers switch to `kinetic_group(m, idx)`).

## 9.1 Expected Code-Reduction Targets

Non-binding indicative targets (exact figures depend on implementation — not acceptance criteria, but signal that the design is carrying its weight):

- `src/types.jl`: eliminate magic-index accessors, shrink `CatSites` from 7 to 2 fields, shrink `RegSites` entry from 5 to 3 fields, delete `graph()`, delete `RegulatorRole` hierarchy, delete `param_constraints(::AllostericEnzymeMechanism)` stub. Expect ≥30% reduction.
- `src/rate_eq_derivation.jl`: the T-state derivation collapses into a single substitution pass over the R-state result (driven by tag lookups), eliminating `_is_tr_equiv_catalytic_param`, `_is_r_only_catalytic_param`, and the parallel R/T control flow in `_dependent_param_exprs` and `_build_allosteric_rate_body`. Expect ≥25% reduction.
- `src/mechanism_enumeration.jl`: `_valid_allosteric_differentiations` disappears (K-type / V-type branches become uniform tag enumeration); `_rewrap_allosteric`, `_tr_equiv_met_delta`, and the `::AllostericMechanismSpec` variants of each expansion move collapse into single methods parametric over tag category. Expect ≥30% reduction.
- `src/dsl.jl`: one grammar with a shared helper set (`_parse_species_field`, `_parse_step_or_group`, atom-balance checker) driving both macros. `_walk_rhs!` / `_parse_constraint_rhs` deleted entirely. Expect small net growth *only* if two macros turn out to justify 2× the parser surface — otherwise expect reduction.
- `src/sym_poly_for_rate_eq_derivation.jl`: `_rs_tr_equiv` / `_rs_r_only` / `_rs_t_only` magic-index helpers deleted outright.

**Total: target large net reduction across `src/`. If the implementation doesn't deliver this, it means the new representation is still redundant and we should iterate on the design.**

---

## 10. Testing Strategy

### 10.1 Unit tests for the new DSL

- Tag-free plain `@enzyme_mechanism` mirrors existing semantics (tests migrated, not rewritten).
- `@allosteric_mechanism` rejects tag-free steps and step-groups.
- `@allosteric_mechanism` rejects iso step with `:: OnlyT`.
- `@enzyme_mechanism` rejects `site(...)` / tags / allosteric fields.
- Step-group validation: cross-metabolite group, iso-step in group, RE+SS mix, substrate-type + inhibitor-type mix in one group, group containing a single step (tautological — single steps need no group).

### 10.2 Hand-verified rate-equation tests (new)

Two realistic enzyme mechanisms exercise the full feature set together. Each is compared against a hand-derived closed-form rate equation at multiple concentration points; agreement to floating-point tolerance is required.

#### 10.2.1 PFK-1 (phosphofructokinase-1)

Reaction: `F6P + ATP ⇌ F16BP + ADP`.

Mechanism:
- Random-order bi-bi at catalytic site; all binding steps RE, iso step SS.
- Oligomeric state 4.
- **Catalytic-site TR modes**:
  - F6P binding groups: `:: OnlyR` (K-type allosteric — F6P absent from T).
  - ATP, F16BP, ADP binding groups: `:: EqualRT`.
  - Iso step: `:: EqualRT` (Vmax equal in R and T).
- **Allosteric regulators**:
  - Pi: `::EqualRT`, competes with ATP at reg site 1.
  - ATP: `::OnlyT`, competes with Pi at reg site 1. (ATP is also a substrate; catalytic and reg-site contexts have independent tags.)
  - ADP: `::OnlyR`, own reg site.
  - Citrate: `::OnlyT`, own reg site.
  - F26BP: `::NonequalRT`, own reg site.

Tests:
- `rate_equation_string(m)` matches hand-derived form. Pi with `::EqualRT` appears in both R and T partition factors at site 1 (not cancelled because ATP is `::OnlyT` at the same site).
- At saturating F6P and ATP, zero other regulators, rate matches Vmax × E_total.
- At zero F6P, rate is zero (OnlyR substrate — but no allosteric shift can rescue catalysis without R-state bound F6P... actually still zero because no forward flux).
- ADP concentration modulates rate via R-only stabilization; Citrate / ATP modulate via T-only stabilization.
- F26BP produces differential kinetics between R and T (NonequalRT).

#### 10.2.2 HK (hexokinase)

Reaction: `Glucose + ATP ⇌ G6P + ADP`.

Mechanism:
- Random-order bi-bi at catalytic site; all binding steps RE, iso step SS.
- Oligomeric state 2 (default — can also run with 4).
- **Catalytic-site TR modes**: all groups `:: EqualRT` (no K-type or V-type refinement).
- **Catalytic inhibitor**: G6P. G6P is also a product, so the same symbol appears in `products:` and `catalytic_inhibitors:`. G6P competes with ATP and ADP at the catalytic site (binding to the enzyme in states where ATP or ADP are bound, blocking turnover). Per the grouping rules, substrate-type G6P release (from `E_G6P_ADP` → `E_ADP` + G6P) and inhibitor-type G6P binding (to `E_ATP`, `E_ADP`) cannot share a kinetic group.
- **Allosteric regulators**:
  - G6P (same symbol, third role): `::OnlyT`, competes with Pi at reg site 1. Catalytic + catalytic-inhibitor + allosteric tags are all independent.
  - Pi: `::EqualRT`, competes with G6P at reg site 1 — allowed because G6P is `::OnlyT` (non-cancelling).

Tests:
- `rate_equation_string(m)` matches hand-derived form. G6P appears in three distinct partition-function contributions: catalytic-product release, catalytic-inhibitor dead-end, reg-site T-state binding.
- G6P-inhibition curve at fixed Glucose, ATP: matches the analytical substrate-inhibition form plus the allosteric T-stabilization contribution.
- Pi at zero G6P shifts nothing at the reg site (Pi alone with no co-ligand — but here G6P is always present in reg site 1's partition function since `::OnlyT` means "binds T only", not "absent from site"; Pi does shift kinetics because at zero G6P the R/T partition functions still differ only by Pi... wait — if Pi is `::EqualRT` and the only OTHER ligand at the site is G6P `::OnlyT`, then at zero G6P, reg_Q_R = 1 + Pi/K_Pi and reg_Q_T = 1 + Pi/K_Pi (G6P term is zero by concentration). The Pi term cancels in this concentration regime. But when G6P > 0, reg_Q_T gets an extra G6P/K_G6P_T term, and Pi's competition with G6P for the site is meaningful.

### 10.2.3 Ancillary narrow tests

Small targeted tests for individual features, each exercising one dimension:

- OnlyT substrate (exotic but well-defined): verify rate goes to zero as K_T → ∞.
- Iso step `:: OnlyR` alone (V-type with every other group EqualRT): verify numerator at T-state is zero.
- `::EqualRT` at a single-ligand reg site → construction-time error.
- Group mixing substrate-type and inhibitor-type bindings of G6P in HK → construction-time error.

Each fixture includes an `analytical_rate_fn(params, concs)` alongside the mechanism definition; test asserts `rate_equation(m, concs, params)` ≈ `analytical_rate_fn(params, concs)` at multiple concentration points and parameter values.

### 10.3 Enumeration invariants

`init_mechanisms(reaction)` and `expand_mechanisms(specs, reaction)` count checks retained from the existing test suite (bi-bi = 11, ter-ter = 283, pyruvate carboxylase = 312, pyruvate dehydrogenase = 334), with updates for any count changes introduced by uniform per-step tag enumeration (previously K-type was hardcoded as "subsets of substrates + subsets of products"; new implementation enumerates per-step tag assignments, which may produce a superset or subset depending on equivalences).

### 10.4 kcat / rescale invariants

Existing kcat analytical-formula and `rescale_parameter_values` tests extend uniformly over the new tag taxonomy: for every refactored test mechanism, verify that `rescale_parameter_values(m, params; kcat=target)` produces `rate_equation(m, sat, rescaled) ≈ E_total * target`.

### 10.5 Aqua / JET

Runs as before; unused-export and type-stability checks cover the new accessors.

---

## 11. Migration Notes

This is a **breaking change** to `EnzymeMechanism`, `AllostericEnzymeMechanism`, `@enzyme_mechanism`, and adds `@allosteric_mechanism` as a new macro.

- `compile_mechanism(spec::MechanismSpec)` continues to produce `EnzymeMechanism`; `compile_mechanism(spec::AllostericMechanismSpec)` continues to produce `AllostericEnzymeMechanism`. Internals change; public API preserved.
- Existing `@enzyme_mechanism` call sites that used the allosteric syntax (`metabolites:` + `site(:catalytic, N):`) must migrate to `@allosteric_mechanism`.
- Existing `constraints:` blocks (which are all simple equalities in current code) translate mechanically into parenthesized step-groups at the step's declaration site. For example, `constraints: begin K2 = K1 end` alongside steps `[E, S] ⇌ [ES]` and `[ES, S] ⇌ [ESS]` becomes `([E, S] ⇌ [ES], [ES, S] ⇌ [ESS])` (followed by `:: Tag` in allosteric).
- Existing `enzymes:` blocks are deleted (forms inferred from steps).
- `AllostericEnzymeMechanism`'s accessor API (`substrates`, `products`, `regulators`, `metabolites`, `n_states`, `n_steps`, `equilibrium_steps`, `reactions`, `regulatory_sites`) becomes the sole read interface; any caller doing `.parameters[k]` indexing is updated.

All migrations are mechanical. Test files carry the bulk of the migration work (one pass per file).

---

## 12. Out of Scope

- Changes to `rate_equation`, `rate_equation_string`, `parameters`, `fit_rate_equation`, `identify_rate_equation`, `FittingProblem`, `IdentifyRateEquationProblem`, `IdentifyRateEquationResults` public signatures.
- Changes to the rate equation's mathematical semantics.
- Changes to the beam-search pipeline (though `mechanism_enumeration.jl` internals simplify).
- First-class `Vmax` parameter (option C from brainstorm); V-type remains expressible via uniform per-iso-step tagging.
- Per-k_f / per-k_r granularity in TR-equivalence. Not expressible in DSL; users needing this construct mechanisms via the low-level `EnzymeMechanism` / `AllostericEnzymeMechanism` constructors.

---

## 13. Sequence

1. **User reviews this design** and requests edits.
2. Implementation plan drafted (separate doc) via `superpowers:writing-plans`.
3. Plan executed. Execution PRs should foreground the code-reduction metric; if a PR doesn't delete code, the design revisits before proceeding.
