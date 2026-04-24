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
    Reactions,            # ((lhs_syms, rhs_syms, is_eq::Bool),) — RE/SS flag merged
    SameKineticsSteps,    # (((step_idx, step_idx, ...), ...),) — same-metabolite-binding groups only
} <: AbstractEnzymeMechanism end
```

Changes vs. current:
- `Species` (4-tuple) → `Metabolites` (3-tuple of substrates/products/regulators). Enzyme forms no longer stored in type parameters; derived from `Reactions` via atomic-balance bookkeeping (already exists).
- `Reactions` + `EquilibriumSteps` merged into a single tuple of `(lhs, rhs, is_eq)` triples. Removes parallel-tuple bookkeeping.
- `ParamConstraints` → `SameKineticsSteps`. Pure equality groups of step indices. Removes support for the unused `K4 = K1 * K3 / K2` monomial form (all user test constraints are simple equalities; Wegscheider/Haldane closure is machinery-internal).

### 4.2 `AllostericEnzymeMechanism`

```julia
struct AllostericEnzymeMechanism{
    CatalyticMech,        # embedded EnzymeMechanism type
    CatSites,             # (multiplicity, species_tags, step_tags) — non-default-only storage
    RegSites,             # (((ligands,), multiplicity, ligand_tags),)  — one entry per reg site
} <: AbstractEnzymeMechanism end
```

`CatSites` fields (3, down from 7):
- `multiplicity::Int` — catalytic subunit count (was `CS[2]`).
- `species_tags::Tuple{Pair{Symbol, Symbol}, ...}` — `(name, tag)` pairs for species with non-default TR tags. A species with no pair has tag `:NonequalRT`.
- `step_tags::Tuple{Pair{Int, Symbol}, ...}` — `(step_idx, tag)` pairs for steps with non-default TR tags. A step with no pair has tag `:NonequalRT`.

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
| `reactions(m)` | Step tuple `((lhs, rhs, is_eq),)` |
| `equilibrium_steps(m)` | `Tuple{Vararg{Bool}}` — extracted from `reactions(m)` |
| `n_steps(m)` | `length(reactions(m))` |
| `enzyme_forms(m)` | Tuple of `(name, atoms)` enzyme forms (derived, not stored) |
| `n_states(m)` | `length(enzyme_forms(m))` |
| `same_kinetics_steps(m)` | Tuple of tuples of step indices |
| `stoich_matrix(m)` | `Matrix{Int}` (mets × steps) |
| `graph(m)` | `(SimpleDiGraph, enzyme_forms_tuple)` |

### 5.2 Allosteric-specific accessors

| Accessor | Description |
|---|---|
| `catalytic_mechanism(m)` | `EnzymeMechanism` singleton (embedded) |
| `catalytic_multiplicity(m)` | `Int` — subunit count |
| `species_tag(m, name)` | `:OnlyR`/`:OnlyT`/`:EqualRT`/`:NonequalRT` |
| `step_tag(m, idx)` | `:OnlyR`/`:OnlyT`/`:EqualRT`/`:NonequalRT` for binding steps; `:OnlyT` omitted for iso steps (forbidden) |
| `allosteric_regulators(m)` | Tuple of `(name, reg_site_tag)` pairs |
| `catalytic_inhibitors(m)` | Tuple of dead-end inhibitor names |
| `regulatory_sites(m)` | Tuple of site descriptors |
| `regulatory_site_ligands(m, i)` | Tuple of ligand names at site `i` |
| `regulatory_site_multiplicity(m, i)` | `Int` |
| `regulatory_ligand_tag(m, i, lig)` | `:OnlyR`/`:OnlyT`/`:NonequalRT` |

Accessors return defaults (`:NonequalRT`) for absent entries in `species_tags` / `step_tags` / `ligand_tags` — callers never need to distinguish "present but default" from "absent."

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
      [E, S]  ⇌   [ES]
      [ES, I] ⇌   [ESI]                               # dead-end
      [ES]   <--> [EP]
      [EP]   ⇌   [E, P]

      same_kinetics: [E,S]⇌[ES], [EP,S]⇌[EPS]         # optional; same-metabolite only
    end
end
```

- `substrates:`, `products:`, `regulators:` replace the current `species: begin ... end` block; atoms use bracket syntax (`S[C]`, `S[C6H12O6]`, `A[C,N]`).
- `enzymes:` block removed — enzyme form names are the species symbols appearing in steps that are neither substrates, products, nor regulators. Atom content inferred from atomic-balance (existing machinery).
- `constraints:` block removed — replaced by `same_kinetics:` inside `steps:`.
- `same_kinetics:` lists steps by their literal chemistry (`[E,S]⇌[ES]`), robust to canonicalization reorder. Multiple `same_kinetics:` lines declare multiple disjoint groups.

### 6.2 `@allosteric_mechanism`

```julia
@allosteric_mechanism begin
    substrates:            S[C]
    products:              P[C]
    allosteric_regulators: I::OnlyT, A::OnlyT, R::NonequalRT
    catalytic_inhibitors:  J

    site(:catalytic, 2): begin
      steps: begin
        [E, S]   ⇌    [ES]       :: EqualRT
        [ES, S]  ⇌    [ESS]      :: EqualRT
        [EP, S]  ⇌    [EPS]      :: NonequalRT
        [ES]     <--> [EP]       :: OnlyR            # iso step, T inactive
        [EP]     ⇌    [E, P]     :: EqualRT
        [ES, J]  ⇌    [ESJ]      :: OnlyR            # dead-end J in R state only

        same_kinetics: [E,S]⇌[ES], [ES,S]⇌[ESS], [EP,S]⇌[EPS]
      end
    end

    site(:regulatory, 2): begin                     # optional — only for competing ligands
      ligands: A, I
    end
    # R tagged in allosteric_regulators: but not in an explicit reg site →
    # its own independent site, multiplicity = catalytic multiplicity (N = 2 here)
end
```

- `substrates:` and `products:` carry **no tags** (no magic defaults). Their T/R binding behavior is fully determined by the step tags on each of their binding steps.
- `allosteric_regulators:` carries **required tag** `::Tag` where `Tag ∈ {OnlyR, OnlyT, NonequalRT}`. No `EqualRT` (cancels at reg sites).
- `catalytic_inhibitors:` carries **no tag**. These regulators appear only in catalytic steps (dead-end binding) and their T/R mode is set by per-step tags there.
- `site(:catalytic, N):` is **required** in `@allosteric_mechanism` and specifies the catalytic multiplicity. Contains only a `steps:` block.
- `site(:regulatory, N):` is **optional**. Used only to declare allosteric regulators that *compete for the same site* (shared partition-function factor). Non-listed allosteric regulators are placed at their own independent site with multiplicity = catalytic multiplicity (the `N` from `site(:catalytic, N):`).
- **Every step has a required `:: Tag`.** Tag vocabulary: `{OnlyR, OnlyT, EqualRT, NonequalRT}`; iso steps forbid `OnlyT`. No site-level aggregate, no defaults.
- `same_kinetics:` inside `steps:` works the same as plain `@enzyme_mechanism`.

### 6.3 Detection: two macros, no heuristics

`@enzyme_mechanism` and `@allosteric_mechanism` are two separate macros. No detection/inference of which type to build from block contents. Each macro rejects blocks that don't match its grammar with clear errors.

---

## 7. Tag Semantics

### 7.1 Species-level (allosteric regulators only)

| Tag | Meaning at reg site |
|---|---|
| `::OnlyR` | Ligand binds only in R state (no T-state K) |
| `::OnlyT` | Ligand binds only in T state (no R-state K) |
| `::NonequalRT` | Independent K_T, K_R |

(`::EqualRT` at a reg site cancels identically in numerator and denominator → disallowed.)

### 7.2 Step-level (binding or iso)

| Tag | On binding step | On iso step |
|---|---|---|
| `:: OnlyR` | Only R-state binds at this step (no K_T for this step) | Only R catalyzes (k_T_f = k_T_r = 0) |
| `:: OnlyT` | Only T-state binds (no K_R) | **Forbidden** (R-inactive is a relabel) |
| `:: EqualRT` | K_T = K_R (or k_T_f = k_f and k_T_r = k_r for SS) at this step | k_T_f = k_f, k_T_r = k_r at this iso step |
| `:: NonequalRT` | Independent T, R params at this step | Independent T, R params at this iso step |

### 7.3 Per-step independence (no metabolite-level aggregation)

Substrates, products, and catalytic inhibitors carry no metabolite-level TR tags — every binding step has its own `:: Tag` and the step tags are independent. A metabolite can bind at one enzyme form with `:: EqualRT` and at another with `:: NonequalRT`, or even be `:: OnlyR` at one step and `:: OnlyT` at another. This is strictly more flexible than per-species tagging. Nonsensical mixtures (e.g., an unreachable T-state enzyme form due to broken T-state binding path) surface as unidentifiable-parameter warnings via `structural_identifiability_deficit`, not as DSL errors.

### 7.4 Reg-site ligand tag

The `allosteric_regulators:` tag governs **reg-site binding only**. If the same regulator also appears in catalytic steps (dead-end binding), the step tag there governs that context independently. Two contexts, two specifications.

---

## 8. Error Cases (DSL-level)

All raised at macro-expansion time with specific diagnostic messages:

- `same_kinetics:` group mixes different metabolites (K for metabolite M isn't interchangeable with K for M').
- `same_kinetics:` group includes an iso step (cycle closure via Haldane is the right mechanism, not manual equality).
- `same_kinetics:` group mixes RE and SS binding steps (different parameter structures).
- `same_kinetics:` references a step not present in `steps:`.
- A step appears in two different `same_kinetics:` groups (ambiguous).
- Iso step tagged `:: OnlyT`.
- `allosteric_regulators:` ligand tagged `::EqualRT`.
- `catalytic_inhibitors:` entry carrying a tag.
- `@enzyme_mechanism` block contains `site(...)` / `::Tag` / `allosteric_regulators:` / `catalytic_inhibitors:` — wrong macro.
- `@allosteric_mechanism` block missing `site(:catalytic, N):`.
- Step in `@allosteric_mechanism` without a `:: Tag` (tags required on every step).
- Atomic-balance inconsistency across steps.

---

## 9. Expected Code-Reduction Targets

Non-binding indicative targets (exact figures depend on implementation — not acceptance criteria, but signal that the design is carrying its weight):

- `src/types.jl`: eliminate magic-index accessors, shrink `CatSites` from 7 to 3 fields, shrink `RegSites` entry from 5 to 3 fields. Accessor code consolidates. Expect ≥30% reduction.
- `src/rate_eq_derivation.jl`: the T-state derivation collapses into a single substitution pass over the R-state result (driven by tag lookups), eliminating `_is_tr_equiv_catalytic_param`, `_is_r_only_catalytic_param`, and the parallel R/T control flow in `_dependent_param_exprs` and `_build_allosteric_rate_body`. Expect ≥25% reduction.
- `src/mechanism_enumeration.jl`: `_valid_allosteric_differentiations` disappears (K-type / V-type branches become uniform tag enumeration); `_rewrap_allosteric`, `_tr_equiv_met_delta`, and the `::AllostericMechanismSpec` variants of each expansion move collapse into single methods parametric over tag category. Expect ≥30% reduction.
- `src/dsl.jl`: one grammar with a shared helper set (`_parse_species_field`, `_parse_step_line`, `_parse_same_kinetics`) driving both macros. Expect small net growth *only* if two macros turn out to justify 2× the parser surface — otherwise expect reduction.
- `src/sym_poly_for_rate_eq_derivation.jl`: `_rs_tr_equiv` / `_rs_r_only` / `_rs_t_only` magic-index helpers deleted outright.

**Total: target large net reduction across `src/`. If the implementation doesn't deliver this, it means the new representation is still redundant and we should iterate on the design.**

---

## 10. Testing Strategy

### 10.1 Unit tests for the new DSL

- Tag-free plain `@enzyme_mechanism` mirrors existing semantics (tests migrated, not rewritten).
- `@allosteric_mechanism` rejects tag-free steps.
- `@allosteric_mechanism` rejects iso step with `:: OnlyT`.
- `@enzyme_mechanism` rejects `site(...)` / tags / allosteric fields.
- `same_kinetics:` rejects cross-metabolite, RE+SS mix, iso-step inclusion, missing-step ref, duplicate-group-membership.

### 10.2 Hand-verified rate-equation tests (new)

Mechanisms constructed via DSL, then `rate_equation_string(m)` / `rate_equation(m, concs, params)` compared to hand-derived closed forms:

1. Monomeric Michaelis-Menten (baseline — no allosteric features).
2. Homodimer MWC with every step `:: EqualRT` (should match existing reference `rate_mwc_dimer_oligo`).
3. Homodimer with substrate binding `:: OnlyR` at every S-binding step (K-type, substrate absent from T).
4. Homodimer with product binding `:: OnlyR` at every P-binding step.
5. Homodimer with substrate binding `:: OnlyT` — exotic but well-defined.
6. Homodimer with iso step `:: OnlyR` (V-type, T catalytically inactive).
7. Homodimer with allosteric regulator `I::OnlyT` at own reg site.
8. Homodimer with allosteric regulator `I::OnlyR`.
9. Homodimer with allosteric regulator `I::NonequalRT` (independent K_I_R, K_I_T).
10. Homodimer where substrate S is also an allosteric regulator: `substrates: S[C]` + `allosteric_regulators: S::OnlyT` + catalytic steps tagged per case. Compare against a hand-derived expression where S appears in both the catalytic binding partition function and the reg-site partition function.
11. Mixed metabolite with step-level per-step overrides: substrate S with `:: EqualRT` at one form and `:: NonequalRT` at another.

Each test fixture includes an `analytical_rate_fn(params, concs)` alongside the mechanism definition; test asserts `rate_equation(m, concs, params)` ≈ `analytical_rate_fn(params, concs)` at multiple concentration points.

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
- Existing `constraints:` blocks translate to `same_kinetics:` entries (all existing user-written constraints are simple equalities).
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
