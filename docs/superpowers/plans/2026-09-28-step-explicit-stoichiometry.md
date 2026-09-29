# Steps with Explicit Stoichiometry Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace `Step`'s single `bound_metabolite` with explicit `consumed` / `released` metabolite lists, so fused, Theorell–Chance and multi-metabolite steps derive correctly and every mechanism the enumerator emits today is unchanged.

**Architecture:** `Step(from_species, to_species, consumed, released, is_equilibrium)` becomes the only constructor. The derivation reads the lists and never guesses a side (`_step_sides` is deleted); the rate numerator becomes an orientation-free telescoped sum over SS steps. Pure bindings keep `K_`/`kon_`/`koff_` names; every other non-isomerization step is named by form pair. The `Mechanism` constructor orients every non-binding step and rejects parameter-name collisions.

**Tech Stack:** Julia 1.12, the EnzymeRates package (`src/`), `Test`, `TestEnv` for focused runs.

**Spec:** `docs/superpowers/specs/2026-09-28-step-explicit-stoichiometry-design.md`. Background and evidence: `docs/superpowers/specs/2026-09-28-identifiable-enumeration-findings-and-plan.md` (point 11).

## Global Constraints

- No backward compatibility: the four-argument `Step(from, to, metabolite_or_nothing, is_equilibrium)` must not survive, in `src/` or in `test/`, and no helper may re-create it (Denis, 2026-09-28).
- Every mechanism the enumerator emits today keeps its fitted parameter names, its Reduced and Full `rate_equation_string`, and its step and group order (CLAUDE.md "Canonical Step Form" is load-bearing).
- Pure bindings keep `K_M_F`, `kon_M_F`, `koff_M_F` (with the `inh` marker for a competitive-inhibitor copy); every other non-isomerization step is named by form pair: `Kiso_F1_to_F2`, `k_F1_to_F2`, `k_F2_to_F1`.
- All `Parameter → Symbol` rendering goes through `name(p, m)`; no `Symbol("K…")`/`Symbol("k…")` literal outside a renderer (guarded by `test/test_types.jl`).
- `rate_equation` stays allocation-free and under 120 ns per call for every spec (`test_rate_equation_performance` in `test/test_rate_eq_derivation.jl`).
- 92-character lines, 4-space indentation; every new file starts with two `# ABOUTME: ` lines; match surrounding style; never change unrelated whitespace.
- TDD: write the failing test, run it and see it fail, implement, see it pass.
- Only one Julia process at a time (7.7 GB RAM, no swap). The full suite takes about 13 minutes: run it in the background with output to a log and poll the log in a bounded shell loop (`for i in $(seq 1 40); do sleep 30; grep -q "Test Summary" LOG && break; done`); never wait on a monitor.
- Full suite: `julia --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`.
- Focused run of one test file: `julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<file>.jl")'` (skips Aqua/JET; an `UndefVarError` for a helper defined in another test file is an artifact of the focused run). If `TestEnv` is missing: `julia -e 'using Pkg; Pkg.add("TestEnv")'`.
- Tests in `test/test_mechanism_enumeration.jl` follow the three "Enumeration-engine tests" rules in CLAUDE.md.
- Commit after each task with a message ending in
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` and
  `Claude-Session: https://claude.ai/code/session_014ErtWiLZFt2gypmQ95qxFM`. Never skip hooks; never `git add -A`.

## Review Focus

1. A pure release written from the bound form (`Step(EA, E, [], [A], eq)`) must be stored exactly like the binding written forward, including for a conformation-changing release (`E*(A) → E + A`); covered in Task 1 Step 1.
2. A competitive-inhibitor copy of a substrate (`CompetitiveInhibitor(:A)`) and the substrate `Substrate(:A)` must stay distinct in `==`, `hash`, sort keys, `ligand` and names (`K_Ainh_E` vs `K_A_E`), while sharing the concentration `A` in the rate law; covered in Task 1 Step 1 and Task 2's dead-end cases.
3. A covalent ping-pong step: a binding onto a residual form (`E_res + B → EB_res`) is a pure binding; a step whose residual changes is a transformation and `_assert_chemistry_is_iso` must reject it as an `expand_mechanisms` parent; covered in Task 1 Step 1.
4. Allosteric mechanisms (`AllostericMechanism` constructor, per-state derivation, dead inactive states) must derive identically; covered by the golden tests in Task 1/2 and the regression in Task 4.
5. A step with metabolites on both sides at rapid equilibrium (RE Theorell–Chance) and a step that consumes two metabolites must derive correctly, not error; covered in Task 3.

---

### Task 1: Step with explicit consumed/released lists

**Files:**
- Modify: `src/types.jl` (Step struct and accessors ~132-188; `_entry_kind` ~474; `_canonical_iso_direction` ~495; `_canonicalize_iso_groups` ~520; `_step_canonical_key` ~538; `_assert_no_re_ss_duplicate` ~568; `_re_segment_extras` ~590; `_bottomless_re_segment` ~656; `_to_sig(s::Step)` ~896; `_step_from_sig` ~945; `_drop_unbound_regulators` ~1030; `reactions` ~1440; `_render_binding` ~1611)
- Modify: `src/rate_eq_derivation.jl` (`_step_sides` ~147-186 deleted; `_compute_alpha` ~269; `_raw_symbolic_rate_polys` ~335; `_segment_graph_terms` ~434; `_compute_numerator` ~518)
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl` (`_free_enz_set` ~71; `_step_priority` ~103; `_step_lex_key` ~115; `_thermodynamic_constraints` ~175)
- Modify: `src/mechanism_enumeration.jl` (every `Step(` call and every `bound_metabolite` reader, listed in Step 5)
- Modify: `src/dsl.jl` (`_build_step_expr` ~732-752)
- Create: `test/test_step_stoichiometry.jl`
- Modify: `test/runtests.jl` (include the new file after `test_types.jl`)
- Modify: `test/test_types.jl`, `test/test_dsl.jl`, `test/test_rate_eq_derivation.jl`, `test/test_mechanism_enumeration.jl`, `test/test_identify_rate_equation.jl` (every four-argument `Step(` and every `bound_metabolite`)

**Interfaces:**
- Produces (used by Tasks 2–4):
  - `Step(from_species::Species, to_species::Species, consumed::Vector{<:Metabolite}, released::Vector{<:Metabolite}, is_equilibrium::Bool)`
  - `consumed(s::Step)::Vector{Metabolite}`, `released(s::Step)::Vector{Metabolite}`
  - `ligand(s::Step)::Union{Metabolite, Nothing}`, `is_binding(s) = ligand(s) !== nothing`, `is_iso(s) = isempty(consumed(s)) && isempty(released(s))`
  - `_met_sort_key(m::Metabolite)`, `_binds_ligand(free::Species, bound_form::Species, m::Metabolite)::Bool`
  - `_re_weight_ratio(s::Step, K::Symbol; inverse::Bool = false)::POLY`

- [ ] **Step 1: Write the failing unit tests**

Create `test/test_step_stoichiometry.jl`:

```julia
# ABOUTME: Tests for Steps that carry explicit consumed/released metabolite lists:
# ABOUTME: construction, canonical orientation, kinds, names, and orientation-free derivation.

const ER = EnzymeRates
_sp(bound, conf = :E) = ER.Species(ER.Metabolite[bound...], conf)
_sp(bound, conf, res) = ER.Species(ER.Metabolite[bound...], conf, res)

@testset "Step: explicit consumed/released lists" begin
    A, B, P, Q = ER.Substrate(:A), ER.Substrate(:B), ER.Product(:P), ER.Product(:Q)
    E, EA, EQ, EAB = _sp([]), _sp([A]), _sp([Q]), _sp([A, B])

    @testset "pure binding keeps its written orientation" begin
        s = ER.Step(E, EA, [A], ER.Metabolite[], true)
        @test ER.from_species(s) == E && ER.to_species(s) == EA
        @test ER.consumed(s) == ER.Metabolite[A] && isempty(ER.released(s))
        @test ER.ligand(s) == A && ER.is_binding(s) && !ER.is_iso(s)
    end

    @testset "a pure release is stored as the binding it reverses" begin
        s = ER.Step(EA, E, ER.Metabolite[], [A], false)
        @test s == ER.Step(E, EA, [A], ER.Metabolite[], false)
        # conformation change allowed: E*(A) → E + A is the binding E + A → E*(A)
        Estar_A = _sp([A], :Estar)
        r = ER.Step(Estar_A, E, ER.Metabolite[], [A], true)
        @test ER.from_species(r) == E && ER.to_species(r) == Estar_A
        @test ER.ligand(r) == A
    end

    @testset "isomerization and transformations" begin
        iso = ER.Step(EAB, _sp([P, Q]), ER.Metabolite[], ER.Metabolite[], false)
        @test ER.is_iso(iso) && ER.ligand(iso) === nothing && !ER.is_binding(iso)
        fused = ER.Step(EAB, EQ, ER.Metabolite[], [P], false)        # chemistry + release
        @test !ER.is_iso(fused) && ER.ligand(fused) === nothing
        @test ER.released(fused) == ER.Metabolite[P]
        tc = ER.Step(EA, EQ, [B], [P], false)                          # Theorell–Chance
        @test ER.consumed(tc) == ER.Metabolite[B] && ER.released(tc) == ER.Metabolite[P]
        @test ER.ligand(tc) === nothing
        two = ER.Step(E, EAB, [B, A], ER.Metabolite[], true)           # lists are sorted
        @test ER.consumed(two) == ER.Metabolite[A, B] && ER.ligand(two) === nothing
    end

    @testset "covalent residual: binding onto a residual form vs chemistry" begin
        res = ER.Residual([A], [P])
        F, FB = _sp([], :E, res), _sp([B], :E, res)
        @test ER.ligand(ER.Step(F, FB, [B], ER.Metabolite[], true)) == B
        chem = ER.Step(EA, F, ER.Metabolite[], [P], false)            # E(A) → F + P
        @test ER.ligand(chem) === nothing
        m = ER.Mechanism(
            @enzyme_reaction(begin
                substrates: A[CX], B[N]
                products: P[C], Q[NX]
            end),
            [[ER.Step(E, EA, [A], ER.Metabolite[], false)], [chem],
             [ER.Step(F, FB, [B], ER.Metabolite[], false)],
             [ER.Step(FB, E, ER.Metabolite[], [Q], false)]])
        @test_throws ErrorException ER._assert_chemistry_is_iso(m)
    end

    @testset "inhibitor copy stays distinct from the substrate" begin
        Ai = ER.CompetitiveInhibitor(:A)
        s_sub = ER.Step(E, EA, [A], ER.Metabolite[], true)
        s_inh = ER.Step(E, _sp([Ai]), [Ai], ER.Metabolite[], true)
        @test s_sub != s_inh && hash(s_sub) != hash(s_inh)
        @test ER.ligand(s_inh) == Ai
    end

    @testset "rejections" begin
        @test_throws ErrorException ER.Step(E, E, [A], ER.Metabolite[], true)
        @test_throws ErrorException ER.Step(EA, EQ, [A], [A], false)
    end

    @testset "sort key reproduces today's order" begin
        s = ER.Step(E, EA, [A], ER.Metabolite[], true)
        @test ER._step_canonical_key(s) == ("E", "EA", "A", "", true)
        iso = ER.Step(EAB, _sp([P, Q]), ER.Metabolite[], ER.Metabolite[], false)
        @test ER._step_canonical_key(iso) == ("EAB", "EPQ", "", "", false)
    end

    @testset "signature round trip" begin
        tc = ER.Step(EA, EQ, [B], [P], false)
        @test ER._step_from_sig(ER._to_sig(tc)) == tc
    end
end

@testset "transformation steps are named by form pair" begin
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end
    names = Set(ER.parameters(m, ER.Full))
    @test :k_EAB_to_EQ in names && :k_EQ_to_EAB in names
    @test !(:kon_P_EAB in names) && !(:koff_P_EAB in names)
    @test :kon_A_E in names && :koff_Q_E in names
end
```

Add to `test/runtests.jl`, after `include("test_types.jl")`:

```julia
    include("test_step_stoichiometry.jl")
```

(`ER.parameters(m, ER.Full)` returns the full parameter-name tuple; check its exact
return type in `src/rate_eq_derivation.jl` ~43 and adapt the `Set(...)` call if it
returns a NamedTuple or Symbols wrapped differently.)

- [ ] **Step 2: Run the new file and confirm it fails**

Run the focused command with `test/test_step_stoichiometry.jl`.
Expected: FAIL with `MethodError: no method matching EnzymeRates.Step(::Species, ::Species, ::Vector{…}, ::Vector{…}, ::Bool)`.

- [ ] **Step 3: Replace the Step definition and its accessors in `src/types.jl`**

Replace the comment block, docstring and `struct Step … end` plus the accessor lines
(`from_species` … `Base.hash(s::Step, …)`, currently ~132-188) with:

```julia
# Sort key for metabolite lists, shared with `Species.bound`: by name, a
# competitive-inhibitor copy after the metabolite of the same name.
_met_sort_key(m::Metabolite) = (name(m), m isa CompetitiveInhibitor)

# Whether `bound_form` is `free` with metabolite `m` added: same residual, and
# `bound(bound_form)` equals `bound(free)` plus `m` as a multiset. The
# conformation may differ (a binding may change the enzyme's conformation).
_binds_ligand(free::Species, bound_form::Species, m::Metabolite) =
    residual(free) == residual(bound_form) &&
    bound(bound_form) == sort(Metabolite[bound(free)..., m]; by = _met_sort_key)

"""
    Step

One elementary transition between two enzyme `Species`, with its stoichiometry:
going `from_species → to_species` it takes up the metabolites in `consumed` from
solution and gives off those in `released`. Both lists may be empty (an
isomerization) or non-empty (for example a Theorell–Chance step EA + B → EQ + P).
`is_equilibrium` flags a rapid-equilibrium step (`true`) versus a steady-state
step (`false`). A pure binding (`ligand`) is stored with its metabolite consumed,
bound on `to_species`; every other step is oriented by the `Mechanism` /
`AllostericMechanism` constructor. See CLAUDE.md "Canonical Step Form".
"""
struct Step
    from_species::Species
    to_species::Species
    consumed::Vector{Metabolite}
    released::Vector{Metabolite}
    is_equilibrium::Bool
    function Step(from_species::Species, to_species::Species,
                  consumed::Vector{<:Metabolite}, released::Vector{<:Metabolite},
                  is_equilibrium::Bool)
        from_species == to_species && error(
            "Step: both ends are the form $(name(from_species)); a step joins two forms")
        c = sort(Vector{Metabolite}(consumed); by = _met_sort_key)
        r = sort(Vector{Metabolite}(released); by = _met_sort_key)
        shared = intersect(c, r)
        isempty(shared) || error(
            "Step $(name(from_species)) → $(name(to_species)): " *
            "$(join(name.(shared), ", ")) both consumed and released")
        # A pure release is stored as the binding it reverses (metabolite on `to`).
        if isempty(c) && length(r) == 1 && _binds_ligand(to_species, from_species, only(r))
            from_species, to_species, c, r = to_species, from_species, r, c
        end
        new(from_species, to_species, c, r, is_equilibrium)
    end
end

from_species(s::Step)   = s.from_species
to_species(s::Step)     = s.to_species
consumed(s::Step)       = s.consumed
released(s::Step)       = s.released
is_equilibrium(s::Step) = s.is_equilibrium

"""The metabolite of a pure binding (it consumes exactly that metabolite, releases
nothing, and `to_species` is `from_species` with it bound), else `nothing`."""
function ligand(s::Step)
    length(s.consumed) == 1 && isempty(s.released) || return nothing
    m = only(s.consumed)
    _binds_ligand(s.from_species, s.to_species, m) ? m : nothing
end
is_binding(s::Step) = ligand(s) !== nothing
is_iso(s::Step)     = isempty(s.consumed) && isempty(s.released)
direction(s::Step)  = is_binding(s) ? :binding : :iso

Base.:(==)(a::Step, b::Step) =
    a.from_species == b.from_species && a.to_species == b.to_species &&
    a.consumed == b.consumed && a.released == b.released &&
    a.is_equilibrium == b.is_equilibrium
Base.hash(s::Step, h::UInt) =
    hash(s.is_equilibrium, hash(s.released, hash(s.consumed,
        hash(s.to_species, hash(s.from_species, hash(:Step, h))))))
```

Check with `grep -n "direction(" src/*.jl` whether `direction` has callers; delete it
if it has none (CLAUDE.md: remove unused features).

- [ ] **Step 4: Update the identity helpers in `src/types.jl` and the thermodynamics file**

```julia
# src/types.jl — canonical sort key (was bound-metabolite name or "")
_met_names(v) = join((String(name(m)) for m in v), "+")
_step_canonical_key(s::Step) =
    (String(name(from_species(s))), String(name(to_species(s))),
     _met_names(consumed(s)), _met_names(released(s)), is_equilibrium(s))
```

In `_assert_no_re_ss_duplicate` change the key to
`(from_species(s), to_species(s), consumed(s), released(s))` and its `Dict` key type
to `Tuple{Species, Species, Vector{Metabolite}, Vector{Metabolite}}`.

```julia
# src/types.jl — signature encoding
_to_sig(s::Step) = (
    _to_sig(from_species(s)),
    _to_sig(to_species(s)),
    Tuple(_to_sig(m) for m in consumed(s)),
    Tuple(_to_sig(m) for m in released(s)),
    is_equilibrium(s),
)

function _step_from_sig(sig::Tuple)
    from_sig, to_sig, consumed_sig, released_sig, is_eq = sig
    Step(_species_from_sig(from_sig), _species_from_sig(to_sig),
         Metabolite[_metabolite_from_sig(t) for t in consumed_sig],
         Metabolite[_metabolite_from_sig(t) for t in released_sig], is_eq)
end
```

```julia
# src/thermodynamic_constr_for_rate_eq_derivation.jl — group-representative tiebreak
_step_lex_key(s::Step) =
    (String(name(from_species(s))), String(name(to_species(s))),
     _met_names(consumed(s)), _met_names(released(s)), is_equilibrium(s))
```

- [ ] **Step 5: Convert every construction site and reader in `src/`**

Conversion rules for a four-argument call `Step(a, b, x, eq)`:
- `x === nothing` → `Step(a, b, Metabolite[], Metabolite[], eq)`.
- `x` a metabolite bound on `b` (a binding written forward) → `Step(a, b, Metabolite[x], Metabolite[], eq)`.
- `x` a metabolite bound on `a` (a release written from the bound form) → `Step(a, b, Metabolite[], Metabolite[x], eq)` (the constructor stores it as the binding).
- A rebuild that copies an existing step `s` (`Step(f, t, bound_metabolite(s), eq)`) → `Step(f, t, consumed(s), released(s), eq)`.

Sites (line numbers are approximate; search for each):
- `src/mechanism_enumeration.jl` `_release_products!` ~249 and `backtrack!` ~379 (product releases: `Step(cur, new, Metabolite[], Metabolite[Product(p)], true)`); substrate bindings ~407, ~431, ~487, ~529, ~553, ~578, ~633, ~662 (inspect each: `Substrate(s)` bound on the new form → consumed; `nothing` → isomerization).
- `_expand_substrate_product_dead_ends` ~1093, ~1116, ~1155 (`_role(met)` bound on `sp2` → consumed) and the mirror rebuild ~1119.
- `_iso`/grouping rebuild ~808-811: `Step(from_species(s), to_species(s), consumed(s), released(s), i != iso_idx)`.
- `_flip_group_to_ss` ~1290: `Step(from_species(s), to_species(s), consumed(s), released(s), false)`.
- dead-end move ~1930: `Step(base, de_species, Metabolite[CompetitiveInhibitor(reg_name)], Metabolite[], true)`; mirrors ~1943: `Step(de_species_map[fn], de_species_map[tn], consumed(s), released(s), is_equilibrium(s))`.
- `src/types.jl` `_canonical_iso_direction` ~505, ~511, ~514: `Step(t, f, Metabolite[], Metabolite[], is_equilibrium(s))`.

Readers of `bound_metabolite` (replace, keeping behaviour for every pure binding and isomerization):
- `_binding_order` / `_release_order` ~266-274: `name(ligand(s))` for `is_binding(s) && ligand(s) isa Substrate` (resp. `Product`).
- ~771-772: `m = ligand(step); m === nothing && continue`.
- `_apply_equivalence_grouping` ~1205: key each step with free metabolites by `(Tuple(name.(consumed(s))), Tuple(name.(released(s))), is_equilibrium(s))`; skip `is_iso(s)`. The `Dict` key type becomes `Tuple{Tuple, Tuple, Bool}`.
- ~1263 (flip units exclude regulator-binding groups): `any(m -> m isa Regulator, consumed(s)) || any(m -> m isa Regulator, released(s))`.
- `_context_bipartitions` ~1657: `own = ligand(first(group))`.
- `_forms_with_binding_step_native` ~1745: `m = ligand(s); m === nothing && continue` then compare `name(m)`.
- dead-end kernel ~1842 and ~1884: regulator/inhibitor detection over `consumed(s)`.
- `_bound_comp_inhibitors` ~2641: over `consumed(s)`.
- `_assert_mechanism_invariants` ~2685-2692: replace the two branches by
  `is_binding(s) || is_iso(s) || error("step $(name(from_species(s))) → $(name(to_species(s))) is neither a pure binding nor an isomerization")`.
- ~2709, ~2722-2723: read `consumed(s)`/`released(s)`; `kinds` keeps only steps with free metabolites, keyed by `(is_equilibrium(s), consumed(s), released(s))`.
- `_assert_step_atom_conserving` ~139:

```julia
function _assert_step_atom_conserving(reaction::EnzymeReaction, s::Step)
    diff = Dict{Symbol,Int}()
    _accumulate_atoms!(diff, _species_atoms(reaction, to_species(s)), 1)
    _accumulate_atoms!(diff, _species_atoms(reaction, from_species(s)), -1)
    for m in consumed(s); _accumulate_atoms!(diff, _atoms_dict(reaction, name(m)), -1); end
    for m in released(s); _accumulate_atoms!(diff, _atoms_dict(reaction, name(m)), 1); end
    diff = _nonzero_atoms(diff)
    isempty(diff) || error(
        "atom-non-conserving step $(name(from_species(s))) → $(name(to_species(s))) " *
        "(consumed $(name.(consumed(s))), released $(name.(released(s)))): " *
        "atoms(from) + atoms(consumed) − atoms(to) − atoms(released) = $diff")
    nothing
end
```

  (`_atoms_dict(reaction, name)` is the function the old code called on the bound metabolite, so inhibitor steps are checked exactly as before.)
- `_assert_chemistry_is_iso` ~163: a parent passed to the moves may contain only pure bindings and isomerizations:

```julia
function _assert_chemistry_is_iso(m::Union{Mechanism, AllostericMechanism})
    for group in steps(m), s in group
        is_binding(s) || is_iso(s) || error(
            "step $(name(from_species(s))) → $(name(to_species(s))) folds chemistry " *
            "into a binding or release; the moves need the chemistry as an " *
            "isomerization and each binding or release as its own step")
    end
    nothing
end
```

  Update its docstring accordingly.
- `src/types.jl` `_entry_kind` ~474: iterate steps with `from_species(s) == sp` and non-empty `consumed(s)`; classify every consumed metabolite name as substrate/product. In `_canonicalize_iso_groups` ~525 pass `filter(s -> !isempty(consumed(s)), flat0)`.
- `_canonical_iso_direction` ~497: replace `is_binding(s) && return s` by `is_iso(s) || return s` (Task 3 extends orientation to transformations).
- `_bottomless_re_segment` ~659: for every RE step, record `side[name(m)] = typeof(m)` for each `m` in `consumed(s)` and `released(s)`.
- `_drop_unbound_regulators` ~1038: add the names of `consumed(s)` and `released(s)`.
- `_render_binding` ~1612: `met = ligand(rep)`.
- `src/thermodynamic_constr_for_rate_eq_derivation.jl`:
  - `_free_enz_set` ~80: exclude `name(to_species(s))` for every RE step with non-empty `consumed(s)`; `bound_in` collects `to_species` of every step with non-empty `consumed(s)`.
  - `_step_priority` ~104: `has_met = !is_iso(s)`.

- [ ] **Step 6: Delete `_step_sides` and read the lists at each caller**

Delete `_step_sides` and its docstring (`src/rate_eq_derivation.jl` ~147-186). At each caller:
- `_raw_symbolic_rate_polys` ~355: `i_form = enz_name_to_form[name(from_species(s))]`, `j_form = enz_name_to_form[name(to_species(s))]`, and pass `Symbol[name(m) for m in consumed(s)]` / `…released(s)` to `_ss_contrib`.
- `_segment_graph_terms` ~443: the two form names.
- `_thermodynamic_constraints` ~207: `-1` per consumed metabolite, `+1` per released one (keep the `haskey(met_idx, …)` guard); delete the comment block that refers to `_step_sides` and write one line saying the stoichiometry comes from the step's lists.
- `_re_segment_extras` (`src/types.jl` ~601): `d[name(m)] += 1` for consumed, `-= 1` for released.
- `reactions` (`src/types.jl` ~1447): `((name(from_species(s)), name.(consumed(s))...), (name(to_species(s)), name.(released(s))...), is_equilibrium(s), g)`; update its docstring.
- `_hyperbolic_catalysis` (`src/mechanism_enumeration.jl` ~1421): `m_lhs = Symbol[name(m) for m in consumed(s)]`, `m_rhs = …released(s)`.
- `_context_form` ~1634: `isempty(consumed(s)) || return from_species(s); isempty(released(s)) || return to_species(s); from_species(s)`; update its docstring.
- `_compute_numerator` (~536-623): replace `bound_metabolite(s)` by `bm = only(vcat(consumed(s), released(s)))` for steps with exactly one free metabolite (all steps before Task 3), take `m_lhs`/`m_rhs` from the lists, and orient a product step by whether the product is consumed. Task 2 replaces this function.

- [ ] **Step 7: One RE weight rule in `_compute_alpha`**

Add to `src/rate_eq_derivation.jl` and use it in `_compute_alpha` for both traversal
directions (replace the `is_iso` / binding branches):

```julia
"""
    _re_weight_ratio(s, K; inverse = false) -> POLY

w(to)/w(from) of rapid-equilibrium step `s` as a monomial (its inverse when
`inverse`): [M]/K for a pure binding of M (K a dissociation constant), and
K·Π[consumed]/Π[released] for every other step (K in the association direction).
"""
function _re_weight_ratio(s::Step, K::Symbol; inverse::Bool = false)
    sgn = inverse ? -1 : 1
    d = Dict{Symbol, Int}(K => sgn * (is_binding(s) ? -1 : 1))
    for m in consumed(s); d[name(m)] = get(d, name(m), 0) + sgn; end
    for m in released(s); d[name(m)] = get(d, name(m), 0) - sgn; end
    filter!(p -> p.second != 0, d)
    POLY(sort!(MONO(collect(d)); by = first) => 1)
end
```

In `_compute_alpha`, for an RE step with `i_f = form of from_species(s)` and
`j_f = form of to_species(s)`:
`alpha[j_f] = poly_mul(alpha[cur], _re_weight_ratio(s, step_to_K[idx]))` when entering
`to` from `from`, and
`alpha[i_f] = poly_mul(alpha[cur], _re_weight_ratio(s, step_to_K[idx]; inverse = true))`
the other way. `_assemble_constraints` already signs an RE constant −1 exactly for
pure bindings (`binding_K_set` uses `is_binding`), which matches this rule.

- [ ] **Step 8: Convert `src/dsl.jl` `_build_step_expr` (one metabolite per side still)**

```julia
function _build_step_expr(lhs::Vector{_StepSideTerm},
                          rhs::Vector{_StepSideTerm},
                          is_eq::Bool,
                          role_of::Dict{Symbol,Symbol})
    lhs_enzyme, lhs_met = _split_side(lhs)
    rhs_enzyme, rhs_met = _split_side(rhs)
    met_exprs(t) = t === nothing ? Expr[] :
        [_metabolite_expr(t.sym, role_of, t.role)]
    from_expr = _species_expr_from_term(lhs_enzyme, role_of)
    to_expr   = _species_expr_from_term(rhs_enzyme, role_of)
    :(EnzymeRates.Step($from_expr, $to_expr,
                       EnzymeRates.Metabolite[$(met_exprs(lhs_met)...)],
                       EnzymeRates.Metabolite[$(met_exprs(rhs_met)...)], $is_eq))
end
```

Update its docstring (the left-hand metabolite is consumed, the right-hand one released).

- [ ] **Step 9: Convert the tests**

In `test/`, convert every four-argument `EnzymeRates.Step(` with the Step 5 rules
(sites: `test_identify_rate_equation.jl` ~1537-1586, `test_rate_eq_derivation.jl`
~1793-1820, `test_mechanism_enumeration.jl` ~233-358, ~2022-2024, ~6029,
`test_types.jl` ~342-397, ~819), and every `bound_metabolite(s)` (60 uses; `grep -n
bound_metabolite test/*.jl`): use `ligand(s)` where the step is a pure binding, and the
lists where the test reads a step's metabolite generally (e.g. the helpers
`_t_reactants`/`_t_products` at the top of `test_mechanism_enumeration.jl` use
`EnzymeRates.is_iso(s)` and `only(vcat(EnzymeRates.consumed(s), EnzymeRates.released(s)))`).
Tests whose expectations name the parameters of a hand-written fused step
(for example `kon_P_EAB`) change to the form-pair names; record every such rename in
the commit message.

- [ ] **Step 10: Run the focused files, then the full suite**

Focused: `test/test_step_stoichiometry.jl`, `test/test_types.jl`, `test/test_dsl.jl`,
`test/test_rate_eq_derivation.jl`, `test/test_mechanism_enumeration.jl` (one at a time).
Expected: all pass. Then the full suite in the background; expected: 0 failures, 0 errors.
Also confirm `grep -rn "bound_metabolite\|_step_sides" src test` prints nothing.

- [ ] **Step 11: Commit**

```bash
git add src/types.jl src/rate_eq_derivation.jl src/thermodynamic_constr_for_rate_eq_derivation.jl \
        src/mechanism_enumeration.jl src/dsl.jl test/test_step_stoichiometry.jl test/runtests.jl \
        test/test_types.jl test/test_dsl.jl test/test_rate_eq_derivation.jl \
        test/test_mechanism_enumeration.jl test/test_identify_rate_equation.jl
git commit   # message: "Give Step explicit consumed/released lists" + renames + trailers
```

---

### Task 2: Orientation-free rate numerator

**Files:**
- Modify: `src/rate_eq_derivation.jl` (`_compute_alpha`: RE-cycle check; `_compute_numerator` rewritten; `_raw_symbolic_rate_polys` loses `subs_species`, `prods_species`, `allow_dead`; `_state_rate_polys` ~1301; docstrings mentioning `allow_dead` ~1158, ~1692)
- Modify: `test/test_step_stoichiometry.jl` (oracle helper + suite)
- Modify: `test/test_rate_eq_derivation.jl` ("Numerator: ambiguous central cut (regulator sibling) raises" ~1517)

**Interfaces:**
- Consumes: Task 1's `Step` accessors and `_re_weight_ratio`.
- Produces: `_compute_numerator(mech, enz_name_to_form, step_params, alpha, form_to_group, D)::POLY`;
  `_raw_symbolic_rate_polys(mech, step_params, rename_map)` (no metabolite lists, no keyword);
  test helpers `_testhelper_consistent_point(em, rng)` and `_testhelper_mass_action_rate(em, point, concs)`.

- [ ] **Step 1: Write the oracle helpers and the failing correctness tests**

Append to `test/test_step_stoichiometry.jl`:

```julia
# ── Brute-force mass-action oracle (no Cha segments) ────────────────────────
# A consistent parameter point: free energies g for every enzyme form and
# metabolite make each step's association constant
#   Ka = exp(g_from + Σ g_consumed − g_to − Σ g_released),
# so every cycle multiplies to Keq^n by construction. Each step is its own
# kinetic group in the mechanisms used here.
function _testhelper_consistent_point(em, rng)
    mech = ER.Mechanism(em)
    flat = ER._flat_steps(mech)
    forms = unique(vcat([[ER.from_species(s), ER.to_species(s)] for (s, _) in flat]...))
    g = Dict{Any, BigFloat}(f => 4 * big(rand(rng)) - 2 for f in forms)
    rxn = ER.reaction(mech)
    mets = vcat([ER.name(x) for x in ER.substrates(rxn)], [ER.name(x) for x in ER.products(rxn)],
                [ER.name(ER.regulator(r)) for r in ER.regulators(rxn)])
    for x in mets; g[x] = 4 * big(rand(rng)) - 2; end
    Keq = exp(sum(g[ER.name(x)] for x in ER.substrates(rxn)) -
              sum(g[ER.name(x)] for x in ER.products(rxn)))
    steps = []
    vals = Dict{Symbol, BigFloat}()
    sp = ER._step_parameters(mech)
    for (idx, (s, _)) in enumerate(flat)
        Ka = exp(g[ER.from_species(s)] + sum((g[ER.name(m)] for m in ER.consumed(s)); init = big(0)) -
                 g[ER.to_species(s)] - sum((g[ER.name(m)] for m in ER.released(s)); init = big(0)))
        if ER.is_equilibrium(s)
            p = sp[idx][1]
            vals[ER.name(p, mech)] = p isa ER.Kd ? 1 / Ka : Ka
            push!(steps, (s, Ka, nothing))
        else
            kf = exp(4 * big(rand(rng)) - 2)
            vals[ER.name(sp[idx][1], mech)] = kf
            vals[ER.name(sp[idx][2], mech)] = kf / Ka
            push!(steps, (s, Ka, kf))
        end
    end
    (; steps, vals, Keq)
end

# v per unit enzyme from the full linear steady state over every form; RE steps
# run at rate scale Λ. v = net consumption of the first substrate.
function _testhelper_mass_action_rate(em, point, concs::Dict{Symbol, BigFloat};
                                      Λ = big(10)^60)
    setprecision(BigFloat, 512) do
        mech = ER.Mechanism(em)
        forms = unique(vcat([[ER.from_species(s), ER.to_species(s)] for (s, _, _) in point.steps]...))
        idx = Dict(f => i for (i, f) in enumerate(forms))
        n = length(forms)
        M = zeros(BigFloat, n, n)
        rates = []
        for (s, Ka, kf) in point.steps
            f_ = kf === nothing ? Λ : kf
            r_ = kf === nothing ? Λ / Ka : kf / Ka
            rf = f_ * prod((concs[ER.name(m)] for m in ER.consumed(s)); init = big(1))
            rr = r_ * prod((concs[ER.name(m)] for m in ER.released(s)); init = big(1))
            a, b = idx[ER.from_species(s)], idx[ER.to_species(s)]
            M[b, a] += rf; M[a, a] -= rf; M[a, b] += rr; M[b, b] -= rr
            push!(rates, (s, a, b, rf, rr))
        end
        M[n, :] .= 1
        rhs = zeros(BigFloat, n); rhs[n] = 1
        x = M \ rhs
        first_sub = ER.name(first(ER.substrates(ER.reaction(mech))))
        sum(((count(m -> ER.name(m) == first_sub, ER.consumed(s)) -
              count(m -> ER.name(m) == first_sub, ER.released(s))) * (rf * x[a] - rr * x[b])
             for (s, a, b, rf, rr) in rates); init = big(0))
    end
end

function _testhelper_check_against_mass_action(em; n_points = 3, seed = 1)
    rng = Random.MersenneTwister(seed)
    fitted = ER.fitted_params(em)
    mets = ER.metabolites(em)
    for _ in 1:n_points
        point = _testhelper_consistent_point(em, rng)
        concs = Dict{Symbol, BigFloat}(x => exp(3 * big(rand(rng)) - 1.5) for x in mets)
        params = NamedTuple{(fitted..., :Keq, :E_total)}(
            (Float64.(getindex.(Ref(point.vals), fitted))..., Float64(point.Keq), 1.0))
        concs_nt = NamedTuple{Tuple(mets)}(Tuple(Float64(concs[x]) for x in mets))
        v_pkg = ER.rate_equation(em, concs_nt, params)
        v_ref = _testhelper_mass_action_rate(em, point, concs)
        @test isapprox(v_pkg, Float64(v_ref); rtol = 1e-8)
    end
end

@testset "fused steps derive the mass-action rate" begin
    cases = [
        # uni-uni: fused release / fused binding, both orientations, both flags
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S <--> E(S); E(S) <--> E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S ⇌ E(S); E(S) <--> E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S <--> E(S); E(S) ⇌ E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S <--> E(S); E + P <--> E(S) end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E + S <--> E(P); E(P) <--> E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E(P) <--> E + S; E + P <--> E(P) end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E(P) ⇌ E + S; E(P) <--> E + P end end),
        @enzyme_mechanism(begin substrates: S; products: P; steps: begin
            E(P) <--> E + S; E(P) ⇌ E + P end end),
        # ordered bi-bi, one central complex (Segel IX-87) and mixes
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A <--> E(A); E(A) + B <--> E(A, B); E(A, B) <--> E(Q) + P
            E(Q) <--> E + Q end end),
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A <--> E(A); E(A) + B <--> E(A, B); E(A, B) ⇌ E(Q) + P
            E(Q) <--> E + Q end end),
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A ⇌ E(A); E(P, Q) <--> E(A) + B; E(P, Q) <--> E(Q) + P
            E(Q) ⇌ E + Q end end),
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A ⇌ E(A); E(P, Q) <--> E(A) + B; E(P, Q) ⇌ E(Q) + P
            E(Q) ⇌ E + Q end end),
        # ping-pong with an RE fused release
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A <--> E(A); E(A) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B <--> E(B; residual = A - P)
            E(B; residual = A - P) <--> E + Q end end),
        # merged random bi-bi: RE bindings and releases (former missing cut),
        # and all SS with a reversed fused binding (former mixed cut)
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A ⇌ E(A); E + B ⇌ E(B); E(A) + B ⇌ E(A, B); E(B) + A ⇌ E(A, B)
            E(A, B) <--> E(Q) + P; E(A, B) <--> E(P) + Q; E(Q) ⇌ E + Q; E(P) ⇌ E + P
        end end),
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; steps: begin
            E + A <--> E(A); E + B <--> E(B); E(A) + B <--> E(P, Q)
            E(P, Q) <--> E(B) + A; E(P, Q) <--> E(Q) + P; E(P, Q) <--> E(P) + Q
            E(Q) <--> E + Q; E(P) <--> E + P end end),
        # dead ends on a merged complex
        @enzyme_mechanism(begin substrates: A, B; products: P, Q; regulators: I; steps: begin
            E + A ⇌ E(A); E(A) + B ⇌ E(A, B); E(A, B) <--> E(Q) + P; E(Q) ⇌ E + Q
            E + I ⇌ E(I); E(Q) + I ⇌ E(I, Q) end end),
    ]
    for em in cases
        _testhelper_check_against_mass_action(em)
    end
end
```

(The DSL's `steps: begin … end` block takes one step per line; if `;` separators are
not accepted inside the macro, put each step on its own line. Check `test/test_dsl.jl`
for the accepted `regulators:` label and adjust the headers accordingly.)

- [ ] **Step 2: Run and confirm the expected failures**

Focused run of `test/test_step_stoichiometry.jl`. Expected: the fused suite fails on
the reversed fused-step cases (the package returns −v or a wrong law), and the merged
random case with RE bindings and releases errors with "no rapid-equilibrium-consistent
reaction cut". Record which cases fail.

- [ ] **Step 3: Replace the numerator**

In `src/rate_eq_derivation.jl`:

```julia
"""
Numerator of the rate: v·den summed over steady-state steps. With u(f) the
first substrate's exponent in form f's RE weight, each SS step e contributes
ω_e·(forward − reverse flux) with ω_e = (copies of the first substrate e
consumes − copies it releases) + u(from_e) − u(to_e). Flux conservation at every
form makes the sum equal the net consumption of the first substrate for every
parameter value, so it needs no choice of reaction cut and does not depend on
the direction a step is written in. Steps with ω_e = 0 contribute nothing.
"""
function _compute_numerator(mech::Mechanism, enz_name_to_form, step_params,
                            alpha, form_to_group, D)
    x = name(first(substrates(reaction(mech))))
    expo(p) = (mono = only(keys(p));
               k = findfirst(q -> q.first == x, mono);
               k === nothing ? 0 : mono[k].second)
    num = poly_zero()
    for (idx, (s, _)) in enumerate(_flat_steps(mech))
        is_equilibrium(s) && continue
        i_form = enz_name_to_form[name(from_species(s))]
        j_form = enz_name_to_form[name(to_species(s))]
        ω = count(m -> name(m) == x, consumed(s)) - count(m -> name(m) == x, released(s)) +
            expo(alpha[i_form]) - expo(alpha[j_form])
        ω == 0 && continue
        fwd = _ss_contrib(poly_sym(name(step_params[idx][1], mech)),
                          Symbol[name(m) for m in consumed(s)], i_form, alpha)
        rev = _ss_contrib(poly_sym(name(step_params[idx][2], mech)),
                          Symbol[name(m) for m in released(s)], j_form, alpha)
        term = poly_sub(poly_mul(fwd, D[form_to_group[i_form]]),
                        poly_mul(rev, D[form_to_group[j_form]]))
        num = poly_add(num, poly_mul(poly_const(ω), term))
    end
    num
end
```

Delete the old cut search, its docstring, the "ambiguous central-complex cut" error
and the `allow_dead` keyword. In `_raw_symbolic_rate_polys` drop the parameters
`subs_species, prods_species` and the `allow_dead` keyword and call
`_compute_numerator(mech, enz_name_to_form, step_params, alpha, form_to_group, D)`;
update both callers (`_raw_symbolic_rate_polys(M::Type{<:EnzymeMechanism})` ~403 and
`_state_rate_polys` ~1301) and the docstrings that mention `allow_dead` (~1158, ~1692):
a dead inactive state now yields the zero polynomial because its fluxes cancel exactly.

- [ ] **Step 4: Detect an all-RE catalytic cycle while building the weights**

At the end of `_compute_alpha`, check every RE step against the weights: the
concentration exponents of `alpha[to] / alpha[from]` must equal the step's
(consumed − released) counts. A mismatch means an RE cycle performs turnover:

```julia
    for (idx, (s, _)) in enumerate(flat)
        is_equilibrium(s) || continue
        a = enz_name_to_form[name(from_species(s))]
        b = enz_name_to_form[name(to_species(s))]
        ratio = poly_mul(alpha[b], _invert_monomial(alpha[a]))
        step = _re_weight_ratio(s, step_to_K[idx])
        conc(p) = Dict(k => v for (k, v) in only(keys(p)) if k in conc_set)
        conc(ratio) == conc(step) || error(
            "rate_equation: the rapid-equilibrium steps close a catalytic cycle " *
            "(through $(name(from_species(s))) ⇌ $(name(to_species(s)))), so the " *
            "mechanism has no finite rate. Make one step of the cycle steady-state.")
    end
```

`conc_set` is `_concentration_symbols(mech)`; pass it in or compute it inside
`_compute_alpha`.

- [ ] **Step 5: Update the ambiguous-cut test**

In `test/test_rate_eq_derivation.jl`, the testset "Numerator: ambiguous central cut
(regulator sibling) raises": keep the mechanism, rename the testset to
"All-RE catalytic cycle raises", replace its comment with one saying the activator
route E(R) + S ⇌ E(S, R) ⇌ E(P, R) ⇌ E(R) + P is an all-RE catalytic cycle, and change
the assertion to `@test occursin("no finite rate", err.msg)`.

- [ ] **Step 6: Run and confirm**

Focused: `test/test_step_stoichiometry.jl` (every fused case passes),
`test/test_rate_eq_derivation.jl`, `test/test_allosteric_golden.jl`,
`test/test_allosteric_collapse.jl`, `test/allosteric_ground_truth.jl`. Then the full
suite in the background. Expected: 0 failures, 0 errors.

- [ ] **Step 7: Commit**

```bash
git add src/rate_eq_derivation.jl test/test_step_stoichiometry.jl test/test_rate_eq_derivation.jl
git commit   # "Derive the rate numerator without a reaction cut" + trailers
```

---

### Task 3: Theorell–Chance and multi-metabolite steps end to end

**Files:**
- Modify: `src/types.jl` (`_canonical_iso_direction` → `_canonical_step_direction`; `_canonicalize_iso_groups` → `_canonicalize_step_directions`; new `_assert_unique_parameter_names`; both constructors)
- Modify: `src/dsl.jl` (`_split_side`, `_build_step_expr`)
- Modify: `test/test_step_stoichiometry.jl`, `test/test_dsl.jl`

**Interfaces:**
- Consumes: Tasks 1–2.
- Produces: `_canonical_step_direction(s, subs, prods, pure_binding_steps)::Step`,
  `_canonicalize_step_directions(reaction, groups)`, `_assert_unique_parameter_names(steps::Vector{Vector{Step}})`.

- [ ] **Step 1: Write the failing tests**

Append to `test/test_step_stoichiometry.jl`:

```julia
@testset "Theorell–Chance and two-metabolite steps" begin
    tc_ss = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end
    tc_re_outer = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(Q) + P
            E(Q) ⇌ E + Q
        end
    end
    tc_re_step = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B ⇌ E(Q) + P
            E(Q) <--> E + Q
        end
    end
    two_in = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A + B <--> E(A, B)
            E(A, B) <--> E(P, Q)
            E(P, Q) <--> E + P + Q
        end
    end
    @test length(ER.fitted_params(tc_ss)) == 5
    @test length(ER.fitted_params(tc_re_outer)) == 3
    for em in (tc_ss, tc_re_outer, tc_re_step, two_in)
        _testhelper_check_against_mass_action(em)
    end
    # canonical orientation: the substrate side is `from`
    s = only(s for g in ER.steps(ER.Mechanism(tc_ss)) for s in g if !ER.is_iso(s) && !ER.is_binding(s))
    @test ER.name(ER.from_species(s)) == :EA && ER.consumed(s) == ER.Metabolite[ER.Substrate(:B)]
    @test :k_EA_to_EQ in ER.parameters(tc_ss, ER.Full)
    @test :k_EQ_to_EA in ER.parameters(tc_ss, ER.Full)
end

@testset "reversing written steps changes nothing" begin
    rev(s) = ER.Step(ER.to_species(s), ER.from_species(s), ER.released(s), ER.consumed(s),
                     ER.is_equilibrium(s))
    rng = Random.MersenneTwister(7)
    ems = Any[spec.mechanism for spec in MECHANISM_TEST_SPECS
              if spec.mechanism isa ER.EnzymeMechanism]
    for em in ems
        m = ER.Mechanism(em)
        flipped = [[rand(rng, Bool) ? rev(s) : s for s in g] for g in ER.steps(m)]
        m2 = ER.Mechanism(ER.reaction(m), flipped)
        @test m2 == m
        @test ER.rate_equation_string(ER.compile_mechanism(m2)) == ER.rate_equation_string(em)
    end
end

@testset "parameter-name collisions are rejected" begin
    A, B, P = ER.Substrate(:A), ER.Substrate(:B), ER.Product(:P)
    E, EA, EstarA, EAB, EP = _sp([]), _sp([A]), _sp([A], :Estar), _sp([A, B]), _sp([P])
    rxn = @enzyme_reaction(begin
        substrates: A[C]
        products: P[C]
    end)
    # A binds E into two different forms from two groups: both would be kon_A_E.
    err = try
        ER.Mechanism(rxn, [
            [ER.Step(E, EA, [A], ER.Metabolite[], false)],
            [ER.Step(E, EstarA, [A], ER.Metabolite[], false)],
            [ER.Step(EA, EP, ER.Metabolite[], ER.Metabolite[], false)],
            [ER.Step(E, EP, [P], ER.Metabolite[], false)]])
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("same parameter names", err.msg)
    # The same isomerization in two groups: both would be k_EA_to_EP.
    @test_throws ErrorException ER.Mechanism(rxn, [
        [ER.Step(E, EA, [A], ER.Metabolite[], false)],
        [ER.Step(EA, EP, ER.Metabolite[], ER.Metabolite[], false)],
        [ER.Step(EP, EA, ER.Metabolite[], ER.Metabolite[], false)],
        [ER.Step(E, EP, [P], ER.Metabolite[], false)]])
end
```

In `test/test_dsl.jl` add a testset asserting that `E(A) + B <--> E(Q) + P` yields a
step with `consumed == [Substrate(:B)]`, `released == [Product(:P)]`, and that
`E + A + B <--> E(A, B)` yields `consumed == [Substrate(:A), Substrate(:B)]`; and that a
side with two enzyme forms still errors with "more than one enzyme-form term".

- [ ] **Step 2: Run and confirm failure**

Expected: the DSL drops or rejects the second metabolite; the orientation and
collision tests fail.

- [ ] **Step 3: DSL sides with several metabolites**

`_split_side` returns `(enzyme_term, met_terms::Vector{_StepSideTerm})`, collecting
every metabolite term and keeping the "no enzyme-form term" and "more than one
enzyme-form term" errors. `_build_step_expr` builds
`EnzymeRates.Metabolite[$(lhs met exprs...)]` and `…[$(rhs met exprs...)]`. Update both
docstrings and remove the "each elementary step binds at most one metabolite" error.

- [ ] **Step 4: Orient every non-binding step**

Rename `_canonical_iso_direction` to `_canonical_step_direction` and
`_canonicalize_iso_groups` to `_canonicalize_step_directions` (update both callers in
the `Mechanism` and `AllostericMechanism` constructors and all comments that name
them). New body:

```julia
function _canonical_step_direction(s::Step, subs::Set{Symbol}, prods::Set{Symbol},
                                   binding_steps::Vector{Step})
    is_binding(s) && return s
    f, t = from_species(s), to_species(s)
    flip() = Step(t, f, released(s), consumed(s), is_equilibrium(s))
    # Tier 1: substrate/product progression over bound plus free metabolites.
    score(sp, free) = (count(m -> name(m) in subs, bound(sp)) + count(m -> name(m) in subs, free),
                       -count(m -> name(m) in prods, bound(sp)) - count(m -> name(m) in prods, free))
    sf, st = score(f, consumed(s)), score(t, released(s))
    sf > st && return s
    sf < st && return flip()
    # Tier 2: 1-hop binding context.
    fk = _entry_kind(f, binding_steps, subs, prods)
    tk = _entry_kind(t, binding_steps, subs, prods)
    fk == :product_only   && tk == :substrate_only && return s
    fk == :substrate_only && tk == :product_only   && return flip()
    # Tier 3: lex fallback.
    string(name(f)) ≤ string(name(t)) ? s : flip()
end
```

Keep the existing comments on each tier (they stay true).

- [ ] **Step 5: Reject parameter-name collisions**

```julia
"""
Error when two kinetic groups would render the same parameter name. A pure
binding's constants are named by (metabolite, free form), every other step's by
its two forms, so two groups collide when they bind the same metabolite to the
same form or join the same two forms.
"""
function _assert_unique_parameter_names(steps::Vector{Vector{Step}})
    seen = Dict{Tuple, Tuple{Int, Step}}()
    for (g, group) in enumerate(steps), s in group
        m = ligand(s)
        key = m === nothing ?
            (:pair, minmax(String(name(from_species(s))), String(name(to_species(s))))...) :
            (:binding, name(m), m isa CompetitiveInhibitor, name(from_species(s)))
        if haskey(seen, key) && first(seen[key]) != g
            other = last(seen[key])
            error("Mechanism: steps $(name(from_species(other))) → $(name(to_species(other))) " *
                  "and $(name(from_species(s))) → $(name(to_species(s))) are in different " *
                  "kinetic groups but would get the same parameter names")
        end
        seen[key] = (g, s)
    end
end
```

Call it in the `Mechanism` constructor after `_assert_no_re_ss_duplicate(steps)` and at
the corresponding point of the `AllostericMechanism` constructor.

- [ ] **Step 6: Run the focused files, then the full suite**

`test/test_step_stoichiometry.jl`, `test/test_dsl.jl`, `test/test_types.jl`,
`test/test_mechanism_enumeration.jl`, `test/test_rate_eq_derivation.jl`; then the full
suite in the background. Expected: 0 failures, 0 errors. If the collision guard fires on
an existing fixture or an enumerated mechanism, stop and report it: it means two
constants were silently tied before.

- [ ] **Step 7: Commit**

```bash
git add src/types.jl src/dsl.jl test/test_step_stoichiometry.jl test/test_dsl.jl
git commit   # "Write and derive Theorell–Chance and multi-metabolite steps" + trailers
```

---

### Task 4: Behaviour preservation, docs and version

**Files:**
- Run: `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/8633b908-eb4b-4a04-bd9a-bfa7049838d2/scratchpad/ident/regress/regress.jl` (outside the repository; the baseline `snap_base_*` was recorded on the unmodified code)
- Modify: `docs/src/developer.md` (~55, the `Step` field list), any docstring or docs page that mentions `bound_metabolite`, `_step_sides` or "at most one metabolite" (`grep -rn` over `src docs/src`)
- Modify: `Project.toml` (version `0.7.0` → `0.8.0`: the `Step` constructor and fused-step parameter names change)
- Create: `docs/superpowers/specs/2026-09-28-step-explicit-stoichiometry-regression.md` (the comparison report)

**Interfaces:**
- Consumes: everything above.

- [ ] **Step 1: Snapshot the new code and compare**

```bash
R=/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/8633b908-eb4b-4a04-bd9a-bfa7049838d2/scratchpad/ident/regress
nohup julia --project $R/regress.jl snap new base > $R/snap_new.log 2>&1 &
# poll: for i in $(seq 1 240); do sleep 30; grep -q "specs done" $R/snap_new.log && break; done
julia --project $R/regress.jl compare base new
```

Expected: every reaction shows the same per-level counts, no missing mechanisms, and
zero differences in `fitted`, `reduced`, `full`, `error` and `okey`; `specs` differs
only for fixtures with fused steps (renamed parameters). Median time within 20% of
the baseline. For any other difference: evaluate the old and new laws against
`_testhelper_mass_action_rate` (load `test/test_step_stoichiometry.jl` helpers in a
REPL) and record whether the new law is correct; stop and report if the new law is
wrong.

- [ ] **Step 2: Write the regression report**

Create `docs/superpowers/specs/2026-09-28-step-explicit-stoichiometry-regression.md`
with: the populations and counts compared, the result of each comparison, the list of
fixtures whose names changed, timings, and any mechanism whose law changed with the
oracle verdict.

- [ ] **Step 3: Docs and version**

Update `docs/src/developer.md` to describe `Step`'s fields (`from_species`,
`to_species`, `consumed`, `released`, `is_equilibrium`) and that the constructors
canonicalize orientation. Fix every other mention found by
`grep -rn "bound_metabolite\|_step_sides\|at most one metabolite" src docs/src README.md`.
Bump `Project.toml` to `version = "0.8.0"`. If `docs/Project.toml` can be instantiated,
run the doctests: `julia --project=docs -e 'using Documenter, EnzymeRates; doctest(EnzymeRates)'`.

- [ ] **Step 4: Full suite**

Run the full suite in the background. Expected: 0 failures, 0 errors, and
`test_rate_equation_performance` green.

- [ ] **Step 5: Commit**

```bash
git add docs/src/developer.md Project.toml docs/superpowers/specs/2026-09-28-step-explicit-stoichiometry-regression.md
git commit   # "Record the Step refactor's regression check; bump to 0.8.0" + trailers
```
