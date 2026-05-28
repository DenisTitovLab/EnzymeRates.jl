# Structural Parameter Names Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace index-based parameter symbols (`:K1`, `:k6f`, `:K_G6P_reg1`) with structural, direction-independent names (`:K_ATP_E`, `:k_ES_to_EP`, `:K_G6Preg`), delete `Step.source_idx` and the index-context naming path, canonicalize SS steps, and rename the allosteric R/T notion to A/I.

**Architecture:** Names become a pure function of a `Step`'s species + bound metabolite + allosteric state, rendered through the existing single chokepoint (`_param_symbol` / `name(p, m)`). The kinetic-group naming representative and the Haldane elimination pivot share one extracted `_step_priority`. The refactor proceeds in build-green phases; numerical physics oracles survive the rename via a permanent positional-remap test helper, while golden output strings are regenerated mechanically.

**Tech Stack:** Julia, `@generated` rate-equation derivation, Test stdlib, Aqua, JET.

**Spec:** `docs/superpowers/specs/2026-05-27-structural-parameter-names-design.md`

---

## Conventions used in every phase

- **Full suite:** `julia --project -e 'using Pkg; Pkg.test()'` (cold, slow — pays precompile + JIT).
- **Faster single-file iteration:** start a session and include one file, e.g.
  `julia --project -e 'using EnzymeRates; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_rate_eq_derivation.jl")'`
  (some test files `include` the shared mechanism-definitions file; check the file head).
- **Golden regeneration is a capture step, not a guess.** Where a step says "regenerate goldens," the literal new strings are obtained by running the failing test and pasting the *actual* value the test prints into the `expected_*` field. These cannot be pre-written because they are outputs of the code being changed. Each such step names the exact command and the file to edit.
- **Commit after every task.** End commit messages with the `Co-Authored-By` trailer per the repo convention.
- Branch: `refactor-to-concrete-types-instead-of-symbols` (already current).

## File-by-file responsibility map

| File | Change |
|---|---|
| `src/types.jl` | `name(Species)` concat; chokepoint `_param_symbol`/`name` structural rewrite; `_inactive_name` helper; delete `name(::Type{P},idx)`, `_rep_idx_for_step`; `Step.source_idx` removal; `Mechanism`/`AllostericMechanism` constructor cleanup (incl. unified RE+SS physical-forward iso canonicalization, Phase 5); remove the RE-iso lex branch from the `Step` constructor; `RegulatorySite` A/I validator |
| `src/rate_eq_derivation.jl` | A/I branch-state + helper renames (`_onlyR_*`→`_onlyA_*`, `_T_rename*`→`_I_rename*`); `@generated` callers switch to value-context `name(p,m)`; `_step_priority` consumers |
| `src/thermodynamic_constr_for_rate_eq_derivation.jl` | extract `_step_priority`; rep = `argmin`; drop `source_idx`-keyed ordering in `_step_parameters` |
| `src/mechanism_enumeration.jl` | SS canonicalization fallout in dedup; simplify canonical-hash token layer; A/I in `_T_rename_parameters` callers |
| `src/dsl.jl` | A/I state-annotation parsing |
| `src/identify_rate_equation.jl` | canonical-hash simplification (cleanup phase) |
| `.claude/CLAUDE.md` | A/I taxonomy + naming-convention docs |
| `test/test_rate_eq_derivation.jl` | oracle positional-shim wiring; golden regen |
| `test/mechanism_definitions_for_test_enzyme_derivation.jl` | golden regen; DSL A/I annotations |
| `test/test_types.jl`, `test/test_accessors.jl` | form-name + param-name golden regen |
| `test/test_dsl.jl`, `test/test_mechanism_enumeration.jl`, `test/test_identify_rate_equation.jl` | golden + name-literal regen |
| `test/test_chokepoint.jl` | drop allowance for deleted index-context entry point |

---

## Phase 0 — Oracle positional-shim scaffolding

Goal: let numerical oracles keep destructuring `k1f, k2f, …` after the rename. The shim is a no-op today (oracle mechanisms are one-step-per-group, so consecutive index == current rep index), which Phase 0 verifies.

### Task 0.1: Add `positional_params` helper and route oracles through it

**Files:**
- Modify: `test/test_rate_eq_derivation.jl` (add helper near line 102; edit `test_analytical_rate` at line 503-508)

- [ ] **Step 1: Write the helper**

In `test/test_rate_eq_derivation.jl`, just above `random_reduced_params` (line 104):

```julia
"""
Re-key a structural-named parameter NamedTuple to the **per-step** positional
names (`K1`, `k1f`, `k1r`, …) that hand-derived analytical oracles destructure.
Walks every step in source order; emits one positional name per RE step or two
per SS step keyed on the step's source position; each value is looked up by the
rep parameter's structural name in `nt` (so group members share the rep's
value, matching how the package's destructure-by-name works today).
Permanent test utility — oracles are inherently positional and source-indexed.
"""
function positional_params(m, nt::NamedTuple)
    mech = m isa EnzymeRates.Mechanism ? m : EnzymeRates.Mechanism(m)
    names = Symbol[]
    vals  = Any[]
    for group in EnzymeRates.steps(mech)
        rep = first(group)
        if EnzymeRates.is_equilibrium(rep)
            rep_name = EnzymeRates.name(EnzymeRates.Kd(rep, :None), mech)
            for s in group
                push!(names, Symbol("K", EnzymeRates.source_idx(s)))
                push!(vals,  nt[rep_name])
            end
        else
            fwd_name = EnzymeRates.name(EnzymeRates.Kfor(rep, :None), mech)
            rev_name = EnzymeRates.name(EnzymeRates.Krev(rep, :None), mech)
            for s in group
                idx = EnzymeRates.source_idx(s)
                push!(names, Symbol("k", idx, "f")); push!(vals, nt[fwd_name])
                push!(names, Symbol("k", idx, "r")); push!(vals, nt[rev_name])
            end
        end
    end
    NamedTuple{Tuple(names)}(Tuple(vals))
end
```

> **Note: shim uses `source_idx` of each step.** This matches today's per-step naming convention exactly (RE iso steps `:K2` etc. → `:K2` value = rep's K-value; multi-step kinetic groups share values via the lookup). After Phase 6 deletes `source_idx`, replace the lookup with the step's position in the flat iteration order (which equals `source_idx` today by Task 6.1's precondition assertion).

- [ ] **Step 2: Route the oracle call through the shim**

In `test_analytical_rate` (line 507), change:

```julia
        p = merge(all_params, (Et=Et,))
```
to:
```julia
        p = merge(positional_params(m, all_params), (Et=Et,))
```

- [ ] **Step 3: Run analytical tests to verify no-op (still green)**

Run: `julia --project -e 'using EnzymeRates, Random, Test; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_rate_eq_derivation.jl")'`
Expected: PASS (the shim is currently an identity remap for oracle mechanisms). If any analytical mechanism is multi-step-per-group, the `@assert` or a `KeyError` fires here — that mechanism's oracle numbering must be reconciled before proceeding.

- [ ] **Step 4: Commit**

```bash
git add test/test_rate_eq_derivation.jl
git commit -m "test: route analytical oracles through positional_params shim"
```

---

## Phase 1 — `name(Species)` full concat

Goal: form names render without internal underscores (`:E_ATP`→`:EATP`). Only form-name goldens change; parameter names (still `:K1`) are unaffected.

### Task 1.1: Concatenate conformation + bound metabolites in `name(Species)`

**Files:**
- Modify: `src/types.jl:77-93`
- Test: form-name assertions in `test/test_types.jl`, `test/test_accessors.jl`

- [ ] **Step 1: Change the renderer**

In `src/types.jl`, replace the body of `name(s::Species)` (lines 77-92). The first bound metabolite is glued to the conformation, and subsequent bounds are glued to each other; the residual block keeps its `_res`/`+`/`-` markers:

```julia
function name(s::Species)
    head = String(conformation(s))
    for m in bound(s)
        head *= m isa CompetitiveInhibitor ?
                String(name(m)) * "inh" : String(name(m))
    end
    parts = String[head]
    if has_residual(s)
        push!(parts, "res")
        for a in added(residual(s))
            push!(parts, "+" * String(name(a)))
        end
        for r in subtracted(residual(s))
            push!(parts, "-" * String(name(r)))
        end
    end
    Symbol(join(parts, "_"))
end
```

Update the comment at lines 71-76 to describe the concat form (`:E`, `:ES`, `:EATP`, `:EstarA_res_+P`) and keep the "metabolite Symbols must not contain `_`" domain note.

- [ ] **Step 2: Run form-name tests to see the new actuals**

Run: `julia --project -e 'using EnzymeRates, Test; include("test/test_types.jl")'`
Expected: FAIL on form-name assertions (e.g. expecting `:E_ATP`, got `:EATP`).

- [ ] **Step 3: Regenerate form-name goldens**

Edit `test/test_types.jl` and `test/test_accessors.jl`: replace each underscored form-name literal (`:E_S`, `:E_A_B`, …) with its concat form (`:ES`, `:EAB`, …) as printed by the failing run. Grep aid: `grep -nE ':E[a-z]*_[A-Z]' test/test_types.jl test/test_accessors.jl`.

- [ ] **Step 4: Run to verify pass**

Run: `julia --project -e 'using EnzymeRates, Test; include("test/test_types.jl"); include("test/test_accessors.jl")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl test/test_accessors.jl
git commit -m "refactor: name(Species) concatenates conformation and bound metabolites"
```

### Task 1.2: Sweep remaining form-name goldens

**Files:** `test/mechanism_definitions_for_test_enzyme_derivation.jl`, `test/test_dsl.jl`, `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Run the full suite to surface every remaining form-name golden**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL only on form-name string assertions (enzyme_forms, stoich labels, `rate_equation_string` enzyme-form references).

- [ ] **Step 2: Regenerate each failing form-name golden from the printed actual.**

- [ ] **Step 3: Re-run full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "test: regenerate form-name goldens for concat Species names"
```

---

## Phase 2 — A/I rename (internal, render-preserving)

Goal: rename the allosteric R/T notion to A/I across taxonomy symbols, branch-state values, helper names, the DSL, and CLAUDE.md — **without changing any rendered parameter symbol** (the chokepoint keeps emitting today's format: `:None`/`:A`→plain, `:I`→`_T` suffix, until Phase 3). Only DSL state annotations in test mechanism definitions change.

### Task 2.1: Rename taxonomy symbols in the `RegulatorySite` validator and `AllostericMechanism`

**Files:** `src/types.jl:109-113` (validator), `AllostericMechanism` validation, `src/dsl.jl` (annotation parser)

- [ ] **Step 1: Write/adjust a failing test** in `test/test_types.jl`: assert `RegulatorySite([AllostericRegulator(:G6P)], 1, [:OnlyA])` constructs and `[:OnlyR]` errors.

```julia
@test RegulatorySite([AllostericRegulator(:G6P)], 1, [:EqualAI]) isa RegulatorySite
@test_throws ErrorException RegulatorySite([AllostericRegulator(:G6P)], 1, [:OnlyR])
```

- [ ] **Step 2: Run to confirm fail** — Run the file; expected FAIL (`:OnlyA` rejected today).

- [ ] **Step 3: Update the validator** at `src/types.jl:110-112`:

```julia
            st in (:OnlyA, :OnlyI, :EqualAI, :NonequalAI) ||
                error("RegulatorySite: allo state $st must be one of " *
                      ":OnlyA, :OnlyI, :EqualAI, :NonequalAI")
```
Apply the same symbol set to the `AllostericMechanism` constructor's allo-state validation (search `:NonequalRT` in `src/types.jl`).

- [ ] **Step 4: Run to confirm pass.**

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "refactor: allosteric-state taxonomy R/T -> A/I in validators"
```

### Task 2.2: Rename taxonomy + branch states throughout `rate_eq_derivation.jl` and `mechanism_enumeration.jl`

**Files:** `src/rate_eq_derivation.jl`, `src/mechanism_enumeration.jl`, `src/sym_poly_for_rate_eq_derivation.jl`

- [ ] **Step 1: Mechanical symbol rename.** Apply these substitutions across the three files (literal symbol tokens and identifier names):
  - taxonomy: `:OnlyR`→`:OnlyA`, `:OnlyT`→`:OnlyI`, `:EqualRT`→`:EqualAI`, `:NonequalRT`→`:NonequalAI`
  - branch-state values in `Parameter` construction and rename maps: `, :R)`→`, :A)`, `, :T)`→`, :I)`, and `=== :T`→`=== :I`, `=== :R`→`=== :A`, `? :R`→`? :A`, `: :T`→`: :I` (review each hit — do not touch unrelated `:R`/`:T` such as matrix `R` variables)
  - helper identifiers: `_onlyR_syms`→`_onlyA_syms`, `_onlyR_parameters`→`_onlyA_parameters`, `_T_rename`→`_I_rename`, `_T_rename_parameters`→`_I_rename_parameters`, and any `K_R`/`K_T` local names → `K_A`/`K_I`.

  Grep to enumerate hits first: `grep -nE ':OnlyR|:OnlyT|:EqualRT|:NonequalRT|_onlyR|_T_rename|, :R\)|, :T\)|=== :[RT]\b' src/rate_eq_derivation.jl src/mechanism_enumeration.jl src/sym_poly_for_rate_eq_derivation.jl`

- [ ] **Step 2: Keep the chokepoint render unchanged for now.** In `src/types.jl`, the index-context `_param_symbol(::Type{P}, idx, state)` (line 1407) currently keys on `state === :T`. Update it to `state === :I` so the `_T` suffix still renders for inactive branches (render output byte-identical):

```julia
_param_symbol(::Type{P}, idx::Int, state::Symbol) where {P<:Parameter} =
    state === :I ? Symbol(_param_symbol(P, idx), "_T") :
                   _param_symbol(P, idx)
```
And `_param_symbol(::Type{Kreg}, …)` (line 1413) keep `_T_reg` for `state === :I`. (These strings die in Phase 3; here we only keep output stable.)

- [ ] **Step 3: Update DSL annotations in test mechanism definitions.** In `test/mechanism_definitions_for_test_enzyme_derivation.jl` and `test/test_dsl.jl`, change user-written allosteric state annotations `OnlyR/OnlyT/EqualRT/NonequalRT` → `OnlyA/OnlyI/EqualAI/NonequalAI`. Grep: `grep -rnE 'OnlyR|OnlyT|EqualRT|NonequalRT' test/`.

- [ ] **Step 4: Run full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS — rendered parameter symbols are unchanged, so no golden churn beyond the DSL annotation edits.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "refactor: rename allosteric R/T to A/I (internal, render-preserving)"
```

### Task 2.3: Update CLAUDE.md allosteric-state taxonomy docs

**Files:** `.claude/CLAUDE.md`

- [ ] **Step 1:** Replace the `:OnlyR/:OnlyT/:EqualRT/:NonequalRT` taxonomy descriptions with `:OnlyA` (active-only), `:OnlyI` (inactive-only), `:EqualAI`, `:NonequalAI`, and update prose ("R-state"→"active state", "T-state"→"inactive state"). Keep the semantics identical.

- [ ] **Step 2: Commit**

```bash
git add .claude/CLAUDE.md
git commit -m "docs: CLAUDE.md allosteric taxonomy R/T -> A/I"
```

---

## Phase 3 — Structural chokepoint rewrite

Goal: `_param_symbol`/`name(p, m)` emit structural names. Delete the index-context companion and `_rep_idx_for_step`; switch `@generated` callers to value-context. Rep stays `first(group)`. Regenerate all parameter-name goldens.

### Task 3.1: Rewrite the chokepoint to render structural names

**Files:** `src/types.jl:1396-1449`

- [ ] **Step 1: Write a failing unit test** in `test/test_types.jl` for a known mechanism:

```julia
@testset "structural parameter names" begin
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E(P)
            E(P) <--> E + P
        end
    end
    ps = EnzymeRates.parameters(m)
    @test :k_ES_to_EP in ps        # SS catalytic forward, structural
    @test !any(p -> occursin(r"^k[0-9]", String(p)), ps)  # no index names
end
```

- [ ] **Step 2: Run to confirm fail** — expected `:k1f`-style names present, `:k_ES_to_EP` absent.

- [ ] **Step 3: Replace the chokepoint bodies** at `src/types.jl:1396-1449`. The value-context `name(p, m)` must still resolve a step to its kinetic-group **representative** (group members share one parameter), so `_rep_idx_for_step` is *replaced* by `_rep_step` (returns the rep `Step` instead of an index) rather than deleted. The render is then a pure function of the rep step + state:

```julia
# Structural parameter-name rendering. Every Parameter → Symbol passes
# through here. The name is a pure function of the kinetic-group rep
# Step's species pair + bound metabolite + allosteric branch state.
#
# State token (placed right after the type prefix): active branch → "A_",
# inactive branch → "I_", shared/none → omitted.
_state_tag(state::Symbol) = state === :A ? "A_" :
                            state === :I ? "I_" : ""

# Render a binding param from its rep step: metabolite + pre-binding form.
_render_binding(prefix::String, rep::Step, state::Symbol) =
    Symbol(prefix, _state_tag(state),
           String(name(bound_metabolite(rep))), "_",
           String(name(from_species(rep))))

# Render an iso param by directed species pair (each rate constant names
# its own direction; storage order is canonical after Phase 5).
_render_iso(prefix::String, from::Species, to::Species, state::Symbol) =
    Symbol(prefix, _state_tag(state),
           String(name(from)), "_to_", String(name(to)))

# Find the kinetic group containing `step`; return its naming rep.
# Phase 3: rep = first(group). Phase 4 routes this through _group_rep.
function _rep_step(step::Step, m::Union{Mechanism,AllostericMechanism})
    for group in steps(m)
        step in group && return first(group)
    end
    error("Step not found in mechanism: $step")
end
_rep_step(step::Step, m::EnzymeMechanism) = _rep_step(step, Mechanism(m))
_rep_step(step::Step, m::AllostericEnzymeMechanism) =
    _rep_step(step, AllostericMechanism(m))

const _AnyMech =
    Union{Mechanism,EnzymeMechanism,AllostericMechanism,AllostericEnzymeMechanism}

name(p::Kd,   m::_AnyMech) = _render_binding("K_",    _rep_step(p.step, m), p.state)
name(p::Kon,  m::_AnyMech) = _render_binding("kon_",  _rep_step(p.step, m), p.state)
name(p::Koff, m::_AnyMech) = _render_binding("koff_", _rep_step(p.step, m), p.state)
function name(p::Kiso, m::_AnyMech)
    rep = _rep_step(p.step, m); _render_iso("Kiso_", from_species(rep), to_species(rep), p.state)
end
function name(p::Kfor, m::_AnyMech)
    rep = _rep_step(p.step, m); _render_iso("k_", from_species(rep), to_species(rep), p.state)
end
function name(p::Krev, m::_AnyMech)
    rep = _rep_step(p.step, m); _render_iso("k_", to_species(rep), from_species(rep), p.state)
end

# Regulator-site parameter: ligand name + "reg", no site index (Risk R2
# accepted: same ligand at two sites is not enumerated). No group lookup
# needed — the ligand is carried on the Parameter.
name(p::Kreg, ::Union{AllostericMechanism,AllostericEnzymeMechanism}) =
    Symbol("K_", _state_tag(p.state), String(name(p.ligand)), "reg")

name(::Keq,   _) = :Keq
name(::Etot,  _) = :E_total   # unchanged — see plan rationale (~15 hard-coded sites)
name(::Lallo, _) = :L
```

Delete `_param_symbol` (all methods), `name(::Type{P}, idx)` / `name(::Type{P}, idx, state)`, the old `_rep_idx_for_step`, and `_site_idx_of` if now unused. Delete the index-context comment block at lines 1396-1425.

> **`Kreg` name injectivity (Risk R2).** The body above never inspects `m` and cannot detect a same-ligand-two-sites collision at name time. Per Denis's note ("same ligand is almost never an allo regulator at two sites"), do NOT add an unreachable error branch in `name(p::Kreg, …)`. Instead, add an assertion in the `AllostericMechanism` constructor that no ligand appears in two distinct `RegulatorySite`s. That gives the same safety with sibling visibility and a clear error site.

> Note on `Krev`/`Kiso`: with unified iso canonicalization not yet in place (Phase 5), iso `from`/`to` follow stored order (RE iso = Step ctor's lex; SS iso = source). `Kfor`/`Krev` name the forward and reverse directed transitions of the rep step, so the pair is direction-independent regardless of storage order. After Phase 5 the stored direction becomes physical-forward for both RE and SS iso.

- [ ] **Step 4: Run the unit test to confirm pass.**

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "refactor: chokepoint renders structural parameter names"
```

### Task 3.1.5: Consolidate the 11 `string(k) * "_T"` synth-dep sites through one `_inactive_name` helper

The chokepoint (Task 3.1) now emits the inactive token mid-name (`:K_I_ATP_E`). But ~11 sites synthesize inactive-state names for *dependent* params that have no `Parameter` struct (Wegscheider/Haldane-eliminated symbols) by appending `_T` to the rendered active name. After Task 3.1 those sites would produce `:K_ATP_E_T` while the chokepoint produces `:K_I_ATP_E` for the same parameter — the allosteric rate body would reference dep symbols that don't match the base param names. This task introduces one structural helper and routes every synth-dep site through it (the DRY simplification Denis asked for).

**Files:** `src/rate_eq_derivation.jl` (sites at lines 813, 1206, 1256, 1300, 1303, 1401, 1454, 1473, 1515), `src/mechanism_enumeration.jl` (sites at lines 2000-2002), `src/types.jl` (helper definition next to chokepoint)

- [ ] **Step 1: Write a failing allosteric integration test** in `test/test_rate_eq_derivation.jl`: pick an existing `MECHANISM_TEST_SPECS` entry that has a `NonequalAI` group (T-state distinct params) and a Haldane-derived dep, build params, call `rate_equation`, compare to its analytical oracle through the positional shim. After Task 3.1 alone this test should FAIL because the dep symbol in the generated body is `:K_…_T` while `parameters(m)` advertises `:K_I_…`. Run to confirm FAIL.

- [ ] **Step 2: Add the helper** in `src/types.jl` next to the chokepoint render helpers:

```julia
# Structural inactive-state name from an active-state name. Synth-dep
# machinery (Wegscheider/Haldane elimination) operates on already-rendered
# symbols of NonequalAI-active dep parameters; those carry an `_A_` state
# token after the type prefix (e.g. :K_A_ATP_E). The inactive variant
# **swaps** that token to `_I_` (not appends/inserts) so the result
# matches what the chokepoint emits for the same parameter with state=:I.
#
# Errors if the input has no `_A_` token: an EqualAI / non-allosteric
# parameter has no separate inactive variant, and the synth-dep machinery
# must not invoke this on such symbols.
function _inactive_name(sym::Symbol)
    s = String(sym)
    occursin("_A_", s) ||
        error("_inactive_name: $sym has no _A_ token to swap to _I_ " *
              "(EqualAI / non-allosteric symbols have no inactive variant)")
    Symbol(replace(s, "_A_" => "_I_"; count=1))
end
```

- [ ] **Step 3: Route every `Symbol(string(k) * "_T")` site through `_inactive_name`.**

Mechanical substitution. For each line:

| File:Line | Old | New |
|---|---|---|
| `rate_eq_derivation.jl:813` | `rename_T[k] = Symbol(string(k) * "_T")` | `rename_T[k] = _inactive_name(k)` |
| `rate_eq_derivation.jl:1206` | `push!(names, Symbol(string(k) * "_T"))` | `push!(names, _inactive_name(k))` |
| `rate_eq_derivation.jl:1256` | `rename_T[k] = Symbol(string(k) * "_T")` | `rename_T[k] = _inactive_name(k)` |
| `rate_eq_derivation.jl:1300` | `push!(indep_T_list, Symbol(string(p) * "_T"))` | `push!(indep_T_list, _inactive_name(p))` |
| `rate_eq_derivation.jl:1303` | `dep_T[Symbol(string(p) * "_T")] = p` | `dep_T[_inactive_name(p)] = p` |
| `rate_eq_derivation.jl:1401` | `rename_T[k] = Symbol(string(k) * "_T")` | `rename_T[k] = _inactive_name(k)` |
| `rate_eq_derivation.jl:1454` | `Expr(:(=), Symbol(string(p) * "_T"), p)` | `Expr(:(=), _inactive_name(p), p)` |
| `rate_eq_derivation.jl:1473` | `t_sym = Symbol(string(sym) * "_T")` | `t_sym = _inactive_name(sym)` |
| `rate_eq_derivation.jl:1515` | `rename_T[k] = Symbol(string(k) * "_T")` | `rename_T[k] = _inactive_name(k)` |
| `mechanism_enumeration.jl:2000` | `t_str = r_str * "_T"` | `t_str = String(_inactive_name(Symbol(r_str)))` |
| `mechanism_enumeration.jl:2002` | `name_map[t_str] = tok * "_T"` | `name_map[t_str] = tok * "_T"` (this one stays — `tok` is the canonical `p_$i` token, not a parameter symbol; document why) |

Grep audit after edits: `grep -rn '\* *"_T"' src/` should return only the `tok * "_T"` line at `mechanism_enumeration.jl:2002` (with a comment justifying it) and zero other hits. If anything else remains, route it.

- [ ] **Step 4: Update `test_chokepoint.jl`** to recognize `_inactive_name` as a chokepoint body. Add `:_inactive_name` to the allowed `fn_name` set so its `Symbol(...)` call doesn't trip the AST walker.

- [ ] **Step 5: Run the allosteric integration test + full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS (the failing test from Step 1 now green; allosteric goldens may shift — regenerate as part of Task 3.3).

- [ ] **Step 6: Commit**

```bash
git add src/types.jl src/rate_eq_derivation.jl src/mechanism_enumeration.jl test/test_chokepoint.jl test/test_rate_eq_derivation.jl
git commit -m "refactor: route synth-dep T-names through _inactive_name helper"
```

### Task 3.2: Delete the now-identity Pass-1 rename; convert index-context callers to value-context

The rate-equation polynomial already labels its symbols via value-context
`_raw_param_symbols` (one rep parameter per kinetic group), so it becomes
structural automatically once the chokepoint is structural (Task 3.1). The
only index-context (`name(::Type{P}, idx)`) callers are in
`_build_kinetic_rename_map` and the `rate_equation_string` display block. Since
value-context `name(p, m)` already collapses a group member to its rep's name,
the **Pass-1 kinetic-group rename entries become identity** (`rep => rep`) and
the Pass-1 loop is deleted. Pass-2 (Wegscheider single-symbol ties) stays.

**Files:** `src/rate_eq_derivation.jl` (`_build_kinetic_rename_map` 119-170; `rate_equation_string` display 609-635)

- [ ] **Step 1: Rewrite `_build_kinetic_rename_map`** (`src/rate_eq_derivation.jl:119-170`). Delete the Pass-1 loop (125-138); start `rename` empty; build `binding_set` from value-context rep names; keep Pass-2 verbatim:

```julia
function _build_kinetic_rename_map(M::Type{<:EnzymeMechanism})
    m = M()
    mech = Mechanism(m)
    rename = Dict{Symbol, Symbol}()
    rxns = reactions(m)
    eq = equilibrium_steps(m)
    enz_set = Set(enzyme_forms(m))
    step_params = _step_parameters(mech)
    # binding-K set: value-context rep name of each RE binding group.
    binding_set = Set{Symbol}()
    for (idx, (lhs, _, _, _)) in enumerate(rxns)
        eq[idx] || continue
        any(s ∉ enz_set for s in lhs) || continue
        push!(binding_set, name(step_params[idx][1], mech))
    end
    # Pass 2: single-symbol Wegscheider RE ties between two binding K's.
    dep_raw, _ = _dependent_param_exprs_kernel(M, rename)
    for (lhs, rhs) in dep_raw
        rhs isa Symbol || continue
        lhs in binding_set && rhs in binding_set || continue
        target = get(rename, rhs, rhs)
        rename[lhs] = target
        for k in collect(keys(rename))
            rename[k] == lhs && (rename[k] = target)
        end
    end
    rename
end
```
(`_step_parameters[idx]` aligns with `rxns[idx]` — both flat source order; `name(p, mech)` renders the rep name.)

- [ ] **Step 2: Delete the Pass-1 display loops** in BOTH non-allosteric and allosteric `rate_equation_string` blocks. (a) `src/rate_eq_derivation.jl:617-635` — the `user_lines` loop interpolating `"K$idx = K$rep"`. (b) `src/rate_eq_derivation.jl:1633-1638` — the allosteric Pass-1 display analog. With structural names, tied group members share one symbol, so there is nothing to annotate at either site. Remove both `user_lines` accumulations and any now-empty `# user-defined … equalities:` headers. Keep the Wegscheider Pass-2 annotation block that follows.

- [ ] **Step 3: Confirm no index-context callers remain**

Run: `grep -rnE 'name\((Kd|Kiso|Kon|Koff|Kfor|Krev|Kreg),\s*[a-z_]*idx' src/`
Expected: no matches. (The reviewer's audit pattern; if this leaves a match, the deletion missed a site.)

- [ ] **Step 4: Run derivation + string tests**

Run: `julia --project -e 'using EnzymeRates, Test; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_rate_eq_derivation.jl")'`
Expected: FAIL only on golden strings (names now structural, fewer annotation lines); numerical oracles PASS (shim).

- [ ] **Step 5: Commit the code change** (goldens regenerated in Task 3.3):

```bash
git add src/rate_eq_derivation.jl
git commit -m "refactor: derivation uses value-context names; drop identity Pass-1 rename"
```

### Task 3.3: Regenerate all parameter-name goldens + name literals

**Files:** `test/mechanism_definitions_for_test_enzyme_derivation.jl`, `test/test_rate_eq_derivation.jl`, `test/test_types.jl`, `test/test_mechanism_enumeration.jl`, `test/test_identify_rate_equation.jl`, `test/test_fitting.jl`

- [ ] **Step 1: Run full suite to surface every param-name golden + `parameters()` tuple assertion.**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL on `expected_factored_num`/`expected_factored_denom`, `rate_equation_string` comparisons, and literal `:K1`/`:k6f` `parameters()` assertions.

- [ ] **Step 2: Regenerate each failing golden / name-literal from the printed actual.** Work file by file. For `parameters()` tuple assertions, replace index names with the structural names printed. (Numerical oracle formulas are NOT touched — the shim handles them.)

- [ ] **Step 3: Re-run full suite to PASS.**

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "test: regenerate parameter-name goldens for structural chokepoint"
```

### Task 3.4: Update the chokepoint AST test

**Files:** `test/test_chokepoint.jl`

- [ ] **Step 1:** The walker forbids raw `Symbol("[KkVL]…")` outside chokepoint bodies. The new render helpers `_render_binding(prefix::String, ...)` and `_render_iso(prefix::String, ...)` construct `Symbol("K_"...)`/`Symbol("k_"...)` but their first arg is a `String`, not a Parameter, so the current `_is_chokepoint_def` matcher (which requires the first arg to mention `Parameter|::K[a-z]|::Type{`) will NOT allow them and the test will fail. Relax the matcher to also accept by function name: change the gate so a def is a chokepoint body if `fn_name in (:name, :_param_symbol, :_render_binding, :_render_iso)` OR (`fn_name == :name` AND first-arg mentions a Parameter). Concretely, add `:_render_binding` and `:_render_iso` to the `fn_name in (...)` set at the line that currently reads `fn_name in (:name, :_param_symbol) || return false`.

- [ ] **Step 2: Run** `julia --project -e 'using EnzymeRates, Test; include("test/test_chokepoint.jl")'` — Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add test/test_chokepoint.jl
git commit -m "test: chokepoint walker recognizes structural-name helpers"
```

---

## Phase 4 — `_step_priority` extraction and group-rep = argmin

Goal: select the kinetic-group naming representative by structural primacy (one priority function shared with the Haldane pivot), instead of `first(group)`.

### Task 4.1: Extract `_step_priority(step, m)`

**Files:** `src/thermodynamic_constr_for_rate_eq_derivation.jl:272-300`

- [ ] **Step 1: Write a failing test** in `test/test_rate_eq_derivation.jl` (or `test_types.jl`) asserting that for a mechanism with a dead-end mirror sharing a group with a free-enzyme binding step, the group's parameter name is the free-enzyme binding step's structural name (not the mirror's). Use an existing allosteric/dead-end fixture from the mechanism-definitions file; assert the expected structural name is in `parameters(m)`.

- [ ] **Step 2: Run to confirm fail** (today rep = `first(group)`, which after canonicalization may be the mirror).

- [ ] **Step 3: Extract the scorer.** Pull the priority computation (currently inline at lines 272-300, keyed on reactions) into a pure function scoring a single `Step`:

```julia
"""
Structural primacy base score for a step (lower = more primary / less
eliminable). Free-enzyme RE binding (-1) < free-enzyme SS binding (0) <
non-free metabolite step (10) < internal isomerization (20). Shared by the
kinetic-group name representative (argmin) and the Haldane elimination pivot
(argmax, which adds a +0/+1 forward/reverse offset per rate constant).
"""
function _step_priority(s::Step, free_enz_set::Set{Symbol})
    has_met = is_binding(s)
    is_free = (name(from_species(s)) in free_enz_set) ||
              (name(to_species(s))   in free_enz_set)
    is_equilibrium(s) && has_met && is_free && return -1
    return !has_met ? 20 : is_free ? 0 : 10
end
```
(Compute `free_enz_set` once via the existing logic at lines 215-223.)

In the pivot loop (lines 276-300), replace the inline base computation with `_step_priority(step, free_enz_set)`, keeping the existing per-parameter offset for SS columns: eq steps → `priority = _step_priority(step, free_enz_set)`; SS steps → `priority = _step_priority(step, free_enz_set) + (offset - 1)`. This reproduces today's numeric ranking exactly (pure refactor — Decision 1).

- [ ] **Step 4: Use it for rep selection.** Add a helper returning the naming rep of a group, with a deterministic lexical tiebreak on the step's species-pair string:

```julia
# Total lexical tiebreak: species pair + bound metabolite + RE/SS so any
# two distinct steps in a group break to a deterministic order.
_step_lex_key(s::Step) =
    (String(name(from_species(s))), String(name(to_species(s))),
     String(bound_metabolite(s) === nothing ? "" : name(bound_metabolite(s))),
     is_equilibrium(s))

_group_rep(group::Vector{Step}, free_enz_set::Set{Symbol}) =
    argmin(s -> (_step_priority(s, free_enz_set), _step_lex_key(s)), group)
```

Route **every** existing "rep = first(group)" naming site through `_group_rep`. Reviewer audit identified at least these (the implementer must grep for the exact set and route all of them in one task — any one site left at `first(group)` while another uses `_group_rep` causes the generated body to reference symbols `parameters(m)` doesn't advertise):

| File:Line | Site |
|---|---|
| `src/types.jl:1462` | `_enumerate_parameters_full` (non-allosteric polynomial leaves) |
| `src/rate_eq_derivation.jl:1045` | `_onlyA_parameters` (allosteric R-branch enumeration) |
| `src/rate_eq_derivation.jl:1072` | `_I_rename_parameters` R-side loop |
| `src/rate_eq_derivation.jl:1105` | `_I_rename_parameters` T-side loop |
| `src/rate_eq_derivation.jl:1155` | `_enumerate_parameters_full_allosteric` |
| chokepoint `_rep_step` (Task 3.1) | non-rep parameter rendering |

Grep audit before commit: `grep -rn 'first(group)\|first(idxs)' src/` should show no remaining naming-context uses. `_rep_step` becomes:

```julia
function _rep_step(step::Step, m::Union{Mechanism,AllostericMechanism})
    fes = _free_enz_set(m)
    for group in steps(m)
        step in group && return _group_rep(group, fes)
    end
    error("Step not found in mechanism: $step")
end
```
where `_free_enz_set(m)` is the extracted free-enzyme-set helper (lines 215-223). Put `_step_priority`/`_group_rep`/`_step_lex_key`/`_free_enz_set` in `src/thermodynamic_constr_for_rate_eq_derivation.jl`. `types.jl`'s `_rep_step` may call them even though `types.jl` is `include`d first — Julia resolves plain function calls at call time, so intra-module definition order does not matter (these run at rate-equation build/runtime, by which point all methods exist).

- [ ] **Step 5: Add a Decision-1 dep-set guard keyed on structural Parameter IDENTITY, not rendered names.** Rendered names change for two unrelated reasons across this refactor (Phase 4 rep change → different rep name; Phase 5 iso flip → different from/to). A name-keyed guard conflates those with a real dep-choice drift. The correct invariant for Decision 1 is "the same *Parameter struct* becomes dependent." Snapshot the dep set as `{(Type, hash(step or site), state)}` tuples — invariant under both rep-change and iso direction flip.

```julia
# test/test_dep_set_invariance.jl  (new file)
"""Structural canonical key for a Parameter — Step.hash already ignores
source_idx, so this is stable under the Phase 4/5/6 transformations."""
_dep_key(p::EnzymeRates.Parameter) = _dep_key(p, EnzymeRates.name)  # see below
function _struct_key(name_sym::Symbol, m)
    # Resolve a dep RHS Symbol back to its Parameter struct via
    # `_enumerate_parameters_full(m)` + name lookup; gives a structural
    # tuple even for synth-dep T-symbols (route through _inactive_name).
    ...
end

@testset "Decision 1: dep-set structurally unchanged" begin
    for spec in MECHANISM_TEST_SPECS
        em = EnzymeMechanism(spec.mechanism)
        dep_exprs, indep = EnzymeRates._dependent_param_exprs(typeof(em))
        @test Set(_struct_key(k, spec.mechanism) for k in keys(dep_exprs)) ==
              spec.expected_dep_struct_keys     # fixture
        @test Set(_struct_key(s, spec.mechanism) for s in indep) ==
              spec.expected_indep_struct_keys
    end
end
```
Capture the `expected_*_struct_keys` fixtures from a pre-Phase-4 run; the test must stay green after **each** of Phases 4, 5, and 6 (re-run as a gate at each phase boundary). If it drifts at Phase 5, a SS iso direction flip silently moved the dep choice between physical kf and physical kr — investigate before regenerating goldens.

- [ ] **Step 6: Run to confirm pass; run full suite; regenerate any goldens whose rep changed.**

Run: `julia --project -e 'using Pkg; Pkg.test()'` → regenerate the (small) set of multi-step-group goldens that shifted; the dep-set snapshot test (Step 5) must stay green. Re-run to PASS.

- [ ] **Step 7: Commit**

```bash
git add -A
git commit -m "refactor: kinetic-group name rep via shared _step_priority (argmin)"
```

---

## Phase 5 — Canonicalize all iso steps (RE and SS) in the physical-forward direction

Goal: iso step storage direction is canonicalized so the `from` side is the side further along the catalytic cycle's substrate→product progression. The direction question is identical for RE iso (one `K`) and SS iso (two `kf`/`kr`); only parameter count differs, so **one canonicalization rule covers both**. Iso canonicalization moves out of `Step` into the `Mechanism` / `AllostericMechanism` constructor (which has the `EnzymeReaction`'s substrate/product sets). RE *binding* canonicalization stays in `Step` (metabolite on `from` side — needs no reaction context). The derivation's "LHS = kf-side" assumption (`rate_eq_derivation.jl:356-360`) becomes physical-forward-aligned.

Algorithm (`_canonical_iso_direction`, per iso step):
- **Tier 1:** `score(sp) = (n_subs_bound, -n_prods_bound)`. Higher = more "from". Atom-balance progression; decides the 95% case.
- **Tier 2 (1-hop graph context):** classify each species by which RE binding steps touch it as the free side — `:substrate_only` / `:product_only` / `:both` / `:neither`. If from=`:product_only` and to=`:substrate_only`, the step is the cycle-closing/regeneration step; forward = exit → entry. Handles Segel Iso Uni Uni `F ⇌ E` fully source-independently.
- **Tier 3 (lex on conformation name):** deterministic fallback for the truly-symmetric tail.

Residual atoms are NOT a tier: by atom conservation, an iso step can't change residual without also changing bound metabolites, so a Tier-1 tie implies a residual tie automatically.

### Task 5.1: Add `_canonical_iso_direction` (unified RE + SS) and apply it in the Mechanism constructor

**Files:** `src/types.jl:140-177` (`Step` constructor — **remove** the RE-iso lex branch at 164-172; RE binding canonicalization at 159-163 stays), `src/types.jl:353-384` (`Mechanism` constructor), `AllostericMechanism` constructor, new helpers placed next to `Mechanism`.

- [ ] **Step 1: Write a failing test** in `test/test_types.jl` covering Tier 1, Tier 2 (the Segel `F ⇌ E` case), and an RE iso reorientation:

```julia
@testset "iso canonicalization (RE + SS, all tiers)" begin
    # Tier 1 (SS iso, score differs): forward = substrate-bound -> product-bound.
    m_fwd = @enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + S <--> E(S)
            E(S)  -->  E(P)
            E(P) <--> E + P
        end
    end
    m_rev = @enzyme_mechanism begin
        substrates: S; products: P
        steps: begin
            E + S <--> E(S)
            E(P)  -->  E(S)        # SS iso written backwards
            E(P) <--> E + P
        end
    end
    @test m_fwd == m_rev

    # Tier 2 (the Segel Iso Uni Uni case): pure conformational F <--> E,
    # tied on Tier 1, decided by entry_kind. Two opposite-source-direction
    # mechanisms canonicalize to the same form.
    s1 = @enzyme_mechanism begin
        substrates: A; products: P
        steps: begin
            E + A <--> E(A); E(A) <--> E(P); E(P) <--> F + P; F <--> E
        end
    end
    s2 = @enzyme_mechanism begin
        substrates: A; products: P
        steps: begin
            E + A <--> E(A); E(A) <--> E(P); E(P) <--> F + P; E <--> F   # last step flipped
        end
    end
    @test s1 == s2
    iso = only(s for grp in EnzymeRates.steps(EnzymeRates.Mechanism(s2))
                   for s in grp
                   if !EnzymeRates.is_binding(s) && !EnzymeRates.is_equilibrium(s) == false  # RE iso
                   && EnzymeRates.bound_metabolite(s) === nothing
                   && EnzymeRates.name(EnzymeRates.from_species(s)) in (:E, :F))
    @test EnzymeRates.name(EnzymeRates.from_species(iso)) == :F  # product-exit
    @test EnzymeRates.name(EnzymeRates.to_species(iso))   == :E  # substrate-entry
end
```

- [ ] **Step 2: Run to confirm fail** — today both equalities fail (SS iso never canonicalized; RE iso canonicalized by raw lex which picks `:E → :F`, not `:F → :E`).

- [ ] **Step 3: Add the canonicalizer + `entry_kind` helper** in `src/types.jl` next to the `Mechanism` constructor:

```julia
# Classify how a species participates in BINDING steps (RE or SS) as the
# FREE side (canonical binding puts the metabolite on the from_species side,
# so the free side IS from_species). Used by `_canonical_iso_direction`
# Tier 2 to decide direction for pure-conformational iso steps where Tier 1
# ties.
#
# Why ALL binding steps (not just RE): the "substrate-entry / product-exit"
# property is a chemistry fact about which forms metabolites enter and
# leave at — it does NOT depend on whether the binding step is rapid-
# equilibrium or steady-state. The DSL parses `<-->` as SS and `⇌` as RE;
# fixtures like Segel Iso Uni Uni (`E + A <--> EA ⇌ EP <--> F + P, F <--> E`)
# use `<-->` throughout, so an RE-only filter would mis-classify both `E`
# and `F` as `:neither` and the F⇌E case would fall through to Tier 3 lex.
function _entry_kind(sp::Species, binding_steps, subs::Set{Symbol}, prods::Set{Symbol})
    has_sub = false; has_prod = false
    for s in binding_steps
        from_species(s) == sp || continue
        b = bound_metabolite(s); b === nothing && continue
        n = name(b)
        n in subs  && (has_sub  = true)
        n in prods && (has_prod = true)
    end
    has_sub && has_prod && return :both
    has_sub  && return :substrate_only
    has_prod && return :product_only
    return :neither
end

# Canonicalize an iso step's storage direction to physical-forward, so
# `from` is further from product-release / closer to substrate-binding.
# Applies to RE iso AND SS iso — the direction question is identical;
# only the parameter count differs. (RE/SS *binding* steps are already
# canonicalized: RE by the Step constructor; SS binding stays as-is —
# kf/koff/kon/koff direction is meaningful per-step regardless.)
function _canonical_iso_direction(s::Step, subs::Set{Symbol}, prods::Set{Symbol},
                                  all_binding_steps::Vector{Step})
    is_binding(s) && return s
    f, t = from_species(s), to_species(s)

    # Tier 1: atom-balance progression.
    score(sp) = (count(m -> name(m) in subs,  bound(sp)),
                -count(m -> name(m) in prods, bound(sp)))
    sf, st = score(f), score(t)
    sf > st && return s
    sf < st && return Step(t, f, nothing, is_equilibrium(s))

    # Tier 2: 1-hop binding (RE+SS) graph context.
    fk = _entry_kind(f, all_binding_steps, subs, prods)
    tk = _entry_kind(t, all_binding_steps, subs, prods)
    fk == :product_only   && tk == :substrate_only && return s
    fk == :substrate_only && tk == :product_only   && return Step(t, f, nothing, is_equilibrium(s))

    # Tier 3: lex fallback (source-independent).
    string(name(f)) ≤ string(name(t)) ? s : Step(t, f, nothing, is_equilibrium(s))
end
```

- [ ] **Step 4: Remove the RE-iso lex branch from `Step`** at `src/types.jl:164-172`. The `elseif is_equilibrium` block goes away; only the RE-binding swap at 159-163 remains in the `Step` constructor. Update the constructor's comment to say "RE binding canonicalizes here; ALL iso steps canonicalize in the Mechanism constructor via `_canonical_iso_direction`."

- [ ] **Step 5: Apply the canonicalizer in `Mechanism`** (`src/types.jl:353-384`):

```julia
function Mechanism(reaction::EnzymeReaction, steps::Vector{Vector{Step}})
    subs  = Set(substrates(reaction))
    prods = Set(products(reaction))
    flat  = Step[s for g in steps for s in g]
    binding_steps = filter(is_binding, flat)
    canon = [[_canonical_iso_direction(s, subs, prods, binding_steps) for s in g] for g in steps]
    # (the existing source_idx auto-assign block operates on `canon` and stays until Phase 6 deletes it.)
    ...
end
```
Mirror in `AllostericMechanism`'s constructor for its `cat_steps`.

- [ ] **Step 6: Run the unit test to PASS.**

- [ ] **Step 7: Audit `src/` AND `test/` for bare `Step(...)` constructions of iso steps.** With RE-iso canonicalization moved out of the Step constructor, any code path that builds opposite-direction iso Steps and compares them (or hashes them into a dedup set) outside a Mechanism context loses determinism. Grep: `grep -rnE 'Step\(' src/ test/`. For each hit verify the construction is followed by `Mechanism(...)`-wrapping before any `==`/`hash`/`in`/`Set` operation. Specifically audit `src/mechanism_enumeration.jl` expansion moves (`_expand_re_to_ss`, `_expand_split_kinetic_group`, dead-end propagation, `init_mechanisms` topology builder) — these construct Steps freely and pass them to `Mechanism(reaction, steps)` at the end, which is fine; just confirm no intermediate dedup or hashing happens on bare iso Steps. For test sites that previously relied on Step `==` for opposite-direction RE iso, rewrite under a Mechanism wrapper.

- [ ] **Step 8: Run the Decision-1 dep-set guard test (Phase 4 Step 5) as a Phase-5 gate.** Iso canonicalization may flip the stored `from`/`to` of some SS iso steps; the pivot scorer's `is_free` check reads from the now-canonical direction. The structural dep-set must stay identical. If it shifts here, the iso flip silently moved the dep choice between physical kf and kr for some mechanism — diagnose before continuing.

```bash
julia --project -e 'using EnzymeRates, Test; include("test/test_dep_set_invariance.jl")'
```

- [ ] **Step 9: Run full suite; regenerate goldens for mechanisms whose iso direction changed; verify numerical oracles.**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
- Golden string churn: regenerate from actuals.
- **Numerical oracle policy (Risk R1).** Physical-forward canonicalization makes oracles written in the natural forward direction already canonical, so slot flips should be rare. If a numerical test fails because the oracle author wrote a step backward, fix the `positional_params` shim mapping for that mechanism (swap f/r slot assignment), never the oracle formula.

- [ ] **Step 10: Commit**

```bash
git add -A
git commit -m "refactor: unified iso canonicalization to physical-forward direction"
```

### Task 5.2: Update comments and CLAUDE.md for the unified iso canonicalization

**Files:** `src/types.jl` (Step ctor comment), `src/rate_eq_derivation.jl:356-360` (derivation source-order comment), `src/mechanism_enumeration.jl:509`, `.claude/CLAUDE.md`

- [ ] **Step 1:** Update the `Step` constructor comment to describe the new split: RE binding canonicalizes here (metabolite on `from`); ALL iso steps (RE and SS) canonicalize in the Mechanism constructor via `_canonical_iso_direction`.

- [ ] **Step 2:** Rewrite the `rate_eq_derivation.jl:356-360` comment. The warning "Using Step's from/to would be wrong" is now stale: iso storage IS canonical and physically forward, so reading direction from Step's `from`/`to` (or equivalently from `rxns`) is correct. The "LHS = kf-side" now means "physical-forward side = kf-side."

- [ ] **Step 3:** Update CLAUDE.md "Canonical Step Form": RE binding canonicalization in the Step constructor is unchanged (metabolite on `from`). ALL iso steps (RE and SS) are canonicalized in the Mechanism constructor via `_canonical_iso_direction` (Tier 1 substrate/product counts → Tier 2 1-hop RE-binding graph context → Tier 3 lex). Replace any text claiming "SS steps are NOT canonicalized" or "RE iso canonicalized by lex in Step."

- [ ] **Step 4: Run full suite to PASS.**

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "docs: unified iso canonicalization (RE + SS) in Canonical Step Form"
```

---

## Phase 6 — Remove `Step.source_idx`

Goal: nothing renders or selects reps from `source_idx` anymore; delete the field and its machinery.

### Task 6.1: Drop `source_idx`-keyed ordering in `_step_parameters`

**Files:** `src/thermodynamic_constr_for_rate_eq_derivation.jl:34-46`

- [ ] **Step 0: Verify the iteration order matches the `source_idx` order today** (proves the rewrite is safe and that `source_idx` was effectively flat-position all along). Add a temporary assertion test that for every `MECHANISM_TEST_SPECS` entry: `[source_idx(s) for (s, _) in _flat_steps(m)] == 1:n_steps(m)`. Run the test — if any spec violates the assertion, the dedup/canonicalization path has produced a `source_idx` ordering that differs from `_flat_steps` and the rewrite below is NOT safe; stop and inspect before proceeding.

```julia
@testset "source_idx == flat iteration order (precondition for removal)" begin
    for spec in MECHANISM_TEST_SPECS
        m = EnzymeRates.Mechanism(spec.mechanism)
        @test [EnzymeRates.source_idx(s) for (s, _) in EnzymeRates._flat_steps(m)] ==
              collect(1:EnzymeRates.n_steps(m))
    end
end
```
Delete the assertion test in Task 6.2 once `source_idx` itself is gone.

- [ ] **Step 1:** `_step_parameters` writes `out[source_idx(s)] = params`. Replace with order-of-iteration over the flat steps (the thermodynamic-constraint columns are built in `reactions(m)` order, which the mechanism already defines without `source_idx`). Rewrite to push in flat-iteration order:

```julia
function _step_parameters(m::Mechanism)
    out = Vector{Vector{Parameter}}()
    for (s, _) in _flat_steps(m)
        push!(out, is_equilibrium(s) ?
            Parameter[is_binding(s) ? Kd(s, :None) : Kiso(s, :None)] :
            Parameter[is_binding(s) ? Kon(s, :None)  : Kfor(s, :None),
                      is_binding(s) ? Koff(s, :None) : Krev(s, :None)])
    end
    out
end
```
Confirm `_flat_steps` order matches the `reactions(m)` / column order used by `_dependent_param_exprs_kernel` (both iterate groups then steps); if not, align them.

- [ ] **Step 2: Run derivation tests to PASS** (no golden change expected — same params, same order).

- [ ] **Step 3: Commit**

```bash
git add src/thermodynamic_constr_for_rate_eq_derivation.jl
git commit -m "refactor: _step_parameters orders by iteration, not source_idx"
```

### Task 6.2: Delete the `source_idx` field and constructor logic

**Files:** `src/types.jl` (`Step` struct 140-177; `source_idx` accessor 183; `Mechanism` constructor 353-384; `AllostericMechanism` constructor), `src/mechanism_enumeration.jl` (expansion moves passing `source_idx=`; dedup density checks 2118-2194), `test/*` (any `source_idx` references)

- [ ] **Step 1: Find every reference**

Run: `grep -rn 'source_idx' src/ test/`

- [ ] **Step 2: Remove the field** from `Step` (line 145), its keyword from the inner constructor (lines 148-149, 174-175), and the `source_idx(s::Step)` accessor (line 183). Update the `Step` doc comment (lines 133-139) to drop the `source_idx` explanation.

- [ ] **Step 3: Simplify the `Mechanism` constructor** (353-384) to the auto-path only — it no longer renumbers anything:

```julia
struct Mechanism
    reaction::EnzymeReaction
    steps::Vector{Vector{Step}}
end
```
Remove the `any_set`/`all_set` validation and renumber loop. Apply the analogous simplification to the `AllostericMechanism` constructor (drop `source_idx` density validation in `mechanism_enumeration.jl:2118-2194`).

- [ ] **Step 4: Remove `source_idx = …` keyword args** in `src/mechanism_enumeration.jl` expansion moves (lines ~1139, 1434, 1448) and any other `Step(...; source_idx=...)` call sites.

- [ ] **Step 5: Run full suite to PASS.** (Numerical + golden unaffected; if a dedup/density test referenced `source_idx`, update or delete that specific assertion.)

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "refactor: delete Step.source_idx and its constructor machinery"
```

---

## Phase 7 — Cleanup sweep

Goal: simplify the now-redundant canonical-hash token layer and remove dead index-era comments. Guarded by the hash-partition tests.

### Task 7.1: Simplify the rate-equation canonical-hash token layer

**Files:** `src/mechanism_enumeration.jl:1856-1900+` (`_parameter_canonical_key` + substitution), `src/identify_rate_equation.jl` (`_canonical_rate_eq_hash_data_impl_struct`)

- [ ] **Step 1: Run the hash-partition test as the safety net first**

Run: `julia --project -e 'using EnzymeRates, Test; include("test/test_canonical_hash_partition.jl")'`
Expected: PASS (baseline).

- [ ] **Step 2: Assess collapse.** With globally structural names, `name(p, m)` is itself a position-independent canonical key. Where `_parameter_canonical_key`/token-substitution exists only to re-derive what the rendered name now encodes, replace the per-Parameter token with the rendered `name(p, m)` Symbol and drop the substitution pass. Make the **smallest** change that keeps `test_canonical_hash_partition.jl` green; do not over-collapse if a test distinguishes a case the name does not encode.

- [ ] **Step 3: Run hash-partition + identify tests to PASS.**

Run: `julia --project -e 'using EnzymeRates, Test; include("test/test_canonical_hash_partition.jl"); include("test/test_identify_rate_equation.jl")'`

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "refactor: canonical hash keys on structural name, drop token layer"
```

### Task 7.1b: Delete the `_step_sides` / `rxns`-walk indirection in the derivation

Round-2 reviewer finding: once iso direction is canonical (Phase 5), the derivation no longer needs to read direction from the `rxns` tuple via `_step_sides` — `from_species(s)` / `to_species(s)` IS the canonical direction, and `rxns` is just a re-projection of those. The "Using Step's from/to would be wrong" warning at `rate_eq_derivation.jl:356-360` is now stale (Task 5.2 already rewrote the comment); this task deletes the underlying indirection.

**Files:** `src/rate_eq_derivation.jl` (`_raw_symbolic_rate_polys`, `_compute_numerator`, `_compute_alpha`, `_compute_re_groups` and any `enumerate(rxns)` walks), `src/thermodynamic_constr_for_rate_eq_derivation.jl` (`_dependent_param_exprs_kernel` line 276+ `enumerate(rxns)` walk), helpers `_step_sides` / `_split_reaction_side`.

- [ ] **Step 1: Map every `enumerate(rxns)` walk in the derivation.** Grep: `grep -rnE 'enumerate\(rxns\)|_step_sides|_split_reaction_side' src/`.

- [ ] **Step 2: Rewrite each walk to iterate `_flat_steps(mech)` and read direction/metabolite directly from the Step.** A typical translation:
```julia
# old:
for (idx, (lhs, rhs, _, _)) in enumerate(rxns)
    e_lhs, m_lhs = _split_reaction_side(lhs, enz_set)
    e_rhs, m_rhs = _split_reaction_side(rhs, enz_set)
    ...
end
# new (after Phase 5 makes Step.from/to canonical):
for (idx, (s, _)) in enumerate(_flat_steps(mech))
    e_lhs = name(from_species(s)); e_rhs = name(to_species(s))
    b = bound_metabolite(s)
    m_lhs = b === nothing ? Symbol[] : [name(b)]   # canonical RE binding: met on from
    m_rhs = Symbol[]                                # SS binding still has met-on-from after Step ctor
    ...
end
```
The exact `m_lhs`/`m_rhs` policy depends on whether the step is RE binding (metabolite on `from`), SS binding (metabolite stays where Step ctor canonicalized it, also `from`), or iso (no metabolite). Verify carefully against each call site.

- [ ] **Step 3: Delete `_step_sides` and `_split_reaction_side` once no callers remain.** Confirm via grep.

- [ ] **Step 4: Run full suite to PASS.** Numerical oracle tests + golden strings are the safety net — if a rewrite gets metabolite-side assignment wrong, numerical oracles fail loudly.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "refactor: derivation reads direction from Step directly; delete _step_sides/rxns walks"
```

### Task 7.1c: Verify Pass-2 Wegscheider absorption never fires, then delete `_build_kinetic_rename_map`

Round-2 reviewer finding + the project memory note (`project_dedup_pass2_dead_code.md`): "Pass-2 single-symbol Wegscheider absorption doesn't fire on real or minimal mechanisms." Task 3.2 already deletes Pass-1. If Pass-2 is also empirically never-fires, the whole `_build_kinetic_rename_map` function — plus the `rename` plumbing through `_dependent_param_exprs_kernel`, the `indep` filter at `thermodynamic_constr…:185`, the `ANNOTATION_SUBSTITUTED` constant (`rate_eq_derivation.jl:12`), and its consumer regex — can be deleted wholesale. ~100-200 LOC of pure deletion.

- [ ] **Step 1: Add an empirical assertion** that runs through every `MECHANISM_TEST_SPECS` entry (and the enumerator's `init_mechanisms` output for a sample of canonical reactions) and asserts `_build_kinetic_rename_map(typeof(em)) == Dict{Symbol,Symbol}()` — empty. Run it.

- [ ] **Step 2: If Step 1 passes,** delete `_build_kinetic_rename_map`. Replace its call sites in `_dependent_param_exprs` and the `rate_equation_string` allosteric block with the empty-Dict equivalent (inlined). Drop the kernel's `rename` parameter; the kernel becomes single-argument. Delete `ANNOTATION_SUBSTITUTED` and the `_canonicalize_for_hash` regex consumer in `identify_rate_equation.jl`. Re-run full suite + hash-partition test.

- [ ] **Step 3: If Step 1 fails** (Pass-2 fires for some mechanism), document the firing mechanism in a comment and keep the function. Note in the commit message.

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "refactor: delete _build_kinetic_rename_map and ANNOTATION_SUBSTITUTED (Pass-2 never fires)"
```

### Task 7.2: Dead-comment + dead-code sweep

**Files:** all of `src/`

- [ ] **Step 1: Grep for stale references**

Run: `grep -rniE 'source_idx|rep_idx|positional|:K1|R-state|T-state|:OnlyR|kNf|index-context' src/`
Expected: only legitimate hits (e.g. `_step_priority` docstring). Remove comments describing the deleted index naming, the deleted index-context companion, and any "future refactor will…" notes now realized.

- [ ] **Step 2:** Per CLAUDE.md code-style, re-read each changed file for dead helpers left unused after the rewrite (e.g. `_site_idx_of`, leftover `_param_symbol` references). Delete unused ones.

- [ ] **Step 3: Run full suite (incl. Aqua/JET) to PASS.**

Run: `julia --project -e 'using Pkg; Pkg.test()'`

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "refactor: remove dead index-era code and comments"
```

### Task 7.3: Final verification

- [ ] **Step 1: Full suite green** — `julia --project -e 'using Pkg; Pkg.test()'`.
- [ ] **Step 2: Confirm deletions** — `grep -rn 'source_idx\|_rep_idx_for_step\|name(::Type{' src/` returns nothing.
- [ ] **Step 3: Confirm structural names** — spot-check `parameters` for a few `MECHANISM_TEST_SPECS` show structural names, no `:k[0-9]`.
- [ ] **Step 4: LOC check** — `git diff --stat main` net source LOC not increased.
- [ ] **Step 5:** Use `superpowers:requesting-code-review` before merge.

---

## Self-review notes (spec coverage)

- Spec Decision 1 (pivot unchanged) → Phase 4 Task 4.1 Step 3 explicitly preserves pivot ordering; only the scorer is extracted.
- Decision 2 (canonicalize SS) → Phase 5.
- Decision 3 (A/I everywhere) → Phase 2 (+ chokepoint token in Phase 3 Task 3.1 `_state_tag`).
- Decision 4 (single `_step_priority`) → Phase 4.
- Decision 5 (hybrid tests) → Phase 0 shim + golden-capture steps throughout.
- Naming table → Phase 3 Task 3.1 chokepoint bodies.
- `name(Species)` concat → Phase 1. `:Et` rename DROPPED per Denis (~15 hard-coded `:E_total` sites; pure churn unrelated to structural naming).
- `_inactive_name` helper (consolidates 11 `string(k) * "_T"` synth-dep sites) → Phase 3 Task 3.1.5. This was the DRY simplification Denis flagged when approving the A/I mid-name placement.
- Phase 5 rewritten to canonicalize **all iso steps (RE and SS)** with one rule — `_canonical_iso_direction` — in the `Mechanism` constructor (which has reaction context). Tier 1 substrate/product counts → Tier 2 1-hop **binding** (RE+SS) graph context (handles the Segel `F ⇌ E` case fully source-independently) → Tier 3 lex fallback. Residual is NOT a tier (atom conservation makes it redundant with Tier 1). The RE-iso lex branch is removed from the `Step` constructor.
- Round-2 reviewer fixes folded in:
  - Tier 2 filter changed from RE-only to all binding steps (the keystone fixture Segel Iso Uni Uni uses `<-->` (SS) throughout, so an RE-only filter would degenerate to lex on its own test). `all_re_steps` → `all_binding_steps`.
  - `_inactive_name` redefined to **swap** `_A_` → `_I_` (not insert), erroring on no-A-token symbols — the original "insert after first underscore" double-encoded state for NonequalAI groups.
  - Phase 0 `positional_params` shim rewritten to walk **flat steps in source order**, emitting per-step `:K{src_idx}`/`:k{src_idx}f/r` names with rep-value lookup, so the 11 multi-step-group analytical oracles (HK, PFK-1, PK, the inhibitor specs, MWC tetramer, etc.) keep working.
  - Decision-1 dep-set guard re-keyed on **structural Parameter identity** (Parameter type + step hash + state), not rendered names — invariant under rep change AND iso flip. Re-runs as a gate at the end of Phase 4, 5, AND 6.
  - Phase 5 audit extended to `src/` Step construction sites, not just `test/`.
- Two-reviewer audits (round 1 + round 2): confirmed perf gate (R3) safe; the round-2 audit caught three blocking bugs in my draft and one large simplification I'd left on the table.
- Phase 7 expanded: in addition to the conservative canonical-hash simplification (Task 7.1), Task 7.1b deletes `_step_sides`/`rxns`-walks in the derivation (Step now is the single source of truth post-canonicalization) and Task 7.1c attempts to delete `_build_kinetic_rename_map` + `ANNOTATION_SUBSTITUTED` entirely (gated on the empirical observation, per project memory, that Pass-2 Wegscheider absorption never fires on real mechanisms). Combined LOC win likely larger than Task 7.1 alone.
- Risk R1 (SS×shim) → Phase 5 Task 5.1 Step 5 explicit policy. R2 (collision) → accepted, noted in `name(p::Kreg)` comment. R3 (perf) → `test_rate_equation_performance` runs in the full suite every phase.
