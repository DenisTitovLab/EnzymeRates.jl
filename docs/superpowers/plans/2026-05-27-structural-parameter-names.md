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
| `src/types.jl` | `name(Species)` concat; chokepoint `_param_symbol`/`name` structural rewrite; delete `name(::Type{P},idx)`, `_rep_idx_for_step`; `Step.source_idx` removal; `Mechanism`/`AllostericMechanism` constructor cleanup; `RegulatorySite` A/I validator; `name(::Etot)=:Et` |
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
Re-key a structural-named parameter NamedTuple to the positional names
(`K1`, `k1f`, `k1r`, …) that hand-derived analytical oracles destructure.
Walks `m`'s kinetic groups in storage order, emitting one positional name
per group (RE) or two (SS, forward then reverse), and pairs them with the
structural names from `parameters(m)` in the same order. Permanent test
utility: oracles are inherently positional, so they bind to slot, not to
the structural identity.
"""
function positional_params(m, nt::NamedTuple)
    mech = m isa EnzymeRates.Mechanism ? m : EnzymeRates.Mechanism(m)
    structural = collect(EnzymeRates.parameters(m))
    positional = Symbol[]
    i = 0
    for group in EnzymeRates.steps(mech)
        i += 1
        rep = first(group)
        if EnzymeRates.is_equilibrium(rep)
            push!(positional, Symbol("K", i))
        else
            push!(positional, Symbol("k", i, "f"))
            push!(positional, Symbol("k", i, "r"))
        end
    end
    @assert length(positional) == length(structural) "positional/structural length mismatch for $(EnzymeRates.name(EnzymeMechanism(mech)))"
    NamedTuple{Tuple(positional)}(Tuple(nt[s] for s in structural))
end
```

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
name(::Etot,  _) = :Et
name(::Lallo, _) = :L
```

Delete `_param_symbol` (all methods), `name(::Type{P}, idx)` / `name(::Type{P}, idx, state)`, the old `_rep_idx_for_step`, and `_site_idx_of` if now unused. Delete the index-context comment block at lines 1396-1425.

> Note on `Krev`: with SS canonicalization not yet in place (Phase 5), `from`/`to` follow stored order. `Kfor`/`Krev` name the forward and reverse directed transitions of the rep step, so the pair is direction-independent regardless of storage order.

- [ ] **Step 4: Run the unit test to confirm pass.**

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/test_types.jl
git commit -m "refactor: chokepoint renders structural parameter names"
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

- [ ] **Step 2: Delete the Pass-1 display loop** in `rate_equation_string` (`src/rate_eq_derivation.jl:617-635`). The `user_lines` block annotated `"K9 = K4"` ties; with structural names, tied members share one symbol, so there is nothing to annotate. Remove the `user_lines` accumulation and any now-empty `# user-defined ... equalities:` header that fed it. Keep the Wegscheider Pass-2 annotation block that follows.

- [ ] **Step 3: Confirm no index-context callers remain**

Run: `grep -rnE 'name\((Kd|Kiso|Kon|Koff|Kfor|Krev|Kreg),\s*[a-z_]*idx' src/`
Expected: no matches.

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
_step_lex_key(s::Step) =
    (String(name(from_species(s))), String(name(to_species(s))))

_group_rep(group::Vector{Step}, free_enz_set::Set{Symbol}) =
    argmin(s -> (_step_priority(s, free_enz_set), _step_lex_key(s)), group)
```
Route the chokepoint's `_rep_step` (Task 3.1) and `_enumerate_parameters_full`/`_onlyA_parameters`/`_I_rename_parameters` "rep = first(group)" through `_group_rep`. `_rep_step` becomes:

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

- [ ] **Step 5: Run to confirm pass; run full suite; regenerate any goldens whose rep changed.**

Run: `julia --project -e 'using Pkg; Pkg.test()'` → regenerate the (small) set of multi-step-group goldens that shifted, re-run to PASS.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "refactor: kinetic-group name rep via shared _step_priority (argmin)"
```

---

## Phase 5 — Canonicalize SS steps

Goal: SS steps canonicalize storage direction like RE steps; delete the "SS not canonicalized" special case. Numerical oracles are the safety net (Risk R1).

### Task 5.1: Canonicalize SS step direction in the `Step` constructor

**Files:** `src/types.jl:146-176`

- [ ] **Step 1: Write a failing test** in `test/test_types.jl`:

```julia
@testset "SS steps canonicalize direction" begin
    a = Species([Substrate(:S)], :E); b = Species([Product(:P)], :E)
    s1 = Step(a, b, nothing, false)   # SS iso, not equilibrium
    s2 = Step(b, a, nothing, false)   # opposite source direction
    @test s1 == s2                    # now structurally equal
end
```

- [ ] **Step 2: Run to confirm fail** (today SS iso steps preserve direction → unequal).

- [ ] **Step 3: Generalize canonicalization.** In the `Step` constructor (lines 150-173), drop the `is_equilibrium &&` guards so binding and iso steps canonicalize regardless of RE/SS:
  - binding: put `bound_metabolite` on the `from_species` side (drop the `is_equilibrium &&` at line 161).
  - iso: order by lex on `name(from_species)` (drop the `elseif is_equilibrium` restriction at line 164 so it also applies to SS iso).

  Update the constructor comment (lines 150-172) to state all steps canonicalize. Update the CLAUDE.md "Canonical Step Form" invariant in a later doc step.

- [ ] **Step 4: Run the unit test to PASS.**

- [ ] **Step 5: Run full suite; regenerate goldens for any mechanism whose SS step storage flipped; VERIFY numerical oracles.**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
- Golden string churn: regenerate from actuals.
- **If a numerical analytical test fails:** a step's canonical direction flipped relative to its oracle's positional assumption. Per Risk R1, fix the **`positional_params` shim mapping** for that mechanism (e.g. swap the f/r assignment or group order to match the oracle), never the oracle formula. Re-run to PASS.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "refactor: canonicalize SS step direction; delete special case"
```

### Task 5.2: Remove SS special-casing in dedup canonical keys + CLAUDE.md

**Files:** `src/mechanism_enumeration.jl` (`_step_canonical_key` and any SS-direction comments), `.claude/CLAUDE.md`

- [ ] **Step 1:** Confirm `_step_canonical_key` (line 1753) needs no SS branch (it already keys on `from`/`to` hashes; with canonicalization these are now direction-stable). Remove now-stale comments at `src/types.jl:152-157,165-169` and `mechanism_enumeration.jl:509` that describe SS as direction-preserving.

- [ ] **Step 2:** Update CLAUDE.md "Canonical Step Form" section: all steps (RE and SS) canonicalize; delete the "SS steps are NOT canonicalized" bullets and the rationale about `:kNf`.

- [ ] **Step 3: Run full suite to PASS.**

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "refactor: drop SS-not-canonicalized special case + docs"
```

---

## Phase 6 — Remove `Step.source_idx`

Goal: nothing renders or selects reps from `source_idx` anymore; delete the field and its machinery.

### Task 6.1: Drop `source_idx`-keyed ordering in `_step_parameters`

**Files:** `src/thermodynamic_constr_for_rate_eq_derivation.jl:34-46`

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
- `name(Species)` concat → Phase 1. `Et` → Phase 3 Task 3.1.
- Risk R1 (SS×shim) → Phase 5 Task 5.1 Step 5 explicit policy. R2 (collision) → accepted, noted in `name(p::Kreg)` comment. R3 (perf) → `test_rate_equation_performance` runs in the full suite every phase.
