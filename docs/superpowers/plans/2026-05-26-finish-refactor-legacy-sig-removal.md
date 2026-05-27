# Finish the Concrete-Types Refactor — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. **This work has one exploratory task (A3 validator port) — do NOT delegate it to an unsupervised subagent.** Drive it directly. The rest is deletion / mechanical conversion / movement.

**Goal:** Close the Symbol→concrete-struct refactor by deleting the dead legacy opaque-Symbol code paths, routing all field access through accessors, and reorganizing file layout — all behavior-preserving.

**Architecture:** Six buckets in locked order: **B** (dsl.jl dead `rxns` chain) → **B′** (opaque-rejection guard) → **A** (legacy Sig path + validator audit) → **C** (stale memory) → **D** (accessor discipline) → **E** (file layout). Deletions first, accessors next, pure movement last. Spec: `docs/superpowers/specs/2026-05-26-finish-refactor-legacy-sig-removal-design.md`.

**Tech Stack:** Julia 1.10+, Test/Aqua/JET. Non-negotiable gates run **per commit**: `test_rate_equation_performance` (0-alloc/<100ns), 3 compile-budget gates, `test_chokepoint.jl`. Plus `scripts/check_test_integrity.sh main` (EXIT=0) at every test-touching commit.

---

## Conventions (apply to every task)

- **No `--amend`.** Stay on branch `refactor-to-concrete-types-instead-of-symbols`.
- **Per-file test recipe** (resolve deps once, iterate fast):
  ```bash
  julia --project=. -e '
    using Pkg; Pkg.activate(temp=true); Pkg.develop(path=".")
    Pkg.add(["Test","Aqua","JET","OptimizationBBO","OptimizationPyCMA","OrdinaryDiffEqFIRK","Tables","DataFrames","Statistics","Optimization","Random","CSV"])
    using Test, EnzymeRates
    include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<file>.jl")'
  ```
- **Full suite:** `julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/out.log 2>&1; echo done` then `grep "Test Summary" /tmp/out.log` — **read the summary line**, don't trust the trailing exit code.
- **Test integrity:** any deletion/weakening of a test requires a `docs/superpowers/refactor-deleted-tests.md` §2.1 entry **in the same commit**. Check: `bash scripts/check_test_integrity.sh main; echo "EXIT=$?"` → 0 (no pipe).
- **No temporal-context comments** in code ("legacy", "previously", "Stage N"). Evergreen only.
- **Commit footer:** `src delta: -X / +Y net Z, cumulative: ±W` (vs main 7,136 — `wc -l src/*.jl`). End with `Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>`.

---

## File Structure

| File | Buckets touching it | Responsibility |
|---|---|---|
| `src/dsl.jl` | B, B′, D, E | macros; delete dead `rxns` chain + opaque guard |
| `src/types.jl` | A, D, E | structs + accessors; delete legacy Sig path |
| `src/mechanism_enumeration.jl` | A, D, E | enumeration; `_assert_mechanism_invariants` gains ported validators |
| `src/rate_eq_derivation.jl` | A(comment), D, E | derivation; doc fix + accessors |
| `src/thermodynamic_constr_for_rate_eq_derivation.jl`, `src/fitting.jl`, `src/identify_rate_equation.jl`, `src/sym_poly_for_rate_eq_derivation.jl` | D, E | accessors + layout |
| `test/test_types.jl` | A | migrate/delete legacy-constructor tests |
| `test/test_mechanism_enumeration.jl` | A | add bi-bi exit-gate test; assert ported validators |
| `test/test_dsl.jl` | B′ | delete opaque-rejection tests |
| `docs/superpowers/refactor-deleted-tests.md` | A, B′ | §2.1 entries |

---

# BUCKET B — dsl.jl dead `rxns` chain

### Task B1: Delete the dead legacy emission chain in dsl.jl

**Files:**
- Modify: `src/dsl.jl` — `_parse_steps_block_with_groups` (825–895), call sites (618, 1204), delete `_parse_single_step` (949–959), `_parse_step_side_symbols` (246–258), `_step_side_term_to_symbol` (264–274), `_synthesize_species_name` (389–391), and the `.sym` synthesis in `_call_form_term_info` (369–381).

- [ ] **Step 1: Confirm the chain is dead.**

Run:
```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
sed -n '617,618p;1204,1206p' src/dsl.jl   # both must destructure `_,` (discard rxns)
```
Expected: line 617 `_, side_terms_per_step =`, line 1204 `_, group_tags, side_terms_per_step =`. Confirms `rxns` (1st return slot) is discarded by both callers.

- [ ] **Step 2: Drop `rxns` from `_parse_steps_block_with_groups`.**

In `src/dsl.jl`, edit `_parse_steps_block_with_groups` (825–895):
- Delete line 828 `rxns = Expr(:tuple)`.
- Delete the four `push!(rxns.args, _parse_single_step(...))` blocks (847–848, 861–862, 882) — keep the adjacent `push!(side_terms_per_step, _step_struct_info(...))` calls.
- Change the returns (890–894) to:
```julia
    if allow_tag
        return tags, side_terms_per_step
    else
        return side_terms_per_step
    end
```

- [ ] **Step 3: Update the two call sites.**

`src/dsl.jl:617`:
```julia
    side_terms_per_step =
        _parse_steps_block_with_groups(steps_block, declared_mets)
```
`src/dsl.jl:1204`:
```julia
    group_tags, side_terms_per_step = _parse_steps_block_with_groups(
        cat_steps_block, declared_mets; allow_tag=true,
    )
```

- [ ] **Step 4: Delete the now-orphaned functions.**

Delete from `src/dsl.jl`: `_parse_single_step` (949–959), `_parse_step_side_symbols` (246–258), `_step_side_term_to_symbol` (264–274), `_synthesize_species_name` (389–391). Delete the `.sym` synthesis block inside `_call_form_term_info` (369–381) and remove `sym` from the `_StepSideTerm` it builds (verify by reading 317–388 first; if `_StepSideTerm.sym` field is now unused everywhere, drop the field from the struct definition too — grep `.sym` across `src/dsl.jl` to confirm zero remaining reads).

- [ ] **Step 5: Update touched docstrings.**

Reword any docstring on the touched functions that references "legacy", "rxns tuple", or "mirrors `_parse_step_side_symbols`" to describe the decomposed-only behavior (evergreen, no temporal words). Specifically check lines ~250, 277–278, 314–316, 818–823.

- [ ] **Step 6: Verify no dangling references.**

```bash
grep -rn "_parse_single_step\|_parse_step_side_symbols\|_step_side_term_to_symbol\|_synthesize_species_name" src/ test/
```
Expected: empty.

- [ ] **Step 7: Run the per-commit gates + dsl tests + full suite.**

Run `test_dsl.jl` via the per-file recipe (Expected: all pass), then the full suite (Expected: `Test Summary` shows 0 failures, 0 errors), then:
```bash
bash scripts/check_test_integrity.sh main; echo "EXIT=$?"   # 0
```
Also run `test_rate_equation_performance`, the 3 compile-budget gates, and `test_chokepoint.jl` (per-file recipe). Expected: all green.

- [ ] **Step 8: Commit.**
```bash
git add src/dsl.jl
git commit   # "refactor: delete dead legacy rxns-emission chain in dsl.jl"
```

---

# BUCKET B′ — remove opaque-rejection guard

### Task B′1: Confirm opaque forms still error without the explicit guard (HARD GATE)

**Files:** none (probe only).

- [ ] **Step 1: Write a throwaway probe** at `/tmp/opaque_probe.jl`:
```julia
using Pkg; Pkg.activate(temp=true); Pkg.develop(path=".")
using EnzymeRates
threw = false
try
    Core.eval(Main, quote
        using EnzymeRates
        EnzymeRates.@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES
                ES <--> E + P
            end
        end
    end)
catch err
    global threw = true
    println("opaque form errored: ", first(sprint(showerror, err), 150))
end
threw || error("REGRESSION: opaque form was SILENTLY ACCEPTED — keep the guard")
println("OK — opaque still errors")
```

- [ ] **Step 2: Run it after temporarily removing the guard call.**

Comment out (do NOT yet delete) the `_assert_no_opaque_terms(...)` calls at `src/dsl.jl:620` and `:1208`, then:
```bash
julia --project=. /tmp/opaque_probe.jl
```
Expected: `OK — opaque still errors`.
**If it prints the REGRESSION error instead:** STOP. Opaque forms become silently accepted without the guard → keep the guard, skip Bucket B′ entirely, and report to Denis. Restore the commented calls.

- [ ] **Step 3: If green, proceed; if red, restore guard and stop.** Record the probe's error message — it goes in the §2.1 entry (Task B′2 Step 3).

### Task B′2: Delete the guard + its tests

**Files:**
- Modify: `src/dsl.jl` — delete `_assert_no_opaque_terms` (435–451), `_is_conformation_shape` (428–429), call sites (620, 1208).
- Modify: `test/test_dsl.jl` — delete opaque-rejection testsets (~321–342; read to confirm exact range).
- Modify: `docs/superpowers/refactor-deleted-tests.md` — §2.1 entry.

- [ ] **Step 1: Delete the guard functions + call sites.**

Remove `_assert_no_opaque_terms` (435–451), `_is_conformation_shape` (428–429), and the two call sites (620, 1208, already commented in B′1).

- [ ] **Step 2: Delete the opaque-rejection tests.**

Read `test/test_dsl.jl` around 300–350 to find the exact testset(s) asserting opaque rejection (`@test_throws` on `:ES`/`E_S` forms). Delete them.

- [ ] **Step 3: Add §2.1 entry** to `docs/superpowers/refactor-deleted-tests.md`:
> **§2.1** — `test_dsl.jl` opaque-rejection tests deleted with `_assert_no_opaque_terms`/`_is_conformation_shape`. Opaque bound-form names (`:ES`) are now rejected by the decomposed grammar's natural parse failure (probe 2026-05-26: error `<paste from B′1 Step 3>`), not an explicit guard. No production path emits opaque forms.

- [ ] **Step 4: Verify + gates + full suite.**
```bash
grep -rn "_assert_no_opaque_terms\|_is_conformation_shape" src/ test/   # empty
bash scripts/check_test_integrity.sh main; echo "EXIT=$?"               # 0
```
Run `test_dsl.jl`, full suite, perf/compile/chokepoint gates. Expected: all green.

- [ ] **Step 5: Commit.**
```bash
git add src/dsl.jl test/test_dsl.jl docs/superpowers/refactor-deleted-tests.md
git commit   # "refactor: remove opaque-rejection guard; grammar parse failure suffices"
```

---

# BUCKET A — legacy Sig path

### Task A1: Commit the bi-bi exit-gate as a permanent test (lock the precondition)

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (add a testset near the existing `bi_bi_rxn` const at line 65).

- [ ] **Step 1: Add the test.** After the existing bi-bi enumeration testsets, add:
```julia
@testset "bi-bi exit gate: every init mechanism derives a rate equation" begin
    mechs = EnzymeRates.init_mechanisms(bi_bi_rxn)
    @test length(mechs) == 77
    for (i, m) in enumerate(mechs)
        cm = EnzymeRates.compile_mechanism(m)
        s = EnzymeRates.rate_equation_string(cm)
        @test s isa AbstractString
        @test !isempty(s)
    end
end
```

- [ ] **Step 2: Run it.** Via per-file recipe on `test/test_mechanism_enumeration.jl`.
Expected: PASS (77 mechanisms, all derive). This was proven green 2026-05-26.

- [ ] **Step 3: Commit.**
```bash
git add test/test_mechanism_enumeration.jl
git commit   # "test: bi-bi exit gate — all init mechanisms derive (Phase 2 precondition)"
```

### Task A2: Collapse the 12 `@generated` accessor branches onto the new-Sig body

**Files:**
- Modify: `src/types.jl` — accessors with `if _is_new_sig(Sig)` at lines 1325, 1338, 1351, 1369, 1470, 1484, 1493, 1500, 1511, 1522, 1545, 1597 (and `Mechanism(em)` at 712–714).

- [ ] **Step 1: For each accessor, read the full function, then collapse it.**

Pattern — each accessor currently looks like:
```julia
@generated function substrates(::EnzymeMechanism{Sig}) where {Sig}
    if _is_new_sig(Sig)
        <NEW BODY>
        return X
    end
    <LEGACY BODY>          # e.g. Sig[1][1]
end
```
Collapse to (drop the guard + the trailing legacy body, keep the new body un-indented):
```julia
@generated function substrates(::EnzymeMechanism{Sig}) where {Sig}
    <NEW BODY>
    return X
end
```
For the ternary form at 1493 (`_is_new_sig(Sig) ? A : B`), keep only `A`.

**Do this one accessor at a time.** After each, run `test_accessors.jl` (per-file recipe) — Expected: PASS (it's on decomposed grammar, exercises the new body). A regression here means a wrong collapse.

- [ ] **Step 2: Simplify `Mechanism(em)` (712–714).**
```julia
Mechanism(em::EnzymeMechanism{Sig}) where {Sig} = _mechanism_from_sig(Sig)
```

- [ ] **Step 3: Run `test_accessors.jl` + full suite + gates.**
Expected: all green. Confirm `test_rate_equation_performance` 0-alloc still passes.

- [ ] **Step 4: Commit.**
```bash
git add src/types.jl
git commit   # "refactor: collapse EnzymeMechanism accessors onto the single decomposed-Sig body"
```

### Task A3: Audit + port the meaningful legacy validators into `_assert_mechanism_invariants` (EXPLORATORY — drive directly)

**Files:**
- Modify: `src/mechanism_enumeration.jl` — `_assert_mechanism_invariants(m::Mechanism)` (2122–2143).
- Modify: `test/test_mechanism_enumeration.jl` — add tests for the ported invariants.

The audit (from spec Bucket A table): **port** "every metabolite appears in some step" and the kinetic-group "RE/SS-mixing + same-metabolite-per-group" checks (the iso/binding `bound_metabolite` consistency is already in `_assert_mechanism_invariants`); **drop** the rest (moot: typed-struct-guaranteed; superseded: connectivity + stoich, per `test_types.jl:322–331` design intent).

- [ ] **Step 1: Write failing tests for the two ported invariants** in `test/test_mechanism_enumeration.jl`:
```julia
@testset "ported invariants: unused metabolite + kinetic-group composition" begin
    # Build a Mechanism whose reaction declares a regulator R that no step uses.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: R
    end
    s1 = EnzymeRates.Step(EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E),
                          EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E_S),
                          EnzymeRates.Substrate(:S), true; source_idx = 1)
    s2 = EnzymeRates.Step(EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E_S),
                          EnzymeRates.Species([EnzymeRates.Product(:P)], :E_P),
                          nothing, false; source_idx = 2)
    s3 = EnzymeRates.Step(EnzymeRates.Species([EnzymeRates.Product(:P)], :E_P),
                          EnzymeRates.Species(EnzymeRates.Metabolite[], :E),
                          EnzymeRates.Product(:P), true; source_idx = 3)
    m_unused = EnzymeRates.Mechanism(rxn, [[s1], [s2], [s3]])
    @test_throws ErrorException EnzymeRates._assert_mechanism_invariants(m_unused)

    # A kinetic group binding two different metabolites → error.
    rxn2 = @enzyme_reaction begin
        substrates: S[C], A[N]
        products: P[CN]
    end
    g1a = EnzymeRates.Step(EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E),
                           EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E_S),
                           EnzymeRates.Substrate(:S), true; source_idx = 1)
    g1b = EnzymeRates.Step(EnzymeRates.Species([EnzymeRates.Substrate(:A)], :E),
                           EnzymeRates.Species([EnzymeRates.Substrate(:A)], :E_A),
                           EnzymeRates.Substrate(:A), true; source_idx = 2)
    g2  = EnzymeRates.Step(EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E_S),
                           EnzymeRates.Species([EnzymeRates.Product(:P)], :E_P),
                           nothing, false; source_idx = 3)
    m_mixed_met = EnzymeRates.Mechanism(rxn2, [[g1a, g1b], [g2]])
    @test_throws ErrorException EnzymeRates._assert_mechanism_invariants(m_mixed_met)
end
```
**Note:** verify the exact `Species`/`Step`/`Metabolite` constructor signatures by reading `src/types.jl:11-160` before writing these — adjust the fixture to the real constructors. The shapes above follow the decomposed grammar (`E(S)` = `Species([Substrate(:S)], :E)`).

- [ ] **Step 2: Run the test — confirm it FAILS** (current `_assert_mechanism_invariants` doesn't check these).
Expected: FAIL (no error thrown → `@test_throws` fails).

- [ ] **Step 3: Port the two checks into `_assert_mechanism_invariants(m::Mechanism)`** (after the existing per-step loop, before `nothing`):
```julia
    # Every declared metabolite/regulator must appear in some step.
    appearing = Set{Symbol}()
    for s in flat
        for sp in (from_species(s), to_species(s))
            for met in bound(sp); push!(appearing, name(met)); end
        end
        bm = bound_metabolite(s)
        bm === nothing || push!(appearing, name(bm))
    end
    for nm in metabolites(reaction(m))
        nm in appearing ||
            error("declared metabolite/regulator $nm appears in no step")
    end

    # Within a kinetic group of size > 1: all binding (no RE/SS mix), same metabolite.
    for group in m.steps
        length(group) == 1 && continue
        kinds = [(s.is_equilibrium, bound_metabolite(s)) for s in group
                 if bound_metabolite(s) !== nothing]
        isempty(kinds) && continue
        first_eq, first_met = kinds[1]
        for (eq, met) in kinds[2:end]
            eq == first_eq ||
                error("kinetic group mixes RE and SS binding steps")
            met == first_met ||
                error("kinetic group binds different metabolites: " *
                      "$(name(first_met)) and $(name(met))")
        end
    end
```
**Note:** confirm `metabolites(reaction(m))` returns the declared substrate∪product∪regulator names (read `src/types.jl` accessors first); if the accessor name differs, use the correct one. `name(::Substrate/Product/CompetitiveInhibitor/AllostericRegulator)` already exists (types.jl:23–26), so `name(met)` is valid here.

- [ ] **Step 4: Run the new test — confirm it PASSES**, then run the **full** `test_mechanism_enumeration.jl` — Expected: PASS, including the ~80 existing `_assert_mechanism_invariants` calls (the ported checks must not reject any currently-valid enumerated mechanism). **If any existing call now errors, STOP** — a real enumerated mechanism violates a ported invariant, meaning the invariant is wrong for the decomposed world; re-examine before proceeding.

- [ ] **Step 5: Commit.**
```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit   # "refactor: port unused-metabolite + kinetic-group-composition checks to _assert_mechanism_invariants"
```

### Task A4: Migrate / delete the legacy-constructor tests in test_types.jl

**Files:**
- Modify: `test/test_types.jl` — testsets at 76–114, 265–279, 281–332, 334–342, 344–360, 362–370.
- Modify: `docs/superpowers/refactor-deleted-tests.md` — §2.1 entries for dropped tests.

- [ ] **Step 1: Migrate the positive construction tests** (76–97, 265–279, 334–342) to decomposed grammar.

Replace the opaque-tuple `EnzymeMechanism(mets, rxns)` constructions with `@enzyme_mechanism` decomposed-grammar fixtures that produce the same structure, preserving the assertions (`n_steps`, `substrates`, `kinetic_group` sharing, `isa EnzymeMechanism`). Example for 76–85:
```julia
@testset "EnzymeMechanism construction (decomposed)" begin
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E(P)
            E(P) <--> E + P
        end
    end
    @test EnzymeRates.n_steps(m) == 3
    @test EnzymeRates.substrates(m) == (:S,)
end
```
For the same-kinetics group test (87–97), use a competitive-inhibitor fixture binding E and E(S) in one group (see `test_dsl.jl`/`test_accessors.jl` for the grammar) and assert `kinetic_group` equality.

- [ ] **Step 2: Re-point the *ported*-validator `@test_throws`** (kinetic-group different-metabolite at 348–352; RE/SS mix at 355–359) onto `_assert_mechanism_invariants` over a decomposed `Mechanism` (reuse the A3 fixtures' shape). These move to `test_mechanism_enumeration.jl` only if cleaner; otherwise keep in `test_types.jl` asserting via `_assert_mechanism_invariants`.

- [ ] **Step 3: Delete the moot/superseded `@test_throws`** with §2.1 entries:
  - 104 (stoich "vanish"), 113 (iso group size>1 — iso singleton is in `_assert_mechanism_invariants`? confirm; if not covered, this is a *ported* check, move to Step 2), 289 (empty reactions — covered by `_assert_mechanism_invariants` `isempty`), 292 (dup substrate names — confirm where enforced now), 300 (zero enzyme — moot), 308 (two metabolites — moot), 315 (unknown metabolite — moot, unrepresentable in decomposed `Step`), 320 (net stoich mismatch — superseded), 362–370 (orphan connectivity — superseded per design note).
  - §2.1 entry text per the spec audit table reasons.

- [ ] **Step 4: Run `test_types.jl` + full suite + integrity.**
```bash
bash scripts/check_test_integrity.sh main; echo "EXIT=$?"   # 0
```
Expected: all green; `@testset` count change accounted for by §2.1 entries.

- [ ] **Step 5: Commit.**
```bash
git add test/test_types.jl docs/superpowers/refactor-deleted-tests.md
git commit   # "test: migrate legacy-constructor tests to decomposed grammar; drop moot/superseded validators (§2.1)"
```

### Task A5: Delete the legacy Sig path itself

**Files:**
- Modify: `src/types.jl` — `_is_new_sig` (1292–1308) + doc block (1280–1287), `_mechanism_from_legacy_sig` (727–~795), 2-arg `EnzymeMechanism` constructor (803–888) + `_validate_kinetic_groups`/`_validate_enzyme_connectivity`/`_validate_stoichiometry`/`_pretty_reaction` (894–971) **if now unused**, stale doc comments (641, 649, 796, 1275, 1291).
- Modify: `src/rate_eq_derivation.jl:189` (doc comment).

- [ ] **Step 1: Confirm the validators are now unused in src/.**
```bash
grep -rn "_validate_kinetic_groups\|_validate_enzyme_connectivity\|_validate_stoichiometry\|_pretty_reaction" src/ test/
```
Expected: only their definitions + the 2-arg constructor body. (Task A4 removed the test callers.) If a test still references them, resolve that first.

- [ ] **Step 2: Delete the 2-arg constructor + its now-orphaned validators.** Remove `EnzymeMechanism(mets, rxns)` (803–888) and `_validate_kinetic_groups`, `_validate_enzyme_connectivity`, `_validate_stoichiometry`, `_pretty_reaction` (894–971) **only if** Step 1 showed them unused.

- [ ] **Step 3: Delete `_mechanism_from_legacy_sig`** (727–~795). It is no longer referenced (Task A2 made `Mechanism(em)` call `_mechanism_from_sig` directly).

- [ ] **Step 4: Delete `_is_new_sig`** (1292–1308) + its doc-comment block (1280–1287).

- [ ] **Step 5: Confirm `_legacy_step_tuple`/`_species_sym` are gone** and rename the misnamed reconstructor.
```bash
grep -rn "_legacy_step_tuple\b\|_species_sym\b" src/   # expect empty (already absent)
```
`_step_tuple_from_sig` (types.jl:1443) is already correctly named — no rename needed (the earlier `_legacy_step_tuple_from_sig` was renamed in commit 388944c). Confirm no `_legacy_*` names remain:
```bash
grep -rn "_legacy" src/
```
Expected: empty (or only the rate_eq_derivation.jl:189 comment, fixed next).

- [ ] **Step 6: Fix stale doc comments.** types.jl lines 641, 649, 796, 1275, 1291 (references to `EnzymeMechanism(metabolites, reactions)`) and rate_eq_derivation.jl:189 — reword to describe the current single Sig shape; no temporal words.

- [ ] **Step 7: Verify the closure of "one struct family".**
```bash
grep -rnE "_is_new_sig|_mechanism_from_legacy_sig|EnzymeMechanism\(metabolites" src/   # empty
```

- [ ] **Step 8: Full suite + integrity + all gates.** Expected: all green. Re-run `test_rate_equation_performance` (0-alloc), 3 compile budgets, `test_chokepoint.jl`.

- [ ] **Step 9: Commit.**
```bash
git add src/types.jl src/rate_eq_derivation.jl
git commit   # "refactor: delete legacy opaque-Sig path — one struct family"
```

---

# BUCKET C — stale memory

### Task C1: Delete the resolved accessor-allocation memory

**Files:**
- Delete: `/home/denis.linux/.claude/projects/-home-denis-linux--julia-dev-EnzymeRates/memory/project-accessor-allocates-on-decomposed-sig.md`
- Modify: that dir's `MEMORY.md` — remove the pointer line.

- [ ] **Step 1: Confirm resolution.** `test_accessors.jl` is on decomposed grammar (`E + S <--> E(S)`) and asserts `== 0` allocs, green in Bucket A. The accessor-allocation concern is resolved (commit `bc1c592`).
- [ ] **Step 2: Delete the memory file** and remove its line from `MEMORY.md`. (Not a git-tracked repo file — no commit needed; this is the assistant's memory store.)

---

# BUCKET D — accessor discipline (full conversion, no enforcement test)

### Task D1: Accessor inventory + add the missing accessors

**Files:**
- Modify: `src/types.jl` — add any missing field accessors.

- [ ] **Step 1: Enumerate every concrete-struct field and its accessor.** Read `src/types.jl` struct definitions (lines ~11–415). For each field, confirm an accessor exists; build the gap list. Known existing: `name` (Metabolite family, 23–26), `bound`/`conformation`/`residual` (Species), `from_species`/`to_species`/`bound_metabolite` (Step), `reaction`/`steps`/`regulatory_sites`. Likely missing (verify): `is_equilibrium`/`source_idx` (Step), `cat_allo_states`/`catalytic_multiplicity` (AllostericMechanism — note `.cat_steps` maps to existing `steps(m)`), `atoms` (ReactantAtoms), `ligands` (RegulatorySite, RegulatorMults).

- [ ] **Step 2: Write a failing test** in `test/test_accessors.jl` asserting each newly-added accessor returns the field value (use an existing fixture `m`/`am`/a constructed `Step`/`Species`). Example:
```julia
@testset "added field accessors" begin
    st = EnzymeRates.Step(EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E),
                          EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E_S),
                          EnzymeRates.Substrate(:S), true; source_idx = 1)
    @test EnzymeRates.is_equilibrium(st) == true
    @test EnzymeRates.source_idx(st) == 1
end
```
Adjust to the actual missing set from Step 1.

- [ ] **Step 3: Run — confirm FAIL** (accessors undefined).

- [ ] **Step 4: Add the accessors** next to the related existing ones in `src/types.jl`, one-liners:
```julia
is_equilibrium(s::Step) = s.is_equilibrium
source_idx(s::Step)     = s.source_idx
name(m::Substrate) = m.name        # + Product / CompetitiveInhibitor / AllostericRegulator
# … etc for the confirmed-missing set
```
For the `Metabolite` family `name`, add one method per concrete subtype (or `name(m::Metabolite) = m.name` if the field is shared via the abstract type's layout — verify). **Do not** collide with `name(p::Parameter, m)` / `name(::Type{P}, idx)` — different signatures, fine.

- [ ] **Step 5: Run — confirm PASS.** Then full suite + perf gate (new accessors inline → 0-alloc unaffected).

- [ ] **Step 6: Commit.**
```bash
git add src/types.jl test/test_accessors.jl
git commit   # "refactor: add field accessors for remaining concrete-struct fields"
```

### Task D2–D8: Convert `.field` → accessor, one file per task

Run once per file in this list (each its own task + commit):
`src/types.jl`, `src/mechanism_enumeration.jl`, `src/rate_eq_derivation.jl`, `src/thermodynamic_constr_for_rate_eq_derivation.jl`, `src/fitting.jl`, `src/identify_rate_equation.jl`, `src/dsl.jl`, `src/sym_poly_for_rate_eq_derivation.jl`.

**For each file:**

- [ ] **Step 1: List the field accesses in this file.**
```bash
grep -nE "\.(reaction|steps|cat_steps|cat_allo_states|catalytic_multiplicity|regulatory_sites|from_species|to_species|bound_metabolite|is_equilibrium|source_idx|bound|conformation|residual|ligands|atoms|name)\b" src/<file>
```

- [ ] **Step 2: Replace each `x.field` with `field(x)`**, **except**:
  - inside accessor definitions themselves (`foo(x) = x.foo`),
  - inside struct **constructors** (`new(...)`, field reads during construction),
  - inside `@generated` function bodies that index the `Sig` tuple (those touch tuple slots, not struct fields),
  - `AllostericMechanism.cat_steps` → use `steps(m)` (not a new `cat_steps` accessor).

  Read each match in context before converting; a mechanical regex replace will wrongly hit the exempt sites.

- [ ] **Step 3: Run this file's test(s) + full suite.** Expected: green. For `src/types.jl` and `src/rate_eq_derivation.jl`, also run `test_rate_equation_performance` + 3 compile budgets + `test_chokepoint.jl` (Expected: 0-alloc/<100ns unchanged).

- [ ] **Step 4: Commit.**
```bash
git add src/<file>
git commit   # "refactor: route field access through accessors in <file>"
```

- [ ] **Final D-step: confirm no stray field access remains** (outside the exempt definitions/constructors/@generated):
```bash
grep -rnE "\.(from_species|to_species|bound_metabolite|conformation|residual)\b" src/ | grep -v "= .*\."   # spot-check; expect only definitions
```

---

# BUCKET E — file layout reorganization (pure movement, last)

### Task E1–E8: Reorder each src file, one file per task

Run once per file (each its own commit). **Behavior-preserving movement only.**

**For each file:**

- [ ] **Step 1: Read the whole file** and sketch a target order: ABOUTME header → major public types (for `types.jl`) → major public functions → grouped helper clusters (functions that do similar things adjacent) → minor/private helpers at the end.

- [ ] **Step 2: Apply the reorder**, respecting Julia constraints:
  - A `struct` must precede any method signature or struct field referencing its type → all `struct`/`abstract type` blocks stay at the top of `types.jl`, in dependency order.
  - A macro must be defined before its use site is parsed → in `dsl.jl`, keep `@enzyme_reaction`/`@enzyme_mechanism`/`@allosteric_mechanism` and their parse-helpers in an order where helpers precede the macro that calls them at expansion.
  - `const`s used at top level precede their use.
  - Plain function call resolution is order-independent — helpers may move freely subject to the above.

- [ ] **Step 3: Run the file's tests + full suite + all gates.** Expected: green. **A failure means a genuine order dependency** (a macro used before definition, a const, a type) — surface it and fix the ordering; do not paper over.

- [ ] **Step 4: Commit.**
```bash
git add src/<file>
git commit   # "refactor: reorganize <file> layout — important definitions first"
```

Suggested file order (lowest-risk first): `sym_poly_for_rate_eq_derivation.jl`, `fitting.jl`, `thermodynamic_constr_for_rate_eq_derivation.jl`, `identify_rate_equation.jl`, `rate_eq_derivation.jl`, `mechanism_enumeration.jl`, `dsl.jl`, `types.jl` (riskiest — struct ordering — last).

---

# FINAL — verification, LOC, docs, PR

### Task F1: Success-criteria verification + docs + PR

- [ ] **Step 1: §3 success criteria.**
```bash
grep -rnE "_is_new_sig|_mechanism_from_legacy_sig|EnzymeMechanism\(metabolites" src/   # empty (§3.1 #1)
bash scripts/check_test_integrity.sh main; echo "EXIT=$?"                              # 0 (§3.1 #3)
```
`test_chokepoint.jl` green (§3.1 #2); `test_rate_equation_performance` 0-alloc/<100ns + 3 compile budgets green (§3.1 #4).

- [ ] **Step 2: Full suite green + @testset count.**
```bash
julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/out.log 2>&1; echo done
grep "Test Summary" /tmp/out.log
```
Expected: 0 failures, 0 errors. Confirm `@testset` count ≥ main minus the documented §2.1 removals.

- [ ] **Step 3: LOC re-baseline.**
```bash
wc -l src/*.jl   # record total
```
Write the honest achieved total into the PR body; record it in the spec / roadmap as the renegotiated measured target (LOC non-gating).

- [ ] **Step 4: Docs final pass.** Re-read README + CLAUDE.md for opaque-form / legacy-Sig references invalidated by the deletions. Update continuation spec §11 status to "done — legacy Sig removed, accessors unified, layout reorganized."

- [ ] **Step 5: Push + PR.**
```bash
git push -u origin refactor-to-concrete-types-instead-of-symbols
gh pr create --fill
```
PR body: honest LOC numbers; celebrate the single-struct-family closure; link the deferred parameter-naming refactor (continuation spec §7) and the deferred `_n_fit_params_estimate` undercount as the next steps. End the body with the Claude Code attribution footer.

---

## Self-Review Notes

- **Spec coverage:** B → Task B1; B′ → Tasks B′1, B′2; A → Tasks A1–A5; C → Task C1; D → Tasks D1, D2–D8; E → Tasks E1–E8; Final → Task F1. All buckets covered.
- **Riskiest tasks:** A3 (exploratory validator port — drive directly, stop if a real mechanism violates a ported invariant), A2 (12 accessor collapses — one at a time, `test_accessors.jl` after each), E on `types.jl` (struct ordering — last, smallest moves).
- **Hard gates with stop-conditions:** B′1 (opaque must still error or keep guard), A3 Step 4 (no enumerated mechanism may fail a ported invariant), A5 Step 1 (validators must be unused before deletion).
- **Don't trust `Pkg.test` exit codes** — read the `Test Summary` line.
- **Line numbers drift** as deletions land — each task re-greps/re-reads before editing rather than trusting the numbers from this plan verbatim.
