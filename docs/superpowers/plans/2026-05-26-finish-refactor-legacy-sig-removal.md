# Finish the Concrete-Types Refactor — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. **This work has one exploratory task (A3 validator port) — do NOT delegate it to an unsupervised subagent.** Drive it directly. The rest is deletion / mechanical conversion / movement.

**Goal:** Delete the dead *legacy* opaque-Symbol code paths, route field access through accessors, and reorganize file layout — all behavior-preserving. **This unifies the FRONT END only.** It does NOT achieve "one struct family / no parallel representations": the derivation back-end still runs on regenerated opaque Symbol tuples (`_species_name_from_sig`/`_step_tuple_from_sig` → `reactions`/`enzyme_forms` → King-Altman/Wegscheider). That back-end unification is a deferred next phase, not this PR (see spec "Deferred").

**Architecture:** Five buckets in locked order: **B** (dsl.jl dead `rxns` chain) → **A** (legacy Sig path + validator audit) → **C** (stale memory) → **D** (accessor discipline) → **E** (file layout). Deletions first, accessors next, pure movement last. **Bucket B′ (remove opaque guard) was DROPPED** — review proved the guard is load-bearing (without it, opaque `ES <--> E + P` is silently accepted, product dropped). Spec: `docs/superpowers/specs/2026-05-26-finish-refactor-legacy-sig-removal-design.md`.

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
| `src/dsl.jl` | B, D, E | macros; delete dead `rxns` chain (keep opaque guard + `.sym` field) |
| `src/types.jl` | A, D, E | structs + accessors; delete legacy Sig path |
| `src/mechanism_enumeration.jl` | A, D, E | enumeration; `_assert_mechanism_invariants` gains ported validators |
| `src/rate_eq_derivation.jl` | A(comment), D, E | derivation; doc fix + accessors |
| `src/thermodynamic_constr_for_rate_eq_derivation.jl`, `src/fitting.jl`, `src/identify_rate_equation.jl`, `src/sym_poly_for_rate_eq_derivation.jl` | D, E | accessors + layout |
| `test/test_types.jl` | A | migrate ALL legacy-Sig constructions (2-arg AND `EnzymeMechanism{(...)}()`) |
| `test/test_mechanism_enumeration.jl` | A | add bi-bi exit-gate (subset) test; assert ported validators |
| `docs/superpowers/refactor-deleted-tests.md` | A | §2.1 entries |

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

- [ ] **Step 4: Delete the now-orphaned functions — but KEEP the `.sym` field.**

Delete from `src/dsl.jl`: `_parse_single_step` (949–959), `_parse_step_side_symbols` (246–258), `_step_side_term_to_symbol` (264–274), `_synthesize_species_name` (389–391).

**Do NOT drop the `_StepSideTerm.sym` field.** Review verified it's read on the live emission path: `_build_step_expr` reads `bound_met_term.sym` (`dsl.jl:705`), `_split_side` reads it (720/726/734), and the opaque guard reads it (445). For `:metabolite`/`:bare_enzyme` terms `sym` is the real name. Only the dead `:call`-branch *synthesis* (the `parts`/`join`/`Symbol(...)` block at 369–378) is removable — replace it so the `:call` `_StepSideTerm`'s first arg is the bare conformation name (read 317–388 first to get the exact constructor call), not a string-joined synthesized symbol. Fix the misleading "legacy synthesized" comment (~407) while there.

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

# BUCKET A — legacy Sig path

> **Ordering is review-corrected.** The accessor collapse (A5) MUST come after
> all legacy-Sig construction is gone — `test_types.jl` calls accessors on
> legacy-Sig mechanisms (built both via the 2-arg constructor AND directly via
> `EnzymeMechanism{(mets,rxns)}()` at lines 9, 33, 67, 413). Order: A1 exit
> gate → A2 port validators → A3 migrate ALL legacy-Sig tests → A4 delete
> 2-arg constructor → A5 collapse accessors → A6 delete `_is_new_sig` +
> `_mechanism_from_legacy_sig`.

### Task A1: Commit the bi-bi exit-gate as a permanent test (cost-bounded subset)

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (add a testset near the existing `bi_bi_rxn` const at line 65).

Deriving all 77 bi-bi rate equations costs ~86 s cold (review-measured) and risks the compile-budget gate. Assert the count (cheap) + derive only a small fixed subset (proves the enumerator routes through the new Sig).

- [ ] **Step 1: Add the test.**
```julia
@testset "bi-bi exit gate: init mechanisms derive (subset)" begin
    mechs = EnzymeRates.init_mechanisms(bi_bi_rxn)
    @test length(mechs) == 77
    # Derive a small subset only — full-77 derivation is ~86s and not worth
    # it on every cold run. Pick the 5 with the fewest enzyme forms.
    by_size = sort(mechs; by = m -> EnzymeRates.n_steps(m))
    for m in by_size[1:5]
        s = EnzymeRates.rate_equation_string(EnzymeRates.compile_mechanism(m))
        @test s isa AbstractString && !isempty(s)
    end
end
```

- [ ] **Step 2: Run it.** Per-file recipe on `test/test_mechanism_enumeration.jl`. Expected: PASS. (`init_mechanisms(bi_bi_rxn)` → 77 was proven 2026-05-26.) Note the added wall-time; if even the 5-subset trips a compile budget, reduce to 2.

- [ ] **Step 3: Commit.**
```bash
git add test/test_mechanism_enumeration.jl
git commit   # "test: bi-bi exit gate — init mechanisms derive (enumerator routes through new Sig)"
```

### Task A2: Port the meaningful legacy validators into `_assert_mechanism_invariants` (EXPLORATORY — drive directly)

**Files:**
- Modify: `src/mechanism_enumeration.jl` — `_assert_mechanism_invariants(m::Mechanism)` (~2122–2143).
- Modify: `test/test_mechanism_enumeration.jl` — add tests for the ported invariants.

Per the spec audit: **port** "every declared **substrate/product** appears in some step" (regulators EXCLUDED — see below) and the kinetic-group "RE/SS-mixing + same-metabolite-per-group" checks; **drop** the rest (moot / superseded). The iso/binding `bound_metabolite` consistency is already checked.

**Why regulators are excluded:** `_drop_unbound_regulators` (types.jl:671) intentionally lets `init_mechanisms` declare a dead-end inhibitor that no step binds yet (`expand_mechanisms` binds it later). Including regulators in the "appears" check would error on every inhibitor `init_mechanisms` output — and `_assert_mechanism_invariants` is called ~80× over enumerated mechanisms. Substrates/products are never dropped, so they are safe to require.

- [ ] **Step 1: Write failing tests** in `test/test_mechanism_enumeration.jl`. First read `src/types.jl:11–160` to confirm the `Species`/`Step`/`Substrate`/`Product` constructor signatures, then:
```julia
@testset "ported invariants: unused substrate/product + kinetic-group composition" begin
    # POSITIVE: an init mechanism with an unbound declared inhibitor must NOT error.
    rxn_inh = @enzyme_reaction begin
        substrates: S[C]; products: P[C]; competitive_inhibitors: R
    end
    for m in EnzymeRates.init_mechanisms(rxn_inh)
        @test EnzymeRates._assert_mechanism_invariants(m) === nothing
    end

    # NEGATIVE 1: a declared SUBSTRATE that no step binds → error.
    rxn_unused = @enzyme_reaction begin
        substrates: S[C], T[C]; products: P[C2]   # T never appears in a step
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
    m_unused = EnzymeRates.Mechanism(rxn_unused, [[s1], [s2], [s3]])
    @test_throws ErrorException EnzymeRates._assert_mechanism_invariants(m_unused)

    # NEGATIVE 2: a kinetic group binding two different metabolites → error.
    rxn2 = @enzyme_reaction begin
        substrates: S[C], A[N]; products: P[CN]
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
    m_mixed = EnzymeRates.Mechanism(rxn2, [[g1a, g1b], [g2]])
    @test_throws ErrorException EnzymeRates._assert_mechanism_invariants(m_mixed)
end
```
**Note:** if `rxn_unused`/the manual `Mechanism` construction itself errors (e.g. an existing constructor invariant rejects an unbound substrate before `_assert_mechanism_invariants` runs), adjust the negative fixture — the point is to reach the ported check.

- [ ] **Step 2: Run — confirm the NEGATIVE cases FAIL** (current `_assert_mechanism_invariants` doesn't check these; the positive case already passes).

- [ ] **Step 3: Port the checks** into `_assert_mechanism_invariants(m::Mechanism)` (after the existing per-step loop, before `nothing`). `substrates(reaction(m))`/`products(reaction(m))` return `Vector{Substrate}`/`Vector{Product}` (types.jl:312–315); `name(met)` gives the Symbol (23–26):
```julia
    # Every declared substrate/product must appear in some step. Regulators
    # are intentionally excluded — _drop_unbound_regulators (types.jl) lets
    # init_mechanisms declare a dead-end inhibitor no step binds yet.
    appearing = Set{Symbol}()
    for s in flat
        for sp in (from_species(s), to_species(s))
            for met in bound(sp); push!(appearing, name(met)); end
        end
        bm = bound_metabolite(s)
        bm === nothing || push!(appearing, name(bm))
    end
    for met in (substrates(reaction(m))..., products(reaction(m))...)
        name(met) in appearing ||
            error("declared substrate/product $(name(met)) appears in no step")
    end

    # Within a kinetic group of size > 1: all binding same metabolite, no RE/SS mix.
    for group in m.steps
        length(group) == 1 && continue
        kinds = [(is_equilibrium(s), bound_metabolite(s)) for s in group
                 if bound_metabolite(s) !== nothing]
        isempty(kinds) && continue
        first_eq, first_met = kinds[1]
        for (eq, met) in kinds[2:end]
            eq == first_eq || error("kinetic group mixes RE and SS binding steps")
            met == first_met ||
                error("kinetic group binds different metabolites: " *
                      "$(name(first_met)) and $(name(met))")
        end
    end
```

- [ ] **Step 4: Run the new test — confirm PASS**, then the **full** `test_mechanism_enumeration.jl` (~80 existing `_assert_mechanism_invariants` calls). Expected: PASS. **If ANY existing call now errors, STOP** — a real enumerated mechanism violates a ported invariant → the invariant is wrong for the decomposed world; re-examine (do not loosen the test).

- [ ] **Step 5: Commit.**
```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit   # "refactor: port substrate/product-coverage + kinetic-group-composition checks to _assert_mechanism_invariants"
```

### Task A3: Migrate ALL legacy-Sig constructions in test_types.jl

**Files:**
- Modify: `test/test_types.jl`.
- Modify: `docs/superpowers/refactor-deleted-tests.md` — §2.1 entries.

- [ ] **Step 1: Inventory every legacy-Sig construction.**
```bash
grep -nE "EnzymeMechanism\{\(|EnzymeMechanism\(\(\(|EnzymeMechanism\(mets|EnzymeMechanism\(base_mets|EnzymeMechanism\(.*_mets" test/test_types.jl
```
Known sites (review-verified): direct `EnzymeMechanism{(mets,rxns)}()` at **9, 33, 67, 413**; 2-arg `EnzymeMechanism(((...)))`/`(mets,rxns)` throughout 76–360. Cover ALL of them — especially the testsets at **2–37** ("EnzymeMechanism struct + accessors"), **60–74** ("stoich_matrix"), and the `AllostericEnzymeMechanism` validator block **372–417** + **449**, which the earlier plan draft omitted.

- [ ] **Step 2: Migrate positive construction/accessor tests** (the 2–37 accessor testset, 60–74 stoich, 76–97, 265–279, 334–342, and any allosteric `cm` fixture) to `@enzyme_mechanism` decomposed grammar, preserving each assertion. Example:
```julia
@testset "EnzymeMechanism construction + accessors (decomposed)" begin
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
    # ...port the remaining accessor assertions (reactions/enzyme_forms/
    #    stoich_matrix/enzyme_row_range/...) onto this decomposed `m`.
end
```
For the same-kinetics group test (87–97), use a competitive-inhibitor fixture binding E and E(S) in one group (grammar in `test_dsl.jl`/`test_accessors.jl`) and assert `kinetic_group` equality. For the allosteric validator block (372–417), rebuild `cm`/`cm_bad` via `@enzyme_mechanism` instead of `EnzymeMechanism{(...)}()`.

- [ ] **Step 3: Re-point ported-validator `@test_throws`** (kinetic-group different-metabolite ~348–352; RE/SS mix ~355–359) onto `_assert_mechanism_invariants` over a decomposed `Mechanism` (reuse the A2 fixture shapes).

- [ ] **Step 4: Delete moot/superseded `@test_throws`** with a §2.1 entry **per assertion** (the audit table reasons): stoich "vanish" (~104), iso-group-size>1 (~113 — if not covered by `_assert_mechanism_invariants`'s bound/iso check, this is a *ported* check → Step 3 instead), empty reactions (~289, covered by `isempty`), dup substrate names (~292 — confirm where enforced; if nowhere, note as dropped), zero enzyme (~300, moot), two metabolites (~308, moot), unknown metabolite (~315, moot), net stoich mismatch (~320, superseded), orphan connectivity (~362–370, superseded).

- [ ] **Step 5: Integrity check — manual, because the script can't see assertion-level deletions.**
```bash
bash scripts/check_test_integrity.sh main; echo "EXIT=$?"   # 0 (necessary, not sufficient)
```
`check_test_integrity.sh` counts `@testset` headings, NOT `@test_throws` inside a surviving testset (script line ~47). So a dropped assertion can pass EXIT=0 silently. **Manually confirm a §2.1 entry exists for each `@test_throws` removed in Step 4.** Then run `test_types.jl` + full suite. Expected: all green.

- [ ] **Step 6: Commit.**
```bash
git add test/test_types.jl docs/superpowers/refactor-deleted-tests.md
git commit   # "test: migrate all legacy-Sig constructions to decomposed grammar; drop moot/superseded validators (§2.1)"
```

### Task A4: Delete the 2-arg constructor + its orphaned validators

**Files:**
- Modify: `src/types.jl` — 2-arg `EnzymeMechanism` constructor (~803–888) + `_validate_kinetic_groups`/`_validate_enzyme_connectivity`/`_validate_stoichiometry`/`_pretty_reaction` (~894–971).

- [ ] **Step 1: Confirm zero remaining callers** (A3 removed the test callers; nothing in src calls the 2-arg form — DSL emits the 1-arg lift):
```bash
grep -rnE "EnzymeMechanism\(\(\(|EnzymeMechanism\(mets|EnzymeMechanism\{\(" src/ test/
grep -rn "_validate_kinetic_groups\|_validate_enzyme_connectivity\|_validate_stoichiometry\|_pretty_reaction" src/ test/
```
Expected: only the definitions themselves. If a test still constructs a legacy Sig, return to A3.

- [ ] **Step 2: Delete** the 2-arg `EnzymeMechanism(mets, rxns)` constructor and the four now-orphaned validators. (Keep `_assert_mechanism_invariants` — it's the live decomposed-world checker, now carrying the ported invariants.)

- [ ] **Step 3: Run full suite + gates.** Expected: green (legacy-Sig construction is gone everywhere). Commit:
```bash
git add src/types.jl
git commit   # "refactor: delete the 2-arg legacy EnzymeMechanism constructor + its orphaned validators"
```

### Task A5: Collapse the ~12 `@generated` accessor branches onto the new-Sig body

**Files:**
- Modify: `src/types.jl` — accessors with `if _is_new_sig(Sig)` (~1325, 1338, 1351, 1369, 1470, 1484, 1493, 1500, 1511, 1522, 1545, 1597) + `Mechanism(em)` (~711–714).

Safe now: A3+A4 removed every legacy-Sig producer, so no `EnzymeMechanism{Sig}` with a legacy Sig exists.

- [ ] **Step 1: Collapse each accessor, one at a time.** Each currently looks like:
```julia
@generated function substrates(::EnzymeMechanism{Sig}) where {Sig}
    if _is_new_sig(Sig)
        <NEW BODY>
        return X
    end
    <LEGACY BODY>          # e.g. Sig[1][1]
end
```
Collapse to (drop guard + trailing legacy body, keep the new body):
```julia
@generated function substrates(::EnzymeMechanism{Sig}) where {Sig}
    <NEW BODY>
    return X
end
```
For the ternary at ~1493 (`_is_new_sig(Sig) ? A : B`), keep only `A`. **After EACH accessor, run `test_accessors.jl`** (per-file recipe) — Expected: PASS. A regression = wrong collapse.

- [ ] **Step 2: Simplify `Mechanism(em)`** to `Mechanism(em::EnzymeMechanism{Sig}) where {Sig} = _mechanism_from_sig(Sig)`.

- [ ] **Step 3: Run `test_accessors.jl` + full suite + perf gate.** Expected: all green; `test_rate_equation_performance` 0-alloc still passes. Commit:
```bash
git add src/types.jl
git commit   # "refactor: collapse EnzymeMechanism accessors onto the single decomposed-Sig body"
```

### Task A6: Delete `_is_new_sig` + `_mechanism_from_legacy_sig` + stale comments

**Files:**
- Modify: `src/types.jl` — `_is_new_sig` (~1292–1308) + doc block (~1280–1287), `_mechanism_from_legacy_sig` (~727–795), stale doc comments (~641, 649, 796, 1275, 1291).
- Modify: `src/rate_eq_derivation.jl:189` (doc comment).

- [ ] **Step 1: Confirm both are now unused.**
```bash
grep -rn "_is_new_sig\|_mechanism_from_legacy_sig" src/   # only their own definitions
```
(`_is_new_sig` lost its callers in A5; `_mechanism_from_legacy_sig` lost its only caller when A5 simplified `Mechanism(em)`.)

- [ ] **Step 2: Delete `_is_new_sig` + its doc block, and `_mechanism_from_legacy_sig`.**

- [ ] **Step 3: Fix stale doc comments** referencing `EnzymeMechanism(metabolites, reactions)` (types.jl ~641, 649, 796, 1275, 1291; rate_eq_derivation.jl ~189) — reword to the current single Sig shape, no temporal words. Confirm no `_legacy_*` names remain (`grep -rn "_legacy" src/` → empty; the `_legacy_step_tuple_from_sig` → `_step_tuple_from_sig` rename already landed in commit 388944c).

- [ ] **Step 4: Verify + full suite + all gates.**
```bash
grep -rnE "_is_new_sig|_mechanism_from_legacy_sig|EnzymeMechanism\(metabolites" src/   # empty
bash scripts/check_test_integrity.sh main; echo "EXIT=$?"                              # 0
```
Full suite green; `test_rate_equation_performance` 0-alloc/<100ns; 3 compile budgets; `test_chokepoint.jl`. Commit:
```bash
git add src/types.jl src/rate_eq_derivation.jl
git commit   # "refactor: delete legacy opaque-Sig discriminator + reconstructor — front-end on one struct family"
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

- [ ] **Step 1: Enumerate every concrete-struct field and its accessor.** Read `src/types.jl` struct definitions (lines ~11–415). Review verified **most accessors already exist** — `name` (Metabolite family, 23–26), `is_equilibrium` (184), `ligands` (120), `atoms` (258), `catalytic_multiplicity` (469), plus `bound`/`conformation`/`residual`/`from_species`/`to_species`/`bound_metabolite`/`reaction`/`steps`/`regulatory_sites`. **Genuinely missing — add only these:** `source_idx` (Step), `cat_allo_states` (AllostericMechanism), and `allo_states` on `RegulatorySite` (read as `site.allo_states` 6× in `rate_eq_derivation.jl`). `.cat_steps` maps to existing `steps(m)` (no new accessor). Re-grep to confirm before adding anything.

- [ ] **Step 2: Write a failing test** in `test/test_accessors.jl` for the genuinely-missing accessors only:
```julia
@testset "added field accessors" begin
    st = EnzymeRates.Step(EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E),
                          EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E_S),
                          EnzymeRates.Substrate(:S), true; source_idx = 7)
    @test EnzymeRates.source_idx(st) == 7
    # + cat_allo_states(am) and allo_states(site) on existing fixtures
end
```
Adjust to the confirmed-missing set from Step 1 (do NOT assert `is_equilibrium`/`name`/etc. — those already exist and would pass immediately, defeating the failing-test step).

- [ ] **Step 3: Run — confirm FAIL** (the genuinely-missing accessors are undefined).

- [ ] **Step 4: Add the missing accessors** next to the related existing ones in `src/types.jl`:
```julia
source_idx(s::Step) = s.source_idx
cat_allo_states(m::AllostericMechanism) = m.cat_allo_states
allo_states(site::RegulatorySite) = site.allo_states
```
(Also add `multiplicity(site::RegulatorySite) = site.multiplicity` if a re-grep shows `.multiplicity` is read without an accessor.) Do NOT collide with `name(p::Parameter, m)` / `name(::Type{P}, idx)` — different signatures, fine.

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

- [ ] **Step 4: Docs final pass.** Re-read README + CLAUDE.md for opaque-form / legacy-Sig references invalidated by the deletions. Update continuation spec §11 status to "legacy Sig removed; front end on one struct family; **derivation back-end opaque-tuple representation deferred** (see spec 'Deferred')."

- [ ] **Step 5: Push + PR — HONEST framing.**
```bash
git push -u origin refactor-to-concrete-types-instead-of-symbols
gh pr create --fill
```
PR body must state plainly: this removes the *legacy* Sig and unifies the **front end** (DSL + enumeration + public mechanism surface) onto one struct family; it does **NOT** close §3.1 #1 — the King-Altman/Wegscheider **derivation back-end still runs on regenerated opaque Symbol tuples** (`_species_name_from_sig`/`_step_tuple_from_sig` → `reactions`/`enzyme_forms`), which remains a deliberate second representation. Honest LOC numbers (branch lands above main — the cost of structured types). Link the three deferred next steps: **(1) derivation back-end unification** (the real remaining work), (2) parameter-naming refactor (continuation spec §7), (3) `_n_fit_params_estimate` undercount. End with the Claude Code attribution footer.

---

## Self-Review Notes

- **Spec coverage:** B → Task B1; A → Tasks A1–A6; C → Task C1; D → Tasks D1, D2–D8; E → Tasks E1–E8; Final → Task F1. Bucket B′ dropped (guard is load-bearing). All remaining buckets covered.
- **Bucket A ordering (review-corrected):** migrate tests (A3) + delete 2-arg constructor (A4) BEFORE collapsing accessors (A5) — `test_types.jl` calls accessors on legacy-Sig mechanisms, so collapsing first would red the suite mid-bucket.
- **Riskiest tasks:** A2 (exploratory validator port — drive directly; STOP if a real enumerated mechanism fails a ported invariant), A5 (12 accessor collapses — one at a time, `test_accessors.jl` after each), E on `types.jl` (struct ordering — last, smallest moves).
- **Hard gates with stop-conditions:** A2 Step 4 (no enumerated mechanism may fail a ported invariant — and regulators are excluded precisely so inhibitor init mechanisms don't), A4 Step 1 (zero legacy-Sig producers before deleting the constructor), A5 (no legacy Sig exists before collapsing accessors).
- **`_StepSideTerm.sym` field stays** (B Task B1 Step 4) — it's live on the emission path; only the `:call` synthesis is dead.
- **Integrity gate is necessary-not-sufficient for A3** — the script counts `@testset`, not `@test_throws` inside one; manually confirm a §2.1 entry per dropped assertion.
- **Honest scope:** this is front-end unification + legacy-Sig removal, NOT full "one struct family." The PR must say so.
- **Don't trust `Pkg.test` exit codes** — read the `Test Summary` line.
- **Line numbers drift** as deletions land — each task re-greps/re-reads before editing.
