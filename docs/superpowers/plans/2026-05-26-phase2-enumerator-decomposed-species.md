# Phase 2 — Enumerator to decomposed Species: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (Denis drives directly per spec §4; do NOT delegate Task 1 to unsupervised subagents — it is exploratory, continuation spec §11.4). Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `init_mechanisms`/`expand_mechanisms` produce decomposed `Species` (`Species([Substrate(:A)], :E)`) instead of opaque (`Species([], :E_A)`), then delete the `_RawSpec` round-trip and the form-name parse/synthesis helpers — closing continuation-spec §3.1 #1 (one struct family) and #2 (no form-name parse-back).

**Architecture:** The topology backtracker (`_catalytic_topologies`/`backtrack!`) already builds decomposed `Vector{Vector{Step}}`; the decomposition is thrown away by a round-trip through `_RawSpec` (`_mechanism_spec_from_steps` at the tail of `_catalytic_topologies`, then dead-end + grouping helpers operate on opaque `_RawSpec`, then `EnzymeMechanism(spec::_RawSpec)` rebuilds Species opaquely). We keep the decomposed Steps flowing through dead-end enumeration and equivalence grouping into `Mechanism` directly. Because `name(decomposed) == name(opaque)` (spec §2 diagnostic 3), form names are byte-identical and the refactor is behavior-preserving, gated by enumeration counts + derivation regression.

**Tech Stack:** Julia 1.10+, Test/Aqua/JET. Source: `src/mechanism_enumeration.jl` (2,799 LOC). Spec: [`../specs/2026-05-26-phase2-enumerator-decomposed-species-design.md`](../specs/2026-05-26-phase2-enumerator-decomposed-species-design.md).

**Base:** branch `refactor-to-concrete-types-instead-of-symbols`, tip `08e54b3` (spec commit). Phase 1 green gate confirmed.

---

## Conventions (every commit)

- **Per-file test recipe** (resolve deps once, then iterate):
  ```bash
  julia --project=. -e '
    using Pkg; Pkg.activate(temp=true); Pkg.develop(path=".")
    Pkg.add(["Test","Aqua","JET","OptimizationBBO","OptimizationPyCMA","OrdinaryDiffEqFIRK","Tables","DataFrames","Statistics","Optimization","Random","CSV"])
    using Test, EnzymeRates
    include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<file>.jl")'
  ```
- **Per-commit gates** (build their own small mechanisms; run every commit):
  - `test/test_rate_eq_derivation.jl` → `test_rate_equation_performance` (0-alloc/<100ns)
  - `test/test_compile_budget.jl` (3 budget gates)
  - `test/test_chokepoint.jl`
- **Full suite:** `julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/out.log 2>&1; echo "EXIT=$?"` then read the `Test Summary` line (the trailing echo swallows `Pkg.test`'s exit — trust the summary).
- **Test integrity:** `bash scripts/check_test_integrity.sh main; echo "EXIT=$?"` must be 0 (no pipe). Any helper-test deletion needs a `docs/superpowers/refactor-deleted-tests.md` §2.1 entry in the same commit.
- **No `--amend`.** Evergreen comments (no "Stage N"/"previously"/"legacy"). Commit footer: `src delta: -X / +Y net Z, cumulative: ±W` (`wc -l src/*.jl`, cumulative vs main 7,136) + `Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>`.

---

## Task 0: Lock the regression baseline

**Why:** every later task is behavior-preserving. The substrate/product dead-end path Task 1 rewrites lives inside `init_mechanisms` and only fires for **multi-substrate** reactions (uni-uni has no substrate/product dead-ends — binding the single sub + single prod is excluded as "all subs and all prods"). The dead-end-bearing mechanisms also tend to be the *larger* forms that hit the documented compile-cost wall (spec §2), so a *derivation*-based golden can't cover them. The right anchor is a **structural** golden: a canonical, derivation-free snapshot of every `init_mechanisms` mechanism's steps (form names + bound metabolite + RE/SS + kinetic grouping) for bi-bi. Because `name(decomposed) == name(opaque)`, this snapshot is byte-identical before/after iff the rewrite preserves behavior — and it is cheap (no `@generated`).

The existing `bi-bi=11`/`ter-ter=283` count testsets guard *topology* generation (`_catalytic_topologies`), which Task 1 does not change; the new structural golden guards the *dead-end + grouping* output, which it does.

**Files:**
- Test: `test/test_mechanism_enumeration.jl` (add one testset)

- [ ] **Step 1: Confirm the existing count invariants are present and green; identify the `init_mechanisms` count assertion.**

Run:
```bash
grep -n "== 11\|== 283\|init_mechanisms\|enumerate_all" test/test_mechanism_enumeration.jl | head -30
```
Expected: testsets asserting the topology counts (`11`/`283`) and at least one asserting an `init_mechanisms`/`enumerate_all` count. Note the line numbers — the topology counts + the `init_mechanisms` count are part of the Task 1/2/3 gate alongside the structural golden.

- [ ] **Step 2: Write a structural (derivation-free) golden test over bi-bi `init_mechanisms`.**

Add to `test/test_mechanism_enumeration.jl` near the integration testsets. `bi_bi_rxn` is already defined at the top of the file. The test asserts the current structure equals a committed fixture (`test/fixtures/phase2_init_golden.txt`); on first run, if the fixture is absent it writes it (bootstrap), then a human commits it. After bootstrap, the `@test keys == golden` line is the permanent gate that fails loudly if Task 1–3 changes structure:

```julia
@testset "Phase 2 baseline: bi-bi init_mechanisms structural golden" begin
    # Canonical, derivation-free key for one Mechanism: per kinetic group,
    # the sorted set of (from form, to form, bound metabolite, RE/SS) tuples;
    # groups themselves sorted so the key is order-independent.
    function _mech_struct_key(m::EnzymeRates.Mechanism)
        grpkeys = String[]
        for grp in m.steps
            stepkeys = sort([
                string((EnzymeRates.name(EnzymeRates.from_species(s)),
                        EnzymeRates.name(EnzymeRates.to_species(s)),
                        EnzymeRates.bound_metabolite(s) === nothing ? :iso :
                            EnzymeRates.name(EnzymeRates.bound_metabolite(s)),
                        EnzymeRates.is_equilibrium(s)))
                for s in grp])
            push!(grpkeys, join(stepkeys, "|"))
        end
        join(sort(grpkeys), " ;; ")
    end

    init = EnzymeRates.init_mechanisms(bi_bi_rxn)
    keys = sort([_mech_struct_key(m) for m in init])
    @test length(unique(keys)) == length(keys)   # no structural duplicates post-dedup

    fixture = joinpath(@__DIR__, "fixtures", "phase2_init_golden.txt")
    if !isfile(fixture)
        mkpath(dirname(fixture)); write(fixture, join(keys, "\n"))
        @warn "bootstrapped phase2 golden fixture; commit it" fixture
    end
    golden = readlines(fixture)
    @test keys == golden          # ← permanent structural regression gate
    @test length(init) == length(golden)   # init_mechanisms count invariant
end
```

- [ ] **Step 3: Bootstrap + verify the fixture.**

Run the per-file recipe with `<file>` = `test_mechanism_enumeration.jl` twice: first run bootstraps `test/fixtures/phase2_init_golden.txt` (warns); second run must PASS the `keys == golden` assertion against the committed-to-be fixture.
```bash
wc -l test/fixtures/phase2_init_golden.txt   # record line count (= bi-bi init count, ~77)
```

- [ ] **Step 4: Commit the test + fixture.**

```bash
git add test/test_mechanism_enumeration.jl test/fixtures/phase2_init_golden.txt
git commit  # "test: Phase 2 structural golden for bi-bi init_mechanisms"
```
This fixture is a permanent regression test — Tasks 1–3 keep it green; it is NOT deleted in Task 5.

---

## Task 1: SPIKE — decomposed dead-end enumeration + equivalence grouping

**This is the exploratory core (spec §5).** It is NOT pre-written code: the dead-end enumeration (`_expand_substrate_product_dead_ends`, `src/mechanism_enumeration.jl:1133-1271`) and `_apply_equivalence_grouping` (`:1433-1454`) are ~190 lines keyed entirely on opaque form-name Symbols (`Dict{Symbol,Set{Symbol}}`, `_dead_end_form_name`, `_RawStep` reactant/product Symbol vectors). The spike rewrites them to operate on decomposed `Vector{Vector{Step}}`, **still converting to `_RawSpec` at the tail of `_init_raw_specs`** so the rest of the pipeline is untouched and the change is isolated. The four §5 gates are the empirical arbiter.

**Transformation recipe (the structural mapping; verified by the gates, not asserted as final code):**

| Opaque today | Decomposed target |
|---|---|
| `_catalytic_topologies` returns `Vector{_RawSpec}` (converts at `:807` via `_mechanism_spec_from_steps`) | returns `Vector{Vector{Step}}` (the `steps` at `:801`, *not* converted) |
| `bound::Dict{Symbol,Set{Symbol}}` via `_parse_bound` | read `bound(from_species(s))` / `bound(to_species(s))` directly off each Step |
| `cat_forms::Set{Symbol}` = `all_form_names(spec)` | set of `from_species`/`to_species` *Species values* (or their `name`s, since `name(decomposed)` is stable) |
| `_dead_end_form_name(f, fb, m)` → `:E_A_I` | construct `Species(sort([bound(f)...; CompetitiveInhibitor(m)]), conformation(f))`; its `name` renders the same Symbol |
| `_RawStep([cat_form, met], [de_name], true, g)` | `Step(cat_species, de_species, <Substrate/Product/CompetitiveInhibitor>(met), true)` in inner group `g` |
| `_apply_equivalence_grouping` merges `_RawStep`s by `(step_metabolite, RE/SS)` into shared `kinetic_group` | merge inner `Vector{Step}`s by `(bound_metabolite name, is_equilibrium)` — same logic `_expand_split_kinetic_group` already uses on Mechanisms |

**Files:**
- Modify: `src/mechanism_enumeration.jl` — `_catalytic_topologies` tail (`:761-811`), `_expand_substrate_product_dead_ends` (`:1133-1271`), `_apply_equivalence_grouping` (`:1433-1454`), `_init_raw_specs` (`:1404-1421`)
- Test: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Read the three functions end-to-end before touching them.**

Run:
```bash
sed -n '300,648p;740,880p;1119,1271p;1404,1490p' src/mechanism_enumeration.jl
```
Confirm: the dead-end metabolite role (substrate/product/inhibitor) is recoverable at construction time — substrates/products from `reactants(reaction)`, inhibitors from `regulators(reaction)` of type `CompetitiveInhibitor`. This determines which `Metabolite` subtype goes in the decomposed bound list.

- [ ] **Step 2: Make `_catalytic_topologies` return `Vector{Vector{Step}}`.**

Change the tail (`:761-811`): replace `result = _RawSpec[]` with `result = Vector{Vector{Step}}[]`... (actually `Vector{Vector{Step}}` is one topology = a `Vector{Step}`; the collection is `Vector{Vector{Step}}`). At `:807`, push `steps` (the sorted `Vector{Step}` from `:801`) instead of `_mechanism_spec_from_steps(reaction, steps)`. Update the function's return-type docstring.

- [ ] **Step 3: Add a temporary conversion at the `_init_raw_specs` boundary.**

In `_init_raw_specs` (`:1404`), the call `_expand_substrate_product_dead_ends(topos, reaction)` now receives `Vector{Vector{Step}}`. After the dead-end + grouping rewrite (Steps 4–5) produces `Vector{Vector{Vector{Step}}}` (a list of mechanisms, each a `Vector{Vector{Step}}` group-list), convert each to `_RawSpec` here via a small local `_rawspec_from_groups(reaction, groups)` that flattens groups to a `Vector{Step}` with kinetic-group indices and calls the existing `_mechanism_spec_from_steps` machinery. **This keeps the `_RawSpec → _mechanism_from_raw → EnzymeMechanism` tail unchanged** — the isolation that makes this a spike.

- [ ] **Step 4: Rewrite `_expand_substrate_product_dead_ends` to consume/produce decomposed groups.**

Signature → `_expand_substrate_product_dead_ends(topos::Vector{Vector{Step}}, reaction)::Vector{Vector{Vector{Step}}}`. Apply the recipe table: replace every `Dict{Symbol,Set{Symbol}}`/`_dead_end_form_name`/`_RawStep` construction with decomposed-Species/`Step` construction. Keep `_substrate_product_dead_end_opportunities`, `_competition_patterns` logic identical (they reason over metabolite-name sets, which are unchanged) — only the form representation changes. Preserve mirror-step kinetic-group inheritance (`:1236-1262`) by carrying group membership in the group-list structure.

- [ ] **Step 5: Rewrite `_apply_equivalence_grouping` to merge decomposed groups.**

Operate on the group-list: merge inner `Vector{Step}` groups whose representative binding step shares `(name(bound_metabolite), is_equilibrium)`. Preserve the `floor_pc` param-count flooring and `_n_fit_params_estimate_from_steps` semantics (compute the estimate from the merged group structure).

- [ ] **Step 6: GATE — run the spike's four checks.**

Run the per-file recipe with `test_mechanism_enumeration.jl`. Expected, ALL of:
1. Topology counts unchanged: bi-bi `== 11`, ter-ter `== 283`; `init_mechanisms` count unchanged.
2. dedup testsets PASS (struct `==`/`hash` canonicalize identically).
3. The structural-golden testset's `@test keys == golden` PASSES — every init mechanism's decomposed steps render to the same form names + bound metabolite + RE/SS + grouping as the committed (opaque-baseline) fixture.
4. `test/test_compile_budget.jl` `init_mechanisms` budget (750) PASS.

To inspect a mismatch (gate #3 fails), diff the live structure against the fixture:
```bash
julia --project=. -e 'include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); ...' # or re-run the testset and read the @test failure
```
The structural golden uses `name(::Species)`, and `name(decomposed) == name(opaque)` is proven (spec §2) — so a mismatch means the decomposed Step you constructed has a *different bound list* than the opaque form encoded (wrong `Metabolite` subtype, missing/extra metabolite, wrong conformation): a real regression. Fix the Step construction; do **not** edit the fixture. **If the opaque rep proves load-bearing in a way Steps can't express → STOP and escalate to Denis (spec §5 clause).**

- [ ] **Step 6b: Derivation cross-check on the derivable subset.**

The structural golden proves the *mechanisms* are identical; this confirms derivation of the derivable ones is unchanged. The existing derivation testsets in the full suite (`test_rate_eq_derivation.jl`, `test_identify_rate_equation.jl`) already pin many mechanisms — running the full suite (Step 7) is the derivation cross-check. No separate fixture needed.

- [ ] **Step 7: Run per-commit gates + full suite.**

Run `test_rate_equation_performance`, the 3 compile budgets, `test_chokepoint.jl`, then the full suite + `check_test_integrity.sh main`. All green.

- [ ] **Step 8: Commit.**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit  # "refactor: decomposed dead-end enumeration + equivalence grouping behind _RawSpec boundary"
```
**← Denis reviews here before Task 2.**

---

## Task 2: Collapse the `_RawSpec` round-trip

**Why:** with the dead-end/grouping rewrite proven, `init_mechanisms` can build `Mechanism` directly from the decomposed group-lists; the entire `_RawSpec` machinery becomes dead.

**Files:**
- Modify: `src/mechanism_enumeration.jl` — `_init_raw_specs`/`init_mechanisms` (`:1404`, `:2716`), delete `_RawSpec`/`_RawStep` (`:13-41`), `_mechanism_spec_from_steps` (`:821`), `_stepspec_from_step` (`:854`), `_mechanism_from_raw` (`:2726`), `_rawspec_from_groups` (Task 1 temp helper), `EnzymeMechanism(spec::_RawSpec)` (`:1308`), `compile_mechanism(spec::_RawSpec)` (`:1286`)
- Test: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Add `Mechanism(reaction, groups::Vector{Vector{Step}})` path (if not already the constructor).**

Confirm `Mechanism(reaction, steps::Vector{Vector{Step}})` exists and accepts the group-list Task 1 produces. Run:
```bash
grep -n "struct Mechanism\|^Mechanism(\|function Mechanism(" src/types.jl
```
If the existing constructor matches, no new code; otherwise add the minimal constructor.

- [ ] **Step 2: Rewrite `_init_raw_specs` → `init_mechanisms` to return `Mechanism` directly.**

Rename/inline `_init_raw_specs` so `init_mechanisms(r::EnzymeReaction)` is: `topos = _catalytic_topologies(r)`; `expanded = _expand_substrate_product_dead_ends(topos, r)`; `[Mechanism(r, _apply_equivalence_grouping(groups, floor_pc)) for groups in expanded]` (using the same `_apply_equivalence_grouping` rewritten in Task 1 Step 5, which now takes and returns group-lists). Remove the `_rawspec_from_groups` temp helper and the `_mechanism_from_raw` hop.

- [ ] **Step 3: GATE — golden + counts.**

Run `test_mechanism_enumeration.jl`; the structural-golden `@test keys == golden` PASSES; topology counts 11/283 + init-count invariant PASS.

- [ ] **Step 4: Delete the now-dead `_RawSpec` machinery.**

Delete `_RawSpec`, `_RawStep` structs; `_mechanism_spec_from_steps`; `_stepspec_from_step`; `_mechanism_from_raw`; `EnzymeMechanism(spec::_RawSpec)`; `compile_mechanism(spec::_RawSpec)`. After each deletion:
```bash
grep -rn "_RawSpec\|_RawStep\|_mechanism_spec_from_steps\|_stepspec_from_step\|_mechanism_from_raw" src/
```
Expected: only `_RawAllostericSpec` references remain (that struct still wraps `_RawSpec`? — if so, defer its deletion to Task 3; verify it no longer references the deleted `_RawSpec`). If any non-allosteric caller remains, it was missed — resolve before continuing.

- [ ] **Step 5: Handle any test that constructed `_RawSpec` directly.**

Run:
```bash
grep -rn "_RawSpec\|_RawStep\|_mechanism_spec_from_steps" test/
```
For each, adapt to the `Mechanism`/`Step` surface (mechanical syntax adaptation per spec §2 — assertions unchanged). If a testset covered *only* a deleted private helper, add a §2.1 deleted-tests log entry in this commit.

- [ ] **Step 6: Per-commit gates + full suite + integrity.** All green.

- [ ] **Step 7: Commit.**

```bash
git add src/mechanism_enumeration.jl test/ docs/superpowers/refactor-deleted-tests.md
git commit  # "refactor: build Mechanism directly in init_mechanisms; delete _RawSpec round-trip"
```

---

## Task 3: Sweep the `_expand_*` moves + delete `_RawAllostericSpec`

**Why:** the `expand_mechanisms` moves still synthesize opaque dead-end names via `_dead_end_form_name` at 4 sites (`:1009`, `:1163`, `:1245-1247`, `:2010`). Once init emits decomposed Mechanisms, the moves must too. `_RawAllostericSpec` + `_raw_from_mechanism` are the allosteric round-trip analog.

**Files:**
- Modify: `src/mechanism_enumeration.jl` — the 4 `_dead_end_form_name` call sites; delete `_RawAllostericSpec` (`:61-98`), `_raw_from_mechanism` (`:1698`), `compile_mechanism(spec::_RawAllostericSpec)` (`:1288`)
- Test: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Inventory the 4 call sites and their context.**

Run:
```bash
grep -n "_dead_end_form_name\|_raw_from_mechanism\|_RawAllostericSpec" src/mechanism_enumeration.jl
```
Three of the four sites (`:1009`, `:1163`, `:1245-1247`) are inside the now-decomposed dead-end enumeration (Task 1 may have already removed them — re-grep). The `:2010` site is inside `_expand_add_dead_end_regulator`'s kernel.

- [ ] **Step 2: Rewrite the `_expand_add_dead_end_regulator` site (`~:2010`) to build a decomposed dead-end Step.**

Replace `de_name = _dead_end_form_name(cf, bound[cf], reg_name)` + the `_RawStep`/opaque-Species construction with a `Step(cat_species, Species(sort([bound(cat_species)...; CompetitiveInhibitor(reg_name)]), conformation(cat_species)), CompetitiveInhibitor(reg_name), true)` in a fresh kinetic group. Read `~:1880-2020` first for the exact surrounding loop variables.

- [ ] **Step 3: GATE — golden + counts + regulator dead-end path.**

Run `test_mechanism_enumeration.jl`; structural-golden `@test keys == golden` PASSES; topology counts 11/283 + init-count invariant; the regulator-bearing enumeration testsets (which exercise `_expand_add_dead_end_regulator`, e.g. `uni_uni_with_reg`) PASS. **Note:** the regulator dead-end path output is *not* in the bi-bi structural golden (that golden covers `init_mechanisms`, which has no regulator dead-ends); its regression guard is the existing `expand_mechanisms`/`uni_uni_with_reg` testsets + the full suite.

- [ ] **Step 4: Delete `_RawAllostericSpec`, `_raw_from_mechanism`, `compile_mechanism(::_RawAllostericSpec)`.**

After deletion:
```bash
grep -rn "_RawAllostericSpec\|_raw_from_mechanism" src/ test/
```
Expected: empty. Adapt/log any test referents per §2/§2.1.

- [ ] **Step 5: Per-commit gates + full suite + integrity.** All green.

- [ ] **Step 6: Commit.**

```bash
git add src/mechanism_enumeration.jl test/
git commit  # "refactor: decomposed dead-end steps in expand moves; delete _RawAllostericSpec round-trip"
```

---

## Task 4: Delete the form-name parse/synthesis helpers

**Why:** with no remaining callers, the opaque-form helpers can go — this is the literal closure of §3.1 #2.

**Files:**
- Modify: `src/mechanism_enumeration.jl` — delete `_dead_end_form_name` (`:962`), `_is_estar_form` (`:972`), `_parse_bound` (nested, `~:912`), `_bound_mets_from_form_name` (`:1755`), `_form_name` (already 0 callers); audit `all_form_names`/`_bound_metabolites_at_forms`

- [ ] **Step 1: Confirm zero callers for each helper.**

```bash
for h in _dead_end_form_name _is_estar_form _parse_bound _bound_mets_from_form_name _form_name; do
  echo "== $h =="; grep -rn "$h" src/ test/ | grep -v "function $h\|^[^:]*:[0-9]*:#"
done
```
Expected: each shows only its own definition (and docstring), no call sites. If a caller remains, it belongs to a prior task that was incompletely swept — fix there first.

- [ ] **Step 2: Delete the helpers + any now-orphaned partners.**

Delete the five helpers. Then check `all_form_names`, `_bound_metabolites_at_forms`, `step_forms`, `step_metabolite` for orphan status:
```bash
for h in all_form_names _bound_metabolites_at_forms; do echo "== $h =="; grep -rn "$h" src/ test/; done
```
Delete any that drop to zero callers; keep those still used by canonicalization/dedup.

- [ ] **Step 3: §2.1 log for any helper that had dedicated unit tests.**

If `test/test_mechanism_enumeration.jl` had a `@testset` targeting a deleted helper by name, add a `docs/superpowers/refactor-deleted-tests.md` entry (test file + testset + commit, deleted helper, replacement path = decomposed Step construction, integration coverage = the golden + count tests).

- [ ] **Step 4: Per-commit gates + full suite + integrity.** All green.

- [ ] **Step 5: Commit.**

```bash
git add src/mechanism_enumeration.jl test/ docs/superpowers/refactor-deleted-tests.md
git commit  # "refactor: delete opaque form-name parse/synthesis helpers (no callers)"
```
**← Denis reviews here before the Final task.**

---

## Task 5: Final — sweep, LOC re-baseline, verification, docs, PR

**Files:**
- Modify: `README.md`, `.claude/CLAUDE.md`, the continuation spec §11 status. (The `test/fixtures/phase2_init_golden.txt` structural golden stays — it is a permanent regression test.)

- [ ] **Step 1: Dead-code re-read of `src/mechanism_enumeration.jl`.** Re-read end-to-end; inline single-use helpers exposed by the deletions; remove now-redundant comments. Confirm no `Species([], :SomeName)` opaque construction remains in `src/`:
```bash
grep -rn "Species(\[\]\|Species(Metabolite\[\]" src/ | grep -iE "_[A-Z]|Estar_"
```
Expected: empty (decomposed construction only).

- [ ] **Step 2: §3.1 success-criteria verification.**
```bash
grep -rn "_RawSpec\|_RawStep\|_RawAllostericSpec\|_mechanism_from_raw\|_dead_end_form_name\|_parse_bound\|_bound_mets_from_form_name\|_is_estar_form\|_form_name" src/  # → empty
bash scripts/check_test_integrity.sh main; echo "EXIT=$?"   # → 0
```
`test_chokepoint.jl` green; `test_rate_equation_performance` 0-alloc/<100ns green; 3 compile budgets green.

- [ ] **Step 3: LOC re-baseline + gate renegotiation (spec §8).**
```bash
wc -l src/*.jl   # record total
```
Record the cumulative delta vs main (7,136). Write the honest achieved number into the PR body and update continuation spec §3.2 / roadmap with the renegotiated measured target (LOC is non-gating). If the number is far above expectation, surface to Denis before the PR.

- [ ] **Step 4: Full suite green; @testset count ≥ main.**
```bash
julia --project=. -e 'using Pkg; Pkg.test()' > /tmp/out.log 2>&1; echo "done"
grep "Test Summary" /tmp/out.log
```

- [ ] **Step 5: Docs.** README DSL examples + CLAUDE.md "Source Layout"/"Key Architecture Decisions" — confirm no opaque-form / `_RawSpec` references survive the deletions. Update continuation spec §11 status to "done — enumerator decomposed".

- [ ] **Step 6: Commit + PR.**
```bash
git status   # confirm no stray scratch files
git add -A
git commit   # "refactor: dead-code sweep + docs; enumerator fully decomposed"
```
Push the branch, open the PR (honest LOC numbers; celebrate the single-struct-family closure; link the deferred parameter-naming refactor as next). End the PR body with the Claude Code attribution footer.

---

## Self-review notes
- **Riskiest task is Task 1** (the spike) — isolated behind the `_RawSpec` boundary, gated by golden-derivation-identity + counts + dedup + compile budget, with the §5 escalation clause.
- **The committed structural golden** (`test/fixtures/phase2_init_golden.txt`) is the cross-task regression anchor; the `@test keys == golden` assertion gates every structural task. It is NOT edited after Task 0 — a mismatch is a regression, not a reason to re-bootstrap.
- **`name(decomposed) == name(opaque)` is proven (spec §2),** so the structural golden should stay byte-identical through Task 1–3; any mismatch points at a wrong bound list in the decomposed Step you constructed.
- **Don't trust `Pkg.test` notification exit codes** — read the summary line.
- **Task 1 is not delegated to unsupervised subagents** (continuation §11.4).
