# Dead-inactive `:OnlyA` catalysis — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `_expand_to_allosteric` emit, for every mechanism, all valid dead-inactive `:OnlyA`-binding combinations — each with **all** chemical (iso) steps `:OnlyA` — and never a partial-catalysis form; so the ping-pong `:OnlyA` sink/crash class is never enumerated.

**Architecture:** All work is in `src/mechanism_enumeration.jl` plus tests. `_expand_to_allosteric` is rewritten to emit dead-inactive combos directly; `_expand_change_allo_state` gains a partial-catalysis filter; `_expand_promote_catalytic_to_onlya` is dropped contingent on a measured coverage diff; `_valid_onlya_completions` is removed if unused. No `@generated` codegen changes, so the `rate_equation` contract is untouched.

**Tech Stack:** Julia 1.12, EnzymeRates.jl. Tests via `TestEnv` per-file drivers.

**Spec:** `docs/superpowers/specs/2026-07-17-dead-inactive-onlya-catalysis-design.md`. Read it before Task 1. This plan builds on the err1 residual-island stranding fix already on this branch's parent `mwc-targeted-fixes` (commit `29eee7e`) — that fix is load-bearing for the dead-inactive form's clean `d_free_I = 1`.

## Global Constraints

- 92-character lines, 4-space indentation. Match surrounding style.
- `rate_equation` MUST stay allocation-free and < 120 ns (`test_rate_equation_performance`). No task changes codegen, so this gate must stay green **unchanged**. If it moves, STOP.
- `test/reference/allosteric_golden_reference.txt` must stay **byte-identical**. If a block moves, STOP and report — do not regenerate. (The static specs contain no allosteric ping-pong-with-`:OnlyA`, so it is expected unchanged.)
- All `Parameter → Symbol` rendering flows through `name(p, m)`. No `Symbol("K…")`/`Symbol("k…")`/`Symbol("V…")`/`Symbol("L…")` literals in `src/` (AST-walker test at `test/test_types.jl`).
- All 12+ `test/allosteric_ground_truth.jl` gates stay green.
- Do NOT edit a test to make it pass. A failing assertion is a finding to report.
- Run Julia **FOREGROUND and WAIT** (600000 ms timeout). Never background-and-yield. ONE julia process at a time (~5 GiB free / 4 cores). Do NOT run full `Pkg.test()` during tasks — it runs once, in the final task.
- A `@test` failure does NOT abort an `include`. **Grep the output for a non-zero `Fail`/`Error` column**, not merely the absence of `ERROR`.
- Commit after each task. Never skip a pre-commit hook.

## Per-file test drivers

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
# ground truth:
julia --project -e 'using TestEnv; TestEnv.activate(); include("test/allosteric_ground_truth.jl")' 2>&1 | tail -25
# enumeration:
julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates; include("test/test_mechanism_enumeration.jl")' 2>&1 | tail -25
# derivation (shared defs first) — golden + perf gate live here:
julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_rate_eq_derivation.jl")' 2>&1 | tail -25
```

---

## Background the implementer needs

**The rule.** An `:OnlyA` catalytic **binding** (the inactive conformation cannot bind that metabolite) means the inactive conformation cannot complete the catalytic cycle, so it is **catalytically dead**: every chemical (iso) step is `:OnlyA`. The current enumeration instead promotes the *fewest* chemical steps to satisfy the Haldane (one), producing a partial form that sinks and crashes `_kcat_forward` for multi-chemical-step (ping-pong) mechanisms.

**Group model.** `steps(m)` / `steps(am)` is a vector of catalytic **kinetic groups**; `cat_allo_states` is one tag per group. A group is a **binding** if its representative step binds a metabolite (`!is_iso(rep_step(m, g))` / `!is_iso(g[1])`), an **iso/chemical** step otherwise (`is_iso`). Dead-inactive = every iso group `:OnlyA`.

**`:OnlyA` ⟹ dead cycle, verified.** Measured over the PFK-P ping-pong: `_valid_onlya_completions` (`:1776`) returns four size-1 completions (2 partial-chem sinks + 2 opposing-binding), and never the dead-inactive form (both chem `:OnlyA`), which is size-2. The dead-inactive form derives clean today: `d_free_I = 1`, `kcat = 1.3`, no sink, v = 0 at equilibrium.

**Combinatorics.** With all chemical steps pinned `:OnlyA`, the only free variable is the binding subset, so the count is `2^B` (B = binding groups): measured 15 for ping-pong bi-bi (all valid), 3 for uni-uni.

---

### Task 1: Pin the dead-inactive ping-pong derivation (ground-truth gate)

The enumeration will produce dead-inactive ping-pong `:OnlyA` mechanisms; this gate pins that the derivation handles one correctly. It **passes on current code** (the derivation is already correct via the merged stranding fix) — it is the acceptance target the enumeration must feed.

**CORRECTED shape (a dead-inactive form needs pure ISO steps).** The 4-step lumped ping-pong (`E+A↔E(A)↔F+P; F+B↔F(B)↔E+Q`) has **no iso steps** — its chemical conversions release a product, so all four groups are binding/release — and the guard rejects every non-empty `:OnlyA` subset there. A dead-inactive form only exists on a shape with **separate pure isomerisation steps**: the 6-step PFK topology
```
E + ATP ⇌ E(ATP)                                                       :: <t>
E(ATP) <--> E(F16BP; residual = ATP - F16BP)                           :: <t>   iso (chem1)
E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP) :: <t>
E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP)     :: <t>
E(F6P; residual = ATP - F16BP) ⇌ E(ADP)                                :: <t>   iso (chem2)
E + ADP ⇌ E(ADP)                                                       :: <t>
```
Dead-inactive = a **substrate** binding (F6P) `:OnlyA` + **both** iso steps (chem1, chem2) `:OnlyA`, other bindings `:EqualAI`. Measured on current code: guard admits it, `d_free_I = 1`, `rate = 0.1045`, `kcat = 1.3`, no sink.

**Files:**
- Modify: `test/allosteric_ground_truth.jl` (append).

- [ ] **Step 1: Add the derivation gate** for the 6-step dead-inactive mechanism above. Assert:
  - `EnzymeRates._onlya_haldane_violation(...) === nothing` (valid),
  - `_state_rate_polys(am, :I)`'s `d_free_I` is the constant `1`,
  - `rate_equation` finite at a normal point **and** as a product → 0 (no sink),
  - `_kcat_forward` finite and `> 0`,
  - `v = 0` at the equilibrium metabolite ratio (`ADP·F16BP = Keq·ATP·F6P`),
  - **L = 0 cross-check:** `rate_equation(L=0)` equals the active-cycle rate. Build the same 6 steps as a non-allosteric `@enzyme_mechanism` and compare its `rate_equation` to the allosteric one at `L = 0` over ~6 random points (rtol 1e-4). If `@enzyme_mechanism` cannot parse the residual syntax, fall back to comparing `rate_equation(L=0)` against `E_total · N_A / D_A` built from `_state_rate_polys(am, :A)` (a within-derivation consistency check). Document the param mapping in a comment.

  A bespoke two-conformation mass-action oracle for this 6-step shape is **optional** — add one only if it self-validates cleanly (`L=0` reduction + `v=0` at equilibrium); do not block on the from-scratch ping-pong Haldane, which has burned prior attempts twice.

- [ ] **Step 2: Run; confirm green on current code.** Ground-truth driver: all gates pass (this one and the existing 12+). If the derivation gate fails, STOP and report (that would contradict the measured `d_free_I = 1`).

- [ ] **Step 5: Commit.**
```
git add test/allosteric_ground_truth.jl
git commit -m "Gate the dead-inactive ping-pong :OnlyA derivation"
```

---

### Task 2: Rewrite `_expand_to_allosteric` to emit all dead-inactive combos

**Files:**
- Modify: `src/mechanism_enumeration.jl` — `_expand_to_allosteric(m::Mechanism, rxn)` (`:1834-1863`) and its docstring (`:1795-1833`).
- Test: `test/test_mechanism_enumeration.jl` (new testset).

**Interfaces:**
- Produces: `_expand_to_allosteric(m::Mechanism, rxn)` returns, per multiplicity, every non-empty binding-subset dead-inactive form (bare) plus the all-chem-`:OnlyA` V-type (regulator-paired). No partial-catalysis form.

- [ ] **Step 1: Write the failing test.** Append to `test/test_mechanism_enumeration.jl`:
```julia
@testset ":OnlyA to_allosteric emits dead-inactive combos only" begin
    # A ping-pong bi-bi (2 chemical steps). Build the non-allosteric mechanism,
    # allostericise it, and require every :OnlyA-tagged child to have BOTH
    # chemical steps :OnlyA (no partial-catalysis form), and the count of bare
    # (regulator-free) :OnlyA children to equal the 15 valid binding subsets.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) <--> F + P
            F + B <--> F(B)
            F(B) <--> E + Q
        end
    end
    rxn = EnzymeRates.reaction(m)
    kids = EnzymeRates._expand_to_allosteric(m, rxn)
    onlya = [k for k in kids if any(==(:OnlyA), EnzymeRates.cat_allo_states(k))]
    isochem(k) = [g for g in eachindex(EnzymeRates.steps(k))
                  if EnzymeRates.is_iso(EnzymeRates.steps(k)[g][1])]
    # every :OnlyA child is dead-inactive: all chemical steps :OnlyA
    for k in onlya
        @test all(EnzymeRates.cat_allo_states(k)[g] === :OnlyA for g in isochem(k))
    end
    bare = [k for k in onlya if isempty(EnzymeRates.regulatory_sites(k))]
    @test length(bare) == 15
end
```
(If this reaction's `allowed_catalytic_multiplicities` yields more than one multiplicity, scope the `== 15` assertion to `catalytic_multiplicity(k) == 1`, or assert `length(bare) == 15 * n_multiplicities`; measure and use the exact number.)

- [ ] **Step 2: Run it; see it fail** (current emits partial forms and a different count). Enumeration driver.

- [ ] **Step 3: Implement.** Replace the body of `_expand_to_allosteric(m::Mechanism, rxn)`:
```julia
function _expand_to_allosteric(m::Mechanism, rxn::EnzymeReaction)
    n_g = length(steps(m))
    iso = [g for g in 1:n_g if is_iso(rep_step(m, g))]
    bind = [g for g in 1:n_g if !is_iso(rep_step(m, g))]
    regs = Symbol[]
    for rm in regulators(rxn)
        reg = regulator(rm)
        reg isa AllostericRegulator && push!(regs, name(reg))
    end
    sort!(regs)
    results = AllostericMechanism[]
    for cn in allowed_catalytic_multiplicities(rxn)
        # K-type: every non-empty subset of binding groups :OnlyA, with all
        # chemical steps :OnlyA — a catalytically-dead inactive conformation. The
        # :OnlyA binding's metabolite reveals L, so each is emitted bare.
        for mask in 1:(2^length(bind) - 1)
            tags = Symbol[:EqualAI for _ in 1:n_g]
            for g in iso
                tags[g] = :OnlyA
            end
            for (i, g) in enumerate(bind)
                (mask >> (i - 1)) & 1 == 1 && (tags[g] = :OnlyA)
            end
            _onlya_haldane_violation(rxn, steps(m), tags) === nothing || continue
            push!(results, AllostericMechanism(
                reaction(m), copy(steps(m)), tags, cn, RegulatorySite[]))
        end
        # V-type: no :OnlyA binding, all chemical steps :OnlyA. L folds into kcat
        # and is unobservable, so emit only paired with a declared regulator.
        if !isempty(iso) && !isempty(regs)
            vtags = Symbol[:EqualAI for _ in 1:n_g]
            for g in iso
                vtags[g] = :OnlyA
            end
            am_cat = AllostericMechanism(
                reaction(m), copy(steps(m)), vtags, cn, RegulatorySite[])
            for reg in regs, tag in (:OnlyA, :OnlyI)
                push!(results, _make_am_with_added_reg(am_cat, reg, tag, 0))
            end
        end
    end
    unique!(results)
end
```
Rewrite the docstring (`:1795-1833`) to describe the present behaviour: dead-inactive combos (all chemical steps `:OnlyA`), one per valid non-empty binding subset (K-type, bare) plus the all-chem-`:OnlyA` V-type (regulator-paired); the all-`:EqualAI` baseline is never emitted. Describe the code as it is — no "previously"/"changed"/temporal language.

- [ ] **Step 4: Run; the new testset passes.** Other enumeration tests WILL break (counts change) — record which; do NOT fix them here (Task 5 re-measures). `_valid_onlya_completions` is now unreferenced by this function but still used by `_expand_promote_catalytic_to_onlya` — leave it until Task 4.

- [ ] **Step 5: Commit.**
```
git commit -am "Emit dead-inactive :OnlyA-binding combos from to_allosteric"
```

---

### Task 3: `_partial_onlya_catalysis` predicate + `_expand_change_allo_state` filter

Without a constructor guard (decision 4), `_expand_change_allo_state` is the only place a relaxation could re-create a partial-catalysis form (relaxing one chemical step of a dead-inactive `:OnlyA`-binding mechanism = the err2 sink). It must filter those.

**Files:**
- Modify: `src/mechanism_enumeration.jl` — add `_partial_onlya_catalysis`; filter in `_expand_change_allo_state` (`:2013-2038`).
- Test: `test/test_mechanism_enumeration.jl`.

**Interfaces:**
- Produces: `_partial_onlya_catalysis(cat_steps, cat_allo_states) → Bool` — true when an `:OnlyA` binding group coexists with a non-`:OnlyA` iso group.

- [ ] **Step 1: Write the failing test.**
```julia
@testset "change_allo_state drops partial-catalysis relaxations" begin
    # Dead-inactive ping-pong: a substrate binding + both chemical steps :OnlyA.
    dead = @allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A <--> E(A)   :: OnlyA
            E(A) <--> F + P   :: OnlyA
            F + B <--> F(B)   :: EqualAI
            F(B) <--> E + Q   :: OnlyA
        end
    end
    am = EnzymeRates.AllostericMechanism(dead)
    kids = EnzymeRates._expand_change_allo_state(am)
    isochem(k) = [g for g in eachindex(EnzymeRates.steps(k))
                  if EnzymeRates.is_iso(EnzymeRates.steps(k)[g][1])]
    hasonlyabind(k) = any(EnzymeRates.cat_allo_states(k)[g] === :OnlyA &&
                          !EnzymeRates.is_iso(EnzymeRates.steps(k)[g][1])
                          for g in eachindex(EnzymeRates.steps(k)))
    # no child may be a partial: :OnlyA binding present but some chem not :OnlyA
    for k in kids
        @test !(hasonlyabind(k) &&
                !all(EnzymeRates.cat_allo_states(k)[g] === :OnlyA for g in isochem(k)))
    end
end
```
(Adjust the mechanism if that exact tag set is not constructable on the current tree — pick a constructable dead-inactive ping-pong; the invariant asserted is what matters.)

- [ ] **Step 2: Run it; see it fail** (a chem-step relaxation currently produces a partial).

- [ ] **Step 3: Implement.** Add above `_expand_change_allo_state`:
```julia
"""
    _partial_onlya_catalysis(cat_steps, cat_allo_states) → Bool

True when an `:OnlyA` catalytic **binding** group coexists with an
isomerisation (chemical) group not tagged `:OnlyA` — the structural signature
of a live-catalysis inactive conformation under an `:OnlyA` binding, which
cannot complete the catalytic cycle and produces a kinetic sink. Used by the
enumeration moves to avoid generating such a form.
"""
function _partial_onlya_catalysis(cat_steps::Vector{Vector{Step}},
                                  cat_allo_states::Vector{Symbol})
    has_onlya_binding = any(cat_allo_states[g] === :OnlyA &&
                            is_binding(cat_steps[g][1]) for g in eachindex(cat_steps))
    has_onlya_binding || return false
    any(is_iso(cat_steps[g][1]) && cat_allo_states[g] !== :OnlyA
        for g in eachindex(cat_steps))
end
```
In `_expand_change_allo_state`, in the catalytic-group loop where a relaxation's `new_states` is validated against `_onlya_haldane_violation` (`:2020-2021`), also skip when `_partial_onlya_catalysis(steps(am), new_states)` is true:
```julia
        _onlya_haldane_violation(reaction(am), steps(am), new_states) ===
            nothing || continue
        _partial_onlya_catalysis(steps(am), new_states) && continue
```

- [ ] **Step 4: Run; the new testset passes; the existing `change_allo_state` tests still pass.** Enumeration driver — grep for non-zero Fail/Error.

- [ ] **Step 5: Commit.**
```
git commit -am "Filter partial-catalysis relaxations in change_allo_state"
```

---

### Task 4: Resolve the promote move (coverage diff) + remove dead completion code

**Files:**
- Modify: `src/mechanism_enumeration.jl` — `_expand_promote_catalytic_to_onlya` (`:2074-2094`), `_add_expansions_mech!` (`:2263-2273`), and possibly `_valid_onlya_completions` / `_each_subset`.
- Create: `scratchpad/coverage_diff.jl` (measurement, not committed to `src`/`test`).

- [ ] **Step 1: Measure the coverage diff.** Write a script that enumerates several reactions (uni-uni; ordered bi-bi; the ping-pong `A,B→P,Q`; and one ter-substrate) to a fixed BFS depth (3–4) **two ways**: (a) the current move set, and (b) the same but with `_expand_promote_catalytic_to_onlya` removed from `_add_expansions_mech!`. For each, collect the set of enumerated `AllostericMechanism`s carrying any `:OnlyA` tag, keyed by `string(k)` (or `eq_hash` if available). Report the set difference **both directions** and the ter-substrate combo count (to check `2^B` tractability). Run it FOREGROUND.

- [ ] **Step 2: Decide from the measurement.**
  - **If (b) loses no `:OnlyA` mechanism** vs (a): remove `_expand_promote_catalytic_to_onlya` (both methods, `:2074-2094`) and its call in `_add_expansions_mech!` (`:2269` area). Record the diff-empty result.
  - **If (b) loses any mechanism:** KEEP `_expand_promote_catalytic_to_onlya`, but change its completion so it never yields a partial — replace its `_valid_onlya_completions` call with the dead-inactive construction (promote a binding ⟹ that binding + all iso groups `:OnlyA`; a promoted iso group ⟹ all iso groups `:OnlyA`). Record which mechanisms needed it.

- [ ] **Step 3: Remove now-dead completion code.** Grep for `_valid_onlya_completions` and `_each_subset`. If no caller remains (to_allosteric no longer uses them after Task 2; promote either removed or switched to the dead-inactive construction), delete both functions and their docstrings. If a caller remains, leave them.

- [ ] **Step 4: Run the enumeration + ground-truth drivers.** Enumeration counts will differ from the pre-change baseline (expected); the ground-truth gates and the Task-1 gate stay green. Grep for non-zero Fail/Error. Do not update `expected_n_*` counts here — Task 5 does the migration.

- [ ] **Step 5: Commit.**
```
git commit -am "Resolve promote move via coverage diff; drop dead completion code"
```
Put the diff outcome (dropped, or kept + which mechanisms) in the commit body.

---

### Task 5: Migration and full-suite verification

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` and/or `test/test_mechanism_enumeration.jl` — updated `expected_n_*` counts for reactions whose enumeration changed.

- [ ] **Step 1: End-to-end no-crash check.** Write and run (FOREGROUND) a scratch script that enumerates the ping-pong `A,B→P,Q` (and a ter-substrate) to depth 3–4 with the full move set and calls `_kcat_forward` on every enumerated `AllostericMechanism`. Assert **zero** "no kcat components" crashes and **zero** `NaN`/non-finite `rate_equation` at a saturating-substrate, zero-product point. This is the design's ultimate goal (the 112-class no longer enumerated) — it must hold. If any mechanism still crashes, characterise it (its tag vector) and STOP: a surviving partial means a move still produces one (likely the `:NonequalAI`-via-V-type-relaxation path — see Risks).

- [ ] **Step 2: Re-measure enumeration counts.** Run the enumeration driver; for each failing count assertion, re-measure the actual value and update `expected_n_*` to it, with a one-line comment/commit-note explaining the change as the dead-inactive rule (ping-pong `:OnlyA` reactions produce the dead-inactive combo family; single-chem reactions unchanged up to rate-equivalent dedup). A count change on a **single-chemical-step** reaction that is NOT explained by rate-equivalent dedup is a finding — STOP and report.

- [ ] **Step 3: Golden byte-identical.** Confirm `git diff --stat <branch-base> HEAD -- test/reference/allosteric_golden_reference.txt` is empty. If it moved, STOP and explain (no static spec should be an allosteric ping-pong-with-`:OnlyA`).

- [ ] **Step 4: Full suite (controller-run, SOLO).**
```
pgrep -af runtests.jl && pkill -9 -f runtests.jl
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -40
```
Expect 0 failures / 0 errors. Confirm explicitly in the output: `test_rate_equation_performance` green (0 alloc, <120 ns) and the allosteric golden byte-identical.

- [ ] **Step 5: Report.** LOC delta (`git diff --stat`), the promote-drop outcome, the enumeration-count changes with the rule as explanation, and the end-to-end no-crash result.

- [ ] **Step 6: Commit.**
```
git commit -am "Migrate enumeration counts; verify dead-inactive :OnlyA end-to-end"
```

## Self-Review

- **Spec coverage:** decision 1/2 (dead-inactive, all chem) → Task 2; decision 2 (all combos) → Task 2; V-type all-chem-`:OnlyA` → Task 2; decision 3 (promote drop, coverage-gated) → Task 4; decision 4 (enumeration-only, no guard) → the `_partial_onlya_catalysis` predicate is used only by the moves (Task 3), never the constructor; change_allo_state filter → Task 3; derivation soundness → Task 1; migration + end-to-end no-crash → Task 5. Covered.
- **Ordering:** pin the derivation target (Task 1, passes now) → core emission rewrite (Task 2) → partial filter (Task 3) → promote decision + dead-code removal (Task 4) → migration + full suite (Task 5). Nothing is deleted while referenced (`_valid_onlya_completions` removed only in Task 4, after Task 2 stops using it and Task 4 resolves promote).
- **Type consistency:** `_partial_onlya_catalysis(cat_steps::Vector{Vector{Step}}, cat_allo_states::Vector{Symbol})` defined and used in Tasks 3–4; `pingpong_dead_inactive_flux` defined and used in Task 1; the rewritten `_expand_to_allosteric` keeps its `(m::Mechanism, rxn)` signature and the `AllostericMechanism(reaction, steps, tags, cn, RegulatorySite[])` / `_make_am_with_added_reg(am_cat, reg, tag, 0)` calls it already uses.
- **Key risk (Task 5 Step 1):** the new all-chem-`:OnlyA` V-type changes how a productive-inactive `:NonequalAI` ping-pong is reached (relaxing one V-type chem step now leaves the other `:OnlyA` — a partial with no `:OnlyA` binding, which the Task-3 filter does not catch because it is keyed on `:OnlyA` bindings). The end-to-end no-crash check is the guard against that producing a crashing mechanism; if it fires, the fix is either to extend the partial filter to the no-`:OnlyA`-binding V-type case or to reach both-chem-`:NonequalAI` without the partial intermediate — decide from the measured failure.
