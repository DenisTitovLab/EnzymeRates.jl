# Catalytic-topology connectivity fix — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `_catalytic_topologies` produce only fully-connected enzyme-form graphs (no dangling single-metabolite forms), and overhaul the enumeration tests to assert the connectivity invariant + specific mechanisms rather than bare counts.

**Architecture:** The backtracker already builds correct, complete catalytic *paths*. The bug is only in the *combiner* that reassembles them: `_steps_for_ordering` cherry-picks individual binding steps by source-accessibility (upper bound only), so free `E` emits a first-binding edge for every metabolite while later second-binding edges are dropped, leaving dangling forms. Fix: build each topology as the **union of whole paths** whose substrate-consumption and product-release orders linearize the chosen weak orderings. Union of complete paths is fully-connected by construction and reads binding history from the path (so ping-pong works, unlike any per-step bound-set rule).

**Tech Stack:** Julia 1.9+; package `EnzymeRates` at `~/.julia/dev/EnzymeRates/`. Run tests with `julia --project -e 'using Pkg; Pkg.test()'` from the package root. Single test files are not run standalone (shared fixtures); use the full suite or a targeted `julia --project -e '…'` snippet.

**Spec:** `docs/superpowers/specs/2026-06-03-catalytic-topology-connectivity-design.md`.

**Invariants you MUST NOT break:** `rate_equation` stays allocation-free and `<100e-9` s/call (`test/test_rate_eq_derivation.jl`) — this fix does not touch derivation. Export count = 18. 92-char lines, 4-space indent. Never skip a pre-commit hook (there are none here).

---

## Key source locations (current tree, `src/mechanism_enumeration.jl`)

- `_catalytic_topologies` — `:164`. Backtracker `:197-504` (correct, untouched). Path dedup + iso-grouping `:506-528`. `_weak_orderings`/`_wo_recurse!` `:531-561` (correct, untouched). **`_steps_for_ordering` `:563-595` (DELETE).** Combiner assembly `:597-674` (REPLACE).
- `is_binding`, `is_iso`, `bound_metabolite`, `from_species`, `to_species`, `name`, `bound`, `conformation`, `residual` are existing accessors.
- `Substrate`, `Product` are `Metabolite` subtypes; `bound_metabolite(step) isa Substrate` classifies a binding step.

---

## Task 0: Baseline + reproduce

**Files:** none (measurement only).

- [ ] **Step 1: Record the suite baseline.**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -5`
Record pass / broken / fail counts and wall time. (Expected at authoring: ~27073 pass / 1 broken / 0 fail — the 1 broken is the deferred Case-B `parameters(Full)` item, out of scope.)

- [ ] **Step 2: Reproduce the bug numbers.** Save `/tmp/repro.jl`:

```julia
using EnzymeRates; const ER=EnzymeRates
boundset(sp)=Set(ER.name(m) for m in ER.bound(sp))
spkey(sp)=(ER.conformation(sp),ER.residual(sp),Tuple(sort(collect(boundset(sp)))))
function viol(steps)
    flat = eltype(steps) <: AbstractVector ? collect(Iterators.flatten(steps)) : steps
    forms=Dict{Any,Any}(); edges=Set{Tuple{Any,Any}}()
    for s in flat
        a,b=ER.from_species(s),ER.to_species(s); forms[spkey(a)]=a; forms[spkey(b)]=b
        push!(edges,(spkey(a),spkey(b))); push!(edges,(spkey(b),spkey(a)))
    end
    v=Tuple{Symbol,Symbol}[]
    for s1 in values(forms), s2 in values(forms)
        (ER.conformation(s1)==ER.conformation(s2) && ER.residual(s1)==ER.residual(s2)) || continue
        b1,b2=boundset(s1),boundset(s2)
        if length(b2)==length(b1)+1 && issubset(b1,b2) && !((spkey(s1),spkey(s2)) in edges)
            push!(v,(ER.name(s1),ER.name(s2)))
        end
    end
    v
end
ldh = ER.@enzyme_reaction begin
    substrates:NADH[N], Pyruvate[C]
    products:Lactate[C], NAD[N]
end
topos = ER._catalytic_topologies(ldh)
println("bi-bi topologies: ", length(topos))
println("violating topologies: ", count(t->!isempty(viol(t)), topos))
```
Run: `julia --project /tmp/repro.jl`
Expected: `bi-bi topologies: 11` and `violating topologies: 8`. This confirms the bug is reproduced and is the RED state the fix must clear.

---

## Task 1: Connectivity-invariant test helper

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (add helper in the helper block near the top, ~`:28-51`).

- [ ] **Step 1: Add the predicate helper.** After the existing topology helpers near the top of the file, add:

```julia
# Connectivity invariant: two enzyme forms identical in conformation+residual
# whose bound-metabolite sets differ by exactly one metabolite MUST be joined by
# a binding step. Returns the list of (formA, formB) pairs that violate it.
# Accepts a flat `Vector{Step}` (a topology) or a `Vector{Vector{Step}}`
# (a Mechanism's kinetic groups).
function _connectivity_violations(steps)
    flat = eltype(steps) <: AbstractVector ?
           collect(Iterators.flatten(steps)) : steps
    _bset(sp) = Set(EnzymeRates.name(m) for m in EnzymeRates.bound(sp))
    _key(sp) = (EnzymeRates.conformation(sp), EnzymeRates.residual(sp),
                Tuple(sort(collect(_bset(sp)))))
    forms = Dict{Any,Any}()
    edges = Set{Tuple{Any,Any}}()
    for s in flat
        a, b = EnzymeRates.from_species(s), EnzymeRates.to_species(s)
        forms[_key(a)] = a; forms[_key(b)] = b
        push!(edges, (_key(a), _key(b))); push!(edges, (_key(b), _key(a)))
    end
    viol = Tuple{Symbol,Symbol}[]
    fv = collect(values(forms))
    for s1 in fv, s2 in fv
        (EnzymeRates.conformation(s1) == EnzymeRates.conformation(s2) &&
         EnzymeRates.residual(s1) == EnzymeRates.residual(s2)) || continue
        b1, b2 = _bset(s1), _bset(s2)
        if length(b2) == length(b1) + 1 && issubset(b1, b2) &&
           !((_key(s1), _key(s2)) in edges)
            push!(viol, (EnzymeRates.name(s1), EnzymeRates.name(s2)))
        end
    end
    viol
end
```

- [ ] **Step 2: Sanity-check the helper compiles** by running the targeted snippet (it's used by later tasks; no standalone test run needed yet). Commit:

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "test: add connectivity-invariant predicate helper"
```

---

## Task 2: RED — assert the invariant across all three enumeration layers

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` — `@testset "_catalytic_topologies"` (`:229`), `@testset "init_mechanisms"` (`:1042`), and the `_expand_substrate_product_dead_ends` testset (locate by `grep -n '_expand_substrate_product_dead_ends' test/test_mechanism_enumeration.jl` — it has its own `@testset`).

- [ ] **Step 1: `_catalytic_topologies` testset.** In EACH sub-testset that calls `_catalytic_topologies` (Uni-Uni `:232`, Uni-Bi `:242`, Bi-Bi `:252`, Bi-Bi Ping-Pong `:262`, Ter-Ter `:272`, and the `:283` case), immediately after the `topos = EnzymeRates._catalytic_topologies(...)` line, add:

```julia
        @test all(isempty(_connectivity_violations(t)) for t in topos)
```

- [ ] **Step 2: `init_mechanisms` testset.** Find each `init_mechanisms(rxn)` call in `@testset "init_mechanisms"` and after it add an assertion over all produced mechanisms, e.g.:

```julia
        @test all(isempty(_connectivity_violations(EnzymeRates.steps(m))) for m in specs)
```
(use the variable the test binds the result to — `specs` / `init_specs` / `mechs`).

- [ ] **Step 3: dead-end testset.** In the `_expand_substrate_product_dead_ends` testset, after the call that produces the expanded `(steps, groups)` list (variable e.g. `expanded`), add:

```julia
        @test all(isempty(_connectivity_violations(steps)) for (steps, _groups) in expanded)
```

- [ ] **Step 4: Run and CONFIRM FAIL.**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep -E "Test Failed|_connectivity|catalytic_topologies|init_mechanisms|dead.end" | head -40`
Expected: the new `@test`s FAIL (bi-bi topologies 8/11 violate, LDH init violates, dead-end violates). This is the intended RED state — do NOT fix tests to pass yet.

- [ ] **Step 5: Commit the (failing) tests.**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "test: assert connectivity invariant on topologies, init, dead-ends (RED)"
```

---

## Task 3: GREEN — path-union combiner

**Files:**
- Modify: `src/mechanism_enumeration.jl` — add 3 helpers before `_catalytic_topologies` (`:164`); delete `_steps_for_ordering` (`:563-595`); replace combiner assembly (`:597-674`).

- [ ] **Step 1: Add the 3 helpers** immediately before `function _catalytic_topologies(` (line ~164):

```julia
# Substrate bound-metabolite names in route (path) order.
_binding_order(path::Vector{Step}) =
    Symbol[name(bound_metabolite(s)) for s in path
           if is_binding(s) && bound_metabolite(s) isa Substrate]

# Product bound-metabolite names in route (path) order.
_release_order(path::Vector{Step}) =
    Symbol[name(bound_metabolite(s)) for s in path
           if is_binding(s) && bound_metabolite(s) isa Product]

# True iff `order` is a linearization of weak ordering `wo` (a vector of
# levels): every metabolite of `wo` appears exactly once in `order`, and the
# level index along `order` is non-decreasing (earlier levels strictly before
# later levels; any order within a level).
function _linearizes(order::Vector{Symbol}, wo::Vector{Vector{Symbol}})
    level = Dict{Symbol,Int}()
    for (i, lvl) in enumerate(wo), m in lvl
        level[m] = i
    end
    length(order) == length(level) || return false
    prev = 0
    for m in order
        haskey(level, m) || return false
        level[m] < prev && return false
        prev = level[m]
    end
    true
end
```

- [ ] **Step 2: Delete `_steps_for_ordering`.** Remove its doc-comment + body (`:563-595`), i.e. the block from the comment starting `# Select binding steps whose source species' bound mets` through the `end` closing `_steps_for_ordering`. Leave `_weak_orderings`/`_wo_recurse!` intact.

- [ ] **Step 3: Replace the combiner assembly.** Replace this exact block (`:597-674`):

```julia
    # --- Build topologies ---
    sub_met_set = Set(sub_names)
    prod_met_set = Set(prod_names)
    is_sub(m::Metabolite) = name(m) in sub_met_set
    is_prod(m::Metabolite) = name(m) in prod_met_set

    # Iterate iso_groups in a deterministic order: prefer
    # smaller iso-step counts first (sequential before
    # ping-pong), then by sorted iso-step names. This stabilizes
    # the topology output order so test ordering invariants hold
    # — `Set{Step}` hashing is not value-stable, so unordered
    # iteration would reorder topologies.
    sorted_iso_pats = sort(collect(keys(iso_groups));
        by = pat -> (
            length(pat),
            sort([
                (string(name(from_species(s))),
                 string(name(to_species(s))))
                for s in pat])))

    result = Vector{Step}[]
    for iso_pat in sorted_iso_pats
        group_paths = iso_groups[iso_pat]
        all_group_steps = Set{Step}()
        for path in group_paths
            union!(all_group_steps, path)
        end

        # Always include all isomerization steps
        iso_steps_set = Set{Step}(
            s for s in all_group_steps if is_iso(s))

        sub_binding_mets = Set{Symbol}()
        prod_binding_mets = Set{Symbol}()
        for step in all_group_steps
            is_binding(step) || continue
            bm = bound_metabolite(step)
            if bm isa Substrate
                push!(sub_binding_mets, name(bm))
            elseif bm isa Product
                push!(prod_binding_mets, name(bm))
            end
        end

        sub_orderings = _weak_orderings(
            sort(collect(sub_binding_mets)))
        prod_orderings = _weak_orderings(
            sort(collect(prod_binding_mets)))

        seen_topos = Set{Set{Step}}()
        for sub_ord in sub_orderings, prod_ord in prod_orderings
            sub_keys = _steps_for_ordering(
                all_group_steps, sub_ord, is_sub)
            prod_keys = _steps_for_ordering(
                all_group_steps, prod_ord, is_prod)
            topo_keys = union(iso_steps_set, sub_keys,
                prod_keys)
            topo_keys ∈ seen_topos && continue
            push!(seen_topos, topo_keys)

            steps = sort(collect(topo_keys); by=s -> (
                is_iso(s) ? 1 : 0,
                string(name(from_species(s))),
                string(name(to_species(s))),
            ))

            iso_idx = findfirst(is_iso, steps)
            push!(result, Step[
                Step(from_species(s), to_species(s),
                     bound_metabolite(s), i != iso_idx)
                for (i, s) in enumerate(steps)])
        end
    end
    result
end
```

with:

```julia
    # --- Build topologies: union whole paths consistent with each
    # (substrate weak-ordering, product weak-ordering). Unioning complete
    # paths (rather than cherry-picking steps) keeps every form connected to
    # the catalytic complex — paths consistent with one weak ordering never
    # carry contradictory binding orders, so no dangling single-metabolite
    # forms arise. Binding history is read from the path, so ping-pong (where
    # a consumed substrate leaves the bound set) is handled correctly.
    #
    # Iterate iso_groups deterministically (smaller iso-step counts first,
    # then by sorted iso-step names) so topology output order is stable —
    # `Set{Step}` hashing is not value-stable.
    sorted_iso_pats = sort(collect(keys(iso_groups));
        by = pat -> (
            length(pat),
            sort([
                (string(name(from_species(s))),
                 string(name(to_species(s))))
                for s in pat])))

    result = Vector{Step}[]
    for iso_pat in sorted_iso_pats
        group_paths = iso_groups[iso_pat]

        sub_binding_mets = Set{Symbol}()
        prod_binding_mets = Set{Symbol}()
        for path in group_paths, step in path
            is_binding(step) || continue
            bm = bound_metabolite(step)
            if bm isa Substrate
                push!(sub_binding_mets, name(bm))
            elseif bm isa Product
                push!(prod_binding_mets, name(bm))
            end
        end

        sub_orderings = _weak_orderings(
            sort(collect(sub_binding_mets)))
        prod_orderings = _weak_orderings(
            sort(collect(prod_binding_mets)))

        seen_topos = Set{Set{Step}}()
        for sub_ord in sub_orderings, prod_ord in prod_orderings
            topo_keys = Set{Step}()
            matched = false
            for path in group_paths
                _linearizes(_binding_order(path), sub_ord) || continue
                _linearizes(_release_order(path), prod_ord) || continue
                union!(topo_keys, path)
                matched = true
            end
            matched || continue
            topo_keys ∈ seen_topos && continue
            push!(seen_topos, topo_keys)

            steps = sort(collect(topo_keys); by=s -> (
                is_iso(s) ? 1 : 0,
                string(name(from_species(s))),
                string(name(to_species(s))),
            ))

            iso_idx = findfirst(is_iso, steps)
            push!(result, Step[
                Step(from_species(s), to_species(s),
                     bound_metabolite(s), i != iso_idx)
                for (i, s) in enumerate(steps)])
        end
    end
    result
end
```

- [ ] **Step 4: Validate the invariant is cleared (the crux gate).** Run `/tmp/repro.jl` (Task 0 Step 2) plus a multi-reaction connectivity check. Save `/tmp/validate.jl`:

```julia
using EnzymeRates; const ER=EnzymeRates
include("/tmp/repro.jl")   # defines viol(), ldh, prints bi-bi numbers
checks = [
  ("uni-uni", ER.@enzyme_reaction(begin; substrates:S[C]; products:P[C]; end)),
  ("uni-bi",  ER.@enzyme_reaction(begin; substrates:S[CN]; products:P[C], Q[N]; end)),
  ("bi-bi",   ldh),
  ("ter-ter", ER.@enzyme_reaction(begin; substrates:A[C],B[N],C[O]; products:P[C],Q[N],R[O]; end)),
]
for (nm, rxn) in checks
    topos = ER._catalytic_topologies(rxn)
    bad = count(t->!isempty(viol(t)), topos)
    println("$nm: $(length(topos)) topologies, $bad violating")
end
```
Run: `julia --project /tmp/validate.jl`
Expected: **0 violating** for every reaction. bi-bi topology count = **11**. If any reaction still has violations, STOP — the fix is incomplete; re-analyze the combiner (do NOT layer another patch). Record the topology counts for each reaction (needed in Task 4).

- [ ] **Step 5: Commit the fix.**

```bash
git add src/mechanism_enumeration.jl
git commit -m "fix: build catalytic topologies by whole-path union (no dangling forms)"
```

---

## Task 4: Re-derive counts + specific-mechanism assertions

**Files:**
- Modify: `test/test_mechanism_enumeration.jl` (count assertions in `_catalytic_topologies`, `init_mechanisms`, dead-end testsets; add specific-mechanism assertions). Possibly other test files whose counts shift.

- [ ] **Step 1: Re-derive every enumeration count the fix touches.** Save `/tmp/counts.jl`:

```julia
using EnzymeRates; const ER=EnzymeRates
rxns = Dict(
  "uni-uni" => ER.@enzyme_reaction(begin; substrates:S[C]; products:P[C]; end),
  "uni-bi"  => ER.@enzyme_reaction(begin; substrates:S[CN]; products:P[C], Q[N]; end),
  "bi-bi"   => ER.@enzyme_reaction(begin; substrates:A[C],B[N]; products:P[C],Q[N]; end),
  "ter-ter" => ER.@enzyme_reaction(begin; substrates:A[C],B[N],C[O]; products:P[C],Q[N],R[O]; end),
)
for (nm, rxn) in sort(collect(rxns); by=first)
    nt = length(ER._catalytic_topologies(rxn))
    ni = length(ER.init_mechanisms(rxn))
    println(rpad(nm,9), " topologies=", nt, "  init_mechanisms=", ni)
end
```
Run: `julia --project /tmp/counts.jl`. Also re-derive the **bi-bi ping-pong** and the `:283` case used in the existing testset (copy those exact reaction definitions from the testset). Record all numbers.

- [ ] **Step 2: Update count assertions** to the re-derived values in `@testset "_catalytic_topologies"` (each `@test length(topos) == N`), `@testset "init_mechanisms"`, and the dead-end testset. bi-bi topologies SHOULD remain `== 11`; if any count differs from the old hardcoded value, that is expected (the old value encoded the bug) — update it.

- [ ] **Step 3: Specific-mechanism assertions (replace count-only with structure).** For **uni-uni** and **uni-bi**, capture the exact topologies from `/tmp/counts.jl`-style output, hand-verify each against the connectivity invariant + biochemistry, and assert the exact edge-sets. Helper to render a topology as a comparable Set of `(from, met, to)` tuples — add near the top of the test file:

```julia
_edge_set(topo) = Set((EnzymeRates.name(EnzymeRates.from_species(s)),
                       EnzymeRates.bound_metabolite(s) === nothing ? :iso :
                           EnzymeRates.name(EnzymeRates.bound_metabolite(s)),
                       EnzymeRates.name(EnzymeRates.to_species(s))) for s in topo)
```

For **bi-bi**, assert landmark topologies exist (encodes "ordered mechanisms are clean, not dangling-form variants"):

```julia
        topos = EnzymeRates._catalytic_topologies(bi_bi_rxn)
        formsets = [Set(EnzymeRates.name(sp) for s in t
                        for sp in (EnzymeRates.from_species(s), EnzymeRates.to_species(s)))
                    for t in topos]
        # A genuinely NADH-first-ordered substrate side has ENADH but NOT EPyruvate
        # (the buggy generator instead produced EPyruvate as a dangling form).
        @test any(fs -> (:EA in fs) && !(:EB in fs), formsets) ||
              any(fs -> (:EB in fs) && !(:EA in fs), formsets)
        # The fully-random topology has both single-substrate forms.
        @test any(fs -> (:EA in fs) && (:EB in fs), formsets)
```
(Use the actual single-substrate form names for `bi_bi_rxn` — derive them from `/tmp/counts.jl` output; e.g. for substrates `A,B` they are `:EA`,`:EB`.)

- [ ] **Step 4: Run the enumeration tests to confirm green.**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep -E "test_mechanism_enumeration|Test Summary|Fail|Error" | head -30`
Expected: enumeration testsets pass (counts updated, invariant holds, specific mechanisms present).

- [ ] **Step 5: Commit.**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "test: re-derive enumeration counts and assert specific mechanisms"
```

---

## Task 5: Full suite, downstream counts, gates, docs

**Files:**
- Modify: any test asserting enumeration-derived counts that shifted (search broadly); `test/test_compile_budget.jl` golden if it moved; `.claude/CLAUDE.md` "Verified topology counts" line.

- [ ] **Step 1: Full suite.**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -40`
Triage every new failure. Failures will cluster in tests that hardcoded enumeration-derived counts (`init_mechanisms` sizes, dead-end counts, compile-budget). For each: confirm the new number is correct (it reflects the now-clean enumeration), update the assertion. Do NOT update a count without confirming via a direct `_catalytic_topologies` / `init_mechanisms` / `fitted_params` probe that the new value is right.

- [ ] **Step 2: Compile-budget gate.** If `test/test_compile_budget.jl` fails (the bi-bi init trace-compile golden, ~750), re-derive the new trace-compile count and update the golden with a one-line comment noting it changed because the clean enumeration produces a different (smaller/cleaner) mechanism set. Confirm the perf test (`test_rate_equation_performance`) still passes (derivation untouched — it should).

- [ ] **Step 3: Update CLAUDE.md** if any topology count changed: `.claude/CLAUDE.md` "Verified topology counts: bi-bi=11, ter-ter=283, …" — set to the re-derived values (bi-bi stays 11). Update the "init_mechanisms … competition-filtered" prose only if a documented count is now wrong.

- [ ] **Step 4: Final full-suite green.**

Run: `julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -5`
Expected: pass count ≥ Task-0 baseline (minus only intentionally-changed assertions), 0 fail, the same 1 pre-existing `@test_broken`.

- [ ] **Step 5: Commit + push.**

```bash
git add -A   # only after a git status confirms only intended files
git status
git commit -m "fix: update downstream enumeration counts + docs for connectivity fix"
git push
```

---

## Self-review checklist (run before declaring done)

- [ ] Connectivity invariant asserted in ALL THREE layers (`_catalytic_topologies`, `init_mechanisms`, `_expand_substrate_product_dead_ends`) and green.
- [ ] `_connectivity_violations` returns empty for uni-uni, uni-bi, bi-bi, bi-bi ping-pong, ter-ter.
- [ ] bi-bi topology count = 11; all other enumeration counts re-derived & asserted (no count left at a value that was never re-verified post-fix).
- [ ] Specific-mechanism assertions present for the small reactions; bi-bi landmark (ordered vs random) assertions present.
- [ ] `_steps_for_ordering` deleted; `_binding_order`/`_release_order`/`_linearizes` added; no dangling references.
- [ ] Full suite green (0 fail); perf gate green; compile-budget gate green (re-baselined w/ justification if it moved).
- [ ] Exports still 18; `rate_equation` derivation untouched.
