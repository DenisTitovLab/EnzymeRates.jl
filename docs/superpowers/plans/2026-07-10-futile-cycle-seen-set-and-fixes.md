# Futile-cycle: seen-set + canonicalization + RE→SS fixes — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the beam search terminate on the LDH four-inhibitor problem and fix the two confirmed correctness bugs behind the futile cycle.

**Architecture:** A structural seen-set in `_process_batch` (the primary termination guarantee) plus two enumeration/canonicalization correctness fixes: idempotent `_canonical_mechanism`, and `_expand_re_to_ss` RE-only invariants. Each ships as its own commit.

**Tech Stack:** Julia 1.12, EnzymeRates.jl. Tests via `Pkg.test()` or focused `TestEnv` runs.

## Global Constraints

- 92-char line length, 4-space indent. Match surrounding style.
- `rate_equation` MUST stay allocation-free and sub-120 ns (`test_rate_eq_derivation.jl` perf gate). None of these changes touch the derived `rate_equation` body, so the gate must stay green unchanged.
- Parameter-name rendering flows through `name(p, m)`; the AST-walker guard (`test_types.jl`) must stay green.
- Machine is memory-limited (8 GB, 4 cores): run Julia **serially**, never two suites at once. Use focused single-file test runs while iterating; full `Pkg.test()` only at task end.
- TDD: failing test → confirm fail → minimal implement → confirm pass → commit.
- A "golden re-baseline" is a **validated** regeneration: regenerate, diff old vs new, confirm only the intended mechanisms changed and each new form passes an independent correctness check (equilibrium-flux oracle `v=0` at `Q=Keq`; parameter-count sanity). If an unexpected mechanism changes or a change can't be validated, STOP.

---

### Task 1: Structural seen-set (termination)

**Files:**
- Modify: `src/identify_rate_equation.jl` — `_process_batch` (524), `_batch_summary` (472), `_beam_search` (655) base call (677) + child call (750) + all-skipped branch (782).
- Test: `test/test_identify_rate_equation.jl`

**Interfaces:**
- Produces: `_process_batch(...; seen::Set{UInt64}=Set{UInt64}(), ...)` now returns a 5-tuple `(entries, failures, n_param_skip, n_complexity_skip, n_seen_skip)`. `_batch_summary(...; n_seen_skipped::Int, ...)`.

- [ ] **Step 1: Write the failing test** — a batch containing a duplicate structure reports it as already-seen and does not re-emit it.

```julia
@testset "seen-set: repeat structures are skipped, not reprocessed" begin
    rxn = @enzyme_reaction begin
        substrates: A[C]
        products: P[C]
    end
    ms = EnzymeRates.init_mechanisms(rxn)
    m = first(ms)
    prob = IdentifyRateEquationProblem(rxn, _uniuni_data(); Keq=1.0)  # helper below
    seen = Set{UInt64}()
    e1, f1, ps1, cs1, ss1 = EnzymeRates._process_batch([m], prob;
        optimizer=_test_optimizer(), max_param_count=13, memo=Dict{UInt64,NamedTuple}(), seen)
    @test ss1 == 0 && length(e1) == 1
    # same structure again in a later batch → seen-skipped, no new entry
    e2, f2, ps2, cs2, ss2 = EnzymeRates._process_batch([m], prob;
        optimizer=_test_optimizer(), max_param_count=13, memo=Dict{UInt64,NamedTuple}(), seen)
    @test ss2 == 1 && isempty(e2) && isempty(f2)
end
```

(Reuse whatever `IdentifyRateEquationProblem` / data / optimizer helpers the neighboring `_process_batch` tests already use in this file; do not invent new ones if a fixture exists.)

- [ ] **Step 2: Run test to verify it fails**

Run (focused): `julia --project=@. -e 'using TestEnv; TestEnv.activate(); include("test/test_identify_rate_equation.jl")'` (or the repo's focused-run recipe).
Expected: FAIL — `_process_batch` currently returns a 4-tuple and takes no `seen`.

- [ ] **Step 3: Add `seen` filtering to `_process_batch`**

Add the keyword and a master-side pre-filter at the very top of the function body (before PASS-1), and add `n_seen_skip` to the return:

```julia
function _process_batch(
    mechs, prob::IdentifyRateEquationProblem;
    optimizer, max_param_count, eq_complexity_filter::Int = typemax(Int),
    memo::Dict{UInt64,NamedTuple}=Dict{UInt64,NamedTuple}(),
    seen::Set{UInt64}=Set{UInt64}(),
    parent_of::AbstractDict = Dict(), kwargs...
)
    # Skip structures already produced in an earlier batch — expand each once.
    # Added to `seen` on first sight regardless of outcome (fit / cap / error).
    fresh = empty(mechs)
    n_seen_skip = 0
    for m in mechs
        h = hash(m)
        if h in seen
            n_seen_skip += 1
        else
            push!(seen, h)
            push!(fresh, m)
        end
    end

    # PASS 1 (workers): complexity-cap + compile + param-cap + render.
    compiled = pmap(fresh) do m
        # … unchanged body …
    end
    # … unchanged PASS-2 and row assembly …
    return entries, failures,
           count(x -> x === nothing, compiled),          # param-count skips
           count(x -> x === :complexity_skip, compiled), # complexity skips
           n_seen_skip
end
```

The only edits: the new `seen` kwarg, the pre-filter block, `pmap(fresh)` instead of `pmap(mechs)`, and the added return element. Everything between is untouched.

- [ ] **Step 4: Add the bucket to `_batch_summary`**

```julia
function _batch_summary(
    entries::Vector{BatchEntry}, failures::Vector{FitFailure};
    n_param_skipped::Int, n_complexity_skipped::Int, n_seen_skipped::Int,
    max_param_count::Int, eq_complexity_filter::Int)
    n_fit   = length(entries)
    n_err   = length(failures)
    n_new   = count(e -> !e.row.fit_inherited, entries)
    n_inh   = n_fit - n_new
    n_succ  = count(e -> e.retcode === :Success, entries)
    n_other = n_fit - n_succ
    pct(x)  = n_fit == 0 ? 0.0 : round(100 * x / n_fit; digits=1)
    string(n_new, " new fits + ", n_inh, " inherited + ",
           n_seen_skipped, " skipped (already seen) + ",
           n_param_skipped, " skipped (>", max_param_count, " params) + ",
           n_complexity_skipped, " skipped (>", eq_complexity_filter,
           " complexity) + ", n_err, " errored | Success ", pct(n_succ),
           "% | non-Success retcode ", pct(n_other), "%")
end
```

Update its docstring "five buckets" → "six buckets".

- [ ] **Step 5: Thread `seen` through `_beam_search`**

After `memo = Dict{UInt64,NamedTuple}()` (≈ line 670) add:

```julia
    # Structures already produced — each is expanded at most once (termination).
    seen = Set{UInt64}()
```

Base call (677): unpack 5 values, pass `seen`, pass the new summary kwarg:

```julia
    base_entries, base_failures, n_base_param_skip, n_base_cx_skip, n_base_seen_skip =
        _process_batch(base, prob;
            optimizer, max_param_count, eq_complexity_filter, memo, seen, kwargs...)
```
```julia
        _batch_summary(base_entries, base_failures;
                       n_param_skipped=n_base_param_skip,
                       n_complexity_skipped=n_base_cx_skip,
                       n_seen_skipped=n_base_seen_skip,
                       max_param_count, eq_complexity_filter),
```

Child call (750): same shape:

```julia
            child_entries, child_failures, n_child_param_skip, n_child_cx_skip, n_child_seen_skip =
                _process_batch(children, prob;
                    optimizer, max_param_count, eq_complexity_filter, memo,
                    seen, parent_of, kwargs...)
```
```julia
                    _batch_summary(child_entries, child_failures;
                                   n_param_skipped=n_child_param_skip,
                                   n_complexity_skipped=n_child_cx_skip,
                                   n_seen_skipped=n_child_seen_skip,
                                   max_param_count, eq_complexity_filter),
```

All-skipped branch (782-790): add seen to the report:

```julia
                _progress(save_dir, show_progress, string(
                    "Expanded ", length(to_expand), " parents → ",
                    length(children), " children | all skipped (",
                    n_child_seen_skip, " already seen, ",
                    n_child_param_skip, " >", max_param_count, " params, ",
                    n_child_cx_skip, " >", eq_complexity_filter, " complexity)"))
```

- [ ] **Step 6: Run the new test — verify it passes.** Focused run as in Step 2. Expected: PASS.

- [ ] **Step 7: Run the identify-rate-equation suite** to confirm no regression in beam behavior or the selection tests.

Run: focused `test/test_identify_rate_equation.jl`. Expected: all PASS (selected model unchanged; existing `_process_batch` call sites in tests that unpack 4 values must be updated to 5 — grep `_process_batch(` in the test file and fix any).

- [ ] **Step 8: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Beam: structural seen-set stops re-expansion (terminates the futile cycle)"
```

---

### Task 2: Idempotent `_canonical_mechanism` (class A)

**Files:**
- Modify: `src/mechanism_enumeration.jl` — `_canonical_mechanism(::Mechanism)` (1401), `_canonical_mechanism(::AllostericMechanism)` (1403).
- Test: `test/test_identify_rate_equation.jl` (near the existing canon tests, ~1576).

**Interfaces:**
- Produces: `_canonical_mechanism` reaches a fixed point (`_canonical_mechanism(_canonical_mechanism(m)) == _canonical_mechanism(m)`).

- [ ] **Step 1: Write the failing test** — idempotency over all specs plus the reproducer parent.

```julia
@testset "_canonical_mechanism is idempotent" begin
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism isa EnzymeRates.AllostericEnzymeMechanism ?
            EnzymeRates.AllostericMechanism(m_from_spec(spec)) : # use the repo's spec→concrete helper
            EnzymeRates.Mechanism(m_from_spec(spec))
        c = EnzymeRates._canonical_mechanism(m)
        @test EnzymeRates._canonical_mechanism(c) == c
    end
    # Reproducer A parent: a split no-op canonicalizes 9→9 on one pass, 9→8 on two.
    amA = EnzymeRates.AllostericMechanism(Core.eval(EnzymeRates, Meta.parse(A_SIG))())
    c = EnzymeRates._canonical_mechanism(amA)
    @test EnzymeRates._canonical_mechanism(c) == c
end
```

(Use whatever spec→concrete-`Mechanism`/`AllostericMechanism` conversion the file already uses; `A_SIG` is the string in `docs/superpowers/specs/2026-07-10-futile-cycle-reproducers.jl`.)

- [ ] **Step 2: Run to verify it fails** (focused). Expected: FAIL on at least one mechanism whose second pass differs.

- [ ] **Step 3: Iterate the merge to a fixed point**

```julia
function _canonical_mechanism(m::Mechanism)
    prev = m
    for _ in 1:8   # convergence is ≤2 passes; the bound guards a pathological loop
        merged = Mechanism(reaction(prev), _merge_tied_kinetic_groups(prev))
        merged == prev && return merged
        prev = merged
    end
    prev
end

function _canonical_mechanism(am::AllostericMechanism)
    prev = am
    for _ in 1:8
        cat_steps, cat_states = _merge_tied_kinetic_groups(prev)
        merged = AllostericMechanism(reaction(prev), cat_steps, cat_states,
                                     catalytic_multiplicity(prev),
                                     copy(regulatory_sites(prev)))
        merged == prev && return merged
        prev = merged
    end
    prev
end
```

The existing `@test_throws ErrorException _canonical_mechanism(m_bad)` still holds — the first iteration's constructor throws exactly as before.

- [ ] **Step 4: Run the new test — verify PASS** (focused).

- [ ] **Step 5: Golden re-baseline check.** Run the full suite; if `test_allosteric_golden.jl` fails, some spec's canonical form changed. For each changed spec, confirm the change is a legitimate additional merge (a Wegscheider-tied pair now fused) and that the equilibrium-flux oracle still holds for it, then regenerate `test/reference/allosteric_golden_reference.txt` from `_allosteric_golden_lines()`. If a change is not a merge or the oracle breaks, STOP and record it.

Run: `julia --project -e 'using Pkg; Pkg.test()'`

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_identify_rate_equation.jl test/reference/allosteric_golden_reference.txt
git commit -m "Enumeration: iterate _canonical_mechanism to a fixed point (drops delta-0 split no-ops)"
```

---

### Task 3: `_expand_re_to_ss` — inhibitor bindings RE-only (invariant 1)

**Files:**
- Modify: `src/mechanism_enumeration.jl` — `_expand_re_to_ss` (1207).
- Test: `test/test_mechanism_enumeration.jl` — new assertion + update the "Substrate-as-dead-end-inhibitor overlap" testset (1719).

**Interfaces:**
- Consumes: `bound_metabolite(::Step)`, `Regulator` abstract type (`CompetitiveInhibitor <: Regulator`).
- Produces: `_expand_re_to_ss` never flips a group whose binding step binds a `Regulator`.

- [ ] **Step 1: Write the failing test** — no SS inhibitor binding is ever produced.

```julia
@testset "_expand_re_to_ss keeps inhibitor bindings RE" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: S
    end
    m = first(EnzymeRates._expand_add_dead_end_regulator(
              first(EnzymeRates.init_mechanisms(rxn)), rxn))
    for r in EnzymeRates._expand_re_to_ss(m)
        for grp in EnzymeRates.steps(r), s in grp
            if EnzymeRates.bound_metabolite(s) isa EnzymeRates.Regulator
                @test EnzymeRates.is_equilibrium(s)   # inhibitor binding stays RE
            end
        end
    end
end
```

- [ ] **Step 2: Run to verify it fails** (focused `test_mechanism_enumeration.jl`). Expected: FAIL — the dead-end-inhibitor group is currently flipped to SS.

- [ ] **Step 3: Add invariant 1 to `_expand_re_to_ss`**

```julia
function _expand_re_to_ss(m::Union{Mechanism, AllostericMechanism})
    results = typeof(m)[]
    for g in kinetic_groups(m)
        grp = steps(m)[g]
        all(is_equilibrium, grp) || continue
        # Inhibitor (dead-end) bindings are rapid-equilibrium only — their speed
        # is never identifiable; never flip them to steady state.
        any(s -> bound_metabolite(s) isa Regulator, grp) && continue
        push!(results, _with_steps(m, _flip_group_to_ss(steps(m), g)))
    end
    results
end
```

- [ ] **Step 4: Update the "Substrate-as-dead-end-inhibitor overlap" testset** (1719). The inhibitor group no longer flips, so there are **2** variants (substrate + product), not 3:

```julia
        result = EnzymeRates._expand_re_to_ss(m)
        # substrate-binding (RE) and product-binding (RE) flip; the
        # dead-end-inhibitor group stays RE (inhibitor bindings are RE-only).
        @test length(result) == 2
```
and change the `@test length(unique(flipped_groups)) == 3` to `== 2`, and drop the "treats … as independent … → 3" comment wording to reflect the RE-only rule.

- [ ] **Step 5: Run to verify PASS** (focused). Then scan the rest of the `_expand_re_to_ss` testsets (1364–1820) for any other assertion that counted an inhibitor-binding flip; update counts to the RE-only rule and confirm each new count by inspection.

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Enumeration: _expand_re_to_ss keeps competitive-inhibitor bindings RE-only"
```

---

### Task 4: `_expand_re_to_ss` — mirror type-lock (invariant 2) — GATED

This is the subtle one. Implement only if it validates cleanly against the B2 reproducer and the full suite; otherwise leave Tasks 1–3 landed and record the blocker for Denis.

**Files:**
- Modify: `src/mechanism_enumeration.jl` — `_expand_re_to_ss` (1207), add a `_step_core` / mirror-class helper.
- Test: `test/test_mechanism_enumeration.jl`.

**Interfaces:**
- Consumes: `bound(::Species)`, `conformation`, `residual`, `Species` ctor, `Regulator`.
- Produces: a catalytic binding group and its inhibitor-bound mirror group(s) always share RE/SS type after `_expand_re_to_ss`.

- [ ] **Step 1: Write the failing test** — starting from a mechanism where a catalytic binding and its inhibitor-bound mirror are BOTH RE in separate groups, `_expand_re_to_ss` never flips one without the other.

```julia
@testset "_expand_re_to_ss flips inhibitor-bound mirrors together" begin
    # Build a bi-substrate mechanism + dead-end inhibitor, split so the
    # inhibitor-bound Pyruvate mirror is its own RE group, then assert no
    # variant leaves the base RE while the mirror is SS (or vice versa).
    # (Construct via init_mechanisms + _expand_add_dead_end_regulator +
    #  _expand_split_kinetic_group on the LDH-like reaction; identify the
    #  base group and its mirror by shared inhibitor-free core.)
    # For every result r and every pair (base g, mirror g') sharing a core:
    #   is_all_ss(r,g) == is_all_ss(r,g')
end
```

Flesh this out concretely during implementation using the same reaction shape as the B2 reproducer; the oracle is: define `core(step)` = (from-species with regulators stripped, to-species with regulators stripped, bound_metabolite); two groups are mirror-linked iff they share a core; assert linked groups have equal all-SS status in every result.

- [ ] **Step 2: Run to verify it fails.**

- [ ] **Step 3: Implement mirror-class flipping.**

```julia
# Inhibitor-free core of a step: the same catalytic binding with all Regulators
# stripped from both species. Two steps with equal cores are the same binding in
# different inhibitor contexts (mirrors).
function _step_core(s::Step)
    strip(sp) = Species(Metabolite[b for b in bound(sp) if !(b isa Regulator)],
                        conformation(sp), residual(sp))
    (strip(from_species(s)), strip(to_species(s)), bound_metabolite(s))
end

# Connected components of the graph groups—share—a—core, restricted to eligible
# groups (all-RE, no regulator binding). Each component flips atomically.
function _re_to_ss_flip_units(m)
    elig = [g for g in kinetic_groups(m)
            if all(is_equilibrium, steps(m)[g]) &&
               !any(s -> bound_metabolite(s) isa Regulator, steps(m)[g])]
    core_of = Dict{Any, Vector{Int}}()
    for g in elig, s in steps(m)[g]
        push!(get!(core_of, _step_core(s), Int[]), g)
    end
    parent = Dict(g => g for g in elig)
    find(x) = (parent[x] == x ? x : (parent[x] = find(parent[x])))
    union!(a, b) = (parent[find(a)] = find(b))
    for gs in values(core_of), i in 2:length(gs)
        union!(gs[1], gs[i])
    end
    comps = Dict{Int, Vector{Int}}()
    for g in elig
        push!(get!(comps, find(g), Int[]), g)
    end
    collect(values(comps))
end

function _expand_re_to_ss(m::Union{Mechanism, AllostericMechanism})
    results = typeof(m)[]
    for unit in _re_to_ss_flip_units(m)
        new_groups = steps(m)
        for g in unit
            new_groups = _flip_group_to_ss(new_groups, g)
        end
        push!(results, _with_steps(m, new_groups))
    end
    results
end
```

Note `_flip_group_to_ss` returns a fresh groups vector; chaining it per group in the unit flips them all. `union!` here is a local helper (shadow is fine; it operates on the `parent` dict).

- [ ] **Step 4: Run the new test — verify PASS.**

- [ ] **Step 5: Validate against B2** — with the both-RE precursor (dead-end SS bindings RE), the fully-flipped form is fully identifiable (`np 10→8`, rank 8) and the split produces no distinct-text duplicate (reuse `docs/superpowers/specs/2026-07-10-futile-cycle-reproducers.jl` and the identifiability check from the seen-set investigation scratchpad). Equilibrium-flux oracle holds.

- [ ] **Step 6: Full suite + enumeration re-baseline.** Run `Pkg.test()`. Update `_expand_re_to_ss` count assertions (1364–1820) to the mirror-linked behavior; confirm each new count by inspection (mirror groups now share a variant). If any count change can't be explained by "mirror groups flip together," STOP and record it.

- [ ] **Step 7: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Enumeration: _expand_re_to_ss flips inhibitor-bound mirrors with their base"
```

---

## Verification (whole change)

- [ ] Full `Pkg.test()` green (incl. golden re-baselines, the `rate_equation` 0-alloc / sub-120 ns gate, the parameter-naming chokepoint).
- [ ] The A and B2 reproducers: A no-op dropped; B2 parent fully identifiable / no split duplicate (after Task 4).
- [ ] Bounded local LDH four-inhibitor beam run terminates by draining the frontier (narrow width, low max-time) — the seen-set's end-to-end proof. Deferred if too heavy locally; note it for an at-scale run.

## Self-review notes

- Spec coverage: seen-set (Task 1), idempotent canon (Task 2), `_expand_re_to_ss` invariants (Tasks 3–4). Defect 2 and the change_allo prune are spec non-goals — no task, by design.
- Task 1 is behavior-preserving for the selected model (no golden re-baseline). Tasks 2–4 each carry their own validated re-baseline.
- Task 4 is gated: land Tasks 1–3 regardless; attempt Task 4 and, if it does not validate cleanly, leave it uncommitted with a written blocker.
