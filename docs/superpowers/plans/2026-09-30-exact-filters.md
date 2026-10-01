# Exact Enumeration Filters Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the enumerator stop emitting mechanisms with a steady-state kinetic group that carries no flux, or a dead-end copy of a substrate or product that creates no new complex, without losing any family the moves reach today.

**Architecture:** One structural flux predicate on the rapid-equilibrium segment graph replaces the `is_iso` block screen and is enforced inside the flip (a zero-flux flip set counts as failed and is extended) and the split (a part with no flux-carrying step is rebuilt at rapid equilibrium). One twin predicate with two keys (composition; RE segment plus offsets) is enforced in the dead-end move (a pattern whose sites are all twins is skipped) and the split (a bipartition isolating twin-only copy sites is not a unit). `expand_mechanisms` asserts both rules on each parent; `seed_mechanisms` errors when no seed binds every required regulator.

**Tech Stack:** Julia 1.12, the EnzymeRates package (`src/`), `Test`, `TestEnv` for focused runs.

**Spec:** `docs/superpowers/specs/2026-09-30-exact-filters-design.md`. Evidence: `docs/superpowers/specs/2026-09-28-identifiable-enumeration-findings-and-plan.md` points 1 and 2, and the track reports `t1_flux_report.md`, `t2_inhdup_report.md` with their skeptic verdicts in `docs/superpowers/specs/2026-09-28-identifiability-evidence/`.

## Global Constraints

- Both rules are emission policies of the moves. The `Mechanism` and `AllostericMechanism` constructors, `_assert_mechanism_invariants`, `init_mechanisms` and every derivation function stay unchanged. A hand-written mechanism with a zero-flux SS group or a duplicate copy still derives.
- The rules are group-level, structural and exact: no floating point, no rank, no random parameters anywhere in `src/`.
- The flux test runs on the catalytic step graph as stored in `steps(m)`, which for an `AllostericMechanism` is its A-state graph. Nothing consults the I-state.
- The flip uses extend semantics (a zero-flux flip set counts as failed inside `gains`). The split uses revert semantics (a part with no flux-carrying step becomes RE) and evaluates its gain test on the reverted child.
- Tests in `test/test_mechanism_enumeration.jl` follow the three "Enumeration-engine tests" rules in CLAUDE.md: every fixture inline with `@enzyme_mechanism` or `@allosteric_mechanism`, exact child sets and counts (property assertions only in addition), chemistry written the enumerator's way (isomerization to a product-bound form, then a release). File-level helpers are prefixed `_testhelper_`. An aggregate pin over a whole enumerated population says so in a comment.
- Existing tests whose expectations change are rewritten to keep their intent, never deleted (Denis, 2026-09-30). Every rewritten exact-children expectation is derived by hand from the rules in a comment and cross-checked with `_testhelper_identifiable_rank`.
- 92-character lines, 4-space indentation; match surrounding style; never change unrelated whitespace; docstrings on every new function; no comment or name that refers to what the code used to do.
- All `Parameter → Symbol` rendering goes through `name(p, m)`; no `Symbol("K…")`/`Symbol("k…")` literal outside a renderer (guarded by `test/test_types.jl`). This plan adds no parameter rendering.
- `rate_equation` stays allocation-free and under 120 ns per call for every spec (`test_rate_equation_performance`). This plan touches no derivation code, so the guard must pass unchanged.
- TDD: write the failing test, run it and see it fail, implement, see it pass.
- Only one Julia process at a time (7.7 GB RAM, no swap). The full suite takes about 14 minutes: run it in the background with output to a log and poll the log in a bounded foreground loop (`for i in $(seq 1 40); do sleep 30; grep -q "Test Summary" LOG && break; done`); never wait on a monitor. A background job started by a subagent dies when the subagent hands back, so the agent that starts the suite polls it to the end. Before reporting a process as running, check it is alive with `pgrep -f bin/julia`.
- Full suite: `julia --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`.
- Focused run of the enumeration file (about 4 minutes, warm): `julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_mechanism_enumeration.jl")'` (skips Aqua/JET; an `UndefVarError` for a helper defined in another test file is an artifact of the focused run). For the red/green cycle of one new testset, copy that testset into a scratch file under the session scratchpad, prefix it with the same `using`/`include` lines plus `const ER = EnzymeRates`, and run that; run the whole enumeration file before committing.
- Line references below are to the branch at commit 3f79e32 and drift as tasks land; find functions by name.
- Commit after each task with a message ending in
  `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>` and
  `Claude-Session: https://claude.ai/code/session_01R4zCpbZSoygD66kRecfrDN`. Never skip hooks; never `git add -A`.

## Review Focus

1. A parent whose SS group is a self-loop of the segment graph (both forms in one RE segment): an unbalanced self-loop (chemistry inside a segment) carries flux and a balanced one (a binding whose two forms an RE route already joins) does not; Task 1 Step 1 pins both.
2. A reactant copy bound at a form whose composition is new but whose weight is proportional to an existing form's (the conformational-isomer and ping-pong cases): the dead-end move must skip it; Task 3 Step 1 pins the predicate on both fixtures and the move on the isomer fixture.
3. A split of a copy group whose one part has only twin sites must not be a unit even when the other part is fine; Task 4 Step 1 pins the exact (empty) child set.
4. A hand-written parent that violates a rule reaches `expand_mechanisms`: the error names the group's steps; Task 5 Step 1 pins both messages.
5. A required competitive inhibitor with no admissible placement (uni-uni with S or P as its own inhibitor) must error from `seed_mechanisms`, not return an empty result; Task 5 Step 1 pins the error and its message.

---

### Task 1: The flux predicate on the segment graph

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_flux_carrying_groups` and its docstring, about lines 1343–1379; the `units` line of `_expand_re_to_ss`, about line 1261; the `_expand_re_to_ss` docstring's sentence on `_flux_carrying_groups`)
- Test: `test/test_mechanism_enumeration.jl` (the `@testset "_flux_carrying_groups"` block, about lines 5434–5543)

**Interfaces:**
- Produces (used by Tasks 2, 4, 5, 6):
  - `_all_steady_state(groups::Vector{Vector{Step}})::Vector{Vector{Step}}` — every step rebuilt with `is_equilibrium = false`.
  - `_flux_carrying_steps(groups::Vector{Vector{Step}}, rxn::EnzymeReaction)::Vector{BitVector}` — one flag per step, parallel to `groups`; RE steps are `false`.
  - `_flux_carrying_groups(groups::Vector{Vector{Step}}, rxn::EnzymeReaction)::BitVector` — `any` over each group's step flags.
  - `_flux_carrying_groups(m::Union{Mechanism, AllostericMechanism})::BitVector` — the same on `steps(m)`, `reaction(m)`.

- [ ] **Step 1: Write the failing tests**

In `test/test_mechanism_enumeration.jl`, inside `@testset "_flux_carrying_groups"`, change every existing call `EnzymeRates._flux_carrying_groups(x)` (six of them: `m`, `leaf`, `pingpong`, `mirror`, `pendant`, `pp`) to the all-steady-state form, for example:

```julia
    @test all(EnzymeRates._flux_carrying_groups(
        EnzymeRates._all_steady_state(EnzymeRates.steps(m)), EnzymeRates.reaction(m)))
```

and likewise `fc = EnzymeRates._flux_carrying_groups(EnzymeRates._all_steady_state(EnzymeRates.steps(leaf)), EnzymeRates.reaction(leaf))`, and the same for `pingpong`, `mirror`, `pendant` (`pfc = ...`) and `pp`. Their assertions stay as they are: the all-steady-state flags equal today's screen on every mechanism the enumerator emits.

Then append these cases before the testset's closing `end`:

```julia
    # A Theorell–Chance step consumes B and releases P in one edge. The cycle
    # E → E(A) → E(Q) → E has weight 1 + 2 + 1 = 4 = one turnover times the
    # four reactants, so it is unbalanced and every step carries flux. The
    # former screen looked for an isomerization step and flagged nothing here.
    tc = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(Q) + P
            E + Q ⇌ E(Q)
        end
    end)
    @test all(EnzymeRates._flux_carrying_groups(
        EnzymeRates._all_steady_state(EnzymeRates.steps(tc)), EnzymeRates.reaction(tc)))
    @test EnzymeRates._flux_carrying_groups(tc) ==
        [EnzymeRates.is_equilibrium(first(g)) ? false : true for g in EnzymeRates.steps(tc)]

    # A merged central complex X with two fused steps into it and no
    # isomerization: E + A ⇌ E(A), E(A) + B → X, E(Q) + P → X, E + Q ⇌ E(Q).
    # Every step lies on the unbalanced cycle E → E(A) → X → E(Q) → E.
    merged = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(P, Q)
            E(Q) + P <--> E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    @test all(EnzymeRates._flux_carrying_groups(
        EnzymeRates._all_steady_state(EnzymeRates.steps(merged)),
        EnzymeRates.reaction(merged)))

    # An RE shunt. E + A and E + Q are steady-state while the abortive route
    # E(Q) + A ⇌ E(A, Q) ⇌ E(A) + Q stays at rapid equilibrium. The segment
    # {E} joins the rest by two parallel edges of weight 0 (A uptake +1 against
    # E(A)'s offset +1; Q uptake −1 against E(Q)'s offset −1): a balanced block,
    # so neither steady-state binding carries flux. The turnover runs inside the
    # big segment, where the chemistry step is an unbalanced self-loop of
    # weight 4 (E(A, B) sits at A + B, E(P, Q) at −P − Q). On the all-steady-state
    # graph every group carries flux, which is what the former screen reported.
    shunt = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(Q) + A ⇌ E(A, Q)
            E + Q <--> E(Q)
            E(A) + Q ⇌ E(A, Q)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
        end
    end)
    binder(grp) = EnzymeRates.bound_metabolite(first(grp))
    source(grp) = EnzymeRates.name(EnzymeRates.from_species(first(grp)))
    ss_binding(grp) = binder(grp) !== nothing && !EnzymeRates.is_equilibrium(first(grp))
    flags = EnzymeRates._flux_carrying_groups(shunt)
    for (g, grp) in enumerate(EnzymeRates.steps(shunt))
        if ss_binding(grp)
            @test source(grp) == :E && !flags[g]
        elseif EnzymeRates.is_iso(first(grp))
            @test flags[g]
        else
            @test EnzymeRates.is_equilibrium(first(grp)) && !flags[g]
        end
    end
    @test count(ss_binding, EnzymeRates.steps(shunt)) == 2
    @test all(EnzymeRates._flux_carrying_groups(
        EnzymeRates._all_steady_state(EnzymeRates.steps(shunt)),
        EnzymeRates.reaction(shunt)))

    # A balanced self-loop: E + A is steady-state, but E and E(A) already share
    # a segment through E ⇌ E(Q) ⇌ E(A, Q) ⇌ E(A), so the step is a self-loop
    # of weight 1 + 0 − 1 = 0 and carries no flux. The chemistry self-loop in
    # the same segment has weight 4 and does.
    loop = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(Q) + A ⇌ E(A, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
        end
    end)
    lflags = EnzymeRates._flux_carrying_groups(loop)
    for (g, grp) in enumerate(EnzymeRates.steps(loop))
        if ss_binding(grp)
            @test !lflags[g]
        elseif EnzymeRates.is_iso(first(grp))
            @test lflags[g]
        end
    end
    @test count(ss_binding, EnzymeRates.steps(loop)) == 1
    @test EnzymeRates._re_segment_count(loop) == 1
```

- [ ] **Step 2: Run the testset and watch it fail**

Copy the whole `@testset "_flux_carrying_groups"` into a scratch file with the focused-run prelude and run it. Expected: `UndefVarError: _all_steady_state not defined` (and `MethodError` on the two-argument `_flux_carrying_groups`).

- [ ] **Step 3: Implement the predicate**

In `src/mechanism_enumeration.jl`, replace the `_flux_carrying_groups` docstring and function (keep `_edge_blocks` above it unchanged) with:

```julia
"""Every step of `groups` rebuilt at steady state. The flip move pre-filters its
units on this graph: a step that carries no flux with every step steady-state
carries none under any assignment (flips only refine segments), so a group with
no such step can never flip usefully."""
_all_steady_state(groups::Vector{Vector{Step}}) =
    [Step[Step(from_species(s), to_species(s), consumed(s), released(s), false)
          for s in group] for group in groups]

"""
    _flux_carrying_steps(groups, rxn) -> Vector{BitVector}

One flag per step of `groups`, parallel to `groups`: whether the step carries net
steady-state flux for generic parameter values. Rapid-equilibrium steps hold
equilibrium mass and are flagged `false`; the rule that uses these flags concerns
steady-state groups only.

Contract each rapid-equilibrium segment to one vertex (`_re_segment_extras`) and
make every steady-state step an edge between its two forms' segments, weighted by
the substrates it takes up minus the products, plus its `from` form's offsets minus
its `to` form's. Around any cycle the offsets telescope and the uptakes sum to the
net turnover times the reactant count, so a cycle has weight zero exactly when it
runs no net reaction. A self-loop (both forms in one segment) carries flux iff its
weight is nonzero. Any other edge carries flux iff its biconnected block
(`_edge_blocks`) holds a cycle of nonzero weight: every edge of such a block lies
on one (join the edge to the cycle by two disjoint paths; one of the two resulting
cycles is unbalanced), while in a balanced block detailed balance holds along every
cycle and each step's flux vanishes. The zero-flux verdict holds for any parameters
and grouping; the flux-carrying verdict for one-way steps, which every step the
enumerator emits is. Reading the metabolite lists, the test covers fused and
Theorell–Chance steps.
"""
function _flux_carrying_steps(groups::Vector{Vector{Step}}, rxn::EnzymeReaction)
    species, segments, extras = _re_segment_extras(groups)
    idx = Dict(sp => i for (i, sp) in enumerate(species))
    seg = zeros(Int, length(species))
    for (k, members) in enumerate(segments), i in members
        seg[i] = k
    end
    rho = Dict{Symbol, Int}()
    for s in substrates(rxn); rho[name(s)] = get(rho, name(s), 0) + 1; end
    for p in products(rxn);   rho[name(p)] = get(rho, name(p), 0) - 1; end
    offset(i) = sum(get(rho, x, 0) * e for (x, e) in extras[i]; init = 0)
    weight(s) = sum(get(rho, name(x), 0) for x in consumed(s); init = 0) -
                sum(get(rho, name(x), 0) for x in released(s); init = 0) +
                offset(idx[from_species(s)]) - offset(idx[to_species(s)])
    flags = [falses(length(group)) for group in groups]
    edges = Tuple{Int, Int}[]; weights = Int[]; owner = Tuple{Int, Int}[]
    for (g, group) in enumerate(groups), (j, s) in enumerate(group)
        is_equilibrium(s) && continue
        a, b = seg[idx[from_species(s)]], seg[idx[to_species(s)]]
        if a == b
            flags[g][j] = weight(s) != 0
        else
            push!(edges, (a, b)); push!(weights, weight(s)); push!(owner, (g, j))
        end
    end
    block = _edge_blocks(length(segments), edges)
    members = Dict{Int, Vector{Int}}()
    for e in eachindex(edges)
        push!(get!(members, block[e], Int[]), e)
    end
    unbalanced = Set(b for (b, es) in members if !_block_balanced(edges, weights, es))
    for (e, (g, j)) in enumerate(owner)
        flags[g][j] = block[e] in unbalanced
    end
    flags
end

"""Whether every cycle of the block made of the edges `es` has weight zero. A BFS
spanning tree gives each vertex a potential, rising by an edge's weight along its
stored direction; the block is balanced iff every edge's weight equals the
potential difference of its ends. A block is connected, so one search from any of
its vertices visits all of them."""
function _block_balanced(edges, weights, es::Vector{Int})
    adj = Dict{Int, Vector{Int}}()
    for e in es
        u, v = edges[e]
        push!(get!(adj, u, Int[]), e); push!(get!(adj, v, Int[]), e)
    end
    root = edges[first(es)][1]
    phi = Dict(root => 0)
    queue = [root]
    while !isempty(queue)
        u = popfirst!(queue)
        for e in adj[u]
            a, b = edges[e]
            v = a == u ? b : a
            w = a == u ? weights[e] : -weights[e]
            if haskey(phi, v)
                phi[v] == phi[u] + w || return false
            else
                phi[v] = phi[u] + w
                push!(queue, v)
            end
        end
    end
    true
end

"""
    _flux_carrying_groups(groups, rxn) -> BitVector
    _flux_carrying_groups(m) -> BitVector

One flag per kinetic group: whether some step of the group carries net
steady-state flux (`_flux_carrying_steps`). A steady-state group with no such step
exposes only the ratio of its two constants, so the moves never emit one; a
zero-flux step inside a group that also holds a flux-carrying step costs nothing,
because the group's shared constants are pinned by the step that carries flux.
The mechanism method reads `steps(m)`, which for an allosteric mechanism is its
A-state catalytic graph.
"""
_flux_carrying_groups(groups::Vector{Vector{Step}}, rxn::EnzymeReaction) =
    BitVector([any(f) for f in _flux_carrying_steps(groups, rxn)])
_flux_carrying_groups(m::Union{Mechanism, AllostericMechanism}) =
    _flux_carrying_groups(steps(m), reaction(m))
```

In `_expand_re_to_ss`, change the first line of the body from `flux = _flux_carrying_groups(m)` to

```julia
    flux = _flux_carrying_groups(_all_steady_state(steps(m)), reaction(m))
```

and in its docstring replace the parenthesis "(`_flux_carrying_groups`; a group with none exposes only equilibrium ratios and would gain a phantom parameter)" with "(`_flux_carrying_groups` on the all-steady-state graph, `_all_steady_state`; a group with none exposes only equilibrium ratios under every assignment and would gain a phantom parameter)".

- [ ] **Step 4: Run the testset, then the whole enumeration file**

Scratch run of the testset: PASS. Then the focused run of `test/test_mechanism_enumeration.jl`: every test passes, including "seed child count is unchanged (220 over bi-bi seeds)" and "a dead-end leaf group never flips" — the all-steady-state flags equal the former screen on every enumerated mechanism.

- [ ] **Step 5: Run the full suite, then commit**

Start the full suite in the background with its output in a log, poll the log in a bounded foreground loop until it prints its summary, and commit only when every test passes:

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Test flux on the rapid-equilibrium segment graph

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R4zCpbZSoygD66kRecfrDN"
```

---

### Task 2: The flip extends a zero-flux set instead of emitting it

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_expand_re_to_ss`: the `gains` closure and the docstring, about lines 1230–1281)
- Test: `test/test_mechanism_enumeration.jl` (new testsets at the end of `@testset "_expand_re_to_ss (group-set flips)"`, about line 6257 onward; the testset "split ter-ter pairs the A and B parts", about lines 6630–6688)

**Interfaces:**
- Consumes: `_flux_carrying_groups(groups, rxn)` from Task 1.
- Produces: nothing new; `_expand_re_to_ss(m)` keeps its signature and now never emits a child whose flipped group carries no flux.

- [ ] **Step 1: Write the failing tests**

Append inside `@testset "_expand_re_to_ss (group-set flips)"`, before its closing `end`:

```julia
    @testset "_expand_re_to_ss: a zero-flux pair is not emitted, its supersets stay reachable" begin
        # Ordered bi-bi with the abortive complex E(A, Q), after both the A group and
        # the Q group were split by context: six rapid-equilibrium binding groups,
        # each holding one step. Singles: E(A) + B and E(Q) + P are bridges of the RE
        # graph (E(A, B) and E(P, Q) are leaves there), so each flips alone; the four
        # bindings on the square E–E(A)–E(A, Q)–E(Q) are not. Pairs of square edges
        # all divide the segment. {E + A, E + Q} isolates {E}, joined to the rest by
        # two parallel edges of weight 0 (A uptake +1 against E(A)'s offset +1; Q
        # uptake −1 against E(Q)'s offset −1): a balanced block, so neither flipped
        # group carries flux and the child would be its parent plus two phantoms.
        # {E(A) + Q, E(Q) + A} isolates {E(A, Q)} the same way (weights −1 + 1 and
        # 1 − 1). The other four pairs each cut a block that also holds the chemistry
        # edge, whose weight 3 or 4 against a binding edge of weight 1 makes the block
        # unbalanced. Each of the two zero-flux pairs has no extension: every
        # three-set containing it also contains a gaining pair, so its supersets are
        # reached by a later flip of those children.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(Q) + A ⇌ E(A, Q)
                E + Q ⇌ E(Q)
                E(A) + Q ⇌ E(A, Q)
                E(A) + B ⇌ E(A, B)
                E(A, B) <--> E(P, Q)
                E(Q) + P ⇌ E(P, Q)
            end
        end)
        expected = [
            EnzymeRates.Mechanism(@enzyme_mechanism begin      # {E(A) + B}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ E(A)
                    E(Q) + A ⇌ E(A, Q)
                    E + Q ⇌ E(Q)
                    E(A) + Q ⇌ E(A, Q)
                    E(A) + B <--> E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(Q) + P ⇌ E(P, Q)
                end
            end),
            EnzymeRates.Mechanism(@enzyme_mechanism begin      # {E(Q) + P}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ E(A)
                    E(Q) + A ⇌ E(A, Q)
                    E + Q ⇌ E(Q)
                    E(A) + Q ⇌ E(A, Q)
                    E(A) + B ⇌ E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(Q) + P <--> E(P, Q)
                end
            end),
            EnzymeRates.Mechanism(@enzyme_mechanism begin      # {E + A, E(A) + Q}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> E(A)
                    E(Q) + A ⇌ E(A, Q)
                    E + Q ⇌ E(Q)
                    E(A) + Q <--> E(A, Q)
                    E(A) + B ⇌ E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(Q) + P ⇌ E(P, Q)
                end
            end),
            EnzymeRates.Mechanism(@enzyme_mechanism begin      # {E + A, E(Q) + A}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> E(A)
                    E(Q) + A <--> E(A, Q)
                    E + Q ⇌ E(Q)
                    E(A) + Q ⇌ E(A, Q)
                    E(A) + B ⇌ E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(Q) + P ⇌ E(P, Q)
                end
            end),
            EnzymeRates.Mechanism(@enzyme_mechanism begin      # {E + Q, E(A) + Q}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ E(A)
                    E(Q) + A ⇌ E(A, Q)
                    E + Q <--> E(Q)
                    E(A) + Q <--> E(A, Q)
                    E(A) + B ⇌ E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(Q) + P ⇌ E(P, Q)
                end
            end),
            EnzymeRates.Mechanism(@enzyme_mechanism begin      # {E + Q, E(Q) + A}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ E(A)
                    E(Q) + A <--> E(A, Q)
                    E + Q <--> E(Q)
                    E(A) + Q ⇌ E(A, Q)
                    E(A) + B ⇌ E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(Q) + P ⇌ E(P, Q)
                end
            end),
        ]
        kids = EnzymeRates._expand_re_to_ss(m)
        @test length(kids) == 6
        @test Set(kids) == Set(expected)
        # The two zero-flux pairs are absent, and each has exactly its parent's rank.
        r0 = _testhelper_identifiable_rank(m)
        for absent in (
            EnzymeRates.Mechanism(@enzyme_mechanism begin  # {E + A, E + Q}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> E(A)
                    E(Q) + A ⇌ E(A, Q)
                    E + Q <--> E(Q)
                    E(A) + Q ⇌ E(A, Q)
                    E(A) + B ⇌ E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(Q) + P ⇌ E(P, Q)
                end
            end),
            EnzymeRates.Mechanism(@enzyme_mechanism begin  # {E(A) + Q, E(Q) + A}
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ E(A)
                    E(Q) + A <--> E(A, Q)
                    E + Q ⇌ E(Q)
                    E(A) + Q <--> E(A, Q)
                    E(A) + B ⇌ E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(Q) + P ⇌ E(P, Q)
                end
            end))
            @test !(absent in kids)
            @test EnzymeRates._re_segment_count(absent) > EnzymeRates._re_segment_count(m)
            @test _testhelper_identifiable_rank(absent) == r0
        end
    end

    @testset "_expand_re_to_ss: a Theorell–Chance parent flips each binding" begin
        # No isomerization step: the former screen found no unit here. Every step
        # lies on the unbalanced cycle E → E(A) → E(Q) → E (weights 1, 2, 1), so both
        # rapid-equilibrium bindings are units; each divides the one segment alone.
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) + B <--> E(Q) + P
                E + Q ⇌ E(Q)
            end
        end)
        expected = [
            EnzymeRates.Mechanism(@enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> E(A)
                    E(A) + B <--> E(Q) + P
                    E + Q ⇌ E(Q)
                end
            end),
            EnzymeRates.Mechanism(@enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ E(A)
                    E(A) + B <--> E(Q) + P
                    E + Q <--> E(Q)
                end
            end),
        ]
        kids = EnzymeRates._expand_re_to_ss(m)
        @test length(kids) == 2
        @test Set(kids) == Set(expected)
    end
```

Then rewrite the testset "_expand_re_to_ss: split ter-ter pairs the A and B parts". Replace its leading comment with:

```julia
        # Ter-ter after one context split: A is split by whether B is bound and B
        # by whether A is bound. Each of the four parts alone leaves E and E(A)
        # joined through its sibling, so nothing flips one part by itself; C, P, Q
        # and R still flip alone. Of the six pairs, {A1, B1} cuts both binding edges
        # at E and E(C) and leaves {E(A), E(B), E(A, B), E(A, C), E(B, C), E(A, B, C)}
        # as a segment with no substrate-free form, which the constructor rejects.
        # {A1, B2} isolates {E(A), E(A, C)}: its four edges to the rest all have
        # weight 1 inward (A uptake) and −1 outward (B uptake against the +2 offset
        # of E(A, B)), so every cycle through them is balanced and the pair carries
        # no flux; the turnover runs through E(B) inside the big segment. {A2, B1}
        # isolates {E(B), E(B, C)} the same way. Both count as failed and have no
        # extension, since every three-set containing them holds a gaining pair.
        # The other three pairs cut a block that also holds the chemistry edge
        # (weight 6 against edges of weight 2) and are emitted.
```

and replace the `expected` list, the count and the loop with:

```julia
        expected = [
            flip(C),        # {C}
            flip(P),        # {P}
            flip(Q),        # {Q}
            flip(R),        # {R}
            flip(A1, A2),   # {A1, A2}
            flip(B1, B2),   # {B1, B2}
            flip(A2, B2),   # {A2, B2}
        ]
        kids = EnzymeRates._expand_re_to_ss(m)
        @test length(kids) == 7
        @test Set(kids) == Set(expected)
        for c in kids
            ss = count((A1, A2, B1, B2)) do key
                grp = only(grp for grp in EnzymeRates.steps(c)
                           if (binds(grp), sources(grp)) == key)
                !any(EnzymeRates.is_equilibrium, grp)
            end
            @test ss != 1
        end
        # The two zero-flux pairs divide a segment yet keep the parent's rank.
        r0 = _testhelper_identifiable_rank(m)
        for absent in (flip(A1, B2), flip(A2, B1))
            @test !(absent in kids)
            @test EnzymeRates._re_segment_count(absent) > EnzymeRates._re_segment_count(m)
            @test _testhelper_identifiable_rank(absent) == r0
        end
```

- [ ] **Step 2: Run the three testsets and watch them fail**

Scratch run (the group-set flips section needs `_testhelper_flip_groups` and `_testhelper_identifiable_rank`, defined earlier in the file at about lines 5941 and 6213; copy them into the scratch file). Expected: the zero-flux pair test fails with 8 children instead of 6; the Theorell–Chance test passes already (Task 1 gave it units); the ter-ter test fails with 9 children instead of 7.

- [ ] **Step 3: Extend zero-flux sets in `gains`**

In `_expand_re_to_ss`, replace the `gains` closure

```julia
    gains(sel) =
        _re_segment_count_after_flip(m, Set(units[u] for u in sel)) > base &&
        _bottomless_re_segment(flipped_groups(sel)) === nothing
```

with

```julia
    gains(sel) = begin
        _re_segment_count_after_flip(m, Set(units[u] for u in sel)) > base || return false
        groups = flipped_groups(sel)
        _bottomless_re_segment(groups) === nothing || return false
        flux = _flux_carrying_groups(groups, reaction(m))
        all(u -> flux[units[u]], sel)
    end
```

In the docstring, after the sentence ending "...so the minimal-set search extends it instead of emitting it;" and before "the non-degenerate mechanisms beyond it stay reachable", insert: "so does a set one of whose flipped groups carries no flux in the child (`_flux_carrying_groups` on the flipped groups): that group's constants enter the rate only as their ratio, the child is its parent plus one phantom per such group, and its flux-carrying supersets are the smallest identifiable relaxations;". Read the resulting sentence once and make it grammatical.

- [ ] **Step 4: Run the testsets, then the whole enumeration file**

Scratch run: all three PASS. Focused run of `test/test_mechanism_enumeration.jl`: every test passes. "seed child count is unchanged (220 over bi-bi seeds)" must still hold (no seed flip is zero-flux). If any other exact-children flip test fails, its parent is a split child with a zero-flux pair: derive the balanced block by hand as in the comments above, remove that child from `expected`, add the rank check, and record the case in the task report.

- [ ] **Step 5: Run the full suite, then commit**

Start the full suite in the background with its output in a log, poll the log in a bounded foreground loop until it prints its summary, and commit only when every test passes:

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Extend a flip set whose flipped group carries no flux

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R4zCpbZSoygD66kRecfrDN"
```

---

### Task 3: The twin predicate and the dead-end move's new-complex rule

> Amended during execution (fix round 1, controller ruling): the twin test is judged in every conformational state where the copy binds. `_copy_twin_test(m)` wraps `_twin_site_test(groups)` per state, `_duplicate_copy_groups` takes the mechanism, and the dead-end move passes the copy's tag `:EqualAI`. Tasks 4–6 below consume the amended interface.

**Files:**
- Modify: `src/mechanism_enumeration.jl` (new functions next to `_bound_at_forms`, about line 1720; `_expand_add_dead_end_regulator_native`, about lines 1829–1960; the `_expand_add_dead_end_regulator` docstring, about lines 1766–1790)
- Test: `test/test_mechanism_enumeration.jl` (new `@testset "_twin_site_test"` after `@testset "_flux_carrying_groups"`, about line 5543; new and rewritten testsets inside `@testset "_expand_add_dead_end_regulator"`, about lines 2343–2720, replacing "Mechanism — Substrate-as-dead-end-inhibitor overlap" at about line 2655; rewritten flip testsets at about lines 1870–1961 and 1963–1978; the `_hyperbolic_catalysis` assertions at about lines 5748–5766)

**Interfaces:**
- Produces (used by Tasks 4 and 5):
  - `_composition(sp::Species, extra::Metabolite...)` — `(sorted bound names plus extra names, conformation, residual)`.
  - `_twin_site_test(groups::Vector{Vector{Step}})` — returns a closure `(site::Species, ligand::Metabolite) -> Bool`.
  - `_duplicate_copy_groups(groups::Vector{Vector{Step}})::Vector{Int}` — indices of kinetic groups that bind a competitive inhibitor only at twin sites.

- [ ] **Step 1: Write the failing tests**

**(a) The predicate.** After `@testset "_flux_carrying_groups"` add:

```julia
@testset "_twin_site_test" begin
    ER = EnzymeRates
    # The form named `nm` among the ends of `m`'s steps.
    form(m, nm) = only(unique(sp for g in ER.steps(m) for s in g
                              for sp in (ER.from_species(s), ER.to_species(s))
                              if ER.name(sp) == nm))
    Ainh, Binh, Sinh, Qinh, I = ER.CompetitiveInhibitor(:A), ER.CompetitiveInhibitor(:B),
        ER.CompetitiveInhibitor(:S), ER.CompetitiveInhibitor(:Q), ER.CompetitiveInhibitor(:I)

    # Ordered bi-bi. A copy of A at E has the composition of E(A); at E(Q) it is new;
    # a copy of B at E is new; a foreign inhibitor is never a twin.
    ordered = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    twin = ER._twin_site_test(ER.steps(ordered))
    @test twin(form(ordered, :E), Ainh)
    @test !twin(form(ordered, :EQ), Ainh)
    @test !twin(form(ordered, :E), Binh)
    @test !twin(form(ordered, :E), I) && !twin(form(ordered, :EQ), I)

    # With the abortive complex E(A, Q) present, the copy at E(Q) is a twin too.
    abortive = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))
        end
    end)
    @test ER._twin_site_test(ER.steps(abortive))(form(abortive, :EQ), Ainh)

    # A copy whose complex already exists is judged against the other forms only.
    with_copy = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    twin_c = ER._twin_site_test(ER.steps(with_copy))
    @test twin_c(form(with_copy, :E), Ainh)
    @test !twin_c(form(with_copy, :EQ), Ainh)

    # The offsets key. E ⇌ E* is a rapid-equilibrium isomerization, so E and E*
    # have equal offsets; a copy of S at E has a composition no form has, yet its
    # weight is proportional to E*(S)'s. With the isomerization at steady state the
    # two forms sit in different segments and the copy is not a twin.
    iso = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E ⇌ Estar
            Estar + S ⇌ Estar(S)
            Estar(S) <--> E(P)
            E + P ⇌ E(P)
        end
    end)
    estar_s = only(ER.to_species(s) for g in ER.steps(iso) for s in g
                   if ER.bound_metabolite(s) !== nothing &&
                      ER.name(ER.bound_metabolite(s)) == :S)
    @test ER._twin_site_test(ER.steps(iso))(form(iso, :E), Sinh)
    @test ER._composition(form(iso, :E), Sinh) != ER._composition(estar_s)
    iso_ss = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E <--> Estar
            Estar + S ⇌ Estar(S)
            Estar(S) <--> E(P)
            E + P ⇌ E(P)
        end
    end)
    @test !ER._twin_site_test(ER.steps(iso_ss))(form(iso_ss, :E), Sinh)

    # Ping-pong with the second chemistry step at rapid equilibrium, an abortive
    # E(B, P; residual) and a copy of P at E. A copy of Q at E(P::Inh) has the new
    # composition {P, Q}, but E(B; residual) ⇌ E(Q) puts E(B, P; residual) at the
    # offsets P + Q of E's segment, the same as E(P::Inh)·Q: a twin. At
    # E(; residual) the copy is not.
    pp = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) ⇌ E(Q)
            E(Q) ⇌ E + Q
            E(B; residual = A - P) + P ⇌ E(B, P; residual = A - P)
            E + P::Inh ⇌ E(P::Inh)
        end
    end)
    copy_form = only(ER.to_species(s) for g in ER.steps(pp) for s in g
                     if ER.bound_metabolite(s) isa ER.CompetitiveInhibitor)
    free_res = only(ER.from_species(s) for g in ER.steps(pp) for s in g
                    if ER.bound_metabolite(s) !== nothing &&
                       ER.name(ER.bound_metabolite(s)) == :B)
    twin_pp = ER._twin_site_test(ER.steps(pp))
    @test twin_pp(copy_form, Qinh)
    @test !twin_pp(free_res, Qinh)

    # Groups that bind a copy only at twin sites.
    @test isempty(ER._duplicate_copy_groups(ER.steps(with_copy)))
    twin_only = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            E + A::Inh ⇌ E(A::Inh)
        end
    end)
    dup = ER._duplicate_copy_groups(ER.steps(twin_only))
    @test length(dup) == 1
    @test ER.bound_metabolite(first(ER.steps(twin_only)[only(dup)])) == Ainh
end
```

**(b) The dead-end move.** Inside `@testset "_expand_add_dead_end_regulator"`, replace the testset "Mechanism — Substrate-as-dead-end-inhibitor overlap" (about line 2655) with these four testsets:

```julia
@testset "Mechanism — uni-uni: a substrate or product as its own inhibitor has no placement" begin
    # The only site that carries neither S nor P is free E, and E·S* has the
    # composition of E(S) (E·P* that of E(P)): the copy's constant would enter the
    # rate only as K_S + K_S,inh with kcat rescaled. No child is emitted.
    rxn_S = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: S
    end
    rxn_P = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: P
    end
    for rxn in (rxn_S, rxn_P)
        m = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E + P ⇌ E(P)
            end
        end)
        @test isempty(EnzymeRates._expand_add_dead_end_regulator(m, rxn))
    end
end

@testset "Mechanism — ordered bi-bi with A as its own inhibitor: three placements" begin
    # A binds productively at E only, so A* competing with P puts it at {E, E(Q)},
    # with Q puts it at {E} alone, and competing with B puts it at {E(A), E(Q)} or
    # {E, E(A)}. E·A* has the composition of E(A); E(Q)·A* and E(A)·A* are new. The
    # placement at E alone is all-twin and is not emitted. Mirrors: E + Q ⇌ E(Q)
    # onto {E, E(Q)} in the Q group; E + A ⇌ E(A) onto {E, E(A)} in the A group.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: A
    end
    # A dead-end child carries the declared reaction (atoms and the inhibitor), so
    # every fixture is rebuilt on it for `==`.
    withrxn(em) = EnzymeRates.Mechanism(rxn, EnzymeRates.steps(EnzymeRates.Mechanism(em)))
    m = withrxn(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    at_E_EQ = withrxn(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    at_EA_EQ = withrxn(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            (E(A) + A::Inh ⇌ E(A, A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    at_E_EA = withrxn(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(A::Inh) + A ⇌ E(A, A::Inh))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            (E + A::Inh ⇌ E(A::Inh), E(A) + A::Inh ⇌ E(A, A::Inh))
        end
    end)
    at_E = withrxn(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            E + A::Inh ⇌ E(A::Inh)
        end
    end)
    kids = EnzymeRates._expand_add_dead_end_regulator(m, rxn)
    @test length(kids) == 3
    @test Set(kids) == Set([at_E_EQ, at_EA_EQ, at_E_EA])
    @test !(at_E in kids)
    # Each child binds the copy as a CompetitiveInhibitor named :A, in one new
    # kinetic group, and adds exactly one fitted constant.
    base_fitted = length(EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(m)))
    for r in kids
        @test any(grp -> (bm = EnzymeRates.bound_metabolite(first(grp));
                          bm isa EnzymeRates.CompetitiveInhibitor && EnzymeRates.name(bm) === :A),
                  EnzymeRates.steps(r))
        @test length(EnzymeRates.steps(r)) == length(EnzymeRates.steps(m)) + 1
        @test length(EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(r))) ==
            base_fitted + 1
    end
    # The absent placement adds a constant the data cannot see.
    @test _testhelper_identifiable_rank(at_E) == _testhelper_identifiable_rank(m)
end

@testset "Mechanism — shared A group with abortive E(A, Q): twin-only placements absent" begin
    # A binds at E and, in the same kinetic group, at E(Q) to form the abortive
    # complex E(A, Q). A copy of A at E duplicates E(A) and at E(Q) duplicates
    # E(A, Q), so the placements {E} and {E, E(Q)} are all-twin. {E, E(Q)} is a
    # phantom. {E} is identifiable, by accident: the shared group pins K_A at E(Q),
    # where no copy binds, so the absorption of the copy's constant is blocked. A
    # non-productive orientation of E(A) is not a distinct hypothesis, so it is
    # dropped with the rest (Denis, 2026-09-30). {E(A), E(Q)} and {E, E(A)} create
    # E(A, A*) and are emitted; the second mirrors E + A ⇌ E(A) in the A group.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: A
    end
    # A dead-end child carries the declared reaction (atoms and the inhibitor), so
    # every fixture is rebuilt on it for `==`.
    withrxn(em) = EnzymeRates.Mechanism(rxn, EnzymeRates.steps(EnzymeRates.Mechanism(em)))
    m = withrxn(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))
        end
    end)
    at_EA_EQ = withrxn(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))
            (E(A) + A::Inh ⇌ E(A, A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    at_E_EA = withrxn(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q), E(A::Inh) + A ⇌ E(A, A::Inh))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))
            (E + A::Inh ⇌ E(A::Inh), E(A) + A::Inh ⇌ E(A, A::Inh))
        end
    end)
    at_E = withrxn(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))
            E + A::Inh ⇌ E(A::Inh)
        end
    end)
    at_E_EQ = withrxn(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    kids = EnzymeRates._expand_add_dead_end_regulator(m, rxn)
    @test length(kids) == 2
    @test Set(kids) == Set([at_EA_EQ, at_E_EA])
    r0 = _testhelper_identifiable_rank(m)
    @test _testhelper_identifiable_rank(at_E_EQ) == r0        # a phantom
    @test _testhelper_identifiable_rank(at_E) == r0 + 1       # tie-only, dropped anyway
end

@testset "Mechanism — a copy that only matches a form through a conformational isomer" begin
    # E ⇌ E* at rapid equilibrium; S binds E* only, P binds E only. S* competing
    # with P is placed at E (where P binds) and E* (where S binds), with the
    # isomerization mirrored between the copy forms. E*·S* has E*(S)'s composition.
    # E·S* has a composition no form has, but E and E* share offsets, so E·S* and
    # E*(S) have proportional weights: the pattern is all-twin and is not emitted;
    # the would-be child has its parent's rank.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: S
    end
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E ⇌ Estar
            Estar + S ⇌ Estar(S)
            Estar(S) <--> E(P)
            E + P ⇌ E(P)
        end
    end)
    @test isempty(EnzymeRates._expand_add_dead_end_regulator(m, rxn))
    absent = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E ⇌ Estar, E(S::Inh) ⇌ Estar(S::Inh))
            Estar + S ⇌ Estar(S)
            Estar(S) <--> E(P)
            E + P ⇌ E(P)
            (E + S::Inh ⇌ E(S::Inh), Estar + S::Inh ⇌ Estar(S::Inh))
        end
    end)
    @test _testhelper_identifiable_rank(absent) == _testhelper_identifiable_rank(m)
end
```

`_testhelper_identifiable_rank` is defined at about line 6213, after the sections that now use it. Move its definition, with its docstring, to the file's top-level helper block before the first `@testset` (next to `_connectivity_violations`), and delete it from its old place. The file already has `LinearAlgebra` and `Random` from `runtests.jl`.

**(c) Flip tests that build the uni-uni Case 1 fixture.** Replace the three testsets "Mechanism — Substrate-as-dead-end-inhibitor overlap" (about line 1870, under `@testset "_expand_re_to_ss"`), "Mechanism — Allosteric substrate-as-dead-end-I overlap" (about line 1922) and "_expand_re_to_ss keeps inhibitor bindings RE" (about line 1963) with:

```julia
@testset "Mechanism — a substrate's inhibitor copy never flips; the catalytic groups do" begin
    # Ordered bi-bi with A also bound as a competitive-inhibitor copy at E and
    # E(Q), the Q binding mirrored onto the copy forms in the Q group, as the
    # dead-end move builds it. The copy group binds a regulator and is never a
    # unit. A, B, P and Q each divide the one segment alone (the Q flip takes its
    # mirror with it and isolates {E(Q), E(P, Q), E(A::Inh, Q)}), and each flipped
    # group shares a block with the chemistry edge, so all four carry flux.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    expected = [
        EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(A) + B ⇌ E(A, B)
                E(A, B) <--> E(P, Q)
                E(Q) + P ⇌ E(P, Q)
                (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
                (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            end
        end),
        EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) + B <--> E(A, B)
                E(A, B) <--> E(P, Q)
                E(Q) + P ⇌ E(P, Q)
                (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
                (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            end
        end),
        EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) + B ⇌ E(A, B)
                E(A, B) <--> E(P, Q)
                E(Q) + P <--> E(P, Q)
                (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
                (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            end
        end),
        EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) + B ⇌ E(A, B)
                E(A, B) <--> E(P, Q)
                E(Q) + P ⇌ E(P, Q)
                (E + Q <--> E(Q), E(A::Inh) + Q <--> E(A::Inh, Q))
                (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
            end
        end),
    ]
    kids = EnzymeRates._expand_re_to_ss(m)
    @test length(kids) == 4
    @test Set(kids) == Set(expected)
    for r in kids, grp in EnzymeRates.steps(r), s in grp
        EnzymeRates.bound_metabolite(s) isa EnzymeRates.Regulator &&
            @test EnzymeRates.is_equilibrium(s)
    end
end

@testset "AllostericMechanism — the inhibitor copy never flips; tags are preserved" begin
    # The same mechanism as an allosteric one with a dead inactive conformation
    # (the chemistry `:OnlyA`, every binding `:EqualAI`) and one catalytic subunit.
    # The same four groups flip; every child keeps one `:OnlyA` group, the
    # chemistry, and the multiplicity and (empty) regulatory sites.
    am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        catalytic_inhibitors: A
        catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A ⇌ E(A)                                            :: EqualAI
            E(A) + B ⇌ E(A, B)                                      :: EqualAI
            E(A, B) <--> E(P, Q)                                    :: OnlyA
            E(Q) + P ⇌ E(P, Q)                                      :: EqualAI
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))            :: EqualAI
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))  :: EqualAI
        end
    end)
    flip(rewrite) = EnzymeRates.AllostericMechanism(rewrite)
    expected = [
        flip(@allosteric_mechanism begin
            substrates: A, B
            products: P, Q
            catalytic_inhibitors: A
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + A <--> E(A)                                         :: EqualAI
                E(A) + B ⇌ E(A, B)                                      :: EqualAI
                E(A, B) <--> E(P, Q)                                    :: OnlyA
                E(Q) + P ⇌ E(P, Q)                                      :: EqualAI
                (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))            :: EqualAI
                (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))  :: EqualAI
            end
        end),
        flip(@allosteric_mechanism begin
            substrates: A, B
            products: P, Q
            catalytic_inhibitors: A
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + A ⇌ E(A)                                            :: EqualAI
                E(A) + B <--> E(A, B)                                   :: EqualAI
                E(A, B) <--> E(P, Q)                                    :: OnlyA
                E(Q) + P ⇌ E(P, Q)                                      :: EqualAI
                (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))            :: EqualAI
                (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))  :: EqualAI
            end
        end),
        flip(@allosteric_mechanism begin
            substrates: A, B
            products: P, Q
            catalytic_inhibitors: A
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + A ⇌ E(A)                                            :: EqualAI
                E(A) + B ⇌ E(A, B)                                      :: EqualAI
                E(A, B) <--> E(P, Q)                                    :: OnlyA
                E(Q) + P <--> E(P, Q)                                   :: EqualAI
                (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))            :: EqualAI
                (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))  :: EqualAI
            end
        end),
        flip(@allosteric_mechanism begin
            substrates: A, B
            products: P, Q
            catalytic_inhibitors: A
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + A ⇌ E(A)                                            :: EqualAI
                E(A) + B ⇌ E(A, B)                                      :: EqualAI
                E(A, B) <--> E(P, Q)                                    :: OnlyA
                E(Q) + P ⇌ E(P, Q)                                      :: EqualAI
                (E + Q <--> E(Q), E(A::Inh) + Q <--> E(A::Inh, Q))      :: EqualAI
                (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))  :: EqualAI
            end
        end),
    ]
    kids = EnzymeRates._expand_re_to_ss(am)
    @test length(kids) == 4
    @test Set(kids) == Set(expected)
    for r in kids
        @test count(==(:OnlyA), r.cat_allo_states) == 1
        @test EnzymeRates.is_iso(first(r.cat_steps[findfirst(==(:OnlyA), r.cat_allo_states)]))
        @test r.catalytic_multiplicity == am.catalytic_multiplicity
        @test r.regulatory_sites == am.regulatory_sites
        @test EnzymeRates.reaction(r) == EnzymeRates.reaction(am)
    end
end
```

If `@allosteric_mechanism` rejects `catalytic_inhibitors: A` for a name that is also a substrate, build each allosteric fixture from the corresponding inline `@enzyme_mechanism` instead: `EnzymeRates.AllostericMechanism(EnzymeRates.reaction(mm), EnzymeRates.steps(mm), tags, 1, EnzymeRates.RegulatorySite[])` with `tags = [EnzymeRates.is_iso(first(g)) ? :OnlyA : :EqualAI for g in EnzymeRates.steps(mm)]`, and note the DSL gap in the task report.

**(d) The `_hyperbolic_catalysis` counts with A as its own inhibitor** (about lines 5751–5766). Replace the block from the comment "A substrate declared as a dead-end inhibitor..." through `@test !any(EnzymeRates._hyperbolic_catalysis, dead_end_with_a)` with:

```julia
    # A substrate declared as a dead-end inhibitor binds its inhibitor site under
    # the substrate's own name. Its bindings put A² and A³ into the derived
    # denominator, but those powers come from the inhibitor site, which the
    # predicate leaves out: the three placements on the ordered scheme pass. The
    # fourth placement, the copy at E alone, duplicates E(A) and is not emitted.
    rxn_a_inhibits = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: A
    end
    copy_sites(k) = Set(EnzymeRates.name(EnzymeRates.from_species(s))
                        for grp in EnzymeRates.steps(k) for s in grp
                        if EnzymeRates.bound_metabolite(s) isa EnzymeRates.CompetitiveInhibitor)
    with_a = EnzymeRates._expand_add_dead_end_regulator(ordered_ss, rxn_a_inhibits)
    den_a_degree(k) = maximum(get(Dict(mono), :A, 0) for mono in keys(
        EnzymeRates._raw_symbolic_rate_polys(k, EnzymeRates._step_parameters(k),
            EnzymeRates._build_wegscheider_rename_map(k))[2]))
    @test Set(copy_sites.(with_a)) ==
        Set([Set([:E, :EQ]), Set([:EA, :EQ]), Set([:E, :EA])])
    @test all(EnzymeRates._hyperbolic_catalysis, with_a)
    @test sort(den_a_degree.(with_a)) == DEGREES

    # Inhibitor-bound forms leave, but catalytic-site powers stay: adding the
    # inhibitor role of A to the abortive-complex scheme keeps it non-hyperbolic.
    # The placements at E and at {E, E(Q)} duplicate E(A) and E(A, Q) and are not
    # emitted.
    dead_end_with_a = EnzymeRates._expand_add_dead_end_regulator(dead_end, rxn_a_inhibits)
    @test Set(copy_sites.(dead_end_with_a)) == Set([Set([:EA, :EQ]), Set([:E, :EA])])
    @test !any(EnzymeRates._hyperbolic_catalysis, dead_end_with_a)
```

`DEGREES` is the sorted list of the three survivors' denominator degrees in A. It is a derivation fact, not an enumeration expectation: after Step 3, evaluate `sort(den_a_degree.(with_a))` once in the scratch session, confirm it is the former `[1, 2, 2, 3]` minus one entry, write the literal in place of `DEGREES`, and say in the comment which placement's degree left.

- [ ] **Step 2: Run the new and rewritten testsets and watch them fail**

Scratch run of the `_twin_site_test` testset: `UndefVarError: _twin_site_test`. Scratch run of the dead-end testsets: the uni-uni test fails (one child each), the ordered test fails (4 children), the shared-group test fails (4 children), the isomer test fails (1 child). The flip rewrites pass already (they never depended on the rule). The hyperbolic block fails on the site sets (4 placements each).

- [ ] **Step 3: Implement the predicate and the rule**

In `src/mechanism_enumeration.jl`, after `_bound_at_forms` (about line 1738), add:

```julia
"""
Composition of a form with `extra` bound as well: the sorted names of its bound
metabolites, a competitive-inhibitor copy counting as the reactant it copies, with
its conformation and residual. Two forms of one composition hold the same ligands
in different ways, which steady-state data cannot tell apart.
"""
_composition(sp::Species, extra::Metabolite...) =
    (sort!(vcat(Symbol[name(b) for b in bound(sp)], Symbol[name(x) for x in extra])),
     conformation(sp), residual(sp))

"""
    _twin_site_test(groups) -> (site::Species, ligand::Metabolite) -> Bool

Whether binding `ligand` at `site` by a rapid-equilibrium step gives a complex that
duplicates another form of `groups`: a form of the same composition, or a form of
the same rapid-equilibrium segment with the same offsets (`_re_segment_extras`),
whose weight is then the complex's weight times a constant. Steady-state data see
only the sum of two such weights, so the binding's constant enters the rate only
through that sum. The complex itself, when `groups` already holds it, is not its
own twin. Both keys are needed: a copy of a metabolite bound at steady state has a
twin by composition alone, and a copy across a rapid-equilibrium isomerization
(a ping-pong second chemistry step, a conformational isomer) by offsets alone.
"""
function _twin_site_test(groups::Vector{Vector{Step}})
    species, segments, extras = _re_segment_extras(groups)
    idx = Dict(sp => i for (i, sp) in enumerate(species))
    seg = zeros(Int, length(species))
    for (k, members) in enumerate(segments), i in members
        seg[i] = k
    end
    compositions = Dict{Any, Int}()
    for sp in species
        key = _composition(sp)
        compositions[key] = get(compositions, key, 0) + 1
    end
    function (site::Species, ligand::Metabolite)
        complex = Species(Metabolite[bound(site)..., ligand], conformation(site),
                          residual(site))
        own = get(idx, complex, 0)
        get(compositions, _composition(complex), 0) - (own == 0 ? 0 : 1) > 0 && return true
        i = idx[site]
        target = mergewith(+, extras[i], Dict(name(ligand) => 1))
        any(j -> j != own && extras[j] == target, segments[seg[i]])
    end
end

"""
Kinetic groups of `groups` that bind a competitive inhibitor only at twin sites
(`_twin_site_test`). Such a group's constant enters the rate only through a sum
with an existing constant, or names a second orientation of an existing complex
that a shared group happens to pin; neither is a hypothesis the moves emit.
"""
function _duplicate_copy_groups(groups::Vector{Vector{Step}})
    twin = _twin_site_test(groups)
    [g for (g, group) in enumerate(groups)
     if bound_metabolite(first(group)) isa CompetitiveInhibitor &&
        all(s -> twin(from_species(s), bound_metabolite(s)), group)]
end
```

In `_expand_add_dead_end_regulator_native`, after `boundmap = _bound_at_forms(m)` add

```julia
    twin = _twin_site_test(steps(m))
```

and inside the pattern loop, right after `isempty(active) && continue`, add

```julia
            all(f -> twin(form_sp[f], CompetitiveInhibitor(reg_name)), active) && continue
```

In the `_expand_add_dead_end_regulator` docstring (the one that starts "Add a dead-end regulator binding step set."), after the sentence ending "(one K_R parameter).", add: "A pattern whose every site is a twin (`_twin_site_test`: the copy's complex has the composition of an existing form, or its segment and offsets) is skipped. Such a copy is a second orientation of a complex the mechanism already has; its constant enters the rate only through a sum with the existing binding's, or, when a shared kinetic group happens to pin that binding, names a hypothesis no different from the existing complex."

- [ ] **Step 4: Run the testsets, fill in `DEGREES`, then the whole enumeration file**

Scratch runs: every new and rewritten testset passes. Focused run of `test/test_mechanism_enumeration.jl`: every test passes. Tests of foreign inhibitors ("Sequential bi-bi + I: 4 variants", "Bi-bi random + I: 9 variants", "Two regulators competition: 17 variants" and the others under `_expand_add_dead_end_regulator`) must be unchanged: a foreign inhibitor's complex carries a name no form has. If one changes, stop and report; the rule is wrong, not the test.

- [ ] **Step 5: Run the full suite, then commit**

Start the full suite in the background with its output in a log, poll the log in a bounded foreground loop until it prints its summary, and commit only when every test passes:

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Skip a dead-end copy that creates no new complex

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R4zCpbZSoygD66kRecfrDN"
```

---

### Task 4: The split reverts a zero-flux part and never isolates twin-only copy sites

**Files:**
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl` (`_partition_independent_count`, about lines 320–365; new `_count_kind`)
- Modify: `src/mechanism_enumeration.jl` (`_expand_split_kinetic_group` and its docstring, about lines 1543–1587; `_split_gain_test(m::Mechanism, units)`, about lines 1589–1607; new `_revert_zero_flux_parts` next to `_apply_bipartitions`)
- Test: `test/test_mechanism_enumeration.jl` (new testsets inside `@testset "_expand_split_kinetic_group"`, about lines 2157–2341, before its closing `end`; a new testset next to the `_partition_independent_count` testset, about line 6195)

**Interfaces:**
- Consumes: `_flux_carrying_steps(groups, rxn)` (Task 1), `_copy_twin_test(m)` (Task 3, as amended in its fix round: a closure `(site, ligand, tag) -> Bool` that judges a copy in every conformational state where it binds).
- Produces:
  - `_count_kind(s::Step)::Symbol` — `:ss`, `:binding_K` or `:iso_K`.
  - `_partition_independent_count(parent)` now returns `counter(group_of_step, kind_of_step = parent kinds)`.
  - `_revert_zero_flux_parts(group::Vector{Step}, bp, flags::BitVector)` — the bipartition with every zero-flux part of a steady-state group rebuilt at rapid equilibrium; returns `bp` itself when nothing changes.

- [ ] **Step 1: Write the failing tests**

Inside `@testset "_expand_split_kinetic_group"`, before its closing `end`, add:

```julia
@testset "Mechanism — a split part with no flux-carrying step is emitted at rapid equilibrium" begin
    # Ordered bi-bi whose steady-state B group holds the catalytic binding
    # E(A) + B and the abortive binding E(Q) + B. Splitting by context leaves
    # E(Q) + B → E(B, Q) alone: E(B, Q) is a dead end, so the step is a pendant
    # bridge of the segment graph and carries no flux. Its two constants would
    # enter the rate only as their ratio, so the part is emitted at rapid
    # equilibrium, one constant fewer than the raw split and the same family.
    # E(B, Q) lies on no cycle, so the new constant raises the independent count: one
    # child. (Amended during execution: the reverted child is 7 fitted / rank 6 against
    # the raw child's 8 / 6 and the parent's 6 / 6; see the assertions below.)
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            (E(A) + B <--> E(A, B), E(Q) + B <--> E(B, Q))
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    reverted = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(A, B)
            E(Q) + B ⇌ E(B, Q)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    raw = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(A, B)
            E(Q) + B <--> E(B, Q)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    kids = EnzymeRates._expand_split_kinetic_group(m)
    @test length(kids) == 1
    @test Set(kids) == Set([reverted])
    @test !(raw in kids)
    @test EnzymeRates._independent_param_count(reverted) ==
        EnzymeRates._independent_param_count(m) + 1
    fitted(k) = length(EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(k)))
    r_rev, r_raw = _testhelper_identifiable_rank(reverted), _testhelper_identifiable_rank(raw)
    @test r_rev == r_raw                                          # the same family
    @test fitted(reverted) == fitted(raw) - 1                     # one constant fewer
    # Both children keep one phantom of another class: once the abortive step no
    # longer pins kf_B/kr_B, E(A) + B → E(A, B) is a steady-state binding into a
    # form with one exit, and its three constants enter the law through two
    # combinations (the chain class sub-project C merges). The revert removes
    # exactly the zero-flux phantom.
    @test fitted(raw) - r_raw == fitted(reverted) - r_rev + 1
end

@testset "Mechanism — a bipartition that isolates twin-only copy sites is not a unit" begin
    # Ordered bi-bi with B also bound as a competitive-inhibitor copy at E(A) and
    # E(Q). E(A, B::Inh) has the composition of E(A, B); E(Q, B::Inh) is new, so
    # the placement stands. Splitting the copy group by context would leave
    # E(A) + B::Inh alone, a group that binds only where it duplicates a form,
    # and the child would be its parent plus a phantom. Every other group holds
    # one step, so the move emits nothing.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            (E(A) + B::Inh ⇌ E(A, B::Inh), E(Q) + B::Inh ⇌ E(B::Inh, Q))
        end
    end)
    @test isempty(EnzymeRates._expand_split_kinetic_group(m))
    absent = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            E(A) + B::Inh ⇌ E(A, B::Inh)
            E(Q) + B::Inh ⇌ E(B::Inh, Q)
        end
    end)
    @test EnzymeRates._independent_param_count(absent) ==
        EnzymeRates._independent_param_count(m) + 1
    @test _testhelper_identifiable_rank(absent) == _testhelper_identifiable_rank(m)
end
```

Also search for a small parent whose reverted child the Wegscheider ties absorb, so that the split emits nothing: try the E1 parent of Task 2 with its two Q bindings made one steady-state group and split by context, and the ordered parent with an abortive complex reached by two rapid-equilibrium routes. When one is found, add a testset pinning its exact (empty or reduced) child set with the tie explained in a comment; when none is found, say so in the task report. The spec asks for this search, not for a result.

Next to the `_partition_independent_count` testset (about line 6195) add:

```julia
@testset "_partition_independent_count counts a reverted part as its RE constant" begin
    # The E(Q) + B step of the steady-state B group, relabelled into a new group
    # and counted as a binding K, gives the count of the built child with that
    # step at rapid equilibrium: the cycle basis does not depend on the flags.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            (E(A) + B <--> E(A, B), E(Q) + B <--> E(B, Q))
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    child = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(A, B)
            E(Q) + B ⇌ E(B, Q)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    counter = EnzymeRates._partition_independent_count(m)
    flat = EnzymeRates._flat_steps(m)
    ids = [g for (_, g) in flat]
    kinds = [EnzymeRates._count_kind(s) for (s, _) in flat]
    j = only(j for (j, (s, _)) in enumerate(flat)
             if EnzymeRates.name(EnzymeRates.from_species(s)) == :EQ &&
                EnzymeRates.bound_metabolite(s) !== nothing &&
                EnzymeRates.name(EnzymeRates.bound_metabolite(s)) == :B)
    @test kinds[j] == :ss
    base = counter(ids)
    ids[j] = length(EnzymeRates.steps(m)) + 1
    kinds[j] = :binding_K
    @test counter(ids, kinds) == EnzymeRates._independent_param_count(child)
    @test counter(ids, kinds) == base + 1
end
```

- [ ] **Step 2: Run the three testsets and watch them fail**

Scratch run. Expected: the revert test fails (the raw split child is emitted, not the reverted one); the copy-split test fails (one child emitted); the counter test fails with `UndefVarError: _count_kind` or a `MethodError` on the two-argument counter.

- [ ] **Step 3: Implement**

In `src/thermodynamic_constr_for_rate_eq_derivation.jl`, replace `_partition_independent_count` and its docstring with:

```julia
"""How a step's constants enter the constraint columns: `:ss` (a forward and a
reverse rate), `:binding_K` (a dissociation constant, whose column carries a sign
flip) or `:iso_K` (an equilibrium constant)."""
_count_kind(s::Step) = is_equilibrium(s) ? (is_binding(s) ? :binding_K : :iso_K) : :ss

"""
    _partition_independent_count(parent::Mechanism) -> counter

Return `counter(group_of_step, kind_of_step = kinds of the parent's steps)`, the
independent-parameter count of the mechanism obtained by regrouping `parent`'s
flat steps (in `_flat_steps` order) into the groups labelled by `group_of_step`,
each step counted under the constant kind `kind_of_step` gives it (`_count_kind`).
Regrouping moves no edges and a flag changes none, so the cycle basis of the step
graph (`_thermodynamic_constraints`) is the same for every call and is computed
once here; each call only merges step columns by group and takes the rank. The
split move passes a step's own kind, or `:binding_K`/`:iso_K` for a steady-state
step it reverts to rapid equilibrium. Equals `_independent_param_count` of the
constructed child: the kernel's independent set is the columns minus the pivots,
and folding a single-symbol Wegscheider tie onto its target removes one column and
one rank together, so the count is invariant to the rename. Column sign
conventions match `_assemble_constraints`: a binding K enters with a sign flip, an
iso K without, an SS step contributes `+kf` and `-kr`.
"""
function _partition_independent_count(parent::Mechanism)
    C, _ = _thermodynamic_constraints(parent)
    kinds = [_count_kind(s) for (s, _) in _flat_steps(parent)]
    function counter(group_of_step::AbstractVector{Int},
                     kind_of_step::AbstractVector{Symbol} = kinds)
        length(group_of_step) == length(kinds) ||
            error("group_of_step must label every flat step of the parent")
        length(kind_of_step) == length(kinds) ||
            error("kind_of_step must label every flat step of the parent")
        column = Dict{Tuple{Int, Int}, Int}()
        for (j, g) in enumerate(group_of_step)
            get!(column, (g, 1), length(column) + 1)
            kind_of_step[j] === :ss && get!(column, (g, 2), length(column) + 1)
        end
        A = zeros(Int, size(C, 1), length(column))
        for (j, g) in enumerate(group_of_step), i in axes(C, 1)
            c = C[i, j]
            c == 0 && continue
            if kind_of_step[j] === :ss
                A[i, column[(g, 1)]] += c
                A[i, column[(g, 2)]] -= c
            else
                A[i, column[(g, 1)]] += kind_of_step[j] === :binding_K ? -c : c
            end
        end
        length(column) - length(_rref_partition(A)[1])
    end
    counter
end
```

In `src/mechanism_enumeration.jl`, replace the body of `_expand_split_kinetic_group` with:

```julia
function _expand_split_kinetic_group(m::Union{Mechanism, AllostericMechanism})
    groups = steps(m)
    flux = _flux_carrying_steps(groups, reaction(m))
    twin = _copy_twin_test(m)
    tag(g) = m isa AllostericMechanism ? cat_allo_state(m, g) : :EqualAI
    duplicate_only(part, g) = bound_metabolite(first(part)) isa CompetitiveInhibitor &&
        all(s -> twin(from_species(s), bound_metabolite(s)::Metabolite, tag(g)), part)
    units = Tuple{Int, Tuple{Vector{Step}, Vector{Step}}}[]
    reverted = Bool[]
    for g in kinetic_groups(m), bp in _context_bipartitions(groups[g])
        (duplicate_only(bp[1], g) || duplicate_only(bp[2], g)) && continue
        parts = _revert_zero_flux_parts(groups[g], bp, flux[g])
        push!(units, (g, parts)); push!(reverted, parts !== bp)
    end
    isempty(units) && return typeof(m)[]
    selection(sel) = [units[u] for u in sel]
    gain = _split_gain_test(m, units)
    gains(sel) =
        (!any(u -> reverted[u], sel) ||
         _bottomless_re_segment(_bipartitioned_groups(groups, selection(sel))[1]) ===
         nothing) && gain(sel)
    segments = _group_re_segments(m)
    partners(sel) = begin
        isempty(sel) && return 1:length(units)
        used = Set(units[u][1] for u in sel)
        touched = reduce(union, (segments[units[u][1]] for u in sel))
        [u for u in 1:length(units)
         if !(units[u][1] in used) &&
            !isempty(intersect(segments[units[u][1]], touched))]
    end
    sets = _minimal_gaining_sets(gains, partners)
    typeof(m)[_apply_bipartitions(m, selection(sel)) for sel in sets]
end
```

Extend its docstring: after the sentence ending "...for an `AllostericMechanism`." add two sentences: "A part of a steady-state group none of whose steps carries flux (`_flux_carrying_steps`, computed once on the parent, since a split moves no edge) is emitted at rapid equilibrium (`_revert_zero_flux_parts`): its two constants would enter the rate only as their ratio, and the rapid-equilibrium part is the same family with one constant fewer; the gain test counts the reverted part under its new kind, so a reverted constant the Wegscheider ties pull back is absorbed like any tied split, and a candidate whose reverted groups leave a rapid-equilibrium segment without a bottom form counts as failed. A bipartition of a competitive-inhibitor group in which one part binds only at twin sites (`_copy_twin_test`, judged in every conformational state where the copy binds) is not a unit: that part's constant would be invisible beside the existing bindings, and every superset of the unit recreates it."

Replace `_split_gain_test(m::Mechanism, units)` with:

```julia
function _split_gain_test(m::Mechanism, units)
    counter = _partition_independent_count(m)
    flat = _flat_steps(m)
    reaction_key(s) = (from_species(s), to_species(s), consumed(s), released(s))
    position = Dict(reaction_key(s) => j for (j, (s, _)) in enumerate(flat))
    parent_ids = [g for (_, g) in flat]
    parent_kinds = [_count_kind(s) for (s, _) in flat]
    base = counter(parent_ids, parent_kinds)
    function gains(sel)
        ids = copy(parent_ids)
        kinds = copy(parent_kinds)
        next_id = length(steps(m))
        for u in sel
            next_id += 1
            for s in units[u][2][1]
                kinds[position[reaction_key(s)]] = _count_kind(s)
            end
            for s in units[u][2][2]
                j = position[reaction_key(s)]
                ids[j] = next_id
                kinds[j] = _count_kind(s)
            end
        end
        counter(ids, kinds) > base
    end
    gains
end
```

and change its docstring to: "Gain test for the split move: `sel -> Bool`, true when applying the selected units raises the independent-parameter count above the parent's. A unit's steps are matched to the parent's flat steps by reaction (forms and metabolite lists), since a reverted part carries the other flag; each step is counted under the kind it has in the child." The `AllostericMechanism` method is unchanged: it builds the child, which carries the reverted flags.

Next to `_apply_bipartitions` add:

```julia
"""
The bipartition `bp` of `group` with every part none of whose steps carries flux
(`flags`, one per step of `group`) rebuilt at rapid equilibrium; `bp` itself when
no part changes, and always for a rapid-equilibrium group. A steady-state group
with no flux-carrying step exposes only the ratio of its constants, and its
rapid-equilibrium form is the same family with one constant fewer, so the split
emits that form instead of the raw part.
"""
function _revert_zero_flux_parts(group::Vector{Step}, bp, flags::BitVector)
    is_equilibrium(first(group)) && return bp
    carries = Dict(s => flags[j] for (j, s) in enumerate(group))
    revert(part) = any(s -> carries[s], part) ? part :
        Step[Step(from_species(s), to_species(s), consumed(s), released(s), true)
             for s in part]
    p1, p2 = revert(bp[1]), revert(bp[2])
    p1 === bp[1] && p2 === bp[2] ? bp : (p1, p2)
end
```

- [ ] **Step 4: Run the testsets, then the whole enumeration file**

Scratch runs: all three PASS. Focused run of `test/test_mechanism_enumeration.jl`: every test passes. In particular "mixed RE/SS multi-step groups: every child gains" is unchanged (both steady-state A bindings of random-order bi-bi lie on the unbalanced block, so nothing reverts), the ter-ter split tests are unchanged, and "expand_mechanisms: ter-ter random-order seed within budget" still emits 81 children under 120 s. If an exact-children split test changes, its parent has a steady-state group with a pendant part: derive the flags by hand, replace the raw child by the reverted one in `expected`, and add the rank checks as in the revert test above.

- [ ] **Step 5: Run the full suite, then commit**

Start the full suite in the background with its output in a log, poll the log in a bounded foreground loop until it prints its summary, and commit only when every test passes:

```bash
git add src/thermodynamic_constr_for_rate_eq_derivation.jl src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Emit a zero-flux split part at rapid equilibrium

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R4zCpbZSoygD66kRecfrDN"
```

---

### Task 5: `expand_mechanisms` asserts both rules on each parent; `seed_mechanisms` errors on an empty seed set

**Files:**
- Modify: `src/mechanism_enumeration.jl` (new `_assert_emission_rules` next to `_assert_chemistry_is_iso`, about line 164; `expand_mechanisms`, about line 2481; `seed_mechanisms`, about lines 2570–2610, and its docstring)
- Test: `test/test_mechanism_enumeration.jl` (new testsets next to "expand_mechanisms rejects chemistry folded into a release step", about line 7330, and next to the `seed_mechanisms` tests, found with `grep -n 'seed_mechanisms' test/test_mechanism_enumeration.jl`)

**Interfaces:**
- Consumes: `_flux_carrying_groups(m)` (Task 1), `_duplicate_copy_groups(m)` (Task 3 as amended: takes the mechanism, reads each group's allosteric tag), `_forward_sides(s)` (types.jl).
- Produces: `_assert_emission_rules(m::Union{Mechanism, AllostericMechanism})::Nothing` (used by Task 6's population test).

- [ ] **Step 1: Write the failing tests**

Next to "expand_mechanisms rejects chemistry folded into a release step" add:

```julia
@testset "expand_mechanisms rejects a parent with a zero-flux steady-state group" begin
    # E + A and E + Q at steady state while E(Q) + A ⇌ E(A, Q) ⇌ E(A) + Q stays at
    # rapid equilibrium: {E} joins the rest by two edges of weight 0, so neither
    # group carries flux and their constants enter the rate only as ratios. The
    # flip never emits this child; a hand-written parent is refused.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    shunt = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(Q) + A ⇌ E(A, Q)
            E + Q <--> E(Q)
            E(A) + Q ⇌ E(A, Q)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
        end
    end)
    err = try
        EnzymeRates.expand_mechanisms([shunt], rxn); nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("carries net flux", err.msg)
    @test occursin("E + A", err.msg) || occursin("E + Q", err.msg)
end

@testset "expand_mechanisms rejects a parent whose inhibitor copy binds only at twin sites" begin
    # A bound as its own competitive inhibitor at E alone: E(A::Inh) has the
    # composition of E(A), so the copy's constant enters the rate only added to
    # K_A's. The dead-end move never emits this child; a hand-written parent is
    # refused.
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: A
    end
    twin_only = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
            E + A::Inh ⇌ E(A::Inh)
        end
    end)
    err = try
        EnzymeRates.expand_mechanisms([twin_only], rxn); nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("duplicates an existing form", err.msg)
    # The same copy placed where it also creates a new complex is a valid parent.
    kept = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            (E + Q ⇌ E(Q), E(A::Inh) + Q ⇌ E(A::Inh, Q))
            (E + A::Inh ⇌ E(A::Inh), E(Q) + A::Inh ⇌ E(A::Inh, Q))
        end
    end)
    @test !isempty(EnzymeRates.expand_mechanisms([kept], rxn))
    @test EnzymeRates._assert_emission_rules(kept) === nothing
end
```

Next to the existing `seed_mechanisms` tests add:

```julia
@testset "seed_mechanisms errors when no mechanism binds every required regulator" begin
    # Uni-uni with S as its own competitive inhibitor: the only site that carries
    # neither S nor P is free E, where the copy duplicates E(S). No seed exists,
    # and the beam must say so rather than return nothing.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: S
    end
    err = try
        EnzymeRates.seed_mechanisms(rxn, Set{Symbol}(), Set([:S])); nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("required regulator", err.msg)
    @test occursin("competitive inhibitors: S", err.msg)
    @test occursin("optional_competitive_inhibitors", err.msg)
    # Bi-bi with A as its own inhibitor has placements that create a new complex.
    bibi = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        dead_end_inhibitors: A
    end
    seeds = EnzymeRates.seed_mechanisms(bibi, Set{Symbol}(), Set([:A]))
    @test !isempty(seeds)
    @test all(m -> EnzymeRates._assert_emission_rules(m) === nothing, seeds)
    @test all(m -> :A in EnzymeRates._bound_comp_inhibitors(m), seeds)
end
```

- [ ] **Step 2: Run the testsets and watch them fail**

Scratch run. Expected: the two `expand_mechanisms` testsets fail (no error is raised; `err` is `nothing`) and the `seed_mechanisms` testset fails (an empty vector is returned).

- [ ] **Step 3: Implement**

In `src/mechanism_enumeration.jl`, after `_assert_chemistry_is_iso`, add:

```julia
"""
    _assert_emission_rules(m)

The two rules every mechanism the moves emit satisfies, checked on a parent before
it is expanded. Every steady-state kinetic group holds a step that carries net flux
(`_flux_carrying_groups`): otherwise its two constants enter the rate only as their
ratio. Every kinetic group that binds a competitive inhibitor binds it somewhere
that creates a new complex (`_duplicate_copy_groups`): otherwise its constant is not
separable from the existing binding's. The flip tests only the groups it flips and
the split only the parts it makes, so a parent must already satisfy both; a move
that emitted a violator fails here at the next expansion instead of propagating it.
"""
function _assert_emission_rules(m::Union{Mechanism, AllostericMechanism})
    label(g) = join((join(_forward_sides(s), " → ") for s in steps(m)[g]), ", ")
    flux = _flux_carrying_groups(m)
    for (g, group) in enumerate(steps(m))
        is_equilibrium(first(group)) || flux[g] || error(
            "expand_mechanisms: steady-state kinetic group {" * label(g) * "} has no " *
            "step that carries net flux, so its two constants enter the rate only as " *
            "their ratio; write the group at rapid equilibrium")
    end
    for g in _duplicate_copy_groups(m)
        error("expand_mechanisms: kinetic group {" * label(g) * "} binds a competitive " *
              "inhibitor only where the complex duplicates an existing form, so its " *
              "constant is not separable from the existing binding's; bind the " *
              "inhibitor where it forms a new complex, or drop it")
    end
    nothing
end
```

In `expand_mechanisms`, after `_assert_chemistry_is_iso(m)` add `_assert_emission_rules(m)`, and in its docstring add the sentence: "Each parent must satisfy the two emission rules (`_assert_emission_rules`): every steady-state group carries flux and every competitive-inhibitor group creates a new complex."

In `seed_mechanisms`, before the final `seeds`, add:

```julia
    if isempty(seeds)
        list(s) = isempty(s) ? "none" : join(sort!(collect(s)), ", ")
        error("seed_mechanisms: no mechanism binds every required regulator " *
              "(competitive inhibitors: " * list(required_comp) *
              "; allosteric regulators: " * list(required_allo) * "). A substrate or " *
              "product declared as a competitive inhibitor binds only where it forms " *
              "a new complex, and a uni-uni mechanism has no such site. Mark a " *
              "regulator optional with `optional_competitive_inhibitors` or " *
              "`optional_allosteric_regulators`, or remove its declaration")
    end
```

and end the docstring with: "Errors when no mechanism binds every required regulator, naming the regulators and the keywords that make one optional; the first case is a uni-uni reaction with its substrate or product declared as a competitive inhibitor, whose only placement duplicates an existing complex."

- [ ] **Step 4: Run the testsets, then the whole enumeration file and the identify file**

Scratch runs: PASS. Focused run of `test/test_mechanism_enumeration.jl` and of `test/test_identify_rate_equation.jl`. Every call of `expand_mechanisms` on a hand-written parent (22 in the enumeration file, 1 in the identify file) must still pass. If one now errors, the fixture violates a rule: rewrite it to satisfy the rule while keeping the test's intent (an all-RE version of a zero-flux group; a copy placed where it creates a new complex), say so in a comment, and list the change in the task report. Do not delete a test.

- [ ] **Step 5: Run the full suite, then commit**

Start the full suite in the background with its output in a log, poll the log in a bounded foreground loop until it prints its summary, and commit only when every test passes:

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Assert the emission rules on parents and refuse an empty seed set

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R4zCpbZSoygD66kRecfrDN"
```

---

### Task 6: Population pin, docs, version, regression record

> Amended 2026-10-01 (Denis): no version bump in B; the branch rises once over `main`, in A.

**Files:**
- Test: `test/test_mechanism_enumeration.jl` (new testset next to "expand_mechanisms: ter-ter random-order seed within budget", about line 7305)
- Modify: `docs/src/developer.md` (the paragraph beginning "The two refinement moves never emit", about lines 86–100)
- Modify: `docs/src/identify/enumeration_engine.md` (about lines 93–96, 123–130, 133–160, 242–246, 262–264, 294–297)
- Modify: `Project.toml` (`version = "0.8.0"` → `"0.9.0"`)
- Modify: `docs/superpowers/specs/2026-09-30-a-plus-b-plus-c-handoff.md` (the "Next" bullet)
- Create: `docs/superpowers/specs/2026-09-30-exact-filters-regression.md`
- Scratch (not committed): `<scratchpad>/regress/regress.jl`, a detached worktree of a2a02b1 at `<scratchpad>/regress/base`

**Interfaces:**
- Consumes: `_assert_emission_rules(m)` (Task 5), `_flux_carrying_groups(m)` (Task 1), `_duplicate_copy_groups(m)` (Task 3 as amended), `_sig_of`/`_mechanism_from_sig` (types.jl).

- [ ] **Step 1: Write the population test and watch it fail**

Next to "expand_mechanisms: ter-ter random-order seed within budget" add:

```julia
@testset "catalytic moves on bi-bi to depth 2: counts and both rules on every child" begin
    # Aggregate regression pin over the whole enumerated population of the
    # findings' reaction R4, whose atoms admit ping-pong: `init_mechanisms` plus two
    # levels of the flip, split and dead-end moves, deduplicated across levels.
    # Before the rules the levels held 62, 369 and 1,388 mechanisms. The rules drop
    # 148 zero-flux flip children at level 2 and turn the 40 zero-flux split
    # children into duplicates of level-2 flip children, so 1,200 remain; no seed
    # or level-1 child changes. Every mechanism satisfies both emission rules.
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    moves(m) = vcat(EnzymeRates._expand_re_to_ss(m),
                    EnzymeRates._expand_split_kinetic_group(m),
                    EnzymeRates._expand_add_dead_end_regulator(m, rxn))
    level = unique!(EnzymeRates.init_mechanisms(rxn))
    seen = Set(level)
    counts = [length(level)]
    for _ in 1:2
        next = EnzymeRates.Mechanism[]
        for m in level, c in moves(m)
            c in seen && continue
            push!(seen, c); push!(next, c)
        end
        level = next
        push!(counts, length(level))
    end
    @test counts == [62, 369, 1200]
    @test all(m -> EnzymeRates._assert_emission_rules(m) === nothing, seen)

    # The same seeds with every reactant also a competitive inhibitor (R6), one
    # level: every dead-end child binds its copy where it creates a new complex.
    rxn6 = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        dead_end_inhibitors: A, B, P, Q
    end
    seeds6 = unique!(EnzymeRates.init_mechanisms(rxn6))
    @test length(seeds6) == 62
    kids6 = unique!(vcat((vcat(EnzymeRates._expand_re_to_ss(m),
                                EnzymeRates._expand_split_kinetic_group(m),
                                EnzymeRates._expand_add_dead_end_regulator(m, rxn6))
                           for m in seeds6)...))
    @test all(m -> EnzymeRates._assert_emission_rules(m) === nothing, kids6)
    @test any(m -> !isempty(EnzymeRates._bound_comp_inhibitors(m)), kids6)
end
```

Run it in a scratch file. It must pass already if Tasks 1–5 are correct; the point of running it first is to confirm the counts. If `counts` differs from `[62, 369, 1200]`, stop: check the convention by running the same loop against the a2a02b1 worktree (Step 4 below), where it must give `[62, 369, 1388]`. A different baseline means the loop's convention differs from the findings' (fix the loop); the same baseline with a different post-change count means a rule is wrong (report, do not adjust the pin).

- [ ] **Step 2: Update the docs**

`docs/src/developer.md`: replace the paragraph beginning "The two refinement moves never emit a child" and ending "nothing in `src/` estimates identifiability numerically." with:

```markdown
The two refinement moves never emit a child that is provably a reparameterization
of its parent. The split move divides a group by binding context and accepts a
set of splits only if the independent-parameter count rises; the count comes from
the thermodynamic constraint solve, and for a `Mechanism` it is evaluated without
building the child, from a cycle basis computed once per parent
(`_partition_independent_count`). The RE→SS move flips whole groups and accepts a
set only if the rapid-equilibrium segment count rises and every flipped group
carries net flux in the child. Flux is decided on the graph of rapid-equilibrium
segments (`_flux_carrying_groups`): each steady-state step is an edge weighted by
its net uptake of substrates minus products, segment offsets included, and a step
carries flux exactly when its biconnected block holds a cycle of nonzero weight. A
steady-state group with no such step exposes only the ratio of its constants, so a
flip set that leaves one counts as failed and is extended, and a split part with
none is emitted at rapid equilibrium. A dead-end copy of a substrate or product
must create a complex no form duplicates, by composition or by segment and
offsets (`_twin_site_test`); the dead-end move skips a pattern whose sites are all
twins, and the split never isolates twin-only copy sites. `expand_mechanisms`
asserts both rules on every parent (`_assert_emission_rules`). Both refinement
moves share one minimal-set search (`_minimal_gaining_sets`). Duplicate equations
that survive these proofs are collapsed at compile time by `eq_hash`. A numerical
identifiability rank exists only in the test suite, as an oracle for the proofs;
nothing in `src/` estimates identifiability numerically.
```

`docs/src/identify/enumeration_engine.md`, six edits:

1. Replace the paragraph "Two kinds of group never flip. ... flipping it would add a parameter the data cannot determine." with:

```markdown
Two kinds of group never flip. Competitive-inhibitor binding stays at rapid
equilibrium by modeling choice. A steady-state group must carry net flux: a group
none of whose steps would carry flux with every step at steady state is never a
unit, and a set of flips that leaves one of its groups without flux in the child,
because an equilibrated route around it carries the turnover, counts as failed and
is extended like a set that divides no segment. Such a group's two constants enter
the equation only as their ratio, which its equilibrium form already has.
```

2. In "### 2. Split a kinetic group by binding context", after the paragraph ending "...before either constant is free." insert:

```markdown
A part of a steady-state group that carries no flux, such as an abortive binding
separated from the catalytic binding whose constants it shared, is emitted at
rapid equilibrium: its two rates would enter the equation only as their ratio, and
the equilibrium form is the same family with one constant fewer. The count test
runs on that form, so a constant the thermodynamic ties pull back is still
rejected. A competitive-inhibitor group is never divided so that one part binds
only where its complex duplicates an existing form (see move 3).
```

3. In "### 3. Add a competitive inhibitor binding site", change "subject to two rules:" to "subject to three rules:" and add after the "Mirror steps" bullet:

```markdown
- **A new complex.** A dead-end complex must differ from every form the mechanism
  has: in composition, or, when a rapid-equilibrium route joins it to an existing
  form, in the metabolites bound along that route. A substrate or product declared
  as its own competitive inhibitor otherwise binds in a second orientation of a
  complex the mechanism already has, and its constant enters the equation only
  added to the existing one. A placement whose every site duplicates a form is not
  emitted; in a uni-uni mechanism that is every placement, and the beam reports
  the unsatisfiable requirement.
```

4. Replace the paragraph "**Dead-end branches stay at equilibrium.** ... so the group never flips." with:

```markdown
**A steady-state group carries flux.** A group whose every step lies off every
cycle that runs the reaction carries no net flux at steady state. Its forward and
reverse rates enter the equation only as their ratio, which the equilibrium form
already has, so the moves never leave such a group at steady state: the flip does
not make one, and the split reverts the one it would create. The test is
structural, on the graph of rapid-equilibrium segments (`_flux_carrying_groups`),
and reads each step's metabolite lists, so it holds for fused and
Theorell–Chance steps as well.
```

5. Change the sentence ending "declaring a substrate as a dead-end inhibitor is how substrate inhibition enters an allosteric mechanism." to end "...is how substrate inhibition enters an allosteric mechanism, provided the inhibitor complex is new: a copy that only duplicates an existing complex is never emitted (move 3)."

6. In "**A child is never a provable copy of its parent.**" replace "a split the constraint solver ties back, a flip that leaves the segment count unchanged, a flip of a dead-end group." with "a split the constraint solver ties back, a flip that leaves the segment count unchanged, a flip that leaves a steady-state group without flux, a dead-end copy that creates no new complex."

Read each edited section once as a whole for flow and for stale sentences elsewhere on the page that still describe the former screen (search the page for "dead-end branch" and "chemistry step" and fix any that now contradict the above).

`Project.toml`: `version = "0.9.0"`.

`docs/superpowers/specs/2026-09-30-a-plus-b-plus-c-handoff.md`: replace the bullet beginning "- Next: B, exact filters" with "- B, exact filters, is implemented (spec `2026-09-30-exact-filters-design.md`, plan `../plans/2026-09-30-exact-filters.md`, regression record `2026-09-30-exact-filters-regression.md`). Next: C, merged and Theorell–Chance seeds and canonical merges (findings plan section C), starting with its spec."

Build the docs to check the pages render: `julia --project=docs -e 'using Pkg; Pkg.develop(path="."); using Documenter, EnzymeRates; include("docs/make.jl")'` — or, if `docs/make.jl` deploys, run only the doctests: `julia --project=docs -e 'using Documenter, EnzymeRates; doctest(EnzymeRates)'`. Either must finish without an error.

- [ ] **Step 3: Write the regression script**

Save as `<scratchpad>/regress/regress.jl` (`<scratchpad>` is the session scratch directory named in the environment; it is not committed):

```julia
# ABOUTME: Regression check for the exact enumeration filters: snapshots the enumerated
# ABOUTME: populations of one checkout and compares two snapshots mechanism by mechanism.
# Usage (one Julia process at a time):
#   julia --project=<checkout> regress.jl snap <outdir>
#   julia --project=<new checkout> regress.jl compare <old outdir> <new outdir> <report.md>
using EnzymeRates, LinearAlgebra, Random, Statistics
const ER = EnzymeRates

function reactions()
    r1 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    r2 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: I
    end
    r3 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: S, P
    end
    r4 = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    r5 = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        dead_end_inhibitors: I
    end
    r6 = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        dead_end_inhibitors: A, B, P, Q
    end
    ["R1" => r1, "R2" => r2, "R3" => r3, "R4" => r4, "R5" => r5, "R6" => r6]
end

key(m::ER.Mechanism) = repr(ER._sig_of(m))
key(am::ER.AllostericMechanism) = repr((
    ER._sig_of(ER.Mechanism(ER.reaction(am), ER.steps(am))),
    Tuple(ER.cat_allo_states(am)), ER.catalytic_multiplicity(am)))
function from_key(k::AbstractString, allosteric::Bool)
    sig = eval(Meta.parse(k))
    allosteric || return ER._mechanism_from_sig(sig)
    cm = ER._mechanism_from_sig(sig[1])
    ER.AllostericMechanism(ER.reaction(cm), ER.steps(cm), collect(sig[2]), sig[3],
                           ER.RegulatorySite[])
end

catalytic_moves(m, rxn) = vcat(ER._expand_re_to_ss(m), ER._expand_split_kinetic_group(m),
                               ER._expand_add_dead_end_regulator(m, rxn))

"Levels 0–2 of the catalytic moves on `rxn`, deduplicated across levels."
function population(rxn)
    level = unique!(ER.init_mechanisms(rxn))
    seen = Set(level); levels = [level]
    for _ in 1:2
        next = ER.Mechanism[]
        for m in level, c in catalytic_moves(m, rxn)
            c in seen && continue
            push!(seen, c); push!(next, c)
        end
        level = next; push!(levels, level)
    end
    levels
end

"Allosteric children of the R4 seeds (level 1) and their flip and split children."
function allosteric_population(rxn)
    seeds = unique!(ER.init_mechanisms(rxn))
    l1 = unique!(vcat((ER._expand_to_allosteric(m, rxn) for m in seeds)...))
    seen = Set(l1); l2 = ER.AllostericMechanism[]
    for am in l1, c in vcat(ER._expand_re_to_ss(am), ER._expand_split_kinetic_group(am))
        c in seen && continue
        push!(seen, c); push!(l2, c)
    end
    [l1, l2]
end

function derive(m)
    try
        em = ER.compile_mechanism(m)
        (join(string.(ER.fitted_params(em)), " "), ER.rate_equation_string(em))
    catch e      # a mechanism over the term budget errors the same way in both runs
        ("ERROR", sprint(showerror, e))
    end
end

"Finite-difference rank of ∂v/∂log θ over the fitted parameters (the test suite's oracle)."
function identifiable_rank(m; npts = 60, ndraws = 3, h = 1e-5)
    em = ER.compile_mechanism(m)
    fp = collect(ER.fitted_params(em))
    cm = m isa ER.Mechanism ? m : ER._state_mechanism(m, :A)
    mets = sort!(collect(ER._concentration_symbols(cm)))
    rng = MersenneTwister(hash(ER.rate_equation_string(em)) % 2^31)
    best = 0
    for _ in 1:ndraws
        θ = exp.(randn(rng, length(fp)))
        keq = exp(randn(rng))
        concs = [NamedTuple{Tuple(mets)}(Tuple(exp.(2 .* randn(rng, length(mets)))))
                 for _ in 1:npts]
        J = zeros(npts, length(fp))
        for j in eachindex(fp), sgn in (1, -1)
            θp = copy(θ); θp[j] *= exp(sgn * h)
            p = NamedTuple{(fp..., :Keq, :E_total)}((θp..., keq, 1.0))
            for (i, c) in enumerate(concs)
                J[i, j] += sgn * ER.rate_equation(em, c, p) / (2h)
            end
        end
        all(isfinite, J) || continue
        sv = svdvals(J)
        best = max(best, count(>(1e-7 * sv[1]), sv))
    end
    best
end

function snap(outdir)
    mkpath(outdir)
    keys_io = open(joinpath(outdir, "keys.tsv"), "w")
    strings_io = open(joinpath(outdir, "strings.tsv"), "w")
    timing_io = open(joinpath(outdir, "timing.tsv"), "w")
    rng = MersenneTwister(20260930)
    for (label, rxn) in reactions()
        t = @elapsed levels = population(rxn)
        println(timing_io, label, '\t', t)
        for (lvl, level) in enumerate(levels), m in level
            println(keys_io, label, '\t', lvl - 1, '\t', key(m))
        end
        all_m = vcat(levels...)
        sample = all_m[randperm(rng, length(all_m))[1:min(100, length(all_m))]]
        for m in sample
            names, str = derive(m)
            println(strings_io, label, '\t', key(m), '\t', names, '\t', str)
        end
        println(stderr, label, ": ", length.(levels), " in ", round(t; digits = 2), " s")
    end
    r4 = last(reactions()[4])
    t = @elapsed alevels = allosteric_population(r4)
    println(timing_io, "ALLO\t", t)
    for (lvl, level) in enumerate(alevels), am in level
        println(keys_io, "ALLO\t", lvl, '\t', key(am))
    end
    all_a = vcat(alevels...)
    for am in all_a[randperm(rng, length(all_a))[1:min(50, length(all_a))]]
        names, str = derive(am)
        println(strings_io, "ALLO\t", key(am), '\t', names, '\t', str)
    end
    println(stderr, "ALLO: ", length.(alevels), " in ", round(t; digits = 2), " s")
    close(keys_io); close(strings_io); close(timing_io)
end

read_keys(dir) = [(split(l, '\t'; limit = 3)...,) for l in eachline(joinpath(dir, "keys.tsv"))]
read_strings(dir) = Dict((r[1], r[2]) => (r[3], r[4]) for r in
    (split(l, '\t'; limit = 4) for l in eachline(joinpath(dir, "strings.tsv"))))
read_timing(dir) = Dict(split(l, '\t')[1] => parse(Float64, split(l, '\t')[2])
                        for l in eachline(joinpath(dir, "timing.tsv")))

zero_flux_groups(m) = [g for (g, grp) in enumerate(ER.steps(m))
                       if !ER.is_equilibrium(first(grp)) && !ER._flux_carrying_groups(m)[g]]
function re_twin(m::ER.Mechanism)
    zf = Set(zero_flux_groups(m))
    groups = [g in zf ? ER.Step[ER.Step(ER.from_species(s), ER.to_species(s), ER.consumed(s),
                                        ER.released(s), true) for s in grp] : grp
              for (g, grp) in enumerate(ER.steps(m))]
    ER.Mechanism(ER.reaction(m), groups)
end
obeys(m) = (ER._assert_emission_rules(m); true)

function compare(old, new, report)
    ko, kn = read_keys(old), read_keys(new)
    sets = unique(first.(ko))
    io = IOBuffer()
    println(io, "# Exact filters: regression check\n")
    println(io, "Populations: levels 0–2 of the flip, split and dead-end moves on R1–R6; ",
            "ALLO = allosteric children of the R4 seeds and their flip and split children.\n")
    println(io, "| Set | Level | Old | New | Kept | Removed | Added |\n|---|---|---|---|---|---|---|")
    removed_all = Tuple{String, String}[]; added_all = Tuple{String, String}[]
    for s in sets
        levels = sort!(unique([r[2] for r in ko if r[1] == s] ∪ [r[2] for r in kn if r[1] == s]))
        for l in levels
            o = Set(r[3] for r in ko if r[1] == s && r[2] == l)
            n = Set(r[3] for r in kn if r[1] == s && r[2] == l)
            println(io, "| $s | $l | $(length(o)) | $(length(n)) | $(length(o ∩ n)) | ",
                    "$(length(setdiff(o, n))) | $(length(setdiff(n, o))) |")
            append!(removed_all, (s, k) for k in setdiff(o, n))
            append!(added_all, (s, k) for k in setdiff(n, o))
        end
    end
    new_keys_by_set = Dict(s => Set(r[3] for r in kn if r[1] == s) for s in sets)
    println(io, "\n## Removed mechanisms\n")
    zf = 0; dup = 0; unexplained = String[]; twin_found = 0; twin_missing = String[]
    rank_ok = 0; rank_bad = String[]
    for (s, k) in removed_all
        m = from_key(k, s == "ALLO")
        z = zero_flux_groups(m)
        d = ER._duplicate_copy_groups(m)
        if !isempty(z)
            zf += 1
            if m isa ER.Mechanism
                key(re_twin(m)) in new_keys_by_set[s] ? (twin_found += 1) : push!(twin_missing, k)
                if s == "R4"
                    fitted = length(ER.fitted_params(ER.compile_mechanism(m)))
                    fitted - identifiable_rank(m) >= length(z) ? (rank_ok += 1) : push!(rank_bad, k)
                end
            end
        elseif !isempty(d)
            dup += 1
        else
            push!(unexplained, "$s $k")
        end
    end
    println(io, "- Removed: $(length(removed_all)); zero-flux group: $zf; twin-only copy group: $dup; ",
            "neither: $(length(unexplained)).")
    println(io, "- RE twin of a removed zero-flux mechanism present in the new population: ",
            "$twin_found of $zf; missing: $(length(twin_missing)).")
    println(io, "- Removed R4 zero-flux mechanisms with fitted − rank ≥ zero-flux groups: ",
            "$rank_ok; failures: $(length(rank_bad)).")
    for u in unexplained; println(io, "  - UNEXPLAINED: ", u); end
    for u in twin_missing; println(io, "  - TWIN MISSING: ", u); end
    for u in rank_bad; println(io, "  - RANK: ", u); end
    println(io, "\n## Added mechanisms\n")
    bad_added = [(s, k) for (s, k) in added_all if !obeys(from_key(k, s == "ALLO"))]
    println(io, "- Added: $(length(added_all)); violating a rule: $(length(bad_added)).")
    for (s, k) in added_all; println(io, "  - $s: $k"); end
    println(io, "\n## Derivation of kept mechanisms\n")
    so, sn = read_strings(old), read_strings(new)
    common = intersect(keys(so), keys(sn))
    diffs = [k for k in common if so[k] != sn[k]]
    println(io, "- Compared: $(length(common)); differing fitted names or Reduced string: ",
            "$(length(diffs)).")
    for k in diffs; println(io, "  - ", k[1], " ", k[2]); end
    println(io, "\n## Enumeration time (s)\n\n| Set | Old | New | Change |\n|---|---|---|---|")
    to, tn = read_timing(old), read_timing(new)
    for s in sort!(collect(keys(to)))
        println(io, "| $s | $(round(to[s]; digits = 2)) | $(round(tn[s]; digits = 2)) | ",
                "$(round(100 * (tn[s] / to[s] - 1); digits = 0))% |")
    end
    write(report, String(take!(io)))
    println(stderr, "report written to ", report)
end

if ARGS[1] == "snap"
    snap(ARGS[2])
elseif ARGS[1] == "compare"
    compare(ARGS[2], ARGS[3], ARGS[4])
end
```

- [ ] **Step 4: Run the regression**

```bash
S=<scratchpad>/regress
mkdir -p $S
git worktree add --detach $S/base a2a02b1
cd $S/base && julia --project=. -e 'using Pkg; Pkg.instantiate()' && cd -
nohup julia --project=$S/base $S/regress.jl snap $S/old > $S/snap_old.log 2>&1 &
# poll: for i in $(seq 1 60); do sleep 30; grep -q "ALLO:" $S/snap_old.log && break; done
nohup julia --project=. $S/regress.jl snap $S/new > $S/snap_new.log 2>&1 &
# poll the same way
julia --project=. $S/regress.jl compare $S/old $S/new $S/report.md
git worktree remove $S/base
```

Run the two snaps one after the other, never together. Read `$S/report.md`. Required: 0 UNEXPLAINED removed mechanisms, 0 TWIN MISSING, 0 RANK failures, 0 added mechanisms violating a rule, 0 derivation differences, and R4 levels 62 / 369 / 1,388 old against 62 / 369 / 1,200 new. Expected shape elsewhere: R1–R3 lose the R3 dead-end children (all uni-uni copies); R5 loses the same zero-flux structures as R4 and gains a few extended flip pairs only at level 3, which this population does not reach; R6 loses dead-end children at levels 1 and 2 and zero-flux flip children at level 2; ALLO loses flip children of allosteric parents whose flipped groups carry no flux, if any. The R6 enumeration time must be within 20% of the old run. Anything outside this shape is investigated before the record is written, not explained away.

- [ ] **Step 5: Write the regression record**

Create `docs/superpowers/specs/2026-09-30-exact-filters-regression.md` in the style of `2026-09-28-step-explicit-stoichiometry-regression.md`: a Method section (the populations, the key, the two commits, the script's checks, where the script lived), a Results table copied from the report, a section on the removed mechanisms (counts by rule, twin and rank results), one on added mechanisms (list them with the reason each is new: an extended flip pair or a reverted split child not present before), one on derivation of kept mechanisms, one on timing, and a Conclusion. Every number in the record comes from the report; state the commits the snaps ran on.

- [ ] **Step 6: Full suite, then commit**

Run the full suite in the background and poll it to the end. It must be green. Then:

```bash
git add test/test_mechanism_enumeration.jl docs/src/developer.md docs/src/identify/enumeration_engine.md Project.toml docs/superpowers/specs/2026-09-30-a-plus-b-plus-c-handoff.md docs/superpowers/specs/2026-09-30-exact-filters-regression.md
git commit -m "Pin the filtered bi-bi population; record the regression; bump to 0.9.0

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R4zCpbZSoygD66kRecfrDN"
```
