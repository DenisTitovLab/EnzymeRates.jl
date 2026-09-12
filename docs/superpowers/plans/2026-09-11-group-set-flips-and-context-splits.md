# Group-Set Flips and Context Splits Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the split move's canonicalization and the per-group RE→SS flip with two moves whose every rejection is a proof: whole-group flips in minimal segment-raising sets, and context-bipartition splits in minimal independent-count-raising sets.

**Architecture:** Both moves are "unit generator + shared minimal-set search + proof-based gain test." The flip move's units are eligible kinetic groups and its gain test is the RE segment count. The split move's units are `(group, context bipartition)` pairs and its gain test is the independent-parameter count, evaluated for `Mechanism` without constructing the child from a once-per-parent cycle basis. `_canonical_mechanism` and `_merge_tied_kinetic_groups` are deleted; compile-time equation dedup is the only dedup.

**Tech Stack:** Julia 1.12, EnzymeRates internals (`Mechanism`, `AllostericMechanism`, `Step`, `_compute_re_groups`, `_thermodynamic_constraints`, `_rref_partition`), `Test`.

**Spec:** `docs/superpowers/specs/2026-09-11-group-set-flips-and-context-splits-design.md`

## Global Constraints

- 92-character line length, 4-space indentation; match surrounding style.
- Every new source or test file starts with two `# ABOUTME:` lines.
- `rate_equation` runtime is untouched: 0 allocations, sub-120 ns (`test_rate_equation_performance`).
- No numerical identifiability test inside `src/`. The finite-difference rank lives only in tests.
- Names describe what code does, never history: no "new", "old", "fast", "legacy" in names or comments.
- Never remove a code comment unless it is provably false; comments in refactored code move with the code.
- Run focused tests via `TestEnv.activate()` + `include` (see *Running tests*). Run the full suite in the background (~13 min, exceeds the 600 s foreground cap) and never concurrently with another Julia process (7.7 GB RAM, no swap).
- Commit after every task with the attribution lines from the session's system reminder.
- `test/test_mechanism_enumeration.jl:4947-4956` pins `expand_mechanisms` to the six named move functions by name; the function names `_expand_re_to_ss` and `_expand_split_kinetic_group` do not change.

## Running tests

Focused (fast, use for TDD):

```bash
julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_flip_and_split_moves.jl")'
```

If `TestEnv` is not installed: `julia -e 'using Pkg; Pkg.add("TestEnv")'` once, in the default environment.

Full suite (background, detached, poll the log):

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates && (nohup julia --project -e 'using Pkg; Pkg.test()' > /tmp/enzymerates_fullsuite.log 2>&1 &)
# poll: tail -5 /tmp/enzymerates_fullsuite.log ; grep -E "Test Summary|Error|FAIL" /tmp/enzymerates_fullsuite.log
```

A `UndefVarError` for a helper defined in another test file (e.g. `uni_uni_rxn`) under a focused run is a focused-run artifact; define the helper locally in the new test file instead of depending on it.

## File structure

- `src/thermodynamic_constr_for_rate_eq_derivation.jl` — gains `_dependent_param_exprs(::Mechanism)` (the concrete form of the existing type method), `_independent_param_count`, and `_partition_independent_count` (the once-per-parent cycle-basis counter). Constraint math lives here already.
- `src/rate_eq_derivation.jl` — `_dependent_param_exprs(::Type{AllostericEnzymeMechanism})` becomes a delegate to a concrete `_dependent_param_exprs(::AllostericMechanism)`.
- `src/mechanism_enumeration.jl` — flux-carrying groups, RE segment count, `_minimal_gaining_sets`, context bipartitions, the two rewritten moves; loses `_re_to_ss_flip_units`, `_step_core`, `_split_one_step`, `_merge_tied_kinetic_groups` (both), `_canonical_mechanism` (both).
- `test/test_flip_and_split_moves.jl` — new file, every new test in this plan; included from `test/runtests.jl` after `test_mechanism_enumeration.jl`.
- `test/test_mechanism_enumeration.jl`, `test/test_types.jl`, `test/test_identify_rate_equation.jl` — existing tests updated or removed as listed per task.
- `docs/src/identify/enumeration_engine.md`, `docs/src/identify/combinatorics.md`, `docs/src/developer.md` — documentation.

---

### Task 1: Independent-parameter count on concrete mechanisms

**Files:**
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl:291-320` (`_dependent_param_exprs(M::Type{<:EnzymeMechanism})`)
- Modify: `src/rate_eq_derivation.jl:1502-1545` (`_dependent_param_exprs(::Type{AllostericEnzymeMechanism{CM,CS,RS}})`)
- Create: `test/test_flip_and_split_moves.jl`
- Modify: `test/runtests.jl:19` (add the include after `test_mechanism_enumeration.jl`)

**Interfaces:**
- Produces: `_dependent_param_exprs(m::Mechanism)`, `_dependent_param_exprs(am::AllostericMechanism)` — same `(dep_exprs, indep)` return as the type methods, which now delegate to them.
- Produces: `_independent_param_count(m::Union{Mechanism, AllostericMechanism})::Int = length(_dependent_param_exprs(m)[2])`.

- [ ] **Step 1: Create the test file and write the failing test**

```julia
# ABOUTME: Tests for the RE→SS group-set flip move and the context-bipartition split move:
# ABOUTME: eligibility, gain proofs, minimal sets, reachability, and the once-per-parent counter.
using Test
using EnzymeRates
using LinearAlgebra
using Random

const _bibi_rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end
const _uni_uni_rxn = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
end

@testset "_independent_param_count matches fitted_params on the test specs" begin
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism isa EnzymeRates.AllostericEnzymeMechanism ?
            EnzymeRates.AllostericMechanism(spec.mechanism) :
            EnzymeRates.Mechanism(spec.mechanism)
        @test EnzymeRates._independent_param_count(m) ==
              length(EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(m)))
    end
end
```

Add to `test/runtests.jl` after the `test_mechanism_enumeration.jl` line:

```julia
    include("test_flip_and_split_moves.jl")
```

- [ ] **Step 2: Run the focused test to verify it fails**

Run the focused command from *Running tests*.
Expected: FAIL with `UndefVarError: _independent_param_count not defined`.

- [ ] **Step 3: Factor the non-allosteric method onto the concrete mechanism**

In `src/thermodynamic_constr_for_rate_eq_derivation.jl`, replace the body of `_dependent_param_exprs(M::Type{<:EnzymeMechanism})` so the existing code and its comments live on the concrete method and the type method delegates:

```julia
function _dependent_param_exprs(mech::Mechanism)
    rename = _build_wegscheider_rename_map(mech)
    dep_exprs, indep = _dependent_param_exprs_kernel(mech, rename)
    # Filter Pass-2-absorbed symbols out of indep. Pass 2 of
    # `_build_wegscheider_rename_map` adds entries like `K_P_E => K_S_E`
    # when a Wegscheider tie collapses two binding-K group reps to the
    # same name. After the merge, the absorbed symbol doesn't appear in
    # the v polynomial — its column has been folded into the target.
    # But the absorbed symbol is still a kinetic-group rep in the
    # mechanism, so `_raw_param_symbols` emits it and the kernel keeps it
    # in `indep`. Without this filter, `fitted_params` exposes a fittable
    # dummy dimension that doesn't affect the loss, and finite-restart
    # convergence suffers (the same rate equation can land at noticeably
    # different fitted losses depending on which absorbed symbol got
    # the dummy slot).
    indep = Tuple(p for p in indep if get(rename, p, p) == p)
    # Sort by name so `fitted_params` / the params destructuring is
    # content-canonical: two mechanisms with the same independent set (e.g.
    # graph-distinct but rate-equivalent ones) produce the identical rate
    # equation string and therefore the same dedup key.
    indep = Tuple(sort(collect(indep); by = string))
    return dep_exprs, indep
end

_dependent_param_exprs(M::Type{<:EnzymeMechanism}) = _dependent_param_exprs(Mechanism(M()))

"""Number of independent (fitted) rate constants of a concrete mechanism, computed
from the thermodynamic constraint solve without compiling the mechanism."""
_independent_param_count(m::Union{Mechanism, AllostericMechanism}) =
    length(_dependent_param_exprs(m)[2])
```

Keep the existing docstring on the type method's position (move it above the concrete method; it describes the same computation).

- [ ] **Step 4: Factor the allosteric method the same way**

In `src/rate_eq_derivation.jl`, the method `_dependent_param_exprs(::Type{AllostericEnzymeMechanism{CM,CS,RS}})` begins with `am = AllostericMechanism(AllostericEnzymeMechanism{CM,CS,RS}())`. Change the signature to `function _dependent_param_exprs(am::AllostericMechanism)`, delete that first line, keep every other line and comment unchanged, and add after the function:

```julia
_dependent_param_exprs(::Type{AllostericEnzymeMechanism{CM,CS,RS}}) where {CM,CS,RS} =
    _dependent_param_exprs(AllostericMechanism(AllostericEnzymeMechanism{CM,CS,RS}()))
```

- [ ] **Step 5: Run the focused test to verify it passes**

Expected: PASS for every spec.

- [ ] **Step 6: Run the derivation and enumeration test files to check nothing else moved**

```bash
julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_rate_eq_derivation.jl")'
```
Expected: all pass.

- [ ] **Step 7: Commit**

```bash
git add src/thermodynamic_constr_for_rate_eq_derivation.jl src/rate_eq_derivation.jl test/test_flip_and_split_moves.jl test/runtests.jl
git commit -m "Independent-parameter count on concrete mechanisms"
```

---

### Task 2: Flux-carrying groups and RE segment count

**Files:**
- Modify: `src/mechanism_enumeration.jl` (add after `_flip_group_to_ss`, before `_expand_split_kinetic_group`)
- Test: `test/test_flip_and_split_moves.jl`

**Interfaces:**
- Produces: `_flux_carrying_groups(m::Union{Mechanism, AllostericMechanism})::BitVector` — one flag per kinetic group: the group holds a step that shares a biconnected block of the step graph with an isomerization step.
- Produces: `_re_segment_count(m::Union{Mechanism, AllostericMechanism})::Int` — number of rapid-equilibrium segments (A-state projection for allosteric).
- Consumes: `_compute_re_groups(::Mechanism)`, `_state_mechanism(am, :A)`.

- [ ] **Step 1: Write the failing tests**

Append to `test/test_flip_and_split_moves.jl`:

```julia
@testset "_flux_carrying_groups" begin
    # Ordered uni-uni: every step is on the catalytic cycle.
    m = first(EnzymeRates.init_mechanisms(_uni_uni_rxn))
    @test all(EnzymeRates._flux_carrying_groups(m))

    # A dead-end leaf hanging off the cycle is a bridge: not flux-carrying.
    leaf = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P ⇌ E(P)
            E(P) + S ⇌ E(P, S)
        end
    end)
    fc = EnzymeRates._flux_carrying_groups(leaf)
    # The leaf group is the one whose step forms the doubly-bound E(P, S).
    leaf_group = only(g for (g, grp) in enumerate(EnzymeRates.steps(leaf))
                      if any(s -> length(EnzymeRates.bound(EnzymeRates.to_species(s))) == 2, grp))
    @test !fc[leaf_group]
    @test count(!, fc) == 1

    # Ping-pong: the second chemistry step starts at rapid equilibrium and lies
    # on the catalytic cycle, so every group is flux-carrying.
    pingpong = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            Estar + B ⇌ Estar(B)
            E + Q ⇌ E(Q)
            Estar + P ⇌ Estar(A, P)
            E(A) <--> Estar(A, P)
            Estar(B) ⇌ E(Q)
        end
    end)
    @test all(EnzymeRates._flux_carrying_groups(pingpong))

    # A pendant binding-only square (inhibitor with a mirror) shares the E→E(S)
    # edge with the catalytic cycle, so its block contains chemistry: flux-carrying.
    mirror = EnzymeRates.AllostericMechanism(EnzymeRates.@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_inhibitors: I
        catalytic_steps: begin
            E + S ⇌ E(S)          :: EqualAI
            E(S) <--> E(P)        :: EqualAI
            E + P ⇌ E(P)          :: EqualAI
            E + I ⇌ E(I)          :: EqualAI
            E(I) + S ⇌ E(I, S)    :: EqualAI
        end
    end)
    @test all(EnzymeRates._flux_carrying_groups(mirror))
end

@testset "_re_segment_count" begin
    m = first(EnzymeRates.init_mechanisms(_uni_uni_rxn))
    @test EnzymeRates._re_segment_count(m) == 1
    groups = EnzymeRates.steps(m)
    g = findfirst(grp -> all(EnzymeRates.is_equilibrium, grp), groups)
    flipped = EnzymeRates.Mechanism(EnzymeRates.reaction(m),
                                    EnzymeRates._flip_group_to_ss(groups, g))
    @test EnzymeRates._re_segment_count(flipped) == 2
end
```

- [ ] **Step 2: Run the focused test to verify it fails**

Expected: FAIL with `UndefVarError: _flux_carrying_groups not defined`.

- [ ] **Step 3: Implement**

Add to `src/mechanism_enumeration.jl` immediately after `_flip_group_to_ss`:

```julia
"""
Biconnected blocks of an undirected multigraph. `edges[e] = (u, v)` with
vertices `1:nv`. Returns the block id of every edge (Tarjan's edge-stack
algorithm); a bridge is a block of its own.
"""
function _edge_blocks(nv::Int, edges::Vector{Tuple{Int, Int}})
    adj = [Int[] for _ in 1:nv]
    for (e, (u, v)) in enumerate(edges)
        push!(adj[u], e); push!(adj[v], e)
    end
    disc = zeros(Int, nv); low = zeros(Int, nv)
    block = zeros(Int, length(edges))
    stack = Int[]; clock = Ref(0); nblocks = Ref(0)
    function visit(u, parent_edge)
        clock[] += 1; disc[u] = low[u] = clock[]
        for e in adj[u]
            e == parent_edge && continue
            w = edges[e][1] == u ? edges[e][2] : edges[e][1]
            if disc[w] == 0
                push!(stack, e)
                visit(w, e)
                low[u] = min(low[u], low[w])
                if low[w] >= disc[u]
                    nblocks[] += 1
                    while true
                        x = pop!(stack); block[x] = nblocks[]
                        x == e && break
                    end
                end
            elseif disc[w] < disc[u]
                push!(stack, e)
                low[u] = min(low[u], disc[w])
            end
        end
    end
    for v in 1:nv
        disc[v] == 0 && visit(v, 0)
    end
    block
end

"""
    _flux_carrying_groups(m) -> BitVector

One flag per kinetic group: the group holds a step that lies on a cycle of the
step graph containing a chemistry (isomerization) step, i.e. shares a
biconnected block with one. A binding-only cycle satisfies detailed balance and
carries no net flux at steady state, so a group whose every step sits in such a
pendant region exposes only equilibrium ratios however it is flagged; flipping
it to steady state adds a phantom parameter. RE and SS steps are both edges
here: flux-carrying-ness depends on the graph, not on the flags.
"""
function _flux_carrying_groups(m::Union{Mechanism, AllostericMechanism})
    forms = Dict{Species, Int}()
    edges = Tuple{Int, Int}[]
    edge_group = Int[]
    edge_is_iso = Bool[]
    vertex(sp) = get!(forms, sp, length(forms) + 1)
    for (g, group) in enumerate(steps(m)), s in group
        push!(edges, (vertex(from_species(s)), vertex(to_species(s))))
        push!(edge_group, g); push!(edge_is_iso, is_iso(s))
    end
    block = _edge_blocks(length(forms), edges)
    chem_blocks = Set(block[e] for e in eachindex(edges) if edge_is_iso[e])
    flags = falses(length(steps(m)))
    for e in eachindex(edges)
        block[e] in chem_blocks && (flags[edge_group[e]] = true)
    end
    flags
end

"""Number of rapid-equilibrium segments (connected components of the RE
subgraph). An allosteric mechanism is measured on its A-state projection, which
holds every catalytic group."""
_re_segment_count(m::Mechanism) = length(_compute_re_groups(m)[2])
_re_segment_count(am::AllostericMechanism) = _re_segment_count(_state_mechanism(am, :A))
```

- [ ] **Step 4: Run the focused test to verify it passes**

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_flip_and_split_moves.jl
git commit -m "Flux-carrying groups via biconnected blocks; RE segment count"
```

---

### Task 3: Shared minimal-set search

**Files:**
- Modify: `src/mechanism_enumeration.jl` (add after `_re_segment_count`)
- Test: `test/test_flip_and_split_moves.jl`

**Interfaces:**
- Produces: `_minimal_gaining_sets(n::Int, gains, partners) -> Vector{Vector{Int}}`. `gains(set::Vector{Int})::Bool`; `partners(set::Vector{Int})` returns the unit indices allowed to extend `set` (an iterable of `Int`). Returns every minimal set (sorted index vectors, in the order found: by level, then lexicographic).

- [ ] **Step 1: Write the failing test**

```julia
@testset "_minimal_gaining_sets" begin
    # Units 1..4. Sets gain iff they contain {1,2} or contain 3.
    gains(set) = (1 in set && 2 in set) || 3 in set
    sets = EnzymeRates._minimal_gaining_sets(4, gains, _ -> 1:4)
    @test sets == [[3], [1, 2]]
    # Partner pruning: unit 2 may never join unit 1, so {1,2} is unreachable.
    partners(set) = (1 in set || 2 in set) ? [3, 4] : 1:4
    @test EnzymeRates._minimal_gaining_sets(4, gains, partners) == [[3]]
    # Nothing gains: empty result, and the search terminates.
    @test isempty(EnzymeRates._minimal_gaining_sets(3, _ -> false, _ -> 1:3))
    # No units at all.
    @test isempty(EnzymeRates._minimal_gaining_sets(0, _ -> true, _ -> 1:0))
end
```

- [ ] **Step 2: Run to verify it fails**

Expected: FAIL with `UndefVarError: _minimal_gaining_sets not defined`.

- [ ] **Step 3: Implement**

```julia
"""
    _minimal_gaining_sets(n, gains, partners) -> Vector{Vector{Int}}

Every minimal subset of units `1:n` for which `gains(set)` holds, found
Apriori-style: level 1 tests each unit; a level-`j` set is tested only if every
`(j−1)`-subset was tested and failed at the previous level, so no superset of a
gaining set is ever tested. `partners(set)` lists the units allowed to extend
`set` (units already in `set` are skipped). The loop ends when a level fails
nothing. Each returned set is sorted; sets are ordered by level, then
lexicographically.
"""
function _minimal_gaining_sets(n::Int, gains, partners)
    out = Vector{Int}[]
    failed = Vector{Int}[Int[]]
    while !isempty(failed)
        failed_keys = Set(failed)
        candidates = Set{Vector{Int}}()
        for set in failed, u in partners(set)
            u in set && continue
            c = sort!(vcat(set, u))
            c in candidates && continue
            all(setdiff(c, [x]) in failed_keys for x in c) || continue
            push!(candidates, c)
        end
        failed = Vector{Int}[]
        for c in sort!(collect(candidates))
            gains(c) ? push!(out, c) : push!(failed, c)
        end
    end
    out
end
```

- [ ] **Step 4: Run to verify it passes**

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_flip_and_split_moves.jl
git commit -m "Minimal gaining-set search shared by the flip and split moves"
```

---

### Task 4: The flip move

**Files:**
- Modify: `src/mechanism_enumeration.jl:1213-1274` (`_step_core`, `_re_to_ss_flip_units`, `_expand_re_to_ss`)
- Modify: `test/test_mechanism_enumeration.jl:1870-1905` (mirror-lock testset)
- Test: `test/test_flip_and_split_moves.jl`

**Interfaces:**
- Consumes: `_flux_carrying_groups`, `_re_segment_count`, `_minimal_gaining_sets`, `_flip_group_to_ss`, `_with_steps`.
- Produces: `_expand_re_to_ss(m)` with the same signature and return type as today. Deletes `_re_to_ss_flip_units` and `_step_core`.

- [ ] **Step 1: Write the failing tests**

```julia
"Whole-group flip of groups `gs` (test helper; production uses _flip_group_to_ss)."
function _flip_groups(m, gs)
    groups = EnzymeRates.steps(m)
    for g in gs
        groups = EnzymeRates._flip_group_to_ss(groups, g)
    end
    EnzymeRates._with_steps(m, groups)
end

@testset "_expand_re_to_ss: every child raises the RE segment count" begin
    for rxn in (_uni_uni_rxn, _bibi_rxn), m in EnzymeRates.init_mechanisms(rxn)
        kids = EnzymeRates._expand_re_to_ss(m)
        @test !isempty(kids)
        for c in kids
            @test EnzymeRates._re_segment_count(c) > EnzymeRates._re_segment_count(m)
        end
    end
end

@testset "_expand_re_to_ss: seed child count is unchanged (220 over bi-bi seeds)" begin
    seeds = EnzymeRates.init_mechanisms(_bibi_rxn)
    @test length(seeds) == 55
    @test sum(length(EnzymeRates._expand_re_to_ss(m)) for m in seeds) == 220
end

@testset "_expand_re_to_ss: a group with no flux-carrying step never flips" begin
    leaf = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P ⇌ E(P)
            E(P) + S ⇌ E(P, S)
        end
    end)
    fc = EnzymeRates._flux_carrying_groups(leaf)
    leaf_group = only(findall(!, fc))
    for c in EnzymeRates._expand_re_to_ss(leaf)
        @test all(EnzymeRates.is_equilibrium, EnzymeRates.steps(c)[leaf_group])
    end
    @test !isempty(EnzymeRates._expand_re_to_ss(leaf))
end

@testset "_expand_re_to_ss: two segment-flat groups flip together" begin
    # Random-order square with A's two binding steps in separate groups. Flipping
    # either A group alone leaves E and E(A) joined through the other A step, so
    # neither is emitted alone; the pair cuts the square and is emitted.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(B) + A ⇌ E(A, B)
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
            (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
            (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
            E(A, B) <--> E(P, Q)
        end
    end)
    groups = EnzymeRates.steps(m)
    a_groups = [g for (g, grp) in enumerate(groups)
                if length(grp) == 1 && EnzymeRates.name(EnzymeRates.bound_metabolite(only(grp))) == :A]
    @test length(a_groups) == 2
    base = EnzymeRates._re_segment_count(m)
    for g in a_groups
        @test EnzymeRates._re_segment_count(_flip_groups(m, [g])) == base
    end
    @test EnzymeRates._re_segment_count(_flip_groups(m, a_groups)) > base
    kids = EnzymeRates._expand_re_to_ss(m)
    # Children re-sort their groups, so match the A steps by content.
    a_steps_ss = [EnzymeRates.Step(EnzymeRates.from_species(s), EnzymeRates.to_species(s),
                                   EnzymeRates.bound_metabolite(s), false)
                  for g in a_groups for s in groups[g]]
    ss_steps(c) = Set(s for grp in EnzymeRates.steps(c) for s in grp if !EnzymeRates.is_equilibrium(s))
    n_a_ss(c) = count(s -> s in ss_steps(c), a_steps_ss)
    @test any(c -> n_a_ss(c) == 2, kids)
    @test !any(c -> n_a_ss(c) == 1, kids)
end

@testset "_expand_re_to_ss: uni-uni flips are emitted (documented no-op)" begin
    m = first(EnzymeRates.init_mechanisms(_uni_uni_rxn))
    @test length(EnzymeRates._expand_re_to_ss(m)) == 2
end

@testset "_expand_re_to_ss: no emitted set is a superset of another" begin
    for m in EnzymeRates.init_mechanisms(_bibi_rxn)[1:10]
        kids = EnzymeRates._expand_re_to_ss(m)
        flipped(c) = Set(g for (g, grp) in enumerate(EnzymeRates.steps(m))
                         if all(EnzymeRates.is_equilibrium, grp) &&
                            !any(EnzymeRates.is_equilibrium, EnzymeRates.steps(c)[g]))
        sets = flipped.(kids)
        for (i, s) in enumerate(sets), (j, t) in enumerate(sets)
            i != j && @test !(s ⊊ t)
        end
    end
end
```

Note: `_expand_re_to_ss` returns children built by `_with_steps`, whose constructor re-sorts groups, which is why the pair test matches the A steps by content rather than by group index. The superset test below uses the same step-content matching if group indices do not line up: compare the sets of parent groups whose steps are all SS in the child.

- [ ] **Step 2: Run to verify the new tests fail**

Expected: the "two segment-flat groups flip together" test fails (today's move emits each A group alone via `_re_to_ss_flip_units`, since the two groups share no inhibitor-free step core); the leaf test may pass or fail; the count test passes today (220 is today's number too).

- [ ] **Step 3: Rewrite the move**

Delete `_step_core` and `_re_to_ss_flip_units` (docstrings included) and replace `_expand_re_to_ss` with:

```julia
"""
    _expand_re_to_ss(m::Union{Mechanism, AllostericMechanism})

RE→SS expansion move. A flip unit is a whole kinetic group that is all-RE, binds
no regulator (competitive-inhibitor binding stays at rapid equilibrium by
modeling choice), and holds a flux-carrying step (`_flux_carrying_groups`; a
group with none exposes only equilibrium ratios and would gain a phantom
parameter). One child is produced per minimal set of units whose joint flip
raises the RE segment count (`_minimal_gaining_sets`): a single group when it
cuts a segment on its own, several groups when each alone is bridged by an RE
route through the others — as happens once a split has separated a
metabolite's binding steps, or a catalytic step from its inhibitor-bound
mirror. A flip that leaves the segment count unchanged adds an SS step whose
endpoints share a segment, which the rate equation never sees. All other
groups, the reaction, and (for allosteric) the catalytic-allo tags,
multiplicity, and regulatory sites are preserved verbatim.
"""
function _expand_re_to_ss(m::Union{Mechanism, AllostericMechanism})
    flux = _flux_carrying_groups(m)
    units = [g for g in kinetic_groups(m)
             if all(is_equilibrium, steps(m)[g]) && flux[g] &&
                !any(s -> bound_metabolite(s) isa Regulator, steps(m)[g])]
    child(sel) = begin
        groups = steps(m)
        for u in sel
            groups = _flip_group_to_ss(groups, units[u])
        end
        _with_steps(m, groups)
    end
    base = _re_segment_count(m)
    sets = _minimal_gaining_sets(length(units), sel -> _re_segment_count(child(sel)) > base,
                                 _ -> 1:length(units))
    typeof(m)[child(sel) for sel in sets]
end
```

- [ ] **Step 4: Update the mirror-lock testset in `test/test_mechanism_enumeration.jl:1870-1905`**

Replace the two lines that reference deleted helpers. The non-vacuity assertion `@test any(u -> length(u) > 1, EnzymeRates._re_to_ss_flip_units(m))` becomes:

```julia
        # Non-vacuity: the A-binding group and its inhibitor-bound mirror are
        # separate all-RE groups, each bridged by the other, so they can only
        # flip as a pair.
        kids = EnzymeRates._expand_re_to_ss(m)
        @test !isempty(kids)
```

and the mirror-lock loop replaces `EnzymeRates._step_core(s)` with a local inhibitor-free core:

```julia
        core(s) = begin
            strip(sp) = EnzymeRates.Species(
                EnzymeRates.Metabolite[b for b in EnzymeRates.bound(sp)
                                       if !(b isa EnzymeRates.Regulator)],
                EnzymeRates.conformation(sp), EnzymeRates.residual(sp))
            (strip(EnzymeRates.from_species(s)), strip(EnzymeRates.to_species(s)),
             EnzymeRates.bound_metabolite(s))
        end
```

using `kids` in place of the second `EnzymeRates._expand_re_to_ss(m)` call and `core(s)` in place of `_step_core(s)`. Update the testset's opening comment: the pair flips together because the minimal-set search finds that neither cuts alone, not because of a mirror class.

- [ ] **Step 5: Run the new file and the enumeration file**

```bash
julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_flip_and_split_moves.jl"); include("test/test_mechanism_enumeration.jl")'
```

Expected: all pass. The existing exact RE→SS counts (lines 1409, 1460, 1508, 1553, 1908 and the allosteric tag tests) are seed counts and must not move; if one does, stop and report the mechanism and both counts rather than editing the expectation.

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_flip_and_split_moves.jl test/test_mechanism_enumeration.jl
git commit -m "RE→SS move: whole-group flips in minimal segment-raising sets"
```

---

### Task 5: Context bipartitions and the split child builder

**Files:**
- Modify: `src/mechanism_enumeration.jl` (replace `_split_one_step`, lines 1488-1509)
- Test: `test/test_flip_and_split_moves.jl`

**Interfaces:**
- Produces: `_context_bipartitions(group::Vector{Step}) -> Vector{Tuple{Vector{Step}, Vector{Step}}}` — ordered by the context ligand; each part nonempty; the first part contains `first(group)`.
- Produces: `_apply_bipartitions(m, selection) -> typeof(m)` where `selection::Vector{Tuple{Int, Tuple{Vector{Step}, Vector{Step}}}}` pairs a group index with one of its bipartitions; each selected group is replaced by its two parts, allosteric tags duplicated.
- Deletes: `_split_one_step`.

- [ ] **Step 1: Write the failing tests**

```julia
@testset "_context_bipartitions" begin
    # A binds E, E(B), and E(P): contexts B and P give two bipartitions.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(P) + A ⇌ E(A, P))
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
            (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q), E(A) + P ⇌ E(A, P))
            (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
            E(A, B) <--> E(P, Q)
        end
    end)
    a_group = only(grp for grp in EnzymeRates.steps(m)
                   if length(grp) == 3 &&
                      EnzymeRates.name(EnzymeRates.bound_metabolite(first(grp))) == :A)
    bps = EnzymeRates._context_bipartitions(a_group)
    @test length(bps) == 2
    for (with, without) in bps
        @test !isempty(with) && !isempty(without)
        @test length(with) + length(without) == 3
        @test first(a_group) in with
        @test isempty(intersect(with, without))
    end
    # Contexts are ordered by ligand name: B before P.
    ctx(part) = Set(EnzymeRates.name(b) for s in part for b in EnzymeRates.bound(EnzymeRates.from_species(s)))
    @test :B in ctx(bps[1][2]) || :B in ctx(bps[1][1])
    # A two-step group with one context has one bipartition; a group whose
    # source forms carry no other ligand has none.
    b_group = only(grp for grp in EnzymeRates.steps(m)
                   if length(grp) == 2 &&
                      EnzymeRates.name(EnzymeRates.bound_metabolite(first(grp))) == :B)
    @test length(EnzymeRates._context_bipartitions(b_group)) == 1
    iso_group = only(grp for grp in EnzymeRates.steps(m) if EnzymeRates.is_iso(first(grp)))
    @test isempty(EnzymeRates._context_bipartitions(iso_group))
end

@testset "_context_bipartitions separates an inhibitor-bound mirror" begin
    am = EnzymeRates.AllostericMechanism(EnzymeRates.@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_inhibitors: I
        catalytic_steps: begin
            (E + S ⇌ E(S), E(I) + S ⇌ E(I, S))   :: EqualAI
            E(S) <--> E(P)                        :: EqualAI
            E + P ⇌ E(P)                          :: EqualAI
            E + I ⇌ E(I)                          :: EqualAI
        end
    end)
    s_group = only(grp for grp in EnzymeRates.steps(am) if length(grp) == 2)
    bps = EnzymeRates._context_bipartitions(s_group)
    @test length(bps) == 1
    @test all(part -> length(part) == 1, bps[1])
end

@testset "_apply_bipartitions" begin
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
            (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
            (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
            E(A, B) <--> E(P, Q)
        end
    end)
    groups = EnzymeRates.steps(m)
    g = findfirst(grp -> length(grp) == 2, groups)
    bp = only(EnzymeRates._context_bipartitions(groups[g]))
    child = EnzymeRates._apply_bipartitions(m, [(g, bp)])
    @test length(EnzymeRates.steps(child)) == length(groups) + 1
    @test EnzymeRates.n_steps(child) == EnzymeRates.n_steps(m)
    @test Set(s for grp in EnzymeRates.steps(child) for s in grp) ==
          Set(s for grp in groups for s in grp)
    @test any(grp -> Set(grp) == Set(bp[1]), EnzymeRates.steps(child))
    @test any(grp -> Set(grp) == Set(bp[2]), EnzymeRates.steps(child))

    am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        catalytic_multiplicity: 2
        catalytic_steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))    :: NonequalAI
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))    :: EqualAI
            E + P ⇌ E(P)             :: EqualAI
            E(P) + Q ⇌ E(P, Q)       :: EqualAI
            E + Q ⇌ E(Q)             :: EqualAI
            E(Q) + P ⇌ E(P, Q)       :: EqualAI
            E(A, B) <--> E(P, Q)     :: EqualAI
        end
    end)
    ga = findfirst(grp -> length(grp) == 2 &&
                   EnzymeRates.name(EnzymeRates.bound_metabolite(first(grp))) == :A,
                   EnzymeRates.steps(am))
    bpa = only(EnzymeRates._context_bipartitions(EnzymeRates.steps(am)[ga]))
    achild = EnzymeRates._apply_bipartitions(am, [(ga, bpa)])
    @test length(EnzymeRates.cat_allo_states(achild)) == length(EnzymeRates.steps(achild))
    for (gi, grp) in enumerate(EnzymeRates.steps(achild))
        Set(grp) ⊆ Set(bpa[1]) || Set(grp) ⊆ Set(bpa[2]) || continue
        @test EnzymeRates.cat_allo_state(achild, gi) == :NonequalAI
    end
    @test EnzymeRates.catalytic_multiplicity(achild) == 2
    @test EnzymeRates.regulatory_sites(achild) == EnzymeRates.regulatory_sites(am)
end
```

- [ ] **Step 2: Run to verify it fails**

Expected: FAIL with `UndefVarError: _context_bipartitions not defined`.

- [ ] **Step 3: Implement, replacing `_split_one_step`**

```julia
"""
The endpoint of `s` that does not carry the step's bound metabolite: the form the
metabolite binds to. Iso steps use `from_species`. Canonical RE binding steps
carry the metabolite on `to_species`; SS release steps may carry it on
`from_species`, so the test is on the bound list rather than on direction.
"""
function _context_form(s::Step)
    bm = bound_metabolite(s)
    bm === nothing && return from_species(s)
    any(b -> name(b) == name(bm), bound(from_species(s))) ? to_species(s) : from_species(s)
end

"""
    _context_bipartitions(group) -> Vector{Tuple{Vector{Step}, Vector{Step}}}

Ways to divide one kinetic group in two by binding context: for each other
ligand Y carried by some step's context form (`_context_form`), the steps whose
context form carries Y against the rest. Encodes "the affinity for this
metabolite may depend on which other ligand is already bound" — including a
competitive inhibitor, which separates a catalytic step from its inhibitor-bound
mirror. The group's own bound metabolite is never a context. Both parts are
nonempty, the first part holds the group's first step, ligands that induce the
same division give one bipartition, and the order follows the ligand's role and
name for deterministic output.
"""
function _context_bipartitions(group::Vector{Step})
    own = bound_metabolite(first(group))
    ligands = Set{Metabolite}()
    for s in group, b in bound(_context_form(s))
        own !== nothing && b == own && continue
        push!(ligands, b)
    end
    carries(s, y) = y in bound(_context_form(s))
    seen = Set{Vector{Step}}()
    out = Tuple{Vector{Step}, Vector{Step}}[]
    for y in sort!(collect(ligands); by = b -> (string(typeof(b)), string(name(b))))
        with = Step[s for s in group if carries(s, y)]
        without = Step[s for s in group if !carries(s, y)]
        (isempty(with) || isempty(without)) && continue
        first_part, second_part = first(group) in with ? (with, without) : (without, with)
        first_part in seen && continue
        push!(seen, first_part)
        push!(out, (first_part, second_part))
    end
    out
end

"""
Replace each selected group by its two bipartition parts. `selection` pairs a
group index with one of that group's `_context_bipartitions`; each group appears
at most once. For an allosteric mechanism both parts inherit the group's
catalytic allo-state tag (splitting is a parameter-relaxation move that must not
change A/I semantics).
"""
function _apply_bipartitions(m::Mechanism, selection)
    _with_steps(m, _bipartitioned_groups(steps(m), selection)[1])
end

function _apply_bipartitions(am::AllostericMechanism, selection)
    groups, origin = _bipartitioned_groups(steps(am), selection)
    _with_steps_and_cat_states(am, groups, cat_allo_states(am)[origin])
end

"""Groups of `groups` with each selected group replaced by its two parts, plus
the index of the original group each new group came from."""
function _bipartitioned_groups(groups::Vector{Vector{Step}}, selection)
    parts = Dict(g => bp for (g, bp) in selection)
    out = Vector{Vector{Step}}()
    origin = Int[]
    for (g, group) in enumerate(groups)
        if haskey(parts, g)
            push!(out, parts[g][1]); push!(origin, g)
            push!(out, parts[g][2]); push!(origin, g)
        else
            push!(out, group); push!(origin, g)
        end
    end
    out, origin
end
```

`Metabolite` values compare by `==` on type and fields (`Substrate(:A) == Substrate(:A)`, and `Product(:Lactate) != CompetitiveInhibitor(:Lactate)`); if `Metabolite` is not `==`-comparable, key the set on `(typeof(b), name(b))` instead and adapt `carries`.

- [ ] **Step 4: Run to verify it passes**

Expected: PASS. Tests in `test_identify_rate_equation.jl:1794-1796` still reference `_split_one_step`; they are removed in Task 7, so the focused run of the new file is the gate here.

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_flip_and_split_moves.jl
git commit -m "Context bipartitions and the split child builder"
```

---

### Task 6: Once-per-parent independent-count from a regrouping

**Files:**
- Modify: `src/thermodynamic_constr_for_rate_eq_derivation.jl` (add after `_independent_param_count`)
- Test: `test/test_flip_and_split_moves.jl`

**Interfaces:**
- Produces: `_partition_independent_count(parent::Mechanism) -> Function`; the returned `count(group_of_step::AbstractVector{Int})::Int` gives `_independent_param_count` of the mechanism obtained by regrouping `parent`'s flat steps (in `_flat_steps(parent)` order) into groups labelled by `group_of_step`, without constructing it.
- Consumes: `_thermodynamic_constraints`, `_rref_partition`, `_flat_steps`.

- [ ] **Step 1: Write the failing test**

```julia
@testset "_partition_independent_count agrees with _independent_param_count" begin
    seeds = EnzymeRates.init_mechanisms(_bibi_rxn)
    checked = 0
    for m in seeds
        count = EnzymeRates._partition_independent_count(m)
        flat = EnzymeRates._flat_steps(m)
        parent_ids = [g for (_, g) in flat]
        @test count(parent_ids) == EnzymeRates._independent_param_count(m)
        pos = Dict(s => j for (j, (s, _)) in enumerate(flat))
        groups = EnzymeRates.steps(m)
        for g in eachindex(groups), bp in EnzymeRates._context_bipartitions(groups[g])
            ids = copy(parent_ids)
            for s in bp[2]
                ids[pos[s]] = length(groups) + 1
            end
            child = EnzymeRates._apply_bipartitions(m, [(g, bp)])
            @test count(ids) == EnzymeRates._independent_param_count(child)
            checked += 1
        end
    end
    @test checked > 100
end
```

- [ ] **Step 2: Run to verify it fails**

Expected: FAIL with `UndefVarError: _partition_independent_count not defined`.

- [ ] **Step 3: Implement**

```julia
"""
    _partition_independent_count(parent::Mechanism) -> count

Return `count(group_of_step)`, the independent-parameter count of the mechanism
obtained by regrouping `parent`'s flat steps (in `_flat_steps` order) into the
groups labelled by `group_of_step`. Regrouping moves no edges, so the cycle basis
of the step graph (`_thermodynamic_constraints`) is the same for every
regrouping and is computed once here; each call only merges step columns by
group and takes the rank. Equals `_independent_param_count` of the constructed
child: the kernel's independent set is the columns minus the pivots, and folding
a single-symbol Wegscheider tie onto its target removes one column and one rank
together, so the count is invariant to the rename. Column sign conventions match
`_assemble_constraints`: a binding K enters with a sign flip, an iso K without,
an SS step contributes `+kf` and `-kr`.
"""
function _partition_independent_count(parent::Mechanism)
    C, _ = _thermodynamic_constraints(parent)
    kinds = [is_equilibrium(s) ? (is_binding(s) ? :binding_K : :iso_K) : :ss
             for (s, _) in _flat_steps(parent)]
    function count(group_of_step::AbstractVector{Int})
        length(group_of_step) == length(kinds) ||
            error("group_of_step must label every flat step of the parent")
        column = Dict{Tuple{Int, Int}, Int}()
        for (j, g) in enumerate(group_of_step)
            get!(column, (g, 1), length(column) + 1)
            kinds[j] === :ss && get!(column, (g, 2), length(column) + 1)
        end
        A = zeros(Int, size(C, 1), length(column))
        for (j, g) in enumerate(group_of_step), i in axes(C, 1)
            c = C[i, j]
            c == 0 && continue
            if kinds[j] === :ss
                A[i, column[(g, 1)]] += c
                A[i, column[(g, 2)]] -= c
            else
                A[i, column[(g, 1)]] += kinds[j] === :binding_K ? -c : c
            end
        end
        length(column) - length(_rref_partition(A)[1])
    end
    count
end
```

- [ ] **Step 4: Run to verify it passes**

Expected: PASS with `checked > 100` (the bi-bi seeds have 172 level-1 context candidates).

- [ ] **Step 5: Commit**

```bash
git add src/thermodynamic_constr_for_rate_eq_derivation.jl test/test_flip_and_split_moves.jl
git commit -m "Independent-parameter count of a regrouping from a once-per-parent cycle basis"
```

---

### Task 7: The split move; delete canonicalization

**Files:**
- Modify: `src/mechanism_enumeration.jl:1297-1346` (`_expand_split_kinetic_group`, both overloads) and delete lines 1348-1486 (`_merge_tied_kinetic_groups` ×2, `_canonical_mechanism` ×2)
- Modify: `test/test_mechanism_enumeration.jl:2018-2230` (split testsets), `:5375-5412` (delete "expand_mechanisms output is canonical")
- Modify: `test/test_types.jl:1199-1229` (delete the `_merge_tied_kinetic_groups` testset)
- Modify: `test/test_identify_rate_equation.jl:1630-1660`, `:1720-1735`, `:1737-1800`, `:1802-1840`
- Test: `test/test_flip_and_split_moves.jl`

**Interfaces:**
- Consumes: `_context_bipartitions`, `_apply_bipartitions`, `_partition_independent_count`, `_independent_param_count`, `_minimal_gaining_sets`, `_compute_re_groups`, `_state_mechanism`.
- Produces: `_expand_split_kinetic_group(m)` with the same signature and return type as today.

- [ ] **Step 1: Write the failing tests**

```julia
"Finite-difference rank of ∂v/∂log θ over the fitted parameters (test oracle only)."
function _identifiable_rank(m; npts = 60, ndraws = 3, h = 1e-5)
    em = EnzymeRates.compile_mechanism(m)
    fp = collect(EnzymeRates.fitted_params(em))
    cm = m isa EnzymeRates.Mechanism ? m : EnzymeRates._state_mechanism(m, :A)
    mets = sort!(collect(EnzymeRates._concentration_symbols(cm)))
    rng = MersenneTwister(hash(m) % 2^31)
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
                J[i, j] += sgn * rate_equation(em, c, p) / (2h)
            end
        end
        all(isfinite, J) || continue
        sv = svdvals(J)
        best = max(best, count(>(1e-7 * sv[1]), sv))
    end
    best
end

const _random_bibi = EnzymeRates.Mechanism(@enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(P) + A ⇌ E(A, P))
        (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
        (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q), E(A) + P ⇌ E(A, P))
        (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q), E(B) + Q ⇌ E(B, Q))
        E(A, B) <--> E(P, Q)
    end
end)

"Closure of `seed` under `gen`, by structural identity."
function _closure(seed, gen; maxn = 5_000)
    seen = Dict{UInt64, Any}(hash(seed) => seed); queue = Any[seed]
    while !isempty(queue)
        m = popfirst!(queue)
        for c in gen(m)
            h = hash(c)
            haskey(seen, h) && continue
            seen[h] = c; push!(queue, c)
            length(seen) > maxn && error("closure exceeded $maxn")
        end
    end
    collect(values(seen))
end

@testset "_expand_split_kinetic_group: random-order bi-bi reaches independent constants" begin
    # Today's canonicalization drops every split of this seed. The split closure
    # must reach the form with every step in its own group and 9 independent
    # parameters (measured: 16 structures).
    @test EnzymeRates._independent_param_count(_random_bibi) == 5
    cl = _closure(_random_bibi, EnzymeRates._expand_split_kinetic_group)
    @test maximum(length(EnzymeRates.steps(m)) for m in cl) == 13
    @test maximum(EnzymeRates._independent_param_count(m) for m in cl) == 9
    @test length(cl) == 16
end

@testset "_expand_split_kinetic_group: every child gains and no set contains another" begin
    for m in EnzymeRates.init_mechanisms(_bibi_rxn)
        base = EnzymeRates._independent_param_count(m)
        kids = EnzymeRates._expand_split_kinetic_group(m)
        for c in kids
            @test EnzymeRates._independent_param_count(c) > base
            @test EnzymeRates.n_steps(c) == EnzymeRates.n_steps(m)
        end
        # Minimality: the set of groups a child splits is never a strict superset
        # of another child's.
        split_groups(c) = Set(g for (g, grp) in enumerate(EnzymeRates.steps(m))
                              if !any(cg -> Set(cg) == Set(grp), EnzymeRates.steps(c)))
        sets = split_groups.(kids)
        for (i, s) in enumerate(sets), (j, t) in enumerate(sets)
            i != j && @test !(s ⊊ t)
        end
    end
end

@testset "_expand_split_kinetic_group: bi-bi seeds emit 102 children" begin
    seeds = EnzymeRates.init_mechanisms(_bibi_rxn)
    @test sum(length(EnzymeRates._expand_split_kinetic_group(m)) for m in seeds) == 102
end

@testset "_expand_split_kinetic_group: a rejected split is the parent's model" begin
    # Level-1 candidates the count test rejects have the parent's identifiable
    # rank; the rendered equation may differ only in which tied name survives.
    m = _random_bibi
    count = EnzymeRates._partition_independent_count(m)
    flat = EnzymeRates._flat_steps(m)
    pos = Dict(s => j for (j, (s, _)) in enumerate(flat))
    ids0 = [g for (_, g) in flat]
    base = count(ids0)
    r0 = _identifiable_rank(m)
    @test r0 == base
    rejected = 0
    for (g, grp) in enumerate(EnzymeRates.steps(m)), bp in EnzymeRates._context_bipartitions(grp)
        ids = copy(ids0)
        for s in bp[2]; ids[pos[s]] = length(EnzymeRates.steps(m)) + 1; end
        count(ids) > base && continue
        rejected += 1
        child = EnzymeRates._apply_bipartitions(m, [(g, bp)])
        @test _identifiable_rank(child) == r0
    end
    @test rejected > 0
end

@testset "_expand_split_kinetic_group: a four-step group splits two and two" begin
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(P) + A ⇌ E(A, P), E(B, P) + A ⇌ E(A, B, P))
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(P) + B ⇌ E(B, P), E(A, P) + B ⇌ E(A, B, P))
            (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q), E(A) + P ⇌ E(A, P), E(B) + P ⇌ E(B, P), E(A, B) + P ⇌ E(A, B, P))
            (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
            E(A, B) <--> E(P, Q)
        end
    end)
    a_group = only(grp for grp in EnzymeRates.steps(m)
                   if length(grp) == 4 && EnzymeRates.name(EnzymeRates.bound_metabolite(first(grp))) == :A)
    bps = EnzymeRates._context_bipartitions(a_group)
    # Context B (and context P) divide the four A steps two and two, which the
    # singleton carve of the old move could never produce.
    @test any(bp -> length(bp[1]) == 2 && length(bp[2]) == 2, bps)
    child = EnzymeRates._apply_bipartitions(m, [(findfirst(==(a_group), EnzymeRates.steps(m)),
                                               only(bp for bp in bps if length(bp[1]) == 2))])
    @test count(grp -> length(grp) == 2 && Set(grp) ⊆ Set(a_group), EnzymeRates.steps(child)) == 2
end

@testset "_expand_split_kinetic_group: allosteric parent" begin
    am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        catalytic_multiplicity: 2
        catalytic_steps: begin
            (E + A <--> E(A), E(B) + A <--> E(A, B))    :: NonequalAI
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))          :: EqualAI
            (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))          :: EqualAI
            (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))          :: EqualAI
            E(A, B) <--> E(P, Q)                        :: EqualAI
        end
    end)
    base = EnzymeRates._independent_param_count(am)
    kids = EnzymeRates._expand_split_kinetic_group(am)
    @test !isempty(kids)
    for c in kids
        @test c isa EnzymeRates.AllostericMechanism
        @test EnzymeRates._independent_param_count(c) > base
        @test length(EnzymeRates.cat_allo_states(c)) == length(EnzymeRates.steps(c))
    end
end

@testset "_expand_split_kinetic_group: ter-ter random-order seed finishes in budget" begin
    terter = @enzyme_reaction begin
        substrates: A[C], B[N], C[O]
        products: P[C], Q[N], R[O]
    end
    seeds = EnzymeRates.init_mechanisms(terter)
    nst = [EnzymeRates.n_steps(m) for m in seeds]
    worst = seeds[argmax(nst)]
    @test EnzymeRates.n_steps(worst) == 55
    t = @elapsed kids = EnzymeRates._expand_split_kinetic_group(worst)
    @test length(kids) == 12
    @test t < 60
end
```

- [ ] **Step 2: Run to verify the new tests fail**

Expected: the random-order reachability test fails (today's move returns nothing for `_random_bibi`); the 102 count fails.

- [ ] **Step 3: Rewrite the move and delete canonicalization**

Replace both `_expand_split_kinetic_group` methods and their docstring, and delete `_merge_tied_kinetic_groups` (both) and `_canonical_mechanism` (both) with their docstrings:

```julia
"""
    _expand_split_kinetic_group(m::Mechanism) → Vector{Mechanism}
    _expand_split_kinetic_group(am::AllostericMechanism) → Vector{AllostericMechanism}

Kinetic-group split move. A unit is one context bipartition of one group
(`_context_bipartitions`: the affinity for a metabolite may depend on which
other ligand is already bound). One child is produced per minimal set of units,
at most one per group, whose joint application raises the independent-parameter
count (`_minimal_gaining_sets`). A single bipartition whose new constant a
Wegscheider cycle ties straight back is not a model: the derivation substitutes
the constant away and the child's equation is the parent's. Such a bipartition
is kept as a seed for pairs and larger sets, because the cycle that absorbs it
is broken by also splitting another group on that cycle — random-order binding
needs one split per substrate before any constant frees up. Partners for a
rejected set are the all-RE groups sharing a rapid-equilibrium segment with its
carved steps: a tie needs an RE cycle through the carve, and every group on that
cycle lies in that segment. The count is evaluated without building the child
for a `Mechanism` (`_partition_independent_count`, cycle basis once per parent)
and by the combined state solve on the built child for an `AllostericMechanism`.
The reaction and (for allosteric) multiplicity and regulatory sites are
preserved; both parts of a split group inherit its catalytic allo-state tag.
"""
function _expand_split_kinetic_group(m::Union{Mechanism, AllostericMechanism})
    groups = steps(m)
    units = [(g, bp) for g in kinetic_groups(m) for bp in _context_bipartitions(groups[g])]
    isempty(units) && return typeof(m)[]
    selection(sel) = [units[u] for u in sel]
    gains = _split_gain_test(m, units)
    segments = _group_re_segments(m)
    partners(sel) = begin
        isempty(sel) && return 1:length(units)
        used = Set(units[u][1] for u in sel)
        touched = reduce(union, (segments[units[u][1]] for u in sel))
        [u for u in 1:length(units)
         if !(units[u][1] in used) && !isempty(intersect(segments[units[u][1]], touched))]
    end
    sets = _minimal_gaining_sets(length(units), gains, partners)
    typeof(m)[_apply_bipartitions(m, selection(sel)) for sel in sets]
end

"""Gain test for the split move: `sel -> Bool`, true when applying the selected
units raises the independent-parameter count above the parent's."""
function _split_gain_test(m::Mechanism, units)
    count = _partition_independent_count(m)
    flat = _flat_steps(m)
    position = Dict(s => j for (j, (s, _)) in enumerate(flat))
    parent_ids = [g for (_, g) in flat]
    base = count(parent_ids)
    function gains(sel)
        ids = copy(parent_ids)
        next_id = length(steps(m))
        for u in sel
            next_id += 1
            for s in units[u][2][2]
                ids[position[s]] = next_id
            end
        end
        count(ids) > base
    end
    gains
end

function _split_gain_test(am::AllostericMechanism, units)
    base = _independent_param_count(am)
    sel -> _independent_param_count(_apply_bipartitions(am, [units[u] for u in sel])) > base
end

"""RE segment ids touched by each kinetic group's steps; empty for a group holding
an SS step, which lies on no RE cycle and is never a split partner."""
function _group_re_segments(m::Union{Mechanism, AllostericMechanism})
    cm = m isa Mechanism ? m : _state_mechanism(m, :A)
    species, _, form_to_segment = _compute_re_groups(cm)
    segment(sp) = form_to_segment[findfirst(==(sp), species)]
    map(steps(m)) do group
        all(is_equilibrium, group) || return Set{Int}()
        Set{Int}(segment(from_species(s)) for s in group)
    end
end
```

`_state_mechanism(am, :A)` carries every group with the same `Species`, so `findfirst(==(sp), species)` resolves; if a species is missing (returns `nothing`), the projection differs from the assumption and the implementer must stop and report rather than guard it.

- [ ] **Step 4: Update the split testsets in `test/test_mechanism_enumeration.jl:2018-2230`**

Replace the four testsets' bodies as follows, keeping their seeds:

"Mechanism — mixed RE/SS multi-step groups split per member" → rename the testset to "Mechanism — mixed RE/SS multi-step groups: every child gains" and replace items 1, 2, and 4 with property assertions:

```julia
        result = EnzymeRates._expand_split_kinetic_group(m)
        @test !isempty(result)
        base = EnzymeRates._independent_param_count(m)
        for r in result
            @test EnzymeRates._independent_param_count(r) > base
            @test r isa EnzymeRates.Mechanism
            EnzymeRates._assert_mechanism_invariants(r)
            @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
            @test EnzymeRates.n_steps(r) == EnzymeRates.n_steps(m)
            @test EnzymeRates.reaction(r) == EnzymeRates.reaction(m)
        end
        # The SS A-binding group has no equilibrium constant to tie, so its
        # context split by B gains alone and is emitted as a single split.
        a_split = [r for r in result if length(EnzymeRates.steps(r)) == length(m.steps) + 1 &&
                   any(grp -> length(grp) == 1 && !EnzymeRates.is_equilibrium(only(grp)) &&
                       EnzymeRates.name(EnzymeRates.bound_metabolite(only(grp))) == :A,
                       EnzymeRates.steps(r))]
        @test length(a_split) == 1
```

"AllostericMechanism — SS multi-step :NonequalAI split" → keep items 3, 4, 5; replace items 1 and 2 with `@test !isempty(result)` and `@test all(r -> EnzymeRates._independent_param_count(r) > EnzymeRates._independent_param_count(am), result)`, and assert a child that splits only the `:NonequalAI` A group exists with `length(EnzymeRates.steps(r)) == length(am.cat_steps) + 1`.

"Mechanism — bi-bi random: all splits Wegscheider-tied (negative)" → rename to "Mechanism — bi-bi random: single splits are tied, the pair of A and B frees a constant" and replace the body's assertion:

```julia
        m = EnzymeRates.Mechanism(m_seed)
        result = EnzymeRates._expand_split_kinetic_group(m)
        @test !isempty(result)
        base = EnzymeRates._independent_param_count(m)
        for r in result
            @test EnzymeRates._independent_param_count(r) > base
            # Every child splits at least two groups: no single split gains here.
            @test length(EnzymeRates.steps(r)) >= length(m.steps) + 2
        end
```

"Mechanism — all singleton groups: empty (negative)" → unchanged.

"AllostericMechanism — bi-bi RE groups all Wegscheider-tied (negative)" → rename to "AllostericMechanism — bi-bi RE groups: every emitted split gains" and replace `@test isempty(...)` with `result = EnzymeRates._expand_split_kinetic_group(am)`, `@test !isempty(result)`, and `@test all(r -> EnzymeRates._independent_param_count(r) > EnzymeRates._independent_param_count(am), result)`. Do not assert how many groups a child splits: whether a single per-state split gains here is not measured.

Delete the testset "expand_mechanisms output is canonical" (`:5375-5412`) entirely.

- [ ] **Step 5: Update `test/test_types.jl:1199-1229`**

Delete the testset "`_merge_tied_kinetic_groups` re-merges tied splits".

- [ ] **Step 6: Update `test/test_identify_rate_equation.jl`**

- Testset "canonical-partition dedup collapses renaming-dups in _process_batch" (`:1630-1650`): keep the two precondition tests (`m1 != m2`, keys differ) and replace the three `_canonical_mechanism` lines with the model-level statement: `@test EnzymeRates._independent_param_count(m1) == EnzymeRates._independent_param_count(m2)`. Rename the testset to "renaming-dup pair: same independent count, different eq_hash" and rewrite its comment: the pair is the same model under the constraint solve; `eq_hash` does not see through the surviving-name choice, which is why the split move rejects by count rather than by equation string.
- Testset "allosteric canonical-partition dedup" (`:1720-1735`): same treatment (keep `am1 != am2`, keys differ; replace the `c1`/`c2` lines with `@test EnzymeRates._independent_param_count(am1) == EnzymeRates._independent_param_count(am2)`), rename to "allosteric renaming-dup pair: same independent count".
- Delete the testsets "_canonical_mechanism is idempotent" and "_canonical_mechanism errors when it does not converge" and the `_canon_a_parent` fixture with its comment block if nothing else uses it (grep first).
- Testset "_process_batch failures report the ORIGINAL mechanism": delete the line `@test_throws ErrorException EnzymeRates._canonical_mechanism(m_bad)` and the line `@test split != EnzymeRates._canonical_mechanism(split)   # merge changes it`; keep the rest. Update the PASS-1/PASS-2 comments to drop the canonicalization references (they describe a stage that no longer exists).

- [ ] **Step 7: Run the new file and the three updated files**

```bash
julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_flip_and_split_moves.jl"); include("test/test_mechanism_enumeration.jl"); include("test/test_types.jl")'
```

Then `test/test_identify_rate_equation.jl` the same way (it is slow; run it alone). Expected: all pass. If the ter-ter budget test exceeds 60 s, report the time; do not raise the budget. If an exact count (102, 16, 12) differs, report the actual number and the mechanism, do not edit the expectation.

- [ ] **Step 8: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_flip_and_split_moves.jl test/test_mechanism_enumeration.jl test/test_types.jl test/test_identify_rate_equation.jl
git commit -m "Split move: context bipartitions in minimal gaining sets; delete canonicalization"
```

---

### Task 8: Full suite, dead code, and the no-op flip measurement

**Files:**
- Modify: `src/mechanism_enumeration.jl` (dead-code pass)
- Test: full suite

- [ ] **Step 1: Grep for dead references**

```bash
grep -rn "_canonical_mechanism\|_merge_tied_kinetic_groups\|_split_one_step\|_re_to_ss_flip_units\|_step_core" src test docs
```
Expected: no hits in `src/` or `test/`. Hits in `docs/superpowers/specs` are history and stay.

- [ ] **Step 2: Re-read the changed source regions for further simplification**

Open `src/mechanism_enumeration.jl` around the two moves and `src/thermodynamic_constr_for_rate_eq_derivation.jl` around the counters. Remove any helper left with no caller. Do not restructure beyond that.

- [ ] **Step 3: Run the full suite in the background and wait for it**

Use the full-suite command from *Running tests*. Poll until `Test Summary` appears. Expected: 0 failures, 0 errors. `test_compile_budget.jl` and `test/fixtures/phase2_init_golden.txt` are untouched by this plan and must pass unchanged.

- [ ] **Step 4: Commit any dead-code removal**

```bash
git add src && git commit -m "Remove helpers left without callers by the move rewrite"
```

---

### Task 9: Documentation

**Files:**
- Modify: `docs/src/identify/enumeration_engine.md:63-96` and add a section after move 7
- Modify: `docs/src/identify/combinatorics.md` (full rewrite)
- Modify: `docs/src/developer.md:48-68`

Apply the `elements-of-style:writing-clearly-and-concisely` skill to every paragraph.

- [ ] **Step 1: Rewrite moves 1 and 2 in `enumeration_engine.md`**

Replace the two subsections with:

```markdown
### 1. Flip rapid-equilibrium groups to steady state

Turns whole kinetic groups from rapid equilibrium (RE) to steady state (SS): every
step in a group converts together, keeping one pair of rate constants per group.
The move emits one child per smallest set of groups whose joint flip divides a
rapid-equilibrium segment in two. On a starting mechanism every group cuts on its
own, so each RE group gives one child. Once splits have separated a metabolite's
binding steps, a single group may be bridged by an RE route through the others;
the move then flips the bridging groups together, because a steady-state step
whose two ends stay in one equilibrated segment never reaches the rate equation.

Two kinds of group never flip. Competitive-inhibitor binding stays at rapid
equilibrium by modeling choice. A group whose steps all lie on dead-end branches
carries no net flux at steady state, so the equation could only ever see its
equilibrium ratio; flipping it would add a parameter the data cannot determine.

**Parameter delta:** +1 per flipped group in most cases. A few flips leave the
equation unchanged even though they divide a segment; uni-uni is the classic
case, where the steady-state and rapid-equilibrium laws have the same form.
Those children are fit once and lose to their parent on parsimony.

### 2. Split a kinetic group by binding context

Divides one kinetic group into two, so the two parts have separate rate
constants. The division is always by context: the steps whose enzyme form already
carries some other ligand Y, against the rest. It encodes the hypothesis that the
affinity for a metabolite depends on what is bound next to it; with Y a
competitive inhibitor it separates a catalytic step from its inhibitor-bound
mirror.

A single split often frees no parameter. In random-order binding, thermodynamics
ties the new constant straight back to the old one, and the equation is
unchanged. The move therefore emits the smallest sets of simultaneous splits, at
most one per group, whose combined effect raises the number of independent
parameters; the thermodynamic constraint solver decides. Random-order bi-bi needs
the A group and the B group split together before either constant is free.

**Parameter delta:** at least +1 by construction. A child that would add nothing
is never emitted.
```

- [ ] **Step 2: Add a "Modeling choices" section after move 7**

```markdown
## Modeling choices

The moves encode a few decisions about which mechanisms are worth fitting. They
are choices, not theorems, and this section states them so they can be revisited.

**Steady-state detail enters by whole groups.** Steps that share a rate constant
flip to steady state together, in the smallest set of groups that separates a
new equilibrated segment. Sharing structure is refined first, by splits, and
steady-state detail follows it. A mechanism in which one enzyme form binds a
metabolite at equilibrium while another binds it at steady state is reachable,
but only after a split has given the two bindings separate constants.

**Dead-end branches stay at equilibrium.** A group whose every step lies off the
catalytic cycle carries no net flux at steady state. Its forward and reverse rates
enter the equation only as their ratio, which the equilibrium form already has,
so the group never flips.

**Competitive-inhibitor binding stays at equilibrium.** An inhibitor bound to two
forms that a catalytic step connects does carry flux through its mirror step, so
this is a modeling choice ("inhibitor binding is fast"), not a consequence of the
previous one.

**Splits follow binding context.** A group is divided only by whether another
ligand is already bound. Arbitrary partitions are not tried, a group is never
split three ways in one move, and groups never merge. A partition that no
sequence of context splits produces is unreachable.

**A child is never a provable copy of its parent.** Both moves reject a child
whose equation can be shown to equal the parent's: a split the constraint solver
ties back, a flip that leaves the segment count unchanged, a flip of a dead-end
group. A few equation-identical children survive, uni-uni flips among them; they
cost one fit each and never win selection.
```

- [ ] **Step 3: Rewrite `combinatorics.md`**

Replace the whole file with:

```markdown
# How many mechanisms? The combinatorics of enumeration

Even a simple enzyme has a surprising number of plausible mechanisms. The reason
is that a mechanism is built from several independent choices, and the counts
multiply. This page walks through those choices for an enzyme with two
substrates and two products, a bi-bi reaction, with concrete numbers. The point
is not the exact totals but how fast they grow, which is why
[`identify_rate_equation`](@ref) filters the search instead of fitting
everything.

## Binding order

Substrates can bind one after another in a fixed order, or in any order. Two
substrates give three arrangements: A first, B first, or either. Products leave
the same way, independently, so a bi-bi has 3 × 3 = 9 binding-and-release
skeletons. A third substrate raises the substrate side from 3 to 13, and a
ping-pong enzyme, which releases a product before the next substrate binds, adds
more skeletons on top.

## Dead-end complexes

An enzyme can bind a substrate and a product at the same time and get stuck: a
dead-end complex that sits off the catalytic cycle. Each pairing of a substrate
with a product may or may not be allowed, and every reactant must take part in
at least one. For a bi-bi that gives 7 patterns of dead ends. Together with
binding order, that is why `init_mechanisms` returns 55 bi-bi starting
mechanisms before any refinement, and 35,665 for a ter-ter reaction.

The `shared_catalytic_site` keyword removes patterns a chemist already knows are
impossible, such as ATP and ADP occupying one site together.

## Which forms share an affinity

A starting mechanism gives each metabolite one binding constant, whatever the
enzyme already holds. The search then asks whether the affinity for A changes
once B is bound, or once a product is bound. Each such question divides A's
binding steps into two groups with separate constants, and the questions
compose: A's affinity may depend on B alone, on P alone, or on both. For the
random-order bi-bi with dead ends, the search reaches 16 different ways of
sharing constants, from one constant per metabolite up to one constant per
binding step.

## Equilibrium or steady state

Each binding can be fast enough to sit at equilibrium, or slow enough that its
on and off rates both matter. The search flips groups of shared steps to steady
state, one group or several at a time, and each flip divides the enzyme's forms
into more regions that are in equilibrium among themselves. The number of ways
to divide the forms is the real count here: the random-order bi-bi with dead
ends has 9 enzyme forms and 1,434 distinct ways of partitioning them into
equilibrated regions. The number of steady-state variants grows with the number
of enzyme forms, which grows with substrates, products, and dead ends together.

## Competitive inhibitors

A declared competitive inhibitor can bind any set of the enzyme forms that hold
the metabolite it competes with. One inhibitor adds a median of about four
variants per mechanism; a second inhibitor binds independently and roughly
doubles that. Each further inhibitor multiplies the count again.

## Allosteric regulation

Making an enzyme allosteric adds a second conformation. Every binding group may
then bind in both conformations equally, in both with different affinities, or
in the active one only, and every regulator takes one of four states. A bi-bi
promotes to about 15 allosteric variants with no regulator declared, and each
regulator adds variants of its own; two regulators may also share a site or sit
at separate sites.

## Putting it together

The choices multiply. Binding order times dead ends gives the 55 starting
mechanisms; sharing structure, steady-state detail, inhibitors, and allosteric
states each multiply that by tens to thousands. Fitting one equation takes
seconds to minutes, so an exhaustive search of a bi-bi with one regulator would
take years on one core. The search therefore keeps only the promising
candidates at each parameter count, drops equations too dense to fit, starts
from fully regulated mechanisms when the regulators are known, and honors
`shared_catalytic_site`. See [Best mechanism selection](@ref) for those filters
and [The enumeration engine](@ref) for the moves that build the space.
```

- [ ] **Step 4: Update `developer.md`**

In *Enumeration engine architecture*, after the paragraph ending "…cost no compilation at all. Only the candidates the search actually fits are lifted to singleton types, one at a time, through `compile_mechanism`.", add:

```markdown
The two refinement moves never emit a child that is provably a reparameterization
of its parent. The split move divides a group by binding context and accepts a
set of splits only if the independent-parameter count rises; the count comes from
the thermodynamic constraint solve, and for a `Mechanism` it is evaluated without
building the child, from a cycle basis computed once per parent
(`_partition_independent_count`). The RE→SS move flips whole groups and accepts a
set only if the rapid-equilibrium segment count rises; groups with no
flux-carrying step (no cycle through a chemistry step, `_flux_carrying_groups`)
never flip. Both moves share one minimal-set search (`_minimal_gaining_sets`).
Duplicate equations that survive these proofs are collapsed at compile time by
`eq_hash`. A numerical identifiability rank exists only in the test suite, as an
oracle for the proofs; nothing in `src/` estimates identifiability numerically.
```

- [ ] **Step 5: Build the docs**

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate(); include("docs/make.jl")'
```
Expected: builds without doctest or cross-reference errors. If the docs environment is not set up, run the same with `Pkg.instantiate()` first and report any failure verbatim.

- [ ] **Step 6: Commit**

```bash
git add docs/src/identify/enumeration_engine.md docs/src/identify/combinatorics.md docs/src/developer.md
git commit -m "Docs: group-set flips, context splits, modeling choices, plain-language combinatorics"
```

---

### Task 10: Bump the version and record the change

**Files:**
- Modify: `Project.toml` (version)

- [ ] **Step 1: Bump the minor version**

The move semantics change what the search reaches, which is a behavior change for every user; bump `version` in `Project.toml` from `0.4.x` to `0.5.0`.

- [ ] **Step 2: Commit**

```bash
git add Project.toml && git commit -m "Bump version to 0.5.0"
```
