# Hyperbolic Catalysis in Conformational Mechanisms Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the enumerator refuse to combine an MWC conformational equilibrium with a catalytic scheme whose own rate equation carries concentration powers, while leaving hand-written mechanisms and dead-end inhibition untouched.

**Architecture:** One structural predicate, `_hyperbolic_catalysis(m)`, decides whether the flux-carrying step graph's rate equation has degree at most 1 in every substrate and product, from the rapid-equilibrium segment graph alone (no derivation). One trait, `_requires_hyperbolic_catalysis(::T)`, names the mechanism types the rule applies to. Two moves consult them: promotion to allosteric refuses a non-hyperbolic parent, and the RE→SS flip drops non-hyperbolic children of a conformational parent. A shared helper `_re_segment_extras` (extracted from `_bottomless_re_segment`) supplies the per-form metabolite exponents both need.

**Tech Stack:** Julia 1.x, EnzymeRates.jl internals (`Mechanism`, `AllostericMechanism`, `Step`, `Species`, `_flux_carrying_groups`, `_step_sides`, `_raw_symbolic_rate_polys`), Test.jl.

**Spec:** `docs/superpowers/specs/2026-09-22-hyperbolic-catalysis-in-conformational-mechanisms-design.md`

## Global Constraints

- 92-character line length, 4-space indentation; match surrounding style exactly.
- Every fixture in `test/test_mechanism_enumeration.jl` is written inline with `@enzyme_mechanism` / `@allosteric_mechanism`; move tests assert exact child counts and `Set(children) == Set(expected)`; file-level helpers are prefixed `_testhelper_`.
- The rule lives in enumeration moves only, never in a constructor.
- Dead-end groups (substrate, product, regulator) are ignored by the predicate.
- No new comments that describe history or a change; docstrings describe the code as it is.
- Focused test run (per file): `julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_mechanism_enumeration.jl")'`. Run ONE Julia process at a time (7.7 GB RAM, no swap). The enumeration file takes ~4 min warm.
- Commit messages end with the two attribution lines given in the session (Co-Authored-By and Claude-Session).

---

### Task 1: `_re_segment_extras` helper and the `_hyperbolic_catalysis` predicate

**Files:**
- Modify: `src/types.jl:585-645` (`_bottomless_re_segment` → extract `_re_segment_extras`)
- Modify: `src/mechanism_enumeration.jl` (add predicate + trait after `_flux_carrying_groups`, which ends near line 1363)
- Test: `test/test_mechanism_enumeration.jl` (new testset after the `_flux_carrying_groups` testset, which starts near line 5310)

**Interfaces:**
- Consumes: `_flux_carrying_groups(m) -> BitVector`, `_step_sides(s::Step) -> (e_lhs, e_rhs, m_lhs::Vector{Symbol}, m_rhs::Vector{Symbol})`, `steps(m)`, `kinetic_groups(m)`, `reaction(m)`, `substrates(rxn)`, `products(rxn)`, `name(::Metabolite)`.
- Produces:
  - `_re_segment_extras(steps::Vector{Vector{Step}}) -> (species::Vector{Species}, segments::Vector{Vector{Int}}, extras::Vector{Dict{Symbol,Int}})` in `src/types.jl`.
  - `_hyperbolic_catalysis(m::Union{Mechanism, AllostericMechanism}) -> Bool` in `src/mechanism_enumeration.jl`.
  - `_requires_hyperbolic_catalysis(::Mechanism) = false`, `_requires_hyperbolic_catalysis(::AllostericMechanism) = true` in `src/mechanism_enumeration.jl`.
  - `_all_reach(n::Int, edges, root::Int, fixed) -> Bool` (private helper of the predicate).

- [ ] **Step 1: Write the failing predicate tests**

Insert immediately after the closing `end` of the `@testset "_flux_carrying_groups"` block (search for `@testset "_re_segment_count"` and insert before it):

```julia
@testset "_hyperbolic_catalysis" begin
    # Ordered SS bi-bi: every metabolite binds on one edge and every segment is
    # a single form, so no denominator term carries a concentration twice.
    ordered = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E + Q <--> E(Q)
            E(Q) + P <--> E(P, Q)
            E(A, B) <--> E(P, Q)
        end
    end)
    @test EnzymeRates._hyperbolic_catalysis(ordered)

    # Random SS bi-bi: A binds on E → E(A) and on E(B) → E(A, B); both edges lie
    # in one tree toward E(A, B), so the denominator carries A².
    random_ss = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> E(A), E(B) + A <--> E(A, B))
            (E + B <--> E(B), E(A) + B <--> E(A, B))
            (E + P <--> E(P), E(Q) + P <--> E(P, Q))
            (E + Q <--> E(Q), E(P) + Q <--> E(P, Q))
            E(A, B) <--> E(P, Q)
        end
    end)
    @test !EnzymeRates._hyperbolic_catalysis(random_ss)

    # Random RE bi-bi with the chemistry step as the only SS step: one segment,
    # so the equation is the rapid-equilibrium law, degree 1 throughout.
    random_re = EnzymeRates.Mechanism(@enzyme_mechanism begin
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
    @test EnzymeRates._hyperbolic_catalysis(random_re)

    # Random RE bi-bi with the A group flipped: the segment {E(A), E(A, B)} has
    # E(A, B) carrying B beyond its bottom E(A), and the edge E(B) → E(A, B) into
    # it leaves a form carrying B, so the term rooted there carries B².
    random_flip_a = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> E(A), E(B) + A <--> E(A, B))
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
            (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
            (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
            E(A, B) <--> E(P, Q)
        end
    end)
    @test !EnzymeRates._hyperbolic_catalysis(random_flip_a)

    # Ordered SS bi-bi with substrate A as a dead end on E(Q): the derived
    # denominator carries A² (substrate inhibition), but the dead-end group is
    # not flux-carrying and the predicate ignores it.
    dead_end = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E + Q <--> E(Q)
            E(Q) + P <--> E(P, Q)
            E(A, B) <--> E(P, Q)
            E(Q) + A ⇌ E(A, Q)
        end
    end)
    @test EnzymeRates._hyperbolic_catalysis(dead_end)

    # Uni-bi with random SS product release: P binds on E → E(P) and on
    # E(Q) → E(P, Q), both toward E(P, Q), so the denominator carries P².
    unibi_ss = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P, Q)
            (E + P <--> E(P), E(Q) + P <--> E(P, Q))
            (E + Q <--> E(Q), E(P) + Q <--> E(P, Q))
        end
    end)
    @test !EnzymeRates._hyperbolic_catalysis(unibi_ss)

    # Ping-pong with one SS chemistry step per half-reaction: two segments,
    # each metabolite carried once per term.
    pingpong = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            E(Q) ⇌ E + Q
        end
    end)
    @test EnzymeRates._hyperbolic_catalysis(pingpong)

    # An allosteric mechanism is scored on its catalytic steps: uni-bi with the
    # P group at steady state carries Q² (the segment {E(P), E(P, Q)} plus the
    # edge E(Q) → E(P, Q)).
    allo_unibi_flip_p = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P, Q
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)                                        :: NonequalAI
            E(S) <--> E(P, Q)                                   :: NonequalAI
            (E + P <--> E(P), E(Q) + P <--> E(P, Q))            :: NonequalAI
            (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))                  :: NonequalAI
        end
    end)
    @test !EnzymeRates._hyperbolic_catalysis(allo_unibi_flip_p)

    @test !EnzymeRates._requires_hyperbolic_catalysis(ordered)
    @test EnzymeRates._requires_hyperbolic_catalysis(allo_unibi_flip_p)
end
```

- [ ] **Step 2: Run the enumeration test file to verify the new testset fails**

Run the focused command from Global Constraints. Expected: the new testset errors with `UndefVarError: _hyperbolic_catalysis not defined`; every other testset in the file passes.

- [ ] **Step 3: Extract `_re_segment_extras` from `_bottomless_re_segment` in `src/types.jl`**

Replace the whole `_bottomless_re_segment` function (keep its docstring verbatim above it) with the two functions below. The new helper's body is the exponent walk that `_bottomless_re_segment` used to hold; `_bottomless_re_segment` keeps only the side bookkeeping and the bottom test.

```julia
"""
Rapid-equilibrium segments of `steps` with each form's metabolite exponents
cleared to the lowest in its segment: `(species, segments, extras)`, where
`segments[k]` indexes `species` and `extras[i][x]` is how many more `x` form `i`
carries than the lowest form of its segment. Within a segment, rapid equilibrium
fixes each form's weight relative to any other as a monomial in concentrations,
one factor per metabolite bound or released along the RE path between them;
chemistry steps pass the exponents through unchanged.
"""
function _re_segment_extras(steps::Vector{Vector{Step}})
    species = Species[]
    for group in steps, s in group
        from_species(s) in species || push!(species, from_species(s))
        to_species(s) in species   || push!(species, to_species(s))
    end
    idx = Dict(sp => i for (i, sp) in enumerate(species))
    # RE adjacency with the exponent change of each metabolite from → to.
    adj = [Tuple{Int, Dict{Symbol, Int}}[] for _ in species]
    for group in steps, s in group
        is_equilibrium(s) || continue
        _, _, m_lhs, m_rhs = _step_sides(s)
        d = Dict{Symbol, Int}()
        for x in m_lhs; d[x] = get(d, x, 0) + 1; end
        for x in m_rhs; d[x] = get(d, x, 0) - 1; end
        a, b = idx[from_species(s)], idx[to_species(s)]
        push!(adj[a], (b, d))
        push!(adj[b], (a, Dict(x => -e for (x, e) in d)))
    end
    expo = Dict{Int, Dict{Symbol, Int}}()
    segments = Vector{Int}[]
    for root in eachindex(species)
        haskey(expo, root) && continue
        expo[root] = Dict{Symbol, Int}()
        segment = [root]
        queue = [root]
        while !isempty(queue)
            u = popfirst!(queue)
            for (v, d) in adj[u]
                haskey(expo, v) && continue
                expo[v] = mergewith(+, expo[u], d)
                push!(segment, v); push!(queue, v)
            end
        end
        for x in union(Set{Symbol}(), (keys(expo[i]) for i in segment)...)
            lowest = minimum(get(expo[i], x, 0) for i in segment)
            for i in segment
                expo[i][x] = get(expo[i], x, 0) - lowest
            end
        end
        push!(segments, segment)
    end
    extras = [filter(p -> p.second > 0, expo[i]) for i in eachindex(species)]
    species, segments, extras
end

function _bottomless_re_segment(steps::Vector{Vector{Step}})
    side = Dict{Symbol, Type}()
    for group in steps, s in group
        is_equilibrium(s) && is_binding(s) || continue
        side[name(bound_metabolite(s))] = typeof(bound_metabolite(s))
    end
    species, segments, extras = _re_segment_extras(steps)
    for segment in segments
        length(segment) < 2 && continue
        for S in (Product, Substrate)
            all(i -> any(x -> get(side, x, Nothing) === S, keys(extras[i])),
                segment) && return species[segment]
        end
    end
    nothing
end
```

Note: `_step_sides` is defined in `src/rate_eq_derivation.jl`, which is included after `src/types.jl`; that is already the case for the existing `_bottomless_re_segment`, so no include-order change is needed.

- [ ] **Step 4: Add the trait and predicate to `src/mechanism_enumeration.jl`**

Insert directly after the `_flux_carrying_groups` function (before the `_re_segment_count` docstring):

```julia
"""
Whether the enumerator may only give this mechanism type a catalytic scheme
that passes `_hyperbolic_catalysis`. True for every type that layers a
conformational equilibrium over its catalytic scheme.
"""
_requires_hyperbolic_catalysis(::Mechanism) = false
_requires_hyperbolic_catalysis(::AllostericMechanism) = true

"""
    _hyperbolic_catalysis(m) -> Bool

Whether the rate equation of `m`'s flux-carrying step graph has degree at most 1
in every substrate and product concentration. Dead-end groups (substrate,
product, or regulator binding off the catalytic cycle, `_flux_carrying_groups`)
are ignored: their inhibition terms are a separate source of concentration
powers that conformational mechanisms keep.

The equation is a sum over rapid-equilibrium (RE) segments and spanning
arborescences of the segment graph toward each segment. A denominator term is
the root segment's weight times the weight of every tree edge, so the exponent
of `X` in a term is the number of `X` the root segment's forms carry beyond the
segment's bottom form, plus, per tree edge, one if the step binds `X` in the
tree direction and one if the edge's source form carries `X` beyond its bottom
(`_re_segment_extras`). The degree exceeds 1 exactly when one edge scores 2, or
the root scores 1 and some arborescence toward it holds a scoring edge, or some
arborescence holds two scoring edges. An arborescence toward `S` containing
given edges exists iff every segment still reaches `S` once each given edge's
source keeps that edge as its only way out (`_all_reach`).
"""
function _hyperbolic_catalysis(m::Union{Mechanism, AllostericMechanism})
    flux = _flux_carrying_groups(m)
    groups = [steps(m)[g] for g in kinetic_groups(m) if flux[g]]
    species, segments, extras = _re_segment_extras(groups)
    idx = Dict(sp => i for (i, sp) in enumerate(species))
    seg_of = zeros(Int, length(species))
    for (k, segment) in enumerate(segments), i in segment
        seg_of[i] = k
    end
    # Directed segment-graph edges: source segment, target segment, source form,
    # metabolites bound in that direction.
    edges = Tuple{Int, Int, Int, Vector{Symbol}}[]
    for group in groups, s in group
        is_equilibrium(s) && continue
        _, _, m_lhs, m_rhs = _step_sides(s)
        a, b = idx[from_species(s)], idx[to_species(s)]
        seg_of[a] == seg_of[b] && continue
        push!(edges, (seg_of[a], seg_of[b], a, m_lhs))
        push!(edges, (seg_of[b], seg_of[a], b, m_rhs))
    end
    rxn = reaction(m)
    mets = vcat(Symbol[name(s) for s in substrates(rxn)],
                Symbol[name(p) for p in products(rxn)])
    n = length(segments)
    for x in mets
        score(e) = count(==(x), e[4]) + get(extras[e[3]], x, 0)
        carrying = [e for e in edges if score(e) > 0]
        any(e -> score(e) > 1, carrying) && return false
        roots = [k for (k, segment) in enumerate(segments)
                 if any(i -> get(extras[i], x, 0) > 0, segment)]
        for e in carrying, k in roots
            k != e[1] && _all_reach(n, edges, k, (e,)) && return false
        end
        for (p, e1) in enumerate(carrying), e2 in carrying[p + 1:end]
            e1[1] == e2[1] && continue
            any(k -> k != e1[1] && k != e2[1] && _all_reach(n, edges, k, (e1, e2)),
                1:n) && return false
        end
    end
    true
end

"""
Whether every segment of the segment graph reaches `root` when each edge in
`fixed` is its source segment's only way out. A digraph has a spanning
arborescence toward `root` iff every vertex reaches `root`, and with the fixed
edges as their sources' only exits every such arborescence contains them.
"""
function _all_reach(n::Int, edges, root::Int, fixed)
    pinned = Dict(e[1] => e[2] for e in fixed)
    into = [Int[] for _ in 1:n]
    for (u, v, _, _) in edges
        get(pinned, u, v) == v && push!(into[v], u)
    end
    seen = falses(n)
    seen[root] = true
    queue = [root]
    while !isempty(queue)
        v = popfirst!(queue)
        for u in into[v]
            seen[u] || (seen[u] = true; push!(queue, u))
        end
    end
    all(seen)
end
```

- [ ] **Step 5: Run the enumeration test file to verify it passes**

Run the focused command. Expected: the new `_hyperbolic_catalysis` testset passes (10 tests), and every existing testset still passes, in particular `"Mechanism — flip past a bottomless RE segment"` and the `_bottomless_re_segment` coverage in `test/test_types.jl`. Also run the types file once: same command with `include("test/test_types.jl")` in place of the enumeration file. Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add src/types.jl src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add _hyperbolic_catalysis: degree of the catalytic scheme's rate equation from its segment graph

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R2a1VuQs4B3Xn8dPRkg8Mk"
```

---

### Task 2: Exactness test against the derived denominator

**Files:**
- Test: `test/test_mechanism_enumeration.jl` (new testset directly after the `_hyperbolic_catalysis` testset from Task 1)

**Interfaces:**
- Consumes: `_hyperbolic_catalysis(m)`, `_flux_carrying_groups(m)`, `_raw_symbolic_rate_polys(m::Mechanism, step_params, rename_map, subs::Vector{Symbol}, prods::Vector{Symbol}) -> (num::POLY, den::POLY, d_free::POLY)`, `_step_parameters(m::Mechanism)`, `_build_wegscheider_rename_map(m::Mechanism)`, `_state_rate_polys(am::AllostericMechanism, state::Symbol) -> (num, den, d_free)`, `_with_steps_and_cat_states(am, steps, tags)`, `cat_allo_states(am)`, `_eq_complexity(m)`, `init_mechanisms(rxn)`, `expand_mechanisms(mechs, rxn)`. A `POLY` is `Dict{Vector{Pair{Symbol,Int}}, Rational{Int}}`: keys are monomials, each a sorted vector of `symbol => exponent`.
- Produces: nothing new in `src/`.

- [ ] **Step 1: Write the exactness test**

Insert directly after the `_hyperbolic_catalysis` testset:

```julia
@testset "_hyperbolic_catalysis matches the derived denominator" begin
    # The structural predicate against the exponents of the derived denominator,
    # over every mechanism reachable from the seeds in a bounded number of
    # expansion levels. Dead-end groups are stripped before deriving, since the
    # predicate ignores them; the reactions declare no regulators, so stripping
    # never leaves the reaction naming a metabolite no step binds. For an
    # allosteric mechanism the predicate is compared with the A-state, and the
    # I-state is checked to be hyperbolic whenever the A-state is.
    _testhelper_poly_hyperbolic(p, mets) =
        all(e <= 1 for mono in keys(p) for (s, e) in mono if s in mets)
    _testhelper_mets(rxn) = vcat(
        Symbol[EnzymeRates.name(s) for s in EnzymeRates.substrates(rxn)],
        Symbol[EnzymeRates.name(p) for p in EnzymeRates.products(rxn)])
    function _testhelper_den_hyperbolic(m::EnzymeRates.Mechanism)
        rxn = EnzymeRates.reaction(m)
        subs = Symbol[EnzymeRates.name(s) for s in EnzymeRates.substrates(rxn)]
        prods = Symbol[EnzymeRates.name(p) for p in EnzymeRates.products(rxn)]
        _, den, _ = EnzymeRates._raw_symbolic_rate_polys(
            m, EnzymeRates._step_parameters(m),
            EnzymeRates._build_wegscheider_rename_map(m), subs, prods)
        _testhelper_poly_hyperbolic(den, _testhelper_mets(rxn))
    end
    function _testhelper_flux_only(m::EnzymeRates.Mechanism)
        flux = EnzymeRates._flux_carrying_groups(m)
        keep = [g for g in EnzymeRates.kinetic_groups(m) if flux[g]]
        EnzymeRates.Mechanism(EnzymeRates.reaction(m), EnzymeRates.steps(m)[keep])
    end
    function _testhelper_flux_only(am::EnzymeRates.AllostericMechanism)
        flux = EnzymeRates._flux_carrying_groups(am)
        keep = [g for g in EnzymeRates.kinetic_groups(am) if flux[g]]
        EnzymeRates._with_steps_and_cat_states(
            am, EnzymeRates.steps(am)[keep], EnzymeRates.cat_allo_states(am)[keep])
    end
    function _testhelper_levels(rxn, depth)
        M = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}
        level = M[m for m in EnzymeRates.init_mechanisms(rxn)]
        mechs = unique(level)
        for _ in 1:depth
            level = unique!(EnzymeRates.expand_mechanisms(level, rxn))
            append!(mechs, level)
        end
        unique!(mechs)
    end

    unibi = @enzyme_reaction begin
        substrates: S[AB]
        products: P[A], Q[B]
        oligomeric_state: 2
    end
    bibi = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        oligomeric_state: 2
    end
    pingpong = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
        oligomeric_state: 2
    end
    n_checked = 0
    n_nonhyperbolic = 0
    for (rxn, depth) in ((unibi, 2), (bibi, 1), (pingpong, 2))
        mets = _testhelper_mets(rxn)
        for m in _testhelper_levels(rxn, depth)
            EnzymeRates._eq_complexity(m) <= 337 || continue
            stripped = _testhelper_flux_only(m)
            structural = EnzymeRates._hyperbolic_catalysis(m)
            if m isa EnzymeRates.Mechanism
                @test structural == _testhelper_den_hyperbolic(stripped)
            else
                _, den_a, _ = EnzymeRates._state_rate_polys(stripped, :A)
                _, den_i, _ = EnzymeRates._state_rate_polys(stripped, :I)
                hyp_a = _testhelper_poly_hyperbolic(den_a, mets)
                @test structural == hyp_a
                @test _testhelper_poly_hyperbolic(den_i, mets) || !hyp_a
            end
            n_checked += 1
            structural || (n_nonhyperbolic += 1)
        end
    end
    # Both classes must be exercised for the comparison to mean anything.
    @test n_checked > 100
    @test n_nonhyperbolic > 0
    @test n_nonhyperbolic < n_checked
end
```

The three reactions are the file's `uni_bi_rxn`, `bi_bi_rxn`, and `bi_bi_pp_rxn` constants (top of the file) with `oligomeric_state: 2` added, so that promotion children exist; they are written inline here because the rule for this file is that a testset's fixtures are readable in place.

- [ ] **Step 2: Run the enumeration test file**

Run the focused command. Expected: the new testset passes with zero failures. Record its wall time from the `@testset` summary. If any `@test structural == ...` fails, the predicate has a defect: do NOT edit the test; report the failing mechanism (print it with `show`) and stop for review. If the testset takes more than 120 s, lower `bibi` to depth 1 with allosteric children excluded (`filter(m -> m isa EnzymeRates.Mechanism, ...)` on the bi-bi level-1 set only) and report the numbers.

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Test _hyperbolic_catalysis against derived denominator degrees over the enumeration

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R2a1VuQs4B3Xn8dPRkg8Mk"
```

---

### Task 3: Promotion refuses a non-hyperbolic parent

**Files:**
- Modify: `src/mechanism_enumeration.jl:1887` (`_expand_to_allosteric(m::Mechanism, rxn)`) and its docstring (starts near line 1854)
- Test: `test/test_mechanism_enumeration.jl` inside `@testset "_expand_to_allosteric"` (starts near line 2952), add two testsets before the closing `end` of that block (search for `@testset "distinguishability invariant on every emitted mechanism"` and insert after that testset's `end`)

**Interfaces:**
- Consumes: `_hyperbolic_catalysis(m)`, `_expand_to_allosteric(m::Mechanism, rxn::EnzymeReaction) -> Vector{AllostericMechanism}`.
- Produces: the same signature; an empty vector for a non-hyperbolic parent.

- [ ] **Step 1: Write the failing move tests**

```julia
@testset "Mechanism — non-hyperbolic catalytic scheme: no children" begin
    # Uni-bi with random SS product release carries P² and Q² in its own
    # denominator; a conformational equilibrium on top would stack a second
    # source of concentration powers, so the promotion emits nothing.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P, Q)
            (E + P <--> E(P), E(Q) + P <--> E(P, Q))
            (E + Q <--> E(Q), E(P) + Q <--> E(P, Q))
        end
    end)
    rxn = @enzyme_reaction begin
        substrates: S[AB]
        products: P[A], Q[B]
        oligomeric_state: 2
    end
    @test isempty(EnzymeRates._expand_to_allosteric(m, rxn))
end

@testset "Mechanism — substrate dead-end keeps the scheme hyperbolic" begin
    # Ordered SS bi-bi with A as a dead end on E(Q): the dead-end's A² term is
    # substrate inhibition, not random-order steady state, so the promotion
    # proceeds. Five binding groups (A, B, Q, P, the dead end), each subset
    # :OnlyA with the chemistry :OnlyA: 2^5 - 1 = 31 K-type children.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E + Q <--> E(Q)
            E(Q) + P <--> E(P, Q)
            E(A, B) <--> E(P, Q)
            E(Q) + A ⇌ E(A, Q)
        end
    end)
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        oligomeric_state: 2
    end
    children = EnzymeRates._expand_to_allosteric(m, rxn)
    @test length(children) == 31
    @test all(EnzymeRates._hyperbolic_catalysis, children)
end
```

If the second count comes out different from 31, do not adjust the number blindly: check whether `_onlya_haldane_violation` drops subsets for this parent (print which masks it rejects) and report; only set the count to what the move provably emits.

- [ ] **Step 2: Run the enumeration test file to verify the first new testset fails**

Expected: `"Mechanism — non-hyperbolic catalytic scheme: no children"` fails (the move still emits children); the dead-end testset passes or reports a count to investigate.

- [ ] **Step 3: Add the check to `_expand_to_allosteric`**

At the top of `function _expand_to_allosteric(m::Mechanism, rxn::EnzymeReaction)`, before `n_g = length(steps(m))`, add:

```julia
    _hyperbolic_catalysis(m) || return AllostericMechanism[]
```

And extend the docstring of that method: after the sentence ending "duplicate variants are removed." add a paragraph:

```
A parent whose catalytic scheme fails `_hyperbolic_catalysis` (random-order
steady-state binding, whose own equation carries concentration powers) emits no
children: a conformational mechanism only carries a hyperbolic catalytic scheme.
```

- [ ] **Step 4: Run the enumeration test file to verify it passes**

Expected: both new testsets pass; every other testset passes. Pay attention to `"expand_mechanisms reaches ≥2 distinct-metabolite catalytic :OnlyA"`, the `_expand_to_allosteric` counts (15 for ordered bi-bi; ping-pong and uni-uni counts), and the `Integration` testsets. If a pre-existing testset fails, report which and why; do not edit it.

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Promotion to allosteric refuses a catalytic scheme with concentration powers

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R2a1VuQs4B3Xn8dPRkg8Mk"
```

---

### Task 4: The flip move keeps a conformational mechanism hyperbolic

**Files:**
- Modify: `src/mechanism_enumeration.jl:1252-1268` (`_expand_re_to_ss`) and its docstring (starts near line 1229)
- Test: `test/test_mechanism_enumeration.jl` inside `@testset "_expand_re_to_ss (group-set flips)"` (starts near line 5788); add one testset after `"_expand_re_to_ss: ping-pong"` and before that block's closing `end`

**Interfaces:**
- Consumes: `_hyperbolic_catalysis(m)`, `_requires_hyperbolic_catalysis(m)`, `_expand_re_to_ss(m) -> Vector{typeof(m)}`.
- Produces: the same signature; for an `AllostericMechanism` parent only hyperbolic children.

- [ ] **Step 1: Write the failing move test**

```julia
    @testset "_expand_re_to_ss: an allosteric parent keeps a hyperbolic scheme" begin
        # Uni-bi with random product release. Flipping the S group alone splits
        # off {E(S)} and the equation stays degree 1. Flipping the P group splits
        # off {E(P), E(P, Q)}, whose term carries Q from E(P, Q) beyond its
        # bottom E(P) and again from the edge E(Q) → E(P, Q), so the equation
        # carries Q²; flipping Q mirrors it with P². A plain Mechanism emits all
        # three; an allosteric parent emits only the S flip.
        plain = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P, Q)
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
            end
        end)
        flipS = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S <--> E(S)
                E(S) <--> E(P, Q)
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
            end
        end)
        flipP = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P, Q)
                (E + P <--> E(P), E(Q) + P <--> E(P, Q))
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
            end
        end)
        flipQ = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P, Q)
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))
                (E + Q <--> E(Q), E(P) + Q <--> E(P, Q))
            end
        end)
        kids = EnzymeRates._expand_re_to_ss(plain)
        @test length(kids) == 3
        @test Set(kids) == Set([flipS, flipP, flipQ])

        allo = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products: P, Q
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S ⇌ E(S)                                :: NonequalAI
                E(S) <--> E(P, Q)                           :: NonequalAI
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))          :: NonequalAI
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))          :: NonequalAI
            end
        end)
        allo_flipS = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products: P, Q
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S <--> E(S)                             :: NonequalAI
                E(S) <--> E(P, Q)                           :: NonequalAI
                (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q))          :: NonequalAI
                (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))          :: NonequalAI
            end
        end)
        allo_kids = EnzymeRates._expand_re_to_ss(allo)
        @test length(allo_kids) == 1
        @test Set(allo_kids) == Set([allo_flipS])
    end
```

- [ ] **Step 2: Run the enumeration test file to verify the new testset fails**

Expected: the plain-Mechanism assertions pass (3 children); the allosteric assertions fail (3 children instead of 1).

- [ ] **Step 3: Filter the emitted children in `_expand_re_to_ss`**

Replace the last line of `_expand_re_to_ss`,

```julia
    typeof(m)[_with_steps(m, flipped_groups(sel)) for sel in sets]
```

with

```julia
    children = typeof(m)[_with_steps(m, flipped_groups(sel)) for sel in sets]
    _requires_hyperbolic_catalysis(m) ? filter(_hyperbolic_catalysis, children) : children
```

Extend the docstring of `_expand_re_to_ss`: after the sentence "…or through the extended set when no other order gains." add:

```
A parent that `_requires_hyperbolic_catalysis` (a conformational mechanism) keeps
only the children whose catalytic scheme passes `_hyperbolic_catalysis`. The
check runs on the emitted minimal sets, not inside the gain predicate: a
failing set would otherwise be extended with more flips, and more steady-state
steps never restore a hyperbolic equation, so its supersets need no visit.
```

- [ ] **Step 4: Run the enumeration test file to verify it passes**

Expected: the new testset passes; every other testset passes, in particular `"AllostericMechanism — :NonequalAI: 2 variants, tags preserved"` (uni-uni: both flips remain hyperbolic) and the 220-child seed pin `"_expand_re_to_ss: seed child count is unchanged (220 over bi-bi seeds)"` (that pin is over plain `Mechanism` seeds, so it must not move). If a pre-existing testset fails, report which and why; do not edit it.

- [ ] **Step 5: Time the move on a ter-ter allosteric parent**

Run once in a fresh Julia process (`julia --project`), pasting:

```julia
using EnzymeRates
am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
    substrates: A, B, C
    products: P, Q, R
    catalytic_multiplicity: 2
    catalytic_steps: begin
        (E + A ⇌ E(A), E(B) + A ⇌ E(A, B), E(C) + A ⇌ E(A, C),
         E(B, C) + A ⇌ E(A, B, C))                                   :: NonequalAI
        (E + B ⇌ E(B), E(A) + B ⇌ E(A, B), E(C) + B ⇌ E(B, C),
         E(A, C) + B ⇌ E(A, B, C))                                   :: NonequalAI
        (E + C ⇌ E(C), E(A) + C ⇌ E(A, C), E(B) + C ⇌ E(B, C),
         E(A, B) + C ⇌ E(A, B, C))                                   :: NonequalAI
        E(A, B, C) <--> E(P, Q, R)                                   :: NonequalAI
        E(Q, R) + P ⇌ E(P, Q, R)                                     :: NonequalAI
        E(R) + Q ⇌ E(Q, R)                                           :: NonequalAI
        E + R ⇌ E(R)                                                 :: NonequalAI
    end
end)
EnzymeRates._expand_re_to_ss(am)
@time EnzymeRates._expand_re_to_ss(am)
@time EnzymeRates._hyperbolic_catalysis(am)
```

Report both timings (second call of each, after compilation). The predicate must be well under a millisecond; the move's time must be dominated by what it did before this task. If the predicate exceeds a millisecond, report it and stop for review rather than optimizing.

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Flip move keeps a conformational mechanism's catalytic scheme hyperbolic

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R2a1VuQs4B3Xn8dPRkg8Mk"
```

---

### Task 5: Docs and version

**Files:**
- Modify: `docs/src/identify/enumeration_engine.md` (sections "1. Flip rapid-equilibrium groups to steady state", "4. Promote a non-allosteric to allosteric mechanism", and "Modeling choices")
- Modify: `docs/src/developer.md` (section "Enumeration engine architecture", the paragraph ending "…nothing in `src/` estimates identifiability numerically.")
- Modify: `Project.toml:3` (`version = "0.6.0"` → `version = "0.7.0"`)

**Interfaces:** none.

- [ ] **Step 1: Document the flip move's rule**

In `docs/src/identify/enumeration_engine.md`, section "1. Flip rapid-equilibrium groups to steady state", after the paragraph beginning "Two kinds of group never flip." add:

```
An allosteric parent also drops every child whose catalytic scheme would carry a
concentration to a power of its own — the random-order steady-state pattern in
which a metabolite binds on two steps that one King–Altman tree can hold. See
**Conformational mechanisms carry hyperbolic catalytic schemes** under
[Modeling choices](@ref).
```

Check that `[Modeling choices](@ref)` resolves: the heading is `## Modeling choices` in the same file. If Documenter's cross-reference form differs elsewhere in the file (search for `](@ref)`), copy that form.

- [ ] **Step 2: Document the promotion move's rule**

In section "4. Promote a non-allosteric to allosteric mechanism", after "No-op on an already allosteric input." add:

```
Also a no-op on a parent whose catalytic scheme already carries concentration
powers (a random-order scheme with steady-state binding); see
**Conformational mechanisms carry hyperbolic catalytic schemes** under
[Modeling choices](@ref).
```

- [ ] **Step 3: State the modeling choice**

In section "Modeling choices", after the paragraph "**Dead-end branches stay at equilibrium.**…" add:

```
**Conformational mechanisms carry hyperbolic catalytic schemes.** An MWC
conformational equilibrium and a random-order steady-state catalytic scheme
each put a metabolite concentration to a power in the rate equation, and
sigmoidal data cannot tell the two sources apart. The moves therefore never
combine them: promotion to allosteric skips a parent whose catalytic scheme
already carries a power, and the flip move on an allosteric parent drops a
child that would introduce one. The test is structural
(`_hyperbolic_catalysis`): the equation's degree in a metabolite is read from the
rapid-equilibrium segment graph without deriving it. Competitive inhibition by a
substrate or product is a separate source of powers and stays allowed inside an
allosteric mechanism; the test ignores dead-end groups. Hand-written mechanisms
are not subject to the rule: an `@allosteric_mechanism` with random-order
steady-state binding still derives and fits.
```

- [ ] **Step 4: Document the predicate on the Developer page**

In `docs/src/developer.md`, after the paragraph ending "…nothing in `src/` estimates identifiability numerically." add:

```
Conformational mechanism types declare `_requires_hyperbolic_catalysis` (true for
`AllostericMechanism`), and the moves that can give such a type a non-hyperbolic
catalytic scheme consult `_hyperbolic_catalysis`: `_expand_to_allosteric` refuses
the parent, `_expand_re_to_ss` filters its emitted children. The predicate scores
the rapid-equilibrium segment graph of the flux-carrying groups: for each
metabolite, a directed steady-state edge scores one if the step binds it in that
direction and one if the edge's source form carries it beyond its segment's
bottom form (`_re_segment_extras`, shared with `_bottomless_re_segment`), and a
segment scores one if any of its forms does. The equation's degree in the
metabolite is the best score over root segments and spanning arborescences
toward the root; the predicate decides "at most 1" with reachability checks
(`_all_reach`) rather than a max-arborescence solve, since a weight-2 pattern is
one edge scoring 2, a scoring root plus a compatible scoring edge, or two
compatible scoring edges. The exactness of the predicate against the derived
denominator is pinned by a test over the uni-bi, bi-bi, and ping-pong
enumerations. A future conformational type (KNF, mnemonic, slow isomerization)
adds one `_requires_hyperbolic_catalysis` method and calls the predicate in its
promotion move.
```

- [ ] **Step 5: Bump the version**

In `Project.toml`, change line 3 to `version = "0.7.0"`.

- [ ] **Step 6: Build the docs to check cross-references**

Run: `julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate(); include("docs/make.jl")'` from the repo root, if a `docs/Project.toml` exists (check with `ls docs`). Expected: no "cross-reference" or "@ref" warnings for the edited pages. If the docs environment is not set up and the build cannot run within a few minutes, skip it and say so in the commit message body; do not leave it unreported.

- [ ] **Step 7: Commit**

```bash
git add docs/src/identify/enumeration_engine.md docs/src/developer.md Project.toml
git commit -m "Document the hyperbolic-catalysis rule for conformational mechanisms (v0.7.0)

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01R2a1VuQs4B3Xn8dPRkg8Mk"
```

---

### Task 6: Full test suite

**Files:** none modified unless the suite fails.

- [ ] **Step 1: Run the full suite in the background and poll**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
nohup julia --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])' > /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/4905eb94-aab8-4f1b-8edf-37cd1f79ec43/scratchpad/fullsuite.log 2>&1 &
```

Poll with a bounded loop (never a background monitor): every 60 s for up to 25 minutes, `tail -3` the log; stop when it prints the `Test Summary` table or `Testing EnzymeRates tests passed`. Make sure no other Julia process is running first (`pgrep -a julia`).

- [ ] **Step 2: Read the result**

Expected: `Pass` count with 0 `Fail` and 0 `Error`, Aqua and JET included. If anything fails, report the testset name and output verbatim; do not delete or weaken any test.

- [ ] **Step 3: Nothing to commit**

If the suite is green, the branch is complete: 5 commits on top of the two spec commits. Report the final `git log --oneline main..HEAD`.
