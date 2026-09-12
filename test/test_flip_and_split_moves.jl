# ABOUTME: Tests for the RE→SS group-set flip move and the context-bipartition split move:
# ABOUTME: eligibility, gain proofs, minimal sets, reachability, once-per-parent counter.
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
                      if any(s -> length(EnzymeRates.bound(
                                       EnzymeRates.to_species(s))) == 2, grp))
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
            E(S) + I ⇌ E(I, S)    :: EqualAI
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

    # Allosteric: measured on the A-state projection.
    am = EnzymeRates.AllostericMechanism(EnzymeRates.@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)   :: EqualAI
            E(S) <--> E(P) :: EqualAI
            E + P ⇌ E(P)   :: EqualAI
        end
    end)
    @test EnzymeRates._re_segment_count(am) == 1
    am_groups = EnzymeRates.steps(am)
    am_g = findfirst(grp -> all(EnzymeRates.is_equilibrium, grp), am_groups)
    am_flipped = EnzymeRates._with_steps(am,
                                        EnzymeRates._flip_group_to_ss(am_groups, am_g))
    @test EnzymeRates._re_segment_count(am_flipped) == 2
end

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
    binds_a(s) = (bm = EnzymeRates.bound_metabolite(s);
                  bm !== nothing && EnzymeRates.name(bm) == :A)
    a_groups = [g for (g, grp) in enumerate(groups)
                if length(grp) == 1 && binds_a(only(grp))]
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
    ss_steps(c) = Set(s for grp in EnzymeRates.steps(c) for s in grp
                      if !EnzymeRates.is_equilibrium(s))
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
    # Contexts are ordered by ligand role (Product before Substrate) then
    # name: bps[1] is the P-context division, bps[2] the B-context division.
    p_step = only(s for s in a_group
                  if any(b -> b isa EnzymeRates.Product,
                         EnzymeRates.bound(EnzymeRates.from_species(s))))
    b_step = only(s for s in a_group
                  if any(b -> b isa EnzymeRates.Substrate,
                         EnzymeRates.bound(EnzymeRates.from_species(s))))
    @test bps[1][2] == [p_step]
    @test bps[2][2] == [b_step]
    # _context_form: canonical RE binding puts the metabolite on to_species,
    # so the context form is from_species.
    @test EnzymeRates._context_form(first(a_group)) ==
          EnzymeRates.from_species(first(a_group))
    # _context_form: an SS dissociation step whose bound metabolite is in
    # neither endpoint's bound list (the Segel ping-pong step shape) puts
    # the metabolite on to_species.
    ping_pong_step = EnzymeRates.Step(
        EnzymeRates.Species([EnzymeRates.Substrate(:A)], :E),
        EnzymeRates.Species(EnzymeRates.Metabolite[], :F),
        EnzymeRates.Product(:P), false)
    @test EnzymeRates._context_form(ping_pong_step) ==
          EnzymeRates.to_species(ping_pong_step)
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

"Finite-difference rank of ∂v/∂log θ over the fitted parameters (test oracle only)."
function _identifiable_rank(m; npts = 60, ndraws = 3, h = 1e-5)
    em = EnzymeRates.compile_mechanism(m)
    fp = collect(EnzymeRates.fitted_params(em))
    cm = m isa EnzymeRates.Mechanism ? m : EnzymeRates._state_mechanism(m, :A)
    mets = sort!(collect(EnzymeRates._concentration_symbols(cm)))
    # Seeded from the rendered equation: a struct hash mixes in objectid and
    # would not reproduce a failure in a later session.
    rng = MersenneTwister(hash(rate_equation_string(em)) % 2^31)
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

@testset "_expand_split_kinetic_group: random-order bi-bi frees independent constants" begin
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
    for (g, grp) in enumerate(EnzymeRates.steps(m)),
        bp in EnzymeRates._context_bipartitions(grp)
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
            (E + P ⇌ E(P), E(Q) + P ⇌ E(P, Q), E(A) + P ⇌ E(A, P),
             E(B) + P ⇌ E(B, P), E(A, B) + P ⇌ E(A, B, P))
            (E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))
            E(A, B) <--> E(P, Q)
        end
    end)
    binder(grp) = EnzymeRates.name(EnzymeRates.bound_metabolite(first(grp)))
    a_group = only(grp for grp in EnzymeRates.steps(m)
                   if length(grp) == 4 && binder(grp) == :A)
    bps = EnzymeRates._context_bipartitions(a_group)
    # Context B and context P each divide the four A steps two and two, a
    # division no carve of a single step can produce.
    @test count(bp -> length(bp[1]) == 2 && length(bp[2]) == 2, bps) == 2
    gi = findfirst(==(a_group), EnzymeRates.steps(m))
    even_bp = first(bp for bp in bps if length(bp[1]) == 2)
    child = EnzymeRates._apply_bipartitions(m, [(gi, even_bp)])
    @test count(grp -> length(grp) == 2 && Set(grp) ⊆ Set(a_group),
                EnzymeRates.steps(child)) == 2
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
