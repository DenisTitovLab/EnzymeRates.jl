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
