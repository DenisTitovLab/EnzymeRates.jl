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
