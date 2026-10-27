# ABOUTME: Tests of the graph accessors of a mechanism (species, steps, kinetic groups,
# ABOUTME: metabolite multisets, allosteric parts): exact values, and no allocation.

const ER = EnzymeRates

"""Bytes the call `f(args...)` allocates, after a first call compiles it."""
function _testhelper_allocated(f, args...)
    f(args...)
    @allocated f(args...)
end

"""Test that the scalar graph accessors of `m`, on form `x` and step `k`, allocate
nothing."""
function _testhelper_test_graph_accessors_allocate_nothing(m, x, k)
    M = ER.Mechanisms
    s, t = M.bound(m, x), M.taken_up(m, k)
    @test _testhelper_allocated(M.species, m) == 0
    @test _testhelper_allocated(M.bound, m, x) == 0
    @test _testhelper_allocated(M.conformation, m, x) == 0
    @test _testhelper_allocated(M.steps, m) == 0
    @test _testhelper_allocated(M.ends, m, k) == 0
    @test _testhelper_allocated(M.taken_up, m, k) == 0
    @test _testhelper_allocated(M.given_off, m, k) == 0
    @test _testhelper_allocated(M.is_rapid_equilibrium, m, k) == 0
    @test _testhelper_allocated(M.step_kind, m, k) == 0
    @test _testhelper_allocated(M.group, m, k) == 0
    @test _testhelper_allocated(M.ngroups, m) == 0
    @test _testhelper_allocated(M.copies, s, :S) == 0
    @test _testhelper_allocated(isempty, s) == 0
    @test _testhelper_allocated(∪, s, t) == 0
    @test _testhelper_allocated(∩, s, t) == 0
end

@testset "graph accessors: random-order bi-bi" begin
    M = ER.Mechanisms
    m = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
            E(A, B) <--> E(P, Q)
            E(P, Q) <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end)
    @test M.catalytic(m) === m

    # Form numbers follow canonical form: free E, then the forms by their bound
    # metabolites read as sorted lists, so E(A, B) comes before E(B).
    @testset "species" begin
        @test M.species(m) === 1:6
        @test [M.metabolite_names(M.bound(m, x)) for x in M.species(m)] ==
              [Symbol[], [:A], [:A, :B], [:B], [:P, :Q], [:Q]]
        @test [M.conformation(m, x) for x in M.species(m)] == fill(:E, 6)
    end

    # Step numbers and directions follow canonical form: each step is stored
    # physical-forward, so a product release is stored as the binding it reverses, and the
    # steps are sorted by their stored ends.
    @testset "steps" begin
        @test M.steps(m) === 1:7
        @test [M.ends(m, k) for k in M.steps(m)] ==
              [(1, 2), (1, 4), (1, 6), (2, 3), (3, 5), (4, 3), (6, 5)]
        @test [M.metabolite_names(M.taken_up(m, k)) for k in M.steps(m)] ==
              [[:A], [:B], [:Q], [:B], Symbol[], [:A], [:P]]
        @test all(isempty(M.given_off(m, k)) for k in M.steps(m))
        @test [M.is_rapid_equilibrium(m, k) for k in M.steps(m)] ==
              [true, true, false, true, false, true, false]
        @test [M.step_kind(m, k) for k in M.steps(m)] ==
              [:binding, :binding, :binding, :binding, :isomerization, :binding, :binding]
    end

    # Kinetic groups are numbered by their first step in step order.
    @testset "groups" begin
        @test M.ngroups(m) == 5
        @test [M.group(m, k) for k in M.steps(m)] == [1, 2, 3, 2, 4, 1, 5]
        @test [M.members(m, g) for g in 1:5] == [[1, 6], [2, 4], [3], [5], [7]]
    end

    @testset "metabolite multisets" begin
        EA, EAB, EB = M.bound(m, 2), M.bound(m, 3), M.bound(m, 4)
        @test EAB isa M.Metabolites
        @test (M.copies(EAB, :A), M.copies(EAB, :B), M.copies(EAB, :P)) == (1, 1, 0)
        @test !isempty(EA) && isempty(M.bound(m, 1))
        @test EA ∪ EB == EAB
        # A multiset union holds each metabolite as often as the fuller side does.
        @test EA ∪ EAB == EAB
        @test EA ∩ EAB == EA
        @test isempty(EA ∩ EB)
        @test M.taken_up(m, 1) == EA
    end

    @testset "scalar accessors allocate nothing" begin
        _testhelper_test_graph_accessors_allocate_nothing(m, 3, 4)
    end
end

@testset "graph accessors: allosteric uni-uni" begin
    M = ER.Mechanisms
    am = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)      :: NonequalAI
            E(S) <--> E(P)    :: OnlyA
            E(P) ⇌ E + P      :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end)
    m = M.catalytic(am)
    @test m isa ER.Mechanism

    # Canonical form numbers the forms E, E(P), E(S) and stores the steps E + P → E(P),
    # E + S → E(S), E(S) → E(P); each step is its own kinetic group, in that order.
    @testset "species" begin
        @test M.species(m) === 1:3
        @test [M.metabolite_names(M.bound(m, x)) for x in M.species(m)] ==
              [Symbol[], [:P], [:S]]
        @test [M.conformation(m, x) for x in M.species(m)] == [:E, :E, :E]
    end

    @testset "steps" begin
        @test M.steps(m) === 1:3
        @test [M.ends(m, k) for k in M.steps(m)] == [(1, 2), (1, 3), (3, 2)]
        @test [M.metabolite_names(M.taken_up(m, k)) for k in M.steps(m)] ==
              [[:P], [:S], Symbol[]]
        @test all(isempty(M.given_off(m, k)) for k in M.steps(m))
        @test [M.is_rapid_equilibrium(m, k) for k in M.steps(m)] == [true, true, false]
        @test [M.step_kind(m, k) for k in M.steps(m)] ==
              [:binding, :binding, :isomerization]
    end

    @testset "groups" begin
        @test M.ngroups(m) == 3
        @test [M.group(m, k) for k in M.steps(m)] == [1, 2, 3]
        @test [M.members(m, g) for g in 1:3] == [[1], [2], [3]]
    end

    @testset "allosteric parts" begin
        # Each tag follows its group to the canonical group number.
        @test [M.tag(am, g) for g in 1:3] == [:EqualAI, :NonequalAI, :OnlyA]
        @test M.catalytic_multiplicity(am) == 2
        @test M.regulatory_sites(am) == [M.RegulatorySite([:R], 2, [:OnlyI])]
    end

    @testset "scalar accessors allocate nothing" begin
        _testhelper_test_graph_accessors_allocate_nothing(m, 3, 2)
        @test _testhelper_allocated(M.catalytic, am) == 0
        @test _testhelper_allocated(M.tag, am, 3) == 0
        @test _testhelper_allocated(M.catalytic_multiplicity, am) == 0
        @test _testhelper_allocated(M.regulatory_sites, am) == 0
    end
end

@testset "graph accessors: conformations and chemistry" begin
    M = ER.Mechanisms
    # Uni-uni through a second free form F: S binds straight into the product-bound form.
    m = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> F(P)
            F(P) ⇌ F + P
            E <--> F
        end
    end)
    # Canonical form numbers the forms E, F, F(P). It stores the isomerization from F,
    # where only products enter or leave, to E, where only substrates do: against its
    # writing.
    @test [M.conformation(m, x) for x in M.species(m)] == [:E, :F, :F]
    @test [M.ends(m, k) for k in M.steps(m)] == [(1, 3), (2, 1), (2, 3)]
    @test [M.step_kind(m, k) for k in M.steps(m)] ==
          [:chemistry, :isomerization, :binding]
    @test [M.metabolite_names(M.taken_up(m, k)) for k in M.steps(m)] ==
          [[:S], Symbol[], [:P]]
    @test _testhelper_allocated(M.step_kind, m, 1) == 0
end

@testset "graph accessors: a metabolite and its modifier copy" begin
    M = ER.Mechanisms
    # Substrate inhibition: a second S binds E(S) as a modifier.
    em = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: S
        steps: begin
            E + S ⇌ E(S)
            E(S) + S::Mod ⇌ E(S, S::Mod)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    m = ER.Mechanism(em)
    # Canonical form numbers the forms E, E(P), E(S), E(S, S::Mod) and sorts the steps
    # E + P → E(P), E + S → E(S), E(S) → E(P), E(S) + S::Mod → E(S, S::Mod).
    ESS = M.bound(m, 4)
    # A reactant and its modifier copy count under one name, as they share a concentration.
    @test M.metabolite_names(ESS) == [:S, :S]
    @test M.copies(ESS, :S) == 2 && M.copies(ESS, :P) == 0
    @test M.ends(m, 4) == (3, 4)
    @test M.metabolite_names(M.taken_up(m, 4)) == [:S]
    @test M.step_kind(m, 4) === :binding
    # The modifier copy is an entry of its own, so it adds to the bound S.
    @test M.copies(M.bound(m, 3) ∪ M.taken_up(m, 4), :S) == 2
    @test isempty(M.bound(m, 3) ∩ M.taken_up(m, 4))

    # Multisets compare by content, so two reads of one compiled mechanism agree.
    other = ER.Mechanism(em)
    @test M.bound(other, 4) == ESS && hash(M.bound(other, 4)) == hash(ESS)
    @test M.bound(other, 3) != ESS
    # They combine only over one metabolite table.
    uniuni = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
    @test_throws "combine the metabolites of one mechanism" ESS ∪ M.bound(uniuni, 3)
    @test_throws "combine the metabolites of one mechanism" ESS ∩ M.bound(uniuni, 3)
end
