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

@testset "conformation_steps" begin
    M = ER.Mechanisms
    am = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)          :: EqualAI
            E(S) <--> E(P)        :: OnlyA
            E(P) ⇌ E + P          :: EqualAI
            E + I ⇌ E(I)          :: OnlyI
            E(I) + P ⇌ E(I, P)    :: EqualAI
        end
    end)
    m = M.catalytic(am)
    # Canonical form numbers the forms E, E(I), E(I, P), E(P), E(S) and stores the steps
    # E + I → E(I), E + P → E(P), E + S → E(S), E(I) + P → E(I, P), E(S) → E(P).
    @test [M.ends(m, k) for k in M.steps(m)] == [(1, 2), (1, 4), (1, 5), (2, 3), (5, 4)]
    # The active conformation drops the :OnlyI binding of I, which strands E(I) and so
    # also the binding of P to it.
    @test M.conformation_steps(am, :active) == BitVector([0, 1, 1, 0, 1])
    # The inactive conformation drops the :OnlyA chemistry and keeps every form.
    @test M.conformation_steps(am, :inactive) == BitVector([1, 1, 1, 1, 0])
    # A plain mechanism has only the active conformation.
    @test M.conformation_steps(m, :active) == trues(5)
    @test M.conformation_steps(m, :inactive) == falses(5)
    @test_throws ArgumentError M.conformation_steps(am, :I)
    @test_throws ArgumentError M.conformation_steps(m, :A)
end

"""Whether a cycle of the steps of graph `m` that `on` keeps takes up some metabolite name
net: for some name, the uptake of the kept steps is no difference of form potentials,
which holds when appending the uptake column to their incidence matrix raises its rank."""
function _testhelper_net_uptake_cycle(m, on::BitVector)
    M = ER.Mechanisms
    kept = filter(k -> on[k], collect(M.steps(m)))
    incidence = zeros(length(kept), length(M.species(m)))
    for (r, k) in enumerate(kept)
        a, b = M.ends(m, k)
        incidence[r, a] -= 1
        incidence[r, b] += 1
    end
    names = unique(reduce(vcat, [[M.metabolite_names(M.taken_up(m, k));
                                  M.metabolite_names(M.given_off(m, k))]
                                 for k in M.steps(m)]))
    any(names) do x
        uptake = Int[M.copies(M.taken_up(m, k), x) - M.copies(M.given_off(m, k), x)
                     for k in kept]
        rank([incidence uptake]) > rank(incidence)
    end
end


@testset "turnover cycle: every subset of steps of five mechanisms" begin
    M = ER.Mechanisms
    # Each mechanism with the cycles of its steps that turn the reaction over, written as
    # step texts (`_testhelper_edge_text`). A set of steps holds a turnover cycle exactly
    # when it holds all the steps of one of these.
    cases = [
        (@enzyme_mechanism(begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E + P
            end
        end), [["E + S ⇌ ES", "E + P <--> ES"]]),
        # Uni-uni through a second free form F.
        (@enzyme_mechanism(begin
            substrates: S
            products: P
            steps: begin
                E + S <--> F(P)
                F(P) ⇌ F + P
                E <--> F
            end
        end), [["E + S <--> FP", "F + P ⇌ FP", "F <--> E"]]),
        # Random-order bi-bi: the square of the two binding orders takes up A and B once
        # each way, so it turns nothing over; each binding order with the chemistry and
        # the releases does.
        (@enzyme_mechanism(begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
                (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                E(A, B) <--> E(P, Q)
                E(P, Q) <--> E(Q) + P
                E(Q) <--> E + Q
            end
        end), [["E + A ⇌ EA", "EA + B ⇌ EAB", "EAB <--> EPQ", "EQ + P <--> EPQ",
                "E + Q <--> EQ"],
               ["E + B ⇌ EB", "EB + A ⇌ EAB", "EAB <--> EPQ", "EQ + P <--> EPQ",
                "E + Q <--> EQ"]]),
        # Substrate inhibition: the dead-end binding of a second S lies on no cycle.
        (@enzyme_mechanism(begin
            substrates: S
            products: P
            modifiers: S
            steps: begin
                E + S ⇌ E(S)
                E(S) + S::Mod ⇌ E(S, S::Mod)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end), [["E + S ⇌ ES", "ES <--> EP", "E + P ⇌ EP"]]),
        # A Theorell–Chance step beside the merged complex E(P, Q): the triangle of the
        # two routes from E(A) to E(Q) takes up B and gives off P once each way.
        (@enzyme_mechanism(begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) + B <--> E(P, Q)
                E(P, Q) <--> E(Q) + P
                E(Q) <--> E + Q
                E(A) + B <--> E(Q) + P
            end
        end), [["E + A ⇌ EA", "EA + B <--> EPQ", "EQ + P <--> EPQ", "E + Q <--> EQ"],
               ["E + A ⇌ EA", "EA + B <--> EQ + P", "E + Q <--> EQ"]]),
    ]
    for (em, cycles) in cases
        m = ER.Mechanism(em)
        n = length(M.steps(m))
        at = [[_testhelper_edge(m, text) for text in cycle] for cycle in cycles]
        for bits in 0:(2^n - 1)
            on = BitVector([isodd(bits >> (k - 1)) for k in 1:n])
            expected = any(ks -> all(k -> on[k], ks), at)
            @test (findall(on), M.has_turnover_cycle(m, on)) == (findall(on), expected)
        end
    end
end

@testset "turnover cycle: MECHANISM_TEST_SPECS" begin
    M = ER.Mechanisms
    for spec in MECHANISM_TEST_SPECS
        mech = _testhelper_lift(spec.mechanism)
        m = M.catalytic(mech)
        n = length(M.steps(m))
        # Every mechanism turns the reaction over, and none at rapid equilibrium alone.
        @test (spec.name, M.has_turnover_cycle(m, trues(n))) == (spec.name, true)
        @test (spec.name, M.has_turnover_cycle(m, M._flags(m))) == (spec.name, false)
        subsets = [[Int[]]; [[i] for i in 1:n]; [[i, j] for i in 1:n for j in (i + 1):n]]
        masks = [BitVector(k in s for k in 1:n) for s in subsets]
        mech isa ER.AllostericMechanism &&
            push!(masks, M.conformation_steps(mech, :active),
                  M.conformation_steps(mech, :inactive))
        for on in masks
            @test (spec.name, findall(on), M.has_turnover_cycle(m, on)) ==
                  (spec.name, findall(on), _testhelper_net_uptake_cycle(m, on))
        end
    end
end

@testset "turnover cycle: malformed input" begin
    M = ER.Mechanisms
    # Hand-written graphs whose cycles do not each run the reaction a whole number of
    # times: a cycle takes up A and Q and gives off nothing, or takes up B alone. The
    # turnover test weighs each step by its uptake of the first substrate, A, so it finds
    # the first cycle and not the second.
    aq = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> E(Q) + P
            E(Q) <--> E + Q
            E(A) + Q <--> E(A, Q)
            E(A, Q) <--> E
        end
    end)
    on = BitVector(k in [_testhelper_edge(aq, t)
                         for t in ("E + A <--> EA", "EA + Q <--> EAQ", "EAQ <--> E")]
                   for k in M.steps(aq))
    @test M.has_turnover_cycle(aq, on)
    b_only = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> E(Q) + P
            E(Q) <--> E + Q
            E + B <--> E(B)
            E(B) <--> E
        end
    end)
    on = BitVector(k in [_testhelper_edge(b_only, t)
                         for t in ("E + B <--> EB", "EB <--> E")]
                   for k in M.steps(b_only))
    @test !M.has_turnover_cycle(b_only, on)
    # A mask with fewer flags than steps is out of bounds.
    @test_throws BoundsError M.has_turnover_cycle(aq, trues(2))
end
