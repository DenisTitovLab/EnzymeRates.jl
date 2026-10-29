# ABOUTME: Tests of what a mechanism's stored graph holds (nodes, steps, kinetic groups,
# ABOUTME: allosteric parts), `step_kind`, `groups` and `conformation_steps`.

const ER = EnzymeRates

"""The names of the metabolites in `list`, in order."""
_testhelper_names(list) = Symbol[x.name for x in list]

"""Bytes the call `f(args...)` allocates, after a first call compiles it."""
function _testhelper_allocated(f, args...)
    f(args...)
    @allocated f(args...)
end

"""Test that `step_kind` and `catalytic` of `m`, on step `k`, allocate nothing."""
function _testhelper_test_graph_reads_allocate_nothing(m, k)
    @test _testhelper_allocated(ER.Mechanisms.step_kind, m, k) == 0
    @test _testhelper_allocated(ER.Mechanisms.catalytic, m) == 0
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
        @test length(m.nodes) == 6
        @test [_testhelper_names(n.bound) for n in m.nodes] ==
              [Symbol[], [:A], [:A, :B], [:B], [:P, :Q], [:Q]]
        @test [n.conformation for n in m.nodes] == fill(:E, 6)
    end

    # Step numbers and directions follow canonical form: each step is stored
    # physical-forward, so a product release is stored as the binding it reverses, and the
    # steps are sorted by their stored ends.
    @testset "steps" begin
        @test length(m.steps) == 7
        @test m.steps[1].taken_up == m.nodes[2].bound
        @test [(e.from, e.to) for e in m.steps] ==
              [(1, 2), (1, 4), (1, 6), (2, 3), (3, 5), (4, 3), (6, 5)]
        @test [_testhelper_names(e.taken_up) for e in m.steps] ==
              [[:A], [:B], [:Q], [:B], Symbol[], [:A], [:P]]
        @test all(isempty(e.given_off) for e in m.steps)
        @test [e.rapid_equilibrium for e in m.steps] ==
              [true, true, false, true, false, true, false]
        @test [M.step_kind(m, k) for k in eachindex(m.steps)] ==
              [:binding, :binding, :binding, :binding, :isomerization, :binding, :binding]
    end

    # Kinetic groups are numbered by their first step in step order.
    @testset "groups" begin
        @test length(M.groups(m)) == 5
        @test [e.group for e in m.steps] == [1, 2, 3, 2, 4, 1, 5]
        @test M.groups(m) == [[1, 6], [2, 4], [3], [5], [7]]
    end

    @testset "step_kind and catalytic allocate nothing" begin
        _testhelper_test_graph_reads_allocate_nothing(m, 4)
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
        @test length(m.nodes) == 3
        @test [_testhelper_names(n.bound) for n in m.nodes] ==
              [Symbol[], [:P], [:S]]
        @test [n.conformation for n in m.nodes] == [:E, :E, :E]
    end

    @testset "steps" begin
        @test length(m.steps) == 3
        @test [(e.from, e.to) for e in m.steps] == [(1, 2), (1, 3), (3, 2)]
        @test [_testhelper_names(e.taken_up) for e in m.steps] ==
              [[:P], [:S], Symbol[]]
        @test all(isempty(e.given_off) for e in m.steps)
        @test [e.rapid_equilibrium for e in m.steps] == [true, true, false]
        @test [M.step_kind(m, k) for k in eachindex(m.steps)] ==
              [:binding, :binding, :isomerization]
    end

    @testset "groups" begin
        @test length(M.groups(m)) == 3
        @test [e.group for e in m.steps] == [1, 2, 3]
        @test M.groups(m) == [[1], [2], [3]]
    end

    @testset "allosteric parts" begin
        # Each tag follows its group to the canonical group number.
        @test [am.mwc_tags[g] for g in 1:3] == [:EqualAI, :NonequalAI, :OnlyA]
        @test am.catalytic_multiplicity == 2
        @test am.regulatory_sites == [M.RegulatorySite([:R], 2, [:OnlyI])]
    end

    @testset "step_kind and catalytic allocate nothing" begin
        _testhelper_test_graph_reads_allocate_nothing(m, 2)
        @test _testhelper_allocated(M.catalytic, am) == 0
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
    @test [n.conformation for n in m.nodes] == [:E, :F, :F]
    @test [(e.from, e.to) for e in m.steps] == [(1, 3), (2, 1), (2, 3)]
    # A fused binding takes up one metabolite and gives off none, into a form that differs
    # from its source by more than that metabolite.
    @test [M.step_kind(m, k) for k in eachindex(m.steps)] ==
          [:fused_binding, :isomerization, :binding]
    @test [_testhelper_names(e.taken_up) for e in m.steps] ==
          [[:S], Symbol[], [:P]]
    @test _testhelper_allocated(M.step_kind, m, 1) == 0

    # A Theorell–Chance step takes up B and gives off P: chemistry, no binding.
    tc = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end)
    # Canonical form numbers the forms E, E(A), E(Q) and stores the steps E + A → E(A),
    # E + Q → E(Q), E(A) + B → E(Q) + P.
    @test [(e.from, e.to) for e in tc.steps] == [(1, 2), (1, 3), (2, 3)]
    @test [M.step_kind(tc, k) for k in eachindex(tc.steps)] ==
          [:binding, :binding, :chemistry]

    # A merged complex entered by the last substrate: a fused binding; the isomerization
    # and the plain binding of the product release keep their kinds.
    fused = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(P, Q)
            E(P, Q) ⇌ E(Q) + P
            E(Q) ⇌ E + Q
        end
    end)
    @test sort([M.step_kind(fused, k) for k in eachindex(fused.steps)]) ==
          [:binding, :binding, :binding, :fused_binding]
end

@testset "graph accessors: a ping-pong's residual" begin
    M = ER.Mechanisms
    m = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(; residual = A - P) + P
            E(; residual = A - P) + B <--> E(B; residual = A - P)
            E(B; residual = A - P) <--> E + Q
        end
    end)
    # Canonical form numbers the forms E, E(; residual = A - P), E(A),
    # E(B; residual = A - P): forms of one label sort by their bound metabolites, the
    # empty list first.
    @test [_testhelper_names(n.bound) for n in m.nodes] ==
          [Symbol[], Symbol[], [:A], [:B]]
    # The residual: the substrates a form adds and the products it subtracts.
    @test [(n.added, n.subtracted) for n in m.nodes] ==
          [(Symbol[], Symbol[]), ([:A], [:P]), (Symbol[], Symbol[]), ([:A], [:P])]
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
    ESS = m.nodes[4].bound
    # A reactant and its modifier copy are two entries, read under one name, as they share
    # a concentration.
    @test ESS == [M.Met(:S, false), M.Met(:S, true)]
    @test _testhelper_names(ESS) == [:S, :S]
    @test (m.steps[4].from, m.steps[4].to) == (3, 4)
    @test m.steps[4].taken_up == [M.Met(:S, true)]
    @test M.step_kind(m, 4) === :binding
    # The modifier copy is an entry of its own, so it adds to the bound S.
    @test count(x -> x.name === :S, [m.nodes[3].bound; m.steps[4].taken_up]) == 2
    @test isempty(intersect(m.nodes[3].bound, m.steps[4].taken_up))
    # The copies held as modifiers: the modifier copy of S, not S itself.
    @test filter(x -> x.modifier, ESS) == [M.Met(:S, true)]
    @test isempty(filter(x -> x.modifier, m.nodes[3].bound))
    @test all(x -> x.modifier, m.steps[4].taken_up)
    @test isempty(filter(x -> x.modifier, m.steps[2].taken_up))

    # Two reads of one mechanism agree, and a form of another mechanism compares by content.
    other = ER.Mechanism(em)
    @test other.nodes[4].bound == ESS && hash(other.nodes[4].bound) == hash(ESS)
    @test other.nodes[3].bound != ESS
    uniuni = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
    @test uniuni.nodes[3].bound == m.nodes[3].bound
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
    @test [(e.from, e.to) for e in m.steps] == [(1, 2), (1, 4), (1, 5), (2, 3), (5, 4)]
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
    kept = filter(k -> on[k], collect(eachindex(m.steps)))
    incidence = zeros(length(kept), length(eachindex(m.nodes)))
    for (r, k) in enumerate(kept)
        a, b = m.steps[k].from, m.steps[k].to
        incidence[r, a] -= 1
        incidence[r, b] += 1
    end
    names = unique(reduce(vcat, [[_testhelper_names(e.taken_up);
                                  _testhelper_names(e.given_off)]
                                 for e in m.steps]))
    any(names) do x
        uptake = Int[count(y -> y.name === x, m.steps[k].taken_up) -
                     count(y -> y.name === x, m.steps[k].given_off) for k in kept]
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
        n = length(eachindex(m.steps))
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
        n = length(eachindex(m.steps))
        # Every mechanism turns the reaction over, and none at rapid equilibrium alone.
        @test (spec.name, M.has_turnover_cycle(m, trues(n))) == (spec.name, true)
        @test (spec.name, M.has_turnover_cycle(m, _testhelper_flags(m))) ==
              (spec.name, false)
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
                   for k in eachindex(aq.steps))
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
                   for k in eachindex(b_only.steps))
    @test !M.has_turnover_cycle(b_only, on)
    # A mask with fewer flags than steps is out of bounds.
    @test_throws BoundsError M.has_turnover_cycle(aq, trues(2))
end

@testset "segment_graph: a steady-state step inside a segment" begin
    M = ER.Mechanisms
    # Random-order bi-bi whose binding of B to free E alone is steady state: the other
    # three bindings put E, E(A), E(B) and E(A, B) in one rapid-equilibrium segment. The
    # guard rejects it (E(B) fills at zero A and B); the view reads any graph.
    m = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(B) + A ⇌ E(A, B)
            E + B <--> E(B)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(P, Q) <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end
    every = trues(length(eachindex(m.steps)))
    text(k) = _testhelper_edge_text(m, k)
    g = M.segment_graph(m, every)
    # Canonical form numbers the forms E, E(A), E(A, B), E(B), E(P, Q), E(Q).
    @test g.segments == [[1, 2, 3, 4], [5], [6]]
    # E + B <--> E(B) joins two forms of the first segment, so it is no step of the graph.
    @test text.(g.steps) == ["E + Q <--> EQ", "EAB <--> EPQ", "EQ + P <--> EPQ"]
    @test g.bottom == Int[]
end

@testset "step_text and form_text: steps and forms as the macros write them" begin
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
    # Canonical form numbers the forms E, E(A), E(A, B), E(B), E(P, Q), E(Q) and stores
    # each step physical-forward, a release as the binding it reverses.
    @test [M.form_text(m, x) for x in eachindex(m.nodes)] ==
          ["E", "E(A)", "E(A, B)", "E(B)", "E(P, Q)", "E(Q)"]
    @test [M.step_text(m, k) for k in eachindex(m.steps)] ==
          ["E + A ⇌ E(A)", "E + B ⇌ E(B)", "E + Q <--> E(Q)", "E(A) + B ⇌ E(A, B)",
           "E(A, B) <--> E(P, Q)", "E(B) + A ⇌ E(A, B)", "E(Q) + P <--> E(P, Q)"]
    # A residual and a modifier copy.
    pingpong = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(; residual = A - P) + P
            E(; residual = A - P) + B <--> E(B; residual = A - P)
            E(B; residual = A - P) <--> E + Q
        end
    end)
    @test [M.form_text(pingpong, x) for x in eachindex(pingpong.nodes)] ==
          ["E", "E(; residual = A - P)", "E(A)", "E(B; residual = A - P)"]
    inhibited = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: S
        steps: begin
            E + S ⇌ E(S)
            E(S) + S::Mod ⇌ E(S, S::Mod)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
    @test M.form_text(inhibited, 4) == "E(S, S::Mod)"
    @test M.step_text(inhibited, 4) == "E(S) + S::Mod ⇌ E(S, S::Mod)"
end

@testset "runs_toward_products: a ping-pong's half-reactions" begin
    M = ER.Mechanisms
    m = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(; residual = A - P) + P
            E(; residual = A - P) + B <--> E(B; residual = A - P)
            E(B; residual = A - P) <--> E + Q
        end
    end)
    # Each step is stored as its uptake: E + A → E(A), E + Q → E(B; residual = A - P),
    # E(; residual = A - P) + P → E(A), E(; residual = A - P) + B → E(B; residual = A - P).
    # The releases of P and Q run toward products against their stored direction; the
    # uptakes of A and B tie (one substrate on each side) and run as stored.
    @test [M.step_text(m, k) for k in eachindex(m.steps)] ==
          ["E + A <--> E(A)", "E + Q <--> E(B; residual = A - P)",
           "E(; residual = A - P) + P <--> E(A)",
           "E(; residual = A - P) + B <--> E(B; residual = A - P)"]
    @test [M.runs_toward_products(m, k) for k in eachindex(m.steps)] ==
          [true, false, false, true]
end

@testset "excludes and exclusion_components: an :OnlyA dead-end binding" begin
    M = ER.Mechanisms
    am = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: S
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)                    :: EqualAI
            E(S) + S::Mod ⇌ E(S, S::Mod)    :: OnlyA
            E(S) <--> E(P)                  :: EqualAI
            E(P) ⇌ E + P                    :: EqualAI
        end
    end)
    m = M.catalytic(am)
    # Canonical form numbers the forms E, E(P), E(S), E(S, S::Mod) and stores the steps
    # E + P → E(P), E + S → E(S), E(S) → E(P), E(S) + S::Mod → E(S, S::Mod). The inactive
    # conformation excludes the :OnlyA binding, which leaves E(S, S::Mod) alone.
    @test M.step_text(m, 4) == "E(S) + S::Mod ⇌ E(S, S::Mod)"
    @test [M.excludes(am, k, :active) for k in eachindex(m.steps)] ==
          [false, false, false, false]
    @test [M.excludes(am, k, :inactive) for k in eachindex(m.steps)] ==
          [false, false, false, true]
    @test M.exclusion_components(am, :active) == [1, 1, 1, 1]
    @test M.exclusion_components(am, :inactive) == [1, 1, 1, 4]
    # A plain mechanism excludes no step from either conformation.
    @test [M.excludes(m, k, c) for k in eachindex(m.steps) for c in (:active, :inactive)] ==
          falses(8)
    # The rule by tag: an :OnlyA step leaves the inactive conformation, an :OnlyI step the
    # active one; an :EqualAI or :NonequalAI step stays in both.
    @test [M.excludes(tag, c) for tag in (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)
           for c in (:active, :inactive)] ==
          [false, true, true, false, false, false, false, false]
    @test_throws ArgumentError M.excludes(:OnlyA, :both)
    @test_throws ArgumentError M.excludes(am, 1, :both)
    @test_throws ArgumentError M.excludes(m, 1, :both)
    @test_throws ArgumentError M.exclusion_components(am, :both)
end

@testset "constant_columns: the joint solve's columns and their names" begin
    M = ER.Mechanisms
    col(edge, slot, inactive) = (; edge, slot, inactive)
    # Uni-uni: E + P → E(P) and E + S → E(S) at rapid equilibrium, E(S) → E(P) steady
    # state, each its own group. The isomerization's columns come first (kind 3), its
    # reverse rate, which runs against the reaction, before its forward one; then the
    # bindings, the larger content rank first.
    m = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
    c = M.constant_columns(m)
    @test c.columns == [col(3, 2, false), col(3, 1, false), col(2, 1, false),
                        col(1, 1, false)]
    @test c.names == [:k_EP_to_ES, :k_ES_to_EP, :K_ES_to_E_S, :K_EP_to_E_P]
    @test c.in_law == trues(4)
    @test (c.flat, c.tags, c.law) == ([1, 2, 3], [:None, :None, :None], falses(3))

    # Allosteric uni-uni: E + P → E(P) :EqualAI, E + S → E(S) :NonequalAI and
    # E(S) → E(P) :OnlyA, one group each. The :NonequalAI and :OnlyA groups get inactive
    # columns; the inactive copies of the :OnlyA step lie outside the inactive law
    # (block 2) and come first, then the inactive binding (block 1), then the active
    # columns (block 0). An :OnlyA step's inactive copy takes its own name.
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
    a = M.constant_columns(am)
    @test a.columns == [col(3, 2, true), col(3, 1, true), col(2, 1, true),
                        col(3, 2, false), col(3, 1, false), col(2, 1, false),
                        col(1, 1, false)]
    @test a.names == [:k_I_EP_to_ES, :k_I_ES_to_EP, :K_I_ES_to_E_S, :k_A_EP_to_ES,
                      :k_A_ES_to_EP, :K_A_ES_to_E_S, :K_EP_to_E_P]
    @test a.in_law == BitVector([0, 0, 1, 1, 1, 1, 1])
    @test (a.flat, a.tags, a.law) ==
          ([1, 2, 3], [:EqualAI, :NonequalAI, :OnlyA], BitVector([1, 1, 0]))
end

@testset "regulator_constant: a regulator ligand's constant per conformation" begin
    M = ER.Mechanisms
    am = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R, X, Y, Z
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)      :: EqualAI
            E(S) <--> E(P)    :: NonequalAI
            E(P) ⇌ E + P      :: EqualAI
        end
        regulatory_site(2): R::OnlyI
        regulatory_site(2): X::NonequalAI
        regulatory_site(2): Y::EqualAI
        regulatory_site(2): Z::OnlyA
    end)
    # Each ligand with its own tag, in the active and then the inactive conformation: an
    # :OnlyI ligand is absent from the active one and an :OnlyA ligand from the inactive
    # one; an :EqualAI ligand has one bare name in both.
    @test [M.regulator_constant(am, lig, tag, inactive)
           for (lig, tag) in ((:R, :OnlyI), (:X, :NonequalAI), (:Y, :EqualAI), (:Z, :OnlyA))
           for inactive in (false, true)] ==
          [nothing, :K_I_Rreg, :K_A_Xreg, :K_I_Xreg, :K_Yreg, :K_Yreg, :K_A_Zreg, nothing]
end

@testset "segment_walk and has_unbalanced_cycle on plain arrays" begin
    M = ER.Mechanisms
    # Nodes 1-3 present, node 4 not. Transitions: 1 → 2 (+1), 2 → 1 (-1), 2 → 3 (+1),
    # 3 → 1 (-2), each taking up one name by `change`.
    present = BitVector([1, 1, 1, 0])
    out = [[1], [2, 3], [4], Int[]]
    to = [2, 1, 3, 1]
    w = M.segment_walk(present, out, to, [1 -1 1 -2])
    @test w.segment_of == [1, 1, 1, 0]
    @test w.segments == [[1, 2, 3]]
    @test w.offsets == [0 1 2 0]
    @test w.parent == [0, 1, 3, 0]
    @test w.disagreeing == 0
    # With 3 → 1 taking up -1, the cycle 1 → 2 → 3 → 1 takes up one copy net: transition 4
    # is the first whose offsets disagree.
    @test M.segment_walk(present, out, to, [1 -1 1 -1]).disagreeing == 4

    # The triangle 1 → 2 → 3 → 1 balances at weights 1, 1, -2 and not at 1, 1, -1; without
    # its third step it closes no cycle.
    from, to3 = [1, 2, 3], [2, 3, 1]
    buffers() = (zeros(Int, 3), zeros(Int, 3))
    @test !M.has_unbalanced_cycle(from, to3, [1, 1, -2], trues(3), buffers()...)
    @test M.has_unbalanced_cycle(from, to3, [1, 1, -1], trues(3), buffers()...)
    @test !M.has_unbalanced_cycle(from, to3, [1, 1, -1], BitVector([1, 1, 0]),
                                  buffers()...)
end
