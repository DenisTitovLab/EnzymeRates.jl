# ABOUTME: Tests of the segments and cycles of a step mask: `re_segments`, `segment_graph`,
# ABOUTME: `segment_walk`, `has_turnover_cycle` and `has_unbalanced_cycle`.

const ER = EnzymeRates

"""The names of the metabolites in `list`, in order."""
_testhelper_met_names(list) = Symbol[x.name for x in list]

"""Whether a cycle of the steps of graph `m` that `on` keeps takes up some metabolite name
net: for some name, the uptake of the kept steps is no difference of form potentials,
which holds when appending the uptake column to their incidence matrix raises its rank."""
function _testhelper_net_uptake_cycle(m, on::BitVector)
    kept = filter(k -> on[k], collect(eachindex(m.steps)))
    incidence = zeros(length(kept), length(eachindex(m.nodes)))
    for (r, k) in enumerate(kept)
        a, b = m.steps[k].from, m.steps[k].to
        incidence[r, a] -= 1
        incidence[r, b] += 1
    end
    names = unique(reduce(vcat, [[_testhelper_met_names(e.taken_up);
                                  _testhelper_met_names(e.given_off)]
                                 for e in m.steps]))
    any(names) do x
        uptake = Int[count(y -> y.name === x, m.steps[k].taken_up) -
                     count(y -> y.name === x, m.steps[k].given_off) for k in kept]
        rank([incidence uptake]) > rank(incidence)
    end
end

@testset "turnover cycle: every subset of steps of five mechanisms" begin
    M = ER.Mechanisms
    D = ER.Derivation
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
            @test (findall(on), D.has_turnover_cycle(m, on)) == (findall(on), expected)
        end
    end
end

@testset "turnover cycle: MECHANISM_TEST_SPECS" begin
    M = ER.Mechanisms
    D = ER.Derivation
    for spec in MECHANISM_TEST_SPECS
        mech = _testhelper_lift(spec.mechanism)
        m = M.catalytic(mech)
        n = length(eachindex(m.steps))
        # Every mechanism turns the reaction over, and none at rapid equilibrium alone.
        @test (spec.name, D.has_turnover_cycle(m, trues(n))) == (spec.name, true)
        @test (spec.name, D.has_turnover_cycle(m, _testhelper_flags(m))) ==
              (spec.name, false)
        subsets = [[Int[]]; [[i] for i in 1:n]; [[i, j] for i in 1:n for j in (i + 1):n]]
        masks = [BitVector(k in s for k in 1:n) for s in subsets]
        mech isa ER.AllostericMechanism &&
            push!(masks, mech.active_steps,
                  mech.inactive_steps)
        for on in masks
            @test (spec.name, findall(on), D.has_turnover_cycle(m, on)) ==
                  (spec.name, findall(on), _testhelper_net_uptake_cycle(m, on))
        end
    end
end

@testset "turnover cycle: malformed input" begin
    M = ER.Mechanisms
    D = ER.Derivation
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
    @test D.has_turnover_cycle(aq, on)
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
    @test !D.has_turnover_cycle(b_only, on)
    # A mask with fewer flags than steps is out of bounds.
    @test_throws BoundsError D.has_turnover_cycle(aq, trues(2))
end

@testset "segment_graph: a steady-state step inside a segment" begin
    M = ER.Mechanisms
    D = ER.Derivation
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
    g = D.segment_graph(m, every)
    # Canonical form numbers the forms E, E(A), E(A, B), E(B), E(P, Q), E(Q).
    @test g.segments == [[1, 2, 3, 4], [5], [6]]
    # E + B <--> E(B) joins two forms of the first segment, so it is no step of the graph.
    @test text.(g.steps) == ["E + Q <--> EQ", "EAB <--> EPQ", "EQ + P <--> EPQ"]
    @test g.bottom == Int[]
end

@testset "segment_walk and has_unbalanced_cycle on plain arrays" begin
    M = ER.Mechanisms
    D = ER.Derivation
    # Nodes 1-3 present, node 4 not. Transitions: 1 → 2 (+1), 2 → 1 (-1), 2 → 3 (+1),
    # 3 → 1 (-2), each taking up one name by `change`.
    present = BitVector([1, 1, 1, 0])
    out = [[1], [2, 3], [4], Int[]]
    to = [2, 1, 3, 1]
    w = D.segment_walk(present, out, to, [1 -1 1 -2])
    @test w.segment_of == [1, 1, 1, 0]
    @test w.segments == [[1, 2, 3]]
    @test w.offsets == [0 1 2 0]
    @test w.parent == [0, 1, 3, 0]
    @test w.disagreeing == 0
    # With 3 → 1 taking up -1, the cycle 1 → 2 → 3 → 1 takes up one copy net: transition 4
    # is the first whose offsets disagree.
    @test D.segment_walk(present, out, to, [1 -1 1 -1]).disagreeing == 4

    # The triangle 1 → 2 → 3 → 1 balances at weights 1, 1, -2 and not at 1, 1, -1; without
    # its third step it closes no cycle.
    from, to3 = [1, 2, 3], [2, 3, 1]
    buffers() = (zeros(Int, 3), zeros(Int, 3))
    @test !D.has_unbalanced_cycle(from, to3, [1, 1, -2], trues(3), buffers()...)
    @test D.has_unbalanced_cycle(from, to3, [1, 1, -1], trues(3), buffers()...)
    @test !D.has_unbalanced_cycle(from, to3, [1, 1, -1], BitVector([1, 1, 0]),
                                  buffers()...)
end

@testset "a sub-graph reads the forms, segments and bottom of the edges a mask keeps" begin
    # Random release of P and Q at rapid equilibrium between steady-state releases: the
    # segment {E(P), E(P, Q), E(Q)} has every form holding a product.
    m = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P, Q)
            E(P, Q) ⇌ E(Q) + P
            E(P, Q) ⇌ E(P) + Q
            E(P) <--> E + P
            E(Q) <--> E + Q
        end
    end
    names(forms) = [ER.Mechanisms._form_name(m, f) for f in forms]
    @test names(eachindex(m.nodes)) == ["E", "EP", "EPQ", "EQ", "ES"]
    # One edge per group, so the edges group by group are the edges in order.
    @test [_testhelper_edge_text(m, k) for k in ER.Derivation._flat_edges(m)] ==
          ["E + P <--> EP", "E + Q <--> EQ", "E + S ⇌ ES", "EP + Q ⇌ EPQ", "EQ + P ⇌ EPQ",
           "ES <--> EPQ"]
    @test _testhelper_flags(m) == BitVector([0, 0, 1, 1, 1, 0])
    @test [ER.Mechanisms._is_free(f) for f in m.nodes] == [true, false, false, false, false]
    every = trues(length(m.steps))
    @test ER.Derivation.present_species(m, every) == trues(5)
    # Each segment breadth-first from its first form; offsets per metabolite name (rows
    # P, Q, S) above the segment's lowest.
    s = ER.Derivation.re_segments(m, every, _testhelper_flags(m))
    @test [names(seg) for seg in s.segments] == [["E", "ES"], ["EP", "EPQ", "EQ"]]
    @test s.segment_of == [1, 2, 2, 2, 1]
    @test s.offsets == [0 1 1 0 0; 0 0 1 1 0; 0 0 0 0 1]
    # Every form of the second segment holds a product, so it has no bottom form.
    @test names(s.bottom) == ["EP", "EPQ", "EQ"]
    # The steady-state steps that join two segments, which are all of them here.
    g = ER.Derivation.segment_graph(m, every)
    @test g.steps == [1, 2, 6]
    @test (g.segment_of, g.segments, g.offsets, g.bottom) ==
          (s.segment_of, s.segments, s.offsets, s.bottom)
    # Without the release of P from E(P, Q), E(Q) stands alone and E(P) is the bottom of
    # {E(P), E(P, Q)}.
    no_p = BitVector([1, 1, 1, 1, 0, 1])
    @test ER.Derivation.present_species(m, no_p) == trues(5)
    s = ER.Derivation.re_segments(m, no_p, _testhelper_flags(m))
    @test [names(seg) for seg in s.segments] == [["E", "ES"], ["EP", "EPQ"], ["EQ"]]
    @test isempty(s.bottom)
    @test ER.Derivation.segment_graph(m, no_p).steps == [1, 2, 6]
    # A mask keeps only the forms its edges touch; a form it leaves is in no segment.
    bindings = BitVector([1, 1, 0, 0, 0, 0])
    @test ER.Derivation.present_species(m, bindings) == BitVector([1, 1, 0, 1, 0])
    s = ER.Derivation.re_segments(m, bindings, _testhelper_flags(m))
    @test [names(seg) for seg in s.segments] == [["E"], ["EP"], ["EQ"]]
    @test s.segment_of == [1, 2, 0, 3, 0]
    @test ER.Derivation.segment_graph(m, bindings).steps == [1, 2]
    # `re` picks the edges that join segments: with none, each present form stands alone.
    @test length(ER.Derivation.re_segments(m, every, falses(6)).segments) == 5

    # The inactive conformation keeps the edges outside :OnlyA groups between forms that
    # such edges reach from free E: here E(A) and E(A, I) lie behind the :OnlyA
    # binding of A.
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A
        products: P
        modifiers: I
        catalytic_multiplicity: 1
        steps: begin
            E + A <--> E(A)                                   :: OnlyA
            E(A) <--> E(P)                                    :: OnlyA
            E + P ⇌ E(P)                                      :: EqualAI
            (E(A) + I ⇌ E(A, I), E(P) + I ⇌ E(I, P))          :: NonequalAI
        end
    end
    inactive = am.inactive_steps
    @test [_testhelper_edge_text(am, k) for k in findall(inactive)] ==
          ["E + P ⇌ EP", "EP + I ⇌ EIP"]
    @test [ER.Mechanisms._form_name(am.catalytic, f)
           for f in findall(ER.Derivation.present_species(am.catalytic, inactive))] ==
          ["E", "EIP", "EP"]
end
