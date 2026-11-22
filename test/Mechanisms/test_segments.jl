# ABOUTME: Tests of the segments and cycles of a step mask: `rapid_equilibrium_segments`,
# ABOUTME: `has_turnover_cycle` and `has_nonzero_weight_cycle`.

const ER = EnzymeRates

"""Whether a cycle of the steps of graph `m` that `on` keeps takes up some metabolite net:
for some metabolite, the uptake of the kept steps is no difference of form potentials,
which holds when appending the uptake column to their incidence matrix raises its rank."""
function _testhelper_net_uptake_cycle(m, on::BitVector)
    kept = filter(k -> on[k], collect(eachindex(m.steps)))
    incidence = zeros(length(kept), length(m.enzyme_forms))
    for (r, k) in enumerate(kept)
        a, b = m.steps[k].from, m.steps[k].to
        incidence[r, a] -= 1
        incidence[r, b] += 1
    end
    metabolites = unique(reduce(vcat, [[e.taken_up; e.given_off] for e in m.steps]))
    any(metabolites) do x
        uptake = Int[count(==(x), m.steps[k].taken_up) -
                     count(==(x), m.steps[k].given_off) for k in kept]
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
                E + S ⇌ ES
                ES <--> E + P
            end
        end), [["E + S ⇌ ES", "E + P <--> ES"]]),
        # Uni-uni through a second free form F.
        (@enzyme_mechanism(begin
            substrates: S
            products: P
            steps: begin
                E + S <--> FP
                FP ⇌ F + P
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
                (E + A ⇌ EA, EB + A ⇌ EAB)
                (E + B ⇌ EB, EA + B ⇌ EAB)
                EAB <--> EPQ
                EPQ <--> EQ + P
                EQ <--> E + Q
            end
        end), [["E + A ⇌ EA", "EA + B ⇌ EAB", "EAB <--> EPQ", "EQ + P <--> EPQ",
                "E + Q <--> EQ"],
               ["E + B ⇌ EB", "EB + A ⇌ EAB", "EAB <--> EPQ", "EQ + P <--> EPQ",
                "E + Q <--> EQ"]]),
        # Substrate inhibition: the dead-end binding of a second S lies on no cycle.
        (@enzyme_mechanism(begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ ES
                ES + S ⇌ ESSmod
                ES <--> EP
                EP ⇌ E + P
            end
        end), [["E + S ⇌ ES", "ES <--> EP", "E + P ⇌ EP"]]),
        # A Theorell–Chance step beside the merged complex EPQ: the triangle of the
        # two routes from EA to EQ takes up B and gives off P once each way.
        (@enzyme_mechanism(begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                EA + B <--> EPQ
                EPQ <--> EQ + P
                EQ <--> E + Q
                EA + B <--> EQ + P
            end
        end), [["E + A ⇌ EA", "EA + B <--> EPQ", "EQ + P <--> EPQ", "E + Q <--> EQ"],
               ["E + A ⇌ EA", "EA + B <--> EQ + P", "E + Q <--> EQ"]]),
    ]
    for (m, cycles) in cases
        n = length(eachindex(m.catalytic.steps))
        at = [[_testhelper_edge(m, text) for text in cycle] for cycle in cycles]
        for bits in 0:(2^n - 1)
            on = BitVector([isodd(bits >> (k - 1)) for k in 1:n])
            expected = any(ks -> all(k -> on[k], ks), at)
            @test (findall(on), M.has_turnover_cycle(m.catalytic, on)) ==
                  (findall(on), expected)
        end
    end
end

@testset "turnover cycle: MECHANISM_TEST_SPECS" begin
    M = ER.Mechanisms
    for spec in MECHANISM_TEST_SPECS
        mech = spec.mechanism
        m = mech.catalytic
        n = length(eachindex(m.steps))
        # Every mechanism turns the reaction over, and none at rapid equilibrium alone.
        @test (spec.name, M.has_turnover_cycle(m, trues(n))) == (spec.name, true)
        @test (spec.name, M.has_turnover_cycle(m, _testhelper_flags(m))) ==
              (spec.name, false)
        subsets = [[Int[]]; [[i] for i in 1:n]; [[i, j] for i in 1:n for j in (i + 1):n]]
        masks = [BitVector(k in s for k in 1:n) for s in subsets]
        _testhelper_is_allosteric(mech) &&
            push!(masks, M.steps_in_active_conformation(mech),
                  M.steps_in_inactive_conformation(mech))
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
    # the first cycle and not the second. The guard rejects both graphs.
    aq = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B <--> EAB
            EAB <--> EQ + P
            EQ <--> E + Q
            EA + Q <--> EAQ
            EAQ <--> E
        end
    end
    on = BitVector(k in [_testhelper_edge(aq, t)
                         for t in ("E + A <--> EA", "EA + Q <--> EAQ", "EAQ <--> E")]
                   for k in eachindex(aq.catalytic.steps))
    @test M.has_turnover_cycle(aq.catalytic, on)
    b_only = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B <--> EAB
            EAB <--> EQ + P
            EQ <--> E + Q
            E + B <--> EB
            EB <--> E
        end
    end
    on = BitVector(k in [_testhelper_edge(b_only, t)
                         for t in ("E + B <--> EB", "E <--> EB")]
                   for k in eachindex(b_only.catalytic.steps))
    @test !M.has_turnover_cycle(b_only.catalytic, on)
    # A mask with fewer flags than steps is out of bounds.
    @test_throws BoundsError M.has_turnover_cycle(aq.catalytic, trues(2))
end

@testset "segments: a steady-state step inside a segment is no step of the graph" begin
    M = ER.Mechanisms
    # Random-order bi-bi whose binding of B to free E alone is steady state: the other
    # three bindings put E, EA, EB and EAB in one rapid-equilibrium segment. The
    # guard rejects it (EB fills at zero A and B); the view reads any graph.
    m = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EB + A ⇌ EAB
            E + B <--> EB
            EA + B ⇌ EAB
            EAB <--> EPQ
            EPQ <--> EQ + P
            EQ <--> E + Q
        end
    end
    every = trues(length(eachindex(m.catalytic.steps)))
    text(k) = _testhelper_edge_text(m, k)
    g = M.rapid_equilibrium_segments(m.catalytic, every)
    # Canonical form numbers the forms E, EA, EAB, EB, EPQ, EQ.
    @test g.segments == [[1, 2, 3, 4], [5], [6]]
    # E + B <--> EB joins two forms of the first segment, so it is no step of the graph.
    @test text.(g.steps) == ["E + Q <--> EQ", "EAB <--> EPQ", "EQ + P <--> EPQ"]
end

@testset "has_nonzero_weight_cycle on plain arrays" begin
    M = ER.Mechanisms
    # The triangle 1 → 2 → 3 → 1 balances at weights 1, 1, -2 and not at 1, 1, -1; without
    # its third step it closes no cycle.
    from, to3 = [1, 2, 3], [2, 3, 1]
    buffers() = (zeros(Int, 3), zeros(Int, 3))
    @test !M.has_nonzero_weight_cycle(from, to3, [1, 1, -2], trues(3), buffers()...)
    @test M.has_nonzero_weight_cycle(from, to3, [1, 1, -1], trues(3), buffers()...)
    @test !M.has_nonzero_weight_cycle(from, to3, [1, 1, -1], BitVector([1, 1, 0]),
                                      buffers()...)
end

@testset "a sub-graph reads the forms and segments of the edges a mask keeps" begin
    # Random release of P and Q at rapid equilibrium between steady-state releases: the
    # segment {EP, EPQ, EQ} has every form holding a product.
    m = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES
            ES <--> EPQ
            EPQ ⇌ EQ + P
            EPQ ⇌ EP + Q
            EP <--> E + P
            EQ <--> E + Q
        end
    end
    names(forms) = String.(m.catalytic.enzyme_forms[forms])
    @test names(eachindex(m.catalytic.enzyme_forms)) == ["E", "EP", "EPQ", "EQ", "ES"]
    # One edge per group, so the edges group by group are the edges in order.
    @test [_testhelper_edge_text(m, k) for k in ER.Derivation._flat_edges(m.catalytic)] ==
          ["E + P <--> EP", "E + Q <--> EQ", "E + S ⇌ ES", "EP + Q ⇌ EPQ", "EQ + P ⇌ EPQ",
           "ES <--> EPQ"]
    @test _testhelper_flags(m) == BitVector([0, 0, 1, 1, 1, 0])
    every = trues(length(m.catalytic.steps))
    @test ER.Mechanisms._present_species(m.catalytic, every) == trues(5)
    # Each segment breadth-first from its first form; offsets per metabolite (rows S, P,
    # Q, the metabolites in stored order) above the segment's lowest.
    s = ER.Mechanisms.rapid_equilibrium_segments(m.catalytic, every, _testhelper_flags(m))
    @test [names(seg) for seg in s.segments] == [["E", "ES"], ["EP", "EPQ", "EQ"]]
    @test s.segment_of == [1, 2, 2, 2, 1]
    @test s.offsets == [0 0 0 0 1; 0 1 1 0 0; 0 0 1 1 0]
    # The steady-state steps that join two segments, which are all of them here.
    g = ER.Mechanisms.rapid_equilibrium_segments(m.catalytic, every)
    @test g.steps == [1, 2, 6]
    @test (g.segment_of, g.segments, g.offsets) == (s.segment_of, s.segments, s.offsets)
    # Without the release of P from EPQ, EQ stands alone.
    no_p = BitVector([1, 1, 1, 1, 0, 1])
    @test ER.Mechanisms._present_species(m.catalytic, no_p) == trues(5)
    s = ER.Mechanisms.rapid_equilibrium_segments(m.catalytic, no_p, _testhelper_flags(m))
    @test [names(seg) for seg in s.segments] == [["E", "ES"], ["EP", "EPQ"], ["EQ"]]
    @test ER.Mechanisms.rapid_equilibrium_segments(m.catalytic, no_p).steps == [1, 2, 6]
    # A mask keeps only the forms its edges touch; a form it leaves is in no segment.
    bindings = BitVector([1, 1, 0, 0, 0, 0])
    @test ER.Mechanisms._present_species(m.catalytic, bindings) ==
          BitVector([1, 1, 0, 1, 0])
    s = ER.Mechanisms.rapid_equilibrium_segments(m.catalytic, bindings,
                                                 _testhelper_flags(m))
    @test [names(seg) for seg in s.segments] == [["E"], ["EP"], ["EQ"]]
    @test s.segment_of == [1, 2, 0, 3, 0]
    @test ER.Mechanisms.rapid_equilibrium_segments(m.catalytic, bindings).steps == [1, 2]
    # `re` picks the edges that join segments: with none, each present form stands alone.
    none = ER.Mechanisms.rapid_equilibrium_segments(m.catalytic, every, falses(6))
    @test length(none.segments) == 5
    @test none.steps == 1:6

    # The inactive conformation keeps the edges outside :OnlyA groups between forms that
    # such edges reach from free E: here EA and EAI lie behind the :OnlyA
    # binding of A.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I
        oligomeric_state: 1
        steps: begin
            E + A <--> EA                                   :: OnlyA
            EA <--> EP                                    :: OnlyA
            E + P ⇌ EP                                      :: EqualAI
            (EA + I ⇌ EAI, EP + I ⇌ EIP)          :: NonequalAI
        end
    end
    inactive = ER.Mechanisms.steps_in_inactive_conformation(am)
    @test [_testhelper_edge_text(am, k) for k in findall(inactive)] ==
          ["E + P ⇌ EP", "EP + I ⇌ EIP"]
    @test am.catalytic.enzyme_forms[
              findall(ER.Mechanisms._present_species(am.catalytic, inactive))] ==
          [:E, :EIP, :EP]
end
