# ABOUTME: Tests of what a mechanism's stored graph holds (nodes, steps, kinetic groups,
# ABOUTME: allosteric parts), `step_kind`, `groups`, the stored step facts, step text.

const ER = EnzymeRates

"""The names of the metabolites in `list`, in order."""
_testhelper_names(list) = Symbol[x.name for x in list]

@testset "graph accessors: random-order bi-bi" begin
    M = ER.Mechanisms
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
            E(A, B) <--> E(P, Q)
            E(P, Q) <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end

    # Form numbers follow canonical form: free E, then the forms by their bound
    # metabolites read as sorted lists, so E(A, B) comes before E(B).
    @testset "species" begin
        @test length(m.catalytic.nodes) == 6
        @test [_testhelper_names(n.bound) for n in m.catalytic.nodes] ==
              [Symbol[], [:A], [:A, :B], [:B], [:P, :Q], [:Q]]
        @test [n.conformation for n in m.catalytic.nodes] == fill(:E, 6)
    end

    # Step numbers and directions follow canonical form: each step is stored the way the
    # graph gives it, so a product release is stored as the binding it reverses, and the
    # steps are sorted by their stored ends.
    @testset "steps" begin
        @test length(m.catalytic.steps) == 7
        @test m.catalytic.steps[1].taken_up == m.catalytic.nodes[2].bound
        @test [(e.from, e.to) for e in m.catalytic.steps] ==
              [(1, 2), (1, 4), (1, 6), (2, 3), (3, 5), (4, 3), (6, 5)]
        @test [_testhelper_names(e.taken_up) for e in m.catalytic.steps] ==
              [[:A], [:B], [:Q], [:B], Symbol[], [:A], [:P]]
        @test all(isempty(e.given_off) for e in m.catalytic.steps)
        @test [e.rapid_equilibrium for e in m.catalytic.steps] ==
              [true, true, false, true, false, true, false]
        @test [m.catalytic.step_kinds[k] for k in eachindex(m.catalytic.steps)] ==
              [:binding, :binding, :binding, :binding, :isomerization, :binding, :binding]
    end

    # Kinetic groups are numbered by their first step in step order.
    @testset "groups" begin
        @test length(m.catalytic.groups) == 5
        @test [e.group for e in m.catalytic.steps] == [1, 2, 3, 2, 4, 1, 5]
        @test m.catalytic.groups == [[1, 6], [2, 4], [3], [5], [7]]
    end
end

@testset "graph accessors: allosteric uni-uni" begin
    M = ER.Mechanisms
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ E(S)      :: NonequalAI
            E(S) <--> E(P)    :: OnlyA
            E(P) ⇌ E + P      :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    m = am.catalytic
    @test m isa ER.Mechanisms.CatalyticScheme

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
        @test [m.step_kinds[k] for k in eachindex(m.steps)] ==
              [:binding, :binding, :isomerization]
    end

    @testset "groups" begin
        @test length(m.groups) == 3
        @test [e.group for e in m.steps] == [1, 2, 3]
        @test m.groups == [[1], [2], [3]]
    end

    @testset "allosteric parts" begin
        # Each tag follows its group to the canonical group number.
        @test [am.allosteric.allosteric_states[g] for g in 1:3] ==
              [:EqualAI, :NonequalAI, :OnlyA]
        @test am.allosteric.oligomeric_state == 2
        @test am.allosteric.regulatory_sites == [M.RegulatorySite([:R], 2, [:OnlyI])]
    end
end

@testset "graph accessors: conformations and chemistry" begin
    M = ER.Mechanisms
    # Uni-uni through a second free form F: S binds straight into the product-bound form.
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> F(P)
            F(P) ⇌ F + P
            E <--> F
        end
    end
    # Canonical form numbers the forms E, F, F(P). It stores the isomerization from F,
    # where only products enter or leave, to E, where only substrates do: against its
    # writing.
    @test [n.conformation for n in m.catalytic.nodes] == [:E, :F, :F]
    @test [(e.from, e.to) for e in m.catalytic.steps] == [(1, 3), (2, 1), (2, 3)]
    # A fused binding takes up one metabolite and gives off none, into a form that differs
    # from its source by more than that metabolite.
    @test [m.catalytic.step_kinds[k] for k in eachindex(m.catalytic.steps)] ==
          [:fused_binding, :isomerization, :binding]
    @test [_testhelper_names(e.taken_up) for e in m.catalytic.steps] ==
          [[:S], Symbol[], [:P]]

    # A Theorell–Chance step takes up B and gives off P: chemistry, no binding.
    tc = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end
    # Canonical form numbers the forms E, E(A), E(Q) and stores the steps E + A → E(A),
    # E + Q → E(Q), E(A) + B → E(Q) + P.
    @test [(e.from, e.to) for e in tc.catalytic.steps] == [(1, 2), (1, 3), (2, 3)]
    @test [tc.catalytic.step_kinds[k] for k in eachindex(tc.catalytic.steps)] ==
          [:binding, :binding, :chemistry]

    # A merged complex entered by the last substrate: a fused binding; the isomerization
    # and the plain binding of the product release keep their kinds.
    fused = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(P, Q)
            E(P, Q) ⇌ E(Q) + P
            E(Q) ⇌ E + Q
        end
    end
    @test sort([fused.catalytic.step_kinds[k] for k in eachindex(fused.catalytic.steps)]) ==
          [:binding, :binding, :binding, :fused_binding]
end

@testset "graph accessors: a ping-pong's residual" begin
    M = ER.Mechanisms
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(; residual = A - P) + P
            E(; residual = A - P) + B <--> E(B; residual = A - P)
            E(B; residual = A - P) <--> E + Q
        end
    end
    # Canonical form numbers the forms E, E(; residual = A - P), E(A),
    # E(B; residual = A - P): forms of one label sort by their bound metabolites, the
    # empty list first.
    @test [_testhelper_names(n.bound) for n in m.catalytic.nodes] ==
          [Symbol[], Symbol[], [:A], [:B]]
    # The residual: the substrates a form adds and the products it subtracts.
    @test [(n.added, n.subtracted) for n in m.catalytic.nodes] ==
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
    m = em
    # Canonical form numbers the forms E, E(P), E(S), E(S, S::Mod) and sorts the steps
    # E + P → E(P), E + S → E(S), E(S) → E(P), E(S) + S::Mod → E(S, S::Mod).
    ESS = m.catalytic.nodes[4].bound
    # A reactant and its modifier copy are two entries, read under one name, as they share
    # a concentration.
    @test ESS == [M.Metabolite(:S, false), M.Metabolite(:S, true)]
    @test _testhelper_names(ESS) == [:S, :S]
    @test (m.catalytic.steps[4].from, m.catalytic.steps[4].to) == (3, 4)
    @test m.catalytic.steps[4].taken_up == [M.Metabolite(:S, true)]
    @test m.catalytic.step_kinds[4] === :binding
    # The modifier copy is an entry of its own, so it adds to the bound S.
    @test count(x -> x.name === :S,
                [m.catalytic.nodes[3].bound; m.catalytic.steps[4].taken_up]) == 2
    @test isempty(intersect(m.catalytic.nodes[3].bound, m.catalytic.steps[4].taken_up))
    # The copies held as modifiers: the modifier copy of S, not S itself.
    @test filter(x -> x.modifier, ESS) == [M.Metabolite(:S, true)]
    @test isempty(filter(x -> x.modifier, m.catalytic.nodes[3].bound))
    @test all(x -> x.modifier, m.catalytic.steps[4].taken_up)
    @test isempty(filter(x -> x.modifier, m.catalytic.steps[2].taken_up))

    # A copy of one mechanism agrees, and a form of another mechanism compares by content.
    other = deepcopy(em)
    @test other.catalytic.nodes[4].bound == ESS &&
          hash(other.catalytic.nodes[4].bound) == hash(ESS)
    @test other.catalytic.nodes[3].bound != ESS
    uniuni = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    @test uniuni.catalytic.nodes[3].bound == m.catalytic.nodes[3].bound
end

@testset "step_text and form_text: steps and forms as the macros write them" begin
    M = ER.Mechanisms
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
            E(A, B) <--> E(P, Q)
            E(P, Q) <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end
    # Canonical form numbers the forms E, E(A), E(A, B), E(B), E(P, Q), E(Q) and stores
    # each step the way the graph gives it, a release as the binding it reverses.
    @test [M.form_text(m.catalytic, x) for x in eachindex(m.catalytic.nodes)] ==
          ["E", "E(A)", "E(A, B)", "E(B)", "E(P, Q)", "E(Q)"]
    @test [M.step_text(m.catalytic, k) for k in eachindex(m.catalytic.steps)] ==
          ["E + A ⇌ E(A)", "E + B ⇌ E(B)", "E + Q <--> E(Q)", "E(A) + B ⇌ E(A, B)",
           "E(A, B) <--> E(P, Q)", "E(B) + A ⇌ E(A, B)", "E(Q) + P <--> E(P, Q)"]
    # A residual and a modifier copy.
    pingpong = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(; residual = A - P) + P
            E(; residual = A - P) + B <--> E(B; residual = A - P)
            E(B; residual = A - P) <--> E + Q
        end
    end
    @test [M.form_text(pingpong.catalytic, x)
           for x in eachindex(pingpong.catalytic.nodes)] ==
          ["E", "E(; residual = A - P)", "E(A)", "E(B; residual = A - P)"]
    inhibited = @enzyme_mechanism begin
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
    @test M.form_text(inhibited.catalytic, 4) == "E(S, S::Mod)"
    @test M.step_text(inhibited.catalytic, 4) == "E(S) + S::Mod ⇌ E(S, S::Mod)"
end

@testset "toward_products: a ping-pong's half-reactions" begin
    M = ER.Mechanisms
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(; residual = A - P) + P
            E(; residual = A - P) + B <--> E(B; residual = A - P)
            E(B; residual = A - P) <--> E + Q
        end
    end
    # Each step is stored as its uptake: E + A → E(A), E + Q → E(B; residual = A - P),
    # E(; residual = A - P) + P → E(A), E(; residual = A - P) + B → E(B; residual = A - P).
    # The releases of P and Q run toward products against their stored direction; the
    # uptakes of A and B tie (one substrate on each side) and run as stored.
    @test [M.step_text(m.catalytic, k) for k in eachindex(m.catalytic.steps)] ==
          ["E + A <--> E(A)", "E + Q <--> E(B; residual = A - P)",
           "E(; residual = A - P) + P <--> E(A)",
           "E(; residual = A - P) + B <--> E(B; residual = A - P)"]
    @test [m.catalytic.toward_products[k] for k in eachindex(m.catalytic.steps)] ==
          [true, false, false, true]
end
