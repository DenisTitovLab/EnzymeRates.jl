# ABOUTME: Tests of the written form of a mechanism and of the writers that rebuild one with
# ABOUTME: new flags, groups or tags: exact values, each expected mechanism written inline.

@testset "written form: forms, steps and mechanisms" begin
    M = EnzymeRates.Mechanisms
    # Substrate inhibition: a second S binds E(S) as a modifier.
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
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
    # Canonical form numbers the forms E, E(P), E(S), E(S, S::Mod) and sorts the steps
    # E + P → E(P), E + S → E(S), E(S) → E(P), E(S) + S::Mod → E(S, S::Mod), each its own
    # kinetic group.
    S, Smod, P, I = M.Met(:S, false), M.Met(:S, true), M.Met(:P, false), M.Met(:I, true)
    E, ES, EP = M.RawForm(:E, M.Met[]), M.RawForm(:E, [S]), M.RawForm(:E, [P])
    ESS = M.RawForm(:E, [S, Smod])
    @testset "forms" begin
        @test M.RawForm(m, 1) == E
        @test M.RawForm(m, 4) == ESS
        @test M.RawForm(m, 3, I) == M.RawForm(:E, [S, I])
        @test M.RawForm(m, 1, I, I) == M.RawForm(:E, [I, I])
    end
    @testset "steps" begin
        @test M.RawStep(m, 1) == M.RawStep(E, EP, [P], M.Met[], true)
        @test M.RawStep(m, 3) == M.RawStep(ES, EP, M.Met[], M.Met[], false)
        @test M.RawStep(m, 4) == M.RawStep(ES, ESS, [Smod], M.Met[], true)
        # With extra metabolites, both ends hold them: the step's mirror.
        @test M.RawStep(m, 3, I) == M.RawStep(M.RawForm(:E, [S, I]),
                                              M.RawForm(:E, [P, I]), M.Met[], M.Met[],
                                              false)
    end
    @testset "mechanism" begin
        raw = M.RawMechanism(m)
        @test raw == M.RawMechanism([:S], [:P], [:S],
                                    [M.RawStep(E, EP, [P], M.Met[], true),
                                     M.RawStep(E, ES, [S], M.Met[], true),
                                     M.RawStep(ES, EP, M.Met[], M.Met[], false),
                                     M.RawStep(ES, ESS, [Smod], M.Met[], true)],
                                    [1, 2, 3, 4])
        @test EnzymeRates.Mechanism(raw) == m
    end
    @testset "a covalent residual" begin
        pingpong = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(A) <--> E(; residual = A - P) + P
                E(; residual = A - P) + B <--> E(B; residual = A - P)
                E(B; residual = A - P) <--> E + Q
            end
        end)
        # Form 2 is E(; residual = A - P).
        @test M.RawForm(pingpong, 2) == M.RawForm(:E, M.Met[], [:A], [:P])
        @test EnzymeRates.Mechanism(M.RawMechanism(pingpong)) == pingpong
    end
    @testset "allosteric mechanism" begin
        am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
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
        raw = M.RawAllostericMechanism(am)
        @test raw == M.RawAllostericMechanism(M.RawMechanism(M.catalytic(am)),
                                              [:EqualAI, :NonequalAI, :OnlyA], 2,
                                              [M.RawSite([:R], 2, [:OnlyI])])
        @test EnzymeRates.AllostericMechanism(raw) == am
        @test M.RawSite(only(am.regulatory_sites)) == M.RawSite([:R], 2, [:OnlyI])
    end
end

@testset "with_flags" begin
    M = EnzymeRates.Mechanisms
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
    # Canonical form stores the steps E + P → E(P), E + S → E(S), E(S) → E(P).
    @test M.with_flags(m, BitVector([1, 0, 0])) ==
          EnzymeRates.Mechanism(@enzyme_mechanism begin
              substrates: S
              products: P
              steps: begin
                  E + S <--> E(S)
                  E(S) <--> E(P)
                  E(P) ⇌ E + P
              end
          end)
    @test M.with_flags(m, BitVector([1, 1, 0])) == m

    am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)      :: NonequalAI
            E(S) <--> E(P)    :: OnlyA
            E(P) ⇌ E + P      :: EqualAI
        end
    end)
    @test M.with_flags(am, BitVector([1, 0, 0])) ==
          EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
              substrates: S
              products: P
              steps: begin
                  E + S <--> E(S)   :: NonequalAI
                  E(S) <--> E(P)    :: OnlyA
                  E(P) ⇌ E + P      :: EqualAI
              end
          end)
end

@testset "regroup" begin
    M = EnzymeRates.Mechanisms
    m = EnzymeRates.Mechanism(@enzyme_mechanism begin
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
    # Canonical form stores the steps E + A → E(A), E + B → E(B), E + Q → E(Q),
    # E(A) + B → E(A, B), E(A, B) → E(P, Q), E(B) + A → E(A, B), E(Q) + P → E(P, Q), in
    # the groups [1, 2, 3, 2, 4, 1, 5]. The groups below put E(B) + A in a group of its own.
    split = [1, 2, 3, 2, 4, 6, 5]
    @test M.regroup(m, split) == EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(B) + A ⇌ E(A, B)
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
            E(A, B) <--> E(P, Q)
            E(P, Q) <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end)
    # New groups and new flags at once: E(B) + A also turns steady state.
    @test M.regroup(m, split, BitVector([1, 1, 0, 1, 0, 0, 0])) ==
          EnzymeRates.Mechanism(@enzyme_mechanism begin
              substrates: A, B
              products: P, Q
              steps: begin
                  E + A ⇌ E(A)
                  E(B) + A <--> E(A, B)
                  (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))
                  E(A, B) <--> E(P, Q)
                  E(P, Q) <--> E(Q) + P
                  E(Q) <--> E + Q
              end
          end)
    @test M.regroup(m, [e.group for e in m.steps]) == m

    # Each new group of an allosteric mechanism takes the tag of its first step's group.
    am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))  :: NonequalAI
            (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))  :: EqualAI
            E(A, B) <--> E(P, Q)                :: OnlyA
            E(P, Q) <--> E(Q) + P               :: EqualAI
            E(Q) <--> E + Q                     :: EqualAI
        end
    end)
    @test M.regroup(am, split) ==
          EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
              substrates: A, B
              products: P, Q
              steps: begin
                  E + A ⇌ E(A)                        :: NonequalAI
                  E(B) + A ⇌ E(A, B)                  :: NonequalAI
                  (E + B ⇌ E(B), E(A) + B ⇌ E(A, B))  :: EqualAI
                  E(A, B) <--> E(P, Q)                :: OnlyA
                  E(P, Q) <--> E(Q) + P               :: EqualAI
                  E(Q) <--> E + Q                     :: EqualAI
              end
          end)
end

@testset "with_tags" begin
    M = EnzymeRates.Mechanisms
    am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
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
    # Canonical form numbers the groups E + P → E(P), E + S → E(S), E(S) → E(P).
    @test M.with_tags(am, [:EqualAI, :NonequalAI, :NonequalAI]) ==
          EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
              substrates: S
              products: P
              allosteric_regulators: R
              catalytic_multiplicity: 2
              steps: begin
                  E + S ⇌ E(S)      :: NonequalAI
                  E(S) <--> E(P)    :: NonequalAI
                  E(P) ⇌ E + P      :: EqualAI
              end
              regulatory_site(2): R::OnlyI
          end)
    @test_throws "must match the number of kinetic groups 3" M.with_tags(am,
                                                                    [:EqualAI, :NonequalAI])
end

@testset "canonical_graph: the canonical graph of raw input, aligned but unchecked" begin
    M = EnzymeRates.Mechanisms
    # Uni-uni written with its release first: the graph is the mechanism the checked
    # constructor builds.
    block = :(begin
        substrates: S
        products: P
        steps: begin
            E(P) ⇌ E + P
            E(S) <--> E(P)
            E + S ⇌ E(S)
        end
    end)
    raw = EnzymeRates.DSL.raw_mechanism(block, false)
    @test M.canonical_graph(raw) == EnzymeRates.Mechanism(raw)
    # A declared modifier that no step binds fails the checked constructor; the graph
    # still builds.
    unused = M.RawMechanism(raw.substrates, raw.products, [:I], raw.steps, raw.group)
    @test_throws("modifier I is declared, but no step binds it",
                 EnzymeRates.Mechanism(unused))
    g = M.canonical_graph(unused)
    @test [M.step_text(g, k) for k in eachindex(g.steps)] ==
          ["E + P ⇌ E(P)", "E + S ⇌ E(S)", "E(S) <--> E(P)"]
    # Tied steps that make different conformational changes have no canonical graph.
    @test_throws(ErrorException(
        "Mechanism: a kinetic group holds E → F and F → G, which make different " *
        "conformational changes; the steps of a kinetic group share their constants, so " *
        "they must make the same conformational change"),
        M.canonical_graph(EnzymeRates.DSL.raw_mechanism(:(begin
            substrates: A
            products: P
            steps: begin
                E + A <--> E(A)
                E(A) <--> E(P)
                E(P) <--> E + P
                (E ⇌ F, F ⇌ G, G ⇌ E)
            end
        end), false)))
end
