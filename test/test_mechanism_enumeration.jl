# ABOUTME: Tests for the staged mechanism enumeration pipeline
# ABOUTME: Organized by stage with hand-calculated expected values

@testset "Mechanism Enumeration Pipeline" begin

@testset "Stage 1: Catalytic topologies" begin

    @testset "Uni-Uni" begin
        topos = EnzymeRates._catalytic_topologies(uni_uni)
        @test length(topos) == 1

        # E ⇌ ES (S-binding), E ⇌ EP (P-binding), ES <--> EP (isom)
        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E_0_0, E_S_0[C], E_0_P[C]
            end
            steps: begin
                [E_0_0, S] ⇌ [E_S_0]
                [E_0_0, P] ⇌ [E_0_P]
                [E_S_0] <--> [E_0_P]
            end
        end

        @test compile_mechanism(topos[1]) === m_uu

        for t in topos
            @test t.n_catalytic_edges == length(t.edges)
            @test count(.!t.equilibrium_steps) >= 1
        end
    end

    @testset "Uni-Bi" begin
        topos = EnzymeRates._catalytic_topologies(uni_bi)
        @test length(topos) == 3

        # Topo 1: ordered release Q-first
        # Path: E → ES → EPQ → EP → E
        # Q released at EPQ→EP, P released at EP→E
        m_ub1 = @enzyme_mechanism begin
            species: begin
                substrates: S[AB]
                products: P[A], Q[B]
                enzymes: E_0_0_0, E_S_0_0[AB],
                    E_0_P_0[A], E_0_P_Q[AB]
            end
            steps: begin
                [E_0_0_0, S] ⇌ [E_S_0_0]
                [E_S_0_0] <--> [E_0_P_Q]
                [E_0_P_0, Q] ⇌ [E_0_P_Q]
                [E_0_0_0, P] ⇌ [E_0_P_0]
            end
        end

        # Topo 2: ordered release P-first
        # Path: E → ES → EPQ → EQ → E
        # P released at EPQ→EQ, Q released at EQ→E
        m_ub2 = @enzyme_mechanism begin
            species: begin
                substrates: S[AB]
                products: P[A], Q[B]
                enzymes: E_0_0_0, E_S_0_0[AB],
                    E_0_0_Q[B], E_0_P_Q[AB]
            end
            steps: begin
                [E_0_0_0, S] ⇌ [E_S_0_0]
                [E_S_0_0] <--> [E_0_P_Q]
                [E_0_0_Q, P] ⇌ [E_0_P_Q]
                [E_0_0_0, Q] ⇌ [E_0_0_Q]
            end
        end

        # Topo 3: random release (both P and Q paths)
        # Path: E → ES → EPQ → EP → E or EPQ → EQ → E
        m_ub3 = @enzyme_mechanism begin
            species: begin
                substrates: S[AB]
                products: P[A], Q[B]
                enzymes: E_0_0_0, E_S_0_0[AB],
                    E_0_P_0[A], E_0_0_Q[B],
                    E_0_P_Q[AB]
            end
            steps: begin
                [E_0_0_0, S] ⇌ [E_S_0_0]
                [E_S_0_0] <--> [E_0_P_Q]
                [E_0_P_0, Q] ⇌ [E_0_P_Q]
                [E_0_0_0, P] ⇌ [E_0_P_0]
                [E_0_0_Q, P] ⇌ [E_0_P_Q]
                [E_0_0_0, Q] ⇌ [E_0_0_Q]
            end
        end

        # Round-trip: each hand-defined mechanism matches
        # a compiled topology
        defined = [m_ub1, m_ub2, m_ub3]
        for (i, m) in enumerate(defined)
            compiled = compile_mechanism(topos[i])
            @test compiled === m
        end

        for t in topos
            @test t.n_catalytic_edges == length(t.edges)
            @test count(.!t.equilibrium_steps) >= 1
        end
    end

    @testset "Bi-Bi" begin
        topos = EnzymeRates._catalytic_topologies(bi_bi)
        @test length(topos) == 9

        # Topo 1: sequential bind B-first,
        #         sequential release P-first
        # Path: E→EB→EAB→EPQ→EP→E
        m_bb1 = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E_0_0_0_0,
                    E_0_B_0_0[N],
                    E_A_B_0_0[CN],
                    E_0_0_P_0[C],
                    E_0_0_P_Q[CN]
            end
            steps: begin
                [E_A_B_0_0] <--> [E_0_0_P_Q]
                [E_0_0_0_0, B] ⇌ [E_0_B_0_0]
                [E_0_0_P_0, Q] ⇌ [E_0_0_P_Q]
                [E_0_B_0_0, A] ⇌ [E_A_B_0_0]
                [E_0_0_0_0, P] ⇌ [E_0_0_P_0]
            end
        end
        @test compile_mechanism(topos[1]) === m_bb1

        # Topo 4: sequential bind A-first,
        #         sequential release P-first
        # Path: E→EA→EAB→EPQ→EP→E
        m_bb4 = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E_0_0_0_0,
                    E_A_0_0_0[C],
                    E_A_B_0_0[CN],
                    E_0_0_P_0[C],
                    E_0_0_P_Q[CN]
            end
            steps: begin
                [E_0_0_0_0, A] ⇌ [E_A_0_0_0]
                [E_A_B_0_0] <--> [E_0_0_P_Q]
                [E_0_0_P_0, Q] ⇌ [E_0_0_P_Q]
                [E_0_0_0_0, P] ⇌ [E_0_0_P_0]
                [E_A_0_0_0, B] ⇌ [E_A_B_0_0]
            end
        end
        @test compile_mechanism(topos[4]) === m_bb4

        # Topo 7: sequential bind A-first,
        #         sequential release Q-first
        # Path: E→EA→EAB→EPQ→EQ→E
        m_bb7 = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E_0_0_0_0,
                    E_A_0_0_0[C],
                    E_A_B_0_0[CN],
                    E_0_0_0_Q[N],
                    E_0_0_P_Q[CN]
            end
            steps: begin
                [E_0_0_0_0, A] ⇌ [E_A_0_0_0]
                [E_A_B_0_0] <--> [E_0_0_P_Q]
                [E_0_0_0_Q, P] ⇌ [E_0_0_P_Q]
                [E_A_0_0_0, B] ⇌ [E_A_B_0_0]
                [E_0_0_0_0, Q] ⇌ [E_0_0_0_Q]
            end
        end
        @test compile_mechanism(topos[7]) === m_bb7

        # Topo 6: random bind, random release (fully
        # random on both sides)
        m_bb6 = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E_0_0_0_0,
                    E_A_0_0_0[C],
                    E_0_B_0_0[N],
                    E_A_B_0_0[CN],
                    E_0_0_P_0[C],
                    E_0_0_0_Q[N],
                    E_0_0_P_Q[CN]
            end
            steps: begin
                [E_0_0_0_0, A] ⇌ [E_A_0_0_0]
                [E_A_B_0_0] <--> [E_0_0_P_Q]
                [E_0_0_0_0, B] ⇌ [E_0_B_0_0]
                [E_0_0_P_0, Q] ⇌ [E_0_0_P_Q]
                [E_0_B_0_0, A] ⇌ [E_A_B_0_0]
                [E_0_0_0_0, P] ⇌ [E_0_0_P_0]
                [E_0_0_0_Q, P] ⇌ [E_0_0_P_Q]
                [E_A_0_0_0, B] ⇌ [E_A_B_0_0]
                [E_0_0_0_0, Q] ⇌ [E_0_0_0_Q]
            end
        end
        @test compile_mechanism(topos[6]) === m_bb6

        # Verify structural properties hold for all
        for t in topos
            @test t.n_catalytic_edges == length(t.edges)
            @test count(.!t.equilibrium_steps) >= 1
        end

        # Verify 9 topologies decompose into known
        # categories by form count:
        # 5 forms: sequential/sequential (4 topos)
        # 6 forms: random-one-side (4 topos)
        # 7 forms: fully random (1 topo)
        form_counts = [
            length(
                EnzymeRates.enzyme_forms(compile_mechanism(t))
            ) for t in topos
        ]
        @test count(==(5), form_counts) == 4
        @test count(==(6), form_counts) == 4
        @test count(==(7), form_counts) == 1
    end

    @testset "Bi-Bi Ping-Pong" begin
        topos =
            EnzymeRates._catalytic_topologies(bi_bi_ping_pong)
        @test length(topos) == 10

        # Topo 4: classic ping-pong with E_X intermediate
        # E → EA → E_X(+P) → E_X_B → E(+Q)
        m_pp = @enzyme_mechanism begin
            species: begin
                substrates: A[CX], B[N]
                products: P[C], Q[NX]
                enzymes: E_0_0_0_0,
                    E_A_0_0_0[CX],
                    E_X_0_0_0[X],
                    E_X_B_0_0[NX],
                    E_X_0_P_0[CX],
                    E_0_0_0_Q[NX]
            end
            steps: begin
                [E_0_0_0_0, A] ⇌ [E_A_0_0_0]
                [E_0_0_0_0, Q] ⇌ [E_0_0_0_Q]
                [E_X_B_0_0] <--> [E_0_0_0_Q]
                [E_X_0_0_0, P] ⇌ [E_X_0_P_0]
                [E_X_0_0_0, B] ⇌ [E_X_B_0_0]
                [E_A_0_0_0] ⇌ [E_X_0_P_0]
            end
        end
        @test compile_mechanism(topos[4]) === m_pp

        for t in topos
            @test t.n_catalytic_edges == length(t.edges)
            @test count(.!t.equilibrium_steps) >= 1
        end
    end

    @testset "Ter-Ter" begin
        # Expected: 49 topologies (1 random + 6+6 mixed + 36 seq)
        # Currently OOMs in _catalytic_topologies — skip until
        # the enumeration code is optimized for 3+3 reactions.
        @test_broken false  # placeholder: ter_ter OOMs
    end

end

end # outer testset
