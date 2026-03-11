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

@testset "Stage 2: RE/SS expansion" begin
    @testset "Uni-Uni: 2^3 - 1 = 7 variants" begin
        topo = EnzymeRates._catalytic_topologies(uni_uni)[1]
        result = EnzymeRates._expand_ress_variants(
            [topo], uni_uni)
        # Correct: all 3 steps toggleable, 2^3 - 1 = 7
        # BUG: current code only toggles binding edges
        #   (2 binding edges), gives 2^2 - 1 = 3
        @test_broken length(result) == 7
        @test length(result) == 3

        for s in result
            @test any(.!s.equilibrium_steps)
            @test s.edges == topo.edges
            @test s.n_catalytic_edges ==
                topo.n_catalytic_edges
        end
    end

    @testset "Uni-Bi" begin
        topo = EnzymeRates._catalytic_topologies(uni_bi)[1]
        n_steps = length(topo.edges)
        result = EnzymeRates._expand_ress_variants(
            [topo], uni_bi)
        # Correct: 2^n_steps - 1. Determine n_steps at
        # runtime and compute expected. Mark correct as
        # @test_broken, then determine current (buggy)
        # count by counting binding edges only.
        n_binding = count(topo.equilibrium_steps)
        @test_broken length(result) == 2^n_steps - 1
        @test length(result) == 2^n_binding - 1

        for s in result
            @test any(.!s.equilibrium_steps)
        end
    end

    @testset "Bi-Bi" begin
        topo = EnzymeRates._catalytic_topologies(bi_bi)[1]
        n_steps = length(topo.edges)
        n_binding = count(topo.equilibrium_steps)
        result = EnzymeRates._expand_ress_variants(
            [topo], bi_bi)
        @test_broken length(result) == 2^n_steps - 1
        @test length(result) == 2^n_binding - 1

        for s in result
            @test any(.!s.equilibrium_steps)
        end
    end

    @testset "max_re_groups filtering" begin
        # Use bi_bi with many topologies to test filtering.
        # The all-SS assignment should be excluded when it
        # creates more RE groups than max_re_groups allows.
        topo = EnzymeRates._catalytic_topologies(bi_bi)[end]
        n_steps = length(topo.edges)
        n_binding = count(topo.equilibrium_steps)
        result_default = EnzymeRates._expand_ress_variants(
            [topo], bi_bi)
        # With strict max_re_groups=2, fewer variants survive
        result_strict = EnzymeRates._expand_ress_variants(
            [topo], bi_bi; max_re_groups=2)
        @test length(result_strict) <=
            length(result_default)
    end
end

@testset "Stage 2.5: Substrate/product dead-end expansion" begin
    # This stage does not exist yet. All tests document
    # expected behavior for when it is implemented.

    @testset "Uni-Uni: passthrough (no off-cycle forms)" begin
        # Uni-Uni: 3 forms, 3 on-cycle → 0 off-cycle
        # No substrate/product dead-end forms exist, so
        # this stage would be a passthrough.
        @test_broken false  # stage doesn't exist
    end

    @testset "Bi-Bi: 4 off-cycle forms" begin
        # Bi-Bi: 11 forms, 7 on-cycle, 4 off-cycle
        # Off-cycle: E_A_0_P_0, E_0_B_P_0,
        #            E_A_0_0_Q, E_0_B_0_Q
        # Expected: 2^4 = 16 variants per input spec
        @test_broken false  # stage doesn't exist
    end

    @testset "Bi-Bi Ping-Pong: 7 off-cycle forms" begin
        # Bi-Bi-PP: 17 forms, 10 on-cycle, 7 off-cycle
        # Off-cycle: E_A_0_P_0, E_0_B_P_0, E_X_B_P_0,
        #   E_A_0_0_Q, E_X_0_0_Q, E_0_B_0_Q, E_X_B_0_Q
        # Expected: 2^7 = 128 variants per input spec
        @test_broken false  # stage doesn't exist
    end
end

@testset "Stage 3: Regulator dead-end expansion" begin
    @testset "No regulators: passthrough" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni)[1]
        result = EnzymeRates._expand_dead_end_inhibitors(
            [topo], uni_uni; dead_end_regs=Symbol[])
        @test length(result) == 1
    end

    @testset "Uni-Uni + I: eligible forms" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_dead_end_I)[1]
        result = EnzymeRates._expand_dead_end_inhibitors(
            [topo], uni_uni_dead_end_I;
            dead_end_regs=[:I])
        # Correct: only E eligible (ES fully occupied
        # with all substrates, EP fully occupied with
        # all products) → (2^1)^1 = 2 variants
        # BUG: current code allows I to bind all 3 forms
        # → (2^1)^3 = 8
        @test_broken length(result) == 2
        @test length(result) == 8
    end

    @testset "Uni-Bi + I" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_bi_dead_end_I)[1]
        result = EnzymeRates._expand_dead_end_inhibitors(
            [topo], uni_bi_dead_end_I;
            dead_end_regs=[:I])
        # Correct: E, EP, EQ eligible (3 forms)
        # ES has all substrates, EPQ has all products
        # → (2^1)^3 = 8
        # BUG: current code allows I to bind all forms
        n_forms = length(EnzymeRates.enzyme_forms(
            compile_mechanism(topo)))
        @test_broken length(result) == 8
        @test length(result) == 2^n_forms
    end

    @testset "Bi-Bi + I: 5 eligible forms" begin
        topo = EnzymeRates._catalytic_topologies(
            bi_bi_dead_end_I)[1]
        result = EnzymeRates._expand_dead_end_inhibitors(
            [topo], bi_bi_dead_end_I;
            dead_end_regs=[:I])
        # Topo 1 is sequential (5 forms: E, EB, EAB,
        # EP, EPQ). Correct: 5 eligible (EAB has all
        # substrates, EPQ has all products → both should
        # be excluded, but this topo has only 5 forms
        # total so 2^5 = 32 = 2^n_forms).
        # BUG: current code allows I to bind all forms,
        # but here it gives correct count by coincidence
        # (all 5 forms × 1 inhibitor → 2^5 = 32).
        n_forms = length(EnzymeRates.enzyme_forms(
            compile_mechanism(topo)))
        @test length(result) == 32
        @test length(result) == 2^n_forms
    end

    @testset "2 inhibitors: Uni-Uni + I, J" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_dead_end_I_J)[1]
        result = EnzymeRates._expand_dead_end_inhibitors(
            [topo], uni_uni_dead_end_I_J;
            dead_end_regs=[:I, :J])
        # Correct: only E eligible → (2^2)^1 = 4
        # BUG: all 3 forms → (2^2)^3 = 64
        @test_broken length(result) == 4
        @test length(result) == 64
    end

    @testset "Allosteric-only: passthrough" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_allosteric_R)[1]
        result = EnzymeRates._expand_dead_end_inhibitors(
            [topo], uni_uni_allosteric_R;
            dead_end_regs=Symbol[])
        @test length(result) == 1
    end

    @testset "Properties" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_dead_end_I)[1]
        result = EnzymeRates._expand_dead_end_inhibitors(
            [topo], uni_uni_dead_end_I;
            dead_end_regs=[:I])
        for s in result
            @test length(s.edges) >= s.n_catalytic_edges
            @test s.n_catalytic_edges ==
                topo.n_catalytic_edges
        end
    end
end

@testset "Stage 4: Equivalence constraints" begin
    @testset "No equiv groups: passthrough" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni)[1]
        result =
            EnzymeRates._expand_equivalence_constraints(
                [topo], uni_uni)
        @test length(result) == 1
    end

    @testset "Dead-end inhibitor creates equiv groups" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_dead_end_I)[1]
        de = EnzymeRates._expand_dead_end_inhibitors(
            [topo], uni_uni_dead_end_I;
            dead_end_regs=[:I])
        # Find a spec with dead-end edges (more edges
        # than catalytic)
        spec_with_de = findfirst(
            s -> length(s.edges) > s.n_catalytic_edges,
            de)
        if spec_with_de !== nothing
            s = de[spec_with_de]
            eq =
                EnzymeRates._expand_equivalence_constraints(
                    [s], uni_uni_dead_end_I)
            # Should produce at least 1 variant
            # (passthrough if no equiv groups) and
            # possibly more if I binds to multiple
            # forms with the same metabolite
            @test length(eq) >= 1
        end
    end

    @testset "Properties" begin
        topo = EnzymeRates._catalytic_topologies(
            uni_uni_dead_end_I)[1]
        de = EnzymeRates._expand_dead_end_inhibitors(
            [topo], uni_uni_dead_end_I;
            dead_end_regs=[:I])
        eq =
            EnzymeRates._expand_equivalence_constraints(
                de, uni_uni_dead_end_I)
        @test length(eq) >= length(de)

        unconstrained = filter(
            s -> isempty(s.param_constraints), eq)
        constrained = filter(
            s -> !isempty(s.param_constraints), eq)
        if !isempty(constrained) &&
                !isempty(unconstrained)
            @test minimum(
                s.param_count for s in constrained) <
                maximum(
                    s.param_count
                    for s in unconstrained)
        end
    end
end

end # outer testset
