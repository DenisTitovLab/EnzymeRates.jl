@testset "Structural identifiability detection" begin
    # Simple Uni-Uni (E + S ⇌ ES ⇌ E + P)
    # 2 steps, 4 k's, 1 Haldane → 3 independent k's
    # Rate: (a*S - b*P)/(1 + c*S + d*P)
    #   Numerator: 2 monomials (S, P)
    #   Denominator: 3 monomials (1, S, P)
    # Identifiable coefficients: (2-1) + (3-1) = 3
    # Deficit: 3 - 3 = 0
    @testset "Uni-Uni is identifiable" begin
        m_uu, _ = make_uni_uni()
        @test structural_identifiability_deficit(m_uu) == 0
        @test is_identifiable(m_uu)
    end

    # 3-step isomerization (E + S ⇌ ES ⇌ ES' ⇌ E + P)
    # 3 steps, 6 k's, 1 Haldane → 5 independent k's
    # Rate: same form as uni-uni → 3 identifiable coefficients
    # Deficit: 5 - 3 = 2
    @testset "Three-step isomerization is non-identifiable" begin
        m_3step, _ = make_three_step_isomerization()
        @test structural_identifiability_deficit(m_3step) == 2
        @test !is_identifiable(m_3step)
    end

    # Sequential Bi-Bi (E + S1 ⇌ ES1 ⇌ ES1S2 ⇌ EP1P2 ⇌ EP2 ⇌ E + P2)
    # 5 steps, 10 k's, 1 Haldane → 9 independent k's
    # Overdetermined: more identifiable coefficients than parameters
    @testset "Sequential Bi-Bi (overdetermined)" begin
        m_bb, _ = make_seq_bibi()
        deficit = structural_identifiability_deficit(m_bb)
        @test deficit < 0  # Overdetermined
        @test is_identifiable(m_bb)
        # The mechanism has 9 independent k's after Haldane
        # Compute number of independent params from parameters(m) minus Keq and E_total
        all_params = parameters(m_bb)
        n_indep = length(all_params) - 2  # minus Keq and E_total
        @test n_indep == 9
    end

    # Random-order Bi-Bi (has Wegscheider condition)
    # 7 steps, 14 k's, 1 Haldane + 1 Wegscheider → 12 independent k's
    # Wegscheider reduces k's but NOT identifiable coefficients
    # Overdetermined: more identifiable coefficients than parameters
    @testset "Random-order Bi-Bi (Wegscheider, overdetermined)" begin
        m_ro, _ = make_random_bibi()
        deficit = structural_identifiability_deficit(m_ro)
        @test deficit < 0  # Overdetermined
        @test is_identifiable(m_ro)
        # Verify the mechanism has Wegscheider (2 constraints total)
        dep_exprs, indep = EnzymeRates._dependent_param_exprs(typeof(m_ro))
        @test length(dep_exprs) == 2
        @test length(indep) == 12
    end

    # Doubly branched mechanism (has Wegscheider cycle)
    # 5 steps, 10 k's, 1 Haldane + 1 Wegscheider → 8 independent k's
    @testset "Doubly branched (Wegscheider cycle)" begin
        m_db, _ = make_doubly_branched()
        deficit = structural_identifiability_deficit(m_db)
        @test deficit > 0  # Underdetermined
        @test !is_identifiable(m_db)
        # Verify the mechanism has 2 constraints (Haldane + Wegscheider)
        dep_exprs, indep = EnzymeRates._dependent_param_exprs(typeof(m_db))
        @test length(dep_exprs) == 2
        @test length(indep) == 8
    end

    # Ping-pong Bi-Bi
    @testset "Ping-pong Bi-Bi" begin
        m_pp, _ = make_pingpong_bibi()
        deficit = structural_identifiability_deficit(m_pp)
        @test deficit > 0  # Underdetermined
        @test !is_identifiable(m_pp)
        @test deficit isa Int
    end

    @testset "Monomial counting" begin
        # Test the internal monomial counting function
        m_uu, _ = make_uni_uni()
        n_num, n_denom = EnzymeRates._count_rate_monomials(typeof(m_uu))
        # Uni-Uni rate: (a*S - b*P)/(1 + c*S + d*P)
        @test n_num == 2  # S and P terms
        @test n_denom == 3  # 1, S, and P terms

        m_3step, _ = make_three_step_isomerization()
        n_num_3, n_denom_3 = EnzymeRates._count_rate_monomials(typeof(m_3step))
        # Same rate equation form
        @test n_num_3 == 2
        @test n_denom_3 == 3
    end

    @testset "Consistency checks" begin
        mechanisms = [
            (make_uni_uni, "Uni-Uni"),
            (make_three_step_isomerization, "Three-step iso"),
            (make_seq_bibi, "Seq Bi-Bi"),
            (make_random_bibi, "Random Bi-Bi"),
            (make_doubly_branched, "Doubly branched"),
            (make_pingpong_bibi, "Ping-pong Bi-Bi"),
        ]

        for (make_fn, label) in mechanisms
            m, _ = make_fn()
            _, indep = EnzymeRates._dependent_param_exprs(typeof(m))
            n_indep = length(indep)
            deficit = structural_identifiability_deficit(m)
            n_num, n_denom = EnzymeRates._count_rate_monomials(typeof(m))
            n_identifiable = (n_num - 1) + (n_denom - 1)

            # Verify the formula: deficit = n_k - n_identifiable (can be negative for overdetermined)
            @test deficit == n_indep - n_identifiable

            # Verify is_identifiable consistency (identifiable when deficit <= 0)
            @test is_identifiable(m) == (deficit <= 0)
        end
    end
end
