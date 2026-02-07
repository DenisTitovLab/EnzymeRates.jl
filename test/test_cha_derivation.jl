# Tests for Cha method (mixed rapid-equilibrium + steady-state steps)

using OrdinaryDiffEqFIRK

# ── RE group computation tests ────────────────────────────────────────────────

@testset "RE Groups" begin
    @testset "All SS = singleton groups" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products:   P[C]
                enzymes:    E, ES[C]
            end
            steps: begin
                [E, S] <--> [ES]
                [ES] <--> [E, P]
            end
        end
        enzs = EnzymeRates.enzyme_forms(m)
        rxns = EnzymeRates.reactions(m)
        eq = EnzymeRates.equilibrium_steps(m)
        enz_names = Tuple(e[1] for e in enzs)
        enz_set = Set(enz_names)
        groups, ftg = EnzymeRates._compute_re_groups(enz_names, enz_set, rxns, eq)
        @test length(groups) == 2
        @test all(length(g) == 1 for g in groups)
    end

    @testset "One RE step = two forms merged" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C]
                products:   P[C]
                enzymes:    E, EA[C], EP[C]
            end
            steps: begin
                [E, A] ⇌ [EA]
                [EA] <--> [EP]
                [EP] <--> [E, P]
            end
        end
        enzs = EnzymeRates.enzyme_forms(m)
        rxns = EnzymeRates.reactions(m)
        eq = EnzymeRates.equilibrium_steps(m)
        enz_names = Tuple(e[1] for e in enzs)
        enz_set = Set(enz_names)
        groups, ftg = EnzymeRates._compute_re_groups(enz_names, enz_set, rxns, eq)
        # E and EA should be in one group (connected by RE step)
        @test length(groups) == 2  # {E, EA} and {EP}
        # Find group containing E (index 1)
        e_idx = findfirst(==(enz_names[1]), enz_names)
        ea_idx = findfirst(==(:EA), enz_names)
        @test ftg[e_idx] == ftg[ea_idx]
    end
end

# ── Simple RE Uni-Uni: E + A ⇌_RE EA <-->_SS E + P ──────────────────────────

@testset "RE Uni-Uni (rapid equilibrium binding)" begin
    m = @enzyme_mechanism begin
        species: begin
            substrates: A[C]
            products:   P[C]
            enzymes:    E, EA[C]
        end
        steps: begin
            [E, A] ⇌ [EA]
            [EA] <--> [E, P]
        end
    end

    @testset "Structure" begin
        @test EnzymeRates.n_states(m) == 2
        @test EnzymeRates.n_steps(m) == 2
        @test EnzymeRates.equilibrium_steps(m) == (true, false)
    end

    @testset "Parameters Raw" begin
        p = parameters(m, EnzymeRates.Raw)
        @test :K1 in p
        @test :k2f in p
        @test :k2r in p
        @test :E_total in p
        @test length(p) == 4  # K1, k2f, k2r, E_total
    end

    @testset "Rapid-equilibrium Michaelis-Menten" begin
        # For E + A ⇌_RE EA <-->_SS E + P:
        # At rapid equilibrium: K1 = [EA]/([E][A])
        # Rate = E_total * (k2f * K1 * A - k2r * P) / (1 + K1 * A)
        # This is the classic rapid-equilibrium Michaelis-Menten
        for _ in 1:20
            K1 = 0.1 + 9.9 * rand()
            k2f = 0.1 + 9.9 * rand()
            k2r = 0.1 + 9.9 * rand()
            Et = 0.1 + 9.9 * rand()
            A = 0.1 + 9.9 * rand()
            P = 0.1 + 9.9 * rand()

            params = (K1=K1, k2f=k2f, k2r=k2r, E_total=Et)
            concs = (A=A, P=P)

            v_cha = rate_equation(m, params, concs, EnzymeRates.Raw)

            # Analytical rapid-equilibrium formula
            num = k2f * K1 * A - k2r * P
            denom = 1.0 + K1 * A
            v_analytical = Et * num / denom

            @test v_cha ≈ v_analytical rtol=1e-10
        end
    end

    @testset "Show method" begin
        s = sprint(show, m)
        @test contains(s, "⇌")
        @test contains(s, "<-->")
    end

    @testset "Zero allocation" begin
        params = (K1=1.0, k2f=2.0, k2r=0.5, E_total=1.0)
        concs = (A=1.0, P=0.5)
        rate_equation(m, params, concs, EnzymeRates.Raw)  # warmup
        allocs = @allocated rate_equation(m, params, concs, EnzymeRates.Raw)
        @test allocs == 0
    end
end

# ── 3-step with RE binding: E + A ⇌_RE EA <-->_SS EP <-->_SS E + P ──────────

@testset "3-step RE binding" begin
    m = @enzyme_mechanism begin
        species: begin
            substrates: A[C]
            products:   P[C]
            enzymes:    E, EA[C], EP[C]
        end
        steps: begin
            [E, A] ⇌ [EA]
            [EA] <--> [EP]
            [EP] <--> [E, P]
        end
    end

    @testset "Structure" begin
        @test EnzymeRates.equilibrium_steps(m) == (true, false, false)
    end

    @testset "Parameters" begin
        p = parameters(m, EnzymeRates.Raw)
        @test :K1 in p
        @test :k2f in p
        @test :k2r in p
        @test :k3f in p
        @test :k3r in p
        @test length(p) == 6  # K1, k2f, k2r, k3f, k3r, E_total
    end

    @testset "Analytical verification" begin
        # E + A ⇌_RE EA <-->_SS EP <-->_SS E + P
        # Groups: {E, EA} (with alpha_EA = K1 * A) and {EP} (singleton)
        # Rate matrix (2x2 over groups):
        #   R[{E,EA}, {EP}] = k2f * K1 * A  (forward through step 2, alpha = K1*A)
        #   R[{EP}, {E,EA}] = k2r             (reverse through step 2, alpha = 1)
        #   R[{EP}, {E,EA}] += k3f            (forward through step 3, alpha = 1)
        #   R[{E,EA}, {EP}] += k3r * P        (reverse through step 3, alpha = 1)
        # Wait, step 3: EP <--> E + P, EP is in group {EP}, E is in group {E, EA}
        # So R[{EP}, {E,EA}] += k3f (forward from EP to E)
        #    R[{E,EA}, {EP}] += k3r * P (reverse from E to EP, alpha_E = 1)

        for _ in 1:20
            K1 = 0.1 + 9.9 * rand()
            k2f = 0.1 + 9.9 * rand()
            k2r = 0.1 + 9.9 * rand()
            k3f = 0.1 + 9.9 * rand()
            k3r = 0.1 + 9.9 * rand()
            Et = 0.1 + 9.9 * rand()
            A = 0.1 + 9.9 * rand()
            P = 0.1 + 9.9 * rand()

            params = (K1=K1, k2f=k2f, k2r=k2r, k3f=k3f, k3r=k3r, E_total=Et)
            concs = (A=A, P=P)

            v_cha = rate_equation(m, params, concs, EnzymeRates.Raw)

            # Analytical: 2 groups, sigma1 = 1 + K1*A, sigma2 = 1
            # R12 = k2f * K1 * A + k3r * P (from group 1 to group 2)
            # R21 = k2r + k3f (from group 2 to group 1)
            # For 2x2: D[1] = R21, D[2] = R12
            # denom = sigma1 * D[1] + sigma2 * D[2] = (1 + K1*A) * (k2r + k3f) + (k2f * K1 * A + k3r * P)
            # Numerator (for substrate A consumed in RE step 1, but A only in RE step):
            # We use P which appears in SS step 3:
            # flux_3 = Et * (k3f * D[{EP}] - k3r * P * D[{E,EA}]) / denom
            # v(P produced) = flux_3, v(A consumed) = -v(P produced) for 1:1 stoichiometry
            # Actually let's just compute numerically

            sigma1 = 1 + K1 * A
            R12 = k2f * K1 * A + k3r * P  # group1 -> group2
            R21 = k2r + k3f                 # group2 -> group1
            D1 = R21  # cofactor for group 1 (delete row/col 1 from 2x2 Laplacian)
            D2 = R12  # cofactor for group 2
            denom_val = sigma1 * D1 + 1 * D2

            # Flux of A: consumed in step 1 (RE), but step 1 is RE so doesn't appear
            # in SS flux. Use step 3 flux involving P:
            # step 3: EP <--> E + P, EP in group 2 (g2), E in group 1 (g1)
            # rf3 = k3f * alpha_EP = k3f (EP is singleton, alpha=1)
            # rr3 = k3r * P * alpha_E = k3r * P (alpha_E = 1)
            # flux_3 = Et * (rf3 * D[g2] - rr3 * D[g1]) / denom
            #        = Et * (k3f * D2 - k3r * P * D1) / denom

            # P is produced in step 3 (on RHS), so v_P = flux_3
            # Net: v_A (consumption) = -v_P = -flux_3 ... wait, that's for 1:1 stoich
            # Actually: A consumed = P produced for this mechanism (same C atoms)
            # v = consumption rate of A = production rate of P

            # Step 3 produces P (on rhs). So P_produced = flux_step3.
            # v = P_produced = flux_step3 (since |nu_A| = |nu_P| = 1)
            # But the function computes v as consumption of A (first substrate)

            # Since A only appears in RE step 1, the fallback kicks in:
            # It finds P in SS step 3, nu_A = -1, nu_P = +1
            # alt_flux = sum over SS steps involving P = flux_3 (P on rhs of step 3: -flux_3)
            # Wait, let me trace through the code more carefully.
            # ref_name = A, not in SS steps -> fallback
            # alt_name = P, nu_ref = -1 (net stoich of A), nu_alt = +1 (net stoich of P)
            # flux computation for P in step 3:
            #   met_r = P (rhs of step 3), so push(-flux_3) (consuming P is negative)
            #   Wait no: step 3 is EP <--> E + P
            #   met_f = nothing (lhs has only EP = enzyme), met_r = P
            #   So met_r === alt_name (P), push(-flux_3)
            # alt_flux = -flux_3 = -(Et * (k3f * D2 - k3r * P * D1) / denom)
            # ratio = nu_ref / nu_alt = -1 / 1 = -1
            # net_expr = alt_flux * ratio = -(-flux_3) * 1 = flux_3... wait
            # Actually ratio = nu_ref // nu_alt = -1 // 1 = -1
            # And the code does: return :(-$alt_flux) for ratio == -1
            # So net_expr = -alt_flux = -(- flux_3) = flux_3
            # Then abs_nu = 1, so result = net_expr = flux_3
            # = Et * (k3f * D2 - k3r * P * D1) / denom_val
            v_analytical = Et * (k3f * D2 - k3r * P * D1) / denom_val

            @test v_cha ≈ v_analytical rtol=1e-10
        end
    end
end

# ── Both bindings RE: E + A ⇌_RE EA <-->_SS EP ⇌_RE E + P ──────────────────

@testset "Both bindings RE" begin
    m = @enzyme_mechanism begin
        species: begin
            substrates: A[C]
            products:   P[C]
            enzymes:    E, EA[C], EP[C]
        end
        steps: begin
            [E, A] ⇌ [EA]
            [EA] <--> [EP]
            [EP] ⇌ [E, P]
        end
    end

    @testset "Structure" begin
        @test EnzymeRates.equilibrium_steps(m) == (true, false, true)
    end

    @testset "Analytical verification" begin
        # Groups: {E, EA} and {EP, E_again}?? No -- E is shared!
        # Wait: step 1 connects E and EA (RE), step 3 connects EP and E (RE)
        # So E, EA, EP are ALL in one RE group!
        # That means G = 1 (single group), which means there's no rate matrix...
        # Actually this can't work: if all forms are in one RE group, there's no SS flux
        # between different groups. We need at least 2 groups for the SS step to connect.

        # Let me re-examine: step 2 is SS connecting EA (form) and EP (form)
        # If E, EA, EP are all in one group via RE steps, then g1 == g2 for step 2
        # and the SS step is a self-loop which doesn't contribute to the rate matrix.

        # Actually the RE groups are computed only from RE steps:
        # Step 1 (RE): E ↔ EA → group {E, EA}
        # Step 3 (RE): EP ↔ E → group {EP, E}
        # So by transitivity: {E, EA, EP} all in one group!

        # With G = 1, rate matrix is 1x1 (zero, since diagonal is excluded)
        # Cofactor D[1] = 1 (0x0 determinant)
        # Denominator = sigma * 1 = (1 + K1*A + P/(K3*A))... hmm wait

        # This is a degenerate case. Let me check what the code does.
        # The alpha computation:
        # ref = E (first form in group), alpha_E = 1
        # Step 1 (RE, forward E→EA): alpha_EA = alpha_E * K1 * A = K1 * A
        # Step 3 (RE, connecting EP and E): EP <--> E + P
        #   Step 3 forward is EP → E + P
        #   Traversing from E to EP (reverse): alpha_EP = alpha_E / K3 * ... hmm

        # Let's just test the numerator fallback: all metabolites in RE steps
        # No metabolites in SS step (step 2 is EA <--> EP, no metabolites)
        # So the "ALL metabolites in RE steps" branch is triggered
        # Net flux through step 2 (the only SS step):
        # flux_2 = Et * (k2f * alpha_EA * D[1] - k2r * alpha_EP * D[1]) / (sigma * D[1])
        #        = Et * (k2f * K1*A - k2r * alpha_EP) / sigma

        # For step 3 (RE): EP ⇌ E + P means K3 = [E][P] / [EP]
        # So [EP] = [E][P]/K3, alpha_EP = alpha_E * P / K3 = P/K3

        # Wait, let me trace the BFS more carefully:
        # Starting from ref=E (alpha=1), looking at step 1 (RE):
        #   Forward (E→EA): alpha_EA = 1 * K1 * A = K1*A  ✓
        # Now from EA, looking at step 3 (RE) which connects EP and E:
        #   EA is visited, so this step doesn't fire from EA.
        # From E (already visited), looking at step 3 (RE):
        #   Step 3: (EP,) -> (E, P)
        #   e_lhs = EP, e_rhs = E
        #   j_form = e_rhs position. i_form = e_lhs position.
        #   current = E index, so j_form (E index) == current → need i_form (EP index) ∉ visited
        #   j_form == current means we go to "reverse traversal" branch
        #   alpha_EP = alpha[current] / K3 * P... wait let me trace the code exactly

        # In the BFS code:
        #   For step 3: lhs = (:EP,), rhs = (:E, :P)
        #   e_lhs = EP, e_rhs = E, m_lhs = [], m_rhs = [P]
        #   i_form = index of EP, j_form = index of E
        #   Current = index of E
        #
        #   Check: i_form == current? No (EP != E index)
        #   Check: j_form == current? Yes (E == E index)
        #   So reverse traversal: i_form (EP) not in visited? Yes!
        #   alpha[EP] = alpha[current] * [Met_rhs] / (K3 * [Met_lhs])
        #   num_factors = [alpha[current]] + m_rhs concs = [1, P]
        #   denom_factors = [K3] + m_lhs concs = [K3] (m_lhs is empty)
        #   alpha_EP = P / K3

        # sigma = 1 + K1*A + P/K3
        # The SS step 2 (EA <--> EP, no metabolites) is a self-loop in the single group
        # g1 == g2, so it doesn't contribute to the rate matrix
        #
        # With G=1, there's only one group, D[1] = 1
        # denom = sigma * 1 = 1 + K1*A + P/K3
        # The fallback for all-RE-metabolites isomerization:
        # step 2 is the only SS step (EA <--> EP, isomerization)
        # g1 = g2 = 1 (same group!), so D[g1] == D[g2]
        # flux = Et * (k2f * alpha_EA - k2r * alpha_EP) * D[1] / denom... wait
        # flux = Et * (k2f*alpha_EA*D[g1] - k2r*alpha_EP*D[g2]) / denom
        # But g1 == g2 so D[g1] == D[g2] = 1
        # flux = Et * (k2f * K1*A - k2r * P/K3) / (1 + K1*A + P/K3)
        # And since nu_ref (A) = -1, sign_factor = sign(-1) = -1
        # result = -flux (before dividing by abs_nu=1)
        # Hmm, this gives a negative sign. Let me think...
        # Actually nu_ref = -1 (consumed), sign_factor = -1
        # return :(-$flux)
        # So the rate = -flux = -Et * (k2f*K1*A - k2r*P/K3) / (sigma)
        # = Et * (k2r*P/K3 - k2f*K1*A) / sigma
        # That's wrong — the rate should be positive when A > P*Keq relation

        # Wait, I think the issue is the sign convention. Let me reconsider.
        # The rate is defined as net consumption of first substrate (A).
        # In step 2 (EA → EP), going forward consumes EA and produces EP.
        # This doesn't directly consume/produce A. But the RE equilibrium means:
        # consuming EA shifts the equilibrium E + A ⇌ EA, consuming A.
        # The Cha method handles this automatically through the group structure.

        # For the unicyclic case with G=1:
        # The net flux through step 2 = k2f * [EA] - k2r * [EP]
        # = k2f * alpha_EA * V - k2r * alpha_EP * V
        # where V = Et / sigma
        # = V * (k2f * K1*A - k2r * P/K3)
        # = Et * (k2f * K1*A - k2r * P/K3) / (1 + K1*A + P/K3)

        # This flux is the rate of conversion EA → EP, which equals the overall reaction rate.
        # d[A]/dt = -v, so v = positive when A is consumed (forward reaction dominant)
        # When k2f*K1*A > k2r*P/K3, forward dominates, A is consumed, v > 0
        # The flux through step 2 IS the overall rate (positive = forward = A consumed)

        # So v = flux_2 = Et * (k2f * K1*A - k2r * P/K3) / sigma
        # But the code returns -flux with sign_factor = -1...
        # Let me re-examine. sign(nu_ref) = sign(-1) = -1
        # The code says: sign_factor == -1 ? :(-$flux) : flux
        # flux = Et * (rf * D1 - rr * D2) / denom
        # For step 2: rf = k2f * alpha_EA = k2f * K1*A, rr = k2r * alpha_EP = k2r * P/K3
        # D1 = D2 = 1 (same group)
        # flux = Et * (k2f*K1*A - k2r*P/K3) / sigma
        # return -flux... that's wrong!

        # Actually, I think the issue is that for isomerization steps in a single group,
        # the sign depends on which direction represents consumption of substrate.
        # For A: the forward direction of step 2 (EA→EP) means A is consumed.
        # But the fallback code just uses sign(nu_ref) = sign(-1) = -1 and returns -flux.
        # This gives -(positive) = negative, which means the rate is negative when forward
        # dominates. That's incorrect.

        # I think for G=1 all-isomerization case, we should return flux directly
        # because the flux already has the correct sign (positive = forward = A consumed).
        # Let me fix this in the fallback code.

        # Actually, let me reconsider. In the standard QSSA case where G=N (all SS),
        # the numerator works by computing consumption/production of the reference substrate.
        # Step 2 is an isomerization, so ref substrate A doesn't appear directly.
        # In the standard case, the rate would be computed from step 1 which involves A.
        # But here step 1 is RE, so we need the fallback.

        # The correct answer for this mechanism is:
        # v = Et * (k2f * K1 * A - k2r * P / K3) / (1 + K1*A + P/K3)

        for _ in 1:20
            K1 = 0.1 + 9.9 * rand()
            k2f = 0.1 + 9.9 * rand()
            k2r = 0.1 + 9.9 * rand()
            K3 = 0.1 + 9.9 * rand()
            Et = 0.1 + 9.9 * rand()
            A = 0.1 + 9.9 * rand()
            P = 0.1 + 9.9 * rand()

            params = (K1=K1, k2f=k2f, k2r=k2r, K3=K3, E_total=Et)
            concs = (A=A, P=P)

            v_cha = rate_equation(m, params, concs, EnzymeRates.Raw)

            sigma = 1.0 + K1 * A + P / K3
            v_analytical = Et * (k2f * K1 * A - k2r * P / K3) / sigma

            @test v_cha ≈ v_analytical rtol=1e-10
        end
    end
end

# ── All SS = full QSSA equivalence ───────────────────────────────────────────

@testset "All-SS equals QSSA" begin
    # Verify that when all steps are SS, Cha gives identical results to standard QSSA
    # Use the Uni-Uni mechanism from the main test suite
    m = @enzyme_mechanism begin
        species: begin
            substrates: S[C]
            products:   P[C]
            enzymes:    E, ES[C]
        end
        steps: begin
            [E, S] <--> [ES]
            [ES] <--> [E, P]
        end
    end

    @test EnzymeRates.equilibrium_steps(m) == (false, false)

    rng = Random.MersenneTwister(42)
    for _ in 1:20
        k1f = 0.1 + 9.9 * rand(rng)
        k1r = 0.1 + 9.9 * rand(rng)
        k2f = 0.1 + 9.9 * rand(rng)
        k2r = 0.1 + 9.9 * rand(rng)
        Et = 0.1 + 9.9 * rand(rng)
        S = 0.1 + 9.9 * rand(rng)
        P = 0.1 + 9.9 * rand(rng)

        params = (k1f=k1f, k1r=k1r, k2f=k2f, k2r=k2r, E_total=Et)
        concs = (S=S, P=P)

        v = rate_equation(m, params, concs, EnzymeRates.Raw)

        # Standard QSSA Uni-Uni formula
        num = k1f * k2f * S - k1r * k2r * P
        denom = k1r + k2f + k1f * S + k2r * P
        v_qssa = Et * num / denom

        @test v ≈ v_qssa rtol=1e-10
    end
end

# ── ODE verification: Cha rate ≈ ODE steady state with large RE rate constants ─

@testset "ODE verification for Cha" begin
    # Build a mechanism with RE binding step, verify against ODE with large RE rates
    m_cha = @enzyme_mechanism begin
        species: begin
            substrates: A[C]
            products:   P[C]
            enzymes:    E, EA[C]
        end
        steps: begin
            [E, A] ⇌ [EA]
            [EA] <--> [E, P]
        end
    end

    # Build equivalent full-SS mechanism for ODE simulation
    m_ss = @enzyme_mechanism begin
        species: begin
            substrates: A[C]
            products:   P[C]
            enzymes:    E, EA[C]
        end
        steps: begin
            [E, A] <--> [EA]
            [EA] <--> [E, P]
        end
    end

    rng = Random.MersenneTwister(123)
    for _ in 1:5
        K1 = 0.1 + 9.9 * rand(rng)
        k2f = 0.1 + 9.9 * rand(rng)
        k2r = 0.1 + 9.9 * rand(rng)
        Et = 0.1 + 9.9 * rand(rng)
        A = 0.1 + 9.9 * rand(rng)
        P = 0.1 + 9.9 * rand(rng)

        # Cha rate
        params_cha = (K1=K1, k2f=k2f, k2r=k2r, E_total=Et)
        concs = (A=A, P=P)
        v_cha = rate_equation(m_cha, params_cha, concs, EnzymeRates.Raw)

        # ODE with large RE rate constants: k1f = 1e6 * K1, k1r = 1e6
        # So K1 = k1f / k1r = 1e6*K1 / 1e6 = K1
        k1f_ode = 1e6 * K1
        k1r_ode = 1e6
        params_ode = (k1f=k1f_ode, k1r=k1r_ode, k2f=k2f, k2r=k2r, E_total=Et)
        v_ode = ode_steady_state_flux(m_ss, params_ode, concs)

        @test v_cha ≈ v_ode rtol=1e-3
    end
end

# ── Equilibrium test: v = 0 at Keq for mechanisms with RE steps ──────────────

@testset "Equilibrium v=0 with RE steps" begin
    m = @enzyme_mechanism begin
        species: begin
            substrates: A[C]
            products:   P[C]
            enzymes:    E, EA[C]
        end
        steps: begin
            [E, A] ⇌ [EA]
            [EA] <--> [E, P]
        end
    end

    for _ in 1:10
        K1 = 0.1 + 9.9 * rand()
        k2f = 0.1 + 9.9 * rand()
        k2r = 0.1 + 9.9 * rand()
        Et = 0.1 + 9.9 * rand()

        # Overall Keq = K1 * k2f / k2r (for A → P)
        Keq = K1 * k2f / k2r
        A_eq = 1.0
        P_eq = Keq * A_eq

        params = (K1=K1, k2f=k2f, k2r=k2r, E_total=Et)
        concs = (A=A_eq, P=P_eq)
        v = rate_equation(m, params, concs, EnzymeRates.Raw)
        @test abs(v) < 1e-10
    end
end

# ── HaldaneWegscheider mode with RE steps ────────────────────────────────────

@testset "HaldaneWegscheider with RE steps" begin
    m = @enzyme_mechanism begin
        species: begin
            substrates: A[C]
            products:   P[C]
            enzymes:    E, EA[C]
        end
        steps: begin
            [E, A] ⇌ [EA]
            [EA] <--> [E, P]
        end
    end

    @testset "Parameters" begin
        p = parameters(m, EnzymeRates.HaldaneWegscheider)
        @test :Keq in p
        @test :E_total in p
    end

    @testset "Consistency Raw vs HW" begin
        for _ in 1:20
            K1 = 0.1 + 9.9 * rand()
            k2f = 0.1 + 9.9 * rand()
            k2r = 0.1 + 9.9 * rand()
            Et = 0.1 + 9.9 * rand()
            A = 0.1 + 9.9 * rand()
            P = 0.1 + 9.9 * rand()

            raw_params = (K1=K1, k2f=k2f, k2r=k2r, E_total=Et)
            concs = (A=A, P=P)
            v_raw = rate_equation(m, raw_params, concs, EnzymeRates.Raw)

            # Build HW params: need to know which params are independent
            hw_p = parameters(m, EnzymeRates.HaldaneWegscheider)
            Keq = K1 * k2f / k2r
            # Construct HW params NamedTuple
            hw_dict = Dict{Symbol,Float64}()
            for sym in hw_p
                if sym == :Keq
                    hw_dict[sym] = Keq
                elseif sym == :E_total
                    hw_dict[sym] = Et
                elseif sym == :K1
                    hw_dict[sym] = K1
                elseif sym == :k2f
                    hw_dict[sym] = k2f
                elseif sym == :k2r
                    hw_dict[sym] = k2r
                end
            end
            hw_params = NamedTuple{Tuple(hw_p)}(Tuple(hw_dict[k] for k in hw_p))
            v_hw = rate_equation(m, hw_params, concs, EnzymeRates.HaldaneWegscheider)

            @test v_raw ≈ v_hw rtol=1e-10
        end
    end
end

# ── Rate equation string with K parameters ──────────────────────────────────

@testset "Rate equation string with RE" begin
    m = @enzyme_mechanism begin
        species: begin
            substrates: A[C]
            products:   P[C]
            enzymes:    E, EA[C]
        end
        steps: begin
            [E, A] ⇌ [EA]
            [EA] <--> [E, P]
        end
    end

    s = rate_equation_string(m, EnzymeRates.Raw)
    @test !occursin("params.", s)
    @test !occursin("concs.", s)
    @test occursin("v = E_total * (", s)
    @test occursin("K1", s)
    @test occursin("k2f", s)
end
