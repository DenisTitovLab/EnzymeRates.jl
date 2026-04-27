@testset "Types" begin
    @testset "EnzymeMechanism struct + accessors (new design)" begin
        mets = ((:S,), (:P,), ())
        rxns = (
            ((:E, :S), (:ES,), true,  1),
            ((:ES,),   (:EP,), false, 2),
            ((:EP,),   (:E, :P), true, 3),
        )
        m = EnzymeRates.EnzymeMechanism{mets, rxns}()

        @test EnzymeRates.substrates(m) == (:S,)
        @test EnzymeRates.products(m) == (:P,)
        @test EnzymeRates.regulators(m) == ()
        @test EnzymeRates.metabolites(m) == (:S, :P)
        @test EnzymeRates.reactions(m) == rxns
        @test EnzymeRates.equilibrium_steps(m) == (true, false, true)
        @test EnzymeRates.n_steps(m) == 3
        @test EnzymeRates.kinetic_group(m, 1) == 1
        @test EnzymeRates.kinetic_group(m, 2) == 2
        @test EnzymeRates.kinetic_group(m, 3) == 3
        @test EnzymeRates.kinetic_groups(m) == (1, 2, 3)
        @test EnzymeRates.steps_in_group(m, 1) == (1,)
        @test EnzymeRates.enzyme_forms(m) == (:E, :ES, :EP)
        @test EnzymeRates.n_states(m) == 3

        # Shared kinetic-group: two steps in group 1
        rxns_shared = (
            ((:E, :S),  (:ES,),  true,  1),
            ((:ES, :S), (:ESS,), true,  1),
            ((:ESS,),   (:EP,),  false, 2),
            ((:EP,),    (:E, :P), true, 3),
        )
        m2 = EnzymeRates.EnzymeMechanism{mets, rxns_shared}()
        @test EnzymeRates.kinetic_group(m2, 1) == 1
        @test EnzymeRates.kinetic_group(m2, 2) == 1
        @test EnzymeRates.steps_in_group(m2, 1) == (1, 2)
    end

    @testset "stoich_matrix has expected enzyme/metabolite rows" begin
        mets = ((:S,), (:P,), ())
        rxns = (
            ((:E, :S), (:ES,), true,  1),
            ((:ES,),   (:EP,), false, 2),
            ((:EP,),   (:E, :P), true, 3),
        )
        m = EnzymeRates.EnzymeMechanism{mets, rxns}()
        S = EnzymeRates.stoich_matrix(m)

        enz_idx = EnzymeRates.enzyme_row_range(m)
        met_idx = EnzymeRates.metabolite_row_range(m)
        @test all(sum(S[enz_idx, j]) == 0 for j in 1:size(S, 2))   # enzyme conservation
        @test size(S[met_idx, :]) == (2, 3)                          # 2 metabolites × 3 steps
    end

    @testset "EnzymeMechanism constructor" begin
        mets = ((:S,), (:P,), ())
        rxns = (
            ((:E, :S), (:ES,), true,  1),
            ((:ES,),   (:EP,), false, 2),
            ((:EP,),   (:E, :P), true, 3),
        )
        m = EnzymeRates.EnzymeMechanism(mets, rxns)
        @test EnzymeRates.n_steps(m) == 3
        @test EnzymeRates.substrates(m) == (:S,)

        # Same-kinetics group test: regulator R binds E and ES sharing one K
        mets_r = ((:S,), (:P,), (:R,))
        rxns_grouped = (
            ((:E, :S),   (:ES,),    true,  1),
            ((:ES,),     (:EP,),    false, 2),
            ((:EP,),     (:E, :P),  true,  3),
            ((:E, :R),   (:ER,),    true,  4),
            ((:ES, :R),  (:ESR,),   true,  4),
        )
        m_g = EnzymeRates.EnzymeMechanism(mets_r, rxns_grouped)
        @test EnzymeRates.kinetic_group(m_g, 4) == EnzymeRates.kinetic_group(m_g, 5)

        # Stoichiometry violation: substrate not actually consumed
        bad_rxns = (
            ((:E, :S), (:ES,), true,  1),
            ((:ES,),   (:E,),  false, 2),    # S "vanishes" — no product
        )
        @test_throws ErrorException EnzymeRates.EnzymeMechanism(((:S,), (:P,), ()), bad_rxns)

        # Iso group with size > 1 should error
        bad_iso = (
            ((:E, :S), (:ES,), true,  1),
            ((:ES,),   (:EP,), false, 99),
            ((:EP,),   (:EQ,), false, 99),    # second iso step in same group → error
            ((:EQ,),   (:E, :P), true, 2),
        )
        @test_throws ErrorException EnzymeRates.EnzymeMechanism(((:S,), (:P,), ()), bad_iso)
    end

    @testset "AllostericEnzymeMechanism (new design)" begin
        cm = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                [E, S] ⇌ [ES]
                [ES] <--> [EP]
                [EP] ⇌ [E, P]
            end
        end
        cat_sites = (2, ((2, :OnlyR),))
        reg_sites = ((((:I,), 2, ((:I, :OnlyT),)),),)
        m = EnzymeRates.AllostericEnzymeMechanism{typeof(cm), cat_sites, reg_sites[1]}()

        @test EnzymeRates.catalytic_mechanism(m) === cm
        @test EnzymeRates.catalytic_multiplicity(m) == 2
        @test EnzymeRates.group_tag(m, 1) == :NonequalRT   # default
        @test EnzymeRates.group_tag(m, 2) == :OnlyR
        @test EnzymeRates.regulatory_sites(m) == reg_sites[1]
    end

    @testset "AllostericEnzymeMechanism constructor + DSL" begin
        cm = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                [E, S] ⇌ [ES]
                [ES] <--> [EP]
                [EP] ⇌ [E, P]
            end
        end

        # Single-ligand :EqualRT reg site → error
        @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, ()), ((((:I,), 2, ((:I, :EqualRT),)),),)[1],
        )

        # Iso group :OnlyT → error
        @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, ((2, :OnlyT),)), (),
        )

        # Build via DSL
        m = @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: I::OnlyT

            site(:catalytic, 2): begin
                steps: begin
                    [E, S] ⇌ [ES]    :: EqualRT
                    [ES] <--> [EP]    :: OnlyR
                    [EP] ⇌ [E, P]    :: EqualRT
                end
            end
        end
        @test EnzymeRates.catalytic_multiplicity(m) == 2
        @test EnzymeRates.group_tag(m, 1) == :EqualRT
        @test EnzymeRates.group_tag(m, 2) == :OnlyR
        @test EnzymeRates.allosteric_regulators(m) == ((:I, :OnlyT),)
    end

    @testset "EnzymeReaction" begin
        r = EnzymeReaction(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
        )
        @test r isa EnzymeReaction
        @test EnzymeRates.substrates(r) == ((:S, ((:C, 1),)),)
        @test EnzymeRates.products(r) == ((:P, ((:C, 1),)),)
        @test EnzymeRates.regulators(r) == ()
    end

    @testset "EnzymeReaction with regulators" begin
        r = EnzymeReaction(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            (:I,),
        )
        @test EnzymeRates.regulators(r) == (:I,)
    end

    @testset "EnzymeReaction with regulator roles" begin
        r = EnzymeReaction(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            ((:I, :dead_end), (:A, :allosteric)),
        )
        @test EnzymeRates.regulators(r) == (:A, :I)
        roles = EnzymeRates.regulator_roles(r)
        @test length(roles) == 2
        role_dict = Dict(r[1] => r[2] for r in roles)
        @test role_dict[:I] == :dead_end
        @test role_dict[:A] == :allosteric
    end

    @testset "EnzymeReaction canonical ordering" begin
        r1 = EnzymeReaction(
            ((:S2, ((:C, 2),)), (:S1, ((:C, 1),))),
            ((:P1, ((:C, 1),)), (:P2, ((:C, 2),))),
        )
        r2 = EnzymeReaction(
            ((:S1, ((:C, 1),)), (:S2, ((:C, 2),))),
            ((:P1, ((:C, 1),)), (:P2, ((:C, 2),))),
        )
        @test typeof(r1) === typeof(r2)
        @test EnzymeRates.substrates(r1) == ((:S1, ((:C, 1),)), (:S2, ((:C, 2),)))
    end

    @testset "EnzymeReaction regulator same as substrate allowed" begin
        r = EnzymeReaction(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            (:S,),  # regulator same as substrate - allowed
        )
        @test r isa EnzymeReaction
    end

    @testset "EnzymeReaction regulator same as product allowed" begin
        r = EnzymeReaction(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            (:P,),  # regulator same as product - allowed
        )
        @test r isa EnzymeReaction
    end

    @testset "EnzymeReaction duplicate substrate names" begin
        @test_throws ErrorException EnzymeReaction(
            ((:S, ((:C, 1),)), (:S, ((:C, 1),))),
            ((:P, ((:C, 1),)),),
        )
    end

    @testset "EnzymeReaction duplicate product names" begin
        @test_throws ErrorException EnzymeReaction(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)), (:P, ((:C, 1),))),
        )
    end

    @testset "EnzymeReaction oligomeric_state" begin
        # Default oligomeric_state is 1
        rxn = EnzymeReaction(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
        )
        @test EnzymeRates.oligomeric_state(rxn) == 1

        # Explicit oligomeric_state
        rxn2 = EnzymeReaction(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            ();
            oligomeric_state=4
        )
        @test EnzymeRates.oligomeric_state(rxn2) == 4

        # Different oligomeric_state = different type
        @test typeof(rxn) !== typeof(rxn2)
    end

    @testset "Pretty printing" begin
        r = EnzymeReaction(((:S, ((:C, 1),)),), ((:P, ((:C, 1),)),))
        @test sprint(show, r) == "EnzymeReaction: S ⇌ P"

        r2 = EnzymeReaction(
            ((:S, ((:C, 1),)), (:ATP, ((:C, 10),))),
            ((:P, ((:C, 1),)), (:ADP, ((:C, 10),))),
            (:I,),
        )
        @test sprint(show, r2) == "EnzymeReaction: ATP + S ⇌ ADP + P | regulators: I"

        # Linear mechanism: compact chain
        species = (
            ((:S, ((:C, 1),)),), ((:P, ((:C, 1),)),), (),
            ((:E, ()), (:ES, ((:C, 1),))),
        )
        m = EnzymeMechanism(
            species,
            (((:E, :S), (:ES,)), ((:ES,), (:E, :P))),
            (false, false),
        )
        @test sprint(show, m) == "EnzymeMechanism: E + S <--> ES <--> E + P"

        # Branched mechanism: multi-line
        species_b = (
            ((:S1, ((:C, 1),)), (:S2, ((:H, 1),))),
            ((:P1, ((:C, 1),)), (:P2, ((:H, 1),))),
            (),
            ((:E, ()), (:ES1, ((:C, 1),)), (:ES2, ((:H, 1),)),
             (:ES1S2, ((:C, 1), (:H, 1))), (:EP1P2, ((:C, 1), (:H, 1))),
             (:EP2, ((:H, 1),))),
        )
        rxns_b = (
            ((:E, :S1), (:ES1,)), ((:E, :S2), (:ES2,)),
            ((:ES1, :S2), (:ES1S2,)), ((:ES2, :S1), (:ES1S2,)),
            ((:ES1S2,), (:EP1P2,)), ((:EP1P2,), (:EP2, :P1)), ((:EP2,), (:E, :P2)),
        )
        m_b = EnzymeMechanism(species_b, rxns_b, ntuple(Returns(false), length(rxns_b)))
        str = sprint(show, m_b)
        @test startswith(str, "EnzymeMechanism (7 steps, 6 enzyme forms):")
        @test contains(str, "E + S1 <--> ES1")
        @test contains(str, "EP2 <--> E + P2")
    end

    @testset "EnzymeMechanism different orderings produce valid mechanisms" begin
        mets = ((:S,), (:P,), ())
        rxns1 = (
            ((:E, :S), (:ES,), false, 1),
            ((:ES,), (:E, :P), false, 2),
        )
        rxns2 = (
            ((:ES,), (:E, :P), false, 1),
            ((:E, :S), (:ES,), false, 2),
        )
        m1 = EnzymeMechanism(mets, rxns1)
        m2 = EnzymeMechanism(mets, rxns2)
        @test m1 isa EnzymeMechanism
        @test m2 isa EnzymeMechanism
    end

    @testset "EnzymeMechanism error cases" begin
        base_mets = ((:S,), (:P,), ())
        base_rxns = (
            ((:E, :S), (:ES,), false, 1),
            ((:ES,), (:E, :P), false, 2),
        )

        # Empty reactions tuple
        @test_throws ErrorException EnzymeMechanism(base_mets, ())

        # Duplicate substrate names in metabolites
        @test_throws ErrorException EnzymeMechanism(
            ((:S, :S), (:P,), ()), base_rxns)

        # Reaction with zero enzymes on LHS (only metabolite)
        no_enz_rxns = (
            ((:S,), (:ES,), false, 1),
            ((:ES,), (:E, :P), false, 2),
        )
        @test_throws ErrorException EnzymeMechanism(base_mets, no_enz_rxns)

        # Reaction with two metabolites on one side
        two_met_mets = ((:S1, :S2), (:P,), ())
        two_met_rxns = (
            ((:E, :S1, :S2), (:ES,), false, 1),
            ((:ES,), (:E, :P), false, 2),
        )
        @test_throws ErrorException EnzymeMechanism(two_met_mets, two_met_rxns)

        # Unknown metabolite in reaction
        unknown_rxns = (
            ((:E, :X), (:ES,), false, 1),
            ((:ES,), (:E, :P), false, 2),
        )
        @test_throws ErrorException EnzymeMechanism(base_mets, unknown_rxns)

        # Net stoichiometry mismatch (substrate consumed but never produced)
        # (Substrate listed but doesn't appear in any step.)
        net_mismatch_mets = ((:S, :S2), (:P,), ())
        @test_throws ErrorException EnzymeMechanism(net_mismatch_mets, base_rxns)

        # Duplicate reactions
        dup_rxns = (
            ((:E, :S), (:ES,), false, 1),
            ((:ES,), (:E, :P), false, 2),
            ((:E, :S), (:ES,), false, 3),
        )
        @test_throws ErrorException EnzymeMechanism(base_mets, dup_rxns)

        # Unreachable enzyme form (EX referenced but never in a reaction)
        unreachable_rxns = (
            ((:E, :S), (:ES,), false, 1),
            ((:ES,), (:E, :P), false, 2),
            ((:EX,), (:E,), false, 3),  # EX has no producing step
        )
        @test_throws ErrorException EnzymeMechanism(base_mets, unreachable_rxns)
    end

    @testset "EnzymeMechanism valid with reachable enzyme forms" begin
        mets = ((:S,), (:P,), ())
        rxns = (
            ((:E, :S), (:ES,), false, 1),
            ((:ES,), (:E, :P), false, 2),
        )
        m = EnzymeMechanism(mets, rxns)
        @test m isa EnzymeMechanism
    end

    # The "Constraint constructor validation errors" testset was deleted —
    # it tested the OLD 4-arg `EnzymeMechanism(species, rxns, eq, constraints)`
    # constructor's `param_constraints`-symbol validators (Invalid target,
    # Self-reference, Duplicate target, Invalid replacement, Zero coefficient).
    # The new design has no symbol-level constraints; kinetic groups are Ints.

    # TODO: Phase 3 — testsets below use cross-metabolite constraints
    # (`K2 = K1` between A-binding and B-binding, `k3r = k1r` between
    # A-binding and P-binding kinetics) that the new step-group syntax
    # cannot express. Restore once richer constraint expressions return.
    @static if false
    @testset "parameters() excludes constrained params" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C]
                products:   P[C]
                enzymes:    E, EA[C], EP[C]
            end
            steps: begin
                [E, A] <--> [EA]
                [EA] <--> [EP]
                [EP] <--> [E, P]
            end
            constraints: begin
                k3r = k1r
            end
        end
        raw_params = parameters(m, Full)
        @test :k3r ∉ raw_params
        @test :k1r ∈ raw_params
        @test :E_total ∈ raw_params
    end

    @testset "Constrained Uni-Uni k3r = k1r: rate equation correctness" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C]
                products:   P[C]
                enzymes:    E, EA[C], EP[C]
            end
            steps: begin
                [E, A] <--> [EA]
                [EA] <--> [EP]
                [EP] <--> [E, P]
            end
            constraints: begin
                k3r = k1r
            end
        end

        function rate_constrained(p, c)
            (; k1f, k1r, k2f, k2r, k3f, E_total) = p
            (; A, P) = c
            k3r = k1r
            num = k1f * k2f * k3f * A - k1r * k2r * k3r * P
            denom = (k1r * k2r + k1r * k3f + k2f * k3f) +
                    k1f * (k2r + k2f + k3f) * A +
                    k3r * (k1r + k2r + k2f) * P
            E_total * num / denom
        end

        rng = Random.MersenneTwister(42)
        @test all(1:20) do _
            params = (k1f = 0.1 + 9.9 * rand(rng),
                      k1r = 0.1 + 9.9 * rand(rng),
                      k2f = 0.1 + 9.9 * rand(rng),
                      k2r = 0.1 + 9.9 * rand(rng),
                      k3f = 0.1 + 9.9 * rand(rng),
                      E_total = 0.1 + 9.9 * rand(rng))
            concs = (A = 0.1 + 9.9 * rand(rng), P = 0.1 + 9.9 * rand(rng))
            isapprox(
                rate_equation(m, concs, params, Full),
                rate_constrained(params, concs);
                rtol=1e-10)
        end
    end

    @testset "Constrained RE K2 = K1: rate equation correctness" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products:   P[C], Q[N]
                enzymes:    E, EA[C], EABEPQ[CN], EQ[N]
            end
            steps: begin
                [E, A] ⇌ [EA]
                [EA, B] ⇌ [EABEPQ]
                [EABEPQ] <--> [EQ, P]
                [EQ] <--> [E, Q]
            end
            constraints: begin
                K2 = K1
            end
        end

        function rate_constrained_re(p, c)
            (; K1, k3f, k3r, k4f, k4r, E_total) = p
            (; A, B, P, Q) = c
            num = k3f * k4f * A * B / (K1 * K1) - k3r * k4r * P * Q
            sigma1 = 1.0 + A / K1 + A * B / (K1 * K1)
            R12 = k3f * A * B / (K1 * K1) + k4r * Q
            R21 = k3r * P + k4f
            denom = sigma1 * R21 + R12
            E_total * num / denom
        end

        raw_params = parameters(m, Full)
        @test :K2 ∉ raw_params
        @test :K1 ∈ raw_params

        rng = Random.MersenneTwister(99)
        @test all(1:20) do _
            params = (K1 = 0.1 + 9.9 * rand(rng),
                      k3f = 0.1 + 9.9 * rand(rng),
                      k3r = 0.1 + 9.9 * rand(rng),
                      k4f = 0.1 + 9.9 * rand(rng),
                      k4r = 0.1 + 9.9 * rand(rng),
                      E_total = 0.1 + 9.9 * rand(rng))
            concs = (A = 0.1 + 9.9 * rand(rng), B = 0.1 + 9.9 * rand(rng),
                     P = 0.1 + 9.9 * rand(rng), Q = 0.1 + 9.9 * rand(rng))
            isapprox(
                rate_equation(m, concs, params, Full),
                rate_constrained_re(params, concs);
                rtol=1e-10)
        end
    end

    @testset "Constrained equilibrium (v=0 at Keq)" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C]
                products:   P[C]
                enzymes:    E, EA[C], EP[C]
            end
            steps: begin
                [E, A] <--> [EA]
                [EA] <--> [EP]
                [EP] <--> [E, P]
            end
            constraints: begin
                k3r = k1r
            end
        end

        hw_params = parameters(m)
        @test :k3r ∉ hw_params
        @test :Keq ∈ hw_params

        rng = Random.MersenneTwister(77)
        hw = parameters(m)
        vals = Tuple(0.1 + 9.9 * rand(rng) for _ in hw)
        p = NamedTuple{hw}(vals)
        Keq = p.Keq

        eq_concs = (A = 1.0, P = Keq)
        v = rate_equation(m, eq_concs, p)
        @test abs(v) < 1e-10
    end

    @testset "Constrained zero allocation" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C]
                products:   P[C]
                enzymes:    E, EA[C], EP[C]
            end
            steps: begin
                [E, A] <--> [EA]
                [EA] <--> [EP]
                [EP] <--> [E, P]
            end
            constraints: begin
                k3r = k1r
            end
        end
        rng = Random.MersenneTwister(55)
        hw = parameters(m)
        vals = Tuple(0.1 + 9.9 * rand(rng) for _ in hw)
        p = NamedTuple{hw}(vals)
        concs = (A = 1.0, P = 2.0)

        rate_equation(m, concs, p)  # warmup
        allocs = @allocated rate_equation(m, concs, p)
        @test allocs == 0
    end

    @testset "rate_equation_string shows constraints" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C]
                products:   P[C]
                enzymes:    E, EA[C], EP[C]
            end
            steps: begin
                [E, A] <--> [EA]
                [EA] <--> [EP]
                [EP] <--> [E, P]
            end
            constraints: begin
                k3r = k1r
            end
        end

        s_raw = rate_equation_string(m, Full)
        @test occursin("k3r = k1r", s_raw)
        @test occursin("v = E_total * (", s_raw)
        @test !occursin("k3r", replace(s_raw, "k3r = k1r" => ""))

        s_hw = rate_equation_string(m)
        @test occursin("k3r = k1r", s_hw)
    end

    @testset "show method displays constraints" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C]
                products:   P[C]
                enzymes:    E, EA[C], EP[C]
            end
            steps: begin
                [E, A] <--> [EA]
                [EA] <--> [EP]
                [EP] <--> [E, P]
            end
            constraints: begin
                k3r = k1r
            end
        end
        s = sprint(show, m)
        @test occursin("constraints:", s)
        @test occursin("k3r = k1r", s)
    end
    end # @static if false (Phase 3: cross-metabolite constraints)
end
