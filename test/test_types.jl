@testset "Types" begin
    @testset "EnzymeMechanism struct + accessors (new design)" begin
        mets = ((:S,), (:P,), ())
        rxns = (
            ((:E, :S), (:ES,), true,  1),
            ((:ES,),   (:EP,), false, 2),
            ((:EP,),   (:E, :P), true, 3),
        )
        m = EnzymeRates.EnzymeMechanism{(mets, rxns)}()

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
        m2 = EnzymeRates.EnzymeMechanism{(mets, rxns_shared)}()
        @test EnzymeRates.kinetic_group(m2, 1) == 1
        @test EnzymeRates.kinetic_group(m2, 2) == 1
        @test EnzymeRates.steps_in_group(m2, 1) == (1, 2)
    end

    @testset "EnzymeMechanism Sig repack" begin
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
            end
        end

        @test m isa EnzymeRates.EnzymeMechanism
        Sig = typeof(m).parameters[1]
        @test Sig isa Tuple
        @test length(Sig) == 2

        @test EnzymeRates.substrates(m) == (:S,)
        @test EnzymeRates.products(m)   == (:P,)
        @test EnzymeRates.n_steps(m)    == 3
    end

    @testset "stoich_matrix has expected enzyme/metabolite rows" begin
        mets = ((:S,), (:P,), ())
        rxns = (
            ((:E, :S), (:ES,), true,  1),
            ((:ES,),   (:EP,), false, 2),
            ((:EP,),   (:E, :P), true, 3),
        )
        m = EnzymeRates.EnzymeMechanism{(mets, rxns)}()
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
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
            end
        end
        # Dense format: (multiplicity, cat_allo_states) and (ligands, mult, reg_allo_states)
        cat_sites = (2, (:NonequalRT, :OnlyR, :NonequalRT))
        reg_sites = ((((:I,), 2, (:OnlyT,)),),)
        m = EnzymeRates.AllostericEnzymeMechanism{typeof(cm), cat_sites, reg_sites[1]}()

        @test EnzymeRates.catalytic_mechanism(m) === cm
        @test EnzymeRates.catalytic_multiplicity(m) == 2
        @test EnzymeRates.cat_allo_state(m, 1) == :NonequalRT
        @test EnzymeRates.cat_allo_state(m, 2) == :OnlyR
        @test EnzymeRates.regulatory_sites(m) == reg_sites[1]
    end

    @testset "AllostericEnzymeMechanism constructor + DSL" begin
        cm = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
            end
        end

        # Single-ligand :EqualRT reg site → error
        @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:NonequalRT, :NonequalRT, :NonequalRT)),
            (((:I,), 2, (:EqualRT,)),),
        )

        # Catalytic group :OnlyT → error
        @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:NonequalRT, :OnlyT, :NonequalRT)), (),
        )

        # Build via DSL
        m = @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: I::OnlyT

            site(:catalytic, 2): begin
                steps: begin
                    E + S ⇌ ES    :: EqualRT
                    ES <--> EP    :: OnlyR
                    EP ⇌ E + P    :: EqualRT
                end
            end
        end
        @test EnzymeRates.catalytic_multiplicity(m) == 2
        @test EnzymeRates.cat_allo_state(m, 1) == :EqualRT
        @test EnzymeRates.cat_allo_state(m, 2) == :OnlyR
        @test EnzymeRates.allosteric_regulators(m) == ((:I, :OnlyT),)
    end

    @testset "EnzymeReactionLegacy" begin
        r = EnzymeReactionLegacy(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
        )
        @test r isa EnzymeReactionLegacy
        @test EnzymeRates.substrates(r) == ((:S, ((:C, 1),)),)
        @test EnzymeRates.products(r) == ((:P, ((:C, 1),)),)
        @test EnzymeRates.regulators(r) == ()
    end

    @testset "EnzymeReactionLegacy with regulators" begin
        r = EnzymeReactionLegacy(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            (:I,),
        )
        @test EnzymeRates.regulators(r) == (:I,)
    end

    @testset "EnzymeReactionLegacy with regulator roles" begin
        r = EnzymeReactionLegacy(
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

    @testset "EnzymeReactionLegacy canonical ordering" begin
        r1 = EnzymeReactionLegacy(
            ((:S2, ((:C, 2),)), (:S1, ((:C, 1),))),
            ((:P1, ((:C, 1),)), (:P2, ((:C, 2),))),
        )
        r2 = EnzymeReactionLegacy(
            ((:S1, ((:C, 1),)), (:S2, ((:C, 2),))),
            ((:P1, ((:C, 1),)), (:P2, ((:C, 2),))),
        )
        @test typeof(r1) === typeof(r2)
        @test EnzymeRates.substrates(r1) == ((:S1, ((:C, 1),)), (:S2, ((:C, 2),)))
    end

    @testset "EnzymeReactionLegacy regulator same as substrate allowed" begin
        r = EnzymeReactionLegacy(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            (:S,),  # regulator same as substrate - allowed
        )
        @test r isa EnzymeReactionLegacy
    end

    @testset "EnzymeReactionLegacy regulator same as product allowed" begin
        r = EnzymeReactionLegacy(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            (:P,),  # regulator same as product - allowed
        )
        @test r isa EnzymeReactionLegacy
    end

    @testset "EnzymeReactionLegacy duplicate substrate names" begin
        @test_throws ErrorException EnzymeReactionLegacy(
            ((:S, ((:C, 1),)), (:S, ((:C, 1),))),
            ((:P, ((:C, 1),)),),
        )
    end

    @testset "EnzymeReactionLegacy duplicate product names" begin
        @test_throws ErrorException EnzymeReactionLegacy(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)), (:P, ((:C, 1),))),
        )
    end

    @testset "EnzymeReactionLegacy oligomeric_state" begin
        # Default oligomeric_state is 1
        rxn = EnzymeReactionLegacy(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
        )
        @test EnzymeRates.oligomeric_state(rxn) == 1

        # Explicit oligomeric_state
        rxn2 = EnzymeReactionLegacy(
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
        r = EnzymeReactionLegacy(((:S, ((:C, 1),)),), ((:P, ((:C, 1),)),))
        @test sprint(show, r) == "EnzymeReaction: S ⇌ P"

        r2 = EnzymeReactionLegacy(
            ((:S, ((:C, 1),)), (:ATP, ((:C, 10),))),
            ((:P, ((:C, 1),)), (:ADP, ((:C, 10),))),
            (:I,),
        )
        @test sprint(show, r2) == "EnzymeReaction: ATP + S ⇌ ADP + P | regulators: I"

        # Linear mechanism: compact chain.
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S <--> ES
                ES <--> E + P
            end
        end
        @test sprint(show, m) ==
            "EnzymeMechanism: E + S <--> ES <--> E + P"

        # Branched mechanism: multi-line with header summary.
        m_b = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                E + A <--> EA
                E + B <--> EB
                EA + B <--> EAB
                EB + A <--> EAB
                EAB <--> EPQ
                EPQ <--> EQ + P
                EQ <--> E + Q
            end
        end
        s = sprint(show, m_b)
        @test startswith(s, "EnzymeMechanism (7 steps, 6 enzyme forms):")
        @test contains(s, "E + A <--> EA")
        @test contains(s, "EQ <--> E + Q")

        # Linear chain with canonical RE binding: chain-walk renders the
        # release step in reverse so the chain stays linear.
        m_re = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
            end
        end
        @test sprint(show, m_re) ==
            "EnzymeMechanism: E + S ⇌ ES <--> EP ⇌ E + P"

        # Mechanism with regulators: appended at end.
        m_reg = @enzyme_mechanism begin
            substrates: S
            products:   P
            regulators: I
            steps: begin
                E + S <--> ES
                ES <--> E + P
                E + I <--> EI
            end
        end
        @test contains(sprint(show, m_reg), "| regulators: I")

        # EnzymeReactionLegacy with oligomeric_state > 1.
        rxn_oligo = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            oligomeric_state: 4
        end
        @test sprint(show, rxn_oligo) ==
            "EnzymeReaction: S ⇌ P | oligomeric_state: 4"

        # AllostericEnzymeMechanism (smoke).
        m_allo = @allosteric_mechanism begin
            substrates: F6P
            products:   F16BP
            allosteric_regulators: I::OnlyT
            site(:catalytic, 2): begin
                steps: begin
                    E + F6P ⇌ E_F6P   :: EqualRT
                    E_F6P <--> E_F16BP :: EqualRT
                    E_F16BP ⇌ E + F16BP :: EqualRT
                end
            end
        end
        s_allo = sprint(show, m_allo)
        @test contains(s_allo, "AllostericEnzymeMechanism (cat_n=2")
        @test contains(s_allo, "reg sites")
        @test contains(s_allo, "I::OnlyT")
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

        # NOTE: "Duplicate reactions" and "Unreachable enzyme form" tests
        # were dropped — the new design accepts both. Two reactions with
        # the same (lhs, rhs) but distinct kinetic_groups are valid (they
        # represent dead-end mirrors with different parameters in the OLD
        # design's terms; in the new design, the kinetic_group integer
        # disambiguates). And the new constructor doesn't enforce a
        # connectivity invariant — enzyme forms are inferred from steps,
        # so an "unreachable" form simply has its own steps in isolation,
        # which is structurally valid (graph connectivity is a downstream
        # concern caught by Wegscheider analysis if it matters).
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

    @testset "Kinetic-group validator error paths" begin
        mets = ((:S, :A), (:P,), ())

        # Group binding different metabolites → error
        @test_throws ErrorException EnzymeMechanism(mets, (
            ((:E, :S), (:ES,), true, 1),
            ((:E, :A), (:EA,), true, 1),     # group 1 also binds A
            ((:EA,), (:E, :P), false, 2),
        ))

        # Group mixing RE and SS → error
        @test_throws ErrorException EnzymeMechanism(((:S,), (:P,), ()), (
            ((:E, :S), (:ES,), true, 1),     # RE
            ((:E, :S), (:ES,), false, 1),    # SS — same group as RE step
            ((:ES,), (:E, :P), true, 2),
        ))
    end

    @testset "Connectivity validator: orphan enzyme form → error" begin
        # Two disjoint subgraphs: E↔ES (catalytic) and EX↔EY (orphan).
        @test_throws ErrorException EnzymeMechanism(((:S,), (:P,), ()), (
            ((:E, :S), (:ES,), true, 1),
            ((:ES,), (:E, :P), true, 2),
            ((:EX,), (:EY,), false, 3),
            ((:EY,), (:EX,), false, 4),
        ))
    end

    @testset "AllostericEnzymeMechanism constructor validators" begin
        cm = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
            end
        end
        # Wrong-length cat_allo_states (4 entries for 3 kinetic groups) → error
        @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:NonequalRT, :NonequalRT, :NonequalRT, :OnlyR)), ())

        # Invalid allo state value → error
        @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:NotAState, :NonequalRT, :NonequalRT)), ())

        # Reg site with no ligands → error
        @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:NonequalRT, :NonequalRT, :NonequalRT)), (((), 2, ()),))

        # Reg site with all-:EqualRT ligands → error (cancels identically)
        @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:NonequalRT, :NonequalRT, :NonequalRT)),
            (((:I, :J), 2, (:EqualRT, :EqualRT)),))

        # Invalid reg-site ligand allo state → error
        @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:NonequalRT, :NonequalRT, :NonequalRT)),
            (((:I,), 2, (:NotAState,)),))

        # Non-consecutive kinetic_group numbers (1, 3 instead of 1, 2) → error
        # The cat_allo_states tuple is indexed by group number, so a hole
        # would cause OOB or wrong-state lookup at runtime.
        bad_mets = ((:S,), (:P,), ())
        bad_rxns = (
            ((:E, :S),  (:ES,),    true,  1),
            ((:ES,),    (:EP,),    false, 3),
            ((:E, :P),  (:EP,),    true,  1),
        )
        cm_bad = EnzymeRates.EnzymeMechanism{(bad_mets, bad_rxns)}()
        @test_throws ErrorException EnzymeRates.AllostericEnzymeMechanism(
            cm_bad, (2, (:NonequalRT, :NonequalRT)), ())
    end

    @testset "Base.show displays all dense states" begin
        m = @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: I::NonequalRT, J::OnlyT
            site(:catalytic, 2): begin
                steps: begin
                    E + S ⇌ ES    :: NonequalRT
                    ES <--> EP   :: OnlyR
                    EP ⇌ E + P   :: EqualRT
                end
            end
        end
        s = sprint(show, m)
        # Every catalytic state appears in cat_allo_states line
        cm_inner = EnzymeRates.catalytic_mechanism(m)
        n_groups = length(unique(EnzymeRates.kinetic_group(cm_inner, i)
                                 for i in 1:EnzymeRates.n_steps(cm_inner)))
        for g in 1:n_groups
            @test occursin(string(EnzymeRates.cat_allo_state(m, g)), s)
        end
        # No :NonequalRT ligand silently hidden from reg-site display
        for (i, _) in enumerate(EnzymeRates.regulatory_sites(m))
            for lig in EnzymeRates.regulatory_site_ligands(m, i)
                state = EnzymeRates.reg_allo_state(m, i, lig)
                @test occursin("$lig::$state", s)
            end
        end
    end

    @testset "EnzymeReactionLegacy: atom mandatory" begin
        @test_throws ErrorException EnzymeReactionLegacy(
            ((:S, ()),),
            ((:P, ((:C, 1),)),)
        )
        @test_throws ErrorException EnzymeReactionLegacy(
            ((:S, ((:C, 1),)),),
            ((:P, ()),)
        )
        @test EnzymeReactionLegacy(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),)
        ) isa EnzymeReactionLegacy
    end

    @testset "EnzymeReactionLegacy: atom balance" begin
        @test_throws ErrorException EnzymeReactionLegacy(
            ((:S, ((:C, 6),)),),
            ((:P, ((:C, 5),)),)
        )
        @test_throws ErrorException EnzymeReactionLegacy(
            ((:S, ((:C, 6), (:H, 12))),),
            ((:P, ((:C, 6),)),)
        )
        @test EnzymeReactionLegacy(
            ((:A, ((:C, 6),)), (:B, ((:N, 1),))),
            ((:P, ((:C, 6),)), (:Q, ((:N, 1),)))
        ) isa EnzymeReactionLegacy
    end

    @testset "EnzymeMechanism: strict regulator binding" begin
        # Regulator :A listed but never bound in any step -> error
        @test_throws ErrorException EnzymeRates.EnzymeMechanism(
            ((:S,), (:P,), (:A,)),
            (((:E, :S), (:E_S,), true, 1),
             ((:E_S,), (:E_P,), false, 2),
             ((:E_P,), (:E, :P), true, 3))
        )
        # All regulators bound -> ok
        @test EnzymeRates.EnzymeMechanism(
            ((:S,), (:P,), (:A,)),
            (((:E, :S), (:E_S,), true, 1),
             ((:E_S, :A), (:E_S_A,), true, 4),
             ((:E_S,), (:E_P,), false, 2),
             ((:E_P,), (:E, :P), true, 3))
        ) isa EnzymeRates.EnzymeMechanism
        # No regulators -> ok
        @test EnzymeRates.EnzymeMechanism(
            ((:S,), (:P,), ()),
            (((:E, :S), (:E_S,), true, 1),
             ((:E_S,), (:E_P,), false, 2),
             ((:E_P,), (:E, :P), true, 3))
        ) isa EnzymeRates.EnzymeMechanism
    end

    @testset "AllostericEnzymeMechanism display format" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: R
            oligomeric_state: 2
        end
        init = EnzymeRates.init_mechanisms(rxn)
        base = first(init)
        used_groups = sort!(collect(
            Set(s.kinetic_group for s in base.steps)))
        g_s = first(s.kinetic_group for s in base.steps
                    if EnzymeRates.step_metabolite(s) === :S)
        g_p = first(s.kinetic_group for s in base.steps
                    if EnzymeRates.step_metabolite(s) === :P)
        group_tags = Dict{Int,Symbol}(
            g => :NonequalRT for g in used_groups)
        group_tags[g_s] = :EqualRT
        group_tags[g_p] = :EqualRT
        spec = EnzymeRates.AllostericMechanismSpec(
            base, 2, [[:R]], [2],
            group_tags,
            Dict(:R => :OnlyT),
            base.n_fit_params_estimate + 1)
        m = EnzymeRates.AllostericEnzymeMechanism(spec)
        s = repr(m)

        # Old summary line gone:
        @test !occursin("cat_allo_states:", s)
        # Inline ::Tag annotations on each step or step group:
        @test occursin(":: EqualRT", s)
        # Multi-line catalytic display (no chain shortcut):
        n_steps_re = count(c -> c == '\n', s)
        @test n_steps_re >= 3   # header + ≥3 step lines
    end

    @testset "AllostericEnzymeMechanism display: shared kinetic group" begin
        cm = EnzymeMechanism(
            ((:S,), (:P,), ()),
            (((:E, :S),    (:E_S,),  true,  1),
             ((:E_P, :S),  (:E_PS,), true,  1),
             ((:E_S,),     (:E_P,),  false, 2),
             ((:E_P,),     (:E, :P), true,  3)))
        am = EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:EqualRT, :EqualRT, :EqualRT)), ())
        s = repr(am)

        # The parenthesized-group branch must execute: look for the
        # exact "(...) :: EqualRT" shape with a comma-separated body.
        paren_group_match = match(
            r"\([^()]*,[^()]*\) :: EqualRT", s)
        @test paren_group_match !== nothing
        # Format: lhs/rhs joined with " + " (no brackets).
        @test occursin("E_S <--> E_P :: EqualRT", s)
        @test occursin(":: EqualRT", s)
        @test !occursin("cat_allo_states:", s)
    end

    # ─── Concrete type hierarchy (spec §5.1–5.7) ──────────────────────

    @testset "Metabolite hierarchy: Substrate / Product / Regulators" begin
        s = EnzymeRates.Substrate(:ATP)
        p = EnzymeRates.Product(:ADP)
        a = EnzymeRates.AllostericRegulator(:cAMP)
        i = EnzymeRates.CompetitiveInhibitor(:I)

        @test s isa EnzymeRates.Substrate
        @test s isa EnzymeRates.Reactant
        @test s isa EnzymeRates.Metabolite
        @test p isa EnzymeRates.Product
        @test p isa EnzymeRates.Reactant
        @test a isa EnzymeRates.AllostericRegulator
        @test a isa EnzymeRates.Regulator
        @test a isa EnzymeRates.Metabolite
        @test i isa EnzymeRates.CompetitiveInhibitor
        @test i isa EnzymeRates.Regulator

        @test EnzymeRates.name(s) === :ATP
        @test EnzymeRates.name(p) === :ADP
        @test EnzymeRates.name(a) === :cAMP
        @test EnzymeRates.name(i) === :I

        @test EnzymeRates.Substrate(:X) == EnzymeRates.Substrate(:X)
        @test EnzymeRates.Substrate(:X) != EnzymeRates.Substrate(:Y)
        # Distinct subtypes with same name are NOT equal (struct identity matters).
        @test EnzymeRates.Substrate(:X) != EnzymeRates.Product(:X)
        @test hash(EnzymeRates.Substrate(:X)) == hash(EnzymeRates.Substrate(:X))
    end

    @testset "Residual: empty default + canonical ordering" begin
        empty_r = EnzymeRates.Residual()
        @test isempty(empty_r)
        @test EnzymeRates.added(empty_r) == EnzymeRates.Substrate[]
        @test EnzymeRates.subtracted(empty_r) == EnzymeRates.Product[]

        r1 = EnzymeRates.Residual(
            [EnzymeRates.Substrate(:B), EnzymeRates.Substrate(:A)],
            [EnzymeRates.Product(:Q), EnzymeRates.Product(:P)],
        )
        r2 = EnzymeRates.Residual(
            [EnzymeRates.Substrate(:A), EnzymeRates.Substrate(:B)],
            [EnzymeRates.Product(:P), EnzymeRates.Product(:Q)],
        )
        @test r1 == r2
        @test hash(r1) == hash(r2)
        @test !isempty(r1)
        @test EnzymeRates.added(r1) ==
              [EnzymeRates.Substrate(:A), EnzymeRates.Substrate(:B)]
        @test EnzymeRates.subtracted(r1) ==
              [EnzymeRates.Product(:P), EnzymeRates.Product(:Q)]
    end

    @testset "Species: canonical bound ordering + accessors + name" begin
        s1 = EnzymeRates.Species(
            EnzymeRates.Metabolite[
                EnzymeRates.Substrate(:B), EnzymeRates.Substrate(:A)],
            :E,
        )
        s2 = EnzymeRates.Species(
            EnzymeRates.Metabolite[
                EnzymeRates.Substrate(:A), EnzymeRates.Substrate(:B)],
            :E,
        )
        @test s1 == s2
        @test hash(s1) == hash(s2)
        @test EnzymeRates.conformation(s1) === :E
        @test EnzymeRates.residual(s1) == EnzymeRates.Residual()
        @test !EnzymeRates.has_residual(s1)
        @test EnzymeRates.bound(s1) ==
              EnzymeRates.Metabolite[
                  EnzymeRates.Substrate(:A), EnzymeRates.Substrate(:B)]

        # Empty bound, :E conformation → :E
        s_e = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        @test EnzymeRates.name(s_e) === :E

        # Bound metabolites are appended in canonical order
        s_es = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Substrate(:A),
                                   EnzymeRates.Substrate(:B)],
            :E,
        )
        @test EnzymeRates.name(s_es) === :E_A_B

        # Estar conformation
        s_estar = EnzymeRates.Species(EnzymeRates.Metabolite[], :Estar)
        @test EnzymeRates.name(s_estar) === :Estar

        # Three-arg constructor exposes residual
        res = EnzymeRates.Residual(
            [EnzymeRates.Substrate(:A)],
            [EnzymeRates.Product(:P)],
        )
        s_res = EnzymeRates.Species(EnzymeRates.Metabolite[], :Estar, res)
        @test EnzymeRates.has_residual(s_res)
        @test EnzymeRates.residual(s_res) == res
    end

    @testset "RegulatorySite: validation + accessors" begin
        lig_a = EnzymeRates.AllostericRegulator(:A)
        lig_b = EnzymeRates.AllostericRegulator(:B)
        site = EnzymeRates.RegulatorySite(
            [lig_a, lig_b], 4, [:OnlyR, :NonequalRT],
        )
        @test EnzymeRates.ligands(site) == [lig_a, lig_b]
        @test EnzymeRates.multiplicity(site) == 4
        @test EnzymeRates.allo_states(site) == [:OnlyR, :NonequalRT]

        # Mismatched ligand / allo_state length → error
        @test_throws ErrorException EnzymeRates.RegulatorySite(
            [lig_a, lig_b], 4, [:OnlyR])

        # Multiplicity < 1 → error
        @test_throws ErrorException EnzymeRates.RegulatorySite(
            [lig_a], 0, [:OnlyR])

        # Invalid allo state → error
        @test_throws ErrorException EnzymeRates.RegulatorySite(
            [lig_a], 1, [:NotAState])

        # All four allowed states accepted
        for st in (:OnlyR, :OnlyT, :EqualRT, :NonequalRT)
            @test EnzymeRates.RegulatorySite([lig_a], 1, [st]) isa
                  EnzymeRates.RegulatorySite
        end

        # Equality / hash
        site2 = EnzymeRates.RegulatorySite(
            [lig_a, lig_b], 4, [:OnlyR, :NonequalRT])
        @test site == site2
        @test hash(site) == hash(site2)
        # Order-sensitive: ligand ordering is parallel to allo_states,
        # so [A,B] / [OnlyR,NonequalRT] != [B,A] / [OnlyR,NonequalRT].
        site_reordered = EnzymeRates.RegulatorySite(
            [lig_b, lig_a], 4, [:OnlyR, :NonequalRT])
        @test site != site_reordered
    end

    @testset "Step has no source_idx field — rep_idx comes from position" begin
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Substrate(:S)], :E)
        e_p = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Product(:P)], :E)

        s1 = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
        s2 = EnzymeRates.Step(e_s, e_p, nothing, false)
        s3 = EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)

        @test fieldnames(EnzymeRates.Step) ==
              (:from_species, :to_species, :bound_metabolite, :is_equilibrium)
        @test !(:source_idx in fieldnames(EnzymeRates.Step))

        @test EnzymeRates.from_species(s1) === e
        @test EnzymeRates.to_species(s1) === e_s
        @test EnzymeRates.bound_metabolite(s1) ==
              EnzymeRates.Substrate(:S)
        @test EnzymeRates.is_equilibrium(s1)
        @test EnzymeRates.is_binding(s1)
        @test !EnzymeRates.is_iso(s1)
        @test EnzymeRates.direction(s1) === :binding

        @test EnzymeRates.bound_metabolite(s2) === nothing
        @test EnzymeRates.is_iso(s2)
        @test !EnzymeRates.is_binding(s2)
        @test EnzymeRates.direction(s2) === :iso

        @test EnzymeRates.is_binding(s3)
    end

    @testset "Step canonicalizes binding direction" begin
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Substrate(:S)], :E)

        # User authored release direction (E_S → E + S, metabolite on RHS).
        # Constructor swaps to binding direction (E + S → E_S).
        released = EnzymeRates.Step(e_s, e, EnzymeRates.Substrate(:S), true)
        bound    = EnzymeRates.Step(e,   e_s, EnzymeRates.Substrate(:S), true)
        @test released == bound
        @test hash(released) == hash(bound)
        @test EnzymeRates.from_species(released) === e
        @test EnzymeRates.to_species(released) === e_s
    end

    @testset "Step canonicalizes iso direction" begin
        e_s = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Substrate(:S)], :E)
        e_p = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Product(:P)], :E)
        fwd = EnzymeRates.Step(e_s, e_p, nothing, false)
        rev = EnzymeRates.Step(e_p, e_s, nothing, false)
        @test fwd == rev
        @test hash(fwd) == hash(rev)
    end

    @testset "Parameter family: step-bound, Kreg, mechanism-level" begin
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Substrate(:S)], :E)
        step = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)

        kd_none = EnzymeRates.Kd(step, :None)
        kd_t    = EnzymeRates.Kd(step, :T)
        @test kd_none isa EnzymeRates.Kd
        @test kd_none isa EnzymeRates.Parameter
        @test EnzymeRates.governing_step(kd_none) === step
        @test !EnzymeRates.is_t_state(kd_none)
        @test EnzymeRates.is_t_state(kd_t)
        @test kd_none == EnzymeRates.Kd(step, :None)
        @test kd_none != kd_t

        for T in (EnzymeRates.Kiso, EnzymeRates.Kon, EnzymeRates.Koff,
                  EnzymeRates.Kfor, EnzymeRates.Krev)
            p = T(step, :None)
            @test p isa EnzymeRates.Parameter
            @test EnzymeRates.governing_step(p) === step
            @test !EnzymeRates.is_t_state(p)
            @test EnzymeRates.is_t_state(T(step, :T))
        end

        lig_a = EnzymeRates.AllostericRegulator(:A)
        site = EnzymeRates.RegulatorySite([lig_a], 2, [:OnlyR])
        kr = EnzymeRates.Kreg(site, lig_a, :R)
        @test kr isa EnzymeRates.Parameter
        @test EnzymeRates.is_t_state(EnzymeRates.Kreg(site, lig_a, :T))
        @test !EnzymeRates.is_t_state(kr)
        @test kr == EnzymeRates.Kreg(site, lig_a, :R)

        # Mechanism-level scalars: singletons
        @test EnzymeRates.Keq() == EnzymeRates.Keq()
        @test EnzymeRates.Etot() == EnzymeRates.Etot()
        @test EnzymeRates.Lallo() == EnzymeRates.Lallo()
        @test EnzymeRates.Keq() isa EnzymeRates.Parameter
        @test EnzymeRates.Etot() isa EnzymeRates.Parameter
        @test EnzymeRates.Lallo() isa EnzymeRates.Parameter
    end

    @testset "ReactantAtoms canonicalizes atom ordering" begin
        ra1 = EnzymeRates.ReactantAtoms(
            EnzymeRates.Substrate(:ATP),
            [:C => 10, :H => 16, :N => 5],
        )
        ra2 = EnzymeRates.ReactantAtoms(
            EnzymeRates.Substrate(:ATP),
            [:N => 5, :H => 16, :C => 10],
        )
        @test ra1 == ra2
        @test hash(ra1) == hash(ra2)
        @test EnzymeRates.metabolite(ra1) == EnzymeRates.Substrate(:ATP)
        @test EnzymeRates.atoms(ra1) == [:C => 10, :H => 16, :N => 5]

        # Distinct metabolite kinds: ATP-Substrate != ATP-Product
        ra_p = EnzymeRates.ReactantAtoms(
            EnzymeRates.Product(:ATP), [:C => 10, :H => 16, :N => 5])
        @test ra1 != ra_p
    end

    @testset "RegulatorMults canonicalizes ordering + validates" begin
        rm1 = EnzymeRates.RegulatorMults(
            EnzymeRates.AllostericRegulator(:A), [4, 1, 2])
        rm2 = EnzymeRates.RegulatorMults(
            EnzymeRates.AllostericRegulator(:A), [1, 2, 4])
        @test rm1 == rm2
        @test hash(rm1) == hash(rm2)
        @test EnzymeRates.regulator(rm1) ==
              EnzymeRates.AllostericRegulator(:A)
        @test EnzymeRates.allowed_multiplicities(rm1) == [1, 2, 4]

        # Multiplicity < 1 → error
        @test_throws Exception EnzymeRates.RegulatorMults(
            EnzymeRates.AllostericRegulator(:A), [0, 1])
        @test_throws Exception EnzymeRates.RegulatorMults(
            EnzymeRates.AllostericRegulator(:A), [-1])

        # CompetitiveInhibitor also accepted
        rm_ci = EnzymeRates.RegulatorMults(
            EnzymeRates.CompetitiveInhibitor(:I), [1])
        @test EnzymeRates.regulator(rm_ci) ==
              EnzymeRates.CompetitiveInhibitor(:I)
    end

    @testset "EnzymeReaction (new concrete)" begin
        r = EnzymeReaction(
            [EnzymeRates.ReactantAtoms(
                 EnzymeRates.Substrate(:ATP), [:C => 10]),
             EnzymeRates.ReactantAtoms(
                 EnzymeRates.Product(:ADP), [:C => 10])],
            [EnzymeRates.RegulatorMults(
                 EnzymeRates.AllostericRegulator(:cAMP), [2])],
            [1, 2],
        )

        @test length(EnzymeRates.reactants(r)) == 2
        @test EnzymeRates.allowed_catalytic_multiplicities(r) == [1, 2]
        @test length(EnzymeRates.regulators(r)) == 1
        @test EnzymeRates.substrates(r) == [EnzymeRates.Substrate(:ATP)]
        @test EnzymeRates.products(r) == [EnzymeRates.Product(:ADP)]
    end

    @testset "EnzymeReaction canonicalizes reactant + regulator ordering" begin
        r1 = EnzymeReaction(
            [EnzymeRates.ReactantAtoms(
                 EnzymeRates.Substrate(:B), [:C => 1]),
             EnzymeRates.ReactantAtoms(
                 EnzymeRates.Substrate(:A), [:C => 1]),
             EnzymeRates.ReactantAtoms(
                 EnzymeRates.Product(:P), [:C => 2])],
            [EnzymeRates.RegulatorMults(
                 EnzymeRates.AllostericRegulator(:Y), [2]),
             EnzymeRates.RegulatorMults(
                 EnzymeRates.AllostericRegulator(:X), [2])],
            [3, 1, 2],
        )
        r2 = EnzymeReaction(
            [EnzymeRates.ReactantAtoms(
                 EnzymeRates.Substrate(:A), [:C => 1]),
             EnzymeRates.ReactantAtoms(
                 EnzymeRates.Substrate(:B), [:C => 1]),
             EnzymeRates.ReactantAtoms(
                 EnzymeRates.Product(:P), [:C => 2])],
            [EnzymeRates.RegulatorMults(
                 EnzymeRates.AllostericRegulator(:X), [2]),
             EnzymeRates.RegulatorMults(
                 EnzymeRates.AllostericRegulator(:Y), [2])],
            [1, 2, 3],
        )
        @test r1 == r2
        @test hash(r1) == hash(r2)
    end

    @testset "EnzymeReaction rejects multiplicity < 1" begin
        @test_throws ErrorException EnzymeReaction(
            EnzymeRates.ReactantAtoms[
                EnzymeRates.ReactantAtoms(
                    EnzymeRates.Substrate(:S), [:C => 1])],
            EnzymeRates.RegulatorMults[],
            [0, 1],
        )
    end

    @testset "Mechanism (non-parametric)" begin
        r = EnzymeReaction(
            [EnzymeRates.ReactantAtoms(
                 EnzymeRates.Substrate(:S), [:C => 1]),
             EnzymeRates.ReactantAtoms(
                 EnzymeRates.Product(:P), [:C => 1])],
            EnzymeRates.RegulatorMults[],
            [1],
        )
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
        e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)

        s_bind = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
        s_iso  = EnzymeRates.Step(e_s, e_p, nothing, false)
        s_rel  = EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)

        m = EnzymeRates.Mechanism(r, [[s_bind], [s_iso], [s_rel]])
        @test EnzymeRates.reaction(m) == r
        @test EnzymeRates.steps(m) == [[s_bind], [s_iso], [s_rel]]
        @test EnzymeRates.kinetic_groups(m) == 1:3
        @test EnzymeRates.n_steps(m) == 3
        @test EnzymeRates.rep_step(m, 2) == s_iso

        m2 = EnzymeRates.Mechanism(r, [[s_bind], [s_iso], [s_rel]])
        @test m == m2
        @test hash(m) == hash(m2)
    end

    @testset "AllostericMechanism (non-parametric)" begin
        r = EnzymeReaction(
            [EnzymeRates.ReactantAtoms(
                 EnzymeRates.Substrate(:S), [:C => 1]),
             EnzymeRates.ReactantAtoms(
                 EnzymeRates.Product(:P), [:C => 1])],
            EnzymeRates.RegulatorMults[],
            [2],
        )
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
        e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)
        s_bind = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
        s_iso  = EnzymeRates.Step(e_s, e_p, nothing, false)
        s_rel  = EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)

        site = EnzymeRates.RegulatorySite(
            [EnzymeRates.AllostericRegulator(:I)], 1, [:OnlyT])
        m = EnzymeRates.AllostericMechanism(
            r, [[s_bind], [s_iso], [s_rel]],
            [:EqualRT, :NonequalRT, :OnlyR], 2,
            [site])

        @test EnzymeRates.reaction(m) == r
        @test EnzymeRates.steps(m) == [[s_bind], [s_iso], [s_rel]]
        @test EnzymeRates.cat_allo_state(m, 1) == :EqualRT
        @test EnzymeRates.cat_allo_state(m, 3) == :OnlyR
        @test EnzymeRates.catalytic_multiplicity(m) == 2
        @test EnzymeRates.regulatory_sites(m) == [site]
        @test EnzymeRates.kinetic_groups(m) == 1:3
        @test EnzymeRates.n_steps(m) == 3
        @test EnzymeRates.rep_step(m, 2) == s_iso
        @test EnzymeRates.allosteric_regulators(m) ==
              [EnzymeRates.AllostericRegulator(:I)]

        m2 = EnzymeRates.AllostericMechanism(
            r, [[s_bind], [s_iso], [s_rel]],
            [:EqualRT, :NonequalRT, :OnlyR], 2,
            [site])
        @test m == m2
        @test hash(m) == hash(m2)
    end

    @testset "AllostericMechanism validation errors" begin
        r = EnzymeReaction(
            [EnzymeRates.ReactantAtoms(
                 EnzymeRates.Substrate(:S), [:C => 1]),
             EnzymeRates.ReactantAtoms(
                 EnzymeRates.Product(:P), [:C => 1])],
            EnzymeRates.RegulatorMults[],
            [1],
        )
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
        s_bind = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
        cat_steps = [[s_bind]]

        # :OnlyT for catalytic group is rejected (R-state-active convention)
        @test_throws ErrorException EnzymeRates.AllostericMechanism(
            r, cat_steps, [:OnlyT], 1,
            EnzymeRates.RegulatorySite[])

        # Length mismatch
        @test_throws ErrorException EnzymeRates.AllostericMechanism(
            r, cat_steps, [:EqualRT, :NonequalRT], 1,
            EnzymeRates.RegulatorySite[])

        # catalytic_multiplicity < 1
        @test_throws ErrorException EnzymeRates.AllostericMechanism(
            r, cat_steps, [:EqualRT], 0,
            EnzymeRates.RegulatorySite[])

        # Unknown allo-state symbol
        @test_throws ErrorException EnzymeRates.AllostericMechanism(
            r, cat_steps, [:Bogus], 1,
            EnzymeRates.RegulatorySite[])
    end
end
