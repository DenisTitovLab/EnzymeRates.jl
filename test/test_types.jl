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
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
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
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
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
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
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

            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S ⇌ E(S)     :: EqualRT
                E(S) <--> E(P)   :: OnlyR
                E(P) ⇌ E + P     :: EqualRT
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
                E + S <--> E(S)
                E(S) <--> E + P
            end
        end
        @test sprint(show, m) ==
            "EnzymeMechanism: E + S <--> E_S <--> E + P"

        # Branched mechanism: multi-line with header summary.
        m_b = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                E + A <--> E(A)
                E + B <--> E(B)
                E(A) + B <--> E(A, B)
                E(B) + A <--> E(A, B)
                E(A, B) <--> E(P, Q)
                E(P, Q) <--> E(Q) + P
                E(Q) <--> E + Q
            end
        end
        s = sprint(show, m_b)
        @test startswith(s, "EnzymeMechanism (7 steps, 6 enzyme forms):")
        @test contains(s, "E + A <--> E_A")
        @test contains(s, "E_Q <--> E + Q")

        # Linear chain with canonical RE binding: chain-walk renders the
        # release step in reverse so the chain stays linear.
        m_re = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        @test sprint(show, m_re) ==
            "EnzymeMechanism: E + S ⇌ E_S <--> E_P ⇌ E + P"

        # Mechanism with regulators: appended at end.
        m_reg = @enzyme_mechanism begin
            substrates: S
            products:   P
            regulators: I
            steps: begin
                E + S <--> E(S)
                E(S) <--> E + P
                E + I <--> E(I)
            end
        end
        @test contains(sprint(show, m_reg), "| regulators: I")

        # EnzymeReaction with oligomeric_state > 1.
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
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + F6P ⇌ E(F6P)         :: EqualRT
                E(F6P) <--> E(F16BP)     :: EqualRT
                E(F16BP) ⇌ E + F16BP     :: EqualRT
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
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
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
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S ⇌ E(S)     :: NonequalRT
                E(S) <--> E(P)   :: OnlyR
                E(P) ⇌ E + P     :: EqualRT
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
        base = first(EnzymeRates.init_mechanisms(rxn))
        cat_allo_states = Symbol[]
        for g in EnzymeRates.kinetic_groups(base)
            rep = EnzymeRates.rep_step(base, g)
            met = EnzymeRates.bound_metabolite(rep)
            tag = (met isa EnzymeRates.Reactant &&
                   EnzymeRates.name(met) in (:S, :P)) ?
                  :EqualRT : :NonequalRT
            push!(cat_allo_states, tag)
        end
        site = EnzymeRates.RegulatorySite(
            [EnzymeRates.AllostericRegulator(:R)], 2, [:OnlyT])
        am = EnzymeRates.AllostericMechanism(
            EnzymeRates.reaction(base), copy(EnzymeRates.steps(base)),
            cat_allo_states, 2, [site])
        m = EnzymeRates.AllostericEnzymeMechanism(am)
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

    @testset "Step carries source_idx presentation metadata" begin
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Substrate(:S)], :E)
        e_p = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Product(:P)], :E)

        s1 = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
        s2 = EnzymeRates.Step(e_s, e_p, nothing, false)
        s3 = EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)

        @test fieldnames(EnzymeRates.Step) ==
              (:from_species, :to_species, :bound_metabolite,
               :is_equilibrium, :source_idx)
        @test :source_idx in fieldnames(EnzymeRates.Step)
        # Default is 0 (unset) when not supplied; the `Mechanism`
        # constructor auto-fills by flat position.
        @test s1.source_idx == 0
        @test s2.source_idx == 0
        @test s3.source_idx == 0

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

    @testset "Step `==` / hash ignore source_idx" begin
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Substrate(:S)], :E)
        s   = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
        s_at_5 = EnzymeRates.Step(
            e, e_s, EnzymeRates.Substrate(:S), true; source_idx = 5)
        @test s == s_at_5
        @test hash(s) == hash(s_at_5)
        @test s.source_idx == 0
        @test s_at_5.source_idx == 5
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

    @testset "Step canonicalizes RE iso direction, preserves SS iso" begin
        e_s = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Substrate(:S)], :E)
        e_p = EnzymeRates.Species(
            EnzymeRates.Metabolite[EnzymeRates.Product(:P)], :E)

        # RE iso: deterministic lex-on-name canonicalization for dedup.
        re_fwd = EnzymeRates.Step(e_s, e_p, nothing, true)
        re_rev = EnzymeRates.Step(e_p, e_s, nothing, true)
        @test re_fwd == re_rev
        @test hash(re_fwd) == hash(re_rev)

        # SS iso: direction preserved (kf/kr labels are direction-sensitive;
        # analytical formulas reference :kNf as source-forward). See
        # CLAUDE.md "Canonical Step Form".
        ss_fwd = EnzymeRates.Step(e_s, e_p, nothing, false)
        ss_rev = EnzymeRates.Step(e_p, e_s, nothing, false)
        @test ss_fwd != ss_rev
        @test EnzymeRates.from_species(ss_fwd) === e_s
        @test EnzymeRates.from_species(ss_rev) === e_p
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

    @testset "Mechanism auto-assigns source_idx by flat position" begin
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

        s1 = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
        s2 = EnzymeRates.Step(e_s, e_p, nothing, false)
        s3 = EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)
        @test s1.source_idx == 0
        @test s2.source_idx == 0
        @test s3.source_idx == 0

        m = EnzymeRates.Mechanism(r, [[s1], [s2], [s3]])
        @test m.steps[1][1].source_idx == 1
        @test m.steps[2][1].source_idx == 2
        @test m.steps[3][1].source_idx == 3
    end

    @testset "Mechanism preserves explicit caller source_idx" begin
        # Legacy lift sets source_idx explicitly based on the original
        # source position (before grouping). The Mechanism ctor must
        # preserve these so `_rep_idx_for_step` reproduces today's
        # source-order rep-idx naming for mechanisms whose group
        # members are non-contiguous in source order.
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

        s1 = EnzymeRates.Step(
            e, e_s, EnzymeRates.Substrate(:S), true; source_idx = 1)
        s2 = EnzymeRates.Step(
            e_s, e_p, nothing, false; source_idx = 7)
        s3 = EnzymeRates.Step(
            e, e_p, EnzymeRates.Product(:P), true; source_idx = 3)
        m = EnzymeRates.Mechanism(r, [[s1], [s2], [s3]])
        @test m.steps[1][1].source_idx == 1
        @test m.steps[2][1].source_idx == 7
        @test m.steps[3][1].source_idx == 3
    end

    @testset "Mechanism rejects mixed set / unset source_idx" begin
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

        s_set   = EnzymeRates.Step(
            e, e_s, EnzymeRates.Substrate(:S), true; source_idx = 1)
        s_unset = EnzymeRates.Step(e_s, e_p, nothing, false)
        @test_throws ErrorException EnzymeRates.Mechanism(
            r, [[s_set], [s_unset]])
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

    @testset "AllostericMechanism(::AllostericEnzymeMechanism) converter" begin
        cm = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        aem = EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:EqualRT, :NonequalRT, :OnlyR)),
            (((:A, :B), 1, (:OnlyR, :NonequalRT)),),
        )

        am = EnzymeRates.AllostericMechanism(aem)
        @test am isa EnzymeRates.AllostericMechanism

        # Catalytic side lifted via Mechanism(CM())
        @test EnzymeRates.steps(am) == EnzymeRates.Mechanism(cm).steps
        @test EnzymeRates.reaction(am) == EnzymeRates.Mechanism(cm).reaction

        # cat_allo_states and multiplicity extracted from CS
        @test EnzymeRates.cat_allo_state(am, 1) === :EqualRT
        @test EnzymeRates.cat_allo_state(am, 2) === :NonequalRT
        @test EnzymeRates.cat_allo_state(am, 3) === :OnlyR
        @test EnzymeRates.catalytic_multiplicity(am) == 2

        # Regulatory sites: ligand Symbols wrapped as AllostericRegulator
        sites = EnzymeRates.regulatory_sites(am)
        @test length(sites) == 1
        @test sites[1].ligands ==
              [EnzymeRates.AllostericRegulator(:A),
               EnzymeRates.AllostericRegulator(:B)]
        @test sites[1].multiplicity == 1
        @test sites[1].allo_states == [:OnlyR, :NonequalRT]

        # Idempotent: two calls give equal results
        @test EnzymeRates.AllostericMechanism(aem) == am

        # _to_mechanism bridge dispatches to AllostericMechanism
        @test EnzymeRates._to_mechanism(aem) == am
        @test EnzymeRates._to_mechanism(aem) isa EnzymeRates.AllostericMechanism
        # And matches the non-allosteric side too
        @test EnzymeRates._to_mechanism(cm) == EnzymeRates.Mechanism(cm)
    end

    @testset "name(p, ::AllostericEnzymeMechanism) chokepoint overloads" begin
        cm = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        aem = EnzymeRates.AllostericEnzymeMechanism(
            cm, (2, (:NonequalRT, :EqualRT, :NonequalRT)),
            (((:R,), 1, (:NonequalRT,)),),
        )
        am = EnzymeRates.AllostericMechanism(aem)

        # Step-bound parameters: AEM dispatch matches AM dispatch
        rep_bind = first(EnzymeRates.steps(am)[1])
        @test EnzymeRates.name(EnzymeRates.Kd(rep_bind, :None), aem) ==
              EnzymeRates.name(EnzymeRates.Kd(rep_bind, :None), am)
        @test EnzymeRates.name(EnzymeRates.Kd(rep_bind, :T), aem) ==
              EnzymeRates.name(EnzymeRates.Kd(rep_bind, :T), am)
        @test EnzymeRates.name(EnzymeRates.Kd(rep_bind, :T), aem) === :K1_T

        rep_iso  = first(EnzymeRates.steps(am)[2])
        @test EnzymeRates.name(EnzymeRates.Kfor(rep_iso, :None), aem) ==
              EnzymeRates.name(EnzymeRates.Kfor(rep_iso, :None), am)
        @test EnzymeRates.name(EnzymeRates.Kfor(rep_iso, :None), aem) === :k2f

        # Kreg: AEM dispatch matches AM dispatch
        site = EnzymeRates.regulatory_sites(am)[1]
        lig  = first(site.ligands)
        @test EnzymeRates.name(EnzymeRates.Kreg(site, lig, :R), aem) ==
              EnzymeRates.name(EnzymeRates.Kreg(site, lig, :R), am)
        @test EnzymeRates.name(EnzymeRates.Kreg(site, lig, :T), aem) ===
              :K_R_T_reg1
    end

    @testset "_sig_of / _mechanism_from_sig roundtrip" begin
        # Use real atom data to catch type-parameter validity issues.
        # Pair{Symbol,Int} is NOT a valid type-parameter value; encoding
        # must use Tuple{Symbol,Int} leaves.
        r = EnzymeRates.EnzymeReaction(
            [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:ATP),
                                        [:C => 10, :H => 16, :N => 5,
                                         :O => 13, :P => 3]),
             EnzymeRates.ReactantAtoms(EnzymeRates.Product(:ADP),
                                        [:C => 10, :H => 15, :N => 5,
                                         :O => 10, :P => 2])],
            EnzymeRates.RegulatorMults[],
            [1],
        )
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species([EnzymeRates.Substrate(:ATP)], :E)
        e_p = EnzymeRates.Species([EnzymeRates.Product(:ADP)], :E)

        m = EnzymeRates.Mechanism(r, [
            [EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:ATP), true)],
            [EnzymeRates.Step(e_s, e_p, nothing, false)],
            [EnzymeRates.Step(e, e_p, EnzymeRates.Product(:ADP), true)],
        ])

        sig = EnzymeRates._sig_of(m)
        @test sig isa Tuple

        m_recon = EnzymeRates._mechanism_from_sig(sig)
        @test m_recon == m   # roundtrip

        # CRITICAL: sig MUST be usable as a type parameter. Throws TypeError
        # if any leaf is invalid (Pair, Vector, DataType inside value-tuple).
        em_type = EnzymeRates.EnzymeMechanism{sig}
        @test em_type <: EnzymeRates.EnzymeMechanism
        em_inst = em_type()
        @test em_inst isa EnzymeRates.EnzymeMechanism

        # Roundtrip through the type-parameter form preserves everything.
        @test EnzymeRates.Mechanism(em_inst) == m
    end

    @testset "Mechanism <-> EnzymeMechanism converters" begin
        r = EnzymeRates.EnzymeReaction(
            [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1]),
             EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])],
            EnzymeRates.RegulatorMults[],
            [1],
        )
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
        e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)

        m = EnzymeRates.Mechanism(r, [
            [EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)],
            [EnzymeRates.Step(e_s, e_p, nothing, false)],
            [EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)],
        ])

        em = EnzymeMechanism(m)
        @test em isa EnzymeMechanism

        m_back = EnzymeRates.Mechanism(em)
        @test m_back == m
    end

    @testset "name(p::Parameter, m) chokepoint" begin
        r = EnzymeRates.EnzymeReaction(
            [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1]),
             EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])],
            EnzymeRates.RegulatorMults[],
            [1],
        )
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
        e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)

        step1 = EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)
        step2 = EnzymeRates.Step(e_s, e_p, nothing, false)
        step3 = EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)

        m = EnzymeRates.Mechanism(r, [[step1], [step2], [step3]])

        # Positional naming: rep_idx for kinetic group g = position of first
        # step in steps(m)[g] within the flattened steps list.
        @test EnzymeRates.name(EnzymeRates.Kd(step1, :None), m) === :K1
        @test EnzymeRates.name(EnzymeRates.Kd(step1, :T),    m) === :K1_T
        @test EnzymeRates.name(EnzymeRates.Kon(step2, :None), m) === :k2f
        @test EnzymeRates.name(EnzymeRates.Koff(step2, :None), m) === :k2r
        @test EnzymeRates.name(EnzymeRates.Kfor(step2, :None), m) === :k2f
        @test EnzymeRates.name(EnzymeRates.Krev(step2, :None), m) === :k2r
        @test EnzymeRates.name(EnzymeRates.Kd(step3, :None), m) === :K3

        # T-suffix on SS step
        @test EnzymeRates.name(EnzymeRates.Kon(step2, :T),  m) === :k2f_T
        @test EnzymeRates.name(EnzymeRates.Koff(step2, :T), m) === :k2r_T

        # Kiso uses K-naming (RE iso)
        @test EnzymeRates.name(EnzymeRates.Kiso(step2, :None), m) === :K2
        @test EnzymeRates.name(EnzymeRates.Kiso(step2, :T),    m) === :K2_T

        # Mechanism-level scalars
        @test EnzymeRates.name(EnzymeRates.Keq(),   m) === :Keq
        @test EnzymeRates.name(EnzymeRates.Etot(),  m) === :E_total
        @test EnzymeRates.name(EnzymeRates.Lallo(), m) === :L

        # Same names resolve via EnzymeMechanism(m) (the parametric form).
        em = EnzymeMechanism(m)
        @test EnzymeRates.name(EnzymeRates.Kd(step1, :None), em) === :K1
        @test EnzymeRates.name(EnzymeRates.Kon(step2, :None), em) === :k2f
        @test EnzymeRates.name(EnzymeRates.Keq(),   em) === :Keq
        @test EnzymeRates.name(EnzymeRates.Etot(),  em) === :E_total
        @test EnzymeRates.name(EnzymeRates.Lallo(), em) === :L
    end

    @testset "name(p::Parameter, m) rep_idx for shared kinetic group" begin
        # Group with 2 steps: rep is the first step's position in the
        # flattened step list. If group 1 contains steps at positions 1
        # and 2, rep_idx is 1; if group 2 starts at position 3 with two
        # steps, rep_idx for group 2 is 3.
        r = EnzymeRates.EnzymeReaction(
            [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1]),
             EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])],
            EnzymeRates.RegulatorMults[],
            [1],
        )
        e    = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s  = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
        e_p  = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)
        e_sp = EnzymeRates.Species(
            EnzymeRates.Metabolite[
                EnzymeRates.Substrate(:S), EnzymeRates.Product(:P)], :E)

        step_a = EnzymeRates.Step(e,   e_s,  EnzymeRates.Substrate(:S), true)
        step_b = EnzymeRates.Step(e_p, e_sp, EnzymeRates.Substrate(:S), true)
        step_c = EnzymeRates.Step(e_s, e_p,  nothing, false)
        step_d = EnzymeRates.Step(e,   e_p,  EnzymeRates.Product(:P), true)

        m = EnzymeRates.Mechanism(r, [[step_a, step_b], [step_c], [step_d]])

        # Group 1: steps 1, 2 — rep_idx = 1 for either step.
        @test EnzymeRates.name(EnzymeRates.Kd(step_a, :None), m) === :K1
        @test EnzymeRates.name(EnzymeRates.Kd(step_b, :None), m) === :K1
        # Group 2: step 3 — rep_idx = 3.
        @test EnzymeRates.name(EnzymeRates.Kon(step_c, :None), m) === :k3f
        # Group 3: step 4 — rep_idx = 4.
        @test EnzymeRates.name(EnzymeRates.Kd(step_d, :None), m) === :K4
    end

    @testset "name(p::Kreg, m) chokepoint" begin
        r = EnzymeRates.EnzymeReaction(
            [EnzymeRates.ReactantAtoms(EnzymeRates.Substrate(:S), [:C => 1]),
             EnzymeRates.ReactantAtoms(EnzymeRates.Product(:P), [:C => 1])],
            [EnzymeRates.RegulatorMults(
                 EnzymeRates.AllostericRegulator(:A), [2])],
            [2],
        )
        e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
        e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
        e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)

        cat_steps = [
            [EnzymeRates.Step(e, e_s, EnzymeRates.Substrate(:S), true)],
            [EnzymeRates.Step(e_s, e_p, nothing, false)],
            [EnzymeRates.Step(e, e_p, EnzymeRates.Product(:P), true)],
        ]
        site_a = EnzymeRates.RegulatorySite(
            [EnzymeRates.AllostericRegulator(:A)], 2, [:NonequalRT])
        am = EnzymeRates.AllostericMechanism(
            r, cat_steps,
            [:EqualRT, :EqualRT, :EqualRT], 2, [site_a])

        @test EnzymeRates.name(
            EnzymeRates.Kreg(site_a, EnzymeRates.AllostericRegulator(:A), :R),
            am) === :K_A_reg1
        @test EnzymeRates.name(
            EnzymeRates.Kreg(site_a, EnzymeRates.AllostericRegulator(:A), :T),
            am) === :K_A_T_reg1

        # Step-bound parameters also resolve via AllostericMechanism.
        rep = first(cat_steps[1])
        @test EnzymeRates.name(EnzymeRates.Kd(rep, :None), am) === :K1
        @test EnzymeRates.name(EnzymeRates.Kd(rep, :T),    am) === :K1_T

        # Iso step in second kinetic group
        iso_step = first(cat_steps[2])
        @test EnzymeRates.name(EnzymeRates.Kiso(iso_step, :None), am) === :K2
        @test EnzymeRates.name(EnzymeRates.Kon(iso_step, :None),  am) === :k2f

        # Scalars also dispatch on AllostericMechanism
        @test EnzymeRates.name(EnzymeRates.Keq(),   am) === :Keq
        @test EnzymeRates.name(EnzymeRates.Etot(),  am) === :E_total
        @test EnzymeRates.name(EnzymeRates.Lallo(), am) === :L
    end
end
