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
                [E, S] ⇌ [ES]
                [ES] <--> [EP]
                [EP] ⇌ [E, P]
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
                    [E, S] ⇌ [ES]    :: EqualRT
                    [ES] <--> [EP]    :: OnlyR
                    [EP] ⇌ [E, P]    :: EqualRT
                end
            end
        end
        @test EnzymeRates.catalytic_multiplicity(m) == 2
        @test EnzymeRates.cat_allo_state(m, 1) == :EqualRT
        @test EnzymeRates.cat_allo_state(m, 2) == :OnlyR
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

        # Linear mechanism: compact chain.
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                [E, S] <--> [ES]
                [ES] <--> [E, P]
            end
        end
        @test sprint(show, m) ==
            "EnzymeMechanism: E + S <--> ES <--> E + P"

        # Branched mechanism: multi-line with header summary.
        m_b = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                [E, A] <--> [EA]
                [E, B] <--> [EB]
                [EA, B] <--> [EAB]
                [EB, A] <--> [EAB]
                [EAB] <--> [EPQ]
                [EPQ] <--> [EQ, P]
                [EQ] <--> [E, Q]
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
                [E, S] ⇌ [ES]
                [ES] <--> [EP]
                [EP] ⇌ [E, P]
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
                [E, S] <--> [ES]
                [ES] <--> [E, P]
                [E, I] <--> [EI]
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
            site(:catalytic, 2): begin
                steps: begin
                    [E, F6P] ⇌ [E_F6P]   :: EqualRT
                    [E_F6P] <--> [E_F16BP] :: EqualRT
                    [E_F16BP] ⇌ [E, F16BP] :: EqualRT
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
                [E, S] ⇌ [ES]
                [ES] <--> [EP]
                [EP] ⇌ [E, P]
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
    end

    @testset "Base.show displays all dense states" begin
        m = @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: I::NonequalRT, J::OnlyT
            site(:catalytic, 2): begin
                steps: begin
                    [E, S] ⇌ [ES]    :: NonequalRT
                    [ES] <--> [EP]   :: OnlyR
                    [EP] ⇌ [E, P]   :: EqualRT
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
end
