_used_forms = EnzymeRates._used_forms

@testset "Mechanism Enumeration" begin

    # ─── Test Reactions ──────────────────────────────────────

    uni_uni = @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
    end

    bi_bi = @enzyme_reaction begin
        substrates: A[C], B[N]
        products:   P[C], Q[N]
    end

    pingpong = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products:   P[C], Q[NX]
    end

    uni_uni_act = @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        regulators: A[N]
    end

    uni_uni_inh = @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        regulators: I[N]
    end

    uni_uni_2reg = @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        regulators: I[N], J[P2]
    end

    # ─── Spec-Driven Pipeline Tests ──────────────────────────

    @testset "Pipeline: $(spec.name)" for spec in ENUMERATION_TEST_SPECS
        rxn = spec.reaction
        mf = spec.max_forms

        # Stage 1: enumerate_enzyme_forms
        forms = enumerate_enzyme_forms(rxn)
        @test length(forms) == spec.expected_n_forms

        # Stage 2: catalytic topologies
        cat_specs = EnzymeRates._enumerate_only_catalytic_mechanisms(
            forms, rxn; max_forms=mf,
        )
        @test length(cat_specs) == spec.expected_n_catalytic

        # Stage 3: activator configs
        act_specs = EnzymeRates.MechanismSpec[
            s for c in cat_specs
            for s in EnzymeRates._generate_activator_configs(
                c, forms, rxn,
            )
        ]
        filter!(
            s -> EnzymeRates._used_form_count(s) <= mf, act_specs,
        )
        @test length(act_specs) == spec.expected_n_cat_with_act

        # Stage 4: dead-end configs
        de_specs = EnzymeRates.MechanismSpec[
            s for a in act_specs
            for s in EnzymeRates._enumerate_dead_end_configs(
                a, forms; max_forms=mf,
            )
        ]
        @test length(de_specs) == spec.expected_n_cat_act_de

        # Independent dead-end verification: for each base spec,
        # check that our independent counter agrees
        for base in act_specs
            independent_de = _compute_expected_dead_end_count(
                base, forms, mf,
            )
            actual_de = length(
                EnzymeRates._enumerate_dead_end_configs(
                    base, forms; max_forms=mf,
                ),
            )
            @test independent_de == actual_de
        end

        # Stage 5+6: RE/SS + constraints
        if !spec.skip_ress_test
            # Independent RE/SS count verification
            for base in de_specs
                independent_ress = _compute_independent_ress_count(
                    base, forms,
                )
                actual_ress = length(
                    EnzymeRates._enumerate_ress_and_constraints(
                        base, forms,
                    ),
                )
                @test independent_ress == actual_ress
            end

            # Full pipeline total
            t = @elapsed begin
                all_mechs = collect(
                    enumerate_mechanisms(rxn; max_forms=mf),
                )
            end
            @test length(all_mechs) == spec.expected_n_total

            # Performance check
            if isfinite(spec.max_enumeration_time)
                @test t < spec.max_enumeration_time
            end
        end

        # All base specs compile to EnzymeMechanism
        @test all(cat_specs) do s
            EnzymeMechanism(s) isa EnzymeMechanism
        end
    end

    # ─── Section 1: enumerate_enzyme_forms ─────────────────

    @testset "enumerate_enzyme_forms" begin
        @testset "Uni-Uni" begin
            forms = enumerate_enzyme_forms(uni_uni)
            names = Set(f.name for f in forms)
            @test names == Set([:E_0_0, :E_S_0, :E_0_P])
            # No form has all subs full AND any prod occupied
            @test all(forms) do f
                all_sub = all(
                    s -> s.role != :sub || s.index != 1 ||
                        s.atoms == s.full_atoms,
                    f.sites,
                )
                any_prod = any(
                    s -> s.role == :prod && s.index == 1 &&
                        s.atoms !== nothing,
                    f.sites,
                )
                !(all_sub && any_prod)
            end
            @test any(
                f -> all(s -> s.atoms === nothing, f.sites),
                forms,
            )
        end

        @testset "Uni-Uni max_sites=2" begin
            r = @enzyme_reaction begin
                substrates: S[C, 2]
                products:   P[C]
            end
            forms = enumerate_enzyme_forms(r)
            names = Set(f.name for f in forms)
            # S1, P1, S2 -> 2^3 - 3 excluded = 5
            @test length(forms) == 5
            @test :E_0_0_0 ∈ names
            @test :E_S_0_0 ∈ names
            @test :E_0_0_S ∈ names
            @test :E_S_0_S ∈ names
            @test :E_0_P_0 ∈ names
        end

        @testset "Bi-Bi" begin
            forms = enumerate_enzyme_forms(bi_bi)
            names = Set(f.name for f in forms)
            @test :E_0_0_0_0 ∈ names
            @test :E_A_0_0_0 ∈ names
            @test :E_0_0_P_Q ∈ names
            # Site ordering: core subs, core prods, extras, regs
            f = first(f for f in forms if f.name == :E_0_0_0_0)
            @test f.sites[1].metabolite == :A
            @test f.sites[1].role == :sub
            @test f.sites[2].metabolite == :B
            @test f.sites[2].role == :sub
            @test f.sites[3].metabolite == :P
            @test f.sites[3].role == :prod
            @test f.sites[4].metabolite == :Q
            @test f.sites[4].role == :prod
        end

        @testset "Ping-Pong" begin
            forms = enumerate_enzyme_forms(pingpong)
            names = Set(f.name for f in forms)
            @test :E_X_0_0_0 ∈ names
            @test :E_X_B_0_0 ∈ names
            # Residual sites: atoms != nothing AND atoms != full
            rf = first(f for f in forms if f.name == :E_X_0_0_0)
            rs = first(
                s for s in rf.sites
                if s.role == :sub && s.index == 1
            )
            @test rs.atoms !== nothing
            @test rs.atoms != rs.full_atoms
        end

        @testset "Regulators" begin
            forms = enumerate_enzyme_forms(uni_uni_inh)
            names = Set(f.name for f in forms)
            @test :E_0_0_0 ∈ names
            @test :E_S_0_I ∈ names
            @test :E_0_P_I ∈ names
            f = first(f for f in forms if f.name == :E_0_0_0)
            @test f.sites[end].role == :reg
        end

        @testset "No-atom species" begin
            r = @enzyme_reaction begin
                substrates: S
                products:   P
            end
            forms = enumerate_enzyme_forms(r)
            @test length(forms) == 3
            @test Set(f.name for f in forms) ==
                Set([:E_0_0, :E_S_0, :E_0_P])
        end

        @testset "SiteState enrichment" begin
            forms = enumerate_enzyme_forms(uni_uni_inh)
            f = first(f for f in forms if f.name == :E_0_0_0)
            @test f.sites[1].role == :sub
            @test f.sites[1].full_atoms == [:C => 1]
            @test f.sites[2].role == :prod
            @test f.sites[2].full_atoms == [:C => 1]
            @test f.sites[3].role == :reg
            @test f.sites[3].full_atoms == [:N => 1]
        end

        @testset "Backward compat: 2-element tuples" begin
            r = EnzymeReaction(
                ((:S, ((:C, 1),)),),
                ((:P, ((:C, 1),)),),
            )
            @test length(enumerate_enzyme_forms(r)) == 3
        end

        @testset "DSL produces 3-element tuples" begin
            subs = EnzymeRates.substrates(uni_uni)
            @test length(subs[1]) == 3
            @test subs[1][3] == 1
            r2 = @enzyme_reaction begin
                substrates: S[C, 3]
                products:   P[C]
            end
            @test EnzymeRates.substrates(r2)[1][3] == 3
        end

        @testset "Site ordering (extra sites)" begin
            r = @enzyme_reaction begin
                substrates: A[C, 2], B[N]
                products:   P[C], Q[N]
            end
            forms = enumerate_enzyme_forms(r)
            f = first(
                f for f in forms if f.name == :E_0_0_0_0_0
            )
            @test length(f.sites) == 5
            @test f.sites[1].metabolite == :A
            @test f.sites[1].index == 1
            @test f.sites[2].metabolite == :B
            @test f.sites[2].index == 1
            @test f.sites[3].metabolite == :P
            @test f.sites[3].index == 1
            @test f.sites[4].metabolite == :Q
            @test f.sites[4].index == 1
            @test f.sites[5].metabolite == :A
            @test f.sites[5].index == 2
        end

        @testset "Multi-product ping-pong residual" begin
            r = @enzyme_reaction begin
                substrates: A[C2X], B[N2]
                products:   P1[C], P2[C], Q[N2X]
            end
            forms = enumerate_enzyme_forms(r)
            names = Set(f.name for f in forms)
            @test :E_CX_0_0_0_0 ∈ names
            @test :E_X_0_0_0_0 ∈ names
        end

        @testset "Compilation time: 7-site reaction" begin
            script = """
            using EnzymeRates
            rxn = @enzyme_reaction begin
                substrates: Glu[C6H12O6], ATP[C10H16N5O13P3]
                products: G6P[C6H13O9P], ADP[C10H15N5O10P2]
                regulators: Phosphate[PO4], G6P[C6H13O9P], G6P[C6H13O9P]
            end
            t1 = @elapsed enumerate_enzyme_forms(rxn)
            t2 = @elapsed enumerate_enzyme_forms(rxn)
            print(t1, " ", t2)
            """
            proj = dirname(@__DIR__)
            out = read(
                `$(Base.julia_cmd()) --project=$proj -e $script`,
                String,
            )
            t_cold, t_warm = parse.(Float64, split(out))
            @test t_cold < 1.0
            @test t_warm < 0.1
        end
    end

    # ─── Section 2: _enumerate_only_catalytic_mechanisms ───────

    @testset "_enumerate_only_catalytic_mechanisms" begin
        @testset "Uni-Uni: single topology" begin
            forms = enumerate_enzyme_forms(uni_uni)
            specs = EnzymeRates._enumerate_only_catalytic_mechanisms(
                forms, uni_uni; max_forms=6,
            )
            spec = specs[1]
            # 3 steps: bind S, isomerize, release P
            @test length(spec.reactions) == 3
            @test EnzymeRates._used_form_count(spec) == 3
            @test all(==(false), spec.equilibrium_steps)
            @test isempty(spec.param_constraints)
            @test Set(spec.forms) == Set(f.name for f in forms)
            @test _used_forms(spec) ==
                Set([:E_0_0, :E_S_0, :E_0_P])
        end

        @testset "Bi-Bi (max_forms=5): ordered cycles" begin
            forms = enumerate_enzyme_forms(bi_bi)
            specs = EnzymeRates._enumerate_only_catalytic_mechanisms(
                forms, bi_bi; max_forms=5,
            )
            @test all(specs) do spec
                EnzymeRates._used_form_count(spec) == 5 &&
                    all(==(false), spec.equilibrium_steps) &&
                    isempty(spec.param_constraints) &&
                    :E_0_0_0_0 ∈ _used_forms(spec)
            end
        end

        @testset "Bi-Bi (max_forms=7): includes multi-cycle" begin
            forms = enumerate_enzyme_forms(bi_bi)
            specs5 = EnzymeRates._enumerate_only_catalytic_mechanisms(
                forms, bi_bi; max_forms=5,
            )
            specs7 = EnzymeRates._enumerate_only_catalytic_mechanisms(
                forms, bi_bi; max_forms=7,
            )
            @test length(specs7) > length(specs5)
            @test any(
                s -> EnzymeRates._used_form_count(s) > 5,
                specs7,
            )
        end

        @testset "Ping-Pong: topologies exist" begin
            forms = enumerate_enzyme_forms(pingpong)
            specs = EnzymeRates._enumerate_only_catalytic_mechanisms(
                forms, pingpong; max_forms=20,
            )
            @test all(
                s -> :E_0_0_0_0 ∈ _used_forms(s), specs,
            )
            @test all(
                s -> all(==(false), s.equilibrium_steps),
                specs,
            )
        end
    end

    # ─── Section 3: _generate_activator_configs ────────────

    @testset "_generate_activator_configs" begin
        @testset "No regulators: returns input only" begin
            forms = enumerate_enzyme_forms(uni_uni)
            specs = EnzymeRates._enumerate_only_catalytic_mechanisms(
                forms, uni_uni; max_forms=6,
            )
            input = specs[1]
            result = EnzymeRates._generate_activator_configs(
                input, forms, uni_uni,
            )
            @test result[1].reactions == input.reactions
        end

        @testset "Activator: absent + non-essential + essential" begin
            forms = enumerate_enzyme_forms(uni_uni_act)
            specs = EnzymeRates._enumerate_only_catalytic_mechanisms(
                forms, uni_uni_act; max_forms=10,
            )
            input = specs[1]
            result = EnzymeRates._generate_activator_configs(
                input, forms, uni_uni_act,
            )
            @test length(result) == 3

            # Config 1: absent (same as input)
            @test result[1].reactions == input.reactions

            # All shadow variants have more reactions than input
            @test all(
                length(r.reactions) > length(input.reactions)
                for r in result[2:end]
            )
            @test all(r.forms == input.forms for r in result)

            # Config 2: non-essential shadow
            @test :E_0_0_A ∈ _used_forms(result[2])
            @test :E_S_0_A ∈ _used_forms(result[2])
            @test :E_0_P_A ∈ _used_forms(result[2])
            base_edges = Set(EnzymeRates._spec_to_edges(input, forms))
            ne_edges = Set(EnzymeRates._spec_to_edges(result[2], forms))
            @test base_edges ⊆ ne_edges  # base cycle retained

            # Config 3: essential shadow
            @test :E_0_0_A ∈ _used_forms(result[3])
            @test :E_S_0_A ∈ _used_forms(result[3])
            @test :E_0_P_A ∈ _used_forms(result[3])
            ess_edges = Set(EnzymeRates._spec_to_edges(result[3], forms))
            @test !issubset(base_edges, ess_edges)  # base cycle removed
        end

        @testset "Essential activator: valid rate equation" begin
            forms = enumerate_enzyme_forms(uni_uni_act)
            specs = EnzymeRates._enumerate_only_catalytic_mechanisms(
                forms, uni_uni_act; max_forms=10,
            )
            result = EnzymeRates._generate_activator_configs(
                specs[1], forms, uni_uni_act,
            )
            # Convert essential variant (config 3) to EnzymeMechanism
            ess_spec = result[3]
            mech = EnzymeMechanism(ess_spec)
            req = rate_equation_string(mech)
            @test occursin("A", req)  # activator in rate equation
        end
    end

    # ─── Section 4: _enumerate_dead_end_configs ────────────

    @testset "_enumerate_dead_end_configs" begin
        @testset "No regulatory sites: input only" begin
            forms = enumerate_enzyme_forms(uni_uni)
            specs = EnzymeRates._enumerate_only_catalytic_mechanisms(
                forms, uni_uni; max_forms=6,
            )
            input = specs[1]
            result = EnzymeRates._enumerate_dead_end_configs(
                input, forms; max_forms=6,
            )
            @test result[1].reactions == input.reactions
        end

        @testset "Uni-Uni + inhibitor: dead-end configs" begin
            forms = enumerate_enzyme_forms(uni_uni_inh)
            specs = EnzymeRates._enumerate_only_catalytic_mechanisms(
                forms, uni_uni_inh; max_forms=10,
            )
            input = specs[1]
            result = EnzymeRates._enumerate_dead_end_configs(
                input, forms; max_forms=10,
            )
            @test result[1].reactions == input.reactions
            @test EnzymeRates._used_form_count(result[1]) == 3
            @test any(
                s -> :E_0_0_I ∈ _used_forms(s), result[2:end],
            )
            # Dead-end forms are reachable from topology
            topo_used = _used_forms(input)
            @test all(result[2:end]) do spec
                used = _used_forms(spec)
                de_forms = setdiff(used, topo_used)
                filter!(f -> f ∈ Set(spec.forms), de_forms)
                !isempty(de_forms)
            end
            @test all(
                EnzymeRates._used_form_count(s) <= 10
                for s in result
            )
        end

        @testset "Uni-Uni + 2 regulators: box rule" begin
            forms = enumerate_enzyme_forms(uni_uni_2reg)
            specs = EnzymeRates._enumerate_only_catalytic_mechanisms(
                forms, uni_uni_2reg; max_forms=10,
            )
            input = specs[1]
            result = EnzymeRates._enumerate_dead_end_configs(
                input, forms; max_forms=10,
            )
            @test result[1].reactions == input.reactions
            @test any(
                s -> :E_0_0_I_J ∈ _used_forms(s), result,
            )
            # Box rule: dual-reg form requires BOTH single-reg parents
            @test all(result) do spec
                used = _used_forms(spec)
                :E_0_0_I_J ∉ used ||
                    (:E_0_0_I_0 ∈ used && :E_0_0_0_J ∈ used)
            end
        end
    end

    # ─── Section 5: _enumerate_ress_and_constraints ────────

    @testset "_enumerate_ress_and_constraints" begin
        @testset "Uni-Uni: no equiv groups" begin
            forms = enumerate_enzyme_forms(uni_uni)
            specs = EnzymeRates._enumerate_only_catalytic_mechanisms(
                forms, uni_uni; max_forms=6,
            )
            result = EnzymeRates._enumerate_ress_and_constraints(
                specs[1], forms,
            )
            @test all(
                isempty(r.param_constraints) for r in result
            )
            @test all(
                r -> any(==(false), r.equilibrium_steps),
                result,
            )
            @test any(
                r -> all(==(false), r.equilibrium_steps),
                result,
            )
        end

        @testset "Dead-ends: equiv groups create constraints" begin
            forms = enumerate_enzyme_forms(uni_uni_inh)
            specs = EnzymeRates._enumerate_only_catalytic_mechanisms(
                forms, uni_uni_inh; max_forms=10,
            )
            de_specs = EnzymeRates._enumerate_dead_end_configs(
                specs[1], forms; max_forms=10,
            )
            # Find config with >= 2 inhibitor binding steps
            multi_inh = first(
                s for s in de_specs
                if length(s.reactions) >= 6
            )
            result = EnzymeRates._enumerate_ress_and_constraints(
                multi_inh, forms,
            )
            @test any(
                r -> !isempty(r.param_constraints), result,
            )
            @test any(
                r -> isempty(r.param_constraints), result,
            )
            @test all(result) do r
                all(r.param_constraints) do (target, coeff, factors)
                    target isa Symbol &&
                        coeff > 0 &&
                        !isempty(factors)
                end
            end
        end
    end

    # ─── Section 6: enumerate_mechanisms (end-to-end) ──────────

    @testset "enumerate_mechanisms (end-to-end)" begin
        @testset "Uni-Uni: length matches iteration" begin
            iter = enumerate_mechanisms(uni_uni)
            @test length(iter) == length(collect(iter))
        end

        @testset "Bi-Bi: ordered and random-order" begin
            mechs6 = collect(
                enumerate_mechanisms(bi_bi; max_forms=6),
            )
            @test any(mechs6) do spec
                EnzymeRates._used_form_count(spec) == 5 &&
                    !any(spec.equilibrium_steps)
            end

            mechs7 = collect(
                enumerate_mechanisms(bi_bi; max_forms=7),
            )
            @test any(mechs7) do spec
                EnzymeRates._used_form_count(spec) >= 7 &&
                    !any(spec.equilibrium_steps)
            end
        end

        @testset "Regulator dead-end inhibitor" begin
            mechs = collect(
                enumerate_mechanisms(uni_uni_inh; max_forms=6),
            )
            @test any(
                s -> :E_0_0_I ∈ _used_forms(s), mechs,
            )
        end

        @testset "Activator mechanism" begin
            mechs = collect(
                enumerate_mechanisms(uni_uni_act; max_forms=8),
            )
            # Activator cycle: shadow forms in reactions
            @test any(mechs) do spec
                rf = _used_forms(spec)
                :E_0_0_A ∈ rf && :E_S_0_A ∈ rf
            end
            # Dead-end-only variants
            @test any(mechs) do spec
                rf = _used_forms(spec)
                :E_0_0_A ∈ rf && :E_S_0_A ∉ rf
            end
        end

        @testset "All Uni-Uni compile rate equations" begin
            @test all(enumerate_mechanisms(uni_uni)) do spec
                m = EnzymeMechanism(spec)
                s = rate_equation_string(m)
                s isa String && !isempty(s)
            end
        end

        @testset "max_forms limit" begin
            small = collect(
                enumerate_mechanisms(uni_uni_inh; max_forms=4),
            )
            large = collect(
                enumerate_mechanisms(uni_uni_inh; max_forms=8),
            )
            @test length(small) <= length(large)
            @test all(
                EnzymeRates._used_form_count(s) <= 4
                for s in small
            )
            @test all(
                EnzymeRates._used_form_count(s) <= 8
                for s in large
            )
        end

        @testset "Equivalent step constraints" begin
            mechs = collect(
                enumerate_mechanisms(uni_uni_inh; max_forms=6),
            )
            @test any(
                s -> !isempty(s.param_constraints), mechs,
            )
            @test any(
                s -> isempty(s.param_constraints), mechs,
            )
        end

        @testset "Ping-pong: no invalid empty->residual edges" begin
            rxn1 = @enzyme_reaction begin
                substrates: Glu[C6H12O6], ATP[C10H16N5O13P3]
                products: G6P[C6H13O9P], ADP[C10H15N5O10P2]
            end
            rxn2 = @enzyme_reaction begin
                substrates: A[C], B[C2]
                products: P[C2], Q[C]
            end
            @test all([rxn1, rxn2]) do rxn
                iter = enumerate_mechanisms(rxn; max_forms=6)
                all(
                    spec -> EnzymeMechanism(spec) isa EnzymeMechanism,
                    iter,
                )
            end
        end

        @testset "Regulators: correctness" begin
            rxn_reg = @enzyme_reaction begin
                substrates: Glu[C6H12O6], ATP[C10H16N5O13P3]
                products: G6P[C6H13O9P], ADP[C10H15N5O10P2]
                regulators: Phosphate[PO4], G6P[C6H13O9P], G6P[C6H13O9P]
            end
            rxn_no_reg = @enzyme_reaction begin
                substrates: Glu[C6H12O6], ATP[C10H16N5O13P3]
                products: G6P[C6H13O9P], ADP[C10H15N5O10P2]
            end
            mechs_reg = collect(
                enumerate_mechanisms(rxn_reg; max_forms=5),
            )
            mechs_no_reg = collect(
                enumerate_mechanisms(rxn_no_reg; max_forms=5),
            )
            @test length(mechs_reg) == length(mechs_no_reg)
        end
    end

    # ─── Section 7: edge_class and helpers ─────────────────

    @testset "edge_class and helpers" begin
        @testset "Core binding -> MustExist" begin
            forms = enumerate_enzyme_forms(uni_uni)
            fe = first(f for f in forms if f.name == :E_0_0)
            fs = first(f for f in forms if f.name == :E_S_0)
            ec, met, etype = EnzymeRates.edge_class(fe, fs)
            @test ec isa EnzymeRates.MustExist
            @test met == :S
            @test etype == :binding
        end

        @testset "Isomerization -> MustExist" begin
            forms = enumerate_enzyme_forms(uni_uni)
            fs = first(f for f in forms if f.name == :E_S_0)
            fp = first(f for f in forms if f.name == :E_0_P)
            ec, met, etype = EnzymeRates.edge_class(fs, fp)
            @test ec isa EnzymeRates.MustExist
            @test met === nothing
            @test etype == :isomerization
        end

        @testset "Same form -> Forbidden" begin
            forms = enumerate_enzyme_forms(uni_uni)
            fe = first(f for f in forms if f.name == :E_0_0)
            ec, _, _ = EnzymeRates.edge_class(fe, fe)
            @test ec isa EnzymeRates.Forbidden
        end

        @testset "Regulator binding -> CouldExist" begin
            forms = enumerate_enzyme_forms(uni_uni_inh)
            fe = first(
                f for f in forms if f.name == :E_0_0_0
            )
            fi = first(
                f for f in forms if f.name == :E_0_0_I
            )
            ec, met, etype = EnzymeRates.edge_class(fe, fi)
            @test ec isa EnzymeRates.CouldExist
            @test met == :I
            @test etype == :binding
        end

        @testset "n_sites helper" begin
            @test n_sites(uni_uni) == 2
            r2 = @enzyme_reaction begin
                substrates: A[C, 2], B[N]
                products:   P[C], Q[N]
                regulators: I[N2, 3]
            end
            @test n_sites(r2) == 2 + 1 + 1 + 1 + 3
        end
    end
end
