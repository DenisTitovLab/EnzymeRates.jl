@testset "Mechanism Enumeration" begin
    @testset "Pipeline: $(spec.name)" for spec in ENUMERATION_TEST_SPECS
        forms = EnzymeRates.enumerate_enzyme_forms(spec.reaction)
        @test length(forms) == spec.expected_n_forms

        catalytic = EnzymeRates.enumerate_mechanisms(
            spec.reaction;
            stage=EnzymeRates.Catalytic(),
        )
        @test length(catalytic) == spec.expected_n_catalytic

        with_de = collect(EnzymeRates.enumerate_mechanisms(
            spec.reaction;
            stage=EnzymeRates.WithDeadEnd(),
        ))
        @test length(with_de) == spec.expected_n_cat_de

        # Independent dead-end verification:
        # (2^r_inh)^n_topo per catalytic topology, summed
        expected_de_total = _compute_expected_dead_end_count(
            catalytic, forms)
        @test expected_de_total == length(with_de)

        # Total mechanism count (O(1) for lazy iterator)
        final = EnzymeRates.enumerate_mechanisms(
            spec.reaction)
        @test length(final) == spec.expected_n_total

        # RE/SS + constraints (closed-form formula verification)
        if !spec.skip_ress_test
            expected_ress_total = sum(with_de) do base
                _compute_expected_n_total(base, forms)
            end
            @test expected_ress_total == length(final)

            if isfinite(spec.max_enumeration_time)
                t = @elapsed EnzymeRates.enumerate_mechanisms(
                    spec.reaction)
                @test t < spec.max_enumeration_time
            end
        end

        # Verify catalytic topologies are subsets of dead-end specs.
        # Dead-end expansion may add binding edges between topology
        # forms, so we check subset (⊆) not exact match.
        @testset "Stage subset" begin
            de_edge_sets = [Set(s.edges) for s in with_de]
            for s in catalytic
                cat_set = Set(s.edges)
                @test any(de -> cat_set ⊆ de, de_edge_sets)
            end
        end

        # Regulator partition verification
        @testset "Regulator partitioning" begin
            regs = collect(Symbol,
                EnzymeRates.regulators(spec.reaction))
            n_reg = length(regs)
            # Every dead-end spec must have valid partition info
            for s in with_de
                @test sort(vcat(s.dead_end_regulators,
                    s.allosteric_regulators)) == sort(regs)
            end
            # Number of distinct partitions = 2^n_reg
            partitions = Set(
                (Tuple(sort(s.dead_end_regulators)),
                 Tuple(sort(s.allosteric_regulators)))
                for s in with_de)
            @test length(partitions) == 1 << n_reg
        end

        # EnzymeMechanism construction from MechanismSpec.
        # Use first 10 (simplest) — large mechanisms compile slowly
        # in @generated rate_equation due to LLVM register allocation.
        mech_test_specs = first(with_de, 10)
        mechanisms = EnzymeMechanism[]
        @testset "EnzymeMechanism construction" begin
            t1 = @elapsed for s in mech_test_specs
                push!(mechanisms, EnzymeMechanism(s))
            end
            @test length(mechanisms) == length(mech_test_specs)
            @test t1 < 10.0

            mechanisms2 = EnzymeMechanism[]
            t2 = @elapsed for s in mech_test_specs
                push!(mechanisms2, EnzymeMechanism(s))
            end
            @test t2 < 1.0
        end

        # Rate equation functions on constructed mechanisms
        @testset "Rate equation functions" begin
            for m in mechanisms
                metabs = metabolites(m)
                @test metabs isa Tuple{Vararg{Symbol}}
                @test length(metabs) > 0

                params_tup = parameters(m)
                @test params_tup isa Tuple{Vararg{Symbol}}
                @test :E_total ∈ params_tup
                @test :Keq ∈ params_tup

                concs = NamedTuple{metabs}(ones(length(metabs)))
                pvals = NamedTuple{params_tup}(
                    ones(length(params_tup)))
                v = rate_equation(m, concs, pvals)
                @test v isa Real
                @test isfinite(v)

                s = rate_equation_string(m)
                @test s isa String
                @test length(s) > 0
            end
        end

        # compile_mechanism dispatches correctly (exported)
        @testset "compile_mechanism" begin
            for s in mech_test_specs
                em = compile_mechanism(s)
                @test em isa EnzymeMechanism
            end
        end

        # Oligomeric expansion (catalytic_n=2)
        @testset "Oligomeric expansion" begin
            cat_n = 2
            final_oligo = EnzymeRates.enumerate_mechanisms(
                spec.reaction; catalytic_n=cat_n)

            # Verify count formula
            if !spec.skip_ress_test
                expected_oligo = _compute_expected_oligomeric_total(
                    with_de, forms, cat_n)
                @test length(final_oligo) == expected_oligo
            end

            # Verify EM + OEM mix: iterate a few specs
            sample = Iterators.take(final_oligo, 20)
            em_count = 0
            oem_count = 0
            for s in sample
                if s.catalytic_n == 0
                    em_count += 1
                else
                    oem_count += 1
                    @test s.catalytic_n == cat_n
                    @test s.n_conf == 2
                    k = length(s.allosteric_regulators)
                    @test length(s.allosteric_multiplicities) == k
                    if k > 0
                        @test all(
                            1 .<= s.allosteric_multiplicities .<= cat_n)
                    end
                end
            end
            @test em_count > 0
            @test oem_count > 0

            # compile_mechanism produces OligomericEnzymeMechanism
            # for the first OEM spec found
            oem_spec = nothing
            for s in Iterators.take(final_oligo, 50)
                if s.catalytic_n > 0
                    oem_spec = s
                    break
                end
            end
            if oem_spec !== nothing
                oem = EnzymeRates.compile_mechanism(oem_spec)
                @test oem isa OligomericEnzymeMechanism

                # Verify type parameters
                metabs = metabolites(oem)
                @test metabs isa Tuple{Vararg{Symbol}}
                @test length(metabs) > 0

                params_tup = parameters(oem)
                @test :E_total ∈ params_tup
                @test :Keq ∈ params_tup
                @test :L ∈ params_tup  # NConf=2 → L param

                # Rate equation evaluates
                concs = NamedTuple{metabs}(
                    ones(length(metabs)))
                pvals = NamedTuple{params_tup}(
                    ones(length(params_tup)))
                v = rate_equation(oem, concs, pvals)
                @test v isa Real
                @test isfinite(v)
            end
        end
    end

    # ── Integration tests ──────────────────────────────────────────
    # Use Uni-Uni (simplest, fast compilation) for detailed checks.
    @testset "Integration: G≤7, SS step, round-trip" begin
        rxn = @enzyme_reaction begin
            substrates:S[C]
            products:P[C]
        end
        iter = EnzymeRates.enumerate_mechanisms(rxn)
        specs = collect(iter)
        @test length(specs) > 0

        for s in specs
            # Every mechanism must have at least 1 SS step
            @test any(.!s.equilibrium_steps)
            # G must be in [2, 7]
            G = _compute_re_group_count_test(
                s.edges, s.equilibrium_steps)
            @test 2 <= G <= 7

            # Round-trip: compile → rate_equation
            m = compile_mechanism(s)
            @test m isa EnzymeMechanism
            metabs = metabolites(m)
            params_tup = parameters(m)
            concs = NamedTuple{metabs}(ones(length(metabs)))
            pvals = NamedTuple{params_tup}(
                ones(length(params_tup)))
            v = rate_equation(m, concs, pvals)
            @test isfinite(v)
        end
    end

    @testset "Integration: OEM round-trip + NConf=2" begin
        rxn = @enzyme_reaction begin
            substrates:S[C]
            products:P[C]
            regulators: R
        end
        iter = EnzymeRates.enumerate_mechanisms(
            rxn; catalytic_n=2)
        # Check first few OEM specs
        oem_checked = 0
        for s in Iterators.take(iter, 100)
            s.catalytic_n == 0 && continue
            @test s.n_conf == 2
            m = compile_mechanism(s)
            @test m isa OligomericEnzymeMechanism
            metabs = metabolites(m)
            params_tup = parameters(m)
            @test :L ∈ params_tup
            concs = NamedTuple{metabs}(ones(length(metabs)))
            pvals = NamedTuple{params_tup}(
                ones(length(params_tup)))
            v = rate_equation(m, concs, pvals)
            @test isfinite(v)
            oem_checked += 1
            oem_checked >= 3 && break
        end
        @test oem_checked >= 1
    end

    @testset "max_re_groups parameter" begin
        rxn = @enzyme_reaction begin
            substrates:S[C]
            products:P[C]
        end
        iter_5 = EnzymeRates.enumerate_mechanisms(
            rxn; max_re_groups=5)
        iter_7 = EnzymeRates.enumerate_mechanisms(rxn)
        # Tighter cap → fewer or equal mechanisms
        @test length(iter_5) <= length(iter_7)
        # All specs from iter_5 satisfy G ≤ 5
        for s in iter_5
            G = _compute_re_group_count_test(
                s.edges, s.equilibrium_steps)
            @test 2 <= G <= 5
        end
    end
end
