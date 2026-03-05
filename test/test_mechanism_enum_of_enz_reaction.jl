@testset "Mechanism Enumeration" begin
    @testset "Pipeline: $(spec.name)" for spec in ENUMERATION_TEST_SPECS
        forms = EnzymeRates.enumerate_enzyme_forms(spec.reaction)
        @test length(forms) == spec.expected_n_forms

        catalytic = EnzymeRates.enumerate_mechanisms(
            spec.reaction;
            stage=EnzymeRates.Catalytic(),
            max_forms=spec.max_forms,
        )
        @test length(catalytic) == spec.expected_n_catalytic

        with_de = collect(EnzymeRates.enumerate_mechanisms(
            spec.reaction;
            stage=EnzymeRates.WithDeadEnd(),
            max_forms=spec.max_forms,
        ))
        @test length(with_de) == spec.expected_n_cat_de

        # Independent dead-end verification:
        # (2^r_inh)^n_topo per catalytic topology, summed
        expected_de_total = _compute_expected_dead_end_count(
            catalytic, forms)
        @test expected_de_total == length(with_de)

        # Total mechanism count (O(1) for lazy iterator)
        final = EnzymeRates.enumerate_mechanisms(
            spec.reaction; max_forms=spec.max_forms)
        @test length(final) == spec.expected_n_total

        # RE/SS + constraints (closed-form formula verification)
        if !spec.skip_ress_test
            expected_ress_total = sum(with_de) do base
                _compute_expected_n_total(base, forms)
            end
            @test expected_ress_total == length(final)

            if isfinite(spec.max_enumeration_time)
                t = @elapsed EnzymeRates.enumerate_mechanisms(
                    spec.reaction; max_forms=spec.max_forms)
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
    end
end
