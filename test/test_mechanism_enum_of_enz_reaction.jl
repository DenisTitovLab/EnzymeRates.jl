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

        with_act = EnzymeRates.enumerate_mechanisms(
            spec.reaction;
            stage=EnzymeRates.WithActivator(),
            max_forms=spec.max_forms,
        )
        @test length(with_act) == spec.expected_n_cat_with_act

        with_de = collect(EnzymeRates.enumerate_mechanisms(
            spec.reaction;
            stage=EnzymeRates.WithDeadEnd(),
            max_forms=spec.max_forms,
        ))
        @test length(with_de) == spec.expected_n_cat_act_de

        # Independent dead-end verification:
        # (2^r_inh)^n_topo per activator config, summed
        expected_de_total = _compute_expected_dead_end_count(
            with_act, forms)
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

        # Rate equation smoke test
        if spec.test_rate_equation
            @test any(catalytic) do s
                m = EnzymeMechanism(s)
                rate_equation_string(m) isa String
            end
        end
    end
end
