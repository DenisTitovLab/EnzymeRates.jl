@testset "Mechanism Enumeration" begin
    @testset "Pipeline: $(spec.name)" for spec in ENUMERATION_TEST_SPECS
        stages = enumerate_mechanism_stages(
            spec.reaction; max_forms=spec.max_forms)

        # Stage counts
        @test length(stages.forms) == spec.expected_n_forms
        @test length(stages.catalytic) == spec.expected_n_catalytic
        @test length(stages.with_activator) == spec.expected_n_cat_with_act
        @test length(stages.with_dead_end) == spec.expected_n_cat_act_de

        # Independent dead-end verification (aggregate):
        # sum of per-base independent counts == total dead-end count
        expected_de_total = sum(stages.with_activator) do base
            _compute_expected_dead_end_count(
                base, stages.forms, spec.max_forms)
        end
        @test expected_de_total == length(stages.with_dead_end)

        # Total mechanism count (O(1) for lazy iterator)
        @test length(stages.final) == spec.expected_n_total

        # RE/SS + constraints (expensive independent verification)
        if !spec.skip_ress_test
            expected_ress_total = sum(stages.with_dead_end) do base
                _compute_independent_ress_count(base, stages.forms)
            end
            @test expected_ress_total == length(stages.final)

            t = @elapsed n_total = length(enumerate_mechanisms(
                spec.reaction; max_forms=spec.max_forms))
            @test n_total == spec.expected_n_total

            if isfinite(spec.max_enumeration_time)
                @test t < spec.max_enumeration_time
            end
        end

        # Rate equation smoke test
        if spec.test_rate_equation
            @test any(stages.catalytic) do s
                m = EnzymeMechanism(s)
                rate_equation_string(m) isa String
            end
        end
    end
end
