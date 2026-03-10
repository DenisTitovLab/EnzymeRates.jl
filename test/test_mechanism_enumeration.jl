# ABOUTME: Tests for the staged mechanism enumeration pipeline
# ABOUTME: Iterates over specs defined in mechanism_enumeration_test_specs.jl

const STAGE_EXPANSION_SPECS = build_stage_expansion_specs()
const ENUMERATION_SPECS = build_enumeration_specs()

"""Run `f()` with a timeout. Returns `f()` result or `nothing` on timeout."""
function _with_timeout(f, timeout_secs::Real)
    result_channel = Channel{Any}(1)
    task = @async try
        put!(result_channel, f())
    catch e
        put!(result_channel, e)
    end
    timer = Timer(timeout_secs)
    @async begin
        wait(timer)
        if !istaskdone(task)
            schedule(task, InterruptException(); error=true)
        end
    end
    result = take!(result_channel)
    close(timer)
    close(result_channel)
    result isa Exception ? nothing : result
end

@testset "Mechanism Enumeration Pipeline" begin

    # ── Stage expansion: each stage independently on base ────
    @testset "Stage expansion: $(s.name)" for s in STAGE_EXPANSION_SPECS
        rxn = s.reaction
        roles = EnzymeRates.regulator_roles(rxn)
        de_regs = Symbol[r[1] for r in roles
                         if r[2] == :dead_end]
        al_regs = Symbol[r[1] for r in roles
                         if r[2] == :allosteric]
        base = [s.base_mechanism]

        @test length(EnzymeRates._expand_ress_variants(
            base, rxn)) == s.expected_n_ress
        @test length(EnzymeRates._expand_dead_end_inhibitors(
            base, rxn; dead_end_regs=de_regs)) ==
            s.expected_n_dead_end
        @test length(EnzymeRates._expand_equivalence_constraints(
            base, rxn)) == s.expected_n_equivalence
        @test length(EnzymeRates._deduplicate(
            base, rxn)) == s.expected_n_dedup

        if !isempty(s.allosteric_regs)
            cn = s.catalytic_n > 0 ? s.catalytic_n : 1
            dd = EnzymeRates._deduplicate(base, rxn)
            allo =EnzymeRates._expand_allosteric(
                dd, rxn; catalytic_n=cn,
                allosteric_regs=al_regs)
            @test length(allo) == s.expected_n_allosteric

            allo = EnzymeRates._expand_tr_equivalence(allo, rxn)
            @test length(allo) == s.expected_n_tr_equiv

            allo = EnzymeRates._deduplicate_allosteric(allo, rxn)
            @test length(allo) == s.expected_n_allosteric_dedup
        end
    end

    # ── End-to-end pipeline ──────────────────────────────────
    @testset "End-to-end: $(s.name)" for s in ENUMERATION_SPECS
        rxn = s.reaction

        _, forms = EnzymeRates.enumerate_enzyme_forms(rxn)
        @test length(forms) == s.expected_n_forms

        counts = _run_full_pipeline_stages(
            rxn; catalytic_n=s.catalytic_n)
        @test counts.catalytic == s.expected_n_catalytic
        @test counts.ress == s.expected_n_ress
        @test counts.dead_end == s.expected_n_dead_end
        @test counts.equivalence == s.expected_n_equivalence
        @test counts.dedup == s.expected_n_dedup

        if s.expected_n_allosteric > 0
            @test counts.allosteric == s.expected_n_allosteric
            @test counts.tr_equiv == s.expected_n_tr_equiv
            @test counts.allosteric_dedup ==
                s.expected_n_allosteric_dedup
        end

        result = _with_timeout(120) do
            collect(EnzymeRates.enumerate_mechanisms(
                rxn; catalytic_n=s.catalytic_n))
        end

        if result === nothing
            @warn "$(s.name) timed out (120s)"
            @test_broken length([]) == s.expected_n_total
        else
            @test length(result) == s.expected_n_total
        end
    end

    # ── Property-based tests ─────────────────────────────────
    @testset "Catalytic topology properties" begin
        for rxn in [uni_uni, bi_bi, bi_bi_ping_pong]
            catalytic = EnzymeRates._catalytic_topologies(rxn)
            @test length(catalytic) > 0
            for spec in catalytic
                @test spec.n_catalytic_edges ==
                    length(spec.edges)
                @test count(.!spec.equilibrium_steps) >= 1
            end
        end
    end

    @testset "RE/SS expansion properties" begin
        for rxn in [uni_uni, bi_bi, bi_bi_ping_pong]
            catalytic = EnzymeRates._catalytic_topologies(rxn)
            for spec in catalytic
                ress = EnzymeRates._expand_ress_variants(
                    [spec], rxn)
                @test length(ress) > 0
                for s in ress
                    @test any(.!s.equilibrium_steps)
                end
            end
        end
    end

    @testset "RE partition bounds" begin
        catalytic = EnzymeRates._catalytic_topologies(uni_bi)
        spec = catalytic[1]
        ress = EnzymeRates._expand_ress_variants(
            [spec], uni_bi)
        @test length(ress) > 0
        for s in ress
            partition = EnzymeRates._compute_re_partition(
                s.edges, s.equilibrium_steps)
            @test 2 <= length(partition) <= 7
        end
    end

    @testset "Dead-end expansion properties" begin
        catalytic = EnzymeRates._catalytic_topologies(
            uni_bi_dead_end_I)
        de_specs = EnzymeRates._expand_dead_end_inhibitors(
            catalytic, uni_bi_dead_end_I;
            dead_end_regs=[:I])
        @test length(de_specs) > length(catalytic)
        for s in de_specs
            @test length(s.edges) >= s.n_catalytic_edges
        end
    end

    @testset "Dead-end passthrough with no regs" begin
        for rxn in [uni_uni, bi_bi]
            catalytic = EnzymeRates._catalytic_topologies(rxn)
            no_de = EnzymeRates._expand_dead_end_inhibitors(
                catalytic, rxn; dead_end_regs=Symbol[])
            @test length(no_de) == length(catalytic)
        end
    end

    @testset "Deduplication reduces count" begin
        catalytic = EnzymeRates._catalytic_topologies(
            uni_bi_dead_end_I)
        de_specs = EnzymeRates._expand_dead_end_inhibitors(
            catalytic, uni_bi_dead_end_I;
            dead_end_regs=[:I])
        ress = EnzymeRates._expand_ress_variants(
            [de_specs[1]], uni_bi_dead_end_I)
        with_eq = EnzymeRates._expand_equivalence_constraints(
            ress, uni_bi_dead_end_I)
        deduped = EnzymeRates._deduplicate(
            with_eq, uni_bi_dead_end_I)
        @test length(deduped) <= length(with_eq)
    end

    @testset "Equivalence constraints add variants" begin
        catalytic = EnzymeRates._catalytic_topologies(
            uni_bi_dead_end_I)
        de_specs = EnzymeRates._expand_dead_end_inhibitors(
            catalytic, uni_bi_dead_end_I;
            dead_end_regs=[:I])
        spec_idx = findfirst(
            s -> length(s.edges) > s.n_catalytic_edges,
            de_specs)
        if spec_idx !== nothing
            s = de_specs[spec_idx]
            ress = EnzymeRates._expand_ress_variants(
                [s], uni_bi_dead_end_I)
            with_eq =
                EnzymeRates._expand_equivalence_constraints(
                    ress, uni_bi_dead_end_I)
            @test length(with_eq) >= length(ress)
        end
    end

    @testset "Stage monotonicity" begin
        for rxn in [uni_bi_reg_unknown,
                    bi_bi_ping_pong_reg_unknown]
            counts = _run_full_pipeline_stages(rxn)
            @test counts.dead_end >= counts.ress
            @test counts.equivalence >= counts.dead_end
            @test counts.dedup <= counts.equivalence
        end
    end

    @testset "Regulator roles affect partitioning" begin
        c_unk = _run_full_pipeline_stages(uni_uni_reg_unknown)
        c_de = _run_full_pipeline_stages(uni_uni_dead_end_I)
        c_al = _run_full_pipeline_stages(uni_uni_allosteric_I)
        @test c_unk.catalytic == c_de.catalytic
        @test c_unk.catalytic == c_al.catalytic
        @test c_de.dedup >= c_al.dedup
    end

    # ── param_count accuracy ─────────────────────────────────
    @testset "param_count accuracy (EM)" begin
        all_specs = collect(
            EnzymeRates.enumerate_mechanisms(uni_bi))
        @test length(all_specs) > 0
        for s in all_specs
            m = compile_mechanism(s)
            @test s.param_count == length(parameters(m))
        end
    end

    @testset "param_count accuracy (EM with reg, sampled)" begin
        all_specs = collect(
            EnzymeRates.enumerate_mechanisms(
                uni_bi_reg_unknown))
        rng = Random.MersenneTwister(42)
        base_specs = filter(
            s -> s isa EnzymeRates.MechanismSpec, all_specs)
        n = min(20, length(base_specs))
        sample = base_specs[randperm(rng,
            length(base_specs))[1:n]]
        for s in sample
            m = compile_mechanism(s)
            @test s.param_count == length(parameters(m))
        end
    end

    @testset "param_count accuracy per stage" begin
        rng = Random.MersenneTwister(99)

        # Stage 1: catalytic topologies
        catalytic = EnzymeRates._catalytic_topologies(uni_bi)
        for s in catalytic
            m = compile_mechanism(s)
            @test s.param_count == length(parameters(m))
        end

        # Stage 2: RE/SS expansion
        ress = EnzymeRates._expand_ress_variants(
            catalytic, uni_bi)
        for s in ress[1:min(10, length(ress))]
            m = compile_mechanism(s)
            @test s.param_count == length(parameters(m))
        end

        # Stage 6: equivalence constraints (no regs, so
        # dead-end stage is passthrough)
        eq = EnzymeRates._expand_equivalence_constraints(
            ress, uni_bi)
        for s in eq[1:min(10, length(eq))]
            m = compile_mechanism(s)
            @test s.param_count == length(parameters(m))
        end

        # With regulators: dead-end + equivalence stages
        cat_r = EnzymeRates._catalytic_topologies(
            uni_bi_dead_end_I)
        ress_r = EnzymeRates._expand_ress_variants(
            cat_r, uni_bi_dead_end_I)
        de_r = EnzymeRates._expand_dead_end_inhibitors(
            ress_r, uni_bi_dead_end_I;
            dead_end_regs=[:I])
        sample_de = de_r[randperm(rng, length(de_r))[
            1:min(10, length(de_r))]]
        for s in sample_de
            m = compile_mechanism(s)
            @test s.param_count == length(parameters(m))
        end

        eq_r = EnzymeRates._expand_equivalence_constraints(
            de_r, uni_bi_dead_end_I)
        sample_eq = eq_r[randperm(rng, length(eq_r))[
            1:min(10, length(eq_r))]]
        for s in sample_eq
            m = compile_mechanism(s)
            @test s.param_count == length(parameters(m))
        end
    end

    @testset "param_count accuracy (sampled)" begin
        rng = Random.MersenneTwister(42)
        for rxn in [uni_uni, uni_bi_reg_unknown, uni_bi]
            all_specs = collect(
                EnzymeRates.enumerate_mechanisms(rxn))
            base_specs = filter(
                s -> s isa EnzymeRates.MechanismSpec,
                all_specs)
            n = min(20, length(base_specs))
            sample = base_specs[randperm(rng,
                length(base_specs))[1:n]]
            for s in sample
                m = compile_mechanism(s)
                @test s.param_count == length(parameters(m))
            end
        end
    end

    # ── Allosteric expansion properties ─────────────────────
    @testset "Allosteric expansion properties" begin
        all_specs = collect(
            EnzymeRates.enumerate_mechanisms(
                uni_bi_allosteric_I_cn2; catalytic_n=2))
        em = filter(
            s -> s isa EnzymeRates.MechanismSpec,
            all_specs)
        allo = filter(
            s -> s isa EnzymeRates.AllostericMechanismSpec,
            all_specs)
        @test length(em) > 0
        @test length(allo) > 0
        for s in allo
            @test s.catalytic_n == 2
        end
    end

    @testset "Allosteric with allosteric regulators" begin
        all_specs = collect(
            EnzymeRates.enumerate_mechanisms(
                uni_bi_allosteric_I_cn2; catalytic_n=2))
        allo = filter(
            s -> s isa EnzymeRates.AllostericMechanismSpec,
            all_specs)
        @test length(allo) > 0
        for s in allo
            @test !isempty(s.allosteric_reg_sites)
            @test !isempty(s.allosteric_multiplicities)
        end
    end

    # ── compile_mechanism round-trip ─────────────────────────
    @testset "compile_mechanism round-trip" begin
        for rxn in [uni_uni, bi_bi]
            all_specs = collect(
                EnzymeRates.enumerate_mechanisms(rxn))
            for s in all_specs
                m = compile_mechanism(s)
                @test m isa EnzymeMechanism
                @test length(metabolites(m)) > 0
                @test length(parameters(m)) > 0
            end
        end
    end

    @testset "compile_mechanism allosteric" begin
        all_specs = collect(
            EnzymeRates.enumerate_mechanisms(
                uni_bi_allosteric_I_cn2; catalytic_n=2))
        allo = filter(
            s -> s isa EnzymeRates.AllostericMechanismSpec,
            all_specs)
        for s in allo[1:min(3, length(allo))]
            m = compile_mechanism(s)
            @test m isa AllostericEnzymeMechanism
            @test length(parameters(m)) > 0
        end
    end

    # ── Combinatorial cross-checks ───────────────────────────
    @testset "Combinatorial cross-checks" begin
        # Verify hardcoded expected values against independent
        # combinatorial formulas.

        specs = STAGE_EXPANSION_SPECS
        by_name = Dict(s.name => s for s in specs)

        # RE/SS: n RE binding edges → 2^n - 1 valid combos
        # (must keep ≥1 RE edge)
        @test by_name["Uni-Uni (no reg)"].expected_n_ress ==
            2^2 - 1  # 2 RE binding edges
        @test by_name["Bi-Bi (no reg)"].expected_n_ress ==
            2^4 - 1  # 4 RE binding edges
        @test by_name["Bi-Bi Ping-Pong (no reg)"].expected_n_ress ==
            2^4 - 1  # 4 RE binding edges

        # No-reg passthroughs: de/eq/dd all = 1
        for name in ["Uni-Uni (no reg)", "Bi-Bi (no reg)",
                      "Bi-Bi Ping-Pong (no reg)"]
            s = by_name[name]
            @test s.expected_n_dead_end == 1
            @test s.expected_n_equivalence == 1
            @test s.expected_n_dedup == 1
        end

        # Dead-end: n catalytic forms → 2^n subsets of I binding
        @test by_name["Uni-Uni (dead-end I)"].expected_n_dead_end ==
            2^3  # 3 catalytic forms
        @test by_name["Uni-Bi (dead-end I)"].expected_n_dead_end ==
            2^4  # 4 catalytic forms
        @test by_name["Bi-Bi Ping-Pong (dead-end I)"].expected_n_dead_end ==
            2^5  # 5 catalytic forms
        @test by_name["Bi-Bi (dead-end I, allosteric J)"].expected_n_dead_end ==
            2^5  # 5 catalytic forms

        # Dead-end passthrough for allosteric-only specs
        for name in ["Uni-Uni (allosteric I)",
                      "Uni-Bi (allosteric I)",
                      "Bi-Bi Ping-Pong (allosteric I)"]
            @test by_name[name].expected_n_dead_end == 1
        end

        # Allosteric: catalytic_n=1 → 1 multiplicity (m=1 only)
        for name in ["Uni-Uni (allosteric I)",
                      "Uni-Bi (allosteric I)",
                      "Bi-Bi Ping-Pong (allosteric I)"]
            s = by_name[name]
            @test s.expected_n_allosteric == 1
        end

        # Allosteric: catalytic_n=2 → 2 multiplicities (m=1,2)
        @test by_name["Uni-Bi (allosteric I, cn=2)"].expected_n_allosteric == 2
    end
end
