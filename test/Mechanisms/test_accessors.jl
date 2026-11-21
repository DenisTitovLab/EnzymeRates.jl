# ABOUTME: Tests for the mechanism/parameter accessor API on plain and allosteric
# ABOUTME: mechanisms and on their compiled rate laws.
@testset "parameters API symmetry" begin
    # A bi-bi seed whose catalytic cycle's Haldane relation makes a constant dependent;
    # B and P also bind as dead ends, at E(Q) and E(A).
    m_mono = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            E + Q ⇌ EQ
            (EA + B ⇌ EAB, EQ + B ⇌ EBQ)
            (EA + P ⇌ EAP, EQ + P ⇌ EPQ)
            EAB <--> EPQ
        end
    end
    params_mono = parameters(m_mono)
    @test isdisjoint(params_mono, (:Keq, :E_total))
    # dependent_parameter_values(m, values, Keq) is exactly the dependent constants to the
    # parameters.
    values_mono = NamedTuple{params_mono}(ntuple(_ -> 1.0, length(params_mono)))
    dep_mono, _ = _testhelper_dependent_param_exprs(m_mono)
    @test !isempty(dep_mono)
    @test Set(keys(dependent_parameter_values(m_mono, values_mono, 1.0))) ==
          Set(keys(dep_mono))
    @test isdisjoint(keys(dependent_parameter_values(m_mono, values_mono, 1.0)),
                     params_mono)

    base = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    allosteric_states = fill(:NonequalAI, length(_testhelper_group_members(base.catalytic)))
    site = EnzymeRates.Mechanisms.RegulatorySite([:R], 2, [:NonequalAI])
    am = EnzymeRates.Mechanism(
        base.catalytic, EnzymeRates.Mechanisms.Allostery(allosteric_states, 2, [site]))
    @test am.allosteric.allosteric_states == allosteric_states
    params_allo = parameters(am)
    @test any(occursin("_I_", string(p)) for p in params_allo)
    @test any(occursin("reg", string(p)) for p in params_allo)
    @test :L in params_allo
    @test isdisjoint(params_allo, (:Keq, :E_total))
    values_allo = NamedTuple{params_allo}(ntuple(_ -> 1.0, length(params_allo)))
    dep_allo, _ = _testhelper_dependent_param_exprs(am)
    @test !isempty(dep_allo)
    @test Set(keys(dependent_parameter_values(am, values_allo, 1.0))) ==
          Set(keys(dep_allo))
    @test isdisjoint(keys(dependent_parameter_values(am, values_allo, 1.0)),
                     params_allo)
end

@testset "fitted parameters of a mechanism equal its compiled rate law's" begin
    for spec in MECHANISM_TEST_SPECS
        @test parameters(spec.mechanism) ===
              parameters(compile_rate_equation(spec.mechanism))
    end
end

@testset "dependent_parameter_values returns the dependent constants" begin
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    p = (K_EP_to_E_P = 1.0, K_ES_to_E_S = 2.0, k_ES_to_EP = 3.0)
    deps = dependent_parameter_values(m, p, 4.0)
    # The Haldane relation: k_EP_to_ES = K_EP_to_E_P · k_ES_to_EP / (K_ES_to_E_S · Keq).
    @test deps == (k_EP_to_ES = 1.0 * 3.0 / (2.0 * 4.0),)
    @test keys(deps) == (:k_EP_to_ES,)
    # A dependent constant already in params is ignored and recomputed.
    @test dependent_parameter_values(m, merge(p, (k_EP_to_ES = 99.0,)), 4.0) == deps
    @test dependent_parameter_values(m, merge((k_EP_to_ES = 99.0,), p), 4.0) == deps
    # The concrete mechanism gives the same values, the argument Keq wins over a Keq in
    # params, and BigFloat values stay BigFloat.
    @test dependent_parameter_values(m, p, 4.0) == deps
    @test dependent_parameter_values(m, merge(p, (Keq = 99.0, E_total = 99.0)), 4.0) ==
          deps
    @test dependent_parameter_values(m, map(big, p), big(4.0)).k_EP_to_ES == big(3) / 8
    @test dependent_parameter_values(m, map(big, p), big(4.0)).k_EP_to_ES isa BigFloat
    @test_throws ErrorException(
        "dependent_parameter_values: params must hold every name of " *
        "parameters(m); it lacks K_ES_to_E_S",
    ) dependent_parameter_values(m, (K_EP_to_E_P = 1.0, k_ES_to_EP = 3.0), 4.0)
end

@testset "dependent_parameter_values reproduces the dependents of a consistent point" begin
    # The point's constants come from free energies, independently of the solve; the
    # fitted ones and Keq determine the rest. `half` has a dependent with
    # half-integer exponents. `cube` has one with third-integer exponents (k_F_to_E
    # through Keq^-1/3), which a Float64 exponent would miss by about 1e-16, far outside
    # the 1e-40 tolerance.
    half = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (ES <--> FS, E <--> F)
            FS <--> E + P
            F + S <--> ES
        end
    end
    cube = @enzyme_mechanism begin
        substrates: S, B
        products: P
        steps: begin
            (E <--> F, ES <--> FS, EP <--> FP)
            F + S <--> ES
            FS + B <--> EP
            FP <--> E + P
        end
    end
    rng = Random.MersenneTwister(12)
    for m in [[s.mechanism for s in MECHANISM_TEST_SPECS]; half; cube]
        point = setprecision(() -> _testhelper_consistent_point(m, rng), BigFloat, 256)
        fitted = parameters(m)
        p = NamedTuple{fitted}(getindex.(Ref(point.vals), fitted))
        deps = setprecision(() -> dependent_parameter_values(m, p, point.Keq), BigFloat,
                            256)
        all_p = merge(p, deps)
        @test all(isapprox(all_p[k], point.vals[k]; rtol = big(10)^-40)
                  for k in keys(all_p))
        @test issorted(collect(keys(all_p))[length(p) + 1:end])
        # Every step constant the law reads, in each conformation, is among the keys; the
        # regulator constants and L are fitted, so p holds them.
        states = _testhelper_is_allosteric(m) ? (:A, :I) : (:None,)
        @test issubset((n for st in states
                        for ns in _testhelper_edge_names(m, st) for n in ns),
                       keys(all_p))
    end
end

@testset "dependent_parameter_values computes rational powers of Keq from the solve" begin
    # The dependent of this conformational square is k_F_to_E = Keq^-1/2 ·
    # k_ES_to_F_S^-1/2 · k_E_P_to_FS^-1/2 · k_E_to_F · k_FS_to_E_P^1/2 · k_F_S_to_ES^1/2
    # (test_joint_solve.jl): 1/4 at Keq = k_ES_to_F_S = 4 and every other constant 1,
    # and 1 once k_F_S_to_ES is 16. BigFloat inputs give the exact value.
    half = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (ES <--> FS, E <--> F)
            FS <--> E + P
            F + S <--> ES
        end
    end
    p = (k_E_to_F = 1.0, k_E_P_to_FS = 1.0, k_FS_to_E_P = 1.0, k_ES_to_F_S = 4.0,
         k_F_S_to_ES = 1.0)
    @test dependent_parameter_values(half, p, 4.0).k_F_to_E == 0.25
    @test dependent_parameter_values(half, merge(p, (k_F_S_to_ES = 16.0,)), 4.0).k_F_to_E ==
          1.0
    @test dependent_parameter_values(half, map(big, p), big(4.0)).k_F_to_E == big(0.25)
end
