# ABOUTME: Tests for the mechanism/parameter accessor API and its symmetry
# ABOUTME: across the concrete and singleton mechanism representations.
@testset "parameters API symmetry" begin
    # A bi-bi seed whose catalytic cycle's Haldane relation makes a constant dependent;
    # B and P also bind as dead ends, at E(Q) and E(A).
    m_mono = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E + Q ⇌ E(Q)
            (E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
            (E(A) + P ⇌ E(A, P), E(Q) + P ⇌ E(P, Q))
            E(A, B) <--> E(P, Q)
        end
    end
    params_mono = parameters(m_mono)
    @test :E_total in params_mono
    @test :Keq in params_mono
    # parameters(m, values) adds exactly the dependent constants to the parameters.
    values_mono = NamedTuple{params_mono}(ntuple(_ -> 1.0, length(params_mono)))
    dep_mono, _ = _testhelper_dependent_param_exprs(EnzymeRates.Mechanism(m_mono))
    @test !isempty(dep_mono)
    @test Set(keys(parameters(m_mono, values_mono))) ==
          union(Set(params_mono), keys(dep_mono))
    @test issubset(keys(parameters(m_mono, values_mono)),
                   parameters(m_mono, AllParameters()))

    base = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
    mwc_tags = fill(:NonequalAI, length(EnzymeRates.Mechanisms.groups(base)))
    site = EnzymeRates.Mechanisms.RegulatorySite([:R], 2, [:NonequalAI])
    am = EnzymeRates.AllostericMechanism(base, mwc_tags, 2, [site])
    @test am.mwc_tags == mwc_tags
    m_allo = EnzymeRates.AllostericEnzymeMechanism(am)
    params_allo = parameters(m_allo)
    @test any(occursin("_I_", string(p)) for p in params_allo)
    @test any(occursin("reg", string(p)) for p in params_allo)
    @test :L in params_allo
    @test :Keq in params_allo
    @test :E_total in params_allo
    values_allo = NamedTuple{params_allo}(ntuple(_ -> 1.0, length(params_allo)))
    dep_allo, _ = _testhelper_dependent_param_exprs(am)
    @test !isempty(dep_allo)
    @test Set(keys(parameters(m_allo, values_allo))) ==
          union(Set(params_allo), keys(dep_allo))
    @test issubset(keys(parameters(m_allo, values_allo)),
                   parameters(m_allo, AllParameters()))
end

@testset "fitted parameters of a concrete mechanism equal its compiled form's" begin
    for spec in MECHANISM_TEST_SPECS
        @test parameters(_testhelper_lift(spec.mechanism), FittedParameters()) ===
              parameters(spec.mechanism, FittedParameters())
    end
end

@testset "parameters(m, values) adds every dependent constant" begin
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    p = (K_EP_to_E_P = 1.0, K_ES_to_E_S = 2.0, k_ES_to_EP = 3.0, Keq = 4.0, E_total = 1.0)
    all_p = parameters(m, p)
    # The Haldane relation: k_EP_to_ES = K_EP_to_E_P · k_ES_to_EP / (K_ES_to_E_S · Keq).
    @test all_p.k_EP_to_ES == 1.0 * 3.0 / (2.0 * 4.0)
    @test keys(all_p) == (keys(p)..., :k_EP_to_ES)
    @test all_p[keys(p)] == p
    # A dependent constant already in params keeps its slot and is recomputed.
    stale = parameters(m, merge(p, (k_EP_to_ES = 99.0,)))
    @test keys(stale) == keys(merge(p, (k_EP_to_ES = 99.0,)))
    @test stale.k_EP_to_ES == all_p.k_EP_to_ES
    # Its slot need not be the last: a dependent constant given first stays first.
    q = merge((k_EP_to_ES = 99.0,), p)
    @test keys(parameters(m, q)) == keys(q)
    @test parameters(m, q).k_EP_to_ES == all_p.k_EP_to_ES
    # The concrete mechanism gives the same values, E_total may be left out, and BigFloat
    # values stay BigFloat.
    @test parameters(EnzymeRates.Mechanism(m), p) == all_p
    @test parameters(m, Base.structdiff(p, (E_total = 0,))).k_EP_to_ES ==
          all_p.k_EP_to_ES
    @test parameters(m, map(big, p)).k_EP_to_ES == big(3) / 8
    @test parameters(m, map(big, p)).k_EP_to_ES isa BigFloat
    @test_throws ErrorException("parameters: params must hold every name of " *
                                "parameters(m) except E_total; it lacks " *
                                "K_ES_to_E_S, Keq") parameters(
        m, (K_EP_to_E_P = 1.0, k_ES_to_EP = 3.0))
end

@testset "parameters(m, values) reproduces the dependents of a consistent point" begin
    # The point's constants come from free energies, independently of the solve; the
    # fitted ones and Keq determine the rest. `half` has a dependent with
    # half-integer exponents. `cube` has one with third-integer exponents (k_F_to_E
    # through Keq^-1/3), which a Float64 exponent would miss by about 1e-16, far outside
    # the 1e-40 tolerance.
    half = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E(S) <--> F(S), E <--> F)
            F(S) <--> E + P
            F + S <--> E(S)
        end
    end
    cube = @enzyme_mechanism begin
        substrates: S, B
        products: P
        steps: begin
            (E <--> F, E(S) <--> F(S), E(P) <--> F(P))
            F + S <--> E(S)
            F(S) + B <--> E(P)
            F(P) <--> E + P
        end
    end
    rng = Random.MersenneTwister(12)
    for m in [[s.mechanism for s in MECHANISM_TEST_SPECS]; half; cube]
        point = setprecision(() -> _testhelper_consistent_point(m, rng), BigFloat, 256)
        fitted = parameters(m, FittedParameters())
        p = NamedTuple{(fitted..., :Keq)}((getindex.(Ref(point.vals), fitted)...,
                                           point.Keq))
        all_p = setprecision(() -> parameters(m, p), BigFloat, 256)
        @test all(isapprox(all_p[k], point.vals[k]; rtol = big(10)^-40)
                  for k in keys(all_p) if k !== :Keq)
        @test issorted(collect(keys(all_p))[length(p) + 1:end])
        # Every step constant the law reads, in each conformation, is among the keys; the
        # regulator constants and L are fitted, so p holds them.
        mech = _testhelper_lift(m)
        states = mech isa EnzymeRates.Mechanism ? (:None,) : (:A, :I)
        @test issubset((n for st in states
                        for ns in EnzymeRates.Derivation._edge_names(mech, st) for n in ns),
                       keys(all_p))
    end
end

@testset "parameters(m, values) computes rational powers of Keq from the solve" begin
    # The dependent of this conformational square is k_F_to_E = Keq^-1/2 ·
    # k_ES_to_F_S^-1/2 · k_E_P_to_FS^-1/2 · k_E_to_F · k_FS_to_E_P^1/2 · k_F_S_to_ES^1/2
    # (test_joint_solve.jl): 1/4 at Keq = k_ES_to_F_S = 4 and every other constant 1,
    # and 1 once k_F_S_to_ES is 16. BigFloat inputs give the exact value.
    half = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E(S) <--> F(S), E <--> F)
            F(S) <--> E + P
            F + S <--> E(S)
        end
    end
    p = (k_E_to_F = 1.0, k_E_P_to_FS = 1.0, k_FS_to_E_P = 1.0, k_ES_to_F_S = 4.0,
         k_F_S_to_ES = 1.0, Keq = 4.0)
    @test parameters(half, p).k_F_to_E == 0.25
    @test parameters(half, merge(p, (k_F_S_to_ES = 16.0,))).k_F_to_E == 1.0
    @test parameters(half, map(big, p)).k_F_to_E == big(0.25)
end
