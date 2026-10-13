# ABOUTME: Tests for the mechanism/parameter accessor API and its symmetry
# ABOUTME: across the concrete and singleton mechanism representations.
@testset "parameters API symmetry" begin
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products:   P[C], Q[N]
    end
    m_mono = EnzymeRates.EnzymeMechanism(
        first(EnzymeRates.init_mechanisms(rxn)))
    params_mono = parameters(m_mono)
    @test :E_total in params_mono
    @test :Keq in params_mono
    # all_parameters adds exactly the dependent constants to the parameters.
    values_mono = NamedTuple{params_mono}(ntuple(_ -> 1.0, length(params_mono)))
    dep_mono, _ = EnzymeRates._dependent_param_exprs(EnzymeRates.Mechanism(m_mono))
    @test !isempty(dep_mono)
    @test Set(keys(all_parameters(m_mono, values_mono))) ==
          union(Set(params_mono), keys(dep_mono))

    rxn_allo = @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        allosteric_regulators: R
        oligomeric_state: 2
    end
    base = first(EnzymeRates.init_mechanisms(rxn_allo))
    cat_allo_states = fill(:NonequalAI, length(EnzymeRates.steps(base)))
    site = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:R)], 2, [:NonequalAI])
    am = EnzymeRates.AllostericMechanism(
        EnzymeRates.reaction(base), copy(EnzymeRates.steps(base)),
        cat_allo_states, 2, [site])
    @test EnzymeRates.cat_allo_states(am) == cat_allo_states
    m_allo = EnzymeRates.AllostericEnzymeMechanism(am)
    params_allo = parameters(m_allo)
    @test any(occursin("_I_", string(p)) for p in params_allo)
    @test any(occursin("reg", string(p)) for p in params_allo)
    @test :L in params_allo
    @test :Keq in params_allo
    @test :E_total in params_allo
    values_allo = NamedTuple{params_allo}(ntuple(_ -> 1.0, length(params_allo)))
    dep_allo, _ = EnzymeRates._dependent_param_exprs(am)
    @test !isempty(dep_allo)
    @test Set(keys(all_parameters(m_allo, values_allo))) ==
          union(Set(params_allo), keys(dep_allo))
end

@testset "all_parameters adds every dependent constant" begin
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
    all_p = all_parameters(m, p)
    # The Haldane relation: k_EP_to_ES = K_EP_to_E_P · k_ES_to_EP / (K_ES_to_E_S · Keq).
    @test all_p.k_EP_to_ES == 1.0 * 3.0 / (2.0 * 4.0)
    @test keys(all_p) == (keys(p)..., :k_EP_to_ES)
    @test all_p[keys(p)] == p
    # A dependent constant already in params keeps its slot and is recomputed.
    stale = all_parameters(m, merge(p, (k_EP_to_ES = 99.0,)))
    @test keys(stale) == keys(merge(p, (k_EP_to_ES = 99.0,)))
    @test stale.k_EP_to_ES == all_p.k_EP_to_ES
    # Its slot need not be the last: a dependent constant given first stays first.
    q = merge((k_EP_to_ES = 99.0,), p)
    @test keys(all_parameters(m, q)) == keys(q)
    @test all_parameters(m, q).k_EP_to_ES == all_p.k_EP_to_ES
    # The concrete mechanism gives the same values, E_total may be left out, and BigFloat
    # values stay BigFloat.
    @test all_parameters(EnzymeRates.Mechanism(m), p) == all_p
    @test all_parameters(m, Base.structdiff(p, (E_total = 0,))).k_EP_to_ES ==
          all_p.k_EP_to_ES
    @test all_parameters(m, map(big, p)).k_EP_to_ES == big(3) / 8
    @test all_parameters(m, map(big, p)).k_EP_to_ES isa BigFloat
    @test_throws ErrorException("all_parameters: params lacks K_ES_to_E_S, Keq, which " *
                                "the dependent constants need") all_parameters(
        m, (K_EP_to_E_P = 1.0, k_ES_to_EP = 3.0))
end

@testset "all_parameters reproduces the dependent constants of a consistent point" begin
    # The point's constants come from free energies, independently of the solve; the
    # fitted ones and Keq determine the rest. `half` has a dependent with
    # half-integer exponents.
    half = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E(S) <--> F(S), E <--> F)
            F(S) <--> E + P
            F + S <--> E(S)
        end
    end
    rng = Random.MersenneTwister(12)
    plain = [s.mechanism for s in MECHANISM_TEST_SPECS
             if s.mechanism isa EnzymeRates.EnzymeMechanism]
    for m in [plain; half]
        point = setprecision(() -> _testhelper_consistent_point(m, rng), BigFloat, 256)
        fitted = EnzymeRates.fitted_params(m)
        p = NamedTuple{(fitted..., :Keq)}((getindex.(Ref(point.vals), fitted)...,
                                           point.Keq))
        all_p = setprecision(() -> all_parameters(m, p), BigFloat, 256)
        @test all(isapprox(all_p[k], point.vals[k]; rtol = big(10)^-40)
                  for k in keys(all_p) if k !== :Keq)
        @test issorted(collect(keys(all_p))[length(p) + 1:end])
    end
end

@testset "_evaluate_dependent reads the solve's power products only" begin
    p = (a = 2.0, Keq = 4.0)
    @test EnzymeRates._evaluate_dependent(:(Keq ^ -0.5 * a / a ^ 2), p) == 0.25
    @test_throws ErrorException EnzymeRates._evaluate_dependent(:(sqrt(Keq)), p)
end
