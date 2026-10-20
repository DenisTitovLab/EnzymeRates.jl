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
    # all_parameters adds exactly the dependent constants to the parameters.
    values_mono = NamedTuple{params_mono}(ntuple(_ -> 1.0, length(params_mono)))
    dep_mono, _ = EnzymeRates._dependent_param_exprs(EnzymeRates.Mechanism(m_mono))
    @test !isempty(dep_mono)
    @test Set(keys(all_parameters(m_mono, values_mono))) ==
          union(Set(params_mono), keys(dep_mono))

    base = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
    cat_allo_states = fill(:NonequalAI, EnzymeRates.ngroups(base))
    site = EnzymeRates.RegulatorySite([:R], 2, [:NonequalAI])
    am = EnzymeRates.AllostericMechanism(base, cat_allo_states, 2, [site])
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
    @test_throws ErrorException("all_parameters: params must hold every name of " *
                                "parameters(m) except E_total; it lacks " *
                                "K_ES_to_E_S, Keq") all_parameters(
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
    for m in [[s.mechanism for s in MECHANISM_TEST_SPECS]; half]
        point = setprecision(() -> _testhelper_consistent_point(m, rng), BigFloat, 256)
        fitted = EnzymeRates.fitted_params(m)
        p = NamedTuple{(fitted..., :Keq)}((getindex.(Ref(point.vals), fitted)...,
                                           point.Keq))
        all_p = setprecision(() -> all_parameters(m, p), BigFloat, 256)
        @test all(isapprox(all_p[k], point.vals[k]; rtol = big(10)^-40)
                  for k in keys(all_p) if k !== :Keq)
        @test issorted(collect(keys(all_p))[length(p) + 1:end])
        # Every step constant the law reads, in each conformation, is among the keys; the
        # regulator constants and L are fitted, so p holds them.
        mech = _testhelper_lift(m)
        states = mech isa EnzymeRates.Mechanism ? (:None,) : (:A, :I)
        @test issubset((n for st in states for ns in EnzymeRates._edge_names(mech, st)
                        for n in ns), keys(all_p))
    end
end

@testset "all_parameters computes rational powers of Keq from the solve" begin
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
    @test all_parameters(half, p).k_F_to_E == 0.25
    @test all_parameters(half, merge(p, (k_F_S_to_ES = 16.0,))).k_F_to_E == 1.0
    @test all_parameters(half, map(big, p)).k_F_to_E == big(0.25)
end
