# ABOUTME: Tests of the exact steady-state oracle, and of the laws against it over
# ABOUTME: regenerated populations: agreement, detailed balance and declared dimension.
using Random

# The ping-pong-capable bi-bi reaction with two catalytic subunits: its seeds hold every
# bi-bi mechanism, and one expansion round adds the allosteric children.
const _testhelper_pingpong_allo_rxn = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    oligomeric_state: 2
end

# The seeds and one expansion round of the reaction above, and a sample of them that is
# the same on every run.
const _testhelper_population_one_round =
    _testhelper_population(_testhelper_pingpong_allo_rxn, 1)
const _testhelper_population_sample = _testhelper_population_one_round[randperm(
    Random.MersenneTwister(2026), length(_testhelper_population_one_round))[1:120]]

@testset "Exact steady-state oracle" begin
    @testset "matches the steady-state Michaelis–Menten law exactly" begin
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E + P
            end
        end
        point = setprecision(() -> _testhelper_consistent_point(
            m, Random.MersenneTwister(1)), BigFloat, 256)
        # Both steps are stored as bindings: E + S → E(S) and E + P → E(S).
        _, Ka_S, k1 = only(x for x in point.active
                           if EnzymeRates.consumed(x[1]) == [EnzymeRates.Substrate(:S)])
        _, Ka_P, km2 = only(x for x in point.active
                            if EnzymeRates.consumed(x[1]) == [EnzymeRates.Product(:P)])
        km1, k2 = k1 / Ka_S, km2 / Ka_P
        S, P = big"0.7", big"1.3"
        v = setprecision(BigFloat, 256) do
            (k1 * k2 * S - km1 * km2 * P) / (km1 + k2 + k1 * S + km2 * P)
        end
        oracle = _testhelper_steady_state_rate(m, point, Dict(:S => S, :P => P))
        @test abs(oracle - v) <= big(10)^-60 * abs(v)
    end

    @testset "a rapid-equilibrium step's error falls as 1/Λ" begin
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        rng = Random.MersenneTwister(3)
        point = setprecision(() -> _testhelper_consistent_point(m, rng), BigFloat, 256)
        concs = Dict(:S => big"0.7", :P => big"1.3")
        params = Dict{Symbol, BigFloat}(k => point.vals[k]
                                        for k in EnzymeRates.fitted_params(m))
        params[:Keq] = point.Keq
        law = _testhelper_law_value(m, params, concs)
        err(Λ) = abs(_testhelper_steady_state_rate(m, point, concs; Λ) - law) / abs(law)
        e10, e20, e30 = err(big(10)^10), err(big(10)^20), err(big(10)^30)
        @test e20 < 1e-9 * e10
        @test e30 < 1e-9 * e20
        @test e30 < big(10)^-20
    end

    @testset "rejects a disconnected network" begin
        two_enzymes = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E + P
                F + S <--> F(S)
                F(S) <--> F + P
            end
        end
        point = _testhelper_consistent_point(two_enzymes, Random.MersenneTwister(1))
        @test_throws ErrorException("_testhelper_steady_state_rate: the network of 4 " *
                                    "enzyme forms is disconnected, so its steady state " *
                                    "is not unique") _testhelper_steady_state_rate(
            two_enzymes, point, Dict(:S => big(1), :P => big(1)))
    end

    @testset "the two conformations combine at the flip ratio L" begin
        # One catalytic subunit, each conformation a rapid-equilibrium Michaelis–Menten
        # enzyme with its own S binding and chemistry: v = Σ_c w_c·(kf_c·S/KS_c −
        # kr_c·P/KP) / Σ_c w_c·(1 + S/KS_c + P/KP), with w = (1, L).
        m = @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + S ⇌ E(S)      :: NonequalAI
                E(S) <--> E(P)    :: NonequalAI
                E(P) ⇌ E + P      :: EqualAI
            end
        end
        am = EnzymeRates.AllostericMechanism(m)
        point = setprecision(() -> _testhelper_consistent_point(
            m, Random.MersenneTwister(5)), BigFloat, 256)
        flat = first.(EnzymeRates._flat_steps(EnzymeRates.steps(am)))
        bind_s = only(s for s in flat if EnzymeRates.consumed(s) ==
                                         [EnzymeRates.Substrate(:S)])
        bind_p = only(s for s in flat if EnzymeRates.consumed(s) ==
                                         [EnzymeRates.Product(:P)])
        chem = only(s for s in flat if EnzymeRates.is_iso(s))
        val(p) = point.vals[EnzymeRates.name(p, am)]
        KP = val(EnzymeRates.Krapid(bind_p, :EqualAI))
        S, P = big"0.7", big"1.3"
        terms = map(((:A, big(1)), (:I, point.L))) do (state, w)
            KS = val(EnzymeRates.Krapid(bind_s, state))
            kf, kr = val(EnzymeRates.Kfor(chem, state)), val(EnzymeRates.Krev(chem, state))
            (w * (kf * S / KS - kr * P / KP), w * (1 + S / KS + P / KP))
        end
        v = setprecision(() -> sum(first, terms) / sum(last, terms), BigFloat, 256)
        oracle = _testhelper_steady_state_rate(m, point, Dict(:S => S, :P => P))
        @test abs(oracle - v) <= big(10)^-25 * abs(v)
    end

    @testset "a dependent with a half-integer exponent matches the oracle" begin
        # The conformational square ties two isomerizations of one group, so the
        # Haldane relation fixes k_F_to_E through square roots (Keq ^ -0.5 · …).
        half = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                (E(S) <--> F(S), E <--> F)
                F(S) <--> E + P
                F + S <--> E(S)
            end
        end
        dep, _ = EnzymeRates._dependent_param_exprs(typeof(half))
        @test any(ex -> occursin("^ -0.5", string(ex)), values(dep))
        _testhelper_check_law_against_oracle(half)
    end
end

@testset "Detailed balance over a regenerated population" begin
    rng = Random.MersenneTwister(7)
    @test count(m -> m isa EnzymeRates.AllostericMechanism,
                _testhelper_population_sample) >= 30
    for m in _testhelper_population_sample
        @test all(<(big(10)^-30), _testhelper_detailed_balance(m, rng))
    end
end
