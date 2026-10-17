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

# The seeds and one expansion round of the reaction above, and 120 of them spaced evenly
# through that order. The pick reads no random stream, so every Julia release tests the
# same sample.
const _testhelper_population_one_round =
    _testhelper_population(_testhelper_pingpong_allo_rxn, 1)
const _testhelper_population_sample = _testhelper_population_one_round[round.(Int,
    range(1, length(_testhelper_population_one_round); length = 120))]

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
        # A 1% error in Keq moves the dependent constant, and the oracle sees it.
        rng = Random.MersenneTwister(3)
        point = setprecision(() -> _testhelper_consistent_point(half, rng), BigFloat, 256)
        fitted = EnzymeRates._dependent_param_exprs(typeof(half))[2]
        params = Dict{Symbol, BigFloat}(k => point.vals[k] for k in fitted)
        params[:Keq] = point.Keq * big"1.01"
        concs = Dict{Symbol, BigFloat}(:S => big"0.7", :P => big"1.3")
        oracle = _testhelper_steady_state_rate(half, point, concs)
        wrong = _testhelper_law_value(half, params, concs)
        @test abs(wrong - oracle) > big(10)^-10 * abs(oracle)
    end
end

@testset "Detailed balance over a regenerated population" begin
    @test count(m -> m isa EnzymeRates.AllostericMechanism,
                _testhelper_population_sample) >= 30
    for (i, m) in enumerate(_testhelper_population_sample)
        @testset "population mechanism $i" begin
            @test all(<(big(10)^-30),
                      _testhelper_detailed_balance(m, Random.MersenneTwister(i)))
        end
    end
end

# The one-round population has no inactive-conformation catalysis, so these specs check it.
@testset "Detailed balance of every spec, per conformation and combined" begin
    n_inactive_checked = 0
    for (i, spec) in enumerate(MECHANISM_TEST_SPECS)
        @testset "$(spec.name)" begin
            ratios = _testhelper_detailed_balance(spec.mechanism, Random.MersenneTwister(i))
            @test all(<(big(10)^-30), ratios)
            n_inactive_checked += length(ratios) == 3
        end
    end
    @test n_inactive_checked == 5
end

@testset "Laws over a regenerated population" begin
    @testset "every sampled law matches the exact oracle" begin
        for (i, m) in enumerate(_testhelper_population_sample)
            _testhelper_check_law_against_oracle(EnzymeRates.compile_mechanism(m);
                                                 n_points = 2, seed = i, compiled = false)
        end
    end

    @testset "fitted counts equal the declared dimension" begin
        # Success criterion 4: a phantom or folded fitted name would raise the count
        # above the dimension the test counts on its own.
        miscounted = [i for (i, m) in enumerate(_testhelper_population_one_round)
                      if EnzymeRates._independent_param_count(m) !=
                         _testhelper_declared_dimension(m)]
        @test isempty(miscounted)
        for m in _testhelper_population_sample[1:10]
            @test length(EnzymeRates.fitted_params(m)) == _testhelper_declared_dimension(m)
        end
    end

    @testset "every fitted catalytic name names a constant of the law" begin
        # Every fitted catalytic constant names a constant of one conformation's law.
        stray = Symbol[]
        for m in _testhelper_population_one_round
            m isa EnzymeRates.AllostericMechanism || continue
            law_names = Set(n for state in (:A, :I)
                            for ns in EnzymeRates._step_names(m, state) for n in ns)
            reg = Set(EnzymeRates.name(p, m) for state in (:A, :I)
                      for p in EnzymeRates._kreg_params(m, state))
            fitted = EnzymeRates._dependent_param_exprs(m)[2]
            append!(stray, setdiff(fitted, law_names, reg, (:L,)))
        end
        @test isempty(stray)
    end

    @testset "a same-model regrouping keeps every fitted name" begin
        # A split that the cycles tie straight back leaves the count, and so the model,
        # unchanged, and it leaves every fitted name unchanged too (E3). An allosteric
        # split is built only when the constructor's :OnlyA check accepts it.
        n_plain, n_allo = 0, 0
        for m in _testhelper_population_sample
            names = Set(EnzymeRates._dependent_param_exprs(m)[2])
            n = EnzymeRates._independent_param_count(m)
            groups = EnzymeRates.steps(m)
            for (g, group) in enumerate(groups),
                bp in EnzymeRates._context_bipartitions(group)
                if m isa EnzymeRates.AllostericMechanism
                    parts, origin = EnzymeRates._bipartitioned_groups(groups, [(g, bp)])
                    EnzymeRates._onlya_haldane_violation(
                        EnzymeRates.reaction(m), parts,
                        EnzymeRates.cat_allo_states(m)[origin]) === nothing || continue
                end
                child = EnzymeRates._apply_bipartitions(m, [(g, bp)])
                EnzymeRates._independent_param_count(child) == n || continue
                m isa EnzymeRates.Mechanism ? (n_plain += 1) : (n_allo += 1)
                @test Set(EnzymeRates._dependent_param_exprs(child)[2]) == names
            end
        end
        @test n_plain >= 20
        @test n_allo >= 50
    end

    @testset "a :NonequalAI group behind an :OnlyA binding takes no inactive name" begin
        # E(S) and E(I, S) lie behind the :OnlyA binding of S, so the inactive
        # conformation never reaches the I binding: its inactive class is outside the
        # law, and no K_I_ constant is fitted.
        island = @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_inhibitors: I
            catalytic_multiplicity: 1
            catalytic_steps: begin
                E + S <--> E(S)          :: OnlyA
                E(S) <--> E(P)           :: OnlyA
                E + P ⇌ E(P)             :: EqualAI
                E(S) + I ⇌ E(I, S)       :: NonequalAI
            end
        end
        @test EnzymeRates.fitted_params(island) ==
              (:K_A_EIS_to_ES_I, :K_EP_to_E_P, :k_A_ES_to_EP, :k_A_ES_to_E_S,
               :k_A_E_S_to_ES, :L)
        @test length(EnzymeRates.fitted_params(island)) ==
              _testhelper_declared_dimension(island)
        _testhelper_check_law_against_oracle(island)
    end

    @testset "an :EqualAI regulator ligand takes one bare name and no mirror" begin
        # R binds both conformations with one constant: one bare name, K_Rreg, which both
        # conformations' site factors read, and no dependent that mirrors it (B10).
        b10 = @allosteric_mechanism begin
            substrates: S
            products: P
            allosteric_regulators: R::EqualAI
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S ⇌ E(S)      :: NonequalAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: NonequalAI
            end
            regulatory_site(multiplicity = 2): begin
                ligands: R
            end
        end
        dep, indep =
            EnzymeRates._dependent_param_exprs(EnzymeRates.AllostericMechanism(b10))
        @test :K_Rreg in indep
        @test !any(k -> occursin("Rreg", string(k)), keys(dep))
        _testhelper_check_law_against_oracle(b10)
    end
end
