# ABOUTME: Tests of the test helpers' mechanism reader and exact steady-state oracle, and
# ABOUTME: of the laws over regenerated populations: agreement, detailed balance, dimension.
using Random

# The ping-pong-capable bi-bi reaction with two catalytic subunits: its seeds hold every
# bi-bi mechanism, and one expansion round adds the allosteric children.
const _testhelper_pingpong_allo_rxn = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    oligomeric_state: 2
end

"""
`m` as text that reads its content alone: each kinetic group as the sorted texts of its
steps (`_testhelper_step_text`), with its tag when `m` is allosteric, the groups sorted,
then an allosteric mechanism's catalytic multiplicity and its regulatory sites, sorted.
"""
function _testhelper_content_text(m)
    g = _testhelper_graph(m)
    allo = g.catalytic_multiplicity !== nothing
    groups = sort([join(sort([_testhelper_step_text(g, s)
                              for s in g.steps if s.group == k]), ", ") *
                   (allo ? " :: $(g.tags[k])" : "")
                   for k in unique(s.group for s in g.steps)])
    allo || return join(groups, " / ")
    sites = sort([string(site.multiplicity, ": ",
                         join(("$l::$t" for (l, t) in zip(site.ligands, site.allo_states)),
                              ", "))
                  for site in g.sites])
    join(groups, " / ") * " | n = $(g.catalytic_multiplicity) | " * join(sites, " / ")
end

# The seeds and one expansion round of the reaction above, and 120 of them spaced evenly
# through their content order (`_testhelper_content_text`). The pick reads neither the
# order the moves emit in nor a random stream, so every Julia release tests the same
# sample.
const _testhelper_population_one_round =
    _testhelper_population(_testhelper_pingpong_allo_rxn, 1)
const _testhelper_population_sample = let ms = _testhelper_population_one_round
    sorted = ms[sortperm(_testhelper_content_text.(ms))]
    sorted[round.(Int, range(1, length(sorted); length = 120))]
end

@testset "the test helpers read a mechanism's forms, steps, tags and sites" begin
    # A plain mechanism: a two-step group binding the modifier I, and the product's
    # modifier copy bound at E(S). Forms come free E first, then by content; steps and
    # groups by content, each binding stored with its metabolite taken up.
    plain = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I, P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
            (E + I ⇌ E(I), E(S) + I ⇌ E(I, S))
            E(S) + P::Mod ⇌ E(P::Mod, S)
        end
    end
    g = _testhelper_graph(plain)
    @test keys(g) == (:forms, :free, :steps, :tags, :inactive, :catalytic_multiplicity,
                      :sites, :substrates, :products, :modifiers)
    @test keys(first(g.steps)) == (:from, :to, :consumed, :released, :consumed_modifier,
                                   :released_modifier, :rapid_equilibrium, :group, :names,
                                   :inactive_names)
    @test g.forms == ["E", "EI", "EIS", "EP", "EPmodS", "ES"]
    @test g.free == [true, false, false, false, false, false]
    @test g.steps == [
        (from = 1, to = 2, consumed = [:I], released = Symbol[], consumed_modifier = [true],
         released_modifier = Bool[], rapid_equilibrium = true, group = 1,
         names = [:K_EI_to_E_I], inactive_names = Symbol[]),
        (from = 6, to = 3, consumed = [:I], released = Symbol[], consumed_modifier = [true],
         released_modifier = Bool[], rapid_equilibrium = true, group = 1,
         names = [:K_EI_to_E_I], inactive_names = Symbol[]),
        (from = 1, to = 4, consumed = [:P], released = Symbol[],
         consumed_modifier = [false], released_modifier = Bool[], rapid_equilibrium = true,
         group = 2, names = [:K_EP_to_E_P], inactive_names = Symbol[]),
        (from = 1, to = 6, consumed = [:S], released = Symbol[],
         consumed_modifier = [false], released_modifier = Bool[], rapid_equilibrium = true,
         group = 3, names = [:K_ES_to_E_S], inactive_names = Symbol[]),
        (from = 6, to = 4, consumed = Symbol[], released = Symbol[],
         consumed_modifier = Bool[], released_modifier = Bool[], rapid_equilibrium = false,
         group = 4, names = [:k_ES_to_EP, :k_EP_to_ES], inactive_names = Symbol[]),
        (from = 6, to = 5, consumed = [:P], released = Symbol[], consumed_modifier = [true],
         released_modifier = Bool[], rapid_equilibrium = true, group = 5,
         names = [:K_EPmodS_to_ES_Pmod], inactive_names = Symbol[])]
    @test (g.tags, g.inactive, g.catalytic_multiplicity, g.sites) ==
          (Symbol[], Int[], nothing, [])
    @test (g.substrates, g.products, g.modifiers) == ([:S], [:P], [:I, :P])
    # A step's text names each side's form and metabolites, as `_step_text` does.
    @test _testhelper_step_text.(Ref(g), g.steps) ==
          ["E + I ⇌ EI", "ES + I ⇌ EIS", "E + P ⇌ EP", "E + S ⇌ ES", "ES <--> EP",
           "ES + P ⇌ EPmodS"]
    @test _testhelper_graph(EnzymeRates.Mechanism(plain)) == g

    # A ping-pong: the residual forms sort by their bound metabolites, the residual-only
    # form before E(A), and each half-reaction's isomerization runs to the product side.
    pingpong = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            E(Q) ⇌ E + Q
        end
    end
    g = _testhelper_graph(pingpong)
    @test g.forms == ["E", "E_res_plusA_minusP", "EA", "EB_res_plusA_minusP",
                      "EP_res_plusA_minusP", "EQ"]
    @test g.free == [true, false, false, false, false, false]
    @test g.steps == [
        (from = 1, to = 3, consumed = [:A], released = Symbol[],
         consumed_modifier = [false], released_modifier = Bool[], rapid_equilibrium = true,
         group = 1, names = [:K_EA_to_E_A], inactive_names = Symbol[]),
        (from = 1, to = 6, consumed = [:Q], released = Symbol[],
         consumed_modifier = [false], released_modifier = Bool[], rapid_equilibrium = true,
         group = 2, names = [:K_EQ_to_E_Q], inactive_names = Symbol[]),
        (from = 2, to = 4, consumed = [:B], released = Symbol[],
         consumed_modifier = [false], released_modifier = Bool[], rapid_equilibrium = true,
         group = 3, names = [:K_EB_res_plusA_minusP_to_E_res_plusA_minusP_B],
         inactive_names = Symbol[]),
        (from = 2, to = 5, consumed = [:P], released = Symbol[],
         consumed_modifier = [false], released_modifier = Bool[], rapid_equilibrium = true,
         group = 4, names = [:K_EP_res_plusA_minusP_to_E_res_plusA_minusP_P],
         inactive_names = Symbol[]),
        (from = 3, to = 5, consumed = Symbol[], released = Symbol[],
         consumed_modifier = Bool[], released_modifier = Bool[], rapid_equilibrium = false,
         group = 5, names = [:k_EA_to_EP_res_plusA_minusP, :k_EP_res_plusA_minusP_to_EA],
         inactive_names = Symbol[]),
        (from = 4, to = 6, consumed = Symbol[], released = Symbol[],
         consumed_modifier = Bool[], released_modifier = Bool[], rapid_equilibrium = false,
         group = 6, names = [:k_EB_res_plusA_minusP_to_EQ, :k_EQ_to_EB_res_plusA_minusP],
         inactive_names = Symbol[])]
    @test (g.tags, g.inactive, g.catalytic_multiplicity, g.sites) ==
          (Symbol[], Int[], nothing, [])
    @test (g.substrates, g.products, g.modifiers) == ([:A, :B], [:P, :Q], Symbol[])
    @test _testhelper_graph(EnzymeRates.Mechanism(pingpong)) == g

    # Two conformations, so two free forms: free E first, then F among the forms of F.
    two_free = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E(S) <--> F(S), E <--> F)
            F(S) <--> E + P
            F + S <--> E(S)
        end
    end
    g = _testhelper_graph(two_free)
    @test (g.forms, g.free) == (["E", "ES", "F", "FS"], [true, false, true, false])

    # An allosteric mechanism: the inactive conformation reaches neither E(P) nor E(I, P),
    # behind the :OnlyA binding of P and the :OnlyA chemistry, so its law keeps the I
    # binding at E and the S binding only. An :EqualAI group's constant takes one bare
    # name in both conformations; an :OnlyA step's inactive copy is named after itself.
    allo = @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        allosteric_regulators: R1, R2, R3
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)                          :: EqualAI
            E(S) <--> E(P)                        :: OnlyA
            E + P ⇌ E(P)                          :: OnlyA
            (E + I ⇌ E(I), E(P) + I ⇌ E(I, P))    :: NonequalAI
        end
        regulatory_site(2): R1::OnlyI, R2::EqualAI
        regulatory_site(1): R3::NonequalAI
    end
    g = _testhelper_graph(allo)
    @test g.forms == ["E", "EI", "EIP", "EP", "ES"]
    @test g.free == [true, false, false, false, false]
    @test g.steps == [
        (from = 1, to = 2, consumed = [:I], released = Symbol[], consumed_modifier = [true],
         released_modifier = Bool[], rapid_equilibrium = true, group = 1,
         names = [:K_A_EI_to_E_I], inactive_names = [:K_I_EI_to_E_I]),
        (from = 4, to = 3, consumed = [:I], released = Symbol[], consumed_modifier = [true],
         released_modifier = Bool[], rapid_equilibrium = true, group = 1,
         names = [:K_A_EI_to_E_I], inactive_names = [:K_I_EI_to_E_I]),
        (from = 1, to = 4, consumed = [:P], released = Symbol[],
         consumed_modifier = [false], released_modifier = Bool[], rapid_equilibrium = true,
         group = 2, names = [:K_A_EP_to_E_P], inactive_names = [:K_I_EP_to_E_P]),
        (from = 1, to = 5, consumed = [:S], released = Symbol[],
         consumed_modifier = [false], released_modifier = Bool[], rapid_equilibrium = true,
         group = 3, names = [:K_ES_to_E_S], inactive_names = [:K_ES_to_E_S]),
        (from = 5, to = 4, consumed = Symbol[], released = Symbol[],
         consumed_modifier = Bool[], released_modifier = Bool[], rapid_equilibrium = false,
         group = 4, names = [:k_A_ES_to_EP, :k_A_EP_to_ES],
         inactive_names = [:k_I_ES_to_EP, :k_I_EP_to_ES])]
    @test g.tags == [:NonequalAI, :OnlyA, :EqualAI, :OnlyA]
    @test g.inactive == [1, 4]
    @test g.catalytic_multiplicity == 2
    @test g.sites ==
          [(multiplicity = 2, ligands = [:R1, :R2], allo_states = [:OnlyI, :EqualAI]),
           (multiplicity = 1, ligands = [:R3], allo_states = [:NonequalAI])]
    @test (g.substrates, g.products, g.modifiers) == ([:S], [:P], [:I])
    @test _testhelper_graph(EnzymeRates.AllostericMechanism(allo)) == g
end

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
        _, Ka_S, k1 = only(x for x in point.active if x[1].consumed == [:S])
        _, Ka_P, km2 = only(x for x in point.active if x[1].consumed == [:P])
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
            steps: begin
                E + S ⇌ E(S)      :: NonequalAI
                E(S) <--> E(P)    :: NonequalAI
                E(P) ⇌ E + P      :: EqualAI
            end
        end
        point = setprecision(() -> _testhelper_consistent_point(
            m, Random.MersenneTwister(5)), BigFloat, 256)
        flat = _testhelper_graph(m).steps
        bind_s = only(s for s in flat if s.consumed == [:S])
        bind_p = only(s for s in flat if s.consumed == [:P])
        chem = only(s for s in flat if _testhelper_is_iso(s))
        KP = point.vals[only(bind_p.names)]
        S, P = big"0.7", big"1.3"
        terms = map(((false, big(1)), (true, point.L))) do (inactive, w)
            names(s) = inactive ? s.inactive_names : s.names
            KS = point.vals[only(names(bind_s))]
            kf, kr = getindex.(Ref(point.vals), names(chem))
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
        # split is built only when the boundary check's :OnlyA test accepts it.
        n_plain, n_allo = 0, 0
        for m in _testhelper_population_sample
            names = Set(EnzymeRates._dependent_param_exprs(m)[2])
            n = EnzymeRates._independent_param_count(m)
            groups = EnzymeRates.steps(m)
            for (g, group) in enumerate(groups),
                bp in EnzymeRates._context_bipartitions(group)
                if m isa EnzymeRates.AllostericMechanism
                    parts, origin = EnzymeRates._bipartitioned_groups(groups, [(g, bp)])
                    isempty(EnzymeRates._onlya_haldane_violation(
                        m, parts, EnzymeRates.cat_allo_states(m)[origin])) || continue
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
            modifiers: I
            catalytic_multiplicity: 1
            steps: begin
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
            allosteric_regulators: R
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: NonequalAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: NonequalAI
            end
            regulatory_site(2): R::EqualAI
        end
        dep, indep =
            EnzymeRates._dependent_param_exprs(EnzymeRates.AllostericMechanism(b10))
        @test :K_Rreg in indep
        @test !any(k -> occursin("Rreg", string(k)), keys(dep))
        _testhelper_check_law_against_oracle(b10)
    end
end
