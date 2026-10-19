# ABOUTME: Tests of the test helpers' mechanism reader and of the exact, RE-speed, ε-path,
# ABOUTME: LP and n-protomer oracles, and of the laws over regenerated populations.
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
    @test keys(g) == (:forms, :free, :steps, :tags, :inactive, :active,
                      :catalytic_multiplicity, :sites, :substrates, :products, :modifiers)
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
    @test g.active == 1:6
    @test (g.substrates, g.products, g.modifiers) == ([:S], [:P], [:I, :P])
    # A step's text names each side's form and metabolites, a modifier copy by its
    # reactant's name.
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

    # A Theorell–Chance step takes up B and gives off P in one step, so its `released`
    # holds a name; it runs physical-forward, B in and P out.
    theorell_chance = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(Q) + P
            E(Q) ⇌ E + Q
        end
    end
    g = _testhelper_graph(theorell_chance)
    @test g.forms == ["E", "EA", "EQ"]
    @test g.steps == [
        (from = 1, to = 2, consumed = [:A], released = Symbol[],
         consumed_modifier = [false], released_modifier = Bool[], rapid_equilibrium = true,
         group = 1, names = [:K_EA_to_E_A], inactive_names = Symbol[]),
        (from = 1, to = 3, consumed = [:Q], released = Symbol[],
         consumed_modifier = [false], released_modifier = Bool[], rapid_equilibrium = true,
         group = 2, names = [:K_EQ_to_E_Q], inactive_names = Symbol[]),
        (from = 2, to = 3, consumed = [:B], released = [:P], consumed_modifier = [false],
         released_modifier = [false], rapid_equilibrium = false, group = 3,
         names = [:k_EA_B_to_EQ_P, :k_EQ_P_to_EA_B], inactive_names = Symbol[])]
    @test _testhelper_step_text.(Ref(g), g.steps) ==
          ["E + A ⇌ EA", "E + Q ⇌ EQ", "EA + B <--> EQ + P"]

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
    # name in both conformations; a step the inactive law does not keep has no inactive
    # names.
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
         names = [:K_A_EI_to_E_I], inactive_names = Symbol[]),
        (from = 1, to = 4, consumed = [:P], released = Symbol[],
         consumed_modifier = [false], released_modifier = Bool[], rapid_equilibrium = true,
         group = 2, names = [:K_A_EP_to_E_P], inactive_names = Symbol[]),
        (from = 1, to = 5, consumed = [:S], released = Symbol[],
         consumed_modifier = [false], released_modifier = Bool[], rapid_equilibrium = true,
         group = 3, names = [:K_ES_to_E_S], inactive_names = [:K_ES_to_E_S]),
        (from = 5, to = 4, consumed = Symbol[], released = Symbol[],
         consumed_modifier = Bool[], released_modifier = Bool[], rapid_equilibrium = false,
         group = 4, names = [:k_A_ES_to_EP, :k_A_EP_to_ES], inactive_names = Symbol[])]
    @test g.tags == [:NonequalAI, :OnlyA, :EqualAI, :OnlyA]
    @test g.inactive == [1, 4]
    @test g.active == 1:5
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
        # Two enzymes that never meet have no single steady state. The constructor
        # rejects their network; the canonical graph of the walker's raw value, which no
        # constructor check has read, holds it, and the oracle rejects it too, concrete
        # or compiled.
        @test_throws(ErrorException("Mechanism: the steps do not form one connected " *
                                    "graph; no steps lead from free E to F"),
                     @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E + P
                F + S <--> F(S)
                F(S) <--> F + P
            end
        end)
        two_enzymes = EnzymeRates._canonicalize(EnzymeRates._raw_mechanism(:(begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E + P
                F + S <--> F(S)
                F(S) <--> F + P
            end
        end), false)).mechanism
        for m in (two_enzymes, EnzymeRates.compile_mechanism(two_enzymes))
            point = _testhelper_consistent_point(m, Random.MersenneTwister(1))
            @test_throws(ErrorException("_testhelper_steady_state_rate: the network of " *
                                        "4 enzyme forms is disconnected, so its steady " *
                                        "state is not unique"),
                         _testhelper_steady_state_rate(m, point,
                                                       Dict(:S => big(1), :P => big(1))))
        end
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
                            for ns in EnzymeRates._edge_names(m, state) for n in ns)
            reg = Set(EnzymeRates.name(p, m) for state in (:A, :I)
                      for p in EnzymeRates._kreg_params(m, state))
            fitted = EnzymeRates._dependent_param_exprs(m)[2]
            append!(stray, setdiff(fitted, law_names, reg, (:L,)))
        end
        @test isempty(stray)
    end

    @testset "a same-model regrouping keeps every fitted name" begin
        # A split that the cycles tie straight back leaves the count, and so the model,
        # unchanged, and it leaves every fitted name unchanged too (E3). Only a split the
        # exit keeps counts: one the guard accepts (`_is_well_posed`).
        n_plain, n_allo = 0, 0
        for m in _testhelper_population_sample
            names = Set(EnzymeRates._dependent_param_exprs(m)[2])
            n = EnzymeRates._independent_param_count(m)
            cm = EnzymeRates.catalytic(m)
            for g in 1:EnzymeRates.ngroups(cm),
                bp in EnzymeRates._context_bipartitions(cm, g)
                child = EnzymeRates._apply_bipartitions(m, [(g, bp)])
                EnzymeRates._is_well_posed(child) || continue
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

@testset "the unguarded macro builds what the guard refuses" begin
    # R2 (spec section 7), which the macro refuses: at zero product its inactive limit
    # depends on how fast the isomerization vanishes. The unguarded macro builds it, as
    # the constructors do from the walker's raw value.
    r2 = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: EqualAI
            F + P ⇌ E(P)      :: EqualAI
            F <--> E          :: OnlyA
        end
    end
    @test r2 isa EnzymeRates.AllostericMechanism
    @test r2 == EnzymeRates.AllostericMechanism(EnzymeRates._raw_mechanism(:(begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: EqualAI
            F + P ⇌ E(P)      :: EqualAI
            F <--> E          :: OnlyA
        end
    end), true))
    @test_throws ErrorException @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: EqualAI
            F + P ⇌ E(P)      :: EqualAI
            F <--> E          :: OnlyA
        end
    end
    mm = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
        end
    end
    @test mm == EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
        end
    end)
end

@testset "the oracle's faces follow the spec's list" begin
    bibi = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(P, Q) ⇌ E(Q) + P
            E(Q) ⇌ E + Q
        end
    end
    @test _testhelper_faces(bibi) == [Set{Symbol}(), Set([:P, :Q]), Set([:A, :B]),
                                      Set([:A]), Set([:B]), Set([:P]), Set([:Q])]
    inhibited = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(P, Q) ⇌ E(Q) + P
            E(Q) ⇌ E + Q
            (E + I ⇌ E(I), E(A) + I ⇌ E(A, I))
        end
    end
    faces = _testhelper_faces(inhibited)
    @test length(faces) == 14
    @test faces[8:14] == [union(z, Set([:I])) for z in faces[1:7]]
    terter = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, C
        products: P, Q, R
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) + C ⇌ E(A, B, C)
            E(A, B, C) <--> E(P, Q, R)
            E(P, Q, R) ⇌ E(Q, R) + P
            E(Q, R) ⇌ E(R) + Q
            E(R) ⇌ E + R
        end
    end
    @test length(_testhelper_faces(terter)) == 9
    # P's modifier copy carries P's name: it is zeroed with P, and P adds no modifier face.
    copy_P = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
            E(S) + P::Mod ⇌ E(P::Mod, S)
        end
    end
    @test _testhelper_faces(copy_P) == [Set{Symbol}(), Set([:P]), Set([:S])]
    # A regulatory-site ligand is never zeroed.
    regulated = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R
        steps: begin
            E + S ⇌ E(S)      :: NonequalAI
            E(S) <--> E(P)    :: EqualAI
            E(P) ⇌ E + P      :: EqualAI
        end
        regulatory_site(1): R::OnlyA
    end
    @test _testhelper_faces(regulated) == [Set{Symbol}(), Set([:P]), Set([:S])]
end

@testset "the RE-speed oracle's error falls as 1/Λ on Michaelis–Menten" begin
    mm = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
        end
    end
    setprecision(BigFloat, 256) do
        rng = Random.MersenneTwister(11)
        point = _testhelper_consistent_point(mm, rng)
        concs = Dict(:S => big"0.7", :P => big"1.3")
        speeds = _testhelper_draw_speeds(mm, rng)
        net = _testhelper_network(mm, point, concs; speeds)
        law = _testhelper_package_law(mm, point, concs)
        # The steady state approaches the law as c/Λ, so the gap times Λ is constant.
        gap(Λ) = (_testhelper_steady_rate(net; Λ).v - law) * Λ
        @test isapprox(gap(big(10.0)^10), gap(big(10.0)^20); rtol = 1e-3)
        @test abs(_testhelper_re_limit(net) - law) <= big(10)^-60 * abs(law)
        @test _testhelper_respeed_verdicts(mm; rng) == [:agree, :agree, :agree]
    end
end

@testset "the RE-speed oracle sees a trap's pattern and an unfed closed class" begin
    # LDH-like: at zero products the forms behind E(P) + Q ⇌ E(P, Q) trap enzyme, and
    # at Q = 0 the segment of free E splits; both limits move with the RE speeds.
    ldh = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> E(P, Q)
            E(P) + Q ⇌ E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + P <--> E(P)
            E + Q ⇌ E(Q)
        end
    end
    @test _testhelper_respeed_verdicts(ldh; rng = Random.MersenneTwister(1)) ==
          [:agree, :pattern, :agree, :agree, :agree, :agree, :pattern]
    # At zero products F(P, P) is a closed set the enzyme never leaves.
    fuzz5 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E ⇌ F
            E + I ⇌ F(I)
            (E + P <--> E(P), E(P) + P <--> E(P, P))
            E + S <--> E(P)
            E + S ⇌ F(S)
            E(P, P) + S <--> F(P, P) + P
        end
    end
    @test _testhelper_respeed_verdicts(fuzz5; rng = Random.MersenneTwister(2))[2] ===
          :nonunique
end

@testset "the RE-speed oracle sees RE turnover beside a steady-state route diverge" begin
    beside = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) ⇌ E(P)
            E + P ⇌ E(P)
            E(S) <--> F(P)
            F(P) <--> E + P
        end
    end
    @test _testhelper_respeed_verdicts(beside; rng = Random.MersenneTwister(3))[1] ===
          :diverges
end

@testset "the RE-speed oracle judges an oligomer at one catalytic subunit" begin
    # The oracle's network is one subunit, so it compares with the one-subunit law, not
    # the law of the two subunits the mechanism declares.
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)      :: NonequalAI
            E(S) <--> E(P)    :: OnlyA
            E(P) ⇌ E + P      :: EqualAI
        end
        regulatory_site(2): R::NonequalAI
    end
    setprecision(BigFloat, 256) do
        rng = Random.MersenneTwister(71)
        point = _testhelper_consistent_point(am, rng)
        concs = Dict{Symbol, BigFloat}(x => big(10)^(2 * big(rand(rng)) - 1)
                                       for x in EnzymeRates.metabolites(am))
        @test _testhelper_respeed(am, point, concs, Set{Symbol}(); rng) === :agree
    end
end

@testset "the ε-path oracle rejects islands, and R2 only at P = 0" begin
    # R3 and R4 strand F behind tagged chemistry in the inactive conformation, an island:
    # its weight follows the potential the path gives it, so along the oracle's five
    # orders of the two limits an island may read direction-dependent or a wrong limit.
    r3 = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)      :: EqualAI
            E(S) <--> E(P)    :: EqualAI
            E + P ⇌ E(P)      :: EqualAI
            E <--> F          :: OnlyA
            F + I ⇌ F(I)      :: EqualAI
        end
    end
    r4 = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)      :: EqualAI
            E(S) <--> E(P)    :: NonequalAI
            E + P ⇌ E(P)      :: EqualAI
            F ⇌ E             :: OnlyA
        end
    end
    @test _testhelper_epath_verdict(r3; rng = Random.MersenneTwister(31)) in
          (:direction_dependent, :wrong_limit)
    @test _testhelper_epath_verdict(r4; rng = Random.MersenneTwister(32)) in
          (:direction_dependent, :wrong_limit)
    r2 = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: EqualAI
            F + P ⇌ E(P)      :: EqualAI
            F <--> E          :: OnlyA
        end
    end
    @test _testhelper_epath_verdict(r2; rng = Random.MersenneTwister(33)) === :converges
    @test _testhelper_epath_verdict(r2; rng = Random.MersenneTwister(34),
                                    face = Set([:P])) === :direction_dependent
    @test _testhelper_epath_verdict(r2; rng = Random.MersenneTwister(35),
                                    face = Set([:S])) === :converges
    # An exclusion that no direction can reach: no path.
    no_path = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: EqualAI
            E + P ⇌ E(P)      :: EqualAI
        end
    end
    @test _testhelper_epath_verdict(no_path; rng = Random.MersenneTwister(36)) === :no_path
end

@testset "an admitted law converges inside the ε window" begin
    # The inactive conformation's chemistry is blocked: the error falls as ε^σ.
    chem = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)      :: EqualAI
            E(S) <--> E(P)    :: OnlyA
            E(P) ⇌ E + P      :: EqualAI
        end
    end
    @test _testhelper_in_epath_window(chem; rng = Random.MersenneTwister(41))
    @test _testhelper_epath_verdict(chem; rng = Random.MersenneTwister(42)) === :converges
end

@testset "the test's certified direction satisfies every cycle row" begin
    grid = _testhelper_uniuni_grid()
    @test length(grid) == 108
    certified = 0
    for am in grid
        ex = _testhelper_exclusion(am, 2)
        @test _testhelper_certified(am, 2) == (ex.failure === :none)
        certified += ex.failure === :none
    end
    # 44 tag sets with an :OnlyA group have a direction, and the 32 without one.
    @test certified == 76
end

@testset "the exact LP agrees with the test's graph test on the uni-uni grid" begin
    grid = _testhelper_uniuni_grid()
    @test count(am -> _testhelper_lp_existence(am, 2) ==
                      (_testhelper_exclusion(am, 2).failure === :none), grid) == 108
    @test count(am -> _testhelper_lp_existence(am, 2), grid) == 76
end

@testset "the corner set holds every corner up to 128" begin
    rng = Random.MersenneTwister(51)
    seven = _testhelper_corners(7, rng)
    @test length(seven) == 128
    @test Set(seven) == Set(digits(x; base = 2, pad = 7) for x in 0:127)
    eight = _testhelper_corners(8, rng)
    @test length(eight) == 128 && length(Set(eight)) == 128
    @test all(c -> length(c) == 8 && all(in((0, 1)), c), eight)
    @test _testhelper_corners(0, rng) == [Int[]]
end

@testset "the n-protomer oracle equals the one-subunit solver at n = 1" begin
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R
        steps: begin
            E + S ⇌ E(S)      :: NonequalAI
            E(S) <--> E(P)    :: OnlyA
            E(P) ⇌ E + P      :: EqualAI
        end
        regulatory_site(2): R::NonequalAI
    end
    setprecision(BigFloat, 256) do
        rng = Random.MersenneTwister(61)
        point = _testhelper_consistent_point(am, rng)
        concs = Dict{Symbol, BigFloat}(x => big(10)^(2 * big(rand(rng)) - 1)
                                       for x in EnzymeRates.metabolites(am))
        exact = _testhelper_steady_state_rate(am, point, concs)
        @test abs(_testhelper_protomer_rate(am, point, concs; n = 1) - exact) <=
              big(10)^-20 * abs(exact)
        # Two subunits, each at one of a conformation's three law forms: nine states per
        # conformation.
        two = _testhelper_protomer_network(am, point, concs; n = 2,
                                           speeds = _testhelper_draw_speeds(am, rng))
        @test length(two.nodes) == 3^2 + 3^2
    end
end

@testset "on a plain mechanism class 2 equals the moves' RE-turnover test" begin
    # A plain mechanism's ε-network is its law network, so class 2 must flag exactly the
    # cycles the moves' test flags: on each plain spec and population member as written,
    # with one more kinetic group at rapid equilibrium, and with every step at it.
    plain = [m for m in [[_testhelper_lift(s.mechanism) for s in MECHANISM_TEST_SPECS];
                         _testhelper_population_sample] if m isa EnzymeRates.Mechanism]
    found = Bool[]
    for m in plain
        flags = EnzymeRates._flags(m)
        for re in [[flags]; [flags .| (m.group .== g) for g in 1:EnzymeRates.ngroups(m)];
                   [trues(length(m.edges))]]
            net = EnzymeRates._eps_network(m, re)
            turnover = !isempty(EnzymeRates._re_turnover_transitions(net))
            @test turnover == EnzymeRates._re_turnover_cycle(m, re)
            push!(found, turnover)
        end
    end
    @test any(found) && !all(found)
end

"""
The iso-mechanism families that test the flip at free E, as `@allosteric_mechanism` text:
Iso Uni Uni (`:uu_ss`, steady-state bindings; `:uu_re`, rapid-equilibrium bindings) and
iso ordered Bi Bi (`:bibi`), the tags in step order, and `extra` steps appended.
"""
function _testhelper_iso_text(family::Symbol, tags; n::Int = 1, extra::String = "")
    steps = family === :bibi ?
        ["E + A ⇌ E(A)", "E(A) + B <--> E(A, B)", "E(A, B) <--> Eprime(P, Q)",
         "Eprime(Q) + P <--> Eprime(P, Q)", "Eprime + Q <--> Eprime(Q)", "Eprime <--> E"] :
        family === :uu_ss ?
        ["E + S <--> E(S)", "E(S) <--> Eprime(P)", "Eprime + P <--> Eprime(P)",
         "Eprime <--> E"] :
        ["E + S ⇌ E(S)", "E(S) <--> Eprime(P)", "Eprime + P ⇌ Eprime(P)",
         "Eprime <--> E"]
    reactants = family === :bibi ? "substrates: A, B\nproducts: P, Q" :
                                   "substrates: S\nproducts: P"
    lines = join(["$s :: $t" for (s, t) in zip(steps, tags)], "\n")
    "@allosteric_mechanism begin\n$reactants\ncatalytic_multiplicity: $n\n" *
        "steps: begin\n$lines\n$extra\nend\nend"
end

@testset "the conformations flip at free E" begin
    E_, N_, O_ = :EqualAI, :NonequalAI, :OnlyA
    uu_tags = ((E_, E_, E_, N_), (E_, N_, E_, E_), (E_, O_, E_, E_), (N_, E_, E_, E_),
               (E_, E_, E_, O_), (E_, O_, E_, N_), (N_, N_, N_, N_), (O_, O_, E_, E_))
    bibi_tags = ((E_, E_, E_, E_, E_, N_), (E_, E_, N_, E_, E_, E_),
                 (E_, E_, O_, E_, E_, E_), (E_, N_, E_, E_, E_, E_),
                 (E_, E_, E_, E_, E_, O_), (N_, E_, E_, E_, E_, E_),
                 (E_, E_, E_, N_, E_, N_), (O_, E_, O_, E_, E_, E_),
                 (N_, N_, N_, N_, N_, N_))
    texts = [[_testhelper_iso_text(:uu_ss, t) for t in uu_tags];
             [_testhelper_iso_text(:uu_re, t) for t in uu_tags];
             [_testhelper_iso_text(:bibi, t) for t in bibi_tags];
             _testhelper_iso_text(:uu_ss, (E_, E_, E_, E_); extra = "D ⇌ E :: NonequalAI");
             _testhelper_iso_text(:uu_ss, (E_, N_, E_, E_); extra = "D ⇌ E :: NonequalAI")]
    # Each family member at one subunit equals the exact two-conformation steady state,
    # the conformations flipping at free E only. The macro builds every one: a second
    # free form needs no rapid-equilibrium link to E.
    for (i, text) in enumerate(texts)
        @testset "family member $i" begin
            _testhelper_check_law_against_oracle(eval(Meta.parse(text)); seed = i)
        end
    end
    # More subunits: the concerted oracle, in which only the fully unliganded oligomer
    # flips, at free E.
    for (i, (family, tags, n)) in enumerate(((:uu_ss, (E_, E_, E_, N_), 2),
                                             (:uu_ss, (E_, O_, E_, N_), 2),
                                             (:bibi, (E_, E_, E_, N_, E_, N_), 2),
                                             (:uu_ss, (E_, E_, E_, N_), 3)))
        m = eval(Meta.parse(_testhelper_iso_text(family, tags; n)))
        am = EnzymeRates.AllostericMechanism(m)
        fitted = EnzymeRates._dependent_param_exprs(typeof(m))[2]
        rng = Random.MersenneTwister(100 + i)
        point, concs, params = _testhelper_law_inputs(am, fitted, rng)
        exact = setprecision(() -> _testhelper_protomer_rate(am, point, concs; n),
                             BigFloat, 256)
        @test abs(_testhelper_law_value(m, params, concs) - exact) <=
              big(10)^-20 * abs(exact)
    end
    # A flip at the other free form gives another law: the law tells the flip form.
    m = eval(Meta.parse(_testhelper_iso_text(:uu_ss, (E_, E_, E_, N_))))
    am = EnzymeRates.AllostericMechanism(m)
    eprime = only(f for f in 2:length(am.catalytic.forms)
                  if EnzymeRates._is_free(am.catalytic.forms[f]))
    fitted = EnzymeRates._dependent_param_exprs(typeof(m))[2]
    rng = Random.MersenneTwister(201)
    for _ in 1:3
        point, concs, params = _testhelper_law_inputs(am, fitted, rng)
        law = _testhelper_law_value(m, params, concs)
        at_e = _testhelper_steady_state_rate(am, point, concs)
        at_eprime = _testhelper_steady_state_rate(am, point, concs; flip = eprime)
        @test abs(law - at_e) <= big(10)^-20 * abs(at_e)
        @test abs(law - at_eprime) > big(10)^-4 * abs(at_e)
    end
    # A free form D in rapid equilibrium with E (D's label sorts before E's) changes no
    # law: flipping at E with ratio L is flipping at D with ratio L·r_I/r_A, where
    # r_c = [D_c]/[E_c]. So L means [E_I]/[E_A].
    m = eval(Meta.parse(_testhelper_iso_text(:uu_ss, (E_, E_, E_, E_);
                                             extra = "D ⇌ E :: NonequalAI")))
    am = EnzymeRates.AllostericMechanism(m)
    cm = am.catalytic
    d = only(f for f in 2:length(cm.forms)
             if EnzymeRates._is_free(cm.forms[f]) &&
                cm.conformations[cm.forms[f].conformation] === :D)
    kd = only(j for (j, e) in enumerate(cm.edges) if Set((e.from, e.to)) == Set((1, d)))
    fitted = EnzymeRates._dependent_param_exprs(typeof(m))[2]
    rng = Random.MersenneTwister(202)
    for _ in 1:3
        point, concs, params = _testhelper_law_inputs(am, fitted, rng)
        ratio(c) = cm.edges[kd].from == d ? 1 / point.K[kd, c] : point.K[kd, c]
        at_d = _testhelper_steady_state_rate(am, merge(point, (L = point.L * ratio(2) /
                                                               ratio(1),)),
                                             concs; flip = d)
        @test abs(_testhelper_law_value(m, params, concs) - at_d) <=
              big(10)^-20 * abs(at_d)
    end
end

@testset "an :OnlyI law is exact, balanced and counts its declared dimension" begin
    # The thirteen mechanisms whose :OnlyI groups bind modifiers alone. Each law, read
    # unchecked from the concrete mechanism, equals the exact limit of the oracle's own
    # deletion network (each conformation without the steps it excludes and the forms
    # that strands, pruned by the oracle itself; the conformations flip at free E),
    # satisfies detailed balance, and fits as many constants as the model declares. The
    # guard rejects two of them; the identities hold for their laws as well.
    for (i, (label, am)) in enumerate(pairs(_testhelper_onlyi_fixtures()))
        @testset "$label" begin
            _testhelper_check_law_against_limit(am; seed = 400 + i)
            rng = Random.MersenneTwister(450 + i)
            @test all(<(big(10)^-30), _testhelper_detailed_balance(am, rng))
            @test length(EnzymeRates._dependent_param_exprs(am)[2]) ==
                  _testhelper_declared_dimension(am)
        end
    end
end

@testset "a 512-bit ε-path of bridge_REcycle2 converges at zero modifiers" begin
    # bridge_REcycle2: at zero modifiers the I bindings vanish with I, and every path of
    # the ε-network tends to the law. At 256 bits the paths lose the limit to rounding, so
    # the check runs at 512. Class 2 rejects the mechanism on the interior ("an RE binding
    # whose affinity vanishes still closes a fast cycle").
    m = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)                      :: EqualAI
            E(S) ⇌ F + P                      :: OnlyA
            F <--> E                          :: EqualAI
            (E + I ⇌ E(I), F + I ⇌ F(I))      :: OnlyI
            F(I) ⇌ E(I)                       :: OnlyA
        end
    end
    @test _testhelper_epath_verdict(m; rng = Random.MersenneTwister(91),
                                    face = Set([:I]), precision = 512) === :converges
end
