# ABOUTME: Tests of the mechanism reader, the exact, RE-speed, ε-path, LP and n-protomer
# ABOUTME: oracles, and the laws and the guard's verdicts over regenerated populations.
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
then an allosteric mechanism's oligomeric state and its regulatory sites, sorted.
"""
function _testhelper_content_text(m)
    g = _testhelper_graph(m)
    allo = g.oligomeric_state !== nothing
    groups = sort([join(sort([_testhelper_step_text(g, s)
                              for s in g.steps if s.group == k]), ", ") *
                   (allo ? " :: $(g.tags[k])" : "")
                   for k in unique(s.group for s in g.steps)])
    allo || return join(groups, " / ")
    sites = sort([string(site.multiplicity, ": ",
                         join(("$l::$t" for (l, t) in zip(site.ligands, site.allosteric_states)),
                              ", "))
                  for site in g.sites])
    join(groups, " / ") * " | n = $(g.oligomeric_state) | " * join(sites, " / ")
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

"""Whether `m`, an enumerated mechanism, holds a covalent form: one the enumerator names
with a residual, `_res` in its name."""
_testhelper_covalent(m) = any(f -> occursin("_res", String(f)), m.catalytic.enzyme_forms)

@testset "the test helpers read a mechanism's forms, steps, tags and sites" begin
    # A plain mechanism: a two-step group binding the modifier I, and the product bound
    # at ES as a modifier. Forms come free E first, then by name; steps and groups in
    # canonical order, each binding stored with its metabolite taken up.
    plain = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
            (E + I ⇌ EI, ES + I ⇌ EIS)
            ES + P ⇌ EPmodS
        end
    end
    g = _testhelper_graph(plain)
    @test keys(g) == (:forms, :steps, :tags, :inactive, :active,
                      :oligomeric_state, :sites, :substrates, :products, :modifiers)
    @test keys(first(g.steps)) == (:from, :to, :consumed, :released, :rapid_equilibrium,
                                   :group, :names, :inactive_names)
    @test g.forms == ["E", "EI", "EIS", "EP", "EPmodS", "ES"]
    @test g.steps == [
        (from = 1, to = 2, consumed = [:I], released = Symbol[], rapid_equilibrium = true,
         group = 1, names = [:K_EI_to_E_I], inactive_names = Symbol[]),
        (from = 6, to = 3, consumed = [:I], released = Symbol[], rapid_equilibrium = true,
         group = 1, names = [:K_EI_to_E_I], inactive_names = Symbol[]),
        (from = 1, to = 4, consumed = [:P], released = Symbol[], rapid_equilibrium = true,
         group = 2, names = [:K_EP_to_E_P], inactive_names = Symbol[]),
        (from = 1, to = 6, consumed = [:S], released = Symbol[], rapid_equilibrium = true,
         group = 3, names = [:K_ES_to_E_S], inactive_names = Symbol[]),
        (from = 6, to = 4, consumed = Symbol[], released = Symbol[],
         rapid_equilibrium = false, group = 4, names = [:k_ES_to_EP, :k_EP_to_ES],
         inactive_names = Symbol[]),
        (from = 6, to = 5, consumed = [:P], released = Symbol[], rapid_equilibrium = true,
         group = 5, names = [:K_EPmodS_to_ES_P], inactive_names = Symbol[])]
    @test (g.tags, g.inactive, g.oligomeric_state, g.sites) ==
          (Symbol[], Int[], nothing, [])
    @test g.active == 1:6
    @test (g.substrates, g.products, g.modifiers) == ([:S], [:P], [:I])
    # A step's text names each side's form and metabolites.
    @test _testhelper_step_text.(Ref(g), g.steps) ==
          ["E + I ⇌ EI", "ES + I ⇌ EIS", "E + P ⇌ EP", "E + S ⇌ ES", "ES <--> EP",
           "ES + P ⇌ EPmodS"]
    @test _testhelper_graph(plain) == g

    # A ping-pong: the forms by name, free E first, and each half-reaction's
    # isomerization runs to the product side.
    pingpong = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP <--> EQ
            EQ ⇌ E + Q
        end
    end
    g = _testhelper_graph(pingpong)
    @test g.forms == ["E", "EA", "EB_res_plusA_minusP", "EP_res_plusA_minusP", "EQ",
                      "E_res_plusA_minusP"]
    @test g.steps == [
        (from = 1, to = 2, consumed = [:A], released = Symbol[], rapid_equilibrium = true,
         group = 1, names = [:K_EA_to_E_A], inactive_names = Symbol[]),
        (from = 1, to = 5, consumed = [:Q], released = Symbol[], rapid_equilibrium = true,
         group = 2, names = [:K_EQ_to_E_Q], inactive_names = Symbol[]),
        (from = 2, to = 4, consumed = Symbol[], released = Symbol[],
         rapid_equilibrium = false, group = 3,
         names = [:k_EA_to_EP_res_plusA_minusP, :k_EP_res_plusA_minusP_to_EA],
         inactive_names = Symbol[]),
        (from = 3, to = 5, consumed = Symbol[], released = Symbol[],
         rapid_equilibrium = false, group = 4,
         names = [:k_EB_res_plusA_minusP_to_EQ, :k_EQ_to_EB_res_plusA_minusP],
         inactive_names = Symbol[]),
        (from = 6, to = 3, consumed = [:B], released = Symbol[], rapid_equilibrium = true,
         group = 5, names = [:K_EB_res_plusA_minusP_to_E_res_plusA_minusP_B],
         inactive_names = Symbol[]),
        (from = 6, to = 4, consumed = [:P], released = Symbol[], rapid_equilibrium = true,
         group = 6, names = [:K_EP_res_plusA_minusP_to_E_res_plusA_minusP_P],
         inactive_names = Symbol[])]
    @test (g.tags, g.inactive, g.oligomeric_state, g.sites) ==
          (Symbol[], Int[], nothing, [])
    @test (g.substrates, g.products, g.modifiers) == ([:A, :B], [:P, :Q], Symbol[])
    @test _testhelper_graph(pingpong) == g

    # A Theorell–Chance step takes up B and gives off P in one step, so its `released`
    # holds a name; it runs the way the reaction does, B in and P out.
    theorell_chance = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EQ + P
            EQ ⇌ E + Q
        end
    end
    g = _testhelper_graph(theorell_chance)
    @test g.forms == ["E", "EA", "EQ"]
    @test g.steps == [
        (from = 1, to = 2, consumed = [:A], released = Symbol[], rapid_equilibrium = true,
         group = 1, names = [:K_EA_to_E_A], inactive_names = Symbol[]),
        (from = 1, to = 3, consumed = [:Q], released = Symbol[], rapid_equilibrium = true,
         group = 2, names = [:K_EQ_to_E_Q], inactive_names = Symbol[]),
        (from = 2, to = 3, consumed = [:B], released = [:P], rapid_equilibrium = false,
         group = 3, names = [:k_EA_B_to_EQ_P, :k_EQ_P_to_EA_B], inactive_names = Symbol[])]
    @test _testhelper_step_text.(Ref(g), g.steps) ==
          ["E + A ⇌ EA", "E + Q ⇌ EQ", "EA + B <--> EQ + P"]

    # Two conformations, so two free forms: free E first, then the others by name.
    two_free = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (ES <--> FS, E <--> F)
            FS <--> E + P
            F + S <--> ES
        end
    end
    g = _testhelper_graph(two_free)
    @test g.forms == ["E", "ES", "F", "FS"]

    # An allosteric mechanism: the inactive conformation reaches neither EP nor EIP,
    # behind the :OnlyA binding of P and the :OnlyA chemistry, so its law keeps the I
    # binding at E and the S binding only. An :EqualAI group's constant takes one bare
    # name in both conformations; a step the inactive law does not keep has no inactive
    # names.
    allo = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                          :: EqualAI
            ES <--> EP                        :: OnlyA
            E + P ⇌ EP                          :: OnlyA
            (E + I ⇌ EI, EP + I ⇌ EIP)    :: NonequalAI
        end
        regulatory_site(2): R1::OnlyI, R2::EqualAI
        regulatory_site(1): R3::NonequalAI
    end
    g = _testhelper_graph(allo)
    @test g.forms == ["E", "EI", "EIP", "EP", "ES"]
    @test g.steps == [
        (from = 1, to = 2, consumed = [:I], released = Symbol[], rapid_equilibrium = true,
         group = 1, names = [:K_A_EI_to_E_I], inactive_names = [:K_I_EI_to_E_I]),
        (from = 4, to = 3, consumed = [:I], released = Symbol[], rapid_equilibrium = true,
         group = 1, names = [:K_A_EI_to_E_I], inactive_names = Symbol[]),
        (from = 1, to = 4, consumed = [:P], released = Symbol[], rapid_equilibrium = true,
         group = 2, names = [:K_A_EP_to_E_P], inactive_names = Symbol[]),
        (from = 1, to = 5, consumed = [:S], released = Symbol[], rapid_equilibrium = true,
         group = 3, names = [:K_ES_to_E_S], inactive_names = [:K_ES_to_E_S]),
        (from = 5, to = 4, consumed = Symbol[], released = Symbol[],
         rapid_equilibrium = false, group = 4, names = [:k_A_ES_to_EP, :k_A_EP_to_ES],
         inactive_names = Symbol[])]
    @test g.tags == [:NonequalAI, :OnlyA, :EqualAI, :OnlyA]
    @test g.inactive == [1, 4]
    @test g.active == 1:5
    @test g.oligomeric_state == 2
    @test g.sites ==
          [(multiplicity = 2, ligands = [:R1, :R2], allosteric_states = [:OnlyI, :EqualAI]),
           (multiplicity = 1, ligands = [:R3], allosteric_states = [:NonequalAI])]
    @test (g.substrates, g.products, g.modifiers) == ([:S], [:P], [:I])
    @test _testhelper_graph(allo) == g
end

@testset "Exact steady-state oracle" begin
    @testset "matches the steady-state Michaelis–Menten law exactly" begin
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES
                ES <--> E + P
            end
        end
        point = setprecision(() -> _testhelper_consistent_point(
            m, Random.MersenneTwister(1)), BigFloat, 256)
        # Both steps are stored as bindings: E + S → ES and E + P → ES.
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
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
            end
        end
        rng = Random.MersenneTwister(3)
        point = setprecision(() -> _testhelper_consistent_point(m, rng), BigFloat, 256)
        concs = Dict(:S => big"0.7", :P => big"1.3")
        params = Dict{Symbol, BigFloat}(k => point.vals[k]
                                        for k in parameters(m))
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
        # rejects their network; the canonical graph of the walker's raw value, built
        # without the constructor's checks, holds it, and the oracle rejects it too.
        @test_throws(ErrorException("Mechanism: the steps do not form one connected " *
                                    "graph; no steps lead from free E to F, FS"),
                     @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES
                ES <--> E + P
                F + S <--> FS
                FS <--> F + P
            end
        end)
        raw = EnzymeRates.DSL.raw_mechanism(:(begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES
                ES <--> E + P
                F + S <--> FS
                FS <--> F + P
            end
        end))
        two_enzymes = _testhelper_plain_of(EnzymeRates.Mechanisms._canonical(
            [raw.substrates; raw.products], raw.forms, [1], [2], raw.steps))
        point = _testhelper_consistent_point(two_enzymes, Random.MersenneTwister(1))
        @test_throws(ErrorException("_testhelper_steady_state_rate: the network of " *
                                    "4 enzyme forms is disconnected, so its steady " *
                                    "state is not unique"),
                     _testhelper_steady_state_rate(two_enzymes, point,
                                                   Dict(:S => big(1), :P => big(1))))
    end

    @testset "the two conformations combine at the flip ratio L" begin
        # One catalytic subunit, each conformation a rapid-equilibrium Michaelis–Menten
        # enzyme with its own S binding and chemistry: v = Σ_c w_c·(kf_c·S/KS_c −
        # kr_c·P/KP) / Σ_c w_c·(1 + S/KS_c + P/KP), with w = (1, L).
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 1
            steps: begin
                E + S ⇌ ES      :: NonequalAI
                ES <--> EP    :: NonequalAI
                EP ⇌ E + P      :: EqualAI
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
        # Haldane relation fixes k_F_to_E through square roots (Keq ^ -1//2 · …).
        half = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                (ES <--> FS, E <--> F)
                FS <--> E + P
                F + S <--> ES
            end
        end
        dep, _ = _testhelper_dependent_param_exprs(half)
        @test any(ex -> occursin("Keq ^ -1//2", string(ex)), values(dep))
        _testhelper_check_law_against_oracle(half)
        # A 1% error in Keq moves the dependent constant, and the oracle sees it.
        rng = Random.MersenneTwister(3)
        point = setprecision(() -> _testhelper_consistent_point(half, rng), BigFloat, 256)
        fitted = _testhelper_dependent_param_exprs(half)[2]
        params = Dict{Symbol, BigFloat}(k => point.vals[k] for k in fitted)
        params[:Keq] = point.Keq * big"1.01"
        concs = Dict{Symbol, BigFloat}(:S => big"0.7", :P => big"1.3")
        oracle = _testhelper_steady_state_rate(half, point, concs)
        wrong = _testhelper_law_value(half, params, concs)
        @test abs(wrong - oracle) > big(10)^-10 * abs(oracle)
    end
end

@testset "the derived laws of the population sample are exact, positive and balanced" begin
    # Spec section 14: each sampled law, derived from the concrete mechanism at one
    # catalytic subunit, equals the exact steady state at two points, holds positive
    # coefficients in each conformation's Q, vanishes at equilibrium, per conformation
    # and combined, and has one unit per denominator, which its normalizer shares when it
    # has one (`_testhelper_check_units`). Both printed forms of each law equal it at a
    # further point, with no printing artifact.
    @test count(m -> _testhelper_is_allosteric(m),
                _testhelper_population_sample) >= 30
    for (i, m) in enumerate(_testhelper_population_sample)
        @testset "population mechanism $i" begin
            mech = _testhelper_single_subunit(m)
            s = EnzymeRates.Derivation._solve(mech)
            L = EnzymeRates.Derivation._law(mech, s)
            rng = Random.MersenneTwister(i)
            for _ in 1:2
                point, concs, params = _testhelper_law_inputs(mech, s.fitted, rng)
                oracle = _testhelper_steady_state_rate(mech, point, concs)
                law = _testhelper_evaluate_law(L, merge(params, concs))
                @test abs(law - oracle) <= big(10)^-20 * abs(oracle)
            end
            Qs = L isa EnzymeRates.Derivation._Law ? [L.Q] : [L.active.Q, L.inactive.Q]
            @test all(Q -> all(>(0), values(Q.terms)), Qs)
            @test all(<(big(10)^-30),
                      _testhelper_detailed_balance(m, Random.MersenneTwister(i)))
            _testhelper_check_units(mech, L)
            # The printed law and its coefficient form of the law above, read back in
            # 256-bit BigFloat, equal it at a consistent point, with no printing artifact.
            full = EnzymeRates.Derivation._law_text(mech, s, L)
            coef = EnzymeRates.Derivation._law_text(mech, s, L; coefficients = true)
            @test isempty(_testhelper_print_artifacts(full)) &&
                  isempty(_testhelper_print_artifacts(coef))
            _, pconcs, pparams = _testhelper_law_inputs(mech, s.fitted,
                                                        Random.MersenneTwister(1000 + i))
            law = _testhelper_evaluate_law(L, merge(pparams, pconcs))
            printed(t) = _testhelper_eval_rate_string(t, pparams, pconcs, pparams[:Keq])
            @test all(t -> abs(printed(t) - law) <= big(10)^-20 * abs(law), (full, coef))
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

"""The regroupings of `m` that keep its model: each split of one group by a context
bipartition (`_context_bipartitions`) that the guard accepts and that keeps the
independent count."""
function _testhelper_same_model_regroupings(m)
    n = EnzymeRates.Derivation.fitted_parameter_count(m)
    cm = m.catalytic
    kept = EnzymeRates.Mechanism[]
    for g in eachindex(_testhelper_group_members(cm)),
        bp in EnzymeRates.Enumeration._context_bipartitions(
            _testhelper_candidate(_testhelper_plain_of(m)), g)
        child = EnzymeRates.Enumeration._apply_bipartitions(_testhelper_candidate(m),
                                                            [(g, bp)]).mechanism
        _testhelper_is_valid(child) &&
            EnzymeRates.Derivation.fitted_parameter_count(child) == n &&
            push!(kept, child)
    end
    kept
end

# The same-model regroupings of each mechanism of `_testhelper_population_sample`, in its
# order, computed once for the two testsets that read them (the guard runs on every child).
const _testhelper_sample_regroupings =
    map(_testhelper_same_model_regroupings, _testhelper_population_sample)

@testset "Laws over a regenerated population" begin
    @testset "fitted counts equal the declared dimension" begin
        # Success criterion 4: a phantom or folded fitted name would raise the count
        # above the dimension the test counts on its own.
        miscounted = [i for (i, m) in enumerate(_testhelper_population_one_round)
                      if EnzymeRates.Derivation.fitted_parameter_count(m) !=
                         _testhelper_declared_dimension(m)]
        @test isempty(miscounted)
        for m in _testhelper_population_sample[1:10]
            @test length(parameters(m)) ==
                  _testhelper_declared_dimension(m)
        end
    end

    @testset "every fitted catalytic name names a constant of the law" begin
        # Every fitted catalytic constant names a constant of one conformation's law.
        stray = Symbol[]
        for m in _testhelper_population_one_round
            _testhelper_is_allosteric(m) || continue
            law_names =
                Set(n for state in (:A, :I)
                    for ns in _testhelper_edge_names(m, state) for n in ns)
            reg = Set(K for state in (:A, :I)
                      for K in _testhelper_regulator_constants(m, state))
            fitted = _testhelper_dependent_param_exprs(m)[2]
            append!(stray, setdiff(fitted, law_names, reg, (:L,)))
        end
        @test isempty(stray)
    end

    @testset "a same-model regrouping keeps every fitted name" begin
        # A split that the cycles tie straight back leaves the count, and so the model,
        # unchanged, and it leaves every fitted name unchanged too (E3). Only a split the
        # exit keeps counts: one the guard accepts (`check_valid_rate_equation`).
        n_plain, n_allo = 0, 0
        for (m, regrouped) in zip(_testhelper_population_sample,
                                  _testhelper_sample_regroupings)
            names = Set(_testhelper_dependent_param_exprs(m)[2])
            for child in regrouped
                _testhelper_is_allosteric(m) ? (n_allo += 1) : (n_plain += 1)
                @test Set(_testhelper_dependent_param_exprs(child)[2]) == names
            end
        end
        @test n_plain >= 20
        @test n_allo >= 50
    end

    @testset "a :NonequalAI group behind an :OnlyA binding takes no inactive name" begin
        # ES and EIS lie behind the :OnlyA binding of S, so the inactive
        # conformation never reaches the I binding: its inactive class is outside the
        # law, and no K_I_ constant is fitted.
        island = @enzyme_mechanism begin
            substrates: S
            products: P
            modifiers: I
            oligomeric_state: 1
            steps: begin
                E + S <--> ES          :: OnlyA
                ES <--> EP           :: OnlyA
                E + P ⇌ EP             :: EqualAI
                ES + I ⇌ EIS       :: NonequalAI
            end
        end
        @test parameters(island) ==
              (:K_A_EIS_to_ES_I, :K_EP_to_E_P, :k_A_ES_to_EP, :k_A_ES_to_E_S,
               :k_A_E_S_to_ES, :L)
        @test length(parameters(island)) ==
              _testhelper_declared_dimension(island)
        _testhelper_check_law_against_oracle(island)
    end

    @testset "an :EqualAI regulator ligand takes one bare name and no mirror" begin
        # R binds both conformations with one constant: one bare name, K_Rreg, which both
        # conformations' site factors read, and no dependent that mirrors it (B10).
        b10 = @enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: NonequalAI
                ES <--> EP    :: EqualAI
                EP ⇌ E + P      :: NonequalAI
            end
            regulatory_site(2): R::EqualAI
        end
        dep, indep =
            _testhelper_dependent_param_exprs(b10)
        @test :K_Rreg in indep
        @test !any(k -> occursin("Rreg", string(k)), keys(dep))
        _testhelper_check_law_against_oracle(b10)
    end
end

@testset "the unguarded macro builds what the guard refuses" begin
    # S binds :OnlyA at rapid equilibrium while ES <--> EP and E + P ⇌ EP keep its
    # two sides joined in the inactive conformation: thermodynamics pins the step's
    # equilibrium constant, and a step keeps its kind, so no reading switches it off and
    # the macro refuses the mechanism. The unguarded macro builds it, as the constructors
    # do from the walker's raw value.
    pinned = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: EqualAI
            E + P ⇌ EP      :: EqualAI
        end
    end
    @test _testhelper_is_allosteric(pinned)
    @test pinned == EnzymeRates.DSL._construct(
        EnzymeRates.DSL.raw_mechanism(:(begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: EqualAI
            E + P ⇌ EP      :: EqualAI
        end
    end)))
    @test_throws ErrorException @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: EqualAI
            E + P ⇌ EP      :: EqualAI
        end
    end
    # R2 (spec section 7): under the guard's built readings, which stop F <--> E, its
    # inactive limit at zero product depends on how fast that step vanishes. The guard
    # admits it under a reading in which F <--> E loses its rate into F instead, and both
    # oracles find the law the limit on every face.
    r2 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: EqualAI
            F + P ⇌ EP      :: EqualAI
            F <--> E          :: OnlyA
        end
    end
    @test EnzymeRates.Derivation._guard(r2) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(r2; rng = Random.MersenneTwister(551),
                                              corners = true))
    _testhelper_check_law_against_limit(r2)
    mm = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
        end
    end
    @test mm == @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
        end
    end
end

@testset "the oracle's faces follow the spec's list" begin
    bibi = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EPQ ⇌ EQ + P
            EQ ⇌ E + Q
        end
    end
    @test _testhelper_faces(bibi) == [Set{Symbol}(), Set([:P, :Q]), Set([:A, :B]),
                                      Set([:A]), Set([:B]), Set([:P]), Set([:Q])]
    inhibited = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EPQ ⇌ EQ + P
            EQ ⇌ E + Q
            (E + I ⇌ EI, EA + I ⇌ EAI)
        end
    end
    faces = _testhelper_faces(inhibited)
    @test length(faces) == 14
    @test faces[8:14] == [union(z, Set([:I])) for z in faces[1:7]]
    terter = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, C
        products: P, Q, R
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB + C ⇌ EABC
            EABC <--> EPQR
            EPQR ⇌ EQR + P
            EQR ⇌ ER + Q
            ER ⇌ E + R
        end
    end
    @test length(_testhelper_faces(terter)) == 9
    # P's modifier copy carries P's name: it is zeroed with P, and P adds no modifier face.
    copy_P = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
            ES + P ⇌ EPmodS
        end
    end
    @test _testhelper_faces(copy_P) == [Set{Symbol}(), Set([:P]), Set([:S])]
    # A regulatory-site ligand is never zeroed.
    regulated = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: NonequalAI
            ES <--> EP    :: EqualAI
            EP ⇌ E + P      :: EqualAI
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
            E + S ⇌ ES
            ES <--> E + P
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
    # LDH-like: at zero products the forms behind EP + Q ⇌ EPQ trap enzyme, and
    # at Q = 0 the segment of free E splits; both limits move with the RE speeds.
    ldh = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EAB
            EAB <--> EPQ
            EP + Q ⇌ EPQ
            EQ + P ⇌ EPQ
            E + P <--> EP
            E + Q ⇌ EQ
        end
    end
    @test _testhelper_respeed_verdicts(ldh; rng = Random.MersenneTwister(1)) ==
          [:agree, :pattern, :agree, :agree, :agree, :agree, :pattern]
    # At zero products FPP is a closed set the enzyme never leaves.
    fuzz5 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E ⇌ F
            E + I ⇌ FI
            (E + P <--> EP, EP + P <--> EPP)
            E + S <--> EP
            E + S ⇌ FS
            EPP + S <--> FPP + P
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
            E + S ⇌ ES
            ES ⇌ EP
            E + P ⇌ EP
            ES <--> FP
            FP <--> E + P
        end
    end
    @test _testhelper_respeed_verdicts(beside; rng = Random.MersenneTwister(3))[1] ===
          :diverges
end

@testset "the RE-speed oracle judges an oligomer at one catalytic subunit" begin
    # The oracle's network is one subunit, so it compares with the one-subunit law, not
    # the law of the two subunits the mechanism declares.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: NonequalAI
            ES <--> EP    :: OnlyA
            EP ⇌ E + P      :: EqualAI
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

@testset "under the built readings the ε-path rejects islands, and R2 only at P = 0" begin
    # Read as the guard's built readings read it, stopped, a steady-state :OnlyA
    # isomerization strands F in the inactive conformation of R3 and of R4 with F <--> E
    # steady-state, an island: its weight follows the potential the path gives it, so along
    # the oracle's five orders of the two limits an island may read direction-dependent or
    # a wrong limit.
    r3 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: EqualAI
            E + P ⇌ EP      :: EqualAI
            E <--> F          :: OnlyA
            F + I ⇌ FI      :: EqualAI
        end
    end
    r4_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: NonequalAI
            E + P ⇌ EP      :: EqualAI
            F <--> E          :: OnlyA
        end
    end
    @test _testhelper_epath_verdict(r3; rng = Random.MersenneTwister(31)) in
          (:direction_dependent, :wrong_limit)
    @test _testhelper_epath_verdict(r4_ss; rng = Random.MersenneTwister(32)) in
          (:direction_dependent, :wrong_limit)
    # R4 itself strands nothing: a step keeps its kind, so no reading stops its
    # rapid-equilibrium F ⇌ E, which loses its rate into F and keeps the fast way back.
    # The guard admits it so, both oracles agreeing on every face.
    r4 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: NonequalAI
            E + P ⇌ EP      :: EqualAI
            F ⇌ E             :: OnlyA
        end
    end
    @test EnzymeRates.Derivation._guard(r4) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(r4; rng = Random.MersenneTwister(32),
                                              corners = true))
    _testhelper_check_law_against_limit(r4)
    # R2 read as the guard's built readings read it, F <--> E stopped: its ε-path converges
    # at the interior and at zero S, and at zero P depends on how fast F <--> E vanishes.
    # The guard admits R2 under another reading ("the unguarded macro builds what the
    # guard refuses").
    r2 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: EqualAI
            F + P ⇌ EP      :: EqualAI
            F <--> E          :: OnlyA
        end
    end
    @test _testhelper_epath_verdict(r2; rng = Random.MersenneTwister(33)) === :converges
    @test _testhelper_epath_verdict(r2; rng = Random.MersenneTwister(34),
                                    face = Set([:P])) === :direction_dependent
    @test _testhelper_epath_verdict(r2; rng = Random.MersenneTwister(35),
                                    face = Set([:S])) === :converges
    # A rapid-equilibrium :OnlyA binding whose sides the other steps keep joined, which no
    # reading switches off: no path.
    no_path = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: EqualAI
            E + P ⇌ EP      :: EqualAI
        end
    end
    @test _testhelper_epath_verdict(no_path; rng = Random.MersenneTwister(36)) === :no_path
end

@testset "an admitted law converges inside the ε window" begin
    # The inactive conformation's chemistry is blocked: the error falls as ε^σ.
    chem = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: OnlyA
            EP ⇌ E + P      :: EqualAI
        end
    end
    @test _testhelper_in_epath_window(chem; rng = Random.MersenneTwister(41))
    @test _testhelper_epath_verdict(chem; rng = Random.MersenneTwister(42)) === :converges
    # R3's island, F stranded behind tagged chemistry in the inactive conformation, takes
    # the weight the path's potential gives it: along some direction its error does not
    # fall as a power of ε.
    r3 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: EqualAI
            E + P ⇌ EP      :: EqualAI
            E <--> F          :: OnlyA
            F + I ⇌ FI      :: EqualAI
        end
    end
    @test !_testhelper_in_epath_window(r3; rng = Random.MersenneTwister(41))
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
    # 60 tag sets with an :OnlyA group have a direction, and the 32 without one. The other
    # 16 tag a rapid-equilibrium binding whose sides the other steps keep joined, a step
    # no reading switches off.
    @test certified == 92
end

@testset "the exact LP agrees with the test's graph test on the uni-uni grid" begin
    # A limit exists exactly when no excluded rapid-equilibrium step is pinned, which the
    # LP does not read, and the LP is feasible. Under the first built reading every grid
    # member's LP is feasible: a binding whose sides stay joined is stopped at steady state
    # and pinned at rapid equilibrium, and neither adds a row. So here the agreement only
    # sorts the pinned members from the rest;
    # test_allosteric_rate_equation_correctness_guard.jl's "the graph test equals
    # the exact LP on the uni-uni grid" reaches the LP's infeasible side, under every
    # reading of the inactive conformation's excluded steps.
    grid = _testhelper_uniuni_grid()
    @test count(am -> (_testhelper_exclusion(am, 2).failure !== :fast &&
                       _testhelper_lp_existence(am, 2)) ==
                      (_testhelper_exclusion(am, 2).failure === :none), grid) == 108
    @test count(am -> _testhelper_lp_existence(am, 2), grid) == 108
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
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: NonequalAI
            ES <--> EP    :: OnlyA
            EP ⇌ E + P      :: EqualAI
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
    plain = [m for m in [[s.mechanism for s in MECHANISM_TEST_SPECS];
                         _testhelper_population_sample] if !_testhelper_is_allosteric(m)]
    found = Bool[]
    for m in plain
        flags = _testhelper_flags(m)
        for re in [[flags];
                   [flags .| ([e.kinetic_group for e in m.catalytic.steps] .== g)
                    for g in eachindex(_testhelper_group_members(m.catalytic))];
                   [trues(length(m.catalytic.steps))]]
            net = EnzymeRates.Derivation._eps_network(
                m, re, EnzymeRates.Derivation._built_reading(m))
            turnover = !isempty(EnzymeRates.Derivation._re_turnover_transitions(net))
            @test turnover == EnzymeRates.Mechanisms.has_turnover_cycle(m.catalytic, re)
            push!(found, turnover)
        end
    end
    @test any(found) && !all(found)
end

"""
The iso-mechanism families that test the flip at free E, as `@enzyme_mechanism` text:
Iso Uni Uni (`:uu_ss`, steady-state bindings; `:uu_re`, rapid-equilibrium bindings) and
iso ordered Bi Bi (`:bibi`), the tags in step order, and `extra` steps appended.
"""
function _testhelper_iso_text(family::Symbol, tags; n::Int = 1, extra::String = "")
    steps = family === :bibi ?
        ["E + A ⇌ EA", "EA + B <--> EAB", "EAB <--> EprimePQ",
         "EprimeQ + P <--> EprimePQ", "Eprime + Q <--> EprimeQ", "Eprime <--> E"] :
        family === :uu_ss ?
        ["E + S <--> ES", "ES <--> EprimeP", "Eprime + P <--> EprimeP",
         "Eprime <--> E"] :
        ["E + S ⇌ ES", "ES <--> EprimeP", "Eprime + P ⇌ EprimeP",
         "Eprime <--> E"]
    reactants = family === :bibi ? "substrates: A, B\nproducts: P, Q" :
                                   "substrates: S\nproducts: P"
    lines = join(["$s :: $t" for (s, t) in zip(steps, tags)], "\n")
    "@enzyme_mechanism begin\n$reactants\noligomeric_state: $n\n" *
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
        fitted = _testhelper_dependent_param_exprs(m)[2]
        rng = Random.MersenneTwister(100 + i)
        point, concs, params = _testhelper_law_inputs(m, fitted, rng)
        exact = setprecision(() -> _testhelper_protomer_rate(m, point, concs; n),
                             BigFloat, 256)
        @test abs(_testhelper_law_value(m, params, concs) - exact) <=
              big(10)^-20 * abs(exact)
    end
    # A flip at the other free form gives another law: the law tells the flip form.
    m = eval(Meta.parse(_testhelper_iso_text(:uu_ss, (E_, E_, E_, N_))))
    eprime = findfirst(==(:Eprime), m.catalytic.enzyme_forms)
    fitted = _testhelper_dependent_param_exprs(m)[2]
    rng = Random.MersenneTwister(201)
    for _ in 1:3
        point, concs, params = _testhelper_law_inputs(m, fitted, rng)
        law = _testhelper_law_value(m, params, concs)
        at_e = _testhelper_steady_state_rate(m, point, concs)
        at_eprime = _testhelper_steady_state_rate(m, point, concs; flip = eprime)
        @test abs(law - at_e) <= big(10)^-20 * abs(at_e)
        @test abs(law - at_eprime) > big(10)^-4 * abs(at_e)
    end
    # A free form D in rapid equilibrium with E (D's label sorts before E's) changes no
    # law: flipping at E with ratio L is flipping at D with ratio L·r_I/r_A, where
    # r_c = [D_c]/[E_c]. So L means [E_I]/[E_A].
    m = eval(Meta.parse(_testhelper_iso_text(:uu_ss, (E_, E_, E_, E_);
                                             extra = "D ⇌ E :: NonequalAI")))
    cm = m.catalytic
    d = findfirst(==(:D), cm.enzyme_forms)
    kd = only(j for (j, e) in enumerate(cm.steps) if Set((e.from, e.to)) == Set((1, d)))
    fitted = _testhelper_dependent_param_exprs(m)[2]
    rng = Random.MersenneTwister(202)
    for _ in 1:3
        point, concs, params = _testhelper_law_inputs(m, fitted, rng)
        ratio(c) = cm.steps[kd].from == d ? 1 / point.K[kd, c] : point.K[kd, c]
        at_d = _testhelper_steady_state_rate(m, merge(point, (L = point.L * ratio(2) /
                                                               ratio(1),)),
                                             concs; flip = d)
        @test abs(_testhelper_law_value(m, params, concs) - at_d) <=
              big(10)^-20 * abs(at_d)
    end
end

@testset "an :OnlyI law is exact, balanced and counts its declared dimension" begin
    # The thirteen mechanisms whose :OnlyI groups bind modifiers alone, read unchecked from
    # the concrete mechanism, each fit as many constants as the model declares. The law of
    # each of twelve, all but the essential activator, whose active conformation runs no
    # step and so has none, equals the exact limit of the oracle's own deletion network
    # (each conformation without the steps it excludes and the forms that strands, pruned
    # by the oracle itself; the conformations flip at free E) and satisfies detailed
    # balance. The guard rejects the inhibitor bridge; the identities hold for its law as
    # well.
    for (i, (label, am)) in enumerate(pairs(_testhelper_onlyi_fixtures()))
        @testset "$label" begin
            if label !== :essential_activator
                _testhelper_check_law_against_limit(am; seed = 300 + i)
                rng = Random.MersenneTwister(450 + i)
                @test all(<(big(10)^-30), _testhelper_detailed_balance(am, rng))
            end
            @test length(_testhelper_dependent_param_exprs(am)[2]) ==
                  _testhelper_declared_dimension(am)
        end
    end
end

@testset "an admitted :OnlyI mechanism compiles to its law" begin
    # The eleven fixtures the guard admits (all but the inhibitor bridge and the essential
    # activator): the mechanism fits the constants its law reads, and its Float64
    # compiled law equals that law at a consistent point.
    for (i, (label, am)) in enumerate(pairs(_testhelper_onlyi_fixtures()))
        label in (:copy_group_bridge, :essential_activator) && continue
        @testset "$label" begin
            fitted = _testhelper_dependent_param_exprs(am)[2]
            @test parameters(am) == fitted
            point, concs, params = _testhelper_law_inputs(am, fitted,
                                                          Random.MersenneTwister(600 + i))
            c64, p64, keq64 = _testhelper_float64_inputs(fitted, params, concs,
                                                         EnzymeRates.metabolites(am))
            @test isapprox(compile_rate_equation(am)(c64, p64, keq64),
                           Float64(_testhelper_law_value(am, params, concs)); rtol = 1e-12)
        end
    end
end

@testset "bridge_REcycle2 has no ε-path; slowed, it converges at zero modifiers" begin
    # bridge_REcycle2: the steps the inactive conformation keeps join the two sides of
    # each rapid-equilibrium :OnlyA step, and a step keeps its kind, so no reading switches
    # them off and the ε-path oracle finds no path. Class 2 rejects the mechanism on the
    # interior ("an RE binding whose affinity vanishes still closes a fast cycle").
    m = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES                      :: EqualAI
            ES ⇌ F + P                      :: OnlyA
            F <--> E                          :: EqualAI
            (E + I ⇌ EI, F + I ⇌ FI)      :: OnlyI
            FI ⇌ EI                       :: OnlyA
        end
    end
    @test EnzymeRates.Derivation._guard(m).fault === EnzymeRates.Derivation._RE_TURNOVER
    @test _testhelper_epath_verdict(m; rng = Random.MersenneTwister(91),
                                    face = Set([:I]), precision = 512) === :no_path
    # With both :OnlyA steps steady-state the inactive conformation stops them, the guard
    # admits the mechanism, and at zero modifiers, where the I bindings vanish with I,
    # every path of the ε-network tends to the law at 512 bits.
    slow = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES                      :: EqualAI
            ES <--> F + P                   :: OnlyA
            F <--> E                          :: EqualAI
            (E + I ⇌ EI, F + I ⇌ FI)      :: OnlyI
            FI <--> EI                    :: OnlyA
        end
    end
    @test EnzymeRates.Derivation._guard(slow) == EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_epath_verdict(slow; rng = Random.MersenneTwister(91),
                                    face = Set([:I]), precision = 512) === :converges
end

# The census ping-pong declaration (census/README.md, "Reactions"): the
# ping-pong-capable bi-bi reaction with the allosteric regulator R, two catalytic
# subunits and A and P on one shared site.
const _testhelper_pingpong_census_rxn = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    allosteric_regulators: R
    oligomeric_state: 2
    shared_catalytic_site: (A, P)
end

@testset "ping-pong seeds' tagging children are admitted, distinct and oracle-checked" begin
    # An aggregate pin over the seed set: the 18 ping-pong seeds (seeds with a covalent
    # form) of the census ping-pong reaction. Each seed gets one base, which stops the
    # inactive conformation's catalysis and passes the guard (test_expansion_moves.jl,
    # "the tagging search extends a tag set whose child the guard rejects", holds a parent
    # with several bases). The move emits 94 children and the guard keeps 53; class 3
    # rejects none of the other 41, admission rejects them all. Every kept child is
    # distinct, its verdict matches the one its tags get at one subunit, on every face its
    # guard verdict agrees with the RE-speed oracle, and, under the reading the guard
    # admits it with, it converges along the ε-path at the interior, the one face the
    # suite's time budget allows it.
    rxn = _testhelper_pingpong_census_rxn
    rng = Random.MersenneTwister(406)
    seeds = [m for m in _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn))
             if _testhelper_covalent(m)]
    @test length(seeds) == 18
    nbases, n_moved, n_kids, n_class3 = Dict{Int, Int}(), 0, 0, 0
    for m in seeds
        tags = EnzymeRates.Enumeration._tagging_bases(_testhelper_candidate(m))
        nbases[length(tags)] = get(nbases, length(tags), 0) + 1
        @test all(t -> EnzymeRates.Enumeration._valid_tagging(
                           _testhelper_candidate(m), t) &&
                       EnzymeRates.Enumeration._stops_inactive_turnover(
                           _testhelper_candidate(m), t), tags)
        moved = _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, m, rxn)
        kids = filter(_testhelper_is_valid, moved)
        n_moved += length(moved)
        n_kids += length(kids)
        n_class3 += count(c -> any(!=(EnzymeRates.Derivation._NO_FAULT),
                                   _testhelper_face_faults(c)), moved)
        @test allunique(kids)
        for c in kids
            @test EnzymeRates.Enumeration._valid_tagging(
                _testhelper_candidate(_testhelper_plain_of(c)),
                c.allosteric.allosteric_states)
            faults = _testhelper_face_faults(c)
            admitted = all(==(EnzymeRates.Derivation._NO_FAULT), faults)
            verdicts = _testhelper_respeed_verdicts(c; rng,
                                                    law = admitted ? :package : :limit)
            @test all(_testhelper_agrees.(faults, verdicts))
            admitted || continue
            g = EnzymeRates.Derivation._guard_with_reading(c)
            @test g.verdict === EnzymeRates.Derivation._WELL_POSED
            reading = _testhelper_reading(g.reading)
            @test _testhelper_reading_is_consistent(c, reading)
            @test _testhelper_epath_verdict(c; rng, reading) === :converges
        end
    end
    @test nbases == Dict(1 => 18)
    @test n_moved == 94
    @test n_kids == 53
    @test n_class3 == 0
end

# ── Population tests of the well-posedness guard ────────────────────────────

# The census bi-bi declaration with A and I as competitive inhibitors (census/README.md,
# "Reactions"), whose seeds' children are the population of the two inhibitor testsets
# below (built in a let block around them).
const _testhelper_inhibitor_rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    competitive_inhibitors: A, I
end

"""`m` built by the constructors, without the guard, from its written form
(`_testhelper_written`) with `rewrite` applied to it."""
_testhelper_rebuilt(m, rewrite) = _testhelper_build(rewrite(_testhelper_written(m)))

"""`m` with every metabolite renamed by `rename`, and every form but `E` renamed so that the
forms' name order reverses (`F000` the last by name, `F001` the one before), built by the
constructors from its written form (`_testhelper_written`), without the guard. The
forms' stored order then reverses but for `E`, which stays first."""
function _testhelper_renamed(m, rename::Dict{Symbol, Symbol})
    r(xs::Vector{Symbol}) = Symbol[get(rename, x, x) for x in xs]
    others = sort(filter(!=(:E), m.catalytic.enzyme_forms))
    form = Dict(f => Symbol("F", lpad(length(others) - i, 3, '0'))
                for (i, f) in enumerate(others))
    site(s) = EnzymeRates.Mechanisms.RegulatorySite(r(s.ligands), s.multiplicity,
                                                    s.allosteric_states)
    _testhelper_rebuilt(m, w -> (; w..., substrates = r(w.substrates),
                                 products = r(w.products), modifiers = r(w.modifiers),
                                 forms = Symbol[get(form, f, f) for f in w.forms],
                                 regulatory_sites = site.(w.regulatory_sites)))
end

"""The renaming of `m`'s metabolites that reverses their name order, so the stored order
of the metabolites and the first substrate change. In sorted order the names become `Z`,
`Y`, `X`, `W`, `U`, `T`, `R`, `M`, `H` and `G`, so `m` holds at most ten metabolites (with
more, the slice throws a `BoundsError`). The letters descend, so the order reverses. None
of them is a name no metabolite may take (`E`, `L`, `N`) or a form's name (`F`), and none
is a metabolite name of the population sample or of the fixtures below (`A`, `B`, `P`,
`Q`, `S`, `ATP`, `F6P`, `ADP` and `F16BP`), so every name changes."""
function _testhelper_reversing_rename(m)
    names = sort(EnzymeRates.metabolites(m))
    Dict(zip(names, [:Z, :Y, :X, :W, :U, :T, :R, :M, :H, :G][1:length(names)]))
end

"""`m` built by the constructors from its written form with its steps in reverse order,
each written the other way round (`_testhelper_turned`), without the guard: canonical form
gives `m` back."""
_testhelper_written_backwards(m) =
    _testhelper_rebuilt(m, w -> (; w..., steps = reverse(_testhelper_turned.(w.steps))))

@testset "the guard equals the RE-speed oracle on every face of the population sample" begin
    # The guard admits each of the 120 sampled mechanisms, as the exits they passed
    # require, and on every face its class-3 verdict agrees with the RE-speed oracle at
    # one point; an allosteric member is judged at one catalytic subunit.
    rng = Random.MersenneTwister(401)
    for m in _testhelper_population_sample
        @test _testhelper_is_valid(m)
        faults = _testhelper_face_faults(m)
        @test all(_testhelper_agrees.(faults, _testhelper_respeed_verdicts(m; rng)))
    end
end

@testset "the population's allosteric members converge along the ε-path on every face" begin
    # Its excluded steps switched off as the guard's admitting reading says, a reading the
    # tests' own rule finds consistent, each allosteric member of the sample tends to its
    # deletion law on every face along the standard directions and every order of the two
    # limits, and at the interior, along the σ = 1 standard directions, its error falls as
    # a power of ε (`_testhelper_in_epath_window`).
    rng = Random.MersenneTwister(402)
    allo = [m for m in _testhelper_population_sample
            if _testhelper_is_allosteric(m)]
    @test length(allo) >= 30
    for am in allo
        g = EnzymeRates.Derivation._guard_with_reading(am)
        @test g.verdict === EnzymeRates.Derivation._WELL_POSED
        reading = _testhelper_reading(g.reading)
        @test _testhelper_reading_is_consistent(am, reading)
        @test all(face -> _testhelper_epath_verdict(am; rng, face, reading) === :converges,
                  _testhelper_faces(am))
        @test _testhelper_in_epath_window(am; rng, reading)
    end
end

@testset "the graph test equals the exact LP over the grid and the population" begin
    # Admission's existence test agrees with the exact LP, and its fault with the test's
    # own exclusion structure, on the 108 uni-uni grid mechanisms and on every allosteric
    # member of the one-round population.
    allo = [m for m in _testhelper_population_one_round
            if _testhelper_is_allosteric(m)]
    @test length(allo) >= 30
    for am in [_testhelper_uniuni_grid(); allo]
        @test _testhelper_existence_agrees(am)
    end
end

"""60 children of `children` that bind I, spaced evenly through their content order
(`_testhelper_content_text`); when fewer than 60 bind I, some repeat."""
function _testhelper_i_binding_sample(children)
    binds_i = unique!([c for c in children if :I in c.catalytic.metabolites])
    sorted = binds_i[sortperm(_testhelper_content_text.(binds_i))]
    sorted[round.(Int, range(1, length(sorted); length = 60))]
end

# The next three testsets read one population, built in a let block around them so that the
# collector frees it after them: the children the moves give the inhibitor declaration's
# seeds before candidate_mechanisms drops the ones the guard rejects.
let children = [c
                for m in _testhelper_mechanisms(
                    EnzymeRates.Enumeration.init_mechanisms(_testhelper_inhibitor_rxn))
                for c in _testhelper_plain_moves(m, _testhelper_inhibitor_rxn)]
@testset "all but 34 inhibitor children pass the guard; none tags catalysis :OnlyI" begin
    # Every child the moves give the census bi-bi seeds with inhibitors A and I passes the
    # guard but 34, which candidate_mechanisms drops, and no allosteric child tags a
    # catalytic group :OnlyI (spec section 13). Each of the 34 tags :OnlyA a fused
    # rapid-equilibrium binding into EPQ, EA + B ⇌ EPQ or EB + A ⇌ EPQ,
    # which the inactive conformation's other steps leave between two groups: at zero A
    # and B it runs one way into a cut-off form that fills faster than it empties
    # (`_FILLED`). An aggregate regression pin over the declaration's whole seed set
    # (CLAUDE.md enumeration rule 1's exception).
    @test length(children) == 5505
    rejected = filter(!_testhelper_is_valid, children)
    @test length(rejected) == 34
    @test all(c -> EnzymeRates.Derivation._guard(c).fault ===
                   EnzymeRates.Derivation._FILLED, rejected)
    allo = [c for c in children if _testhelper_is_allosteric(c)]
    @test !isempty(allo)
    @test all(c -> :OnlyI ∉ c.allosteric.allosteric_states, allo)
end

@testset "the guard equals the RE-speed oracle on inhibitor children's modifier faces" begin
    # 60 of the children above that bind I, spaced evenly through their content order; on
    # every face that zeroes I the guard's class-3 verdict agrees with the RE-speed oracle
    # at one point.
    rng = Random.MersenneTwister(405)
    sample = _testhelper_i_binding_sample(children)
    # At least 60 children bind I, so the 60 spaced picks are distinct.
    @test allunique(sample)
    for c in sample
        modifier_faces = findall(face -> :I in face, _testhelper_faces(c))
        @test !isempty(modifier_faces)
        faults = _testhelper_face_faults(c)
        verdicts = _testhelper_respeed_verdicts(c; rng)
        @test all(f -> _testhelper_agrees(faults[f], verdicts[f]), modifier_faces)
    end
end

@testset "laws that share a key are one law under the canonical map" begin
    # Spec section 14: the sample, the 60 I-binding inhibitor children, and every
    # tagging child of ten catalytic graphs of the one-round population, whose siblings
    # often share a law (tag sets that differ only on groups the inactive conformation
    # never reaches, or only in names). Members of one key evaluate alike with constants
    # mapped by canonical position, and swapping two of them breaks the map.
    by_graph = Dict{EnzymeRates.Mechanisms.CatalyticScheme,
                    Vector{EnzymeRates.Mechanism}}()
    for m in _testhelper_population_one_round
        _testhelper_is_allosteric(m) || continue
        push!(get!(by_graph, m.catalytic, EnzymeRates.Mechanism[]), m)
    end
    graphs = sort!(collect(keys(by_graph));
                   by = g -> _testhelper_content_text(_testhelper_plain_of(g)))
    siblings = reduce(vcat, [by_graph[g] for g in
                             graphs[round.(Int, range(1, length(graphs); length = 10))]])
    # A mechanism listed twice would form a merged group with itself.
    pool = unique!([_testhelper_population_sample; _testhelper_i_binding_sample(children);
                    siblings])
    keyed = Dict{String, Vector{Any}}()
    for m in pool
        L = EnzymeRates.Derivation._memo_law(m)
        c = EnzymeRates.Derivation._canonical(L, m)
        @test !c.capped
        push!(get!(keyed, c.key, Any[]), (L, c))
    end
    merged = [v for v in values(keyed) if length(v) > 1]
    @test !isempty(merged)
    rng = Random.MersenneTwister(613)
    for group in merged
        (L1, c1) = group[1]
        for (L2, c2) in group[2:end]
            @test _testhelper_same_law_under_map(L1, L2, c1.order, c2.order; rng)
        end
        length(c1.order) >= 2 || continue
        (L2, c2) = group[2]
        @test !_testhelper_same_law_under_map(L1, L2, c1.order, c2.order; rng, swap = true)
    end
end
end

@testset "the guard's verdict survives reordering, reversal, regrouping and renames" begin
    # The sample, R2, which the guard admits under a reading it finds by search, and five
    # mechanisms the guard rejects: a lumpable split rejected at Q zero only (limit test
    # `differ`), an unfed trap rejected at Q zero only (`speed`), F3 with A's binding and
    # its rapid-equilibrium Theorell–Chance step :OnlyA, whose EA fills faster than it
    # empties at zero substrates (`_FILLED`), S1's second closed class at zero substrates,
    # and a ping-pong whose inactive conformation runs neither half-reaction, pumped at
    # zero products (`_PUMPED`). Reversing the step order and every step's direction gives
    # the same mechanism. A same-model regrouping (a split the cycles tie back, the count
    # unchanged) keeps the fault. A renaming that reverses the name order moves the faces'
    # indices and the first substrate, and keeps the fault, the limit verdict, the
    # conformation and the zeroed names.
    lumpable_732 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> EA, EB + A <--> EAB, EP + A <--> EAP)
            (EA + B ⇌ EAB, EQ + B ⇌ EBQ)
            (E + P <--> EP, EA + P <--> EAP, EQ + P <--> EPQ)
            (E + Q ⇌ EQ, EB + Q ⇌ EBQ, EP + Q ⇌ EPQ)
            EAB <--> EPQ
        end
    end
    unfed_trap = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EAB
            EAB <--> EQ + P
            E + Q ⇌ EQ
            EQ + P ⇌ EPQ
            EP + Q ⇌ EPQ
            E + P <--> EP
        end
    end
    r2 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: EqualAI
            F + P ⇌ EP      :: EqualAI
            F <--> E          :: OnlyA
        end
    end
    @test EnzymeRates.Derivation._guard(r2) == EnzymeRates.Derivation._WELL_POSED
    tc = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA            :: OnlyA
            EA + B ⇌ EQ + P        :: OnlyA
            EQ <--> E + Q            :: EqualAI
        end
    end
    s1 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A <--> EA                 :: OnlyA
            EA + B <--> EAB           :: EqualAI
            EAB <--> EP               :: OnlyA
            EP <--> E + P                 :: EqualAI
            E + A ⇌ EAmod          :: OnlyA
            EAmod + B ⇌ EAmodB    :: EqualAI
            EB + A ⇌ EAmodB    :: EqualAI
            F + B <--> EB                 :: EqualAI
        end
    end
    dead = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        steps: begin
            E + ATP ⇌ EATP                                                      :: EqualAI
            EATP <--> EF16BP_res_plusATP_minusF16BP                          :: OnlyA
            E_res_plusATP_minusF16BP + F16BP ⇌ EF16BP_res_plusATP_minusF16BP:: EqualAI
            E_res_plusATP_minusF16BP + F6P ⇌ EF6P_res_plusATP_minusF16BP    :: OnlyA
            EF6P_res_plusATP_minusF16BP ⇌ EADP                               :: OnlyA
            E + ADP ⇌ EADP                                                      :: EqualAI
        end
    end
    rejected = [lumpable_732, unfed_trap, tc, s1, dead]
    verdicts = EnzymeRates.Derivation._guard.(rejected)
    @test all(v -> v.fault !== EnzymeRates.Derivation._NO_FAULT, verdicts)
    @test [(v.fault, v.limit, v.face) for v in verdicts[1:2]] ==
          [(EnzymeRates.Derivation._SPLIT, EnzymeRates.Derivation._LIMIT_DIFFER, 7),
           (EnzymeRates.Derivation._SPLIT, EnzymeRates.Derivation._LIMIT_SPEED, 7)]
    @test (verdicts[3].fault, verdicts[3].face, verdicts[3].conformation) ==
          (EnzymeRates.Derivation._FILLED, 3, 2)
    @test (verdicts[4].fault, verdicts[4].face) == (EnzymeRates.Derivation._NONUNIQUE, 3)
    @test (verdicts[5].fault, verdicts[5].face, verdicts[5].conformation) ==
          (EnzymeRates.Derivation._PUMPED, 2, 2)
    _testhelper_zeroed(m, v) = v.face == 0 ? Set{Symbol}() :
        Set(EnzymeRates.Derivation._face_names(
            m.catalytic,
            EnzymeRates.Derivation._faces(m.catalytic)[v.face]))
    for m in [_testhelper_population_sample; r2; rejected]
        v = EnzymeRates.Derivation._guard(m)
        rewritten = _testhelper_written_backwards(m)
        @test rewritten == m && EnzymeRates.Derivation._guard(rewritten) == v
        rename = _testhelper_reversing_rename(m)
        renamed = _testhelper_renamed(m, rename)
        w = EnzymeRates.Derivation._guard(renamed)
        @test (w.fault, w.limit, w.conformation) == (v.fault, v.limit, v.conformation)
        @test _testhelper_zeroed(renamed, w) ==
              Set(rename[x] for x in _testhelper_zeroed(m, v))
    end
    # The renaming moves Q's face: it is the last face of the lumpable split, and its
    # renamed image W sorts before X, P's image, so the fault moves to face 6.
    rename_732 = _testhelper_reversing_rename(lumpable_732)
    @test EnzymeRates.Derivation._guard(
        _testhelper_renamed(lumpable_732, rename_732)).face == 6
    n_regrouped = 0
    for m in _testhelper_population_sample
        fault, n = EnzymeRates.Derivation._guard(m).fault,
            EnzymeRates.Derivation.fitted_parameter_count(m)
        cm = m.catalytic
        for g in eachindex(_testhelper_group_members(cm)),
            bp in EnzymeRates.Enumeration._context_bipartitions(_testhelper_candidate(cm),
                                                                g)
            child = EnzymeRates.Enumeration._apply_bipartitions(_testhelper_candidate(m),
                                                                [(g, bp)]).mechanism
            EnzymeRates.Derivation.fitted_parameter_count(child) == n || continue
            n_regrouped += 1
            @test EnzymeRates.Derivation._guard(child).fault === fault
        end
    end
    @test n_regrouped >= 20
end

@testset "rewrites, regroupings and constant renames keep every sampled key" begin
    # Spec section 14: writing a mechanism backwards, regrouping it without changing
    # its model, and renaming its fitted constants change no law and no key, and the
    # canonical order maps one labelling onto the other. A same-model regrouping keeps
    # every fitted name too ("a same-model regrouping keeps every fitted name").
    rng = Random.MersenneTwister(612)
    n_regrouped = 0
    for (m, regrouped) in zip(_testhelper_population_sample, _testhelper_sample_regroupings)
        L = EnzymeRates.Derivation._memo_law(m)
        c = EnzymeRates.Derivation._canonical(L, m)
        @test !c.capped
        backwards = _testhelper_written_backwards(m)
        back = EnzymeRates.Derivation._memo_law(backwards)
        @test EnzymeRates.Derivation._canonical(back, backwards).key == c.key
        for child in regrouped
            n_regrouped += 1
            Lc = EnzymeRates.Derivation._memo_law(child)
            cc = EnzymeRates.Derivation._canonical(Lc, child)
            @test cc.key == c.key
            @test _testhelper_same_law_under_map(L, Lc, c.order, cc.order; rng)
        end
        constants = filter(x -> x ∉ (:Keq, :L), L.vars.names[.!L.vars.conc])
        rename = Dict(zip(constants, shuffle(rng, constants)))
        R = _testhelper_renamed_law(L, rename)
        r = EnzymeRates.Derivation._canonical(R, m)
        @test r.key == c.key
        @test Set(r.order) == Set(rename[x] for x in c.order)
        @test _testhelper_same_law_under_map(L, R, c.order, r.order; rng)
    end
    @test n_regrouped >= 70
end

# ── L phantoms over regenerated populations ─────────────────────────────────

# The census bi-bi reaction at one catalytic subunit, whose init mechanisms' tagging
# children make up the bi-bi L-phantom population tested below.
const _testhelper_phantom_bibi_rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end

# The census ping-pong declaration at one catalytic subunit, whose ping-pong parents'
# tagging children make up the ping-pong L-phantom population tested below.
const _testhelper_phantom_pingpong_rxn = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    allosteric_regulators: R
    shared_catalytic_site: (A, P)
end

"""`n` of the mechanisms `ms`, spaced evenly through their content order
(`_testhelper_content_text`), or all of them when there are at most `n`."""
function _testhelper_spaced(ms, n)
    sorted = ms[sortperm(_testhelper_content_text.(ms))]
    length(sorted) <= n ? sorted :
        sorted[round.(Int, range(1, length(sorted); length = n))]
end

"""The well-posed tagging children (`_expand_to_allosteric`) of `parents` over `rxn`,
without duplicates, in content order (`_testhelper_content_text`)."""
function _testhelper_tagging_children(parents, rxn)
    kids = unique!(
        [c for p in parents
         for c in _testhelper_move(EnzymeRates.Enumeration._expand_to_allosteric, p, rxn)
         if _testhelper_is_valid(c)])
    kids[sortperm(_testhelper_content_text.(kids))]
end

@testset "every flagged one-subunit bi-bi variant is its L→0 law's family" begin
    # Spec section 13: the tagging children of 40 bi-bi init mechanisms spaced
    # through their content order, 200 of them spaced the same way. Each variant
    # `_is_phantom` flags reproduces its law at L = 0 under the certificate's rescaling, so
    # its family is that law's and no flagged variant gains rank; the rank oracle confirms
    # it on the first four. The counts are aggregate regression pins over the
    # declaration's seed set (CLAUDE.md enumeration rule 1's exception).
    rxn = _testhelper_phantom_bibi_rxn
    parents = _testhelper_spaced(
        collect(_testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn))), 40)
    kids = _testhelper_spaced(_testhelper_tagging_children(parents, rxn), 200)
    @test length(kids) == 200
    flagged = [am for am in kids
               if EnzymeRates.Derivation._is_phantom(EnzymeRates.Derivation._memo_law(am))]
    @test length(flagged) == 76
    rng = Random.MersenneTwister(505)
    @test all(am -> _testhelper_phantom_map_holds(am; rng), flagged)
    for am in first(flagged, 4)
        @test _testhelper_identifiable_rank(am) ==
              _testhelper_identifiable_rank(_testhelper_plain_of(am))
    end
end

@testset "every flagged one-subunit ping-pong tagging child is its L→0 law's family" begin
    # The tagging children at one catalytic subunit of 30 ping-pong parents,
    # the covalent init mechanisms and their well-posed covalent flips, spaced through
    # their content order. In 93 children the inactive conformation keeps no step, and 111
    # children, kept steps or not, have a certificate. Each of these reproduces its law at
    # L = 0 under the certificate's rescaling, so its family is that law's and no flagged
    # child gains rank; the rank oracle confirms it on the first four. The counts are
    # aggregate regression pins over the declaration's seed set (CLAUDE.md enumeration
    # rule 1's exception).
    rxn = _testhelper_phantom_pingpong_rxn
    covalent = _testhelper_covalent
    inits = filter(covalent, collect(
        _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn))))
    flips = [c for m in inits
             for c in _testhelper_move(EnzymeRates.Enumeration._expand_re_to_ss, m)
             if _testhelper_is_valid(c) && covalent(c)]
    parents = _testhelper_spaced(unique!([inits; flips]), 30)
    @test length(parents) == 30
    kids = _testhelper_tagging_children(parents, rxn)
    @test length(kids) == 247
    @test count(kids) do am
        !any(EnzymeRates.Mechanisms.steps_running_in_allosteric_conformations(am).inactive)
    end == 93
    flagged = [am for am in kids
               if EnzymeRates.Derivation._is_phantom(EnzymeRates.Derivation._memo_law(am))]
    @test length(flagged) == 111
    rng = Random.MersenneTwister(505)
    @test all(am -> _testhelper_phantom_map_holds(am; rng), flagged)
    for am in first(flagged, 4)
        @test _testhelper_identifiable_rank(am) ==
              _testhelper_identifiable_rank(_testhelper_plain_of(am))
    end
end

# ── The invariant checker ─────────────────────────────────────────────────────────────

@testset "_testhelper_assert_mechanism_invariants: coverage of every reactant" begin
    # POSITIVE: the init mechanisms of a reaction with an inhibitor bind no inhibitor
    # yet; the coverage check reads substrates and products only.
    rxn_inh = @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        competitive_inhibitors: R
    end
    for m in _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn_inh))
        @test _testhelper_assert_mechanism_invariants(m) === nothing
    end

    # NEGATIVE: a declared substrate T that no step takes up. The macro rejects such a
    # mechanism; the canonical graph of the walker's raw value, which no constructor
    # check has read, holds it, and the helper must fail on it. This half checks that the
    # oracle can fail.
    @test_throws(ErrorException("Mechanism: substrate T is declared, but no step uses it"),
                 @enzyme_mechanism begin
        substrates: S, T
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end)
    raw_unused = EnzymeRates.DSL.raw_mechanism(:(begin
            substrates: S, T
            products: P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
            end
        end))
    m_unused = _testhelper_plain_of(EnzymeRates.Mechanisms._canonical(
        [raw_unused.substrates; raw_unused.products], raw_unused.forms, [1, 2], [3],
        raw_unused.steps))
    @test_throws(ErrorException("declared substrate/product T appears in no step"),
                 _testhelper_assert_mechanism_invariants(m_unused))
end
