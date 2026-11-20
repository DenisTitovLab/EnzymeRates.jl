# ABOUTME: Tests of the well-posedness guard and admission: networks, faces, classes,
# ABOUTME: verdicts, messages, the oracles' agreement, and the macros, parser and copy rule.
using Random

@testset "the guard's faces follow the spec's list and the oracle's" begin
    for spec in MECHANISM_TEST_SPECS
        cm = EnzymeRates.Mechanisms.catalytic(_testhelper_lift(spec.mechanism))
        faces = EnzymeRates.Derivation._faces(cm)
        @test faces[1] == 0
        @test [Set(EnzymeRates.Derivation._face_names(cm, f)) for f in faces] ==
              _testhelper_faces(cm)
        # No face zeroes a substrate with a product: the guard lists no mixed corner.
        @test !any(f -> (n = EnzymeRates.Derivation._face_names(cm, f);
                         any(in(cm.substrates), n) && any(in(cm.products), n)), faces)
    end
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
    @test length.(EnzymeRates.Derivation._faces.((bibi, inhibited, terter))) == (7, 14, 9)
    # P's modifier copy shares P's name and bit: it is zeroed with P, and P adds no
    # modifier face.
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
    @test EnzymeRates.Derivation._faces(copy_P) == UInt64[0, 1, 2]
    @test EnzymeRates.Derivation._face_names(copy_P, UInt64(1)) == [:P]
end

@testset "a law the guard admits can be 0/0 at a mixed corner, scored as a 0 in a fit" begin
    # Census mechanism d7a09031d5dd2fc6 (ter-bi). At B = P = Q = 0 the rapid-equilibrium
    # segment {E(B, D), E(B, D, P), E(D, P)} has no form free of B and P, so the law is 0/0
    # there; the guard's faces never zero one substrate with every product (spec section
    # 8). At A = P = Q = 0, another such corner, the law is 0.
    am = @allosteric_mechanism begin
        substrates: A, B, D
        products: P, Q
        steps: begin
            (E + B ⇌ E(B), E(D, P) + B ⇌ E(B, D, P), E(P) + B ⇌ E(B, P)) :: OnlyA
            (E + P ⇌ E(P), E(B) + P ⇌ E(B, P), E(B, D) + P ⇌ E(B, D, P),
             E(Q) + P ⇌ E(P, Q)) :: EqualAI
            (E + Q <--> E(Q), E(P) + Q <--> E(P, Q)) :: EqualAI
            (E(A, B) + D <--> E(P, Q), E(B) + D <--> E(B, D),
             E(B, P) + D <--> E(B, D, P), E(P) + D <--> E(D, P)) :: OnlyA
            (E(B) + A <--> E(A, B), E(B, D) + A <--> E(P, Q)) :: OnlyA
            E(Q) + A ⇌ E(A, Q) :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._guard(EnzymeRates.AllostericMechanism(am)) ===
          EnzymeRates.Derivation._WELL_POSED
    params = _testhelper_random_independent_params(am; rng = Random.MersenneTwister(1))
    @test isnan(rate_equation(am, (A = 1.0, B = 0.0, D = 1.0, P = 0.0, Q = 0.0), params))
    @test rate_equation(am, (A = 0.0, B = 1.0, D = 1.0, P = 0.0, Q = 0.0), params) == 0.0
    # A fit scores the NaN as it scores the 0 a well-defined law predicts there: a zero rate
    # is refused, and a nonzero rate costs the same at the corner as at A = P = Q = 0, where
    # the law is 0.
    _testhelper_table(rate, A, B) =
        (group = [1, 1], Rate = [1.0, rate], A = [1.0, A], B = [1.0, B], D = [1.0, 1.0],
         P = [0.1, 0.0], Q = [0.1, 0.0])
    reference = _testhelper_row_reference(_testhelper_table(0.5, 1.0, 0.0))
    @test_throws(ErrorException("Zero rate at row 2: log(0) is undefined"),
                 FittingProblem(am, _testhelper_table(0.0, 1.0, 0.0); Keq = 1.0,
                                reference))
    x = [log(params[p]) for p in parameters(am, FittedParameters())]
    at_nan = FittingProblem(am, _testhelper_table(0.5, 1.0, 0.0); Keq = 1.0, reference)
    at_zero = FittingProblem(am, _testhelper_table(0.5, 0.0, 1.0); Keq = 1.0, reference)
    @test EnzymeRates.Fitting.loss!(x, at_nan) == EnzymeRates.Fitting.loss!(x, at_zero)
end

@testset "the law and ε networks: a node per form and copy, two transitions per edge" begin
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
    bind, iso = _testhelper_edge(r2, "E + S ⇌ ES"), _testhelper_edge(r2, "F <--> E")
    @test r2.excluded_inactive[bind] &&
          !r2.excluded_active[bind]
    @test r2.active_steps == trues(4)
    # The inactive law deletes both tagged steps, which strands every form but E.
    @test r2.inactive_steps == falses(4)
    law = EnzymeRates.Derivation._law_network(r2)
    @test (law.nforms, law.ncopies, length(law.from)) == (4, 2, 10)
    @test law.present == [true, true, true, true, true, false, false, false]
    @test law.edge[9:10] == [0, 0] && law.conformation[9:10] == [0, 0]
    @test (law.from[9], law.to[9]) ==
          (EnzymeRates.Derivation._node(law, 1, 1), EnzymeRates.Derivation._node(law, 1, 2))
    @test law.fast[9] && law.fast[10] && law.forward[9] && !law.forward[10]
    @test !any(law.vanishing)
    eps = EnzymeRates.Derivation._eps_network(
        r2, EnzymeRates.Derivation._re_flags(r2.catalytic),
        EnzymeRates.Derivation._built_reading(r2))
    @test (length(eps.from), eps.present) == (18, trues(8))
    _testhelper_transition(net, k, c, forward) =
        only(t for t in eachindex(net.from)
             if net.edge[t] == k && net.conformation[t] == c && net.forward[t] == forward)
    uptake, release = (_testhelper_transition(eps, bind, 2, d) for d in (true, false))
    blocked = [_testhelper_transition(eps, iso, 2, d) for d in (true, false)]
    # The tagged rapid-equilibrium binding stays fast while its affinity vanishes: its
    # uptake vanishes, its release keeps its rate.
    @test eps.fast[uptake] && eps.fast[release]
    @test eps.vanishing[uptake] && !eps.vanishing[release]
    # The tagged isomerization is blocked: slow, vanishing both ways.
    @test !any(eps.fast[blocked]) && all(eps.vanishing[blocked])
    # Made rapid-equilibrium, the tagged isomerization stays fast where its group is kept
    # and turns slow where it is excluded, since it is not a plain binding.
    eps_re = EnzymeRates.Derivation._eps_network(r2, trues(4),
                                                 EnzymeRates.Derivation._built_reading(r2))
    iso_re = [_testhelper_transition(eps_re, iso, c, d) for c in 1:2, d in (true, false)]
    @test all(eps_re.fast[iso_re[1, :]]) && !any(eps_re.fast[iso_re[2, :]])
    @test !any(eps.vanishing[t] for t in eachindex(eps.from) if eps.conformation[t] == 1)
    # Names index the table [P, S]: the binding takes up one S, its release gives it off.
    active_uptake = _testhelper_transition(eps, bind, 1, true)
    active_release = _testhelper_transition(eps, bind, 1, false)
    @test eps.change[:, active_uptake] == [0, 1]
    @test (eps.takes[active_uptake], eps.takes[active_release]) == (UInt64(2), UInt64(0))
    mm = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
        end
    end
    law_mm, eps_mm = EnzymeRates.Derivation._law_network(mm),
        EnzymeRates.Derivation._eps_network(mm, EnzymeRates.Derivation._re_flags(mm),
                                            EnzymeRates.Derivation._built_reading(mm))
    @test (law_mm.nforms, law_mm.ncopies, length(law_mm.from)) == (2, 1, 4)
    @test all(f -> getfield(law_mm, f) == getfield(eps_mm, f),
              fieldnames(EnzymeRates.Derivation._Network))
    s_bind = _testhelper_edge(mm, "E + S ⇌ ES")
    @test [law_mm.fast[t] for t in eachindex(law_mm.from) if law_mm.edge[t] == s_bind] ==
          [true, true]
    @test count(law_mm.fast) == 2
end

@testset "class 1 rejects a mechanism whose cycles never turn the reaction over" begin
    no_turnover = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) <--> F + P
        end
    end
    @test !EnzymeRates.Mechanisms.has_turnover_cycle(no_turnover, trues(3))
    v = EnzymeRates.Derivation._guard(no_turnover)
    @test v === EnzymeRates.Derivation._Verdict(
        EnzymeRates.Derivation._NO_CYCLE, EnzymeRates.Derivation._LIMIT_UNTESTED, 0, 0, 0)
    @test EnzymeRates.Derivation._guard_message(no_turnover, v) ==
          "Mechanism: no cycle of steps turns substrates into products, so the rate is " *
          "zero everywhere."
    @test_throws "holds no fault" EnzymeRates.Derivation._guard_message(
        no_turnover, EnzymeRates.Derivation._WELL_POSED)
    # The oracle's steady state carries no flux.
    setprecision(BigFloat, 256) do
        rng = Random.MersenneTwister(1)
        point = _testhelper_consistent_point(no_turnover, rng)
        net = _testhelper_network(no_turnover, point, Dict(:S => big"0.7", :P => big"1.3");
                                  speeds = _testhelper_draw_speeds(no_turnover, rng))
        @test abs(_testhelper_steady_rate(net; Λ = big(10.0)^20).v) < big(10)^-60
    end
end

@testset "class 2 rejects a rapid-equilibrium cycle that turns the reaction over" begin
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
    v = EnzymeRates.Derivation._guard(beside)
    @test (v.fault, v.conformation) == (EnzymeRates.Derivation._RE_TURNOVER, 1)
    @test EnzymeRates.Derivation._guard_message(beside, v) ==
          "Mechanism: the rapid-equilibrium steps E + P ⇌ E(P), E + S ⇌ E(S) and " *
          "E(S) ⇌ E(P) close a cycle that turns the reaction over, so the mechanism has " *
          "no finite rate. Make one step of the cycle steady-state."
    beside_allo = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)      :: EqualAI
            E(S) ⇌ E(P)       :: EqualAI
            E + P ⇌ E(P)      :: EqualAI
            E(S) <--> F(P)    :: EqualAI
            F(P) <--> E + P   :: EqualAI
        end
    end
    allo_message = EnzymeRates.Derivation._guard_message(
        beside_allo, EnzymeRates.Derivation._guard(beside_allo))
    @test startswith(allo_message,
                     "AllostericMechanism: the rapid-equilibrium steps E + P ⇌ E(P), ")
    flip2 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E + P ⇌ E(P)
            E + Q ⇌ E(Q)
            E(A) + B ⇌ E(P, Q)
            E(P) + Q ⇌ E(P, Q)
            E(Q) + P <--> E(P, Q)
        end
    end
    # test_rate_eq_derivation.jl "All-RE catalytic cycle raises": its two mechanisms.
    allre = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S1, S2
        products: P1, P2
        steps: begin
            E + S1 ⇌ E(S1)
            E + S2 ⇌ E(S2)
            E(S1) + S2 <--> E(S1, S2)
            E(S2) + S1 ⇌ E(S1, S2)
            E(S1, S2) ⇌ E(P1, P2)
            E(P1, P2) ⇌ E(P1) + P2
            E(P1, P2) <--> E(P2) + P1
            E(P1) ⇌ E + P1
            E(P2) ⇌ E + P2
        end
    end
    sib = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: R
        steps: begin
            E + S ⇌ E(S)
            E(S) ⇌ E(P)
            E + P <--> E(P)
            E(R) + S ⇌ E(S, R)
            E(S, R) ⇌ E(P, R)
            E(R) + P ⇌ E(P, R)
            (E + R ⇌ E(R), E(S) + R ⇌ E(S, R), E(P) + R ⇌ E(P, R))
        end
    end
    fix = " Make one step of the cycle steady-state."
    for (i, m) in enumerate((beside, flip2, allre, sib))
        net = EnzymeRates.Derivation._eps_network(
            m, EnzymeRates.Derivation._re_flags(EnzymeRates.Mechanisms.catalytic(m)),
            EnzymeRates.Derivation._built_reading(m))
        cycle = EnzymeRates.Derivation._re_turnover_transitions(net)
        @test !isempty(cycle) && all(t -> net.fast[t], cycle)
        v = EnzymeRates.Derivation._guard(m)
        @test v.fault === EnzymeRates.Derivation._RE_TURNOVER
        msg = EnzymeRates.Derivation._guard_message(m, v)
        @test all(
            t -> occursin(EnzymeRates.Mechanisms.step_text(m, net.edge[t]), msg), cycle)
        # The fix is printed only where each step of the cycle, made steady-state, cures
        # every fast turnover cycle (`_every_step_cures`); on sib some of its steps leave a
        # second fast cycle.
        @test endswith(msg, "so the mechanism has no finite rate." *
                            (m === sib ? "" : fix))
        # The RE-speed oracle at the interior. On sib's faces with R zero no fast cycle
        # holds enzyme, so the oracle would compare with the package law there, and a
        # mechanism class 2 rejects has no law.
        verdict = setprecision(BigFloat, 256) do
            rng = Random.MersenneTwister(i)
            point = _testhelper_consistent_point(m, rng)
            concs = Dict{Symbol, BigFloat}(x => big(10)^(2 * big(rand(rng)) - 1)
                                           for x in EnzymeRates.metabolites(m))
            _testhelper_respeed(m, point, concs, Set{Symbol}(); rng)
        end
        @test verdict === :diverges
    end
    # The fix: one step of the cycle steady-state. Each variant is well-posed.
    fixed = (
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E + P ⇌ E(P)
                E(S) <--> F(P)
                F(P) <--> E + P
            end
        end),
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E + P ⇌ E(P)
                E + Q ⇌ E(Q)
                E(A) + B <--> E(P, Q)
                E(P) + Q ⇌ E(P, Q)
                E(Q) + P <--> E(P, Q)
            end
        end),
        @_testhelper_unguarded(@enzyme_mechanism begin
            substrates: S1, S2
            products: P1, P2
            steps: begin
                E + S1 ⇌ E(S1)
                E + S2 ⇌ E(S2)
                E(S1) + S2 <--> E(S1, S2)
                E(S2) + S1 ⇌ E(S1, S2)
                E(S1, S2) <--> E(P1, P2)
                E(P1, P2) ⇌ E(P1) + P2
                E(P1, P2) <--> E(P2) + P1
                E(P1) ⇌ E + P1
                E(P2) ⇌ E + P2
            end
        end))
    @test all(
        m -> EnzymeRates.Derivation._guard(m) === EnzymeRates.Derivation._WELL_POSED, fixed)
    # sib turns over through two fast cycles, and some steps of the reported one leave
    # the other when made steady-state, so its message offers no fix: with both
    # chemistry steps steady-state no fast cycle turns over, but at R = 0 that variant
    # splits the segment of free E, which class 3 rejects.
    sib_fixed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: R
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P <--> E(P)
            E(R) + S ⇌ E(S, R)
            E(S, R) <--> E(P, R)
            E(R) + P ⇌ E(P, R)
            (E + R ⇌ E(R), E(S) + R ⇌ E(S, R), E(P) + R ⇌ E(P, R))
        end
    end
    @test isempty(EnzymeRates.Derivation._re_turnover_transitions(
        EnzymeRates.Derivation._eps_network(
            sib_fixed, EnzymeRates.Derivation._re_flags(sib_fixed),
            EnzymeRates.Derivation._built_reading(sib_fixed))))
end

@testset "the fitted classes pass classes 1 and 2" begin
    michaelis_menten = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
        end
    end
    abortive_complex = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P <--> E(P, Q)
            E + Q <--> E(Q)
            E(Q) + A ⇌ E(A, Q)
        end
    end
    product_activated = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ E(P)
            E(P) + S ⇌ E(P, S)
            E(P, S) <--> E(P, P)
            E(P) + P ⇌ E(P, P)
        end
    end
    ss_flux_free_dead_end = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
            E + I <--> E(I)
            E(S) <--> F(S)
        end
    end
    for m in (michaelis_menten, abortive_complex, product_activated, ss_flux_free_dead_end)
        @test EnzymeRates.Derivation._guard(m) === EnzymeRates.Derivation._WELL_POSED
    end
end

@testset "each step a message names reads back as that step" begin
    pingpong = @_testhelper_unguarded @enzyme_mechanism begin
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
    theorell_chance = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(Q) + P
            E(Q) ⇌ E + Q
        end
    end
    merged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(P, Q)
            E(P, Q) ⇌ E(Q) + P
            E(Q) ⇌ E + Q
        end
    end
    copies = @_testhelper_unguarded @enzyme_mechanism begin
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
    mechs = [[_testhelper_lift(s.mechanism) for s in MECHANISM_TEST_SPECS];
             [pingpong, theorell_chance, merged, copies]]
    for mech in mechs
        cm = EnzymeRates.Mechanisms.catalytic(mech)
        names = (substrates = cm.substrates, products = cm.products,
                 modifiers = cm.modifiers,
                 regulators = mech isa EnzymeRates.AllostericMechanism ?
                              EnzymeRates.Mechanisms.allosteric_regulators(mech) : Symbol[])
        for (k, e) in enumerate(cm.steps)
            op, lhs, rhs = EnzymeRates.DSL._step_parts(
                Meta.parse(EnzymeRates.Mechanisms.step_text(cm, k)), "@enzyme_mechanism")
            nodes = EnzymeRates.Mechanisms.EnzymeForm[]
            s = EnzymeRates.DSL._raw_step!(nodes, op, lhs, rhs, e.group, names,
                                           "@enzyme_mechanism")
            @test nodes == cm.nodes[[e.from, e.to]] &&
                  s == EnzymeRates.Mechanisms.Step(1, 2, e.taken_up, e.given_off,
                                                   e.rapid_equilibrium, e.group)
        end
    end
    line(m, text) = EnzymeRates.Mechanisms.step_text(m, _testhelper_edge(m, text))
    @test line(pingpong, "E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP") ==
          "E(; residual = A - P) + B ⇌ E(B; residual = A - P)"
    @test line(theorell_chance, "EA + B <--> EQ + P") == "E(A) + B <--> E(Q) + P"
    @test line(merged, "EA + B <--> EPQ") == "E(A) + B <--> E(P, Q)"
    @test line(copies, "ES + Pmod ⇌ EPmodS") == "E(S) + P::Mod ⇌ E(P::Mod, S)"
    @test line(copies, "E + I ⇌ EI") == "E + I ⇌ E(I)"
end

@testset "class 3 accepts the fitted classes and lumpable splits on every face" begin
    michaelis_menten = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
        end
    end
    abortive_complex = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P <--> E(P, Q)
            E + Q <--> E(Q)
            E(Q) + A ⇌ E(A, Q)
        end
    end
    product_activated = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ E(P)
            E(P) + S ⇌ E(P, S)
            E(P, S) <--> E(P, P)
            E(P) + P ⇌ E(P, P)
        end
    end
    ss_flux_free_dead_end = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
            E + I <--> E(I)
            E(S) <--> F(S)
        end
    end
    # R2's active conformation: a second free form, and no trap (spec section 8).
    r2_active = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            F + P ⇌ E(P)
            F <--> E
        end
    end
    tied_isoforms = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E + S <--> E(S), F + S <--> F(S))
            (E(S) <--> E + P, F(S) <--> F + P)
            (E + P ⇌ G(P), F + P ⇌ G(P))
            E <--> F
        end
    end
    untied_isoforms = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(S)
            F + S <--> F(S)
            (E(S) <--> E + P, F(S) <--> F + P)
            (E + P ⇌ G(P), F + P ⇌ G(P))
            E <--> F
        end
    end
    copy_trap = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
            E(P) + P::Mod ⇌ E(P, P::Mod)
            E(P, P::Mod) <--> E(P::Mod) + P
            E + P::Mod <--> E(P::Mod)
        end
    end
    lumpable_172 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> E(A), E(B) + A <--> E(P, Q))
            (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))
            (E + Q ⇌ E(Q), E(B) + Q ⇌ E(B, Q))
            (E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q))
        end
    end
    lumpable_532 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            (E(B) + A ⇌ E(P, Q), E(Q) + A ⇌ E(A, Q))
            (E + B <--> E(B), E(A) + B <--> E(P, Q), E(I) + B <--> E(B, I))
            (E + I ⇌ E(I), E(B) + I ⇌ E(B, I), E(Q) + I ⇌ E(I, Q))
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q), E(I) + Q ⇌ E(I, Q))
            E(Q) + P <--> E(P, Q)
        end
    end
    deadend_split = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
            E + P ⇌ E(P)
            F + P ⇌ E(P)
            F <--> E
        end
    end
    accepted = (michaelis_menten, abortive_complex, product_activated,
                ss_flux_free_dead_end, r2_active, tied_isoforms, untied_isoforms, copy_trap,
                lumpable_172, lumpable_532, deadend_split)
    for (i, m) in enumerate(accepted)
        @test all(==(EnzymeRates.Derivation._NO_FAULT), _testhelper_face_faults(m))
        @test all(==(:agree),
                  _testhelper_respeed_verdicts(m; rng = Random.MersenneTwister(i)))
        @test EnzymeRates.Derivation._guard(m) === EnzymeRates.Derivation._WELL_POSED
    end
    # Where the structure leaves a split, the limit test finds the law exact.
    limits(m) = [v.limit for v in _testhelper_face_verdicts(m)]
    @test limits(lumpable_172)[[2, 7]] ==
          [EnzymeRates.Derivation._LIMIT_EQUAL, EnzymeRates.Derivation._LIMIT_EQUAL]
    @test EnzymeRates.Derivation._LIMIT_EQUAL in limits(lumpable_532)
    @test limits(deadend_split)[2] === EnzymeRates.Derivation._LIMIT_EQUAL
end

@testset "class 3 rejects a trap and names the step that feeds it" begin
    ok, trap, split = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._TRAP,
        EnzymeRates.Derivation._SPLIT
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
    faults = _testhelper_face_faults(ldh)
    @test faults == [ok, trap, ok, ok, ok, ok, split]
    v = EnzymeRates.Derivation._guard(ldh)
    @test (v.fault, v.limit, v.face, v.conformation, v.edge) ==
          (trap, EnzymeRates.Derivation._LIMIT_SPEED, 2, 1,
           _testhelper_edge(ldh, "EP + Q ⇌ EPQ"))
    @test EnzymeRates.Derivation._guard_message(ldh, v) ==
          "Mechanism: at zero P and Q, the rapid-equilibrium step E(P) + Q ⇌ E(P, Q) " *
          "fills the forms {E(P)}, which no rapid-equilibrium step running there " *
          "empties, so the rate law misweights them. Make that step's kinetic group " *
          "steady-state."
    verdicts = _testhelper_respeed_verdicts(ldh; rng = Random.MersenneTwister(1))
    @test verdicts[[2, 7]] == [:pattern, :pattern]
    @test all(_testhelper_agrees.(faults, verdicts))
    # The fix: the feeding step steady-state.
    ldh_fed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> E(P, Q)
            E(P) + Q <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + P <--> E(P)
            E + Q ⇌ E(Q)
        end
    end
    @test EnzymeRates.Derivation._guard(ldh_fed) === EnzymeRates.Derivation._WELL_POSED
    @test all(==(:agree),
              _testhelper_respeed_verdicts(ldh_fed; rng = Random.MersenneTwister(2)))
    child9 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E + P ⇌ E(P)
            E + Q <--> E(Q)
            E(A) + B ⇌ E(A, B)
            E(A) + Q <--> E(A, Q)
            E(P) + B ⇌ E(B, P)
            E(P) + Q ⇌ E(P, Q)
            E(Q) + A <--> E(A, Q)
            E(Q) + P ⇌ E(P, Q)
            E(A, B) <--> E(P, Q)
        end
    end
    faults9 = _testhelper_face_faults(child9)
    @test faults9 == [ok, trap, ok, ok, ok, split, ok]
    v9 = _testhelper_face_verdicts(child9)[2]
    @test v9.edge == _testhelper_edge(child9, "EQ + P ⇌ EPQ")
    @test occursin("the rapid-equilibrium step E(Q) + P ⇌ E(P, Q) fills the forms",
                   EnzymeRates.Derivation._guard_message(child9, v9))
    verdicts9 = _testhelper_respeed_verdicts(child9; rng = Random.MersenneTwister(3))
    @test verdicts9[[2, 6]] == [:pattern, :pattern]
    @test all(_testhelper_agrees.(faults9, verdicts9))
    child9_fed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E + P ⇌ E(P)
            E + Q <--> E(Q)
            E(A) + B ⇌ E(A, B)
            E(A) + Q <--> E(A, Q)
            E(P) + B ⇌ E(B, P)
            E(P) + Q ⇌ E(P, Q)
            E(Q) + A <--> E(A, Q)
            E(Q) + P <--> E(P, Q)
            E(A, B) <--> E(P, Q)
        end
    end
    @test EnzymeRates.Derivation._guard(child9_fed) === EnzymeRates.Derivation._WELL_POSED
    @test all(==(:agree),
              _testhelper_respeed_verdicts(child9_fed; rng = Random.MersenneTwister(4)))
end

@testset "class 3 rejects a split, with no fix" begin
    ok, split = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._SPLIT
    split_child1 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E + B <--> E(B)
            E + Q ⇌ E(Q)
            E(A) + B <--> E(P, Q)
            E(B) + Q ⇌ E(B, Q)
            E(B) + A <--> E(P, Q)
            E(Q) + B ⇌ E(B, Q)
            E(Q) + P <--> E(P, Q)
        end
    end
    faults = _testhelper_face_faults(split_child1)
    @test faults == [ok, split, ok, ok, ok, ok, split]
    at_Q = _testhelper_face_verdicts(split_child1)[7]
    msg = EnzymeRates.Derivation._guard_message(split_child1, at_Q)
    @test startswith(msg, "Mechanism: at zero Q, the forms ")
    @test endswith(msg, "so the rate law cannot weight one against the other.")
    @test !occursin("Make", msg)
    verdicts = _testhelper_respeed_verdicts(split_child1; rng = Random.MersenneTwister(5))
    @test verdicts[[2, 7]] == [:pattern, :pattern]
    @test all(_testhelper_agrees.(faults, verdicts))
    # LDH-like at Q = 0: the lowest forms of free E's segment fall into two classes.
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
    @test EnzymeRates.Derivation._guard_message(ldh, _testhelper_face_verdicts(ldh)[7]) ==
          "Mechanism: at zero Q, the forms {E, E(A)} and {E(P)} of one rapid-equilibrium " *
          "segment are joined only by steps that take up a zeroed metabolite " *
          "(E + Q ⇌ E(Q) and E(P) + Q ⇌ E(P, Q)), so the rate law cannot weight one " *
          "against the other."
end

@testset "a split across both conformations names each form's conformation" begin
    # Both halves of the ping-pong :NonequalAI run chemistry in both conformations on one
    # rapid-equilibrium segment, joined through the flip at free E. At zero ATP and F6P
    # the covalent forms of each conformation make a part of their own, joined to the
    # rest only by steps that take up a zeroed metabolite. The two parts hold the same
    # forms, one part in each conformation, so the message names each form's
    # conformation and none for the whole split.
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: ATP, F6P
        products: ADP, F16BP
        catalytic_multiplicity: 1
        steps: begin
            E + ATP ⇌ E(ATP) :: EqualAI
            E(ATP) <--> E(F16BP; residual = ATP - F16BP) :: NonequalAI
            E(; residual = ATP - F16BP) + F16BP ⇌ E(F16BP; residual = ATP - F16BP):: EqualAI
            E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP) :: EqualAI
            E(F6P; residual = ATP - F16BP) ⇌ E(ADP) :: NonequalAI
            E + ADP ⇌ E(ADP) :: EqualAI
        end
    end
    v = EnzymeRates.Derivation._guard(am)
    @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._SPLIT, 3, 0)
    @test EnzymeRates.Derivation._guard_message(am, v) ==
          "AllostericMechanism: at zero ATP and F6P, the forms " *
          "{E(; residual = ATP - F16BP) (active), " *
          "E(F16BP; residual = ATP - F16BP) (active)} and " *
          "{E(; residual = ATP - F16BP) (inactive), " *
          "E(F16BP; residual = ATP - F16BP) (inactive)} of one rapid-equilibrium segment " *
          "are joined only by steps that take up a zeroed metabolite " *
          "(E(; residual = ATP - F16BP) + F6P ⇌ E(F6P; residual = ATP - F16BP)), so " *
          "the rate law cannot weight one against the other."
    # The LDH-like split at zero Q, every group :EqualAI: one part holds free E and E(A)
    # of both conformations, joined through the flip, so the labels go on the forms.
    ldh = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        catalytic_multiplicity: 1
        steps: begin
            E + A ⇌ E(A) :: EqualAI
            E(A) + B <--> E(A, B) :: EqualAI
            E(A, B) <--> E(P, Q) :: EqualAI
            E(P) + Q ⇌ E(P, Q) :: EqualAI
            E(Q) + P ⇌ E(P, Q) :: EqualAI
            E + P <--> E(P) :: EqualAI
            E + Q ⇌ E(Q) :: EqualAI
        end
    end
    at_Q = _testhelper_face_verdicts(ldh)[7]
    @test (at_Q.fault, at_Q.conformation) == (EnzymeRates.Derivation._SPLIT, 0)
    @test EnzymeRates.Derivation._guard_message(ldh, at_Q) ==
          "AllostericMechanism: at zero Q, the forms {E (active), E(A) (active), " *
          "E (inactive), E(A) (inactive)} and {E(P) (active)} and {E(P) (inactive)} of " *
          "one rapid-equilibrium segment are joined only by steps that take up a zeroed " *
          "metabolite (E + Q ⇌ E(Q) and E(P) + Q ⇌ E(P, Q)), so the rate law cannot " *
          "weight one against the other."
end

@testset "a split whose parts lie in the inactive conformation is placed there" begin
    # A ping-pong whose X- and Y-bound routes run only in the inactive conformation. Free
    # E's segment starts at active E, but at zero P and Q its lowest parts, the X- and
    # Y-bound covalent forms, lie in the inactive conformation alone.
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: X, Y
        steps: begin
            E + A ⇌ E(A) :: EqualAI
            E(A) <--> E(P; residual = A - P) :: EqualAI
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P :: EqualAI
            E(; residual = A - P) + B ⇌ E(B; residual = A - P) :: EqualAI
            E(B; residual = A - P) <--> E(Q) :: EqualAI
            E(Q) ⇌ E + Q :: EqualAI
            E(A) + X ⇌ E(A, X) :: OnlyI
            E(A, X) ⇌ E(P, X; residual = A - P) :: EqualAI
            E(P, X; residual = A - P) ⇌ E(X; residual = A - P) + P :: EqualAI
            E(X; residual = A - P) + B ⇌ E(B, X; residual = A - P) :: EqualAI
            E(B, X; residual = A - P) <--> E(Q, X) :: EqualAI
            E(Q) + X ⇌ E(Q, X) :: OnlyI
            E(A) + Y ⇌ E(A, Y) :: OnlyI
            E(A, Y) ⇌ E(P, Y; residual = A - P) :: EqualAI
            E(P, Y; residual = A - P) ⇌ E(Y; residual = A - P) + P :: EqualAI
            E(Y; residual = A - P) + B ⇌ E(B, Y; residual = A - P) :: EqualAI
            E(B, Y; residual = A - P) <--> E(Q, Y) :: EqualAI
            E(Q) + Y ⇌ E(Q, Y) :: OnlyI
        end
    end
    v = EnzymeRates.Derivation._guard(am)
    @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._SPLIT, 2, 2)
    @test EnzymeRates.Derivation._guard_message(am, v) ==
          "AllostericMechanism: at zero P and Q in the inactive conformation, the forms " *
          "{E(B, X; residual = A - P), E(X; residual = A - P)} and " *
          "{E(B, Y; residual = A - P), E(Y; residual = A - P)} of one rapid-equilibrium " *
          "segment are joined only by steps that take up a zeroed metabolite " *
          "(E(X; residual = A - P) + P ⇌ E(P, X; residual = A - P) and " *
          "E(Y; residual = A - P) + P ⇌ E(P, Y; residual = A - P)), so the rate law " *
          "cannot weight one against the other."
end

@testset "a bottomless segment across both conformations labels each form" begin
    # At zero P and X the rapid-equilibrium segment of free E holds active E and, through
    # the flip, inactive E and the X-bound forms only the inactive conformation reaches.
    tc_re = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: X
        steps: begin
            E + S <--> E(S)          :: EqualAI
            E(S) <--> E(P)           :: EqualAI
            E(P) <--> E + P          :: EqualAI
            E + X ⇌ E(X)             :: OnlyI
            E(X) + S ⇌ F(X) + P      :: EqualAI
            F(X) <--> E(X)           :: EqualAI
        end
    end
    v = EnzymeRates.Derivation._guard(tc_re)
    @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._BOTTOMLESS, 5, 0)
    @test EnzymeRates.Derivation._guard_message(tc_re, v) ==
          "AllostericMechanism: at zero P and X, the rapid-equilibrium segment " *
          "{E (active), E (inactive), E(X) (inactive), F(X) (inactive)} has no form " *
          "whose weight keeps the segment's lowest power of each zeroed metabolite, so " *
          "the forms' shares depend on how those concentrations vanish and the rate is " *
          "0/0 there."
    # No fix is printed: the segment's one binding made steady-state leaves a law with no
    # value on that face.
    tc_ss = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: X
        steps: begin
            E + S <--> E(S)          :: EqualAI
            E(S) <--> E(P)           :: EqualAI
            E(P) <--> E + P          :: EqualAI
            E + X <--> E(X)          :: OnlyI
            E(X) + S ⇌ F(X) + P      :: EqualAI
            F(X) <--> E(X)           :: EqualAI
        end
    end
    w = EnzymeRates.Derivation._guard(tc_ss)
    @test (w.fault, w.limit, w.face) ==
          (EnzymeRates.Derivation._LAW_MULTI, EnzymeRates.Derivation._LIMIT_LAW_NAN, 5)
end

@testset "a trap across both conformations labels each form" begin
    # Every step but the binding of B at E(P) runs in the active conformation alone. At
    # zero P and Q the fast class of free E holds E of both conformations, joined by the
    # flip, and the step E + P ⇌ E(P) that joins it to the rest of its segment takes up P.
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + B <--> E(B)                                         :: OnlyA
            (E + P ⇌ E(P), E(B) + P ⇌ E(B, P), E(Q) + P ⇌ E(P, Q))  :: OnlyA
            (E + Q <--> E(Q), E(P) + Q <--> E(P, Q))                :: OnlyA
            E(B) + A ⇌ E(P, Q)                                      :: OnlyA
            E(P) + B ⇌ E(B, P)                                      :: EqualAI
        end
    end
    # E(B) + A ⇌ E(P, Q) takes up A into E(P, Q), which the inactive conformation cuts
    # off, so it loses its affinity there and stays fast. The bridge its cut-off forms
    # open at the interior is set aside, and the harsh ε-path oracle converges there, so
    # the guard reports class 3's verdict at zero P and Q: the trap.
    v = EnzymeRates.Derivation._guard(am)
    @test v == _testhelper_face_verdicts(am)[2]
    @test (v.fault, v.limit, v.face, v.conformation, v.edge) ==
          (EnzymeRates.Derivation._TRAP, EnzymeRates.Derivation._LIMIT_DIFFER, 2, 0,
           _testhelper_edge(am, "E + P ⇌ EP"))
    @test _testhelper_epath_verdict(am; rng = Random.MersenneTwister(33),
                                    corners = true) === :converges
    @test all(_testhelper_admission_agreement(am; rng = Random.MersenneTwister(33)))
    @test EnzymeRates.Derivation._guard_message(am, v) ==
          "AllostericMechanism: at zero P and Q, the forms {E (active), E (inactive)} " *
          "hold enzyme, but the rapid-equilibrium step E + P ⇌ E(P) that joins them to " *
          "the rest of their segment takes up a zeroed metabolite, so the rate law " *
          "misweights them. Make that step's kinetic group steady-state."
    # The face condition alone refuses the interior for that bridge; its message names
    # both tagged steps the cheapest tree needs.
    D = EnzymeRates.Derivation
    eps = D._eps_network(am, D._re_flags(am.catalytic), D._built_reading(am))
    w = D._face_condition(eps, D._law_network(am).present, UInt64(0),
                          D._turnover_name(am.catalytic))
    @test (w.fault, w.conformation, w.edge) ==
          (D._BRIDGE, 2, _testhelper_edge(am, "EB + P ⇌ EBP"))
    @test D._guard_message(am, D._Verdict(w.fault, w.limit, 1, w.conformation, w.edge)) ==
          "AllostericMechanism: the :OnlyA steps E(B) + P ⇌ E(B, P); E(Q) + P ⇌ E(P, Q) " *
          "bridge forms of the inactive conformation that its other steps keep apart, so " *
          "admission cannot show that the limit does not depend on how fast those steps " *
          "vanish. Make the group of E(B) + P ⇌ E(B, P) steady-state."
    faults = _testhelper_face_faults(am)
    verdicts = _testhelper_respeed_verdicts(am; rng = Random.MersenneTwister(33))
    @test verdicts[[2, 6]] == [:mismatch, :mismatch]
    @test all(_testhelper_agrees.(faults, verdicts))
    # The fix: the binding of P steady-state.
    fixed = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + B <--> E(B)                                                 :: OnlyA
            (E + P <--> E(P), E(B) + P <--> E(B, P), E(Q) + P <--> E(P, Q)) :: OnlyA
            (E + Q <--> E(Q), E(P) + Q <--> E(P, Q))                        :: OnlyA
            E(B) + A ⇌ E(P, Q)                                              :: OnlyA
            E(P) + B ⇌ E(B, P)                                              :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._guard(fixed) === EnzymeRates.Derivation._WELL_POSED
    @test all(==(:agree),
              _testhelper_respeed_verdicts(fixed; rng = Random.MersenneTwister(34)))
end

@testset "a bottomless segment no binding's steady-state variant cures prints no fix" begin
    # Random release of three products at rapid equilibrium, the last releases steady
    # state: at zero products every form of the segment holds a product, and making any
    # one of its bindings steady-state leaves the segment bottomless.
    terbot = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A
        products: P, Q, R
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P, Q, R)
            E(Q, R) + P ⇌ E(P, Q, R)
            E(P, R) + Q ⇌ E(P, Q, R)
            E(P, Q) + R ⇌ E(P, Q, R)
            E(R) + Q ⇌ E(Q, R)
            E(Q) + R ⇌ E(Q, R)
            E(R) + P ⇌ E(P, R)
            E(P) + R ⇌ E(P, R)
            E(Q) + P ⇌ E(P, Q)
            E(P) + Q ⇌ E(P, Q)
            E + P <--> E(P)
            E + Q <--> E(Q)
            E + R <--> E(R)
        end
    end
    v = EnzymeRates.Derivation._guard(terbot)
    @test (v.fault, v.face) == (EnzymeRates.Derivation._BOTTOMLESS, 2)
    @test EnzymeRates.Derivation._guard_message(terbot, v) ==
          "Mechanism: at zero P, Q and R, the rapid-equilibrium segment " *
          "{E(P), E(P, Q), E(P, R), E(P, Q, R), E(Q), E(R), E(Q, R)} has no form whose " *
          "weight keeps the segment's lowest power of each zeroed metabolite, so the " *
          "forms' shares depend on how those concentrations vanish and the rate is 0/0 " *
          "there."
    terbot_ss = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A
        products: P, Q, R
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P, Q, R)
            E(Q, R) + P ⇌ E(P, Q, R)
            E(P, R) + Q ⇌ E(P, Q, R)
            E(P, Q) + R ⇌ E(P, Q, R)
            E(R) + Q ⇌ E(Q, R)
            E(Q) + R ⇌ E(Q, R)
            E(R) + P ⇌ E(P, R)
            E(P) + R ⇌ E(P, R)
            E(Q) + P ⇌ E(P, Q)
            E(P) + Q <--> E(P, Q)
            E + P <--> E(P)
            E + Q <--> E(Q)
            E + R <--> E(R)
        end
    end
    w = EnzymeRates.Derivation._guard(terbot_ss)
    @test (w.fault, w.face) == (EnzymeRates.Derivation._BOTTOMLESS, 2)
end

@testset "class 3 rejects a bottomless segment" begin
    ok, bottomless = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._BOTTOMLESS
    entered = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(P) + Q ⇌ E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + P <--> E(P)
            E + Q <--> E(Q)
        end
    end
    faults = _testhelper_face_faults(entered)
    @test faults == [ok, bottomless, ok, ok, ok, ok, ok]
    v = EnzymeRates.Derivation._guard(entered)
    @test (v.fault, v.face) == (bottomless, 2)
    @test EnzymeRates.Derivation._guard_message(entered, v) ==
          "Mechanism: at zero P and Q, the rapid-equilibrium segment " *
          "{E(P), E(P, Q), E(Q)} has no form whose weight keeps the segment's lowest " *
          "power of each zeroed metabolite, so the forms' shares depend on how those " *
          "concentrations vanish and the rate is 0/0 there. Make the group of " *
          "E(P) + Q ⇌ E(P, Q) steady-state."
    verdicts = _testhelper_respeed_verdicts(entered; rng = Random.MersenneTwister(6),
                                            law = :limit)
    @test verdicts[2] === :pattern
    @test all(_testhelper_agrees.(faults, verdicts))
    # The printed fix passes the guard.
    entered_named = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(P) + Q <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + P <--> E(P)
            E + Q <--> E(Q)
        end
    end
    @test EnzymeRates.Derivation._guard(entered_named) ===
          EnzymeRates.Derivation._WELL_POSED
    # So does the segment's other binding made steady-state.
    entered_fixed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(P) + Q ⇌ E(P, Q)
            E(Q) + P <--> E(P, Q)
            E + P <--> E(P)
            E + Q <--> E(Q)
        end
    end
    @test EnzymeRates.Derivation._guard(entered_fixed) ===
          EnzymeRates.Derivation._WELL_POSED
    @test all(==(:agree),
              _testhelper_respeed_verdicts(entered_fixed; rng = Random.MersenneTwister(7)))
    # Catalysis never feeds this segment, so the face's steady state has the limit the
    # oracle's own law finds; the package law is 0/0 there all the same, and the guard
    # rejects the face on purpose.
    unfed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> E + P + Q
            E + P <--> E(P)
            E(P) + Q ⇌ E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q <--> E(Q)
        end
    end
    @test _testhelper_face_faults(unfed)[2] === bottomless
    @test _testhelper_respeed_verdicts(unfed; rng = Random.MersenneTwister(8),
                                       law = :limit)[2] === :agree
    # The package law is 0/0 at zero products: the guard rejects the face as bottomless.
    @test !_testhelper_is_valid(unfed)
    v = EnzymeRates.Derivation._guard(unfed)
    @test (v.fault, v.face) == (bottomless, 2)
end

@testset "class 3 rejects a second closed class on a face" begin
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
    # At P = 0 the dead end F(P, S) is closed and unreachable from E.
    unweighted = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
            E + P <--> E(P)
            E(P) + P <--> E(P, P)
            E(P, P) + S <--> F(P, S) + P
        end
    end
    weighted = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
            E + P <--> E(P)
            E(P) + S <--> F(S) + P
        end
    end
    detour = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: X
        steps: begin
            E + S ⇌ E(S)              :: EqualAI
            E(S) <--> E + P           :: EqualAI
            E + P ⇌ E(P)              :: EqualAI
            E(P) + P ⇌ E(P, P)        :: EqualAI
            F(P) + P <--> E(P, P)     :: EqualAI
            E + X ⇌ E(X)              :: OnlyA
            F(P) + X <--> F(P, X)     :: OnlyA
            E(X) + S ⇌ E(S, X)        :: EqualAI
            E(S, X) <--> F(P, X)      :: EqualAI
        end
    end
    cases = ((fuzz5, [2, 5]), (unweighted, [2]), (weighted, [2]), (detour, [2, 5]))
    for (i, (m, zero_products)) in enumerate(cases)
        faults = _testhelper_face_faults(m)
        @test all(==(EnzymeRates.Derivation._NONUNIQUE), faults[zero_products])
        # detour's package law is 0/0 at face 5, where the package oracle reads NaN.
        verdicts = _testhelper_respeed_verdicts(m; rng = Random.MersenneTwister(10 + i),
                                                law = m === detour ? :limit : :package)
        @test all(==(:nonunique), verdicts[zero_products])
        @test all(_testhelper_agrees.(faults, verdicts))
        @test all(v -> _testhelper_names_its_step(m, v),
                  _testhelper_face_verdicts(m)[zero_products])
    end
    # The message names the closed set and the steps out of it, with no fix.
    v = _testhelper_face_verdicts(weighted)[2]
    leaving = only(k for (k, e) in enumerate(weighted.steps)
                   if weighted.nodes[e.from].conformation !=
                      weighted.nodes[e.to].conformation)
    @test v.edge == leaving
    @test EnzymeRates.Derivation._guard_message(weighted, v) ==
          "Mechanism: at zero P, the forms {F(S)} keep enzyme away from free E: every " *
          "step that leaves them takes up a zeroed metabolite " *
          "($(EnzymeRates.Mechanisms.step_text(weighted, leaving))), so the steady " *
          "state is not unique."
    # An allosteric message names the conformation the closed set lies in.
    at_P, at_PX = _testhelper_face_verdicts(detour)[[2, 5]]
    @test EnzymeRates.Derivation._guard_message(detour, at_P) ==
          "AllostericMechanism: at zero P in the inactive conformation, the forms {F(P)} " *
          "keep enzyme away from free E: every step that leaves them takes up a zeroed " *
          "metabolite (F(P) + P <--> E(P, P)), so the steady state is not unique."
    @test EnzymeRates.Derivation._guard_message(detour, at_PX) ==
          "AllostericMechanism: at zero P and X in the active conformation, the forms " *
          "{F(P)} keep enzyme away from free E: every step that leaves them takes up a " *
          "zeroed metabolite (F(P) + P <--> E(P, P) and F(P) + X <--> F(P, X)), so the " *
          "steady state is not unique."
end

@testset "the limit test decides what structure cannot" begin
    unlumped = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(B) + A <--> E(P, Q)
            (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))
            (E + Q ⇌ E(Q), E(B) + Q ⇌ E(B, Q))
            (E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q))
        end
    end
    lumpable_732 = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> E(A), E(B) + A <--> E(A, B), E(P) + A <--> E(A, P))
            (E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
            (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q))
            (E + Q ⇌ E(Q), E(B) + Q ⇌ E(B, Q), E(P) + Q ⇌ E(P, Q))
            E(A, B) <--> E(P, Q)
        end
    end
    unfed_trap = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> E(Q) + P
            E + Q ⇌ E(Q)
            E(Q) + P ⇌ E(P, Q)
            E(P) + Q ⇌ E(P, Q)
            E + P <--> E(P)
        end
    end
    dead_end = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q <--> E(Q)
            E(P) + Q <--> E(P, Q)
        end
    end
    # At zero P the aggregated chain has two closed classes, and so does the law's segment
    # chain, so the face stops at `_LAW_MULTI`.
    law_multi = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P ⇌ E(P)
            E(P) <--> F + P
            F + S ⇌ F(S)
            F(S) <--> F(P)
            F + P ⇌ F(P)
            F(P) <--> E + P
        end
    end
    two_transients = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> E(P, Q)
            E(P, Q) ⇌ G(P, Q)
            E(Q) + P ⇌ G(P, Q)
            E(P) + Q ⇌ G(P, Q)
            E + P <--> E(P)
            E + Q ⇌ E(Q)
        end
    end
    split, trap, multi = EnzymeRates.Derivation._SPLIT, EnzymeRates.Derivation._TRAP,
        EnzymeRates.Derivation._LAW_MULTI
    rows =
        ((unlumped, [2, 7], [split, split], EnzymeRates.Derivation._LIMIT_SPEED, :pattern),
         (lumpable_732, [7], [split], EnzymeRates.Derivation._LIMIT_DIFFER, :mismatch),
         (unfed_trap, [7], [split], EnzymeRates.Derivation._LIMIT_SPEED, :pattern),
         (dead_end, [2], [multi], EnzymeRates.Derivation._LIMIT_LAW_NAN, :law_nan),
         (law_multi, [2], [multi], EnzymeRates.Derivation._LIMIT_SPEED, :pattern),
         (two_transients, [2, 7], [trap, split], EnzymeRates.Derivation._LIMIT_SPEED,
          :pattern))
    for (i, (m, faces, faults, limit, oracle)) in enumerate(rows)
        verdicts = _testhelper_face_verdicts(m)
        @test [verdicts[f].fault for f in faces] == faults
        @test all(f -> verdicts[f].limit === limit, faces)
        # dead_end's package law is 0/0 at zero products, where the package oracle reads
        # NaN; the oracle's own deletion limit has no value there either.
        oracles = _testhelper_respeed_verdicts(m; rng = Random.MersenneTwister(20 + i),
                                               law = m === dead_end ? :limit : :package)
        @test all(==(oracle), oracles[faces])
        @test all(_testhelper_agrees.([v.fault for v in verdicts], oracles))
        @test all(f -> _testhelper_names_its_step(m, verdicts[f]), faces)
    end
    @test !_testhelper_is_valid(dead_end)
    v = EnzymeRates.Derivation._guard(dead_end)
    @test (v.fault, v.face) == (multi, 2)
    # Where the law's segment chain has two closed classes, the message names the face and
    # what the limit test found, with no fix.
    message(m, f) =
        EnzymeRates.Derivation._guard_message(m, _testhelper_face_verdicts(m)[f])
    @test message(dead_end, 2) ==
          "Mechanism: at zero P and Q, the rate law has no value (0/0), while the steady " *
          "state does."
    @test message(law_multi, 2) ==
          "Mechanism: at zero P, the steady state depends on how fast the " *
          "rapid-equilibrium steps run, which the rate law does not carry."
    # unfed_trap at zero products: no transition feeds the trap there, and the limit test
    # finds the law exact.
    @test _testhelper_face_verdicts(unfed_trap)[2] === EnzymeRates.Derivation._Verdict(
        EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._LIMIT_EQUAL, 2, 0, 0)
end

@testset "a law chain whose fast limit differs from the law is rejected" begin
    # At zero P and Q the law's segment chain has two closed classes, and the steady
    # state's fast-equilibrium limit, the same at new speeds, differs from the law's value
    # as P and Q vanish. The RE-speed oracle compares the limit with its own law's limit.
    differ = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(P; residual = A - P)                          :: EqualAI
            E + Q <--> E(Q)                                         :: NonequalAI
            E(; residual = A - P) + P <--> E(P; residual = A - P)   :: NonequalAI
            E(; residual = A - P) + B ⇌ E(Q)                        :: OnlyA
            E(P; residual = A - P) + B ⇌ E(B, P; residual = A - P)  :: NonequalAI
            E(Q) + P ⇌ E(B, P; residual = A - P)                    :: NonequalAI
        end
    end
    # The :OnlyA binding of B into E(Q), which the inactive conformation's Q binding still
    # reaches, cannot be switched off at rapid equilibrium, so the guard rejects the
    # interior first; class 3 sees the law chain at zero P and Q.
    v = EnzymeRates.Derivation._guard(differ)
    k = _testhelper_edge(differ, "E_res_plusA_minusP + B ⇌ EQ")
    @test v == EnzymeRates.Derivation._Verdict(
        EnzymeRates.Derivation._FAST_STEP_PINNED, EnzymeRates.Derivation._LIMIT_UNTESTED,
        0, 2, k)
    @test EnzymeRates.Derivation._guard_message(differ, v) ==
          _testhelper_pinned_cause(differ, k, 2)
    law = _testhelper_face_verdicts(differ)[2]
    @test law == EnzymeRates.Derivation._Verdict(
        EnzymeRates.Derivation._LAW_MULTI, EnzymeRates.Derivation._LIMIT_DIFFER, 2, 0, 0)
    @test EnzymeRates.Derivation._guard_message(differ, law) ==
          "AllostericMechanism: at zero P and Q, the steady state's fast-equilibrium " *
          "limit differs from the rate law."
    @test _testhelper_respeed_verdicts(differ; rng = Random.MersenneTwister(1),
                                       law = :limit)[2] === :mismatch
    # With that binding steady-state it is stopped, the interior passes, and the guard
    # itself reaches the law chain at zero P and Q; both oracles agree on every face.
    differ_ss = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(P; residual = A - P)                          :: EqualAI
            E + Q <--> E(Q)                                         :: NonequalAI
            E(; residual = A - P) + P <--> E(P; residual = A - P)   :: NonequalAI
            E(; residual = A - P) + B <--> E(Q)                     :: OnlyA
            E(P; residual = A - P) + B ⇌ E(B, P; residual = A - P)  :: NonequalAI
            E(Q) + P ⇌ E(B, P; residual = A - P)                    :: NonequalAI
        end
    end
    v = EnzymeRates.Derivation._guard(differ_ss)
    @test v == EnzymeRates.Derivation._Verdict(
        EnzymeRates.Derivation._LAW_MULTI, EnzymeRates.Derivation._LIMIT_DIFFER, 2, 0, 0)
    @test EnzymeRates.Derivation._guard_message(differ_ss, v) ==
          "AllostericMechanism: at zero P and Q, the steady state's fast-equilibrium " *
          "limit differs from the rate law."
    @test all(_testhelper_admission_agreement(differ_ss; rng = Random.MersenneTwister(2)))
end

@testset "a bridge's message names the tagged steps it needs" begin
    # E + A ⇌ E(P; residual = A - P) binds A into a form the inactive conformation cuts
    # off, so it loses its affinity and stays fast: the cheapest spanning in-tree of that
    # conformation's fed forms needs only the vanishing :OnlyA binding of A at E(Q), and
    # the message names it alone (the trap test's mechanism names a bridge of two). At zero
    # P and Q the ε-path converges: the rejection is conservative, as the message allows.
    bridge2 = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(P; residual = A - P)                          :: OnlyA
            E + Q <--> E(Q)                                         :: EqualAI
            E(; residual = A - P) + P <--> E(P; residual = A - P)   :: OnlyA
            E(; residual = A - P) + B <--> E(Q)                     :: NonequalAI
            E(P; residual = A - P) + Q ⇌ E(A, Q)                    :: EqualAI
            E(Q) + A ⇌ E(A, Q)                                      :: OnlyA
        end
    end
    v = EnzymeRates.Derivation._guard(bridge2)
    @test (v.fault, v.face, v.conformation, v.edge) ==
          (EnzymeRates.Derivation._BRIDGE, 1, 2, _testhelper_edge(bridge2, "EQ + A ⇌ EAQ"))
    @test EnzymeRates.Derivation._guard_message(bridge2, v) ==
          "AllostericMechanism: the :OnlyA step E(Q) + A ⇌ E(A, Q) bridges forms of the " *
          "inactive conformation that its other steps keep apart, so admission cannot " *
          "show that the limit does not depend on how fast that step vanishes. Make the " *
          "group of E(Q) + A ⇌ E(A, Q) steady-state."
    @test _testhelper_epath_verdict(bridge2; rng = Random.MersenneTwister(2),
                                    face = Set([:P, :Q])) === :converges
    # With the binding of A into E(P; residual = A - P) steady-state, the guard's first
    # built reading still finds a bridge at zero P and Q; the search admits the mechanism
    # under another reading, and the oracles agree on every face.
    leading_ss = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(P; residual = A - P)                       :: OnlyA
            E + Q <--> E(Q)                                         :: EqualAI
            E(; residual = A - P) + P <--> E(P; residual = A - P)   :: OnlyA
            E(; residual = A - P) + B <--> E(Q)                     :: NonequalAI
            E(P; residual = A - P) + Q ⇌ E(A, Q)                    :: EqualAI
            E(Q) + A ⇌ E(A, Q)                                      :: OnlyA
        end
    end
    w = EnzymeRates.Derivation._admission(
        leading_ss, EnzymeRates.Derivation._built_reading(leading_ss))
    @test (w.fault, w.face, w.conformation) == (EnzymeRates.Derivation._BRIDGE, 2, 2)
    @test EnzymeRates.Derivation._guard(leading_ss) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(leading_ss; rng = Random.MersenneTwister(3),
                                              corners = true))
    _testhelper_check_law_against_limit(EnzymeRates.compile_mechanism(leading_ss))
end

@testset "the modifier faces catch a steady-state inhibitor binding" begin
    # At I = 0 the rapid-equilibrium square through E(I) and E(I, S) no longer joins E
    # and E(S), whose ratio the steady-state binding of S then sets.
    mixed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S <--> E(S)
            E(S) <--> E + P
            (E + I ⇌ E(I), E(S) + I ⇌ E(I, S))
            E(I) + S ⇌ E(I, S)
        end
    end
    faults = _testhelper_face_faults(mixed)
    @test all(==(EnzymeRates.Derivation._NO_FAULT), faults[1:3])
    @test all(!=(EnzymeRates.Derivation._NO_FAULT), faults[4:6])
    verdicts = _testhelper_respeed_verdicts(mixed; rng = Random.MersenneTwister(31))
    @test verdicts[4:6] == [:mismatch, :mismatch, :mismatch]
    @test all(_testhelper_agrees.(faults, verdicts))
    @test all(v -> _testhelper_names_its_step(mixed, v),
              _testhelper_face_verdicts(mixed)[4:6])
    # At zero I and S the slow reverse catalysis fills E(S), and the rapid-equilibrium
    # binding of I that joins it to free E's segment takes up I.
    v = _testhelper_face_verdicts(mixed)[6]
    @test (v.fault, v.limit, v.conformation, v.edge) ==
          (EnzymeRates.Derivation._TRAP, EnzymeRates.Derivation._LIMIT_DIFFER, 1,
           _testhelper_edge(mixed, "ES + I ⇌ EIS"))
    @test EnzymeRates.Derivation._guard_message(mixed, v) ==
          "Mechanism: at zero I and S, the forms {E(S)} hold enzyme, but the " *
          "rapid-equilibrium step E(S) + I ⇌ E(I, S) that joins them to the rest of " *
          "their segment takes up a zeroed metabolite, so the rate law misweights them. " *
          "Make that step's kinetic group steady-state."
    # The fix: the binding of I steady-state.
    mixed_fixed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S <--> E(S)
            E(S) <--> E + P
            (E + I <--> E(I), E(S) + I <--> E(I, S))
            E(I) + S ⇌ E(I, S)
        end
    end
    @test EnzymeRates.Derivation._guard(mixed_fixed) === EnzymeRates.Derivation._WELL_POSED
    @test all(==(:agree),
              _testhelper_respeed_verdicts(mixed_fixed; rng = Random.MersenneTwister(32)))
end

@testset "the limit test draws one verdict for every writing of a mechanism" begin
    written = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> E(A), E(B) + A <--> E(A, B), E(P) + A <--> E(A, P))
            (E(A) + B ⇌ E(A, B), E(Q) + B ⇌ E(B, Q))
            (E + P <--> E(P), E(A) + P <--> E(A, P), E(Q) + P <--> E(P, Q))
            (E + Q ⇌ E(Q), E(B) + Q ⇌ E(B, Q), E(P) + Q ⇌ E(P, Q))
            E(A, B) <--> E(P, Q)
        end
    end
    # The same mechanism with its groups and steps in reverse order, each binding
    # written as its release.
    rewritten = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E(A, B) <--> E(P, Q)
            (E(P, Q) ⇌ E(P) + Q, E(B, Q) ⇌ E(B) + Q, E(Q) ⇌ E + Q)
            (E(P, Q) <--> E(Q) + P, E(A, P) <--> E(A) + P, E(P) <--> E + P)
            (E(B, Q) ⇌ E(Q) + B, E(A, B) ⇌ E(A) + B)
            (E(A, P) <--> E(P) + A, E(A, B) <--> E(B) + A, E(A) <--> E + A)
        end
    end
    @test written == rewritten
    @test hash(written) == hash(rewritten)
    @test _testhelper_face_verdicts(written) == _testhelper_face_verdicts(rewritten)
    @test _testhelper_face_verdicts(written)[7].limit ===
          EnzymeRates.Derivation._LIMIT_DIFFER
    # The residue stream is pinned so that every process and Julia version draws the same
    # residues.
    r = EnzymeRates.Derivation._ResidueStream(UInt64(0))
    @test [EnzymeRates.Derivation._residue!(r) for _ in 1:3] ==
          [153307352162749901, 1042757494553273858, 487617019471545681]
end

@testset "the guard on LDH's flip candidates: the trap, its fix and a fast cycle" begin
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
    candidate(re) = EnzymeRates.Enumeration._with_flags(ldh, re)
    flags = _testhelper_flags(ldh)
    @test EnzymeRates.Derivation._guard(candidate(flags)).fault ===
          EnzymeRates.Derivation._TRAP
    # The feeding step steady-state through the flags: the graph of the trap test's fix.
    fed = copy(flags)
    fed[_testhelper_edge(ldh, "EP + Q ⇌ EPQ")] = false
    @test _testhelper_is_valid(candidate(fed))
    # Every step at rapid equilibrium: a fast cycle turns the reaction over.
    @test EnzymeRates.Derivation._guard(candidate(trues(length(ldh.steps)))).fault ===
          EnzymeRates.Derivation._RE_TURNOVER
end

# ── Tag sets of the admission tests ─────────────────────────────────────────

"""
Tag sets by name, each built without the guard (`@_testhelper_unguarded`): those that pass
every face under the guard's first built reading, those a face rejects under it
(`_testhelper_face_rejected`), and those whose inactive conformation excludes a
rapid-equilibrium step that reading cannot switch off, so admission's interior rejects it
(`_testhelper_pinned_tag_sets`); the search admits some of the last two kinds. The
ping-pong tag sets write their
covalent form `E(; residual = A - P)`, which gives the graph a free form `F` would;
AD2a_iso_chemNonequal is R1, AD2e is R2 and AD2c is R3 written again, and AD6 holds one
reaction in two groups, which the constructor rejects.
"""
const _testhelper_tag_sets = Dict{String, EnzymeRates.AllostericMechanism}()

_testhelper_tag_sets["R1"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> E(P)      :: NonequalAI
        F + P ⇌ E(P)        :: EqualAI
        F <--> E            :: NonequalAI
    end
end
_testhelper_tag_sets["A5"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> E(P)      :: NonequalAI
        E + P ⇌ E(P)        :: EqualAI
        F ⇌ E               :: NonequalAI
    end
end
_testhelper_tag_sets["AD2a"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> E(P)      :: OnlyA
        F + P ⇌ E(P)        :: EqualAI
        F <--> E            :: EqualAI
    end
end
_testhelper_tag_sets["AD2d"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> E(P)      :: EqualAI
        F + P ⇌ E(P)        :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["AD3a"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> E(P)      :: NonequalAI
        E + P ⇌ E(P)        :: EqualAI
        E + I ⇌ E(I)        :: OnlyA
    end
end
_testhelper_tag_sets["AD3b"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: S
    steps: begin
        E + S ⇌ E(S)                                       :: EqualAI
        E(S) <--> E(P)                                     :: EqualAI
        E + P ⇌ E(P)                                       :: EqualAI
        (E + S::Mod ⇌ E(S::Mod), E(S) + S::Mod ⇌ E(S, S::Mod)) :: OnlyA
        E(S::Mod) + S ⇌ E(S, S::Mod)                       :: EqualAI
    end
end
_testhelper_tag_sets["AD3d"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: S
    steps: begin
        E + S ⇌ E(S)                    :: EqualAI
        E(S) <--> E(P)                  :: EqualAI
        E + P ⇌ E(P)                    :: EqualAI
        E + S::Mod <--> E(S::Mod)       :: OnlyA
        E(S) + S::Mod <--> E(S, S::Mod) :: OnlyA
        E(S::Mod) + S ⇌ E(S, S::Mod)    :: EqualAI
    end
end
_testhelper_tag_sets["AD4"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P, Q
    steps: begin
        E + S ⇌ E(S)            :: EqualAI
        E(S) <--> E(P, Q)       :: EqualAI
        E(P, Q) <--> E + P + Q  :: OnlyA
    end
end
_testhelper_tag_sets["AD4b"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P, Q
    steps: begin
        E + S ⇌ E(S)            :: EqualAI
        E(S) <--> E(P, Q)       :: OnlyA
        E(P, Q) ⇌ E + P + Q     :: EqualAI
    end
end
_testhelper_tag_sets["AD5"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E <--> F                        :: EqualAI
        (E + S ⇌ E(S), F + S ⇌ F(S))    :: OnlyA
        E(S) <--> F(S)                  :: EqualAI
        F(S) <--> F(P)                  :: OnlyA
        F + P ⇌ F(P)                    :: EqualAI
    end
end
_testhelper_tag_sets["AD7a"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A
    products: P
    modifiers: I
    steps: begin
        E + A <--> E(A)                 :: EqualAI
        (E + I ⇌ E(I), F + I ⇌ E(I))    :: OnlyA
        E + P <--> E(P)                 :: EqualAI
        E ⇌ F                           :: EqualAI
        E(A) <--> E(P)                  :: EqualAI
    end
end
_testhelper_tag_sets["AD7b"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A
    products: P
    modifiers: I
    steps: begin
        E + A <--> E(A)                 :: EqualAI
        (E + I ⇌ E(I), F + I ⇌ E(I))    :: EqualAI
        E + P <--> E(P)                 :: EqualAI
        E ⇌ F                           :: OnlyA
        E(A) <--> E(P)                  :: EqualAI
    end
end
_testhelper_tag_sets["AD7c"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A
    products: P
    modifiers: I
    steps: begin
        E + A <--> E(A)                 :: EqualAI
        (E + I ⇌ E(I), F + I ⇌ E(I))    :: EqualAI
        E + P <--> E(P)                 :: EqualAI
        E ⇌ F                           :: EqualAI
        E(A) <--> E(P)                  :: OnlyA
    end
end
_testhelper_tag_sets["AD8"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ E(A)            :: EqualAI
        E(A) + B <--> E(Q) + P  :: OnlyA
        E + Q ⇌ E(Q)            :: EqualAI
        E(Q) + B ⇌ E(B, Q)      :: EqualAI
    end
end
_testhelper_tag_sets["AD9"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ E(A)            :: EqualAI
        E(A) + B ⇌ E(A, B)      :: EqualAI
        E(A, B) ⇌ E(P, Q)       :: OnlyA
        E(Q) + P <--> E(P, Q)   :: EqualAI
        E + Q <--> E(Q)         :: EqualAI
        E(A) + Q <--> E(A, Q)   :: OnlyA
    end
end
_testhelper_tag_sets["AD10"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)         :: OnlyA
        E + B ⇌ E(B)            :: EqualAI
        E(A) + B ⇌ E(A, B)      :: EqualAI
        E(B) + A <--> E(A, B)   :: OnlyA
        E(A, B) <--> E(P, Q)    :: OnlyA
        E(Q) + P ⇌ E(P, Q)      :: EqualAI
        E + Q ⇌ E(Q)            :: EqualAI
    end
end
# AD13: a second free form joined to E at rapid equilibrium, each label sorting before or
# after E, with the chemistry :NonequalAI or :OnlyA.
for label in (:D, :Ca, :Ea, :F, :Z), chem in (:NonequalAI, :OnlyA)
    block = quote
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: EqualAI
            E(S) <--> E(P)      :: $chem
            E + P ⇌ E(P)        :: EqualAI
            $label ⇌ E          :: NonequalAI
        end
    end
    _testhelper_tag_sets["AD13_$(chem)_$label"] =
        EnzymeRates.DSL._construct(EnzymeRates.DSL.raw_mechanism(block, true))
end
_testhelper_tag_sets["AD14"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ E(S)                            :: EqualAI
        E + P ⇌ E(P)                            :: EqualAI
        (E(S) <--> E(P), E(S, I) <--> E(P, I))  :: OnlyA
        E(I) + S ⇌ E(S, I)                      :: OnlyA
        E(I) + P ⇌ E(P, I)                      :: EqualAI
        E + I ⇌ E(I)                            :: EqualAI
    end
end
_testhelper_tag_sets["AD14b"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ E(S)                            :: EqualAI
        E + P ⇌ E(P)                            :: EqualAI
        (E(S) <--> E(P), E(S, I) <--> E(P, I))  :: OnlyA
        E(I) + S <--> E(S, I)                   :: OnlyA
        E(I) + P ⇌ E(P, I)                      :: EqualAI
        E + I ⇌ E(I)                            :: EqualAI
    end
end
_testhelper_tag_sets["AD14c"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ E(S)                            :: EqualAI
        E + P ⇌ E(P)                            :: EqualAI
        (E(S) <--> E(P), E(S, I) <--> E(P, I))  :: OnlyA
        E(I) + S ⇌ E(S, I)                      :: EqualAI
        E(I) + P ⇌ E(P, I)                      :: EqualAI
        E + I ⇌ E(I)                            :: EqualAI
    end
end
_testhelper_tag_sets["AD15"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        (E + A ⇌ E(A), E(B) + A ⇌ E(A, B))  :: OnlyA
        E(A) + B ⇌ E(A, B)                  :: OnlyA
        E + B ⇌ E(B)                        :: EqualAI
        E(A, B) <--> E(P, Q)                :: OnlyA
        E(Q) + P ⇌ E(P, Q)                  :: EqualAI
        E + Q ⇌ E(Q)                        :: EqualAI
    end
end
_testhelper_tag_sets["R2P_sealed"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> E + P     :: EqualAI
        E + P ⇌ E(P)        :: OnlyA
        F + P ⇌ E(P)        :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["B_uphill"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: X
    steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> E(P)      :: EqualAI
        E + P ⇌ E(P)        :: EqualAI
        E + X ⇌ E(X)        :: OnlyA
        F + X <--> E(X)     :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["TC1"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ E(A)            :: EqualAI
        E(A) + B ⇌ E(Q) + P     :: OnlyA
        E(Q) <--> E + Q         :: EqualAI
        E(A) + B <--> E(A, B)   :: EqualAI
        E(A, B) <--> E(P, Q)    :: EqualAI
        E(P, Q) ⇌ E(Q) + P      :: EqualAI
    end
end
_testhelper_tag_sets["pp_V_half1"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)                                       :: EqualAI
        E(A) <--> E(; residual = A - P) + P                   :: OnlyA
        E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
        E(B; residual = A - P) <--> E + Q                     :: EqualAI
    end
end
_testhelper_tag_sets["pp_V_half2"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)                                       :: EqualAI
        E(A) <--> E(; residual = A - P) + P                   :: EqualAI
        E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
        E(B; residual = A - P) <--> E + Q                     :: OnlyA
    end
end
_testhelper_tag_sets["pp_K_Abind_half1"] =
    @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A <--> E(A)                                       :: OnlyA
        E(A) <--> E(; residual = A - P) + P                   :: OnlyA
        E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
        E(B; residual = A - P) <--> E + Q                     :: EqualAI
    end
end
_testhelper_tag_sets["PP1"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ E(A)                                          :: EqualAI
        E(A) <--> E(; residual = A - P) + P                   :: OnlyA
        E(; residual = A - P) + B ⇌ E(B; residual = A - P)    :: EqualAI
        E(B; residual = A - P) <--> E + Q                     :: EqualAI
    end
end
_testhelper_tag_sets["PP2"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    modifiers: I
    steps: begin
        E + A <--> E(A)                                       :: EqualAI
        E(A) <--> E(; residual = A - P) + P                   :: EqualAI
        E(; residual = A - P) + B ⇌ E(B; residual = A - P)    :: EqualAI
        E(B; residual = A - P) <--> E + Q                     :: OnlyA
        E(; residual = A - P) + I ⇌ E(I; residual = A - P)    :: OnlyA
    end
end
# Census ter-bi d75d2415361c4dae: its faces with P zero pass by T0+ alone.
_testhelper_tag_sets["terbi_d75d2415361c4dae"] =
    @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B, D
    products: P, Q
    steps: begin
        (E + A ⇌ E(A), E(B, P) + A ⇌ E(A, B, P), E(P) + A ⇌ E(A, P))       :: OnlyA
        (E + P ⇌ E(P), E(A) + P ⇌ E(A, P), E(A, B) + P ⇌ E(A, B, P),
         E(Q) + P ⇌ E(P, Q))                                               :: EqualAI
        (E + Q <--> E(Q), E(P) + Q <--> E(P, Q))                           :: OnlyA
        (E(A) + B <--> E(A, B), E(A, P) + B <--> E(A, B, P),
         E(P) + B <--> E(B, P))                                            :: EqualAI
        (E(A, B) + D <--> E(P, Q), E(Q) + D <--> E(D, Q))                  :: OnlyA
    end
end
# An :OnlyA group mixing B's fused binding into E(P, Q) with its plain binding at E(Q),
# which no untagged step joins to E(B, Q): its limit exists (spec section 7, "Mixed groups
# stay allowed").
_testhelper_tag_sets["mixed_group_limit"] =
    @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    steps: begin
        E + A ⇌ E(A)                                    :: EqualAI
        (E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q))  :: OnlyA
        E(Q) + P ⇌ E(P, Q)                              :: EqualAI
        E + Q ⇌ E(Q)                                    :: EqualAI
    end
end
_testhelper_tag_sets["R2"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)        :: OnlyA
        E(S) <--> E(P)      :: EqualAI
        F + P ⇌ E(P)        :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["BS1"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> E(P)      :: EqualAI
        E(P) ⇌ E + P        :: EqualAI
        E <--> F            :: OnlyA
        E + I ⇌ F(I)        :: OnlyA
        F + I ⇌ F(I)        :: EqualAI
        F + S ⇌ F(S)        :: EqualAI
        F(S) <--> F(P)      :: EqualAI
        F(P) ⇌ F + P        :: EqualAI
    end
end
_testhelper_tag_sets["BS2"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)        :: OnlyA
        E(S) <--> E(P)      :: EqualAI
        F + P ⇌ E(P)        :: EqualAI
        F ⇌ E               :: OnlyA
    end
end
_testhelper_tag_sets["BS3"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S <--> E(S)     :: OnlyA
        E(S) <--> E(P)      :: EqualAI
        F + P ⇌ E(P)        :: EqualAI
        F <--> E            :: OnlyA
    end
end
_testhelper_tag_sets["BS4"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ E(S)            :: EqualAI
        E(S) <--> E + P         :: EqualAI
        E + I ⇌ E(I)            :: OnlyA
        E(I) + S <--> E(I, S)   :: EqualAI
        E(I, S) <--> F(I) + P   :: EqualAI
        F(I) ⇌ F + I            :: EqualAI
        F <--> E                :: OnlyA
    end
end
_testhelper_tag_sets["ADbridge"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S <--> E(S)                     :: EqualAI
        E(S) <--> E + P                     :: EqualAI
        (E + I ⇌ E(I), E(S) + I ⇌ E(I, S))  :: OnlyA
        E(I) + S ⇌ E(I, S)                  :: EqualAI
    end
end
_testhelper_tag_sets["ADbridge_noturnover"] =
    @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S <--> E(S)                     :: EqualAI
        E(S) <--> E + P                     :: OnlyA
        (E + I ⇌ E(I), E(S) + I ⇌ E(I, S))  :: OnlyA
        E(I) + S ⇌ E(I, S)                  :: EqualAI
    end
end
_testhelper_tag_sets["bridge_drain"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S ⇌ E(S)        :: EqualAI
        E(S) <--> F + P     :: EqualAI
        F <--> E            :: OnlyA
        E + I ⇌ E(I)        :: OnlyA
        F + I ⇌ F(I)        :: OnlyA
        F(I) ⇌ E(I)         :: EqualAI
    end
end
_testhelper_tag_sets["bridge_deadend"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I, J
    steps: begin
        E + S <--> E(S)         :: EqualAI
        E(S) <--> E + P         :: EqualAI
        E + J <--> E(J)         :: EqualAI
        E + I ⇌ E(I)            :: OnlyA
        E(J) + I ⇌ E(I, J)      :: OnlyA
        E(I) + J ⇌ E(I, J)      :: EqualAI
    end
end
_testhelper_tag_sets["R2funnel"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)        :: OnlyA
        E(S) <--> E(P)      :: EqualAI
        F + P ⇌ E(P)        :: EqualAI
        F <--> E            :: OnlyA
        F + S ⇌ F(S)        :: EqualAI
        F(S) <--> F(P)      :: EqualAI
        F + P ⇌ F(P)        :: EqualAI
        F(P) ⇌ E(P)         :: EqualAI
    end
end
_testhelper_tag_sets["R2funnel_SS"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S <--> E(S)     :: OnlyA
        E(S) <--> E(P)      :: EqualAI
        F + P ⇌ E(P)        :: EqualAI
        F <--> E            :: OnlyA
        F + S ⇌ F(S)        :: EqualAI
        F(S) <--> F(P)      :: EqualAI
        F + P ⇌ F(P)        :: EqualAI
        F(P) ⇌ E(P)         :: EqualAI
    end
end
_testhelper_tag_sets["S1"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: A, B
    products: P
    modifiers: A
    steps: begin
        E + A <--> E(A)                 :: OnlyA
        E(A) + B <--> E(A, B)           :: EqualAI
        E(A, B) <--> E(P)               :: OnlyA
        E(P) <--> E + P                 :: EqualAI
        E + A::Mod ⇌ E(A::Mod)          :: OnlyA
        E(A::Mod) + B ⇌ E(A::Mod, B)    :: EqualAI
        E(B) + A::Mod ⇌ E(A::Mod, B)    :: EqualAI
        F + B <--> E(B)                 :: EqualAI
    end
end
# Fuzz mechanism fz101_10232 (seed 101).
_testhelper_tag_sets["fz101_10232"] = @_testhelper_unguarded @allosteric_mechanism begin
    substrates: S
    products: P
    modifiers: I
    steps: begin
        E + S <--> E(S)                             :: OnlyA
        E(S) <--> E(P)                              :: OnlyA
        E(P) <--> F + P                             :: EqualAI
        F <--> E                                    :: OnlyA
        (E(S) + I <--> E(I, S), F + I <--> F(I))    :: OnlyA
        E(P) + S ⇌ E(S) + P                         :: EqualAI
    end
end

@testset "how an excluded step is switched off" begin
    D = EnzymeRates.Derivation
    rng = Random.MersenneTwister(2611)
    # The S binding takes up one metabolite into E(S), which the inactive conformation's
    # P binding no longer joins to E: it loses its affinity. The chemistry takes up
    # nothing, so it has no affinity to lose; at steady state it is stopped.
    k_type = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: OnlyA
            E(S) <--> E(P)      :: OnlyA
            E + P ⇌ E(P)        :: EqualAI
        end
    end
    off = D._built_reading(k_type)[:, 2]
    @test off[_testhelper_edge(k_type, "E + S ⇌ ES")] === D._LOSES_FORWARD
    @test off[_testhelper_edge(k_type, "ES <--> EP")] === D._STOPPED
    @test off[_testhelper_edge(k_type, "E + P ⇌ EP")] === D._RUNS
    @test all(==(D._RUNS), D._built_reading(k_type)[:, 1])
    @test D._guard(k_type).fault === D._NO_FAULT
    @test all(_testhelper_admission_agreement(k_type; rng))
    # With the chemistry untagged, E(S) stays, made from P: thermodynamics pins the S
    # binding's equilibrium constant, and a rapid-equilibrium step cannot turn slow.
    pinned = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: OnlyA
            E(S) <--> E(P)      :: EqualAI
            E + P ⇌ E(P)        :: EqualAI
        end
    end
    v = D._guard(pinned)
    @test (v.fault, v.conformation, v.edge) ==
          (D._FAST_STEP_PINNED, 2, _testhelper_edge(pinned, "E + S ⇌ ES"))
    @test _testhelper_epath_verdict(pinned; rng) === :no_path
    # The same binding at steady state is stopped: S neither binds nor leaves in the
    # inactive conformation, and E(S) is made from P there.
    shut_door = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(S)     :: OnlyA
            E(S) <--> E(P)      :: EqualAI
            E + P ⇌ E(P)        :: EqualAI
        end
    end
    @test D._STOPPED in D._built_reading(shut_door)[:, 2]
    @test D._guard(shut_door).fault === D._NO_FAULT
    @test all(_testhelper_admission_agreement(shut_door; rng))
    # A rapid-equilibrium step that binds B and reacts at once, into E(P, Q), which the
    # P binding still makes: it cannot be switched off.
    fused_re = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)         :: EqualAI
            E + Q ⇌ E(Q)            :: EqualAI
            E(A) + B ⇌ E(P, Q)      :: OnlyA
            E(Q) + P <--> E(P, Q)   :: EqualAI
        end
    end
    @test D._guard(fused_re).fault === D._FAST_STEP_PINNED
    # At steady state the same step is stopped.
    fused_ss = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)         :: EqualAI
            E + Q ⇌ E(Q)            :: EqualAI
            E(A) + B <--> E(P, Q)   :: OnlyA
            E(Q) + P <--> E(P, Q)   :: EqualAI
        end
    end
    @test D._guard(fused_ss).fault === D._NO_FAULT
    @test all(_testhelper_admission_agreement(fused_ss; rng))
    # The package's reading agrees with the test's own statement of the rule.
    for x in (k_type, pinned, shut_door, fused_re, fused_ss)
        @test _testhelper_reading(D._built_reading(x)) ==
              _testhelper_reading(x, :all_lose_affinity)
    end
end

@testset "a fast isomerization excluded between two groups loses a direction" begin
    # The isomerization takes up nothing, so no built reading switches it off. But the
    # steps the inactive conformation keeps do not join its two sides: E(P) lies in a group
    # of its own, and the isomerization loses its rate into E(P), staying fast. The search
    # admits the mechanism, and the oracles agree on every face.
    D = EnzymeRates.Derivation
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: EqualAI
            E(S) ⇌ E(P)         :: OnlyA
            E + P <--> E(P)     :: OnlyA
        end
    end
    iso = _testhelper_edge(am, "ES ⇌ EP")
    @test D._admission_interior(am, D._built_reading(am)) ==
          D._Verdict(D._FAST_STEP_PINNED, D._LIMIT_UNTESTED, 0, 2, iso)
    @test D._guard(am) == D._WELL_POSED
    @test D._admitting_reading(am)[iso, 2] === D._LOSES_FORWARD
    @test all(_testhelper_admission_agreement(am; rng = Random.MersenneTwister(2615),
                                              corners = true))
    _testhelper_check_law_against_limit(EnzymeRates.compile_mechanism(am))
end

@testset "a group that loses one affinity names its step that cannot be switched off" begin
    # E(A) is a dead end, cut off in the inactive conformation, so the A binding at E
    # loses its affinity; its groupmate binds A into E(P, Q), which the P binding still
    # makes, at rapid equilibrium: the verdict names the groupmate.
    D = EnzymeRates.Derivation
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(B) + A ⇌ E(P, Q))  :: OnlyA
            E + B ⇌ E(B)                        :: EqualAI
            E(Q) + P <--> E(P, Q)               :: EqualAI
            E + Q ⇌ E(Q)                        :: EqualAI
        end
    end
    lost, pinned = _testhelper_edge(am, "E + A ⇌ EA"), _testhelper_edge(am, "EB + A ⇌ EPQ")
    off = D._built_reading(am)[:, 2]
    @test (off[lost], off[pinned]) == (D._LOSES_FORWARD, D._CANNOT_SWITCH_OFF)
    v = D._guard(am)
    @test (v.fault, v.conformation, v.edge) == (D._FAST_STEP_PINNED, 2, pinned)
end

@testset "an uptake from one cut-off form into another loses its affinity" begin
    # Every step :OnlyA: the inactive conformation keeps only E. The B binding takes up
    # one metabolite from E(A), cut off, into E(A, B), cut off too, so it loses its
    # affinity, and free E reaches both along the bindings.
    D = EnzymeRates.Derivation
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)            :: OnlyA
            E(A) + B ⇌ E(A, B)      :: OnlyA
            E(A, B) <--> E(P, Q)    :: OnlyA
            E(P, Q) ⇌ E(Q) + P      :: OnlyA
            E(Q) ⇌ E + Q            :: OnlyA
        end
    end
    off = D._built_reading(am)[:, 2]
    @test off[_testhelper_edge(am, "E + A ⇌ EA")] === D._LOSES_FORWARD
    @test off[_testhelper_edge(am, "EA + B ⇌ EAB")] === D._LOSES_FORWARD
    @test D._guard(am).fault === D._NO_FAULT
    @test all(_testhelper_admission_agreement(am; rng = Random.MersenneTwister(2612)))
end

@testset "how an excluded step is switched off depends on neither subunits nor sites" begin
    # The pinned S binding and its steady-state variant of "how an excluded step is
    # switched off", with two catalytic subunits and a regulatory site.
    D = EnzymeRates.Derivation
    pinned = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)        :: OnlyA
            E(S) <--> E(P)      :: EqualAI
            E + P ⇌ E(P)        :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    shut_door = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R
        catalytic_multiplicity: 2
        steps: begin
            E + S <--> E(S)     :: OnlyA
            E(S) <--> E(P)      :: EqualAI
            E + P ⇌ E(P)        :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    v = D._guard(pinned)
    @test (v.fault, v.conformation, v.edge) ==
          (D._FAST_STEP_PINNED, 2, _testhelper_edge(pinned, "E + S ⇌ ES"))
    @test D._guard(shut_door).fault === D._NO_FAULT
end

@testset "a steady-state uptake into a cut-off form may be stopped" begin
    # The guard first reads every excluded steady-state step that takes up one metabolite
    # into a cut-off form as losing its affinity; when that leaves an admission fault, it
    # checks again with such a step losing its affinity only when it leads away from free
    # E, then with every such step stopped. Four ter-ter mechanisms of the gate's
    # populations, built unguarded since the first reading rejects each: one that the
    # second reading admits, one that only the third admits, and two whose bridges are
    # rapid-equilibrium steps, which no reading stops.
    D = EnzymeRates.Derivation
    rng = Random.MersenneTwister(2613)
    second_try_a = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + D ⇌ E(D), E(B, Q) + D ⇌ E(B, D, Q), E(P, Q) + D ⇌ E(D, P, Q),
             E(Q) + D ⇌ E(D, Q)) :: EqualAI
            (E + Q ⇌ E(Q), E(B, D) + Q ⇌ E(B, D, Q), E(D) + Q ⇌ E(D, Q),
             E(D, P) + Q ⇌ E(D, P, Q)) :: OnlyA
            E(B, D) + A <--> E(P, Q, R) :: OnlyA
            (E(D) + B ⇌ E(B, D), E(D, Q) + B ⇌ E(B, D, Q), E(Q) + B ⇌ E(B, Q)) :: OnlyA
            (E(D) + P <--> E(D, P), E(D, Q) + P <--> E(D, P, Q),
             E(Q) + P <--> E(P, Q), E(Q, R) + P <--> E(P, Q, R)) :: EqualAI
            (E(P, Q) + R ⇌ E(P, Q, R), E(Q) + R ⇌ E(Q, R)) :: EqualAI
        end
    end
    third_try = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + D ⇌ E(D), E(A, Q) + D ⇌ E(A, D, Q), E(P, Q) + D ⇌ E(D, P, Q),
             E(Q) + D ⇌ E(D, Q)) :: OnlyA
            (E + Q ⇌ E(Q), E(A, D) + Q ⇌ E(A, D, Q), E(D) + Q ⇌ E(D, Q),
             E(D, P) + Q ⇌ E(D, P, Q)) :: EqualAI
            (E(A, D) + B ⇌ E(P, Q, R), E(D) + B ⇌ E(B, D)) :: OnlyA
            (E(B, D) + A <--> E(P, Q, R), E(D) + A <--> E(A, D),
             E(D, Q) + A <--> E(A, D, Q), E(Q) + A <--> E(A, Q)) :: EqualAI
            (E(D) + P ⇌ E(D, P), E(D, Q) + P ⇌ E(D, P, Q), E(Q) + P ⇌ E(P, Q)) :: EqualAI
            E(P, Q) + R <--> E(P, Q, R) :: OnlyA
        end
    end
    # second_try_a, inactive conformation: the chemistry step joins two cut-off forms.
    k = _testhelper_edge(second_try_a, "EBD + A <--> EPQR")
    losing, stopped = D._built_reading(second_try_a)[:, 2],
                      D._built_reading(second_try_a, D._ALL_STOPPED)[:, 2]
    @test losing[k] === D._LOSES_FORWARD
    @test stopped[k] === D._STOPPED
    @test [losing[j] for j in eachindex(losing) if j != k] ==
          [stopped[j] for j in eachindex(stopped) if j != k]
    @test _testhelper_reading(second_try_a, :all_lose_affinity)[k, 2] === :loses_forward
    @test _testhelper_reading(second_try_a, :all_stopped)[k, 2] === :stopped
    # The two conformations are read alike.
    for r in (D._ALL_LOSE_AFFINITY, D._ALL_STOPPED)
        @test D._built_reading(_testhelper_mirror(second_try_a), r)[:, 1] ==
              D._built_reading(second_try_a, r)[:, 2]
    end
    # The first reading alone finds a bridge at zero substrates in the inactive
    # conformation.
    v = D._admission(second_try_a, D._built_reading(second_try_a))
    @test (_testhelper_faces(second_try_a)[v.face], v.conformation, v.edge) ==
          (Set([:A, :B, :D]), 2, _testhelper_edge(second_try_a, "EDP + Q ⇌ EDPQ"))
    for am in (second_try_a, third_try)
        # The first try alone.
        @test D._admission(am, D._built_reading(am)).fault === D._BRIDGE
        @test D._admission(am, D._built_reading(am, D._ALL_STOPPED)) == D._WELL_POSED
        @test D._guard(am) == D._WELL_POSED
        @test all(_testhelper_admission_agreement(am; rng = Random.MersenneTwister(1)))
        # The compiled law, derived through the guard, against the oracle's own deletion
        # limit.
        _testhelper_check_law_against_limit(EnzymeRates.compile_mechanism(am))
    end
    # The outward reading admits second_try_a but not third_try, so each later try is
    # needed.
    @test D._admission(second_try_a,
                       D._built_reading(second_try_a, D._OUTWARD_LOSE_AFFINITY)) ==
          D._WELL_POSED
    @test D._admission(third_try,
                       D._built_reading(third_try, D._OUTWARD_LOSE_AFFINITY)).fault ===
          D._BRIDGE
    # bridge_fast_a: its bridge, the rapid-equilibrium P binding E(D, Q) + P ⇌ E(D, P, Q),
    # loses its affinity under every reading, and stopping its steady-state uptakes
    # leaves forms that free E does not reach. At zero P the ε-path oracle finds that the
    # limit depends on the direction: the rejection is real.
    bridge_fast_a =@_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + D ⇌ E(D), E(A, P) + D ⇌ E(A, D, P), E(P) + D ⇌ E(D, P),
             E(P, Q) + D ⇌ E(D, P, Q)) :: EqualAI
            (E + P ⇌ E(P), E(D) + P ⇌ E(D, P), E(D, Q) + P ⇌ E(D, P, Q),
             E(R; residual = A + D - R) + P ⇌ E(A, D, P)) :: OnlyA
            E(; residual = A + D - R) + B ⇌ E(P, Q) :: OnlyA
            E(; residual = A + D - R) + R ⇌ E(R; residual = A + D - R) :: EqualAI
            (E(A, P) + Q <--> E(A, P, Q), E(D) + Q <--> E(D, Q),
             E(D, P) + Q <--> E(D, P, Q), E(P) + Q <--> E(P, Q),
             E(R; residual = A + D - R) + Q <--> E(A, D, Q)) :: OnlyA
            (E(D) + A <--> E(R; residual = A + D - R), E(D, P) + A <--> E(A, D, P),
             E(D, Q) + A <--> E(A, D, Q), E(P) + A <--> E(A, P),
             E(P, Q) + A <--> E(A, P, Q)) :: EqualAI
        end
    end
    v = D._guard(bridge_fast_a)
    @test v == D._admission(bridge_fast_a, D._built_reading(bridge_fast_a))
    @test (v.fault, v.conformation, v.edge) ==
          (D._BRIDGE, 2, _testhelper_edge(bridge_fast_a, "EDQ + P ⇌ EDPQ"))
    @test D._admission(bridge_fast_a,
                       D._built_reading(bridge_fast_a, D._ALL_STOPPED)).fault === D._ISLAND
    @test D._admission(bridge_fast_a,
                       D._built_reading(bridge_fast_a, D._OUTWARD_LOSE_AFFINITY)).fault ===
          D._BRIDGE
    @test _testhelper_epath_verdict(bridge_fast_a; rng, face = Set([:P])) ===
          :direction_dependent
    # bridge_fast_b: every excluded uptake into a cut-off form is a rapid-equilibrium step,
    # so there is no steady-state one to read again.
    bridge_fast_b = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + A ⇌ E(A), E(B, Q) + A ⇌ E(A, B, Q),
             E(P; residual = A + B - P) + A ⇌ E(A, P; residual = A + B - P),
             E(Q) + A ⇌ E(A, Q)) :: OnlyA
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q),
             E(P; residual = A + B - P) + Q ⇌ E(A, B, Q)) :: OnlyA
            (E(; residual = A + B - P) + P ⇌ E(P; residual = A + B - P),
             E(A) + P ⇌ E(A, P),
             E(P; residual = A + B - P) + P ⇌ E(A, B, P)) :: EqualAI
            E(; residual = A + B - P) + D ⇌ E(Q, R) :: EqualAI
            (E(A) + B <--> E(P; residual = A + B - P), E(A, P) + B <--> E(A, B, P),
             E(A, Q) + B <--> E(A, B, Q),
             E(P; residual = A + B - P) + B <--> E(B, P; residual = A + B - P),
             E(Q) + B <--> E(B, Q)) :: EqualAI
            (E(Q) + R <--> E(Q, R),
             E(Q, R) + R <--> E(D, R; residual = A + B - P)) :: EqualAI
        end
    end
    @test D._guard(bridge_fast_b).fault === D._BRIDGE
    for c in 1:2
        @test D._built_reading(bridge_fast_b)[:, c] ==
              D._built_reading(bridge_fast_b, D._ALL_STOPPED)[:, c]
    end
    # The first try decides when it passes: with its S binding stopped, the K-type
    # uni-uni mechanism would leave E(S) out of free E's reach.
    k_type = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(S)     :: OnlyA
            E(S) <--> E(P)      :: OnlyA
            E(P) ⇌ E + P        :: EqualAI
        end
    end
    @test D._guard(k_type) == D._WELL_POSED
    @test D._admission_interior(k_type, D._built_reading(k_type, D._ALL_STOPPED)).fault ===
          D._ISLAND
end

@testset "a mechanism may need one uptake to lose its affinity and another to stop" begin
    # The guard's second try (`_OUTWARD_LOSE_AFFINITY`): an excluded steady-state step that
    # takes up one metabolite into a cut-off form loses its affinity only when it leads
    # away from free E, its bound side's part farther from free E's than its unbound
    # side's; otherwise it is stopped. A part's distance counts the excluded uptakes into
    # cut-off forms on the way to it from free E's part. The fixtures come from the gate's
    # populations.
    D = EnzymeRates.Derivation
    # both_readings, inactive conformation: the steps it keeps join E(A) to E(A, P) and
    # E(P) to E(P, Q), each part one uptake from free E's. The B binding, from one part to
    # the other, must stop, or the steps that lose their affinity close a cycle with the
    # A binding at E(P); the P binding at E must then lose its affinity, or no such step
    # leads to E(P) and E(P, Q).
    both_readings = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(P) + A ⇌ E(A, P))  :: OnlyA
            E + P <--> E(P)                     :: OnlyA
            E(A) + P <--> E(A, P)               :: EqualAI
            E(A) + B <--> E(P, Q)               :: OnlyA
            E(P) + Q ⇌ E(P, Q)                  :: EqualAI
        end
    end
    am = both_readings
    p = _testhelper_edge(am, "E + P <--> EP")              # leads away: distance 0 to 1
    b = _testhelper_edge(am, "EA + B <--> EPQ")            # distance 1 to 1
    a = _testhelper_edge(am, "EP + A ⇌ EAP")               # rapid equilibrium, 1 to 1
    L, S = D._LOSES_FORWARD, D._STOPPED
    for (reading, name, at_p, at_b) in (
            (D._ALL_LOSE_AFFINITY, :all_lose_affinity, L, L),
            (D._OUTWARD_LOSE_AFFINITY, :outward_lose_affinity, L, S),
            (D._ALL_STOPPED, :all_stopped, S, S))
        off = D._built_reading(am, reading)[:, 2]
        @test off[p] === at_p && off[b] === at_b
        @test off[a] === D._LOSES_FORWARD                     # a fast step is read one way
        @test _testhelper_reading(D._built_reading(am, reading))[:, 2] ==
              _testhelper_reading(am, name)[:, 2]
    end
    # The same need with the substrates' roles swapped and an :OnlyI P binding, which the
    # guard rejected before this rule.
    both_readings_only_i = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + B ⇌ E(B), E(P) + B ⇌ E(B, P))  :: OnlyA
            E + P <--> E(P)                     :: OnlyA
            E(B) + P <--> E(B, P)               :: OnlyI
            E(B) + A <--> E(P, Q)               :: OnlyA
            E(P) + Q ⇌ E(P, Q)                  :: EqualAI
        end
    end
    for am in (both_readings, both_readings_only_i)
        @test D._admission(am, D._built_reading(am)).fault === D._ARC_CYCLE
        @test D._admission(am, D._built_reading(am, D._ALL_STOPPED)).fault === D._ISLAND
        @test D._admission(am, D._built_reading(am, D._OUTWARD_LOSE_AFFINITY)) ==
              D._WELL_POSED
        @test D._guard(am) == D._WELL_POSED
        @test all(_testhelper_admission_agreement(am; rng = Random.MersenneTwister(1)))
        # The compiled law, derived through the guard, against the oracle's own deletion
        # limit.
        _testhelper_check_law_against_limit(EnzymeRates.compile_mechanism(am))
    end
    # outward_pumped: the first try admits it, and the outward reading alone would lose
    # it.
    outward_pumped = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            E + A ⇌ E(A)                                    :: OnlyA
            E + R <--> E(R)                                 :: OnlyA
            (E(A) + B <--> E(A, B), E(R) + B <--> E(B, R))  :: OnlyA
            (E(A) + P ⇌ E(A, P), E(R) + P ⇌ E(P, R))        :: EqualAI
            E(A, B) + D <--> E(P, R) + Q                    :: EqualAI
        end
    end
    @test D._admission(outward_pumped, D._built_reading(outward_pumped)) == D._WELL_POSED
    @test D._admission(outward_pumped,
                       D._built_reading(outward_pumped, D._OUTWARD_LOSE_AFFINITY)).fault ===
          D._PUMPED
    @test D._guard(outward_pumped) == D._WELL_POSED
    @test all(_testhelper_admission_agreement(outward_pumped;
                                              rng = Random.MersenneTwister(1)))
    # In the inactive conformation no uptake leads from free E to E(; residual = A - P),
    # which only the chemistry and the P release reach: its B binding leaves a part no
    # uptake reaches, so the outward reading stops it.
    unreached = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                            :: EqualAI
            E(A) <--> E(P; residual = A - P)                        :: OnlyA
            E(P; residual = A - P) ⇌ E(; residual = A - P) + P      :: OnlyA
            E(; residual = A - P) + B <--> E(B; residual = A - P)   :: OnlyA
            E(B; residual = A - P) <--> E(Q)                        :: OnlyA
            E(Q) ⇌ E + Q                                            :: OnlyA
        end
    end
    k = _testhelper_edge(unreached, "E_res_plusA_minusP + B <--> EB_res_plusA_minusP")
    @test D._built_reading(unreached)[k, 2] === D._LOSES_FORWARD
    @test D._built_reading(unreached, D._OUTWARD_LOSE_AFFINITY)[k, 2] === D._STOPPED
    # The package's reading agrees with the test's own statement of the rule, and the two
    # conformations are read alike.
    for x in (both_readings, both_readings_only_i, outward_pumped, unreached),
        (name, r) in pairs(_testhelper_guard_readings)
        @test _testhelper_reading(D._built_reading(x, r)) == _testhelper_reading(x, name)
        @test D._built_reading(_testhelper_mirror(x), r)[:, 1] ==
              D._built_reading(x, r)[:, 2]
    end
end

@testset "a reading may lose a step in either direction" begin
    # The ping-pong of the spec's fourth example, both chemistry steps :OnlyA. In the
    # inactive conformation the covalent forms F, F(B) and F(P) form one cut-off group
    # (F is E(; residual = A - P)). E(A) <--> F(P) loses its stored direction, into that
    # group, and F(B) ⇌ E(Q) its reverse, also into that group: both point uphill from
    # free E's group, and the cut-off group empties through both halves. No built reading
    # loses a step in reverse.
    D = EnzymeRates.Derivation
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                            :: EqualAI
            E + Q ⇌ E(Q)                                            :: EqualAI
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)      :: EqualAI
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)      :: EqualAI
            E(A) <--> E(P; residual = A - P)                        :: OnlyA
            E(B; residual = A - P) ⇌ E(Q)                           :: OnlyA
        end
    end
    chem_a = _testhelper_edge(am, "EA <--> EP_res_plusA_minusP")
    chem_b = _testhelper_edge(am, "EB_res_plusA_minusP ⇌ EQ")
    r = fill(D._RUNS, length(am.catalytic.steps), 2)
    r[chem_a, 2], r[chem_b, 2] = D._LOSES_FORWARD, D._LOSES_REVERSE
    @test D._admission_interior(am, r) == D._WELL_POSED
    eps = D._eps_network(am, D._re_flags(am.catalytic), r)
    fwd, rev = (only(t for t in eachindex(eps.from)
                     if eps.edge[t] == chem_b && eps.conformation[t] == 2 &&
                        eps.forward[t] == forward)
                for forward in (true, false))
    # F(B) ⇌ E(Q) in copy 2: the reverse (E(Q) → F(B)) vanishes, the forward does not,
    # both fast.
    @test [eps.vanishing[t] for t in (fwd, rev)] == [false, true]
    @test eps.fast[fwd] && eps.fast[rev]
    # Exchanged, both steps point into free E's group, below which nothing lies: no
    # energies exist.
    flipped = copy(r)
    flipped[chem_a, 2], flipped[chem_b, 2] = D._LOSES_REVERSE, D._LOSES_FORWARD
    # The verdict names the first of them in edge order.
    @test D._admission_interior(am, flipped) ==
          D._Verdict(D._ARC_CYCLE, D._LIMIT_UNTESTED, 0, 2, chem_a)
    # A step keeps its kind, so a reading that stops the rapid-equilibrium F(B) ⇌ E(Q) is
    # refused as one that cannot switch it off.
    stopped_fast = copy(r)
    stopped_fast[chem_b, 2] = D._STOPPED
    @test D._admission_interior(am, stopped_fast) ==
          D._Verdict(D._FAST_STEP_PINNED, D._LIMIT_UNTESTED, 0, 2, chem_b)
    stopped_names = _testhelper_reading(stopped_fast)
    @test !_testhelper_reading_is_consistent(am, stopped_names)
    @test _testhelper_exclusion(am, 2; reading = stopped_names).failure === :fast
    # A step the conformation keeps joins its two sides into one group, so a direction it
    # loses points uphill inside that group and closes a cycle alone.
    b_binding = _testhelper_edge(am, "E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP")
    inside = copy(r)
    inside[b_binding, 2] = D._LOSES_FORWARD
    @test D._admission_interior(am, inside) ==
          D._Verdict(D._ARC_CYCLE, D._LIMIT_UNTESTED, 0, 2, b_binding)
    # The tests' own names for the reading's entries.
    expected = fill(:runs, length(am.catalytic.steps), 2)
    expected[chem_a, 2], expected[chem_b, 2] = :loses_forward, :loses_reverse
    @test _testhelper_reading(r) == expected
    # The ε-path oracle, fed the same reading, finds the limit on every face; under the
    # guard's first built reading F(B) ⇌ E(Q) takes up nothing, so it cannot be switched
    # off, and there is no path.
    rng = Random.MersenneTwister(2614)
    for face in _testhelper_faces(am)
        @test _testhelper_epath_verdict(am; rng, face, corners = true,
                                        reading = _testhelper_reading(r)) === :converges
    end
    @test _testhelper_epath_verdict(
        am; rng, reading = _testhelper_reading(am, :all_lose_affinity)) === :no_path
end

@testset "no form fills faster than it empties" begin
    D = EnzymeRates.Derivation
    NO, FI, U = D._NO_FAULT, D._FILLED, D._LIMIT_UNTESTED
    rng = Random.MersenneTwister(2650)
    # `filled`: a ping-pong whose B binding at the covalent form F (E(; residual = A - P))
    # and first chemistry step are :OnlyA. In the inactive conformation F and F(P) form a
    # cut-off group: E(A) <--> F(P) loses its stored direction, into it, and F + B ⇌ F(B)
    # its reverse, the release of B, also into it. With B zero the uptake of B dies, and
    # the release, fast though it vanishes, runs one way from the law's form F(B) into F,
    # whose group empties only through the slow chemistry F(P) → E(A): F fills faster
    # than it empties.
    filled = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                            :: EqualAI
            E + Q ⇌ E(Q)                                            :: EqualAI
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)      :: OnlyA
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)      :: EqualAI
            E(A) <--> E(P; residual = A - P)                        :: OnlyA
            E(B; residual = A - P) ⇌ E(Q)                           :: EqualAI
        end
    end
    chem_a = _testhelper_edge(filled, "EA <--> EP_res_plusA_minusP")
    b_binding = _testhelper_edge(filled, "E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP")
    r = fill(D._RUNS, length(filled.catalytic.steps), 2)
    r[chem_a, 2], r[b_binding, 2] = D._LOSES_FORWARD, D._LOSES_REVERSE
    faces = _testhelper_faces(filled)
    # The first face with B zero zeroes A as well; B alone is refused the same way.
    v = D._admission(filled, r)
    @test v == D._Verdict(FI, U, 3, 2, b_binding)
    @test faces[3] == Set([:A, :B])
    @test faces[5] == Set([:B])
    @test _testhelper_admission_faults(filled; reading = _testhelper_reading(r)) ==
          [NO, NO, FI, NO, FI, NO, NO]
    # The message names the face, the step, the form it fills and the group's slow way
    # out, and offers the fix the guard admits: F + B <--> F(B) steady-state.
    @test D._steady_state_passes(filled, b_binding)
    @test D._admission_message(filled, v, r) ==
          "AllostericMechanism: at zero A and B, in the inactive conformation the :OnlyA " *
          "rapid-equilibrium step E(; residual = A - P) + B ⇌ E(B; residual = A - P) " *
          "runs one way into E(; residual = A - P), at a rate that vanishes yet grows " *
          "with the rapid-equilibrium speed, and the forms E(; residual = A - P), " *
          "E(P; residual = A - P) leave only through the slower step " *
          "E(A) <--> E(P; residual = A - P), so E(; residual = A - P) fills faster than " *
          "it empties and its weight depends on how fast the tagged step vanishes. Make " *
          "the group of E(; residual = A - P) + B ⇌ E(B; residual = A - P) steady-state."
    filled_ss = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                            :: EqualAI
            E + Q ⇌ E(Q)                                            :: EqualAI
            E(; residual = A - P) + B <--> E(B; residual = A - P)   :: OnlyA
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)      :: EqualAI
            E(A) <--> E(P; residual = A - P)                        :: OnlyA
            E(B; residual = A - P) ⇌ E(Q)                           :: EqualAI
        end
    end
    @test D._guard(filled_ss) == D._WELL_POSED
    @test all(_testhelper_admission_agreement(filled_ss; rng = Random.MersenneTwister(2653),
                                              corners = true))
    _testhelper_check_law_against_limit(EnzymeRates.compile_mechanism(filled_ss))
    # The harsh ε-path oracle under the same reading: the limit depends on the direction
    # on both faces with B zero, and exists on the interior.
    for i in (3, 5)
        @test _testhelper_epath_verdict(filled; rng, face = faces[i], corners = true,
                                        reading = _testhelper_reading(r)) ===
              :direction_dependent
    end
    @test _testhelper_epath_verdict(filled; rng, corners = true,
                                    reading = _testhelper_reading(r)) === :converges

    # The six-step ping-pong with both chemistry steps :OnlyA, under the reading that
    # admits it: neither excluded step takes up a metabolite, so no face kills the
    # direction it keeps, and nothing fills one way.
    pingpong = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                            :: EqualAI
            E + Q ⇌ E(Q)                                            :: EqualAI
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)      :: EqualAI
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)      :: EqualAI
            E(A) <--> E(P; residual = A - P)                        :: OnlyA
            E(B; residual = A - P) ⇌ E(Q)                           :: OnlyA
        end
    end
    p = fill(D._RUNS, length(pingpong.catalytic.steps), 2)
    p[_testhelper_edge(pingpong, "EA <--> EP_res_plusA_minusP"), 2] = D._LOSES_FORWARD
    p[_testhelper_edge(pingpong, "EB_res_plusA_minusP ⇌ EQ"), 2] = D._LOSES_REVERSE
    @test D._admission(pingpong, p) == D._WELL_POSED

    # `slow_source`: a ping-pong whose B binding at F gives E(Q) at once, every step
    # :OnlyA but F + P ⇌ F(P). In the inactive conformation free E stands alone and the
    # steps that leave it are steady-state uptakes, which lose their uptake. With B zero
    # F + B ⇌ E(Q) runs one way, fast though it vanishes, from E(Q) into F; but enzyme
    # reaches E(Q) only through steady-state steps, so the flow into F vanishes, and the
    # faces with B zero pass.
    slow_source = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> E(P; residual = A - P), E(Q) + A <--> E(A, Q))    :: OnlyA
            (E + Q <--> E(Q), E(P; residual = A - P) + Q <--> E(A, Q))    :: OnlyA
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)            :: EqualAI
            E(; residual = A - P) + B ⇌ E(Q)                              :: OnlyA
        end
    end
    b_step = _testhelper_edge(slow_source, "E_res_plusA_minusP + B ⇌ EQ")
    s = fill(D._RUNS, length(slow_source.catalytic.steps), 2)
    for text in ("E + A <--> EP_res_plusA_minusP", "E + Q <--> EQ",
                 "EP_res_plusA_minusP + Q <--> EAQ", "EQ + A <--> EAQ")
        s[_testhelper_edge(slow_source, text), 2] = D._LOSES_FORWARD
    end
    s[b_step, 2] = D._LOSES_REVERSE
    @test D._admission(slow_source, s) == D._WELL_POSED
    eps = D._eps_network(slow_source, D._re_flags(slow_source.catalytic), s)
    law_nodes = D._law_network(slow_source).present
    sets = D._face_sets(eps, law_nodes, D._faces(slow_source.catalytic)[3])
    fast_from_law = D._reached(D._successors(eps, sets.alive .& eps.fast),
                               findall(law_nodes))
    uptake, release = (only(t for t in eachindex(eps.from)
                            if eps.edge[t] == b_step && eps.conformation[t] == 2 &&
                               eps.forward[t] == forward)
                       for forward in (true, false))
    # At zero A and B, the release E(Q) → F + B is fast, vanishing and alive while the
    # uptake is dead; its source is fed, but not along fast transitions from the law's
    # forms.
    @test eps.fast[release] && eps.vanishing[release] && !eps.vanishing[uptake]
    @test sets.alive[release] && !sets.alive[uptake]
    @test sets.fed[eps.from[release]] && !fast_from_law[eps.from[release]]
    for face in _testhelper_faces(slow_source)
        @test _testhelper_epath_verdict(slow_source; rng, face, corners = true,
                                        reading = _testhelper_reading(s)) === :converges
    end

    # No built reading raises the fault. Under each of them `filled`'s excluded steps are
    # stopped or cannot be switched off, so none keeps a direction; `slow_source`'s
    # F + B ⇌ E(Q) loses its affinity, and the direction it keeps is the release of B,
    # which takes up nothing, so no face kills it.
    for am in (filled, slow_source)
        law, built_faces = D._law_network(am), D._faces(am.catalytic)
        for b in instances(D._SlowUptakeReading)
            eps = D._eps_network(am, D._re_flags(am.catalytic), D._built_reading(am, b))
            @test all(i -> D._admission_face(am, eps, law, built_faces, i).fault !== FI,
                      eachindex(built_faces))
        end
    end
    @test D._built_reading(slow_source)[b_step, 2] === D._LOSES_FORWARD
end

@testset "cut-off forms that cannot matter are set aside" begin
    D = EnzymeRates.Derivation
    NO, BR, U = D._NO_FAULT, D._BRIDGE, D._LIMIT_UNTESTED
    rng = Random.MersenneTwister(2651)
    # Whether face `i` of `am` under reading `r` is refused by the face condition alone and
    # passed by admission's face check, which sets groups aside.
    function _testhelper_passes_only_set_aside(am, r, i)
        eps = D._eps_network(am, D._re_flags(am.catalytic), r)
        law, faces = D._law_network(am), D._faces(am.catalytic)
        v = D._face_condition(eps, law.present, faces[i], D._turnover_name(am.catalytic))
        v.fault === BR && D._admission_face(am, eps, law, faces, i) == D._WELL_POSED
    end
    # `set_aside`: a ping-pong whose chemistry steps and the bindings of A and of P are
    # :OnlyA. In the inactive conformation the law's forms are free E and E(Q), and the
    # cut-off forms E(A), F(P), F and F(B) (F is E(; residual = A - P)), which the excluded
    # steps join, are entered from free E by the A binding and from E(Q) by
    # F(B) ⇌ E(Q), both losing their direction into them. Their vanishing steps open a
    # route from free E to E(Q) on which the reaction turns over, so the face condition
    # refuses the interior as a bridge. But the law's own rapid-equilibrium step
    # E + Q ⇌ E(Q) holds the two forms they touch in equilibrium, their fast classes drain,
    # and none fills faster than it empties: they are set aside.
    set_aside = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                            :: OnlyA
            E + Q ⇌ E(Q)                                            :: EqualAI
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)      :: EqualAI
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)      :: OnlyA
            E(A) <--> E(P; residual = A - P)                        :: OnlyA
            E(B; residual = A - P) ⇌ E(Q)                           :: OnlyA
        end
    end
    r = fill(D._RUNS, length(set_aside.catalytic.steps), 2)
    r[_testhelper_edge(set_aside, "E + A ⇌ EA"), 2] = D._LOSES_FORWARD
    r[_testhelper_edge(set_aside, "EA <--> EP_res_plusA_minusP"), 2] = D._LOSES_FORWARD
    r[_testhelper_edge(set_aside, "E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP"), 2] =
        D._LOSES_REVERSE
    r[_testhelper_edge(set_aside, "EB_res_plusA_minusP ⇌ EQ"), 2] = D._LOSES_REVERSE
    @test D._admission(set_aside, r) == D._WELL_POSED
    @test _testhelper_admission_faults(set_aside; reading = _testhelper_reading(r)) ==
          fill(NO, 7)
    # The face condition alone refuses the interior and the face with A zero. There the
    # uptake E + A → E(A) is dead, so that pair runs one way, through its release, which
    # does not vanish and only empties E(A): it is safe, and the cut-off forms go again.
    @test _testhelper_faces(set_aside)[4] == Set([:A])
    let eps = D._eps_network(set_aside, D._re_flags(set_aside.catalytic), r),
        law = D._law_network(set_aside), faces = D._faces(set_aside.catalytic)
        x = D._turnover_name(set_aside.catalytic)
        @test [D._face_condition(eps, law.present, f, x).fault for f in faces] ==
              [BR, NO, NO, BR, NO, NO, NO]
    end
    @test _testhelper_passes_only_set_aside(set_aside, r, 1)
    @test _testhelper_passes_only_set_aside(set_aside, r, 4)
    for face in _testhelper_faces(set_aside)
        @test _testhelper_epath_verdict(set_aside; rng, face, corners = true,
                                        reading = _testhelper_reading(r)) === :converges
    end
    # `kept`: the same ping-pong with E + Q <--> E(Q) steady-state, so free E and E(Q) are
    # no longer held in equilibrium, and E(A) ⇌ F(P) at rapid equilibrium, so the route
    # through the cut-off forms runs fast, its rates vanishing only as the reading says.
    # In the interior every transition at the cut-off forms runs both ways, and each of
    # them empties into the law's forms along fast steps that do not vanish, so none lies
    # in a fast class that could fail to drain; but the route joins two law forms whose
    # ratio the law sets otherwise: nothing is set aside, the interior stays a bridge, and
    # the limit there depends on how the two limits are taken.
    kept = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                            :: OnlyA
            E + Q <--> E(Q)                                         :: EqualAI
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)      :: EqualAI
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)      :: OnlyA
            E(A) ⇌ E(P; residual = A - P)                           :: OnlyA
            E(B; residual = A - P) ⇌ E(Q)                           :: OnlyA
        end
    end
    k = fill(D._RUNS, length(kept.catalytic.steps), 2)
    k[_testhelper_edge(kept, "E + A ⇌ EA"), 2] = D._LOSES_FORWARD
    k[_testhelper_edge(kept, "EA ⇌ EP_res_plusA_minusP"), 2] = D._LOSES_FORWARD
    k[_testhelper_edge(kept, "E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP"), 2] =
        D._LOSES_REVERSE
    k[_testhelper_edge(kept, "EB_res_plusA_minusP ⇌ EQ"), 2] = D._LOSES_REVERSE
    eps = D._eps_network(kept, D._re_flags(kept.catalytic), k)
    law = D._law_network(kept)
    sets = D._face_sets(eps, law.present, UInt64(0))
    @test all(D._fast_law_reach(eps, law.present, sets).to_law[.!law.present])
    @test D._without_set_aside_groups(eps, law.present, UInt64(0)) === eps
    @test D._admission(kept, k) ==
          D._Verdict(BR, U, 1, 2, _testhelper_edge(kept, "EB_res_plusA_minusP ⇌ EQ"))
    @test _testhelper_epath_verdict(kept; rng, corners = true,
                                    reading = _testhelper_reading(k)) ===
          :direction_dependent

    # One touching law form counts: in this ping-pong's inactive conformation free E stands
    # alone, and the cut-off forms E(A), F(P), F, F(B) and E(Q) touch nothing else.
    one_touch = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                            :: OnlyA
            E + Q ⇌ E(Q)                                            :: OnlyA
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)      :: EqualAI
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)      :: OnlyA
            E(A) <--> E(P; residual = A - P)                        :: EqualAI
            E(B; residual = A - P) ⇌ E(Q)                           :: EqualAI
        end
    end
    o = fill(D._RUNS, length(one_touch.catalytic.steps), 2)
    o[_testhelper_edge(one_touch, "E + A ⇌ EA"), 2] = D._LOSES_FORWARD
    o[_testhelper_edge(one_touch, "E + Q ⇌ EQ"), 2] = D._LOSES_FORWARD
    o[_testhelper_edge(one_touch, "E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP"), 2] =
        D._LOSES_REVERSE
    @test D._admission(one_touch, o) == D._WELL_POSED
    @test _testhelper_admission_faults(one_touch; reading = _testhelper_reading(o)) ==
          fill(NO, 7)
    # The face condition alone refuses five faces. With P zero the uptake F + P → F(P) is
    # dead, and F + P ⇌ F(P) runs one way between two cut-off forms, through its release,
    # which vanishes and is fast; it is safe because no fast transition from the law's
    # forms reaches its source F(P), so the flow through it vanishes.
    @test _testhelper_faces(one_touch)[6] == Set([:P])
    let eps = D._eps_network(one_touch, D._re_flags(one_touch.catalytic), o),
        law = D._law_network(one_touch), faces = D._faces(one_touch.catalytic)
        x = D._turnover_name(one_touch.catalytic)
        @test [D._face_condition(eps, law.present, f, x).fault for f in faces] ==
              [BR, BR, NO, BR, NO, BR, BR]
    end
    @test _testhelper_passes_only_set_aside(one_touch, o, 1)
    @test _testhelper_passes_only_set_aside(one_touch, o, 6)
    for face in _testhelper_faces(one_touch)
        @test _testhelper_epath_verdict(one_touch; rng, face, corners = true,
                                        reading = _testhelper_reading(o)) === :converges
    end

    # `first_try`, a ter-ter mechanism of the gate's populations: under the guard's first
    # reading the face condition alone refuses five faces as bridges of the inactive
    # conformation, and admission passes every one by setting cut-off forms aside, so the
    # first try admits it. The oracles agree on every face under that reading, and the
    # compiled law equals the oracle's deletion limit.
    first_try = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + A ⇌ E(A), E(D, P) + A ⇌ E(A, D, P), E(P) + A ⇌ E(A, P),
             E(P, R) + A ⇌ E(A, P, R)) :: OnlyA
            (E + P ⇌ E(P), E(A) + P ⇌ E(A, P), E(A, D) + P ⇌ E(A, D, P),
             E(A, R) + P ⇌ E(A, P, R)) :: OnlyA
            (E(A) + D ⇌ E(A, D), E(A, P) + D ⇌ E(A, D, P), E(A, R) + D ⇌ E(A, D, R),
             E(P) + D ⇌ E(D, P), E(P, R) + D ⇌ E(D, P, R)) :: EqualAI
            (E(A) + R <--> E(A, R), E(A, D) + R <--> E(A, D, R),
             E(A, P) + R <--> E(A, P, R), E(D, P) + R <--> E(D, P, R),
             E(P) + R <--> E(P, R)) :: EqualAI
            E(A, D) + B <--> E(P, Q, R) :: OnlyA
            E(P, R) + Q ⇌ E(P, Q, R) :: OnlyA
        end
    end
    let reading = D._built_reading(first_try)
        eps = D._eps_network(first_try, D._re_flags(first_try.catalytic), reading)
        law, faces = D._law_network(first_try), D._faces(first_try.catalytic)
        x = D._turnover_name(first_try.catalytic)
        @test [D._face_condition(eps, law.present, f, x).fault for f in faces] ==
              [BR, NO, BR, BR, BR, BR, NO, NO, NO]
        @test D._admission(first_try, reading) == D._WELL_POSED
    end
    @test D._guard(first_try) == D._WELL_POSED
    @test all(_testhelper_admission_agreement(first_try; rng = Random.MersenneTwister(1)))
    _testhelper_check_law_against_limit(EnzymeRates.compile_mechanism(first_try))

    # R2funnel and R2funnel_SS under the guard's first built reading, at zero P. Their
    # cut-off forms touch free E alone, and each reaches it along transitions that do not
    # vanish; but F and F(S) form a fast class whose one slow way out, to F(P), leads back
    # into it along fast steps, while the slow E(P) <--> E(S) leaves E(P), which a fast
    # step empties into F. The class does not drain, nothing is set aside, the face stays
    # a bridge, and the limit there depends on the direction.
    for name in ("R2funnel", "R2funnel_SS")
        am = _testhelper_tag_sets[name]
        @test _testhelper_faces(am)[2] == Set([:P])
        eps = D._eps_network(am, D._re_flags(am.catalytic), D._built_reading(am))
        law, faces = D._law_network(am), D._faces(am.catalytic)
        @test !any(D._face_sets(eps, law.present, faces[2]).stranded)
        @test D._without_set_aside_groups(eps, law.present, faces[2]) === eps
        @test D._admission_face(am, eps, law, faces, 2).fault === BR
        @test _testhelper_epath_verdict(am; rng, face = _testhelper_faces(am)[2],
                                        corners = true) === :direction_dependent
    end

    # The safety clause deciding. `stopped_exit` (a gate entry; the mechanism of
    # `slow_source` in "no form fills faster than it empties") under a reading that stops
    # E + A <--> F(P), F being E(; residual = A - P), rather than letting it lose its
    # uptake. In the inactive conformation free E stands alone, so the cut-off forms
    # touch one law form, and their fast classes drain through F(P) + Q <--> E(A, Q),
    # E(A, Q) → E(Q) + A and E(Q) → E + Q. With A and B zero the stopped step runs one
    # way, from F(P) into free E, along its vanishing direction: a pair at a law form with
    # only its vanishing direction alive, which is not safe. It alone keeps them: nothing
    # is set aside, the face stays a bridge, and the limit there depends on the direction.
    stopped_exit = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A <--> E(P; residual = A - P), E(Q) + A <--> E(A, Q))    :: OnlyA
            (E + Q <--> E(Q), E(P; residual = A - P) + Q <--> E(A, Q))    :: OnlyA
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)            :: EqualAI
            E(; residual = A - P) + B ⇌ E(Q)                              :: OnlyA
        end
    end
    se = fill(D._RUNS, length(stopped_exit.catalytic.steps), 2)
    se[_testhelper_edge(stopped_exit, "E + A <--> EP_res_plusA_minusP"), 2] = D._STOPPED
    for text in ("E + Q <--> EQ", "EQ + A <--> EAQ")
        se[_testhelper_edge(stopped_exit, text), 2] = D._LOSES_FORWARD
    end
    for text in ("E_res_plusA_minusP + B ⇌ EQ", "EP_res_plusA_minusP + Q <--> EAQ")
        se[_testhelper_edge(stopped_exit, text), 2] = D._LOSES_REVERSE
    end
    @test _testhelper_faces(stopped_exit)[3] == Set([:A, :B])
    let eps = D._eps_network(stopped_exit, D._re_flags(stopped_exit.catalytic), se),
        law = D._law_network(stopped_exit), faces = D._faces(stopped_exit.catalytic)
        @test count(law.present[(length(stopped_exit.catalytic.nodes) + 1):end]) == 1
        @test !any(D._face_sets(eps, law.present, faces[3]).stranded)
        @test D._without_set_aside_groups(eps, law.present, faces[3]) === eps
    end
    @test D._admission(stopped_exit, se) ==
          D._Verdict(BR, U, 3, 2, _testhelper_edge(stopped_exit, "E + Q <--> EQ"))
    @test _testhelper_epath_verdict(stopped_exit; rng,
                                    face = _testhelper_faces(stopped_exit)[3],
                                    corners = true,
                                    reading = _testhelper_reading(se)) ===
          :direction_dependent

    # The conservative rejections the suite pins stay rejected under the guard's first
    # built reading. R-B1's cut-off forms E(I) and E(I, J) touch free E and the dead end
    # E(J), which only the steady-state E + J <--> E(J) joins; R-B2's basins reach the
    # law's forms only through vanishing steps, so they are stranded. (The search admits
    # R-B2 under other readings: "conservative rejections the oracle admits stay rejected
    # under their reading".)
    for (name, fault) in (("bridge_deadend", BR), ("BS1", D._PUMPED),
                          ("fz101_10232", D._PUMPED))
        am = _testhelper_tag_sets[name]
        @test (name, D._admission(am, D._built_reading(am)).fault) == (name, fault)
    end
    @test D._guard(_testhelper_tag_sets["bridge_deadend"]).fault === BR
end

@testset "the guard searches for a reading" begin
    D = EnzymeRates.Derivation
    rng = Random.MersenneTwister(2652)
    # The six-step ping-pong of the spec's fourth example, both chemistry steps :OnlyA (F is
    # E(; residual = A - P)). No built reading switches off F(B) ⇌ E(Q), a
    # rapid-equilibrium step that takes up nothing, so the guard searches. The first
    # reading it finds loses E(A) <--> F(P)'s stored direction and F(B) ⇌ E(Q)'s reverse,
    # both into the covalent group: the reading of "a reading may lose a step in either
    # direction".
    pingpong = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                            :: EqualAI
            E + Q ⇌ E(Q)                                            :: EqualAI
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)      :: EqualAI
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)      :: EqualAI
            E(A) <--> E(P; residual = A - P)                        :: OnlyA
            E(B; residual = A - P) ⇌ E(Q)                           :: OnlyA
        end
    end
    pingpong_dimer = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        catalytic_multiplicity: 2
        steps: begin
            E + A ⇌ E(A)                                            :: EqualAI
            E + Q ⇌ E(Q)                                            :: EqualAI
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)      :: EqualAI
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)      :: EqualAI
            E(A) <--> E(P; residual = A - P)                        :: OnlyA
            E(B; residual = A - P) ⇌ E(Q)                           :: OnlyA
        end
    end
    chem_a = _testhelper_edge(pingpong, "EA <--> EP_res_plusA_minusP")
    chem_b = _testhelper_edge(pingpong, "EB_res_plusA_minusP ⇌ EQ")
    r = fill(D._RUNS, length(pingpong.catalytic.steps), 2)
    r[chem_a, 2], r[chem_b, 2] = D._LOSES_FORWARD, D._LOSES_REVERSE
    @test D._built_reading(pingpong)[chem_b, 2] === D._CANNOT_SWITCH_OFF
    for am in (pingpong, pingpong_dimer)
        @test D._guard(am) == D._WELL_POSED
        @test D._admitting_reading(am) == r
        _testhelper_check_law_against_limit(EnzymeRates.compile_mechanism(am))
        @test all(_testhelper_admission_agreement(am; rng, corners = true))
    end
    # Its mirror image tags both chemistry steps :OnlyI, so the active conformation keeps
    # no chemistry: no reading helps.
    @test D._guard(_testhelper_mirror(pingpong)).fault === D._NO_ACTIVE_TURNOVER
    # The verdict does not depend on how the steps are written: the macro stores them in
    # canonical order. Nor on the order in which they are stored, which only orders the
    # search: stored backward, the same reading admits the mechanism.
    pingpong_rewritten = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E(B; residual = A - P) ⇌ E(Q)                           :: OnlyA
            E(A) <--> E(P; residual = A - P)                        :: OnlyA
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)      :: EqualAI
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)      :: EqualAI
            E + Q ⇌ E(Q)                                            :: EqualAI
            E + A ⇌ E(A)                                            :: EqualAI
        end
    end
    @test pingpong_rewritten == pingpong
    @test D._guard(pingpong_rewritten) == D._WELL_POSED
    let cm = pingpong.catalytic, back = reverse(eachindex(cm.steps))
        backward = EnzymeRates.AllostericMechanism(
            EnzymeRates.Mechanisms._unchecked_mechanism(cm.substrates, cm.products,
                                                        cm.nodes, cm.steps[back]),
            pingpong.mwc_tags, 1, pingpong.regulatory_sites)
        @test D._guard(backward) == D._WELL_POSED
        @test D._admitting_reading(backward) == r[back, :]
    end

    # The three-step ping-pong: the inactive conformation keeps only E + Q <--> E(Q). The
    # built readings stop both chemistry steps, which leaves F out of free E's reach; the
    # first reading the search finds loses the first half's stored direction and the
    # second half's reverse, both into F.
    three_step = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> F + P        :: OnlyA
            E + Q <--> E(Q)         :: EqualAI
            F + B <--> E(Q)         :: OnlyA
        end
    end
    half_a = _testhelper_edge(three_step, "E + A <--> F + P")
    half_b = _testhelper_edge(three_step, "F + B <--> EQ")
    t = fill(D._RUNS, length(three_step.catalytic.steps), 2)
    t[half_a, 2], t[half_b, 2] = D._LOSES_FORWARD, D._LOSES_REVERSE
    @test D._admission_interior(three_step, D._built_reading(three_step)).fault ===
          D._ISLAND
    @test D._guard(three_step) == D._WELL_POSED
    @test D._admitting_reading(three_step) == t
    _testhelper_check_law_against_limit(EnzymeRates.compile_mechanism(three_step))
    @test all(_testhelper_admission_agreement(three_step; rng, corners = true))

    # A ping-pong with no limit under any reading: the covalent forms empty only by
    # binding B. The guard finds no reading, and the ε-path oracle converges on every face
    # under none of the consistent readings.
    no_limit = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(P; residual = A - P)                       :: OnlyA
            E + Q <--> E(Q)                                         :: EqualAI
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)      :: EqualAI
            E(; residual = A - P) + B ⇌ E(Q)                        :: OnlyA
        end
    end
    @test D._admitting_reading(no_limit) === nothing
    @test D._guard(no_limit).fault !== D._NO_FAULT
    consistent = _testhelper_consistent_readings(no_limit)
    @test length(consistent) == 2
    for reading in consistent
        @test any(face -> _testhelper_epath_verdict(no_limit; rng, face, corners = true,
                                                    reading) !== :converges,
                  _testhelper_faces(no_limit))
    end

    # `filled` of "no form fills faster than it empties": its other readings fail too.
    filled = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                            :: EqualAI
            E + Q ⇌ E(Q)                                            :: EqualAI
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)      :: OnlyA
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)      :: EqualAI
            E(A) <--> E(P; residual = A - P)                        :: OnlyA
            E(B; residual = A - P) ⇌ E(Q)                           :: EqualAI
        end
    end
    @test D._admitting_reading(filled) === nothing
    @test D._guard(filled).fault !== D._NO_FAULT

    # A mechanism a built reading admits is admitted with it: `both_readings` of "a
    # mechanism may need one uptake to lose its affinity and another to stop".
    both_readings = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(P) + A ⇌ E(A, P))  :: OnlyA
            E + P <--> E(P)                     :: OnlyA
            E(A) + P <--> E(A, P)               :: EqualAI
            E(A) + B <--> E(P, Q)               :: OnlyA
            E(P) + Q ⇌ E(P, Q)                  :: EqualAI
        end
    end
    @test D._admitting_reading(both_readings) ==
          D._built_reading(both_readings, D._OUTWARD_LOSE_AFFINITY)

    # A rapid-equilibrium step inside one group is the fault no reading changes: the spec's
    # second example. The message offers the two fixes, each of whose variants passes.
    pinned = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: OnlyA
            E(S) <--> E(P)      :: EqualAI
            E + P ⇌ E(P)        :: EqualAI
        end
    end
    v = D._guard(pinned)
    @test v == D._Verdict(D._FAST_STEP_PINNED, D._LIMIT_UNTESTED, 0, 2,
                          _testhelper_edge(pinned, "E + S ⇌ ES"))
    @test D._guard_message(pinned, v) ==
          "AllostericMechanism: in the inactive conformation the rapid-equilibrium step " *
          "E + S ⇌ E(S) is tagged :OnlyA, but the steps that conformation keeps still " *
          "join its two sides. Thermodynamics then fixes its equilibrium constant, and a " *
          "step keeps its kind in both conformations, so it cannot be switched off. Make " *
          "its group steady-state, or tag :OnlyA the steps that join its sides: " *
          "E + P ⇌ E(P); E(S) <--> E(P)."
    @test _testhelper_consistent_readings(pinned) == Matrix{Symbol}[]

    # A cycle of rapid-equilibrium steps that turns the reaction over, each conformation
    # excluding two of them between its two groups. The built readings cannot switch off
    # the isomerizations, so under them the cycle is slow in each conformation; every
    # consistent reading lets each excluded step lose one direction, so the cycle runs
    # fast in both, and class 2 rejects it under the first, naming how each excluded step
    # of the cycle is switched off. The ε-path oracle finds the rate diverging.
    fast_cycle = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: OnlyI
            E(S) ⇌ Estar(S)     :: OnlyA
            Estar(S) ⇌ E(P)     :: OnlyI
            E + P ⇌ E(P)        :: OnlyA
        end
    end
    @test D._verdict_under(fast_cycle, D._built_reading(fast_cycle)).fault ===
          D._FAST_STEP_PINNED
    v = D._guard(fast_cycle)
    @test v.fault === D._RE_TURNOVER
    reported = _testhelper_reading(D._guard_with_reading(fast_cycle).reading)
    @test _testhelper_reading_is_consistent(fast_cycle, reported)
    msg = D._guard_message(fast_cycle, v)
    @test occursin("(whose affinity vanishes in the ", msg)
    @test occursin("(whose rate one way vanishes in the ", msg)
    @test _testhelper_epath_verdict(fast_cycle; rng, reading = reported) === :diverges

    # A mechanism and its mirror image, every :OnlyA and :OnlyI swapped, get one verdict
    # when both conformations catalyze: the six-step ping-pong with a ternary-complex path
    # beside it, which both conformations keep. The search admits both, the mirror with
    # the same reading in the other conformation, and the oracles agree on each.
    pp_ternary = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                            :: EqualAI
            E + Q ⇌ E(Q)                                            :: EqualAI
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)      :: EqualAI
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)      :: EqualAI
            E(A) <--> E(P; residual = A - P)                        :: OnlyA
            E(B; residual = A - P) ⇌ E(Q)                           :: OnlyA
            E(A) + B ⇌ E(A, B)                                      :: EqualAI
            E(A, B) <--> E(P, Q)                                    :: EqualAI
            E(Q) + P ⇌ E(P, Q)                                      :: EqualAI
        end
    end
    mirror = _testhelper_mirror(pp_ternary)
    @test D._admission_interior(pp_ternary, D._built_reading(pp_ternary)).fault ===
          D._FAST_STEP_PINNED
    p = fill(D._RUNS, length(pp_ternary.catalytic.steps), 2)
    p[_testhelper_edge(pp_ternary, "EA <--> EP_res_plusA_minusP"), 2] = D._LOSES_FORWARD
    p[_testhelper_edge(pp_ternary, "EB_res_plusA_minusP ⇌ EQ"), 2] = D._LOSES_REVERSE
    @test D._admitting_reading(pp_ternary) == p
    @test D._admitting_reading(mirror) == p[:, [2, 1]]
    for am in (pp_ternary, mirror)
        @test D._guard(am) == D._WELL_POSED
        _testhelper_check_law_against_limit(EnzymeRates.compile_mechanism(am))
        @test all(_testhelper_admission_agreement(am; rng, corners = true))
    end

    # The search is skipped whole for a mechanism too large to search: three ter-ter
    # taggings of the gate's populations whose inactive conformation has more ways to
    # switch off its excluded steps than the guard searches, counted by the tests' own
    # rule. The built readings admit `large_admitted`. They reject `large_island` and
    # `large_pumped`, and so does the guard, though a reading it never tries passes
    # every check: under no built reading are `large_island`'s forms all reached, so the
    # guard reports that there are too many readings; `large_pumped`'s first built reading
    # is consistent, so the guard reports its verdict.
    # The ways conformation `c` of `am` may switch off its excluded steps between two
    # groups, by the tests' own rule: each loses toward either side but free E's group,
    # or, at steady state, stops.
    function _testhelper_choice_count(am, c = 2)
        group = _testhelper_kept_components(am, c)
        prod((group[e.to] != group[1]) + (group[e.from] != group[1]) + !e.rapid_equilibrium
             for (k, e) in enumerate(am.catalytic.steps)
             if _testhelper_excludes(am, k, c) && group[e.from] != group[e.to]; init = 1)
    end
    large_admitted = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            E + D ⇌ E(D) :: EqualAI
            E + R ⇌ E(R) :: OnlyA
            (E(A, R) + P <--> E(A, P, R), E(Q, R) + P <--> E(P, Q, R),
             E(R) + P <--> E(P, R)) :: OnlyA
            (E(B, D) + A <--> E(P, Q, R), E(P, R) + A <--> E(A, P, R),
             E(R) + A <--> E(A, R)) :: OnlyA
            (E(B, R) + Q ⇌ E(B, Q, R), E(P, R) + Q ⇌ E(P, Q, R),
             E(R) + Q ⇌ E(Q, R)) :: OnlyA
            (E(D) + B ⇌ E(B, D), E(Q, R) + B ⇌ E(B, Q, R), E(R) + B ⇌ E(B, R)) :: OnlyA
        end
    end
    large_island = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + B <--> E(B), E(D) + B <--> E(Q; residual = B + D - Q)) :: OnlyA
            (E + D <--> E(D), E(B) + D <--> E(Q; residual = B + D - Q),
             E(P) + D <--> E(D, P),
             E(Q; residual = B + D - Q) + D <--> E(D, Q; residual = B + D - Q)) :: OnlyA
            (E + P <--> E(P), E(D) + P <--> E(D, P), E(R) + P <--> E(P, R)) :: OnlyA
            (E + R <--> E(R), E(P) + R <--> E(P, R)) :: EqualAI
            E(; residual = B + D - Q) + A ⇌ E(P, R) :: OnlyA
            (E(; residual = B + D - Q) + Q ⇌ E(Q; residual = B + D - Q),
             E(D) + Q ⇌ E(D, Q)) :: OnlyA
        end
    end
    large_pumped = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + A <--> E(A), E(B) + A <--> E(P; residual = A + B - P),
             E(R) + A <--> E(A, R)) :: OnlyA
            (E + B <--> E(B), E(A) + B <--> E(P; residual = A + B - P)) :: OnlyA
            (E + R <--> E(R), E(A) + R <--> E(A, R),
             E(Q, R) + R <--> E(D, R; residual = A + B - P)) :: OnlyA
            (E(; residual = A + B - P) + P ⇌ E(P; residual = A + B - P),
             E(Q, R) + P ⇌ E(D, P; residual = A + B - P)) :: EqualAI
            (E(; residual = A + B - P) + D <--> E(Q, R),
             E(P; residual = A + B - P) + D <--> E(D, P; residual = A + B - P),
             E(R) + D <--> E(D, R)) :: OnlyA
            E(R) + Q ⇌ E(Q, R) :: OnlyA
        end
    end
    for am in (large_admitted, large_island, large_pumped)
        @test _testhelper_choice_count(am) > D._MAX_CHOICES_PER_CONFORMATION
    end
    @test D._admitting_reading(large_admitted) == D._built_reading(large_admitted)
    @test D._guard(large_island) ==
          D._Verdict(D._TOO_MANY_READINGS, D._LIMIT_UNTESTED, 0, 2, 0)
    @test D._guard_message(large_island, D._guard(large_island)) ==
          "AllostericMechanism: the inactive conformation can switch off the tagged " *
          "steps in more than 20000 ways, too many for the guard to search, and none of " *
          "the three ways it tries first fits the rule for switching them off, so it " *
          "cannot show that the rate law is the limit."
    @test D._guard(large_pumped).fault === D._PUMPED
    @test D._guard(large_pumped) ==
          D._admission(large_pumped, D._built_reading(large_pumped))
    names = Dict(name => o for (o, name) in _testhelper_switch_off_names)
    LF, LR = :loses_forward, :loses_reverse
    for (am, column) in ((large_island, [LF, LF, LF, :runs, LR, LR, LF, LF, LF, LF, LF,
                                         :runs, LF, LF]),
                         (large_pumped, [LF, LF, LF, :runs, LR, LF, LF, LF, LR, :runs, LF,
                                         LF, LF, LF]))
        reading = [fill(:runs, length(column)) column]
        @test _testhelper_reading_is_consistent(am, reading)
        @test D._verdict_under(am, map(name -> names[name], reading)) == D._WELL_POSED
        @test D._admitting_reading(am) === nothing
    end

    # The two other ways a search can be too large, each first in a ter-ter tagging that
    # mixes :OnlyA and :OnlyI, found by sampling the gate's ter-ter init mechanisms. In
    # `too_many_pairs` each conformation stays within `_MAX_CHOICES_PER_CONFORMATION`, but
    # its consistent columns, counted by the tests' own rule and equal to the package's,
    # pair up more than `_MAX_READING_PAIRS` times; in `too_many_active` it is the active
    # conformation that has too many choices. No built reading of either is consistent.
    # Neither active conformation turns the reaction over, which no reading changes, so
    # the guard reports that and not the cap: the cap does not hide it.
    function _testhelper_consistent_column_count(am, c)
        n = length(am.catalytic.steps)
        group = _testhelper_kept_components(am, c)
        between = [k for (k, e) in enumerate(am.catalytic.steps)
                   if _testhelper_excludes(am, k, c) && group[e.from] != group[e.to]]
        r = fill(:runs, n, 2)
        for (k, e) in enumerate(am.catalytic.steps)
            _testhelper_excludes(am, k, c) && group[e.from] == group[e.to] &&
                (r[k, c] = e.rapid_equilibrium ? :cannot_switch_off : :stopped)
        end
        ways(k) = [w for w in (:loses_forward, :loses_reverse, :stopped)
                   if !(w === :loses_forward &&
                        group[am.catalytic.steps[k].to] == group[1]) &&
                      !(w === :loses_reverse &&
                        group[am.catalytic.steps[k].from] == group[1]) &&
                      !(w === :stopped && am.catalytic.steps[k].rapid_equilibrium)]
        count(Iterators.product((ways(k) for k in between)...)) do choice
            r[between, c] .= collect(choice)
            ex = _testhelper_exclusion(am, c; reading = r)
            ex.failure === :none && isempty(ex.lost)
        end
    end
    too_many_pairs = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + A <--> E(A), E(D) + A <--> E(P; residual = A + D - P),
             E(Q) + A <--> E(A, Q)) :: OnlyI
            (E + D <--> E(D), E(A) + D <--> E(P; residual = A + D - P),
             E(R) + D <--> E(D, R)) :: OnlyI
            (E + Q <--> E(Q), E(A) + Q <--> E(A, Q),
             E(Q, R) + Q <--> E(B, Q; residual = A + D - P), E(R) + Q <--> E(Q, R)) :: OnlyI
            (E + R <--> E(R), E(D) + R <--> E(D, R), E(Q) + R <--> E(Q, R)) :: OnlyA
            (E(; residual = A + D - P) + P ⇌ E(P; residual = A + D - P),
             E(Q, R) + P ⇌ E(B, P; residual = A + D - P)) :: OnlyA
            (E(; residual = A + D - P) + B ⇌ E(Q, R),
             E(P; residual = A + D - P) + B ⇌ E(B, P; residual = A + D - P),
             E(Q) + B ⇌ E(B, Q)) :: EqualAI
        end
    end
    too_many_active = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            (E + A ⇌ E(A), E(B, R) + A ⇌ E(A, B, R), E(R) + A ⇌ E(A, R)) :: OnlyI
            (E + R <--> E(R), E(A) + R <--> E(A, R),
             E(Q; residual = A + B - Q) + R <--> E(A, B, R)) :: OnlyI
            E(; residual = A + B - Q) + D <--> E(P, R) :: OnlyI
            E(; residual = A + B - Q) + Q ⇌ E(Q; residual = A + B - Q) :: OnlyI
            (E(A) + B <--> E(Q; residual = A + B - Q), E(A, R) + B <--> E(A, B, R),
             E(R) + B <--> E(B, R)) :: OnlyI
            E(R) + P ⇌ E(P, R) :: OnlyI
        end
    end
    pair_columns = [_testhelper_consistent_column_count(too_many_pairs, c) for c in 1:2]
    @test pair_columns == [length(D._consistent_columns(too_many_pairs, c).columns)
                           for c in 1:2]
    @test all(c -> _testhelper_choice_count(too_many_pairs, c) <=
                   D._MAX_CHOICES_PER_CONFORMATION, 1:2)
    @test prod(pair_columns) > D._MAX_READING_PAIRS
    @test _testhelper_choice_count(too_many_active, 1) > D._MAX_CHOICES_PER_CONFORMATION
    for am in (too_many_pairs, too_many_active)
        @test !any(r -> _testhelper_reading_is_consistent(
                       am, _testhelper_reading(am, r)), keys(_testhelper_guard_readings))
    end
    for am in (too_many_pairs, too_many_active)
        @test !EnzymeRates.Mechanisms.has_turnover_cycle(am.catalytic, am.active_steps)
        @test D._guard(am) == D._Verdict(D._NO_ACTIVE_TURNOVER, D._LIMIT_UNTESTED, 0, 1, 1)
    end
    @test D._guard_message(too_many_pairs, D._guard(too_many_pairs)) ==
          "AllostericMechanism: without its :OnlyI steps E + A <--> E(A); " *
          "E + D <--> E(D); E + Q <--> E(Q); E(A) + Q <--> E(A, Q); " *
          "E(A) + D <--> E(P; residual = A + D - P); " *
          "E(D) + A <--> E(P; residual = A + D - P); E(Q) + A <--> E(A, Q); " *
          "E(Q, R) + Q <--> E(B, Q; residual = A + D - P); E(R) + D <--> E(D, R); " *
          "E(R) + Q <--> E(Q, R) the active conformation turns no substrate into " *
          "product, but it must catalyze the reaction."
    @test D._guard_message(too_many_active, D._guard(too_many_active)) ==
          "AllostericMechanism: without its :OnlyI steps E + A ⇌ E(A); E + R <--> E(R); " *
          "E(; residual = A + B - Q) + D <--> E(P, R); " *
          "E(; residual = A + B - Q) + Q ⇌ E(Q; residual = A + B - Q); " *
          "E(A) + R <--> E(A, R); E(A) + B <--> E(Q; residual = A + B - Q); " *
          "E(A, R) + B <--> E(A, B, R); E(B, R) + A ⇌ E(A, B, R); " *
          "E(Q; residual = A + B - Q) + R <--> E(A, B, R); E(R) + A ⇌ E(A, R); " *
          "E(R) + B <--> E(B, R); E(R) + P ⇌ E(P, R) the active conformation turns no " *
          "substrate into product, but it must catalyze the reaction. Tag their groups " *
          ":EqualAI or :NonequalAI."

    # The same two ways with an active conformation that turns the reaction over, built
    # on random-order ter-ter: one full path of bindings, chemistry and releases kept by
    # both conformations; the other nine bindings, a second chemistry through Estar(A, B, D)
    # and a second release order of Q and R tagged. In `cap_active` all three are :OnlyI,
    # which leaves the active conformation more than `_MAX_CHOICES_PER_CONFORMATION` ways;
    # in `cap_pairs` the bindings are :OnlyI and the rest :OnlyA, and the two
    # conformations' consistent columns pair up more than `_MAX_READING_PAIRS` times. Each
    # built reading leaves a group that no switched-off step reaches (Estar(A, B, D) under
    # the two that lose affinity, E(B) or E(Q) when all stop), so none is consistent, and
    # the guard reports the cap.
    cap_active = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            E + A <--> E(A)                 :: EqualAI
            E(A) + B <--> E(A, B)           :: EqualAI
            E(A, B) + D <--> E(A, B, D)     :: EqualAI
            E + B <--> E(B)                 :: OnlyI
            E + D <--> E(D)                 :: OnlyI
            E(A) + D <--> E(A, D)           :: OnlyI
            E(B) + A <--> E(A, B)           :: OnlyI
            E(B) + D <--> E(B, D)           :: OnlyI
            E(D) + A <--> E(A, D)           :: OnlyI
            E(D) + B <--> E(B, D)           :: OnlyI
            E(A, D) + B <--> E(A, B, D)     :: OnlyI
            E(B, D) + A <--> E(A, B, D)     :: OnlyI
            E(A, B, D) <--> E(P, Q, R)      :: EqualAI
            E(A, B, D) <--> Estar(A, B, D)  :: OnlyI
            Estar(A, B, D) <--> E(P, Q, R)  :: OnlyI
            E(P, Q, R) <--> E(Q, R) + P     :: EqualAI
            E(Q, R) <--> E(R) + Q           :: EqualAI
            E(R) <--> E + R                 :: EqualAI
            E(Q, R) <--> E(Q) + R           :: OnlyI
            E(Q) <--> E + Q                 :: OnlyI
        end
    end
    cap_pairs = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, D
        products: P, Q, R
        steps: begin
            E + A <--> E(A)                 :: EqualAI
            E(A) + B <--> E(A, B)           :: EqualAI
            E(A, B) + D <--> E(A, B, D)     :: EqualAI
            E + B <--> E(B)                 :: OnlyI
            E + D <--> E(D)                 :: OnlyI
            E(A) + D <--> E(A, D)           :: OnlyI
            E(B) + A <--> E(A, B)           :: OnlyI
            E(B) + D <--> E(B, D)           :: OnlyI
            E(D) + A <--> E(A, D)           :: OnlyI
            E(D) + B <--> E(B, D)           :: OnlyI
            E(A, D) + B <--> E(A, B, D)     :: OnlyI
            E(B, D) + A <--> E(A, B, D)     :: OnlyI
            E(A, B, D) <--> E(P, Q, R)      :: EqualAI
            E(A, B, D) <--> Estar(A, B, D)  :: OnlyA
            Estar(A, B, D) <--> E(P, Q, R)  :: OnlyA
            E(P, Q, R) <--> E(Q, R) + P     :: EqualAI
            E(Q, R) <--> E(R) + Q           :: EqualAI
            E(R) <--> E + R                 :: EqualAI
            E(Q, R) <--> E(Q) + R           :: OnlyA
            E(Q) <--> E + Q                 :: OnlyA
        end
    end
    @test _testhelper_choice_count(cap_active, 1) > D._MAX_CHOICES_PER_CONFORMATION
    cap_columns = [_testhelper_consistent_column_count(cap_pairs, c) for c in 1:2]
    @test cap_columns ==
          [length(D._consistent_columns(cap_pairs, c).columns) for c in 1:2]
    @test all(c -> _testhelper_choice_count(cap_pairs, c) <=
                   D._MAX_CHOICES_PER_CONFORMATION, 1:2)
    @test prod(cap_columns) > D._MAX_READING_PAIRS
    for am in (cap_active, cap_pairs)
        @test EnzymeRates.Mechanisms.has_turnover_cycle(am.catalytic, am.active_steps)
        @test !any(r -> _testhelper_reading_is_consistent(
                       am, _testhelper_reading(am, r)), keys(_testhelper_guard_readings))
    end
    @test D._guard(cap_active) ==
          D._Verdict(D._TOO_MANY_READINGS, D._LIMIT_UNTESTED, 0, 1, 0)
    @test D._guard_message(cap_active, D._guard(cap_active)) ==
          "AllostericMechanism: the active conformation can switch off the tagged steps " *
          "in more than 20000 ways, too many for the guard to search, and none of the " *
          "three ways it tries first fits the rule for switching them off, so it cannot " *
          "show that the rate law is the limit."
    @test D._guard(cap_pairs) ==
          D._Verdict(D._TOO_MANY_READINGS, D._LIMIT_UNTESTED, 0, 0, 0)
    @test D._guard_message(cap_pairs, D._guard(cap_pairs)) ==
          "AllostericMechanism: the two conformations together can switch off the tagged " *
          "steps in more than 5000 ways, too many for the guard to search, and none of " *
          "the three ways it tries first fits the rule for switching them off, so it " *
          "cannot show that the rate law is the limit."
end

@testset "the graph test equals the exact LP on the uni-uni grid" begin
    # A limit exists exactly when no excluded rapid-equilibrium step is pinned (the reading
    # cannot switch it off, or stops it, though a step keeps its kind) and the
    # steps that lose a direction, each pointing uphill, close no cycle and never point
    # into free E's component; the test oracle's own reading of the rule and its exact
    # simplex over the cycle rows decide the same question, and its exclusion structure
    # names the same failure. So under the guard's first built reading, and under every
    # reading of the inactive conformation's excluded steps, the search's among them.
    grid = _testhelper_uniuni_grid()
    @test length(grid) == 108
    for am in grid
        @test _testhelper_existence_agrees(am)
        n = length(am.catalytic.steps)
        excluded = [k for k in 1:n if _testhelper_excludes(am, k, 2)]
        ways = (:loses_forward, :loses_reverse, :stopped, :cannot_switch_off)
        readings = map(Iterators.product(fill(ways, length(excluded))...)) do o
            r = fill(:runs, n, 2)
            r[excluded, 2] .= collect(o)
            r
        end
        @test all(r -> _testhelper_existence_agrees(am; reading = r), readings)
    end
end

@testset "a rapid-equilibrium step whose sides stay joined cannot be switched off" begin
    rng = Random.MersenneTwister(31)
    FP, U = EnzymeRates.Derivation._FAST_STEP_PINNED, EnzymeRates.Derivation._LIMIT_UNTESTED
    # The uni-uni S binding :OnlyA: the chemistry and the P binding join E to E(S) in the
    # inactive conformation, so the S binding cannot lose its affinity there, and at rapid
    # equilibrium it cannot be stopped either.
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: OnlyA
            E(S) <--> E(P)      :: EqualAI
            E + P ⇌ E(P)        :: EqualAI
        end
    end
    m = am.catalytic
    bind, release, chem = (_testhelper_edge(am, t) for t in
                           ("E + S ⇌ ES", "E + P ⇌ EP", "ES <--> EP"))
    # The chemistry is stored substrate side first, so it runs toward products as stored.
    @test m.toward_products[chem]
    v = EnzymeRates.Derivation._guard(am)
    @test v == EnzymeRates.Derivation._Verdict(FP, U, 0, 2, bind)
    @test EnzymeRates.Derivation._admission_interior(
        am, EnzymeRates.Derivation._built_reading(am)) == v
    @test _testhelper_epath_verdict(am; rng) === :no_path
    line(k) = EnzymeRates.Mechanisms.step_text(m, k)
    @test EnzymeRates.Derivation._guard_message(am, v) ==
          _testhelper_pinned_cause(am, bind, 2) * " Make its group steady-state, or tag " *
          ":OnlyA the steps that join its sides: $(line(release)); $(line(chem))."
    # The fixes pass the guard: the steps that join its sides tagged, and the binding
    # steady-state ("how an excluded step is switched off"); so does the binding untagged,
    # which the message does not offer.
    untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: EqualAI
            E(S) <--> E(P)      :: EqualAI
            E + P ⇌ E(P)        :: EqualAI
        end
    end
    joins_tagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: OnlyA
            E(S) <--> E(P)      :: OnlyA
            E + P ⇌ E(P)        :: OnlyA
        end
    end
    @test EnzymeRates.Derivation._guard(untagged) == EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(joins_tagged) == EnzymeRates.Derivation._WELL_POSED
    # AD1 and AD11: a mixed :OnlyA group that binds B into E(P, Q) and into E(B, Q), forms
    # the untagged steps still make. At steady state the group is stopped, and the law is
    # the limit on every face; at rapid equilibrium (AD1 RE) it cannot be switched off.
    # AD3c: a ::Mod copy's rapid-equilibrium binding whose sides untagged steps join.
    ad1_ss = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                    :: EqualAI
            (E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q))  :: OnlyA
            E(Q) + P ⇌ E(P, Q)                              :: EqualAI
            E + Q ⇌ E(Q)                                    :: EqualAI
            E + B ⇌ E(B)                                    :: EqualAI
            E(B) + Q ⇌ E(B, Q)                              :: EqualAI
        end
    end
    ad1_re = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                :: EqualAI
            (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))    :: OnlyA
            E(Q) + P <--> E(P, Q)                       :: EqualAI
            E + Q ⇌ E(Q)                                :: EqualAI
            E + B ⇌ E(B)                                :: EqualAI
            E(B) + Q ⇌ E(B, Q)                          :: EqualAI
        end
    end
    ad11 = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                            :: EqualAI
            (E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q),
             E(P) + B <--> E(B, P))                                 :: OnlyA
            E(Q) + P ⇌ E(P, Q)                                      :: EqualAI
            E + Q ⇌ E(Q)                                            :: EqualAI
            E + P ⇌ E(P)                                            :: EqualAI
            E + B ⇌ E(B)                                            :: EqualAI
            E(B) + Q ⇌ E(B, Q)                                      :: EqualAI
        end
    end
    ad3c = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: S
        steps: begin
            E + S ⇌ E(S)                    :: EqualAI
            E(S) <--> E(P)                  :: OnlyA
            E + P ⇌ E(P)                    :: EqualAI
            E + S::Mod ⇌ E(S::Mod)          :: OnlyA
            E(S) + S::Mod ⇌ E(S, S::Mod)    :: EqualAI
            E(S::Mod) + S ⇌ E(S, S::Mod)    :: EqualAI
        end
    end
    for x in (ad1_ss, ad11)
        @test EnzymeRates.Derivation._guard(x) == EnzymeRates.Derivation._WELL_POSED
        @test all(_testhelper_admission_agreement(x; rng))
    end
    for (x, text) in ((ad1_re, "EA + B ⇌ EPQ"), (ad3c, "E + Smod ⇌ ESmod"))
        @test EnzymeRates.Derivation._guard(x) ==
              EnzymeRates.Derivation._Verdict(FP, U, 0, 2, _testhelper_edge(x, text))
        @test _testhelper_exclusion(x, 2).failure === :fast
        @test _testhelper_epath_verdict(x; rng) === :no_path
    end
    # On AD3c the message's steps that join the sides run from the binding's unbound side
    # to its bound side, and that fix passes there. On AD1 RE both fixes pass.
    message(x) = EnzymeRates.Derivation._guard_message(x, EnzymeRates.Derivation._guard(x))
    joins = join((EnzymeRates.Mechanisms.step_text(ad3c.catalytic,
                                                   _testhelper_edge(ad3c, t))
                  for t in ("E + S ⇌ ES", "ES + Smod ⇌ ESSmod", "ESmod + S ⇌ ESSmod")),
                 "; ")
    @test endswith(message(ad3c), "Tag :OnlyA the steps that join its sides: $joins.")
    ad1_re_joins = join((EnzymeRates.Mechanisms.step_text(ad1_re.catalytic,
                                                         _testhelper_edge(ad1_re, t))
                         for t in ("E + A ⇌ EA", "E + Q ⇌ EQ", "EQ + P <--> EPQ")), "; ")
    @test endswith(message(ad1_re), "Make its group steady-state, or tag :OnlyA the " *
                                    "steps that join its sides: $ad1_re_joins.")
    ad1_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                    :: EqualAI
            (E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q))  :: EqualAI
            E(Q) + P ⇌ E(P, Q)                              :: EqualAI
            E + Q ⇌ E(Q)                                    :: EqualAI
            E + B ⇌ E(B)                                    :: EqualAI
            E(B) + Q ⇌ E(B, Q)                              :: EqualAI
        end
    end
    ad1_joins_tagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                    :: EqualAI
            (E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q))  :: OnlyA
            E(Q) + P ⇌ E(P, Q)                              :: EqualAI
            E + Q ⇌ E(Q)                                    :: OnlyA
            E + B ⇌ E(B)                                    :: OnlyA
            E(B) + Q ⇌ E(B, Q)                              :: OnlyA
        end
    end
    ad1_re_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                :: EqualAI
            (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))    :: EqualAI
            E(Q) + P <--> E(P, Q)                       :: EqualAI
            E + Q ⇌ E(Q)                                :: EqualAI
            E + B ⇌ E(B)                                :: EqualAI
            E(B) + Q ⇌ E(B, Q)                          :: EqualAI
        end
    end
    ad1_re_joins_tagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                :: EqualAI
            (E(A) + B ⇌ E(P, Q), E(Q) + B ⇌ E(B, Q))    :: OnlyA
            E(Q) + P <--> E(P, Q)                       :: EqualAI
            E + Q ⇌ E(Q)                                :: OnlyA
            E + B ⇌ E(B)                                :: OnlyA
            E(B) + Q ⇌ E(B, Q)                          :: OnlyA
        end
    end
    ad3c_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: S
        steps: begin
            E + S ⇌ E(S)                    :: EqualAI
            E(S) <--> E(P)                  :: OnlyA
            E + P ⇌ E(P)                    :: EqualAI
            E + S::Mod ⇌ E(S::Mod)          :: EqualAI
            E(S) + S::Mod ⇌ E(S, S::Mod)    :: EqualAI
            E(S::Mod) + S ⇌ E(S, S::Mod)    :: EqualAI
        end
    end
    ad3c_joins_tagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: S
        steps: begin
            E + S ⇌ E(S)                    :: OnlyA
            E(S) <--> E(P)                  :: OnlyA
            E + P ⇌ E(P)                    :: EqualAI
            E + S::Mod ⇌ E(S::Mod)          :: OnlyA
            E(S) + S::Mod ⇌ E(S, S::Mod)    :: OnlyA
            E(S::Mod) + S ⇌ E(S, S::Mod)    :: OnlyA
        end
    end
    # Each variant of AD1, AD1 RE, AD3c and AD11 with the group untagged, or with the steps
    # that join the sides of its B or S::Mod binding tagged, passes the guard.
    ad11_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                            :: EqualAI
            (E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q),
             E(P) + B <--> E(B, P))                                 :: EqualAI
            E(Q) + P ⇌ E(P, Q)                                      :: EqualAI
            E + Q ⇌ E(Q)                                            :: EqualAI
            E + P ⇌ E(P)                                            :: EqualAI
            E + B ⇌ E(B)                                            :: EqualAI
            E(B) + Q ⇌ E(B, Q)                                      :: EqualAI
        end
    end
    ad11_joins_tagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                            :: EqualAI
            (E(A) + B <--> E(P, Q), E(Q) + B <--> E(B, Q),
             E(P) + B <--> E(B, P))                                 :: OnlyA
            E(Q) + P ⇌ E(P, Q)                                      :: EqualAI
            E + Q ⇌ E(Q)                                            :: OnlyA
            E + P ⇌ E(P)                                            :: EqualAI
            E + B ⇌ E(B)                                            :: OnlyA
            E(B) + Q ⇌ E(B, Q)                                      :: OnlyA
        end
    end
    for y in (ad1_untagged, ad1_joins_tagged, ad1_re_untagged, ad1_re_joins_tagged,
              ad3c_untagged, ad3c_joins_tagged, ad11_untagged, ad11_joins_tagged)
        @test EnzymeRates.Derivation._guard(y) == EnzymeRates.Derivation._WELL_POSED
    end
    # The message prints a fix only when its variant passes. A dead-end I joins the S
    # binding's sides a second way, through E(I) and E(I, S), so tagging the steps of the
    # path the message would name leaves that binding inside one component, and the guard
    # rejects its steady-state variant too: the message offers no fix. (Untagging the
    # binding passes, but is not one of the two fixes the message offers.)
    deadend = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)        :: OnlyA
            E(S) <--> E(P)      :: EqualAI
            E + P ⇌ E(P)        :: EqualAI
            E + I ⇌ E(I)        :: EqualAI
            E(S) + I ⇌ E(I, S)  :: EqualAI
            E(I) + S ⇌ E(I, S)  :: EqualAI
        end
    end
    deadend_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)        :: EqualAI
            E(S) <--> E(P)      :: EqualAI
            E + P ⇌ E(P)        :: EqualAI
            E + I ⇌ E(I)        :: EqualAI
            E(S) + I ⇌ E(I, S)  :: EqualAI
            E(I) + S ⇌ E(I, S)  :: EqualAI
        end
    end
    deadend_joins_tagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)        :: OnlyA
            E(S) <--> E(P)      :: OnlyA
            E + P ⇌ E(P)        :: OnlyA
            E + I ⇌ E(I)        :: EqualAI
            E(S) + I ⇌ E(I, S)  :: EqualAI
            E(I) + S ⇌ E(I, S)  :: EqualAI
        end
    end
    v = EnzymeRates.Derivation._guard(deadend)
    @test v == EnzymeRates.Derivation._Verdict(FP, U, 0, 2,
                                               _testhelper_edge(deadend, "E + S ⇌ ES"))
    @test _testhelper_epath_verdict(deadend; rng) === :no_path
    @test EnzymeRates.Derivation._guard_message(deadend, v) ==
          _testhelper_pinned_cause(deadend, v.edge, 2)
    @test EnzymeRates.Derivation._guard(deadend_untagged) ==
          EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(deadend_joins_tagged).fault === FP
    # I's binding at E :OnlyA too: untagging it leaves the S binding inside one component,
    # while tagging the steps that join its sides separates every binding, so the message
    # offers only the joins fix.
    two_bindings = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)        :: OnlyA
            E(S) <--> E(P)      :: EqualAI
            E + P ⇌ E(P)        :: EqualAI
            E + I ⇌ E(I)        :: OnlyA
            E(S) + I ⇌ E(I, S)  :: EqualAI
            E(I) + S ⇌ E(I, S)  :: EqualAI
        end
    end
    two_bindings_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)        :: OnlyA
            E(S) <--> E(P)      :: EqualAI
            E + P ⇌ E(P)        :: EqualAI
            E + I ⇌ E(I)        :: EqualAI
            E(S) + I ⇌ E(I, S)  :: EqualAI
            E(I) + S ⇌ E(I, S)  :: EqualAI
        end
    end
    two_bindings_joins_tagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)        :: OnlyA
            E(S) <--> E(P)      :: OnlyA
            E + P ⇌ E(P)        :: OnlyA
            E + I ⇌ E(I)        :: OnlyA
            E(S) + I ⇌ E(I, S)  :: OnlyA
            E(I) + S ⇌ E(I, S)  :: OnlyA
        end
    end
    v = EnzymeRates.Derivation._guard(two_bindings)
    @test v == EnzymeRates.Derivation._Verdict(FP, U, 0, 2,
                                               _testhelper_edge(two_bindings, "E + I ⇌ EI"))
    @test _testhelper_epath_verdict(two_bindings; rng) === :no_path
    joins = join((EnzymeRates.Mechanisms.step_text(two_bindings.catalytic,
                                                   _testhelper_edge(two_bindings, t))
                  for t in ("E + P ⇌ EP", "ES <--> EP", "ES + I ⇌ EIS", "EI + S ⇌ EIS")),
                 "; ")
    @test EnzymeRates.Derivation._guard_message(two_bindings, v) ==
          _testhelper_pinned_cause(two_bindings, v.edge, 2) *
          " Tag :OnlyA the steps that join its sides: $joins."
    @test EnzymeRates.Derivation._guard(two_bindings_untagged).fault === FP
    @test EnzymeRates.Derivation._guard(two_bindings_joins_tagged) ==
          EnzymeRates.Derivation._WELL_POSED
end

@testset "a tagged binding pinned inside one group, or closing a cycle of a reading" begin
    rng = Random.MersenneTwister(37)
    FP, U = EnzymeRates.Derivation._FAST_STEP_PINNED, EnzymeRates.Derivation._LIMIT_UNTESTED
    # The ter cube of seven :OnlyA edges: the untagged E(A) + C, E(C) + B and E(B) + A
    # bindings join E(A) to E(A, B), so its rapid-equilibrium :OnlyA binding of B stays
    # inside one component and cannot be switched off.
    cube = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, C
        products: P
        steps: begin
            E + A ⇌ E(A)             :: OnlyA
            E + B ⇌ E(B)             :: OnlyA
            E + C ⇌ E(C)             :: OnlyA
            E(A) + B ⇌ E(A, B)       :: OnlyA
            E(A) + C ⇌ E(A, C)       :: EqualAI
            E(B) + A ⇌ E(A, B)       :: EqualAI
            E(B) + C ⇌ E(B, C)       :: EqualAI
            E(C) + A ⇌ E(A, C)       :: EqualAI
            E(C) + B ⇌ E(B, C)       :: EqualAI
            E(A, B) + C ⇌ E(A, B, C) :: OnlyA
            E(A, C) + B ⇌ E(A, B, C) :: OnlyA
            E(B, C) + A ⇌ E(A, B, C) :: OnlyA
            E(A, B, C) <--> E(P)     :: OnlyA
            E + P ⇌ E(P)             :: EqualAI
        end
    end
    @test _testhelper_exclusion(cube, 2).failure === :fast
    @test EnzymeRates.Derivation._guard(cube) ==
          EnzymeRates.Derivation._Verdict(FP, U, 0, 2,
                                          _testhelper_edge(cube, "EA + B ⇌ EAB"))
    @test _testhelper_epath_verdict(cube; rng) === :no_path
    # The message offers three fixes, each passing the guard: the B binding at E(A)
    # untagged, the five untagged bindings that join its sides (E(A) + C, E(C) + A,
    # E(C) + B, E(B) + C, E(B) + A) tagged, and the B binding steady-state.
    cube_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, C
        products: P
        steps: begin
            E + A ⇌ E(A)             :: OnlyA
            E + B ⇌ E(B)             :: OnlyA
            E + C ⇌ E(C)             :: OnlyA
            E(A) + B ⇌ E(A, B)       :: EqualAI
            E(A) + C ⇌ E(A, C)       :: EqualAI
            E(B) + A ⇌ E(A, B)       :: EqualAI
            E(B) + C ⇌ E(B, C)       :: EqualAI
            E(C) + A ⇌ E(A, C)       :: EqualAI
            E(C) + B ⇌ E(B, C)       :: EqualAI
            E(A, B) + C ⇌ E(A, B, C) :: OnlyA
            E(A, C) + B ⇌ E(A, B, C) :: OnlyA
            E(B, C) + A ⇌ E(A, B, C) :: OnlyA
            E(A, B, C) <--> E(P)     :: OnlyA
            E + P ⇌ E(P)             :: EqualAI
        end
    end
    cube_joins_tagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, C
        products: P
        steps: begin
            E + A ⇌ E(A)             :: OnlyA
            E + B ⇌ E(B)             :: OnlyA
            E + C ⇌ E(C)             :: OnlyA
            E(A) + B ⇌ E(A, B)       :: OnlyA
            E(A) + C ⇌ E(A, C)       :: OnlyA
            E(B) + A ⇌ E(A, B)       :: OnlyA
            E(B) + C ⇌ E(B, C)       :: OnlyA
            E(C) + A ⇌ E(A, C)       :: OnlyA
            E(C) + B ⇌ E(B, C)       :: OnlyA
            E(A, B) + C ⇌ E(A, B, C) :: OnlyA
            E(A, C) + B ⇌ E(A, B, C) :: OnlyA
            E(B, C) + A ⇌ E(A, B, C) :: OnlyA
            E(A, B, C) <--> E(P)     :: OnlyA
            E + P ⇌ E(P)             :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._guard(cube_untagged) == EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(cube_joins_tagged) ==
          EnzymeRates.Derivation._WELL_POSED
    joins = join((EnzymeRates.Mechanisms.step_text(cube.catalytic,
                                                   _testhelper_edge(cube, t))
                  for t in ("EA + C ⇌ EAC", "EC + A ⇌ EAC", "EC + B ⇌ EBC",
                            "EB + C ⇌ EBC", "EB + A ⇌ EAB")), "; ")
    @test endswith(
        EnzymeRates.Derivation._guard_message(cube, EnzymeRates.Derivation._guard(cube)),
        "Make its group steady-state, or tag :OnlyA the steps that join its sides: " *
        "$joins.")
    # Both bindings of the ordered bi-uni :OnlyA: the chemistry and the P binding join
    # E(A, B) back to E, so E(A) alone is cut off. No built reading switches off the B
    # binding, which takes up B into a form free E's group holds; the one consistent
    # reading loses the A binding's uptake and the B binding's release, both into E(A).
    # Their vanishing directions then bridge E and E(A, B), and the limit depends on the
    # direction.
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A ⇌ E(A)            :: OnlyA
            E(A) + B ⇌ E(A, B)      :: OnlyA
            E(A, B) <--> E(P)       :: EqualAI
            E + P ⇌ E(P)            :: EqualAI
        end
    end
    b = _testhelper_edge(am, "EA + B ⇌ EAB")
    @test _testhelper_exclusion(am, 2).failure === :fast
    @test EnzymeRates.Derivation.exclusion_components(am, :inactive) ==
          [1, 2, 1, 1]  # E, E(A), E(A, B), E(P)
    @test EnzymeRates.Derivation.exclusion_components(am, :active) == [1, 1, 1, 1]
    v = EnzymeRates.Derivation._guard(am)
    @test v == EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._BRIDGE, U, 1, 2, b)
    @test _testhelper_epath_verdict(am; rng) === :no_path
    reported = _testhelper_reading(EnzymeRates.Derivation._guard_with_reading(am).reading)
    @test _testhelper_reading_is_consistent(am, reported)
    @test _testhelper_epath_verdict(am; rng, corners = true, reading = reported) ===
          :direction_dependent
    @test EnzymeRates.Derivation._guard_message(am, v) ==
          "AllostericMechanism: the :OnlyA step E(A) + B ⇌ E(A, B) bridges forms of the " *
          "inactive conformation that its other steps keep apart, so admission cannot " *
          "show that the limit does not depend on how fast that step vanishes. Make the " *
          "group of E(A) + B ⇌ E(A, B) steady-state."
    # Both bindings untagged, or the chemistry and the P binding tagged too, the mechanism
    # passes the guard; untagging the B binding alone leaves the A binding's sides joined.
    untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A ⇌ E(A)            :: EqualAI
            E(A) + B ⇌ E(A, B)      :: EqualAI
            E(A, B) <--> E(P)       :: EqualAI
            E + P ⇌ E(P)            :: EqualAI
        end
    end
    joins_tagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A ⇌ E(A)            :: OnlyA
            E(A) + B ⇌ E(A, B)      :: OnlyA
            E(A, B) <--> E(P)       :: OnlyA
            E + P ⇌ E(P)            :: OnlyA
        end
    end
    one_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P
        steps: begin
            E + A ⇌ E(A)            :: OnlyA
            E(A) + B ⇌ E(A, B)      :: EqualAI
            E(A, B) <--> E(P)       :: EqualAI
            E + P ⇌ E(P)            :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._guard(untagged) == EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(joins_tagged) == EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(one_untagged).fault === FP
    # With I forms, E(A) is still cut off alone, and the one consistent reading loses both
    # bindings into it, as above. Admission passes the faces with I present, and class 3
    # rejects zero I, where the law splits a rapid-equilibrium segment, as the RE-speed
    # oracle finds. Tagging the chemistry and the P binding too leaves the law a second
    # closed class on a face.
    with_i = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P
        modifiers: I
        steps: begin
            E + A ⇌ E(A)                :: OnlyA
            E(A) + B ⇌ E(A, B)          :: OnlyA
            E(A, B) <--> E(P)           :: EqualAI
            E + P ⇌ E(P)                :: EqualAI
            E + I ⇌ E(I)                :: EqualAI
            E(A, B) + I ⇌ E(A, B, I)    :: EqualAI
            E(I) + A ⇌ E(A, I)          :: EqualAI
            E(A, I) + B ⇌ E(A, B, I)    :: EqualAI
        end
    end
    with_i_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P
        modifiers: I
        steps: begin
            E + A ⇌ E(A)                :: EqualAI
            E(A) + B ⇌ E(A, B)          :: EqualAI
            E(A, B) <--> E(P)           :: EqualAI
            E + P ⇌ E(P)                :: EqualAI
            E + I ⇌ E(I)                :: EqualAI
            E(A, B) + I ⇌ E(A, B, I)    :: EqualAI
            E(I) + A ⇌ E(A, I)          :: EqualAI
            E(A, I) + B ⇌ E(A, B, I)    :: EqualAI
        end
    end
    with_i_joins_tagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P
        modifiers: I
        steps: begin
            E + A ⇌ E(A)                :: OnlyA
            E(A) + B ⇌ E(A, B)          :: OnlyA
            E(A, B) <--> E(P)           :: OnlyA
            E + P ⇌ E(P)                :: OnlyA
            E + I ⇌ E(I)                :: EqualAI
            E(A, B) + I ⇌ E(A, B, I)    :: EqualAI
            E(I) + A ⇌ E(A, I)          :: EqualAI
            E(A, I) + B ⇌ E(A, B, I)    :: EqualAI
        end
    end
    v = EnzymeRates.Derivation._guard(with_i)
    @test (v.fault, _testhelper_faces(with_i)[v.face]) ==
          (EnzymeRates.Derivation._SPLIT, Set([:I]))
    @test _testhelper_epath_verdict(with_i; rng) === :no_path
    @test _testhelper_agrees(v.fault, _testhelper_respeed_verdicts(
        with_i; rng = Random.MersenneTwister(38), law = :limit)[v.face])
    @test EnzymeRates.Derivation._guard_message(with_i, v) ==
          "AllostericMechanism: at zero I, the forms {E (active), E(A) (active), " *
          "E(A, B) (active), E(P) (active), E (inactive), E(P) (inactive)} and " *
          "{E(A, B) (inactive)} of one rapid-equilibrium segment are joined only by " *
          "steps that take up a zeroed metabolite (E + I ⇌ E(I) and " *
          "E(A, B) + I ⇌ E(A, B, I)), so the rate law cannot weight one against the other."
    @test EnzymeRates.Derivation._guard(with_i_untagged) ==
          EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(with_i_joins_tagged).fault ===
          EnzymeRates.Derivation._NONUNIQUE
    # Under the guard's built readings the arcs close a cycle between two cut-off parts:
    # E(A) with E(A, Q), and E(Q) with E(A, B) and E(P, Q). The B binding at E(A) enters
    # the second, the A binding at E(Q) the first, so the two affinities cannot vanish
    # together. A reading in which one of them loses its release instead has no cycle,
    # and the search admits the mechanism under it, the oracles agreeing.
    cycle = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)              :: OnlyA
            E(A) + Q ⇌ E(A, Q)        :: EqualAI
            E(A) + B ⇌ E(A, B)        :: OnlyA
            E(A, B) <--> E(P, Q)      :: EqualAI
            E(Q) + P ⇌ E(P, Q)        :: EqualAI
            E(Q) + A ⇌ E(A, Q)        :: OnlyA
            E + Q ⇌ E(Q)              :: OnlyA
        end
    end
    line_c(k) = EnzymeRates.Mechanisms.step_text(cycle.catalytic, k)
    b_c, a_c = (_testhelper_edge(cycle, t) for t in ("EA + B ⇌ EAB", "EQ + A ⇌ EAQ"))
    @test _testhelper_exclusion(cycle, 2).failure === :cycle
    @test EnzymeRates.Derivation._admission_interior(
        cycle, EnzymeRates.Derivation._built_reading(cycle)) ==
          EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._ARC_CYCLE, U, 0, 2, a_c)
    @test _testhelper_epath_verdict(cycle; rng) === :no_path
    @test EnzymeRates.Derivation._guard(cycle) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(cycle; rng, corners = true))
    _testhelper_check_law_against_limit(EnzymeRates.compile_mechanism(cycle))
end

@testset "forms the built readings reach only through tagged chemistry" begin
    rng = Random.MersenneTwister(41)
    D = EnzymeRates.Derivation
    IS, U = EnzymeRates.Derivation._ISLAND, EnzymeRates.Derivation._LIMIT_UNTESTED
    FP, PU = EnzymeRates.Derivation._FAST_STEP_PINNED, EnzymeRates.Derivation._PUMPED
    # Each mechanism here the guard's built readings reject, with an island, a pumped
    # basin or a step they cannot switch off; most the search admits, the oracles agreeing
    # on every face under the reading it admits with.
    function _testhelper_admitted_by_search(x)
        @test D._guard(x) == D._WELL_POSED
        @test !(D._admitting_reading(x) in
                [D._built_reading(x, r) for r in instances(D._SlowUptakeReading)])
        @test all(_testhelper_admission_agreement(x; rng, corners = true))
        _testhelper_check_law_against_limit(EnzymeRates.compile_mechanism(x))
    end
    _testhelper_edge_between(x, a, b) = only(k for (k, e) in enumerate(x.catalytic.steps)
        if Set(EnzymeRates.Mechanisms._form_name.(Ref(x.catalytic), (e.from, e.to))) ==
           Set((a, b)))
    _testhelper_island_verdict(x) =
        _testhelper_epath_verdict(x; rng) in (:direction_dependent, :wrong_limit)
    # Both half-reactions of the ping-pong :OnlyA. Under the first built reading half 2
    # takes up Q into E(B; residual = A - P), which the inactive conformation cuts off, so
    # it loses its affinity and free E reaches the covalent forms along it: no island. At
    # zero A and B, B's binding feeds E(; residual = A - P), whose only way out, half 1,
    # is stopped: a pumped basin, as the ε-path oracle finds under that reading. The
    # search admits the ping-pong under a reading in which half 1 loses its rate into the
    # covalent forms instead.
    both = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)                                       :: EqualAI
            E(A) <--> E(; residual = A - P) + P                   :: OnlyA
            E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
            E(B; residual = A - P) <--> E + Q                     :: OnlyA
        end
    end
    half1 = _testhelper_edge(both, "E_res_plusA_minusP + P <--> EA")
    feed = _testhelper_edge(both, "E_res_plusA_minusP + B <--> EB_res_plusA_minusP")
    @test !both.catalytic.toward_products[half1]
    built = D._built_reading(both)
    v = D._admission(both, built)
    @test v == EnzymeRates.Derivation._Verdict(PU, U, 3, 2, feed)
    @test isempty(_testhelper_exclusion(both, 2).lost)
    @test _testhelper_admission_agrees(
        PU, _testhelper_epath_verdict(both; rng, face = _testhelper_faces(both)[3]))
    @test D._admission_message(both, v, built) ==
          "AllostericMechanism: at zero A and B, " *
          "$(EnzymeRates.Mechanisms.step_text(both.catalytic, feed)) feeds the forms " *
          "E(; residual = A - P) of the inactive conformation, which leave only through " *
          "the vanishing step $(EnzymeRates.Mechanisms.step_text(both.catalytic, half1)) " *
          "while nothing runs back the way it came, so admission cannot show that their " *
          "weight does not depend on how fast that step vanishes."
    _testhelper_admitted_by_search(both)
    @test EnzymeRates.Derivation._guard(_testhelper_tag_sets["pp_V_half2"]) ==
          EnzymeRates.Derivation._WELL_POSED
    # R3 and R4: a second free form behind a tagged isomerization; untagging it passes.
    # Under the built readings R3's isomerization is steady-state, so it is stopped and
    # F, F(I) are an island; R4's is at rapid equilibrium, which they cannot switch off.
    # Each isomerization joins F's group to free E's, and the search lets it lose its rate
    # into F.
    r3 = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)        :: EqualAI
            E(S) <--> E(P)      :: EqualAI
            E + P ⇌ E(P)        :: EqualAI
            E <--> F            :: OnlyA
            F + I ⇌ F(I)        :: EqualAI
        end
    end
    r3_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)        :: EqualAI
            E(S) <--> E(P)      :: EqualAI
            E + P ⇌ E(P)        :: EqualAI
            E <--> F            :: EqualAI
            F + I ⇌ F(I)        :: EqualAI
        end
    end
    r4 = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: EqualAI
            E(S) <--> E(P)      :: NonequalAI
            E + P ⇌ E(P)        :: EqualAI
            F ⇌ E               :: OnlyA
        end
    end
    r4_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: EqualAI
            E(S) <--> E(P)      :: NonequalAI
            E + P ⇌ E(P)        :: EqualAI
            F ⇌ E               :: EqualAI
        end
    end
    iso = _testhelper_edge_between(r3, "E", "F")
    # An isomerization of free forms ties on both sides, so it runs as stored.
    @test r3.catalytic.toward_products[iso]
    @test D._admission_interior(r3, D._built_reading(r3)) ==
          EnzymeRates.Derivation._Verdict(IS, U, 0, 2, iso)
    @test _testhelper_island_verdict(r3)
    _testhelper_admitted_by_search(r3)
    iso = _testhelper_edge_between(r4, "E", "F")
    @test r4.catalytic.toward_products[iso]
    @test D._admission_interior(r4, D._built_reading(r4)) ==
          EnzymeRates.Derivation._Verdict(FP, U, 0, 2, iso)
    @test _testhelper_epath_verdict(r4; rng) === :no_path
    _testhelper_admitted_by_search(r4)
    for fixed in (r3_untagged, r4_untagged)
        @test EnzymeRates.Derivation._guard(fixed) == EnzymeRates.Derivation._WELL_POSED
    end
    # AD2b's P binding at F is at rapid equilibrium and joins F to E(P), in free E's group,
    # so no built reading switches it off. Its consistent readings let it lose its release
    # into F; with P zero its uptake then dies, and it runs one way into F, which empties
    # only through the slow F <--> E: F fills faster than it empties, as the ε-path oracle
    # finds, and the message offers the steady-state fix, which the search admits. In the
    # ping-pong with B's binding and half 1 :OnlyA, the covalent form is an island behind
    # half 1 under the built readings; untagged, half 1 joins it to E, and B's steady-state
    # binding, its sides joined, is stopped, a fix that passes every face, the oracles
    # agreeing.
    ad2b = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: EqualAI
            E(S) <--> E(P)      :: EqualAI
            F + P ⇌ E(P)        :: OnlyA
            F <--> E            :: OnlyA
        end
    end
    bbind = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)                                       :: EqualAI
            E(A) <--> E(; residual = A - P) + P                   :: OnlyA
            E(; residual = A - P) + B <--> E(B; residual = A - P) :: OnlyA
            E(B; residual = A - P) <--> E + Q                     :: EqualAI
        end
    end
    bbind_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)                                       :: EqualAI
            E(A) <--> E(; residual = A - P) + P                   :: EqualAI
            E(; residual = A - P) + B <--> E(B; residual = A - P) :: OnlyA
            E(B; residual = A - P) <--> E + Q                     :: EqualAI
        end
    end
    ad2b_ss = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: EqualAI
            E(S) <--> E(P)      :: EqualAI
            F + P <--> E(P)     :: OnlyA
            F <--> E            :: OnlyA
        end
    end
    pbind = _testhelper_edge(ad2b, "F + P ⇌ EP")
    @test D._admission_interior(ad2b, D._built_reading(ad2b)) ==
          EnzymeRates.Derivation._Verdict(FP, U, 0, 2, pbind)
    @test _testhelper_epath_verdict(ad2b; rng) === :no_path
    v = EnzymeRates.Derivation._guard(ad2b)
    @test v == EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._FILLED, U, 2, 2,
                                               pbind)
    reported = _testhelper_reading(D._guard_with_reading(ad2b).reading)
    @test _testhelper_reading_is_consistent(ad2b, reported)
    @test _testhelper_admission_agrees(
        v.fault, _testhelper_epath_verdict(ad2b; rng, face = _testhelper_faces(ad2b)[2],
                                           corners = true, reading = reported))
    @test EnzymeRates.Derivation._guard_message(ad2b, v) ==
          "AllostericMechanism: at zero P, in the inactive conformation the :OnlyA " *
          "rapid-equilibrium step F + P ⇌ E(P) runs one way into F, at a rate that " *
          "vanishes yet grows with the rapid-equilibrium speed, and F leaves only " *
          "through the slower step F <--> E, so F fills faster than it empties and its " *
          "weight depends on how fast the tagged step vanishes. Make the group of " *
          "F + P ⇌ E(P) steady-state."
    _testhelper_admitted_by_search(ad2b_ss)
    half1 = _testhelper_edge(bbind, "E_res_plusA_minusP + P <--> EA")
    @test D._admission_interior(bbind, D._built_reading(bbind)) ==
          EnzymeRates.Derivation._Verdict(IS, U, 0, 2, half1)
    @test _testhelper_island_verdict(bbind)
    _testhelper_admitted_by_search(bbind)
    @test EnzymeRates.Derivation._guard(bbind_untagged) ==
          EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(bbind_untagged; rng))
    # Half 1 and its copy at E(A, I) share one :OnlyA group, and I's binding at E is
    # :OnlyA too. As in the ping-pong above, under the first built reading, at zero A and
    # B, B's binding feeds the covalent forms, whose only ways out, the two copies of half
    # 1, vanish; the search admits it as it admits that ping-pong. Untagging that group
    # joins E(I) to E, so I's rapid-equilibrium binding cannot be switched off.
    shared_half1 = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A <--> E(A)                                       :: EqualAI
            (E(A) <--> E(; residual = A - P) + P,
             E(A, I) <--> E(I; residual = A - P) + P)             :: OnlyA
            E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
            E(B; residual = A - P) <--> E + Q                     :: OnlyA
            E + I ⇌ E(I)                                          :: OnlyA
            E(I) + A <--> E(A, I)                                 :: EqualAI
            E(; residual = A - P) + I ⇌ E(I; residual = A - P)    :: EqualAI
        end
    end
    shared_half1_untagged = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A <--> E(A)                                       :: EqualAI
            (E(A) <--> E(; residual = A - P) + P,
             E(A, I) <--> E(I; residual = A - P) + P)             :: EqualAI
            E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
            E(B; residual = A - P) <--> E + Q                     :: OnlyA
            E + I ⇌ E(I)                                          :: OnlyA
            E(I) + A <--> E(A, I)                                 :: EqualAI
            E(; residual = A - P) + I ⇌ E(I; residual = A - P)    :: EqualAI
        end
    end
    shared_line(t) = EnzymeRates.Mechanisms.step_text(shared_half1.catalytic,
                                                      _testhelper_edge(shared_half1, t))
    built = D._built_reading(shared_half1)
    v = D._admission(shared_half1, built)
    @test v == EnzymeRates.Derivation._Verdict(
        PU, U, 3, 2,
        _testhelper_edge(shared_half1, "E_res_plusA_minusP + B <--> EB_res_plusA_minusP"))
    @test _testhelper_admission_agrees(
        PU, _testhelper_epath_verdict(shared_half1; rng,
                                      face = _testhelper_faces(shared_half1)[3]))
    @test D._admission_message(shared_half1, v, built) ==
          "AllostericMechanism: at zero A and B, " *
          "$(shared_line("E_res_plusA_minusP + B <--> EB_res_plusA_minusP")) feeds the " *
          "forms E(; residual = A - P), E(I; residual = A - P) of the inactive " *
          "conformation, which leave only through the vanishing steps " *
          "$(shared_line("E_res_plusA_minusP + P <--> EA")); " *
          "$(shared_line("EI_res_plusA_minusP + P <--> EAI")) while nothing runs back " *
          "the way it came, so admission cannot show that their weight does not depend " *
          "on how fast those steps vanish."
    _testhelper_admitted_by_search(shared_half1)
    @test EnzymeRates.Derivation._guard(shared_half1_untagged).fault === FP
end

"""The tag sets of `_testhelper_tag_sets` whose inactive conformation excludes a
rapid-equilibrium step that the guard's first built reading cannot switch off, so
admission's interior rejects that reading (`_FAST_STEP_PINNED`). In AD7b, AD9 and TC1 the
step lies inside one group, so no reading switches it off; BS2's isomerization joins two
groups, and the search admits BS2."""
const _testhelper_pinned_tag_sets = Set(["AD7b", "AD9", "BS2", "TC1"])

@testset "admission's interior rejects only the tag sets with a pinned fast step" begin
    rng = Random.MersenneTwister(2613)
    @testset "$name" for name in sort!(collect(keys(_testhelper_tag_sets)))
        am = _testhelper_tag_sets[name]
        v = EnzymeRates.Derivation._admission_interior(
            am, EnzymeRates.Derivation._built_reading(am))
        if name in _testhelper_pinned_tag_sets
            @test (v.fault, v.conformation) == (EnzymeRates.Derivation._FAST_STEP_PINNED, 2)
            @test _testhelper_epath_verdict(am; rng) === :no_path
        else
            @test v == EnzymeRates.Derivation._WELL_POSED
        end
    end
end

@testset "a mixed :OnlyA group whose limit exists is admitted" begin
    # B's fused binding into E(P, Q) and its plain binding at E(Q) share one :OnlyA group;
    # in the inactive conformation no untagged step joins E(Q) to E(B, Q), so the plain
    # binding's affinity can vanish, and the limit is deletion plus pruning.
    am = _testhelper_tag_sets["mixed_group_limit"]
    @test _testhelper_exclusion(am, 2).failure === :none
    @test _testhelper_lp_existence(am, 2)
    @test EnzymeRates.Derivation._admission_interior(
        am, EnzymeRates.Derivation._built_reading(am)) ==
          EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(am) == EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_epath_verdict(am; rng = Random.MersenneTwister(43)) === :converges
end

"""The tag sets of `_testhelper_tag_sets` that a face rejects under the guard's first built
reading; every other one passes every face under it. The search admits R2, BS1, BS2, BS3,
BS4, R2funnel, R2funnel_SS and fz101_10232 under other readings."""
const _testhelper_face_rejected = Set(["R2", "BS1", "BS2", "BS3", "BS4", "ADbridge",
    "ADbridge_noturnover", "bridge_drain", "bridge_deadend", "R2funnel", "R2funnel_SS",
    "S1", "AD7b", "fz101_10232"])

@testset "R2's limit at zero product depends on how fast its isomerization vanishes" begin
    rng = Random.MersenneTwister(43)
    NO, PU = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._PUMPED
    # Classes 1-3 and admission's interior accept R2, so under the guard's first built
    # reading the face condition alone rejects it (spec section 7): at P zero, F's exits
    # vanish (F <--> E, stopped) or take up P (F + P ⇌ E(P)), and E(P) feeds it by a
    # release whose reverse is dead. The search admits R2 under a reading in which
    # F <--> E loses its rate into F instead, so F drains into free E, and the oracles
    # agree on every face under it.
    am = _testhelper_tag_sets["R2"]
    m = am.catalytic
    @test _testhelper_face_faults(am) == fill(NO, 3)
    @test EnzymeRates.Derivation._admission_interior(
        am, EnzymeRates.Derivation._built_reading(am)) ==
          EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_admission_faults(am) == [NO, PU, NO]
    @test _testhelper_epath_verdict(am; rng, face = _testhelper_faces(am)[2]) ===
          :direction_dependent
    k, iso = _testhelper_edge(am, "F + P ⇌ EP"), _testhelper_edge(am, "F <--> E")
    built = EnzymeRates.Derivation._built_reading(am)
    v = EnzymeRates.Derivation._admission(am, built)
    @test v == EnzymeRates.Derivation._Verdict(
        PU, EnzymeRates.Derivation._LIMIT_UNTESTED, 2, 2, k)
    line(j) = EnzymeRates.Mechanisms.step_text(m, j)
    @test EnzymeRates.Derivation._admission_message(am, v, built) ==
          "AllostericMechanism: at zero P, $(line(k)) feeds the forms F of the inactive " *
          "conformation, which leave only through the vanishing step $(line(iso)) while " *
          "nothing runs back the way it came, so admission cannot show that their weight " *
          "does not depend on how fast that step vanishes."
    @test EnzymeRates.Derivation._guard(am) == EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._admitting_reading(am)[iso, 2] !==
          EnzymeRates.Derivation._STOPPED
    @test all(_testhelper_admission_agreement(am; rng, corners = true))
    _testhelper_check_law_against_limit(EnzymeRates.compile_mechanism(am))
end

@testset "a stranded basin with two vanishing exits names both steps in its message" begin
    # Under the guard's first built reading, at P zero, E(S) feeds F by a release whose
    # reverse takes up P; F and F(I) leave the inactive conformation's basin only through
    # the :OnlyA steps E(I) ⇌ F(I) and F <--> E, which vanish there, so admission's face
    # condition rejects that face. But that reading cannot switch off E(I) ⇌ F(I), at
    # rapid equilibrium, which consistent readings let lose its rate into F(I); under the
    # first of them class 3 rejects zero I, where the law splits a rapid-equilibrium
    # segment, as the RE-speed oracle finds.
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + I ⇌ E(I)        :: NonequalAI
            E + S ⇌ E(S)        :: OnlyA
            E(I) ⇌ F(I)         :: OnlyA
            F <--> E            :: OnlyA
            F + P <--> E(S)     :: NonequalAI
            F + I ⇌ F(I)        :: OnlyA
        end
    end
    rng = Random.MersenneTwister(43)
    k = _testhelper_edge(am, "EI ⇌ FI")
    @test EnzymeRates.Derivation._admission_interior(
        am, EnzymeRates.Derivation._built_reading(am)) ==
          EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._FAST_STEP_PINNED,
                                          EnzymeRates.Derivation._LIMIT_UNTESTED, 0, 2, k)
    @test _testhelper_epath_verdict(am; rng, face = _testhelper_faces(am)[2]) === :no_path
    v = EnzymeRates.Derivation._guard(am)
    @test (v.fault, _testhelper_faces(am)[v.face], v.conformation) ==
          (EnzymeRates.Derivation._SPLIT, Set([:I]), 0)
    @test _testhelper_agrees(v.fault, _testhelper_respeed_verdicts(
        am; rng = Random.MersenneTwister(44), law = :limit)[v.face])
    @test EnzymeRates.Derivation._guard_message(am, v) ==
          "AllostericMechanism: at zero I, the forms {E (active), E(S) (active), " *
          "E (inactive)} and {F (active)} of one rapid-equilibrium segment are joined " *
          "only by steps that take up a zeroed metabolite (E + I ⇌ E(I) and " *
          "F + I ⇌ F(I)), so the rate law cannot weight one against the other."
    faces = EnzymeRates.Derivation._faces(am.catalytic)
    built = EnzymeRates.Derivation._built_reading(am)
    eps = EnzymeRates.Derivation._eps_network(
        am, EnzymeRates.Derivation._re_flags(am.catalytic), built)
    w = EnzymeRates.Derivation._admission_face(am, eps,
                                               EnzymeRates.Derivation._law_network(am),
                                               faces, 2)
    @test (w.fault, w.face, w.conformation) == (EnzymeRates.Derivation._PUMPED, 2, 2)
    @test EnzymeRates.Derivation._admission_message(am, w, built) ==
          "AllostericMechanism: at zero P, F + P <--> E(S) feeds the forms F, F(I) of " *
          "the inactive conformation, which leave only through the vanishing steps " *
          "E(I) ⇌ F(I); F <--> E while nothing runs back the way it came, so admission " *
          "cannot show that their weight does not depend on how fast those steps vanish."
    # With E(I) <--> F(I) steady-state the first built reading stops it, the interior
    # passes, and the face condition rejects the pumped basin at zero P, the ε-path
    # oracle agreeing. The search admits it under a reading in which the :OnlyA steps
    # lose their rates into F's group instead, both oracles agreeing on every face.
    stranded_ss = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + I ⇌ E(I)        :: NonequalAI
            E + S ⇌ E(S)        :: OnlyA
            E(I) <--> F(I)      :: OnlyA
            F <--> E            :: OnlyA
            F + P <--> E(S)     :: NonequalAI
            F + I ⇌ F(I)        :: OnlyA
        end
    end
    built_ss = EnzymeRates.Derivation._built_reading(stranded_ss)
    @test EnzymeRates.Derivation._admission_interior(stranded_ss, built_ss) ==
          EnzymeRates.Derivation._WELL_POSED
    v = EnzymeRates.Derivation._admission(stranded_ss, built_ss)
    @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._PUMPED, 2, 2)
    @test _testhelper_admission_agrees(
        v.fault, _testhelper_epath_verdict(stranded_ss; rng,
                                           face = _testhelper_faces(stranded_ss)[2]))
    @test EnzymeRates.Derivation._admission_message(stranded_ss, v, built_ss) ==
          "AllostericMechanism: at zero P, F + P <--> E(S) feeds the forms F, F(I) of " *
          "the inactive conformation, which leave only through the vanishing steps " *
          "E(I) <--> F(I); F <--> E while nothing runs back the way it came, so " *
          "admission cannot show that their weight does not depend on how fast those " *
          "steps vanish."
    @test EnzymeRates.Derivation._guard(stranded_ss) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(stranded_ss; rng, corners = true))
    _testhelper_check_law_against_limit(EnzymeRates.compile_mechanism(stranded_ss))
end

@testset "admitted tag sets converge on every face" begin
    rng = Random.MersenneTwister(47)
    names = sort!([n for n in keys(_testhelper_tag_sets)
                   if n ∉ _testhelper_face_rejected && n ∉ _testhelper_pinned_tag_sets])
    @test length(names) == 37
    @testset "$name" for name in names
        am = _testhelper_tag_sets[name]
        @test EnzymeRates.Derivation._guard(am) == EnzymeRates.Derivation._WELL_POSED
        @test all(==(EnzymeRates.Derivation._NO_FAULT), _testhelper_admission_faults(am))
        @test all(_testhelper_admission_agreement(am; rng))
    end
    # The search admits these of the tag sets the first built reading rejects.
    @test Set(n for n in union(_testhelper_face_rejected, _testhelper_pinned_tag_sets)
              if EnzymeRates.Derivation._guard(_testhelper_tag_sets[n]) ==
                 EnzymeRates.Derivation._WELL_POSED) ==
          Set(["R2", "BS1", "BS2", "BS3", "BS4", "R2funnel", "R2funnel_SS", "fz101_10232"])
    # AD9's chemistry and TC1's exchange of B for P, both at rapid equilibrium and :OnlyA,
    # join forms the inactive conformation keeps joined, so it cannot switch them off.
    @testset "$name" for name in ("AD9", "TC1")
        am = _testhelper_tag_sets[name]
        @test EnzymeRates.Derivation._guard(am).fault ===
              EnzymeRates.Derivation._FAST_STEP_PINNED
        @test _testhelper_epath_verdict(am; rng) === :no_path
    end
end

@testset "a sealed basin drains, and a pumped one does not" begin
    rng = Random.MersenneTwister(53)
    NO, PU = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._PUMPED
    function _testhelper_stranded(am, i)
        face = EnzymeRates.Derivation._faces(am.catalytic)[i]
        s = EnzymeRates.Derivation._face_sets(
            EnzymeRates.Derivation._eps_network(
                am, EnzymeRates.Derivation._re_flags(am.catalytic),
                EnzymeRates.Derivation._built_reading(am)),
            EnzymeRates.Derivation._law_network(am).present, face)
        any(s.stranded)
    end
    # A stranded basin entered only reversibly from the law network is sealed: its weight
    # vanishes with the tagged constants, and the face passes. R2P_sealed at P zero,
    # B_uphill and BS4 with X or I zero.
    for (name, sealed) in (("R2P_sealed", [2]), ("B_uphill", [4, 5, 6]), ("BS4", [4, 5, 6]))
        am = _testhelper_tag_sets[name]
        @test all(i -> _testhelper_stranded(am, i), sealed)
        @test all(==(NO), _testhelper_admission_faults(am)[sealed])
    end
    # Under the guard's first built reading, which stops F <--> E, a basin fed by a step
    # whose reverse is dead is pumped: BS2, BS3 and BS4 at P zero. The search admits BS3
    # and BS4 under a reading in which F <--> E loses its rate into F, so F drains into
    # free E; the oracles agree on every face under it.
    @test _testhelper_admission_faults(_testhelper_tag_sets["BS2"]) == [NO, PU, NO]
    @test _testhelper_admission_faults(_testhelper_tag_sets["BS3"]) == [NO, PU, NO]
    @test _testhelper_admission_faults(_testhelper_tag_sets["BS4"]) ==
          [NO, PU, NO, NO, NO, NO]
    for name in ("R2P_sealed", "B_uphill")
        @test all(_testhelper_admission_agreement(_testhelper_tag_sets[name]; rng))
    end
    for name in ("BS3", "BS4")
        am = _testhelper_tag_sets[name]
        @test (name, EnzymeRates.Derivation._guard(am)) ==
              (name, EnzymeRates.Derivation._WELL_POSED)
        @test all(_testhelper_admission_agreement(am; rng, corners = true))
        _testhelper_check_law_against_limit(EnzymeRates.compile_mechanism(am))
    end
    # BS2's isomerization F ⇌ E is at rapid equilibrium, and the built readings cannot
    # switch it off. It joins F's group to free E's, and the search lets it lose its rate
    # into F.
    bs2 = _testhelper_tag_sets["BS2"]
    @test EnzymeRates.Derivation._admission_interior(
        bs2, EnzymeRates.Derivation._built_reading(bs2)).fault ===
          EnzymeRates.Derivation._FAST_STEP_PINNED
    @test _testhelper_epath_verdict(bs2; rng) === :no_path
    @test EnzymeRates.Derivation._guard(bs2) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(bs2; rng, corners = true))
    _testhelper_check_law_against_limit(EnzymeRates.compile_mechanism(bs2))
    # With F <--> E steady-state (R2's graph) the first built reading stops it, the
    # interior passes, and the face condition rejects the pumped basin at zero P; the
    # search admits it as it admits R2, both oracles agreeing on every face.
    bs2_ss = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: OnlyA
            E(S) <--> E(P)      :: EqualAI
            F + P ⇌ E(P)        :: EqualAI
            F <--> E            :: OnlyA
        end
    end
    v = EnzymeRates.Derivation._admission(bs2_ss,
                                          EnzymeRates.Derivation._built_reading(bs2_ss))
    @test (v.fault, v.face, v.conformation) == (PU, 2, 2)
    @test _testhelper_admission_faults(bs2_ss) == [NO, PU, NO]
    @test EnzymeRates.Derivation._guard(bs2_ss) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(bs2_ss; rng, corners = true))
    _testhelper_check_law_against_limit(EnzymeRates.compile_mechanism(bs2_ss))
end

@testset "a tagged binding that bridges two basins under a reading" begin
    rng = Random.MersenneTwister(59)
    NO, BR = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._BRIDGE
    SP, TR = EnzymeRates.Derivation._SPLIT, EnzymeRates.Derivation._TRAP
    # ADbridge: the :OnlyA bindings of I join E and E(S) of the inactive conformation, which
    # turns over through its steady-state S binding. With them the fed forms' cheapest
    # in-tree needs one slow step fewer (R1 fails), and the turnover is no potential (T0+
    # fails). With I zero the bridge is dead, and class 3 finds the inhibitor square split.
    am = _testhelper_tag_sets["ADbridge"]
    m = am.catalytic
    @test _testhelper_admission_faults(am) == [BR, BR, BR, SP, SP, TR]
    @test all(_testhelper_admission_agreement(am; rng))
    k = _testhelper_edge(am, "ES + I ⇌ EIS")
    v = EnzymeRates.Derivation._guard(am)
    @test v == EnzymeRates.Derivation._Verdict(
        BR, EnzymeRates.Derivation._LIMIT_UNTESTED, 1, 2, k)
    line = EnzymeRates.Mechanisms.step_text(m, k)
    @test EnzymeRates.Derivation._guard_message(am, v) ==
          "AllostericMechanism: the :OnlyA step $line bridges forms of the inactive " *
          "conformation that its other steps keep apart, so admission cannot show that " *
          "the limit does not depend on how fast that step vanishes. Make the group of " *
          "$line steady-state."
    fixed = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S <--> E(S)                         :: EqualAI
            E(S) <--> E + P                         :: EqualAI
            (E + I <--> E(I), E(S) + I <--> E(I, S)) :: OnlyA
            E(I) + S ⇌ E(I, S)                      :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._guard(fixed) == EnzymeRates.Derivation._WELL_POSED
    # The same bridge with the chemistry :OnlyA turns nothing over in the inactive
    # conformation: R1 still fails, but the first substrate's uptake is a potential on
    # every strongly connected part, so T0+ admits the first three faces.
    noturn = _testhelper_tag_sets["ADbridge_noturnover"]
    @test _testhelper_admission_faults(noturn) == [NO, NO, NO, SP, SP, TR]
    @test all(_testhelper_admission_agreement(noturn; rng))
    for (x, potential) in ((am, false), (noturn, true))
        eps = EnzymeRates.Derivation._eps_network(
            x, EnzymeRates.Derivation._re_flags(x.catalytic),
            EnzymeRates.Derivation._built_reading(x))
        law = EnzymeRates.Derivation._law_network(x)
        row = EnzymeRates.Derivation._turnover_name(x.catalytic)
        for face in EnzymeRates.Derivation._faces(x.catalytic)[1:3]
            s = EnzymeRates.Derivation._face_sets(eps, law.present, face)
            nodes = s.fed .& .!s.stranded
            (; van, held) = EnzymeRates.Derivation._copy_transitions(eps, s.alive, nodes, 2)
            @test EnzymeRates.Derivation._shortens_in_tree(eps, nodes, s.n0, van)
            @test EnzymeRates.Derivation._potential_on_sccs(eps, held, row) == potential
        end
    end
    # bridge_drain: :OnlyA bindings of I let the inactive F return to E.
    @test _testhelper_admission_faults(_testhelper_tag_sets["bridge_drain"]) ==
          [BR, BR, BR, SP, NO, SP]
    @test all(_testhelper_admission_agreement(_testhelper_tag_sets["bridge_drain"]; rng))
    drain = _testhelper_tag_sets["bridge_drain"]
    k = _testhelper_edge(drain, "F + I ⇌ FI")
    v = EnzymeRates.Derivation._guard(drain)
    @test v == EnzymeRates.Derivation._Verdict(
        BR, EnzymeRates.Derivation._LIMIT_UNTESTED, 1, 2, k)
    @test EnzymeRates.Derivation._guard_message(drain, v) ==
          "AllostericMechanism: the :OnlyA step F + I ⇌ F(I) bridges forms of the " *
          "inactive conformation that its other steps keep apart, so admission cannot " *
          "show that the limit does not depend on how fast that step vanishes. Make the " *
          "group of F + I ⇌ F(I) steady-state."
    drain_fixed = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)        :: EqualAI
            E(S) <--> F + P     :: EqualAI
            F <--> E            :: OnlyA
            E + I ⇌ E(I)        :: OnlyA
            F + I <--> F(I)     :: OnlyA
            F(I) ⇌ E(I)         :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._guard(drain_fixed) == EnzymeRates.Derivation._WELL_POSED
    # R2funnel: under the guard's first built reading, which stops F <--> E, at P zero the
    # trap {F, F(S)} drains only through a slow step out of a fast transient, and the
    # bridge stays with the S binding steady-state (R2funnel_SS); the limit there depends
    # on the direction. The search admits both under a reading in which F <--> E loses
    # its rate into F instead, the oracles agreeing on every face.
    for name in ("R2funnel", "R2funnel_SS")
        x = _testhelper_tag_sets[name]
        @test _testhelper_admission_faults(x) == [NO, BR, NO]
        w = EnzymeRates.Derivation._admission(x, EnzymeRates.Derivation._built_reading(x))
        @test (w.fault, w.face, w.conformation) === (BR, 2, 2)
        @test _testhelper_epath_verdict(x; rng, face = _testhelper_faces(x)[2],
                                        corners = true) === :direction_dependent
        @test EnzymeRates.Derivation._guard(x) == EnzymeRates.Derivation._WELL_POSED
        @test all(_testhelper_admission_agreement(x; rng, corners = true))
        _testhelper_check_law_against_limit(EnzymeRates.compile_mechanism(x))
    end
end

@testset "the census ter-bi mechanism passes by a potential on every component" begin
    rng = Random.MersenneTwister(61)
    am = _testhelper_tag_sets["terbi_d75d2415361c4dae"]
    m = am.catalytic
    @test EnzymeRates.Derivation._guard(am) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(am; rng))
    # With P zero its :OnlyA bindings shorten the inactive conformation's in-tree (R1
    # fails), yet the first substrate's uptake is a potential on each strongly connected
    # part of its kept and fast transitions (T0+): the census's case for T0+.
    eps = EnzymeRates.Derivation._eps_network(am, EnzymeRates.Derivation._re_flags(m),
                                              EnzymeRates.Derivation._built_reading(am))
    law = EnzymeRates.Derivation._law_network(am)
    row = EnzymeRates.Derivation._turnover_name(m)
    faces = EnzymeRates.Derivation._faces(m)
    by_potential = map(faces) do face
        s = EnzymeRates.Derivation._face_sets(eps, law.present, face)
        nodes = s.fed .& .!s.stranded
        (; van, held) = EnzymeRates.Derivation._copy_transitions(eps, s.alive, nodes, 2)
        EnzymeRates.Derivation._shortens_in_tree(eps, nodes, s.n0, van) &&
            EnzymeRates.Derivation._potential_on_sccs(eps, held, row)
    end
    @test [Set(EnzymeRates.Derivation._face_names(m, f)) for f in faces[by_potential]] ==
          [Set([:P, :Q]), Set([:P])]
end

@testset "conservative rejections the oracle admits stay rejected under their reading" begin
    rng = Random.MersenneTwister(67)
    NO, BR, PU = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._BRIDGE,
        EnzymeRates.Derivation._PUMPED
    _testhelper_corner_verdict(am, i) =
        _testhelper_epath_verdict(am; rng, face = _testhelper_faces(am)[i], corners = true)
    # R-B1, an equilibrium bridge: the bridge joins E and the dead end E(J), which deletion
    # already holds at its equilibrium ratio, so the limit is deletion's, but no potential
    # proves it. The steady-state fix passes and the message offers it.
    am = _testhelper_tag_sets["bridge_deadend"]
    @test _testhelper_admission_faults(am) == [BR, BR, BR, NO, NO, NO]
    @test all(i -> _testhelper_corner_verdict(am, i) === :converges, 1:3)
    @test all(_testhelper_admission_agreement(am; rng)[4:6])
    k = _testhelper_edge(am, "EJ + I ⇌ EIJ")
    v = EnzymeRates.Derivation._guard(am)
    @test v == EnzymeRates.Derivation._Verdict(
        BR, EnzymeRates.Derivation._LIMIT_UNTESTED, 1, 2, k)
    @test endswith(EnzymeRates.Derivation._guard_message(am, v),
                   "Make the group of " *
                   "$(EnzymeRates.Mechanisms.step_text(am.catalytic, k)) steady-state.")
    fixed = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S <--> E(S)         :: EqualAI
            E(S) <--> E + P         :: EqualAI
            E + J <--> E(J)         :: EqualAI
            E + I ⇌ E(I)            :: OnlyA
            E(J) + I <--> E(I, J)   :: OnlyA
            E(I) + J ⇌ E(I, J)      :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._guard(fixed) == EnzymeRates.Derivation._WELL_POSED
    # R-B2, under the guard's first built reading, which stops E <--> F: a stranded basin
    # fed only through vanishing sources, BS1 with I zero and P or S zero, fz101_10232 on
    # its two faces with S zero. Every corner direction converges. The search admits both
    # under a reading in which E <--> F loses its rate into F, the oracles agreeing on
    # every face.
    bs1 = _testhelper_tag_sets["BS1"]
    @test _testhelper_admission_faults(bs1) == [NO, NO, NO, NO, PU, PU]
    @test all(i -> _testhelper_corner_verdict(bs1, i) === :converges, 5:6)
    # The step that rejects BS1 under that reading has both its forms in the basin: it
    # runs one way among them rather than feeding them.
    built = EnzymeRates.Derivation._built_reading(bs1)
    v = EnzymeRates.Derivation._admission(bs1, built)
    @test (v.fault, v.face, v.conformation) == (PU, 5, 2)
    @test EnzymeRates.Derivation._admission_message(bs1, v, built) ==
          "AllostericMechanism: at zero I and P, F + P ⇌ F(P) runs one way among the " *
          "forms F, F(P), F(S) of the inactive conformation, which leave only through " *
          "the vanishing step E <--> F while nothing runs back the way it came, so " *
          "admission cannot show that their weight does not depend on how fast that " *
          "step vanishes."
    fz = _testhelper_tag_sets["fz101_10232"]
    @test _testhelper_admission_faults(fz) == [NO, PU, PU, NO, PU, PU]
    @test _testhelper_corner_verdict(fz, 3) === :converges
    # On both faces with P zero, I zero or not (faces 2 and 5), that reading's rejection
    # is right: at face 2 F's only alive exit, the blocked F <--> E, can vanish faster
    # than its inflow, which one corner direction shows.
    @test _testhelper_corner_verdict(fz, 2) !== :converges
    # The faces with I zero as well come last, so the earlier draws of `rng` stay put.
    @test _testhelper_corner_verdict(fz, 5) !== :converges
    @test _testhelper_corner_verdict(fz, 6) === :converges
    for (x, corners) in ((bs1, false), (fz, true))
        @test EnzymeRates.Derivation._guard(x) == EnzymeRates.Derivation._WELL_POSED
        @test all(_testhelper_admission_agreement(x; rng, corners))
        _testhelper_check_law_against_limit(EnzymeRates.compile_mechanism(x))
    end
    # R-B3, the set-aside rule's equilibrium condition. In this ping-pong's inactive
    # conformation the cut-off forms E(A), F(P), F and F(B) (F is E(; residual = A - P))
    # touch free E and E(Q), which only the steady-state E + Q <--> E(Q) joins, so they are
    # not set aside and the interior stays a bridge. Yet the route through them crosses the
    # slow chemistry E(A) <--> F(P), and every corner direction converges there: a
    # conservative rejection.
    D = EnzymeRates.Derivation
    rb3 = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                            :: OnlyA
            E + Q <--> E(Q)                                         :: EqualAI
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)      :: EqualAI
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)      :: OnlyA
            E(A) <--> E(P; residual = A - P)                        :: OnlyA
            E(B; residual = A - P) ⇌ E(Q)                           :: OnlyA
        end
    end
    r = fill(D._RUNS, length(rb3.catalytic.steps), 2)
    r[_testhelper_edge(rb3, "E + A ⇌ EA"), 2] = D._LOSES_FORWARD
    r[_testhelper_edge(rb3, "EA <--> EP_res_plusA_minusP"), 2] = D._LOSES_FORWARD
    r[_testhelper_edge(rb3, "E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP"), 2] =
        D._LOSES_REVERSE
    r[_testhelper_edge(rb3, "EB_res_plusA_minusP ⇌ EQ"), 2] = D._LOSES_REVERSE
    @test D._admission(rb3, r) ==
          D._Verdict(BR, D._LIMIT_UNTESTED, 1, 2,
                     _testhelper_edge(rb3, "EB_res_plusA_minusP ⇌ EQ"))
    @test _testhelper_epath_verdict(rb3; rng, corners = true,
                                    reading = _testhelper_reading(r)) === :converges
    # R-B4, the set-aside rule's safety clause. In this random-order mechanism's inactive
    # conformation, with P and Q zero, the stopped E + P <--> E(P) runs one way from E(P)
    # into free E along its vanishing direction, a pair at a law form that is not safe;
    # the cut-off forms touch free E alone and their fast classes drain, so that pair
    # alone keeps them, and the face stays a bridge. Every corner direction converges
    # there: a conservative rejection.
    rb4 = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ E(A), E(P) + A ⇌ E(A, P))        :: OnlyA
            (E + P <--> E(P), E(A) + P <--> E(A, P))  :: OnlyA
            E(A) + B <--> E(P, Q)                      :: EqualAI
            E(P) + Q ⇌ E(P, Q)                         :: OnlyA
        end
    end
    r = fill(D._RUNS, length(rb4.catalytic.steps), 2)
    for text in ("E + A ⇌ EA", "EA + P <--> EAP")
        r[_testhelper_edge(rb4, text), 2] = D._LOSES_FORWARD
    end
    r[_testhelper_edge(rb4, "E + P <--> EP"), 2] = D._STOPPED
    for text in ("EP + A ⇌ EAP", "EP + Q ⇌ EPQ")
        r[_testhelper_edge(rb4, text), 2] = D._LOSES_REVERSE
    end
    @test _testhelper_faces(rb4)[2] == Set([:P, :Q])
    let eps = D._eps_network(rb4, D._re_flags(rb4.catalytic), r),
        law = D._law_network(rb4), faces = D._faces(rb4.catalytic)
        @test D._without_set_aside_groups(eps, law.present, faces[2]) === eps
    end
    @test D._admission(rb4, r) ==
          D._Verdict(BR, D._LIMIT_UNTESTED, 2, 2, _testhelper_edge(rb4, "EP + Q ⇌ EPQ"))
    @test _testhelper_epath_verdict(rb4; rng, face = _testhelper_faces(rb4)[2],
                                    corners = true, reading = _testhelper_reading(r)) ===
          :converges
end

@testset "S1 passes admission's interior and fails class 3, as AD7b's faces do" begin
    rng = Random.MersenneTwister(71)
    NO, NU = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._NONUNIQUE
    # S1's stray free form F binds B at steady state: with A or B zero the law network has
    # a second closed class. The ε-network's closed-class count sees the same split, so
    # admission alone rejects those faces too; the guard reports class 3, which runs first.
    s1 = _testhelper_tag_sets["S1"]
    @test EnzymeRates.Derivation._admission_interior(
        s1, EnzymeRates.Derivation._built_reading(s1)) ==
          EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_admission_faults(s1) == [NO, NO, NU, NU, NU]
    @test EnzymeRates.Derivation._guard(s1).fault === NU
    @test all(_testhelper_admission_agreement(s1; rng))
    v = EnzymeRates.Derivation._admission(s1, EnzymeRates.Derivation._built_reading(s1))
    @test (v.fault, v.face, v.conformation) ===
          (EnzymeRates.Derivation._EPS_NONUNIQUE, 3, 2)
    @test EnzymeRates.Derivation._guard_message(s1, v) ==
          "AllostericMechanism: at zero A and B, the forms F (inactive) hold the enzyme " *
          "they receive and are reached only through tagged steps, so the steady state " *
          "depends on where the enzyme starts."
    # AD7b: with I zero no spanning in-tree reaches the fed forms, so the face condition
    # passes and class 3 decides. Its isomerization E ⇌ F is at rapid equilibrium, and
    # the inactive conformation cannot switch it off, so admission rejects the interior.
    ad7b = _testhelper_tag_sets["AD7b"]
    @test EnzymeRates.Derivation._admission(
        ad7b, EnzymeRates.Derivation._built_reading(ad7b)).fault ===
          EnzymeRates.Derivation._FAST_STEP_PINNED
    @test _testhelper_admission_faults(ad7b) == [NO, NO, NO, NU, NU, NU]
    @test _testhelper_epath_verdict(ad7b; rng) === :no_path
    # Class 3's faces agree with the RE-speed oracle.
    @test all(_testhelper_admission_agreement(ad7b; rng)[4:6])
    # With E <--> F steady-state it is stopped: admission passes every face, and the guard
    # itself reaches class 3's second closed class with I zero; both oracles agree on
    # every face.
    ad7b_ss = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A
        products: P
        modifiers: I
        steps: begin
            E + A <--> E(A)                 :: EqualAI
            (E + I ⇌ E(I), F + I ⇌ E(I))    :: EqualAI
            E + P <--> E(P)                 :: EqualAI
            E <--> F                        :: OnlyA
            E(A) <--> E(P)                  :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._admission(
        ad7b_ss, EnzymeRates.Derivation._built_reading(ad7b_ss)) ==
          EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_admission_faults(ad7b_ss) == [NO, NO, NO, NU, NU, NU]
    v = EnzymeRates.Derivation._guard(ad7b_ss)
    @test (v.fault, v.face, v.conformation) == (NU, 4, 2)
    @test all(_testhelper_admission_agreement(ad7b_ss; rng))
end

@testset "joint R1 catches a tree that only both conformations' bridges shorten" begin
    # The in-tree costs where each conformation passes R1 and both together fail it,
    # written in the guard's own type: the smallest network that shows the pattern, for
    # the in-tree arithmetic; "joint R1 rejects a mechanism whose bridges cut the tree
    # only together" tests the check on a mechanism.
    # Copy A holds E and X (nodes 1, 2), copy I holds E, Y and Z (nodes 4, 5, 6), and the
    # flip joins the two copies of E. Alone, copy A's vanishing E → X and copy I's
    # vanishing Y → E each leave the cheapest spanning in-tree at one slow transition;
    # together they make it free, rooted at X. Transitions in order: E→X slow, X→E slow,
    # E→X fast and vanishing (copy A); E→Y fast, Y→Z slow, Z→Y fast, Z→E slow, Y→E fast and
    # vanishing (copy I); the flip. Only the transitions these costs need are present, so
    # most lack the reverse partner a mechanism's would have; the face condition reads the
    # partners only on a stranded basin, which this network has none of.
    net = EnzymeRates.Derivation._Network(3, 2, BitVector([1, 1, 0, 1, 1, 1]),
        [1, 2, 1, 4, 5, 6, 6, 5, 1, 4],                  # from
        [2, 1, 2, 5, 6, 5, 4, 4, 4, 1],                  # to
        [1, 1, 2, 3, 4, 5, 6, 7, 0, 0],                  # edge
        [1, 1, 1, 2, 2, 2, 2, 2, 0, 0],                  # conformation
        BitVector([1, 0, 1, 1, 1, 1, 1, 1, 1, 0]),       # forward
        BitVector([0, 0, 1, 1, 0, 1, 0, 1, 1, 1]),       # fast
        BitVector([0, 0, 1, 0, 0, 0, 0, 1, 0, 0]),       # vanishing
        zeros(Int, 1, 10), zeros(UInt64, 10))
    s = EnzymeRates.Derivation._face_sets(net, net.present, UInt64(0))
    nodes = s.fed .& .!s.stranded
    @test nodes == net.present
    @test EnzymeRates.Derivation._fewest_slow(net, nodes, s.n0) == 1
    for c in 1:2
        (; van) = EnzymeRates.Derivation._copy_transitions(net, s.alive, nodes, c)
        @test count(van) == 1
        @test !EnzymeRates.Derivation._shortens_in_tree(net, nodes, s.n0, van)
    end
    @test EnzymeRates.Derivation._fewest_slow(net, nodes, s.alive) == 0
    @test EnzymeRates.Derivation._face_condition(net, net.present, UInt64(0), 1) ==
          EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._BRIDGE,
                                          EnzymeRates.Derivation._LIMIT_UNTESTED, 0, 0, 0)
    # The in-tree test reports -1 when no in-tree spans the nodes: nodes 2 and 6 share no
    # transition.
    @test EnzymeRates.Derivation._fewest_slow(
        net, BitVector([0, 1, 0, 0, 0, 1]), trues(10)) == -1
    # T0+'s potential: the copy-I cycle E → Y → Z → E taking up one first substrate is no
    # potential; without its first transition the cycle is broken and one exists.
    change = zeros(Int, 1, 10)
    change[1, 4] = 1
    turning = EnzymeRates.Derivation._Network(
        net.nforms, net.ncopies, net.present, net.from, net.to, net.edge, net.conformation,
        net.forward, net.fast, net.vanishing, change, net.takes)
    @test EnzymeRates.Derivation._potential_on_sccs(net, trues(10), 1)
    @test !EnzymeRates.Derivation._potential_on_sccs(turning, trues(10), 1)
    broken = BitVector([1, 1, 1, 0, 1, 1, 1, 1, 1, 1])
    @test EnzymeRates.Derivation._potential_on_sccs(turning, broken, 1)
end

@testset "the guard rejects a steady-state dead end whose law is 0/0 at zero products" begin
    # E(P) leaves only by binding Q, so every other form's weight holds Q, and E(P) is fed
    # only from E(P, Q), which E(Q) + P ⇌ E(P, Q) ties to P: each term of the denominator
    # holds P or Q, so the law is 0/0 with both zero (face 2), where the law's segment
    # chain has two closed classes.
    m = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q <--> E(Q)
            E(P) + Q <--> E(P, Q)
        end
    end
    @test !_testhelper_is_valid(m)
    v = EnzymeRates.Derivation._guard(m)
    @test (v.fault, v.face) == (EnzymeRates.Derivation._LAW_MULTI, 2)
    # The negative control of `_testhelper_finite_on_every_face`, which cannot derive
    # `rate_equation` for a mechanism the guard rejects: at the concentrations the check
    # reads, the law read without the guard is NaN exactly where P and Q are zero.
    params = Dict{Symbol, Float64}(k => 1.0
                                   for k in EnzymeRates.Derivation._solve(m).fitted)
    params[:Keq] = 2.0
    faces = _testhelper_face_concentrations(m; rng = Random.MersenneTwister(5))
    nan = [isnan(_testhelper_law_value(m, params, Dict(pairs(c)))) for c in faces]
    @test count(nan) == 1
    @test nan == [c.P == 0 && c.Q == 0 for c in faces]
    # The same dead end with E(P) + Q <--> E(P, Q) in the active conformation only. The
    # active denominator holds P or Q in every term, as above, and the inactive one is
    # weighted by the active free-enzyme weight, which holds Q in every term: the law is
    # 0/0 with P and Q zero, though the inactive denominator alone keeps a constant term.
    # In the inactive conformation E(P) joins free E only through that :OnlyA step, which
    # the guard's built readings stop, leaving an island; the consistent reading lets it
    # lose its rate into E(P), and class 3 rejects zero P and Q, where the law's segment
    # chain has two closed classes, as for the plain mechanism.
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ E(S)          :: EqualAI
            E(S) <--> E(P, Q)     :: EqualAI
            E(Q) + P ⇌ E(P, Q)    :: EqualAI
            E + Q <--> E(Q)       :: EqualAI
            E(P) + Q <--> E(P, Q) :: OnlyA
        end
    end
    @test !_testhelper_is_valid(am)
    @test EnzymeRates.Derivation._admission_interior(
        am, EnzymeRates.Derivation._built_reading(am)).fault ===
          EnzymeRates.Derivation._ISLAND
    w = EnzymeRates.Derivation._guard(am)
    @test (w.fault, w.face) == (EnzymeRates.Derivation._LAW_MULTI, 2)
end

@testset "every spec's rate equation is finite on every face" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test _testhelper_finite_on_every_face(spec.mechanism;
                                               rng = Random.MersenneTwister(5))
    end
end

@testset "every sampled census mechanism is finite on every face" begin
    # An aggregate pin over a seed set: 300 mechanisms sampled from init_mechanisms and one
    # expand_mechanisms round of the bi-bi and ping-pong reactions.
    bibi = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    pingpong = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    mechs = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[]
    for rxn in (bibi, pingpong)
        seeds = EnzymeRates.Enumeration.init_mechanisms(rxn)
        append!(mechs, seeds, EnzymeRates.Enumeration.expand_mechanisms(seeds, rxn))
    end
    sample = mechs[randperm(Random.MersenneTwister(6), length(mechs))[1:300]]
    @test isempty(filter(m -> !_testhelper_finite_on_every_face(
                             m; rng = Random.MersenneTwister(5)), sample))
end

@testset "a redundant steady-state binding traps S1-bound enzyme at zero substrates" begin
    rng = Random.MersenneTwister(73)
    NO, TR, SP = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._TRAP,
        EnzymeRates.Derivation._SPLIT
    # With S1 and S2 zero, E(S1, S2) feeds E(S1) by the fast release of S2, and E(S1)
    # leaves only by its slow release of S1: a trap, whose limit depends on the RE speeds.
    # With S2 zero, the lowest forms are joined only through E(S1, S2): a split.
    m = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S1, S2
        products: P
        steps: begin
            E + S1 <--> E(S1)
            E + S2 ⇌ E(S2)
            E(S1) + S2 ⇌ E(S1, S2)
            E(S2) + S1 ⇌ E(S1, S2)
            E(S1, S2) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    @test _testhelper_face_faults(m) == [NO, NO, TR, NO, SP]
    @test _testhelper_respeed_verdicts(m; rng, law = :limit) ==
          [:agree, :agree, :pattern, :agree, :pattern]
    feed = _testhelper_edge(m, "ES1 + S2 ⇌ ES1S2")
    v = EnzymeRates.Derivation._guard(m)
    @test (v.fault, v.limit, v.face, v.edge) ===
          (TR, EnzymeRates.Derivation._LIMIT_SPEED, 3, feed)
    @test occursin(EnzymeRates.Mechanisms.step_text(m, feed),
                   EnzymeRates.Derivation._guard_message(m, v))
    # The variant moves the steady-state binding of S1 from free E to E(S2), where it
    # stays redundant inside the one RE segment: every face passes, and its law is the
    # exact steady state, before any pin.
    spec = only(s for s in MECHANISM_TEST_SPECS
                if s.name == "Numerator: redundant SS-bind")
    x = _testhelper_lift(spec.mechanism)
    _testhelper_check_law_against_oracle(spec.mechanism)
    @test _testhelper_face_faults(x) == fill(NO, 5)
    @test _testhelper_respeed_verdicts(spec.mechanism; rng, law = :package) ==
          fill(:agree, 5)
    @test EnzymeRates.Derivation._guard(x) == EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_finite_on_every_face(spec.mechanism; rng = Random.MersenneTwister(5))
    @test length(parameters(spec.mechanism, FittedParameters())) == 6
    # The fix the message prints, the feeding step steady-state, passes as well.
    fed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S1, S2
        products: P
        steps: begin
            E + S1 <--> E(S1)
            E + S2 ⇌ E(S2)
            E(S1) + S2 <--> E(S1, S2)
            E(S2) + S1 ⇌ E(S1, S2)
            E(S1, S2) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    @test EnzymeRates.Derivation._guard(fed) === EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_respeed_verdicts(fed; rng) == fill(:agree, 5)
end

@testset "the LDH i-state 6-group mechanism traps active lactate-bound enzyme" begin
    rng = Random.MersenneTwister(79)
    NO, TR, SP = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._TRAP,
        EnzymeRates.Derivation._SPLIT
    # With Lactate and NAD zero, E(Lactate, NAD) feeds E(Lactate) by the fast release of
    # NAD, and E(Lactate) leaves only slowly: a trap. With NAD zero, a split.
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: NADH, Pyruvate
        products: Lactate, NAD
        catalytic_multiplicity: 4
        steps: begin
            (E + Lactate <--> E(Lactate),
             E(NADH) + Lactate <--> E(Lactate, NADH)) :: NonequalAI
            (E + NAD ⇌ E(NAD),
             E(Lactate) + NAD ⇌ E(Lactate, NAD),
             E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate)) :: EqualAI
            (E + NADH <--> E(NADH),
             E(Lactate) + NADH <--> E(Lactate, NADH),
             E(Pyruvate) + NADH <--> E(NADH, Pyruvate)) :: OnlyA
            (E + Pyruvate ⇌ E(Pyruvate),
             E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate),
             E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate)) :: EqualAI
            E(NAD) + Lactate ⇌ E(Lactate, NAD) :: EqualAI
            E(NADH, Pyruvate) <--> E(Lactate, NAD) :: OnlyA
        end
    end
    @test _testhelper_admission_faults(am) == [NO, TR, NO, NO, NO, NO, SP]
    @test _testhelper_respeed_verdicts(am; rng, law = :limit) ==
          [:agree, :pattern, :agree, :agree, :agree, :agree, :pattern]
    feed = _testhelper_edge(am, "ELactate + NAD ⇌ ELactateNAD")
    v = EnzymeRates.Derivation._guard(am)
    @test (v.fault, v.limit, v.face, v.edge) ===
          (TR, EnzymeRates.Derivation._LIMIT_SPEED, 2, feed)
    @test occursin(EnzymeRates.Mechanisms.step_text(am.catalytic, feed),
                   EnzymeRates.Derivation._guard_message(am, v))
    # The variant makes the :EqualAI NAD group steady-state, keeping its six groups: every
    # face passes the guard and agrees with both oracles, and its law is exact.
    spec = only(s for s in MECHANISM_TEST_SPECS
                if s.name == "LDH i-state NonequalAI 6-group, steady-state NAD group")
    x = _testhelper_lift(spec.mechanism)
    _testhelper_check_law_against_oracle(spec.mechanism)
    @test _testhelper_admission_faults(x) == fill(NO, 7)
    @test all(_testhelper_admission_agreement(x; rng))
    @test _testhelper_respeed_verdicts(spec.mechanism; rng, law = :package) ==
          fill(:agree, 7)
    @test EnzymeRates.Derivation._guard(x) == EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_finite_on_every_face(spec.mechanism; rng = Random.MersenneTwister(5))
    @test length(parameters(spec.mechanism, FittedParameters())) == 10
end

# ── The guard wired: macros, exits, flip move and parser ──────────────────────

@testset "the macros reject what the guard rejects when the expanded code runs" begin
    # One mechanism per class: no cycle turns the reaction over (class 1); a cycle of
    # rapid-equilibrium steps does (class 2); a trap at zero products (class 3); a
    # rapid-equilibrium :OnlyA binding whose two sides untagged steps join, which the
    # inactive conformation cannot switch off (class 4). Each `@test_throws`
    # holds the macro call itself, not an `eval`: had a macro thrown at expansion, this
    # file would not load.
    no_cycle = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) <--> F + P
        end
    end
    re_cycle = @_testhelper_unguarded @enzyme_mechanism begin
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
    trap = @_testhelper_unguarded @enzyme_mechanism begin
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
    no_limit = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: EqualAI
            E + P ⇌ E(P)      :: EqualAI
        end
    end
    faults = [EnzymeRates.Derivation._guard(m).fault
              for m in (no_cycle, re_cycle, trap, no_limit)]
    @test faults == [EnzymeRates.Derivation._NO_CYCLE, EnzymeRates.Derivation._RE_TURNOVER,
                     EnzymeRates.Derivation._TRAP, EnzymeRates.Derivation._FAST_STEP_PINNED]
    thrown(m) = ErrorException(
        EnzymeRates.Derivation._guard_message(m, EnzymeRates.Derivation._guard(m)))
    @test_throws thrown(no_cycle) @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) <--> F + P
        end
    end
    @test_throws thrown(re_cycle) @enzyme_mechanism begin
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
    @test_throws thrown(trap) @enzyme_mechanism begin
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
    @test_throws thrown(no_limit) @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)      :: OnlyA
            E(S) <--> E(P)    :: EqualAI
            E + P ⇌ E(P)      :: EqualAI
        end
    end
end

@testset "the parser rejects a mechanism the guard rejects, naming its step" begin
    # An identify CSV row or `show` text can hold an allosteric ping-pong with both
    # half-reactions :OnlyA. The six-step one, which the search admits, reads back as
    # itself. Read back, the text of one the guard rejects raises the guard's message, not
    # a parse error, and the message names the step at fault: in this four-step one, with
    # no limit under any reading ("the guard searches for a reading"), the release of B
    # fills the covalent form faster than it empties at zero A and B.
    both = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)                                          :: EqualAI
            E(A) <--> E(P; residual = A - P)                      :: OnlyA
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)    :: EqualAI
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)    :: EqualAI
            E(B; residual = A - P) ⇌ E(Q)                         :: OnlyA
            E + Q ⇌ E(Q)                                          :: EqualAI
        end
    end
    @test EnzymeRates.Derivation._guard(both) == EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.DSL.parse_mechanism(sprint(show, both)) == both
    no_limit = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(P; residual = A - P)                       :: OnlyA
            E + Q <--> E(Q)                                         :: EqualAI
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)      :: EqualAI
            E(; residual = A - P) + B ⇌ E(Q)                        :: OnlyA
        end
    end
    v = EnzymeRates.Derivation._guard(no_limit)
    half2 = _testhelper_edge(no_limit, "E_res_plusA_minusP + B ⇌ EQ")
    @test (v.fault, v.edge) == (EnzymeRates.Derivation._FILLED, half2)
    msg = EnzymeRates.Derivation._guard_message(no_limit, v)
    @test occursin(EnzymeRates.Mechanisms.step_text(no_limit.catalytic, half2), msg)
    @test_throws ErrorException(msg) EnzymeRates.DSL.parse_mechanism(sprint(show, no_limit))
end

@testset "a second free form keeps its law exact once the conformations flip at E" begin
    # Every admitted tag set of `_testhelper_tag_sets` that holds a free form besides E:
    # R1, AD2a, AD2d, AD5, R2P_sealed and B_uphill join it to E by a steady-state or
    # tagged step, A5, AD7a, AD7c and the AD13 family in rapid equilibrium (its label
    # sorting before or after E). The conformations flip at E alone, so every one
    # derives through its derivation entry, and each law equals the exact limit of the
    # oracle's own network, which prunes from free E by itself (R2P_sealed and B_uphill
    # lose their second free form in the inactive conformation).
    tag_sets = [["R1", "AD2a", "AD2d", "AD5", "R2P_sealed", "B_uphill", "A5", "AD7a",
                 "AD7c"];
                ["AD13_$(chem)_$label" for chem in (:NonequalAI, :OnlyA)
                 for label in (:D, :Ca, :Ea, :F, :Z)]]
    for (i, name) in enumerate(tag_sets)
        @testset "$name" begin
            _testhelper_check_law_against_limit(
                EnzymeRates.compile_mechanism(_testhelper_tag_sets[name]); seed = i)
        end
    end
end

# ── :OnlyI on modifier-binding groups ─────────────────────────────────────────

@testset "an :OnlyI modifier binding mirrors :OnlyA in the active conformation" begin
    fx = _testhelper_onlyi_fixtures()
    # The fixtures but the inhibitor bridge and the essential activator pass the guard,
    # and the oracles agree on every face of four of them: steady-state and
    # rapid-equilibrium steps, a substrate's modifier copy and an activator route.
    admitted = [k for k in keys(fx) if !(k in (:copy_group_bridge, :essential_activator))]
    @test length(admitted) == 11
    for k in admitted
        @test (k, EnzymeRates.Derivation._guard(fx[k]).fault) ==
              (k, EnzymeRates.Derivation._NO_FAULT)
    end
    for (i, k) in enumerate((:competitive_SS, :uncompetitive_RE, :mod_copy_substrate,
                             :nonessential_activator))
        agree = _testhelper_admission_agreement(
            fx[k]; rng = Random.MersenneTwister(500 + i))
        @test (k, all(agree)) == (k, true)
    end
    # The competitive inhibitor tagged :OnlyI is the :OnlyA one with the conformations
    # swapped: each conformation's law keeps what the other's drops.
    onlya = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)            :: EqualAI
            E(S) <--> E + P         :: EqualAI
            E + I ⇌ E(I)            :: OnlyA
        end
    end
    onlyi = fx.competitive_RE
    @test onlyi.active_steps ==
          onlya.inactive_steps
    @test onlyi.inactive_steps ==
          onlya.active_steps
    @test !all(onlyi.active_steps)
end

@testset "an RE binding whose affinity vanishes still closes a fast cycle" begin
    # The B1 set: an :OnlyI binding at rapid equilibrium loses its affinity in the active
    # conformation, yet on the ε-network it stays fast and closes, with the other
    # rapid-equilibrium steps, a cycle that turns the reaction over: class 2 rejects on
    # the interior. The law networks hold no such cycle, so class 3 and the RE-speed
    # oracle agree on every face. The ε-path finds no path (checked on B1a): the inactive
    # conformation cannot switch off B1a's chemistry of the inhibitor complex, an :OnlyA
    # step at rapid equilibrium.
    b1a = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)            :: EqualAI
            E(S) <--> E(P)          :: EqualAI
            E + P ⇌ E(P)            :: EqualAI
            E + I ⇌ E(I)            :: OnlyI
            E(I) + S ⇌ E(I, S)      :: EqualAI
            E(I, S) ⇌ E(I, P)       :: OnlyA
            E(P) + I ⇌ E(I, P)      :: OnlyI
        end
    end
    b1b = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)            :: OnlyA
            E + P ⇌ E(P)            :: OnlyA
            E(S) <--> E(P)          :: EqualAI
            E(S) + I ⇌ E(I, S)      :: OnlyI
            E(I, S) ⇌ E(I, P)       :: EqualAI
            E(P) + I ⇌ E(I, P)      :: OnlyI
        end
    end
    b1d = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)                          :: EqualAI
            E(S) <--> E(P)                        :: EqualAI
            E + P ⇌ E(P)                          :: EqualAI
            (E + I ⇌ E(I), E(P) + I ⇌ E(I, P))    :: OnlyI
            E(I) + S ⇌ E(I, S)                    :: EqualAI
            E(I, S) ⇌ E(I, P)                     :: OnlyA
        end
    end
    b1e = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: P
        steps: begin
            E + S ⇌ E(S)                    :: EqualAI
            E(S) <--> E(P)                  :: EqualAI
            E + P ⇌ E(P)                    :: EqualAI
            E + P::Mod ⇌ E(P::Mod)          :: OnlyI
            E(P::Mod) + S ⇌ E(P::Mod, S)    :: EqualAI
            E(P::Mod, S) ⇌ E(P::Mod, P)     :: OnlyA
            E(P) + P::Mod ⇌ E(P::Mod, P)    :: OnlyI
        end
    end
    b1f = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ E(A)                    :: EqualAI
            E(A) + B ⇌ E(A, B)              :: EqualAI
            E(A, B) <--> E(P, Q)            :: EqualAI
            E(Q) + P ⇌ E(P, Q)              :: EqualAI
            E + Q ⇌ E(Q)                    :: EqualAI
            E(A, B) + I ⇌ E(A, B, I)        :: OnlyI
            E(A, B, I) ⇌ E(I, P, Q)         :: OnlyA
            E(P, Q) + I ⇌ E(I, P, Q)        :: OnlyI
        end
    end
    bridge_REcycle = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)                      :: EqualAI
            E(S) ⇌ F + P                      :: EqualAI
            F <--> E                          :: EqualAI
            (E + I ⇌ E(I), F + I ⇌ F(I))      :: OnlyI
            F(I) ⇌ E(I)                       :: OnlyA
        end
    end
    bridge_REcycle2 = @_testhelper_unguarded @allosteric_mechanism begin
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
    for (i, (label, m)) in enumerate(pairs((; b1a, b1b, b1d, b1e, b1f, bridge_REcycle,
                                               bridge_REcycle2)))
        v = EnzymeRates.Derivation._guard(m)
        @test (label, v.fault, v.face) == (label, EnzymeRates.Derivation._RE_TURNOVER, 0)
        msg = EnzymeRates.Derivation._guard_message(m, v)
        @test (label, occursin("(whose affinity vanishes in the ", msg) &&
                      occursin("turns the reaction over", msg)) == (label, true)
        # The law networks agree with the RE-speed oracle face by face: class 3 accepts
        # every face but bridge_REcycle's zero products, where both see the steady state
        # move with the speeds.
        respeed = _testhelper_respeed_verdicts(m; rng = Random.MersenneTwister(60 + i),
                                               law = :limit)
        @test (label, all(_testhelper_agrees.(_testhelper_face_faults(m), respeed))) ==
              (label, true)
    end
    @test all(==(:agree), _testhelper_respeed_verdicts(
        b1a; rng = Random.MersenneTwister(61), law = :limit))
    @test _testhelper_epath_verdict(b1a; rng = Random.MersenneTwister(62)) === :no_path
    # B1c: the inhibitor complex's chemistry at steady state breaks the fast cycle.
    b1c = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)            :: EqualAI
            E(S) <--> E(P)          :: EqualAI
            E + P ⇌ E(P)            :: EqualAI
            E + I ⇌ E(I)            :: OnlyI
            E(I) + S ⇌ E(I, S)      :: EqualAI
            E(I, S) <--> E(I, P)    :: OnlyA
            E(P) + I ⇌ E(I, P)      :: OnlyI
        end
    end)
    @test all(_testhelper_admission_agreement(b1c; rng = Random.MersenneTwister(63)))
    # With bridge_REcycle's I bindings at steady state no fast cycle is left; class 3
    # then rejects zero products, where the steady state depends on how fast the
    # rapid-equilibrium steps run, and passes the other faces. The guard rejects the
    # interior first: the inactive conformation cannot switch off F(I) ⇌ E(I).
    ssbind = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)                          :: EqualAI
            E(S) ⇌ F + P                          :: EqualAI
            F <--> E                              :: EqualAI
            (E + I <--> E(I), F + I <--> F(I))    :: OnlyI
            F(I) ⇌ E(I)                           :: OnlyA
        end
    end
    v = EnzymeRates.Derivation._guard(ssbind)
    @test (v.fault, v.conformation, v.edge) ==
          (EnzymeRates.Derivation._FAST_STEP_PINNED, 2, _testhelper_edge(ssbind, "EI ⇌ FI"))
    w = _testhelper_face_verdicts(ssbind)[2]
    @test w.face == 2 && w.limit === EnzymeRates.Derivation._LIMIT_SPEED
    @test _testhelper_respeed_verdicts(ssbind; rng = Random.MersenneTwister(64),
                                       law = :limit) ==
          [:agree, :pattern, :agree, :agree, :agree, :agree]
    # With F(I) <--> E(I) steady-state as well it is stopped, the interior passes, and the
    # guard itself rejects zero products, where the steady state depends on the speeds;
    # both oracles agree on every face.
    ssbind_ss = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)                          :: EqualAI
            E(S) ⇌ F + P                          :: EqualAI
            F <--> E                              :: EqualAI
            (E + I <--> E(I), F + I <--> F(I))    :: OnlyI
            F(I) <--> E(I)                        :: OnlyA
        end
    end
    v = EnzymeRates.Derivation._guard(ssbind_ss)
    @test v.face == 2 && v.limit === EnzymeRates.Derivation._LIMIT_SPEED
    @test _testhelper_respeed_verdicts(ssbind_ss; rng = Random.MersenneTwister(65),
                                       law = :limit) ==
          [:agree, :pattern, :agree, :agree, :agree, :agree]
    @test all(_testhelper_admission_agreement(ssbind_ss; rng = Random.MersenneTwister(66)))
end

@testset "a fast cycle that only the inactive conformation keeps is named there" begin
    # The two :OnlyI steps each take up I and J at once, so no built reading switches them
    # off; but they join the active conformation's group of E(I, J, S) and E(I, J, P) to
    # free E's, and consistent readings let them lose their rates into it. With the
    # product binding :OnlyA, the inactive conformation keeps the I, J steps, which still
    # join E(P) to free E: the product binding lies inside one group, the fault no reading
    # changes. No fast cycle is left to class 2.
    m = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ E(S)                    :: EqualAI
            E(S) + I + J ⇌ E(I, J, S)       :: OnlyI
            E(I, J, S) ⇌ E(I, J, P)         :: EqualAI
            E(P) + I + J ⇌ E(I, J, P)       :: OnlyI
            E + P ⇌ E(P)                    :: OnlyA
        end
    end
    k = _testhelper_edge(m, "EIJP ⇌ EP + I + J")
    @test EnzymeRates.Derivation._admission_interior(
        m, EnzymeRates.Derivation._built_reading(m)) ==
          EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._FAST_STEP_PINNED,
                                          EnzymeRates.Derivation._LIMIT_UNTESTED, 0, 1, k)
    p = _testhelper_edge(m, "E + P ⇌ EP")
    v = EnzymeRates.Derivation._guard(m)
    @test (v.fault, v.face, v.conformation, v.edge) ==
          (EnzymeRates.Derivation._FAST_STEP_PINNED, 0, 2, p)
    @test EnzymeRates.Derivation._guard_message(m, v) == _testhelper_pinned_cause(m, p, 2)
    # The same with the chemistry steady-state.
    iso_ss = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ E(S)                    :: EqualAI
            E(S) + I + J ⇌ E(I, J, S)       :: OnlyI
            E(I, J, S) <--> E(I, J, P)      :: EqualAI
            E(P) + I + J ⇌ E(I, J, P)       :: OnlyI
            E + P ⇌ E(P)                    :: OnlyA
        end
    end
    @test EnzymeRates.Derivation._guard(iso_ss).fault ===
          EnzymeRates.Derivation._FAST_STEP_PINNED
    # With the product binding untagged, the inactive conformation keeps every step: the
    # only cycle of rapid-equilibrium steps, which turns the reaction over, runs there, and
    # class 2, which runs before admission, names that conformation.
    kept = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ E(S)                    :: EqualAI
            E(S) + I + J ⇌ E(I, J, S)       :: OnlyI
            E(I, J, S) ⇌ E(I, J, P)         :: EqualAI
            E(P) + I + J ⇌ E(I, J, P)       :: OnlyI
            E + P ⇌ E(P)                    :: EqualAI
        end
    end
    v = EnzymeRates.Derivation._guard(kept)
    @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._RE_TURNOVER, 0, 2)
    @test EnzymeRates.Derivation._guard_message(kept, v) ==
          "AllostericMechanism: the rapid-equilibrium steps E(I, J, P) ⇌ E(P) + I + J, " *
          "E + P ⇌ E(P), E + S ⇌ E(S), E(I, J, S) ⇌ E(S) + I + J and " *
          "E(I, J, S) ⇌ E(I, J, P) close a cycle in the inactive conformation that turns " *
          "the reaction over, so the mechanism has no finite rate. Make one step of the " *
          "cycle steady-state."
end

@testset "the ε-network's second closed class rejects the face" begin
    # fz101_1964: at zero products the :OnlyA steps vanish in the inactive conformation
    # and the :OnlyI binding in the active one, and the alive transitions of the
    # ε-network leave two closed classes, so the steady state depends on where the enzyme
    # starts; the same at zero products and modifiers. Class 3 passes every face: the law
    # networks keep one class each.
    fz = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)            :: OnlyA
            E(S) <--> E + P         :: OnlyA
            E(S) + I ⇌ E(I, S)      :: OnlyI
            E(I, S) ⇌ E(I) + P      :: NonequalAI
        end
    end
    v = EnzymeRates.Derivation._guard(fz)
    @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._EPS_NONUNIQUE, 2, 1)
    # The class the message names, E(I) of the active conformation, which only the
    # vanishing I binding at E(S) feeds, holds no form of the law networks.
    @test EnzymeRates.Derivation._guard_message(fz, v) ==
          "AllostericMechanism: at zero P, the forms E(I) (active) hold the enzyme they " *
          "receive and are reached only through tagged steps, so the steady state " *
          "depends on where the enzyme starts."
    eps = EnzymeRates.Derivation._eps_network(
        fz, EnzymeRates.Derivation._re_flags(fz.catalytic),
        EnzymeRates.Derivation._built_reading(fz))
    law = EnzymeRates.Derivation._law_network(fz)
    s = EnzymeRates.Derivation._face_sets(eps, law.present,
                                          EnzymeRates.Derivation._faces(fz.catalytic)[2])
    named = EnzymeRates.Derivation._class_outside_law(
        EnzymeRates.Derivation._closed_classes(eps, s.alive), law.present)
    @test !any(view(law.present, named))
    NO, NU = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._EPS_NONUNIQUE
    @test _testhelper_admission_faults(fz) == [NO, NU, NO, NO, NU, NO]
    @test all(==(NO), _testhelper_face_faults(fz))
    @test _testhelper_epath_verdict(fz; rng = Random.MersenneTwister(65),
                                    face = Set([:P])) === :nonunique
    @test _testhelper_epath_verdict(fz; rng = Random.MersenneTwister(66),
                                    face = Set([:P, :I])) === :nonunique
end

@testset "an :OnlyI bridge is rejected" begin
    # copy_group_bridge and B2b: the :OnlyI bindings of I at E and E(S) bridge, in the
    # active conformation, forms its other steps keep apart, so the limit depends on how
    # fast the vanishing affinities fall; the ε-path settles on another value. With the
    # bridging group at steady state the guard admits the mechanism, the fix the message
    # prints. At zero modifiers class 3 also rejects copy_group_bridge.
    bridge = _testhelper_onlyi_fixtures().copy_group_bridge
    b2b = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S <--> E(S)                       :: EqualAI
            E(S) <--> E + P                       :: NonequalAI
            (E + I ⇌ E(I), E(S) + I ⇌ E(I, S))    :: OnlyI
            E(I) + S ⇌ E(I, S)                    :: EqualAI
            E + J ⇌ E(J)                          :: OnlyA
        end
    end
    for m in (bridge, b2b)
        v = EnzymeRates.Derivation._guard(m)
        @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._BRIDGE, 1, 1)
        @test _testhelper_admission_faults(m)[1:3] ==
              fill(EnzymeRates.Derivation._BRIDGE, 3)
        @test all(!=(EnzymeRates.Derivation._NO_FAULT), _testhelper_face_faults(m)[4:6])
    end
    msg =
        EnzymeRates.Derivation._guard_message(bridge, EnzymeRates.Derivation._guard(bridge))
    @test startswith(msg, "AllostericMechanism: the :OnlyI step ") &&
          occursin(" bridges forms of the active conformation that its other steps " *
                   "keep apart, so admission cannot show that the limit does not depend " *
                   "on how fast that step vanishes. Make the group of ", msg) &&
          endswith(msg, " steady-state.")
    # Along the oracle's five orders the ε-path settles on another value or depends on
    # the order.
    @test _testhelper_epath_verdict(bridge; rng = Random.MersenneTwister(71)) in
          (:wrong_limit, :direction_dependent)
    fixed = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S <--> E(S)                           :: EqualAI
            E(S) <--> E + P                           :: EqualAI
            (E + I <--> E(I), E(S) + I <--> E(S, I))  :: OnlyI
            E(I) + S ⇌ E(S, I)                        :: EqualAI
        end
    end)
    @test all(_testhelper_admission_agreement(fixed; rng = Random.MersenneTwister(72)))
    # B2a bridges in both conformations, yet admission passes every face; class 3
    # rejects the three faces with zero modifiers, where the limit differs from the law.
    b2a = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S <--> E(S)                       :: EqualAI
            E(S) <--> E + P                       :: EqualAI
            (E + I ⇌ E(I), E(S) + I ⇌ E(I, S))    :: OnlyA
            E(I) + S ⇌ E(I, S)                    :: EqualAI
            (E + J ⇌ E(J), E(S) + J ⇌ E(J, S))    :: OnlyI
            E(J) + S ⇌ E(J, S)                    :: EqualAI
        end
    end
    v = EnzymeRates.Derivation._guard(b2a)
    @test v.face == 4 && v.limit === EnzymeRates.Derivation._LIMIT_DIFFER
    @test _testhelper_admission_faults(b2a)[1:3] ==
          fill(EnzymeRates.Derivation._NO_FAULT, 3)
    @test _testhelper_respeed_verdicts(b2a; rng = Random.MersenneTwister(73),
                                       law = :limit) ==
          [:agree, :agree, :agree, :mismatch, :mismatch, :mismatch]
end

@testset "an essential activator tagged :OnlyI stops catalysis" begin
    # X binds only the inactive conformation, and every turnover passes through E(X):
    # the active conformation, which must catalyze, turns nothing over. Tagging the X
    # binding :EqualAI or :NonequalAI, the message's fix, passes.
    ea = _testhelper_onlyi_fixtures().essential_activator
    v = EnzymeRates.Derivation._guard(ea)
    @test (v.fault, v.face, v.conformation) ==
          (EnzymeRates.Derivation._NO_ACTIVE_TURNOVER, 0, 1)
    @test EnzymeRates.Derivation._guard_message(ea, v) ==
          "AllostericMechanism: without its :OnlyI steps E + X ⇌ E(X) the active " *
          "conformation turns no substrate into product, but it must catalyze the " *
          "reaction. Tag their groups :EqualAI or :NonequalAI."
    for (i, tag) in enumerate((:EqualAI, :NonequalAI))
        fixed = EnzymeRates.AllostericMechanism(eval(:(@allosteric_mechanism begin
            substrates: S
            products: P
            modifiers: X
            steps: begin
                E + X ⇌ E(X)            :: $tag
                E(X) + S ⇌ E(S, X)      :: EqualAI
                E(S, X) <--> E(X) + P   :: EqualAI
            end
        end)))
        @test all(_testhelper_admission_agreement(fixed;
                                                  rng = Random.MersenneTwister(80 + i)))
    end
end

@testset "an :OnlyI activator whose retags fail is offered no fix" begin
    # X binds E and E(S) in the inactive conformation only, and P leaves only from
    # E(S, X), so the active conformation turns nothing over. Tagging the X bindings
    # :EqualAI or :NonequalAI joins E, E(X), E(S, X) and E(S) in one rapid-equilibrium
    # segment that only the X bindings hold together, so each variant splits it at zero
    # X and the message offers no fix.
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: X
        steps: begin
            E + S <--> E(S)         :: EqualAI
            E(S) + X ⇌ E(S, X)      :: OnlyI
            E(S, X) ⇌ E(X) + P      :: EqualAI
            E + X ⇌ E(X)            :: OnlyI
        end
    end
    v = EnzymeRates.Derivation._guard(am)
    @test (v.fault, v.face, v.conformation) ==
          (EnzymeRates.Derivation._NO_ACTIVE_TURNOVER, 0, 1)
    ks = [k for k in eachindex(am.catalytic.steps)
          if am.excluded_active[k]]
    @test !EnzymeRates.Derivation._retag_passes(am, ks, :EqualAI)
    @test !EnzymeRates.Derivation._retag_passes(am, ks, :NonequalAI)
    @test EnzymeRates.Derivation._guard_message(am, v) ==
          "AllostericMechanism: without its :OnlyI steps E + X ⇌ E(X); E(S) + X ⇌ " *
          "E(S, X) the active conformation turns no substrate into product, but it " *
          "must catalyze the reaction."
end

@testset "the active conformation's interior faults name the :OnlyI step" begin
    # The mirror of the :OnlyA interior faults, in conformation 1. Each message offers the
    # fixes whose variants pass. The oracles agree with each verdict.
    rng = Random.MersenneTwister(91)
    FP, U = EnzymeRates.Derivation._FAST_STEP_PINNED, EnzymeRates.Derivation._LIMIT_UNTESTED
    # I binds E and E(S); without E + I ⇌ E(I), the untagged I binding at E(S) and S
    # binding at E(I) still join E to E(I), so the I binding cannot lose its affinity,
    # and at rapid equilibrium the active conformation cannot switch it off. Neither fix
    # the message offers passes; untagging it would.
    no_limit = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)            :: EqualAI
            E(S) <--> E + P         :: EqualAI
            E + I ⇌ E(I)            :: OnlyI
            E(I) + S ⇌ E(I, S)      :: EqualAI
            E(S) + I ⇌ E(I, S)      :: EqualAI
        end
    end
    k = _testhelper_edge(no_limit, "E + I ⇌ EI")
    v = EnzymeRates.Derivation._guard(no_limit)
    @test v == EnzymeRates.Derivation._Verdict(FP, U, 0, 1, k)
    @test EnzymeRates.Derivation._guard_message(no_limit, v) ==
          _testhelper_pinned_cause(no_limit, k, 1)
    @test all(t -> EnzymeRates.Derivation._retag_passes(no_limit, [k], t),
              (:EqualAI, :NonequalAI))
    @test _testhelper_epath_verdict(no_limit; rng) === :no_path
    # I binds E and J binds E(I), both :OnlyI; the untagged J binding at E and I binding
    # at E(J) join E(I, J) back to E, so the J binding at E(I) cannot lose its affinity,
    # and no built reading switches it off. But E(I) is cut off alone, and the search lets
    # both :OnlyI bindings lose their rates into it, the oracles agreeing.
    square = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ E(S)            :: EqualAI
            E(S) <--> E + P         :: EqualAI
            E + I ⇌ E(I)            :: OnlyI
            E(I) + J ⇌ E(I, J)      :: OnlyI
            E + J ⇌ E(J)            :: EqualAI
            E(J) + I ⇌ E(I, J)      :: EqualAI
        end
    end
    k = _testhelper_edge(square, "EI + J ⇌ EIJ")
    @test EnzymeRates.Derivation._admission_interior(
        square, EnzymeRates.Derivation._built_reading(square)) ==
          EnzymeRates.Derivation._Verdict(FP, U, 0, 1, k)
    @test _testhelper_epath_verdict(square; rng) === :no_path
    @test EnzymeRates.Derivation._guard(square) == EnzymeRates.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(square; rng, corners = true))
    _testhelper_check_law_against_limit(EnzymeRates.compile_mechanism(square))
    # E(I) is reached only as E(I, J) gives off J, and the active conformation's other
    # steps make E(I, J): no built reading switches that J binding off, though the search
    # lets it lose its rate into E(I). But the inactive conformation keeps it, and at zero
    # I and J E(I) there holds the enzyme it receives: class 3 finds a second closed
    # class, as the RE-speed oracle does.
    j_release = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ E(S)            :: EqualAI
            E(S) <--> E + P         :: EqualAI
            E + J ⇌ E(J)            :: EqualAI
            E(J) + I ⇌ E(I, J)      :: EqualAI
            E(I) + J ⇌ E(I, J)      :: OnlyI
        end
    end
    k = _testhelper_edge(j_release, "EI + J ⇌ EIJ")
    @test EnzymeRates.Derivation._admission_interior(
        j_release, EnzymeRates.Derivation._built_reading(j_release)) ==
          EnzymeRates.Derivation._Verdict(FP, U, 0, 1, k)
    @test _testhelper_epath_verdict(j_release; rng) === :no_path
    v = EnzymeRates.Derivation._guard(j_release)
    @test (v.fault, _testhelper_faces(j_release)[v.face], v.conformation) ==
          (EnzymeRates.Derivation._NONUNIQUE, Set([:I, :J]), 2)
    @test _testhelper_agrees(v.fault, _testhelper_respeed_verdicts(
        j_release; rng = Random.MersenneTwister(92), law = :limit)[v.face])
    @test EnzymeRates.Derivation._guard_message(j_release, v) ==
          "AllostericMechanism: at zero I and J in the inactive conformation, the forms " *
          "{E(I)} keep enzyme away from free E: every step that leaves them takes up a " *
          "zeroed metabolite (E(I) + J ⇌ E(I, J)), so the steady state is not unique."
    # At steady state the built readings stop that J binding, so E(I) is reached only
    # through a step that vanishes both ways: an island on the :OnlyI binding's unbound
    # side, which collects the active enzyme. The inactive conformation's second closed
    # class at zero I and J rejects it under every reading, as above.
    island = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ E(S)            :: EqualAI
            E(S) <--> E + P         :: EqualAI
            E + J ⇌ E(J)            :: EqualAI
            E(J) + I ⇌ E(I, J)      :: EqualAI
            E(I) + J <--> E(I, J)   :: OnlyI
        end
    end
    k = _testhelper_edge(island, "EI + J <--> EIJ")
    @test EnzymeRates.Derivation._admission_interior(
        island, EnzymeRates.Derivation._built_reading(island)) ==
          EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._ISLAND, U, 0, 1, k)
    @test _testhelper_epath_verdict(island; rng) in (:direction_dependent, :wrong_limit)
    v = EnzymeRates.Derivation._guard(island)
    @test (v.fault, _testhelper_faces(island)[v.face], v.conformation) ==
          (EnzymeRates.Derivation._NONUNIQUE, Set([:I, :J]), 2)
    @test _testhelper_agrees(v.fault, _testhelper_respeed_verdicts(
        island; rng = Random.MersenneTwister(93), law = :limit)[v.face])
    @test EnzymeRates.Derivation._guard_message(island, v) ==
          "AllostericMechanism: at zero I and J in the inactive conformation, the forms " *
          "{E(I)} keep enzyme away from free E: every step that leaves them takes up a " *
          "zeroed metabolite (E(I) + J <--> E(I, J)), so the steady state is not unique."
end

@testset "a pinned step and its mirror get one message, but for conformation and tag" begin
    # Either conformation is offered the same fixes, each with that conformation's tag,
    # when both conformations catalyze in every fix's variant, as here, where the S path
    # stays untagged. I binds E :OnlyI, and the J bindings join E(I) back to E, so the
    # active conformation cannot switch the I binding off. Tagging :OnlyI the steps that
    # join its sides passes (and so would untagging it, which the message does not
    # offer). The mirror image, the I binding :OnlyA, gets the same message in the
    # inactive conformation, with :OnlyA.
    D = EnzymeRates.Derivation
    swapped(text) = replace(text, "the active conformation" => "the inactive conformation",
                            ":OnlyI" => ":OnlyA")
    square = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I, J
        steps: begin
            E + S ⇌ E(S)            :: EqualAI
            E(S) <--> E + P         :: EqualAI
            E + I ⇌ E(I)            :: OnlyI
            E(I) + J ⇌ E(I, J)      :: EqualAI
            E + J ⇌ E(J)            :: EqualAI
            E(J) + I ⇌ E(I, J)      :: EqualAI
        end
    end
    k = _testhelper_edge(square, "E + I ⇌ EI")
    joins = [_testhelper_edge(square, t) for t in ("E + J ⇌ EJ", "EJ + I ⇌ EIJ",
                                                   "EI + J ⇌ EIJ")]
    v = D._guard(square)
    @test (v.fault, v.conformation, v.edge) == (D._FAST_STEP_PINNED, 1, k)
    msg = D._guard_message(square, v)
    @test msg == _testhelper_pinned_cause(square, k, 1) * " Tag :OnlyI the steps that " *
                 "join its sides: E + J ⇌ E(J); E(J) + I ⇌ E(I, J); E(I) + J ⇌ E(I, J)."
    @test all(t -> D._retag_passes(square, [k], t), (:EqualAI, :NonequalAI))
    @test D._retag_passes(square, joins, :OnlyI)
    mirror = _testhelper_mirror(square)
    w = D._guard(mirror)
    @test (w.fault, w.conformation, w.edge) == (D._FAST_STEP_PINNED, 2, k)
    @test D._guard_message(mirror, w) == swapped(msg)
    # The one rule that tells the states apart: the active conformation must catalyze.
    # The steps that join the sides of this I binding include the P release, so tagging
    # them :OnlyI leaves the active conformation no turnover: only the mirror image is
    # offered that fix.
    no_limit = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)            :: EqualAI
            E(S) <--> E + P         :: EqualAI
            E + I ⇌ E(I)            :: OnlyI
            E(I) + S ⇌ E(I, S)      :: EqualAI
            E(S) + I ⇌ E(I, S)      :: EqualAI
        end
    end
    k = _testhelper_edge(no_limit, "E + I ⇌ EI")
    e = no_limit.catalytic.steps[k]
    joins = D._kept_path(no_limit, 1, e.from, e.to)
    retagged = copy(no_limit.mwc_tags)
    retagged[[no_limit.catalytic.steps[j].group for j in joins]] .= :OnlyI
    @test D._guard(EnzymeRates.AllostericMechanism(
              no_limit.catalytic, retagged, 1,
              EnzymeRates.Mechanisms.RegulatorySite[])).fault === D._NO_ACTIVE_TURNOVER
    msg = D._guard_message(no_limit, D._guard(no_limit))
    mirror = _testhelper_mirror(no_limit)
    @test D._guard_message(mirror, D._guard(mirror)) ==
          swapped(msg) * " Tag :OnlyA the steps that join its sides: " *
          "E + P <--> E(S); E(S) + I ⇌ E(I, S); E(I) + S ⇌ E(I, S)."
end

@testset ":OnlyI on a step that does more than bind a modifier" begin
    # Any group may be :OnlyI while the active conformation still catalyzes the reaction.
    # Each mechanism the macro admits has a law equal to the limit on every face, by the
    # ε-path oracle and by the RE-speed oracle against the package's law.
    D = EnzymeRates.Derivation
    rng = Random.MersenneTwister(101)
    # Two isomerization paths from E(S) to E(P); the direct one runs in the inactive
    # conformation alone. Its mirror image, the direct path :OnlyA, is the same model
    # with L inverted, and the guard admits it too.
    two_paths = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)              :: EqualAI
            E(S) <--> E(P)            :: OnlyI
            E(S) ⇌ Estar(S)           :: EqualAI
            Estar(S) <--> Estar(P)    :: EqualAI
            Estar(P) ⇌ E(P)           :: EqualAI
            E + P ⇌ E(P)              :: EqualAI
        end
    end)
    @test all(_testhelper_admission_agreement(two_paths; rng))
    @test all(==(:agree), _testhelper_respeed_verdicts(two_paths; rng))
    mirror = _testhelper_mirror(two_paths)
    @test D._guard(mirror) == D._WELL_POSED
    @test all(_testhelper_admission_agreement(mirror; rng))
    @test all(==(:agree), _testhelper_respeed_verdicts(mirror; rng))
    # An exchange step beside a ternary-complex path, the exchange :OnlyI.
    exchange = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)              :: EqualAI
            E(A) + B <--> E(Q) + P    :: OnlyI
            E(A) + B ⇌ E(A, B)        :: EqualAI
            E(A, B) <--> E(P, Q)      :: EqualAI
            E(Q) + P ⇌ E(P, Q)        :: EqualAI
            E + Q ⇌ E(Q)              :: EqualAI
        end
    end)
    @test all(_testhelper_admission_agreement(exchange; rng))
    @test all(==(:agree), _testhelper_respeed_verdicts(exchange; rng))
    # The only catalytic step :OnlyI: the active conformation cannot catalyze.
    no_turnover = "the active conformation turns no substrate into product"
    @test_throws no_turnover @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)      :: EqualAI
            E(S) <--> E(P)    :: OnlyI
            E + P ⇌ E(P)      :: EqualAI
        end
    end
    # A rapid-equilibrium S binding :OnlyI whose bound form the chemistry and the P
    # binding keep joined to free E: the active conformation cannot switch it off. The
    # message offers no fix: each would leave the active conformation no turnover
    # (untagging the binding passes, but is not one of the fixes it offers).
    pinned_active = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)      :: OnlyI
            E(S) <--> E(P)    :: EqualAI
            E + P ⇌ E(P)      :: EqualAI
        end
    end
    k = _testhelper_edge(pinned_active, "E + S ⇌ ES")
    v = D._guard(pinned_active)
    @test (v.fault, v.conformation, v.edge) == (D._FAST_STEP_PINNED, 1, k)
    msg = D._guard_message(pinned_active, v)
    @test occursin("in the active conformation", msg)
    @test occursin(":OnlyI", msg)
    @test msg == _testhelper_pinned_cause(pinned_active, k, 1)
    @test all(t -> D._retag_passes(pinned_active, [k], t), (:EqualAI, :NonequalAI))
    @test _testhelper_epath_verdict(pinned_active; rng) === :no_path
end

@testset "sampled taggings: the rule agrees with the mass-action oracles" begin
    # An aggregate check over a population: 60 taggings of the plain init mechanisms of
    # the bi-bi and ping-pong reactions, sorted by their text so that the sample does not
    # depend on the enumeration order, each built without the guard. 40 tag one random
    # nonempty set of groups :OnlyA; 20 draw a state for each group, at least one :OnlyI.
    # An admitted law agrees with both oracles on every face. A mechanism rejected for a
    # pinned step has a step that the tests' own statement of the rule cannot switch off,
    # and the ε-path oracle finds no path to the limit.
    D = EnzymeRates.Derivation
    bibi = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    pingpong = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    plain = sort!(unique([m for rxn in (bibi, pingpong)
                          for m in EnzymeRates.Enumeration.init_mechanisms(rxn)
                          if m isa EnzymeRates.Mechanism]); by = m -> sprint(show, m))
    rng = Random.MersenneTwister(20261119)
    function _testhelper_random_tagging(holds_only_i::Bool)
        m = rand(rng, plain)
        G = length(m.groups)
        states = if holds_only_i
            s = Symbol[]
            while !(:OnlyI in s)
                s = rand(rng, [:EqualAI, :OnlyA, :OnlyI], G)
            end
            s
        else
            bits = rand(rng, 1:(2^G - 1))
            [isodd(bits >> (g - 1)) ? :OnlyA : :EqualAI for g in 1:G]
        end
        EnzymeRates.AllostericMechanism(m, states, 1,
                                        EnzymeRates.Mechanisms.RegulatorySite[])
    end
    sample = [[_testhelper_random_tagging(false) for _ in 1:40];
              [_testhelper_random_tagging(true) for _ in 1:20]]
    faults = [D._guard(am).fault for am in sample]
    for (am, fault) in zip(sample, faults)
        if fault === D._NO_FAULT
            @test all(_testhelper_admission_agreement(am; rng))
        elseif fault === D._FAST_STEP_PINNED
            @test :cannot_switch_off in _testhelper_reading(am, :all_lose_affinity)
            @test _testhelper_epath_verdict(am; rng) === :no_path
        end
    end
    @test count(==(D._NO_FAULT), faults) >= 10
    @test count(==(D._FAST_STEP_PINNED), faults) >= 3
end

@testset "sampled taggings: the energy-level rule agrees with the oracles" begin
    # An aggregate check over a population: the nonempty :OnlyA sets of the plain init
    # mechanisms of 8 steps or fewer of the bi-bi and ping-pong reactions, sorted by their
    # text so that the sample does not depend on the enumeration order, each built without
    # the guard. 40 random ones that the guard admits only by a reading it searched, none
    # of its three built readings passing, agree with both oracles on every face, the
    # ε-path oracle along the corner directions under the reading the guard admitted with.
    # At least 5 of them are admitted with a step that loses its reverse direction, which
    # no built reading does.
    D = EnzymeRates.Derivation
    bibi = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    pingpong = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    plain = sort!(unique([m for rxn in (bibi, pingpong)
                          for m in EnzymeRates.Enumeration.init_mechanisms(rxn)
                          if m isa EnzymeRates.Mechanism && length(m.steps) <= 8]);
                  by = m -> sprint(show, m))
    built(am) = [D._built_reading(am, r) for r in instances(D._SlowUptakeReading)]
    rng = Random.MersenneTwister(20261120)
    sample = EnzymeRates.AllostericMechanism[]
    # The 40 take 1,448 draws; the bound only stops a guard that admits fewer from
    # drawing forever.
    for _ in 1:10_000
        length(sample) == 40 && break
        m = rand(rng, plain)
        G = length(m.groups)
        bits = rand(rng, 1:(2^G - 1))
        am = EnzymeRates.AllostericMechanism(
            m, [isodd(bits >> (g - 1)) ? :OnlyA : :EqualAI for g in 1:G], 1,
            EnzymeRates.Mechanisms.RegulatorySite[])
        reading = D._admitting_reading(am)
        searched = reading !== nothing && !(reading in built(am))
        searched && !(am in sample) && push!(sample, am)
    end
    @test length(sample) == 40
    for am in sample
        @test all(_testhelper_admission_agreement(am; rng, corners = true))
    end
    @test count(am -> :loses_reverse in _testhelper_reading(D._admitting_reading(am)),
                sample) >= 5
end

@testset "joint R1 rejects a mechanism whose bridges cut the tree only together" begin
    # The hand-built joint-R1 network's pattern on a mechanism. Half 2 :OnlyA strands the
    # inactive conformation's covalent forms at zero P, so no in-tree of that
    # conformation reaches its free E without half 2's vanishing steps. M binds E(A) and
    # the covalent form :OnlyI, and M-bound half 1 runs at rapid equilibrium in the active
    # conformation, where the vanishing M binding at E(A) lets free E's fast class leave,
    # one way, for the M-bound covalent form. Each conformation's vanishing transitions
    # alone leave the cheapest spanning in-tree at 4 slow transitions; together they cut
    # it to 3. The guard rejects the interior first: the M bindings, which the inactive
    # conformation keeps, join the two sides of M-bound half 1, an :OnlyA step at rapid
    # equilibrium, so that conformation cannot switch it off.
    am = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: M
        steps: begin
            E + A ⇌ E(A)                                             :: EqualAI
            E(A) <--> E(; residual = A - P) + P                      :: EqualAI
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)       :: EqualAI
            E(B; residual = A - P) <--> E + Q                        :: OnlyA
            E(A) + M ⇌ E(A, M)                                       :: OnlyI
            E(; residual = A - P) + M <--> E(M; residual = A - P)    :: OnlyI
            E(A, M) ⇌ E(M; residual = A - P) + P                     :: OnlyA
        end
    end
    NO, BR, U = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._BRIDGE,
        EnzymeRates.Derivation._LIMIT_UNTESTED
    k = _testhelper_edge(am, "EM_res_plusA_minusP + P ⇌ EAM")
    v = EnzymeRates.Derivation._guard(am)
    @test v == EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._FAST_STEP_PINNED, U,
                                               0, 2, k)
    @test EnzymeRates.Derivation._guard_message(am, v) ==
          _testhelper_pinned_cause(am, k, 2) * " Make its group steady-state."
    @test _testhelper_epath_verdict(am; rng = Random.MersenneTwister(95)) === :no_path
    @test _testhelper_admission_faults(am) == [[NO, BR, NO, NO, NO, BR, NO]; fill(NO, 7)]
    eps = EnzymeRates.Derivation._eps_network(
        am, EnzymeRates.Derivation._re_flags(am.catalytic),
        EnzymeRates.Derivation._built_reading(am))
    law = EnzymeRates.Derivation._law_network(am)
    faces = EnzymeRates.Derivation._faces(am.catalytic)
    for i in (2, 6)   # zero P and Q, zero P
        s = EnzymeRates.Derivation._face_sets(eps, law.present, faces[i])
        nodes = s.fed .& .!s.stranded
        @test EnzymeRates.Derivation._fewest_slow(eps, nodes, s.n0) == 4
        for c in 1:2
            (; van) = EnzymeRates.Derivation._copy_transitions(eps, s.alive, nodes, c)
            @test EnzymeRates.Derivation._fewest_slow(eps, nodes, s.n0 .| van) == 4
        end
        @test EnzymeRates.Derivation._fewest_slow(eps, nodes, s.alive) == 3
        @test EnzymeRates.Derivation._face_condition(
            eps, law.present, faces[i],
            EnzymeRates.Derivation._turnover_name(am.catalytic)) ==
              EnzymeRates.Derivation._Verdict(BR, U, 0, 0, 0)
    end
    # The face condition's message names the vanishing steps a cheapest in-tree needs: the
    # M binding at E(A) of the active conformation and half 2 of the inactive one.
    w = EnzymeRates.Derivation._admission_face(am, eps, law, faces, 2)
    @test w == EnzymeRates.Derivation._Verdict(BR, U, 2, 0, 0)
    @test EnzymeRates.Derivation._guard_message(am, w) ==
          "AllostericMechanism: at zero P and Q, the tagged steps E(A) + M ⇌ E(A, M); " *
          "E + Q <--> E(B; residual = A - P) of both conformations together bridge forms " *
          "that each conformation's other steps keep apart, so admission cannot show " *
          "that the limit does not depend on how fast those steps vanish."
end

@testset "_redundant_copy_groups builds each conformation from its own law" begin
    # An aggregate regression pin over the thirteen :OnlyI fixtures: the copy rule reads
    # the active conformation without its :OnlyI groups and the inactive one without its
    # :OnlyA groups, and it flags none of them.
    for (label, am) in pairs(_testhelper_onlyi_fixtures())
        @test (label, EnzymeRates.Enumeration._redundant_copy_groups(am)) == (label, Int[])
    end
end

@testset "check_valid_rate_equation(Bool, m) returns false on a fault" begin
    # A mechanism whose law, on a face (face 2), has several closed classes in its
    # segment chain (`_LAW_MULTI`) fails the guard; the `Bool` form gives the verdict
    # instead of throwing.
    good = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
    bad = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q <--> E(Q)
            E(P) + Q <--> E(P, Q)
        end
    end
    check = EnzymeRates.Derivation.check_valid_rate_equation
    @test check(good) === true
    @test check(Bool, good) === true
    @test check(Bool, bad) === false
    @test @inferred(check(Bool, good)) === true
    @test @inferred(check(Bool, bad)) === false
    v = EnzymeRates.Derivation._guard(bad)
    @test_throws ErrorException(EnzymeRates.Derivation._guard_message(bad, v)) check(bad)
end
