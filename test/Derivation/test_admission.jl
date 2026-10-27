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
    x = [log(params[p]) for p in EnzymeRates.fitted_params(am)]
    at_nan = FittingProblem(am, _testhelper_table(0.5, 1.0, 0.0); Keq = 1.0, reference)
    at_zero = FittingProblem(am, _testhelper_table(0.5, 0.0, 1.0); Keq = 1.0, reference)
    @test EnzymeRates.loss!(x, at_nan) == EnzymeRates.loss!(x, at_zero)
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
    @test EnzymeRates.Mechanisms._excludes(r2, bind, 2) &&
          !EnzymeRates.Mechanisms._excludes(r2, bind, 1)
    @test EnzymeRates.Mechanisms._law_mask(r2, 1) == trues(4)
    # The inactive law deletes both tagged steps, which strands every form but E.
    @test EnzymeRates.Mechanisms._law_mask(r2, 2) == falses(4)
    law = EnzymeRates.Derivation._law_network(r2)
    @test (law.nforms, law.ncopies, length(law.from)) == (4, 2, 10)
    @test law.present == [true, true, true, true, true, false, false, false]
    @test law.edge[9:10] == [0, 0] && law.conformation[9:10] == [0, 0]
    @test (law.from[9], law.to[9]) ==
          (EnzymeRates.Derivation._node(law, 1, 1), EnzymeRates.Derivation._node(law, 1, 2))
    @test law.fast[9] && law.fast[10] && law.forward[9] && !law.forward[10]
    @test !any(law.vanishing)
    eps = EnzymeRates.Derivation._eps_network(r2)
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
    eps_re = EnzymeRates.Derivation._eps_network(r2, trues(4))
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
        EnzymeRates.Derivation._eps_network(mm)
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
    @test !EnzymeRates.Derivation._has_turnover_cycle(no_turnover, trues(3))
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
        net = EnzymeRates.Derivation._eps_network(m)
        cycle = EnzymeRates.Derivation._re_turnover_transitions(net)
        @test !isempty(cycle) && all(t -> net.fast[t], cycle)
        v = EnzymeRates.Derivation._guard(m)
        @test v.fault === EnzymeRates.Derivation._RE_TURNOVER
        msg = EnzymeRates.Derivation._guard_message(m, v)
        @test all(
            t -> occursin(EnzymeRates.Mechanisms._step_text(m, net.edge[t]), msg), cycle)
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
        EnzymeRates.Derivation._eps_network(sib_fixed)))
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
                 modifiers = EnzymeRates.Mechanisms.modifiers(cm),
                 regulators = mech isa EnzymeRates.AllostericMechanism ?
                              EnzymeRates.Mechanisms.allosteric_regulators(mech) : Symbol[])
        raw = EnzymeRates.Mechanisms._raw(cm)
        for k in eachindex(cm.edges)
            op, lhs, rhs = EnzymeRates.Mechanisms._step_parts(
                Meta.parse(EnzymeRates.Mechanisms._step_text(cm, k)), "@enzyme_mechanism")
            @test EnzymeRates.Mechanisms._raw_step(op, lhs, rhs, names,
                                                   "@enzyme_mechanism") == raw.steps[k]
        end
    end
    line(m, text) = EnzymeRates.Mechanisms._step_text(m, _testhelper_edge(m, text))
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
    v = EnzymeRates.Derivation._guard(am)
    @test (v.fault, v.limit, v.face, v.conformation, v.edge) ==
          (EnzymeRates.Derivation._TRAP, EnzymeRates.Derivation._LIMIT_DIFFER, 2, 0,
           _testhelper_edge(am, "E + P ⇌ EP"))
    @test EnzymeRates.Derivation._guard_message(am, v) ==
          "AllostericMechanism: at zero P and Q, the forms {E (active), E (inactive)} " *
          "hold enzyme, but the rapid-equilibrium step E + P ⇌ E(P) that joins them to " *
          "the rest of their segment takes up a zeroed metabolite, so the rate law " *
          "misweights them. Make that step's kinetic group steady-state."
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
    @test !EnzymeRates.Derivation._is_well_posed(unfed)
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
    leaving = only(k for (k, e) in enumerate(weighted.edges)
                   if weighted.forms[e.from].conformation !=
                      weighted.forms[e.to].conformation)
    @test v.edge == leaving
    @test EnzymeRates.Derivation._guard_message(weighted, v) ==
          "Mechanism: at zero P, the forms {F(S)} keep enzyme away from free E: every " *
          "step that leaves them takes up a zeroed metabolite " *
          "($(EnzymeRates.Mechanisms._step_text(weighted, leaving))), so the steady " *
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
    @test !EnzymeRates.Derivation._is_well_posed(dead_end)
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
    v = EnzymeRates.Derivation._guard(differ)
    @test v == EnzymeRates.Derivation._Verdict(
        EnzymeRates.Derivation._LAW_MULTI, EnzymeRates.Derivation._LIMIT_DIFFER, 2, 0, 0)
    @test EnzymeRates.Derivation._guard_message(differ, v) ==
          "AllostericMechanism: at zero P and Q, the steady state's fast-equilibrium " *
          "limit differs from the rate law."
    @test _testhelper_respeed_verdicts(differ; rng = Random.MersenneTwister(1),
                                       law = :limit)[2] === :mismatch
end

@testset "a bridge of several tagged steps names each one" begin
    # At zero P and Q the cheapest spanning in-tree of the inactive conformation's fed
    # forms uses two of its vanishing :OnlyA steps, so the message names both. The ε-path
    # converges on that face: the rejection is conservative, as the message allows.
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
    @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._BRIDGE, 2, 2)
    @test EnzymeRates.Derivation._guard_message(bridge2, v) ==
          "AllostericMechanism: at zero P and Q, the :OnlyA steps " *
          "E + A ⇌ E(P; residual = A - P); E(Q) + A ⇌ E(A, Q) bridge forms of the " *
          "inactive conformation that its other steps keep apart, so admission cannot " *
          "show that the limit does not depend on how fast those steps vanish."
    @test _testhelper_epath_verdict(bridge2; rng = Random.MersenneTwister(2),
                                    face = Set([:P, :Q])) === :converges
    # No fix is printed: with the leading step's group steady-state the guard still
    # rejects the mechanism.
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
    w = EnzymeRates.Derivation._guard(leading_ss)
    @test (w.fault, w.face, w.conformation) == (EnzymeRates.Derivation._BRIDGE, 2, 2)
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

@testset "the flip move's predicate reads class 3 through the rapid-equilibrium flags" begin
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
    flags = EnzymeRates.Mechanisms._flags(ldh)
    @test !EnzymeRates.Derivation._faces_pass(ldh, flags)
    # The feeding step steady-state through the flags: the graph of the trap test's fix.
    fed = copy(flags)
    fed[_testhelper_edge(ldh, "EP + Q ⇌ EPQ")] = false
    @test EnzymeRates.Derivation._faces_pass(ldh, fed)
    # Every step at rapid equilibrium: a fast cycle turns the reaction over.
    @test !EnzymeRates.Derivation._faces_pass(ldh, trues(length(ldh.edges)))
end

# ── Tag sets admission's interior accepts ───────────────────────────────────

"""
The tag sets that admission's interior accepts, by name, each built without the guard
(`@_testhelper_unguarded`): those that pass every face, and those a face rejects. The
ping-pong tag sets write their covalent form `E(; residual = A - P)`, which gives the graph
a free form `F` would; AD2a_iso_chemNonequal is R1, AD2e is R2 and AD2c is R3 written
again, and AD6 holds one reaction in two groups, which the constructor rejects.
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
        EnzymeRates.AllostericMechanism(EnzymeRates.Mechanisms._raw_mechanism(block, true))
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

@testset "the graph test equals the exact LP on the uni-uni grid" begin
    # Existence holds exactly when no arc stays inside a component and the arcs close no
    # cycle (spec section 7); the test oracle's exact simplex decides the same question
    # from the cycle rows, and its own exclusion structure names the same failure.
    grid = _testhelper_uniuni_grid()
    @test length(grid) == 108
    for am in grid
        @test _testhelper_existence_agrees(am)
    end
end

@testset "a tagged binding that untagged steps join on both sides has no limit" begin
    rng = Random.MersenneTwister(31)
    NL, U = EnzymeRates.Derivation._NO_LIMIT, EnzymeRates.Derivation._LIMIT_UNTESTED
    # The uni-uni S binding :OnlyA: the chemistry and the P binding join E to E(S) in the
    # inactive conformation, so the S affinity cannot vanish there.
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
    @test EnzymeRates.Derivation._runs_toward_products(m, m.edges[chem])
    v = EnzymeRates.Derivation._guard(am)
    @test v == EnzymeRates.Derivation._Verdict(NL, U, 0, 2, bind)
    @test EnzymeRates.Derivation._admission_interior(am) == v
    @test _testhelper_epath_verdict(am; rng) === :no_path
    line(k) = EnzymeRates.Mechanisms._step_text(m, k)
    @test EnzymeRates.Derivation._guard_message(am, v) ==
          "AllostericMechanism: in the inactive conformation, steps outside :OnlyA " *
          "groups still join the two sides of the :OnlyA binding $(line(bind)), so " *
          "thermodynamics keeps its affinity finite and the binding has no limit " *
          "there. Tag its group :EqualAI or :NonequalAI, or tag :OnlyA the steps that " *
          "join its sides: $(line(release)); $(line(chem))."
    # Both fixes pass the guard: the binding untagged, and the steps that join its sides
    # tagged.
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
    # AD1 and AD11: a mixed :OnlyA group whose plain member's two sides untagged steps
    # still join (spec section 7: AD1 and AD11 are rejected). AD3c: a ::Mod copy's binding.
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
    for (x, text) in ((ad1_ss, "EQ + B <--> EBQ"), (ad1_re, "EQ + B ⇌ EBQ"),
                      (ad11, "EQ + B <--> EBQ"), (ad3c, "E + Smod ⇌ ESmod"))
        @test EnzymeRates.Derivation._guard(x) ==
              EnzymeRates.Derivation._Verdict(NL, U, 0, 2, _testhelper_edge(x, text))
        @test _testhelper_exclusion(x, 2).failure === :inside
        @test _testhelper_epath_verdict(x; rng) === :no_path
    end
    # On AD1 SS, AD1 RE and AD3c the message's steps that join the sides run from the
    # binding's unbound side to its bound side, and both fixes pass there too.
    for (x, path) in ((ad1_ss, ("E + Q ⇌ EQ", "E + B ⇌ EB", "EB + Q ⇌ EBQ")),
                      (ad1_re, ("E + Q ⇌ EQ", "E + B ⇌ EB", "EB + Q ⇌ EBQ")),
                      (ad3c, ("E + S ⇌ ES", "ES + Smod ⇌ ESSmod", "ESmod + S ⇌ ESSmod")))
        joins = join((EnzymeRates.Mechanisms._step_text(x.catalytic, _testhelper_edge(x, t))
                      for t in path), "; ")
        @test endswith(
            EnzymeRates.Derivation._guard_message(x, EnzymeRates.Derivation._guard(x)),
            "tag :OnlyA the steps that join its sides: $joins.")
    end
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
    # AD11's message names the steps AD1's does (E + Q ⇌ E(Q); E + B ⇌ E(B);
    # E(B) + Q ⇌ E(B, Q)); both its fixes pass every face, the oracles agreeing.
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
    joins = join((EnzymeRates.Mechanisms._step_text(ad11.catalytic,
                                                     _testhelper_edge(ad11, t))
                  for t in ("E + Q ⇌ EQ", "E + B ⇌ EB", "EB + Q ⇌ EBQ")), "; ")
    @test endswith(
        EnzymeRates.Derivation._guard_message(ad11, EnzymeRates.Derivation._guard(ad11)),
        "tag :OnlyA the steps that join its sides: $joins.")
    for y in (ad1_untagged, ad1_joins_tagged, ad1_re_untagged, ad1_re_joins_tagged,
              ad3c_untagged, ad3c_joins_tagged, ad11_untagged, ad11_joins_tagged)
        @test EnzymeRates.Derivation._guard(y) == EnzymeRates.Derivation._WELL_POSED
    end
    # The message prints a fix only when its variant passes. A dead-end I joins the S
    # binding's sides a second way, through E(I) and E(I, S), so tagging the steps of the
    # path the message would name leaves that binding inside one component: the message
    # offers only the untag fix.
    _testhelper_no_limit_cause(x, k) =
        "AllostericMechanism: in the inactive conformation, steps outside :OnlyA " *
        "groups still join the two sides of the :OnlyA binding " *
        "$(EnzymeRates.Mechanisms._step_text(x.catalytic, k)), so thermodynamics keeps " *
        "its affinity finite and the binding has no limit there."
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
    @test v == EnzymeRates.Derivation._Verdict(NL, U, 0, 2,
                                               _testhelper_edge(deadend, "E + S ⇌ ES"))
    @test _testhelper_epath_verdict(deadend; rng) === :no_path
    @test EnzymeRates.Derivation._guard_message(deadend, v) ==
          _testhelper_no_limit_cause(deadend, v.edge) *
          " Tag its group :EqualAI or :NonequalAI."
    @test EnzymeRates.Derivation._guard(deadend_untagged) ==
          EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(deadend_joins_tagged).fault === NL
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
    @test v == EnzymeRates.Derivation._Verdict(NL, U, 0, 2,
                                               _testhelper_edge(two_bindings, "E + I ⇌ EI"))
    @test _testhelper_epath_verdict(two_bindings; rng) === :no_path
    joins = join((EnzymeRates.Mechanisms._step_text(two_bindings.catalytic,
                                                    _testhelper_edge(two_bindings, t))
                  for t in ("E + P ⇌ EP", "ES <--> EP", "ES + I ⇌ EIS", "EI + S ⇌ EIS")),
                 "; ")
    @test EnzymeRates.Derivation._guard_message(two_bindings, v) ==
          _testhelper_no_limit_cause(two_bindings, v.edge) *
          " Tag :OnlyA the steps that join its sides: $joins."
    @test EnzymeRates.Derivation._guard(two_bindings_untagged).fault === NL
    @test EnzymeRates.Derivation._guard(two_bindings_joins_tagged) ==
          EnzymeRates.Derivation._WELL_POSED
end

@testset "tagged bindings that close a directed cycle have no limit" begin
    rng = Random.MersenneTwister(37)
    U = EnzymeRates.Derivation._LIMIT_UNTESTED
    # The ter cube of seven :OnlyA edges: the untagged E(A) + C, E(C) + B and E(B) + A
    # bindings join E(A) to E(A, B), so its :OnlyA binding of B stays inside one component.
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
    @test _testhelper_exclusion(cube, 2).failure === :inside
    @test EnzymeRates.Derivation._guard(cube) ==
          EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._NO_LIMIT, U, 0, 2,
                                          _testhelper_edge(cube, "EA + B ⇌ EAB"))
    @test _testhelper_epath_verdict(cube; rng) === :no_path
    # Both fixes the message offers pass every face, the oracles agreeing:
    # the B binding at E(A) untagged, and the five untagged bindings
    # that join its sides (E(A) + C, E(C) + A, E(C) + B, E(B) + C, E(B) + A) tagged.
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
    joins = join((EnzymeRates.Mechanisms._step_text(cube.catalytic,
                                                     _testhelper_edge(cube, t))
                  for t in ("EA + C ⇌ EAC", "EC + A ⇌ EAC", "EC + B ⇌ EBC",
                            "EB + C ⇌ EBC", "EB + A ⇌ EAB")), "; ")
    @test endswith(
        EnzymeRates.Derivation._guard_message(cube, EnzymeRates.Derivation._guard(cube)),
        "tag :OnlyA the steps that join its sides: $joins.")
    # Both bindings of the ordered bi-uni :OnlyA: the chemistry and the P binding join
    # E(A, B) back to E, so the arcs E → E(A) → {E, E(A, B), E(P)} close a cycle and the
    # two affinities cannot vanish together.
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
    m = am.catalytic
    a, b, chem, release = (_testhelper_edge(am, t) for t in
                           ("E + A ⇌ EA", "EA + B ⇌ EAB", "EAB <--> EP", "E + P ⇌ EP"))
    @test _testhelper_exclusion(am, 2).failure === :cycle
    @test EnzymeRates.Mechanisms._exclusion_components(am, 2) ==
          [1, 2, 1, 1]  # E, E(A), E(A, B), E(P)
    @test EnzymeRates.Mechanisms._exclusion_components(am, 1) == [1, 1, 1, 1]
    v = EnzymeRates.Derivation._guard(am)
    @test v ==
          EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._ARC_CYCLE, U, 0, 2, b)
    @test _testhelper_epath_verdict(am; rng) === :no_path
    line(k) = EnzymeRates.Mechanisms._step_text(m, k)
    @test EnzymeRates.Derivation._guard_message(am, v) ==
          "AllostericMechanism: the :OnlyA bindings $(line(b)); $(line(a)) close a cycle " *
          "in the inactive conformation: steps outside :OnlyA groups join each one's " *
          "bound side to the next one's unbound side, so thermodynamics forbids their " *
          "affinities to vanish together. Tag them :EqualAI or :NonequalAI, or tag " *
          ":OnlyA the steps that join them: $(line(chem)); $(line(release))."
    # Both fixes pass the guard; untagging one binding alone leaves the other inside a
    # component, so the message does not offer it.
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
    @test EnzymeRates.Derivation._guard(one_untagged).fault ===
          EnzymeRates.Derivation._NO_LIMIT
    # With I forms, I's bindings join E(A, B) back to E a second way, so tagging the
    # chemistry and the P binding leaves the cycle closed: the message offers only the
    # untag fix.
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
    line_i(k) = EnzymeRates.Mechanisms._step_text(with_i.catalytic, k)
    a_i, b_i = (_testhelper_edge(with_i, t) for t in ("E + A ⇌ EA", "EA + B ⇌ EAB"))
    v = EnzymeRates.Derivation._guard(with_i)
    @test v ==
          EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._ARC_CYCLE, U, 0, 2, a_i)
    @test _testhelper_epath_verdict(with_i; rng) === :no_path
    @test EnzymeRates.Derivation._guard_message(with_i, v) ==
          "AllostericMechanism: the :OnlyA bindings $(line_i(a_i)); $(line_i(b_i)) " *
          "close a cycle in the inactive conformation: steps outside :OnlyA groups join " *
          "each one's bound side to the next one's unbound side, so thermodynamics " *
          "forbids their affinities to vanish together. Tag them :EqualAI or :NonequalAI."
    @test EnzymeRates.Derivation._guard(with_i_untagged) ==
          EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(with_i_joins_tagged).fault ===
          EnzymeRates.Derivation._ARC_CYCLE
end

@testset "forms reached only through tagged chemistry are an island" begin
    rng = Random.MersenneTwister(41)
    IS, U = EnzymeRates.Derivation._ISLAND, EnzymeRates.Derivation._LIMIT_UNTESTED
    _testhelper_edge_between(x, a, b) = only(k for (k, e) in enumerate(x.catalytic.edges)
        if Set(EnzymeRates.Mechanisms._form_name.(Ref(x.catalytic), (e.from, e.to))) ==
           Set((a, b)))
    _testhelper_island_verdict(x) =
        _testhelper_epath_verdict(x; rng) in (:direction_dependent, :wrong_limit)
    # Both half-reactions of the ping-pong :OnlyA: the covalent forms join free E only
    # through chemistry. The message names half 1, which runs toward products into the
    # island (it is stored as the uptake of P), and its untagging leaves half 2 :OnlyA, a
    # well-posed variant.
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
    @test !EnzymeRates.Derivation._runs_toward_products(both.catalytic,
                                                        both.catalytic.edges[half1])
    v = EnzymeRates.Derivation._guard(both)
    @test v == EnzymeRates.Derivation._Verdict(IS, U, 0, 2, half1)
    @test !isempty(_testhelper_exclusion(both, 2).lost)
    @test _testhelper_island_verdict(both)
    line = EnzymeRates.Mechanisms._step_text(both.catalytic, half1)
    @test EnzymeRates.Derivation._guard_message(both, v) ==
          "AllostericMechanism: in the inactive conformation the forms " *
          "E(; residual = A - P), E(B; residual = A - P) join free E only through " *
          ":OnlyA steps that are not bindings toward them, such as $line, so their " *
          "weight depends on how fast those steps vanish. Tag the group of $line " *
          ":EqualAI or :NonequalAI. An inactive conformation that runs neither " *
          "half-reaction of a ping-pong is not supported yet (see the Roadmap page)."
    @test EnzymeRates.Derivation._guard(_testhelper_tag_sets["pp_V_half2"]) ==
          EnzymeRates.Derivation._WELL_POSED
    # R3 and R4: a second free form behind a tagged isomerization; untagging it passes.
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
    for (x, fixed, forms) in ((r3, r3_untagged, "F, F(I)"), (r4, r4_untagged, "F"))
        iso = _testhelper_edge_between(x, "E", "F")
        # An isomerization of free forms ties on both sides, so it runs as stored.
        @test EnzymeRates.Derivation._runs_toward_products(x.catalytic,
                                                           x.catalytic.edges[iso])
        v = EnzymeRates.Derivation._guard(x)
        @test v == EnzymeRates.Derivation._Verdict(IS, U, 0, 2, iso)
        @test _testhelper_island_verdict(x)
        msg = EnzymeRates.Derivation._guard_message(x, v)
        line = EnzymeRates.Mechanisms._step_text(x.catalytic, iso)
        @test occursin("the forms $forms join free E", msg)
        @test endswith(msg, "Tag the group of $line :EqualAI or :NonequalAI.")
        @test EnzymeRates.Derivation._guard(fixed) == EnzymeRates.Derivation._WELL_POSED
    end
    # AD2b and the ping-pong with B's binding and half 1 :OnlyA: a tagged binding borders
    # the island too, so untagging the named step leaves that binding inside a component,
    # and the message offers no fix.
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
    for (x, named) in ((ad2b, _testhelper_edge_between(ad2b, "E", "F")),
                       (bbind, _testhelper_edge(bbind, "E_res_plusA_minusP + P <--> EA")))
        v = EnzymeRates.Derivation._guard(x)
        @test v == EnzymeRates.Derivation._Verdict(IS, U, 0, 2, named)
        @test _testhelper_island_verdict(x)
        @test !occursin("Tag the group", EnzymeRates.Derivation._guard_message(x, v))
    end
    # Half 1 and its copy at E(A, I) share one :OnlyA group, and I's binding at E is
    # :OnlyA too. Untagging that group joins E(I) to E through the island, so I's binding
    # falls inside one component and the message offers no fix.
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
    v = EnzymeRates.Derivation._guard(shared_half1)
    @test v == EnzymeRates.Derivation._Verdict(
        IS, U, 0, 2, _testhelper_edge(shared_half1, "E_res_plusA_minusP + P <--> EA"))
    @test _testhelper_island_verdict(shared_half1)
    @test !occursin("Tag the group", EnzymeRates.Derivation._guard_message(shared_half1, v))
    @test EnzymeRates.Derivation._guard(shared_half1_untagged).fault ===
          EnzymeRates.Derivation._NO_LIMIT
end

@testset "the admitted tag sets pass admission's interior" begin
    @testset "$name" for name in sort!(collect(keys(_testhelper_tag_sets)))
        @test EnzymeRates.Derivation._admission_interior(_testhelper_tag_sets[name]) ==
              EnzymeRates.Derivation._WELL_POSED
    end
end

@testset "a mixed :OnlyA group whose limit exists is admitted" begin
    # B's fused binding into E(P, Q) and its plain binding at E(Q) share one :OnlyA group;
    # in the inactive conformation no untagged step joins E(Q) to E(B, Q), so the plain
    # binding's affinity can vanish, and the limit is deletion plus pruning.
    am = _testhelper_tag_sets["mixed_group_limit"]
    @test _testhelper_exclusion(am, 2).failure === :none
    @test _testhelper_lp_existence(am, 2)
    @test EnzymeRates.Derivation._admission_interior(am) ==
          EnzymeRates.Derivation._WELL_POSED
    @test EnzymeRates.Derivation._guard(am) == EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_epath_verdict(am; rng = Random.MersenneTwister(43)) === :converges
end

"""The tag sets of `_testhelper_tag_sets` that a face rejects; every other one passes every
face."""
const _testhelper_face_rejected = Set(["R2", "BS1", "BS2", "BS3", "BS4", "ADbridge",
    "ADbridge_noturnover", "bridge_drain", "bridge_deadend", "R2funnel", "R2funnel_SS",
    "S1", "AD7b", "fz101_10232"])

@testset "R2's limit at zero product depends on how fast its isomerization vanishes" begin
    rng = Random.MersenneTwister(43)
    NO, PU = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._PUMPED
    # Classes 1-3 and admission's interior accept R2, so the face condition alone rejects
    # it (spec section 7): at P zero, F's exits vanish (F <--> E) or take up P (F + P ⇌
    # E(P)), and E(P) feeds it by a release whose reverse is dead.
    am = _testhelper_tag_sets["R2"]
    m = am.catalytic
    @test _testhelper_face_faults(am) == fill(NO, 3)
    @test EnzymeRates.Derivation._admission_interior(am) ==
          EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_admission_faults(am) == [NO, PU, NO]
    @test _testhelper_epath_verdict(am; rng, face = _testhelper_faces(am)[2]) ===
          :direction_dependent
    @test all(_testhelper_admission_agreement(am; rng))
    k, iso = _testhelper_edge(am, "F + P ⇌ EP"), _testhelper_edge(am, "F <--> E")
    v = EnzymeRates.Derivation._guard(am)
    @test v == EnzymeRates.Derivation._Verdict(
        PU, EnzymeRates.Derivation._LIMIT_UNTESTED, 2, 2, k)
    line(j) = EnzymeRates.Mechanisms._step_text(m, j)
    @test EnzymeRates.Derivation._guard_message(am, v) ==
          "AllostericMechanism: at zero P, $(line(k)) feeds the forms F of the inactive " *
          "conformation, which leave only through the vanishing step $(line(iso)) while " *
          "nothing runs back the way it came, so admission cannot show that their weight " *
          "does not depend on how fast that step vanishes."
end

@testset "a stranded basin with two vanishing exits names both steps in its message" begin
    # At P zero, E(S) feeds F by a release whose reverse takes up P; F and F(I) leave the
    # inactive conformation's basin only through the :OnlyA steps E(I) ⇌ F(I) and
    # F <--> E, which vanish there. The ε-path oracle finds that limit direction-dependent.
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
    @test EnzymeRates.Derivation._admission_interior(am) ==
          EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_epath_verdict(am; rng, face = _testhelper_faces(am)[2]) ===
          :direction_dependent
    v = EnzymeRates.Derivation._guard(am)
    @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._PUMPED, 2, 2)
    @test EnzymeRates.Derivation._guard_message(am, v) ==
          "AllostericMechanism: at zero P, F + P <--> E(S) feeds the forms F, F(I) of " *
          "the inactive conformation, which leave only through the vanishing steps " *
          "E(I) ⇌ F(I); F <--> E while nothing runs back the way it came, so admission " *
          "cannot show that their weight does not depend on how fast those steps vanish."
end

@testset "admitted tag sets converge on every face" begin
    rng = Random.MersenneTwister(47)
    names = sort!([n for n in keys(_testhelper_tag_sets) if n ∉ _testhelper_face_rejected])
    @test length(names) == 39
    @testset "$name" for name in names
        am = _testhelper_tag_sets[name]
        @test EnzymeRates.Derivation._guard(am) == EnzymeRates.Derivation._WELL_POSED
        @test all(==(EnzymeRates.Derivation._NO_FAULT), _testhelper_admission_faults(am))
        @test all(_testhelper_admission_agreement(am; rng))
    end
end

@testset "a sealed basin drains, and a pumped one does not" begin
    rng = Random.MersenneTwister(53)
    NO, PU = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._PUMPED
    function _testhelper_stranded(am, i)
        face = EnzymeRates.Derivation._faces(am.catalytic)[i]
        s = EnzymeRates.Derivation._face_sets(
            EnzymeRates.Derivation._eps_network(am),
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
    # A basin fed by a step whose reverse is dead is pumped: BS2, BS3 and BS4 at P zero.
    @test _testhelper_admission_faults(_testhelper_tag_sets["BS2"]) == [NO, PU, NO]
    @test _testhelper_admission_faults(_testhelper_tag_sets["BS3"]) == [NO, PU, NO]
    @test _testhelper_admission_faults(_testhelper_tag_sets["BS4"]) ==
          [NO, PU, NO, NO, NO, NO]
    for name in ("R2P_sealed", "B_uphill", "BS2", "BS3", "BS4")
        @test all(_testhelper_admission_agreement(_testhelper_tag_sets[name]; rng))
    end
end

@testset "a tagged binding that bridges two basins is rejected" begin
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
    line = EnzymeRates.Mechanisms._step_text(m, k)
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
        eps, law = EnzymeRates.Derivation._eps_network(x),
            EnzymeRates.Derivation._law_network(x)
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
    # R2funnel: at P zero the trap {F, F(S)} drains only through a slow step out of a fast
    # transient; its bridge stays with the S binding steady-state (R2funnel_SS), so the
    # message offers no fix.
    for name in ("R2funnel", "R2funnel_SS")
        x = _testhelper_tag_sets[name]
        @test _testhelper_admission_faults(x) == [NO, BR, NO]
        @test all(_testhelper_admission_agreement(x; rng))
        v = EnzymeRates.Derivation._guard(x)
        @test (v.fault, v.face, v.conformation) === (BR, 2, 2)
        @test !occursin("Make the group", EnzymeRates.Derivation._guard_message(x, v))
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
    eps, law = EnzymeRates.Derivation._eps_network(am),
        EnzymeRates.Derivation._law_network(am)
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

@testset "conservative rejections the oracle admits stay rejected" begin
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
                   "$(EnzymeRates.Mechanisms._step_text(am.catalytic, k)) steady-state.")
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
    # R-B2, a stranded basin fed only through vanishing sources: BS1 with I zero and P or
    # S zero, fz101_10232 on its two faces with S zero. Every corner direction converges.
    bs1 = _testhelper_tag_sets["BS1"]
    @test _testhelper_admission_faults(bs1) == [NO, NO, NO, NO, PU, PU]
    @test all(i -> _testhelper_corner_verdict(bs1, i) === :converges, 5:6)
    @test all(_testhelper_admission_agreement(bs1; rng)[1:4])
    # The step that rejects BS1 has both its forms in the basin: it runs one way among
    # them rather than feeding them.
    v = EnzymeRates.Derivation._guard(bs1)
    @test (v.fault, v.face, v.conformation) == (PU, 5, 2)
    @test EnzymeRates.Derivation._guard_message(bs1, v) ==
          "AllostericMechanism: at zero I and P, F + P ⇌ F(P) runs one way among the " *
          "forms F, F(P), F(S) of the inactive conformation, which leave only through " *
          "the vanishing step E <--> F while nothing runs back the way it came, so " *
          "admission cannot show that their weight does not depend on how fast that " *
          "step vanishes."
    fz = _testhelper_tag_sets["fz101_10232"]
    @test _testhelper_admission_faults(fz) == [NO, PU, PU, NO, PU, PU]
    @test _testhelper_corner_verdict(fz, 3) === :converges
    # On both faces with P zero, I zero or not (faces 2 and 5), the rejection is right: at
    # face 2 F's only alive exit, the blocked F <--> E, can vanish faster than its inflow,
    # which one corner direction shows.
    @test _testhelper_corner_verdict(fz, 2) !== :converges
    @test all(_testhelper_admission_agreement(fz; rng)[[1, 4]])
    # The faces with I zero as well come last, so the earlier draws of `rng` stay put.
    @test _testhelper_corner_verdict(fz, 5) !== :converges
    @test _testhelper_corner_verdict(fz, 6) === :converges
end

@testset "S1 and AD7b pass admission's interior and fail class 3" begin
    rng = Random.MersenneTwister(71)
    NO, NU = EnzymeRates.Derivation._NO_FAULT, EnzymeRates.Derivation._NONUNIQUE
    # S1's stray free form F binds B at steady state: with A or B zero the law network has
    # a second closed class. The ε-network's closed-class count sees the same split, so
    # admission alone rejects those faces too; the guard reports class 3, which runs first.
    s1 = _testhelper_tag_sets["S1"]
    @test EnzymeRates.Derivation._admission_interior(s1) ==
          EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_admission_faults(s1) == [NO, NO, NU, NU, NU]
    @test EnzymeRates.Derivation._guard(s1).fault === NU
    @test all(_testhelper_admission_agreement(s1; rng))
    v = EnzymeRates.Derivation._admission(s1)
    @test (v.fault, v.face, v.conformation) ===
          (EnzymeRates.Derivation._EPS_NONUNIQUE, 3, 2)
    @test EnzymeRates.Derivation._guard_message(s1, v) ==
          "AllostericMechanism: at zero A and B, the forms F (inactive) hold the enzyme " *
          "they receive and are reached only through tagged steps, so the steady state " *
          "depends on where the enzyme starts."
    # AD7b: with I zero no spanning in-tree reaches the fed forms, so the face condition
    # passes and class 3 decides.
    ad7b = _testhelper_tag_sets["AD7b"]
    @test EnzymeRates.Derivation._admission(ad7b) == EnzymeRates.Derivation._WELL_POSED
    @test _testhelper_admission_faults(ad7b) == [NO, NO, NO, NU, NU, NU]
    @test all(_testhelper_admission_agreement(ad7b; rng))
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

@testset "the guard rejects a steady-state dead end whose law vanishes at zero products" begin
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
    @test !EnzymeRates.Derivation._is_well_posed(m)
    v = EnzymeRates.Derivation._guard(m)
    @test (v.fault, v.face) == (EnzymeRates.Derivation._LAW_MULTI, 2)
    # The same dead end with E(P) + Q <--> E(P, Q) in the active conformation only. The
    # active denominator holds P or Q in every term, as above, and the inactive one is
    # weighted by the active free-enzyme weight, which holds Q in every term: the law is
    # 0/0 with P and Q zero, though the inactive denominator alone keeps a constant term.
    # In the inactive conformation E(P) joins free E only through that :OnlyA step, an
    # island the guard rejects before it checks any face.
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
    @test !EnzymeRates.Derivation._is_well_posed(am)
    w = EnzymeRates.Derivation._guard(am)
    @test (w.fault, w.face) == (EnzymeRates.Derivation._ISLAND, 0)
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
        seeds = EnzymeRates.init_mechanisms(rxn)
        append!(mechs, seeds, EnzymeRates.expand_mechanisms(seeds, rxn))
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
    @test occursin(EnzymeRates.Mechanisms._step_text(m, feed),
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
    @test length(EnzymeRates.fitted_params(spec.mechanism)) == 6
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
    @test occursin(EnzymeRates.Mechanisms._step_text(am.catalytic, feed),
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
    @test length(EnzymeRates.fitted_params(spec.mechanism)) == 10
end

# ── The guard wired: macros, exits, flip move and parser ──────────────────────

@testset "the macros reject what the guard rejects when the expanded code runs" begin
    # One mechanism per class: no cycle turns the reaction over (class 1); a cycle of
    # rapid-equilibrium steps does (class 2); a trap at zero products (class 3); an
    # :OnlyA binding whose two sides untagged steps join (class 4). Each `@test_throws`
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
                     EnzymeRates.Derivation._TRAP, EnzymeRates.Derivation._NO_LIMIT]
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
    # half-reactions :OnlyA. Read back, such text raises the guard's message, not a parse
    # error, and the message names the half-reaction whose untagging fixes it.
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
    v = EnzymeRates.Derivation._guard(both)
    @test v.fault === EnzymeRates.Derivation._ISLAND
    msg = EnzymeRates.Derivation._guard_message(both, v)
    @test occursin("Tag the group of E(A) <--> E(P; residual = A - P) :EqualAI or " *
                   ":NonequalAI.", msg)
    @test occursin("is not supported yet (see the Roadmap page)", msg)
    @test_throws ErrorException(msg) EnzymeRates._parse_mechanism(sprint(show, both))
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
    @test EnzymeRates.Mechanisms._active_mask(onlyi) ==
          EnzymeRates.Mechanisms._inactive_mask(onlya)
    @test EnzymeRates.Mechanisms._inactive_mask(onlyi) ==
          EnzymeRates.Mechanisms._active_mask(onlya)
    @test !all(EnzymeRates.Mechanisms._active_mask(onlyi))
end

@testset "an RE binding whose affinity vanishes still closes a fast cycle" begin
    # The B1 set: an :OnlyI binding at rapid equilibrium loses its affinity in the active
    # conformation, yet on the ε-network it stays fast and closes, with the other
    # rapid-equilibrium steps, a cycle that turns the reaction over: class 2 rejects on
    # the interior. The law networks hold no such cycle, so class 3 and the RE-speed
    # oracle agree on every face, while the ε-path does not converge (checked on B1a).
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
    @test _testhelper_epath_verdict(b1a; rng = Random.MersenneTwister(62)) in
          (:diverges, :direction_dependent, :wrong_limit)
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
    # rapid-equilibrium steps run, and passes the other faces.
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
    @test v.face == 2 && v.limit === EnzymeRates.Derivation._LIMIT_SPEED
    @test _testhelper_respeed_verdicts(ssbind; rng = Random.MersenneTwister(64),
                                       law = :limit) ==
          [:agree, :pattern, :agree, :agree, :agree, :agree]
end

@testset "a fast cycle that only the inactive conformation keeps is named there" begin
    # The two :OnlyI steps each take up I and J at once, so the active conformation blocks
    # them, while the inactive one keeps them fast: the only cycle of rapid-equilibrium
    # steps, which turns the reaction over, runs in the inactive conformation, through
    # the :OnlyA product binding, whose affinity vanishes there but stays fast.
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
    v = EnzymeRates.Derivation._guard(m)
    @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._RE_TURNOVER, 0, 2)
    @test EnzymeRates.Derivation._guard_message(m, v) ==
          "AllostericMechanism: the rapid-equilibrium steps E(I, J, P) ⇌ E(P) + I + J, " *
          "E + P ⇌ E(P) (whose affinity vanishes in the inactive conformation), " *
          "E + S ⇌ E(S), E(I, J, S) ⇌ E(S) + I + J and E(I, J, S) ⇌ E(I, J, P) close a " *
          "cycle in the inactive conformation that turns the reaction over, so the " *
          "mechanism has no finite rate. Make one step of the cycle steady-state."
    # The mechanism also fails admission: its I, J complexes form an island in the active
    # conformation, so the printed fix only breaks the cycle.
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
    @test EnzymeRates.Derivation._guard(iso_ss).fault === EnzymeRates.Derivation._ISLAND
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
    eps, law = EnzymeRates.Derivation._eps_network(fz),
        EnzymeRates.Derivation._law_network(fz)
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
    # the active conformation, which alone the convention lets catalyze, turns nothing
    # over. Tagging the X binding :EqualAI or :NonequalAI, the message's fix, passes.
    ea = _testhelper_onlyi_fixtures().essential_activator
    v = EnzymeRates.Derivation._guard(ea)
    @test (v.fault, v.face, v.conformation) ==
          (EnzymeRates.Derivation._NO_ACTIVE_TURNOVER, 0, 1)
    @test EnzymeRates.Derivation._guard_message(ea, v) ==
          "AllostericMechanism: without its :OnlyI steps E + X ⇌ E(X) the active " *
          "conformation turns no substrate into product, and only the active " *
          "conformation may catalyze. Tag their groups :EqualAI or :NonequalAI."
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
    ks = [k for k in eachindex(am.catalytic.edges)
          if EnzymeRates.Mechanisms._excludes(am, k, 1)]
    @test !EnzymeRates.Derivation._retag_passes(am, ks, :EqualAI)
    @test !EnzymeRates.Derivation._retag_passes(am, ks, :NonequalAI)
    @test EnzymeRates.Derivation._guard_message(am, v) ==
          "AllostericMechanism: without its :OnlyI steps E + X ⇌ E(X); E(S) + X ⇌ " *
          "E(S, X) the active conformation turns no substrate into product, and only " *
          "the active conformation may catalyze."
end

@testset "the active conformation's interior faults name the :OnlyI step" begin
    # The mirror of the :OnlyA interior faults, in conformation 1. The messages offer no
    # fix: a step tagged :OnlyI binds modifiers alone, and no test shows an untagged
    # variant well-posed. The ε-path oracle agrees with each verdict.
    rng = Random.MersenneTwister(91)
    U = EnzymeRates.Derivation._LIMIT_UNTESTED
    # I binds E and E(S); without E + I ⇌ E(I), the untagged I binding at E(S) and S
    # binding at E(I) still join E to E(I), so the I affinity cannot vanish.
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
    @test v == EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._NO_LIMIT, U, 0, 1, k)
    @test EnzymeRates.Derivation._guard_message(no_limit, v) ==
          "AllostericMechanism: in the active conformation, steps outside :OnlyI " *
          "groups still join the two sides of the :OnlyI binding " *
          "$(EnzymeRates.Mechanisms._step_text(no_limit.catalytic, k)), so " *
          "thermodynamics keeps its affinity finite and the binding has no limit there."
    @test _testhelper_epath_verdict(no_limit; rng) === :no_path
    # I binds E and J binds E(I), both :OnlyI; the untagged J binding at E and I binding
    # at E(J) join E(I, J) back to E, so the two arcs close a cycle.
    arc_cycle = @_testhelper_unguarded @allosteric_mechanism begin
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
    arcs = [_testhelper_edge(arc_cycle, t) for t in ("E + I ⇌ EI", "EI + J ⇌ EIJ")]
    v = EnzymeRates.Derivation._guard(arc_cycle)
    @test (v.fault, v.face, v.conformation) == (EnzymeRates.Derivation._ARC_CYCLE, 0, 1)
    @test v.edge in arcs
    msg = EnzymeRates.Derivation._guard_message(arc_cycle, v)
    @test startswith(msg, "AllostericMechanism: the :OnlyI bindings ")
    @test all(j -> occursin(EnzymeRates.Mechanisms._step_text(arc_cycle.catalytic, j), msg),
              arcs)
    @test endswith(msg, " close a cycle in the active conformation: steps outside " *
                        ":OnlyI groups join each one's bound side to the next one's " *
                        "unbound side, so thermodynamics forbids their affinities to " *
                        "vanish together.")
    @test _testhelper_epath_verdict(arc_cycle; rng) === :no_path
    # E(I) is reached only as E(I, J) gives off J, and in the active conformation J's
    # affinity for E(I) vanishes: an island on the :OnlyI binding's unbound side, which
    # collects the active enzyme.
    island = @_testhelper_unguarded @allosteric_mechanism begin
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
    k = _testhelper_edge(island, "EI + J ⇌ EIJ")
    v = EnzymeRates.Derivation._guard(island)
    @test v == EnzymeRates.Derivation._Verdict(EnzymeRates.Derivation._ISLAND, U, 0, 1, k)
    @test EnzymeRates.Derivation._guard_message(island, v) ==
          "AllostericMechanism: in the active conformation the forms E(I) join free E " *
          "only through :OnlyI steps that are not bindings toward them, such as " *
          "$(EnzymeRates.Mechanisms._step_text(island.catalytic, k)), so their weight " *
          "depends on how fast those steps vanish."
    @test _testhelper_epath_verdict(island; rng) in (:direction_dependent, :wrong_limit)
end

@testset "joint R1 rejects a mechanism whose bridges cut the tree only together" begin
    # The hand-built joint-R1 network's pattern on a mechanism. Half 2 :OnlyA strands the
    # inactive conformation's covalent forms at zero P, so no in-tree of that
    # conformation reaches its free E without half 2's vanishing steps. M binds E(A) and
    # the covalent form :OnlyI, and M-bound half 1 runs at rapid equilibrium in the active
    # conformation, where the vanishing M binding at E(A) lets free E's fast class leave,
    # one way, for the M-bound covalent form. Each conformation's vanishing transitions
    # alone leave the cheapest spanning in-tree at 4 slow transitions; together they cut
    # it to 3. Every other check passes, and the ε-path does not converge on those faces,
    # so the rejection is right; both oracles agree on every face.
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
    v = EnzymeRates.Derivation._guard(am)
    @test v == EnzymeRates.Derivation._Verdict(BR, U, 2, 0, 0)
    @test _testhelper_admission_faults(am) == [[NO, BR, NO, NO, NO, BR, NO]; fill(NO, 7)]
    eps, law = EnzymeRates.Derivation._eps_network(am),
        EnzymeRates.Derivation._law_network(am)
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
    # The message names the vanishing steps a cheapest in-tree needs: the M binding at
    # E(A) of the active conformation and half 2 of the inactive one.
    @test EnzymeRates.Derivation._guard_message(am, v) ==
          "AllostericMechanism: at zero P and Q, the tagged steps E(A) + M ⇌ E(A, M); " *
          "E + Q <--> E(B; residual = A - P) of both conformations together bridge forms " *
          "that each conformation's other steps keep apart, so admission cannot show " *
          "that the limit does not depend on how fast those steps vanish."
    @test all(_testhelper_admission_agreement(am; rng = Random.MersenneTwister(95)))
end

@testset "_redundant_copy_groups builds each conformation from its own law" begin
    # An aggregate regression pin over the thirteen :OnlyI fixtures: the copy rule reads
    # the active conformation without its :OnlyI groups and the inactive one without its
    # :OnlyA groups, and it flags none of them.
    for (label, am) in pairs(_testhelper_onlyi_fixtures())
        @test (label, EnzymeRates._redundant_copy_groups(am)) == (label, Int[])
    end
end
